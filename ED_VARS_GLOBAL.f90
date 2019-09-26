MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE ED_SPARSE_MATRIX
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type effective_bath_component
     real(8)                                       :: v     !spin-keep hyb.
     real(8),dimension(:,:,:,:,:,:),allocatable    :: h     !Replica hamilt [Nlat][Nlat][Nspin][Nspin][Norb][Norb]
  end type effective_bath_component

  type effective_bath
     integer                                                 :: Nmask
     type(effective_bath_component),dimension(:),allocatable :: item
     logical(8),dimension(:,:,:,:,:,:),allocatable           :: mask  !impHloc mask   [Nlat][Nlat][Nspin][Nspin][Norb][Norb]
     logical                                                 :: status=.false.
  end type effective_bath




  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!
  type sector_map
     integer,dimension(:),allocatable :: map
     logical                          :: status=.false.
  end type sector_map

  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate



  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  !dbleMat*dbleVec
  abstract interface
     subroutine dd_sparse_HxV(Nloc,v,Hv)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: v
       real(8),dimension(Nloc) :: Hv
     end subroutine dd_sparse_HxV
  end interface


  !-------------- GMATRIX FOR FAST EVALUATION OF GF ------------------!
  !note that we use a single Qmatrix here which must be intended as
  !component corresponding to the poles. 
  type GFspectrum
     complex(8),dimension(:),allocatable :: weight
     complex(8),dimension(:),allocatable :: poles
  end type GFspectrum
  type GFchannel
     type(GFspectrum),dimension(:),allocatable :: channel !N_channel = 2 (c,cdag), 4 (c,cdag,c pm cdag)
  end type GFchannel
  type GFmatrix
     type(GFchannel),dimension(:),allocatable :: state !state_list%size = # of state in the spectrum 
  end type GFmatrix


  interface GFmatrix_allocate
     module procedure :: allocate_GFmatrix_Nstate
     module procedure :: allocate_GFmatrix_Nchan
     module procedure :: allocate_GFmatrix_Nexc
  end interface GFmatrix_allocate



  !-------------------------- ED  VARIABLES --------------------------!

  !SIZE OF THE PROBLEM
  !=========================================================
  integer,save                                       :: Ns       !Number of levels per spin
  integer,save                                       :: Nsectors !Number of sectors
  integer,save                                       :: Ns_orb
  integer,save                                       :: Ns_ud
  integer                                            :: Nlat     !# size of cluster
  !
  integer                                            :: Nlso

  !local part of the Hamiltonian
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  real(8),dimension(:,:,:,:,:,:),allocatable         :: impHloc           !local hamiltonian [Nlat][Nlat][Nspin][Nspin][Norb][Norb]


  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:)                   :: getDim             ! [Nsectors]
  integer,allocatable,dimension(:,:,:)               :: getCsector         ! [1/Norb,2,NSectors]
  integer,allocatable,dimension(:,:,:)               :: getCDGsector       ! [1/Norb,2,NSectors]
  integer,allocatable,dimension(:,:,:)               :: getBathStride
  ! integer,allocatable,dimension(:,:)                 :: impIndex
  logical,allocatable,dimension(:)                   :: twin_mask
  logical,allocatable,dimension(:)                   :: sectors_mask

  !Effective Bath used in the ED code (this is opaque to user)
  !PRIVATE
  !=========================================================
  type(effective_bath)                               :: dmft_bath


  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  type(sparse_matrix_csr)                            :: spH0d !diagonal part
  type(sparse_matrix_csr)                            :: spH0nd !non-diagonal part
  type(sparse_matrix_csr),dimension(:),allocatable   :: spH0ups,spH0dws !reduced UP and DW parts
  !
  procedure(dd_sparse_HxV),pointer                   :: spHtimesV_p=>null()


  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  integer,allocatable,dimension(:)                   :: neigen_sector
  !--------------- LATTICE WRAP VARIABLES -----------------!  
  integer,allocatable,dimension(:,:)                 :: neigen_sectorii
  integer,allocatable,dimension(:)                   :: neigen_totalii
  logical                                            :: trim_state_list=.false.

  !Partition function
  !PRIVATE
  !=========================================================
  real(8)                                            :: zeta_function
  real(8)                                            :: gs_energy



  !Impurity Green's function and Self-Energies: (Nlat,Nlat,Nspin,Nspin,Norb,Norb,:)
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impGmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impGreal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  !
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impG0mats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impG0real ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  !
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impSmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impSreal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  !
  type(GFmatrix),allocatable,dimension(:,:,:,:,:,:) :: impGmatrix

  ! !--------------- LATTICE WRAP VARIABLES -----------------!
  ! complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Smatsii,Srealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  ! complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Gmatsii,Grealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  ! complex(8),dimension(:,:,:,:,:,:),allocatable,save :: G0matsii,G0realii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  ! complex(8),dimension(:,:,:,:,:)  ,allocatable,save :: imp_density_matrix_ii    ![Nlat][Nspin][Nspin][Norb][Norb]



  ! !Spin Susceptibilities
  ! !=========================================================
  ! real(8),allocatable,dimension(:,:)                 :: spinChi_tau
  ! complex(8),allocatable,dimension(:,:)              :: spinChi_w
  ! complex(8),allocatable,dimension(:,:)              :: spinChi_iv


  ! !Diagonal/Off-diagonal charge-charge Susceptibilities
  ! !=========================================================  
  ! real(8),allocatable,dimension(:,:,:)               :: densChi_tau
  ! complex(8),allocatable,dimension(:,:,:)            :: densChi_w
  ! complex(8),allocatable,dimension(:,:,:)            :: densChi_iv

  ! !Mixed inter-orbital charge-charge Susceptibilities
  ! !=========================================================
  ! real(8),allocatable,dimension(:,:,:)               :: densChi_mix_tau
  ! complex(8),allocatable,dimension(:,:,:)            :: densChi_mix_w
  ! complex(8),allocatable,dimension(:,:,:)            :: densChi_mix_iv

  ! !Total (orbital-sum) Density-density Susceptibilities
  ! !=========================================================
  ! real(8),allocatable,dimension(:)                   :: densChi_tot_tau
  ! complex(8),allocatable,dimension(:)                :: densChi_tot_w
  ! complex(8),allocatable,dimension(:)                :: densChi_tot_iv

  ! !Pair-Pair Susceptibilities
  ! !=========================================================
  ! real(8),allocatable,dimension(:,:)                 :: pairChi_tau
  ! complex(8),allocatable,dimension(:,:)              :: pairChi_w
  ! complex(8),allocatable,dimension(:,:)              :: pairChi_iv


  !Density and double occupancy
  !PRIVATE (now public but accessible thru routines)
  !=========================================================
  real(8),dimension(:,:),allocatable                   ::  ed_dens
  real(8),dimension(:,:),allocatable                   ::  ed_dens_up,ed_dens_dw
  real(8),dimension(:,:),allocatable                   ::  ed_docc
  !--------------- LATTICE WRAP VARIABLES -----------------!
  ! real(8),dimension(:,:),allocatable,save            ::  nii,dii,mii



  !Local energies and generalized double occupancies
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  real(8),dimension(:),allocatable                     :: ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot
  real(8),dimension(:),allocatable                     :: ed_Dust,ed_Dund,ed_Dse,ed_Dph
  ! !--------------- LATTICE WRAP VARIABLES -----------------!
  ! real(8),dimension(:,:),allocatable,save            :: ddii,eii




  ! !Impurity operators
  ! !PRIVATE (now public but accessible thru routine)
  ! !=========================================================
  ! complex(8),allocatable,dimension(:,:,:,:)          :: imp_density_matrix




  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                  :: ed_file_suffix=""       !suffix string attached to the output files.
  character(len=10)                                  :: ineq_site_suffix="_ineq"
  integer                                            :: site_indx_padding=4
  logical                                            :: Jhflag              !spin-exchange and pair-hopping flag.
  logical                                            :: offdiag_gf_flag=.false.
  character(len=200)                                 :: ed_input_file=""


  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                            :: MpiComm_Global=MPI_COMM_NULL
  integer                                            :: MpiComm=MPI_COMM_NULL
#endif
  integer                                            :: MpiGroup_Global=MPI_GROUP_NULL
  integer                                            :: MpiGroup=MPI_GROUP_NULL
  logical                                            :: MpiStatus=.false.
  logical                                            :: MpiMaster=.true.
  integer                                            :: MpiRank=0
  integer                                            :: MpiSize=1
  integer,allocatable,dimension(:)                   :: MpiMembers
  integer                                            :: mpiQup=0
  integer                                            :: mpiRup=0
  integer                                            :: mpiQdw=0
  integer                                            :: mpiRdw=0
  integer                                            :: mpiQ=0
  integer                                            :: mpiR=0
  integer                                            :: mpiIstart
  integer                                            :: mpiIend
  integer                                            :: mpiIshift
  logical                                            :: mpiAllThreads=.true.



contains




  !=========================================================
  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer          :: N
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:)       :: H
    integer,dimension(size(H))          :: N
    integer                             :: i
    do i=1,size(H)
       call map_allocate_scalar(H(i),N(i))
    enddo
  end subroutine map_allocate_vector


  !=========================================================
  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector




  !=========================================================
  subroutine ed_set_MpiComm(comm)
#ifdef _MPI
    integer :: comm,ierr
    ! call MPI_Comm_dup(Comm,MpiComm_Global,ierr)
    ! call MPI_Comm_dup(Comm,MpiComm,ierr)
    MpiComm_Global = comm
    MpiComm        = comm
    call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
    MpiStatus      = .true.
    MpiSize        = get_Size_MPI(MpiComm_Global)
    MpiRank        = get_Rank_MPI(MpiComm_Global)
    MpiMaster      = get_Master_MPI(MpiComm_Global)
#else
    integer,optional :: comm
#endif
  end subroutine ed_set_MpiComm

  subroutine ed_del_MpiComm()
#ifdef _MPI    
    MpiComm_Global = MPI_UNDEFINED
    MpiComm        = MPI_UNDEFINED
    MpiGroup_Global= MPI_GROUP_NULL
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#endif
  end subroutine ed_del_MpiComm





  !Allocate the channels in GFmatrix structure
  subroutine allocate_gfmatrix_Nstate(self,Nstate)
    type(GFmatrix) :: self
    integer        :: Nstate
    if(allocated(self%state))deallocate(self%state)
    allocate(self%state(Nstate))
  end subroutine allocate_gfmatrix_Nstate

  subroutine allocate_gfmatrix_Nchan(self,istate,Nchan)
    type(GFmatrix) :: self
    integer        :: istate,Nchan
    if(allocated(self%state(istate)%channel))deallocate(self%state(istate)%channel)
    allocate(self%state(istate)%channel(Nchan))
  end subroutine allocate_gfmatrix_Nchan

  !Allocate the Excitations spectrum at a given channel
  subroutine allocate_gfmatrix_Nexc(self,istate,ichan,Nexc)
    type(GFmatrix) :: self
    integer        :: istate,ichan
    integer        :: Nexc
    if(allocated(self%state(istate)%channel(ichan)%weight))deallocate(self%state(istate)%channel(ichan)%weight)
    if(allocated(self%state(istate)%channel(ichan)%poles))deallocate(self%state(istate)%channel(ichan)%poles)
    allocate(self%state(istate)%channel(ichan)%weight(Nexc))
    allocate(self%state(istate)%channel(ichan)%poles(Nexc))
  end subroutine allocate_gfmatrix_Nexc


END MODULE ED_VARS_GLOBAL
