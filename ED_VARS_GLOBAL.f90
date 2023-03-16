MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE ED_SPARSE_MATRIX
  USE ED_SPARSE_MAP
  USE ED_INPUT_VARS
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  !-------------------- H EXPANSION STRUCTURE ----------------------!
  type H_operator
     complex(8),dimension(:,:,:,:,:,:),allocatable   :: O  !Bath hamiltonian (replica/general)
  end type H_operator

  type(H_operator),dimension(:),allocatable             :: Hbath_basis  ![Nsym]
  real(8),dimension(:,:),allocatable                    :: Hbath_lambda ![Nbath,Nsym]
  logical                                               :: Hbath_status=.false.
  

  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type effective_bath_component
     integer                          :: N_dec
     real(8),dimension(:),allocatable :: v ![1] for "replica" bath; [Nbath] for "general" bath
     real(8),dimension(:),allocatable :: lambda
  end type effective_bath_component

  type effective_bath
     type(effective_bath_component),dimension(:),allocatable :: item ![Nbath]
     logical                                                 :: status=.false.
  end type effective_bath


  !-------------------- CUSTOM OBSERVABLE STRUCTURE ----------------------!
  type observable
    complex(8),dimension(:,:,:),allocatable :: sij
    character(len=32)                       :: o_name
    real(8)                                 :: o_value
  end type observable

  type custom_observables
     type(observable),dimension(:),allocatable               :: item
     complex(8),dimension(:,:,:),allocatable                 :: Hk
     integer                                                 :: N_asked
     integer                                                 :: N_filled
     logical                                                 :: init=.false.
  end type custom_observables


  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!
  type sector_map
     integer,dimension(:),allocatable :: map
     type(sparse_map)                 :: sp
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
     subroutine cc_sparse_HxV(Nloc,v,Hv)
       integer                    :: Nloc
       complex(8),dimension(Nloc) :: v
       complex(8),dimension(Nloc) :: Hv
     end subroutine cc_sparse_HxV
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
  !
  integer                                            :: Nimp     !Total number of levels in the impurity cluster: Nlat*Norb
  integer                                            :: Nlso     !Nlat*Nspin*Norb

  !local part of the Hamiltonian
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  complex(8),dimension(:,:,:,:,:,:),allocatable         :: impHloc           !local hamiltonian [Nlat][Nlat][Nspin][Nspin][Norb][Norb]


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
  procedure(cc_sparse_HxV),pointer                   :: spHtimesV_p=>null()


  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  integer,allocatable,dimension(:)                   :: neigen_sector
  !--------------- LATTICE WRAP VARIABLES -----------------! 
#if __GFORTRAN__ &&  __GNUC__ > 8      
  integer,allocatable,dimension(:,:)                 :: neigen_sector_ineq
  integer,allocatable,dimension(:)                   :: neigen_total_ineq
#endif
  logical                                            :: trim_state_list=.false.

  !Partition function
  !PRIVATE
  !=========================================================
  real(8)                                            :: zeta_function
  real(8)                                            :: gs_energy
  real(8)                                            :: max_exc



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
#if __GFORTRAN__ &&  __GNUC__ > 8     
   complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: Smatsii,Srealii                   ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
   complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: Gmatsii,Grealii                   ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
   complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: G0matsii,G0realii                 ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
   complex(8),dimension(:,:,:,:,:,:,:)  ,allocatable,save :: single_particle_density_matrix_ii ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb]
   complex(8),dimension(:,:,:)          ,allocatable,save :: cluster_density_matrix_ii         ![Nineq][4**(Nlat*Norb)][4**(Nlat*Norb)]
#endif


   !Spin Susceptibilities
  ! !=========================================================
      real(8),allocatable,dimension(:,:,:,:,:)               :: spinChi_tau
      complex(8),allocatable,dimension(:,:,:,:,:)            :: spinChi_w
      complex(8),allocatable,dimension(:,:,:,:,:)            :: spinChi_iv


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
  real(8),dimension(:,:),allocatable                   ::  ed_docc,ed_mag
  !--------------- LATTICE WRAP VARIABLES -----------------!
  type(custom_observables)                             ::  custom_o

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable                   :: wm,tau,wr,vm,vr




  !Local energies and generalized double occupancies
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  !real(8),dimension(:),allocatable                     :: ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot
  real(8)                                               :: ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot
  !real(8),dimension(:),allocatable                     :: ed_Dust,ed_Dund,ed_Dse,ed_Dph
  real(8)                                               :: ed_Dust,ed_Dund,ed_Dse,ed_Dph
  ! !--------------- LATTICE WRAP VARIABLES -----------------!
   real(8),dimension(:,:),allocatable,save            :: ddii,eii




  ! !Impurity operators
  ! !PRIVATE (now public but accessible thru routine)
  ! !=========================================================
   complex(8),allocatable,dimension(:,:,:,:,:,:)          :: single_particle_density_matrix ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
   complex(8),allocatable,dimension(:,:)                  :: cluster_density_matrix         ![4**(Nlat*Norb),4**(Nlat*Norb)]


#if __GFORTRAN__ &&  __GNUC__ > 8     
  !--------------- LATTICE WRAP VARIABLES -----------------!
  complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: Smats_ineq,Sreal_ineq    ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L]
  !complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: SAmats_ineq,SAreal_ineq
  complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: Gmats_ineq,Greal_ineq
  !complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: Fmats_ineq,Freal_ineq
  complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: G0mats_ineq,G0real_ineq
  !complex(8),dimension(:,:,:,:,:,:,:,:),allocatable,save :: F0mats_ineq,F0real_ineq
  !complex(8),dimension(:,:),allocatable,save         :: Dmats_ph_ineq,Dreal_ph_ineq
  complex(8),dimension(:,:,:,:,:,:,:),allocatable,save   :: single_particle_density_matrix_ineq ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  complex(8),dimension(:,:,:),allocatable,save       :: cluster_density_matrix_ineq ![Nineq,4**(Nlat*Norb),4**(Nlat*Norb)]
  real(8),dimension(:,:,:),allocatable,save          :: dens_ineq 
  real(8),dimension(:,:,:),allocatable,save          :: docc_ineq
  real(8),dimension(:,:,:),allocatable,save          :: mag_ineq
  !real(8),dimension(:,:,:,:),allocatable,save       :: phisc_ineq
  real(8),dimension(:,:),allocatable,save            :: dd_ineq,e_ineq
  real(8),dimension(:,:,:),allocatable               :: Hbath_lambda_ineq ![Nineq,Nbath,Nsym]
#endif


  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                  :: ed_file_suffix=""       !suffix string attached to the output files.
  character(len=10)                                  :: ineq_site_suffix="_ineq"
  integer                                            :: site_indx_padding=4
  logical                                            :: Jhflag              !spin-exchange and pair-hopping flag.
  logical                                            :: offdiag_gf_flag=.false.


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
  subroutine map_allocate_scalar(H,N,Nsp)
    type(sector_map) :: H
    integer          :: N
    integer,optional :: Nsp
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    if(present(Nsp))call sp_init_map(H%sp,Nsp)
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N,Nsp)
    type(sector_map),dimension(:) :: H
    integer,dimension(size(H))    :: N
    integer,optional              :: Nsp
    integer                       :: i
    do i=1,size(H)
       if(present(Nsp))then
          call map_allocate_scalar(H(i),N(i),Nsp)
       else
          call map_allocate_scalar(H(i),N(i))
       endif
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
    call sp_delete_map(H%sp)
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
