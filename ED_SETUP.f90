!|imp_up>|bath_up> * |imp_dw>|bath_dw>
!
!|imp_sigma>=
!| (1..Norb)_1 ... (1..Norb)_Nlat;         <-- Nlat*Norb == stride, indices here: imp_state_index
! [(1..Norb)_1 ... (1..Norb)_Nlat]_1       <-- Nlat*Norb + Nlat*Norb, "      "  : imp_state_index
!  ...
! [(1..Norb)_1 ... (1..Norb)_Nlat]_Nbath>; <-- Nbath*Norb*Nlat
!
!
MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only: free_unit,reg,file_length
  USE SF_MISC,    only: assert_shape
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  interface set_Himpurity
     module procedure :: set_Himpurity_lso_c
     module procedure :: set_Himpurity_nnn_c
  end interface set_Himpurity

  ! interface build_sector
  !    module procedure :: build_sector_full !A unique map for the full sector, separating spins
  !    module procedure :: build_sector_spin !Two spin-resolved maps, separating imp and bath states
  ! end interface build_sector

  ! interface delete_sector
  !    module procedure :: delete_sector_full
  !    module procedure :: delete_sector_spin
  ! end interface delete_sector

  public :: init_ed_structure
  public :: setup_global
  !
  public :: set_Himpurity
  !
  public :: build_sector
  public :: delete_sector
  !
  public :: get_Sector
  public :: get_Indices
  public :: get_Nup
  public :: get_Ndw
  public :: get_DimUp
  public :: get_DimDw
  !
  public :: indices2state
  public :: state2indices
  public :: iup_index
  public :: idw_index
  !
  public :: bdecomp
  public :: breorder
  public :: bjoin
  !
  public :: c,cdg
  !    
  public :: twin_sector_order
  public :: get_twin_sector
  public :: flip_state
  !
  public :: imp_state_index
  !
  public :: binary_search
  public :: sequential_search

#ifdef _MPI
  public :: scatter_vector_MPI
  public :: scatter_basis_MPI
  public :: gather_vector_MPI
  public :: allgather_vector_MPI
#endif



contains

  subroutine ed_checks_global
    !
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "ED ERROR: Nspin > 2 is currently not supported"
    if(Norb>5)stop "ED ERROR: Norb > 5 is currently not supported"
    !
    if(Nspin>1.AND.ed_twin.eqv..true.)then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
       call sleep(1)
    end if
    !
    if(lanc_method=="lanczos")then
       if(lanc_nstates_total>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_total==1, T=0"
       if(lanc_nstates_sector>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_sector==1, T=0"
    endif
    !
  end subroutine ed_checks_global


  !+------------------------------------------------------------------+
  !PURPOSE  : Setup Dimensions of the problem
  ! Norb    = # of impurity orbitals
  ! Nbath   = # of bath levels (depending on bath_type)
  ! Ns      = # of levels (per spin)
  ! Nlevels = 2*Ns = Total # of levels (counting spin degeneracy 2) 
  !+------------------------------------------------------------------+
  subroutine ed_setup_dimensions()
    !
    Nimp = Nlat*Norb    !Total number of levels in the impurity cluster
    Ns = Nimp*(Nbath+1) !Total number of levels per spin
    !
    Ns_Orb = Ns
    Ns_Ud  = 1
    !
    Nsectors = ((Ns_Orb+1)*(Ns_Orb+1))**Ns_Ud
  end subroutine ed_setup_dimensions



  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure()
    logical                          :: control
    integer                          :: i,iud,iorb,jorb,ispin,jspin
    integer,dimension(:),allocatable :: DimUps,DimDws
    !
    call ed_checks_global
    !
    call ed_setup_dimensions
    !
    !
    allocate(DimUps(Ns_Ud))
    allocate(DimDws(Ns_Ud))
    do iud=1,Ns_Ud
       DimUps(iud) = get_sector_dimension(Ns_Orb,Ns_Orb/2)
       DimDws(iud) = get_sector_dimension(Ns_Orb,Ns_Orb-Ns_Orb/2)
    enddo
    if(MpiMaster)then
       write(LOGfile,"(A)")"Summary:"
       write(LOGfile,"(A)")"--------------------------------------------"
       write(LOGfile,"(A,I15)")'# of levels/spin      = ',Ns
       write(LOGfile,"(A,I15)")'Total size            = ',2*Ns
       write(LOGfile,"(A,I15)")'# of sites            = ',Nlat
       write(LOGfile,"(A,I15)")'# of impurities       = ',Norb
       write(LOGfile,"(A,I15)")'# of bath/impurity    = ',Nbath
       write(LOGfile,"(A,I15)")'Number of sectors     = ',Nsectors
       write(LOGfile,"(A,"//str(Ns_Ud)//"I8,2X,"//str(Ns_Ud)//"I8,I20)")&
            'Largest Sector(s)     = ',DimUps,DimDws,product(DimUps)*product(DimDws)
       write(LOGfile,"(A)")"--------------------------------------------"
    endif
    call sleep(1)
    !
    !
    allocate(spH0ups(Ns_Ud))
    allocate(spH0dws(Ns_Ud))
    !
    !Allocate indexing arrays
    allocate(getCsector(Ns_Ud,2,Nsectors))  ;getCsector  =0
    allocate(getCDGsector(Ns_Ud,2,Nsectors));getCDGsector=0
    !
    allocate(getDim(Nsectors));getDim=0
    !
    allocate(getBathStride(Nlat,Norb,Nbath));getBathStride=0
    allocate(twin_mask(Nsectors))
    allocate(sectors_mask(Nsectors))
    allocate(neigen_sector(Nsectors))
    !
    !check finiteT
    finiteT=.true.                    !assume doing finite T per default
    if(lanc_nstates_total==1)then     !is you only want to keep 1 state
       finiteT=.false.                !set to do zero temperature calculations
       write(LOGfile,"(A)")"Required Lanc_nstates_total=1 => set T=0 calculation"
    endif
    !
    !check whether lanc_nstates_sector and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_nstates_sector,2)/=0)then
          lanc_nstates_sector=lanc_nstates_sector+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
       endif
       if(mod(lanc_nstates_total,2)/=0)then
          lanc_nstates_total=lanc_nstates_total+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
       endif
       write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
       !
       write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
       write(LOGfile,"(A,I3)")"Nstates   Total  = ", lanc_nstates_total
       call sleep(1)
    else
       write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
       call sleep(1)
    endif
    !
    Jhflag=.FALSE.
    if(Norb>1.AND.(Jx/=0d0.OR.Jp/=0d0))Jhflag=.TRUE.
    !
    !
    offdiag_gf_flag=.true.
    !
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
    endif

    !ALLOCATE impHloc
    if(.not.allocated(impHloc))then
       allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
       impHloc=zero
    else
       call assert_shape(impHloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"init_ed_structure","impHloc")
    endif
    !
    !
    !allocate functions
    !allocate functions
    allocate(impSmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impSmats=zero
    impSreal=zero
    !
    allocate(impGmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impGmats=zero
    impGreal=zero
    !
    !allocate functions
    allocate(impG0mats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impG0real(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impG0mats=zero
    impG0real=zero
    !

    !
    !allocate observables
    allocate(ed_dens(Nlat,Norb),ed_docc(Nlat,Norb),ed_mag(Nlat,Norb),ed_dens_up(Nlat,Norb),ed_dens_dw(Nlat,Norb))
    ed_dens=0d0
    ed_docc=0d0
    ed_mag=0d0
    ed_dens_up=0d0
    ed_dens_dw=0d0
    !
     if(chiflag)then
      allocate(spinChi_tau(Nlat,Nlat,Norb,Norb,0:Ltau))
      allocate(spinChi_w(Nlat,Nlat,Norb,Norb,Lreal))
      allocate(spinChi_iv(Nlat,Nlat,Norb,Norb,0:Lmats))
    !    !
    !    ! allocate(densChi_tau(Norb,Norb,0:Ltau))
    !    ! allocate(densChi_w(Norb,Norb,Lreal))
    !    ! allocate(densChi_iv(Norb,Norb,0:Lmats))
    !    ! allocate(densChi_mix_tau(Norb,Norb,0:Ltau))
    !    ! allocate(densChi_mix_w(Norb,Norb,Lreal))
    !    ! allocate(densChi_mix_iv(Norb,Norb,0:Lmats))
    !    ! allocate(densChi_tot_tau(0:Ltau))
    !    ! allocate(densChi_tot_w(Lreal))
    !    ! allocate(densChi_tot_iv(0:Lmats))
    !    ! !
    !    ! allocate(pairChi_tau(Norb,0:Ltau))
    !    ! allocate(pairChi_w(Norb,Lreal))
    !    ! allocate(pairChi_iv(Norb,0:Lmats))
     endif
    !
    !allocate density matrices
    !
    allocate(single_particle_density_matrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    single_particle_density_matrix=zero
    !
    allocate(cluster_density_matrix(4**Nimp,4**Nimp))
    cluster_density_matrix=zero
    !
  end subroutine init_ed_structure


  !+------------------------------------------------------------------+
  !PURPOSE  : Setup Himpurity, the local part of the non-interacting Hamiltonian
  !+------------------------------------------------------------------+
  subroutine set_Himpurity_nnn_c(hloc)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hloc
    if(allocated(impHloc))deallocate(impHloc)
    allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb));impHloc=zero
    impHloc = Hloc
    if(ed_verbose>2)call print_hloc(impHloc)
  end subroutine set_Himpurity_nnn_c

  subroutine set_Himpurity_lso_c(hloc)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: hloc
    if(allocated(impHloc))deallocate(impHloc)
    allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb));impHloc=zero
    impHloc = lso2nnn_reshape(Hloc,Nlat,Nspin,Norb)
    if(ed_verbose>2)call print_hloc(impHloc)
  end subroutine set_Himpurity_lso_c




  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  subroutine setup_global
    integer                          :: DimUp,DimDw
    integer                          :: DimUps(Ns_Ud),DimDws(Ns_Ud)
    integer                          :: Indices(2*Ns_Ud),Jndices(2*Ns_Ud)
    integer                          :: Nups(Ns_ud),Ndws(Ns_ud)
    integer                          :: Jups(Ns_ud),Jdws(Ns_ud)
    integer                          :: i,iud,iorb,ilat,stride,ibath
    integer                          :: isector,jsector
    integer                          :: unit,status,istate,ishift,isign
    logical                          :: IOfile
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    !
    !Store full dimension of the sectors:
    do isector=1,Nsectors
       call get_DimUp(isector,DimUps)
       call get_DimDw(isector,DimDws)
       DimUp = product(DimUps)
       DimDw = product(DimDws)       
       getDim(isector)  = DimUp*DimDw
    enddo
    !
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Restarting from a state_list file:"
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       !
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       status=0
       do while(status>=0)
          read(unit,*,iostat=status)istate,isector,indices
          list_sector(istate)=isector
          call get_Nup(isector,Nups)
          call get_Ndw(isector,Ndws)
          if(any(Indices /= [Nups,Ndws]))&
               stop "setup_global error: nups!=nups(isector).OR.ndws!=ndws(isector)"
       enddo
       close(unit)
       !
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getDim(isector),lanc_nstates_sector) !init every sector to required eigenstates
       enddo
    endif
    !
    twin_mask=.true.
    if(ed_twin)then
       do isector=1,Nsectors
          call get_Nup(isector,Nups)
          call get_Ndw(isector,Ndws)
          if(any(Nups .ne. Ndws))then
            call get_Sector([Ndws,Nups],Ns_Orb,jsector)
            if (twin_mask(jsector))twin_mask(isector)=.false.
          endif
       enddo
       write(LOGfile,"(A,I6,A,I9)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
       call sleep(1)
    endif
    !
    stride=Nlat*Norb
    do ibath=1,Nbath
       do ilat=1,Nlat
          do iorb=1,Norb
             getBathStride(ilat,iorb,ibath) = stride + imp_state_index(ilat,iorb) + (ibath-1)*Nlat*Norb
             ! getBathStride(ilat,iorb,ibath) = stride + iorb + (ilat-1)*Norb + (ibath-1)*Nlat*Norb
          enddo
       enddo
    enddo
    !
    getCsector  = 0
    getCDGsector= 0
    do isector=1,Nsectors
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       !
       !UPs:
       !DEL:
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jups(iud)=Jups(iud)-1; if(Jups(iud) < 0)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,1,isector)=jsector
       enddo
       !ADD
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws
          Jups(iud)=Jups(iud)+1; if(Jups(iud) > Ns_Orb)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,1,isector)=jsector
       enddo
       !
       !DWs:
       !DEL
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)-1; if(Jdws(iud) < 0)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,2,isector)=jsector
       enddo
       !DEL
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)+1; if(Jdws(iud) > Ns_Orb)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,2,isector)=jsector
       enddo
    enddo
    return
  end subroutine setup_global













  !##################################################################
  !##################################################################
  !AUXILIARY PROCEDURES - Sectors,Nup,Ndw,DimUp,DimDw,...
  !##################################################################
  !##################################################################
  elemental function get_sector_dimension(n,np) result(dim)
    integer,intent(in) :: n,np
    integer            :: dim
    dim = binomial(n,np)
  end function get_sector_dimension


  subroutine get_Sector(indices,N,isector)
    integer,dimension(:) :: indices
    integer              :: N
    integer              :: isector
    integer              :: i,Nind,factor
    Nind = size(indices)
    Factor = N+1
    isector = 1
    do i=Nind,1,-1
       isector = isector + indices(i)*(Factor)**(Nind-i)
    enddo
  end subroutine get_Sector


  subroutine get_Indices(isector,N,indices)
    integer                          :: isector,N
    integer,dimension(:)             :: indices
    integer                          :: i,count,Dim
    integer,dimension(size(indices)) :: indices_
    !
    Dim = size(indices)
    if(mod(Dim,2)/=0)stop "get_Indices_main error: Dim%2 != 0"
    count=isector-1
    do i=1,Dim
       indices_(i) = mod(count,N+1)
       count      = count/(N+1)
    enddo
    indices = indices_(Dim:1:-1)
  end subroutine get_Indices


  subroutine get_Nup(isector,Nup)
    integer                   :: isector,Nup(Ns_Ud)
    integer                   :: i,count
    integer,dimension(2*Ns_Ud)  :: indices_
    count=isector-1
    do i=1,2*Ns_Ud
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Nup = indices_(2*Ns_Ud:Ns_Ud+1:-1)
  end subroutine get_Nup


  subroutine get_Ndw(isector,Ndw)
    integer                   :: isector,Ndw(Ns_Ud)
    integer                   :: i,count
    integer,dimension(2*Ns_Ud) :: indices_
    count=isector-1
    do i=1,2*Ns_Ud
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Ndw = indices_(Ns_Ud:1:-1)
  end subroutine get_Ndw


  subroutine  get_DimUp(isector,DimUps)
    integer                :: isector,DimUps(Ns_Ud)
    integer                :: Nups(Ns_Ud),iud
    call get_Nup(isector,Nups)
    do iud=1,Ns_Ud
       DimUps(iud) = binomial(Ns_Orb,Nups(iud))
    enddo
  end subroutine get_DimUp


  subroutine get_DimDw(isector,DimDws)
    integer                :: isector,DimDws(Ns_Ud)
    integer                :: Ndws(Ns_Ud),iud
    call get_Ndw(isector,Ndws)
    do iud=1,Ns_Ud
       DimDws(iud) = binomial(Ns_Orb,Ndws(iud))
    enddo
  end subroutine get_DimDw


  subroutine indices2state(ivec,Nvec,istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end subroutine indices2state

  subroutine state2indices(istate,Nvec,ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end subroutine state2indices


  function iup_index(i,DimUp) result(iup)
    integer :: i
    integer :: DimUp
    integer :: iup
    iup = mod(i,DimUp);if(iup==0)iup=DimUp
  end function iup_index


  function idw_index(i,DimUp) result(idw)
    integer :: i
    integer :: DimUp
    integer :: idw
    idw = (i-1)/DimUp+1
  end function idw_index

  !> Find position in the state vector for a given lattice-spin-orbital position for the cluster (no bath considered)
  function imp_state_index(ilat,iorb) result(indx)  
    integer :: ilat
    integer :: iorb
    integer :: indx
    indx = iorb + (ilat-1)*Norb
  end function imp_state_index




#ifdef _MPI
  !! Scatter V into the arrays Vloc on each thread: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine scatter_vector_MPI(MpiComm,v,vloc)
    integer                          :: MpiComm
    complex(8),dimension(:)          :: v    !size[N]
    complex(8),dimension(:)          :: vloc !size[Nloc]
    integer                          :: i,irank,Nloc,N
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
    if( MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null )return
    ! stop "scatter_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "scatter_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    Vloc=0
    call MPI_Scatterv(V,Counts,Offset,MPI_DOUBLE_COMPLEX,Vloc,Nloc,MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
    !
    return
  end subroutine scatter_vector_MPI


  subroutine scatter_basis_MPI(MpiComm,v,vloc)
    integer                   :: MpiComm
    complex(8),dimension(:,:) :: v    !size[N,N]
    complex(8),dimension(:,:) :: vloc !size[Nloc,Neigen]
    integer                   :: N,Nloc,Neigen,i
    N      = size(v,1)
    Nloc   = size(vloc,1)
    Neigen = size(vloc,2)
    if( size(v,2) < Neigen ) stop "error scatter_basis_MPI: size(v,2) < Neigen"
    !
    do i=1,Neigen
       call scatter_vector_MPI(MpiComm,v(:,i),vloc(:,i))
    end do
    !
    return
  end subroutine scatter_basis_MPI


  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine gather_vector_MPI(MpiComm,vloc,v)
    integer                             :: MpiComm
    complex(8),dimension(:)             :: vloc !size[Nloc]
    complex(8),dimension(:)             :: v    !size[N]
    integer                             :: i,irank,Nloc,N
    integer,dimension(:),allocatable    :: Counts,Offset
    integer                             :: MpiSize,MpiIerr
    logical                             :: MpiMaster
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    !stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    call MPI_Gatherv(Vloc,Nloc,MPI_DOUBLE_COMPLEX,V,Counts,Offset,MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
    !
    return
  end subroutine gather_vector_MPI


  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine allgather_vector_MPI(MpiComm,vloc,v)
    integer                             :: MpiComm
    complex(8),dimension(:)             :: vloc !size[Nloc]
    complex(8),dimension(:)             :: v    !size[N]
    integer                             :: i,irank,Nloc,N
    integer,dimension(:),allocatable    :: Counts,Offset
    integer                             :: MpiSize,MpiIerr
    logical                             :: MpiMaster
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    ! stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N    = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "allgather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    V = 0d0
    call MPI_AllGatherv(Vloc,Nloc,MPI_DOUBLE_COMPLEX,V,Counts,Offset,MPI_DOUBLE_COMPLEX,MpiComm,MpiIerr)
    !
    return
  end subroutine Allgather_vector_MPI
#endif




  !######################################################################
  !######################################################################
  !BUILD SECTORS: given isector build / delete the corresponding map(s)
  !######################################################################
  !######################################################################
  !
  subroutine build_sector(isector,H,Itrace)
    integer                             :: isector
    type(sector_map),dimension(2*Ns_Ud) :: H
    logical,optional                    :: itrace
    logical                             :: itrace_
    integer,dimension(Ns_Ud)            :: Nups,Ndws
    integer,dimension(Ns_Ud)            :: DimUps,DimDws
    integer                             :: DimUp,DimDw,impDIM
    integer                             :: iup,idw
    integer                             :: nup_,ndw_
    integer                             :: imap,iud
    integer                             :: iIMP,iBATH
    !
    itrace_=.false. ; if(present(itrace))itrace_=itrace
    !
    impDIM = 2**(Nimp/Ns_ud)
    !
    call get_Nup(isector,Nups)
    call get_Ndw(isector,Ndws)
    call get_DimUp(isector,DimUps)
    call get_DimDw(isector,DimDws)
    !
    if(itrace_)then
       call map_allocate(H,[DimUps,DimDws],impDIM)
    else
       call map_allocate(H,[DimUps,DimDws])
    endif
    !
    do iud=1,Ns_Ud
       !UP    
       imap=0
       do iup=0,2**Ns_Orb-1
          nup_ = popcnt(iup)  
          if(nup_ /= Nups(iud))cycle
          imap = imap+1
          H(iud)%map(imap) = iup
          if(.not.itrace_)cycle
          iIMP  = ibits(iup,0,Nimp)
          iBATH = ibits(iup,Nimp,Nimp*Nbath)
          call sp_insert_state(H(iud)%sp,iIMP,iBATH,imap) 
       enddo
       !DW
       imap=0
       do idw=0,2**Ns_Orb-1
          ndw_= popcnt(idw)
          if(ndw_ /= Ndws(iud))cycle
          imap = imap+1
          H(iud+Ns_Ud)%map(imap) = idw
          if(.not.itrace_)cycle
          iIMP  = ibits(idw,0,Nimp)
          iBATH = ibits(idw,Nimp,Nimp*Nbath)
          call sp_insert_state(H(iud+Ns_Ud)%sp,iIMP,iBATH,imap) 
       enddo
    enddo
    !
  end subroutine build_sector


  subroutine delete_sector(isector,H)
    integer                   :: isector
    type(sector_map)          :: H(:)
    call map_deallocate(H)
  end subroutine delete_sector

  ! subroutine delete_sector_spin(isector,HUP,HDW)
  !   integer                   :: isector
  !   type(sector_map)          :: HUP
  !   type(sector_map)          :: HDW
  !   call map_deallocate(HUP)
  !   call map_deallocate(HDW)
  ! end subroutine delete_sector_spin






  !##################################################################
  !##################################################################
  !CREATION / DESTRUCTION OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg






  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw.
  !
  !- build the map from the A-sector to \HHH
  !- get the list of states in \HHH corresponding to sector B twin of A
  !- return the ordering of B-states in \HHH with respect to those of A
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                             :: isector
    integer,dimension(:)                :: order
    type(sector_map),dimension(2*Ns_Ud) :: H
    integer,dimension(2*Ns_Ud)          :: Indices,Istates
    integer,dimension(Ns_Ud)            :: DimUps,DimDws
    integer                             :: Dim
    integer                             :: i,iud
    !
    Dim = GetDim(isector)
    if(size(Order)/=Dim)stop "twin_sector_order error: wrong dimensions of *order* array"
    call get_DimUp(isector,DimUps)
    call get_DimDw(isector,DimDws)
    !
    call build_sector(isector,H)
    do i=1,Dim
       call state2indices(i,[DimUps,DimDws],Indices)
       forall(iud=1:2*Ns_Ud)Istates(iud) = H(iud)%map(Indices(iud))
       Order(i) = flip_state( Istates )
    enddo
    call delete_sector(isector,H)
    !
    call sort_array(Order)
    !
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  !+------------------------------------------------------------------+
  function flip_state(istate) result(j)
    integer,dimension(2*Ns_Ud) :: istate
    integer                    :: j
    integer,dimension(Ns_Ud)   :: jups,jdws
    integer,dimension(2*Ns_Ud) :: dims
    !
    jups = istate(Ns_Ud+1:2*Ns_Ud)
    jdws = istate(1:Ns_Ud)
    dims = 2**Ns_Orb
    call indices2state([jups,jdws],Dims,j)
    !
  end function flip_state


  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite 
  ! quantum numbers): 
  ! nup,ndw ==> ndw,nup (spin-exchange)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in)       :: isector
    integer                  :: jsector
    integer,dimension(Ns_Ud) :: Iups,Idws
    call get_Nup(isector,iups)
    call get_Ndw(isector,idws)
    call get_Sector([idws,iups],Ns_Orb,jsector)
  end function get_twin_sector










  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp


  !+------------------------------------------------------------------+
  ! Reorder a binary decomposition so to have a state of the form:
  !
  ! BATH_TYPE: |binary layout of the spin states>_spin
  !
  ! normal:  |(1:Norb),([1:Nbath]_1, [1:Nbath]_2, ... ,[1:Nbath]_Norb)>_spin !> NOT IMPLEMENTED HERE <!
  ! hybrid:  |(1:Norb),([1:Nbath])>_spin                                     !> NOT IMPLEMENTED HERE <!
  ! replica: |(1:Norb),([1:Norb]_1, [1:Norb]_2, ...  , [1:Norb]_Nbath)>_spin
  ! general: |(1:Norb),([1:Norb]_1, [1:Norb]_2, ...  , [1:Norb]_Nbath)>_spin
  !
  !> case (ed_total_ud):
  !   (T): Ns_Ud=1, Ns_Orb=Ns.
  !        bdecomp is already of the form above [1:Ns]
  !   (F): Ns_Ud=Norb, Ns_Orb=Ns/Norb==1+Nbath
  !        bdecomp is
  !        |( [1:1+Nbath]_1,...,[1:1+Nbath]_Norb)>_spin
  !+------------------------------------------------------------------+
  function breorder(Nins) result(Ivec)
    integer,intent(in),dimension(Ns_Ud,Ns_Orb) :: Nins ![1,Ns] - [Norb,1+Nbath]
    integer,dimension(Ns)                      :: Ivec ![Ns]
    integer                                    :: iud,ibath,indx
    ! select case (ed_total_ud) !!! this is true by default in the CDMFT code
    ! case (.true.)
    Ivec = Nins(1,:)
    ! case (.false.)
    !    do iud=1,Ns_Ud           ![1:Norb]
    !       Ivec(iud) = Nins(iud,1)
    !       do ibath=1,Nbath
    !          indx = getBathStride(iud,ibath) !take care of normal/
    !          Ivec(indx) = Nins(iud,1+ibath)
    !       enddo
    !    enddo
    ! end select
  end function breorder


  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a integer value in a *sorted* array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !write(*,*)  "WARNING: binary_search: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search




  !+------------------------------------------------------------------+
  !PURPOSE : sequential search of a integer value in a generic array
  !+------------------------------------------------------------------+
  function sequential_search(a,value) result(ssresult)
    integer,intent(in) :: a(:), value
    integer            :: ssresult,i
    do i=1,size(a)
       !print*,"i="//str(i)
       if(a(i)==value)then
          !print*,"here we set to i and exit"
          ssresult=i
          exit
       else
          ! print*,"here we set to zero and loop"
          ssresult=0
       endif
    enddo
    !print*,"result="//str(ssresult)
  end function sequential_search




  !+------------------------------------------------------------------+
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array




end MODULE ED_SETUP













!   ! SPIN-RESOLVED ALGORITHM: two map objects Hsigma, with a *sparse* 
!   !                          structure capable to store separately
!   !                          the impurity and the bath spin-states
!   subroutine build_sector_spin(isector,Hup,Hdw)
!     integer                            :: isector
!     type(sector_map)                   :: HUP !Map for the UPs
!     type(sector_map)                   :: HDW !Map for the Dws
!     integer,dimension(Ns_Ud)           :: Nups,Ndws
!     integer,dimension(Ns_Ud)           :: DimUps,DimDws
!     integer                            :: DimUp,DimDw,impDIM
!     integer                            :: iup,idw
!     integer                            :: nup_,ndw_
!     integer                            :: imap,iud
!     integer                            :: iIMP,iBATH
!     !
!     impDIM = 2**(Nimp/Ns_ud) !Number of states for the impurity
!     !
!     !Init UP sub-sector:
!     ! > allocate UP-map to binomial(Ns_Orb Nup)
!     call get_Nup(isector,Nups)
!     call get_DimUp(isector,DimUps); DimUp = product(DimUps)
!     call map_allocate(HUP,DimUp,impDIM)
!     !
!     !Init DW sub-sector:
!     ! > allocate DW-map to binomial(Ns_Orb Ndw)
!     call get_Ndw(isector,Ndws)
!     call get_DimDw(isector,DimDws); DimDw = product(DimDws)
!     call map_allocate(HDW,DimDw,impDIM)
!     !
!     !Formally dealing with ed_total_ud == .false.
!     !->iud=1:Ns_Ud (formally because Ns_Ud = 1 in CDMFT code...)
!     ![STILL PROBLEMS WITH Ns_Ud, don't know how to deal with it in density_matrix_impurity()] 
!     do iud=1,Ns_Ud 
! #ifdef _DEBUG
!        if(ed_verbose>3)write(Logfile,"(A)")&
!             "  DEBUG build_sector_spin(): working UP sub-sector"
! #endif
!        imap=0
!        do iup=0,2**Ns_Orb-1
!           nup_ = popcnt(iup) !equivalent to sum(binary_decomposition(iup))
!           if(nup_ /= Nups(iud))cycle !the state does not have the required number of UPs
!           imap = imap+1
!           !HUP(iud)%map(imap) = iup
!           HUP%map(imap) = iup
!           !
!           iIMP  = ibits(iup,0,Nimp)
!           iBATH = ibits(iup,Nimp,Nimp*Nbath) !check: Nimp+Nimp*Nbath=Nimp(1*Nbath)=Ns
!           !call sp_insert_state(HUP(iud)%sp,iIMP,iBATH,imap)
! #ifdef _DEBUG
!           if(ed_verbose>4)then 
!              write(Logfile,"(A)")&
!                   "    DEBUG build_sector_spin(): inserting UP state"
!              write(Logfile,"(A)")&
!                   "      Iup: "//str(iup)//" | IimpUp: "//str(iIMP)//" | IbathUp: "//str(iBATH)
!           endif
! #endif
!           call sp_insert_state(HUP%sp,iIMP,iBATH,imap) 
!           !
!        enddo
!        !
! #ifdef _DEBUG
!        if(ed_verbose>3)write(Logfile,"(A)")&
!             "  DEBUG build_sector_spin(): working DW sub-sector"
! #endif
!        imap=0
!        do idw=0,2**Ns_Orb-1
!           ndw_=popcnt(idw)   !equivalent to sum(binary_decomposition(idw))
!           if(ndw_ /= Ndws(iud))cycle !the state does not have the required number of DWs
!           imap = imap+1
!           !HDW(iud)%map(imap) = idw
!           HDW%map(imap) = idw
!           !
!           iIMP  = ibits(idw,0,Nimp)
!           iBATH = ibits(idw,Nimp,Nimp*Nbath) !check: Nimp+Nimp*Nbath=Nimp(1*Nbath)=Ns
!           !call sp_insert_state(HDW(iud)%sp,iIMP,iBATH,imap)
! #ifdef _DEBUG
!           if(ed_verbose>4)then 
!              write(LOGfile,"(A)")&
!                   "    DEBUG build_sector_spin(): inserting DW state"
!              write(LOGfile,"(A)")&
!                   "      Idw: "//str(idw)//" | IimpDw: "//str(iIMP)//" | IbathDw: "//str(iBATH)
!           endif
! #endif
!           call sp_insert_state(HDW%sp,iIMP,iBATH,imap)
!           !
!        enddo
!     enddo
!     !
!   end subroutine build_sector_spin




! getCsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup-1; jdw=ndw; if(jup < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw-1;if(jdw < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,2,isector)=jsector
! enddo
! !
! !
! !
! getCDGsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup+1;jdw=ndw;if(jup > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,2,isector)=jsector
! enddo
