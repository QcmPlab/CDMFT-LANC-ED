module ED_MAIN
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace,delete_eigenspace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG
  USE SF_IOTOOLS, only: str,reg
  USE SF_TIMER,only: start_timer,stop_timer
  implicit none
  private



  interface ed_init_solver
     module procedure :: ed_init_solver_single
#ifdef _MPI
     module procedure :: ed_init_solver_single_mpi
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_init_solver_lattice_mpi   
#endif  
#endif
  end interface ed_init_solver

  interface ed_solve
     module procedure :: ed_solve_single
#ifdef _MPI
     module procedure :: ed_solve_single_mpi
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_solve_lattice_mpi
#endif
#endif
  end interface ed_solve


  public :: ed_init_solver
  public :: ed_solve




contains





  ! PURPOSE: allocate and initialize one or multiple baths -+!
  subroutine ed_init_solver_single(bath)
    real(8),dimension(:),intent(inout) :: bath
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    !Init bath:
    !call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath()
    !if( (Nspin>1) .AND. &
    !     any(Hloc(:,:,1,Nspin,:,:).ne.0d0) )stop "ED ERROR: impHloc.mask(s,s`) /= 0. Spin-Flip terms are not allowed"
    call init_dmft_bath()
    call get_dmft_bath(bath)    !dmft_bath --> user_bath
    !
    if(isetup)call setup_global
    call deallocate_dmft_bath()
    isetup=.false.
    !
  end subroutine ed_init_solver_single

  
  
#ifdef _MPI
  subroutine ed_init_solver_single_mpi(MpiComm,bath)
    integer                            :: MpiComm
    real(8),dimension(:),intent(inout) :: bath
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    !
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm(MpiComm)
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath()
    call init_dmft_bath()
    call get_dmft_bath(bath)
    if(isetup)call setup_global
    call deallocate_dmft_bath()
    isetup=.false.
    !
    call ed_del_MpiComm()
    !
  end subroutine ed_init_solver_single_mpi
#endif


#ifdef _MPI
#if __GFORTRAN__ &&  __GNUC__ > 8     
  subroutine ed_init_solver_lattice_mpi(MpiComm,bath)
    integer                        :: MpiComm
    real(8),dimension(:,:)         :: bath ![Nineq][:]
    integer                        :: iineq,Nineq
    logical                        :: check
    integer                        :: MPI_ERR
    !
    !
    Nineq = size(bath,1)
    if(.not.allocated(Hreplica_lambda_ineq))stop "ERROR ed_init_solver: replica parameters lambda not defined for all sites"
    !
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(dd_ineq))deallocate(dd_ineq)
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    if(allocated(Gmats_ineq))deallocate(Gmats_ineq)
    if(allocated(Greal_ineq))deallocate(Greal_ineq)
    if(allocated(single_particle_density_matrix_ineq))deallocate(single_particle_density_matrix_ineq)
    if(allocated(cluster_density_matrix_ineq))deallocate(cluster_density_matrix_ineq)

    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
    !
    allocate(dens_ineq(Nineq,Nlat,Norb))
    allocate(docc_ineq(Nineq,Nlat,Norb))
    allocate(mag_ineq(Nineq,3,Norb))
    allocate(e_ineq(Nineq,4))
    allocate(dd_ineq(Nineq,4))
    !
    allocate(Smats_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sreal_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    !
    allocate(Gmats_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Greal_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    !
    allocate(single_particle_density_matrix_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(cluster_density_matrix_ineq(Nineq,4**(Nlat*Norb),4**(Nlat*Norb)))

    !
    do iineq=1,Nineq             !all nodes check the bath, u never know...
      !
      ed_file_suffix=reg(ineq_site_suffix)//str(iineq,site_indx_padding)
      call Hreplica_site(iineq)
      call ed_init_solver_single(bath(iineq,:)) !Here we init the solver using impHloc, whatever its origin.
      !
   end do
   !
   call MPI_Barrier(MpiComm,MPI_ERR)
   !
   ed_file_suffix=""
   !
   allocate(neigen_sector_ineq(Nineq,Nsectors))
   allocate(neigen_total_ineq(Nineq))
   do iineq=1,Nineq       
     neigen_sector_ineq(iineq,:) = neigen_sector(:)
     neigen_total_ineq(iineq)    = lanc_nstates_total
   end do
   !
  end subroutine ed_init_solver_lattice_mpi
#endif
#endif

  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single(bath,Hloc)
    real(8),dimension(:),intent(in) :: bath
    complex(8),intent(in)           :: Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    logical                         :: check
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    call set_Himpurity(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE error: wrong bath dimensions"
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath)    !user_bath --> dmft_bath
    call write_dmft_bath(LOGfile)
    if(MpiMaster)call save_dmft_bath(used=.true.)
    !
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()                   !find target states by digonalization of Hamiltonian
    call observables_impurity()                   !obtain impurity observables as thermal averages.
    call get_custom_observables()                 !obtain custom user-defined observables (if initialized)
    call local_energy_impurity()                  !obtain the local energy of the effective impurity problem
    !
    !OPTIONAL (HEAVY) CALCULATIONS 
    if(ed_get_dm) call density_matrix_impurity()  !build the cluster density matrix (\rho_IMP = Tr_BATH(\rho))
    !if(chiflag) call buildchi_impurity()         !build the local susceptibilities (todo)
    !
    !GET IMPURITY GFs and RELATED QUANTITIES
    call buildgf_impurity()                       !build the one-particle impurity Green's functions & Self-energy
    !
    call deallocate_dmft_bath()
    call es_delete_espace(state_list)
    !
    nullify(spHtimesV_p)
  end subroutine ed_solve_single
#ifdef _MPI
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single_mpi(MpiComm,bath,Hloc)
    integer                         :: MpiComm
    real(8),dimension(:),intent(in) :: bath
    complex(8),intent(in)           :: Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    logical                         :: check
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm(MpiComm)
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    call set_Himpurity(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE error: wrong bath dimensions"
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath)    !user_bath --> dmft_bath
    call write_dmft_bath(LOGfile)
    if(MpiMaster)call save_dmft_bath(used=.true.)
    !
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()                   !find target states by digonalization of Hamiltonian
    call observables_impurity()                   !obtain impurity observables as thermal averages.
    call get_custom_observables()                 !obtain custom user-defined observables(if initialized)
    call local_energy_impurity()                  !obtain the local energy of the effective impurity problem
    !
    !OPTIONAL (HEAVY) CALCULATIONS 
    if(ed_get_dm) call density_matrix_impurity()  !build the cluster density matrix (\rho_IMP = Tr_BATH(\rho))
    !if(chiflag) call buildchi_impurity()         !build the local susceptibilities (todo)
    !
    !GET IMPURITY GFs and RELATED QUANTITIES
    call buildgf_impurity()                       !build the one-particle impurity Green's functions  & Self-energy
    !
    call deallocate_dmft_bath()
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    !
    nullify(spHtimesV_p)
  end subroutine ed_solve_single_mpi
#endif

  !FALL BACK: DO A VERSION THAT DOES THE SITES IN PARALLEL USING SERIAL ED CODE
#ifdef _MPI
#if __GFORTRAN__ &&  __GNUC__ > 8     
  subroutine ed_solve_lattice_mpi(MpiComm,bath,Hloc,mpi_lanc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii)
    integer          :: MpiComm
    !inputs
    real(8)          :: bath(:,:) ![Nineq][Nb]
    complex(8)       :: Hloc(size(bath,1),Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    logical,optional :: mpi_lanc
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    real(8),optional :: Jp_ii(size(bath,1))
    real(8),optional :: Jx_ii(size(bath,1))
    !MPI  auxiliary vars
    ! 
    integer          :: i,j,iineq,iorb,jorb,ispin,jspin,ilat,jlat
    integer          :: Nineq
    logical          :: check_dim,mpi_lanc_
    !
    integer          :: MPI_ID=0
    integer          :: MPI_SIZE=1
    logical          :: MPI_MASTER=.true.
    !
    integer          :: mpi_err 
    !
    MPI_ID     = get_Rank_MPI(MpiComm)
    MPI_SIZE   = get_Size_MPI(MpiComm)
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    mpi_lanc_=.false.;if(present(mpi_lanc))mpi_lanc_=mpi_lanc
    !
    ! Check dimensions !
    Nineq=size(bath,1)
    !
    if(size(neigen_sector_ineq,1)<Nineq)stop "ed_solve_lattice error: size(neigen_sectorii,1)<Nineq"
    if(size(neigen_total_ineq)<Nineq)stop "ed_solve_lattice error: size(neigen_totalii,1)<Nineq"
    !
    !Check the dimensions of the bath are ok.
    !This can always be done in parallel no issues with mpi_lanc
    do iineq=1+MPI_ID,Nineq,MPI_SIZE
       check_dim = check_bath_dimension(bath(iineq,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    Smats_ineq    = zero ; Sreal_ineq    = zero
    Gmats_ineq    = zero ; Greal_ineq    = zero
    dens_ineq     = 0d0  ; docc_ineq     = 0d0
    mag_ineq      = 0d0
    e_ineq        = 0d0  ; dd_ineq       = 0d0 
    single_particle_density_matrix_ineq = zero
    cluster_density_matrix_ineq = zero
    !
    !solve sites serial, Lanczos with MPI
    if(MPI_MASTER)call start_timer
    do iineq = 1, Nineq
       write(LOGfile,*)" SOLVING INEQ SITE: "//str(iineq,Npad=4)
       ed_file_suffix=reg(ineq_site_suffix)//str(iineq,site_indx_padding)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(iineq,1:Norb)
       if(present(Ust_ii)) Ust = Ust_ii(iineq)
       if(present(Jh_ii))  Jh  = Jh_ii(iineq)
       if(present(Jp_ii))  Jp  = Jp_ii(iineq)
       if(present(Jx_ii))  Jx  = Jx_ii(iineq)
       !
       !Solve the impurity problem for the iineq-th site
       neigen_sector(:)   = neigen_sector_ineq(iineq,:)
       lanc_nstates_total = neigen_total_ineq(iineq)
       !
       call ed_solve_single_mpi(MpiComm,bath(iineq,:),Hloc(iineq,:,:,:,:,:,:))
       !
       neigen_sector_ineq(iineq,:)  = neigen_sector(:)
       neigen_total_ineq(iineq)     = lanc_nstates_total
       Smats_ineq(iineq,:,:,:,:,:,:,:)  = impSmats(:,:,:,:,:,:,:)
       Sreal_ineq(iineq,:,:,:,:,:,:,:)  = impSreal(:,:,:,:,:,:,:)
       Gmats_ineq(iineq,:,:,:,:,:,:,:)  = impGmats(:,:,:,:,:,:,:)
       Greal_ineq(iineq,:,:,:,:,:,:,:)  = impGreal(:,:,:,:,:,:,:)
       dens_ineq(iineq,1:Ilat,1:Norb)      = ed_dens(1:Ilat,1:Norb)
       docc_ineq(iineq,1:Ilat,1:Norb)      = ed_docc(1:Ilat,1:Norb)
       !mag_ineq(iineq,:,1:Norb)     = ed_mag(:,1:Norb)
       e_ineq(iineq,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       dd_ineq(iineq,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
       single_particle_density_matrix_ineq(iineq,:,:,:,:,:,:) = single_particle_density_matrix(:,:,:,:,:,:)
       cluster_density_matrix_ineq(iineq,:,:) = cluster_density_matrix(:,:)

    enddo
    if(MPI_MASTER)call stop_timer(unit=LOGfile)
    ed_file_suffix=""
    !
  end subroutine ed_solve_lattice_mpi
#endif
#endif














end module ED_MAIN








