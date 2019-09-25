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
  !
  !>INIT ED SOLVER
  !
  interface ed_init_solver
     module procedure :: ed_init_solver_single
#ifdef _MPI
     module procedure :: ed_init_solver_single_mpi
     \#endif
  end interface ed_init_solver
  !>
  public :: ed_init_solver


  !
  !> ED SOLVER
  !
  interface ed_solve
     module procedure :: ed_solve_single
#ifdef _MPI
     module procedure :: ed_solve_single_mpi
#endif
  end interface ed_solve
  !>
  public :: ed_solve


  real(8),dimension(:),allocatable                   :: wr,wm
  character(len=64)                                  :: suffix



contains





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_single(bath,Hloc)
    real(8),dimension(:),intent(inout) :: bath
    complex(8),intent(in)              :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    logical                            :: MPI_MASTER=.true.
    integer                            :: MPI_ERR
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    !Init bath:
    call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath,dreal(Hloc))
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath()
    call init_dmft_bath()
    call get_dmft_bath(bath)    !dmft_bath --> user_bath
    !
    if(isetup)call setup_global
    call deallocate_dmft_bath()
    isetup=.false.
    !
  end subroutine ed_init_solver_single
#ifdef _MPI
  subroutine ed_init_solver_single_mpi(MpiComm,bath,Hloc)
    integer                            :: MpiComm
    real(8),dimension(:),intent(inout) :: bath
    complex(8),intent(in)              :: Hloc(Nspin,Nspin,Norb,Norb)
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
    !Init bath:
    call set_hloc(Hloc)
    !
    check = check_bath_dimension(bath,dreal(Hloc))
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath()
    call init_dmft_bath()
    call get_dmft_bath(bath)
    !
    if(isetup)call setup_global
    call deallocate_dmft_bath()
    isetup=.false.
    !
    call ed_del_MpiComm()
    !
  end subroutine ed_init_solver_single_mpi
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
    complex(8),optional,intent(in)  :: Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    logical                         :: check
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    if(present(Hloc))call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath)    !user_bath --> dmft_bath
    call write_dmft_bath(LOGfile)
    if(MpiMaster)call save_dmft_bath(dmft_bath,used=.true.)
    !
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()         !find target states by digonalization of Hamiltonian
    call buildgf_impurity()             !build the one-particle impurity Green's functions  & Self-energy
    if(chiflag)call buildchi_impurity() !build the local susceptibilities (spin [todo charge])
    call observables_impurity()         !obtain impurity observables as thermal averages.          
    call local_energy_impurity()        !obtain the local energy of the effective impurity problem
    !
    call deallocate_dmft_bath(dmft_bath)
    select case(ed_diag_type)
    case default
       call es_delete_espace(state_list)
    case ("full")
       call delete_eigenspace()
    end select
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
    complex(8),optional,intent(in)  :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                         :: check
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm(MpiComm)
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    if(present(Hloc))call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    call write_dmft_bath(dmft_bath,LOGfile)
    if(MpiMaster)call save_dmft_bath(dmft_bath,used=.true.)
    !
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()         !find target states by digonalization of Hamiltonian
    call buildgf_impurity()             !build the one-particle impurity Green's functions  & Self-energy
    if(chiflag)call buildchi_impurity() !build the local susceptibilities (spin [todo charge])    
    call observables_impurity()         !obtain impurity observables as thermal averages.
    call local_energy_impurity()        !obtain the local energy of the effective impurity problem
    !
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    !
    nullify(spHtimesV_p)
  end subroutine ed_solve_single_mpi
#endif















end module ED_MAIN








