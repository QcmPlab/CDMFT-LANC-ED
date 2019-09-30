program ed_hm_1dchain
  USE CDMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  !
  implicit none
  integer                                                                :: iloop,Nb,Lk,Nx
  logical                                                                :: converged
  real(8)                                                                :: wband,ts,tsp,wmixing,Eout(2),dens
  !Bath:
  real(8),allocatable                                                    :: Bath(:),BathOld(:)
  !The local hybridization function:
  complex(8),allocatable                                                 :: Hloc(:,:,:,:,:,:),Sig(:,:,:),SigA(:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Gmats,Greal,Smats,Sreal,Delta
  character(len=16)                                                      :: finput,fhloc
  real(8),allocatable                                                    :: wt(:),kxgrid(:),kygrid(:)
  complex(8),allocatable                                                 :: wm(:),wr(:)
  complex(8),allocatable                                                 :: Hk(:,:,:)
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer                                                                :: comm
  integer                                                                :: rank
  integer                                                                :: mpi_size
  logical                                                                :: master

  !Init MPI: use of MPI overloaded functions in SciFor
  call init_MPI(comm,.true.) !init MPI and print check message on screen
  rank   = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
  call parse_input_variable(tsp,"TSP",finput,default=0.0d0,comment="hopping parameter t prime")
  call parse_input_variable(Nx,"Nx",finput,default=10,comment="Number of kx point for 2d BZ integration")
  !
  call ed_read_input(trim(finput),comm)
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if (Nspin/=1.or.Norb/=1) stop "You are using too many spin-orbitals"
  Nlso=Nlat*Nspin*Norb
  if(.not.allocated(wm))allocate(wm(Lmats))
  wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)


  !Allocate Weiss Field:
  allocate(delta(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))


  !Build Hk
  call TB_set_bk(bkx=[pi2,0d0])
  Lk = Nx
  allocate(Hk(Nlso,Nlso,Lk),Wt(Lk),Hloc(1,1,1,1,1,1))
  call TB_build_model(Hk(:,:,:),hk_model,Nlso,[Nx,1])
  Wt = 1d0/Lk
  Hloc   = zero
  call TB_write_hk(Hk(:,:,:),"Hk2d_square_2nn.dat",1,&
       Nd=1,Np=0,Nineq=1,&
       Nkvec=[Nx,1])


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bathold(Nb))
  call ed_init_solver(comm,bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath) 
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)


     !Compute the local gfs:
     !call dmft_gloc_matsubara(comm,Hk,Wt,Gmats,Smats)
     !if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     !Get the Weiss field/Delta function to be fitted
     if(master)delta=ed_get_delta_matsubara(wm,bath)

     call Bcast_MPI(comm,Delta)


     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(master)then
        call ed_chi2_fitgf(delta,bath,ispin=1)
        call ed_chi2_fitgf(delta,bath,ispin=2)
        !
        !MIXING:
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*BathOld
        BathOld=Bath
        !
        !Check convergence (if required change chemical potential)
        converged = check_convergence(delta(:,:,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
        !if(nread/=0.d0)then
           !call ed_get_dens(dens,iorb=1)
           !call search_chemical_potential(xmu,dens,converged)
        !endif
     endif
     !
     call Bcast_MPI(comm,bath)
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo

  !Compute the local gfs:
  !call dmft_gloc_realaxis(comm,Hk,Wt,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)

  !Compute the Kinetic Energy:
  !call dmft_kinetic_energy(comm,Hk(:,:,:),Wt,Smats)


  call finalize_MPI()


contains

  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky
    complex(8)           :: hk(N,N)
    kx=kpoint(1)
    Hk = -2d0*ts*(cos(kx))
  end function hk_model


end program ed_hm_1dchain


