program ed_hm_square_2nn
  USE DMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  !
  implicit none
  integer                                       :: iloop,Nb,Lk,Nx,Nso
  logical                                       :: converged
  real(8)                                       :: wband,ts,tsp,wmixing,Eout(2),dens
  !Bath:
  real(8),allocatable                           :: Bath(:),BathOld(:)
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:,:,:),Sig(:,:,:),SigA(:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats,Greal,Smats,Sreal,Delta
  character(len=16)                             :: finput,fhloc
  real(8),allocatable                           :: wt(:),kxgrid(:),kygrid(:),wm(:),wr(:)
  complex(8),allocatable                        :: Hk(:,:,:)
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer :: MpiComm
  integer :: MpiRank
  integer :: MpiSize
  logical :: MpiMaster

  !Init MPI: use of MPI overloaded functions in SciFor
  call init_MPI(MpiComm,.true.) !init MPI and print check message on screen
  MpiRank   = get_Rank_MPI(MpiComm)
  MpiMaster = get_Master_MPI(MpiComm)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
  call parse_input_variable(tsp,"TSP",finput,default=0.0d0,comment="hopping parameter t prime")
  call parse_input_variable(Nx,"Nx",finput,default=10,comment="Number of kx point for 2d BZ integration")
  !
  call ed_read_input(trim(finput),MpiComm)
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if (Nspin/=1.or.Norb/=1) stop "You are using too many spin-orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))


  !Build Hk
  call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
  Lk = Nx*Nx
  allocate(Hk(Nso,Nso,Lk),Wt(Lk),Hloc(1,1,1,1))
  call TB_build_model(Hk(:,:,:),hk_model,Nso,[Nx,Nx])
  Wt = 1d0/Lk
  Hloc   = zero
  call TB_write_hk(Hk(:,:,:),"Hk2d_square_2nn.dat",1,&
       Nd=1,Np=0,Nineq=1,&
       Nkvec=[Nx,Nx])


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bathold(Nb))
  call ed_init_solver(MpiComm,bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(MpiMaster)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(MpiComm,bath) 
     call ed_get_sigma_matsubara(Smats(:,:,:,:,:))
     call ed_get_sigma_real(Sreal(:,:,:,:,:))


     !Compute the local gfs:
     call dmft_gloc_matsubara(MpiComm,Hk,Wt,Gmats,Smats)
     if(MpiMaster)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     !Get the Weiss field/Delta function to be fitted
     if(MpiMaster)call dmft_self_consistency(Gmats,Smats,Delta,Hloc,cg_scheme)
     call Bcast_MPI(MpiComm,Delta)


     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(MpiMaster)then
        call ed_chi2_fitgf(delta,bath,ispin=1)
        !
        !MIXING:
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*BathOld
        BathOld=Bath
        !
        !Check convergence (if required change chemical potential)
        converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
        if(nread/=0.d0)then
           call ed_get_dens(dens,iorb=1)
           call search_chemical_potential(xmu,dens,converged)
        endif
     endif
     !
     call Bcast_MPI(MpiComm,bath)
     call Bcast_MPI(MpiComm,converged)
     call Bcast_MPI(MpiComm,xmu)
     !
     if(MpiMaster)call end_loop
  enddo

  !Compute the local gfs:
  call dmft_gloc_realaxis(MpiComm,Hk,Wt,Greal,Sreal)
  if(MpiMaster)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)

  !Compute the Kinetic Energy:
  call dmft_kinetic_energy(MpiComm,Hk(:,:,:),Wt,Smats)


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
    ky=kpoint(2)
    Hk = -one*2d0*ts*(cos(kx)+cos(ky))-one*4d0*tsp*(cos(kx)*cos(ky))
  end function hk_model


end program ed_hm_square_2nn



