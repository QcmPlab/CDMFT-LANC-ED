program hm_Nbands_bethe
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Nb,Le,Nso,iorb
  logical                                     :: converged
  real(8),dimension(5)                        :: Wbethe,Dbethe
  !Bath:
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !  
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: H0
  real(8),dimension(:),allocatable            :: de,dens
  !
  real(8),dimension(:),allocatable            :: Wband
  !
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal,Weiss_
  complex(8),allocatable,dimension(:) :: Gtest
  character(len=16)                           :: finput
  real(8)                                     :: wmixing,dw,dt,tmax
  !
  integer                                     :: comm,rank
  logical                                     :: master
  logical :: betheSC,wGimp,mixG0,symOrbs

  real(8),dimension(:),allocatable :: Aw,wr
  complex(8),dimension(:),allocatable :: Gw_gtr,Gt_gtr


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wbethe,"WBETHE",finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Dbethe,"DBETHE",finput,default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(betheSC,"BETHESC",finput,default=.false.)
  call parse_input_variable(wGimp,"wGimp",finput,default=.false.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.)
  !
  call ed_read_input(trim(finput),comm)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  if(Nspin/=1.OR.Norb>5)stop "Wrong setup from input file: Nspin=1 OR Norb>5"
  Nso=Nspin*Norb


  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(Wband(Nso))
  allocate(H0(Nso))
  allocate(de(Nso))
  !
  Wband = Wbethe(:Norb)
  H0    = Dbethe(:Norb)
  do iorb=1,Norb
     Ebands(iorb,:) = linspace(-Wband(iorb),Wband(iorb),Le,mesh=de(iorb))
     Dbands(iorb,:) = dens_bethe(Ebands(iorb,:),Wband(iorb))*de(iorb)
  enddo
  if(master)call TB_write_Hloc(one*diag(H0))


  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  Hloc(1,1,:,:)=diag(H0)


  !setup solver
  Nb=get_bath_dimension()
  allocate(bath(Nb))
  allocate(bath_(Nb))
  call ed_init_solver(comm,bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(comm,Ebands,Dbands,H0,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     !
     call dmft_gloc_realaxis(Comm,Ebands,Dbands,H0,Greal,Sreal)
     if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
     !
     !Get the Weiss field/Delta function to be fitted
     if(.not.betheSC)then
        call dmft_self_consistency(comm,Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)
     else
        if(wGimp)call ed_get_gimp_matsubara(Gmats)
        call dmft_self_consistency(comm,Gmats,Weiss,Hloc,SCtype=cg_scheme,wbands=Wband)
     endif
     if(master)call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
     !
     !
     !
     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(symOrbs)then
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1,iorb=1)
        call orb_equality_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
     endif
     !
     !MIXING:
     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     !
     !Check convergence (if required change chemical potential)
     if(master)then
        Gtest=zero
        do iorb=1,Norb
           Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
        enddo
        converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
        if(nread/=0d0)then
           call ed_get_dens(dens)
           call ed_search_variable(xmu,sum(dens),converged)
        endif
     endif
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo


  call dmft_gloc_realaxis(Comm,Ebands,Dbands,H0,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
  call dmft_kinetic_energy(comm,Ebands,Dbands,H0,Smats(1,1,:,:,:))



  if(.true.)then
     print*,"Get FFT of impG_realw"
     !1. Get the spectrum = ImGimp
     call ed_get_gimp_real(Greal)
     !
     allocate(Aw(Lreal),wr(Lreal))
     allocate(Gw_gtr(Lreal),Gt_gtr(0:Lreal))
     Aw = -dimag(Greal(1,1,1,1,:))/pi
     wr = linspace(wini,wfin,Lreal,mesh=dw)
     !2. Obtain  G>(w)
     Gw_gtr = pi2*xi*(fermi(wr,beta)-1d0)*Aw
     call splot("Gw_gtr.out",wr,Gw_gtr)
     !3. FFt G>(w) --> G>(t)
     Gt_gtr  =  fft_rw2rt(Gw_gtr)*dw/pi2
     dt      = pi/wfin
     tmax    = dt*Lreal/2
     call splot("Gt_gtr.out",linspace(-tmax,tmax,Lreal+1),Gt_gtr)
  endif


  call finalize_MPI()


contains



  !+----------------------------------------------------------------+
  !PURPOSE  : Modified real-axis FFT to deal with the special 0:L 
  ! form of the time-axis Keldysh GF.
  !+----------------------------------------------------------------+
  function fft_rw2rt(func_in) result(func_out)
    complex(8),dimension(:)               :: func_in
    complex(8),dimension(0:size(func_in)) :: func_out
    complex(8),dimension(size(func_in))   :: ftmp
    ftmp = func_in
    call fft(ftmp)
    call fftex(ftmp)
    func_out = fftshift(ftmp)*size(ftmp)
  end function fft_rw2rt
  !
  function fft_rt2rw(func_in) result(func_out)
    complex(8),dimension(:)               :: func_in
    complex(8),dimension(size(func_in)-1) :: func_out
    complex(8),dimension(size(func_in)-1)  :: ftmp
    ftmp = func_in(:size(func_in)-1)
    call ifft(ftmp)
    call fftex(ftmp)
    func_out = ifftshift(ftmp)
  end function fft_rt2rw

end program hm_Nbands_bethe



