program cdn_weyl
  USE CDMFT_ED !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  !
  implicit none
  integer                                         :: Nx,Ny,Nz,Nlso,iloop,Nb,Nkx,Nky,Nkz,iw,iii,jjj,kkk,Ntr
  integer,dimension(2)                            :: recover
  logical                                         :: converged
  real(8)                                         :: ts,Mh,lambda,wmixing,observable
  !Bath:
  real(8),allocatable                             :: Bath(:),Bath_prev(:)
  !The local hybridization function:
  complex(8),allocatable                          :: Hloc(:,:)
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal,Smats,Sreal,Weiss,Weiss_old
  character(len=16)                               :: finput
 ! real(8),allocatable                             :: wt(:)
  complex(8),allocatable                          :: wm(:),wr(:)
  complex(8),allocatable                          :: Hk(:,:,:),Smats_lso(:,:,:)
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: observable_matrix
  !Density matrices:
  complex(8),allocatable,dimension(:,:)           :: reduced_density_matrix
  complex(8),allocatable,dimension(:,:)           :: pure_cvec
  real(8),allocatable,dimension(:)                :: pure_prob
  !SYMMETRIES TEST
  real(8),dimension(:),allocatable                :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Hsym_basis
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer                                         :: comm
  integer                                         :: rank
  integer                                         :: mpi_size
  logical                                         :: master,hermiticize,dmexit

  !Init MPI: use of MPI overloaded functions in SciFor
  call init_MPI(comm,.true.)
  rank   = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  !

  !Parse input variables
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.5d0,comment="hopping parameter")
  call parse_input_variable(Mh,"Mh",finput,default=1.d0,comment="crystal field splitting")
  call parse_input_variable(lambda,"lambda",finput,default=0.3d0,comment="spin-orbit coupling")
  call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
  call parse_input_variable(Nky,"Nky",finput,default=10,comment="Number of ky point for BZ integration")
  call parse_input_variable(Nkz,"Nkz",finput,default=10,comment="Number of kz point for BZ integration")
  call parse_input_variable(hermiticize,"HERMITICIZE",finput,default=.true.,comment="are bath replicas hermitian")
  !call parse_input_variable(dmexit,"DMEXIT",finput,default=.true.,comment="Exit after evaluating DM, behave as post-processing")
  !
  call ed_read_input(trim(finput),comm)
  !
  !Add dmft control variables
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  !set global variables
  !if (Nspin/=1.or.Norb/=1) stop "You are using too many spin-orbitals"
  print*,"Forcing Nx=1, Ny=1, Nz=2, Nspin=1, Norb=2"
  Nx=1
  Ny=1
  Nz=2
  Nspin=1
  Norb=2
  Nlat=Nx*Ny*Nz
  Nlso=Nlat*Nspin*Norb
  if (Norb/=2) stop "Norb must be 2!"

  if(.not.allocated(wm))allocate(wm(Lmats))
  wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)

  !Allocate Weiss Field:
  allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats_lso(Nlso,Nlso,Lmats))

  !Build Hk and Hloc
  call generate_hk_hloc()

  !SETUP BATH STEP 1
  allocate(lambdasym_vector(2))
  allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,2))
  !
  lambdasym_vector(1)=Mh
  Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(hloc_model(Nlso,1.d0,0.d0,0.d0))
  !
  lambdasym_vector(2)=ts
  Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(hloc_model(Nlso,0.d0,1.d0,0.d0))
  !
  !
  !SETUP BATH STEP 2 and SETUP SOLVER
  call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
  Nb=ed_get_bath_dimension(Hsym_basis)
  allocate(bath(Nb))
  allocate(bath_prev(Nb))
  call ed_init_solver(comm,bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)

     call ed_solve(comm,bath,lso2nnn(Hloc))
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)

     !Compute the local gfs:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     !
     !Get the Weiss field/Delta function to be fitted
     call dmft_self_consistency(Gmats,Smats,Weiss,lso2nnn(Hloc),cg_scheme)
     call Bcast_MPI(comm,Weiss)
     !
     !MIXING:
     !if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_Old
     !Weiss_old=Weiss
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(master)then
        call ed_chi2_fitgf(Weiss,bath)
        !
        !MIXING:
        !call adaptive_mix(Bath(Nbath+1:),Bath_fitted(Nbath+1:)-Bath(Nbath+1:),wmixing,iloop)
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
        Bath_Prev=Bath
        !
        !Check convergence (if required change chemical potential)
        converged = check_convergence(Weiss(:,:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     endif
     !
     call Bcast_MPI(comm,bath)
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo

  !Compute the local gfs:
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)

  !Compute the Kinetic Energy:
  do iw=1,Lmats
     Smats_lso(:,:,iw)=nnn2lso(Smats(:,:,:,:,:,:,iw))
  enddo
  call dmft_kinetic_energy(Hk(:,:,:),Smats_lso)


  call finalize_MPI()


contains

  !+------------------------------------------------------------------+
  !PURPOSE  : Hloc for the 2d BHZ model
  !+------------------------------------------------------------------+


  function Hloc_model(N,Mh_,ts_,lambda_) result (H0)
    integer                                               :: N,ilat,jlat,ispin,iorb,jorb,ind1,ind2
    real(8)                                               :: Mh_,ts_,lambda_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
    complex(8),dimension(N,N)                             :: H0
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
       hopping_matrix(1,1,ispin,ispin,:,:)= t_m(mh_)
       hopping_matrix(2,2,ispin,ispin,:,:)= t_m(mh_)
       hopping_matrix(1,2,ispin,ispin,:,:)= dconjg(transpose(t_z(ts_)))
       hopping_matrix(2,1,ispin,ispin,:,:)= t_z(ts_)
    enddo
    !
    H0=nnn2lso(hopping_matrix)
    !
  end function hloc_model


  function hk_model(kpoint,N) result(Hk)
    integer                                                      :: N,ilat,jlat,ispin,iorb,jorb,i,j,ind1,ind2
    real(8),dimension(:)                                         :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
    complex(8),dimension(N,N)                                    :: hk
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
      do ilat=1,2
          hopping_matrix(ilat,ilat,ispin,ispin,:,:)=hopping_matrix(ilat,ilat,ispin,ispin,:,:) +&
                                                    2*ts*pauli_sigma_z*(-cos(kpoint(1))-cos(kpoint(2))) + &
                                                    pauli_sigma_x*lambda*sin(kpoint(1)) + &
                                                    pauli_sigma_y*lambda*sin(kpoint(2))
      enddo
      hopping_matrix(1,2,ispin,ispin,:,:)= dconjg(transpose(t_z(ts)))*exp(xi*kpoint(3)*Nz)
      hopping_matrix(2,1,ispin,ispin,:,:)= t_z(ts)*exp(-xi*kpoint(3)*Nz)
    enddo
    
    
    Hk=nnn2lso(hopping_matrix)+hloc_model(N,Mh,ts,lambda)
    !
  end function hk_model


  !AUXILLIARY HOPPING MATRIX CONSTRUCTORS

  function t_m(mass) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: mass
    !
    tmpmat=zero
    tmpmat=mass*pauli_sigma_z
    !
  end function t_m


  function t_x(hop1) result(tmpmat)
   complex(8),dimension(Norb,Norb) :: tmpmat
   real(8)                         :: hop1

   !
   tmpmat=zero
   tmpmat=-hop1*pauli_sigma_x
   !
 end function t_x

 function t_y(hop1) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: hop1
    !
    tmpmat=zero
    tmpmat=-hop1*pauli_sigma_y
    !
  end function t_y
  
  
  function t_z(hop1) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: hop1

    !
    tmpmat=-hop1*pauli_sigma_z
    !
  end function t_z

  !-------------------------------------------------------------------------------------------
  !PURPOSE: generate Hloc and Hk
  !-------------------------------------------------------------------------------------------

  subroutine generate_hk_hloc()
    integer                                     :: ik
    real(8),dimension(Nkx*Nky*Nkz,3)            :: kgrid
    real(8),dimension(3)                        :: e1,e2,e3,bk1,bk2,bk3
    real(8)                                     :: bklen
    !
    e1 = [1d0, 0d0, 0d0]
    e2 = [0d0, 1d0, 0d0]
    e3 = [0d0, 0d0, 2d0]
    call TB_set_ei(eix=e1,eiy=e2,eiz=e3)
    bklen=2d0*pi
    bk1=bklen*[1d0, 0d0, 0d0]
    bk2=bklen*[0d0, 1d0, 0d0]
    bk3=bklen*[0d0, 0d0, 0.5d0]
    call TB_set_bk(bkx=bk1,bky=bk2,bkz=bk3)
    !
    call TB_build_kgrid([Nkx,Nky,Nkz],kgrid)
    !
    if(allocated(hk))deallocate(hk)
    if(allocated(hloc))deallocate(Hloc)
    !
    allocate(Hk(Nlso,Nlso,Nkx*Nky*Nkz),Hloc(Nlso,Nlso))
    hk=zero
    hloc=zero
    !
    call TB_build_model(Hk,hk_model,Nlso,kgrid)
    Hloc=hloc_model(Nlso,Mh,ts,lambda)
    !
  end subroutine generate_hk_hloc




  !+------------------------------------------------------------------+
  !PURPOSE  : Auxilliary reshape functions
  !+------------------------------------------------------------------+

  function lso2nnn(Hlso) result(Hnnn)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                               :: ilat,jlat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hnnn=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn


  function nnn2lso(Hnnn) result(Hlso)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: ilat,jlat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function nnn2lso


end program cdn_weyl








