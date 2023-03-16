program cdn_bhz_2d
   USE CDMFT_ED
   !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   !USE MPI
   !
   implicit none
   integer                                                                :: Nx,Ny,Nz,Nlso,Nso,iloop,Nb,Nkx,Nky,Nkz,iw,iii,jjj,Nkpath
   integer                                                                :: il,jl,is,js,io,jo
   integer,dimension(2):: recover
   logical                                                                :: converged
   real(8)                                                                :: ts,Mh,lambda,wmixing,observable,repart,impart
   !Bath:
   real(8),allocatable                                                    :: Bath(:),Bath_fitted(:)
   !The local hybridization function:
   complex(8),allocatable                                                 :: Hloc(:,:)
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Gmats,Greal,Smats,Sreal,Weiss,Weiss_old
   character(len=16)                                                      :: finput
   real(8),allocatable                                                    :: wt(:),wr(:)
   complex(8),allocatable                                                 :: wm(:)
   complex(8),allocatable                                                 :: Hk(:,:,:)
   complex(8),dimension(:,:,:,:,:,:),allocatable                          :: observable_matrix
   !SYMMETRIES TEST
   real(8),dimension(:),allocatable                                       :: lambdasym_vector
   complex(8),dimension(:,:,:,:,:,:,:),allocatable                        :: Hsym_basis
   !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
   !integer                                                                :: comm
   !integer                                                                :: rank
   !integer                                                                :: mpi_size
   !logical                                                                :: master
   character(len=6)                                                       :: scheme
   type(finter_type)                                                      :: finter_func
   complex(8),dimension(:,:,:,:,:),allocatable                        :: dummy_real
   complex(8),dimension(2,2)                    :: test
   complex(8),dimension(:,:,:,:,:),allocatable                        :: dummy_mats

   !Init MPI: use of MPI overloaded functions in SciFor
   !call init_MPI(comm,.true.)
   !rank   = get_Rank_MPI(comm)
   !master = get_Master_MPI(comm)

   !

   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
   call parse_input_variable(Mh,"Mh",finput,default=1.d0,comment="crystal field splitting")
   call parse_input_variable(lambda,"lambda",finput,default=0.3d0,comment="spin-orbit coupling")
   call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
   call parse_input_variable(Nky,"Nky",finput,default=10,comment="Number of ku point for BZ integration")
  call parse_input_variable(Nkz,"Nkz",finput,default=10,comment="Number of kz point for BZ integration")
   call parse_input_variable(scheme,"SCHEME",finput,default="sigma")
   call parse_input_variable(nkpath,"NKPATH",finput,default=100)
   !
   call ed_read_input(trim(finput))
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

  print*,"Forcing Nx=1, Ny=1, Nz=2, Nspin=1, Norb=2"
  Nx=1
  Ny=1
  Nz=2
  Nspin=1
  Norb=2
  Nso=2
  Nlat=Nx*Ny*Nz
  Nlso=Nlat*Nspin*Norb
  if (Norb/=2) stop "Norb must be 2!"
  
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   !Allocate Weiss Field:
   allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(dummy_mats(Nspin,Nspin,Norb,Norb,Lmats),dummy_real(Nspin,Nspin,Norb,Norb,Lreal))

   !Build Hk and Hloc
   call generate_hk_hloc()
   allocate(lambdasym_vector(2))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,2))

   !SETUP SYMMETRIES (EXPERIMENTAL)
   lambdasym_vector(1)=Mh
   Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(hloc_model(Nlso,1.d0,0.d0,0.d0))
   !
   lambdasym_vector(2)=ts
   Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(hloc_model(Nlso,0.d0,1.d0,0.d0))
   !
   !
   !setup solver
   call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
   Nb=ed_get_bath_dimension(Hsym_basis)
   allocate(bath(Nb))
   allocate(bath_fitted(Nb))
   call ed_init_solver(bath)
   !
   call ed_read_impG
   call ed_read_impSigma
   !
   call ed_get_gimp_matsubara(Gmats)
   call ed_get_gimp_realaxis(Greal)
   call ed_get_sigma_matsubara(Smats)
   call ed_get_sigma_realaxis(Sreal)
   !
   !RETRIEVE AND PERIODIZE
   !call   print_hk_periodized_path()
   call   get_zeros()
   !


contains


   !+------------------------------------------------------------------+
   !PURPOSE  : Hloc for the 2d BHZ model
   !+------------------------------------------------------------------+


  function Hloc_model(N,Mh_,ts_,lambda_) result (H0)
    integer                                               :: N,ilat,jlat,ispin,iorb,jorb
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
    integer                                                      :: N,ilat,jlat,ispin,iorb,jorb,i,j
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
   !-------------------------------------------------------------------------------------------
   !PURPOSE: auxiliaries
   !-------------------------------------------------------------------------------------------

   function hk_periodized(kpoint,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kpoint
    complex(8),dimension(N,N) :: hk
    !
    Hk=(Mh-2*ts*(cos(kpoint(1))+cos(kpoint(2))+cos(kpoint(3))))*pauli_sigma_z+lambda*(sin(kpoint(1))*pauli_sigma_x+sin(kpoint(2))*pauli_sigma_y)
    !
   end function hk_periodized
   !

   !-------------------------------------------------------------------------------------------
   !PURPOSE: periodization M scheme
   !-------------------------------------------------------------------------------------------
   function periodize_sigma_Mscheme_real(kpoint,wprint) result(sreal_periodized)
      integer                                                     :: ilat,jlat,ispin,iorb,jorb,ii
      real(8),dimension(:)                                        :: kpoint
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: tmpmat
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: sreal_periodized,greal_periodized
      complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal)           :: g_lso,s_lso
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal) :: greal_unperiodized
      logical,optional                                            :: wprint
      logical                                                     :: wprint_
      !
      real(8)                                                     :: Mh_tmp,ts_tmp,lambda_tmp
      !
      wprint_=.true.;if(present(wprint))wprint_=wprint
      Mh_tmp=Mh
      ts_tmp=ts
      lambda_tmp=lambda
      !
      sreal_periodized=zero
      greal_periodized=zero
      !
      ts=0d0
      lambda=0d0
      !
      do ii=1,Lreal
         tmpmat=(dcmplx(wr(ii),eps)+xmu)*eye(Nlat*Nspin*Norb) - hk_model(kpoint,Nlat*Nspin*Norb) - nnn2lso(Sreal(:,:,:,:,:,:,ii))
         call inv(tmpmat)
         greal_unperiodized(:,:,:,:,:,:,ii)=lso2nnn(tmpmat)
      enddo
      !
      ts=ts_tmp
      lambda=lambda_tmp
      !
      do ii=1,Lreal
         do ilat=1,2  
            do jlat=1,2
               greal_periodized(:,:,:,:,ii)=greal_periodized(:,:,:,:,ii)+exp(-xi*kpoint(3)*(ilat-jlat))*greal_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
            enddo
         enddo
      enddo
      !
      do ii=1,Lreal
         g_lso(:,:,ii)=nn2so(greal_periodized(:,:,:,:,ii))
         call inv(g_lso(:,:,ii))
         Mh=0
         g_lso(:,:,ii)=g_lso(:,:,ii)-Hk_periodized(kpoint,Nspin*Norb)
         Mh=Mh_tmp
         s_lso(:,:,ii)=(dcmplx(wr(ii),eps)+xmu)*eye(Nspin*Norb)-Hk_periodized(kpoint,Nspin*Norb)-g_lso(:,:,ii)
         call inv(g_lso(:,:,ii))
         greal_periodized(:,:,:,:,ii)=so2nn(g_lso(:,:,ii))
         sreal_periodized(:,:,:,:,ii)=so2nn(s_lso(:,:,ii))
      enddo
      !
      if(wprint_)then
        do iorb=1,Nso
           do jorb=1,Nso
              call splot("perSigma_mscheme_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,s_lso(iorb,jorb,:))
           enddo
        enddo
        do iorb=1,Nso
           do jorb=1,Nso
              call splot("perG_mscheme_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,g_lso(iorb,jorb,:))
           enddo
        enddo
      endif
      !   
   end function periodize_sigma_Mscheme_real
   !
   function periodize_sigma_Mscheme_mats(kpoint,wprint) result(smats_periodized)
      integer                                                     :: ilat,jlat,ispin,iorb,jorb,ii
      real(8),dimension(:)                                        :: kpoint
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: tmpmat
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: smats_periodized,gmats_periodized
      complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats)           :: g_lso,s_lso
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: gmats_unperiodized
      logical,optional                                            :: wprint
      logical                                                     :: wprint_
      !
      real(8)                                                     :: Mh_tmp,ts_tmp,lambda_tmp
      !
      wprint_=.true.;if(present(wprint))wprint_=wprint
      Mh_tmp=Mh
      ts_tmp=ts
      lambda_tmp=lambda
      !
      smats_periodized=zero
      gmats_periodized=zero
      !
      ts=0d0
      lambda=0d0
      !
      do ii=1,Lmats
         tmpmat=(wm(ii)+xmu)*eye(Nlat*Nspin*Norb) - hk_model(kpoint,Nlat*Nspin*Norb) - nnn2lso(Smats(:,:,:,:,:,:,ii))
         call inv(tmpmat)
         gmats_unperiodized(:,:,:,:,:,:,ii)=lso2nnn(tmpmat)
      enddo
      !
      ts=ts_tmp
      lambda=lambda_tmp
      !
      do ii=1,Lmats
         do ilat=1,Nlat
            do jlat=1,Nlat
               gmats_periodized(:,:,:,:,ii)=gmats_periodized(:,:,:,:,ii)+exp(-xi*kpoint(3)*(ilat-jlat))*gmats_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
            enddo
         enddo
      enddo
      !
      do ii=1,Lmats
         g_lso(:,:,ii)=nn2so(gmats_periodized(:,:,:,:,ii))
         call inv(g_lso(:,:,ii))
         Mh=0
         g_lso(:,:,ii)=g_lso(:,:,ii)-Hk_periodized(kpoint,Nspin*Norb)
         Mh=Mh_tmp
         s_lso(:,:,ii)=(wm(ii)+xmu)*eye(Nspin*Norb)-Hk_periodized(kpoint,Nspin*Norb)-g_lso(:,:,ii)
         call inv(g_lso(:,:,ii))
         gmats_periodized(:,:,:,:,ii)=so2nn(g_lso(:,:,ii))
         smats_periodized(:,:,:,:,ii)=so2nn(s_lso(:,:,ii))
      enddo
      !
      if(wprint_)then
        do iorb=1,Nso
           do jorb=1,Nso
              call splot("perSigma_mscheme_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",imag(wm),s_lso(iorb,jorb,:))
           enddo
        enddo
        do iorb=1,Nso
           do jorb=1,Nso
              call splot("perG_mscheme_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",imag(wm),g_lso(iorb,jorb,:))
           enddo
        enddo
      endif
      !   
   end function periodize_sigma_Mscheme_mats
 
   !-------------------------------------------------------------------------------------------
   !PURPOSE: print routines
   !-------------------------------------------------------------------------------------------

   subroutine print_hk_periodized_path()
      integer                                :: i,j,Lk
      integer                                :: Npts
      real(8),dimension(:,:),allocatable     :: kpath
      character(len=64)                      :: file
      !
      write(LOGfile,*)"Build H(k) BHZ along path"
      !
      Npts = 5
      Lk=(Npts-1)*Nkpath
      allocate(kpath(Npts,3))
      kpath(1,:)=[0,0,0]*pi
      kpath(2,:)=[0,0,1]*pi
      kpath(3,:)=[0,1,1]*pi
      kpath(4,:)=[0,1,0]*pi
      kpath(5,:)=[0,0,0]*pi
      file="Eigenbands.nint"
      !
      if(allocated(Hk))deallocate(Hk)
      allocate(Hk(Nspin*Norb,Nspin*Norb,Lk))
      !
      call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
        call TB_Solve_model(hk_periodized,Nspin*Norb,kpath,Nkpath,&
         colors_name=[red1,blue1],&
         points_name=[character(len=20) :: 'G', 'Z', 'R', 'Y', 'G'],&
         file=reg(file))
      !
   end subroutine print_hk_periodized_path

  !
  !---------------------------------------------------------------------
  !PURPOSE: GET ZEROS ON THE REAL AXIS
  !---------------------------------------------------------------------
  subroutine get_zeros()
    integer                                       :: i,j,ik,ix,iy,Nso,Nktot,Npts
    integer                                       :: iorb,jorb
    integer                                       :: isporb,jsporb
    integer                                       :: ispin,jspin
    integer                                       :: iso,unit
    real(8),dimension(Nspin*Norb)                 :: dzeta
    complex(8),allocatable,dimension(:,:,:)       :: Hk_bare,Hk_topo
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: zeta,z_adj,fg,gdelta,fgk
    complex(8),dimension(:,:,:,:,:),allocatable   :: gloc
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Sigmareal,Sigmamats
    complex(8),dimension(:,:,:,:),allocatable     :: gk,gfoo,ReSmat
    complex(8)                                    :: iw
    complex(8),dimension(:,:),allocatable         :: detGiw
    real(8),dimension(Lreal)                      :: Den
    real(8),dimension(:),allocatable              :: Ipoles,Xcsign,Iweight
    real(8),dimension(:,:),allocatable            :: Ipoles3d,kpoints
    real(8),dimension(:,:),allocatable            :: Mpoles,Mweight
    real(8),dimension(:,:,:),allocatable          :: Mpoles3d
    integer                                       :: Linterval
    integer                                       :: count,Ninterval,maxNinterval,int
    real(8)                                       :: sign,sign_old
    real(8),dimension(:,:),allocatable :: kpath
    !
    Nso=Nspin*Norb
    !
    allocate(kpath(3,3))
    kpath(1,:)=[0.0,0.0,0.0]!G-e<-R
    kpath(2,:)=[0.0,0.0,1.0]!G
    kpath(3,:)=[0.0,0.0,2.0]!G
    !kpath(3,:)=[1.0,1.0]!G+e->R
    !kpath(4,:)=[0.0,1.0]!G+e->R
    !kpath(5,:)=[0.0,0.0]!G+e->R
    !kpath(4,:)=[1.0-0.35,0.0]!G-e<-R
    !kpath(5,:)=[1.0,0.0]!G
    !kpath(6,:)=[1.0+0.35,0.0]!G+e->R
    kpath=kpath*pi
    Npts  = size(kpath,1)
    Nktot = (Npts-1)*Nkpath
    !
    if(allocated(Hk_bare))deallocate(Hk_bare)
    allocate(Hk_bare(Nspin*Norb,Nspin*Norb,Nktot));Hk_bare=zero
    call TB_build_model(hk_bare,hk_periodized,Nspin*Norb,kpath,Nkpath)
    !
    !
    allocate(kpoints(Nktot,3))
    call TB_build_kgrid(kpath,Nkpath,kpoints)
    if(allocated(Sigmamats))deallocate(Sigmamats)
    if(allocated(Sigmareal))deallocate(Sigmareal)
    allocate(Sigmamats(Nktot,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sigmareal(Nktot,Nspin,Nspin,Norb,Norb,Lreal))
    do ik=1,Nktot
      Sigmamats(ik,:,:,:,:,:)=periodize_sigma_mscheme_mats(kpoints(ik,:),wprint=.false.)
      Sigmareal(ik,:,:,:,:,:)=periodize_sigma_mscheme_real(kpoints(ik,:),wprint=.false.)
    enddo
    !
    Linterval = 50000 !Maximum number of allowed intervals to look for zeros&poles
    !
    allocate(Xcsign(0:Linterval))
    allocate(Ipoles(Nktot),Iweight(Nktot))
    allocate(Mpoles(Nktot,Linterval),Mweight(Nktot,Linterval))
    !
    !FINDING THE POLES:
    !assume \eps=0.d0 ==> the ImSigma(poles)=0 this condition should be automatically
    !verified at the pole from definition of the pole (the ImSigma at the pole is just
    !an artificial broadening of an otherwise delta function, whose position should be 
    !determined by ReSigma only.
    Ipoles=0.d0   
    Mpoles=0.d0
    write(LOGfile,*)"Solving for the zeros..."
    maxNinterval=-1
    do ik=1,Nktot
       do i=1,Lreal
          zeta(:,:) = (wr(i)+xmu)*eye(Nso) - Hk_bare(:,:,ik) - nn2so(Sigmareal(ik,:,:,:,:,i))
          call inv(zeta)
          Den(i) = dreal(zeta(1,1))*dreal(zeta(2,2)) - dreal(zeta(1,2)*zeta(2,1))
       enddo
       Xcsign(0)=wr(1)
       count=0
       sign_old=sgn(Den(1))
       do i=2,Lreal
          sign=sgn(Den(i))
          if(sign*sign_old<1)then
             count=count+1
             if(count>Linterval)stop "Allocate Xcsign to a larger array."
             Xcsign(count)=wr(i)
          endif
          sign_old=sign
       enddo
       Ninterval=count
       if(count>maxNinterval)maxNinterval=count
       call init_finter(finter_func,wr,Den,3)
       do int=1,Ninterval
          Mpoles(ik,int) = brentq(det_poles,Xcsign(int-1),Xcsign(int))
          Mweight(ik,int)= get_weight(hk_bare(:,:,ik)-nn2so(Sigmamats(ik,:,:,:,:,1)))
       enddo
       !ipoles(ik) = brentq(det_poles,0.d0,wr(Lreal))
       !iweight(ik)= get_weight(hk_bare(:,:,ik)-nn2so(Sigmamats(ik,:,:,:,:,1)))
       call delete_finter(finter_func)
    enddo
    !call splot("BHZzeros.ed",ipoles,iweight)
    do int=1,maxNinterval
       unit=free_unit()
       open(unit,file="BHZzeros_int"//reg(txtfy(int))//".ed")
       if(any((Mpoles(:,int)/=0.d0)))then
          do ik=1,Nktot
             if(Mpoles(ik,int)/=0.d0)write(unit,*)kpoints(ik,1),kpoints(ik,2),kpoints(ik,3),Mpoles(ik,int),Mweight(ik,int)
          enddo
       endif
       close(unit)
    enddo
    !
  end subroutine get_zeros




  function det_poles(w) result(det)
    real(8),intent(in) :: w
    real(8)            :: det
    det = finter(finter_func,w)
  end function det_poles

  function get_weight(hk) result(wt)
    complex(8),dimension(4,4) :: hk,foo
    real(8),dimension(4)      :: eigv
    real(8) :: wt
    foo = hk
    call eigh(foo,eigv)
    wt = sum(foo(:,1))
  end function Get_Weight



 !+------------------------------------------------------------------+
 !PURPOSE  : Auxilliary reshape functions
 !+------------------------------------------------------------------+

   function so2nn(Hso) result(Hnn)
      complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
      complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
      integer                                     :: iorb,ispin,is
      integer                                     :: jorb,jspin,js
      Hnn=zero
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Norb
               do jorb=1,Norb
                  is = iorb + (ispin-1)*Norb  !spin-orbit stride
                  js = jorb + (jspin-1)*Norb  !spin-orbit stride
                  Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
               enddo
            enddo
         enddo
      enddo
   end function so2nn
   !
   function nn2so(Hnn) result(Hso)
      complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
      complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
      integer                                     :: iorb,ispin,is
      integer                                     :: jorb,jspin,js
      Hso=zero
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Norb
               do jorb=1,Norb
                  is = iorb + (ispin-1)*Norb  !spin-orbit stride
                  js = jorb + (jspin-1)*Norb  !spin-orbit stride
                  Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
               enddo
            enddo
         enddo
      enddo
   end function nn2so

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


end program cdn_bhz_2d








