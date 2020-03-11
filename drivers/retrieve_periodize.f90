program cdn_bhz_2d
   USE CDMFT_ED
   !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nx,Ny,Nlso,iloop,Nb,Nkx,Nky,iw,iii,jjj,kkk
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
   real(8),allocatable                                                    :: wt(:)
   complex(8),allocatable                                                 :: wm(:),wr(:)
   complex(8),allocatable                                                 :: Hk(:,:,:),Smats_lso(:,:,:)
   complex(8),dimension(:,:,:,:,:,:),allocatable                          :: observable_matrix
   !SYMMETRIES TEST
   real(8),dimension(:),allocatable                                       :: lambdasym_vector
   complex(8),dimension(:,:,:,:,:,:,:),allocatable                        :: Hsym_basis
   !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
   integer                                                                :: comm
   integer                                                                :: rank
   integer                                                                :: mpi_size
   logical                                                                :: master,hermiticize
   character(len=6)                                                       :: scheme

   !Init MPI: use of MPI overloaded functions in SciFor
   call init_MPI(comm,.true.)
   rank   = get_Rank_MPI(comm)
   master = get_Master_MPI(comm)

   !

   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
   call parse_input_variable(Mh,"Mh",finput,default=1.d0,comment="crystal field splitting")
   call parse_input_variable(lambda,"lambda",finput,default=0.3d0,comment="spin-orbit coupling")
   call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
   call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of cluster sites in y direction")
   call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
   call parse_input_variable(Nky,"Nky",finput,default=10,comment="Number of ku point for BZ integration")
   call parse_input_variable(hermiticize,"HERMITICIZE",finput,default=.true.,comment="are bath replicas hermitian")
   call parse_input_variable(scheme,"SCHEME",finput,default="g",comment="Periodization scheme: possible g or sigma")
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
   !Ny=Nx
   !Nky=Nkx
   Nlat=Nx*Ny
   Nlso=Nlat*Nspin*Norb
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   do iii=1,Lreal
      wr(iii)=dcmplx(DREAL(wr(iii)),eps)
   enddo

   if(ED_VERBOSE > 0)call naming_convention()
   !Allocate Weiss Field:
   allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats_lso(Nlso,Nlso,Lmats))

   !Build Hk and Hloc
   call generate_hk_hloc()
   allocate(lambdasym_vector(3))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,3))
   
   !SETUP SYMMETRIES (EXPERIMENTAL)
   lambdasym_vector(1)=Mh
   Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(hloc_model(Nlso,1.d0,0.d0,0.d0))
   !
   lambdasym_vector(2)=ts
   Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(hloc_model(Nlso,0.d0,1.d0,0.d0))
   !
   lambdasym_vector(3)=lambda
   Hsym_basis(:,:,:,:,:,:,3)=lso2nnn(hloc_model(Nlso,0.d0,0.d0,1.d0))
   !
   !setup solver
   call set_Hloc(Hsym_basis,lambdasym_vector)
   Nb=get_bath_dimension(Hsym_basis)
   allocate(bath(Nb))
   allocate(bath_fitted(Nb))
   call ed_init_solver(comm,bath)
   !
   call ed_read_impG
   call ed_read_impSigma
   !
   call ed_get_gimp_matsubara(Gmats)
   call ed_get_gimp_realaxis(Greal)
   call ed_get_sigma_matsubara(Smats)
   call ed_get_sigma_realaxis(Sreal)
   !
   do iii=1,Lmats
      do il=1,Nlat
         do jl=1,Nlat
            do is=1,Nspin
               do js=1,Nspin
                  do io=1,Norb
                     do jo=1,Norb
                        repart=real(Smats(il,jl,is,js,io,jo,iii))
                        impart=imag(Smats(il,jl,is,js,io,jo,iii))
                        if(abs(repart).le.1d-6)repart=0.d0
                        if(abs(impart).le.1d-6)impart=0.d0
                        Smats(il,jl,is,js,io,jo,iii)=repart+xi*impart
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !
   !RETRIEVE AND PERIODIZE
   !call   print_hk_topological()
   call   print_hk_topological_path()
   !call dmft_gloc_matsubara(comm,Hk,Wt,Gmats,Smats)
   !if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
   !call dmft_gloc_realaxis(comm,Hk,Wt,Greal,Sreal)
   !if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)
   !call print_periodized([Nkx,Nky],hk_model,hk_periodized,scheme)

   call finalize_MPI()


contains

   !include "auxiliary_routines.f90"
   

   !+------------------------------------------------------------------+
   !PURPOSE  : Hloc for the 2d BHZ model
   !+------------------------------------------------------------------+


   function Hloc_model(N,Mh_,ts_,lambda_) result (H0)
      integer                                               :: N,ilat,jlat,ispin,iorb,jorb,ind1,ind2
      real(8)                                               :: Mh_,ts_,lambda_
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
      complex(8),dimension(N,N)                             :: H0
      character(len=64)                                     :: file_
      integer                                               :: unit
      file_ = "tcluster_matrix.dat"
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do ilat=1,Nx
            do jlat=1,Ny
               ind1=indices2N([ilat,jlat])
               hopping_matrix(ind1,ind1,ispin,ispin,:,:)= t_m(mh_)
               if(ilat<Nx)then
                  ind2=indices2N([ilat+1,jlat])
                  hopping_matrix(ind1,ind2,ispin,ispin,:,:)= t_x(ts_,lambda_,ispin)
               endif
               if(ilat>1)then
                  ind2=indices2N([ilat-1,jlat])
                  hopping_matrix(ind1,ind2,ispin,ispin,:,:)= dconjg(transpose(t_x(ts_,lambda_,ispin)))
               endif
               if(jlat<Ny)then
                  ind2=indices2N([ilat,jlat+1])
                  hopping_matrix(ind1,ind2,ispin,ispin,:,:)= t_y(ts_,lambda_)
               endif
               if(jlat>1)then
                  ind2=indices2N([ilat,jlat-1])
                  hopping_matrix(ind1,ind2,ispin,ispin,:,:)= transpose(t_y(ts_,lambda_))
               endif
            enddo
         enddo
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
         do ilat=1,Ny
            ind1=indices2N([1,ilat])
            ind2=indices2N([Nx,ilat])
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + dconjg(transpose(t_x(ts,lambda,ispin)))*exp(xi*kpoint(1)*Nx)
            hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + t_x(ts,lambda,ispin)*exp(-xi*kpoint(1)*Nx)
         enddo
         do ilat =1,Nx
            ind1=indices2N([ilat,1])
            ind2=indices2N([ilat,Ny])
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + transpose(t_y(ts,lambda))*exp(xi*kpoint(2)*Ny)
            hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + t_y(ts,lambda)*exp(-xi*kpoint(2)*Ny)
         enddo
      enddo
      !
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

   function t_x(hop1,hop2,spinsign) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: hop1,hop2,sz
      integer                         :: spinsign
      !
      tmpmat=zero
      sz=(-1.d0)**(spinsign+1)
      tmpmat=-hop1*pauli_sigma_z+0.5d0*sz*xi*hop2*pauli_sigma_x
      !
   end function t_x

   function t_y(hop1,hop2) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: hop1,hop2
      !
      tmpmat=zero
      tmpmat=-hop1*pauli_sigma_z
      tmpmat(1,2)=-hop2*0.5d0
      tmpmat(2,1)=hop2*0.5d0
      !
   end function t_y
   
   !-------------------------------------------------------------------------------------------
   !PURPOSE: periodization & print_hk_topological
   !-------------------------------------------------------------------------------------------

   function hk_periodized(kpoint,N) result(Hk)
      real(8),dimension(:)                          :: kpoint
      integer                                       :: Nlat_,Nx_,Ny_,N
      complex(8),dimension(N,N)                     :: Hk
      complex(8),dimension(Nspin,Nspin,Norb,Norb)   :: tmpmat
      !
      if(scheme=="g")then
         tmpmat=periodize_g(kpoint)
         Hk=nn2so(tmpmat)
      else
         tmpmat=periodize_sigma(kpoint)
         Nlat_=Nlat
         Nx_=Nx
         Ny_=Ny
         Nlat=1
         Nx=1
         Ny=1
         !
         Hk=hk_model(kpoint,Nspin*Norb)+nn2so(tmpmat)
         !
         Nlat=Nlat_
         Nx=Nx_
         Ny=Ny_
      endif
      !
   end function hk_periodized
   
   !

   function periodize_sigma(kpoint) result(smats_periodized_omegazero)
      integer                                                     :: ilat,jlat,ispin,iorb,ii
      real(8),dimension(:)                                        :: kpoint
      integer,dimension(:),allocatable                            :: ind1,ind2
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: smats_periodized
      complex(8),dimension(Nspin,Nspin,Norb,Norb)                 :: smats_periodized_omegazero
      !
      !
      if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
      if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
      !
      smats_periodized=zero
      !
      do ii=1,1!Lmats
         do ilat=1,Nlat
            ind1=N2indices(ilat)        
            do jlat=1,Nlat
               ind2=N2indices(jlat)
               smats_periodized(:,:,:,:,ii)=smats_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*Smats(ilat,jlat,:,:,:,:,ii)/Nlat
            enddo
         enddo
      enddo
      !
      smats_periodized_omegazero=smats_periodized(:,:,:,:,1)
      !
      !   
   end function periodize_sigma


   function periodize_g(kpoint) result(smats_periodized_omegazero)
      integer                                                     :: ilat,jlat,ispin,iorb,ii,itest,Ntest
      integer                                                     :: Nlat_,Nx_,Ny_,N
      real(8),dimension(:)                                        :: kpoint
      integer,dimension(:),allocatable                            :: ind1,ind2
      complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats)           :: invG0mats,invGmats
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: tmpmat,Hk_unper
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: gmats_unperiodized![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: smats_periodized,gmats_periodized
      complex(8),dimension(Nspin,Nspin,Norb,Norb)                 :: smats_periodized_omegazero
      !
      !
      if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
      if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
      !
      gmats_unperiodized=zero
      gmats_periodized=zero
      tmpmat=zero
      !
      do ii=1,1!Lmats
         tmpmat=(xi*wm(ii)+xmu)*eye(Nlat*Nspin*Norb) - hk_model(kpoint,Nlat*Nspin*Norb) - nnn2lso(Smats(:,:,:,:,:,:,ii))
         call inv(tmpmat)
         gmats_unperiodized(:,:,:,:,:,:,ii)=lso2nnn(tmpmat)
      enddo
      !
      Ntest=1
      !
      do ii=1,1!Lmats
         do ilat=1,Nlat
            ind1=N2indices(ilat)        
            do jlat=1,Nlat
               ind2=N2indices(jlat)
               do itest=1,Ntest
                  gmats_periodized(:,:,:,:,ii)=gmats_periodized(:,:,:,:,ii)&
                     +exp(-xi*dot_product(kpoint,(ind1-ind2))*itest/Ntest)*(itest/Ntest)*gmats_unperiodized(ilat,jlat,:,:,:,:,ii)/(Nlat*Ntest)
               enddo
            enddo
         enddo
      enddo
      !
      !
      !Get G0^-1
      !
      !Nlat_=Nlat
      !Nx_=Nx
      !Ny_=Ny
      !Nlat=1
      !Nx=1
      !Ny=1
      !
      !do ii=1,1!Lmats
      !   invG0mats(:,:,ii) = (xi*wm(ii)+xmu)*eye(Nspin*Norb)  - Hk_model(kpoint,Nspin*Norb)            
      !enddo
      !
      !Nlat=Nlat_
      !Nx=Nx_
      !Ny=Ny_
      !
      !Get G^-1
      !
      do ii=1,1!Lmats
         invGmats(:,:,ii) = nn2so(gmats_periodized(:,:,:,:,ii))
         call inv(invGmats(:,:,ii))
      enddo
      !
      !Get Sigma functions: Sigma= G0^-1 - G^-1
      !Smats_periodized=zero
      !
      !do ii=1,1!Lmats
      !   Smats_periodized(:,:,:,:,ii) = so2nn(invG0mats(:,:,ii) - invGmats(:,:,ii))
      !enddo
      !
      !smats_periodized_omegazero=smats_periodized(:,:,:,:,1)
      smats_periodized_omegazero=so2nn(-invGmats(:,:,1))
      !   
   end function periodize_g


   subroutine print_hk_topological()
      integer                                     :: ikx,iky,unit
      real(8)                                     :: kx,ky
      real(8),dimension(Nspin*Norb)               :: eigvals
      complex(8),dimension(Nspin*Norb,Nspin*Norb) :: H_kpoint
      !
      !
      unit=free_unit()
      open(unit,file="bands2d.ed")
      do ikx=1,Nkx
         kx=-pi+2*pi/Nkx*ikx
         do iky=1,Nky
            ky=-pi+2*pi/Nky*iky
            H_kpoint=Hk_periodized([kx,ky],Nspin*Norb)
            call Eigh(H_kpoint,eigvals)
            write(unit,*)kx,ky,eigvals
         enddo
            write(unit,*)" "
      enddo
      close(unit)
      !
   end subroutine print_hk_topological

   subroutine print_hk_topological_path()
    integer                                :: i,j,Lk,Nkpath
    integer                                :: Npts
    real(8),dimension(:,:),allocatable     :: kpath
    character(len=64)                      :: file
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is constructed along this path.
       if(master)write(LOGfile,*)"Build H(k) haldane_square along the path GXMG:"
       Npts = 7
       Nkpath=500
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,2))
       kpath(1,:)=-kpoint_X2(1:2)
       kpath(2,:)=kpoint_Gamma(1:2)
       kpath(3,:)=kpoint_X2(1:2)
       kpath(4,:)=kpoint_M1(1:2)
       kpath(5,:)=kpoint_X1(1:2)
       kpath(6,:)=kpoint_Gamma(1:2)
       kpath(7,:)=-kpoint_X1(1:2)
       file="Eigenbands.nint"
       !
       if(allocated(Hk))deallocate(Hk)
       allocate(Hk(Nspin*Norb,Nspin*Norb,Lk))
       !
       call TB_set_bk([pi2,0d0],[0d0,pi2])
       if(master)  call TB_Solve_model(hk_periodized,Nspin*Norb,kpath,Nkpath,&
         colors_name=[red1,blue1],&
         points_name=[character(len=20) :: '-Y', 'G', 'Y', 'M', 'X', 'G', '-X'],&
         file=reg(file))
      !
   end subroutine print_hk_topological_path

   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate Hloc and Hk
   !-------------------------------------------------------------------------------------------

   subroutine generate_hk_hloc()
      integer                                     :: ik
      real(8),dimension(Nkx*Nky,2)                :: kgrid
      !
      call TB_build_kgrid([Nkx,Nky],kgrid)
      kgrid(:,1)=kgrid(:,1)/Nx
      kgrid(:,2)=kgrid(:,2)/Ny
      !
      if(allocated(hk))deallocate(hk)
      if(allocated(Wt))deallocate(Wt)
      if(allocated(hloc))deallocate(Hloc)
      !
      allocate(Hk(Nlso,Nlso,Nkx*Nky),Wt(Nkx*Nky),Hloc(Nlso,Nlso))
      hk=zero
      wt=zero
      hloc=zero
      !
      call TB_build_model(Hk,hk_model,Nlso,kgrid)
      Wt = 1d0/(Nkx*Nky)
      Hloc=hloc_model(Nlso,Mh,ts,lambda)
      where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
         !
   end subroutine generate_hk_hloc

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


   function indices2N(indices) result(N)
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      !
      N=Nx*(indices(2)-1)+indices(1)
   end function indices2N

   function N2indices(N) result(indices) 
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      indices(1)=mod(N,Nx)
      if(indices(1)==0)then
         indices(1)=Nx
         indices(2)=(N-Nx)/Nx+1
      else
         indices(2)=N/Nx+1
      endif
   end function N2indices

   subroutine naming_convention()
      integer                       :: i,j
      integer,dimension(Nx,Ny)      :: matrix
      !
      do j=1,Ny
         do i=1,Nx
            matrix(i,j)=indices2N([i,j])
         enddo
      enddo
      !
      write(LOGfile,"(A)")"The unique index of each site (on the cartesian plane) is as follows:"
      write(LOGfile,"(A)")" "
      do j=1,Ny
         write(LOGfile,"(20(I2,2x))")(matrix(i,Ny+1-j),i =1,Nx)
      enddo
      write(LOGfile,"(A)")" "
   end subroutine naming_convention
   !
end program cdn_bhz_2d








