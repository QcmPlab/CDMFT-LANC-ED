program cdn_bhz_2d
   USE CDMFT_ED
   !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nx,Ny,Nlso,iloop,Nb,Nkx,Nky,iw,iii,jjj,Nkpath
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
   integer                                                                :: comm
   integer                                                                :: rank
   integer                                                                :: mpi_size
   logical                                                                :: master
   type(finter_type)                                                      :: finter_func

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
   call parse_input_variable(nkpath,"NKPATH",finput,default=100)
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
   !
   if(ED_VERBOSE > 0)call naming_convention()
   !Allocate Weiss Field:
   allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))

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
   !
   !RETRIEVE AND PERIODIZE
   !
   call   print_hk_periodized_path()
   call   print_hk_topological_path()
   call   get_Akw()
   call   get_poles()
   !
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

   function hk_periodized(kvec,N) result(hk)
    integer                   :: N,ii
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky
    kx=kvec(1)
    ky=kvec(2)
    Hk  = zero
    Hk(1:Norb,1:Norb) = hk_bhz2x2(kx,ky)
    if(Nspin .eq. 2)then
      Hk(Nspin+1:Nspin*Norb,Nspin+1:Nspin*Norb) = conjg(hk_bhz2x2(-kx,-ky))
    endif
   end function hk_periodized

   function hk_bhz2x2(kx,ky) result(hk)
     real(8)                   :: kx,ky,epsik
     complex(8),dimension(2,2) :: hk
     epsik   = 2.d0*ts*(cos(kx)+cos(ky))
     hk(1,1) = mh - epsik
     hk(2,2) =-mh + epsik
     hk(1,2) = lambda*(sin(kx)-xi*sin(ky))
     hk(2,1) = lambda*(sin(kx)+xi*sin(ky))
   end function hk_bhz2x2
   !
   !
   function hk_topological(kpoint,N) result(Hk)
      real(8),dimension(:)                                         :: kpoint
      integer                                                      :: Nlat_,Nx_,Ny_,N,i
      complex(8),dimension(N,N)                                    :: Hk
      complex(8),dimension(Nspin*Norb,Nspin*Norb)                  :: Zmats,Sigma_lso
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)            :: sigmamat
      !
      sigmamat=periodize_sigma_mats(kpoint)
      sigma_lso=nn2so(sigmamat(:,:,:,:,1))
      Hk=hk_periodized(kpoint,Nspin*Norb)+sigma_lso
      !
      Zmats=abs( zeye(Nspin*Norb) - IMAG(sigma_lso)/(pi/beta))
      call inv(Zmats)
      Hk = matmul(Zmats,Hk)
      !print*,"Z11",Zmats(1,1)
      !print*,"Z22",Zmats(2,2)
      !print*,"Z12",Zmats(1,2)
      !print*,"Z21",Zmats(2,1)
      !print*,""
      !
   end function hk_topological
   !
   !
   !
   function periodize_sigma_mats(kpoint) result(smats_periodized)
      integer                                                     :: ilat,jlat,ispin,iorb,ii
      real(8),dimension(:)                                        :: kpoint
      integer,dimension(:),allocatable                            :: ind1,ind2
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: smats_periodized
      !
      !
      if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
      if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
      !
      smats_periodized=zero
      !
      do ii=1,Lmats
         do ilat=1,Nlat
            ind1=N2indices(ilat)        
            do jlat=1,Nlat
               ind2=N2indices(jlat)
               smats_periodized(:,:,:,:,ii)=smats_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*Smats(ilat,jlat,:,:,:,:,ii)/Nlat
            enddo
         enddo
      enddo
      !   
   end function periodize_sigma_mats
   !
   !
   !
   function periodize_sigma_real(kpoint) result(sreal_periodized)
      integer                                                     :: ilat,jlat,ispin,iorb,ii
      real(8),dimension(:)                                        :: kpoint
      integer,dimension(:),allocatable                            :: ind1,ind2
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: sreal_periodized
      !
      !
      if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
      if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
      !
      sreal_periodized=zero
      !
      do ii=1,Lreal
         do ilat=1,Nlat
            ind1=N2indices(ilat)        
            do jlat=1,Nlat
               ind2=N2indices(jlat)
               sreal_periodized(:,:,:,:,ii)=sreal_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*Sreal(ilat,jlat,:,:,:,:,ii)/Nlat
            enddo
         enddo
      enddo
      !
      !   
   end function periodize_sigma_real
   !
   !
   subroutine print_hk_periodized_path()
      integer                                :: i,j,Lk
      integer                                :: Npts
      real(8),dimension(:,:),allocatable     :: kpath
      character(len=64)                      :: file
      !
      if(master)write(LOGfile,*)"Build H(k) BHZ along path"
      !
      Npts = 7
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
   end subroutine print_hk_periodized_path
   !
   subroutine print_hk_topological_path()
      integer                                :: i,j,Lk
      integer                                :: Npts
      real(8),dimension(:,:),allocatable     :: kpath
      character(len=64)                      :: file
      !
      if(master)write(LOGfile,*)"Build H(k) BHZ along path"
      !
      Npts = 7
      Lk=(Npts-1)*Nkpath
      allocate(kpath(Npts,2))
      kpath(1,:)=-kpoint_X2(1:2)
      kpath(2,:)=kpoint_Gamma(1:2)
      kpath(3,:)=kpoint_X2(1:2)
      kpath(4,:)=kpoint_M1(1:2)
      kpath(5,:)=kpoint_X1(1:2)
      kpath(6,:)=kpoint_Gamma(1:2)
      kpath(7,:)=-kpoint_X1(1:2)
      file="Eig_Htop.ed"
      !
      if(allocated(Hk))deallocate(Hk)
      allocate(Hk(Nspin*Norb,Nspin*Norb,Lk))
      !
      call TB_set_bk([pi2,0d0],[0d0,pi2])
      if(master)  call TB_Solve_model(hk_topological,Nspin*Norb,kpath,Nkpath,&
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
   !PURPOSE  : Get A(k,w) along the Y-G-X direction
   !+------------------------------------------------------------------+
   !
   !
   subroutine get_Akw()
      integer                                       :: ik,ispin,iorb,ilat
      integer                                       :: Npts,Nktot,unit
      !
      complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sigma
      complex(8),allocatable,dimension(:,:,:)       :: Hk_bare,Hk_topo
      real(8),allocatable,dimension(:,:)            :: Eigval_bare,Eigval_topo
      complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gkreal
      real(8),allocatable,dimension(:,:)            :: Akreal
      real(8),dimension(:),allocatable              :: knumber
      real(8),dimension(:,:),allocatable            :: kpath,kpoints
      !
      if(master)write(LOGfile,*)"Build A(k,w) BHZ along path"
      !
      !path: Y G X
      allocate(kpath(3,2))
      kpath(1,:)=[0d0, pi]
      kpath(2,:)=[0d0,0d0]
      kpath(3,:)=[ pi,0d0]
      !
      !get kpoints
      Npts  = size(kpath,1)
      Nktot = (Npts-1)*Nkpath
      allocate(knumber(Nktot))
      !
      !Explicitly write kpoints (needed for Sigma)
      allocate(kpoints(Nktot,2))
      call TB_build_kgrid(kpath,Nkpath,kpoints)
      !
      !Generate Sigma(k,w) along path
      allocate(Sigma(Nktot,Nspin,Nspin,Norb,Norb,Lreal))
      do ik=1,Nktot
         Sigma(ik,:,:,:,:,:)=periodize_sigma_real(kpoints(ik,:))
      enddo
      !
      !allocate Hamiltonian and build model along path
      allocate(Hk_bare(Nspin*Norb,Nspin*Norb,Nktot));Hk_bare=zero
      allocate(Hk_topo(Nspin*Norb,Nspin*Norb,Nktot));Hk_topo=zero
      allocate(eigval_bare(Nspin*Norb,Nktot));eigval_bare=zero
      allocate(eigval_topo(Nspin*Norb,Nktot));eigval_topo=zero
      call TB_build_model(hk_bare,hk_periodized,Nspin*Norb,kpath,Nkpath)
      call TB_build_model(hk_topo,hk_topological,Nspin*Norb,kpath,Nkpath)
      !
      !allocate and compute Gkw
      allocate(Gkreal(Nktot,Nspin,Nspin,Norb,Norb,Lreal))
      do ik=1,Nktot
         knumber(ik)=ik
         call eigh(hk_topo(:,:,ik),eigval_topo(:,ik))
         call dmft_gk_realaxis(Hk_bare(:,:,ik),1d0,Gkreal(ik,:,:,:,:,:),Sigma(ik,:,:,:,:,:)) 
         call eigh(hk_bare(:,:,ik),eigval_bare(:,ik)) !THIS AFTER GK!!!!!
      enddo
      !
      !get Akw
      allocate(Akreal(Nktot,Lreal))
      Akreal=zero
      do ispin=1,Nspin
         do iorb=1,Norb
            Akreal = Akreal - dimag(Gkreal(:,ispin,ispin,iorb,iorb,:))/pi/Nspin/Norb
         enddo
      enddo
      !
      !print
      !
      call splot3d("Akw.dat",knumber,wr,Akreal(:,:))
      unit=free_unit()
      open(unit,file="Bands_topo.dat")
      do ik=1,Nktot
        write(unit,*)knumber(ik),Eigval_topo(:,ik)
      enddo
      close(unit)
      unit=free_unit()
      open(unit,file="Bands_bare.dat")
      do ik=1,Nktot
        write(unit,*)knumber(ik),Eigval_bare(:,ik)
      enddo
      close(unit)
      !
      !
      ! 
   end subroutine get_Akw




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
  !---------------------------------------------------------------------
  !PURPOSE: GET POLES ON THE REAL AXIS
  !---------------------------------------------------------------------
  subroutine get_poles()
    integer                                       :: i,j,ik,ix,iy,Nso,Nktot,Npts
    integer                                       :: iorb,jorb
    integer                                       :: isporb,jsporb
    integer                                       :: ispin,jspin
    integer                                       :: iso,unit
    real(8),dimension(Nspin*Norb)                 :: dzeta
    complex(8),allocatable,dimension(:,:,:)       :: Hk_bare,Hk_topo
    complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: zeta,fg,gdelta,fgk
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
    allocate(kpath(3,2))
    kpath(1,:)=[0.0,1.0]!G-e<-R
    kpath(2,:)=[0.0,0.0]!G
    kpath(3,:)=[1.0,0.0]!G+e->R
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
    allocate(kpoints(Nktot,2))
    call TB_build_kgrid(kpath,Nkpath,kpoints)
    if(allocated(Sigmamats))deallocate(Sigmamats)
    if(allocated(Sigmareal))deallocate(Sigmareal)
    allocate(Sigmamats(Nktot,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sigmareal(Nktot,Nspin,Nspin,Norb,Norb,Lreal))
    do ik=1,Nktot
      Sigmamats(ik,:,:,:,:,:)=periodize_sigma_mats(kpoints(ik,:))
      Sigmareal(ik,:,:,:,:,:)=periodize_sigma_real(kpoints(ik,:))
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
    write(LOGfile,*)"Solving for the poles..."
    maxNinterval=-1
    do ik=1,Nktot
       do i=1,Lreal
          zeta(:,:) = (wr(i)+xmu)*eye(Nso) - dreal(nn2so(Sigmareal(ik,:,:,:,:,i)))
          Den(i) = dreal((zeta(1,1) - Hk_bare(1,1,ik))*(zeta(2,2) - Hk_bare(2,2,ik))) - Hk_bare(1,2,ik)*Hk_bare(2,1,ik)
       enddo
       Xcsign(0)=0.d0
       count=0
       sign_old=sgn(Den(Lreal/2+1))
       do i=Lreal/2+1,Lreal
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
       ipoles(ik) = brentq(det_poles,0.d0,wr(Lreal))
       iweight(ik)= get_weight(hk_bare(:,:,ik)-nn2so(Sigmamats(ik,:,:,:,:,1)))
       call delete_finter(finter_func)
    enddo
    call splot("BHZpoles.ed",ipoles,iweight)
    do int=1,maxNinterval
       unit=free_unit()
       open(unit,file="BHZpoles_int"//reg(txtfy(int))//".ed")
       if(any((Mpoles(:,int)/=0.d0)))then
          do ik=1,Nktot
             if(Mpoles(ik,int)/=0.d0)write(unit,*)ik-1,Mpoles(ik,int),Mweight(ik,int)
          enddo
       endif
       close(unit)
    enddo
    !
  end subroutine get_poles

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

end program cdn_bhz_2d








