program cdn_ssh_postprocessing
   USE CDMFT_ED
   !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   !USE MPI
   !
   implicit none
   integer                                                                :: Nk,Nlso,Nso,iloop,Nb,Nkx,Nky,iw,iii,jjj,Nkpath,Ndimer
   integer                                                                :: il,jl,is,js,io,jo
   integer,dimension(2):: recover
   logical                                                                :: converged
   real(8)                                                                :: vhop,whop,wmixing,observable,repart,impart,energy_offset
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

   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(vhop,"vhop",finput,default=0.25d0,comment="intra-dimer hopping")
   call parse_input_variable(whop,"whop",finput,default=0.25d0,comment="inter-dimer hopping")
   call parse_input_variable(Ndimer,"Ndimer",finput,default=1,comment="number of dimers")
   call parse_input_variable(Nk,"Nk",finput,default=10,comment="Number of k point for BZ integration")
   call parse_input_variable(Nkpath,"Nkpath",finput,default=100,comment="Number of k point for BZ path")
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

   !set global variables
   !if (Nspin/=1.or.Norb/=1) stop "You are using too many spin-orbitals"
   !
   Nlat=2*Ndimer
   Nlso=Nlat*Nspin*Norb
   Nso=Nspin*Norb
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   !Allocate Weiss Field:
   allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))

   !Build Hk and Hloc
   allocate(lambdasym_vector(3))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,3))
   !
   lambdasym_vector(1)=vhop
   Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(hloc_model(Nlso,1.d0,0.d0),Nlat,Nspin,Norb)
   !
   lambdasym_vector(2)=whop
   Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(hloc_model(Nlso,0.d0,1.d0),Nlat,Nspin,Norb)
   !
   lambdasym_vector(3)=energy_offset
   Hsym_basis(:,:,:,:,:,:,3)=lso2nnn(zeye(Nlso),Nlat,Nspin,Norb)
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
   call ed_get_gimp_realaxis(Greal)
   call ed_get_sigma_realaxis(Sreal)
   !
   !
   !RETRIEVE AND PERIODIZE
   !
   call   get_det_G()
   call   get_local_Sigma()
   call   get_local_g()



contains


   !+------------------------------------------------------------------+
   !PURPOSE  : Hloc for the 2d BHZ model
   !+------------------------------------------------------------------+


   function Hloc_model(N,vhop_,whop_) result (H0)
      integer                                               :: N,ispin,idimer
      real(8)                                               :: vhop_,whop_
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
      complex(8),dimension(N,N)                             :: H0
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do idimer=1,Ndimer
          hopping_matrix(2*idimer-1,2*idimer  ,ispin,ispin,:,:) = -vhop_
          hopping_matrix(2*idimer  ,2*idimer-1,ispin,ispin,:,:) = -vhop_
          !
          if(idimer < Ndimer) then
             hopping_matrix(2*idimer  ,2*idimer+1,ispin,ispin,:,:) = -whop_
             hopping_matrix(2*idimer+1,2*idimer  ,ispin,ispin,:,:) = -whop_
          endif
          !
          if(idimer > Ndimer) then
             hopping_matrix(2*idimer-1,2*idimer-2,ispin,ispin,:,:) = -whop_
             hopping_matrix(2*idimer-2,2*idimer-1,ispin,ispin,:,:) = -whop_
          endif
         enddo
      enddo
      !
      H0=nnn2lso(hopping_matrix,Nlat,Nspin,Norb)
      !
   end function hloc_model


   function hk_model(kpoint,N) result(Hk)
      integer                                                      :: Nispin,ispin,N
      real(8),dimension(:)                                         :: kpoint
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
      complex(8),dimension(N,N)                                    :: hk
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
            hopping_matrix(1,Nlat,ispin,ispin,:,:) = hopping_matrix(1,Nlat,ispin,ispin,:,:) - whop*exp(-xi*kpoint(1)*Ndimer)
            hopping_matrix(Nlat,1,ispin,ispin,:,:) = hopping_matrix(Nlat,1,ispin,ispin,:,:) - whop*exp( xi*kpoint(1)*Ndimer)
      enddo
      !
      Hk=nnn2lso(hopping_matrix,Nlat,Nspin,Norb)+hloc_model(N,vhop,whop)
      !
   end function hk_model

   subroutine generate_hk_hloc()
      integer                                     :: ik
      real(8),dimension(Nk,1)                     :: kgrid
      real(8),dimension(2)                        :: e1,e2,bk1,bk2
      real(8)                                     :: bklen
      !
      e1 = [1d0, 0d0]
      call TB_set_ei(eix=e1)
      bklen=2d0*pi
      bk1=bklen*[1d0, 0d0]
      call TB_set_bk(bkx=bk1)
      !
      call TB_build_kgrid([Nk],kgrid)
      kgrid(:,1)=kgrid(:,1)/Ndimer
      !
      if(allocated(hk))deallocate(hk)
      if(allocated(hloc))deallocate(Hloc)
      !
      allocate(Hk(Nlso,Nlso,Nk),Hloc(Nlso,Nlso))
      hk=zero
      hloc=zero
      !
      call TB_build_model(Hk,hk_model,Nlso,kgrid)
      Hloc=hloc_model(Nlso,vhop,whop)
      !
   end subroutine generate_hk_hloc


   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate the Hamiltonian with the minimal unit cell
   !-------------------------------------------------------------------------------------------

   function hk_periodized(kvec,N) result(hk)
    integer                   :: N,ii,Nlat_tmp,Ndimer_tmp
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    !
    Nlat_tmp=Nlat
    Ndimer_tmp=Ndimer
    Nlat=2
    Ndimer=1
    !
    Hk = hk_model(kvec,Nlat*Nso)
    !
    Nlat=Nlat_tmp
    Ndimer=Ndimer_tmp
    !
   end function hk_periodized
   !
   !
   !
   !
   !-------------------------------------------------------------------------------------------
   !PURPOSE: periodization M scheme
   !-------------------------------------------------------------------------------------------
   function periodize_sigma_Mscheme_real(kpoint,iw) result(s_lso)
      integer                                                     :: ilat,jlat,ispin,iorb,jorb,iw
      real(8),dimension(:)                                        :: kpoint
      integer,dimension(2)                                        :: indices1,indices2
      integer                                                     :: idimer_1,idimer_2,isite_1,isite_2
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: tmpmat
      complex(8),dimension(2,2,Nspin,Nspin,Norb,Norb)             :: sreal_periodized,greal_periodized
      complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb)             :: g_lso,s_lso
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)       :: greal_unperiodized
      !
      real(8)                                                     :: vhop_tmp, whop_tmp
      !
      vhop_tmp=vhop
      whop_tmp=whop
      !
      !
      sreal_periodized=zero
      greal_periodized=zero
      !
      vhop=0d0
      whop=0d0
      !
      tmpmat=(dcmplx(wr(iw),eps)+xmu)*eye(Nlat*Nspin*Norb) - hk_model(kpoint,Nlat*Nspin*Norb) - nnn2lso(Sreal(:,:,:,:,:,:,iw),Nlat,Nspin,Norb)
      call inv(tmpmat)
      greal_unperiodized(:,:,:,:,:,:)=lso2nnn(tmpmat,Nlat,Nspin,Norb)
      !
      vhop=vhop_tmp
      whop=whop_tmp
      !
      do ilat=1,Nlat
         indices1=N2indices(ilat)
         idimer_1=indices1(1)
         isite_1=indices1(2)        
         do jlat=1,Nlat
            indices2=N2indices(jlat)
            idimer_2=indices2(1)
            isite_2=indices2(2)   
            greal_periodized(isite_1,isite_2,:,:,:,:)=greal_periodized(isite_1,isite_2,:,:,:,:)+&
                                                       exp(-xi*kpoint(1)*(idimer_1-idimer_2))*greal_unperiodized(ilat,jlat,:,:,:,:)/Ndimer
         enddo
      enddo
      !
      !
      g_lso(:,:)=nnn2lso(greal_periodized(:,:,:,:,:,:),2,Nspin,Norb)
      call inv(g_lso(:,:))
      g_lso(:,:)=g_lso(:,:)-Hk_periodized(kpoint,2*Nspin*Norb)
      s_lso(:,:)=(dcmplx(wr(iw),eps)+xmu)*eye(2*Nspin*Norb)-Hk_periodized(kpoint,2*Nspin*Norb)-g_lso(:,:)
      call inv(g_lso(:,:))
      !   
   end function periodize_sigma_Mscheme_real

   function periodize_G_Mscheme_real(kpoint,iw) result(g_lso)
      integer                                                     :: ilat,jlat,ispin,iorb,jorb,iw
      real(8),dimension(:)                                        :: kpoint
      integer,dimension(2)                                        :: indices1,indices2
      integer                                                     :: idimer_1,idimer_2,isite_1,isite_2
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: tmpmat
      complex(8),dimension(2,2,Nspin,Nspin,Norb,Norb)             :: sreal_periodized,greal_periodized
      complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb)             :: g_lso
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)       :: greal_unperiodized
      !
      real(8)                                                     :: vhop_tmp, whop_tmp
      !
      vhop_tmp=vhop
      whop_tmp=whop
      !
      !
      sreal_periodized=zero
      greal_periodized=zero
      !
      vhop=0d0
      whop=0d0
      !
      tmpmat=(dcmplx(wr(iw),eps)+xmu)*eye(Nlat*Nspin*Norb) - hk_model(kpoint,Nlat*Nspin*Norb) - nnn2lso(Sreal(:,:,:,:,:,:,iw),Nlat,Nspin,Norb)
      call inv(tmpmat)
      greal_unperiodized(:,:,:,:,:,:)=lso2nnn(tmpmat,Nlat,Nspin,Norb)
      !
      vhop=vhop_tmp
      whop=whop_tmp
      !
      do ilat=1,Nlat
         indices1=N2indices(ilat)
         idimer_1=indices1(1)
         isite_1=indices1(2)        
         do jlat=1,Nlat
            indices2=N2indices(jlat)
            idimer_2=indices2(1)
            isite_2=indices2(2)   
            greal_periodized(isite_1,isite_2,:,:,:,:)=greal_periodized(isite_1,isite_2,:,:,:,:)+&
                                                       exp(-xi*kpoint(1)*(idimer_1-idimer_2))*greal_unperiodized(ilat,jlat,:,:,:,:)/Ndimer
         enddo
      enddo
      !
      !
      g_lso(:,:)=nnn2lso(greal_periodized(:,:,:,:,:,:),2,Nspin,Norb)
      !   
   end function periodize_G_Mscheme_real   
   
   
   !
   !+------------------------------------------------------------------+
   !PURPOSE  : Auxilliary reshape functions
   !+------------------------------------------------------------------+

   function lso2nnn(Hlso,Nlat_,Nspin_,Norb_) result(Hnnn)
      integer                                                     :: ilat,jlat,Nlat_,Nspin_,Norb_
      complex(8),dimension(Nlat_*Nspin_*Norb_,Nlat_*Nspin_*Norb_) :: Hlso
      complex(8),dimension(Nlat_,Nlat_,Nspin_,Nspin_,Norb_,Norb_) :: Hnnn
      integer                                                     :: iorb,jorb
      integer                                                     :: ispin,jspin
      integer                                                     :: is,js
      Hnnn=zero
      do ilat=1,Nlat_
         do jlat=1,Nlat_
            do ispin=1,Nspin_
               do jspin=1,Nspin_
                  do iorb=1,Norb_
                     do jorb=1,Norb_
                        is = iorb + (ilat-1)*Norb_ + (ispin-1)*Norb_*Nlat_
                        js = jorb + (jlat-1)*Norb_ + (jspin-1)*Norb_*Nlat_
                        Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function lso2nnn


   function nnn2lso(Hnnn,Nlat_,Nspin_,Norb_) result(Hlso)
      integer                                                     :: ilat,jlat,Nlat_,Nspin_,Norb_
      complex(8),dimension(Nlat_,Nlat_,Nspin_,Nspin_,Norb_,Norb_) :: Hnnn
      complex(8),dimension(Nlat_*Nspin_*Norb_,Nlat_*Nspin_*Norb_) :: Hlso
      integer                                                     :: iorb,jorb
      integer                                                     :: ispin,jspin
      integer                                                     :: is,js
      Hlso=zero
      do ilat=1,Nlat_
         do jlat=1,Nlat_
            do ispin=1,Nspin_
               do jspin=1,Nspin_
                  do iorb=1,Norb_
                     do jorb=1,Norb_
                        is = iorb + (ilat-1)*Norb_ + (ispin-1)*Norb_*Nlat_
                        js = jorb + (jlat-1)*Norb_ + (jspin-1)*Norb_*Nlat_
                        Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function nnn2lso


   function indices2N(indices) result(N)   !idimer, index in dimer -> ilat
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      !
      N=2*(indices(1)-1)+indices(2)
   end function indices2N

   function N2indices(N) result(indices)    !ilat -> idimer, index in dimer
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      if(mod(N,2)==0)then
         indices(1)=N/2
         indices(2)=2
      else
         indices(1)=(N+1)/2
         indices(2)=1
      endif
   end function N2indices

   !-------------------------------------------------------------------------------------------
   ! PURPOSE:Get determinant of G(k,w)
   !-------------------------------------------------------------------------------------------   
  
  subroutine get_det_G()
    integer                                       :: Lk,Nso,Niso,Npts
    integer                                       :: ik,iw
    integer                                       :: ispin,jspin
    complex(8),dimension(:,:),allocatable         :: Sigmareal
    real(8),dimension(:,:),allocatable            :: Akreal,kpoints
    complex(8),dimension(:,:),allocatable         :: Gkreal
    complex(8),allocatable,dimension(:,:,:)       :: Hk_bare
    real(8),dimension(:,:),allocatable            :: Kpath
    character(len=30)                             :: suffix
    !
    !
    !
    allocate(kpath(3,1))
    !
    kpath(1,:)=[-1.0]!G-e<-R
    kpath(2,:)=[0.0]!G
    kpath(3,:)=[1.0]!G+e->R
    !
    kpath=kpath*pi
    Npts  = size(kpath,1)
    Lk = (Npts-1)*Nkpath
    !   
    Nso=Nspin*Norb
    Niso=2*Nspin*Norb
    !
    !
    if(allocated(Hk_bare))deallocate(Hk_bare)
    allocate(Hk_bare(Niso,Niso,Lk));Hk_bare=zero
    !
    !
    call TB_build_model(hk_bare,hk_periodized,Niso,kpath,Nkpath)
    !
    !    
    allocate(kpoints(Lk,1))
    call TB_build_kgrid(kpath,Nkpath,kpoints)
    !
    !
    if(allocated(Sigmareal))deallocate(Sigmareal)
    if(allocated(Gkreal))deallocate(Gkreal)
    allocate(Sigmareal(Niso,Niso))
    allocate(Gkreal(Niso,Niso))
    Gkreal=zero      
    Sigmareal=zero
    allocate(Akreal(Lk,Lreal));Akreal=zero
    !
    !
    do ik=1,Lk
     do iw=1,Lreal
      Sigmareal=periodize_sigma_Mscheme_real(kpoints(ik,:),iw)
      Gkreal(:,:)=(dcmplx(wr(iw),0d0)+xmu)*eye(Niso) - Hk_bare(:,:,ik) - Sigmareal(:,:)
      call inv(Gkreal(:,:))
      Akreal(ik,iw) = log(abs(det(Gkreal(:,:)))/pi/Niso)
     enddo
    enddo
    !
    call splot3d("det_G_real_nso.dat",kpoints(:,1),wr,Akreal) 
    !
  end subroutine get_det_G


   !-------------------------------------------------------------------------------------------
   ! PURPOSE:Get determinant of G(k,w)
   !-------------------------------------------------------------------------------------------   
  
   subroutine get_local_sigma()
      integer                                                :: ik,iw
      real(8),dimension(Nk,1)                                :: kgrid
      real(8),dimension(2)                                   :: e1,e2,bk1,bk2
      real(8)                                                :: bklen
      complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb,Lreal)  :: s_lso
      complex(8),dimension(2,2,Nspin,Nspin,Norb,Norb,Lreal)  :: s_nnn
      !
      e1 = [1d0, 0d0]
      call TB_set_ei(eix=e1)
      bklen=2d0*pi
      bk1=bklen*[1d0, 0d0]
      call TB_set_bk(bkx=bk1)
      !
      call TB_build_kgrid([Nk],kgrid)
      s_lso=zero
      s_nnn=zero
      !
      do iw=1,Lreal
        do ik=1,Nk
          s_lso(:,:,iw)=s_lso(:,:,iw) + periodize_sigma_Mscheme_real(kgrid(ik,:),iw)/Nk
        enddo
        s_nnn(:,:,:,:,:,:,iw)=lso2nnn(s_lso(:,:,iw),2,Nspin,Norb)
      enddo
      call dmft_print_gf_realaxis(s_nnn,"perSigma",iprint=4)   
      !
   end subroutine get_local_sigma


   subroutine get_local_g()
      integer                                                :: ik,iw
      real(8),dimension(Nk,1)                                :: kgrid
      real(8),dimension(2)                                   :: e1,e2,bk1,bk2
      real(8)                                                :: bklen
      complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb,Lreal)  :: G_lso
      complex(8),dimension(2,2,Nspin,Nspin,Norb,Norb,Lreal)  :: G_nnn
      !
      e1 = [1d0, 0d0]
      call TB_set_ei(eix=e1)
      bklen=2d0*pi
      bk1=bklen*[1d0, 0d0]
      call TB_set_bk(bkx=bk1)
      !
      call TB_build_kgrid([Nk],kgrid)
      G_lso=zero
      G_nnn=zero
      !
      do iw=1,Lreal
        do ik=1,Nk
          G_lso(:,:,iw)=G_lso(:,:,iw) + periodize_G_Mscheme_real(kgrid(ik,:),iw)/Nk
        enddo
        G_nnn(:,:,:,:,:,:,iw)=lso2nnn(G_lso(:,:,iw),2,Nspin,Norb)
      enddo
      call dmft_print_gf_realaxis(G_nnn,"perG",iprint=4)   
      !
   end subroutine get_local_g


end program cdn_ssh_postprocessing








