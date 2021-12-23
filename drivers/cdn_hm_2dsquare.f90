program cdn_hm_2dsquare
   USE CDMFT_ED
   !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nx,Ny,Nso,Nlo,Nlso,iloop,Nb,Nkx,Nky,iw,iii,jjj,kkk,Ntr
   integer,dimension(2):: recover
   logical                                                                :: converged
   real(8)                                                                :: ts,wmixing,delta
   !Bath:
   real(8),allocatable                                                    :: bath(:),bath_prev(:)
   !The local hybridization function:
   complex(8),allocatable                                                 :: Hloc(:,:)
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Gmats,Greal,Smats,Sreal,Weiss
   character(len=16)                                                      :: finput
   complex(8),allocatable                                                 :: wm(:),wr(:)
   complex(8),allocatable                                                 :: Hk(:,:,:),Smats_lso(:,:,:)
   !Density matrices
   complex(8),allocatable,dimension(:,:)                                  :: cluster_density_matrix
   complex(8),allocatable,dimension(:,:)                                  :: reduced_density_matrix
   complex(8),allocatable,dimension(:,:)                                  :: small_dm, big_dm
   complex(8),allocatable,dimension(:,:)                                  :: local_density_matrix
   !SYMMETRIES TEST
   real(8),dimension(:),allocatable                                       :: lambdasym_vector
   complex(8),dimension(:,:,:,:,:,:,:),allocatable                        :: Hsym_basis
   !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
   integer                                                                :: comm
   integer                                                                :: rank
   integer                                                                :: mpi_size
   logical                                                                :: master
   logical                                                                :: periodize
   character(len=6)                                                       :: scheme
   real(8),dimension(:,:),allocatable                                     :: dens, dens_up, dens_dw, docc, mag
   !Init MPI: use of MPI overloaded functions in SciFor
   call init_MPI(comm,.true.)
   rank   = get_Rank_MPI(comm)
   master = get_Master_MPI(comm)

   !

   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputHM2D.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(ts,"TS",finput,default=1.d0,comment="hopping parameter")
   call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
   call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of cluster sites in y direction")
   call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
   call parse_input_variable(Nky,"Nky",finput,default=10,comment="Number of ky point for BZ integration")
   call parse_input_variable(periodize,"PERIODIZE",finput,default=.false.,comment="Periodization: T or F")
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
   !Ny=Nx:w

   !Nky=Nkx
   Nlat=Nx*Ny
   Nso=Nspin*Norb
   Nlo=Nlat*Norb
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
   allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats_lso(Nlso,Nlso,Lmats))

   allocate(dens(Nlat,Norb),dens_up(Nlat,Norb),dens_dw(Nlat,Norb),docc(Nlat,Norb),mag(Nlat,Norb))

   !Build Hk and Hloc
   call generate_hk_hloc()

   !Build Hsym_basis and lambdasym_vector
   allocate(lambdasym_vector(1))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,1))
   Hsym_basis(:,:,:,:,:,:,1)=abs(lso2nnn(Hloc))
   lambdasym_vector=[-1.d0] !not propto TS, since TS is contained in Hloc
   
   !SETUP BATH & SOLVER
   call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
   Nb=ed_get_bath_dimension(Hsym_basis)
   allocate(bath(Nb))
   allocate(bath_prev(Nb))
   bath_prev=zero
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

      !Retrieve ALL REDUCED DENSITY MATRICES DOWN TO THE LOCAL one
      if(master)then
         call ed_get_cluster_dm(cluster_density_matrix,doprint=.false.)
         big_dm = cluster_density_matrix
         do Ntr=1,Nlat-1 ! Ntr: number of sites we want to trace out
            call ed_get_reduced_dm(reduced_density_matrix,Nlat-Ntr,doprint=.true.)
            call subtrace(small_dm,big_dm,Nlat-Ntr+1,Nlat-Ntr)
            !>>>[Whatever you need to do with reduced_density_matrix...]<<<
            !e.g. CROSS-CHECK:
            print*,"*************************************************************************"
            print*,str(Nlat-Ntr)//"-site REDUCED DENSITY MATRIX"
            delta = maxval(abs(reduced_density_matrix-small_dm))
            print*,"direct and iterative subtracings match up to",delta
            print*,"*************************************************************************"
            deallocate(big_dm)
            big_dm = small_dm
         enddo
         !
         if(Norb==1)then
            ! BENCHMARK: build the local-dm according to Eq.4 in Mod.Phys.Lett.B.2013.27:05
            call ed_get_reduced_dm(local_density_matrix,1,doprint=.false.)
            call ed_get_mag(mag)
            call ed_get_dens(dens)
            call ed_get_docc(docc)
            dens_up = 0.5d0*(dens + mag)
            dens_dw = 0.5d0*(dens - mag)
            write(*,*)
            write(*,*) "LOCAL-DM BENCHMARK [Mod.Phys.Lett.B.2013.27:05]"
            write(*,*) "Semi-Analytical Estimate  |  Error"
            write(*,*) 1-dens_up(1,1)-dens_dw(1,1)+docc(1,1), "|", abs(1-dens_up(1,1)-dens_dw(1,1)+docc(1,1)-local_density_matrix(1,1))
            write(*,*) dens_up(1,1)-docc(1,1),                "|", abs(dens_up(1,1)-docc(1,1)-local_density_matrix(2,2))
            write(*,*) dens_dw(1,1)-docc(1,1),                "|", abs(dens_dw(1,1)-docc(1,1)-local_density_matrix(3,3))
            write(*,*) docc(1,1),                             "|", abs(docc(1,1)-local_density_matrix(4,4))
            write(*,*)
         endif
      endif


      !Compute the local gfs:
      call dmft_gloc_matsubara(Hk,Gmats,Smats)
      if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
      !
      !Get the Weiss field/Delta function to be fitted
      call dmft_self_consistency(Gmats,Smats,Weiss,lso2nnn(Hloc),cg_scheme)
      call Bcast_MPI(comm,Weiss)
      !
      !
      !Perform the SELF-CONSISTENCY by fitting the new bath
      if(master)then
         call ed_chi2_fitgf(Weiss,bath)
         !
         !MIXING:
         if(iloop>1)bath = wmixing*bath + (1.d0-wmixing)*bath_prev
         bath_prev=bath
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
   
   !PERIODIZE
   if(periodize)then
      call print_periodized([Nkx,Nky],hk_model,hk_periodized,scheme)
   endif

   call finalize_MPI()


contains


   !-------------------------------------------------------------------------------------------
   !PURPOSE:  Hk model for the 2d square lattice
   !-------------------------------------------------------------------------------------------
   
   function hloc_model(N) result (hloc)
      integer                                               :: ilat,jlat,ispin,iorb,ind1,ind2,N
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
      complex(8),dimension(N,N)                             :: hloc
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do iorb=1,Norb
            do ilat=1,Nx
               do jlat=1,Ny
                  ind1=indices2N([ilat,jlat])
                  hopping_matrix(ind1,ind1,ispin,ispin,iorb,iorb)= 0.d0!-mu_var
                  if(ilat<Nx)then
                     ind2=indices2N([ilat+1,jlat])
                     hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)= -ts
                  endif
                  if(ilat>1)then
                     ind2=indices2N([ilat-1,jlat])
                     hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)= -ts
                  endif
                  if(jlat<Ny)then
                     ind2=indices2N([ilat,jlat+1])
                     hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)= -ts
                  endif
                  if(jlat>1)then
                     ind2=indices2N([ilat,jlat-1])
                     hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)= -ts
                  endif
               enddo
            enddo
         enddo
      enddo
      !
      Hloc=nnn2lso(hopping_matrix)
      !
   end function hloc_model
   !
   !
   !
   function hk_model(kpoint,N) result(hk)
      integer                                               :: n,ilat,ispin,iorb,ind1,ind2
      real(8),dimension(:)                                  :: kpoint
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
      complex(8),dimension(N,N)                             :: hk
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do iorb=1,Norb
            do ilat=1,Nx
               ind1=indices2N([ilat,1])
               ind2=indices2N([ilat,Ny])
               hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)=hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb) -ts*exp(xi*kpoint(2)*Ny)
               hopping_matrix(ind2,ind1,ispin,ispin,iorb,iorb)=hopping_matrix(ind2,ind1,ispin,ispin,iorb,iorb) -ts*exp(-xi*kpoint(2)*Ny)
            enddo
            do ilat=1,Ny
               ind1=indices2N([1,ilat])
               ind2=indices2N([Nx,ilat])
               hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb)=hopping_matrix(ind1,ind2,ispin,ispin,iorb,iorb) -ts*exp(xi*kpoint(1)*Nx)
               hopping_matrix(ind2,ind1,ispin,ispin,iorb,iorb)=hopping_matrix(ind2,ind1,ispin,ispin,iorb,iorb) -ts*exp(-xi*kpoint(1)*Nx)
            enddo
         enddo
      enddo
      !
      Hk=nnn2lso(hopping_matrix)+hloc_model(N)
   end function hk_model
   !
   !
   !
   function hk_periodized(kpoint,N) result(Hk)
      real(8),dimension(:)                          :: kpoint
      integer                                       :: Nlat_,Nx_,Ny_,N
      complex(8),dimension(N,N)                     :: Hk
      !
      Nlat_=Nlat
      Nx_=Nx
      Ny_=Ny
      Nlat=1
      Nx=1
      Ny=1
      !
      Hk=hk_model(kpoint,Nspin*Norb)
      !
      Nlat=Nlat_
      Nx=Nx_
      Ny=Ny_
      !
   end function hk_periodized

   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate Hloc and Hk
   !-------------------------------------------------------------------------------------------

   subroutine generate_hk_hloc()
      integer                                     :: ik
      real(8),dimension(Nkx*Nky,2)                :: kgrid
      real(8),dimension(2)                        :: e1,e2,bk1,bk2
      real(8)                                     :: bklen
      !
      e1 = [1d0, 0d0]
      e2 = [0d0, 1d0]
      call TB_set_ei(eix=e1,eiy=e2)
      bklen=2d0*pi
      bk1=bklen*[1d0, 0d0]
      bk2=bklen*[0d0, 1d0]
      call TB_set_bk(bkx=bk1,bky=bk2)
      !
      call TB_build_kgrid([Nkx,Nky],kgrid)
      kgrid(:,1)=kgrid(:,1)/Nx
      kgrid(:,2)=kgrid(:,2)/Ny
      !
      if(allocated(hk))deallocate(hk)
      if(allocated(hloc))deallocate(Hloc)
      !
      allocate(Hk(Nlso,Nlso,Nkx*Nky),Hloc(Nlso,Nlso))
      hk=zero
      hloc=zero
      !
      call TB_build_model(Hk,hk_model,Nlso,kgrid)
      Hloc=hloc_model(Nlso)
      where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
      !
   end subroutine generate_hk_hloc

   !-------------------------------------------------------------------------------------------
   !PURPOSE: auxiliary reshape functions
   !-------------------------------------------------------------------------------------------

   ! These two functions are totally equivalent to:
   !
   ! > function lso2nnn_reshape(Hlso,Nlat,Nspin,Norb) result(Hnnn)
   !     --> from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]
   ! > function nnn2lso_reshape(Hnnn,Nlat,Nspin,Norb) result(Hlso)
   !     --> from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
   !
   ! which are included in CDMFT_ED (specifically ED_AUX_FUNX).
   !
   ! For some reason they are not retrieved ("no IMPLICIT type"),
   ! so for now we keep these local copies. 

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
   !
   !
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

   !-------------------------------------------------------------------------------------------
   !PURPOSE: managing verbose LOG files
   !-------------------------------------------------------------------------------------------

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



   
   !-------------------------------------------------------------------------------------------
   !PURPOSE: periodization routines
   !-------------------------------------------------------------------------------------------

   !include "auxiliary_routines.f90" 
   !--> NO, including here all this stuff is a nightmare and we don't really need to do so. 
   !    Let's cherrypick what we need instead.

   !--> USER-INTERFACE
   subroutine print_periodized(Nkpts,func1,func2,scheme)
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)  :: gmats_periodized, Smats_periodized, Gloc_per_iw, Sloc_per_iw
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)  :: greal_periodized, Sreal_periodized, Gloc_per_rw, Sloc_per_rw
      integer,dimension(:)                               :: nkpts
      real(8),dimension(product(Nkpts),size(Nkpts))      :: kgrid
      integer                                            :: i
      character(len=6)                                   :: scheme
      interface
      function func1(x,N)
         integer                   :: N
         real(8),dimension(:)      :: x
         complex(8),dimension(N,N) :: func1
      end function func1
      end interface
      interface
      function func2(x,N)
         integer                   :: N
         real(8),dimension(:)      :: x
         complex(8),dimension(N,N) :: func2
      end function func2
      end interface
      !
      Gloc_per_iw=zero
      Sloc_per_iw=zero
      Gloc_per_rw=zero
      Sloc_per_rw=zero
      !
      call TB_build_kgrid(nkpts,kgrid)
      do i=1,size(Nkpts)
         kgrid(:,i)=kgrid(:,i)/Nkpts(i)
      enddo
      !
      if(ED_VERBOSE .ge. 1)write(LOGfile,*)"Computing periodized quantities using    ",scheme," scheme"
      !
      if(MASTER)call start_timer
      do i=1,product(Nkpts)
         if(scheme=="g") then
            call build_sigma_g_scheme(kgrid(i,:),gmats_periodized,greal_periodized,smats_periodized,&
               sreal_periodized,func1(kgrid(i,:),Nlat*Nspin*Norb),func2(kgrid(i,:),Nspin*Norb))
         elseif(scheme=="sigma")then
            call build_g_sigma_scheme(kgrid(i,:),gmats_periodized,greal_periodized,smats_periodized,&
               sreal_periodized,func2(kgrid(i,:),Nspin*Norb))
         else
            STOP "Nonexistent periodization scheme"
         endif
         Gloc_per_iw=Gloc_per_iw+(gmats_periodized/product(Nkpts))
         Sloc_per_iw=Sloc_per_iw+(smats_periodized/product(Nkpts))
         Gloc_per_rw=Gloc_per_rw+(greal_periodized/product(Nkpts))
         Sloc_per_rw=Sloc_per_rw+(sreal_periodized/product(Nkpts))
         if(ED_VERBOSE .ge. 1)call eta(i,product(Nkpts))
      enddo 
      if(MASTER)call stop_timer
      !
      if(master)call dmft_print_gf_matsubara(gloc_per_iw,"Gloc_periodized",iprint=4)
      if(master)call dmft_print_gf_matsubara(sloc_per_iw,"Sigma_periodized",iprint=4)
      !
      if(master)call dmft_print_gf_realaxis(gloc_per_rw,"Gloc_periodized",iprint=4)
      if(master)call dmft_print_gf_realaxis(sloc_per_rw,"Sigma_periodized",iprint=4)
      !
   end subroutine print_periodized
 
   

   !--> SIGMA-SCHEME
   subroutine periodize_sigma_scheme(kpoint,smats_periodized,sreal_periodized)
      integer                                                     :: ilat,jlat,ispin,iorb,ii
      real(8),dimension(:)                                        :: kpoint
      integer,dimension(:),allocatable                            :: ind1,ind2
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: smats_periodized
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: sreal_periodized
      !
      !
      if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
      if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
      !
      smats_periodized=zero
      sreal_periodized=zero
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
      deallocate(ind1,ind2)
      !   
   end subroutine periodize_sigma_scheme
   !
   !
   !
   subroutine build_g_sigma_scheme(kpoint,gmats_periodized,greal_periodized,smats_periodized,sreal_periodized,Hk_per)
      integer                                                     :: i,ispin,iorb,ii
      real(8),dimension(:)                                        :: kpoint
      complex(8),dimension(Nspin*Norb,Nspin*Norb)                 :: Hk_per,tmpmat
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: gmats_periodized, Smats_periodized
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: greal_periodized, Sreal_periodized
      !
      !
      !Get Gimp^-1
      call periodize_sigma_scheme(kpoint,smats_periodized,sreal_periodized)
      !
      do ii=1,Lmats
         tmpmat=(xi*wm(ii)+xmu)*eye(Nspin*Norb) - hk_per - nn2so(Smats_periodized(:,:,:,:,ii))
         call inv(tmpmat)
         gmats_periodized(:,:,:,:,ii)=so2nn(tmpmat)
      enddo
      !
      do ii=1,Lreal
         tmpmat=(wr(ii)+xmu)*eye(Nspin*Norb) - hk_per - nn2so(Sreal_periodized(:,:,:,:,ii))
         call inv(tmpmat)
         greal_periodized(:,:,:,:,ii)=so2nn(tmpmat)
      enddo
      !
   end subroutine build_g_sigma_scheme


   
   !--> G-SCHEME
   subroutine periodize_g_scheme(kpoint,gmats_periodized,greal_periodized,hk_unper)
      integer                                                     :: ilat,jlat,ispin,iorb,ii
      real(8),dimension(:)                                        :: kpoint
      integer,dimension(:),allocatable                            :: ind1,ind2
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: tmpmat,Hk_unper
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: gmats_unperiodized![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal) :: greal_unperiodized ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: gmats_periodized
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: greal_periodized
      !
      !
      if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
      if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
      !
      gmats_unperiodized=zero
      greal_unperiodized=zero
      gmats_periodized=zero
      greal_periodized=zero
      tmpmat=zero
      !
      !
      do ii=1,Lmats
         tmpmat=(xi*wm(ii)+xmu)*eye(Nlat*Nspin*Norb) - hk_unper - nnn2lso(Smats(:,:,:,:,:,:,ii))
         call inv(tmpmat)
         gmats_unperiodized(:,:,:,:,:,:,ii)=lso2nnn(tmpmat)
      enddo
      !
      do ii=1,Lreal
         tmpmat=(wr(ii)+xmu)*eye(Nlat*Nspin*Norb) - hk_unper - nnn2lso(Sreal(:,:,:,:,:,:,ii))
         call inv(tmpmat)
         greal_unperiodized(:,:,:,:,:,:,ii)=lso2nnn(tmpmat)
      enddo
      !
      do ii=1,Lmats
         do ilat=1,Nlat
            ind1=N2indices(ilat)        
            do jlat=1,Nlat
               ind2=N2indices(jlat)
               gmats_periodized(:,:,:,:,ii)=gmats_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*gmats_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
            enddo
         enddo
      enddo
      !
      do ii=1,Lreal   
         do ilat=1,Nlat
            ind1=N2indices(ilat)        
            do jlat=1,Nlat
               ind2=N2indices(jlat)
               greal_periodized(:,:,:,:,ii)=greal_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*greal_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
            enddo
         enddo
      enddo
      !
      deallocate(ind1,ind2)
      !   
   end subroutine periodize_g_scheme
   !
   !
   ! 
   subroutine build_sigma_g_scheme(kpoint,gmats_periodized,greal_periodized,smats_periodized,sreal_periodized,Hk_unper,Hk_per)
      integer                                                     :: i,ispin,iorb,ii
      real(8),dimension(:)                                        :: kpoint
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: Hk_unper
      complex(8),dimension(Nspin*Norb,Nspin*Norb)                 :: Hk_per
      complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats)           :: invG0mats,invGmats
      complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal)           :: invG0real,invGreal
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: gmats_periodized, Smats_periodized
      complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: greal_periodized, Sreal_periodized
      !
      !
      invG0mats = zero
      invGmats  = zero
      invG0real = zero
      invGreal  = zero
      Smats_periodized  = zero
      Sreal_periodized  = zero
      !
      !Get G0^-1
      do ii=1,Lmats
         invG0mats(:,:,ii) = (xi*wm(ii)+xmu)*eye(Nspin*Norb)  - Hk_per            
      enddo
      do ii=1,Lreal
         invG0real(:,:,ii) = (wr(ii)+xmu)*eye(Nspin*Norb)   - Hk_per   
      enddo
      !
      !Get Gimp^-1
      call periodize_g_scheme(kpoint,gmats_periodized,greal_periodized,Hk_unper)
      !
      do ii=1,Lmats
         invGmats(:,:,ii) = nn2so(gmats_periodized(:,:,:,:,ii))
         call inv(invGmats(:,:,ii))
      enddo
      do ii=1,Lreal
         invGreal(:,:,ii) = nn2so(greal_periodized(:,:,:,:,ii))
         call inv(invGreal(:,:,ii))
      enddo
      !
      !Get Sigma functions: Sigma= G0^-1 - G^-1
      Smats_periodized=zero
      Sreal_periodized=zero
      !
      do ii=1,Lmats
         Smats_periodized(:,:,:,:,ii) = so2nn(invG0mats(:,:,ii) - invGmats(:,:,ii))
      enddo
      do ii=1,Lreal
         Sreal_periodized(:,:,:,:,ii) = so2nn(invG0real(:,:,ii) - invGreal(:,:,ii))
      enddo
      !
      !
   end subroutine build_sigma_g_scheme
    

   !--> AUXILIARY-TO-PERIODIZATION
   function indices2N(indices) result(N)
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      !
      N=Nx*(indices(2)-1)+indices(1)
   end function indices2N
   !
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
   !
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


    !==> TO BE REMOVED, here just for testing 
    !+---------------------------------------------------------------------+
    !PURPOSE :  reduce a generic dm by tracing out a given number of sites
    !+---------------------------------------------------------------------+
    subroutine subtrace(red_dm,big_dm,bigNsites,redNsites)
      complex(8),dimension(:,:),allocatable,intent(out) :: red_dm
      complex(8),dimension(:,:),allocatable,intent(in)  :: big_dm
      integer                              ,intent(in)  :: bigNsites
      integer                              ,intent(in)  :: redNsites
      integer         :: i,j,io,jo,iUP,iDW,jUP,jDW
      integer         :: iIMPup,iIMPdw,jIMPup,jIMPdw
      integer         :: iREDup,iREDdw,jREDup,jREDdw
      integer         :: iTrUP,iTrDW,jTrUP,jTrDW
      integer         :: Nbig,Nred,dimBIG,dimRED
      !
      Nbig=Norb*bigNsites
      Nred=Norb*redNsites
      !
      dimBIG = 4**Nbig
      dimRED = 4**Nred
      !
      if(size(big_dm(:,1))/=dimBIG)stop "ERROR: Nsites is not consistent with the given big_dm"
      !
      allocate(red_dm(dimRED,dimRED)); red_dm=0.d0
      !
      do iUP = 1,2**Nbig
         do iDW = 1,2**Nbig
              i = iUP + (iDW-1)*2**Nbig
              iIMPup = iup-1
              iIMPdw = idw-1
              iREDup = Ibits(iIMPup,0,Nred)
              iREDdw = Ibits(iIMPdw,0,Nred)
              iTrUP  = Ibits(iIMPup,Nred,Nbig)
              iTrDW  = Ibits(iIMPdw,Nred,Nbig)
              do jUP = 1,2**Nbig
                 do jDW = 1,2**Nbig
                       j = jUP + (jDW-1)*2**Nbig
                       jIMPup = jup-1
                       jIMPdw = jdw-1
                       jREDup = Ibits(jIMPup,0,Nred)
                       jREDdw = Ibits(jIMPdw,0,Nred)
                       jTrUP  = Ibits(jIMPup,Nred,Nbig)
                       jTrDW  = Ibits(jIMPdw,Nred,Nbig)
                       if(jTrUP/=iTrUP.or.jTrDW/=iTrDW)cycle
                       io = (iREDup+1) + iREDdw*2**Nred
                       jo = (jREDup+1) + jREDdw*2**Nred
                       red_dm(io,jo) = red_dm(io,jo) + big_dm(i,j)
                 enddo
              enddo
         enddo
      enddo
      !
    end subroutine subtrace
    !<== TO BE REMOVED, here just for testing


end program cdn_hm_2dsquare











