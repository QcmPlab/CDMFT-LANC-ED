program cdn_kagome
   USE CDMFT_ED !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nx,Ny,Nlso,iloop,Nb,Nkx,Nky,iw,iii,jjj,kkk
   integer,dimension(2):: recover
   logical                                                                :: converged
   real(8)                                                                :: ts,wmixing,observable
   !Bath:
   real(8),allocatable                                                    :: Bath(:),Bath_prev(:)
   !The local hybridization function:
   complex(8),allocatable                                                 :: Hloc(:,:)
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Gmats,Greal,Smats,Sreal,Weiss,Weiss_old
   character(len=16)                                                      :: finput
   real(8),allocatable                                                    :: wt(:)
   complex(8),allocatable                                                 :: wm(:),wr(:)
   complex(8),allocatable                                                 :: Hk(:,:,:),Smats_lso(:,:,:),G0(:,:,:),G0latt(:,:,:,:)
   complex(8),dimension(:,:,:,:,:,:),allocatable                          :: observable_matrix
   !SYMMETRIES TEST
   real(8),dimension(:),allocatable                                       :: lambdasym_vector
   complex(8),dimension(:,:,:,:,:,:,:),allocatable                        :: Hsym_basis
   !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
   integer                                                                :: comm
   integer                                                                :: rank
   integer                                                                :: mpi_size
   logical                                                                :: master,hermiticize
   complex(8)                                                             :: dummy

   !Init MPI: use of MPI overloaded functions in SciFor
   call init_MPI(comm,.true.)
   rank   = get_Rank_MPI(comm)
   master = get_Master_MPI(comm)

   !

   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(ts,"TS",finput,default=10.d0,comment="hopping parameter")
   call parse_input_variable(Nkx,"Nkx",finput,default=30,comment="Number of kx point for BZ integration")
   call parse_input_variable(Nky,"Nky",finput,default=30,comment="Number of ku point for BZ integration")
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
   if (Norb/=1) stop "You are using too many orbitals. Only 1 allowed"
   if (Nspin/=2) stop "You are using too many spin-orbitals. Only 2 allowed"
   Nlat=3
   Nlso=Nlat*Nspin*Norb
   if(.not.allocated(wm))allocate(wm(Lmats))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)

   !if(ED_VERBOSE > 0)call naming_convention()
   !Allocate Weiss Field:
   allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats_lso(Nlso,Nlso,Lmats))
   allocate(G0(Nlso,Nlso,Lmats))
   allocate(G0latt(Nlso,Nlso,Nkx*Nky,Lmats))

   !Build Hk and Hloc
   call generate_hk_hloc()
   ! Check whether k average gives local hamiltonian
   !do iii=1,12
   !    do jjj=1,12
   !        dummy = 1d0/(Nkx*Nky)*sum(hk(iii,jjj,:))
   !        dummy = dummy - hloc(iii,jjj)
   !        write(LOGFile,"(F0.0,A,F0.0,A)",advance='no') real(dummy), '+i' , aimag(dummy), ' '
   !    enddo
   !    write(LOGFILE,*) "//"
   !enddo

   !CUSTOM OBSERVABLES: n_tot, n_mag, Ekin
   allocate(observable_matrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
   call init_custom_observables(3,Hk)
   observable_matrix=zero
   do iii=1,Nlat
      observable_matrix(iii,iii,1,1,1,1)=one/Nlat
      observable_matrix(iii,iii,Nspin,Nspin,1,1)=one/Nlat
   enddo
   call add_custom_observable("n_tot",nnn2lso(observable_matrix))
   observable_matrix=zero
   do iii=1,Nlat
      observable_matrix(iii,iii,1,1,1,1)=one/Nlat
      observable_matrix(iii,iii,Nspin,Nspin,1,1)=-one/Nlat
   enddo
   call add_custom_observable("n_mag",nnn2lso(observable_matrix))
   call add_custom_observable("Ekin",Hk)

   !SETUP BATH STEP 1
   allocate(lambdasym_vector(8))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,8))
   !
   lambdasym_vector(1)=ts
   Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(bath_model(1.d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0))
   !
   lambdasym_vector(2)=0d0
   Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(bath_model(0.d0, 1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0))
   !
   lambdasym_vector(3)=0d0
   Hsym_basis(:,:,:,:,:,:,3)=lso2nnn(bath_model(0d0,  0d0, 1d0, 0d0, 0d0, 0d0, 0d0, 0d0))
   !
   lambdasym_vector(4)=0d0
   Hsym_basis(:,:,:,:,:,:,4)=lso2nnn(bath_model(0d0,  0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0))
   !
   lambdasym_vector(5)=0d0
   Hsym_basis(:,:,:,:,:,:,5)=lso2nnn(bath_model(0d0,  0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0))
   !
   lambdasym_vector(6)=0d0
   Hsym_basis(:,:,:,:,:,:,6)=lso2nnn(bath_model(0d0,  0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0))
   !
   lambdasym_vector(7)=0d0
   Hsym_basis(:,:,:,:,:,:,7)=lso2nnn(bath_model(0d0,  0d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0))
   !
   lambdasym_vector(8)=0d0
   Hsym_basis(:,:,:,:,:,:,8)=lso2nnn(bath_model(0d0,  0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0))
   !
   !SETUP BATH STEP 2 and SETUP SOLVER
   write(*,*) "setting hreplica"
   call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
   Nb=ed_get_bath_dimension(Hsym_basis)
   allocate(bath(Nb))
   allocate(bath_prev(Nb))
   write(*,*) "initializing solver"
   call ed_init_solver(comm,bath)
   !
   !build G0 and store it in Gmats
   G0 = zero
   do iii=1, Lmats
      do jjj=1, Nkx*Nky
          G0latt(:,:,jjj,iii) = (wm(iii)+xmu)*eye(Nlso) - Hk(:,:,jjj)
          call inv(G0latt(:,:,jjj,iii))
          G0(:,:,iii) = G0(:,:,iii) + G0latt(:,:,jjj,iii)
      enddo
      G0(:,:,iii) = 1./(Nkx*Nky) * G0(:,:,iii)
      Gmats(:,:,:,:,:,:,iii) = lso2nnn(G0(:,:,iii))
   enddo
   call print_G
   !Smats = zero
   write(*,*) "fitting G0"
   call ed_chi2_fitgf(Gmats,bath)
   !call dmft_self_consistency(Gmats,Smats,Weiss,lso2nnn(Hloc),cg_scheme)

   !!DMFT loop
   !iloop=0;converged=.false.
   !do while(.not.converged.AND.iloop<nloop)
      !iloop=iloop+1
      !if(master)call start_loop(iloop,nloop,"DMFT-loop")

      !!Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
      !call ed_solve(comm,bath,lso2nnn(hloc))
      !!!! GO ON HERE !!!
      !call ed_get_sigma_matsubara(Smats)
      !call ed_get_sigma_realaxis(Sreal)

      !!Compute the local gfs:
      !call dmft_gloc_matsubara(Hk,Gmats,Smats)
      !if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
      !!
      !!Get the Weiss field/Delta function to be fitted
      !call dmft_self_consistency(Gmats,Smats,Weiss,lso2nnn(Hloc),cg_scheme)
      !call Bcast_MPI(comm,Weiss)
      !!
      !!MIXING:
      !!if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_Old
      !!Weiss_old=Weiss
      !!
      !!Perform the SELF-CONSISTENCY by fitting the new bath
      !if(master)then
         !call ed_chi2_fitgf(Weiss,bath)
         !!
         !!MIXING:
         !!call adaptive_mix(Bath(Nbath+1:),Bath_fitted(Nbath+1:)-Bath(Nbath+1:),wmixing,iloop)
         !if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
         !Bath_Prev=Bath
         !!
         !!Check convergence (if required change chemical potential)
         !converged = check_convergence(Weiss(:,:,1,1,1,1,:),dmft_error,nsuccess,nloop)
      !endif
      !!
      !call Bcast_MPI(comm,bath)
      !call Bcast_MPI(comm,converged)
      !call Bcast_MPI(comm,xmu)
      !!
      !if(master)call end_loop
   !enddo

   !!Compute the local gfs:
   !call dmft_gloc_realaxis(Hk,Greal,Sreal)
   !if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)

   !!Compute the Kinetic Energy:
   !do iw=1,Lmats
      !Smats_lso(:,:,iw)=nnn2lso(Smats(:,:,:,:,:,:,iw))
   !enddo
   !call dmft_kinetic_energy(Hk(:,:,:),Smats_lso)


   call finalize_MPI()


contains

   !+------------------------------------------------------------------+
   !PURPOSE  : Hloc for the 2d Kagome model
   !+------------------------------------------------------------------+


   function hloc_model(ts_) result (hloc_)
      ! Directly created in lso basis
      integer                                               :: ispin, lowerbound, upperbound, Nlo
      real(8),intent(in)                                    :: ts_
      complex(8),dimension(Nlso,Nlso)                       :: hloc_
      !
      Nlo = Nlat*Norb
      hloc_  = zero
      !
      do ispin=1,Nspin
         lowerbound = (ispin-1)*Nlo+1
         upperbound = ispin*Nlo
         hloc_(lowerbound:upperbound,lowerbound:upperbound) = hloc_matrix(ts_)
      enddo
      !
   end function hloc_model

   function bath_model(tur_, tui_, tdr_, tdi_, e_, s1_, s2_, s3_) result (bath_)
      ! Directly created in lso basis
      integer                                               :: ispin, iorb, lowerbound, upperbound
      integer                                               :: Nlo
      real(8),intent(in)                                    :: tur_, tui_, tdr_, tdi_, e_, s1_, s2_, s3_
      complex(8),dimension(Nlso,Nlso)                       :: bath_
      !
      Nlo = Nlat*Norb
      bath_  = zero
      !
      do ispin=1,Nspin
         lowerbound = (ispin-1)*Nlo+1
         upperbound = ispin*Nlo
         bath_(lowerbound:upperbound,lowerbound:upperbound) = bath_matrix(tur_, tui_, tdr_, tdi_, e_, ispin)
      enddo
      !
      do iorb=1,Norb
         lowerbound = iorb
         upperbound = iorb+Nlo
         bath_(lowerbound:upperbound:Nlo,lowerbound:upperbound:Nlo) = bath_(lowerbound:upperbound:Nlo,lowerbound:upperbound:Nlo)&
                                                                      + bath_matrix_spin(s1_, s2_, s3_)
      enddo
      !
   end function bath_model

   function hk_model(kpoint,N) result(Hk_)
      ! Directly created in lso basis, kpoint MUST be in direct coordinates
      ! the N input is only needed or the TB_build_model routine....
      integer                         :: N, ispin, lowerbound, upperbound, Nlo
      real(8),dimension(:)            :: kpoint
      complex(8),dimension(N,N)       :: Hk_, hloc_, h1_, h2_, h3_, h4_, h5_, h6_
      !
      if (N/=Nlso) stop "dimensionality of hk wrong!"
      Nlo   = Nlat*Norb
      Hk_   = zero
      hloc_ = zero
      h1_   = zero
      h2_   = zero
      h3_   = zero
      h4_   = zero
      h5_   = zero
      h6_   = zero
      !
      do ispin=1,Nspin
         lowerbound = (ispin-1)*Nlo+1
         upperbound = ispin*Nlo
         hloc_(lowerbound:upperbound,lowerbound:upperbound) = hloc_matrix(ts)
         h1_(lowerbound:upperbound,lowerbound:upperbound)   = hhop1_matrix(ts)
         h2_(lowerbound:upperbound,lowerbound:upperbound)   = hhop2_matrix(ts)
         h3_(lowerbound:upperbound,lowerbound:upperbound)   = hhop3_matrix(ts)
         h4_(lowerbound:upperbound,lowerbound:upperbound)   = hhop4_matrix(ts)
         h5_(lowerbound:upperbound,lowerbound:upperbound)   = hhop5_matrix(ts)
         h6_(lowerbound:upperbound,lowerbound:upperbound)   = hhop6_matrix(ts)
      enddo
      !
      Hk_ = hloc_ + h1_*exp(2*pi*(0,1)*kpoint(1))              + h2_*exp(-2*pi*(0,1)*(kpoint(1))) &
                  + h3_*exp(2*pi*(0,1)*kpoint(2))              + h4_*exp(-2*pi*(0,1)*kpoint(2)) &
                  + h5_*exp(-2*pi*(0,1)*(kpoint(1)-kpoint(2))) + h6_*exp(2*pi*(0,1)*(kpoint(1)-kpoint(2)))
      !
   end function hk_model


   !AUXILLIARY MATRIX CONSTRUCTORS for spin Up block, spindown = spinup

   function hloc_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 1, 1 ]
      tempmat(2,:) = [ 1, 0, 1 ]
      tempmat(3,:) = [ 1, 1, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hloc_matrix

   function hhop1_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 1, 0 ]
      tempmat(2,:) = [ 0, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop1_matrix

   function hhop2_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0 ]
      tempmat(2,:) = [ 1, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop2_matrix

   function hhop3_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 1 ]
      tempmat(2,:) = [ 0, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop3_matrix

   function hhop4_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 0 ]
      tempmat(3,:) = [ 1, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop4_matrix

   function hhop5_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 1 ]
      tempmat(3,:) = [ 0, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop5_matrix

   function hhop6_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 0 ]
      tempmat(3,:) = [ 0, 1, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop6_matrix

   function bath_matrix(tur, tui, tdr, tdi, e, spinflav) result(bmat)
      complex(8),dimension(3,3) :: bmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: tur, tui, tdr, tdi, e
      integer,intent(in)        :: spinflav
      !
      bmat    = zero
      tempmat = zero
      !
      if(spinflav==1)then
          tempmat(1,:) = [ 0, 1, 1]
          tempmat(2,:) = [ 1, 0, 1]
          tempmat(3,:) = [ 1, 1, 0]
          bmat         = bmat + tur*tempmat + tui*(0,1)*tempmat
      elseif(spinflav==2)then
          tempmat(1,:) = [ 0, 1, 1]
          tempmat(2,:) = [ 1, 0, 1]
          tempmat(3,:) = [ 1, 1, 0]
          bmat         = bmat + tdr*tempmat + tdi*(0,1)*tempmat
      else
          stop "invalid spinflav intput in bath creation"
      endif
      !
      tempmat(1,:) = [ 1, 0, 0]
      tempmat(2,:) = [ 0, 1, 0]
      tempmat(3,:) = [ 0, 0, 1]
      bmat         = bmat + e*tempmat
      !
   end function bath_matrix

   function bath_matrix_spin(s1, s2, s3) result(bmat)
      complex(8),dimension(2,2) :: bmat
      real(8),dimension(2,2)    :: tempmat
      real(8),intent(in)        :: s1, s2, s3
      !
      bmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 1]
      tempmat(2,:) = [ 1, 0]
      bmat         = bmat + s1*tempmat
      !
      tempmat(1,:) = [ 0, -1]
      tempmat(2,:) = [ 1,  0]
      bmat         = bmat + s2*(0,1)*tempmat
      !
      tempmat(1,:) = [ 1,  0]
      tempmat(2,:) = [ 0, -1]
      bmat         = bmat + s3*tempmat
      !
   end function bath_matrix_spin

   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate Hloc and Hk
   !-------------------------------------------------------------------------------------------

   subroutine generate_hk_hloc()
      integer                      :: ik,i,j
      real(8),dimension(Nkx*Nky,2) :: kgrid
      real(8),dimension(Nkx)       :: gridx
      real(8),dimension(Nky)       :: gridy
      !
      !kmesh in direct coordinates
      gridx = linspace(0d0,1d0,Nkx,iend=.false.)
      gridy = linspace(0d0,1d0,Nky,iend=.false.)
      do i=1,Nkx
          do j=1,Nky
              ik=(i-1)*Nky+j
              kgrid(ik,1)=gridx(i)
              kgrid(ik,2)=gridy(j)
          enddo
      enddo
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
      hloc=hloc_model(ts)
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
   !
  subroutine print_G
    character(len=64) :: suffix
    integer           :: ilat,jlat,iorb,ispin,jorb
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb  
             do jorb=1,Norb
                do ispin=1,Nspin
                   suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                   call splot("G"//reg(suffix)//"_iw.ed",aimag(wm),Gmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
  end subroutine print_G

end program cdn_kagome
