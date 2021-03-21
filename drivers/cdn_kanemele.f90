program cdn_kanemele
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
   real(8)                                                                :: ts,Mh,lambda,wmixing,observable
   !Bath:
   real(8),allocatable                                                    :: Bath(:),Bath_prev(:)
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

   !Init MPI: use of MPI overloaded functions in SciFor
   call init_MPI(comm,.true.)
   rank   = get_Rank_MPI(comm)
   master = get_Master_MPI(comm)

   !

   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(ts,"TS",finput,default=10.d0,comment="hopping parameter")
   call parse_input_variable(Mh,"Mh",finput,default=0.d0,comment="crystal field splitting")
   call parse_input_variable(lambda,"lambda",finput,default=0.3d0,comment="spin-orbit coupling")
   call parse_input_variable(Nkx,"Nkx",finput,default=30,comment="Number of kx point for BZ integration")
   call parse_input_variable(Nky,"Nky",finput,default=30,comment="Number of ku point for BZ integration")
   call parse_input_variable(hermiticize,"HERMITICIZE",finput,default=.true.,comment="are bath replicas hermitian")
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
   Nlat=6
   Nlso=Nlat*Nspin*Norb
   if(.not.allocated(wm))allocate(wm(Lmats))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)

   !if(ED_VERBOSE > 0)call naming_convention()
   !Allocate Weiss Field:
   allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats_lso(Nlso,Nlso,Lmats))

   !Build Hk and Hloc
   call generate_hk_hloc()

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
   allocate(lambdasym_vector(3))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,3))
   !
   lambdasym_vector(1)=Mh
   Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(hloc_model(1.d0,0.d0,0.d0))
   !
   lambdasym_vector(2)=ts
   Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(hloc_model(0.d0,1.d0,0.d0))
   !
   lambdasym_vector(3)=lambda
   Hsym_basis(:,:,:,:,:,:,3)=lso2nnn(hloc_model(0.d0,0.d0,1.d0))
   !
   !SETUP BATH STEP 2 and SETUP SOLVER
   call set_Hloc(Hsym_basis,lambdasym_vector)
   Nb=get_bath_dimension(Hsym_basis)
   allocate(bath(Nb))
   allocate(bath_prev(Nb))
   call ed_init_solver(comm,bath)

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


   function Hloc_model(Mh_,ts_,lambda_) result (H0)
      ! Directly created in lso basis
      integer                                               :: ispin, spsign, lowerbound, upperbound, Nlo
      real(8)                                               :: Mh_,ts_,lambda_
      complex(8),dimension(Nlso,Nlso)                       :: H0
      !
      Nlo = Nlat*Norb
      H0  = zero
      !
      do ispin=1,Nspin
         lowerbound = (ispin-1)*Nlo
         upperbound = ispin*Nlo
         spsign = -2*ispin+3
         H0(lowerbound:upperbound,lowerbound:upperbound) = hloc_matrix(ts_, Mh_, lambda_, spsign)
      enddo
      !
   end function hloc_model


   function hk_model(kpoint,N) result(Hk)
      ! Directly created in lso basis, kpoint MUST be in direct coordinates
      ! the N input is only needed or the TB_build_model routine....
      integer                         :: N, ispin, spsign, lowerbound, upperbound, Nlo
      real(8),dimension(:)            :: kpoint
      complex(8),dimension(N,N)       :: hk, hloc_, h1_, h2_, h3_, h4_, h5_, h6_
      !
      if (N/=Nlso) stop "dimensionality of hk wrong!"
      Nlo = Nlat*Norb
      hk  = zero
      !
      do ispin=1,Nspin
         lowerbound = (ispin-1)*Nlo
         upperbound = ispin*Nlo
         spsign     = -2*ispin+3
         hloc_(lowerbound:upperbound,lowerbound:upperbound) = hloc_matrix(ts, Mh, lambda, spsign)
         h1_(lowerbound:upperbound,lowerbound:upperbound)   = hhop1_matrix(ts, lambda, spsign)
         h2_(lowerbound:upperbound,lowerbound:upperbound)   = hhop2_matrix(ts, lambda, spsign)
         h3_(lowerbound:upperbound,lowerbound:upperbound)   = hhop3_matrix(ts, lambda, spsign)
         h4_(lowerbound:upperbound,lowerbound:upperbound)   = hhop4_matrix(ts, lambda, spsign)
         h5_(lowerbound:upperbound,lowerbound:upperbound)   = hhop5_matrix(ts, lambda, spsign)
         h6_(lowerbound:upperbound,lowerbound:upperbound)   = hhop6_matrix(ts, lambda, spsign)
      enddo
      !
      Hk = hloc_ + h1_*exp(-2*pi*(0,1)*kpoint(2))            + h2_*exp(2*pi*(0,1)*(kpoint(1)-kpoint(2))) &
                 + h3_*exp(2*pi*(0,1)*kpoint(1))             + h4_*exp(2*pi*(0,1)*kpoint(2)) &
                 + h5_*exp(2*pi*(0,1)*(kpoint(2)-kpoint(1))) + h6_*exp(-2*pi*(0,1)*kpoint(1))
      !
   end function hk_model


   !AUXILLIARY MATRIX CONSTRUCTORS for spin Up block, spindown = spinup(-lSOC)

   function hloc_matrix(t,M,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),dimension(6,6)    :: tempmat
      real(8)                   :: t,M,lSOC
      integer(4)                :: spinsign
      !
      if (spinsign==1) then
          continue
      else if (spinsign==-1) then
          lSOC = -1.*lSOC
      else
          stop "Invalid spinsign passed!"
      end if
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = (/ 1, 0, 0, 0, 0, 0 /)
      tempmat(2,:) = (/ 0,-1, 0, 0, 0, 0 /)
      tempmat(3,:) = (/ 0, 0, 1, 0, 0, 0 /)
      tempmat(4,:) = (/ 0, 0, 0,-1, 0, 0 /)
      tempmat(5,:) = (/ 0, 0, 0, 0, 1, 0 /)
      tempmat(6,:) = (/ 0, 0, 0, 0, 0,-1 /)
      hmat         = hmat + M*tempmat
      !
      tempmat(1,:) = (/ 0, 1, 0, 0, 0, 1 /)
      tempmat(2,:) = (/ 1, 0, 1, 0, 0, 0 /)
      tempmat(3,:) = (/ 0, 1, 0, 1, 0, 0 /)
      tempmat(4,:) = (/ 0, 0, 1, 0, 1, 0 /)
      tempmat(5,:) = (/ 0, 0, 0, 1, 0, 1 /)
      tempmat(6,:) = (/ 1, 0, 0, 0, 1, 0 /)
      hmat         = hmat + t*tempmat
      !
      tempmat(1,:) = (/ 0, 0, 1, 0,-1, 0 /)
      tempmat(2,:) = (/ 0, 0, 0, 1, 0,-1 /)
      tempmat(3,:) = (/-1, 0, 0, 0, 1, 0 /)
      tempmat(4,:) = (/ 0,-1, 0, 0, 0, 1 /)
      tempmat(5,:) = (/ 1, 0,-1, 0, 0, 0 /)
      tempmat(6,:) = (/ 0, 1, 0,-1, 0, 0 /)
      hmat         = hmat + lSOC*(0,1)*tempmat
      !
   end function hloc_matrix

   function hhop1_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),dimension(6,6)    :: tempmat
      real(8)                   :: t,lSOC
      integer(4)                :: spinsign
      !
      if (spinsign==1) then
          continue
      else if (spinsign==-1) then
          lSOC = -1.*lSOC
      else
          stop "Invalid spinsign passed!"
      end if
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = (/ 0, 0, 0, 1, 0, 0 /)
      tempmat(2,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(3,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(4,:) = (/ 1, 0, 0, 0, 0, 0 /)
      tempmat(5,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(6,:) = (/ 0, 0, 0, 0, 0, 0 /)
      hmat         = hmat + t*tempmat
      !
      tempmat(1,:) = (/ 0, 0, 1, 0,-1, 0 /)
      tempmat(2,:) = (/ 0, 0, 0, 1, 0, 0 /)
      tempmat(3,:) = (/-1, 0, 0, 0, 0, 0 /)
      tempmat(4,:) = (/ 0,-1, 0, 0, 0, 1 /)
      tempmat(5,:) = (/ 1, 0, 0, 0, 0, 0 /)
      tempmat(6,:) = (/ 0, 0, 0,-1, 0, 0 /)
      hmat         = hmat + lSOC*(0,1)*tempmat
      !
   end function hhop1_matrix

   function hhop2_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),dimension(6,6)    :: tempmat
      real(8)                   :: t,M,lSOC
      integer(4)                :: spinsign
      !
      if (spinsign==1) then
          continue
      else if (spinsign==-1) then
          lSOC = -1.*lSOC
      else
          stop "Invalid spinsign passed!"
      end if
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(2,:) = (/ 0, 0, 0, 0, 1, 0 /)
      tempmat(3,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(4,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(5,:) = (/ 0, 1, 0, 0, 0, 0 /)
      tempmat(6,:) = (/ 0, 0, 0, 0, 0, 0 /)
      hmat         = hmat + t*tempmat
      !
      tempmat(1,:) = (/ 0, 0, 0, 0,-1, 0 /)
      tempmat(2,:) = (/ 0, 0, 0, 1, 0,-1 /)
      tempmat(3,:) = (/ 0, 0, 0, 0, 1, 0 /)
      tempmat(4,:) = (/ 0,-1, 0, 0, 0, 1 /)
      tempmat(5,:) = (/ 1, 0,-1, 0, 0, 0 /)
      tempmat(6,:) = (/ 0, 1, 0, 0, 0, 0 /)
      hmat         = hmat + lSOC*(0,1)*tempmat
      !
   end function hhop2_matrix

   function hhop3_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),dimension(6,6)    :: tempmat
      real(8)                   :: t,lSOC
      integer(4)                :: spinsign
      !
      if (spinsign==1) then
          continue
      else if (spinsign==-1) then
          lSOC = -1.*lSOC
      else
          stop "Invalid spinsign passed!"
      end if
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(2,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(3,:) = (/ 0, 0, 0, 0, 0, 1 /)
      tempmat(4,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(5,:) = (/ 0, 0, 0, 0, 0, 0 /)
      tempmat(6,:) = (/ 0, 0, 1, 0, 0, 0 /)
      hmat         = hmat + t*tempmat
      !
      tempmat(1,:) = (/ 0, 0, 1, 0, 0, 0 /)
      tempmat(2,:) = (/ 0, 0, 0, 0, 0,-1 /)
      tempmat(3,:) = (/-1, 0, 0, 0, 1, 0 /)
      tempmat(4,:) = (/ 0, 0, 0, 0, 0, 1 /)
      tempmat(5,:) = (/ 0, 0,-1, 0, 0, 0 /)
      tempmat(6,:) = (/ 0, 1, 0,-1, 0, 0 /)
      hmat         = hmat + lSOC*(0,1)*tempmat
      !
   end function hhop3_matrix

   function hhop4_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8)                   :: t,lSOC
      integer(4)                :: spinsign
      !
      hmat(:,:) = hhop1_matrix(t,lSOC,spinsign)
   end function hhop4_matrix

   function hhop5_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8)                   :: t,lSOC
      integer(4)                :: spinsign
      !
      hmat(:,:) = hhop2_matrix(t,lSOC,spinsign)
   end function hhop5_matrix

   function hhop6_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8)                   :: t,lSOC
      integer(4)                :: spinsign
      !
      hmat(:,:) = hhop3_matrix(t,lSOC,spinsign)
   end function hhop6_matrix

   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate Hloc and Hk
   !-------------------------------------------------------------------------------------------

   subroutine generate_hk_hloc()
      integer                          :: ik,i,j
      real(8),dimension(Nkx*Nky,2)     :: kgrid
      real(8),dimension(:),allocatable :: gridx,gridy
      !
      !kmesh in direct coordinates
      gridx = linspace(0d0,1d0,Nkx,iend=.false.)
      gridy = linspace(0d0,1d0,Nky,iend=.false.)
      do i=1,Nkx
          do j=1,Nky
              ik=(i-1)*Nkx+j
              kgrid(ik,1)=gridx(i)
              kgrid(ik,2)=gridy(j)
              write(LOGFile,"(2F8.5,A,2F9.5)") gridx(i),gridy(j), " / ", kgrid(ik,:)
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
      Hloc=Hloc_model(Mh,ts,lambda)
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
      stop " naming_convention() NOT IMPLEMENTED!"
      !
      write(LOGfile,"(A)")"The unique index of each site (on the cartesian plane) is as follows:"
      write(LOGfile,"(A)")" "
      do j=1,Ny
         write(LOGfile,"(20(I2,2x))")(matrix(i,Ny+1-j),i =1,Nx)
      enddo
      write(LOGfile,"(A)")" "
   end subroutine naming_convention
   !
end program cdn_kanemele
