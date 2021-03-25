program cdn_kanemele
   USE CDMFT_ED !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nlso, Nkpath
   real(8)                                                                :: ts,Mh,lambda
   !
   character(len=16)                                                      :: finput
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
   call parse_input_variable(ts,"TS",finput,default=10.d0,comment="hopping parameter")
   call parse_input_variable(Mh,"Mh",finput,default=0.d0,comment="crystal field splitting")
   call parse_input_variable(lambda,"lambda",finput,default=0.3d0,comment="spin-orbit coupling")
   call parse_input_variable(Nkpath,"Nkpath",finput,default=30,comment="Number of k points along kpath")
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

   call generate_bands()

   call finalize_MPI()


contains

   !+------------------------------------------------------------------+
   !PURPOSE  : Hloc for the 2d BHZ model
   !+------------------------------------------------------------------+


   function hloc_model(ts_,Mh_,lambda_) result (hloc_)
      ! Directly created in lso basis
      integer                                               :: ispin, spsign, lowerbound, upperbound, Nlo
      real(8),intent(in)                                    :: Mh_,ts_,lambda_
      complex(8),dimension(Nlso,Nlso)                       :: hloc_
      !
      Nlo = Nlat*Norb
      hloc_  = zero
      !
      do ispin=1,Nspin
         lowerbound = (ispin-1)*Nlo+1
         upperbound = ispin*Nlo
         spsign = -2*ispin+3
         hloc_(lowerbound:upperbound,lowerbound:upperbound) = hloc_matrix(ts_, Mh_, lambda_, spsign)
      enddo
      !
   end function hloc_model


   function hk_model(kpoint,N) result(Hk_)
      ! Directly created in lso basis, kpoint MUST be in direct coordinates
      ! the N input is only needed or the TB_build_model routine....
      integer                         :: N, ispin, spsign, lowerbound, upperbound, Nlo
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
      Hk_ = hloc_ + h1_*exp(-2*pi*(0,1)*kpoint(2))            + h2_*exp(2*pi*(0,1)*(kpoint(1)-kpoint(2))) &
                 + h3_*exp(2*pi*(0,1)*kpoint(1))             + h4_*exp(2*pi*(0,1)*kpoint(2)) &
                 + h5_*exp(2*pi*(0,1)*(kpoint(2)-kpoint(1))) + h6_*exp(-2*pi*(0,1)*kpoint(1))
      !
   end function hk_model


   !AUXILLIARY MATRIX CONSTRUCTORS for spin Up block, spindown = spinup(-lSOC)

   function hloc_matrix(t,M,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),dimension(6,6)    :: tempmat
      real(8),intent(in)        :: t,M,lSOC
      integer(4),intent(in)     :: spinsign
      !
      if (.not.(spinsign==1 .or. spinsign==-1)) stop "Invalid spinsign passed!"
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 1, 0, 0, 0, 0, 0 ]
      tempmat(2,:) = [ 0,-1, 0, 0, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 1, 0, 0, 0 ]
      tempmat(4,:) = [ 0, 0, 0,-1, 0, 0 ]
      tempmat(5,:) = [ 0, 0, 0, 0, 1, 0 ]
      tempmat(6,:) = [ 0, 0, 0, 0, 0,-1 ]
      hmat         = hmat + M*tempmat
      !
      tempmat(1,:) = [ 0, 1, 0, 0, 0, 1 ]
      tempmat(2,:) = [ 1, 0, 1, 0, 0, 0 ]
      tempmat(3,:) = [ 0, 1, 0, 1, 0, 0 ]
      tempmat(4,:) = [ 0, 0, 1, 0, 1, 0 ]
      tempmat(5,:) = [ 0, 0, 0, 1, 0, 1 ]
      tempmat(6,:) = [ 1, 0, 0, 0, 1, 0 ]
      hmat         = hmat + t*tempmat
      !
      tempmat(1,:) = [ 0, 0, 1, 0,-1, 0 ]
      tempmat(2,:) = [ 0, 0, 0, 1, 0,-1 ]
      tempmat(3,:) = [-1, 0, 0, 0, 1, 0 ]
      tempmat(4,:) = [ 0,-1, 0, 0, 0, 1 ]
      tempmat(5,:) = [ 1, 0,-1, 0, 0, 0 ]
      tempmat(6,:) = [ 0, 1, 0,-1, 0, 0 ]
      hmat         = hmat + spinsign*lSOC*(0,1)*tempmat
      !
   end function hloc_matrix

   function hhop1_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),dimension(6,6)    :: tempmat
      real(8),intent(in)        :: t,lSOC
      integer(4),intent(in)     :: spinsign
      !
      if (.not.(spinsign==1 .or. spinsign==-1)) stop "Invalid spinsign passed!"
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0, 1, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(4,:) = [ 1, 0, 0, 0, 0, 0 ]
      tempmat(5,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(6,:) = [ 0, 0, 0, 0, 0, 0 ]
      hmat         = hmat + t*tempmat
      !
      tempmat(1,:) = [ 0, 0, 1, 0,-1, 0 ]
      tempmat(2,:) = [ 0, 0, 0, 1, 0, 0 ]
      tempmat(3,:) = [-1, 0, 0, 0, 0, 0 ]
      tempmat(4,:) = [ 0,-1, 0, 0, 0, 1 ]
      tempmat(5,:) = [ 1, 0, 0, 0, 0, 0 ]
      tempmat(6,:) = [ 0, 0, 0,-1, 0, 0 ]
      hmat         = hmat + spinsign*lSOC*(0,1)*tempmat
      !
   end function hhop1_matrix

   function hhop2_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),dimension(6,6)    :: tempmat
      real(8),intent(in)        :: t,lSOC
      integer(4),intent(in)     :: spinsign
      !
      if (.not.(spinsign==1 .or. spinsign==-1)) stop "Invalid spinsign passed!"
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 0, 0, 1, 0 ]
      tempmat(3,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(4,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(5,:) = [ 0, 1, 0, 0, 0, 0 ]
      tempmat(6,:) = [ 0, 0, 0, 0, 0, 0 ]
      hmat         = hmat + t*tempmat
      !
      tempmat(1,:) = [ 0, 0, 0, 0,-1, 0 ]
      tempmat(2,:) = [ 0, 0, 0, 1, 0,-1 ]
      tempmat(3,:) = [ 0, 0, 0, 0, 1, 0 ]
      tempmat(4,:) = [ 0,-1, 0, 0, 0, 1 ]
      tempmat(5,:) = [ 1, 0,-1, 0, 0, 0 ]
      tempmat(6,:) = [ 0, 1, 0, 0, 0, 0 ]
      hmat         = hmat + spinsign*lSOC*(0,1)*tempmat
      !
   end function hhop2_matrix

   function hhop3_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),dimension(6,6)    :: tempmat
      real(8),intent(in)        :: t,lSOC
      integer(4),intent(in)     :: spinsign
      !
      if (.not.(spinsign==1 .or. spinsign==-1)) stop "Invalid spinsign passed!"
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 0, 0, 0, 1 ]
      tempmat(4,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(5,:) = [ 0, 0, 0, 0, 0, 0 ]
      tempmat(6,:) = [ 0, 0, 1, 0, 0, 0 ]
      hmat         = hmat + t*tempmat
      !
      tempmat(1,:) = [ 0, 0, 1, 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 0, 0, 0,-1 ]
      tempmat(3,:) = [-1, 0, 0, 0, 1, 0 ]
      tempmat(4,:) = [ 0, 0, 0, 0, 0, 1 ]
      tempmat(5,:) = [ 0, 0,-1, 0, 0, 0 ]
      tempmat(6,:) = [ 0, 1, 0,-1, 0, 0 ]
      hmat         = hmat + spinsign*lSOC*(0,1)*tempmat
      !
   end function hhop3_matrix

   function hhop4_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),intent(in)        :: t,lSOC
      integer(4),intent(in)     :: spinsign
      !
      hmat(:,:) = hhop1_matrix(t,lSOC,spinsign)
   end function hhop4_matrix

   function hhop5_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),intent(in)        :: t,lSOC
      integer(4),intent(in)     :: spinsign
      !
      hmat(:,:) = hhop2_matrix(t,lSOC,spinsign)
   end function hhop5_matrix

   function hhop6_matrix(t,lSOC,spinsign) result(hmat)
      complex(8),dimension(6,6) :: hmat
      real(8),intent(in)        :: t,lSOC
      integer(4),intent(in)     :: spinsign
      !
      hmat(:,:) = hhop3_matrix(t,lSOC,spinsign)
   end function hhop6_matrix

   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate Hloc and Hk
   !-------------------------------------------------------------------------------------------

   subroutine generate_bands()
      integer                      :: ik,i,j
      real(8),dimension(4,2)       :: kpath
      real(8),dimension(2)         :: pointK, pointKp, bk1, bk2
      !
      ! These are only to make the TB_Solve_routine run
      bk1 = 2*pi*[1,0]
      bk2 = 2*pi*[0,1]
      !
      pointK  = [1d0/3d0, -1d0/3d0]
      pointKp = [2d0/3d0, 1d0/3d0]
      write(LOGFile,*) pointK
      write(LOGFile,*) pointKp
      KPath(1,:)=[0,0]
      KPath(2,:)=pointK
      Kpath(3,:)=pointKp
      KPath(4,:)=[0d0,0d0]
      write(LOGFile,*) KPath
      call TB_set_bk(bkx=bk1,bky=bk2)
      !call TB_set_ei(eix=[1d0,0d0],eiy=[0d0,1d0])
      call TB_Solve_model(hk_model,Nlso,KPath,Nkpath,&
           colors_name=[red1,blue1,red1,blue1,red1,blue1,red1,blue1,red1,blue1,red1,blue1],&
           points_name=[character(len=10) :: "G","K","K`","G"],&
           file="Eigenbands.nint",iproject=.false.)
   end subroutine generate_bands




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
end program cdn_kanemele
