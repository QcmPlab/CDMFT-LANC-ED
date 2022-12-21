MODULE ED_FIT_CHI2
   USE SF_CONSTANTS
   USE SF_OPTIMIZE, only:fmin_cg,fmin_cgminimize
   USE SF_LINALG,   only:eye,zeye,inv,inv_her,trace,operator(.x.) !BLAS xgemm operator overloading
   USE SF_IOTOOLS,  only:reg,free_unit,str
   USE SF_ARRAYS,   only:arange
   USE SF_MISC,     only:assert_shape
   USE ED_INPUT_VARS
   USE ED_VARS_GLOBAL
   USE ED_AUX_FUNX
   USE ED_BATH
   USE ED_BATH_FUNCTIONS
   USE ED_FIT_REPLICA
   USE ED_FIT_GENERAL


   implicit none
   private

   interface ed_chi2_fitgf
      module procedure chi2_fitgf_generic_normal
#if __GFORTRAN__ &&  __GNUC__ > 8
      !RDMFT_WRAPPER
      module procedure chi2_fitgf_lattice_normal
#endif
   end interface ed_chi2_fitgf


   public :: ed_chi2_fitgf


   integer                                               :: Ldelta
   complex(8),dimension(:,:,:,:,:,:,:),allocatable       :: FGmatrix
   logical(8),dimension(:,:,:,:,:,:),allocatable         :: Hmask
   complex(8),dimension(:,:),allocatable                 :: Fdelta
   real(8),dimension(:),allocatable                      :: Xdelta,Wdelta
   integer                                               :: totNorb,totNspin
   integer,dimension(:),allocatable                      :: getIorb,getJorb,getIspin,getJspin,getIlat,getJlat
   integer                                               :: Orb_indx,Spin_indx,Spin_mask
   !location of the maximum of the chisquare over Nlso.
   integer                                               :: maxchi_loc
   !
   type nsymm_vector
      real(8),dimension(:),allocatable                   :: element
   end type nsymm_vector
   !

contains

   !+----------------------------------------------------------------------+
   !PURPOSE  : Chi^2 fit of the G0/Delta
   !+----------------------------------------------------------------------+
   subroutine chi2_fitgf_generic_normal(fg,bath)
      complex(8),dimension(:,:,:,:,:,:,:) :: fg ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Niw]
      real(8),dimension(:)                :: bath
      !
      call assert_shape(fg,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(fg,7)],"chi2_fitgf_generic_normal","fg")
      !
      select case(cg_method)
       case default
         stop "ED Error: cg_method > 1"
       case (0)
         if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
       case (1)
         if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
      end select
      !
      select case(bath_type)
       case('replica')
         call chi2_fitgf_replica(fg,bath)
       case('general')
         call chi2_fitgf_general(fg,bath)
      end select
      !
      !set trim_state_list to true after the first fit has been done: this
      !marks the ends of the cycle of the 1st DMFT loop.
      trim_state_list=.true.
      !
   end subroutine chi2_fitgf_generic_normal


#if __GFORTRAN__ &&  __GNUC__ > 8
   !+----------------------------------------------------------------------!
   ! PURPOSE: given a number of independent baths, evaluate N independent
   ! Delta/G0 functions and fit them to update the effective baths for ED.
   !+----------------------------------------------------------------------!
   !RDMFT WRAPPER:
   subroutine chi2_fitgf_lattice_normal(fg,bath)
      real(8),dimension(:,:)                    :: bath
      complex(8),dimension(:,:,:,:,:,:,:,:)     :: fg
      !MPI auxiliary vars
      integer                                   :: isites
      integer                                   :: Nsites
      character(len=5)                          :: tmp_suffix
      !
      ! Check dimensions !
      Nsites=size(bath,1)
      call assert_shape(fg,[Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(fg,8)],"chi2_fitgf_generic_normal","fg")
      !
      !
      do isites = 1, Nsites
         !
         ed_file_suffix=reg(ineq_site_suffix)//str(isites,site_indx_padding)
         !
         call chi2_fitgf_generic_normal(fg(isites,:,:,:,:,:,:,:),bath(isites,:))
         !
      end do
      !
      !
      ed_file_suffix=""
   end subroutine chi2_fitgf_lattice_normal

#endif

end MODULE ED_FIT_CHI2
