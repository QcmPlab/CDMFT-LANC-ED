MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  private


  interface break_symmetry_bath
     module procedure :: break_symmetry_bath_site
  end interface break_symmetry_bath

  interface spin_symmetrize_bath
     module procedure ::  spin_symmetrize_bath_site
  end interface spin_symmetrize_bath


  interface orb_equality_bath
     module procedure ::  orb_equality_bath_site
  end interface orb_equality_bath

  interface ph_symmetrize_bath
     module procedure ::  ph_symmetrize_bath_site
  end interface ph_symmetrize_bath





  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  public :: get_bath_dimension
  public :: check_bath_dimension
  !explicit symmetries:
  public :: break_symmetry_bath
  public :: spin_symmetrize_bath
  public :: orb_equality_bath
  public :: ph_symmetrize_bath




  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  !PUBLIC   = transparent to the final user
  !INTERNAL = opaque to the user but available for internal use in the code.
  !
  !DMFT BATH procedures:
  public :: allocate_dmft_bath               !INTERNAL (for effective_bath)
  public :: deallocate_dmft_bath             !INTERNAL (for effective_bath)
  public :: init_dmft_bath                   !INTERNAL (for effective_bath)
  public :: init_dmft_bath_mask              !INTERNAL (for effective_bath)
  public :: write_dmft_bath                  !INTERNAL (for effective_bath)
  public :: save_dmft_bath                   !INTERNAL (for effective_bath)
  public :: set_dmft_bath                    !INTERNAL (for effective_bath)
  public :: get_dmft_bath                    !INTERNAL (for effective_bath)
  !







  integer              :: ibath,ilat,iorb



contains


  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  include 'ED_BATH/user_aux.f90'



  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  include 'ED_BATH/dmft_aux.f90'



  !##################################################################
  !
  !     USER BATH CHECKS:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_,Hloc_nn) result(bool)
    real(8),dimension(:)        :: bath_
    integer                     :: Ntrue
    logical                     :: bool
    real(8),optional,intent(in) :: Hloc_nn(:,:,:,:,:,:)![Nlat][:][Nspin][:][Norb][:]
    if (present(Hloc_nn))then
       Ntrue = get_bath_dimension(one*Hloc_nn)
    else
       Ntrue = get_bath_dimension()
    endif
    bool  = ( size(bath_) == Ntrue )
  end function check_bath_dimension


END MODULE ED_BATH
