MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv,trace
  USE SF_MISC, only: assert_shape
  USE SF_ARRAYS, only: linspace
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  private

  

  interface is_identity
     module procedure ::  is_identity_lso
     module procedure ::  is_identity_nnn
  end interface is_identity

  interface is_diagonal
     module procedure ::  is_diagonal_lso
     module procedure ::  is_diagonal_nnn
  end interface is_diagonal

  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  
  interface impose_equal_lambda
     module procedure ::  impose_equal_lambda
  end interface impose_equal_lambda

  interface get_bath_dimension
     module procedure ::  get_bath_dimension_direct
     module procedure ::  get_bath_dimension_symmetries
  end interface get_bath_dimension
  
  interface set_Hreplica
     module procedure init_Hreplica_direct_lso
     module procedure init_Hreplica_direct_nnn
     module procedure init_Hreplica_symmetries_site
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure init_Hreplica_symmetries_lattice
#endif
  end interface set_Hreplica
  
  
  public :: get_bath_dimension
  public :: check_bath_dimension
  !explicit symmetries:
  public :: impose_equal_lambda
  public :: impose_bath_offset
  !
  public :: set_Hreplica




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
  public :: write_dmft_bath                  !INTERNAL (for effective_bath)
  public :: save_dmft_bath                   !INTERNAL (for effective_bath)
  public :: set_dmft_bath                    !INTERNAL (for effective_bath)
  public :: get_dmft_bath                    !INTERNAL (for effective_bath)
  !

  public :: hreplica_build                   !INTERNAL (for effective_bath)
  public :: hreplica_mask                    !INTERNAL (for effective_bath)
#if __GFORTRAN__ &&  __GNUC__ > 8     
  public :: hreplica_site                    !INTERNAL (for effective_bath)
#endif




  integer              :: ibath,ilat,iorb



contains

!+-------------------------------------------------------------------+
!PURPOSE  : Check if a matrix is diagonal
!+-------------------------------------------------------------------+

   function is_diagonal_nnn(mnnn) result(flag)
     complex(8),dimension(nlat,nlat,nspin,nspin,norb,norb)                    :: mnnn
     real(8),dimension(nlat*nspin*norb,nlat*nspin*norb)                       :: mtmp
     integer                                                                  :: i,j
     logical                                                                  :: flag
     !
     flag=.true.
     !
     if ( ANY( abs(dimag(mnnn)) .gt. 1d-6 ) ) flag=.false.
     !
     mtmp=abs(dreal(nnn2lso_reshape(mnnn,nlat,nspin,norb)))
     !
     do i=1,nlat*nspin*norb
        do j=1,nlat*nspin*norb
           if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
         enddo
     enddo
     !
   end function is_diagonal_nnn


   function is_diagonal_lso(mlso) result(flag)
     complex(8),dimension(nlat*nspin*norb,nlat*nspin*norb)                    :: mlso
     real(8),dimension(nlat*nspin*norb,nlat*nspin*norb)                       :: mtmp
     integer                                                                  :: i,j
     logical                                                                  :: flag
     !
     flag=.true.
     !
     if ( ANY( abs(dimag(mlso)) .gt. 1d-6 ) ) flag=.false.
     !
     mtmp=abs(dreal(mlso))
     !
     do i=1,nlat*nspin*norb
        do j=1,nlat*nspin*norb
           if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
         enddo
     enddo
     !
   end function is_diagonal_lso
!+-------------------------------------------------------------------+
!PURPOSE  : Check if a matrix is the identity
!+-------------------------------------------------------------------+

   function is_identity_nnn(mnnn) result(flag)
     complex(8),dimension(nlat,nlat,nspin,nspin,norb,norb)                    :: mnnn
     real(8),dimension(nlat*nspin*norb,nlat*nspin*norb)                       :: mtmp
     integer                                                                  :: i,j
     logical                                                                  :: flag
     !
     flag=.true.
     !
     if ( ANY( abs(dimag(mnnn)) .gt. 1d-6 ) ) flag=.false.
     !
     mtmp=abs(dreal(nnn2lso_reshape(mnnn,nlat,nspin,norb)))
     !
     do i=1,nlat*nspin*norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
     enddo
     !
     do i=1,nlat*nspin*norb
        do j=1,nlat*nspin*norb
           if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
         enddo
     enddo
     !
   end function is_identity_nnn


   function is_identity_lso(mlso) result(flag)
     complex(8),dimension(nlat*nspin*norb,nlat*nspin*norb)                    :: mlso
     real(8),dimension(nlat*nspin*norb,nlat*nspin*norb)                       :: mtmp
     integer                                                                  :: i,j
     logical                                                                  :: flag
     !
     flag=.true.
     !
     if ( ANY( abs(dimag(mlso)) .gt. 1d-6 ) ) flag=.false.
     !
     mtmp=abs(dreal(mlso))
     !
     do i=1,nlat*nspin*norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
     enddo
     !
     do i=1,nlat*nspin*norb
        do j=1,nlat*nspin*norb
           if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
         enddo
     enddo
     !
   end function is_identity_lso



  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  include 'ED_BATH/user_aux.f90'


  !##################################################################
  !
  !     H_REPLICA ROUTINES:
  !
  !##################################################################
  include 'ED_BATH/hreplica_setup.f90'
#if __GFORTRAN__ &&  __GNUC__ > 8     
  include 'ED_BATH/hreplica_setup_lattice.f90'
#endif


  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  include 'ED_BATH/dmft_aux.f90'




END MODULE ED_BATH
