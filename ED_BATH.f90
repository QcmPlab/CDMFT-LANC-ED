MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv,trace
  USE SF_MISC, only: assert_shape
  USE SF_ARRAYS, only: linspace
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_HLOC_DECOMPOSITION
  implicit none

  private


  !interface break_symmetry_bath
     !module procedure :: break_symmetry_bath_site
  !end interface break_symmetry_bath

  !interface spin_symmetrize_bath
     !module procedure ::  spin_symmetrize_bath_site
  !end interface spin_symmetrize_bath

  !interface hermiticize_bath
     !module procedure ::  hermiticize_bath_main
  !end interface hermiticize_bath

  !interface orb_equality_bath
     !module procedure ::  orb_equality_bath_site
  !end interface orb_equality_bath

  !interface ph_symmetrize_bath
     !module procedure ::  ph_symmetrize_bath_site
  !end interface ph_symmetrize_bath

  interface get_bath_dimension
     module procedure ::  get_bath_dimension_direct
     module procedure ::  get_bath_dimension_symmetries
  end interface get_bath_dimension

  interface is_identity
     module procedure ::  is_identity_lso
     module procedure ::  is_identity_nnn
  end interface is_identity


  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  public :: get_bath_dimension
  public :: check_bath_dimension
  !explicit symmetries:
  !public :: hermiticize_bath
  !public :: break_symmetry_bath
  !public :: spin_symmetrize_bath
  !public :: orb_equality_bath
  !public :: ph_symmetrize_bath




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
  public :: bath_from_sym                    !INTERNAL (for effective_bath)
  public :: mask_hloc
  !







  integer              :: ibath,ilat,iorb



contains

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
     mtmp=dreal(nnn2lso_reshape(mnnn,nlat,nspin,norb))
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
     mtmp=dreal(mlso)
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

!+-------------------------------------------------------------------+
!PURPOSE  : Create bath mask
!+-------------------------------------------------------------------+

function mask_hloc(hloc,wdiag,uplo) result(Hmask)
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hloc
    logical,optional                                   :: wdiag,uplo
    logical                                            :: wdiag_,uplo_
    logical,dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hmask
    integer                                            :: ilat,jlat,iorb,jorb,ispin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    Hmask=.false.
    where(abs(Hloc)>1d-6)Hmask=.true.
    !
    !
    if(wdiag_)then
       do ispin=1,Nspin
          do ilat=1,Nlat
             do iorb=1,Norb
                Hmask(ilat,ilat,ispin,ispin,iorb,iorb)=.true.
             enddo
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nspin
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = index_stride_lso(ilat,ispin,iorb)
                      jo = index_stride_lso(jlat,ispin,jorb)
                      if(io>jo)Hmask(ilat,jlat,ispin,ispin,iorb,jorb)=.false.
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    !
  end function mask_hloc



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




END MODULE ED_BATH
