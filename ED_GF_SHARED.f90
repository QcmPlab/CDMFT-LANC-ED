MODULE ED_GF_SHARED
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_LINALG,  only: inv,eigh,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  implicit none



  !Lanczos shared variables
  !=========================================================
  complex(8),dimension(:),allocatable   :: state_cvec
  real(8)                               :: state_e


  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:) :: auxGmats,auxGreal



contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2*i
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
  end subroutine deallocate_grids




END MODULE ED_GF_SHARED
