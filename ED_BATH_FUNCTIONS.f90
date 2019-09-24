MODULE ED_BATH_FUNCTIONS
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv,zeye,inv_her
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none

  private


  !##################################################################
  !
  !\DELTA HYBRIDIZATION FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface delta_bath_mats
     module procedure delta_bath_mats_main
     module procedure delta_bath_mats_ispin_jspin
     module procedure delta_bath_mats_ispin_jspin_iorb_jorb
     module procedure delta_bath_mats_main_
     module procedure delta_bath_mats_ispin_jspin_
     module procedure delta_bath_mats_ispin_jspin_iorb_jorb_
  end interface delta_bath_mats
  !
  interface ed_delta_matsubara
     module procedure delta_bath_mats_main_
     module procedure delta_bath_mats_ispin_jspin_
     module procedure delta_bath_mats_ispin_jspin_iorb_jorb_
  end interface ed_delta_matsubara



  !##################################################################
  !
  !\DELTA HYBRIDIZATION FUNCTION REAL
  !
  !##################################################################
  !NORMAL
  interface delta_bath_real
     module procedure delta_bath_real_main
     module procedure delta_bath_real_ispin_jspin
     module procedure delta_bath_real_ispin_jspin_iorb_jorb
     module procedure delta_bath_real_main_
     module procedure delta_bath_real_ispin_jspin_
     module procedure delta_bath_real_ispin_jspin_iorb_jorb_
  end interface delta_bath_real
  !
  interface ed_delta_realaxis
     module procedure delta_bath_real_main_
     module procedure delta_bath_real_ispin_jspin_
     module procedure delta_bath_real_ispin_jspin_iorb_jorb_
  end interface ed_delta_realaxis


  !##################################################################
  !
  !NON-INTERACTING GREEN'S FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface g0and_bath_mats
     module procedure g0and_bath_mats_main
     module procedure g0and_bath_mats_ispin_jspin
     module procedure g0and_bath_mats_ispin_jspin_iorb_jorb
     module procedure g0and_bath_mats_main_
     module procedure g0and_bath_mats_ispin_jspin_
     module procedure g0and_bath_mats_ispin_jspin_iorb_jorb_
  end interface g0and_bath_mats
  !
  interface ed_g0and_matsubara
     module procedure g0and_bath_mats_main_
     module procedure g0and_bath_mats_ispin_jspin_
     module procedure g0and_bath_mats_ispin_jspin_iorb_jorb_
  end interface ed_g0and_matsubara
  !



  !##################################################################
  !
  !INVERSE NON-INTERACTING GREEN'S FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface invg0_bath_mats
     module procedure invg0_bath_mats_main
     module procedure invg0_bath_mats_ispin_jspin
     module procedure invg0_bath_mats_ispin_jspin_iorb_jorb
     module procedure invg0_bath_mats_main_
     module procedure invg0_bath_mats_ispin_jspin_
     module procedure invg0_bath_mats_ispin_jspin_iorb_jorb_
  end interface invg0_bath_mats
  !
  interface ed_invg0_matsubara
     module procedure invg0_bath_mats_main_
     module procedure invg0_bath_mats_ispin_jspin_
     module procedure invg0_bath_mats_ispin_jspin_iorb_jorb_
  end interface ed_invg0_matsubara




  !##################################################################
  !
  !NON-INTERACTING GREEN'S FUNCTION REAL-AXIS
  !
  !##################################################################
  !NORMAL
  interface g0and_bath_real
     module procedure g0and_bath_real_main
     module procedure g0and_bath_real_ispin_jspin
     module procedure g0and_bath_real_ispin_jspin_iorb_jorb
     module procedure g0and_bath_real_main_
     module procedure g0and_bath_real_ispin_jspin_
     module procedure g0and_bath_real_ispin_jspin_iorb_jorb_
  end interface g0and_bath_real
  !
  interface ed_g0and_realaxis
     module procedure g0and_bath_real_main_
     module procedure g0and_bath_real_ispin_jspin_
     module procedure g0and_bath_real_ispin_jspin_iorb_jorb_
  end interface ed_g0and_realaxis
  !



  !##################################################################
  !
  !INVERSE NON-INTERACTING GREEN'S FUNCTION REAL-AXIS
  !
  !##################################################################
  !NORMAL
  interface invg0_bath_real
     module procedure invg0_bath_real_main
     module procedure invg0_bath_real_ispin_jspin
     module procedure invg0_bath_real_ispin_jspin_iorb_jorb
     module procedure invg0_bath_real_main_
     module procedure invg0_bath_real_ispin_jspin_
     module procedure invg0_bath_real_ispin_jspin_iorb_jorb_
  end interface invg0_bath_real
  !
  interface ed_invg0_realaxis
     module procedure invg0_bath_real_main_
     module procedure invg0_bath_real_ispin_jspin_
     module procedure invg0_bath_real_ispin_jspin_iorb_jorb_
  end interface ed_invg0_realaxis
  !


  !INTERNAL USE:
  public :: delta_bath_mats
  public :: g0and_bath_mats
  public :: invg0_bath_mats
  !
  public :: g0and_bath_real
  public :: delta_bath_real
  public :: invg0_bath_real


  !EXTERNAL USE, per USER:
  public :: ed_delta_matsubara
  public :: ed_g0and_matsubara
  public :: ed_invg0_matsubara
  !
  public :: ed_delta_realaxis
  public :: ed_g0and_realaxis
  public :: ed_invg0_realaxis



contains




  !##################################################################
  !
  !     DELTA FUNCTIONS
  !     G0 FUNCTIONS
  !     G0^{-1} FUNCTIONS
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  Delta and Fdelta functions on the Matsubara axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function delta_bath_mats_main(x,dmft_bath_) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                             :: i,ih,L
    integer                                             :: iorb,jorb,ispin,jspin,ibath
    integer                                             :: io,jo
    real(8),dimension(Nbath)                            :: eps,dps,vps
    real(8),dimension(Norb,Nbath)                       :: vops
    !
    real(8),dimension(Nspin,Nbath)                      :: ehel
    real(8),dimension(Nspin,Nspin,Nbath)                :: whel
    real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wohel
    !
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(x)) :: invH_k
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)   :: invH_knn
    !
    !
    Delta=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
       do ispin=1,Nspin
          do iorb=1,Norb
             eps = dmft_bath_%e(ispin,iorb,1:Nbath)
             vps = dmft_bath_%v(ispin,iorb,1:Nbath)
             do i=1,L
                Delta(ispin,ispin,iorb,iorb,i) = sum( vps(:)*vps(:)/(x(i) - eps(:)) )
             enddo
          enddo
       enddo
       !

       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !

       !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(iw_n - E(k)) ]
       do ispin=1,Nspin
          eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
          vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
          do iorb=1,Norb
             do jorb=1,Norb
                do i=1,L
                   Delta(ispin,ispin,iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(x(i) - eps(:)) )
                enddo
             enddo
          enddo
       enddo

       !
    case ("replica")
       !

       invH_k=zero
       do i=1,L
          invH_knn=zero
          do ibath=1,Nbath
             !
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1) * Norb
                         jo = jorb + (jspin-1) * Norb
                         invH_k(io,jo,i)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)
                      enddo
                   enddo
                enddo
             enddo
             !
             invH_k(:,:,i) = zeye(Nspin*Norb) * x(i) - invH_k(:,:,i)
             call inv(invH_k(:,:,i))
             invH_knn(:,:,:,:,ibath)=so2nn_reshape(invH_k(:,:,i),Nspin,Norb)
             !
          enddo
          !
          do ibath=1,Nbath
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         Delta(ispin,jspin,iorb,jorb,i)=Delta(ispin,jspin,iorb,jorb,i)+ &
                              (dmft_bath_%vr(ibath)) * invH_knn(ispin,jspin,iorb,jorb,ibath) * dmft_bath_%vr(ibath)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       !
    end select
  end function delta_bath_mats_main


  function delta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    Delta = delta_bath_mats_main(x,dmft_bath_)
    G0out = Delta(ispin,jspin,:,:,:)
  end function delta_bath_mats_ispin_jspin


  function delta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(size(x))                       :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    Delta = delta_bath_mats_main(x,dmft_bath_)
    G0out = Delta(ispin,jspin,iorb,jorb,:)
  end function delta_bath_mats_ispin_jspin_iorb_jorb


  function delta_bath_mats_main_(x,bath_) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    Delta = delta_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_mats_main_


  function delta_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(:),intent(in)      :: x
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = delta_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_mats_ispin_jspin_

  function delta_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    type(effective_bath)               :: dmft_bath_
    complex(8),dimension(:),intent(in) :: x
    complex(8),dimension(size(x))      :: G0out
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_mats_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = delta_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_mats_ispin_jspin_iorb_jorb_







  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  Delta and Fdelta functions on the Real axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  function delta_bath_real_main(x,dmft_bath_) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                             :: i,ih,L
    integer                                             :: iorb,jorb,ispin,jspin,ibath
    integer                                             :: io,jo
    real(8),dimension(Nbath)                            :: eps,dps,vps
    real(8),dimension(Norb,Nbath)                       :: vops
    !
    real(8),dimension(Nspin,Nbath)                      :: ehel
    real(8),dimension(Nspin,Nspin,Nbath)                :: whel
    real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wohel
    !
    complex(8),dimension(Nspin*Norb,Nspin*Norb,size(x)) :: invH_k
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)   :: invH_knn
    !
    !
    Delta=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default
       !
       !
       !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(w+i\h - E_{a}(k)) ]
       do ispin=1,Nspin
          do iorb=1,Norb
             eps = dmft_bath_%e(ispin,iorb,1:Nbath)
             vps = dmft_bath_%v(ispin,iorb,1:Nbath)
             do i=1,L
                Delta(ispin,ispin,iorb,iorb,i) = sum( vps(:)*vps(:)/(x(i) - eps(:)) )
             enddo
          enddo
       enddo
       !

       !
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       !

       !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(w+i\h  - E(k)) ]
       do ispin=1,Nspin
          eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
          vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
          do iorb=1,Norb
             do jorb=1,Norb
                do i=1,L
                   Delta(ispin,ispin,iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(x(i) - eps(:)) )
                enddo
             enddo
          enddo
       enddo

    case ("replica")
       !

       invH_k=zero
       do i=1,L
          invH_knn=zero
          do ibath=1,Nbath
             !
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1) * Norb
                         jo = jorb + (jspin-1) * Norb
                         invH_k(io,jo,i)=dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)
                      enddo
                   enddo
                enddo
             enddo
             !
             invH_k(:,:,i) = zeye(Nspin*Norb) * x(i) - invH_k(:,:,i)
             call inv(invH_k(:,:,i))
             invH_knn(:,:,:,:,ibath)=so2nn_reshape(invH_k(:,:,i),Nspin,Norb)
             !
          enddo
          !
          do ibath=1,Nbath
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         Delta(ispin,jspin,iorb,jorb,i)=Delta(ispin,jspin,iorb,jorb,i)+ &
                              (dmft_bath_%vr(ibath)) * invH_knn(ispin,jspin,iorb,jorb,ibath) * dmft_bath_%vr(ibath)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end function delta_bath_real_main


  function delta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    Delta = delta_bath_real_main(x,dmft_bath_)
    G0out = Delta(ispin,jspin,:,:,:)
  end function delta_bath_real_ispin_jspin


  function delta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(size(x))                       :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    Delta = delta_bath_real_main(x,dmft_bath_)
    G0out = Delta(ispin,jspin,iorb,jorb,:)
  end function delta_bath_real_ispin_jspin_iorb_jorb


  function delta_bath_real_main_(x,bath_) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    Delta = delta_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_real_main_


  function delta_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = delta_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_real_ispin_jspin_

  function delta_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8),dimension(size(x))      :: G0out
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = delta_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function delta_bath_real_ispin_jspin_iorb_jorb_



















  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the G0 function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0 and F0 non-interacting Green's functions on the Matsubara axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function g0and_bath_mats_main(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta,Fdelta
    integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    real(8),dimension(size(x))                          :: det
    complex(8),dimension(size(x))                       :: fg,ff
    complex(8),dimension(:,:),allocatable               :: fgorb,zeta
    !
    G0and = zero
    !
    L=size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       !
       Delta = delta_bath_mats(x,dmft_bath_)
       do ispin=1,Nspin
          do iorb=1,Norb
             fg(:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
             G0and(ispin,ispin,iorb,iorb,:) = one/fg(:)
          enddo
       enddo
       !
       !
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       !
       allocate(fgorb(Norb,Norb),zeta(Norb,Norb))
       Delta = delta_bath_mats(x,dmft_bath_)
       do ispin=1,Nspin         !Spin diagonal
          do i=1,L
             fgorb= zero
             zeta = (x(i)+xmu)*eye(Norb)
             do iorb=1,Norb
                do jorb=1,Norb
                   fgorb(iorb,jorb) = zeta(iorb,jorb)-impHloc(ispin,ispin,iorb,jorb)-Delta(ispin,ispin,iorb,jorb,i)
                enddo
             enddo
             call inv(fgorb)
             G0and(ispin,ispin,:,:,i)=fgorb
          enddo
       enddo
       deallocate(fgorb,zeta)
       !
    end select
  end function g0and_bath_mats_main


  function g0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = g0and_bath_mats_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,:,:,:)
  end function g0and_bath_mats_ispin_jspin


  function g0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: G0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = g0and_bath_mats_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,iorb,jorb,:)
  end function g0and_bath_mats_ispin_jspin_iorb_jorb


  function g0and_bath_mats_main_(x,bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0and = g0and_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_mats_main_


  function g0and_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_mats_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = g0and_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_mats_ispin_jspin_

  function g0and_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: G0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_mats_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = g0and_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_mats_ispin_jspin_iorb_jorb_










  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0 and F0 non-interacting Green's functions on the real-axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function g0and_bath_real_main(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta,Fdelta
    integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    complex(8),dimension(size(x))                       :: det,fg,ff
    complex(8),dimension(:,:),allocatable               :: fgorb,zeta
    !
    G0and = zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       Delta = delta_bath_real(x,dmft_bath_)
       do ispin=1,Nspin
          do iorb=1,Norb
             fg(:)    = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
             G0and(ispin,ispin,iorb,iorb,:) = one/fg(:)
          enddo
       enddo
       !
       !
       !
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       !
       allocate(fgorb(Norb,Norb),zeta(Norb,Norb))
       Delta = delta_bath_real(x,dmft_bath_)
       do ispin=1,Nspin
          do i=1,L
             fgorb= zero
             zeta = (x(i)+xmu)*eye(Norb)
             do iorb=1,Norb
                do jorb=1,Norb
                   fgorb(iorb,jorb) = zeta(iorb,jorb)-impHloc(ispin,ispin,iorb,jorb)-Delta(ispin,ispin,iorb,jorb,i)
                enddo
             enddo
             call inv(fgorb)
             G0and(ispin,ispin,:,:,i)=fgorb
          enddo
       enddo
       deallocate(fgorb,zeta)
    end select
  end function g0and_bath_real_main


  function g0and_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = g0and_bath_real_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,:,:,:)
  end function g0and_bath_real_ispin_jspin


  function g0and_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: G0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = g0and_bath_real_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,iorb,jorb,:)
  end function g0and_bath_real_ispin_jspin_iorb_jorb


  function g0and_bath_real_main_(x,bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0and = g0and_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_real_main_


  function g0and_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = g0and_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_real_ispin_jspin_


  function g0and_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: G0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = g0and_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function g0and_bath_real_ispin_jspin_iorb_jorb_









  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the inverse G0 function at a point x from
  ! type(effective_bath) :: dmft_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0^{-1} and F0{-1} non-interacting Green's functions on the Matsubara axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function invg0_bath_mats_main(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta,Fdelta
    integer                                             :: i,iorb,jorb,ispin,jspin,io,jo,Nso,L
    !complex(8),dimension(size(x))                       :: fg,ff
    complex(8),dimension(:,:),allocatable               :: zeta!,fgorb
    !
    G0and = zero
    !
    L=size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       Delta = delta_bath_mats(x,dmft_bath_)
       do ispin=1,Nspin
          do iorb=1,Norb
             G0and(ispin,ispin,iorb,iorb,:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !
       !
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       allocate(zeta(Norb,Norb))
       Delta = delta_bath_mats(x,dmft_bath_)
       do ispin=1,Nspin
          do i=1,L
             zeta = (x(i)+xmu)*eye(Norb)
             do iorb=1,Norb
                do jorb=1,Norb
                   G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb)-impHloc(ispin,ispin,iorb,jorb)-Delta(ispin,ispin,iorb,jorb,i)
                enddo
             enddo
          enddo
       enddo
       deallocate(zeta)
       !
    end select
    !
  end function invg0_bath_mats_main


  function invg0_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = invg0_bath_mats_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,:,:,:)
  end function invg0_bath_mats_ispin_jspin


  function invg0_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: G0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = invg0_bath_mats_main(x,dmft_bath_)
    G0out = G0and(ispin,jspin,iorb,jorb,:)
  end function invg0_bath_mats_ispin_jspin_iorb_jorb


  function invg0_bath_mats_main_(x,bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0and = invg0_bath_mats_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_mats_main_

  function invg0_bath_mats_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_mats_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = invg0_bath_mats_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_mats_ispin_jspin_

  function invg0_bath_mats_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: G0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_mats_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = invg0_bath_mats_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_mats_ispin_jspin_iorb_jorb_










  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0 and F0 non-interacting Green's functions on the real-axis:
  ! _1 : input type(effective_bath) dmft_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function invg0_bath_real_main(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    !complex(8),dimension(size(x))                       :: fg,ff
    complex(8),dimension(:,:),allocatable               :: zeta!,fgorb
    !
    G0and = zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       Delta = delta_bath_real(x,dmft_bath_)
       do ispin=1,Nspin
          do iorb=1,Norb
             G0and(ispin,ispin,iorb,iorb,:) =  x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !

       allocate(zeta(Norb,Norb))
       Delta = delta_bath_real(x,dmft_bath_)
       do ispin=1,Nspin
          do i=1,L
             zeta = (x(i)+xmu)*eye(Norb)
             do iorb=1,Norb
                do jorb=1,Norb
                   G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb) - impHloc(ispin,ispin,iorb,jorb) - Delta(ispin,ispin,iorb,jorb,i)
                enddo
             enddo
          enddo
       enddo
       deallocate(zeta)
       !
    end select
    !
  end function invg0_bath_real_main


  function invg0_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x))             :: G0out
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = invg0_bath_real_main(x,dmft_bath_)
    G0out = zero
    G0out = G0and(ispin,jspin,:,:,:)
  end function invg0_bath_real_ispin_jspin

  function invg0_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_) result(G0out)
    integer,intent(in)                                  :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8)                                          :: G0out(size(x))
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    G0and = invg0_bath_real_main(x,dmft_bath_)
    G0out = zero
    G0out = G0and(ispin,jspin,iorb,jorb,:)
  end function invg0_bath_real_ispin_jspin_iorb_jorb


  function invg0_bath_real_main_(x,bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_real_main_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0and = invg0_bath_real_main(x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_real_main_


  function invg0_bath_real_ispin_jspin_(ispin,jspin,x,bath_) result(G0out)
    integer,intent(in)                      :: ispin,jspin
    complex(8),dimension(:),intent(in)      :: x
    type(effective_bath)                    :: dmft_bath_
    complex(8),dimension(Norb,Norb,size(x)) :: G0out
    real(8),dimension(:)                    :: bath_
    logical                                 :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = zero
    G0out = invg0_bath_real_ispin_jspin(ispin,jspin,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_real_ispin_jspin_

  function invg0_bath_real_ispin_jspin_iorb_jorb_(ispin,jspin,iorb,jorb,x,bath_) result(G0out)
    integer,intent(in)                 :: iorb,jorb,ispin,jspin
    complex(8),dimension(:),intent(in) :: x
    type(effective_bath)               :: dmft_bath_
    complex(8)                         :: G0out(size(x))
    real(8),dimension(:)               :: bath_
    logical                            :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "invg0_bath_real_ error: wrong bath dimensions"
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    G0out = zero
    G0out = invg0_bath_real_ispin_jspin_iorb_jorb(ispin,jspin,iorb,jorb,x,dmft_bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end function invg0_bath_real_ispin_jspin_iorb_jorb_











END MODULE ED_BATH_FUNCTIONS
