MODULE ED_FIT_CHI2
  USE SF_CONSTANTS
  USE SF_OPTIMIZE, only:fmin_cg,fmin_cgminimize
  USE SF_LINALG,   only:eye,zeye,inv,inv_her,operator(.x.)
  USE SF_IOTOOLS,  only:reg,free_unit,txtfy
  USE SF_ARRAYS,   only:arange
  USE SF_MISC,     only:assert_shape 
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_BATH_FUNCTIONS


  implicit none
  private

  interface ed_chi2_fitgf
     module procedure chi2_fitgf_generic_normal
  end interface ed_chi2_fitgf


  public :: ed_chi2_fitgf


  integer                               :: Ldelta
  complex(8),dimension(:,:),allocatable :: Gdelta
  complex(8),dimension(:,:),allocatable :: Fdelta
  real(8),dimension(:),allocatable      :: Xdelta,Wdelta
  integer                               :: totNorb,totNspin,totNlso
  integer,dimension(:),allocatable      :: getIorb,getJorb,getIspin,getJspin,getIlat,getJlat
  integer                               :: Orb_indx,Spin_indx,Spin_mask
  integer,dimension(:),allocatable      :: Nlambdas
  !location of the maximum of the chisquare over Nlso.
  integer                               :: maxchi_loc
  !
  type nsymm_vector
     real(8),dimension(:),allocatable   :: element          
  end type nsymm_vector
  !

contains


  !##################################################################
  ! THE CALCULATION OF THE \chi^2 FUNCTIONS USE PROCEDURES FURTHER 
  ! BELOW TO EVALUATE INDEPENDENTLY THE ANDERSON MODEL:
  !  - DELTA, 
  !  -\GRAD DELTA
  !  - G0
  ! THE LATTER ARE ADAPTED FROM THE PROCEDURES:
  ! DELTA_BATH_MATS
  ! GRAD_DELTA_BATH_MATS
  ! G0 BATH_MATS
  ! FOR, YOU NEED TO DECOMPOSE THE a INPUT ARRAY INTO ELEMENTS.
  !##################################################################

  !+----------------------------------------------------------------------+
  !PURPOSE  : Chi^2 fit of the G0/Delta 
  !+----------------------------------------------------------------------+
  subroutine chi2_fitgf_generic_normal(fg,bath)
    complex(8),dimension(:,:,:,:,:,:,:) :: fg ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Niw] 
    real(8),dimension(:)                :: bath
    !integer,optional                   :: ispin,iorb
    !integer                            :: ispin_
    !
    !ispin_=1;if(present(ispin))ispin_=ispin
    !
    call assert_shape(fg,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(fg,7)],"chi2_fitgf_generic_normal","fg")
    allocate(Nlambdas(Nbath))
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
    call chi2_fitgf_replica(fg,bath)
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
    deallocate(Nlambdas)
  end subroutine chi2_fitgf_generic_normal



  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_replica(fg,bath_)
    complex(8),dimension(:,:,:,:,:,:,:)                   :: fg ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    logical(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hmask
    real(8),dimension(:),intent(inout)                    :: bath_
    real(8),dimension(:),allocatable                      :: array_bath
    integer                                               :: i,j,ilat,jlat,iorb,jorb,ispin,jspin,ibath,io,jo
    integer                                               :: iter,stride,counter,Asize
    real(8)                                               :: chi
    logical                                               :: check
    character(len=256)                                    :: suffix
    integer                                               :: unit
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_replica error: wrong bath dimensions"
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    allocate(array_bath(size(bath_)-Nbath))
    counter=0
    do ibath=1,Nbath
       counter=counter+1
       Nlambdas(ibath)=bath_(counter)
    enddo
    array_bath=bath_(Nbath+1:size(bath_))
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,7))Ldelta=size(fg,7)
    !
    Hmask=mask_hloc(impHloc,wdiag=.true.,uplo=.true.)
    totNlso=count(Hmask)
    !
    allocate(getIlat(totNlso) ,getJlat(totNlso))
    allocate(getIspin(totNlso),getJspin(totNlso))
    allocate(getIorb(totNlso) ,getJorb(totNlso))
    counter=0
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (Hmask(ilat,jlat,ispin,jspin,iorb,jorb))then
                         counter=counter+1
                         getIlat(counter)  = ilat
                         getIspin(counter) = ispin
                         getIorb(counter)  = iorb
                         getJlat(counter)  = jlat
                         getJspin(counter) = jspin
                         getJorb(counter)  = jorb
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    allocate(Gdelta(totNlso,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    write(LOGfile,*)"  fitted functions",totNlso
    do i=1,totNlso
       Gdelta(i,1:Ldelta) = fg(getIlat(i),getJlat(i),getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
    enddo
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default
       if(cg_grad==0)then
          write(LOGfile,*)"  Using analytic gradient"
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,&
                  chi2_weiss_replica,&
                  grad_chi2_weiss_replica,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,&
                  chi2_delta_replica,&
                  grad_chi2_delta_replica,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
          end select
       else
          write(LOGfile,*)"  Using numerical gradient"
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_replica,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_replica,&
                  iter,&
                  chi, &
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
          end select
       endif
       !
    case (1)
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,chi2_weiss_replica,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_replica,&
               iter,&
               chi, &
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,      &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5,A)")"chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  All Orbs, All Spins"
    !
    suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    bath_(Nbath+1:size(bath_))=array_bath
    call set_dmft_bath(bath_)           ! *** bath_ --> dmft_bath ***    (per write fit result)
    call write_dmft_bath(LOGfile)
    call save_dmft_bath()
    !
    call write_fit_result()
    !
    call get_dmft_bath(bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath()
    deallocate(Gdelta,Xdelta,Wdelta)
    deallocate(getIlat,getJlat)
    deallocate(getIspin,getJspin)
    deallocate(getIorb,getJorb)
    deallocate(array_bath)
    !
  contains
    !
    subroutine write_fit_result()
      complex(8)        :: fgand(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)
      integer           :: i,j,s,l,ilat,jlat,iorb,jorb,ispin,jspin
      real(8)           :: w
      if(cg_scheme=='weiss')then
         fgand(:,:,:,:,:,:,:) = g0and_bath(xi*Xdelta(:))
      else
         fgand(:,:,:,:,:,:,:) = delta_bath(xi*Xdelta(:))
      endif
      !
      do l=1,totNlso
         ilat = getIlat(l)
         jlat = getJlat(l)
         iorb = getIorb(l)
         jorb = getJorb(l)
         ispin = getIspin(l)
         jspin = getJspin(l)
         suffix="_i"//reg(txtfy(ilat))//&
              "_j"//reg(txtfy(jlat))//&
              "_l"//reg(txtfy(iorb))//&
              "_m"//reg(txtfy(jorb))//&
              "_s"//reg(txtfy(ispin))//&
              "_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
         unit=free_unit()
         if(cg_scheme=='weiss')then
            open(unit,file="fit_weiss"//reg(suffix)//".ed")
         else
            open(unit,file="fit_delta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            w = Xdelta(i)
            write(unit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(ilat,jlat,ispin,jspin,iorb,jorb,i)),dimag(fgand(ilat,jlat,ispin,jspin,iorb,jorb,i)),&
                 dreal(fg(ilat,jlat,ispin,jspin,iorb,jorb,i)),dreal(fgand(ilat,jlat,ispin,jspin,iorb,jorb,i))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_replica


  !##################################################################
  ! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
  !##################################################################
  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  !+-------------------------------------------------------------+

  function chi2_delta_replica(a) result(chi2)
    real(8),dimension(:)                                         :: a
    real(8)                                                      :: chi2
    real(8),dimension(totNlso)                                   :: chi2_lso
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    real(8),dimension(Ldelta)                                    :: Ctmp
    integer                                                      :: i,l,ilat,jlat,iorb,jorb,ispin,jspin
    !
    Delta = delta_replica(a)
    !
    do l=1,totNlso
       ilat = getIlat(l)
       jlat = getJlat(l)
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ctmp =  abs(Gdelta(l,:)-Delta(ilat,jlat,ispin,jspin,iorb,jorb,:))
       chi2_lso(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_lso)
    chi2=chi2/Ldelta
    !
  end function chi2_delta_replica

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of 
  ! \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function grad_chi2_delta_replica(a) result(dchi2)
    real(8),dimension(:)                                                 :: a
    real(8),dimension(size(a))                                           :: dchi2
    real(8),dimension(totNlso,size(a))                                   :: df
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)         :: Delta
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    complex(8),dimension(Ldelta)                                         :: Ftmp
    real(8),dimension(Ldelta)                                            :: Ctmp
    integer                                                              :: i,j,l,ilat,jlat,iorb,jorb,ispin,jspin
    !
    Delta  = delta_replica(a)
    dDelta = grad_delta_replica(a)
    !
    do l=1,totNlso
       ilat = getIlat(l)
       jlat = getJlat(l)
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ftmp = Gdelta(l,:)-Delta(ilat,jlat,ispin,jspin,iorb,jorb,:)
       Ctmp = abs(Ftmp)**(cg_pow-2)
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Ftmp)*dreal(dDelta(ilat,jlat,ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dDelta(ilat,jlat,ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta )
       enddo
    enddo
    !
    dchi2 = -cg_pow*sum(df,1)/Ldelta     !sum over all orbital indices
    print*,dchi2
    !
  end function grad_chi2_delta_replica


  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient. 
  !+-------------------------------------------------------------+
  function chi2_weiss_replica(a) result(chi2)
    real(8),dimension(:)                                         :: a
    real(8)                                                      :: chi2
    real(8),dimension(totNlso)                                   :: chi2_lso
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
    real(8),dimension(Ldelta)                                    :: Ctmp
    integer                                                      :: i,l,ilat,jlat,iorb,jorb,ispin,jspin
    !
    g0and = g0and_replica(a)
    !
    do l=1,totNlso
       ilat = getIlat(l)
       jlat = getJlat(l)
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ctmp = abs(Gdelta(l,:)-g0and(ilat,jlat,ispin,jspin,iorb,jorb,:))
       chi2_lso(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    !FIXME:THIS NEEDS A THOROUGH DISCUSSION
    !
    !chi2=sum(chi2_lso)
    chi2=maxval(chi2_lso)
    maxchi_loc=maxloc(chi2_lso,1)
    chi2=chi2/Ldelta
    !
  end function chi2_weiss_replica

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of 
  ! \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function grad_chi2_weiss_replica(a) result(dchi2)
    real(8),dimension(:)                                                 :: a
    real(8),dimension(size(a))                                           :: dchi2
    real(8),dimension(totNlso,size(a))                                   :: df
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)         :: g0and
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dg0and
    complex(8),dimension(Ldelta)                                         :: Ftmp
    real(8),dimension(Ldelta)                                            :: Ctmp
    integer                                                              :: i,j,l,ilat,jlat,iorb,jorb,ispin,jspin
    !
    g0and  = g0and_replica(a)
    dg0and = grad_g0and_replica(a)
    !
    do l=1,totNlso
       ilat = getIlat(l)
       jlat = getJlat(l)
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ftmp = Gdelta(l,:)-g0and(ilat,jlat,ispin,jspin,iorb,jorb,:)
       Ctmp = abs(Ftmp)**(cg_pow-2)
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Ftmp)*dreal(dg0and(ilat,jlat,ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dg0and(ilat,jlat,ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta )
       enddo
    enddo
    !
    !dchi2 = -cg_pow*sum(df,1)     !sum over all orbital indices
    dchi2 = -cg_pow*df(maxchi_loc,:)
    dchi2 = dchi2/Ldelta
    !
  end function grad_chi2_weiss_replica

  !##################################################################
  ! THESE PROCEDURES EVALUATE THE 
  ! - \delta
  ! - g0
  ! FUNCTIONS. 
  !##################################################################
  function delta_replica(a) result(Delta)
    real(8),dimension(:)                                          :: a
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)  :: Delta
    integer                                                       :: ilat,jlat,ispin,jspin,iorb,jorb,ibath,isym
    integer                                                       :: i,io,jo,ndx,stride
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: V_k
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: Haux
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: invH_knnn
    real(8),dimension(Nbath)                                      :: dummy_Vbath
    type(nsymm_vector),dimension(Nbath)                           :: dummy_lambda
    complex(8)                                                    :: iw
    !
    !ACHTUNG! here the bath was a temporary one, since we removed the possibility to act on other baths we need to replicate the
    !function behaviour. Rather ugly...
    !Get Hs
    stride = 0
    do ibath=1,Nbath
       allocate(dummy_lambda(ibath)%element(Nlambdas(ibath)))
       !Get Vs
       stride = stride + 1
       dummy_vbath(ibath) = a(stride)
       !get Lambdas
       dummy_lambda(ibath)%element=a(stride+1:stride+Nlambdas(ibath))
       stride=stride+Nlambdas(ibath)
    enddo
    !
    !
    Delta=zero
    do i=1,Ldelta
       iw = xi*Xdelta(i)+xmu
       do ibath=1,Nbath
          invH_knnn  = bath_from_sym(dummy_lambda(ibath)%element)
          Haux      = zeye(Nlat*Nspin*Norb)*iw - nnn2lso_reshape(invH_knnn,Nlat,Nspin,Norb)
          call inv(Haux)
          invH_knnn = lso2nnn_reshape(Haux,Nlat,Nspin,Norb)
          Delta(:,:,:,:,:,:,i)=Delta(:,:,:,:,:,:,i)+ dummy_Vbath(ibath)*invH_knnn(:,:,:,:,:,:)*dummy_Vbath(ibath)
       enddo
    enddo
    !
    do ibath=1,Nbath
       deallocate(dummy_lambda(ibath)%element)
    enddo
    !
  end function delta_replica

  function g0and_replica(a) result(G0and)
    real(8),dimension(:)                                         :: a
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
    complex(8),dimension(Nlat*Norb*Nspin,Nlat*Norb*Nspin)        :: zeta,fgorb
    integer                                                      :: i,Nlso
    integer                                                      :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    Nlso = Nlat*Norb*Nspin
    !
    Delta = delta_replica(a)
    !
    do i=1,Ldelta
       zeta  = (xi*Xdelta(i)+xmu)*zeye(Nlso)
       FGorb = zeta - nnn2lso_reshape(impHloc + Delta(:,:,:,:,:,:,i), Nlat,Nspin,Norb)
       call inv(FGorb)
       G0and(:,:,:,:,:,:,i) = lso2nnn_reshape(FGorb,Nlat,Nspin,Norb)
    enddo
    !
  end function g0and_replica

  !##################################################################
  ! THESE PROCEDURES EVALUATE GRADIENT OF THE 
  ! - \delta
  ! - g0
  ! FUNCTIONS. 
  !##################################################################

  function grad_delta_replica(a) result(dDelta)
    real(8),dimension(:)                                                 :: a
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    integer                                                              :: ilat,jlat,ispin,iorb,jorb,ibath
    integer                                                              :: k,l,io,counter
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                :: H_reconstructed, Htmp,Hbasis_lso
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Ldelta)         :: Haux
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)         :: invH_knn
    real(8),dimension(1,Nbath)                                           :: dummy_Vbath !FIXME: TO EXTEND: 1->NSPIN
    type(nsymm_vector),dimension(Nbath)                                  :: dummy_lambda
    !
    !
    !Get Hs
    counter = 0
    do ibath=1,Nbath
       allocate(dummy_lambda(ibath)%element(Nlambdas(ibath)))
       !
       !FIXME: to extend uncomment Nspin and 1->NSPIN
       !do ispin=1,Nspin
       counter = counter + 1
       dummy_vbath(1,ibath) = a(counter)
       !enddo
       !
       dummy_lambda(ibath)%element=a(counter+1:counter+Nlambdas(ibath))
       counter=counter+Nlambdas(ibath)
    enddo
    !
    dDelta=zero
    counter=0
    !
    do ibath=1,Nbath
       H_reconstructed = nnn2lso_reshape(bath_from_sym(dummy_lambda(ibath)%element),Nlat,Nspin,Norb)
       do l=1,Ldelta
          Haux(:,:,l) = zeye(Nlat*Nspin*Norb)*(xi*Xdelta(l)+xmu) - H_reconstructed
          call inv(Haux(:,:,l))
          invH_knn(:,:,:,:,:,:,l) = lso2nnn_reshape(Haux(:,:,l),Nlat,Nspin,Norb)
       enddo
       !Derivate_Vp
       do ispin=1,Nspin
          counter = counter + 1
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      !FIXME: to extend, 1->ISPIN
                      dDelta(ilat,jlat,ispin,ispin,iorb,jorb,:,counter)=2d0*dummy_Vbath(1,ibath)*invH_knn(ilat,jlat,ispin,ispin,iorb,jorb,:)
                   enddo
                enddo
             enddo
          enddo
       enddo
       !Derivate_lambda_p
       do k=1,size(dummy_lambda(ibath)%element)
          counter = counter + 1
          Hbasis_lso=nnn2lso_reshape(H_basis(k)%O,Nlat,Nspin,Norb)
          do l=1,Ldelta
             ! Hbasis_lso=nnn2lso_reshape(H_basis(k)%O,Nlat,Nspin,Norb)
             ! Htmp=matmul(Haux(:,:,l),Hbasis_lso)
             ! Htmp=matmul(Htmp,Haux(:,:,l))
             Htmp = ((Haux(:,:,l) .x. Hbasis_lso)) .x. Haux(:,:,l)
             dDelta(:,:,:,:,:,:,l,counter)=lso2nnn_reshape((dummy_Vbath(1,ibath)**2)*Htmp,Nlat,Nspin,Norb)
          enddo
       enddo
    enddo
  end function grad_delta_replica


  function grad_g0and_replica(a) result(dG0and)
    real(8),dimension(:)                                                 :: a
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dG0and,dDelta
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)         :: G0and
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                :: dDelta_lso,dG0and_lso,G0and_lso
    integer                                                              :: ilat,jlat,ispin,iorb,jorb
    integer                                                              :: ik,l
    !
    G0and  = g0and_replica(a)
    dDelta = grad_delta_replica(a)
    !
    dG0and = zero
    !
    do l=1,Ldelta
       G0and_lso=nnn2lso_reshape(g0and(:,:,:,:,:,:,l),Nlat,Nspin,Norb)
       do ik=1,size(a)
          dDelta_lso=nnn2lso_reshape(dDelta(:,:,:,:,:,:,l,ik),Nlat,Nspin,Norb)
          ! dG0and_lso=matmul(-G0and_lso,dDelta_lso)
          ! dG0and_lso=matmul(dG0and_lso,G0and_lso)
          dG0and_lso = (G0and_lso .x. dDelta_lso) .x. G0and_lso
          dG0and(:,:,:,:,:,:,l,ik)=-lso2nnn_reshape(dG0and_lso,Nlat,Nspin,Norb)
       enddo
    enddo
    !
  end function grad_g0and_replica


end MODULE ED_FIT_CHI2
