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
   integer,dimension(:),allocatable                      :: Nlambdas
   !location of the maximum of the chisquare over Nlso.
   integer                                               :: maxchi_loc
   !
   type nsymm_vector
      real(8),dimension(:),allocatable                   :: element
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



   !+-------------------------------------------------------------+
   !PURPOSE  : Chi^2 interface for REPLICA BATH
   !+-------------------------------------------------------------+
   subroutine chi2_fitgf_replica(fg,bath_)
      complex(8),dimension(:,:,:,:,:,:,:)                   :: fg ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
      real(8),dimension(:),intent(inout)                    :: bath_
      real(8),dimension(:),allocatable                      :: array_bath
      integer                                               :: i,j,ilat,jlat,iorb,jorb,ispin,jspin,ibath,io,jo
      integer                                               :: iter,stride,counter,Asize
      real(8)                                               :: chi
      logical                                               :: check
      logical                                               :: ed_all_g ! TO BE MOVED TO INPUT VARS (or pruned)
      character(len=256)                                    :: suffix
      integer                                               :: unit
      !
      if(size(fg,1)/=Nlat) stop "chi2_fitgf_replica error: size[fg,1]!=Nlat"
      if(size(fg,2)/=Nlat) stop "chi2_fitgf_replica error: size[fg,2]!=Nlat"
      if(size(fg,3)/=Nspin)stop "chi2_fitgf_replica error: size[fg,3]!=Nspin"
      if(size(fg,4)/=Nspin)stop "chi2_fitgf_replica error: size[fg,4]!=Nspin"
      if(size(fg,5)/=Norb) stop "chi2_fitgf_replica error: size[fg,5]!=Norb"
      if(size(fg,6)/=Norb) stop "chi2_fitgf_replica error: size[fg,6]!=Norb"
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
         Nlambdas(ibath)=NINT(bath_(counter))
      enddo
      array_bath=bath_(Nbath+1:size(bath_))
      !
      Ldelta = Lfit ; if(Ldelta>size(fg,7))Ldelta=size(fg,7)
      !
      !
      allocate(FGmatrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta))
      allocate(Hmask(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
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
      !
      !----------------------------------------------------------------------------------------
      !BORROWED FROM LIB_DMFT_ED
      !> https://github.com/QcmPlab/LIB_DMFT_ED/commit/0e5c272b45eda6b7ff652e2473b9ecda09e5ba8b
      ed_all_g = .true. !(dummy, we assume ed_all_g always set here, for now)
      if(ed_all_g)then
         Hmask=.true.
         ! Aren't we sure about hermiticity?
         ! -> Hmask=Hreplica_mask(wdiag=.false.,uplo=.true.)
      else
         Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
      endif
      ! For now I'd say that we'd better dump everything inside the \chi^2, hence Hmask=.true.,
      ! but we might want to consider exploiting hermiticity at least (probably not looking at
      ! Hreplica though: who guarantees zeros therein imply zeros here?). Discussion needed...
      !----------------------------------------------------------------------------------------
      !
      !
      FGmatrix = fg
      !
      !
      select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
       case default
         if(cg_grad==0)then
#if __GNUC__ >= 8 || __INTEL_COMPILER >= 1500
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
#else
            STOP "analytic gradient not supported for gfortran < 8"
#endif
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
         if(cg_grad==0)then
            write(*,*) "                                                                                "
            write(*,*) "WARNING: analytic gradient not available with cg-method=1 (minimize f77 routine)"
            write(*,*) "         > we will force cg_grad=1 (so let the routine estimate the gradient)   "
            write(*,*) "                                                                                "
         endif
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
      call set_dmft_bath(bath_)               ! *** bath_ --> dmft_bath ***    (per write fit result)
      call write_dmft_bath(LOGfile)
      call save_dmft_bath()
      !
      call write_fit_result()
      !
      call get_dmft_bath(bath_)               ! ***  dmft_bath --> bath_ ***    (bath in output)
      call deallocate_dmft_bath()
      deallocate(FGmatrix,Xdelta,Wdelta)
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
         do ilat=1,Nlat
            do jlat=1,Nlat
               do ispin=1,Nspin
                  do jspin=1,Nspin
                     do iorb=1,Norb
                        do jorb=1,Norb
                           suffix="_i"//reg(str(ilat))//&
                              "_j"//reg(str(jlat))//&
                              "_l"//reg(str(iorb))//&
                              "_m"//reg(str(jorb))//&
                              "_s"//reg(str(ispin))//&
                              "_r"//reg(str(jspin))//reg(ed_file_suffix)
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
                     enddo
                  enddo
               enddo
            enddo
         enddo
      end subroutine write_fit_result
      !
   end subroutine chi2_fitgf_replica












   !##################################################################
   ! THESE PROCEDURES EVALUATE THE \chi^2 FUNCTIONS TO MINIMIZE.
   !##################################################################
   !
   !
   !+-----------------------------------------------------------------+
   !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
   !+-----------------------------------------------------------------+
   function chi2_delta_replica(a) result(chi2)
      real(8),dimension(:)                                         :: a
      real(8)                                                      :: chi2
      !
      select case(cg_norm)
       case ("elemental")
         chi2 = chi2_delta_replica_elemental(a)
       case ("frobenius")
         chi2 = chi2_delta_replica_frobenius(a)
       case default
         stop "chi2_fitgf_replica error: cg_norm != [elemental,frobenius]"
      end select
      !
   end function chi2_delta_replica
   !
   !
   !> ELEMENTAL NORM: weighted sum over i\omega for each matrix element, then weighted sum over elements
   function chi2_delta_replica_elemental(a) result(chi2)
      real(8),dimension(:)                                          :: a
      real(8)                                                       :: chi2
      real(8),dimension(Ldelta)                                     :: chi2_freq
      real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)            :: chi2_mtrx,Wmat
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)  :: Delta
      integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin
      !
#ifdef _DEBUG
      if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_delta_replica_elemental. a:",a
#endif
      !
      Delta = delta_replica(a)
      !
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        chi2_freq = abs(Delta(ilat,jlat,ispin,jspin,iorb,jorb,1:Ldelta) - FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,1:Ldelta))
                        chi2_mtrx(ilat,jlat,ispin,jspin,iorb,jorb) = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
                        select case(cg_matrix)
                         case(0) !FLAT (all matrix elements weighted equal)
                           Wmat(ilat,jlat,ispin,jspin,iorb,jorb) = 1d0!0.25d0 !needs to depend on the hopping I think (\Delta=(D/2)^2*Gloc…)
                         case(1) !SPECTRAL (normalization through A(iw), element by element)
                           Wmat(ilat,jlat,ispin,jspin,iorb,jorb) = -sum(dimag(FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,1:Lmats)))/beta
#ifdef _DEBUG
                           if(ed_verbose>4)then
                              print*, "ilat: "//str(ilat)//" ispin: "//str(ispin)//" iorb: "//str(iorb)
                              print*, "jlat: "//str(jlat)//" jspin: "//str(jspin)//" jorb: "//str(jorb)
                              print*, "> Wmat: ", Wmat(ilat,jlat,ispin,jspin,iorb,jorb)/0.25d0
                           endif
#endif
                        end select
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      chi2 = sum(chi2_mtrx / Wmat, Hmask) !Weighted sum over matrix elements
      chi2 = chi2 / Ldelta / count(Hmask) !Normalization over {iw} and Hmask
      !
#ifdef _DEBUG
      if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_replica_elemental. chi2:",chi2
#endif
      !
   end function chi2_delta_replica_elemental
   !
   !
   !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
   function chi2_delta_replica_frobenius(a) result(chi2)
      real(8),dimension(:)                                         :: a
      real(8)                                                      :: chi2
      real(8),dimension(Ldelta)                                    :: chi2_freq
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)        :: Delta_lso
      integer                                                      :: l
      !
#ifdef _DEBUG
      if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_delta_replica_frobenius. a:",a
#endif
      !
      Delta = delta_replica(a)
      !
      do l=1,Ldelta
         Delta_lso    =  nnn2lso_reshape(delta(:,:,:,:,:,:,l) - FGmatrix(:,:,:,:,:,:,l),Nlat,Nspin,Norb)
         chi2_freq(l) =  sqrt(trace(matmul(Delta_lso,conjg(transpose(Delta_lso)))))
      enddo
      !
      chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
      chi2 = chi2/Ldelta/(Nlat*Nspin*Norb) !Normalization over {iw} and Nlso
      !
#ifdef _DEBUG
      if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_replica_frobenius. chi2:",chi2
#endif
      !
   end function chi2_delta_replica_frobenius
   !
   !
   !+-------------------------------------------------------------+
   !PURPOSE: Evaluate the \chi^2 distance of G0_Anderson function.
   !+-------------------------------------------------------------+
   function chi2_weiss_replica(a) result(chi2)
      real(8),dimension(:)                                         :: a
      real(8)                                                      :: chi2
      !
      select case(cg_norm)
       case ("elemental")
         chi2 = chi2_weiss_replica_elemental(a)
       case ("frobenius")
         chi2 = chi2_weiss_replica_frobenius(a)
       case default
         stop "chi2_fitgf_replica error: cg_norm != [elemental,frobenius]"
      end select
      !
   end function chi2_weiss_replica
   !
   !
   !> ELEMENTAL NORM: weighted sum over i\omega for each matrix element, then weighted sum over elements
   function chi2_weiss_replica_elemental(a) result(chi2)
      real(8),dimension(:)                                          :: a
      real(8)                                                       :: chi2
      real(8),dimension(Ldelta)                                     :: chi2_freq
      real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)            :: chi2_mtrx,Wmat
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)  :: g0and
      integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin
      !
#ifdef _DEBUG
      if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_replica_elemental. a:",a
#endif
      !
      g0and = g0and_replica(a)
      !
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        chi2_freq = abs(g0and(ilat,jlat,ispin,jspin,iorb,jorb,1:Ldelta) - FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,1:Ldelta))
                        chi2_mtrx(ilat,jlat,ispin,jspin,iorb,jorb) = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
                        select case(cg_matrix)
                         case(0) !FLAT (all matrix elements weighted equal)
                           Wmat(ilat,jlat,ispin,jspin,iorb,jorb) = 1d0
                         case(1) !SPECTRAL (normalization through A(iw), element by element)
                           Wmat(ilat,jlat,ispin,jspin,iorb,jorb) = -sum(dimag(FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,1:Lmats)))/beta
#ifdef _DEBUG
                           if(ed_verbose>4)then
                              print*, "ilat: "//str(ilat)//"ispin: "//str(ispin)//"iorb: "//str(iorb)
                              print*, "jlat: "//str(jlat)//"jspin: "//str(jspin)//"jorb: "//str(jorb)
                              print*, "> Wmat = ", Wmat(ilat,jlat,ispin,jspin,iorb,jorb)
                           endif
#endif
                        end select
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      chi2 = sum(chi2_mtrx / Wmat, Hmask) !Weighted sum over matrix elements
      chi2 = chi2 / Ldelta / count(Hmask) !Normalization over {iw} and Hmask
      !
#ifdef _DEBUG
      if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_replica_elemental. chi2:",chi2
#endif
      !
   end function chi2_weiss_replica_elemental
   !
   !
   !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
   function chi2_weiss_replica_frobenius(a) result(chi2)
      real(8),dimension(:)                                         :: a
      real(8)                                                      :: chi2
      real(8),dimension(Ldelta)                                    :: chi2_freq
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)        :: Delta_lso
      integer                                                      :: l
      !
#ifdef _DEBUG
      if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_replica_frobenius. a:",a
#endif
      !
      g0and = g0and_replica(a)
      !
      do l=1,Ldelta
         Delta_lso    =  nnn2lso_reshape(g0and(:,:,:,:,:,:,l) - FGmatrix(:,:,:,:,:,:,l),Nlat,Nspin,Norb)
         chi2_freq(l) =  sqrt(trace(matmul(Delta_lso,conjg(transpose(Delta_lso)))))
      enddo
      !
      chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
      chi2 = chi2/Ldelta/(Nlat*Nspin*Norb) !Normalization over {iw} and Nlso
      !
#ifdef _DEBUG
      if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_replica_frobenius. chi2:",chi2
#endif
      !
   end function chi2_weiss_replica_frobenius







   !######################################################################
   ! THESE PROCEDURES EVALUATE THE >GRADIENTS< OF THE \chi^2 TO MINIMIZE.
   !######################################################################
   !
   !
#if __GNUC__ >= 8 || __INTEL_COMPILER >= 1500
   !+----------------------------------------------------------------------+
   !PURPOSE: Evaluate the gradient \Grad\chi^2 of \Delta_Anderson function.
   !+----------------------------------------------------------------------+
   function grad_chi2_delta_replica(a) result(dchi2)
      real(8),dimension(:)                                         :: a
      real(8),dimension(size(a))                                   :: dchi2
      !
      select case(cg_norm)
       case ("elemental")
         dchi2 = grad_chi2_delta_replica_elemental(a)
       case ("frobenius")
         dchi2 = grad_chi2_delta_replica_frobenius(a)
       case default
         stop "chi2_fitgf_replica error: cg_norm != [elemental,frobenius]"
      end select
      !
   end function grad_chi2_delta_replica
   !
   !
   !> ELEMENTAL NORM: weighted sum over i\omega for each matrix element, then weighted sum over elements
   function grad_chi2_delta_replica_elemental(a) result(dchi2)
      real(8),dimension(:)                                                  :: a
      real(8),dimension(size(a))                                            :: dchi2
      real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(a))            :: df
      real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                    :: Wmat
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)          :: Delta
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a))  :: dDelta
      complex(8),dimension(Ldelta)                                          :: Ftmp
      real(8),dimension(Ldelta)                                             :: Ctmp
      integer                                                               :: ia,ilat,jlat,iorb,jorb,ispin,jspin
      !
#ifdef _DEBUG
      if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_replica_elemental. a:",a
#endif
      !
      Delta  = delta_replica(a)
      dDelta = grad_delta_replica(a)
      !
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        Ftmp = delta(ilat,jlat,ispin,jspin,iorb,jorb,1:Ldelta) - FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,1:Ldelta)
                        Ctmp = abs(Ftmp)**(cg_pow-2)
                        do ia = 1,size(a)
                           df(ilat,jlat,ispin,jspin,iorb,jorb,ia) = & !Weighted sum over matsubara frqs
                              sum( dreal(Ftmp) * dreal(dDelta(ilat,jlat,ispin,jspin,iorb,jorb,:,ia)) * Ctmp/Wdelta ) + &
                              sum( dimag(Ftmp) * dimag(dDelta(ilat,jlat,ispin,jspin,iorb,jorb,:,ia)) * Ctmp/Wdelta )
                        enddo
                        select case(cg_matrix)
                         case(0) !FLAT (all matrix elements weighted equal)
                           Wmat(ilat,jlat,ispin,jspin,iorb,jorb) = 1d0!0.25d0 !needs to depend on the hopping I think (\Delta=(D/2)^2*Gloc…)
                         case(1) !SPECTRAL (normalization through A(iw), element by element)
                           Wmat(ilat,jlat,ispin,jspin,iorb,jorb) = -sum(dimag(FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,1:Lmats)))/beta
                        end select
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      do ia=1,size(a)
         dchi2(ia) = + cg_pow * sum( df(:,:,:,:,:,:,ia) / Wmat, Hmask) !Weighted sum over matrix elements
         dchi2(ia) = dchi2(ia) / Ldelta / count(Hmask) !Normalization over {iw} and Hmask
      enddo
      !
#ifdef _DEBUG
      if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_replica_elemental. dChi2:",dchi2
#endif
      !
   end function grad_chi2_delta_replica_elemental
   !
   !
   !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
   function grad_chi2_delta_replica_frobenius(a) result(dchi2)
      real(8),dimension(:)                                                 :: a
      real(8),dimension(size(a))                                           :: dchi2
      real(8),dimension(Ldelta,size(a))                                    :: df
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)         :: Delta
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
      complex(8),dimension(Ldelta)                                         :: Ftmp
      real(8),dimension(Ldelta,size(a))                                    :: dChi_freq
      integer                                                              :: i,j,idelta,ilat,jlat,iorb,jorb,ispin,jspin
      !
      Delta  = delta_replica(a)
      dDelta = grad_delta_replica(a)
      Ftmp=zero
      df=zero
      !
      do idelta=1,Ldelta
         do ilat=1,Nlat
            do jlat=1,Nlat
               do ispin=1,Nspin
                  do jspin=1,Nspin
                     do iorb=1,Norb
                        do jorb=1,Norb
                           !
                           Ftmp(idelta) = Ftmp(idelta) + abs(Delta(ilat,jlat,ispin,jspin,iorb,jorb,idelta)-FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,idelta))**2
                           do j=1,size(a)
                              df(idelta,j) = df(idelta,j) + &
                                 real(Delta(ilat,jlat,ispin,jspin,iorb,jorb,idelta) - FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,idelta)) * &
                                 real(dDelta(ilat,jlat,ispin,jspin,iorb,jorb,idelta,j)) + &
                                 imag(Delta(ilat,jlat,ispin,jspin,iorb,jorb,idelta) - FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,idelta)) * &
                                 imag(dDelta(ilat,jlat,ispin,jspin,iorb,jorb,idelta,j))
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         Ftmp(idelta) = cg_pow * (sqrt(Ftmp(idelta))**(cg_pow-2)) / Wdelta(idelta)
         dchi_freq(idelta,:) = Ftmp(idelta) * df(idelta,:)
      enddo
      !
      dchi2 = sum(dchi_freq,1)/Ldelta/(Nlat*Nspin*Norb)
      !
#ifdef _DEBUG
      if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_replica_frobenius. dChi2:",dchi2
#endif
      !
   end function grad_chi2_delta_replica_frobenius
   !
   !
   !+------------------------------------------------------------------+
   !PURPOSE: Evaluate the gradient \Grad\chi^2 of G0_Anderson function.
   !+------------------------------------------------------------------+
   function grad_chi2_weiss_replica(a) result(dchi2)
      real(8),dimension(:)                                         :: a
      real(8),dimension(size(a))                                   :: dchi2
      !
      select case(cg_norm)
       case ("elemental")
         dchi2 = grad_chi2_weiss_replica_elemental(a)
       case ("frobenius")
         dchi2 = grad_chi2_weiss_replica_frobenius(a)
       case default
         stop "chi2_fitgf_replica error: cg_norm != [elemental,frobenius]"
      end select
      !
   end function grad_chi2_weiss_replica
   !
   !
   !> ELEMENTAL NORM: weighted sum over i\omega for each matrix element, then weighted sum over elements
   function grad_chi2_weiss_replica_elemental(a) result(dchi2)
      real(8),dimension(:)                                                  :: a
      real(8),dimension(size(a))                                            :: dchi2
      real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(a))            :: df
      real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                    :: Wmat
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)          :: g0and
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a))  :: dg0and
      complex(8),dimension(Ldelta)                                          :: Ftmp
      real(8),dimension(Ldelta)                                             :: Ctmp
      integer                                                               :: ia,ilat,jlat,iorb,jorb,ispin,jspin
      !
#ifdef _DEBUG
      if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_replica_elemental. a:",a
#endif
      !
      g0and  = g0and_replica(a)
      dg0and = grad_g0and_replica(a)
      !
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        Ftmp = g0and(ilat,jlat,ispin,jspin,iorb,jorb,1:Ldelta) - FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,1:Ldelta)
                        Ctmp = abs(Ftmp)**(cg_pow-2)
                        do ia = 1,size(a)
                           df(ilat,jlat,ispin,jspin,iorb,jorb,ia) = & !Weighted sum over matsubara frqs
                              sum( dreal(Ftmp) * dreal(dg0and(ilat,jlat,ispin,jspin,iorb,jorb,:,ia)) * Ctmp/Wdelta ) + &
                              sum( dimag(Ftmp) * dimag(dg0and(ilat,jlat,ispin,jspin,iorb,jorb,:,ia)) * Ctmp/Wdelta )
                        enddo
                        select case(cg_matrix)
                         case(0) !FLAT (all matrix elements weighted equal)
                           Wmat(ilat,jlat,ispin,jspin,iorb,jorb) = 1d0
                         case(1) !SPECTRAL (normalization through A(iw), element by element)
                           Wmat(ilat,jlat,ispin,jspin,iorb,jorb) = -sum(dimag(FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,1:Lmats)))/beta
                        end select
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      do ia=1,size(a)
         dchi2(ia) = + cg_pow * sum( df(:,:,:,:,:,:,ia) / Wmat, Hmask) !Weighted sum over matrix elements
         dchi2(ia) = dchi2(ia) / Ldelta / count(Hmask) !Normalization over {iw} and Hmask
      enddo
      !
#ifdef _DEBUG
      if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_replica_elemental. dChi2:",dchi2
#endif
      !
   end function grad_chi2_weiss_replica_elemental
   !
   !
   !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
   function grad_chi2_weiss_replica_frobenius(a) result(dchi2)
      real(8),dimension(:)                                                 :: a
      real(8),dimension(size(a))                                           :: dchi2
      real(8),dimension(Ldelta,size(a))                                    :: df
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)         :: g0and
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dg0and
      complex(8),dimension(Ldelta)                                         :: Ftmp
      real(8),dimension(Ldelta,size(a))                                    :: dChi_freq
      integer                                                              :: i,j,idelta,ilat,jlat,iorb,jorb,ispin,jspin
      !
      !
      print*, "                                                                     "
      print*, "WARNING: the analytic gradient of the Weiss field Frobenius distance "
      print*, "         has been found giving /QUALITATIVELY WRONG/ fitted spectra! "
      print*, "         > the issue has emerged in some Nlat=Nspin=Norb=1 test runs "
      print*, "         > while the bug is investigated please switch cg_scheme to  "
      print*, "           'delta' or cg_norm to 'elemental' (or give up analytic cg)"
      print*, "                                                                     "
      !
      !
      g0and  = g0and_replica(a)
      dg0and = grad_g0and_replica(a)
      Ftmp=zero
      df=zero
      !
      do idelta=1,Ldelta
         do ilat=1,Nlat
            do jlat=1,Nlat
               do ispin=1,Nspin
                  do jspin=1,Nspin
                     do iorb=1,Norb
                        do jorb=1,Norb
                           !
                           Ftmp(idelta) = Ftmp(idelta) + abs(g0and(ilat,jlat,ispin,jspin,iorb,jorb,idelta)-FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,idelta))**2
                           do j=1,size(a)
                              df(idelta,j) = df(idelta,j) + &
                                 real(g0and(ilat,jlat,ispin,jspin,iorb,jorb,idelta) - FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,idelta)) * &
                                 real(dg0and(ilat,jlat,ispin,jspin,iorb,jorb,idelta,j)) + &
                                 imag(g0and(ilat,jlat,ispin,jspin,iorb,jorb,idelta) - FGmatrix(ilat,jlat,ispin,jspin,iorb,jorb,idelta)) * &
                                 imag(dg0and(ilat,jlat,ispin,jspin,iorb,jorb,idelta,j))
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         Ftmp(idelta) = cg_pow * (sqrt(Ftmp(idelta))**(cg_pow-2)) / Wdelta(idelta)
         dchi_freq(idelta,:) = Ftmp(idelta) * df(idelta,:)
      enddo
      !
      dchi2 = sum(dchi_freq,1)/Ldelta/(Nlat*Nspin*Norb)
      !
#ifdef _DEBUG
      if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_replica_frobenius. dChi2:",dchi2
#endif
      !
   end function grad_chi2_weiss_replica_frobenius
#endif












   !##################################################################
   ! THESE PROCEDURES EVALUATE THE ANDERSON FUNCTIONS:
   ! - \Delta (hybridization)
   ! - g0     (weiss field)
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
            invH_knnn  = Hreplica_build(dummy_lambda(ibath)%element)
            Haux      = zeye(Nlat*Nspin*Norb)*iw - nnn2lso_reshape(invH_knnn,Nlat,Nspin,Norb)
            call inv(Haux) !GUARDA MAIL: FORSE BASTA INVERTIRE UNA VOLTA -> U(iwI-D)U†
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
   !
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
   !
   !
   !
   !
   !
   !
#if __GNUC__ >= 8 || __INTEL_COMPILER >= 1500
   !##################################################################
   ! THESE PROCEDURES EVALUATE THE GRADIENT OF ANDERSON FUNCTIONS:
   ! - \Delta (hybridization)
   ! - g0     (weiss field)
   !##################################################################
   function grad_delta_replica(a) result(dDelta)
      real(8),dimension(:)                                                 :: a
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
      integer                                                              :: ilat,jlat,ispin,iorb,jorb,ibath
      integer                                                              :: k,l,io,counter
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                :: H_reconstructed, Htmp,Hbasis_lso
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Ldelta)         :: Haux
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Ldelta)         :: invH_knn
      real(8),dimension(1,Nbath)                                           :: dummy_Vbath !TODO) to extend to Vup != Vdw: 1->NSPIN
      type(nsymm_vector),dimension(Nbath)                                  :: dummy_lambda
      !
      !
      !Get Hs
      counter = 0
      do ibath=1,Nbath
         if(allocated(dummy_lambda(ibath)%element))deallocate(dummy_lambda(ibath)%element)
         allocate(dummy_lambda(ibath)%element(Nlambdas(ibath)))
         !
         counter = counter + 1
         do ispin=1,Nspin
            dummy_vbath(1,ibath) = a(counter)
            !           ^ TODO) to extend to Vup != Vdw: 1->NSPIN
         enddo
         !
         dummy_lambda(ibath)%element=a(counter+1:counter+Nlambdas(ibath))
         counter=counter+Nlambdas(ibath)
      enddo
      !
      dDelta=zero
      counter=0
      !
      do ibath=1,Nbath
         H_reconstructed = nnn2lso_reshape(Hreplica_build(dummy_lambda(ibath)%element),Nlat,Nspin,Norb)
         do l=1,Ldelta
            Haux(:,:,l) = zeye(Nlat*Nspin*Norb)*(xi*Xdelta(l)+xmu) - H_reconstructed
            call inv(Haux(:,:,l))
            invH_knn(:,:,:,:,:,:,l) = lso2nnn_reshape(Haux(:,:,l),Nlat,Nspin,Norb)
         enddo
         !Derivate_Vp
         counter=counter+1
         do ispin=1,Nspin
            dDelta(:,:,ispin,ispin,:,:,:,counter)=2d0*dummy_Vbath(1,ibath)*invH_knn(:,:,ispin,ispin,:,:,:)
            !                                                     ^ TODO) to extend to Vup != Vdw: 1->ispin
         enddo
         !Derivate_lambda_p
         do k=1,Nlambdas(ibath)
            counter = counter + 1
            Hbasis_lso=nnn2lso_reshape(Hreplica_basis(k)%O,Nlat,Nspin,Norb)
            do l=1,Ldelta
               Htmp = ((Haux(:,:,l) .x. Hbasis_lso)) .x. Haux(:,:,l)
               dDelta(:,:,:,:,:,:,l,counter)=lso2nnn_reshape((dummy_Vbath(1,ibath)**2)*Htmp,Nlat,Nspin,Norb)
            enddo
         enddo
      enddo
      !
      do ibath=1,Nbath
         if(allocated(dummy_lambda(ibath)%element))deallocate(dummy_lambda(ibath)%element)
      enddo
      !
   end function grad_delta_replica
   !
   !
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
            dG0and_lso = (G0and_lso .x. dDelta_lso) .x. G0and_lso
            dG0and(:,:,:,:,:,:,l,ik)=lso2nnn_reshape(dG0and_lso,Nlat,Nspin,Norb)
         enddo
      enddo
      !
   end function grad_g0and_replica

#endif

end MODULE ED_FIT_CHI2
