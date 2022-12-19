!-------------------------------------------------------------------!
! PURPOSE: INITIALIZE INTERNAL Hbath STRUCTURES
!-------------------------------------------------------------------!

   subroutine allocate_Hbath(Nsym)
      integer              :: Nsym,isym
      !
      allocate(Hbath_basis(Nsym))
      allocate(Hbath_lambda(Nbath,Nsym))
      do isym=1,Nsym
         allocate(Hbath_basis(isym)%O(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
         Hbath_basis(isym)%O=zero
         Hbath_lambda(:,isym)=zero
      enddo
   end subroutine allocate_Hbath


   subroutine deallocate_Hbath()
      integer              :: isym
      !
      do isym=1,size(Hbath_basis)
         deallocate(Hbath_basis(isym)%O)
      enddo
      deallocate(Hbath_basis)
      deallocate(Hbath_lambda)
   end subroutine deallocate_Hbath


!+------------------------------------------------------------------+
!PURPOSE  : Set Hbath from user defined Hloc
!1: [Nspin,Nspin,Norb,Norb]
!2: [Nspin*Norb,Nspin*Norb]
!+------------------------------------------------------------------+
   subroutine init_Hbath_direct_nnn(Hloc)
      integer                                               :: ilat,jlat,ispin,jspin,iorb,jorb,counter,io,jo,Nsym
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hloc
      logical(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hmask
      !
      counter=0
      !
      Hmask=.false.
      !SPIN DIAGONAL
      do ispin=1,Nspin
         do jspin=1,Nspin
            do ilat=1,Nlat
               do jlat=1,Nlat
                  do iorb=1,Norb
                     do jorb=1,Norb
                        io=index_stride_lso(ilat,ispin,iorb)
                        jo=index_stride_lso(jlat,jspin,jorb)
                        if((Hloc(ilat,jlat,ispin,jspin,iorb,jorb).ne.zero).and.(io.le.jo))then
                           if(DREAL(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)counter=counter+1
                           if(DIMAG(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)counter=counter+1
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      call allocate_Hbath(counter)
      !
      counter=0
      !
      do ispin=1,Nspin
         do jspin=1,Nspin
            do ilat=1,Nlat
               do jlat=1,Nlat
                  do iorb=1,Norb
                     do jorb=1,Norb
                        io=index_stride_lso(ilat,ispin,iorb)
                        jo=index_stride_lso(jlat,jspin,jorb)
                        if((Hloc(ilat,jlat,ispin,jspin,iorb,jorb).ne.zero).and.(io.le.jo))then
                           if(DREAL(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)then
                              counter=counter+1
                              Hbath_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=one
                              Hbath_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=one
                              Hbath_lambda(:,counter)=DREAL(Hloc(ilat,jlat,ispin,ispin,iorb,jorb))
                           endif
                           !
                           if(DIMAG(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)then
                              counter=counter+1
                              Hbath_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=dcmplx(0.d0,1.d0)
                              Hbath_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=dcmplx(0.d0,1.d0)
                              Hbath_lambda(:,counter)=DIMAG(Hloc(ilat,jlat,ispin,ispin,iorb,jorb))
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine init_Hbath_direct_nnn



   subroutine init_Hbath_direct_lso(Hloc)
      integer                                               :: ilat,jlat,ispin,jspin,iorb,jorb,counter,io,jo,Nsym
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hloc
      logical(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hmask
      !
      !
      Hmask=.false.
      !SPIN DIAGONAL
      do ispin=1,Nspin
         do jspin=1,Nspin
            do ilat=1,Nlat
               do jlat=1,Nlat
                  do iorb=1,Norb
                     do jorb=1,Norb
                        io=index_stride_lso(ilat,ispin,iorb)
                        jo=index_stride_lso(jlat,jspin,jorb)
                        if((Hloc(io,jo).ne.zero).and.(io.le.jo))then
                           if(DREAL(Hloc(io,jo)).ne.0.d0)counter=counter+1
                           if(DIMAG(Hloc(io,jo)).ne.0.d0)counter=counter+1
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      call allocate_Hbath(counter)
      !
      counter=0
      !
      do ispin=1,Nspin
         do jspin=1,Nspin
            do ilat=1,Nlat
               do jlat=1,Nlat
                  do iorb=1,Norb
                     do jorb=1,Norb
                        io=index_stride_lso(ilat,ispin,iorb)
                        jo=index_stride_lso(jlat,jspin,jorb)
                        if((Hloc(io,jo).ne.zero).and.(io.le.jo))then
                           if(DREAL(Hloc(io,jo)).ne.0.d0)then
                              counter=counter+1
                              Hbath_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=one
                              Hbath_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=one
                              Hbath_lambda(:,counter)=DREAL(Hloc(io,jo))
                           endif
                           !
                           if(DIMAG(Hloc(io,jo)).ne.0.d0)then
                              counter=counter+1
                              Hbath_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=dcmplx(0.d0,1.d0)
                              Hbath_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=dcmplx(0.d0,1.d0)
                              Hbath_lambda(:,counter)=DIMAG(Hloc(io,jo))
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine init_Hbath_direct_lso



subroutine init_Hbath_symmetries_site(Hvec,lambdavecs)
  complex(8),dimension(:,:,:,:,:,:,:) :: Hvec       ![size(Hloc),Nsym]
  real(8),dimension(:,:)              :: lambdavecs ![Nbath,Nsym]
  integer                             :: isym,Nsym,ibath
  !
  if(size(lambdavecs(:,1))/=Nbath)then
      write(*,*) "                                                                               "
      write(*,*) "ERROR: if you are trying to init Hbath for inequivalent clusters please note"
      write(*,*) "       that the lambdasym array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
      write(*,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
      write(*,*) "       the new recommended [Nbath]x[Nsym] shape for the single cluster case.   "
      write(*,*) "                                                                               "
      stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
  else
      Nsym=size(lambdavecs(Nbath,:))
  endif
  !
  call assert_shape(Hvec,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nsym],"init_Hbath_symmetries","Hvec")
  !
  call allocate_Hbath(Nsym)
  !
  do isym=1,Nsym
     Hbath_lambda(:,isym) = lambdavecs(:,isym)
     Hbath_basis(isym)%O  = Hvec(:,:,:,:,:,:,isym)
  enddo
  !
  if(ed_verbose>2)then
   do ibath=1,Nbath
      write(*,*) "Hbath #"//str(ibath)//":"
      call print_hloc(Hbath_build(Hbath_lambda(ibath,:)))
   enddo
  endif
  !
end subroutine init_Hbath_symmetries_site

subroutine init_Hbath_symmetries_LEGACY(Hvec,lambdavec)
   complex(8),dimension(:,:,:,:,:,:,:) :: Hvec      ![size(Hloc),Nsym]
   real(8),dimension(:)                :: lambdavec ![Nsym]
   integer                             :: isym,Nsym,ibath
   !
   Nsym=size(lambdavec)
   call assert_shape(Hvec,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nsym],"init_Hbath_symmetries","Hvec")
   !
   call allocate_Hbath(Nsym)
   !
   do isym=1,Nsym
      do ibath=1,Nbath
      !> BACK-COMPATIBILITY PATCH (cfr. init_dmft_bath in dmft_aux.f90)
         Hbath_lambda(ibath,isym) = lambdavec(isym)
      enddo
      Hbath_basis(isym)%O = Hvec(:,:,:,:,:,:,isym)
   enddo
   !
   ! PRINT DEPRECATION MESSAGE TO LOG
   write(*,*) "                                                                               "
   write(*,*) "WARNING: Passing a single lambdasym vector to ed_set_Hbath is /deprecated/. "
   write(*,*) "         You should instead define a different lambda for each bath component, "
   write(*,*) "         namely passing a [Nbath]x[Nsym] array instead of a [Nsym] vector.     "
   write(*,*) "         Your single lambda vector has been internally copied into the required"
   write(*,*) "         higher-rank array, so giving each replica the same set of lambdas.    "
   write(*,*) "         >>> This back-compatibility patch might be removed in a future update."
   write(*,*) "                                                                               "
   !
   if(ed_verbose>2)then
    do ibath=1,Nbath
       write(*,*) "Hbath #"//str(ibath)//":"
       call print_hloc(Hbath_build(Hbath_lambda(ibath,:)))
    enddo
   endif
   !
 end subroutine init_Hbath_symmetries_LEGACY

!+-------------------------------------------------------------------+
!PURPOSE  : Reconstruct H_replica from symmetry vector+matrices
!+-------------------------------------------------------------------+


 function Hbath_build(lambdavec) result (H)
    real(8),dimension(:)                                         :: lambdavec ![Nsym]
    integer                                                      :: isym
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: H
    !
    if(size(lambdavec).ne.size(Hbath_basis)) STOP "H_from_sym: Wrong coefficient vector size"
    H=zero
    do isym=1,size(lambdavec)
       H=H+lambdavec(isym)*Hbath_basis(isym)%O
    enddo
 end function Hbath_build



!+-------------------------------------------------------------------+
!PURPOSE  : Create bath mask
!+-------------------------------------------------------------------+

function Hbath_mask(wdiag,uplo) result(Hmask)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hloc
    logical,optional                                      :: wdiag,uplo
    logical                                               :: wdiag_,uplo_
    logical,dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)    :: Hmask
    integer                                               :: ilat,jlat,iorb,jorb,ispin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    Hloc = Hbath_build(Hbath_lambda(Nbath,:)) !The mask should be replica-independent
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
  end function Hbath_mask






