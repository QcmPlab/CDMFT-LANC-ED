!-------------------------------------------------------------------!
! PURPOSE: INITIALIZE INTERNAL Hreplica STRUCTURES
!-------------------------------------------------------------------!

   subroutine allocate_hreplica(N)
      integer              :: N,isym
      !
      allocate(hreplica_basis(N))
      allocate(Hreplica_lambda(N))
      do isym=1,N
         allocate(hreplica_basis(isym)%O(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
         hreplica_basis(isym)%O=zero
         Hreplica_lambda(isym)=0.d0
      enddo
   end subroutine allocate_hreplica


   subroutine deallocate_hreplica()
      integer              :: isym
      !
      do isym=1,size(hreplica_basis)
         deallocate(hreplica_basis(isym)%O)
      enddo
      deallocate(hreplica_basis)
   end subroutine deallocate_hreplica


!+------------------------------------------------------------------+
!PURPOSE  : Set Hreplica from user defined Hloc
!1: [Nspin,Nspin,Norb,Norb]
!2: [Nspin*Norb,Nspin*Norb]
!+------------------------------------------------------------------+
   subroutine init_hreplica_direct_nnn(Hloc)
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
      call allocate_hreplica(counter)
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
                              hreplica_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=one
                              hreplica_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=one
                              Hreplica_lambda(counter)=DREAL(Hloc(ilat,jlat,ispin,ispin,iorb,jorb))
                           endif
                           !
                           if(DIMAG(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)then
                              counter=counter+1
                              hreplica_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=dcmplx(0.d0,1.d0)
                              hreplica_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=dcmplx(0.d0,1.d0)
                              Hreplica_lambda(counter)=DIMAG(Hloc(ilat,jlat,ispin,ispin,iorb,jorb))
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine init_hreplica_direct_nnn



   subroutine init_hreplica_direct_lso(Hloc)
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
      call allocate_hreplica(counter)
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
                              hreplica_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=one
                              hreplica_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=one
                              Hreplica_lambda(counter)=DREAL(Hloc(io,jo))
                           endif
                           !
                           if(DIMAG(Hloc(io,jo)).ne.0.d0)then
                              counter=counter+1
                              hreplica_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=dcmplx(0.d0,1.d0)
                              hreplica_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=dcmplx(0.d0,1.d0)
                              Hreplica_lambda(counter)=DIMAG(Hloc(io,jo))
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine init_hreplica_direct_lso



subroutine init_Hreplica_symmetries_site(Hvec,lambdavec)
  complex(8),dimension(:,:,:,:,:,:,:) :: Hvec
  real(8),dimension(:)                :: lambdavec
  integer                             :: isym,N
  !
  N=size(lambdavec)
  call assert_shape(Hvec,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,N],"init_Hreplica_symmetries","Hvec")
  !
  call allocate_hreplica(N)
  !
  do isym=1,N
     Hreplica_lambda(isym)  = lambdavec(isym)
     Hreplica_basis(isym)%O = Hvec(:,:,:,:,:,:,isym)
  enddo
  !
  if(ed_verbose>2)call print_hloc(Hreplica_build(Hreplica_lambda))
end subroutine init_Hreplica_symmetries_site


!+-------------------------------------------------------------------+
!PURPOSE  : Reconstruct H_replica from symmetry vector+matrices
!+-------------------------------------------------------------------+


 function Hreplica_build(lambdavec) result (H)
    real(8),dimension(:)                                         :: lambdavec
    integer                                                      :: isym
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: H
    !
    if(size(lambdavec).ne.size(hreplica_basis)) STOP "H_from_sym: Wrong coefficient vector size"
    H=zero
    do isym=1,size(lambdavec)
       H=H+lambdavec(isym)*hreplica_basis(isym)%O
    enddo
 end function Hreplica_build



!+-------------------------------------------------------------------+
!PURPOSE  : Create bath mask
!+-------------------------------------------------------------------+

function Hreplica_mask(hloc,wdiag,uplo) result(Hmask)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hloc
    logical,optional                                      :: wdiag,uplo
    logical                                               :: wdiag_,uplo_
    logical,dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)    :: Hmask
    integer                                               :: ilat,jlat,iorb,jorb,ispin,io,jo
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
  end function Hreplica_mask






