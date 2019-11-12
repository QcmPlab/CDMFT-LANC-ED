MODULE ED_HLOC_DECOMPOSITION
   USE ED_VARS_GLOBAL
   USE ED_AUX_FUNX
   implicit none
   private

   interface set_Hloc
     module procedure init_Hloc_direct_lso
     module procedure init_Hloc_direct_nnn
     module procedure init_Hloc_symmetries
   end interface set_Hloc

   public:: set_Hloc
   public:: allocate_h_basis
   public:: deallocate_h_basis

   contains

   !-------------------------------------------------------------------!
   ! PURPOSE: INITIALIZE INTERNAL HLOC STRUCTURES
   !-------------------------------------------------------------------!

   !allocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient

   subroutine allocate_h_basis(N)
      integer              :: N,isym
      !
      allocate(H_basis(N))
      allocate(lambda_impHloc(N))
      do isym=1,N
         allocate(H_basis(isym)%O(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
         H_basis(isym)%O=zero
         lambda_impHloc(isym)=0.d0
      enddo
   end subroutine allocate_h_basis


   !deallocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient

   subroutine deallocate_h_basis()
      integer              :: isym
      !
      do isym=1,size(H_basis)
         deallocate(H_basis(isym)%O)
      enddo
      deallocate(H_basis)
   end subroutine deallocate_h_basis

   !reconstruct [Nlat][][Nspin][][Norb][] hamiltonian from basis expansion given [lambda]

   function H_from_sym(lambdavec) result (H)
      real(8),dimension(:)                                         :: lambdavec
      integer                                                      :: isym
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: H
      !
      if(size(lambdavec).ne.size(H_basis)) STOP "H_from_sym: Wrong coefficient vector size"
      H=zero
      do isym=1,size(lambdavec)
         H=H+lambdavec(isym)*H_basis(isym)%O
      enddo
   end function H_from_sym

   !initialize impHloc and the set [H_basis,lambda_impHloc]

   subroutine init_hloc_direct_lso(Hloc)
      integer                                               :: ilat,jlat,ispin,jspin,iorb,jorb,counter,io,jo,Nsym
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hloc
      logical(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hmask
      !
      allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
      impHloc=zero
      !
      impHloc=lso2nnn_reshape(Hloc,Nlat,Nspin,Norb)
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
                        if((impHloc(ilat,jlat,ispin,jspin,iorb,jorb).ne.zero).and.(io.le.jo))then
                           if(DREAL(impHloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)counter=counter+1
                           if(DIMAG(impHloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)counter=counter+1
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      call allocate_h_basis(counter)
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
                        if((impHloc(ilat,jlat,ispin,jspin,iorb,jorb).ne.zero).and.(io.le.jo))then
                           if(DREAL(impHloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)then
                              counter=counter+1
                              H_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=one
                              H_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=one
                              lambda_impHloc(counter)=DREAL(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))
                           endif
                           !
                           if(DIMAG(impHloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)then
                              counter=counter+1
                              H_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=dcmplx(0.d0,1.d0)
                              H_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=dcmplx(0.d0,1.d0)
                              lambda_impHloc(counter)=DIMAG(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine init_hloc_direct_lso


   subroutine init_hloc_direct_nnn(Hloc)
      integer                                               :: ilat,jlat,ispin,jspin,iorb,jorb,counter,io,jo,Nsym
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hloc
      logical(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hmask
      !
      allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
      impHloc=zero
      !
      impHloc=Hloc
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
                        if((impHloc(ilat,jlat,ispin,jspin,iorb,jorb).ne.zero).and.(io.le.jo))then
                           if(DREAL(impHloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)counter=counter+1
                           if(DIMAG(impHloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)counter=counter+1
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      call allocate_h_basis(counter)
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
                        if((impHloc(ilat,jlat,ispin,jspin,iorb,jorb).ne.zero).and.(io.le.jo))then
                           if(DREAL(impHloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)then
                              counter=counter+1
                              H_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=one
                              H_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=one
                              lambda_impHloc(counter)=DREAL(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))
                           endif
                           !
                           if(DIMAG(impHloc(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)then
                              counter=counter+1
                              H_basis(counter)%O(ilat,jlat,ispin,jspin,iorb,jorb)=dcmplx(0.d0,1.d0)
                              H_basis(counter)%O(jlat,ilat,ispin,jspin,jorb,iorb)=dcmplx(0.d0,1.d0)
                              lambda_impHloc(counter)=DIMAG(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine init_hloc_direct_nnn
   

   subroutine init_hloc_symmetries(Hvec,lambdavec)
      integer                                   :: isym,N
      complex(8),dimension(:,:,:,:,:,:,:)       :: Hvec
      real(8),dimension(:)                      :: lambdavec
      !
      if(size(lambdavec).ne.size(Hvec,7)) STOP "Init_hloc: Wrong coefficient vector size"
      !
      allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
      impHloc=zero
      !
      N=size(lambdavec)
      !
      call allocate_h_basis(N)
      !
      do isym=1,N
         lambda_impHloc(isym)= lambdavec(isym)
         H_basis(isym)%O     = Hvec(:,:,:,:,:,:,isym)
      enddo
      !
      impHloc=H_from_sym(lambda_impHloc)
      !
   end subroutine init_hloc_symmetries





END MODULE ED_HLOC_DECOMPOSITION
