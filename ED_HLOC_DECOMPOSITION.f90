MODULE ED_HLOC_DECOMPOSITION
   USE ED_VARS_GLOBAL
   USE ED_AUX_FUNX
   implicit none
   private

   interface set_Hloc
     module procedure init_Hloc_direct
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
         H_basis(isym)%O=0.d0
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

   subroutine init_hloc_direct(Hloc)
      integer                                               :: ilat,jlat,ispin,iorb,jorb,counter,io,jo,Nsym
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hloc
      logical(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hmask
      !
      !
      impHloc=lso2nnn_reshape(Hloc,Nlat,Nspin,Norb)
      Hmask=.false.
      do ispin=1,Nspin
         do ilat=1,Nlat
            do jlat=1,Nlat
               do iorb=1,Norb
                  do jorb=1,Norb
                     io=imp_state_index(ilat,iorb)
                     jo=imp_state_index(jlat,jorb)
                     if((impHloc(ilat,jlat,ispin,ispin,iorb,jorb).ne.zero).and.(io.le.jo))then
                        counter=counter+1
                        !COMPLEX
                        !if(io.ne.jo)counter=counter+1
                     endif
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
         do ilat=1,Nlat
            do jlat=1,Nlat
               do iorb=1,Norb
                  do jorb=1,Norb
                     io=imp_state_index(ilat,iorb)
                     jo=imp_state_index(jlat,jorb)
                     if((impHloc(ilat,jlat,ispin,ispin,iorb,jorb).ne.zero).and.(io.le.jo))then
                        counter=counter+1
                        H_basis(counter)%O(ilat,jlat,ispin,ispin,iorb,jorb)=one
                        H_basis(counter)%O(jlat,ilat,ispin,ispin,jorb,iorb)=one
                        !REAL
                        lambda_impHloc(counter)=impHloc(ilat,jlat,ispin,ispin,iorb,jorb)
                        !COMPLEX
                        !lambda_impHloc(counter)=DREAL(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))
                        !
                        !if(io.ne.jo)then
                        !counter=counter+1
                           !H_basis(counter)%O(ilat,jlat,ispin,ispin,iorb,jorb)=xi
                           !H_basis(counter)%O(jlat,ilat,ispin,ispin,jorb,iorb)=-xi
                           !lambda_impHloc(counter)=DIMAG(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))
                        !endif
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine init_hloc_direct


   subroutine init_hloc_symmetries(Hvec,lambdavec)
      integer                                   :: isym,N
      complex(8),dimension(:,:,:,:,:,:,:)       :: Hvec
      real(8),dimension(:)                      :: lambdavec
      !
      if(size(lambdavec).ne.size(H_basis)) STOP "Init_hloc: Wrong coefficient vector size"
      if(size(Hvec(1,1,1,1,1,1,:)).ne.size(H_basis)) STOP "Init_hloc: Wrong H_basis size"
      !
      N=size(lambdavec)
      !
      call allocate_h_basis(N)
      do isym=1,N
         lambda_impHloc(isym)= lambdavec(isym)
         H_basis(isym)%O     = Hvec(:,:,:,:,:,:,isym)
      enddo
      !
      impHloc=H_from_sym(lambda_impHloc)
      !
   end subroutine init_hloc_symmetries





END MODULE ED_HLOC_DECOMPOSITION
