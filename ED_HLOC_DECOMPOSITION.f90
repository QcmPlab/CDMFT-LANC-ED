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
   public:: allocate_h_repr
   public:: deallocate_h_repr

   contains

   !-------------------------------------------------------------------!
   ! PURPOSE: INITIALIZE INTERNAL HLOC STRUCTURES
   !-------------------------------------------------------------------!

   subroutine allocate_h_repr(H,N)
      integer              :: N,isym
      type(H_repr)         :: H
      !
      H%N_dec=N
      allocate(H%decomposition(N))
      do isym=1,N
         H%decomposition(isym)%lambda=0.d0
         allocate(H%decomposition(isym)%O(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
      enddo
   end subroutine allocate_h_repr

   subroutine deallocate_h_repr(H)
      integer              :: isym
      type(H_repr)         :: H
      !
      do isym=1,H%N_dec
         H%decomposition(isym)%lambda=0.d0
         deallocate(H%decomposition(isym)%O)
      enddo
      deallocate(H%decomposition)
   end subroutine deallocate_h_repr

   function H_from_sym(H_sym) result (H)
      type(H_repr)                                                 :: H_sym
      integer                                                      :: isym
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: H
      !
      H=zero
      do isym=1,H_sym%N_dec
         H=H+H_sym%decomposition(isym)%lambda*H_sym%decomposition(isym)%O
      enddo
   end function H_from_sym

   subroutine init_hloc_direct(Hloc)
      integer                                               :: ilat,jlat,ispin,iorb,jorb,counter,io,jo
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hloc
      logical(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hmask
      !
      !
      impHloc=lso2nnn_reshape(Hloc,Nlat,Nspin,Norb)
      Hmask=.false.
      where(abs(impHloc)>1d-6)Hmask=.true.
      !
      !REAL
      call allocate_h_repr(impHloc_sym,2*count(Hmask))
      !COMPLEX
      !call allocate_h_repr(impHloc_sym,2*count(Hmask))
      counter=0
      !
      do ispin=1,Nspin
         do ilat=1,Nlat
            do jlat=1,Nlat
               do iorb=1,Norb
                  do jorb=1,Norb
                     io=imp_state_index(ilat,iorb)
                     jo=imp_state_index(jlat,jorb)
                     if((Hmask(ilat,jlat,ispin,ispin,iorb,jorb)).and.(io.le.jo))then
                        counter=counter+1
                        !
                        impHloc_sym%decomposition(counter)%O=zero
                        !
                        impHloc_sym%decomposition(counter)%O(ilat,jlat,ispin,ispin,iorb,jorb)=one
                        impHloc_sym%decomposition(counter)%O(jlat,ilat,ispin,ispin,jorb,iorb)=one
                        !REAL
                        impHloc_sym%decomposition(counter)%lambda=impHloc(ilat,jlat,ispin,ispin,iorb,jorb)
                        !COMPLEX
                        !impHloc_sym%decomposition(counter)%lambda=DREAL(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))
                        !
                        !counter=counter+1
                        !impHloc_sym%decomposition(counter)%O(ilat,jlat,ispin,ispin,iorb,jorb)=xi
                        !impHloc_sym%decomposition(counter)%O(jlat,ilat,ispin,ispin,jorb,iorb)=-xi
                        !impHloc_sym%decomposition(counter)%lambda=DIMAG(impHloc(ilat,jlat,ispin,ispin,iorb,jorb))
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine init_hloc_direct


   subroutine init_hloc_symmetries(Hvec,lambdavec,N)
      integer                                                       :: isym,N
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,N)       :: Hvec
      real(8),dimension(N)                                          :: lambdavec
      !
      call allocate_h_repr(impHloc_sym,N)
      do isym=1,N
         impHloc_sym%decomposition(isym)%lambda= lambdavec(isym)
         impHloc_sym%decomposition(isym)%O     = Hvec(:,:,:,:,:,:,isym)
      enddo
      !
      impHloc=H_from_sym(impHloc_sym)
      !
   end subroutine init_hloc_symmetries





END MODULE ED_HLOC_DECOMPOSITION
