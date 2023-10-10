MODULE ED_BATH_FUNCTIONS
   USE SF_CONSTANTS, only: zero
   USE SF_IOTOOLS, only:free_unit,reg,file_length,str
   USE SF_LINALG, only: eye,inv,zeye,inv_her
   USE ED_INPUT_VARS
   USE ED_VARS_GLOBAL
   USE ED_BATH
   USE ED_AUX_FUNX
   implicit none

   private

   !##################################################################
   !
   !     DELTA FUNCTIONS
   !     G0 FUNCTIONS
   !     G0^{-1} FUNCTIONS
   !
   !##################################################################
   interface delta_bath
      module procedure :: delta_bath_single
      module procedure :: delta_bath_array
   end interface delta_bath

   interface invg0_bath
      module procedure :: invg0_bath_single
      module procedure :: invg0_bath_array
   end interface invg0_bath

   public :: delta_bath
   public :: g0and_bath
   public :: invg0_bath



contains


   function delta_bath_single(x) result(Delta)
      complex(8),intent(in)                                         :: x
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: Delta
      integer                                                       :: i,ih,L
      integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin,ibath
      integer                                                       :: io,jo
      !
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: Vk
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: invH_k
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: invH_knn
      !
      !
      Delta=zero
      !
      !
      invH_k=zero
      do ibath=1,Nbath
         Vk = dzdiag(dmft_bath%item(ibath)%v(:))
         invH_knn=Hbath_build(dmft_Bath%item(ibath)%lambda)
         invH_k = nnn2lso_reshape(invH_knn,Nlat,Nspin,Norb)
         invH_k = zeye(Nlat*Nspin*Norb)*x - invH_k
         call inv(invH_k)
         invH_k = matmul(matmul(Vk,invH_k),Vk)
         invH_knn = lso2nnn_reshape(invH_k,Nlat,Nspin,Norb)
         Delta = Delta + invH_knn
      enddo
      !
   end function delta_bath_single


   function delta_bath_array(x) result(Delta)
      complex(8),dimension(:),intent(in)                            :: x
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
      integer                                                       :: i,ih,L
      integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin,ibath
      integer                                                       :: io,jo
      !
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: Vk
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: invH_k
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: invH_knn
      !
      !
      Delta=zero
      !
      L = size(x)
      !
      invH_k=zero
      do i=1,L
         do ibath=1,Nbath
            Vk = dzdiag(dmft_bath%item(ibath)%v(:))
            invH_knn=Hbath_build(dmft_Bath%item(ibath)%lambda)
            invH_k = nnn2lso_reshape(invH_knn,Nlat,Nspin,Norb)
            invH_k = zeye(Nlat*Nspin*Norb)*x(i) - invH_k
            call inv(invH_k)
            invH_k = matmul(matmul(Vk,invH_k),Vk)
            invH_knn = lso2nnn_reshape(invH_k,Nlat,Nspin,Norb)
            Delta(:,:,:,:,:,:,i)=Delta(:,:,:,:,:,:,i) + invH_knn
         enddo
      enddo
      !
   end function delta_bath_array


   function g0and_bath(x) result(G0and)
      complex(8),dimension(:),intent(in)                            :: x
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
      integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nso,i,L
      real(8),dimension(size(x))                                    :: det
      complex(8),dimension(size(x))                                 :: fg,ff
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: fgorb,zeta
      !
      G0and = zero
      !
      L=size(x)
      !
      g0and=invg0_bath(x)
      do i=1,L
         fgorb=nnn2lso_reshape(g0and(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
         call inv(fgorb)
         g0and(:,:,:,:,:,:,i)=lso2nnn_reshape(fgorb,Nlat,Nspin,Norb)
      enddo
      !
   end function g0and_bath



   function invg0_bath_single(x) result(G0and)
      complex(8),intent(in)                                         :: x
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: G0and,Delta
      integer                                                       :: i,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nso,L
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: zeta
      !
      G0and = zero
      !
      !
      Delta = delta_bath(x)
      zeta = (x+xmu)*eye(Nlat*Nspin*Norb)
      G0and = lso2nnn_reshape(zeta,Nlat,Nspin,Norb)-impHloc-Delta
   end function invg0_bath_single


   function invg0_bath_array(x) result(G0and)
      complex(8),dimension(:),intent(in)                            :: x
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
      integer                                                       :: i,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nso,L
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: zeta
      !
      G0and = zero
      !
      L=size(x)
      !
      Delta = delta_bath(x)
      do i=1,L
         zeta = (x(i)+xmu)*eye(Nlat*Nspin*Norb)
         G0and(:,:,:,:,:,:,i) = lso2nnn_reshape(zeta,Nlat,Nspin,Norb)-impHloc-Delta(:,:,:,:,:,:,i)
      enddo
   end function invg0_bath_array

   ! Auxiliary ℝ -> ℂ vector-to-diagonal-matrix constructor
   function dzdiag(x) result(A)
      real(8),dimension(:)                   :: x
      complex(8),dimension(:,:),allocatable  :: A
      integer                                :: N,i
      N=size(x,1)
      allocate(A(N,N))
      A=0.d0
      do i=1,N
         A(i,i)=x(i)
      enddo
   end function dzdiag

END MODULE ED_BATH_FUNCTIONS
