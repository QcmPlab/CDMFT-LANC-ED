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
  !     DELTA FUNCTIONS
  !     G0 FUNCTIONS
  !     G0^{-1} FUNCTIONS
  !
  !##################################################################

  public :: delta_bath
  public :: g0and_bath
  public :: invg0_bath



contains




  function delta_bath(x) result(Delta)
    complex(8),dimension(:),intent(in)                            :: x
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                                       :: i,ih,L
    integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin,ibath
    integer                                                       :: io,jo
    real(8),dimension(Nbath)                                      :: eps,dps,vps
    real(8),dimension(Norb,Nbath)                                 :: vops
    !
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
          invH_k = nnn2lso_reshape(dmft_bath%item(ibath)%h,Nlat,Nspin,Norb)
          invH_k = zeye(Nlat*Nspin*Norb)*x(i) - invH_k
          call inv(invH_k)
          Delta(:,:,:,:,:,:,i)=Delta(:,:,:,:,:,:,i) + &
               (dmft_bath%item(ibath)%v**2)*lso2nnn_reshape(invH_k,Nlat,Nspin,Norb) !remember module of V for complex type
       enddo
       !
    enddo
    !
  end function delta_bath


  function g0and_bath(x) result(G0and)
    complex(8),dimension(:),intent(in)                            :: x
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                                       :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    real(8),dimension(size(x))                                    :: det
    complex(8),dimension(size(x))                                 :: fg,ff
    complex(8),dimension(Nlat*Norb,Nlat*Norb)                     :: fgorb,zeta
    !
    G0and = zero
    !
    L=size(x)
    !
    !> workout a spin-block at the time
    Delta = delta_bath(x)
    do ispin=1,Nspin
       do i=1,L
          fgorb= zero
          zeta = (x(i)+xmu)*eye(Nlat*Norb)
          fgorb= zeta-nn2so_reshape(impHloc(:,:,ispin,ispin,:,:)-Delta(:,:,ispin,ispin,:,:,i),Nlat,Norb)
          call inv(fgorb)
          G0and(:,:,ispin,ispin,:,:,i)=so2nn_reshape(fgorb,Nlat,Norb)
       enddo
    enddo
    !
  end function g0and_bath



  function invg0_bath(x) result(G0and)
    complex(8),dimension(:),intent(in)                            :: x
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                                       :: i,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,Nso,L
    complex(8),dimension(Nlat*Norb,Nlat*Norb)                     :: zeta
    !
    G0and = zero
    !
    L=size(x)
    !
    Delta = delta_bath(x)
    do i=1,L
       zeta = (x(i)+xmu)*eye(Nlat**Nspin*Norb)
       G0and(:,:,:,:,:,:,i) = lso2nnn_reshape(zeta,Nlat,Nspin,Norb)-impHloc-Delta(:,:,:,:,:,:,i)
    enddo
  end function invg0_bath





END MODULE ED_BATH_FUNCTIONS
