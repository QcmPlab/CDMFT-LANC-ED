! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private

  integer    :: iiup,iidw
  integer    :: iud,jj
  integer    :: i,iup,idw
  integer    :: j,jup,jdw
  integer    :: m,mup,mdw
  integer    :: ishift
  integer    :: isector,jsector
  integer    :: ms
  integer    :: impi
  integer    :: ilat,jlat,iorb,jorb,ispin,jspin,is,js,ibath
  integer    :: k1,k2,k3,k4
  integer    :: ialfa,ibeta
  real(8)    :: sg1,sg2,sg3,sg4
  complex(8) :: htmp,htmpup,htmpdw
  logical    :: Jcondition



  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
#ifdef _MPI
  public  :: directMatVec_MPI_main
#endif



contains


  subroutine directMatVec_main(Nloc,vin,Hv)
    integer                                                     :: Nloc
    complex(8),dimension(Nloc)                                  :: vin
    complex(8),dimension(Nloc)                                  :: Hv
    complex(8),dimension(:),allocatable                         :: vt,Hvt
    integer,dimension(Ns)                                       :: ibup,ibdw
    integer,dimension(2*Ns_Ud)                                  :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                             :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Nlat,Norb)                                :: Nup,Ndw
    real(8),dimension(Nlat,Nspin,Norb,Nbath)                    :: diag_hybr
    real(8),dimension(Nlat,Nspin,Norb,Nbath)                    :: bath_diag
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_reconstructed

    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    diag_hybr=zero
    bath_diag=zero
    do ibath=1,Nbath
      Hbath_reconstructed(:,:,:,:,:,:,ibath)=bath_from_sym(dmft_bath%item(ibath)%lambda)
      do ilat=1,Nlat
        do ispin=1,Nspin
          do iorb=1,Norb
            diag_hybr(ilat,ispin,iorb,ibath)=dmft_bath%item(ibath)%v
            bath_diag(ilat,ispin,iorb,ibath)=DREAL(Hbath_Reconstructed(ilat,ilat,ispin,ispin,iorb,iorb,ibath))
          enddo
        enddo
      enddo
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "ED_HAMILTONIAN/direct/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS
    include "ED_HAMILTONIAN/direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "ED_HAMILTONIAN/direct/HxV_dw.f90"
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       include "ED_HAMILTONIAN/direct/HxV_non_local.f90"
    endif
    !-----------------------------------------------!
    !
    return
  end subroutine directMatVec_main


#ifdef _MPI
  subroutine directMatVec_MPI_main(Nloc,vin,Hv)
    integer                                                     :: Nloc,N
    complex(8),dimension(Nloc)                                  :: Vin
    complex(8),dimension(Nloc)                                  :: Hv
    complex(8),dimension(:),allocatable                         :: vt,Hvt
    integer,dimension(Ns)                                       :: ibup,ibdw
    integer,dimension(2*Ns_Ud)                                  :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                             :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Nlat,Norb)                                :: Nup,Ndw
    real(8),dimension(Nlat,Nspin,Norb,Nbath)                    :: diag_hybr
    real(8),dimension(Nlat,Nspin,Norb,Nbath)                    :: bath_diag
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_reconstructed

    !
    integer                                  :: MpiIerr
    integer,allocatable,dimension(:)         :: Counts
    integer,allocatable,dimension(:)         :: Offset
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    ! if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    !Get diagonal hybridization, bath energy
    diag_hybr=zero
    bath_diag=zero
    do ibath=1,Nbath
      Hbath_reconstructed(:,:,:,:,:,:,ibath)=bath_from_sym(dmft_bath%item(ibath)%lambda)
      do ilat=1,Nlat
        do ispin=1,Nspin
          do iorb=1,Norb
            diag_hybr(ilat,ispin,iorb,ibath)=dmft_bath%item(ibath)%v
            bath_diag(ilat,ispin,iorb,ibath)=Hbath_Reconstructed(ilat,ilat,ispin,ispin,iorb,iorb,ibath)
          enddo
        enddo
      enddo
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "ED_HAMILTONIAN/direct_mpi/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "ED_HAMILTONIAN/direct_mpi/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    allocate(vt(mpiQup*DimDw)) ;vt=zero
    allocate(Hvt(mpiQup*DimDw));Hvt=zero
    call vector_transpose_MPI(DimUp,MpiQdw,Vin,DimDw,MpiQup,vt) !Vin^T --> Vt
    include "ED_HAMILTONIAN/direct_mpi/HxV_dw.f90"
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=zero        !reallocate Vt
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
    Hv = Hv + Vt
    deallocate(vt)
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       !
       allocate(vt(N)) ; vt = zero
       call allgather_vector_MPI(MpiComm,vin,vt)
       !
       include "ED_HAMILTONIAN/direct_mpi/HxV_non_local.f90"
       !
       deallocate(Vt)
    endif
    !-----------------------------------------------!
    !
    return
  end subroutine directMatVec_MPI_main

#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV
