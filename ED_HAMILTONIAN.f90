MODULE ED_HAMILTONIAN
   USE ED_HAMILTONIAN_COMMON
   USE ED_HAMILTONIAN_SPARSE_HxV
   USE ED_HAMILTONIAN_DIRECT_HxV
   USE SF_SP_LINALG, only: sp_lanc_tridiag
   !
   implicit none
   private


   !>Build sparse hamiltonian of the sector
   public  :: build_Hv_sector
   public  :: delete_Hv_sector
   public  :: vecDim_Hv_sector

   !> Tridiag sparse Hamiltonian of the sector
   public  :: tridiag_Hv_sector

   !>Sparse Mat-Vec product using stored sparse matrix
   public  :: spMatVec_main
#ifdef _MPI
   public  :: spMatVec_MPI_main
#endif


   !>Sparse Mat-Vec direct on-the-fly product
   public  :: directMatVec_main
#ifdef _MPI
   public  :: directMatVec_MPI_main
#endif




contains




   !####################################################################
   !                 MAIN ROUTINES: BUILD/DELETE SECTOR
   !####################################################################
   subroutine build_Hv_sector(isector,Hmat)
      integer                            :: isector,SectorDim
      complex(8),dimension(:,:),optional :: Hmat
      integer                            :: irank,ierr
      integer                            :: i,iup,idw
      integer                            :: j,jup,jdw
      !
      allocate(Hs(2*Ns_Ud))
      allocate(DimUps(Ns_Ud))
      allocate(DimDws(Ns_Ud))
      !LEGACY:
      !call build_sector(isector,Hs)
      !BECOMES:
      call build_sector(isector,global_sector)
      Hs = global_sector%H
      !
      !LEGACY:
      !call get_DimUp(isector,DimUps)
      !call get_DimDw(isector,DimDws)
      !DimUp = product(DimUps)
      !DimDw = product(DimDws)
      !Dim   = getDim(isector)
      !BECOMES:
      DimUp = global_sector%DimUp
      DimDw = global_sector%DimDw
      DimUps = global_sector%DimUps
      DimDws = global_sector%DimDws
      Dim = global_sector%Dim
      !
      mpiAllThreads=.true.
      !>PREAMBLE: check that split of the DW is performed with the minimum #cpu: no idle cpus allowed (with zero elements)
#ifdef _MPI
      if(MpiStatus)then
         if(DimDw < MpiSize)then
            if(MpiMaster.AND.ed_verbose>4)write(*,*)"Reducing N_cpu to DimDw:",DimDw,MpiSize-DimDw
            allocate(MpiMembers(0:DimDw-1))
            forall(irank=0:DimDw-1)MpiMembers(irank)=irank
            call Mpi_Group_Incl(MpiGroup_Global,DimDw,MpiMembers,MpiGroup,ierr)
            call Mpi_Comm_create(MpiComm_Global,MpiGroup,MpiComm,ierr)
            deallocate(MpiMembers)
            mpiAllThreads=.false.
            call Barrier_MPI(MpiComm_Global)
#ifdef _DEBUG
            if(ed_verbose>4)then
               if(MpiMaster)write(LOGfile,*)&
                  "       mpiRank,   MpiComm, Comm_Global, Comm_World, Comm_Null, Undefined"
               do i=0,MpiSize-1
                  call Barrier_MPI(MpiComm_Global)
                  if(MpiRank==i)write(*,*)i,MpiComm,MpiComm_Global,Mpi_Comm_World,Mpi_comm_null,Mpi_Undefined
               enddo
               call Barrier_MPI(MpiComm_Global)
            endif
#endif
         endif
         if( MpiComm /= MPI_COMM_NULL )then
            MpiRank = Get_Rank_MPI(MpiComm)
            MpiSize = Get_Size_MPI(MpiComm)
         endif
      endif
#endif
      !
      !Dw split:
      mpiQdw = DimDw/MpiSize
      mpiRdw = mod(DimDw,MpiSize)
      if(MpiRank < mod(DimDw,MpiSize) ) then
         mpiRdw = 0
         MpiQdw = MpiQdw+1
      endif
      !
      !Total split: split DW \times UP
      mpiQ = DimUp*mpiQdw
      mpiR = DimUp*mpiRdw
      mpiIstart = 1 + MpiRank*mpiQ+mpiR
      mpiIend   = (MpiRank+1)*mpiQ+mpiR
      mpiIshift = MpiRank*mpiQ+mpiR
      !
      !
#ifdef _MPI
#ifdef _DEBUG
      if(MpiStatus.AND.ed_verbose>4.AND.(MpiComm/=Mpi_Comm_Null).AND.MpiSize>=1)then
         if(MpiMaster)write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,      mpi_Qdw,      mpiR_dw,  mpi_Istart,  mpi_Iend,  Iend-Istart,  Comm, Comm_Global"
         do irank=0,MpiSize-1
            call Barrier_MPI(MpiComm)
            if(MpiRank==irank)write(*,*)MpiRank,MpiQ,MpiR,mpiQdw,MpiRdw,MpiIstart,MpiIend,MpiIend-MpiIstart+1,MpiComm,MpiComm_Global
         enddo
         call Barrier_MPI(MpiComm)
      endif
#endif
#endif
      !
      !
      if(present(Hmat))then
         spHtimesV_p => null()
         call ed_buildh_main(isector,Hmat)
         return
      endif
      !
      select case (ed_sparse_H)
       case (.true.)
         spHtimesV_p => spMatVec_main
#ifdef _MPI
         if(MpiStatus)spHtimesV_p => spMatVec_MPI_main
#endif
         call ed_buildh_main(isector)
       case (.false.)
         spHtimesV_p => directMatVec_main
#ifdef _MPI
         if(MpiStatus)spHtimesV_p => directMatVec_MPI_main
#endif
      end select
      !
   end subroutine build_Hv_sector





   subroutine delete_Hv_sector()
      integer :: iud,ierr,i
      call delete_sector(global_sector)
      deallocate(Hs)
      deallocate(DimUps)
      deallocate(DimDws)
      !
      !There is no difference here between Mpi and serial version, as Local part was removed.
#ifdef _MPI
      if(MpiStatus)then
         call sp_delete_matrix(MpiComm,spH0d)
         if(Jhflag)call sp_delete_matrix(MpiComm,spH0nd)
      else
         call sp_delete_matrix(spH0d)
         if(Jhflag)call sp_delete_matrix(spH0nd)
      endif
#else
      call sp_delete_matrix(spH0d)
      if(Jhflag)call sp_delete_matrix(spH0nd)
#endif
      do iud=1,Ns_Ud
         call sp_delete_matrix(spH0ups(iud))
         call sp_delete_matrix(spH0dws(iud))
      enddo
      !
      spHtimesV_p => null()
      !
#ifdef _MPI
      if(MpiStatus)then
         if(MpiGroup/=Mpi_Group_Null)call Mpi_Group_free(MpiGroup,ierr)
         if(MpiComm/=Mpi_Comm_Null.AND.MpiComm/=Mpi_Comm_World)call Mpi_Comm_Free(MpiComm,ierr)
         MpiComm = MpiComm_Global
         MpiSize = get_Size_MPI(MpiComm_Global)
         MpiRank = get_Rank_MPI(MpiComm_Global)
         !call Mpi_Comm_Dup(MpiComm_Global,MpiComm,ierr)
      endif
#endif
      iter=0
      !
   end subroutine delete_Hv_sector






   function vecDim_Hv_sector(isector) result(vecDim)
      integer :: isector
      integer :: vecDim
      integer :: mpiQdw
      integer :: DimUps(Ns_Ud),DimUp
      integer :: DimDws(Ns_Ud),DimDw
      !
      call get_DimUp(isector,DimUps) ; DimUp = product(DimUps)
      call get_DimDw(isector,DimDws) ; DimDw = product(DimDws)
      !
#ifdef _MPI
      if(MpiStatus)then
         !Dw split:
         mpiQdw = DimDw/MpiSize
         if(MpiRank < mod(DimDw,MpiSize) ) MpiQdw = MpiQdw+1
      else
         mpiQdw = DimDw
      endif
#else
      mpiQdw = DimDw
#endif
      !
      vecDim=DimUp*mpiQdw
      !
   end function vecDim_Hv_sector



   subroutine tridiag_Hv_sector(isector,vvinit,alanc,blanc,norm2)
      integer                            :: isector
      complex(8),dimension(:)            :: vvinit
      real(8),dimension(:),allocatable   :: alanc,blanc
      real(8)                            :: norm2
      !
      complex(8),dimension(:),allocatable:: vvloc
      integer                            :: vecDim
      !
      if(MpiMaster)then
         norm2=dot_product(vvinit,vvinit)
         vvinit=vvinit/sqrt(norm2)
      endif
#ifdef _MPI
      if(MpiStatus)call bcast_MPI(MpiComm,norm2)
#endif
      call build_Hv_sector(isector)
      allocate(alanc(global_sector%Nlanc),blanc(global_sector%Nlanc))
      alanc=0d0 ; blanc=0d0
      if(norm2/=0d0)then
#ifdef _MPI
         if(MpiStatus)then
            vecDim = vecDim_Hv_sector(isector)
            allocate(vvloc(vecDim))
            call scatter_vector_MPI(MpiComm,vvinit,vvloc)
            call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alanc,blanc)
         else
            call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
         endif
#else
         call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
#endif
      endif
      call delete_Hv_sector()
   end subroutine tridiag_Hv_sector


end MODULE ED_HAMILTONIAN
