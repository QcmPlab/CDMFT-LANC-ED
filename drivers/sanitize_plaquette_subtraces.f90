program sanitize_plaquette_subtraces
  USE CDMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  !
  implicit none
  integer                                                                :: Nimp,Nx,Ny,ilat,jlat
  !Density matrices:
  complex(8),allocatable,dimension(:,:)                                  :: reduced_density_matrix,cdm
  logical,allocatable,dimension(:,:)                                     :: orbital_mask
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer                                                                :: comm
  logical                                                                :: master

  !Init MPI: use of MPI overloaded functions in SciFor
  call init_MPI(comm,.true.)
  master = get_Master_MPI(comm)

  Nx = 2
  Ny = 2
  Nlat = Nx*Ny ! Standard 2x2 plaquette
  Norb = 1 ! Single band Hubbard model
  Nimp = Nlat*Norb ! needed by the subtrace code
  

  call naming_convention()

  !Load from file the cluster density matrix (cdm)
  call ed_read_dm("reduced_density_matrix_4sites.dat",4**(Nlat*Norb),cdm)

  !Retrieve ALL TWO-SITE REDUCED DENSITY MATRICES
  if(.not.allocated(orbital_mask))allocate(orbital_mask(Nlat,Norb))
  !All independent two-site RDMs (to check NNs and NNNs are equal)
  do ilat=1,Nlat-1
    orbital_mask = .false.
    orbital_mask(ilat,:) = .true.
    do jlat=ilat+1,Nlat
        orbital_mask(ilat+1:Nlat,:) = .false.
        orbital_mask(jlat,:) = .true.
        call ed_get_reduced_dm(reduced_density_matrix,orbital_mask,doprint=.true.,source_dm=cdm)
    enddo
  enddo


contains

  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Conventional indices for the cluster sites:
  !             y ^
  !             3 |    007 008 009
  !             2 |    004 005 006
  !             1 |    001 002 003
  !             0 |___ ___ ___ ___ ___ >
  !                 0   1   2   3   4  x
  !-------------------------------------------------------------------------------------------
  function indices2N(indices) result(N)
    integer,dimension(2)         :: indices
    integer                      :: N,i
    !
    N=Nx*(indices(2)-1)+indices(1)
    !
  end function indices2N
  !
  !
  function N2indices(N) result(indices)
    integer,dimension(2)         :: indices
    integer                      :: N,i
    !
    indices(1)=mod(N,Nx)
    if(indices(1)==0)then
       indices(1)=Nx
       indices(2)=(N-Nx)/Nx+1
    else
       indices(2)=N/Nx+1
    endif
    !
  end function N2indices


  !-------------------------------------------------------------------------------------------
  !PURPOSE: explicitate the cluster-indices convention in LOG files
  !-------------------------------------------------------------------------------------------
  subroutine naming_convention()
    integer                       :: i,j
    integer,dimension(Nx,Ny)      :: matrix
    !
    do j=1,Ny
       do i=1,Nx
          matrix(i,j)=indices2N([i,j])
       enddo
    enddo
    !
    write(LOGfile,"(A)")"The unique index of each site (on the cartesian plane) is as follows:"
    write(LOGfile,"(A)")" "
    do j=1,Ny
       write(LOGfile,"(20(I2,2x))")(matrix(i,Ny+1-j),i =1,Nx)
    enddo
    write(LOGfile,"(A)")" "
  end subroutine naming_convention

end program sanitize_plaquette_subtraces

