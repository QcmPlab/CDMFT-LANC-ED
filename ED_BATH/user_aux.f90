!+-------------------------------------------------------------------+
!PURPOSE  : Inquire the correct bath size to allocate the 
! the bath array in the calling program.
!+-------------------------------------------------------------------+
function get_bath_dimension_direct(Hloc_nn) result(bath_size)
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb),intent(in)          :: Hloc_nn
  integer                                                                   :: bath_size,ndx,ilat,jlat,ispin,jspin,iorb,jorb,io,jo,counter
  !
  counter=0
  !
  !Real part of nonzero elements
  do ispin=1,Nspin
     do jspin=1,Nspin
        do ilat=1,Nlat
           do jlat=1,Nlat
              do iorb=1,Norb
                 do jorb=1,Norb
                    io=index_stride_lso(ilat,ispin,iorb)
                    jo=index_stride_lso(jlat,jspin,jorb)
                    if((Hloc_nn(ilat,jlat,ispin,jspin,iorb,jorb).ne.zero).and.(io.le.jo))then
                       if(DREAL(Hloc_nn(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)counter=counter+1
                       if(DIMAG(Hloc_nn(ilat,jlat,ispin,jspin,iorb,jorb)).ne.0.d0)counter=counter+1
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  ndx   = counter
  ndx   = ndx + 1             !we also print n_Dec
  !
  !number of non vanishing elements for each replica
  ndx = ndx * Nbath
  !diagonal hybridizations: Vs
  ndx = ndx + Nbath
  !
  bath_size = ndx
  !
end function get_bath_dimension_direct


function get_bath_dimension_symmetries(Hloc_nn) result(bath_size)
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)                     :: Hloc_nn
  integer                                                            :: bath_size,ndx,isym,Nsym
  !
  !number of symmetries
  Nsym=size(Hloc_nn(1,1,1,1,1,1,:))
  !
  ndx=Nsym
  !
  !for each replica we also print N_dec
  ndx=ndx+1
  !
  !number of replicas
  ndx = ndx * Nbath
  !diagonal hybridizations: Vs
  ndx = ndx + Nbath
  !
  bath_size = ndx
  !
end function get_bath_dimension_symmetries

!+-------------------------------------------------------------------+
!PURPOSE  : Check if the dimension of the bath array are consistent
!+-------------------------------------------------------------------+
function check_bath_dimension(bath_) result(bool)
  real(8),dimension(:)           :: bath_
  integer                        :: Ntrue,i
  logical                        :: bool
  complex(8),allocatable         :: Hbath(:,:,:,:,:,:,:)![Nlat][:][Nspin][:][Norb][:][Nsym]
  !
  if(.not.allocated(Hbath_basis))STOP "check_bath_dimension: Hbasis not allocated"
  !
  if(.not.allocated(Hbath))allocate(Hbath(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Hbath_basis)))
  !
  do i=1,size(Hbath_basis)
    Hbath(:,:,:,:,:,:,i)=Hbath_basis(i)%O
  enddo
  !
  Ntrue = get_bath_dimension(Hbath)
  bool  = ( size(bath_) == Ntrue )
end function check_bath_dimension


!##################################################################
!
!     USER BATH PREDEFINED SYMMETRIES:
!
!##################################################################

!+-------------------------------------------------------------------+
!PURPOSE  : given a bath array apply a specific transformation or 
! impose a given symmetry:
! - break spin symmetry by applying a symmetry breaking field
! - given a bath array set both spin components to have 
!    the same bath, i.e. impose non-magnetic solution
! - given a bath array enforces the particle-hole symmetry 
!    by setting the positive energies in modulo identical to the negative
!    ones.
!+-------------------------------------------------------------------+
subroutine impose_equal_lambda(bath_,ibath,lambdaindex_vec)
  real(8),dimension(:)    :: bath_
  real(8)                 :: val
  integer,dimension(:)    :: lambdaindex_vec
  integer                 :: i,N,ibath
  !
  call allocate_dmft_bath()
  call set_dmft_bath(bath_)
  !
  N=size(lambdaindex_vec)
  val=0.d0
  do i=1,N
    val=val+dmft_bath%item(ibath)%lambda(lambdaindex_vec(i))/N
  enddo
  !
  do i=1,N
    dmft_bath%item(ibath)%lambda(lambdaindex_vec(i))=val
  enddo
  !
  call get_dmft_bath(bath_)
  call deallocate_dmft_bath()
end subroutine impose_equal_lambda


subroutine impose_bath_offset(bath_,ibath,offset)
  real(8),dimension(:)    :: bath_
  real(8)                 :: offset
  integer                 :: isym,N,ibath
  !
  call allocate_dmft_bath()
  call set_dmft_bath(bath_)
  !
  if(size(Hbath_lambda) .ne. dmft_bath%item(ibath)%N_dec)then
    dmft_bath%item(ibath)%lambda(dmft_bath%item(ibath)%N_dec)=offset
  else
    do isym=1,size(Hbath_lambda)
      if(is_identity(Hbath_basis(isym)%O)) dmft_bath%item(ibath)%lambda(isym)=offset
      return
    enddo
  endif
  !
  !
  call get_dmft_bath(bath_)
  call deallocate_dmft_bath()
  !
end subroutine impose_bath_offset

