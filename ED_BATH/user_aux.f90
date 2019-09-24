!+-------------------------------------------------------------------+
!PURPOSE  : Inquire the correct bath size to allocate the 
! the bath array in the calling program.
!+-------------------------------------------------------------------+
function get_bath_dimension(Hloc_nn) result(bath_size)
  complex(8),optional,intent(in) :: Hloc_nn(:,:,:,:,:,:)
  integer                        :: bath_size,ndx,ilat,jlat,ispin,iorb,jorb,io,jo
  real(8),allocatable            :: Hloc(:,:,:,:,:,:)
  !
  !off-diagonal non-vanishing elements
  if(present(Hloc_nn))then
     allocate(Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb));Hloc=dreal(Hloc_nn)
  elseif(allocated(impHloc))then
     allocate(Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb));Hloc=impHloc
  else
     stop "ERROR: get_bath_dimension: neither Hloc_nn present nor impHloc allocated"
  endif
  !
  ndx=0

  do ispin=1,Nspin
     do ilat=1,Nlat
        do jlat=1,Nlat
           do iorb=1,Norb
              do jorb=1,Norb
                 io = index_stride_lso(ilat,ispin,iorb)
                 jo = index_stride_lso(jlat,ispin,jorb)
                 if( (io<jo) .AND. (abs(Hloc(ilat,jlat,ispin,ispin,iorb,jorb))>1d-6) )ndx=ndx+1
              endif
           enddo
        enddo
     enddo
  enddo
  !Real diagonal elements (always assumed)
  ndx = ndx + Nspin*Nlat*Norb
  !number of non vanishing elements for each replica
  ndx = ndx * Nbath
  !diagonal hybridizations: Vs
  ndx = ndx + Nbath
  !
  bath_size = ndx
  !
end function get_bath_dimension





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
subroutine break_symmetry_bath_site(bath_,field,sign,save)
  real(8),dimension(:) :: bath_
  real(8)              :: field
  real(8)              :: sign
  logical,optional     :: save
  logical              :: save_
  save_=.true.;if(present(save))save_=save
  call allocate_dmft_bath()
  call set_dmft_bath(bath_)
  do ibath=1,Nbath
     do ilat=1,Nlat
        do iorb=1,Norb
           dmft_bath(ibath)%h(ilat,ilat,1,1,iorb,iorb)        = &
                dmft_bath(ibath)%h(ilat,ilat,1,1,iorb,iorb)         + sign*field
           dmft_bath(ibath)%h(ilat,ilat,Nspin,Nspin,iorb,iorb)= &
                dmft_bath(ibath)%h(ilat,ilat,Nspin,Nspin,iorb,iorb) + sign*field
        enddo
     enddo
  enddo
  if(save_)call save_dmft_bath()
  call get_dmft_bath(bath_)
  call deallocate_dmft_bath()
end subroutine break_symmetry_bath_site

!---------------------------------------------------------!

subroutine spin_symmetrize_bath_site(bath_,save)
  real(8),dimension(:)   :: bath_
  logical,optional       :: save
  logical                :: save_
  integer
  save_=.true.;if(present(save))save_=save
  if(Nspin==1)then
     write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
     return
  endif
  !
  call allocate_dmft_bath()
  call init_dmft_bathmask()
  call set_dmft_bath(bath_)

  do ibath=1,Nbath
     do ilat=1,Nlat
        do iorb=1,Norb
           dmft_bath(ibath)%h(ilat,ilat,Nspin,Nspin,iorb,iorb) = dmft_bath(ibath)%h(ilat,ilat,1,1,iorb,iorb)
           dmft_bath(ibath)%v(ilat,Nspin,iorb)                 = dmft_bath(ibath)%v(     ilat,  1,     iorb)
        enddo
     enddo
  enddo
  if(save_)call save_dmft_bath()
  call get_dmft_bath(bath_)
  call deallocate_dmft_bath()
end subroutine spin_symmetrize_bath_site

!---------------------------------------------------------!


subroutine orb_equality_bath_site(bath_,indx,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath
  integer,optional       :: indx
  logical,optional       :: save
  integer                :: indx_
  logical                :: save_
  indx_=1     ;if(present(indx))indx_=indx
  save_=.true.;if(present(save))save_=save
  if(Norb==1)then
     write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
     return
  endif
  !
  call allocate_dmft_bath()
  call init_dmft_bathmask()
  call set_dmft_bath(bath_)
  !
  do ibath=1,Nbath
     do iorb=1,Norb
        if(iorb==indx_)cycle
        dmft_bath(ibath)%h(:,:,:,:,iorb,iorb)=dmft_bath(ibath)%h(:,:,:,:,indx_,indx_)
        dmft_bath(ibath)%v(:,iorb,:)         =dmft_bath(ibath)%v(:,indx_,:)
     enddo
  enddo
  !
  if(save_)call save_dmft_bath()
  call get_dmft_bath(bath_)
  call deallocate_dmft_bath()
end subroutine orb_equality_bath_site



subroutine ph_symmetrize_bath_site(bath_,save)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath
  integer                :: ibath
  logical,optional       :: save
  logical                :: save_
  save_=.true.;if(present(save))save_=save
  call allocate_dmft_bath()
  call set_dmft_bath(bath_)
  if(Nbath==1)return
  !
  do ilat=1,Nlat
     do ispin=1,Nspin
        do iorb=1,Norb
           !
           if(mod(Nbath,2)==0)then
              do ibath=1,Nbath/2
                 dmft_bath(Nbath+1-ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)=-dmft_bath(ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)
                 dmft_bath(Nbath+1-ibath)%v(ilat,ispin,iorb)= dmft_bath(ibath)%v(ilat,ispin,iorb)
              enddo
           else
              do ibath=1,(Nbath-1)/2
                 dmft_bath(Nbath+1-ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)=-dmft_bath(ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)
                 dmft_bath(Nbath+1-ibath)%v(ilat,ispin,iorb)= dmft_bath(ibath)%v(ilat,ispin,iorb)
              enddo
              dmft_bath((Nbath-1)/2+1)%h(ilat,ilat,ispin,ispin,iorb,iorb)=0d0
           endif
           !
        enddo
     enddo
  enddo
  if(save_)call save_dmft_bath()
  call get_dmft_bath(bath_)
  call deallocate_dmft_bath()
end subroutine ph_symmetrize_bath_site
