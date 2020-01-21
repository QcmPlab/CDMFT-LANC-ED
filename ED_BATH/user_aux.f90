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
  complex(8),allocatable         :: Hbasis_reshuffled(:,:,:,:,:,:,:)![Nlat][:][Nspin][:][Norb][:][Nsym]
  !
  if(.not.allocated(H_basis))STOP "check_bath_dimension: Hbasis not allocated"
  !
  if(.not.allocated(Hbasis_reshuffled))allocate(Hbasis_reshuffled(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(H_basis)))
  !
  do i=1,size(H_basis)
    Hbasis_reshuffled(:,:,:,:,:,:,i)=H_basis(i)%O
  enddo
  !
  Ntrue = get_bath_dimension(Hbasis_reshuffled)
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
  if(size(lambda_impHloc) .ne. dmft_bath%item(ibath)%N_dec)then
    dmft_bath%item(ibath)%lambda(dmft_bath%item(ibath)%N_dec)=offset
  else
    do isym=1,size(lambda_impHloc)
      if(is_identity(H_basis(isym)%O)) dmft_bath%item(ibath)%lambda(isym)=offset
      return
    enddo
  endif
  !
  !
  call get_dmft_bath(bath_)
  call deallocate_dmft_bath()
  !
end subroutine impose_bath_offset

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
!subroutine break_symmetry_bath_site(bath_,field,sign,save)
  !real(8),dimension(:) :: bath_
  !real(8)              :: field
  !real(8)              :: sign
  !logical,optional     :: save
  !logical              :: save_
  !save_=.true.;if(present(save))save_=save
  !call allocate_dmft_bath()
  !call set_dmft_bath(bath_)
  !do ibath=1,Nbath
     !do ilat=1,Nlat
        !do iorb=1,Norb
           !dmft_bath%item(ibath)%h(ilat,ilat,1,1,iorb,iorb)        = &
                !dmft_bath%item(ibath)%h(ilat,ilat,1,1,iorb,iorb)         + sign*field
           !dmft_bath%item(ibath)%h(ilat,ilat,Nspin,Nspin,iorb,iorb)= &
                !dmft_bath%item(ibath)%h(ilat,ilat,Nspin,Nspin,iorb,iorb) + sign*field
        !enddo
     !enddo
  !enddo
  !if(save_)call save_dmft_bath()
  !call get_dmft_bath(bath_)
  !call deallocate_dmft_bath()
!end subroutine break_symmetry_bath_site

!!---------------------------------------------------------!

!subroutine spin_symmetrize_bath_site(bath_,save)
  !real(8),dimension(:)   :: bath_
  !logical,optional       :: save
  !logical                :: save_
  !save_=.true.;if(present(save))save_=save
  !if(Nspin==1)then
     !write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
     !return
  !endif
  !!
  !call allocate_dmft_bath()
  !!call init_dmft_bathmask()
  !call set_dmft_bath(bath_)

  !do ibath=1,Nbath
     !do ilat=1,Nlat
        !do iorb=1,Norb
           !dmft_bath%item(ibath)%h(ilat,ilat,Nspin,Nspin,iorb,iorb) = dmft_bath%item(ibath)%h(ilat,ilat,1,1,iorb,iorb)
        !enddo
     !enddo
  !enddo
  !if(save_)call save_dmft_bath()
  !call get_dmft_bath(bath_)
  !call deallocate_dmft_bath()
!end subroutine spin_symmetrize_bath_site

!!---------------------------------------------------------!

!subroutine hermiticize_bath_main(bath_,save)
  !real(8),dimension(:)                               :: bath_
  !real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: h_aux
  !real(8)                                            :: trax
  !logical,optional                                   :: save
  !logical                                            :: save_
  !integer                                            :: io,jo
  !save_=.true.;if(present(save))save_=save
  !write(LOGfile,"(A)")"Hermiticizing bath:"
  !!
  !call allocate_dmft_bath()
  !!call init_dmft_bathmask()
  !call set_dmft_bath(bath_)
  !!
  !do ibath=1,Nbath
    !h_aux=nnn2lso_reshape(dmft_bath%item(ibath)%h,Nlat,Nspin,Norb)
    !do io=1,Nlat*Nspin*Norb
      !do jo=io+1,Nlat*Nspin*Norb
        !h_aux(io,jo)=(h_aux(io,jo)+h_aux(jo,io))*0.5d0
      !enddo
    !enddo
    !trax=trace(h_aux)/(Nlat*Nspin*Norb)
    !do io=1,Nlat*Nspin*Norb
      !h_aux(io,io)=trax
    !enddo
    !dmft_bath%item(ibath)%h=lso2nnn_reshape(h_aux,Nlat,Nspin,Norb)
  !enddo
  !!
  !if(save_)call save_dmft_bath()
  !call get_dmft_bath(bath_)
  !call deallocate_dmft_bath()
!end subroutine hermiticize_bath_main
!!---------------------------------------------------------!


!subroutine orb_equality_bath_site(bath_,indx,save)
  !real(8),dimension(:)   :: bath_
  !integer,optional       :: indx
  !logical,optional       :: save
  !integer                :: indx_,ibath,iorb
  !logical                :: save_
  !indx_=1     ;if(present(indx))indx_=indx
  !save_=.true.;if(present(save))save_=save
  !if(Norb==1)then
     !write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
     !return
  !endif
  !!
  !call allocate_dmft_bath()
!!  call init_dmft_bathmask()
  !call set_dmft_bath(bath_)
  !!
  !do ibath=1,Nbath
     !do iorb=1,Norb
        !if(iorb==indx_)cycle
        !dmft_bath%item(ibath)%h(:,:,:,:,iorb,iorb)=dmft_bath%item(ibath)%h(:,:,:,:,indx_,indx_)
     !enddo
  !enddo
  !!
  !if(save_)call save_dmft_bath()
  !call get_dmft_bath(bath_)
  !call deallocate_dmft_bath()
!end subroutine orb_equality_bath_site



!subroutine ph_symmetrize_bath_site(bath_,save)
  !real(8),dimension(:)   :: bath_
  !integer                :: ibath,ilat,ispin,iorb
  !logical,optional       :: save
  !logical                :: save_
  !save_=.true.;if(present(save))save_=save
  !call allocate_dmft_bath()
  !call set_dmft_bath(bath_)
  !if(Nbath==1)return
  !!

  !!
  !if(mod(Nbath,2)==0)then
     !do ibath=1,Nbath/2
        !forall(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)&
             !dmft_bath%item(Nbath+1-ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)=-dmft_bath%item(ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)
        !dmft_bath%item(Nbath+1-ibath)%v= dmft_bath%item(ibath)%v
     !enddo
  !else
     !do ibath=1,(Nbath-1)/2
        !forall(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)&
             !dmft_bath%item(Nbath+1-ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)=-dmft_bath%item(ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)
        !dmft_bath%item(Nbath+1-ibath)%v= dmft_bath%item(ibath)%v
     !enddo
     !dmft_bath%item((Nbath-1)/2+1)%h(ilat,ilat,ispin,ispin,iorb,iorb)=0d0
  !endif
  !!
  !if(save_)call save_dmft_bath()
  !call get_dmft_bath(bath_)
  !call deallocate_dmft_bath()
!end subroutine ph_symmetrize_bath_site
