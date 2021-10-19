  subroutine ed_get_single_particle_density_matrix_lattice(dm_,custom_rot,dm_eig_,dm_rot_)
    implicit none
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:,:)
    complex(8),allocatable,intent(in) ,optional  :: custom_rot(:,:)
    real(8),allocatable,intent(out)   ,optional  :: dm_eig_(:,:)
    complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:,:)
    !internal
    integer                                      :: unit
    integer                                      :: iorb,jorb,ispin,jspin,io,jo,ilat,jlat,isite,Nlat,Nsites
    complex(8)                                   :: Tr
    complex(8),allocatable                       :: dm_custom_rot(:,:,:)
    real(8)                                      :: soc
    complex(8),allocatable                       :: dm_tmp(:,:)
    real(8),allocatable                          :: dm_eig_tmp(:)
    complex(8),allocatable                       :: dm_rot_tmp(:,:)
    complex(8),allocatable                       :: dm_custom_rot_tmp(:,:)
    !
    Nsites=size(single_particle_density_matrix_ii,1)
    !
    if (.not.allocated(single_particle_density_matrix)) then
       write(LOGfile,"(A)") "single_particle_density_matrix is not allocated"
       stop
    endif
    !
    if(allocated(dm_))deallocate(dm_)                             ;allocate(dm_(Nsites,Nlat*Nspin*Norb,Nlat*Nspin*Norb))          ;dm_ = zero
    if(allocated(dm_custom_rot))deallocate(dm_custom_rot)         ;allocate(dm_custom_rot(Nsites,Nlat*Nspin*Norb,Nlat*Nspin*Norb));dm_custom_rot = zero
    if(present(dm_eig_).and.allocated(dm_eig_))deallocate(dm_eig_);allocate(dm_eig_(Nsites,Nlat*Nspin*Norb))                      ;dm_eig_ = 0.0d0
    if(present(dm_rot_).and.allocated(dm_rot_))deallocate(dm_rot_);allocate(dm_rot_(Nsites,Nlat*Nspin*Norb,Nlat*Nspin*Norb))      ;dm_rot_ = zero
    !
    if(allocated(dm_tmp))deallocate(dm_)                          ;allocate(dm_tmp(Nlat*Nspin*Norb,Nlat*Nspin*Norb))              ;dm_tmp = zero
    if(allocated(dm_custom_rot_tmp))deallocate(dm_custom_rot_tmp) ;allocate(dm_custom_rot_tmp(Nlat*Nspin*Norb,Nlat*Nspin*Norb))   ;dm_custom_rot_tmp = zero
    if(allocated(dm_eig_tmp))deallocate(dm_eig_tmp)               ;allocate(dm_eig_tmp(Nlat*Nspin*Norb))                          ;dm_eig_tmp = 0.0d0
    if(allocated(dm_rot_tmp))deallocate(dm_rot_tmp)               ;allocate(dm_rot_tmp(Nlat*Nspin*Norb,Nlat*Nspin*Norb))          ;dm_rot_tmp = zero
    !
    do isite=1,Nsites
       !
       ! dm in the impurity problem basis
       dm_(isite,:,:) = nnn2lso_reshape(single_particle_density_matrix_ii(isite,:,:,:,:,:,:),Nlat,Nspin,Norb)
       !
       !
     ! dm in her diagonal basis
       if(present(dm_eig_).and.present(dm_rot_))then
         dm_rot_(isite,:,:)=dm_(isite,:,:)
         call eigh(dm_rot_(isite,:,:),dm_eig_(isite,:),jobz='V',uplo='U')
       endif
       !
       ! dm in the basis defined by custom_rot
       dm_custom_rot(isite,:,:)=matmul(transpose(conjg(custom_rot)),matmul(dm_(isite,:,:),custom_rot))
       !
       dm_tmp            = zero ; dm_tmp            = dm_(isite,:,:)
       dm_rot_tmp        = zero ; dm_rot_tmp        = dm_rot_(isite,:,:)
       dm_eig_tmp        = zero ; dm_eig_tmp        = dm_eig_(isite,:)
       dm_custom_rot_tmp = zero ; dm_custom_rot_tmp = dm_custom_rot(isite,:,:)
       call print_dm(dm_tmp,dm_rot_tmp,dm_eig_tmp,dm_custom_rot_tmp,isite)
       !
    enddo
  end subroutine ed_get_single_particle_density_matrix_lattice

