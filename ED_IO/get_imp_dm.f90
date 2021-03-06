  subroutine ed_get_density_matrix_single(dm_,custom_rot,dm_eig_,dm_rot_)
    implicit none
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:)
    complex(8),allocatable,intent(in) ,optional  :: custom_rot(:,:)
    real(8),allocatable,intent(out)   ,optional  :: dm_eig_(:)
    complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:)
    !internal
    integer                                      :: unit
    integer                                      :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    complex(8)                                   :: Tr
    complex(8),allocatable                       :: dm_custom_rot(:,:)
    real(8)                                      :: soc
    !
    if (.not.allocated(imp_density_matrix)) then
       write(LOGfile,"(A)") "imp_density_matrix is not allocated"
       stop
    endif
    !
    if(allocated(dm_))                         deallocate(dm_)          ;allocate(dm_(Nlat*Nspin*Norb,Nlat*Nspin*Norb))          ;dm_ = zero
    if(allocated(dm_custom_rot))               deallocate(dm_custom_rot);allocate(dm_custom_rot(Nlat*Nspin*Norb,Nlat*Nspin*Norb));dm_custom_rot = zero
    if(present(dm_eig_).and.allocated(dm_eig_))deallocate(dm_eig_)      ;allocate(dm_eig_(Nlat*Nspin*Norb))                      ;dm_eig_ = 0.0d0
    if(present(dm_rot_).and.allocated(dm_rot_))deallocate(dm_rot_)      ;allocate(dm_rot_(Nlat*Nspin*Norb,Nlat*Nspin*Norb))      ;dm_rot_ = zero
    !
    ! dm in the impurity problem basis
    dm_ = nnn2lso_reshape(imp_density_matrix,Nlat,Nspin,Norb)
    !
    !
    ! dm in her diagonal basis
    if(present(dm_eig_).and.present(dm_rot_))then
      dm_rot_=dm_
      call eigh(dm_rot_,dm_eig_,jobz='V',uplo='U')
    endif
    !
    ! dm in the basis defined by custom_rot
    dm_custom_rot=matmul(transpose(conjg(custom_rot)),matmul(dm_,custom_rot))
    !
    !
    call print_dm(dm_,dm_rot_,dm_eig_,dm_custom_rot,1)
    !
  end subroutine ed_get_density_matrix_single



  subroutine ed_get_density_matrix_lattice(dm_,custom_rot,dm_eig_,dm_rot_)
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
    Nsites=size(imp_density_matrix_ii,1)
    !
    if (.not.allocated(imp_density_matrix)) then
       write(LOGfile,"(A)") "imp_density_matrix is not allocated"
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
       dm_(isite,:,:) = nnn2lso_reshape(imp_density_matrix_ii(isite,:,:,:,:,:,:),Nlat,Nspin,Norb)
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
  end subroutine ed_get_density_matrix_lattice


  subroutine print_dm(dm_,dm_rot_,dm_eig_,dm_custom_rot,ndx)
    implicit none
    integer               ,intent(in)            :: ndx
    complex(8),allocatable,intent(in)            :: dm_(:,:)
    complex(8),allocatable,intent(in)            :: dm_custom_rot(:,:)
    real(8),allocatable   ,intent(in),optional   :: dm_eig_(:)
    complex(8),allocatable,intent(in),optional   :: dm_rot_(:,:)
    !internal
    integer                                      :: unit
    character(len=24)                            :: suffix
    integer                                      :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    suffix="imp_density_matrix_"//reg(str(ndx))//".dat"
    !
    unit = free_unit()
    open(unit,file=suffix,action="write",position="rewind",status='unknown')
    !
    write(unit,"(A90)")"# density matrix in the impurity problem basis REAL part:"
    do io=1,Nlat*Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)
    !
    write(unit,"(A90)")"# density matrix in the impurity problem basis IMAGINARY part:"
    do io=1,Nlat*Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)
    !
    if(present(dm_eig_).and.present(dm_rot_))then
       write(unit,"(A90)")"# eigenvalues of density matrix"
       write(unit,'(10F22.12)') dm_eig_
       write(unit,*)
       !
       write(unit,"(A90)")"# density matrix eigenvector matrix REAL part:"
       do io=1,Nlat*Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_rot_(io,jo)),jo=1,Nlat*Nspin*Norb)
       enddo
       write(unit,*)
       !
       write(unit,"(A90)")"# density matrix eigenvector matrix IMAGINARY part:"
       do io=1,Nlat*Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_rot_(io,jo)),jo=1,Nlat*Nspin*Norb)
       enddo
       write(unit,*)
    endif
    !
    write(unit,"(A90)")"# density matrix in the basis defined by custom_rot REAL part:"
    do io=1,Nlat*Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_custom_rot(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)
    !
    write(unit,"(A90)")"# density matrix in the basis defined by custom_rot IMAGINARY part:"
    do io=1,Nlat*Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_custom_rot(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)
    write(unit,"(A30)")"# J basis densities"
    write(unit,"(90(F15.9,1X))") (real(dm_custom_rot(io,io)),io=1,Nlat*Nspin*Norb)
    !
    close(unit)
    !
  end subroutine print_dm
