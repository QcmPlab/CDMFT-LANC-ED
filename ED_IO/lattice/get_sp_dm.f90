  subroutine ed_get_single_particle_density_matrix_lattice(dm,doprint)
    complex(8),allocatable,intent(out)           :: dm(:,:,:)
    logical               ,intent(in),optional   :: doprint
    logical                                      :: doprint_
    integer                                      :: ii,Nineq
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    if(.not.allocated(single_particle_density_matrix_ii))then
       stop "ERROR: single_particle_density_matrix_ii is not allocated"
    endif
    !
    Nineq=size(single_particle_density_matrix_ii,1)
    allocate(dm(Nineq,Nlat*Nspin*Norb,Nlat*Nspin*Norb)); dm=zero
    !
    do ii=1,Nineq
       !
       !Impurity problem basis
       dm(ii,:,:) = nnn2lso_reshape(single_particle_density_matrix_ii(ii,:,:,:,:,:,:),Nlat,Nspin,Norb)
       !
       !Print to file (if requested)
       if(doprint_)then
          call ed_print_dm(dm(ii,:,:),Nlat*Nspin*Norb,ineq=ii)
       endif
       !
    enddo
    !
    deallocate(dm)
    !
  end subroutine ed_get_single_particle_density_matrix_lattice


  
