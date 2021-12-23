  subroutine ed_get_cluster_density_matrix_lattice(dm,doprint)
    complex(8),allocatable,intent(out)           :: dm(:,:,:)
    logical               ,intent(in) ,optional  :: doprint
    logical                                      :: doprint_
    integer                                      :: ii,Nineq
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    if(.not.allocated(cluster_density_matrix_ii))then
       stop "ERROR: cluster_density_matrix_ii is not allocated"
    endif
    !
    Nineq=size(cluster_density_matrix_ii,1)
    allocate(dm(Nineq,4**Nimp,4**Nimp)); dm=zero
    !
    do ii=1,Nineq
       !
       dm(ii,:,:) = cluster_density_matrix_ii(ii,:,:)
       !
       !Print to file (if requested)
       if(doprint_)then
          call ed_print_dm(dm(ii,:,:),4**Nimp,ineq=ii)
       endif
       !
    enddo
    !
    deallocate(dm)
    !
  end subroutine ed_get_cluster_density_matrix_lattice



  
