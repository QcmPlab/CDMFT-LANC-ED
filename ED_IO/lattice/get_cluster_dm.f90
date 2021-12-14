  subroutine ed_get_cluster_density_matrix_lattice(dm,doprint)
    complex(8),allocatable,intent(out)           :: dm(:,:,:)
    logical               ,intent(in) ,optional  :: doprint
    logical                                      :: doprint_
    integer                                      :: isite,Nsites
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    Nsites=size(cluster_density_matrix_ii,1)
    !
    if(.not.allocated(cluster_density_matrix_ii))then
       write(LOGfile,"(A)") "cluster_density_matrix_ii is not allocated"
       stop
    endif
    !
    allocate(dm(Nsites,4**Nimp,4**Nimp)); dm=zero
    !
    do isite=1,Nsites
       !
       dm(isite,:,:) = cluster_density_matrix_ii(isite,:,:)
       !
       !Print to file [print_cluster_dm() is defined in ED_IO.f90]
       if(doprint_)then
          call print_cluster_dm(dm(isite,:,:),ineq=isite)
       endif
       !
    enddo
    !
    deallocate(dm)
    !
  end subroutine ed_get_cluster_density_matrix_lattice



