  subroutine ed_get_cluster_density_matrix_single(dm,doprint)
    complex(8),allocatable,intent(out)           :: dm(:,:)
    logical               ,intent(in) ,optional  :: doprint
    logical                                      :: doprint_
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    if(.not.allocated(cluster_density_matrix))then
       stop "ERROR: cluster_density_matrix is not allocated"
    endif
    !
    dm = cluster_density_matrix
    !
    !Print to file (if requested)
    if(doprint_)then
       call ed_print_dm(dm,4**Nimp)
    endif
    !
  end subroutine ed_get_cluster_density_matrix_single






