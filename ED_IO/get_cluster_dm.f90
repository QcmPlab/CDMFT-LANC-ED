  subroutine ed_get_cluster_density_matrix_single(dm,doprint)
    complex(8),dimension(4**Nimp,4**Nimp),intent(out)           :: dm
    logical                              ,intent(in) ,optional  :: doprint
    logical                                                     :: doprint_
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    if(.not.allocated(cluster_density_matrix))then
       write(LOGfile,"(A)") "cluster_density_matrix is not allocated"
       stop
    endif
    !
    dm = cluster_density_matrix
    !
    !Print to file (if requested)
    if(doprint_)then
       call print_cluster_dm(dm)
    endif
    !
  end subroutine ed_get_cluster_density_matrix_single



  
