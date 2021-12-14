  subroutine ed_get_single_particle_density_matrix_single(dm,doprint)
    complex(8),dimension(Nlat*Norb*Nspin,Nlat*Norb*Nspin),intent(out)            :: dm
    logical                                              ,intent(in) ,optional   :: doprint
    logical                                                                      :: doprint_
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    if(.not.allocated(single_particle_density_matrix))then
       write(LOGfile,"(A)") "single_particle_density_matrix is not allocated"
       stop
    endif
    !
    !Impurity problem basis
    dm = nnn2lso_reshape(single_particle_density_matrix,Nlat,Nspin,Norb)
    !
    !Print to file (if requested)
    if(doprint_)then
       call ed_print_dm(dm,Nlat*Norb*Nspin)
    endif
    !
  end subroutine ed_get_single_particle_density_matrix_single






