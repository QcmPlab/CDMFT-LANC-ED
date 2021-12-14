  subroutine ed_get_single_particle_density_matrix_lattice(dm,doprint)
    complex(8),allocatable,intent(out)           :: dm(:,:,:)
    logical               ,intent(in),optional   :: doprint
    logical                                      :: doprint_
    integer                                      :: isite,Nsites
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    Nsites=size(single_particle_density_matrix_ii,1)
    !
    if(.not.allocated(single_particle_density_matrix_ii))then
       write(LOGfile,"(A)") "single_particle_density_matrix_ii is not allocated"
       stop
    endif
    !
    allocate(dm(Nsites,Nlat*Nspin*Norb,Nlat*Nspin*Norb)); dm=zero
    !
    do isite=1,Nsites
       !
       !Impurity problem basis
       dm(isite,:,:) = nnn2lso_reshape(single_particle_density_matrix_ii(isite,:,:,:,:,:,:),Nlat,Nspin,Norb)
       !
       !Print to file [ed_print_dm() is defined in ED_IO.f90]
       if(doprint_)then
          call ed_print_dm(dm(isite,:,:),Nlat*Nspin*Norb,ineq=isite)
       endif
       !
    enddo
    !
    deallocate(dm)
    !
  end subroutine ed_get_single_particle_density_matrix_lattice


  
