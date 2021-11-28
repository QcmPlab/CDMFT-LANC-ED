  subroutine ed_get_cluster_density_matrix_lattice(dm,dm_eig,dm_rot,doprint)
    complex(8),allocatable,intent(out)           :: dm(:,:,:)
    complex(8),allocatable,intent(out),optional  :: dm_rot(:,:,:)
    real(8)   ,allocatable,intent(out),optional  :: dm_eig(:,:)
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
    if(allocated(dm))deallocate(dm)                             ;allocate(dm(Nsites,4**Nimp,4**Nimp))    ;dm     = zero
    if(present(dm_eig).and.allocated(dm_eig))deallocate(dm_eig);allocate(dm_eig(Nsites,4**Nimp))        ;dm_eig = zero
    if(present(dm_rot).and.allocated(dm_rot))deallocate(dm_rot);allocate(dm_rot(Nsites,4**Nimp,4**Nimp));dm_rot = zero
    !
    do isite=1,Nsites
       !
       !Impurity problem basis
       dm(isite,:,:) = cluster_density_matrix_ii(isite,:,:)
       !
       !Diagonal (pure-state) basis
       if(present(dm_eig).and.present(dm_rot))then
         dm_rot(isite,:,:)= dm(isite,:,:)
         call eigh(dm_rot(isite,:,:),dm_eig(isite,:),jobz='V',uplo='U')
       endif
       !
       !Print to file [print_sp_dm() is defined in ED_IO/get_sp_dm]
       if(doprint_)then
          if(present(dm_eig).and.present(dm_rot))then
            call print_cluster_dm(dm(isite,:,:),dm_eig(isite,:),dm_rot(isite,:,:),ineq=isite)
          else
            call print_cluster_dm(dm(isite,:,:),ineq=isite)
          endif
       endif
       !
    enddo
    !
  end subroutine ed_get_cluster_density_matrix_lattice



