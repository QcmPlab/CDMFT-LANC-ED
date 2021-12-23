  !+---------------------------------------------------------------------+
  !PURPOSE :   further reduce the cdms by tracing out (Nlat-Nsites) sites
  !+---------------------------------------------------------------------+
  subroutine ed_get_reduced_density_matrix_lattice(rdm,Nsites,doprint)
    complex(8),allocatable                         :: cdm(:,:,:)
    complex(8),allocatable,intent(out)             :: rdm(:,:,:)
    integer               ,intent(in)              :: Nsites
    logical               ,intent(in),optional     :: doprint
    logical                                        :: doprint_
    integer                                        :: ii,Nineq
    integer    :: i,j,io,jo,iUP,iDW,jUP,jDW
    integer    :: iIMPup,iIMPdw,jIMPup,jIMPdw
    integer    :: iREDup,iREDdw,jREDup,jREDdw
    integer    :: iTrUP,iTrDW,jTrUP,jTrDW,Nred
    !
    if(Nsites>Nlat)then
      stop "ERROR: cannot request a density matrix reduced to more sites then Nlat"
    endif
    !
    doprint_=.false.; if(present(doprint)) doprint_=doprint
    !
    if(.not.allocated(cluster_density_matrix_ii))then
      stop "ERROR: cluster_density_matrix_ii is not allocated"
    endif
    !
    Nineq=size(cluster_density_matrix_ii,1)
    allocate(cdm(Nineq,4**Nimp,4**Nimp)); cdm=zero
    Nred=Norb*Nsites
    allocate(rdm(Nineq,4**Nred,4**Nred)); rdm=zero
    !
    do ii=1,Nineq
       !
       cdm(ii,:,:) = cluster_density_matrix_ii(ii,:,:)
       !
       !Trace the cdm(ii,:,:) to the requested Nsites-rdm(ii,:,:)
       do iUP = 1,2**Nimp
          do iDW = 1,2**Nimp
             i = iUP + (iDW-1)*2**Nimp
             iIMPup = iup-1
             iIMPdw = idw-1
             iREDup = Ibits(iIMPup,0,Nred)
             iREDdw = Ibits(iIMPdw,0,Nred)
             iTrUP  = Ibits(iIMPup,Nred,Nimp)
             iTrDW  = Ibits(iIMPdw,Nred,Nimp)
             do jUP = 1,2**Nimp
                do jDW = 1,2**Nimp
                   j = jUP + (jDW-1)*2**Nimp
                   jIMPup = jup-1
                   jIMPdw = jdw-1
                   jREDup = Ibits(jIMPup,0,Nred)
                   jREDdw = Ibits(jIMPdw,0,Nred)
                   jTrUP  = Ibits(jIMPup,Nred,Nimp)
                   jTrDW  = Ibits(jIMPdw,Nred,Nimp)
                   if(jTrUP/=iTrUP.or.jTrDW/=iTrDW)cycle
                   io = (iREDup+1) + iREDdw*2**Nred
                   jo = (jREDup+1) + jREDdw*2**Nred
                   rdm(ii,io,jo) = rdm(ii,io,jo) + cdm(ii,i,j)
                enddo
             enddo
          enddo
       enddo
       !
       !Print to file (if requested)
       if(doprint_)then
         call ed_print_dm(rdm(ii,:,:),4**Nred,ineq=ii)
       endif
       !
    enddo
    !
    deallocate(cdm)
    !
  end subroutine ed_get_reduced_density_matrix_lattice

! NB: the spin-factorization of the loops is paramount,
!     since we build in such way the original cdm and the
!     ordering of the basis states has to be preserved 
!     in order to trace along the proper matrix elements 



