  subroutine ed_get_cluster_density_matrix_single(dm,dm_eig,dm_rot,doprint)
    complex(8),dimension(4**Nimp,4**Nimp),intent(out)           :: dm
    complex(8),dimension(4**Nimp,4**Nimp),intent(out),optional  :: dm_rot
    real(8)   ,dimension(4**Nimp)        ,intent(out),optional  :: dm_eig
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
    !Cluster (impurity-problem) basis
    dm = cluster_density_matrix
    !
    !Diagonal (pure-state) basis
    if(present(dm_eig).and.present(dm_rot))then
      dm_rot=dm
      write(*,*) "DIAGONALIZING CLUSTER DENSITY MATRIX"
      call eigh(dm_rot,dm_eig,jobz='V',uplo='U')
    endif
    !
    !Print to file (if requested)
    if(doprint_)then
      if(present(dm_eig).and.present(dm_rot))then
         call print_cluster_dm(dm,dm_eig,dm_rot)
      else
         write(*,*) "CLUSTER DENSITY MATRIX IS ALREADY DIAGONAL"
         call print_cluster_dm(dm)
      endif
    endif
    !
  end subroutine ed_get_cluster_density_matrix_single


  subroutine print_cluster_dm(dm,dm_eig,dm_rot,ineq)
    !Passed
    complex(8),dimension(4**Nimp,4**Nimp),intent(in)            :: dm
    complex(8),dimension(4**Nimp,4**Nimp),intent(in),optional   :: dm_rot
    real(8)   ,dimension(4**Nimp)        ,intent(in),optional   :: dm_eig
    integer                              ,intent(in),optional   :: ineq
    !Internal
    integer                              :: unit_dm,unit_eig,unit_rot
    character(len=64)                    :: suffix_dm,suffix_eig,suffix_rot
    integer                              :: io,jo
    !
    if(present(ineq))then
      suffix_dm  = "cluster_density_matrix_ineq"//reg(str(ineq))//".dat"
      suffix_eig = "cluster_probabilities_ineq"//reg(str(ineq))//".dat"
      suffix_rot = "cluster_pure_states_ineq"//reg(str(ineq))//".dat"
    else
      suffix_dm  = "cluster_density_matrix.dat"
      suffix_eig = "cluster_probabilities.dat"
      suffix_rot = "cluster_pure_states.dat"
    endif
    !
    unit_dm = free_unit()
    open(unit_dm,file=suffix_dm,action="write",position="rewind",status='unknown')
    !
    write(unit_dm,"(A90)")"# cluster density matrix [REAL part]:"
    do io=1,4**Nimp
       write(unit_dm,"(90(F15.9,1X))") (real(dm(io,jo)),jo=1,4**Nimp)
    enddo
    write(unit_dm,*)
    !
    write(unit_dm,"(A90)")"# cluster density matrix [IMAG part]:"
    do io=1,4**Nimp
       write(unit_dm,"(90(F15.9,1X))") (aimag(dm(io,jo)),jo=1,4**Nimp)
    enddo
    write(unit_dm,*)
    !
    close(unit_dm)
    !
    if(present(dm_eig).and.present(dm_rot))then
       unit_eig = free_unit()
       open(unit_eig,file=suffix_eig,action="write",position="rewind",status='unknown') 
       !
       do io=1,4**Nimp
          write(unit_eig,"(90(F15.9,1X))") dm_eig(io)
       enddo
       write(unit_eig,*)
       !
       close(unit_eig)
       !
       unit_rot = free_unit()
       open(unit_rot,file=suffix_rot,action="write",position="rewind",status='unknown')
       !
       write(unit_rot,"(A90)")"# cluster pure states [REAL part]:"
       do io=1,4**Nimp
          write(unit_rot,"(90(F15.9,1X))") (real(dm_rot(io,jo)),jo=1,4**Nimp)
       enddo
       write(unit_rot,*)
       !
       write(unit_rot,"(A90)")"# cluster pure states [IMAG part]:"
       do io=1,4**Nimp
          write(unit_rot,"(90(F15.9,1X))") (aimag(dm_rot(io,jo)),jo=1,4**Nimp)
       enddo
       write(unit_rot,*)
       !
       close(unit_rot)
       !
    endif
    !
  end subroutine print_cluster_dm
