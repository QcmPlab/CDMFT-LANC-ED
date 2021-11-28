  subroutine ed_get_single_particle_density_matrix_single(dm,dm_eig,dm_rot,doprint)
    complex(8),dimension(Nlat*Norb*Nspin,Nlat*Norb*Nspin),intent(out)            :: dm
    complex(8),dimension(Nlat*Norb*Nspin,Nlat*Norb*Nspin),intent(out),optional   :: dm_rot
    real(8)   ,dimension(Nlat*Norb*Nspin)                ,intent(out),optional   :: dm_eig
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
    !Diagonal (pure-state) basis
    if(present(dm_eig).and.present(dm_rot))then
      dm_rot=dm
      call eigh(dm_rot,dm_eig,jobz='V',uplo='U')
    endif
    !
    !Print to file (if requested)
    if(doprint_)then
      if(present(dm_eig).and.present(dm_rot))then
         call print_sp_dm(dm,dm_eig,dm_rot)
      else
         call print_sp_dm(dm)
      endif
    endif
    !
  end subroutine ed_get_single_particle_density_matrix_single


  subroutine print_sp_dm(dm,dm_eig,dm_rot,ineq)
    !Passed
    complex(8),dimension(Nlat*Norb*Nspin,Nlat*Norb*Nspin),intent(in)            :: dm
    complex(8),dimension(Nlat*Norb*Nspin,Nlat*Norb*Nspin),intent(in),optional   :: dm_rot
    real(8)   ,dimension(Nlat*Norb*Nspin)                ,intent(in),optional   :: dm_eig
    integer                                              ,intent(in),optional   :: ineq
    !Internal
    integer                                              :: unit
    character(len=64)                                    :: suffix
    integer                                              :: io,jo
    !
    if(present(ineq))then
      suffix="sp_density_matrix_ineq"//reg(str(ineq))//"_info.dat"
    else
      suffix="sp_density_matrix_info.dat"
    endif
    !
    unit = free_unit()
    open(unit,file=suffix,action="write",position="rewind",status='unknown')
    !
    write(unit,"(A90)")"# single-particle density matrix in the impurity problem basis [REAL part]:"
    do io=1,Nlat*Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)
    !
    write(unit,"(A90)")"# single-particle density matrix in the impurity problem basis [IMAG part]:"
    do io=1,Nlat*Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm(io,jo)),jo=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)
    !
    if(present(dm_eig).and.present(dm_rot))then
       write(unit,"(A90)")"# eigenvalues of single-particle density matrix"
       do io=1,Nlat*Nspin*Norb
          write(unit,"(90(F15.9,1X))") dm_eig(io)
       enddo
       write(unit,*)
       !
       write(unit,"(A90)")"# eigenvector matrix [REAL part]:"
       do io=1,Nlat*Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_rot(io,jo)),jo=1,Nlat*Nspin*Norb)
       enddo
       write(unit,*)
       !
       write(unit,"(A90)")"# eigenvector matrix [IMAG part]:"
       do io=1,Nlat*Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_rot(io,jo)),jo=1,Nlat*Nspin*Norb)
       enddo
       write(unit,*)
    endif
    !
    close(unit)
    !
  end subroutine print_sp_dm
