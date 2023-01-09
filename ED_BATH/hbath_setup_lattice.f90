subroutine init_Hbath_symmetries_lattice(Hvec,lambdavec)
  complex(8),dimension(:,:,:,:,:,:,:) :: Hvec      ![size(Hloc),Nsym]
  real(8),dimension(:,:,:)            :: lambdavec ![Nsites,Nbath,Nsym]
  integer                             :: isym,isites,ibath,Nsym,Nsites
  !
  Nsites = size(lambdavec,1)
  Nsym   = size(lambdavec,3)
  call assert_shape(Hvec,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nsym],"init_Hbath_symmetries","Hvec")
  !
  if(allocated(Hbath_lambda_ineq))deallocate(Hbath_lambda_ineq)
  allocate(Hbath_lambda_ineq(Nsites,Nbath,Nsym))
  call allocate_Hbath(Nsym)
  !
  do isym=1,Nsym
     Hbath_lambda_ineq(:,:,isym)  = lambdavec(:,:,isym)
     Hbath_basis(isym)%O = Hvec(:,:,:,:,:,:,isym)
  enddo
  !
  if(ed_verbose>2)then
     do isites=1,Nsites
        write(*,*) "Inequivalent #"//str(isites)//":"
        do ibath=1,Nbath
           write(*,*) "> Hbath #"//str(ibath)//":"
           call print_hloc(Hbath_build(Hbath_lambda_ineq(isites,ibath,:)))
        enddo
     enddo
  endif
  !
end subroutine init_Hbath_symmetries_lattice


!+-------------------------------------------------------------------+
!PURPOSE  : take Hbath for i-th site in the real-space case
!+-------------------------------------------------------------------+

  subroutine Hbath_site(site)
    integer :: site
    if(site<1.OR.site>size(Hbath_lambda_ineq,1))stop "ERROR Hbath_site: site not in [1,Nlat]"
    if(.not.allocated(Hbath_lambda_ineq))stop "ERROR Hbath_site: Hbath_lambda_ineq not allocated"
    Hbath_lambda(:,:) = Hbath_lambda_ineq(site,:,:)
  end subroutine Hbath_site


