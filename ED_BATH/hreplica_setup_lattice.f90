subroutine init_Hreplica_symmetries_lattice(Hvec,lambdavec)
  complex(8),dimension(:,:,:,:,:,:,:) :: Hvec      ![size(Hloc),Nsym]
  real(8),dimension(:,:,:)            :: lambdavec ![Nsites,Nbath,Nsym]
  integer                             :: isym,isites,ibath,Nsym,Nsites
  !
  Nsites = size(lambdavec,1)
  Nsym   = size(lambdavec,3)
  call assert_shape(Hvec,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
  !
  if(allocated(Hreplica_lambda_ineq))deallocate(Hreplica_lambda_ineq)
  allocate(Hreplica_lambda_ineq(Nsites,Nbath,Nsym))
  call allocate_hreplica(Nsym)
  !
  do isym=1,Nsym
     Hreplica_lambda_ineq(:,:,isym)  = lambdavec(:,:,isym)
     Hreplica_basis(isym)%O = Hvec(:,:,:,:,:,:,isym)
  enddo
  !
  if(ed_verbose>2)then
     do isites=1,Nsites
        write(*,*) "Inequivalent #"//str(isites)//":"
        do ibath=1,Nbath
           write(*,*) "> Hreplica #"//str(ibath)//":"
           call print_hloc(Hreplica_build(Hreplica_lambda_ineq(isites,ibath,:)))
        enddo
     enddo
  endif
  !
end subroutine init_Hreplica_symmetries_lattice


!+-------------------------------------------------------------------+
!PURPOSE  : take Hreplica for i-th site in the real-space case
!+-------------------------------------------------------------------+

  subroutine Hreplica_site(site)
    integer :: site
    if(site<1.OR.site>size(Hreplica_lambda_ineq,1))stop "ERROR Hreplica_site: site not in [1,Nlat]"
    if(.not.allocated(Hreplica_lambda_ineq))stop "ERROR Hreplica_site: Hreplica_lambda_ineq not allocated"
    Hreplica_lambda(:,:) = Hreplica_lambda_ineq(site,:,:)
  end subroutine Hreplica_site


