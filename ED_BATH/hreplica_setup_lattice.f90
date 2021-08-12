subroutine init_Hreplica_symmetries_lattice(Hvec,lambdavec)
  complex(8),dimension(:,:,:,:,:,:,:) :: Hvec
  real(8),dimension(:,:)              :: lambdavec ![Nsites,Nsym]
  integer                             :: isym,isites,N,Nsites
  !
  Nsites=size(lambdavec,1)
  N     =size(lambdavec,2)
  call assert_shape(Hvec,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,N],"init_Hreplica_symmetries","Hvec")
  !
  if(allocated(Hreplica_lambda_ineq))deallocate(Hreplica_lambda_ineq)
  allocate(Hreplica_lambda_ineq(Nsites,N))
  call allocate_hreplica(N)
  !
  do isym=1,N
     Hreplica_lambda_ineq(:,isym)  = lambdavec(:,isym)
     Hreplica_basis(isym)%O = Hvec(:,:,:,:,:,:,isym)
  enddo
  !
  if(ed_verbose>2)then
     do isites=1,Nsites
        call print_hloc(Hreplica_build(Hreplica_lambda_ineq(isites,:)))
     enddo
  endif
end subroutine init_Hreplica_symmetries_lattice

