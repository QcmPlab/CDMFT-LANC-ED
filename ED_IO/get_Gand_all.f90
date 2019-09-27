
function delta_bath_main_(x,bath_) result(Delta)
  complex(8),dimension(:),intent(in)                            :: x
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  real(8),dimension(:)                                          :: bath_
  logical                                                       :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath()
  call set_dmft_bath(bath_)
  Delta = delta_bath(x)
  call deallocate_dmft_bath()
end function delta_bath_main_



function g0and_bath_main_(x,bath_) result(G0and)
  complex(8),dimension(:),intent(in)                            :: x
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  real(8),dimension(:)                                          :: bath_
  logical                                                       :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath()
  call set_dmft_bath(bath_)
  G0and = g0and_bath(x)
  call deallocate_dmft_bath()
end function g0and_bath_main_


function invg0_bath_main_(x,bath_) result(G0and)
  complex(8),dimension(:),intent(in)                            :: x
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  real(8),dimension(:)                                          :: bath_
  logical                                                       :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "invg0_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath()
  call set_dmft_bath(bath_)
  G0and = invg0_bath(x)
  call deallocate_dmft_bath()
end function invg0_bath_main_
