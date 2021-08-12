subroutine ed_get_docc_1(docc) 
  real(8),dimension(Nlat,Norb) :: docc
  docc = ed_docc
end subroutine ed_get_docc_1

subroutine ed_get_docc_2(docc,ilat,iorb) 
  real(8)   :: docc
  integer   :: ilat,iorb
  if(ilat>Nlat)stop "ed_get_docc error: lattice index > N_lattice"
  if(iorb>Norb)stop "ed_get_docc error: orbital index > N_orbital"
  docc = ed_docc(ilat,iorb)
end subroutine ed_get_docc_2

