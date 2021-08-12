subroutine ed_get_dens_1(dens)
  real(8),dimension(Nlat,Norb) :: dens
  dens = ed_dens
end subroutine ed_get_dens_1

subroutine ed_get_dens_2(dens,ilat,iorb)
  real(8)   :: dens
  integer   :: ilat,iorb
  if(ilat>Nlat)stop "ed_get_dens error: lattice index > N_lattice"
  if(iorb>Norb)stop "ed_get_dens error: orbital index > N_orbital"
  dens = ed_dens(ilat,iorb)
end subroutine ed_get_dens_2

