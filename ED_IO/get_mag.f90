subroutine ed_get_mag_1(mag) 
  real(8),dimension(Nlat,Norb) :: mag
  mag = (ed_dens_up-ed_dens_dw)
end subroutine ed_get_mag_1

subroutine ed_get_mag_2(mag,ilat,iorb) 
  real(8)   :: mag
  integer   :: ilat,iorb
  if(ilat>Nlat)stop "ed_get_mag error: lattice index > N_lattice"
  if(iorb>Norb)stop "ed_get_mag error: orbital index > N_orbital"
  mag = ed_dens_up(ilat,iorb)-ed_dens_dw(ilat,iorb)
end subroutine ed_get_mag_2


