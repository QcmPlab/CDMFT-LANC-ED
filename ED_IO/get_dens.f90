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

 subroutine ed_get_dens_lattice_1(yii,Nineq)
   integer                            :: Nineq
   real(8),dimension(Nineq,Nlat,Norb) :: yii
   yii=0d0    
   if(allocated(dens_ineq))then
      if(Nineq>size(dens_ineq,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
      yii=dens_ineq
   end if
 end subroutine ed_get_dens_lattice_1

