subroutine ed_get_mag_lattice_1(yii,Nsites)
  integer                             :: Nsites
  real(8),dimension(Nsites,Nlat,Norb) :: yii
  yii=0d0
  if(allocated(mag_ineq))then
     if(Nlat>size(mag_ineq,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
     yii=mag_ineq
  endif
end subroutine ed_get_mag_lattice_1

