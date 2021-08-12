subroutine ed_get_epot_(eimp)
  real(8) :: eimp(Nlat)
  eimp = ed_Epot
end subroutine ed_get_epot_

subroutine ed_get_eint_(eimp)
  real(8) :: eimp(Nlat)
  eimp = ed_Eint
end subroutine ed_get_eint_

subroutine ed_get_ehartree_(eimp)
  real(8) :: eimp(Nlat)
  eimp = ed_Ehartree
end subroutine ed_get_ehartree_

subroutine ed_get_eknot_(eimp)
  real(8) :: eimp(Nlat)
  eimp = ed_Eknot
end subroutine ed_get_eknot_

