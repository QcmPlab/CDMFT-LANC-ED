subroutine ed_get_dust_(docc)
  real(8) :: docc(Nlat)
  docc = ed_Dust
end subroutine ed_get_dust_

subroutine ed_get_dund_(docc)
  real(8) :: docc(Nlat)
  docc = ed_Dund
end subroutine ed_get_dund_

subroutine ed_get_dse_(docc)
  real(8) :: docc(Nlat)
  docc = ed_Dse
end subroutine ed_get_dse_

subroutine ed_get_dph_(docc)
  real(8) :: docc(Nlat)
  docc = ed_Dph
end subroutine ed_get_dph_

