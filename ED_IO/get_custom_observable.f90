subroutine ed_get_custom_observable_1(value,obs_name)
  real(8)                                                 :: value
  character(len=*)                                        :: obs_name
  integer                                                 :: counter
  !
  if(MpiMaster .and. custom_o%init)then
    do counter=1,custom_o%N_filled
      if (custom_o%item(counter)%o_name .eq. obs_name) value=custom_o%item(counter)%o_value
    enddo
  else
    write(*,*)"Can't get custom observable value"
  endif
end subroutine ed_get_custom_observable_1


