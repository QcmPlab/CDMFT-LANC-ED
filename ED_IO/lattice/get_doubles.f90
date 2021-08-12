 subroutine ed_get_doubles_lattice(yii,Nineq)
   integer                      :: Nineq
   real(8),dimension(Nineq,4)    :: yii
   yii=0d0
   if(allocated(ddii))then
      if(Nineq>size(ddii,1)) stop "ed_get_doubles error: required N_sites > evaluated N_sites"
      yii=ddii(:,:)
   endif
 end subroutine ed_get_doubles_lattice

 subroutine ed_get_dust_lattice(yii,Nineq)
   integer                 :: Nineq
   real(8),dimension(Nineq) :: yii
   yii=0d0
   if(allocated(ddii))then
      if(Nineq>size(ddii,1)) stop "ed_get_dust error: required N_sites > evaluated N_sites"
      yii=ddii(:,1)
   endif
 end subroutine ed_get_dust_lattice

 subroutine ed_get_dund_lattice(yii,Nineq)
   integer                 :: Nineq
   real(8),dimension(Nineq) :: yii
   yii=0d0
   if(allocated(ddii))then
      if(Nineq>size(ddii,1)) stop "ed_get_dund error: required N_sites > evaluated N_sites"
      yii=ddii(:,2)
   endif
 end subroutine ed_get_dund_lattice

 subroutine ed_get_dse_lattice(yii,Nineq)
   integer                 :: Nineq
   real(8),dimension(Nineq) :: yii
   yii=0d0
   if(allocated(ddii))then
      if(Nineq>size(ddii,1)) stop "ed_get_dse error: required N_sites > evaluated N_sites"
      yii=ddii(:,3)
   endif
 end subroutine ed_get_dse_lattice

 subroutine ed_get_dph_lattice(yii,Nineq)
   integer                 :: Nineq
   real(8),dimension(Nineq) :: yii
   yii=0d0
   if(allocated(ddii))then
      if(Nineq>size(ddii,1)) stop "ed_get_dph error: required N_sites > evaluated N_sites"
      yii=ddii(:,4)
   endif
 end subroutine ed_get_dph_lattice
