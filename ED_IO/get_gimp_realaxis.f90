!NORMAL, REAL GREEN'S FUNCTION
subroutine ed_get_gimp_real_1(Greal)
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:,:,:) = impGreal(:,:,:,:,:,:,:)
end subroutine ed_get_gimp_real_1

subroutine ed_get_gimp_real_2(Greal)
  complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal),intent(inout) :: Greal
  integer  :: io,jo,ilat,jlat,iorb,jorb,ispin,jspin
  do ilat=1,Nlat
    do jlat=1,Nlat
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Norb
               do jorb=1,Norb
                  io = index_stride_lso(ilat,ispin,iorb)
                  jo = index_stride_lso(jlat,jspin,jorb)
                  Greal(io,jo,:) = impGreal(ilat,jlat,ispin,jspin,iorb,jorb,:)
               enddo
            enddo
         enddo
      enddo
    enddo
  enddo
end subroutine ed_get_gimp_real_2

subroutine ed_get_gimp_real_3(Greal,ilat,jlat,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lreal),intent(inout) :: Greal
  integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin
  Greal(:) = impGreal(ilat,jlat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_gimp_real_3


