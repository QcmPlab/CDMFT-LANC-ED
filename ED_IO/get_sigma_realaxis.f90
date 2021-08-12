!NORMAL, REAL SELF-ENERGY
subroutine ed_get_sigma_real_1(Sreal)
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(:,:,:,:,:,:,:) = impSreal(:,:,:,:,:,:,:)
end subroutine ed_get_sigma_real_1

subroutine ed_get_sigma_real_2(Sreal)
  complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal),intent(inout) :: Sreal
  integer  :: io,jo,ilat,jlat,iorb,jorb,ispin,jspin
  do ilat=1,Nlat
    do jlat=1,Nlat
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Norb
               do jorb=1,Norb
                  io = index_stride_lso(ilat,ispin,iorb) 
                  jo = index_stride_lso(jlat,jspin,jorb)
                  Sreal(io,jo,:) = impSreal(ilat,jlat,ispin,jspin,iorb,jorb,:)
               enddo
            enddo
         enddo
      enddo
    enddo
  enddo
end subroutine ed_get_sigma_real_2

subroutine ed_get_sigma_real_3(Sreal,ilat,jlat,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lreal),intent(inout) :: Sreal
  integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin
  Sreal(:) = impSreal(ilat,jlat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_sigma_real_3



