!NORMAL, MATSUBARA GREEN'S FUNCTIONS
subroutine ed_get_g0imp_matsubara_1(Gmats)
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:,:,:) = impG0mats(:,:,:,:,:,:,:)
end subroutine ed_get_g0imp_matsubara_1

subroutine ed_get_g0imp_matsubara_2(Gmats)
  complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats),intent(inout) :: Gmats
  integer  :: io,jo,ilat,jlat,iorb,jorb,ispin,jspin
  do ilat=1,Nlat
    do jlat=1,Nlat
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Norb
               do jorb=1,Norb
                  io = index_stride_lso(ilat,ispin,iorb)
                  jo = index_stride_lso(jlat,jspin,jorb)
                  Gmats(io,jo,:) = impG0mats(ilat,jlat,ispin,jspin,iorb,jorb,:)
               enddo
            enddo
         enddo
      enddo
    enddo
  enddo
end subroutine ed_get_g0imp_matsubara_2

subroutine ed_get_g0imp_matsubara_3(Gmats,ilat,jlat,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lmats),intent(inout) :: Gmats
  integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin
  Gmats(:) = impG0mats(ilat,jlat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_g0imp_matsubara_3


