!NORMAL, REAL GREEN'S FUNCTION
subroutine ed_get_g0imp_real_1(Greal)
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:,:,:) = impG0real(:,:,:,:,:,:,:)
end subroutine ed_get_g0imp_real_1

subroutine ed_get_g0imp_real_2(Greal)
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
                  Greal(io,jo,:) = impG0real(ilat,jlat,ispin,jspin,iorb,jorb,:)
               enddo
            enddo
         enddo
      enddo
    enddo
  enddo
end subroutine ed_get_g0imp_real_2

subroutine ed_get_g0imp_real_3(Greal,ilat,jlat,ispin,jspin,iorb,jorb)
  complex(8),dimension(Lreal),intent(inout) :: Greal
  integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin
  Greal(:) = impG0real(ilat,jlat,ispin,jspin,iorb,jorb,:)
end subroutine ed_get_g0imp_real_3


subroutine ed_get_g0imp_real_lattice_1(Greal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(1:Nsites,:,:,:,:,:,:,:) = G0real_ineq(1:Nsites,:,:,:,:,:,:,:)
end subroutine ed_get_g0imp_real_lattice_1

subroutine ed_get_g0imp_real_lattice_2(Greal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
  integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat,jlat,isite
  do isite=1,Nsites
   do ilat=1,Nlat
    do jlat=1,Nlat
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = index_stride_lso(ilat,ispin,iorb)
                 jo = index_stride_lso(jlat,jspin,jorb)
                 Greal(isite,io,jo,:) = G0real_ineq(isite,ilat,jlat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
       enddo
     enddo
    enddo
  enddo
end subroutine ed_get_g0imp_real_lattice_2



