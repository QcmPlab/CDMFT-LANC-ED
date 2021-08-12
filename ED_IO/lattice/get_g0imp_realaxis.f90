!NORMAL, REAL GREEN'S FUNCTION
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



