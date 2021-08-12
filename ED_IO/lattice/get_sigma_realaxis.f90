!NORMAL, REAL SELF-ENERGY
subroutine ed_get_sigma_real_lattice_1(Sreal,Nsites)
  integer                                                                          :: Nsites
  complex(8),dimension(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(1:Nsites,:,:,:,:,:,:,:) = Sreal_ineq(1:Nsites,:,:,:,:,:,:,:)
end subroutine ed_get_sigma_real_lattice_1

subroutine ed_get_sigma_real_lattice_2(Sreal,Nsites)
  integer                                                                          :: Nsites
  complex(8),dimension(Nsites,Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal),intent(inout) :: Sreal
  integer                                                                          :: io,jo,iorb,jorb,ispin,jspin,isite,ilat,jlat
  do isite=1,Nsites
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = index_stride_lso(ilat,ispin,iorb)
                 jo = index_stride_lso(jlat,jspin,jorb)
                 Sreal(isite,io,jo,:) = Sreal_ineq(isite,ilat,jlat,ispin,jspin,iorb,jorb,:)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_sigma_real_lattice_2

