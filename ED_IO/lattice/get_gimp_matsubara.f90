!NORMAL, MATSUBARA GREEN'S FUNCTIONS
subroutine ed_get_gimp_matsubara_lattice_1(Gmats,Nsites)
  integer                                                                          :: Nsites
  complex(8),dimension(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(1:Nsites,:,:,:,:,:,:,:) = Gmats_ineq(1:Nsites,:,:,:,:,:,:,:)
end subroutine ed_get_gimp_matsubara_lattice_1

subroutine ed_get_gimp_matsubara_lattice_2(Gmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lmats),intent(inout) :: Gmats
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
                 Gmats(ilat,io,jo,:) = Gmats_ineq(isite,ilat,jlat,ispin,jspin,iorb,jorb,:)
                enddo
              enddo
             enddo
           enddo
        enddo
     enddo
  enddo
end subroutine ed_get_gimp_matsubara_lattice_2


