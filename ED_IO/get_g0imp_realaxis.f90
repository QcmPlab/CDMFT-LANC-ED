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







!subroutine ed_get_g0imp_real_lattice_1(Greal,Nsites)
  !integer                                                                :: Nsites
  !complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  !Greal(1:Nsites,:,:,:,:,:) = G0realii(1:Nsites,:,:,:,:,:)
!end subroutine ed_get_g0imp_real_lattice_1

!subroutine ed_get_g0imp_real_lattice_2(Greal,Nsites)
  !integer                                                                :: Nsites
  !complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
  !integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  !do ilat=1,Nsites
     !do ispin=1,Nspin
        !do jspin=1,Nspin
           !do iorb=1,Norb
              !do jorb=1,Norb
                 !io = iorb + (ispin-1)*Norb
                 !jo = jorb + (jspin-1)*Norb
                 !Greal(ilat,io,jo,:) = G0realii(ilat,ispin,jspin,iorb,jorb,:)
              !enddo
           !enddo
        !enddo
     !enddo
  !enddo
!end subroutine ed_get_g0imp_real_lattice_2

!subroutine ed_get_g0imp_real_lattice_3(Greal,Nsites,ispin,jspin,iorb,jorb)
  !integer                                          :: Nsites
  !complex(8),dimension(Nsites,Lreal),intent(inout) :: Greal
  !integer                                          :: iorb,jorb,ispin,jspin
  !Greal(1:Nsites,:) = G0realii(1:Nsites,ispin,jspin,iorb,jorb,:)
!end subroutine ed_get_g0imp_real_lattice_3








!subroutine ed_get_g0imp_real_lattice_11(Greal,ilat)
  !integer                                                         :: ilat
  !complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  !Greal(:,:,:,:,:) = G0realii(ilat,:,:,:,:,:)
!end subroutine ed_get_g0imp_real_lattice_11

!subroutine ed_get_g0imp_real_lattice_21(Greal,ilat)
  !integer                                                         :: ilat
  !complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Greal
  !integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  !do ispin=1,Nspin
     !do jspin=1,Nspin
        !do iorb=1,Norb
           !do jorb=1,Norb
              !io = iorb + (ispin-1)*Norb
              !jo = jorb + (jspin-1)*Norb
              !Greal(io,jo,:) = G0realii(ilat,ispin,jspin,iorb,jorb,:)
           !enddo
        !enddo
     !enddo
  !enddo
!end subroutine ed_get_g0imp_real_lattice_21

!subroutine ed_get_g0imp_real_lattice_31(Greal,ilat,ispin,jspin,iorb,jorb)
  !integer                                   :: ilat
  !complex(8),dimension(Lreal),intent(inout) :: Greal
  !integer                                   :: iorb,jorb,ispin,jspin
  !Greal(:) = G0realii(ilat,ispin,jspin,iorb,jorb,:)
!end subroutine ed_get_g0imp_real_lattice_31
