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





!subroutine ed_get_sigma_real_lattice_1(Sreal,Nsites)
  !integer                                                                :: Nsites
  !complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  !Sreal(1:Nsites,:,:,:,:,:) = Srealii(1:Nsites,:,:,:,:,:)
!end subroutine ed_get_sigma_real_lattice_1

!subroutine ed_get_sigma_real_lattice_2(Sreal,Nsites)
  !integer                                                                :: Nsites
  !complex(8),dimension(Nsites,Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
  !integer                                                                :: io,jo,iorb,jorb,ispin,jspin,ilat
  !do ilat=1,Nsites
     !do ispin=1,Nspin
        !do jspin=1,Nspin
           !do iorb=1,Norb
              !do jorb=1,Norb
                 !io = iorb + (ispin-1)*Norb
                 !jo = jorb + (jspin-1)*Norb
                 !Sreal(ilat,io,jo,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
              !enddo
           !enddo
        !enddo
     !enddo
  !enddo
!end subroutine ed_get_sigma_real_lattice_2

!subroutine ed_get_sigma_real_lattice_3(Sreal,Nsites,ispin,jspin,iorb,jorb)
  !integer                                          :: Nsites
  !complex(8),dimension(Nsites,Lreal),intent(inout) :: Sreal
  !integer                                          :: iorb,jorb,ispin,jspin
  !Sreal(1:Nsites,:) = Srealii(1:Nsites,ispin,jspin,iorb,jorb,:)
!end subroutine ed_get_sigma_real_lattice_3





!subroutine ed_get_sigma_real_lattice_11(Sreal,ilat)
  !integer                                                         :: ilat
  !complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  !Sreal(:,:,:,:,:) = Srealii(ilat,:,:,:,:,:)
!end subroutine ed_get_sigma_real_lattice_11

!subroutine ed_get_sigma_real_lattice_21(Sreal,ilat)
  !integer                                                         :: ilat
  !complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal),intent(inout) :: Sreal
  !integer                                                         :: io,jo,iorb,jorb,ispin,jspin
  !do ispin=1,Nspin
     !do jspin=1,Nspin
        !do iorb=1,Norb
           !do jorb=1,Norb
              !io = iorb + (ispin-1)*Norb
              !jo = jorb + (jspin-1)*Norb
              !Sreal(io,jo,:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
           !enddo
        !enddo
     !enddo
  !enddo
!end subroutine ed_get_sigma_real_lattice_21

!subroutine ed_get_sigma_real_lattice_31(Sreal,ilat,ispin,jspin,iorb,jorb)
  !integer                                   :: ilat
  !complex(8),dimension(Lreal),intent(inout) :: Sreal
  !integer                                   :: iorb,jorb,ispin,jspin
  !Sreal(:) = Srealii(ilat,ispin,jspin,iorb,jorb,:)
!end subroutine ed_get_sigma_real_lattice_31
