  
  if(allocated(diag_hybr))deallocate(diag_hybr)
  if(allocated(bath_diag))deallocate(bath_diag)

  Nfoo = Norb
  allocate(diag_hybr(Nlat,Nspin,Norb,Nbath));diag_hybr=0d0
  allocate(bath_diag(Nlat,Nspin,Nfoo,Nbath));bath_diag=0d0
  do ibath=1,Nbath
    do ispin=1,Nspin
      do iorb=1,Norb
        diag_hybr(ilat,ispin,iorb,ibath)=dmft_bath%item(ibath)%v
        bath_diag(ilat,ispin,iorb,ibath)=dmft_bath%item(ibath)%h(ilat,ilat,ispin,ispin,iorb,iorb)
      enddo
    enddo
  enddo

