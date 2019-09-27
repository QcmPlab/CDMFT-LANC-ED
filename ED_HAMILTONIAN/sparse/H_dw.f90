  do idw=1,DimDw
     mdw  = Hs(2)%map(idw)
     ibdw  = bdecomp(mdw,Ns)
     !
     !
     !> H_imp: Off-diagonal elements, i.e. non-local part. 
     !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
     do ilat=1,Nlat
        do jlat=1,Nlat
           do iorb=1,Norb
              do jorb=1,Norb
                 is = imp_state_index(ilat,iorb) 
                 js = imp_state_index(jlat,jorb)
                 Jcondition = (impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)/=0d0)&
                      .AND.(ibdw(js)==1).AND.(ibdw(is)==0)
                 if (Jcondition) then
                    call c(js,mdw,k1,sg1)
                    call cdg(is,k1,k2,sg2)
                    jdw = binary_search(Hs(2)%map,k2)
                    htmp = impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)*sg1*sg2
                    !
                    call sp_insert_element(spH0dws(1),htmp,idw,jdw)
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
     !
     !> H_Bath: inter-orbital bath hopping contribution.
     do ibath=1,Nbath
        do ilat=1,Nlat 
           do jlat=1,Nlat
              do iorb=1,Norb
                 do jorb=1,Norb
                    !
                    ialfa = getBathStride(ilat,iorb,ibath)
                    ibeta = getBathStride(jlat,jorb,ibath)
                    Jcondition = &
                         (dmft_bath%item(ibath)%h(ilat,jlat,Nspin,Nspin,iorb,jorb)/=0d0) &
                         .AND. (ibdw(ibeta)==1) .AND. (ibdw(ialfa)==0)
                    !
                    if (Jcondition)then
                       call c(ibeta,mdw,k1,sg1)
                       call cdg(ialfa,k1,k2,sg2)
                       jdw = binary_search(Hs(2)%map,k2)
                       htmp = dmft_bath%item(ibath)%h(ilat,jlat,Nspin,Nspin,iorb,jorb)*sg1*sg2
                       !
                       call sp_insert_element(spH0dws(1),htmp,idw,jdw)
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
     !
     !
     !>H_hyb: hopping terms for a given spin (imp <--> bath)
     do ilat=1,Nlat
        do iorb=1,Norb
           do ibath=1,Nbath
              ialfa=getBathStride(ilat,iorb,ibath)
              is = imp_state_index(ilat,iorb) !imp site
              if( (diag_hybr(ilat,Nspin,iorb,ibath)/=0d0) &
                   .AND. (ibdw(is)==1) .AND. (ibdw(ialfa)==0) )then
                 call c(is,mdw,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 jdw=binary_search(Hs(2)%map,k2)
                 htmp=diag_hybr(ilat,Nspin,iorb,ibath)*sg1*sg2
                 !
                 call sp_insert_element(spH0dws(1),htmp,idw,jdw)
                 !
              endif
              if( (diag_hybr(ilat,Nspin,iorb,ibath)/=0d0) &
                   .AND. (ibdw(is)==0) .AND. (ibdw(ialfa)==1) )then
                 call c(ialfa,mdw,k1,sg1)
                 call cdg(is,k1,k2,sg2)
                 jdw=binary_search(Hs(2)%map,k2)
                 htmp=diag_hybr(ilat,Nspin,iorb,ibath)*sg1*sg2
                 !
                 call sp_insert_element(spH0dws(1),htmp,idw,jdw)
                 !
              endif
           enddo
        enddo
     enddo
     !
  enddo

