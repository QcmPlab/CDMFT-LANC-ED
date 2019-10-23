  do idw=1,DimDw
     do iup=1,DimUp
        mup = Hs(1)%map(iup)
        ibup = bdecomp(mup,Ns)
        i    = iup + (idw-1)*dimUp
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
                    Jcondition = (impHloc(ilat,jlat,1,1,iorb,jorb)/=0d0)&
                         .AND.(ibup(js)==1).AND.(ibup(is)==0)
                    if (Jcondition) then
                       call c(js,mup,k1,sg1)
                       call cdg(is,k1,k2,sg2)
                       jup = binary_search(Hs(1)%map,k2)
                       j   = jup + (idw-1)*DimUp
                       htmp = impHloc(ilat,jlat,1,1,iorb,jorb)*sg1*sg2

                       !
                       Hv(i) = Hv(i) + htmp*vin(j)
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
                            (hbath_reconstructed(ilat,jlat,1,1,iorb,jorb,ibath)/=0d0) &
                            .AND. (ibup(ibeta)==1) .AND. (ibup(ialfa)==0)
                       !
                       if (Jcondition)then
                          call c(ibeta,mup,k1,sg1)
                          call cdg(ialfa,k1,k2,sg2)
                          jup = binary_search(Hs(1)%map,k2)
                          j   = jup + (idw-1)*DimUp
                          htmp = hbath_Reconstructed(ilat,jlat,1,1,iorb,jorb,ibath)*sg1*sg2
                          !
                          hv(i) = hv(i) + htmp*vin(j)
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
                 is = imp_state_index(ilat,iorb)
                 if( (diag_hybr(ilat,1,iorb,ibath)/=0d0) .AND. &
                      (ibup(is)==1) .AND. (ibup(ialfa)==0) )then              
                    call c(is,mup,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    jup = binary_search(Hs(1)%map,k2)
                    j   = jup + (idw-1)*DimUp
                    htmp = diag_hybr(ilat,1,iorb,ibath)*sg1*sg2
                    !
                    hv(i) = hv(i) + htmp*vin(j)
                    !
                 endif
                 if( (diag_hybr(ilat,1,iorb,ibath)/=0d0) .AND. &
                      (ibup(is)==0) .AND. (ibup(ialfa)==1) )then
                    call c(ialfa,mup,k1,sg1)
                    call cdg(is,k1,k2,sg2)
                    jup = binary_search(Hs(1)%map,k2)
                    j   = jup + (idw-1)*DimUp
                    htmp = diag_hybr(ilat,1,iorb,ibath)*sg1*sg2
                    !
                    hv(i) = hv(i) + htmp*vin(j)
                    !
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo





