  do i=MpiIstart,MpiIend
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hs(1)%map(iup)
     mdw = Hs(2)%map(idw)
     !
     ibup = bdecomp(mup,Ns)
     ibdw = bdecomp(mdw,Ns)
     !
     !
     do ilat=1,Nlat
        do iorb=1,Norb
           nup(ilat,iorb)=dble(ibup(imp_state_index(ilat,iorb)))
           ndw(ilat,iorb)=dble(ibdw(imp_state_index(ilat,iorb)))
        enddo
     enddo
     !
     !
     !SPIN-EXCHANGE IN IMPHLOC
     !
     !
     if(Jhflag.AND.Nspin>1)then
        do ilat=1,Nlat
          do jlat=1,Nlat
            do iorb=1,Norb
              do jorb=1,Norb
                is = imp_state_index(ilat,iorb) 
                js = imp_state_index(jlat,jorb)
                !UP->DOWN
                Jcondition=(&
                  (is/=js).AND.&
                  (nup(ilat,jorb)==1).AND.&
                  (ndw(jlat,jorb)==0))
                if(Jcondition)then
                  call c(is,mup,k1,sg1)  !DESTROY UP
                  jup=binary_search(Hs(1)%map,k1)
                  call cdg(js,k1,k2,sg2) !CREATE  DW
                  jdw=binary_search(Hs(2)%map,k2)
                  htmp = impHloc(ilat,jlat,1,2,iorb,jorb)*sg1*sg2
                  j = jup + (jdw-1)*DimUp
                  !
                  ! if(j==0)cycle
                  select case(MpiStatus)
                  case (.true.)
                    call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                  case (.false.)
                    call sp_insert_element(spH0nd,htmp,i,j)
                  end select
                  !
                endif
                !DOWN->UP
                Jcondition=(&
                  (is/=js).AND.&
                  (ndw(ilat,jorb)==1).AND.&
                  (nup(jlat,jorb)==0))
                if(Jcondition)then
                  call c(is,mdw,k1,sg1)  !DESTROY DW
                  jdw=binary_search(Hs(2)%map,k1)
                  call cdg(js,k1,k2,sg2) !CREATE  UP
                  jup=binary_search(Hs(1)%map,k2)
                  htmp = impHloc(ilat,jlat,2,1,iorb,jorb)*sg1*sg2
                  j = jup + (jdw-1)*DimUp
                  !
                  ! if(j==0)cycle
                  select case(MpiStatus)
                  case (.true.)
                    call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                  case (.false.)
                    call sp_insert_element(spH0nd,htmp,i,j)
                  end select
                  !
                endif
              enddo
            enddo
          enddo
        enddo
      endif
      !
      !
      !SPIN-EXCHANGE IN BATH
      !
      !
      if(Jhflag.AND.Nspin>1)then
        do ibath=1,Nbath
          do ilat=1,Nlat
            do jlat=1,Nlat
              do iorb=1,Norb
                do jorb=1,Norb
                  is = getbathstride(ilat,iorb,ibath) 
                  js = getBathStride(jlat,jorb,ibath)
                  !UP->DOWN
                  Jcondition=(&
                    (is/=js).AND.&
                    (ibup(is)==1).AND.&
                    (ibdw(js)==0))
                  if(Jcondition)then
                    call c(is,mup,k1,sg1)  !DESTROY UP
                    jup=binary_search(Hs(1)%map,k1)
                    call cdg(js,k1,k2,sg2) !CREATE  DW
                    jdw=binary_search(Hs(2)%map,k2)
                    htmp = dmft_bath%item(ibath)%h(ilat,jlat,1,2,iorb,jorb)*sg1*sg2
                    j = jup + (jdw-1)*DimUp
                    !
                    ! if(j==0)cycle
                    select case(MpiStatus)
                    case (.true.)
                      call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                    case (.false.)
                      call sp_insert_element(spH0nd,htmp,i,j)
                    end select
                    !
                  endif
                  !DOWN->UP
                  Jcondition=(&
                    (is/=js).AND.&
                    (ibdw(is)==1).AND.&
                    (ibup(js)==0))
                  if(Jcondition)then
                    call c(is,mup,k1,sg1)  !DESTROY DW
                    jdw=binary_search(Hs(2)%map,k1)
                    call cdg(js,k1,k2,sg2) !CREATE  UP
                    jup=binary_search(Hs(1)%map,k2)
                    htmp = dmft_bath%item(ibath)%h(ilat,jlat,2,1,iorb,jorb)*sg1*sg2
                    j = jup + (jdw-1)*DimUp
                    !
                    ! if(j==0)cycle
                    select case(MpiStatus)
                    case (.true.)
                      call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                    case (.false.)
                      call sp_insert_element(spH0nd,htmp,i,j)
                    end select
                    !
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      ! SPIN-EXCHANGE (S-E) TERMS
      !    S-E: J c^+_a_up c^+_b_dw c_a_dw c_b_up
      !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
      !
      !    S-E: J  [c^+_b_dw c_a_dw] [c^+_a_up c_b_up]
      !    S-E: J  [c^+_{jorb} c_{iorb}]_dw [c^+_iorb c_jorb]_up
      if(Jhflag.AND.Jx/=0d0)then
        do ilat=1,Nlat
          do iorb=1,Norb
            do jorb=1,Norb
              is = imp_state_index(ilat,iorb) 
              js = imp_state_index(ilat,jorb)
              Jcondition=(&
                (is/=js).AND.&
                (nup(ilat,jorb)==1).AND.&
                (ndw(ilat,iorb)==1).AND.&
                (ndw(ilat,jorb)==0).AND.&
                (nup(ilat,iorb)==0))
              if(Jcondition)then
                call c(is,mdw,k1,sg1)  !DW
                call cdg(js,k1,k2,sg2) !DW
                jdw=binary_search(Hs(2)%map,k2)
                call c(js,mup,k3,sg3)  !UP
                call cdg(is,k3,k4,sg4) !UP
                jup=binary_search(Hs(1)%map,k4)
                htmp = Jx*sg1*sg2*sg3*sg4
                j = jup + (jdw-1)*DimUp
                !
                ! if(j==0)cycle
                select case(MpiStatus)
                case (.true.)
                  call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                case (.false.)
                  call sp_insert_element(spH0nd,htmp,i,j)
                end select
                !
              endif
            enddo
          enddo
        enddo
      endif
      !
      ! PAIR-HOPPING (P-H) TERMS
      !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
      !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
      if(Jhflag.AND.Jp/=0d0)then
        do ilat=1,Nlat
          do iorb=1,Norb
            do jorb=1,Norb
              is = imp_state_index(ilat,iorb) 
              js = imp_state_index(ilat,jorb)
              Jcondition=(&
                (nup(ilat,jorb)==1).AND.&
                (ndw(ilat,jorb)==1).AND.&
                (ndw(ilat,iorb)==0).AND.&
                (nup(ilat,iorb)==0))
              if(Jcondition)then
                call c(js,mdw,k1,sg1)       !c_jorb_dw
                call cdg(is,k1,k2,sg2)      !c^+_iorb_dw
                jdw = binary_search(Hs(2)%map,k2)
                call c(js,mup,k3,sg3)       !c_jorb_up
                call cdg(is,k3,k4,sg4)      !c^+_iorb_up
                jup = binary_search(Hs(1)%map,k4)
                htmp = Jp*sg1*sg2*sg3*sg4
                j = jup + (jdw-1)*dimup
                !
                ! if(j==0)cycle
                select case(MpiStatus)
                case (.true.)
                  call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
                case (.false.)
                  call sp_insert_element(spH0nd,htmp,i,j)
                end select
                !
              endif
            enddo
          enddo
        enddo
      endif
      !
    enddo
