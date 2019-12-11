!We build the transposed H_non_local here (symmetric)
!to comply with the MPI decomposition of the matrix.
!A better MPI handling might be necessary here... 
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
     ! SPIN-EXCHANGE (S-E) TERMS
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
