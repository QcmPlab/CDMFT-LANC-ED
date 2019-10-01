MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none
  private
  !
  public :: observables_impurity
  public :: local_energy_impurity


  logical,save                           :: iolegend=.true.
  real(8),dimension(:,:),allocatable     :: dens,dens_up,dens_dw
  real(8),dimension(:,:),allocatable     :: docc
  real(8),dimension(:,:),allocatable     :: magz
  real(8),dimension(:,:,:,:),allocatable :: sz2,n2
  real(8),dimension(:,:,:),allocatable   :: zimp,simp
  real(8),dimension(:),allocatable       :: s2tot
  real(8)                                :: Egs
  real(8)                                :: Ei
  real(8)                                :: integrationR
  !
  integer                                :: iorb,jorb,iorb1,jorb1
  integer                                :: ispin,jspin
  integer                                :: ilat,jlat
  integer                                :: ibath
  integer                                :: r,m,k,k1,k2
  integer                                :: iup,idw
  integer                                :: jup,jdw
  integer                                :: mup,mdw
  real(8)                                :: sgn,sgn1,sgn2,sg1,sg2
  real(8)                                :: gs_weight
  !
  real(8)                                :: peso
  real(8)                                :: norm
  !
  integer                                :: i,j,ii
  integer                                :: isector,jsector
  integer                                :: idim,idimUP,idimDW
  !
  real(8),dimension(:),pointer           :: state_cvec
  logical                                :: Jcondition
  !



contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_impurity()
    call lanc_observables()
  end subroutine observables_impurity



  subroutine local_energy_impurity()
    call lanc_local_energy()
  end subroutine local_energy_impurity




  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine lanc_observables()
    integer                             :: istate,Nud(2,Ns),iud(2),jud(2),is,js
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer,dimension(Ns)               :: IbUp,IbDw  ![Ns]
    real(8),dimension(Nlat,Norb)        :: nup,ndw,Sz,nt
    type(sector_map),dimension(2*Ns_Ud) :: HI
    !
    !LOCAL OBSERVABLES:
    allocate(dens(Nlat,Norb),dens_up(Nlat,Norb),dens_dw(Nlat,Norb))
    allocate(docc(Nlat,Norb),s2tot(Nlat))
    allocate(magz(Nlat,Norb),sz2(Nlat,Nlat,Norb,Norb),n2(Nlat,Nlat,Norb,Norb))
    allocate(simp(Nlat,Norb,Nspin),zimp(Nlat,Norb,Nspin))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       iDim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       if(MpiMaster)then
          call build_sector(isector,HI)
          do i = 1,iDim
             call state2indices(i,[iDimUps,iDimDws],Indices)
             do ii=1,Ns_Ud
                mup = HI(ii)%map(Indices(ii))
                mdw = HI(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                IbUp = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                IbDw = Bdecomp(mdw,Ns_Orb)
             enddo
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !Get operators:
             do ilat=1,Nlat
                do iorb=1,Norb
                   nup(ilat,iorb)= ibup(imp_state_index(ilat,iorb))
                   ndw(ilat,iorb)= ibdw(imp_state_index(ilat,iorb))
                   sz(ilat,iorb) = (nup(ilat,iorb) - ndw(ilat,iorb))/2d0
                   nt(ilat,iorb) =  nup(ilat,iorb) + ndw(ilat,iorb)
                enddo
             enddo
             !
             !Evaluate averages of observables:
             !>TODO:
             !add non-local averges like spin-spin, density-density etc...
             !<TODO
             do ilat=1,Nlat
                do iorb=1,Norb
                   dens(ilat,iorb)     = dens(ilat,iorb)      +  nt(ilat,iorb)*gs_weight
                   dens_up(ilat,iorb)  = dens_up(ilat,iorb)   +  nup(ilat,iorb)*gs_weight
                   dens_dw(ilat,iorb)  = dens_dw(ilat,iorb)   +  ndw(ilat,iorb)*gs_weight
                   docc(ilat,iorb)     = docc(ilat,iorb)      +  nup(ilat,iorb)*ndw(ilat,iorb)*gs_weight
                   magz(ilat,iorb)     = magz(ilat,iorb)      +  (nup(ilat,iorb)-ndw(ilat,iorb))*gs_weight
                enddo
                s2tot(ilat) = s2tot(ilat)  + (sum(sz(ilat,:)))**2*gs_weight
             enddo

             do ilat=1,Nlat
                do iorb=1,Norb
                   sz2(ilat,ilat,iorb,iorb) = sz2(ilat,ilat,iorb,iorb)  +  (sz(ilat,iorb)*sz(ilat,iorb))*gs_weight
                   n2(ilat,ilat,iorb,iorb)  = n2(ilat,ilat,iorb,iorb)   +  (nt(ilat,iorb)*nt(ilat,iorb))*gs_weight
                   do jlat=1,Nlat
                      do jorb=iorb+1,Norb
                         sz2(ilat,jlat,iorb,jorb) = sz2(ilat,jlat,iorb,jorb)  +  (sz(ilat,iorb)*sz(jlat,jorb))*gs_weight
                         sz2(ilat,jlat,jorb,iorb) = sz2(ilat,jlat,jorb,iorb)  +  (sz(ilat,jorb)*sz(jlat,iorb))*gs_weight
                         n2(ilat,jlat,iorb,jorb)  = n2(ilat,jlat,iorb,jorb)   +  (nt(ilat,iorb)*nt(jlat,jorb))*gs_weight
                         n2(ilat,jlat,jorb,iorb)  = n2(ilat,jlat,jorb,iorb)   +  (nt(jlat,jorb)*nt(ilat,iorb))*gs_weight
                      enddo
                   enddo
                enddo
             enddo
          enddo
          call delete_sector(isector,HI)
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    !
    !
    !
    if(MPIMASTER)then
       call get_szr
       if(iolegend)call write_legend
       call write_observables()
       !
       do ilat=1,Nlat
          write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(ilat,iorb),iorb=1,Norb),sum(dens(ilat,:))
       enddo
       write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(sum(dens(:,iorb))/Nlat,iorb=1,Norb),sum(dens)/Nlat
       !
       write(LOGfile,"(A,10f18.12,A)")"docc"//reg(ed_file_suffix)//"=",(sum(docc(:,iorb))/Nlat,iorb=1,Norb)
       if(Nspin==2)write(LOGfile,"(A,10f18.12,A)") "mag "//reg(ed_file_suffix)//"=",(sum(magz(:,iorb))/Nlat,iorb=1,Norb)
       !
       ed_dens_up=dens_up
       ed_dens_dw=dens_dw
       ed_dens   =dens
       ed_docc   =docc
    endif
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
    endif
#endif
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2,s2tot)
    deallocate(simp,zimp)
  end subroutine lanc_observables





  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine lanc_local_energy()
    integer                             :: istate,iud(2),jud(2),is,js
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    real(8),dimension(Ns)               :: Nup,Ndw
    type(sector_map),dimension(2*Ns_Ud) :: H
    !
    Egs     = state_list%emin
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    !
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       iDim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       !Master:
       if(MpiMaster)then
          call build_sector(isector,H)
          do i=1,iDim
             call state2indices(i,[iDimUps,iDimDws],Indices)
             do ii=1,Ns_Ud
                mup = H(ii)%map(Indices(ii))
                mdw = H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             Nup = Breorder(Nups)
             Ndw = Breorder(Ndws)
             !
             gs_weight=peso*abs(state_cvec(i))**2
             !
             !> H_Imp: Diagonal Elements, i.e. local part
             do ilat=1,Nlat
                do iorb=1,Norb
                   is = imp_state_index(ilat,iorb)
                   ed_Eknot = ed_Eknot + impHloc(ilat,ilat,1,1,iorb,iorb)*Nup(is)*gs_weight
                   ed_Eknot = ed_Eknot + impHloc(ilat,ilat,Nspin,Nspin,iorb,iorb)*Ndw(is)*gs_weight
                enddo
             enddo
             ! !> H_imp: Off-diagonal elements, i.e. non-local part. 
             iup = Indices(1)  ; idw = Indices(2)
             mup = H(1)%map(iup) ; mdw = H(2)%map(idw)
             do ilat=1,Nlat
                do jlat=1,Nlat
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = imp_state_index(ilat,iorb)
                         js = imp_state_index(jlat,jorb)
                         !UP
                         Jcondition = &
                              (impHloc(ilat,jlat,1,1,iorb,jorb)/=0d0) .AND. &
                              (Nup(js)==1) .AND. (Nup(is)==0)
                         if (Jcondition) then
                            call c(js,mup,k1,sg1)
                            call cdg(is,k1,k2,sg2)
                            jup = binary_search(H(1)%map,k2)
                            j   = jup + (idw-1)*iDimUp
                            ed_Eknot = ed_Eknot + &
                                 impHloc(ilat,jlat,1,1,iorb,jorb)*sg1*sg2*state_cvec(i)*(state_cvec(j))
                         endif
                         !
                         !DW
                         Jcondition = &
                              (impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)/=0d0) .AND. &
                              (ndw(js)==1) .AND. (ndw(is)==0)
                         if (Jcondition) then
                            call c(js,mdw,k1,sg1)
                            call cdg(is,k1,k2,sg2)
                            jdw = binary_search(H(2)%map,k2)
                            j   = iup + (jdw-1)*iDimUp
                            ed_Eknot = ed_Eknot + &
                                 impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)*sg1*sg2*state_cvec(i)*(state_cvec(j))
                         endif
                      enddo
                   enddo
                enddo
             enddo
             !
             !
             !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
             do ilat=1,Nlat
                do iorb=1,Norb
                   is = imp_state_index(ilat,iorb)
                   ed_Epot = ed_Epot + Uloc(iorb)*nup(is)*ndw(is)*gs_weight
                enddo
             enddo
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             if(Norb>1)then
                do ilat=1,Nlat
                   do jlat=ilat+1,Nlat
                      do iorb=1,Norb
                         do jorb=iorb+1,Norb
                            is = imp_state_index(ilat,iorb)
                            js = imp_state_index(jlat,jorb)
                            ed_Epot = ed_Epot + Ust*(nup(is)*ndw(js) + nup(js)*ndw(is))*gs_weight
                            ed_Dust = ed_Dust + (nup(is)*ndw(js) + nup(js)*ndw(is))*gs_weight
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do ilat=1,Nlat
                   do jlat=ilat+1,Nlat
                      do iorb=1,Norb
                         do jorb=iorb+1,Norb
                            is = imp_state_index(ilat,iorb)
                            js = imp_state_index(jlat,jorb)
                            ed_Epot = ed_Epot + (Ust-Jh)*(nup(is)*nup(js) + ndw(is)*ndw(js))*gs_weight
                            ed_Dund = ed_Dund + (nup(is)*nup(js) + ndw(is)*ndw(js))*gs_weight
                         enddo
                      enddo
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                !ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
                do ilat=1,Nlat
                   do iorb=1,Norb
                      is =imp_state_index(ilat,iorb)
                      ed_Ehartree=ed_Ehartree - 0.5d0*uloc(iorb)*(nup(is)+ndw(is))*gs_weight + 0.25d0*uloc(is)*gs_weight
                   enddo
                enddo
                if(Norb>1)then
                   do ilat=1,Nlat
                      do jlat=ilat+1,Nlat
                         do iorb=1,Norb
                            do jorb=iorb+1,Norb
                               is=imp_state_index(ilat,iorb)
                               js=imp_state_index(jlat,jorb)
                               ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(is)+ndw(is)+nup(js)+ndw(js))*gs_weight + 0.25d0*Ust*gs_weight
                               ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(is)+ndw(is)+nup(js)+ndw(js))*gs_weight + 0.25d0*(Ust-Jh)*gs_weight
                            enddo
                         enddo
                      enddo
                   enddo
                endif
             endif
          enddo
          call delete_sector(isector,H)         
       endif
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_Epot)
       call Bcast_MPI(MpiComm,ed_Eknot)
       call Bcast_MPI(MpiComm,ed_Ehartree)
       call Bcast_MPI(MpiComm,ed_Dust)
       call Bcast_MPI(MpiComm,ed_Dund)
    endif
#endif
    !
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose==3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
    endif
    !
    if(MPIMASTER)then
       call write_energy_info()
       call write_energy()
    endif
    !
    !
  end subroutine lanc_local_energy









  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: ilat,ispin,iorb
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3d0*pi/beta
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             simp(ilat,iorb,ispin) = dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1)) - &
                  wm1*(dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,2))-dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1)))/(wm2-wm1)
             zimp(ilat,iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1))/wm1 ))
          enddo
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(5*Norb+1))//"s2",&
         reg(txtfy(5*Norb+2))//"egs",&
         ((reg(txtfy(5*Norb+3+(ispin-1)*Norb+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         ((reg(txtfy((5+Nspin)*Norb+3+(ispin-1)*Norb+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend

  subroutine write_energy_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>"
    close(unit)
  end subroutine write_energy_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: iorb,jorb,ispin,ilat
    unit = free_unit()
    do ilat=1,Nlat
       open(unit,file="observables_all"//reg(ed_file_suffix)//"_site"//str(ilat,3)//".ed",position='append')
       write(unit,"(90(F15.9,1X))")&
            (dens(ilat,iorb),iorb=1,Norb),&
            (docc(ilat,iorb),iorb=1,Norb),&
            (dens_up(ilat,iorb),iorb=1,Norb),&
            (dens_dw(ilat,iorb),iorb=1,Norb),&
            (magz(ilat,iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((zimp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
    enddo
    !
    unit = free_unit()
    do ilat=1,Nlat
       open(unit,file="parameters_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
       close(unit)
    enddo
    !
    unit = free_unit()
    do ilat=1,Nlat
       open(unit,file="observables_last"//reg(ed_file_suffix)//"_site"//str(ilat,3)//".ed")
       write(unit,"(90(F15.9,1X))")&
            (dens(ilat,iorb),iorb=1,Norb),&
            (docc(ilat,iorb),iorb=1,Norb),&
            (dens_up(ilat,iorb),iorb=1,Norb),&
            (dens_dw(ilat,iorb),iorb=1,Norb),&
            (magz(ilat,iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((zimp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
    enddo

    unit = free_unit()
    open(unit,file="Sz_ij_ab_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"#I, J, a, b, Sz(I,J,a,b)"
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                write(unit,"(4I15,F15.9)")ilat,jlat,iorb,jorb,sz2(ilat,jlat,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
    close(unit)

    open(unit,file="N2_ij_ab_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"#I, J, a, b, N2(I,J,a,b)"
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                write(unit,"(4I15,F15.9)")ilat,jlat,iorb,jorb,n2(ilat,jlat,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
    close(unit)
  end subroutine write_observables

  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy



END MODULE ED_OBSERVABLES

















