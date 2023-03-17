MODULE ED_CHI_SPIN
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  implicit none
  private


  public :: build_chi_spin_normal

  integer                          :: istate
  integer                          :: isector
  integer                          :: idim,idimUP,idimDW
  integer                          :: jdim,jdimUP,jdimDW
  complex(8),allocatable           :: vvinit(:),vvloc(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: ialfa,ibeta
  integer                          :: jalfa,jbeta
  integer                          :: r
  integer                          :: i,iup,idw
  integer                          :: j,jup,jdw  
  integer                          :: m,mup,mdw
  integer                          :: ipos,jpos
  real(8)                          :: sgn,norm2
  complex(8),dimension(:),allocatable :: state_cvec
  real(8)                          :: state_e
  integer                          :: Nitermax,Nlanc,vecDim




contains


  !+------------------------------------------------------------------+
  !                            SPIN
  !PURPOSE  : Evaluate the Spin susceptibility \Chi_spin for a 
  ! single orbital: \chi = <S_a(\tau)S_a(0)>
  ! note: as S_a is hermitian particle and holes contributions
  ! are identical so work out only one lanczos tridiag. work out the 
  ! reduction for both values of isign in the same call.
  !+------------------------------------------------------------------+
  subroutine build_chi_spin_normal()
  real(8) :: chan4
  integer :: isite,jsite,iorb,jorb
  !
  if(ed_gf_symmetric)then
    chan4=0.d0
  else
    chan4=1.d0
  endif
    write(LOGfile,"(A)")"Get impurity spin Chi:"
    do isite=1,Nlat
      do iorb=1,Norb
         write(LOGfile,"(A)")"Get Chi_spin_l"//reg(txtfy(iorb))
         call lanc_ed_build_spinChi_diag(isite,iorb)
      enddo
    enddo
    !
    do isite=1,Nlat
      do iorb=1,Norb
        do jsite=1,Nlat
          do jorb=1,Norb
            if(isite==jsite .and. iorb==jorb)cycle
            if(ed_gf_symmetric)then
              call lanc_ed_build_spinChi_mix_chan2(isite,jsite,iorb,jorb)
            else
              call lanc_ed_build_spinChi_mix_chan4(isite,jsite,iorb,jorb)
            endif
          end do
        enddo
      enddo
    enddo
       !
       !
    do isite=1,Nlat
      do jsite=1,Nlat
        do iorb=1,Norb
          do jorb=1,Norb
             spinChi_w(isite,jsite,iorb,jorb,:)   = 0.5d0*(spinChi_w(isite,jsite,iorb,jorb,:) & 
                                        - (one-chan4*xi)*spinChi_w(isite,isite,iorb,iorb,:) - (one-chan4*xi)*spinChi_w(jsite,jsite,jorb,jorb,:))
             spinChi_tau(isite,jsite,iorb,jorb,:) = 0.5d0*(spinChi_tau(isite,jsite,iorb,jorb,:) & 
                                        - (one-chan4*xi)*spinChi_tau(isite,isite,iorb,iorb,:) - (one-chan4*xi)*spinChi_tau(jsite,jsite,jorb,jorb,:))
             spinChi_iv(isite,jsite,iorb,jorb,:)  = 0.5d0*(spinChi_iv(isite,jsite,iorb,jorb,:) & 
                                        - (one-chan4*xi)*spinChi_iv(isite,isite,iorb,iorb,:) - (one-chan4*xi)*spinChi_iv(jsite,jsite,jorb,jorb,:))
          enddo
        enddo
      enddo
    enddo
    !
  end subroutine build_chi_spin_normal
  
  


  !################################################################
  !################################################################
  !################################################################
  !################################################################

  subroutine lanc_ed_build_spinChi_diag(isite,iorb)
    integer                     :: isite,iorb,ipos
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer,dimension(2)        :: Iud
    type(sector_map)            :: HI(2*Ns_Ud)
    !
    ialfa = 1
    ipos = imp_state_index(isite,iorb)
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_cvector(MpiComm,state_list,istate,state_cvec)
       else
          call es_return_cvector(state_list,istate,state_cvec)
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec)
#endif
       idim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,HI)
       !
       if(MpiMaster)then
          if(ed_verbose>=3)write(LOGfile,"(A20,I6)")&
               'Apply Sz',isector
          allocate(vvinit(idim)) ; vvinit=zero
          do i=1,idim
            sgn=0d0
            call state2indices(i,[iDimUps,iDimDws],Indices)
            iud(1)   = HI(ialfa)%map(Indices(ialfa))
            iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
            nud(1,:) = Bdecomp(iud(1),Ns_Orb)
            nud(2,:) = Bdecomp(iud(2),Ns_Orb)
            !
            sgn = dble(nud(1,ipos))-dble(nud(2,ipos))
            sgn = sgn/2d0
            vvinit(i) = sgn*state_cvec(i)
          enddo
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
       alfa_=0.d0
       beta_=0.d0
       call build_Hv_sector(isector)
#ifdef _MPI
       if(MpiStatus)then
         if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,norm2)
         vecDim = vecDim_Hv_sector(isector)
         allocate(vvloc(vecDim))
         if(MpiComm /= MPI_COMM_NULL) call scatter_vector_MPI(MpiComm,vvinit,vvloc)
         call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alfa_,beta_)
       else
         call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
#endif
       call delete_Hv_sector()
       call add_to_lanczos_spinChi(one*norm2,state_e,alfa_,beta_,isite,isite,iorb,iorb)
          !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)          
       if(allocated(vvloc))deallocate(vvloc)
       if(allocated(state_cvec))deallocate(state_cvec)
       call delete_sector(isector,HI)
    enddo
    return
  end subroutine lanc_ed_build_spinChi_diag




  !################################################################




  subroutine lanc_ed_build_spinChi_mix_chan2(isite,jsite,iorb,jorb)
    integer                     :: isite,iorb,ipos,jsite,jorb,jpos
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer,dimension(2)        :: Iud
    type(sector_map)            :: HI(2*Ns_Ud)
    real(8)                     :: Siorb,Sjorb
    !
    !  
    ialfa = 1
    jalfa = 1
    ipos = imp_state_index(isite,iorb)
    jpos = imp_state_index(jsite,jorb)
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_cvector(MpiComm,state_list,istate,state_cvec)
       else
          call es_return_cvector(state_list,istate,state_cvec)
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec)
#endif
       idim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,HI)
       !
       !EVALUATE (Sz_jorb + Sz_iorb)|gs> = Sz_jorb|gs> + Sz_iorb|gs>
       if(MpiMaster)then
          if(ed_verbose>=3)write(LOGfile,"(A20,I6)")&
               'From sector',isector
          if(ed_verbose>=3)write(LOGfile,"(A20,I15)")'Apply (Sz_a + Sz_b)',isector
          allocate(vvinit(idim)) ; vvinit=zero
          do i=1,idim
             sgn=0d0
             call state2indices(i,[iDimUps,iDimDws],Indices)
             iud(1)   = HI(ialfa)%map(Indices(ialfa))
             iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
             nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             Siorb = (dble(nud(1,ipos))-dble(nud(2,ipos)))/2d0
             Sjorb = (dble(nud(1,jpos))-dble(nud(2,jpos)))/2d0
             sgn       = Siorb + Sjorb
             vvinit(i) = sgn*state_cvec(i)
          enddo
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
       alfa_=0.d0
       beta_=0.d0
       call build_Hv_sector(isector)
#ifdef _MPI
       if(MpiStatus)then
         if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,norm2)
         vecDim = vecDim_Hv_sector(isector)
         allocate(vvloc(vecDim))
         if(MpiComm /= MPI_COMM_NULL) call scatter_vector_MPI(MpiComm,vvinit,vvloc)
         call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alfa_,beta_)
       else
         call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
#endif
       call delete_Hv_sector()
       call add_to_lanczos_spinChi(one*norm2,state_e,alfa_,beta_,isite,jsite,iorb,jorb)
          !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)          
       if(allocated(vvloc))deallocate(vvloc)
       if(allocated(state_cvec))deallocate(state_cvec)
       call delete_sector(isector,HI)
    enddo
    return
  end subroutine lanc_ed_build_spinChi_mix_chan2

!################################################################

  subroutine lanc_ed_build_spinChi_mix_chan4(isite,jsite,iorb,jorb)
    integer                     :: isite,iorb,ipos,jsite,jorb,jpos
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer,dimension(2)        :: Iud
    type(sector_map)            :: HI(2*Ns_Ud)
    real(8)                     :: Siorb,Sjorb
    !
    !  
    ialfa = 1
    jalfa = 1
    ipos = imp_state_index(isite,iorb)
    jpos = imp_state_index(jsite,jorb)
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_cvector(MpiComm,state_list,istate,state_cvec)
       else
          call es_return_cvector(state_list,istate,state_cvec)
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec)
#endif
       idim  = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       call build_sector(isector,HI)
       !
       !EVALUATE (Sz_jorb + Sz_iorb)|gs> = Sz_jorb|gs> + i*Sz_iorb|gs>
       if(MpiMaster)then
          if(ed_verbose>=3)write(LOGfile,"(A20,I6)")&
               'From sector',isector
          if(ed_verbose>=3)write(LOGfile,"(A20,I15)")'Apply (Sz_a + Sz_b)',isector
          allocate(vvinit(idim)) ; vvinit=zero
          do i=1,idim
             sgn=0d0
             call state2indices(i,[iDimUps,iDimDws],Indices)
             iud(1)   = HI(ialfa)%map(Indices(ialfa))
             iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
             nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             Siorb = (dble(nud(1,ipos))-dble(nud(2,ipos)))/2d0
             Sjorb = (dble(nud(1,jpos))-dble(nud(2,jpos)))/2d0
             sgn       = Siorb + Sjorb
             vvinit(i) = sgn*state_cvec(i)
          enddo
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
       alfa_=0.d0
       beta_=0.d0
       call build_Hv_sector(isector)
#ifdef _MPI
       if(MpiStatus)then
         if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,norm2)
         vecDim = vecDim_Hv_sector(isector)
         allocate(vvloc(vecDim))
         if(MpiComm /= MPI_COMM_NULL) call scatter_vector_MPI(MpiComm,vvinit,vvloc)
         call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alfa_,beta_)
       else
         call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
#endif
       call delete_Hv_sector()
       call add_to_lanczos_spinChi(one*norm2,state_e,alfa_,beta_,isite,jsite,iorb,jorb)
          !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)          
       if(allocated(vvloc))deallocate(vvloc)
       
       !EVALUATE (Sz_jorb + i*Sz_iorb)|gs> = Sz_jorb|gs> + i*Sz_iorb|gs>
       if(MpiMaster)then
          if(ed_verbose>=3)write(LOGfile,"(A20,I6)")&
               'From sector',isector
          if(ed_verbose>=3)write(LOGfile,"(A20,I15)")'Apply (Sz_a + Sz_b)',isector
          allocate(vvinit(idim)) ; vvinit=zero
          do i=1,idim
             sgn=0d0
             call state2indices(i,[iDimUps,iDimDws],Indices)
             iud(1)   = HI(ialfa)%map(Indices(ialfa))
             iud(2)   = HI(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
             nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             nud(2,:) = Bdecomp(iud(2),Ns_Orb)
             Siorb = (dble(nud(1,ipos))-dble(nud(2,ipos)))/2d0
             Sjorb = (dble(nud(1,jpos))-dble(nud(2,jpos)))/2d0
             sgn       = Siorb + xi*Sjorb
             vvinit(i) = sgn*state_cvec(i)
          enddo
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       nlanc=min(idim,lanc_nGFiter)
       allocate(alfa_(nlanc),beta_(nlanc))
       alfa_=0.d0
       beta_=0.d0
       call build_Hv_sector(isector)
#ifdef _MPI
       if(MpiStatus)then
         if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,norm2)
         vecDim = vecDim_Hv_sector(isector)
         allocate(vvloc(vecDim))
         if(MpiComm /= MPI_COMM_NULL) call scatter_vector_MPI(MpiComm,vvinit,vvloc)
         call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alfa_,beta_)
       else
         call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alfa_,beta_)
#endif
       call delete_Hv_sector()
       call add_to_lanczos_spinChi(one*norm2,state_e,alfa_,beta_,isite,jsite,iorb,jorb)
          !
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)          
       if(allocated(vvloc))deallocate(vvloc)
       if(allocated(state_cvec))deallocate(state_cvec)
       call delete_sector(isector,HI)
    enddo
    return
  end subroutine lanc_ed_build_spinChi_mix_chan4



  !################################################################


  subroutine add_to_lanczos_spinChi(vnorm2,Ei,alanc,blanc,isite,jsite,iorb,jorb)
    real(8)                                    :: Ei,Ej,Egs,de
    complex(8)                                 :: vnorm2,pesoF,pesoAB,pesoBZ,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: iorb,jorb,isite,jsite
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_spinChi: add-up to GF"
#endif
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
#ifdef _MPI
    if(MpiStatus)then
       if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,alanc)
       if(MpiComm /= MPI_COMM_NULL)call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)

    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! This ensures that the correct null contribution is obtained.
       ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! we do not include the contribution (because we are in the situation described above).
       ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       if(beta*dE > 1d-3)spinChi_iv(isite,jsite,iorb,jorb,0)=spinChi_iv(isite,jsite,iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          spinChi_iv(isite,jsite,iorb,jorb,i)=spinChi_iv(isite,jsite,iorb,jorb,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          spinChi_tau(isite,jsite,iorb,jorb,i)=spinChi_tau(isite,jsite,iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          spinChi_w(isite,jsite,iorb,jorb,i)=spinChi_w(isite,jsite,iorb,jorb,i) - peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
  end subroutine add_to_lanczos_spinChi


END MODULE ED_CHI_SPIN

