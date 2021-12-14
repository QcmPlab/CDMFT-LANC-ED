MODULE ED_GREENS_FUNCTIONS
  USE ED_GF_SHARED
  USE ED_GF_NORMAL
  !USE ED_GF_CHISPIN
  !
  implicit none
  private 

  public :: buildGf_impurity
  !public :: buildChi_impurity


  real(8),dimension(:,:,:),allocatable            :: zimp,simp


contains



  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildGF_impurity()
    !
    call allocate_grids
    !
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    call build_sigma_normal()
    !
    if(MPIMASTER)then
       if(ed_print_Sigma)call ed_print_impSigma()
       if(ed_print_G)call ed_print_impG()
       if(ed_print_G0)call ed_print_impG0()
    endif
    !
    allocate(simp(Nlat,Norb,Nspin),zimp(Nlat,Norb,Nspin))
    call get_szr
    call write_szr()
    deallocate(simp,zimp)

    !
    call deallocate_grids
    !
  end subroutine buildgf_impurity








  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  !subroutine buildChi_impurity()
  !!
  !call allocate_grids
  !!
  !!
  !!BUILD SPIN SUSCEPTIBILITY
  !spinChi_tau=zero
  !spinChi_w=zero
  !spinChi_iv=zero
  !call build_chi_spin()
  !!
  !!
  !! !BUILD CHARGE SUSCEPTIBILITY
  !! densChi_tau=zero
  !! densChi_w=zero
  !! densChi_iv=zero
  !! densChi_mix_tau=zero
  !! densChi_mix_w=zero
  !! densChi_mix_iv=zero
  !! densChi_tot_tau=zero
  !! densChi_tot_w=zero
  !! densChi_tot_iv=zero
  !! call build_chi_dens()
  !!
  !!
  !! !BUILD PAIR SUSCEPTIBILITY
  !! pairChi_tau=zero
  !! pairChi_w=zero
  !! pairChi_iv=zero
  !! call build_chi_pair()
  !!
  !!
  !!PRINTING:
  !if(MPIMASTER)call ed_print_impChi()
  !!
  !!
  !call deallocate_grids
  !!
  !end subroutine buildChi_impurity




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
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_szr()
    integer :: unit
    integer :: iorb,jorb,ispin,ilat
    !
    open(free_unit(unit),file="zeta_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         ((reg(txtfy(iorb+(ispin-1)*Norb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    open(free_unit(unit),file="sig_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         ((reg(txtfy(iorb+(ispin-1)*Norb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin) 
    close(unit)
    !
    do ilat=1,Nlat
       open(free_unit(unit),file="zeta_all"//reg(ed_file_suffix)//"_site"//str(ilat,3)//".ed",position='append')
       write(unit,"(90(F15.9,1X))")&
            ((zimp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
       open(free_unit(unit),file="zeta_last"//reg(ed_file_suffix)//"_site"//str(ilat,3)//".ed")
       write(unit,"(90(F15.9,1X))")&
            ((zimp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
       !
       open(free_unit(unit),file="sig_all"//reg(ed_file_suffix)//"_site"//str(ilat,3)//".ed",position='append')
       write(unit,"(90(F15.9,1X))")&
            ((simp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin)        
       close(unit)
       open(free_unit(unit),file="sig_last"//reg(ed_file_suffix)//"_site"//str(ilat,3)//".ed")
       write(unit,"(90(F15.9,1X))")&
            ((simp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin)    
       close(unit)
    enddo
    !
  end subroutine write_szr





end MODULE ED_GREENS_FUNCTIONS
