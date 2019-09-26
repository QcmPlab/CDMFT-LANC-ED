MODULE ED_GREENS_FUNCTIONS
  USE ED_GF_SHARED
  USE ED_GF_NORMAL
  !USE ED_GF_CHISPIN
  ! USE ED_CHI_DENS
  ! USE ED_CHI_PAIR
  !
  implicit none
  private 

  public :: buildGf_impurity
  !public :: buildChi_impurity

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




end MODULE ED_GREENS_FUNCTIONS
