MODULE ED_IO
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  implicit none
  private

  !Retrieve self-energy through routines:
  interface ed_get_sigma_matsubara
     module procedure ed_get_sigma_matsubara_1
     module procedure ed_get_sigma_matsubara_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure ed_get_sigma_matsubara_lattice_1
#endif
     !     module procedure ed_get_sigma_matsubara_lattice_2
  end interface ed_get_sigma_matsubara

  interface ed_get_sigma_realaxis
     module procedure ed_get_sigma_real_1
     module procedure ed_get_sigma_real_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure ed_get_sigma_real_lattice_1
#endif
     !     module procedure ed_get_sigma_real_lattice_2
  end interface ed_get_sigma_realaxis

  !Retrieve imp GF through routines.
  interface ed_get_gimp_matsubara
     module procedure ed_get_gimp_matsubara_1
     module procedure ed_get_gimp_matsubara_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure ed_get_gimp_matsubara_lattice_1
#endif
     !     module procedure ed_get_gimp_matsubara_lattice_2
  end interface ed_get_gimp_matsubara


  interface ed_get_gimp_realaxis
     module procedure ed_get_gimp_real_1
     module procedure ed_get_gimp_real_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     procedure ed_get_gimp_real_lattice_1
#endif
     !module procedure ed_get_gimp_real_lattice_2
  end interface ed_get_gimp_realaxis




  !Retrieve imp GF_0 (G0_and) through routines.
  interface ed_get_g0imp_matsubara
     module procedure ed_get_g0imp_matsubara_1
     module procedure ed_get_g0imp_matsubara_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure ed_get_g0imp_matsubara_lattice_1
#endif
     !     module procedure ed_get_g0imp_matsubara_lattice_2
  end interface ed_get_g0imp_matsubara


  interface ed_get_g0imp_realaxis
     module procedure ed_get_g0imp_real_1
     module procedure ed_get_g0imp_real_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure ed_get_g0imp_real_lattice_1
#endif
     !     module procedure ed_get_g0imp_real_lattice_2
  end interface ed_get_g0imp_realaxis



  interface ed_get_delta_matsubara
     module procedure delta_bath_main_
  end interface ed_get_delta_matsubara

  interface ed_get_delta_realaxis
     module procedure delta_bath_main_ !mats is the same as real
  end interface ed_get_delta_realaxis


  interface ed_get_g0and_matsubara
     module procedure g0and_bath_main_
  end interface ed_get_g0and_matsubara

  interface ed_get_g0and_realaxis
     module procedure g0and_bath_main_
  end interface ed_get_g0and_realaxis



  interface ed_get_invg0and_matsubara
     module procedure invg0_bath_main_
  end interface ed_get_invg0and_matsubara

  interface ed_get_invg0and_realaxis
     module procedure invg0_bath_main_
  end interface ed_get_invg0and_realaxis


  !Retrieve static common observables  
  interface ed_get_dens
     module procedure ed_get_dens_1
     module procedure ed_get_dens_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure ed_get_dens_lattice_1
#endif
     !     module procedure ed_get_dens_lattice_2
  end interface ed_get_dens

  interface ed_get_mag
     module procedure ed_get_mag_1
     module procedure ed_get_mag_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure ed_get_mag_lattice_1
#endif
     !     module procedure ed_get_mag_lattice_2
  end interface ed_get_mag

  interface ed_get_docc
     module procedure ed_get_docc_1
     module procedure ed_get_docc_2
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure ed_get_docc_lattice_1
#endif
     !     module procedure ed_get_docc_lattice_2
  end interface ed_get_docc

  interface ed_get_epot
     module procedure :: ed_get_epot_
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_get_epot_lattice
#endif
  end interface ed_get_epot

  interface ed_get_eint
     module procedure :: ed_get_eint_
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_get_eint_lattice
#endif
  end interface ed_get_eint

  interface ed_get_ehartree
     module procedure :: ed_get_ehartree_
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_get_ehartree_lattice
#endif
  end interface ed_get_ehartree

  interface ed_get_eknot
     module procedure :: ed_get_eknot_
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_get_eknot_lattice
#endif
  end interface ed_get_eknot


  interface ed_get_dust
     module procedure :: ed_get_dust_
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_get_dust_lattice
#endif
  end interface ed_get_dust

  interface ed_get_dund
     module procedure :: ed_get_dund_
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_get_dund_lattice
#endif
  end interface ed_get_dund

  interface ed_get_dse
     module procedure :: ed_get_dse_
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_get_dse_lattice
#endif
  end interface ed_get_dse

  interface ed_get_dph
     module procedure :: ed_get_dph_
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_get_dph_lattice
#endif
  end interface ed_get_dph



  interface ed_get_cluster_dm
      module procedure :: ed_get_cluster_density_matrix_single
#if __GFORTRAN__ &&  __GNUC__ > 8     
      module procedure :: ed_get_cluster_density_matrix_lattice
#endif
   end interface ed_get_cluster_dm

   interface ed_get_reduced_dm
      module procedure :: ed_get_reduced_density_matrix_single
#if __GFORTRAN__ &&  __GNUC__ > 8     
      module procedure :: ed_get_reduced_density_matrix_lattice
#endif
   end interface ed_get_reduced_dm

   interface ed_get_sp_dm
      module procedure :: ed_get_single_particle_density_matrix_single
#if __GFORTRAN__ &&  __GNUC__ > 8     
      module procedure :: ed_get_single_particle_density_matrix_lattice
#endif
   end interface ed_get_sp_dm



  interface ed_gf_cluster
     module procedure :: ed_gf_cluster_scalar
     module procedure :: ed_gf_cluster_array
  end interface ed_gf_cluster

  interface ed_read_impsigma
     module procedure :: ed_read_impsigma_single
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_read_impsigma_lattice
#endif
  end interface ed_read_impsigma

  interface ed_read_impG
     module procedure :: ed_read_impG_single
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: ed_read_impG_lattice
#endif
  end interface ed_read_impG

  public :: ed_get_sigma_matsubara
  public :: ed_get_sigma_realaxis

  public :: ed_get_gimp_matsubara
  public :: ed_get_gimp_realaxis

  public :: ed_get_g0imp_matsubara
  public :: ed_get_g0imp_realaxis

  public :: ed_get_delta_matsubara
  public :: ed_get_delta_realaxis

  public :: ed_get_g0and_matsubara
  public :: ed_get_g0and_realaxis

  public :: ed_get_invG0and_matsubara
  public :: ed_get_invG0and_realaxis


  public :: ed_gf_cluster

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc

  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot

  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph

  public :: ed_get_cluster_dm
  public :: ed_get_reduced_dm
  public :: ed_get_sp_dm

  public :: ed_read_impSigma
  public :: ed_read_impG

  !****************************************************************************************!
  !****************************************************************************************!

  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_dm
  public :: ed_print_impChi


  !****************************************************************************************!
  !****************************************************************************************!




  character(len=64)                :: suffix




contains



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity functions
  !+--------------------------------------------------------------------------+!
  include "ED_IO/get_sigma_matsubara.f90"
  include "ED_IO/get_sigma_realaxis.f90"
  include "ED_IO/get_gimp_matsubara.f90"
  include "ED_IO/get_gimp_realaxis.f90"
  include "ED_IO/get_g0imp_matsubara.f90"
  include "ED_IO/get_g0imp_realaxis.f90"
  include "ED_IO/get_Gand_all.f90"
  include "ED_IO/gf_cluster.f90"

#if __GFORTRAN__ &&  __GNUC__ > 8    
  include "ED_IO/lattice/get_sigma_matsubara.f90"
  include "ED_IO/lattice/get_sigma_realaxis.f90"
  include "ED_IO/lattice/get_gimp_matsubara.f90"
  include "ED_IO/lattice/get_gimp_realaxis.f90"
  include "ED_IO/lattice/get_g0imp_matsubara.f90"
  include "ED_IO/lattice/get_g0imp_realaxis.f90"
#endif


  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+--------------------------------------------------------------------------+!
  include "ED_IO/get_dens.f90"
  include "ED_IO/get_mag.f90"
  include "ED_IO/get_docc.f90"
  include "ED_IO/get_eimp.f90"
  include "ED_IO/get_doubles.f90"
  include "ED_IO/get_cluster_dm.f90"
  include "ED_IO/get_reduced_dm.f90"
  include "ED_IO/get_sp_dm.f90"

#if __GFORTRAN__ &&  __GNUC__ > 8    
  include "ED_IO/lattice/get_dens.f90"
  include "ED_IO/lattice/get_mag.f90"
  include "ED_IO/lattice/get_docc.f90"
  include "ED_IO/lattice/get_eimp.f90"
  include "ED_IO/lattice/get_doubles.f90"
  include "ED_IO/lattice/get_cluster_dm.f90"
  include "ED_IO/lattice/get_reduced_dm.f90"
  include "ED_IO/lattice/get_sp_dm.f90"
#endif

  !+------------------------------------------------------------------+
  !                         PRINT SIGMA:
  !+------------------------------------------------------------------+  
  subroutine ed_print_impSigma
    character(len=64) :: suffix
    integer           :: ilat,jlat,iorb,jorb,ispin
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !Print the impurity Sigma:
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                do ispin=1,Nspin
                   suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                   call splot("impSigma"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impSmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                   call splot("impSigma"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impSreal(ilat,jlat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
  end subroutine ed_print_impSigma




  !+------------------------------------------------------------------+
  !                         PRINT G
  !+------------------------------------------------------------------+  
  subroutine ed_print_impG
    character(len=64) :: suffix
    integer           :: ilat,jlat,iorb,ispin,jorb
    !
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb  
             do jorb=1,Norb
                do ispin=1,Nspin
                   suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                   call splot("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impGmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                   call splot("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impGreal(ilat,jlat,ispin,ispin,iorb,jorb,:))                
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine ed_print_impG



  !+------------------------------------------------------------------+
  !                         PRINT G0
  !+------------------------------------------------------------------+  
  subroutine ed_print_impG0
    character(len=64) :: suffix
    integer           :: ilat,jlat,iorb,ispin,jorb
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                do ispin=1,Nspin
                   suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
                   call splot("impG0"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impG0mats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                   call splot("impG0"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impG0real(ilat,jlat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine ed_print_impG0


  !+------------------------------------------------------------------+
  !                      PRINT DENSITY MATRICES
  !+------------------------------------------------------------------+
  subroutine ed_print_dm(dm,Nrdm,ineq)
   integer                  ,intent(in)            :: Nrdm
   complex(8),dimension(:,:),intent(in)            :: dm
   integer                  ,intent(in),optional   :: ineq
   integer                                         :: unit,Nsites
   character(len=64)                               :: suffix
   integer                                         :: io,jo
   !
   if(size(dm,1)/=Nrdm.OR.size(dm,2)/=Nrdm)then
      stop "ERROR: actual dm argument has incogruent size wrt explicitly passed Nrdm"
   endif
   !
   Nsites = nint( 1/Norb * log(real(Nrdm,kind=8)) / log(4d0) ) !Nrdm = 4**(Nsites*Norb)
   !
   if(present(ineq))then
     suffix = "reduced_density_matrix_"//reg(str(Nsites))//"sites_ineq"//reg(str(ineq))//".dat"
   else
     suffix = "reduced_density_matrix_"//reg(str(Nsites))//"sites.dat"
   endif
   !
   unit = free_unit()
   open(unit,file=suffix,action="write",position="rewind",status='unknown')
   !
   do io=1,Nrdm
      write(unit,"(*(F20.16,1X))") (dreal(dm(io,jo)),jo=1,Nrdm)
   enddo
   write(unit,*)
   !
   if(any(dimag(dm)/=0d0))then
      do io=1,Nrdm
         write(unit,"(*(F20.16,1X))") (dimag(dm(io,jo)),jo=1,Nrdm)
      enddo
      write(unit,*)
   endif
   !
   close(unit)
   !
 end subroutine ed_print_dm

  !+------------------------------------------------------------------+
  ! !                         PRINT CHI:
  ! !+------------------------------------------------------------------+  
   subroutine ed_print_impChi
     call print_chi_spin
  !   call print_chi_dens
  !   call print_chi_dens_mix
  !   call print_chi_dens_tot
  !   call print_chi_pair
   end subroutine ed_print_impChi

  ! !                         SPIN-SPIN
  subroutine print_chi_spin
    integer                               :: ilat,jlat,iorb,jorb
    do ilat=1,Nlat
      do jlat=1,Nlat
        do iorb=1,Norb
           do jorb=1,Norb
              call splot("spinChi_i"//str(ilat)//str(jlat)//"_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,spinChi_tau(ilat,jlat,iorb,jorb,0:))
              call splot("spinChi_i"//str(ilat)//str(jlat)//"_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,spinChi_w(ilat,jlat,iorb,jorb,:))
              call splot("spinChi_i"//str(ilat)//str(jlat)//"_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,spinChi_iv(ilat,jlat,iorb,jorb,:))
           enddo
        enddo
      enddo
    enddo
  end subroutine print_chi_spin

  ! !                     DENSITY-DENSITY
  ! subroutine print_chi_dens
  !   integer                               :: i,j,iorb,jorb
  !   integer                               :: unit(3),unit_mix
  !   do iorb=1,Norb
  !      do jorb=iorb,Norb
  !         call splot("densChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,densChi_tau(iorb,jorb,0:))
  !         call splot("densChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",wr,densChi_w(iorb,jorb,:))
  !         call splot("densChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,densChi_iv(iorb,jorb,:))
  !      enddo
  !   enddo
  ! end subroutine print_chi_dens

  ! subroutine print_chi_dens_mix
  !   integer                               :: i,j,iorb,jorb
  !   integer                               :: unit(3),unit_mix
  !   do iorb=1,Norb
  !      do jorb=1,Norb
  !         call splot("densChi_mix_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,densChi_mix_tau(iorb,jorb,0:))
  !         call splot("densChi_mix_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",wr,densChi_mix_w(iorb,jorb,:))
  !         call splot("densChi_mix_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,densChi_mix_iv(iorb,jorb,:))
  !      enddo
  !   enddo
  ! end subroutine print_chi_dens_mix

  ! subroutine print_chi_dens_tot
  !   integer                               :: i,j,iorb,jorb
  !   integer                               :: unit(3),unit_mix
  !   call splot("densChi_tot_tau"//reg(ed_file_suffix)//".ed",tau,densChi_tot_tau(0:))
  !   call splot("densChi_tot_realw"//reg(ed_file_suffix)//".ed",wr,densChi_tot_w(:))
  !   call splot("densChi_tot_iw"//reg(ed_file_suffix)//".ed",vm,densChi_tot_iv(:))
  ! end subroutine print_chi_dens_tot


  ! !                             PAIR
  ! subroutine print_chi_pair
  !   integer                               :: i,iorb
  !   integer                               :: unit(3)
  !   do iorb=1,Norb
  !      call splot("pairChi_orb"//str(iorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,pairChi_tau(iorb,0:))
  !      call splot("pairChi_orb"//str(iorb)//"_realw"//reg(ed_file_suffix)//".ed",wr,pairChi_w(iorb,:))
  !      call splot("pairChi_orb"//str(iorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,pairChi_iv(iorb,:))
  !   enddo
  ! end subroutine print_chi_pair










   ! PURPOSE: Read self-energy function(s) - also for inequivalent sites.
   !+-----------------------------------------------------------------------------+!
   subroutine ed_read_impSigma_single
     integer                                           :: i,ispin,isign,unit(2),iorb,jorb,ilat,jlat
     character(len=30)                                 :: suffix
     integer,dimension(:),allocatable                  :: getIorb,getJorb
     integer                                           :: totNorb,l
     !
     if(.not.allocated(wm))allocate(wm(Lmats))
     if(.not.allocated(wr))allocate(wr(Lreal))
     wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
     wr     = linspace(wini,wfin,Lreal)
     !
     !!
     !Print the impurity functions:
     do ispin=1,Nspin
      do ilat=1,Nlat
        do jlat=1,Nlat
          do iorb=1,Norb
            do jorb=1,Norb
              suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
              call sread("impSigma"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impSmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
              call sread("impSigma"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impSreal(ilat,jlat,ispin,ispin,iorb,jorb,:))
            enddo
          enddo
        enddo
      enddo
     enddo
     !
     if(allocated(wm))deallocate(wm)
     if(allocated(wr))deallocate(wr)
     !
   end subroutine ed_read_impSigma_single
 
#if __GFORTRAN__ &&  __GNUC__ > 8      
   
   subroutine ed_read_impSigma_lattice(Nineq)
    integer :: Nineq
    integer :: isites
    !
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    !
    allocate(Smats_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sreal_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    !
    Smats_ineq  = zero 
    Sreal_ineq  = zero
    !
    do isites=1,Nineq
       ed_file_suffix=reg(ineq_site_suffix)//str(isites,site_indx_padding)
       call ed_read_impSigma_single
       Smats_ineq(isites,:,:,:,:,:,:,:)  = impSmats
       Sreal_ineq(isites,:,:,:,:,:,:,:)  = impSreal
    enddo
    ed_file_suffix=""
  end subroutine ed_read_impSigma_lattice
   
#endif   

   subroutine ed_read_impG_single
     integer                                           :: i,ispin,isign,unit(2),iorb,jorb,ilat,jlat
     character(len=30)                                 :: suffix
     integer,dimension(:),allocatable                  :: getIorb,getJorb
     integer                                           :: totNorb,l
     !
     if(.not.allocated(wm))allocate(wm(Lmats))
     if(.not.allocated(wr))allocate(wr(Lreal))
     wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
     wr     = linspace(wini,wfin,Lreal)
     !
     !Print the impurity functions:
     do ispin=1,Nspin
      do ilat=1,Nlat
        do jlat=1,Nlat
          do iorb=1,Norb
            do jorb=1,Norb
              suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
              call sread("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impGmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
              call sread("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impGreal(ilat,jlat,ispin,ispin,iorb,jorb,:))
            enddo
          enddo
        enddo
      enddo
     enddo
     !
     if(allocated(wm))deallocate(wm)
     if(allocated(wr))deallocate(wr)
     !
   end subroutine ed_read_impG_single

#if __GFORTRAN__ &&  __GNUC__ > 8      
   
   subroutine ed_read_impG_lattice(Nineq)
    integer :: Nineq
    integer :: isites
    !
    if(allocated(Gmats_ineq))deallocate(Gmats_ineq)
    if(allocated(Greal_ineq))deallocate(Greal_ineq)
    !
    allocate(Gmats_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Greal_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    !
    Gmats_ineq  = zero 
    Greal_ineq  = zero
    !
    do isites=1,Nineq
       ed_file_suffix=reg(ineq_site_suffix)//str(isites,site_indx_padding)
       call ed_read_impG_single
       Gmats_ineq(isites,:,:,:,:,:,:,:)  = impGmats
       Greal_ineq(isites,:,:,:,:,:,:,:)  = impGreal
    enddo
    ed_file_suffix=""
  end subroutine ed_read_impG_lattice
  
#endif

END MODULE ED_IO







