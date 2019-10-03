!This is a collection of auxiliary routines which are often used in 
!the drivers, but do not find space in the ED code:



!-------------------------------------------------------------------------------------------
!PURPOSE: periodization routines
!-------------------------------------------------------------------------------------------


subroutine periodize_g_scheme(kpoint)
   integer                                                     :: ilat,jlat,ispin,iorb,ii
   real(8),dimension(Ndim)                                     :: kpoint,ind1,ind2
   complex(8),allocatable,dimension(:,:,:,:,:,:)               :: gfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
   complex(8),allocatable,dimension(:,:)                       :: gfprime_lso ![Nlso][Nlso]
   complex(8),allocatable,dimension(:,:)                       :: Vk_lso ![Nlso][Nlso]
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: gfreal_unperiodized![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: gfmats_unperiodized ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
   !
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   allocate(gfmats_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(gfreal_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
   allocate(gfprime_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
   allocate(Vk_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
   if(.not.allocated(gfmats_periodized))allocate(gfmats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
   if(.not.allocated(gfreal_periodized))allocate(gfreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
   !
   gfmats_unperiodized=zero
   gfreal_unperiodized=zero
   gfprime=zero
   gfmats_periodized=zero
   gfreal_periodized=zero
   Vk_lso=zero
   gfprime_lso=zero
   !
   !
   Vk_lso=vca_nnn2lso_reshape(tk(kpoint)-t_prime,Nlat,Nspin,Norb)
   !
   do ii=1,Lmats    
      call vca_gf_cluster(xi*wm(ii),gfprime)
      gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
      call inv(gfprime_lso)
      gfprime_lso=gfprime_lso-Vk_lso
      call inv(gfprime_lso)
      gfmats_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(gfprime_lso,Nlat,Nspin,Norb)
   enddo
   !
   do ii=1,Lreal    
      call vca_gf_cluster(dcmplx(wr(ii),eps),gfprime)
      gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
      call inv(gfprime_lso)
      gfprime_lso=gfprime_lso-Vk_lso
      call inv(gfprime_lso)
      gfreal_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(gfprime_lso,Nlat,Nspin,Norb)
   enddo
   !
   do ii=1,Lmats
      do ilat=1,Nlat
         ind1=N2indices(ilat)        
         do jlat=1,Nlat
            ind2=N2indices(jlat)
            gfmats_periodized(:,:,:,:,ii)=gfmats_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*gfmats_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
         enddo
      enddo
   enddo
   !
   do ii=1,Lreal   
      do ilat=1,Nlat
         ind1=N2indices(ilat)        
         do jlat=1,Nlat
            ind2=N2indices(jlat)
            gfreal_periodized(:,:,:,:,ii)=gfreal_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*gfreal_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
         enddo
      enddo
   enddo
   !
   !if(allocated(wm))deallocate(wm)
   !if(allocated(wr))deallocate(wr) 
   !   
end subroutine periodize_g_scheme



subroutine build_sigma_g_scheme(kpoint)
   integer                                                     :: i,ispin,iorb,ii
   real(8),dimension(Ndim)                                     :: kpoint
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: invG0mats,invGmats
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: invG0real,invGreal
   !
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   if(.not.allocated(Smats_periodized))allocate(Smats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
   if(.not.allocated(Sreal_periodized))allocate(Sreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
   invG0mats = zero
   invGmats  = zero
   invG0real = zero
   invGreal  = zero
   Smats_periodized  = zero
   Sreal_periodized  = zero
   !
   !Get G0^-1
   !invG0mats = invg0_bath_mats(dcmplx(0d0,wm(:)),vca_bath)
   !invG0real = invg0_bath_real(dcmplx(wr(:),eps),vca_bath)
   !
   !Get G0^-1
   do ispin=1,Nspin
      do iorb=1,Norb
         do ii=1,Lmats
            invG0mats(ispin,ispin,iorb,iorb,ii) = (xi*wm(ii)+xmu)  + 2.d0*t_var*(cos(kpoint(1)) + cos(kpoint(2)))             ! FIXME: ad-hoc solution
         enddo
         do ii=1,Lreal
            invG0real(ispin,ispin,iorb,iorb,ii) = (wr(ii)+xmu)   + 2.d0*t_var*(cos(kpoint(1)) + cos(kpoint(2)))               ! FIXME: ad-hoc solution
         enddo
      enddo
   enddo
   !
   !Get Gimp^-1
   call periodize_g_scheme(kpoint)
   do ispin=1,Nspin
      do iorb=1,Norb
         invGmats(ispin,ispin,iorb,iorb,:) = one/gfmats_periodized(ispin,ispin,iorb,iorb,:)
         invGreal(ispin,ispin,iorb,iorb,:) = one/gfreal_periodized(ispin,ispin,iorb,iorb,:)
      enddo
   enddo
   !Get Sigma functions: Sigma= G0^-1 - G^-1
   Smats_periodized=zero
   Sreal_periodized=zero
   do ispin=1,Nspin
      do iorb=1,Norb
         Smats_periodized(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
         Sreal_periodized(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
      enddo
   enddo
   !
   !
   !if(allocated(wm))deallocate(wm)
   !if(allocated(wr))deallocate(wr)
   !
end subroutine build_sigma_g_scheme



subroutine print_periodized()
   character(len=64) :: suffix
   integer           :: iorb,ispin
   !
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   do iorb=1,Norb
      do ispin=1,Nspin
         suffix="_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//"_"//str(scheme)//"_scheme"
         call splot("perG"//reg(suffix)//"_iw.vca"   ,wm,gfmats_periodized(ispin,ispin,iorb,iorb,:))
         call splot("perG"//reg(suffix)//"_realw.vca",wr,gfreal_periodized(ispin,ispin,iorb,iorb,:))
         call splot("perSigma"//reg(suffix)//"_iw.vca"   ,wm,Smats_periodized(ispin,ispin,iorb,iorb,:))
         call splot("perSigma"//reg(suffix)//"_realw.vca",wr,Sreal_periodized(ispin,ispin,iorb,iorb,:))
      enddo
   enddo
   !
   if(allocated(wm))deallocate(wm)
   if(allocated(wr))deallocate(wr)
   !
end subroutine print_periodized
