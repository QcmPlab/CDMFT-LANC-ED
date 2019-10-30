!This is a collection of auxiliary routines which are often used in 
!the drivers, but do not find space in the ED code:

!-------------------------------------------------------------------------------------------
!PURPOSE: periodization routines, G-scheme
!-------------------------------------------------------------------------------------------

subroutine periodize_g_scheme(kpoint,gmats_periodized,greal_periodized,hk_unper)
   integer                                                     :: ilat,jlat,ispin,iorb,ii
   real(8),dimension(:)                                        :: kpoint
   integer,dimension(:),allocatable                            :: ind1,ind2
   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: tmpmat,Hk_unper
   complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: gmats_unperiodized![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
   complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal) :: greal_unperiodized ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: gmats_periodized
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: greal_periodized
   !
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
   if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
   !
   gmats_unperiodized=zero
   greal_unperiodized=zero
   gmats_periodized=zero
   greal_periodized=zero
   tmpmat=zero
   !
   do ii=1,Lmats
      tmpmat=(xi*wm(ii)+xmu)*eye(Nlat*Nspin*Norb) - hk_unper - nnn2lso(Smats(:,:,:,:,:,:,ii))
      call inv(tmpmat)
      gmats_unperiodized(:,:,:,:,:,:,ii)=lso2nnn(tmpmat)
   enddo
   !
   do ii=1,Lreal
      tmpmat=(wr(ii)+xmu)*eye(Nlat*Nspin*Norb) - hk_unper - nnn2lso(Sreal(:,:,:,:,:,:,ii))
      call inv(tmpmat)
      greal_unperiodized(:,:,:,:,:,:,ii)=lso2nnn(tmpmat)
   enddo
   !
   do ii=1,Lmats
      do ilat=1,Nlat
         ind1=N2indices(ilat)        
         do jlat=1,Nlat
            ind2=N2indices(jlat)
            gmats_periodized(:,:,:,:,ii)=gmats_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*gmats_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
         enddo
      enddo
   enddo
   !
   do ii=1,Lreal   
      do ilat=1,Nlat
         ind1=N2indices(ilat)        
         do jlat=1,Nlat
            ind2=N2indices(jlat)
            greal_periodized(:,:,:,:,ii)=greal_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*greal_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
         enddo
      enddo
   enddo
   !
   deallocate(ind1,ind2)
   !if(allocated(wm))deallocate(wm)
   !if(allocated(wr))deallocate(wr) 
   !   
end subroutine periodize_g_scheme



subroutine build_sigma_g_scheme(kpoint,gmats_periodized,greal_periodized,smats_periodized,sreal_periodized,Hk_unper,Hk_per)
   integer                                                     :: i,ispin,iorb,ii
   real(8),dimension(:)                                        :: kpoint
   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: Hk_unper
   complex(8),dimension(Nspin*Norb,Nspin*Norb)                 :: Hk_per
   complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats)           :: invG0mats,invGmats
   complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal)           :: invG0real,invGreal
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: gmats_periodized, Smats_periodized
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: greal_periodized, Sreal_periodized
   !
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   invG0mats = zero
   invGmats  = zero
   invG0real = zero
   invGreal  = zero
   Smats_periodized  = zero
   Sreal_periodized  = zero
   !
   !Get G0^-1
   do ii=1,Lmats
      invG0mats(:,:,ii) = (xi*wm(ii)+xmu)  - Hk_per            
   enddo
   do ii=1,Lreal
      invG0real(:,:,ii) = (wr(ii)+xmu)   - Hk_per   
   enddo
   !
   !Get Gimp^-1
   call periodize_g_scheme(kpoint,gmats_periodized,greal_periodized,Hk_unper)
   do ii=1,Lmats
      invGmats(:,:,ii) = nn2so(gmats_periodized(:,:,:,:,ii))
      call inv(invGmats(:,:,ii))
   enddo
   do ii=1,Lmats
      invGreal(:,:,ii) = nn2so(greal_periodized(:,:,:,:,ii))
      call inv(invGreal(:,:,ii))
   enddo
   !
   !Get Sigma functions: Sigma= G0^-1 - G^-1
   Smats_periodized=zero
   Sreal_periodized=zero
   !
   do ii=1,Lmats
      Smats_periodized(:,:,:,:,ii) = so2nn(invG0mats(:,:,ii) - invGmats(:,:,ii))
   enddo
   do ii=1,Lreal
      Sreal_periodized(:,:,:,:,ii) = so2nn(invG0real(:,:,ii) - invGreal(:,:,ii))
   enddo
   !
   !if(allocated(wm))deallocate(wm)
   !if(allocated(wr))deallocate(wr)
   !
end subroutine build_sigma_g_scheme


!-------------------------------------------------------------------------------------------
!PURPOSE: periodization routines, SIGMA-scheme
!-------------------------------------------------------------------------------------------

subroutine periodize_sigma_scheme(kpoint,smats_periodized,sreal_periodized)
   integer                                                     :: ilat,jlat,ispin,iorb,ii
   real(8),dimension(:)                                        :: kpoint
   integer,dimension(:),allocatable                            :: ind1,ind2
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: smats_periodized
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: sreal_periodized
   !
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
   if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
   !
   smats_periodized=zero
   sreal_periodized=zero
   !
   do ii=1,Lmats
      do ilat=1,Nlat
         ind1=N2indices(ilat)        
         do jlat=1,Nlat
            ind2=N2indices(jlat)
            smats_periodized(:,:,:,:,ii)=smats_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*Smats(ilat,jlat,:,:,:,:,ii)/Nlat
         enddo
      enddo
   enddo
   !
   do ii=1,Lreal   
      do ilat=1,Nlat
         ind1=N2indices(ilat)        
         do jlat=1,Nlat
            ind2=N2indices(jlat)
            sreal_periodized(:,:,:,:,ii)=sreal_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*Sreal(ilat,jlat,:,:,:,:,ii)/Nlat
         enddo
      enddo
   enddo
   !
   deallocate(ind1,ind2)
   !if(allocated(wm))deallocate(wm)
   !if(allocated(wr))deallocate(wr) 
   !   
end subroutine periodize_sigma_scheme



subroutine build_g_sigma_scheme(kpoint,gmats_periodized,greal_periodized,smats_periodized,sreal_periodized,Hk_per)
   integer                                                     :: i,ispin,iorb,ii
   real(8),dimension(:)                                        :: kpoint
   complex(8),dimension(Nspin*Norb,Nspin*Norb)                 :: Hk_per,tmpmat
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: gmats_periodized, Smats_periodized
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: greal_periodized, Sreal_periodized
   !
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   !Get Gimp^-1
   call periodize_sigma_scheme(kpoint,smats_periodized,sreal_periodized)
   !
   do ii=1,Lmats
      tmpmat=(xi*wm(ii)+xmu)*eye(Nspin*Norb) - hk_per - nn2so(Smats_periodized(:,:,:,:,ii))
      call inv(tmpmat)
      gmats_periodized(:,:,:,:,ii)=so2nn(tmpmat)
   enddo
   !
   do ii=1,Lreal
      tmpmat=(wr(ii)+xmu)*eye(Nspin*Norb) - hk_per - nn2so(Sreal_periodized(:,:,:,:,ii))
      call inv(tmpmat)
      greal_periodized(:,:,:,:,ii)=so2nn(tmpmat)
   enddo
   !
   !if(allocated(wm))deallocate(wm)
   !if(allocated(wr))deallocate(wr)
   !
end subroutine build_g_sigma_scheme


!-------------------------------------------------------------------------------------------
!PURPOSE: user interface
!-------------------------------------------------------------------------------------------

subroutine print_periodized(Nkpts,func1,func2,scheme)
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)  :: gmats_periodized, Smats_periodized, Gloc_per_iw, Sloc_per_iw
   complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)  :: greal_periodized, Sreal_periodized, Gloc_per_rw, Sloc_per_rw
   integer,dimension(:)                               :: nkpts
   real(8),dimension(product(Nkpts),size(Nkpts))      :: kgrid
   integer                                            :: i
   character(len=6)                                   :: scheme
   interface
     function func1(x,N)
       integer                   :: N
       real(8),dimension(:)      :: x
       complex(8),dimension(N,N) :: func1
     end function func1
   end interface
   interface
     function func2(x,N)
       integer                   :: N
       real(8),dimension(:)      :: x
       complex(8),dimension(N,N) :: func2
     end function func2
   end interface
   !
   Gloc_per_iw=zero
   Sloc_per_iw=zero
   Gloc_per_rw=zero
   Sloc_per_rw=zero
   !
   call TB_build_kgrid(nkpts,kgrid)
   do i=1,size(Nkpts)
      kgrid(:,i)=kgrid(:,i)/Nkpts(i)
   enddo
   !
   if(ED_VERBOSE .ge. 1)write(LOGfile,*)"Computing periodized quantities using    ",scheme," scheme"
   !
   if(MASTER)call start_timer
   do i=1,product(Nkpts)
      if(scheme=="g") then
         call build_sigma_g_scheme(kgrid(i,:),gmats_periodized,greal_periodized,smats_periodized,&
            sreal_periodized,func1(kgrid(i,:),Nlat*Nspin*Norb),func2(kgrid(i,:),Nspin*Norb))
      elseif(scheme=="sigma")then
         call build_g_sigma_scheme(kgrid(i,:),gmats_periodized,greal_periodized,smats_periodized,&
            sreal_periodized,func2(kgrid(i,:),Nspin*Norb))
      else
         STOP "Nonexistent periodization scheme"
      endif
      Gloc_per_iw=Gloc_per_iw+(gmats_periodized/product(Nkpts))
      Sloc_per_iw=Sloc_per_iw+(smats_periodized/product(Nkpts))
      Gloc_per_rw=Gloc_per_rw+(greal_periodized/product(Nkpts))
      Sloc_per_rw=Sloc_per_rw+(sreal_periodized/product(Nkpts))
      if(ED_VERBOSE .ge. 1)call eta(i,product(Nkpts))
   enddo 
   if(MASTER)call stop_timer
   !
   if(master)call dmft_print_gf_matsubara(gloc_per_iw,"Gloc_periodized",iprint=4)
   if(master)call dmft_print_gf_matsubara(sloc_per_iw,"Sigma_periodized",iprint=4)
   !
   if(master)call dmft_print_gf_realaxis(gloc_per_rw,"Gloc_periodized",iprint=4)
   if(master)call dmft_print_gf_realaxis(sloc_per_rw,"Sigma_periodized",iprint=4)
   !
end subroutine print_periodized


!-------------------------------------------------------------------------------------------
!PURPOSE: auxiliary
!-------------------------------------------------------------------------------------------

function so2nn(Hso) result(Hnn)
  complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
  integer                                     :: iorb,ispin,is
  integer                                     :: jorb,jspin,js
  Hnn=zero
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
               is = iorb + (ispin-1)*Norb  !spin-orbit stride
              js = jorb + (jspin-1)*Norb  !spin-orbit stride
              Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
           enddo
        enddo
     enddo
  enddo
end function so2nn
!
function nn2so(Hnn) result(Hso)
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
  complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
  integer                                     :: iorb,ispin,is
  integer                                     :: jorb,jspin,js
  Hso=zero
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              is = iorb + (ispin-1)*Norb  !spin-orbit stride
              js = jorb + (jspin-1)*Norb  !spin-orbit stride
              Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
           enddo
        enddo
     enddo
  enddo
end function nn2so
   !
