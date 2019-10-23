!+-------------------------------------------------------------------+
!PURPOSE  : Deallocate the ED bath
!+-------------------------------------------------------------------+
subroutine deallocate_dmft_bath()
   integer :: ibath,isym
   if(.not.dmft_bath%status)return
   do ibath=1,Nbath
      dmft_bath%item(ibath)%v = 0d0
      call deallocate_h_repr(dmft_bath%item(ibath)%h)
   enddo
   deallocate(dmft_bath%item)
   dmft_bath%status=.false.
end subroutine deallocate_dmft_bath



!+-------------------------------------------------------------------+
!PURPOSE  : Allocate the ED bath
!+-------------------------------------------------------------------+
subroutine allocate_dmft_bath()
   integer :: ibath,isym,Nsym
   real(8) :: maxdiff
   !
   if(.not.allocated(impHloc))stop "impHloc not allocated in allocate_dmft_bath" !FIXME

   call deallocate_dmft_bath()
   !
   allocate(dmft_bath%item(Nbath))
   !
   !CHECK IF IDENDITY IS ONE OF THE SYMMETRIES, IF NOT ADD IT
   Nsym=impHloc_sym%N_dec
   !
   do isym=1,impHloc_sym%N_dec
      maxdiff=maxval(impHloc_sym%decomposition(isym)%O-lso2nnn_reshape(eye(Nlat*Nspin*Norb)))
      if(maxdiff .lt. 1d-6) Nsym=Nsym+1
   enddo
   !
   !ALLOCATE H MATRICES
   !
   do ibath=1,Nbath
      call allocate_h_repr(dmft_bath%item(ibath)%h,Nsym)
   enddo
   !
   dmft_bath%status=.true.
end subroutine allocate_dmft_bath






!+------------------------------------------------------------------+
!PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
!reading previous (converged) solution
!+------------------------------------------------------------------+
subroutine init_dmft_bath()
   real(8)              :: hrep_aux(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
   real(8)              :: re,im,maxdiff
   integer              :: i,ibath,isym,unit,flen,Nh,Nsym,Nsym_
   integer              :: io,jo,iorb,ispin,jorb,jspin
   logical              :: IOfile
   real(8)              :: de
   real(8)              :: offset_b(Nbath),noise_b(Nlat*Nspin*Norb)
   character(len=21)    :: space
   !  
   if(.not.dmft_bath%status)stop "init_dmft_bath error: bath not allocated"
   !
   if(Nbath>1)then
      offset_b=linspace(-HWBAND,HWBAND,Nbath)
   else
      offset_b(1)=0.d0
   endif
   !
   !BATH V INITIALIZATION
   do ibath=1,Nbath
      dmft_bath%item(ibath)%v=max(0.1d0,1.d0/sqrt(dble(Nbath)))
   enddo
   !
   !BATH MATRICES INITIALIZATION
   do ibath=1,Nbath
      Nsym = dmft_bath%item(ibath)%h%N_dec
      Nsym_= impHloc_sym%N_dec
      if(Nsym .ne. Nsym_)then
         dmft_bath%item(ibath)%h%decomposition(1)%O      =  lso2nnn_reshape(eye(Nlat*Nspin*Norb))
         dmft_bath%item(ibath)%h%decomposition(1)%lambda =  xmu+offset_b(ibath)
         do isym=2,Nsym
            dmft_bath%item(ibath)%h%decomposition(isym)%O      =  impHloc_sym%decomposition(isym)%O
            dmft_bath%item(ibath)%h%decomposition(isym)%lambda =  impHloc_sym%decomposition(isym)%lambda
         enddo
      else
         do isym=1,Nsym
            dmft_bath%item(ibath)%h%decomposition(isym)%O      =  impHloc_sym%decomposition(isym)%O
            dmft_bath%item(ibath)%h%decomposition(isym)%lambda =  impHloc_sym%decomposition(isym)%lambda
            maxdiff=maxval(impHloc_sym%decomposition(isym)%O-lso2nnn_reshape(eye(Nlat*Nspin*Norb)))
            if(maxdiff<1.d-6) dmft_bath%item(ibath)%h%decomposition(isym)%lambda =&
               dmft_bath%item(ibath)%h%decomposition(isym)%lambda + (xmu+offset_b(ibath))
         enddo
      endif
   enddo
   !
   !Read from file if exist:
   inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
   if(IOfile)then
      write(LOGfile,"(A)")"Reading bath from file "//trim(Hfile)//trim(ed_file_suffix)//".restart"
      unit = free_unit()
      flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
      !
      open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
      !
      !read V
      !
      do ibath=1,Nbath
         read(unit,"(F21.12,1X)")dmft_bath%item(ibath)%v
      enddo
      !
      !read number of lambdas
      read(unit,"(I3)")Nsym
      !read lambdas and O
      do isym=1,Nsym
         do ibath=1,Nbath
            read(unit,"(F21.12,1X)")dmft_bath%item(ibath)%h%decomposition(isym)%lambda
         enddo
         do io=1,Nlat*Nspin*Norb
            read(unit,*)(hrep_aux(io,jo),jo=1,Nlat*Nspin*Norb)
         enddo
         dmft_bath%item(ibath)%h%decomposition(isym)%O=lso2nnn_reshape(hrep_aux,Nlat,Nspin,Norb)
      enddo
      close(unit)
   endif
end subroutine init_dmft_bath








!+-------------------------------------------------------------------+
!PURPOSE  : write out the bath to a given unit with 
! the following column formatting: 
! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
!+-------------------------------------------------------------------+
subroutine write_dmft_bath(unit)
   integer,optional     :: unit
   integer              :: unit_
   integer              :: ibath
   integer              :: io,jo,iorb,ispin,isym
   complex(8)           :: hrep_aux_nnn(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
   real(8)              :: hrep_aux(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
   unit_=LOGfile;if(present(unit))unit_=unit
   if(.not.dmft_bath%status)stop "write_dmft_bath error: bath not allocated"
   !
   if(unit_==LOGfile)write(unit_,"(A9,a5,90(A9,1X))")"V"," ","H"        
   if(unit_==LOGfile)then
      do ibath=1,Nbath
         hrep_aux=0d0
         hrep_aux_nnn=H_from_symm(dmft_bath%item(ibath)%h)
         Hrep_aux=DREAL(nnn2lso_reshape(hrep_aux_nnn,Nlat,Nspin,Norb))
         write(unit_,"(F9.4,a5,90(F9.4,1X))")dmft_bath%item(ibath)%v,"|",( hrep_aux(1,jo),jo=1,Nlat*Nspin*Norb)        
         do io=2,Nlat*Nspin*Norb
            write(unit_,"(A9,a5,90(F9.4,1X))") "  "  ,"|",(hrep_aux(io,jo),jo=1,Nlat*Nspin*Norb)
         enddo
         write(unit_,"(A9)")" "
      enddo
   else
      do ibath=1,Nbath
         write(unit,"(90(F21.12,1X))")dmft_bath%item(ibath)%v
      enddo
      !
      !write number of lambdas
      write(unit,"(I3)")dmft_bath%item(ibath)%h%N_dec
      !write lambdas and O
      do isym=1,dmft_bath%item(ibath)%h%N_dec
         do ibath=1,Nbath
            write(unit,"(90(F21.12,1X))")dmft_bath%item(ibath)%h%decomposition(isym)%lambda
         enddo
         !O are the same for each replica
         hrep_aux=dmft_bath%item(1)%h%decomposition(isym)%O
         do io=1,Nlat*Nspin*Norb
            write(unit,*)(hrep_aux(io,jo),jo=1,Nlat*Nspin*Norb)
         enddo
      enddo
   endif
   !
end subroutine write_dmft_bath





!+-------------------------------------------------------------------+
!PURPOSE  : save the bath to a given file using the write bath
! procedure and formatting: 
!+-------------------------------------------------------------------+
subroutine save_dmft_bath(file,used)
   character(len=*),optional :: file
   character(len=256)        :: file_
   logical,optional          :: used
   logical                   :: used_
   character(len=16)         :: extension
   integer                   :: unit_
   if(.not.dmft_bath%status)stop "save_dmft_bath error: bath is not allocated"
   used_=.false.;if(present(used))used_=used
   extension=".restart";if(used_)extension=".used"
   file_=str(str(Hfile)//str(ed_file_suffix)//str(extension))
   if(present(file))file_=str(file)
   unit_=free_unit()
   open(unit_,file=str(file_))
   call write_dmft_bath(unit_)
   close(unit_)
end subroutine save_dmft_bath




!+-------------------------------------------------------------------+
!PURPOSE  : copy the bath components back to a 1-dim array 
!+-------------------------------------------------------------------+
subroutine set_dmft_bath(bath_)
   real(8),dimension(:)   :: bath_
   integer                :: stride,ibath,Nmask,io,jo,isym
   logical                :: check
   real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Haux
   !
   if(.not.dmft_bath%status)stop "get_dmft_bath error: bath not allocated"
   !
   !check=check_bath_dimension(bath_)
   !if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
   !
   do ibath=1,Nbath
      do isym=1,dmft_bath%item(ibath)%h%N_dec
         dmft_bath%item(ibath)%v=0d0
         dmft_bath%item(ibath)%h%decomposition(isym)%lambda=0d0
      enddo
   enddo
   !
   stride = 0
   !Get Vs
   do ibath=1,Nbath
      stride = stride + 1
      dmft_bath%item(ibath)%v = bath_(stride)
   enddo
   !Get lambdas
   do ibath=1,Nbath
      do isym=1,dmft_bath%item(ibath)%h%N_dec
         stride=stride+1
         dmft_bath%item(ibath)%h%decomposition(isym)%lambda=bath_(stride)
      enddo
   enddo
   !
end subroutine set_dmft_bath





!+-------------------------------------------------------------------+
!PURPOSE  : set the bath components from a given user provided 
! bath-array 
!+-------------------------------------------------------------------+
subroutine get_dmft_bath(bath_)
   real(8),dimension(:)   :: bath_
   integer                :: stride,ibath,Nmask,isym
   logical                :: check
   !
   if(.not.dmft_bath%status)stop "set_dmft_bath error: bath not allocated"
   !
   !check = check_bath_dimension(bath_)
   !if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
   !
   bath_=0.d0
   !
   stride = 0
   !Set Vs
   do ibath=1,Nbath
      stride = stride + 1
      bath_(stride)=dmft_bath%item(ibath)%v
   enddo
   !Set Lambdas
   do ibath=1,Nbath
      do isym=1,dmft_bath%item(ibath)%h%N_dec
         stride = stride + 1
         bath_(stride)=dmft_bath%item(ibath)%h%decomposition(isym)%lambda
      enddo
   enddo
end subroutine get_dmft_bath


























