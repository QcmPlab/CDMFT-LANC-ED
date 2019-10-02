!+-------------------------------------------------------------------+
!PURPOSE  : Deallocate the ED bath
!+-------------------------------------------------------------------+
subroutine deallocate_dmft_bath()
  integer :: ibath
  if(.not.dmft_bath%status)return
  do ibath=1,Nbath
     dmft_bath%item(ibath)%v = 0d0
     if(allocated(dmft_bath%item(ibath)%h))   deallocate(dmft_bath%item(ibath)%h)
  enddo
  if(allocated(dmft_bath%mask))deallocate(dmft_bath%mask)
  deallocate(dmft_bath%item)
  dmft_bath%Nmask = 0
  dmft_bath%status=.false.
end subroutine deallocate_dmft_bath



!+-------------------------------------------------------------------+
!PURPOSE  : Allocate the ED bath
!+-------------------------------------------------------------------+
subroutine allocate_dmft_bath()
  integer :: ibath
  !
  if(.not.allocated(impHloc))stop "impHloc not allocated in allocate_dmft_bath"

  call deallocate_dmft_bath()
  !
  allocate(dmft_bath%item(Nbath))
  !
  do ibath=1,Nbath
     allocate(dmft_bath%item(ibath)%h(Nlat,Nlat,Nspin,Nspin,Norb,Norb))    !replica hamilt of the bath
  enddo
  allocate(dmft_bath%mask(Nlat,Nlat,Nspin,Nspin,Norb,Norb)) !mask on components
  !
  dmft_bath%mask = mask_hloc(impHloc,wdiag=.true.,uplo=.true.)
  dmft_bath%Nmask= count(dmft_bath%mask)
  dmft_bath%status=.true.
end subroutine allocate_dmft_bath






!+------------------------------------------------------------------+
!PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
!reading previous (converged) solution
!+------------------------------------------------------------------+
subroutine init_dmft_bath()
  real(8)              :: hrep_aux(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
  real(8)              :: re,im
  integer              :: i,ibath,unit,flen,Nh
  integer              :: io,jo,iorb,ispin,jorb,jspin
  logical              :: IOfile
  real(8)              :: de
  real(8),allocatable  :: noise_b(:)
  character(len=21)    :: space
  !  
  if(.not.dmft_bath%status)stop "init_dmft_bath error: bath not allocated"
  !
  allocate(noise_b(Nbath));noise_b=0.d0 
  call random_number(noise_b)
  !
  !BATH INITIALIZATION
  do ibath=1,Nbath
     dmft_bath%item(ibath)%h=impHloc - (xmu+noise_b(ibath))*lso2nnn_reshape(eye(Nlat*Nspin*Norb),Nlat,Nspin,Norb)
     dmft_bath%item(ibath)%v=max(0.1d0,1.d0/sqrt(dble(Nbath)))
  enddo
  !
  !
  !Read from file if exist:
  inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
  if(IOfile)then
     write(LOGfile,"(A)")"Reading bath from file "//trim(Hfile)//trim(ed_file_suffix)//".restart"
     unit = free_unit()
     flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
     !
     open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
     do ibath=1,Nbath
        read(unit,"(F21.12,1X)")dmft_bath%item(ibath)%v
        do io=1,Nlat*Nspin*Norb
           read(unit,*)(hrep_aux(io,jo),jo=1,Nlat*Nspin*Norb)
        enddo
        dmft_bath%item(ibath)%h=lso2nnn_reshape(hrep_aux,Nlat,Nspin,Norb)
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
  integer              :: io,jo,iorb,ispin
  real(8)              :: hrep_aux(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
  unit_=LOGfile;if(present(unit))unit_=unit
  if(.not.dmft_bath%status)stop "write_dmft_bath error: bath not allocated"
  !
  if(unit_==LOGfile)write(unit_,"(A9,a5,90(A9,1X))")"V"," ","H"        
  do ibath=1,Nbath
     hrep_aux=0d0
     hrep_aux=nnn2lso_reshape(dmft_bath%item(ibath)%h,Nlat,Nspin,Norb)
     if(unit_==LOGfile)then
        write(unit_,"(F9.4,a5,90(F9.4,1X))")dmft_bath%item(ibath)%v,"|",( hrep_aux(1,jo),jo=1,Nlat*Nspin*Norb)        
        do io=2,Nlat*Nspin*Norb
           write(unit_,"(A9,a5,90(F9.4,1X))") "  "  ,"|",(hrep_aux(io,jo),jo=1,Nlat*Nspin*Norb)
        enddo
           write(unit_,"(A9)")" "
     else
        write(unit_,"(F21.12)")dmft_bath%item(ibath)%v
        do io=1,Nlat*Nspin*Norb
           write(unit_,"(90(F21.12,1X))")(hrep_aux(io,jo),jo=1,Nlat*Nspin*Norb)
        enddo
     endif
  enddo
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
  integer                :: stride,ibath,Nmask,io,jo
  logical                :: check
  real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Haux
  !
  if(.not.dmft_bath%status)stop "get_dmft_bath error: bath not allocated"
  !
  check=check_bath_dimension(bath_)
  if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
  !
  do ibath=1,Nbath
     dmft_bath%item(ibath)%h=0d0
     dmft_bath%item(ibath)%v=0d0
  enddo
  !
  Nmask  = dmft_bath%Nmask
  stride = 0
  !Get Hs
  do ibath=1,Nbath
     dmft_bath%item(ibath)%h = 0d0
     Haux=0.d0
     dmft_bath%item(ibath)%h = unpack(bath_(stride+1:stride+Nmask), dmft_bath%mask, dmft_bath%item(ibath)%h)
     Haux=nnn2lso_reshape(dmft_bath%item(ibath)%h,Nlat,Nspin,Norb)
     forall(io=1:Nlat*Nspin*Norb,jo=1:Nlat*Nspin*Norb,jo>io) Haux(jo,io)=Haux(io,jo) !conjugate if complex
     dmft_bath%item(ibath)%h=lso2nnn_reshape(Haux,Nlat,Nspin,Norb)
     stride = stride + Nmask
  enddo
  !
  !Get Vs
  do ibath=1,Nbath
     stride = stride + 1
     dmft_bath%item(ibath)%v = bath_(stride)
  enddo
end subroutine set_dmft_bath





!+-------------------------------------------------------------------+
!PURPOSE  : set the bath components from a given user provided 
! bath-array 
!+-------------------------------------------------------------------+
subroutine get_dmft_bath(bath_)
  real(8),dimension(:)   :: bath_
  integer                :: stride,ibath,Nmask
  logical                :: check
  !
  if(.not.dmft_bath%status)stop "set_dmft_bath error: bath not allocated"
  !
  check = check_bath_dimension(bath_)
  if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
  !
  bath_=0.d0
  !
  Nmask  = dmft_bath%Nmask
  stride = 0
  !Set Hs
  do ibath=1,Nbath
     bath_(stride+1:stride+Nmask) = pack(dmft_bath%item(ibath)%h,dmft_bath%mask)
     stride = stride + Nmask
  enddo
  !
  !Set Vs
  do ibath=1,Nbath
     stride = stride + 1
     bath_(stride)=dmft_bath%item(ibath)%v
  enddo
  !
end subroutine get_dmft_bath


























