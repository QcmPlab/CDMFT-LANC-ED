MODULE ED_AUX_FUNX
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_TIMER
  USE SF_LINALG
  USE SF_MISC, only: assert_shape
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  implicit none
  private




  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape

  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape

  interface print_Hloc
     module procedure :: print_Hloc_nnn
     module procedure :: print_Hloc_lso
  end interface print_Hloc

  !interface set_Hloc
     !module procedure set_Hloc_lso
     !module procedure set_Hloc_nnn
  !end interface set_Hloc

#if __GNUC__ > 6
  interface read(unformatted)
     module procedure :: read_unformatted
  end interface read(unformatted)

  interface write(unformatted)
     module procedure :: write_unformatted
  end interface write(unformatted)

  interface read(formatted)
     module procedure :: read_formatted
  end interface read(formatted)

  interface write(formatted)
     module procedure :: write_formatted
  end interface write(formatted)
#endif


  public :: index_stride_lso
  !
  !public :: set_Hloc
  public :: print_Hloc
  !
  public :: save_gfprime
  public :: read_gfprime
  !
  public :: lso2nnn_reshape
  public :: so2nn_reshape
  public :: nnn2lso_reshape
  public :: nn2so_reshape
  !
  public :: ed_search_variable
  !



contains



  !> Get stride position in the one-particle many-body space 
  function index_stride_lso(ilat,ispin,iorb) result(indx)
    integer :: ilat
    integer :: iorb
    integer :: ispin
    integer :: indx
    indx = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
  end function index_stride_lso



  !##################################################################
  !                   HLOC ROUTINES
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : Print Hloc
  !+------------------------------------------------------------------+
  subroutine print_Hloc_nnn(hloc,file)![Nspin][Nspin][Norb][Norb]
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hloc
    character(len=*),optional                          :: file
    integer                                            :: ilat,jlat,iorb,jorb,ispin,jspin
    integer                                            :: unit
    !
    unit=LOGfile
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    do ispin=1,Nspin
       do ilat=1,Nlat
          do iorb=1,Norb
             write(unit,"(20(F7.3,2x))")&
                  (((Hloc(ilat,jlat,ispin,jspin,iorb,jorb),jorb =1,Norb),jlat=1,Nlat),jspin=1,Nspin)
          enddo
       enddo
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hloc_nnn

  subroutine print_Hloc_lso(hloc,file) ![Nlso][Nlso]
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: hloc
    character(len=*),optional                          :: file
    integer                                            :: ilat,iorb,jorb,unit
    unit=LOGfile
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    !
    Nlso = Nlat*Nspin*Norb
    do iorb=1,Nlso
       write(unit,"(20(F7.3,2x))")(Hloc(iorb,jorb),jorb =1,Nlso)
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hloc_lso



  !##################################################################
  !                   RESHAPE ROUTINES
  !##################################################################
  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hnnn=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  !
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                               :: ilat,jlat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hnnn=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn
  !
  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
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
  end function d_nso2nn
  !
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                     :: Nspin,Norb
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
  end function c_nso2nn


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso
  !
  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: ilat,jlat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso
  !
  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
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
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                     :: Nspin,Norb
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
  end function c_nn2nso







#if __GNUC__ > 6
  !##################################################################
  !##################################################################
  ! ROUTINES TO READ AND WRITE CLUSTER GREEN FUNCTION
  ! unformatted and formatted I/O
  !##################################################################
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : write overload for GFmatrix type (formatted)
  !+-------------------------------------------------------------------+

  subroutine write_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
    class(GFmatrix), intent(in)         :: dtv
    integer, intent(in)                 :: unit
    integer, intent(out)                :: iostat
    character(*), intent(in)            :: iotype
    integer, intent(in)                 :: v_list(:)
    integer                             :: Nexc,iexc,Ichan,ilat,jlat,iorb,ispin,istate
    integer                             :: Nchan,Nstates
    character(*), intent(inout)         :: iomsg
    !
    !
    Nstates = size(dtv%state)
    write (unit, *,IOSTAT=iostat, IOMSG=iomsg) Nstates
    do istate=1,Nstates
       Nchan = size(dtv%state(istate)%channel)
       write (unit, *,IOSTAT=iostat, IOMSG=iomsg) Nchan
       do ichan=1,Nchan
          write (unit, *,IOSTAT=iostat, IOMSG=iomsg) size(dtv%state(istate)%channel(ichan)%poles)
          write (unit, *,IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%poles
          write (unit, *,IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%weight
       enddo
       write (unit, *,IOSTAT=iostat, IOMSG=iomsg) "\n"
    enddo
    !
  end subroutine write_formatted

  !+-------------------------------------------------------------------+
  !PURPOSE  : read overload for GFmatrix type (formatted)
  !+-------------------------------------------------------------------+

  subroutine read_formatted(dtv, unit,iotype, v_list, iostat, iomsg)
    class(GFmatrix), intent(inout)                :: dtv
    integer, intent(in)                           :: unit
    integer, intent(out)                          :: iostat
    character(*), intent(in)                      :: iotype
    integer, intent(in)                           :: v_list(:)
    character(*), intent(inout)                   :: iomsg
    logical                                       :: alloc
    integer                                       :: ichan,Nchan,Nlanc,istate,Nstates
    !
    read (unit,*,IOSTAT=iostat, IOMSG=iomsg) Nstates
    call GFmatrix_allocate(dtv,Nstate=Nstates)
    do istate=1,Nstates
       read (unit,*,IOSTAT=iostat, IOMSG=iomsg) Nchan
       call GFmatrix_allocate(dtv,istate=istate,Nchan=Nchan)
       do ichan=1,Nchan
          read (unit,*, IOSTAT=iostat, IOMSG=iomsg) Nlanc
          call GFmatrix_allocate(dtv,istate=istate,ichan=ichan,Nexc=Nlanc)
          read (unit, *, IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%poles
          read (unit, *, IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%weight
       enddo
    enddo
    !
  end subroutine read_formatted



  !+-------------------------------------------------------------------+
  !PURPOSE  : write overload for GFmatrix type (unformatted)
  !+-------------------------------------------------------------------+

  subroutine write_unformatted(dtv, unit, iostat, iomsg)
    class(GFmatrix), intent(in)         :: dtv
    integer, intent(in)                 :: unit
    integer, intent(out)                :: iostat
    integer                             :: Nexc,iexc,Ichan,ilat,jlat,iorb,ispin,istate
    integer                             :: Nchan,Nstates
    character(*), intent(inout)         :: iomsg
    !
    !
    Nstates = size(dtv%state)
    write (unit, IOSTAT=iostat, IOMSG=iomsg) Nstates
    do istate=1,Nstates
       Nchan = size(dtv%state(istate)%channel)
       write (unit, IOSTAT=iostat, IOMSG=iomsg) Nchan
       do ichan=1,Nchan
          write (unit, IOSTAT=iostat, IOMSG=iomsg) size(dtv%state(istate)%channel(ichan)%poles), dtv%state(istate)%channel(ichan)%poles, dtv%state(istate)%channel(ichan)%weight
       enddo
    enddo
    !
  end subroutine write_unformatted

  !+-------------------------------------------------------------------+
  !PURPOSE  : read overload for GFmatrix type (unformatted)
  !+-------------------------------------------------------------------+

  subroutine read_unformatted(dtv, unit, iostat, iomsg)
    class(GFmatrix), intent(inout)                :: dtv
    integer, intent(in)                           :: unit
    integer, intent(out)                          :: iostat
    character(*), intent(inout)                   :: iomsg
    logical                                       :: alloc
    integer                                       :: ichan,Nchan,Nlanc,istate,Nstates
    !
    read (unit, IOSTAT=iostat, IOMSG=iomsg) Nstates
    call GFmatrix_allocate(dtv,Nstate=Nstates)
    do istate=1,Nstates
       read (unit, IOSTAT=iostat, IOMSG=iomsg) Nchan
       call GFmatrix_allocate(dtv,istate=istate,Nchan=Nchan)
       do ichan=1,Nchan
          read (unit, IOSTAT=iostat, IOMSG=iomsg) Nlanc
          call GFmatrix_allocate(dtv,istate=istate,ichan=ichan,Nexc=Nlanc)
          read (unit, IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%poles
          read (unit, IOSTAT=iostat, IOMSG=iomsg) dtv%state(istate)%channel(ichan)%weight
       enddo
    enddo
    !
  end subroutine read_unformatted

#endif


  !+-------------------------------------------------------------------+
  !PURPOSE  : Save cluster GF to file
  !+-------------------------------------------------------------------+

  subroutine save_gfprime(file,used,use_formatted)

    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_
    logical,optional          :: use_formatted
    logical                   :: use_formatted_
    character(len=16)         :: extension
    integer                   :: unit_,Nchannel,Nexc,ichan,iexc,ilat,jlat,ispin,iorb,jorb
    !
#if __GNUC__ > 6
    if(.not.allocated(impGmatrix))stop "ed_gf_cluster ERROR: impGmatrix not allocated!"
    used_=.false.;if(present(used))used_=used
    use_formatted_=.false.;if(present(use_formatted))use_formatted_=use_formatted
    extension=".restart";if(used_)extension=".used"
    file_=str(str(file)//str(ed_file_suffix)//str(extension))
    unit_=free_unit()
    !
    if(use_formatted_)then
       open(unit_,file=str(file_),access='sequential')
    else
       open(unit_,file=str(file_),form='unformatted',access='sequential')
    endif
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if(use_formatted_)then
                      write(unit_,*)impGmatrix(ilat,jlat,ispin,ispin,iorb,jorb)
                   else
                      write(unit_)impGmatrix(ilat,jlat,ispin,ispin,iorb,jorb)
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    close(unit_)
#else
    print*,"Rear/write overloading requires Gfortran 6+"
#endif
  end subroutine save_gfprime

  !+-------------------------------------------------------------------+
  !PURPOSE  : Read cluster GF from file
  !+-------------------------------------------------------------------+

  subroutine read_gfprime(file,used,use_formatted)
    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_
    logical,optional          :: use_formatted
    logical                   :: use_formatted_
    character(len=16)         :: extension
    integer                   :: unit_,Nchannel,Nexc,ichan,iexc,ilat,jlat,ispin,iorb,jorb
    !
#if __GNUC__ > 6
    if(allocated(impGmatrix))deallocate(impGmatrix)
    allocate(impGmatrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    used_=.false.;if(present(used))used_=used
    use_formatted_=.false.;if(present(use_formatted))use_formatted_=use_formatted
    extension=".restart";if(used_)extension=".used"
    file_=str(str(file)//str(ed_file_suffix)//str(extension))
    unit_=free_unit()
    !
    if(use_formatted_)then
       open(unit_,file=str(file_),access='sequential')
    else
       open(unit_,file=str(file_),form='unformatted',access='sequential')
    endif
    !
    rewind(unit_)
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if(use_formatted_)then
                      read(unit_,*)impGmatrix(ilat,jlat,ispin,ispin,iorb,jorb)
                   else
                      read(unit_)impGmatrix(ilat,jlat,ispin,ispin,iorb,jorb)
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    close(unit_)
#else
    print*,"Rear/write overloading requires Gfortran 6+"
#endif
  end subroutine read_gfprime


  !##################################################################
  !##################################################################
  ! ROUTINES TO SEARCH CHEMICAL POTENTIAL UP TO SOME ACCURACY
  ! can be used to fix any other *var so that  *ntmp == nread
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_search_variable(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8),save          :: chich
    real(8),save          :: nold
    real(8),save          :: var_new
    real(8),save          :: var_old
    real(8)               :: var_sign
    !
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer               :: unit
    !
    !check actual value of the density *ntmp* with respect to goal value *nread*
    count=count+1
    totcount=totcount+1
    !  
    if(count==1)then
       chich = ndelta        !~0.2
       inquire(file="var_compressibility.restart",EXIST=bool)
       if(bool)then
          open(free_unit(unit),file="var_compressibility.restart")
          read(unit,*)chich
          close(unit)
       endif
       var_old = var
    endif
    !
    ndiff=ntmp-nread
    !
    !Get 'charge compressibility"
    if(count>1)chich = (ntmp-nold)/(var-var_old)
    !
    !Add here controls on chich: not to be too small....
    if(chich<0.1d0)chich=0.1d0*chich/abs(chich)
    !
    !update chemical potential
    var_new = var - ndiff/chich
    !
    !
    !re-define variables:
    nold    = ntmp
    var_old = var
    var     = var_new
    !
    !Print information
    write(LOGfile,"(A9,F16.9,A,F15.9)")  "n    = ",ntmp,"| instead of",nread
    write(LOGfile,"(A9,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nerr
    var_sign = (var-var_old)/abs(var-var_old)
    if(var_sign>0d0)then
       write(LOGfile,"(A9,ES16.9,A4)")"shift = ",ndiff/chich," ==>"
    else
       write(LOGfile,"(A9,ES16.9,A4)")"shift = ",ndiff/chich," <=="
    end if
    write(LOGfile,"(A9,F16.9)")"var  = ",var
    !
    !Save info about search variable iteration:
    open(free_unit(unit),file="search_variable_iteration_info"//reg(ed_file_suffix)//".ed",position="append")
    if(count==1)write(unit,*)"#var,ntmp,ndiff"
    write(unit,*)var,ntmp,ndiff
    close(unit)
    !
    !If density is not converged set convergence to .false.
    if(abs(ndiff)>nerr)converged=.false.
    !
    write(LOGfile,"(A18,I5)")"Search var count= ",count
    write(LOGfile,"(A19,L2)")"Converged       = ",converged
    print*,""
    !
    open(free_unit(unit),file="var_compressibility.used")
    write(unit,*)chich
    close(unit)
    !
  end subroutine ed_search_variable






END MODULE ED_AUX_FUNX
