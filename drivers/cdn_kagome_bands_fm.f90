program cdn_kanemele
   USE CDMFT_ED !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nlso, Nkpath
   real(8)                                                                :: ts,ts_perturb
   !
   character(len=16)                                                      :: finput
   !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
   integer                                                                :: comm
   integer                                                                :: rank
   integer                                                                :: mpi_size
   logical                                                                :: master,hermiticize

   !Init MPI: use of MPI overloaded functions in SciFor
   call init_MPI(comm,.true.)
   rank   = get_Rank_MPI(comm)
   master = get_Master_MPI(comm)
   !
   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(ts,"TS",finput,default=10.d0,comment="hopping parameter")
   call parse_input_variable(ts_perturb,"TS_PERTURB",finput,default=0.d0,comment="perturbation parameter")
   call parse_input_variable(Nkpath,"Nkpath",finput,default=30,comment="Number of k points along kpath")
   !
   call ed_read_input(trim(finput),comm)
   !
   !Add dmft control variables
   !
   call add_ctrl_var(beta,"BETA")
   call add_ctrl_var(Norb,"Norb")
   call add_ctrl_var(Nspin,"Nspin")
   call add_ctrl_var(xmu,"xmu")
   call add_ctrl_var(wini,"wini")
   call add_ctrl_var(wfin,"wfin")
   call add_ctrl_var(eps,"eps")

   !set global variables
   Nlat=3
   Nlso=Nlat*Nspin*Norb

   call generate_bands()

   call finalize_MPI()


contains

   !+------------------------------------------------------------------+
   !PURPOSE  : Hloc for the 2d Kagome model
   !+------------------------------------------------------------------+


   function hloc_model(ts_) result (hloc_)
      ! Directly created in lso basis
      integer                                               :: ispin, lowerbound, upperbound, Nlo
      real(8),intent(in)                                    :: ts_
      complex(8),dimension(Nlso,Nlso)                       :: hloc_
      !
      Nlo = Nlat*Norb
      hloc_  = zero
      !
      do ispin=1,Nspin
         lowerbound = (ispin-1)*Nlo+1
         upperbound = ispin*Nlo
         hloc_(lowerbound:upperbound,lowerbound:upperbound) = hloc_matrix(ts_)
      enddo
      !
   end function hloc_model


   function hk_model(kpoint,N) result(Hk_)
      ! Directly created in lso basis, kpoint MUST be in direct coordinates
      ! the N input is only needed or the TB_build_model routine....
      integer                         :: N, ispin, lowerbound, upperbound, Nlo
      real(8),dimension(:)            :: kpoint
      complex(8),dimension(N,N)       :: Hk_, hloc_, h1_, h2_, h3_, h4_, h5_, h6_,h_perturb,h_nonhermitian
      !
      if (N/=Nlso) stop "dimensionality of hk wrong!"
      Nlo   = Nlat*Norb
      Hk_   = zero
      hloc_ = zero
      h_perturb=zero
      h1_   = zero
      h2_   = zero
      h3_   = zero
      h4_   = zero
      h5_   = zero
      h6_   = zero
      !
      do ispin=1,Nspin
         lowerbound = (ispin-1)*Nlo+1
         upperbound = ispin*Nlo
         hloc_(lowerbound:upperbound,lowerbound:upperbound) = hloc_matrix(ts)
         h_perturb(lowerbound:upperbound,lowerbound:upperbound) = hloc_matrix(ts_perturb)
         h1_(lowerbound:upperbound,lowerbound:upperbound)   = hhop1_matrix(ts)
         h2_(lowerbound:upperbound,lowerbound:upperbound)   = hhop2_matrix(ts)
         h3_(lowerbound:upperbound,lowerbound:upperbound)   = hhop3_matrix(ts)
         h4_(lowerbound:upperbound,lowerbound:upperbound)   = hhop4_matrix(ts)
         h5_(lowerbound:upperbound,lowerbound:upperbound)   = hhop5_matrix(ts)
         h6_(lowerbound:upperbound,lowerbound:upperbound)   = hhop6_matrix(ts)
         h_nonhermitian(lowerbound:upperbound,lowerbound:upperbound) = xi*((-1.d0)**ispin)*zeye(Nlo)
      enddo
      !
      Hk_ = hloc_ + h1_*exp(2*pi*(0,1)*kpoint(1))              + h2_*exp(-2*pi*(0,1)*(kpoint(1))) &
                  + h3_*exp(2*pi*(0,1)*kpoint(2))              + h4_*exp(-2*pi*(0,1)*kpoint(2)) &
                  + h5_*exp(-2*pi*(0,1)*(kpoint(1)-kpoint(2))) + h6_*exp(2*pi*(0,1)*(kpoint(1)-kpoint(2))) &
                  + h_perturb &
                  + h_nonhermitian

      !
   end function hk_model


   !AUXILLIARY MATRIX CONSTRUCTORS for spin Up block, spindown = spinup

   function hloc_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 1, 1 ]
      tempmat(2,:) = [ 1, 0, 1 ]
      tempmat(3,:) = [ 1, 1, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hloc_matrix

   function hhop1_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 1, 0 ]
      tempmat(2,:) = [ 0, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop1_matrix

   function hhop2_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0 ]
      tempmat(2,:) = [ 1, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop2_matrix

   function hhop3_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 1 ]
      tempmat(2,:) = [ 0, 0, 0 ]
      tempmat(3,:) = [ 0, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop3_matrix

   function hhop4_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 0 ]
      tempmat(3,:) = [ 1, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop4_matrix

   function hhop5_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 1 ]
      tempmat(3,:) = [ 0, 0, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop5_matrix

   function hhop6_matrix(ts_) result(hmat)
      complex(8),dimension(3,3) :: hmat
      real(8),dimension(3,3)    :: tempmat
      real(8),intent(in)        :: ts_
      !
      hmat    = zero
      tempmat = zero
      !
      tempmat(1,:) = [ 0, 0, 0 ]
      tempmat(2,:) = [ 0, 0, 0 ]
      tempmat(3,:) = [ 0, 1, 0 ]
      hmat         = hmat - ts_*tempmat
      !
   end function hhop6_matrix

   subroutine generate_bands()
      integer                      :: ik,i,j
      real(8),dimension(4,2)       :: kpath
      real(8),dimension(2)         :: pointK, pointKp, pointM, bk1, bk2
      !
      bk1 = 2*pi/sqrt(3d0)*[1d0/2d0, sqrt(3d0)/2d0]
      bk2 = 2*pi*[0d0,1d0]
      !
      pointK  = [-1d0/3d0, 1d0/3d0]
      pointM  = [0d0, 1d0/2d0]
      pointKp = [-1d0/3d0, 2d0/3d0]
      write(LOGFile,*) pointK
      write(LOGFile,*) pointKp
      KPath(1,:)=[0d0,0d0]
      KPath(2,:)=pointK
      Kpath(3,:)=pointM
      KPath(4,:)=[0d0,0d0]
      write(LOGFile,*) KPath
      call TB_set_bk(bkx=bk1,bky=bk2)
      call solve_nh_model(hk_model,Nlso,KPath,Nkpath,&
           colors_name=[red1,blue1,red1,blue1,red1,blue1],&
           points_name=[character(len=10) :: "G","K","M","G"],&
           file="Eigenbands.nint",iproject=.false.)
      ! The distances on the kpath are not correctly scaled because the hamiltonian
      ! is implemented in direct coordinates
   end subroutine generate_bands


   !+------------------------------------------------------------------+
   !PURPOSE  : Auxilliary reshape functions
   !+------------------------------------------------------------------+

   function lso2nnn(Hlso) result(Hnnn)
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
                        is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                        js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                        Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function lso2nnn


   function nnn2lso(Hnnn) result(Hlso)
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
                        is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                        js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                        Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function nnn2lso
   !
  subroutine solve_nh_model(hk_model,Nlso,kpath,Nkpath,colors_name,points_name,file,iproject)
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         integer                   :: N
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
    integer                                   :: Nlso,io
    real(8),dimension(:,:)                    :: kpath
    integer                                   :: Nkpath
    type(rgb_color),dimension(Nlso)           :: colors_name
    character(len=*),dimension(size(kpath,1)) :: points_name
    character(len=*),optional                 :: file
    logical,optional                          :: iproject
    character(len=256)                        :: file_,file_real,file_imag
    logical                                   :: iproject_
    character(len=256)                        :: xtics
    integer                                   :: Npts,Ndim,Nktot
    integer                                   :: ipts,ik,ic,unit,iorb
    real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff,bk_x,bk_y,bk_z
    real(8)                                   :: coeff(Nlso),klen,ktics(size(Kpath,1))
    complex(8)                                :: h(Nlso,Nlso),evec(Nlso,Nlso),eval(Nlso)
    type(rgb_color)                           :: corb(Nlso),c(Nlso)
    character(len=10)                         :: chpoint
    character(len=32)                         :: fmt
    real(8),allocatable                       :: kseg(:)
    complex(8),allocatable                    :: Ekval(:,:)
    real(8),allocatable                       :: Ekval_sorted(:,:)
    integer,allocatable                       :: Ekcol(:,:)
    !
    master=.true.
    !
    file_    = "Eigenbands.tb";if(present(file))file_=file
    file_real=reg(file_)//".real"
    file_imag=reg(file_)//".imag"
    iproject_= .false.
    !
    Npts = size(kpath,1)
    Ndim = size(kpath,2)
    Nktot= (Npts-1)*Nkpath
    do iorb=1,Nlso
       corb(iorb) = colors_name(iorb)
    enddo
    !
    bk_x=[pi2,0d0,0d0]
    bk_y=[0d0,pi2,0d0]
    bk_z=[0d0,0d0,pi2]
    !  
    if(iproject_)then
       select case(Ndim)
       case (1)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
       case(2)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
       case (3)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
       end select
    endif
    !
    !
    if(master)then
       write(*,*)"Solving model along the path:"
       write(fmt,"(A3,I0,A)")"(A,",size(kpath,2),"F7.4,A1)"
       do ipts=1,Npts
          write(*,fmt)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
       enddo
    endif
    !
    ic = 0
    allocate(kseg(Nktot))
    allocate(ekval(Nktot,Nlso))
    allocate(ekval_sorted(Nktot,Nlso))
    allocate(ekcol(Nktot,Nlso))
    klen=0d0  
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/dble(Nkpath)
       ktics(ipts)  = klen
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hk_model(kpoint,Nlso)
          call eig(h,Eval,Evec)
          do iorb=1,Nlso
             coeff(:)=h(:,iorb)*conjg(h(:,iorb))
             c(iorb) = coeff.dot.corb
             Ekval(ic,iorb) = Eval(iorb)
             Ekcol(ic,iorb) = rgb(c(iorb))
          enddo
          kseg(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))
       enddo
    enddo
    ktics(Npts) = kseg(ic-1)
    !
    if(master)then
       open(free_unit(unit),file=str(file_real))
       do ic=1,Nktot
        Ekval_sorted(ic,:)=lazy_sort(REAL(Ekval(ic,:)))
       enddo
       do io=1,Nlso
          do ic=1,Nktot
             write(unit,*)kseg(ic),Ekval_sorted(ic,io),Ekcol(ic,io)
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       xtics=""
       xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
       do ipts=2,Npts-1
          xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
       enddo
       xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
       !
       open(unit,file=reg(file_real)//".gp")
       write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
       write(unit,*)"#set out '"//reg(file_real)//".png'"
       write(unit,*)""
       write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
       write(unit,*)"#set out '"//reg(file_real)//".svg'"
       write(unit,*)""
       write(unit,*)"#set term postscript eps enhanced color 'Times'"
       write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_real)//".pdf'"
       write(unit,*)"unset key"
       write(unit,*)"set xtics ("//reg(xtics)//")"
       write(unit,*)"set grid ytics xtics"
       !
       write(unit,*)"plot '"//reg(file_real)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_real)//".gp")
    endif
    !
    if(master)then
       open(free_unit(unit),file=str(file_imag))
       do ic=1,Nktot
        Ekval_sorted(ic,:)=lazy_sort(IMAG(Ekval(ic,:)))
       enddo
       do io=1,Nlso
          do ic=1,Nktot
             write(unit,*)kseg(ic),Ekval_sorted(ic,io),Ekcol(ic,io)
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       xtics=""
       xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
       do ipts=2,Npts-1
          xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
       enddo
       xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
       !
       open(unit,file=reg(file_imag)//".gp")
       write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
       write(unit,*)"#set out '"//reg(file_imag)//".png'"
       write(unit,*)""
       write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
       write(unit,*)"#set out '"//reg(file_imag)//".svg'"
       write(unit,*)""
       write(unit,*)"#set term postscript eps enhanced color 'Times'"
       write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_imag)//".pdf'"
       write(unit,*)"unset key"
       write(unit,*)"set xtics ("//reg(xtics)//")"
       write(unit,*)"set grid ytics xtics"
       !
       write(unit,*)"plot '"//reg(file_imag)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_imag)//".gp")
    endif
  end subroutine solve_nh_model
  
  
  
  subroutine bands_3d(hk_model,Nlso,Nk_)
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         integer                   :: N
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
    integer                                   :: Nlso,io,Nk_
    real(8),dimension(Nk_**2,2)               :: kpath
    type(rgb_color),dimension(Nlso)           :: colors_name
    character(len=256)                        :: file_,file_real,file_imag
    integer                                   :: Ndim,Nktot
    integer                                   :: ik,jk,unit,iorb,ipoint
    real(8),dimension(size(kpath,2))          :: bk_x,bk_y,kpoint
    complex(8)                                :: h(Nlso,Nlso),evec(Nlso,Nlso),eval(Nlso)
    character(len=32)                         :: fmt
    complex(8),allocatable                    :: Ekval(:,:)
    real(8),allocatable                       :: Ekval_sorted(:,:)
    !
    master=.true.
    !
    file_    = "bands_3d.ed"
    file_real=reg(file_)//".real"
    file_imag=reg(file_)//".imag"
    !
    Ndim = size(kpath,2)
    Nktot= Nk_**Ndim
    !
    bk_x=[pi2,0d0]
    bk_y=[0d0,pi2]
    !  
    !
    if(master)then
       write(*,*)"3d bands model:"
    endif
    !
    allocate(ekval(Nktot,Nlso))
    allocate(ekval_sorted(Nktot,Nlso))
    !
    call TB_build_kgrid([Nk_,Nk_],kpath)
    
    do ik=1,Nktot
      kpath(ik,:)=kpath(ik,:)-[pi,pi]
    enddo
    !
    do ik=1,Nktot
       kpoint = kpath(ik,:)
       h = hk_model(kpoint,Nlso)
       call eig(h,Eval,Evec)
       do iorb=1,Nlso
          Ekval(ik,iorb) = Eval(iorb)
       enddo
    enddo
    !
    if(master)then
       open(free_unit(unit),file=str(file_real))
       do ik=1,Nktot
        Ekval_sorted(ik,:)=lazy_sort(REAL(Ekval(ik,:)))
       enddo
       do io=1,2
          do ik=1,Nk_
             do jk=1,Nk_
               ipoint=ik*(Nk_-1)+jk
               write(unit,*)kpath(ipoint,1),kpath(ipoint,2),Ekval_sorted(ipoint,io)
             enddo
             write(unit,*)""
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       !
       open(unit,file=reg(file_real)//".gp")
       write(unit,*)"unset key"
       !
       write(unit,*)"splot '"//reg(file_real)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_real)//".gp")
    endif
    !
    if(master)then
       open(free_unit(unit),file=str(file_imag))
       do ik=1,Nktot
        Ekval_sorted(ik,:)=lazy_sort(IMAG(Ekval(ik,:)))
       enddo
       do io=1,Nlso
          do ik=1,Nk_
             do jk=1,Nk_
               ipoint=ik*(Nk_-1)+jk
               write(unit,*)kpath(ipoint,1),kpath(ipoint,2),Ekval_sorted(ipoint,io)
             enddo
             write(unit,*)""
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       !
       open(unit,file=reg(file_imag)//".gp")
       write(unit,*)"unset key"
       !
       write(unit,*)"plot '"//reg(file_imag)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_imag)//".gp")
    endif
  end subroutine bands_3d
  
  
  


  function lazy_sort(arr_in) result (arr)
    REAL(8), DIMENSION(:)                :: arr_in
    REAL(8), DIMENSION(:), allocatable   :: arr
    INTEGER                              :: i,j,inc,n
    REAL(8)                              :: v
    n=size(arr_in)
    allocate(arr(n))
    arr=arr_in
    inc=1
    do
      inc=3*inc+1
      if (inc > n) exit
    end do
    do
      inc=inc/3
      do i=inc+1,n
        v=arr(i)
        j=i
        do
          if (arr(j-inc) <= v) exit
          arr(j)=arr(j-inc)
          j=j-inc
          if (j <= inc) exit
        end do
        arr(j)=v
      end do
      if (inc <= 1) exit
    end do
  end function lazy_sort
   !
end program cdn_kanemele
