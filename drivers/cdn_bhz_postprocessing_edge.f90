program cdn_bhz_postprocessing_edge
   USE CDMFT_ED
   !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   !USE MPI
   !
   implicit none
   integer                                                                :: Nx,Nso,Nlso,Nilso,iloop,Nb,Nkx,iw,Ly,Nkpath
   integer                                                                :: Nineq,Nsites,ineq,isites,icounter,icounter2
   integer,dimension(2)                                                   :: recover
   logical,allocatable,dimension(:)                                       :: converged_sites
   logical                                                                :: converged,lrsym
   real(8)                                                                :: ts,Mh,lambda,wmixing,observable,repart,impart
   !Bath:
   real(8),allocatable                                                    :: Bath(:,:),Bath_prev(:,:)
   !The local hybridization function:
   complex(8),allocatable                                                 :: Hloc(:,:)
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Hloc_nnnn,Hloc_ineq
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Gmats,Greal,Smats,Sreal
   complex(8),allocatable,dimension(:,:,:,:,:,:,:,:)                      :: Gmats_all,Greal_all,Smats_all,Sreal_all
   complex(8),allocatable,dimension(:,:,:,:,:,:,:,:)                      :: Gmats_ineq,Greal_ineq,Smats_ineq,Sreal_ineq
   character(len=16)                                                      :: finput
   real(8),allocatable                                                    :: wr(:)
   complex(8),allocatable                                                 :: wm(:)
   complex(8),allocatable                                                 :: Hk(:,:,:),Smats_lso(:,:,:)
   !SYMMETRIES TEST
   real(8),dimension(:,:,:),allocatable                                   :: lambdasym_vectors
   complex(8),dimension(:,:,:,:,:,:,:),allocatable                        :: Hsym_basis
   character(len=6)                                                       :: scheme
   type(finter_type)                                                      :: finter_func
   complex(8),dimension(:,:,:),allocatable                                :: dummy_real
   complex(8),dimension(:,:,:),allocatable                                :: dummy_mats
   !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
   integer                                                                :: comm
   integer                                                                :: rank
   integer                                                                :: mpi_size
   logical                                                                :: master

   !Init MPI: use of MPI overloaded functions in SciFor
   call init_MPI(comm,.true.)
   rank   = get_Rank_MPI(comm)
   master = get_Master_MPI(comm)

   !

   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
   call parse_input_variable(Mh,"Mh",finput,default=1.d0,comment="crystal field splitting")
   call parse_input_variable(lambda,"lambda",finput,default=0.3d0,comment="spin-orbit coupling")
   call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
   call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
   call parse_input_variable(Ly,"Ly",finput,default=2,comment="Number of sites in the finite y direction")
   call parse_input_variable(lrsym,"LRSYM",finput,default=.true.)
   call parse_input_variable(scheme,"SCHEME",finput,default="sigma")
   call parse_input_variable(nkpath,"NKPATH",finput,default=100)
   !
   call ed_read_input(trim(finput))
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
   Nsites = Ly
   Nineq= Ly
   if(lrsym)then
     if(mod(Ly,2)/=0)stop "Wrong setup from input file: Ly%2 > 0 (odd number of sites)"
     Nineq=Ly/2
     print*,"Using L-R Symmetry. Solve",Nineq," of",Nsites," sites."
     call sleep(2)
   endif
   Nlat=Nx
   Nso=Nspin*Norb
   Nlso=Nlat*Nspin*Norb
   Nilso=Nsites*Nlat*Nspin*Norb
   !
   if(.not.allocated(wm))allocate(wm(Lmats))
   if(.not.allocated(wr))allocate(wr(Lreal))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)
   wr     = linspace(wini,wfin,Lreal)
   !
   !Allocate:
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
   allocate(Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
   allocate(Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
   !
   allocate(Smats_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats_ineq=zero
   allocate(Sreal_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal_ineq=zero
   allocate(Gmats_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats_ineq=zero
   allocate(Greal_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal_ineq=zero
   !
   allocate(Smats_all(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats_all=zero
   allocate(Sreal_all(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal_all=zero
   allocate(Gmats_all(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats_all=zero
   allocate(Greal_all(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal_ineq=zero
   !
   allocate(Hloc_nnnn(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb));Hloc_nnnn=zero
   allocate(Hloc_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
   allocate(dummy_mats(Nsites*Nso,Nsites*Nso,Lmats),dummy_real(Nsites*Nso,Nsites*Nso,Lreal))

   !Build Hk and Hloc
   call generate_hk_hloc()
   Hloc_nnnn = ilso2nnnn_blocks(Hloc)
   do ineq=1,Nineq
     isites = ineq2isites(ineq)
     Hloc_ineq(ineq,:,:,:,:,:,:) = Hloc_nnnn(isites,:,:,:,:,:,:)
   enddo
   
   !SETUP SYMMETRIE
   allocate(lambdasym_vectors(Nineq,Nbath,3))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,3))
   !
   do icounter=1,Nineq
    lambdasym_vectors(icounter,:,1)=Mh
   enddo
   Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(zeye(Nlso))
   !
   do icounter=1,Nineq
    lambdasym_vectors(icounter,:,2)=ts
   enddo
   Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(zeye(Nlso))
   !
   do icounter=1,Nineq
    lambdasym_vectors(icounter,:,3)=lambda
   enddo
   Hsym_basis(:,:,:,:,:,:,3)=lso2nnn(zeye(Nlso))
   !
   !setup solver
   call ed_set_Hreplica(Hsym_basis,lambdasym_vectors)
   Nb=ed_get_bath_dimension(Hsym_basis)
   allocate(bath(Nineq,Nb))
   allocate(bath_prev(Nineq,Nb))
   call ed_init_solver(comm,bath)
   !
   !call ed_read_impG(Nineq)
   call ed_read_impSigma(Nineq)
   !
   !call ed_get_gimp_matsubara(Gmats_ineq,Nineq)
   !call ed_get_gimp_realaxis(Greal_ineq,Nineq)
   !call ed_get_sigma_matsubara(Smats_ineq,Nineq)
   call ed_get_sigma_realaxis(Sreal_ineq,Nineq)
   !
   do isites=1,Nsites
     ineq = isites2ineq(isites)
     !Smats_all(isites,:,:,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:,:,:)
     Sreal_all(isites,:,:,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:,:,:)
     !Gmats_all(isites,:,:,:,:,:,:,:) = Gmats_ineq(ineq,:,:,:,:,:,:,:)
     !Greal_all(isites,:,:,:,:,:,:,:) = Greal_ineq(ineq,:,:,:,:,:,:,:)
   enddo
   !
   !
   !RETRIEVE AND PERIODIZE
   !call get_zeros()
   call get_akw()
   !

contains


   !+------------------------------------------------------------------+
   !PURPOSE  : the BHZ-edge model hamiltonian
   !+------------------------------------------------------------------+
   
 function bhz_edge_model(kpoint,Nsites,N,pbc) result(Hrk)
    real(8),dimension(:)                    :: kpoint
    real(8)                                 :: kx
    integer                                 :: Nsites,N
    complex(8),dimension(N,N)               :: Hmat,Tmat,TmatH
    complex(8),dimension(Nsites*N,Nsites*N) :: Hrk
    integer                                 :: i,Idmin,Idmax,Itmin,Itmax
    logical                                 :: pbc
    kx=kpoint(1)
    Hrk=zero
    Hmat=hk_model(kx,N)
    Tmat=t0_rk_bhz(N)
    TmatH=conjg(transpose(Tmat))
    do i=1,Nsites
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat 
    enddo
    if(Nsites>1)then
      do i=1,Nsites-1
         Idmin=1 + (i-1)*N
         Idmax=        i*N
         Itmin=1 +     i*N
         Itmax=    (i+1)*N
         Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
         Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
      enddo
    endif
    if(pbc)then
       Itmin=1+(Nsites-1)*N
       Itmax=0+Nsites*N
       Hrk(1:N,Itmin:Itmax)=TmatH
       Hrk(Itmin:Itmax,1:N)=Tmat
    endif
  end function bhz_edge_model


  function t0_rk_bhz(N) result(H)
    integer                    :: N,i,Idmin,Idmax
    complex(8),dimension(N,N)  :: H
    H=zero
    do i=1,N/2
       Idmin=2*i-1
       Idmax=2*i
       H(Idmin:Idmax,Idmin:Idmax)=t_y()
    enddo
  end function t0_rk_bhz


   function hk_model(kpoint,N) result(Hk)
      integer                                                      :: N,ilat,jlat,ispin,iorb,jorb,i,j,ind1,ind2
      real(8)                                                      :: kpoint
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
      complex(8),dimension(N,N)                                    :: hk
      !
      hopping_matrix=zero
      !
      !
      do ispin=1,Nspin
         do ilat=1,Nx
           hopping_matrix(ilat,ilat,ispin,ispin,:,:)= t_m()
           if(ilat<Nx)then
              hopping_matrix(ilat+1,ilat,ispin,ispin,:,:)= t_x(ispin)
           endif
           if(ilat>1)then
              hopping_matrix(ilat-1,ilat,ispin,ispin,:,:)= dconjg(transpose(t_x(ispin)))
           endif
         enddo 
         hopping_matrix(Nx,1,ispin,ispin,:,:)=hopping_matrix(Nx,1,ispin,ispin,:,:) + dconjg(transpose(t_x(ispin)))*exp(xi*kpoint*Nx)
         hopping_matrix(1,Nx,ispin,ispin,:,:)=hopping_matrix(1,Nx,ispin,ispin,:,:) + t_x(ispin)*exp(-xi*kpoint*Nx)
      enddo
      !
      Hk=nnn2lso(hopping_matrix)
      !
   end function hk_model


   !AUXILLIARY HOPPING MATRIX CONSTRUCTORS

   function t_m() result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      !
      tmpmat=zero
      tmpmat=Mh*pauli_sigma_z
      !
   end function t_m

   function t_x(spinsign) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: sz
      integer                         :: spinsign
      !
      tmpmat=zero
      sz=(-1.d0)**(spinsign+1)
      tmpmat=-ts*pauli_sigma_z+0.5d0*sz*xi*lambda*pauli_sigma_x
      !
   end function t_x

   function t_y() result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      !
      tmpmat=zero
      tmpmat=-ts*pauli_sigma_z
      tmpmat(1,2)=-lambda*0.5d0
      tmpmat(2,1)=lambda*0.5d0
      !
   end function t_y

   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate Hloc and Hk
   !-------------------------------------------------------------------------------------------

   subroutine generate_hk_hloc()
      integer                                     :: ik
      real(8),dimension(Nkx,3)                    :: kgrid
      real(8),dimension(3)                        :: e1,e2,e3,bk1,bk2,bk3
      real(8)                                     :: bklen
      !
      e1 = [1d0, 0d0, 0d0]
      e2 = [0d0, 1d0, 0d0]
      e3 = [0d0, 0d0, 1d0]
      call TB_set_ei(eix=e1,eiy=e2,eiz=e3)
      bklen=2d0*pi
      bk1=bklen*[1d0, 0d0, 0d0]
      bk2=bklen*[0d0, 1d0, 0d0]
      bk3=bklen*[0d0, 0d0, 1d0]
      call TB_set_bk(bkx=bk1,bky=bk2,bkz=bk3)
      !
      call TB_build_kgrid([Nkx,1,1],kgrid)
      kgrid(:,1)=kgrid(:,1)/Nx
      !
      if(allocated(hk))deallocate(hk)
      if(allocated(hloc))deallocate(Hloc)
      !
      allocate(Hk(Nilso,Nilso,Nkx),Hloc(Nilso,Nilso))
      hk=zero
      hloc=zero
      !
      ! SEVER !
      !do ik=1,size(kgrid,1)
      !  hk(:,:,ik)=bhz_edge_model(kgrid(ik,:),Nsites,Nlso,.false.)
      !enddo
      call TB_build_model(Hk,bhz_edge_model,Ly,Nlso,kgrid,pbc=.false.,wdos=.false.)
      Hloc = sum(Hk,dim=3)/Nkx
      where(abs(Hloc)<1d-6)Hloc=zero
      !
   end subroutine generate_hk_hloc


   !+------------------------------------------------------------------+
   !PURPOSE  : Auxilliary reshape functions
   !+------------------------------------------------------------------+
   
  function isites2ineq(isites) result(ineq)
    integer,intent(in) :: isites
    integer            :: ineq
    ineq=isites
    if( lrsym .AND. (isites>Nineq) )ineq=Nsites-isites+1
  end function isites2ineq

  function ineq2isites(ineq) result(isites)
    integer,intent(in) :: ineq
    integer            :: isites
    isites=ineq
    if(ineq>Nineq)stop "ineq2ilat error: called with ineq > Nineq"
  end function ineq2isites




  function select_block(ip,Matrix) result(Vblock)
    integer                                                            :: ip
    complex(8),dimension(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb)       :: Matrix
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)              :: Vblock
    integer                                                            :: is,js,ispin,jspin,iorb,jorb,ilat,jlat
    Vblock=zero
    do ilat=1,Nlat
      do jlat=1,Nlat
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                    js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                    Vblock(is,js) = Matrix(ip,ilat,jlat,ispin,jspin,iorb,jorb)
                 enddo
              enddo
           enddo
        enddo
      enddo
    enddo
  end function select_block
  
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
   
   function ilso2nnnn_blocks(Hilso) result(Hnnnn)
      complex(8),dimension(Nsites*Nlat*Nspin*Norb,Nsites*Nlat*Nspin*Norb) :: Hilso
      complex(8),dimension(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: Hnnnn
      integer                                                             :: ilat,jlat,isites
      integer                                                             :: iorb,jorb
      integer                                                             :: ispin,jspin
      integer                                                             :: is,js
      Hnnnn=zero
      do isites=1,Nsites
        do ilat=1,Nlat
           do jlat=1,Nlat
              do ispin=1,Nspin
                 do jspin=1,Nspin
                    do iorb=1,Norb
                       do jorb=1,Norb
                          is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat + (isites-1)*Nlat*Nspin*Norb
                          js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat + (isites-1)*Nlat*Nspin*Norb
                          Hnnnn(isites,ilat,jlat,ispin,jspin,iorb,jorb) = Hilso(is,js)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
      enddo
   end function ilso2nnnn_blocks
   
   
   
   
  subroutine print_Hloc_ilso(hloc,file) ![Nlso][Nlso]
    complex(8),dimension(Nsites*Nlat*Nspin*Norb,Nsites*Nlat*Nspin*Norb) :: hloc
    character(len=*),optional                                           :: file
    integer                                                             :: ilat,is,js,unit
    unit=LOGfile
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    !
    write(unit,*)"Re(H)"
    do is=1,Nilso
       write(unit,"(20(F8.4,2x))")(real(Hloc(is,js)),js =1,Nilso)
    enddo
    write(unit,*)""
    write(unit,*)"Im(H)"
    do is=1,Nilso
       write(unit,"(20(F8.4,2x))")(imag(Hloc(is,js)),js =1,Nilso)
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hloc_ilso
  
   !-------------------------------------------------------------------------------------------
   !
   !
   ! FROM HERE START THE POSTPROCESSING ROUTINES
   !
   !
   !------------------------------------------------------------------------------------------- 
  
   !-------------------------------------------------------------------------------------------
   !PURPOSE: periodized Hamiltonian, for periodization schemes
   !-------------------------------------------------------------------------------------------

   
   function bhz_edge_model_periodized(kpoint,Nsites,N,pbc) result(hk)
    integer                                 :: Nsites,N,ii,Nlat_tmp,Nx_tmp
    real(8),dimension(:)                    :: kpoint
    complex(8),dimension(Nsites*N,Nsites*N) :: hk
    logical                                 :: pbc
    !
    !pbc=.false.
    !
    Nlat_tmp=Nlat
    Nx_tmp=Nx
    Nlat=1
    Nx=1
    !
    !
    Hk = bhz_edge_model(kpoint,Nsites,N,pbc)
    !
    Nlat=Nlat_tmp
    Nx=Nx_tmp
    !
   end function bhz_edge_model_periodized
  
   
   !-------------------------------------------------------------------------------------------
   !PURPOSE: periodization M scheme
   !-------------------------------------------------------------------------------------------
   
   function periodize_sigma_block_real(kpoint_,isite,ifreq,wprint) result(s_lso)
      integer                                                     :: isite,ilat,jlat,ispin,iorb,jorb,ii,ifreq
      real(8),dimension(1)                                        :: kpoint
      real(8),dimension(:)                                        :: kpoint_
      integer,dimension(:),allocatable                            :: ind1,ind2
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: tmpmat
      complex(8),dimension(Nspin,Nspin,Norb,Norb)           :: sreal_periodized,greal_periodized
      complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: g_lso,s_lso
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: greal_unperiodized
      logical,optional                                            :: wprint
      logical                                                     :: wprint_,pbc
      !
      real(8)                                                     :: Mh_tmp,ts_tmp,lambda_tmp
      !
      wprint_=.true.;if(present(wprint))wprint_=wprint
      pbc=.false.
      !
      kpoint(1)=kpoint_(1)
      !
      Mh_tmp=Mh
      ts_tmp=ts
      lambda_tmp=lambda
      !
      if(.not.allocated(ind1))allocate(ind1(size(kpoint)))
      if(.not.allocated(ind2))allocate(ind2(size(kpoint)))
      !
      sreal_periodized=zero
      greal_periodized=zero
      !
      ts=0d0
      lambda=0d0
      !
     tmpmat=(dcmplx(wr(ifreq),eps)+xmu)*eye(Nlat*Nspin*Norb) - bhz_edge_model(kpoint,1,Nlat*Nspin*Norb,pbc) - nnn2lso(Sreal_all(isite,:,:,:,:,:,:,ifreq))
     call inv(tmpmat)
     greal_unperiodized(:,:,:,:,:,:)=lso2nnn(tmpmat)
     !
     ts=ts_tmp
     lambda=lambda_tmp
     !
     do ilat=1,Nlat
        ind1=ilat        
        do jlat=1,Nlat
           ind2=jlat
           greal_periodized(:,:,:,:)=greal_periodized(:,:,:,:)+exp(-xi*dot_product(kpoint,ind1-ind2))*greal_unperiodized(ilat,jlat,:,:,:,:)/Nlat
        enddo
     enddo
     !
     g_lso(:,:)=nn2so(greal_periodized(:,:,:,:))
     call inv(g_lso(:,:))
     Mh=0
     g_lso(:,:)=g_lso(:,:)-bhz_edge_model_periodized(kpoint,1,Nspin*Norb,pbc)
     Mh=Mh_tmp
     s_lso(:,:)=(dcmplx(wr(ifreq),eps)+xmu)*eye(Nspin*Norb)-bhz_edge_model_periodized(kpoint,1,Nspin*Norb,pbc)-g_lso(:,:)
     !   
   end function periodize_sigma_block_real
   
   !-------------------------------------------------------------------------------------------
   ! PURPOSE:Get determinant of G(k,w)
   !-------------------------------------------------------------------------------------------   
  
  subroutine get_Akw()
    integer                                       :: Lk,Nso,Niso,Npts
    integer                                       :: ik,iw
    integer                                       :: iorb,jorb,idmax,idmin
    integer                                       :: ispin,jspin
    complex(8),dimension(:,:),allocatable         :: Sigmareal
    real(8),dimension(:,:),allocatable            :: Akreal,kpoints
    complex(8),dimension(:,:),allocatable         :: Gkreal
    complex(8),allocatable,dimension(:,:,:)       :: Hk_bare
    real(8),dimension(:,:),allocatable            :: Kpath
    character(len=30)                             :: suffix
    !
    !
    !
    allocate(kpath(3,2))
    !
    kpath(1,:)=[0.0,0.0]!G-e<-R
    kpath(2,:)=[1.0,0.0]!G
    kpath(3,:)=[2.0,0.0]!G+e->R
    !
    kpath=kpath*pi
    Npts  = size(kpath,1)
    Lk = (Npts-1)*Nkpath
    !   
    Nso=Nspin*Norb
    Niso=Nsites*Nspin*Norb
    !
    !
    if(allocated(Hk_bare))deallocate(Hk_bare)
    allocate(Hk_bare(Niso,Niso,Lk));Hk_bare=zero
    !
    !
    call TB_build_model(hk_bare,bhz_edge_model_periodized,Nsites,Nso,kpath,Nkpath,pbc=.false.)
    !
    !    
    allocate(kpoints(Lk,2))
    call TB_build_kgrid(kpath,Nkpath,kpoints)
    !
    !
    if(allocated(Sigmareal))deallocate(Sigmareal)
    if(allocated(Gkreal))deallocate(Gkreal)
    allocate(Sigmareal(Niso,Niso))
    allocate(Gkreal(Niso,Niso))
    Gkreal=zero      
    Sigmareal=zero
    allocate(Akreal(Lk,Lreal));Akreal=zero
    !
    !
    do ik=1,Lk
     do iw=1,Lreal
      do isites=1,Nsites
        idmin=1+(isites-1)*Nso
        idmax=isites*Nso
        Sigmareal(idmin:idmax,idmin:idmax)=periodize_sigma_block_real(kpoints(ik,:),isites,iw,wprint=.false.)
      enddo
      !
      Gkreal(:,:)=(dcmplx(wr(iw),0d0)+xmu)*eye(Niso) - Hk_bare(:,:,ik) - Sigmareal(:,:)
      call inv(Gkreal(:,:))
      Akreal(ik,iw) = log(abs(det(Gkreal(:,:)))/pi/Niso)
     enddo
    enddo

    call splot3d("Akw_real_nso.dat",kpoints(:,1),wr,Akreal) 

  end subroutine get_Akw
  
  
  

end program cdn_bhz_postprocessing_edge
