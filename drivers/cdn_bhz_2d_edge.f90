program cdn_bhz_2d
   USE CDMFT_ED !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nx,Nlso,Nilso,iloop,Nb,Nkx,iw,Ly
   integer                                                                :: Nineq,Nsites,ineq,isites,icounter,icounter2
   integer,dimension(2)                                                   :: recover
   logical                                                                :: converged,lrsym
   real(8)                                                                :: ts,Mh,lambda,wmixing,observable
   !Bath:
   real(8),allocatable                                                    :: Bath(:,:),Bath_prev(:,:)
   !The local hybridization function:
   complex(8),allocatable                                                 :: Hloc(:,:)
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Hloc_nnnn,Hloc_ineq
   complex(8),allocatable,dimension(:,:,:,:,:,:,:,:)                      :: Gmats,Greal,Smats,Sreal
   complex(8),allocatable,dimension(:,:,:,:,:,:,:,:)                      :: Gmats_ineq,Smats_ineq,Sreal_ineq,Weiss_ineq,Weiss_old_ineq
   character(len=16)                                                      :: finput
   complex(8),allocatable                                                 :: wm(:),wr(:)
   complex(8),allocatable                                                 :: Hk(:,:,:),Smats_lso(:,:,:)
   !SYMMETRIES TEST
   real(8),dimension(:,:,:),allocatable                                   :: lambdasym_vectors
   complex(8),dimension(:,:,:,:,:,:,:),allocatable                        :: Hsym_basis
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
   call parse_input_variable(Ly,"Ly",finput,default=2,comment="Number of sites in the finite y direction")
   call parse_input_variable(lrsym,"LRSYM",finput,default=.true.)
   call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
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
   !if (Nspin/=1.or.Norb/=1) stop "You are using too many spin-orbitals"
   
   Nsites = Ly
   Nineq= Ly
   if(lrsym)then
     if(mod(Ly,2)/=0)stop "Wrong setup from input file: Ly%2 > 0 (odd number of sites)"
     Nineq=Ly/2
     print*,"Using L-R Symmetry. Solve",Nineq," of",Nsites," sites."
     call sleep(2)
   endif
   Nlat=Nx
   Nlso=Nlat*Nspin*Norb
   Nilso=Nsites*Nlat*Nspin*Norb
   
   if(.not.allocated(wm))allocate(wm(Lmats))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)

   !Allocate Weiss Field:
   allocate(Gmats(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats_lso(Nilso,Nilso,Lmats))
   
   allocate(Weiss_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss_ineq=zero
   allocate(Weiss_old_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss_old_ineq=zero
   allocate(Smats_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats_ineq=zero
   allocate(Sreal_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal_ineq=zero
   allocate(Gmats_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats_ineq=zero
   allocate(Hloc_nnnn(Nsites,Nlat,Nlat,Nspin,Nspin,Norb,Norb));Hloc_nnnn=zero
   allocate(Hloc_ineq(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero

   !Build Hk and Hloc
   call generate_hk_hloc()
   if(ed_verbose>2)call print_hloc_ilso(Hloc)
   Hloc_nnnn = ilso2nnnn_blocks(Hloc)
   do ineq=1,Nineq
     isites = ineq2isites(ineq)
     Hloc_ineq(ineq,:,:,:,:,:,:) = Hloc_nnnn(isites,:,:,:,:,:,:)
   enddo


   !SETUP BATH STEP 1
   allocate(lambdasym_vectors(Nineq,Nbath,3))
   allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,3))
   !
   do icounter=1,Nineq
    lambdasym_vectors(icounter,:,1)=Mh
   enddo
   Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(hloc_model(Nlso,1.d0,0.d0,0.d0))
   !
   do icounter=1,Nineq
    lambdasym_vectors(icounter,:,2)=ts
   enddo
   Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(hloc_model(Nlso,0.d0,1.d0,0.d0))
   !
   do icounter=1,Nineq
    lambdasym_vectors(icounter,:,3)=lambda
   enddo
   Hsym_basis(:,:,:,:,:,:,3)=lso2nnn(hloc_model(Nlso,0.d0,0.d0,1.d0))
   !
   !SETUP BATH STEP 2 and SETUP SOLVER
   call ed_set_Hreplica(Hsym_basis,lambdasym_vectors)
   Nb=ed_get_bath_dimension(Hsym_basis)
   allocate(bath(Nineq,Nb))
   allocate(bath_prev(Nineq,Nb))
   call ed_init_solver(comm,bath)

   !DMFT loop
   iloop=0;converged=.false.
   do while(.not.converged.AND.iloop<nloop)
      iloop=iloop+1
      if(master)call start_loop(iloop,nloop,"DMFT-loop")

      !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
      call ed_solve(comm,bath,Hloc_ineq)
      call ed_get_sigma_matsubara(Smats_ineq,Nineq)
      call ed_get_sigma_realaxis(Sreal_ineq,Nineq)
      !
      do isites=1,Nsites
        ineq = isites2ineq(isites)
        Smats(isites,:,:,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:,:,:)
      enddo
      !
      !Compute the local gfs:
      call dmft_gloc_matsubara(Hk,Gmats,Smats)
      !if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
      do ineq=1,Nineq
        isites = ineq2isites(ineq)
        Gmats_ineq(ineq,:,:,:,:,:,:,:) = Gmats(isites,:,:,:,:,:,:,:)
      enddo
      !
      !Get the Weiss field/Delta function to be fitted
      call dmft_self_consistency(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,cg_scheme)
      call Bcast_MPI(comm,Weiss_ineq)
      !
      !
      !
      if(master)then
         call ed_chi2_fitgf(Weiss_ineq,bath)
         !
         if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
         Bath_Prev=Bath
         !
         converged = check_convergence(Weiss_ineq(:,:,:,1,1,1,1,:),dmft_error,nsuccess,nloop)
         !
      endif
      !
      call Bcast_MPI(comm,bath)
      call Bcast_MPI(comm,converged)
      call Bcast_MPI(comm,xmu)
      !
      if(master)call end_loop
   enddo

   !Compute the local gfs:
   call dmft_gloc_realaxis(Hk,Greal,Sreal)
   !if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)


   call finalize_MPI()


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
    do i=1,Nsites-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 +     i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
    enddo
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
    if (N/=Nlso) stop "t0_rk_bhz: wrong dimension, the block has to be Nlso"
    H=zero
    do i=1,N/2
       Idmin=2*i-1
       Idmax=2*i
       H(Idmin:Idmax,Idmin:Idmax)=t_y(ts,lambda)
    enddo
  end function t0_rk_bhz



   !+------------------------------------------------------------------+
   !PURPOSE  : Hcluster for the 2d BHZ model
   !+------------------------------------------------------------------+


   function Hloc_model(N,Mh_,ts_,lambda_) result (H0)
      integer                                               :: N,ilat,jlat,ispin,iorb,jorb,ind1,ind2
      real(8)                                               :: Mh_,ts_,lambda_
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
      complex(8),dimension(N,N)                             :: H0
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do ilat=1,Nx
           ind1=ilat
           hopping_matrix(ind1,ind1,ispin,ispin,:,:)= t_m(mh_)
           if(ilat<Nx)then
              ind2=ilat+1
              hopping_matrix(ind2,ind1,ispin,ispin,:,:)= t_x(ts_,lambda_,ispin)
           endif
           if(ilat>1)then
              ind2=ilat-1
              hopping_matrix(ind2,ind1,ispin,ispin,:,:)= dconjg(transpose(t_x(ts_,lambda_,ispin)))
           endif
         enddo
      enddo
      !
      H0=nnn2lso(hopping_matrix)
      !
   end function hloc_model


   function hk_model(kpoint,N) result(Hk)
      integer                                                      :: N,ilat,jlat,ispin,iorb,jorb,i,j,ind1,ind2
      real(8)                                                      :: kpoint
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
      complex(8),dimension(N,N)                                    :: hk
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
        ind1=1
        ind2=Nx
        hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + dconjg(transpose(t_x(ts,lambda,ispin)))*exp(xi*kpoint*Nx)
        hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + t_x(ts,lambda,ispin)*exp(-xi*kpoint*Nx)
      enddo
      !
      Hk=nnn2lso(hopping_matrix)+hloc_model(N,Mh,ts,lambda)
      !
   end function hk_model


   !AUXILLIARY HOPPING MATRIX CONSTRUCTORS

   function t_m(mass) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: mass
      !
      tmpmat=zero
      tmpmat=mass*pauli_sigma_z
      !
   end function t_m

   function t_x(hop1,hop2,spinsign) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: hop1,hop2,sz
      integer                         :: spinsign
      !
      tmpmat=zero
      sz=(-1.d0)**(spinsign+1)
      tmpmat=-hop1*pauli_sigma_z+0.5d0*sz*xi*hop2*pauli_sigma_x
      !
   end function t_x

   function t_y(hop1,hop2) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: hop1,hop2
      !
      tmpmat=zero
      tmpmat=-hop1*pauli_sigma_z
      tmpmat(1,2)=-hop2*0.5d0
      tmpmat(2,1)=hop2*0.5d0
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
   !
end program cdn_bhz_2d








