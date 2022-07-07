program cdn_sg77
  USE SCIFOR
  USE DMFT_TOOLS
  USE CDMFT_ED
  USE MPI
  !
  implicit none
  !
  integer                                         :: Nx,Nlso,iloop,Nb,Nk
  logical                                         :: converged
  real(8)                                         :: ts,wmixing
  real(8),allocatable                             :: Bath(:),Bath_prev(:)
  complex(8),allocatable                          :: Hloc(:,:)
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal,Smats,Sreal,Weiss,Weiss_old
  character(len=16)                               :: finput
  complex(8),allocatable                          :: wm(:),wr(:)
  complex(8),allocatable                          :: Hk(:,:,:)
  real(8),dimension(:),allocatable                :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Hsym_basis
  !
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer                                         :: comm
  integer                                         :: rank
  integer                                         :: mpi_size
  logical                                         :: master
  !
  !Init MPI: use of MPI overloaded functions in SciFor
  call init_MPI(comm,.true.)
  rank   = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !
  !Parse input variables
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=1d0,comment="hopping parameter")
  call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
  call parse_input_variable(Nk,"Nk",finput,default=11,comment="Number of k point for BZ integration")
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
  !
  !set global variables
  Nlat=Nx
  Nlso=Nlat*Nspin*Norb
  if(.not.allocated(wm))allocate(wm(Lmats))
  wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)
  !
  !Allocate Weiss Field:
  allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  
  !Build Hk and Hloc
  call print_hk()
  call generate_hk_hloc()
  
  
  !SETUP BATH STEP 1
  allocate(lambdasym_vector(1))
  allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,1))
  !
  lambdasym_vector(1)=ts
  Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(hloc_model(Nlso,1.d0))
  !
  !SETUP BATH STEP 2 and SETUP SOLVER
  call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
  Nb=ed_get_bath_dimension(Hsym_basis)
  allocate(bath(Nb))
  allocate(bath_prev(Nb))
  call ed_init_solver(comm,bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath,lso2nnn(Hloc))
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)

     !Compute the local gfs:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     !
     !Get the Weiss field/Delta function to be fitted
     call dmft_self_consistency(Gmats,Smats,Weiss,lso2nnn(Hloc),cg_scheme)
     call Bcast_MPI(comm,Weiss)
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(master)then
        call ed_chi2_fitgf(Weiss,bath)
        !
        !MIXING:
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
        Bath_Prev=Bath
        !
        !Check convergence (if required change chemical potential)
        converged = check_convergence(Weiss(:,:,1,1,1,1,:),dmft_error,nsuccess,nloop)
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
  if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)

  call finalize_MPI()

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : Hloc and H(k) for the 2d BHZ model
  !+------------------------------------------------------------------+

  function Hloc_model(N,ts_) result (H0)
    integer                                               :: N,ispin
    real(8)                                               :: ts_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
    complex(8),dimension(N,N)                             :: H0
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
        if(Nx>1)then
          !Hop 1 (4 terms are k-dependent)
          hopping_matrix(Nx,1,ispin,ispin,1,1) = ts_/2d0
          hopping_matrix(1,Nx,ispin,ispin,1,1) = ts_/2d0
          hopping_matrix(Nx,1,ispin,ispin,2,2) = -ts_/2d0
          hopping_matrix(1,Nx,ispin,ispin,2,2) = -ts_/2d0
          !Hop 5 (4-terms are k-dependent)
          hopping_matrix(Nx,1,ispin,ispin,1,2) = ts_/4d0
          hopping_matrix(1,Nx,ispin,ispin,1,2) = ts_/4d0
          hopping_matrix(Nx,1,ispin,ispin,2,1) = ts_/4d0
          hopping_matrix(1,Nx,ispin,ispin,2,1) = ts_/4d0
        endif
        !Hop 8 (only 4 terms because on-site)
        hopping_matrix(1,1,ispin,ispin,1,2) = ts_
        hopping_matrix(Nx,Nx,ispin,ispin,1,2) = ts_
        hopping_matrix(1,1,ispin,ispin,2,1) = ts_
        hopping_matrix(Nx,Nx,ispin,ispin,2,1) = ts_
    enddo
    !
    H0=nnn2lso(hopping_matrix)
    !
  end function hloc_model

  function hk_model(kpoint,N) result(Hk)
    integer                                                      :: N,ispin
    real(8),dimension(:)                                         :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
    complex(8),dimension(N,N)                                    :: hk,hktest
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
    !Hop 1 (other four terms are in Hloc)
      hopping_matrix(Nx,1,ispin,ispin,1,1) = hopping_matrix(Nx,1,ispin,ispin,1,1) + (ts/2d0)*exp(xi*kpoint(1)*Nx)
      hopping_matrix(1,Nx,ispin,ispin,1,1) = hopping_matrix(1,Nx,ispin,ispin,1,1) + (ts/2d0)*exp(-xi*kpoint(1)*Nx)
      hopping_matrix(Nx,1,ispin,ispin,2,2) = hopping_matrix(Nx,1,ispin,ispin,2,2) - (ts/2d0)*exp(xi*kpoint(1)*Nx)
      hopping_matrix(1,Nx,ispin,ispin,2,2) = hopping_matrix(1,Nx,ispin,ispin,2,2) - (ts/2d0)*exp(-xi*kpoint(1)*Nx)
    !Hop 2 (cosine)
      hopping_matrix(1,1,ispin,ispin,1,1)   = hopping_matrix(1,1,ispin,ispin,1,1) - (ts/2d0)*2*cos(kpoint(2))
      hopping_matrix(Nx,Nx,ispin,ispin,1,1) = hopping_matrix(Nx,Nx,ispin,ispin,1,1) - (ts/2d0)*2*cos(kpoint(2))
      hopping_matrix(1,1,ispin,ispin,2,2)   = hopping_matrix(1,1,ispin,ispin,2,2) + (ts/2d0)*2*cos(kpoint(2))
      hopping_matrix(Nx,Nx,ispin,ispin,2,2) = hopping_matrix(Nx,Nx,ispin,ispin,2,2) + (ts/2d0)*2*cos(kpoint(2))
    !Hop 3
      hopping_matrix(Nx,1,ispin,ispin,1,1) = hopping_matrix(Nx,1,ispin,ispin,1,1) - (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,1d0,0d0]))
      hopping_matrix(Nx,1,ispin,ispin,1,1) = hopping_matrix(Nx,1,ispin,ispin,1,1) - (ts/4d0)*exp(-xi*dot_product(kpoint,[-Nx*1d0,-1d0,0d0]))
      hopping_matrix(1,Nx,ispin,ispin,1,1) = hopping_matrix(1,Nx,ispin,ispin,1,1) - (ts/4d0)*exp(-xi*dot_product(kpoint,[Nx*1d0,1d0,0d0]))
      hopping_matrix(1,Nx,ispin,ispin,1,1) = hopping_matrix(1,Nx,ispin,ispin,1,1) - (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,-1d0,0d0]))
      hopping_matrix(Nx,1,ispin,ispin,2,2) = hopping_matrix(Nx,1,ispin,ispin,2,2) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,1d0,0d0]))
      hopping_matrix(Nx,1,ispin,ispin,2,2) = hopping_matrix(Nx,1,ispin,ispin,2,2) + (ts/4d0)*exp(-xi*dot_product(kpoint,[-Nx*1d0,-1d0,0d0]))
      hopping_matrix(1,Nx,ispin,ispin,2,2) = hopping_matrix(1,Nx,ispin,ispin,2,2) + (ts/4d0)*exp(-xi*dot_product(kpoint,[Nx*1d0,1d0,0d0]))
      hopping_matrix(1,Nx,ispin,ispin,2,2) = hopping_matrix(1,Nx,ispin,ispin,2,2) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,-1d0,0d0]))
    !Hop 4
      hopping_matrix(Nx,1,ispin,ispin,1,1) = hopping_matrix(Nx,1,ispin,ispin,1,1) + (ts/4d0)*exp(-xi*dot_product(kpoint,[-Nx*1d0,1d0,0d0]))
      hopping_matrix(Nx,1,ispin,ispin,1,1) = hopping_matrix(Nx,1,ispin,ispin,1,1) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,-1d0,0d0]))
      hopping_matrix(1,Nx,ispin,ispin,1,1) = hopping_matrix(1,Nx,ispin,ispin,1,1) + (ts/4d0)*exp(-xi*dot_product(kpoint,[Nx*1d0,-1d0,0d0]))
      hopping_matrix(1,Nx,ispin,ispin,1,1) = hopping_matrix(1,Nx,ispin,ispin,1,1) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,1d0,0d0]))
      hopping_matrix(Nx,1,ispin,ispin,2,2) = hopping_matrix(Nx,1,ispin,ispin,2,2) - (ts/4d0)*exp(-xi*dot_product(kpoint,[-Nx*1d0,1d0,0d0]))
      hopping_matrix(Nx,1,ispin,ispin,2,2) = hopping_matrix(Nx,1,ispin,ispin,2,2) - (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,-1d0,0d0]))
      hopping_matrix(1,Nx,ispin,ispin,2,2) = hopping_matrix(1,Nx,ispin,ispin,2,2) - (ts/4d0)*exp(-xi*dot_product(kpoint,[Nx*1d0,-1d0,0d0]))
      hopping_matrix(1,Nx,ispin,ispin,2,2) = hopping_matrix(1,Nx,ispin,ispin,2,2) - (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,1d0,0d0]))
    !Hop 5 (other 4 terms are in Hloc)
      hopping_matrix(Nx,1,ispin,ispin,1,2) = hopping_matrix(Nx,1,ispin,ispin,1,2) + (ts/4d0)*exp(xi*kpoint(1)*Nx)
      hopping_matrix(1,Nx,ispin,ispin,1,2) = hopping_matrix(1,Nx,ispin,ispin,1,2) + (ts/4d0)*exp(-xi*kpoint(1)*Nx)   
      hopping_matrix(Nx,1,ispin,ispin,2,1) = hopping_matrix(Nx,1,ispin,ispin,2,1) + (ts/4d0)*exp(xi*kpoint(1)*Nx) 
      hopping_matrix(1,Nx,ispin,ispin,2,1) = hopping_matrix(1,Nx,ispin,ispin,2,1) + (ts/4d0)*exp(-xi*kpoint(1)*Nx)
    !Hop 6
      hopping_matrix(1,1,ispin,ispin,1,2) = hopping_matrix(1,1,ispin,ispin,1,2) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,-1d0,-1d0]))
      hopping_matrix(Nx,Nx,ispin,ispin,1,2) = hopping_matrix(Nx,Nx,ispin,ispin,1,2) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,-1d0,-1d0]))
      hopping_matrix(1,1,ispin,ispin,2,1) = hopping_matrix(1,1,ispin,ispin,2,1) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,+1d0,+1d0]))
      hopping_matrix(Nx,Nx,ispin,ispin,2,1) = hopping_matrix(Nx,Nx,ispin,ispin,2,1) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,+1d0,+1d0]))
    !Hop 7
      hopping_matrix(1,1,ispin,ispin,1,2) = hopping_matrix(1,1,ispin,ispin,1,2) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,+1d0,-1d0]))
      hopping_matrix(Nx,Nx,ispin,ispin,1,2) = hopping_matrix(Nx,Nx,ispin,ispin,1,2) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,+1d0,-1d0]))
      hopping_matrix(1,1,ispin,ispin,2,1) = hopping_matrix(1,1,ispin,ispin,2,1) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,-1d0,+1d0]))
      hopping_matrix(Nx,Nx,ispin,ispin,2,1) = hopping_matrix(Nx,Nx,ispin,ispin,2,1) + (ts/4d0)*exp(-xi*dot_product(kpoint,[0d0,-1d0,+1d0]))
    !Hop 8
      !local
    !Hop 9 (cosine)
      hopping_matrix(1,1,ispin,ispin,1,2)   = hopping_matrix(1,1,ispin,ispin,1,2) + ts*(exp(xi*kpoint(3)))
      hopping_matrix(Nx,Nx,ispin,ispin,1,2)   = hopping_matrix(Nx,Nx,ispin,ispin,1,2) + ts*(exp(xi*kpoint(3)))
      hopping_matrix(1,1,ispin,ispin,2,1)   = hopping_matrix(1,1,ispin,ispin,2,1) + ts*(exp(-xi*kpoint(3)))
      hopping_matrix(Nx,Nx,ispin,ispin,2,1)   = hopping_matrix(Nx,Nx,ispin,ispin,2,1) + ts*(exp(-xi*kpoint(3)))
    enddo
    !
    Hk=nnn2lso(hopping_matrix)+hloc_model(N,ts)
    !
    hktest=conjg(transpose(Hk))-Hk
    if(any(abs(Real(hktest)).ge.0.00001))STOP "Not Hermitian (Re)"
    if(any(abs(Imag(hktest)).ge.0.00001))STOP "Not Hermitian (Im)"
    !
  end function hk_model
  
  !-------------------------------------------------------------------------------------------
  !PURPOSE: generate Hloc and Hk
  !-------------------------------------------------------------------------------------------

  subroutine generate_hk_hloc()
    real(8),dimension(Nk**3,3)                  :: kgrid
    real(8),dimension(3)                        :: e1,e2,e3,bk1,bk2,bk3
    real(8)                                     :: bklen
    !
    e1 = [2d0, 0d0, 0d0]
    e2 = [0d0, 1d0, 0d0]
    e3 = [0d0, 0d0, 1d0]
    call TB_set_ei(eix=e1,eiy=e2,eiz=e3)
    bklen=2d0*pi
    bk1=bklen*[0.5d0, 0d0, 0d0]
    bk2=bklen*[0d0, 1d0, 0d0]
    bk2=bklen*[0d0, 0d0, 1d0]
    call TB_set_bk(bkx=bk1,bky=bk2,bkz=bk3)
    !
    call TB_build_kgrid([Nk,Nk,Nk],kgrid)
    kgrid(:,1)=kgrid(:,1)/Nx
    !
    if(allocated(hk))deallocate(hk)
    if(allocated(hloc))deallocate(Hloc)
    !
    allocate(Hk(Nlso,Nlso,Nk**3),Hloc(Nlso,Nlso))
    hk=zero
    hloc=zero
    !
    call TB_build_model(Hk,hk_model,Nlso,kgrid)
    Hloc=hloc_model(Nlso,ts)
    !
  end subroutine generate_hk_hloc
  
  
  subroutine print_hk()
      integer                                :: i,j,Lk
      integer                                :: Npts
      real(8),dimension(:,:),allocatable     :: kpath
      character(len=64)                      :: file
      !
      write(LOGfile,*)"Build H(k) BHZ along path"
      !
      Npts = 12
      Lk=(Npts-1)*500
      allocate(kpath(Npts,3))
      kpath(1,:)=[0,0,0]
      kpath(2,:)=[1,0,0]
      kpath(3,:)=[1,1,0]
      kpath(4,:)=[0,0,0]
      kpath(5,:)=[0,0,1]
      kpath(6,:)=[1,0,1]
      kpath(7,:)=[1,1,1]
      kpath(8,:)=[0,0,1]
      kpath(9,:)=[1,0,0]
      kpath(10,:)=[1,0,1]
      kpath(11,:)=[1,1,0]
      kpath(12,:)=[1,1,1]
      !
      kpath=pi*kpath
      kpath(:,1)=kpath(:,1)/Nx
      file="Eigenbands.ed"
      !
      if(allocated(Hk))deallocate(Hk)
      allocate(Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lk))
      !
      call TB_set_bk([pi,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
        call TB_Solve_model(hk_model,Nlat*Nspin*Norb,kpath,500,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=20) :: 'G', 'X', 'M', 'G', 'Z', 'R', 'A', 'Z', 'X', 'R', 'M', 'A'],&
         file=reg(file))
      !
  end subroutine print_hk

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

end program cdn_sg77
