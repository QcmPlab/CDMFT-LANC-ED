program cdn_bhz_2d
  USE CDMFT_ED !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  !
  implicit none
  integer                                         :: Nx,Ny,Nlso,iloop,Nb,Nkx,Nky,iw,iii,jjj,kkk,Ntr
  integer,dimension(2)                            :: recover
  logical                                         :: converged
  real(8)                                         :: ts,Mh,lambda,wmixing,observable
  !Bath:
  real(8),allocatable                             :: Bath(:),Bath_prev(:)
  !The local hybridization function:
  complex(8),allocatable                          :: Hloc(:,:)
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal,Smats,Sreal,Weiss,Weiss_old
  character(len=16)                               :: finput
  real(8),allocatable                             :: wt(:)
  complex(8),allocatable                          :: wm(:),wr(:)
  complex(8),allocatable                          :: Hk(:,:,:),Smats_lso(:,:,:)
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: observable_matrix
  !Density matrices:
  complex(8),allocatable,dimension(:,:)           :: reduced_density_matrix
  complex(8),allocatable,dimension(:,:)           :: pure_cvec
  real(8),allocatable,dimension(:)                :: pure_prob
  !SYMMETRIES TEST
  real(8),dimension(:),allocatable                :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Hsym_basis
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer                                         :: comm
  integer                                         :: rank
  integer                                         :: mpi_size
  logical                                         :: master,hermiticize,dmexit

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
  call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of cluster sites in y direction")
  call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
  call parse_input_variable(Nky,"Nky",finput,default=10,comment="Number of ku point for BZ integration")
  call parse_input_variable(hermiticize,"HERMITICIZE",finput,default=.true.,comment="are bath replicas hermitian")
  call parse_input_variable(dmexit,"DMEXIT",finput,default=.true.,comment="Exit after evaluating DM, behave as post-processing")
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
  !Ny=Nx
  !Nky=Nkx
  Nlat=Nx*Ny
  Nlso=Nlat*Nspin*Norb
  if(.not.allocated(wm))allocate(wm(Lmats))
  wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)

  if(ED_VERBOSE > 0)call naming_convention()
  !Allocate Weiss Field:
  allocate(Weiss(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats_lso(Nlso,Nlso,Lmats))

  !Build Hk and Hloc
  call generate_hk_hloc()

  !CUSTOM OBSERVABLES: n1, n2, Ekin
  allocate(observable_matrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
  call init_custom_observables(3,Hk)
  observable_matrix=zero
  do iii=1,Nlat
     observable_matrix(iii,iii,1,1,1,1)=one/Nlat
     observable_matrix(iii,iii,Nspin,Nspin,1,1)=one/Nlat
  enddo
  call add_custom_observable("n1",nnn2lso(observable_matrix))
  observable_matrix=zero
  do iii=1,Nlat
     observable_matrix(iii,iii,1,1,2,2)=one/Nlat
     observable_matrix(iii,iii,Nspin,Nspin,2,2)=one/Nlat
  enddo
  call add_custom_observable("n2",nnn2lso(observable_matrix))
  call add_custom_observable("Ekin",Hk)

  !SETUP BATH STEP 1
  allocate(lambdasym_vector(3))
  allocate(Hsym_basis(Nlat,Nlat,Nspin,Nspin,Norb,Norb,3))
  !
  lambdasym_vector(1)=Mh
  Hsym_basis(:,:,:,:,:,:,1)=lso2nnn(hloc_model(Nlso,1.d0,0.d0,0.d0))
  !
  lambdasym_vector(2)=ts
  Hsym_basis(:,:,:,:,:,:,2)=lso2nnn(hloc_model(Nlso,0.d0,1.d0,0.d0))
  !
  lambdasym_vector(3)=lambda
  Hsym_basis(:,:,:,:,:,:,3)=lso2nnn(hloc_model(Nlso,0.d0,0.d0,1.d0))
  !
  !SETUP BATH STEP 2 and SETUP SOLVER
  call ed_set_Hbath(Hsym_basis,lambdasym_vector)
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


     !Retrieve ALL REDUCED DENSITY MATRICES DOWN TO THE LOCAL one
     if(dm_flag.AND.master)then
        do Ntr=0,Nlat-1 ! Ntr: number of cluster sites we want to trace out
           call ed_get_reduced_dm(reduced_density_matrix,Nlat-Ntr,doprint=.true.)
           !
           !DIAGONALIZATION
           write(*,*) " "
           write(*,*) "Diagonalizing "//str(Nlat-Ntr)//"-site REDUCED DENSITY MATRIX"
           allocate(pure_prob(4**((Nlat-Ntr)*Norb)))
           pure_cvec = reduced_density_matrix
           call eigh(pure_cvec,pure_prob,jobz='V',uplo='U')
           !
           !PRINT-TO-FILE (and log)
           call print_pure_states(pure_cvec,pure_prob,Nlat-Ntr)
           !
           deallocate(pure_cvec,pure_prob,reduced_density_matrix)
        enddo
        !
        ! !Semi-analytical crosscheck of the local density matrix TO BE DONE
        ! if(Norb==1)call local_dm_benchmark() 
        !
        if(dmexit)stop "Got DM.. ciao"
     endif

     !Compute the local gfs:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     if(master)call dmft_write_gf(Gmats,"Gloc",axis="matsubara",iprint=4)
     !
     !Get the Weiss field/Delta function to be fitted
     call legacy_dmft_self_consistency(Gmats,Smats,Weiss,lso2nnn(Hloc),cg_scheme)
     !call dmft_self_consistency(Gmats,Smats,Weiss,lso2nnn(Hloc),cg_scheme)
     call Bcast_MPI(comm,Weiss)
     !
     !MIXING:
     !if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_Old
     !Weiss_old=Weiss
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(master)then
        call ed_chi2_fitgf(Weiss,bath)
        !
        !MIXING:
        !call adaptive_mix(Bath(Nbath+1:),Bath_fitted(Nbath+1:)-Bath(Nbath+1:),wmixing,iloop)
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
  if(master)call dmft_write_gf(Greal,"Gloc",axis="realaxis",iprint=4)

  !Compute the Kinetic Energy:
  do iw=1,Lmats
     Smats_lso(:,:,iw)=nnn2lso(Smats(:,:,:,:,:,:,iw))
  enddo
  call dmft_kinetic_energy(Hk(:,:,:),Smats_lso)


  call finalize_MPI()


contains

  !+------------------------------------------------------------------+
  !PURPOSE  : Hloc for the 2d BHZ model
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
          do jlat=1,Ny
             ind1=indices2N([ilat,jlat])
             hopping_matrix(ind1,ind1,ispin,ispin,:,:)= t_m(mh_)
             if(ilat<Nx)then
                ind2=indices2N([ilat+1,jlat])
                hopping_matrix(ind2,ind1,ispin,ispin,:,:)= t_x(ts_,lambda_,ispin)
             endif
             if(ilat>1)then
                ind2=indices2N([ilat-1,jlat])
                hopping_matrix(ind2,ind1,ispin,ispin,:,:)= dconjg(transpose(t_x(ts_,lambda_,ispin)))
             endif
             if(jlat<Ny)then
                ind2=indices2N([ilat,jlat+1])
                hopping_matrix(ind2,ind1,ispin,ispin,:,:)= t_y(ts_,lambda_)
             endif
             if(jlat>1)then
                ind2=indices2N([ilat,jlat-1])
                hopping_matrix(ind2,ind1,ispin,ispin,:,:)= transpose(t_y(ts_,lambda_))
             endif
          enddo
       enddo
    enddo
    !
    H0=nnn2lso(hopping_matrix)
    !
  end function hloc_model


  function hk_model(kpoint,N) result(Hk)
    integer                                                      :: N,ilat,jlat,ispin,iorb,jorb,i,j,ind1,ind2
    real(8),dimension(:)                                         :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
    complex(8),dimension(N,N)                                    :: hk
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
       do ilat=1,Ny
          ind1=indices2N([1,ilat])
          ind2=indices2N([Nx,ilat])
          hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + dconjg(transpose(t_x(ts,lambda,ispin)))*exp(xi*kpoint(1)*Nx)
          hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + t_x(ts,lambda,ispin)*exp(-xi*kpoint(1)*Nx)
       enddo
       do ilat =1,Nx
          ind1=indices2N([ilat,1])
          ind2=indices2N([ilat,Ny])
          hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + transpose(t_y(ts,lambda))*exp(xi*kpoint(2)*Ny)
          hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + t_y(ts,lambda)*exp(-xi*kpoint(2)*Ny)
       enddo
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
    real(8),dimension(Nkx*Nky,2)                :: kgrid
    real(8),dimension(2)                        :: e1,e2,bk1,bk2
    real(8)                                     :: bklen
    !
    e1 = [1d0, 0d0]
    e2 = [0d0, 1d0]
    call TB_set_ei(eix=e1,eiy=e2)
    bklen=2d0*pi
    bk1=bklen*[1d0, 0d0]
    bk2=bklen*[0d0, 1d0]
    call TB_set_bk(bkx=bk1,bky=bk2)
    !
    call TB_build_kgrid([Nkx,Nky],kgrid)
    kgrid(:,1)=kgrid(:,1)/Nx
    kgrid(:,2)=kgrid(:,2)/Ny
    !
    if(allocated(hk))deallocate(hk)
    if(allocated(Wt))deallocate(Wt)
    if(allocated(hloc))deallocate(Hloc)
    !
    allocate(Hk(Nlso,Nlso,Nkx*Nky),Wt(Nkx*Nky),Hloc(Nlso,Nlso))
    hk=zero
    wt=zero
    hloc=zero
    !
    ! SEVER !
    call TB_build_model(Hk,hk_model,Nlso,kgrid)
    Wt = 1d0/(Nkx*Nky)
    Hloc=hloc_model(Nlso,Mh,ts,lambda)
    !
  end subroutine generate_hk_hloc




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


  function indices2N(indices) result(N)
    integer,dimension(2)         :: indices
    integer                      :: N,i
    !
    !
    N=Nx*(indices(2)-1)+indices(1)
  end function indices2N

  function N2indices(N) result(indices)
    integer,dimension(2)         :: indices
    integer                      :: N,i
    !
    indices(1)=mod(N,Nx)
    if(indices(1)==0)then
       indices(1)=Nx
       indices(2)=(N-Nx)/Nx+1
    else
       indices(2)=N/Nx+1
    endif
  end function N2indices

  subroutine naming_convention()
    integer                       :: i,j
    integer,dimension(Nx,Ny)      :: matrix
    !
    do j=1,Ny
       do i=1,Nx
          matrix(i,j)=indices2N([i,j])
       enddo
    enddo
    !
    write(LOGfile,"(A)")"The unique index of each site (on the cartesian plane) is as follows:"
    write(LOGfile,"(A)")" "
    do j=1,Ny
       write(LOGfile,"(20(I2,2x))")(matrix(i,Ny+1-j),i =1,Nx)
    enddo
    write(LOGfile,"(A)")" "
  end subroutine naming_convention
  !







  ! !+---------------------------------------------------------------------------+
  ! !PURPOSE : check the local-dm comparing to Eq.4 in Mod.Phys.Lett.B.2013.27:05
  ! !+---------------------------------------------------------------------------+
  ! subroutine local_dm_benchmark()
  !   complex(8),allocatable,dimension(:,:)  :: local_density_matrix
  !   real(8),allocatable,dimension(:,:)     :: dens,dens_up,dens_dw,docc,mag
  !   !
  !   if(Norb>1) stop "ERROR: local_dm_benchmark available for Norb=1 only."
  !   !
  !   allocate(dens(Nlat,Norb),dens_up(Nlat,Norb),dens_dw(Nlat,Norb),docc(Nlat,Norb),mag(Nlat,Norb))
  !   !
  !   call ed_get_reduced_dm(local_density_matrix,1,doprint=.false.)
  !   call ed_get_mag(mag)
  !   call ed_get_dens(dens)
  !   call ed_get_docc(docc)
  !   dens_up = 0.5d0*(dens + mag)
  !   dens_dw = 0.5d0*(dens - mag)
  !   write(*,*)
  !   write(*,*) "LOCAL-DM BENCHMARK [Mod.Phys.Lett.B.2013.27:05]"
  !   write(*,*) "Semi-Analytical Estimate  |  Error"
  !   write(*,*) 1-dens_up(1,1)-dens_dw(1,1)+docc(1,1), "|", abs(1-dens_up(1,1)-dens_dw(1,1)+docc(1,1)-local_density_matrix(1,1))
  !   write(*,*) dens_up(1,1)-docc(1,1),                "|", abs(dens_up(1,1)-docc(1,1)-local_density_matrix(2,2))
  !   write(*,*) dens_dw(1,1)-docc(1,1),                "|", abs(dens_dw(1,1)-docc(1,1)-local_density_matrix(3,3))
  !   write(*,*) docc(1,1),                             "|", abs(docc(1,1)-local_density_matrix(4,4))
  !   write(*,*)
  !   !
  ! end subroutine local_dm_benchmark

  !+---------------------------------------------------------------------------+
  !PURPOSE : print to file (and stdout) the reduced-dm eigenstuff
  !+---------------------------------------------------------------------------+
  subroutine print_pure_states(pure_cvec,pure_prob,Nsites)
    complex(8),allocatable,dimension(:,:)     :: pure_cvec
    real(8),allocatable,dimension(:)          :: pure_prob
    integer                                   :: Nsites
    integer                                   :: unit      
    !
    if(ed_verbose>3) write(*,*) "> PURE-STATE probabilities:"
    open(free_unit(unit),file="probabilities_"//str(Nsites)//"sites"//".dat",action="write",position="rewind",status='unknown')
    do iii=1,4**(Nsites*Norb)
       if(ed_verbose>3) write(*,"(90(F15.5,1X))") abs(pure_prob(iii))
       write(unit,*) abs(pure_prob(iii)) !Machine zero could be negative.
    enddo
    close(unit)
    !
    if(ed_verbose>4) write(*,*) "> PURE-STATE components:"
    open(free_unit(unit),file="pure-states_"//str(Nsites)//"sites"//".dat",action="write",position="rewind",status='unknown')
    do iii=1,4**(Nsites*Norb)
       if(ed_verbose>4) write(*,"(180(F15.5,1X))") (pure_cvec(jjj,iii), jjj=1,4**(Nsites*Norb))
       write(unit,*) (pure_cvec(jjj,iii), jjj=1,4**(Nsites*Norb))
    enddo
    close(unit)
    write(*,*) " "
    !
  end subroutine print_pure_states





end program cdn_bhz_2d








