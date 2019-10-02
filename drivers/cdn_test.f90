program ed_hm_1dchain
   USE CDMFT_ED
   !
   USE SCIFOR
   USE DMFT_TOOLS
   !
   USE MPI
   !
   implicit none
   integer                                                                :: Nx
   integer                                                                :: iloop,Nb,Nkx,iw
   logical                                                                :: converged
   real(8)                                                                :: ts,tsp,wmixing
   !Bath:
   real(8),allocatable                                                    :: Bath(:),BathOld(:)
   !The local hybridization function:
   complex(8),allocatable                                                 :: Hloc(:,:)
   complex(8),allocatable,dimension(:,:,:,:,:,:,:)                        :: Gmats,Greal,Smats,Sreal,Delta
   character(len=16)                                                      :: finput
   real(8),allocatable                                                    :: wt(:)
   complex(8),allocatable                                                 :: wm(:),wr(:)
   complex(8),allocatable                                                 :: Hk(:,:,:),Smats_lso(:,:,:)
   !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
   integer                                                                :: comm
   integer                                                                :: rank
   integer                                                                :: mpi_size
   logical                                                                :: master

   !Init MPI: use of MPI overloaded functions in SciFor
   call init_MPI(comm,.true.)
   rank   = get_Rank_MPI(comm)
   master = get_Master_MPI(comm)
   
   !Parse input variables
   call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
   call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
   call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
   call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
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
   if (Nspin/=1.or.Norb/=1) stop "You are using too many spin-orbitals"
   Nlat=Nx**2
   Nlso=Nlat*Nspin*Norb
   if(.not.allocated(wm))allocate(wm(Lmats))
   wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)


   !Allocate Weiss Field:
   allocate(delta(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
   allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
   allocate(Smats_lso(Nlso,Nlso,Lmats))

   !Build Hk and Hloc
   call generate_hk_hloc()
   
   !setup solver
   Nb=get_bath_dimension(lso2nnn(hloc))
   allocate(bath(Nb))
   allocate(bathold(Nb))
   call ed_init_solver(comm,bath,lso2nnn(Hloc))


   !DMFT loop
   iloop=0;converged=.false.
   do while(.not.converged.AND.iloop<nloop)
      iloop=iloop+1
      if(master)call start_loop(iloop,nloop,"DMFT-loop")

      !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
      bath(Nb)=0.d0
      call ed_solve(comm,bath) 
      call ed_get_sigma_matsubara(Smats)
      call ed_get_sigma_realaxis(Sreal)


      !Compute the local gfs:
      call dmft_gloc_matsubara(comm,Hk,Wt,Gmats,Smats)
      !if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=3)

      !Get the Weiss field/Delta function to be fitted
      call dmft_self_consistency(comm,Gmats,Smats,Delta,lso2nnn(Hloc),cg_scheme)
      call Bcast_MPI(comm,Delta)


      !Perform the SELF-CONSISTENCY by fitting the new bath
      if(master)then
         call ed_chi2_fitgf(delta,bath)
         !
         !MIXING:
         if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*BathOld
         BathOld=Bath
         !
         !Check convergence (if required change chemical potential)
         converged = check_convergence(delta(1,1,1,1,1,1,:),dmft_error,nsuccess,nloop)
      endif
      !
      call Bcast_MPI(comm,bath)
      call Bcast_MPI(comm,converged)
      call Bcast_MPI(comm,xmu)
      !
      if(master)call end_loop
   enddo

   !Compute the local gfs:
   call dmft_gloc_realaxis(comm,Hk,Wt,Greal,Sreal)
   if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=3)

   !Compute the Kinetic Energy:
   do iw=1,Lmats
      Smats_lso(:,:,iw)=nnn2lso(Smats(:,:,:,:,:,:,iw))
   enddo
   call dmft_kinetic_energy(comm,Hk(:,:,:),Wt,Smats_lso)


   call finalize_MPI()


contains

   !-------------------------------------------------------------------------------------------
   !PURPOSE:  Hk model for the 1d Hubbard chain
   !-------------------------------------------------------------------------------------------
   function hk_model(kpoint,N) result(hk)
      integer                                               :: N
      real(8),dimension(:)                                  :: kpoint
      complex(8),dimension(N,N)                             :: Hk
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: t_k
      !
      t_k=zero
      !
      t_k(1,1,1,1,1,1)=-xmu
      t_k(2,2,1,1,1,1)=-xmu
      t_k(3,3,1,1,1,1)=-xmu
      t_k(4,4,1,1,1,1)=-xmu
      !
      t_k(1,2,1,1,1,1)= -ts
      t_k(2,1,1,1,1,1)= -ts
      !
      t_k(2,3,1,1,1,1)= -ts
      t_k(3,2,1,1,1,1)= -ts
      !
      t_k(3,4,1,1,1,1)= -ts
      t_k(4,3,1,1,1,1)= -ts
      !
      t_k(4,1,1,1,1,1)= -ts
      t_k(1,4,1,1,1,1)= -ts
      !
      t_k(1,4,1,1,1,1)=t_k(1,4,1,1,1,1)-ts*exp(xi*kpoint(2))
      t_k(4,1,1,1,1,1)=t_k(4,1,1,1,1,1)-ts*exp(-xi*kpoint(2))
      t_k(1,2,1,1,1,1)=t_k(1,2,1,1,1,1)-ts*exp(xi*kpoint(1))
      t_k(2,1,1,1,1,1)=t_k(2,1,1,1,1,1)-ts*exp(-xi*kpoint(1))
      !
      t_k(2,3,1,1,1,1)=t_k(2,3,1,1,1,1)-ts*exp(-xi*kpoint(2))
      t_k(3,2,1,1,1,1)=t_k(3,2,1,1,1,1)-ts*exp(xi*kpoint(2))
      t_k(3,4,1,1,1,1)=t_k(3,4,1,1,1,1)-ts*exp(xi*kpoint(1))
      t_k(4,3,1,1,1,1)=t_k(4,3,1,1,1,1)-ts*exp(-xi*kpoint(1))
      !
      Hk=nnn2lso(t_k)
      !
   end function hk_model


   !-------------------------------------------------------------------------------------------
   !PURPOSE: generate Hloc and Hk
   !-------------------------------------------------------------------------------------------

   subroutine generate_hk_hloc()
      if(allocated(hk))deallocate(hk)
      if(allocated(Wt))deallocate(Wt)
      if(allocated(hloc))deallocate(Hloc)
      !
      allocate(Hk(Nlso,Nlso,Nkx**2),Wt(Nkx**2),Hloc(Nlso,Nlso))
      hk=zero
      wt=zero
      hloc=zero
      !
      call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
      call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nkx])
      Wt = 1d0/Nkx**2
      Hloc   = zero
      Hloc=sum(Hk(:,:,:),dim=3)/Nkx**2
      !
   end subroutine generate_hk_hloc

   !-------------------------------------------------------------------------------------------
   !PURPOSE: auxilliary reshape functions
   !-------------------------------------------------------------------------------------------

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

end program ed_hm_1dchain


