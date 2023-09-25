subroutine ed_get_reduced_density_matrix_single_LEGACY(rdm,Nsites,doprint)
   !! further reduce the cdm by tracing out (Nlat-Nsites) sites
   complex(8),dimension(4**Nimp,4**Nimp)                      :: cdm
   complex(8),dimension(:,:),allocatable,intent(out)          :: rdm
   integer                              ,intent(in)           :: Nsites
   logical                              ,intent(in),optional  :: doprint
   logical                                                    :: doprint_
   integer    :: i,j,io,jo,iUP,iDW,jUP,jDW
   integer    :: iIMPup,iIMPdw,jIMPup,jIMPdw
   integer    :: iREDup,iREDdw,jREDup,jREDdw
   integer    :: iTrUP,iTrDW,jTrUP,jTrDW,Nred
   !
   if(Nsites>Nlat)then
      stop "ERROR: cannot request a density matrix reduced to more sites then Nlat"
   endif
   !
   doprint_=.false.; if(present(doprint)) doprint_=doprint
   !
   !Retrieve cdm from the global scope
   if(.not.allocated(cluster_density_matrix))then
      stop "ERROR: cluster_density_matrix is not allocated"
   endif
   cdm = cluster_density_matrix
   !
   Nred=Norb*Nsites !Number of levels associated to the Nsites-reduced system
   allocate(rdm(4**Nred,4**Nred))
   rdm=zero
   !
   !Trace the cdm to the requested Nsites-rdm
   do iUP = 1,2**Nimp
      do iDW = 1,2**Nimp
         i = iUP + (iDW-1)*2**Nimp
         iIMPup = iup-1
         iIMPdw = idw-1
         iREDup = Ibits(iIMPup,0,Nred)
         iREDdw = Ibits(iIMPdw,0,Nred)
         iTrUP  = Ibits(iIMPup,Nred,Nimp)
         iTrDW  = Ibits(iIMPdw,Nred,Nimp)
         do jUP = 1,2**Nimp
            do jDW = 1,2**Nimp
               j = jUP + (jDW-1)*2**Nimp
               jIMPup = jup-1
               jIMPdw = jdw-1
               jREDup = Ibits(jIMPup,0,Nred)
               jREDdw = Ibits(jIMPdw,0,Nred)
               jTrUP  = Ibits(jIMPup,Nred,Nimp)
               jTrDW  = Ibits(jIMPdw,Nred,Nimp)
               if(jTrUP/=iTrUP.or.jTrDW/=iTrDW)cycle
               io = (iREDup+1) + iREDdw*2**Nred
               jo = (jREDup+1) + jREDdw*2**Nred
               rdm(io,jo) = rdm(io,jo) + cdm(i,j)
            enddo
         enddo
         ! NB: the spin-factorization of the loops is paramount,
         !     since we build in such way the original cdm and the
         !     ordering of the basis states has to be preserved
         !     in order to trace along the proper matrix elements
      enddo
   enddo
   !
   !Print to file (if requested)
   if(doprint_)then
      call ed_print_dm_LEGACY(rdm,4**Nred)
   endif
   !
end subroutine ed_get_reduced_density_matrix_single_LEGACY

subroutine ed_get_reduced_density_matrix_single(rdm,orbital_mask,doprint)
   !! further reduce the cdm by tracing out selected orbitals (via mask)
   complex(8),dimension(:,:),allocatable,intent(out)          :: rdm
   logical,dimension(Nlat,Norb),intent(in)                    :: orbital_mask
   logical,intent(in),optional                                :: doprint
   logical                                                    :: doprint_
   logical                                                    :: dotrace_
   integer,dimension(:),allocatable                           :: red_indices
   integer,dimension(:),allocatable                           :: trace_indices
   integer,dimension(:),allocatable :: IbUP,IbDW,JbUP,JbDW
   integer    :: i,j,io,jo,iUP,iDW,jUP,jDW,ilat,jorb,red_count,trace_count
   integer    :: iIMPup,iIMPdw,jIMPup,jIMPdw
   integer    :: iREDup,iREDdw,jREDup,jREDdw
   integer    :: iTrUP,iTrDW,jTrUP,jTrDW,Nred
   real(8)    :: IsignUP,IsignDW,JsignUP,JsignDW,sign

   ! Input handling and checks
   Nred = count(orbital_mask)
   if(Nred<1)then
      stop "ERROR: invalid orbital mask, the reduced system must consist of at least one orbital"
   elseif(Nred==Nimp)then
      dotrace_ = .FALSE.
   else
      dotrace_ = .TRUE.
      allocate(red_indices(Nred),trace_indices(Nimp-Nred))
   endif
   doprint_=.false.; if(present(doprint)) doprint_=doprint

   ! Retrieve cdm from the global scope
   if(.not.allocated(cluster_density_matrix))then
      stop "ERROR: cluster_density_matrix is not allocated"
   endif
   associate(cdm => cluster_density_matrix)
      
      if(.not.dotrace_)then

         ! Pour the whole cdm into the rdm
         rdm = cdm ! and that's all we do…

      else

         ! Retrieve the requested bit-indices for the reduced/traced system
         red_count = 0; trace_count = 0
         do ilat = 1,Nlat
            do jorb=1,Norb
               if(orbital_mask(ilat,jorb))then
                  red_count = red_count + 1
                  red_indices(red_count) = jorb + (ilat-1)*Nlat
               else
                  trace_count = trace_count + 1
                  trace_indices(trace_count) = jorb + (ilat-1)*Nlat
               endif
            enddo
         enddo

         allocate(rdm(4**Nred,4**Nred))
         rdm = zero

         ! Trace the cdm to the requested subsystem, and store into rdm
         do iUP = 1,2**Nimp
            IbUP  = bdecomp(iUP,Nimp)
            call get_sign(IsignUP,IbUP,red_indices)
            call split_state(IbUp,red_indices,trace_indices,iREDup,iTrUP)
            do iDW = 1,2**Nimp
               i = iUP + (iDW-1)*2**Nimp
               IbDW = bdecomp(iDW,Nimp)
               call get_sign(IsignDW,IbDW,red_indices)
               call split_state(IbDw,red_indices,trace_indices,iREDdw,iTrDW)
               do JUP = 1,2**Nimp
                  JbUP  = bdecomp(Jup,Nimp)
                  call get_sign(JsignUP,JbUP,red_indices)
                  call split_state(JbUp,red_indices,trace_indices,jREDup,jTrUP)
                  do jDW = 1,2**Nimp
                     j = jUP + (jDW-1)*2**Nimp
                     JbDW = bdecomp(jDW,Nimp)
                     call get_sign(JsignDW,JbDW,red_indices)
                     call split_state(JbDw,red_indices,trace_indices,jREDdw,jTrDW)
                     if(jTrUP/=iTrUP.or.jTrDW/=iTrDW)cycle
                     io = (iREDup+1) + iREDdw*2**Nred
                     jo = (jREDup+1) + jREDdw*2**Nred
                     sign = IsignUP * IsignDW * JsignUP * JsignDW
                     rdm(io,jo) = rdm(io,jo) + cdm(i,j) * sign
                  enddo
               enddo
               ! NB: the spin-factorization of the loops is paramount,
               !     since we build in such way the original cdm and the
               !     ordering of the basis states has to be preserved
               !     in order to trace along the proper matrix elements
            enddo
         enddo

      end if

   end associate

   !Print to file (if requested)
   if(doprint_)then
      call ed_print_dm(rdm,orbital_mask)
   endif

contains

   subroutine get_sign(sign,state,indices)
      !! Compute the Fermionic sign associated to the required swipes
      real(8), intent(out) :: sign
      integer, intent(in)  :: state(Nimp)
      integer, intent(in)  :: indices(:)
      integer :: filtered(Nimp)
      integer :: N
      integer :: r
      ! FILTER THE STATE TO CONSTRAIN THE SUM
      filtered = state; filtered(indices)=0
      ! PERFORM THE SUM (count permutations)
      N = 0
      do r=1,size(indices)
         N = N + sum(filtered(1:indices(r)))
      enddo
      ! ASSIGN THE SIGN: (-1)^N
      if(mod(N,2)==0)then
         sign = +1
      else
         sign = -1
      endif
   end subroutine get_sign

   subroutine split_state(state,reduced_indices,tracing_indices,reduced_state,tracing_state)
      !! Extract the reduced and tracing subsystems, given the appropriate bit-indices
      integer,intent(in),allocatable  :: state(:)
      integer,intent(in),allocatable  :: reduced_indices(:)
      integer,intent(in),allocatable  :: tracing_indices(:)
      integer,intent(out)             :: reduced_state
      integer,intent(out)             :: tracing_state
      !
      integer,allocatable             :: reduced_ibits(:)
      integer,allocatable             :: tracing_ibits(:)
      !
      reduced_ibits = state(reduced_indices)
      tracing_ibits = state(tracing_indices)
      !
      reduced_state = bjoin(reduced_ibits,Nred)
      tracing_state = bjoin(tracing_ibits,Nimp-Nred)
      !
   end subroutine split_state

end subroutine ed_get_reduced_density_matrix_single

! function ed_get_orbital_mask(site1,site2,site3,site4) result(omask)
!    !! Tentative ergonomic API to build the mask (which is needed by
!    !! both the getter and the printer routines, so let's at least
!    !! do all this mess only in one place, please)
!    logical,dimension(Nlat,Norb)                               :: omask
!    integer,dimension(:),intent(in),optional                   :: site1
!    integer,dimension(:),intent(in),optional                   :: site2
!    integer,dimension(:),intent(in),optional                   :: site3
!    integer,dimension(:),intent(in),optional                   :: site4
!    integer,allocatable,dimension(:)                           :: s1,s2,s3,s4
!    integer                                                    :: counter
!    !
!    ! INPUT CHECKS (so freaking ugly, but for a good cause (shiny interface...)...)
!    if(any(site1>Norb).or.size(site1)>Norb)stop "ERROR: invalid orbital indices for site1"
!    Nred = size(site1)
!    s1 = site1
!    if(present(site2))then
!       if(any(site2>Norb).or.size(site2)>Norb)stop "ERROR: invalid orbital indices for site2"
!       Nred = Nred + size(site2)
!       s2 = site2
!    else
!       allocate(s2(0))
!    endif
!    if(present(site3))then
!       if(any(site3>Norb).or.size(site3)>Norb)stop "ERROR: invalid orbital indices for site3"
!       Nred = Nred + size(site3)
!       s3 = site3
!    else
!       allocate(s3(0))
!    endif
!    if(present(site4))then
!       if(any(site4>Norb).or.size(site4)>Norb)stop "ERROR: invalid orbital indices for site4"
!       Nred = Nred + size(site4)
!       s4 = site4
!    else
!       allocate(s4(0))
!    endif
!    ! BUILD THE LOGICAL MASK
!    counter = 0
!    red_indices(counter+1:counter+size(s1)) = s1 + 0*Norb
!    counter = counter + size(s1)
!    red_indices(counter+1:counter+size(s2)) = s2 + 1*Norb
!    counter = counter + size(s2)
!    red_indices(counter+1:counter+size(s3)) = s3 + 2*Norb
!    counter = counter + size(s3)
!    red_indices(counter+1:counter+size(s4)) = s4 + 3*Norb
! end function ed_get_orbital_mask
