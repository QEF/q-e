module traj_object
   use kinds, only : dp


   implicit none

   type timestep ! it never owns the pointer
       real(dp), pointer :: tps=>null()
       integer,  pointer  :: nstep=>null()
       real(dp), pointer :: tau(:,:)=>null()
       real(dp), pointer :: vel(:,:)=>null()
   end type

   type trajectory ! it always owns the pointers! (there are no other references around)
       integer :: nsteps=0,nsteps_max=0,natoms=0
       real(dp), pointer :: tau(:,:,:) => null(), vel(:,:,:) => null(), tps(:)=>null()
       integer,  pointer :: nstep(:)=>null()
   end type

   contains

   subroutine trajectory_deallocate(t)
       type(trajectory), intent(inout) :: t 

       if (associated(t%tau)) then
           deallocate(t%tau)
           deallocate(t%vel)
           deallocate(t%tps)
           deallocate(t%nstep)
       end if
   end subroutine

   subroutine trajectory_remove_back(t)
       type(trajectory), intent(inout) :: t
       if (t%nsteps > 0) t%nsteps = t%nsteps - 1
   end subroutine

   subroutine trajectory_reallocate_if_necessary(t)
       type(trajectory), intent(inout) :: t
       integer :: newsize
       if (t%nsteps_max == t%nsteps) then
           newsize=(t%nsteps_max*3)/2
           if (newsize < t%nsteps_max ) &
               call errore('trajectory_reallocate_if_necessary','!!overflow!!',1)
           call trajectory_allocate(t,t%natoms,newsize)
       end if
   end subroutine

   subroutine trajectory_push_back(t, natoms, tps, nstep, tau, vel)
       type(trajectory), intent(inout) :: t
       real(dp), intent(in) :: tps, tau(:,:), vel(:,:)
       integer,  intent(in) :: nstep, natoms

       integer :: newsize

       if (natoms /= t%natoms) &
           call errore('trajectory_push_back','trying to push back a frame with a different number of atoms',1)

       !check if we must expand the array
       call trajectory_reallocate_if_necessary(t)

       t%nsteps = t%nsteps + 1
       !copy stuff
       t%tps(t%nsteps)   = tps
       t%nstep(t%nsteps) = nstep
       t%tau(:,:,t%nsteps) = tau
       t%vel(:,:,t%nsteps) = vel
   end subroutine

   subroutine trajectory_get_temporary(t, natoms,tstep)
       type(trajectory), intent(inout) :: t
       type(timestep), intent(out) :: tstep
       integer, intent(in) :: natoms
       integer :: idx 
       if (natoms /= t%natoms) &
           call errore('trajectory_get_temporary','trying to get a temporary frame with a different number of atoms',1)
       
       !check if we must expand the array
       call trajectory_reallocate_if_necessary(t)
       idx=t%nsteps+1 ! 1 past the end of the data
       call trajectory_get(t,idx,tstep,.true.) 
   end subroutine

   subroutine trajectory_push_back_last_temporary(t)
       type(trajectory), intent(inout) :: t
       if (t%nsteps == t%nsteps_max) &
           call errore('trajectory_push_back_last_temporary','there is no last temporary!',1)
       t%nsteps = t%nsteps + 1
   end subroutine

   subroutine trajectory_get(t, idx, tstep, dont_check)
       type(trajectory), intent(in) :: t
       integer, intent(in) :: idx 
       type(timestep), intent(out) :: tstep
       logical, optional, intent(in) :: dont_check
       if (present(dont_check)) then
           if (.not. dont_check .and. &
               .not. (idx>0 .and. idx <= t%nsteps) ) &
                   call errore('get_step','idx out of range',1)
       else
           if (.not. (idx>0 .and. idx <= t%nsteps) ) &
               call errore('get_step','idx out of range',1)
       end if
       tstep%tps   => t%tps(idx)
       tstep%nstep => t%nstep(idx)
       tstep%tau   => t%tau(:,:,idx)
       tstep%vel   => t%vel(:,:,idx)
   end subroutine

   subroutine trajectory_allocate(t,natoms,nsteps_to_allocate)
       !reallocate (to expand the trajectory or shrink) a trajectory type
       !or initialize from scratch

       type(trajectory), intent(inout) :: t
       integer, intent(in) :: natoms,nsteps_to_allocate
       integer :: n

       ! arrays that may be stolen by t
       real(dp), allocatable, target :: tau(:,:,:), vel(:,:,:), tps(:)
       integer, allocatable, target :: nstep(:)                       


       if (associated(t%tau) .and. t%nsteps_max /= nsteps_to_allocate) then
           if (t%natoms /= natoms) &
               call errore('trajectory_allocate','trying to reallocate with different number of atoms',1)
           !reallocate stuff and copy
           n=min(nsteps_to_allocate,t%nsteps)

           allocate(tps(nsteps_to_allocate))
           tps(1:n)=t%tps(1:n)
           deallocate(t%tps)
           t%tps => tps

           allocate(tau(3,natoms,nsteps_to_allocate))
           tau(:,:,1:n)=t%tau(:,:,1:n)
           deallocate(t%tau)
           t%tau => tau

           allocate(vel(3,natoms,nsteps_to_allocate))
           vel(:,:,1:n)=t%vel(:,:,1:n)
           deallocate(t%vel)
           t%vel => vel

           allocate(nstep(nsteps_to_allocate))
           nstep(1:n)=t%nstep(1:n)
           deallocate(t%nstep)
           t%nstep => nstep
           
           t%nsteps=n
           t%nsteps_max=nsteps_to_allocate
           ! now t is extended (or reduced)
       else if (.not. associated(t%tau)) then
           ! initialize structure with nothing inside
           allocate(t%tau(3,natoms,nsteps_to_allocate))
           allocate(t%tps(nsteps_to_allocate))
           allocate(t%vel(3,natoms,nsteps_to_allocate))
           allocate(t%nstep(nsteps_to_allocate))
           t%nsteps=0
           t%natoms=natoms
           t%nsteps_max=nsteps_to_allocate
       else
           if (t%natoms /= natoms) &
               call errore('trajectory_allocate','trying to reallocate with different number of atoms',1)
           !since the dimension is the same, I have nothing to do
       end if

   end subroutine


end module
