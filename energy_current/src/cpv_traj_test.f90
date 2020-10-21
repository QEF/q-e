program cpv_test
   use cpv_traj, only: cpv_trajectory, cpv_trajectory_initialize, cpv_trajectory_deallocate, &
                       cpv_trajectory_read_step, cpv_trajectory_get_step
   use traj_object, only: timestep ! type for timestep data
   use kinds, only: dp

   implicit none

   type(cpv_trajectory) :: test
   type(timestep) :: ts
   character(len=256) :: fpath = ''
   integer :: natoms, step_idx, iatom, i
   read (*, *) natoms, fpath
   call cpv_trajectory_initialize(test, fpath, natoms, 1.0_dp, 1.0_dp, 1.0_dp)
   step_idx = 0
   do while (cpv_trajectory_read_step(test))
      step_idx = step_idx + 1
      call cpv_trajectory_get_step(test, step_idx, ts)
      write (*, *) ts%nstep, ts%tps
      do iatom = 1, natoms
         write (*, *) iatom, (ts%tau(i, iatom), i=1, 3), (ts%vel(i, iatom), i=1, 3)
      end do
   end do

   call cpv_trajectory_deallocate(test)

   write (*, *) ''
   write (*, *) 'total steps', step_idx

end program
