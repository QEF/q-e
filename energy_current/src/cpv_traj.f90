!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
module cpv_traj
   use kinds, only: dp
   use traj_object ! timestep,trajectory and all subroutines that starts with trajectory_

   implicit none

   type cpv_trajectory
      type(trajectory) :: traj
      character(len=256) :: fname = ''
      integer :: iounit_pos = -1, iounit_vel = -1
      logical :: is_open = .false.
      real(dp) :: tau_fac = 1.0_dp, vel_fac = 1.0_dp, tps_fac = 1.0_dp
   end type

contains

   subroutine cpv_trajectory_initialize(t, fname, natoms, tau_fac, vel_fac, tps_fac, circular, ios)
      implicit none
      type(cpv_trajectory), intent(inout) :: t
      character(len=256), intent(in) :: fname
      integer, intent(in) :: natoms
      real(dp), intent(in) :: tau_fac, vel_fac, tps_fac
      logical, intent(in), optional :: circular
      integer, intent(out), optional :: ios
      integer :: iostat
      integer, external :: find_free_unit
      ! try to open fname, allocate traj
      t%iounit_pos = find_free_unit()
      open (unit=t%iounit_pos, file=trim(fname)//'.pos', iostat=iostat, action='read')
      if (.not. present(ios) .and. iostat /= 0) &
         call errore('cpv_trajectory_initialize', 'error opening file "'//trim(fname)//'.pos"', 1)
      if (present(ios)) &
         ios = iostat
      if (iostat /= 0) return
      t%iounit_vel = find_free_unit()
      open (unit=t%iounit_vel, file=trim(fname)//'.vel', iostat=iostat, action='read')
      if (.not. present(ios) .and. iostat /= 0) &
         call errore('cpv_trajectory_initialize', 'error opening file "'//trim(fname)//'.vel"', 1)
      if (present(ios)) &
         ios = iostat
      if (iostat /= 0) return
      t%is_open = .true.
      t%tau_fac = tau_fac
      t%vel_fac = vel_fac
      t%tps_fac = tps_fac
      t%fname = fname
      !allocate traj
      if (present(circular)) then
         call trajectory_allocate(t%traj, natoms, 10, circular) !allocate space for 10 steps
      else
         call trajectory_allocate(t%traj, natoms, 50) !start allocating space for 50 steps
      end if
   end subroutine

   subroutine cpv_trajectory_close(t)
      implicit none
      type(cpv_trajectory), intent(inout) :: t
      logical itsopen
      if (.not. t%is_open) return
      inquire (unit=t%iounit_pos, opened=itsopen)
      if (itsopen) &
         close (t%iounit_pos)
      inquire (unit=t%iounit_vel, opened=itsopen)
      if (itsopen) &
         close (t%iounit_vel)
      t%is_open = .false.
   end subroutine

   subroutine read_header(iounit, nstep, tps, err, eof, tps_not_set)
      implicit none
      logical, intent(out) :: err, eof, tps_not_set
      character(len=256) :: line_read, tps_str
      integer, intent(in) :: iounit
      integer, intent(out) :: nstep
      real(dp), intent(out) :: tps

      !default values
      err = .false.
      eof = .false.
      tps_not_set = .false.
      tps = 0.0_dp
      nstep = 0

      !read line as string
      read (unit=iounit, fmt='(A)', err=200, end=210) line_read
      !read nstep and tps as string
      read (unit=line_read, fmt=*, err=200, end=210) nstep, tps_str
      !read tps
      read (unit=tps_str, fmt=*, err=190, end=190) tps

      return
190   continue
      write (*, *) '!WARNING! Error reading time. Header was'
      write (*, *) trim(line_read)
      tps_not_set = .true.
      return
200   continue
      write (*, *) '!ERROR! Error reading header. Header was'
      write (*, *) trim(line_read)
      err = .true.
      return
210   continue
      eof = .true.
      return
   end subroutine

   function cpv_trajectory_read_step(t) result(res)
      implicit none
      type(cpv_trajectory), intent(inout) :: t
      logical :: res, eof, err, tps_not_set
      type(timestep) :: tstep !this internally is only a pointer to a bigger allocated array
      real(dp) :: tps_
      integer :: nstep_, iatom, i
      if (.not. t%is_open) then
         res = .false.
         return
      endif

      !try to read a new step, if success return true
      !get some space (that can be used again
      call trajectory_get_temporary(t%traj, t%traj%natoms, tstep) ! first get a temporary, then eventually confirm it as valid

      !read headers
      call read_header(t%iounit_pos, tstep%nstep, tstep%tps, err, eof, tps_not_set)
      if (err .or. eof) goto 100
      call read_header(t%iounit_vel, nstep_, tps_, err, eof, tps_not_set)
      if (err .or. eof) goto 100
      if (tstep%nstep /= nstep_ .or. tstep%tps /= tps_) then
         write (*, *) 'error in reading: inconsistent timestep headers in .pos and .vel files!', &
            nstep_, tps_, tstep%nstep, tstep%tps
         goto 100
      end if

      ! read trajectory
      do iatom = 1, t%traj%natoms
         read (t%iounit_pos, *, err=100, end=100) (tstep%tau(i, iatom), i=1, 3)
         read (t%iounit_vel, *, err=100, end=100) (tstep%vel(i, iatom), i=1, 3)
      enddo

      !apply factors
      tstep%tau = tstep%tau*t%tau_fac
      tstep%vel = tstep%vel*t%vel_fac
      tstep%tps = tstep%tps*t%tps_fac

      !store temporary in trajectory
      call trajectory_push_back_last_temporary(t%traj)

      res = .true.

      return

      !reading error
100   res = .false.
      call cpv_trajectory_close(t)

   end function

   subroutine cpv_trajectory_get_last_step(t, tstep)
      implicit none
      type(cpv_trajectory), intent(in) :: t
      type(timestep), intent(out) :: tstep

      ! get the timestep data
      call trajectory_get(t%traj, t%traj%nsteps, tstep)

   end subroutine
   subroutine cpv_trajectory_get_step(t, idx, tstep)
      implicit none
      type(cpv_trajectory), intent(in) :: t
      type(timestep), intent(out) :: tstep
      integer, intent(in) :: idx

      ! get the timestep data
      call trajectory_get(t%traj, idx, tstep)

   end subroutine

   subroutine cpv_trajectory_deallocate(t)
      implicit none
      type(cpv_trajectory), intent(inout) :: t

      !deallocate stuff
      call trajectory_deallocate(t%traj)
      call cpv_trajectory_close(t)

   end subroutine

end module
