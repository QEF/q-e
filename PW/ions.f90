!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ions
  !-----------------------------------------------------------------------
  !
  USE control_flags, ONLY: conv_ions, restart, lscf, lmd, lbfgs, loldbfgs
  USE force_mod, ONLY: lforce, lstres
  !
  implicit none
  !
  call start_clock ('ions')
  conv_ions = .true.
  !
  ! recover from a previous run, if appropriate
  !
  if (restart.and.lscf) call restart_in_ions
  !
  if (lforce) call forces
  !
  if (lstres) call stress
  !
  if (lmd .OR. lbfgs .OR. loldbfgs) then
     call move_ions
      !
      ! save restart information
      !
      call write_config_to_file
      call save_in_ions
  end if
  call stop_clock ('ions')
  return
end subroutine ions

