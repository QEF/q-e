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
  use pwcom
  implicit none
  !
  call start_clock ('ions')
  conv_ions = .true.
  !
  ! recover from a previous run, if appropriate
  !

  if (restart.and.iswitch.ge.0) call restart_in_ions

  if (lforce) call forces

  if (lstres) call stress

  if (iswitch.gt.0) then
     call move_ions
#ifdef PARA
      call check (3 * nat, tau)
#endif
      !
      ! save restart information
      !
      call write_config_to_file
      call save_in_ions
  end if
  call stop_clock ('ions')
  return
end subroutine ions

