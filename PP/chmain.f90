!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
program chdens

  character :: nodenumber * 3
  call start_postproc (nodenumber)
  !
  ! Works on parallel machines but only with one node !!!
  !
  call do_chdens
  call stop_pp
  stop
end program chdens
