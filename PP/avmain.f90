!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
program macroscopic

  character :: nodenumber * 3
  call start_postproc (nodenumber)
  !
  !      work in parallel machines but only with one node!!!
  !
  call average
end program macroscopic
