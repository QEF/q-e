!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!



subroutine show_memory ()
#include "machine.h"
  implicit none
  !      WRITE( stdout,'(5x,"Current number of allocated pointers:",i8)') nptr

  !WRITE( stdout, '(5x,"Dynamical memory: ",f6.2,"Mb current, ", &
  !     &            f6.2,"Mb maximum")') real (totsize)  / 1000000, &
  !     &     real (maxsize)  / 1000000
  return
end subroutine show_memory

