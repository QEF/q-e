
! Copyright (C) 2010-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE set_engine_output()
  !-----------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode
  !
  ! ... assign a unit number different from 6 to stdout,
  ! ... or else the output NEB routines will be redirected as well!
  ! ... only processors performing I/O need to be redirected
  !
  if(ionode) stdout = 66
  !
END SUBROUTINE set_engine_output
!
