
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE set_engine_io_units()
  !-----------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode
  INTEGER, EXTERNAL :: find_free_unit
  !
  if(ionode) stdout = find_free_unit()
  !
END SUBROUTINE set_engine_io_units
!
