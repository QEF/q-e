!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE ring_io_units_module
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the I/O units and files used by  "path"-routines
  ! 
  ! ... Written by Carlo Sbraccia ( 2003-2006 )
  IMPLICIT NONE
  !
  !
  SAVE
  !
  INTEGER :: iunpath = 6
  
  INTEGER :: iunnewimage = 28 ! unit for parallelization among images
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  CONTAINS
  !
  SUBROUTINE set_io_units()
  !
  IMPLICIT NONE

  iunnewimage = find_free_unit()
  !
  END SUBROUTINE set_io_units
  ! 
END MODULE ring_io_units_module
