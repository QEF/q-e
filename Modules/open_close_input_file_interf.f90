!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE open_close_input_file_interf
  !-----------------------------------------------------------------------------
  !
  ! ...  this module contains interface routines for open and close
  ! ... input file xml or not
  ! ...  ---------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: open_input_file
  PUBLIC :: close_input_file
  !
  !
  INTERFACE open_input_file
  SUBROUTINE open_input_file_x(lxmlinput,attr,unit)
  !
  ! ...  this subroutine opens the input file standard input ( unit 5 )
  ! ...  Use "-input filename" to read input from file "filename":
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  or xml input
  ! ...  ---------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  !
  USE io_global,     ONLY : stdout, xmlinputunit
  !
  USE iotk_module,           ONLY : iotk_open_read, iotk_close_read, iotk_attlenx
  !
  IMPLICIT NONE
  !
  LOGICAL, intent(inout), optional :: lxmlinput
  CHARACTER (len=iotk_attlenx), intent(inout), optional :: attr
  INTEGER, intent(in), optional :: unit
  !
END SUBROUTINE open_input_file_x
END INTERFACE
  !
  !
INTERFACE close_input_file
SUBROUTINE close_input_file_x(lxmlinput,unit)
  !
  ! ...  this subroutine close the input file for the specified unit
  ! ...  ( default is unit 5 )
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  or xml input
  ! ...  ---------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  !
  USE iotk_module,           ONLY : iotk_close_read
  !
  IMPLICIT NONE
  !
  LOGICAL, intent(inout), optional :: lxmlinput
  INTEGER, intent(in), optional :: unit
  !
END SUBROUTINE close_input_file_x
  !
END INTERFACE
  !
END MODULE open_close_input_file_interf
