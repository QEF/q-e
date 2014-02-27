!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE plugin_variables
  !=--------------------------------------------------------------------------=!
  !
  ! ... this module contains all basic variables possibly
  ! ... used by plug ins
  !----------------------------------------------
  !
  USE kinds
  USE parameters
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: plugin_etot
  !
  REAL(DP) :: plugin_etot
  !
END MODULE plugin_variables
