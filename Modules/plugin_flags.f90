!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE plugin_flags
  !=--------------------------------------------------------------------------=!
  !! This module contains all basic variables that controls
  !! the use or not use of plugins.
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
  !
  ! ...   declare execution control variables
  !
  CHARACTER(LEN=256), PUBLIC :: plugin_name
  LOGICAL, PUBLIC :: use_plumed
  LOGICAL, PUBLIC :: use_pw2casino
  LOGICAL, PUBLIC :: use_environ
  LOGICAL, PUBLIC :: use_partn
  LOGICAL, PUBLIC :: use_oscdft
  !
END MODULE plugin_flags
