!
! Copyright (C) 2003-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE global_version
  !! QE version
  IMPLICIT NONE
  !
  SAVE
  !
#include "qe_version.h"
  !
END MODULE global_version
