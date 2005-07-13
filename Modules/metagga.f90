!
! Copyright (C) 2002-2004 quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
MODULE metaflag
!------------------------------------------------------------------------------!
  !
  IMPLICIT NONE
  SAVE
  LOGICAL :: ismeta = .FALSE.   ! ismeta = .TRUE. if metagga is used
			        ! default value is .FALSE.
  !
!------------------------------------------------------------------------------!
END MODULE metaflag
!------------------------------------------------------------------------------!
