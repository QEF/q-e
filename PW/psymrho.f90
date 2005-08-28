!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE psymrho( rho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
  !----------------------------------------------------------------------------
  !
  ! ...  p-symmetrize the charge density.
  !
  USE kinds,     ONLY : DP
  USE pfft,      ONLY : nxx
  USE mp_global, ONLY : me_pool
  !
  IMPLICIT NONE
  !
  INTEGER        :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau
  REAL (DP) :: rho(nxx)
  !
#if defined  (__PARA)
  !
  REAL (DP), ALLOCATABLE :: rrho(:)
  !
  !
  ALLOCATE (rrho( nrx1 * nrx2 * nrx3))    
  !
  CALL gather( rho, rrho )
  !
  IF ( me_pool == 0 ) &
     CALL symrho( rrho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
  !
  CALL scatter( rrho, rho )
  !
  DEALLOCATE( rrho )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE psymrho

