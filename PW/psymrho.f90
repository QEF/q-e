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
  USE mp_global, ONLY : me_pool
  USE fft_base,  ONLY : dfftp, grid_gather, grid_scatter
  !
  IMPLICIT NONE
  !
  INTEGER        :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau
  REAL (DP) :: rho( dfftp%nnr )
  !
#if defined  (__PARA)
  !
  REAL (DP), ALLOCATABLE :: rrho(:)
  !
  !
  ALLOCATE (rrho( nrx1 * nrx2 * nrx3))    
  !
  CALL grid_gather( rho, rrho )
  !
  IF ( me_pool == 0 ) &
     CALL symrho( rrho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
  !
  CALL grid_scatter( rrho, rho )
  !
  DEALLOCATE( rrho )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE psymrho

