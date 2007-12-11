!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE ccalbec( nkb, npwx, npw, nbnd, bec, vkb, psi )
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine computes the dot product of the beta functions
  ! ... and the wavefunctions, and save them in the array bec.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... here the dummy variables
  !
  INTEGER :: nkb, npwx, npw, nbnd
    ! input: the total number of beta functions
    ! input: the maximum number of plane waves
    ! input: the length of the vectors
    ! input: the number of bands
  COMPLEX(DP) ::  vkb(npwx,nkb), psi(npwx,nbnd), bec(nkb,nbnd)
    ! input: the FT of the beta functions
    ! input: the wavefunctions
    ! output: dot product of the beta and the wavefunctions
  !
  !
  IF ( nkb == 0 ) RETURN
  !
  CALL start_clock( 'calbec' )
  !
  IF ( nbnd == 1 ) THEN
     !
     CALL ZGEMV( 'C', npw, nkb, ONE, vkb, npwx, psi, 1, ZERO, bec, 1 )
     !
  ELSE
     !
     CALL ZGEMM( 'C', 'N', nkb, nbnd, npw, ONE, &
                 vkb, npwx, psi, npwx, ZERO, bec, nkb )
     !
  END IF
  !
  CALL reduce( 2 * nkb * nbnd, bec )
  !
  CALL stop_clock( 'calbec' )
  !
  RETURN
  !
END SUBROUTINE ccalbec
