!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt( npwx, npw, nbnd, npol, psi, hpsi, spsi, e, &
                         uspp, eigen, reorder, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization.
  !
  USE kinds,         ONLY : DP
  USE control_flags, ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi (npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: spsi(npwx*npol,nbnd)
  REAL(DP),    INTENT(OUT)   :: e(nbnd)
  LOGICAL,     INTENT(IN)    :: uspp
  LOGICAL,     INTENT(IN)    :: eigen
  LOGICAL,     INTENT(IN)    :: reorder
  INTEGER,     INTENT(IN)    :: nbsize
  !
  CALL start_clock( 'gsorth' )
  !
  IF ( gamma_only ) THEN
     !
     CALL gram_schmidt_gamma( npwx, npw, nbnd, psi, hpsi, spsi, e, &
                              uspp, eigen, reorder, nbsize )
     !
  ELSE
     !
     CALL gram_schmidt_k( npwx, npw, nbnd, npol, psi, hpsi, spsi, e, &
                          uspp, eigen, reorder, nbsize )
     !
  END IF
  !
  CALL stop_clock( 'gsorth' )
  !
END SUBROUTINE gram_schmidt
