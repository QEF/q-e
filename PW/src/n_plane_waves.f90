!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
INTEGER FUNCTION n_plane_waves (gcutw, nks, xk, g, ngm) RESULT(npwx)
  !-----------------------------------------------------------------------
  !
  ! Find maximum number of plane waves over all k-points
  !
  USE kinds, ONLY: DP
  USE mp,       ONLY : mp_max
  USE mp_pools, ONLY : inter_pool_comm
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nks, ngm
  real(DP),INTENT(in) :: gcutw, xk (3, nks), g (3, ngm)
  !
  INTEGER :: nk, ng, npw
  real(DP) :: q2
  !
  npwx = 0
  DO nk = 1, nks
     npw = 0
     DO ng = 1, ngm
        q2 = (xk (1, nk) + g (1, ng) ) **2 + (xk (2, nk) + g (2, ng) ) ** &
             2 + (xk (3, nk) + g (3, ng) ) **2
        IF (q2 <= gcutw) THEN
           !
           ! here if |k+G|^2 <= Ecut increase the number of G inside the sphere
           !
           npw = npw + 1
        ELSE
           IF ( sqrt (g (1, ng) **2 + g (2, ng) **2 + g (3, ng) **2) > &
                sqrt (xk(1, nk) **2 + xk(2, nk) **2 + xk(3, nk) **2) + &
                sqrt (gcutw) ) GOTO 100
           !
           ! if |G| > |k| + sqrt(Ecut)  stop search
           !
        ENDIF
     ENDDO
100  npwx = max (npwx, npw )
  ENDDO
  IF (npwx <= 0) CALL errore ('n_plane_waves', &
                'No plane waves found: running on too many processors?', 1)
  !
  ! when using pools, set npwx to the maximum value across pools
  ! (you may run into trouble at restart otherwise)
  !
  CALL mp_max ( npwx, inter_pool_comm )
  !
END FUNCTION n_plane_waves
