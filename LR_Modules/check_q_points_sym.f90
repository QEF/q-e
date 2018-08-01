!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION check_q_points_sym(nqs, q, at, bg, nsym, s, invs, &
                                    nq1, nq2, nq3)
!
!  This function returns .true. if the mesh of q points given as input
!  is compatible with the FFT mesh. It returns .false. if a rotation of
!  the point group gives a q point that is not in the FFT mesh.
!
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nqs, nsym
INTEGER, INTENT(IN) :: nq1, nq2, nq3
INTEGER, INTENT(IN) :: s(3,3,48), invs(48)
REAL(DP), INTENT(IN) :: q(3,nqs), at(3,3), bg(3,3)

INTEGER :: nq, ipol, icar, iq, jq
INTEGER :: nr(3), isq (48), imq
LOGICAL :: lq
REAL(DP) :: xq, sxq(3,48)
REAL(DP) :: eps=1.d-5

nr(1)=nq1
nr(2)=nq2
nr(3)=nq3
lq = .TRUE.
DO iq = 1,nqs
   call star_q (q(:,iq), at, bg, nsym, s, invs, nq, sxq, isq, imq, .FALSE. )
   DO jq=1,nq
      DO ipol=1,3
         xq = 0.0d0
         DO icar=1,3
            xq = xq + at(icar,ipol) * sxq(icar,jq) * nr(ipol)
         END DO
         lq = lq .AND. (ABS(NINT(xq) - xq) .LT. eps)
      ENDDO
   ENDDO
ENDDO
check_q_points_sym=lq

RETURN
END FUNCTION check_q_points_sym
