!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------=
   MODULE stick_base
!=----------------------------------------------------------------------=

      USE kinds
      USE io_global, ONLY: ionode

        IMPLICIT NONE
        PRIVATE
        SAVE

        PUBLIC :: sticks_maps


!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

      SUBROUTINE  sticks_maps( tk, ub, lb, b1, b2, b3, gcut, gcutw, gcuts, st, stw, sts )

          USE mp, ONLY: mp_sum
          USE mp_global, ONLY: mpime, nproc, group

          LOGICAL, INTENT(IN) :: tk    !  if true use the full space grid
          INTEGER, INTENT(IN) :: ub(:) !  upper bounds for i-th grid dimension
          INTEGER, INTENT(IN) :: lb(:) !  lower bounds for i-th grid dimension
          REAL(dbl) , INTENT(IN) :: b1(:), b2(:), b3(:) ! reciprocal space base vectors
          REAL(dbl) , INTENT(IN) :: gcut   ! cut-off for potentials
          REAL(dbl) , INTENT(IN) :: gcutw  ! cut-off for plane waves
          REAL(dbl) , INTENT(IN) :: gcuts  ! cut-off for smooth mesh
          INTEGER, INTENT(OUT) :: st( lb(1): ub(1), lb(2):ub(2) ) ! stick map for potential
          INTEGER, INTENT(OUT) :: stw(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions
          INTEGER, INTENT(OUT) :: sts(lb(1): ub(1), lb(2):ub(2) ) ! stick map for smooth mesh

          INTEGER :: i, j, k, kip
          REAL(dbl) :: gsq

          stw  = 0
          st   = 0
          sts  = 0

! ...     Here find the basic maps of sticks st and stw for the potential
! ...     cut-off gcut and the wavefunction cut-off gcutw
! ...     st(i,j) will contain the number of G vectors of the stick whose
! ...     indices are (i,j). 
 
          IF( .NOT. tk ) THEN

            kip = 0 + ABS(lb(3)) + 1
            IF( MOD( kip, nproc ) == mpime ) THEN
              st (0,0) = st (0,0) + 1
              stw(0,0) = stw(0,0) + 1
            END IF

            DO i= 0, 0 
              DO j= 0, 0 
                DO k= 1, ub(3)
                  kip = k + ABS(lb(3)) + 1
                  IF( MOD( kip, nproc ) == mpime ) THEN
                    gsq=    (REAL(i)*b1(1)+REAL(j)*b2(1)+REAL(k)*b3(1) )**2
                    gsq=gsq+(REAL(i)*b1(2)+REAL(j)*b2(2)+REAL(k)*b3(2) )**2
                    gsq=gsq+(REAL(i)*b1(3)+REAL(j)*b2(3)+REAL(k)*b3(3) )**2
                    IF(gsq.LE.gcut ) THEN
                      st(i,j) = st(i,j) + 1
                      IF(gsq.LE.gcutw) THEN
                        stw(i,j) = stw(i,j) + 1
                      END IF
                      IF(gsq.LE.gcuts) THEN
                        sts(i,j) = sts(i,j) + 1
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END DO

            DO i = 0, 0
              DO j = 1, ub(2)
                DO k = lb(3), ub(3)
                  kip = k + ABS(lb(3)) + 1
                  IF( MOD( kip, nproc) == mpime ) THEN
                    gsq=    (REAL(i)*b1(1)+REAL(j)*b2(1)+REAL(k)*b3(1) )**2
                    gsq=gsq+(REAL(i)*b1(2)+REAL(j)*b2(2)+REAL(k)*b3(2) )**2
                    gsq=gsq+(REAL(i)*b1(3)+REAL(j)*b2(3)+REAL(k)*b3(3) )**2
                    IF(gsq.LE.gcut ) THEN
                      st(i,j) = st(i,j) + 1
                      IF(gsq.LE.gcutw) THEN
                        stw(i,j) = stw(i,j) + 1
                      END IF
                      IF(gsq.LE.gcuts) THEN
                        sts(i,j) = sts(i,j) + 1
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END DO

            DO i = 1, ub(1)
              DO j = lb(2), ub(2)
                DO k = lb(3), ub(3)
                  kip = k + ABS(lb(3)) + 1
                  IF( MOD( kip, nproc) == mpime ) THEN
                    gsq=    (REAL(i)*b1(1)+REAL(j)*b2(1)+REAL(k)*b3(1) )**2
                    gsq=gsq+(REAL(i)*b1(2)+REAL(j)*b2(2)+REAL(k)*b3(2) )**2
                    gsq=gsq+(REAL(i)*b1(3)+REAL(j)*b2(3)+REAL(k)*b3(3) )**2
                    IF(gsq.LE.gcut ) THEN
                      st(i,j) = st(i,j) + 1
                      IF(gsq.LE.gcutw) THEN
                        stw(i,j) = stw(i,j) + 1
                      END IF
                      IF(gsq.LE.gcuts) THEN
                        sts(i,j) = sts(i,j) + 1
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END DO

          ELSE

            DO i= lb(1), ub(1)
              DO j= lb(2), ub(2)
                DO k= lb(3), ub(3)
                  kip = k + ABS(lb(3)) + 1
                  IF( MOD( kip, nproc ) == mpime ) THEN
                    gsq=    (REAL(i)*b1(1)+REAL(j)*b2(1)+REAL(k)*b3(1) )**2
                    gsq=gsq+(REAL(i)*b1(2)+REAL(j)*b2(2)+REAL(k)*b3(2) )**2
                    gsq=gsq+(REAL(i)*b1(3)+REAL(j)*b2(3)+REAL(k)*b3(3) )**2
                    IF(gsq.LE.gcut ) THEN
                      st(i,j) = st(i,j) + 1
                      IF(gsq.LE.gcutw) THEN
                        stw(i,j) = stw(i,j) + 1
                      END IF
                      IF(gsq.LE.gcuts) THEN
                        sts(i,j) = sts(i,j) + 1
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END DO

          END IF

          CALL mp_sum(st  ,group)
          CALL mp_sum(stw ,group)
          CALL mp_sum(sts ,group)

! Test sticks
!          write(6,*) 'testtesttesttesttesttesttesttesttesttest'
!          DO i = lb(1), ub(1)
!            DO j = lb(2), ub(2)
!              write(6,'(2I4,I6)') i,j,stw(i,j)
!            END DO
!          END DO
!          write(6,*) 'testtesttesttesttesttesttesttesttesttest'
! Test sticks

        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------=
   END MODULE stick_base
!=----------------------------------------------------------------------=
