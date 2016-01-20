!
! Copyright (C) 2009 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ws_base
  !============================================================================
  !
  !   Module containing type definitions and auxiliary routines to deal with
  !   basic operations on the Wigner-Seitz cell associated to a given set
  !   of Bravais fundamental lattice vectors.
  !
  !   Should contain low level routines and no reference to other modules
  !   (with the possible exception of kinds and parameters) so as to be 
  !   call-able from any other module.
  !
  ! content:
  !
  ! - ws_type   : derived type definition used to encoded the auxiliary 
  !               quantities needed by the other WS functions or routines
  !
  ! - ws_init(a,ws) 
  !             : a routine that initializes a ws_type variable
  !
  ! - ws_clear(ws) 
  !             : a routine that un-sets a ws_type variable
  !
  ! - ws_test(ws)
  !             : a routine that tests whether a ws_type variable has been
  !               initialized
  !
  ! - ws_vect(r,ws,r_ws) 
  !             : a routine that given a vector returns an equivalent 
  !               vector inside the WS cell
  !
  ! - ws_dist(r,ws)
  !             : a routine that, given a vector, returns the shortest 
  !               distance from any point in the Bravais lattice
  !
  ! - ws_weight(r,ws)
  !             : a routine that given a vector 
  !               returns 1.0      if the vector is inside the WS cell
  !               returns 0.0      if the vector is outside the WS cell
  !               returns 1/(1+NR) if the vector is on the frontier of the 
  !                                WS cell and NR is the number of Bravais 
  !                                lattice points whose distance is the same
  !                                as the one from the origin
  !
  !============================================================================
  !
  USE kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  TYPE ws_type

    PRIVATE ! this means (I hope) that internal variables can only
            ! be accessed through calls of routines inside the module.
    REAL(DP) ::  &
        a(3,3),  & ! the fundamental Bravais lattice vectors
        aa(3,3), & ! a^T*a
        b(3,3),  & ! the inverse of a, i.e. the transponse of the fundamental 
                   ! reciprocal lattice vectors
        norm_b(3)  ! the norm of the fundamental reciprocal lattice vectors
    LOGICAL  ::  &
        initialized = .FALSE. ! .TRUE. when  initialized

  END TYPE ws_type

  PRIVATE
  PUBLIC :: ws_type, ws_init, ws_clean, ws_test, ws_vect, ws_dist, ws_weight, ws_dist_stupid

  !============================================================================
  !
 CONTAINS
!---------------------------------------------------------------
    SUBROUTINE ws_init(a,ws)
!---------------------------------------------------------------
      USE matrix_inversion
      REAL(DP), INTENT(IN) :: a(3,3)
      TYPE(ws_type), INTENT(OUT) :: ws
      INTEGER :: i
      !
      ws%a = a
      CALL invmat( 3, ws%a, ws%b )
      ws%aa = MATMUL(TRANSPOSE(a),a)
      do i=1,3
         ws%norm_b(i) =  DSQRT(SUM(ws%b(i,:)*ws%b(i,:)))
      end do
      ws%initialized = .TRUE.

      RETURN
    END SUBROUTINE ws_init
!
!---------------------------------------------------------------
    SUBROUTINE ws_clean(ws)
!---------------------------------------------------------------
      TYPE(ws_type), INTENT(OUT) :: ws

      ws%initialized = .FALSE.

      RETURN
    END SUBROUTINE ws_clean
!
!---------------------------------------------------------------
    SUBROUTINE ws_test(ws)
!---------------------------------------------------------------
      TYPE(ws_type), INTENT(IN) :: ws

      IF (.NOT.ws%initialized) CALL errore &
               ('ws_test','trying to use an uninitialized ws_type variable',1)

      RETURN
    END SUBROUTINE ws_test

!---------------------------------------------------------------
    SUBROUTINE ws_vect(r,ws,r_ws)
!---------------------------------------------------------------
      REAL(DP), INTENT(IN) :: r(3)
      TYPE(ws_type), INTENT(IN) :: ws
      REAL(DP), INTENT(OUT) :: r_ws(3)
      REAL(DP) :: x(3), y(3), c, ctest
      INTEGER :: lb(3), ub(3), i1, i2, i3, m(3)

      CALL ws_test(ws)

      x = MATMUL(ws%b,r)
      x(:) = x(:) - NINT(x(:))
      c  = SUM(x*MATMUL(ws%aa,x))
      m = 0

      lb(:) =  NINT ( x(:) - DSQRT (c) * ws%norm_b(:) )
            ! CEILING should be enough for lb but NINT might be safer
      ub(:) =  NINT ( x(:) + DSQRT (c) * ws%norm_b(:) )
            ! FLOOR should be enough for ub but NINT might be safer

      DO i1 = lb(1), ub(1)
         DO i2 = lb(2), ub(2)
            DO i3 = lb(3), ub(3)
               y = x - (/i1,i2,i3/)
               ctest = SUM(y*MATMUL(ws%aa,y))
               IF (ctest < c) THEN
                  c = ctest
                  m = (/i1,i2,i3/)
               END IF
            END DO
         END DO
      END DO

      y = x-m
      r_ws =  MATMUL(ws%a,y)

      RETURN
    END SUBROUTINE ws_vect
!
!---------------------------------------------------------------
    FUNCTION ws_dist_stupid(r,ws)
!---------------------------------------------------------------
      REAL(DP), INTENT(IN) :: r(3)
      TYPE(ws_type), INTENT(IN) :: ws
      REAL(DP) :: ws_dist_stupid
      REAL(DP) :: r_ws(3)
 
      integer :: i1,i2,i3
      real(DP) :: rr, rmin, rtest(3)

      CALL ws_test(ws)

      rmin = 1.d+9
      
      do i1=-3,3
      do i2=-3,3
      do i3=-3,3
         rtest(:) = r(:) + ws%a(:,1)*i1 + ws%a(:,2)*i2 + ws%a(:,3)*i3
         rr = sum(rtest(:)**2)
         if (rr < rmin) rmin = rr
      end do
      end do
      end do

      ws_dist_stupid = DSQRT(rmin)

      RETURN
    END FUNCTION ws_dist_stupid
!
!---------------------------------------------------------------
    FUNCTION ws_dist(r,ws)
!---------------------------------------------------------------
      REAL(DP), INTENT(IN) :: r(3)
      TYPE(ws_type), INTENT(IN) :: ws
      REAL(DP) :: ws_dist
      REAL(DP) :: r_ws(3)

      CALL ws_test(ws)

      CALL ws_vect(r,ws,r_ws)

      ws_dist = DSQRT(SUM(r_ws**2))

      RETURN
    END FUNCTION ws_dist
!
!---------------------------------------------------------------
    FUNCTION ws_weight(r,ws)
!---------------------------------------------------------------
      REAL(DP), INTENT(IN) :: r(3)
      TYPE(ws_type), INTENT(IN) :: ws
      REAL(DP) :: ws_weight

      REAL(DP) :: x(3), y(3), c, ctest
      INTEGER :: lb(3), ub(3), i1, i2, i3, m(3)

      REAL(DP), PARAMETER  :: eps6  = 1.0E-6_DP

      ws_weight = 0.0_DP

      CALL ws_test(ws)

      x = MATMUL(ws%b,r)
      c  = SUM(x*MATMUL(ws%aa,x))

      lb(:) =  NINT ( x(:) - DSQRT (c) * ws%norm_b(:) )
            ! CEILING should be enough for lb but NINT might be safer
      ub(:) =  NINT ( x(:) + DSQRT (c) * ws%norm_b(:) )
            ! FLOOR should be enough for ub but NINT might be safer

      DO i1 = lb(1), ub(1)
         DO i2 = lb(2), ub(2)
            DO i3 = lb(3), ub(3)
               y = x - (/i1,i2,i3/)
               ctest = SUM(y*MATMUL(ws%aa,y))
               IF (ctest < c - eps6 ) THEN
                  ws_weight = 0.0_DP
                  RETURN
               END IF
               IF (ctest < c + eps6 ) THEN
                  ws_weight = ws_weight + 1.0_DP
               END IF
            END DO
         END DO
      END DO

      IF (ws_weight == 0.0_DP) CALL errore ('ws_weight','unexpected error',1)

      ws_weight = 1.0_dp / ws_weight

      RETURN
    END FUNCTION ws_weight
!
END MODULE ws_base
