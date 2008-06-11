!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
!--------------------------------------------------------------------
      SUBROUTINE multiscale( f, n1, n2, n3, d1n, d2n, d3n,      &
            g, m1, m2, m3, d1m, d2m, d3m, orig1, orig2, orig3,  &
            do_per )
!--------------------------------------------------------------------
      !
      ! ...
      !
      USE kinds,         ONLY : DP
      !
      IMPLICIT NONE
      !
      INTEGER                :: n1, n2, n3
      INTEGER                :: m1, m2, m3
      REAL( DP )             :: d1n, d2n, d3n
      REAL( DP )             :: d1m, d2m, d3m
      REAL( DP )             :: f( n1 * n2 * n3 )
      REAL( DP )             :: g( m1 * m2 * m3 )
      LOGICAL                :: do_per
      !
      REAL( DP )             :: t1, t2, t3
      REAL( DP )             :: x1, x2, x3
      REAL( DP )             :: ff
      REAL( DP )             :: dff1, dff2, dff3
      REAL( DP )             :: orig1, orig2, orig3
      REAL( DP )             :: l1, l2, l3
      !
      INTEGER                :: i1, i2, i3, i
      INTEGER                :: j1, j2, j3, j
      INTEGER                :: a, b, c
      INTEGER                :: bound1, bound2, bound3
      !
      INTEGER, EXTERNAL      :: compindex
      REAL( DP ), EXTERNAL   :: pinterp
      REAL( DP ), EXTERNAL   :: qinterp
      INTEGER, EXTERNAL      :: bound
      !
      REAL( DP ), PARAMETER  :: epsx = 0.D-1
      !
      l1 = d1n * DBLE( n1 )
      l2 = d2n * DBLE( n2 )
      l3 = d3n * DBLE( n3 )
      !
      g = 0.D0
      bound1 = 0
      bound2 = 0
      bound3 = 0
      !
      DO j1 = 1, m1
       DO j2 = 1, m2
        DO j3 = 1, m3
         !
         x1 = DBLE( j1 - 1 ) * d1m + orig1
         x2 = DBLE( j2 - 1 ) * d2m + orig2
         x3 = DBLE( j3 - 1 ) * d3m + orig3
         !
         IF(    x1 .GE. - epsx .AND. x1 .LE. l1 + epsx                 &
          .AND. x2 .GE. - epsx .AND. x2 .LE. l2 + epsx                 &
          .AND. x3 .GE. - epsx .AND. x3 .LE. l3 + epsx ) THEN
         !
         j = compindex( j1, j2, j3, m1, m2, m3 )
         !
         t1 = x1 / d1n
         t2 = x2 / d2n 
         t3 = x3 / d3n
         !
         i1 = MIN( INT( t1 ) + 1, n1 - 1 )
         i2 = MIN( INT( t2 ) + 1, n2 - 1 )
         i3 = MIN( INT( t3 ) + 1, n3 - 1 )
         !
         t1 = t1 - DBLE( i1 - 1 )
         t2 = t2 - DBLE( i2 - 1 )
         t3 = t3 - DBLE( i3 - 1 )
         !
         IF( .NOT. do_per ) THEN
           bound1 = bound( i1, n1 )
           bound2 = bound( i2, n2 )
           bound3 = bound( i3, n3 )
         END IF
         !
         DO a = 0, 1
          DO b = 0, 1
           DO c = 0, 1
            !
            ff = f( compindex( i1 + a, i2 + b, i3 + c, n1, n2, n3 ) )
            g( j ) = g( j ) + ff * pinterp( t1, a, bound1 )            &
                                 * pinterp( t2, b, bound2 )            &
                                 * pinterp( t3, c, bound3 )
            dff1 = 0.5D0 * (                                           &
              f( compindex( i1 + a + 1, i2 + b, i3 + c, n1, n2, n3 ) ) &
            - f( compindex( i1 + a - 1, i2 + b, i3 + c, n1, n2, n3 ) ) )
            dff2 = 0.5D0 * (                                           &
              f( compindex( i1 + a, i2 + b + 1, i3 + c, n1, n2, n3 ) ) &
            - f( compindex( i1 + a, i2 + b - 1, i3 + c, n1, n2, n3 ) ) )
            dff3 = 0.5D0 * (                                           &
              f( compindex( i1 + a, i2 + b, i3 + c + 1, n1, n2, n3 ) ) &
            - f( compindex( i1 + a, i2 + b, i3 + c - 1, n1, n2, n3 ) ) )
            g( j ) = g( j ) + dff1 * qinterp( t1, a, bound1 )          &
                                   * pinterp( t2, b, bound2 )          &
                                   * pinterp( t3, c, bound3 )          &
                            + dff2 * pinterp( t1, a, bound1 )          &
                                   * qinterp( t2, b, bound2 )          &
                                   * pinterp( t3, c, bound3 )          &
                            + dff3 * pinterp( t1, a, bound1 )          &
                                   * pinterp( t2, b, bound2 )          &
                                   * qinterp( t3, c, bound3 )
            ! 
           END DO
          END DO
         END DO
         !
         END IF
         !
        END DO
       END DO
      END DO
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE multiscale
!--------------------------------------------------------------------

!--------------------------------------------------------------------
      FUNCTION pinterp( x, side, bound )
!--------------------------------------------------------------------
      !
      ! ... Calculates the P(x) where Pinterp
      ! ... is P the interpolation polynomial satifying
      ! ... P'(0)=0,P'(1)=0,P(0)=1,P(1)=0 for side=0,bound=0
      ! ... P'(0)=0,P'(1)=0,P(0)=0,P(1)=1 for side=1,bound=0
      ! ... P'(1)=0,P(0)=0,P(1)=1 for side=1,bound=-1
      ! ... P'(1)=0,P(0)=1,P(1)=0 for side=0,bound=-1
      ! ... P'(0)=0,P(0)=0,P(1)=1 for side=1,bound=1
      ! ... P'(0)=0,P(0)=1,P(1)=0 for side=0,bound=1
      !
      USE kinds,         ONLY : DP
      !
      IMPLICIT NONE
      !
      REAL( DP ) :: pinterp
      REAL( DP ) :: x
      INTEGER    :: side
      INTEGER    :: bound
      !
      IF( bound == 0 .AND. side == 1 ) THEN
        pinterp = 3.D0 * x * x - 2.D0 * x * x * x
      ELSE IF( bound == 0 .AND. side == 0 ) THEN
        pinterp = 1.D0 - 3.D0 * x * x + 2.D0 * x * x * x
      ELSE IF( bound == - 1 .AND. side == 0 ) THEN
        pinterp = 1.D0 - 2.D0 * x + x * x
      ELSE IF( bound == - 1 .AND. side ==  1 ) THEN
        pinterp = 2.D0 * x - x * x
      ELSE IF( bound == 1 .AND. side == 1 ) THEN
        pinterp = x * x
      ELSE IF( bound == 1 .AND. side == 0 ) THEN
        pinterp =  1 - x * x
      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION pinterp
!--------------------------------------------------------------------

!--------------------------------------------------------------------
      FUNCTION dpinterp( x, side, bound )
!--------------------------------------------------------------------
      !
      ! ... Calculates the **DERIVATIVE** OF P(x) where Pinterp
      !
      USE kinds,         ONLY : DP
      !
      IMPLICIT NONE
      !
      REAL( DP ) :: dpinterp
      REAL( DP ) :: x
      INTEGER    :: side
      INTEGER    :: bound
      !
      IF( bound == 0 .AND. side == 1 ) THEN
        dpinterp = 6.D0 * x  - 6.D0 * x * x
      ELSE IF( bound == 0 .AND. side == 0 ) THEN
        dpinterp = - 6.D0 * x + 6.D0 * x * x
      ELSE IF( bound == - 1 .AND. side == 0 ) THEN
        dpinterp = - 2.D0  + 2 * x
      ELSE IF( bound == - 1 .AND. side ==  1 ) THEN
        dpinterp = 2.D0 - 2 * x
      ELSE IF( bound == 1 .AND. side == 1 ) THEN
        dpinterp = 2 * x
      ELSE IF( bound == 1 .AND. side == 0 ) THEN
        dpinterp =  - 2 * x
      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION dpinterp
!--------------------------------------------------------------------


!--------------------------------------------------------------------
      FUNCTION qinterp( x, side, bound )
!--------------------------------------------------------------------
      !
      ! ... Calculates the Q(x) where Qinterp
      ! ... is Q the interpolation polynomial satifying
      ! ... Q'(0)=1,Q'(1)=0,Q(0)=0,Q(1)=0 for side=0,bound=0
      ! ... Q'(0)=0,Q'(1)=1,Q(0)=0,Q(1)=0 for side=1,bound=0
      ! ... Q'(1)=1,Q(0)=0,Q(1)=0 for side=1,bound=-1
      ! ... Q'(0)=1,Q(0)=0,Q(1)=0 for side=0,bound=-1
      ! ... Q'(1)=1,Q(0)=0,Q(1)=0 for side=1,bound=1
      ! ... Q'(0)=1,Q(0)=0,Q(1)=0 for side=0,bound=1
      !
      USE kinds,         ONLY : DP
      !
      IMPLICIT NONE
      !
      REAL( DP ) :: qinterp
      REAL( DP ) :: x
      INTEGER    :: side
      INTEGER    :: bound
      !
      IF( bound == 0 .AND. side == 1 ) THEN
        qinterp = - x * x + x * x * x
      ELSE IF( bound == 0 .AND. side == 0 ) THEN
        qinterp = x - 2.D0 * x * x + x * x * x
      ELSE IF( bound == - 1 .AND. side == 0 ) THEN
        qinterp = 0.D0
      ELSE IF( bound == 1 .AND. side == 1 ) THEN
        qinterp =  0.D0
      ELSE IF( bound == - 1 .AND. side == 1 ) THEN
        qinterp = - x + x * x
      ELSE IF( bound == 1 .AND. side == 0 ) THEN
        qinterp = x - x * x
      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION qinterp
!--------------------------------------------------------------------

!--------------------------------------------------------------------
      FUNCTION dqinterp( x, side, bound )
!--------------------------------------------------------------------
      !
      ! ... Calculates the **DERIVATIVE** of Q(x) where Qinterp
      USE kinds,         ONLY : DP
      !
      IMPLICIT NONE
      !
      REAL( DP ) :: dqinterp
      REAL( DP ) :: x
      INTEGER    :: side
      INTEGER    :: bound
      !
      IF( bound == 0 .AND. side == 1 ) THEN
        dqinterp = - 2 * x  + 3 * x * x
      ELSE IF( bound == 0 .AND. side == 0 ) THEN
        dqinterp = 1 - 4.D0 * x  + 3 * x * x
      ELSE IF( bound == - 1 .AND. side == 0 ) THEN
        dqinterp = 0.D0
      ELSE IF( bound == 1 .AND. side == 1 ) THEN
        dqinterp =  0.D0
      ELSE IF( bound == - 1 .AND. side == 1 ) THEN
        dqinterp = - 1 + 2 * x
      ELSE IF( bound == 1 .AND. side == 0 ) THEN
        dqinterp = 1 - 2 * x
      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION dqinterp
!--------------------------------------------------------------------


!--------------------------------------------------------------------
      FUNCTION bound( i, n )
!--------------------------------------------------------------------
      !
      ! ... Returns -1 if i = 0, 1 if i = n - 1, and 0 otherwise
      !
      IMPLICIT NONE
      !
      INTEGER :: i
      INTEGER :: n
      INTEGER :: bound
      !
      IF( i == 1 ) THEN
        bound = -1
      ELSE IF( i == n - 1 ) THEN
        bound = 1
      ELSE
        bound = 0
      END IF
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION bound
!--------------------------------------------------------------------
