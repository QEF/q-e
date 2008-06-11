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
      SUBROUTINE writetofile( f, filename, mr1, mr2, mr3,              &
                              delta1m, delta2m, delta3m, which_print )
!--------------------------------------------------------------------
      !
      USE kinds,         ONLY : DP
      !
      IMPLICIT NONE
      !
      INTEGER                 :: mr1, mr2, mr3
      REAL( DP )              :: delta1m, delta2m, delta3m
      REAL( DP )              :: f( mr1 * mr2 * mr3 )
      !
      CHARACTER( LEN = 256 )  :: filename
      CHARACTER( LEN = 256 )  :: which_print
      !
      INTEGER                 :: ir, ir1, ir2, ir3, num
      !
      REAL( DP )              :: total
      !
      INTEGER, EXTERNAL       :: compindex
      !

      OPEN( 300, file = TRIM( filename ), status = 'unknown' )
      !
      SELECT CASE( TRIM( which_print ) )
        !
      CASE( 'all' )
        !
        WRITE( 300, * ) f( : )
        !
      CASE( 'x' )
        !
        DO ir1 = 1, mr1
          total = 0.D0
          num = 0
          DO ir2 = mr2 / 2, mr2 / 2
            DO ir3 = mr3 / 2, mr3 / 2
              ir = compindex( ir1, ir2, ir3, mr1, mr2, mr3 )
              total = total + f( ir )
              num = num + 1
            END DO
          END DO
          total = total / DBLE( num )
          WRITE( 300, '(2E30.10)' ) DBLE( ir1 - 1 )   &
               * delta1m, total
        END DO
        !
       CASE( 'y' )
        !
        DO ir2 = 1, mr2
          total = 0.D0
          num = 0
          DO ir1 = mr1 / 2, mr1 / 2
            DO ir3 = mr3 / 2, mr3 / 2
              ir = compindex( ir1, ir2, ir3, mr1, mr2, mr3 )
              total = total + f( ir )
              num = num + 1
            END DO
          END DO
          total = total / DBLE( num )
          WRITE( 300, '(2E30.10)' ) DBLE( ir2 - 1 )   &
               * delta2m, total
        END DO
        !
      CASE( 'z' )
        !
        DO ir3 = 1, mr3
          total = 0.D0
          num = 0
          DO ir2 = mr2 / 2, mr2 / 2
            DO ir1 = mr1 / 2, mr1 / 2
              ir = compindex( ir1, ir2, ir3, mr1, mr2, mr3 )
              total = total + f( ir )
              num = num + 1
            END DO
          END DO
          total = total / DBLE( num )
          WRITE( 300, '(2E30.10)' ) DBLE( ir3 - 1 )   &
               * delta3m, total
        END DO
        !
      CASE( 'ax' )
        !
        DO ir1 = 1, mr1
          total = 0.D0
          num = 0
          DO ir2 = 1, mr2
            DO ir3 = 1, mr3
              ir = compindex( ir1, ir2, ir3, mr1, mr2, mr3 )
              total = total + f( ir )
              num = num + 1
            END DO
          END DO
          total = total / DBLE( num )
          WRITE( 300, '(2E30.10)' ) DBLE( ir1 - 1 )   &
               * delta1m, total
        END DO
        !
       CASE( 'ay' )
        !
        DO ir2 = 1, mr2
          total = 0.D0
          num = 0
          DO ir1 = 1, mr1
            DO ir3 = 1, mr3
              ir = compindex( ir1, ir2, ir3, mr1, mr2, mr3 )
              total = total + f( ir )
              num = num + 1
            END DO
          END DO
          total = total / DBLE( num )
          WRITE( 300, '(2E30.10)' ) DBLE( ir2 - 1 )   &
               * delta2m, total
        END DO
        !
      CASE( 'az' )
        !
        DO ir3 = 1, mr3
          total = 0.D0
          num = 0
          DO ir2 = 1, mr2
            DO ir1 = 1, mr1
              ir = compindex( ir1, ir2, ir3, mr1, mr2, mr3 )
              total = total + f( ir )
              num = num + 1
            END DO
          END DO
          total = total / DBLE( num )
          WRITE( 300, '(2E30.10)' ) DBLE( ir3 - 1 )   &
               * delta3m, total
        END DO
        !
      END SELECT
      !
      CLOSE( 300 )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE writetofile
!--------------------------------------------------------------------
