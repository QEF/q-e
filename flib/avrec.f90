!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

    SUBROUTINE avrec( n, alpha, v, av )

! ... This subroutine try to use fast library to
! ... calculate 
! ...             av(i) = alpha / v(i)
! ...

      USE kinds, ONLY : DP
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      INTEGER :: i
      REAL(DP), INTENT(IN) :: alpha 
      REAL(DP), INTENT(IN) :: v(*)
      REAL(DP), INTENT(OUT) :: av(*)

#if defined __BENCHLIB

      CALL oneover_v( n, v, av )
      IF( alpha /= 1.0d0 ) THEN
        CALL dscal( n, alpha, av, 1 )
      END IF

#elif defined __MASS

      CALL vrec( av, v, n )
      IF( alpha /= 1.0d0 ) THEN
        CALL dscal( n, alpha, av, 1 )
      END IF

#else

      DO i = 1, n
        av(i) = alpha / v(i)
      END DO 

#endif

      RETURN
    END SUBROUTINE avrec
