!
! Copyright (C) 2001-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE sph_bes( msh, r, q, l, jl )
  !----------------------------------------------------------------------------
  !
  ! ... input:
  ! ...   msh      = number of grid points points
  ! ...   r(1:msh) = radial grid
  ! ...   q        = q
  ! ...   l        = angular momentum ( -1 <= l <= 6 )
  ! ... output:
  ! ...   jl(1:msh) = j_l(q*r(i))  (j_l = spherical bessel function)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  !
  IMPLICIT NONE
  !
  INTEGER       :: msh, l
  REAL(kind=DP) :: r(msh), q, jl(msh)
  !
  INTEGER :: ir0
  !
#if defined (__MASS)
  !
  REAL(kind=DP) :: qr(msh), sin_qr(msh), cos_qr(msh)
  !
#endif
  !
  IF ( ABS( q ) < eps8 ) THEN
     !
     IF ( l == -1 ) THEN
        !
        CALL errore( 'sph_bes', 'j_{-1}(0) ?!?', 1 )
        !
     ELSE IF ( l == 0 ) THEN
        !
        jl(:) = 1.D0
        !
     ELSE
        !
        jl(:) = 0.D0
        !
     END IF
     !
  ELSE
     !
     IF ( ABS( q*r(1) ) > eps8 ) THEN
        !
        ir0 = 1
        !
     ELSE
        !
        IF (l == -1) THEN
           !
           CALL errore( 'sph_bes', 'j_{-1}(0) ?!?', 2 )
           !
        ELSE IF ( l == 0 ) THEN
           !
           jl(1) = 1.D0
           !
        ELSE
           !
           jl(1) = 0.D0
           !
        ENDIF
        !
        ir0 = 2
        !
     END IF
     !
     IF ( l == - 1 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh )
        !
        jl(ir0:) = cos_qr(ir0:) / ( q * r(ir0:) )
         !
#else
        !
        jl(ir0:) = COS( q * r(ir0:) ) / ( q * r(ir0:) )
        !
#endif
        !
     ELSE IF ( l == 0 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vsin( sin_qr, qr, msh )
        !
        jl(ir0:) = sin_qr(ir0:) / ( q * r(ir0:) )
        !
#else
        !
        jl(ir0:) = SIN( q * r(ir0:) ) / ( q * r(ir0:) )
        !
#endif
        !
     ELSE IF ( l == 1 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh )
        CALL vsin( sin_qr, qr, msh )
        !
        jl(ir0:) = ( sin_qr(ir) / ( q * r(ir0:) ) - &
                     cos_qr(ir0:) ) / ( q * r(ir0:) )
        !
#else
        !
        jl(ir0:) = ( SIN( q * r(ir0:) ) / ( q * r(ir0:) ) - &
                     COS( q * r(ir0:) ) ) / ( q * r(ir0:) )
        !
#endif
        !
     ELSE IF ( l == 2 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh )
        CALL vsin( sin_qr, qr, msh )
        !
        jl(ir0:) = ( ( 3.D0 / ( q * r(ir0:) ) - ( q * r(ir0:) ) ) * &
                     sin_qr(ir0:) - 3.D0 * cos_qr(ir0:) ) / (q * r(ir0:) ) ** 2
        !
#else
        !
        jl(ir0:) = ( ( 3.D0 / ( q * r(ir0:) ) - ( q * r(ir0:) ) ) * &
                     SIN( q * r(ir0:) ) - 3.D0 * COS( q * r(ir0:) ) ) / &
                   ( q * r(ir0:) ) ** 2
        !
#endif
        !
     ELSE IF ( l == 3 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh )
        CALL vsin( sin_qr, qr, msh )
        !
        jl(ir0:) = ( sin_qr(ir0:) * ( 15.D0 / ( q * r(ir0:) ) - &
                   6.D0 * ( q * r(ir0:) ) ) + cos_qr(ir0:) * &
                   ( ( q * r(ir0:) ) ** 2 - 15.D0 ) ) / ( q * r(ir0:) ) ** 3
        !
#else
        !
        jl(ir0:) = ( SIN( q * r(ir0:) ) * ( 15.D0 / ( q * r(ir0:) ) - &
                   6.D0 * ( q * r(ir0:) ) ) + COS( q * r(ir0:) ) * &
                   ( ( q * r(ir0:) ) ** 2 - 15.D0 ) ) / ( q * r(ir0:) ) ** 3
        !
#endif
        !
     ELSE IF ( l == 4 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh )
        CALL vsin( sin_qr, qr, msh )
        !
        jl(ir0:) = ( sin_qr(ir0:) * ( 105.D0 - 45.D0 * ( q * r(ir0:) ) ** 2 + &
                   ( q * r(ir0:) ) ** 4 ) + cos_qr(ir0:) * ( 10.D0 * &
                   ( q * r(ir0:) ) ** 3 - 105.D0 * ( q * r(ir0:) ) ) ) / &
                   ( q * r(ir0:) ) ** 5
        !
#else
        !
        jl(ir0:) = ( SIN( q * r(ir0:) ) * ( 105.D0 - 45.D0 * &
                   ( q * r(ir0:) ) ** 2 + ( q * r(ir0:) ) **4 ) + &
                   COS( q * r(ir0:) ) * ( 10.D0 * ( q * r(ir0:) ) ** 3 - &
                   105.D0 * ( q * r(ir0:) ) ) ) / ( q * r(ir0:) ) ** 5
        !
#endif
        !
     ELSE IF ( l == 5 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh )
        CALL vsin( sin_qr, qr, msh )
        !
        jl(ir0:) = ( - cos_qr(ir0:) - ( 945.D0 * cos_qr(ir0:) ) / &
                   ( q * r(ir0:) ) ** 4 + ( 105.D0 * cos_qr(ir0:) ) / &
                   ( q * r(ir0:) ) ** 2 + ( 945.D0 * sin_qr(ir0:) ) / &
                   ( q * r(ir0:) ) ** 5 - ( 420.D0 * sin_qr(ir0:) ) / &
                   ( q * r(ir0:) ) ** 3 + (  15.D0 * sin_qr(ir0:) ) / &
                   ( q * r(ir0:) ) ) / ( q * r(ir0:) )
        !
#else
        !
        jl(ir0:) = ( - COS( q * r(ir0:) ) - ( 945.D0 * COS( q * r(ir0:) ) ) / &
                   ( q * r(ir0:) ) ** 4 + ( 105.D0 * COS( q * r(ir0:) ) ) / &
                   ( q * r(ir0:) ) ** 2 + ( 945.D0 * SIN( q * r(ir0:) ) ) / &
                   ( q * r(ir0:) ) ** 5 - ( 420.D0 * SIN( q * r(ir0:) ) ) / &
                   ( q * r(ir0:) ) ** 3 + (  15.D0 * SIN( q * r(ir0:) ) ) / &
                   ( q * r(ir0:) ) ) / ( q * r(ir0:) )
        !
#endif
        !
     ELSE IF ( l == 6 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh )
        CALL vsin( sin_qr, qr, msh )
        !
        jl(ir0:) = ( - sin_qr(ir0:) + &
                     (     21.D0 * cos_qr(ir0:) ) / ( q * r(ir0:) ) + &
                     ( -10395.D0 * cos_qr(ir0:) ) / ( q * r(ir0:) ) ** 5 + &
                     (   1260.D0 * cos_qr(ir0:) ) / ( q * r(ir0:) ) ** 3 - &
                     (  10395.D0 * sin_qr(ir0:) ) / ( q * r(ir0:) ) ** 6 - &
                     (   4725.D0 * sin_qr(ir0:) ) / ( q * r(ir0:) ) ** 4 + &
                     (    210.D0 * sin_qr(ir0:) ) / ( q * r(ir0:) ) ** 2 ) / &
                   ( q * r(ir0:) )
        !
#else
        !
        jl(ir0:) = ( - SIN( q * r(ir0:) ) + &
                     (     21.D0 * COS( q * r(ir0:) ) ) / ( q * r(ir0:) ) + &
                     ( -10395.D0 * COS( q * r(ir0:) ) ) / ( q * r(ir0:) )**5 + &
                     (   1260.D0 * COS( q * r(ir0:) ) ) / ( q * r(ir0:) )**3 - &
                     (  10395.D0 * SIN( q * r(ir0:) ) ) / ( q * r(ir0:) )**6 - &
                     (   4725.D0 * SIN( q * r(ir0:) ) ) / ( q * r(ir0:) )**4 + &
                     (    210.D0 * SIN( q * r(ir0:) ) ) / ( q * r(ir0:) )**2 ) / &
                   ( q * r(ir0:) )
        !
#endif
        !
     ELSE
        !
        CALL errore( 'sph_bes', 'not implemented', l )
        !
     END IF
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE sph_bes
