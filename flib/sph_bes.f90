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
  ! ...   msh     = number of grid points points
  ! ...   r(1:msh)= radial grid
  ! ...   q       = q
  ! ...   l       = angular momentum (-1 <= l <= 6)
  ! ... output:
  ! ...   jl(1:msh) = j_l(q*r(i))  (j_l = spherical bessel function)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  !
  IMPLICIT NONE
  !
  INTEGER :: msh, l
  REAL(kind=DP) :: r(msh), q, jl(msh)
  !
  INTEGER :: ir, ir0
  !
#if defined (__MASS)
  REAL(kind=DP) :: qr(msh), sin_qr(msh), cos_qr(msh)
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
        CALL vcos( cos_qr, qr, msh)
        !
        DO ir = ir0, msh
           jl (ir) = cos_qr(ir) / (q * r (ir) )
        END DO
         !
#else
        !
        DO ir = ir0, msh
           jl (ir) = COS (q * r (ir) ) / (q * r (ir) )
        END DO
        !
#endif
        !
     ELSE IF ( l == 0 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vsin( sin_qr, qr, msh)
        !
        DO ir = ir0, msh
           jl (ir) = sin_qr(ir) / (q * r (ir) )
        END DO
        !
#else
        !
        DO ir = ir0, msh
           jl (ir) = SIN (q * r (ir) ) / (q * r (ir) )
        END DO
        !
#endif
        !
     ELSE IF ( l == 1 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh)
        CALL vsin( sin_qr, qr, msh)
        !
        DO ir = ir0, msh
           jl (ir) = ( sin_qr(ir) / (q * r (ir) ) - cos_qr(ir) ) / (q * r (ir) )
        END DO
        !
#else
        !
        DO ir = ir0, msh
           jl (ir) = (SIN (q * r (ir) ) / (q * r (ir) ) - COS (q * r ( &
                ir) ) ) / (q * r (ir) )
        END DO
        !
#endif
        !
     ELSE IF ( l == 2 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh)
        CALL vsin( sin_qr, qr, msh)
        !
        DO ir = ir0, msh
           jl (ir) = ( (3.d0 / (q * r (ir) ) - (q * r (ir) ) ) * sin_qr( &
                ir ) - 3.d0 * cos_qr(ir) ) / (q * r (ir) ) ** 2
        END DO
        !
#else
        !
        DO ir = ir0, msh
           jl (ir) = ( (3.d0 / (q * r (ir) ) - (q * r (ir) ) ) * SIN ( &
                q * r (ir) ) - 3.d0 * COS (q * r (ir) ) ) / (q * r (ir) ) ** 2
        END DO
        !
#endif
        !
     ELSE IF ( l == 3 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh)
        CALL vsin( sin_qr, qr, msh)
        !
        DO ir = ir0, msh
           jl (ir) = (sin_qr (ir) * (15.d0 / (q * r (ir) ) &
                - 6.d0 * (q * r (ir) ) ) + cos_qr (ir) * ( (q * r (ir) &
                ) **2 - 15.d0) ) / (q * r (ir) ) **3
        ENDDO
        !
#else
        !
        DO ir = ir0, msh
           jl (ir) = (SIN (q * r (ir) ) * (15.d0 / (q * r (ir) ) &
                - 6.d0 * (q * r (ir) ) ) + COS (q * r (ir) ) * ( (q * r (ir) &
                ) **2 - 15.d0) ) / (q * r (ir) ) ** 3
        END DO
        !
#endif
        !
     ELSE IF ( l == 4 ) THEN
        !
#if defined (__MASS)
        !
        qr = q * r
        !
        CALL vcos( cos_qr, qr, msh)
        CALL vsin( sin_qr, qr, msh)
        !
        DO ir = ir0, msh
           jl (ir) = (sin_qr (ir) * (105.d0 - 45.d0 * (q * r (ir) &
                ) **2 + (q * r (ir) ) **4) + cos_qr (ir) * (10.d0 * &
                (q * r (ir) ) **3 - 105.d0 * (q * r (ir) ) ) ) / (q * r (ir) &
                ) **5
        END DO
        !
#else
        !
        DO ir = ir0, msh
           jl (ir) = (SIN (q * r (ir) ) * (105.d0 - 45.d0 * (q * r (ir) &
                ) **2 + (q * r (ir) ) **4) + COS (q * r (ir) ) * (10.d0 * &
                (q * r (ir) ) **3 - 105.d0 * (q * r (ir) ) ) ) / (q * r (ir) &
                ) **5
        END DO
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
        DO ir = ir0, msh
           jl (ir) = (-cos_qr(ir) - (945.d0*cos_qr(ir))/(q*r(ir)) ** 4 + &
                (105.d0*cos_qr(ir))/ (q*r(ir)) ** 2 + &
                (945.d0*sin_qr(ir))/ (q*r(ir)) ** 5 - &
                (420.d0*sin_qr(ir))/(q*r(ir)) ** 3 + &
                (15.d0*sin_qr(ir))/ (q*r(ir)))/(q*r(ir))
        END DO
        !
#else
        !
        DO ir = ir0, msh
           jl (ir) = (-COS(q*r(ir)) - (945.d0*COS(q*r(ir)))/(q*r(ir)) ** 4 + &
                (105.d0*COS(q*r(ir)))/ (q*r(ir)) ** 2 + &
                (945.d0*SIN(q*r(ir)))/ (q*r(ir)) ** 5 - &
                (420.d0*SIN(q*r(ir)))/(q*r(ir)) ** 3 + &
                (15.d0*SIN(q*r(ir)))/ (q*r(ir)))/(q*r(ir))
        END DO
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
        DO ir = ir0, msh
           jl (ir) = ((-10395.d0*cos_qr(ir))/(q*r(ir)) ** 5 + &
                (1260.d0*cos_qr(ir))/(q*r(ir)) ** 3 - &
                (21.d0*cos_qr(ir))/ (q*r(ir)) - sin_qr(ir) + &
                (10395.d0*sin_qr(ir))/(q*r(ir)) ** 6 - &
                (4725.d0*sin_qr(ir))/ (q*r(ir)) ** 4 + &
                (210.d0*sin_qr(ir))/(q*r(ir)) ** 2)/(q*r(ir))
        END DO
        !
#else
        !
        DO ir = ir0, msh
           jl (ir) = ((-10395.d0*COS(q*r(ir)))/(q*r(ir)) ** 5 + &
                (1260.d0*COS(q*r(ir)))/(q*r(ir)) ** 3 - &
                (21.d0*COS(q*r(ir)))/ (q*r(ir)) - SIN(q*r(ir)) + &
                (10395.d0*SIN(q*r(ir)))/(q*r(ir)) ** 6 - &
                (4725.d0*SIN(q*r(ir)))/ (q*r(ir)) ** 4 + &
                (210.d0*SIN(q*r(ir)))/(q*r(ir)) ** 2)/(q*r(ir))
        END DO
        !
#endif
        !
     ELSE
        !
        CALL errore ('sph_bes', 'not implemented', l)
        !
     END IF
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE sph_bes
