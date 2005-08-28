!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
SUBROUTINE sph_besr( msh, r, q, l, jl )  
  !--------------------------------------------------------------------
  !
  ! input:
  !   msh     = number of grid points points
  !   r(1:msh)= radial grid
  !   q       = q
  !   l       = angular momentum (-1 <= l <= 6)
  ! output:
  !   jl(1:msh) = j_l(q*r(i))  (j_l = spherical bessel function)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  !
  IMPLICIT NONE  
  !
  INTEGER :: msh, l  
  REAL(DP) :: r (msh), q, jl (msh)  
  !
  INTEGER :: ir, ir0  
  !
  !
  IF (ABS (q) < eps8) THEN  
     IF (l == -1) THEN  
        CALL errore ('sph_besr', 'j_{-1}(0) ?!?', 1)  
     ELSEIF (l == 0) THEN  
        DO ir = 1, msh  
           jl (ir) = r(ir)   
        ENDDO
     ELSE  
        jl(:) = 0.d0
     ENDIF
  ELSE  
     IF (ABS (q * r (1) ) > eps8) THEN  
        ir0 = 1  
     ELSE  
        IF (l == -1) THEN  
           CALL errore ('sph_besr', 'j_{-1}(0) ?!?', 2)  
        ELSEIF (l == 0) THEN  
           jl (1) = r(1)  
        ELSE  
           jl (1) = 0.d0  
        ENDIF
        ir0 = 2 
     ENDIF
     IF (l == - 1) THEN  
        DO ir = ir0, msh  
           jl (ir) = COS (q * r (ir) ) / q   
        ENDDO
     ELSEIF (l == 0) THEN  
        DO ir = ir0, msh  
           jl (ir) = SIN (q * r (ir) ) / q  
        ENDDO
     ELSEIF (l == 1) THEN  
        DO ir = ir0, msh  
           IF (ABS(q*r(ir)).LT.6.d-2) THEN
              jl(ir)=(q*r(ir))/3.d0-(q*r(ir))**3/30.d0+ &
                     (q*r(ir))**5/840.d0-(q*r(ir))**7/45360.d0
              jl(ir) =  jl(ir)*r(ir)
           ELSE
              jl(ir)=(SIN(q*r(ir))/(q*r(ir))-COS(q*r(ir)))/q
           ENDIF
        ENDDO
     ELSEIF (l == 2) THEN  
        DO ir = ir0, msh  
           IF (ABS(q*r(ir)).LT.6.d-2) THEN
              jl(ir) = (q*r(ir))**2/15.d0-(q*r(ir))**4/210.d0 &
                     +(q*r(ir))**6/7560.d0
              jl(ir) =  jl(ir)*r(ir)
           ELSE
              jl(ir) = ( ( 3.d0-(q*r(ir))**2)*SIN(q*r(ir)) &
               - 3.d0 * q*r(ir)*COS(q*r(ir)) )/(q*r(ir))**2/q
           ENDIF
        ENDDO
     ELSEIF (l == 3) THEN  
        DO ir = ir0, msh  
           IF (ABS(q*r(ir)).LT.6.d-2) THEN
               jl(ir) = 48.d0*(q*r(ir))**3/5040.d0 &
                      - 192.d0*(q*r(ir))**5/362880.d0
               jl(ir) = jl(ir) * r(ir)
           ELSE
              jl(ir) = (SIN(q*r(ir)) * (15.d0/(q*r(ir))      &
                - 6.d0*(q*r(ir))) + COS(q*r(ir))*((q*r(ir))**2 &
                - 15.d0) ) / (q * r (ir) ) **2/q
           ENDIF
        ENDDO
     ELSEIF (l == 4) THEN  
        DO ir = ir0, msh  
           jl (ir) = (SIN (q * r (ir) ) * (105.d0 - 45.d0 * (q * r (ir) &
                ) **2 + (q * r (ir) ) **4) + COS (q * r (ir) ) * (10.d0 * &
                (q * r (ir) ) **3 - 105.d0 * (q * r (ir) ) ) ) / (q * r (ir) &
                ) **4/q
        ENDDO
     ELSEIF (l == 5) THEN  
        DO ir = ir0, msh  
           jl (ir) = (-COS(q*r(ir)) - (945.d0*COS(q*r(ir)))/(q*r(ir)) ** 4 + &
                (105.d0*COS(q*r(ir)))/ (q*r(ir)) ** 2 + &
                (945.d0*SIN(q*r(ir)))/ (q*r(ir)) ** 5 - &
                (420.d0*SIN(q*r(ir)))/(q*r(ir)) ** 3 +  &
                (15.d0*SIN(q*r(ir)))/ (q*r(ir)))/q
        END DO
     ELSEIF (l == 6) THEN  
        DO ir = ir0, msh  
           jl (ir) = ((-10395.d0*COS(q*r(ir)))/(q*r(ir)) ** 5 + &
                (1260.d0*COS(q*r(ir)))/(q*r(ir)) ** 3 - &
                (21.d0*COS(q*r(ir)))/ (q*r(ir)) - SIN(q*r(ir)) + &
                (10395.d0*SIN(q*r(ir)))/(q*r(ir)) ** 6 - &
                (4725.d0*SIN(q*r(ir)))/ (q*r(ir)) ** 4 + &
                (210.d0*SIN(q*r(ir)))/(q*r(ir)) ** 2)/q
        END DO
     ELSE  
        CALL errore ('sph_besr', 'not implemented', l)  
     ENDIF
  ENDIF
  !
  RETURN  
END SUBROUTINE sph_besr
