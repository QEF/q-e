!
! Copyright (C) 2001-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE sph_dbes( MMAX, R, XG, L, DJL )
  !----------------------------------------------------------------------------
  !
  ! ... calculates derivatives of spherical bessel functions  j_l(Gr)
  ! ... with respect to h_alpha,beta (without the factor GAGK(KK,IG)*HTM1)
  ! ... i.e. -x * D(jl(x))/dx
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  !
  IMPLICIT NONE
  !
  INTEGER       :: MMAX, L
  REAL(DP) :: XG
  REAL(DP) :: DJL(MMAX), R(MMAX)
  !
  INTEGER       :: IR
  REAL(DP) :: XRG, XRG2
  !
  !
  IF ( L == 1 ) THEN                      !   S  PART
     IF( XG < eps8 ) THEN
        DO IR=1,MMAX
           DJL(IR) = 0.D0
        END DO
     ELSE
        DJL(1) = 0.D0
        DO IR=2,MMAX
           XRG=R(IR)*XG
           DJL(IR) = SIN(XRG)/XRG-COS(XRG)
        END DO
     ENDIF
  ENDIF
  !
  IF ( L == 2 ) THEN                      !   P  PART
     IF( XG < eps8 ) THEN
        DO IR=1,MMAX
           DJL(IR) = 0.D0
        END DO
     ELSE
        DJL(1) = 0.D0
        DO IR=2,MMAX
           XRG=R(IR)*XG
           DJL(IR) = 2.D0*(SIN(XRG)/XRG-COS(XRG))/XRG - SIN(XRG)
        END DO
     ENDIF
  ENDIF
  !
  IF ( L == 3 ) THEN                      !   D  PART
     IF ( XG < eps8 ) THEN
        DO IR=1,MMAX
           DJL(IR) = 0.D0
        END DO
     ELSE
        DJL(1) = 0.D0
        DO IR=2,MMAX
           XRG=R(IR)*XG
           DJL(IR) = ( SIN(XRG)*(9.D0/(XRG*XRG)-4.D0) - &
                       9.D0*COS(XRG)/XRG ) /XRG + COS(XRG)
        END DO
     END IF
  END IF
  !
  IF ( L == 4 ) THEN                      !   F  PART
     IF ( XG < eps8 ) THEN
        DO IR=1,MMAX
           DJL(IR) = 0.D0
        END DO
     ELSE
        DJL(1) = 0.D0
        DO IR=2,MMAX
           XRG=R(IR)*XG
           XRG2=XRG*XRG
           DJL(IR) = SIN(XRG)*(60.D0/(XRG2*XRG2)-27.D0/XRG2+1.d0) - &
                     COS(XRG)*(60.D0/XRG2-7.D0)/XRG
        END DO
     END IF
  END IF
  !
  IF ( L == 5 ) THEN                      !   G  PART
     IF ( XG < eps8 ) THEN
        DO IR=1,MMAX
           DJL(IR) = 0.D0
        END DO
     ELSE
        DJL(1) = 0.D0
        DO IR=2,MMAX
           XRG=R(IR)*XG
           XRG2=XRG*XRG
           DJL(IR) = SIN(XRG)*(525.D0/(XRG2*XRG2)-240.D0/XRG2+11.D0)/XRG - &
                     COS(XRG)*(525.D0/(XRG2*XRG2)-65.D0/XRG2+1.D0)
        END DO
     END IF
  END IF
  !
  IF ( L <= 0 .OR. L >= 6 ) &
     CALL errore( 'sph_dbes', ' L NOT PROGRAMMED, L= ',L )
  !
  RETURN
  !
END SUBROUTINE sph_dbes

!
SUBROUTINE sph_dbes1 ( nr, r, xg, l, jl, djl )
  !
  ! calculates x*dj_l(x)/dx using the recursion formula
  !    dj_l(x)/dx = l/x*j_l(x) -  j_(l+1)(x)
  ! for l=0, and for l>0 :
  !    dj_l(x)/dx = j_(l-1)(x) - (l+1)/x *  j_l(x)
  ! requires j_l(r) in input
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l, nr
  REAL (DP), INTENT(IN) :: xg, jl(nr), r(nr)
  REAL (DP), INTENT(OUT):: djl(nr)
  !
  if ( xg < eps8 ) then
     !
     ! special case q=0
     ! note that x*dj_l(x)/dx = 0 for x = 0
     !
     djl(:) = 0.0d0
  else
     !
     if ( l > 0 ) then
        call sph_bes ( nr, r, xg, l-1, djl )
        djl(:) = djl(:) * (xg * r(:) ) - (l+1) * jl(:)
     else if ( l == 0 ) then
        call sph_bes ( nr, r, xg, l+1, djl )
        djl(:) = - djl(:) * (xg * r(:) )
     else
        call errore('sph_dbes','l < 0 not implemented', abs(l) )
     end if
  end if
  !
end SUBROUTINE sph_dbes1
