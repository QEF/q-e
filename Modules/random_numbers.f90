!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE random_numbers
  !----------------------------------------------------------------------------
  !! Module for random numbers generation.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTERFACE gauss_dist
     !
     MODULE PROCEDURE gauss_dist_scal, gauss_dist_vect
     !
  END INTERFACE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    FUNCTION randy ( irand )
      !------------------------------------------------------------------------
      !! * x=randy(n): reseed with initial seed \(\text{idum}=n\) ( 0 <= n <= ic, see below)
      !!               if randy is not explicitly initialized, it will be
      !!               initialized with seed \(\text{idum}=0\) the first time it is called;
      !! * x=randy() : generate uniform REAL(DP) numbers x in [0,1].
      !
      REAL(DP) :: randy
      INTEGER, optional    :: irand
      !
      INTEGER , PARAMETER  :: m    = 714025, &
                              ia   = 1366, &
                              ic   = 150889, &
                              ntab = 97
      REAL(DP), PARAMETER  :: rm = 1.0_DP / m
      INTEGER              :: j
      INTEGER, SAVE        :: ir(ntab), iy, idum=0
      LOGICAL, SAVE        :: first=.true.
      !
      IF ( present(irand) ) THEN
         idum = MIN( ABS(irand), ic) 
         first=.true.
      END IF

      IF ( first ) THEN
         !
         first = .false.
         idum = MOD( ic - idum, m )
         !
         DO j=1,ntab
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         END DO
         idum=mod(ia*idum+ic,m)
         iy=idum
      END IF
      j=1+(ntab*iy)/m
      IF( j > ntab .OR. j <  1 ) call errore('randy','j out of range',ABS(j)+1)
      iy=ir(j)
      randy=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      !
      RETURN
      !
    END FUNCTION randy
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_random_seed ( )
      !------------------------------------------------------------------------
      !! poor-man random seed for \(\texttt{randy}\).
      !
      INTEGER, DIMENSION (8) :: itime
      INTEGER :: iseed, irand
      !
      CALL date_and_time ( values = itime ) 
      ! itime contains: year, month, day, time difference (minutes) from UTC,
      !                 hours, minutes, seconds and milliseconds. 
      ! The following rather arbitrary choice, modified as suggested by
      ! Han Hsu, appears to yield sufficiently randomized initial seeds
      !
      iseed = ( itime(8) + itime(6) ) * ( itime(7) + itime(5) )
      irand = randy ( iseed )
      !
    END SUBROUTINE set_random_seed
    !
    !-----------------------------------------------------------------------
    FUNCTION gauss_dist_scal( mu, sigma )
      !-----------------------------------------------------------------------
      !! This function generates a number taken from a normal distribution of
      !! mean value \(\text{mu}\) and variance \(\text{sigma}.
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: mu
      REAL(DP), INTENT(IN) :: sigma
      REAL(DP)             :: gauss_dist_scal
      !
      REAL(DP) :: x1, x2, w
      !
      !
      gaussian_loop: DO
         !
         x1 = 2.0_DP * randy() - 1.0_DP
         x2 = 2.0_DP * randy() - 1.0_DP
         !
         w = x1 * x1 + x2 * x2
         !
         IF ( w < 1.0_DP ) EXIT gaussian_loop
         !
      END DO gaussian_loop
      !
      w = SQRT( ( - 2.0_DP * LOG( w ) ) / w )
      !
      gauss_dist_scal = x1 * w * sigma + mu
      !
      RETURN
      !
    END FUNCTION gauss_dist_scal
    !
    !-----------------------------------------------------------------------
    FUNCTION gauss_dist_cmplx( mu, sigma )
      !-----------------------------------------------------------------------
      !! This function generates a number taken from a normal distribution of
      !! mean value \(\text{mu}\) and variance \(\text{sigma}\).
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: mu
      REAL(DP), INTENT(IN) :: sigma
      COMPLEX(DP)          :: gauss_dist_cmplx
      !
      REAL(DP) :: x1, x2, w
      !
      !
      gaussian_loop: DO
         !
         x1 = 2.0_DP * randy() - 1.0_DP
         x2 = 2.0_DP * randy() - 1.0_DP
         !
         w = x1 * x1 + x2 * x2
         !
         IF ( w < 1.0_DP ) EXIT gaussian_loop
         !
      END DO gaussian_loop
      !
      w = SQRT( ( - 2.0_DP * LOG( w ) ) / w )
      !
      gauss_dist_cmplx = CMPLX( x1 * w * sigma + mu, x2 * w * sigma + mu, kind=DP)
      !
      RETURN
      !
    END FUNCTION gauss_dist_cmplx
    !    
    !-----------------------------------------------------------------------
    FUNCTION gauss_dist_vect( mu, sigma, dim )
      !-----------------------------------------------------------------------
      !! This function generates an array of numbers taken from a normal
      !! distribution of mean value \(\text{mu}\) and variance \(\text{sigma}\).
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: mu
      REAL(DP), INTENT(IN) :: sigma
      INTEGER,  INTENT(IN) :: dim
      REAL(DP)             :: gauss_dist_vect( dim )
      !
      REAL(DP) :: x1, x2, w
      INTEGER  :: i
      !
      !
      DO i = 1, dim, 2
         !
         gaussian_loop: DO
            !
            x1 = 2.0_DP * randy() - 1.0_DP
            x2 = 2.0_DP * randy() - 1.0_DP
            !
            w = x1 * x1 + x2 * x2
            !
            IF ( w < 1.0_DP ) EXIT gaussian_loop
            !
         END DO gaussian_loop
         !
         w = SQRT( ( - 2.0_DP * LOG( w ) ) / w )
         !
         gauss_dist_vect(i) = x1 * w * sigma
         !
         IF ( i >= dim ) EXIT
         !
         gauss_dist_vect(i+1) = x2 * w * sigma
         !
      END DO
      !
      gauss_dist_vect(:) = gauss_dist_vect(:) + mu
      !
      RETURN
      !
    END FUNCTION gauss_dist_vect
    !
    !-----------------------------------------------------------------------
    FUNCTION gamma_dist (ialpha)
      !-----------------------------------------------------------------------
      !! Gamma-distributed random number, implemented as described in
      !! Numerical recipes (Press, Teukolsky, Vetterling, Flannery).
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ialpha
      REAL(DP) gamma_dist
      INTEGER j
      REAL(DP) am,e,s,v1,v2,x,y
      REAL(DP), external :: ran1
      !
      IF ( ialpha < 1 ) CALL errore('gamma_dist',  'bad alpha in gamma_dist', 1)
      !
      ! For small  alpha, it is more efficient to calculate Gamma as the waiting time
      ! to the alpha-th event oin a Poisson random process of unit mean.
      ! Define alpha as small for 0 < alpha < 6:
      IF ( ialpha < 6 ) THEN
        !
        x = 1.0D0
        DO j=1,ialpha
          x = x * randy()
        ENDDO
        x = -LOG(x)
      ELSE
        DO
          v1 = 2.0D0*randy()-1.0D0
          v2 = 2.0D0*randy()-1.0D0
          !
          ! need to get this condition met:
          IF ( v1**2+v2**2 > 1.0D0) CYCLE
          !
          y = v2 / v1
          am = ialpha - 1
          s = sqrt(2.0D0 * am + 1.0D0)
          x = s * y + am
          !
          IF ( x <= 0.) CYCLE
          !
          e = (1.0D0+y**2)* exp( am * log( x / am ) - s * y)
          !
          IF (randy() > e) THEN
            CYCLE
          ELSE
            EXIT
          ENDIF
        ENDDO
      ENDIF
    !
    gamma_dist=x
    !
  ENDFUNCTION gamma_dist
  !
  !-----------------------------------------------------------------------
  FUNCTION sum_of_gaussians2(inum_gaussians)
    !-----------------------------------------------------------------------
    !! Returns the sum of inum independent gaussian noises squared, i.e. the result
    !! is equivalent to summing the square of the return values of inum calls
    !! to \(\texttt{gauss_dist}\).
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: inum_gaussians
    !
    REAL(DP) sum_of_gaussians2
    !
    IF ( inum_gaussians < 0 ) THEN
      CALL errore('sum_of_gaussians2',  'negative number of gaussians', 1)
    ELSEIF ( inum_gaussians == 0 ) THEN
      sum_of_gaussians2 = 0.0D0
    ELSEIF ( inum_gaussians == 1 ) THEN
      sum_of_gaussians2 = gauss_dist( 0.0D0, 1.0D0 )**2
    ELSEIF( MODULO(inum_gaussians,2) == 0 ) THEN
      sum_of_gaussians2 = 2.0 * gamma_dist( inum_gaussians/2 )
    ELSE
      sum_of_gaussians2 = 2.0 * gamma_dist((inum_gaussians-1)/2) + &
                gauss_dist( 0.0D0, 1.0D0 )**2
    ENDIF
    !
  ENDFUNCTION sum_of_gaussians2
  !
END MODULE random_numbers
