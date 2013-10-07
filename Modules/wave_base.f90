!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  BEGIN manual

!==----------------------------------------------==!
        MODULE wave_base
!==----------------------------------------------==!


!  (describe briefly what this module does...)
!  ----------------------------------------------

!  END manual

          USE kinds

          IMPLICIT NONE
          SAVE
          PRIVATE

          REAL(DP) :: frice  = 0.0_DP   !  friction parameter for electronic 
                                        !  damped dynamics
          REAL(DP) :: grease = 0.0_DP   !  friction parameter for electronic 
                                        !  damped dynamics

          PUBLIC :: dotp, hpsi, rande_base, gram_kp_base, gram_gamma_base
          PUBLIC :: converg_base, rande_base_s, scalw

          PUBLIC :: wave_steepest
          PUBLIC :: wave_verlet
          PUBLIC :: wave_speed2

          PUBLIC :: frice, grease

          INTERFACE dotp
            MODULE PROCEDURE dotp_gamma, dotp_kp, dotp_gamma_n, dotp_kp_n
          END INTERFACE

          INTERFACE hpsi
            MODULE PROCEDURE hpsi_gamma, hpsi_kp
          END INTERFACE

          INTERFACE converg_base
            MODULE PROCEDURE converg_base_gamma, converg_base_kp
          END INTERFACE

!==----------------------------------------------==!
        CONTAINS
!==----------------------------------------------==!

      SUBROUTINE gram_kp_base(wf, gid)
        USE mp, ONLY: mp_sum
        COMPLEX(DP) :: wf(:,:)
        INTEGER, INTENT(IN) :: gid
        COMPLEX(DP), PARAMETER :: one  = ( 1.0_DP,0.0_DP)
        COMPLEX(DP), PARAMETER :: onem = (-1.0_DP,0.0_DP)
        COMPLEX(DP), PARAMETER :: zero = ( 0.0_DP,0.0_DP)
        REAL(DP), PARAMETER :: small = 1.e-16_DP
        COMPLEX(DP), ALLOCATABLE :: s(:)
        REAL(DP)    :: anorm
        INTEGER      :: ib, ngw, nb
        ngw = SIZE(wf, 1)
        nb  = SIZE(wf, 2)
        ALLOCATE( s(nb) )
        DO ib = 1, nb
          IF(ib > 1)THEN
             s = zero
             CALL ZGEMV &
               ('C', ngw, ib-1, one, wf(1,1), ngw, wf(1,ib), 1, zero, s(1), 1)
             CALL mp_sum(s,gid)
             CALL ZGEMV &
               ('N', ngw, ib-1, onem, wf(1,1), ngw, s(1), 1, one, wf(1,ib), 1)
          END IF
          anorm = SUM( DBLE( wf(:,ib) * CONJG(wf(:,ib)) ) )
          CALL mp_sum(anorm, gid)
          anorm = 1.0_DP / MAX( SQRT(anorm), small )
          CALL zdscal(ngw, anorm, wf(1,ib), 1)
        END DO
        DEALLOCATE( s )
        RETURN
      END SUBROUTINE gram_kp_base

!==----------------------------------------------==!
!==----------------------------------------------==!
!  BEGIN manual
      SUBROUTINE gram_gamma_base(wf, gzero, gid)

! Gram-Schmidt ortogonalization procedure
! input: cp(2,ngik,n) = ( <g(1 )|psi(1)>..<g(1 )|psi(k)>..<g(1 )|psi(n)> )
!                       ( <g(2 )|psi(1)>..<g(2 )|psi(k)>..<g(2 )|psi(n)> )
!                       ( ...............................................)
!                       ( <g(ng)|psi(1)>..<g(ng)|psi(k)>..<g(ng)|psi(n)> )
! output: the same orthogonalized
!  ----------------------------------------------
! line 7&8   : s(k) = -<psi(k)|g(1)><g(1)|psi(i)>  k=1,..,i-1 (orthonormal)
!                                                  i          (non-orthogonal)
! line   9   : s(k) = 2*sum_g{<psi(k)|g><g|psi(i)>} + s(k)
! line  10   : <g|psi(i)> = <g|psi(i)> - sum_k {s(k) <g|psi(k)>}
! lines 12-15: normalize |psi(i)>
! note: line 2 com. out due to im(<g(1)|psi(k)>)=0 for all k (gam. p. is ass.)
!       s(k) is added in 9 to av. doub. count. of <psi(k)|g(1)><g(1)|psi(i)>
!       |psi(i)> after line 10 is orthogonal to |psi(k)> k=1,...,i-1
!  ----------------------------------------------
!  END manual

        USE mp, ONLY: mp_sum

        COMPLEX(DP), INTENT(INOUT) :: wf(:,:)
        INTEGER, INTENT(IN) :: gid
        LOGICAL, INTENT(IN) :: gzero

        REAL(DP), PARAMETER :: one  =  1.0_DP
        REAL(DP), PARAMETER :: two  =  2.0_DP
        REAL(DP), PARAMETER :: onem = -1.0_DP
        REAL(DP), PARAMETER :: zero =  0.0_DP
        REAL(DP), PARAMETER :: small = 1.e-16_DP
        REAL(DP)  :: dnrm2
        REAL(DP), ALLOCATABLE  :: s(:)
        REAL(DP)  :: anorm, wftmp
        INTEGER    :: ib, nwfr, ngw, nb

        ngw  = SIZE(wf, 1)
        nb   = SIZE(wf, 2)
        nwfr = SIZE(wf, 1) * 2
        ALLOCATE( s(nb) )
        DO ib = 1, nb
          IF(ib.GT.1)THEN
             s = zero
! ...        only the processor that own G=0 
             IF(gzero) THEN
               wftmp = -DBLE(wf(1,ib))
               CALL daxpy(ib-1, wftmp, wf(1,1), nwfr, s(1), 1)
             END IF

             CALL DGEMV('T', nwfr, ib-1, two, wf(1,1), nwfr, wf(1,ib), 1, one, s(1), 1)
             CALL mp_sum(s, gid)
             !WRITE( stdout, fmt = '(I3, 16F8.2)' ) mpime, s(1:nb)
             CALL DGEMV('N', nwfr, ib-1, onem, wf(1,1), nwfr, s(1), 1, one, wf(1,ib), 1)
          END IF
          IF(gzero) THEN
            anorm = dnrm2( 2*(ngw-1), wf(2,ib), 1)
            anorm = 2.0_DP * anorm**2 + DBLE( wf(1,ib) * CONJG(wf(1,ib)) )
          ELSE
            anorm = dnrm2( 2*ngw, wf(1,ib), 1)
            anorm = 2.0_DP * anorm**2
          END IF
          CALL mp_sum(anorm, gid)
          anorm = 1.0_DP / MAX( small, SQRT(anorm) )
          CALL dscal( 2*ngw, anorm, wf(1,ib), 1)
        END DO
        DEALLOCATE( s )

        RETURN
      END SUBROUTINE gram_gamma_base


!==----------------------------------------------==!
!==----------------------------------------------==!

      FUNCTION hpsi_kp( c, dc )

!  (describe briefly what this routine does...)
!  ----------------------------------------------

      IMPLICIT NONE

      COMPLEX(DP) :: zdotc

      COMPLEX(DP) :: c(:,:)
      COMPLEX(DP) :: dc(:)

      COMPLEX(DP), DIMENSION( SIZE( c, 2 ) ) :: hpsi_kp

      INTEGER :: jb, ngw, nx

! ... end of declarations
!  ----------------------------------------------

      IF( SIZE( c, 1 ) /= SIZE( dc ) ) &
        CALL errore(' hpsi_kp ', ' wrong sizes ', 1 )

      ngw = SIZE( c, 1 )
      nx  = SIZE( c, 2 )

      DO jb = 1, nx
        hpsi_kp( jb ) = - zdotc( ngw, c(1,jb), 1, dc(1), 1)
      END DO

      RETURN
      END FUNCTION hpsi_kp

!==----------------------------------------------==!
!==----------------------------------------------==!

      FUNCTION hpsi_gamma( gzero, c, ngw, dc, n, noff )

      IMPLICIT NONE

      COMPLEX(DP) :: c(:,:)
      COMPLEX(DP) :: dc(:)
      LOGICAL, INTENT(IN) :: gzero
      INTEGER, INTENT(IN) :: n, noff, ngw

      REAL(DP), DIMENSION( n ) :: hpsi_gamma

      COMPLEX(DP) :: zdotc

      INTEGER :: j

      IF(gzero) THEN
        DO j = 1, n
          hpsi_gamma(j) = &
            - DBLE( (2.0_DP * zdotc(ngw-1, c(2,j+noff-1), 1, dc(2), 1) + c(1,j+noff-1)*dc(1)) )
        END DO
      ELSE
        DO j = 1, n
          hpsi_gamma(j) = - DBLE( (2.0_DP * zdotc(ngw, c(1,j+noff-1), 1, dc(1), 1)) )
        END DO
      END IF
      RETURN
      END FUNCTION hpsi_gamma

!==----------------------------------------------==!
!==----------------------------------------------==!


!  BEGIN manual

      SUBROUTINE converg_base_gamma(gzero, cgrad, gemax, cnorm, comm)

!  this routine checks for convergence, by computing the norm of the
!  gradients of wavefunctions
!  version for the Gamma point
!  ----------------------------------------------
!  END manual

        USE mp, ONLY: mp_sum, mp_max

        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP) :: cgrad(:,:,:)
        LOGICAL, INTENT(IN) :: gzero
        INTEGER, INTENT(IN) :: comm
        REAL(DP), INTENT(OUT) :: gemax, cnorm

! ...   declare other variables
        INTEGER    :: imx, IZAMAX, i, nb, ngw
        REAL(DP)  :: gemax_l

! ...   end of declarations
!  ----------------------------------------------

        ngw     = SIZE( cgrad, 1)
        nb      = SIZE( cgrad, 2)

        gemax_l = 0.0_DP
        cnorm   = 0.0_DP

        DO i = 1, nb
          imx = IZAMAX( ngw, cgrad(1, i, 1), 1 )
          IF ( gemax_l < ABS( cgrad(imx, i, 1) ) ) THEN
            gemax_l = ABS ( cgrad(imx, i, 1) )
          END IF
          cnorm = cnorm + dotp(gzero, cgrad(:,i,1), cgrad(:,i,1), comm)
        END DO

        CALL mp_max(gemax_l, comm)
        CALL mp_sum(nb, comm)
        CALL mp_sum(ngw, comm)

        gemax = gemax_l
        cnorm = SQRT( cnorm / (nb * ngw) )

        RETURN
      END SUBROUTINE converg_base_gamma

!  ----------------------------------------------
!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE converg_base_kp(weight, cgrad, gemax, cnorm, comm)


!  this routine checks for convergence, by computing the norm of the
!  gradients of wavefunctions
!  version for generic k-points
!  ----------------------------------------------
!  END manual

        USE mp, ONLY: mp_sum, mp_max

        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP) :: cgrad(:,:,:)
        REAL(DP), INTENT(IN)  :: weight(:)
        REAL(DP), INTENT(OUT) :: gemax, cnorm
        INTEGER, INTENT(IN) :: comm

! ...   declare other variables
        INTEGER    :: nb, ngw, nk, iabs, IZAMAX, i, ik
        REAL(DP)  :: gemax_l, cnormk
        COMPLEX(DP) :: zdotc

! ...   end of declarations
!  ----------------------------------------------

        ngw = SIZE( cgrad, 1)
        nb  = SIZE( cgrad, 2)
        nk  = SIZE( cgrad, 3)
 
        gemax_l = 0.0_DP
        cnorm   = 0.0_DP
 
        DO ik = 1, nk
          cnormk  = 0.0_DP
          DO i = 1, nb
            iabs = IZAMAX( ngw, cgrad(1,i,ik), 1)
            IF( gemax_l < ABS( cgrad(iabs,i,ik) ) ) THEN
              gemax_l = ABS( cgrad(iabs,i,ik) )
            END IF
            cnormk = cnormk + DBLE( zdotc(ngw, cgrad(1,i,ik), 1, cgrad(1,i,ik), 1))
          END DO
          cnormk = cnormk * weight(ik)
          cnorm = cnorm + cnormk
        END DO

        CALL mp_max(gemax_l, comm)
        CALL mp_sum(cnorm, comm)
        CALL mp_sum(nb, comm)
        CALL mp_sum(ngw, comm)

        gemax = gemax_l
        cnorm = SQRT( cnorm / ( nb * ngw ) )

        RETURN
      END SUBROUTINE converg_base_kp



!==----------------------------------------------==!
!==----------------------------------------------==!

          REAL(DP) FUNCTION wdot_gamma(gzero, ng, a, b)

            LOGICAL, INTENT(IN) :: gzero
            COMPLEX(DP) :: a(:), b(:)
            INTEGER, OPTIONAL, INTENT(IN) :: ng

            REAL(DP) :: ddot
            INTEGER :: n

            n = MIN( SIZE(a), SIZE(b) )
            IF ( PRESENT (ng) ) n = MIN( n, ng )

            IF ( n < 1 ) &
              CALL errore( ' wdot_gamma ', ' wrong dimension ', 1 )

            IF (gzero) THEN
              wdot_gamma = ddot( 2*(n-1), a(2), 1, b(2), 1)
              wdot_gamma = 2.0_DP * wdot_gamma + DBLE( a(1) ) * DBLE( b(1) ) 
            ELSE
              wdot_gamma = 2.0_DP * ddot( 2*n, a(1), 1, b(1), 1)
            END IF 

            RETURN
          END FUNCTION wdot_gamma

!==----------------------------------------------==!
!==----------------------------------------------==!

          REAL(DP) FUNCTION dotp_gamma(gzero, ng, a, b, comm)

! ... Compute the dot product between distributed complex vectors "a" and "b"
! ... representing HALF-SPACE complex wave functions, with the G-point symmetry
! ... a( -G ) = CONJG( a( G ) ). Only half of the values plus G=0 are really
! ... stored in the array.
!
! ... dotp = < a | b >
!

            USE mp, ONLY: mp_sum

            REAL(DP) :: ddot
            REAL(DP) :: dot_tmp
            INTEGER, INTENT(IN) :: ng
            LOGICAL, INTENT(IN) :: gzero
            INTEGER, INTENT(IN) :: comm

            COMPLEX(DP) :: a(:), b(:)
            INTEGER :: n

            n = MIN( SIZE(a), SIZE(b) )
            n = MIN( n, ng )

            IF ( n < 1 ) &
              CALL errore( ' dotp_gamma ', ' wrong dimension ', 1 )

! ...       gzero is true on the processor where the first element of the
! ...       input arrays is the coefficient of the G=0 plane wave
!
            IF (gzero) THEN
              dot_tmp = ddot( 2*(n-1), a(2), 1, b(2), 1)
              dot_tmp = 2.0_DP * dot_tmp + DBLE( a(1) ) * DBLE( b(1) ) 
            ELSE
              dot_tmp = ddot( 2*n, a(1), 1, b(1), 1)
              dot_tmp = 2.0_DP*dot_tmp
            END IF 

            CALL mp_sum( dot_tmp, comm )
            dotp_gamma = dot_tmp

            RETURN
          END FUNCTION dotp_gamma

!==----------------------------------------------==!
!==----------------------------------------------==!

          REAL(DP) FUNCTION dotp_gamma_n(gzero, a, b, comm)

! ...  Compute the dot product between distributed complex vectors "a" and "b"
! ...  representing HALF-SPACE complex wave functions, with the G-point symmetry
! ...  a( -G ) = CONJG( a( G ) ). Only half of the values plus G=0 are really
! ...  stored in the array.

            USE mp, ONLY: mp_sum

            LOGICAL, INTENT(IN) :: gzero
            INTEGER, INTENT(IN) :: comm

            COMPLEX(DP) :: a(:), b(:)
            INTEGER :: n

            n = MIN( SIZE(a), SIZE(b) )

            IF ( n < 1 ) &
              CALL errore( ' dotp_gamma_n ', ' wrong dimension ', 1 )

            dotp_gamma_n = dotp_gamma(gzero, n, a, b, comm)

            RETURN
          END FUNCTION 


!==----------------------------------------------==!
!==----------------------------------------------==!

          COMPLEX(DP) FUNCTION dotp_kp(ng, a, b, comm)

! ...  Compute the dot product between distributed complex vectors "a" and "b"
! ...  representing FULL-SPACE complex wave functions 

            USE mp, ONLY: mp_sum

            COMPLEX(DP) :: zdotc
            INTEGER, INTENT(IN) :: ng
            COMPLEX(DP) :: a(:),b(:)
            INTEGER, INTENT(IN) :: comm

            COMPLEX(DP) :: dot_tmp
            INTEGER      :: n

            n = MIN( SIZE(a), SIZE(b) )
            n = MIN( n, ng )

            IF ( n < 1 ) &
              CALL errore( ' dotp_kp ', ' wrong dimension ', 1 )

            dot_tmp = zdotc(ng, a(1), 1, b(1), 1)

            CALL mp_sum(dot_tmp, comm)
            dotp_kp = dot_tmp

            RETURN
          END FUNCTION dotp_kp

!==----------------------------------------------==!
!==----------------------------------------------==!

          COMPLEX(DP) FUNCTION dotp_kp_n(a, b, comm)

! ...  Compute the dot product between distributed complex vectors "a" and "b"
! ...  representing FULL-SPACE complex wave functions 

            USE mp, ONLY: mp_sum

            COMPLEX(DP) zdotc
            COMPLEX(DP), INTENT(IN) :: a(:),b(:)
            INTEGER, INTENT(IN) :: comm

            COMPLEX(DP) :: dot_tmp
            INTEGER :: n

            n = MIN( SIZE(a), SIZE(b) )

            IF ( n < 1 ) &
              CALL errore( ' dotp_kp_n ', ' wrong dimension ', 1 )

            dot_tmp = zdotc( n, a(1), 1, b(1), 1)

            CALL mp_sum( dot_tmp, comm )
            dotp_kp_n = dot_tmp

            RETURN
          END FUNCTION dotp_kp_n

!==----------------------------------------------==!
!==----------------------------------------------==!

          COMPLEX(DP) FUNCTION wdot_kp(ng, a, b)

! ...  Compute the dot product between complex vectors "a" and "b"
! ...  representing FULL-SPACE complex wave functions 
! ...  Note this is a _SCALAR_ subroutine

            COMPLEX(DP) :: a(:), b(:)
            INTEGER, INTENT(IN), OPTIONAL :: ng

            COMPLEX(DP) :: zdotc
            INTEGER :: n

            n = MIN( SIZE(a), SIZE(b) )
            IF ( PRESENT (ng) ) n = MIN( n, ng )

            IF ( n < 1 ) &
              CALL errore( ' dotp_kp_n ', ' wrong dimension ', 1 )

            wdot_kp = zdotc(n, a(1), 1, b(1), 1)

            RETURN
          END FUNCTION wdot_kp

!==----------------------------------------------==!
!==----------------------------------------------==!

      SUBROUTINE rande_base(wf,ampre)

!  randomize wave functions coefficients
!  ----------------------------------------------
      USE random_numbers, ONLY : randy
      IMPLICIT NONE
! ... declare subroutine arguments
      COMPLEX(DP)          :: wf(:,:)
      REAL(DP), INTENT(IN) :: ampre

! ... declare other variables
      INTEGER i, j
      REAL(DP)  rranf1, rranf2
! ... end of declarations
!  ----------------------------------------------
      DO i = 1, SIZE(wf, 2)
        DO j = 1, SIZE( wf, 1)
          rranf1 = 0.5_DP - randy()
          rranf2 = 0.5_DP - randy()
          wf(j,i) = wf(j,i) + ampre * CMPLX(rranf1, rranf2, KIND=DP)
        END DO
      END DO
      RETURN
      END SUBROUTINE rande_base

!==----------------------------------------------==!

      SUBROUTINE rande_base_s(wf,ampre)

!  randomize wave functions coefficients
!  ----------------------------------------------
      USE random_numbers, ONLY : randy
      IMPLICIT NONE
! ... declare subroutine arguments
      COMPLEX(DP)          :: wf(:)
      REAL(DP), INTENT(IN) :: ampre
! ... declare other variables
      INTEGER j
      REAL(DP)  rranf1, rranf2
! ... end of declarations
!  ----------------------------------------------
      DO j = 1, SIZE( wf )
        rranf1 = 0.5_DP - randy()
        rranf2 = 0.5_DP - randy()
        wf(j) = wf(j) + ampre * CMPLX(rranf1, rranf2, KIND=DP)
      END DO
      RETURN
      END SUBROUTINE rande_base_s

!==----------------------------------------------==!
!==----------------------------------------------==!


       REAL(DP) FUNCTION scalw(gzero, RR1, RR2, metric, comm)

         USE mp, ONLY: mp_sum

         IMPLICIT NONE

         COMPLEX(DP), INTENT(IN) :: rr1(:), rr2(:), metric(:)
         LOGICAL, INTENT(IN) :: gzero
         INTEGER, INTENT(IN) :: comm
         INTEGER :: ig, gstart, ngw
         REAL(DP) :: rsc

         ngw = MIN( SIZE(rr1), SIZE(rr2), SIZE(metric) )
         rsc = 0.0_DP

         gstart = 1
         IF (gzero) gstart = 2

         DO ig = gstart, ngw
           rsc = rsc + rr1( ig ) * CONJG( rr2( ig ) ) * metric( ig )
         END DO

         CALL mp_sum(rsc, comm)

         scalw = rsc

         RETURN
       END FUNCTION scalw

!==----------------------------------------------==!
!==----------------------------------------------==!

   SUBROUTINE wave_steepest( CP, C0, dt2m, grad, ngw, idx )
      IMPLICIT NONE
      COMPLEX(DP), INTENT(OUT) :: CP(:)
      COMPLEX(DP), INTENT(IN) :: C0(:)
      COMPLEX(DP), INTENT(IN) :: grad(:)
      REAL(DP), INTENT(IN) ::  dt2m(:)
      INTEGER, OPTIONAL, INTENT(IN) :: ngw, idx
      !
      IF( PRESENT( ngw ) .AND. PRESENT( idx ) ) THEN
         CP( : )  = C0( : )  + dt2m(:) * grad( (idx-1)*ngw+1 : idx*ngw )
      ELSE
         CP( : )  = C0( : )  + dt2m(:) * grad(:)
      END IF
      !
      RETURN
   END SUBROUTINE wave_steepest

!==----------------------------------------------==!
!==----------------------------------------------==!

   SUBROUTINE wave_verlet( cm, c0, ver1, ver2, ver3, grad, ngw, idx )
      IMPLICIT NONE
      COMPLEX(DP), INTENT(INOUT) :: cm(:)
      COMPLEX(DP), INTENT(IN) :: c0(:)
      COMPLEX(DP), INTENT(IN) :: grad(:)
      REAL(DP), INTENT(IN) ::  ver1, ver2, ver3(:)
      INTEGER, OPTIONAL, INTENT(IN) :: ngw, idx
      !
      IF( PRESENT( ngw ) .AND. PRESENT( idx ) ) THEN
         cm( : )  = ver1 * c0( : ) + ver2 * cm( : ) + ver3( : ) * grad( (idx-1)*ngw+1:idx*ngw)
      ELSE
         cm( : )  = ver1 * c0( : ) + ver2 * cm( : ) + ver3( : ) * grad( : )
      END IF
      !
      RETURN
   END SUBROUTINE wave_verlet

!==----------------------------------------------==!
!==----------------------------------------------==!

   FUNCTION wave_speed2( cp, cm, wmss, fact )
     IMPLICIT NONE
     COMPLEX(DP), INTENT(IN) :: cp(:)
     COMPLEX(DP), INTENT(IN) :: cm(:)
     REAL(DP) :: wmss(:), fact
     REAL(DP) :: wave_speed2
     REAL(DP) :: ekinc
     COMPLEX(DP) :: speed
     INTEGER :: j
     speed  = ( cp(1) - cm(1) )
     ekinc  = fact * wmss(1) * CONJG( speed ) * speed
     DO j = 2, SIZE( cp )
       speed  = ( cp(j) - cm(j) )
       ekinc  = ekinc + wmss(j) * CONJG( speed ) * speed
     END DO
     wave_speed2 = ekinc
     RETURN
   END FUNCTION wave_speed2


!==----------------------------------------------==!
       END MODULE wave_base
!==----------------------------------------------==!
