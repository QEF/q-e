!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sun Nov 21 13:29:38 MET 1999
!  ----------------------------------------------
!  BEGIN manual

!=----------------------------------------------------------------------------=!
      MODULE charge_mix
!=----------------------------------------------------------------------------=!

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE allocate_charge_mix(ng_l)
!  SUBROUTINE deallocate_charge_mix
!  SUBROUTINE charge_mix_print_info(unit)
!  SUBROUTINE newrho(rhoe,nfi,tcel,drho,gv)
!  SUBROUTINE invgen(aa,dimaa)
!  ----------------------------------------------
!  END manual

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

! ...   declare module-scope variables
        REAL(dbl)  :: achmix
        REAL(dbl)  :: g1met2
        REAL(dbl)  :: g0chmix2
        INTEGER    :: daamax

        REAL(dbl), ALLOCATABLE :: aa_save(:,:)
        COMPLEX(dbl), ALLOCATABLE :: rho(:,:)
        COMPLEX(dbl), ALLOCATABLE :: rr(:,:)
        COMPLEX(dbl), ALLOCATABLE :: chmix(:)
        COMPLEX(dbl), ALLOCATABLE :: metric(:)

! ...   end of module-scope declarations
!  ----------------------------------------------

        PUBLIC :: charge_mix_setup
        PUBLIC :: allocate_charge_mix, deallocate_charge_mix
        PUBLIC :: charge_mix_print_info
        PUBLIC :: newrho

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE charge_mix_setup(achmix_inp, g0chmix2_inp, daamax_inp, g1met2_inp)
          REAL(dbl), INTENT(IN) ::  achmix_inp, g0chmix2_inp
          REAL(dbl), INTENT(IN) ::  g1met2_inp
          INTEGER, INTENT(IN) :: daamax_inp
          achmix = achmix_inp
          g0chmix2 = g0chmix2_inp
          daamax = daamax_inp
          g1met2 = g1met2_inp
          RETURN
        END SUBROUTINE charge_mix_setup

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE allocate_charge_mix(ng)
          INTEGER, INTENT(IN) :: ng
          INTEGER :: ierr
          ALLOCATE( rho(ng, daamax), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating rho ', ierr)
          ALLOCATE( rr(ng, daamax), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating rr ', ierr)
          ALLOCATE( aa_save(daamax, daamax), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating aa_save ', ierr)
          ALLOCATE( chmix(ng), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating chmix ', ierr)
          ALLOCATE( metric(ng), STAT=ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_charge_mix ', ' allocating metric ', ierr)
          RETURN
        END SUBROUTINE allocate_charge_mix

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE deallocate_charge_mix
          INTEGER :: ierr
          IF( ALLOCATED(rho) ) THEN
            DEALLOCATE(rho, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating rho ', ierr)
          END IF
          IF( ALLOCATED(rr) ) THEN
            DEALLOCATE(rr, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating rr ', ierr)
          END IF
          IF( ALLOCATED(aa_save) ) THEN
            DEALLOCATE(aa_save, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating aa_save ', ierr)
          END IF
          IF( ALLOCATED(chmix) ) THEN
            DEALLOCATE(chmix, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating chmix ', ierr)
          END IF
          IF( ALLOCATED(metric) ) THEN
            DEALLOCATE(metric, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_charge_mix ', ' deallocating metric ', ierr)
          END IF
          RETURN
        END SUBROUTINE deallocate_charge_mix

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE charge_mix_print_info(unit)
          INTEGER, INTENT(IN) :: unit
            WRITE(unit,300)
            WRITE(unit,310) achmix, g0chmix2, g1met2
            WRITE(unit,320) daamax
300         FORMAT(/,3X,'Charge mixing parameters:')
310         FORMAT(  3X,'A = ', D14.6, ' G0^2 = ', D14.6,' G1^2 = ',D14.6)
320         FORMAT(  3X,'charge mixing matrix maximum size = ', I5)
          RETURN
        END SUBROUTINE charge_mix_print_info

!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE newrho(rhoe, drho, gv, nfi)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ... declare modules
      USE fft, ONLY: pfwfft, pinvfft
      USE ions_base, ONLY: nsp
      USE cell_base, ONLY: tpiba2
      USE cp_types, ONLY: recvecs
      USE wave_base, ONLY: scalw
      USE mp_global, ONLY: group
      USE io_global, ONLY: stdout
      USE mp, ONLY: mp_sum

      IMPLICIT NONE

! ... declare subroutine arguments
      TYPE (recvecs), INTENT(IN) :: gv
      REAL(dbl), INTENT(INOUT) :: rhoe(:,:,:)
      REAL(dbl), INTENT(OUT) ::  drho
      INTEGER, INTENT(IN) :: nfi

! ... declare other variables
      COMPLEX(dbl) :: dr
      COMPLEX(dbl) :: rhoout(gv%ng_l)
      REAL(dbl) :: g02, g12, ar, den, num, rsc
      REAL(dbl) :: alpha(daamax)
      REAL(dbl), ALLOCATABLE :: aa(:,:)
      REAL(dbl), ALLOCATABLE :: rho_old(:,:,:)
      INTEGER :: ns, sp, is, ism, i, ig, gstart
      LOGICAL, SAVE :: tfirst = .TRUE.
      INTEGER, SAVE :: dimaa, dimaaold, nrho_t, ierr

! ... end of declarations
!  ----------------------------------------------

      gstart = 1
      IF( gv%gzero ) gstart = 2

      IF( nfi /= 0 .AND. tfirst ) THEN

        CALL errore(' newrho ', ' not initialized ', nfi )

      ELSE IF( nfi == 0 )THEN

        IF( tfirst ) THEN
          CALL allocate_charge_mix( gv%ng_l )
        END IF

! ...   define array chmix = A * G^2 / (G^2 + G_0^2) and metric = (G^2 + G_1^2) / G^2
        g02 = g0chmix2 / tpiba2
        g12 = g1met2 / tpiba2
        IF(gv%gzero) THEN
          chmix(1)  = 0.0d0
          metric(1) = 0.0d0
        END IF
        DO ig = gstart, gv%ng_l
          chmix(ig)  = achmix * gv%hg_l(ig) / (gv%hg_l(ig)+g02)
          metric(ig) = (gv%hg_l(ig)+g12)    /  gv%hg_l(ig)
        END DO
        tfirst = .FALSE.

      END IF

! ... Reset matrix dimension for the first iteration / initialization
      IF( nfi <= 1 )THEN
        dimaa  = 0
        nrho_t = 0
      END IF

! ... Now update matrix dimension and counter
      nrho_t = nrho_t + 1

      dimaaold = dimaa                    ! save the previous matrix dimension 
      dimaa    = MIN( daamax, nrho_t-1 )  ! number of densities and rr saved up to now

      ism      = MOD( nrho_t-1, daamax )
      if( ism == 0 ) ism = daamax
      is       = MOD( nrho_t  , daamax )
      if( is  == 0 ) is  = daamax

! ... Fourier tranform of rhoe
      CALL pfwfft(rhoout,rhoe)

      IF( nrho_t == 1 )THEN

        rho(:,1) = rhoout
        RETURN

      ELSE IF( nrho_t.EQ.2 .OR. (daamax.EQ.1 .AND. nrho_t.GT.1) )THEN

        WRITE( stdout, fmt='( 3X,"charge mixing of order  1")' )

        DO ig = gstart, gv%ng_l
          dr = rhoout(ig) - rho(ig,1)
          rr(ig,1) = dr
          rhoout(ig) = rho(ig,1) + chmix(ig) * dr
          rho(ig,is) = rhoout(ig)
        END DO
        IF( gv%gzero ) THEN
          rhoout(1) = rho(1,1)
          rr(1,1)   = (0.d0,0.d0)
        END IF
        IF( daamax /= 1 )THEN
          rsc = scalw(gv%gzero, rr(:,1), rr(:,1), metric)
          aa_save(1, 1) =  rsc
        END IF

      ELSE

        IF( dimaa < 1 .OR. dimaa > daamax ) THEN
          CALL errore(' newrho ', ' dimaa out of range ', dimaa )
        END IF
        IF( dimaaold < 1 .OR. dimaaold > daamax ) THEN
          CALL errore(' newrho ', ' dimaaold out of range ', dimaaold )
        END IF

        WRITE( stdout, fmt='( 3X,"charge mixing of order ",I2)' ) dimaa

        DO ig = gstart, gv%ng_l
          rr(ig,ism) = rhoout(ig) - rho(ig,ism)
        END DO
        IF(gv%gzero) THEN
          rr(1,ism) = (0.d0, 0.d0)
        END IF

! ...   Allocate the new A matrix
        ALLOCATE( aa ( dimaa, dimaa ), STAT=ierr )
        IF( ierr /= 0 ) CALL errore(' newrho ', ' allocating aa ', ierr)

! ...   Fill in new A with the content of the old a
        aa( 1:dimaaold, 1:dimaaold ) = aa_save( 1:dimaaold, 1:dimaaold )

! ...   Compute new matrix A
        DO i = 1, dimaa
          rsc = scalw(gv%gzero,rr(:,i),rr(:,ism),metric)
          aa(i,ism)=  rsc
          aa(ism,i)=  rsc
        END DO

! ...   Save the content of A for the next iteration
        aa_save( 1:dimaa, 1:dimaa ) = aa( 1:dimaa, 1:dimaa )

! ...   Compute alphas
        CALL invgen( aa )
        den = SUM( aa )
        DO i = 1, dimaa
          alpha(i) = SUM( aa(:,i) ) / den
        END DO

        DEALLOCATE( aa, STAT=ierr )
        IF( ierr /= 0 ) CALL errore(' newrho ', ' deallocating aa ', ierr)

        DO ig = gstart, gv%ng_l
          rhoout(ig) = (0.d0,0.d0)
          DO i = 1, dimaa
            rhoout(ig) = rhoout(ig) + alpha(i) * ( rho(ig,i) + chmix(ig) * rr(ig,i) )
          END DO
          rho(ig,is) = rhoout(ig)
        END DO
        IF(gv%gzero) THEN
          rhoout(1) = rho(1,1)
        END IF

      END IF

      ALLOCATE( rho_old( SIZE(rhoe, 1), SIZE(rhoe, 2), SIZE(rhoe, 3) ), STAT=ierr )
      IF( ierr /= 0 ) CALL errore(' newrho ', ' allocating rho_old ', ierr)
      rho_old = rhoe 

! ... rhoe back to real space rhoe = FFT( rhoout )
      CALL pinvfft(rhoe, rhoout)
      drho = SUM( (rho_old - rhoe)**2 )

      DEALLOCATE(rho_old, STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' newrho ', ' deallocating rho_old ', ierr)
      CALL mp_sum(drho, group)

      RETURN
      END SUBROUTINE newrho

!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE invgen( aa )

!  (describe briefly what this routine does...)
!  ----------------------------------------------

          IMPLICIT NONE

! ...     declare subroutine arguments
          INTEGER dimaa
          REAL(dbl) :: aa(:,:)

! ...     declare other variables
          REAL(dbl) ::  scr1(SIZE(aa,1),SIZE(aa,2))
          REAL(dbl) ::  scr2(SIZE(aa,1),SIZE(aa,2))
          REAL(dbl) ::  scr3(4*SIZE(aa,1))
          REAL(dbl) ::  cond, toleig
          INTEGER   ::  info, iopt, mrank

! ... end of declarations
!  ----------------------------------------------

          toleig = 1.d-10
          iopt   = 10
          CALL geninv(aa, SIZE(aa,1), SIZE(aa,2), mrank, cond, scr1, scr2, scr3, toleig, info, iopt)
          RETURN
        END SUBROUTINE invgen

!  ----------------------------------------------
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
      END MODULE charge_mix
!=----------------------------------------------------------------------------=!

