!
! Copyright (C) 2002-2005 FPMD-CPV groups
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
!  SUBROUTINE newrho(rhor,nfi,tcel,drho)
!  SUBROUTINE invgen(aa,dimaa)
!  ----------------------------------------------
!  END manual

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

! ...   declare module-scope variables
        REAL(DP)  :: achmix
        REAL(DP)  :: g1met2
        REAL(DP)  :: g0chmix2
        INTEGER    :: daamax

        REAL(DP), ALLOCATABLE :: aa_save(:,:)
        COMPLEX(DP), ALLOCATABLE :: rho(:,:)
        COMPLEX(DP), ALLOCATABLE :: rr(:,:)
        COMPLEX(DP), ALLOCATABLE :: chmix(:)
        COMPLEX(DP), ALLOCATABLE :: metric(:)

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
          REAL(DP), INTENT(IN) ::  achmix_inp, g0chmix2_inp
          REAL(DP), INTENT(IN) ::  g1met2_inp
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

      SUBROUTINE newrho(rhor, drho, nfi)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ... declare modules
      USE fft_base, ONLY: dfftp
      USE fft_module, ONLY: fwfft, invfft
      USE ions_base, ONLY: nsp
      USE cell_base, ONLY: tpiba2
      USE reciprocal_vectors, ONLY: gstart, gzero, g
      USE gvecp, ONLY: ngm
      USE wave_base, ONLY: scalw
      USE mp_global, ONLY: intra_image_comm
      USE io_global, ONLY: stdout
      USE mp, ONLY: mp_sum

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(DP), INTENT(INOUT) :: rhor(:)
      REAL(DP), INTENT(OUT) ::  drho
      INTEGER, INTENT(IN) :: nfi

! ... declare other variables
      COMPLEX(DP) :: dr
      COMPLEX(DP) :: rhoout(ngm)
      REAL(DP) :: g02, g12, ar, den, num, rsc
      REAL(DP) :: alpha(daamax)
      REAL(DP), ALLOCATABLE :: aa(:,:)
      REAL(DP), ALLOCATABLE :: rho_old(:)
      INTEGER :: ns, sp, is, ism, i, ig
      LOGICAL, SAVE :: tfirst = .TRUE.
      INTEGER, SAVE :: dimaa, dimaaold, nrho_t, ierr
      COMPLEX(DP), ALLOCATABLE :: psi(:)

! ... end of declarations
!  ----------------------------------------------

      IF( nfi /= 0 .AND. tfirst ) THEN

        CALL errore(' newrho ', ' not initialized ', nfi )

      ELSE IF( nfi == 0 )THEN

        IF( tfirst ) THEN
          CALL allocate_charge_mix( ngm )
        END IF

! ...   define array chmix = A * G^2 / (G^2 + G_0^2) and metric = (G^2 + G_1^2) / G^2
        g02 = g0chmix2 / tpiba2
        g12 = g1met2 / tpiba2
        IF(gzero) THEN
          chmix(1)  = 0.0d0
          metric(1) = 0.0d0
        END IF
        DO ig = gstart, ngm
          chmix(ig)  = achmix * g(ig) / (g(ig)+g02)
          metric(ig) = (g(ig)+g12)    /  g(ig)
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

! ... Fourier tranform of rhor

      ALLOCATE( psi( SIZE( rhor ) ) )

      psi = rhor

      CALL fwfft(   'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
      CALL psi2rho( 'Dense', psi, dfftp%nnr, rhoout, ngm )

      DEALLOCATE( psi )
 

      IF( nrho_t == 1 )THEN

        rho(:,1) = rhoout
        RETURN

      ELSE IF( nrho_t.EQ.2 .OR. (daamax.EQ.1 .AND. nrho_t.GT.1) )THEN

        WRITE( stdout, fmt='( 3X,"charge mixing of order  1")' )

        DO ig = gstart, ngm
          dr = rhoout(ig) - rho(ig,1)
          rr(ig,1) = dr
          rhoout(ig) = rho(ig,1) + chmix(ig) * dr
          rho(ig,is) = rhoout(ig)
        END DO
        IF( gzero ) THEN
          rhoout(1) = rho(1,1)
          rr(1,1)   = (0.d0,0.d0)
        END IF
        IF( daamax /= 1 )THEN
          rsc = scalw(gzero, rr(:,1), rr(:,1), metric)
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

        DO ig = gstart, ngm
          rr(ig,ism) = rhoout(ig) - rho(ig,ism)
        END DO
        IF(gzero) THEN
          rr(1,ism) = (0.d0, 0.d0)
        END IF

! ...   Allocate the new A matrix
        ALLOCATE( aa ( dimaa, dimaa ), STAT=ierr )
        IF( ierr /= 0 ) CALL errore(' newrho ', ' allocating aa ', ierr)

! ...   Fill in new A with the content of the old a
        aa( 1:dimaaold, 1:dimaaold ) = aa_save( 1:dimaaold, 1:dimaaold )

! ...   Compute new matrix A
        DO i = 1, dimaa
          rsc = scalw(gzero,rr(:,i),rr(:,ism),metric)
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

        DO ig = gstart, ngm
          rhoout(ig) = (0.d0,0.d0)
          DO i = 1, dimaa
            rhoout(ig) = rhoout(ig) + alpha(i) * ( rho(ig,i) + chmix(ig) * rr(ig,i) )
          END DO
          rho(ig,is) = rhoout(ig)
        END DO
        IF(gzero) THEN
          rhoout(1) = rho(1,1)
        END IF

      END IF

      ALLOCATE( rho_old( SIZE(rhor) ), STAT=ierr )
      IF( ierr /= 0 ) CALL errore(' newrho ', ' allocating rho_old ', ierr)
      rho_old = rhor 

      ! ... rhor back to real space rhor = FFT( rhoout )
      ! CALL pinvfft(rhor, rhoout)

      ALLOCATE( psi( SIZE( rhor ) ) )

      CALL rho2psi( 'Dense', psi, dfftp%nnr, rhoout, ngm )
      CALL invfft(  'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

      rhor = DBLE( psi )

      drho = SUM( (rho_old - rhor)**2 )

      DEALLOCATE(psi)
      DEALLOCATE(rho_old, STAT=ierr)
      IF( ierr /= 0 ) CALL errore(' newrho ', ' deallocating rho_old ', ierr)

      CALL mp_sum(drho, intra_image_comm)

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
          REAL(DP) :: aa(:,:)

! ...     declare other variables
          REAL(DP) ::  scr1(SIZE(aa,1),SIZE(aa,2))
          REAL(DP) ::  scr2(SIZE(aa,1),SIZE(aa,2))
          REAL(DP) ::  scr3(4*SIZE(aa,1))
          REAL(DP) ::  cond, toleig
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

