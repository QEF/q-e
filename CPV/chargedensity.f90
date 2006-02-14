!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  ----------------------------------------------
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------

#include "f_defs.h"

!=----------------------------------------------------------------------=!
      MODULE charge_density
!=----------------------------------------------------------------------=!

!  include modules
        USE kinds
        USE io_files, ONLY: rho_name, rho_name_up, rho_name_down, rho_name_avg
        USE io_files, ONLY: rhounit

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(DP), PARAMETER :: zero = 0.0_DP
        REAL(DP), PARAMETER :: one = 1.0_DP

!  end of module-scope declarations
!  ----------------------------------------------

        PUBLIC :: rhoofr, gradrho

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!


        SUBROUTINE charge_density_closeup()
          INTEGER :: ierr
          RETURN
        END SUBROUTINE charge_density_closeup
!

      REAL(DP) FUNCTION dft_total_charge( ispin, c, cdesc, fi )

!  This subroutine compute the Total Charge in reciprocal space
!  ------------------------------------------------------------

        USE wave_types, ONLY: wave_descriptor

        IMPLICIT NONE

        COMPLEX(DP), INTENT(IN) :: c(:,:)
        INTEGER, INTENT(IN) :: ispin
        TYPE (wave_descriptor), INTENT(IN) :: cdesc
        REAL (DP),  INTENT(IN) :: fi(:)
        INTEGER   :: ib, igs
        REAL(DP) :: rsum
        COMPLEX(DP) :: wdot
        COMPLEX(DP) :: ZDOTC
        EXTERNAL ZDOTC

! ... end of declarations

        IF( ( cdesc%nbl( ispin ) > SIZE( c, 2 ) ) .OR. &
            ( cdesc%nbl( ispin ) > SIZE( fi )     )    ) &
          CALL errore( ' dft_total_charge ', ' wrong sizes ', 1 )

        rsum = 0.0d0

        IF( cdesc%gamma .AND. cdesc%gzero ) THEN

          DO ib = 1, cdesc%nbl( ispin )
            wdot = ZDOTC( ( cdesc%ngwl - 1 ), c(2,ib), 1, c(2,ib), 1 )
            wdot = wdot + DBLE( c(1,ib) )**2 / 2.0d0
            rsum = rsum + fi(ib) * DBLE( wdot )
          END DO

        ELSE

          DO ib = 1, cdesc%nbl( ispin )
            wdot = ZDOTC( cdesc%ngwl, c(1,ib), 1, c(1,ib), 1 )
            rsum = rsum + fi(ib) * DBLE( wdot )
          END DO

        END IF

        dft_total_charge = rsum

        RETURN
      END FUNCTION dft_total_charge



!=----------------------------------------------------------------------=!
!  BEGIN manual

   SUBROUTINE rhoofr (nfi, c0, cdesc, fi, rhoe, box)

!  this routine computes:
!  rhoe = normalized electron density in real space
!
!    rhoe(r) = (sum over ik) weight(ik)
!              (sum over ib) f%s(ib,ik) |psi(r,ib,ik)|^2
!
!    Using quantities in scaled space
!    rhoe(r) = rhoe(s) / Omega
!    rhoe(s) = (sum over ik) weight(ik)
!              (sum over ib) f%s(ib,ik) |psi(s,ib,ik)|^2 
!
!    f%s(ib,ik) = occupation numbers
!    psi(r,ib,ik) = psi(s,ib,ik) / SQRT( Omega ) 
!    psi(s,ib,ik) = INV_FFT (  c0(ik)%w(ig,ib)  )
!
!    ik = index of k point
!    ib = index of band
!    ig = index of G vector
!  ----------------------------------------------
!  END manual

! ... declare modules

    USE fft_base,        ONLY: dfftp, dffts
    USE mp_global,       ONLY: mpime
    USE mp,              ONLY: mp_sum
    USE turbo,           ONLY: tturbo, nturbo, turbo_states, allocate_turbo
    USE cell_module,     ONLY: boxdimensions
    USE wave_types,      ONLY: wave_descriptor
    USE io_global,       ONLY: stdout, ionode
    USE control_flags,   ONLY: force_pairing, iprint
    USE parameters,      ONLY: nspinx
    USE brillouin,       ONLY: kpoints, kp
    USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nnrx
    USE fft_module,      ONLY: invfft



    IMPLICIT NONE

! ... declare subroutine arguments

    INTEGER,              INTENT(IN) :: nfi
    COMPLEX(DP)                     :: c0(:,:,:,:)
    TYPE (boxdimensions), INTENT(IN) :: box
    REAL(DP),          INTENT(IN) :: fi(:,:,:)
    REAL(DP),            INTENT(OUT) :: rhoe(:,:)
    TYPE (wave_descriptor), INTENT(IN) :: cdesc

! ... declare other variables

    INTEGER :: i, is1, is2, j, k, ib, ik, nb, ispin
    INTEGER :: nspin, nbnd, nnr
    REAL(DP)  :: r2, r1, coef3, coef4, omega, rsumg( nspinx ), rsumgs
    REAL(DP)  :: fact, rsumr( nspinx )
    COMPLEX(DP), ALLOCATABLE :: psi2(:)
    INTEGER :: ierr, ispin_wfc
    LOGICAL :: ttprint

! ... end of declarations
!  ----------------------------------------------

    nnr   =  dfftp%nr1x * dfftp%nr2x * dfftp%npl

    omega = box%deth

    nspin = cdesc%nspin

    rsumg = 0.0d0
    rsumr = 0.0d0

    ttprint = .FALSE.
    IF( nfi == 0 .or. mod( nfi, iprint ) == 0 ) ttprint = .TRUE.

    ALLOCATE( psi2( nnrx ), STAT=ierr )
    IF( ierr /= 0 ) CALL errore(' rhoofr ', ' allocating psi2 ', ABS(ierr) )

    IF( tturbo ) THEN
      !
      ! ... if tturbo=.TRUE. some data is stored in memory instead of being
      ! ... recalculated (see card 'TURBO')
      !
      CALL allocate_turbo( nnrx )

    END IF

    rhoe  = zero

    DO ispin = 1, nspin

      ! ... arrange for FFT of wave functions

      ispin_wfc = ispin
      IF( force_pairing ) ispin_wfc = 1

      IF( kp % gamma_only ) THEN

        ! ...  Gamma-point calculation: wave functions are real and can be
        ! ...  Fourier-transformed two at a time as a complex vector

        psi2 = zero

        nbnd = cdesc%nbl( ispin )
        nb   = ( nbnd - MOD( nbnd, 2 ) )

        DO ib = 1, nb / 2

          is1 = 2*ib - 1       ! band index of the first wave function
          is2 = is1  + 1       ! band index of the second wave function

          ! ...  Fourier-transform wave functions to real-scaled space
          ! ...  psi(s,ib,ik) = INV_FFT (  c0(ig,ib,ik)  )

          CALL c2psi( psi2, dffts%nnr, c0( 1, is1, 1, ispin_wfc ), c0( 1, is2, 1, ispin_wfc ), cdesc%ngwl, 2 )
          CALL invfft( 'Wave',psi2, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )

          IF( tturbo .AND. ( ib <= nturbo ) ) THEN
            ! ...  store real-space wave functions to be used in force 
            turbo_states( :, ib ) = psi2( : )
          END IF

          ! ...  occupation numbers divided by cell volume
          ! ...  Remember: rhoe( r ) =  rhoe( s ) / omega

          coef3 = fi( is1, 1, ispin ) / omega
          coef4 = fi( is2, 1, ispin ) / omega 

          ! ...  compute charge density from wave functions

          DO i = 1, nnr

                ! ...  extract wave functions from psi2

                r1 =  DBLE( psi2(i) ) 
                r2 = AIMAG( psi2(i) ) 

                ! ...  add squared moduli to charge density

                rhoe(i,ispin) = rhoe(i,ispin) + coef3 * r1 * r1 + coef4 * r2 * r2

          END DO

        END DO

        IF( MOD( nbnd, 2 ) /= 0 ) THEN

          nb = nbnd

          ! ...  Fourier-transform wave functions to real-scaled space

          CALL c2psi( psi2, dffts%nnr, c0( 1, nb, 1, ispin_wfc ), c0( 1, nb, 1, ispin_wfc ), cdesc%ngwl, 1 )
          CALL invfft( 'Wave', psi2, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )

          ! ...  occupation numbers divided by cell volume

          coef3 = fi( nb, 1, ispin ) / omega

          ! ...  compute charge density from wave functions

          DO i = 1, nnr

             ! ...  extract wave functions from psi2

             r1 = DBLE( psi2(i) )

             ! ...  add squared moduli to charge density

             rhoe(i,ispin) = rhoe(i,ispin) + coef3 * r1 * r1

          END DO

        END IF

      ELSE

        ! ...  calculation with generic k points: wave functions are complex

        psi2 = zero

        DO ik = 1, cdesc%nkl

          DO ib = 1, cdesc%nbl( ispin )

            ! ...  Fourier-transform wave function to real space

            CALL c2psi( psi2, dffts%nnr, c0( 1, nb, 1, ispin_wfc ), c0( 1, nb, 1, ispin_wfc ), cdesc%ngwl, 0 )
            CALL invfft( 'Wave', psi2, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )

            ! ...  occupation numbers divided by cell volume
            ! ...  times the weight of this k point

            coef3 = kp%weight(ik) * fi( ib, ik, ispin ) / omega

            ! ...  compute charge density

            DO i = 1, nnr

                ! ...  add squared modulus to charge density

                rhoe(i,ispin) = rhoe(i,ispin) + coef3 * DBLE( psi2(i) * CONJG(psi2(i)) )

            END DO
          END DO
        END DO

      END IF

      IF( ttprint ) rsumr( ispin ) = SUM( rhoe( :, ispin ) ) * omega / ( nr1 * nr2 * nr3 )

    END DO


    IF( ttprint ) THEN
      !
      DO ispin = 1, nspin
        ispin_wfc = ispin
        IF( force_pairing ) ispin_wfc = 1
        DO ik = 1, cdesc%nkl
          fact = kp%weight(ik)
          IF( cdesc%gamma ) fact = fact * 2.d0
          rsumgs = dft_total_charge( ispin, c0(:,:,ik,ispin_wfc), cdesc, fi(:,ik,ispin) )
          rsumg( ispin ) = rsumg( ispin ) + fact * rsumgs
        END DO
      END DO
      !
      CALL mp_sum( rsumg( 1:nspin ) )
      CALL mp_sum( rsumr( 1:nspin ) )
      !
      if ( nspin == 1 ) then
        WRITE( stdout, 10) rsumg(1), rsumr(1)
      else
        WRITE( stdout, 20) rsumg(1), rsumr(1), rsumg(2), rsumr(2)
      endif

10    FORMAT( /, 3X, 'from rhoofr: total integrated electronic density', &
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 )
20    FORMAT( /, 3X, 'from rhoofr: total integrated electronic density', &
            & /, 3X, 'spin up', & 
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 , &
            & /, 3X, 'spin down', & 
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 )


    END IF

    DEALLOCATE(psi2, STAT=ierr)
    IF( ierr /= 0 ) CALL errore(' rhoofr ', ' deallocating psi2 ', ABS(ierr) )


    RETURN
!=----------------------------------------------------------------------=!
  END SUBROUTINE rhoofr
!=----------------------------------------------------------------------=!


!=----------------------------------------------------------------------=!
   SUBROUTINE gradrho(rhoeg, grho, gx)
!=----------------------------------------------------------------------=!

!     This subroutine calculate the gradient of the charge
!     density in reciprocal space and transforms it to the
!     real space.

     USE fft_base,   ONLY: dfftp
     USE cell_base,  ONLY: tpiba
     USE gvecp,      ONLY: ngm
     USE fft_module, ONLY: invfft

     IMPLICIT NONE

     COMPLEX(DP), INTENT(IN)  :: rhoeg(:)    ! charge density (Reciprocal Space)
     REAL(DP), INTENT(IN)  :: gx(:,:)        ! cartesian components of G-vectors
     REAL(DP),  INTENT(OUT) :: grho(:,:) ! charge density gradient

     INTEGER :: ig, ipol, ierr
     COMPLEX(DP), ALLOCATABLE :: tgrho(:)
     COMPLEX(DP)              :: rg
     COMPLEX(DP), ALLOCATABLE :: psi(:)

     ! ...

     ALLOCATE( tgrho( SIZE( rhoeg ) ), STAT=ierr)
     IF( ierr /= 0 ) CALL errore(' gradrho ', ' allocating tgrho ', ABS(ierr) )
     !
     ALLOCATE( psi( SIZE( grho, 1 ) ), STAT=ierr)
     IF( ierr /= 0 ) CALL errore(' gradrho ', ' allocating psi ', ABS(ierr) )

     DO ipol = 1, 3
       DO ig = 1, SIZE( rhoeg )
         rg        = rhoeg(ig) * gx( ipol, ig )
         tgrho(ig) = tpiba * CMPLX( - AIMAG(rg), DBLE(rg) ) 
       END DO
       CALL rho2psi( 'Dense', psi, dfftp%nnr, tgrho, ngm )
       CALL invfft(  'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
       grho( :, ipol ) = DBLE( psi )
     END DO

     DEALLOCATE(tgrho, STAT=ierr)
     IF( ierr /= 0 ) CALL errore(' gradrho ', ' deallocating tgrho ', ABS(ierr) )

     DEALLOCATE(psi, STAT=ierr)
     IF( ierr /= 0 ) CALL errore(' gradrho ', ' deallocating psi ', ABS(ierr) )

     RETURN
  END SUBROUTINE gradrho


!=----------------------------------------------------------------------=!
      END MODULE charge_density
!=----------------------------------------------------------------------=!


!=----------------------------------------------------------------------=!
!   CP subroutine to compute gradient
!=----------------------------------------------------------------------=!

      subroutine fillgrad(nspin,rhog,gradr)
!     _________________________________________________________________
!
!     calculates gradient of charge density for gradient corrections
!     in: charge density on G-space    out: gradient in R-space
!
      use reciprocal_vectors, only: gx
      use recvecs_indexes, only: np, nm
      use gvecp, only: ng => ngm
      use grid_dimensions, only: nr1, nr2, nr3, &
            nr1x, nr2x, nr3x, nnr => nnrx
      use cell_base, only: tpiba
      USE fft_module, ONLY: invfft
!
      implicit none
! input
      integer nspin
      complex(8) rhog(ng,nspin)
! output
      real(8)    gradr(nnr,3,nspin)
! local
      complex(8), allocatable :: v(:)
      complex(8) ci
      integer iss, ig, ir
!
!
      allocate( v( nnr ) ) 
      !
      ci=(0.0,1.0)
      do iss=1,nspin
         do ig=1,nnr
            v(ig)=(0.0,0.0)
         end do
         do ig=1,ng
            v(np(ig))=      ci*tpiba*gx(1,ig)*rhog(ig,iss)
            v(nm(ig))=CONJG(ci*tpiba*gx(1,ig)*rhog(ig,iss))
         end do
         call invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            gradr(ir,1,iss)=DBLE(v(ir))
         end do
         do ig=1,nnr
            v(ig)=(0.0,0.0)
         end do
         do ig=1,ng
            v(np(ig))= tpiba*(      ci*gx(2,ig)*rhog(ig,iss)-           &
     &                                 gx(3,ig)*rhog(ig,iss) )
            v(nm(ig))= tpiba*(CONJG(ci*gx(2,ig)*rhog(ig,iss)+           &
     &                                 gx(3,ig)*rhog(ig,iss)))
         end do
         call invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnr
            gradr(ir,2,iss)= DBLE(v(ir))
            gradr(ir,3,iss)=AIMAG(v(ir))
         end do
      end do
      !
      deallocate( v )
!
      return
      end subroutine fillgrad

