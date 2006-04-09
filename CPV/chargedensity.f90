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

        PUBLIC :: rhoofr, fillgrad

        INTERFACE rhoofr
           MODULE PROCEDURE rhoofr_fpmd, rhoofr_cp
        END INTERFACE

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

   SUBROUTINE rhoofr_fpmd (nfi, c0, cdesc, fi, rhor, box)

!  this routine computes:
!  rhor = normalized electron density in real space
!
!    rhor(r) = (sum over ik) weight(ik)
!              (sum over ib) f%s(ib,ik) |psi(r,ib,ik)|^2
!
!    Using quantities in scaled space
!    rhor(r) = rhor(s) / Omega
!    rhor(s) = (sum over ik) weight(ik)
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
    USE mp_global,       ONLY: intra_image_comm
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
    REAL(DP),            INTENT(OUT) :: rhor(:,:)
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

    rhor  = zero

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
          ! ...  Remember: rhor( r ) =  rhor( s ) / omega

          coef3 = fi( is1, 1, ispin ) / omega
          coef4 = fi( is2, 1, ispin ) / omega 

          ! ...  compute charge density from wave functions

          DO i = 1, nnr

                ! ...  extract wave functions from psi2

                r1 =  DBLE( psi2(i) ) 
                r2 = AIMAG( psi2(i) ) 

                ! ...  add squared moduli to charge density

                rhor(i,ispin) = rhor(i,ispin) + coef3 * r1 * r1 + coef4 * r2 * r2

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

             rhor(i,ispin) = rhor(i,ispin) + coef3 * r1 * r1

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

                rhor(i,ispin) = rhor(i,ispin) + coef3 * DBLE( psi2(i) * CONJG(psi2(i)) )

            END DO
          END DO
        END DO

      END IF

      IF( ttprint ) rsumr( ispin ) = SUM( rhor( :, ispin ) ) * omega / ( nr1 * nr2 * nr3 )

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
      CALL mp_sum( rsumg( 1:nspin ), intra_image_comm )
      CALL mp_sum( rsumr( 1:nspin ), intra_image_comm )
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
  END SUBROUTINE rhoofr_fpmd
!=----------------------------------------------------------------------=!

!-----------------------------------------------------------------------
   SUBROUTINE rhoofr_cp (nfi,c,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin)
!-----------------------------------------------------------------------
!     the normalized electron density rhor in real space
!     the kinetic energy ekin
!     subroutine uses complex fft so it computes two ft's
!     simultaneously
!
!     rho_i,ij = sum_n < beta_i,i | psi_n >< psi_n | beta_i,j >
!     < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
!                   2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
!
!     e_v = sum_i,ij rho_i,ij d^ion_is,ji
!
      USE kinds,              ONLY: DP
      USE control_flags,      ONLY: iprint, tbuff, iprsta, thdyn, tpre, trhor
      USE ions_base,          ONLY: nat
      USE gvecp,              ONLY: ngm
      USE gvecs,              ONLY: ngs, nps, nms
      USE gvecb,              ONLY: ngb
      USE gvecw,              ONLY: ngw
      USE recvecs_indexes,    ONLY: np, nm
      USE reciprocal_vectors, ONLY: gstart
      USE uspp,               ONLY: nkb
      USE uspp_param,         ONLY: nh, nhm
      USE grid_dimensions,    ONLY: nr1, nr2, nr3, &
                                    nr1x, nr2x, nr3x, nnrx
      USE cell_base,          ONLY: omega
      USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, &
                                        nr1sx, nr2sx, nr3sx, nnrsx
      USE electrons_base,     ONLY: nx => nbspx, n => nbsp, f, ispin, nspin
      USE constants,          ONLY: pi, fpi
      USE mp,                 ONLY: mp_sum
      USE dener,              ONLY: denl, dekin
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_image_comm
      USE funct,              ONLY: dft_is_meta
      USE cg_module,          ONLY: tcg
      USE cp_main_variables,  ONLY: rhopr
      USE fft_module,         ONLY: fwfft, invfft
!
      IMPLICIT NONE
      INTEGER nfi
      REAL(DP) bec(:,:)
      REAL(DP) rhovan(:, :, : )
      REAL(DP) rhor(:,:)
      REAL(DP) rhos(:,:)
      REAL(DP) enl, ekin
      COMPLEX(DP) eigrb( :, : )
      COMPLEX(DP) rhog( :, : )
      COMPLEX(DP) c( :, : )
      INTEGER irb( :, : )

      ! local variables

      INTEGER iss, isup, isdw, iss1, iss2, ios, i, ir, ig
      REAL(DP) rsumr(2), rsumg(2), sa1, sa2
      REAL(DP) rnegsum, rmin, rmax, rsum
      REAL(DP), EXTERNAL :: enkin, ennl
      COMPLEX(DP) ci,fp,fm
      COMPLEX(DP), ALLOCATABLE :: psi(:), psis(:)
      LOGICAL, SAVE :: first = .TRUE.
      
      !

      CALL start_clock( 'rhoofr' )

      ALLOCATE( psi( nnrx ) ) 
      ALLOCATE( psis( nnrsx ) ) 

      ci=(0.0,1.0)

      DO iss=1,nspin
         rhor(:,iss) = 0.d0
         rhos(:,iss) = 0.d0
         rhog(:,iss) = (0.d0, 0.d0)
      END DO
      !
      !  calculation of kinetic energy ekin
      !
      ekin = enkin( c, ngw, f, n )
      !
      IF( tpre ) CALL denkin( c, dekin )
      !
      !     calculation of non-local energy
      !
      enl = ennl( rhovan, bec )
      !
      IF( tpre ) CALL dennl( bec, denl )
      !    
      !    warning! trhor and thdyn are not compatible yet!   
      !
      IF( trhor .AND. ( .NOT. thdyn ) ) THEN
         !
         !   non self-consistent calculation  
         !   charge density is read from unit 47
         !
         IF( first ) THEN
            CALL read_rho( nspin, rhor )
            rhopr = rhor
            first = .FALSE.
         ELSE
            rhor = rhopr
         END IF
!
         IF(nspin.EQ.1)THEN
            iss=1
            DO ir=1,nnrx
               psi(ir)=CMPLX(rhor(ir,iss),0.d0)
            END DO
            CALL fwfft('Dense', psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            DO ig=1,ngm
               rhog(ig,iss)=psi(np(ig))
            END DO
         ELSE
            isup=1
            isdw=2
            DO ir=1,nnrx
               psi(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw))
            END DO
            CALL fwfft('Dense', psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            DO ig=1,ngm
               fp=psi(np(ig))+psi(nm(ig))
               fm=psi(np(ig))-psi(nm(ig))
               rhog(ig,isup)=0.5*CMPLX( DBLE(fp),AIMAG(fm))
               rhog(ig,isdw)=0.5*CMPLX(AIMAG(fp),-DBLE(fm))
            END DO
         ENDIF
!
      ELSE

         !     ==================================================================
         !     self-consistent charge
         !     ==================================================================
         !
         !     important: if n is odd then nx must be .ge.n+1 and c(*,n+1)=0.
         ! 
         IF (MOD(n,2).NE.0) THEN
            DO ig=1,ngw
               c(ig,n+1)=(0.,0.)
            END DO
         ENDIF
         !
         DO i=1,n,2
            psis (:) = (0.d0, 0.d0)
            DO ig=1,ngw
               psis(nms(ig))=CONJG(c(ig,i))+ci*CONJG(c(ig,i+1))
               psis(nps(ig))=c(ig,i)+ci*c(ig,i+1)
               ! write(6,'(I6,4F15.10)') ig, psis(nms(ig)), psis(nps(ig))
            END DO

            CALL invfft('Wave',psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)

            !     wavefunctions in unit 21
            !
            IF(tbuff) WRITE(21,iostat=ios) psis
            iss1=ispin(i)
            sa1=f(i)/omega
            IF (i.NE.n) THEN
               iss2=ispin(i+1)
               sa2=f(i+1)/omega
            ELSE
               iss2=iss1
               sa2=0.0
            END IF
            DO ir=1,nnrsx
               rhos(ir,iss1)=rhos(ir,iss1) + sa1*( DBLE(psis(ir)))**2
               rhos(ir,iss2)=rhos(ir,iss2) + sa2*(AIMAG(psis(ir)))**2
            END DO

            !
            !       buffer 21
            !     
            IF(tbuff) THEN
               IF(ios.NE.0) CALL errore(' rhoofr',' error in writing unit 21',ios)
            ENDIF
            !
         END DO
         !
         IF(tbuff) REWIND 21
         !
         !     smooth charge in g-space is put into rhog(ig)
         !
         IF(nspin.EQ.1)THEN
            iss=1
            DO ir=1,nnrsx
               psis(ir)=CMPLX(rhos(ir,iss),0.d0)
            END DO
            CALL fwfft('Smooth', psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
            DO ig=1,ngs
               rhog(ig,iss)=psis(nps(ig))
            END DO
         ELSE
            isup=1
            isdw=2
             DO ir=1,nnrsx
               psis(ir)=CMPLX(rhos(ir,isup),rhos(ir,isdw))
            END DO
            CALL fwfft('Smooth',psis,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
            DO ig=1,ngs
               fp= psis(nps(ig)) + psis(nms(ig))
               fm= psis(nps(ig)) - psis(nms(ig))
               rhog(ig,isup)=0.5*CMPLX( DBLE(fp),AIMAG(fm))
               rhog(ig,isdw)=0.5*CMPLX(AIMAG(fp),-DBLE(fm))
            END DO
         ENDIF
!
         IF(nspin.EQ.1) THEN
            ! 
            !     case nspin=1
            ! 
            iss=1
            psi (:) = (0.d0, 0.d0)
            DO ig=1,ngs
               psi(nm(ig))=CONJG(rhog(ig,iss))
               psi(np(ig))=      rhog(ig,iss)
            END DO
            CALL invfft('Dense',psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            DO ir=1,nnrx
               rhor(ir,iss)=DBLE(psi(ir))
            END DO
         ELSE 
            !
            !     case nspin=2
            !
            isup=1
            isdw=2
            psi (:) = (0.d0, 0.d0)
            DO ig=1,ngs
               psi(nm(ig))=CONJG(rhog(ig,isup))+ci*CONJG(rhog(ig,isdw))
               psi(np(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
            END DO
            CALL invfft('Dense',psi,nr1,nr2,nr3,nr1x,nr2x,nr3x)
            DO ir=1,nnrx
               rhor(ir,isup)= DBLE(psi(ir))
               rhor(ir,isdw)=AIMAG(psi(ir))
            END DO
         ENDIF
         IF (dft_is_meta()) CALL kedtauofr_meta(c, psi, psis) ! METAGGA
!
         IF(iprsta.GE.3)THEN
            DO iss=1,nspin
               rsumg(iss)=omega*DBLE(rhog(1,iss))
               rsumr(iss)=SUM(rhor(:,iss))*omega/DBLE(nr1*nr2*nr3)
            END DO

            IF ( gstart /= 2 ) THEN
               !
               !    in the parallel case, only one processor has G=0 ! 
               !
               DO iss=1,nspin
                  rsumg(iss)=0.0
               END DO
            END IF
            CALL mp_sum( rsumg( 1:nspin ), intra_image_comm )
            CALL mp_sum( rsumr( 1:nspin ), intra_image_comm )

            IF ( nspin == 1 ) THEN
              WRITE( stdout, 10) rsumg(1), rsumr(1)
            ELSE
              WRITE( stdout, 20) rsumg(1), rsumr(1), rsumg(2), rsumr(2)
            ENDIF

         ENDIF
         !
         !     add vanderbilt contribution to the charge density
         !     drhov called before rhov because input rho must be the smooth part
         !
         IF (tpre) CALL drhov(irb,eigrb,rhovan,rhog,rhor)
         !
         CALL rhov(irb,eigrb,rhovan,rhog,rhor)

         rhopr = rhor

      ENDIF

!     ======================================endif for trhor=============
!
!     here to check the integral of the charge density
!
!
      IF( iprsta .GE. 2 ) THEN
         CALL checkrho(nnrx,nspin,rhor,rmin,rmax,rsum,rnegsum)
         rnegsum=rnegsum*omega/DBLE(nr1*nr2*nr3)
         rsum=rsum*omega/DBLE(nr1*nr2*nr3)
         WRITE( stdout,'(a,4(1x,f12.6))')                                     &
     &     ' rhoofr: rmin rmax rnegsum rsum  ',rmin,rmax,rnegsum,rsum
      END IF
!
      IF( nfi == 0 .OR. MOD(nfi, iprint) == 0 .AND. .NOT. tcg ) THEN

         DO iss=1,nspin
            rsumg(iss)=omega*DBLE(rhog(1,iss))
            rsumr(iss)=SUM(rhor(:,iss),1)*omega/DBLE(nr1*nr2*nr3)
         END DO

         IF (gstart.NE.2) THEN
            ! in the parallel case, only one processor has G=0 ! 
            DO iss=1,nspin
               rsumg(iss)=0.0
            END DO
         END IF

         CALL mp_sum( rsumg( 1:nspin ), intra_image_comm )
         CALL mp_sum( rsumr( 1:nspin ), intra_image_comm )

         IF ( nspin == 1 ) THEN
           WRITE( stdout, 10) rsumg(1), rsumr(1)
         ELSE
           WRITE( stdout, 20) rsumg(1), rsumr(1), rsumg(2), rsumr(2)
         ENDIF

      ENDIF

      DEALLOCATE( psi ) 
      DEALLOCATE( psis ) 

10    FORMAT( /, 3X, 'from rhoofr: total integrated electronic density', &
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 )
20    FORMAT( /, 3X, 'from rhoofr: total integrated electronic density', &
            & /, 3X, 'spin up', &
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 , &
            & /, 3X, 'spin down', &
            & /, 3X, 'in g-space = ', f11.6, 3x, 'in r-space =', f11.6 )
!
      CALL stop_clock( 'rhoofr' )

!
      RETURN
      END SUBROUTINE rhoofr_cp



!=----------------------------------------------------------------------=!

   SUBROUTINE fillgrad( nspin, rhog, gradr )

      !
      !     calculates gradient of charge density for gradient corrections
      !     in: charge density on G-space    out: gradient in R-space
      !
      USE kinds,              ONLY: DP
      use reciprocal_vectors, only: gx
      use recvecs_indexes,    only: np, nm
      use gvecp,              only: ngm
      use grid_dimensions,    only: nr1, nr2, nr3, &
                                    nr1x, nr2x, nr3x, nnrx
      use cell_base,          only: tpiba
      USE fft_module,         ONLY: invfft
!
      implicit none
! input
      integer nspin
      complex(DP) :: rhog( ngm, nspin )
! output
      real(DP) ::    gradr( nnrx, 3, nspin )
! local
      complex(DP), allocatable :: v(:)
      complex(DP) :: ci
      integer     :: iss, ig, ir
!
!
      allocate( v( nnrx ) ) 
      !
      ci = ( 0.0d0, 1.0d0 )
      do iss = 1, nspin
         do ig = 1, nnrx
            v( ig ) = ( 0.0d0, 0.0d0 )
         end do
         do ig=1,ngm
            v(np(ig))=      ci*tpiba*gx(1,ig)*rhog(ig,iss)
            v(nm(ig))=CONJG(ci*tpiba*gx(1,ig)*rhog(ig,iss))
         end do
         call invfft( 'Dense', v, nr1, nr2, nr3, nr1x, nr2x, nr3x )
         do ir=1,nnrx
            gradr(ir,1,iss)=DBLE(v(ir))
         end do
         do ig=1,nnrx
            v(ig)=(0.0,0.0)
         end do
         do ig=1,ngm
            v(np(ig))= tpiba*(      ci*gx(2,ig)*rhog(ig,iss)-           &
     &                                 gx(3,ig)*rhog(ig,iss) )
            v(nm(ig))= tpiba*(CONJG(ci*gx(2,ig)*rhog(ig,iss)+           &
     &                                 gx(3,ig)*rhog(ig,iss)))
         end do
         call invfft('Dense',v,nr1,nr2,nr3,nr1x,nr2x,nr3x)
         do ir=1,nnrx
            gradr(ir,2,iss)= DBLE(v(ir))
            gradr(ir,3,iss)=AIMAG(v(ir))
         end do
      end do
      !
      deallocate( v )
!
      RETURN
    END SUBROUTINE fillgrad


!=----------------------------------------------------------------------=!
      END MODULE charge_density
!=----------------------------------------------------------------------=!
