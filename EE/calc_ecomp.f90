!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
!--------------------------------------------------------------------
      SUBROUTINE calc_ecomp( rho, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin )
!--------------------------------------------------------------------
      !
      ! ... Calculates the periodic-image-correction
      ! ... (charge-compensation energy).
      ! ...
      ! ... Note: we adopt the convention that electrons are positively charged
      ! ...
      ! 
      USE ions_base,     ONLY : nat, ityp, zv, tau
      USE cell_base,     ONLY : alat, omega, at
      USE ee_mod,        ONLY : vcomp, ecomp, which_compensation
      USE kinds,         ONLY : DP
      USE mp,            ONLY : mp_sum
      USE mp_global,     ONLY : intra_pool_comm
      USE fft_base,      ONLY : grid_gather

      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      INTEGER, INTENT(IN)    :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin 
      REAL(DP), INTENT(IN)   :: rho( nrxx, nspin )
      !
      INTEGER                :: ir,                                    &
                                ir1,                                   &
                                ir2,                                   &
                                ir3,                                   &
                                na,                                    &
                                bound1,                                &
                                bound2,                                &
                                bound3,                                &
                                a,                                     &
                                b,                                     &
                                c
      REAL(DP)               :: delta1,                                &
                                delta2,                                &
                                delta3,                                &
                                t1,                                    &
                                t2,                                    &
                                t3,                                    &
                                f,                                     &
                                g,                                     &
                                df1,                                   &
                                df2,                                   &
                                df3
      !
      INTEGER, EXTERNAL      :: COMPINDEX
      INTEGER, EXTERNAL      :: COMPMOD
      REAL( DP ), EXTERNAL   :: PINTERP
      REAL( DP ), EXTERNAL   :: QINTERP
      INTEGER, EXTERNAL      :: BOUND

      REAL( DP ), allocatable :: vaux (:)

      allocate ( vaux(nrx1*nrx2*nrx3) )
#ifdef __PARA
      vaux(:) = 0.d0
      call grid_gather(vcomp,vaux)
      call mp_sum(vaux, intra_pool_comm)
#else
      vaux = vcomp
#endif
      !
      ! ... Initializes the variables
      !
      delta1 = alat * at( 1, 1 ) / DBLE( nr1 )
      delta2 = alat * at( 2, 2 ) / DBLE( nr2 )
      delta3 = alat * at( 3, 3 ) / DBLE( nr3 )
      !
      ecomp = 0.D0
      !
      ! ... Calculates the energy correction
      !
      !  electronic term
      !
      ecomp = ecomp +                                                  &
             0.5D0 * SUM( vcomp( 1:nrxx ) * rho( 1:nrxx, 1 ) )         &
             * omega / ( nr1 * nr2 * nr3 ) 
      IF ( nspin == 2 ) THEN
        ecomp = ecomp +                                                &
             0.5D0 * SUM( vcomp( 1:nrxx ) * rho( 1:nrxx, 2 ) )         &
             * omega / ( nr1 * nr2 * nr3 )
      END IF
      call mp_sum(ecomp, intra_pool_comm)
      !
      ! ionic term
      !
      DO na = 1, nat
        !
        t1 = tau( 1, na ) * alat / delta1
        t2 = tau( 2, na ) * alat / delta2
        t3 = tau( 3, na ) * alat / delta3
        !
        ir1 = INT( t1 ) + 1
        ir2 = INT( t2 ) + 1
        ir3 = INT( t3 ) + 1
        !
        t1 = t1 - DBLE( ir1 - 1 )
        t2 = t2 - DBLE( ir2 - 1 )
        t3 = t3 - DBLE( ir3 - 1 )
        !
        ir1 = COMPMOD( ir1, nr1 )
        ir2 = COMPMOD( ir2, nr2 )
        ir3 = COMPMOD( ir3, nr3 )
        !
        bound1 = BOUND( ir1, nr1 )
        bound2 = BOUND( ir2, nr2 )
        bound3 = BOUND( ir3, nr3 )
        !
        g = 0.D0
        !
        DO a = 0, 1
         DO b = 0, 1
          DO c = 0, 1
           !
           f = vaux( COMPINDEX( ir1+a,ir2+b,ir3+c,nr1,nr2,nr3 ) )
           g = g + f * PINTERP( t1, a, bound1 )                        &
                     * PINTERP( t2, b, bound2 )                        &
                     * PINTERP( t3, c, bound3 )
           df1 = 0.5D0 * (                                             &
             vaux( COMPINDEX( ir1+a+1,ir2+b,ir3+c,nrx1,nrx2,nrx3 ) )     &
           - vaux( COMPINDEX( ir1+a-1,ir2+b,ir3+c,nrx1,nrx2,nrx3 ) ) )
           df2 = 0.5D0 * (                                             &
             vaux( COMPINDEX( ir1+a,ir2+b+1,ir3+c,nrx1,nrx2,nrx3 ) )     &
           - vaux( COMPINDEX( ir1+a,ir2+b-1,ir3+c,nrx1,nrx2,nrx3 ) ) )
           df3 = 0.5D0 * (                                             &
             vaux( COMPINDEX( ir1+a,ir2+b,ir3+c+1,nrx1,nrx2,nrx3 ) )     &
           - vaux( COMPINDEX( ir1+a,ir2+b,ir3+c-1,nrx1,nrx2,nrx3 ) ) )
           g = g + df1 * QINTERP( t1, a, bound1 )                      &
                       * PINTERP( t2, b, bound2 )                      &
                       * PINTERP( t3, c, bound3 )                      &
                 + df2 * PINTERP( t1, a, bound1 )                      &
                       * QINTERP( t2, b, bound2 )                      &
                       * PINTERP( t3, c, bound3 )                      &
                 + df3 * PINTERP( t1, a, bound1 )                      &
                       * PINTERP( t2, b, bound2 )                      &
                       * QINTERP( t3, c, bound3 )
!
          END DO
         END DO
        END DO
        !
        ecomp = ecomp + 0.5D0 * g * zv( ityp( na ) )
        !
      END DO
      !
      !

      ecomp = - ecomp

      deallocate ( vaux )

      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE calc_ecomp
!--------------------------------------------------------------------
