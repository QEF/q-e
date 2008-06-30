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
!============================================================================
      SUBROUTINE dvdr_tao(dv_dtao, nr1, nr2, nr3, nrx1, nrx2, nrx3 )
! Estimates the derivative of the corrective potential
! at the ionic positions
      
      USE ee_mod,        ONLY : vcomp                       
      USE ions_base,     ONLY : nat, ityp, zv, tau
      USE cell_base,     ONLY : alat, omega, at
      USE kinds,         ONLY : DP
      USE mp,            ONLY : mp_sum
      USE mp_global,     ONLY : intra_pool_comm
      USE fft_base,      ONLY : grid_gather
      !
      IMPLICIT NONE
      !
      REAL(DP)               :: dv_dtao(3, nat)
      INTEGER, INTENT(IN)    :: nr1, nr2, nr3
      INTEGER, INTENT(IN)    :: nrx1, nrx2, nrx3
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
                                df1,                                   &
                                df2,                                   &
                                df3,                                   &
                                f,                                     &
                                g1,                                    &
                                g2,                                    &
                                g3                                     
      !
      INTEGER, EXTERNAL      :: COMPINDEX
      INTEGER, EXTERNAL      :: COMPMOD
      REAL( DP ), EXTERNAL   :: PINTERP
      REAL( DP ), EXTERNAL   :: QINTERP
      REAL( DP ), EXTERNAL   :: DPINTERP
      REAL( DP ), EXTERNAL   :: DQINTERP
      INTEGER, EXTERNAL      :: BOUND
 
!      REAL( DP ), allocatable :: vaux (:)
!
!      allocate ( vaux(nrx1*nrx2*nrx3) )
!#ifdef __PARA
!      vaux(:) = 0.d0
!      call grid_gather(vcomp,vaux)
!      call mp_sum(vaux,intra_pool_comm)
!#else
!      vaux = vcomp
!#endif
      !
      ! ... Initializes the variables
      !
      delta1 = alat * at( 1, 1 ) / DBLE( nr1 )
      delta2 = alat * at( 2, 2 ) / DBLE( nr2 )
      delta3 = alat * at( 3, 3 ) / DBLE( nr3 )
      !
      dv_dtao(:,:) = 0.D0 
      df1 = 0.D0
      df2 = 0.D0
      df3 = 0.D0
      !
      DO na = 1, nat
        !
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
! t1 t2 t3 contains (now) the rest "what was lost" in the discretization process
! and will be recovered by interpolation
!
        !
        ir1 = COMPMOD( ir1, nr1 )
        ir2 = COMPMOD( ir2, nr2 )
        ir3 = COMPMOD( ir3, nr3 )
!
! ir1 ir2 ir3 are the components of the ionic position in grid indexes
!
        bound1 = BOUND( ir1, nr1 )
        bound2 = BOUND( ir2, nr2 )
        bound3 = BOUND( ir3, nr3 )
        !
        f = 0 
        g1 = 0
        g2 = 0
        g3 = 0
        !
        DO a = 0, 1
         DO b = 0, 1
          DO c = 0, 1
           !
           f = vcomp( COMPINDEX( ir1+a,ir2+b,ir3+c,nrx1,nrx2,nrx3 ) )
           g1 = f * DPINTERP( t1, a, bound1 )  &
                  *  PINTERP( t2, b, bound2 )  &
                  *  PINTERP( t3, c, bound3 )
           g2 = f *  PINTERP( t1, a, bound1 )  &
                  * DPINTERP( t2, b, bound2 )  &
                  *  PINTERP( t3, c, bound3 )
           g3 = f *  PINTERP( t1, a, bound1 )  &
                  *  PINTERP( t2, b, bound2 )  &
                  * DPINTERP( t3, c, bound3 )

           df1 = 0.5D0 * (                                             &
             vcomp( COMPINDEX( ir1+a+1,ir2+b,ir3+c,nrx1,nrx2,nrx3 ) )     &
           - vcomp( COMPINDEX( ir1+a-1,ir2+b,ir3+c,nrx1,nrx2,nrx3 ) ) )   
                
           df2 = 0.5D0 * (                                             &
             vcomp( COMPINDEX( ir1+a,ir2+b+1,ir3+c,nrx1,nrx2,nrx3 ) )     &
           - vcomp( COMPINDEX( ir1+a,ir2+b-1,ir3+c,nrx1,nrx2,nrx3 ) ) )   
                
           df3 = 0.5D0 * (                                             &
             vcomp( COMPINDEX( ir1+a,ir2+b,ir3+c+1,nrx1,nrx2,nrx3 ) )     &
           - vcomp( COMPINDEX( ir1+a,ir2+b,ir3+c-1,nrx1,nrx2,nrx3 ) ) )   
                
           dv_dtao(1, na) =  dv_dtao(1, na) + g1 +                     &
                       df1  * DQINTERP( t1, a, bound1 )                &
                       * PINTERP( t2, b, bound2 )                      &
                       * PINTERP( t3, c, bound3 )                       
           dv_dtao(2, na) =  dv_dtao(2, na) + g2 +                     &
                       df2 * PINTERP( t1, a, bound1 )                  &
                       * DQINTERP( t2, b, bound2 )                     &
                       * PINTERP( t3, c, bound3 )                      
           dv_dtao(3, na) =  dv_dtao(3, na) + g3 +                     &
                       df3 * PINTERP( t1, a, bound1 )                  &
                       * PINTERP( t2, b, bound2 )                      &
                       * DQINTERP( t3, c, bound3 )
           !
          END DO
         END DO
        END DO
        dv_dtao(1, na) = dv_dtao(1, na) / delta1
        dv_dtao(2, na) = dv_dtao(2, na) / delta2
        dv_dtao(3, na) = dv_dtao(3, na) / delta3

      END DO 

      END SUBROUTINE dvdr_tao
