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
      SUBROUTINE add_dcc_field( vpot, vcorr, rhotot,nelec, &
                                nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                nl, nlm, g, gg, ngm, gstart, nspin  )

!--------------------------------------------------------------------
      !
      ! ... Adds the density-countercharge (DCC) correction vcorr
      ! ... to the periodic-boundary-condition (PBC) potential
      ! ... to obtain the open-boundary-condition (OBC) potential
      !
      USE kinds,         ONLY : DP
      USE io_global,     ONLY : stdout
      USE ions_base,     ONLY : nat,                                   &
                                ityp,                                  &
                                zv
      USE cell_base,     ONLY : alat,                                  &
                                omega,                                 &
                                at
      USE io_global,     ONLY : stdout
      USE ee_mod,        ONLY : mcc => mixing_charge_compensation,     &
                                vloccoul,                              &
                                cellmin,                               &
                                cellmax,                               &
                                mr1, mr2, mr3
      USE mp_global, ONLY: intra_pool_comm, me_pool
      USE mp       , ONLY: mp_sum
      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      INTEGER, INTENT(IN)     :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin, &
                                 ngm, gstart, nl( ngm ), nlm( ngm)

      REAL( DP ), INTENT(IN)  :: g( ngm ), gg( ngm ), nelec

      REAL( DP ), INTENT(IN)  :: rhotot( nrx1*nrx2*nrx3 )
      REAL( DP )              :: vpot( nrx1*nrx2*nrx3 ), &
                                 vcorr( nrx1*nrx2*nrx3 )
      !
      REAL( DP )              :: t1, t2, t3
      REAL( DP )              :: delta1n, delta2n, delta3n
      REAL( DP )              :: delta1m, delta2m, delta3m
      REAL( DP )              :: rhoavg
      REAL( DP )              :: rhoavgcoarse
      REAL( DP)               :: rhocharge, ehart
      !
      INTEGER                 :: mrxx
      !
      CHARACTER( LEN = 256 )  :: which_print
      CHARACTER( LEN = 256 )  :: fileName
      !
      REAL( DP ), ALLOCATABLE :: aux( : )
      REAL( DP ), ALLOCATABLE :: vfft( : )
      REAL( DP ), ALLOCATABLE :: old_vcorr( : )
      REAL( DP ), ALLOCATABLE :: rhocoarse( : )
      REAL( DP ), ALLOCATABLE :: vcorrcoarse( : )
      REAL( DP ), ALLOCATABLE :: eps( : )
      REAL( DP ), ALLOCATABLE :: kap2( : )

      INTEGER :: i

      !
      ! ... Initializes variables
      !
      delta1n = alat * at( 1, 1 ) / DBLE( nr1 )
      delta2n = alat * at( 2, 2 ) / DBLE( nr2 )
      delta3n = alat * at( 3, 3 ) / DBLE( nr3 )
      !
      delta1m = alat * at( 1, 1 ) / DBLE( mr1 )
      delta2m = alat * at( 2, 2 ) / DBLE( mr2 )
      delta3m = alat * at( 3, 3 ) / DBLE( mr3 )

      t1 = alat * at( 1, 1 ) * cellmin( 1 ) / ( cellmax( 1 ) - cellmin( 1 ) )
      t2 = alat * at( 2, 2 ) * cellmin( 2 ) / ( cellmax( 2 ) - cellmin( 2 ) )
      t3 = alat * at( 3, 3 ) * cellmin( 3 ) / ( cellmax( 3 ) - cellmin( 3 ) )
      !
      mrxx = mr1 * mr2 * mr3
      !
      ! ... Saves the old compensating potential
      ! ... taking into account the fact that the
      ! ... average of the local potential is not zero
      !
      ALLOCATE( old_vcorr( nrx1 * nrx2 * nrx3 ) )
      old_vcorr = vcorr
      !
      ! ... the electronic charge density adopts 
      ! ... the convention that electrons are positively charged
      !
      ! ... Interpolates the electronic charge density 
      ! ... on a coarse grid using tricubic interpolation
      !
      ALLOCATE( rhocoarse( mrxx ) )

      CALL multiscale(  rhotot, nr1, nr2, nr3,                         &
                        delta1n, delta2n, delta3n,                     &
                        rhocoarse , mr1, mr2, mr3,                     &
                        delta1m, delta2m, delta3m,                     &
                        - t1, - t2, - t3, .TRUE. )
      !
      ! ... Calculates the average density and the total charge
      !
      rhoavg = SUM( rhotot( : ) ) / DBLE( nr1 * nr2 * nr3 )

      rhoavgcoarse = SUM( rhocoarse( : ) ) / DBLE( mrxx )
      rhocoarse = rhocoarse * rhoavg / rhoavgcoarse
      !
      ! ... Calculates the average charge density and
      ! ... charge to check charge conservation
      !
      rhoavg = ( nelec - SUM( zv( ityp ( 1:nat ) ) ) ) / omega
      !
      ! ... Solves for the potential in Fourier space:
      ! ... nabla^2 vfft = - 4 * pi * e * e * ( rhotot - < rhotot > )
      ! ... with periodic boundary conditions 
      !

      ALLOCATE( aux(nrx1*nrx2*nrx3) )
      aux(:)=vpot

      CALL v_h_from_rho_r( rhotot, nr1, nr2, nr3, nrx1, nrx2, nrx3,  &
                 nrxx, nl,nlm, ngm, gg, gstart, alat, omega, ehart,  &
                 rhocharge, aux )


      ALLOCATE( vfft( mrxx ) )
      CALL multiscale( aux, nr1, nr2, nr3,                            &
                        delta1n, delta2n, delta3n,                     &
                        vfft , mr1, mr2, mr3,                          &
                        delta1m, delta2m, delta3m,                     &
                        - t1, - t2, - t3, .TRUE. )
      !

      DEALLOCATE(aux)

      ALLOCATE( vcorrcoarse( mrxx ) )
      !
      ! ... Defines the source term and the dielectric
      ! ... constant for the Poisson-Boltzmann equation
      !
      ALLOCATE( aux( mrxx ) )
      ALLOCATE( eps( mrxx ) )
      ALLOCATE( kap2( mrxx ) )
      eps( : ) = 1.D0
      kap2( : ) = 0.D0
      aux( : ) = rhoavg
      !
      ! ... Calculates the Dirichlet boundary conditions
      ! ... for the compensating potential
      !
       
      vcorrcoarse = vloccoul - vfft

      CALL add_boundary( rhocoarse, vcorrcoarse, mr1, mr2, mr3,        &
                         delta1m, delta2m, delta3m )
      !
      ! ... Solves the Poisson-Boltzmann
      ! ... equation in real-space with the calculated
      ! ... boundary conditions
      !

      CALL mg_pb_solver( aux, vcorrcoarse, eps, kap2, mr1, mr2, mr3,   &
                         delta1m, delta2m, delta3m )

      !
      DEALLOCATE( aux )
      DEALLOCATE( eps )
      DEALLOCATE( kap2 )
      !
      ! ... Interpolates the corrective potential from the coarse grid
      ! ... to the fine grid using finite-element methods
      !
      CALL multiscale(  vcorrcoarse, mr1, mr2, mr3,                    &
                        delta1m, delta2m, delta3m,                     &
                        vcorr , nr1, nr2, nr3,                         &
                        delta1n, delta2n, delta3n,                     &
                        t1, t2, t3, .FALSE. )

      !
      ! ... Mixes the old and new corrective potentials
      ! ... and adds the resulting correction to vpot
      !
      vcorr = mcc * vcorr + ( 1.D0 - mcc ) * old_vcorr
      !
      ! ... Adds the DCC correction to the PBC potential
      !
      vpot( : ) = vpot( : ) + vcorr( : )
      !
      ! ... Outputs the calculated potentials
      !
      which_print = 'z'
      CALL writetofile( vcorr, 'vcorrz.dat', nr1, nr2, nr3,            &
                        delta1n, delta2n, delta3n, which_print )
      !
      DEALLOCATE( vfft )
      DEALLOCATE( old_vcorr )
      DEALLOCATE( rhocoarse )
      DEALLOCATE( vcorrcoarse )

      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE add_dcc_field
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION compindex( ir1, ir2, ir3, nr1, nr2, nr3 )
!--------------------------------------------------------------------
      ! ... Calculates the composite grid index corresponding
      ! ... to ir1, ir2, ir3
      !
      IMPLICIT NONE
      !
      INTEGER :: compindex
      INTEGER, INTENT(IN) :: ir1, ir2, ir3
      INTEGER, INTENT(IN) :: nr1, nr2, nr3
      INTEGER             :: jr1, jr2, jr3
      !
      jr1 = MODULO( ir1 - 1 , nr1 ) + 1
      jr2 = MODULO( ir2 - 1 , nr2 ) + 1
      jr3 = MODULO( ir3 - 1 , nr3 ) + 1
      !
      compindex = jr1 + ( jr2 -1 ) * nr1 + ( jr3 - 1 ) * nr1 * nr2
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION compindex
!--------------------------------------------------------------------
!--------------------------------------------------------------------
      FUNCTION compmod( ir, nr )
!--------------------------------------------------------------------
      ! ... Calculates the composite grid index corresponding to ir
      !
      IMPLICIT NONE
      !
      INTEGER :: compmod
      INTEGER, INTENT(IN) :: ir
      INTEGER, INTENT(IN) :: nr
      !
      compmod = MODULO( ir - 1, nr ) + 1
      !
      RETURN
      !
!--------------------------------------------------------------------
      END FUNCTION compmod
!--------------------------------------------------------------------
