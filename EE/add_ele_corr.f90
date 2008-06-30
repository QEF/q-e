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
      SUBROUTINE add_ele_corr( vltot, rho, nelec,                     &
                               nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                               nl, nlm, g, gg, ngm, gstart, nspin     )
!--------------------------------------------------------------------
      !
      ! ... Adds the appropriate countercharge correction to 
      ! ... the periodic potential
      !
      USE io_global,     ONLY : stdout
      USE kinds,         ONLY : DP
      USE ee_mod,        ONLY : which_compensation, vcomp
      USE cell_base,     ONLY : alat, omega
      USE mp,            ONLY : mp_sum
      USE mp_global,     ONLY : intra_pool_comm, me_pool
      USE fft_base,      ONLY : grid_gather, grid_scatter

      !
      IMPLICIT NONE
      !
      ! ... Declares variables
      !
      INTEGER, INTENT(IN)     :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin, &
                                 ngm, gstart, nl( ngm ), nlm( ngm)
      REAL( DP ), INTENT(IN)  :: rho ( nrxx, nspin ), g( ngm ), gg( ngm ), nelec
      REAL( DP )              :: vltot( nrxx )
      REAL( DP)               :: rhocharge, ehart
      !
      ! ... Local variables
      !
      REAL( DP ), ALLOCATABLE :: rhotot( : ), v_dcc (:), vcorr(:), &
                                 aux(:)
      INTEGER :: i
      !

      CALL start_clock( 'correction' ) 
      !
      ALLOCATE( rhotot( nrx1 * nrx2 * nrx3 ), v_dcc( nrx1 * nrx2 *nrx3), &
                vcorr( nrx1 * nrx2 * nrx3))
#ifdef __PARA
      ALLOCATE( aux( nrxx ) )
      aux(1:nrxx) = rho(1:nrxx, 1)
      IF( nspin==2 ) aux(1:nrxx) = aux(1:nrxx) + rho(1:nrxx, 2)

      rhotot(:) = 0.d0
      CALL grid_gather( aux, rhotot)
      CALL mp_sum( rhotot, intra_pool_comm )

      DEALLOCATE( aux )
#ifdef DEBUG
      ! test
      v_dcc(:) = 0.d0
      CALL grid_gather( vltot, v_dcc)
      CALL mp_sum( v_dcc, intra_pool_comm )
      do i=1, nrx1*nrx2*nrx3
          write(40+me_pool,*) v_dcc(i),0,0
      end do
      write(40+me_pool,*) "END "
      ! end test
#endif
      v_dcc(:) = 0.d0
      CALL grid_gather( vltot, v_dcc)
      CALL mp_sum( v_dcc, intra_pool_comm )
#else
      rhotot(1:nrxx) = rho(1:nrxx, 1)
      IF( nspin==2 ) rhotot(1:nrxx) = rhotot(1:nrxx) + rho(1:nrxx, 2)
      v_dcc (:) = vltot(:)
#endif
      !
      vcorr (:) = vcomp(:)
      !
#ifdef DEBUG
      do i=1, nrx1*nrx2*nrx3
          write(50+me_pool,*) v_dcc(i), vcorr(i), rhotot(i)
      end do
      write(50+me_pool,*) "END "
#endif
      !
      SELECT CASE( TRIM( which_compensation ) )
        !
 
      CASE( 'dcc' )                           
        !
        CALL add_dcc_field( v_dcc, vcorr, rhotot, nelec,           &
                            nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                               nl, nlm, g, gg, ngm, gstart, nspin )
        !
!#define SOLVATION
#ifdef SOLVATION
      CASE( 'solvation' )
        !
        CALL add_solvation_field( v_dcc, vcorr, rhotot, nelec,     &
                            nr1, nr2, nr3, nrx1, nrx2, nrx3,nrxx,  &
                             nl, nlm, g, gg, ngm, gstart )
#endif
        !
      CASE DEFAULT
        !
        WRITE( *, * ) TRIM(which_compensation), ' not implemented'
        WRITE( *, * ) 'Warning: No Charge Compensation'
        !
      END SELECT
      !

      vcomp (:) = vcorr(:)

#ifdef __PARA
      CALL grid_scatter ( v_dcc, vltot )
#else
      vltot (:) = v_dcc(:)
#endif

      DEALLOCATE( rhotot, vcorr, v_dcc )
      !
      CALL stop_clock( 'correction' )
      !
      RETURN
      !
!--------------------------------------------------------------------
      END SUBROUTINE add_ele_corr
!--------------------------------------------------------------------
