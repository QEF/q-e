!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE h_psi( lda, n, m, psi, hpsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the Hamiltonian
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    hpsi  H*psi
  !
  USE kinds,    ONLY : DP
  USE wvfct,    ONLY : gamma_only 
  USE uspp,     ONLY : vkb, nkb
  USE wvfct,    ONLY : igk, g2kin
  USE gsmooth,  ONLY : nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE ldaU,     ONLY : lda_plus_u
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs  
  USE gvect,    ONLY : gstart
#ifdef EXX
  USE exx,      ONLY : vexx
  USE funct,    ONLY : exx_is_active
#endif

  !
  IMPLICIT NONE
  !
  ! ... input/output arguments
  !
  INTEGER     :: lda, n, m
  COMPLEX(DP) :: psi(lda,m) 
  COMPLEX(DP) :: hpsi(lda,m)   
  !
  !
  CALL start_clock( 'h_psi' )
  !  
  IF ( gamma_only ) THEN
     !
     CALL h_psi_gamma( )
     !
  ELSE  
     !
     CALL h_psi_k( )
     !
  END IF  
  !
  CALL stop_clock( 'h_psi' )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE h_psi_gamma( )
       !-----------------------------------------------------------------------
       ! 
       ! ... gamma version
       !
       USE becmod,  ONLY : rbecp
       !
       IMPLICIT NONE
       !
       INTEGER :: ibnd, j
       !
       !
       CALL start_clock( 'init' )
       !
       ! ... Here we apply the kinetic energy (k+G)^2 psi
       !
       DO ibnd = 1, m
          !
          ! ... set to zero the imaginary part of psi at G=0
          ! ... absolutely needed for numerical stability
          !
          IF ( gstart == 2 ) psi(1,ibnd) = CMPLX( DBLE( psi(1,ibnd) ), 0.D0 )
          !
          hpsi(1:n,ibnd) = g2kin(1:n) * psi(1:n,ibnd)
          !
       END DO
       !
       CALL stop_clock( 'init' )
       !
       ! ... Here we add the Hubbard potential times psi
       !
       IF ( lda_plus_u ) CALL vhpsi( lda, n, m, psi, hpsi )
       !
       ! ... the local potential V_Loc psi
       !
       CALL vloc_psi( lda, n, m, psi, vrs(1,current_spin), hpsi )
       !
       ! ... Here the product with the non local potential V_NL psi
       !
       IF ( nkb > 0 ) THEN
          !
          CALL pw_gemm( 'Y', nkb, m, n, vkb, lda, psi, lda, rbecp, nkb )
          !
          CALL add_vuspsi( lda, n, m, psi, hpsi )
          !
       END IF
       !
#ifdef EXX
       IF ( exx_is_active() ) CALL vexx( lda, n, m, psi, hpsi )
#endif
       !
       RETURN
       !
     END SUBROUTINE h_psi_gamma
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE h_psi_k( )
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE wavefunctions_module, ONLY : psic
       USE becmod,               ONLY : becp
       !
       IMPLICIT NONE
       !
       INTEGER :: ibnd, j
       ! counters
       !
       !
       CALL start_clock( 'init' )
       !
       ! ... Here we apply the kinetic energy (k+G)^2 psi
       !
       DO ibnd = 1, m
          !
          hpsi(1:n,ibnd) = g2kin(1:n) * psi(1:n,ibnd)
          !
       END DO
       !
       CALL stop_clock( 'init' )
       !
       ! ... Here we add the Hubbard potential times psi
       !
       IF ( lda_plus_u ) CALL vhpsi( lda, n, m, psi, hpsi )
       !
       ! ... the local potential V_Loc psi. First the psi in real space
       !
       DO ibnd = 1, m
          !
          CALL start_clock( 'firstfft' )
          !
          psic(1:nrxxs) = ( 0.D0, 0.D0 )
          !
          psic(nls(igk(1:n))) = psi(1:n,ibnd)
          !
          CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
          !
          CALL stop_clock( 'firstfft' )
          !
          ! ... product with the potential vrs = (vltot+vr) on the smooth grid
          !
          psic(1:nrxxs) = psic(1:nrxxs) * vrs(1:nrxxs,current_spin)
          !
          ! ... back to reciprocal space
          !
          CALL start_clock( 'secondfft' )
          !
          CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
          !
          ! ... addition to the total product
          !
          hpsi(1:n,ibnd) = hpsi(1:n,ibnd) + psic(nls(igk(1:n)))
          !
          CALL stop_clock( 'secondfft' )
          !
       END DO
       !
       ! ... Here the product with the non local potential V_NL psi
       !
       IF ( nkb > 0 ) THEN
          !
          CALL ccalbec( nkb, lda, n, m, becp, vkb, psi )
          !
          CALL add_vuspsi( lda, n, m, psi, hpsi )
          !
       END IF
       !
#ifdef EXX
       IF ( exx_is_active() ) CALL vexx( lda, n, m, psi, hpsi )
#endif
       !
       RETURN
       !
     END SUBROUTINE h_psi_k     
     !
END SUBROUTINE h_psi
