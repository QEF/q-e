!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!==========================================================================
!                      metaGGA DRIVERS for E and V
!==========================================================================
!
!--------------------------------------------------------------------------
MODULE qe_drivers_mgga
  !------------------------------------------------------------------------
  !! Contains the mGGA drivers of QE that calculate XC energy and potential.
  !
  USE kind_l,      ONLY: DP
  USE dft_par_mod, ONLY: imeta, imetac, rho_threshold_mgga, grho2_threshold_mgga,&
                         tau_threshold_mgga
  USE metagga
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: tau_xc, tau_xc_spin
  !
  !
CONTAINS
!
!---------------------------------------------------------------------------------
SUBROUTINE tau_xc( length, rho, grho2, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !-------------------------------------------------------------------------------
  !! Meta gradient corrections for exchange and correlation - Hartree a.u.
  !! Available cases: M06L and TPSS. Other mGGA functionals can be used
  !! through Libxc.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Number of k-points
  REAL(DP), DIMENSION(length) :: rho
  !! Charge density
  REAL(DP), DIMENSION(length) :: grho2
  !! Square modulus of the density gradient
  REAL(DP), DIMENSION(length) :: tau
  !! Laplacian of the density
  REAL(DP), DIMENSION(length) :: ex
  !! \(E_x = \int e_x(\text{rho},\text{grho}) dr \)
  REAL(DP), DIMENSION(length) :: ec
  !! \(E_x = \int e_x(\text{rho},\text{grho}) dr \)
  REAL(DP), DIMENSION(length) :: v1x
  !! \( D\ E_x\ /\ D\ \text{rho} \)
  REAL(DP), DIMENSION(length) :: v2x
  !! \( D\ E_x\ /\ D( D\ \text{rho}/D\ r_\alpha )/|\nabla\text{rho}| \)
  REAL(DP), DIMENSION(length) :: v3x
  !! \( D\ E_x\ /\ D\ \text{tau} \)
  REAL(DP), DIMENSION(length) :: v1c
  !! \( D\ E_c\ /\ D\ \text{rho} \)
  REAL(DP), DIMENSION(length) :: v2c
  !! \( D\ E_c\ /\ D( D\ \text{rho}/D\ r_\alpha )/|\nabla\text{rho}| \)
  REAL(DP), DIMENSION(length) :: v3c
  !! \( D\ E_c\ /\ D\ \text{tau} \)
  !
  ! ... local variables
  !
  INTEGER :: k
  REAL(DP) :: arho
  !
  v1x=0.d0 ; v2x=0.d0 ; v3x=0.d0 ; ex=0.d0
  v1c=0.d0 ; v2c=0.d0 ; v3c=0.d0 ; ec=0.d0
  !
  DO k = 1, length
    !
    arho = ABS(rho(k))
    !
    IF ( (arho<=rho_threshold_mgga).OR.(grho2(k)<=grho2_threshold_mgga).OR. &
         (ABS(tau(k))<=rho_threshold_mgga) ) CYCLE
    !
    SELECT CASE( imeta )
    CASE( 1 )
       !
       CALL tpsscxc( arho, grho2(k), tau(k), ex(k), ec(k), v1x(k), v2x(k), &
                     v3x(k), v1c(k), v2c(k), v3c(k) )
       !
    CASE( 2 )
       !
       CALL m06lxc(  arho, grho2(k), tau(k), ex(k), ec(k), v1x(k), v2x(k), &
                     v3x(k), v1c(k), v2c(k), v3c(k) )
       !
    CASE DEFAULT
       !
       CALL xclib_error( 'tau_xc', 'This case is not implemented', imeta )
       !
    END SELECT
    !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE tau_xc
!
!
!-------------------------------------------------------------------------------------
SUBROUTINE tau_xc_spin( length, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !-----------------------------------------------------------------------------------
  !! Meta gradient corrections for exchange and correlation - Hartree a.u. Spin 
  !! polarized case.  
  !! Available cases: M06L and TPSS. Other mGGA functionals can be used
  !! through Libxc.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Number of k-points
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: rho
  !! Charge density
  REAL(DP), INTENT(IN), DIMENSION(3,length,2) :: grho
  !! The density gradient
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: tau
  !! Laplacian of the density
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ex
  !! \(E_x = \int e_x(\text{rho},\text{grho}) dr \)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec
  !! \(E_x = \int e_x(\text{rho},\text{grho}) dr \)
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1x
  !! \( D\ E_x\ /\ D\ \text{rho} \)
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v2x
  !! \( D\ E_x\ /\ D( D\ \text{rho}/D\ r_\alpha )/|\nabla\text{rho}| \)
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v3x
  !! \( D\ E_x\ /\ D\ \text{tau} \)
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1c
  !! \( D\ E_c\ /\ D\ \text{rho} \)
  REAL(DP), INTENT(OUT), DIMENSION(3,length,2) :: v2c
  !! \( D\ E_c\ /\ D( D\ \text{rho}/D\ r_\alpha )/|\nabla\text{rho}| \)
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v3c
  !! \( D\ E_c\ /\ D\ \text{tau} \)
  !
  !  ... local variables
  !
  INTEGER :: k, ipol
  REAL(DP) :: rh, zeta, atau, grho2(2), ggrho2
  REAL(DP) :: v2cup, v2cdw
  !
  ex=0.0_DP ; v1x=0.0_DP ; v2x=0.0_DP ; v3x=0.0_DP
  ec=0.0_DP ; v1c=0.0_DP ; v2c=0.0_DP ; v3c=0.0_DP
  !
  ! FIXME: for SCAN, this will be calculated later
  !
  DO k = 1, length
     !
     rh   = rho(k,1) + rho(k,2)
     atau = tau(k,1) + tau(k,2)             ! KE-density in Hartree
     grho2(1) = SUM( grho(:,k,1)**2 )
     grho2(2) = SUM( grho(:,k,2)**2 )
     ggrho2 = ( grho2(1) + grho2(2) ) * 4.0_DP
     !
     IF ( (rh <= rho_threshold_mgga).OR.(ggrho2 <= grho2_threshold_mgga).OR.&
          (ABS(atau) <= tau_threshold_mgga) ) CYCLE
     !
     SELECT CASE( imeta )
     CASE( 1 )
        !
        CALL tpsscx_spin( rho(k,1), rho(k,2), grho2(1), grho2(2), tau(k,1), &
                          tau(k,2), ex(k), v1x(k,1), v1x(k,2), v2x(k,1),    &
                          v2x(k,2), v3x(k,1), v3x(k,2) )
        !
        zeta = (rho(k,1) - rho(k,2)) / rh
        zeta = MAX( MIN( 0.99999999_DP, zeta ), -0.99999999_DP )
        !
        CALL tpsscc_spin( rh, zeta, grho(:,k,1), grho(:,k,2), atau, ec(k), &
                          v1c(k,1), v1c(k,2), v2c(:,k,1), v2c(:,k,2),      &
                          v3c(k,1), v3c(k,2) )
        !
     CASE( 2 )
        !
        CALL m06lxc_spin( rho(k,1), rho(k,2), grho2(1), grho2(2), tau(k,1), &
                          tau(k,2), ex(k), ec(k), v1x(k,1), v1x(k,2),       &
                          v2x(k,1), v2x(k,2), v3x(k,1), v3x(k,2), v1c(k,1), &
                          v1c(k,2), v2cup, v2cdw, v3c(k,1), v3c(k,2) )
        !
        v2c(:,k,1) = v2cup*grho(:,k,1)
        v2c(:,k,2) = v2cdw*grho(:,k,2)
        !
     CASE DEFAULT
        !
        CALL xclib_error( 'tau_xc_spin', 'This case not implemented', imeta )
        !
     END SELECT
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE tau_xc_spin
!
END MODULE qe_drivers_mgga
