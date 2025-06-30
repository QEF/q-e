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
  USE kind_l,               ONLY: DP
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
  USE dft_setting_params,   ONLY: imeta, rho_threshold_mgga, &
                                  grho2_threshold_mgga
  USE metagga
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Number of k-points
  REAL(DP), INTENT(IN) :: rho(length)
  !! Charge density
  REAL(DP), INTENT(IN) :: grho2(length)
  !! Square modulus of the density gradient
  REAL(DP), INTENT(IN) :: tau(length)
  !! Laplacian of the density
  REAL(DP), INTENT(OUT) :: ex(length)
  !! \(E_x = \int e_x(\text{rho},\text{grho}) dr \)
  REAL(DP), INTENT(OUT) :: ec(length)
  !! \(E_x = \int e_x(\text{rho},\text{grho}) dr \)
  REAL(DP), INTENT(OUT) :: v1x(length)
  !! \( D\ E_x\ /\ D\ \text{rho} \)
  REAL(DP), INTENT(OUT) :: v2x(length)
  !! \( D\ E_x\ /\ D( D\ \text{rho}/D\ r_\alpha )/|\nabla\text{rho}| \)
  REAL(DP), INTENT(OUT) :: v3x(length)
  !! \( D\ E_x\ /\ D\ \text{tau} \)
  REAL(DP), INTENT(OUT) :: v1c(length)
  !! \( D\ E_c\ /\ D\ \text{rho} \)
  REAL(DP), INTENT(OUT) :: v2c(1,length,1)
  !! \( D\ E_c\ /\ D( D\ \text{rho}/D\ r_\alpha )/|\nabla\text{rho}| \)
  REAL(DP), INTENT(OUT) :: v3c(length)
  !! \( D\ E_c\ /\ D\ \text{tau} \)
  !
  ! ... local variables
  !
  INTEGER :: k
  REAL(DP) :: arho
  !
#if defined(_OPENACC)
  !$acc data present( rho, grho2, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !$acc parallel loop
#endif
  DO k = 1, length
    !
    arho = ABS(rho(k))
    !
    IF ( (arho<=rho_threshold_mgga).OR.(grho2(k)<=grho2_threshold_mgga).OR. &
         (ABS(tau(k))<=rho_threshold_mgga) ) THEN
      v1x(k)=0.d0 ; v2x(k)=0.d0 ; v3x(k)=0.d0 ; ex(k)=0.d0
      v1c(k)=0.d0 ; v2c(1,k,1)=0.d0 ; v3c(k)=0.d0 ; ec(k)=0.d0
      CYCLE
    ENDIF
    !
    ! ...libxc-like threshold management
    !grho2(k) = MIN( grho2(k), (8.d0*rho(k)*tau(k))**2 )
    !
    SELECT CASE( imeta )
    CASE( 1 )
       !
       CALL tpsscxc( arho, grho2(k), tau(k), ex(k), ec(k), v1x(k), v2x(k), &
                     v3x(k), v1c(k), v2c(1,k,1), v3c(k) )
       !
    CASE( 2 )
       !
       CALL m06lxc( arho, grho2(k), tau(k), ex(k), ec(k), v1x(k), v2x(k), &
                    v3x(k), v1c(k), v2c(1,k,1), v3c(k) )
       !
    CASE DEFAULT
       !
       v1x(k)=0.d0 ; v2x(k)=0.d0 ; v3x(k)=0.d0 ; ex(k)=0.d0
       v1c(k)=0.d0 ; v2c(1,k,1)=0.d0 ; v3c(k)=0.d0 ; ec(k)=0.d0
       !
    END SELECT
    !
  ENDDO
  !
#if defined(_OPENACC)
  !$acc end data
#endif
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
  USE dft_setting_params,   ONLY: imeta, rho_threshold_mgga, &
                                  grho2_threshold_mgga, tau_threshold_mgga
  USE metagga
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Number of k-points
  REAL(DP), INTENT(IN) :: rho(length,2)
  !! Charge density
  REAL(DP), INTENT(IN) :: grho(3,length,2)
  !! The density gradient
  REAL(DP), INTENT(IN) :: tau(length,2)
  !! Laplacian of the density
  REAL(DP), INTENT(OUT) :: ex(length)
  !! \(E_x = \int e_x(\text{rho},\text{grho}) dr \)
  REAL(DP), INTENT(OUT) :: ec(length)
  !! \(E_x = \int e_x(\text{rho},\text{grho}) dr \)
  REAL(DP), INTENT(OUT) :: v1x(length,2)
  !! \( D\ E_x\ /\ D\ \text{rho} \)
  REAL(DP), INTENT(OUT) :: v2x(length,2)
  !! \( D\ E_x\ /\ D( D\ \text{rho}/D\ r_\alpha )/|\nabla\text{rho}| \)
  REAL(DP), INTENT(OUT) :: v3x(length,2)
  !! \( D\ E_x\ /\ D\ \text{tau} \)
  REAL(DP), INTENT(OUT) :: v1c(length,2)
  !! \( D\ E_c\ /\ D\ \text{rho} \)
  REAL(DP), INTENT(OUT) :: v2c(3,length,2)
  !! \( D\ E_c\ /\ D( D\ \text{rho}/D\ r_\alpha )/|\nabla\text{rho}| \)
  REAL(DP), INTENT(OUT) :: v3c(length,2)
  !! \( D\ E_c\ /\ D\ \text{tau} \)
  !
  !  ... local variables
  !
  INTEGER :: k
  REAL(DP) :: rh, zeta, atau, grho2up, grho2dw, ggrho2
  REAL(DP) :: v2cup, v2cdw
  !
#if defined(_OPENACC)
  !$acc data present( rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !$acc parallel loop
#endif
  DO k = 1, length
     !
     rh   = rho(k,1) + rho(k,2)
     atau = tau(k,1) + tau(k,2)             ! KE-density in Hartree
     ! ...libxc-like threshold management
     !grho2up = MIN( SUM(grho(:,k,1)**2), (8.d0*rho(k,1)*tau(k,1))**2 )
     !grho2dw = MIN( SUM(grho(:,k,2)**2), (8.d0*rho(k,2)*tau(k,2))**2 )
     grho2up = SUM(grho(:,k,1)**2) 
     grho2dw = SUM(grho(:,k,2)**2)
     ggrho2 = ( grho2up + grho2dw ) * 4.0_DP
     !
     IF ( (rh <= rho_threshold_mgga).OR.(ggrho2 <= grho2_threshold_mgga).OR.&
          (ABS(atau) <= tau_threshold_mgga) ) THEN
       v1x(k,:)=0.d0 ; v2x(k,:)=0.d0   ; v3x(k,:)=0.d0 ; ex(k)=0.d0
       v1c(k,:)=0.d0 ; v2c(:,k,:)=0.d0 ; v3c(k,:)=0.d0 ; ec(k)=0.d0
       CYCLE
     ENDIF
     !
     SELECT CASE( imeta )
     CASE( 1 )
        !
        CALL tpsscx_spin( rho(k,1), rho(k,2), grho2up, grho2dw, tau(k,1), &
                          tau(k,2), ex(k), v1x(k,1), v1x(k,2), v2x(k,1),    &
                          v2x(k,2), v3x(k,1), v3x(k,2) )
        !
        zeta = MAX( MIN( 0.99999999_DP, (rho(k,1)-rho(k,2))/rh ), -0.99999999_DP )
        !
        CALL tpsscc_spin( rh, zeta, grho(:,k,1), grho(:,k,2), atau, ec(k), &
                          v1c(k,1), v1c(k,2), v2c(:,k,1), v2c(:,k,2),      &
                          v3c(k,1), v3c(k,2) )
        !
     CASE( 2 )
        !
        CALL m06lxc_spin( rho(k,1), rho(k,2), grho2up, grho2dw, tau(k,1),   &
                          tau(k,2), ex(k), ec(k), v1x(k,1), v1x(k,2),       &
                          v2x(k,1), v2x(k,2), v3x(k,1), v3x(k,2), v1c(k,1), &
                          v1c(k,2), v2cup, v2cdw, v3c(k,1), v3c(k,2) )
        !
        v2c(:,k,1) = v2cup*grho(:,k,1)
        v2c(:,k,2) = v2cdw*grho(:,k,2)
        !
     CASE DEFAULT
        !
        v1x(k,:)=0.d0 ; v2x(k,:)=0.d0   ; v3x(k,:)=0.d0 ; ex(k)=0.d0
        v1c(k,:)=0.d0 ; v2c(:,k,:)=0.d0 ; v3c(k,:)=0.d0 ; ec(k)=0.d0
        !
     END SELECT
     !
  ENDDO
  !
#if defined(_OPENACC)
  !$acc end data
#endif
  !
  RETURN
  !
END SUBROUTINE tau_xc_spin
!
END MODULE qe_drivers_mgga
