!---------------------------------------------------------------------------------
SUBROUTINE tau_xc_l( length, rho, grho2, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !-------------------------------------------------------------------------------
  !  gradient corrections for exchange and correlation - Hartree a.u.
  !  See comments at the beginning of module for implemented cases
  !
  !  input:  rho, grho=|\nabla rho|^2
  !
  !  definition:  E_x = \int e_x(rho,grho) dr
  !
  !  output: sx = e_x(rho,grho) = grad corr
  !          v1x= D(E_x)/D(rho)
  !          v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !          v3x= D(E_x)/D(tau)
  !
  !          sc, v1c, v2c as above for correlation
  !
  USE kind_l
  USE dft_par_mod
  USE metagga_l
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !
  INTEGER :: k
  REAL(DP) :: arho
  REAL(DP), DIMENSION(length) :: rho, grho2, tau, &
                                 ex, ec, v1x, v2x, v3x, v1c, v2c, v3c
  !
  v1x=0.d0 ; v2x=0.d0 ; v3x=0.d0 ; ex=0.d0
  v1c=0.d0 ; v2c=0.d0 ; v3c=0.d0 ; ec=0.d0
  !
  DO k = 1, length
    !
    arho = ABS(rho(k))
    !
    IF ( (arho<=rho_threshold_mgga).OR.(grho2(k)<=grho2_threshold_mgga).OR.(ABS(tau(k))<=rho_threshold_mgga) ) CYCLE
    !
    SELECT CASE( imeta )
    CASE( 1 )
       CALL tpsscxc_l( arho, grho2(k), tau(k), ex(k), ec(k), v1x(k), v2x(k), v3x(k), v1c(k), v2c(k), v3c(k) )
    CASE( 2 )
       CALL m06lxc_l(  arho, grho2(k), tau(k), ex(k), ec(k), v1x(k), v2x(k), v3x(k), v1c(k), v2c(k), v3c(k) )
    CASE DEFAULT
       CALL errore( 'tau_xc', 'This case is not implemented', imeta )
    END SELECT
    !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE tau_xc_l
!
!
!----------------------------------------------------------------------------------------
SUBROUTINE tau_xc_spin_l( length, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !------------------------------------------------------------------------------------
  !
  USE kind_l
  USE dft_par_mod
  USE metagga_l
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN) :: rho(length,2), tau(length,2)
  REAL(DP), INTENT(IN) :: grho(3,length,2)
  !
  REAL(DP), INTENT(OUT) :: ex(length), ec(length), v1x(length,2), v2x(length,2), &
                           v3x(length,2), v1c(length,2), v3c(length,2)
  REAL(DP), INTENT(OUT) :: v2c(3,length,2)
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
     IF ((rh <= rho_threshold_mgga) .OR. (ggrho2 <= grho2_threshold_mgga) .OR. (ABS(atau) <= tau_threshold_mgga)) CYCLE
     !
     SELECT CASE( imeta )
     CASE( 1 )
        !
        CALL tpsscx_spin_l( rho(k,1), rho(k,2), grho2(1), grho2(2), tau(k,1), &
                          tau(k,2), ex(k), v1x(k,1), v1x(k,2), v2x(k,1), v2x(k,2), v3x(k,1), v3x(k,2) )
        !
        zeta = (rho(k,1) - rho(k,2)) / rh
        zeta = MAX( MIN( 0.99999999_DP, zeta ), -0.99999999_DP )
        !
        CALL tpsscc_spin_l( rh, zeta, grho(:,k,1), grho(:,k,2), atau, ec(k), &
                          v1c(k,1), v1c(k,2), v2c(:,k,1), v2c(:,k,2), v3c(k,1), v3c(k,2) )
        !
     CASE( 2 )
        !
        CALL m06lxc_spin_l( rho(k,1), rho(k,2), grho2(1), grho2(2), tau(k,1), tau(k,2), ex(k), ec(k), &
                          v1x(k,1), v1x(k,2), v2x(k,1), v2x(k,2), v3x(k,1), v3x(k,2), &
                          v1c(k,1), v1c(k,2), v2cup   , v2cdw   , v3c(k,1), v3c(k,2)  )
        !
        v2c(:,k,1) = v2cup*grho(:,k,1)
        v2c(:,k,2) = v2cdw*grho(:,k,2)
        !
     CASE DEFAULT
        !
        CALL errore( 'tau_xc_spin', 'This case not implemented', imeta )
        !
     END SELECT
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE tau_xc_spin_l
