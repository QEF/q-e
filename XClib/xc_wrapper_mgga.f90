!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------------------
SUBROUTINE xc_metagcx( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !-------------------------------------------------------------------------------------
  !! Wrapper routine. Calls internal metaGGA drivers or the Libxc ones,
  !! depending on the input choice.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_mod,       ONLY: get_libxc_flags_exc
  USE dft_par_mod,   ONLY: xc_func, xc_info
#endif 
  !
  USE kind_l,        ONLY: DP
  USE dft_par_mod,   ONLY: imeta, imetac, is_libxc, rho_threshold_mgga,        &
                           grho2_threshold_mgga, tau_threshold_mgga, scan_exx, &
                           exx_started, exx_fraction
  USE qe_drivers_mgga
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER, INTENT(IN) :: ns
  !! spin components
  INTEGER, INTENT(IN) :: np
  !! first dimension of v2c
  REAL(DP), INTENT(IN) :: rho(length,ns)
  !! the charge density
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  !! grho = \nabla rho
  REAL(DP), INTENT(IN) :: tau(length,ns)
  !! kinetic energy density
  REAL(DP), INTENT(OUT) :: ex(length)
  !! sx = E_x(rho,grho)
  REAL(DP), INTENT(OUT) :: ec(length)
  !! sc = E_c(rho,grho)
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  !! v1x = D(E_x)/D(rho)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  !! v2x = D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3x(length,ns)
  !! v3x = D(E_x)/D(tau)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  !! v1c = D(E_c)/D(rho)
  REAL(DP), INTENT(OUT) :: v2c(np,length,ns)
  !! v2c = D(E_c)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3c(length,ns)
  !! v3c = D(E_c)/D(tau)
  !
  ! ... local variables
  !
  INTEGER :: is
  REAL(DP), ALLOCATABLE :: grho2(:,:)
  !
#if defined(__LIBXC)
  REAL(DP), ALLOCATABLE :: rho_lxc(:), sigma(:), tau_lxc(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vx_sigma(:), vx_tau(:)
  REAL(DP), ALLOCATABLE :: vc_rho(:), vc_sigma(:), vc_tau(:)
  REAL(DP), ALLOCATABLE :: lapl_rho(:), vlapl_rho(:) ! not used in TPSS
  !
  INTEGER :: k, ipol, pol_unpol, eflag
  LOGICAL :: POLARIZED
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
  !
  lengthxc = length
  !
  ex = 0.0_DP ;  v1x = 0.0_DP ;  v2x = 0.0_DP ;  v3x = 0.0_DP
  ec = 0.0_DP ;  v1c = 0.0_DP ;  v2c = 0.0_DP ;  v3c = 0.0_DP
  !
  POLARIZED = .FALSE.
  IF (ns == 2) THEN
     POLARIZED = .TRUE.
  ENDIF
  !
  pol_unpol = ns
  !
  ALLOCATE( rho_lxc(length*ns), sigma(length*np), tau_lxc(length*ns) )
  ALLOCATE( lapl_rho(length*ns) )
  !
  ALLOCATE( ex_lxc(length)    , ec_lxc(length) )
  ALLOCATE( vx_rho(length*ns) , vx_sigma(length*np), vx_tau(length*ns) )
  ALLOCATE( vc_rho(length*ns) , vc_sigma(length*np), vc_tau(length*ns) )
  ALLOCATE( vlapl_rho(length*ns) )
  !
  !
  IF ( ns == 1 ) THEN
    !
    DO k = 1, length
      rho_lxc(k) = ABS( rho(k,1) )
      sigma(k) = MAX( grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2, &
                      grho2_threshold_mgga )
      tau_lxc(k) = MAX( tau(k,1), tau_threshold_mgga )
    ENDDO
    !
  ELSE
    !
    DO k = 1, length
       rho_lxc(2*k-1) = ABS( rho(k,1) )
       rho_lxc(2*k)   = ABS( rho(k,2) )
       !
       sigma(3*k-2) = MAX( grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2, &
                           grho2_threshold_mgga )
       sigma(3*k-1) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                      grho(3,k,1) * grho(3,k,2)
       sigma(3*k)   = MAX( grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2, &
                           grho2_threshold_mgga )
       !
       tau_lxc(2*k-1) = MAX( tau(k,1), tau_threshold_mgga )
       tau_lxc(2*k)   = MAX( tau(k,2), tau_threshold_mgga )
    ENDDO
    !
  ENDIF
  !
  IF ( .NOT.is_libxc(5) .AND. imetac==0 ) THEN
    !
    ALLOCATE( grho2(length,ns) )
    !
    DO is = 1, ns
       grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
    ENDDO
    !
    IF (ns == 1) THEN
       CALL tau_xc( length, rho(:,1), grho2(:,1), tau(:,1), ex, ec, v1x(:,1), &
                    v2x(:,1), v3x(:,1), v1c(:,1), v2c(1,:,1), v3c(:,1) )
    ELSEIF (ns == 2) THEN
       CALL tau_xc_spin( length, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, &
                         v2c, v3c )
    ENDIF
    !
    DEALLOCATE( grho2 )
    !
  ENDIF
  !
  ! META EXCHANGE
  !
  IF ( is_libxc(5) ) THEN
    CALL xc_f03_func_set_dens_threshold( xc_func(5), rho_threshold_mgga )
    CALL get_libxc_flags_exc( xc_info(5), eflag )
    IF (eflag==1) THEN
      CALL xc_f03_mgga_exc_vxc( xc_func(5), lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                                ex_lxc(1), vx_rho(1), vx_sigma(1), vlapl_rho(1), vx_tau(1) )
    ELSE
      CALL xc_f03_mgga_vxc( xc_func(5), lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                            vx_rho(1), vx_sigma(1), vlapl_rho(1), vx_tau(1) )
    ENDIF
    !
    IF (.NOT. POLARIZED) THEN
      DO k = 1, length
        ex(k) = ex_lxc(k) * rho_lxc(k)
        v1x(k,1) = vx_rho(k)
        v2x(k,1) = vx_sigma(k) * 2.0_DP
        v3x(k,1) = vx_tau(k)
      ENDDO
    ELSE
      DO k = 1, length
        ex(k) = ex_lxc(k) * (rho_lxc(2*k-1)+rho_lxc(2*k))
        v1x(k,1) = vx_rho(2*k-1)
        v1x(k,2) = vx_rho(2*k)
        v2x(k,1) = vx_sigma(3*k-2)*2.d0
        v2x(k,2) = vx_sigma(3*k)*2.d0
        v3x(k,1) = vx_tau(2*k-1)
        v3x(k,2) = vx_tau(2*k)
      ENDDO
    ENDIF
    !
    ! ... only for HK/MCA: SCAN0 (used in CPV)
    IF ( scan_exx ) THEN
       IF (exx_started) THEN
         ex  = (1.0_DP - exx_fraction) * ex
         v1x = (1.0_DP - exx_fraction) * v1x
         v2x = (1.0_DP - exx_fraction) * v2x
         v3x = (1.0_DP - exx_fraction) * v3x
       ENDIF
    ENDIF
    !
  ENDIF
  !
  ! META CORRELATION
  !
  IF ( is_libxc(6) ) THEN
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(6), rho_threshold_mgga )
    CALL xc_f03_mgga_exc_vxc( xc_func(6), lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                              ec_lxc(1), vc_rho(1), vc_sigma(1), vlapl_rho(1), vc_tau(1) )
    !
    IF (.NOT. POLARIZED) THEN
       DO k = 1, length
          ec(k) = ec_lxc(k) * rho_lxc(k) 
          v1c(k,1) = vc_rho(k)
          v2c(1,k,1) = vc_sigma(k) * 2.0_DP
          v3c(k,1) = vc_tau(k)
       ENDDO
    ELSE
       DO k = 1, length
          ec(k) = ec_lxc(k) * (rho_lxc(2*k-1)+rho_lxc(2*k))
          v1c(k,1) = vc_rho(2*k-1)
          v1c(k,2) = vc_rho(2*k)
          DO ipol = 1, 3
            v2c(ipol,k,1) = vc_sigma(3*k-2)*grho(ipol,k,1)*2.D0 + vc_sigma(3*k-1)*grho(ipol,k,2)
            v2c(ipol,k,2) = vc_sigma(3*k)  *grho(ipol,k,2)*2.D0 + vc_sigma(3*k-1)*grho(ipol,k,1)
          ENDDO
          v3c(k,1) = vc_tau(2*k-1)
          v3c(k,2) = vc_tau(2*k)
       ENDDO
    ENDIF
    !
  ENDIF
  !
  DEALLOCATE( rho_lxc, sigma, tau_lxc, lapl_rho )
  DEALLOCATE( ex_lxc , ec_lxc )
  DEALLOCATE( vx_rho , vx_sigma, vx_tau )
  DEALLOCATE( vc_rho , vc_sigma, vc_tau, vlapl_rho )
  !
#else
  !
  ALLOCATE( grho2(length,ns) )
  !
  DO is = 1, ns
     grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
  ENDDO
  !
  IF (ns == 1) THEN
     !
     CALL tau_xc( length, rho(:,1), grho2(:,1), tau(:,1), ex, ec, v1x(:,1), &
                  v2x(:,1), v3x(:,1), v1c(:,1), v2c(1,:,1), v3c(:,1) )
     !
  ELSEIF (ns == 2) THEN
     !
     CALL tau_xc_spin( length, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, &
                       v2c, v3c )
     !
  ENDIF
  !
  DEALLOCATE( grho2 )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE xc_metagcx
