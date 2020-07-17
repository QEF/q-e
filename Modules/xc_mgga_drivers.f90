MODULE xc_mgga
!
USE kinds,     ONLY: DP
USE funct,     ONLY: get_meta, get_metac, is_libxc, &
                     exx_is_active, scan_exx, get_exx_fraction
!
IMPLICIT NONE
!
PRIVATE
SAVE
!
!  GGA exchange-correlation drivers
PUBLIC :: xc_metagcx
PUBLIC :: tau_xc, tau_xc_spin
PUBLIC :: change_threshold_mgga
!
!
!  input thresholds (default values)
REAL(DP) :: rho_threshold   = 1.0E-12_DP
REAL(DP) :: grho2_threshold = 1.0E-24_DP
REAL(DP) :: tau_threshold   = 1.0E-12_DP
!
!
 CONTAINS
!
!
!-------------------------------------------------------------------------------------
SUBROUTINE change_threshold_mgga( rho_thr_in, grho2_thr_in, tau_thr_in )
  !------------------------------------------------------------------------------------
  !! Change rho, grho and tau thresholds.
  ! 
  REAL(DP), INTENT(IN) :: rho_thr_in
  REAL(DP), INTENT(IN), OPTIONAL :: grho2_thr_in
  REAL(DP), INTENT(IN), OPTIONAL :: tau_thr_in
  !
  rho_threshold = rho_thr_in
  IF ( PRESENT(grho2_thr_in) ) grho2_threshold = grho2_thr_in
  IF ( PRESENT(tau_thr_in)  ) tau_threshold  = tau_thr_in
  !
  RETURN
  !
END SUBROUTINE change_threshold_mgga
!
!
!----------------------------------------------------------------------------------------
SUBROUTINE xc_metagcx( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !-------------------------------------------------------------------------------------
  !! Wrapper routine. Calls metaGGA drivers from internal libraries
  !! of q-e or from the external libxc, depending on the input choice.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE funct,            ONLY : get_libxc_flags_exc
  USE xc_f03_lib_m
#endif 
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
  INTEGER :: is, imeta, imetac
  REAL(DP), ALLOCATABLE :: grho2(:,:)
  !
#if defined(__LIBXC)
  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info1, xc_info2
  !
  REAL(DP), ALLOCATABLE :: rho_lxc(:), sigma(:), tau_lxc(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vx_sigma(:), vx_tau(:)
  REAL(DP), ALLOCATABLE :: vc_rho(:), vc_sigma(:), vc_tau(:)
  REAL(DP), ALLOCATABLE :: lapl_rho(:), vlapl_rho(:) ! not used in TPSS
  !
  REAL(DP) :: exx_fraction
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
  imeta  = get_meta()
  imetac = get_metac()
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
                      grho2_threshold )
      tau_lxc(k) = MAX( tau(k,1), tau_threshold )
    ENDDO
    !
  ELSE
    !
    DO k = 1, length
       rho_lxc(2*k-1) = ABS( rho(k,1) )
       rho_lxc(2*k)   = ABS( rho(k,2) )
       !
       sigma(3*k-2) = MAX( grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2, &
                           grho2_threshold )
       sigma(3*k-1) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                      grho(3,k,1) * grho(3,k,2)
       sigma(3*k)   = MAX( grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2, &
                           grho2_threshold )
       !
       tau_lxc(2*k-1) = MAX( tau(k,1), tau_threshold )
       tau_lxc(2*k)   = MAX( tau(k,2), tau_threshold )
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
     CALL xc_f03_func_init( xc_func, imeta, pol_unpol )
     xc_info1 = xc_f03_func_get_info( xc_func )
     CALL xc_f03_func_set_dens_threshold( xc_func, rho_threshold )
     CALL get_libxc_flags_exc( xc_info1, eflag )
     IF (eflag==1) THEN
       CALL xc_f03_mgga_exc_vxc( xc_func, lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                                 ex_lxc(1), vx_rho(1), vx_sigma(1), vlapl_rho(1), vx_tau(1) )
     ELSE
       CALL xc_f03_mgga_vxc( xc_func, lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                             vx_rho(1), vx_sigma(1), vlapl_rho(1), vx_tau(1) )
     ENDIF
    CALL xc_f03_func_end( xc_func )
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
       exx_fraction = get_exx_fraction()
       IF (exx_is_active()) THEN
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
    CALL xc_f03_func_init( xc_func, imetac, pol_unpol )
    xc_info1 = xc_f03_func_get_info( xc_func )
    CALL xc_f03_func_set_dens_threshold( xc_func, rho_threshold )
    CALL xc_f03_mgga_exc_vxc( xc_func, lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                               ec_lxc(1), vc_rho(1), vc_sigma(1), vlapl_rho(1), vc_tau(1) )
    CALL xc_f03_func_end( xc_func )
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
!
!
!---------------------------------------------------------------------------------
SUBROUTINE tau_xc( length, rho, grho2, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
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
  USE metagga
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !
  INTEGER :: k, imeta
  REAL(DP) :: arho
  REAL(DP), DIMENSION(length) :: rho, grho2, tau, &
                                 ex, ec, v1x, v2x, v3x, v1c, v2c, v3c
  !
  imeta = get_meta()
  !
  v1x=0.d0 ; v2x=0.d0 ; v3x=0.d0 ; ex=0.d0
  v1c=0.d0 ; v2c=0.d0 ; v3c=0.d0 ; ec=0.d0
  !
  DO k = 1, length
    !
    arho = ABS(rho(k))
    !
    IF ( (arho<=rho_threshold).OR.(grho2(k)<=grho2_threshold).OR.(ABS(tau(k))<=rho_threshold) ) CYCLE
    !
    SELECT CASE( imeta )
    CASE( 1 )
       CALL tpsscxc( arho, grho2(k), tau(k), ex(k), ec(k), v1x(k), v2x(k), v3x(k), v1c(k), v2c(k), v3c(k) )
    CASE( 2 )
       CALL m06lxc(  arho, grho2(k), tau(k), ex(k), ec(k), v1x(k), v2x(k), v3x(k), v1c(k), v2c(k), v3c(k) )
    CASE DEFAULT
       CALL errore( 'tau_xc', 'This case is not implemented', imeta )
    END SELECT
    !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE tau_xc
!
!
!----------------------------------------------------------------------------------------
SUBROUTINE tau_xc_spin( length, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !------------------------------------------------------------------------------------
  !
  USE metagga
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
  INTEGER :: k, ipol, imeta
  REAL(DP) :: rh, zeta, atau, grho2(2), ggrho2
  REAL(DP) :: v2cup, v2cdw
  !
  imeta = get_meta()
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
     IF ((rh <= rho_threshold) .OR. (ggrho2 <= grho2_threshold) .OR. (ABS(atau) <= tau_threshold)) CYCLE
     !
     SELECT CASE( imeta )
     CASE( 1 )
        !
        CALL tpsscx_spin( rho(k,1), rho(k,2), grho2(1), grho2(2), tau(k,1), &
                          tau(k,2), ex(k), v1x(k,1), v1x(k,2), v2x(k,1), v2x(k,2), v3x(k,1), v3x(k,2) )
        !
        zeta = (rho(k,1) - rho(k,2)) / rh
        zeta = MAX( MIN( 0.99999999_DP, zeta ), -0.99999999_DP )
        !
        CALL tpsscc_spin( rh, zeta, grho(:,k,1), grho(:,k,2), atau, ec(k), &
                          v1c(k,1), v1c(k,2), v2c(:,k,1), v2c(:,k,2), v3c(k,1), v3c(k,2) )
        !
     CASE( 2 )
        !
        CALL m06lxc_spin( rho(k,1), rho(k,2), grho2(1), grho2(2), tau(k,1), tau(k,2), ex(k), ec(k), &
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
END SUBROUTINE tau_xc_spin
!
!
END MODULE xc_mgga
