!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------------------
SUBROUTINE xc_metagcx( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, &
                       v2c, v3c, gpu_args_ )
  !----------------------------------------------------------------------------------
  !! Wrapper to gpu or non gpu version of \(\texttt{xc_metagcx}\).
  !
  USE kind_l,               ONLY: DP
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
  LOGICAL, INTENT(IN), OPTIONAL :: gpu_args_
  !! whether you wish to run on gpu in case use_gpu is true
  !
  LOGICAL :: gpu_args
  !
  gpu_args = .FALSE.
  IF ( PRESENT(gpu_args_) ) gpu_args = gpu_args_
  !
  IF ( gpu_args ) THEN
    !
    !$acc data present( rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
    CALL xc_metagcx_( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, &
                      v2c, v3c )
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( rho, grho, tau ), copyout( ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
    CALL xc_metagcx_( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, &
                      v2c, v3c )
    !$acc end data
    !
  ENDIF  
  !
  RETURN
  !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------------------
SUBROUTINE xc_metagcx_( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !-------------------------------------------------------------------------------------
  !! Wrapper routine. Calls internal metaGGA drivers or the Libxc ones,
  !! depending on the input choice.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_setting_params,   ONLY: xc_func, libxc_flags
#endif 
  !
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: imeta, imetac, is_libxc, rho_threshold_mgga,&
                                  grho2_threshold_mgga, tau_threshold_mgga,   &
                                  ishybrid, exx_started, exx_fraction
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
  INTEGER :: k, is, ipol
  REAL(DP), ALLOCATABLE :: grho2(:,:)
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
#if defined(__LIBXC)
  REAL(DP), ALLOCATABLE :: rho_lxc(:), sigma(:), tau_lxc(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vx_sigma(:), vx_tau(:)
  REAL(DP), ALLOCATABLE :: vc_rho(:), vc_sigma(:), vc_tau(:)
  REAL(DP), ALLOCATABLE :: lapl_rho(:), vlapl_rho(:) ! not used in QE
  !
  REAL(DP) :: rh, ggrho2, atau, xcoef
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  !$acc data present( rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !
#if defined(__LIBXC)
  lengthxc = length
  !
  ALLOCATE( rho_lxc(length*ns), sigma(length*np) )
  ALLOCATE( tau_lxc(length*ns), lapl_rho(length*ns) )
  !$acc data create( rho_lxc, sigma, tau_lxc, lapl_rho )
  !
  IF ( is_libxc(5) ) THEN
    ALLOCATE( ex_lxc(length) )
    ALLOCATE( vx_rho(length*ns), vx_sigma(length*np), vx_tau(length*ns) )
    IF ( imetac==0 ) THEN
      !$acc parallel loop
      DO k = 1, length
        ec(k) = 0.d0
        DO is = 1, ns
          v1c(k,is) = 0.d0 ; v3c(k,is) = 0.d0
          DO ipol = 1, np
            v2c(ipol,k,is) = 0.d0
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  IF ( is_libxc(6) ) THEN
    ALLOCATE( ec_lxc(length) )
    ALLOCATE( vc_rho(length*ns), vc_sigma(length*np), vc_tau(length*ns) )
    IF ( imeta==0 ) THEN
      !$acc parallel loop
      DO k = 1, length
        ex(k) = 0.d0
        DO is = 1, ns
          v1x(k,is) = 0.d0 ; v2x(k,is) = 0.d0 ; v3x(k,is) = 0.d0
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  IF ( ANY(is_libxc(5:6)) ) ALLOCATE( vlapl_rho(length*ns) )
  !
  IF ( ns == 1 ) THEN
    !
    !$acc parallel loop
    DO k = 1, length
      rho_lxc(k) = ABS( rho(k,1) )
      sigma(k) = MAX( grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2, &
                      grho2_threshold_mgga )
      tau_lxc(k) = MAX( tau(k,1), tau_threshold_mgga )
      lapl_rho(k) = 0.d0
    ENDDO
    !
  ELSE
    !
    !$acc parallel loop
    DO k = 1, length
       rho_lxc(2*k-1) = ABS( rho(k,1) )
       rho_lxc(2*k)   = ABS( rho(k,2) )
       sigma(3*k-2) = MAX( grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2, &
                           grho2_threshold_mgga )
       sigma(3*k-1) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) +&
                      grho(3,k,1) * grho(3,k,2)
       sigma(3*k)   = MAX( grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2, &
                           grho2_threshold_mgga )
       tau_lxc(2*k-1) = MAX( tau(k,1), small )
       tau_lxc(2*k)   = MAX( tau(k,2), small )
       lapl_rho(2*k-1) = 0.d0
       lapl_rho(2*k)   = 0.d0
    ENDDO
    !
  ENDIF
  !
  !$acc update self( rho_lxc, sigma, tau_lxc, lapl_rho )
  !
#endif
  !
  IF ( .NOT.is_libxc(5) .AND. imetac==0 ) THEN
    IF (ns == 1) THEN
      !
      ALLOCATE( grho2(length,ns) )
      !$acc data create( grho2 )
      !
      !$acc parallel loop
      DO k = 1, length
        grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
      ENDDO
      !
      CALL tau_xc( length, rho(:,1), grho2, tau(:,1), ex, ec, v1x(:,1), &
                   v2x(:,1), v3x(:,1), v1c(:,1), v2c, v3c(:,1) )
      !
      !$acc end data
      DEALLOCATE( grho2 )
      !
    ELSEIF (ns == 2) THEN
      !
      CALL tau_xc_spin( length, rho, grho, tau, ex, ec, v1x, v2x, v3x, &
                        v1c, v2c, v3c )
      !
    ENDIF
  ENDIF
  !
#if defined(__LIBXC)
  !
  ! ... META EXCHANGE
  !
  IF ( is_libxc(5) ) THEN
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(5), rho_threshold_mgga )
    IF (libxc_flags(5,0)==1) THEN
      CALL xc_f03_mgga_exc_vxc( xc_func(5), lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                                ex_lxc(1), vx_rho(1), vx_sigma(1), vlapl_rho(1), vx_tau(1) )
    ELSE
      CALL xc_f03_mgga_vxc( xc_func(5), lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                            vx_rho(1), vx_sigma(1), vlapl_rho(1), vx_tau(1) )
      ex_lxc = 0.d0 
    ENDIF
    !
    xcoef = 1.d0
    IF ( ishybrid .AND. exx_started .AND. exx_fraction>0.d0) xcoef = 1.d0-exx_fraction
    !
    !$acc data copyin( ex_lxc, vx_rho, vx_sigma, vx_tau )
    IF ( ns==1 ) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( ABS(rho_lxc(k))<=rho_threshold_mgga .OR. &
             sigma(k)<=grho2_threshold_mgga      .OR. &
             ABS(tau_lxc(k))<=rho_threshold_mgga ) THEN
          ex(k) = 0.d0    ; v1x(k,1) = 0.d0
          v2x(k,1) = 0.d0 ; v3x(k,1) = 0.d0
          CYCLE
        ENDIF  
        ex(k) = xcoef * ex_lxc(k) * rho_lxc(k)
        v1x(k,1) = xcoef * vx_rho(k)
        v2x(k,1) = xcoef * vx_sigma(k) * 2.0_DP
        v3x(k,1) = xcoef * vx_tau(k)
      ENDDO
    ELSE
      !$acc parallel loop
      DO k = 1, length
        IF (rho_lxc(2*k-1)+rho_lxc(2*k) <= rho_threshold_mgga) THEN
          ex(k) = 0.d0    
          v1x(k,1) = 0.d0 ; v2x(k,1) = 0.d0 ; v3x(k,1) = 0.d0
          v1x(k,2) = 0.d0 ; v2x(k,2) = 0.d0 ; v3x(k,2) = 0.d0
          CYCLE
        ENDIF
        ex(k) = xcoef * ex_lxc(k) * (rho_lxc(2*k-1)+rho_lxc(2*k))
        IF ( ABS(rho_lxc(2*k-1))>rho_threshold_mgga .AND. &
             sigma(3*k-2)>grho2_threshold_mgga      .AND. &
             ABS(tau_lxc(2*k-1))>tau_threshold_mgga ) THEN
          v1x(k,1) = xcoef * vx_rho(2*k-1)
          v2x(k,1) = xcoef * vx_sigma(3*k-2)*2.d0
          v3x(k,1) = xcoef * vx_tau(2*k-1)
        ELSE
          v1x(k,1) = 0.d0 ; v2x(k,1) = 0.d0 ; v3x(k,1) = 0.d0
        ENDIF
        IF ( ABS(rho_lxc(2*k))>rho_threshold_mgga .AND. &
             sigma(3*k)>grho2_threshold_mgga      .AND. &
             ABS(tau_lxc(2*k))>tau_threshold_mgga ) THEN
          v1x(k,2) = xcoef * vx_rho(2*k)
          v2x(k,2) = xcoef * vx_sigma(3*k)*2.d0
          v3x(k,2) = xcoef * vx_tau(2*k)
        ELSE
          v1x(k,2) = 0.d0 ; v2x(k,2) = 0.d0 ; v3x(k,2) = 0.d0
        ENDIF
      ENDDO
    ENDIF
    !
    !$acc end data
    DEALLOCATE( ex_lxc, vx_rho, vx_sigma, vx_tau )
    !
  ENDIF
  !
  ! ... META CORRELATION
  !
  IF ( is_libxc(6) ) THEN
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(6), rho_threshold_mgga )
    IF (libxc_flags(6,0)==1) THEN
      CALL xc_f03_mgga_exc_vxc( xc_func(6), lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                                ec_lxc(1), vc_rho(1), vc_sigma(1), vlapl_rho(1), vc_tau(1) )
    ELSE
      CALL xc_f03_mgga_vxc( xc_func(6), lengthxc, rho_lxc(1), sigma(1), lapl_rho(1), tau_lxc(1), &
                            vc_rho(1), vc_sigma(1), vlapl_rho(1), vc_tau(1) )
      ec_lxc = 0.d0
    ENDIF
    !
    !$acc data copyin( ec_lxc, vc_rho, vc_sigma, vc_tau )
    IF ( ns==1 ) THEN
       !$acc parallel loop
       DO k = 1, length
         IF ( ABS(rho_lxc(k))<=rho_threshold_mgga    .OR. &
                     sigma(k)<=grho2_threshold_mgga  .OR. &
              ABS(tau_lxc(k))<=rho_threshold_mgga  ) THEN
           ec(k) = 0.d0      ; v1c(k,1) = 0.d0
           v2c(1,k,1) = 0.d0 ; v3c(k,1) = 0.d0
           CYCLE
         ENDIF  
         ec(k) = ec_lxc(k) * rho_lxc(k) 
         v1c(k,1) = vc_rho(k)
         v2c(1,k,1) = vc_sigma(k) * 2.0_DP
         v3c(k,1) = vc_tau(k)
       ENDDO
    ELSE
       !$acc parallel loop
       DO k = 1, length
          rh   = rho_lxc(2*k-1) + rho_lxc(2*k)
          atau = ABS(tau_lxc(2*k-1) + tau_lxc(2*k))
          ggrho2 = (sigma(3*k-2) + sigma(3*k))*4.0_DP
          IF ( rh <= rho_threshold_mgga   .OR. &
           ggrho2 <= grho2_threshold_mgga .OR. &
             atau <= tau_threshold_mgga  ) THEN
            ec(k) = 0.d0    
            v1c(k,1) = 0.d0 ; v3c(k,1) = 0.d0
            v1c(k,2) = 0.d0 ; v3c(k,2) = 0.d0
            DO ipol = 1, 3
              v2c(ipol,k,1) = 0.d0
              v2c(ipol,k,2) = 0.d0
            ENDDO
            CYCLE
          ENDIF   
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
    !$acc end data
    !
    DEALLOCATE( ec_lxc, vc_rho, vc_sigma, vc_tau )
  ENDIF
  !
  IF ( ANY(is_libxc(5:6)) ) DEALLOCATE( vlapl_rho )
  !
  !$acc end data
  DEALLOCATE( rho_lxc, sigma, tau_lxc, lapl_rho )
  !
#endif
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE xc_metagcx_
