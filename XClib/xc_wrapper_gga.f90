!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud, &
                   gpu_args_ )
  !-------------------------------------------------------------------------
  !! Wrapper to gpu or non gpu version of \(\texttt{xc_gcx}\).
  !
  USE kind_l,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: ns
  !! spin dimension for input
  REAL(DP), INTENT(IN) :: rho(length,ns)
  !! Charge density
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  !! gradient
  REAL(DP), INTENT(OUT) :: ex(length)
  !! exchange energy
  REAL(DP), INTENT(OUT) :: ec(length)
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  !! exchange potential (density part)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  !! exchange potential (gradient part)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  !! correlation potential (density part)
  REAL(DP), INTENT(OUT) :: v2c(length,ns)
  !! correlation potential (gradient part)
  REAL(DP), INTENT(OUT), OPTIONAL :: v2c_ud(length)
  !! correlation potential, cross term
  LOGICAL, INTENT(IN), OPTIONAL :: gpu_args_
  !! whether you wish to run on gpu in case use_gpu is true
  !
  LOGICAL :: gpu_args
  REAL(DP), ALLOCATABLE :: v2c_dummy(:)
  !
  gpu_args = .FALSE.
  !
  IF ( PRESENT(gpu_args_) ) gpu_args = gpu_args_
  !
  IF (ns==2 .AND. .NOT. PRESENT(v2c_ud)) CALL xclib_infomsg( 'xc_gcx', 'WARNING: cross &
                &term v2c_ud not found xc_gcx (gga) call with polarized case' )
  !
  IF ( gpu_args ) THEN
    !
    !$acc data present( rho, grho, ex, ec, v1x, v2x, v1c, v2c )
    IF (PRESENT(v2c_ud)) THEN
      !$acc data present( v2c_ud )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
      !$acc end data
    ELSE
      ALLOCATE( v2c_dummy(length) )
      !$acc data create( v2c_dummy )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_dummy )
      !$acc end data
      DEALLOCATE( v2c_dummy )
    ENDIF
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( rho, grho ), copyout( ex, ec, v1x, v2x, v1c, v2c )
    IF (PRESENT(v2c_ud)) THEN
      !$acc data copyout( v2c_ud )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
      !$acc end data
    ELSE
      ALLOCATE( v2c_dummy(length) )
      !$acc data create( v2c_dummy )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_dummy )
      !$acc end data
      DEALLOCATE( v2c_dummy )
    ENDIF
    !$acc end data
    !
  ENDIF  
  !
  RETURN
  !
END SUBROUTINE
!
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
  !-------------------------------------------------------------------------
  !! GGA wrapper routine - gpu double.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_setting_params,   ONLY: xc_func, xc_info, libxc_flags
#endif
  !
  USE kind_l,               ONLY: DP
  USE xclib_utils_and_para, ONLY: error_msg, nowarning
  USE dft_setting_params,   ONLY: igcx, igcc, is_libxc, rho_threshold_gga, &
                                  grho_threshold_gga, rho_threshold_lda,   &
                                  ishybrid, exx_started, exx_fraction
  USE qe_drivers_gga
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  INTEGER,  INTENT(IN) :: ns
  REAL(DP), INTENT(IN) :: rho(length,ns)
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  REAL(DP), INTENT(OUT) :: ex(length)
  REAL(DP), INTENT(OUT) :: ec(length)
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  REAL(DP), INTENT(OUT) :: v2c(length,ns)
  REAL(DP), INTENT(OUT) :: v2c_ud(length)
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  REAL(DP), ALLOCATABLE :: rho_lxc(:), sigma(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vx_sigma(:)
  REAL(DP), ALLOCATABLE :: vc_rho(:), vc_sigma(:)
  !
  INTEGER :: np
  REAL(DP) :: rs, rtot, zet, vc_2(2), arho_k, xcoef
  REAL(DP), PARAMETER :: pi34 = 0.6203504908994_DP
  !
  LOGICAL :: POLARIZED
  INTEGER :: ildax, ildac, pol_unpol
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grho2(:,:), grho_ud(:)
  !
  LOGICAL :: fkind_is_XC
  INTEGER :: k, is, ierr, fkind_x
  REAL(DP) :: rho_up, rho_dw, grho_up, grho_dw, sgn1
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  !$acc data present( rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
  !
  ierr = 0
  fkind_x = -1
  fkind_is_XC = .FALSE.
  !
#if defined(__LIBXC)
  !
  lengthxc = length
  !
  POLARIZED = .FALSE.
  IF (ns == 2) THEN
     POLARIZED = .TRUE.
  ENDIF
  !
  pol_unpol = 1
  np = 1
  IF ( ns == 2 ) THEN
     pol_unpol = 2
     np = 3
  ENDIF
  !
  IF (ANY(is_libxc(3:4))) THEN
    !
    ALLOCATE( rho_lxc(length*ns), sigma(length*np) )
    !$acc enter data create( rho_lxc, sigma )
    !
    IF ( is_libxc(3) ) THEN
      ALLOCATE( ex_lxc(length) )
      ALLOCATE( vx_rho(length*ns), vx_sigma(length*np) )
    ENDIF
    IF ( is_libxc(4) ) THEN
      ALLOCATE( ec_lxc(length) )
      ALLOCATE( vc_rho(length*ns), vc_sigma(length*np) )
    ENDIF
    !
    IF ( ns == 1 ) THEN
      !$acc parallel loop
      DO k = 1, length
        arho_k = ABS(rho(k,1))
        rho_lxc(k) = arho_k
        sigma(k) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
      ENDDO
    ELSE
      !$acc parallel loop
      DO k = 1, length
        rho_lxc(2*k-1) = rho(k,1)
        rho_lxc(2*k)   = rho(k,2)
        !
        sigma(3*k-2) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
        sigma(3*k-1) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                       grho(3,k,1) * grho(3,k,2)
        sigma(3*k)   = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
      ENDDO
    ENDIF
    !$acc update self( rho_lxc, sigma )
    !
  ENDIF
  !
#endif
  !
  IF (ANY(.NOT.is_libxc(3:4))) THEN
    !
    ALLOCATE( rh(length), grho2(length,ns) )
    !$acc enter data create( rh, grho2 )
    !$acc parallel loop
    DO k = 1, length
       rh(k) = ABS(rho(k,1))
       grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
    ENDDO
    !
    IF ( ns==1 ) THEN
      CALL gcxc( length, rh, grho2(:,1), ex, ec, v1x(:,1), v2x(:,1), v1c(:,1), &
                 v2c(:,1), ierr )
      !
      !$acc parallel loop
      DO k = 1, length
         sgn1 = SIGN(1._DP, rho(k,1))
         ex(k) = ex(k) * sgn1
         ec(k) = ec(k) * sgn1
      ENDDO
    ENDIF
    !
  ENDIF
  !
  ! ---- GGA CORRELATION
  !
  IF ( is_libxc(4) ) THEN  !lda part of LYP not present in libxc (still so? - check)
    !
#if defined(__LIBXC)
    !
    fkind_x = xc_f03_func_info_get_kind( xc_info(4) )
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(4), small )!rho_threshold_gga )
    IF (libxc_flags(4,0)==1) THEN
      CALL xc_f03_gga_exc_vxc( xc_func(4), lengthxc, rho_lxc(1), sigma(1), &
                               ec_lxc(1), vc_rho(1), vc_sigma(1) )
    ELSE
      CALL xc_f03_gga_vxc( xc_func(4), lengthxc, rho_lxc(1), sigma(1), &
                           vc_rho(1), vc_sigma(1) )
      ec_lxc = 0.d0
    ENDIF
    !
    !$acc data copyin( ec_lxc, vc_rho, vc_sigma )
    IF (.NOT. POLARIZED) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( rho_lxc(k) <= rho_threshold_lda ) THEN
          ec(k)=0.d0 ;  v1c(k,1)=0.d0 ;  v2c(k,1)=0.d0
          CYCLE
        ENDIF  
        ec(k) = ec_lxc(k) * rho_lxc(k) * SIGN(1.0_DP, rho(k,1))
        v1c(k,1) = vc_rho(k)
        IF ( rho_lxc(k) <= rho_threshold_gga .OR. &
             SQRT(ABS(sigma(k))) <= grho_threshold_gga) THEN
           v2c(k,1) = 0.d0
           CYCLE
        ENDIF
        v2c(k,1) = vc_sigma(k)*2.d0
      ENDDO
    ELSE
      !$acc parallel loop
      DO k = 1, length
        rho_up = rho_lxc(2*k-1)
        rho_dw = rho_lxc(2*k)
        grho_up = SQRT(ABS(sigma(3*k-2)))
        grho_dw = SQRT(ABS(sigma(3*k)))
        IF ( rho_up <= rho_threshold_lda .OR. rho_dw <= rho_threshold_lda ) THEN
          ec(k) = 0.d0    ;  v1c(k,1) = 0.d0 ;  v1c(k,2) = 0.d0
          v2c(k,1) = 0.d0 ;  v2c_ud(k)= 0.d0 ;  v2c(k,2) = 0.d0
          CYCLE
        ENDIF
        ec(k) = ec_lxc(k) * (rho_up+rho_dw)
        v1c(k,1) = vc_rho(2*k-1)
        v1c(k,2) = vc_rho(2*k)
        IF ( rho_up <= rho_threshold_gga .OR. rho_dw <= rho_threshold_gga .OR. &
             grho_up<=grho_threshold_gga .OR. grho_dw<=grho_threshold_gga ) THEN
          v2c(k,1) = 0.d0 ; v2c_ud(k)= 0.d0 ; v2c(k,2) = 0.d0
          CYCLE
        ENDIF
        v2c(k,1) = vc_sigma(3*k-2)*2.d0
        v2c_ud(k)= vc_sigma(3*k-1)
        v2c(k,2) = vc_sigma(3*k)*2.d0
      ENDDO
    ENDIF
    !
    !$acc end data
    DEALLOCATE( ec_lxc, vc_rho, vc_sigma )
    !
    fkind_is_XC = (fkind_x==XC_EXCHANGE_CORRELATION)
    !
#endif
    !
  ELSEIF ( (.NOT.is_libxc(4)) .AND. (.NOT.fkind_is_XC) ) THEN
    !
    IF ( ns /= 1 ) THEN
       !
       IF (igcc==3 .OR. igcc==7 .OR. igcc==13 ) THEN
          !
          ALLOCATE( grho_ud(length) )
          !$acc data create(grho_ud)
          !
          !$acc parallel loop
          DO k = 1, length
            grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
            grho_ud(k) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                         grho(3,k,1) * grho(3,k,2)
            grho2(k,2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
          ENDDO
          !
          CALL gcc_spin_more( length, rho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
          !
          !$acc end data
          DEALLOCATE( grho_ud )
          !
       ELSE
          !
          ALLOCATE( zeta(length) )
          !$acc data create( zeta )
          !
          !$acc parallel loop
          DO k = 1, length
            rh(k) = rho(k,1) + rho(k,2)
            IF ( rh(k) > rho_threshold_gga ) THEN
              zeta(k) = (rho(k,1)-rho(k,2)) / rh(k)
            ELSE
              zeta(k) = 2.0_DP ! trash value, gcc-routines get rid of it when present
            ENDIF
            grho2(k,1) = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                         ( grho(2,k,1) + grho(2,k,2) )**2 + &
                         ( grho(3,k,1) + grho(3,k,2) )**2
            grho2(k,2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
          ENDDO
          !
          CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
          !
          !$acc parallel loop
          DO k = 1, length
            v2c(k,2) = v2c(k,1)
            IF ( ns==2 ) v2c_ud(k) = v2c(k,1)
          ENDDO
          !
          !$acc end data
          DEALLOCATE( zeta )
          !
       ENDIF
       !
    ENDIF
    !
  ENDIF
  !
  ! --- GGA EXCHANGE
  !
  IF ( is_libxc(3) ) THEN
    !
#if defined(__LIBXC)
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(3), grho_threshold_gga )
    IF (libxc_flags(3,0)==1) THEN
      CALL xc_f03_gga_exc_vxc( xc_func(3), lengthxc, rho_lxc(1), sigma(1), ex_lxc(1), vx_rho(1), vx_sigma(1) )
    ELSE
      CALL xc_f03_gga_vxc( xc_func(3), lengthxc, rho_lxc(1), sigma(1), vx_rho(1), vx_sigma(1) )
      ex_lxc = 0.d0
    ENDIF
    !
    xcoef = 1.d0
    IF ( ishybrid .AND. exx_started .AND. exx_fraction>0.d0) xcoef = 1.d0-exx_fraction
    !
    !$acc data copyin( ex_lxc, vx_rho, vx_sigma )
    IF (.NOT. POLARIZED) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( rho_lxc(k) <= rho_threshold_lda ) THEN
          ex(k) = 0.d0 ;  v1x(k,1) = 0.d0 ;  v2x(k,1) = 0.d0
          CYCLE
        ENDIF
        ex(k) = xcoef * ex_lxc(k) * rho_lxc(k) * SIGN(1.0_DP, rho(k,1))
        v1x(k,1) = xcoef * vx_rho(k)
        IF ( rho_lxc(k) <= rho_threshold_gga .OR. &
             SQRT(ABS(sigma(k))) <= grho_threshold_gga) THEN
          v2x(k,1) = 0.d0
          CYCLE
        ENDIF
        v2x(k,1) = xcoef * vx_sigma(k)*2.d0
      ENDDO
    ELSE
      !$acc parallel loop
      DO k = 1, length
        rho_up = rho_lxc(2*k-1)
        rho_dw = rho_lxc(2*k)
        grho_up = SQRT(ABS(sigma(3*k-2)))
        grho_dw = SQRT(ABS(sigma(3*k)))
        IF ( rho_up <= rho_threshold_lda .OR. rho_dw <= rho_threshold_lda ) THEN
          ex(k) = 0.d0
          v1x(k,1) = 0.d0 ;  v1x(k,2) = 0.d0
          v2x(k,1) = 0.d0 ;  v2x(k,2) = 0.d0
          CYCLE
        ENDIF
        ex(k) = xcoef * ex_lxc(k) * (rho_up+rho_dw)
        v1x(k,1) = xcoef * vx_rho(2*k-1)
        v1x(k,2) = xcoef * vx_rho(2*k)
        IF ( rho_up <= rho_threshold_gga .OR. rho_dw <= rho_threshold_gga .OR. &
             grho_up<=grho_threshold_gga .OR. grho_dw<=grho_threshold_gga ) THEN
          v2x(k,1) = 0.d0 ;  v2x(k,2) = 0.d0
          CYCLE
        ENDIF
        v2x(k,1) = xcoef * vx_sigma(3*k-2)*2.d0
        v2x(k,2) = xcoef * vx_sigma(3*k)*2.d0
      ENDDO
    ENDIF
    !
    !$acc end data
    DEALLOCATE( ex_lxc, vx_rho, vx_sigma )
    !
#endif
    !
  ELSE
    !
    IF ( ns > 1 ) THEN
      !$acc parallel loop collapse(2)
      DO is = 1, ns
        DO k = 1, length
          grho2(k,is) = grho(1,k,is)**2 + grho(2,k,is)**2 + grho(3,k,is)**2
        ENDDO
      ENDDO
      !
      CALL gcx_spin( length, rho, grho2, ex, v1x, v2x, ierr )
    ENDIF
    !
  ENDIF
  !
  IF (ANY(.NOT.is_libxc(3:4))) THEN
    !$acc exit data delete( rh, grho2 )
    DEALLOCATE( rh, grho2 )
  ENDIF
#if defined(__LIBXC)
  IF (ANY(is_libxc(3:4))) THEN
    !$acc exit data delete(rho_lxc,sigma)
    DEALLOCATE( rho_lxc, sigma )
  ENDIF
#endif
  !
  !$acc end data
  !
  IF (ierr/=0 .AND. .NOT.nowarning) CALL xclib_error( 'xc_gcx_', error_msg(ierr), 1 )
  !
  RETURN
  !
END SUBROUTINE xc_gcx_
