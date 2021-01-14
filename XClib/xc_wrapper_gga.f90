!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
  !-------------------------------------------------------------------------
  !! Wrapper routine. Calls xc_gga-driver internal routines or the external
  !! ones from libxc, depending on the input choice.
  !
  !! NOTE: differently from 'xc_lda_drivers', here the input rho is in (up,down)
  !!       form (in the LSDA case).
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_par_mod,   ONLY: xc_func, xc_info 
#endif
  !
  USE kind_l,        ONLY: DP
  USE dft_par_mod,   ONLY: igcx, igcc, is_libxc, rho_threshold_gga, &
                           grho_threshold_gga
  USE qe_drivers_gga
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: ns
  !! spin dimension for input
  REAL(DP), INTENT(IN) :: rho(:,:)
  !! Charge density
  REAL(DP), INTENT(IN) :: grho(:,:,:)
  !! gradient
  REAL(DP), INTENT(OUT) :: ex(:)
  !! exchange energy
  REAL(DP), INTENT(OUT) :: ec(:)
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1x(:,:)
  !! exchange potential (density part)
  REAL(DP), INTENT(OUT) :: v2x(:,:)
  !! exchange potential (gradient part)
  REAL(DP), INTENT(OUT) :: v1c(:,:)
  !! correlation potential (density part)
  REAL(DP), INTENT(OUT) :: v2c(:,:)
  !! correlation potential (gradient part)
  REAL(DP), INTENT(OUT), OPTIONAL :: v2c_ud(:)
  !! correlation potential, cross term
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  REAL(DP), ALLOCATABLE :: rho_lxc(:), sigma(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vx_sigma(:)
  REAL(DP), ALLOCATABLE :: vc_rho(:), vc_sigma(:)
  !
  INTEGER :: fkind_x, np
  REAL(DP) :: rs, rtot, zet, vc_2(2)
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
  REAL(DP), ALLOCATABLE :: arho(:,:)
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grho2(:,:), grho_ud(:)
  !
  INTEGER :: k, is
  REAL(DP) :: sgn(2)
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  !
  IF (ns==2 .AND. .NOT. PRESENT(v2c_ud)) CALL xclib_infomsg( 'xc_gcx', 'WARNING: cross &
                &term v2c_ud not found xc_gcx (gga) call with polarized case' )
  !
  ex = 0.0_DP ;  v1x = 0.0_DP ;  v2x = 0.0_DP
  ec = 0.0_DP ;  v1c = 0.0_DP ;  v2c = 0.0_DP
  IF ( PRESENT(v2c_ud) ) v2c_ud = 0.0_DP
  !
  !
#if defined(__LIBXC)
  !
  fkind_x = -1
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
  ALLOCATE( rho_lxc(length*ns) )
  ALLOCATE( sigma(length*np) )
  !
  ALLOCATE( ex_lxc(length)    , ec_lxc(length)      )
  ALLOCATE( vx_rho(length*ns) , vx_sigma(length*np) )
  ALLOCATE( vc_rho(length*ns) , vc_sigma(length*np) )
  !
  !
  IF ( ns == 1 ) THEN
    !
    DO k = 1, length
      rho_lxc(k) = ABS( rho(k,1) )
      IF ( rho_lxc(k) > rho_threshold_gga ) &
         sigma(k) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
    ENDDO
    !
  ELSE
    !
    DO k = 1, length
       rho_lxc(2*k-1) = rho(k,1)
       rho_lxc(2*k)   = rho(k,2)
       !
       sigma(3*k-2) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
       sigma(3*k-1) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                      grho(3,k,1) * grho(3,k,2)
       sigma(3*k)   = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
    ENDDO
    !
  ENDIF
  !
  IF ( ns==1 .AND. ANY(.NOT.is_libxc(3:4)) ) THEN
     !
     CALL gcxc( length, ABS(rho(:,1)), sigma, ex, ec, v1x(:,1), v2x(:,1), v1c(:,1), v2c(:,1) )  
     !
     DO k = 1, length
        sgn(1) = SIGN(1._DP, rho(k,1))
        ex(k) = ex(k) * sgn(1)
        ec(k) = ec(k) * sgn(1)
     ENDDO
     !
  ENDIF
  !
  ! ---- GGA CORRELATION
  !
  IF ( is_libxc(4) ) THEN  !lda part of LYP not present in libxc (still so? - check)
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(4), rho_threshold_gga )
    fkind_x  = xc_f03_func_info_get_kind( xc_info(4) )
    CALL xc_f03_gga_exc_vxc( xc_func(4), lengthxc, rho_lxc(1), sigma(1), ec_lxc(1), vc_rho(1), vc_sigma(1) )
    !
    IF (.NOT. POLARIZED) THEN
      DO k = 1, length
        ec(k) = ec_lxc(k) * rho_lxc(k) * SIGN(1.0_DP, rho(k,1))
        v1c(k,1) = vc_rho(k)
        v2c(k,1) = vc_sigma(k)*2.d0
      ENDDO
    ELSE
      DO k = 1, length
        sgn(:) = 1.d0
        IF (rho_lxc(2*k-1)<rho_threshold_gga .OR. SQRT(ABS(sigma(3*k-2)))<grho_threshold_gga) sgn(1)=0.d0
        IF (rho_lxc(2*k)  <rho_threshold_gga .OR. SQRT(ABS(sigma(3*k)))  <grho_threshold_gga) sgn(2)=0.d0
        ec(k) = ec_lxc(k) * (rho_lxc(2*k-1)*sgn(1)+rho_lxc(2*k)*sgn(2))
        v1c(k,1) = vc_rho(2*k-1) * sgn(1)
        v1c(k,2) = vc_rho(2*k) * sgn(2)
        v2c(k,1) = vc_sigma(3*k-2)*2.d0 * sgn(1)
        v2c_ud(k)= vc_sigma(3*k-1) * sgn(1)*sgn(2)
        v2c(k,2) = vc_sigma(3*k)*2.d0 * sgn(2)
      ENDDO
    ENDIF
    !  
  ELSEIF ( (.NOT.is_libxc(4)) .AND. fkind_x/=XC_EXCHANGE_CORRELATION ) THEN
    !
    ALLOCATE( arho(length,ns), grho2(length,ns) )
    !
    IF ( ns /= 1 ) THEN
       !
       DO is = 1, 2
          grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
       ENDDO
       !
       IF (igcc==3 .OR. igcc==7 .OR. igcc==13 ) THEN
          !
          ALLOCATE( grho_ud(length) )
          !
          grho_ud = grho(1,:,1) * grho(1,:,2) + grho(2,:,1) * grho(2,:,2) + &
                    grho(3,:,1) * grho(3,:,2)
          !
          arho = rho
          !
          WHERE ( rho(:,1)+rho(:,2) < rho_threshold_gga )
             arho(:,1) = 0.0_DP
             arho(:,2) = 0.0_DP
          ENDWHERE
          !
          CALL gcc_spin_more( length, arho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
          !
          DEALLOCATE( grho_ud )
          !
       ELSE
          !
          ALLOCATE( rh(length), zeta(length) )
          !
          rh = rho(:,1) + rho(:,2)
          !
          zeta = 2.0_DP ! trash value, gcc-routines get rid of it when present
          WHERE ( rh > rho_threshold_gga ) zeta = ( rho(:,1) - rho(:,2) ) / rh(:)
          !
          grho2(:,1) = ( grho(1,:,1) + grho(1,:,2) )**2 + &
                       ( grho(2,:,1) + grho(2,:,2) )**2 + &
                       ( grho(3,:,1) + grho(3,:,2) )**2
          !
          CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
          !
          v2c(:,2)  = v2c(:,1)
          IF ( PRESENT(v2c_ud) ) v2c_ud(:) = v2c(:,1)
          !
          DEALLOCATE( rh, zeta )
          !
       ENDIF
       !   
    ENDIF
    !
    DEALLOCATE( arho, grho2 )
    !
  ENDIF  
  !
  ! --- GGA EXCHANGE
  !
  IF ( is_libxc(3) ) THEN
    !
    CALL xc_f03_func_set_dens_threshold( xc_func(3), rho_threshold_gga )
    CALL xc_f03_gga_exc_vxc( xc_func(3), lengthxc, rho_lxc(1), sigma(1), ex_lxc(1), vx_rho(1), vx_sigma(1) )
    !
    IF (.NOT. POLARIZED) THEN
      DO k = 1, length
        ex(k) = ex_lxc(k) * rho_lxc(k) * SIGN(1.0_DP, rho(k,1))
        v1x(k,1) = vx_rho(k)
        v2x(k,1) = vx_sigma(k)*2.d0
      ENDDO
    ELSE
      DO k = 1, length
        ex(k) = ex_lxc(k) * (rho_lxc(2*k-1)+rho_lxc(2*k))
        v1x(k,1) = vx_rho(2*k-1)
        v1x(k,2) = vx_rho(2*k)
        v2x(k,1) = vx_sigma(3*k-2)*2.d0
        v2x(k,2) = vx_sigma(3*k)*2.d0
      ENDDO
    ENDIF
    !
  ELSE
    !
    ALLOCATE( grho2(length,ns) )
    !
    IF ( ns /= 1 ) THEN
       !
       DO is = 1, 2
          grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
       ENDDO
       !
       CALL gcx_spin( length, rho, grho2, ex, v1x, v2x )
       !
    ENDIF
    !
    DEALLOCATE( grho2 )
    !
  ENDIF
  !
  DEALLOCATE( rho_lxc, sigma )
  DEALLOCATE( ex_lxc , ec_lxc   )
  DEALLOCATE( vx_rho , vx_sigma )
  DEALLOCATE( vc_rho , vc_sigma )
  !
#else
  !
  ALLOCATE( arho(length,ns), grho2(length,ns) )
  arho  = 0.0_DP
  grho2 = 0.0_DP
  !
  IF ( ns == 1 ) THEN
     !
     ! ... This is the spin-unpolarised case
     DO k = 1, length
        IF ( ABS(rho(k,1)) > rho_threshold_gga ) &
          grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
     ENDDO
     !
     !
     CALL gcxc( length, ABS(rho(:,1)), grho2(:,1), ex, ec, v1x(:,1), v2x(:,1), v1c(:,1), v2c(:,1) )
     !
     DO k = 1, length
        sgn(1) = SIGN(1._DP, rho(k,1))
        ex(k) = ex(k) * sgn(1)
        ec(k) = ec(k) * sgn(1)
     ENDDO
     !
  ELSE
     !
     DO is = 1, 2
        grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
     ENDDO
     !
     CALL gcx_spin( length, rho, grho2, ex, v1x, v2x )
     !
     IF (igcc==3 .OR. igcc==7 .OR. igcc==13 ) THEN
        !
        ALLOCATE( grho_ud(length) )
        !
        grho_ud = grho(1,:,1) * grho(1,:,2) + grho(2,:,1) * grho(2,:,2) + &
                  grho(3,:,1) * grho(3,:,2)
        !
        arho = rho
        !
        CALL gcc_spin_more( length, arho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
        !
        DEALLOCATE( grho_ud )
        !
     ELSE
        !
        ALLOCATE( rh(length), zeta(length) )
        !
        rh = rho(:,1) + rho(:,2)
        !
        zeta = 2.0_DP ! trash value, gcc-routines get rid of it when present
        WHERE ( rh > rho_threshold_gga ) zeta = ( rho(:,1) - rho(:,2) ) / rh(:)
        !
        grho2(:,1) = ( grho(1,:,1) + grho(1,:,2) )**2 + &
                     ( grho(2,:,1) + grho(2,:,2) )**2 + &
                     ( grho(3,:,1) + grho(3,:,2) )**2
        !
        CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
        !
        v2c(:,2)  = v2c(:,1)
        IF ( PRESENT(v2c_ud) ) v2c_ud(:) = v2c(:,1)
        !
        DEALLOCATE( rh, zeta )
        !
     ENDIF
     !   
  ENDIF
  !
  DEALLOCATE( arho, grho2 )
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE xc_gcx_
