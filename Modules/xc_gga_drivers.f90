MODULE xc_gga
!
USE kinds,     ONLY: DP
USE funct,     ONLY: get_igcx, get_igcc, is_libxc,    &
                     exx_is_active, get_exx_fraction, &
                     get_screening_parameter, get_gau_parameter
!
IMPLICIT NONE
!
PRIVATE
SAVE
!
!  GGA exchange-correlation drivers
PUBLIC :: xc_gcx, gcxc, gcx_spin, gcc_spin, gcc_spin_more, &
          change_threshold_gga
!
!  input thresholds (default values)
REAL(DP) :: rho_threshold = 1.D-6
REAL(DP) :: grho_threshold = 1.D-10
!
!
 CONTAINS
!
!
!-----------------------------------------------------------------------
SUBROUTINE change_threshold_gga( rho_thr_in, grho_thr_in )
  !--------------------------------------------------------------------
  !! Change rho and grho thresholds.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho_thr_in
  REAL(DP), INTENT(IN), OPTIONAL :: grho_thr_in
  !
  rho_threshold = rho_thr_in
  IF (PRESENT(grho_thr_in)) grho_threshold = grho_thr_in
  !
  RETURN
  !
END SUBROUTINE
!
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
  !-------------------------------------------------------------------------
  !! Wrapper routine. Calls xc_gga-driver routines from internal libraries
  !! of q-e or from the external libxc, depending on the input choice.
  !
  !! NOTE: differently from 'xc_lda_drivers', here the input rho is in (up,down)
  !!       form (in the LSDA case).
  !
  !! NOTE: look at 'PP/src/benchmark_libxc.f90' to test and see the differences
  !!       between q-e and libxc libraries.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
#endif
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
  !! correlation (gradient part)
  REAL(DP), INTENT(OUT), OPTIONAL :: v2c_ud(:)
  !! correlation
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info1, xc_info2
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
  INTEGER :: igcx, igcc
  INTEGER :: k, is
  REAL(DP) :: sgn(2)
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  !
  IF (ns==2 .AND. .NOT. PRESENT(v2c_ud)) CALL errore( 'xc_gga', 'cross &
                                             &term v2c_ud not found', 1 )
  !
  igcx = get_igcx()
  igcc = get_igcc()
  !
  ex = 0.0_DP ;  v1x = 0.0_DP ;  v2x = 0.0_DP
  ec = 0.0_DP ;  v1c = 0.0_DP ;  v2c = 0.0_DP
  IF ( PRESENT(v2c_ud) ) v2c_ud = 0.0_DP
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
      IF ( rho_lxc(k) > rho_threshold ) &
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
  ! --- GGA EXCHANGE
  !
  IF ( is_libxc(3) ) THEN
    !
    CALL xc_f03_func_init( xc_func, igcx, pol_unpol )
     xc_info1 = xc_f03_func_get_info( xc_func )
     CALL xc_f03_func_set_dens_threshold( xc_func, rho_threshold )
     fkind_x  = xc_f03_func_info_get_kind( xc_info1 )
     CALL xc_f03_gga_exc_vxc( xc_func, lengthxc, rho_lxc(1), sigma(1), ex_lxc(1), vx_rho(1), vx_sigma(1) )
    CALL xc_f03_func_end( xc_func )
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
  ! ---- GGA CORRELATION
  !
  IF ( is_libxc(4) ) THEN  !lda part of LYP not present in libxc
    !
    CALL xc_f03_func_init( xc_func, igcc, pol_unpol )
     xc_info2 = xc_f03_func_get_info( xc_func )
     CALL xc_f03_func_set_dens_threshold( xc_func, rho_threshold )
     CALL xc_f03_gga_exc_vxc( xc_func, lengthxc, rho_lxc(1), sigma(1), ec_lxc(1), vc_rho(1), vc_sigma(1) )
    CALL xc_f03_func_end( xc_func )
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
        IF (rho_lxc(2*k-1)<rho_threshold .OR. SQRT(ABS(sigma(3*k-2)))<grho_threshold) sgn(1)=0.d0
        IF (rho_lxc(2*k)  <rho_threshold .OR. SQRT(ABS(sigma(3*k)))  <grho_threshold) sgn(2)=0.d0
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
          WHERE ( rho(:,1)+rho(:,2) < rho_threshold )
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
          WHERE ( rh > rho_threshold ) zeta = ( rho(:,1) - rho(:,2) ) / rh(:)
          !
          grho2(:,1) = ( grho(1,:,1) + grho(1,:,2) )**2 + &
                       ( grho(2,:,1) + grho(2,:,2) )**2 + &
                       ( grho(3,:,1) + grho(3,:,2) )**2
          !
          CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
          !
          v2c(:,2)  = v2c(:,1)
          v2c_ud(:) = v2c(:,1)
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
        IF ( ABS(rho(k,1)) > rho_threshold ) &
          grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
     ENDDO
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
        WHERE ( rh > rho_threshold ) zeta = ( rho(:,1) - rho(:,2) ) / rh(:)
        !
        grho2(:,1) = ( grho(1,:,1) + grho(1,:,2) )**2 + &
                     ( grho(2,:,1) + grho(2,:,2) )**2 + &
                     ( grho(3,:,1) + grho(3,:,2) )**2
        !
        CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
        !
        v2c(:,2)  = v2c(:,1)
        v2c_ud(:) = v2c(:,1)
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
END SUBROUTINE xc_gcx
!
!
!-----------------------------------------------------------------------
!------- GRADIENT CORRECTION DRIVERS ----------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE gcxc( length, rho_in, grho_in, sx_out, sc_out, v1x_out, &
                                          v2x_out, v1c_out, v2c_out )
  !---------------------------------------------------------------------
  !! Gradient corrections for exchange and correlation - Hartree a.u. 
  !! See comments at the beginning of module for implemented cases
  !
  ! Input:  rho, grho=|\nabla rho|^2
  ! Definition:  E_x = \int E_x(rho,grho) dr
  ! Output: sx = E_x(rho,grho)
  !         v1x= D(E_x)/D(rho)
  !         v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !         sc, v1c, v2c as above for correlation
  !
  USE exch_gga
  USE corr_gga
  !
#if defined(use_beef)
  USE beef_interface, ONLY: beefx, beeflocalcorr
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in, grho_in
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sx_out, sc_out, v1x_out, &
                                              v2x_out, v1c_out, v2c_out
  !
  ! ... local variables
  !
  INTEGER :: ir, igcx, igcc
  REAL(DP) :: rho, grho
  REAL(DP) :: sx, v1x, v2x
  REAL(DP) :: sx_, v1x_, v2x_
  REAL(DP) :: sxsr, v1xsr, v2xsr
  REAL(DP) :: sc, v1c, v2c
  REAL(DP) :: screening_parameter, gau_parameter
  REAL(DP) :: exx_fraction
  LOGICAL  :: exx_started
  
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
  igcx = get_igcx()
  igcc = get_igcc()
  exx_started = exx_is_active()
  exx_fraction = get_exx_fraction()
  IF (igcx == 12) screening_parameter = get_screening_parameter()
  IF (igcx == 20) gau_parameter = get_gau_parameter()
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, grho, sx, sx_, sxsr, v1x, v1x_, v1xsr, &
!$omp             v2x, v2x_, v2xsr, sc, v1c, v2c )
  DO ir = 1, length  
     !
     grho = grho_in(ir)
     !
     IF ( rho_in(ir) <= rho_threshold .OR. grho <= grho_threshold ) THEN
        sx_out(ir)  = 0.0_DP ;   sc_out(ir)  = 0.0_DP
        v1x_out(ir) = 0.0_DP ;   v1c_out(ir) = 0.0_DP
        v2x_out(ir) = 0.0_DP ;   v2c_out(ir) = 0.0_DP
        CYCLE
     ENDIF
     !
     rho  = ABS(rho_in(ir))
     !
     ! ... EXCHANGE
     !  
     SELECT CASE( igcx )
     CASE( 1 )
        !
        CALL becke88( rho, grho, sx, v1x, v2x )
        !
     CASE( 2 )
        !
        CALL ggax( rho, grho, sx, v1x, v2x )
        !
     CASE( 3 )
        !
        CALL pbex( rho, grho, 1, sx, v1x, v2x )
        !
     CASE( 4 )
        !
        CALL pbex( rho, grho, 2, sx, v1x, v2x )
        !
     CASE( 5 )
        !
        IF (igcc == 5) CALL hcth( rho, grho, sx, v1x, v2x )
        !
     CASE( 6 )
        !
        CALL optx( rho, grho, sx, v1x, v2x )
        !
     ! case igcx == 7 (meta-GGA) must be treated in a separate call to another
     ! routine: needs kinetic energy density in addition to rho and grad rho
     CASE( 8 ) ! 'PBE0'
        !
        CALL pbex( rho, grho, 1, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 9 ) ! 'B3LYP'
        !
        CALL becke88( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = 0.72_DP * sx
           v1x = 0.72_DP * v1x
           v2x = 0.72_DP * v2x
        ENDIF
        !
     CASE( 10 ) ! 'pbesol'
        !
        CALL pbex( rho, grho, 3, sx, v1x, v2x )
        !
     CASE( 11 ) ! 'wc'
        !
        CALL wcx( rho, grho, sx, v1x, v2x )
        !
     CASE( 12 ) ! 'pbexsr'
        !
        CALL pbex( rho, grho, 1, sx, v1x, v2x )
        !
        IF (exx_started) THEN
          CALL pbexsr( rho, grho, sxsr, v1xsr, v2xsr, screening_parameter )
          sx  = sx  - exx_fraction * sxsr
          v1x = v1x - exx_fraction * v1xsr
          v2x = v2x - exx_fraction * v2xsr
        ENDIF
        !
     CASE( 13 ) ! 'rPW86'
        !
        CALL rPW86( rho, grho, sx, v1x, v2x )
        !
     CASE( 16 ) ! 'C09x'
        !
        CALL c09x( rho, grho, sx, v1x, v2x )
        !
     CASE( 17 ) ! 'sogga'
        !
        CALL sogga( rho, grho, sx, v1x, v2x )
        !
     CASE( 19 ) ! 'pbeq2d'
        !
        CALL pbex( rho, grho, 4, sx, v1x, v2x )
        !
     CASE( 20 ) ! 'gau-pbe'
        !
        CALL pbex( rho, grho, 1, sx, v1x, v2x )
        IF (exx_started) THEN
          CALL pbexgau( rho, grho, sxsr, v1xsr, v2xsr, gau_parameter )
          sx  = sx  - exx_fraction * sxsr
          v1x = v1x - exx_fraction * v1xsr
          v2x = v2x - exx_fraction * v2xsr
        ENDIF
        !
     CASE( 21 ) ! 'pw86'
        !
        CALL pw86( rho, grho, sx, v1x, v2x )
        !
     CASE( 22 ) ! 'b86b'
        !
        CALL becke86b( rho, grho, sx, v1x, v2x )
        ! CALL b86b( rho, grho, 1, sx, v1x, v2x )
        !
     CASE( 23 ) ! 'optB88'
        !
        CALL pbex( rho, grho, 5, sx, v1x, v2x )
        !
     CASE( 24 ) ! 'optB86b'
        !
        CALL pbex( rho, grho, 6, sx, v1x, v2x )
        ! CALL b86b (rho, grho, 2, sx, v1x, v2x)
        !
     CASE( 25 ) ! 'ev93'
        !
        CALL pbex( rho, grho, 7, sx, v1x, v2x )
        !
     CASE( 26 ) ! 'b86r'
        !
        CALL b86b( rho, grho, 3, sx, v1x, v2x )
        !
     CASE( 27 ) ! 'cx13'
        !
        CALL cx13( rho, grho, sx, v1x, v2x )
        !
     CASE( 28 ) ! 'X3LYP'
        !
        CALL becke88( rho, grho, sx, v1x, v2x )
        CALL pbex( rho, grho, 1, sx_, v1x_, v2x_ )
        IF (exx_started) THEN
           sx  = REAL(0.765*0.709,DP) * sx
           v1x = REAL(0.765*0.709,DP) * v1x
           v2x = REAL(0.765*0.709,DP) * v2x
           sx  = sx  + REAL(0.235*0.709,DP) * sx_
           v1x = v1x + REAL(0.235*0.709,DP) * v1x_
           v2x = v2x + REAL(0.235*0.709,DP) * v2x_
        ENDIF
        !
     CASE( 29, 31 ) ! 'cx0'or `cx0p'
        !
        CALL cx13( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 30 ) ! 'r860'
        !
        CALL rPW86( rho, grho, sx, v1x, v2x )
        !
        IF (exx_started) then
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 38 ) ! 'BR0'
        !
        CALL b86b( rho, grho, 3, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 40 ) ! 'c090'
        !
        CALL c09x( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 41 ) ! 'B86BPBEX'
        !
        CALL becke86b( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 42 ) ! 'BHANDHLYP'
        !
        CALL becke88( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
#ifdef use_beef
     CASE( 43 ) ! 'beefx'
        ! last parameter = 0 means do not add LDA (=Slater) exchange
        ! (espresso) will add it itself
        CALL beefx(rho, grho, sx, v1x, v2x, 0)
#endif
        !
     CASE( 44 ) ! 'RPBE'
        !
        CALL pbex( rho, grho, 8, sx, v1x, v2x )
        !
     CASE( 45 ) ! 'W31X'
        !
        CALL pbex( rho, grho, 9, sx, v1x, v2x )
        !
     CASE( 46 ) ! 'W32X'
        !
        CALL b86b( rho, grho, 4, sx, v1x, v2x )
        !
     CASE DEFAULT
        !
        sx  = 0.0_DP
        v1x = 0.0_DP
        v2x = 0.0_DP
        !
     END SELECT
     !
     !
     ! ... CORRELATION
     !
     SELECT CASE( igcc )
     CASE( 1 )
        !
        CALL perdew86( rho, grho, sc, v1c, v2c )
        !
     CASE( 2 )
        !
        CALL ggac( rho, grho, sc, v1c, v2c )
        !
     CASE( 3 )
        !
        CALL glyp( rho, grho, sc, v1c, v2c )
        !
     CASE( 4 )
        !
        CALL pbec( rho, grho, 1, sc, v1c, v2c )
        !
     ! igcc == 5 (HCTH) is calculated together with case igcx=5
     ! igcc == 6 (meta-GGA) is treated in a different routine
     CASE( 7 ) !'B3LYP'
        !
        CALL glyp( rho, grho, sc, v1c, v2c )
        IF (exx_started) THEN
           sc  = 0.81_DP * sc
           v1c = 0.81_DP * v1c
           v2c = 0.81_DP * v2c
        ENDIF
        !
     CASE( 8 ) ! 'PBEsol'
        !
        CALL pbec( rho, grho, 2, sc, v1c, v2c )
        !
     ! igcc ==  9 set to 5, back-compatibility
     ! igcc == 10 set to 6, back-compatibility
     ! igcc == 11 M06L calculated in another routine
     CASE( 12 ) ! 'Q2D'
        !
        CALL pbec( rho, grho, 3, sc, v1c, v2c )
        !
     CASE( 13 ) !'X3LYP'
        !
        CALL glyp( rho, grho, sc, v1c, v2c )
        IF (exx_started) THEN
           sc  = 0.871_DP * sc
           v1c = 0.871_DP * v1c
           v2c = 0.871_DP * v2c
        ENDIF
#ifdef use_beef
     CASE( 14 ) ! 'BEEF'
        ! last parameter 0 means: do not add lda contributions
        ! espresso will do that itself
        call beeflocalcorr(rho, grho, sc, v1c, v2c, 0)
#endif
        !
     CASE DEFAULT
        !
        sc = 0.0_DP
        v1c = 0.0_DP
        v2c = 0.0_DP
        !
     END SELECT
     !
     sx_out(ir)  = sx    ;  sc_out(ir)  = sc
     v1x_out(ir) = v1x   ;  v1c_out(ir) = v1c
     v2x_out(ir) = v2x   ;  v2c_out(ir) = v2c
     !
  ENDDO 
!$omp end do
!$omp end parallel
  !
  !
  RETURN
  !
END SUBROUTINE gcxc
!
!
!===============> SPIN <===============!
!
!-------------------------------------------------------------------------
SUBROUTINE gcx_spin( length, rho_in, grho2_in, sx_tot, v1x_out, v2x_out )
  !-----------------------------------------------------------------------
  !! Gradient corrections for exchange - Hartree a.u.
  !
  USE exch_gga
#if defined(use_beef)
  USE beef_interface, ONLY: beefx
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Length of the input/output arrays
  REAL(DP), INTENT(IN),  DIMENSION(length,2) :: rho_in
  !! Up and down charge density
  REAL(DP), INTENT(IN),  DIMENSION(length,2) :: grho2_in
  !! Up and down gradient of the charge
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sx_tot
  !! Energy exchange GGA
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1x_out
  !! Exchange potential (density part)
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v2x_out
  !! Exchange potantial (gradient part)
  !
  ! ... local variables
  !
  INTEGER :: ir, is, iflag, igcx, igcc
  REAL(DP) :: rho(2), grho2(2)
  REAL(DP) :: v1x(2), v2x(2)
  REAL(DP) :: sx(2), rnull(2)
  REAL(DP) :: sxsr(2)
  REAL(DP) :: v1xsr(2), v2xsr(2)
  REAL(DP) :: screening_parameter, gau_parameter
  REAL(DP) :: exx_fraction
  LOGICAL  :: exx_started
  !
  REAL(DP), PARAMETER :: small=1.D-10
  REAL(DP), PARAMETER :: rho_trash=0.5_DP, grho2_trash=0.2_DP
  ! temporary values assigned to rho and grho when they
  ! are too small in order to avoid numerical problems.
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
  sx_tot = 0.0_DP
  !
  igcx = get_igcx()
  igcc = get_igcc()
  exx_started = exx_is_active()
  exx_fraction = get_exx_fraction()
  IF (igcx == 12 .AND. exx_started) screening_parameter = get_screening_parameter()
  IF (igcx == 20 .AND. exx_started) gau_parameter = get_gau_parameter()
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, grho2, rnull, sx, sxsr, v1x, v1xsr, &
!$omp             v2x, v2xsr )
  DO ir = 1, length  
     !
     rho(:) = rho_in(ir,:)
     grho2(:) = grho2_in(ir,:)
     rnull(:) = 1.0_DP
     !
     IF ( rho(1)+rho(2) <= small ) THEN
        sx_tot(ir) = 0.0_DP
        v1x_out(ir,:) = 0.0_DP ; v2x_out(ir,:) = 0.0_DP
        CYCLE
     ELSE
        DO is = 1, 2
           IF ( rho(is)<=small .OR. SQRT(ABS(grho2(is)))<=small ) THEN
             rho(is) = rho_trash
             grho2(is) = grho2_trash
             rnull(is) = 0.0_DP
           ENDIF
        ENDDO
     ENDIF
     !
     !
     ! ... exchange
     !
     SELECT CASE( igcx )
     CASE( 0 )
        !
        sx_tot(ir) = 0.0_DP
        v1x = 0.0_DP
        v2x = 0.0_DP
        !
     CASE( 1 )
        !
        CALL becke88_spin( rho(1), rho(2), grho2(1), grho2(2), sx(1), sx(2), v1x(1), v1x(2), v2x(1), v2x(2) )
        !
        sx_tot(ir) = sx(1)*rnull(1) + sx(2)*rnull(2)
        !
     CASE( 2 )
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL ggax( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL ggax( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 3, 4, 8, 10, 12, 20, 23, 24, 25, 44, 45 )
        ! igcx=3:  PBE,  igcx=4:  revised PBE, igcx=8:  PBE0, igcx=10: PBEsol
        ! igcx=12: HSE,  igcx=20: gau-pbe,     igcx=23: obk8, igcx=24: ob86,
        ! igcx=25: ev93, igcx=44: RPBE,        igcx=45: W31X
        !
        iflag = 1
        IF ( igcx== 4 ) iflag = 2
        IF ( igcx==10 ) iflag = 3
        IF ( igcx==23 ) iflag = 5
        IF ( igcx==24 ) iflag = 6
        IF ( igcx==25 ) iflag = 7
        IF ( igcx==44 ) iflag = 8
        IF ( igcx==45 ) iflag = 9
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL pbex( rho(1), grho2(1), iflag, sx(1), v1x(1), v2x(1) )
        CALL pbex( rho(2), grho2(2), iflag, sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( igcx == 8 .AND. exx_started ) THEN
           !
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
           !
        ELSEIF ( igcx == 12 .AND. exx_started ) THEN
           !
           CALL pbexsr( rho(1), grho2(1), sxsr(1), v1xsr(1), &
                                          v2xsr(1), screening_parameter )
           CALL pbexsr( rho(2), grho2(2), sxsr(2), v1xsr(2), &
                                          v2xsr(2), screening_parameter )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction*0.5_DP * ( sxsr(1)*rnull(1) + &
                                                               sxsr(2)*rnull(2) )
           v1x = v1x - exx_fraction * v1xsr
           v2x = v2x - exx_fraction * v2xsr * 2.0_DP
           !
        ELSEIF ( igcx == 20 .AND. exx_started ) THEN
           ! gau-pbe
           !CALL pbexgau_lsd( rho, grho2, sxsr, v1xsr, v2xsr, gau_parameter_l )
           CALL pbexgau( rho(1), grho2(1), sxsr(1), v1xsr(1), v2xsr(1), gau_parameter )
           CALL pbexgau( rho(2), grho2(2), sxsr(2), v1xsr(2), v2xsr(2), gau_parameter )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction*0.5_DP * ( sxsr(1)*rnull(1) + &
                                                               sxsr(2)*rnull(2) )
           v1x = v1x - exx_fraction * v1xsr
           v2x = v2x - exx_fraction * v2xsr * 2.0_DP
           !
        ENDIF
        !
     CASE( 9 )                    ! B3LYP
        !
        CALL becke88_spin( rho(1), rho(2), grho2(1), grho2(2), sx(1), sx(2), v1x(1), v1x(2), v2x(1), v2x(2) )
        !
        sx_tot(ir) = sx(1)*rnull(1) + sx(2)*rnull(2)
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = 0.72_DP * sx_tot(ir)
           v1x = 0.72_DP * v1x
           v2x = 0.72_DP * v2x
        ENDIF
        !
     CASE( 11 )                   ! 'Wu-Cohen'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL wcx( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL wcx( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 13 )                   ! 'revised PW86 for vdw-df2'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL rPW86( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL rPW86( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 16 )                   ! 'c09x for vdw-df-c09.'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL c09x( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL c09x( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 21 )                   ! 'PW86'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL pw86( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL pw86( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 22 )                   ! 'B86B'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL becke86b( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL becke86b( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
      CASE( 26, 46 )                  ! 'B86R for rev-vdW-DF2'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        IF ( igcx==26 ) iflag = 3 ! B86R for rev-vdW-DF2
        IF ( igcx==46 ) iflag = 4 ! W32X for vdW-DF3-opt2
        CALL b86b( rho(1), grho2(1), iflag, sx(1), v1x(1), v2x(1) )
        CALL b86b( rho(2), grho2(2), iflag, sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 27 )                   ! 'cx13 for vdw-df-cx'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL cx13( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL cx13( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 28 )                   ! X3LYP
        !
        CALL becke88_spin( rho(1), rho(2), grho2(1), grho2(2), sx(1), sx(2), v1x(1), v1x(2), v2x(1), v2x(2) )
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL pbex( rho(1), grho2(1), 1, sxsr(1), v1xsr(1), v2xsr(1) )
        CALL pbex( rho(2), grho2(2), 1, sxsr(2), v1xsr(2), v2xsr(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sxsr(1)*rnull(1) + sxsr(2)*rnull(2) ) * 0.235_DP + &
                              (   sx(1)*rnull(1) +   sx(2)*rnull(2) ) * 0.765_DP
        v1x = v1xsr * 0.235_DP + v1x * 0.765_DP
        v2x = v2xsr * 0.235_DP * 2.0_DP + v2x * 0.765_DP
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = 0.709_DP * sx_tot(ir)
           v1x = 0.709_DP * v1x
           v2x = 0.709_DP * v2x
        ENDIF
        !
     CASE( 29, 31 )               ! 'cx0 for vdw-df-cx0' or `cx0p for vdW-DF-cx0p'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL cx13( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL cx13( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 30 )                   ! 'R860' = 'rPW86-0' for vdw-df2-0'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL rPW86( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL rPW86( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 38 )                  ! 'br0 for vdw-df2-BR0' etc
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL b86b( rho(1), grho2(1), 3, sx(1), v1x(1), v2x(1) )
        CALL b86b( rho(2), grho2(2), 3, sx(2), v1x(2), v2x(2) )     
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF  
        !
     CASE( 40 )                  ! 'c090 for vdw-df-c090' etc
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL c09x( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL c09x( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )  
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 41 )                 ! B86X for B86BPBEX hybrid
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL becke86b( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL becke86b( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started ) THEN
           sx_tot = (1.0_DP - exx_fraction) * sx_tot
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 42 )                ! B88X for BHANDHLYP
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL becke88( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL becke88( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot = 0.5_DP * ( sx(1)*rnull(1) + sx(2)*rnull(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started ) THEN
           sx_tot = (1.0_DP - exx_fraction) * sx_tot
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
#ifdef use_beef
     CASE( 43 ) ! 'beefx'
        IF (rho(1) > small .AND. SQRT (ABS (grho2(1)) ) > small) THEN
           call beefx(2.0_DP * rho(1), 4.0_DP * grho2(1), sx(1), v1x(1), v2x(1), 0)
        ELSE
           sx(1) = 0.0_DP
           v1x(1) = 0.0_DP
           v2x(1) = 0.0_DP
        ENDIF
        IF (rho(2) > small .AND. SQRT (ABS (grho2(2)) ) > small) THEN
           CALL beefx(2.0_DP * rho(2), 4.0_DP * grho2(2), sx(2), v1x(2), v2x(2), 0)
           CALL beefx(2.0_DP * rho(2), 4.0_DP * grho2(2), sx(2), v1x(2), v2x(2), 0)
        ELSE
           sx(2) = 0.0_DP
           v1x(2) = 0.0_DP
           v2x(2) = 0.0_DP
        ENDIF
        sx_tot(ir) = 0.5_DP * (sx(1) + sx(2))
        v2x  = 2.0_DP * v2x
#endif
     !
     ! case igcx == 5 (HCTH) and 6 (OPTX) not implemented
     ! case igcx == 7 (meta-GGA) must be treated in a separate call to another
     ! routine: needs kinetic energy density in addition to rho and grad rho
     !
     CASE DEFAULT
        !
        sx = 0.0_DP
        v1x = 0.0_DP
        v2x = 0.0_DP
        !
     END SELECT
     !
     v1x_out(ir,:) = v1x(:) * rnull(:)
     v2x_out(ir,:) = v2x(:) * rnull(:)
     !
  ENDDO
!$omp end do
!$omp end parallel
  !
  !
  RETURN
  !
END SUBROUTINE gcx_spin
!
!
!--------------------------------------------------------------------------------
SUBROUTINE gcc_spin( length, rho_in, zeta_io, grho_in, sc_out, v1c_out, v2c_out )
  !-------------------------------------------------------------------------------
  !! Gradient corrections for correlations - Hartree a.u.  
  !! Implemented: Perdew86, GGA (PW91), PBE
  !
  USE corr_gga
  !
#if defined(use_beef)
  USE beef_interface, ONLY: beeflocalcorrspin
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! the length of the I/O arrays
  REAL(DP), INTENT(IN), DIMENSION(length) :: rho_in
  !! the total charge
  REAL(DP), INTENT(INOUT), DIMENSION(length) :: zeta_io
  !! the magnetization
  REAL(DP), INTENT(IN), DIMENSION(length) :: grho_in
  !! the gradient of the charge squared
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sc_out
  !! correlation energies
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1c_out
  !! correlation potential (density part)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c_out
  !! correlation potential (gradient part)
  !
  ! ... local variables
  !
  INTEGER :: ir, igcc
  REAL(DP) :: rho, zeta, grho
  REAL(DP) :: sc, v1c(2), v2c
  !REAL(DP), PARAMETER :: small=1.E-10_DP !, epsr=1.E-6_DP
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
  igcc = get_igcc()
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, zeta, grho, sc, v1c, v2c )
  DO ir = 1, length
    !
    rho  = rho_in(ir)
    grho = grho_in(ir)
    IF ( ABS(zeta_io(ir))<=1.0_DP ) zeta_io(ir) = SIGN( MIN(ABS(zeta_io(ir)), &
                                    (1.0_DP-rho_threshold)), zeta_io(ir) )
    zeta = zeta_io(ir)
    !
    IF ( ABS(zeta)>1.0_DP .OR. rho<=rho_threshold .OR. SQRT(ABS(grho))<=rho_threshold ) THEN
       sc_out(ir) = 0.0_DP
       v1c_out(ir,:) = 0.0_DP ; v2c_out(ir) = 0.0_DP
       CYCLE
    ENDIF
    !
    SELECT CASE( igcc )
    CASE( 0 )
       !
       sc  = 0.0_DP
       v1c = 0.0_DP
       v2c = 0.0_DP
       !
    CASE( 1 )
       !
       CALL perdew86_spin( rho, zeta, grho, sc, v1c(1), v1c(2), v2c )
       !
    CASE( 2 )
       !
       CALL ggac_spin( rho, zeta, grho, sc, v1c(1), v1c(2), v2c )
       !
    CASE( 4 )
       !
       CALL pbec_spin( rho, zeta, grho, 1, sc, v1c(1), v1c(2), v2c )
       !
    CASE( 8 )
       !
       CALL pbec_spin( rho, zeta, grho, 2, sc, v1c(1), v1c(2), v2c )
       !
#ifdef use_beef
    CASE( 14 )
       !
       call beeflocalcorrspin(rho, zeta, grho, sc, v1c(1), v1c(2), v2c, 0)
#endif
    CASE DEFAULT
       !
       sc = 0.0_DP
       v1c = 0.0_DP
       v2c = 0.0_DP
       !
    END SELECT
    !
    sc_out(ir)  = sc
    v1c_out(ir,:) = v1c(:)
    v2c_out(ir) = v2c
    !
  ENDDO
!$omp end do
!$omp end parallel
  !
  RETURN
  !
END SUBROUTINE gcc_spin
!
!
!---------------------------------------------------------------------------
SUBROUTINE gcc_spin_more( length, rho_in, grho_in, grho_ud_in, &
                                               sc, v1c, v2c, v2c_ud )
  !-------------------------------------------------------------------------
  !! Gradient corrections for exchange and correlation.
  !
  !! * Exchange:
  !!    * Becke88;
  !!    * GGAX.
  !! * Correlation:
  !!    * Perdew86;
  !!    * Lee, Yang & Parr;
  !!    * GGAC.
  !
  USE corr_gga
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: rho_in
  !! the total charge
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: grho_in
  !! the gradient of the charge squared
  REAL(DP), INTENT(IN), DIMENSION(length) :: grho_ud_in
  !! gradient off-diagonal term up-down
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sc
  !! correlation energies
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1c
  !! correlation potential (density part)
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v2c
  !! correlation potential (gradient part)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c_ud
  !!correlation potential (off-diag. term)
  !
  ! ... local variables
  !
  INTEGER :: ir, igcc
  REAL(DP) :: rho(2), grho(2)
  REAL(DP) :: grho_ud
  LOGICAL  :: exx_started
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
#endif    
  !
  igcc = get_igcc()
  sc  = 0.0_DP
  v1c = 0.0_DP
  v2c = 0.0_DP
  v2c_ud = 0.0_DP
  exx_started = exx_is_active()
  !
#if defined(_OPENMP)
  ntids = omp_get_num_threads()
#endif
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, grho, grho_ud )
  DO ir = 1, length
    !
    rho(:) = rho_in(ir,:)
    grho(:) = grho_in(ir,:)
    grho_ud = grho_ud_in(ir)
    !
    IF ( rho(1)+rho(2) < rho_threshold ) THEN
       sc(ir) = 0.0_DP
       v1c(ir,:) = 0.0_DP
       v2c(ir,:) = 0.0_DP ; v2c_ud(ir) = 0.0_DP
       CYCLE
    ENDIF
    !
    CALL lsd_glyp( rho(1), rho(2), grho(1), grho(2), grho_ud, &
                   sc(ir), v1c(ir,1), v1c(ir,2), v2c(ir,1),   &
                   v2c(ir,2), v2c_ud(ir) )
    !
    SELECT CASE( igcc )
    CASE( 3 )
       !
       ! ... void
       !
    CASE( 7 )
       !
       IF ( exx_started ) THEN
          sc(ir) = 0.81_DP * sc(ir)
          v1c(ir,:) = 0.81_DP * v1c(ir,:)
          v2c(ir,:) = 0.81_DP * v2c(ir,:)
          v2c_ud(ir) = 0.81_DP * v2c_ud(ir)
       ENDIF
       !
    CASE( 13 )
       !
       IF ( exx_started ) THEN
          sc(ir) = 0.871_DP * sc(ir)
          v1c(ir,:) = 0.871_DP * v1c(ir,:)
          v2c(ir,:) = 0.871_DP * v2c(ir,:)
          v2c_ud(ir) = 0.871_DP * v2c_ud(ir)
       ENDIF
       !
    CASE DEFAULT
       !
       CALL errore( " gcc_spin_more ", " gradient correction not implemented ", 1 )
       !
    END SELECT
    !
  ENDDO
!$omp end do
!$omp end parallel
  !
  RETURN
  !
END SUBROUTINE gcc_spin_more
!
!
END MODULE xc_gga
