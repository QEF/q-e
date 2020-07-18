!-----------------------------------------------------------------------
!------- DRIVERS FOR DERIVATIVES OF XC POTENTIAL (GGA CASE) ------------
!-----------------------------------------------------------------------
!
!---------------------------------------------------------------------
SUBROUTINE dgcxc( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
  !---------------------------------------------------------------------
  !! Wrapper routine. Calls dgcx-driver routines from internal libraries
  !! or from the external libxc, depending on the input choice.
  !
  USE constants,        ONLY: e2
  USE kinds,            ONLY: DP
  USE funct,            ONLY: get_igcx, get_igcc, is_libxc
  USE xc_gga,           ONLY: gcxc, gcx_spin
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: sp
  !! number of spin components
  REAL(DP), INTENT(IN) :: r_in(length,sp)
  !! charge density
  REAL(DP), INTENT(IN) :: g_in(length,3,sp)
  !! gradient
  REAL(DP), INTENT(OUT) :: dvxc_rr(length,sp,sp), dvxc_sr(length,sp,sp), &
                           dvxc_ss(length,sp,sp)
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: vrrx(:,:), vsrx(:,:), vssx(:,:)
  REAL(DP), ALLOCATABLE :: vrrc(:,:), vsrc(:,:), vssc(:), vrzc(:,:) 
  !
#if defined(__LIBXC)
  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info1, xc_info2
  INTEGER :: fkind
  REAL(DP), ALLOCATABLE :: rho_lbxc(:)
  REAL(DP), ALLOCATABLE :: v2rho2_x(:), v2rhosigma_x(:), v2sigma2_x(:) 
  REAL(DP), ALLOCATABLE :: v2rho2_c(:), v2rhosigma_c(:), v2sigma2_c(:) 
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  INTEGER :: igcx, igcc
  INTEGER :: k, ir, length_lxc, length_dlxc
  REAL(DP) :: rht, zeta
  REAL(DP), ALLOCATABLE :: sigma(:)
  REAL(DP), PARAMETER :: small = 1.E-10_DP, rho_trash = 0.5_DP
  REAL(DP), PARAMETER :: epsr=1.0d-6, epsg=1.0d-10
  !
  igcx = get_igcx()
  igcc = get_igcc()
  !
#if defined(__LIBXC)
  !
  fkind = -1
  lengthxc = length
  !
  IF ( (is_libxc(3) .OR. igcx==0) .AND. (is_libxc(4) .OR. igcc==0)) THEN
    !
    length_lxc = length*sp
    length_dlxc = length
    IF (sp == 2) length_dlxc = length*3
    !
    ALLOCATE( rho_lbxc(length_lxc), sigma(length_dlxc) )
    !
    ! ... set libxc input
    SELECT CASE( sp )
    CASE( 1 )
      !
      DO k = 1, length
        rho_lbxc(k) = r_in(k,1) 
        sigma(k) = g_in(k,1,1)**2 + g_in(k,2,1)**2 + g_in(k,3,1)**2
      ENDDO
      !
    CASE( 2 )
      !
      DO k = 1, length
        rho_lbxc(2*k-1) = r_in(k,1)
        rho_lbxc(2*k)   = r_in(k,2)
        !
        sigma(3*k-2) = g_in(k,1,1)**2 + g_in(k,2,1)**2 + g_in(k,3,1)**2
        sigma(3*k-1) = g_in(k,1,1) * g_in(k,1,2) + g_in(k,2,1) * g_in(k,2,2) + &
                       g_in(k,3,1) * g_in(k,3,2)
        sigma(3*k)   = g_in(k,1,2)**2 + g_in(k,2,2)**2 + g_in(k,3,2)**2
      ENDDO
      !
    CASE( 4 )
      !
      CALL errore( 'dgcxc', 'The derivative of the xc potential with libxc &
                            &is not available for noncollinear case', 1 )
      !
    CASE DEFAULT
      !
      CALL errore( 'dgcxc', 'Wrong number of spin dimensions', 2 )
      !
    END SELECT
    !
    ALLOCATE( v2rho2_x(length_dlxc), v2rhosigma_x(length_dlxc*sp), v2sigma2_x(length_dlxc*sp) )
    ALLOCATE( v2rho2_c(length_dlxc), v2rhosigma_c(length_dlxc*sp), v2sigma2_c(length_dlxc*sp) )
    !
    ! ... DERIVATIVE FOR EXCHANGE
    v2rho2_x = 0.d0 ;  v2rhosigma_x = 0.d0 ;  v2sigma2_x = 0.d0
    IF (igcx /= 0) THEN 
      CALL xc_f03_func_init( xc_func, igcx, sp )
       xc_info1 = xc_f03_func_get_info( xc_func )
       fkind  = xc_f03_func_info_get_kind( xc_info1 )
       CALL xc_f03_gga_fxc( xc_func, lengthxc, rho_lbxc(1), sigma(1), v2rho2_x(1), v2rhosigma_x(1), v2sigma2_x(1) )
      CALL xc_f03_func_end( xc_func )
    ENDIF
    !
    ! ... DERIVATIVE FOR CORRELATION
    v2rho2_c = 0.d0 ;  v2rhosigma_c = 0.d0 ;  v2sigma2_c = 0.d0
    IF (igcc /= 0) THEN 
      CALL xc_f03_func_init( xc_func, igcc, sp )
       xc_info2 = xc_f03_func_get_info( xc_func )
       CALL xc_f03_gga_fxc( xc_func, lengthxc, rho_lbxc(1), sigma(1), v2rho2_c(1), v2rhosigma_c(1), v2sigma2_c(1) )
      CALL xc_f03_func_end( xc_func )
    ENDIF
    !
    dvxc_rr = 0.d0
    dvxc_sr = 0.d0
    dvxc_ss = 0.d0
    !
    IF ( sp == 1 ) THEN
       !
       dvxc_rr(:,1,1) = e2 * (v2rho2_x(:) + v2rho2_c(:))
       dvxc_sr(:,1,1) = e2 * (v2rhosigma_x(:) + v2rhosigma_c(:))*2.d0
       dvxc_ss(:,1,1) = e2 * (v2sigma2_x(:) + v2sigma2_c(:))*4.d0
       !
    ELSEIF ( sp == 2 ) THEN
       !
       DO k = 1, length
         rht = r_in(k,1) + r_in(k,2)
         IF (rht > epsr) THEN
           dvxc_rr(k,1,1) = e2 * (v2rho2_x(3*k-2) + v2rho2_c(3*k-2))
           dvxc_rr(k,1,2) = e2 * (v2rho2_x(3*k-1) + v2rho2_c(3*k-1))
           dvxc_rr(k,2,1) = e2 * (v2rho2_x(3*k-1) + v2rho2_c(3*k-1))
           dvxc_rr(k,2,2) = e2 * (v2rho2_x(3*k)   + v2rho2_c(3*k))
         ENDIF
         !
         dvxc_sr(k,1,1) = e2 * (v2rhosigma_x(6*k-5) + v2rhosigma_c(6*k-5))*2.d0
         dvxc_ss(k,1,1) = e2 * (v2sigma2_x(6*k-5) + v2sigma2_c(6*k))*4.d0
         IF ( fkind==XC_EXCHANGE_CORRELATION ) THEN
            dvxc_sr(k,1,2) = e2 * v2rhosigma_x(6*k-4)
            dvxc_sr(k,2,1) = e2 * v2rhosigma_x(6*k-1)
            dvxc_ss(k,1,2) = e2 * v2sigma2_x(6*k-2)
            dvxc_ss(k,2,1) = e2 * v2sigma2_x(6*k-2)
         ELSE
            dvxc_sr(k,1,2) = e2 * v2rhosigma_c(6*k-4)
            dvxc_sr(k,2,1) = e2 * v2rhosigma_c(6*k-1)
            dvxc_ss(k,1,2) = e2 * v2sigma2_c(6*k-2)
            dvxc_ss(k,2,1) = e2 * v2sigma2_c(6*k-2)
         ENDIF
         dvxc_sr(k,2,2) = e2 * (v2rhosigma_x(6*k) + v2rhosigma_c(6*k))*2.d0
         dvxc_ss(k,2,2) = e2 * (v2sigma2_x(6*k) + v2sigma2_c(6*k))*4.d0
         !
       ENDDO
       !
    ENDIF
    !
    DEALLOCATE( rho_lbxc, sigma )
    DEALLOCATE( v2rho2_x, v2rhosigma_x, v2sigma2_x )
    DEALLOCATE( v2rho2_c, v2rhosigma_c, v2sigma2_c )
    !
  ELSEIF ((.NOT.is_libxc(3)) .AND. (.NOT.is_libxc(4))) THEN
    !
    ALLOCATE( vrrx(length,sp), vsrx(length,sp), vssx(length,sp) )
    ALLOCATE( vrrc(length,sp), vsrc(length,sp), vssc(length) )
    !
    IF ( sp == 1 ) THEN
       !
       ALLOCATE( sigma(length) )
       sigma(:) = g_in(:,1,1)**2 + g_in(:,2,1)**2 + g_in(:,3,1)**2
       CALL dgcxc_unpol( length, r_in, sigma, vrrx, vsrx, vssx, vrrc, vsrc, vssc )
       DEALLOCATE( sigma )
       !
       dvxc_rr(:,1,1) = e2 * (vrrx(:,1) + vrrc(:,1))
       dvxc_sr(:,1,1) = e2 * (vsrx(:,1) + vsrc(:,1))
       dvxc_ss(:,1,1) = e2 * (vssx(:,1) + vssc(:)  )
       !
    ELSEIF ( sp == 2 ) THEN
       !
       ALLOCATE( vrzc(length,sp) )
       !
       CALL dgcxc_spin( length, r_in, g_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc, vrzc )
       !
       DO k = 1, length
         !
         rht = r_in(k,1) + r_in(k,2)
         IF (rht > epsr) THEN
           zeta = (r_in(k,1) - r_in(k,2)) / rht
           !
           dvxc_rr(k,1,1) = e2 * (vrrx(k,1) + vrrc(k,1) + vrzc(k,1) * (1.d0 - zeta) / rht)
           dvxc_rr(k,1,2) = e2 * (vrrc(k,1) - vrzc(k,1) * (1.d0 + zeta) / rht)
           dvxc_rr(k,2,1) = e2 * (vrrc(k,2) + vrzc(k,2) * (1.d0 - zeta) / rht)
           dvxc_rr(k,2,2) = e2 * (vrrx(k,2) + vrrc(k,2) - vrzc(k,2) * (1.d0 + zeta) / rht)
         ENDIF
         !
         dvxc_sr(k,1,1) = e2 * (vsrx(k,1) + vsrc(k,1))   
         dvxc_sr(k,1,2) = e2 * vsrc(k,1)   
         dvxc_sr(k,2,1) = e2 * vsrc(k,2)   
         dvxc_sr(k,2,2) = e2 * (vsrx(k,2) + vsrc(k,2))   
         !
         dvxc_ss(k,1,1) = e2 * (vssx(k,1) + vssc(k))   
         dvxc_ss(k,1,2) = e2 * vssc(k)   
         dvxc_ss(k,2,1) = e2 * vssc(k)   
         dvxc_ss(k,2,2) = e2 * (vssx(k,2) + vssc(k))
         !
       ENDDO
       !
       DEALLOCATE( vrzc )
       !
    ENDIF
    !
    DEALLOCATE( vrrx, vsrx, vssx )
    DEALLOCATE( vrrc, vsrc, vssc )
    !
  ELSE
    !
    CALL errore( 'dgcxc', 'Derivatives of exchange and correlation terms, &
                         & at present, must be both qe or both libxc.', 3 )
    !
  ENDIF
  !
#else
  !
  ALLOCATE( vrrx(length,sp), vsrx(length,sp), vssx(length,sp) )
  ALLOCATE( vrrc(length,sp), vsrc(length,sp), vssc(length) )
  !
  SELECT CASE( sp )
  CASE( 1 )
     !
     ALLOCATE( sigma(length) )
     sigma(:) = g_in(:,1,1)**2 + g_in(:,2,1)**2 + g_in(:,3,1)**2
     CALL dgcxc_unpol( length, r_in, sigma, vrrx, vsrx, vssx, vrrc, vsrc, vssc )
     DEALLOCATE( sigma )
     !
     dvxc_rr(:,1,1) = e2 * (vrrx(:,1) + vrrc(:,1))
     dvxc_sr(:,1,1) = e2 * (vsrx(:,1) + vsrc(:,1))
     dvxc_ss(:,1,1) = e2 * (vssx(:,1) + vssc(:)  )     
     !
  CASE( 2 )
     !
     ALLOCATE( vrzc(length,sp) )
     !
     CALL dgcxc_spin( length, r_in, g_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc, vrzc )
     !
     DO k = 1, length
        !
        rht = r_in(k,1) + r_in(k,2)
        IF (rht > epsr) THEN
           zeta = (r_in(k,1) - r_in(k,2)) / rht
           !
           dvxc_rr(k,1,1) = e2 * (vrrx(k,1) + vrrc(k,1) + vrzc(k,1) * (1.d0 - zeta) / rht)
           dvxc_rr(k,1,2) = e2 * (vrrc(k,1) - vrzc(k,1) * (1.d0 + zeta) / rht)
           dvxc_rr(k,2,1) = e2 * (vrrc(k,2) + vrzc(k,2) * (1.d0 - zeta) / rht)
           dvxc_rr(k,2,2) = e2 * (vrrx(k,2) + vrrc(k,2) - vrzc(k,2) * (1.d0 + zeta) / rht)
        ENDIF
        !
        dvxc_sr(k,1,1) = e2 * (vsrx(k,1) + vsrc(k,1))   
        dvxc_sr(k,1,2) = e2 * vsrc(k,1)   
        dvxc_sr(k,2,1) = e2 * vsrc(k,2)   
        dvxc_sr(k,2,2) = e2 * (vsrx(k,2) + vsrc(k,2))   
        !
        dvxc_ss(k,1,1) = e2 * (vssx(k,1) + vssc(k))   
        dvxc_ss(k,1,2) = e2 * vssc(k)   
        dvxc_ss(k,2,1) = e2 * vssc(k)   
        dvxc_ss(k,2,2) = e2 * (vssx(k,2) + vssc(k))
     ENDDO
     !
     DEALLOCATE( vrzc )
     !
  CASE DEFAULT
     !
     CALL errore( 'dgcxc', 'Wrong ns input', 4 )
     !
  END SELECT
  !
  DEALLOCATE( vrrx, vsrx, vssx )
  DEALLOCATE( vrrc, vsrc, vssc )
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE
!
!
!---------------------------------------------------------------------------
SUBROUTINE dgcxc_unpol( length, r_in, s2_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc )
  !-------------------------------------------------------------------------
  !! This routine computes the derivative of the exchange and correlation
  !! potentials.
  !
  USE kinds,        ONLY: DP
  USE xc_gga,       ONLY: gcxc
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN), DIMENSION(length) :: r_in, s2_in
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrx, vsrx, vssx
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrc, vsrc, vssc
  !
  ! ... local variables
  !
  INTEGER :: i1, i2, i3, i4, f1, f2, f3, f4
  REAL(DP), DIMENSION(length) :: dr, s, ds
  REAL(DP), DIMENSION(4*length) :: raux, s2aux
  REAL(DP), ALLOCATABLE :: v1x(:), v2x(:), v1c(:), v2c(:)
  REAL(DP), ALLOCATABLE :: sx(:), sc(:)
  REAL(DP), PARAMETER :: small = 1.E-30_DP
  !
  ALLOCATE( v1x(4*length), v2x(4*length), sx(4*length) )
  ALLOCATE( v1c(4*length), v2c(4*length), sc(4*length) )
  !
  i1 = 1     ;   f1 = length     !4 blocks:  [ rho+dr ,    grho2    ]
  i2 = f1+1  ;   f2 = 2*length   !           [ rho-dr ,    grho2    ]
  i3 = f2+1  ;   f3 = 3*length   !           [ rho    , (grho+ds)^2 ]
  i4 = f3+1  ;   f4 = 4*length   !           [ rho    , (grho-ds)^2 ]
  !
  s  = SQRT(s2_in)
  dr = MIN(1.d-4, 1.d-2*r_in)
  ds = MIN(1.d-4, 1.d-2*s)
  !
  raux(i1:f1) = r_in+dr  ;   s2aux(i1:f1) = s2_in
  raux(i2:f2) = r_in-dr  ;   s2aux(i2:f2) = s2_in
  raux(i3:f3) = r_in     ;   s2aux(i3:f3) = (s+ds)**2
  raux(i4:f4) = r_in     ;   s2aux(i4:f4) = (s-ds)**2
  !
  CALL gcxc( length*4, raux, s2aux, sx, sc, v1x, v2x, v1c, v2c )
  !
  ! ... to avoid NaN in the next operations
  WHERE( r_in<=small .OR. s2_in<=small )
    dr = 1._DP ; ds = 1._DP ; s = 1._DP
  END WHERE
  !
  vrrx = 0.5_DP * (v1x(i1:f1) - v1x(i2:f2)) / dr
  vrrc = 0.5_DP * (v1c(i1:f1) - v1c(i2:f2)) / dr
  !
  vsrx = 0.25_DP * ((v2x(i1:f1) - v2x(i2:f2)) / dr + &
                    (v1x(i3:f3) - v1x(i4:f4)) / ds / s)
  vsrc = 0.25_DP * ((v2c(i1:f1) - v2c(i2:f2)) / dr + &
                    (v1c(i3:f3) - v1c(i4:f4)) / ds / s)
  !
  vssx = 0.5_DP * (v2x(i3:f3) - v2x(i4:f4)) / ds / s
  vssc = 0.5_DP * (v2c(i3:f3) - v2c(i4:f4)) / ds / s
  !
  DEALLOCATE( v1x, v2x, sx )
  DEALLOCATE( v1c, v2c, sc )
  !
  RETURN
  !
END SUBROUTINE dgcxc_unpol
!
!
!--------------------------------------------------------------------------
SUBROUTINE dgcxc_spin( length, r_in, g_in, vrrx, vrsx, vssx, vrrc, vrsc, &
                       vssc, vrzc )
  !------------------------------------------------------------------------
  !! This routine computes the derivative of the exchange and correlation
  !! potentials in the spin-polarized case.
  !
  USE xc_gga,       ONLY: gcx_spin, gcc_spin
  USE kinds,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: r_in
  REAL(DP), INTENT(IN), DIMENSION(length,3,2) :: g_in
  ! input: the charges and the gradient
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrrx, vrsx, vssx
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrrc, vrsc, vrzc
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vssc
  ! output: derivatives of the exchange and of the correlation
  !
  ! ... local variables
  !
  INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8
  INTEGER :: f1, f2, f3, f4, f5, f6, f7, f8
  ! block delimiters
  REAL(DP), DIMENSION(length,2) :: r, s, s2
  REAL(DP), DIMENSION(length,2) :: drup, drdw, dsup, dsdw
  ! deltas for rho and gradient
  REAL(DP), ALLOCATABLE :: sx(:), v1x(:,:), v2x(:,:)
  ! exchange energy and potentials for each block
  REAL(DP), ALLOCATABLE :: sc(:), v1c(:,:), v2c(:)
  ! correlation energy and potentials for each block
  REAL(DP), DIMENSION(length) :: rt, zeta, st, s2t
  ! rho tot, zeta, gradient, square tot gradient
  REAL(DP), DIMENSION(length) :: dr, ds, dz
  ! deltas for rho tot, gradient and zeta
  REAL(DP), DIMENSION(length,2) :: null_v
  ! used to set output values to zero when input values 
  ! are too small (e.g. rho<eps)
  ! 
  REAL(DP), ALLOCATABLE :: raux(:,:), s2aux(:,:)
  REAL(DP), ALLOCATABLE :: rtaux(:), s2taux(:), zetaux(:)
  ! auxiliary arrays for gcx- and gcc- routines input
  !
  REAL(DP), PARAMETER :: eps = 1.D-6
  REAL(DP), PARAMETER :: rho_trash = 0.4_DP, zeta_trash = 0.2_DP, &
                         s2_trash = 0.1_DP
  !
  vrrx = 0.0_DP ; vrsx = 0.0_DP ; vssx = 0.0_DP
  vrrc = 0.0_DP ; vrsc = 0.0_DP ; vrzc = 0.0_DP
  vssc = 0.0_DP
  !
  ! ... EXCHANGE
  !
  i1 = 1     ;   f1 = length     !8 blocks(x2): [ rho+drup , grho2         ]
  i2 = f1+1  ;   f2 = 2*length   !              [ rho-drup , grho2         ]
  i3 = f2+1  ;   f3 = 3*length   !              [ rho      , (grho+dsup)^2 ]
  i4 = f3+1  ;   f4 = 4*length   !              [ rho      , (grho-dsup)^2 ]
  i5 = f4+1  ;   f5 = 5*length   !              [ rho+drdw , grho2         ]
  i6 = f5+1  ;   f6 = 6*length   !              [ rho-drdw , grho2         ]  
  i7 = f6+1  ;   f7 = 7*length   !              [ rho      , (grho+dsdw)^2 ]
  i8 = f7+1  ;   f8 = 8*length   !              [ rho      , (grho-dsdw)^2 ]    
  !
  ALLOCATE( raux(length*8,2), s2aux(length*8,2) )
  ALLOCATE( sx(length*8), v1x(length*8,2), v2x(length*8,2) )
  !
  s2(:,1) = g_in(:,1,1)**2 + g_in(:,2,1)**2 + g_in(:,3,1)**2 !up
  s2(:,2) = g_in(:,1,2)**2 + g_in(:,2,2)**2 + g_in(:,3,2)**2 !down
  s = SQRT(s2)
  r = r_in
  !
  null_v = 1.0_DP
  !
  ! ... thresholds
  WHERE ( r_in(:,1)<=eps .OR. s(:,1)<=eps )
     r(:,1) = rho_trash
     s2(:,1) = s2_trash ; s(:,1) = SQRT(s2_trash)
     null_v(:,1) = 0.0_DP
  END WHERE
  !
  WHERE ( r_in(:,2)<=eps .OR. s(:,2)<=eps )
     r(:,2) = rho_trash
     s2(:,2) = s2_trash ; s(:,2) = SQRT(s2_trash)
     null_v(:,2) = 0.0_DP
  END WHERE
  !
  drup = 0.0_DP ;  drdw = 0.0_DP
  dsup = 0.0_DP ;  dsdw = 0.0_DP
  !
  drup(:,1) = MIN(1.D-4, 1.D-2*r(:,1)) ; dsup(:,1) = MIN(1.D-4, 1.D-2*s(:,1))
  drdw(:,2) = MIN(1.D-4, 1.D-2*r(:,2)) ; dsdw(:,2) = MIN(1.D-4, 1.D-2*s(:,2))
  !
  ! ... up
  raux(i1:f1,:) = r+drup ;  s2aux(i1:f1,:) = s2
  raux(i2:f2,:) = r-drup ;  s2aux(i2:f2,:) = s2
  raux(i3:f3,:) = r      ;  s2aux(i3:f3,:) = (s+dsup)**2
  raux(i4:f4,:) = r      ;  s2aux(i4:f4,:) = (s-dsup)**2
  ! ... down
  raux(i5:f5,:) = r+drdw ;  s2aux(i5:f5,:) = s2
  raux(i6:f6,:) = r-drdw ;  s2aux(i6:f6,:) = s2
  raux(i7:f7,:) = r      ;  s2aux(i7:f7,:) = (s+dsdw)**2
  raux(i8:f8,:) = r      ;  s2aux(i8:f8,:) = (s-dsdw)**2
  !
  !
  CALL gcx_spin( length*8, raux, s2aux, sx, v1x, v2x )
  !
  ! ... up
  vrrx(:,1) = 0.5_DP  *  (v1x(i1:f1,1) - v1x(i2:f2,1)) / drup(:,1)
  vrsx(:,1) = 0.25_DP * ((v2x(i1:f1,1) - v2x(i2:f2,1)) / drup(:,1) + &
                         (v1x(i3:f3,1) - v1x(i4:f4,1)) / dsup(:,1) / s(:,1))
  vssx(:,1) = 0.5_DP  *  (v2x(i3:f3,1) - v2x(i4:f4,1)) / dsup(:,1) / s(:,1)
  ! ... down
  vrrx(:,2) = 0.5_DP  *  (v1x(i5:f5,2) - v1x(i6:f6,2)) / drdw(:,2)
  vrsx(:,2) = 0.25_DP * ((v2x(i5:f5,2) - v2x(i6:f6,2)) / drdw(:,2) + &
                                     (v1x(i7:f7,2) - v1x(i8:f8,2)) / dsdw(:,2) / s(:,2))
  vssx(:,2) = 0.5_DP  *  (v2x(i7:f7,2) - v2x(i8:f8,2)) / dsdw(:,2) / s(:,2)
  !
  vrrx(:,1) = vrrx(:,1)*null_v(:,1) ;  vrrx(:,2) = vrrx(:,2)*null_v(:,2)
  vrsx(:,1) = vrsx(:,1)*null_v(:,1) ;  vrsx(:,2) = vrsx(:,2)*null_v(:,2)
  vssx(:,1) = vssx(:,1)*null_v(:,1) ;  vssx(:,2) = vssx(:,2)*null_v(:,2)
  !
  DEALLOCATE( raux, s2aux  )
  DEALLOCATE( sx, v1x, v2x )
  !
  ! ... CORRELATION
  !
  i1 = 1     ;   f1 = length     !6 blocks: [ rt+dr , s2t       , zeta    ]
  i2 = f1+1  ;   f2 = 2*length   !          [ rt-dr , s2t       , zeta    ]
  i3 = f2+1  ;   f3 = 3*length   !          [ rt    , (st+ds)^2 , zeta    ]
  i4 = f3+1  ;   f4 = 4*length   !          [ rt    , (st-ds)^2 , zeta    ]
  i5 = f4+1  ;   f5 = 5*length   !          [ rt    , grho2     , zeta+dz ]
  i6 = f5+1  ;   f6 = 6*length   !          [ rt    , grho2     , zeta-dz ]  
  !
  ALLOCATE( rtaux(length*6), s2taux(length*6), zetaux(length*6) )
  ALLOCATE( v1c(length*6,2), v2c(length*6), sc(length*6) )
  !
  rt(:) = r_in(:,1) + r_in(:,2)
  !
  null_v = 1.0_DP
  !
  WHERE (rt > eps)
     zeta = (r_in(:,1) - r_in(:,2)) / rt(:)
  ELSEWHERE
     zeta = zeta_trash
     null_v(:,1) = 0.0_DP
  END WHERE
  !
  s2t = (g_in(:,1,1) + g_in(:,1,2))**2 + &
        (g_in(:,2,1) + g_in(:,2,2))**2 + &
        (g_in(:,3,1) + g_in(:,3,2))**2
  st = SQRT(s2t)
  !
  WHERE (rt<eps .OR. ABS(zeta)>1._DP .OR. st<eps)
     rt(:)  = rho_trash
     s2t(:) = s2_trash ; st = SQRT(s2_trash)
     zeta(:) = zeta_trash
     null_v(:,1) = 0.0_DP
  END WHERE
  !
  dr = MIN(1.D-4, 1.D-2 * rt)
  ds = MIN(1.D-4, 1.D-2 * st)
  !dz = MIN(1.D-4, 1.D-2 * ABS(zeta) )
  dz = 1.D-6
  !
  ! ... If zeta is too close to +-1 the derivative is evaluated at a
  ! slightly smaller value.
  zeta = SIGN( MIN(ABS(zeta), (1.0_DP - 2.0_DP*dz)), zeta )
  !
  rtaux(i1:f1) = rt+dr ;  s2taux(i1:f1) = s2t        ;  zetaux(i1:f1) = zeta
  rtaux(i2:f2) = rt-dr ;  s2taux(i2:f2) = s2t        ;  zetaux(i2:f2) = zeta
  rtaux(i3:f3) = rt    ;  s2taux(i3:f3) = (st+ds)**2 ;  zetaux(i3:f3) = zeta
  rtaux(i4:f4) = rt    ;  s2taux(i4:f4) = (st-ds)**2 ;  zetaux(i4:f4) = zeta
  rtaux(i5:f5) = rt    ;  s2taux(i5:f5) = s2t        ;  zetaux(i5:f5) = zeta+dz
  rtaux(i6:f6) = rt    ;  s2taux(i6:f6) = s2t        ;  zetaux(i6:f6) = zeta-dz
  !
  CALL gcc_spin( length*6, rtaux, zetaux, s2taux, sc, v1c, v2c )
  !
  vrrc(:,1) = 0.5_DP * (v1c(i1:f1,1) - v1c(i2:f2,1)) / dr    * null_v(:,1)
  vrrc(:,2) = 0.5_DP * (v1c(i1:f1,2) - v1c(i2:f2,2)) / dr    * null_v(:,1)
  vrsc(:,1) = 0.5_DP * (v1c(i3:f3,1) - v1c(i4:f4,1)) / ds/st * null_v(:,1)
  vrsc(:,2) = 0.5_DP * (v1c(i3:f3,2) - v1c(i4:f4,2)) / ds/st * null_v(:,1)
  vssc(:)   = 0.5_DP * (v2c(i3:f3)   - v2c(i4:f4)  ) / ds/st * null_v(:,1)
  vrzc(:,1) = 0.5_DP * (v1c(i5:f5,1) - v1c(i6:f6,1)) / dz    * null_v(:,1)
  vrzc(:,2) = 0.5_DP * (v1c(i5:f5,2) - v1c(i6:f6,2)) / dz    * null_v(:,1)
  !
  RETURN
  !
END SUBROUTINE dgcxc_spin
!
!
!-----------------------------------------------------------------------
SUBROUTINE d3gcxc( r, s2, vrrrx, vsrrx, vssrx, vsssx, &
                   vrrrc, vsrrc, vssrc, vsssc )
  !-----------------------------------------------------------------------
  !    wat20101006: Calculates all derivatives of the exchange (x) and
  !                 correlation (c) potential in third order.
  !                 of the Exc.
  !
  !    input:       r = rho, s2=|\nabla rho|^2
  !    definition:  E_xc = \int ( f_x(r,s2) + f_c(r,s2) ) dr
  !    output:      vrrrx = d^3(f_x)/d(r)^3
  !                 vsrrx = d^3(f_x)/d(|\nabla r|)d(r)^2 / |\nabla r|
  !                 vssrx = d/d(|\nabla r|) [ &
  !                           d^2(f_x)/d(|\nabla r|)d(r) / |\nabla r| ] &
  !                                                           / |\nabla r|
  !                 vsssx = d/d(|\nabla r|) [ &
  !                           d/d(|\nabla r|) [ &
  !                           d(f_x)/d(|\nabla r|) / |\nabla r| ] &
  !                                                    / |\nabla r| ] &
  !                                                        / |\nabla r|
  !                 same for (c)
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: r, s2
  REAL(DP) :: vrrrx, vsrrx, vssrx, vsssx, vrrrc, vsrrc, vssrc, vsssc
  !
  ! ... local variables
  !
  REAL(DP) :: dr, s, ds
  REAL(DP), DIMENSION(4) :: raux, s2aux
  REAL(DP), DIMENSION(4) :: vrrx_rs, vsrx_rs, vssx_rs, vrrc_rs, &
                            vsrc_rs, vssc_rs
  !
  s = SQRT(s2)
  dr = MIN(1.d-4, 1.d-2 * r)
  ds = MIN(1.d-4, 1.d-2 * s)
  !
  raux(1) = r+dr  ; s2aux(1) = s2           !4 blocks:  [ rho+dr ,    grho2    ]
  raux(2) = r-dr  ; s2aux(2) = s2           !           [ rho-dr ,    grho2    ]
  raux(3) = r     ; s2aux(3) = (s+ds)**2    !           [ rho    , (grho+ds)^2 ]
  raux(4) = r     ; s2aux(4) = (s-ds)**2    !           [ rho    , (grho-ds)^2 ]
  !
  CALL dgcxc_unpol( 4, raux, s2aux, vrrx_rs, vsrx_rs, vssx_rs, vrrc_rs, vsrc_rs, vssc_rs )
  !
  vrrrx = 0.5d0  *  (vrrx_rs(1) - vrrx_rs(2)) / dr
  vsrrx = 0.25d0 * ((vsrx_rs(1) - vsrx_rs(2)) / dr &
                   +(vrrx_rs(3) - vrrx_rs(4)) / ds / s)
  vssrx = 0.25d0 * ((vssx_rs(1) - vssx_rs(2)) / dr &
                   +(vsrx_rs(3) - vsrx_rs(4)) / ds / s)
  vsssx = 0.5d0  *  (vssx_rs(3) - vssx_rs(4)) / ds / s
  !
  vrrrc = 0.5d0  *  (vrrc_rs(1) - vrrc_rs(2)) / dr
  vsrrc = 0.25d0 * ((vsrc_rs(1) - vsrc_rs(2)) / dr &
                   +(vrrc_rs(3) - vrrc_rs(4)) / ds / s)
  vssrc = 0.25d0 * ((vssc_rs(1) - vssc_rs(2)) / dr &
                   +(vsrc_rs(3) - vsrc_rs(4)) / ds / s)
  vsssc = 0.5d0  *  (vssc_rs(3) - vssc_rs(4)) / ds / s
  !
  RETURN
  !
END SUBROUTINE d3gcxc
