MODULE xc_gga
!
USE kinds,     ONLY: DP 
!
IMPLICIT NONE
!
PRIVATE
SAVE
!
!  GGA exchange-correlation drivers
PUBLIC :: xc_gcx, gcxc, gcx_spin, gcc_spin, gcc_spin_more, &
          select_gga_functionals, change_threshold_gga
!
PUBLIC :: libxc_switches_gga
PUBLIC :: igcx_l, igcc_l
PUBLIC :: exx_started_g, exx_fraction_g
PUBLIC :: screening_parameter_l, gau_parameter_l
!
!  libxc on/off
INTEGER  :: libxc_switches_gga(2)
!
!  indexes defining xc functionals
INTEGER  :: igcx_l, igcc_l
!
!  input thresholds (default values)
REAL(DP) :: rho_threshold = 1.D-6
REAL(DP) :: grho_threshold = 1.D-10
!
!  variables for hybrid exchange
LOGICAL  :: exx_started_g
REAL(DP) :: exx_fraction_g
!
!  screening_ and gau_parameters
REAL(DP) :: screening_parameter_l, gau_parameter_l
!
!
 CONTAINS
!
!
!----------------------------------------------------------------------------
!----- Select functionals by the corresponding indexes ----------------------
!----------------------------------------------------------------------------
SUBROUTINE select_gga_functionals( igcx, igcc, exx_fraction, screening_parameter, &
                                   gau_parameter )
   !-----------------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   INTEGER,  INTENT(IN) :: igcx, igcc
   REAL(DP), INTENT(IN), OPTIONAL :: exx_fraction
   REAL(DP), INTENT(IN), OPTIONAL :: screening_parameter
   REAL(DP), INTENT(IN), OPTIONAL :: gau_parameter
   !
   ! exchange-correlation indexes
   igcx_l = igcx
   igcc_l = igcc
   !
   ! hybrid exchange vars
   exx_started_g  = .FALSE.
   exx_fraction_g = 0._DP
   IF ( PRESENT(exx_fraction) ) THEN
      exx_started_g  = .TRUE.
      exx_fraction_g = exx_fraction
   ENDIF
   !
   ! screening_ and gau_parameter
   screening_parameter_l = 0.0_DP
   gau_parameter_l = 0.0_DP
   !
   IF ( PRESENT(screening_parameter) ) THEN
      screening_parameter_l = screening_parameter
   ENDIF
   !
   IF ( PRESENT(gau_parameter) ) THEN
      gau_parameter_l = gau_parameter
   ENDIF
   !
   RETURN
   !
END SUBROUTINE select_gga_functionals
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
  !! NOTE: look at 'PP/src/benchmark_libxc.f90' to test and see the differences
  !!       between q-e and libxc libraries.
  !
#if defined(__LIBXC)
  USE xc_f90_types_m
  USE xc_f90_lib_m
#endif
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
  !! exchange potential
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  !! exchange
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  !! correlation potential
  REAL(DP), INTENT(OUT) :: v2c(length,ns)
  !! correlation
  REAL(DP), INTENT(OUT), OPTIONAL :: v2c_ud(length)
  !! correlation
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info1, xc_info2
  REAL(DP), ALLOCATABLE :: rho_lxc(:), sigma(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vx_sigma(:)
  REAL(DP), ALLOCATABLE :: vc_rho(:), vc_sigma(:)
  !
  INTEGER :: np
  REAL(DP) :: rs, rtot, zet, sgn(2), vc_2(2)
  REAL(DP), PARAMETER :: pi34 = 0.6203504908994_DP
  !
  LOGICAL :: POLARIZED
  INTEGER :: ildax, ildac, pol_unpol
#endif
  REAL(DP), ALLOCATABLE :: arho(:,:), sign_v(:)
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grho2(:,:), grho_ud(:)
  !
  INTEGER :: k, is
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  !
  IF (ns==2 .AND. .NOT. PRESENT(v2c_ud)) CALL errore( 'xc_gga', 'cross &
                                             &term v2c_ud not found', 1 )
  !
#if defined(__LIBXC)
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
  IF ( ns == 1 ) ALLOCATE( sign_v(length) )
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
       IF ( rho_lxc(k) > rho_threshold ) THEN
          sigma(k) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
          IF ( sigma(k) > grho_threshold ) THEN
             sign_v(k) = SIGN( 1._DP, rho(k,1) )
          ELSE
             rho_lxc(k) = 0.5_DP
             sigma(k) = 0.1_DP
             sign_v(k) = 0.0_DP
          ENDIF
       ELSE
          rho_lxc(k) = 0.5_DP
          sigma(k) = 0.1_DP
          sign_v(k) = 0.0_DP
       ENDIF
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
  IF ( ns==1 .AND. SUM(libxc_switches_gga(:))/=2) &
              CALL gcxc( length, rho_lxc, sigma, ex, ec, v1x(:,1), v2x(:,1), v1c(:,1), v2c(:,1) )  
  !
  ! --- GGA EXCHANGE
  !
  IF ( libxc_switches_gga(1) == 1 ) THEN
    !
    CALL xc_f90_func_init( xc_func, xc_info1, igcx_l, pol_unpol )
    CALL xc_f90_func_set_dens_threshold( xc_func, rho_threshold )
     CALL xc_f90_gga_exc_vxc( xc_func, length, rho_lxc(1), sigma(1), ex_lxc(1), vx_rho(1), vx_sigma(1) )
    CALL xc_f90_func_end( xc_func )
    !
    IF (.NOT. POLARIZED) THEN
      DO k = 1, length
        ex(k) = ex_lxc(k) * rho_lxc(k) * sign_v(k)
        sgn(1) = ABS(sign_v(k))
        v1x(k,1) = vx_rho(k)*sgn(1)
        v2x(k,1) = vx_sigma(k)*sgn(1)
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
    IF ( ns == 1 ) THEN
       !
       ! ... This is the spin-unpolarised case
       !
       ex = ex*sign_v
       DO k = 1, length
         sgn(1) = ABS(sign_v(k))
         v1x(k,1) = v1x(k,1)*sgn(1)  ;  v2x(k,1) = v2x(k,1)*sgn(1)
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
    ENDIF
    !
    DEALLOCATE( grho2 )
    !
  ENDIF
  !
  ! ---- GGA CORRELATION
  !
  IF ( libxc_switches_gga(2) == 1 ) THEN  !lda part of LYP not present in libxc
    !
    CALL xc_f90_func_init( xc_func, xc_info2, igcc_l, pol_unpol )
    CALL xc_f90_func_set_dens_threshold( xc_func, rho_threshold )
     CALL xc_f90_gga_exc_vxc( xc_func, length, rho_lxc(1), sigma(1), ec_lxc(1), vc_rho(1), vc_sigma(1) )
    CALL xc_f90_func_end( xc_func )
    !
    IF (.NOT. POLARIZED) THEN
      DO k = 1, length
        ec(k) = ec_lxc(k) * rho_lxc(k) * sign_v(k)
        sgn(1) = ABS(sign_v(k))
        v1c(k,1) = vc_rho(k) * sgn(1)
        v2c(k,1) = vc_sigma(k) * sgn(1)
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
  ELSE
    !
    ALLOCATE( arho(length,ns), grho2(length,ns) )
    !
    IF ( ns == 1 ) THEN
       !
       ! ... This is the spin-unpolarised case
       !
       ec = ec*sign_v
       sign_v = ABS(sign_v)
       v1c(:,1) = v1c(:,1)*sign_v  ;  v2c(:,1) = v2c(:,1)*sign_v
       !
    ELSE
       !
       DO is = 1, 2
          grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
       ENDDO
       !
       IF (igcc_l==3 .OR. igcc_l==7 .OR. igcc_l==13 ) THEN
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
  IF (ns == 1) DEALLOCATE( sign_v )
  DEALLOCATE( ex_lxc , ec_lxc   )
  DEALLOCATE( vx_rho , vx_sigma )
  DEALLOCATE( vc_rho , vc_sigma )
  !
#else
  !
  ALLOCATE( arho(length,ns), grho2(length,ns) )
  !
  IF ( ns == 1 ) THEN
     !
     ! ... This is the spin-unpolarised case
     ALLOCATE( sign_v(length) )
     !
     DO k = 1, length
        arho(k,1) = ABS( rho(k,1) )
        IF ( arho(k,1) > rho_threshold ) THEN
           grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
           IF ( grho2(k,1) > grho_threshold ) THEN
              sign_v(k) = SIGN( 1._DP, rho(k,1) )
           ELSE
              arho(k,1)  = 0.5_DP
              grho2(k,1) = 0.1_DP
              sign_v(k)  = 0.0_DP
           ENDIF
        ELSE
           arho(k,1)  = 0.5_DP
           grho2(k,1) = 0.1_DP
           sign_v(k)  = 0.0_DP
        ENDIF
     ENDDO
     !
     CALL gcxc( length, arho(:,1), grho2(:,1), ex, ec, v1x(:,1), v2x(:,1), v1c(:,1), v2c(:,1) )
     !
     ex = ex*sign_v              ;  ec = ec * sign_v
     sign_v = ABS(sign_v)
     v1x(:,1) = v1x(:,1)*sign_v  ;  v2x(:,1) = v2x(:,1)*sign_v
     v1c(:,1) = v1c(:,1)*sign_v  ;  v2c(:,1) = v2c(:,1)*sign_v
     !
     DEALLOCATE( sign_v )
     !
  ELSE
     !
     DO is = 1, 2
        grho2(:,is) = grho(1,:,is)**2 + grho(2,:,is)**2 + grho(3,:,is)**2
     ENDDO
     !
     CALL gcx_spin( length, rho, grho2, ex, v1x, v2x )
     !
     IF (igcc_l==3 .OR. igcc_l==7 .OR. igcc_l==13 ) THEN
        !
        ALLOCATE( grho_ud(length) )
        !
        grho_ud = grho(1,:,1) * grho(1,:,2) + grho(2,:,1) * grho(2,:,2) + &
                  grho(3,:,1) * grho(3,:,2)
        !
        arho = rho
        WHERE ( rho(:,1)+rho(:,2) < rho_threshold )
           arho(:,1) = 0.0_DP !trash value
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
#endif
  !
  !
  RETURN
  !
END SUBROUTINE xc_gcx
!
!
!-----------------------------------------------------------------------
!------- GRADIENT CORRECTIONS DRIVERS ----------------------------------
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
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in, grho_in
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sx_out, sc_out, v1x_out, &
                                              v2x_out, v1c_out, v2c_out
  !
  ! ... local variables
  !
  INTEGER :: ir
  REAL(DP) :: rho, grho
  REAL(DP) :: sx, v1x, v2x
  REAL(DP) :: sx_, v1x_, v2x_
  REAL(DP) :: sxsr, v1xsr, v2xsr
  REAL(DP) :: sc, v1c, v2c
  REAL(DP), PARAMETER :: small = 1.E-10_DP
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
#endif
  !
#if defined(_OPENMP)
  ntids = omp_get_num_threads()
#endif
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, grho, sx, sx_, sxsr, v1x, v1x_, v1xsr, &
!$omp             v2x, v2x_, v2xsr, sc, v1c, v2c )
  DO ir = 1, length  
     !
     rho  = rho_in(ir)
     grho = grho_in(ir)
     !
     IF ( rho <= small ) THEN
        sx_out(ir)  = 0.0_DP ;   sc_out(ir)  = 0.0_DP
        v1x_out(ir) = 0.0_DP ;   v1c_out(ir) = 0.0_DP
        v2x_out(ir) = 0.0_DP ;   v2c_out(ir) = 0.0_DP
        CYCLE
     ENDIF
     !
     !
     ! ... EXCHANGE
     !  
     SELECT CASE( igcx_l )
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
        IF (igcc_l == 5) CALL hcth( rho, grho, sx, v1x, v2x )
        !
     CASE( 6 )
        !
        CALL optx( rho, grho, sx, v1x, v2x )
        !
     ! case igcx_l == 7 (meta-GGA) must be treated in a separate call to another
     ! routine: needs kinetic energy density in addition to rho and grad rho
     CASE( 8 ) ! 'PBE0'
        !
        CALL pbex( rho, grho, 1, sx, v1x, v2x )
        IF (exx_started_g) THEN
           sx  = (1.0_DP - exx_fraction_g) * sx
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
        ENDIF
        !
     CASE( 9 ) ! 'B3LYP'
        !
        CALL becke88( rho, grho, sx, v1x, v2x )
        IF (exx_started_g) THEN
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
        IF (exx_started_g) THEN
          CALL pbexsr( rho, grho, sxsr, v1xsr, v2xsr, screening_parameter_l )
          sx  = sx  - exx_fraction_g * sxsr
          v1x = v1x - exx_fraction_g * v1xsr
          v2x = v2x - exx_fraction_g * v2xsr
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
        IF (exx_started_g) THEN
          CALL pbexgau( rho, grho, sxsr, v1xsr, v2xsr, gau_parameter_l )
          sx  = sx  - exx_fraction_g * sxsr
          v1x = v1x - exx_fraction_g * v1xsr
          v2x = v2x - exx_fraction_g * v2xsr
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
        IF (exx_started_g) THEN
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
        IF (exx_started_g) THEN
           sx  = (1.0_DP - exx_fraction_g) * sx
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
        ENDIF
        !
     CASE( 30 ) ! 'r860'
        !
        CALL rPW86( rho, grho, sx, v1x, v2x )
        !
        IF (exx_started_g) then
           sx  = (1.0_DP - exx_fraction_g) * sx
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
        ENDIF
        !
     CASE( 38 ) ! 'BR0'
        !
        CALL b86b( rho, grho, 3, sx, v1x, v2x )
        IF (exx_started_g) THEN
           sx  = (1.0_DP - exx_fraction_g) * sx
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
        ENDIF
        !
     CASE( 40 ) ! 'c090'
        !
        CALL c09x( rho, grho, sx, v1x, v2x )
        IF (exx_started_g) THEN
           sx  = (1.0_DP - exx_fraction_g) * sx
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
        ENDIF
        !
     CASE( 41 ) ! 'B86BPBEX'
        !
        CALL becke86b( rho, grho, sx, v1x, v2x )
        IF (exx_started_g) THEN
           sx  = (1.0_DP - exx_fraction_g) * sx
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
        ENDIF
        !
     CASE( 42 ) ! 'BHANDHLYP'
        !
        CALL becke88( rho, grho, sx, v1x, v2x )
        IF (exx_started_g) THEN
           sx  = (1.0_DP - exx_fraction_g) * sx
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
        ENDIF
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
     SELECT CASE( igcc_l )
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
     ! igcc_l == 5 (HCTH) is calculated together with case igcx_l=5
     ! igcc_l == 6 (meta-GGA) is treated in a different routine
     CASE( 7 ) !'B3LYP'
        !
        CALL glyp( rho, grho, sc, v1c, v2c )
        IF (exx_started_g) THEN
           sc  = 0.81_DP * sc
           v1c = 0.81_DP * v1c
           v2c = 0.81_DP * v2c
        ENDIF
        !
     CASE( 8 ) ! 'PBEsol'
        !
        CALL pbec( rho, grho, 2, sc, v1c, v2c )
        !
     ! igcc_l ==  9 set to 5, back-compatibility
     ! igcc_l == 10 set to 6, back-compatibility
     ! igcc_l == 11 M06L calculated in another routine
     CASE( 12 ) ! 'Q2D'
        !
        CALL pbec( rho, grho, 3, sc, v1c, v2c )
        !
     CASE( 13 ) !'X3LYP'
        !
        CALL glyp( rho, grho, sc, v1c, v2c )
        IF (exx_started_g) THEN
           sc  = 0.871_DP * sc
           v1c = 0.871_DP * v1c
           v2c = 0.871_DP * v2c
        ENDIF
        !
     CASE DEFAULT
        !
        sc = 0.0_DP
        v1c = 0.0_DP
        v2c = 0.0_DP
        !
     END SELECT
     !
     sx_out(ir)  = sx  ;  sc_out(ir)  = sc
     v1x_out(ir) = v1x ;  v1c_out(ir) = v1c
     v2x_out(ir) = v2x ;  v2c_out(ir) = v2c
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
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Length of the input/output arrays
  REAL(DP), INTENT(IN),  DIMENSION(length,2) :: rho_in
  !! Up and down charge density
  REAL(DP), INTENT(IN),  DIMENSION(length,2) :: grho2_in
  !! Up and down gradient of the charge
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sx_tot
  !! Energy
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1x_out
  !! Derivatives of exchange wr. rho
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v2x_out
  !! Derivatives of exchange wr. grho
  !
  ! ... local variables
  !
  INTEGER :: ir, is, iflag
  REAL(DP) :: rho(2), grho2(2)
  REAL(DP) :: v1x(2), v2x(2)
  REAL(DP) :: sx(2), null_v(2)
  REAL(DP) :: sxsr(2)
  REAL(DP) :: v1xsr(2), v2xsr(2)
  !
  REAL(DP), PARAMETER :: small=1.E-10_DP
  REAL(DP), PARAMETER :: rho_trash=0.5_DP, grho2_trash=0.2_DP
  ! temporary values assigned to rho and grho when they
  ! are too small in order to avoid numerical problems.
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
#endif    
  !
  sx_tot = 0.0_DP
  !
#if defined(_OPENMP)
  ntids = omp_get_num_threads()
#endif
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, grho2, null_v, sx, sxsr, v1x, v1xsr, &
!$omp             v2x, v2xsr )
  DO ir = 1, length  
     !
     rho(:) = rho_in(ir,:)
     grho2(:) = grho2_in(ir,:)
     null_v(:) = 1.0_DP
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
             null_v(is) = 0.0_DP
           ENDIF
        ENDDO
     ENDIF
     !
     !
     ! ... exchange
     !
     SELECT CASE( igcx_l )
     CASE( 0 )
        !
        sx_tot(ir) = 0.0_DP
        v1x = 0.0_DP
        v2x = 0.0_DP
        !
     CASE( 1 )
        !
        CALL becke88_spin( rho, grho2, sx, v1x, v2x )
        !
        sx_tot(ir) = sx(1)*null_v(1) + sx(2)*null_v(2)
        !
     CASE( 2 )
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL ggax( rho(1), grho2(1), sx(1), v1x(1), v2x(1) )
        CALL ggax( rho(2), grho2(2), sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 3, 4, 8, 10, 12, 20, 23, 24, 25 )
        ! igcx_l=3: PBE, igcx_l=4: revised PBE, igcx_l=8: PBE0, igcx_l=10: PBEsol
        ! igcx_l=12: HSE,  igcx_l=20: gau-pbe, igcx_l=23: obk8, igcx_l=24: ob86, igcx_l=25: ev93
        !
        iflag = 1
        IF ( igcx_l== 4 ) iflag = 2
        IF ( igcx_l==10 ) iflag = 3
        IF ( igcx_l==23 ) iflag = 5
        IF ( igcx_l==24 ) iflag = 6
        IF ( igcx_l==25 ) iflag = 7
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL pbex( rho(1), grho2(1), iflag, sx(1), v1x(1), v2x(1) )
        CALL pbex( rho(2), grho2(2), iflag, sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( igcx_l == 8 .AND. exx_started_g ) THEN
           !
           sx_tot(ir) = (1.0_DP - exx_fraction_g) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
           !
        ELSEIF ( igcx_l == 12 .AND. exx_started_g ) THEN
           !
           CALL pbexsr( rho(1), grho2(1), sxsr(1), v1xsr(1), &
                                          v2xsr(1), screening_parameter_l )
           CALL pbexsr( rho(2), grho2(2), sxsr(2), v1xsr(2), &
                                          v2xsr(2), screening_parameter_l )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction_g*0.5_DP * ( sxsr(1)*null_v(1) + &
                                                               sxsr(2)*null_v(2) )
           v1x = v1x - exx_fraction_g * v1xsr
           v2x = v2x - exx_fraction_g * v2xsr * 2.0_DP
           !
        ELSEIF ( igcx_l == 20 .AND. exx_started_g ) THEN
           ! gau-pbe
           !CALL pbexgau_lsd( rho, grho2, sxsr, v1xsr, v2xsr, gau_parameter_l )
           CALL pbexgau( rho(1), grho2(1), sxsr(1), v1xsr(1), v2xsr(1), gau_parameter_l )
           CALL pbexgau( rho(2), grho2(2), sxsr(2), v1xsr(2), v2xsr(2), gau_parameter_l )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction_g*0.5_DP * ( sxsr(1)*null_v(1) + &
                                                               sxsr(2)*null_v(2) )
           v1x = v1x - exx_fraction_g * v1xsr
           v2x = v2x - exx_fraction_g * v2xsr * 2.0_DP
           !
        ENDIF
        !
     CASE( 9 )                    ! B3LYP
        !
        CALL becke88_spin( rho, grho2, sx, v1x, v2x )
        !
        sx_tot(ir) = sx(1)*null_v(1) + sx(2)*null_v(2)
        !
        IF ( exx_started_g ) THEN
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
      CASE( 26 )                  ! 'B86R for rev-vdW-DF2'
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL b86b( rho(1), grho2(1), 3, sx(1), v1x(1), v2x(1) )
        CALL b86b( rho(2), grho2(2), 3, sx(2), v1x(2), v2x(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
     CASE( 28 )                   ! X3LYP
        !
        CALL becke88_spin( rho, grho2, sx, v1x, v2x )
        !
        rho = 2.0_DP * rho
        grho2 = 4.0_DP * grho2
        !
        CALL pbex( rho(1), grho2(1), 1, sxsr(1), v1xsr(1), v2xsr(1) )
        CALL pbex( rho(2), grho2(2), 1, sxsr(2), v1xsr(2), v2xsr(2) )
        !
        sx_tot(ir) = 0.5_DP * ( sxsr(1)*null_v(1) + sxsr(2)*null_v(2) ) * 0.235_DP + &
                              (   sx(1)*null_v(1) +   sx(2)*null_v(2) ) * 0.765_DP
        v1x = v1xsr * 0.235_DP + v1x * 0.765_DP
        v2x = v2xsr * 0.235_DP * 2.0_DP + v2x * 0.765_DP
        !
        IF ( exx_started_g ) THEN
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started_g ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction_g) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started_g ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction_g) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started_g ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction_g) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
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
        sx_tot(ir) = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started_g ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction_g) * sx_tot(ir)
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
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
        sx_tot = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started_g ) THEN
           sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
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
        sx_tot = 0.5_DP * ( sx(1)*null_v(1) + sx(2)*null_v(2) )
        v2x = 2.0_DP * v2x
        !
        IF ( exx_started_g ) THEN
           sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
           v1x = (1.0_DP - exx_fraction_g) * v1x
           v2x = (1.0_DP - exx_fraction_g) * v2x
        ENDIF
        !
     ! case igcx_l == 5 (HCTH) and 6 (OPTX) not implemented
     ! case igcx_l == 7 (meta-GGA) must be treated in a separate call to another
     ! routine: needs kinetic energy density in addition to rho and grad rho
     CASE DEFAULT
        !
        CALL errore( 'gcx_spin', 'not implemented', igcx_l )
        !
     END SELECT
     !
     v1x_out(ir,:) = v1x(:) * null_v(:)
     v2x_out(ir,:) = v2x(:) * null_v(:)
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
  !! Implemented:  Perdew86, GGA (PW91), PBE
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
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c_out
  !! derivatives of correlation wr. grho
  !
  ! ... local variables
  !
  INTEGER :: ir
  REAL(DP) :: rho, zeta, grho
  REAL(DP) :: sc, v1c(2), v2c
  REAL(DP), PARAMETER :: small=1.E-10_DP !, epsr=1.E-6_DP
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
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
    IF ( ABS(zeta)>1.0_DP .OR. rho<=small .OR. SQRT(ABS(grho))<=small ) THEN
       sc_out(ir) = 0.d0
       v1c_out(ir,:) = 0.0_DP ; v2c_out(ir) = 0.0_DP
       CYCLE
    ENDIF
    !
    SELECT CASE( igcc_l )
    CASE( 0 )
       !
       sc  = 0.0_DP
       v1c = 0.0_DP
       v2c = 0.0_DP
       !
    CASE( 1 )
       !
       CALL perdew86_spin( rho, zeta, grho, sc, v1c, v2c )
       !
    CASE( 2 )
       !
       CALL ggac_spin( rho, zeta, grho, sc, v1c, v2c )
       !
    CASE( 4 )
       !
       CALL pbec_spin( rho, zeta, grho, 1, sc, v1c, v2c )
       !
    CASE( 8 )
       !
       CALL pbec_spin( rho, zeta, grho, 2, sc, v1c, v2c )
       !
    CASE DEFAULT
       !
       CALL errore( 'xc_gga_drivers (gcc_spin)', 'not implemented', igcc_l )
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
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v2c
  !! derivative of correlation wr. grho
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c_ud
  !! derivative of correlation wr. grho, off-diag. term
  !
  ! ... local variables
  !
  INTEGER :: ir
  REAL(DP) :: rho(2), grho(2)
  REAL(DP) :: grho_ud
  REAL(DP), PARAMETER :: small=1.E-20_DP
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
#endif    
  !
  sc  = 0.0_DP
  v1c = 0.0_DP
  v2c = 0.0_DP
  v2c_ud = 0.0_DP
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
    IF ( rho(1)+rho(2) < small ) THEN
       sc(ir) = 0.0_DP
       v1c(ir,:) = 0.0_DP
       v2c(ir,:) = 0.0_DP ; v2c_ud(ir) = 0.0_DP
       CYCLE
    ENDIF
    !
    CALL lsd_glyp( rho, grho, grho_ud, sc(ir), v1c(ir,:), v2c(ir,:), v2c_ud(ir) )
    !
    SELECT CASE( igcc_l )
    CASE( 3 )
       !
       ! ... void
       !
    CASE( 7 )
       !
       IF ( exx_started_g ) THEN
          sc(ir) = 0.81_DP * sc(ir)
          v1c(ir,:) = 0.81_DP * v1c(ir,:)
          v2c(ir,:) = 0.81_DP * v2c(ir,:)
          v2c_ud(ir) = 0.81_DP * v2c_ud(ir)
       ENDIF
       !
    CASE( 13 )
       !
       IF ( exx_started_g ) THEN
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
