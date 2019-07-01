!
! Copyright (C) 2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
PROGRAM benchmark_libxc
  !
  !------------------------------------------------------------------------------------!
  !  To be run on a single processor
  !------------------------------------------------------------------------------------!
  !
#if defined(__LIBXC)
  !
  USE xc_f90_types_m
  USE xc_f90_lib_m
  !
  USE xc_lda_lsda
  USE xc_gga
  !
  IMPLICIT NONE
  !
  !-------- Common vars ----------------------
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
  INTEGER, PARAMETER :: nnr = 6
  CHARACTER(LEN=120) :: aprx, e_q, f_q
  INTEGER :: ii, ns, np, quit, i_sub, family
  REAL(DP) :: exx_frctn
  LOGICAL :: LDA, GGA, POLARIZED, ENERGY_ONLY, DF_OK
  REAL(DP), PARAMETER :: null = 0.0_DP, pi34 = 0.6203504908994_DP
  !
  !----------QE vars --------------------------
  INTEGER :: iexch_qe, icorr_qe
  REAL(DP) :: rs(nnr), ec_qe2(nnr), vc_qe2(nnr,2)
  REAL(DP), ALLOCATABLE :: rho_qe(:,:)
  REAL(DP), ALLOCATABLE :: rho_tot(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grho(:,:,:), grho_ud(:), grho2(:,:), grh2(:)
  REAL(DP), ALLOCATABLE :: ex_qe(:), ec_qe(:)
  REAL(DP), ALLOCATABLE :: vx_qe(:,:), vc_qe(:,:)
  REAL(DP), ALLOCATABLE :: dmuxc(:,:,:)
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:), v2c_ud(:)
  REAL(DP), ALLOCATABLE :: vrrx(:,:), vsrx(:,:), vssx(:,:)
  REAL(DP), ALLOCATABLE :: vrrc(:,:), vsrc(:,:), vssc(:), vrzc(:,:)
  !
  !--------- LIBXC vars -----------------------
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info1, xc_info2, xc_info3, xc_info4
  CHARACTER(LEN=120) :: name1, name2
  INTEGER :: iexch_lxc, icorr_lxc
  INTEGER :: pol_unpol
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: sigma(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_lxc(:), vc_lxc(:)
  REAL(DP), ALLOCATABLE :: dex_lxc(:), dcr_lxc(:), df_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_rho(:), vc_rho(:)
  REAL(DP), ALLOCATABLE :: vx_sigma(:), vc_sigma(:)
  REAL(DP), ALLOCATABLE :: vx_lxc2(:), vc_lxc2(:)
  REAL(DP), ALLOCATABLE :: ex_lxc2(:), ec_lxc2(:)
  REAL(DP), ALLOCATABLE :: v2rho2_x(:), v2rhosigma_x(:), v2sigma2_x(:)
  REAL(DP), ALLOCATABLE :: v2rho2_c(:), v2rhosigma_c(:), v2sigma2_c(:)
  !
  !
  ! *******************************************************************************
  ! *-----------------------------------------------------------------------------*
  ! * libxc funct. indexes: http://bigdft.org/Wiki/index.php?title=XC_codes       *
  ! * qe      "       "   : see comments in Modules/funct.f90                     *
  ! *-----------------------------------------------------------------------------*
  ! *                                                                             *
  ! *  ... some examples:                                                         *
  ! *                                                                             *
  ! *                         LDA                                 GGA             *       
  ! *              |   q-e     |   libxc  |              |   q-e   |    libxc   | *
  ! *              |___________|__________|              |_________|____________| *
  ! *  slater (x)  |    1      |     1    |  Becke88 (x) |    1    |     106    | *
  ! *  pz (c)      |    1      |     9    |  PW86(c)-POL |    1    |     132    | *
  ! *              |           |          |  PW86(c)-UNP |   21    |      "     | *
  ! *  wigner (c)  |    5      |     2    |  PBE (c)     |    4    |     130    | *
  ! *  vwn (c)     |    2      |     7    |  LYP (c)     |    3    |     131    | *
  ! *  pw (c)      |    4      |    12    |  PW91 (c)    |    2    |     134    | *
  ! *  ...         |   ...     |    ...   |  ...         |   ...   |     ...    | *
  ! *                                                                             *
  ! *******************************************************************************
  !
  !
  !*******************
  !  compatibility for GGA:
  !
  !  - lyp =>    qe:     ec = ec_lyp*rho + ec_glyp  / vc = vc_lyp + v1c_glyp
  !              libxc:  ec = ec_glyp (icorr=131)   / vc = vc_glyp (131)
  !                     ... same for polarized case
  !
  !  - pbe =>  differences of the order of 1-5% only in the polarized case by
  !            adding thw pw lda part (in qe)
  !
  !
  PRINT *, CHAR(10)//" --- BENCHMARK TEST BETWEEN QE AND LIBXC ---"//CHAR(10)//" "
  !
  WRITE (*,'(/,1x,a)', ADVANCE='no') "Derivative of xc?(y/n) "
  READ(*,*) f_q
  DF_OK = .FALSE.
  IF ( TRIM(f_q) == 'y' ) DF_OK = .TRUE.
  IF ( TRIM(f_q) /= 'y' .AND. TRIM(f_q) /= 'n' ) THEN
     PRINT *, CHAR(10)//"ERROR: it is yes (y) or no (n)"//CHAR(10)
     GO TO 10
  ENDIF
  !
  !
  IF ( .NOT.DF_OK ) THEN
    WRITE (*,'(/,1x,a)', ADVANCE='no') "Energy only (y/n)? "
    READ(*,*) e_q
    ENERGY_ONLY = .FALSE.
    IF ( TRIM(e_q) == 'y' ) ENERGY_ONLY = .TRUE.
    IF ( TRIM(e_q) /= 'y' .AND. TRIM(e_q) /= 'n' ) THEN
       PRINT *, CHAR(10)//"ERROR: it is yes (y) or no (n)"//CHAR(10)
       GO TO 10
    ENDIF
  ENDIF
  !
  WRITE (*,'(/,1x,a)', ADVANCE='no') "lda or gga ?  "
  READ(*,*) aprx
  IF ( TRIM(aprx) /= 'lda' .AND. TRIM(aprx) /= 'gga' ) THEN
     PRINT *, CHAR(10)//"ERROR: you can only choose lda or gga"//CHAR(10)
     GO TO 10
  ENDIF
  WRITE (*,'(/,1x,a)', ADVANCE='no') "Polarization switch (1 unpolarized,  & 
                                                         & 2 polarized):  "
  READ(*,*) ns
  IF ( ns/=1 .AND. ns/=2 ) THEN
     PRINT *, CHAR(10)//"ERROR: you can only choose 1 or 2"//CHAR(10)
     GO TO 10
  ENDIF
  WRITE (*,'(/,1x,a)') "-- Functional indexes "
  WRITE (*,'(/,1x,a)', ADVANCE='no') "iexch_libxc  icorr_libxc: "
  READ(*,*) iexch_lxc, icorr_lxc
  WRITE (*,'(/,1x,a)', ADVANCE='no') "iexch_qe  icorr_qe: "
  READ(*,*) iexch_qe, icorr_qe
  IF (ns == 2 .AND. icorr_qe/=0 .AND. icorr_qe/=1 .AND. icorr_qe/=2 .AND. &
                    icorr_qe/=4 .AND. icorr_qe/=8 .AND. icorr_qe/=3 .AND. &
                    icorr_qe/=7 .AND. icorr_qe/=13) THEN
     PRINT *, CHAR(10)//" ERROR: icorr_qe not available at these conditions"//CHAR(10)
     GO TO 10
  ENDIF
  !
  !
  IF ( TRIM(aprx) == 'lda' ) THEN
     LDA = .TRUE.
     GGA = .FALSE.
  ELSE 
     LDA = .FALSE.
     GGA = .TRUE.
  ENDIF
  !
  POLARIZED = .FALSE.
  IF (ns == 2) THEN
     POLARIZED = .TRUE.
  ENDIF
  !
  pol_unpol = XC_UNPOLARIZED
  np = 1
  !
  IF ( ns == 2 ) THEN
     pol_unpol = XC_POLARIZED
     np = 3
  ENDIF
  !
  !
  ! *******************************************************************************
  ! *                     Polarized case:                                         *
  ! *                                                                             *
  ! *  qe =>  rho_qe(:,1) = up    |      libxc     =>        rho_lxc(2n+1) = up   *
  ! *         rho_qe(:,2) = down  |            (dim=2*nnr)   rho_lxc(2n+2) = down *
  ! *                             |                                               *
  ! *         grho(:,1)  = uu     |            (dim=3*nnr)   sigma(3n+1) = uu     *
  ! *         grho(:,2)  = dd     |                          sigma(3n+2) = ud     *
  ! *         grho_ud(:) = ud     |                          sigma(3n+3) = dd     *
  ! *                                                                             *
  ! *******************************************************************************
  !
  !
  ! ------ Allocations --------
  !
  ! ... qe
  !
  ALLOCATE( rho_qe(nnr,ns) )
  IF ( POLARIZED ) ALLOCATE( rho_tot(nnr), zeta(nnr) )
  ALLOCATE( ex_qe(nnr), ec_qe(nnr) )
  IF ( LDA ) THEN
     ALLOCATE( vx_qe(nnr,ns), vc_qe(nnr,ns) )
     IF ( DF_OK ) THEN
       IF ( .NOT.POLARIZED ) ALLOCATE( dmuxc(nnr,1,1) )
       IF ( POLARIZED ) ALLOCATE( dmuxc(nnr,2,2) )
     ENDIF
  ELSEIF ( GGA ) THEN
     ALLOCATE( grho(nnr,3,ns), grho2(nnr,ns), grho_ud(nnr), grh2(nnr) )
     ALLOCATE( v1x(nnr,ns), v2x(nnr,ns) )
     ALLOCATE( v1c(nnr,ns), v2c(nnr,ns), v2c_ud(nnr) )
     IF (DF_OK) THEN
       ALLOCATE( vrrx(nnr,ns), vsrx(nnr,ns), vssx(nnr,ns) )
       ALLOCATE( vrrc(nnr,ns), vsrc(nnr,ns), vssc(nnr), vrzc(nnr,ns) )
     ENDIF
  ENDIF
  !
  ! ... libxc
  !
  ALLOCATE( rho_lxc(nnr*ns) )
  ALLOCATE( ex_lxc(nnr), ec_lxc(nnr) )
  IF ( LDA ) THEN
    ALLOCATE( vx_lxc(nnr*ns), vc_lxc(nnr*ns) )
    IF ( DF_OK ) THEN
      IF ( .NOT.POLARIZED ) ALLOCATE( dex_lxc(nnr), dcr_lxc(nnr), df_lxc(nnr) )
      IF ( POLARIZED ) ALLOCATE( dex_lxc(nnr*3), dcr_lxc(nnr*3), df_lxc(nnr*3) )
    ENDIF
  ELSEIF ( GGA ) THEN
    ALLOCATE( sigma(nnr*np) )
    ALLOCATE( vx_rho(nnr*ns), vx_sigma(nnr*np) )
    ALLOCATE( vc_rho(nnr*ns), vc_sigma(nnr*np) )
    ALLOCATE( vx_lxc2(nnr*ns), vc_lxc2(nnr*ns) )
    ALLOCATE( ex_lxc2(nnr), ec_lxc2(nnr) )
    IF ( DF_OK ) THEN
      IF ( .NOT.POLARIZED ) THEN
        ALLOCATE( dex_lxc(nnr), dcr_lxc(nnr), df_lxc(nnr) )
        ALLOCATE( v2rho2_x(nnr), v2rhosigma_x(nnr), v2sigma2_x(nnr) )
        ALLOCATE( v2rho2_c(nnr), v2rhosigma_c(nnr), v2sigma2_c(nnr) )
      ELSEIF ( POLARIZED ) THEN
        ALLOCATE( dex_lxc(nnr*3), dcr_lxc(nnr*3), df_lxc(nnr*3) )
        ALLOCATE( v2rho2_x(3*nnr), v2rhosigma_x(6*nnr), v2sigma2_x(6*nnr) )
        ALLOCATE( v2rho2_c(3*nnr), v2rhosigma_c(6*nnr), v2sigma2_c(6*nnr) )
      ENDIF
    ENDIF
  ENDIF
  !
  !
  !----- Initializations -----
  !
  ! ... qe
  !
  rho_qe = 0.0_DP
  IF ( POLARIZED ) THEN
     rho_tot = 0.0_DP
     zeta = 0.0_DP
  ENDIF
  IF ( GGA ) THEN
     grho = 0.0_DP
     grho2 = 0.0_DP
     grho_ud = 0.0_DP
  ENDIF
  !
  ! ... libcx
  !
  rho_lxc = 0.0_DP
  IF ( GGA ) sigma = 0.0_DP
  !
  ! -------- Setting up an arbitrary input for both qe and libxc -----
  !
  ! ... qe
  !
  DO ii = 1, nnr
     !
     rho_qe(ii,1) = DBLE(ii)/DBLE(nnr+2)
     !
     IF ( GGA ) THEN
        grho(ii,1,1) = ABS( 0.05_DP + 0.8_DP*SIN(DBLE(ii)) )
        grho(ii,2,1) = ABS( 0.05_DP + 0.7_DP*SIN(DBLE(ii)) )
        grho(ii,3,1) = ABS( 0.05_DP + 0.6_DP*SIN(DBLE(ii)) )
     ENDIF
     grho2(ii,1) = grho(ii,1,1)**2 + grho(ii,2,1)**2 + grho(ii,3,1)**2
     !
     IF ( POLARIZED ) THEN
        !
        rho_qe(ii,2) = (1.0_DP - rho_qe(ii,1))*0.7_DP
        rho_tot(ii) = rho_qe(ii,1) + rho_qe(ii,2)
        zeta(ii) = (rho_qe(ii,1) - rho_qe(ii,2)) / rho_tot(ii)
        !
        IF ( GGA ) THEN
           grho(ii,1,2) = ABS( (1.0_DP - grho(ii,1,1))*0.7_DP )
           grho(ii,2,2) = ABS( (1.0_DP - grho(ii,2,1))*0.6_DP )
           grho(ii,3,2) = ABS( (1.0_DP - grho(ii,3,1))*0.5_DP )
           !
           grh2(ii)= ( grho(ii,1,1) + grho(ii,1,2) )**2 + &
                     ( grho(ii,2,1) + grho(ii,2,2) )**2 + &
                     ( grho(ii,3,1) + grho(ii,3,2) )**2
           !
           grho2(ii,2) = ( grho(ii,1,2)**2 + grho(ii,2,2)**2 + grho(ii,3,2)**2 )
           !
           grho_ud(ii) = grho(ii,1,1) * grho(ii,1,2) + &
                         grho(ii,2,1) * grho(ii,2,2) + &
                         grho(ii,3,1) * grho(ii,3,2)
           !
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  ! ... libxc
  !
  DO ii = 1, nnr
     !
     IF ( .NOT. POLARIZED ) THEN
        !
        rho_lxc(ii) = rho_qe(ii,1)
        !
        IF ( GGA ) sigma(ii) = grho2(ii,1)
        !
     ELSE
        !
        rho_lxc(2*ii-1) = rho_qe(ii,1)
        rho_lxc(2*ii) = rho_qe(ii,2)
        !
        IF ( GGA ) THEN
          sigma(3*ii-2) = grho2(ii,1)
          sigma(3*ii-1) = grho_ud(ii)
          sigma(3*ii) = grho2(ii,2)
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  !
  !-------- Calculation of energies and potential -------------------
  !
  IF ( LDA ) THEN
     !
     !------ LIBXC ------
     !
     CALL xc_f90_func_init( xc_func, xc_info1, iexch_lxc, pol_unpol )   ! ... EXCHANGE
      CALL xc_f90_lda_exc_vxc( xc_func, nnr, rho_lxc(1), ex_lxc(1), vx_lxc(1) )
      !
      IF ( DF_OK ) CALL xc_f90_lda_fxc( xc_func, nnr, rho_lxc(1), dex_lxc(1) )
     CALL xc_f90_func_end( xc_func )
     ! 
     CALL xc_f90_func_init( xc_func, xc_info2, icorr_lxc, pol_unpol )  ! ... CORRELATION
      CALL xc_f90_lda_exc_vxc( xc_func, nnr, rho_lxc(1), ec_lxc(1), vc_lxc(1) )
      !
      IF ( DF_OK ) THEN   
        CALL xc_f90_lda_fxc( xc_func, nnr, rho_lxc(1), dcr_lxc(1) )   
        df_lxc = (dex_lxc + dcr_lxc) * 2.0_DP   
      ENDIF
     CALL xc_f90_func_end( xc_func )
     !
     !----- QE ----------
     !
     CALL select_lda_functionals( iexch_qe, icorr_qe )   ! ... EXCHANGE and CORRELATION
     !
     IF ( DF_OK ) THEN
       IF ( .NOT.POLARIZED ) THEN
         CALL dmxc_lda( nnr, rho_qe(:,1), dmuxc(:,1,1) )
       ELSE
         CALL dmxc_lsda( nnr, rho_qe, dmuxc )
       ENDIF
     ENDIF
     !  
     IF ( .NOT. POLARIZED ) THEN
        CALL xc_lda( nnr, rho_qe(:,1), ex_qe, ec_qe, vx_qe(:,1), vc_qe(:,1) )
     ELSE
        CALL xc_lsda( nnr, rho_tot, zeta, ex_qe, ec_qe, vx_qe, vc_qe )
     ENDIF
     !
     !
  ELSEIF ( GGA ) THEN
     !
     !------ LIBXC ------
     !
     exx_frctn = 0.0_DP
     !
     CALL xc_f90_func_init( xc_func, xc_info1, iexch_lxc, pol_unpol )
      family = xc_f90_info_family( xc_info1 )
      IF (family == XC_FAMILY_HYB_GGA) CALL xc_f90_hyb_exx_coef( xc_func, exx_frctn )
      !
      CALL xc_f90_gga_exc_vxc( xc_func, nnr, rho_lxc(1), sigma(1), ex_lxc(1), vx_rho(1), vx_sigma(1) )
      !
      IF ( DF_OK ) CALL xc_f90_gga_fxc( xc_func, nnr, rho_lxc(1), sigma(1), v2rho2_x(1), v2rhosigma_x(1), v2sigma2_x(1) )
     CALL xc_f90_func_end( xc_func )
     !
     ! remove Slater term for compatibility with QE
     CALL xc_f90_func_init( xc_func, xc_info3, 1, pol_unpol )
      CALL xc_f90_lda_exc_vxc( xc_func, nnr, rho_lxc(1), ex_lxc2(1), vx_lxc2(1) )
      !
      IF ( DF_OK ) CALL xc_f90_lda_fxc( xc_func, nnr, rho_lxc(1), dex_lxc(1) )
     CALL xc_f90_func_end( xc_func )
     !
     IF (DF_OK) THEN
       v2rho2_x = v2rho2_x - dex_lxc
     ENDIF
     !
     !
     DO ii = 1, nnr
       IF ( .NOT. POLARIZED) THEN
          ex_lxc(ii) = (ex_lxc(ii) - ex_lxc2(ii)) * rho_lxc(ii)
       ELSE
          ex_lxc(ii) = (ex_lxc(ii) - ex_lxc2(ii)) * rho_tot(ii)
       ENDIF
     ENDDO
     vx_rho = vx_rho - vx_lxc2
     !
     ! ----
     !
     CALL xc_f90_func_init( xc_func, xc_info2, icorr_lxc, pol_unpol )
      CALL xc_f90_gga_exc_vxc( xc_func, nnr, rho_lxc(1), sigma(1), ec_lxc(1), vc_rho(1), vc_sigma(1) )
      !
      IF ( DF_OK ) CALL xc_f90_gga_fxc( xc_func, nnr, rho_lxc(1), sigma(1), v2rho2_c(1), v2rhosigma_c(1), v2sigma2_c(1) )
     CALL xc_f90_func_end( xc_func )
     !
     ec_lxc2 = 0.d0
     vc_lxc2 = 0.d0
     IF (icorr_lxc /= 131 ) then ! .AND. .NOT.(icorr_lxc == 130 .AND. POLARIZED)) THEN 
        ! remove LDA correlation for compatibility with QE
        i_sub=12 !(pw)
        !IF (icorr_lxc == 132) i_sub=9 !(pz)
        CALL xc_f90_func_init( xc_func, xc_info4, i_sub, pol_unpol )
         CALL xc_f90_lda_exc_vxc( xc_func, nnr, rho_lxc(1), ec_lxc2(1), vc_lxc2(1) )
         !
         IF ( DF_OK ) CALL xc_f90_lda_fxc( xc_func, nnr, rho_lxc(1), dex_lxc(1) )
        CALL xc_f90_func_end( xc_func )
        !
        IF (DF_OK) THEN
          v2rho2_c = v2rho2_c - dex_lxc
        ENDIF
        !
     ENDIF
     !
     !
     DO ii = 1, nnr
        IF ( .NOT. POLARIZED) THEN
           ec_lxc(ii) = (ec_lxc(ii) - ec_lxc2(ii)) * rho_lxc(ii)
        ELSE
           ec_lxc(ii) = (ec_lxc(ii) - ec_lxc2(ii)) * rho_tot(ii)
        ENDIF
     ENDDO
     vc_rho(:) = vc_rho(:) - vc_lxc2(:)
     !
     vx_sigma = vx_sigma*2.0_DP
     vc_sigma = vc_sigma*2.0_DP
     IF ( POLARIZED ) THEN
        DO ii = 1, nnr
           vc_sigma(3*ii-1) = vc_sigma(3*ii-1) / 2.0_DP
        ENDDO
     ENDIF
     !
     !
     !----- QE ----------
     !
     CALL select_gga_functionals( iexch_qe, icorr_qe, exx_fraction=exx_frctn )
     !
     IF ( DF_OK ) THEN
        !
        !CALL select_lda_functionals( iexch_qe, icorr_qe )
        !
        IF (.NOT. POLARIZED) THEN
           CALL dgcxc( nnr, rho_qe(:,1), grho2(:,1), vrrx(:,1), vsrx(:,1), vssx(:,1), &
                       vrrc(:,1), vsrc(:,1), vssc )
           !
        ELSE
           CALL dgcxc_spin( nnr, rho_qe, grho, vrrx, vsrx, vssx, vrrc, vsrc, &
                            vssc, vrzc )
        ENDIF
        !
     ELSE
        ! 
        IF ( .NOT. POLARIZED ) THEN
          !
          CALL gcxc( nnr, rho_qe(:,1), grho2(:,1), ex_qe, ec_qe, v1x(:,1), v2x(:,1), v1c(:,1), v2c(:,1) )
          !
          IF ( icorr_qe == 3 ) THEN     ! LDA part of Lee-Yang-Parr is not available in libxc.
            rs(:) = pi34 / rho_qe(:,1)**(1.d0/3.d0)
            CALL lyp( nnr, rs, ec_qe2, vc_qe2(:,1) )
            ec_qe(:) = ec_qe(:) + ec_qe2(:)*rho_qe(:,1)
            v1c(:,1) = v1c(:,1) + vc_qe2(:,1)
          ENDIF
         !
       ELSE
         !
         CALL gcx_spin( nnr, rho_qe, grho2, ex_qe, v1x, v2x )
         !
         IF (icorr_qe/=3 .AND. icorr_qe/=7 .AND. icorr_qe/=13 ) THEN
           CALL gcc_spin( nnr, rho_tot, zeta, grh2, ec_qe, v1c, v2c(:,1) )
           v2c(:,2) = v2c(:,1)
           v2c_ud(:) = v2c(:,1)
           !
           !
           !IF (icorr_qe == 4) THEN
           !   rs(:) = pi34 / rho_tot(:)**(1.d0/3.d0)   !  .... trovare le combinazioni e mettere ordine
           !   CALL pw_spin( nnr, rs, zeta, ec_qe2, vc_qe2 )
           !   ec_qe(:) = ec_qe(:) + ec_qe2(:)*rho_tot(:)
           !   v1c(:,1) = v1c(:,1) + vc_qe2(:,1)
           !   v1c(:,2) = v1c(:,2) + vc_qe2(:,2)
           !ENDIF
           !
         ELSE
           CALL gcc_spin_more( nnr, rho_qe, grho2, grho_ud, ec_qe, v1c, v2c, v2c_ud )
           CALL lsd_lyp( nnr, rho_tot, zeta, ec_qe2, vc_qe2 )
           ec_qe(:) = ec_qe(:) + ec_qe2(:)*rho_tot(:)
           v1c(:,1) = v1c(:,1) + vc_qe2(:,1)
           v1c(:,2) = v1c(:,2) + vc_qe2(:,2)
         ENDIF
         !
       ENDIF
       !
       !
     ENDIF
     !
     !ec_qe=ec_qe*2.d0
     !v1c = v1c * 2.d0
     !
  ENDIF 
  !
  !------------------
  !
  CALL xc_f90_info_name( xc_info1, name1 )
  CALL xc_f90_info_name( xc_info2, name2 )
  !
  PRINT *, "=================================== "//CHAR(10)//" "
  PRINT *, "libxc functionals:"//CHAR(10)//" "
  PRINT *, "Exchange: ", TRIM(name1)
  PRINT *, "Correlation: ", TRIM(name2)
  PRINT *, " "
  !
  IF ( LDA ) THEN
     !
     DO ii = 1, nnr, nnr-1
        WRITE(*,909) ii, nnr
        IF ( .NOT. POLARIZED ) THEN
           WRITE (*, 401 ) rho_qe(ii,1)
        ELSE
           WRITE (*, 402 ) rho_qe(ii,1), rho_qe(ii,2)
        ENDIF
        PRINT *, " "
        IF (.NOT. DF_OK) THEN  
          PRINT *, "=== Exchange and correlation energies: ==="
          WRITE (*,102) ex_qe(ii),  ec_qe(ii)
          WRITE (*,202) ex_lxc(ii), ec_lxc(ii)
          PRINT *, " --- "
          WRITE (*,302) ex_qe(ii)-ex_lxc(ii), ec_qe(ii)-ec_lxc(ii)
          !
          IF (.NOT. ENERGY_ONLY) THEN
             PRINT *, " "
             PRINT *, "=== Exchange potential ==="
             IF ( .NOT. POLARIZED ) THEN
                WRITE (*,101) vx_qe(ii,1)
                WRITE (*,201) vx_lxc(ii)
                PRINT *, " --- "
                WRITE (*,301) vx_qe(ii,1)-vx_lxc(ii)
             ELSEIF ( POLARIZED ) THEN
                WRITE (*,102) vx_qe(ii,1), vx_qe(ii,2)
                WRITE (*,202) vx_lxc(2*ii-1), vx_lxc(2*ii)
                PRINT *, " --- "
                WRITE (*,302) vx_qe(ii,1)-vx_lxc(2*ii-1), vx_qe(ii,2)-vx_lxc(2*ii)
             ENDIF
             PRINT *, " "
             PRINT *, "=== Correlation potential ==="
             IF ( .NOT. POLARIZED ) THEN
                WRITE (*,101) vc_qe(ii,1)
                WRITE (*,201) vc_lxc(ii)
                PRINT *, " --- "
                WRITE (*,301) vc_qe(ii,1)-vc_lxc(ii)
             ELSEIF ( POLARIZED ) THEN
                WRITE (*,102) vc_qe(ii,1), vc_qe(ii,2)
                WRITE (*,202) vc_lxc(2*ii-1), vc_lxc(2*ii)
                PRINT *, " --- "
                WRITE (*,302) vc_qe(ii,1)-vc_lxc(2*ii-1), vc_qe(ii,2)-vc_lxc(2*ii)
             ENDIF
             PRINT *, "--- ---"
          ENDIF
          !
        ELSE
          !
          PRINT *, " "   
          PRINT *, "=== First derivative of xc functional (exch+corr): ==="   
          IF ( .NOT. POLARIZED ) THEN   
            WRITE (*,101) dmuxc(ii,1,1)   
            WRITE (*,201) df_lxc(ii)   
            PRINT *, " --- "   
            WRITE (*,301) dmuxc(ii,1,1)-df_lxc(ii)   
          ELSE   
            WRITE (*,104) dmuxc(ii,1,1), dmuxc(ii,2,1), dmuxc(ii,2,2), dmuxc(ii,1,2)   
            WRITE (*,203) df_lxc(3*ii-2), df_lxc(3*ii-1), df_lxc(3*ii)   
            PRINT *, " --- "   
            WRITE (*,303) dmuxc(ii,1,1)-df_lxc(3*ii-2), dmuxc(ii,2,1)-df_lxc(3*ii-1), &   
                          dmuxc(ii,2,2)-df_lxc(3*ii)   
          ENDIF   
        ENDIF
        !
     ENDDO
     !
  ELSEIF ( GGA ) THEN
     !
     DO ii = 1, nnr !, nnr-1
        WRITE(*,*) ' '
        WRITE(*,*) ' '
        WRITE(*,909) ii, nnr
        IF (.NOT. POLARIZED ) THEN
           WRITE (*,401) rho_qe(ii,1)
           WRITE (*,501) grho2(ii,1)
        ELSE
           WRITE (*,402) rho_qe(ii,1), rho_qe(ii,2)
           WRITE (*,503) grho2(ii,1), grho_ud(ii), grho2(ii,2)
        ENDIF
        PRINT *, " "
        IF (.NOT. DF_OK) THEN
          PRINT *, "=== Exchange and correlation energies: ==="
          WRITE (*,102) ex_qe(ii),  ec_qe(ii)
          WRITE (*,202) ex_lxc(ii), ec_lxc(ii)
          PRINT *, " --- "
          WRITE (*,302) ex_qe(ii)-ex_lxc(ii), ec_qe(ii)-ec_lxc(ii)
          !WRITE (*,302) ex_qe(ii)/ex_lxc(ii), ec_qe(ii)/ec_lxc(ii)
          !
          IF (.NOT. ENERGY_ONLY) THEN
             !
             PRINT *, " "
             PRINT *, "=== Exchange potential vrho ==="
             IF ( .NOT. POLARIZED ) THEN
                WRITE (*,101) v1x(ii,1)
                WRITE (*,201) vx_rho(ii)
                PRINT *, " --- "
                WRITE (*,301) v1x(ii,1)-vx_rho(ii)
             ELSEIF ( POLARIZED ) THEN
                WRITE (*,102) v1x(ii,1), v1x(ii,2)
                WRITE (*,202) vx_rho(2*ii-1), vx_rho(2*ii)
                PRINT *, " --- "
                WRITE (*,302) v1x(ii,1)-vx_rho(2*ii-1), v1x(ii,2)-vx_rho(2*ii)
             ENDIF        
             !
             PRINT *, " "
             PRINT *, "=== Exchange potential vsigma ==="
             IF ( .NOT. POLARIZED ) THEN
                WRITE (*,101) v2x(ii,1)
                WRITE (*,201) vx_sigma(ii)
                PRINT *, " --- "
                WRITE (*,301) v2x(ii,1)-vx_sigma(ii)
             ELSEIF ( POLARIZED ) THEN
                WRITE (*,103) v2x(ii,1), null, v2x(ii,2)
                WRITE (*,203) vx_sigma(3*ii-2), vx_sigma(3*ii-1), vx_sigma(3*ii)
                PRINT *, " --- "
                WRITE (*,303) v2x(ii,1)-vx_sigma(3*ii-2), null-vx_sigma(3*ii-1), &
                                                       v2x(ii,2)-vx_sigma(3*ii)
             ENDIF   
             !
             PRINT *, " "
             PRINT *, "=== Correlation potential vrho ==="
             IF ( .NOT. POLARIZED ) THEN
                WRITE (*,101) v1c(ii,1)
                WRITE (*,201) vc_rho(ii)
                WRITE (*,301) v1c(ii,1)-vc_rho(ii)
             ELSEIF ( POLARIZED ) THEN
                WRITE (*,102) v1c(ii,1), v1c(ii,2)
                WRITE (*,202) vc_rho(2*ii-1), vc_rho(2*ii)
                PRINT *, " --- "
                WRITE (*,302) v1c(ii,1)-vc_rho(2*ii-1), v1c(ii,2)-vc_rho(2*ii)
             ENDIF        
             !
             PRINT *, " "
             PRINT *, "=== Correlation potential vsigma ==="
             IF ( .NOT. POLARIZED ) THEN
                WRITE (*,101) v2c(ii,1)
                WRITE (*,201) vc_sigma(ii)
                PRINT *, " --- "
                WRITE (*,301) v2c(ii,1)-vc_sigma(ii)
             ELSEIF ( POLARIZED ) THEN
                WRITE (*,103) v2c(ii,1), v2c_ud(ii), v2c(ii,2)
                WRITE (*,203) vc_sigma(3*ii-2), vc_sigma(3*ii-1), vc_sigma(3*ii)
                PRINT *, " --- "
                WRITE (*,303) v2c(ii,1)-vc_sigma(3*ii-2), v2c_ud(ii)-vc_sigma(3*ii-1), &
                                                          v2c(ii,2)-vc_sigma(3*ii)
             ENDIF
             !
          ENDIF
        ELSE
          PRINT *, " "   
          PRINT *, "====== First derivative of xc functional: ==="  
          !
          PRINT *, " "
          PRINT *, "=== Exchange part ==="
          !
          IF ( .NOT. POLARIZED ) THEN
            WRITE (*,103) vrrx(ii,1), vsrx(ii,1), vssx(ii,1)
            WRITE (*,203) v2rho2_x(ii), v2rhosigma_x(ii)*2.d0, v2sigma2_x(ii)*4.d0
            PRINT *, " --- "
            WRITE (*,303) vrrx(ii,1)-v2rho2_x(ii), vsrx(ii,1)-v2rhosigma_x(ii)*2.d0, vssx(ii,1)-v2sigma2_x(ii)*4.d0
          ELSE   
            PRINT *, 'vrrx'
            WRITE (*,102) vrrx(ii,1), vrrx(ii,2)
            WRITE (*,202) v2rho2_x(3*ii-2), v2rho2_x(3*ii)
            PRINT *, " --- "
            WRITE (*,302) vrrx(ii,1)-v2rho2_x(3*ii-2), vrrx(ii,2)-v2rho2_x(3*ii)
            PRINT *, 'vsrx'
            WRITE (*,102) vsrx(ii,1), vsrx(ii,2)
            WRITE (*,202) v2rhosigma_x(6*ii-5)*2.d0, v2rhosigma_x(6*ii)*2.d0
            PRINT *, " --- "
            WRITE (*,303) vsrx(ii,1)-v2rhosigma_x(6*ii-5)*2, vsrx(ii,2)-v2rhosigma_x(6*ii)*2.d0
            PRINT *, 'vssx'
            WRITE (*,102) vssx(ii,1), vssx(ii,2)
            !
            WRITE (*,202) v2sigma2_x(6*ii-5)*4, v2sigma2_x(6*ii)*2
            !WRITE (*,203) v2sigma2_x(6*ii-5)*2, v2sigma2_x(6*ii-4), v2sigma2_x(6*ii-3)
            !WRITE (*,203) v2sigma2_x(6*ii-2), v2sigma2_x(6*ii-1), v2sigma2_x(6*ii)
            PRINT *, " --- "
            WRITE (*,302) vssx(ii,1)-v2sigma2_x(6*ii-5)*4, vssx(ii,2)-v2sigma2_x(6*ii)*4
          ENDIF
          !
          PRINT *, " "
          PRINT *, "=== Corr part ==="
          !
          IF ( .NOT. POLARIZED ) THEN
            WRITE (*,103) vrrc(ii,1), vsrc(ii,1), vssc(ii)
            WRITE (*,203) v2rho2_c(ii), v2rhosigma_c(ii)*2.d0, v2sigma2_c(ii)*4.d0
            PRINT *, " --- "
            WRITE (*,303) vrrc(ii,1)-v2rho2_c(ii), vsrc(ii,1)-v2rhosigma_c(ii)*2.d0, vssc(ii)-v2sigma2_c(ii)*4.d0
          ELSE
            PRINT *, 'vrrc'
            WRITE (*,102) vrrc(ii,1), vrrc(ii,2)
            WRITE (*,203) v2rho2_c(3*ii-1), v2rho2_c(3*ii-2), v2rho2_c(3*ii)
            PRINT *, " --- "
            !WRITE (*,302) vrrc(ii,1)-v2rho2_c(3*ii-2), vrrc(ii,2)-v2rho2_c(3*ii)
            PRINT *, 'vsrc'
            WRITE (*,102) vsrc(ii,1)+4.d0*vsrc(ii,2), vsrc(ii,2)
            WRITE (*,202) v2rhosigma_c(6*ii-5)*2.d0, v2rhosigma_c(6*ii)*2.d0
            !WRITE (*,203) v2rhosigma_c(6*ii-5), v2rhosigma_c(6*ii-4), v2rhosigma_c(6*ii-3)
            !WRITE (*,203) v2rhosigma_c(6*ii-2), v2rhosigma_c(6*ii-1), v2rhosigma_c(6*ii)
            PRINT *, " --- "
            WRITE (*,303) vsrc(ii,1)-v2rhosigma_c(6*ii-5)*2, vsrc(ii,2)-v2rhosigma_c(6*ii)*2.d0
            PRINT *, 'vssc'
            WRITE (*,103) vssc(ii), vrzc(ii,1), vrzc(ii,2)
            !WRITE (*,202) v2sigma2_c(6*ii-5)*4, v2sigma2_c(6*ii)*2
            !
            WRITE (*,203) v2sigma2_c(6*ii-5), v2sigma2_c(6*ii-4), v2sigma2_c(6*ii-3)
            WRITE (*,203) v2sigma2_c(6*ii-2), v2sigma2_c(6*ii-1), v2sigma2_c(6*ii)
            PRINT *, " --- "
            !WRITE (*,302) vssc(ii,1)-v2sigma2_c(6*ii-5)*4, vssc(ii,2)-v2sigma2_c(6*ii)*4
            !
          ENDIF
          !
        ENDIF
        !
     ENDDO
     !
  ENDIF
  !
  101 FORMAT('qe: ',3x,F17.14)
  102 FORMAT('qe: ',3x,F17.14,4x,F17.14)
  103 FORMAT('qe: ',3x,F17.14,4x,F17.14,4x,F17.14)
  104 FORMAT('qe: ',3x,F17.14,4x,F17.14,4x,F17.14,4x,F17.14)
  !
  201 FORMAT('libxc: ',F17.14)
  202 FORMAT('libxc: ',F17.14,4x,F17.14)
  203 FORMAT('libxc: ',F17.14,4x,F17.14,4x,F17.14)
  !
  301 FORMAT('diffs: ',F17.14)
  302 FORMAT('diffs: ',F17.14,4x,F17.14)
  303 FORMAT('diffs: ',F17.14,4x,F17.14,4x,F17.14)
  !
  401 FORMAT('rho: ',F17.14)
  402 FORMAT('rho(up,down): ',F17.14,4x,F17.14)
  !
  501 FORMAT('grho2: ',F17.14)
  503 FORMAT('grho2(uu,ud,dd): ',F17.14,4x,F17.14,4x,F17.14)
  !
  909 FORMAT('grid-point: ',I4,' of',I4)
  !
  DEALLOCATE( rho_qe )
  IF ( POLARIZED ) DEALLOCATE( rho_tot, zeta )
  DEALLOCATE( ex_qe, ec_qe )
  IF ( LDA ) THEN
     DEALLOCATE( vx_qe, vc_qe )
     IF ( DF_OK ) THEN
       DEALLOCATE( dmuxc )
     ENDIF
  ELSEIF ( GGA ) THEN
     DEALLOCATE( grho, grho2, grho_ud, grh2 )
     DEALLOCATE( v1x, v2x )
     DEALLOCATE( v1c, v2c, v2c_ud )
     IF (DF_OK) THEN
       DEALLOCATE( vrrx, vsrx, vssx )
       DEALLOCATE( vrrc, vsrc, vssc, vrzc )
     ENDIF
  ENDIF
  !
  DEALLOCATE( rho_lxc )
  DEALLOCATE( ex_lxc, ec_lxc )
  IF ( LDA ) THEN
    DEALLOCATE( vx_lxc, vc_lxc )
    IF ( DF_OK ) THEN
      DEALLOCATE( dex_lxc, dcr_lxc, df_lxc )
    ENDIF
  ELSEIF ( GGA ) THEN
    DEALLOCATE( sigma )
    DEALLOCATE( vx_rho, vx_sigma )
    DEALLOCATE( vc_rho, vc_sigma )
    DEALLOCATE( vx_lxc2, vc_lxc2 )
    DEALLOCATE( ex_lxc2, ec_lxc2 )
    IF (DF_OK) THEN
      DEALLOCATE( v2rho2_x, v2rhosigma_x, v2sigma2_x )
      DEALLOCATE( v2rho2_c, v2rhosigma_c, v2sigma2_c )
    ENDIF
  ENDIF
  !
  PRINT *, " "
  !
#else
  !
  PRINT *, "ERROR: library libxc not included."
  !
#endif
10 STOP
  !
END PROGRAM benchmark_libxc
