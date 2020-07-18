!
! Copyright (C) 2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------------
PROGRAM xctest_qe_libxc
  !--------------------------------------------------------------------------------
  !! This program tests the output results (energies and potentials) of the 
  !! xc-driver routines by comparing the results obtained with q-e internal library
  !! of xc functionals with those obtained with libxc library.
  !! It also tests the xc-derivative part (mainly used in PHonon) by comparing q-e
  !! derivative of xc (for LDA and GGA) with the libxc one.
  !
  !! * LDA ;
  !! * derivative of LDA (dmxc) ;
  !! * GGA ;
  !! * derivative of GGA (dgcxc) ;
  !! * metaGGA.
  !
  !------------------------------------------------------------------------------------!
  !  To be run on a single processor
  !------------------------------------------------------------------------------------!
  !
#if defined(__LIBXC)
  !
  USE xc_f03_lib_m
  !
  USE kinds,          ONLY: DP
  USE funct,          ONLY: set_dft_from_name, set_exx_fraction
  USE funct,          ONLY: get_iexch, get_icorr, get_igcx, get_igcc, &
                            get_meta, get_metac, reset_dft
  USE corr_lda,       ONLY: lyp, lsd_lyp
  USE xc_lda_lsda,    ONLY: xc
  USE xc_gga,         ONLY: xc_gcx
  USE xc_mgga,        ONLY: xc_metagcx
  !
  IMPLICIT NONE
  !
  !-------- Common vars ----------------------
  INTEGER :: nnr
  CHARACTER(LEN=120) :: aprx, e_q, f_q
  INTEGER :: ii, ns, np, ipol, quit, i_sub, family, ithr, nthr
  REAL(DP) :: exx_frctn
  LOGICAL :: LDA, GGA, MGGA, POLARIZED, ENERGY_ONLY, DF_OK
  !LOGICAL :: is_it_out, is_dit_out
  REAL(DP), PARAMETER :: null=0.0_DP, pi34=0.6203504908994_DP
  REAL(DP), PARAMETER :: thresh_lda  = 0.d0, & !1.E-6_DP, &
                         thresh_gga  = 0.d0, & !1.E-6_DP, &
                         thresh_mgga = 0.d0 !1.E-6_DP
  !
  INTEGER :: iexch_qe, icorr_qe, igcx_qe, igcc_qe, imeta_qe, imetac_qe
  INTEGER :: iexch_lxc, icorr_lxc, igcx_lxc, igcc_lxc, imeta_lxc, imetac_lxc
  !
  !----------QE vars --------------------------
  REAL(DP), ALLOCATABLE :: rho(:,:), rho_tz(:,:)
  REAL(DP), ALLOCATABLE :: grho(:,:,:), grh(:,:,:)
  REAL(DP), ALLOCATABLE :: ex_qe(:), ec_qe(:)
  REAL(DP), ALLOCATABLE :: exg_qe(:), ecg_qe(:)
  REAL(DP), ALLOCATABLE :: vx_qe(:,:), vc_qe(:,:)
  REAL(DP), ALLOCATABLE :: dmuxc_qe(:,:,:)
  REAL(DP), ALLOCATABLE :: dvxcrr_qe(:,:,:), dvxcsr_qe(:,:,:), &
                           dvxcss_qe(:,:,:)
  !
  REAL(DP), ALLOCATABLE :: v1x_qe(:,:), v2x_qe(:,:), v3x_qe(:,:)
  REAL(DP), ALLOCATABLE :: v1c_qe(:,:), v2c_ud_qe(:)
  REAL(DP), ALLOCATABLE :: v2c_qe(:,:), v3c_qe(:,:)
  !
  !--------- LIBXC vars -----------------------
  CHARACTER(LEN=120) :: name1, name2
  !
  CHARACTER(LEN=30) :: qe_dft
  CHARACTER(LEN=30) :: lxc_dft
  !
  INTEGER :: pol_unpol
  !
  REAL(DP), ALLOCATABLE :: tau(:,:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: exg_lxc(:), ecg_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_lxc(:,:), vc_lxc(:,:)
  REAL(DP), ALLOCATABLE :: dmuxc_lxc(:,:,:)
  REAL(DP), ALLOCATABLE :: dvxcrr_lxc(:,:,:), dvxcsr_lxc(:,:,:), &
                           dvxcss_lxc(:,:,:)
  !
  REAL(DP), ALLOCATABLE :: v1x_lxc(:,:), v2x_lxc(:,:), v3x_lxc(:,:)
  REAL(DP), ALLOCATABLE :: v1c_lxc(:,:), v2c_lxc(:,:), v2c_ud_lxc(:)
  REAL(DP), ALLOCATABLE :: v2cm(:,:,:), v3c_lxc(:,:)
  !
  !----------diff vars --------------------------
  LOGICAL :: ex_is_out, ec_is_out, vx_is_out, vc_is_out, &
             something_out, dmuxc_is_out
  LOGICAL :: v1x_is_out, v2x_is_out, v1c_is_out, v2c_is_out, &
             v3x_is_out, v3c_is_out, &
             dvxcrr_is_out, dvxcsr_is_out, dvxcss_is_out, dvgga_is_out
  !
  REAL(DP) :: diff_thr_v_lda, diff_thr_e_lda, diff_thr_dmuxc
  REAL(DP) :: diff_thr_e_gga, diff_thr_v_gga, diff_thr_dvxc
  REAL(DP) :: diff_thr_e_mgga, diff_thr_v_mgga
  !
  REAL(DP), ALLOCATABLE :: grho2(:)
  REAL(DP) :: grho_ud
  !
  !
  REAL(DP) :: ex_aver(2), ec_aver(2), ex_max(2), ec_max(2), ex_min(2), ec_min(2)
  !
  REAL(DP) :: vx_aver(2,2), vx_max(2,2), vx_min(2,2), &
              vc_aver(2,2), vc_max(2,2), vc_min(2,2)
  REAL(DP) :: dmuxc_aver(2,3), dmuxc_max(2,3), dmuxc_min(2,3)
              
  REAL(DP) :: vaver(2), vmax(2), vmin(2)
  REAL(DP) :: v1x_aver(2,2), v1x_max(2,2), v1x_min(2,2), &
              v1c_aver(2,2), v1c_max(2,2), v1c_min(2,2)
  REAL(DP) :: v2x_aver(2,2), v2x_max(2,2), v2x_min(2,2), &
              v2c_aver(2,3), v2c_max(2,3), v2c_min(2,3), &
              v3x_aver(2,2), v3x_max(2,2), v3x_min(2,2), &
              v3c_aver(2,2), v3c_max(2,2), v3c_min(2,2)
  !
  REAL(DP) :: dvxcrr_aver(2,3), dvxcrr_max(2,3), dvxcrr_min(2,3), &
              dvxcsr_aver(2,3), dvxcsr_max(2,3), dvxcsr_min(2,3), &
              dvxcss_aver(2,3), dvxcss_max(2,3), dvxcss_min(2,3)
  !
  ! *******************************************************************************
  ! *-----------------------------------------------------------------------------*
  ! * To find the names of the Libxc functionals look at:                         *
  ! *                                                                             *
  ! *        https://tddft.org/programs/libxc/functionals/                        *
  ! *                                                                             *
  ! * For QE names see the comments in Modules/funct.f90                          *
  ! *******************************************************************************
  !
  !
  ! *******************************************************************************
  ! * Compatibility for GGA with lyp:                                             *
  ! *                                                                             *
  ! *  qe:     ec = ec_lyp*rho + ec_glyp  / vc = vc_lyp + v1c_glyp                *
  ! *  libxc:  ec = ec_glyp (icorr=131)   / vc = vc_glyp (131)                    *
  ! *         ... same for polarized case                                         *
  ! *******************************************************************************
  !
  !
  PRINT *, CHAR(10)//" --- XC TEST BETWEEN QE AND LIBXC ---"//CHAR(10)//" "
  !
  WRITE (*,'(/,1x,a)', ADVANCE='no') "# points? "
  READ(*,*) nnr
  !
  WRITE (*,'(/,1x,a)', ADVANCE='no') "Derivative of xc (y/n)? "
  READ(*,*) f_q
  DF_OK = .FALSE.
  IF ( TRIM(f_q) == 'y' ) DF_OK = .TRUE.
  IF ( TRIM(f_q) /= 'y' .AND. TRIM(f_q) /= 'n' ) THEN
     PRINT *, CHAR(10)//"Wrong answer"//CHAR(10)
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
       PRINT *, CHAR(10)//"Wrong answer"//CHAR(10)
       GO TO 10
    ENDIF
  ENDIF
  !
  !
  WRITE (*,'(/,1x,a)', ADVANCE='no') "QE-input_dft:     "
  READ(*,'(A)') qe_dft
  WRITE (*,'(/,1x,a)', ADVANCE='no') "LIBXC-input_dft:  "
  READ(*,'(A)') lxc_dft
  WRITE (*,'(/,1x,a)', ADVANCE='no') "Polarization switch (1 unpolarized,  & 
                                                         & 2 polarized):  "
  READ(*,*) ns
  IF ( ns/=1 .AND. ns/=2 ) THEN
     PRINT *, CHAR(10)//"ERROR: you can only choose 1 or 2"//CHAR(10)
     GO TO 10
  ENDIF
  !
  CALL set_dft_from_name( qe_dft )
  !
  LDA = .FALSE.
  GGA = .FALSE.
  MGGA= .FALSE.
  !
  iexch_qe = get_iexch()
  icorr_qe = get_icorr()
  IF (iexch_qe+icorr_qe/=0)  LDA = .TRUE.
  igcx_qe = get_igcx()
  igcc_qe = get_igcc()
  IF (igcx_qe+igcc_qe/=0)    GGA = .TRUE.
  imeta_qe  = get_meta()
  imetac_qe = get_metac()
  IF (imeta_qe+imetac_qe/=0) MGGA = .TRUE.
  !
  IF ( MGGA .AND. DF_OK ) THEN
    PRINT *, " "
    PRINT *, "ERROR: No derivative of V with MGGA."
    GO TO 10
  ENDIF
  !
  IF (ns == 2 .AND. icorr_qe/=0 .AND. icorr_qe/=1 .AND. icorr_qe/=2 .AND. &
                    icorr_qe/=4 .AND. icorr_qe/=8 .AND. icorr_qe/=3 .AND. &
                    icorr_qe/=7 .AND. icorr_qe/=13) THEN
     PRINT *, CHAR(10)//" ERROR: icorr_qe not available at these conditions"//CHAR(10)
     GO TO 10
  ENDIF
  !
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
  ! *  qe =>  rho_qe(n,1) -> up    |      libxc     =>      rho_lxc(2n+1) -> up   *
  ! *         rho_qe(n,2) -> down  |            (dim=2*nnr) rho_lxc(2n+2) -> down *
  ! *                              |                                              *
  ! *         grho(n,1)  -> uu     |            (dim=3*nnr) sigma(3n+1) -> uu     *
  ! *         grho(n,2)  -> dd     |                        sigma(3n+2) -> ud     *
  ! *         grho_ud(n) -> ud     |                        sigma(3n+3) -> dd     *
  ! *                                                                             *
  ! *******************************************************************************
  !
  !
  ! ------ ALLOCATIONS --------
  !
  IF ( LDA )  nthr = 1
  IF ( GGA )  nthr = 2
  IF ( MGGA ) nthr = 3
  !
  ! ... input
  !
  ALLOCATE( rho(nnr+nthr,ns) )
  ALLOCATE( rho_tz(nnr+nthr,ns) )
  !
  IF ( GGA .OR. MGGA ) ALLOCATE( grho(3,nnr+nthr,ns) )
  IF ( MGGA ) ALLOCATE( tau(nnr+nthr,ns) )
  !
  ! ... q-e output
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( LDA ) THEN
       IF ( .NOT. DF_OK ) THEN
         ALLOCATE( ex_qe(nnr+nthr), ec_qe(nnr+nthr) )
         ALLOCATE( vx_qe(nnr+nthr,ns), vc_qe(nnr+nthr,ns) )
       ELSE
         ALLOCATE( dmuxc_qe(nnr+nthr,ns,ns) )
       ENDIF 
     ENDIF
     !
     IF ( GGA ) THEN
       IF ( .NOT. DF_OK ) THEN
         ALLOCATE( exg_qe(nnr+nthr), ecg_qe(nnr+nthr) )
         ALLOCATE( v1x_qe(nnr+nthr,ns), v2x_qe(nnr+nthr,ns) )
         ALLOCATE( v1c_qe(nnr+nthr,ns) )
         ALLOCATE( v2c_qe(nnr+nthr,ns), v2c_ud_qe(nnr+nthr) )
       ELSE
         ALLOCATE( grh(nnr+nthr,3,ns) )
         ALLOCATE( dvxcrr_qe(nnr+nthr,ns,ns), dvxcsr_qe(nnr+nthr,ns,ns), &
                   dvxcss_qe(nnr+nthr,ns,ns) )
       ENDIF
     ENDIF
  ELSEIF ( MGGA ) THEN
     ALLOCATE( ex_qe(nnr+nthr), ec_qe(nnr+nthr) )
     ALLOCATE( v1x_qe(nnr+nthr,ns), v2x_qe(nnr+nthr,ns) )
     ALLOCATE( v1c_qe(nnr+nthr,ns) )
     ALLOCATE( v2c_qe(nnr+nthr,ns) )
     ALLOCATE( v3x_qe(nnr+nthr,ns), v3c_qe(nnr+nthr,ns) )
  ENDIF
  !
  ! ... libxc output
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( LDA ) THEN
       IF ( .NOT. DF_OK ) THEN
         ALLOCATE( ex_lxc(nnr+nthr), ec_lxc(nnr+nthr) )
         ALLOCATE( vx_lxc(nnr+nthr,ns), vc_lxc(nnr+nthr,ns) )
       ELSE
         ALLOCATE( dmuxc_lxc(nnr+nthr,ns,ns) )
       ENDIF 
     ENDIF
     !
     IF ( GGA ) THEN
       IF ( .NOT. DF_OK ) THEN
         ALLOCATE( exg_lxc(nnr+nthr), ecg_lxc(nnr+nthr) )
         ALLOCATE( v1x_lxc(nnr+nthr,ns), v2x_lxc(nnr+nthr,ns) )
         ALLOCATE( v1c_lxc(nnr+nthr,ns) )
         ALLOCATE( v2c_lxc(nnr+nthr,ns), v2c_ud_lxc(nnr+nthr) )
       ELSE
         ALLOCATE( dvxcrr_lxc(nnr+nthr,ns,ns), dvxcsr_lxc(nnr+nthr,ns,ns), &
                   dvxcss_lxc(nnr+nthr,ns,ns) )
       ENDIF
     ENDIF
  ELSEIF ( MGGA ) THEN
     ALLOCATE( ex_lxc(nnr+nthr), ec_lxc(nnr+nthr) )
     ALLOCATE( v1x_lxc(nnr+nthr,ns), v2x_lxc(nnr+nthr,ns) )
     ALLOCATE( v1c_lxc(nnr+nthr,ns) )
     ALLOCATE( v2c_lxc(nnr+nthr,ns) )
     ALLOCATE( v3x_lxc(nnr+nthr,ns), v3c_lxc(nnr+nthr,ns) )
  ENDIF
  !
  ! ... diff output
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( GGA ) THEN
       ALLOCATE( grho2(ns) )
     ENDIF
  ELSEIF ( MGGA ) THEN
     ALLOCATE( grho2(ns) )
  ENDIF
  !
  !----- BUILD INPUT -----  
  !  
  rho  = 0.0_DP  
  IF ( GGA .OR. MGGA ) grho = 0.0_DP  
  IF ( MGGA ) tau = 0.0_DP  
  !  
  DO ii = 1, nnr
     !  
     ! LDA: tot,diff  
     ! GGA: up,dw  
     rho(ii,1) = DBLE(ii)/DBLE(nnr+2)  
     IF (.NOT. POLARIZED) rho_tz(ii,1) = rho(ii,1)  
     !  
     IF ( GGA .OR. MGGA ) THEN  
       grho(1,ii,1) = ABS( 0.05_DP + 0.8_DP*SIN(DBLE(ii)) ) !0.1d0  
       grho(2,ii,1) = ABS( 0.05_DP + 0.7_DP*SIN(DBLE(ii)) ) !0.1d0  
       grho(3,ii,1) = ABS( 0.05_DP + 0.6_DP*SIN(DBLE(ii)) ) !0.1d0  
     ENDIF  
     !  
     IF ( MGGA ) tau(ii,1) = ABS( 0.05_DP + 0.8_DP*SIN(DBLE(ii)) )*0.5_DP  
     !  
     IF ( POLARIZED ) THEN  
        !  
        rho(ii,2) = (1.0_DP - rho(ii,1))*0.7_DP  
        rho_tz(ii,1) = rho(ii,1) + rho(ii,2)  
        rho_tz(ii,2) = rho(ii,1) - rho(ii,2)  
        !  
        IF ( GGA .OR. MGGA ) THEN  
           grho(1,ii,2) = ABS( (1.0_DP - grho(1,ii,1))*0.7_DP ) !0.1d0   
           grho(2,ii,2) = ABS( (1.0_DP - grho(2,ii,1))*0.6_DP ) !0.1d0   
           grho(3,ii,2) = ABS( (1.0_DP - grho(3,ii,1))*0.5_DP ) !0.1d0   
        ENDIF  
        !  
        IF ( MGGA ) tau(ii,2) = ABS( 0.05_DP + 0.8_DP*SIN(DBLE(ii)) )*0.2_DP  
        !  
     ENDIF  
     !  
  ENDDO  
  !  
  !  
  !--- THRESHOLD POINTS ---  
  rho(nnr+1,1) = thresh_lda/3.0_DP  
  IF (.NOT. POLARIZED) rho_tz(nnr+1,1) = rho(nnr+1,1)
  IF ( POLARIZED ) THEN
     rho(nnr+1,2) = rho(nnr,1)  
     rho_tz(nnr+1,1) = rho(nnr+1,1) + rho(nnr+1,2)
     rho_tz(nnr+1,2) = rho(nnr+1,1) - rho(nnr+1,2)
  ENDIF
  !  
  IF ( GGA .OR. MGGA ) THEN  
    grho(:,nnr+1,1) = grho(:,nnr,1)  
    !  
    rho(nnr+2,:) = rho(nnr,:)
    IF (.NOT. POLARIZED) rho_tz(nnr+2,1) = rho(nnr+2,1)
    !
    grho(:,nnr+2,1) = thresh_gga/10
    IF ( POLARIZED ) THEN
      rho_tz(nnr+2,1) = rho(nnr+2,1) + rho(nnr+2,2)
      rho_tz(nnr+2,2) = rho(nnr+2,1) - rho(nnr+2,2)
      grho(:,nnr+2,2) = grho(:,nnr,2)
    ENDIF
    !
    IF ( MGGA ) THEN  
      tau(nnr+1,:) = tau(nnr,:)
      tau(nnr+2,:) = tau(nnr,:)
      !
      rho(nnr+3,:) = rho(nnr,:)
      grho(:,nnr+3,:) = grho(:,nnr,:)
      !
      tau(nnr+3,1) = thresh_mgga/10
      IF ( POLARIZED ) tau(nnr+3,2) = tau(nnr,2)
    ENDIF
    !
  ENDIF
  !
  !
  diff_thr_e_lda = 1.0E-6_DP
  diff_thr_v_lda = 1.0E-6_DP
  diff_thr_dmuxc = 1.0E-6_DP
  !
  diff_thr_e_gga = 1.0E-12_DP
  diff_thr_v_gga = 1.0E-12_DP
  diff_thr_dvxc  = 1.0E-16_DP
  !
  diff_thr_e_mgga = 1.0E-12_DP
  diff_thr_v_mgga = 1.0E-12_DP
  !
  !
  !-------- Calculation of energies and potential ------------------
  !
  IF ( LDA ) THEN
     !
     IF (.NOT. DF_OK ) THEN 
       CALL xc( nnr+nthr, ns, ns, rho_tz, ex_qe, ec_qe, vx_qe, vc_qe )
     ELSE
       CALL dmxc( nnr+nthr, ns, rho, dmuxc_qe )
     ENDIF
     !
  ENDIF
  !
  IF ( GGA ) THEN
    !
    IF ( .NOT. DF_OK ) THEN
      !
      IF ( .NOT. LDA ) THEN
        ex_qe = 0.d0  ;  ec_qe = 0.d0
        vx_qe = 0.d0  ;  vc_qe = 0.d0
      ENDIF
      !
      CALL xc_gcx( nnr+nthr, ns, rho, grho, exg_qe, ecg_qe, v1x_qe, v2x_qe, v1c_qe, v2c_qe, v2c_ud_qe )
      !
      ex_qe = ex_qe*rho_tz(:,1) + exg_qe
      ec_qe = ec_qe*rho_tz(:,1) + ecg_qe
      !
      v1x_qe = v1x_qe + vx_qe
      v1c_qe = v1c_qe + vc_qe
      !
    ELSE
      !
      DO ii = 1, nnr+nthr
        grh(ii,1:3,1) = grho(1:3,ii,1)
        IF (ns==2) grh(ii,1:3,2) = grho(1:3,ii,2)
      ENDDO
      !
      CALL dgcxc( nnr+nthr, ns, rho, grh, dvxcrr_qe, dvxcsr_qe, dvxcss_qe )
      !
      dvxcrr_qe = dvxcrr_qe + dmuxc_qe
      !
    ENDIF
     !
  ENDIF
  !
  !
  IF ( MGGA ) THEN
     ALLOCATE( v2cm(np,nnr+nthr,ns) )
     CALL xc_metagcx( nnr+nthr, ns, np, rho, grho, tau, ex_qe, &
                      ec_qe, v1x_qe, v2x_qe, v3x_qe,        &
                      v1c_qe, v2cm, v3c_qe )
     !
     v2c_qe = v2cm(1,:,:)
     DEALLOCATE( v2cm )
  ENDIF
  !
  !-----------------------------------------------LXC
  !
  CALL reset_dft()
  !
  CALL set_dft_from_name( lxc_dft )
  !
  LDA = .FALSE.
  GGA = .FALSE.
  MGGA= .FALSE.
  iexch_lxc = get_iexch()
  icorr_lxc = get_icorr()
  IF (iexch_lxc+icorr_lxc/=0)  LDA = .TRUE.
  igcx_lxc = get_igcx()
  igcc_lxc = get_igcc()
  IF (igcx_lxc+igcc_lxc/=0)    GGA = .TRUE.
  imeta_lxc  = get_meta()
  imetac_lxc = get_metac()
  IF (imeta_lxc+imetac_lxc/=0) MGGA = .TRUE.
  !
  !
  IF ( LDA ) THEN
    !
    IF ( .NOT. DF_OK ) THEN
      !
      CALL xc( nnr+nthr, ns, ns, rho_tz, ex_lxc, ec_lxc, vx_lxc, vc_lxc )
      !
    ELSE
      !
      CALL dmxc( nnr+nthr, ns, rho, dmuxc_lxc )
      !
    ENDIF
    !
  ENDIF
  !
  IF ( GGA ) THEN
    !
    IF ( .NOT. DF_OK ) THEN
      !
      IF ( .NOT. LDA ) THEN
        ex_lxc = 0.d0  ;  ec_lxc = 0.d0
        vx_lxc = 0.d0  ;  vc_lxc = 0.d0
      ENDIF
      !
      CALL xc_gcx( nnr+nthr, ns, rho, grho, exg_lxc, ecg_lxc, v1x_lxc, v2x_lxc, &
                   v1c_lxc, v2c_lxc, v2c_ud_lxc )
      !
      ex_lxc = ex_lxc*rho_tz(:,1) + exg_lxc
      ec_lxc = ec_lxc*rho_tz(:,1) + ecg_lxc
      !
      v1x_lxc = v1x_lxc + vx_lxc
      v1c_lxc = v1c_lxc + vc_lxc
      !
    ELSE
      !
      CALL dgcxc( nnr+nthr, ns, rho, grh, dvxcrr_lxc, dvxcsr_lxc, dvxcss_lxc )
      !
    ENDIF
    !
  ENDIF
  !
  !
  IF ( MGGA ) THEN
    ALLOCATE( v2cm(np,nnr+nthr,ns) )
    CALL xc_metagcx( nnr+nthr, ns, np, rho, grho, tau, ex_lxc, &
                     ec_lxc, v1x_lxc, v2x_lxc, v3x_lxc,   &
                     v1c_lxc, v2cm, v3c_lxc )
    v2c_lxc = v2cm(1,:,:)
    DEALLOCATE( v2cm )
  ENDIF
  !
  ! ... testing the name-reading part
  !
  PRINT *, " "
  PRINT *, "=================================== "//CHAR(10)//" "
  PRINT *, "QE functional IDs:"//CHAR(10)//" "
  !
  ! ... q-e
  !
  PRINT *, "- LDA IDs: ",  iexch_qe, icorr_qe
  !
  IF ( GGA )  PRINT *, "- GGA IDs: ",  igcx_qe,  igcc_qe
  !
  IF ( MGGA ) PRINT *, "- MGGA IDs: ", imeta_qe, imetac_qe
  !
  ! ... libxc
  !
  PRINT *, " "
  PRINT *, "LIBXC functional IDs and info:"//CHAR(10)//" "
  IF ( LDA ) THEN
    !
    name1 = xc_f03_functional_get_name( iexch_lxc )
    name2 = xc_f03_functional_get_name( icorr_lxc )
    !
    PRINT *, "- LDA IDs: ", iexch_lxc, icorr_lxc
    PRINT *, "- LDA exch: " , TRIM(name1)
    PRINT *, "- LDA corr: " , TRIM(name2)
    PRINT *, " "
  ENDIF
  !
  IF ( GGA ) THEN
    !
    name1 = xc_f03_functional_get_name( igcx_lxc )
    name2 = xc_f03_functional_get_name( igcc_lxc )
    !
    PRINT *, "- GGA IDs: ", igcx_lxc, igcc_lxc
    PRINT *, "- GGA exch: " , TRIM(name1)
    PRINT *, "- GGA corr: " , TRIM(name2)
    PRINT *, " "
  ENDIF
  !
  IF ( MGGA ) THEN
    !
    name1 = xc_f03_functional_get_name( imeta_lxc )
    name2 = xc_f03_functional_get_name( imetac_lxc )
    !
    PRINT *, "- MGGA IDs: ", imeta_lxc, imetac_lxc
    PRINT *, "- MGGA exch: " , TRIM(name1)
    PRINT *, "- MGGA corr: " , TRIM(name2)
    PRINT *, " "
  ENDIF
  !
  !
  IF ( .NOT. DF_OK ) THEN
    !
    CALL evxc_stats( 'Ex', diff_thr_e_lda, ex_qe, ex_lxc )
    CALL evxc_stats( 'Ec', diff_thr_e_lda, ec_qe, ec_lxc )
    !
  ENDIF
  !
  IF ( LDA .AND. .NOT. GGA ) THEN
     !
     IF ( .NOT. DF_OK ) THEN
       !
       ! V exchange
       !
       CALL evxc_stats( 'Vx', diff_thr_v_lda, vx_qe, vx_lxc )
       !
       ! V correlation
       !
       CALL evxc_stats( 'Vc', diff_thr_v_lda, vc_qe, vc_lxc )
       !
     ELSE
       !
       ! Derivative of Vxc
       !
       CALL derivxc_stats( 'dmuxc', diff_thr_v_lda, dmuxc_qe, dmuxc_lxc )
       !
     ENDIF
     !
     vx_is_out = .FALSE.
     vc_is_out = .FALSE.
     something_out = .FALSE.
     dmuxc_is_out = .FALSE.
     !
     DO ii = 1, nnr
        !
        IF ( .NOT. DF_OK ) THEN
          ex_is_out = is_it_out( diff_thr_e_lda, 1, ex_qe(ii:ii), ex_lxc(ii:ii) )
          ec_is_out = is_it_out( diff_thr_e_lda, 1, ec_qe(ii:ii), ec_lxc(ii:ii) )
          !
          IF (.NOT. energy_only) THEN
            vx_is_out = is_it_out( diff_thr_v_lda,ns, vx_qe(ii,:), vx_lxc(ii,:) )
            vc_is_out = is_it_out( diff_thr_v_lda,ns, vc_qe(ii,:), vc_lxc(ii,:) )
          ENDIF
          something_out = ANY( (/ex_is_out, ec_is_out, vx_is_out, vc_is_out/) )
        ELSE
          dmuxc_is_out = is_dit_out( diff_thr_dmuxc, dmuxc_qe(ii,:,:), dmuxc_lxc(ii,:,:) )
        ENDIF
        !
        PRINT *, " "
        WRITE(*,909) ii, nnr
        !  
        IF ( .NOT. POLARIZED ) THEN
           WRITE (*, 401 ) rho(ii,1)
        ELSE
           WRITE (*, 402 ) rho(ii,1), rho(ii,2)
        ENDIF
        PRINT *, " "
        !
        IF ( something_out .OR. dmuxc_is_out ) THEN
          !
          IF (.NOT. DF_OK) THEN
            !
            IF ( ex_is_out ) CALL print_diff( 'Ex', ex_qe(ii:ii), ex_lxc(ii:ii) )
            IF ( ec_is_out ) CALL print_diff( 'Ec', ec_qe(ii:ii), ec_lxc(ii:ii) )
            !
            IF (.NOT. ENERGY_ONLY) THEN
               IF ( vx_is_out ) CALL print_diff( 'Vx', vx_qe(ii,:), vx_lxc(ii,:) )
               IF ( vc_is_out ) CALL print_diff( 'Vc', vc_qe(ii,:), vc_lxc(ii,:) )
            ENDIF
            !
          ELSE  
            !
            CALL print_diff2( 'dmuxc', dmuxc_qe(ii,:,:), dmuxc_lxc(ii,:,:) )
            !
          ENDIF !df_ok
          !
        ELSE !something out
          !
          PRINT *, "... match "
          !
        ENDIF
        !
     ENDDO
     !
     !
     ! ============ THRESHOLD TEST ==================
     PRINT *, " "
     PRINT *, "--- INPUT THRESHOLD CHECK ---"
     PRINT *, " "
     !
     DO ithr = 1, nthr
       PRINT *, " "
       WRITE(*,910) ithr, nthr
       !  
       IF ( .NOT. POLARIZED ) THEN
         WRITE (*, 401 ) rho(nnr+ithr,1)
       ELSE
         WRITE (*, 402 ) rho(nnr+ithr,1), rho(nnr+ithr,2)         
       ENDIF
       PRINT *, " "
       !
       IF (.NOT. DF_OK) THEN
         CALL print_diff( 'Ex', ex_qe(nnr+ithr:nnr+ithr), ex_lxc(nnr+ithr:nnr+ithr) )
         CALL print_diff( 'Ec', ec_qe(nnr+ithr:nnr+ithr), ec_lxc(nnr+ithr:nnr+ithr) )
         IF (.NOT. ENERGY_ONLY) THEN
           CALL print_diff( 'Vx', vx_qe(nnr+ithr,:), vx_lxc(nnr+ithr,:) )
           CALL print_diff( 'Vc', vc_qe(nnr+ithr,:), vc_lxc(nnr+ithr,:) )
         ENDIF
       ELSE
         CALL print_diff2( 'dmuxc', dmuxc_qe(nnr+ithr,:,:), dmuxc_lxc(nnr+ithr,:,:) )
       ENDIF !df_ok
     ENDDO
     !
     ! ===============================================
     !
   ELSEIF ( GGA ) THEN
      !
      IF ( .NOT. DF_OK ) THEN
        !
        ! V1 exchange
        !
        CALL evxc_stats( 'V1x', diff_thr_v_gga, v1x_qe, v1x_lxc )
        !
        ! V2 exchange
        !
        CALL evxc_stats( 'V2x', diff_thr_v_gga, v2x_qe, v2x_lxc )
        !
        ! V1 correlation
        !
        CALL evxc_stats( 'V1c', diff_thr_v_gga, v1c_qe, v1c_lxc )
        !
        ! V2 correlation
        !
        CALL evxc_stats( 'V2c', diff_thr_v_gga, v2c_qe, v2c_lxc )
        !
        IF ( POLARIZED ) THEN
          CALL diff_average( diff_thr_v_gga, v2c_ud_qe, v2c_ud_lxc, vaver )
          CALL diff_max( diff_thr_v_gga, v2c_ud_qe, v2c_ud_lxc, vmax )
          CALL diff_min( diff_thr_v_gga, v2c_ud_qe, v2c_ud_lxc, vmin )
          !
          CALL print_stat( 'cross', vaver, vmax, vmin )
        ENDIF
        !
      ELSE
        !
        CALL derivxc_stats( 'dvxcrr', diff_thr_v_gga, dvxcrr_qe, dvxcrr_lxc )
        !
        CALL derivxc_stats( 'dvxcsr', diff_thr_v_gga, dvxcsr_qe, dvxcsr_lxc )
        !
        CALL derivxc_stats( 'dvxcss', diff_thr_v_gga, dvxcss_qe, dvxcss_lxc )
        ! 
      ENDIF
      !
      !
      v1x_is_out = .FALSE.
      v2x_is_out = .FALSE.
      v1c_is_out = .FALSE.
      v2c_is_out = .FALSE.
      something_out = .FALSE.
      dvgga_is_out = .FALSE.
      !
      DO ii = 1, nnr
         !
         IF ( .NOT. DF_OK ) THEN
           ex_is_out = is_it_out( diff_thr_e_gga, 1, ex_qe(ii:ii), ex_lxc(ii:ii) )
           ec_is_out = is_it_out( diff_thr_e_gga, 1, ec_qe(ii:ii), ec_lxc(ii:ii) )
           !
           IF (.NOT. energy_only) THEN
             v1x_is_out = is_it_out( diff_thr_v_gga,ns, v1x_qe(ii,:), v1x_lxc(ii,:) )
             v2x_is_out = is_it_out( diff_thr_v_gga,ns, v2x_qe(ii,:), v2x_lxc(ii,:) )
             v1c_is_out = is_it_out( diff_thr_v_gga,ns, v1c_qe(ii,:), v1c_lxc(ii,:) )
             v2c_is_out = is_it_out( diff_thr_v_gga,np, v2c_qe(ii,:), v2c_lxc(ii,:), v2c_ud_qe(ii), v2c_ud_lxc(ii) )
           ENDIF
           something_out = ANY( (/ex_is_out, ec_is_out, v1x_is_out, v2x_is_out, &
                                  v1c_is_out, v2c_is_out/) )
         ELSE           
           dvxcrr_is_out = is_dit_out( diff_thr_dvxc, dvxcrr_qe(ii,:,:), dvxcrr_lxc(ii,:,:) )
           dvxcsr_is_out = is_dit_out( diff_thr_dvxc, dvxcsr_qe(ii,:,:), dvxcsr_lxc(ii,:,:) )
           dvxcss_is_out = is_dit_out( diff_thr_dvxc, dvxcss_qe(ii,:,:), dvxcss_lxc(ii,:,:) )
           !
           dvgga_is_out  = ANY((/dvxcrr_is_out, dvxcsr_is_out, dvxcss_is_out/))
         ENDIF
         !
         WRITE(*,*) " "
         WRITE(*,909) ii, nnr
         IF (.NOT. POLARIZED ) THEN
            WRITE (*,401) rho(ii,1)
            grho2(1) = grho(1,ii,1)**2 + grho(2,ii,1)**2 + grho(3,ii,1)**2
            WRITE (*,501) grho2(1)
         ELSE
            WRITE (*,402) rho(ii,1), rho(ii,2)
            grho2(1) = grho(1,ii,1)**2 + grho(2,ii,1)**2 + grho(3,ii,1)**2
            grho2(2) = grho(1,ii,2)**2 + grho(2,ii,2)**2 + grho(3,ii,2)**2
            grho_ud  = grho(1,ii,1) * grho(1,ii,2) + &
                       grho(2,ii,1) * grho(2,ii,2) + &
                       grho(3,ii,1) * grho(3,ii,2)
            WRITE (*,503) grho2(1), grho_ud, grho2(2)
         ENDIF
         PRINT *, " "
         !
         IF ( something_out .OR. dvgga_is_out ) THEN
           !
           IF (.NOT. DF_OK) THEN
             !
             IF ( ex_is_out ) CALL print_diff( 'Ex', ex_qe(ii:ii), ex_lxc(ii:ii) )
             IF ( ec_is_out ) CALL print_diff( 'Ec', ec_qe(ii:ii), ec_lxc(ii:ii) )
             !
             IF (.NOT. ENERGY_ONLY) THEN
               !
               IF ( v1x_is_out ) CALL print_diff( 'V1x', v1x_qe(ii,:), v1x_lxc(ii,:) )
               IF ( v2x_is_out ) CALL print_diff( 'V2x', v2x_qe(ii,:), v2x_lxc(ii,:) )
               IF ( v1c_is_out ) CALL print_diff( 'V1c', v1c_qe(ii,:), v1c_lxc(ii,:) )
               IF ( v2c_is_out ) CALL print_diff( 'V2c', v2c_qe(ii,:), v2c_lxc(ii,:),  &
                                                         v2c_ud_qe(ii), v2c_ud_lxc(ii) )
               !
             ENDIF
             !
           ELSE
             !
             PRINT *, " "
             !
             IF ( dvxcrr_is_out ) CALL print_diff2( 'dvxcrr', dvxcrr_qe(ii,:,:), dvxcrr_lxc(ii,:,:) )
             IF ( dvxcsr_is_out ) CALL print_diff2( 'dvxcsr', dvxcsr_qe(ii,:,:), dvxcsr_lxc(ii,:,:) )
             IF ( dvxcss_is_out ) CALL print_diff2( 'dvxcss', dvxcss_qe(ii,:,:), dvxcss_lxc(ii,:,:) )
             !
           ENDIF
           !
         ELSE !something_out
           !
           PRINT *, "... match "
           !
         ENDIF
         !
      ENDDO
      !
      ! ============ THRESHOLD TEST ==================
      PRINT *, " "
      PRINT *, "--- INPUT THRESHOLD CHECK ---"
      PRINT *, " "
      !
      DO ithr = 1, nthr
        !
        IF (.NOT. DF_OK) THEN
          WRITE(*,*) " "
          WRITE(*,910) ithr, nthr
          IF (.NOT. POLARIZED ) THEN
            WRITE (*,401) rho(nnr+ithr,1)
            grho2(1) = grho(nnr+ithr,1,1)**2 + grho(nnr+ithr,2,1)**2 + grho(nnr+ithr,3,1)**2
            WRITE (*,501) grho2(1)
          ELSE
            WRITE (*,402) rho(nnr+ithr,1), rho(nnr+ithr,2)
            grho2(1) = grho(1,nnr+ithr,1)**2 + grho(2,nnr+ithr,1)**2 + grho(3,nnr+ithr,1)**2
            grho2(2) = grho(1,nnr+ithr,2)**2 + grho(2,nnr+ithr,2)**2 + grho(3,nnr+ithr,2)**2
            grho_ud  = grho(1,nnr+ithr,1) * grho(1,nnr+ithr,2) + &
                       grho(2,nnr+ithr,1) * grho(2,nnr+ithr,2) + &
                       grho(3,nnr+ithr,1) * grho(3,nnr+ithr,2)
            WRITE (*,503) grho2(1), grho_ud, grho2(2)
          ENDIF
          PRINT *, " "
          !
          CALL print_diff( 'Ex', ex_qe(nnr+ithr:nnr+ithr), ex_lxc(nnr+ithr:nnr+ithr) )
          CALL print_diff( 'Ec', ec_qe(nnr+ithr:nnr+ithr), ec_lxc(nnr+ithr:nnr+ithr) )
          IF (.NOT. ENERGY_ONLY) THEN
            CALL print_diff( 'V1x', v1x_qe(nnr+ithr,:), v1x_lxc(nnr+ithr,:) )
            CALL print_diff( 'V2x', v2x_qe(nnr+ithr,:), v2x_lxc(nnr+ithr,:) )
            CALL print_diff( 'V1c', v1c_qe(nnr+ithr,:), v1c_lxc(nnr+ithr,:) )
            CALL print_diff( 'V2c', v2c_qe(nnr+ithr,:), v2c_lxc(nnr+ithr,:),  &
                                    v2c_ud_qe(nnr+ithr), v2c_ud_lxc(nnr+ithr) )
          ENDIF
        ELSE  
          CALL print_diff2( 'dvxcrr', dvxcrr_qe(nnr+ithr,:,:), dvxcrr_lxc(nnr+ithr,:,:) )
          CALL print_diff2( 'dvxcsr', dvxcsr_qe(nnr+ithr,:,:), dvxcsr_lxc(nnr+ithr,:,:) )
          CALL print_diff2( 'dvxcss', dvxcss_qe(nnr+ithr,:,:), dvxcss_lxc(nnr+ithr,:,:) )
        ENDIF
        !
      ENDDO !df_ok
      ! ===============================================
      !
      !
   ELSEIF ( MGGA ) THEN
     !
     ! V1 exchange
     !
     CALL evxc_stats( 'V1x', diff_thr_v_mgga, v1x_qe, v1x_lxc )
     !   
     ! V2 exchange   
     !   
     CALL evxc_stats( 'V2x', diff_thr_v_mgga, v2x_qe, v2x_lxc )
     !
     ! V3 exchange
     !
     CALL evxc_stats( 'V3x', diff_thr_v_mgga, v3x_qe, v3x_lxc )
     !
     ! V1 correlation
     !
     CALL evxc_stats( 'V1c', diff_thr_v_mgga, v1c_qe, v1c_lxc )
     !
     ! V2 correlation
     !
     CALL evxc_stats( 'V2c', diff_thr_v_mgga, v2c_qe, v2c_lxc )
     !
     ! V3 correlation
     !
     CALL evxc_stats( 'V3c', diff_thr_v_mgga, v3c_qe, v3c_lxc )
     !
     v1x_is_out = .FALSE.
     v2x_is_out = .FALSE.
     v3x_is_out = .FALSE.
     v1c_is_out = .FALSE.
     v2c_is_out = .FALSE.
     v3c_is_out = .FALSE.
     something_out = .FALSE.
     !
     DO ii = 1, nnr
        !
        ex_is_out = is_it_out( diff_thr_e_mgga, 1, ex_qe(ii:ii), ex_lxc(ii:ii) )
        ec_is_out = is_it_out( diff_thr_e_mgga, 1, ec_qe(ii:ii), ec_lxc(ii:ii) )
        !
        IF (.NOT. energy_only) THEN
          v1x_is_out = is_it_out( diff_thr_v_mgga,ns, v1x_qe(ii,:), v1x_lxc(ii,:) )
          v2x_is_out = is_it_out( diff_thr_v_mgga,ns, v2x_qe(ii,:), v2x_lxc(ii,:) )
          v3x_is_out = is_it_out( diff_thr_v_mgga,ns, v3x_qe(ii,:), v3x_lxc(ii,:) )
          v1c_is_out = is_it_out( diff_thr_v_mgga,ns, v1c_qe(ii,:), v1c_lxc(ii,:) )
          v2c_is_out = is_it_out( diff_thr_v_mgga,ns, v2c_qe(ii,:), v2c_lxc(ii,:) )
          v3c_is_out = is_it_out( diff_thr_v_mgga,ns, v3c_qe(ii,:), v3c_lxc(ii,:) )
        ENDIF
        !
        something_out = ANY( (/ex_is_out, ec_is_out, v1x_is_out, v2x_is_out, &
                               v3x_is_out, v1c_is_out, v2c_is_out, &
                               v3c_is_out/) )
        !
        WRITE(*,*) " "
        WRITE(*,909) ii, nnr
        IF (.NOT. POLARIZED ) THEN
           WRITE (*,401) rho(ii,1)
           grho2(1) = grho(1,ii,1)**2 + grho(2,ii,1)**2 + grho(3,ii,1)**2
           WRITE (*,501) grho2(1)
           WRITE (*,601) tau(ii,1)
        ELSE
           WRITE (*,402) rho(ii,1), rho(ii,2)
           grho2(1) = grho(1,ii,1)**2 + grho(2,ii,1)**2 + grho(3,ii,1)**2
           grho2(2) = grho(1,ii,2)**2 + grho(2,ii,2)**2 + grho(3,ii,2)**2
           grho_ud  = grho(1,ii,1) * grho(1,ii,2) + &
                      grho(2,ii,1) * grho(2,ii,2) + &
                      grho(3,ii,1) * grho(3,ii,2)
           WRITE (*,503) grho2(1), grho_ud, grho2(2)
           WRITE (*,602) tau(ii,1), tau(ii,2)
        ENDIF
        !
        IF ( something_out ) THEN
           !
           IF ( ex_is_out ) CALL print_diff( 'Ex', ex_qe(ii:ii), ex_lxc(ii:ii) )
           IF ( ec_is_out ) CALL print_diff( 'Ec', ec_qe(ii:ii), ec_lxc(ii:ii) )
           !
           IF (.NOT. ENERGY_ONLY) THEN
             IF ( v1x_is_out ) CALL print_diff( 'V1x', v1x_qe(ii,:), v1x_lxc(ii,:) )
             IF ( v2x_is_out ) CALL print_diff( 'V2x', v2x_qe(ii,:), v2x_lxc(ii,:) )
             IF ( v3x_is_out ) CALL print_diff( 'V3x', v3x_qe(ii,:), v3x_lxc(ii,:) )
             IF ( v1c_is_out ) CALL print_diff( 'V1c', v1c_qe(ii,:), v1c_lxc(ii,:) )
             IF ( v2c_is_out ) CALL print_diff( 'V2c', v2c_qe(ii,:), v2c_lxc(ii,:) )
             IF ( v3c_is_out ) CALL print_diff( 'V3c', v3c_qe(ii,:), v3c_lxc(ii,:) )
           ENDIF !not energy only
           !  
        ELSE ! something out
           !
           PRINT *, "... match "
           !
        ENDIF
        !
     ENDDO
     !
     ! ============ THRESHOLD TEST ================== 
     PRINT *, " " 
     PRINT *, "--- INPUT THRESHOLD CHECK ---" 
     PRINT *, " " 
     ! 
     DO ithr = 1, nthr 
       !
       WRITE(*,*) " "
       WRITE(*,910) ithr, nthr   
       IF (.NOT. POLARIZED ) THEN   
         WRITE (*,401) rho(nnr+ithr,1)
         grho2(1) = grho(1,nnr+ithr,1)**2 + grho(2,nnr+ithr,1)**2 + grho(3,nnr+ithr,1)**2   
         WRITE (*,501) grho2(1) 
         WRITE (*,601) tau(nnr+ithr,1)
       ELSE   
         WRITE (*,402) rho(nnr+ithr,1), rho(nnr+ithr,2)   
         grho2(1) = grho(1,nnr+ithr,1)**2 + grho(2,nnr+ithr,1)**2 + grho(3,nnr+ithr,1)**2   
         grho2(2) = grho(1,nnr+ithr,2)**2 + grho(2,nnr+ithr,2)**2 + grho(3,nnr+ithr,2)**2   
         grho_ud  = grho(1,nnr+ithr,1) * grho(1,nnr+ithr,2) + &   
                    grho(2,nnr+ithr,1) * grho(2,nnr+ithr,2) + &   
                    grho(3,nnr+ithr,1) * grho(3,nnr+ithr,2)   
         WRITE (*,503) grho2(1), grho_ud, grho2(2)
         WRITE (*,602) tau(nnr+ithr,1), tau(nnr+ithr,2)
       ENDIF   
       PRINT *, " "   
       !   
       CALL print_diff( 'Ex', ex_qe(nnr+ithr:nnr+ithr), ex_lxc(nnr+ithr:nnr+ithr) )   
       CALL print_diff( 'Ec', ec_qe(nnr+ithr:nnr+ithr), ec_lxc(nnr+ithr:nnr+ithr) )   
       IF (.NOT. ENERGY_ONLY) THEN   
         CALL print_diff( 'V1x', v1x_qe(nnr+ithr,:), v1x_lxc(nnr+ithr,:) )   
         CALL print_diff( 'V2x', v2x_qe(nnr+ithr,:), v2x_lxc(nnr+ithr,:) )   
         CALL print_diff( 'V3x', v3x_qe(nnr+ithr,:), v3x_lxc(nnr+ithr,:) )   
         CALL print_diff( 'V1c', v1c_qe(nnr+ithr,:), v1c_lxc(nnr+ithr,:) )   
         CALL print_diff( 'V2c', v2c_qe(nnr+ithr,:), v2c_lxc(nnr+ithr,:) )   
         CALL print_diff( 'V3c', v3c_qe(nnr+ithr,:), v3c_lxc(nnr+ithr,:) )   
       ENDIF   
       ! 
     ENDDO 
     ! =============================================== 
     !
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
  502 FORMAT('grho2(uu,dd): ',F17.14,4x,F17.14)
  503 FORMAT('grho2(uu,ud,dd): ',F17.14,4x,F17.14,4x,F17.14)
  !
  601 FORMAT('tau: ',F17.14)
  602 FORMAT('tau(up,down): ',F17.14,4x,F17.14)
  !
  909 FORMAT('grid-point: ',I4,' of',I4)
  910 FORMAT('threshold-point: ',I4,' of',I4)
  !
  !
  !
  ! ------ DEALLOCATIONS --------
  !
  DEALLOCATE( rho, rho_tz )
  !
  IF ( GGA .OR. MGGA ) DEALLOCATE( grho )
  IF ( MGGA ) DEALLOCATE( tau )
  !
  ! ... q-e output
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( LDA ) THEN
       IF ( .NOT. DF_OK ) THEN
         DEALLOCATE( ex_qe, ec_qe )
         DEALLOCATE( vx_qe, vc_qe )
       ELSE
         DEALLOCATE( dmuxc_qe )
       ENDIF 
     ENDIF
     !
     IF ( GGA ) THEN
       IF ( .NOT. DF_OK ) THEN
         DEALLOCATE( exg_qe, ecg_qe )
         DEALLOCATE( v1x_qe, v2x_qe )
         DEALLOCATE( v1c_qe )
         DEALLOCATE( v2c_qe, v2c_ud_qe )
       ELSE
         DEALLOCATE( grh )
         DEALLOCATE( dvxcrr_qe, dvxcsr_qe, dvxcss_qe )
       ENDIF
     ENDIF
  ELSEIF ( MGGA ) THEN
     DEALLOCATE( ex_qe, ec_qe )
     DEALLOCATE( v1x_qe, v2x_qe )
     DEALLOCATE( v1c_qe )
     DEALLOCATE( v2c_qe )
     DEALLOCATE( v3x_qe, v3c_qe )
  ENDIF
  !
  ! ... libxc output
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( LDA ) THEN
       IF ( .NOT. DF_OK ) THEN
         DEALLOCATE( ex_lxc, ec_lxc )
         DEALLOCATE( vx_lxc, vc_lxc )
       ELSE
         DEALLOCATE( dmuxc_lxc )
       ENDIF 
     ENDIF
     !
     IF ( GGA ) THEN
       IF ( .NOT. DF_OK ) THEN
         DEALLOCATE( exg_lxc, ecg_lxc )
         DEALLOCATE( v1x_lxc, v2x_lxc )
         DEALLOCATE( v1c_lxc )
         DEALLOCATE( v2c_lxc, v2c_ud_lxc )
       ELSE
         DEALLOCATE( dvxcrr_lxc, dvxcsr_lxc, dvxcss_lxc )
       ENDIF
     ENDIF
  ELSEIF ( MGGA ) THEN
     DEALLOCATE( ex_lxc, ec_lxc )
     DEALLOCATE( v1x_lxc, v2x_lxc )
     DEALLOCATE( v1c_lxc )
     DEALLOCATE( v2c_lxc )
     DEALLOCATE( v3x_lxc, v3c_lxc )
  ENDIF
  !
  PRINT *, " "
  !
#else
  !
  PRINT *, "ERROR: library libxc not linked."
  !
#endif
10 STOP
  !
#if defined(__LIBXC)
  !
 CONTAINS
!
!
!------------------------------------------------------------------------
SUBROUTINE diff_average( thr, x_qe, x_lxc, aver_abs_perc )
  !----------------------------------------------------------------------
  !! Calculates average difference (both absolute and percentage) between
  !! qe and libxc quantities.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: thr
  REAL(DP), INTENT(IN)  :: x_qe(nnr), x_lxc(nnr)
  REAL(DP), INTENT(OUT) :: aver_abs_perc(2)
  !! 1: absolute difference;  2: percentage difference
  !
  INTEGER  :: i, nnr_int
  REAL(DP) :: calc_perc_diff
  REAL(DP) :: abs_diff, perc_diff
  !
  nnr_int = 0
  aver_abs_perc = 0.d0
  !
  DO i = 1, nnr
    !
    abs_diff  = ABS(x_qe(i) - x_lxc(i))
    perc_diff = calc_perc_diff( thr, x_qe(i), x_lxc(i) )
    !
    aver_abs_perc(1) = aver_abs_perc(1) + abs_diff
    !
    IF ( perc_diff < 0.d0 ) CYCLE
    !
    nnr_int = nnr_int+1
    aver_abs_perc(2) = aver_abs_perc(2) + perc_diff
    !
  ENDDO
  !
  aver_abs_perc(1) = aver_abs_perc(1) / DBLE(nnr)
  aver_abs_perc(2) = aver_abs_perc(2) / DBLE(nnr_int)
  !
  RETURN
  !
END SUBROUTINE diff_average
!
!
!-------------------------------------------------------------------------
SUBROUTINE diff_max( thr, x_qe, x_lxc, max_abs_perc )
  !-----------------------------------------------------------------------
  !! Finds the max difference (both absolute and percentage) between qe
  !! and libxc quantities.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: thr
  REAL(DP), INTENT(IN)  :: x_qe(nnr), x_lxc(nnr)
  REAL(DP), INTENT(OUT) :: max_abs_perc(2)
  !
  INTEGER :: i
  REAL(DP) :: calc_perc_diff
  REAL(DP) :: abs_diff, perc_diff
  REAL(DP) :: abs_diff_prev, perc_diff_prev
  !
  abs_diff_prev  = 0.d0
  perc_diff_prev = 0.d0
  !
  DO i = 1, nnr
    !
    abs_diff  = ABS(x_qe(i) - x_lxc(i))
    perc_diff = calc_perc_diff( thr, x_qe(i), x_lxc(i) )
    !
    IF ( abs_diff > abs_diff_prev ) THEN
      max_abs_perc(1) = abs_diff
      abs_diff_prev = abs_diff
    ENDIF
    !
    IF ( perc_diff > perc_diff_prev ) THEN
      max_abs_perc(2) = perc_diff
      perc_diff_prev = perc_diff
    ENDIF
    !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE diff_max
!
!
!--------------------------------------------------------------------------
SUBROUTINE diff_min( thr, x_qe, x_lxc, min_abs_perc )
  !------------------------------------------------------------------------
  !! Finds the max difference (both absoulte and percentage) between qe
  !! and libxc quantities.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: thr
  REAL(DP), INTENT(IN)  :: x_qe(nnr), x_lxc(nnr)
  REAL(DP), INTENT(OUT) :: min_abs_perc(2)
  !
  INTEGER  :: i
  REAL(DP) :: calc_perc_diff
  REAL(DP) :: abs_diff, perc_diff
  REAL(DP) :: perc_diff_prev, abs_diff_prev
  !
  abs_diff_prev  = 1000.d0
  perc_diff_prev = 1000.d0
  !
  DO i = 1, nnr
    !
    abs_diff  = ABS(x_qe(i) - x_lxc(i))
    perc_diff = calc_perc_diff( thr, x_qe(i), x_lxc(i) )
    !
    IF ( abs_diff < abs_diff_prev .and. abs_diff>0.d0 ) THEN
      min_abs_perc(1) = abs_diff
      abs_diff_prev = abs_diff
    ENDIF
    !
    IF ( perc_diff < 0.d0 ) CYCLE
    !
    IF ( perc_diff < perc_diff_prev ) THEN
      min_abs_perc(2) = perc_diff
      perc_diff_prev = perc_diff
    ENDIF
    !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE diff_min
!
!--------------------------------------------------------------------
SUBROUTINE print_stat( what, vaver, vmax, vmin )
  !------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: what
  REAL(DP), DIMENSION(2), INTENT(IN) :: vaver, vmax, vmin
  !
  PRINT *, " "   
  PRINT *, " ", TRIM(what)   
  PRINT *, "AVR abs diff: ", vaver(1), "   AVR % diff: ", vaver(2)   
  PRINT *, "MAX abs diff: ", vmax(1),  "   MAX % diff: ", vmax(2)   
  PRINT *, "MIN abs diff: ", vmin(1),  "   MIN % diff: ", vmin(2)
  !
END SUBROUTINE
!
!------------------------------------------------------------------
SUBROUTINE print_diff( what, x_qe, x_lxc, x_ud_qe, x_ud_lxc )
  !-----------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: x_qe(ns), x_lxc(ns)
  REAL(DP), INTENT(IN), OPTIONAL :: x_ud_qe, x_ud_lxc
  !
  PRINT *," "
  PRINT *, what
  !
  IF ( .NOT. POLARIZED .OR. what(1:1)=='E'  ) THEN
    WRITE (*,101) x_qe(1)
    WRITE (*,201) x_lxc(1)
    PRINT *, " --- "
    WRITE (*,301) ABS(x_qe(1)-x_lxc(1))
  ELSEIF ( POLARIZED ) THEN
    IF ( .NOT. PRESENT(x_ud_qe) ) THEN
      WRITE (*,102) x_qe(1), x_qe(2)
      WRITE (*,202) x_lxc(1), x_lxc(2)
      PRINT *, " --- "
      WRITE (*,302) ABS(x_qe(1)-x_lxc(1)), &
                    ABS(x_qe(2)-x_lxc(2))
    ELSE
      WRITE (*,103) x_qe(1), x_ud_qe, x_qe(2)
      WRITE (*,203) x_lxc(1), x_ud_lxc, x_lxc(2)
      PRINT *, " --- "
      WRITE (*,303) ABS(x_qe(1)-x_lxc(1)), &
                    ABS(x_ud_qe-x_ud_lxc), &
                    ABS(x_qe(2)-x_lxc(2))
    ENDIF
  ENDIF
  !
  101 FORMAT('qe: ',3x,F17.14)
  102 FORMAT('qe: ',3x,F17.14,4x,F17.14)
  103 FORMAT('qe: ',3x,F17.14,4x,F17.14,4x,F17.14)
  201 FORMAT('libxc: ',F17.14)
  202 FORMAT('libxc: ',F17.14,4x,F17.14)
  203 FORMAT('libxc: ',F17.14,4x,F17.14,4x,F17.14)
  301 FORMAT('diffs: ',F17.14)
  302 FORMAT('diffs: ',F17.14,4x,F17.14)
  303 FORMAT('diffs: ',F17.14,4x,F17.14,4x,F17.14)
  !
END SUBROUTINE print_diff
!
!-----------------------------------------------------------------------
SUBROUTINE print_diff2( what, dxc_qe, dxc_lxc )
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: dxc_qe(ns,ns), dxc_lxc(ns,ns)
  !
  PRINT *, " "
  !
  IF ( POLARIZED ) PRINT *, what
  !
  IF ( .NOT. POLARIZED ) THEN   
    WRITE (*,101) dxc_qe(1,1)
    WRITE (*,201) dxc_lxc(1,1)
    PRINT *, " --- "
    WRITE (*,301) dxc_qe(1,1)-dxc_lxc(1,1)
  ELSE
    WRITE (*,103) dxc_qe(1,1), dxc_qe(2,1), dxc_qe(2,2) !, dxc_qe(1,2)
    WRITE (*,203) dxc_lxc(1,1), dxc_lxc(2,1), dxc_lxc(2,2)
    PRINT *, " --- "
    WRITE (*,303) dxc_qe(1,1)-dxc_lxc(1,1), &
                  dxc_qe(2,1)-dxc_lxc(2,1), &
                  dxc_qe(2,2)-dxc_lxc(2,2)
  ENDIF
  !
  101 FORMAT('qe: ',3x,F17.14)
  103 FORMAT('qe: ',3x,F17.14,4x,F17.14,4x,F17.14)
  201 FORMAT('libxc: ',F17.14)
  203 FORMAT('libxc: ',F17.14,4x,F17.14,4x,F17.14)
  301 FORMAT('diffs: ',F17.14)
  303 FORMAT('diffs: ',F17.14,4x,F17.14,4x,F17.14)
  !
END SUBROUTINE print_diff2
!
!---------------------------------------------------------------------
SUBROUTINE calc_stats( thr, x_qe, x_lxc, x_aver, x_max, x_min )
  !-------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: thr
  REAL(DP), INTENT(IN)  :: x_qe(nnr), x_lxc(nnr)
  REAL(DP), INTENT(OUT) :: x_aver(2), x_max(2), x_min(2)
  !
  CALL diff_average( thr, x_qe, x_lxc, x_aver )
  CALL diff_max( thr, x_qe, x_lxc, x_max )
  CALL diff_min( thr, x_qe, x_lxc, x_min )
  !
  RETURN
  !
END SUBROUTINE calc_stats
!
!
SUBROUTINE evxc_stats( what, thr, xc_qe, xc_lxc )
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: what
  REAL(DP), INTENT(IN) :: thr, xc_qe(nnr+nthr,ns), xc_lxc(nnr+nthr,ns)
  REAL(DP) :: xc_aver(2,ns), xc_max(2,ns), xc_min(2,ns)
  !
  PRINT *," "
  !
  IF ( POLARIZED .AND. what(1:1)/='E' ) PRINT *, what
  !
  CALL diff_average( thr, xc_qe(1:nnr,1), xc_lxc(1:nnr,1), xc_aver(1:nnr,1) )
  CALL diff_max( thr, xc_qe(1:nnr,1), xc_lxc(1:nnr,1), xc_max(1:nnr,1) )
  CALL diff_min( thr, xc_qe(1:nnr,1), xc_lxc(1:nnr,1), xc_min(1:nnr,1) )
  !
  IF ( .NOT. POLARIZED .OR. what(1:1)=='E' ) THEN
    CALL print_stat( what, xc_aver(1:nnr,1), xc_max(1:nnr,1), xc_min(1:nnr,1) )
  ELSE
    CALL diff_average( thr, xc_qe(1:nnr,2), xc_lxc(1:nnr,2), xc_aver(1:nnr,2) )
    CALL diff_max( thr, xc_qe(1:nnr,2), xc_lxc(1:nnr,2), xc_max(1:nnr,2) )
    CALL diff_min( thr, xc_qe(1:nnr,2), xc_lxc(1:nnr,2), xc_min(1:nnr,2) )
    !
    CALL print_stat( 'up', xc_aver(1:nnr,1), xc_max(1:nnr,1), xc_min(1:nnr,1) )
    CALL print_stat( 'down', xc_aver(1:nnr,2), xc_max(1:nnr,2), xc_min(1:nnr,2) )
  ENDIF
  !
END SUBROUTINE evxc_stats

!
!
!---------------------------------------------------------------------
SUBROUTINE derivxc_stats( what, thr, dxc_qe, dxc_lxc )
  !-------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: thr
  REAL(DP), INTENT(IN)  :: dxc_qe(nnr+nthr,ns,ns), dxc_lxc(nnr+nthr,ns,ns)
  REAL(DP) :: dxc_aver(2,np), dxc_max(2,np), dxc_min(2,np)
  !
  PRINT *," "
  PRINT *, what
  !
  CALL diff_average( thr, dxc_qe(1:nnr,1,1), dxc_lxc(1:nnr,1,1), dxc_aver(1:nnr,1) )
  CALL diff_max( thr, dxc_qe(1:nnr,1,1), dxc_lxc(1:nnr,1,1), dxc_max(1:nnr,1) )
  CALL diff_min( thr, dxc_qe(1:nnr,1,1), dxc_lxc(1:nnr,1,1), dxc_min(1:nnr,1) )
  !
  IF ( .NOT. POLARIZED ) THEN
    CALL print_stat( what, dxc_aver(1:nnr,1), dxc_max(1:nnr,1), dxc_min(1:nnr,1) )
  ELSE
    CALL diff_average( thr, dxc_qe(1:nnr,1,2), dxc_lxc(1:nnr,1,2), dxc_aver(1:nnr,2) )
    CALL diff_max( thr, dxc_qe(1:nnr,1,2), dxc_lxc(1:nnr,1,2), dxc_max(1:nnr,2) )
    CALL diff_min( thr, dxc_qe(1:nnr,1,2), dxc_lxc(1:nnr,1,2), dxc_min(1:nnr,2) )
    !
    CALL diff_average( thr, dxc_qe(1:nnr,2,2), dxc_lxc(1:nnr,2,2), dxc_aver(1:nnr,3) )
    CALL diff_max( thr, dxc_qe(1:nnr,2,2), dxc_lxc(1:nnr,2,2), dxc_max(1:nnr,3) )
    CALL diff_min( thr, dxc_qe(1:nnr,2,2), dxc_lxc(1:nnr,2,2), dxc_min(1:nnr,3) )
    !
    CALL print_stat( 'up-up', dxc_aver(1:nnr,1), dxc_max(1:nnr,1), dxc_min(1:nnr,1) )
    CALL print_stat( 'up-down', dxc_aver(1:nnr,2), dxc_max(1:nnr,2), dxc_min(1:nnr,2) )
    CALL print_stat( 'down-down', dxc_aver(1:nnr,3), dxc_max(1:nnr,3), dxc_min(1:nnr,3) )
  ENDIF
  !
END SUBROUTINE derivxc_stats
!
!------------------------------------------------------------------------
FUNCTION is_it_out( diff_thr, dm, x_qe, x_lxc, x_ud_qe, x_ud_lxc )
  !----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: dm
  REAL(DP), INTENT(IN) :: diff_thr, x_qe(dm), x_lxc(dm)
  REAL(DP), INTENT(IN), OPTIONAL :: x_ud_qe, x_ud_lxc
  LOGICAL :: is_it_out, is_it_out_ud
  !
  is_it_out = ANY(ABS(x_qe(1:dm)-x_lxc(1:dm)) > diff_thr)
  !
  IF (PRESENT(x_ud_qe)) THEN
    is_it_out_ud =  ABS(x_ud_qe-x_ud_lxc) > diff_thr
    is_it_out = ANY( (/ is_it_out, is_it_out_ud /) )
  ENDIF
  !
END FUNCTION
!
!------------------------------------------------------------------------
FUNCTION is_dit_out( diff_thr, dx_qe, dx_lxc )
  !----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: diff_thr, dx_qe(ns,ns), dx_lxc(ns,ns)
  REAL(DP) :: dxc_diff(np)
  LOGICAL :: is_dit_out
  !
  dxc_diff(1) = dx_qe(1,1)-dx_lxc(1,1)
  IF ( POLARIZED ) THEN
    dxc_diff(2) = dx_qe(2,1)-dx_lxc(2,1)
    dxc_diff(3) = dx_qe(2,2)-dx_lxc(2,2)
  ENDIF
  !
  is_dit_out = ANY(dxc_diff(:) > diff_thr)
  !
END FUNCTION
!
#endif
!
END PROGRAM xctest_qe_libxc
!
!
#if defined(__LIBXC)
!------------------------------------------------------------------------
FUNCTION calc_perc_diff( thr, x_qe, x_lxc )
  !----------------------------------------------------------------------
  !! Calculates difference between qe and libxc quantities in percentage.
  !
  IMPLICIT NONE
  !
  REAL(8), INTENT(IN) :: thr
  REAL(8), INTENT(IN) :: x_qe, x_lxc
  REAL(8) :: calc_perc_diff
  !
  REAL(8) :: perc_diff
  !
  perc_diff = -1.d0
  !
  IF ( ABS(x_qe)<10.d0*thr .AND. ABS(x_qe-x_lxc)<10.d0*thr ) RETURN
  IF ( ABS(x_qe)==0.d0 .AND. ABS(x_qe-x_lxc)>thr ) calc_perc_diff = 100.d0
  IF ( ABS(x_qe)>thr ) calc_perc_diff = ABS( (x_qe-x_lxc)/x_qe )*100.d0
  !
  RETURN
  !
END FUNCTION calc_perc_diff
#endif
