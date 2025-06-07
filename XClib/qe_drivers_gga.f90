!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!========================================================================
!                GRADIENT CORRECTION DRIVERS for E and V
!========================================================================
!
!------------------------------------------------------------------------
MODULE qe_drivers_gga
  !----------------------------------------------------------------------
  !! Contains the GGA drivers that calculate the XC energy and potential.
  !
  USE kind_l,               ONLY: DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: gcxc, gcx_spin, gcc_spin, gcc_spin_more
  !
  !
CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE gcxc( length, rho_in, grho_in, sx_out, sc_out, v1x_out, &
                                        v2x_out, v1c_out, v2c_out, err_out )
  !---------------------------------------------------------------------
  !! Gradient corrections for exchange and correlation - Hartree a.u. 
  !! See comments at the beginning of module for implemented cases
  !
  USE dft_setting_params,   ONLY: igcx, igcc, rho_threshold_gga,     &
                                  grho_threshold_gga, exx_started,   &
                                  exx_fraction, screening_parameter, &
                                  gau_parameter
  USE exch_gga
  USE corr_gga
  USE beef_interface, ONLY: beefx, beeflocalcorr
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! Length of the input/output arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! Charge density
  REAL(DP), INTENT(IN),  DIMENSION(length) :: grho_in
  !! \(\text{grho}=|\nabla rho|^2\)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sx_out
  !! Exchange energy: \(s_x = \int e_x(\text{rho},\text{grho}) dr\)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sc_out
  !! Correlation energy: \(s_c = \int e_c(\text{rho},\text{grho}) dr\)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v1x_out
  !! Exchange potential: \(D\ E_x\ /\ D\ \text{rho} \)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v2x_out
  !! Exchange potential: \(D\ E_x\ /\ D(D\text{rho}/D r_\alpha)\ /
  !! \ |\nabla\text{rho}| \)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v1c_out
  !! Correlation potential (density term)
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c_out
  !! Correlation potential (gradient term)
  INTEGER, INTENT(OUT) :: err_out
  !! error index
  !
  ! ... local variables
  !
  INTEGER :: ir, in_err, iflag ! Added iflag for AH series
  REAL(DP) :: rho, grho
  REAL(DP) :: sx, v1x, v2x
  REAL(DP) :: sx_, v1x_, v2x_
  REAL(DP) :: sxsr, v1xsr, v2xsr
  REAL(DP) :: sc, v1c, v2c
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
  err_out = 0
  !
#if defined(_OPENACC)
! ... workaround for atomic-related bug in hpc_sdk 21.5 and older
#if defined(__PGI) && (__PGIC__ < 21 || (__PGIC__ == 21 && __PGIC_MINOR__ < 7))
!$acc data present( rho_in, grho_in, sx_out, sc_out, v1x_out, v2x_out, v1c_out, v2c_out )
#else
!$acc data present( rho_in, grho_in, sx_out, sc_out, v1x_out, v2x_out, v1c_out, v2c_out ) copy( err_out )
#endif
!$acc parallel loop  
#else
!$omp parallel if(ntids==1) default(none) &
!$omp private( rho, grho, sx, sx_, sxsr, v1x, v1x_, v1xsr, &
!$omp          v2x, v2x_, v2xsr, sc, v1c, v2c, iflag, in_err ) &
!$omp shared( rho_in, grho_in, length, igcx, exx_started, &
!$omp         grho_threshold_gga, rho_threshold_gga, gau_parameter, &
!$omp         screening_parameter, exx_fraction, igcc, v1x_out, v2x_out, &
!$omp         v1c_out, v2c_out, sx_out, sc_out, err_out )
!$omp do
#endif
  DO ir = 1, length  
     !
     in_err = 0
     grho = grho_in(ir)
     !
     IF ( rho_in(ir) <= rho_threshold_gga .OR. grho <= grho_threshold_gga ) THEN
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
          CALL pbexsr( rho, grho, sxsr, v1xsr, v2xsr, screening_parameter, in_err )
          sx  = sx  - exx_fraction * sxsr
          v1x = v1x - exx_fraction * v1xsr
          v2x = v2x - exx_fraction * v2xsr
        ENDIF
        !
     CASE( 34, 35 ) ! ' AH series for GGA cross checks
        !
        iflag = 0
        IF ( igcx== 34 ) THEN ! PBE-AH cross check
           CALL pbex( rho, grho, 1, sx, v1x, v2x )
           iflag = 1 ! AHPB for PBE cross check
        ELSEIF ( igcx== 35 ) THEN ! PBESOL-AH cross check
           CALL pbex( rho, grho, 3, sx, v1x, v2x )
           iflag = 2 ! AHPS for PBEsol-based cross check
        ENDIF
        !
        IF ( iflag == 0) in_err = 4 ! Sorting GGA-AHs failed
        !
        IF (exx_started) THEN
          CALL axsr( iflag, rho, grho, sxsr, v1xsr, v2xsr, screening_parameter, in_err )
          sx  = sx  - exx_fraction * sxsr
          v1x = v1x - exx_fraction * v1xsr
          v2x = v2x - exx_fraction * v2xsr
        ENDIF
        !
     CASE( 32, 33, 47 ) ! 'AH series for vdW-DFs, JPCM 34, 025902 (2022)
        !
        iflag = 0
        IF ( igcx == 32) THEN ! vdW-DF-ahcx
           CALL cx13( rho, grho, sx, v1x, v2x )
           iflag = 3 ! for cx13 - analytical sr hole
        ELSEIF ( igcx == 33) THEN ! vdW-DF2-ah
           CALL rPW86( rho, grho, sx, v1x, v2x )
           iflag = 4 ! for rPW86 - analytical sr hole
        ELSEIF ( igcx == 47) THEN ! vdW-DF2-ahbr
           CALL b86b( rho, grho, 3, sx, v1x, v2x ) 
           iflag = 6 ! for test-reserve - analytical sr hole
        ENDIF
        !
        IF ( iflag == 0) in_err = 5  ! Sorting vdW-DF-AHs failed
        !
        IF (exx_started) THEN
          CALL axsr( iflag, rho, grho, sxsr, v1xsr, v2xsr, screening_parameter, in_err )
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
        !
     CASE( 43 ) ! 'BEEX'
        !
        CALL beefx( rho, grho, sx, v1x, v2x, 0 )
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
        !
     CASE( 14 ) !'BEEC'
        ! last parameter 0 means: do not add lda contributions
        ! espresso will do that itself
        CALL beeflocalcorr( rho, grho, sc, v1c, v2c, 0 )
        !
     CASE DEFAULT
        !
        sc = 0.0_DP
        v1c = 0.0_DP
        v2c = 0.0_DP
        !
     END SELECT
     !
     IF (in_err/=0) THEN
#if defined(_OPENACC)
!$acc atomic write
#else
!$omp atomic write
#endif
       err_out = in_err
     ENDIF
     !
     sx_out(ir)  = sx    ;  sc_out(ir)  = sc
     v1x_out(ir) = v1x   ;  v1c_out(ir) = v1c
     v2x_out(ir) = v2x   ;  v2c_out(ir) = v2c
     !
  ENDDO
#if defined(_OPENACC)
!$acc end data
#else
!$omp end do
!$omp end parallel
#endif
  !
  RETURN
  !
END SUBROUTINE gcxc
!
!
!===============> SPIN <===============!
!
!-------------------------------------------------------------------------
SUBROUTINE gcx_spin( length, rho_in, grho2_in, sx_tot, v1x_out, v2x_out, err_out )
  !-----------------------------------------------------------------------
  !! Gradient corrections for exchange - Hartree a.u.
  !
  USE dft_setting_params,   ONLY: igcx, igcc, exx_started,   &
                                  exx_fraction, screening_parameter, &
                                  gau_parameter
  USE exch_gga
  USE beef_interface, ONLY: beefx
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Length of the input/output arrays
  REAL(DP), INTENT(IN),  DIMENSION(length,2) :: rho_in
  !! Up and down charge density
  REAL(DP), INTENT(IN),  DIMENSION(length,2) :: grho2_in
  !! Up and down gradient of the charge
  REAL(DP), INTENT(OUT), DIMENSION(length)   :: sx_tot
  !! Energy exchange GGA
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1x_out
  !! Exchange potential (density part)
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v2x_out
  !! Exchange potantial (gradient part)
  INTEGER, INTENT(OUT) :: err_out
  !! error index
  !
  ! ... local variables
  !
  INTEGER :: ir, iflag, in_err
  REAL(DP) :: rho_up, rho_dw, grho2_up, grho2_dw
  REAL(DP) :: v1x_up, v1x_dw, v2x_up, v2x_dw
  REAL(DP) :: sx_up, sx_dw, rnull_up, rnull_dw
  REAL(DP) :: sxsr_up, sxsr_dw
  REAL(DP) :: v1xsr_up, v1xsr_dw, v2xsr_up, v2xsr_dw
  !
  REAL(DP), PARAMETER :: small=1.D-10
  REAL(DP), PARAMETER :: rho_trash=0.5_DP, grho2_trash=0.2_DP
  ! temporary values assigned to rho and grho when they
  ! are too small in order to avoid numerical problems.
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
  err_out = 0
  !
#if defined(_OPENACC)
#if defined(__PGI) && (__PGIC__ < 21 || (__PGIC__ == 21 && __PGIC_MINOR__ < 7))
!$acc data present( rho_in, grho2_in, sx_tot, v1x_out, v2x_out )
#else
!$acc data present( rho_in, grho2_in, sx_tot, v1x_out, v2x_out ) copy( err_out )
#endif
!$acc parallel loop
#else
!$omp parallel if(ntids==1) default(none) &
!$omp private( rho_up, rho_dw, grho2_up, grho2_dw, rnull_up, rnull_dw, &
!$omp          sx_up, sx_dw, sxsr_up, sxsr_dw, v1xsr_up, v1xsr_dw, &
!$omp          v1x_up, v1x_dw, v2x_up, v2x_dw, v2xsr_up, v2xsr_dw, &
!$omp          iflag, in_err ) &
!$omp  shared( rho_in, length, grho2_in, sx_tot, v1x_out, v2x_out,  &
!$omp          igcx, exx_started, exx_fraction, screening_parameter,&
!$omp          gau_parameter, err_out )
!$omp do
#endif
  DO ir = 1, length  
     !
     in_err = 0
     !
     rho_up = rho_in(ir,1)
     rho_dw = rho_in(ir,2)
     grho2_up = grho2_in(ir,1)
     grho2_dw = grho2_in(ir,2)
     rnull_up = 1.0_DP
     rnull_dw = 1.0_DP
     !
     IF ( rho_up+rho_dw <= small ) THEN
        sx_tot(ir) = 0.0_DP
        v1x_out(ir,1) = 0.0_DP
        v2x_out(ir,1) = 0.0_DP
        v1x_out(ir,2) = 0.0_DP
        v2x_out(ir,2) = 0.0_DP
        CYCLE
     ELSE
        IF ( rho_up<=small .OR. SQRT(ABS(grho2_up))<=small ) THEN
          rho_up = rho_trash
          grho2_up = grho2_trash
          rnull_up = 0.0_DP
        ENDIF
        IF ( rho_dw<=small .OR. SQRT(ABS(grho2_dw))<=small ) THEN
          rho_dw = rho_trash
          grho2_dw = grho2_trash
          rnull_dw = 0.0_DP
        ENDIF
     ENDIF
     !
     ! ... exchange
     !
     SELECT CASE( igcx )
     CASE( 1 )
        !
        CALL becke88_spin( rho_up, rho_dw, grho2_up, grho2_dw, sx_up, sx_dw, &
                           v1x_up, v1x_dw, v2x_up, v2x_dw )
        !
        sx_tot(ir) = sx_up*rnull_up + sx_dw*rnull_dw
        !
     CASE( 2 )
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL ggax( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL ggax( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
     CASE( 3, 4, 8, 10, 12, 20, 23, 24, 25, 34, 35, 44, 45 )
        ! igcx=3:  PBE,  igcx=4:  revised PBE, igcx=8:  PBE0, igcx=10: PBEsol
        ! igcx=12: HSE,  igcx=20: gau-pbe,     igcx=23: obk8, igcx=24: ob86,
        ! igcx=25: ev93, igcx=34: PBE-AH, igcx=35: PBESOL-AH,
        ! igcx=44: RPBE,        igcx=45: W31X
        !
        iflag = 1
        IF ( igcx== 4 ) iflag = 2
        IF ( igcx==10 ) iflag = 3
        IF ( igcx==23 ) iflag = 5
        IF ( igcx==24 ) iflag = 6
        IF ( igcx==25 ) iflag = 7
        IF ( igcx==44 ) iflag = 8
        IF ( igcx==45 ) iflag = 9
        IF ( igcx==34 ) iflag = 1
        IF ( igcx==35 ) iflag = 3
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL pbex( rho_up, grho2_up, iflag, sx_up, v1x_up, v2x_up )
        CALL pbex( rho_dw, grho2_dw, iflag, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
        IF ( igcx == 8 .AND. exx_started ) THEN
           !
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x_up = (1.0_DP - exx_fraction) * v1x_up
           v1x_dw = (1.0_DP - exx_fraction) * v1x_dw
           v2x_up = (1.0_DP - exx_fraction) * v2x_up
           v2x_dw = (1.0_DP - exx_fraction) * v2x_dw
           !
        ELSEIF ( igcx == 12 .AND. exx_started ) THEN
           !
           CALL pbexsr( rho_up, grho2_up, sxsr_up, v1xsr_up, &
                                          v2xsr_up, screening_parameter, in_err )
           CALL pbexsr( rho_dw, grho2_dw, sxsr_dw, v1xsr_dw, &
                                          v2xsr_dw, screening_parameter, in_err )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction*0.5_DP*( sxsr_up*rnull_up + &
                                                           sxsr_dw*rnull_dw )
           v1x_up = v1x_up - exx_fraction * v1xsr_up
           v1x_dw = v1x_dw - exx_fraction * v1xsr_dw
           v2x_up = v2x_up - exx_fraction * v2xsr_up * 2.0_DP
           v2x_dw = v2x_dw - exx_fraction * v2xsr_dw * 2.0_DP
           !
        ELSEIF ( igcx == 34 .AND. exx_started ) THEN
           !
           CALL axsr( 1, rho_up, grho2_up, sxsr_up, v1xsr_up, &
                                          v2xsr_up, screening_parameter, in_err )
           CALL axsr( 1, rho_dw, grho2_dw, sxsr_dw, v1xsr_dw, &
                                          v2xsr_dw, screening_parameter, in_err )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction*0.5_DP * ( sxsr_up*rnull_up + &
                                                             sxsr_dw*rnull_dw )
           v1x_up = v1x_up - exx_fraction * v1xsr_up
           v1x_dw = v1x_dw - exx_fraction * v1xsr_dw
           v2x_up = v2x_up - exx_fraction * v2xsr_up * 2.0_DP
           v2x_dw = v2x_dw - exx_fraction * v2xsr_dw * 2.0_DP
           !
        ELSEIF ( igcx == 35 .AND. exx_started ) THEN
           !
           CALL axsr( 2, rho_up, grho2_up, sxsr_up, v1xsr_up, &
                                          v2xsr_up, screening_parameter, in_err )
           CALL axsr( 2, rho_dw, grho2_dw, sxsr_dw, v1xsr_dw, &
                                          v2xsr_dw, screening_parameter, in_err )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction*0.5_DP * ( sxsr_up*rnull_up + &
                                                             sxsr_dw*rnull_dw )
           v1x_up = v1x_up - exx_fraction * v1xsr_up
           v1x_dw = v1x_dw - exx_fraction * v1xsr_dw
           v2x_up = v2x_up - exx_fraction * v2xsr_up * 2.0_DP
           v2x_dw = v2x_dw - exx_fraction * v2xsr_dw * 2.0_DP
           !
        ELSEIF ( igcx == 20 .AND. exx_started ) THEN
           ! gau-pbe
           !CALL pbexgau_lsd( rho, grho2, sxsr, v1xsr, v2xsr, gau_parameter )
           CALL pbexgau( rho_up,grho2_up, sxsr_up, v1xsr_up,v2xsr_up, gau_parameter )
           CALL pbexgau( rho_dw,grho2_dw, sxsr_dw, v1xsr_dw,v2xsr_dw, gau_parameter )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction*0.5_DP * ( sxsr_up*rnull_up + &
                                                             sxsr_dw*rnull_dw )
           v1x_up = v1x_up - exx_fraction * v1xsr_up
           v1x_dw = v1x_dw - exx_fraction * v1xsr_dw
           v2x_up = v2x_up - exx_fraction * v2xsr_up * 2.0_DP
           v2x_dw = v2x_dw - exx_fraction * v2xsr_dw * 2.0_DP
           !
        ENDIF
        !
     CASE( 9 )                    ! B3LYP
        !
        CALL becke88_spin( rho_up, rho_dw, grho2_up, grho2_dw, sx_up, sx_dw, &
                           v1x_up, v1x_dw, v2x_up, v2x_dw )
        !
        sx_tot(ir) = sx_up*rnull_up + sx_dw*rnull_dw
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = 0.72_DP * sx_tot(ir)
           v1x_up = 0.72_DP * v1x_up ; v1x_dw = 0.72_DP * v1x_dw
           v2x_up = 0.72_DP * v2x_up ; v2x_dw = 0.72_DP * v2x_dw
        ENDIF
        !
     CASE( 11 )                   ! 'Wu-Cohen'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL wcx( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL wcx( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
     CASE( 13 )                   ! 'revised PW86 for vdw-df2'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL rPW86( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL rPW86( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
     CASE( 16 )                   ! 'c09x for vdw-df-c09.'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL c09x( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL c09x( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
     CASE( 21 )                   ! 'PW86'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL pw86( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL pw86( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
     CASE( 22 )                   ! 'B86B'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL becke86b( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL becke86b( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
      CASE( 26, 46 )                  ! 'B86R for rev-vdW-DF2'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        IF ( igcx==26 ) iflag = 3 ! B86R for rev-vdW-DF2
        IF ( igcx==46 ) iflag = 4 ! W32X for vdW-DF3-opt2
        CALL b86b( rho_up, grho2_up, iflag, sx_up, v1x_up, v2x_up )
        CALL b86b( rho_dw, grho2_dw, iflag, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
     CASE( 27 )                   ! 'cx13 for vdw-df-cx'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL cx13( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL cx13( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
     CASE( 28 )                   ! X3LYP
        !
        CALL becke88_spin( rho_up, rho_dw, grho2_up, grho2_dw, sx_up, sx_dw, &
                           v1x_up, v1x_dw, v2x_up, v2x_dw )
        !
        rho_up = 2.0_DP * rho_up
        rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up
        grho2_dw = 4.0_DP * grho2_dw
        !
        CALL pbex( rho_up, grho2_up, 1, sxsr_up, v1xsr_up, v2xsr_up )
        CALL pbex( rho_dw, grho2_dw, 1, sxsr_dw, v1xsr_dw, v2xsr_dw )
        !
        sx_tot(ir) = 0.5_DP*( sxsr_up*rnull_up + sxsr_dw*rnull_dw )*0.235_DP + &
                            (   sx_up*rnull_up +   sx_dw*rnull_dw )*0.765_DP
        v1x_up = v1xsr_up * 0.235_DP + v1x_up * 0.765_DP
        v1x_dw = v1xsr_dw * 0.235_DP + v1x_dw * 0.765_DP
        v2x_up = v2xsr_up * 0.235_DP * 2.0_DP + v2x_up * 0.765_DP
        v2x_dw = v2xsr_dw * 0.235_DP * 2.0_DP + v2x_dw * 0.765_DP
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = 0.709_DP * sx_tot(ir)
           v1x_up = 0.709_DP * v1x_up
           v1x_dw = 0.709_DP * v1x_dw
           v2x_up = 0.709_DP * v2x_up
           v2x_dw = 0.709_DP * v2x_dw
        ENDIF
        !
     CASE( 29, 31 )               ! 'cx0 for vdw-df-cx0' or `cx0p for vdW-DF-cx0p'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL cx13( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL cx13( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x_up = (1.0_DP - exx_fraction) * v1x_up
           v1x_dw = (1.0_DP - exx_fraction) * v1x_dw
           v2x_up = (1.0_DP - exx_fraction) * v2x_up
           v2x_dw = (1.0_DP - exx_fraction) * v2x_dw
        ENDIF
        !
     CASE( 30 )                   ! 'R860' = 'rPW86-0' for vdw-df2-0'
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL rPW86( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL rPW86( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x_up = (1.0_DP - exx_fraction) * v1x_up
           v1x_dw = (1.0_DP - exx_fraction) * v1x_dw
           v2x_up = (1.0_DP - exx_fraction) * v2x_up
           v2x_dw = (1.0_DP - exx_fraction) * v2x_dw
        ENDIF
        !
     CASE( 38 )                  ! 'br0 for vdw-df2-BR0' etc
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL b86b( rho_up, grho2_up, 3, sx_up, v1x_up, v2x_up )
        CALL b86b( rho_dw, grho2_dw, 3, sx_dw, v1x_dw, v2x_dw )     
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x_up = (1.0_DP - exx_fraction) * v1x_up
           v1x_dw = (1.0_DP - exx_fraction) * v1x_dw
           v2x_up = (1.0_DP - exx_fraction) * v2x_up
           v2x_dw = (1.0_DP - exx_fraction) * v2x_dw
        ENDIF  
        !
     CASE( 32, 33, 47 ) ! ! 'AH series for vdW-DFs, JPCM 34, 025902 (2022)
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        ! igcx=32:  vdw-df-ahcx
        ! igcx=33:  vdw-df2-AH
        ! igcx=47:  vdw-df2-ahbr
        !
        iflag = 0
        IF ( igcx == 32) THEN ! vdW-DF-ahcx
           CALL cx13( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
           CALL cx13( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
           iflag = 3 ! for cx13 - sr hole
        ELSEIF ( igcx == 33) THEN ! vdW-DF2-ah
           CALL rPW86( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
           CALL rPW86( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
           iflag = 4 ! for rPW86 - sr hole
        ELSEIF ( igcx == 47) THEN ! vdW-DF2-ahbr
           CALL b86b( rho_up, grho2_up, 3, sx_up, v1x_up, v2x_up ) 
           CALL b86b( rho_dw, grho2_dw, 3, sx_dw, v1x_dw, v2x_dw ) 
           iflag = 6 ! for test-reserve - sr hole
        ENDIF
        !
        IF ( iflag == 0) THEN
          in_err = 5   ! Sorting vdW-DF-AHs failed
        ELSE
          sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
          v2x_up = 2.0_DP * v2x_up
          v2x_dw = 2.0_DP * v2x_dw
        ENDIF
        !
        IF ( exx_started ) THEN
           !
           CALL axsr( iflag, rho_up, grho2_up, sxsr_up, v1xsr_up, &
                                          v2xsr_up, screening_parameter, in_err )
           CALL axsr( iflag, rho_dw, grho2_dw, sxsr_dw, v1xsr_dw, &
                                          v2xsr_dw, screening_parameter, in_err )
           !
           sx_tot(ir) = sx_tot(ir) - exx_fraction*0.5_DP * ( sxsr_up*rnull_up + &
                                                             sxsr_dw*rnull_dw )
           v1x_up = v1x_up - exx_fraction * v1xsr_up
           v1x_dw = v1x_dw - exx_fraction * v1xsr_dw
           v2x_up = v2x_up - exx_fraction * v2xsr_up * 2.0_DP
           v2x_dw = v2x_dw - exx_fraction * v2xsr_dw * 2.0_DP
        END IF
        !
     CASE( 40 )                  ! 'c090 for vdw-df-c090' etc
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL c09x( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL c09x( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )  
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x_up = (1.0_DP - exx_fraction) * v1x_up
           v1x_dw = (1.0_DP - exx_fraction) * v1x_dw
           v2x_up = (1.0_DP - exx_fraction) * v2x_up
           v2x_dw = (1.0_DP - exx_fraction) * v2x_dw
        ENDIF
        !
     CASE( 41 )                 ! B86X for B86BPBEX hybrid
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL becke86b( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL becke86b( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x_up = (1.0_DP - exx_fraction) * v1x_up
           v1x_dw = (1.0_DP - exx_fraction) * v1x_dw
           v2x_up = (1.0_DP - exx_fraction) * v2x_up
           v2x_dw = (1.0_DP - exx_fraction) * v2x_dw
        ENDIF
        !
     CASE( 42 )                ! B88X for BHANDHLYP
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL becke88( rho_up, grho2_up, sx_up, v1x_up, v2x_up )
        CALL becke88( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw )
        !
        sx_tot(ir) = 0.5_DP * ( sx_up*rnull_up + sx_dw*rnull_dw )
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
        IF ( exx_started ) THEN
           sx_tot(ir) = (1.0_DP - exx_fraction) * sx_tot(ir)
           v1x_up = (1.0_DP - exx_fraction) * v1x_up
           v1x_dw = (1.0_DP - exx_fraction) * v1x_dw
           v2x_up = (1.0_DP - exx_fraction) * v2x_up
           v2x_dw = (1.0_DP - exx_fraction) * v2x_dw
        ENDIF
        !
        ! case igcx == 5 (HCTH) and 6 (OPTX) not implemented
        ! case igcx == 7 (meta-GGA) must be treated in a separate call to another
        ! routine: needs kinetic energy density in addition to rho and grad rho
        !
     CASE( 43 )                ! BEEX
        !
        rho_up = 2.0_DP * rho_up     ; rho_dw = 2.0_DP * rho_dw
        grho2_up = 4.0_DP * grho2_up ; grho2_dw = 4.0_DP * grho2_dw
        !
        CALL beefx( rho_up, grho2_up, sx_up, v1x_up, v2x_up, 0 )
        CALL beefx( rho_dw, grho2_dw, sx_dw, v1x_dw, v2x_dw, 0 )
        !
        sx_tot(ir) = 0.5_DP * (sx_up*rnull_up + sx_dw*rnull_dw)
        v2x_up = 2.0_DP * v2x_up
        v2x_dw = 2.0_DP * v2x_dw
        !
     CASE DEFAULT
        !
        sx_tot(ir) = 0.0_DP
        v1x_up = 0.0_DP ; v1x_dw = 0.0_DP
        v2x_up = 0.0_DP ; v2x_dw = 0.0_DP
        !
     END SELECT
     !
     IF (in_err/=0) THEN
#if defined(_OPENACC)
!$acc atomic write
#else
!$omp atomic write
#endif
       err_out = in_err
     ENDIF
     !
     v1x_out(ir,1) = v1x_up * rnull_up
     v1x_out(ir,2) = v1x_dw * rnull_dw
     v2x_out(ir,1) = v2x_up * rnull_up
     v2x_out(ir,2) = v2x_dw * rnull_dw
     !
  ENDDO
#if defined(_OPENACC)
!$acc end data
#else
!$omp end do
!$omp end parallel
#endif
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
  USE dft_setting_params,   ONLY: igcx, igcc, rho_threshold_gga
  USE corr_gga
  USE beef_interface, ONLY: beeflocalcorrspin
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
  INTEGER :: ir
  REAL(DP) :: rho, zeta, grho
  REAL(DP) :: sc, v1c_up, v1c_dw, v2c
  !REAL(DP), PARAMETER :: small=1.E-10_DP !, epsr=1.E-6_DP
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
#if defined(_OPENACC)
!$acc data present( rho_in, zeta_io, grho_in, sc_out, v1c_out, v2c_out )
!$acc parallel loop
#else
!$omp parallel if(ntids==1) default(none) &
!$omp private( rho, zeta, grho, sc, v1c_up, v1c_dw, v2c ) &
!$omp shared( igcc, sc_out, v1c_out, v2c_out, &
!$omp         rho_threshold_gga, zeta_io, length, &
!$omp         grho_in, rho_in )
!$omp do
#endif
  DO ir = 1, length
    !
    rho  = rho_in(ir)
    grho = grho_in(ir)
    IF ( ABS(zeta_io(ir))<=1.0_DP ) zeta_io(ir) = SIGN( MIN(ABS(zeta_io(ir)), &
                                    (1.0_DP-rho_threshold_gga)), zeta_io(ir) )
    zeta = zeta_io(ir)
    !
    IF ( ABS(zeta)>1.0_DP .OR. rho<=rho_threshold_gga .OR. &
         SQRT(ABS(grho))<=rho_threshold_gga ) THEN
       sc_out(ir) = 0.0_DP
       v1c_out(ir,1) = 0.0_DP ; v2c_out(ir) = 0.0_DP
       v1c_out(ir,2) = 0.0_DP
       CYCLE
    ENDIF
    !
    SELECT CASE( igcc )
    CASE( 1 )
       !
       CALL perdew86_spin( rho, zeta, grho, sc, v1c_up, v1c_dw, v2c )
       !
    CASE( 2 )
       !
       CALL ggac_spin( rho, zeta, grho, sc, v1c_up, v1c_dw, v2c )
       !
    CASE( 4 )
       !
       CALL pbec_spin( rho, zeta, grho, 1, sc, v1c_up, v1c_dw, v2c )
       !
    CASE( 8 )
       !
       CALL pbec_spin( rho, zeta, grho, 2, sc, v1c_up, v1c_dw, v2c )
       !
    CASE( 14 )
       !  
       CALL beeflocalcorrspin( rho, zeta, grho, sc, v1c_up, v1c_dw, v2c, 0 )
       !
    CASE DEFAULT
       !
       sc = 0.0_DP
       v1c_up = 0.0_DP
       v1c_dw = 0.0_DP
       v2c = 0.0_DP
       !
    END SELECT
    !
    sc_out(ir)  = sc
    v1c_out(ir,1) = v1c_up
    v1c_out(ir,2) = v1c_dw
    v2c_out(ir) = v2c
    !
  ENDDO
#if defined(_OPENACC)
!$acc end data
#else
!$omp end do
!$omp end parallel
#endif
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
  USE dft_setting_params,   ONLY: igcx, igcc, rho_threshold_gga,     &
                                  exx_started
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
  !! correlation potential (off-diag. term)
  !
  ! ... local variables
  !
  INTEGER :: ir
  REAL(DP) :: rho_up, rho_dw, grho_up, grho_dw
  REAL(DP) :: grho_ud
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif    
  !
#if defined(_OPENACC) 
!$acc data present( rho_in, grho_in, grho_ud_in, sc, v1c, v2c, v2c_ud )
!$acc parallel loop
#else 
!$omp parallel if(ntids==1) default(none) &
!$omp private( rho_up, rho_dw, grho_up, grho_dw, grho_ud ) &
!$omp shared( length, rho_in, grho_in, grho_ud_in, &
!$omp         rho_threshold_gga, sc, exx_started, &
!$omp         igcc, v1c, v2c, v2c_ud)
!$omp do
#endif
  DO ir = 1, length
    !
    rho_up = rho_in(ir,1)
    rho_dw = rho_in(ir,2)
    grho_up = grho_in(ir,1)
    grho_dw = grho_in(ir,2)
    grho_ud = grho_ud_in(ir)
    !
    IF ( rho_up+rho_dw < rho_threshold_gga ) THEN
       sc(ir) = 0.0_DP
       v1c(ir,1) = 0.0_DP ; v1c(ir,2) = 0.0_DP
       v2c(ir,1) = 0.0_DP ; v2c_ud(ir) = 0.0_DP
       v2c(ir,2) = 0.0_DP
       CYCLE
    ENDIF
    !
    CALL lsd_glyp( rho_up, rho_dw, grho_up, grho_dw, grho_ud, &
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
          v1c(ir,1) = 0.81_DP * v1c(ir,1)
          v1c(ir,2) = 0.81_DP * v1c(ir,2)
          v2c(ir,1) = 0.81_DP * v2c(ir,1)
          v2c(ir,2) = 0.81_DP * v2c(ir,2)
          v2c_ud(ir) = 0.81_DP * v2c_ud(ir)
       ENDIF
       !
    CASE( 13 )
       !
       IF ( exx_started ) THEN
          sc(ir) = 0.871_DP * sc(ir)
          v1c(ir,1) = 0.871_DP * v1c(ir,1)
          v1c(ir,2) = 0.871_DP * v1c(ir,2)
          v2c(ir,1) = 0.871_DP * v2c(ir,1)
          v2c(ir,2) = 0.871_DP * v2c(ir,2)
          v2c_ud(ir) = 0.871_DP * v2c_ud(ir)
       ENDIF
       !
    CASE DEFAULT
       !
       ! ... void
       !
    END SELECT
    !
  ENDDO
#if defined(_OPENACC)
!$acc end data
#else
!$omp end do
!$omp end parallel
#endif
  !
  RETURN
  !
END SUBROUTINE gcc_spin_more
!
!
END MODULE qe_drivers_gga
