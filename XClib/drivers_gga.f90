!-----------------------------------------------------------------------
!------- GRADIENT CORRECTION DRIVERS ----------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE gcxc_l( length, rho_in, grho_in, sx_out, sc_out, v1x_out, &
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
  USE kind_l, ONLY: DP
  USE dft_par_mod
  USE exch_gga_l
  USE corr_gga_l
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
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, grho, sx, sx_, sxsr, v1x, v1x_, v1xsr, &
!$omp             v2x, v2x_, v2xsr, sc, v1c, v2c )
  DO ir = 1, length  
     !
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
        CALL becke88_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 2 )
        !
        CALL ggax_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 3 )
        !
        CALL pbex_l( rho, grho, 1, sx, v1x, v2x )
        !
     CASE( 4 )
        !
        CALL pbex_l( rho, grho, 2, sx, v1x, v2x )
        !
     CASE( 5 )
        !
        IF (igcc == 5) CALL hcth_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 6 )
        !
        CALL optx_l( rho, grho, sx, v1x, v2x )
        !
     ! case igcx == 7 (meta-GGA) must be treated in a separate call to another
     ! routine: needs kinetic energy density in addition to rho and grad rho
     CASE( 8 ) ! 'PBE0'
        !
        CALL pbex_l( rho, grho, 1, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 9 ) ! 'B3LYP'
        !
        CALL becke88_l( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = 0.72_DP * sx
           v1x = 0.72_DP * v1x
           v2x = 0.72_DP * v2x
        ENDIF
        !
     CASE( 10 ) ! 'pbesol'
        !
        CALL pbex_l( rho, grho, 3, sx, v1x, v2x )
        !
     CASE( 11 ) ! 'wc'
        !
        CALL wcx_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 12 ) ! 'pbexsr'
        !
        CALL pbex_l( rho, grho, 1, sx, v1x, v2x )
        !
        IF (exx_started) THEN
          CALL pbexsr_l( rho, grho, sxsr, v1xsr, v2xsr, screening_parameter )
          sx  = sx  - exx_fraction * sxsr
          v1x = v1x - exx_fraction * v1xsr
          v2x = v2x - exx_fraction * v2xsr
        ENDIF
        !
     CASE( 13 ) ! 'rPW86'
        !
        CALL rPW86_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 16 ) ! 'C09x'
        !
        CALL c09x_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 17 ) ! 'sogga'
        !
        CALL sogga_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 19 ) ! 'pbeq2d'
        !
        CALL pbex_l( rho, grho, 4, sx, v1x, v2x )
        !
     CASE( 20 ) ! 'gau-pbe'
        !
        CALL pbex_l( rho, grho, 1, sx, v1x, v2x )
        IF (exx_started) THEN
          CALL pbexgau_l( rho, grho, sxsr, v1xsr, v2xsr, gau_parameter )
          sx  = sx  - exx_fraction * sxsr
          v1x = v1x - exx_fraction * v1xsr
          v2x = v2x - exx_fraction * v2xsr
        ENDIF
        !
     CASE( 21 ) ! 'pw86'
        !
        CALL pw86_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 22 ) ! 'b86b'
        !
        CALL becke86b_l( rho, grho, sx, v1x, v2x )
        ! CALL b86b( rho, grho, 1, sx, v1x, v2x )
        !
     CASE( 23 ) ! 'optB88'
        !
        CALL pbex_l( rho, grho, 5, sx, v1x, v2x )
        !
     CASE( 24 ) ! 'optB86b'
        !
        CALL pbex_l( rho, grho, 6, sx, v1x, v2x )
        ! CALL b86b (rho, grho, 2, sx, v1x, v2x)
        !
     CASE( 25 ) ! 'ev93'
        !
        CALL pbex_l( rho, grho, 7, sx, v1x, v2x )
        !
     CASE( 26 ) ! 'b86r'
        !
        CALL b86b_l( rho, grho, 3, sx, v1x, v2x )
        !
     CASE( 27 ) ! 'cx13'
        !
        CALL cx13_l( rho, grho, sx, v1x, v2x )
        !
     CASE( 28 ) ! 'X3LYP'
        !
        CALL becke88_l( rho, grho, sx, v1x, v2x )
        CALL pbex_l( rho, grho, 1, sx_, v1x_, v2x_ )
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
        CALL cx13_l( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 30 ) ! 'r860'
        !
        CALL rPW86_l( rho, grho, sx, v1x, v2x )
        !
        IF (exx_started) then
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 38 ) ! 'BR0'
        !
        CALL b86b_l( rho, grho, 3, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 40 ) ! 'c090'
        !
        CALL c09x_l( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 41 ) ! 'B86BPBEX'
        !
        CALL becke86b_l( rho, grho, sx, v1x, v2x )
        IF (exx_started) THEN
           sx  = (1.0_DP - exx_fraction) * sx
           v1x = (1.0_DP - exx_fraction) * v1x
           v2x = (1.0_DP - exx_fraction) * v2x
        ENDIF
        !
     CASE( 42 ) ! 'BHANDHLYP'
        !
        CALL becke88_l( rho, grho, sx, v1x, v2x )
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
        CALL pbex_l( rho, grho, 8, sx, v1x, v2x )
        !
     CASE( 45 ) ! 'W31X'
        !
        CALL pbex_l( rho, grho, 9, sx, v1x, v2x )
        !
     CASE( 46 ) ! 'W32X'
        !
        CALL b86b_l( rho, grho, 4, sx, v1x, v2x )
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
        CALL perdew86_l( rho, grho, sc, v1c, v2c )
        !
     CASE( 2 )
        !
        CALL ggac_l( rho, grho, sc, v1c, v2c )
        !
     CASE( 3 )
        !
        CALL glyp_l( rho, grho, sc, v1c, v2c )
        !
     CASE( 4 )
        !
        CALL pbec_l( rho, grho, 1, sc, v1c, v2c )
        !
     ! igcc == 5 (HCTH) is calculated together with case igcx=5
     ! igcc == 6 (meta-GGA) is treated in a different routine
     CASE( 7 ) !'B3LYP'
        !
        CALL glyp_l( rho, grho, sc, v1c, v2c )
        IF (exx_started) THEN
           sc  = 0.81_DP * sc
           v1c = 0.81_DP * v1c
           v2c = 0.81_DP * v2c
        ENDIF
        !
     CASE( 8 ) ! 'PBEsol'
        !
        CALL pbec_l( rho, grho, 2, sc, v1c, v2c )
        !
     ! igcc ==  9 set to 5, back-compatibility
     ! igcc == 10 set to 6, back-compatibility
     ! igcc == 11 M06L calculated in another routine
     CASE( 12 ) ! 'Q2D'
        !
        CALL pbec_l( rho, grho, 3, sc, v1c, v2c )
        !
     CASE( 13 ) !'X3LYP'
        !
        CALL glyp_l( rho, grho, sc, v1c, v2c )
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
END SUBROUTINE gcxc_l
