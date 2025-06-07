!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!========================================================================
!                      LDA-LSDA DRIVERS for E and V
!========================================================================
!
!-------------------------------------------------------------------------
MODULE qe_drivers_lda_lsda
  !-----------------------------------------------------------------------
  !! Contains the LDA drivers of QE that calculate XC energy and potential.
  !
  USE kind_l,               ONLY: DP
  USE exch_lda
  USE corr_lda
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: xc_lda, xc_lsda
  !
  !
CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE xc_lda( length, rho_in, ex_out, ec_out, vx_out, vc_out )
  !--------------------------------------------------------------------------
  !! LDA exchange and correlation functionals - Hartree a.u.
  !
  !! * Exchange:
  !!    * Slater;
  !!    * relativistic Slater.
  !! * Correlation:
  !!    * Ceperley-Alder (Perdew-Zunger parameters);
  !!    * Vosko-Wilk-Nusair;
  !!    * Lee-Yang-Parr;
  !!    * Perdew-Wang;
  !!    * Wigner;
  !!    * Hedin-Lundqvist;
  !!    * Ortiz-Ballone (Perdew-Zunger formula);
  !!    * Ortiz-Ballone (Perdew-Wang formula);
  !!    * Gunnarsson-Lundqvist.
  !
  !! NOTE:
  !! $$ E_x = \int E_x(\text{rho}) dr, E_x(\text{rho}) = 
  !!               \text{rho}\epsilon_c(\text{rho})\ . $$
  !! Same for correlation.
  !
  USE dft_setting_params,   ONLY: iexch, icorr, rho_threshold_lda, exx_started, &
                                  exx_fraction, finite_size_cell_volume
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! Charge density
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ex_out
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vx_out
  !! \(dE_x(\text{rho})/d\text{rho}\)  ( NOT 
  !! \(d\epsilon_x(\text{rho})/d\text{rho}\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec_out
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vc_out
  !! \(dE_c(\text{rho})/d\text{rho}\)  ( NOT 
  !! \(d\epsilon_c(\text{rho})/d\text{rho}\) )
  !
  ! ... local variables
  !
  INTEGER  :: ir
  REAL(DP) :: rho, rs
  REAL(DP) :: ex, ec, ec_
  REAL(DP) :: vx, vc, vc_
  REAL(DP), PARAMETER :: third = 1.0_DP/3.0_DP, &
                         pi34 = 0.6203504908994_DP, e2 = 2.0_DP
  !                      pi34 = (3/4pi)^(1/3)
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
#if defined(_OPENACC)
!$acc data present( rho_in, ex_out, vx_out, ec_out, vc_out )
!$acc parallel loop
#else
!$omp parallel if(ntids==1) default(none) &
!$omp private( rho, rs, ex, ec, ec_, vx, vc, vc_ ) &
!$omp shared( rho_in, length, iexch, icorr, ex_out, ec_out, vx_out, vc_out, &
!$omp         finite_size_cell_volume, exx_fraction, exx_started, &
!$omp         rho_threshold_lda )
!$omp do
#endif
  DO ir = 1, length
     !
     rho = ABS(rho_in(ir))
     !
     ! ... RHO THRESHOLD
     !
     IF ( rho > rho_threshold_lda ) THEN
        rs = pi34 / rho**third
     ELSE
        ex_out(ir) = 0.0_DP  ;  ec_out(ir) = 0.0_DP
        vx_out(ir) = 0.0_DP  ;  vc_out(ir) = 0.0_DP
        CYCLE
     ENDIF
     !
     ! ... EXCHANGE
     !
     SELECT CASE( iexch )
     CASE( 1 )                      ! 'sla'
        !
        CALL slater( rs, ex, vx )
        !
     CASE( 2 )                      ! 'sl1'
        !
        CALL slater1( rs, ex, vx )
        !
     CASE( 3 )                      ! 'rxc'
        !
        CALL slater_rxc( rs, ex, vx )
        !
     CASE( 4, 5 )                   ! 'oep','hf'
        !
        IF ( exx_started ) THEN
           ex = 0.0_DP
           vx = 0.0_DP
        ELSE
           CALL slater( rs, ex, vx )
        ENDIF
        !
     CASE( 6, 7 )                   ! 'pb0x' or 'DF-cx-0', or 'DF2-0',
        !                           ! 'B3LYP'
        CALL slater( rs, ex, vx )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE( 8 )                      ! 'sla+kzk'
        !
        CALL slaterKZK( rs, ex, vx, finite_size_cell_volume )
        !
     CASE( 9 )                      ! 'X3LYP'
        !
        CALL slater( rs, ex, vx )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE DEFAULT
        !
        ex = 0.0_DP
        vx = 0.0_DP
        !
     END SELECT
     !
     !
     ! ... CORRELATION
     !
     SELECT CASE( icorr )
     CASE( 1 )
        !
        CALL pz( rs, 1, ec, vc )
        !
     CASE( 2 )
        !
        CALL vwn( rs, ec, vc )
        !
     CASE( 3 )
        !
        CALL lyp( rs, ec, vc )
        !
     CASE( 4 )
        !
        CALL pw( rs, 1, ec, vc )
        !
     CASE( 5 )
        !
        CALL wignerc( rs, ec, vc )
        !
     CASE( 6 )
        !
        CALL hl( rs, ec, vc )
        !
     CASE( 7 )
        !
        CALL pz( rs, 2, ec, vc )
        ! 
     CASE( 8 )
        !
        CALL pw( rs, 2, ec, vc )
        !
     CASE( 9 )
        !
        CALL gl( rs, ec, vc )
        !
     CASE( 10 )
        !
        CALL pzKZK( rs, ec, vc, finite_size_cell_volume )
        !
     CASE( 11 )
        !
        CALL vwn1_rpa( rs, ec, vc )
        !
     CASE( 12 )                ! 'B3LYP'
        !
        CALL vwn( rs, ec, vc )
        ec = 0.19_DP * ec
        vc = 0.19_DP * vc
        !
        CALL lyp( rs, ec_, vc_ )
        ec = ec + 0.81_DP * ec_
        vc = vc + 0.81_DP * vc_
        !
     CASE( 13 )                ! 'B3LYP-V1R'
        !
        CALL vwn1_rpa( rs, ec, vc )
        ec = 0.19_DP * ec
        vc = 0.19_DP * vc
        !
        CALL lyp( rs, ec_, vc_ )
        ec = ec + 0.81_DP * ec_
        vc = vc + 0.81_DP * vc_
        !
     CASE( 14 )                ! 'X3LYP'
        !
        CALL vwn1_rpa( rs, ec, vc )
        ec = 0.129_DP * ec
        vc = 0.129_DP * vc
        !
        CALL lyp( rs, ec_, vc_ )
        ec = ec + 0.871_DP * ec_
        vc = vc + 0.871_DP * vc_
        !
     CASE DEFAULT
        !
        ec = 0.0_DP
        vc = 0.0_DP
        !
     END SELECT
     !
     ex_out(ir) = ex  ;  ec_out(ir) = ec
     vx_out(ir) = vx  ;  vc_out(ir) = vc
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
END SUBROUTINE xc_lda
!
!
!-----------------------------------------------------------------------------
SUBROUTINE xc_lsda( length, rho_in, zeta_in, ex_out, ec_out, vx_out, vc_out )
  !-----------------------------------------------------------------------------
  !! LSD exchange and correlation functionals - Hartree a.u.
  !
  !! * Exchange:
  !!    * Slater (alpha=2/3).
  !! * Correlation:
  !!    * Ceperley & Alder (Perdew-Zunger parameters);
  !!    * Perdew & Wang.
  !
  USE dft_setting_params,   ONLY: iexch, icorr, rho_threshold_lda, exx_started, &
                                  exx_fraction
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! Total charge density
  REAL(DP), INTENT(IN),  DIMENSION(length) :: zeta_in
  !! zeta = mag / rho_tot
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ex_out
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec_out
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vx_out
  !! \(dE_x(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_x(\text{rho})/d\text{rho}\) )
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vc_out
  !! \(dE_c(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_c(\text{rho})/d\text{rho}\) )
  !
  ! ...  local variables
  !
  INTEGER  :: ir
  REAL(DP) :: rho, rs, zeta
  REAL(DP) :: ex, ec, ec_
  REAL(DP) :: vx_up, vc_up, vc_up_
  REAL(DP) :: vx_dw, vc_dw, vc_dw_
  !
  REAL(DP), PARAMETER :: third = 1.0_DP/3.0_DP, &
                         pi34 = 0.6203504908994_DP
  !                      pi34 = (3/4pi)^(1/3)
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
#if defined(_OPENACC)  
!$acc data present( rho_in, zeta_in, ex_out, vx_out, ec_out, vc_out )
!$acc parallel loop  
#else
!$omp parallel if(ntids==1) default(none) &
!$omp private( rho, rs, zeta, ex, ec, ec_, vx_up, vx_dw, vc_up, &
!$omp          vc_dw, vc_up_, vc_dw_ ) &
!$omp shared( length, iexch, icorr, exx_fraction, &
!$omp         vx_out, vc_out, ex_out, ec_out, &
!$omp         zeta_in, exx_started, rho_in, rho_threshold_lda )
!$omp do
#endif
  DO ir = 1, length
     !
     zeta = zeta_in(ir)
     IF (ABS(zeta) > 1.D0) zeta = SIGN( 1.D0, zeta )
     !
     rho = ABS(rho_in(ir))
     !
     IF ( rho > rho_threshold_lda ) THEN
        rs = pi34 / rho**third
     ELSE
        ex_out(ir) = 0.0_DP  ;  vx_out(ir,1) = 0.0_DP  ;  vx_out(ir,2) = 0.0_DP
        ec_out(ir) = 0.0_DP  ;  vc_out(ir,1) = 0.0_DP  ;  vc_out(ir,2) = 0.0_DP
        CYCLE
     ENDIF
     !
     !
     ! ... EXCHANGE
     !
     SELECT CASE( iexch )
     CASE( 1 )                                      ! 'sla'
        !
        CALL slater_spin( rho, zeta, ex, vx_up, vx_dw )
        !
     CASE( 2 )                                      ! 'sl1'
        !
        CALL slater1_spin( rho, zeta, ex, vx_up, vx_dw )
        !
     CASE( 3 )                                      ! 'rxc'
        !
        CALL slater_rxc_spin( rho, zeta, ex, vx_up, vx_dw )
        !
     CASE( 4, 5 )                                   ! 'oep','hf'
        !
        IF ( exx_started ) THEN
           ex = 0.0_DP
           vx_up = 0.0_DP ; vx_dw = 0.0_DP
        ELSE
           CALL slater_spin( rho, zeta, ex, vx_up, vx_dw )
        ENDIF
        !
     CASE( 6 )                                      ! 'pb0x'
        !
        CALL slater_spin( rho, zeta, ex, vx_up, vx_dw )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx_up = (1.0_DP - exx_fraction) * vx_up
           vx_dw = (1.0_DP - exx_fraction) * vx_dw
        ENDIF
        !
     CASE( 7 )                                      ! 'B3LYP'
        !
        CALL slater_spin( rho, zeta, ex, vx_up, vx_dw )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx_up = (1.0_DP - exx_fraction) * vx_up
           vx_dw = (1.0_DP - exx_fraction) * vx_dw
        ENDIF
        !
     CASE( 9 )                                      ! 'X3LYP'
        !
        CALL slater_spin( rho, zeta, ex, vx_up, vx_dw )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx_up = (1.0_DP - exx_fraction) * vx_up
           vx_dw = (1.0_DP - exx_fraction) * vx_dw
        ENDIF
        !
     CASE DEFAULT
        !
        ex = 0.0_DP
        vx_up = 0.0_DP
        vx_dw = 0.0_DP
        !
     END SELECT
     !
     !
     ! ... CORRELATION
     !
     SELECT CASE( icorr )
     CASE( 1 )
        !
        CALL pz_spin( rs, zeta, ec, vc_up, vc_dw )
        !
     CASE( 2 )
        !
        CALL vwn_spin( rs, zeta, ec, vc_up, vc_dw )
        !
     CASE( 3 )
        !
        CALL lsd_lyp( rho, zeta, ec, vc_up, vc_dw )       ! from CP/FPMD
        !
     CASE( 4 )
        !
        CALL pw_spin( rs, zeta, ec, vc_up, vc_dw )
        !
     CASE( 12 )                                           ! 'B3LYP'
        !
        CALL vwn_spin( rs, zeta, ec, vc_up, vc_dw )
        ec = 0.19_DP * ec
        vc_up = 0.19_DP * vc_up
        vc_dw = 0.19_DP * vc_dw
        !
        CALL lsd_lyp( rho, zeta, ec_, vc_up_, vc_dw_ )    ! from CP/FPMD
        ec = ec + 0.81_DP * ec_
        vc_up = vc_up + 0.81_DP * vc_up_
        vc_dw = vc_dw + 0.81_DP * vc_dw_
        !     
     CASE( 13 )                                           ! 'B3LYP-V1R'
        !
        CALL vwn1_rpa_spin( rs, zeta, ec, vc_up, vc_dw )
        ec = 0.19_DP * ec
        vc_up = 0.19_DP * vc_up
        vc_dw = 0.19_DP * vc_dw
        !
        CALL lsd_lyp( rho, zeta, ec_, vc_up_, vc_dw_ )    ! from CP/FPMD
        ec = ec + 0.81_DP * ec_
        vc_up = vc_up + 0.81_DP * vc_up_
        vc_dw = vc_dw + 0.81_DP * vc_dw_
        !
     CASE( 14 )                                           ! 'X3LYP
        !
        CALL vwn1_rpa_spin( rs, zeta, ec, vc_up, vc_dw )
        ec = 0.129_DP * ec
        vc_up = 0.129_DP * vc_up
        vc_dw = 0.129_DP * vc_dw
        !
        CALL lsd_lyp( rho, zeta, ec_, vc_up_, vc_dw_ )    ! from CP/FPMD
        ec = ec + 0.871_DP * ec_
        vc_up = vc_up + 0.871_DP * vc_up_
        vc_dw = vc_dw + 0.871_DP * vc_dw_
        !
     CASE DEFAULT
        !
        ec = 0.0_DP
        vc_up = 0.0_DP
        vc_dw = 0.0_DP
        !
     END SELECT
     !
     ex_out(ir) = ex  ;  vx_out(ir,1) = vx_up  ;  vx_out(ir,2) = vx_dw 
     ec_out(ir) = ec  ;  vc_out(ir,1) = vc_up  ;  vc_out(ir,2) = vc_dw
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
END SUBROUTINE xc_lsda
!
END MODULE qe_drivers_lda_lsda
