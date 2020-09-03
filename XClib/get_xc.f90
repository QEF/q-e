!
!
SUBROUTINE get_IDs( iexch_, icorr_, igcx_, igcc_, imeta_, imetac_, is_libxc_ )
   !
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iexch_, icorr_
   INTEGER, INTENT(IN) :: igcx_, igcc_
   INTEGER, INTENT(IN) :: imeta_, imetac_
   LOGICAL, OPTIONAL   :: is_libxc_(6)
   !
   iexch = iexch_
   icorr = icorr_
   igcx  = igcx_
   igcc  = igcc_
   imeta = imeta_
   !
   IF ( PRESENT(is_libxc_) ) is_libxc_l = is_libxc_
   !
   RETURN
   !
END SUBROUTINE
!
!
SUBROUTINE get_exx_started_l( exx_started_ )
   !
   USE kind_l,  ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: exx_started_
   !
   exx_started = exx_started_
   !
   RETURN
   !
END SUBROUTINE
!
SUBROUTINE get_exx_fraction_l( exx_fraction_ )
   !
   USE kind_l,  ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: exx_fraction_
   !
   exx_fraction = exx_fraction_
   !
   RETURN
   !
END SUBROUTINE
!
!  LDA 
!
SUBROUTINE get_finite_size_cell_l( finite_size_cell_volume_ )
   !
   USE kind_l, ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: finite_size_cell_volume_
   !
   finite_size_cell_volume = finite_size_cell_volume_
   !
   RETURN
   !
END SUBROUTINE
!
!
SUBROUTINE get_gau_scr_par_l( gau_scr_par_ )
   !
   USE kind_l,  ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: gau_scr_par_
   !
   IF ( igcx == 12 ) &
      screening_parameter = gau_scr_par_
   IF ( igcx == 20 ) &
      gau_parameter = gau_scr_par_
   !
   RETURN
   !
END SUBROUTINE
!
!
SUBROUTINE set_threshold( fkind, rho_threshold_, grho_threshold_, tau_threshold_ )
   !
   USE kind_l, ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   CHARACTER(len=4), INTENT(IN) :: fkind
   REAL(DP), INTENT(IN) :: rho_threshold_
   REAL(DP), INTENT(IN), OPTIONAL :: grho_threshold_
   REAL(DP), INTENT(IN), OPTIONAL :: tau_threshold_
   !
   SELECT CASE( TRIM(fkind) )
   CASE( 'lda' )
     rho_threshold_lda = rho_threshold_
   CASE( 'gga' )
     rho_threshold_gga = rho_threshold_
     IF ( PRESENT(grho_threshold_) ) grho_threshold_gga = grho_threshold_
   CASE( 'mgga' )
     rho_threshold_gga = rho_threshold_
     IF ( PRESENT(grho_threshold_) ) grho2_threshold_mgga = grho_threshold_
     IF ( PRESENT(tau_threshold_)  ) tau_threshold_mgga   = tau_threshold_
   END SELECT
   !
  RETURN
   !
END SUBROUTINE set_threshold


 SUBROUTINE get_lda_threshold( rho_threshold_ )                           !---------cancella
    !
    USE kind_l, ONLY: DP
    USE dft_par_mod
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: rho_threshold_
    !
    rho_threshold_lda = rho_threshold_
    !
    RETURN
    !
 END SUBROUTINE


SUBROUTINE get_gga_threshold( rho_threshold_, grho_threshold_ )
   !
   USE kind_l, ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: rho_threshold_
   REAL(DP), INTENT(IN) :: grho_threshold_
   !
   rho_threshold_gga = rho_threshold_
   rho_threshold_gga = grho_threshold_
   !
   RETURN
   !
END SUBROUTINE


SUBROUTINE get_mgga_threshold( rho_threshold_, grho2_threshold_, tau_threshold_ )
   !
   USE kind_l, ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: rho_threshold_
   REAL(DP), INTENT(IN) :: grho2_threshold_
   REAL(DP), INTENT(IN) :: tau_threshold_
   !
   rho_threshold_mgga = rho_threshold_
   grho2_threshold_mgga = grho2_threshold_
   tau_threshold_mgga  = tau_threshold_
   !
   RETURN
   !
END SUBROUTINE
