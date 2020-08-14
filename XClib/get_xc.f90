!
!
SUBROUTINE get_xclib_IDs( iexch_, icorr_, igcx_, igcc_, imeta_, is_libxc_ )
   !
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iexch_, icorr_
   INTEGER, INTENT(IN) :: igcx_, igcc_
   INTEGER, INTENT(IN) :: imeta_
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
!  LDA 
!
SUBROUTINE get_ldaxcparlib( finite_size_cell_volume_, exx_started_, exx_fraction_ )
   !
   USE kind_l, ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   REAL(DP), OPTIONAL, INTENT(IN) :: finite_size_cell_volume_
   LOGICAL,  OPTIONAL, INTENT(IN) :: exx_started_
   REAL(DP), OPTIONAL, INTENT(IN) :: exx_fraction_
   !
   finite_size_cell_volume = -1.d0
   IF ( present(finite_size_cell_volume_) ) &
      finite_size_cell_volume = finite_size_cell_volume_
   exx_started = .FALSE.
   IF ( present(exx_started_) ) &
      exx_started = exx_started_
   exx_fraction = 0.d0
   IF ( present(exx_fraction_) ) &
      exx_fraction = exx_fraction_
   !
   RETURN
   !
END SUBROUTINE


SUBROUTINE get_ggaxcparlib( gau_scr_par_, exx_started_, exx_fraction_ )
   !
   USE kind_l, ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: gau_scr_par_
   LOGICAL, OPTIONAL, INTENT(IN) :: exx_started_
   REAL(DP), OPTIONAL, INTENT(IN) :: exx_fraction_
   !
   IF ( igcx == 12 ) &
      screening_parameter = gau_scr_par_
   IF ( igcx == 20 ) &
      gau_parameter = gau_scr_par_
   !
   !exx_started = .FALSE.
   IF ( present(exx_started_) ) &
      exx_started = exx_started_
   !exx_fraction = 0.d0
   IF ( present(exx_fraction_) ) THEN
      exx_fraction = 0.d0
      IF ( exx_started ) exx_fraction = exx_fraction_ 
   ENDIF
   !
   RETURN
   !
END SUBROUTINE


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
