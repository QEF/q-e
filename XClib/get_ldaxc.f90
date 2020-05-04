!
! ... LDAlib
!
SUBROUTINE get_ldaxclib( iexch_, icorr_ )
   !
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iexch_, icorr_
   !
   iexch = iexch_
   icorr = icorr_
   !
   RETURN
   !
END SUBROUTINE


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


SUBROUTINE get_threshold( rho_threshold_ )
   !
   USE kind_l, ONLY: DP
   USE dft_par_mod
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: rho_threshold_
   !
   rho_threshold = rho_threshold_
   !
   RETURN
   !
END SUBROUTINE
