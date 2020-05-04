!
!  ... LDAlib
!
MODULE dft_par_mod
    !
    USE kind_l, ONLY: DP
    !
    IMPLICIT NONE
    !
    INTEGER  :: iexch, icorr
    LOGICAL  :: exx_started !, is_there_finite_size_corr
    REAL(DP) :: finite_size_cell_volume
    REAL(DP) :: exx_fraction = 0.0_DP
    REAL(DP) :: rho_threshold = 1.E-10_DP
    !
END MODULE
