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
    REAL(DP) :: finite_size_cell_volume, exx_fraction
    REAL(DP) :: rho_threshold
    !
END MODULE
