!
!
MODULE dft_par_mod
    !
    USE kind_l, ONLY: DP
    !
    IMPLICIT NONE
    !
    LOGICAL  :: exx_started
    REAL(DP) :: exx_fraction = 0.0_DP
    !
    ! LDA
    INTEGER  :: iexch, icorr
    REAL(DP) :: finite_size_cell_volume
    REAL(DP) :: rho_threshold_lda = 1.E-10_DP
    !
    ! GGA
    INTEGER  :: igcx, igcc
    REAL(DP) :: rho_threshold_gga = 1.E-10_DP
    REAL(DP) :: grho_threshold_gga = 1.E-10_DP
    REAL(DP) :: screening_parameter, gau_parameter
    !
    ! MGGA
    INTEGER  :: imeta
    REAL(DP) :: rho_threshold_mgga = 1.E-12_DP
    REAL(DP) :: grho2_threshold_mgga = 1.E-24_DP
    REAL(DP) :: tau_threshold_mgga = 1.0E-12_DP
    !
END MODULE
