!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE dft_setting_params
    !----------------------------------------------------------------
    !! Parameters that define the XC functionals.
    !
    USE kind_l, ONLY: DP
#if defined(__LIBXC)
    USE xc_f03_lib_m
#endif
    !
    IMPLICIT NONE
    !
    SAVE
    !
    CHARACTER(LEN=32) :: dft = 'not set'
    !! Full name of the XC functional
    INTEGER, PARAMETER :: notset = -1
    !! Value of indexes that have not been set yet
    INTEGER, PARAMETER :: max_flags = 29
    !! maximum number of libxc flags per functional
    !
    LOGICAL :: is_libxc(6) = .FALSE.
    !! \(\text{is_libxc(i)}=TRUE\) if the i-th term of the input 
    !! functional is from Libxc
    !
    LOGICAL :: libxc_initialized(6) = .FALSE.
    !! TRUE if libxc functionals have been initialized
    !
#if defined(__LIBXC)
    TYPE(xc_f03_func_t) :: xc_func(6)
    !! pointers to libxc functional structs
    TYPE(xc_f03_func_info_t) :: xc_info(6)
    !! pointers to libxc info structs
    INTEGER :: libxc_flags(6,0:max_flags) = notset
    !! flags related to specific features of libxc functionals 
    !! (see Libxc doc)
    INTEGER :: n_ext_params(6) = 0
    !! number of external parameters for each functional
    REAL(DP) :: par_list(6,50)
    !! list of external parameters
    LOGICAL :: xc_kind_error = .FALSE.
    !! some Libxc dfts are currently not usable in QE (for example
    !! kinetic dfts). If one of them is detected this variable is
    !! set to TRUE.
    CHARACTER(LEN=100) :: lxc_exx_desc(5)
    DATA lxc_exx_desc / 'Fraction of HF exchange',            &
                        'Fraction of exact exchange',         &
                        'Fraction of Hartree-Fock exchange ', &
                        'Coefficient for exact exchange ',    &
                        'Mixing parameter '/
    !
    CHARACTER(LEN=100) :: lxc_scr_desc(5)
    DATA lxc_scr_desc / 'Range separation parameter', &
                        'Range-separation parameter', &
                        'Range separation constant',  &
                        'Range-separation constant',  &
                        'Screening parameter for HF' /
#endif
    !
    LOGICAL  :: exx_started = .FALSE.
    !! TRUE if Exact Exchange is active
    REAL(DP) :: exx_fraction = 0.0_DP
    !! Exact Exchange fraction parameter
    !
    INTEGER  :: iexch = notset
    !! LDA exchange index
    INTEGER  :: icorr = notset
    !! LDA correlation index
    REAL(DP) :: finite_size_cell_volume = -1._DP
    !! Cell volume parameter for finite size correction
    REAL(DP) :: rho_threshold_lda = 1.E-10_DP
    !! Threshold value for the density in LDA
    !
    INTEGER  :: igcx = notset
    !! GGA exchange index
    INTEGER  :: igcc = notset
    !! GGA correlation index
    REAL(DP) :: rho_threshold_gga = 1.E-6_DP
    !! Threshold value for the density in GGA
    REAL(DP) :: grho_threshold_gga = 1.E-10_DP
    !! Threshold value for the density gradient in GGA
    REAL(DP) :: screening_parameter
    !! Screening parameter for PBE exchange
    REAL(DP) :: gau_parameter
    !! Gaussian parameter for PBE exchange
    !
    INTEGER  :: imeta = notset
    !! MGGA exchange index
    INTEGER  :: imetac = notset
    !! MGGA correlation index
    REAL(DP) :: rho_threshold_mgga = 1.E-12_DP
    !! Threshold value for the density in MGGA
    REAL(DP) :: grho2_threshold_mgga = 1.E-24_DP
    !! Threshold value for the density gradient (square modulus) in MGGA
    REAL(DP) :: tau_threshold_mgga = 1.0E-12_DP
    !! Threshold value for the density laplacian in MGGA
    !
    ! 
    LOGICAL :: islda       = .FALSE.
    !! TRUE if the functional is LDA only
    LOGICAL :: isgradient  = .FALSE.
    !! TRUE if the functional is GGA
    LOGICAL :: has_finite_size_correction = .FALSE.
    !! TRUE if finite size correction is present
    LOGICAL :: finite_size_cell_volume_set = .FALSE.
    !! TRUE if the cell volume has been set for finite size correction.
    LOGICAL :: ismeta      = .FALSE.
    !! TRUE if the functional is MGGA
    LOGICAL :: ishybrid    = .FALSE.
    !! TRUE if the functional is hybrid
    INTEGER :: exx_term    = 0
    !! term of the functional that is hybrid
    !
    LOGICAL :: discard_input_dft = .FALSE.
    !! TRUE if input DFT can be overwritten
    INTEGER :: beeftype = -1
    !! Index for BEEF functional
    INTEGER :: beefvdw = 0
    !! Index for vdw term of BEEF
    !
END MODULE dft_setting_params
