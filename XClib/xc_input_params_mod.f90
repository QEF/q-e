!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE dft_par_mod
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
    CHARACTER(LEN=32) :: dft = 'not set'
    !! Full name of the XC functional
    INTEGER, PARAMETER :: notset = -1
    !! Value of indexes that have not been set yet
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
    INTEGER  :: n_ext_params(6) = 0._DP
    !! number of external parameters for each functional
    REAL(DP) :: par_list(6,10)
    !! list of external parameters
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
    LOGICAL :: scan_exx    = .FALSE.
    !! TRUE if SCAN0 functional is active
    !
    LOGICAL :: discard_input_dft = .FALSE.
    !! TRUE if input DFT can be overwritten
    INTEGER :: beeftype = -1
    !! Index for BEEF functional
    INTEGER :: beefvdw = 0
    !! Index for vdw term of BEEF
    !
    INTEGER, PARAMETER :: nxc=8, ncc=10, ngcx=46, ngcc=13, nmeta=6
    CHARACTER(LEN=4) :: exc, corr, gradx, gradc, meta
    DIMENSION :: exc(0:nxc), corr(0:ncc), gradx(0:ngcx), gradc(0:ngcc), &
                 meta(0:nmeta)
    !
    DATA exc  / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', 'B3LP', 'KZK' /
    DATA corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
                'OBW', 'GL' , 'KZK' /
    !
    DATA gradx / 'NOGX', 'B88',  'GGX',  'PBX',  'RPB',  'HCTH', 'OPTX', &
                 'xxxx', 'PB0X', 'B3LP', 'PSX',  'WCX',  'HSE',  'RW86', 'PBE', &
                 'xxxx', 'C09X', 'SOX',  'xxxx', 'Q2DX', 'GAUP', 'PW86', 'B86B', &
                 'OBK8', 'OB86', 'EVX',  'B86R', 'CX13', 'X3LP', &
                 'CX0',  'R860', 'CX0P', 'AHCX', 'AHF2', &
                 'AHPB', 'AHPS', 'CX14', 'CX15', 'BR0',  'CX16', 'C090', &
                 'B86X', 'B88X', 'BEEX', 'RPBX', 'W31X', 'W32X' /
    !
    DATA gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH', 'NONE',&
                 'B3LP', 'PSC', 'PBE', 'xxxx', 'xxxx', 'Q2DC', 'BEEC' /
    !
    DATA meta  / 'NONE', 'TPSS', 'M06L', 'TB09', 'META', 'SCAN', 'SCA0' /
    !
    !
END MODULE
