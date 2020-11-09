!
!
MODULE dft_par_mod
    !
    USE kind_l, ONLY: DP
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=25) :: dft = 'not set'
    !
    INTEGER, PARAMETER :: notset = -1
    !
    LOGICAL  :: is_libxc(6) = .FALSE.
    !
    LOGICAL  :: exx_started = .FALSE.
    REAL(DP) :: exx_fraction = 0.0_DP
    !
    ! LDA
    INTEGER  :: iexch = notset
    INTEGER  :: icorr = notset
    REAL(DP) :: finite_size_cell_volume = -1._DP
    REAL(DP) :: rho_threshold_lda = 1.E-10_DP
    !
    ! GGA
    INTEGER  :: igcx = notset
    INTEGER  :: igcc = notset
    REAL(DP) :: rho_threshold_gga = 1.E-6_DP
    REAL(DP) :: grho_threshold_gga = 1.E-10_DP
    REAL(DP) :: screening_parameter, gau_parameter
    !
    ! MGGA
    INTEGER  :: imeta = notset
    INTEGER  :: imetac = notset
    REAL(DP) :: rho_threshold_mgga = 1.E-12_DP
    REAL(DP) :: grho2_threshold_mgga = 1.E-24_DP
    REAL(DP) :: tau_threshold_mgga = 1.0E-12_DP
    !
    
    
    LOGICAL :: islda       = .FALSE.
    LOGICAL :: isgradient  = .FALSE.
    LOGICAL :: has_finite_size_correction = .FALSE.
    LOGICAL :: finite_size_cell_volume_set = .FALSE.
    LOGICAL :: ismeta      = .FALSE.
    LOGICAL :: ishybrid    = .FALSE.
    LOGICAL :: scan_exx    = .FALSE.
    
    
    
    !
    LOGICAL :: discard_input_dft = .FALSE.
    
    INTEGER :: beeftype = -1
    INTEGER :: beefvdw = 0
    
    
    
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
