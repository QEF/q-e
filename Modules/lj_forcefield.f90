!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE lj_uff(atmn, eps, sig, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... get Lennard-Jones parameters of UFF.
  ! ... (A.K.Rappe et al., J. Am. Chem. Soc. 1992, 144, 10024-10035)
  !
  ! ... Variables:
  ! ...   atmn: atomic number
  ! ...   eps:  L.J. parameter epsilon (kcal/mol)
  ! ...   sig:  L.J. parameter sigma   (angstrom)
  !
  USE err_rism, ONLY : IERR_RISM_NULL, IERR_RISM_LJ_UNSUPPORTED
  USE kinds,    ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: atmn
  REAL(DP), INTENT(OUT) :: eps
  REAL(DP), INTENT(OUT) :: sig
  INTEGER,  INTENT(OUT) :: ierr
  !
  INTEGER,  PARAMETER :: MAX_ATMN = 103
  REAL(DP)            :: eps_uff(MAX_ATMN)
  REAL(DP)            :: r0_uff( MAX_ATMN)
  !
  DATA eps_uff / &
& 0.044_DP, 0.056_DP, 0.025_DP, 0.085_DP, 0.180_DP, 0.105_DP, 0.069_DP, 0.060_DP, &
& 0.050_DP, 0.042_DP, 0.030_DP, 0.111_DP, 0.505_DP, 0.402_DP, 0.305_DP, 0.274_DP, &
& 0.227_DP, 0.185_DP, 0.035_DP, 0.238_DP, 0.019_DP, 0.017_DP, 0.016_DP, 0.015_DP, &
& 0.013_DP, 0.013_DP, 0.014_DP, 0.015_DP, 0.005_DP, 0.124_DP, 0.415_DP, 0.379_DP, &
& 0.309_DP, 0.291_DP, 0.251_DP, 0.220_DP, 0.040_DP, 0.235_DP, 0.072_DP, 0.069_DP, &
& 0.059_DP, 0.056_DP, 0.048_DP, 0.056_DP, 0.053_DP, 0.048_DP, 0.036_DP, 0.228_DP, &
& 0.599_DP, 0.567_DP, 0.449_DP, 0.398_DP, 0.339_DP, 0.332_DP, 0.045_DP, 0.364_DP, &
& 0.017_DP, 0.013_DP, 0.010_DP, 0.010_DP, 0.009_DP, 0.008_DP, 0.008_DP, 0.009_DP, &
& 0.007_DP, 0.007_DP, 0.007_DP, 0.007_DP, 0.006_DP, 0.228_DP, 0.041_DP, 0.072_DP, &
& 0.081_DP, 0.067_DP, 0.066_DP, 0.037_DP, 0.073_DP, 0.080_DP, 0.039_DP, 0.385_DP, &
& 0.680_DP, 0.663_DP, 0.518_DP, 0.325_DP, 0.284_DP, 0.248_DP, 0.050_DP, 0.404_DP, &
& 0.033_DP, 0.026_DP, 0.022_DP, 0.022_DP, 0.019_DP, 0.016_DP, 0.014_DP, 0.013_DP, &
& 0.013_DP, 0.013_DP, 0.012_DP, 0.012_DP, 0.011_DP, 0.011_DP, 0.011_DP /
  !
  DATA r0_uff / &
& 2.886_DP, 2.362_DP, 2.451_DP, 2.745_DP, 4.083_DP, 3.851_DP, 3.660_DP, 3.500_DP, &
& 3.364_DP, 3.243_DP, 2.983_DP, 3.021_DP, 4.499_DP, 4.295_DP, 4.147_DP, 4.035_DP, &
& 3.947_DP, 3.868_DP, 3.812_DP, 3.399_DP, 3.295_DP, 3.175_DP, 3.144_DP, 3.023_DP, &
& 2.961_DP, 2.912_DP, 2.872_DP, 2.834_DP, 3.495_DP, 2.763_DP, 4.383_DP, 4.280_DP, &
& 4.230_DP, 4.205_DP, 4.189_DP, 4.141_DP, 4.114_DP, 3.641_DP, 3.345_DP, 3.124_DP, &
& 3.165_DP, 3.052_DP, 2.998_DP, 2.963_DP, 2.929_DP, 2.899_DP, 3.148_DP, 2.848_DP, &
& 4.463_DP, 4.392_DP, 4.420_DP, 4.470_DP, 4.500_DP, 4.404_DP, 4.517_DP, 3.703_DP, &
& 3.522_DP, 3.556_DP, 3.606_DP, 3.575_DP, 3.547_DP, 3.520_DP, 3.493_DP, 3.368_DP, &
& 3.451_DP, 3.428_DP, 3.409_DP, 3.391_DP, 3.374_DP, 3.355_DP, 3.640_DP, 3.141_DP, &
& 3.170_DP, 3.069_DP, 2.954_DP, 3.120_DP, 2.840_DP, 2.754_DP, 3.293_DP, 2.705_DP, &
& 4.347_DP, 4.297_DP, 4.370_DP, 4.709_DP, 4.750_DP, 4.765_DP, 4.900_DP, 3.677_DP, &
& 3.478_DP, 3.396_DP, 3.424_DP, 3.395_DP, 3.424_DP, 3.424_DP, 3.381_DP, 3.326_DP, &
& 3.339_DP, 3.313_DP, 3.299_DP, 3.286_DP, 3.274_DP, 3.248_DP, 3.236_DP /
  !
  IF (1 <= atmn .AND. atmn <= MAX_ATMN) THEN
    eps  = eps_uff(atmn)
    sig  = r0_uff(atmn) / (2.0_DP ** (1.0_DP / 6.0_DP))
    ierr = IERR_RISM_NULL
    !
  ELSE
    eps  = 0.0_DP
    sig  = 0.0_DP
    ierr = IERR_RISM_LJ_UNSUPPORTED
    !
  END IF
  !
END SUBROUTINE lj_uff
!
!--------------------------------------------------------------------------
SUBROUTINE lj_clayff(atmn, crdn, eps, sig, label, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... get Lennard-Jones parameters of ClayFF.
  ! ... (R.T.Cygan et al., J. Phys. Chem. B 2004, 108, 1255-1266)
  !
  ! ... Variables:
  ! ...   atmn:  atomic number
  ! ...   crdn:  coordination number
  ! ...   eps:   L.J. parameter epsilon (kcal/mol)
  ! ...   sig:   L.J. parameter sigma   (angstrom)
  ! ...   label: optional name of element
  !
  USE err_rism, ONLY : IERR_RISM_NULL, IERR_RISM_LJ_UNSUPPORTED
  USE kinds,    ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,          INTENT(IN)  :: atmn
  INTEGER,          INTENT(IN)  :: crdn
  REAL(DP),         INTENT(OUT) :: eps
  REAL(DP),         INTENT(OUT) :: sig
  CHARACTER(LEN=5), INTENT(OUT) :: label
  INTEGER,          INTENT(OUT) :: ierr
  !
  REAL(DP) :: r0
  !
  eps   = 0.0_DP
  r0    = 0.0_DP
  label = ''
  ierr  = IERR_RISM_LJ_UNSUPPORTED
  !
  SELECT CASE(atmn)
  CASE( 1) ! H
    eps  = 0.0460_DP
    r0   = 1.1225_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 8) ! O
    eps  = 0.1554_DP
    r0   = 3.5532_DP
    ierr = IERR_RISM_NULL
    !
  CASE(14) ! Si
    IF (3 <= crdn .AND. crdn <= 4) THEN ! tetrahedral
      eps   = 1.8405E-6_DP
      r0    = 3.7064_DP
      label = ' [Td]'
      ierr  = IERR_RISM_NULL
    END IF
    !
  CASE(13) ! Al
    IF (3 <= crdn .AND. crdn <= 4) THEN ! tetrahedral
      eps   = 1.8405E-6_DP
      r0    = 4.7943_DP
      label = ' [Td]'
      ierr  = IERR_RISM_NULL
      !
    ELSE IF (5 <= crdn .AND. crdn <= 6) THEN ! octahedral
      eps   = 1.3298E-6_DP
      r0    = 4.7943_DP
      label = ' [Oh]'
      ierr  = IERR_RISM_NULL
      !
    END IF
    !
  CASE(12) ! Mg
    IF (5 <= crdn .AND. crdn <= 6) THEN ! octahedral
      eps   = 9.0298E-7_DP
      r0    = 5.9090_DP
      label = ' [Oh]'
      ierr  = IERR_RISM_NULL
    END IF
    !
  CASE(20) ! Ca
    IF (5 <= crdn .AND. crdn <= 6) THEN ! octahedral
      eps   = 5.0298E-6_DP
      r0    = 6.2484_DP
      label = ' [Oh]'
      ierr  = IERR_RISM_NULL
    END IF
    !
  CASE(26) ! Fe
    IF (5 <= crdn .AND. crdn <= 6) THEN ! octahedral
      eps   = 9.0298E-6_DP
      r0    = 5.5070_DP
      label = ' [Oh]'
      ierr  = IERR_RISM_NULL
    END IF
    !
  CASE( 3) ! Li
    IF (5 <= crdn .AND. crdn <= 6) THEN ! octahedral
      eps   = 9.0298E-6_DP
      r0    = 4.7257_DP
      label = ' [Oh]'
      ierr  = IERR_RISM_NULL
    END IF
    !
  END SELECT
  !
  sig = r0 / (2.0_DP ** (1.0_DP / 6.0_DP))
  !
END SUBROUTINE lj_clayff
!
!--------------------------------------------------------------------------
SUBROUTINE lj_oplsaa(atmn, eps, sig, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... get Lennard-Jones parameters of OPLS-AA.
  ! ... (generic parameters for QM/MM)
  !
  ! ... Variables:
  ! ...   atmn: atomic number
  ! ...   eps:  L.J. parameter epsilon (kcal/mol)
  ! ...   sig:  L.J. parameter sigma   (angstrom)
  !
  USE err_rism, ONLY : IERR_RISM_NULL, IERR_RISM_LJ_UNSUPPORTED
  USE kinds,    ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: atmn
  REAL(DP), INTENT(OUT) :: eps
  REAL(DP), INTENT(OUT) :: sig
  INTEGER,  INTENT(OUT) :: ierr
  !
  eps  = 0.0_DP
  sig  = 0.0_DP
  ierr = IERR_RISM_LJ_UNSUPPORTED
  !
  SELECT CASE(atmn)
  CASE( 1) ! H
    eps  = 0.030_DP
    sig  = 2.460_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 2) ! He
    eps  = 0.020_DP
    sig  = 2.556_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 3) ! Li
    eps  = 0.018_DP
    sig  = 2.126_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 4) ! Be
    eps  = 0.050_DP
    sig  = 3.250_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 5) ! B
    eps  = 0.050_DP
    sig  = 3.600_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 6) ! C
    eps  = 0.070_DP
    sig  = 3.550_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 7) ! N
    eps  = 0.170_DP
    sig  = 3.250_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 8) ! O
    eps  = 0.170_DP
    sig  = 3.000_DP
    ierr = IERR_RISM_NULL
    !
  CASE( 9) ! F
    eps  = 0.060_DP
    sig  = 2.900_DP
    ierr = IERR_RISM_NULL
    !
  CASE(10) ! Ne
    eps  = 0.069_DP
    sig  = 2.780_DP
    ierr = IERR_RISM_NULL
    !
  CASE(11) ! Na
    eps  = 0.003_DP
    sig  = 3.330_DP
    ierr = IERR_RISM_NULL
    !
  CASE(12) ! Mg
    eps  = 0.050_DP
    sig  = 3.400_DP
    ierr = IERR_RISM_NULL
    !
  CASE(13) ! Al
    eps  = 0.100_DP
    sig  = 4.050_DP
    ierr = IERR_RISM_NULL
    !
  CASE(14) ! Si
    eps  = 0.100_DP
    sig  = 4.000_DP
    ierr = IERR_RISM_NULL
    !
  CASE(15) ! P
    eps  = 0.200_DP
    sig  = 3.740_DP
    ierr = IERR_RISM_NULL
    !
  CASE(16) ! S
    eps  = 0.250_DP
    sig  = 3.550_DP
    ierr = IERR_RISM_NULL
    !
  CASE(17) ! Cl
    eps  = 0.300_DP
    sig  = 3.400_DP
    ierr = IERR_RISM_NULL
    !
  CASE(18) ! Ar
    eps  = 0.234_DP
    sig  = 3.401_DP
    ierr = IERR_RISM_NULL
    !
  CASE(35) ! Br
    eps  = 0.470_DP
    sig  = 3.470_DP
    ierr = IERR_RISM_NULL
    !
  CASE(53) ! I
    eps  = 0.580_DP
    sig  = 3.550_DP
    ierr = IERR_RISM_NULL
    !
  END SELECT
  !
END SUBROUTINE lj_oplsaa
