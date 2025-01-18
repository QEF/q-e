!
! Copyright (C) 2020-2020 Quantum ESPRESSO group
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!------------------------------------------------------------------------------
PROGRAM postahc
!------------------------------------------------------------------------------
!!
!! This program
!!   - Reads the electron-phonon quantities calculated by ph.x with the
!!     electron_phonon='ahc' option.
!!   - Calculate the phonon-induced electron self-energy in the full matrix
!!     form at a given temperature.
!!
!! Input data (namelist "input") is described in Doc/INPUT_POSTAHC.
!!
!------------------------------------------------------------------------------
  USE kinds,       ONLY : DP
  USE constants,   ONLY : ry_to_kelvin, AMU_RY, RY_TO_CMM1
  USE mp,          ONLY : mp_bcast, mp_sum
  USE mp_world,    ONLY : world_comm
  USE mp_global,   ONLY : mp_startup, mp_global_end
  USE io_global,   ONLY : ionode_id, ionode, stdout
  USE environment, ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: nat_max = 1000
  !! Max number of atoms.
  !
  ! Input variables
  !
  LOGICAL :: skip_upperfan
  !! Skip calculation of upper Fan self-energy
  LOGICAL :: skip_dw
  !! Skip calculation of Debye-Waller self-energy
  INTEGER :: nk
  !! Number of k points
  INTEGER :: nbnd
  !! Number of bands in the NSCF calculation
  INTEGER :: ahc_nbnd
  !! Number of bands for which the electron self-energy is to be computed.
  INTEGER :: ahc_nbndskip
  !! Number of bands to exclude when computing the self-energy. The
  !! self-energy is computed for ibnd from ahc_nbndskip + 1
  !! to ahc_nbndskip + ahc_nbnd.
  INTEGER :: nat
  !! Number of atoms
  INTEGER :: nq
  !! Number of q points
  CHARACTER(LEN=256) :: filename
  !! Name of binary files
  CHARACTER(LEN=256) :: ahc_dir
  !! Directory where the output binary files for AHC e-ph coupling are located
  CHARACTER(LEN=256) :: flvec
  !! output file for normalized phonon displacements generated by matdyn.x
  !! The normalized phonon displacements are the eigenvectors divided by the
  !! square root of the mass then normalized. As such they are not orthogonal.
  REAL(DP) :: eta
  !! Infinitesimal to add in the denominator for self-energy, in Ry
  REAL(DP) :: temp_kelvin
  !! Temperature at which the electron self-energy is calculated, in Kelvins.
  REAL(DP) :: efermi
  !! Fermi level energy, in Ry
  REAL(DP) :: amass_amu(nat_max)
  !! Mass of atoms one for each atom (not type), in atomic mass unit.
  !
  ! Local variables
  !
  LOGICAL :: lgamma
  !! .true. if q == Gamma
  INTEGER :: ios
  !! io status
  INTEGER :: iat
  !! Counter for atoms
  INTEGER :: ibnd
  !! Counter for bands
  INTEGER :: jbnd
  !! Counter for bands
  INTEGER :: ik
  !! Counter for k points
  INTEGER :: iq
  !! Counter for q points
  INTEGER :: imode
  !! Counter for modes
  INTEGER :: iun
  !! Unit for reading file
  INTEGER :: recl
  !! Record length for reading file
  INTEGER :: count
  !! Counter for degeneracy
  INTEGER :: nmodes
  !! Number of modes. 3 * nat
  REAL(DP) :: temperature
  !! temp_kelvin transformed from Kelvin to Ry
  REAL(DP) :: unorm
  !! Norm of u multiplied with amass
  REAL(DP) :: rval
  !! Temporary real variables
  REAL(DP) :: omega_zero_cutoff = 1.d-4
  !! Cutoff of phonon frequency. Modes with omega smaller than
  !! omega_zero_cutoff is neglected.
  REAL(DP) :: e_degen_cutoff = 2.d-5
  !! degeneracy cutoff. Ignore couping between degenerate states with energy
  !! difference less than e_degen_cutoff.
  COMPLEX(DP) :: selfen_avg_temp(5)
  !! Diagonal self-energy averaged over degenerate states
  REAL(DP), ALLOCATABLE :: inv_omega(:)
  !! (nmodes) 1 / omega
  REAL(DP), ALLOCATABLE :: occph(:)
  !! (nmodes) Bose-Einstein occupation of phonon
  REAL(DP), ALLOCATABLE :: wtq(:)
  !! (nq) Weight of q points. Set to 1/nq.
  REAL(DP), ALLOCATABLE :: amass(:)
  !! Mass of atoms in Ry
  REAL(DP), ALLOCATABLE :: xq(:, :)
  !! (3, nq) q point vectors in Cartesian basis
  REAL(DP), ALLOCATABLE :: omega(:, :)
  !! (nmodes, nq) Phonon frequency
  REAL(DP), ALLOCATABLE :: etk(:, :)
  !! (nbnd, nk) Energy at k
  REAL(DP), ALLOCATABLE :: etk_all(:, :)
  !! (ahc_nbnd, nk) Energy at k
  COMPLEX(DP), ALLOCATABLE :: u(:, :, :)
  !! (nmodes, nmodes, nq) Phonon modes
  COMPLEX(DP), ALLOCATABLE :: ahc_dw(:, :, :, :, :)
  !! Debye-Waller matrix element
  COMPLEX(DP), ALLOCATABLE :: selfen_dw(:, :, :)
  !! Debye-Waller self-energy
  COMPLEX(DP), ALLOCATABLE :: selfen_upfan(:, :, :)
  !! Upper Fan self-energy
  COMPLEX(DP), ALLOCATABLE :: selfen_lofan(:, :, :)
  !! Lower Fan self-energy
  COMPLEX(DP), ALLOCATABLE :: selfen_fan(:, :, :)
  !! Fan self-energy (lower + upper)
  COMPLEX(DP), ALLOCATABLE :: selfen_tot(:, :, :)
  !! Total self-energy
  COMPLEX(DP), ALLOCATABLE :: selfen_diag(:, :)
  !! Diagonal self-energy
  COMPLEX(DP), ALLOCATABLE :: selfen_diag_avg(:, :)
  !! Diagonal self-energy averaged over degenerate states
  !
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(len=6), EXTERNAL :: int_to_char
  REAL(DP), EXTERNAL :: wgauss
  !
  ! ---------------------------------------------------------------------------
  !
  NAMELIST / input / skip_upperfan, skip_dw, nk, nbnd, nat, nq, ahc_nbnd, &
      ahc_nbndskip, ahc_dir, flvec, eta, temp_kelvin, efermi, amass_amu
  !
  CALL mp_startup()
  CALL environment_start('POSTAHC')
  !
  ! ---------------------------------------------------------------------------
  !
  ! Reading input
  !
  IF (ionode) CALL input_from_file()
  !
  ! Default values of input arguments
  !
  skip_upperfan = .FALSE.
  skip_dw = .FALSE.
  nk = -1
  nbnd = -1
  nat = -1
  nq = -1
  ahc_nbnd = -1
  ahc_nbndskip = 0
  ahc_dir = ' '
  flvec = ' '
  eta = -9999.d0
  temp_kelvin = -9999.d0
  efermi = -9999.d0
  amass_amu(:) = -9999.d0
  !
  ! Read input file
  !
  IF (ionode) READ (5, input, IOSTAT=ios)
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('postahc','error reading input namelist', ABS(ios))
  !
  ! Broadcast input arguments
  !
  CALL mp_bcast(skip_upperfan, ionode_id, world_comm)
  CALL mp_bcast(skip_dw, ionode_id, world_comm)
  CALL mp_bcast(nk, ionode_id, world_comm)
  CALL mp_bcast(nbnd, ionode_id, world_comm)
  CALL mp_bcast(nat, ionode_id, world_comm)
  CALL mp_bcast(nq, ionode_id, world_comm)
  CALL mp_bcast(ahc_nbnd, ionode_id, world_comm)
  CALL mp_bcast(ahc_nbndskip, ionode_id, world_comm)
  CALL mp_bcast(ahc_dir, ionode_id, world_comm)
  CALL mp_bcast(flvec, ionode_id, world_comm)
  CALL mp_bcast(eta, ionode_id, world_comm)
  CALL mp_bcast(temp_kelvin, ionode_id, world_comm)
  CALL mp_bcast(efermi, ionode_id, world_comm)
  CALL mp_bcast(amass_amu, ionode_id, world_comm)
  !
  ! Check input argument validity
  !
  IF (nk < 0) CALL errore('postahc', 'nk must be specified', 1)
  IF (nbnd < 0) CALL errore('postahc', 'nbnd must be specified', 1)
  IF (nat < 0) CALL errore('postahc', 'nat must be specified', 1)
  IF (nq < 0) CALL errore('postahc', 'nq must be specified', 1)
  IF (ahc_nbnd < 0) CALL errore('postahc', 'ahc_nbnd must be specified', 1)
  IF (ahc_dir == '') CALL errore('postahc', 'ahc_dir must be specified', 1)
  IF (flvec == '') CALL errore('postahc', 'flvec must be specified', 1)
  IF (eta < -9990.d0) CALL errore('postahc', 'eta must be specified', 1)
  IF (efermi < -9990.d0) CALL errore('postahc', 'efermi must be specified', 1)
  IF (temp_kelvin < -9990.d0) CALL errore('postahc', &
      'temp_kelvin must be specified', 1)
  DO iat = 1, nat
    IF (amass_amu(iat) < -9990.d0) CALL errore('postahc', &
      'amass_amu(iat) must be specified for iat = 1 to nat', 1)
  ENDDO
  !
  ! Parallelization not implemented. Other processors goto the end of the code.
  !
  IF (.NOT. ionode) GOTO 1001
  !
  ! Setup local variables, unit converstion
  !
  ahc_dir = trimcheck(ahc_dir)
  nmodes = 3 * nat
  temperature = temp_kelvin / ry_to_kelvin
  !
  ALLOCATE(amass(nmodes))
  ALLOCATE(wtq(nq))
  ALLOCATE(xq(3, nq))
  ALLOCATE(omega(nmodes, nq))
  ALLOCATE(u(nmodes, nmodes, nq))
  ALLOCATE(occph(nmodes))
  ALLOCATE(inv_omega(nmodes))
  ALLOCATE(etk_all(nbnd, nk))
  ALLOCATE(etk(ahc_nbnd, nk))
  ALLOCATE(ahc_dw(ahc_nbnd, ahc_nbnd, nmodes, 3, nk))
  ALLOCATE(selfen_dw(ahc_nbnd, ahc_nbnd, nk))
  ALLOCATE(selfen_upfan(ahc_nbnd, ahc_nbnd, nk))
  ALLOCATE(selfen_lofan(ahc_nbnd, ahc_nbnd, nk))
  ALLOCATE(selfen_fan(ahc_nbnd, ahc_nbnd, nk))
  ALLOCATE(selfen_tot(ahc_nbnd, ahc_nbnd, nk))
  ALLOCATE(selfen_diag(ahc_nbnd, 5))
  ALLOCATE(selfen_diag_avg(ahc_nbnd, 5))
  !
  etk_all = 0.d0
  etk = 0.d0
  selfen_dw = (0.d0, 0.d0)
  selfen_upfan = (0.d0, 0.d0)
  selfen_lofan = (0.d0, 0.d0)
  selfen_fan = (0.d0, 0.d0)
  selfen_tot = (0.d0, 0.d0)
  !
  DO iat = 1, nat
    amass((iat-1) * 3 + 1:(iat-1) * 3 + 3) = amass_amu(iat) * AMU_RY
  ENDDO
  wtq = 1.d0 / REAL(nq, KIND=DP)
  !
  ! ---------------------------------------------------------------------------
  !
  ! Read flvec file
  !
  CALL read_flvec(flvec, nq, nmodes, xq, omega, u)
  !
  omega = omega / RY_TO_CMM1
  !
  DO iq = 1, nq
    DO imode = 1, nmodes
      unorm = SUM(CONJG(u(:, imode, iq)) * u(:, imode, iq) * amass(:))
      u(:, imode, iq) = u(:, imode, iq) / SQRT(unorm)
    ENDDO
  ENDDO
  !
  ! ---------------------------------------------------------------------------
  !
  ! Read ahc_dw which does not depend on iq.
  !
  iun = find_free_unit()
  !
  filename = TRIM(ahc_dir) // 'ahc_dw.bin'
  !
  INQUIRE(IOLENGTH=recl) ahc_dw
  OPEN(iun, FILE=TRIM(filename), STATUS='OLD', FORM='UNFORMATTED', &
    ACCESS='DIRECT', RECL=recl, IOSTAT=ios)
  IF (ios /= 0) CALL errore('postahc', 'Error opening ' // TRIM(filename), 1)
  READ(iun, REC=1) ahc_dw
  CLOSE(iun)
  !
  ! Read ahc_etk_iq
  !
  DO iq = 1, nq
    filename = TRIM(ahc_dir) // 'ahc_etk_iq' // TRIM(int_to_char(iq)) // '.bin'
    !
    INQUIRE(IOLENGTH=recl) etk_all
    OPEN(iun, FILE=TRIM(filename), STATUS='OLD', FORM='UNFORMATTED', &
      ACCESS='DIRECT', RECL=recl, IOSTAT=ios)
    IF (ios /= 0) CALL errore('postahc', 'Error opening ' // TRIM(filename), 1)
    READ(iun, REC=1) etk_all
    CLOSE(iun)
  ENDDO
  !
  ! Skip ahc_nbndskip bands in etk
  !
  etk = etk_all(ahc_nbndskip+1:ahc_nbndskip+ahc_nbnd, :)
  !
  ! ---------------------------------------------------------------------------
  !
  ! Loop over q points and calculate self-energy
  !
  WRITE(stdout, *)
  WRITE(stdout,'(5x,a)') 'Calculating electron self-energy. Loop over q points'
  !
  DO iq = 1, nq
    !
    WRITE(stdout, '(i8)', ADVANCE='no') iq
    IF(MOD(iq, 10) == 0 ) WRITE(stdout,*)
    FLUSH(stdout)
    !
    lgamma = .FALSE.
    IF ( ALL( ABS(xq(:, iq)) < 1.d-4 ) ) lgamma = .TRUE.
    !
    ! Set Bose-Einstein occupation of phonons and set inv_omega
    !
    DO imode = 1, nmodes
      IF (omega(imode, iq) < omega_zero_cutoff) THEN
        inv_omega(imode) = 0.d0
        occph(imode) = 0.d0
      ELSE
        inv_omega(imode) = 1.d0 / omega(imode, iq)
        IF (temperature < 1.d-4) THEN
          occph = 0.d0
        ELSE
          rval = wgauss( omega(imode, iq) / temperature, -99 )
          occph(imode) = 1.d0 / (4.d0 * rval - 2.d0) - 0.5d0
        ENDIF
      ENDIF
    ENDDO
    !
    IF (.NOT. skip_dw) CALL calc_debye_waller(iq, selfen_dw)
    !
    IF (.NOT. skip_upperfan) CALL calc_upper_fan(iq, selfen_upfan)
    !
    CALL calc_lower_fan(iq, selfen_lofan)
    !
  ENDDO ! iq
  !
  ! ---------------------------------------------------------------------------
  !
  selfen_fan = selfen_lofan + selfen_upfan
  selfen_tot = selfen_fan + selfen_dw
  !
  ! Write self-energy to stdout
  !
  WRITE(stdout, *)
  WRITE(stdout, *)
  IF (skip_dw) THEN
    WRITE(stdout, '(5x,a)') 'Skip Debye-Waller: Debye-Waller self-energy &
                            &is set to zero'
  ENDIF
  IF (skip_upperfan) THEN
    WRITE(stdout, '(5x,a)') 'Skip upper Fan: upper Fan self-energy is set &
                            &to zero'
  ENDIF
  !
  WRITE(stdout, *)
  WRITE(stdout, '(5x,a)') 'Real part of diagonal electron self-energy in Ry'
  WRITE(stdout, '(5x,a)') 'Self-energy of degenerate states are averaged.'
  WRITE(stdout, '(5x,a)') 'Total_Fan = Upper_Fan + Lower_Fan'
  WRITE(stdout, '(5x,a)') 'Total = Total_Fan + DW'
  WRITE(stdout, *)
  WRITE(stdout, '(5x,a)') 'Begin postahc output'
  WRITE(stdout, '(5x,a)') '    ik  ibnd       Total          DW   Total_Fan&
                          &   Upper_Fan   Lower_Fan'
  DO ik = 1, nk
    DO ibnd = 1, ahc_nbnd
      selfen_diag(ibnd, 1) = selfen_tot(ibnd, ibnd, ik)
      selfen_diag(ibnd, 2) = selfen_dw(ibnd, ibnd, ik)
      selfen_diag(ibnd, 3) = selfen_fan(ibnd, ibnd, ik)
      selfen_diag(ibnd, 4) = selfen_upfan(ibnd, ibnd, ik)
      selfen_diag(ibnd, 5) = selfen_lofan(ibnd, ibnd, ik)
    ENDDO
    !
    ! Average over degenerate states
    !
    DO ibnd = 1, ahc_nbnd
      !
      selfen_avg_temp = (0.d0, 0.d0)
      count = 0
      !
      DO jbnd = 1, ahc_nbnd
        !
        IF (ABS(etk(ibnd, ik) - etk(jbnd, ik)) < e_degen_cutoff) THEN
          count = count + 1
          selfen_avg_temp = selfen_avg_temp + selfen_diag(jbnd, :)
        ENDIF
        !
      ENDDO
      !
      selfen_diag_avg(ibnd, :) = selfen_avg_temp / REAL(count, DP)
      !
    ENDDO ! ibnd
    !
    ! Write averaged self-energy
    !
    DO ibnd = 1, ahc_nbnd
      WRITE(stdout, '(5x, 2I6, 5F12.7)') ik, ibnd, REAL(selfen_diag_avg(ibnd, :))
    ENDDO
    !
  ENDDO ! ik
  !
  WRITE(stdout, '(5x,a)') 'End postahc output'
  WRITE(stdout, *)
  WRITE(stdout, '(5x,a)') 'Full off-diagonal complex self-energy &
                          &matrix is written in files'
  WRITE(stdout, '(5x,a)') 'selfen_real.dat and selfen_imag.dat. &
                          &These data can differ from'
  WRITE(stdout, '(5x,a)') 'the output above because the self-energy &
                          &of degenerate states are'
  WRITE(stdout, '(5x,a)') 'NOT averaged in the selfen_*.dat output.'
  !
  ! Write self-energy to selfen.dat
  !
  iun = find_free_unit()
  !
  OPEN(iun, FILE='selfen_real.dat', FORM='FORMATTED')
  WRITE(iun, '(a)') '# Real part of diagonal electron self-energy &
                    &Re[sigma(ibnd, jbnd, ik)] in Ry'
  WRITE(iun, '(a)') '#   ik  ibnd  jbnd       Total          DW   Total_Fan&
                    &   Upper_Fan   Lower_Fan'
  DO ik = 1, nk
    DO jbnd = 1, ahc_nbnd
      DO ibnd = 1, ahc_nbnd
        WRITE(iun, '(3I6, 5F12.7)') ik, ibnd, jbnd, &
          REAL(selfen_tot(ibnd, jbnd, ik)), REAL(selfen_dw(ibnd, jbnd, ik)), &
          REAL(selfen_fan(ibnd, jbnd, ik)), REAL(selfen_upfan(ibnd, jbnd, ik)), &
          REAL(selfen_lofan(ibnd, jbnd, ik))
      ENDDO
    ENDDO
  ENDDO
  CLOSE(iun)
  !
  OPEN(iun, FILE='selfen_imag.dat', FORM='FORMATTED')
  WRITE(iun, '(a)') '# Imaginary part of diagonal electron self-energy &
                    &Im[sigma(ibnd, jbnd, ik)] in Ry'
  WRITE(iun, '(a)') '#   ik  ibnd  jbnd       Total          DW   Total_Fan&
                    &   Upper_Fan   Lower_Fan'
  DO ik = 1, nk
    DO jbnd = 1, ahc_nbnd
      DO ibnd = 1, ahc_nbnd
        WRITE(iun, '(3I6, 5F12.7)') ik, ibnd, jbnd, &
          AIMAG(selfen_tot(ibnd, jbnd, ik)), AIMAG(selfen_dw(ibnd, jbnd, ik)), &
          AIMAG(selfen_fan(ibnd, jbnd, ik)), AIMAG(selfen_upfan(ibnd, jbnd, ik)), &
          AIMAG(selfen_lofan(ibnd, jbnd, ik))
      ENDDO
    ENDDO
  ENDDO
  CLOSE(iun)
  !
  ! ---------------------------------------------------------------------------
  !
  DEALLOCATE(amass)
  DEALLOCATE(wtq)
  DEALLOCATE(xq)
  DEALLOCATE(omega)
  DEALLOCATE(u)
  DEALLOCATE(occph)
  DEALLOCATE(inv_omega)
  DEALLOCATE(etk_all)
  DEALLOCATE(etk)
  DEALLOCATE(ahc_dw)
  DEALLOCATE(selfen_dw)
  DEALLOCATE(selfen_upfan)
  DEALLOCATE(selfen_lofan)
  DEALLOCATE(selfen_fan)
  DEALLOCATE(selfen_tot)
  DEALLOCATE(selfen_diag)
  DEALLOCATE(selfen_diag_avg)
  !
  ! ---------------------------------------------------------------------------
  !
  WRITE(stdout, *)
  WRITE(stdout, *)
  CALL print_clock('debye_waller')
  CALL print_clock('lower_fan')
  CALL print_clock('upper_fan')
  !
1001 CONTINUE
  !
  CALL environment_end('POSTAHC')
  CALL mp_global_end()
  !
CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE calc_debye_waller(iq, selfen_dw)
    !----------------------------------------------------------------------------
    !!
    !! Compute Debye-Waller self-energy at iq
    !!
    !! Implements Eq.(8) of PHonon/Doc/dfpt_self_energy.pdf
    !!
    !! Here, the "operator-generalized acoustic sum rule" is used to represent
    !! Debye-Waller self-energy as a simple matrix element.
    !! See Eq.(13) of the following reference:
    !! Jae-Mo Lihm and Cheol-Hwan Park, Phys. Rev. B, 101, 121102(R) (2020).
    !!
    !----------------------------------------------------------------------------
    USE kinds,       ONLY : DP
    !
    INTEGER, INTENT(IN) :: iq
    COMPLEX(DP), INTENT(INOUT) :: selfen_dw(ahc_nbnd, ahc_nbnd, nk)
    !
    INTEGER :: imode, jmode, kmode
    !! Counter for modes
    INTEGER :: idir, jdir
    !! Counter for directions
    COMPLEX(DP), ALLOCATABLE :: selfen_dw_iq(:, :, :)
    !! Debye-Waller self-energy at iq
    COMPLEX(DP), ALLOCATABLE :: coeff_dw(:, :, :)
    !! Coefficients for Debye-Waller
    !
    CALL start_clock('debye_waller')
    !
    ALLOCATE(selfen_dw_iq(ahc_nbnd, ahc_nbnd, nk))
    ALLOCATE(coeff_dw(nmodes, 3, nmodes))
    coeff_dw = (0.d0, 0.d0)
    selfen_dw_iq = (0.d0, 0.d0)
    !
    ! coeff_dw(3*iat + idir, jdir, kmode)
    ! = 0.5 * ( CONJG(u(iat*3+idir, kmode)) * u(iat*3+jdir, kmode) ).real
    !   * (occph(kmode) + 0.5) * inv_omega(kmode)
    !
    DO kmode = 1, nmodes
      DO jdir = 1, 3
        DO idir = 1, 3
          DO iat = 1, nat
            imode = 3 * (iat - 1) + idir
            jmode = 3 * (iat - 1) + jdir
            coeff_dw(imode, jdir, kmode) = coeff_dw(imode, jdir, kmode) &
              + 0.5d0 * CONJG(u(imode, kmode, iq)) * u(jmode, kmode, iq) &
                * (occph(kmode) + 0.5d0) * inv_omega(kmode)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    coeff_dw = REAL(coeff_dw, KIND=DP)
    !
    DO kmode = 1, nmodes
      DO jdir = 1, 3
        DO imode = 1, nmodes
          selfen_dw_iq = selfen_dw_iq + ahc_dw(:, :, imode, jdir, :) &
                                      * coeff_dw(imode, jdir, kmode)
        ENDDO
      ENDDO
    ENDDO
    !
    selfen_dw = selfen_dw + selfen_dw_iq * wtq(iq)
    !
    DEALLOCATE(coeff_dw)
    DEALLOCATE(selfen_dw_iq)
    !
    CALL stop_clock('debye_waller')
    !
  !------------------------------------------------------------------------------
  END SUBROUTINE calc_debye_waller
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  SUBROUTINE calc_upper_fan(iq, selfen_upfan)
    !----------------------------------------------------------------------------
    !!
    !! Compute upper Fan (high-energy band contribution from Sternheimer
    !! equation) self-energy at iq
    !!
    !! Implements Eq.(11) of PHonon/Doc/dfpt_self_energy.pdf
    !!
    !----------------------------------------------------------------------------
    USE kinds,       ONLY : DP
    !
    INTEGER, INTENT(IN) :: iq
    COMPLEX(DP), INTENT(INOUT) :: selfen_upfan(ahc_nbnd, ahc_nbnd, nk)
    !
    CHARACTER(LEN=256) :: filename
    !! Name of ahc_upfan_iq*.bin file
    INTEGER :: iun
    !! Unit for reading file
    INTEGER :: recl
    !! Record length for reading file
    INTEGER :: ibnd, jbnd
    !! Counter for bands
    INTEGER :: imode, jmode, kmode
    !! Counter for modes
    REAL(DP) :: coeff
    !! real coefficient
    COMPLEX(DP), ALLOCATABLE :: selfen_upfan_iq(:, :, :)
    !! Upper Fan self-energy at iq
    COMPLEX(DP), ALLOCATABLE :: ahc_upfan(:, :, :, :, :)
    !! Upper Fan matrix element
    COMPLEX(DP), ALLOCATABLE :: ahc_upfan_mode(:, :, :, :)
    !! Upper Fan matrix element in mode basis
    INTEGER, EXTERNAL :: find_free_unit
    !
    CALL start_clock('upper_fan')
    !
    ALLOCATE(selfen_upfan_iq(ahc_nbnd, ahc_nbnd, nk))
    ALLOCATE(ahc_upfan(ahc_nbnd, ahc_nbnd, nmodes, nmodes, nk))
    ALLOCATE(ahc_upfan_mode(ahc_nbnd, ahc_nbnd, nmodes, nk))
    selfen_upfan_iq = (0.d0, 0.d0)
    ahc_upfan_mode = (0.d0, 0.d0)
    !
    ! Read ahc_upfan_iq*.bin
    !
    filename = TRIM(ahc_dir) // 'ahc_upfan_iq' // TRIM(int_to_char(iq)) // '.bin'
    !
    iun = find_free_unit()
    !
    INQUIRE(IOLENGTH=recl) ahc_upfan
    OPEN(iun, FILE=TRIM(filename), STATUS='OLD', FORM='UNFORMATTED', &
      ACCESS='DIRECT', RECL=recl, IOSTAT=ios)
    IF (ios /= 0) CALL errore('postahc', 'Error reading ' // TRIM(filename), 1)
    !
    READ(iun, REC=1) ahc_upfan
    CLOSE(iun)
    !
    ! rotate ahc_upfan from Cartesian to eigenmode basis
    !
    DO ik = 1, nk
      DO imode = 1, nmodes
        DO kmode = 1, nmodes
          DO jmode = 1, nmodes
            ahc_upfan_mode(:, :, imode, ik) = ahc_upfan_mode(:, :, imode, ik) &
            + ahc_upfan(:, :, jmode, kmode, ik) &
            * CONJG(u(jmode, imode, iq)) * u(kmode, imode, iq)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    ! Compute selfen_upfan_iq
    !
    DO ik = 1, nk
      DO imode = 1, nmodes
        !
        coeff = 0.5d0 * inv_omega(imode) * (occph(imode) + 0.5d0)
        !
        DO jbnd = 1, ahc_nbnd
          DO ibnd = 1, ahc_nbnd
            selfen_upfan_iq(ibnd, jbnd, ik) = selfen_upfan_iq(ibnd, jbnd, ik) &
              + ( ahc_upfan_mode(ibnd, jbnd, imode, ik) &
                + CONJG(ahc_upfan_mode(jbnd, ibnd, imode, ik)) ) &
              * coeff
          ENDDO
        ENDDO
        !
      ENDDO
    ENDDO
    !
    selfen_upfan = selfen_upfan + selfen_upfan_iq * wtq(iq)
    !
    DEALLOCATE(ahc_upfan)
    DEALLOCATE(selfen_upfan_iq)
    !
    CALL stop_clock('upper_fan')
    !
  !------------------------------------------------------------------------------
  END SUBROUTINE calc_upper_fan
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  SUBROUTINE calc_lower_fan(iq, selfen_lofan)
    !----------------------------------------------------------------------------
    !!
    !! Compute lower Fan (low-energy band contribution) self-energy at iq
    !!
    !! Implements Eq.(10) of PHonon/Doc/dfpt_self_energy.pdf
    !!
    !----------------------------------------------------------------------------
    USE kinds,       ONLY : DP
    !
    INTEGER, INTENT(IN) :: iq
    COMPLEX(DP), INTENT(INOUT) :: selfen_lofan(ahc_nbnd, ahc_nbnd, nk)
    !
    CHARACTER(LEN=256) :: filename
    !! Name of ahc_gkk_iq*.bin file
    INTEGER :: iun
    !! Unit for reading file
    INTEGER :: recl
    !! Record length for reading file
    INTEGER :: ibq, ibk, jbk
    !! Counter for bands
    INTEGER :: imode, jmode
    !! Counter for modes
    REAL(DP) :: sign
    !! coefficients
    COMPLEX(DP) :: de1, de2
    !! coefficients
    REAL(DP), ALLOCATABLE :: etq(:, :)
    !! (nbnd, nk) Energy at k+q
    REAL(DP), ALLOCATABLE :: occq(:, :)
    !! (nbnd, nk) Fermi-Dirac occupation of energy at k+q
    COMPLEX(DP), ALLOCATABLE :: selfen_lofan_iq(:, :, :)
    !! Lower Fan self-energy at iq
    COMPLEX(DP), ALLOCATABLE :: ahc_gkk(:, :, :, :)
    !! electron-phonon matrix element
    COMPLEX(DP), ALLOCATABLE :: ahc_gkk_mode(:, :, :, :)
    !! electron-phonon matrix element in mode basis
    COMPLEX(DP), ALLOCATABLE :: coeff(:, :, :, :)
    !! coefficient for lower Fan self-energy
    COMPLEX(DP), ALLOCATABLE :: coeff1(:, :, :, :)
    !! coefficient for lower Fan self-energy
    COMPLEX(DP), ALLOCATABLE :: coeff2(:, :, :, :)
    !! coefficient for lower Fan self-energy
    INTEGER, EXTERNAL :: find_free_unit
    !
    CALL start_clock('lower_fan')
    !
    ALLOCATE(selfen_lofan_iq(ahc_nbnd, ahc_nbnd, nk))
    ALLOCATE(ahc_gkk(nbnd, ahc_nbnd, nmodes, nk))
    ALLOCATE(ahc_gkk_mode(nbnd, ahc_nbnd, nmodes, nk))
    ALLOCATE(coeff(nbnd, ahc_nbnd, nmodes, nk))
    ALLOCATE(coeff1(nbnd, ahc_nbnd, nmodes, nk))
    ALLOCATE(coeff2(nbnd, ahc_nbnd, nmodes, nk))
    ALLOCATE(etq(nbnd, nk))
    ALLOCATE(occq(nbnd, nk))
    !
    selfen_lofan_iq = (0.d0, 0.d0)
    ahc_gkk = (0.d0, 0.d0)
    ahc_gkk_mode = (0.d0, 0.d0)
    coeff = (0.d0, 0.d0)
    coeff1 = (0.d0, 0.d0)
    coeff2 = (0.d0, 0.d0)
    etq = 0.d0
    occq = 0.d0
    !
    ! Read files: ahc_etq, ahc_gkk
    !
    iun = find_free_unit()
    !
    filename = TRIM(ahc_dir) // 'ahc_gkk_iq' // TRIM(int_to_char(iq)) // '.bin'
    !
    INQUIRE(IOLENGTH=recl) ahc_gkk
    OPEN(iun, FILE=TRIM(filename), STATUS='OLD', FORM='UNFORMATTED', &
      ACCESS='DIRECT', RECL=recl, IOSTAT=ios)
    IF (ios /= 0) CALL errore('postahc', 'Error opening ' // TRIM(filename), 1)
    READ(iun, REC=1) ahc_gkk
    CLOSE(iun)
    !
    filename = TRIM(ahc_dir) // 'ahc_etq_iq' // TRIM(int_to_char(iq)) // '.bin'
    !
    INQUIRE(IOLENGTH=recl) etq
    OPEN(iun, FILE=TRIM(filename), STATUS='OLD', FORM='UNFORMATTED', &
      ACCESS='DIRECT', RECL=recl, IOSTAT=ios)
    IF (ios /= 0) CALL errore('postahc', 'Error opening ' // TRIM(filename), 1)
    READ(iun, REC=1) etq
    CLOSE(iun)
    !
    ! Fermi-Dirac occupation at k+q
    !
    DO ik = 1, nk
      DO ibnd = 1, nbnd
        IF (temperature < 1.d-4) THEN
          IF (etq(ibnd, ik) < efermi) THEN
            occq(ibnd, ik) = 1.0d0
          ELSEIF (etq(ibnd, ik) > efermi) THEN
            occq(ibnd, ik) = 0.0d0
          ELSE
            occq(ibnd, ik) = 0.5d0
          ENDIF
        ELSE
          occq(ibnd, ik) = wgauss( (efermi - etq(ibnd, ik)) / temperature, -99 )
        ENDIF
      ENDDO
    ENDDO
    !
    ! rotate ahc_gkk from Cartesian to eigenmode basis
    !
    DO ik = 1, nk
      DO imode = 1, nmodes
        DO jmode = 1, nmodes
          ahc_gkk_mode(:, :, imode, ik) = ahc_gkk_mode(:, :, imode, ik) &
          + ahc_gkk(:, :, jmode, ik) * u(jmode, imode, iq)
        ENDDO
      ENDDO
    ENDDO
    !
    ! Compute coefficients
    !
    ! sign = +1 if etk(ibk, ik) > efermi
    !      = -1 otherwise
    !
    ! coeff1(ibq, ibk, imode, ik)
    ! = ( 1 - occq(ibq, ik) + occph(imode) )
    ! / ( etk(ibk, ik) - etq(ibq, ik) - omega(imode) + 1j * eta * sign )
    !
    ! coeff2(ibq, ibk, imode, ik)
    ! = ( occq(ibq, ik) + occph(imode) )
    ! / ( etk(ibk, ik) - etq(ibq, ik) + omega(imode) + 1j * eta * sign )
    !
    ! coeff = (coeff1 + coeff2) * 0.5 * inv_omega(imode)
    !
    DO ik = 1, nk
      DO imode = 1, nmodes
        DO ibk = 1, ahc_nbnd
          !
          sign = 1.d0
          IF (etk(ibk, ik) < efermi) sign = -1.d0
          !
          DO ibq = 1, nbnd
            de1 = CMPLX(etk(ibk, ik) - etq(ibq, ik) - omega(imode, iq),&
                        eta * sign, KIND=DP)
            coeff1(ibq, ibk, imode, ik) = (1.d0 - occq(ibq, ik) + occph(imode)) &
                                        / de1
            !
            de2 = CMPLX(etk(ibk, ik) - etq(ibq, ik) + omega(imode, iq),&
                        eta * sign, KIND=DP)
            coeff2(ibq, ibk, imode, ik) = ( occq(ibq, ik) + occph(imode) ) &
                                        / de2
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    coeff = coeff1 + coeff2
    !
    DO ik = 1, nk
      DO imode = 1, nmodes
        coeff(:, :, imode, ik) = coeff(:, :, imode, ik) * 0.5d0 * inv_omega(imode)
      ENDDO
    ENDDO
    !
    ! Remove coupling with oneself
    !
    IF (lgamma) THEN
      DO ibk = 1, ahc_nbnd
        coeff(ibk + ahc_nbndskip, ibk, :, :) = (0.d0, 0.d0)
      ENDDO
    ENDIF
    !
    ! Remove coupling between degenerate states
    !
    IF (lgamma) THEN
      DO ik = 1, nk
        DO ibk = 1, ahc_nbnd
          DO ibq = 1, nbnd
            IF ( ABS(etk(ibk, ik) - etq(ibq, ik)) < e_degen_cutoff ) THEN
              coeff(ibq, ibk, :, ik) = (0.d0, 0.d0)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ! Compute selfen_lofan_iq
    !
    DO ik = 1, nk
      DO imode = 1, nmodes
        !
        DO jbk = 1, ahc_nbnd
          DO ibk = 1, ahc_nbnd
            DO ibq = 1, nbnd
              !
              selfen_lofan_iq(ibk, jbk, ik) = selfen_lofan_iq(ibk, jbk, ik) &
                  + CONJG(ahc_gkk_mode(ibq, ibk, imode, ik)) &
                  * ahc_gkk_mode(ibq, jbk, imode, ik) &
                  * ( coeff(ibq, ibk, imode, ik) &
                    + coeff(ibq, jbk, imode, ik) ) * 0.5d0
              !
            ENDDO
          ENDDO
        ENDDO
        !
      ENDDO
    ENDDO
    !
    selfen_lofan = selfen_lofan + selfen_lofan_iq * wtq(iq)
    !
    DEALLOCATE(coeff)
    DEALLOCATE(coeff1)
    DEALLOCATE(coeff2)
    DEALLOCATE(ahc_gkk)
    DEALLOCATE(ahc_gkk_mode)
    DEALLOCATE(selfen_lofan_iq)
    DEALLOCATE(etq)
    DEALLOCATE(occq)
    !
    CALL stop_clock('lower_fan')
    !
  !------------------------------------------------------------------------------
  END SUBROUTINE calc_lower_fan
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  SUBROUTINE read_flvec(flvec, nq, nmodes, xq, omega, u)
    !----------------------------------------------------------------------------
    !!
    !! Read flvec (.modes) file generated by matdyn.x
    !!
    !! Output
    !!   - xq : list of q point vectors in Cartesian coordinate
    !!   - omega : phonon frequency in Ry
    !!   - u : phonon modes, normalized to satisfy u^dagger * M * u = identity
    !!         (Eq.(1) of PHonon/Doc/dfpt_self_energy.pdf)
    !!
    !----------------------------------------------------------------------------
    USE kinds,       ONLY : DP
    !
    CHARACTER(LEN=256), INTENT(IN) :: flvec
    !! Name of the modes file
    INTEGER, INTENT(IN) :: nq
    !! Number of q points
    INTEGER, INTENT(IN) :: nmodes
    !! Number of modes
    REAL(DP), INTENT(INOUT) :: xq(3, nq)
    !! q point vectors
    REAL(DP), INTENT(INOUT) :: omega(nmodes, nq)
    !! Phonon frequency
    COMPLEX(DP), INTENT(INOUT) :: u(nmodes, nmodes, nq)
    !! Phonon modes
    !
    INTEGER :: i
    !! dummy variable for reading flvec
    REAL(DP) :: omega_
    !! dummy variable for reading flvec
    INTEGER :: iq
    !! Counter for q points
    INTEGER :: imode
    !! Counter for modes
    INTEGER :: iat
    !! Counter for atoms
    INTEGER :: nat
    !! number of atoms. nmodes / 3
    INTEGER :: iun
    !! Unit for reading flvec
    INTEGER :: ios
    !! io status
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    nat = nmodes / 3
    !
    iun = find_free_unit()
    OPEN(UNIT=iun, FILE=TRIM(flvec), STATUS='OLD', FORM='FORMATTED', IOSTAT=ios)
    IF (ios /= 0) CALL errore('postahc', &
        'problem reading flvec file ' // TRIM(flvec), 1)
    !
    DO iq = 1, nq
      READ(iun, '(a)')
      READ(iun, '(a)')
      READ(iun, '(1x,6x,3f12.4)') (xq(i, iq), i=1, 3)
      READ(iun, '(a)')
      DO imode = 1, nmodes
        READ(iun, 9010) i, omega_, omega(imode, iq)
        DO iat = 1, nat
          READ(iun, 9020) (u(i, imode, iq), i=(iat-1)*3+1, (iat-1)*3+3)
        ENDDO
      ENDDO
      READ(iun, '(a)')
    ENDDO
  9010 format(5x, 6x, i5, 3x, f15.6, 8x, f15.6, 7x)
  9020 format(1x, 1x, 3(f10.6, 1x, f10.6, 3x), 1x)
    !
    CLOSE(iun, IOSTAT=ios)
    IF (ios /= 0) CALL errore('postahc', &
        'problem closing flvec file ' // TRIM(flvec), 1)
    !
  !------------------------------------------------------------------------------
  END SUBROUTINE read_flvec
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  SUBROUTINE postahc_read_unformatted_file(filename, irec, array)
    !----------------------------------------------------------------------------
    !! Read a unformatted file to an array.
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename
    !! Name of the unformatted file to read
    INTEGER, INTENT(IN) :: irec
    !! Record index to read
    COMPLEX(DP), INTENT(INOUT) :: array(:, :, :, :)
    !! Output. Data read from the file.
    !
    INTEGER :: iun
    !! Unit for reading file
    INTEGER :: recl
    !! Record length for reading file
    INTEGER :: ios
    !! io status
    !
    INQUIRE(IOLENGTH=recl) array
    OPEN(NEWUNIT=iun, FILE=TRIM(filename), STATUS='OLD', FORM='UNFORMATTED', &
      ACCESS='DIRECT', RECL=recl, IOSTAT=ios)
    IF (ios /= 0) CALL errore('postahc', 'Error opening ' // TRIM(filename), 1)
    READ(iun, REC=irec) array
    CLOSE(iun, STATUS='KEEP')
  END SUBROUTINE
  !------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
END PROGRAM postahc
!------------------------------------------------------------------------------
