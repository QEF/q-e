!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE iosys_fcp()
  !--------------------------------------------------------------------------
  !
  ! ...  Copy data read from input file (in subroutine "read_input_file") and
  ! ...  stored in modules input_parameters into internal modules of FCP
  !
  USE cell_base,             ONLY : alat, at
  USE constants,             ONLY : RYTOEV
  USE control_flags,         ONLY : lbfgs, lmd
  USE fcp_dynamics,          ONLY : fcpdyn_init, fcpdyn_prm_mass, &
                                  & fcpdyn_prm_velocity, fcpdyn_prm_temp
  USE fcp_module,            ONLY : fcp_mu_ => fcp_mu, &
                                  & fcp_eps, fcp_eps0, fcp_calc, fcp_check, &
                                  & fcp_is_dynamics
  USE fcp_relaxation,        ONLY : fcprlx_init, fcprlx_prm
  USE ions_base,             ONLY : if_pos
  USE kinds,                 ONLY : DP
  USE read_namelists_module, ONLY : fcp_not_set
  USE rism_module,           ONLY : lrism
  !
  ! ... CONTROL namelist
  !
  USE input_parameters,      ONLY : calculation
  !
  ! ... FCP namelist
  !
  USE input_parameters,      ONLY : fcp_mu, fcp_dynamics_ => fcp_dynamics, fcp_conv_thr, &
                                  & fcp_ndiis, fcp_rdiis, fcp_mass, fcp_velocity, fcp_temperature, &
                                  & fcp_tempw, fcp_tolp, fcp_delta_t, fcp_nraise, &
                                  & freeze_all_atoms
  !
  IMPLICIT NONE
  !
  REAL(DP) :: area_xy
  !
  REAL(DP), PARAMETER :: MASS_DEF   = 5.0E+6_DP
  REAL(DP), PARAMETER :: SCALE_RISM = 100.0_DP
  !
  ! ... modify fcp_mass
  !
  IF (fcp_mass <= 0.0_DP) THEN
     !
     area_xy = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
     !
     fcp_mass = MASS_DEF / area_xy
     !
     IF (lrism) THEN
        !
        fcp_mass = fcp_mass / SCALE_RISM
        !
     END IF
     !
  END IF
  !
  ! ... resolve fcp_dynamics -> fcp_calc
  !
  SELECT CASE (TRIM(calculation))
     !
  CASE ('relax', 'vc-relax')
     !
     SELECT CASE(TRIM(fcp_dynamics_))
        !
     CASE('lm', 'line-min', 'line-minimization', 'line-minimisation')
        !
        fcp_calc = 'lm'
        !
     CASE('damp')
        !
        fcp_calc = 'damp'
        !
     CASE('newton')
        !
        fcp_calc = 'newton'
        !
     CASE('bfgs')
        !
        fcp_calc = 'bfgs'
        !
     CASE DEFAULT
        !
        CALL errore('iosys', 'calculation=' // TRIM(calculation) // &
                  & ': fcp_dynamics=' // TRIM(fcp_dynamics_) // &
                  & ' not supported', 1)
        !
     END SELECT
     !
     IF (lbfgs .AND. TRIM(fcp_calc) /= 'bfgs') THEN
        !
        fcp_calc = 'bfgs'
        !
        CALL infomsg('iosys', 'calculation='// TRIM(calculation) // &
                   & ': fcp_dynamics=' // TRIM(fcp_dynamics_) // &
                   & " ignored, 'bfgs' assumed")
        !
     END IF
     !
     IF (lmd .AND. TRIM(fcp_calc) == 'bfgs') THEN
        !
        fcp_calc = 'lm'
        !
        CALL infomsg('iosys', 'calculation='// TRIM(calculation) // &
                   & ': fcp_dynamics=' // TRIM(fcp_dynamics_) // &
                   & " ignored, 'lm' assumed")
        !
     END IF
     !
  CASE ('md')
     !
     SELECT CASE(TRIM(fcp_dynamics_))
        !
     CASE('verlet')
        !
        fcp_calc = 'verlet'
        !
     CASE('vv', 'vverlet', 'velocityverlet', 'velocity-verlet')
        !
        fcp_calc = 'velocity-verlet'
        !
     CASE DEFAULT
        !
        CALL errore('iosys', 'calculation=' // TRIM(calculation) // &
                  & ': fcp_dynamics=' // TRIM(fcp_dynamics_) // &
                  & ' not supported', 1 )
        !
     END SELECT
     !
  CASE DEFAULT
     !
     CALL errore('iosys', 'calculation=' // TRIM(calculation) // &
               & ' not supported, for FCP', 1)
     !
  END SELECT
  !
  ! ... set variables from namelist
  !
  fcp_mu_  = fcp_mu / RYTOEV
  fcp_eps  = fcp_conv_thr / RYTOEV
  fcp_eps0 = fcp_eps
  !
  IF (fcp_is_dynamics()) THEN
     !
     ! ... initialize fcp_dynamics
     !
     CALL fcpdyn_init()
     !
     ! ... set parameters of fcp_dynamics
     !
     CALL fcpdyn_prm_mass(fcp_mass)
     !
     IF (fcp_velocity /= fcp_not_set) THEN
        !
        CALL fcpdyn_prm_velocity(fcp_velocity)
        !
     END IF
     !
     CALL fcpdyn_prm_temp(fcp_temperature, &
                        & fcp_tempw, fcp_tolp, fcp_delta_t, fcp_nraise)
     !
  ELSE
     !
     ! ... initialize fcp_relaxation
     !
     CALL fcprlx_init()
     !
     ! ... set parameters of fcp_relaxation
     !
     CALL fcprlx_prm(fcp_ndiis, fcp_rdiis)
     !
  END IF
  !
  ! ... freeze all atoms
  !
  IF (freeze_all_atoms) THEN
     !
     if_pos(:, :) = 0
     !
  END IF
  !
  ! ... check condition
  !
  CALL fcp_check()
  !
END SUBROUTINE iosys_fcp
