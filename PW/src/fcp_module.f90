!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_module
  !--------------------------------------------------------------------------
  !
  ! ... This module controls the Fictitious Charge Particle (FCP) for constant-mu
  ! ... method developed by N. Bonnet, T. Morishita, O. Sugino, and M. Otani,
  ! ... Phys. Rev. Lett. 109, 266101 (2012).
  ! ...
  ! ... Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
  ! ... description of the system connected to a potentiostat which preserves
  ! ... the Fermi energy of the system as the target Fermi energy (mu).
  ! ...
  ! ... Original version is implemented by Minoru Otani (AIST) and Nicephore
  ! ... Bonnet (AIST), and merged by S. Hagiwara (AIST).
  ! ... Newton and BFGS algorithms are implemented by S. Nishihara (2016-2017)
  ! ...
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ! ...   This module is the facade of FCP calculations.
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  !
  USE constants,       ONLY : RYTOEV
  USE control_flags,   ONLY : lbfgs, lmd, tr2
  USE dynamics_module, ONLY : dt
  USE esm,             ONLY : do_comp_esm, esm_bc
  USE ener,            ONLY : ef
  USE exx_base,        ONLY : x_gamma_extrapolation
  USE fcp_dynamics,    ONLY : fcpdyn_final, fcpdyn_update, fcpdyn_set_verlet, &
                            & fcpdyn_set_velocity_verlet, fcpdyn_set_proj_verlet
  USE fcp_relaxation,  ONLY : fcprlx_final, fcprlx_update, &
                            & fcprlx_set_line_min, fcprlx_set_newton
  USE fixed_occ,       ONLY : tfixed_occ
  USE xc_lib,          ONLY : exx_is_active
  USE io_global,       ONLY : stdout
  USE kinds,           ONLY : DP
  USE klist,           ONLY : tot_charge, lgauss, degauss, ltetra, two_fermi_energies
  USE relax,           ONLY : starting_scf_threshold
  USE rism_module,     ONLY : lrism
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... variables for FCP method
  LOGICAL           :: lfcp     = .FALSE.  ! to calculate FCP method, or not
  REAL(DP)          :: fcp_mu   = 0.0_DP   ! target Fermi energy (in Ry)
  REAL(DP)          :: fcp_eps  = 0.0_DP   ! convergence threshold (in Ry)
  REAL(DP)          :: fcp_eps0 = 0.0_DP   ! initial convergence threshold (in Ry)
  CHARACTER(LEN=16) :: fcp_calc = ''       ! type of calculation
                                           ! {bfgs|lm|newton|damp|verlet|velocity-verlet}
  !
  ! ... public components
  PUBLIC :: lfcp
  PUBLIC :: fcp_mu
  PUBLIC :: fcp_eps
  PUBLIC :: fcp_eps0
  PUBLIC :: fcp_calc
  !
  PUBLIC :: fcp_check
  PUBLIC :: fcp_iosys
  PUBLIC :: fcp_summary
  PUBLIC :: fcp_relax
  PUBLIC :: fcp_verlet
  PUBLIC :: fcp_terminate
  PUBLIC :: fcp_new_conv_thr
  PUBLIC :: fcp_is_dynamics
  PUBLIC :: output_fcp
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_check(lneb)
    !----------------------------------------------------------------------------
    !
    ! ... check conditions,
    ! ... and stop program if incorrect conditions.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN), OPTIONAL :: lneb
    LOGICAL                       :: lneb_
    !
    lneb_ = .FALSE.
    IF (PRESENT(lneb)) THEN
       lneb_ = lneb
    END IF
    !
    ! ... only ESM
    !
    IF (.NOT. do_comp_esm) THEN
       CALL errore('fcp_check', 'please set assume_isolated = "esm", for FCP', 1)
    END IF
    !
    ! ... cannot use PBC
    !
    IF (TRIM(esm_bc) == 'pbc') THEN
       CALL errore('fcp_check', 'please do not set esm_bc = "pbc", for FCP', 1)
    END IF
    !
    ! ... cannot use Vacuum/Slab/Vacuum
    !
    IF (TRIM(esm_bc) == 'bc1' .AND. (.NOT. lrism)) THEN
       CALL errore('fcp_check', 'cannot use ESM-BC1 without RISM, for FCP', 1)
    END IF
    !
    ! ... correct Vexx(G=0) ?
    !
    IF (exx_is_active() .AND. (.NOT. x_gamma_extrapolation)) THEN
       CALL errore('fcp_check', 'FCP calculation requires Vexx(G=0)', 1)
    END IF
    !
    ! ... only for metallic system
    !
    IF (tfixed_occ .OR. ltetra .OR. (.NOT. lgauss) .OR. (degauss <= 0.0_DP)) THEN
       CALL errore('fcp_check', 'please set occupations = "smearing", for FCP', 1)
    END IF
    !
    ! ... only one Fermi energy
    !
    IF (two_fermi_energies) THEN
       CALL errore('fcp_check', 'please do not set tot_magnetization, for FCP', 1)
    END IF
    !
    ! ... must be relax or md or NEB
    !
    IF (.NOT. (lbfgs .OR. lmd .OR. lneb_)) THEN
       CALL errore('fcp_check', 'calculation has to be relax or md, for FCP', 1)
    END IF
    !
    ! ... cannot use FCP of PWscf, if NEB
    !
    IF (lneb_ .AND. lfcp) THEN
       CALL errore('fcp_check', 'cannot use FCP of PWscf, if NEB', 1)
    END IF
    !
  END SUBROUTINE fcp_check
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_iosys(tfcp)
    !----------------------------------------------------------------------------
    !
    ! ... set variables from input file
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: tfcp
    !
    lfcp = tfcp
    !
    IF (lfcp) THEN
       !
       CALL iosys_fcp()
       !
    END IF
    !
  END SUBROUTINE fcp_iosys
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_summary()
    !----------------------------------------------------------------------------
    !
    ! ... This routine writes on output the FCP's information.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lfcp) RETURN
    !
    IF (fcp_is_dynamics()) THEN
       !
       WRITE(stdout, '(/,5X,">>>>> FCP Dynamics is activated <<<<<<")' )
       !
    ELSE
       !
       WRITE(stdout, '(/,5X,">>>> FCP Relaxation is activated <<<<<")' )
       !
    END IF
    !
    WRITE(stdout, '(5X,"Initial Total Charge = ",F12.6," e"   )') tot_charge
    WRITE(stdout, '(5X,"Target Fermi Energy  = ",F12.6," Ry"  )') fcp_mu
    WRITE(stdout, '(5X,"                     = ",F12.6," eV"  )') fcp_mu * RYTOEV
    !
    FLUSH(stdout)
    !
  END SUBROUTINE fcp_summary
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_relax(conv)
    !----------------------------------------------------------------------------
    !
    ! ... relaxation of FCP. in the current version,
    ! ... Projected-Verlet, Line-Minimization and MDIIS are implemented.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(INOUT) :: conv
    !
    REAL(DP) :: force
    REAL(DP) :: capacitance
    REAL(DP) :: step_max
    !
    IF (.NOT. lfcp) RETURN
    !
    CALL fcp_check()
    !
    ! ... evaluate maximum step
    !
    !force = fcp_mu - ef
    force = 0.1_DP ! max step is 0.1Ry
    !
    CALL fcp_capacitance(capacitance)
    !
    step_max = ABS(capacitance * force)
    !
    ! ... perform each algorithm
    !
    IF (TRIM(fcp_calc) == 'lm') THEN
       !
       ! ... update nelec by Line-Minimization
       !
       CALL fcprlx_set_line_min(fcp_eps, step_max)
       !
       CALL fcprlx_update(fcp_mu, conv)
       !
    ELSE IF (TRIM(fcp_calc) == 'newton') THEN
       !
       ! ... update nelec by Newton-Raphson
       !
       CALL fcprlx_set_newton(fcp_eps, step_max)
       !
       CALL fcprlx_update(fcp_mu, conv)
       !
       !
    ELSE IF (TRIM(fcp_calc) == 'damp') THEN
       !
       ! ... update nelec by Projected-Verlet
       !
       CALL fcpdyn_set_proj_verlet(fcp_eps, step_max)
       !
       CALL fcpdyn_update(fcp_mu, dt, conv)
       !
    ELSE
       !
       CALL errore('fcp_relax', 'incorrect calculation: ' // TRIM(fcp_calc), 1)
       !
    END IF
    !
  END SUBROUTINE fcp_relax
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_verlet()
    !----------------------------------------------------------------------------
    !
    ! ... dynamics of FCP, using Verlet or Velocity-Verlet algorithm.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lfcp) RETURN
    !
    CALL fcp_check()
    !
    IF (TRIM(fcp_calc) == 'verlet') THEN
       !
       ! ... update nelec by Verlet
       !
       CALL fcpdyn_set_verlet()
       !
       CALL fcpdyn_update(fcp_mu, dt)
       !
    ELSE IF (TRIM(fcp_calc) == 'velocity-verlet') THEN
       !
       ! ... update nelec by Velocity-Verlet
       !
       CALL fcpdyn_set_velocity_verlet()
       !
       CALL fcpdyn_update(fcp_mu, dt)
       !
    ELSE
       !
       CALL errore('fcp_verlet', 'incorrect calculation: ' // TRIM(fcp_calc), 1)
       !
    END IF
    !
  END SUBROUTINE fcp_verlet
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_terminate()
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    IF (.NOT. lfcp) RETURN
    !
    IF (fcp_is_dynamics()) THEN
       !
       ! ... finalize fcp_dynamics
       !
       CALL fcpdyn_final()
       !
    ELSE
       !
       ! ... finalize fcp_relaxation
       !
       CALL fcprlx_final()
       !
    END IF
    !
  END SUBROUTINE fcp_terminate
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_new_conv_thr()
    !----------------------------------------------------------------------------
    !
    ! ... update convergence threshold.
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER :: TR2_EXPON = 0.50_DP
    !
    IF (.NOT. lfcp) RETURN
    !
    IF (fcp_eps > 0.0_DP .AND. fcp_eps0               > 0.0_DP .AND. &
        tr2     > 0.0_DP .AND. starting_scf_threshold > 0.0_DP) THEN
       !
       fcp_eps = fcp_eps0 * ((tr2 / starting_scf_threshold) ** TR2_EXPON)
       !
    ELSE
       !
       fcp_eps = fcp_eps0
       !
    END IF
    !
  END SUBROUTINE fcp_new_conv_thr
  !
  !----------------------------------------------------------------------------
  LOGICAL FUNCTION fcp_is_dynamics()
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    IF (TRIM(fcp_calc) == 'damp'   .OR. &
    &   TRIM(fcp_calc) == 'verlet' .OR. &
    &   TRIM(fcp_calc) == 'velocity-verlet') THEN
       !
       fcp_is_dynamics = .TRUE.
       !
    ELSE
       !
       fcp_is_dynamics = .FALSE.
       !
    END IF
    !
  END FUNCTION fcp_is_dynamics
  !
  !----------------------------------------------------------------------------
  SUBROUTINE output_fcp(tot_charge_, conv)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: tot_charge_
    LOGICAL,  INTENT(IN) :: conv
    !
    REAL(DP) :: force
    !
    IF (.NOT. lfcp) RETURN
    !
    IF (.NOT. conv) THEN
       WRITE(stdout, '(5X,"FCP: Total Charge = ",F12.6,"  -> ",F12.6)') tot_charge_, tot_charge
    ELSE
       WRITE(stdout, '(5X,"FCP: Total Charge = ",F12.6)') tot_charge
    END IF
    !
    force = fcp_mu - ef
    !
    WRITE(stdout, '(5X,"FCP: Fermi Energy = ",F12.6," Ry (",F12.6," eV)")') ef,      ef      * RYTOEV
    WRITE(stdout, '(5X,"FCP: Target Level = ",F12.6," Ry (",F12.6," eV)")') fcp_mu,  fcp_mu  * RYTOEV
    WRITE(stdout, '(5X,"FCP: Force on FCP = ",F12.6," Ry (",F12.6," eV)")') force,   force   * RYTOEV
    WRITE(stdout, '(5X,"FCP: Force Thr.   = ",F12.6," Ry (",F12.6," eV)")') fcp_eps, fcp_eps * RYTOEV
    WRITE(stdout, '(/)')
    !
  END SUBROUTINE output_fcp
  !
END MODULE fcp_module
