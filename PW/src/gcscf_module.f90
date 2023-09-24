!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE gcscf_module
  !--------------------------------------------------------------------------
  !
  ! ... This module controls Grand-Canonical SCF for constant-mu method developed
  ! ... by R. Sundararaman, W. A. Goddard-III, T. A. Arias, J. Chem. Phys. 
  ! ... 146, 114104 (2017).
  ! ...
  ! ... Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
  ! ... description of the system connected to a potentiostat which preserves
  ! ... the Fermi energy of the system as the target Fermi energy (mu).
  ! ...
  ! ... Original version is implemented by S. Nishihara (2016-2017) and merged by
  ! ... S. Hagiwara (AIST).
  ! ...
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ! ...   This module is the facade of GCSCF calculations.
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  !
  USE constants,       ONLY : RYTOEV
  USE control_flags,   ONLY : imix, lscf
  USE ener,            ONLY : egrand
  USE esm,             ONLY : do_comp_esm, esm_bc
  USE exx_base,        ONLY : x_gamma_extrapolation
  USE fcp_module,      ONLY : lfcp
  USE fixed_occ,       ONLY : tfixed_occ
  USE xc_lib,          ONLY : exx_is_active
  USE io_global,       ONLY : stdout
  USE ions_base,       ONLY : nat, ityp, zv
  USE kinds,           ONLY : DP
  USE klist,           ONLY : nks, wk, nelec, tot_charge, &
                            & lgauss, degauss, ltetra, two_fermi_energies
  USE mp,              ONLY : mp_sum
  USE mp_pools,        ONLY : inter_pool_comm
  USE rism_module,     ONLY : lrism
  USE wvfct,           ONLY : nbnd, wg
#if defined (__ENVIRON)
  USE plugin_flags,          ONLY : use_environ
#endif
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... variables for GC-SCF method
  LOGICAL  :: lgcscf           = .FALSE.  ! to calculate GC-SCF method, or not
  LOGICAL  :: gcscf_ignore_mun = .FALSE.  ! ignore -mu * N, or not
  LOGICAL  :: gcscf_environ    = .FALSE.
  REAL(DP) :: gcscf_mu         = 0.0_DP   ! target Fermi energy (in Ry)
  REAL(DP) :: gcscf_eps        = 0.0_DP   ! convergence threshold (in Ry)
  REAL(DP) :: gcscf_gk         = 0.0_DP   ! wavenumber shift for Kerker operator (in 1/bohr)
  REAL(DP) :: gcscf_gh         = 0.0_DP   ! wavenumber shift for Hartree metric (in 1/bohr)
  REAL(DP) :: gcscf_beta       = 0.0_DP   ! mixing rate of Fermi energy
  !
  ! ... public components
  PUBLIC :: lgcscf
  PUBLIC :: gcscf_ignore_mun
  PUBLIC :: gcscf_mu
  PUBLIC :: gcscf_eps
  PUBLIC :: gcscf_gk
  PUBLIC :: gcscf_gh
  PUBLIC :: gcscf_beta
  !
  PUBLIC :: gcscf_check
  PUBLIC :: gcscf_iosys
  PUBLIC :: gcscf_summary
  PUBLIC :: gcscf_set_nelec
  PUBLIC :: gcscf_calc_nelec
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_check()
    !----------------------------------------------------------------------------
    !
    ! ... check conditions,
    ! ... and stop program if incorrect conditions.
    !
    IMPLICIT NONE
    !
    ! ... only ESM
    !
#if defined (__ENVIRON)
    IF (use_environ) THEN
      gcscf_environ = .TRUE.
    ENDIF
#endif
    IF (.NOT. gcscf_environ) THEN
         IF (.NOT. do_comp_esm) THEN
            CALL errore('gcscf_check', 'please set assume_isolated = "esm", for GC-SCF', 1)
         END IF
         !
         ! ... cannot use PBC
         !
         IF (TRIM(esm_bc) == 'pbc') THEN
            CALL errore('gcscf_check', 'please do not set esm_bc = "pbc", for GC-SCF', 1)
         END IF
         !
         ! ... cannot use Vacuum/Slab/Vacuum
         !
         IF (TRIM(esm_bc) == 'bc1' .AND. (.NOT. lrism)) THEN
            CALL errore('gcscf_check', 'cannot use ESM-BC1 without RISM, for GC-SCF', 1)
         END IF
    END IF
    !
    ! ... correct Vexx(G=0) ?
    !
    IF (exx_is_active() .AND. (.NOT. x_gamma_extrapolation)) THEN
       CALL errore('gcscf_check', 'GC-SCF calculation requires Vexx(G=0)', 1)
    END IF
    !
    ! ... cannot use FCP
    !
    IF (lfcp) THEN
       CALL errore('gcscf_check', 'cannot use FCP with GC-SCF', 1)
    END IF
    !
    ! ... only for metallic system
    !
    IF (tfixed_occ .OR. ltetra .OR. (.NOT. lgauss) .OR. (degauss <= 0.0_DP)) THEN
       CALL errore('gcscf_check', 'please set occupations = "smearing", for GC-SCF', 1)
    END IF
    !
    ! ... only one Fermi energy
    !
    IF (two_fermi_energies) THEN
       CALL errore('gcscf_check', 'please do not set tot_magnetization, for GC-SCF', 1)
    END IF
    !
    ! ... only TF or local-TF mixing
    !
    IF (imix /= 1 .AND. imix /= 2) THEN
       CALL errore('gcscf_check', 'please set mixing_mode = "TF" or "local-TF", for GC-SCF', 1)
    END IF
    !
    ! ... cannot non-SCF
    !
    IF (.NOT. lscf) THEN
       CALL infomsg('gcscf_check', 'cannot use calculation=nscf for GC-SCF, lgcscf is ignored')
    END IF
   
    !
  END SUBROUTINE gcscf_check
  !
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_iosys(tgcscf)
    !----------------------------------------------------------------------------
    !
    ! ... set variables from input file
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: tgcscf
    !
    lgcscf = tgcscf
    !
    IF (lgcscf) THEN
       !
       CALL iosys_gcscf()
       !
    END IF
    !
  END SUBROUTINE gcscf_iosys
  !
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_summary()
    !----------------------------------------------------------------------------
    !
    ! ... This routine writes on output the GC-SCF's information.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lgcscf) RETURN
    !
    WRITE(stdout, '(/,5X,">>>>> Grand-Canonical SCF is activated <<<<<")' )
    WRITE(stdout, '(5X,"Initial Total Charge  = ",F12.6," e"      )') tot_charge
    WRITE(stdout, '(5X,"Target Fermi Energy   = ",F12.6," eV"     )') gcscf_mu  * RYTOEV
    WRITE(stdout, '(5X,"Thr. of Fermi Energy  = ",F12.6," eV"     )') gcscf_eps * RYTOEV
    WRITE(stdout, '(5X,"Wave-shift of Kerker  = ",F12.6," bohr^-1")') gcscf_gk
    WRITE(stdout, '(5X,"Wave-shift of Hartree = ",F12.6," bohr^-1")') gcscf_gh
    WRITE(stdout, '(5X,"Mixing Rate of Fermi  = ",F12.6           )') gcscf_beta
    !
    FLUSH(stdout)
    !
  END SUBROUTINE gcscf_summary
  !
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_set_nelec(charge)
    !----------------------------------------------------------------------------
    !
    ! ... set number of electrons,
    ! ... also total charge is evaluated.
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: charge
    !
    IF (.NOT. lgcscf) RETURN
    !
    nelec = charge
    !
    tot_charge = SUM(zv(ityp(1:nat))) - nelec
    !
  END SUBROUTINE gcscf_set_nelec
  !
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_calc_nelec()
    !----------------------------------------------------------------------------
    !
    ! ... calculate number of electrons from weights,
    ! ... also total charge and -mu*N are evaluated.
    !
    IMPLICIT NONE
    !
    INTEGER :: ik
    INTEGER :: ibnd
    !
    IF (.NOT. lgcscf) RETURN
    !
    nelec = 0.0_DP
    !
    DO ik = 1, nks
       !
       DO ibnd = 1, nbnd
          !
          nelec = nelec + wg(ibnd, ik)
          !
       END DO
       !
    END DO
    !
    CALL mp_sum(nelec, inter_pool_comm)
    !
    tot_charge = SUM(zv(ityp(1:nat))) - nelec
    !
    egrand = gcscf_mu * tot_charge
    !
  END SUBROUTINE gcscf_calc_nelec
  !
END MODULE gcscf_module
