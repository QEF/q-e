!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_force(rismt, force, vloc, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate force of 3D-RISM.
  ! ...
  ! ... F_sol = F_ele + F_ion + F_lj
  ! ...
  ! ... (1) contribution from electrons, which is included in Hellmann-Feynman force,
  ! ...     is set to zero now.
  ! ...     F_ele = 0
  ! ...
  ! ... (2) contribution from ions (local part)
  ! ...     F_ion = -Omega \Sum_G -n_sol*(G) d V_loc(G)/d R_i
  ! ...
  ! ... (3) contribution from Lennard-Jones potential
  ! ...     F_lj = -Omega \Sum_r n_sol*(r) d V_lj(r)/d R_i
  ! ...
  ! ... (Force of solute: T.Luchko et al., J. Chem. Theory Comput. 2010, 6, 607-624)
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE gvect,     ONLY : ngl
  USE ions_base, ONLY : nsp, nat
  USE kinds,     ONLY : DP
  USE solute,    ONLY : get_solU_LJ_force
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: force(3, nat)
  REAL(DP),        INTENT(IN)  :: vloc(ngl, nsp)
  INTEGER,         INTENT(OUT) :: ierr
  !
  REAL(DP), ALLOCATABLE :: forceion(:,:)
  REAL(DP), ALLOCATABLE :: forcelj (:,:)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  ALLOCATE(forceion(3, nat))
  ALLOCATE(forcelj (3, nat))
  !
  ! ... ions term (local part)
  forceion = 0.0_DP
  CALL solvation_force_ion(rismt, forceion, vloc, ierr)
  !
  IF (ierr /= IERR_RISM_NULL) THEN
    GOTO 1
  END IF
  !
  ! ... Lennad-Jones term
  forcelj = 0.0_DP
  CALL get_solU_LJ_force(rismt, forcelj, ierr)
  !
  IF (ierr /= IERR_RISM_NULL) THEN
    GOTO 1
  END IF
  !
  ! ... total solvation force
  force = forceion + forcelj
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
  ! ... deallocate memory
1 CONTINUE
  DEALLOCATE(forceion)
  DEALLOCATE(forcelj)
  !
END SUBROUTINE solvation_force
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_force_ion(rismt, force, vloc, ierr)
  !---------------------------------------------------------------------------
  !
  USE cell_base,     ONLY : alat, omega
  USE constants,     ONLY : tpi
  USE control_flags, ONLY : gamma_only
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE gvect,         ONLY : ngl, igtongl
  USE ions_base,     ONLY : nsp, nat, ityp, tau
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE rism,          ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: force(3, nat)
  REAL(DP),        INTENT(IN)  :: vloc(ngl, nsp)
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER                  :: ipol
  INTEGER                  :: ig
  INTEGER                  :: na
  REAL(DP)                 :: arg
  REAL(DP)                 :: rhogr
  REAL(DP)                 :: rhogi
  REAL(DP)                 :: vrho
  REAL(DP)                 :: fact
  REAL(DP)                 :: forc1(3)
  REAL(DP),    ALLOCATABLE :: forcesm(:,:)
  COMPLEX(DP), POINTER     :: rhogt(:)
  LOGICAL                  :: laue
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ng < rismt%gvec%ngm) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... allocate memory
  IF (laue .AND. (nat > 0)) THEN
    ALLOCATE(forcesm(3, nat))
  END IF
  !
  ! ... factor of Gamma-trick
  IF (gamma_only) THEN
    fact = 2.0_DP
  ELSE
    fact = 1.0_DP
  END IF
  !
  ! ... set solvent density
  IF (laue) THEN
    rhogt => rismt%rhog_pbc
  ELSE
    rhogt => rismt%rhog
  END IF
  !
  ! ... calculate forces for each atom
  DO na = 1, nat
    DO ipol = 1, 3
      forc1(ipol) = 0.0_DP
    END DO
    !
    ! ... contribution from G=0 is zero
!$omp parallel do default(shared) private(ig, ipol, arg, rhogr, rhogi, vrho) reduction(+:forc1)
    DO ig = rismt%gvec%gstart, rismt%gvec%ngm
      arg = (rismt%gvec%g(1, ig) * tau(1, na) &
        & +  rismt%gvec%g(2, ig) * tau(2, na) &
        & +  rismt%gvec%g(3, ig) * tau(3, na)) * tpi
      rhogr = -DBLE( rhogt(ig))
      rhogi = -AIMAG(rhogt(ig))
      vrho = vloc(igtongl(ig), ityp(na)) * (SIN(arg) * rhogr + COS(arg) * rhogi)
      !
      DO ipol = 1, 3
        forc1(ipol) = forc1(ipol) + rismt%gvec%g(ipol, ig) * vrho
      END DO
    END DO
!$omp end parallel do
    !
    DO ipol = 1, 3
      force(ipol, na) = fact * forc1(ipol) * omega * tpi / alat
    END DO
  END DO
  !
  CALL mp_sum(force, rismt%mp_site%intra_sitg_comm)
  !
  ! ... clear solvent density
  NULLIFY(rhogt)
  !
  ! ... add long-range from ESM(BC1)
  IF (laue .AND. (nat > 0)) THEN
    forcesm = 0.0_DP
    CALL solvation_esm_force(rismt, 1.0_DP / alat, forcesm, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 1
    END IF
    !
    force = force + forcesm
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
  ! ... deallocate memory
1 CONTINUE
  IF (laue .AND. (nat > 0)) THEN
    DEALLOCATE(forcesm)
  END IF
  !
END SUBROUTINE solvation_force_ion
