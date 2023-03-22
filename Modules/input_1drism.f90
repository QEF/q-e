!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE iosys_1drism(laue)
  !-----------------------------------------------------------------------------
  !
  ! ...  Copy data read from input file (in subroutine "read_input_file") and
  ! ...  stored in modules input_parameters into internal modules
  ! ...  Note: this subroutine requires pseudo_dir(io_files), omega(cell_base),
  ! ...        ecutrho(gvect), dual(gvecs).
  !
  USE cell_base,        ONLY : alat, at, omega
  USE constants,        ONLY : pi, tpi, eps8, eps12
  USE gvecs,            ONLY : dual
  USE gvect,            ONLY : ecutrho
  USE io_files,         ONLY : molfile
  USE kinds,            ONLY : DP
  USE molecule_const,   ONLY : BOHRm3_TO_MOLCMm3, BOHRm3_TO_MOLLm1
  USE read_solv_module, ONLY : read_solvents
  USE rism,             ONLY : CLOSURE_HNC, CLOSURE_KH
  USE rism1d_facade,    ONLY : nproc_sub, nproc_switch, starting_corr, niter, epsv, &
                             & bond_width, dielectric, molesize, &
                             & mdiis_size, mdiis_step, rism1t, rism1d_initialize, &
                             & rism1d_activate_right, rism1d_activate_left
  USE solvmol,          ONLY : nsolV_ => nsolV, solVs, get_nsite_in_solVs
  !
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : restart_mode
  !
  ! ... RISM namelist
  !
  USE input_parameters, ONLY : nsolv, closure, starting1d, tempv, rmax1d, &
                               smear1d, rism1d_maxstep, rism1d_conv_thr, &
                               rism1d_bond_width, rism1d_dielectric, rism1d_molesize, &
                               rism1d_nproc, rism1d_nproc_switch, mdiis1d_size, mdiis1d_step, &
                               laue_expand_right, laue_expand_left, laue_both_hands
  !
  ! ... SOLVENTS card
  !
  USE input_parameters, ONLY : solv_label, solv_mfile, solv_dens1, solv_dens2, solvents_unit
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: laue
  !
  INTEGER                        :: ngrid
  INTEGER                        :: nsite
  INTEGER                        :: isolV
  INTEGER                        :: iatom
  REAL(DP)                       :: z0
  REAL(DP)                       :: zright
  REAL(DP)                       :: zleft
  CHARACTER(LEN=10), ALLOCATABLE :: slabel(:)
  REAL(DP),          ALLOCATABLE :: sdens1(:)
  REAL(DP),          ALLOCATABLE :: sdens2(:)
  REAL(DP)                       :: tdens1
  REAL(DP)                       :: tdens2
  REAL(DP)                       :: rsol(3)
  REAL(DP)                       :: qsol
  REAL(DP)                       :: dsol(3)
  REAL(DP)                       :: dtot
  INTEGER                        :: ihand
  INTEGER                        :: nhand
  !
  REAL(DP),          PARAMETER   :: RMAX1D_SCALE    = 1.5_DP
  REAL(DP),          PARAMETER   :: BOND_SCALE      = 4.0_DP
  INTEGER,           PARAMETER   :: MDIIS_SWITCH    = 6
  REAL(DP),          PARAMETER   :: MDIIS_STEP_DEF1 = 0.5_DP
  REAL(DP),          PARAMETER   :: MDIIS_STEP_DEF2 = 0.1_DP
  !
  ! ... allocate memory
  ALLOCATE(slabel(nsolv))
  ALLOCATE(sdens1(nsolv))
  ALLOCATE(sdens2(nsolv))
  !
  ! ... check starting condition.
  IF (TRIM(restart_mode) == 'restart') THEN
    IF (TRIM(starting1d) /= 'file' .AND. TRIM(starting1d) /= 'fix') THEN
      CALL infomsg('input','WARNING: "starting1d" set to '//TRIM(starting1d)//' may spoil restart')
      starting1d = 'file'
    END IF
  END IF
  !
  ! ... check both-hand Laue-RISM w/ DRISM.
  IF (laue .AND. laue_both_hands) THEN
    IF (dielectric > 0.0_DP) THEN
      CALL errore('iosys_1drism', 'Cannot use both-hand Laue-RISM with DRISM', 1)
    END IF
  END IF
  !
  ! ... modify rmax1d
  IF (laue) THEN
    z0 = 0.5_DP * alat * at(3, 3)
    zright =  z0 + MAX(0.0_DP, laue_expand_right)
    zleft  = -z0 - MAX(0.0_DP, laue_expand_left )
    rmax1d = MAX(rmax1d, RMAX1D_SCALE * (zright - zleft))
  ENDIF
  !
  ! ... modify rism1d_bond_width
  IF (laue .AND. rism1d_bond_width <= 0.0_DP) THEN
    rism1d_bond_width = BOND_SCALE / SQRT(ecutrho * 4.0_DP / dual)
  ENDIF
  !
  ! ... evaluate #grid
  ngrid = number_of_grids(rmax1d)
  !
  ! ... set from namelist. these data are already checked.
  nsolV_        = nsolv
  starting_corr = starting1d
  nproc_sub     = rism1d_nproc
  nproc_switch  = rism1d_nproc_switch
  niter         = rism1d_maxstep
  epsv          = rism1d_conv_thr
  bond_width    = rism1d_bond_width
  dielectric    = rism1d_dielectric
  molesize      = rism1d_molesize
  mdiis_size    = mdiis1d_size
  mdiis_step    = mdiis1d_step
  !
  ! ... set from card
  DO isolV = 1, nsolV_
    slabel( isolV) = solv_label(isolV)
    molfile(isolV) = solv_mfile(isolV)
    sdens1( isolV) = solv_dens1(isolV)
    sdens2( isolV) = solv_dens2(isolV)
  END DO
  !
  ! ... read solvents from molecule files
  CALL read_solvents()
  !
  ! ... set variables for solVs
  tdens1 = 0.0_DP
  tdens2 = 0.0_DP
  !
  DO isolV = 1, nsolV_
    !
    ! ... name
    IF (LEN_TRIM(slabel(isolV)) > 0) THEN
      solVs(isolV)%name = slabel(isolv)
    ELSE
      CALL infomsg('iosys_1drism', &
      & 'default molecular name(formula) from MOL file('//TRIM(molfile(isolV))//') is used')
    END IF
    IF (LEN_TRIM(solVs(isolV)%name) <= 0) THEN
      CALL errore('iosys_1drism', 'invalid name', isolV)
    END IF
    !
    ! ... density
    IF (sdens1(isolV) >= 0.0_DP) THEN
      CALL convert_dens(TRIM(ADJUSTL(solvents_unit)), isolV, sdens1(isolV))
      solVs(isolV)%density = sdens1(isolv)
    ELSE
      IF (laue .AND. laue_both_hands) THEN
        CALL infomsg('iosys_1drism', &
        & 'default density(right) from MOL file('//TRIM(molfile(isolV))//') is used')
      ELSE
        CALL infomsg('iosys_1drism', &
        & 'default density from MOL file('//TRIM(molfile(isolV))//') is used')
      END IF
    END IF
    !
    IF (solVs(isolV)%density < 0.0_DP) THEN
      IF (laue .AND. laue_both_hands) THEN
        CALL errore('iosys_1drism', 'invalid density(right)', isolV)
      ELSE
        CALL errore('iosys_1drism', 'invalid density', isolV)
      END IF
    END IF
    !
    ! ... subdensity
    IF (laue .AND. laue_both_hands) THEN
      IF (sdens2(isolV) >= 0.0_DP) THEN
        CALL convert_dens(TRIM(ADJUSTL(solvents_unit)), isolV, sdens2(isolV))
        solVs(isolV)%subdensity = sdens2(isolv)
      ELSE
        CALL infomsg('iosys_1drism', &
        & 'default density(left) from MOL file('//TRIM(molfile(isolV))//') is used')
      END IF
      !
      IF (solVs(isolV)%subdensity < 0.0_DP) THEN
        CALL errore('iosys_1drism', 'invalid density(left)', isolV)
      END IF
      !
    ELSE
      solVs(isolV)%subdensity = solVs(isolV)%density
    END IF
    !
    IF (solVs(isolV)%density <= 0.0_DP .AND. solVs(isolV)%subdensity <= 0.0_DP) THEN
      CALL errore('iosys_1drism', 'solvent density is zero', isolV)
    END IF
    !
    tdens1 = tdens1 + solVs(isolV)%density
    tdens2 = tdens2 + solVs(isolV)%subdensity
    !
    ! ... dipole moment, and rotation of coordinate, if DRISM
    solVs(isolV)%dipole   = 0.0_DP
    solVs(isolV)%is_polar = .FALSE.
    !
    IF (dielectric > 0.0_DP) THEN
      !
      rsol = 0.0_DP
      DO iatom = 1, solVs(isolV)%natom
        rsol = rsol + solVs(isolV)%coord(:, iatom)
      END DO
      !
      rsol = rsol / DBLE(solVs(isolV)%natom)
      !
      qsol = 0.0_DP
      dsol = 0.0_DP
      DO iatom = 1, solVs(isolV)%natom
        qsol = qsol + solVs(isolV)%charge(iatom)
        dsol = dsol + solVs(isolV)%charge(iatom) * (solVs(isolV)%coord(:, iatom) - rsol)
      END DO
      !
      dtot = SQRT(dsol(1) * dsol(1) + dsol(2) * dsol(2) + dsol(3) * dsol(3))
      !
      IF (ABS(qsol) <= eps8 .AND. dtot > eps8) THEN
        DO iatom = 1, solVs(isolV)%natom
          solVs(isolV)%coord(:, iatom) = solVs(isolV)%coord(:, iatom) - rsol
        END DO
        !
        CALL rotate_coord(isolV, dsol)
        !
        dsol = 0.0_DP
        DO iatom = 1, solVs(isolV)%natom
          dsol = dsol + solVs(isolV)%charge(iatom) * solVs(isolV)%coord(:, iatom)
        END DO
        !
        IF (ABS(dsol(1)) > eps8 .OR. ABS(dsol(2)) > eps8 .OR. ABS(dtot - dsol(3)) > eps8) THEN
          CALL errore('iosys_1drism', 'error in rotation of molecule', isolV)
        END IF
        !
        solVs(isolV)%dipole   = dsol(3)
        solVs(isolV)%is_polar = .TRUE.
      END IF
      !
    END IF
    !
  END DO
  !
  IF (tdens1 < eps8 .OR. tdens2 < eps8) THEN
    CALL errore('iosys_1drism', 'all of solvent densities are zero', 1)
  END IF
  !
  ! ... modify rism1d_dielectric
  IF (rism1d_dielectric > 0.0_DP) THEN
    IF(.NOT. ANY(solVs(:)%is_polar)) THEN
      rism1d_dielectric = -1.0_DP
      dielectric = rism1d_dielectric
      CALL infomsg('iosys_1drism', 'there are no polar solvents, DRISM is not used.')
    END IF
  END IF
  !
  ! ... modify mdiis1d_step (this operation must be after read_solvents)
  IF (mdiis1d_step < 0.0_DP) THEN
    nsite = get_nsite_in_solVs()
    IF (nsite <= MDIIS_SWITCH) THEN
      mdiis1d_step = MDIIS_STEP_DEF1
    ELSE
      mdiis1d_step = MDIIS_STEP_DEF2
    END IF
    mdiis_step = mdiis1d_step
  END IF
  !
  ! ... initialize rism1d_facade
  IF (laue .AND. laue_both_hands) THEN
    nhand = 2
  ELSE
    nhand = 1
  END IF
  !
  DO ihand = 1, nhand
    IF (ihand == 1) THEN
      CALL rism1d_activate_right()
    ELSE
      CALL rism1d_activate_left()
    END IF
    !
    IF (TRIM(closure) == 'hnc') THEN
      rism1t%closure = CLOSURE_HNC
    ELSE IF (TRIM(closure) == 'kh') THEN
      rism1t%closure = CLOSURE_KH
    END IF
    rism1t%temp = tempv
    rism1t%tau  = smear1d
  END DO
  !
  CALL rism1d_initialize(ngrid, rmax1d, nhand > 1)
  !
  ! ... deallocate memory
  DEALLOCATE(slabel)
  DEALLOCATE(sdens1)
  DEALLOCATE(sdens2)
  !
CONTAINS
  !
  SUBROUTINE convert_dens(dens_format, isolV, dens)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)    :: dens_format
    INTEGER,          INTENT(IN)    :: isolV
    REAL(DP),         INTENT(INOUT) :: dens
    !
    SELECT CASE (dens_format)
    CASE ('1/cell')
      dens = dens / omega
      !
    CASE ('mol/L')
      dens = dens / BOHRm3_TO_MOLLm1
      !
    CASE ('g/cm^3')
      dens = (dens / solVs(isolV)%mass) / BOHRm3_TO_MOLCMm3
      !
    CASE DEFAULT
      CALL errore('iosys_1drism', 'dens_format='// &
                & TRIM(dens_format)//' not implemented', isolV)
      !
    END SELECT
  END SUBROUTINE convert_dens
  !
  SUBROUTINE rotate_coord(isolV, ax)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: isolV
    REAL(DP), INTENT(IN) :: ax(3)
    !
    INTEGER  :: iatom
    REAL(DP) :: er
    REAL(DP) :: e0(3)
    REAL(DP) :: e1(3), r1, s1
    REAL(DP) :: e2(3), r2, s2
    REAL(DP) :: e3(3), r3, s3
    REAL(DP) :: cos0, sin0, tan0
    REAL(DP) :: x, y, z
    REAL(DP) :: xx, yy, xy
    REAL(DP) :: phi
    !
    e3 = ax
    er = SQRT(MAX(0.0_DP, e3(1) * e3(1) + e3(2) * e3(2) + e3(3) * e3(3)))
    IF (er <= 0.0_DP) CALL errore('iosys_1drism', 'e3 = 0', isolV)
    e3 = e3 / er
    !
    IF (ABS(e3(1)) > eps12 .OR. ABS(e3(2)) > eps12) THEN
      !
      e0(1) = 0.0_DP
      e0(2) = 0.0_DP
      e0(3) = 1.0_DP
      !
      cos0 = e0(1) * e3(1) + e0(2) * e3(2) + e0(3) * e3(3)
      sin0 = -SQRT(MAX(0.0_DP, 1.0_DP - cos0 * cos0))
      !
      e1(1) = e0(2) * e3(3) - e0(3) * e3(2)
      e1(2) = e0(3) * e3(1) - e0(1) * e3(3)
      e1(3) = e0(1) * e3(2) - e0(2) * e3(1)
      er = SQRT(MAX(0.0_DP, e1(1) * e1(1) + e1(2) * e1(2) + e1(3) * e1(3)))
      IF (er <= 0.0_DP) CALL errore('iosys_1drism', 'e1 = 0', isolV)
      e1 = e1 / er
      !
      e2(1) = e3(2) * e1(3) - e3(3) * e1(2)
      e2(2) = e3(3) * e1(1) - e3(1) * e1(3)
      e2(3) = e3(1) * e1(2) - e3(2) * e1(1)
      er = SQRT(MAX(0.0_DP, e2(1) * e2(1) + e2(2) * e2(2) + e2(3) * e2(3)))
      IF (er <= 0.0_DP) CALL errore('iosys_1drism', 'e2 = 0', isolV)
      e2 = e2 / er
      !
      DO iatom = 1, solVs(isolV)%natom
        !
        x = solVs(isolV)%coord(1, iatom)
        y = solVs(isolV)%coord(2, iatom)
        z = solVs(isolV)%coord(3, iatom)
        !
        r1 = x * e1(1) + y * e1(2) + z * e1(3)
        r2 = x * e2(1) + y * e2(2) + z * e2(3)
        r3 = x * e3(1) + y * e3(2) + z * e3(3)
        !
        s1 = r1
        s2 = r2 * cos0 - r3 * sin0
        s3 = r2 * sin0 + r3 * cos0
        !
        x = s1 * e1(1) + s2 * e2(1) + s3 * e3(1)
        y = s1 * e1(2) + s2 * e2(2) + s3 * e3(2)
        z = s1 * e1(3) + s2 * e2(3) + s3 * e3(3)
        !
        solVs(isolV)%coord(1, iatom) = x
        solVs(isolV)%coord(2, iatom) = y
        solVs(isolV)%coord(3, iatom) = z
        !
      END DO
      !
    ELSE IF (e3(3) < 0.0_DP) THEN
      !
      DO iatom = 1, solVs(isolV)%natom
        !
        x = solVs(isolV)%coord(1, iatom)
        y = solVs(isolV)%coord(2, iatom)
        z = solVs(isolV)%coord(3, iatom)
        !
        solVs(isolV)%coord(1, iatom) = +x ! rotate with x-axis
        solVs(isolV)%coord(2, iatom) = -y
        solVs(isolV)%coord(3, iatom) = -z
        !
      END DO
      !
    END IF
    !
    xx = 0.0_DP
    yy = 0.0_DP
    xy = 0.0_DP
    !
    DO iatom = 1, solVs(isolV)%natom
      !
      x = solVs(isolV)%coord(1, iatom)
      y = solVs(isolV)%coord(2, iatom)
      !
      xx = xx + x * x
      yy = yy + y * y
      xy = xy + x * y
      !
    END DO
    !
    IF (ABS(xx - yy) > eps12) THEN
      tan0 = -2.0_DP * xy / (xx - yy)
      phi  = ATAN(tan0) / 2.0_DP
    ELSE
      phi  = pi / 4.0_DP
    END IF
    !
    cos0 = COS(phi)
    sin0 = SIN(phi)
    !
    DO iatom = 1, solVs(isolV)%natom
      !
      x = solVs(isolV)%coord(1, iatom)
      y = solVs(isolV)%coord(2, iatom)
      !
      solVs(isolV)%coord(1, iatom) = x * cos0 - y * sin0
      solVs(isolV)%coord(2, iatom) = x * sin0 + y * cos0
      !
    END DO
    !
  END SUBROUTINE rotate_coord
  !
  FUNCTION number_of_grids(rmax) RESULT(nr)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: rmax
    INTEGER              :: nr
    !
    REAL(DP) :: tpibr
    REAL(DP) :: tpibr2
    REAL(DP) :: ecutrism
    REAL(DP) :: gcutrism
    !
    REAL(DP), PARAMETER :: ECUT_SCALE = 1.1_DP
    !
    tpibr  = tpi / (2.0_DP * rmax)
    tpibr2 = tpibr * tpibr
    ecutrism = ECUT_SCALE * MAX(ecutrho, ecutrho * 4.0_DP / dual)
    !
    gcutrism = ecutrism / tpibr2
    nr = INT(SQRT(gcutrism)) + 1
  END FUNCTION number_of_grids
  !
END SUBROUTINE iosys_1drism
