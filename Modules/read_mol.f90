!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE read_mol_module
  !--------------------------------------------------------------------------
  !
  ! ... this module handles the reading of molecular data
  !
  USE constants,      ONLY : BOHR_RADIUS_ANGS, RYTOEV
  USE kinds,          ONLY : DP
  USE molecule_const, ONLY : RY_TO_KJMOLm1, RY_TO_KCALMOLm1, BOHRm3_TO_MOLCMm3, BOHRm3_TO_MOLLm1
  USE molecule_types, ONLY : molecule, deallocate_molecule
  USE parser,         ONLY : version_compare
  USE iotk_module
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... data to avoid reading
  LOGICAL :: without_mass         = .FALSE.
  LOGICAL :: without_density      = .FALSE.
  LOGICAL :: without_permittivity = .FALSE.
  LOGICAL :: without_element      = .FALSE.
  LOGICAL :: without_xyz          = .FALSE.
  LOGICAL :: without_charge       = .FALSE.
  LOGICAL :: without_lj           = .FALSE.
  !
  ! ... public components
  PUBLIC :: without_mass
  PUBLIC :: without_density
  PUBLIC :: without_permittivity
  PUBLIC :: without_element
  PUBLIC :: without_xyz
  PUBLIC :: without_charge
  PUBLIC :: without_lj
  PUBLIC :: clean_mol_readables
  PUBLIC :: read_mol
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE clean_mol_readables()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    without_mass         = .FALSE.
    without_density      = .FALSE.
    without_permittivity = .FALSE.
    without_element      = .FALSE.
    without_xyz          = .FALSE.
    without_charge       = .FALSE.
    without_lj           = .FALSE.
  END SUBROUTINE clean_mol_readables
  !
  !--------------------------------------------------------------------------
  SUBROUTINE read_mol(mol, ierr, unit, filename)
    !--------------------------------------------------------------------------
    !
    ! ... read molecule in MOL format
    ! ... ierr =  0 : read MOL v.1
    ! ... ierr =  1 : not an MOL file, or error while reading
    !
    IMPLICIT NONE
    !
    TYPE(molecule),             INTENT(INOUT) :: mol       ! the molecular data
    INTEGER,                    INTENT(OUT)   :: ierr
    INTEGER,          OPTIONAL, INTENT(IN)    :: unit      ! i/o unit
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: filename  ! i/o filename
    !
    INTEGER :: u  ! i/o unit
    !
    ierr = 0
    !
    IF (.NOT. PRESENT(unit)) THEN
      IF (.NOT. PRESENT(filename)) THEN
        CALL errore('read_mol', 'You have to specify at least one between filename and unit', 1)
      END IF
      CALL iotk_free_unit(u)
    ELSE
      u = unit
    END IF
    !
    IF (PRESENT(filename)) THEN
      OPEN(unit=u, file=filename, status='old', form='formatted', iostat=ierr)
    END IF
    IF (ierr > 0) THEN
      CALL errore('read_mol', 'Cannot open file: ' // TRIM(filename), 1)
    END IF
    !
    CALL read_mol_v1(u, mol, ierr)
    !
    IF (ierr > 0) THEN
      REWIND(u)
      CALL deallocate_molecule(mol)
    END IF
    !
  END SUBROUTINE read_mol
  !
  !--------------------------------------------------------------------------
  SUBROUTINE read_mol_v1(u, mol, ierr)
    !--------------------------------------------------------------------------
    !
    ! ... read molecule in MOL format(v.1), using iotk
    !
    IMPLICIT NONE
    !
    INTEGER,           INTENT(IN)    :: u    ! i/o unit
    TYPE(molecule),    INTENT(INOUT) :: mol  ! the molecular data
    INTEGER, OPTIONAL, INTENT(OUT)   :: ierr ! /= 0 if something went wrong
    !
    CHARACTER(LEN=iotk_namlenx) :: root
    CHARACTER(LEN=iotk_attlenx) :: attr
    INTEGER                     :: ierr_
    LOGICAL                     :: found
    CHARACTER(len=6), PARAMETER :: max_version = '1.0.0'
    !
    LOGICAL, EXTERNAL :: matches
    !
    ! ... initialize the file
    CALL iotk_open_read(u, attr=attr, root=root, ierr=ierr_)
    !
    IF((ABS(ierr_) > 0) .OR. (.NOT. matches('MOL', root))) THEN
      !
      CALL iotk_close_read(u, ierr=ierr_)
      IF(.NOT. PRESENT(ierr)) THEN
        CALL errore('read_mol_v1', 'Cannot open MOL file.', 1)
      END IF
      ierr = 1
      RETURN
    END IF
    !
    CALL iotk_scan_attr(attr, 'version', mol%nv)
    IF (version_compare(mol%nv, max_version) == 'newer') THEN
      CALL errore('read_mol_v1', 'Unknown MOL format version: ' // TRIM(mol%nv), 1)
    END IF
    !
    ! ... skip human-readable header
    CALL iotk_scan_begin(u, 'MOL_INFO', found=found)
    if(found) CALL iotk_scan_end(u, 'MOL_INFO')
    !
    ! ... read machine-readable header
    CALL read_mol_header(u, mol)
    !
    ! ... read mass of molecule
    IF (.NOT. without_mass) THEN
      CALL read_mol_mass(u, mol)
    END IF
    !
    ! ... read density of molecule
    IF (.NOT. without_density) THEN
      CALL read_mol_density(u, mol)
    END IF
    !
    ! ... read permittivity of molecule
    IF (.NOT. without_permittivity) THEN
      CALL read_mol_permittivity(u, mol)
    END IF
    !
    ! ... read element of every atom
    IF (.NOT. without_element) THEN
      CALL read_mol_element(u, mol)
    END IF
    !
    ! ... read xyz-coordinate of every atom
    IF (.NOT. without_xyz) THEN
       CALL read_mol_xyz(u, mol)
    END IF
    !
    ! ... read charge of every atom
    IF (.NOT. without_charge) THEN
      CALL read_mol_charge(u, mol)
    END IF
    !
    ! ... read Lennard-Jones of every atom
    IF (.NOT. without_lj) THEN
      CALL read_mol_lj(u, mol)
    END IF
    !
    ! ... close the file (not the unit!)
    CALL iotk_close_read(u)
    !
    IF (PRESENT(ierr)) ierr = 0
    !
  CONTAINS
    !
    SUBROUTINE read_mol_header(u, mol)
      IMPLICIT NONE
      INTEGER,        INTENT(IN)    :: u    ! i/o unit
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      CHARACTER(LEN=iotk_attlenx) :: attr
      !
      CALL iotk_scan_empty(u, 'MOL_HEADER', attr=attr)
      !
      CALL iotk_scan_attr(attr, 'author',          mol%author,     default='anonymous')
      CALL iotk_scan_attr(attr, 'date',            mol%date,       default=' ')
      CALL iotk_scan_attr(attr, 'comment',         mol%comment,    default=' ')
      !
      CALL iotk_scan_attr(attr, 'formula',         mol%formula)
      IF (LEN_TRIM(mol%formula) < 1) THEN
        CALL errore('read_mol_v1', 'formula is empty @MOL_HEADER', 1)
      END IF
      mol%name = ADJUSTL(mol%formula)
      !
      CALL iotk_scan_attr(attr, 'number_of_atoms', mol%natom)
      IF (mol%natom < 1) THEN
        CALL errore('read_mol_v1', 'number_of_atoms is not positive @MOL_HEADER', 1)
      END IF
      !
      CALL iotk_scan_attr(attr, 'has_charge',      mol%has_charge, default=.FALSE.)
      CALL iotk_scan_attr(attr, 'has_lj',          mol%has_lj,     default=.FALSE.)
    END SUBROUTINE read_mol_header
    !
    SUBROUTINE read_mol_mass(u, mol)
      IMPLICIT NONE
      INTEGER,        INTENT(IN)    :: u    ! i/o unit
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      CHARACTER(LEN=iotk_attlenx) :: attr
      CHARACTER(LEN=16)           :: units
      INTEGER                     :: i
      !
      CHARACTER(LEN=1), EXTERNAL :: capital
      !
      CALL iotk_scan_dat(u, 'MOL_MASS', mol%mass, attr=attr)
      IF (mol%mass <= 0.0_DP) THEN
        CALL errore('read_mol_v1', 'molecular mass is not positive @MOL_MASS', 1)
      END IF
      !
      CALL iotk_scan_attr(attr, 'UNITS', units, default='a.m.u.')
      DO i = 1, LEN_TRIM(units)
        units(i:i) = capital(units(i:i))
      END DO
      !
      SELECT CASE (TRIM(units))
      CASE ('A.M.U.')
        ! NOP
      CASE ('G/MOL')
        ! NOP
      CASE DEFAULT
        CALL errore('read_mol_v1', 'incorrect units @MOL_MASS', 1)
      END SELECT
    END SUBROUTINE read_mol_mass
    !
    SUBROUTINE read_mol_density(u, mol)
      IMPLICIT NONE
      INTEGER,        INTENT(IN)    :: u    ! i/o unit
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      CHARACTER(LEN=iotk_attlenx) :: attr
      CHARACTER(LEN=16)           :: units
      INTEGER                     :: i
      !
      CHARACTER(LEN=1), EXTERNAL :: capital
      !
      CALL iotk_scan_dat(u, 'MOL_DENSITY', mol%density, attr=attr)
      IF (mol%density <= 0.0_DP) THEN
        CALL errore('read_mol_v1', 'molecular density is not positive @MOL_DENSITY', 1)
      END IF
      !
      CALL iotk_scan_attr(attr, 'UNITS', units, default='1/bohr^3')
      DO i = 1, LEN_TRIM(units)
        units(i:i) = capital(units(i:i))
      END DO
      !
      SELECT CASE (TRIM(units))
      CASE ('1/BOHR^3')
        ! NOP
      CASE ('G/CM^3')
        mol%density = (mol%density / mol%mass) / BOHRm3_TO_MOLCMm3
      CASE ('MOL/L')
        mol%density = mol%density / BOHRm3_TO_MOLLm1
      CASE DEFAULT
        CALL errore('read_mol_v1', 'incorrect units @MOL_DENSITY', 1)
      END SELECT
      !
      mol%subdensity = mol%density
      !
    END SUBROUTINE read_mol_density
    !
    SUBROUTINE read_mol_permittivity(u, mol)
      IMPLICIT NONE
      INTEGER,        INTENT(IN)    :: u    ! i/o unit
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      LOGICAL :: found
      !
      CALL iotk_scan_dat(u, 'MOL_PERMITTIVITY', mol%permittivity, FOUND=found)
      IF (.NOT. found) THEN
        mol%permittivity = 0.0_DP
      END IF
      !
      IF (mol%permittivity < 0.0_DP) THEN
        CALL infomsg('read_mol_v1', 'molecular permittivity is negative @MOL_PERMITTIVITY')
        mol%permittivity = 0.0_DP
      END IF
    END SUBROUTINE read_mol_permittivity
    !
    SUBROUTINE read_mol_element(u, mol)
      IMPLICIT NONE
      INTEGER,        INTENT(IN)    :: u    ! i/o unit
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      INTEGER :: iatom
      !
      IF (ASSOCIATED(mol%aname)) THEN
        DEALLOCATE(mol%aname)
      END IF
      !
      ALLOCATE(mol%aname(mol%natom))
      CALL iotk_scan_dat(u, 'MOL_ELEMENT', mol%aname)
      !
      DO iatom = 1, mol%natom
        mol%aname(iatom) = ADJUSTL(mol%aname(iatom))
        IF (LEN_TRIM(mol%aname(iatom)) < 1) THEN
          CALL errore('read_mol_v1', 'atomic name is empty @MOL_ELEMENT', iatom)
        END IF
      END DO
    END SUBROUTINE read_mol_element
    !
    SUBROUTINE read_mol_xyz(u, mol)
      IMPLICIT NONE
      INTEGER,        INTENT(IN)    :: u    ! i/o unit
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      CHARACTER(LEN=iotk_attlenx) :: attr
      CHARACTER(LEN=16)           :: units
      INTEGER                     :: i
      !
      CHARACTER(LEN=1), EXTERNAL :: capital
      !
      IF (ASSOCIATED(mol%coord)) THEN
        DEALLOCATE(mol%coord)
      END IF
      !
      ALLOCATE(mol%coord(3, mol%natom))
      CALL iotk_scan_dat(u, 'MOL_XYZ', mol%coord, attr=attr)
      !
      CALL iotk_scan_attr(attr, 'UNITS', units, default='bohr')
      DO i = 1, LEN_TRIM(units)
         units(i:i) = capital(units(i:i))
      END DO
      !
      SELECT CASE (TRIM(units))
      CASE ('BOHR')
        ! NOP
      CASE ('ANGSTROM')
        mol%coord = mol%coord / BOHR_RADIUS_ANGS
      CASE DEFAULT
        CALL errore('read_mol_v1', 'incorrect units @MOL_XYZ', 1)
      END SELECT
    END SUBROUTINE read_mol_xyz
    !
    SUBROUTINE read_mol_charge(u, mol)
      IMPLICIT NONE
      INTEGER,        INTENT(IN)    :: u    ! i/o unit
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      IF (.NOT. mol%has_charge) THEN
        RETURN
      END IF
      !
      IF (ASSOCIATED(mol%charge)) THEN
        DEALLOCATE(mol%charge)
      END IF
      !
      ALLOCATE(mol%charge(mol%natom))
      CALL iotk_scan_dat(u, 'MOL_CHARGE', mol%charge)
    END SUBROUTINE read_mol_charge
    !
    SUBROUTINE read_mol_lj(u, mol)
      IMPLICIT NONE
      INTEGER,        INTENT(IN)    :: u    ! i/o unit
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      CHARACTER(LEN=iotk_attlenx) :: attr
      CHARACTER(LEN=16)           :: units
      INTEGER                     :: i
      INTEGER                     :: iatom
      !
      CHARACTER(LEN=1), EXTERNAL :: capital
      !
      IF (.NOT. mol%has_lj) THEN
        RETURN
      END IF
      !
      CALL iotk_scan_begin(u, 'MOL_LJ')
      !
      ! ... read epsilon parameter
      IF (ASSOCIATED(mol%ljeps)) THEN
        DEALLOCATE(mol%ljeps)
      END IF
      !
      ALLOCATE(mol%ljeps(mol%natom))
      CALL iotk_scan_dat(u, 'MOL_EPSILON', mol%ljeps, attr=attr)
      DO iatom = 1, mol%natom
        IF (mol%ljeps(iatom) < 0.0_DP) THEN
          CALL errore('read_mol_v1', 'L.J.-epsilon is negative @MOL_EPSILON', iatom)
        END IF
      END DO
      !
      CALL iotk_scan_attr(attr, 'UNITS', units, default='rydberg')
      DO i = 1, LEN_TRIM(units)
        units(i:i) = capital(units(i:i))
      END DO
      !
      SELECT CASE (TRIM(units))
      CASE ('RYDBERG')
        ! NOP
      CASE ('RY')
        ! NOP
      CASE ('HA')
        mol%ljeps = mol%ljeps * 2.0_DP
      CASE ('HARTREE')
        mol%ljeps = mol%ljeps * 2.0_DP
      CASE ('EV')
        mol%ljeps = mol%ljeps / RYTOEV
      CASE ('KJ/MOL')
        mol%ljeps = mol%ljeps / RY_TO_KJMOLm1
      CASE ('KCAL/MOL')
        mol%ljeps = mol%ljeps / RY_TO_KCALMOLm1
      CASE DEFAULT
        CALL errore('read_mol_v1', 'incorrect units @MOL_EPSILON', 1)
      END SELECT
      !
      ! ... read sigma parameter
      IF (ASSOCIATED(mol%ljsig)) THEN
        DEALLOCATE(mol%ljsig)
      END IF
      !
      ALLOCATE(mol%ljsig(mol%natom))
      CALL iotk_scan_dat(u, 'MOL_SIGMA', mol%ljsig, attr=attr)
      DO iatom = 1, mol%natom
        IF (mol%ljsig(iatom) < 0.0_DP) THEN
          CALL errore('read_mol_v1', 'L.J.-sigma is negative @MOL_SIGMA', iatom)
        END IF
      END DO
      !
      CALL iotk_scan_attr(attr, 'UNITS', units, default='bohr')
      DO i = 1, LEN_TRIM(units)
        units(i:i) = capital(units(i:i))
      END DO
      !
      SELECT CASE (TRIM(units))
      CASE ('BOHR')
        ! NOP
      CASE ('ANGSTROM')
        mol%ljsig = mol%ljsig / BOHR_RADIUS_ANGS
      CASE DEFAULT
        CALL errore('read_mol_v1', 'incorrect units @MOL_SIGMA', 1)
      END SELECT
      !
      CALL iotk_scan_end(u, 'MOL_LJ')
    END SUBROUTINE read_mol_lj
    !
  END SUBROUTINE read_mol_v1
  !
END MODULE read_mol_module
