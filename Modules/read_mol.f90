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
  USE upf_utils,      ONLY : version_compare
#if defined(__fox)
  USE FoX_dom
#else
  USE dom
#endif
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
    TYPE(Node)      , POINTER                 :: doc
    TYPE(Node)      , POINTER                 :: u         ! pointer to root DOM node
    TYPE(DOMException)                        :: ex
    CHARACTER(LEN=10240)                      :: strbuf
    CHARACTER(LEN=256)                        :: linebuf
    !
    ierr = 0
    !
    IF (.NOT. PRESENT(unit)) THEN
      IF (.NOT. PRESENT(filename)) THEN
        CALL errore('read_mol', 'You have to specify at least one between filename and unit', 1)
      ELSE
        doc => parseFile(filename, EX=ex)
        ierr = getExceptionCode(ex)
      END IF
    ELSE
      strbuf = ''
      DO
        READ(unit, '(a)', iostat=ierr) linebuf
        IF (ierr /= 0) EXIT
        strbuf = TRIM(strbuf) // new_line('a') // TRIM(linebuf)
      END DO
      IF (ierr < 0) THEN
        doc => parseString(TRIM(strbuf), EX=ex)
        ierr = getExceptionCode(ex)
      END IF
    END IF
    !
    IF (ierr > 0) THEN
      CALL errore('read_mol', 'Cannot open file: ' // TRIM(filename), 1)
    ELSE
      u => getFirstChild(doc)
      CALL read_mol_v1(u, mol, ierr)
    END IF
    !
    IF (ierr > 0) THEN
      CALL deallocate_molecule(mol)
    END IF
    !
    CALL destroy(doc)
    !
  END SUBROUTINE read_mol
  !
  !--------------------------------------------------------------------------
  SUBROUTINE read_mol_v1(u, mol, ierr)
    !--------------------------------------------------------------------------
    !
    ! ... read molecule in MOL format(v.1), using FoX_dom
    !
    IMPLICIT NONE
    !
    TYPE(Node),POINTER,INTENT(IN)    :: u    ! pointer to root DOM node
    TYPE(molecule),    INTENT(INOUT) :: mol  ! the molecular data
    INTEGER, OPTIONAL, INTENT(OUT)   :: ierr ! /= 0 if something went wrong
    !
    CHARACTER(LEN=64)           :: root
    CHARACTER(LEN=256)          :: attr
    INTEGER                     :: ierr_
    TYPE(DOMException)          :: ex
    LOGICAL                     :: found
    CHARACTER(len=6), PARAMETER :: max_version = '1.0.0'
    !
    LOGICAL, EXTERNAL :: matches
    !
    ! ... check DOM
    root = getTagname(u, EX=ex)
    ierr_ = getExceptionCode(ex)
    !
    IF((ABS(ierr_) > 0) .OR. (.NOT. matches('MOL', root))) THEN
      !
      IF(.NOT. PRESENT(ierr)) THEN
        CALL errore('read_mol_v1', 'Cannot open MOL file.', 1)
      END IF
      ierr = 1
      RETURN
    END IF
    !
    CALL extractDataAttribute(u, 'version', mol%nv)
    IF (version_compare(mol%nv, max_version) == 'newer') THEN
      CALL errore('read_mol_v1', 'Unknown MOL format version: ' // TRIM(mol%nv), 1)
    END IF
    !
    ! ... read header
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
    IF (PRESENT(ierr)) ierr = 0
    !
  CONTAINS
    !
    SUBROUTINE read_mol_header(u, mol)
      IMPLICIT NONE
      TYPE(Node),POINTER,INTENT(IN) :: u    ! parent node pointer
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      TYPE(Node), POINTER         :: hdrNode
      CHARACTER(LEN=256)          :: attr
      TYPE(DOMException)          :: ex
      INTEGER                     :: ios
      !
      hdrNode => item(getElementsByTagname(u, 'MOL_HEADER'), 0)
      !
      IF (hasAttribute(hdrNode, 'author')) THEN
        CALL extractDataAttribute(hdrNode, 'author', mol%author)
      ELSE
        mol%author = 'anonymous'
      END IF
      IF (hasAttribute( hdrNode, 'date')) THEN
        CALL extractDataAttribute(hdrNode, 'date', mol%date)
      ELSE
        mol%date = ' '
      END IF
      IF (hasAttribute( hdrNode, 'comment')) THEN
        CALL extractDataAttribute(hdrNode, 'comment', mol%comment)
      ELSE
        mol%comment = ' '
      END IF
      !
      CALL extractDataAttribute(hdrNode, 'formula', mol%formula)
      IF (LEN_TRIM(mol%formula) < 1) THEN
        CALL errore('read_mol_v1', 'formula is empty @MOL_HEADER', 1)
      END IF
      mol%name = ADJUSTL(mol%formula)
      !
      CALL extractDataAttribute(hdrNode, 'number_of_atoms', mol%natom)
      IF (mol%natom < 1) THEN
        CALL errore('read_mol_v1', 'number_of_atoms is not positive @MOL_HEADER', 1)
      END IF
      !
      IF (hasAttribute(hdrNode, 'has_charge')) THEN
        CALL extractDataAttribute(hdrNode, 'has_charge', mol%has_charge, iostat=ios)
        IF (ios /= 0) THEN
          CALL extractDataAttribute (hdrNode, 'has_charge', attr)
          mol%has_charge = (INDEX(attr, 'T') > 0)
        END IF
      ELSE
        mol%has_charge = .FALSE.
      END IF
      IF (hasAttribute(hdrNode, 'has_lj')) THEN
        CALL extractDataAttribute(hdrNode, 'has_lj', mol%has_lj, iostat=ios)
        IF (ios /= 0) THEN
          CALL extractDataAttribute (hdrNode, 'has_lj', attr)
          mol%has_lj = (INDEX(attr, 'T') > 0)
        END IF
      ELSE
        mol%has_lj = .FALSE.
      END IF
    END SUBROUTINE read_mol_header
    !
    SUBROUTINE read_mol_mass(u, mol)
      IMPLICIT NONE
      TYPE(Node),POINTER,INTENT(IN) :: u    ! parent node pointer
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      TYPE(Node), POINTER         :: masNode
      CHARACTER(LEN=16)           :: units
      INTEGER                     :: i
      !
      CHARACTER(LEN=1), EXTERNAL :: capital
      !
      masNode => item(getElementsByTagname(u, 'MOL_MASS'), 0)
      !
      CALL extractDataContent(masNode, mol%mass)
      IF (mol%mass <= 0.0_DP) THEN
        CALL errore('read_mol_v1', 'molecular mass is not positive @MOL_MASS', 1)
      END IF
      !
      IF (hasAttribute(masNode, 'UNITS')) THEN
        CALL extractDataAttribute(masNode, 'UNITS', units)
      ELSE
        units = 'a.m.u.'
      END IF
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
      TYPE(Node),POINTER,INTENT(IN) :: u    ! parent node pointer
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      TYPE(Node), POINTER         :: denNode
      CHARACTER(LEN=16)           :: units
      INTEGER                     :: i
      !
      CHARACTER(LEN=1), EXTERNAL :: capital
      !
      denNode => item(getElementsByTagname(u, 'MOL_DENSITY'), 0)
      !
      CALL extractDataContent(denNode, mol%density)
      IF (mol%density <= 0.0_DP) THEN
        CALL errore('read_mol_v1', 'molecular density is not positive @MOL_DENSITY', 1)
      END IF
      !
      IF (hasAttribute(denNode, 'UNITS')) THEN
        CALL extractDataAttribute(denNode, 'UNITS', units)
      ELSE
        units = '1/bohr^3'
      END IF
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
      TYPE(Node),POINTER,INTENT(IN) :: u    ! parent node pointer
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      TYPE(NodeList), POINTER     :: perNodes
      !
      perNodes => getElementsByTagname(u, 'MOL_PERMITTIVITY')
      !
      IF (getLength(perNodes) > 0) THEN
        CALL extractDataContent(item(perNodes, 0), mol%permittivity)
      ELSE
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
      TYPE(Node),POINTER,INTENT(IN) :: u    ! parent node pointer
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      TYPE(Node), POINTER         :: eleNode
      INTEGER :: iatom
      !
      eleNode => item(getElementsByTagname(u, 'MOL_ELEMENT'), 0)
      !
      IF (ASSOCIATED(mol%aname)) THEN
        DEALLOCATE(mol%aname)
      END IF
      !
      ALLOCATE(mol%aname(mol%natom))
      CALL extractDataContent(eleNode, mol%aname)
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
      TYPE(Node),POINTER,INTENT(IN) :: u    ! parent node pointer
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      TYPE(Node), POINTER         :: xyzNode
      CHARACTER(LEN=16)           :: units
      INTEGER                     :: i
      !
      CHARACTER(LEN=1), EXTERNAL :: capital
      !
      xyzNode => item(getElementsByTagname(u, 'MOL_XYZ'), 0)
      !
      IF (ASSOCIATED(mol%coord)) THEN
        DEALLOCATE(mol%coord)
      END IF
      !
      ALLOCATE(mol%coord(3, mol%natom))
      CALL extractDataContent(xyzNode, mol%coord)
      !
      IF (hasAttribute(xyzNode, 'UNITS')) THEN
        CALL extractDataAttribute(xyzNode, 'UNITS', units)
      ELSE
        units = 'bohr'
      END IF
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
      TYPE(Node),POINTER,INTENT(IN) :: u    ! parent node pointer
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      TYPE(Node), POINTER         :: chrNode
      !
      IF (.NOT. mol%has_charge) THEN
        RETURN
      END IF
      !
      chrNode => item(getElementsByTagname(u, 'MOL_CHARGE'), 0)
      !
      IF (ASSOCIATED(mol%charge)) THEN
        DEALLOCATE(mol%charge)
      END IF
      !
      ALLOCATE(mol%charge(mol%natom))
      CALL extractDataContent(chrNode, mol%charge)
    END SUBROUTINE read_mol_charge
    !
    SUBROUTINE read_mol_lj(u, mol)
      IMPLICIT NONE
      TYPE(Node),POINTER,INTENT(IN) :: u    ! parent node pointer
      TYPE(molecule), INTENT(INOUT) :: mol  ! the molecular data
      !
      TYPE(Node), POINTER         :: ljNode
      TYPE(Node), POINTER         :: epsNode
      TYPE(Node), POINTER         :: sigNode
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
      ljNode => item(getElementsByTagname(u, 'MOL_LJ'), 0)
      epsNode => item(getElementsByTagname(ljNode, 'MOL_EPSILON'), 0)
      sigNode => item(getElementsByTagname(ljNode, 'MOL_SIGMA'), 0)
      !
      ! ... read epsilon parameter
      IF (ASSOCIATED(mol%ljeps)) THEN
        DEALLOCATE(mol%ljeps)
      END IF
      !
      ALLOCATE(mol%ljeps(mol%natom))
      CALL extractDataContent(epsNode, mol%ljeps)
      DO iatom = 1, mol%natom
        IF (mol%ljeps(iatom) < 0.0_DP) THEN
          CALL errore('read_mol_v1', 'L.J.-epsilon is negative @MOL_EPSILON', iatom)
        END IF
      END DO
      !
      IF (hasAttribute(epsNode, 'UNITS')) THEN
        CALL extractDataAttribute(epsNode, 'UNITS', units)
      ELSE
        units = 'rydberg'
      END IF
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
      CALL extractDataContent(sigNode, mol%ljsig)
      DO iatom = 1, mol%natom
        IF (mol%ljsig(iatom) < 0.0_DP) THEN
          CALL errore('read_mol_v1', 'L.J.-sigma is negative @MOL_SIGMA', iatom)
        END IF
      END DO
      !
      IF (hasAttribute(sigNode, 'UNITS')) THEN
        CALL extractDataAttribute(sigNode, 'UNITS', units)
      ELSE
        units = 'bohr'
      END IF
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
    END SUBROUTINE read_mol_lj
    !
  END SUBROUTINE read_mol_v1
  !
END MODULE read_mol_module
