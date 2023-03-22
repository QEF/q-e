!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE molecule_types
  !--------------------------------------------------------------------------
  !
  ! ... this module defines types of molecule.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define constants
  INTEGER, PARAMETER :: LEN_NV      = 11
  INTEGER, PARAMETER :: LEN_AUTHOR  = 80
  INTEGER, PARAMETER :: LEN_DATE    = 80
  INTEGER, PARAMETER :: LEN_COMMENT = 80
  INTEGER, PARAMETER :: LEN_FORMULA = 256
  INTEGER, PARAMETER :: LEN_NAME    = 16
  INTEGER, PARAMETER :: LEN_ANAME   = 8
  !
  ! ... define type of molecule
  TYPE molecule
    CHARACTER(LEN=LEN_NV)             :: nv           ! MOL file three-digit version i.e. 1.0.0
    CHARACTER(LEN=LEN_AUTHOR)         :: author       ! author of molecular file
    CHARACTER(LEN=LEN_DATE)           :: date         ! generation date
    CHARACTER(LEN=LEN_COMMENT)        :: comment      ! author's comment
    CHARACTER(LEN=LEN_FORMULA)        :: formula      ! chemical formula of molecule
    CHARACTER(LEN=LEN_NAME)           :: name         ! name of molecule
    INTEGER                           :: natom        ! number of atoms in molecule
    REAL(DP)                          :: mass         ! mass of molecule           (in a.m.u.)
    REAL(DP)                          :: density      ! density of molecule        (in 1/bohr^3)
    REAL(DP)                          :: subdensity   ! second density of molecule (in 1/bohr^3)
    REAL(DP)                          :: permittivity ! relative permittivity
    REAL(DP)                          :: dipole       ! dipole moment              (in e*bohr)
    LOGICAL                           :: has_charge   ! if .true. includes charge
    LOGICAL                           :: has_lj       ! if .true. includes Lennard-Jones
    LOGICAL                           :: is_polar     ! if .true. this is a polar molecule
    CHARACTER(LEN=LEN_ANAME), POINTER :: aname(:)     ! name of atoms
    REAL(DP),                 POINTER :: coord(:,:)   ! xyz-coordinate of atoms (in bohr)
    REAL(DP),                 POINTER :: charge(:)    ! charge of atoms         (in e)
    REAL(DP),                 POINTER :: ljeps(:)     ! L.J.-epsilon of atoms   (in Ry)
    REAL(DP),                 POINTER :: ljsig(:)     ! L.J.-sgmma of atoms     (in bohr)
  END TYPE molecule
  !
  ! ... public components
  PUBLIC :: LEN_NAME
  PUBLIC :: LEN_ANAME
  PUBLIC :: molecule
  PUBLIC :: deallocate_molecule
  PUBLIC :: nullify_molecule
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE nullify_molecule(mol)
    !--------------------------------------------------------------------------
    !
    ! ... clear data of molecule
    !
    IMPLICIT NONE
    !
    TYPE(molecule), INTENT(INOUT) :: mol
    !
    mol%nv           = ''
    mol%author       = ''
    mol%date         = ''
    mol%comment      = ''
    mol%formula      = ''
    mol%name         = ''
    mol%natom        = 0
    mol%mass         = 0.0_DP
    mol%density      = 0.0_DP
    mol%subdensity   = 0.0_DP
    mol%permittivity = 0.0_DP
    mol%dipole       = 0.0_DP
    mol%has_charge   = .FALSE.
    mol%has_lj       = .FALSE.
    mol%is_polar     = .FALSE.
    NULLIFY(mol%aname)
    NULLIFY(mol%coord)
    NULLIFY(mol%charge)
    NULLIFY(mol%ljeps)
    NULLIFY(mol%ljsig)
    !
  END SUBROUTINE nullify_molecule
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_molecule(mol)
    !--------------------------------------------------------------------------
    !
    ! ... free memory
    !
    IMPLICIT NONE
    !
    TYPE(molecule), INTENT(INOUT) :: mol
    !
    IF (ASSOCIATED(mol%aname )) DEALLOCATE(mol%aname)
    IF (ASSOCIATED(mol%coord )) DEALLOCATE(mol%coord)
    IF (ASSOCIATED(mol%charge)) DEALLOCATE(mol%charge)
    IF (ASSOCIATED(mol%ljeps )) DEALLOCATE(mol%ljeps)
    IF (ASSOCIATED(mol%ljsig )) DEALLOCATE(mol%ljsig)
    !
  END SUBROUTINE deallocate_molecule
  !
END MODULE molecule_types
