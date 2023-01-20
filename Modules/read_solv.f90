!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE read_solv_module
  !--------------------------------------------------------------------------
  !
  ! ... read molecular files as solvents
  !
  USE io_files,        ONLY : pseudo_dir, pseudo_dir_cur, molfile
  USE io_global,       ONLY : stdout, ionode
  USE kinds,           ONLY : DP
  USE mp,              ONLY : mp_sum
  USE mp_images,       ONLY : intra_image_comm
  USE read_mol_module, ONLY : read_mol, clean_mol_readables, &
                            & without_mass_         => without_mass,         &
                            & without_density_      => without_density,      &
                            & without_permittivity_ => without_permittivity, &
                            & without_element_      => without_element,      &
                            & without_xyz_          => without_xyz,          &
                            & without_charge_       => without_charge,       &
                            & without_lj_           => without_lj
  USE solvmol,         ONLY : nsolV, solVs, allocate_solVs, deallocate_solVs, set_index_of_solVs
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  PUBLIC :: read_solvents
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE read_solvents(printout, without_mass, without_density, without_permittivity, &
                         & without_element, without_xyz, without_charge, without_lj)
    !--------------------------------------------------------------------------
    !
    ! ... read molecular files to put into `molecule' structure,
    ! ... and initialize `solvent' module.
    !
    IMPLICIT NONE
    !
    LOGICAL, OPTIONAL, INTENT(IN) :: printout
    LOGICAL, OPTIONAL, INTENT(IN) :: without_mass
    LOGICAL, OPTIONAL, INTENT(IN) :: without_density
    LOGICAL, OPTIONAL, INTENT(IN) :: without_permittivity
    LOGICAL, OPTIONAL, INTENT(IN) :: without_element
    LOGICAL, OPTIONAL, INTENT(IN) :: without_xyz
    LOGICAL, OPTIONAL, INTENT(IN) :: without_charge
    LOGICAL, OPTIONAL, INTENT(IN) :: without_lj
    !
    INTEGER :: isolV
    INTEGER :: ios
    INTEGER :: iunsv
    INTEGER :: ismol
    INTEGER :: nsolV_
    CHARACTER(len=256) :: file_solvent
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    LOGICAL, SAVE :: printout_ = .FALSE.
    !
    CALL clean_mol_readables()
    !
    ! ... set parameters to avoid reading
    IF (PRESENT(without_mass)) THEN
      without_mass_ = without_mass
    END IF
    !
    IF (PRESENT(without_density)) THEN
      without_density_ = without_density
    END IF
    !
    IF (PRESENT(without_permittivity)) THEN
      without_permittivity_ = without_permittivity
    END IF
    !
    IF (PRESENT(without_element)) THEN
      without_element_ = without_element
    END IF
    !
    IF (PRESENT(without_xyz)) THEN
      without_xyz_ = without_xyz
    END IF
    !
    IF (PRESENT(without_charge)) THEN
      without_charge_ = without_charge
    END IF
    !
    IF (PRESENT(without_lj)) THEN
      without_lj_ = without_lj
    END IF
    !
    ! ... initialize solvent module
    iunsv = find_free_unit()
    !
    IF (ALLOCATED(solVs)) THEN
      IF (SIZE(solVs) /= nsolV) THEN
        nsolV_ = nsolV
        CALL deallocate_solVs()
        nsolV = nsolV_
      END IF
    END IF
    !
    IF (.NOT. ALLOCATED(solVs)) THEN
      CALL allocate_solVs()
    END IF
    !
    IF (PRESENT(printout)) THEN
      printout_ = printout
    END IF
    !
    IF (ionode .AND. printout_) THEN
      WRITE(stdout, "(//,3X,'Solvent Molecular Parameters',/, &
                    &    3X,'----------------------------' )" )
    END IF
    !
    ! ... read every solvent
    DO isolV = 1, nsolV
      !
      ! try first pseudo_dir_cur if set: in case of restart from file,
      ! this is where PP files should be located
      !
      ios = 1
      IF (pseudo_dir_cur /= ' ') THEN
        file_solvent = TRIM(pseudo_dir_cur) // TRIM(molfile(isolV))
        OPEN(unit=iunsv, file=file_solvent, status='old', form='formatted', action='read', iostat=ios)
        CALL mp_sum(ios, intra_image_comm)
        IF (ios /= 0) THEN
          CALL infomsg('read_solvents', 'file ' // TRIM(file_solvent) // ' not found')
        END IF
        !
        ! file not found? no panic (yet): if the restart file is not visible
        ! to all processors, this may happen. Try the original location
      END IF
      !
      ! try the original location pseudo_dir, as set in input
      ! (it should already contain a slash at the end)
      !
      IF (ios /= 0) THEN
        file_solvent = TRIM(pseudo_dir) // TRIM(molfile(isolV))
        OPEN(unit=iunsv, file=file_solvent, status='old', form='formatted', action='read', iostat=ios)
        CALL mp_sum(ios, intra_image_comm)
        CALL errore('read_solvents', 'file ' // TRIM(file_solvent) // ' not found', ABS(ios))
      END IF
      !
      ! start reading - MOL first: the MOL format is detected via the
      ! presence of the keyword '<MOL>' at the beginning of the file
      !
      IF (ionode .AND. printout_) THEN
        WRITE(stdout, "(/,3X,'Reading molecule of # ',I2, &
                      & ' from file :',/,3X,A)") isolV, TRIM(file_solvent)
      END IF
      !
      ! reading molecule
      !
      CALL read_mol(solVs(isolV), ismol, unit=iunsv)
      !
      IF (ismol == 0) THEN
        IF (ionode .AND. printout_) THEN
          WRITE(stdout, "(3X,'file type is MOL v.',i1)") (ismol + 1)
        END IF
        !
      ELSE
        CALL errore('read_solvents', 'cannot read file ' // TRIM(file_solvent), MAX(ABS(ismol), 1))
      END IF
      !
      ! end of reading
      !
      CLOSE(iunsv)
      !
    END DO
    !
    ! more initializations
    !
    CALL set_index_of_solVs()
    !
    CALL clean_mol_readables()
    !
  END SUBROUTINE read_solvents
  !
END MODULE read_solv_module
