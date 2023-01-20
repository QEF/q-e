!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE solvmol
  !--------------------------------------------------------------------------
  !
  ! ... this module keeps data of solvents.
  !
  USE kinds,          ONLY : DP
  USE molecule_types, ONLY : LEN_ANAME, molecule, nullify_molecule, deallocate_molecule
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... number of solvents
  INTEGER :: nsolV = 0
  !
  ! ... total number of sites in solvents
  INTEGER :: nsite_save = 0
  !
  ! ... total number of unique sites in solvents
  INTEGER :: nuniq_save = 0
  !
  ! ... array of solvents
  TYPE(molecule), ALLOCATABLE :: solVs(:)
  !
  ! ... index of sites in solvents
  INTEGER, ALLOCATABLE :: isite_to_isolV(:)    ! index of site -> index of solvent
  INTEGER, ALLOCATABLE :: isite_to_iatom(:)    ! index of site -> index of atom in a solvent
  INTEGER, ALLOCATABLE :: iuniq_to_nsite(:)    ! index of unique site -> number of same sites
  INTEGER, ALLOCATABLE :: iuniq_to_isite(:,:)  ! index of unique site -> index of site
  !
  ! ... public components
  PUBLIC :: nsolV
  PUBLIC :: solVs
  PUBLIC :: isite_to_isolV
  PUBLIC :: isite_to_iatom
  PUBLIC :: iuniq_to_nsite
  PUBLIC :: iuniq_to_isite
  PUBLIC :: allocate_solVs
  PUBLIC :: deallocate_solVs
  PUBLIC :: get_nsite_in_solVs
  PUBLIC :: get_nuniq_in_solVs
  PUBLIC :: set_index_of_solVs
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_solVs(ns)
    !--------------------------------------------------------------------------
    !
    ! ... initialize solVs
    !
    IMPLICIT NONE
    !
    INTEGER, OPTIONAL, INTENT(IN) :: ns
    !
    INTEGER :: isolV
    !
    IF (PRESENT(ns)) THEN
      nsolV = ns
    END IF
    !
    ALLOCATE(solVs(nsolV))
    !
    DO isolV = 1, nsolV
      call nullify_molecule(solVs(isolV))
    END DO
    !
    nsite_save = 0
    nuniq_save = 0
    !
  END SUBROUTINE allocate_solVs
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_solVs()
    !--------------------------------------------------------------------------
    !
    ! ... finalize solVs
    !
    IMPLICIT NONE
    !
    INTEGER :: isolV
    !
    IF (ALLOCATED(solVs)) THEN
      ! deallocate each solvent
      DO isolV = 1, nsolV
        CALL deallocate_molecule(solVs(isolV))
        CALL nullify_molecule(solVs(isolV))
      END DO
      ! deallocate solvent array
      DEALLOCATE(solVs)
    END IF
    !
    ! deallocate indexes of solvents
    IF (ALLOCATED(isite_to_isolV)) DEALLOCATE(isite_to_isolV)
    IF (ALLOCATED(isite_to_iatom)) DEALLOCATE(isite_to_iatom)
    IF (ALLOCATED(iuniq_to_nsite)) DEALLOCATE(iuniq_to_nsite)
    IF (ALLOCATED(iuniq_to_isite)) DEALLOCATE(iuniq_to_isite)
    !
    nsolV = 0
    !
    nsite_save = 0
    nuniq_save = 0
    !
  END SUBROUTINE deallocate_solVs
  !
  !--------------------------------------------------------------------------
  FUNCTION get_nsite_in_solVs() RESULT(nsite)
    !--------------------------------------------------------------------------
    !
    ! ... get total number of sites in solvents
    !
    IMPLICIT NONE
    !
    INTEGER :: nsite
    INTEGER :: isolV
    !
    IF (nsite_save > 0) THEN
      nsite = nsite_save
      RETURN
    END IF
    !
    nsite = 0
    DO isolV = 1, nsolV
      nsite = nsite + solVs(isolV)%natom
    END DO
    !
    nsite_save = nsite
    !
  END FUNCTION get_nsite_in_solVs
  !
  !--------------------------------------------------------------------------
  FUNCTION get_nuniq_in_solVs() RESULT(nuniq)
    !--------------------------------------------------------------------------
    !
    ! ... get total number of unique sites in solvents
    !
    IMPLICIT NONE
    !
    INTEGER :: nuniq
    INTEGER :: isolV
    INTEGER :: iatom1
    INTEGER :: iatom2
    INTEGER :: ipso
    CHARACTER(LEN=LEN_ANAME) :: aname1
    CHARACTER(LEN=LEN_ANAME) :: aname2
    !
    IF (nuniq_save > 0) THEN
      nuniq = nuniq_save
      RETURN
    END IF
    !
    nuniq = 0
    DO isolV = 1, nsolV
      DO iatom1 = 1, solVs(isolV)%natom
        aname1 = solVs(isolV)%aname(iatom1)
        ipso = 0
        DO iatom2 = 1, (iatom1 - 1)
          aname2 = solVs(isolV)%aname(iatom2)
          IF (aname1 == aname2) ipso = ipso + 1
        END DO
        IF (ipso < 1) THEN
          nuniq = nuniq + 1
        END IF
      END DO
    END DO
    !
    nuniq_save = nuniq
    !
  END FUNCTION get_nuniq_in_solVs
  !
  !--------------------------------------------------------------------------
  FUNCTION get_max_same_sites() RESULT(msite)
    !--------------------------------------------------------------------------
    !
    ! ... get maximum number of identical sites
    !
    IMPLICIT NONE
    !
    INTEGER :: msite
    !
    INTEGER :: nsite
    INTEGER :: isolV
    INTEGER :: iatom1
    INTEGER :: iatom2
    INTEGER :: ipso
    CHARACTER(LEN=LEN_ANAME) :: aname1
    CHARACTER(LEN=LEN_ANAME) :: aname2
    !
    msite = 1
    !
    DO isolV = 1, nsolV
      DO iatom1 = 1, solVs(isolV)%natom
        !
        aname1 = solVs(isolV)%aname(iatom1)
        !
        ! ... count same sites (in iatom2 < iatom1)
        ipso = 0
        DO iatom2 = 1, (iatom1 - 1)
          aname2 = solVs(isolV)%aname(iatom2)
          IF (aname1 == aname2) ipso = ipso + 1
        END DO
        !
        ! ... no same sites (i.e. to be unique site)
        IF (ipso < 1) THEN
          !
          ! ... count number of sites named aname1
          nsite = 1
          DO iatom2 = (iatom1 + 1), solVs(isolV)%natom
            aname2 = solVs(isolV)%aname(iatom2)
            IF (aname1 == aname2) THEN
              nsite = nsite + 1
            END IF
          END DO
          !
          ! ... set maximum of nsite
          msite = MAX(msite, nsite)
          !
        END IF
        !
      END DO
    END DO
    !
  END FUNCTION get_max_same_sites
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_index_of_solVs()
    !--------------------------------------------------------------------------
    !
    ! ... create indexes `isite_to_isolV' and `isite_to_iatom'
    !
    IMPLICIT NONE
    !
    INTEGER :: nsite
    INTEGER :: msite
    INTEGER :: isite
    INTEGER :: jsite
    INTEGER :: nuniq
    INTEGER :: iuniq
    INTEGER :: isolV
    INTEGER :: iatom1
    INTEGER :: iatom2
    INTEGER :: ipso
    CHARACTER(LEN=LEN_ANAME) :: aname1
    CHARACTER(LEN=LEN_ANAME) :: aname2
    !
    ! ... deallocate, if already allocated
    IF (ALLOCATED(isite_to_isolV)) DEALLOCATE(isite_to_isolV)
    IF (ALLOCATED(isite_to_iatom)) DEALLOCATE(isite_to_iatom)
    IF (ALLOCATED(iuniq_to_nsite)) DEALLOCATE(iuniq_to_nsite)
    IF (ALLOCATED(iuniq_to_isite)) DEALLOCATE(iuniq_to_isite)
    !
    ! ... allocate array
    nsite = get_nsite_in_solVs()
    nuniq = get_nuniq_in_solVs()
    msite = get_max_same_sites()
    ALLOCATE(isite_to_isolV(nsite))
    ALLOCATE(isite_to_iatom(nsite))
    ALLOCATE(iuniq_to_nsite(nuniq))
    ALLOCATE(iuniq_to_isite(msite, nuniq))
    !
    ! ... create indexes
    isite = 0
    iuniq = 0
    DO isolV = 1, nsolV
      DO iatom1 = 1, solVs(isolV)%natom
        !
        aname1 = solVs(isolV)%aname(iatom1)
        !
        ! ... set site's index
        isite = isite + 1
        isite_to_isolV(isite) = isolV
        isite_to_iatom(isite) = iatom1
        !
        ! ... count same sites (in iatom2 < iatom1)
        ipso = 0
        DO iatom2 = 1, (iatom1 - 1)
          aname2 = solVs(isolV)%aname(iatom2)
          IF (aname1 == aname2) ipso = ipso + 1
        END DO
        !
        ! ... no same sites (i.e. to be unique site)
        IF (ipso < 1) THEN
          !
          ! ... set unique site's index
          iuniq = iuniq + 1
          iuniq_to_nsite(iuniq) = 1
          iuniq_to_isite(1, iuniq) = isite
          jsite = isite
          DO iatom2 = (iatom1 + 1), solVs(isolV)%natom
            jsite = jsite + 1
            aname2 = solVs(isolV)%aname(iatom2)
            IF (aname1 == aname2) THEN
              iuniq_to_nsite(iuniq) = iuniq_to_nsite(iuniq) + 1
              iuniq_to_isite(iuniq_to_nsite(iuniq), iuniq) = jsite
            END IF
          END DO
          !
        END IF
        !
      END DO
    END DO
    !
  END SUBROUTINE set_index_of_solVs
  !
END MODULE solvmol
