!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE hp_read_chi()
!------------------------------------------------------------
  !
  ! Read chi0 and chi from file.
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix, tmp_dir
  USE ldaU_hp,       ONLY : nath_sc, nath, nah_pert, chi0, chi, &
                            todo_atom, tmp_dir_hp
  !
  IMPLICIT NONE
  !
  INTEGER :: nb
  INTEGER :: iunitchi ! unit number
  CHARACTER(len=50) :: filenamechi
  CHARACTER(len=256) :: tempfile
  CHARACTER(len=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  INTEGER, EXTERNAL :: find_free_unit
  !
  tmp_dir = tmp_dir_hp
  !
  chi0(:,:) = 0.0d0
  chi(:,:)  = 0.0d0
  !
  DO nb = 1, nat
     !
     IF ( todo_atom(nb) ) THEN 
        !
        iunitchi = find_free_unit()
        filenamechi = TRIM(prefix) // TRIM(".chi.pert_") // TRIM(int_to_char(nb)) // TRIM(".dat")
        tempfile = TRIM(tmp_dir) // TRIM(filenamechi)
        !
        INQUIRE (file = tempfile, exist = exst)
        IF (.NOT.exst) THEN
           WRITE( stdout, * ) "    WARNING: " // TRIM(filenamechi) // " does not exist !!!"
           WRITE( stdout, * ) "    Check the folder: ", TRIM(tmp_dir)
           CALL errore('hp_read_chi','Missing file',1)
        ENDIF
        !
        OPEN(iunitchi, file = tempfile, form = 'formatted', status = 'unknown')
        !
        CALL read_chi(nb)
        !
        CLOSE(iunitchi)
        !
     ENDIF
     !
  ENDDO
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE read_chi(nb)
  ! 
  ! Read nb-th column of chi0 and chi
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nb
  INTEGER :: na, na_, nb_
  !
  ! Read chi0
  READ(iunitchi, *)
  DO na = 1, nath_sc
     READ(iunitchi,*) na_, nb_, chi0(na, nb)
  ENDDO
  READ(iunitchi,*)
  !
  ! Read chi
  READ(iunitchi,*)
  DO na = 1, nath_sc
     READ(iunitchi,*) na_, nb_, chi(na, nb)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE read_chi
  !
END SUBROUTINE hp_read_chi
