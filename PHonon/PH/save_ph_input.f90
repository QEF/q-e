!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE save_ph
  !----------------------------------------------------------------------------
  !
  ! ... this module contains methods to read and write data saved by the
  !     phonon code to restart smoothly
  !
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: save_ph_input_variables, restore_ph_input_variables, &
            clean_input_variables
  !
  INTEGER, PRIVATE :: nat_todo_save
  INTEGER, ALLOCATABLE, PRIVATE :: atomo_save(:)
  CHARACTER(LEN=256), PUBLIC :: tmp_dir_save
  !
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE save_ph_input_variables()
      !------------------------------------------------------------------------
      !
      USE ions_base,  ONLY : nat
      USE partial,    ONLY : atomo, nat_todo
      USE control_ph, ONLY : search_sym_save, search_sym
      !
      IMPLICIT NONE
      !
      ALLOCATE(atomo_save(nat))
      nat_todo_save=nat_todo
      atomo_save=atomo
      search_sym_save=search_sym

      RETURN
    END SUBROUTINE save_ph_input_variables
    !
    SUBROUTINE restore_ph_input_variables(  )
      !------------------------------------------------------------------------
      !
      USE io_files,   ONLY : tmp_dir
      USE ions_base,  ONLY : nat
      USE partial,    ONLY : atomo, nat_todo
      USE control_ph, ONLY : search_sym_save, search_sym
      !
      IMPLICIT NONE
      !
      nat_todo=nat_todo_save
      atomo=atomo_save
      tmp_dir=tmp_dir_save
      search_sym = search_sym_save

      RETURN
    END SUBROUTINE restore_ph_input_variables

    SUBROUTINE clean_input_variables()
    IMPLICIT NONE

    DEALLOCATE(atomo_save)

    RETURN
    END SUBROUTINE clean_input_variables
    !
END MODULE save_ph
