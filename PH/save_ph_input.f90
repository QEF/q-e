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
  INTEGER, PRIVATE :: nat_todo_save, nrapp_save
  INTEGER, ALLOCATABLE, PRIVATE :: list_save(:), atomo_save(:) 
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
      USE partial,    ONLY : list, atomo, nat_todo, nrapp
      !
      IMPLICIT NONE
      !
      ALLOCATE(list_save(3*nat))
      ALLOCATE(atomo_save(nat))
      nat_todo_save=nat_todo
      nrapp_save=nrapp
      list_save=list
      atomo_save=atomo

      RETURN
    END SUBROUTINE save_ph_input_variables
    !
    SUBROUTINE restore_ph_input_variables(  )
      !------------------------------------------------------------------------
      !
      USE io_files,   ONLY : tmp_dir
      USE ions_base,  ONLY : nat
      USE partial,    ONLY : list, atomo, nat_todo, nrapp
      !
      IMPLICIT NONE
      !
      nat_todo=nat_todo_save
      nrapp=nrapp_save
      list=list_save
      atomo=atomo_save
      tmp_dir=tmp_dir_save

      RETURN
    END SUBROUTINE restore_ph_input_variables

    SUBROUTINE clean_input_variables()
    IMPLICIT NONE

    DEALLOCATE(list_save)
    DEALLOCATE(atomo_save)    

    RETURN
    END SUBROUTINE clean_input_variables
    !
END MODULE save_ph
