!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE startup( nd_nmbr, code, version )
  !----------------------------------------------------------------------------
  !
  ! ... Wrapper routine, initializes MPI and various other things
  !
  USE control_flags, ONLY: use_task_groups, ortho_para
  USE mp_global,     ONLY: mp_startup
  USE environment,   ONLY: environment_start
  !
  CHARACTER (LEN=6)  :: nd_nmbr, version
  CHARACTER (LEN=9)  :: code
  !
#ifdef __PARA
  CALL mp_startup ( use_task_groups, ortho_para )
#endif
  CALL environment_start ( code )
  !
  RETURN
  !     
END SUBROUTINE startup
