!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_arguments()
  !-----------------------------------------------------------------------------
  !
  ! the name of each new plugin has to be set here. Its use or not
  ! is controlled by its presence as argument
  !
  !
  USE kinds,         ONLY : DP
  !
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp_global,     ONLY : intra_image_comm
  !
  USE plugin_flags
  !
  USE mp,            ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER  :: iiarg, nargs, iargc, ierr
  CHARACTER (len=50) :: arg
  CHARACTER (len=256) :: np
  !
  !
#if defined(__ABSOFT)
#   define getarg getarg_
#   define iargc  iargc_
#endif
  !
  nargs = iargc()
  !
  !
  DO iiarg = 1, ( nargs )
    CALL getarg( iiarg, plugin_name)
    IF((TRIM(plugin_name)=='-plumed').or.(TRIM(plugin_name)=='-PLUMED')) THEN
      use_plumed = .true.
      CALL mp_bcast( use_plumed, intra_image_comm )
    ENDIF
  ENDDO
  !
  !
  RETURN
  !
END SUBROUTINE plugin_arguments
