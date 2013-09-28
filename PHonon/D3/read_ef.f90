!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE read_ef()
  !-----------------------------------------------------------------------
  !
  ! Reads the shift of the Fermi Energy
  !
  USE pwcom
  USE d3com
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  !
  IF (degauss == 0.d0 ) RETURN
  !
  IF ( ionode ) THEN
     !
     REWIND (unit = iuef)
     READ (iuef, err = 100, iostat = ios) ef_sh
     !
     !
  END IF

100 CALL mp_bcast(ios, ionode_id, world_comm)

  CALL errore ('d3_valence', 'reading iuef', ABS (ios) )

  CALL mp_bcast( ef_sh, ionode_id, world_comm )

  RETURN
END SUBROUTINE read_ef
