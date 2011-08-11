!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE read_input
   !---------------------------------------------------------------------------
   !
   USE kinds,     ONLY: DP
   !
   IMPLICIT NONE
   SAVE
   !
   PRIVATE
   PUBLIC :: read_input_file, has_been_read
   !
   LOGICAL :: has_been_read = .FALSE.
   !
   CONTAINS
   !
   !-------------------------------------------------------------------------
   SUBROUTINE read_input_file ( prog )
     !-------------------------------------------------------------------------
     !
     USE read_namelists_module, ONLY : read_namelists
     USE read_cards_module,     ONLY : read_cards
     USE io_global,             ONLY : ionode, ionode_id
     USE xml_input,             ONLY : xml_input_dump
     USE read_xml_module,       ONLY : read_xml
     USE mp,                    ONLY : mp_bcast
     USE mp_global,             ONLY : intra_image_comm
     USE iotk_module,           ONLY : iotk_attlenx
     USE open_close_input_file_interf, ONLY : open_input_file, close_input_file
     !
     IMPLICIT NONE
     !
     CHARACTER(LEN=2), INTENT (IN) :: prog
     CHARACTER(LEN=iotk_attlenx) :: attr
     LOGICAL :: xmlinput
     !
     !
     IF ( ionode ) THEN
        IF ( prog == 'CP' ) CALL xml_input_dump()
        CALL open_input_file( xmlinput, attr) 
     END IF
     !
     ! bcast of xmlinput and attr is needed because
     ! it is only the open statement inside read_cards 
     ! and read_namelist (in Modules) that has if(ionode)
     !
     CALL mp_bcast( xmlinput, ionode_id, intra_image_comm )
     CALL mp_bcast( attr, ionode_id, intra_image_comm )
     !
     IF ( xmlinput ) THEN
        !
        CALL read_xml ( prog, attr )
        !
     ELSE
        !
        ! ... Read NAMELISTS 
        !
        CALL read_namelists( prog )
        !
        ! ... Read CARDS 
        !
        CALL read_cards ( prog )
        !
     END IF
     IF ( ionode) CALL close_input_file( xmlinput )
     !
     has_been_read = .TRUE.
     !
     RETURN
     !
   END SUBROUTINE read_input_file
  !
END MODULE read_input
