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
   SUBROUTINE read_input_file ( prog, input_file_ )
     !-------------------------------------------------------------------------
     !
     USE input_parameters,      ONLY : reset_input_checks
     USE read_namelists_module, ONLY : read_namelists
     USE read_cards_module,     ONLY : read_cards
     USE io_global,             ONLY : ionode, ionode_id, qestdin
     USE xml_input,             ONLY : xml_input_dump
     USE read_xml_module,       ONLY : read_xml
     USE mp,                    ONLY : mp_bcast
     USE mp_images,             ONLY : intra_image_comm
     USE iotk_module,           ONLY : iotk_attlenx
     USE open_close_input_file, ONLY : open_input_file, close_input_file
     !
     IMPLICIT NONE
     !
     CHARACTER(LEN=2), INTENT (IN) :: prog
     CHARACTER(LEN=*), INTENT (IN) :: input_file_
     !
     CHARACTER(LEN=iotk_attlenx) :: attr
     LOGICAL :: xmlinput
     INTEGER :: ierr
     !
     IF ( ionode ) THEN
        IF ( prog == 'CP' ) CALL xml_input_dump()
        ierr = open_input_file( input_file_, xmlinput, attr) 
     END IF
     !
     CALL mp_bcast( ierr, ionode_id, intra_image_comm )
     IF ( ierr > 0 ) CALL errore('read_input', 'opening input file',ierr)
     CALL mp_bcast( xmlinput, ionode_id, intra_image_comm )
     CALL mp_bcast( attr, ionode_id, intra_image_comm )
     !
     CALL reset_input_checks () 
     !
     IF ( xmlinput ) THEN
        !
        CALL read_xml ( prog, attr )
        !
     ELSE
        !
        ! ... Read NAMELISTS 
        !
        CALL read_namelists( prog, qestdin )
        !
        ! ... Read CARDS 
        !
        CALL read_cards ( prog, qestdin )
        !
     END IF
     IF ( ionode) ierr = close_input_file( )
     !
     has_been_read = .TRUE.
     !
     RETURN
     !
   END SUBROUTINE read_input_file
  !
END MODULE read_input
