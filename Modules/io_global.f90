!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE io_global

        IMPLICIT NONE

        PRIVATE
        SAVE

        PUBLIC :: io_global_start, io_global_getionode
        PUBLIC :: ionode, ionode_id

        INTEGER :: ionode_id = 0
        LOGICAL :: ionode = .TRUE.

        LOGICAL :: first = .TRUE.
      
      CONTAINS

        SUBROUTINE io_global_start(mpime,ionode_set)
          INTEGER, INTENT(IN) :: mpime,ionode_set
          IF(mpime.EQ.ionode_set) THEN
            ionode = .TRUE.
          ELSE
            ionode = .FALSE.
          END IF
          ionode_id = ionode_set
          first = .FALSE.
          RETURN
        END SUBROUTINE io_global_start


        SUBROUTINE io_global_getionode( ionode_out, ionode_id_out )
          LOGICAL, INTENT(OUT) :: ionode_out
          INTEGER, INTENT(OUT) :: ionode_id_out
          IF( first ) &
            CALL error( ' get_ionode ', ' ionode not yet defined ', 1 )
          ionode_out = ionode
          ionode_id_out = ionode_id
          RETURN
        END SUBROUTINE io_global_getionode
       
    
      END MODULE io_global
