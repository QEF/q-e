!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE io_global

        USE kinds
        IMPLICIT NONE
        PRIVATE
        SAVE

        PUBLIC :: io_global_start, ionode, ionode_id

        INTEGER :: ionode_id = 0
        LOGICAL :: ionode = .TRUE.
      
      CONTAINS

        SUBROUTINE io_global_start(mpime,ionode_set)
          INTEGER, INTENT(IN) :: mpime,ionode_set
          IF(mpime.EQ.ionode_set) THEN
            ionode = .TRUE.
          ELSE
            ionode = .FALSE.
          END IF
          ionode_id = ionode_set
          RETURN
        END SUBROUTINE io_global_start
    
      END MODULE io_global
