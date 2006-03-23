!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      MODULE turbo

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

        LOGICAL :: TTURBO
        INTEGER :: NTURBO
        COMPLEX(DP), ALLOCATABLE :: turbo_states(:,:) 

        PUBLIC :: tturbo, nturbo, turbo_states, turbo_init, allocate_turbo
        PUBLIC :: deallocate_turbo

      CONTAINS

        SUBROUTINE turbo_init(tturbo_inp, nturbo_inp)
          USE io_global, ONLY: ionode
          USE io_global, ONLY: stdout
          LOGICAL, INTENT(IN) :: tturbo_inp
          INTEGER, INTENT(IN) :: nturbo_inp
          tturbo = tturbo_inp
          nturbo = nturbo_inp
          IF( ionode .AND. tturbo ) THEN
            WRITE( stdout,fmt='(/,3X,"TURBO setup, nturbo = ",I10)') nturbo
          END IF
          RETURN
        END SUBROUTINE turbo_init

        SUBROUTINE allocate_turbo( nnr )
          USE io_global, ONLY: ionode, stdout
          USE mp_global, ONLY: intra_image_comm
          USE mp, ONLY: mp_sum
          INTEGER :: nnr
          INTEGER :: ierr
          IF( ionode ) THEN
            WRITE( stdout,fmt='(/,3X,"TURBO: allocating ",I10," bytes ")') &
              16*nnr*nturbo
          END IF
          IF( .NOT. ALLOCATED( turbo_states ) ) THEN
            ALLOCATE( turbo_states( nnr, nturbo ), STAT = ierr)
            CALL mp_sum( ierr, intra_image_comm )
            IF( ierr /= 0 ) THEN 
              IF( ionode ) THEN
                WRITE( stdout,fmt='(3X,"TURBO: insufficient memory, turbo is switched off ")')
              END IF
              tturbo = .FALSE.
              nturbo = 0
            END IF
          END IF
          RETURN
        END SUBROUTINE allocate_turbo 

        SUBROUTINE deallocate_turbo
          INTEGER :: ierr
          IF( ALLOCATED(turbo_states) ) THEN
            DEALLOCATE(turbo_states, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_turbo ', ' deallocating turbo_states ', ierr)
          END IF
          RETURN
        END SUBROUTINE deallocate_turbo

      END MODULE turbo
