!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE blacs
!=----------------------------------------------------------------------------=!

        USE kinds

        IMPLICIT NONE
        PRIVATE
        SAVE

        PUBLIC :: start_blacs, stop_blacs

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

        SUBROUTINE start_blacs(iam_out, nprocs_out)

          INTEGER, INTENT(OUT), OPTIONAL :: iam_out, nprocs_out
          INTEGER :: iam, nprocs

#if defined __SCALAPACK

          CALL BLACS_PINFO( iam, nprocs )
          IF(PRESENT(iam_out)) THEN
            iam_out = iam
          END IF
          IF(PRESENT(nprocs_out)) THEN
            nprocs_out = nprocs
          END IF
#else

          iam_out = 0
          nprocs_out = 0

#endif

          RETURN
        END SUBROUTINE

!=----------------------------------------------------------------------------=!

        SUBROUTINE stop_blacs()

          INTEGER :: cont = 1  ! continue with the message passing after exiting from blacs

#if defined __SCALAPACK
          CALL BLACS_EXIT( cont )
#endif

          RETURN
        END SUBROUTINE

!=----------------------------------------------------------------------------=!
   END MODULE blacs
!=----------------------------------------------------------------------------=!
