!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

       SUBROUTINE HANGUP
       IMPLICIT NONE
#if defined __MPI
       include 'mpif.h'
       INTEGER IERR
       CALL MPI_FINALIZE(IERR)
#endif
       RETURN
       END SUBROUTINE hangup
