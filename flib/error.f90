!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Mon Nov 15 13:26:42 MET 1999
!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE error(a,b,n)

!  this routine prints an error and warning message and 
!  if necessary terminates the program.
!  INPUT: a, b, n
!    a   (character)   subroutine name
!    b   (character)   error message
!    n   (integer)     error code
!                      if n > 0 write the error message and terminate the execution
!                      if n = 0 do nothing
!                      if n < 0 print the error message and return
!  OUTPUT: none
!  ----------------------------------------------
!  END manual

      IMPLICIT NONE

#if defined __MPI
      include 'mpif.h'
#endif

! ... declare subroutine arguments
      CHARACTER(LEN=*)    :: a, b
      INTEGER, INTENT(IN) :: n

      INTEGER :: ip, nproc, mpime, ierr

! ... declare function

!  end of declarations
!  ----------------------------------------------

#if defined __MPI
      CALL mpi_comm_size(mpi_comm_world,nproc,ierr)
      CALL mpi_comm_rank(mpi_comm_world,mpime,ierr)
#else
      MPIME = 0
      NPROC = 1
#endif

! ... print the error message
!
      DO ip = 0, nproc-1
        IF( n > 0 ) THEN
          WRITE (6,100) mpime, a, b, n
          OPEN(UNIT=15, FILE='CRASH', POSITION='append', STATUS='unknown')
          WRITE (15,100) mpime, a, b, n
          CLOSE(UNIT=15)
        ELSE IF ( n < 0 ) THEN
          IF( mpime == 0 ) WRITE (6,200) a, b
        END IF
#if defined __MPI
        CALL mpi_barrier(mpi_comm_world,ierr)
#endif
      END DO

! ... terminate the program
!
      CALL cpflush  ! flush output streams

      IF( n > 0 ) THEN
#if defined __MPI
        CALL mpi_finalize(ierr)
        IF ( ierr/=0 ) THEN
          CALL mpi_abort(mpi_comm_world, ierr)
        END IF
#endif
      END IF

100   FORMAT (/,' *** from PE    : ',I5, &
              /,' *** in routine : ',A, &
              /,' *** error msg. : ',A, &
              /,' *** error code : ',I5, &
              /,' *** aborting ***', /)

200   FORMAT ('   Warning (', A, ') : ', A)

      IF( n > 0 ) THEN
        STOP 'CRASH'
      ELSE IF( n == 0 ) THEN
        WRITE(6,*) ' ERROR DEBUG ', a, b
      END IF
 
      RETURN

      END SUBROUTINE


      SUBROUTINE warning(a, b, n)
! ... declare subroutine arguments
        CHARACTER(LEN=*), INTENT(IN) :: a, b
        INTEGER, INTENT(IN) :: n
        WRITE (6,100) a, b, n
        RETURN
100     FORMAT (/,'*WARNING* from ',A,', ',A,/)
      END SUBROUTINE
