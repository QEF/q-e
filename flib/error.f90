!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine errore (routin, messag, ierr)
  !----------------------------------------------------------------------
  !
  !    This is a simple routine which writes an error message to
  !    output. If ierr = 0 it does nothing, if ierr < 0 it
  !    writes the message but does not stop, if ierr > 0 stops.
  !
  !    **** Important note for parallel execution ***
  !    in parallel execution unit 6 is written only by the first node;
  !    all other nodes have unit 6 redirected to nothing (/dev/null).
  !    As a consequence an error not occurring on the first node
  !    will be invisible. For t3e and origin machines, this problem
  !    is solved by writing an error message to unit * instead of 6.
  !    For ibm sp machines, we write to the standard error, unit 0
  !    (this will appear in the error files produced by loadleveler).
  !
  USE io_global,  ONLY : stdout
  USE io_files, ONLY: crashunit, crash_file
  USE parallel_include
  !
  implicit none
  !
  character (len=*) :: routin, messag
    ! the name of the calling routine
    ! the output message
  integer :: ierr
    ! the error flag
  integer :: mpime
    ! the task id
  !
  !
  if (ierr == 0) return

  WRITE( * , * ) ' '
  WRITE( * , '(1x,78("%"))' )
  WRITE( * , '(5x,"from ",a," : error #",i10)' ) routin, ierr
  WRITE( * , '(5x,a)' ) messag
  WRITE( * , '(1x,78("%"))' )

#ifdef __PARA
#ifdef __AIX
  write (0, * ) ' '
  write (0, '(1x,78("%"))')
  write (0, '(5x,"from ",a," : error #",i10)') routin, ierr
  write (0, '(5x,a)') messag
  write (0, '(1x,78("%"))')
#endif
#endif

  if (ierr.gt.0) then

     WRITE( stdout , '("     stopping ...")')

#ifdef FLUSH
     call flush (6)
#endif

      mpime = 0

#ifdef __PARA
#ifdef __MPI
      CALL mpi_comm_rank( mpi_comm_world, mpime, ierr )
#endif
#endif

      !  .. write the message to a file and close it before exiting
      !  .. this will prevent loss of information on systems that
      !  .. do not flush the open streams
      !  .. added by C.C.

      OPEN( UNIT = crashunit, FILE = crash_file, POSITION = 'append', STATUS = 'unknown' )
      write (crashunit, * ) ' '
      write (crashunit, '(1x,78("%"))')
      write (crashunit, '(5x,"task #",i10)') mpime
      write (crashunit, '(5x,"from ",a," : error #",i10)') routin, ierr
      write (crashunit, '(5x,a)') messag
      write (crashunit, '(1x,78("%"))')
      CLOSE( UNIT = crashunit )

#ifdef __PARA
#ifdef __MPI
     CALL mpi_finalize( ierr )  !  try to exit in a smooth way
     IF ( ierr /= 0 ) THEN
       CALL MPI_ABORT( MPI_COMM_WORLD, ierr )
     END IF
#endif
#endif

     stop 2

  else

     WRITE( stdout, * ) ' '

  endif

end subroutine errore


!----------------------------------------------------------------------
subroutine infomsg (routin, messag, info)
  !----------------------------------------------------------------------
  !
  !    This is a simple routine which writes an info message 
  !    from a given routine to output. 
  !
  USE io_global,  ONLY : stdout, ionode
  USE parallel_include
  !
  implicit none
  !
  character (len=*) :: routin, messag
    ! the name of the calling routine
    ! the output message
  integer :: info
    ! the info code
  !
  IF( ionode ) THEN
    WRITE( stdout , '(5x,"from ",a," : info #",i10)' ) routin, info
    WRITE( stdout , '(5x,a)' ) messag
  END IF

  return

end subroutine infomsg
