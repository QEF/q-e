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
  use parameters
  implicit none
#ifdef __PARA
  include 'mpif.h'
#endif
  character (len=*) :: routin, messag
  ! the name of the calling routine
  ! the output message

  integer :: ierr
  ! the error flag
  if (ierr.eq.0) return
  write (6, * ) ' '
  write (6, '(1x,78("%"))')
  write ( * , '(5x,"from ",a," : error #",i10)') routin, ierr
  write ( * , '(5x,a)') messag
  write (6, '(1x,78("%"))')
#ifdef __PARA
#ifdef AIX
  write (0, * ) ' '
  write (0, '(1x,78("%"))')
  write (0, '(5x,"from ",a," : error #",i10)') routin, ierr
  write (0, '(5x,a)') messag
  write (0, '(1x,78("%"))')
#endif
#endif
  if (ierr.gt.0) then
     write ( * , '("     stopping ...")')
#ifdef FLUSH
     call flush (6)
#endif
#ifdef __PARA
     call mpi_abort (MPI_COMM_WORLD, ierr, ierr)
#endif
     stop 2
  else
     write (6, * ) ' '
     return
  endif
end subroutine errore

