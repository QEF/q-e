!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine start_postproc (nodenumber)
  !
  !  Usage: [mpirun, mpprun, whatever] postproc [-npool N] [filename]
  !
  !  Reads input data from filename if present
  !  from standard input otherwise
  !
#include "machine.h"
  implicit none
  integer :: nargs, ierr, ilen
  integer, external :: iargc
  logical :: exst


  character :: filin * 80, nodenumber * 3, version * 12
  version = 'POSTPROC-121'
  filin = ' '
  nodenumber = '   '
  !
  !  Read the number of arguments of the command
  !
  nargs = iargc ()
  if (nargs.gt.3) call errore ('postproc', 'wrong no. of arguments', &
       nargs)
  !
  if (nargs.eq.1.or.nargs.eq.3) then
     !
     !  Read the input file name (if any). It must be the last argument
     !
     call getarg (nargs, filin)
     if (filin.ne.' ') then
        inquire (file = filin, exist = exst)
        if (.not.exst) call errore ('postproc', 'file '//filin//' not found', 1)
        open (unit = 5, form = 'formatted', status = 'old', file = filin)
     endif
  endif

  call startup (nodenumber, version)

#ifdef __PARA
  call init_pool
#endif

  return
end subroutine start_postproc
