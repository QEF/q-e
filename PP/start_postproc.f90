!
! Copyright (C) 2001 PWSCF group
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
  !  Reads input data from filename if present, from standard input otherw
  !  On parallel machines filename is needed
  !
#include "machine.h"
  implicit none  
  integer :: nargs  
#ifdef T3D
  integer :: ipxfargc, ierr, ilen  
#else
  integer :: iargc  
#endif
  logical :: exst  


  character :: filin * 14, nodenumber * 3, version * 12  
  version = 'POSTPROC1.2-'  
  filin = ' '  
  nodenumber = '   '  
  !
  !  Read the number of arguments of the command
  !
#ifdef T3D
  nargs = ipxfargc ()  
#else
  nargs = iargc ()  
#endif
  if (nargs.gt.3) call error ('postproc', 'wrong no. of arguments', &
       nargs)
  !
  if (nargs.eq.1.or.nargs.eq.3) then  
     !
     !  Read the input file name (if any). It must be the last argument
     !
#ifdef T3D
     call pxfgetarg (nargs, filin, ilen, ierr)  
#else
     call getarg (nargs, filin)  
#endif
     if (filin.ne.' ') then  
        inquire (file = filin, exist = exst)  
        if (.not.exst) call error ('postproc', 'file '//filin//' not found', 1)
        open (unit = 5, form = 'formatted', status = 'old', file = filin)
     endif
  endif
#ifdef PARA

  if (filin.eq.' ') call error ('postproc','Usage: postproc [-npool N] filename',2)

  call startup (nodenumber, version)  

  call init_pool  
#endif
  return  
end subroutine start_postproc
