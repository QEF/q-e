!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine diropn (unit, filename, recl, exst)
  !-----------------------------------------------------------------------
  !
  !     this routine opens a file in tmp_dir for direct I/O access
  !     If appropriate, the node number is added to the file name
  !
#include "machine.h"
  use parameters
  use io
  use mp_global, only: mpime
  implicit none

  !
  !    first the input variables
  !
  character :: filename * ( * )
  ! input: name of the file to ope
  integer :: unit, recl
  ! input: unit of the file to ope
  ! input: length of the records
  logical :: exst
  ! output: if true the file exist
  !
  !    local variables
  !
  character(len=80) :: assstr, tempfile
  ! complete file name
  integer :: ios, unf_recl, ierr
  ! used to check I/O operations
  ! length of the record
  ! error code
  logical :: opnd
  ! if true the file is already opened


  if (unit.le.0) call errore ('diropn', 'wrong unit', 1)
  !
  !    we first check that the file is not already openend
  !
  ios = 0
  inquire (unit = unit, opened = opnd)
  if (opnd) call errore ('diropn', 'can"t open a connected unit', abs(unit))
  !
  !      then we check the filename
  !

  if (filename.eq.' ') call errore ('diropn', 'filename not given', 2)
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
  ! debug
  !write(200+mpime,*) trim(tmp_dir)
  !write(200+mpime,*) trim(filename)
  !write(200+mpime,*) nd_nmbr
  !write(200+mpime,*) tempfile
  !close(200+mpime)
  !return
  ! end debug
  inquire (file = tempfile, exist = exst)
  !
  !      the unit for record length is unfortunately machine-dependent
  !
  unf_recl = DIRECT_IO_FACTOR * recl
  if (unf_recl.le.0) call errore ('diropn', 'wrong record length', 3)
  !
  !     on T3E reduce the size of the buffer if it is too large
  !
#ifdef __T3E
  if (unf_recl.gt.5000000) then
     if (unit.lt.10) then
        write (assstr, '("assign -b 1 u:",i1)') unit
     else if(unit.lt.100) then
        write (assstr, '("assign -b 1 u:",i2)') unit
     else
        call errore ('diropn', 'unit too large', 1)
     endif
     call assign (assstr, ierr)
  endif
#endif

  open (unit, file = tempfile, iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)

  if (ios.ne.0) call errore ('diropn', 'error opening '//filename, unit)
  return
end subroutine diropn

