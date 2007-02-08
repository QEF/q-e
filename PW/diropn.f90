!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine diropn (unit, extension, recl, exst)
  !-----------------------------------------------------------------------
  !
  !     this routine opens a file named "prefix"."extension" in tmp_dir 
  !     for direct I/O access
  !     If appropriate, the node number is added to the file name
  !
#if defined(__ALPHA)
#  define DIRECT_IO_FACTOR 2
#else
#  define DIRECT_IO_FACTOR 8 
#endif
  !
  ! the  record length in direct-access I/O is given by the number of
  ! real*8 words times DIRECT_IO_FACTOR (may depend on the compiler)
  !
  USE kinds
  use io_files, only: prefix, tmp_dir, nd_nmbr
  use mp_global, only: mpime
  implicit none
  !
  !    first the input variables
  !
  character(len=*) :: extension
  ! input: name of the file to open
  integer :: unit, recl
  ! input: unit of the file to open
  ! input: length of the records
  logical :: exst
  ! output: if true the file exists
  !
  !    local variables
  !
  character(len=256) :: tempfile, filename
  ! complete file name
  integer :: ios, unf_recl
  ! used to check I/O operations
  ! length of the record
  logical :: opnd
  ! if true the file is already opened
  !
  if (unit < 0) call errore ('diropn', 'wrong unit', 1)
  !
  !    we first check that the file is not already openend
  !
  ios = 0
  inquire (unit = unit, opened = opnd)
  if (opnd) call errore ('diropn', "can't open a connected unit", abs(unit))
  !
  !      then we check the filename extension
  !
  if (extension == ' ') call errore ('diropn','filename extension not given',2)
  filename = trim(prefix) // "." // trim(extension)
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
  inquire (file = tempfile, exist = exst)
  !
  !      the unit for record length is unfortunately machine-dependent
  !
  unf_recl = DIRECT_IO_FACTOR * recl
  if (unf_recl <= 0) call errore ('diropn', 'wrong record length', 3)
  !
  open (unit, file = tempfile, iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)

  if (ios /= 0) call errore ('diropn', 'error opening '//filename, unit)
  return
end subroutine diropn

