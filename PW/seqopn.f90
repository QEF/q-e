!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine seqopn (unit, filename, formatt, exst)
  !-----------------------------------------------------------------------
  !
  !     this routine opens a file in tmp_dir for sequential I/O access
  !     If appropriate, the node number is added to the file name
  !
  use parameters
  use io_files
  implicit none
  !
  !    first the dummy variables
  !
  character(len=*) :: filename, formatt
  ! input: name of the file to connect
  ! input: 'formatted' or 'unformatted'
  integer :: unit
  ! input: unit to connect
  logical :: exst
  ! output: true if the file already exist
  !
  !    here the local variables
  !
  character(len=256) :: tempfile
  ! complete file name
  integer :: ios
  ! integer variable to test I/O status
  logical :: opnd
  ! true if the file is already opened


  if (unit < 1) call errore ('seqopn', 'wrong unit', 1)
  !
  !    test if the file is already opened
  !
  ios = 0
  inquire (unit = unit, opened = opnd)
  if (opnd) call errore ('seqopn', "can't open a connected unit", &
       abs (unit) )
  !
  !      then we check the filename
  !

  if (filename.eq.' ') call errore ('seqopn', 'filename not given', 2)
  if ( trim(nd_nmbr) == '1' .or. trim(nd_nmbr) == '01'.or. &
       trim(nd_nmbr) == '001') then
     !
     ! do not add processor number to files opened by processor 1
     ! in parallel execution: if only the first processor writes,
     ! we do not want the filename to be dependent on the number
     ! of processors
     !
     tempfile = trim(tmp_dir) // trim(filename)
  else
     tempfile = trim(tmp_dir) // trim(filename) // nd_nmbr
  end if
  inquire (file = tempfile, exist = exst)
  !
  !    Open the file
  !

  open (unit = unit, file = tempfile, form = formatt, status = &
       'unknown', iostat = ios)

  if (ios /= 0) call errore ('seqopn', 'error opening '//filename, &
       unit)
  return
end subroutine seqopn
