! FOR GWW
!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Modified by P. Umari
! Modified by G. Stenuit
!
!-----------------------------------------------------------------------
subroutine dirdel ( extension)
  !-----------------------------------------------------------------------
! #ifdef __GWW
  !
  !     this routine opens a file named "prefix"."extension" in tmp_dir
  !     for direct I/O access
  !     If appropriate, the node number is added to the file name
  !
  USE kinds
  use io_files, only: prefix, tmp_dir, nd_nmbr, find_free_unit, diropn
  use mp_global, only: mpime
  use io_global, only: stdout
  implicit none

  !
  !    first the input variables
  !
  character(len=*) :: extension
  ! input: name of the file to open
  integer :: unit
  logical :: exst
  ! output: if true the file exists
  !
  !    local variables
  !
  character(len=256) :: tempfile, filename
  ! complete file name
  character(len=80) :: assstr
  integer :: ios, unf_recl, ierr,iunit
  ! used to check I/O operations
  ! length of the record
  ! error code
  ! if true the file is already opened



  if (extension == ' ') call errore ('diropn','filename extension not given',2)
  filename = trim(prefix) // "." // trim(extension)
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr

  INQUIRE( FILE = tempfile, EXIST = exst )
    !
  IF ( exst ) THEN
       !
       iunit = find_free_unit()

       OPEN(  UNIT = iunit, FILE = tempfile , STATUS = 'OLD' )
       CLOSE( UNIT = iunit, STATUS = 'DELETE' )
       !

       WRITE( UNIT = stdout, FMT = '(/,5X,"WARNING: ",A, &
            & " file was present; old file deleted")' ) tempfile

    ENDIF


! #endif __GWW
  return
end subroutine dirdel
