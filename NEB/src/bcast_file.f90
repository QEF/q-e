!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE bcast_file ( filin, root, comm, ios ) 
  !-----------------------------------------------------------------------------
  !! Check whether file "filin", that must be present on processor "root",
  !! is also visible to all processes in "comm"; make a local copy, if not.
  !! On output, ios contains the return status:
  !! ios=0 nothing done, ios=-1 successfully copied, ios=1 failure
  !---------------------------------------------------------------
  !
  USE mp, only: mp_rank, mp_size, mp_bcast, mp_sum
  !
  IMPLICIT NONE
  !
  CHARACTER (len=*), intent(in) :: filin
  INTEGER, intent(in) :: root
  INTEGER, intent(in) :: comm
  INTEGER, intent(out):: ios
  !
  character(len=512) :: line
  integer :: root_unit, localunit
  integer :: filesize, goodsize
  integer :: nlines, n
  logical :: ishere, tohere
  !
  ! true if the original file is on this processor
  ishere = ( mp_rank(comm) == root )
  INQUIRE( FILE=filin, SIZE=filesize )
  IF ( ishere ) goodsize = filesize
  CALL mp_bcast( goodsize, root, comm )
  ! check: can all images see a file with the same size?
  ios = abs(filesize - goodsize)
  CALL mp_sum( ios, comm )
  IF ( ios == 0 ) RETURN
  !
  ! true if the original file must be copied to this processor
  tohere = ( filesize /= goodsize )
  IF ( ishere ) THEN
     OPEN( NEWUNIT = root_unit, FILE = filin, STATUS = 'old', &
          FORM='formatted', iostat = ios )
  ELSE IF ( tohere ) THEN
     OPEN( NEWUNIT = localunit, FILE = filin, STATUS='unknown',&
          FORM='formatted', iostat = ios )
  END IF
  CALL mp_sum( ios, comm )
  IF ( ios > 0 ) RETURN
  !
  ! count lines: not smart but I haven't found a smarter way
  ! (and no, you cannot just use "END=", you end up with a deadlock)
  !
  nlines = 0 
  IF ( ishere ) THEN
     DO
        READ(root_unit,'(A512)', END=10) line
        nlines = nlines+1
     END DO
10   REWIND(root_unit)
  END IF
  CALL mp_bcast( nlines, root, comm )
  DO n = 1, nlines
     IF ( ishere ) READ(root_unit,'(A512)') line
     CALL mp_bcast( line, root, comm )
     IF ( tohere ) WRITE(localunit,'(A)') trim(line)
  END DO
  IF ( ishere ) CLOSE ( unit=root_unit )
  IF ( tohere ) CLOSE ( unit=localunit )
  ios = -1
  !
END SUBROUTINE bcast_file
