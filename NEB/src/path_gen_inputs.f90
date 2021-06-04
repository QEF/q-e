!
! Copyright (C) 2011-2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
LOGICAL FUNCTION skip_line(dummy)
  !
  ! True if line is empty or starts with # or !
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), intent(in) :: dummy 
  CHARACTER(len=1) :: dum1
  !
  IF (len_trim(dummy) < 1) THEN
    ! empty line
    skip_line = .TRUE.
    RETURN
  END IF
  !
  dum1 = trim(adjustl(dummy))
  ! Comment line (starts with ! or #)
  skip_line = dum1 .EQ. "!" .OR. dum1 .EQ. "#"
  !
END FUNCTION skip_line
!
SUBROUTINE path_gen_inputs(parse_file_name,engine_prefix,nimage,root,comm)
  !
  USE mp, ONLY : mp_rank, mp_bcast
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: parse_file_name
  CHARACTER(len=*), INTENT(in) :: engine_prefix
  INTEGER, INTENT(out) :: nimage
  INTEGER, INTENT(in) :: root
  INTEGER, INTENT(in) :: comm
  !
  CHARACTER(len=512) :: dummy
  CHARACTER(LEN=80)  :: engine_file, tmp_file
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  INTEGER :: iimage, status
  INTEGER :: parse_unit, engine_unit, tmp_unit
  LOGICAL :: skip_line
  !
  IF( mp_rank(comm) == root) THEN
     !
     OPEN( newunit=parse_unit, file=trim(parse_file_name), status='old')
     nimage= 0
     status=-1
     !
     read_loop: DO
        !
        READ(parse_unit,fmt='(A512)',END=10) dummy
        IF ( skip_line(dummy) ) CYCLE
        IF( dummy == "END" ) THEN
           EXIT read_loop
        ELSE IF( trim(adjustl(dummy)) == "BEGIN" ) THEN
           status=0
        ELSE IF( trim(adjustl(dummy)) == "BEGIN_PATH_INPUT"  ) THEN
           status=1
           OPEN( newunit=tmp_unit, file='neb.dat', status='unknown' )
           ! ---------------------------------------------------
           ! NEB input data
           ! ---------------------------------------------------
        ELSE IF( trim(adjustl(dummy)) == "END_PATH_INPUT"    ) THEN
           status=0
           CLOSE( tmp_unit )
        ELSE IF( trim(adjustl(dummy)) == "BEGIN_ENGINE_INPUT") THEN
           status=2
           tmp_file = 'tmp.in' // trim(int_to_char(nimage))
           OPEN( newunit=tmp_unit, file=tmp_file, status='unknown' )
           ! ---------------------------------------------------
           ! Engine generic input data
           ! ---------------------------------------------------
        ELSE IF( trim(adjustl(dummy)) == "END_ENGINE_INPUT"  ) THEN
           status=0
           CLOSE( tmp_unit )
        ELSE IF( trim(adjustl(dummy)) == "BEGIN_POSITIONS"   ) THEN
           status=3
           ! ---------------------------------------------------
           ! Engine image-specific data here
           ! ---------------------------------------------------
        ELSE IF( trim(adjustl(dummy)) == "END_POSITIONS"      ) THEN
           status=2
           CLOSE (tmp_unit)
           tmp_file = 'tmp.in' // trim(int_to_char(nimage+1))
           OPEN( newunit=tmp_unit, file=tmp_file, status='unknown' )
        ELSE IF( trim(adjustl(dummy)) == "FIRST_IMAGE"       .OR. &
                trim(adjustl(dummy)) == "INTERMEDIATE_IMAGE" .OR. &
                trim(adjustl(dummy)) == "LAST_IMAGE"        ) THEN
           nimage = nimage + 1
           IF ( nimage > 1 ) CLOSE (tmp_unit)
           tmp_file = 'tmp.in' // trim(int_to_char(nimage))
           OPEN( newunit=tmp_unit, file=tmp_file, status='unknown' )
        ELSE IF( status > 0 ) THEN
           !
           WRITE(tmp_unit,*) trim(dummy)
           !
        END IF
        !
     END DO read_loop
     !
10   CONTINUE
     CLOSE (parse_unit)
     CLOSE (tmp_unit)
     !
  END IF
  !
  CALL mp_bcast(nimage, root, comm)
  !
  IF( mp_rank(comm) == root) THEN
     !
     DO iimage = 1, nimage
        !
        engine_file = trim(engine_prefix) // trim(int_to_char(iimage)) //".in"
        OPEN( newunit=engine_unit, file=engine_file, status='unknown' )
        !
        tmp_file = 'tmp.in0'
        OPEN( newunit=tmp_unit, file=tmp_file, status='old' )
        DO
           READ ( tmp_unit, fmt='(A512)', END=20) dummy
           WRITE( engine_unit, '(A)' ) trim(dummy)
        END DO
20      CLOSE (tmp_unit)
        !
        tmp_file = 'tmp.in' // trim(int_to_char(iimage))
        OPEN( newunit=tmp_unit, file=tmp_file, status='old' )
        DO
           READ ( tmp_unit, fmt='(A512)', END=30) dummy
           WRITE( engine_unit, '(A)' ) trim(dummy)
        END DO
30      CLOSE (tmp_unit) 
        !
        tmp_file = 'tmp.in' // trim(int_to_char(nimage+1))
        OPEN( newunit=tmp_unit, file=tmp_file, status='old' )
        DO
           READ ( tmp_unit, fmt='(A512)', END=40) dummy
           WRITE( engine_unit, '(A)' ) trim(dummy)
        END DO
40      CLOSE (tmp_unit)
        !
        CLOSE (engine_unit) 
        !
     END DO
     !
   END IF
   !
   END SUBROUTINE path_gen_inputs
