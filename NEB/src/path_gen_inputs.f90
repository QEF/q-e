!
! Copyright (C) 2011-2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
SUBROUTINE path_gen_inputs (parse_file_name, engine_prefix, nimage)
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: parse_file_name
  CHARACTER(len=*), INTENT(in) :: engine_prefix
  INTEGER, INTENT(out) :: nimage
  !
  CHARACTER(LEN=80)  :: neb_file
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  INTEGER :: iimage
  INTEGER :: parse_unit, neb_unit
  !
  !
  OPEN( newunit=parse_unit, file=parse_file_name, status='old')
  iimage = 0
  neb_file = 'neb.dat'
  OPEN( newunit=neb_unit, file=neb_file, status='unknown')
  CALL parse_and_split (parse_unit, neb_unit, iimage, nimage)
  CLOSE( neb_unit )
  !
  image_loop: DO iimage = 1, nimage
     !
     REWIND (parse_unit)
     neb_file = trim(engine_prefix) // trim(int_to_char(iimage)) //".in"
     OPEN( newunit=neb_unit, file=neb_file, status='unknown' )
     CALL parse_and_split (parse_unit, neb_unit, iimage, nimage)
     CLOSE (neb_unit)
     !
  END DO image_loop
  !
  CLOSE (parse_unit)
  !
END SUBROUTINE path_gen_inputs
!
SUBROUTINE parse_and_split ( parse_unit, neb_unit, iimage, nimage )
  !
  ! Poor man parsing of NEB input file
  ! Data is read from "parse_unit", written to "neb_unit", for image "iimage"
  ! On output, "nimage" is the number of images read from "parse_unit"
  ! Must be called first with iimage = 0, later with iimage = 1,...,nimage
  ! parse_unit  must connected to the NEB input file
  ! neb_unit    must be connected to the output file, either neb input
  !             for iimage=0 or engine input for iimage > 0
  ! iimage      neb data is written if iimage = 0
  !             data for image iimage=1,..., nimage is written otherwise
  !  
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)  :: parse_unit, neb_unit, iimage
  INTEGER, INTENT(out) :: nimage
  !
  CHARACTER(len=512) :: dummy
  INTEGER :: status
  LOGICAL :: skip_line
  !
  nimage= 0
  status=-1
  !
  read_loop: DO
     !
     READ(parse_unit,fmt='(A512)',END=10) dummy
     IF( skip_line(dummy) ) CYCLE
     IF( trim(adjustl(dummy)) == "END" ) THEN
        EXIT read_loop
     ELSE IF( trim(adjustl(dummy)) == "BEGIN" ) THEN
        status=0
     ELSE IF( trim(adjustl(dummy)) == "BEGIN_PATH_INPUT"  ) THEN
        status=1
        ! ---------------------------------------------------
        ! NEB input data
        ! ---------------------------------------------------
     ELSE IF( trim(adjustl(dummy)) == "END_PATH_INPUT"    ) THEN
        status=0
     ELSE IF( trim(adjustl(dummy)) == "BEGIN_ENGINE_INPUT") THEN
        status=2
        ! ---------------------------------------------------
        ! Engine generic input data
        ! ---------------------------------------------------
     ELSE IF( trim(adjustl(dummy)) == "END_ENGINE_INPUT"  ) THEN
        status=0
     ELSE IF( trim(adjustl(dummy)) == "BEGIN_POSITIONS"   ) THEN
        status=3
        ! ---------------------------------------------------
        ! Engine image-specific data
        ! ---------------------------------------------------
     ELSE IF( trim(adjustl(dummy)) == "END_POSITIONS"      ) THEN
        status=2
     ELSE IF( trim(adjustl(dummy)) == "FIRST_IMAGE"       .OR. &
          trim(adjustl(dummy)) == "INTERMEDIATE_IMAGE" .OR. &
          trim(adjustl(dummy)) == "LAST_IMAGE"        ) THEN
        !
        nimage = nimage + 1
        !
     ELSE IF( status == 1 .AND. iimage == 0 ) THEN
        !
        WRITE(neb_unit,'(A)') trim(dummy)
        !
     ELSE IF( status == 2 .AND. iimage > 0 ) THEN
        !
        WRITE(neb_unit,'(A)') trim(dummy)
        !
     ELSE IF( status == 3 .AND. iimage > 0 .AND. iimage == nimage ) THEN
        !
        WRITE(neb_unit,'(A)') trim(dummy)
        !
     END IF
     !
  END DO read_loop
  !
  10 IF ( status /= 0 ) CALL infomsg('parse_and_split','parsing error?')
  !
END SUBROUTINE parse_and_split
!
LOGICAL FUNCTION skip_line(dummy)
  !
  ! True if line is empty or starts with # or !
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), intent(in) :: dummy 
  INTEGER :: i
  !
  IF (len_trim(dummy) < 1) THEN
    ! empty line
    skip_line = .TRUE.
  ELSE
    DO i=1,len_trim(dummy)
       IF ( dummy(i:i) /= ' ' ) EXIT
    END DO 
    ! Comment line (starts with ! or #)
    skip_line = dummy(i:i) .EQ. '!' .OR. dummy(i:i) .EQ. '#'
  END IF
  !
END FUNCTION skip_line
