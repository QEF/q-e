!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE path_gen_inputs(parse_file_name,engine_prefix,nimage,root,comm)
  !
  USE mp_global, ONLY : mp_rank
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
  INTEGER :: i, j
  INTEGER :: parse_unit, neb_unit
  INTEGER, ALLOCATABLE :: unit_tmp(:)
  INTEGER :: unit_tmp_i
  CHARACTER(len=10) :: a_tmp
  INTEGER :: myrank
  INTEGER, EXTERNAL :: find_free_unit
  !
  myrank =  mp_rank(comm)
  parse_unit = find_free_unit()
  OPEN(unit=parse_unit,file=trim(parse_file_name),status="old")
  ! ---------------------------------------------------
  ! NEB INPUT PART
  ! ---------------------------------------------------  
  i=0
  nimage = 0
  neb_unit = find_free_unit()
  OPEN(unit=neb_unit,file='neb.dat',status="unknown")
  dummy=""
  DO WHILE (len_trim(dummy)<1)
     READ(parse_unit,fmt='(A512)',END=10) dummy
  ENDDO
  
  IF(trim(ADJUSTL(dummy)) == "BEGIN") THEN
     DO WHILE (trim(ADJUSTL(dummy)) /= "END")
        READ(parse_unit,*) dummy
        IF(trim(ADJUSTL(dummy)) == "BEGIN_PATH_INPUT") THEN
           
           READ(parse_unit,'(A512)') dummy
           
           DO WHILE (trim(ADJUSTL(dummy)) /= "END_PATH_INPUT")
              IF(myrank==root) WRITE(neb_unit,*) trim(ADJUSTL(dummy))
              READ(parse_unit,'(A512)') dummy
           ENDDO
        ENDIF
        IF(trim(ADJUSTL(dummy)) == "FIRST_IMAGE") THEN
           nimage = nimage + 1
        ENDIF
        IF(trim(ADJUSTL(dummy)) == "INTERMEDIATE_IMAGE") THEN
           nimage = nimage + 1
        ENDIF
        IF(trim(ADJUSTL(dummy)) == "LAST_IMAGE") THEN
           nimage=nimage+1
        ENDIF
     ENDDO
  ELSE
     CALL infomsg('path_gen_inputs','key word BEGIN missing')
  ENDIF
  CLOSE(neb_unit)
  !------------------------------------------------
  !
  !
  ! ------------------------------------------------
  ! ENGINE INPUT PART
  ! ------------------------------------------------
  
  ALLOCATE(unit_tmp(1:nimage))
  unit_tmp(:) = 0
  
  DO i=1,nimage
     unit_tmp(i) = find_free_unit()
  ENDDO
  
  DO i=1,nimage
     IF(i>=1.and.i<10) THEN
        WRITE(a_tmp,'(i1)') i
     ELSEIF(i>=10.and.i<100) THEN
        WRITE(a_tmp,'(i2)') i
     ELSEIF(i>=100.and.i<1000) THEN
        WRITE(a_tmp,'(i3)')
     ENDIF
     unit_tmp_i = unit_tmp(i)
     OPEN(unit=unit_tmp_i,file=trim(engine_prefix)//trim(a_tmp)//".in")
     REWIND(parse_unit)
     dummy=""
     DO WHILE (len_trim(dummy)<1)
        READ(parse_unit,fmt='(A512)',END=10) dummy
     ENDDO
     
     IF(trim(ADJUSTL(dummy)) == "BEGIN") THEN
        DO WHILE (trim(ADJUSTL(dummy)) /= "END")
           dummy=""
           DO WHILE (len_trim(dummy)<1)
              READ(parse_unit,fmt='(A512)',END=10) dummy
           ENDDO
           
           IF(trim(ADJUSTL(dummy)) == "BEGIN_ENGINE_INPUT") THEN
              dummy=""
              DO WHILE (len_trim(dummy)<1)
                 READ(parse_unit,fmt='(A512)',END=10) dummy
              ENDDO
              
              DO WHILE (trim(ADJUSTL(dummy)) /= "BEGIN_POSITIONS")
                 IF(myrank==root) WRITE(unit_tmp_i,'(A)') trim(dummy)
                 READ(parse_unit,'(A512)') dummy
              ENDDO
              IF(i==1) THEN
                 DO WHILE (trim(ADJUSTL(dummy)) /= "FIRST_IMAGE")
                    READ(parse_unit,'(A512)') dummy
                 ENDDO
                 IF(trim(ADJUSTL(dummy)) == "FIRST_IMAGE") THEN
                    READ(parse_unit,'(A512)') dummy
                    DO WHILE (trim(ADJUSTL(dummy)) /= "INTERMEDIATE_IMAGE" &
                        .and.(trim(ADJUSTL(dummy)) /= "LAST_IMAGE"))
                       IF(myrank==root) WRITE(unit_tmp_i,'(A)') trim(ADJUSTL(dummy))
                       READ(parse_unit,'(A512)') dummy
                    ENDDO
                    DO WHILE (trim(ADJUSTL(dummy)) /= "END_POSITIONS")
                       READ(parse_unit,'(A512)') dummy
                    ENDDO
                    READ(parse_unit,'(A512)') dummy
                    DO WHILE (trim(ADJUSTL(dummy)) /= "END_ENGINE_INPUT")
                       IF(myrank==root) WRITE(unit_tmp_i,'(A)') trim(ADJUSTL(dummy))
                       READ(parse_unit,'(A512)') dummy
                    ENDDO
                 ENDIF
              ENDIF
              !
              IF(i==nimage) THEN
                 DO WHILE (trim(ADJUSTL(dummy)) /= "LAST_IMAGE")
                    READ(parse_unit,'(A512)') dummy
                 ENDDO
                 IF(trim(ADJUSTL(dummy)) == "LAST_IMAGE") THEN
                    READ(parse_unit,'(A512)') dummy
                    DO WHILE (trim(ADJUSTL(dummy)) /= "END_POSITIONS")
                       IF(myrank==root) WRITE(unit_tmp_i,'(A)') trim(ADJUSTL(dummy))
                       READ(parse_unit,'(A512)') dummy
                    ENDDO
                    READ(parse_unit,'(A512)') dummy
                    DO WHILE (trim(ADJUSTL(dummy)) /= "END_ENGINE_INPUT")
                       IF(myrank==root) WRITE(unit_tmp_i,'(A)') trim(ADJUSTL(dummy))
                       READ(parse_unit,'(A512)') dummy
                    ENDDO
                 ENDIF
              ENDIF
              !
              IF(i/=nimage.and.i/=1) THEN
                 DO j=2,i
                    dummy=""
                    DO WHILE (trim(ADJUSTL(dummy)) /= "INTERMEDIATE_IMAGE")
                       READ(parse_unit,'(A512)') dummy
!!! write(0,*) i,j,trim(dummy)
                    ENDDO
                 ENDDO
                 IF(trim(ADJUSTL(dummy)) == "INTERMEDIATE_IMAGE") THEN
                    READ(parse_unit,'(A512)') dummy
                    DO WHILE ((trim(ADJUSTL(dummy)) /= "LAST_IMAGE") .and. &
                               trim(ADJUSTL(dummy)) /= "INTERMEDIATE_IMAGE")
                       IF(myrank==root) WRITE(unit_tmp_i,'(A)') trim(ADJUSTL(dummy))
                       READ(parse_unit,'(A512)') dummy
                    ENDDO
                    DO WHILE (trim(ADJUSTL(dummy)) /= "END_POSITIONS")
                       READ(parse_unit,'(A512)') dummy
                    ENDDO
                    READ(parse_unit,'(A512)') dummy
                    DO WHILE (trim(ADJUSTL(dummy)) /= "END_ENGINE_INPUT")
                       IF(myrank==root) WRITE(unit_tmp_i,'(A)') trim(ADJUSTL(dummy))
                       READ(parse_unit,'(A512)') dummy
                    ENDDO
                 ENDIF
              ENDIF
              !
           ENDIF
        ENDDO
     ENDIF
     
     CLOSE(unit_tmp_i)
  ENDDO
  
  DEALLOCATE(unit_tmp)
  
  CLOSE(parse_unit)
  !
10 CONTINUE
  !
END SUBROUTINE path_gen_inputs
