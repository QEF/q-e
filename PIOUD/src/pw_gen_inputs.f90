!
! Copyright (C) 2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Aadhityan A, Lorenzo Paulatto, Michele Casula, Tommaso Morresi
!
!----------------------------------------------------------------------------
SUBROUTINE pw_gen_inputs(parse_file_name,engine_prefix,nimage,root,comm)
  !
  USE mp, ONLY : mp_rank
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
  INTEGER :: parse_unit
  INTEGER, ALLOCATABLE :: unit_tmp(:)
  INTEGER :: unit_tmp_i
  CHARACTER(len=10) :: a_tmp
  INTEGER :: myrank
  INTEGER, EXTERNAL :: find_free_unit
  nimage=1
  !
  myrank =  mp_rank(comm)
  parse_unit = find_free_unit()
  OPEN(unit=parse_unit,file=trim(parse_file_name),status="old")
  !
  ! ------------------------------------------------
  ! ENGINE INPUT PART
  ! ------------------------------------------------
  
  ALLOCATE(unit_tmp(1:nimage))
  unit_tmp(:) = 0
  
  !WRITE(9999,*) "Reading images", nimage

!   DO i=1,nimage
!      unit_tmp(i) = find_free_unit()
!   ENDDO
  
  DO i=1,nimage
     IF(i>=1.and.i<10) THEN
        WRITE(a_tmp,'(i1)') i
     ELSEIF(i>=10.and.i<100) THEN
        WRITE(a_tmp,'(i2)') i
     ELSEIF(i>=100.and.i<1000) THEN
        WRITE(a_tmp,'(i3)')
     ENDIF 
     OPEN(newunit=unit_tmp_i,file=trim(engine_prefix)//trim(a_tmp)//".in")
     unit_tmp(i) = unit_tmp_i 
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

               IF(trim(ADJUSTL(dummy)) == "BEGIN_POSITIONS") THEN
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
END SUBROUTINE pw_gen_inputs
