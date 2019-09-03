  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  SUBROUTINE system_mem_usage(valueRSS)
  !----------------------------------------------------------------------
  ! 
  !! Report the memory usage ( VIRT and REAL ) from the current PID 
  !! process ( so it will be master only in case of MPI ).
  !! Memory is reported from the /proc/PID_NUMBER/status file
  ! 
  ! ---------------------------------------------------------------------
#ifdef __INTEL_COMPILER
  USE ifport !if on intel compiler
#endif  
  USE io_global,   ONLY : stdout
  USE io_epw,      ONLY : iunimem
  ! 
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(inout) :: valueRSS(2)
  !! Contains the value of the memory in kB
  ! 
  CHARACTER(LEN = 200) :: filename = ' '
  !! Name of the file
  CHARACTER(LEN = 80) :: line
  !! Line in the file
  CHARACTER(LEN = 8) :: pid_char = ' '
#if defined(__PGI) || defined(__CRAY) || defined(__XLF)
  INTEGER, EXTERNAL :: getpid
  !! PID of the process
#endif
  INTEGER :: pid
  !! PID of the process
  LOGICAL :: ifxst
  !! Does the file exists
  ! 
  valueRSS = -1    ! return negative number if not found
  ! 
  ! Get process ID
  ! 
  pid = getpid()
  WRITE(pid_char, '(I8)') pid
  filename = '/proc/' // TRIM(ADJUSTL(pid_char)) // '/status'
  ! 
  ! Read system file
  !
  INQUIRE(FILE = filename, EXIST = ifxst)
  IF (.NOT. ifxst) THEN
    WRITE(stdout, '(a)') 'System file does not exist'
    RETURN
  ENDIF
  ! 
  OPEN(UNIT = iunimem, FILE = filename, ACTION = 'read')
  ! 
  DO
    READ(iunimem, '(a)', END = 120) line
    ! Peak virtual memory usage
    IF (line(1:7) == 'VmPeak:') THEN
      READ(line(8:), *) valueRSS(1)
    ENDIF
    ! Peak resident set size
    IF (line(1:6) == 'VmHWM:') THEN
      READ(line(7:), *) valueRSS(2)
      CLOSE(UNIT = iunimem, STATUS = 'keep')
      EXIT
    ENDIF
  ENDDO
  120 CONTINUE
  ! 
  !--------------------------------------------------------------------------
  END SUBROUTINE system_mem_usage
  !--------------------------------------------------------------------------
