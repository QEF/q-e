  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  SUBROUTINE system_mem_usage (valueRSS)
  !----------------------------------------------------------------------
  !! Report the memory usage ( VIRT and REAL ) from the current PID 
  !! process ( so it will be master only in case of MPI ).
  !! Memory is reported from the /proc/PID_NUMBER/status file
  ! ---------------------------------------------------------------------
#ifdef __INTEL_COMPILER
  USE ifport !if on intel compiler
#endif  
  USE io_global,   ONLY : stdout
  USE io_epw,      ONLY : iunimem
  ! 
  implicit none
  ! 
  INTEGER, INTENT (inout) :: valueRSS(2)
  !! Contains the value of the memory in kB
  ! 
  CHARACTER(len=200):: filename=' '
  CHARACTER(len=80) :: line
  CHARACTER(len=8)  :: pid_char=' '
#if defined(__PGI)
  INTEGER, EXTERNAL :: getpid
#endif
  INTEGER :: pid
  LOGICAL :: ifxst

  valueRSS = -1    ! return negative number if not found

  !--- get process ID

  pid=getpid()
  WRITE(pid_char,'(I8)') pid
  filename='/proc/'//trim(adjustl(pid_char))//'/status'

  !--- read system file

  inquire (file=filename,exist=ifxst)
  if (.not.ifxst) then
    write (stdout,'(a)') 'System file does not exist'
    return
  endif

  open(unit=iunimem, file=filename, action='read')
  do
    read (iunimem,'(a)',end=120) line
    ! Peak virtual memory usage
    if (line(1:7).eq.'VmPeak:') then
       read (line(8:),*) valueRSS(1)
    endif
    ! Peak resident set size
    if (line(1:6).eq.'VmHWM:') then
       read (line(7:),*) valueRSS(2)
       CLOSE (unit = iunimem, status = 'keep')
       exit
    endif        
  enddo
  120 continue

  END SUBROUTINE system_mem_usage
  !--------------------------------------------------------------------------
