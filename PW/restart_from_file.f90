!
! Copyright (C) 2001,2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine restart_from_file
  !-----------------------------------------------------------------------
  USE io_global,  ONLY : stdout
  USE io_files,  ONLY : iunres
  USE control_flags, ONLY: restart
  implicit none

  character :: where * 20  ! parameter indicating from where to restart
  logical :: exst
  !
  !     check if restart file is present
  !
  iunres = 1
  if (.not.restart) then
     ! WRITE( stdout, '(/5x,"RECOVER from restart file has been switched off on input")')
     call seqopn (iunres, 'restart', 'unformatted', exst)
     ! if (exst) WRITE( stdout,'(/5x,"Existing restart file has been removed")')
     close (unit = iunres, status = 'delete')
     return
  endif

  call seqopn (iunres, 'restart', 'unformatted', restart)
  if (.not.restart) then
     WRITE( stdout, '(/5x,"RECOVER from restart file failed: file not found")')
     close (unit = iunres, status = 'delete')
     return
  endif
  !
  WRITE( stdout, '(/5x,"read information from restart file")')
  read (iunres, err = 10, end = 10) where
  WRITE( stdout, '(5x,"Restarting in ",a)') where
  if (where.ne.'ELECTRONS'.and.where.ne.'IONS') then
     WRITE( stdout,*) where, '......?'
     call errore ('readin', ' wrong recover file ', 1)
  endif
  !
  !  close the file for later use
  !
  close (unit = iunres, status = 'keep')

  return

10 call errore ('readin', 'problems in reading recover file', 1)
end subroutine restart_from_file
