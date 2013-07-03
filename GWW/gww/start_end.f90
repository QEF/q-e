!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


MODULE start_end
!this module contains routines to initialize the MPI environment

#ifdef __OPENMP
  INTEGER, SAVE :: ntids
#endif

CONTAINS

  SUBROUTINE startup

  !
  USE io_global,  ONLY : stdout, ionode
  USE mp_global,  ONLY : nproc, mpime
  USE mp_global,  ONLY : mp_startup
  USE mp,         ONLY : mp_barrier

  IMPLICIT NONE
  CHARACTER(5) :: name_proc
  INTEGER :: gid

#ifdef __PARA


!  CALL mp_start( nproc, mpime, gid )
!  CALL mp_barrier()
!  CALL io_global_start( mpime, 0 )
!  CALL mp_global_start_new( 0, mpime, gid, nproc )

  CALL mp_startup()

  if(ionode) then
     write(stdout,*) 'MPI PARALLEL VERSION'
     write(stdout,*) 'Number of procs: ', nproc
  endif
!  write(name_proc,'(5i1)') &
!       & (mpime+1)/10000,mod(mpime+1,10000)/1000,mod(mpime+1,1000)/100,mod(mpime+1,100)/10,mod(mpime+1,10)
!ATTENZIONE
!  OPEN( UNIT = stdout, FILE = './out_'//name_proc, STATUS = 'UNKNOWN' )
#else
!  OPEN( UNIT = stdout, FILE = './out_00', STATUS = 'UNKNOWN' )
#endif

  return

  END SUBROUTINE startup

  SUBROUTINE stop_run
!this subroutine kills the MPI environment

    USE io_global,         ONLY : stdout, ionode
    USE mp,                ONLY : mp_barrier, mp_end

    IMPLICIT NONE

#ifdef __PARA

    if(ionode) write(stdout,*) 'Stopping MPI environment'
    call mp_barrier()
    call mp_end()
#endif

    return
  END SUBROUTINE stop_run


END MODULE start_end
