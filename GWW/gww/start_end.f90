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
    IMPLICIT NONE
    CHARACTER (len=10), PARAMETER :: code = 'GWW'
#if defined(__OPENMP)
    INTEGER, SAVE :: ntids
#endif

CONTAINS

  SUBROUTINE startup

  !
  USE io_global,  ONLY : stdout, ionode
  USE mp_world,   ONLY : nproc
  USE mp_global,  ONLY : mp_startup
  USE environment,           ONLY: environment_start
  
  IMPLICIT NONE

#if defined(__MPI)
  CALL mp_startup()
#endif
  
  CALL environment_start ( code )
  
#if defined(__MPI)
  if(ionode) then
     write(stdout,*) 'MPI PARALLEL VERSION'
     write(stdout,*) 'Number of procs: ', nproc
     write(stdout,*)  'GWL: Version 1.00'
  endif
#else 
   write(stdout,*)  'GWL: Version 1.00'
#endif
  return

  END SUBROUTINE startup

  SUBROUTINE stop_run
!this subroutine kills the MPI environment

    USE io_global,         ONLY : stdout, ionode
    USE mp_global,         ONLY : mp_global_end

    IMPLICIT NONE

#if defined(__MPI)

    if(ionode) write(stdout,*) 'Stopping MPI environment'
    call mp_global_end( )
#endif

    return
  END SUBROUTINE stop_run


END MODULE start_end
