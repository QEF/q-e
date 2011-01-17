!program GWW P.Umari

MODULE start_end
!this module contains routines to initialize the MPI environment

CONTAINS

  SUBROUTINE startup

  !
  USE io_global,  ONLY : stdout, io_global_start, ionode
  USE mp_global,  ONLY : nproc, mpime, mp_global_start
  USE mp,         ONLY : mp_start

  IMPLICIT NONE
  CHARACTER(5) :: name_proc
  INTEGER :: gid

#ifdef __PARA


  CALL mp_start( nproc, mpime, gid )
  CALL io_global_start( mpime, 0 )
  CALL mp_global_start( 0, mpime, gid, nproc )

  if(ionode) then
     write(stdout,*) 'MPI PARALLEL VERSION'
     write(stdout,*) 'Number of procs: ', nproc
  endif
  write(name_proc,'(5i1)') &
       & (mpime+1)/10000,mod(mpime+1,10000)/1000,mod(mpime+1,1000)/100,mod(mpime+1,100)/10,mod(mpime+1,10)
!ATTENZIONE
  OPEN( UNIT = stdout, FILE = './out_'//name_proc, STATUS = 'UNKNOWN' )
#else
  OPEN( UNIT = stdout, FILE = './out_00', STATUS = 'UNKNOWN' )
#endif

  return

  END SUBROUTINE startup

  SUBROUTINE stop_run
!this subroutine kills the MPI environment

    USE io_global,         ONLY : stdout, ionode
    USE mp_global,         ONLY : mp_global_end

    IMPLICIT NONE

#ifdef __PARA

    if(ionode) write(stdout,*) 'Stopping MPI environment'
    call mp_global_end()
#endif

    return
  END SUBROUTINE stop_run


END MODULE start_end
