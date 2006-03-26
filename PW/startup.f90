!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE startup( nd_nmbr, code, version )
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine initializes MPI
  !
  ! ... Processes are organized in NIMAGE images each dealing with a subset of
  ! ... images used to discretize the "path" (this only in "path" optimizations)
  ! ... Within each image processes are organized in NPOOL pools each dealing 
  ! ... with a subset of kpoints.
  ! ... Within each pool R & G space distribution is performed.
  ! ... NPROC is read from command line or can be set with the appropriate
  ! ... environment variable ( for example use 'setenv MP_PROCS 8' on IBM SP
  ! ... machine to run on NPROC=8 processors ); NIMAGE and NPOOL are read from 
  ! ... command line.
  ! ... NPOOL must be a whole divisor of NPROC
  !
  ! ... An example without any environment variable set is the following:
  !
  ! ... T3E :
  ! ...      mpprun -n 16 pw.x -npool 8 < input
  !
  ! ... IBM SP :
  ! ...      poe pw.x -procs 16 -npool 8 < input
  !
  ! ... ORIGIN /PC clusters using "mpirun" :
  ! ...      mpirun -np 16 pw.x -npool 8 < input
  !
  ! ... COMPAQ :
  ! ...      prun -n 16 sh -c 'pw.x -npool 8 < input'
  !
  ! ... PC clusters using "mpiexec" :
  ! ...      mpiexec -n 16 pw.x -npool 8 < input 
  ! 
  ! ... In this example you will use 16 processors divided into 8 pools
  ! ... of 2 processors each (in this case you must have at least 8 k-points)
  !
  ! ... The following two modules hold global information about processors
  ! ... number, IDs and communicators
  !
  USE io_global,  ONLY : stdout, io_global_start, meta_ionode, meta_ionode_id
  USE mp_global,  ONLY : nproc, nproc_image, nimage, mpime, me_image, &
                         my_image_id, root_image, npool, nproc_pool
  USE mp_global,  ONLY : mp_global_start, init_pool
  USE mp,         ONLY : mp_start, mp_env, mp_barrier, mp_bcast
  USE para_const, ONLY : maxproc
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=3)  :: nd_nmbr
  CHARACTER (LEN=6)  :: version
  CHARACTER (LEN=9)  :: code, cdate, ctime
  CHARACTER (LEN=80) :: np
  INTEGER            :: gid, node_number
  INTEGER            :: ierr = 0, ilen, nargs, iiarg
  INTEGER, EXTERNAL  :: iargc
  !
  !
#if defined (__PARA)
  !
  ! ... parallel case setup :  MPI environment is initialized
  !  
#  if defined (__T3E)
  !
  ! ... set streambuffers on
  !
  CALL set_d_stream( 1 )
  !
#  endif
  !
  CALL mp_start()
  !
  CALL mp_env( nproc, mpime, gid )
  !
  ! ... Set the I/O node
  !
  CALL io_global_start( mpime, 0 )
  !
  ! ... Set global coordinate for this processor
  !
  CALL mp_global_start( 0, mpime, gid, nproc )  
  !
  IF ( meta_ionode ) THEN
     !
     ! ... How many pools ?
     !
     npool = 1
     nargs = iargc() 
     !
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, np )
        !
        IF ( TRIM( np ) == '-npool' .OR. TRIM( np ) == '-npools' ) THEN
          !
          CALL getarg( ( iiarg + 1 ), np )  
          READ( np, * ) npool  
          !
        END IF
        !
     END DO
     !
     npool = MAX( npool, 1 )
     npool = MIN( npool, nproc )
     !
     ! ... How many parallel images ?
     !
     nimage = 1 
     nargs  = iargc() 
     !
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, np )
        !
        IF ( TRIM( np ) == '-nimage' .OR. TRIM( np ) == '-nimages' ) THEN
          !
          CALL getarg( ( iiarg + 1 ), np )  
          READ( np, * ) nimage 
          !
        END IF
        !
     END DO
     !
     nimage = MAX( nimage, 1 )
     nimage = MIN( nimage, nproc )
     !          
  END IF
  !
  CALL mp_barrier() 
  !
  ! ... transmit npool and nimage
  !
  CALL mp_bcast( npool,  meta_ionode_id )
  CALL mp_bcast( nimage, meta_ionode_id )
  !
  IF ( nproc > maxproc ) &
     CALL errore( 'startup', ' too many processors', nproc )
  !
  ! ... all pools are initialized here
  !
  CALL init_pool()
  !
  ! ... set the processor label for files ( remember that 
  ! ... me_image = 0 : ( nproc_image - 1 ) )
  !
  nd_nmbr = '   '
  !
  node_number = ( me_image + 1 )
  !
  IF ( nproc_image < 10 ) THEN
     !
     WRITE( nd_nmbr(1:1), '(I1)' ) node_number
     !
  ELSE IF ( nproc_image < 100 ) THEN
     !
     IF ( node_number < 10 ) THEN
        !
        nd_nmbr = '0'
        !
        WRITE( nd_nmbr(2:2), '(I1)' ) node_number
        !
     ELSE
        !
        WRITE( nd_nmbr(1:2), '(I2)' ) node_number
        !
     END IF
     !
  ELSE
     !
     IF ( node_number < 10 ) THEN
        !
        nd_nmbr = '00'
        !     
        WRITE( nd_nmbr(3:3), '(I1)' ) node_number
        !
     ELSE IF ( node_number < 100 ) THEN
        !
        nd_nmbr = '0'
        !
        WRITE( nd_nmbr(2:3), '(I2)' ) node_number
        !
     ELSE
        !
        WRITE( nd_nmbr, '(I3)' ) node_number
        !
     END IF
     !
  END IF    
  !
  ! ... stdout is printed only by the root_image ( set in init_pool() )
  !
#  if defined (DEBUG)
  !
  IF ( me_image /= root_image ) &
     OPEN( UNIT = stdout, FILE = './out_'//nd_nmbr, STATUS = 'UNKNOWN' )
  !   
#  else
  !
  IF ( me_image /= root_image ) &
     OPEN( UNIT = stdout, FILE = '/dev/null', STATUS = 'UNKNOWN' )
  !   
#  endif
  !
  ! ... information printout
  !  
  IF ( meta_ionode ) THEN
     !
     CALL date_and_tim( cdate, ctime )
     !
     WRITE( stdout, '(/5X,"Program ",A9," v.",A6," starts ...",&
                     &/5X,"Today is ",A9," at ",A9)' ) &
         code, version, cdate, ctime
     !
     WRITE( stdout, '(/5X,"Parallel version (MPI)",/)' )
     !
     WRITE( stdout, '(5X,"Number of processors in use:    ",I4)' ) nproc
     !
     IF ( nimage > 1 ) &
        WRITE( stdout, &
               '(5X,"path-images division:  nimage    = ",I4)' ) nimage
     IF ( npool > 1 ) &
        WRITE( stdout, &
               '(5X,"K-points division:     npool     = ",I4)' ) npool
     IF ( nproc_pool > 1 ) &
        WRITE( stdout, &
               '(5X,"R & G space division:  proc/pool = ",I4)' ) nproc_pool
     !
  END IF   
  !
#else
  !
  ! ... serial case setup :  only information printout
  !
  nd_nmbr = '   '
  !
  CALL date_and_tim( cdate, ctime )
  !
  WRITE( stdout, '(/5X,"Program ",A9," v.",A6," starts ...",&
                  &/5X,"Today is ",A9," at ",A9)' ) code, version, cdate, ctime
  !
#endif
  !
  RETURN
  !     
END SUBROUTINE startup
