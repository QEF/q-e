!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE startup( nd_nmbr, code, version )
  !----------------------------------------------------------------------------
  !
  !  ... This subroutine initializes MPI
  !
  !  ... Processes are organized in NPOOL pools each dealing with a subset of
  !  ... kpoints. Within each pool R & G space distribution is performed.
  !  ... NPROC is read from command line or can be set with the appropriate
  !  ... environment variable ( for example use 'setenv MP_PROCS 8' on IBM SP
  !  ... machine to run on NPROC=8 processors ) and NPOOL is read from command
  !  ... line.
  !  ... NPOOL must be a whole divisor of NPROC
  !
  !  ... An example without any environment variable set is the following:
  !
  !  ... T3E :
  !  ...      mpprun -n 16 pw.x -npool 8 < input
  !
  !  ... IBM SP :
  !  ...      poe pw.x -procs 16 -npool 8 < input
  !
  !  ... ORIGIN /PC clusters using "mpirun" :
  !  ...      mpirun -np 16 pw.x -npool 8 < input
  !
  !  ... COMPAQ :
  !  ...      prun -n 16 sh -c 'pw.x -npool 8 < input'
  !
  !  ... PC clusters using "mpiexec" :
  !  ...      mpiexec -n 16 pw.x -npool 8 < input 
  ! 
  !  ... In this example you will use 16 processors divided into 8 pools
  !  ... of 2 processors each (in this case you must have at least 8 k-points)
  !
  ! ... The following two modules hold global information about processors
  ! ... number, IDs and communicators
  !
  USE io_global,  ONLY :  stdout, io_global_start, ionode, ionode_id
  USE mp_global,  ONLY :  nproc, nimage, mpime, me_image, root, root_image
  USE mp_global,  ONLY :  mp_global_start
  USE mp,         ONLY :  mp_start, mp_env, mp_barrier, mp_bcast
  USE para_const, ONLY :  maxproc
  USE para,       ONLY :  me, npool, nprocp 
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=3)  :: nd_nmbr
  CHARACTER (LEN=6)  :: version
  CHARACTER (LEN=9)  :: code, cdate, ctime
  CHARACTER (LEN=80) :: np
  INTEGER            :: gid
  EXTERNAL              date_and_tim
  INTEGER            :: ierr = 0, ilen, iargc, nargs, iiarg
  !
  !
#if defined (__PARA)
  !
  ! ... prallel case setup :  MPI environment is initialized
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
  CALL mp_env( nproc, me, gid )
  !
  ! ... Set the I/O node
  !
  CALL io_global_start( me, 0 )
  !
  ! ... Set global coordinate for this processor
  !
  CALL mp_global_start( 0, me, gid, nproc )  
  !
  IF ( ionode ) THEN
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
  CALL mp_barrier( gid ) 
  !
  ! ... transmit npool and nimage
  !
  CALL mp_bcast( npool,  ionode_id, gid )
  CALL mp_bcast( nimage, ionode_id, gid )
  !
  IF ( nproc > maxproc ) &
     CALL errore( 'startup', ' too many processors', nproc )
  !
  ! ... all pools are initialized here
  !
  CALL init_pool()  
  !
  ! ... set the processor label for files
  !
  nd_nmbr = '   '
  !
  IF ( nproc < 10 ) THEN
     !
     WRITE( nd_nmbr(1:1) , '(I1)' ) ( me_image + 1 )
     !
  ELSE IF ( nproc < 100 ) THEN
     !
     IF ( me < 10 ) THEN
        nd_nmbr = '0'
        WRITE( nd_nmbr(2:2) , '(I1)' ) ( me_image + 1 )
     ELSE
        WRITE( nd_nmbr(1:2) , '(I2)' ) ( me_image + 1 )
     END IF
     !
  ELSE
     !
     IF ( me < 10 ) THEN
        nd_nmbr = '00'
        WRITE( nd_nmbr(3:3) , '(I1)' ) ( me_image + 1 )
     ELSE IF ( me < 100 ) THEN
        nd_nmbr = '0'
        WRITE( nd_nmbr(2:3) , '(I2)' ) ( me_image + 1 )
     ELSE
        WRITE( nd_nmbr, '(I3)' ) ( me_image + 1 )
     END IF
     !
  END IF    
  !
  ! ... stdout is printed only by the root_image (set in init_pool())
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
  IF ( mpime == root ) THEN
     !
     CALL date_and_tim( cdate, ctime )
     !
     WRITE( stdout, '(/5X,"Program ",A9," v.",A6," starts ...",&
                     &/5X,"Today is ",A9," at ",A9)' ) &
         code, version, cdate, ctime
     WRITE( stdout, '(/5X,"Parallel version (MPI)")' )
     WRITE( stdout, '(5X,"Number of processors in use:   ",I4)' ) nproc
     IF ( nimage > 1 ) &
        WRITE( stdout, '(5X,"NEB images division:  nimage = ",i4)' ) nimage
     IF ( npool > 1 ) &
        WRITE( stdout, '(5X,"K-points division:    npool  = ",i4)' ) npool
     IF ( nprocp > 1 ) &
        WRITE( stdout, '(5X,"R & G space division: nprocp = ",i4/)' ) nprocp
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
