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
  !  This subroutine initializes MPI
  !  Processes are organized in NPOOL pools each dealing with a subset of
  !  kpoints. Within each pool R & G space distribution is performed.
  !  NPROC is read from command line or can be set with the appropriate
  !  environment variable ( for example use 'setenv MP_PROCS 8' on IBM SP
  !  machine to run on NPROC=8 processors ) and NPOOL is read from command
  !  line.
  !  NPOOL must be a whole divisor of NPROC
  !  An example without any environment variable set is the following:
  !  T3E :
  !  mpprun -n 16 pw.x -npool 8 < input
  !  IBM SP :
  !  poe pw.x -procs 16 -npool 8 < input
  !  ORIGIN /PC clusters using "mpirun" :
  !  mpirun -np 16 pw.x -npool 8 < input
  !  COMPAQ :
  !  prun -n 16 sh -c 'pw.x -npool 8 < input'
  !  PC clusters using "mpiexec" :
  !  mpiexec -n 16 pw.x -npool 8 < input  
  !  In this example you will use 16 processors divided into 8 pools
  !  of 2 processors each (in this case you must have at least 8 k-points)
  !
  !  Use "-input filename" to read input from file "filename":
  !  may be useful if you have trouble reading from standard input
  !-----------------------------------------------------------------------
  !
  ! ... The following two modules hold global information about processors
  ! ... number, IDs and communicators
  !
  USE io_global,  ONLY :  stdout, io_global_start, ionode_id
  USE mp_global,  ONLY :  mp_global_start
  USE mp,         ONLY :  mp_start, mp_env, mp_barrier, mp_bcast
#ifdef __PARA
  USE para_const, ONLY :  maxproc
  USE para,       ONLY :  me, mypool, nproc, npool, nprocp 
#endif  
  !
  IMPLICIT NONE
  !
  CHARACTER :: nd_nmbr*3, code*9, version*6
  INTEGER   :: gid
#if ! defined __PARA
  INTEGER   :: me, nproc
#endif
  CHARACTER :: np*80, cdate * 9, ctime * 9
  EXTERNAL     date_and_tim
#ifdef __PARA
  INTEGER   :: ierr, ilen, iargc, nargs, iiarg
  !
  !
#  ifdef __T3E
  !
  ! ... set streambuffers on
  !
  CALL set_d_stream( 1 )
#  endif
  !
#endif
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
  ! ... This is added for compatibility with PVM notations
  ! ... parent process (source) will have me=1 - child process me=2,...,NPROC
  !
  me = me + 1
  !
#ifdef __PARA
  !
  IF ( me == 1 ) THEN
     !
     ! ... How many pools ?
     !
     npool = 1
     nargs = IARGC() 
     !
     DO iiarg = 1, ( nargs - 1 )
        CALL getarg( iiarg, np )  
        IF ( TRIM( np ) == '-npool' .OR. TRIM( np ) == '-npools' ) THEN
          CALL GETARG( ( iiarg + 1 ), np )  
          READ(np,*) npool  
        END IF
        IF ( TRIM( np ) == '-input' .OR. TRIM( np ) == '-inp' .OR. &
             TRIM( np ) == '-in' ) THEN
           CALL GETARG( ( iiarg + 1 ) , np )  
           OPEN ( UNIT = 5, FILE = np, FORM = 'formatted', &
                STATUS = 'old', IOSTAT = ierr)
           CALL errore( 'startup', 'input file '//TRIM(np)//' not found',&
                ierr )
        end if

     END DO
     npool = MAX( npool, 1 )
     npool = MIN( npool, nproc )
     !
     ! ... set number of processes per pool (must be equal for all pools)
     !
     nprocp = nproc / npool
     IF ( nproc /= ( nprocp * npool ) ) &
        CALL errore( 'startup', 'nproc /= nprocp*npool', 1 )
     !
  END IF
  !
  CALL mp_barrier( gid ) 
  !
  ! ... transmit nprocp and npool
  !
  CALL mp_bcast( nprocp, ionode_id, gid )
  CALL mp_bcast( npool, ionode_id, gid )
  !
  ! ... set the processor label for files
  !
  IF ( nproc > maxproc ) &
     call errore( 'startup', ' too many processors', nproc )
  nd_nmbr = '   '
  IF ( nproc < 10 ) THEN
     WRITE( nd_nmbr(1:1) , '(I1)' ) me
  ELSE IF ( nproc < 100 ) THEN
     IF ( me < 10 ) THEN
        nd_nmbr = '0'
        WRITE( nd_nmbr(2:2) , '(I1)' ) me
     ELSE
        WRITE( nd_nmbr(1:2) , '(I2)' ) me
     END IF
  ELSE
     IF ( me < 10 ) THEN
        nd_nmbr = '00'
        WRITE( nd_nmbr(3:3) , '(I1)' ) me
     ELSE IF ( me < 100 ) THEN
        nd_nmbr = '0'
        WRITE( nd_nmbr(2:3) , '(I2)' ) me
     ELSE
        WRITE( nd_nmbr, '(I3)' ) me
     END IF
  END IF
  !
  ! ... pools are initialized
  !
  CALL init_pool()  
  !
#  ifdef DEBUG
  IF ( me /= 1 .OR. mypool /= 1 ) &
     OPEN( UNIT = stdout, FILE = './out_'//nd_nmbr, STATUS = 'UNKNOWN' )
#  else
  IF ( me /= 1 .OR. mypool /= 1 ) &
     OPEN( UNIT = stdout, FILE = '/dev/null', STATUS = 'UNKNOWN' )
#  endif
  !
  IF ( me == 1 ) THEN
     !
     CALL date_and_tim( cdate, ctime )
     !
     WRITE( stdout, 9000 ) code, version, cdate, ctime
     WRITE( stdout, '(/5X,"Parallel version (MPI)")' )
     WRITE( stdout, '(5X,"Number of processors in use:   ",I4)' ) nproc
     IF ( npool /= 1 ) &
        WRITE( stdout, '(5X,"K-points division:    npool  = ",i4)' ) npool
     IF ( nprocp /= 1 ) &
        WRITE( stdout, '(5X,"R & G space division: nprocp = ",i4/)' ) nprocp
     !
  END IF
  !
#else
  !
  nd_nmbr = '   '
  CALL date_and_tim( cdate, ctime )
  WRITE( stdout, 9000 ) code, version, cdate, ctime
  !
#endif
  !
  RETURN
  !
9000 FORMAT( /5X,'Program ',A9,' v.',A6,' starts ...',/5X, &
           &     'Today is ',A9,' at ',A9)
  !     
END SUBROUTINE startup

