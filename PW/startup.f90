!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine startup (nd_nmbr, version)  
  !-----------------------------------------------------------------------
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
  !  pw.x -procs 16 -npool 8 < input
  !  ORIGIN (valid also on T3E):
  !  mpirun -np 16 pw.x -npool 8 < input
  !  COMPAQ :
  !  prun -n 16 sh -c 'pw.x -npool 8 < input'
  !  In this example you will use 16 processors divided into 8 pools
  !  of 2 processors each (in this case you must have at least 8 k-points)
  !-----------------------------------------------------------------------
  !
#ifdef PARA
  use para
#endif

  ! The following two modules hold global informations about processors
  ! nuber, IDs and communicators 
  use io_global, only: io_global_start
  use mp_global, only: mp_global_start

  implicit none  
  character :: nd_nmbr * 3, version * 12  
#ifdef PARA
  include 'mpif.h'  
  character :: np * 2, cdate * 9, ctime * 9  
  external date_and_tim  
  integer :: ierr, ilen, iargc, nargs 
#ifdef T3D
  integer ipxfargc
  !
  ! set streambuffers on
  !

  call set_d_stream (1)  
#endif
  call mpi_init (ierr)  
  call error ('startup', 'mpi_init', ierr)  
  call mpi_comm_size (MPI_COMM_WORLD, nproc, ierr)  
  call error ('startup', 'mpi_comm_size', ierr)  
  call mpi_comm_rank (MPI_COMM_WORLD, me, ierr)  

  call error ('startup', 'mpi_comm_rank', ierr)  
  !
  ! This is added for compatibility with PVM notations
  ! parent process (source) will have me=1 - child process me=2,...,NPROC
  !
  me = me+1  
  if (me.eq.1) then  
     !
     ! How many pools?
     !
#ifdef T3D
     nargs = ipxfargc ()
     if (nargs.ge.2) then
        call pxfgetarg (2, np, ilen, ierr)  
#else
#ifdef HITACHI
     nargs = iargc () 
     if (nargs.ge.3) then
        call getarg (3, np)  
#else
     nargs = iargc ()
     if (nargs.ge.2) then
        call getarg (2, np)  
#endif
#endif
        read (np, '(i2)') npool  
        if (npool.le.0) npool = 1  
        if (npool.gt.nproc) npool = nproc  
     else
        npool = 1
     end if
     !
     ! set number of processes per pool (must be equal for all pool)
     !
     nprocp = nproc / npool  
     if (nproc.ne.nprocp * npool) &
          &call error ('startup','nproc.ne.nprocp*npool', 1)

  endif
  call mpi_barrier (MPI_COMM_WORLD, ierr)  

  call error ('startup', 'mpi_barrier', ierr)  
  !
  ! transmit  nprocp and npool
  !
  call mpi_bcast (nprocp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)  
  call error ('startup', 'mpi_bcast 1', ierr)  
  call mpi_bcast (npool, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)  
  call error ('startup', 'mpi_bcast 2', ierr)  
  !
  ! set the processor label for files
  !
  if (nproc.gt.maxproc) call error ('startup', ' too many processors', nproc)
  nd_nmbr = '   '  
  if (nproc.lt.10) then  
     write (nd_nmbr (1:1) , '(i1)') me  
  elseif (nproc.lt.100) then  
     if (me.lt.10) then  
        nd_nmbr = '0'  
        write (nd_nmbr (2:2) , '(i1)') me  
     else  
        write (nd_nmbr (1:2) , '(i2)') me  
     endif
  else  
     if (me.lt.10) then  
        nd_nmbr = '00'  
        write (nd_nmbr (3:3) , '(i1)') me  
     elseif (me.lt.100) then  
        nd_nmbr = '0'  
        write (nd_nmbr (2:3) , '(i2)') me  
     else  
        write (nd_nmbr, '(i3)') me  
     endif
  endif
#ifdef DEBUG
  if (me.ne.1) open (6, file = './out_'//nd_nmbr, status = 'unknown')
#else
  if (me.ne.1) open (6, file = '/dev/null', status = 'unknown')  
#endif
  if (me.eq.1) then  
     call date_and_tim (cdate, ctime)  
     write (6, 9000) version, cdate, ctime  
     write (6, '(/5x,''Parallel version (MPI)'')')  
     write (6, '(5x,''Number of processors in use:   '',i4)') nproc  
     if (npool.ne.1) &
          write (6, '(5x,''K-points division:    npool  = '',i4)') npool
     if (nprocp.ne.1)& 
          write (6, '(5x,''R & G space division: nprocp = '',i4/)') nprocp
  endif

  !
  ! Set the I/O node 
  !
  call io_global_start( (me-1), 0 )

  !
  ! Set global coordinate for this processor
  !
  call mp_global_start( 0, (me-1), MPI_COMM_WORLD, nproc )

9000 format (/5x,'Program ',a12,' starts ...'/5x, &
       &            'Today is ',a9,' at ',a9)
#endif
  return  
end subroutine startup

