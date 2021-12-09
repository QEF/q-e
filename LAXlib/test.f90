 PROGRAM lax_test
  use laxlib_descriptor
  USE laxlib_parallel_include
  use dspev_module
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  include 'laxlib_param.fh'
  include 'laxlib_hi.fh'
  include 'laxlib_low.fh'
#if defined(__MPI)
  INTEGER    STATUS(MPI_STATUS_SIZE)
#else
#define MPI_MAX_PROCESSOR_NAME 64
#endif
  INTEGER :: mype, npes, comm, ntgs, root
  LOGICAL :: iope
  INTEGER :: ierr
  INTEGER :: stdout
  !
  INTEGER :: np_ortho(2) = 1  ! size of the processor grid used in ortho
  INTEGER :: me_ortho(2) = 0  ! coordinates of the processors
  INTEGER :: me_ortho1   = 0  ! task id for the ortho group
  INTEGER :: nproc_ortho = 1  ! size of the ortho group:
                              ! of two neighbour processors in ortho_comm
  INTEGER :: ortho_comm  = 0  ! communicator for the ortho group
  INTEGER :: ortho_row_comm  = 0  ! communicator for the ortho row group
  INTEGER :: ortho_col_comm  = 0  ! communicator for the ortho col group
  INTEGER :: ortho_comm_id= 0 ! id of the ortho_comm
  INTEGER :: ortho_parent_comm  = 0  ! parent communicator from which ortho group has been created
  !
  LOGICAL :: do_distr_diag_inside_bgrp = .FALSE.
  LOGICAL :: la_proc
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  !
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
  ! the following is a workaround for Intel 12.1 bug
#if __INTEL_COMPILER  < 9999
!dir$ attributes align: 4096 :: a, s, c, d
#endif
#endif
#endif
  REAL(DP), ALLOCATABLE :: a(:,:)
  REAL(DP), ALLOCATABLE :: s(:,:)
  REAL(DP), ALLOCATABLE :: c(:,:)
  REAL(DP), ALLOCATABLE :: d(:)
  !
  REAL(DP) :: time1, time2 
  REAL*8  :: tempo(100)
  REAL*8, allocatable :: tempo_tutti(:)
  REAL*8, allocatable :: perf_matrix(:,:)
  REAL*8, allocatable :: latency_matrix(:,:)
  integer, allocatable :: perf_count(:,:)
  REAL*8  :: tempo_mio(100)
  REAL*8  :: tempo_min(100)
  REAL*8  :: tempo_max(100)
  REAL*8  :: tempo_avg(100)

  TYPE(la_descriptor) :: desc
  INTEGER :: idesc(LAX_DESC_SIZE)
  INTEGER :: i, ir, ic, nx, n, nr, nc  ! size of the matrix
  INTEGER :: n_in, nlen, dest, sour, tag, ii
  INTEGER :: n_diag ! number of MPI processes participated in diagonalization
  INTEGER :: nnodes ! number of nodes by hostname
  !
  integer :: nargs
  CHARACTER(LEN=80) :: arg
  CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME), allocatable :: proc_name(:)
  CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME), allocatable :: node_name(:)
  INTEGER, allocatable :: proc2node(:)
#if defined(_OPENMP)
  INTEGER, EXTERNAL :: omp_get_max_threads
#endif
  !
#if defined(_OPENMP)
  INTEGER :: PROVIDED
#endif
  !
  !   ........
  !
  !   default parameter 
  !
  n_in = 512
  !
  nargs = command_argument_count()
  do i = 1, nargs - 1
     CALL get_command_argument(i, arg)
     IF( TRIM( arg ) == '-n' ) THEN
        CALL get_command_argument(i+1, arg)
        READ( arg, * ) n_in
     END IF
  end do

#if defined(__MPI)

#if defined(_OPENMP)
  CALL MPI_Init_thread(MPI_THREAD_FUNNELED, PROVIDED, ierr)
#else
  CALL MPI_Init(ierr)
#endif

  CALL mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  comm = MPI_COMM_WORLD
  ntgs = 1
  root = 0
  IF(mype==root) THEN
     iope = .true.
  ELSE
     iope = .false.
  ENDIF

#else

  mype = 0
  npes = 1
  comm = 0
  ntgs = 1
  root = 0
  iope = .true.

#endif


  !OPEN ( unit = 6, file = TRIM('test.out'), status='unknown' )

  !
  !write(6,*) 'mype = ', mype, ' npes = ', npes
  !
  !
  !  Broadcast input parameter first
  !
#if defined(__MPI)
  CALL MPI_BCAST(n_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
#endif

  n = n_in

  if( mype == 0 ) then

    write(6,*) '+-----------------------------------+'
    write(6,*) '|         QE Linear Algebra         |'
    write(6,*) '|          testing & timing         |'
    write(6,*) '|         by Carlo Cavazzoni        |'
    write(6,*) '+-----------------------------------+'
    write(6,*)
    write(6,*) 'matrix size = ', n, ' x ', n
    write(6,*) 'num. procs  = ', npes
#if defined(_OPENMP)
    write(6,*) 'thr x proc  = ', omp_get_max_threads()
#endif
    write(6,*)

  endif

  allocate( proc_name( npes ) )
  allocate( node_name( npes ) )
  allocate( proc2node( npes ) )
  nnodes = 0

  do i = 1, npes
#if defined(__MPI)
     if( mype == i-1 ) then
         call MPI_Get_processor_name( proc_name(i), nlen, ierr )
     end if
     CALL MPI_BCAST( nlen, 1, MPI_INT, i-1, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST( proc_name(i), MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, i-1, MPI_COMM_WORLD, ierr )
#else
     proc_name(i) = 'localhost'
#endif
     if( mype == 0 ) then
!        write(6,310)  i, proc_name(i)
     end if
310 FORMAT('pe = ',I5,' name = ', A20) 
    do ii = 1, nnodes
       if( proc_name(i) == node_name(ii) ) then
          exit
       end if 
    end do
    if( ii > nnodes ) then
       nnodes = nnodes + 1
       node_name( nnodes ) = proc_name( i )
    end if
    proc2node( i ) = ii
  end do
  !
  if( mype == 0 ) then
    write(6,*) '+-----------------------------------+'
    write(6,*) '|  node list                        |'
    write(6,*) '+-----------------------------------+'
    write(6,*) 
    do ii = 1, nnodes
       write(6,310)  ii, node_name(ii)
    end do
  end if
  
  allocate( perf_matrix( nnodes, nnodes ) )
  allocate( latency_matrix( nnodes, nnodes ) )
  allocate( perf_count( nnodes, nnodes ) )
  perf_matrix = 0.0d0
  latency_matrix = 0.0d0
  perf_count = 0

  ! Check core speed
  !
#if defined(__MPI)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif
  nx = 1024
  ALLOCATE( s( nx, nx ) )
  ALLOCATE( a( nx, nx ) )
  ALLOCATE( c( nx, nx ) )
  ALLOCATE( tempo_tutti( npes ) )
  tempo_tutti = 0.0d0
  a = 1.0d0
  s = 1.0d0
  c = 1.0d0
  CALL dgemm('n', 'n', nx, nx, nx, 1.0d0, A, nx, s, nx, 1.0d0, C, nx)
  tempo(1) = mpi_wall_time()
  CALL dgemm('n', 'n', nx, nx, nx, 1.0d0, A, nx, s, nx, 1.0d0, C, nx)
  CALL dgemm('n', 'n', nx, nx, nx, 1.0d0, A, nx, s, nx, 1.0d0, C, nx)
  CALL dgemm('n', 'n', nx, nx, nx, 1.0d0, A, nx, s, nx, 1.0d0, C, nx)
  CALL dgemm('n', 'n', nx, nx, nx, 1.0d0, A, nx, s, nx, 1.0d0, C, nx)
  CALL dgemm('n', 'n', nx, nx, nx, 1.0d0, A, nx, s, nx, 1.0d0, C, nx)
  tempo(2) = mpi_wall_time()
  DEALLOCATE( s )
  DEALLOCATE( a )
  DEALLOCATE( c )
  tempo_tutti(mype+1) = tempo(2)-tempo(1)
#if defined(__MPI)
  CALL MPI_ALLREDUCE( MPI_IN_PLACE, tempo_tutti, npes, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
#endif
  if( mype == 0 ) then
     write(6,*)
     write(6,*)
     write(6,*)
     write(6,*) '+-----------------------------------+'
     write(6,*) '|    measured task performances     |'
     write(6,*) '+-----------------------------------+'
     do i = 1, npes
        write(6,300)  i, 5.0d0*DBLE(nx*nx*nx)*2.0d0/tempo_tutti(i)/1.0D+9, proc_name(i)
     end do
  end if
300 FORMAT('pe = ',I5,',', F8.3, ' GFlops', ',  node: ', A20) 
  !
  ! Check network speed
  !
  nx = 2048
  ALLOCATE( s( nx, nx ) )
  tempo_tutti = 0.0d0

  do ii = 0, npes-1
  do i = 0, ii-1
     sour = ii
     dest = i
     tag = i + ii * npes
#if defined(__MPI)
     CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
     if( ( mype == sour ) .or. ( mype == dest ) ) THEN
        tempo(1) = mpi_wall_time()
        if( mype == dest ) then
           CALL MPI_SEND(s, nx*nx, MPI_DOUBLE_PRECISION, sour, TAG, MPI_COMM_WORLD, ierr)
           CALL MPI_RECV(s, nx*nx, MPI_DOUBLE_PRECISION, sour, TAG+NPES*NPES, MPI_COMM_WORLD, status, ierr)
        else if( mype == sour ) then
           CALL MPI_RECV(s, nx*nx, MPI_DOUBLE_PRECISION, dest, TAG, MPI_COMM_WORLD, status, ierr)
           CALL MPI_SEND(s, nx*nx, MPI_DOUBLE_PRECISION, dest, TAG+NPES*NPES, MPI_COMM_WORLD, ierr)
        endif
        tempo(2) = mpi_wall_time()
        perf_matrix( proc2node( ii+1 ), proc2node( i+1 ) ) = perf_matrix( proc2node( ii+1 ), proc2node( i+1 ) ) + &
           2.0d0*DBLE(nx*nx)*8.0d0/(tempo(2)-tempo(1))/1.0D+9
        perf_count( proc2node( ii+1 ), proc2node( i+1 ) ) = perf_count( proc2node( ii+1 ), proc2node( i+1 ) ) + 1
     END IF
     CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
     if( ( mype == sour ) .or. ( mype == dest ) ) THEN
        tempo(1) = mpi_wall_time()
        if( mype == dest ) then
           CALL MPI_SEND(ii, 1, MPI_BYTE, sour, TAG, MPI_COMM_WORLD, ierr)
           CALL MPI_RECV(ii, 1, MPI_BYTE, sour, TAG+NPES, MPI_COMM_WORLD, status, ierr)
        else if( mype == sour ) then
           CALL MPI_RECV(ii, 1, MPI_BYTE, dest, TAG, MPI_COMM_WORLD, status, ierr)
           CALL MPI_SEND(ii, 1, MPI_BYTE, dest, TAG+NPES, MPI_COMM_WORLD, ierr)
        endif
        tempo(2) = mpi_wall_time()
        latency_matrix( proc2node( ii+1 ), proc2node( i+1 ) ) = latency_matrix( proc2node( ii+1 ), proc2node( i+1 ) ) + &
           (tempo(2)-tempo(1))
     END IF
#endif
  end do
  end do
#if defined(__MPI)
  CALL MPI_ALLREDUCE( MPI_IN_PLACE, perf_matrix, SIZE(perf_matrix), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( MPI_IN_PLACE, latency_matrix, SIZE(latency_matrix), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( MPI_IN_PLACE, perf_count, SIZE(perf_count), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
#endif
  if( mype == 0 ) then
     write(6,*)
     write(6,*)
     write(6,*)
     write(6,*) '+-----------------------------------+'
     write(6,*) '|    ping-pong network bandwidth    |'
     write(6,*) '+-----------------------------------+'
     write(6,*) 
     do ii = 1, nnodes
        do i = 1, nnodes
           if( perf_count(i,ii) > 0 ) then
             perf_matrix(i,ii) = perf_matrix(i,ii) / perf_count(i,ii)
             write( 6, 314 ) node_name(i), node_name(ii), perf_count(i,ii), perf_matrix(i,ii)
           end if
        end do
     end do
314 FORMAT( A20, A20, I5, ':', F8.3, 'GBytes') 
     write(6,*)
     write(6,*) '+-----------------------------------+'
     write(6,*) '|    ping-pong network latency      |'
     write(6,*) '+-----------------------------------+'
     write(6,*) 
     do ii = 1, nnodes
        do i = 1, nnodes
           if( perf_count(i,ii) > 0 ) then
             latency_matrix(i,ii) = latency_matrix(i,ii) / perf_count(i,ii)
             write( 6, 315 ) node_name(i), node_name(ii), perf_count(i,ii), latency_matrix(i,ii)*1000000.0d0
           end if
        end do
     end do
315 FORMAT( A20, A20, I5, ':', F10.3, 'usec') 
  end if
  DEALLOCATE( s )
  DEALLOCATE( tempo_tutti )
  !
  !
  n_diag = n
  CALL laxlib_start(n_diag, mpi_comm_world, do_distr_diag_inside_bgrp)
  CALL laxlib_getval( np_ortho = np_ortho, ortho_comm = ortho_comm, &
    do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp )
  !
  ALLOCATE(rank_ip( np_ortho(1), np_ortho(2) ))
  ALLOCATE( irc_ip( np_ortho(1) ), nrc_ip (np_ortho(1) ) )
  CALL desc_init( n, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip  )
  !
  CALL laxlib_intarray_to_desc( desc, idesc )
  !
  nx = 1
  IF( desc%active_node > 0 ) nx = desc%nrcx
  !
  ALLOCATE( d( n ) )
  ALLOCATE( s( nx, nx ) )
  ALLOCATE( a( nx, nx ) )
  ALLOCATE( c( nx, nx ) )

  nr = desc%nr
  nc = desc%nc
  ir = desc%ir
  ic = desc%ic
  !
  ! do not take the time of the first execution, it may be biased from MPI
  ! initialization stuff
  !
  CALL set_a()
  !
  CALL diagonalize_parallel_x( n, a, d, s, idesc )
  !
  tempo = 0.0d0
  tempo_mio = 0.0d0
  tempo_min = 0.0d0
  tempo_max = 0.0d0
  tempo_avg = 0.0d0


  CALL set_a()
  !
#if defined(__MPI)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif
  tempo(1) = mpi_wall_time()
  !
  CALL diagonalize_parallel_x( n, a, d, s, idesc )
  !
#if defined(__MPI)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif
  tempo(2) = mpi_wall_time()
  !
  CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, a, nx, s, nx, 0.0d0, c, nr, idesc)
  !
#if defined(__MPI)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif
  tempo(3) = mpi_wall_time()
  !
  do i = 2, 10
     tempo_mio(i) = tempo(i)-tempo(i-1)
  end do
#if defined(__MPI)
  CALL MPI_ALLREDUCE( tempo_mio, tempo_min, 100, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( tempo_mio, tempo_max, 100, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( tempo_mio, tempo_avg, 100, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
#else
  tempo_min = tempo_mio
  tempo_max = tempo_mio
  tempo_avg = tempo_mio
#endif

  tempo_avg = tempo_avg / npes
  !
  IF( mype == 0 ) THEN
     write(6,*) 
     write(6,*) ' Matrix eigenvalues '
     write(6,*) 
     IF ( n <= 16 ) THEN
       DO i = 1, n
         write(6,*) ' D(',i,')=',d(i)
       END DO
     ELSE
       DO i = 1, 8
         write(6,*) ' D(',i,')=',d(i)
       END DO
       write(6,*) ' ... '
       DO i = n-8, n 
         write(6,*) ' D(',i,')=',d(i)
       END DO
     END IF
     write(6,*) 
  ENDIF

  if( mype == 0 ) then

    write(6,*) '**** LA Timing ****'
    write(6,*) 
    write(6,200) 2.0d0*n*n*n / 1.D9 / tempo_avg(3)
    write(6,*) 

    write(6,100)
    write(6,1)
    write(6,100)
    write(6,2) tempo_min(2), tempo_max(2), tempo_avg(2)
    write(6,100)
    write(6,3) tempo_min(3), tempo_max(3), tempo_avg(3)
    write(6,100)

200 FORMAT(' GFlops = ', F14.2 )
100 FORMAT(' +--------------------+----------------+-----------------+----------------+' )
1   FORMAT(' |LAX subroutine      |  sec. min      | sec. max        | sec.  avg      |' )
2   FORMAT(' |diagonalize_parallel| ',    D14.3, ' | ',   D14.3,  '  | ', D14.3,    ' |' )
3   FORMAT(' |sqr_mm_cannon       | ',    D14.3, ' | ',   D14.3,  '  | ', D14.3,    ' |' )


  end if

  deallocate( proc_name )
  deallocate( node_name )
  deallocate( proc2node )
  deallocate( perf_matrix )
  deallocate( latency_matrix )
  deallocate( perf_count )
  deallocate( a, s, d, c )
  deallocate( rank_ip, irc_ip, nrc_ip )

#if defined(__MPI)
  CALL mpi_finalize(ierr)
#endif


contains

  SUBROUTINE set_a()
     INTEGER :: i, j, ii, jj
     IF( desc%active_node < 0 ) RETURN
     DO j = 1, nc
         DO i = 1, nr
            ii = i + ir - 1
            jj = j + ic - 1
            IF( ii == jj ) THEN
               a(i,j) = ( DBLE( n-ii+1 ) ) / DBLE( n ) + 1.0d0 / ( DBLE( ii+jj ) - 1.0d0 )
            ELSE
               a(i,j) = 1.0d0 / ( DBLE( ii+jj ) - 1.0d0 )
            END IF
         END DO
      END DO
      RETURN
   END SUBROUTINE set_a

  function mpi_wall_time ()
    real*8 :: mpi_wall_time
#if defined(__MPI)
    mpi_wall_time = MPI_WTIME()
#else
    ! standard way to get the wall time, sometimes with very low precision
    integer :: cr, nc
    real*8, save :: t0 = -1.0
    !
    call system_clock(count_rate=cr)
    call system_clock(count=nc)
    if ( t0 < 0.0 ) t0 = dble(nc)/cr
    mpi_wall_time = dble(nc)/cr - t0
#endif
  end function mpi_wall_time
end program lax_test
