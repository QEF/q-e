program lax_test
  use descriptors
  use dspev_module
  IMPLICIT NONE
#ifdef __MPI
  include 'mpif.h'
#endif
#include "la_param.f90"
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
#if defined __SCALAPACK
  INTEGER :: me_blacs   =  0  ! BLACS processor index starting from 0
  INTEGER :: np_blacs   =  1  ! BLACS number of processor
#endif
  !
  INTEGER :: world_cntx = -1  ! BLACS context of all processor
  INTEGER :: ortho_cntx = -1  ! BLACS context for ortho_comm
  !
  REAL(DP), ALLOCATABLE :: a(:,:)
  REAL(DP), ALLOCATABLE :: s(:,:)
  REAL(DP), ALLOCATABLE :: c(:,:)
  REAL(DP), ALLOCATABLE :: d(:)
  !
  REAL(DP) :: time1, time2 
  REAL*8  :: tempo(100)
  REAL*8  :: tempo_mio(100)
  REAL*8  :: tempo_min(100)
  REAL*8  :: tempo_max(100)
  REAL*8  :: tempo_avg(100)

  TYPE(la_descriptor) :: desc
  INTEGER :: i, ir, ic, nx, n, nr, nc  ! size of the matrix
  INTEGER :: n_in
  !
  integer :: nargs
  CHARACTER(LEN=80) :: arg
  !
#if defined(__OPENMP)
  INTEGER :: PROVIDED
#endif
  !
  !   ........
  !
  !   default parameter 
  !
  n_in = 1024
  !
  nargs = command_argument_count()
  do i = 1, nargs - 1
     CALL get_command_argument(i, arg)
     IF( TRIM( arg ) == '-n' ) THEN
        CALL get_command_argument(i+1, arg)
        READ( arg, * ) n_in
     END IF
  end do

#ifdef __MPI

#if defined(__OPENMP)
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
  !
  !write(*,*) 'mype = ', mype, ' npes = ', npes
  !
  !
  !  Broadcast input parameter first
  !
  CALL MPI_BCAST(n_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

  n = n_in

  if( mype == 0 ) then

    write(*,*) '+-----------------------------------+'
    write(*,*) '|         QE Linear Algebra         |'
    write(*,*) '|          testing & timing         |'
    write(*,*) '|         by Carlo Cavazzoni        |'
    write(*,*) '+-----------------------------------+'
    write(*,*)
    write(*,*) 'matrix size = ', n, ' x ', n
    write(*,*) 'num. procs  = ', npes
    write(*,*)

  endif

  call mp_start_diag()
  !
  CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
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
  CALL diagonalize_parallel( n, a, d, s, desc )
  !
  tempo = 0.0d0
  tempo_mio = 0.0d0
  tempo_min = 0.0d0
  tempo_max = 0.0d0
  tempo_avg = 0.0d0

  CALL set_a()
  !
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  tempo(1) = MPI_WTIME()
  !
  CALL diagonalize_parallel( n, a, d, s, desc )
  !
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  tempo(2) = MPI_WTIME()
  !
  CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, a, nx, s, nx, 0.0d0, c, nr, desc)
  !
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  tempo(3) = MPI_WTIME()
  !
  do i = 2, 10
     tempo_mio(i) = tempo(i)-tempo(i-1)
  end do
  !
#ifdef __MPI
  CALL MPI_ALLREDUCE( tempo_mio, tempo_min, 100, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( tempo_mio, tempo_max, 100, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  CALL MPI_ALLREDUCE( tempo_mio, tempo_avg, 100, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
#else
  tempo_min = tempo
  tempo_max = tempo
#endif

  tempo_avg = tempo_avg / npes
  !
  IF( mype == 0 ) THEN
     write(*,*) 
     write(*,*) ' Matrix eigenvalues '
     write(*,*) 
     IF ( n <= 16 ) THEN
       DO i = 1, n
         write(*,*) ' D(',i,')=',d(i)
       END DO
     ELSE
       DO i = 1, 8
         write(*,*) ' D(',i,')=',d(i)
       END DO
       write(*,*) ' ... '
       DO i = n-8, n 
         write(*,*) ' D(',i,')=',d(i)
       END DO
     END IF
     write(*,*) 
  ENDIF

  if( mype == 0 ) then

    write(*,*) '**** LA Timing ****'
    write(*,*) 
    write(*,200) 2.0d0*n*n*n / 1.D9 / tempo_avg(3)
    write(*,*) 

    write(*,100)
    write(*,1)
    write(*,100)
    write(*,2) tempo_min(2), tempo_max(2), tempo_avg(2)
    write(*,100)
    write(*,3) tempo_min(3), tempo_max(3), tempo_avg(3)
    write(*,100)

200 FORMAT(' GFlops = ', F14.2 )
100 FORMAT(' +--------------------+----------------+-----------------+----------------+' )
1   FORMAT(' |LAX subroutine      |  sec. min      | sec. max        | sec.  avg      |' )
2   FORMAT(' |diagonalize_parallel| ',    D14.3, ' | ',   D14.3,  '  | ', D14.3,    ' |' )
3   FORMAT(' |sqr_mm_cannon       | ',    D14.3, ' | ',   D14.3,  '  | ', D14.3,    ' |' )


  end if



#ifdef __MPI
  CALL mpi_finalize(ierr)
#endif


contains

  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_diag( )
    !---------------------------------------------------------------------------
    !
    ! ... Ortho/diag/linear algebra group initialization
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr = 0
    INTEGER :: color, key, nproc_try

#if defined __SCALAPACK
    INTEGER, ALLOCATABLE :: blacsmap(:,:)
    INTEGER :: ortho_cntx_pe
    INTEGER :: nprow, npcol, myrow, mycol, i, j, k
    INTEGER, EXTERNAL :: BLACS_PNUM
    !
    INTEGER :: nparent=1
    INTEGER :: total_nproc=1
    INTEGER :: total_mype=0
    INTEGER :: nproc_parent=1
    INTEGER :: my_parent_id=0
#endif

    !
#if defined __SCALAPACK
    !
    CALL mpi_comm_rank( MPI_COMM_WORLD, me_blacs, ierr)
    CALL mpi_comm_size( MPI_COMM_WORLD, np_blacs, ierr)
    !
    ! define a 1D grid containing all MPI tasks of the global communicator
    ! NOTE: world_cntx has the MPI communicator on entry and the BLACS context
    ! on exit
    !       BLACS_GRID_INIT() will create a copy of the communicator, which can
    !       be
    !       later retrieved using CALL BLACS_GET(world_cntx, 10, comm_copy)
    !
    world_cntx = MPI_COMM_WORLD
    CALL BLACS_GRIDINIT( world_cntx, 'Row', 1, np_blacs )
    !
#endif
    !
    ! the ortho group for parallel linear algebra is a sub-group of the pool,
    ! then there are as many ortho groups as pools.
    !
#if defined __MPI

    nproc_try = MAX( npes, 1 )

    !  find the square closer (but lower) to nproc_try
    !
    CALL grid2d_dims( 'S', nproc_try, np_ortho(1), np_ortho(2) )
    !
    !  now, and only now, it is possible to define the number of tasks
    !  in the ortho group for parallel linear algebra
    !
    nproc_ortho = np_ortho(1) * np_ortho(2)
    !
    !  here we choose the first "nproc_ortho" processors
    !
    color = 0
    IF( mype < nproc_ortho ) color = 1
    !
    key   = mype
    !
    !  initialize the communicator for the new group by splitting the input
    !  communicator
    !
    CALL mpi_comm_split( MPI_COMM_WORLD , color, key, ortho_comm, ierr )
    !
    ! Computes coordinates of the processors, in row maior order
    !
    CALL mpi_comm_rank( ortho_comm, me_ortho1, ierr)
    !
    IF( mype == 0 .AND. me_ortho1 /= 0 ) &
         CALL lax_error__( " init_ortho_group ", " wrong root task in ortho group ", ierr )
    !
    if( color == 1 ) then
       ! this task belong to the ortho_group compute its coordinates
       ortho_comm_id = 1
       CALL GRID2D_COORDS( 'R', me_ortho1, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
       CALL GRID2D_RANK( 'R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr )
       IF( ierr /= me_ortho1 ) &
            CALL lax_error__( " init_ortho_group ", " wrong task coordinates in ortho group ", ierr )
       IF( me_ortho1 /= mype ) &
            CALL lax_error__( " init_ortho_group ", " wrong rank assignment in ortho group ", ierr )
       CALL mpi_comm_split( ortho_comm , me_ortho(2), me_ortho(1), ortho_col_comm, ierr )
       CALL mpi_comm_split( ortho_comm , me_ortho(1), me_ortho(2), ortho_row_comm, ierr )
    else
       ! this task does NOT belong to the ortho_group set dummy values
       ortho_comm_id = 0
       me_ortho(1) = me_ortho1
       me_ortho(2) = me_ortho1
    endif

#if defined __SCALAPACK
    !
    !  This part is used to eliminate the image dependency from ortho groups
    !  SCALAPACK is now independent of whatever level of parallelization
    !  is present on top of pool parallelization
    !
    total_nproc = npes
    total_mype = mype
    !
    ALLOCATE( blacsmap( np_ortho(1), np_ortho(2) ) )

    CALL BLACS_GET( world_cntx, 10, ortho_cntx_pe ) ! retrieve communicator of world context
    blacsmap = 0
    nprow = np_ortho(1)
    npcol = np_ortho(2)

    IF( ortho_comm_id > 0  ) THEN

       blacsmap( me_ortho(1) + 1, me_ortho(2) + 1 ) = BLACS_PNUM( world_cntx, 0, me_blacs )

    END IF

   ! All MPI tasks defined in the global communicator take part in the definition of the BLACS grid

   CALL MPI_ALLREDUCE( MPI_IN_PLACE, blacsmap, SIZE( blacsmap ), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

   CALL BLACS_GRIDMAP( ortho_cntx_pe, blacsmap, nprow, nprow, npcol)
   CALL BLACS_GRIDINFO( ortho_cntx_pe, nprow, npcol, myrow, mycol )

   IF( ortho_comm_id > 0) THEN

        IF(  np_ortho(1) /= nprow ) &
           CALL lax_error__( ' init_ortho_group ', ' problem with SCALAPACK, wrong no. of task rows ', 1 )
        IF(  np_ortho(2) /= npcol ) &
           CALL lax_error__( ' init_ortho_group ', ' problem with SCALAPACK, wrong no. of task columns ', 1 )
        IF(  me_ortho(1) /= myrow ) &
           CALL lax_error__( ' init_ortho_group ', ' problem with SCALAPACK, wrong task row ID ', 1 )
        IF(  me_ortho(2) /= mycol ) &
           CALL lax_error__( ' init_ortho_group ', ' problem with SCALAPACK, wrong task columns ID ', 1 )

        ortho_cntx = ortho_cntx_pe

    END IF

    DEALLOCATE( blacsmap )


#endif

#else

    ortho_comm_id = 1

#endif

    RETURN
  END SUBROUTINE 



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



end program lax_test
