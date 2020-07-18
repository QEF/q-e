!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

SUBROUTINE laxlib_end()
  use laxlib_processors_grid
  CALL laxlib_end_drv ( )
END SUBROUTINE laxlib_end


SUBROUTINE laxlib_getval_ ( nproc_ortho, leg_ortho, np_ortho, me_ortho, ortho_comm, ortho_row_comm, ortho_col_comm, &
  ortho_comm_id, ortho_parent_comm, me_blacs, np_blacs, ortho_cntx, world_cntx, do_distr_diag_inside_bgrp  )
  use laxlib_processors_grid, ONLY : &
    nproc_ortho_ => nproc_ortho, &
    leg_ortho_   => leg_ortho, &
    np_ortho_    => np_ortho, &
    me_ortho_    => me_ortho, &
    ortho_comm_  => ortho_comm, & 
    ortho_row_comm_ => ortho_row_comm, &
    ortho_col_comm_ => ortho_col_comm, & 
    ortho_comm_id_  => ortho_comm_id, &
    ortho_parent_comm_ => ortho_parent_comm, &
    me_blacs_    => me_blacs,  &
    np_blacs_    => np_blacs, &
    ortho_cntx_  => ortho_cntx, &
    world_cntx_  => world_cntx, &
    do_distr_diag_inside_bgrp_ => do_distr_diag_inside_bgrp
  IMPLICIT NONE
  INTEGER, OPTIONAL, INTENT(OUT) :: nproc_ortho
  INTEGER, OPTIONAL, INTENT(OUT) :: leg_ortho
  INTEGER, OPTIONAL, INTENT(OUT) :: np_ortho(2)
  INTEGER, OPTIONAL, INTENT(OUT) :: me_ortho(2)
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_row_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_col_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_comm_id
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_parent_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: me_blacs
  INTEGER, OPTIONAL, INTENT(OUT) :: np_blacs
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_cntx
  INTEGER, OPTIONAL, INTENT(OUT) :: world_cntx
  LOGICAL, OPTIONAL, INTENT(OUT) :: do_distr_diag_inside_bgrp
  IF( PRESENT(nproc_ortho) ) nproc_ortho = nproc_ortho_
  IF( PRESENT(leg_ortho) ) leg_ortho = leg_ortho_
  IF( PRESENT(np_ortho) ) np_ortho = np_ortho_
  IF( PRESENT(me_ortho) ) me_ortho = me_ortho_
  IF( PRESENT(ortho_comm) ) ortho_comm = ortho_comm_
  IF( PRESENT(ortho_row_comm) ) ortho_row_comm = ortho_row_comm_
  IF( PRESENT(ortho_col_comm) ) ortho_col_comm = ortho_col_comm_
  IF( PRESENT(ortho_comm_id) ) ortho_comm_id = ortho_comm_id_
  IF( PRESENT(ortho_parent_comm) ) ortho_parent_comm = ortho_parent_comm_
  IF( PRESENT(me_blacs) ) me_blacs = me_blacs_
  IF( PRESENT(np_blacs) ) np_blacs = np_blacs_
  IF( PRESENT(ortho_cntx) ) ortho_cntx = ortho_cntx_
  IF( PRESENT(world_cntx) ) world_cntx = world_cntx_
  IF( PRESENT(do_distr_diag_inside_bgrp) ) do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp_
END SUBROUTINE
!
SUBROUTINE laxlib_get_status_x ( lax_status )
  use laxlib_processors_grid, ONLY : &
    nproc_ortho_ => nproc_ortho, &
    leg_ortho_   => leg_ortho, &
    np_ortho_    => np_ortho, &
    me_ortho_    => me_ortho, &
    ortho_comm_  => ortho_comm, & 
    ortho_row_comm_ => ortho_row_comm, &
    ortho_col_comm_ => ortho_col_comm, & 
    ortho_comm_id_  => ortho_comm_id, &
    ortho_parent_comm_ => ortho_parent_comm, &
    me_blacs_    => me_blacs,  &
    np_blacs_    => np_blacs, &
    ortho_cntx_  => ortho_cntx, &
    world_cntx_  => world_cntx, &
    do_distr_diag_inside_bgrp_ => do_distr_diag_inside_bgrp
  IMPLICIT NONE
  include 'laxlib_param.fh'
  INTEGER, INTENT(OUT) :: LAX_STATUS(:)
  lax_status(LAX_STATUS_NPROC)= nproc_ortho_
  lax_status(LAX_STATUS_LEG)= leg_ortho_
  lax_status(LAX_STATUS_NP1)= np_ortho_( 1 )
  lax_status(LAX_STATUS_NP2)= np_ortho_( 2 )
  lax_status(LAX_STATUS_ME1)= me_ortho_( 1 )
  lax_status(LAX_STATUS_ME2)= me_ortho_( 2 )
  lax_status(LAX_STATUS_COMM)= ortho_comm_
  lax_status(LAX_STATUS_ROWCOMM)= ortho_row_comm_
  lax_status(LAX_STATUS_COLCOMM)= ortho_col_comm_
  lax_status(LAX_STATUS_COMMID)= ortho_comm_id_
  lax_status(LAX_STATUS_PARENTCOMM)= ortho_parent_comm_
  lax_status(LAX_STATUS_MEBLACS)= me_blacs_
  lax_status(LAX_STATUS_NPBLACS)= np_blacs_
  lax_status(LAX_STATUS_ORTHOCNTX)= ortho_cntx_
  lax_status(LAX_STATUS_WORLDCNTX)= world_cntx_
  IF( do_distr_diag_inside_bgrp_ ) THEN
     lax_status(LAX_STATUS_DISTDIAG)= 1
  ELSE
     lax_status(LAX_STATUS_DISTDIAG)= 2
  END IF
END SUBROUTINE

!----------------------------------------------------------------------------

SUBROUTINE laxlib_start_drv( ndiag_, my_world_comm, parent_comm, do_distr_diag_inside_bgrp_  )
    !
    use laxlib_processors_grid
    USE laxlib_parallel_include
    !
    !
    ! ... Ortho/diag/linear algebra group initialization
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(INOUT) :: ndiag_  ! (IN) input number of procs in the diag group, (OUT) actual number
    INTEGER, INTENT(IN) :: my_world_comm ! parallel communicator of the "local" world
    INTEGER, INTENT(IN) :: parent_comm ! parallel communicator inside which the distributed linear algebra group
                                       ! communicators are created
    LOGICAL, INTENT(IN) :: do_distr_diag_inside_bgrp_  ! comme son nom l'indique
    !
    INTEGER :: mpime      =  0  ! the global MPI task index (used in clocks) can be set with a laxlib_rank call
    !
    INTEGER :: nproc_ortho_try
    INTEGER :: parent_nproc ! nproc of the parent group
    INTEGER :: world_nproc  ! nproc of the world group
    INTEGER :: my_parent_id ! id of the parent communicator 
    INTEGER :: nparent_comm ! mumber of parent communicators
    INTEGER :: ierr = 0
    !
    IF( lax_is_initialized ) &
       CALL laxlib_end_drv ( ) 

    world_nproc  = laxlib_size( my_world_comm ) ! the global number of processors in world_comm
    mpime        = laxlib_rank( my_world_comm ) ! set the global MPI task index  (used in clocks)
    parent_nproc = laxlib_size( parent_comm )! the number of processors in the current parent communicator
    my_parent_id = mpime / parent_nproc  ! set the index of the current parent communicator
    nparent_comm = world_nproc/parent_nproc ! number of paren communicators

    ! save input value inside the module
    do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp_ 

    !
#if defined __SCALAPACK
    np_blacs     = laxlib_size( my_world_comm )
    me_blacs     = laxlib_rank( my_world_comm )
    !
    ! define a 1D grid containing all MPI tasks of the global communicator
    ! NOTE: world_cntx has the MPI communicator on entry and the BLACS context on exit
    !       BLACS_GRIDINIT() will create a copy of the communicator, which can be
    !       later retrieved using CALL BLACS_GET(world_cntx, 10, comm_copy)
    !
    world_cntx = my_world_comm 
    CALL BLACS_GRIDINIT( world_cntx, 'Row', 1, np_blacs )
    !
#endif
    !
    IF( ndiag_ > 0 ) THEN
       ! command-line argument -ndiag N or -northo N set to a value N
       ! use the command line value ensuring that it falls in the proper range
       nproc_ortho_try = MIN( ndiag_ , parent_nproc )
    ELSE 
       ! no command-line argument -ndiag N or -northo N is present
       ! insert here custom architecture specific default definitions
#if defined __SCALAPACK
       nproc_ortho_try = MAX( parent_nproc, 1 )
#else
       nproc_ortho_try = 1
#endif
    END IF
    !
    ! the ortho group for parallel linear algebra is a sub-group of the pool,
    ! then there are as many ortho groups as pools.
    !
    CALL init_ortho_group ( nproc_ortho_try, my_world_comm, parent_comm, nparent_comm, my_parent_id )
    !
    ! set the number of processors in the diag group to the actual number used
    !
    ndiag_ = nproc_ortho
    !
    lax_is_initialized = .true.
    !  
    RETURN
    !
CONTAINS

  SUBROUTINE init_ortho_group ( nproc_try_in, my_world_comm, comm_all, nparent_comm, my_parent_id )
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nproc_try_in, comm_all
    INTEGER, INTENT(IN) :: my_world_comm ! parallel communicator of the "local" world
    INTEGER, INTENT(IN) :: nparent_comm
    INTEGER, INTENT(IN) :: my_parent_id ! id of the parent communicator 

    INTEGER :: ierr, color, key, me_all, nproc_all, nproc_try

#if defined __SCALAPACK
    INTEGER, ALLOCATABLE :: blacsmap(:,:)
    INTEGER, ALLOCATABLE :: ortho_cntx_pe(:)
    INTEGER :: nprow, npcol, myrow, mycol, i, j, k
    INTEGER, EXTERNAL :: BLACS_PNUM
#endif

#if defined __MPI

    me_all    = laxlib_rank( comm_all )
    !
    nproc_all = laxlib_size( comm_all )
    !
    nproc_try = MIN( nproc_try_in, nproc_all )
    nproc_try = MAX( nproc_try, 1 )

    !  find the square closer (but lower) to nproc_try
    !
    CALL grid2d_dims( 'S', nproc_try, np_ortho(1), np_ortho(2) )
    !
    !  now, and only now, it is possible to define the number of tasks
    !  in the ortho group for parallel linear algebra
    !
    nproc_ortho = np_ortho(1) * np_ortho(2)
    !
    IF( nproc_all >= 4*nproc_ortho ) THEN
       !
       !  here we choose a processor every 4, in order not to stress memory BW
       !  on multi core procs, for which further performance enhancements are
       !  possible using OpenMP BLAS inside regter/cegter/rdiaghg/cdiaghg
       !  (to be implemented)
       !
       color = 0
       IF( me_all < 4*nproc_ortho .AND. MOD( me_all, 4 ) == 0 ) color = 1
       !
       leg_ortho = 4
       !
    ELSE IF( nproc_all >= 2*nproc_ortho ) THEN
       !
       !  here we choose a processor every 2, in order not to stress memory BW
       !
       color = 0
       IF( me_all < 2*nproc_ortho .AND. MOD( me_all, 2 ) == 0 ) color = 1
       !
       leg_ortho = 2
       !
    ELSE
       !
       !  here we choose the first processors
       !
       color = 0
       IF( me_all < nproc_ortho ) color = 1
       !
       leg_ortho = 1
       !
    END IF
    !
    key   = me_all
    !
    !  initialize the communicator for the new group by splitting the input communicator
    !
    CALL laxlib_comm_split ( comm_all, color, key, ortho_comm )
    !
    ! and remember where it comes from
    !
    ortho_parent_comm = comm_all
    !
    !  Computes coordinates of the processors, in row maior order
    !
    me_ortho1   = laxlib_rank( ortho_comm )
    !
    IF( me_all == 0 .AND. me_ortho1 /= 0 ) &
         CALL lax_error__( " init_ortho_group ", " wrong root task in ortho group ", ierr )
    !
    if( color == 1 ) then
       ortho_comm_id = 1
       CALL GRID2D_COORDS( 'R', me_ortho1, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
       CALL GRID2D_RANK( 'R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr )
       IF( ierr /= me_ortho1 ) &
            CALL lax_error__( " init_ortho_group ", " wrong task coordinates in ortho group ", ierr )
       IF( me_ortho1*leg_ortho /= me_all ) &
            CALL lax_error__( " init_ortho_group ", " wrong rank assignment in ortho group ", ierr )

       CALL laxlib_comm_split( ortho_comm, me_ortho(2), me_ortho(1), ortho_col_comm)
       CALL laxlib_comm_split( ortho_comm, me_ortho(1), me_ortho(2), ortho_row_comm)

    else
       ortho_comm_id = 0
       me_ortho(1) = me_ortho1
       me_ortho(2) = me_ortho1
    endif

#if defined __SCALAPACK
    !
    !  This part is used to eliminate the image dependency from ortho groups
    !  SCALAPACK is now independent from whatever level of parallelization
    !  is present on top of pool parallelization
    !
    ALLOCATE( ortho_cntx_pe( nparent_comm ) )
    ALLOCATE( blacsmap( np_ortho(1), np_ortho(2) ) )

    DO j = 1, nparent_comm

         CALL BLACS_GET(world_cntx, 10, ortho_cntx_pe( j ) ) ! retrieve communicator of world context
         blacsmap = 0
         nprow = np_ortho(1)
         npcol = np_ortho(2)

         IF( ( j == ( my_parent_id + 1 ) ) .and. ( ortho_comm_id > 0 ) ) THEN

           blacsmap( me_ortho(1) + 1, me_ortho(2) + 1 ) = BLACS_PNUM( world_cntx, 0, me_blacs )

         END IF

         ! All MPI tasks defined in the global communicator take part in the definition of the BLACS grid

         CALL MPI_ALLREDUCE( MPI_IN_PLACE, blacsmap,  SIZE(blacsmap), MPI_INTEGER, MPI_SUM, my_world_comm, ierr )
         IF( ierr /= 0 ) &
            CALL lax_error__( ' init_ortho_group ', ' problem in MPI_ALLREDUCE of blacsmap ', ierr )

         CALL BLACS_GRIDMAP( ortho_cntx_pe( j ), blacsmap, nprow, nprow, npcol )

         CALL BLACS_GRIDINFO( ortho_cntx_pe( j ), nprow, npcol, myrow, mycol )

         IF( ( j == ( my_parent_id + 1 ) ) .and. ( ortho_comm_id > 0 ) ) THEN

            IF(  np_ortho(1) /= nprow ) &
               CALL lax_error__( ' init_ortho_group ', ' problem with SCALAPACK, wrong no. of task rows ', 1 )
            IF(  np_ortho(2) /= npcol ) &
               CALL lax_error__( ' init_ortho_group ', ' problem with SCALAPACK, wrong no. of task columns ', 1 )
            IF(  me_ortho(1) /= myrow ) &
               CALL lax_error__( ' init_ortho_group ', ' problem with SCALAPACK, wrong task row ID ', 1 )
            IF(  me_ortho(2) /= mycol ) &
               CALL lax_error__( ' init_ortho_group ', ' problem with SCALAPACK, wrong task columns ID ', 1 )

            ortho_cntx = ortho_cntx_pe( j )

         END IF

    END DO 

    DEALLOCATE( blacsmap )
    DEALLOCATE( ortho_cntx_pe )

    !  end SCALAPACK code block

#endif

#else

    ortho_comm_id = 1

#endif

    RETURN
  END SUBROUTINE init_ortho_group

END SUBROUTINE laxlib_start_drv

!------------------------------------------------------------------------------!

SUBROUTINE print_lambda_x( lambda, idesc, n, nshow, nudx, ccc, ionode, iunit )
    IMPLICIT NONE
    include 'laxlib_low.fh'
    include 'laxlib_kinds.fh'
    real(DP), intent(in) :: lambda(:,:,:), ccc
    INTEGER, INTENT(IN) :: idesc(:,:)
    integer, intent(in) :: n, nshow, nudx
    logical, intent(in) :: ionode
    integer, intent(in) :: iunit
    !
    integer :: nnn, j, i, is
    real(DP), allocatable :: lambda_repl(:,:)
    nnn = min( nudx, nshow )
    ALLOCATE( lambda_repl( nudx, nudx ) )
    IF( ionode ) WRITE( iunit,*)
    DO is = 1, SIZE( lambda, 3 )
       CALL collect_lambda( lambda_repl, lambda(:,:,is), idesc(:,is) )
       IF( ionode ) THEN
          WRITE( iunit,3370) '    lambda   nudx, spin = ', nudx, is
          IF( nnn < n ) WRITE( iunit,3370) '    print only first ', nnn
          DO i=1,nnn
             WRITE( iunit,3380) (lambda_repl(i,j)*ccc,j=1,nnn)
          END DO
       END IF
    END DO
    DEALLOCATE( lambda_repl )
3370   FORMAT(26x,a,2i4)
3380   FORMAT(9f8.4)
    RETURN
END SUBROUTINE print_lambda_x

  SUBROUTINE laxlib_desc_init1( nsiz, nx, la_proc, idesc, rank_ip, idesc_ip )
     !
     IMPLICIT NONE
     include 'laxlib_low.fh'
     include 'laxlib_param.fh'
     include 'laxlib_kinds.fh'
     !
     INTEGER, INTENT(IN)  :: nsiz
     INTEGER, INTENT(OUT) :: nx
     LOGICAL, INTENT(OUT) :: la_proc
     INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
     INTEGER, INTENT(OUT), ALLOCATABLE :: rank_ip( :, : )
     INTEGER, INTENT(OUT), ALLOCATABLE :: idesc_ip(:,:,:)
     !
     INTEGER :: ortho_comm, np_ortho(2), me_ortho(2), ortho_comm_id, &
          leg_ortho
     !
     CALL laxlib_getval( np_ortho = np_ortho, me_ortho = me_ortho, &
          ortho_comm = ortho_comm, leg_ortho = leg_ortho, &
          ortho_comm_id = ortho_comm_id )
     !
     IF ( .NOT. ALLOCATED (idesc_ip) ) THEN
        ALLOCATE( idesc_ip( LAX_DESC_SIZE, np_ortho(1), np_ortho(2) ) )
     ELSE
        IF ( SIZE (idesc_ip,2) /= np_ortho(1) .OR. &
             SIZE (idesc_ip,3) /= np_ortho(2) ) &
             CALL lax_error__( " desc_init ", " inconsistent dimension ", 2 )
     END IF
     IF ( .NOT. ALLOCATED (rank_ip) ) &
          ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
     !
     CALL laxlib_init_desc( idesc, idesc_ip, rank_ip, nsiz, nsiz )
     ! 
     nx = idesc(LAX_DESC_NRCX)
     !
     la_proc = .FALSE.
     IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) la_proc = .TRUE.
     !
     RETURN
   END SUBROUTINE laxlib_desc_init1
   !
   SUBROUTINE laxlib_desc_init2( nsiz, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
     !
     IMPLICIT NONE
     include 'laxlib_low.fh'
     include 'laxlib_param.fh'
     include 'laxlib_kinds.fh'
     !
     INTEGER, INTENT(IN)  :: nsiz
     INTEGER, INTENT(OUT) :: nx
     LOGICAL, INTENT(OUT) :: la_proc
     INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
     INTEGER, INTENT(OUT), ALLOCATABLE :: rank_ip(:,:)
     INTEGER, INTENT(OUT), ALLOCATABLE :: irc_ip(:)
     INTEGER, INTENT(OUT), ALLOCATABLE :: nrc_ip(:)

     INTEGER :: i, j, rank
     INTEGER :: ortho_comm, np_ortho(2), me_ortho(2), ortho_comm_id, &
          leg_ortho, ortho_cntx
     !
     CALL laxlib_getval( np_ortho = np_ortho, me_ortho = me_ortho, &
          ortho_comm = ortho_comm, leg_ortho = leg_ortho, &
          ortho_comm_id = ortho_comm_id, ortho_cntx = ortho_cntx )
     !
     CALL laxlib_init_desc( idesc, nsiz, nsiz, np_ortho, me_ortho, &
          ortho_comm, ortho_cntx, ortho_comm_id )
     !
     nx = idesc(LAX_DESC_NRCX)
     !
     IF ( .NOT. ALLOCATED (rank_ip) ) THEN
        ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
        ALLOCATE( irc_ip( np_ortho(1) ), nrc_ip (np_ortho(1) ) )
     ELSE
        IF ( SIZE (rank_ip,1) /= np_ortho(1) .OR. &
             SIZE (rank_ip,2) /= np_ortho(2) ) &
             CALL lax_error__( " desc_init ", " inconsistent dimension ", 1 )
     END IF
     DO j = 0, idesc(LAX_DESC_NPC) - 1
        CALL laxlib_local_dims( irc_ip( j + 1 ), nrc_ip( j + 1 ), &
             idesc(LAX_DESC_N), idesc(LAX_DESC_NX), np_ortho(1), j )
        DO i = 0, idesc(LAX_DESC_NPR) - 1
           CALL GRID2D_RANK( 'R', idesc(LAX_DESC_NPR), idesc(LAX_DESC_NPC), i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) la_proc = .TRUE.
     !
     RETURN
   END SUBROUTINE laxlib_desc_init2
  !
  !
SUBROUTINE laxlib_init_desc_x( idesc, n, nx, np, me, comm, cntx, comm_id )
    USE laxlib_descriptor,       ONLY: la_descriptor, descla_init, laxlib_desc_to_intarray
    IMPLICIT NONE
    include 'laxlib_param.fh'
    INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
    INTEGER, INTENT(IN)  :: n   !  the size of this matrix
    INTEGER, INTENT(IN)  :: nx  !  the max among different matrixes sharing this descriptor or the same data distribution
    INTEGER, INTENT(IN)  :: np(2), me(2), comm, cntx
    INTEGER, INTENT(IN)  :: comm_id
    !
    TYPE(la_descriptor) :: descla
    !
    CALL descla_init( descla, n, nx, np, me, comm, cntx, comm_id )
    CALL laxlib_desc_to_intarray( idesc, descla )
    RETURN
END SUBROUTINE laxlib_init_desc_x


SUBROUTINE laxlib_multi_init_desc_x( idesc, idesc_ip, rank_ip, n, nx  )
    USE laxlib_descriptor,       ONLY: la_descriptor, descla_init, laxlib_desc_to_intarray
    use laxlib_processors_grid,  ONLY: leg_ortho, np_ortho, me_ortho, ortho_comm, ortho_comm_id, ortho_cntx
    IMPLICIT NONE
    include 'laxlib_param.fh'
    INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
    INTEGER, INTENT(OUT) :: idesc_ip(:,:,:)
    INTEGER, INTENT(OUT) :: rank_ip(:,:)
    INTEGER, INTENT(IN)  :: n   !  the size of this matrix
    INTEGER, INTENT(IN)  :: nx  !  the max among different matrixes sharing this descriptor or the same data distribution

    INTEGER :: i, j, rank, includeme
    INTEGER :: coor_ip( 2 )
    !
    TYPE(la_descriptor) :: descla
    !
    CALL descla_init( descla, n, nx, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
    !
    CALL laxlib_desc_to_intarray( idesc, descla )
    !
    includeme = 1
    !
    DO j = 0, idesc(LAX_DESC_NPC) - 1
       DO i = 0, idesc(LAX_DESC_NPR) - 1
          coor_ip( 1 ) = i
          coor_ip( 2 ) = j
          CALL descla_init( descla, idesc(LAX_DESC_N), idesc(LAX_DESC_NX), &
               np_ortho, coor_ip, ortho_comm, ortho_cntx, includeme )
          CALL laxlib_desc_to_intarray( idesc_ip(:,i+1,j+1), descla )
          CALL GRID2D_RANK( 'R', idesc(LAX_DESC_NPR), idesc(LAX_DESC_NPC), i, j, rank )
          rank_ip( i+1, j+1 ) = rank * leg_ortho
       END DO
    END DO
    !
    RETURN
END SUBROUTINE laxlib_multi_init_desc_x

   SUBROUTINE descla_local_dims( i2g, nl, n, nx, np, me )
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: i2g  !  global index of the first local element
      INTEGER, INTENT(OUT) :: nl   !  local number of elements
      INTEGER, INTENT(IN)  :: n    !  number of actual element in the global array
      INTEGER, INTENT(IN)  :: nx   !  dimension of the global array (nx>=n) to be distributed
      INTEGER, INTENT(IN)  :: np   !  number of processors
      INTEGER, INTENT(IN)  :: me   !  taskid for which i2g and nl are computed
      !
      !  note that we can distribute a global array larger than the
      !  number of actual elements. This could be required for performance
      !  reasons, and to have an equal partition of matrix having different size
      !  like matrixes of spin-up and spin-down
      !
      INTEGER, EXTERNAL ::  ldim_block, ldim_cyclic, ldim_block_sca
      INTEGER, EXTERNAL ::  gind_block, gind_cyclic, gind_block_sca
      !
#if __SCALAPACK
      nl  = ldim_block_sca( nx, np, me )
      i2g = gind_block_sca( 1, nx, np, me )
#else
      nl  = ldim_block( nx, np, me )
      i2g = gind_block( 1, nx, np, me )
#endif
      ! This is to try to keep a matrix N * N into the same
      ! distribution of a matrix NX * NX, useful to have
      ! the matrix of spin-up distributed in the same way
      ! of the matrix of spin-down
      !
      IF( i2g + nl - 1 > n ) nl = n - i2g + 1
      IF( nl < 0 ) nl = 0
      RETURN
      !
   END SUBROUTINE descla_local_dims


!   ----------------------------------------------
!   Simplified driver 

   SUBROUTINE diagonalize_parallel_x( n, rhos, rhod, s, idesc )

      USE dspev_module

      IMPLICIT NONE
      include 'laxlib_kinds.fh'
      include 'laxlib_param.fh'
      include 'laxlib_mid.fh'
      include 'laxlib_low.fh'
      REAL(DP), INTENT(IN)  :: rhos(:,:) !  input symmetric matrix
      REAL(DP)              :: rhod(:)   !  output eigenvalues
      REAL(DP)              :: s(:,:)    !  output eigenvectors
      INTEGER,  INTENT(IN) :: n         !  size of the global matrix
      INTEGER,  INTENT(IN) :: idesc(LAX_DESC_SIZE)

      IF( n < 1 ) RETURN

      !  Matrix is distributed on the same processors group
      !  used for parallel matrix multiplication
      !
      IF( SIZE(s,1) /= SIZE(rhos,1) .OR. SIZE(s,2) /= SIZE(rhos,2) ) &
         CALL lax_error__( " diagonalize_parallel ", " inconsistent dimension for s and rhos ", 1 )

      IF ( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         !
         IF( SIZE(s,1) /= idesc(LAX_DESC_NRCX) ) &
            CALL lax_error__( " diagonalize_parallel ", " inconsistent dimension ", 1)
         !
         !  Compute local dimension of the cyclically distributed matrix
         !
         s = rhos
         !
#if defined(__SCALAPACK)
         CALL pdsyevd_drv( .true. , n, idesc(LAX_DESC_NRCX), s, SIZE(s,1), rhod, idesc(LAX_DESC_CNTX), idesc(LAX_DESC_COMM) )
#else
         CALL laxlib_pdsyevd( .true., n, idesc, s, SIZE(s,1), rhod )
#endif
         !
      END IF

      RETURN

   END SUBROUTINE diagonalize_parallel_x


   SUBROUTINE diagonalize_serial_x( n, rhos, rhod )
      IMPLICIT NONE
      include 'laxlib_kinds.fh'
      include 'laxlib_low.fh'
      INTEGER,  INTENT(IN)  :: n
      REAL(DP)              :: rhos(:,:)
      REAL(DP)              :: rhod(:)
      !
      ! inputs:
      ! n     size of the eigenproblem
      ! rhos  the symmetric matrix
      ! outputs:
      ! rhos  eigenvectors
      ! rhod  eigenvalues
      !
      REAL(DP), ALLOCATABLE :: aux(:)
      INTEGER :: i, j, k

      IF( n < 1 ) RETURN

      ALLOCATE( aux( n * ( n + 1 ) / 2 ) )

      !  pack lower triangle of rho into aux
      !
      k = 0
      DO j = 1, n
         DO i = j, n
            k = k + 1
            aux( k ) = rhos( i, j )
         END DO
      END DO

      CALL dspev_drv( 'V', 'L', n, aux, rhod, rhos, SIZE(rhos,1) )

      DEALLOCATE( aux )

      RETURN
   END SUBROUTINE diagonalize_serial_x

   SUBROUTINE diagonalize_serial_gpu( m, rhos, rhod, s, info )
#if defined(__CUDA)
      use cudafor
#if defined ( __USE_CUSOLVER )
      USE cusolverDn
#else
      use eigsolve_vars
      use nvtx_inters
      use dsyevd_gpu
#endif
      IMPLICIT NONE
      include 'laxlib_kinds.fh'
      INTEGER, INTENT(IN) :: m
      REAL(DP), DEVICE, INTENT(IN) :: rhos(:,:)
      REAL(DP), DEVICE, INTENT(OUT) :: rhod(:)
      REAL(DP), DEVICE, INTENT(OUT) :: s(:,:)
      INTEGER, INTENT(OUT) :: info
      !
      INTEGER :: lwork_d
      INTEGER :: i, j, lda
      !
#if defined (__USE_CUSOLVER)
      !
      INTEGER, DEVICE        :: devInfo
      TYPE(cusolverDnHandle) :: cuSolverHandle
      REAL(DP), ALLOCATABLE, DEVICE :: work_d(:)
      !
#else
      !
      REAL(DP), ALLOCATABLE :: work_d(:), a(:,:)
      ATTRIBUTES( DEVICE ) :: work_d, a
      REAL(DP), ALLOCATABLE :: b(:,:)
      REAL(DP), ALLOCATABLE :: work_h(:), w_h(:), z_h(:,:)
      ATTRIBUTES( PINNED ) :: work_h, w_h, z_h
      INTEGER, ALLOCATABLE :: iwork_h(:)
      ATTRIBUTES( PINNED ) :: iwork_h
      !
      INTEGER :: lwork_h, liwork_h
      !
#endif
      ! .... Subroutine Body
      !
#if defined (__USE_CUSOLVER)
      !
      s = rhos
      lda = SIZE( rhos, 1 )
      !
      info = cusolverDnCreate(cuSolverHandle)
      IF ( info /= CUSOLVER_STATUS_SUCCESS ) &
         CALL lax_error__( ' diagonalize_serial_gpu ', 'cusolverDnCreate',  ABS( info ) )

      info = cusolverDnDsyevd_bufferSize( &
             cuSolverHandle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, m, s, lda, rhod, lwork_d)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' error in solver 1 ', ABS( info ) )

      ALLOCATE( work_d ( lwork_d ), STAT=info )
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' allocate work_d ', ABS( info ) )

      info = cusolverDnDsyevd( &
             cuSolverHandle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, m, s, lda, rhod, work_d, lwork_d, devInfo)
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' error in solver 2 ', ABS( info ) )

      info = cudaDeviceSynchronize()

      info = cusolverDnDestroy(cuSolverHandle)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' diagonalize_serial_gpu ', ' cusolverDnDestroy failed ', ABS( info ) )

      DEALLOCATE( work_d )

      !
#else
      !
      info = 0
      lwork_d  = 2*64*64 + 66*SIZE(rhos,1)
      lwork_h = 1 + 6*SIZE(rhos,1) + 2*SIZE(rhos,1)*SIZE(rhos,1)
      liwork_h = 3 + 5*SIZE(rhos,1)
      ALLOCATE(work_d(lwork_d),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' allocate work_d ', ABS( info ) )
      ALLOCATE(a(SIZE(rhos,1),SIZE(rhos,2)),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' allocate a ', ABS( info ) )
      ALLOCATE(work_h(lwork_h),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' allocate work_h ', ABS( info ) )
      ALLOCATE(iwork_h(liwork_h),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' allocate iwork_h ', ABS( info ) )
      !
      ALLOCATE(w_h(SIZE(rhod)),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' allocate w_h ', ABS( info ) )
      ALLOCATE(z_h(SIZE(s,1),SIZE(s,2)),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' allocate z_h ', ABS( info ) )

      if(initialized == 0) call init_eigsolve_gpu
      
      info = cudaMemcpy(a, rhos, SIZE(rhos,1)*SIZE(rhos,2), cudaMemcpyDeviceToDevice)
      lda = SIZE(rhos,1)
      !$cuf kernel do(2) <<<*,*>>>
      do j = 1,m
        do i = 1,m
          if (i > j) then
            s(i,j) = a(i,j)
          endif
        end do
      end do
      
      call dsyevd_gpu('V', 'U', 1, m, m, a, lda, s, lda, rhod, work_d, lwork_d, &
                      work_h, lwork_h, iwork_h, liwork_h, z_h, lda, w_h, info)

      DEALLOCATE(z_h)
      DEALLOCATE(w_h)
      DEALLOCATE(iwork_h)
      DEALLOCATE(work_h)
      DEALLOCATE(a)
      DEALLOCATE(work_d)
#endif
#else
      CALL lax_error__( ' laxlib diagonalize_serial_gpu ', ' not compiled in this version ', 0 )
#endif
   END SUBROUTINE

