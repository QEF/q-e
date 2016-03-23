!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_diag
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_size, mp_rank, mp_sum, mp_comm_free, mp_comm_split
  !
  ! The following variables are needed in order to set up the communicator
  ! for scalapack
  !
  USE mp_world, ONLY : world_comm
  !
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... linear-algebra group (also known as "ortho" or "diag" group). 
  ! ... Used for parallelization of dense-matrix diagonalization used in
  ! ... iterative diagonalization/orthonormalization, matrix-matrix products
  !
  INTEGER :: np_ortho(2) = 1  ! size of the processor grid used in ortho
  INTEGER :: me_ortho(2) = 0  ! coordinates of the processors
  INTEGER :: me_ortho1   = 0  ! task id for the ortho group
  INTEGER :: nproc_ortho = 1  ! size of the ortho group:
  INTEGER :: leg_ortho   = 1  ! the distance in the father communicator
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
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_diag( ndiag_, parent_comm, nparent_comm, my_parent_id )
    !---------------------------------------------------------------------------
    !
    ! ... Ortho/diag/linear algebra group initialization
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ndiag_, parent_comm
    INTEGER, INTENT(IN) :: nparent_comm ! number of parent communicators
    INTEGER, INTENT(IN) :: my_parent_id ! id of the parent communicator 
    !
    INTEGER :: nproc_ortho_try
    INTEGER :: parent_nproc  ! nproc of the parent group
    INTEGER :: ierr = 0
    !
    parent_nproc = mp_size( parent_comm )
    !
#if defined __SCALAPACK
    np_blacs     = mp_size( world_comm )
    me_blacs     = mp_rank( world_comm )
    !
    ! define a 1D grid containing all MPI tasks of the global communicator
    ! NOTE: world_cntx has the MPI communicator on entry and the BLACS context on exit
    !       BLACS_GRID_INIT() will create a copy of the communicator, which can be
    !       later retrieved using CALL BLACS_GET(world_cntx, 10, comm_copy)
    !
    world_cntx = world_comm
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
       nproc_ortho_try = MAX( parent_nproc/2, 1 )
#else
       nproc_ortho_try = 1
#endif
    END IF
    !
    ! the ortho group for parallel linear algebra is a sub-group of the pool,
    ! then there are as many ortho groups as pools.
    !
    CALL init_ortho_group( nproc_ortho_try, parent_comm, nparent_comm, my_parent_id )
    !  
    RETURN
    !
  END SUBROUTINE mp_start_diag
  !
  !
  SUBROUTINE init_ortho_group( nproc_try_in, comm_all, nparent_comm, my_parent_id )
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nproc_try_in, comm_all
    INTEGER, INTENT(IN) :: nparent_comm
    INTEGER, INTENT(IN) :: my_parent_id ! id of the parent communicator 

    LOGICAL, SAVE :: first = .true.
    INTEGER :: ierr, color, key, me_all, nproc_all, nproc_try

#if defined __SCALAPACK
    INTEGER, ALLOCATABLE :: blacsmap(:,:)
    INTEGER, ALLOCATABLE :: ortho_cntx_pe(:)
    INTEGER :: nprow, npcol, myrow, mycol, i, j, k
    INTEGER, EXTERNAL :: BLACS_PNUM
#endif


#if defined __MPI

    me_all    = mp_rank( comm_all )
    !
    nproc_all = mp_size( comm_all )
    !
    nproc_try = MIN( nproc_try_in, nproc_all )
    nproc_try = MAX( nproc_try, 1 )

    IF( .NOT. first ) CALL clean_ortho_group ( ) 

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
    CALL mp_comm_split ( comm_all, color, key, ortho_comm )
    !
    ! and remember where it comes from
    !
    ortho_parent_comm = comm_all
    !
    !  Computes coordinates of the processors, in row maior order
    !
    me_ortho1   = mp_rank( ortho_comm )
    !
    IF( me_all == 0 .AND. me_ortho1 /= 0 ) &
         CALL errore( " init_ortho_group ", " wrong root task in ortho group ", ierr )
    !
    if( color == 1 ) then
       ortho_comm_id = 1
       CALL GRID2D_COORDS( 'R', me_ortho1, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
       CALL GRID2D_RANK( 'R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr )
       IF( ierr /= me_ortho1 ) &
            CALL errore( " init_ortho_group ", " wrong task coordinates in ortho group ", ierr )
       IF( me_ortho1*leg_ortho /= me_all ) &
            CALL errore( " init_ortho_group ", " wrong rank assignment in ortho group ", ierr )

       CALL mp_comm_split( ortho_comm, me_ortho(2), me_ortho(1), ortho_col_comm)
       CALL mp_comm_split( ortho_comm, me_ortho(1), me_ortho(2), ortho_row_comm)

    else
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

         CALL mp_sum( blacsmap, world_comm ) 

         CALL BLACS_GRIDMAP( ortho_cntx_pe( j ), blacsmap, nprow, nprow, npcol )

         CALL BLACS_GRIDINFO( ortho_cntx_pe( j ), nprow, npcol, myrow, mycol )

         IF( ( j == ( my_parent_id + 1 ) ) .and. ( ortho_comm_id > 0 ) ) THEN

            IF(  np_ortho(1) /= nprow ) &
               CALL errore( ' init_ortho_group ', ' problem with SCALAPACK, wrong no. of task rows ', 1 )
            IF(  np_ortho(2) /= npcol ) &
               CALL errore( ' init_ortho_group ', ' problem with SCALAPACK, wrong no. of task columns ', 1 )
            IF(  me_ortho(1) /= myrow ) &
               CALL errore( ' init_ortho_group ', ' problem with SCALAPACK, wrong task row ID ', 1 )
            IF(  me_ortho(2) /= mycol ) &
               CALL errore( ' init_ortho_group ', ' problem with SCALAPACK, wrong task columns ID ', 1 )

            ortho_cntx = ortho_cntx_pe( j )

         END IF

    END DO 

    DEALLOCATE( blacsmap )
    DEALLOCATE( ortho_cntx_pe )


#endif

#else

    ortho_comm_id = 1

#endif

    first = .false.

    RETURN
  END SUBROUTINE init_ortho_group
  !
  SUBROUTINE clean_ortho_group ( )
    !  
    !  free resources associated to the communicator
    !
    CALL mp_comm_free( ortho_comm )
    IF(  ortho_comm_id > 0  ) THEN
       CALL mp_comm_free( ortho_col_comm )
       CALL mp_comm_free( ortho_row_comm )
    ENDIF
#if defined __SCALAPACK
    IF(  ortho_cntx /= -1 ) CALL BLACS_GRIDEXIT( ortho_cntx )
    ortho_cntx = -1
#endif
    !
  END SUBROUTINE clean_ortho_group
  !
END MODULE mp_diag
