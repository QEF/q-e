!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_global
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_comm_free, mp_size, mp_rank, mp_sum, mp_barrier, &
       mp_bcast, mp_start, mp_end, mp_env
  USE io_global, ONLY : stdout, io_global_start, io_global_getmeta
  USE parallel_include
  !
  IMPLICIT NONE 
  !
  SAVE
  !
  INTEGER :: mpime = 0  ! absolute processor index starting from 0
  INTEGER :: root  = 0  ! index of the absolute root processor
  INTEGER :: nproc = 1  ! absolute number of processor
  INTEGER :: nproc_file = 1  ! absolute number of processor written in the 
  ! xml punch file
  INTEGER :: world_comm = 0  ! communicator of all processor
#if defined __SCALAPACK
  INTEGER :: me_blacs   =  0  ! BLACS processor index starting from 0
  INTEGER :: np_blacs   =  1  ! BLACS number of processor
  INTEGER :: world_cntx = -1  ! BLACS context of all processor 
#endif

  INTEGER :: kunit = 1  ! granularity of k-point distribution
  !
  ! ... indeces ( all starting from 0 !!! )
  !
  INTEGER :: me_pool     = 0  ! index of the processor within a pool 
  INTEGER :: me_image    = 0  ! index of the processor within an image
  INTEGER :: me_bgrp     = 0  ! index of the processor within a band group
  INTEGER :: root_pool   = 0  ! index of the root processor within a pool
  INTEGER :: root_image  = 0  ! index of the root processor within an image
  INTEGER :: root_bgrp   = 0  ! index of the root processor within a band group
  INTEGER :: my_pool_id  = 0  ! index of my pool
  INTEGER :: my_image_id = 0  ! index of my image
  INTEGER :: my_bgrp_id  = 0  ! index of my band group
  INTEGER :: me_ortho(2) = 0  ! coordinates of the processors
  INTEGER :: me_ortho1   = 0  ! task id for the ortho group
  INTEGER :: me_pgrp     = 0  ! task id for plane wave task group
  !
  INTEGER :: npool       = 1  ! number of "k-points"-pools
  INTEGER :: nimage      = 1  ! number of "path-images"-pools
  INTEGER :: nbgrp       = 1  ! number of band groups
  INTEGER :: nogrp       = 1  ! number of proc. in an orbital "task group" 
  INTEGER :: npgrp       = 1  ! number of proc. in a plane-wave "task group" 
  INTEGER :: nproc_pool  = 1  ! number of processor within a pool
  INTEGER :: nproc_pool_file  = 1  ! number of processor within a pool of
  !   written in the xml punch file
  INTEGER :: nproc_image = 1  ! number of processor within an image
  INTEGER :: nproc_image_file  = 1  ! number of processor within a image
  INTEGER :: nproc_bgrp  = 1  ! number of processor within a band group
  INTEGER :: np_ortho(2) = 1  ! size of the processor grid used in ortho
  INTEGER :: nproc_ortho = 1  ! size of the ortho group:
  INTEGER :: leg_ortho   = 1  ! the distance in the father communicator
  ! of two neighbour processors in ortho_comm
  INTEGER, ALLOCATABLE :: nolist(:) ! list of processors in my orbital task group 
  INTEGER, ALLOCATABLE :: nplist(:) ! list of processors in my plane wave task group 
  !
  ! ... communicators
  !
  INTEGER :: inter_pool_comm  = 0  ! inter pool communicator
  INTEGER :: intra_pool_comm  = 0  ! intra pool communicator
  INTEGER :: inter_image_comm = 0  ! inter image communicator
  INTEGER :: intra_image_comm = 0  ! intra image communicator  
  INTEGER :: inter_bgrp_comm  = 0  ! inter band group communicator
  INTEGER :: intra_bgrp_comm  = 0  ! intra band group communicator  
  INTEGER :: pgrp_comm        = 0  ! plane-wave group communicator (task grouping)
  INTEGER :: ogrp_comm        = 0  ! orbital group communicarot (task grouping)
  INTEGER :: ortho_comm       = 0  ! communicator used for fast and memory saving ortho
  INTEGER :: ortho_comm_id    = 0  ! id of the ortho_comm
#if defined __SCALAPACK
  INTEGER :: ortho_cntx       = -1 ! BLACS context for ortho_comm
#endif
  !
  ! ... Task Groups parallelization
  !
  LOGICAL :: &
    use_task_groups = .FALSE.  ! if TRUE task groups parallelization is used
  !
  ! ... Module Private stuff
  !
  LOGICAL, PRIVATE :: user_nproc_ortho = .FALSE.
  !
  PRIVATE :: init_pool
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_startup ( ) 
    !-----------------------------------------------------------------------
    ! ... This subroutine initializes MPI
    ! ... Processes are organized in NIMAGE images each dealing with a subset of
    ! ... images used to discretize the "path" (only in "path" optimizations)
    ! ... Within each image processes are organized in NPOOL pools each dealing
    ! ... with a subset of kpoints.
    ! ... Within each pool R & G space distribution is performed.
    ! ... NPROC is read from command line or can be set with the appropriate
    ! ... environment variable ( for example use 'setenv MP_PROCS 8' on IBM SP
    ! ... machine to run on NPROC=8 processors ); NIMAGE and NPOOL are read from
    ! ... command line.
    ! ... NPOOL must be a whole divisor of NPROC
    !
    IMPLICIT NONE
    INTEGER :: world, ntask_groups, nproc_ortho_in, meta_ionode_id 
    INTEGER :: root = 0
    LOGICAL :: meta_ionode
    !
    !
    CALL mp_start()
    !
    ! ... get the basic parameters from communications sub-system
    ! ... to handle processors
    ! ... mpime = processor number, starting from 0
    ! ... nproc = number of processors
    ! ... world = group index of all processors
    !
    CALL mp_env( nproc, mpime, world )
    !
    !
    ! ... now initialize processors and groups variables
    ! ... set global coordinate for this processor
    ! ... root  = index of the root processor
    !
    CALL mp_global_start( root, mpime, world, nproc )
    !
    ! ... initialize input/output, set the I/O node
    !
    CALL io_global_start( mpime, root )
    !
    ! ... get the "meta" I/O node
    !
    CALL io_global_getmeta ( meta_ionode, meta_ionode_id )
    !
    IF ( meta_ionode ) THEN
       !
       ! ... How many parallel images ?
       !
       CALL get_arg_nimage( nimage )
       !
       nimage = MAX( nimage, 1 )
       nimage = MIN( nimage, nproc )
       !
       ! ... How many parallel images ?
       !
       CALL get_arg_nbgrp( nbgrp )
       !
       nbgrp = MAX( nbgrp, 1 )
       nbgrp = MIN( nbgrp, nproc )
       !
       ! ... How many pools ?
       !
       CALL get_arg_npool( npool )
       !
       npool = MAX( npool, 1 )
       npool = MIN( npool, nproc )
       !
       ! ... How many task groups ?
       !
       CALL get_arg_ntg( ntask_groups )
       !
       ! ... How many processors involved in diagonalization of Hamiltonian ?
       !
       CALL get_arg_northo( nproc_ortho_in )
       !
       IF( nproc_ortho_in < 1 ) THEN
          !  any invalid value means use the default
          user_nproc_ortho = .FALSE.
       ELSE
          user_nproc_ortho = .TRUE.
       END IF
       !
       nproc_ortho_in = MAX( nproc_ortho_in, 1 )
       nproc_ortho_in = MIN( nproc_ortho_in, nproc )
       !
    END IF
    !
    CALL mp_barrier()
    !
    ! ... transmit npool and nimage
    !
    CALL mp_bcast( npool,  meta_ionode_id )
    CALL mp_bcast( nimage, meta_ionode_id )
    CALL mp_bcast( nbgrp, meta_ionode_id )
    CALL mp_bcast( ntask_groups, meta_ionode_id )
    CALL mp_bcast( nproc_ortho_in, meta_ionode_id )
    CALL mp_bcast( user_nproc_ortho, meta_ionode_id )
    !
    use_task_groups = ( ntask_groups > 1 )
    !
    ! ... all pools are initialized here
    !
    CALL init_pool( nimage, ntask_groups, nproc_ortho_in )
    !
    RETURN
    !
  END SUBROUTINE mp_startup
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_global_start( root_i, mpime_i, group_i, nproc_i )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: root_i, mpime_i, group_i, nproc_i
    !
    root             = root_i
    mpime            = mpime_i
    world_comm       = group_i
    nproc            = nproc_i
    nproc_pool       = nproc_i
    nproc_image      = nproc_i
    nproc_bgrp       = nproc_i
    my_pool_id       = 0
    my_image_id      = 0
    my_bgrp_id       = 0
    me_pool          = mpime
    me_image         = mpime
    me_bgrp          = mpime
    me_pgrp          = me_pool
    root_pool        = root
    root_image       = root
    root_bgrp        = root
    inter_pool_comm  = group_i
    intra_pool_comm  = group_i
    inter_image_comm = group_i
    intra_image_comm = group_i
    inter_bgrp_comm  = group_i
    intra_bgrp_comm  = group_i
    ortho_comm       = group_i
    ALLOCATE( nolist( nproc_i ) )
    ALLOCATE( nplist( nproc_i ) )
    nolist = 0
    nplist = 0
    !
    RETURN
    !
  END SUBROUTINE mp_global_start
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_global_end ( )
    !-----------------------------------------------------------------------
    !
    CALL mp_barrier()
    CALL mp_end ()
    IF (ALLOCATED (nolist) ) DEALLOCATE ( nolist )
    IF (ALLOCATED (nplist) ) DEALLOCATE ( nplist )
    !
  END SUBROUTINE mp_global_end
  !
  !-----------------------------------------------------------------------     
  SUBROUTINE mp_global_group_start( mep, myp, nprocp, num_of_pools )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !     
    INTEGER, INTENT(IN) :: mep, myp, nprocp, num_of_pools
    !
    me_pool    = mep
    my_pool_id = myp
    nproc_pool = nprocp
    npool      = num_of_pools
    !
    RETURN
    !
  END SUBROUTINE mp_global_group_start
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE init_pool( nimage_ , ntask_groups_ , nproc_ortho_in )
    !----------------------------------------------------------------------------
    !
    ! ... This routine initialize the pool :  MPI division in pools and images
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nimage_
    INTEGER, INTENT(IN) :: ntask_groups_
    INTEGER, INTENT(IN) :: nproc_ortho_in
    !
    INTEGER :: nproc_ortho_try
    INTEGER :: ierr = 0
    !
#if defined (__PARA)
    ! 
    !
    IF( nimage < 1 ) &
       CALL errore( 'init_pool', 'invalid number of images, less than one', 1 )
    !
    nimage = nimage_
    !  
    ! ... here we set all parallel indeces (defined in mp_global): 
    !
    !
    ! ... number of cpus per image
    !
    nproc_image = nproc / nimage
    !
    IF ( nproc < nimage ) &
       CALL errore( 'init_pool', 'invalid number of images, nimage > nproc', 1 )
    !
    IF ( MOD( nproc, nimage ) /= 0 ) &
       CALL errore( 'init_pool', 'invalid number of images, nproc /= nproc_image * nimage', 1 ) 
    !
    ! ... my_image_id  =  image index for this processor   ( 0 : nimage - 1 )
    ! ... me_image     =  processor index within the image ( 0 : nproc_image - 1 )
    !
    my_image_id = mpime / nproc_image
    me_image    = MOD( mpime, nproc_image )
    !
    CALL mp_barrier()
    !
    ! ... the intra_image_comm communicator is created
    !
    CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, my_image_id, mpime, intra_image_comm, ierr )
    !
    IF ( ierr /= 0 ) &
       CALL errore( 'init_pool', 'intra image communicator initialization', ABS(ierr) )
    !
    CALL mp_barrier()
    !
    ! ... the inter_image_comm communicator is created                     
    !     
    CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, me_image, mpime, inter_image_comm, ierr )  
    !
    IF ( ierr /= 0 ) &
       CALL errore( 'init_pool', 'inter image communicator initialization', ABS(ierr) )
    !
    ! ... Now the band group communicator
    !
    nproc_bgrp = nproc_image / nbgrp
    !
    IF ( MOD( nproc_image, nbgrp ) /= 0 ) &
         CALL errore( 'init_pool', 'invalid number of band group, nproc_image /= nproc_bgrp * nbgrp', 1 )  
    !
    ! ... my_bgrp_id  =  band group index for this processor   ( 0 : nbgrp - 1 )
    ! ... me_bgrp     =  processor index within the band group ( 0 : nproc_bgrp - 1 )
    !
    my_bgrp_id = me_image / nproc_bgrp
    me_bgrp    = MOD( me_image, nproc_bgrp )
    !
    CALL mp_barrier()
    !
    ! ... the intra_bgrp_comm communicator is created
    !
    CALL MPI_COMM_SPLIT( intra_image_comm, my_bgrp_id, me_image, intra_bgrp_comm, ierr )
    !
    IF ( ierr /= 0 ) &
       CALL errore( 'init_pool', 'intra band group communicator initialization', ABS(ierr) )
    !
    CALL mp_barrier()
    !
    ! ... the inter_bgrp_comm communicator is created                     
    !     
    CALL MPI_COMM_SPLIT( intra_image_comm, me_bgrp, me_image, inter_bgrp_comm, ierr )  
    !
    IF ( ierr /= 0 ) &
       CALL errore( 'init_pool', 'inter band group communicator initialization', ABS(ierr) )
    !
    ! ... number of cpus per pool of k-points (they are created inside each image)
    !
    nproc_pool = nproc_bgrp / npool
    !
    IF ( MOD( nproc_bgrp, npool ) /= 0 ) &
         CALL errore( 'init_pool', 'invalid number of pools, nproc_bgrp /= nproc_pool * npool', 1 )  
    !
    ! ... my_pool_id  =  pool index for this processor    ( 0 : npool - 1 )
    ! ... me_pool     =  processor index within the pool  ( 0 : nproc_pool - 1 )
    !
    my_pool_id = me_bgrp / nproc_pool    
    me_pool    = MOD( me_bgrp, nproc_pool )
    !
    CALL mp_barrier( intra_bgrp_comm )
    !
    ! ... the intra_pool_comm communicator is created
    !
    CALL MPI_COMM_SPLIT( intra_bgrp_comm, my_pool_id, me_bgrp, intra_pool_comm, ierr )
    !
    IF ( ierr /= 0 ) &
       CALL errore( 'init_pool', 'intra pool communicator initialization', ABS(ierr) )
    !
    CALL mp_barrier( intra_bgrp_comm )
    !
    ! ... the inter_pool_comm communicator is created
    !
    CALL MPI_COMM_SPLIT( intra_bgrp_comm, me_pool, me_bgrp, inter_pool_comm, ierr )
    !
    IF ( ierr /= 0 ) &
       CALL errore( 'init_pool', 'inter pool communicator initialization', ABS(ierr) )
    !
#endif
    !
    !
#if defined __SCALAPACK

    ! define a 1D grid containing all MPI task of MPI_COMM_WORLD communicator
    !
    CALL BLACS_PINFO( me_blacs, np_blacs )
    CALL BLACS_GET( -1, 0, world_cntx )
    CALL BLACS_GRIDINIT( world_cntx, 'Row', 1, np_blacs )
    !
#endif
    !
    IF( user_nproc_ortho ) THEN
       ! use the command line value ensuring that it falls in the proper range.
       nproc_ortho_try = MIN( nproc_ortho_in , nproc_pool )
       nproc_ortho_try = MAX( nproc_ortho_try , 1 )
    ELSE
       ! here we can play with custom architecture specific default definitions
#if defined __SCALAPACK
       nproc_ortho_try = MAX( nproc_pool/2, 1 )
#else
       nproc_ortho_try = 1
#endif
    END IF
    !
    ! the ortho group for parallel linear algebra is a sub-group of the pool,
    ! then there are as many ortho groups as pools.
    !
    CALL init_ortho_group( nproc_ortho_try, intra_pool_comm )
    !  
    IF( ntask_groups_ > 1 ) THEN
       nogrp = ntask_groups_
       CALL init_task_groups( )
    END IF
    !
    RETURN
    !
  END SUBROUTINE init_pool
  !
  !
  SUBROUTINE init_task_groups( )
    !
    INTEGER :: i, n1, ipos, color, key, ierr, itsk, ntsk
    INTEGER :: pgroup( nproc_pool )
    !
    !SUBDIVIDE THE PROCESSORS IN GROUPS
    !
    !THE NUMBER OF GROUPS HAS TO BE A DIVISOR OF THE NUMBER
    !OF PROCESSORS
    !
    IF( MOD( nproc_pool, nogrp ) /= 0 ) &
         CALL errore( " init_task_groups ", "the number of task groups should be a divisor of nproc_pool ", 1 )
    !
    npgrp = nproc_pool / nogrp

    DO i = 1, nproc_pool
       pgroup( i ) = i - 1
    ENDDO
    !
    !LIST OF PROCESSORS IN MY ORBITAL GROUP
    !
    !  processors in these group have contiguous indexes
    !
    N1 = ( me_pool / NOGRP ) * NOGRP - 1
    DO i = 1, nogrp
       nolist( I ) = pgroup( N1 + I + 1 )
       IF( me_pool == nolist( I ) ) ipos = i - 1
    ENDDO
    !
    !LIST OF PROCESSORS IN MY PLANE WAVE GROUP
    !
    DO I = 1, npgrp
       nplist( I ) = pgroup( ipos + ( i - 1 ) * nogrp + 1 )
    ENDDO

    !
    !SET UP THE GROUPS
    !
    !
    !CREATE ORBITAL GROUPS
    !
#if defined __MPI
    color = me_pool / nogrp
    key   = MOD( me_pool , nogrp )
    CALL MPI_COMM_SPLIT( intra_pool_comm, color, key, ogrp_comm, ierr )
    if( ierr /= 0 ) &
         CALL errore( ' init_task_groups ', ' creating ogrp_comm ', ABS(ierr) )
    CALL MPI_COMM_RANK( ogrp_comm, itsk, IERR )
    CALL MPI_COMM_SIZE( ogrp_comm, ntsk, IERR )
    IF( nogrp /= ntsk ) CALL errore( ' init_task_groups ', ' ogrp_comm size ', ntsk )
    DO i = 1, nogrp
       IF( me_pool == nolist( i ) ) THEN
          IF( (i-1) /= itsk ) CALL errore( ' init_task_groups ', ' ogrp_comm rank ', itsk )
       END IF
    END DO
#endif
    !
    !CREATE PLANEWAVE GROUPS
    !
#if defined __MPI
    color = MOD( me_pool , nogrp )
    key   = me_pool / nogrp
    CALL MPI_COMM_SPLIT( intra_pool_comm, color, key, pgrp_comm, ierr )
    if( ierr /= 0 ) &
         CALL errore( ' init_task_groups ', ' creating pgrp_comm ', ABS(ierr) )
    CALL MPI_COMM_RANK( pgrp_comm, itsk, IERR )
    CALL MPI_COMM_SIZE( pgrp_comm, ntsk, IERR )
    IF( npgrp /= ntsk ) CALL errore( ' init_task_groups ', ' pgrp_comm size ', ntsk )
    DO i = 1, npgrp
       IF( me_pool == nplist( i ) ) THEN
          IF( (i-1) /= itsk ) CALL errore( ' init_task_groups ', ' pgrp_comm rank ', itsk )
       END IF
    END DO
    me_pgrp = itsk
#endif


    RETURN
  END SUBROUTINE init_task_groups
  !
  !
  SUBROUTINE init_ortho_group( nproc_try_in, comm_all )
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nproc_try_in, comm_all

    LOGICAL, SAVE :: first = .true.
    INTEGER :: ierr, color, key, me_all, nproc_all, nproc_try

#if defined __SCALAPACK
    INTEGER, ALLOCATABLE :: blacsmap(:,:)
    INTEGER, ALLOCATABLE :: ortho_cntx_pe(:,:,:)
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

    IF( .NOT. first ) THEN
       !  
       !  free resources associated to the communicator
       !
       CALL mp_comm_free( ortho_comm )
       !
#if defined __SCALAPACK
       IF(  ortho_comm_id > 0  ) THEN
          CALL BLACS_GRIDEXIT( ortho_cntx )
       ENDIF
       ortho_cntx = -1
#endif
       !
    END IF

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
    CALL MPI_COMM_SPLIT( comm_all, color, key, ortho_comm, ierr )
    IF( ierr /= 0 ) &
         CALL errore( " init_ortho_group ", " initializing ortho group communicator ", ierr )
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
    else
       ortho_comm_id = 0
       me_ortho(1) = me_ortho1
       me_ortho(2) = me_ortho1
    endif

#if defined __SCALAPACK

    ALLOCATE( ortho_cntx_pe( npool, nbgrp, nimage ) )
    ALLOCATE( blacsmap( np_ortho(1), np_ortho(2) ) )

    DO j = 1, nimage

     DO k = 1, nbgrp

       DO i = 1, npool

         CALL BLACS_GET( -1, 0, ortho_cntx_pe( i, k, j ) ) ! take a default value 

         blacsmap = 0
         nprow = np_ortho(1)
         npcol = np_ortho(2)

         IF( ( j == ( my_image_id + 1 ) ) .and. ( k == ( my_bgrp_id + 1 ) ) .and.  &
             ( i == ( my_pool_id  + 1 ) ) .and. ( ortho_comm_id > 0 ) ) THEN

           blacsmap( me_ortho(1) + 1, me_ortho(2) + 1 ) = BLACS_PNUM( world_cntx, 0, me_blacs )

         END IF

         ! All MPI tasks defined in world comm take part in the definition of the BLACS grid

         CALL mp_sum( blacsmap ) 

         CALL BLACS_GRIDMAP( ortho_cntx_pe(i,k,j), blacsmap, nprow, nprow, npcol )

         CALL BLACS_GRIDINFO( ortho_cntx_pe(i,k,j), nprow, npcol, myrow, mycol )

         IF( ( j == ( my_image_id + 1 ) ) .and. ( k == ( my_bgrp_id + 1 ) ) .and. &
             ( i == ( my_pool_id  + 1 ) ) .and. ( ortho_comm_id > 0 ) ) THEN

            IF(  np_ortho(1) /= nprow ) &
               CALL errore( ' init_ortho_group ', ' problem with SCALAPACK, wrong no. of task rows ', 1 )
            IF(  np_ortho(2) /= npcol ) &
               CALL errore( ' init_ortho_group ', ' problem with SCALAPACK, wrong no. of task columns ', 1 )
            IF(  me_ortho(1) /= myrow ) &
               CALL errore( ' init_ortho_group ', ' problem with SCALAPACK, wrong task row ID ', 1 )
            IF(  me_ortho(2) /= mycol ) &
               CALL errore( ' init_ortho_group ', ' problem with SCALAPACK, wrong task columns ID ', 1 )

            ortho_cntx = ortho_cntx_pe(i,k,j)

         END IF

       END DO

     END DO

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
  !
END MODULE mp_global
