!
! Copyright (C) 2002-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_global
  !----------------------------------------------------------------------------
  !
#if defined (__SHMEM)
  USE shmem_include
#endif
  !
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
  INTEGER :: kunit = 1  ! granularity of k-point distribution
  !
  ! ... indeces ( all starting from 0 !!! )
  !
  INTEGER :: me_pool     = 0  ! index of the processor within a pool 
  INTEGER :: me_image    = 0  ! index of the processor within an image
  INTEGER :: root_pool   = 0  ! index of the root processor within a pool
  INTEGER :: root_image  = 0  ! index of the root processor within an image
  INTEGER :: my_pool_id  = 0  ! index of my pool
  INTEGER :: my_image_id = 0  ! index of my image
  INTEGER :: me_ortho(2) = 0  ! coordinates of the processors
  !
  INTEGER :: npool       = 1  ! number of "k-points"-pools
  INTEGER :: nimage      = 1  ! number of "path-images"-pools
  INTEGER :: nogrp       = 1  ! number of proc. in an orbital "task group" 
  INTEGER :: npgrp       = 1  ! number of proc. in a plane-wave "task group" 
  INTEGER :: nproc_pool  = 1  ! number of processor within a pool
  INTEGER :: nproc_pool_file  = 1  ! number of processor within a pool of
                              !   written in the xml punch file
  INTEGER :: nproc_image = 1  ! number of processor within an image
  INTEGER :: np_ortho(2) = 1  ! size of the processor grid used in ortho
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
  INTEGER :: pgrp_comm        = 0  ! plane-wave group communicator (task grouping)
  INTEGER :: ogrp_comm        = 0  ! orbital group communicarot (task grouping)
  INTEGER :: ortho_comm       = 0  ! communicator used for fast and memory saving ortho
  INTEGER :: ortho_comm_id    = 0  ! id of the ortho_comm
  !
  CONTAINS
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
       my_pool_id       = 0
       my_image_id      = 0
       me_pool          = mpime
       me_image         = mpime
       root_pool        = root
       root_image       = root
       inter_pool_comm  = group_i
       intra_pool_comm  = group_i
       inter_image_comm = group_i
       intra_image_comm = group_i
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
SUBROUTINE init_pool( nimage_ , ntask_groups_ , nproc_ortho_ )
  !----------------------------------------------------------------------------
  !
  ! ... This routine initialize the pool :  MPI division in pools and images
  !
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER, OPTIONAL, INTENT(IN) :: nimage_
  INTEGER, OPTIONAL, INTENT(IN) :: ntask_groups_
  INTEGER, OPTIONAL, INTENT(IN) :: nproc_ortho_
  !
  INTEGER :: ierr = 0
  INTEGER :: nproc_ortho
  !
#if defined (__PARA)
  ! 
  !
  IF( PRESENT( nimage_ ) ) THEN
     nimage = nimage_
  END IF
  !  
  ! ... here we set all parallel indeces (defined in mp_global): 
  !
  !
  ! ... number of cpus per image
  !
  nproc_image = nproc / nimage
  !
  IF ( nproc < nimage ) &
     CALL errore( 'startup', 'nproc < nimage', 1 )
  !
  IF ( MOD( nproc, nimage ) /= 0 ) &
     CALL errore( 'startup', 'nproc /= nproc_image * nimage', 1 ) 
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
  CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, &
                       my_image_id, mpime, intra_image_comm, ierr )
  !
  CALL errore( 'init_pool', 'intra_image_comm is wrong', ierr )
  !
  CALL mp_barrier()
  !
  ! ... the inter_image_comm communicator is created                     
  !     
  CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, &
                       me_image, mpime, inter_image_comm, ierr )  
  !
  CALL errore( 'init_pool', 'inter_image_comm is wrong', ierr )
  !
  ! ... number of cpus per pool of k-points (they are created inside each image)
  !
  nproc_pool = nproc_image / npool
  !
  IF ( MOD( nproc, npool ) /= 0 ) &
     CALL errore( 'startup', 'nproc /= nproc_pool * npool', 1 )  
  !
  ! ... my_pool_id  =  pool index for this processor    ( 0 : npool - 1 )
  ! ... me_pool     =  processor index within the pool  ( 0 : nproc_pool - 1 )
  !
  my_pool_id = me_image / nproc_pool    
  me_pool    = MOD( me_image, nproc_pool )
  !
  CALL mp_barrier( intra_image_comm )
  !
  ! ... the intra_pool_comm communicator is created
  !
  CALL MPI_COMM_SPLIT( intra_image_comm, &
                       my_pool_id, me_image, intra_pool_comm, ierr )
  !
  CALL errore( 'init_pool', 'intra_pool_comm is wrong', ierr )
  !
  CALL mp_barrier( intra_image_comm )
  !
  ! ... the inter_pool_comm communicator is created
  !
  CALL MPI_COMM_SPLIT( intra_image_comm, &
                       me_pool, me_image, inter_pool_comm, ierr )
  !
  call errore( 'init_pool', 'inter_pool_comm is wrong', ierr )
  !
#endif
  !
  nproc_ortho = nproc_pool
  !
  IF( PRESENT( nproc_ortho_ ) ) THEN
     IF( nproc_ortho_ < nproc_pool ) nproc_ortho = nproc_ortho_
  END IF
  !
  CALL init_ortho_group( nproc_ortho, intra_pool_comm )
  !  
  IF( PRESENT( ntask_groups_ ) ) THEN
     IF( ntask_groups_ > 0 ) THEN
        nogrp = ntask_groups_
        CALL init_task_groups( )
     END IF
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
      CALL errore( " init_pool ", " nogrp should be a divisor of nproc_pool ", 1 )
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
     CALL errore( ' task_groups_init ', ' creating ogrp_comm ', ABS(ierr) )
  CALL MPI_COMM_RANK( ogrp_comm, itsk, IERR )
  CALL MPI_COMM_SIZE( ogrp_comm, ntsk, IERR )
  IF( nogrp /= ntsk ) CALL errore( ' task_groups_init ', ' ogrp_comm size ', ntsk )
  DO i = 1, nogrp
     IF( me_pool == nolist( i ) ) THEN
        IF( (i-1) /= itsk ) CALL errore( ' task_groups_init ', ' ogrp_comm rank ', itsk )
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
     CALL errore( ' task_groups_init ', ' creating pgrp_comm ', ABS(ierr) )
  CALL MPI_COMM_RANK( pgrp_comm, itsk, IERR )
  CALL MPI_COMM_SIZE( pgrp_comm, ntsk, IERR )
  IF( npgrp /= ntsk ) CALL errore( ' task_groups_init ', ' pgrp_comm size ', ntsk )
  DO i = 1, npgrp
     IF( me_pool == nplist( i ) ) THEN
        IF( (i-1) /= itsk ) CALL errore( ' task_groups_init ', ' pgrp_comm rank ', itsk )
     END IF
  END DO
#endif

  
  RETURN
END SUBROUTINE init_task_groups
!
!
SUBROUTINE init_ortho_group( nproc_try, comm_all )
   !
   USE mp, ONLY : mp_comm_free, mp_size, mp_rank
   !
   IMPLICIT NONE
    
   INTEGER, INTENT(IN) :: nproc_try, comm_all
    
   LOGICAL, SAVE :: first = .true.
   INTEGER :: ierr, color, key, me_all, newid, nproc_all
   INTEGER :: i, nproc_ortho
    
#if defined __MPI

   me_all    = mp_rank( comm_all )
   nproc_all = mp_size( comm_all )

   IF( nproc_try > nproc_all ) THEN
      CALL errore( " init_ortho_group ", " argument 1 out of range ", nproc_try )
   END IF

   IF( .NOT. first ) THEN
      !  
      !  free resources associated to the communicator
      !
      CALL mp_comm_free( ortho_comm )
      !
   END IF

   !  find the square closer (but lower) to nproc_try
   !
   CALL grid2d_dims( 'S', nproc_try, np_ortho(1), np_ortho(2) )
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
   !  initialize the communicator for the new group
   !
   CALL MPI_COMM_SPLIT( comm_all, color, key, ortho_comm, ierr )
   IF( ierr /= 0 ) &
      CALL errore( " init_ortho_group ", " error splitting communicator ", ierr )
   !
   !  Computes coordinates of the processors, in row maior order
   !
   newid       = mp_rank( ortho_comm )
   nproc_ortho = mp_size( ortho_comm )
   IF( color == 1 .AND. nproc_ortho /= np_ortho(1) * np_ortho(2) ) &
      CALL errore( " init_ortho_group ", " wrong number of proc in ortho_comm ", ierr )
   !
   IF( me_all == 0 .AND. newid /= 0 ) &
      CALL errore( " init_ortho_group ", " wrong root in ortho_comm ", ierr )
   !
   if( color == 1 ) then
      ortho_comm_id = 1
      CALL GRID2D_COORDS( 'R', newid, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
      CALL GRID2D_RANK( 'R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr )
      IF( ierr /= newid ) &
         CALL errore( " init_ortho_group ", " wrong coordinates in ortho_comm ", ierr )
      IF( newid*leg_ortho /= me_all ) &
         CALL errore( " init_ortho_group ", " wrong rank assignment in ortho_comm ", ierr )
   else
      ortho_comm_id = 0
      me_ortho(1) = newid
      me_ortho(2) = newid
   endif

#else

   ortho_comm_id = 1

#endif

   first = .false.
    
   RETURN
END SUBROUTINE init_ortho_group
     !
     !
END MODULE mp_global
