!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------=
   MODULE task_groups
!=----------------------------------------------------------------------=

!  ... Distribute G-vectors across processors as sticks and planes,
!  ... initialize FFT descriptors for both dense and smooth grids

      IMPLICIT NONE

      INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

      TYPE task_groups_descriptor
        !
        !  parent group parallelization
        !
        INTEGER :: mype     = 0          ! my processor id (starting from 0) in the parent fft group
        INTEGER :: comm     = 0          ! communicator of the parent fft gruop 
        INTEGER :: nproc    = 1          ! number of processor in the parent fft group
        INTEGER :: root     = 0          ! root processor
        !
        !  task groups
        !
        LOGICAL :: have_task_groups = .FALSE.
        !
        INTEGER :: me_pgrp   = 0          ! task id for plane wave task group
        INTEGER :: nogrp     = 1          ! number of proc. in an orbital "task group"
        INTEGER :: npgrp     = 1          ! number of proc. in a plane-wave "task group"
        INTEGER :: ogrp_comm = 0          ! orbital group communicator
        INTEGER :: pgrp_comm = 0          ! plane-wave group communicator
        INTEGER, ALLOCATABLE :: nolist(:) ! list of pes in orbital group
        INTEGER, ALLOCATABLE :: nplist(:) ! list of pes in pw group
        !
        INTEGER :: tg_nnr = 0            ! maximum among nnr
        INTEGER, ALLOCATABLE :: tg_nsw(:) ! number of sticks per task group ( wave func )
        INTEGER, ALLOCATABLE :: tg_npp(:) ! number of "Z" planes per task group
        INTEGER, ALLOCATABLE :: tg_snd(:) ! number of element to be sent in group redist
        INTEGER, ALLOCATABLE :: tg_rcv(:) ! number of element to be received in group redist
        INTEGER, ALLOCATABLE :: tg_psdsp(:)! send displacement for all to all (pack)
        INTEGER, ALLOCATABLE :: tg_usdsp(:)! send displacement for all to all (unpack)
        INTEGER, ALLOCATABLE :: tg_rdsp(:)! receive displacement for all to all
        INTEGER :: tg_nppx = 0  ! max of tg_npp
        INTEGER :: tg_ncpx = 0  ! max of tg_ncpx
        !
      END TYPE

      PRIVATE
      SAVE

      PUBLIC :: task_groups_init, task_groups_descriptor, task_groups_deallocate

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

!-----------------------------------------
! Task groups Contributed by C. Bekas, October 2005
! Revised by C. Cavazzoni
!--------------------------------------------

SUBROUTINE task_groups_init( dffts, dtgs, nogrp )

   !
   USE fft_types,      ONLY : fft_type_descriptor

   ! T.G.
   ! NPGRP:      Number of processors per group
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE
#if defined(__MPI)
   INCLUDE 'mpif.h'
#endif


   TYPE(fft_type_descriptor), INTENT(inout) :: dffts
   TYPE(task_groups_descriptor), INTENT(inout) :: dtgs
   INTEGER, INTENT(in) :: nogrp   ! number of task groups

   INTEGER :: stdout = 6

   !----------------------------------
   !Local Variables declaration
   !----------------------------------

   INTEGER  :: I
   INTEGER  :: IERR
   INTEGER  :: num_planes, num_sticks
   INTEGER  :: nnrsx_vec ( dffts%nproc )
   INTEGER  :: pgroup( dffts%nproc )
   INTEGER  :: strd

   CALL task_groups_deallocate( dtgs )

   dtgs%mype  = dffts%mype
   dtgs%comm  = dffts%comm
   dtgs%nproc = dffts%nproc
   dtgs%root  = dffts%root

   CALL task_groups_init_first( dffts, dtgs, nogrp )
   !
   !Find maximum chunk of local data concerning coefficients of eigenfunctions in g-space

#if defined(__MPI)
   CALL MPI_Allgather( dffts%nnr, 1, MPI_INTEGER, nnrsx_vec, 1, MPI_INTEGER, dffts%comm, IERR)
   strd = maxval( nnrsx_vec( 1:dffts%nproc ) )
#else
   strd = dffts%nnr
#endif

   IF( strd /= dtgs%tg_nnr ) CALL fftx_error__( ' task_groups_init ', ' inconsistent nnr ', 1 )

   !-------------------------------------------------------------------------------------
   !C. Bekas...TASK GROUP RELATED. FFT DATA STRUCTURES ARE ALREADY DEFINED ABOVE
   !-------------------------------------------------------------------------------------
   !dfft%nsw(me) holds the number of z-sticks for the current processor per wave-function
   !We can either send these in the group with an mpi_allgather...or put them
   !in the PSIS vector (in special positions) and send them with them.
   !Otherwise we can do this once at the beginning, before the loop.
   !we choose to do the latter one.
   !-------------------------------------------------------------------------------------
   !
   !
   ALLOCATE( dtgs%tg_nsw(dtgs%nproc))
   ALLOCATE( dtgs%tg_npp(dtgs%nproc))

   num_sticks = 0
   num_planes = 0
   DO i = 1, dtgs%nogrp
      num_sticks = num_sticks + dffts%nsw( dtgs%nolist(i) + 1 )
      num_planes = num_planes + dffts%npp( dtgs%nolist(i) + 1 )
   ENDDO

#if defined(__MPI)
   CALL MPI_ALLGATHER(num_sticks, 1, MPI_INTEGER, dtgs%tg_nsw(1), 1, MPI_INTEGER, dtgs%comm, IERR)
   CALL MPI_ALLGATHER(num_planes, 1, MPI_INTEGER, dtgs%tg_npp(1), 1, MPI_INTEGER, dtgs%comm, IERR)
#else
   dtgs%tg_nsw(1) = num_sticks
   dtgs%tg_npp(1) = num_planes
#endif

#if defined(__TASK_PRINTOUT)
   write (6,*) 'TASK GROUP stick AND plane DISTRIBUTION '
   do i=1,dtgs%nproc
      write (6,*) i, dtgs%tg_nsw(i), dtgs%tg_npp(i)
   end do
#endif

   ALLOCATE( dtgs%tg_snd( dtgs%nogrp ) )
   ALLOCATE( dtgs%tg_rcv( dtgs%nogrp ) )
   ALLOCATE( dtgs%tg_psdsp( dtgs%nogrp ) )
   ALLOCATE( dtgs%tg_usdsp( dtgs%nogrp ) )
   ALLOCATE( dtgs%tg_rdsp( dtgs%nogrp ) )

   dtgs%tg_snd(1)   = dffts%nr3x * dffts%nsw( dffts%mype + 1 )
   IF( dffts%nr3x * dffts%nsw( dffts%mype + 1 ) > dtgs%tg_nnr ) THEN
      CALL fftx_error__( ' task_groups_init ', ' inconsistent dtgs%tg_nnr ', 1 )
   ENDIF
   dtgs%tg_psdsp(1) = 0
   dtgs%tg_usdsp(1) = 0
   dtgs%tg_rcv(1)  = dffts%nr3x * dffts%nsw( dtgs%nolist(1) + 1 )
   dtgs%tg_rdsp(1) = 0
   DO i = 2, dtgs%nogrp
      dtgs%tg_snd(i)  = dffts%nr3x * dffts%nsw( dffts%mype + 1 )
      dtgs%tg_psdsp(i) = dtgs%tg_psdsp(i-1) + dtgs%tg_nnr
      dtgs%tg_usdsp(i) = dtgs%tg_usdsp(i-1) + dtgs%tg_snd(i-1)
      dtgs%tg_rcv(i)  = dffts%nr3x * dffts%nsw( dtgs%nolist(i) + 1 )
      dtgs%tg_rdsp(i) = dtgs%tg_rdsp(i-1) + dtgs%tg_rcv(i-1)
   ENDDO

   dtgs%tg_ncpx = 0
   dtgs%tg_nppx = 0
   DO i = 1, dtgs%npgrp
      dtgs%tg_ncpx = max( dtgs%tg_ncpx, dtgs%tg_nsw ( dtgs%nplist(i) + 1 ) )
      dtgs%tg_nppx = max( dtgs%tg_nppx, dtgs%tg_npp ( dtgs%nplist(i) + 1 ) )
   ENDDO


   RETURN

END SUBROUTINE task_groups_init


  !
SUBROUTINE task_groups_init_first( dffts, dtgs, nogrp )
   !
   USE fft_types,      ONLY : fft_type_descriptor
   !
   IMPLICIT NONE
#if defined(__MPI)
   INCLUDE 'mpif.h'
#endif
   !
   TYPE(fft_type_descriptor), INTENT(inout) :: dffts
   TYPE(task_groups_descriptor), INTENT(inout) :: dtgs

   INTEGER, INTENT(in) :: nogrp   ! number of task groups
    !
    INTEGER :: i, nlen, n1, ipos, color, key, ierr, itsk, ntsk
    INTEGER :: nppx, ncpx
#if defined(__TASK_MAPPING)
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME), ALLOCATABLE :: proc_name(:)
    INTEGER, ALLOCATABLE :: proc_id(:)
#endif
    !
    !SUBDIVIDE THE PROCESSORS IN GROUPS
    !
    dtgs%have_task_groups = ( nogrp > 1 )
    dtgs%me_pgrp = 0

    IF( MOD( dtgs%nproc, MAX( 1, nogrp ) ) /= 0 ) &
       CALL fftx_error__("task_groups_init_first","the number of task groups should be a divisor of the number of MPI task",1)
    IF( nogrp > dtgs%nproc ) &
       CALL fftx_error__( "task_groups_init_first","the number of task groups should be less than the number of MPI task",1)

    dtgs%nogrp = MAX( 1, nogrp )
    dtgs%npgrp = dtgs%nproc / MAX( 1, nogrp )
    dtgs%ogrp_comm = 0
    dtgs%pgrp_comm = 0
    ALLOCATE( dtgs%nolist( dtgs%nogrp ) )
    ALLOCATE( dtgs%nplist( dtgs%npgrp ) )
    dtgs%nolist = 0
    dtgs%nplist = 0


    nppx = 0
    ncpx = 0
    DO i = 1, dffts%nproc
       nppx = MAX( nppx, dffts%npp( i ) )  ! maximum number of planes per processor
       ncpx = MAX( ncpx, dffts%nsp( i ) )  ! maximum number of columns per processor
    END DO

    IF ( dtgs%nproc == 1 ) THEN
      dtgs%tg_nnr = dffts%nnr
    ELSE
      dtgs%tg_nnr = dffts%nnr
      ! this is required to contain the local data in G space (should be already granted!)
      dtgs%tg_nnr = max( dtgs%tg_nnr, dffts%nr3x * ncpx ) 
      ! this is required to contain the local data in R space (should be already granted!)
      dtgs%tg_nnr = max( dtgs%tg_nnr, dffts%nr1x * dffts%nr2x * nppx ) 
      dtgs%tg_nnr = max( 1, dtgs%tg_nnr ) ! ensure that dffts%nrr > 0 ( for extreme parallelism )
    ENDIF

    !
    !SET UP THE GROUPS
    !
    !CREATE ORBITAL GROUPS
    !LIST OF PROCESSORS IN MY ORBITAL GROUP 
    !     (processors dealing with my same pw's of different orbitals)
    !
#if defined(__MPI)
    ! processes with the same color are in the same new communicator

#if defined(__TASK_GROUP_WAVE_ORDER)
    color = MOD( dtgs%mype , dtgs%npgrp )
    key   = dtgs%mype / dtgs%npgrp
#else
    color = dtgs%mype / dtgs%nogrp
    key   = MOD( dtgs%mype , dtgs%nogrp )
#endif

    CALL MPI_COMM_SPLIT( dtgs%comm, color, key, dtgs%ogrp_comm, ierr )
    if( ierr /= 0 ) &
         CALL fftx_error__( ' task_groups_init_first ', ' creating ogrp_comm ', ABS(ierr) )
    CALL MPI_COMM_RANK( dtgs%ogrp_comm, itsk, IERR )
    CALL MPI_COMM_SIZE( dtgs%ogrp_comm, ntsk, IERR )
    IF( dtgs%nogrp /= ntsk ) CALL fftx_error__( ' task_groups_init_first ', ' ogrp_comm size ', ntsk )
    dtgs%nolist = 0
    dtgs%nolist( itsk + 1 ) = dtgs%mype
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, dtgs%nolist, dtgs%nogrp, MPI_INTEGER, MPI_SUM, dtgs%ogrp_comm, ierr)
#endif
    !
    !CREATE PLANEWAVE GROUPS
    !LIST OF PROCESSORS IN MY PLANE WAVE GROUP
    !     (processors dealing with different pw's of my same orbital)
    !
    !
#if defined(__MPI)
    ! processes with the same color are in the same new communicator

    CALL MPI_COMM_SPLIT( dtgs%comm, key, color, dtgs%pgrp_comm, ierr )
    if( ierr /= 0 ) &
         CALL fftx_error__( ' task_groups_init_first ', ' creating pgrp_comm ', ABS(ierr) )
    CALL MPI_COMM_RANK( dtgs%pgrp_comm, itsk, IERR )
    CALL MPI_COMM_SIZE( dtgs%pgrp_comm, ntsk, IERR )
    IF( dtgs%npgrp /= ntsk ) CALL fftx_error__( ' task_groups_init_first ', ' pgrp_comm size ', ntsk )
    dtgs%nplist = 0
    dtgs%nplist( itsk + 1 ) = dtgs%mype
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, dtgs%nplist, dtgs%npgrp, MPI_INTEGER, MPI_SUM, dtgs%pgrp_comm, ierr)
    dtgs%me_pgrp = itsk
#endif

#if defined(__TASK_MAPPING)

    allocate( proc_name( dtgs%nproc ) )
    allocate( proc_id( dtgs%nproc ) )

    do i = 1, dtgs%nproc
#if defined(__MPI)
       if( dtgs%mype == i-1 ) then
          call MPI_Get_processor_name( proc_name(i), nlen, ierr )
       end if
       CALL MPI_BCAST( nlen, 1, MPI_INT, i-1, dtgs%comm, ierr )
       CALL MPI_BCAST( proc_name(i), MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, i-1, dtgs%comm, ierr )
#else
       proc_name(i) = 'localhost'
#endif
    end do

    IF( dtgs%mype == 0 ) THEN
       WRITE(6, 70) 
       DO i = 1, dtgs%nproc
          WRITE(6,100)  i-1, proc_name( i )
       END DO
    END IF
 70 FORMAT( 'MPI Task to node MAP' )
100 FORMAT( I5, '   ', A20 )

    proc_id = 0
    CALL MPI_COMM_RANK( dtgs%ogrp_comm, itsk, IERR )
    proc_id( dtgs%mype + 1 ) = itsk
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, proc_id, dtgs%nproc, MPI_INTEGER, MPI_SUM, dtgs%comm, ierr)
    IF( dtgs%mype == 0 ) THEN
       WRITE(6, 90) 
       DO i = 1, dtgs%nproc
          WRITE(6,100) proc_id( i ), proc_name( i )
       END DO
    END IF
 90 FORMAT( 'Task Group: Orbital Groups' )
    proc_id = 0
    CALL MPI_COMM_RANK( dtgs%pgrp_comm, itsk, IERR )
    proc_id( dtgs%mype + 1 ) = itsk
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, proc_id, dtgs%nproc, MPI_INTEGER, MPI_SUM, dtgs%comm, ierr)
    IF( dtgs%mype == 0 ) THEN
       WRITE(6, 80) 
       DO i = 1, dtgs%nproc
          WRITE(6,100) proc_id( i ), proc_name( i )
       END DO
    END IF
 80 FORMAT( 'Task Group: Wave Groups' )
    
    deallocate( proc_id )
    deallocate( proc_name )

#endif

    RETURN
  END SUBROUTINE task_groups_init_first
  !

  SUBROUTINE task_groups_deallocate ( desc )
    TYPE (task_groups_descriptor), INTENT(inout) :: desc
    IF ( ALLOCATED( desc%nolist ) )   DEALLOCATE( desc%nolist )
    IF ( ALLOCATED( desc%nplist ) )   DEALLOCATE( desc%nplist )
    IF ( ALLOCATED( desc%tg_nsw ) )   DEALLOCATE( desc%tg_nsw )
    IF ( ALLOCATED( desc%tg_npp ) )   DEALLOCATE( desc%tg_npp )
    IF ( ALLOCATED( desc%tg_snd ) )   DEALLOCATE( desc%tg_snd )
    IF ( ALLOCATED( desc%tg_rcv ) )   DEALLOCATE( desc%tg_rcv )
    IF ( ALLOCATED( desc%tg_psdsp ) )   DEALLOCATE( desc%tg_psdsp )
    IF ( ALLOCATED( desc%tg_usdsp ) )   DEALLOCATE( desc%tg_usdsp )
    IF ( ALLOCATED( desc%tg_rdsp ) )   DEALLOCATE( desc%tg_rdsp )
    desc%have_task_groups = .FALSE.
    RETURN
  END SUBROUTINE task_groups_deallocate

!=----------------------------------------------------------------------=
   END MODULE task_groups
!=----------------------------------------------------------------------=
