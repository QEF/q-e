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
    LOGICAL :: have_task_groups
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

      PUBLIC :: task_groups_init

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

!-----------------------------------------
! Task groups Contributed by C. Bekas, October 2005
! Revised by C. Cavazzoni
!--------------------------------------------

SUBROUTINE task_groups_init( dffts )

   !
   USE fft_types,      ONLY : fft_dlay_descriptor

   ! T.G.
   ! NPGRP:      Number of processors per group
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE
#if defined(__MPI)
   INCLUDE 'mpif.h'
#endif


   TYPE(fft_dlay_descriptor), INTENT(inout) :: dffts

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

   CALL task_groups_init_first( dffts )
   !
#if defined(DEBUG)
   IF ( dffts%nogrp > 1 ) WRITE( stdout, 100 ) dffts%nogrp, dffts%npgrp

100 FORMAT( /,3X,'Task Groups are in USE',/,3X,'groups and procs/group : ',I5,I5 )
#endif

   !Find maximum chunk of local data concerning coefficients of eigenfunctions in g-space

#if defined(__MPI)
   CALL MPI_Allgather( dffts%nnr, 1, MPI_INTEGER, nnrsx_vec, 1, MPI_INTEGER, dffts%comm, IERR)
   strd = maxval( nnrsx_vec( 1:dffts%nproc ) )
#else
   strd = dffts%nnr
#endif

   IF( strd /= dffts%tg_nnr ) CALL fftx_error__( ' task_groups_init ', ' inconsistent nnr ', 1 )

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
   ALLOCATE( dffts%tg_nsw(dffts%nproc))
   ALLOCATE( dffts%tg_npp(dffts%nproc))

   num_sticks = 0
   num_planes = 0
   DO i = 1, dffts%nogrp
      num_sticks = num_sticks + dffts%nsw( dffts%nolist(i) + 1 )
      num_planes = num_planes + dffts%npp( dffts%nolist(i) + 1 )
   ENDDO

#if defined(__MPI)
   CALL MPI_ALLGATHER(num_sticks, 1, MPI_INTEGER, dffts%tg_nsw(1), 1, MPI_INTEGER, dffts%comm, IERR)
   CALL MPI_ALLGATHER(num_planes, 1, MPI_INTEGER, dffts%tg_npp(1), 1, MPI_INTEGER, dffts%comm, IERR)
#else
   dffts%tg_nsw(1) = num_sticks
   dffts%tg_npp(1) = num_planes
#endif

   ALLOCATE( dffts%tg_snd( dffts%nogrp ) )
   ALLOCATE( dffts%tg_rcv( dffts%nogrp ) )
   ALLOCATE( dffts%tg_psdsp( dffts%nogrp ) )
   ALLOCATE( dffts%tg_usdsp( dffts%nogrp ) )
   ALLOCATE( dffts%tg_rdsp( dffts%nogrp ) )

   dffts%tg_snd(1)   = dffts%nr3x * dffts%nsw( dffts%mype + 1 )
   IF( dffts%nr3x * dffts%nsw( dffts%mype + 1 ) > dffts%tg_nnr ) THEN
      CALL fftx_error__( ' task_groups_init ', ' inconsistent dffts%tg_nnr ', 1 )
   ENDIF
   dffts%tg_psdsp(1) = 0
   dffts%tg_usdsp(1) = 0
   dffts%tg_rcv(1)  = dffts%nr3x * dffts%nsw( dffts%nolist(1) + 1 )
   dffts%tg_rdsp(1) = 0
   DO i = 2, dffts%nogrp
      dffts%tg_snd(i)  = dffts%nr3x * dffts%nsw( dffts%mype + 1 )
      dffts%tg_psdsp(i) = dffts%tg_psdsp(i-1) + dffts%tg_nnr
      dffts%tg_usdsp(i) = dffts%tg_usdsp(i-1) + dffts%tg_snd(i-1)
      dffts%tg_rcv(i)  = dffts%nr3x * dffts%nsw( dffts%nolist(i) + 1 )
      dffts%tg_rdsp(i) = dffts%tg_rdsp(i-1) + dffts%tg_rcv(i-1)
   ENDDO

   dffts%tg_ncpx = 0
   dffts%tg_nppx = 0
   DO i = 1, dffts%npgrp
      dffts%tg_ncpx = max( dffts%tg_ncpx, dffts%tg_nsw ( dffts%nplist(i) + 1 ) )
      dffts%tg_nppx = max( dffts%tg_nppx, dffts%tg_npp ( dffts%nplist(i) + 1 ) )
   ENDDO


   RETURN

END SUBROUTINE task_groups_init


  !
SUBROUTINE task_groups_init_first( dffts )
   !
   USE fft_types,      ONLY : fft_dlay_descriptor
   !
   IMPLICIT NONE
#if defined(__MPI)
   INCLUDE 'mpif.h'
#endif
   !
   TYPE(fft_dlay_descriptor), INTENT(inout) :: dffts
    !
    INTEGER :: i, n1, ipos, color, key, ierr, itsk, ntsk
    INTEGER :: pgroup( dffts%nproc )
    !
    !SUBDIVIDE THE PROCESSORS IN GROUPS
    !
    DO i = 1, dffts%nproc
       pgroup( i ) = i - 1
    ENDDO
    !
#if defined(__TASK_GROUP_WAVE_ORDER)
    n1 = ( dffts%mype / dffts%npgrp ) * dffts%npgrp 
    ipos = dffts%mype - n1
#else
    n1 = ( dffts%mype / dffts%nogrp ) * dffts%nogrp 
    ipos = dffts%mype - n1
#endif
    !
    !LIST OF PROCESSORS IN MY ORBITAL GROUP 
    !     (processors dealing with my same pw's of different orbitals)
    !
    !  processors in these group have contiguous indexes
    !
    DO i = 1, dffts%nogrp
#if defined(__TASK_GROUP_WAVE_ORDER)
       dffts%nolist( i ) = pgroup( ipos + ( i - 1 ) * dffts%npgrp + 1 )
#else
       dffts%nolist( i ) = pgroup( n1 + i )
#endif
    ENDDO
    !
    !LIST OF PROCESSORS IN MY PLANE WAVE GROUP
    !     (processors dealing with different pw's of my same orbital)
    !
    DO i = 1, dffts%npgrp
#if defined(__TASK_GROUP_WAVE_ORDER)
       dffts%nplist( i ) = pgroup( n1 + i )
#else
       dffts%nplist( i ) = pgroup( ipos + ( i - 1 ) * dffts%nogrp + 1 )
#endif
    ENDDO
    !
    !SET UP THE GROUPS
    !
    !CREATE ORBITAL GROUPS
    !
#if defined(__MPI)
    ! processes with the same color are in the same new communicator

#if defined(__TASK_GROUP_WAVE_ORDER)
    color = MOD( dffts%mype , dffts%npgrp )
    key   = dffts%mype / dffts%npgrp
#else
    color = dffts%mype / dffts%nogrp
    key   = MOD( dffts%mype , dffts%nogrp )
#endif


    CALL MPI_COMM_SPLIT( dffts%comm, color, key, dffts%ogrp_comm, ierr )
    if( ierr /= 0 ) &
         CALL fftx_error__( ' task_groups_init_first ', ' creating ogrp_comm ', ABS(ierr) )
    CALL MPI_COMM_RANK( dffts%ogrp_comm, itsk, IERR )
    CALL MPI_COMM_SIZE( dffts%ogrp_comm, ntsk, IERR )
    IF( dffts%nogrp /= ntsk ) CALL fftx_error__( ' task_groups_init_first ', ' ogrp_comm size ', ntsk )
    DO i = 1, dffts%nogrp
       IF( dffts%mype == dffts%nolist( i ) ) THEN
          IF( (i-1) /= itsk ) CALL fftx_error__( ' task_groups_init_first ', ' ogrp_comm rank ', itsk )
       END IF
    END DO
#endif
    !
    !CREATE PLANEWAVE GROUPS
    !
#if defined(__MPI)
    ! processes with the same color are in the same new communicator

#if defined(__TASK_GROUP_WAVE_ORDER)
    color = dffts%mype / dffts%npgrp
    key   = MOD( dffts%mype , dffts%npgrp )
#else
    color = MOD( dffts%mype , dffts%nogrp )
    key   = dffts%mype / dffts%nogrp
#endif

    CALL MPI_COMM_SPLIT( dffts%comm, color, key, dffts%pgrp_comm, ierr )
    if( ierr /= 0 ) &
         CALL fftx_error__( ' task_groups_init_first ', ' creating pgrp_comm ', ABS(ierr) )
    CALL MPI_COMM_RANK( dffts%pgrp_comm, itsk, IERR )
    CALL MPI_COMM_SIZE( dffts%pgrp_comm, ntsk, IERR )
    IF( dffts%npgrp /= ntsk ) CALL fftx_error__( ' task_groups_init_first ', ' pgrp_comm size ', ntsk )
    DO i = 1, dffts%npgrp
       IF( dffts%mype == dffts%nplist( i ) ) THEN
          IF( (i-1) /= itsk ) CALL fftx_error__( ' task_groups_init_first ', ' pgrp_comm rank ', itsk )
       END IF
    END DO
    dffts%me_pgrp = itsk
#endif

    RETURN
  END SUBROUTINE task_groups_init_first
  !

!=----------------------------------------------------------------------=
   END MODULE task_groups
!=----------------------------------------------------------------------=
