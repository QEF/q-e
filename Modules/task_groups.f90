!
! Copyright (C) 2002-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------
! Contributed by C. Bekas, October 2005
! Revised by C. Cavazzoni
!--------------------------------------------

MODULE task_groups

   USE kinds,      ONLY: DP

   IMPLICIT NONE
   SAVE

   INTEGER,  ALLOCATABLE :: nolist(:), nplist(:), pgroup(:)
   INTEGER,  ALLOCATABLE :: tmp_nsw(:), tmp_npp(:)
   INTEGER   :: strd
   INTEGER   :: nswx  ! maximum number of stick per processor

CONTAINS

SUBROUTINE DEALLOCATE_GROUPS

   IMPLICIT NONE

   !  ... Deallocate groups related arrays

   IF (ALLOCATED(nolist))       DEALLOCATE(nolist)
   IF (ALLOCATED(nplist))       DEALLOCATE(nplist)
   IF (ALLOCATED(pgroup))       DEALLOCATE(pgroup)
   IF (ALLOCATED(tmp_nsw))      DEALLOCATE(tmp_nsw)
   IF (ALLOCATED(tmp_npp))      DEALLOCATE(tmp_npp)

END SUBROUTINE DEALLOCATE_GROUPS


!========================================================================================
! ADDED SUBROUTINEs FOR TASK GROUP PARALLIZATION
! C. Bekas, IBM Research, Zurich
!        - GROUPS: Define and initialize Task Groups
!        - tg_ivfftw: Inverse FFT driver for Task Groups
!=======================================================================================


!-----------------------------------------------------------------------
!      SUBROUTINE GROUPS (added by C. Bekas)
!      Define groups for task group parallilization
!-----------------------------------------------------------------------

SUBROUTINE task_groups_init( dffts )

   USE parallel_include
   !
   USE mp_global,      ONLY : me_image, nproc_image, intra_image_comm, root
   USE mp_global,      ONLY : NOGRP, NPGRP, ogrp_comm, pgrp_comm  
   USE mp,             ONLY : mp_bcast
   USE io_global,      only : stdout
   USE fft_types,      only : fft_dlay_descriptor

   ! T.G. 
   ! NPGRP:      Number of processors per group
   ! NOGRP:      Number of group

   IMPLICIT NONE 

   TYPE(fft_dlay_descriptor), INTENT(IN) :: dffts

   !----------------------------------
   !Local Variables declaration
   !----------------------------------

   INTEGER  ::  MSGLEN, I, J, N1, IPOS, WORLD, NEWGROUP
   INTEGER  ::  IERR
   INTEGER  ::  itsk, ntsk, color, key
   INTEGER  ::  num_planes, num_sticks
   INTEGER,  DIMENSION(:),   ALLOCATABLE  ::  nnrsx_vec 

   !
   WRITE( stdout, 100 ) nogrp, npgrp

100 FORMAT( /,3X,'Task Groups are in use',/,3X,'groups and procs/group : ',I5,I5 )

   !--------------------------------------------------------------
   !SUBDIVIDE THE PROCESSORS IN GROUPS
   !
   !THE NUMBER OF GROUPS HAS TO BE A DIVISOR OF THE NUMBER
   !OF PROCESSORS
   !--------------------------------------------------------------

   IF( MOD( nproc_image, nogrp ) /= 0 ) &
      CALL errore( " groups ", " nogrp should be a divisor of nproc_image ", 1 )
 
   ALLOCATE( pgroup( nproc_image ) )
   !
   DO i = 1, nproc_image
      pgroup( i ) = i - 1
   ENDDO
   !
   ALLOCATE( nplist( npgrp ) )
   !
   ALLOCATE( nolist( nogrp ) )
   !
   !--------------------------------------
   !LIST OF PROCESSORS IN MY ORBITAL GROUP
   !--------------------------------------
   !
   !  processors in these group have contiguous indexes
   !
   N1 = ( me_image / NOGRP ) * NOGRP - 1
   DO i = 1, nogrp
      nolist( I ) = pgroup( N1 + I + 1 )
      IF( me_image == nolist( I ) ) ipos = i - 1
   ENDDO

   !-----------------------------------------
   !LIST OF PROCESSORS IN MY PLANE WAVE GROUP
   !-----------------------------------------
   !
   DO I = 1, npgrp
      nplist( I ) = pgroup( ipos + ( i - 1 ) * nogrp + 1 )
   ENDDO

   !-----------------
   !SET UP THE GROUPS
   !-----------------
   !

   !---------------------------------------
   !CREATE ORBITAL GROUPS
   !---------------------------------------
   !
#if defined __MPI
   color = me_image / nogrp
   key   = MOD( me_image , nogrp )
   CALL MPI_COMM_SPLIT( intra_image_comm, color, key, ogrp_comm, ierr )
   if( ierr /= 0 ) &
      CALL errore( ' task_groups_init ', ' creating ogrp_comm ', ABS(ierr) )
   CALL MPI_COMM_RANK( ogrp_comm, itsk, IERR )
   CALL MPI_COMM_SIZE( ogrp_comm, ntsk, IERR )
   IF( nogrp /= ntsk ) CALL errore( ' task_groups_init ', ' ogrp_comm size ', ntsk )
   DO i = 1, nogrp
      IF( me_image == nolist( i ) ) THEN
         IF( (i-1) /= itsk ) CALL errore( ' task_groups_init ', ' ogrp_comm rank ', itsk )
      END IF
   END DO
#endif

   !---------------------------------------
   !CREATE PLANEWAVE GROUPS
   !---------------------------------------
   !
#if defined __MPI
   color = MOD( me_image , nogrp )
   key   = me_image / nogrp
   CALL MPI_COMM_SPLIT( intra_image_comm, color, key, pgrp_comm, ierr )
   if( ierr /= 0 ) &
      CALL errore( ' task_groups_init ', ' creating pgrp_comm ', ABS(ierr) )
   CALL MPI_COMM_RANK( pgrp_comm, itsk, IERR )
   CALL MPI_COMM_SIZE( pgrp_comm, ntsk, IERR )
   IF( npgrp /= ntsk ) CALL errore( ' task_groups_init ', ' pgrp_comm size ', ntsk )
   DO i = 1, npgrp
      IF( me_image == nplist( i ) ) THEN
         IF( (i-1) /= itsk ) CALL errore( ' task_groups_init ', ' pgrp_comm rank ', itsk )
      END IF
   END DO
#endif


   ALLOCATE( nnrsx_vec( nproc_image ) )

   !Find maximum chunk of local data concerning coefficients of eigenfunctions in g-space

#if defined __MPI
   CALL MPI_Allgather( dffts%nnr, 1, MPI_INTEGER, nnrsx_vec, 1, MPI_INTEGER, intra_image_comm, IERR)
   strd = MAXVAL( nnrsx_vec( 1:nproc_image ) )
   nswx = MAXVAL( dffts%nsw( 1:nproc_image ) )
#else
   strd = dffts%nnr 
   nswx = dffts%nsw(1)
#endif

   DEALLOCATE( nnrsx_vec )

   !-------------------------------------------------------------------------------------
   !C. Bekas...TASK GROUP RELATED. FFT DATA STRUCTURES ARE ALREADY DEFINED ABOVE
   !-------------------------------------------------------------------------------------
   !dfft%nsw(me) holds the number of z-sticks for the current processor per wave-function
   !We can either send these in the group with an mpi_allgather...or put the
   !in the PSIS vector (in special positions) and send them with them.
   !Otherwise we can do this once at the beginning, before the loop.
   !we choose to do the latter one.
   !-------------------------------------------------------------------------------------
   !
   ALLOCATE(tmp_nsw(nproc_image))
   ALLOCATE(tmp_npp(nproc_image))

   num_sticks = 0
   num_planes = 0
   DO i = 1, nogrp
      num_sticks = num_sticks + dffts%nsw( nolist(i) + 1 )
      num_planes = num_planes + dffts%npp( nolist(i) + 1 )
   ENDDO

#if defined __MPI
   CALL MPI_ALLGATHER(num_sticks, 1, MPI_INTEGER, tmp_nsw, 1, MPI_INTEGER, intra_image_comm, IERR)
   CALL MPI_ALLGATHER(num_planes, 1, MPI_INTEGER, tmp_npp, 1, MPI_INTEGER, intra_image_comm, IERR)
#else
   tmp_nsw(1) = num_sticks
   tmp_npp(1) = num_planes
#endif

   RETURN

   END SUBROUTINE task_groups_init


END MODULE task_groups
