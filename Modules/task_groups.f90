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
   USE parameters, ONLY: maxcpu

   IMPLICIT NONE
   SAVE

   INTEGER, DIMENSION(:), ALLOCATABLE :: ALL_Z_STICKS
   INTEGER, DIMENSION(:), ALLOCATABLE  :: NOLIST, nplist, PGROUP
   INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_nsw, tmp_npp, tmp_planes, tmp_revs, recvs, tmp_ismap
   INTEGER, DIMENSION(:), ALLOCATABLE :: ngw_vec !GLOBAL VECTOR OF ALL NGW VECTORS
   COMPLEX(DP), DIMENSION(:,:),  ALLOCATABLE :: tg_betae
   COMPLEX(DP), DIMENSION(:),    ALLOCATABLE :: tg_c2, tg_c3
   REAL(DP),  DIMENSION(:,:), ALLOCATABLE :: tg_rhos
   REAL(DP),  DIMENSION(:),   ALLOCATABLE :: tg_ggp
   INTEGER,  DIMENSION(:),   ALLOCATABLE  ::  nnrsx_vec !INCREASE THIS TO THE MAXIMUM NUMBER OF PROCS IF NEEDED
   REAL(DP), DIMENSION(:,:), ALLOCATABLE  ::  tmp_rhos, local_rhos
   INTEGER  :: SZ, CLOCK1, CLOCK2, CLOCK3, CLOCK4
   INTEGER  :: sticks_index, eig_offset, strd       
   REAL(DP) :: tm_tg, tm_rhoofr

CONTAINS

SUBROUTINE DEALLOCATE_GROUPS

   IMPLICIT NONE

   !DEALLOCATE GROUPS RELATED ARRAYS

   IF (ALLOCATED(ALL_Z_STICKS)) DEALLOCATE(ALL_Z_STICKS)
   IF (ALLOCATED(NOLIST))       DEALLOCATE(NOLIST)
   IF (ALLOCATED(nplist))       DEALLOCATE(nplist)
   IF (ALLOCATED(PGROUP))       DEALLOCATE(PGROUP)
   IF (ALLOCATED(tmp_nsw))      DEALLOCATE(tmp_nsw)
   IF (ALLOCATED(tmp_npp))      DEALLOCATE(tmp_npp)
   IF (ALLOCATED(tmp_planes))   DEALLOCATE(tmp_planes)
   IF (ALLOCATED(tmp_revs))     DEALLOCATE(tmp_revs)
   IF (ALLOCATED(recvs))        DEALLOCATE(recvs)
   IF (ALLOCATED(tmp_ismap))    DEALLOCATE(tmp_ismap)
   IF (ALLOCATED(ngw_vec))      DEALLOCATE(ngw_vec)
   IF (ALLOCATED(tg_betae))     DEALLOCATE(tg_betae)
   IF (ALLOCATED(tg_c2))        DEALLOCATE(tg_c2)
   IF (ALLOCATED(tg_c3))        DEALLOCATE(tg_c3)
   IF (ALLOCATED(tg_rhos))      DEALLOCATE(tg_rhos)
   IF (ALLOCATED(tg_ggp))       DEALLOCATE(tg_ggp)
   IF (ALLOCATED(nnrsx_vec))    DEALLOCATE(nnrsx_vec)
   IF (ALLOCATED(tmp_rhos))     DEALLOCATE(tmp_rhos)
   IF (ALLOCATED(local_rhos))   DEALLOCATE(local_rhos)

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

   USE mp_global,  ONLY : me_image, nproc_image, intra_image_comm, root
   USE mp_global,  ONLY : NOGRP, NPGRP, ME_OGRP, ME_PGRP  
   USE mp,         ONLY : mp_bcast
   USE io_global,  only : stdout
   USE fft_types,  only : fft_dlay_descriptor
   USE electrons_base, only: nspin
   USE parallel_include

   ! T.G. 
   ! NPGRP:      Number of processors per group
   ! NOGRP:      Number of group

   IMPLICIT NONE 

   TYPE(fft_dlay_descriptor), INTENT(IN) :: dffts

   !----------------------------------
   !Local Variables declaration
   !----------------------------------

   INTEGER  ::  MSGLEN, I, J, N1, LABEL, IPOS, WORLD, NEWGROUP
   INTEGER  ::  IERR

   !--------------------------------------------------------------
   !

   WRITE( stdout, 100 ) nogrp, npgrp

100 FORMAT( /,3X,'Task Groups are in use',/,3X,'groups and procs/group : ',I5,I5 )

   tm_tg     = 0D0
   tm_rhoofr = 0D0

   !  Find the number of processors and my rank
   !
   sz = nproc_image

   !--------------------------------------------------------------
   !SUBDIVIDE THE PROCESSORS IN GROUPS
   !
   !THE NUMBER OF GROUPS HAS TO BE A DIVISOR OF THE NUMBER
   !OF PROCESSORS
   !--------------------------------------------------------------

   IF( MOD( nproc_image, nogrp ) /= 0 ) &
      CALL errore( " groups ", " nogrp should be a divisor of nproc_image ", 1 )

 
   ALLOCATE( pgroup( nproc_image ) )
   DO i = 1, nproc_image
      pgroup( i ) = i - 1
   ENDDO

   ALLOCATE( nnrsx_vec( sz ) )

   !Find maximum chunk of local data concerning coefficients of eigenfunctions in g-space

#if defined __MPI
   CALL MPI_Allgather(dffts%nnr, 1, MPI_INTEGER, nnrsx_vec, 1, MPI_INTEGER, intra_image_comm, IERR)
   strd = MAXVAL( nnrsx_vec( 1:sz ) )
#else
   strd = dffts%nnr 
#endif


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
   ALLOCATE( all_z_sticks( sz ) )
   !
   !ALL-Gather number of Z-sticks from all processors
   !
#if defined __MPI
   CALL MPI_Allgather(dffts%nsw(me_image+1), 1, MPI_INTEGER, all_z_sticks, 1, MPI_INTEGER, intra_image_comm, ierr)
#else
   all_z_sticks( 1 ) = dffts%nsw( 1 )
#endif
   IF (.NOT.ALLOCATED(tmp_nsw)) ALLOCATE(tmp_nsw(sz))
   IF (.NOT.ALLOCATED(tmp_npp)) THEN
      ALLOCATE(tmp_npp(sz))
      tmp_npp(1)=-1
   ENDIF


   ALLOCATE( nplist( npgrp ) )
   !
   ALLOCATE( nolist( nogrp ) )


   !--------------------------------------
   !LIST OF PROCESSORS IN MY ORBITAL GROUP
   !--------------------------------------
   N1 = ( me_image / NOGRP ) * NOGRP - 1
   DO I = 1, NOGRP
      NOLIST( I ) = PGROUP( N1 + I + 1 )
      IF( me_image .EQ. NOLIST( I ) ) IPOS = I - 1
   ENDDO

   !-----------------------------------------
   !LIST OF PROCESSORS IN MY PLANE WAVE GROUP
   !-----------------------------------------
   DO I = 1, NPGRP
      nplist( I ) = PGROUP( IPOS + ( I - 1 ) * NOGRP + 1 )
   ENDDO

   !-----------------
   !SET UP THE GROUPS
   !-----------------
   DO I = 1, NPGRP
      IF( me_image == nplist( I ) ) LABEL = I
   ENDDO

   !---------------------------------------
   !CREATE ORBITAL GROUPS
   !---------------------------------------
   !
#if defined __MPI
   ! get the handle world to the main group
   CALL MPI_COMM_GROUP( intra_image_comm, WORLD, IERR )  
   ! create a new group handle containing processor whose indices are
   ! contained in array nogrp()
   CALL MPI_GROUP_INCL( WORLD, NOGRP, NOLIST, NEWGROUP, IERR )
   ! build a communicator for the newgroup and store it in ME_OGRP
   CALL MPI_COMM_CREATE( intra_image_comm, NEWGROUP, ME_OGRP, IERR )
#endif

   DO I = 1, NOGRP
      IF( me_image .EQ. NOLIST( I ) ) LABEL = I + MAXCPU
   ENDDO

   !---------------------------------------
   !CREATE PLANEWAVE GROUPS
   !---------------------------------------
   !
#if defined __MPI
   CALL MPI_COMM_GROUP( intra_image_comm, WORLD, IERR )
   CALL MPI_GROUP_INCL( WORLD, NPGRP, nplist, NEWGROUP, IERR )
   CALL MPI_COMM_CREATE( intra_image_comm, NEWGROUP, ME_PGRP, IERR )
#endif

   RETURN

   END SUBROUTINE task_groups_init


END MODULE task_groups
