!-----------------------------------------
! VARIABLES FOR TASKGROUPS
! C. Bekas, October 2005
!-----------------------------------------
!
!Variable description
!--------------------
!MAXGRP:	Maximum number of task-groups
!--------------------------------------------

MODULE GROUPS_MODULE

   USE kinds, ONLY: DP
   USE parameters, ONLY: MAXCPU, MAXGRP

   IMPLICIT NONE
   SAVE

   INTEGER, DIMENSION(:), ALLOCATABLE :: ALL_Z_STICKS
   INTEGER, DIMENSION(:), ALLOCATABLE  :: NOLIST, NPLIST, PGROUP
   INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_nsw, tmp_npp, tmp_planes, tmp_revs, recvs, tmp_ismap
   INTEGER, DIMENSION(:), ALLOCATABLE :: ngw_vec !GLOBAL VECTOR OF ALL NGW VECTORS
   COMPLEX(DP), DIMENSION(:,:),  ALLOCATABLE :: tg_betae
   COMPLEX(DP), DIMENSION(:),    ALLOCATABLE :: tg_c2, tg_c3
   REAL(DP),  DIMENSION(:,:), ALLOCATABLE :: tg_rhos
   REAL(DP),  DIMENSION(:),   ALLOCATABLE :: tg_ggp
   INTEGER,  DIMENSION(:),   ALLOCATABLE  ::  nnrsx_vec !INCREASE THIS TO THE MAXIMUM NUMBER OF PROCS IF NEEDED
   REAL(DP), DIMENSION(:,:), ALLOCATABLE  ::  tmp_rhos, local_rhos
   INTEGER  :: recv_cnt(MAXGRP), recv_displ(MAXGRP), send_cnt(MAXGRP), send_displ(MAXGRP)
   INTEGER  :: SZ, CLOCK1, CLOCK2, CLOCK3, CLOCK4
   INTEGER  :: sticks_index, eig_offset, strd       
   REAL(DP) :: tm_tg, tm_rhoofr

CONTAINS

SUBROUTINE DEALLOCATE_GROUPS

   IMPLICIT NONE

   !DEALLOCATE GROUPS RELATED ARRAYS

   IF (ALLOCATED(ALL_Z_STICKS)) DEALLOCATE(ALL_Z_STICKS)
   IF (ALLOCATED(NOLIST))       DEALLOCATE(NOLIST)
   IF (ALLOCATED(NPLIST))       DEALLOCATE(NPLIST)
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
!	- GROUPS: Define and initialize Task Groups
!	- tg_ivfftw: Inverse FFT driver for Task Groups
!=======================================================================================


!-----------------------------------------------------------------------
!=======================================================================
!      SUBROUTINE GROUPS (added by C. Bekas)
!      Define groups for task group parallilization
!=======================================================================
!-----------------------------------------------------------------------
SUBROUTINE GROUPS( nogrp_ , dffts )
   !------------
   !Modules used
   !------------

   USE mp_global,  ONLY : mpime, nproc, group, root
   USE mp_global,  ONLY : NOGRP, ME_OGRP, ME_PGRP  !Variables: NOGRP, MAXGRP, ME_OGRP, ME_PGRP
   USE mp,         ONLY : mp_bcast
   USE parameters, ONLY : MAXGRP
   USE io_global,  only : stdout
   USE fft_type,   only : fft_dlay_descriptor
   USE electrons_base, only: nspin
   USE parallel_include

   IMPLICIT NONE 

   INTEGER, INTENT(IN) :: nogrp_
   TYPE(fft_dlay_descriptor), INTENT(IN) :: dffts

#if defined (__MPI)

   !----------------------------------
   !Local Variables declaration
   !----------------------------------
   !NPROC:      Total number of processors
   !NPGRP:      Number of processors per group
   INTEGER  ::  MSGLEN, I, J, N1, LABEL, IPOS, WORLD, NEWGROUP
   INTEGER  ::  NPGRP, ios, IERR

   !--------------------------------------------------------------
   !Allocations
   !--------------------------------------------------------------
   ALLOCATE(NOLIST(MAXGRP))
   ALLOCATE(NPLIST(MAXGRP))
   

   tm_tg     = 0D0
   tm_rhoofr = 0D0

   !Find the number of processors and my rank
   SZ = NPROC

   !--------------------------------------------------------------
   !SUBDIVIDE THE PROCESSORS IN GROUPS
   !
   !THE NUMBER OF GROUPS HAS TO BE A DIVISOR OF THE NUMBER
   !OF PROCESSORS
   !--------------------------------------------------------------

   IF( MOD( nproc, nogrp_ ) /= 0 ) &
      CALL errore( " groups ", " nogrp should be a divisor of nproc ", 1 )
 
   ALLOCATE( PGROUP( NPROC ) )
   DO I=1, NPROC
      PGROUP(I) = I-1
   ENDDO

   nogrp = nogrp_

   allocate( nnrsx_vec( SZ ) )

   !Find maximum chunk of local data concerning coefficients of eigenfunctions in g-space

#if defined __MPI
   CALL MPI_Allgather(dffts%nnr, 1, MPI_INTEGER, nnrsx_vec, 1, MPI_INTEGER, group, IERR)
   strd = maxval( nnrsx_vec( 1:SZ ) )
#else
   strd = dffts%nnr 
#endif


   !---------------------------------------------------------------
   !Broadcast the number of groups: NOGRP
   !---------------------------------------------------------------

   CALL mp_bcast( nogrp, root, group )

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
   ALLOCATE(ALL_Z_STICKS(SZ))
   !
   !ALL-Gather number of Z-sticks from all processors
   !
#if defined __MPI
   CALL MPI_Allgather(dffts%nsw(mpime+1), 1, MPI_INTEGER, ALL_Z_STICKS, 1, MPI_INTEGER, group, IERR)
#else
   all_z_sticks( 1 ) = dffts%nsw( 1 )
#endif
   IF (.NOT.ALLOCATED(tmp_nsw)) ALLOCATE(tmp_nsw(SZ))
   IF (.NOT.ALLOCATED(tmp_npp)) THEN
      ALLOCATE(tmp_npp(SZ))
      tmp_npp(1)=-1
   ENDIF


   IF( NOGRP == 1 ) RETURN

   NPGRP = NPROC / NOGRP

   IF( NPGRP > MAXGRP ) THEN
      CALL errore( "groups", "too many npgrp", 1 )
   ENDIF

   !--------------------------------------
   !LIST OF PROCESSORS IN MY ORBITAL GROUP
   !--------------------------------------
   N1 = ( mpime / NOGRP ) * NOGRP - 1
   DO I = 1, NOGRP
      NOLIST( I ) = PGROUP( N1 + I + 1 )
      IF( mpime .EQ. NOLIST( I ) ) IPOS = I - 1
   ENDDO

   !-----------------------------------------
   !LIST OF PROCESSORS IN MY PLANE WAVE GROUP
   !-----------------------------------------
   DO I = 1, NPGRP
      NPLIST( I ) = PGROUP( IPOS + ( I - 1 ) * NOGRP + 1 )
   ENDDO

   !-----------------
   !SET UP THE GROUPS
   !-----------------
   DO I = 1, NPGRP
      IF( mpime .EQ. NPLIST( I ) ) LABEL = I
   ENDDO


   !---------------------------------------
   !CREATE ORBITAL GROUPS
   !---------------------------------------
#if defined __MPI
   CALL MPI_COMM_GROUP( group, WORLD, IERR )
   CALL MPI_GROUP_INCL( WORLD, NOGRP, NOLIST, NEWGROUP, IERR )
   CALL MPI_COMM_CREATE( group, NEWGROUP, ME_OGRP, IERR )
#endif


   DO I = 1, NOGRP
      IF( mpime .EQ. NOLIST( I ) ) LABEL = I + MAXCPU
   ENDDO


   !---------------------------------------
   !CREATE PLANEWAVE GROUPS
   !---------------------------------------
#if defined __MPI
   CALL MPI_COMM_GROUP(group, WORLD, IERR)
   CALL MPI_GROUP_INCL(WORLD, NPGRP, NPLIST, NEWGROUP, IERR)
   CALL MPI_COMM_CREATE(group, NEWGROUP, ME_PGRP, IERR)
#endif

   !--------
   !END
   !--------

#endif

   RETURN

   END SUBROUTINE GROUPS

SUBROUTINE GROUPS_NEW( nogrp_ , dffts )
   !------------
   !Modules used
   !------------
   USE mp_global,  ONLY : mpime, nproc, group
   USE mp_global,  ONLY : NOGRP, ME_OGRP, ME_PGRP  !Variables: NOGRP, MAXGRP, ME_OGRP, ME_PGRP
   USE parameters, ONLY : MAXGRP
   USE io_global,  only : stdout
   USE fft_type,   only : fft_dlay_descriptor
   USE electrons_base, only: nspin
   USE parallel_include

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: nogrp_
   TYPE(fft_dlay_descriptor), INTENT(IN) :: dffts

#if defined (__MPI)

   !----------------------------------
   !Local Variables declaration
   !----------------------------------
   !NPROC:      Total number of processors
   !NPGRP:      Number of processors per group
   INTEGER  ::  MSGLEN, I, J, N1, LABEL, IPOS, WORLD, NEWGROUP
   INTEGER  ::  NPGRP, ios, IERR
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: NOLIST_MATRIX, T_NOLIST_MATRIX
   INTEGER, DIMENSION(MPI_STATUS_SIZE)  :: STATUS

   INTEGER      tmp1(128)
   !--------------------------------------------------------------
   !Allocations
   !--------------------------------------------------------------
   ALLOCATE(NOLIST(MAXGRP))
   ALLOCATE(NPLIST(MAXGRP))


   tm_tg = 0D0
   tm_rhoofr = 0D0

   !Find the number of processors and my rank
   SZ = NPROC

   ALLOCATE(PGROUP(NPROC))
   DO I=1, NPROC
      PGROUP(I) = I-1
   ENDDO

   nogrp = nogrp_

   allocate(nnrsx_vec(SZ))
   !Find maximum chunk of local data concerning coefficients of eigenfunctions in g-space
   CALL MPI_Allgather(dffts%nnr, 1, MPI_INTEGER, nnrsx_vec, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
   strd = maxval(nnrsx_vec(1:SZ))

   !--------------------------------------------------------------
   !SUBDIVIDE THE PROCESSORS IN GROUPS
   !
   !THE NUMBER OF GROUPS HAS TO BE A DIVISOR OF THE NUMBER
   !OF PROCESSORS
   !--------------------------------------------------------------

   !---------------------------------------------------------------
   !Broadcast the number of groups: NOGRP
   !---------------------------------------------------------------
   CALL MPI_BCAST(NOGRP ,1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
   !Error check for broadcast
   IF (mpime.EQ.0) THEN
      IF (IERR.NE.0) THEN
         !Abort: Broadcast has failed
         CALL MPI_Abort !To be replaced by a proper exit routine
      ENDIF
   ENDIF

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


   ALLOCATE(ALL_Z_STICKS(SZ))
   !ALL-Gather number of Z-sticks from all processors
   CALL MPI_Allgather(dffts%nsw(mpime+1), 1, MPI_INTEGER, ALL_Z_STICKS, 1, MPI_INTEGER,  MPI_COMM_WORLD, IERR)
   IF (.NOT.ALLOCATED(tmp_nsw)) ALLOCATE(tmp_nsw(SZ))
   IF (.NOT.ALLOCATED(tmp_npp)) THEN
      ALLOCATE(tmp_npp(SZ))
      tmp_npp(1)=-1
   ENDIF



   IF(NOGRP.EQ.1) RETURN
   IF(NOGRP.GT.MAXGRP) THEN
      IF(mpime.EQ.0) THEN
         WRITE(stdout,*) ' MAXIMUM NUMBER OF GROUPS IS ',MAXGRP
         WRITE(stdout,*) ' NUMBER OF GROUPS     :',NOGRP
         CALL MPI_Abort !To be replaced by a proper exit routine
      ENDIF
   ENDIF
   IF(MOD(NPROC,NOGRP).NE.0) THEN
      IF(mpime.EQ.0) THEN
         WRITE(stdout,*) ' THE NUMBER OF GROUPS HAS TO BE A DIVISOR OF THE NUMBER OF PROCESSORS'
         WRITE(stdout,*) ' NUMBER OF PROCESSORS :',NPROC
         WRITE(stdout,*) ' NUMBER OF GROUPS     :',NOGRP
         CALL MPI_Abort !To be replaced by a proper exit routine
      ENDIF
   ENDIF
   NPGRP=NPROC/NOGRP
   IF(NPGRP.GT.MAXGRP) THEN
      IF(mpime.EQ.0) THEN
         WRITE(stdout,*) ' MINIMUM NUMBER OF GROUPS IS ',NPROC/MAXGRP
         WRITE(stdout,*) ' NUMBER OF GROUPS     :',NOGRP
         CALL MPI_Abort !To be replaced by a proper exit routine
      ENDIF
   ENDIF


   !--------------------------------------
   !LIST OF PROCESSORS IN MY ORBITAL GROUP
   !--------------------------------------
   !N1 = (mpime/NOGRP)*NOGRP-1
   !DO I=1,NOGRP
   !   NOLIST(I)=PGROUP(N1+I+1)
   !   IF(mpime.EQ.NOLIST(I)) IPOS=I-1
   !ENDDO

   IF (mpime.EQ.0) THEN
      ALLOCATE(NOLIST_MATRIX(NOGRP,NPROC/NOGRP), T_NOLIST_MATRIX(NPROC/NOGRP,NOGRP))
      DO I=1,NPROC/NOGRP
         DO J=1, NOGRP
            NOLIST_MATRIX(J,I) = PGROUP(I+(J-1)*NPROC/NOGRP)
            !PRINT *, NOLIST_MATRIX(J,I)
         ENDDO
      ENDDO
      DO I=1,NPROC/NOGRP 
         IF (I.EQ.1) THEN
            NOLIST(1:NOGRP) = NOLIST_MATRIX(1:NOGRP,I)
            DO J=2,NOGRP
               !PRINT *, NOLIST_MATRIX(J,I)
               CALL MPI_Send(NOLIST_MATRIX(1,I), NOGRP, MPI_INTEGER, NOLIST_MATRIX(J,I), 0, MPI_COMM_WORLD, IERR)
            ENDDO
         ELSE
            DO J=1,NOGRP
               CALL MPI_Send(NOLIST_MATRIX(1,I), NOGRP, MPI_INTEGER, NOLIST_MATRIX(J,I), 0, MPI_COMM_WORLD, IERR)
            ENDDO
         ENDIF
      ENDDO
   ELSE
      CALL MPI_Recv(NOLIST, NOGRP, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, STATUS, IERR)
   ENDIF
   T_NOLIST_MATRIX = transpose(NOLIST_MATRIX)
   IF (mpime.EQ.0) THEN
      NPLIST(1:NPROC/NOGRP) = NOLIST_MATRIX(1,1:NPROC/NOGRP)
      DO J=1, NOGRP
         DO I=1, NPROC/NOGRP
            IF ((I.NE.1).OR.(J.NE.1))  THEN
               CALL MPI_Send(T_NOLIST_MATRIX(1,J), NPROC/NOGRP, MPI_INTEGER, T_NOLIST_MATRIX(I,J), 0, MPI_COMM_WORLD, IERR)
            ENDIF
         ENDDO
      ENDDO
   ELSE
      CALL MPI_Recv(NPLIST, NPROC/NOGRP, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, STATUS, IERR)
   ENDIF


   IF (mpime.EQ.0) DEALLOCATE(NOLIST_MATRIX, T_NOLIST_MATRIX)

   CALL MPI_Barrier(MPI_COMM_WORLD, IERR)

   CALL MPI_gather(NOLIST, 2, MPI_INTEGER, tmp1, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
   DO I=1,NPROC*2
      PRINT *, tmp1(I)
   ENDDO
   !CALL MPI_ABORT

   !-----------------------------------------
   !LIST OF PROCESSORS IN MY PLANE WAVE GROUP
   !-----------------------------------------
   !DO I=1,NPGRP
   !   NPLIST(I)=PGROUP(IPOS+(I-1)*NOGRP+1)
   !ENDDO
   !-----------------
   !SET UP THE GROUPS
   !-----------------
   DO I=1,NPGRP
      IF(mpime.EQ.NPLIST(I)) LABEL=I
   ENDDO


   !---------------------------------------
   !CREATE ORBITAL GROUPS
   !---------------------------------------
   CALL MPI_COMM_GROUP(MPI_COMM_WORLD,WORLD,IERR)
   CALL MPI_GROUP_INCL(WORLD, NOGRP, NOLIST, NEWGROUP, IERR)
   CALL MPI_COMM_CREATE(MPI_COMM_WORLD,NEWGROUP, ME_OGRP, IERR)


   DO I=1,NOGRP
      IF(mpime.EQ.NOLIST(I)) LABEL=I+MAXCPU
   ENDDO


   !---------------------------------------
   !CREATE PLANEWAVE GROUPS
   !---------------------------------------
   CALL MPI_COMM_GROUP(MPI_COMM_WORLD, WORLD, IERR)
   CALL MPI_GROUP_INCL(WORLD, NPGRP, NPLIST, NEWGROUP, IERR)
   CALL MPI_COMM_CREATE(MPI_COMM_WORLD,NEWGROUP, ME_PGRP, IERR)

   !--------
   !END
   !--------

   PRINT *, "AFTER GROUPS"

#endif

   RETURN

   END SUBROUTINE GROUPS_NEW

END MODULE GROUPS_MODULE
