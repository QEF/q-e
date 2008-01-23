!
! Copyright (C) 2002-2003 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined __HPM
#  include "/cineca/prod/hpm/include/f_hpm.h"
#endif



!------------------------------------------------------------------------------!
    MODULE mp
!------------------------------------------------------------------------------!
      USE kinds,     ONLY : DP, i4b
      USE io_global, ONLY : stdout
      USE parallel_include
      !
      IMPLICIT NONE

      PUBLIC :: mp_start, mp_end, mp_env, &
        mp_bcast, mp_stop, mp_sum, mp_max, mp_min, mp_rank, mp_size, &
        mp_gather, mp_get, mp_put, mp_barrier, mp_report, mp_group_free, &
        mp_root_sum, mp_comm_free, mp_comm_create, mp_comm_group, mp_group_create, &
        mp_comm_split
!
      INTERFACE mp_bcast
#if defined __T3E
        MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, &
          mp_bcast_z, mp_bcast_zv, &
          mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
          mp_bcast_cm, mp_bcast_im, mp_bcast_it, mp_bcast_rt, mp_bcast_lv, &
          mp_bcast_lm, mp_bcast_r4d, mp_bcast_r5d, mp_bcast_ct, mp_bcast_c4d, &
          mp_bcast_i4b
#else
        MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, &
          mp_bcast_z, mp_bcast_zv, &
          mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
          mp_bcast_cm, mp_bcast_im, mp_bcast_it, mp_bcast_rt, mp_bcast_lv, &
          mp_bcast_lm, mp_bcast_r4d, mp_bcast_r5d, mp_bcast_ct,  mp_bcast_c4d
#endif
      END INTERFACE

      INTERFACE mp_sum
        MODULE PROCEDURE mp_sum_i1, mp_sum_iv, mp_sum_im, mp_sum_it, & 
          mp_sum_r1, mp_sum_rv, mp_sum_rm, mp_sum_rt, mp_sum_r4d, &
          mp_sum_c1, mp_sum_cv, mp_sum_cm, mp_sum_ct, mp_sum_c4d, &
          mp_sum_c6d, mp_sum_rmm, mp_sum_cmm
      END INTERFACE

      INTERFACE mp_root_sum
        MODULE PROCEDURE mp_root_sum_rm, mp_root_sum_cm
      END INTERFACE

      INTERFACE mp_get
        MODULE PROCEDURE mp_get_rv, mp_get_cv, mp_get_i1, mp_get_iv, &
          mp_get_rm
      END INTERFACE

      INTERFACE mp_put
        MODULE PROCEDURE mp_put_rv, mp_put_cv, mp_put_i1, mp_put_iv, &
          mp_put_rm
      END INTERFACE

      INTERFACE mp_max
        MODULE PROCEDURE mp_max_i, mp_max_r, mp_max_rv, mp_max_iv
      END INTERFACE
      INTERFACE mp_min
        MODULE PROCEDURE mp_min_i, mp_min_r, mp_min_rv, mp_min_iv
      END INTERFACE
      INTERFACE mp_gather
        MODULE PROCEDURE mp_gather_iv
      END INTERFACE

      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: mp_call_count(:)
      INTEGER, ALLOCATABLE, PRIVATE, SAVE :: mp_call_sizex(:)


      CHARACTER(LEN=80), PRIVATE :: err_msg = ' '

!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
!..mp_gather_iv
!..Carlo Cavazzoni
      SUBROUTINE mp_gather_iv(mydata, alldata, root, gid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: mydata(:), root
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER, INTENT(OUT) :: alldata(:,:)
        INTEGER :: msglen, ierr


#if defined (__MPI)
        msglen = SIZE(mydata)
        IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop( 8002 )
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL MPI_GATHER(mydata, msglen, MPI_INTEGER, alldata, msglen, MPI_INTEGER, root, group, IERR)
        IF (ierr/=0) CALL mp_stop( 8003 )
#else
        msglen = SIZE(mydata)
        IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop( 8004 )
        alldata(:,1) = mydata(:)
#endif
        RETURN
      END SUBROUTINE mp_gather_iv

!
!------------------------------------------------------------------------------!
!..mp_start
      SUBROUTINE mp_start

! ...
        IMPLICIT NONE
        INTEGER :: ierr, taskid
! ...
        ierr = 0
        taskid = 0

        ALLOCATE( mp_call_count( 1000 ) )
        mp_call_count = 0
        ALLOCATE( mp_call_sizex( 1000 ) )
        mp_call_sizex = 0

#if defined(__MPI) || defined (__SHMEM)
        CALL MPI_INIT(ierr)
        IF (ierr/=0) CALL mp_stop( 8005 )
#endif

#if defined __HPM 

        !   initialize the IBM Harware performance monitor
      
#  if defined(__MPI) || defined (__SHMEM)
        CALL mpi_comm_rank( mpi_comm_world, taskid, ierr)
#  endif
        CALL f_hpminit( taskid, 'profiling' )
#endif
! ...

      END SUBROUTINE mp_start
!
!------------------------------------------------------------------------------!
!..mp_end

      SUBROUTINE mp_end
        IMPLICIT NONE
        INTEGER :: ierr, taskid

        ierr = 0
        taskid = 0

        IF ( ALLOCATED ( mp_call_count ) ) DEALLOCATE( mp_call_count )
        IF ( ALLOCATED ( mp_call_sizex ) ) DEALLOCATE( mp_call_sizex )

#if defined __HPM 

        !   terminate the IBM Harware performance monitor

#  if defined(__MPI)
        CALL mpi_comm_rank( mpi_comm_world, taskid, ierr)
#  endif
        CALL f_hpmterminate( taskid )
#endif

#if defined(__MPI)
        CALL mpi_finalize(ierr)
        IF (ierr/=0) CALL mp_stop( 8006 )
#endif
        RETURN
      END SUBROUTINE mp_end
!
!------------------------------------------------------------------------------!
!..mp_env

      SUBROUTINE mp_env(numtask, taskid, groupid)
        IMPLICIT NONE
        INTEGER, INTENT (OUT) :: numtask, taskid, groupid
        INTEGER :: ierr

        ierr = 0
        numtask = 1
        taskid = 0
        groupid = 0

#if defined(__MPI)

        CALL mpi_comm_rank(mpi_comm_world,taskid,ierr)
        IF (ierr/=0) CALL mp_stop( 8007 )
        CALL mpi_comm_size(mpi_comm_world,numtask,ierr)
        groupid = mpi_comm_world
        IF (ierr/=0) CALL mp_stop( 8008 )

#endif

        RETURN
      END SUBROUTINE mp_env

!------------------------------------------------------------------------------!
!..mp_group

      SUBROUTINE mp_comm_group( comm, group )
         IMPLICIT NONE
         INTEGER, INTENT (IN) :: comm
         INTEGER, INTENT (OUT) :: group
         INTEGER :: ierr
         ierr = 0
#if defined(__MPI)
         CALL mpi_comm_group( comm, group, ierr )
         IF (ierr/=0) CALL mp_stop( 8009 )
#else
         group = 0
#endif
      END SUBROUTINE  mp_comm_group

      SUBROUTINE mp_comm_split( old_comm, color, key, new_comm )
         IMPLICIT NONE
         INTEGER, INTENT (IN) :: old_comm
         INTEGER, INTENT (IN) :: color, key
         INTEGER, INTENT (OUT) :: new_comm
         INTEGER :: ierr
         ierr = 0
#if defined(__MPI)
         CALL MPI_COMM_SPLIT( old_comm, color, key, new_comm, ierr )
         IF (ierr/=0) CALL mp_stop( 8009 )
#else
         new_comm = old_comm
#endif  
      END SUBROUTINE  mp_comm_split


      SUBROUTINE mp_group_create( group_list, group_size, old_grp, new_grp )
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: group_list(:), group_size, old_grp
        INTEGER, INTENT (OUT) :: new_grp
        INTEGER :: ierr

        ierr = 0
        new_grp = old_grp
#if defined(__MPI)
        CALL mpi_group_incl( old_grp, group_size, group_list, new_grp, ierr )
        IF (ierr/=0) CALL mp_stop( 8010 )
#endif
      END SUBROUTINE mp_group_create

!------------------------------------------------------------------------------!
      SUBROUTINE mp_comm_create( old_comm, new_grp, new_comm )
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: old_comm
        INTEGER, INTENT (IN) :: new_grp
        INTEGER, INTENT (OUT) :: new_comm
        INTEGER :: ierr

        ierr = 0
        new_comm = old_comm
#if defined(__MPI)
        CALL mpi_comm_create( old_comm, new_grp, new_comm, ierr )
        IF (ierr/=0) CALL mp_stop( 8011 )
#endif
      END SUBROUTINE mp_comm_create

!------------------------------------------------------------------------------!
!..mp_group_free
      SUBROUTINE mp_group_free( group )
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: group
        INTEGER :: ierr
        ierr = 0
#if defined(__MPI)
        CALL mpi_group_free( group, ierr )
        IF (ierr/=0) CALL mp_stop( 8014 )
#endif
      END SUBROUTINE mp_group_free
!------------------------------------------------------------------------------!

      SUBROUTINE mp_comm_free( comm )
         IMPLICIT NONE
         INTEGER, INTENT (INOUT) :: comm
         INTEGER :: ierr
         ierr = 0
#if defined(__MPI)
         IF( comm /= MPI_COMM_NULL ) THEN
            CALL mpi_comm_free( comm, ierr )
            IF (ierr/=0) CALL mp_stop( 8013 )
         END IF
#endif
         RETURN
      END SUBROUTINE mp_comm_free

!------------------------------------------------------------------------------!
!..mp_bcast

#if defined (__T3E)

      SUBROUTINE mp_bcast_i4b(msg,source,gid)
        IMPLICIT NONE
        INTEGER(i4b) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: source
        INTEGER :: msglen, ierr, imsg

#if defined(__MPI)
        ierr = 0
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        imsg = msg
        CALL mpi_bcast(imsg, msglen, mpi_integer, source, group, ierr)
        msg = imsg
        IF (ierr/=0) CALL mp_stop( 8019 )
#endif
      END SUBROUTINE mp_bcast_i4b

#endif


!------------------------------------------------------------------------------!
!..mp_bcast

      SUBROUTINE mp_bcast_i1(msg,source,gid)
        IMPLICIT NONE
        INTEGER :: msg
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen

#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL BCAST_INTEGER( msg, msglen, source, group )
        mp_call_count( 1 ) = mp_call_count( 1 ) + 1
        mp_call_sizex( 1 ) = MAX( mp_call_sizex( 1 ), msglen )
#endif
      END SUBROUTINE mp_bcast_i1
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_iv(msg,source,gid)
        IMPLICIT NONE
        INTEGER :: msg(:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL BCAST_INTEGER( msg, msglen, source, group )
        mp_call_count( 2 ) = mp_call_count( 2 ) + 1
        mp_call_sizex( 2 ) = MAX( mp_call_sizex( 2 ), msglen )
#endif
      END SUBROUTINE mp_bcast_iv
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_im( msg, source, gid )
        IMPLICIT NONE
        INTEGER :: msg(:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL BCAST_INTEGER( msg, msglen, source, group )
        mp_call_count( 3 ) = mp_call_count( 3 ) + 1
        mp_call_sizex( 3 ) = MAX( mp_call_sizex( 3 ), msglen )
#endif
      END SUBROUTINE mp_bcast_im
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_it(msg,source,gid)
        IMPLICIT NONE
        INTEGER :: msg(:,:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL BCAST_INTEGER( msg, msglen, source, group )
        mp_call_count( 4 ) = mp_call_count( 4 ) + 1
        mp_call_sizex( 4 ) = MAX( mp_call_sizex( 4 ), msglen )
#endif
      END SUBROUTINE mp_bcast_it
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_r1(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg
        INTEGER :: msglen, source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, msglen, source, group )
        mp_call_count( 5 ) = mp_call_count( 5 ) + 1
        mp_call_sizex( 5 ) = MAX( mp_call_sizex( 5 ), msglen )
#endif
      END SUBROUTINE mp_bcast_r1
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rv(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen

#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, msglen, source, group )
        mp_call_count( 6 ) = mp_call_count( 6 ) + 1
        mp_call_sizex( 6 ) = MAX( mp_call_sizex( 6 ), msglen )
#endif
      END SUBROUTINE mp_bcast_rv
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rm(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, msglen, source, group )
        mp_call_count( 7 ) = mp_call_count( 7 ) + 1
        mp_call_sizex( 7 ) = MAX( mp_call_sizex( 7 ), msglen )
#endif
      END SUBROUTINE mp_bcast_rm
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_rt(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, msglen, source, group )
        mp_call_count( 8 ) = mp_call_count( 8 ) + 1
        mp_call_sizex( 8 ) = MAX( mp_call_sizex( 8 ), msglen )
#endif
      END SUBROUTINE mp_bcast_rt
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_r4d(msg, source, gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:,:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, msglen, source, group )
        mp_call_count( 9 ) = mp_call_count( 9 ) + 1
        mp_call_sizex( 9 ) = MAX( mp_call_sizex( 9 ), msglen )
#endif
      END SUBROUTINE mp_bcast_r4d

!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_r5d(msg, source, gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:,:,:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, msglen, source, group )
        mp_call_count( 9 ) = mp_call_count( 9 ) + 1
        mp_call_sizex( 9 ) = MAX( mp_call_sizex( 9 ), msglen )
#endif
      END SUBROUTINE mp_bcast_r5d

!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_c1(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, 2 * msglen, source, group )
        mp_call_count( 10 ) = mp_call_count( 10 ) + 1
        mp_call_sizex( 10 ) = MAX( mp_call_sizex( 10 ), msglen )
#endif
      END SUBROUTINE mp_bcast_c1
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_cv(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, 2 * msglen, source, group )
        mp_call_count( 11 ) = mp_call_count( 11 ) + 1
        mp_call_sizex( 11 ) = MAX( mp_call_sizex( 11 ), msglen )
#endif
      END SUBROUTINE mp_bcast_cv
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_cm(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, 2 * msglen, source, group )
        mp_call_count( 12 ) = mp_call_count( 12 ) + 1
        mp_call_sizex( 12 ) = MAX( mp_call_sizex( 12 ), msglen )
#endif
      END SUBROUTINE mp_bcast_cm
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_ct(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, 2 * msglen, source, group )
        mp_call_count( 13 ) = mp_call_count( 13 ) + 1
        mp_call_sizex( 13 ) = MAX( mp_call_sizex( 13 ), msglen )
#endif
      END SUBROUTINE mp_bcast_ct

!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_c4d(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_real( msg, 2 * msglen, source, group )
        mp_call_count( 14 ) = mp_call_count( 14 ) + 1
        mp_call_sizex( 14 ) = MAX( mp_call_sizex( 14 ), msglen )
#endif
      END SUBROUTINE mp_bcast_c4d

!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_bcast_l(msg,source,gid)
        IMPLICIT NONE
        LOGICAL :: msg
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_logical( msg, msglen, source, group )
        mp_call_count( 15 ) = mp_call_count( 15 ) + 1
        mp_call_sizex( 15 ) = MAX( mp_call_sizex( 15 ), msglen )
#endif
      END SUBROUTINE mp_bcast_l
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_lv(msg,source,gid)
        IMPLICIT NONE
        LOGICAL :: msg(:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_logical( msg, msglen, source, group )
        mp_call_count( 16 ) = mp_call_count( 16 ) + 1
        mp_call_sizex( 16 ) = MAX( mp_call_sizex( 16 ), msglen )
#endif
      END SUBROUTINE mp_bcast_lv

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_lm(msg,source,gid)
        IMPLICIT NONE
        LOGICAL :: msg(:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL bcast_logical( msg, msglen, source, group )
        mp_call_count( 17 ) = mp_call_count( 17 ) + 1
        mp_call_sizex( 17 ) = MAX( mp_call_sizex( 17 ), msglen )
#endif
      END SUBROUTINE mp_bcast_lm


!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_z(msg,source,gid)
        IMPLICIT NONE
        CHARACTER (len=*) :: msg
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr, i
        INTEGER, ALLOCATABLE :: imsg(:)
#if defined(__MPI)
        ierr = 0
        msglen = len(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        IF (ierr/=0) CALL mp_stop( 8051 )
        ALLOCATE (imsg(1:msglen), STAT=ierr)
        IF (ierr/=0) CALL mp_stop( 8052 )
        DO i = 1, msglen
          imsg(i) = ichar(msg(i:i))
        END DO
        CALL bcast_integer( imsg, msglen, source, group )
        DO i = 1, msglen
          msg(i:i) = char(imsg(i))
        END DO
        DEALLOCATE (imsg, STAT=ierr)
        IF (ierr/=0) CALL mp_stop( 8054 )
        mp_call_count( 18 ) = mp_call_count( 18 ) + 1
        mp_call_sizex( 18 ) = MAX( mp_call_sizex( 18 ), msglen )
#endif
      END SUBROUTINE mp_bcast_z
!
!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_zv(msg,source,gid)
        IMPLICIT NONE
        CHARACTER (len=*) :: msg(:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, m1, m2, ierr, i, j
        INTEGER, ALLOCATABLE :: imsg(:,:)
#if defined(__MPI)
        ierr = 0
        m1 = LEN(msg)
        m2 = SIZE(msg)
        msglen = LEN(msg)*SIZE(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        ALLOCATE (imsg(1:m1,1:m2), STAT=ierr)
        IF (ierr/=0) CALL mp_stop( 8057 )
        DO j = 1, m2
          DO i = 1, m1
            imsg(i,j) = ichar(msg(j)(i:i))
          END DO
        END DO
        CALL bcast_integer( imsg, msglen, source, group )
        DO j = 1, m2
          DO i = 1, m1
            msg(j)(i:i) = char(imsg(i,j))
          END DO
        END DO
        DEALLOCATE (imsg, STAT=ierr)
        IF (ierr/=0) CALL mp_stop( 8059 )
        mp_call_count( 19 ) = mp_call_count( 19 ) + 1
        mp_call_sizex( 19 ) = MAX( mp_call_sizex( 19 ), msglen )
#endif
      END SUBROUTINE mp_bcast_zv
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_i1(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        INTEGER :: msg_dest, msg_sour
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen = 1

#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif

        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF(dest .NE. sour) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, msglen, MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8060 )
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, msglen, MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8061 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8062 )
           END IF
#endif
        ELSE
          msg_dest = msg_sour
        END IF

#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8063 )
#endif

        mp_call_count( 20 ) = mp_call_count( 20 ) + 1
        mp_call_sizex( 20 ) = MAX( mp_call_sizex( 20 ), msglen )

        RETURN
      END SUBROUTINE mp_get_i1

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_iv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        INTEGER :: msg_dest(:), msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif

        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen = SIZE(msg_sour)
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8064 )
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8065 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8066 )
             msglen = nrcv
           END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8068 )
#endif
        mp_call_count( 21 ) = mp_call_count( 21 ) + 1
        mp_call_sizex( 21 ) = MAX( mp_call_sizex( 21 ), msglen )
        RETURN
      END SUBROUTINE mp_get_iv

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_rv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        REAL (DP) :: msg_dest(:), msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif

        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen = SIZE(msg_sour) 
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8069 )
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8070 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8071 )
             msglen = nrcv
           END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8073 )
#endif
        mp_call_count( 22 ) = mp_call_count( 22 ) + 1
        mp_call_sizex( 22 ) = MAX( mp_call_sizex( 22 ), msglen )
        RETURN
      END SUBROUTINE mp_get_rv

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_rm(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        REAL (DP) :: msg_dest(:,:), msg_sour(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif

        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8074 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8075 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8076 )
             msglen = nrcv
           END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour,1), 1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE( msg_sour )
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8078 )
#endif
        mp_call_count( 23 ) = mp_call_count( 23 ) + 1
        mp_call_sizex( 23 ) = MAX( mp_call_sizex( 23 ), msglen )
        RETURN
      END SUBROUTINE mp_get_rm


!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_cv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        COMPLEX (DP) :: msg_dest(:), msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif

        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF( dest .NE. sour ) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8079 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8080 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8081 )
             msglen = nrcv
           END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8083 )
#endif
        mp_call_count( 24 ) = mp_call_count( 24 ) + 1
        mp_call_sizex( 24 ) = MAX( mp_call_sizex( 24 ), msglen )
        RETURN
      END SUBROUTINE mp_get_cv
!------------------------------------------------------------------------------!
!
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_put_i1(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        INTEGER :: msg_dest, msg_sour
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif

        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF(dest .NE. sour) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, 1, MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8084 )
             msglen = 1
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, 1, MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8085 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8086 )
             msglen = 1
           END IF
#endif
        ELSE
          msg_dest = msg_sour
          msglen = 1
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8088 )
#endif
        mp_call_count( 25 ) = mp_call_count( 25 ) + 1
        mp_call_sizex( 25 ) = MAX( mp_call_sizex( 25 ), msglen )
        RETURN
      END SUBROUTINE mp_put_i1

!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_iv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        INTEGER :: msg_dest(:), msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif
        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8089 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8090 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8091 )
             msglen = nrcv
           END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8093 )
#endif
        mp_call_count( 26 ) = mp_call_count( 26 ) + 1
        mp_call_sizex( 26 ) = MAX( mp_call_sizex( 26 ), msglen )
        RETURN
      END SUBROUTINE mp_put_iv

!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_rv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        REAL (DP) :: msg_dest(:), msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif
        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8094 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8095 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8096 )
             msglen = nrcv
           END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8098 )
#endif
        mp_call_count( 27 ) = mp_call_count( 27 ) + 1
        mp_call_sizex( 27 ) = MAX( mp_call_sizex( 27 ), msglen )
        RETURN
      END SUBROUTINE mp_put_rv

!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_rm(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        REAL (DP) :: msg_dest(:,:), msg_sour(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif
        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8099 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8100 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8101 )
             msglen = nrcv
           END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour,1),1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8103 )
#endif
        mp_call_count( 28 ) = mp_call_count( 28 ) + 1
        mp_call_sizex( 28 ) = MAX( mp_call_sizex( 28 ), msglen )
        RETURN
      END SUBROUTINE mp_put_rm


!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_cv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        COMPLEX (DP) :: msg_dest(:), msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
#endif
        ! processors not taking part in the communication have 0 lenght message

        msglen = 0

        IF( dest .NE. sour ) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8104 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8105 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8106 )
             msglen = nrcv
           END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8108 )
#endif
        mp_call_count( 29 ) = mp_call_count( 29 ) + 1
        mp_call_sizex( 29 ) = MAX( mp_call_sizex( 29 ), msglen )
        RETURN
      END SUBROUTINE mp_put_cv

!
!------------------------------------------------------------------------------!
!
!..mp_stop
!
      SUBROUTINE mp_stop(code)
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: code
        WRITE( stdout, fmt='( "*** error in Message Passing (mp) module ***")' )
        WRITE( stdout, fmt='( "*** error msg:  ",A60)' ) TRIM( err_msg )
        WRITE( stdout, fmt='( "*** error code: ",I5)' ) code
#if defined(__MPI)
        CALL mpi_abort(mpi_comm_world,code)
#endif
        STOP
      END SUBROUTINE mp_stop
!------------------------------------------------------------------------------!
!
!..mp_sum
      SUBROUTINE mp_sum_i1(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_integer( msglen, msg, group, -1 )
        mp_call_count( 30 ) = mp_call_count( 30 ) + 1
        mp_call_sizex( 30 ) = MAX( mp_call_sizex( 30 ), msglen )
#endif
      END SUBROUTINE mp_sum_i1
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, group, -1 )
        mp_call_count( 31 ) = mp_call_count( 31 ) + 1
        mp_call_sizex( 31 ) = MAX( mp_call_sizex( 31 ), msglen )
#endif
      END SUBROUTINE mp_sum_iv
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_im(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, group, -1 )
        mp_call_count( 32 ) = mp_call_count( 32 ) + 1
        mp_call_sizex( 32 ) = MAX( mp_call_sizex( 32 ), msglen )
#endif
      END SUBROUTINE mp_sum_im
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_it(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, group, -1 )
        mp_call_count( 33 ) = mp_call_count( 33 ) + 1
        mp_call_sizex( 33 ) = MAX( mp_call_sizex( 33 ), msglen )
#endif
      END SUBROUTINE mp_sum_it

!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_r1(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( msglen, msg, group, -1 )
        mp_call_count( 34 ) = mp_call_count( 34 ) + 1
        mp_call_sizex( 34 ) = MAX( mp_call_sizex( 34 ), msglen )
#endif
      END SUBROUTINE mp_sum_r1

!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( msglen, msg, group, -1 )
        mp_call_count( 35 ) = mp_call_count( 35 ) + 1
        mp_call_sizex( 35 ) = MAX( mp_call_sizex( 35 ), msglen )
#endif
      END SUBROUTINE mp_sum_rv
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_rm(msg, gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( msglen, msg, group, -1 )
        mp_call_count( 36 ) = mp_call_count( 36 ) + 1
        mp_call_sizex( 36 ) = MAX( mp_call_sizex( 36 ), msglen )
#endif
      END SUBROUTINE mp_sum_rm


      SUBROUTINE mp_root_sum_rm( msg, res, root, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (IN)  :: msg(:,:)
        REAL (DP), INTENT (OUT) :: res(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr, taskid

#if defined(__MPI)

        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid

        CALL mpi_comm_rank( group, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 8129 )
        !
        IF( taskid == root ) THEN
           IF( msglen > size(res) ) CALL mp_stop( 8129 )
        END IF

        CALL reduce_base_real_to( msglen, msg, res, group, root )

        mp_call_count( 36 ) = mp_call_count( 36 ) + 1
        mp_call_sizex( 36 ) = MAX( mp_call_sizex( 36 ), msglen )

#else

        res = msg

#endif

      END SUBROUTINE mp_root_sum_rm


      SUBROUTINE mp_root_sum_cm( msg, res, root, gid )
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN)  :: msg(:,:)
        COMPLEX (DP), INTENT (OUT) :: res(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr, taskid

#if defined(__MPI)

        msglen = size(msg)

        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid

        CALL mpi_comm_rank( group, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 8129 )

        IF( taskid == root ) THEN
           IF( msglen > size(res) ) CALL mp_stop( 8129 )
        END IF

        CALL reduce_base_real_to( 2 * msglen, msg, res, group, root )

        mp_call_count( 36 ) = mp_call_count( 36 ) + 1
        mp_call_sizex( 36 ) = MAX( mp_call_sizex( 36 ), msglen )

#else

        res = msg

#endif

      END SUBROUTINE mp_root_sum_cm

!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!

      SUBROUTINE mp_sum_rmm( msg, res, root, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (IN) :: msg(:,:)
        REAL (DP), INTENT (OUT) :: res(:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: root
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
        INTEGER :: taskid, ierr

        msglen = size(msg)

#if defined(__MPI)

        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid

        IF( PRESENT( root ) ) THEN
           !
           CALL mpi_comm_rank( group, taskid, ierr)
           IF( ierr /= 0 ) CALL mp_stop( 8129 )

           IF( taskid == root ) THEN
              IF( msglen > size(res) ) CALL mp_stop( 8129 )
           END IF
           !
           CALL reduce_base_real_to( msglen, msg, res, group, root )
           !
        ELSE
           !
           IF( msglen > size(res) ) CALL mp_stop( 8129 )
           !
           CALL reduce_base_real_to( msglen, msg, res, group, -1 )
           !
        END IF

        mp_call_count( 37 ) = mp_call_count( 37 ) + 1
        mp_call_sizex( 37 ) = MAX( mp_call_sizex( 37 ), msglen )

#else
        res = msg
#endif

      END SUBROUTINE mp_sum_rmm


!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_rt( msg, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( msglen, msg, group, -1 )
        mp_call_count( 38 ) = mp_call_count( 38 ) + 1
        mp_call_sizex( 38 ) = MAX( mp_call_sizex( 38 ), msglen )
#endif
      END SUBROUTINE mp_sum_rt

!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_r4d(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( msglen, msg, group, -1 )
        mp_call_count( 39 ) = mp_call_count( 39 ) + 1
        mp_call_sizex( 39 ) = MAX( mp_call_sizex( 39 ), msglen )
#endif
      END SUBROUTINE mp_sum_r4d



!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_c1(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen

#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( 2 * msglen, msg, group, -1 )
        mp_call_count( 40 ) = mp_call_count( 40 ) + 1
        mp_call_sizex( 40 ) = MAX( mp_call_sizex( 40 ), msglen )
#endif
      END SUBROUTINE mp_sum_c1
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_cv(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( 2 * msglen, msg, group, -1 )
        mp_call_count( 41 ) = mp_call_count( 41 ) + 1
        mp_call_sizex( 41 ) = MAX( mp_call_sizex( 41 ), msglen )
#endif
      END SUBROUTINE mp_sum_cv
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_cm(msg, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( 2 * msglen, msg, group, -1 )
        mp_call_count( 42 ) = mp_call_count( 42 ) + 1
        mp_call_sizex( 42 ) = MAX( mp_call_sizex( 42 ), msglen )
#endif
      END SUBROUTINE mp_sum_cm
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_cmm(msg, res, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN) :: msg(:,:)
        COMPLEX (DP), INTENT (OUT) :: res(:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real_to( 2 * msglen, msg, res, group, -1 )
        mp_call_count( 43 ) = mp_call_count( 43 ) + 1
        mp_call_sizex( 43 ) = MAX( mp_call_sizex( 43 ), msglen )
#else
        res = msg
#endif
      END SUBROUTINE mp_sum_cmm


!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_ct(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = SIZE(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( 2 * msglen, msg, group, -1 )
        mp_call_count( 44 ) = mp_call_count( 44 ) + 1
        mp_call_sizex( 44 ) = MAX( mp_call_sizex( 44 ), msglen )
#endif
      END SUBROUTINE mp_sum_ct

!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c4d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( 2 * msglen, msg, group, -1 )
        mp_call_count( 45 ) = mp_call_count( 45 ) + 1
        mp_call_sizex( 45 ) = MAX( mp_call_sizex( 45 ), msglen )
#endif
      END SUBROUTINE mp_sum_c4d

!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c6d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL reduce_base_real( 2 * msglen, msg, group, -1 )
        mp_call_count( 45 ) = mp_call_count( 45 ) + 1
        mp_call_sizex( 45 ) = MAX( mp_call_sizex( 45 ), msglen )
#endif
      END SUBROUTINE mp_sum_c6d



!------------------------------------------------------------------------------!
      SUBROUTINE mp_max_i(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL parallel_max_integer( msglen, msg, group, -1 )
        mp_call_count( 46 ) = mp_call_count( 46 ) + 1
        mp_call_sizex( 46 ) = MAX( mp_call_sizex( 46 ), msglen )
#endif
      END SUBROUTINE mp_max_i
!
!------------------------------------------------------------------------------!
!
!..mp_max_iv
!..Carlo Cavazzoni
!
      SUBROUTINE mp_max_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        CALL parallel_max_integer( msglen, msg, group, -1 )
        mp_call_count( 47 ) = mp_call_count( 47 ) + 1
        mp_call_sizex( 47 ) = MAX( mp_call_sizex( 47 ), msglen )
#endif
      END SUBROUTINE mp_max_iv
!
!----------------------------------------------------------------------

      SUBROUTINE mp_max_r(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = 1
        CALL parallel_max_real( msglen, msg, group, -1 )
        mp_call_count( 48 ) = mp_call_count( 48 ) + 1
        mp_call_sizex( 48 ) = MAX( mp_call_sizex( 48 ), msglen )
#endif
      END SUBROUTINE mp_max_r
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_max_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        CALL parallel_max_real( msglen, msg, group, -1 )
        mp_call_count( 49 ) = mp_call_count( 49 ) + 1
        mp_call_sizex( 49 ) = MAX( mp_call_sizex( 49 ), msglen )
#endif
      END SUBROUTINE mp_max_rv
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_i(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = 1
        CALL parallel_min_integer( msglen, msg, group, -1 )
        mp_call_count( 50 ) = mp_call_count( 50 ) + 1
        mp_call_sizex( 50 ) = MAX( mp_call_sizex( 50 ), msglen )
#endif
      END SUBROUTINE mp_min_i
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = SIZE(msg)
        CALL parallel_min_integer( msglen, msg, group, -1 )
        mp_call_count( 51 ) = mp_call_count( 51 ) + 1
        mp_call_sizex( 51 ) = MAX( mp_call_sizex( 51 ), msglen )
#endif
      END SUBROUTINE mp_min_iv
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_r(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = 1
        CALL parallel_min_real( msglen, msg, group, -1 )
        mp_call_count( 52 ) = mp_call_count( 52 ) + 1
        mp_call_sizex( 52 ) = MAX( mp_call_sizex( 52 ), msglen )
#endif
      END SUBROUTINE mp_min_r
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        CALL parallel_min_real( msglen, msg, group, -1 )
        mp_call_count( 53 ) = mp_call_count( 53 ) + 1
        mp_call_sizex( 53 ) = MAX( mp_call_sizex( 53 ), msglen )
#endif
      END SUBROUTINE mp_min_rv

!------------------------------------------------------------------------------!

      SUBROUTINE mp_barrier(gid)
        IMPLICIT NONE
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL MPI_BARRIER(group,IERR)
        IF (ierr/=0) CALL mp_stop( 8180 )
#endif
      END SUBROUTINE mp_barrier

!------------------------------------------------------------------------------!
!.. Carlo Cavazzoni
!..mp_rank
      FUNCTION mp_rank( comm )
        IMPLICIT NONE
        INTEGER :: mp_rank
        INTEGER, OPTIONAL, INTENT(IN) :: comm
        INTEGER :: ierr, taskid

        ierr = 0
        taskid = 0
#if defined(__MPI)
        IF( PRESENT( comm ) ) THEN
           CALL mpi_comm_rank(comm,taskid,ierr)
        ELSE
           CALL mpi_comm_rank(mpi_comm_world,taskid,ierr)
        END IF
        IF (ierr/=0) CALL mp_stop( 8181 )
#endif
        mp_rank = taskid
      END FUNCTION mp_rank

!------------------------------------------------------------------------------!
!.. Carlo Cavazzoni
!..mp_size
      FUNCTION mp_size( comm )
        IMPLICIT NONE
        INTEGER :: mp_size
        INTEGER, OPTIONAL, INTENT(IN) :: comm
        INTEGER :: ierr, numtask

        ierr = 0
        numtask = 1
#if defined(__MPI)
        IF( PRESENT( comm ) ) THEN
           CALL mpi_comm_size(comm,numtask,ierr)
        ELSE
           CALL mpi_comm_size(mpi_comm_world,numtask,ierr)
        END IF
        IF (ierr/=0) CALL mp_stop( 8182 )
#endif
        mp_size = numtask
      END FUNCTION mp_size

      SUBROUTINE mp_report
        INTEGER :: i
        WRITE( stdout, *)
#if defined(__MPI)
#  if defined (__MP_STAT)
        WRITE( stdout, 20 )
        DO i = 1, SIZE( mp_call_count )
           IF( mp_call_count( i ) > 0 ) THEN
              WRITE( stdout, 30 ) i, mp_call_count( i ),  mp_call_sizex( i )
           END IF
        END DO
#  endif
10      FORMAT(3X,'Message Passing, maximum message size (bytes) : ',I15)
20      FORMAT(3X,'Sub.   calls   maxsize')
30      FORMAT(3X,I4,I8,I10)
#else
        WRITE( stdout, *) 
#endif
        RETURN
      END SUBROUTINE mp_report


!------------------------------------------------------------------------------!
    END MODULE mp
!------------------------------------------------------------------------------!

