!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains interfaces to most low-level MPI operations:
! initialization and stopping, broadcast, parallel sum, etc.
!
!------------------------------------------------------------------------------!
    MODULE mp
!------------------------------------------------------------------------------!
      USE util_param,     ONLY : DP, stdout
      USE parallel_include
#if defined(__CUDA)
      USE cudafor,        ONLY: cudamemcpy, cudamemcpy2d, &
                                & cudaMemcpyDeviceToDevice, &
                                & cudaDeviceSynchronize
#endif
      !
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: mp_start, mp_abort, mp_stop, mp_end, &
        mp_bcast, mp_sum, mp_max, mp_min, mp_rank, mp_size, &
        mp_gather, mp_alltoall, mp_get, mp_put, &
        mp_barrier, mp_report, mp_group_free, &
        mp_root_sum, mp_comm_free, mp_comm_create, mp_comm_group, &
        mp_group_create, mp_comm_split, mp_set_displs, &
        mp_circular_shift_left, &
        mp_get_comm_null, mp_get_comm_self, mp_count_nodes, &
        mp_type_create_column_section, mp_type_free, &
        mp_allgather

!
      INTERFACE mp_bcast
        MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, &
          mp_bcast_z, mp_bcast_zv, &
          mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
          mp_bcast_cm, mp_bcast_im, mp_bcast_it, mp_bcast_i4d, mp_bcast_rt, mp_bcast_lv, &
          mp_bcast_lm, mp_bcast_r4d, mp_bcast_r5d, mp_bcast_ct,  mp_bcast_c4d,&
          mp_bcast_c5d
#if defined(__CUDA)
        ! this code duplication can be done with regex: rexp = re.compile(r"(\b\w+\b)") ; rexp.sub(r'\1_gpu',code_here)
        MODULE PROCEDURE mp_bcast_i1_gpu, mp_bcast_r1_gpu, mp_bcast_c1_gpu, &
          !mp_bcast_z_gpu, mp_bcast_zv_gpu, &
          mp_bcast_iv_gpu, mp_bcast_rv_gpu, mp_bcast_cv_gpu, mp_bcast_l_gpu, mp_bcast_rm_gpu, &
          mp_bcast_cm_gpu, mp_bcast_im_gpu, mp_bcast_it_gpu, mp_bcast_i4d_gpu, mp_bcast_rt_gpu, mp_bcast_lv_gpu, &
          mp_bcast_lm_gpu, mp_bcast_r4d_gpu, mp_bcast_r5d_gpu, mp_bcast_ct_gpu,  mp_bcast_c4d_gpu,&
          mp_bcast_c5d_gpu
#endif
      END INTERFACE

      INTERFACE mp_sum
        MODULE PROCEDURE mp_sum_i1, mp_sum_iv, mp_sum_im, mp_sum_it, &
          mp_sum_r1, mp_sum_rv, mp_sum_rm, mp_sum_rt, mp_sum_r4d, &
          mp_sum_c1, mp_sum_cv, mp_sum_cm, mp_sum_ct, mp_sum_c4d, &
          mp_sum_c5d, mp_sum_c6d, mp_sum_rmm, mp_sum_cmm, mp_sum_r5d
#if defined(__CUDA)
        MODULE PROCEDURE  mp_sum_i1_gpu, mp_sum_iv_gpu, mp_sum_im_gpu, mp_sum_it_gpu, &
          mp_sum_r1_gpu, mp_sum_rv_gpu, mp_sum_rm_gpu, mp_sum_rt_gpu, mp_sum_r4d_gpu, &
          mp_sum_c1_gpu, mp_sum_cv_gpu, mp_sum_cm_gpu, mp_sum_ct_gpu, mp_sum_c4d_gpu, &
          mp_sum_c5d_gpu, mp_sum_c6d_gpu, mp_sum_rmm_gpu, mp_sum_cmm_gpu, mp_sum_r5d_gpu
#endif
      END INTERFACE

      INTERFACE mp_root_sum
        MODULE PROCEDURE mp_root_sum_rm, mp_root_sum_cm
#if defined(__CUDA)
        MODULE PROCEDURE mp_root_sum_rm_gpu, mp_root_sum_cm_gpu
#endif
      END INTERFACE

      INTERFACE mp_get
        MODULE PROCEDURE mp_get_r1, mp_get_rv, mp_get_cv, mp_get_i1, mp_get_iv, &
          mp_get_rm, mp_get_cm
#if defined(__CUDA)
        MODULE PROCEDURE mp_get_r1_gpu, mp_get_rv_gpu, mp_get_cv_gpu, mp_get_i1_gpu, mp_get_iv_gpu, &
          mp_get_rm_gpu, mp_get_cm_gpu
#endif
      END INTERFACE

      INTERFACE mp_put
        MODULE PROCEDURE mp_put_rv, mp_put_cv, mp_put_i1, mp_put_iv, &
          mp_put_rm
#if defined(__CUDA)
        MODULE PROCEDURE mp_put_rv_gpu, mp_put_cv_gpu, mp_put_i1_gpu, mp_put_iv_gpu, &
          mp_put_rm_gpu
#endif
      END INTERFACE

      INTERFACE mp_max
        MODULE PROCEDURE mp_max_i, mp_max_r, mp_max_rv, mp_max_iv
#if defined(__CUDA)
        MODULE PROCEDURE mp_max_i_gpu, mp_max_r_gpu, mp_max_rv_gpu, mp_max_iv_gpu
#endif
      END INTERFACE

      INTERFACE mp_min
        MODULE PROCEDURE mp_min_i, mp_min_r, mp_min_rv, mp_min_iv
#if defined(__CUDA)
        MODULE PROCEDURE mp_min_i_gpu, mp_min_r_gpu, mp_min_rv_gpu, mp_min_iv_gpu
#endif
      END INTERFACE

      INTERFACE mp_gather
        MODULE PROCEDURE mp_gather_i1, mp_gather_iv, mp_gatherv_rv, mp_gatherv_iv, &
                         mp_gatherv_rm, mp_gatherv_im, mp_gatherv_cv, &
                         mp_gatherv_inplace_cplx_array
#if defined(__CUDA)
        MODULE PROCEDURE mp_gather_i1_gpu, mp_gather_iv_gpu, mp_gatherv_rv_gpu, mp_gatherv_iv_gpu, &
          mp_gatherv_rm_gpu, mp_gatherv_im_gpu, mp_gatherv_cv_gpu, mp_gatherv_inplace_cplx_array_gpu
#endif
      END INTERFACE

      INTERFACE mp_allgather
        MODULE PROCEDURE mp_allgatherv_inplace_cplx_array
#if defined(__CUDA)
        MODULE PROCEDURE mp_allgatherv_inplace_cplx_array_gpu
#endif
      END INTERFACE

      INTERFACE mp_alltoall
        MODULE PROCEDURE mp_alltoall_c3d, mp_alltoall_i3d
#if defined(__CUDA)
        MODULE PROCEDURE mp_alltoall_c3d_gpu, mp_alltoall_i3d_gpu
#endif
      END INTERFACE

      INTERFACE mp_circular_shift_left
        MODULE PROCEDURE mp_circular_shift_left_i0, &
          mp_circular_shift_left_i1, &
          mp_circular_shift_left_i2, &
          mp_circular_shift_left_r2d, &
          mp_circular_shift_left_c2d
#if defined(__CUDA)
        MODULE PROCEDURE mp_circular_shift_left_i0_gpu, &
          mp_circular_shift_left_i1_gpu, &
          mp_circular_shift_left_i2_gpu, &
          mp_circular_shift_left_r2d_gpu, &
          mp_circular_shift_left_c2d_gpu
#endif
      END INTERFACE

      INTERFACE mp_type_create_column_section
        MODULE PROCEDURE mp_type_create_cplx_column_section
#if defined(__CUDA)
        MODULE PROCEDURE mp_type_create_cplx_column_section_gpu
#endif
      END INTERFACE

!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
!..mp_gather_i1
      SUBROUTINE mp_gather_i1(mydata, alldata, root, gid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: mydata, root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER, INTENT(OUT) :: alldata(:)
        INTEGER :: ierr


#if defined (__MPI)
        group = gid
        CALL MPI_GATHER(mydata, 1, MPI_INTEGER, alldata, 1, MPI_INTEGER, root, group, IERR)
        IF (ierr/=0) CALL mp_stop( 8013 )
#else
        alldata(1) = mydata
#endif
        RETURN
      END SUBROUTINE mp_gather_i1

!------------------------------------------------------------------------------!
!..mp_gather_iv
!..Carlo Cavazzoni
      SUBROUTINE mp_gather_iv(mydata, alldata, root, gid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: mydata(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER, INTENT(OUT) :: alldata(:,:)
        INTEGER :: msglen, ierr


#if defined (__MPI)
        msglen = SIZE(mydata)
        IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop( 8014 )
        group = gid
        CALL MPI_GATHER(mydata, msglen, MPI_INTEGER, alldata, msglen, MPI_INTEGER, root, group, IERR)
        IF (ierr/=0) CALL mp_stop( 8014 )
#else
        msglen = SIZE(mydata)
        IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop( 8014 )
        alldata(:,1) = mydata(:)
#endif
        RETURN
      END SUBROUTINE mp_gather_iv

!
!------------------------------------------------------------------------------!
!..mp_start
      SUBROUTINE mp_start(numtask, taskid, group)

! ...
        IMPLICIT NONE
        INTEGER, INTENT (OUT) :: numtask, taskid
        INTEGER, INTENT (IN)  :: group
        INTEGER :: ierr
! ...
        ierr = 0
        numtask = 1
        taskid = 0

#  if defined(__MPI)
        IF (ierr/=0) CALL mp_stop( 8004 )
        CALL mpi_comm_rank(group,taskid,ierr)
        IF (ierr/=0) CALL mp_stop( 8005 )
        CALL mpi_comm_size(group,numtask,ierr)
        IF (ierr/=0) CALL mp_stop( 8006 )
! ...
        CALL allocate_buffers()
#if defined(__CUDA)
        CALL allocate_buffers_gpu()
#endif
#endif
        RETURN
      END SUBROUTINE mp_start
!
!------------------------------------------------------------------------------!
!..mp_abort

      SUBROUTINE mp_abort(errorcode,gid)
        IMPLICIT NONE
        INTEGER :: ierr
        INTEGER, INTENT(IN):: errorcode, gid
#if defined(__MPI)
        CALL deallocate_buffers()
#if defined(__CUDA)
        CALL deallocate_buffers_gpu()
#endif
        CALL mpi_abort(gid, errorcode, ierr)
#endif
      END SUBROUTINE mp_abort
!
!------------------------------------------------------------------------------!
!..mp_end

      SUBROUTINE mp_end(groupid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: groupid
        INTEGER :: ierr, taskid

        ierr = 0
        taskid = 0

#if defined(__MPI)
        CALL mpi_comm_rank( groupid, taskid, ierr)
        CALL deallocate_buffers()
#if defined(__CUDA)
        CALL deallocate_buffers_gpu()
#endif
#endif
        RETURN
      END SUBROUTINE mp_end

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
         IF (ierr/=0) CALL mp_stop( 8007 )
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
         IF (ierr/=0) CALL mp_stop( 8008 )
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
        IF (ierr/=0) CALL mp_stop( 8009 )
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
        IF (ierr/=0) CALL mp_stop( 8010 )
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
        IF (ierr/=0) CALL mp_stop( 8011 )
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
            IF (ierr/=0) CALL mp_stop( 8012 )
         END IF
#endif
         RETURN
      END SUBROUTINE mp_comm_free

!------------------------------------------------------------------------------!
!..mp_bcast

      SUBROUTINE mp_bcast_i1(msg,source,gid)
        IMPLICIT NONE
        INTEGER :: msg
        INTEGER :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen

#if defined(__MPI)
        msglen = 1
        group = gid
        CALL bcast_integer( msg, msglen, source, group )
#endif
      END SUBROUTINE mp_bcast_i1
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_iv(msg,source,gid)
        IMPLICIT NONE
        INTEGER :: msg(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_integer( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_iv
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_im( msg, source, gid )
        IMPLICIT NONE
        INTEGER :: msg(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_integer( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_im
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_it( msg, source, gid )
        IMPLICIT NONE
        INTEGER :: msg(:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_integer( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_it
!
!------------------------------------------------------------------------------!
!
! Samuel Ponce
!
      SUBROUTINE mp_bcast_i4d(msg, source, gid)
        IMPLICIT NONE
        INTEGER :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_integer( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_i4d
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_r1( msg, source, gid )
        IMPLICIT NONE
        REAL (DP) :: msg
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL bcast_real( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_r1
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rv(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_rv
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rm(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, msglen, source, gid )
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
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, msglen, source, gid )
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
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, msglen, source, gid )
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
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_r5d

!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_c1(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL bcast_real( msg, 2 * msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_c1
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_cv(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, 2 * msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_cv
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_cm(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, 2 * msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_cm
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_ct(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, 2 * msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_ct

!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_c4d(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, 2 * msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_c4d

      SUBROUTINE mp_bcast_c5d(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_real( msg, 2 * msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_c5d

!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_bcast_l(msg,source,gid)
        IMPLICIT NONE
        LOGICAL :: msg
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL bcast_logical( msg, msglen, source, gid )
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
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_logical( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_lv

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_lm(msg,source,gid)
        IMPLICIT NONE
        LOGICAL :: msg(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL bcast_logical( msg, msglen, source, gid )
#endif
      END SUBROUTINE mp_bcast_lm


!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_z(msg,source,gid)
        IMPLICIT NONE
        CHARACTER (len=*) :: msg
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr, i
        INTEGER, ALLOCATABLE :: imsg(:)
#if defined(__MPI)
        ierr = 0
        msglen = len(msg)
        group = gid
        ALLOCATE (imsg(1:msglen), STAT=ierr)
        IF (ierr/=0) CALL mp_stop( 8015 )
        DO i = 1, msglen
          imsg(i) = ichar(msg(i:i))
        END DO
        CALL bcast_integer( imsg, msglen, source, group )
        DO i = 1, msglen
          msg(i:i) = char(imsg(i))
        END DO
        DEALLOCATE (imsg, STAT=ierr)
        IF (ierr/=0) CALL mp_stop( 8016 )
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
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, m1, m2, ierr, i, j
        INTEGER, ALLOCATABLE :: imsg(:,:)
#if defined(__MPI)
        ierr = 0
        m1 = LEN(msg)
        m2 = SIZE(msg)
        msglen = LEN(msg)*SIZE(msg)
        group = gid
        ALLOCATE (imsg(1:m1,1:m2), STAT=ierr)
        IF (ierr/=0) CALL mp_stop( 8017 )
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
        IF (ierr/=0) CALL mp_stop( 8018 )
#endif
      END SUBROUTINE mp_bcast_zv
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_i1(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        INTEGER :: msg_dest
        INTEGER, INTENT(IN) :: msg_sour
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen = 1

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(dest .NE. sour) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen=1
             CALL MPI_SEND( msg_sour, msglen, MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8019 )
           ELSE IF(mpime .EQ. dest) THEN
             msglen=1
             CALL MPI_RECV( msg_dest, msglen, MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8020 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8021 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest = msg_sour
          msglen = 1
        END IF

#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8022 )
#endif


        RETURN
      END SUBROUTINE mp_get_i1

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_iv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        INTEGER :: msg_dest(:)
        INTEGER, INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen = SIZE(msg_sour)
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8023 )
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8024 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8025 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8026 )
#endif
        RETURN
      END SUBROUTINE mp_get_iv

!------------------------------------------------------------------------------!

      SUBROUTINE mp_get_r1(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        REAL (DP)             :: msg_dest
        REAL (DP), INTENT(IN) :: msg_sour
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen = 1
             CALL MPI_SEND( msg_sour, msglen, MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8027 )
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, 1, MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8028 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8029 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest = msg_sour
          msglen = 1
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8030 )
#endif
        RETURN
      END SUBROUTINE mp_get_r1

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_rv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        REAL (DP)             :: msg_dest(:)
        REAL (DP), INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen = SIZE(msg_sour)
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8027 )
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8028 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8029 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8030 )
#endif
        RETURN
      END SUBROUTINE mp_get_rv

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_rm(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        REAL (DP) :: msg_dest(:,:)
        REAL (DP), INTENT(IN) :: msg_sour(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8031 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8032 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8033 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour,1), 1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE( msg_sour )
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8034 )
#endif
        RETURN
      END SUBROUTINE mp_get_rm


!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_cv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        COMPLEX (DP)             :: msg_dest(:)
        COMPLEX (DP), INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF( dest .NE. sour ) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8035 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8036 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8037 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8038 )
#endif
        RETURN
      END SUBROUTINE mp_get_cv



!------------------------------------------------------------------------------!
!
! Marco Govoni
!
      SUBROUTINE mp_get_cm(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        COMPLEX (DP)              :: msg_dest(:,:)
        COMPLEX (DP), INTENT(IN)  :: msg_sour(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8031 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8032 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8033 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour,1), 1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE( msg_sour )
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8034 )
#endif
        RETURN
      END SUBROUTINE mp_get_cm
!------------------------------------------------------------------------------!
!
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_put_i1(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        INTEGER :: msg_dest
        INTEGER, INTENT(IN) :: msg_sour
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(dest .NE. sour) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, 1, MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8039 )
             msglen = 1
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, 1, MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8040 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8041 )
             msglen = 1
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest = msg_sour
          msglen = 1
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8042 )
#endif
        RETURN
      END SUBROUTINE mp_put_i1

!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_iv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        INTEGER             :: msg_dest(:)
        INTEGER, INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if defined(__MPI)
        group = gid
#endif
        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8043 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8044 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8045 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8046 )
#endif
        RETURN
      END SUBROUTINE mp_put_iv

!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_rv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        REAL (DP)             :: msg_dest(:)
        REAL (DP), INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if defined(__MPI)
        group = gid
#endif
        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8047 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8048 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8049 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8050 )
#endif
        RETURN
      END SUBROUTINE mp_put_rv

!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_rm(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        REAL (DP)             :: msg_dest(:,:)
        REAL (DP), INTENT(IN) :: msg_sour(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if defined(__MPI)
        group = gid
#endif
        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8051 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8052 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8053 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour,1),1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8054 )
#endif
        RETURN
      END SUBROUTINE mp_put_rm


!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_cv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        COMPLEX (DP)             :: msg_dest(:)
        COMPLEX (DP), INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if defined(__MPI)
        group = gid
#endif
        ! processors not taking part in the communication have 0 length message

        msglen = 0

        IF( dest .NE. sour ) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 8055 )
             msglen = SIZE(msg_sour)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 8056 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 8057 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 8058 )
#endif
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
        INTEGER :: ierr
        WRITE( stdout, fmt='( "*** error in Message Passing (mp) module ***")' )
        WRITE( stdout, fmt='( "*** error code: ",I5)' ) code
#if defined(__MPI)
        ! abort with extreme prejudice across the entire MPI set of tasks
        CALL mpi_abort(MPI_COMM_WORLD,code,ierr)
#endif
        STOP
      END SUBROUTINE mp_stop
!------------------------------------------------------------------------------!
!
!..mp_sum
      SUBROUTINE mp_sum_i1(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_i1
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_iv
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_im(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_im
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_it(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_it

!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_r1(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_r1

!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_rv
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_rm(msg, gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_rm


      SUBROUTINE mp_root_sum_rm( msg, res, root, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (IN)  :: msg(:,:)
        REAL (DP), INTENT (OUT) :: res(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER,   INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen, ierr, taskid

        msglen = size(msg)

        CALL mpi_comm_rank( gid, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 8059 )
        !
        IF( taskid == root ) THEN
           IF( msglen > size(res) ) CALL mp_stop( 8060 )
        END IF

        CALL reduce_base_real_to( msglen, msg, res, gid, root )

#else

        res = msg

#endif

      END SUBROUTINE mp_root_sum_rm


      SUBROUTINE mp_root_sum_cm( msg, res, root, gid )
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN)  :: msg(:,:)
        COMPLEX (DP), INTENT (OUT) :: res(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER,  INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen, ierr, taskid

        msglen = size(msg)

        CALL mpi_comm_rank( gid, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 8061 )

        IF( taskid == root ) THEN
           IF( msglen > size(res) ) CALL mp_stop( 8062 )
        END IF

        CALL reduce_base_real_to( 2 * msglen, msg, res, gid, root )

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
        INTEGER, INTENT (IN) :: root
        INTEGER, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
        INTEGER :: taskid, ierr

        msglen = size(msg)

#if defined(__MPI)

        group = gid
        !
        CALL mpi_comm_rank( group, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 8063 )

        IF( taskid == root ) THEN
           IF( msglen > size(res) ) CALL mp_stop( 8064 )
        END IF
        !
        CALL reduce_base_real_to( msglen, msg, res, group, root )
        !

#else
        res = msg
#endif

      END SUBROUTINE mp_sum_rmm


!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_rt( msg, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
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
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_r4d



!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_c1(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_c1
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_cv(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_cv
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_cm(msg, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_cm
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_cmm(msg, res, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN) :: msg(:,:)
        COMPLEX (DP), INTENT (OUT) :: res(:,:)
        INTEGER, INTENT (IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real_to( 2 * msglen, msg, res, gid, -1 )
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
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = SIZE(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
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
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_c4d
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c5d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_c5d

!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_r5d(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_r5d



!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c6d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL reduce_base_real( 2 * msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_sum_c6d



!------------------------------------------------------------------------------!
      SUBROUTINE mp_max_i(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL parallel_max_integer( msglen, msg, gid, -1 )
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
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL parallel_max_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_max_iv
!
!----------------------------------------------------------------------

      SUBROUTINE mp_max_r(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL parallel_max_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_max_r
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_max_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL parallel_max_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_max_rv
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_i(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL parallel_min_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_min_i
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = SIZE(msg)
        CALL parallel_min_integer( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_min_iv
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_r(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = 1
        CALL parallel_min_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_min_r
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
#if defined(__MPI)
        INTEGER :: msglen
        msglen = size(msg)
        CALL parallel_min_real( msglen, msg, gid, -1 )
#endif
      END SUBROUTINE mp_min_rv

!------------------------------------------------------------------------------!

      SUBROUTINE mp_barrier(gid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: gid
        INTEGER :: ierr
#if defined(__MPI)
        CALL MPI_BARRIER(gid,IERR)
        IF (ierr/=0) CALL mp_stop( 8066 )
#endif
      END SUBROUTINE mp_barrier

!------------------------------------------------------------------------------!
!.. Carlo Cavazzoni
!..mp_rank
      FUNCTION mp_rank( comm )
        IMPLICIT NONE
        INTEGER :: mp_rank
        INTEGER, INTENT(IN) :: comm
        INTEGER :: ierr, taskid

        ierr = 0
        taskid = 0
#if defined(__MPI)
        CALL mpi_comm_rank(comm,taskid,ierr)
        IF (ierr/=0) CALL mp_stop( 8067 )
#endif
        mp_rank = taskid
      END FUNCTION mp_rank

!------------------------------------------------------------------------------!
!.. Carlo Cavazzoni
!..mp_size
      FUNCTION mp_size( comm )
        IMPLICIT NONE
        INTEGER :: mp_size
        INTEGER, INTENT(IN) :: comm
        INTEGER :: ierr, numtask

        ierr = 0
        numtask = 1
#if defined(__MPI)
        CALL mpi_comm_size(comm,numtask,ierr)
        IF (ierr/=0) CALL mp_stop( 8068 )
#endif
        mp_size = numtask
      END FUNCTION mp_size

      SUBROUTINE mp_report
        INTEGER :: i
        WRITE( stdout, *)
#if defined(__MPI)
#  if defined (__MP_STAT)
        WRITE( stdout, 20 )
20      FORMAT(3X,'please use an MPI profiler to analyze communications ')
#  endif
#else
        WRITE( stdout, *)
#endif
        RETURN
      END SUBROUTINE mp_report


!------------------------------------------------------------------------------!
!..mp_gatherv_rv
!..Carlo Cavazzoni

      SUBROUTINE mp_gatherv_rv( mydata, alldata, recvcount, displs, root, gid)
        IMPLICIT NONE
        REAL(DP) :: mydata(:)
        REAL(DP) :: alldata(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid

#if defined (__MPI)
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 8069 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 8070 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 8072 )
        END IF
        IF ( SIZE( mydata ) < recvcount( myid + 1 ) ) CALL mp_stop( 8073 )
        !
        CALL MPI_GATHERV( mydata, recvcount( myid + 1 ), MPI_DOUBLE_PRECISION, &
                         alldata, recvcount, displs, MPI_DOUBLE_PRECISION, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 8074 )
#else
        IF ( SIZE( alldata ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
#endif
        RETURN
      END SUBROUTINE mp_gatherv_rv

!------------------------------------------------------------------------------!
!..mp_gatherv_cv
!..Carlo Cavazzoni

      SUBROUTINE mp_gatherv_cv( mydata, alldata, recvcount, displs, root, gid)
        IMPLICIT NONE
        COMPLEX(DP) :: mydata(:)
        COMPLEX(DP) :: alldata(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid

#if defined (__MPI)
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 8069 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 8070 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 8072 )
        END IF
        IF ( SIZE( mydata ) < recvcount( myid + 1 ) ) CALL mp_stop( 8073 )
        !
        CALL MPI_GATHERV( mydata, recvcount( myid + 1 ), MPI_DOUBLE_COMPLEX, &
                         alldata, recvcount, displs, MPI_DOUBLE_COMPLEX, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 8074 )
#else
        IF ( SIZE( alldata ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
#endif
        RETURN
      END SUBROUTINE mp_gatherv_cv

!------------------------------------------------------------------------------!
!..mp_gatherv_rv
!..Carlo Cavazzoni

      SUBROUTINE mp_gatherv_iv( mydata, alldata, recvcount, displs, root, gid)
        IMPLICIT NONE
        INTEGER :: mydata(:)
        INTEGER :: alldata(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid

#if defined (__MPI)
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 8069 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 8070 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 8072 )
        END IF
        IF ( SIZE( mydata ) < recvcount( myid + 1 ) ) CALL mp_stop( 8073 )
        !
        CALL MPI_GATHERV( mydata, recvcount( myid + 1 ), MPI_INTEGER, &
                         alldata, recvcount, displs, MPI_INTEGER, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 8074 )
#else
        IF ( SIZE( alldata ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
#endif
        RETURN
      END SUBROUTINE mp_gatherv_iv


!------------------------------------------------------------------------------!
!..mp_gatherv_rm
!..Carlo Cavazzoni

      SUBROUTINE mp_gatherv_rm( mydata, alldata, recvcount, displs, root, gid)
        IMPLICIT NONE
        REAL(DP) :: mydata(:,:)  ! Warning first dimension is supposed constant!
        REAL(DP) :: alldata(:,:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid, nsiz
        INTEGER, ALLOCATABLE :: nrecv(:), ndisp(:)


#if defined (__MPI)
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 8069 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 8070 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata, 2 ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 8072 )
           IF ( SIZE( alldata, 1 ) /= SIZE( mydata, 1 ) ) CALL mp_stop( 8072 )
        END IF
        IF ( SIZE( mydata, 2 ) < recvcount( myid + 1 ) ) CALL mp_stop( 8073 )
        !
        ALLOCATE( nrecv( npe ), ndisp( npe ) )
        !
        nrecv( 1:npe ) = recvcount( 1:npe ) * SIZE( mydata, 1 )
        ndisp( 1:npe ) = displs( 1:npe ) * SIZE( mydata, 1 )
        !
        CALL MPI_GATHERV( mydata, nrecv( myid + 1 ), MPI_DOUBLE_PRECISION, &
                         alldata, nrecv, ndisp, MPI_DOUBLE_PRECISION, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 8074 )
        !
        DEALLOCATE( nrecv, ndisp )
        !
#else
        IF ( SIZE( alldata, 1 ) /= SIZE( mydata, 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( alldata, 2 ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata, 2  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( :, 1:recvcount( 1 ) ) = mydata( :, 1:recvcount( 1 ) )
#endif
        RETURN
      END SUBROUTINE mp_gatherv_rm

!------------------------------------------------------------------------------!
!..mp_gatherv_im
!..Carlo Cavazzoni

      SUBROUTINE mp_gatherv_im( mydata, alldata, recvcount, displs, root, gid)
        IMPLICIT NONE
        INTEGER :: mydata(:,:)  ! Warning first dimension is supposed constant!
        INTEGER :: alldata(:,:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid, nsiz
        INTEGER, ALLOCATABLE :: nrecv(:), ndisp(:)


#if defined (__MPI)
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 8069 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 8070 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata, 2 ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 8072 )
           IF ( SIZE( alldata, 1 ) /= SIZE( mydata, 1 ) ) CALL mp_stop( 8072 )
        END IF
        IF ( SIZE( mydata, 2 ) < recvcount( myid + 1 ) ) CALL mp_stop( 8073 )
        !
        ALLOCATE( nrecv( npe ), ndisp( npe ) )
        !
        nrecv( 1:npe ) = recvcount( 1:npe ) * SIZE( mydata, 1 )
        ndisp( 1:npe ) = displs( 1:npe ) * SIZE( mydata, 1 )
        !
        CALL MPI_GATHERV( mydata, nrecv( myid + 1 ), MPI_INTEGER, &
                         alldata, nrecv, ndisp, MPI_INTEGER, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 8074 )
        !
        DEALLOCATE( nrecv, ndisp )
        !
#else
        IF ( SIZE( alldata, 1 ) /= SIZE( mydata, 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( alldata, 2 ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata, 2  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( :, 1:recvcount( 1 ) ) = mydata( :, 1:recvcount( 1 ) )
#endif
        RETURN
      END SUBROUTINE mp_gatherv_im


!------------------------------------------------------------------------------!
!..mp_gatherv_inplace_cplx_array
!..Ye Luo

      SUBROUTINE mp_gatherv_inplace_cplx_array(alldata, my_column_type, recvcount, displs, root, gid)
        IMPLICIT NONE
        COMPLEX(DP) :: alldata(:,:)
        INTEGER, INTENT(IN) :: my_column_type
        INTEGER, INTENT(IN) :: recvcount(:), displs(:)
        INTEGER, INTENT(IN) :: root, gid
        INTEGER :: ierr, npe, myid

#if defined (__MPI)
        CALL mpi_comm_size( gid, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 8069 )
        CALL mpi_comm_rank( gid, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 8070 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
        !
        IF (myid==root) THEN
           CALL MPI_GATHERV( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                             alldata, recvcount, displs, my_column_type, root, gid, ierr )
        ELSE
           CALL MPI_GATHERV( alldata(1,displs(myid+1)+1), recvcount(myid+1), my_column_type, &
                             MPI_IN_PLACE, recvcount, displs, MPI_DATATYPE_NULL, root, gid, ierr )
        ENDIF
        IF (ierr/=0) CALL mp_stop( 8074 )
#endif
        RETURN
      END SUBROUTINE mp_gatherv_inplace_cplx_array

!------------------------------------------------------------------------------!
!..mp_allgatherv_inplace_cplx_array
!..Ye Luo

      SUBROUTINE mp_allgatherv_inplace_cplx_array(alldata, my_element_type, recvcount, displs, gid)
        IMPLICIT NONE
        COMPLEX(DP) :: alldata(:,:)
        INTEGER, INTENT(IN) :: my_element_type
        INTEGER, INTENT(IN) :: recvcount(:), displs(:)
        INTEGER, INTENT(IN) :: gid
        INTEGER :: ierr, npe, myid

#if defined (__MPI)
        CALL mpi_comm_size( gid, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 8069 )
        CALL mpi_comm_rank( gid, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 8070 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
        !
        CALL MPI_ALLGATHERV( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                             alldata, recvcount, displs, my_element_type, gid, ierr )
        IF (ierr/=0) CALL mp_stop( 8074 )
#endif
        RETURN
      END SUBROUTINE mp_allgatherv_inplace_cplx_array

!------------------------------------------------------------------------------!

      SUBROUTINE mp_set_displs( recvcount, displs, ntot, nproc )
        !  Given the number of elements on each processor (recvcount), this subroutine
        !  sets the correct offsets (displs) to collect them on a single
        !  array with contiguous elemets
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: recvcount(:) ! number of elements on each processor
        INTEGER, INTENT(OUT) :: displs(:)   ! offsets/displacements
        INTEGER, INTENT(OUT) :: ntot
        INTEGER, INTENT(IN) :: nproc
        INTEGER :: i

        displs( 1 ) = 0
        !
#if defined (__MPI)
        IF( nproc < 1 ) CALL mp_stop( 8090 )
        DO i = 2, nproc
           displs( i ) = displs( i - 1 ) + recvcount( i - 1 )
        END DO
        ntot = displs( nproc ) + recvcount( nproc )
#else
        ntot = recvcount( 1 )
#endif
        RETURN
      END SUBROUTINE mp_set_displs

!------------------------------------------------------------------------------!


SUBROUTINE mp_alltoall_c3d( sndbuf, rcvbuf, gid )
   IMPLICIT NONE
   COMPLEX(DP) :: sndbuf( :, :, : )
   COMPLEX(DP) :: rcvbuf( :, :, : )
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe

#if defined (__MPI)

   group = gid

   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8069 )

   IF ( SIZE( sndbuf, 3 ) < npe ) CALL mp_stop( 8069 )
   IF ( SIZE( rcvbuf, 3 ) < npe ) CALL mp_stop( 8069 )

   nsiz = SIZE( sndbuf, 1 ) * SIZE( sndbuf, 2 )

   CALL MPI_ALLTOALL( sndbuf, nsiz, MPI_DOUBLE_COMPLEX, &
                      rcvbuf, nsiz, MPI_DOUBLE_COMPLEX, group, ierr )

   IF (ierr/=0) CALL mp_stop( 8074 )

#else

   rcvbuf = sndbuf

#endif

   RETURN
END SUBROUTINE mp_alltoall_c3d


!------------------------------------------------------------------------------!

SUBROUTINE mp_alltoall_i3d( sndbuf, rcvbuf, gid )
   IMPLICIT NONE
   INTEGER :: sndbuf( :, :, : )
   INTEGER :: rcvbuf( :, :, : )
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe

#if defined (__MPI)

   group = gid

   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8069 )

   IF ( SIZE( sndbuf, 3 ) < npe ) CALL mp_stop( 8069 )
   IF ( SIZE( rcvbuf, 3 ) < npe ) CALL mp_stop( 8069 )

   nsiz = SIZE( sndbuf, 1 ) * SIZE( sndbuf, 2 )

   CALL MPI_ALLTOALL( sndbuf, nsiz, MPI_INTEGER, &
                      rcvbuf, nsiz, MPI_INTEGER, group, ierr )

   IF (ierr/=0) CALL mp_stop( 8074 )

#else

   rcvbuf = sndbuf

#endif

   RETURN
END SUBROUTINE mp_alltoall_i3d

SUBROUTINE mp_circular_shift_left_i0( buf, itag, gid )
   IMPLICIT NONE
   INTEGER :: buf
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)

   INTEGER :: istatus( mpi_status_size )
   !
   group = gid
   !
   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8100 )
   CALL mpi_comm_rank( group, mype, ierr )
   IF (ierr/=0) CALL mp_stop( 8101 )
   !
   sour = mype + 1
   IF( sour == npe ) sour = 0
   dest = mype - 1
   IF( dest == -1 ) dest = npe - 1
   !
   CALL MPI_Sendrecv_replace( buf, 1, MPI_INTEGER, &
        dest, itag, sour, itag, group, istatus, ierr)
   !
   IF (ierr/=0) CALL mp_stop( 8102 )
   !
#else
   ! do nothing
#endif
   RETURN
END SUBROUTINE mp_circular_shift_left_i0


SUBROUTINE mp_circular_shift_left_i1( buf, itag, gid )
   IMPLICIT NONE
   INTEGER :: buf(:)
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)

   INTEGER :: istatus( mpi_status_size )
   !
   group = gid
   !
   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8100 )
   CALL mpi_comm_rank( group, mype, ierr )
   IF (ierr/=0) CALL mp_stop( 8101 )
   !
   sour = mype + 1
   IF( sour == npe ) sour = 0
   dest = mype - 1
   IF( dest == -1 ) dest = npe - 1
   !
   CALL MPI_Sendrecv_replace( buf, SIZE(buf), MPI_INTEGER, &
        dest, itag, sour, itag, group, istatus, ierr)
   !
   IF (ierr/=0) CALL mp_stop( 8102 )
   !
#else
   ! do nothing
#endif
   RETURN
END SUBROUTINE mp_circular_shift_left_i1


SUBROUTINE mp_circular_shift_left_i2( buf, itag, gid )
   IMPLICIT NONE
   INTEGER :: buf(:,:)
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)

   INTEGER :: istatus( mpi_status_size )
   !
   group = gid
   !
   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8100 )
   CALL mpi_comm_rank( group, mype, ierr )
   IF (ierr/=0) CALL mp_stop( 8101 )
   !
   sour = mype + 1
   IF( sour == npe ) sour = 0
   dest = mype - 1
   IF( dest == -1 ) dest = npe - 1
   !
   CALL MPI_Sendrecv_replace( buf, SIZE(buf), MPI_INTEGER, &
        dest, itag, sour, itag, group, istatus, ierr)
   !
   IF (ierr/=0) CALL mp_stop( 8102 )
   !
#else
   ! do nothing
#endif
   RETURN
END SUBROUTINE mp_circular_shift_left_i2


SUBROUTINE mp_circular_shift_left_r2d( buf, itag, gid )
   IMPLICIT NONE
   REAL(DP) :: buf( :, : )
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)

   INTEGER :: istatus( mpi_status_size )
   !
   group = gid
   !
   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8100 )
   CALL mpi_comm_rank( group, mype, ierr )
   IF (ierr/=0) CALL mp_stop( 8101 )
   !
   sour = mype + 1
   IF( sour == npe ) sour = 0
   dest = mype - 1
   IF( dest == -1 ) dest = npe - 1
   !
   CALL MPI_Sendrecv_replace( buf, SIZE(buf), MPI_DOUBLE_PRECISION, &
        dest, itag, sour, itag, group, istatus, ierr)
   !
   IF (ierr/=0) CALL mp_stop( 8102 )
   !
#else
   ! do nothing
#endif
   RETURN
END SUBROUTINE mp_circular_shift_left_r2d

SUBROUTINE mp_circular_shift_left_c2d( buf, itag, gid )
   IMPLICIT NONE
   COMPLEX(DP) :: buf( :, : )
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)

   INTEGER :: istatus( mpi_status_size )
   !
   group = gid
   !
   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8100 )
   CALL mpi_comm_rank( group, mype, ierr )
   IF (ierr/=0) CALL mp_stop( 8101 )
   !
   sour = mype + 1
   IF( sour == npe ) sour = 0
   dest = mype - 1
   IF( dest == -1 ) dest = npe - 1
   !
   CALL MPI_Sendrecv_replace( buf, SIZE(buf), MPI_DOUBLE_COMPLEX, &
        dest, itag, sour, itag, group, istatus, ierr)
   !
   IF (ierr/=0) CALL mp_stop( 8102 )
   !
#else
   ! do nothing
#endif
   RETURN
END SUBROUTINE mp_circular_shift_left_c2d
!
!
!------------------------------------------------------------------------------!
!..mp_count_nodes
SUBROUTINE mp_count_nodes(num_nodes, color, key, group)
  !
  ! ... This routine counts the number of nodes using
  ! ...  MPI_GET_PROCESSOR_NAME in the group specified by `group`.
  ! ...  It returns colors and keys to be used in MPI_COMM_SPLIT.
  ! ...  When running in parallel, the evaluation of color and key
  ! ...  is done by all processors.
  ! ...
  ! ...
  ! ... input:
  ! ...    group      Communicator used to count nodes.
  !
  ! ... output:
  ! ...    num_nodes  Number of unique nodes in the communicator
  ! ...    color      Integer (positive), same for all processes residing on a node.
  ! ...    key        Integer, unique number identifying each process on the same node.
  ! ...
  IMPLICIT NONE
  INTEGER, INTENT (OUT) :: num_nodes
  INTEGER, INTENT (OUT) :: color
  INTEGER, INTENT (OUT) :: key
  INTEGER, INTENT (IN)  :: group
#if defined (__MPI)
  CHARACTER(len=MPI_MAX_PROCESSOR_NAME) :: hostname
  CHARACTER(len=MPI_MAX_PROCESSOR_NAME), ALLOCATABLE :: host_list(:)
#endif
  
  LOGICAL, ALLOCATABLE   :: found_list(:)
  INTEGER, ALLOCATABLE   :: color_list(:)
  INTEGER, ALLOCATABLE   :: key_list(:)
  !
  INTEGER :: hostname_len, max_hostname_len, numtask, me, ierr
  !
  ! Loops variables
  INTEGER :: i, j, e, s, c, k
  ! ...
  ierr      = 0
  num_nodes = 1
  color     = 1
  key       = 0
  !
#if defined(__MPI)
  !
  CALL MPI_GET_PROCESSOR_NAME(hostname, hostname_len, ierr)
  IF (ierr/=0)  CALL mp_stop( 8103 )

  ! find total number of ranks and my rank in communicator
  CALL MPI_COMM_SIZE(group, numtask, ierr)
  IF (ierr/=0) CALL mp_stop( 8104 )
  !
  CALL MPI_COMM_RANK(group, me, ierr)
  IF (ierr/=0) CALL mp_stop( 8105 )
  !
  ALLOCATE(host_list(0:numtask-1))
  !
  host_list(me) = hostname(1:hostname_len)
  !
  ! Each process broadcast its name to the others
  DO i=0,numtask-1
    CALL MPI_BCAST(host_list(i), MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER,&
                     i, group, ierr)
    IF (ierr/=0) CALL mp_stop( 8106 )
  END DO
  !
  ! Simple algorithm to count unique entries.
  !
  ALLOCATE(found_list(0:numtask-1),color_list(0:numtask-1))
  ALLOCATE(key_list(0:numtask-1))
  found_list(:) = .false.
  color_list(:) = -1
  key_list(:)   = -1
  !
  ! c is the counter for colors
  ! k is the counter for keys
  !
  c = 0
  DO i=0,numtask-1
    ! if node_counter == .true., this element has already been found,
    ! so skip it.
    IF (found_list(i)) CYCLE
    ! else increment color counter and reset key counter
    c = c + 1; k = 0
    color_list(i) = c
    key_list(i)   = k
    !
    DO j=i+1,numtask-1
      !
      IF ( LLE(host_list(i),host_list(j)) .and. &
           LGE(host_list(i),host_list(j))        ) THEN
        ! increment the key, key=0 is the one we are comparing to
        k = k + 1
        ! element should not be already found
        IF ( found_list(j) ) CALL mp_stop( 8107 )
        found_list(j) = .true.
        color_list(j) = c
        key_list(j)   = k
      END IF
    END DO
  END DO
  ! Sanity checks
  IF ( MINVAL(color_list) < 0 ) CALL mp_stop( 8108 )
  IF ( MINVAL(key_list)   < 0 ) CALL mp_stop( 8109 )
  !
  color     = color_list(me)
  key       = key_list(me)
  num_nodes = MAXVAL(color_list)
  DEALLOCATE(host_list,found_list,color_list,key_list)
!
#endif
  RETURN
END SUBROUTINE mp_count_nodes
!
FUNCTION mp_get_comm_null( )
  IMPLICIT NONE
  INTEGER :: mp_get_comm_null
  mp_get_comm_null = MPI_COMM_NULL
END FUNCTION mp_get_comm_null

FUNCTION mp_get_comm_self( )
  IMPLICIT NONE
  INTEGER :: mp_get_comm_self
  mp_get_comm_self = MPI_COMM_SELF
END FUNCTION mp_get_comm_self

SUBROUTINE mp_type_create_cplx_column_section(dummy, start, length, stride, mytype)
  IMPLICIT NONE
  !
  COMPLEX (DP), INTENT(IN) :: dummy
  INTEGER, INTENT(IN) :: start, length, stride
  INTEGER, INTENT(OUT) :: mytype
  !
#if defined(__MPI)
  INTEGER :: ierr
  !
  CALL MPI_TYPE_CREATE_SUBARRAY(1, stride, length, start, MPI_ORDER_FORTRAN,&
                                MPI_DOUBLE_COMPLEX, mytype, ierr)
  IF (ierr/=0) CALL mp_stop( 8081 )
  CALL MPI_Type_commit(mytype, ierr)
  IF (ierr/=0) CALL mp_stop( 8082 )
#else
  mytype = 0;
#endif
  !
  RETURN
END SUBROUTINE mp_type_create_cplx_column_section

SUBROUTINE mp_type_free(mytype)
  IMPLICIT NONE
  INTEGER :: mytype, ierr
  !
#if defined(__MPI)
  CALL MPI_TYPE_FREE(mytype, ierr)
  IF (ierr/=0) CALL mp_stop( 8083 )
#endif
  !
  RETURN
END SUBROUTINE mp_type_free
!------------------------------------------------------------------------------!
!  GPU specific subroutines (Pietro Bonfa')
!------------------------------------------------------------------------------!
! Before hacking on the CUDA part remember that:
!
! 1. all mp_* interface should be blocking with respect to both MPI and CUDA.
!    MPI will only wait for completion on the default stream therefore device
!    synchronization must be enforced.
! 2. Host -> device memory copies of a memory block of 64 KB or less are
!    asynchronous in the sense that they may return before the data is actually
!    available on the GPU. However, the user is still free to change the buffer
!    as soon as those calls return with no ill effects.
!    (https://devtalk.nvidia.com/default/topic/471866/cuda-programming-and-performance/host-device-memory-copies-up-to-64-kb-are-asynchronous/)
! 3. For transfers from device to either pageable or pinned host memory,
!    the function returns only once the copy has completed.
! 4. GPU synchronization is always enforced even if no communication takes place.
!------------------------------------------------------------------------------!
#ifdef __CUDA

!------------------------------------------------------------------------------!
!..mp_bcast

      SUBROUTINE mp_bcast_i1_gpu(msg_d,source,gid)
        IMPLICIT NONE
        INTEGER, DEVICE :: msg_d
        INTEGER :: msg_h
        INTEGER :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        !
#if defined(__MPI)
        msglen = 1
        group = gid
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()  ! This syncs __GPU_MPI case
        CALL bcast_integer_gpu( msg_d, msglen, source, group )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        msg_h = msg_d                   ! This syncs __MPI case
        CALL bcast_integer( msg_h, msglen, source, group )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_i1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_iv_gpu(msg_d,source,gid)
        IMPLICIT NONE
        INTEGER, DEVICE :: msg_d(:)
        INTEGER, ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
        !
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_integer_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_integer( msg_h, msglen, source, gid )
        msg_d = msg_h; DEALLOCATE( msg_h )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_iv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_im_gpu( msg_d, source, gid )
        IMPLICIT NONE
        INTEGER, DEVICE :: msg_d(:,:)
        INTEGER, ALLOCATABLE :: msg_h(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_integer_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_integer( msg_h, msglen, source, gid )
        msg_d = msg_h; DEALLOCATE( msg_h )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_im_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_it_gpu( msg_d, source, gid )
        IMPLICIT NONE
        INTEGER, DEVICE :: msg_d(:,:,:)

        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_integer_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        INTEGER, ALLOCATABLE :: msg_h(:,:,:)
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_integer( msg_h, msglen, source, gid )
        msg_d = msg_h; DEALLOCATE( msg_h )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_it_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_i4d_gpu(msg_d, source, gid)
        IMPLICIT NONE
        INTEGER, DEVICE :: msg_d(:,:,:,:)
        INTEGER, ALLOCATABLE :: msg_h(:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_integer_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_integer( msg_h, msglen, source, gid )
        msg_d = msg_h; DEALLOCATE( msg_h )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_i4d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_r1_gpu( msg_d, source, gid )
        IMPLICIT NONE
        REAL (DP), DEVICE :: msg_d
        REAL (DP) :: msg_h
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        msg_h=msg_d                         ! This syncs __MPI case
        CALL bcast_real( msg_h, msglen, source, gid )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_r1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rv_gpu(msg_d,source,gid)
        IMPLICIT NONE
        REAL (DP), DEVICE :: msg_d(:)
        REAL (DP), ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_rv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rm_gpu(msg_d,source,gid)
        IMPLICIT NONE
        REAL (DP), DEVICE :: msg_d(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        REAL (DP), ALLOCATABLE :: msg_h(:,:)
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_rm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rt_gpu(msg_d,source,gid)
        IMPLICIT NONE
        REAL (DP), DEVICE :: msg_d(:,:,:)
        REAL (DP), ALLOCATABLE :: msg_h(:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_rt_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_r4d_gpu(msg_d, source, gid)
        IMPLICIT NONE
        REAL (DP), DEVICE :: msg_d(:,:,:,:)
        REAL (DP), ALLOCATABLE :: msg_h(:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_r4d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_r5d_gpu(msg_d, source, gid)
        IMPLICIT NONE
        REAL (DP), DEVICE :: msg_d(:,:,:,:,:)
        REAL (DP), ALLOCATABLE :: msg_h(:,:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_r5d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_c1_gpu(msg_d,source,gid)
        IMPLICIT NONE
        COMPLEX (DP), DEVICE :: msg_d
        COMPLEX (DP) :: msg_h
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, 2 * msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        msg_h=msg_d                         ! This syncs __MPI case
        CALL bcast_real( msg_h, 2 * msglen, source, gid )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_c1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_cv_gpu(msg_d,source,gid)
        IMPLICIT NONE
        COMPLEX (DP), DEVICE :: msg_d(:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, 2 * msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, 2 * msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_cv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_cm_gpu(msg_d,source,gid)
        IMPLICIT NONE
        COMPLEX (DP), DEVICE :: msg_d(:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, 2 * msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, 2 * msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_cm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_ct_gpu(msg_d,source,gid)
        IMPLICIT NONE
        COMPLEX (DP), DEVICE :: msg_d(:,:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, 2 * msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, 2 * msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_ct_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_c4d_gpu(msg_d,source,gid)
        IMPLICIT NONE
        COMPLEX (DP), DEVICE :: msg_d(:,:,:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, 2 * msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, 2 * msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_c4d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_c5d_gpu(msg_d,source,gid)
        IMPLICIT NONE
        COMPLEX (DP), DEVICE :: msg_d(:,:,:,:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_real_gpu( msg_d, 2 * msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        ALLOCATE( msg_h, source=msg_d )     ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_real( msg_h, 2 * msglen, source, gid )
        msg_d = msg_h ; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_c5d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_l_gpu(msg_d,source,gid)
        IMPLICIT NONE
        LOGICAL, DEVICE :: msg_d
        LOGICAL         :: msg_h
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_logical_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        msg_h = msg_d                       ! This syncs __MPI case
        CALL bcast_logical( msg_h, msglen, source, gid )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_l_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_lv_gpu(msg_d,source,gid)
        IMPLICIT NONE
        LOGICAL, DEVICE :: msg_d(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_logical_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        LOGICAL, ALLOCATABLE :: msg_h(:)
        ALLOCATE(msg_h, source=msg_d)       ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_logical( msg_h, msglen, source, gid )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_lv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_lm_gpu(msg_d,source,gid)
        IMPLICIT NONE
        LOGICAL, DEVICE :: msg_d(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()      ! This syncs __GPU_MPI case
        CALL bcast_logical_gpu( msg_d, msglen, source, gid )
        RETURN ! Sync done by MPI call (or inside bcast_xxx_gpu)
#else
        LOGICAL, ALLOCATABLE :: msg_h(:,:)
        ALLOCATE(msg_h, source=msg_d)       ! This syncs __MPI case
        msglen = size(msg_h)
        CALL bcast_logical( msg_h, msglen, source, gid )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_bcast_lm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_get_i1_gpu(msg_dest_d, msg_sour_d, mpime, dest, sour, ip, gid)
        INTEGER, DEVICE             :: msg_dest_d
        INTEGER, INTENT(IN), DEVICE :: msg_sour_d
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group, ierr
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: nrcv
        INTEGER :: msglen = 1

#if ! defined(__GPU_MPI)
        ! Call CPU implementation
        INTEGER :: msg_dest_h, msg_sour_h
        !
        msg_dest_h = msg_dest_d; msg_sour_h = msg_sour_d      ! This syncs __MPI case
        CALL mp_get_i1(msg_dest_h, msg_sour_h, mpime, dest, sour, ip, gid)
        msg_dest_d = msg_dest_h
#else

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(dest .NE. sour) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen=1
             CALL MPI_SEND( msg_sour_d, msglen, MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( -8001 )
           ELSE IF(mpime .EQ. dest) THEN
             msglen=1
             CALL MPI_RECV( msg_dest_d, msglen, MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( -8002 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( -8003 )
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest_d = msg_sour_d
          msglen = 1
        END IF

#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( -8004 )
#endif

#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_get_i1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_get_iv_gpu(msg_dest_d, msg_sour_d, mpime, dest, sour, ip, gid)
        INTEGER, DEVICE             :: msg_dest_d(:)
        INTEGER, INTENT(IN), DEVICE :: msg_sour_d(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        INTEGER, ALLOCATABLE :: msg_dest_h(:), msg_sour_h(:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );       ! This syncs __MPI case
        CALL mp_get_iv(msg_dest_h, msg_sour_h, mpime, dest, sour, ip, gid)
        msg_dest_d = msg_dest_h
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen = SIZE(msg_sour_d)
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9001 )
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9002 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9003 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d)) = msg_sour_d(:)
          ierr = cudaMemcpy(msg_dest_d(1) , msg_sour_d(1), SIZE(msg_sour_d), cudaMemcpyDeviceToDevice )
          msglen = SIZE(msg_sour_d)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9004 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_get_iv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_get_r1_gpu(msg_dest_d, msg_sour_d, mpime, dest, sour, ip, gid)
        REAL (DP), DEVICE             :: msg_dest_d
        REAL (DP), INTENT(IN), DEVICE :: msg_sour_d
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
#if ! defined(__GPU_MPI)
        REAL(DP) :: msg_dest_h, msg_sour_h
        !
        msg_dest_h=msg_dest_d; msg_sour_h=msg_sour_d      ! This syncs __MPI case
        CALL mp_get_r1(msg_dest_h, msg_sour_h, mpime, dest, sour, ip, gid)
        msg_dest_d = msg_dest_h
#else
#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen = 1
             CALL MPI_SEND( msg_sour_d, msglen, MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9005 )
           ELSE IF(mpime .EQ. dest) THEN
             msglen = 1
             CALL MPI_RECV( msg_dest_d, msglen, MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9006 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9007 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest_d = msg_sour_d
          msglen = 1
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9008 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_get_r1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_get_rv_gpu(msg_dest_d, msg_sour_d, mpime, dest, sour, ip, gid)
        REAL (DP), DEVICE             :: msg_dest_d(:)
        REAL (DP), INTENT(IN), DEVICE :: msg_sour_d(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        REAL (DP), ALLOCATABLE :: msg_dest_h(:), msg_sour_h(:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );       ! This syncs __MPI case
        CALL mp_get_rv(msg_dest_h, msg_sour_h, mpime, dest, sour, ip, gid)
        msg_dest_d = msg_dest_h
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else
        !
#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             msglen = SIZE(msg_sour_d)
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9009 )
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9010 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9011 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d)) = msg_sour_d(:)
          ierr = cudaMemcpy(msg_dest_d(1) , msg_sour_d(1), SIZE(msg_sour_d), cudaMemcpyDeviceToDevice )
          msglen = SIZE(msg_sour_d)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9012 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_get_rv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_get_rm_gpu(msg_dest_d, msg_sour_d, mpime, dest, sour, ip, gid)
        REAL (DP), DEVICE             :: msg_dest_d(:,:)
        REAL (DP), INTENT(IN), DEVICE :: msg_sour_d(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        REAL (DP), ALLOCATABLE :: msg_dest_h(:,:), msg_sour_h(:,:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );        ! This syncs __MPI case
        CALL mp_get_rm(msg_dest_h, msg_sour_h, mpime, dest, sour, ip, gid)
        msg_dest_d = msg_dest_h
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9013 )
             msglen = SIZE(msg_sour_d)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9014 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9015 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d,1), 1:SIZE(msg_sour_d,2)) = msg_sour_d(:,:)
          ! function cudaMemcpy2D(dst, dpitch, src, spitch, width, height, kdir)
          ierr = cudaMemcpy2D(msg_dest_d, SIZE(msg_dest_d,1),&
                              msg_sour_d, SIZE(msg_sour_d,1),&
                              SIZE(msg_sour_d,1), SIZE(msg_sour_d,2), &
                              cudaMemcpyDeviceToDevice )
          msglen = SIZE( msg_sour_d )
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9016 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_get_rm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_get_cv_gpu(msg_dest_d, msg_sour_d, mpime, dest, sour, ip, gid)
        COMPLEX (DP), DEVICE             :: msg_dest_d(:)
        COMPLEX (DP), INTENT(IN), DEVICE :: msg_sour_d(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        COMPLEX (DP), ALLOCATABLE :: msg_dest_h(:), msg_sour_h(:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );         ! This syncs __MPI case
        CALL mp_get_cv(msg_dest_h, msg_sour_h, mpime, dest, sour, ip, gid)
        msg_dest_d = msg_dest_h;
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else
        !
#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF( dest .NE. sour ) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9017 )
             msglen = SIZE(msg_sour_d)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9018 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9019 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d)) = msg_sour_d(:)
          ierr = cudaMemcpy(msg_dest_d(1) , msg_sour_d(1), SIZE(msg_sour_d), cudaMemcpyDeviceToDevice )
          msglen = SIZE(msg_sour_d)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9020 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_get_cv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_get_cm_gpu(msg_dest_d, msg_sour_d, mpime, dest, sour, ip, gid)
        COMPLEX (DP), INTENT(IN), DEVICE :: msg_sour_d(:,:)
        COMPLEX (DP), DEVICE             :: msg_dest_d(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        COMPLEX (DP), ALLOCATABLE :: msg_dest_h(:,:), msg_sour_h(:,:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );          ! This syncs __MPI case
        CALL mp_get_cm(msg_dest_h, msg_sour_h, mpime, dest, sour, ip, gid)
        msg_dest_d = msg_dest_h;
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else
        !
#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9021 )
             msglen = SIZE(msg_sour_d)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9022 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9023 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d,1), 1:SIZE(msg_sour_d,2)) = msg_sour_d(:,:)
          ierr = cudaMemcpy2D(msg_dest_d, SIZE(msg_dest_d,1),&
                              msg_sour_d, SIZE(msg_sour_d,1),&
                              SIZE(msg_sour_d,1), SIZE(msg_sour_d,2), &
                              cudaMemcpyDeviceToDevice )
          msglen = SIZE( msg_sour_d )
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9024 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_get_cm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_put_i1_gpu(msg_dest_d, msg_sour_d, mpime, sour, dest, ip, gid)
        INTEGER, DEVICE             :: msg_dest_d
        INTEGER, INTENT(IN), DEVICE :: msg_sour_d
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        INTEGER :: msg_dest_h, msg_sour_h
        !
        msg_dest_h=msg_dest_d ; msg_sour_h=msg_sour_d           ! This syncs __MPI case
        CALL mp_put_i1(msg_dest_h, msg_sour_h, mpime, sour, dest, ip, gid)
        msg_dest_d = msg_dest_h
#else

#if defined(__MPI)
        group = gid
#endif

        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(dest .NE. sour) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour_d, 1, MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9025 )
             msglen = 1
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, 1, MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9026 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9027 )
             msglen = 1
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest_d = msg_sour_d
          msglen = 1
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9028 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_put_i1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_put_iv_gpu(msg_dest_d, msg_sour_d, mpime, sour, dest, ip, gid)
        INTEGER, DEVICE             :: msg_dest_d(:)
        INTEGER, INTENT(IN), DEVICE :: msg_sour_d(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        !
#if ! defined(__GPU_MPI)
        INTEGER, ALLOCATABLE :: msg_dest_h(:), msg_sour_h(:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );           ! This syncs __MPI case
        CALL mp_put_iv(msg_dest_h, msg_sour_h, mpime, sour, dest, ip, gid)
        msg_dest_d = msg_dest_h
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else
        !
#if defined(__MPI)
        group = gid
#endif
        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_INTEGER, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9029 )
             msglen = SIZE(msg_sour_d)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_INTEGER, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9030 )
             CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9031 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d)) = msg_sour_d(:)
          ierr = cudaMemcpy(msg_dest_d(1) , msg_sour_d(1), SIZE(msg_sour_d), cudaMemcpyDeviceToDevice )
          msglen = SIZE(msg_sour_d)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9032 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_put_iv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_put_rv_gpu(msg_dest_d, msg_sour_d, mpime, sour, dest, ip, gid)
        REAL (DP), DEVICE             :: msg_dest_d(:)
        REAL (DP), INTENT(IN), DEVICE :: msg_sour_d(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        REAL (DP), ALLOCATABLE :: msg_dest_h(:), msg_sour_h(:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );           ! This syncs __MPI case
        CALL mp_put_rv(msg_dest_h, msg_sour_h, mpime, sour, dest, ip, gid)
        msg_dest_d = msg_dest_h
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else
        !
#if defined(__MPI)
        group = gid
#endif
        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9033 )
             msglen = SIZE(msg_sour_d)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9034 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9035 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d)) = msg_sour_d(:)
          ierr = cudaMemcpy(msg_dest_d(1) , msg_sour_d(1), SIZE(msg_sour_d), cudaMemcpyDeviceToDevice )
          msglen = SIZE(msg_sour_d)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9036 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_put_rv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_put_rm_gpu(msg_dest_d, msg_sour_d, mpime, sour, dest, ip, gid)
        REAL (DP), DEVICE             :: msg_dest_d(:,:)
        REAL (DP), INTENT(IN), DEVICE :: msg_sour_d(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        REAL (DP), ALLOCATABLE :: msg_dest_h(:,:), msg_sour_h(:,:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );           ! This syncs __MPI case
        CALL mp_put_rm(msg_dest_h, msg_sour_h, mpime, sour, dest, ip, gid)
        msg_dest_d = msg_dest_h
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else
        !
#if defined(__MPI)
        group = gid
#endif
        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF(sour .NE. dest) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9037 )
             msglen = SIZE(msg_sour_d)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9038 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9039 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d,1),1:SIZE(msg_sour_d,2)) = msg_sour_d(:,:)
          ierr = cudaMemcpy2D(msg_dest_d, SIZE(msg_dest_d,1),&
                              msg_sour_d, SIZE(msg_sour_d,1),&
                              SIZE(msg_sour_d,1), SIZE(msg_sour_d,2), &
                              cudaMemcpyDeviceToDevice )
          msglen = SIZE(msg_sour_d)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9040 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_put_rm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_put_cv_gpu(msg_dest_d, msg_sour_d, mpime, sour, dest, ip, gid)
        COMPLEX (DP),             DEVICE :: msg_dest_d(:)
        COMPLEX (DP), INTENT(IN), DEVICE :: msg_sour_d(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
#if defined(__MPI)
        INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
        ! 
#if ! defined(__GPU_MPI)
        COMPLEX (DP), ALLOCATABLE :: msg_dest_h(:), msg_sour_h(:)
        !
        ALLOCATE( msg_dest_h, source=msg_dest_d ); ALLOCATE( msg_sour_h, source=msg_sour_d );           ! This syncs __MPI case
        CALL mp_put_cv(msg_dest_h, msg_sour_h, mpime, sour, dest, ip, gid)
        msg_dest_d = msg_dest_h
        DEALLOCATE(msg_dest_h, msg_sour_h)
#else
        !
#if defined(__MPI)
        group = gid
#endif
        ! processors not taking part in the communication have 0 length message

        msglen = 0
        !
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL and __GPU_MPI
        !
        IF( dest .NE. sour ) THEN
#if defined(__MPI)
           IF(mpime .EQ. sour) THEN
             CALL MPI_SEND( msg_sour_d, SIZE(msg_sour_d), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
             IF (ierr/=0) CALL mp_stop( 9041 )
             msglen = SIZE(msg_sour_d)
           ELSE IF(mpime .EQ. dest) THEN
             CALL MPI_RECV( msg_dest_d, SIZE(msg_dest_d), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
             IF (ierr/=0) CALL mp_stop( 9042 )
             CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
             IF (ierr/=0) CALL mp_stop( 9043 )
             msglen = nrcv
           END IF
#endif
        ELSEIF(mpime .EQ. sour)THEN
          !msg_dest_d(1:SIZE(msg_sour_d)) = msg_sour_d(:)
          ierr = cudaMemcpy(msg_dest_d(1) , msg_sour_d(1), SIZE(msg_sour_d), cudaMemcpyDeviceToDevice )
          msglen = SIZE(msg_sour_d)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop( 9044 )
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI and __GPU_MPI
        RETURN
      END SUBROUTINE mp_put_cv_gpu
!
!------------------------------------------------------------------------------!
!
!..mp_sum
      SUBROUTINE mp_sum_i1_gpu(msg_d,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT), DEVICE :: msg_d
        INTEGER, msg_h
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()  ! This syncs __GPU_MPI
        CALL reduce_base_integer_gpu( msglen, msg_d, gid, -1 )
        RETURN ! No need for final syncronization
#else
        !
        msg_h = msg_d                   ! This syncs __MPI case
        CALL reduce_base_integer( msglen, msg_h, gid, -1 )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_i1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_iv_gpu(msg_d,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT), DEVICE :: msg_d(:)
        INTEGER, ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()  ! This syncs __GPU_MPI
        CALL reduce_base_integer_gpu( msglen, msg_d, gid, -1 )
        RETURN ! No need for final syncronization
#else
        ALLOCATE( msg_h, source=msg_d ) ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_integer( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_iv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_im_gpu(msg_d,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT), DEVICE :: msg_d(:,:)
        INTEGER, ALLOCATABLE :: msg_h(:,:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()  ! This syncs __GPU_MPI
        CALL reduce_base_integer_gpu( msglen, msg_d, gid, -1 )
        RETURN ! No need for final syncronization
#else
        ALLOCATE( msg_h, source=msg_d ) ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_integer( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_im_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_it_gpu(msg_d,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT), DEVICE :: msg_d(:,:,:)
        INTEGER, ALLOCATABLE :: msg_h(:,:,:)
        INTEGER, INTENT (IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()  ! This syncs __GPU_MPI
        CALL reduce_base_integer_gpu( msglen, msg_d, gid, -1 )
        RETURN ! No need for final syncronization
#else
        ALLOCATE( msg_h, source=msg_d ) ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_integer( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_it_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_r1_gpu(msg_d,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d
        REAL(DP) :: msg_h
        INTEGER, INTENT (IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()  ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! No need for final syncronization
#else
        msg_h=msg_d                     ! This syncs __MPI case
        CALL reduce_base_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_r1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_rv_gpu(msg_d,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d(:)
        REAL(DP), ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT (IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()  ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! No need for final syncronization
#else
        ALLOCATE( msg_h, source=msg_d ) ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_rv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_rm_gpu(msg_d, gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d(:,:)
        REAL (DP), ALLOCATABLE :: msg_h(:,:)
        INTEGER, INTENT (IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()  ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! No need for final syncronization
#else
        ALLOCATE( msg_h, source=msg_d ) ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_rm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_root_sum_rm_gpu( msg_d, res_d, root, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (IN) , DEVICE :: msg_d(:,:)
        REAL (DP), INTENT (OUT), DEVICE :: res_d(:,:)
        REAL (DP), ALLOCATABLE :: res_h(:,:), msg_h(:,:)
        INTEGER,   INTENT (IN) :: root
        INTEGER,   INTENT (IN) :: gid
        !
        INTEGER :: msglen, ierr, taskid
#if defined(__MPI)
        !
        CALL mpi_comm_rank( gid, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 9045 )
        !
        msglen = size(msg_d)
        IF( taskid == root ) THEN
           IF( msglen > size(res_d) ) CALL mp_stop( 9046 )
        END IF
#if  defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_to_gpu( msglen, msg_d, res_d, gid, root )
        RETURN ! Sync not needed in this case
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        ALLOCATE( res_h(lbound(msg_h,1):ubound(msg_h,1), lbound(msg_h,2):ubound(msg_h,2)));
        CALL reduce_base_real_to( msglen, msg_h, res_h, gid, root )
        res_d = res_h; DEALLOCATE(msg_h, res_h)
#endif

#else
        res_d = msg_d
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_root_sum_rm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_root_sum_cm_gpu( msg_d, res_d, root, gid )
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN) , DEVICE :: msg_d(:,:)
        COMPLEX (DP), INTENT (OUT), DEVICE :: res_d(:,:)
        COMPLEX (DP), ALLOCATABLE          :: res_h(:,:), msg_h(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER,  INTENT (IN) :: gid
        !
        INTEGER :: msglen, ierr, taskid
#if defined(__MPI)
        msglen = size(msg_d)

        CALL mpi_comm_rank( gid, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 9047 )

        IF( taskid == root ) THEN
           IF( msglen > size(res_d) ) CALL mp_stop( 9048 )
        END IF
#if  defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_to_gpu( 2 * msglen, msg_d, res_d, gid, root )
        RETURN ! Sync not needed in this case
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        ALLOCATE( res_h(lbound(msg_h,1):ubound(msg_h,1), lbound(msg_h,2):ubound(msg_h,2)));
        CALL reduce_base_real_to( 2 * msglen, msg_h, res_h, gid, root )
        res_d = res_h; DEALLOCATE(msg_h, res_h)
#endif
#else
        res_d = msg_d
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_root_sum_cm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_rmm_gpu( msg_d, res_d, root, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (IN), DEVICE :: msg_d(:,:)
        REAL (DP), INTENT (OUT),DEVICE :: res_d(:,:)
        REAL (DP), ALLOCATABLE         :: res_h(:,:), msg_h(:,:)
        INTEGER, INTENT (IN) :: root
        INTEGER, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
        INTEGER :: taskid, ierr



#if defined(__MPI)

        msglen = size(msg_d)
        !
        group = gid
        !
        CALL mpi_comm_rank( group, taskid, ierr)
        IF( ierr /= 0 ) CALL mp_stop( 9049 )

        IF( taskid == root ) THEN
           IF( msglen > size(res_d) ) CALL mp_stop( 9050 )
        END IF
        !
#if  defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_to_gpu( msglen, msg_d, res_d, group, root )
        RETURN ! Sync not needed in this case
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        ALLOCATE( res_h(lbound(msg_h,1):ubound(msg_h,1), lbound(msg_h,2):ubound(msg_h,2)));
        CALL reduce_base_real_to( msglen, msg_h, res_h, gid, root )
        res_d = res_h; DEALLOCATE(msg_h, res_h)
#endif
        !
#else
        res_d = msg_d
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_rmm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_rt_gpu( msg_d, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d(:,:,:)
        REAL (DP), ALLOCATABLE :: msg_h(:,:,:)
        INTEGER, INTENT(IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_rt_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_r4d_gpu(msg_d,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d(:,:,:,:)
        REAL (DP), ALLOCATABLE :: msg_h(:,:,:,:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_r4d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_c1_gpu(msg_d,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT), DEVICE :: msg_d
        COMPLEX (DP) :: msg_h
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if  defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( 2 * msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        msg_h=msg_d                               ! This syncs __MPI case
        CALL reduce_base_real( 2 * msglen, msg_h, gid, -1 )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_c1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_cv_gpu(msg_d,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT), DEVICE :: msg_d(:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( 2 * msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( 2 * msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_cv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_cm_gpu(msg_d, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT), DEVICE :: msg_d(:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:)
        INTEGER, INTENT (IN) :: gid
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( 2 * msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( 2 * msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_cm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_cmm_gpu(msg_d, res_d, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN), DEVICE :: msg_d(:,:)
        COMPLEX (DP), INTENT (OUT), DEVICE :: res_d(:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:), res_h(:,:)
        INTEGER, INTENT (IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_to_gpu( 2 * msglen, msg_d, res_h, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        ALLOCATE( res_h(lbound(msg_h,1):ubound(msg_h,1), lbound(msg_h,2):ubound(msg_h,2)));
        CALL reduce_base_real_to( 2 * msglen, msg_h, res_h, gid, -1 )
        res_d = res_h; DEALLOCATE(msg_h, res_h)
#endif
#else
        res_d = msg_d
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_cmm_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_ct_gpu(msg_d,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT), DEVICE :: msg_d(:,:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:,:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)        
#if  defined(__GPU_MPI)
        msglen = SIZE(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( 2 * msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( 2 * msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_ct_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_c4d_gpu(msg_d,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT), DEVICE :: msg_d(:,:,:,:)
        COMPLEX (DP),ALLOCATABLE :: msg_h(:,:,:,:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( 2 * msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( 2 * msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_c4d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_c5d_gpu(msg_d,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT), DEVICE :: msg_d(:,:,:,:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( 2 * msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( 2 * msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_c5d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_r5d_gpu(msg_d,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d(:,:,:,:,:)
        REAL (DP), ALLOCATABLE :: msg_h(:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_r5d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_sum_c6d_gpu(msg_d,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT), DEVICE :: msg_d(:,:,:,:,:,:)
        COMPLEX (DP), ALLOCATABLE :: msg_h(:,:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL reduce_base_real_gpu( 2 * msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL reduce_base_real( 2 * msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_sum_c6d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_max_i_gpu(msg_d,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT), DEVICE :: msg_d
        INTEGER :: msg_h
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if  defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL parallel_max_integer_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        msg_h = msg_d                             ! This syncs __MPI case
        CALL parallel_max_integer( msglen, msg_h, gid, -1 )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_max_i_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_max_iv_gpu(msg_d,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT), DEVICE :: msg_d(:)
        INTEGER, ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()            ! This syncs __GPU_MPI
        CALL parallel_max_integer_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )           ! This syncs __MPI case
        msglen = size(msg_h)
        CALL parallel_max_integer( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_max_iv_gpu
!
!----------------------------------------------------------------------
!
      SUBROUTINE mp_max_r_gpu(msg_d,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d
        REAL (DP) :: msg_h
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL parallel_max_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        msg_h = msg_d                    ! This syncs __MPI case
        CALL parallel_max_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()  ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_max_r_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_max_rv_gpu(msg_d,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d(:)
        REAL (DP), ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL parallel_max_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )  ! This syncs __MPI case
        msglen = size(msg_h)
        CALL parallel_max_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_max_rv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_min_i_gpu(msg_d,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT), DEVICE :: msg_d
        INTEGER  :: msg_h
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if  defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL parallel_min_integer_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        msg_h = msg_d                    ! This syncs __MPI case
        CALL parallel_min_integer( msglen, msg_h, gid, -1 )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_min_i_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_min_iv_gpu(msg_d,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT), DEVICE :: msg_d(:)
        INTEGER, ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)
        msglen = SIZE(msg_d)
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL parallel_min_integer_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )  ! This syncs __MPI case
        msglen = size(msg_h)
        CALL parallel_min_integer( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_min_iv_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_min_r_gpu(msg_d,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d
        REAL (DP) :: msg_h
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
        msglen = 1
#if  defined(__GPU_MPI)
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL parallel_min_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        msg_h = msg_d                    ! This syncs __MPI case
        CALL parallel_min_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h
#endif
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_min_r_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_min_rv_gpu(msg_d,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT), DEVICE :: msg_d(:)
        REAL (DP), ALLOCATABLE :: msg_h(:)
        INTEGER, INTENT(IN) :: gid
        !
        INTEGER :: msglen, ierr
#if defined(__MPI)
#if  defined(__GPU_MPI)   
        msglen = size(msg_d)
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL parallel_min_real_gpu( msglen, msg_d, gid, -1 )
        RETURN ! Sync not needed after MPI call
#else
        ALLOCATE( msg_h, source=msg_d )  ! This syncs __MPI case
        msglen = size(msg_h)
        CALL parallel_min_real( msglen, msg_h, gid, -1 )
        msg_d = msg_h; DEALLOCATE(msg_h)
#endif
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_min_rv_gpu
!
!------------------------------------------------------------------------------!
!..mp_gather

      SUBROUTINE mp_gather_i1_gpu(mydata_d, alldata_d, root, gid)
        IMPLICIT NONE
        INTEGER, DEVICE :: mydata_d
        INTEGER, INTENT(IN) :: gid, root
        INTEGER :: group
        INTEGER, INTENT(OUT), DEVICE :: alldata_d(:)
        INTEGER :: ierr


#if defined (__MPI)
#if ! defined(__GPU_MPI)
        INTEGER :: mydata_h
        INTEGER, ALLOCATABLE :: alldata_h(:)
        ALLOCATE( alldata_h, source=alldata_d ) ! This syncs __MPI
        mydata_h = mydata_d
        CALL mp_gather_i1(mydata_h, alldata_h, root, gid)
        mydata_d = mydata_h; alldata_d = alldata_h
        DEALLOCATE(alldata_h)
#else
        group = gid
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL MPI_GATHER(mydata_d, 1, MPI_INTEGER, alldata_d, 1, MPI_INTEGER, root, group, IERR)
        IF (ierr/=0) CALL mp_stop( 9051 )
        RETURN ! Sync not needed after MPI call
#endif
#else
        !alldata_d(1) = mydata_d
        ierr = cudaMemcpy( alldata_d(1), mydata_d, 1, &
                                            & cudaMemcpyDeviceToDevice )
        IF (ierr/=0) CALL mp_stop( 9052 )
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_gather_i1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_gather_iv_gpu(mydata_d, alldata_d, root, gid)
        IMPLICIT NONE
        INTEGER, DEVICE :: mydata_d(:)
        INTEGER, INTENT(IN) :: gid, root
        INTEGER :: group
        INTEGER, INTENT(OUT), DEVICE :: alldata_d(:,:)
        INTEGER :: msglen, ierr, i


#if defined (__MPI)
#if ! defined(__GPU_MPI)
        INTEGER, ALLOCATABLE :: mydata_h(:)
        INTEGER, ALLOCATABLE :: alldata_h(:,:)
        ALLOCATE( mydata_h, source=mydata_d )     ! This syncs __MPI
        ALLOCATE( alldata_h, source=alldata_d )
        
        CALL mp_gather_iv(mydata_h, alldata_h, root, gid)
        mydata_d = mydata_h; alldata_d = alldata_h
        DEALLOCATE(alldata_h, mydata_h)
#else
        msglen = SIZE(mydata_d)
        IF( msglen .NE. SIZE(alldata_d, 1) ) CALL mp_stop( 9053 )
        group = gid
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL MPI_GATHER(mydata_d, msglen, MPI_INTEGER, alldata_d, msglen, MPI_INTEGER, root, group, IERR)
        IF (ierr/=0) CALL mp_stop( 9054 )
        RETURN ! Sync not needed after MPI call
#endif
#else
        msglen = SIZE(mydata_d)
        IF( msglen .NE. SIZE(alldata_d, 1) ) CALL mp_stop( 9055 )
        !alldata_d(:,1) = mydata_d(:)
        ierr = cudaMemcpy(alldata_d(:,1) , mydata_d(1), msglen, cudaMemcpyDeviceToDevice )
        IF (ierr/=0) CALL mp_stop( 9056 )
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_gather_iv_gpu
!
!------------------------------------------------------------------------------!
!..mp_gatherv_rv
!
      SUBROUTINE mp_gatherv_rv_gpu( mydata_d, alldata_d, recvcount, displs, root, gid)
        IMPLICIT NONE
        REAL(DP), DEVICE :: mydata_d(:)
        REAL(DP), DEVICE :: alldata_d(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid
#if defined (__MPI)
#if ! defined(__GPU_MPI)
        REAL(DP), ALLOCATABLE :: mydata_h(:)
        REAL(DP), ALLOCATABLE :: alldata_h(:)
        
        ALLOCATE(mydata_h, source=mydata_d)     ! This syncs __MPI
        ALLOCATE(alldata_h, source=alldata_d)
        CALL mp_gatherv_rv( mydata_h, alldata_h, recvcount, displs, root, gid)
        alldata_d = alldata_h ; mydata_d = mydata_h
        DEALLOCATE(alldata_h , mydata_h)
#else

        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 9057 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 9058 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 9059 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata_d ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 9060 )
        END IF
        IF ( SIZE( mydata_d ) < recvcount( myid + 1 ) ) CALL mp_stop( 9061 )
        !
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL MPI_GATHERV( mydata_d, recvcount( myid + 1 ), MPI_DOUBLE_PRECISION, &
                         alldata_d, recvcount, displs, MPI_DOUBLE_PRECISION, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 9062 )
        RETURN ! Sync not needed after MPI call
#endif
#else
        IF ( SIZE( alldata_d ) < recvcount( 1 ) ) CALL mp_stop( 9063 )
        IF ( SIZE( mydata_d  ) < recvcount( 1 ) ) CALL mp_stop( 9064 )
        !
        !alldata_d( 1:recvcount( 1 ) ) = mydata_d( 1:recvcount( 1 ) )
        ierr = cudaMemcpy(alldata_d(1) , mydata_d(1), recvcount( 1 ), cudaMemcpyDeviceToDevice )
        IF (ierr/=0) CALL mp_stop( 9065 )
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_gatherv_rv_gpu
!
!------------------------------------------------------------------------------!
!..mp_gatherv_cv
!
      SUBROUTINE mp_gatherv_cv_gpu( mydata_d, alldata_d, recvcount, displs, root, gid)
        IMPLICIT NONE
        COMPLEX(DP), DEVICE :: mydata_d(:)
        COMPLEX(DP), DEVICE :: alldata_d(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid

#if defined (__MPI)
#if ! defined(__GPU_MPI)
        COMPLEX(DP), ALLOCATABLE :: mydata_h(:)
        COMPLEX(DP), ALLOCATABLE :: alldata_h(:)
        
        ALLOCATE(mydata_h, source=mydata_d)     ! This syncs __MPI
        ALLOCATE(alldata_h, source=alldata_d)
        CALL mp_gatherv_cv( mydata_h, alldata_h, recvcount, displs, root, gid)
        alldata_d = alldata_h ; mydata_d = mydata_h
        DEALLOCATE(alldata_h , mydata_h)
#else
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 9066 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 9067 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 9068 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata_d ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 9069 )
        END IF
        IF ( SIZE( mydata_d ) < recvcount( myid + 1 ) ) CALL mp_stop( 9070 )
        !
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL MPI_GATHERV( mydata_d, recvcount( myid + 1 ), MPI_DOUBLE_COMPLEX, &
                         alldata_d, recvcount, displs, MPI_DOUBLE_COMPLEX, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 9071 )
        RETURN ! Sync not needed after MPI call
#endif
#else
        IF ( SIZE( alldata_d ) < recvcount( 1 ) ) CALL mp_stop( 9072 )
        IF ( SIZE( mydata_d  ) < recvcount( 1 ) ) CALL mp_stop( 9073 )
        !
        !alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
        ierr = cudaMemcpy(alldata_d(1) , mydata_d(1), recvcount( 1 ), cudaMemcpyDeviceToDevice )
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_gatherv_cv_gpu
!
!------------------------------------------------------------------------------!
!..mp_gatherv_rv_gpu
!

      SUBROUTINE mp_gatherv_iv_gpu( mydata_d, alldata_d, recvcount, displs, root, gid)
        IMPLICIT NONE
        INTEGER, DEVICE :: mydata_d(:)
        INTEGER, DEVICE :: alldata_d(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid

#if defined (__MPI)
#if ! defined(__GPU_MPI)
        INTEGER, ALLOCATABLE :: mydata_h(:)
        INTEGER, ALLOCATABLE :: alldata_h(:)
        
        ALLOCATE(mydata_h, source=mydata_d)     ! This syncs __MPI
        ALLOCATE(alldata_h, source=alldata_d)
        CALL mp_gatherv_iv( mydata_h, alldata_h, recvcount, displs, root, gid)
        alldata_d = alldata_h ; mydata_d = mydata_h
        DEALLOCATE(alldata_h , mydata_h)
#else
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 9074 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 9075 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 9076 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata_d ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 9077 )
        END IF
        IF ( SIZE( mydata_d ) < recvcount( myid + 1 ) ) CALL mp_stop( 9078 )
        !
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL MPI_GATHERV( mydata_d, recvcount( myid + 1 ), MPI_INTEGER, &
                         alldata_d, recvcount, displs, MPI_INTEGER, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 9079 )
        RETURN ! Sync not needed after MPI call
#endif
#else
        IF ( SIZE( alldata_d ) < recvcount( 1 ) ) CALL mp_stop( 9080 )
        IF ( SIZE( mydata_d  ) < recvcount( 1 ) ) CALL mp_stop( 9081 )
        !
        !alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
        ierr = cudaMemcpy(alldata_d(1) , mydata_d(1), recvcount( 1 ), cudaMemcpyDeviceToDevice )
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_gatherv_iv_gpu
!
!------------------------------------------------------------------------------!
!..mp_gatherv_rm
!

      SUBROUTINE mp_gatherv_rm_gpu( mydata_d, alldata_d, recvcount, displs, root, gid)
        IMPLICIT NONE
        REAL(DP), DEVICE :: mydata_d(:,:)  ! Warning first dimension is supposed constant!
        REAL(DP), DEVICE :: alldata_d(:,:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid, nsiz
        INTEGER, ALLOCATABLE :: nrecv(:), ndisp(:)


#if defined (__MPI)
#if ! defined(__GPU_MPI)
        REAL(DP), ALLOCATABLE :: mydata_h(:,:)
        REAL(DP), ALLOCATABLE :: alldata_h(:,:)
        
        ALLOCATE(mydata_h, source=mydata_d)     ! This syncs __MPI
        ALLOCATE(alldata_h, source=alldata_d)
        CALL mp_gatherv_rm( mydata_h, alldata_h, recvcount, displs, root, gid)
        alldata_d = alldata_h ; mydata_d = mydata_h
        DEALLOCATE(alldata_h , mydata_h)
#else
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 9082 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 9083 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 9084 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata_d, 2 ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 9085 )
           IF ( SIZE( alldata_d, 1 ) /= SIZE( mydata_d, 1 ) ) CALL mp_stop( 9086 )
        END IF
        IF ( SIZE( mydata_d, 2 ) < recvcount( myid + 1 ) ) CALL mp_stop( 9087 )
        !
        ALLOCATE( nrecv( npe ), ndisp( npe ) )
        !
        nrecv( 1:npe ) = recvcount( 1:npe ) * SIZE( mydata_d, 1 )
        ndisp( 1:npe ) = displs( 1:npe ) * SIZE( mydata_d, 1 )
        !
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL MPI_GATHERV( mydata_d, nrecv( myid + 1 ), MPI_DOUBLE_PRECISION, &
                         alldata_d, nrecv, ndisp, MPI_DOUBLE_PRECISION, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 9088 )
        !
        DEALLOCATE( nrecv, ndisp )
        !
        RETURN ! Sync not needed after MPI call
#endif
#else
        IF ( SIZE( alldata_d, 1 ) /= SIZE( mydata_d, 1 ) ) CALL mp_stop( 9089 )
        IF ( SIZE( alldata_d, 2 ) < recvcount( 1 ) ) CALL mp_stop( 9090 )
        IF ( SIZE( mydata_d, 2  ) < recvcount( 1 ) ) CALL mp_stop( 9091 )
        !
        !alldata( :, 1:recvcount( 1 ) ) = mydata( :, 1:recvcount( 1 ) )
        
        ierr = cudaMemcpy2D(alldata_d, SIZE(alldata_d,1),&
                              mydata_d, SIZE(mydata_d,1),&
                              SIZE(mydata_d,1), recvcount( 1 ), &
                              cudaMemcpyDeviceToDevice )

        IF (ierr/=0) CALL mp_stop( 9092 )
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_gatherv_rm_gpu
!
!------------------------------------------------------------------------------!
!..mp_gatherv_im
!
      SUBROUTINE mp_gatherv_im_gpu( mydata_d, alldata_d, recvcount, displs, root, gid)
        IMPLICIT NONE
        INTEGER, DEVICE :: mydata_d(:,:)  ! Warning first dimension is supposed constant!
        INTEGER, DEVICE :: alldata_d(:,:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid, nsiz
        INTEGER, ALLOCATABLE :: nrecv(:), ndisp(:)


#if defined (__MPI)
#if ! defined(__GPU_MPI)
        INTEGER, ALLOCATABLE :: mydata_h(:,:)
        INTEGER, ALLOCATABLE :: alldata_h(:,:)
        
        ALLOCATE(mydata_h, source=mydata_d)     ! This syncs __MPI
        ALLOCATE(alldata_h, source=alldata_d)
        CALL mp_gatherv_im( mydata_h, alldata_h, recvcount, displs, root, gid)
        alldata_d = alldata_h ; mydata_d = mydata_h
        DEALLOCATE(alldata_h , mydata_h)
#else
        group = gid
        CALL mpi_comm_size( group, npe, ierr )
        IF (ierr/=0) CALL mp_stop( 9093 )
        CALL mpi_comm_rank( group, myid, ierr )
        IF (ierr/=0) CALL mp_stop( 9094 )
        !
        IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 9095 )
        IF ( myid == root ) THEN
           IF ( SIZE( alldata_d, 2 ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 9096 )
           IF ( SIZE( alldata_d, 1 ) /= SIZE( mydata_d, 1 ) ) CALL mp_stop( 9097 )
        END IF
        IF ( SIZE( mydata_d, 2 ) < recvcount( myid + 1 ) ) CALL mp_stop( 9098 )
        !
        ALLOCATE( nrecv( npe ), ndisp( npe ) )
        !
        nrecv( 1:npe ) = recvcount( 1:npe ) * SIZE( mydata_d, 1 )
        ndisp( 1:npe ) = displs( 1:npe ) * SIZE( mydata_d, 1 )
        !
        ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
        CALL MPI_GATHERV( mydata_d, nrecv( myid + 1 ), MPI_INTEGER, &
                         alldata_d, nrecv, ndisp, MPI_INTEGER, root, group, ierr )
        IF (ierr/=0) CALL mp_stop( 9099 )
        !
        DEALLOCATE( nrecv, ndisp )
        !
        RETURN ! Sync not needed after MPI call
#endif
#else
        IF ( SIZE( alldata_d, 1 ) /= SIZE( mydata_d, 1 ) ) CALL mp_stop( 9100 )
        IF ( SIZE( alldata_d, 2 ) < recvcount( 1 ) ) CALL mp_stop( 9101 )
        IF ( SIZE( mydata_d, 2  ) < recvcount( 1 ) ) CALL mp_stop( 9102 )
        !
        !alldata( :, 1:recvcount( 1 ) ) = mydata( :, 1:recvcount( 1 ) )
        
        ierr = cudaMemcpy2D(alldata_d, SIZE(alldata_d,1),&
                              mydata_d, SIZE(mydata_d,1),&
                              SIZE(mydata_d,1), recvcount( 1 ), &
                              cudaMemcpyDeviceToDevice )
        
        IF (ierr/=0) CALL mp_stop( 9103 )
#endif
        ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_gatherv_im_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_alltoall_c3d_gpu( sndbuf_d, rcvbuf_d, gid )
         IMPLICIT NONE
         COMPLEX(DP), DEVICE :: sndbuf_d( :, :, : )
         COMPLEX(DP), DEVICE :: rcvbuf_d( :, :, : )
         INTEGER, INTENT(IN) :: gid
         INTEGER :: nsiz, group, ierr, npe

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         COMPLEX(DP), ALLOCATABLE :: sndbuf_h(:,:,:)
         COMPLEX(DP), ALLOCATABLE :: rcvbuf_h(:,:,:)
         
         ALLOCATE(sndbuf_h, source=sndbuf_d)     ! This syncs __MPI
         ALLOCATE(rcvbuf_h, source=rcvbuf_d)
         CALL mp_alltoall_c3d( sndbuf_h, rcvbuf_h, gid )
         sndbuf_d = sndbuf_h ; rcvbuf_d = rcvbuf_h
         DEALLOCATE(sndbuf_h , rcvbuf_h)
#else
         group = gid
      
         CALL mpi_comm_size( group, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9104 )
      
         IF ( SIZE( sndbuf_d, 3 ) < npe ) CALL mp_stop( 9105 )
         IF ( SIZE( rcvbuf_d, 3 ) < npe ) CALL mp_stop( 9106 )
      
         nsiz = SIZE( sndbuf_d, 1 ) * SIZE( sndbuf_d, 2 )
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         !
         CALL MPI_ALLTOALL( sndbuf_d, nsiz, MPI_DOUBLE_COMPLEX, &
                            rcvbuf_d, nsiz, MPI_DOUBLE_COMPLEX, group, ierr )
      
         IF (ierr/=0) CALL mp_stop( 9107 )
         RETURN ! Sync not needed after MPI call
#endif
#else
         rcvbuf_d = sndbuf_d
#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_alltoall_c3d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_alltoall_i3d_gpu( sndbuf_d, rcvbuf_d, gid )
         IMPLICIT NONE
         INTEGER, DEVICE :: sndbuf_d( :, :, : )
         INTEGER, DEVICE :: rcvbuf_d( :, :, : )
         INTEGER, INTENT(IN) :: gid
         INTEGER :: nsiz, group, ierr, npe

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         INTEGER, ALLOCATABLE :: sndbuf_h(:,:,:)
         INTEGER, ALLOCATABLE :: rcvbuf_h(:,:,:)
         
         ALLOCATE(sndbuf_h, source=sndbuf_d)     ! This syncs __MPI
         ALLOCATE(rcvbuf_h, source=rcvbuf_d)
         CALL mp_alltoall_i3d( sndbuf_h, rcvbuf_h, gid )
         sndbuf_d = sndbuf_h ; rcvbuf_d = rcvbuf_h
         DEALLOCATE(sndbuf_h , rcvbuf_h)
#else
         group = gid
         
         CALL mpi_comm_size( group, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9108 )
         
         IF ( SIZE( sndbuf_d, 3 ) < npe ) CALL mp_stop( 9109 )
         IF ( SIZE( rcvbuf_d, 3 ) < npe ) CALL mp_stop( 9110 )
         
         nsiz = SIZE( sndbuf_d, 1 ) * SIZE( sndbuf_d, 2 )
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         !
         CALL MPI_ALLTOALL( sndbuf_d, nsiz, MPI_INTEGER, &
                            rcvbuf_d, nsiz, MPI_INTEGER, group, ierr )
         
         IF (ierr/=0) CALL mp_stop( 9111 )
         RETURN ! Sync not needed after MPI call
#endif
#else

         rcvbuf_d = sndbuf_d

#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_alltoall_i3d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_circular_shift_left_i0_gpu( buf_d, itag, gid )
         IMPLICIT NONE
         INTEGER, DEVICE :: buf_d
         INTEGER, INTENT(IN) :: itag
         INTEGER, INTENT(IN) :: gid
         INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         INTEGER :: buf_h
         buf_h = buf_d     ! This syncs __MPI
         CALL mp_circular_shift_left_i0( buf_h, itag, gid )
         buf_d = buf_h
#else
         INTEGER :: istatus( mpi_status_size )
         !
         group = gid
         !
         CALL mpi_comm_size( group, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9112 )
         CALL mpi_comm_rank( group, mype, ierr )
         IF (ierr/=0) CALL mp_stop( 9113 )
         !
         sour = mype + 1
         IF( sour == npe ) sour = 0
         dest = mype - 1
         IF( dest == -1 ) dest = npe - 1
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         CALL MPI_Sendrecv_replace( buf_d, 1, MPI_INTEGER, &
              dest, itag, sour, itag, group, istatus, ierr)
         !
         IF (ierr/=0) CALL mp_stop( 9114 )
         !
         RETURN ! Sync not needed after MPI call
#endif
#else
         ! do nothing
#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_circular_shift_left_i0_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_circular_shift_left_i1_gpu( buf_d, itag, gid )
         IMPLICIT NONE
         INTEGER, DEVICE :: buf_d(:)
         INTEGER, INTENT(IN) :: itag
         INTEGER, INTENT(IN) :: gid
         INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         INTEGER, ALLOCATABLE :: buf_h(:)
         ALLOCATE(buf_h, source=buf_d)    ! This syncs __MPI
         CALL mp_circular_shift_left_i1( buf_h, itag, gid )
         buf_d = buf_h; DEALLOCATE(buf_h)
#else
         INTEGER :: istatus( mpi_status_size )
         !
         group = gid
         !
         CALL mpi_comm_size( group, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9115 )
         CALL mpi_comm_rank( group, mype, ierr )
         IF (ierr/=0) CALL mp_stop( 9116 )
         !
         sour = mype + 1
         IF( sour == npe ) sour = 0
         dest = mype - 1
         IF( dest == -1 ) dest = npe - 1
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         CALL MPI_Sendrecv_replace( buf_d, SIZE(buf_d), MPI_INTEGER, &
              dest, itag, sour, itag, group, istatus, ierr)
         !
         IF (ierr/=0) CALL mp_stop( 9117 )
         !
         RETURN ! Sync not needed after MPI call
#endif
#else
         ! do nothing
#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_circular_shift_left_i1_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_circular_shift_left_i2_gpu( buf_d, itag, gid )
         IMPLICIT NONE
         INTEGER, DEVICE :: buf_d(:,:)
         INTEGER, INTENT(IN) :: itag
         INTEGER, INTENT(IN) :: gid
         INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         INTEGER, ALLOCATABLE :: buf_h(:,:)
         ALLOCATE(buf_h, source=buf_d)    ! This syncs __MPI
         CALL mp_circular_shift_left_i2( buf_h, itag, gid )
         buf_d = buf_h; DEALLOCATE(buf_h)
#else
         INTEGER :: istatus( mpi_status_size )
         !
         group = gid
         !
         CALL mpi_comm_size( group, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9118 )
         CALL mpi_comm_rank( group, mype, ierr )
         IF (ierr/=0) CALL mp_stop( 9119 )
         !
         sour = mype + 1
         IF( sour == npe ) sour = 0
         dest = mype - 1
         IF( dest == -1 ) dest = npe - 1
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         CALL MPI_Sendrecv_replace( buf_d, SIZE(buf_d), MPI_INTEGER, &
              dest, itag, sour, itag, group, istatus, ierr)
         !
         IF (ierr/=0) CALL mp_stop( 9120 )
         !
         RETURN ! Sync not needed after MPI call
#endif
#else
         ! do nothing
#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_circular_shift_left_i2_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_circular_shift_left_r2d_gpu( buf_d, itag, gid )
         IMPLICIT NONE
         REAL(DP), DEVICE :: buf_d( :, : )
         INTEGER, INTENT(IN) :: itag
         INTEGER, INTENT(IN) :: gid
         INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         REAL(DP), ALLOCATABLE :: buf_h(:, :)
         ALLOCATE(buf_h, source=buf_d)    ! This syncs __MPI
         CALL mp_circular_shift_left_r2d( buf_h, itag, gid )
         buf_d = buf_h; DEALLOCATE(buf_h)
#else
         INTEGER :: istatus( mpi_status_size )
         !
         group = gid
         !
         CALL mpi_comm_size( group, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9121 )
         CALL mpi_comm_rank( group, mype, ierr )
         IF (ierr/=0) CALL mp_stop( 9122 )
         !
         sour = mype + 1
         IF( sour == npe ) sour = 0
         dest = mype - 1
         IF( dest == -1 ) dest = npe - 1
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         CALL MPI_Sendrecv_replace( buf_d, SIZE(buf_d), MPI_DOUBLE_PRECISION, &
              dest, itag, sour, itag, group, istatus, ierr)
         !
         IF (ierr/=0) CALL mp_stop( 9123 )
         !
         RETURN ! Sync not needed after MPI call
#endif
#else
         ! do nothing
#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_circular_shift_left_r2d_gpu
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_circular_shift_left_c2d_gpu( buf_d, itag, gid )
         IMPLICIT NONE
         COMPLEX(DP), DEVICE :: buf_d( :, : )
         INTEGER, INTENT(IN) :: itag
         INTEGER, INTENT(IN) :: gid
         INTEGER :: nsiz, group, ierr, npe, sour, dest, mype

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         COMPLEX(DP), ALLOCATABLE :: buf_h(:, :)
         ALLOCATE(buf_h, source=buf_d)    ! This syncs __MPI
         CALL mp_circular_shift_left_c2d( buf_h, itag, gid )
         buf_d = buf_h; DEALLOCATE(buf_h)
#else
         INTEGER :: istatus( mpi_status_size )
         !
         group = gid
         !
         CALL mpi_comm_size( group, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9124 )
         CALL mpi_comm_rank( group, mype, ierr )
         IF (ierr/=0) CALL mp_stop( 9125 )
         !
         sour = mype + 1
         IF( sour == npe ) sour = 0
         dest = mype - 1
         IF( dest == -1 ) dest = npe - 1
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         CALL MPI_Sendrecv_replace( buf_d, SIZE(buf_d), MPI_DOUBLE_COMPLEX, &
              dest, itag, sour, itag, group, istatus, ierr)
         !
         IF (ierr/=0) CALL mp_stop( 9126 )
         !
         RETURN ! Sync not needed after MPI call
#endif
#else
         ! do nothing
#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_circular_shift_left_c2d_gpu

!------------------------------------------------------------------------------!
!..mp_gatherv_inplace_cplx_array
!

      SUBROUTINE mp_gatherv_inplace_cplx_array_gpu(alldata_d, my_column_type, recvcount, displs, root, gid)
         IMPLICIT NONE
         COMPLEX(DP), DEVICE :: alldata_d(:,:)
         INTEGER, INTENT(IN) :: my_column_type
         INTEGER, INTENT(IN) :: recvcount(:), displs(:)
         INTEGER, INTENT(IN) :: root, gid
         INTEGER :: ierr, npe, myid

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         COMPLEX(DP), ALLOCATABLE :: alldata_h(:, :)
         ALLOCATE(alldata_h, source=alldata_d)    ! This syncs __MPI
         CALL mp_gatherv_inplace_cplx_array(alldata_h, my_column_type, recvcount, displs, root, gid)
         alldata_d = alldata_h; DEALLOCATE(alldata_h)
#else
         CALL mpi_comm_size( gid, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9127 )
         CALL mpi_comm_rank( gid, myid, ierr )
         IF (ierr/=0) CALL mp_stop( 9128 )
         !
         IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 9129 )
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         IF (myid==root) THEN
            CALL MPI_GATHERV( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                              alldata_d, recvcount, displs, my_column_type, root, gid, ierr )
         ELSE
            CALL MPI_GATHERV( alldata_d(1,displs(myid+1)+1), recvcount(myid+1), my_column_type, &
                              MPI_IN_PLACE, recvcount, displs, MPI_DATATYPE_NULL, root, gid, ierr )
         ENDIF
         !
         IF (ierr/=0) CALL mp_stop( 9130 )
         !
         RETURN ! Sync not needed after MPI call
#endif
#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_gatherv_inplace_cplx_array_gpu

!------------------------------------------------------------------------------!
!..mp_allgatherv_inplace_cplx_array
!

      SUBROUTINE mp_allgatherv_inplace_cplx_array_gpu(alldata_d, my_element_type, recvcount, displs, gid)
         IMPLICIT NONE
         COMPLEX(DP), DEVICE :: alldata_d(:,:)
         INTEGER, INTENT(IN) :: my_element_type
         INTEGER, INTENT(IN) :: recvcount(:), displs(:)
         INTEGER, INTENT(IN) :: gid
         INTEGER :: ierr, npe, myid

#if defined (__MPI)
#if ! defined(__GPU_MPI)
         COMPLEX(DP), ALLOCATABLE :: alldata_h(:, :)
         ALLOCATE(alldata_h, source=alldata_d)! This syncs __MPI
         CALL mp_allgatherv_inplace_cplx_array(alldata_h, my_element_type, recvcount, displs, gid)
         alldata_d = alldata_h; DEALLOCATE(alldata_h)
#else
         CALL mpi_comm_size( gid, npe, ierr )
         IF (ierr/=0) CALL mp_stop( 9131 )
         CALL mpi_comm_rank( gid, myid, ierr )
         IF (ierr/=0) CALL mp_stop( 9132 )
         !
         IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 9133 )
         !
         ierr = cudaDeviceSynchronize()   ! This syncs __GPU_MPI
         CALL MPI_ALLGATHERV( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                              alldata_d, recvcount, displs, my_element_type, gid, ierr )
         IF (ierr/=0) CALL mp_stop( 9134 )
         RETURN ! Sync not needed after MPI call
#endif
#endif
         ierr = cudaDeviceSynchronize()   ! This syncs SERIAL, __MPI
      END SUBROUTINE mp_allgatherv_inplace_cplx_array_gpu

      SUBROUTINE mp_type_create_cplx_column_section_gpu(dummy, start, length, stride, mytype)
         IMPLICIT NONE
         !
         COMPLEX (DP), DEVICE, INTENT(IN) :: dummy
         INTEGER, INTENT(IN) :: start, length, stride
         INTEGER, INTENT(OUT) :: mytype
         !
#if defined(__MPI)
         INTEGER :: ierr
         !
         CALL MPI_TYPE_CREATE_SUBARRAY(1, stride, length, start, MPI_ORDER_FORTRAN,&
                                       MPI_DOUBLE_COMPLEX, mytype, ierr)
         IF (ierr/=0) CALL mp_stop( 8081 )
         CALL MPI_Type_commit(mytype, ierr)
         IF (ierr/=0) CALL mp_stop( 8082 )
#else
         mytype = 0;
#endif
         !
         RETURN
      END SUBROUTINE mp_type_create_cplx_column_section_gpu

!------------------------------------------------------------------------------!

#endif
!------------------------------------------------------------------------------!
    END MODULE mp
!------------------------------------------------------------------------------!

! Script to generate stop messages:
!   # coding: utf-8
!   import re
!   import sys
!   i = 8000
!   def replace(match):
!       global i
!       i += 1
!       return 'mp_stop( {0} )'.format(i)
!   
!   with open(sys.argv[1],'r') as f:
!       data = re.sub(r"mp_stop\(\s?\d+\s?\)", replace, f.read())
!       with open(sys.argv[1]+'.new','w') as fo:
!           fo.write(data)
