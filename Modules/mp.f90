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

!      PRIVATE
      PUBLIC :: mp_start, mp_end, mp_env, mp_group, mp_cart_create, &
        mp_bcast, mp_stop, mp_sum, mp_max, mp_min, mp_rank, mp_size, &
        mp_excng, mp_gather, mp_get, mp_put, mp_barrier, mp_report
!
      INTERFACE mp_excng    ! Carlo Cavazzoni
        MODULE PROCEDURE mp_excng_i
      END INTERFACE

      INTERFACE mp_bcast
#if defined __T3E
        MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, &
          mp_bcast_z, mp_bcast_zv, &
          mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
          mp_bcast_cm, mp_bcast_im, mp_bcast_it, mp_bcast_rt, mp_bcast_lv, &
          mp_bcast_lm, mp_bcast_r4d, mp_bcast_ct, mp_bcast_c4d, &
          mp_bcast_i4b
#else
        MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, &
          mp_bcast_z, mp_bcast_zv, &
          mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
          mp_bcast_cm, mp_bcast_im, mp_bcast_it, mp_bcast_rt, mp_bcast_lv, &
          mp_bcast_lm, mp_bcast_r4d, mp_bcast_ct,  mp_bcast_c4d
#endif
      END INTERFACE

      INTERFACE mp_sum
        MODULE PROCEDURE mp_sum_i1, mp_sum_iv, mp_sum_im, mp_sum_it, & 
          mp_sum_r1, mp_sum_rv, mp_sum_rm, mp_sum_rt, &
          mp_sum_c1, mp_sum_cv, mp_sum_cm, mp_sum_ct, mp_sum_c4d, &
          mp_sum_rmm, mp_sum_cmm
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
        MODULE PROCEDURE mp_gather_iv, mp_gather_cvv, mp_gather_rvv
      END INTERFACE

      INTEGER, PRIVATE, SAVE :: mp_high_watermark = 0

      INTEGER, PRIVATE, PARAMETER :: mp_msgsiz_max = 900000000

      CHARACTER(LEN=80), PRIVATE :: err_msg = ' '

!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
!..mp_excng
!..Carlo Cavazzoni
!     THIS SUBROUTINE performs the following operation :
!
!               ARRAY                      ARRAY
!     P0   [D0][  ][  ][  ]           [D0][D1][D2][D3]
!     P1   [  ][D1][  ][  ]    --\    [D0][D1][D2][D3]
!     P2   [  ][  ][D2][  ]    --/    [D0][D1][D2][D3]
!     P3   [  ][  ][  ][D3]           [D0][D1][D2][D3]

      SUBROUTINE mp_excng_i(mydata, alldata, gid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: mydata
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER, INTENT(OUT) :: alldata(:)
        INTEGER :: taskid, ierr
        INTEGER :: msglen = 1

#if defined (__MPI)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8900)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_comm_rank(group, taskid, ierr)
        IF (ierr/=0) CALL mp_stop(8001)
        alldata(taskid+1) = mydata
        CALL MPI_ALLGATHER(mydata, 1, MPI_INTEGER, alldata, 1, MPI_INTEGER, group, IERR)
        IF (ierr/=0) CALL mp_stop(8001)
#else
        alldata(1) = mydata
#endif
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
        RETURN
      END SUBROUTINE mp_excng_i

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
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8901)
        IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop(8002)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL MPI_GATHER(mydata, msglen, MPI_INTEGER, alldata, msglen, MPI_INTEGER, root, group, IERR)
        IF (ierr/=0) CALL mp_stop(8001)
#else
        msglen = SIZE(mydata)
        IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop(8002)
        alldata(:,1) = mydata(:)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
        RETURN
      END SUBROUTINE mp_gather_iv

!------------------------------------------------------------------------------!
!..mp_gather_cvv
!..Carlo Cavazzoni
      SUBROUTINE mp_gather_cvv(mydata, alldata, root, gid)
        IMPLICIT NONE
        COMPLEX(DP), INTENT(IN) :: mydata(:)
        COMPLEX(DP), INTENT(OUT) :: alldata(:)
        INTEGER, INTENT(IN) :: root
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
#if defined (__MPI)
        group = MPI_COMM_WORLD
        msglen = SIZE(mydata)
        IF( msglen*16 > mp_msgsiz_max ) CALL mp_stop(8902)
        IF( PRESENT( gid ) ) group = gid
        CALL MPI_GATHER(mydata, msglen, MPI_DOUBLE_COMPLEX, alldata, msglen, &
          MPI_DOUBLE_COMPLEX, root, group, IERR)
        IF (ierr/=0) CALL mp_stop(8001)
#else
        msglen = SIZE(mydata)
        alldata = mydata
#endif
        mp_high_watermark = MAX( mp_high_watermark, 16 * msglen ) 
        RETURN
      END SUBROUTINE mp_gather_cvv

!------------------------------------------------------------------------------!
!..mp_gather_rvv
!..Carlo Cavazzoni
      SUBROUTINE mp_gather_rvv(mydata, alldata, root, gid)
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: mydata(:)
        REAL(DP), INTENT(OUT) :: alldata(:)
        INTEGER, INTENT(IN) :: root
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
#if defined (__MPI)
        group = mpi_comm_world
        msglen = SIZE(mydata)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8903)
        IF( PRESENT( gid ) ) group = gid
        CALL MPI_GATHER(mydata, msglen, MPI_DOUBLE_PRECISION, alldata, msglen, &
          MPI_DOUBLE_PRECISION, root, group, IERR)
        IF (ierr/=0) CALL mp_stop(8001)
#else
        msglen = SIZE(mydata)
        alldata = mydata
#endif
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
        return
      END SUBROUTINE mp_gather_rvv


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

#if defined(__MPI) || defined (__SHMEM)
        CALL MPI_INIT(ierr)
        IF (ierr/=0) CALL mp_stop(8000)
#endif

#if defined __HPM && defined __AIX

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

#if defined __HPM && defined __AIX

        !   terminate the IBM Harware performance monitor

#  if defined(__MPI)
        CALL mpi_comm_rank( mpi_comm_world, taskid, ierr)
#  endif
        CALL f_hpmterminate( taskid )
#endif

#if defined(__MPI)
        CALL mpi_finalize(ierr)
        IF (ierr/=0) CALL mp_stop(8904)
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
        IF (ierr/=0) CALL mp_stop(8001)
        CALL mpi_comm_size(mpi_comm_world,numtask,ierr)
        groupid = mpi_comm_world
        IF (ierr/=0) CALL mp_stop(8002)

#endif

        RETURN
      END SUBROUTINE mp_env
!------------------------------------------------------------------------------!
!..mp_group
      SUBROUTINE mp_group(group_list, group_size, base_group, groupid)
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: group_list(:), group_size, base_group
        INTEGER, INTENT (OUT) :: groupid
        INTEGER :: base, newgroup, ierr

        ierr = 0
        groupid = base_group
#if defined(__MPI)
        CALL mpi_comm_group(base_group,base,ierr)
        IF (ierr/=0) CALL mp_stop(8010)
        CALL mpi_group_incl(base,group_size,group_list,newgroup,ierr)
        IF (ierr/=0) CALL mp_stop(8011)
        CALL mpi_comm_create(base_group,newgroup,groupid,ierr)
        IF (ierr/=0) CALL mp_stop(8012)
#endif
      END SUBROUTINE mp_group
!------------------------------------------------------------------------------!
!..mp_cart_create
      SUBROUTINE mp_cart_create(comm_old,ndims,dims,pos,comm_cart)
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: comm_old, ndims
        INTEGER, INTENT (OUT) :: dims(:), pos(:), comm_cart
        INTEGER :: ierr, nodes
        LOGICAL :: period(1:ndims), reorder

        ierr = 0
        dims(1:ndims) = 1
        pos(1:ndims) = 1
        comm_cart = comm_old
#if defined(__MPI)
        dims(1:ndims) = 0
        CALL mpi_comm_size(comm_old,nodes,ierr)
        IF (ierr/=0) CALL mp_stop(8020)
        CALL mpi_dims_create(nodes,ndims,dims,ierr)
        IF (ierr/=0) CALL mp_stop(8021)
        reorder = .TRUE.
        period = .TRUE.
        CALL mpi_cart_create(comm_old,ndims,dims,period,reorder,comm_cart, ierr)
        IF (ierr/=0) CALL mp_stop(8022)
        CALL mpi_cart_get(comm_cart,ndims,dims,period,pos,ierr)
        IF (ierr/=0) CALL mp_stop(8023)
#endif
      END SUBROUTINE mp_cart_create
!------------------------------------------------------------------------------!
!..mp_bcast
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
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8905)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        imsg = msg
        CALL mpi_bcast(imsg, msglen, mpi_integer, source, group, ierr)
        msg = imsg
        IF (ierr/=0) CALL mp_stop(8101)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_bcast_i4b


!------------------------------------------------------------------------------!
!..mp_bcast
      SUBROUTINE mp_bcast_i1(msg,source,gid)
        IMPLICIT NONE
        INTEGER :: msg
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr

#if defined(__MPI)
        ierr = 0
        msglen = 1
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8906)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_integer,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8101)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8907)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_integer,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8102)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_bcast_iv
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_im(msg,source,gid)
        IMPLICIT NONE
        INTEGER :: msg(:,:)
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8908)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_integer,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8102)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8909)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_integer,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8102)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        INTEGER :: ierr
#if defined(__MPI)
        ierr = 0
        msglen = 1
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8910)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_precision,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8111)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        INTEGER :: msglen, ierr

#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8911)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_precision,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8112)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8912)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_precision,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8113)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8913)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_precision,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8113)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8914)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg, msglen, mpi_double_precision, source, group, ierr)
        IF (ierr/=0) CALL mp_stop(8113)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_bcast_r4d

!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_c1(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg
        INTEGER :: source
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_complex,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8121)
        mp_high_watermark = MAX( mp_high_watermark, 16 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*16 > mp_msgsiz_max ) CALL mp_stop(8916)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_complex,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8122)
        mp_high_watermark = MAX( mp_high_watermark, 16 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*16 > mp_msgsiz_max ) CALL mp_stop(8915)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_complex,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8123)
        mp_high_watermark = MAX( mp_high_watermark, 16 * msglen )
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*16 > mp_msgsiz_max ) CALL mp_stop(8915)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_complex,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8123)
        mp_high_watermark = MAX( mp_high_watermark, 16 * msglen )
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        IF( msglen*16 > mp_msgsiz_max ) CALL mp_stop(8915)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_double_complex,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8123)
        mp_high_watermark = MAX( mp_high_watermark, 16 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_bcast(msg,msglen,mpi_logical,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8130)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8916)
        CALL mpi_bcast(msg,msglen,mpi_logical,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8130)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        INTEGER :: msglen, ierr
#if defined(__MPI)
        ierr = 0
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8916)
        CALL mpi_bcast(msg,msglen,mpi_logical,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8130)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen )
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
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8917)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
! this is a workaround to avoid problems on the T3E
! at the moment we have a data alignment error when trying to
! broadcast characters on the T3E (not always!)
! JH 3/19/99 on galileo
!       CALL mpi_bcast(msg,msglen,mpi_character,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8140)
        ALLOCATE (imsg(1:msglen), STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8140)
        DO i = 1, msglen
          imsg(i) = ichar(msg(i:i))
        END DO
        CALL mpi_bcast(imsg,msglen,mpi_integer,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8140)
        DO i = 1, msglen
          msg(i:i) = char(imsg(i))
        END DO
        DEALLOCATE (imsg, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8140)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8917)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
! ...      CALL mpi_bcast(msg,msglen,mpi_character,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8140)
        ALLOCATE (imsg(1:m1,1:m2), STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8140)
        DO j = 1, m2
          DO i = 1, m1
            imsg(i,j) = ichar(msg(j)(i:i))
          END DO
        END DO
        CALL mpi_bcast(imsg,msglen,mpi_integer,source,group,ierr)
        IF (ierr/=0) CALL mp_stop(8140)
        DO j = 1, m2
          DO i = 1, m1
            msg(j)(i:i) = char(imsg(i,j))
          END DO
        END DO
        DEALLOCATE (imsg, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8140)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen )
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

        IF(dest .NE. sour) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          CALL MPI_SEND( msg_sour, msglen, MPI_INTEGER, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, msglen, MPI_INTEGER, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
        END IF
#endif
        ELSE
          msg_dest = msg_sour
        END IF

#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 

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

        IF(sour .NE. dest) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          msglen = SIZE(msg_sour)
          CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_INTEGER, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_INTEGER, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = nrcv
        END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8918)
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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

        IF(sour .NE. dest) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          msglen = SIZE(msg_sour) 
          CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = nrcv
        END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
        IF( msglen*8 > mp_msgsiz_max ) THEN
          WRITE( err_msg, * ) "Message too long, ", msglen
          CALL mp_stop(8919)
        END IF
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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

        IF(sour .NE. dest) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = SIZE(msg_sour)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = nrcv
        END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour,1), 1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE( msg_sour )
        END IF
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8920)
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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

        IF( dest .NE. sour ) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = SIZE(msg_sour)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = nrcv
        END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
        IF( msglen*16 > mp_msgsiz_max ) CALL mp_stop(8921)
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 16 * msglen ) 
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

        IF(dest .NE. sour) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          CALL MPI_SEND( msg_sour, 1, MPI_INTEGER, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = 1
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, 1, MPI_INTEGER, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = 1
        END IF
#endif
        ELSE
          msg_dest = msg_sour
          msglen = 1
        END IF
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8922)
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        IF(sour .NE. dest) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_INTEGER, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = SIZE(msg_sour)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_INTEGER, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_INTEGER, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = nrcv
        END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8923)
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        IF(sour .NE. dest) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = SIZE(msg_sour)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = nrcv
        END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8924)
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        IF(sour .NE. dest) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_PRECISION, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = SIZE(msg_sour)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_PRECISION, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_PRECISION, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = nrcv
        END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour,1),1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE(msg_sour)
        END IF
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8925)
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        IF( dest .NE. sour ) THEN
#if defined(__MPI)
        IF(mpime .EQ. sour) THEN
          CALL MPI_SEND( msg_sour, SIZE(msg_sour), MPI_DOUBLE_COMPLEX, dest, ip, group, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = SIZE(msg_sour)
        END IF
        IF(mpime .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, SIZE(msg_dest), MPI_DOUBLE_COMPLEX, sour, ip, group, istatus, IERR )
          IF (ierr/=0) CALL mp_stop(8140)
          CALL MPI_GET_COUNT(istatus, MPI_DOUBLE_COMPLEX, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(8140)
          msglen = nrcv
        END IF
#endif
        ELSE
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
        IF( msglen*16 > mp_msgsiz_max ) CALL mp_stop(8926)
#if defined(__MPI)
        CALL MPI_BARRIER(group, IERR)
        IF (ierr/=0) CALL mp_stop(8140)
#endif
        mp_high_watermark = MAX( mp_high_watermark, 16 * msglen ) 
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
        INTEGER :: msglen, res, ierr

        msglen = 1
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_allreduce(msg,res,msglen,mpi_integer,mpi_sum,group,ierr)
        IF (ierr/=0) CALL mp_stop(8200)
        msg = res
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_sum_i1
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        INTEGER, ALLOCATABLE :: res(:)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8927)
        ALLOCATE (res(1:msglen),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8201)
        CALL mpi_allreduce(msg,res,msglen,mpi_integer,mpi_sum,group,ierr)
        IF (ierr/=0) CALL mp_stop(8200)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8202)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_sum_iv
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_im(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, m1, m2, ierr
        INTEGER, ALLOCATABLE :: res(:,:)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8928)
        m1 = size(msg(:,1))
        m2 = size(msg(1,:))
        ALLOCATE (res(m1,m2),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8204)
        CALL mpi_allreduce(msg,res,msglen,mpi_integer,mpi_sum,group,ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_sum_im
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_it(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, m1, m2, m3, ierr
        INTEGER, ALLOCATABLE :: res(:,:,:)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8929)
        m1 = size(msg,1)
        m2 = size(msg,2)
        m3 = size(msg,3)
        ALLOCATE (res(m1,m2,m3),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8204)
        CALL mpi_allreduce(msg,res,msglen,mpi_integer,mpi_sum,group,ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_sum_it

!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_r1(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        REAL (DP) :: res
#if defined(__MPI)
        msglen = 1
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,group,ierr)
        IF (ierr/=0) CALL mp_stop(8203)
        msg = res
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_sum_r1

!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        REAL (DP), ALLOCATABLE :: res(:)
#if defined(__MPI)
        msglen = size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8930)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        ALLOCATE (res(1:msglen),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8204)
        CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,group, ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_sum_rv
!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_rm(msg, gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, m1, m2, ierr, i, j, k
        REAL (DP), ALLOCATABLE :: res(:,:)
        REAL (DP), ALLOCATABLE :: resv(:)
#if defined(__MPI)
        msglen = size(msg)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid

        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8931)
        m1 = size(msg(:,1))
        m2 = size(msg(1,:))
        ALLOCATE (res(m1,m2),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8204)
        CALL mpi_allreduce(msg, res, msglen, mpi_double_precision, mpi_sum, group, ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8205)

        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 

#endif

      END SUBROUTINE mp_sum_rm
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_rmm(msg, res, root, gid)
        IMPLICIT NONE
        REAL (DP), INTENT (IN) :: msg(:,:)
        REAL (DP), INTENT (OUT) :: res(:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: root
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr, j

        msglen = size(msg)

#if defined(__MPI)

        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid

        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8932)

        IF( PRESENT( root ) ) THEN
          CALL mpi_reduce(msg, res, msglen, mpi_double_precision, mpi_sum, root, group, ierr)
        ELSE
          CALL mpi_allreduce(msg, res, msglen, mpi_double_precision, mpi_sum, group, ierr)
        END IF
        IF (ierr/=0) CALL mp_stop(8205)

#else
        res = msg
#endif
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 

      END SUBROUTINE mp_sum_rmm


!
!------------------------------------------------------------------------------!


      SUBROUTINE mp_sum_rt(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, m1, m2, m3, ierr
        REAL (DP), ALLOCATABLE :: res(:,:,:)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8933)
        m1 = size(msg,1)
        m2 = size(msg,2)
        m3 = size(msg,3)
        ALLOCATE (res(m1,m2,m3),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8204)
        CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,group, ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8205)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_sum_rt

!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_c1(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        COMPLEX (DP) :: res

#if defined(__MPI)
        msglen = 2
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,group,ierr)
        msg = res
        IF (ierr/=0) CALL mp_stop(8205)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_sum_c1
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_cv(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        COMPLEX (DP), ALLOCATABLE :: res(:)
#if defined(__MPI)
        msglen = 2*size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8934)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        ALLOCATE (res(1:size(msg)),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8206)
        CALL mpi_allreduce(msg(1),res(1),msglen,mpi_double_precision,mpi_sum,group, ierr)
        IF (ierr/=0) CALL mp_stop(8207)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8207)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_sum_cv
!
!------------------------------------------------------------------------------!

      SUBROUTINE mp_sum_cm(msg, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, OPTIONAL, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, m1, m2, ierr
        COMPLEX (DP), ALLOCATABLE :: res(:,:)
#if defined(__MPI)
        msglen = 2*size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8935)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        m1 = size(msg(:,1))
        m2 = size(msg(1,:))
        ALLOCATE (res(m1,m2),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8204)
        CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,group, ierr)
        IF (ierr/=0) CALL mp_stop(8208)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8208)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        INTEGER :: msglen, ierr
        msglen = 2*size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8936)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        CALL mpi_allreduce(msg, res, msglen, mpi_double_precision, mpi_sum, group, ierr)
        IF (ierr/=0) CALL mp_stop(8208)
#else
        res = msg
#endif
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        INTEGER :: i, msglen, ierr
        COMPLEX (DP), ALLOCATABLE :: res(:,:,:)
#if defined(__MPI)
        msglen = 2 * SIZE(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8937)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        ALLOCATE (res(size(msg,1),size(msg,2),size(msg,3)),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8204)
        CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,group,ierr)
        IF (ierr/=0) CALL mp_stop(8208)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8208)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        INTEGER :: i, msglen, ierr
        COMPLEX (DP), ALLOCATABLE :: res(:,:,:,:)
#if defined(__MPI)
        msglen = 2*size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8938)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        ALLOCATE (res(size(msg,1),size(msg,2),size(msg,3),size(msg,4)), STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8204)
        CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,group, ierr)
        IF (ierr/=0) CALL mp_stop(8208)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8208)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_sum_c4d


!------------------------------------------------------------------------------!
      SUBROUTINE mp_max_i(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        INTEGER :: res
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = 1
        CALL MPI_ALLREDUCE(MSG,res,msglen,MPI_INTEGER,MPI_MAX,group,IERR)
        IF (ierr/=0) CALL mp_stop(8300)
        msg = res
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
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
        INTEGER :: msglen, ierr
        INTEGER, ALLOCATABLE :: res(:)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8939)
        ALLOCATE (res(1:msglen),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8302)
        CALL MPI_ALLREDUCE(MSG,res,msglen,MPI_INTEGER, MPI_MAX,group,IERR)
        IF (ierr/=0) CALL mp_stop(8303)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8303)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_max_iv
!
!----------------------------------------------------------------------

      SUBROUTINE mp_max_r(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        REAL (DP) :: res
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = 1
        CALL MPI_ALLREDUCE(MSG,res,msglen,MPI_DOUBLE_PRECISION, MPI_MAX,group,IERR)
        IF (ierr/=0) CALL mp_stop(8301)
        msg = res
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_max_r
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_max_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        REAL (DP), ALLOCATABLE :: res(:)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8940)
        ALLOCATE (res(1:msglen),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8302)
        CALL MPI_ALLREDUCE(MSG,res,msglen,MPI_DOUBLE_PRECISION, MPI_MAX,group,IERR)
        IF (ierr/=0) CALL mp_stop(8303)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8303)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_max_rv
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_i(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        INTEGER :: res
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = 1
        CALL MPI_ALLREDUCE(MSG,res,msglen,MPI_INTEGER,MPI_MIN,group,IERR)
        IF (ierr/=0) CALL mp_stop(8310)
        msg = res
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_min_i
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        INTEGER, ALLOCATABLE :: res(:)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = SIZE(msg)
        IF( msglen*4 > mp_msgsiz_max ) CALL mp_stop(8941)
        ALLOCATE (res(1:msglen),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8312)
        CALL MPI_ALLREDUCE(MSG,res,msglen,MPI_INTEGER,MPI_MIN,group,IERR)
        IF (ierr/=0) CALL mp_stop(8313)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8313)
        mp_high_watermark = MAX( mp_high_watermark, 4 * msglen ) 
#endif
      END SUBROUTINE mp_min_iv
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_r(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        REAL (DP) :: res
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = 1
        CALL MPI_ALLREDUCE(MSG,res,msglen,MPI_DOUBLE_PRECISION, MPI_MIN,group,IERR)
        IF (ierr/=0) CALL mp_stop(8311)
        msg = res
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
#endif
      END SUBROUTINE mp_min_r
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, OPTIONAL, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr
        REAL (DP), ALLOCATABLE :: res(:)
#if defined(__MPI)
        group = mpi_comm_world
        IF( PRESENT( gid ) ) group = gid
        msglen = size(msg)
        IF( msglen*8 > mp_msgsiz_max ) CALL mp_stop(8942)
        ALLOCATE (res(1:msglen),STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8312)
        CALL MPI_ALLREDUCE(MSG,res,msglen,MPI_DOUBLE_PRECISION, MPI_MIN,group,IERR)
        IF (ierr/=0) CALL mp_stop(8313)
        msg = res
        DEALLOCATE (res, STAT=ierr)
        IF (ierr/=0) CALL mp_stop(8313)
        mp_high_watermark = MAX( mp_high_watermark, 8 * msglen ) 
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
        IF (ierr/=0) CALL mp_stop(8313)
#endif
      END SUBROUTINE mp_barrier

!------------------------------------------------------------------------------!
!.. Carlo Cavazzoni
!..mp_rank
      FUNCTION mp_rank()
        IMPLICIT NONE
        INTEGER :: mp_rank
        INTEGER :: ierr, taskid

        ierr = 0
        taskid = 0
#if defined(__MPI)
        CALL mpi_comm_rank(mpi_comm_world,taskid,ierr)
        IF (ierr/=0) CALL mp_stop(8003)
#endif
        mp_rank = taskid
      END FUNCTION mp_rank

!------------------------------------------------------------------------------!
!.. Carlo Cavazzoni
!..mp_size
      FUNCTION mp_size()
        IMPLICIT NONE
        INTEGER :: mp_size
        INTEGER :: ierr, numtask

        ierr = 0
        numtask = 1
#if defined(__MPI)
        CALL mpi_comm_size(mpi_comm_world,numtask,ierr)
        IF (ierr/=0) CALL mp_stop(8004)
#endif
        mp_size = numtask
      END FUNCTION mp_size

      SUBROUTINE mp_report
        WRITE( stdout, *) 
        WRITE( stdout, *) '  mp: high_watermark (bytes): ', mp_high_watermark
        RETURN
      END SUBROUTINE mp_report


!------------------------------------------------------------------------------!
    END MODULE mp
!------------------------------------------------------------------------------!
