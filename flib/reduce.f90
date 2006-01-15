!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! FIXME: this seems like a workaround for the MPICH defaults.
!  NOTE: that this can worked around via environment variables.
! AK 2005/10/18
!
#if defined(__LINUX) || defined(__LINUX64)
#  define MESSAGE_MAX_SIZE         262144
#  define INTEGER_MESSAGE_MAX_SIZE  65536
#  define REAL_MESSAGE_MAX_SIZE     32768
#  define COMPLEX_MESSAGE_MAX_SIZE  16384
#else
#  define MESSAGE_MAX_SIZE          20000000
#  define INTEGER_MESSAGE_MAX_SIZE   5000000
#  define REAL_MESSAGE_MAX_SIZE      2500000
#  define COMPLEX_MESSAGE_MAX_SIZE   1250000
#endif


      SUBROUTINE PARALLEL_SUM_REAL(ARRAY,N, gid)
        IMPLICIT NONE
        INTEGER N, ERR, I, gid
        REAL*8 ARRAY(N)

#if defined __MPI 

        INCLUDE 'mpif.h'
        INTEGER nblk, blksiz, msgsiz, iblk, istart
        INTEGER msgsiz_max
        REAL*8 ARRAY_TMP(N)
        msgsiz_max = REAL_MESSAGE_MAX_SIZE

        IF( n .LT. msgsiz_max ) THEN
          CALL MPI_ALLREDUCE(ARRAY, ARRAY_TMP, N, MPI_DOUBLE_PRECISION,  &
     &       MPI_SUM, gid, ERR)
        ELSE
          nblk   = n / msgsiz_max
          blksiz = msgsiz_max
          DO iblk = 1, nblk
            istart = (iblk-1)*msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY(istart), ARRAY_TMP(istart), blksiz, &
     &         MPI_DOUBLE_PRECISION, MPI_SUM, gid, ERR)
          END DO
          blksiz = MOD( n, msgsiz_max )
          IF( blksiz .GT. 0 ) THEN
            istart = nblk * msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY(istart), ARRAY_TMP(istart), blksiz, &
     &         MPI_DOUBLE_PRECISION, MPI_SUM, gid, ERR)
          END IF
        END IF
        ARRAY = ARRAY_TMP
#endif
        RETURN
      END SUBROUTINE PARALLEL_SUM_REAL

!=----------------------------------------------------------------------------=!

      SUBROUTINE PARALLEL_SUM_COMPLEX(ARRAY,N,gid)
        IMPLICIT NONE
        INTEGER N, ERR, I, gid
        COMPLEX*16 ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart
        COMPLEX*16 ARRAY_TMP(N)
        msgsiz_max = COMPLEX_MESSAGE_MAX_SIZE 

        IF( n .LT. msgsiz_max ) THEN
          CALL MPI_ALLREDUCE(ARRAY, ARRAY_TMP, N, MPI_DOUBLE_COMPLEX,  &
     &       MPI_SUM, gid, ERR)
        ELSE
          nblk   = n / msgsiz_max
          blksiz = msgsiz_max
          DO iblk = 1, nblk
            istart = (iblk-1)*msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY(istart), ARRAY_TMP(istart), blksiz, &
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, gid, ERR)
          END DO
          blksiz = MOD( n, msgsiz_max )
          IF( blksiz .GT. 0 ) THEN
            istart = nblk * msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY(istart), ARRAY_TMP(istart), blksiz, &
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, gid, ERR)
          END IF
        END IF

        ARRAY = ARRAY_TMP
#endif
        RETURN
      END SUBROUTINE PARALLEL_SUM_COMPLEX

!=----------------------------------------------------------------------------=!
!
      SUBROUTINE PARALLEL_SUM_INTEGER(ARRAY,N,gid)
        IMPLICIT NONE
        INTEGER N, ERR, I, gid
        INTEGER ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart
        INTEGER ARRAY_TMP(N)
        msgsiz_max = INTEGER_MESSAGE_MAX_SIZE 
        IF( n .LT. msgsiz_max ) THEN
          CALL MPI_ALLREDUCE(ARRAY, ARRAY_TMP, N, MPI_INTEGER,  &
     &       MPI_SUM, gid, ERR)
        ELSE
          nblk   = n / msgsiz_max
          blksiz = msgsiz_max
          DO iblk = 1, nblk
            istart = (iblk-1)*msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY(istart), ARRAY_TMP(istart), blksiz, &
     &         MPI_INTEGER, MPI_SUM, gid, ERR)
          END DO
          blksiz = MOD( n, msgsiz_max )
          IF( blksiz .GT. 0 ) THEN
            istart = nblk * msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY(istart), ARRAY_TMP(istart), blksiz, &
     &         MPI_INTEGER, MPI_SUM, gid, ERR)
          END IF
        END IF
        ARRAY = ARRAY_TMP
#endif
      RETURN
      END SUBROUTINE PARALLEL_SUM_INTEGER
!
!=----------------------------------------------------------------------------=!
!
      SUBROUTINE PARALLEL_SUM_REAL_TO(ARRAY_IN,ARRAY_OUT,N, gid)
        IMPLICIT NONE
        INTEGER N,ERR, gid
        REAL*8 ARRAY_IN(N)
        REAL*8 ARRAY_OUT(N)
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart
#if defined __MPI 
        INCLUDE 'mpif.h'
        msgsiz_max = REAL_MESSAGE_MAX_SIZE 

        IF( n .LT. msgsiz_max) THEN
          CALL MPI_ALLREDUCE(ARRAY_IN, ARRAY_OUT, N, MPI_DOUBLE_PRECISION,  &
     &       MPI_SUM, gid, ERR)
        ELSE
          nblk   = n / msgsiz_max
          blksiz = msgsiz_max
          DO iblk = 1, nblk
            istart = (iblk-1)*msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY_IN(istart), ARRAY_OUT(istart), blksiz, &
     &         MPI_DOUBLE_PRECISION, MPI_SUM, gid, ERR)
          END DO
          blksiz = MOD( n, msgsiz_max )
          IF( blksiz .GT. 0 ) THEN
            istart = nblk * msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY_IN(istart), ARRAY_OUT(istart), blksiz, &
     &         MPI_DOUBLE_PRECISION, MPI_SUM, gid, ERR)
          END IF
        END IF
#else
        ARRAY_OUT = ARRAY_IN
#endif
        RETURN
      END SUBROUTINE PARALLEL_SUM_REAL_TO

!=----------------------------------------------------------------------------=!
!
      SUBROUTINE PARALLEL_SUM_REAL_TO_ROOT(ARRAY_IN,ARRAY_OUT,N, root, gid)
        IMPLICIT NONE
        INTEGER N,ERR, gid, root
        REAL*8 ARRAY_IN(N)
        REAL*8 ARRAY_OUT(N)
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart

#if defined __MPI 

        INCLUDE 'mpif.h'
        msgsiz_max = REAL_MESSAGE_MAX_SIZE 

        IF( n .LT. msgsiz_max) THEN
          CALL MPI_REDUCE(ARRAY_IN, ARRAY_OUT, N, MPI_DOUBLE_PRECISION,  &
     &       MPI_SUM, root, gid, ERR)
        ELSE
          nblk   = n / msgsiz_max
          blksiz = msgsiz_max
          DO iblk = 1, nblk
            istart = (iblk-1)*msgsiz_max + 1
            CALL MPI_REDUCE(ARRAY_IN(istart), ARRAY_OUT(istart), blksiz, &
     &         MPI_DOUBLE_PRECISION, MPI_SUM, root, gid, ERR)
          END DO
          blksiz = MOD( n, msgsiz_max )
          IF( blksiz .GT. 0 ) THEN
            istart = nblk * msgsiz_max + 1
            CALL MPI_REDUCE(ARRAY_IN(istart), ARRAY_OUT(istart), blksiz, &
     &         MPI_DOUBLE_PRECISION, MPI_SUM, root, gid, ERR)
          END IF
        END IF

#else

        ARRAY_OUT = ARRAY_IN

#endif

        RETURN
      END SUBROUTINE PARALLEL_SUM_REAL_TO_ROOT

!=----------------------------------------------------------------------------=!

      SUBROUTINE PARALLEL_SUM_COMPLEX_TO(ARRAY_IN,ARRAY_OUT,N, gid)
        IMPLICIT NONE
        INTEGER N,ERR, gid
        COMPLEX*16 ARRAY_IN(N)
        COMPLEX*16 ARRAY_OUT(N)
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart
#if defined __MPI
        INCLUDE 'mpif.h'
        msgsiz_max = COMPLEX_MESSAGE_MAX_SIZE

        IF( n .LT. msgsiz_max ) THEN
          CALL MPI_ALLREDUCE(ARRAY_IN, ARRAY_OUT, N, MPI_DOUBLE_COMPLEX,  &
     &       MPI_SUM, gid, ERR)
        ELSE
          nblk   = n / msgsiz_max
          blksiz = msgsiz_max
          DO iblk = 1, nblk
            istart = (iblk-1)*msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY_IN(istart), ARRAY_OUT(istart), blksiz, &
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, gid, ERR)
          END DO 
          blksiz = MOD( n, msgsiz_max )
          IF( blksiz .GT. 0 ) THEN
            istart = nblk * msgsiz_max + 1
            CALL MPI_ALLREDUCE(ARRAY_IN(istart), ARRAY_OUT(istart), blksiz, &
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, gid, ERR)
          END IF
        END IF
#else   
        ARRAY_OUT = ARRAY_IN
#endif  
        RETURN
      END SUBROUTINE PARALLEL_SUM_COMPLEX_TO

!=----------------------------------------------------------------------------=!

      SUBROUTINE PARALLEL_MAX_INTEGER(ARRAY,N, gid)
        IMPLICIT NONE
        INTEGER N, I, IERR, gid
        INTEGER ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart
        INTEGER IWORK(N)
        msgsiz_max = INTEGER_MESSAGE_MAX_SIZE 
        IF( n .GT. msgsiz_max ) THEN
          WRITE(6,*) ' *** WARNING PARALLEL_MAX_INTEGER: MSGSIZ > MSGSIZ_MAX '
        END IF
        CALL MPI_ALLREDUCE(ARRAY,IWORK,N,MPI_INTEGER, MPI_MAX,gid,IERR)
        ARRAY = IWORK
#endif
        RETURN
      END SUBROUTINE PARALLEL_MAX_INTEGER

!=----------------------------------------------------------------------------=!

      SUBROUTINE PARALLEL_MIN_INTEGER(ARRAY,N, gid)
        IMPLICIT NONE
        INTEGER N, I, IERR, gid
        INTEGER ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart
        INTEGER IWORK(N)
        msgsiz_max = INTEGER_MESSAGE_MAX_SIZE 
        IF( n .GT. msgsiz_max ) THEN
          WRITE(6,*) ' *** WARNING PARALLEL_MIN_INTEGER: MSGSIZ > MSGSIZ_MAX '
        END IF
        CALL MPI_ALLREDUCE(ARRAY,IWORK,N,MPI_INTEGER, MPI_MIN,gid,IERR)
        ARRAY = IWORK
#endif
        RETURN
      END SUBROUTINE PARALLEL_MIN_INTEGER

!=----------------------------------------------------------------------------=!

      SUBROUTINE PARALLEL_MAX_REAL(ARRAY,N, gid)
        IMPLICIT NONE
        INTEGER N, IERR, gid
        REAL*8  ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart
        REAL*8  RWORK(N)
        msgsiz_max = REAL_MESSAGE_MAX_SIZE
        IF( n .GT. msgsiz_max ) THEN
          WRITE(6,*) ' *** WARNING PARALLEL_MAX_REAL: MSGSIZ > MSGSIZ_MAX '
        END IF
        CALL MPI_ALLREDUCE(ARRAY,RWORK,N,MPI_DOUBLE_PRECISION, &
     &       MPI_MAX,gid,IERR)
        ARRAY = RWORK
#endif
        RETURN
      END SUBROUTINE PARALLEL_MAX_REAL

!=----------------------------------------------------------------------------=!

      SUBROUTINE PARALLEL_MIN_REAL(ARRAY,N, gid)
        IMPLICIT NONE
        INTEGER N, IERR, gid
        REAL*8  ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INTEGER msgsiz_max
        INTEGER nblk, blksiz, msgsiz, iblk, istart
        REAL*8  RWORK(N)
        msgsiz_max = REAL_MESSAGE_MAX_SIZE
        IF( n .GT. msgsiz_max ) THEN
          WRITE(6,*) ' *** WARNING PARALLEL_MIN_REAL: MSGSIZ > MSGSIZ_MAX '
        END IF
        CALL MPI_ALLREDUCE(ARRAY,RWORK,N,MPI_DOUBLE_PRECISION, &
     &       MPI_MIN,gid,IERR)
        ARRAY = RWORK
#endif
        RETURN
      END SUBROUTINE PARALLEL_MIN_REAL


#if defined __T3E


      SUBROUTINE SHMEM_SUM_REAL(ARRAY, N)

        IMPLICIT NONE
        INTEGER   :: n, nproc
        REAL*8    :: array(n)
        INTEGER   :: err, ib, nblock, indx, i, block_size
        INTEGER   :: num_pes

#if defined __SHMEM

        include 'reduce.h'
!
        
        nproc = NUM_PES()
        IF ( n <= bufsiz ) THEN
          buffer(1:n) = array(1:n)
          CALL shmem_barrier_all
          CALL SHMEM_REAL8_SUM_TO_ALL(buffer, buffer, n, 0, 0, nproc, work, pSync_sta)
          CALL shmem_barrier_all
          array(1:n) = buffer(1:n)
        ELSE
          nblock = n / bufsiz
          DO ib = 1, nblock
            indx = (ib-1) * bufsiz
            DO i = 1, bufsiz
              buffer(i) = array(indx+i)
            END DO
            CALL shmem_barrier_all
            CALL SHMEM_REAL8_SUM_TO_ALL(buffer, buffer, bufsiz, 0, 0, nproc, work, pSync_sta)
            CALL shmem_barrier_all
            DO i = 1, bufsiz
              array(indx+i) = buffer(i) 
            END DO
          END DO
          block_size = MOD(n, bufsiz)
          IF(block_size > 0 ) THEN
            indx = nblock * bufsiz 
            DO i = 1, block_size
              buffer(i) = array(indx+i)
            END DO
            CALL shmem_barrier_all
            CALL SHMEM_REAL8_SUM_TO_ALL(buffer, buffer, block_size, 0, 0, nproc, work, pSync_sta)
            CALL shmem_barrier_all
            DO i = 1, block_size
              array(indx+i) = buffer(i)
            END DO
          END IF
        END IF
#endif
        RETURN
      END SUBROUTINE SHMEM_SUM_REAL



      SUBROUTINE SHMEM_SUM_REAL_TO(ARRAY, ARRAY_SUM, N)

        IMPLICIT NONE
        INTEGER   :: n, nproc
        REAL*8    :: array(n)
        REAL*8    :: array_sum(n)
        INTEGER   :: err, ib, nblock, indx, i, block_size
        INTEGER   :: num_pes

#if defined __SHMEM

!
!     Include file for SHMEM Library
!
        include 'reduce.h'

        nproc = NUM_PES()
        IF( n .LT. bufsiz ) THEN
          buffer(1:n) = array(1:n)
          CALL shmem_barrier_all
          CALL SHMEM_REAL8_SUM_TO_ALL(buffer, buffer, n, 0, 0, nproc, work, pSync_sta)
          CALL shmem_barrier_all
          array_sum(1:n) = buffer(1:n)
        ELSE
          nblock = n / bufsiz
          DO ib = 1, nblock
            indx = (ib-1) * bufsiz
            DO i = 1, bufsiz
              buffer(i) = array(indx+i)
            END DO
            CALL shmem_barrier_all
            CALL SHMEM_REAL8_SUM_TO_ALL(buffer, buffer, bufsiz, 0, 0, nproc, work, pSync_sta)
            CALL shmem_barrier_all
            DO i = 1, bufsiz
              array_sum(indx+i) = buffer(i) 
            END DO
          END DO
          block_size = n - nblock * bufsiz
          IF(block_size.GT.0) THEN
            indx = nblock * bufsiz 
            DO i = 1, block_size
              buffer(i) = array(indx+i)
            END DO
            CALL shmem_barrier_all
            CALL SHMEM_REAL8_SUM_TO_ALL(buffer, buffer, block_size, 0, 0, nproc, work, pSync_sta)
            CALL shmem_barrier_all
            DO i = 1, block_size
              array_sum(indx+i) = buffer(i)
            END DO
          END IF
        END IF
#endif
        RETURN
      END SUBROUTINE SHMEM_SUM_REAL_TO

#endif
