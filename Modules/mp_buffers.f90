!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
!   Carlo Cavazzoni
!   Last update 19 May 2001
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!
    MODULE mp_buffers

! This module is used to implement, when possible, high efficient buffered
! communications among processors.
! In particular two buffers are defined:
!
!   mp_snd_buffer   Send buffer
!
!   mp_rcv_buffer   Receive buffer
!
! together with the buffers the module contains initialization and
! communication functions, that may depend on the particular hardware
!
!------------------------------------------------------------------------------!

      USE kinds, ONLY : dbl
      USE parallel_include
      USE shmem_include

      PRIVATE
      PUBLIC :: mp_sendrecv_buffers, mp_allocate_buffers, mp_deallocate_buffers, &
        mp_barrier_buffers, mp_alltoall_buffers,mp_excng, mp_sum_buffers, mp_report_buffers

      SAVE

#if defined __SHMEM
          pointer (mp_p_snd_buffer,mp_snd_buffer)
          pointer (mp_p_rcv_buffer,mp_rcv_buffer)
          complex (dbl) :: mp_snd_buffer(1)
          complex (dbl) :: mp_rcv_buffer(1)
#else
          integer :: mp_p_snd_buffer = 0
          integer :: mp_p_rcv_buffer = 0
          complex (dbl), pointer :: mp_snd_buffer(:)
          complex (dbl), pointer :: mp_rcv_buffer(:)
#endif
          PUBLIC :: mp_snd_buffer, mp_rcv_buffer, &
                    mp_p_snd_buffer, mp_p_rcv_buffer

          integer :: mp_bufsize
          integer :: mp_high_watermark = 0

#if defined COMPLEX_MESSAGE_MAX_SIZE
          INTEGER, PARAMETER :: mp_bufsize_msgmax = COMPLEX_MESSAGE_MAX_SIZE
#else
          INTEGER, PARAMETER :: mp_bufsize_msgmax = 2**20  ! 1Mb 2^20
#endif

          PUBLIC :: mp_bufsize_msgmax

!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!..mp_allocate_buffers
!..Carlo Cavazzoni
      SUBROUTINE mp_allocate_buffers(bufsize)
        IMPLICIT NONE
          INTEGER, INTENT(IN) :: bufsize
          INTEGER :: ERR

#if defined __SHMEM
          CALL SHPALLOC(mp_p_snd_buffer,2*bufsize, err, 0)
#else
          ALLOCATE(mp_snd_buffer(bufsize), STAT = err)
#endif
          IF(ERR /= 0) CALL errore(' mp_allocate_buffers ', ' allocating mp_snd_buffer ',err)

#if defined __SHMEM
          CALL SHPALLOC(mp_p_rcv_buffer,2*bufsize, err, 0)
#else
          ALLOCATE(mp_rcv_buffer(bufsize), STAT = err)
#endif
          IF(ERR /= 0) CALL errore(' mp_allocate_buffers ', ' allocating mp_rcv_buffer ',err)

          mp_bufsize = bufsize

        RETURN
       END SUBROUTINE mp_allocate_buffers      

!
!------------------------------------------------------------------------------!
!..mp_deallocate_buffers
!..Carlo Cavazzoni
      SUBROUTINE mp_deallocate_buffers
        IMPLICIT NONE
          integer err
#if defined __SHMEM
          CALL shmem_barrier_all
          CALL SHPDEALLC(mp_p_rcv_buffer, err, 0)
#else
          DEALLOCATE(mp_snd_buffer, STAT = err)
#endif
          IF(ERR /= 0) CALL errore(' mp_deallocate_buffers ', ' deallocating mp_rcv_buffer ',err)

#if defined __SHMEM
          CALL SHPDEALLC(mp_p_snd_buffer, err, 0)
#else
          DEALLOCATE(mp_rcv_buffer, STAT = err)
#endif
          IF(ERR /= 0) CALL errore(' mp_deallocate_buffers ', ' deallocating mp_snd_buffer ',err)

        RETURN
      END SUBROUTINE mp_deallocate_buffers
!
!------------------------------------------------------------------------------!
!..mp_barrier_buffers
!..Carlo Cavazzoni
      SUBROUTINE mp_barrier_buffers
        IMPLICIT NONE
#if defined __SHMEM
          call shmem_barrier_all
#else
#endif
        RETURN
      END SUBROUTINE mp_barrier_buffers
!
!------------------------------------------------------------------------------!
!..mp_sum_buffers
!..Carlo Cavazzoni

      SUBROUTINE mp_sum_buffers

      IMPLICIT NONE

      INTEGER ierr, pwrksize

#if defined __PARA

#  if defined __SHMEM

      pointer (p_pWrk,pWrk)
      REAL(dbl)  pWrk(1)
      INTEGER :: nproc, num_pes

      nproc = num_pes()
      pwrksize = MAX(2*mp_bufsize,SHMEM_REDUCE_MIN_WRKDATA_SIZE)
      CALL SHPALLOC(p_pWrk, pwrksize, ierr, 0)
      IF(IERR /= 0) THEN
        CALL errore(' mp_sum_buffers ', ' allocating p_pWrk ',ierr)
      END IF

      call shmem_barrier_all
      CALL SHMEM_REAL8_SUM_TO_ALL(mp_rcv_buffer, mp_snd_buffer, &
           2*mp_bufsize, 0, 0, nproc, pWrk, pSync_sta)

      call shmem_barrier_all
      CALL SHPDEALLC(p_pwrk, ierr, 0)
      IF(IERR /= 0) call errore(' mp_sum_buffers ', &
        ' deallocating p_pWrk ',ierr)

#  elif defined __MPI

      CALL MPI_ALLREDUCE(mp_snd_buffer,mp_rcv_buffer,mp_bufsize, &
           MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
      IF(IERR /= 0) call errore(' mp_sum_buffers ', ' mpi_allreduce ',ierr)

#  endif

#else

      call ZCOPY(mp_bufsize,mp_snd_buffer,1,mp_rcv_buffer,1)

#endif

        mp_high_watermark = MAX( mp_high_watermark, 16 * mp_bufsize )

        return
      END SUBROUTINE mp_sum_buffers




!
!------------------------------------------------------------------------------!
!..mp_alltoall_buffers
!..Carlo Cavazzoni

      SUBROUTINE mp_alltoall_buffers(mp_snd_buffer, mp_rcv_buffer)


      IMPLICIT NONE

      COMPLEX(dbl) :: mp_snd_buffer(:)
      COMPLEX(dbl) :: mp_rcv_buffer(:)
      INTEGER :: ierr, nproc, i
      INTEGER :: msg_size


#if defined __PARA



#  if defined __SHMEM

      integer ip, isour, mpime, nproc
      integer my_pe, num_pes
      call shmem_barrier_all
      mpime = my_pe()
      nproc = num_pes()
      msg_size = mp_bufsize/nproc
      IF( (msg_size + 1) > mp_bufsize_msgmax ) THEN
         CALL errore(' mp_alltoall_buffers ', ' bufsize too large ', msg_size)
      END IF

      do ip =1,nproc
        ISOUR = MOD(MPIME-IP+NPROC,NPROC)
        call shmem_get64(mp_rcv_buffer( 1 + isour*msg_size), &
             mp_snd_buffer(1+mpime*msg_size), msg_size*2, isour)
      end do

#  elif defined __MPI


      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      msg_size = mp_bufsize/nproc
      IF( (msg_size + 1) > mp_bufsize_msgmax ) THEN
         CALL errore(' mp_alltoall_buffers ', ' bufsize too large ', msg_size)
      END IF

      !WRITE(6,*) ' MP_BUFFERS DEBUG ', msg_size
      !WRITE(6,*) ' MP_BUFFERS DEBUG ', mp_snd_buffer(1)
      !WRITE(6,*) ' MP_BUFFERS DEBUG ', mp_snd_buffer(1+msg_size)

      call MPI_ALLTOALL(mp_snd_buffer,msg_size,MPI_DOUBLE_COMPLEX, &
                        mp_rcv_buffer,msg_size,MPI_DOUBLE_COMPLEX, &
                        MPI_COMM_WORLD,IERR)

      !WRITE(6, fmt='(10D8.2)' ) mp_rcv_buffer(1:mp_bufsize)
      !WRITE(6,*) ' MP_BUFFERS DEBUG ', mp_rcv_buffer(1)
      !WRITE(6,*) ' MP_BUFFERS DEBUG ', mp_rcv_buffer(1+msg_size)

      IF(IERR /= 0) call errore(' mp_alltoall_buffers ', ' mpi_alltoall ',ierr)

#  endif

#else

      msg_size = mp_bufsize
      CALL ZCOPY(msg_size, mp_snd_buffer, 1, mp_rcv_buffer, 1)

#endif

        mp_high_watermark = MAX( mp_high_watermark, 16 * msg_size )

        return
      END SUBROUTINE mp_alltoall_buffers

!
!------------------------------------------------------------------------------!
!..mp_sendrecv_buffers
!..Carlo Cavazzoni
      SUBROUTINE mp_sendrecv_buffers(isour, idest, ip)


      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  isour, idest, ip

#if defined __PARA

#  if defined __MPI

      INTEGER :: istatus(MPI_STATUS_SIZE)
      INTEGER :: ierr

#  endif

#  if defined __SHMEM

      call shmem_barrier_all
      call shmem_get64(mp_rcv_buffer, mp_snd_buffer, mp_bufsize*2, isour-1)

#  elif defined __MPI

      CALL MPI_SENDRECV(mp_snd_buffer, mp_bufsize, MPI_DOUBLE_COMPLEX, &
          IDEST-1, ip, mp_rcv_buffer, mp_bufsize, MPI_DOUBLE_COMPLEX, &
          ISOUR-1, ip, MPI_COMM_WORLD, ISTATUS, ierr)
      IF(ierr /= 0) call errore(' mp_sendrecv_buffers ', ' MPI_SENDRECV ', ierr)

#  endif

#else 

      CALL ZCOPY(mp_bufsize, mp_snd_buffer, 1, mp_rcv_buffer, 1)

#endif

        mp_high_watermark = MAX( mp_high_watermark, 16 * mp_bufsize )

        RETURN
      END SUBROUTINE mp_sendrecv_buffers

      SUBROUTINE mp_report_buffers
        WRITE(6, *) 
        WRITE(6, *) '  mp_buffers: high_watermark (bytes): ', mp_high_watermark
        RETURN
      END SUBROUTINE mp_report_buffers

    END MODULE mp_buffers
