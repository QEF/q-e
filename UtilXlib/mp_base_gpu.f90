!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!  Wrapper for MPI implementations that have problems with large messages
!


!  In some MPI implementation the communication subsystem
!  crashes when message exceeds a given size, so we need
!  to break down MPI communications in smaller pieces
!
#define __MSGSIZ_MAX 100000
#define __BCAST_MSGSIZ_MAX 100000

!  Some implementation of MPI (OpenMPI) if it is not well tuned for the given
!  network hardware (InfiniBand) tend to lose performance or get stuck inside
!  collective routines if processors are not well synchronized
!  A barrier fixes the problem
!
#define __USE_BARRIER


!=----------------------------------------------------------------------------=!
!
! These routines allocate buffer spaces used in reduce_base_real_gpu.
! These should be in data_buffer.f90 but need to be here becouse size is
! depends on the __MSGSIZ_MAX definition

#if defined (__CUDA)
   SUBROUTINE allocate_buffers_gpu
       USE data_buffer
       IMPLICIT NONE
       INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
       !    
       ALLOCATE(mp_buff_r_d(maxb), mp_buff_c_d(maxb), mp_buff_i_d(maxb))
       !
   END SUBROUTINE allocate_buffers_gpu
   
   SUBROUTINE deallocate_buffers_gpu
       USE data_buffer
       IMPLICIT NONE
       !
       DEALLOCATE(mp_buff_r_d, mp_buff_c_d, mp_buff_i_d)
       !
   END SUBROUTINE deallocate_buffers_gpu


!=----------------------------------------------------------------------------=!
!

   SUBROUTINE bcast_real_gpu( array_d, n, root, gid )
        USE util_param, ONLY: DP
        USE parallel_include  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n, root, gid
        REAL(DP), DEVICE :: array_d( n )
#if defined __MPI
        INTEGER :: msgsiz_max = __BCAST_MSGSIZ_MAX
        INTEGER :: nblk, blksiz, iblk, istart, ierr

#if defined __TRACE
        write(*,*) 'BCAST_REAL_GPU IN'
#endif
        IF( n <= 0 ) GO TO 1

#if defined __USE_BARRIER
        CALL mp_synchronize( gid )
#endif

        IF( n <= msgsiz_max ) THEN
           CALL MPI_BCAST( array_d, n, MPI_DOUBLE_PRECISION, root, gid, ierr )
           IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 1 ', ierr )
        ELSE
           nblk   = n / msgsiz_max
           blksiz = msgsiz_max
           DO iblk = 1, nblk
              istart = (iblk-1)*msgsiz_max + 1
              CALL MPI_BCAST( array_d( istart ), blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 2 ', ierr )
           END DO
           blksiz = MOD( n, msgsiz_max )
           IF( blksiz > 0 ) THEN
              istart = nblk * msgsiz_max + 1
              CALL MPI_BCAST( array_d( istart ), blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 3 ', ierr )
           END IF
        ENDIF

1       CONTINUE
#if defined __TRACE
        write(*,*) 'BCAST_REAL_GPU OUT'
#endif

#endif

        RETURN
   END SUBROUTINE bcast_real_gpu

   SUBROUTINE bcast_integer_gpu( array_d, n, root, gid )
        USE parallel_include  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n, root, gid
        INTEGER, DEVICE :: array_d( n )
#if defined __MPI
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, iblk, istart, ierr

#if defined __TRACE
        write(*,*) 'BCAST_INTEGER_GPU IN'
#endif

        IF( n <= 0 ) GO TO 1

#if defined __USE_BARRIER
        CALL mp_synchronize( gid )
#endif

        IF( n <= msgsiz_max ) THEN
           CALL MPI_BCAST( array_d, n, MPI_INTEGER, root, gid, ierr )
           IF( ierr /= 0 ) CALL errore( ' bcast_integer_gpu ', ' error in mpi_bcast 1 ', ierr )
        ELSE
           nblk   = n / msgsiz_max
           blksiz = msgsiz_max
           DO iblk = 1, nblk
              istart = (iblk-1)*msgsiz_max + 1
              CALL MPI_BCAST( array_d( istart ), blksiz, MPI_INTEGER, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_integer_gpu ', ' error in mpi_bcast 2 ', ierr )
           END DO
           blksiz = MOD( n, msgsiz_max )
           IF( blksiz > 0 ) THEN
              istart = nblk * msgsiz_max + 1
              CALL MPI_BCAST( array_d( istart ), blksiz, MPI_INTEGER, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_integer_gpu ', ' error in mpi_bcast 3 ', ierr )
           END IF
        END IF
1       CONTINUE
#if defined __TRACE
        write(*,*) 'BCAST_INTEGER_GPU OUT'
#endif
#endif
        RETURN
   END SUBROUTINE bcast_integer_gpu


   SUBROUTINE bcast_logical_gpu( array_d, n, root, gid )
        USE parallel_include  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n, root, gid
        LOGICAL, DEVICE :: array_d( n )
#if defined __MPI
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, iblk, istart, ierr

#if defined __TRACE
        write(*,*) 'BCAST_LOGICAL_GPU IN'
#endif

        IF( n <= 0 ) GO TO 1

#if defined __USE_BARRIER
        CALL mp_synchronize( gid )
#endif

        IF( n <= msgsiz_max ) THEN
           CALL MPI_BCAST( array_d, n, MPI_LOGICAL, root, gid, ierr )
           IF( ierr /= 0 ) CALL errore( ' bcast_logical_gpu ', ' error in mpi_bcast 1 ', ierr )
        ELSE
           nblk   = n / msgsiz_max
           blksiz = msgsiz_max
           DO iblk = 1, nblk
              istart = (iblk-1)*msgsiz_max + 1
              CALL MPI_BCAST( array_d( istart ), blksiz, MPI_LOGICAL, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_logical_gpu ', ' error in mpi_bcast 2 ', ierr )
           END DO
           blksiz = MOD( n, msgsiz_max )
           IF( blksiz > 0 ) THEN
              istart = nblk * msgsiz_max + 1
              CALL MPI_BCAST( array_d( istart ), blksiz, MPI_LOGICAL, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_logical_gpu ', ' error in mpi_bcast 3 ', ierr )
           END IF
        END IF

1       CONTINUE
#if defined __TRACE
        write(*,*) 'BCAST_LOGICAL_GPU OUT'
#endif
#endif
        RETURN
   END SUBROUTINE bcast_logical_gpu


!
! ... "reduce"-like subroutines
!
#if defined (__USE_INPLACE_MPI)
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_real_gpu( dim, ps_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY: DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim     ! size of the array
  REAL(DP), DEVICE        :: ps_d(dim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_gpu IN'
#endif
  !
  IF ( dim <= 0 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  IF( root >= 0 ) THEN
     CALL MPI_REDUCE( MPI_IN_PLACE, ps_d, dim, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
     IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in mpi_reduce 1', info )
  ELSE
     CALL MPI_ALLREDUCE( MPI_IN_PLACE, ps_d, dim, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
     IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in mpi_allreduce 1', info )
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_gpu OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real_gpu
!
#else
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_real_gpu( dim, ps_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer,    ONLY : mp_buff_r_d
  USE parallel_include  
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim     ! size of the array
  REAL(DP), DEVICE        :: ps_d(dim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_gpu IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+(n-1)*maxb), mp_buff_r_d, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+(n-1)*maxb), mp_buff_r_d, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in mpi_allreduce 1', info )
     END IF
     !                    
     IF( root < 0 ) THEN
        info = cudaMemcpy( ps_d((1+(n-1)*maxb)), mp_buff_r_d(1), maxb, cudaMemcpyDeviceToDevice )
        !ps((1+(n-1)*maxb):(n*maxb)) = mp_buff_r(1:maxb)
        IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        info = cudaMemcpy( ps_d((1+(n-1)*maxb)), mp_buff_r_d(1), maxb, cudaMemcpyDeviceToDevice )
        !ps((1+(n-1)*maxb):(n*maxb)) = mp_buff_r(1:maxb)
        IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+nbuf*maxb), mp_buff_r_d, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+nbuf*maxb), mp_buff_r_d, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), mp_buff_r_d(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = mp_buff_r(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), mp_buff_r_d(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = mp_buff_r(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'reduce_base_real_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_gpu OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real_gpu
!
#endif
!
!
#if defined (__USE_INPLACE_MPI)
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_integer_gpu( dim, ps_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY: DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim     ! size of the array
  INTEGER,  DEVICE        :: ps_d(dim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_gpu IN'
#endif
  !
  IF ( dim <= 0 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  IF( root >= 0 ) THEN
     CALL MPI_REDUCE( MPI_IN_PLACE, ps_d, dim, MPI_INTEGER, MPI_SUM, root, comm, info )
     IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in mpi_reduce 1', info )
  ELSE
     CALL MPI_ALLREDUCE( MPI_IN_PLACE, ps_d, dim, MPI_INTEGER, MPI_SUM, comm, info )
     IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in mpi_allreduce 1', info )
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_gpu OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_integer_gpu
!
#else
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_integer_gpu( dim, ps_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer,    ONLY : mp_buff_i_d
  USE parallel_include  
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim     ! size of the array
  INTEGER,  DEVICE        :: ps_d(dim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_gpu IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1  ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+(n-1)*maxb), mp_buff_i_d, maxb, MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+(n-1)*maxb), mp_buff_i_d, maxb, MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in mpi_allreduce 1', info )
     END IF
     !                    
     IF( root < 0 ) THEN
        info = cudaMemcpy( ps_d((1+(n-1)*maxb)), mp_buff_i_d(1), maxb, cudaMemcpyDeviceToDevice )
        !ps((1+(n-1)*maxb):(n*maxb)) = mp_buff_r(1:maxb)
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        info = cudaMemcpy( ps_d((1+(n-1)*maxb)), mp_buff_i_d(1), maxb, cudaMemcpyDeviceToDevice )
        !ps((1+(n-1)*maxb):(n*maxb)) = mp_buff_r(1:maxb)
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+nbuf*maxb), mp_buff_i_d, (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+nbuf*maxb), mp_buff_i_d, (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), mp_buff_i_d(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = mp_buff_r(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), mp_buff_i_d(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = mp_buff_r(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_gpu OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_integer_gpu
!
#endif
!
! ... "reduce"-like subroutines
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_real_to_gpu( dim, ps_d, psout_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors,
  ! ... and store the results in variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE util_param, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)          :: dim
  REAL(DP), DEVICE, INTENT(IN)  :: ps_d(dim)
  REAL(DP), DEVICE              :: psout_d(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                   ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to_gpu IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to_gpu', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to_gpu', 'error in mpi_comm_rank', info )
  !
  IF ( dim > 0 .AND. nproc <= 1 ) THEN
     psout_d = ps_d
  END IF
  IF( dim <= 0 .OR. nproc <= 1 ) GO TO 1 ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+(n-1)*maxb), psout_d(1+(n-1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to_gpu', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+(n-1)*maxb), psout_d(1+(n-1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to_gpu', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+nbuf*maxb), psout_d(1+nbuf*maxb), (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to_gpu', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+nbuf*maxb), psout_d(1+nbuf*maxb), (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to_gpu', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to_gpu OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real_to_gpu
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_integer_to_gpu( dim, ps_d, psout_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed integer variable ps(dim) over the processors, and
  ! ... saves the result on the output variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE util_param, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)          :: dim
  INTEGER, DEVICE, INTENT(IN)  :: ps_d(dim)
  INTEGER, DEVICE              :: psout_d(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_to_gpu IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_to_gpu', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_to_gpu', 'error in mpi_comm_rank', info )
  !
  IF ( dim > 0 .AND. nproc <= 1 ) THEN
     psout_d = ps_d
  END IF
  IF( dim <= 0 .OR. nproc <= 1 ) GO TO 1 ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+(n-1)*maxb), psout_d( 1+(n-1)*maxb ), maxb, MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to_gpu', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+(n-1)*maxb), psout_d( 1+(n-1)*maxb ), maxb, MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to_gpu', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+nbuf*maxb), psout_d(1+nbuf*maxb), (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to_gpu', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+nbuf*maxb), psout_d(1+nbuf*maxb), (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to_gpu', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_to_gpu OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_integer_to_gpu
!
!
!  Parallel MIN and MAX
!

!----------------------------------------------------------------------------
SUBROUTINE parallel_min_integer_gpu( dim, ps_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the minimum of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_i_d
  USE parallel_include  
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER,  DEVICE        :: ps_d(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'parallel_min_integer_gpu IN'
#endif
  !
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        !ps_d((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
        info = cudaMemcpy( ps_d(1+(n-1)*maxb), buff(1), maxb, cudaMemcpyDeviceToDevice )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        !ps_d((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
        info = cudaMemcpy( ps_d((1+(n-1)*maxb)), buff(1), maxb, cudaMemcpyDeviceToDevice )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), buff(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), buff(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'parallel_min_integer_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'parallel_min_integer_gpu OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE parallel_min_integer_gpu

!
!----------------------------------------------------------------------------
SUBROUTINE parallel_max_integer_gpu( dim, ps_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the maximum of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_i_d
  USE parallel_include  
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER,  DEVICE        :: ps_d(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'parallel_max_integer_gpu IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        !ps_d((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
        info = cudaMemcpy( ps_d(1+(n-1)*maxb), buff(1), maxb, cudaMemcpyDeviceToDevice )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        !ps_d((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
        info = cudaMemcpy( ps_d((1+(n-1)*maxb)), buff(1), maxb, cudaMemcpyDeviceToDevice )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), buff(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), buff(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'parallel_max_integer_gpu OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE parallel_max_integer_gpu


!----------------------------------------------------------------------------
SUBROUTINE parallel_min_real_gpu( dim, ps_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the minimum of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_r_d
  USE parallel_include  
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP), DEVICE        :: ps_d(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'parallel_min_real_gpu IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_real_gpu', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_real_gpu', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_real_gpu', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_real_gpu', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        !ps_d((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
        info = cudaMemcpy( ps_d(1+(n-1)*maxb), buff(1), maxb, cudaMemcpyDeviceToDevice )
     ELSE IF( root == myid ) THEN
        !ps_d((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
        info = cudaMemcpy( ps_d((1+(n-1)*maxb)), buff(1), maxb, cudaMemcpyDeviceToDevice )
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_real_gpu', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_real_gpu', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), buff(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), buff(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'parallel_max_integer_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'parallel_min_real_gpu OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE parallel_min_real_gpu

!
!----------------------------------------------------------------------------
SUBROUTINE parallel_max_real_gpu( dim, ps_d, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the maximum of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_r_d
  USE parallel_include
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP), DEVICE        :: ps_d(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'parallel_max_real_gpu IN'
#endif
  !
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) GO TO 1
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_max, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_max, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        !ps_d((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
        info = cudaMemcpy( ps_d(1+(n-1)*maxb), buff(1), maxb, cudaMemcpyDeviceToDevice )
        IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        !ps_d((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
        info = cudaMemcpy( ps_d((1+(n-1)*maxb)), buff(1), maxb, cudaMemcpyDeviceToDevice )
        IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps_d(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_max, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps_d(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_max, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), buff(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in cudaMemcpy ', info )
     ELSE IF( root == myid ) THEN
        info = cudaMemcpy( ps_d((1+nbuf*maxb)), buff(1), dim-nbuf*maxb, cudaMemcpyDeviceToDevice )
        !ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
        IF( info /= 0 ) CALL errore( 'parallel_max_real_gpu', 'error in cudaMemcpy ', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'parallel_max_real_gpu OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE parallel_max_real_gpu
#endif
