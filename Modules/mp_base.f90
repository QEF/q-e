!
! Copyright (C) 2002-2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#define __MSGSIZ_MAX 10000


      !
      !  Wrapper for MPI implementations that have problems with large messages
      !
      !  In some Cray XD1 systems the communicaction subsystem
      !  crashes when message exceed a given size, so we need
      !  to break down MPI communications in smaller pieces
      !

!=----------------------------------------------------------------------------=!
!

   SUBROUTINE BCAST_REAL( array, n, root, gid )
        USE kinds, ONLY: DP
        USE parallel_include  
        IMPLICIT NONE
        INTEGER  :: n, root, gid, ierr
        REAL(DP) :: array( n )
#if defined __MPI
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, msgsiz, iblk, istart, i
#if defined __TRACE
        write(*,*) 'BCAST_REAL IN'
#endif
        IF( n <= msgsiz_max ) THEN
           CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, root, gid, ierr )
           IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 1 ', ierr )
        ELSE
           nblk   = n / msgsiz_max
           blksiz = msgsiz_max
           DO iblk = 1, nblk
              istart = (iblk-1)*msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 2 ', ierr )
           END DO
           blksiz = MOD( n, msgsiz_max )
           IF( blksiz > 0 ) THEN
              istart = nblk * msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_DOUBLE_PRECISION, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_real ', ' error in mpi_bcast 3 ', ierr )
           END IF
        END IF
#if defined __TRACE
        write(*,*) 'BCAST_REAL OUT'
#endif
#endif
        RETURN
   END SUBROUTINE BCAST_REAL 


   SUBROUTINE BCAST_INTEGER( array, n, root, gid )
        USE parallel_include  
        IMPLICIT NONE
        INTEGER :: n, root, gid, ierr
        INTEGER :: array( n )
#if defined __MPI
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, msgsiz, iblk, istart, i
#if defined __TRACE
        write(*,*) 'BCAST_INTEGER IN'
#endif
        IF( n <= msgsiz_max ) THEN
           CALL MPI_BCAST( array, n, MPI_INTEGER, root, gid, ierr )
           IF( ierr /= 0 ) CALL errore( ' bcast_integer ', ' error in mpi_bcast 1 ', ierr )
        ELSE
           nblk   = n / msgsiz_max
           blksiz = msgsiz_max
           DO iblk = 1, nblk
              istart = (iblk-1)*msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_INTEGER, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_integer ', ' error in mpi_bcast 2 ', ierr )
           END DO
           blksiz = MOD( n, msgsiz_max )
           IF( blksiz > 0 ) THEN
              istart = nblk * msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_INTEGER, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_integer ', ' error in mpi_bcast 3 ', ierr )
           END IF
        END IF
#if defined __TRACE
        write(*,*) 'BCAST_INTEGER OUT'
#endif
#endif
        RETURN
   END SUBROUTINE BCAST_INTEGER


   SUBROUTINE BCAST_LOGICAL( array, n, root, gid )
        USE parallel_include  
        IMPLICIT NONE
        INTEGER :: n, root, gid, ierr
        LOGICAL :: array( n )
#if defined __MPI
        INTEGER :: msgsiz_max = __MSGSIZ_MAX
        INTEGER :: nblk, blksiz, msgsiz, iblk, istart, i
#if defined __TRACE
        write(*,*) 'BCAST_LOGICAL IN'
#endif
        IF( n <= msgsiz_max ) THEN
           CALL MPI_BCAST( array, n, MPI_LOGICAL, root, gid, ierr )
           IF( ierr /= 0 ) CALL errore( ' bcast_logical ', ' error in mpi_bcast 1 ', ierr )
        ELSE
           nblk   = n / msgsiz_max
           blksiz = msgsiz_max
           DO iblk = 1, nblk
              istart = (iblk-1)*msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_LOGICAL, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_logical ', ' error in mpi_bcast 2 ', ierr )
           END DO
           blksiz = MOD( n, msgsiz_max )
           IF( blksiz > 0 ) THEN
              istart = nblk * msgsiz_max + 1
              CALL MPI_BCAST( array( istart ), blksiz, MPI_LOGICAL, root, gid, ierr )
              IF( ierr /= 0 ) CALL errore( ' bcast_logical ', ' error in mpi_bcast 3 ', ierr )
           END IF
        END IF
#if defined __TRACE
        write(*,*) 'BCAST_LOGICAL OUT'
#endif
#endif
        RETURN
   END SUBROUTINE BCAST_LOGICAL


!
! ... "reduce"-like subroutines
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE kinds, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP), INTENT(INOUT) :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined (__SHMEM) && (defined __ALTIX || defined __ORIGIN)
  INTEGER  :: sym_len
  LOGICAL  :: first
  REAL(DP) :: buff(*), snd_buff(*)
  POINTER     (buff_p, buff), (snd_buff_p, snd_buff)
  COMMON /sym_heap1/ buff_p, snd_buff_p, sym_len, first
#else
  REAL(DP) :: buff(maxb)  
#endif
  !
#if defined (__SHMEM)
  !
  ! ... SHMEM specific 
  !
  INCLUDE 'mpp/shmem.fh'
#if defined (__ALTIX) || defined (__ORIGIN)
  INTEGER    :: pWrkSync(SHMEM_REDUCE_SYNC_SIZE), &
                pWrkData(1024*1024), start
  DATA pWrkSync /SHMEM_REDUCE_SYNC_SIZE*SHMEM_SYNC_VALUE/
  DATA pWrkData / 1048576 * 0 /
#else
  ! T3E ? likely obsolete
  INTEGER :: pWrkSync, pWrkData, start
  COMMON / SH_SYNC / pWrkSync(SHMEM_BARRIER_SYNC_dim)
  COMMON / SH_DATA / pWrkData(1024*1024)
  DATA pWrkData / 1048576 * 0 /
  DATA pWrkSync / SHMEM_BARRIER_SYNC_dim * SHMEM_SYNC_VALUE /
!DIR$ CACHE_ALIGN /SH_SYNC/
!DIR$ CACHE_ALIGN /SH_DATA/
#endif
  !
#endif
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_comm_rank', info )

  !
  IF ( dim <= 0 .OR. nproc <= 1 ) RETURN
  !
  ! ... synchronize processes
  !
  CALL mpi_barrier( comm, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_barrier', info )
  !
  nbuf = dim / maxb
  !
#if defined (__SHMEM)

#if defined (__ALTIX) || defined (__ORIGIN)
  IF (dim .GT. sym_len) THEN
     IF (sym_len .NE. 0) THEN
        CALL shpdeallc( snd_buff_p, info, -1 )
     END IF
     sym_len = dim
     CALL shpalloc( snd_buff_p, 2*sym_len, info, -1 )
  END IF
  IF (first .NE. .TRUE.) THEN
     CALL shpalloc( buff_p, 2*maxb, info, -1 )
     first = .TRUE.
  END IF
  snd_buff(1:dim) = ps(1:dim)
#endif
  !
  start = myid * nproc
  !
#endif
  !
  DO n = 1, nbuf
     !
#if defined (__SHMEM)
     !
#if defined (__ALTIX) || defined (__ORIGIN)
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, snd_buff(1+(n-1)*maxb), maxb, &
                                  start, 0, nproc, pWrkData, pWrkSync )
#else
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+(n-1)*maxb), maxb, &
                                  start, 0, nproc, pWrkData, pWrkSync )
#endif
     !                             
#else
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 1', info )
     END IF
     !                    
#endif
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
#if defined (__SHMEM)
     !
#if defined (__ALTIX) || defined (__ORIGIN)
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, snd_buff(1+nbuf*maxb),          &
     &                            (dim-nbuf*maxb), start, 0, nproc,&
     &                            pWrkData, pWrkSync )
#else
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+nbuf*maxb), (dim-nbuf*maxb), &
                                  start, 0, nproc, pWrkData, pWrkSync )
#endif
     !                             
#else
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real', 'error in mpi_allreduce 2', info )
     END IF
     !
#endif
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_integer( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE kinds, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER,  INTENT(INOUT) :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  INTEGER :: buff(maxb)  
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) RETURN
  !
  ! ... synchronize processes
  !
  CALL mpi_barrier( comm, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer', 'error in mpi_barrier', info )
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_integer

!
! ... "reduce"-like subroutines
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_real_to( dim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE kinds, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: dim
  REAL(DP), INTENT(IN)  :: ps(dim)
  REAL(DP), INTENT(OUT) :: psout(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined (__SHMEM) && (defined __ALTIX || defined __ORIGIN)
  INTEGER  :: sym_len
  LOGICAL  :: first
  REAL(DP) :: buff(*), snd_buff(*)
  POINTER     (buff_p, buff), (snd_buff_p, snd_buff)
  COMMON /sym_heap1/ buff_p, snd_buff_p, sym_len, first
#else
  REAL(DP) :: buff(maxb)  
#endif
  !
#if defined (__SHMEM)
  !
  ! ... SHMEM specific 
  !
  INCLUDE 'mpp/shmem.fh'
#if defined (__ALTIX) || defined (__ORIGIN)
  INTEGER    :: pWrkSync(SHMEM_REDUCE_SYNC_SIZE), &
                pWrkData(1024*1024), start
  DATA pWrkSync /SHMEM_REDUCE_SYNC_SIZE*SHMEM_SYNC_VALUE/
  DATA pWrkData / 1048576 * 0 /
#else
  ! T3E ? likely obsolete
  INTEGER :: pWrkSync, pWrkData, start
  COMMON / SH_SYNC / pWrkSync(SHMEM_BARRIER_SYNC_dim)
  COMMON / SH_DATA / pWrkData(1024*1024)
  DATA pWrkData / 1048576 * 0 /
  DATA pWrkSync / SHMEM_BARRIER_SYNC_dim * SHMEM_SYNC_VALUE /
!DIR$ CACHE_ALIGN /SH_SYNC/
!DIR$ CACHE_ALIGN /SH_DATA/
#endif
  !
#endif
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_comm_rank', info )

  !
  IF ( dim <= 0 .OR. nproc < 1 ) RETURN
  !
  IF ( nproc == 1 ) THEN
     psout = ps
     RETURN
  END IF
  !
  ! ... synchronize processes
  !
  CALL mpi_barrier( comm, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_barrier', info )
  !
  nbuf = dim / maxb
  !
#if defined (__SHMEM)

#if defined (__ALTIX) || defined (__ORIGIN)
  IF (dim .GT. sym_len) THEN
     IF (sym_len .NE. 0) THEN
        CALL shpdeallc( snd_buff_p, info, -1 )
     END IF
     sym_len = dim
     CALL shpalloc( snd_buff_p, 2*sym_len, info, -1 )
  END IF
  IF (first .NE. .TRUE.) THEN
     CALL shpalloc( buff_p, 2*maxb, info, -1 )
     first = .TRUE.
  END IF
  snd_buff(1:dim) = ps(1:dim)
#endif
  !
  start = myid * nproc
  !
#endif
  !
  DO n = 1, nbuf
     !
#if defined (__SHMEM)
     !
#if defined (__ALTIX) || defined (__ORIGIN)
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, snd_buff(1+(n-1)*maxb), maxb, &
                                  start, 0, nproc, pWrkData, pWrkSync )
#else
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+(n-1)*maxb), maxb, &
                                  start, 0, nproc, pWrkData, pWrkSync )
#endif
     !                             
     IF( root < 0 ) THEN
        psout((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        psout((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
#else
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_allreduce 1', info )
     END IF
     !                    
#endif
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
#if defined (__SHMEM)
     !
#if defined (__ALTIX) || defined (__ORIGIN)
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, snd_buff(1+nbuf*maxb),          &
     &                            (dim-nbuf*maxb), start, 0, nproc,&
     &                            pWrkData, pWrkSync )
#else
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+nbuf*maxb), (dim-nbuf*maxb), &
                                  start, 0, nproc, pWrkData, pWrkSync )
#endif
     !                             
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
#else
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_allreduce 2', info )
     END IF
     !
#endif
     !
  END IF
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real_to
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_integer_to( dim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE kinds, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: dim
  INTEGER,  INTENT(IN)  :: ps(dim)
  INTEGER,  INTENT(OUT) :: psout(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  INTEGER :: buff(maxb)  
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_to IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc < 1 ) RETURN
  !
  IF ( nproc == 1 ) THEN
     psout = ps
     RETURN
  END IF
  !
  ! ... synchronize processes
  !
  CALL mpi_barrier( comm, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_barrier', info )
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), psout( 1+(n-1)*maxb ), maxb, MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), psout( 1+(n-1)*maxb ), maxb, MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), (dim-nbuf*maxb), MPI_INTEGER, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_integer_to', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
#if defined __TRACE
  write(*,*) 'reduce_base_integer_to OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_integer_to
!
!
!  Parallel MIN and MAX
!

!----------------------------------------------------------------------------
SUBROUTINE parallel_min_integer( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE kinds, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER,  INTENT(INOUT) :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  INTEGER :: buff(maxb)  
  !
#if defined __TRACE
  write(*,*) 'parallel_min_integer IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_integer', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_integer', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) RETURN
  !
  ! ... synchronize processes
  !
  CALL mpi_barrier( comm, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_integer', 'error in mpi_barrier', info )
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_integer', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
#if defined __TRACE
  write(*,*) 'parallel_min_integer OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE parallel_min_integer

!
!----------------------------------------------------------------------------
SUBROUTINE parallel_max_integer( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE kinds, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER,  INTENT(INOUT) :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  INTEGER :: buff(maxb)  
  !
#if defined __TRACE
  write(*,*) 'parallel_max_integer IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_integer', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_integer', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) RETURN
  !
  ! ... synchronize processes
  !
  CALL mpi_barrier( comm, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_integer', 'error in mpi_barrier', info )
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_INTEGER, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_INTEGER, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_integer', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
#if defined __TRACE
  write(*,*) 'parallel_max_integer OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE parallel_max_integer


!----------------------------------------------------------------------------
SUBROUTINE parallel_min_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE kinds, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP), INTENT(INOUT) :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  REAL(DP) :: buff(maxb)  
  !
#if defined __TRACE
  write(*,*) 'parallel_min_real IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_real', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_real', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) RETURN
  !
  ! ... synchronize processes
  !
  CALL mpi_barrier( comm, info )
  IF( info /= 0 ) CALL errore( 'parallel_min_real', 'error in mpi_barrier', info )
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_real', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MIN, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MIN, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_min_real', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
#if defined __TRACE
  write(*,*) 'parallel_min_real OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE parallel_min_real

!
!----------------------------------------------------------------------------
SUBROUTINE parallel_max_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE kinds, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP), INTENT(INOUT) :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
  REAL(DP) :: buff(maxb)  
  !
#if defined __TRACE
  write(*,*) 'parallel_max_real IN'
#endif
  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_real', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_real', 'error in mpi_comm_rank', info )
  !
  IF ( dim <= 0 .OR. nproc <= 1 ) RETURN
  !
  ! ... synchronize processes
  !
  CALL mpi_barrier( comm, info )
  IF( info /= 0 ) CALL errore( 'parallel_max_real', 'error in mpi_barrier', info )
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_real', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_real', 'error in mpi_allreduce 1', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     ELSE IF( root == myid ) THEN
        ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     END IF
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MAX, root, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_real', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_MAX, comm, info )
        IF( info /= 0 ) CALL errore( 'parallel_max_real', 'error in mpi_allreduce 2', info )
     END IF
     !
     IF( root < 0 ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     ELSE IF( root == myid ) THEN
        ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     END IF
     !
  END IF
  !
#if defined __TRACE
  write(*,*) 'parallel_max_real OUT'
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE parallel_max_real


#if defined (__MPI)  
SUBROUTINE hangup()
  INCLUDE 'mpif.h'
  INTEGER IERR
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  CALL MPI_FINALIZE( ierr )
  STOP 'hangup'
END SUBROUTINE hangup
#endif
