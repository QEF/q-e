!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined __TEST_MYMATMUL
PROGRAM matmul
  IMPLICIT NONE
  REAL*8, ALLOCATABLE :: ax(:,:), bx(:,:), cx(:,:), dx(:,:), ex(:,:), fx(:,:)
  REAL*8, ALLOCATABLE :: a(:,:), b(:,:), c(:,:), d(:,:), e(:,:), f(:,:)
  INTEGER :: n, nloc, mpime, nproc, root, group
  INTEGER :: i, j, k, iloc, jloc
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic

  CALL PARALLEL_STARTUP(NPROC,MPIME,ROOT,GROUP)

  n = 5
  nloc = ldim_cyclic(n, nproc, mpime)
  ALLOCATE( a(n,n), b(n,n), c(n,n), d(n,n), e(n,n), f(n,n) )
  ALLOCATE( ax(nloc,n), bx(nloc,n), cx(nloc,n), dx(n,nloc), ex(nloc,n), fx(n,nloc) )

  CALL RANDOM_NUMBER( a )
  CALL RANDOM_NUMBER( b )

  do j = 1, n
    do i = 1, n
      IF( mpime == owner_cyclic(i, n, nproc) ) THEN
        iloc = lind_cyclic(i, n, nproc, mpime)
        ax( iloc, j ) = a( i, j ) 
        bx( iloc, j ) = b( i, j ) 
      END IF
      IF( mpime == owner_cyclic(j, n, nproc) ) THEN
        jloc = lind_cyclic(j, n, nproc, mpime)
        dx( i, jloc ) = a( i, j )
        fx( i, jloc ) = b( i, j )
      END IF
    end do
  end do

  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0d0
      do k = 1, n
        !c(i,j) = c(i,j) + a(i,k) * b(k,j)
        !c(i,j) = c(i,j) + a(i,k) * b(j,k)
        !c(i,j) = c(i,j) + a(k,i) * b(k,j)
        c(i,j) = c(i,j) + a(k,i) * b(j,k)
      end do
    end do
  end do

  IF( mpime == root ) THEN
    do i = 1, n
      write(6,100) (c(i,j),j=1,n), (a(i,j),j=1,n), (b(i,j),j=1,n)
    end do
100 FORMAT( 5F7.4, ' |', 5F7.4, ' |', 5F7.4 )
  END IF

  !CALL mymatmul(ax, nloc, 'T', 'R', bx, nloc, 'T', 'R', cx, nloc, 'R', n, mpime, nproc)
  !CALL mymatmul(ax, nloc, 'N', 'R', bx, nloc, 'T', 'R', cx, nloc, 'R', n, mpime, nproc)
  !CALL mymatmul(ax, nloc, 'T', 'R', bx, nloc, 'N', 'R', cx, nloc, 'R', n, mpime, nproc)
  CALL mymatmul(ax, nloc, 'T', 'R', bx, nloc, 'T', 'R', cx, nloc, 'R', n, mpime, nproc)

!  CALL mymatmul(ax, nloc, 'N', 'R', bx, nloc, 'N', 'R', cx, nloc, 'R', n, mpime, nproc)
!  CALL mymatmul(ax, nloc, 'N', 'R', fx, n, 'N', 'C', cx, nloc, 'R', n, mpime, nproc)
!  CALL mymatmul(dx, n, 'N', 'C', bx, nloc, 'N', 'R', cx, nloc, 'R', n, mpime, nproc)
!  CALL mymatmul(dx, n, 'N', 'C', fx, n, 'N', 'C', cx, nloc, 'R', n, mpime, nproc)

  CALL mytrasp_dist(cx, nloc, 'R', dx, n, 'C', n, mpime, nproc)
  CALL mytrasp_dati(dx, n, 'C', ex, nloc, 'R', n, mpime, nproc)
  CALL mytrasp_dati(dx, n, 'C', fx, n, 'C', n, mpime, nproc)

  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0d0
      b(i,j) = 0.0d0
      a(i,j) = 0.0d0
      d(i,j) = 0.0d0
      e(i,j) = 0.0d0
      f(i,j) = 0.0d0
      IF( mpime == owner_cyclic(i, n, nproc) ) THEN
        iloc = lind_cyclic(i, n, nproc, mpime)
        c(i, j) = cx( iloc, j ) 
        a(i, j) = ax( iloc, j ) 
        b(i, j) = bx( iloc, j ) 
        e(i, j) = ex( iloc, j ) 
      END IF
      IF( mpime == owner_cyclic(j, n, nproc) ) THEN
        jloc = lind_cyclic(j, n, nproc, mpime)
        d(i, j) = dx( i, jloc ) 
        f(i, j) = fx( i, jloc ) 
      END IF
    END DO
  END DO

  CALL PARALLEL_SUM_REAL( c, SIZE(c) )
  CALL PARALLEL_SUM_REAL( a, SIZE(a) )
  CALL PARALLEL_SUM_REAL( b, SIZE(b) )
  CALL PARALLEL_SUM_REAL( d, SIZE(d) )
  CALL PARALLEL_SUM_REAL( e, SIZE(e) )
  CALL PARALLEL_SUM_REAL( f, SIZE(f) )

  IF( mpime == root ) THEN
    write(6,*) &
    '--------------------------------------------------------------------------------------------------' 
    do i = 1, n
      write(6,110) (c(i,j),j=1,n), (a(i,j),j=1,n), (b(i,j),j=1,n)
    end do
110 FORMAT( 5F7.4, ' |', 5F7.4, ' |', 5F7.4 )
  END IF
  IF( mpime == root ) THEN
    write(6,*) &
    '--------------------------------------------------------------------------------------------------' 
    do i = 1, n
      write(6,120) (d(i,j),j=1,n), (e(i,j),j=1,n), (f(i,j),j=1,n)
    end do
120 FORMAT( 5F7.4, ' |', 5F7.4, ' |', 5F7.4 )
  END IF

  DEALLOCATE( a, b, c, d, e, f )
  DEALLOCATE( ax, bx, cx, dx, ex, fx )

  CALL PARALLEL_HANGUP
END PROGRAM

#endif


SUBROUTINE mymatmul(ax, lda, tac, dista, bx, ldb, tbc, distb, cx, ldc, distc, n, mpime, nproc)
  IMPLICIT NONE
  INTEGER :: lda, ldc, ldb, n, mpime, nproc
  CHARACTER(LEN=1) :: tac, tbc
  CHARACTER(LEN=1) :: dista, distb, distc
  REAL*8 :: ax(lda,*), bx(ldb,*), cx(ldc,*)
  REAL*8, ALLOCATABLE :: bs(:,:), br(:,:)
  REAL*8, ALLOCATABLE :: as(:,:), ar(:,:)
  REAL*8, ALLOCATABLE :: cs(:,:), cr(:,:), cc(:,:)
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic
  INTEGER :: ip, i, ii, j, jj, k, kk, nloc_ip, nloc
  INTEGER :: nloc_src, nloc_dst, idest, isour
  LOGICAL :: ta, tb 
  CHARACTER(LEN=1) :: da, db, dc

  ta = .FALSE.
  tb = .FALSE.
  da = 'R'
  db = 'R'
  dc = 'R'
  IF( tac == 't' .OR. tac == 'T') THEN
    ta = .TRUE.
  END IF
  IF( tbc == 't' .OR. tbc == 'T') THEN
    tb = .TRUE.
  END IF
  IF( dista == 'c' .OR. dista == 'C' ) THEN
    da = 'C'
  END IF
  IF( distb == 'c' .OR. distb == 'C' ) THEN
    db = 'C'
  END IF
  IF( distc == 'c' .OR. distc == 'C' ) THEN
    dc = 'C'
  END IF
    
  nloc = ldim_cyclic(n, nproc, mpime)

  IF( .NOT. ta .AND. .NOT. tb ) THEN

    IF( dc == 'R' .AND. da == 'R' .AND. db == 'R' ) THEN

      CALL mymatmul_rrrnn(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)

    ELSE IF( dc == 'R' .AND. da == 'R' .AND. db == 'C' ) THEN

      CALL mymatmul_rrcnn(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)

    ELSE IF( dc == 'R' .AND. da == 'C' .AND. db == 'R' ) THEN

      ALLOCATE( ar(nloc, n) )
      CALL mytrasp_dist(ax, lda, 'C', ar, nloc, 'R', n, mpime, nproc)
      CALL mymatmul_rrrnn(ar, nloc, bx, ldb, cx, ldc, n, mpime, nproc)
      DEALLOCATE( ar )

    ELSE IF( dc == 'R' .AND. da == 'C' .AND. db == 'C' ) THEN

      ALLOCATE( ar(nloc, n) )
      CALL mytrasp_dist(ax, lda, 'C', ar, nloc, 'R', n, mpime, nproc)
      CALL mymatmul_rrcnn(ar, nloc, bx, ldb, cx, ldc, n, mpime, nproc)
      DEALLOCATE( ar )

    ELSE IF( dc == 'C' .AND. da == 'R' .AND. db == 'R' ) THEN 

      ALLOCATE( cc(nloc, n) )
      CALL mymatmul_rrrnn(ax, lda, bx, ldb, cc, nloc, n, mpime, nproc)
      CALL mytrasp_dist(cc, nloc, 'R', cx, ldc, 'C', n, mpime, nproc)
      DEALLOCATE( cc )

    ELSE IF( dc == 'C' .AND. da == 'R' .AND. db == 'C' ) THEN

      ALLOCATE( cc(nloc, n) )
      CALL mymatmul_rrcnn(ax, lda, bx, ldb, cc, nloc, n, mpime, nproc)
      CALL mytrasp_dist(cc, nloc, 'R', cx, ldc, 'C', n, mpime, nproc)
      DEALLOCATE( cc )

    ELSE IF( dc == 'C' .AND. da == 'C' .AND. db == 'R' ) THEN

      ALLOCATE( ar(nloc, n) )
      ALLOCATE( cc(nloc, n) )
      CALL mytrasp_dist(ax, lda, 'C', ar, nloc, 'R', n, mpime, nproc)
      CALL mymatmul_rrrnn(ar, nloc, bx, ldb, cc, nloc, n, mpime, nproc)
      CALL mytrasp_dist(cc, nloc, 'R', cx, ldc, 'C', n, mpime, nproc)
      DEALLOCATE( ar )
      DEALLOCATE( cc )

    ELSE IF( dc == 'C' .AND. da == 'C' .AND. db == 'C' ) THEN

      ALLOCATE( ar(nloc, n) )
      ALLOCATE( cc(nloc, n) )
      CALL mytrasp_dist(ax, lda, 'C', ar, nloc, 'R', n, mpime, nproc)
      CALL mymatmul_rrcnn(ar, nloc, bx, ldb, cc, nloc, n, mpime, nproc)
      CALL mytrasp_dist(cc, nloc, 'R', cx, ldc, 'C', n, mpime, nproc)
      DEALLOCATE( ar )
      DEALLOCATE( cc )

    ELSE

          WRITE(6,*) ' ** ERROR (1) in mymatmul ** '
          STOP

    END IF

  ELSE IF(        ta  .AND. (.NOT. tb) ) THEN

    IF( dc == 'R' .AND. da == 'R' .AND. db == 'R' ) THEN

      CALL mymatmul_rrrtn(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)

    ELSE IF( dc == 'C' .AND. da == 'R' .AND. db == 'R' ) THEN

      ALLOCATE( cc(nloc, n) )
      CALL mymatmul_rrrtn(ax, lda, bx, ldb, cc, nloc, n, mpime, nproc)
      CALL mytrasp_dist(cc, nloc, 'R', cx, ldc, 'C', n, mpime, nproc)
      DEALLOCATE( cc )

    ELSE

          WRITE(6,*) ' ** ERROR (2) in mymatmul ** '
          STOP

    END IF

  ELSE IF( (.NOT. ta) .AND.        tb  ) THEN

    IF( dc == 'R' .AND. da == 'R' .AND. db == 'R' ) THEN

      CALL mymatmul_rrrnt(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)

    ELSE IF( dc == 'C' .AND. da == 'R' .AND. db == 'R' ) THEN

      ALLOCATE( cc(nloc, n) )
      CALL mymatmul_rrrnt(ax, lda, bx, ldb, cc, nloc, n, mpime, nproc)
      CALL mytrasp_dist(cc, nloc, 'R', cx, ldc, 'C', n, mpime, nproc)
      DEALLOCATE( cc )

    ELSE

          WRITE(6,*) ' ** ERROR (3) in mymatmul ** '
          STOP

    END IF

  ELSE IF(        ta  .AND.        tb  ) THEN

    IF( dc == 'R' .AND. da == 'R' .AND. db == 'R' ) THEN

      CALL mymatmul_rrrtt(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)

    ELSE IF( dc == 'C' .AND. da == 'R' .AND. db == 'R' ) THEN

      ALLOCATE( cc(nloc, n) )
      CALL mymatmul_rrrtt(ax, lda, bx, ldb, cc, nloc, n, mpime, nproc)
      CALL mytrasp_dist(cc, nloc, 'R', cx, ldc, 'C', n, mpime, nproc)
      DEALLOCATE( cc )

    ELSE

          WRITE(6,*) ' ** ERROR (4) in mymatmul ** '
          STOP

    END IF

  ELSE

          WRITE(6,*) ' ** ERROR (5) in mymatmul ** '
          STOP

  END IF

  RETURN
END SUBROUTINE


SUBROUTINE mymatmul_rrrnn(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)
  IMPLICIT NONE
#if defined __MPI
  INCLUDE 'mpif.h'
  INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
#endif
  INTEGER :: lda, ldc, ldb, n, mpime, nproc
  REAL*8 :: ax(lda,*), bx(ldb,*), cx(ldc,*)
  REAL*8, ALLOCATABLE :: bs(:,:), br(:,:)
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic
  INTEGER :: ip, i, ii, j, jj, k, kk, nloc_ip, nloc
  INTEGER :: nloc_src, nloc_dst, idest, isour
  nloc = ldim_cyclic(n, nproc, mpime)
  DO j = 1, n
    DO i = 1, nloc
      cx(i,j) = 0.0d0
    END DO
  END DO

#if defined __MPI

  DO ip = 0, nproc - 1
    isour = MOD(mpime-ip+nproc, nproc)
    idest = MOD(mpime+ip      , nproc)
    nloc_src = ldim_cyclic(n, nproc, isour)
    ALLOCATE( br( nloc_src, n ) )
    IF( nloc /= ldb ) THEN
      ALLOCATE( bs( nloc, n ) )
      DO j = 1, n
        DO i = 1, nloc
          bs(i,j) = bx(i,j)
        END DO
      END DO
      CALL MPI_SENDRECV(bs, SIZE(bs), MPI_DOUBLE_PRECISION,  &
         IDEST, ip,  br, SIZE(br), MPI_DOUBLE_PRECISION, &
         ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
      DEALLOCATE( bs )
    ELSE
      CALL MPI_SENDRECV(bx, nloc*n, MPI_DOUBLE_PRECISION,  &
         IDEST, ip,  br, SIZE(br), MPI_DOUBLE_PRECISION, &
         ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
    END IF
    IF(ierr .NE. 0) THEN
       WRITE(6,*) ' ** ERROR in sendrecv ** '
       STOP
    END IF

    DO j = 1, n
      kk = isour+1
      DO k = 1, nloc_src 
        DO i = 1, nloc
          cx(i,j) = cx(i,j) + ax(i,kk) * br(k,j)
        END DO
        kk = kk + nproc
      END DO
    END DO
    DEALLOCATE( br )
  END DO

#else

  DO j = 1, n
    DO k = 1, n
      DO i = 1, n
        cx(i,j) = cx(i,j) + ax(i,k) * bx(k,j)
      END DO
    END DO
  END DO

#endif

  RETURN
END SUBROUTINE

SUBROUTINE mymatmul_rrcnn(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)
  IMPLICIT NONE
#if defined __MPI
  INCLUDE 'mpif.h'
  INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
#endif
  INTEGER :: lda, ldc, ldb, n, mpime, nproc
  REAL*8 :: ax(lda,*), bx(ldb,*), cx(ldc,*)
  REAL*8, ALLOCATABLE :: bs(:,:), br(:,:)
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic
  INTEGER :: ip, i, ii, j, jj, k, kk, nloc_ip, nloc
  INTEGER :: nloc_src, nloc_dst, idest, isour
  nloc = ldim_cyclic(n, nproc, mpime)
  DO j = 1, n
    DO i = 1, nloc
      cx(i,j) = 0.0d0
    END DO
  END DO
  DO ip = 0, nproc - 1
    isour = MOD(mpime-ip+nproc, nproc)
    idest = MOD(mpime+ip      , nproc)
    nloc_src = ldim_cyclic(n, nproc, isour)
    ALLOCATE( bs( n, nloc) )
    ALLOCATE( br( n, nloc_src ) )
    DO j = 1, nloc
      DO i = 1, n
        bs(i,j) = bx(i,j)
      END DO
    END DO
    !CALL sendrecv_real(bs, SIZE(bs), idest, br, SIZE(br), isour, ip)
#if defined __MPI
    CALL MPI_SENDRECV(bs, SIZE(bs), MPI_DOUBLE_PRECISION,  &
         IDEST, ip,  br, SIZE(br), MPI_DOUBLE_PRECISION, &
         ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
    IF(ierr .NE. 0) THEN
       WRITE(6,*) ' ** ERROR in sendrecv ** '
       STOP
    END IF
#else
    br = bs
#endif
    jj = isour+1
    DO j = 1, nloc_src
      DO k = 1, n
        DO i = 1, nloc
          cx(i,jj) = cx(i,jj) + ax(i,k) * br(k,j)
        END DO
      END DO
      jj = jj + nproc
    END DO
    DEALLOCATE( bs )
    DEALLOCATE( br )
  END DO
  RETURN
END SUBROUTINE


SUBROUTINE mymatmul_rrrtn(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)
  IMPLICIT NONE
#if defined __MPI
  INCLUDE 'mpif.h'
  INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
#endif
  INTEGER :: lda, ldc, ldb, n, mpime, nproc
  REAL*8 :: ax(lda,*), bx(ldb,*), cx(ldc,*)
  REAL*8, ALLOCATABLE :: bs(:,:), br(:,:)
  REAL*8, ALLOCATABLE :: as(:,:), ar(:,:)
  REAL*8, ALLOCATABLE :: cs(:,:), cr(:,:), cc(:,:)
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic
  INTEGER :: ip, i, ii, j, jj, k, kk, nloc_ip, nloc
  INTEGER :: nloc_src, nloc_dst, idest, isour
  nloc = ldim_cyclic(n, nproc, mpime)

#if defined __MPI

  DO ip = 0, nproc - 1
    nloc_ip = ldim_cyclic(n, nproc, ip)
    ALLOCATE( cs( nloc_ip, n ) )
    DO j = 1, n
      ii = ip+1
      DO i = 1, nloc_ip
        cs(i,j) = 0.0d0
        DO k = 1, nloc
          cs(i,j) = cs(i,j) + ax(k, ii) * bx(k, j)
        END DO
        ii = ii + nproc
      END DO
    END DO
    IF( ldc /= nloc_ip ) THEN
      ALLOCATE( cr( nloc_ip, n ) )
      CALL MPI_REDUCE(cs, cr, SIZE(cs), MPI_DOUBLE_PRECISION, MPI_SUM, ip, MPI_COMM_WORLD, ierr)
      IF( mpime == ip ) THEN
        DO j = 1, n
          DO i = 1, nloc_ip
            cx(i,j) = cr(i,j)
          END DO
        END DO
      END IF
      DEALLOCATE( cr )
    ELSE
      CALL MPI_REDUCE(cs, cx, n*nloc_ip, MPI_DOUBLE_PRECISION, MPI_SUM, ip, MPI_COMM_WORLD, ierr)
    END IF
    DEALLOCATE( cs )
  END DO

#else

  DO j = 1, n
    DO i = 1, n
      cx(i,j) = 0.0d0
      DO k = 1, n
        cx(i,j) = cx(i,j) + ax(k, i) * bx(k, j)
      END DO
    END DO
  END DO

#endif

  RETURN
END SUBROUTINE


SUBROUTINE mymatmul_rrrnt(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)
  IMPLICIT NONE
#if defined __MPI
  INCLUDE 'mpif.h'
  INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
#endif
  INTEGER :: lda, ldc, ldb, n, mpime, nproc
  REAL*8 :: ax(lda,*), bx(ldb,*), cx(ldc,*)
  REAL*8, ALLOCATABLE :: bs(:,:), br(:,:)
  REAL*8, ALLOCATABLE :: as(:,:), ar(:,:)
  REAL*8, ALLOCATABLE :: cs(:,:), cr(:,:), cc(:,:)
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic
  INTEGER :: ip, i, ii, j, jj, k, kk, nloc_ip, nloc
  INTEGER :: nloc_src, nloc_dst, idest, isour
  nloc = ldim_cyclic(n, nproc, mpime)
  DO j = 1, n
    DO i = 1, nloc
      cx(i,j) = 0.0d0
    END DO
  END DO

#if defined __MPI

  DO ip = 0, nproc - 1
    isour = MOD(mpime-ip+nproc, nproc)
    idest = MOD(mpime+ip      , nproc)
    nloc_src = ldim_cyclic(n, nproc, isour)
    ALLOCATE( br( nloc_src, n ) )

    IF( nloc /= ldb ) THEN
      ALLOCATE( bs( nloc, n ) )
      DO j = 1, n
        DO i = 1, nloc
          bs(i,j) = bx(i,j)
        END DO
      END DO
      CALL MPI_SENDRECV(bs, n*nloc, MPI_DOUBLE_PRECISION,  &
         IDEST, ip,  br, SIZE(br), MPI_DOUBLE_PRECISION, &
         ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
      DEALLOCATE( bs )
    ELSE
      CALL MPI_SENDRECV(bx, n*nloc, MPI_DOUBLE_PRECISION,  &
         IDEST, ip,  br, SIZE(br), MPI_DOUBLE_PRECISION, &
         ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
    ENDIF
    IF(ierr .NE. 0) THEN
       WRITE(6,*) ' ** ERROR in sendrecv ** '
       STOP
    END IF
    jj = isour+1
    DO j = 1, nloc_src
      kk = isour+1
      DO k = 1, n
        DO i = 1, nloc
          cx(i,jj) = cx(i,jj) + ax(i,k) * br(j,k)
        END DO
      END DO
      jj = jj + nproc
    END DO
    DEALLOCATE( br )
  END DO

#else

    DO j = 1, n
      DO k = 1, n
        DO i = 1, n
          cx(i,j) = cx(i,j) + ax(i,k) * bx(j,k)
        END DO
      END DO
    END DO

#endif

  RETURN
END SUBROUTINE



SUBROUTINE mymatmul_rrrtt(ax, lda, bx, ldb, cx, ldc, n, mpime, nproc)
  IMPLICIT NONE
#if defined __MPI
  INCLUDE 'mpif.h'
  INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
#endif
  INTEGER :: lda, ldc, ldb, n, mpime, nproc
  REAL*8 :: ax(lda,*), bx(ldb,*), cx(ldc,*)
  REAL*8, ALLOCATABLE :: bs(:,:), br(:,:)
  REAL*8, ALLOCATABLE :: as(:,:), ar(:,:)
  REAL*8, ALLOCATABLE :: cs(:,:), cr(:,:), cc(:,:)
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic
  INTEGER :: ip, i, ii, j, jj, k, kk, nloc_ip, nloc
  INTEGER :: nloc_src, nloc_dst, idest, isour
  nloc = ldim_cyclic(n, nproc, mpime)
  ALLOCATE( cc( nloc, n ) )
  DO j = 1, n
    DO i = 1, nloc
      cc(i,j) = 0.0d0
    END DO
  END DO
  DO ip = 0, nproc - 1
    isour = MOD(mpime-ip+nproc, nproc)
    idest = MOD(mpime+ip      , nproc)
    nloc_src = ldim_cyclic(n, nproc, isour)
    ALLOCATE( as( nloc, n ) )
    ALLOCATE( ar( nloc_src, n ) )
    DO j = 1, n
      DO i = 1, nloc
        as(i,j) = ax(i,j)
      END DO
    END DO
    !CALL sendrecv_real(as, SIZE(as), idest, ar, SIZE(ar), isour, ip)
#if defined __MPI
    CALL MPI_SENDRECV(as, SIZE(as), MPI_DOUBLE_PRECISION,  &
         IDEST, ip,  ar, SIZE(ar), MPI_DOUBLE_PRECISION, &
         ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
    IF(ierr .NE. 0) THEN
       WRITE(6,*) ' ** ERROR in sendrecv ** '
       STOP
    END IF
#else
    ar = as
#endif
    DO j = 1, n
      kk = isour+1
      DO k = 1, nloc_src
        DO i = 1, nloc
          cc(i,j) = cc(i,j) + bx(i,kk) * ar(k,j)
        END DO
        kk = kk + nproc
      END DO
    END DO
    DEALLOCATE( as )
    DEALLOCATE( ar )
  END DO
  CALL mytrasp_dati(cc, nloc, 'R', cx, nloc, 'R', n, mpime, nproc)
  DEALLOCATE( cc )
  RETURN
END SUBROUTINE



SUBROUTINE mytrasp_dati(ax, lda, dista, bx, ldb, distb, n, mpime, nproc)
  IMPLICIT NONE
#if defined __MPI
  INCLUDE 'mpif.h'
  INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
#endif
  INTEGER :: lda, ldb, n, mpime, nproc
  CHARACTER(LEN=1) :: dista, distb
  REAL*8 :: ax(lda,*), bx(ldb,*)
  REAL*8, ALLOCATABLE :: bs(:,:), br(:,:)
  REAL*8, ALLOCATABLE :: as(:,:), ar(:,:)
  REAL*8, ALLOCATABLE :: cs(:,:), cr(:,:), cc(:,:)
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic
  INTEGER :: ip, i, ii, j, jj, k, kk, nloc_ip, nloc
  INTEGER :: nloc_src, nloc_dst, idest, isour
  CHARACTER(LEN=1) :: da, db

  da = 'R'
  db = 'R'
  IF( dista == 'c' .OR. dista == 'C' ) THEN
    da = 'C'
  END IF
  IF( distb == 'c' .OR. distb == 'C' ) THEN
    db = 'C'
  END IF

  nloc = ldim_cyclic(n, nproc, mpime)

    IF( da == 'R' .AND. db == 'R' ) THEN

      DO ip = 0, nproc - 1
        isour = MOD(mpime-ip+nproc, nproc)
        idest = MOD(mpime+ip      , nproc)
        nloc_src = ldim_cyclic(n, nproc, isour)
        nloc_dst = ldim_cyclic(n, nproc, idest)
        ALLOCATE( cs( nloc, nloc_dst ) )
        ALLOCATE( cr( nloc_src, nloc ) )
        jj = idest+1
        DO j = 1, nloc_dst
          DO i = 1, nloc
            cs(i,j) = ax(i,jj)
          END DO
          jj = jj + nproc
        END DO
        ! CALL sendrecv_real(cs, SIZE(cs), idest, cr, SIZE(cr), isour, ip)
#if defined __MPI
        CALL MPI_SENDRECV(cs, SIZE(cs), MPI_DOUBLE_PRECISION,  &
           IDEST, ip,  cr, SIZE(cr), MPI_DOUBLE_PRECISION, &
           ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
        IF(ierr .NE. 0) THEN
           WRITE(6,*) ' ** ERROR in sendrecv ** '
           STOP
        END IF
#else
        cr = cs
#endif
        jj = isour+1
        DO j = 1, nloc_src
          DO i = 1, nloc
            bx(i,jj) = cr(j,i)
          END DO
          jj = jj + nproc
        END DO
        DEALLOCATE( cs )
        DEALLOCATE( cr )
      END DO

    ELSE IF( da == 'R' .AND. db == 'C' ) THEN

      DO j = 1, n
        DO i = 1, nloc
          bx(j,i) = ax(i,j)
        END DO
      END DO

    ELSE IF( da == 'C' .AND. db == 'R' ) THEN

      DO j = 1, n
        DO i = 1, nloc
          bx(i,j) = ax(j,i)
        END DO
      END DO

    ELSE IF( da == 'C' .AND. db == 'C' ) THEN

      DO ip = 0, nproc - 1
        isour = MOD(mpime-ip+nproc, nproc)
        idest = MOD(mpime+ip      , nproc)
        nloc_src = ldim_cyclic(n, nproc, isour)
        nloc_dst = ldim_cyclic(n, nproc, idest)
        ALLOCATE( cs( nloc_dst, nloc ) )
        ALLOCATE( cr( nloc, nloc_src ) )
        DO i = 1, nloc
          jj = idest+1
          DO j = 1, nloc_dst
            cs(j,i) = ax(jj,i)
            jj = jj + nproc
          END DO
        END DO
        !CALL sendrecv_real(cs, SIZE(cs), idest, cr, SIZE(cr), isour, ip)
#if defined __MPI
        CALL MPI_SENDRECV(cs, SIZE(cs), MPI_DOUBLE_PRECISION,  &
           IDEST, ip,  cr, SIZE(cr), MPI_DOUBLE_PRECISION, &
           ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
        IF(ierr .NE. 0) THEN
           WRITE(6,*) ' ** ERROR in sendrecv ** '
           STOP
        END IF
#else
        cr = cs
#endif
        DO i = 1, nloc
          jj = isour+1
          DO j = 1, nloc_src
            bx(jj,i) = cr(i,j)
            jj = jj + nproc
          END DO
        END DO
        DEALLOCATE( cs )
        DEALLOCATE( cr )
      END DO

    ELSE
      ! ERROR 
    END IF

  RETURN
END SUBROUTINE


SUBROUTINE mytrasp_dist(ax, lda, dista, bx, ldb, distb, n, mpime, nproc)
  IMPLICIT NONE
#if defined __MPI
  INCLUDE 'mpif.h'
  INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
#endif
  INTEGER :: lda, ldb, n, mpime, nproc
  CHARACTER(LEN=1) :: dista, distb
  REAL*8 :: ax(lda,*), bx(ldb,*)
  REAL*8, ALLOCATABLE :: bs(:,:), br(:,:)
  REAL*8, ALLOCATABLE :: as(:,:), ar(:,:)
  REAL*8, ALLOCATABLE :: cs(:,:), cr(:,:), cc(:,:)
  INTEGER :: ldim_cyclic, owner_cyclic, lind_cyclic
  INTEGER :: ip, i, ii, j, jj, k, kk, nloc_ip, nloc
  INTEGER :: nloc_src, nloc_dst, idest, isour
  CHARACTER(LEN=1) :: da, db

  da = 'R'
  db = 'R'
  IF( dista == 'c' .OR. dista == 'C' ) THEN
    da = 'C'
  END IF
  IF( distb == 'c' .OR. distb == 'C' ) THEN
    db = 'C'
  END IF

  nloc = ldim_cyclic(n, nproc, mpime)

    IF( da == 'R' .AND. db == 'R' ) THEN

        DO j = 1, n
          DO i = 1, nloc
            bx(i,j) = ax(i,j)
          END DO
        END DO

    ELSE IF( da == 'R' .AND. db == 'C' ) THEN

      DO ip = 0, nproc - 1
        isour = MOD(mpime-ip+nproc, nproc)
        idest = MOD(mpime+ip      , nproc)
        nloc_src = ldim_cyclic(n, nproc, isour)
        nloc_dst = ldim_cyclic(n, nproc, idest)
        ALLOCATE( cs( nloc, nloc_dst ) )
        ALLOCATE( cr( nloc_src, nloc ) )
        jj = idest+1
        DO j = 1, nloc_dst
          DO i = 1, nloc
            cs(i,j) = ax(i,jj)
          END DO
          jj = jj + nproc
        END DO
        !CALL sendrecv_real(cs, SIZE(cs), idest, cr, SIZE(cr), isour, ip)
#if defined __MPI
        CALL MPI_SENDRECV(cs, SIZE(cs), MPI_DOUBLE_PRECISION,  &
           IDEST, ip,  cr, SIZE(cr), MPI_DOUBLE_PRECISION, &
           ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
        IF(ierr .NE. 0) THEN
           WRITE(6,*) ' ** ERROR in sendrecv ** '
           STOP
        END IF
#else
        cr = cs
#endif
        DO j = 1, nloc
          ii = isour+1
          DO i = 1, nloc_src
            bx(ii,j) = cr(i,j)
            ii = ii + nproc
          END DO
        END DO
        DEALLOCATE( cs )
        DEALLOCATE( cr )
      END DO

    ELSE IF( da == 'C' .AND. db == 'R' ) THEN

      DO ip = 0, nproc - 1
        isour = MOD(mpime-ip+nproc, nproc)
        idest = MOD(mpime+ip      , nproc)
        nloc_src = ldim_cyclic(n, nproc, isour)
        nloc_dst = ldim_cyclic(n, nproc, idest)
        ALLOCATE( cs( nloc_dst, nloc ) )
        ALLOCATE( cr( nloc, nloc_src ) )
        DO i = 1, nloc
          jj = idest+1
          DO j = 1, nloc_dst
            cs(j, i) = ax(jj, i)
            jj = jj + nproc
          END DO
        END DO
        !CALL sendrecv_real(cs, SIZE(cs), idest, cr, SIZE(cr), isour, ip)
#if defined __MPI
        CALL MPI_SENDRECV(cs, SIZE(cs), MPI_DOUBLE_PRECISION,  &
           IDEST, ip,  cr, SIZE(cr), MPI_DOUBLE_PRECISION, &
           ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
        IF(ierr .NE. 0) THEN
           WRITE(6,*) ' ** ERROR in sendrecv ** '
           STOP
        END IF
#else
        cr = cs
#endif
        ii = isour+1
        DO i = 1, nloc_src
          DO j = 1, nloc
            bx(j,ii) = cr(j,i)
          END DO
          ii = ii + nproc
        END DO
        DEALLOCATE( cs )
        DEALLOCATE( cr )
      END DO

    ELSE IF( da == 'C' .AND. db == 'C' ) THEN

        DO j = 1, nloc
          DO i = 1, n
            bx(i,j) = ax(i,j)
          END DO
        END DO

    ELSE
      ! ERROR 
    END IF


  RETURN
END SUBROUTINE

