!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!==----------------------------------------------==!
    MODULE parallel_toolkit
!==----------------------------------------------==!

    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE parallel_include

    IMPLICIT NONE
    SAVE
    PRIVATE

    PUBLIC :: rep_matmul_drv
    PUBLIC :: zrep_matmul_drv

    CONTAINS

!==----------------------------------------------==!
!
!  My parallel blas
!
!==----------------------------------------------==!

SUBROUTINE mattr_drv( m, k, a, lda, b, ldb, nb, dims, coor, comm )
  !
  !  Compute B as the transpose of matrix A 
  !  A and B are distributed on a 2D cartesian processor
  !  grid in a block cyclic way (as in scalapack),
  !  using a block size of NB
  !
  !     B :=  A'
  !
  !  A is a K by M matrix
  !  B is an M by K matrix
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: m, k
  INTEGER, INTENT(IN) :: lda, ldb
  REAL(DP)            :: a(lda,*), b(ldb,*)
  INTEGER, INTENT(IN) :: nb, dims(2), coor(2), comm
  !
#if defined __MPI
  !
  integer ierr
  integer ndims, rowid, colid
  integer coosrc(2), coodst(2), ipsrc, ipdst, mpime
  integer ihsnd, ihrcv
  logical periods(2)
  !
  integer :: iu
  integer :: i, j, nk, nm
  integer :: ii, jj
  integer :: isrc, jsrc
  integer :: idst, jdst
  integer :: itag
  integer :: nmb, nkb 
  integer :: istatus( MPI_STATUS_SIZE )
  real(DP), allocatable :: abuf(:,:)
  !
  integer :: numroc
  integer :: indxg2l
  external :: numroc, indxg2l

  !
  CALL GRID2D_RANK( 'C', dims(1), dims(2), coor(1), coor(2), mpime )
  !
  iu = 200 + mpime
  !
  !  Compute the global number of blocks for matrix dimension
  !
  nmb = ( m + nb - 1 ) / nb  
  nkb = ( k + nb - 1 ) / nb
  !
  ALLOCATE( abuf( nb, nb ) )
  !
  DO i = 1, nmb
     DO j = 1, nkb
        !
        itag = j + nkb * (i-1)
        !
        coosrc(1) = MOD( (j-1), dims(1) )
        coosrc(2) = MOD( (i-1), dims(2) ) 
        !
        coodst(1) = MOD( (i-1), dims(1) )
        coodst(2) = MOD( (j-1), dims(2) ) 
        !
        CALL GRID2D_RANK( 'C', dims(1), dims(2), coosrc(1), coosrc(2), ipsrc )
        CALL GRID2D_RANK( 'C', dims(1), dims(2), coodst(1), coodst(2), ipdst )
        !
        jsrc = INDXG2L( 1 + (j-1)*nb, nb, coor(1), 0, dims(1) )
        isrc = INDXG2L( 1 + (i-1)*nb, nb, coor(2), 0, dims(2) )
        !
        jdst = INDXG2L( 1 + (j-1)*nb, nb, coor(2), 0, dims(2) )
        idst = INDXG2L( 1 + (i-1)*nb, nb, coor(1), 0, dims(1) )
        !
        nk = MIN( nb, k - (j-1)*nb ) !  number of element in the block
        nm = MIN( nb, m - (i-1)*nb ) !  number of element in the block
        !
        IF( ipsrc == ipdst ) THEN
          IF( ipsrc == mpime ) THEN
            DO ii = 1, nm
              DO jj = 1, nk
                b( idst + ii - 1, jdst + jj - 1 ) = a( jsrc + jj - 1, isrc + ii - 1 )
              END DO
            END DO
          END IF
        ELSE
          IF( ipsrc == mpime ) THEN
            DO ii = 1, nm
              DO jj = 1, nk
                abuf( ii, jj ) = a( jsrc + jj - 1, isrc + ii - 1 )
                !
              END DO
            END DO
            CALL MPI_ISEND( abuf, nb*nb, MPI_DOUBLE_PRECISION, ipdst, itag, comm, ihsnd, ierr )
            CALL mpi_wait(ihsnd, istatus, ierr)
          ELSE IF( ipdst == mpime ) THEN
            CALL MPI_IRECV( abuf, nb*nb, MPI_DOUBLE_PRECISION, ipsrc, itag, comm, ihrcv, ierr )
            CALL mpi_wait(ihrcv, istatus, ierr)
            DO jj = 1, nk
              DO ii = 1, nm
                !
                b( idst + ii - 1, jdst + jj - 1 ) = abuf( ii, jj )
              END DO
            END DO
          END IF
        END IF 
        !
     END DO
  END DO

#else

  INTEGER :: i, j

  DO j = 1, k
     DO i = 1, m
        B( i, j ) = A( j, i )
     END DO
  END DO

#endif

  RETURN

END SUBROUTINE mattr_drv

! ---------------------------------------------------------------------------------

SUBROUTINE matsplit_drv( m, k, ar, ldar, a, lda, nb, dims, coor, comm )
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: m, k
  INTEGER, INTENT(IN) :: ldar
  REAL(DP)            :: ar(ldar,*)  !  matrix to be splitted, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  REAL(DP)            :: a(lda,*)
  INTEGER, INTENT(IN) :: nb, coor(2), dims(2), comm
  !
  INTEGER :: i, j, nra, nca, ii, jj
  !
  INTEGER  :: numroc, INDXL2G
  EXTERNAL :: numroc, INDXL2G

  nra = NUMROC( m, nb, coor(1), 0, dims(1) )  !  total number of local row for matrix A, C
  nca = NUMROC( k, nb, coor(2), 0, dims(2) )  !  total number of local columns of A

  do j = 1, nca
     jj = INDXL2G( j, NB, coor(2), 0, dims(2) )
     do i = 1, nra
        ii = INDXL2G( i, NB, coor(1), 0, dims(1) )
        a( i, j ) = ar( ii, jj )
     end do
  end do

  RETURN

END SUBROUTINE matsplit_drv

! ---------------------------------------------------------------------------------

SUBROUTINE matmerge_drv( m, k, a, lda, ar, ldar, nb, dims, coor, comm )
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: m, k
  INTEGER, INTENT(IN) :: ldar
  REAL(DP)            :: ar(ldar,*)  !  matrix to be merged, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  REAL(DP)            :: a(lda,*)
  INTEGER, INTENT(IN) :: nb, coor(2), dims(2), comm
  !
  INTEGER :: i, j, ii, jj, ierr

#if defined __MPI
  !

  INTEGER :: jsrc, isrc, ipsrc, coosrc(2)
  INTEGER :: nmb, nkb, nk, nm, mpime

  REAL(DP), ALLOCATABLE :: buf(:,:)
  !
  INTEGER  :: INDXG2L
  EXTERNAL :: INDXG2L
  !
  CALL GRID2D_RANK( 'C', dims(1), dims(2), coor(1), coor(2), mpime )

  nmb = ( m + nb - 1 ) / nb  
  nkb = ( k + nb - 1 ) / nb

  ALLOCATE( buf( nb, nb ) )

  DO j = 1, nkb
     DO i = 1, nmb
        !
        coosrc(1) = MOD( (i-1), dims(1) )
        coosrc(2) = MOD( (j-1), dims(2) )
        !
        CALL GRID2D_RANK( 'C', dims(1), dims(2), coosrc(1), coosrc(2), ipsrc )
        !
        isrc = INDXG2L( 1 + (i-1)*nb, nb, coor(1), 0, dims(1) )
        jsrc = INDXG2L( 1 + (j-1)*nb, nb, coor(2), 0, dims(2) )
        !
        nm = MIN( nb, m - (i-1)*nb ) !  number of element in the block
        nk = MIN( nb, k - (j-1)*nb ) !  number of element in the block

        IF( ipsrc == mpime ) THEN
           DO jj = 1, nk
              DO ii = 1, nm
                 buf( ii, jj ) = a( isrc + ii - 1, jsrc + jj - 1 )
              END DO
           END DO
        ENDIF
        !
        CALL MPI_BCAST( buf, nb*nb, MPI_DOUBLE_PRECISION, ipsrc, comm, ierr )
        !
        do jj = 1, nk
           do ii = 1, nm
              ar( ii + (i-1)*nb, jj + (j-1)*nb ) = buf( ii, jj )
           end do
        end do
        !
     END DO
  END DO
  !
  DEALLOCATE( buf )

#else

  DO j = 1, k
     DO i = 1, m
        ar( i, j ) = a( i, j )
     END DO
  END DO

#endif

  RETURN
END SUBROUTINE matmerge_drv

! ---------------------------------------------------------------------------------

SUBROUTINE matscal_drv( m, n, beta, c, ldc, nb, dims, coor, comm )
  !
  implicit none
  !
  INTEGER,  INTENT(IN) :: m, n
  REAL(DP), INTENT(IN) :: beta
  INTEGER,  INTENT(IN) :: ldc
  REAL(DP)             :: c(ldc,*)
  INTEGER,  INTENT(IN) :: nb, coor(2), dims(2), comm
  !
  INTEGER :: i, j, nr, nc, ierr
  !
  INTEGER  :: numroc
  EXTERNAL :: numroc

  nr  = NUMROC( m, nb, coor(1), 0, dims(1) )  ! local row of C
  nc  = NUMROC( n, nb, coor(2), 0, dims(2) )  ! local colum of C 

  IF( beta == 0.0_DP ) THEN
    do j = 1, nc
      do i = 1, nr
        c(i,j) = 0.0_DP
      end do
    end do 
  ELSE
    do j = 1, nc
      do i = 1, nr
        c(i,j) = beta * c(i,j)
      end do
    end do 
  END IF

  RETURN

END SUBROUTINE

! ---------------------------------------------------------------------------------

SUBROUTINE matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, nb, dims, coor, comm )
  !
  implicit none
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa, transb
  INTEGER,  INTENT(IN) :: m, n, k
  REAL(DP), INTENT(IN) :: alpha, beta
  INTEGER,  INTENT(IN) :: lda, ldb, ldc
  REAL(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER, INTENT(IN) :: nb, dims(2), coor(2), comm
  !
  !  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !
  !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !
  !  WHERE  OP( X ) IS ONE OF
  !
  !     OP( X ) = X   OR   OP( X ) = X',
  !
  !  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !
  !
  !

#if defined __MPI
  !
  integer ierr
  integer ndims, rowid, colid
  integer comm_row, comm_col
  !
  integer :: ib, jb, kb, ibl, kbl, jbl
  integer :: i, j, kk, ni, nj, nk, nm, il, jl
  integer :: nnb, nmb, nkb 
  integer :: nr, nra, nca, nc, nrb, ncb, ii, jj
  integer :: nrt, ncat, nct, nrbt
  real(DP), allocatable :: abuf(:,:), bbuf(:,:)
  real(DP), allocatable :: at(:,:)
  real(DP), allocatable :: bt(:,:)
  !
  integer :: numroc
  integer :: indxg2l
  external :: numroc, indxg2l
  !
  IF( dims(1) * dims(2) == 1 ) THEN

     !  if there is only one proc no need of using parallel alg.

     call dgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )

     RETURN

  END IF
  !

  CALL MPI_COMM_SPLIT( COMM, coor(2), coor(1), COMM_COL, IERR )
  CALL MPI_COMM_RANK( COMM_COL, rowid, IERR )
  !
  CALL MPI_COMM_SPLIT( COMM, coor(1), coor(2), COMM_ROW, IERR )
  CALL MPI_COMM_RANK( COMM_ROW, colid, IERR )
  !
  !  Compute the global number of blocks for matrix dimension
  !
  nmb = ( m + nb - 1 ) / nb  
  !
  nnb = ( n + nb - 1 ) / nb
  !
  nkb = ( k + nb - 1 ) / nb
  !
  !  Compute the total number of local row for matrix A, C
  !
  nr  = NUMROC( m, nb, coor(1), 0, dims(1) )  ! local row of C
  !
  nra = NUMROC( m, nb, coor(1), 0, dims(1) )  ! local row of OP( A )
  nca = NUMROC( k, nb, coor(2), 0, dims(2) )  ! local columns of OP( A )
  !
  nrb = NUMROC( k, nb, coor(1), 0, dims(1) )  ! local row of OP( B )
  ncb = NUMROC( n, nb, coor(2), 0, dims(2) )  ! local colum of OP( B )
  !
  nc  = NUMROC( n, nb, coor(2), 0, dims(2) )  ! local colum of C 
  !
  IF( transa == 'T' .OR. transa == 't' ) THEN
    !
    ALLOCATE( at( nra, nca ) )
    !
    CALL mattr_drv( m, k, a, lda, at, nra, nb, dims, coor, comm )
    !
  END IF
  !
  IF( transb == 'T' .OR. transb == 't' ) THEN
    !
    ALLOCATE( bt( nrb, ncb ) )
    !
    CALL mattr_drv( k, n, b, ldb, bt, nrb, nb, dims, coor, comm )
    !
  END IF
  !
  !  Scale matrix C
  !
  CALL matscal_drv( m, n, beta, c, ldc, nb, dims, coor, comm )
  !
  !  loop over the rows/columns blocks of matrix OP(A)/OP(B)
  !
  do kb = 1, nkb
    !
    kk  = ( kb - 1 ) * nb + 1  !  first element of the block (global index)
    nk = MIN( nb, k - kk + 1 ) !  number of element in the block

    colid = MOD( (kb-1), dims(2) )  ! processor owning the block
    rowid = MOD( (kb-1), dims(1) )

    allocate( abuf( nr, nk ) )

    if( colid == coor(2) ) then
      nrt = 0
      ibl = 0
      kbl = INDXG2L( 1 + (kb-1)*nb, nb, coor(2), 0, dims(2) )
      do ib = 1 + coor(1), nmb, dims(1)
        i = ( ib - 1 ) * nb + 1
        ni = MIN( nb, m - i + 1 )
        IF( transa == 'T' .OR. transa == 't' ) THEN
          do jj = 1, nk
            do ii = 1, ni
              abuf( ii + nrt, jj ) = at( ii + ibl*nb, jj + kbl - 1 )
            end do
          end do
        ELSE
          do jj = 1, nk
            do ii = 1, ni
              abuf( ii + nrt, jj ) = a( ii + ibl*nb, jj + kbl - 1 )
            end do
          end do
        END IF
        nrt = nrt + ni
        ibl = ibl + 1
      end do
    end if
    CALL MPI_BCAST( abuf(1,1), nr*nk, MPI_DOUBLE_PRECISION, colid, COMM_ROW, IERR )

    allocate( bbuf( nk, nc ) )

    if( rowid == coor(1) ) then
      nct = 0 
      jbl = 0
      kbl = INDXG2L( 1 + (kb-1)*nb, nb, coor(1), 0, dims(1) )
      do jb = 1 + coor(2), nnb, dims(2)
        j = ( jb - 1 ) * nb + 1
        nj = MIN( nb, n - j + 1 )
        IF( transb == 'T' .OR. transb == 't' ) THEN
          do jj = 1, nj
            do ii = 1, nk
              bbuf( ii, jj + nct ) = bt( ii + kbl - 1, jj + jbl*nb )
            end do
          end do
        ELSE
          do jj = 1, nj
            do ii = 1, nk
              bbuf( ii, jj + nct ) = b( ii + kbl - 1, jj + jbl*nb )
            end do
          end do
        END IF
        nct = nct + nj
        jbl = jbl + 1
      end do
    end if

    CALL MPI_BCAST( bbuf(1,1), nk*nc, MPI_DOUBLE_PRECISION, rowid, COMM_COL, IERR )

    ii = 1
    do ib = 1 + coor(1), nmb, dims(1)
      i = ( ib - 1 ) * nb + 1
      il = INDXG2L( i, nb, coor(1), 0, dims(1) )
      ni = MIN( nb, m - i + 1 )
      jj = 1
      do jb = 1 + coor(2), nnb, dims(2)
        j = ( jb - 1 ) * nb + 1
        jl = INDXG2L( j, nb, coor(2), 0, dims(2) )
        nj = MIN( nb, n - j + 1 )
        call dgemm( 'n', 'n', ni, nj, nk, alpha, abuf( ii, 1 ), nra, bbuf( 1, jj ), nk, 1.0_DP, c( il, jl ), ldc )
        jj = jj + nj
      end do
      ii = ii + ni
    end do

    deallocate( abuf )
    deallocate( bbuf )

  end do

  IF( ALLOCATED( at ) ) DEALLOCATE( at )
  IF( ALLOCATED( bt ) ) DEALLOCATE( bt )

  CALL MPI_COMM_FREE(COMM_ROW, ierr)
  CALL MPI_COMM_FREE(COMM_COL, ierr)


#else

     !  if we are not compiling with __MPI this is equivalent to a blas call

     call dgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )

#endif

  RETURN

END SUBROUTINE

!==----------------------------------------------==!
!
! Copyright (C) 2005 Carlo Cavazzoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE rep_matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, comm )
  !
  !  Parallel matrix multiplication with replicated matrix
  !
  implicit none
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa, transb
  INTEGER, INTENT(IN) :: m, n, k
  REAL(DP), INTENT(IN) :: alpha, beta
  INTEGER, INTENT(IN) :: lda, ldb, ldc
  REAL(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER, INTENT(IN) :: comm
  !
  !  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !
  !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !
  !  WHERE  OP( X ) IS ONE OF
  !
  !     OP( X ) = X   OR   OP( X ) = X',
  !
  !  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !
  !
  !

#if defined __MPI

  !

  INTEGER :: ME, I, II, J, JJ, IP, SOUR, DEST, INFO, IERR, ioff, ldx
  INTEGER :: NB, IB_S, NB_SOUR, IB_SOUR, IBUF
  INTEGER :: nproc, mpime, q, r

  REAL(8), ALLOCATABLE :: auxa( : )
  REAL(8), ALLOCATABLE :: auxc( : )

  !
  ! ... BODY
  !

  CALL MPI_COMM_SIZE(comm, NPROC, IERR)
  CALL MPI_COMM_RANK(comm, MPIME, IERR)

  IF ( NPROC == 1 ) THEN

     !  if there is only one proc no need of using parallel alg.

     CALL DGEMM(TRANSA, TRANSB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

     RETURN

  END IF

  ME = MPIME + 1
  Q = INT( m / NPROC )
  R = MOD( m , NPROC )

  ! ... Find out the number of elements in the local block
  !     along "M" first dimension os matrix A

  NB = Q
  IF( ME <= R ) NB = NB + 1

  ! ... Find out the global index of the local first row

  IF( ME <= R ) THEN
     ib_s = (Q+1)*(ME-1) + 1
  ELSE
     ib_s = Q*(ME-1) + R + 1
  END IF

  ldx = m / nproc + 1

  ALLOCATE( auxa( MAX( n, m ) * ldx ) )
  ALLOCATE( auxc( MAX( n, m ) * ldx ) )

  IF( TRANSA == 'N' .OR. TRANSA == 'n' ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = A( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  ELSE
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = A( J, I + ioff )
        END DO
        ibuf = ibuf + ldx
     END DO
     !ioff = ib_s - 1
     !call mytranspose( A( 1, ioff + 1 ), lda, auxa(1), ldx, m, nb)
  END IF

  IF( beta /= 0.0_DP ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, n
        DO I = 1, NB
           auxc( ibuf + I ) = C( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  END IF

  CALL DGEMM( 'N', transb, nb, n, k, alpha, auxa(1), ldx, B, ldb, beta, auxc(1), ldx )

  ! ... Here processors exchange blocks

  DO IP = 0, NPROC-1

     ! ...    Find out the number of elements in the block of processor SOUR

     NB_SOUR = q
     IF( (IP+1) .LE. r ) NB_SOUR = NB_SOUR+1

     ! ...    Find out the global index of the first row owned by SOUR

     IF( (IP+1) .LE. r ) THEN
        ib_sour = (Q+1)*IP + 1
     ELSE
        ib_sour = Q*IP + R + 1
     END IF

     IF( mpime == ip ) auxa(1:n*ldx) = auxc(1:n*ldx)

     CALL MPI_BCAST( auxa(1), ldx*n, mpi_double_precision, ip, comm, IERR)

     IBUF = 0
     ioff = IB_SOUR - 1
     DO J = 1, N
        DO I = 1, NB_SOUR
           C( I + ioff, J ) = AUXA( IBUF + I )
        END DO
        IBUF = IBUF + ldx
     END DO

  END DO

  DEALLOCATE( auxa, auxc )

#else

     !  if we are not compiling with __MPI this is equivalent to a blas call

     CALL DGEMM(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)

#endif

  RETURN

END SUBROUTINE rep_matmul_drv


SUBROUTINE zrep_matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, comm )
  !
  !  Parallel matrix multiplication with replicated matrix
  !
  implicit none
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa, transb
  INTEGER, INTENT(IN) :: m, n, k
  COMPLEX(DP), INTENT(IN) :: alpha, beta
  INTEGER, INTENT(IN) :: lda, ldb, ldc
  COMPLEX(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER, INTENT(IN) :: comm
  !
  !  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !
  !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !
  !  WHERE  OP( X ) IS ONE OF
  !
  !     OP( X ) = X   OR   OP( X ) = X',
  !
  !  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !
  !
  !

#if defined __MPI

  !

  INTEGER :: ME, I, II, J, JJ, IP, SOUR, DEST, INFO, IERR, ioff, ldx
  INTEGER :: NB, IB_S, NB_SOUR, IB_SOUR, IBUF
  INTEGER :: nproc, mpime, q, r

  COMPLEX(DP), ALLOCATABLE :: auxa( : )
  COMPLEX(DP), ALLOCATABLE :: auxc( : )

  !
  ! ... BODY
  !

  CALL MPI_COMM_SIZE(comm, NPROC, IERR)
  CALL MPI_COMM_RANK(comm, MPIME, IERR)

  IF ( NPROC == 1 ) THEN

     !  if there is only one proc no need of using parallel alg.

     CALL ZGEMM(TRANSA, TRANSB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

     RETURN

  END IF

  ME = MPIME + 1
  Q = INT( m / NPROC )
  R = MOD( m , NPROC )

  ! ... Find out the number of elements in the local block
  !     along "M" first dimension os matrix A

  NB = Q
  IF( ME <= R ) NB = NB + 1

  ! ... Find out the global index of the local first row

  IF( ME <= R ) THEN
     ib_s = (Q+1)*(ME-1) + 1
  ELSE
     ib_s = Q*(ME-1) + R + 1
  END IF

  ldx = m / nproc + 1

  ALLOCATE( auxa( MAX( n, m ) * ldx ) )
  ALLOCATE( auxc( MAX( n, m ) * ldx ) )

  IF( TRANSA == 'N' .OR. TRANSA == 'n' ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = A( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  ELSE
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = CONJG( A( J, I + ioff ) )
        END DO
        ibuf = ibuf + ldx
     END DO
     !ioff = ib_s - 1
     !call mytranspose( A( 1, ioff + 1 ), lda, auxa(1), ldx, m, nb)
  END IF

  IF( beta /= 0.0_DP ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, n
        DO I = 1, NB
           auxc( ibuf + I ) = C( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  END IF

  CALL ZGEMM( 'N', transb, nb, n, k, alpha, auxa(1), ldx, B, ldb, beta, auxc(1), ldx )

  ! ... Here processors exchange blocks

  DO IP = 0, NPROC-1

     ! ...    Find out the number of elements in the block of processor SOUR

     NB_SOUR = q
     IF( (IP+1) .LE. r ) NB_SOUR = NB_SOUR+1

     ! ...    Find out the global index of the first row owned by SOUR

     IF( (IP+1) .LE. r ) THEN
        ib_sour = (Q+1)*IP + 1
     ELSE
        ib_sour = Q*IP + R + 1
     END IF

     IF( mpime == ip ) auxa(1:n*ldx) = auxc(1:n*ldx)

     CALL MPI_BCAST( auxa(1), ldx*n, mpi_double_complex, ip, comm, IERR)

     IBUF = 0
     ioff = IB_SOUR - 1
     DO J = 1, N
        DO I = 1, NB_SOUR
           C( I + ioff, J ) = AUXA( IBUF + I )
        END DO
        IBUF = IBUF + ldx
     END DO

  END DO

  DEALLOCATE( auxa, auxc )

#else

     !  if we are not compiling with __MPI this is equivalent to a blas call

     CALL ZGEMM(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)

#endif

  RETURN

END SUBROUTINE zrep_matmul_drv


!==----------------------------------------------==!
END MODULE parallel_toolkit
!==----------------------------------------------==!

!
! ... some simple routines for parallel linear algebra (the matrices are
! ... always replicated on all the cpus)
!
! ... written by carlo sbraccia ( 2006 )
!
!----------------------------------------------------------------------------
SUBROUTINE para_dgemm( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (splitting matrix B by columns) of DGEMM 
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_bcast
  USE parallel_include
  USE parallel_toolkit
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1), INTENT(IN)    :: transa, transb
  INTEGER,          INTENT(IN)    :: m, n, k
  REAL(DP),         INTENT(IN)    :: alpha, beta
  INTEGER,          INTENT(IN)    :: lda, ldb, ldc
  REAL(DP),         INTENT(INOUT) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER,          INTENT(IN)    :: comm
  !
  INTEGER              :: i, mpime, nproc, ierr
  INTEGER              :: ncol, i0, i1
  INTEGER, ALLOCATABLE :: i0a(:), i1a(:)
  !
  ! ... quick return if possible
  !
  IF ( m == 0 .OR. n == 0 .OR. &
       ( ( alpha == 0.0_DP .OR. k == 0 ) .AND. beta == 1.0_DP ) ) RETURN
  !
#if defined (__XD1)
  !
  CALL rep_matmul_drv( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  RETURN
  !
#endif

#if defined (__MPI)
  !
  CALL MPI_COMM_SIZE( comm, nproc, ierr )
  CALL MPI_COMM_RANK( comm, mpime, ierr )
  !
#else
  !
  nproc = 1
  mpime = 0
  !
#endif
  !
  ncol = n / nproc
  !
  ALLOCATE( i0a( 0:nproc-1 ), i1a( 0:nproc-1 ) )
  !
  DO i = 0, nproc -1
     !
     i0a(i) = i*ncol + 1
     i1a(i) = ( i + 1 )*ncol
     !
  END DO
  !
  i1a(nproc-1) = n
  !
  i0 = i0a(mpime)
  i1 = i1a(mpime)
  !
  IF ( transb == 'n' .OR. transb == 'N' ) THEN
     !
     CALL DGEMM( transa, transb, m, i1 - i0 + 1, k, &
                 alpha, a, lda, b(1,i0), ldb, beta, c(1,i0), ldc )
     !
  ELSE
     !
     CALL DGEMM( transa, transb, m, i1 - i0 + 1, k, &
                 alpha, a, lda, b(i0,1), ldb, beta, c(1,i0), ldc )
     !
  END IF
  !
  DO i = 0 , nproc - 1
     !
     CALL mp_bcast( c(1:ldc,i0a(i):i1a(i)), i, comm )
     !
  END DO
  !
  DEALLOCATE( i0a, i1a )
  !
  RETURN
  !
END SUBROUTINE para_dgemm
!
!----------------------------------------------------------------------------
SUBROUTINE para_zgemm( transa, transb, m, n, k, &
                       alpha, a, lda, b, ldb, beta, c, ldc, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (splitting matrix B by columns) of ZGEMM
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_bcast
  USE parallel_include
  USE parallel_toolkit
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1), INTENT(IN)    :: transa, transb
  INTEGER,          INTENT(IN)    :: m, n, k
  COMPLEX(DP),      INTENT(IN)    :: alpha, beta
  INTEGER,          INTENT(IN)    :: lda, ldb, ldc
  COMPLEX(DP),      INTENT(INOUT) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER,          INTENT(IN)    :: comm
  !
  INTEGER                :: i, mpime, nproc, ierr
  INTEGER                :: ncol, i0, i1
  INTEGER, ALLOCATABLE   :: i0a(:), i1a(:)
  !
  COMPLEX(DP), PARAMETER :: ONE = (1.0_DP,0.0_DP), ZERO = ( 0.0_DP, 0.0_DP )
  !
  ! ... quick return if possible
  !
  IF ( m == 0 .OR. n == 0 .OR. &
       ( ( alpha == 0.0_DP .OR. k == 0 ) .AND. beta == ONE ) ) RETURN
  !
#if defined (__XD1)
  !
  CALL zrep_matmul_drv( transa, transb, m, n, k, &
                        alpha, a, lda, b, ldb, beta, c, ldc, comm )
  RETURN
  !
#endif
  !
#if defined (__MPI)
  !
  CALL MPI_COMM_SIZE( comm, nproc, ierr )
  CALL MPI_COMM_RANK( comm, mpime, ierr )
  !
#else
  !
  nproc = 1
  mpime = 0
  !
#endif
  !
  ncol = n / nproc
  !
  ALLOCATE( i0a( 0:nproc-1 ), i1a( 0:nproc-1 ) )
  !
  DO i = 0, nproc -1
     !
     i0a(i) = i*ncol + 1
     i1a(i) = ( i + 1 )*ncol
     !
  END DO
  !
  i1a(nproc-1) = n
  !
  i0 = i0a(mpime)
  i1 = i1a(mpime)
  !
  IF ( transb == 'n' .OR. transb == 'N' ) THEN
     !
     CALL ZGEMM( transa, transb, m, i1 - i0 + 1, k, &
                 alpha, a, lda, b(1,i0), ldb, beta, c(1,i0), ldc )
     !
  ELSE
     !
     CALL ZGEMM( transa, transb, m, i1 - i0 + 1, k, &
                 alpha, a, lda, b(i0,1), ldb, beta, c(1,i0), ldc )
     !
  END IF
  !
  DO i = 0 , nproc - 1
     !
     CALL mp_bcast( c(1:ldc,i0a(i):i1a(i)), i, comm )
     !
  END DO
  !
  DEALLOCATE( i0a, i1a )
  !
  RETURN
  !
END SUBROUTINE para_zgemm
!
!----------------------------------------------------------------------------
SUBROUTINE para_dgemv( trans, m, n, alpha, &
                       a, lda, x, incx, beta, y, incy, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (splitting matrix A by rows) of DGEMV
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_bcast
  USE parallel_include
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1), INTENT(IN)    :: trans
  INTEGER,          INTENT(IN)    :: m, n
  REAL(DP),         INTENT(IN)    :: alpha, beta
  INTEGER,          INTENT(IN)    :: lda, incx, incy
  REAL(DP),         INTENT(INOUT) :: a(lda,*)
  REAL(DP),         INTENT(INOUT) :: x(*), y(*)
  INTEGER,          INTENT(IN)    :: comm
  !
  INTEGER               :: i, j, mpime, nproc, ierr
  INTEGER               :: nrow, i0, i1, dim, ydum_size
  INTEGER,  ALLOCATABLE :: i0a(:), i1a(:)
  REAL(DP), ALLOCATABLE :: ydum(:)
  !
  ! ... quick return if possible
  !
  IF ( m == 0 .OR. n == 0 .OR. &
       ( alpha == 0.0_DP .AND. beta == 1.0_DP ) ) RETURN
  !
#if defined (__MPI)
  !
  CALL MPI_COMM_SIZE( comm, nproc, ierr )
  CALL MPI_COMM_RANK( comm, mpime, ierr )
  !
#else
  !
  nproc = 1
  mpime = 0
  !
#endif
  !
  nrow = m / nproc
  !
  ALLOCATE( i0a( 0:nproc-1 ), i1a( 0:nproc-1 ) )
  !
  DO i = 0, nproc -1
     !
     i0a(i) = i*nrow + 1
     i1a(i) = ( i + 1 )*nrow
     !
  END DO
  !
  i1a(nproc-1) = m
  !
  i0 = i0a(mpime)
  i1 = i1a(mpime)
  !
  IF ( trans == 'n' .OR. trans == 'N' ) THEN
     !
     dim = ( 1 + ( m - 1 )*ABS( incy ) )
     !
     ydum_size = m
     !
  ELSE
     !
     dim = ( 1 + ( n - 1 )*ABS( incy ) )
     !
     ydum_size = n
     !
  END IF
  !
  ALLOCATE( ydum( ydum_size ) )
  !
  i = 0
  !
  DO j = 1, dim, incy
     !
     i = i + 1
     !
     IF ( i < i0 .OR. i > i1 ) CYCLE
     !
     ydum(i) = y(j)
     !
  END DO
  !
  IF ( trans == 'n' .OR. trans == 'N' ) THEN
     !
     CALL DGEMV( trans, i1 - i0 + 1, n, &
                 alpha, a(i0,1), lda, x, incx, beta, ydum(i0), 1 )
     !
  ELSE
     !
     CALL DGEMV( trans, i1 - i0 + 1, n, &
                 alpha, a(1,i0), lda, x, incx, beta, ydum(i0), 1 )
     !
  END IF
  !
  DO i = 0 , nproc - 1
     !
     CALL mp_bcast( ydum(i0a(i):i1a(i)), i, comm )
     !
  END DO
  !
  i = 0
  !
  DO j = 1, dim, incy
     !
     i = i + 1
     !
     y(j) = ydum(i)
     !
  END DO
  !
  DEALLOCATE( ydum, i0a, i1a )
  !
  RETURN
  !
END SUBROUTINE para_dgemv
!
!----------------------------------------------------------------------------
SUBROUTINE para_zgemv( trans, m, n, alpha, &
                       a, lda, x, incx, beta, y, incy, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (splitting matrix A by rows) of ZGEMV
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_bcast
  USE parallel_include
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1), INTENT(IN)    :: trans
  INTEGER,          INTENT(IN)    :: m, n
  COMPLEX(DP),      INTENT(IN)    :: alpha, beta
  INTEGER,          INTENT(IN)    :: lda, incx, incy
  COMPLEX(DP),      INTENT(INOUT) :: a(lda,*)
  COMPLEX(DP),      INTENT(INOUT) :: x(*), y(*)
  INTEGER,          INTENT(IN)    :: comm
  !
  INTEGER                  :: i, j, mpime, nproc, ierr
  INTEGER                  :: nrow, i0, i1, dim, ydum_size
  INTEGER,     ALLOCATABLE :: i0a(:), i1a(:)
  COMPLEX(DP), ALLOCATABLE :: ydum(:)
  !
  ! ... quick return if possible
  !
  IF ( m == 0 .OR. n == 0 .OR. &
       ( alpha == 0.0_DP .AND. beta == 1.0_DP ) ) RETURN
  !
#if defined (__MPI)
  !
  CALL MPI_COMM_SIZE( comm, nproc, ierr )
  CALL MPI_COMM_RANK( comm, mpime, ierr )
  !
#else
  !
  nproc = 1
  mpime = 0
  !
#endif
  !
  nrow = m / nproc
  !
  ALLOCATE( i0a( 0:nproc-1 ), i1a( 0:nproc-1 ) )
  !
  DO i = 0, nproc -1
     !
     i0a(i) = i*nrow + 1
     i1a(i) = ( i + 1 )*nrow
     !
  END DO
  !
  i1a(nproc-1) = m
  !
  i0 = i0a(mpime)
  i1 = i1a(mpime)
  !
  IF ( trans == 'n' .OR. trans == 'N' ) THEN
     !
     dim = ( 1 + ( m - 1 )*ABS( incy ) )
     !
     ydum_size = m
     !
  ELSE
     !
     dim = ( 1 + ( n - 1 )*ABS( incy ) )
     !
     ydum_size = n
     !
  END IF
  !
  ALLOCATE( ydum( ydum_size ) )
  !
  i = 0
  !
  DO j = 1, dim, incy
     !
     i = i + 1
     !
     IF ( i < i0 .OR. i > i1 ) CYCLE
     !
     ydum(i) = y(j)
     !
  END DO
  !
  IF ( trans == 'n' .OR. trans == 'N' ) THEN
     !
     CALL ZGEMV( trans, i1 - i0 + 1, n, &
                 alpha, a(i0,1), lda, x, incx, beta, ydum(i0), 1 )
     !
  ELSE
     !
     CALL ZGEMV( trans, i1 - i0 + 1, n, &
                 alpha, a(1,i0), lda, x, incx, beta, ydum(i0), 1 )
     !
  END IF
  !
  DO i = 0 , nproc - 1
     !
     CALL mp_bcast( ydum(i0a(i):i1a(i)), i, comm )
     !
  END DO
  !
  i = 0
  !
  DO j = 1, dim, incy
     !
     i = i + 1
     !
     y(j) = ydum(i)
     !
  END DO  
  !
  DEALLOCATE( ydum, i0a, i1a )
  !
  RETURN
  !
END SUBROUTINE para_zgemv
!
!----------------------------------------------------------------------------
SUBROUTINE para_dcholdc( n, a, lda, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (using a parallel version of DGEMV) of
  ! ... the Cholesky decomposition (equivalent to DPOTF2)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: n
  INTEGER,  INTENT(IN)    :: lda
  REAL(DP), INTENT(INOUT) :: a(lda,*)
  INTEGER,  INTENT(IN)    :: comm
  !
  INTEGER            :: i, j
  REAL(DP)           :: aii
  REAL(DP), EXTERNAL :: DDOT
  !
  !
  DO i = 1, n
     !
     aii = a(i,i) - DDOT( i-1, a(i,1), lda, a(i,1), lda )
     !
     IF ( aii < 0.0_DP ) &
        CALL errore( 'para_dcholdc', 'a is not positive definite', i )
     !
     aii = SQRT( aii )
     !
     a(i,i) = aii
     !
     IF ( i < n ) THEN
        !
        CALL para_dgemv( 'N', n-i, i-1, -1.0_DP, a(i+1,1), &
                         lda, a(i,1), lda, 1.0_DP, a(i+1,i), 1, comm )
        !
        CALL DSCAL( n-i, 1.0_DP / aii, a(i+1,i), 1 )
        !
     END IF
     !
  END DO
  !
  FORALL( i = 1:n, j = 1:n, j > i ) a(i,j) = 0.0_DP
  !
  RETURN
  !
END SUBROUTINE para_dcholdc
!
!----------------------------------------------------------------------------
SUBROUTINE para_zcholdc( n, a, lda, comm )
  !----------------------------------------------------------------------------
  !
  ! ... trivial parallelization (using a parallel version of ZGEMV) of
  ! ... the Cholesky decomposition (equivalent to ZPOTF2)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: n
  INTEGER,     INTENT(IN)    :: lda
  COMPLEX(DP), INTENT(INOUT) :: a(lda,*)
  INTEGER,     INTENT(IN)    :: comm
  !
  INTEGER               :: i, j
  REAL(DP)              :: aii
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  COMPLEX(DP), PARAMETER :: ONE = (1.0_DP,0.0_DP), ZERO = ( 0.0_DP, 0.0_DP )
  !
  !
  DO i = 1, n
     !
     aii = REAL( a(i,i) ) - ZDOTC( i-1, a(i,1), lda, a(i,1), lda )
     !
     IF ( aii < 0.0_DP ) &
        CALL errore( 'para_zcholdc', 'a is not positive definite', i )
     !
     aii = SQRT( aii )
     !
     a(i,i) = aii
     !
     IF ( i < n ) THEN
        !
        CALL ZLACGV( i-1, a(i,1), lda )
        !
        CALL para_zgemv( 'N', n-i, i-1, -ONE, a(i+1,1), &
                         lda, a(i,1), lda, ONE, a(i+1,i), 1, comm )
        !
        CALL ZLACGV( i-1, a(i,1), lda )
        !
        CALL ZDSCAL( n-i, 1.0_DP / aii, a(i+1,i), 1 )
        !
     END IF
     !
  END DO
  !
  FORALL( i = 1:n, j = 1:n, j > i ) a(i,j) = ZERO
  !
  RETURN
  !
END SUBROUTINE para_zcholdc
!
!----------------------------------------------------------------------------
SUBROUTINE para_dtrtri( n, a, lda, comm )
  !----------------------------------------------------------------------------
  !
  ! ... parallel inversion of a lower triangular matrix done distributing
  ! ... by columns ( the number of columns assigned to each processor are
  ! ... chosen to optimize the load balance )
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_bcast
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: n
  INTEGER,  INTENT(IN)    :: lda
  REAL(DP), INTENT(INOUT) :: a(lda,*)
  INTEGER,  INTENT(IN)    :: comm
  !
  INTEGER               :: i, j, k
  INTEGER               :: i0, i1, mpime, nproc, ierr
  INTEGER,  ALLOCATABLE :: i0a(:), i1a(:)
  REAL(DP)              :: an, xfrac
  REAL(DP)              :: sum
  REAL(DP), ALLOCATABLE :: inva(:,:)
  !
  !
#if defined (__MPI)
  !
  CALL MPI_COMM_SIZE( comm, nproc, ierr )
  CALL MPI_COMM_RANK( comm, mpime, ierr )
  !
#else
  !
  nproc = 1
  mpime = 0
  !
#endif
  !
  ALLOCATE( i0a( 0:nproc-1 ), i1a( 0:nproc-1 ) )
  !
  an = 1.0_DP / DBLE( nproc )
  !
  i0a(0) = 1
  !
  DO i = 0, nproc - 2
     !
     xfrac = 1.0_DP - SQRT( 1.0_DP - DBLE( i+1 )*an )
     !
     i1a(i)   = ANINT( xfrac*n )
     i0a(i+1) = i1a(i) + 1
     !
  END DO
  !
  i1a(nproc-1) = n
  !
  i0 = i0a(mpime)
  i1 = i1a(mpime)
  !
  ALLOCATE( inva( n, i0:i1 ) )
  !
  inva(:,:) = 0.0_DP
  !
  DO i = i0, i1
     !
     inva(i,i) = 1.0_DP / a(i,i)
     !
     DO j = i + 1, n
        !
        sum = 0.0_DP
        !
        DO k = i, j - 1
           !
           sum = sum + inva(k,i)*a(j,k)
           !
        END DO
        !
        inva(j,i) = - sum / a(j,j)
        !
     END DO
     !
  END DO
  !
  a(1:lda,1:n) = 0.0_DP
  a(1:n,i0:i1) = inva(:,:)
  !
  DEALLOCATE( inva )
  !
  DO i = 0 , nproc - 1
     !
#if defined __XD1
     CALL BCAST_REAL( a(1,i0a(i)), i1a(i)-i0a(i)+1, i, comm, ierr )
#else
     CALL mp_bcast( a(1:n,i0a(i):i1a(i)), i, comm )
#endif
     !
  END DO
  !
  DEALLOCATE( i0a, i1a )
  !
  RETURN
  !
END SUBROUTINE para_dtrtri
!
!----------------------------------------------------------------------------
SUBROUTINE para_ztrtri( n, a, lda, comm )
  !----------------------------------------------------------------------------
  !
  ! ... parallel inversion of a lower triangular matrix done distributing
  ! ... by columns ( the number of columns assigned to each processor are
  ! ... chosen to optimize the load balance in the limit of large matrices )
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_bcast
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: n
  INTEGER,     INTENT(IN)    :: lda
  COMPLEX(DP), INTENT(INOUT) :: a(lda,*)
  INTEGER,     INTENT(IN)    :: comm
  !
  INTEGER                  :: i, j, k
  INTEGER                  :: i0, i1, mpime, nproc, ierr
  INTEGER,     ALLOCATABLE :: i0a(:), i1a(:)
  REAL(DP)                 :: an, xfrac
  COMPLEX(DP)              :: sum
  COMPLEX(DP), ALLOCATABLE :: inva(:,:)
  !
  COMPLEX(DP), PARAMETER :: ONE = (1.0_DP,0.0_DP), ZERO = ( 0.0_DP, 0.0_DP )
  !
  !
#if defined (__MPI)
  !
  CALL MPI_COMM_SIZE( comm, nproc, ierr )
  CALL MPI_COMM_RANK( comm, mpime, ierr )
  !
#else
  !
  nproc = 1
  mpime = 0
  !
#endif
  !
  ALLOCATE( i0a( 0:nproc-1 ), i1a( 0:nproc-1 ) )
  !
  an = 1.0_DP / DBLE( nproc )
  !
  i0a(0) = 1
  !
  DO i = 0, nproc - 2
     !
     xfrac = 1.0_DP - SQRT( 1.0_DP - DBLE( i+1 )*an )
     !
     i1a(i)   = ANINT( xfrac*n )
     i0a(i+1) = i1a(i) + 1
     !
  END DO
  !
  i1a(nproc-1) = n
  !
  i0 = i0a(mpime)
  i1 = i1a(mpime)
  !
  ALLOCATE( inva( n, i0:i1 ) )
  !
  inva(:,:) = ZERO
  !
  DO i = i0, i1
     !
     inva(i,i) = ONE / a(i,i)
     !
     DO j = i + 1, n
        !
        sum = ZERO
        !
        DO k = i, j - 1
           !
           sum = sum + inva(k,i)*a(j,k)
           !
        END DO
        !
        inva(j,i) = - sum / a(j,j)
        !
     END DO
     !
  END DO
  !
  a(1:lda,1:n) = ZERO
  a(1:n,i0:i1) = inva(:,:)
  !
  DEALLOCATE( inva )
  !
  DO i = 0 , nproc - 1
     !
     CALL mp_bcast( a(1:n,i0a(i):i1a(i)), i, comm )
     !
  END DO
  !
  DEALLOCATE( i0a, i1a )
  !
  RETURN
  !
END SUBROUTINE para_ztrtri

!=----------------------------------------------------------------------------=!
!
!
!  Cannon's algorithms for parallel matrix multiplication
!  written by Carlo Cavazzoni
!  
!
!

SUBROUTINE sqr_mm_cannon( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, desc )
   !
   !  Parallel square matrix multiplication with Cannon's algorithm
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , la_npc_ , &
                           la_comm_ , lambda_node_ , la_npr_ , la_npc_ , la_myr_ , la_myc_
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   REAL(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
   INTEGER, INTENT(IN) :: desc(*)
   !
   !  performs one of the matrix-matrix operations
   !
   !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
   !
   !  where  op( x ) is one of
   !
   !     OP( X ) = X   OR   OP( X ) = X',
   !
   !  alpha and beta are scalars, and a, b and c are square matrices
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr
   integer :: np
   integer :: i, j, nr, nc, nb, iter, rowid, colid
   logical :: ta, tb
   INTEGER :: comm
   !
   !
   real(DP), allocatable :: bblk(:,:), ablk(:,:)
   !
#if defined (__MPI)
   !
   integer :: istatus( MPI_STATUS_SIZE )
   !
#endif
   !
   IF( desc( lambda_node_ ) < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF

   IF( desc( la_npr_ ) == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL dgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) &
      CALL errore( ' sqr_mm_cannon ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' sqr_mm_cannon ', ' n less or equal zero ', 1 )
   IF( n < desc( la_npr_ ) ) &
      CALL errore( ' sqr_mm_cannon ', ' n less than the mesh size ', 1 )
   IF( n < desc( nlax_ ) ) &
      CALL errore( ' sqr_mm_cannon ', ' n less than the block size ', 1 )

   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc( la_npr_ )
   comm  = desc( la_comm_ )
   rowid = desc( la_myr_  )
   colid = desc( la_myc_  )
   !
   !  Retrieve the size of the local block
   !
   nr    = desc( nlar_ ) 
   nc    = desc( nlac_ ) 
   nb    = desc( nlax_ )
   !
#if defined (__MPI)
   CALL MPI_BARRIER( comm, ierr )
#endif
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   !
   allocate( bblk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         bblk( i, j ) = b( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         bblk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         bblk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   !
   ta = ( TRANSA == 'T' .OR. TRANSA == 't' )
   tb = ( TRANSB == 'T' .OR. TRANSB == 't' )
   !
   !  Shift A rowid+1 places to the west
   ! 
   IF( ta ) THEN
      CALL shift_exch_block( ablk, 'W', 1 )
   ELSE
      CALL shift_block( ablk, 'W', rowid+1, 1 )
   END IF
   !
   !  Shift B colid+1 places to the north
   ! 
   IF( tb ) THEN
      CALL shift_exch_block( bblk, 'N', np+1 )
   ELSE
      CALL shift_block( bblk, 'N', colid+1, np+1 )
   END IF
   !
   !  Accumulate on C
   !
   CALL dgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, beta, c, ldc)
   !
   DO iter = 2, np
      !
      !  Shift A 1 places to the east
      ! 
      CALL shift_block( ablk, 'E', 1, iter )
      !
      !  Shift B 1 places to the south
      ! 
      CALL shift_block( bblk, 'S', 1, np+iter )
      !
      !  Accumulate on C
      !
      CALL dgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, 1.0_DP, c, ldc)
      !
   END DO

   deallocate( ablk, bblk )
   
   RETURN

CONTAINS

   SUBROUTINE shift_block( blk, dir, ln, tag )
      !
      !   Block shift 
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir      ! shift direction
      INTEGER,          INTENT(IN) :: ln       ! shift lenght
      INTEGER,          INTENT(IN) :: tag      ! communication tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      IF( dir == 'W' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid - ln + np, np )
         icsrc = MOD( colid + ln + np, np )
         !
      ELSE IF( dir == 'E' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid + ln + np, np )
         icsrc = MOD( colid - ln + np, np )
         !
      ELSE IF( dir == 'N' ) THEN

         irdst = MOD( rowid - ln + np, np )
         irsrc = MOD( rowid + ln + np, np )
         icdst = colid
         icsrc = colid

      ELSE IF( dir == 'S' ) THEN

         irdst = MOD( rowid + ln + np, np )
         irsrc = MOD( rowid - ln + np, np )
         icdst = colid
         icsrc = colid

      ELSE

         CALL errore( ' sqr_mm_cannon ', ' unknown shift direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_PRECISION, &
           idest, tag, isour, tag, comm, istatus, ierr)
      !
#endif
      RETURN
   END SUBROUTINE shift_block

   SUBROUTINE shift_exch_block( blk, dir, tag )
      !
      !   Combined block shift and exchange
      !   only used for the first step
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir
      INTEGER,          INTENT(IN) :: tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      INTEGER :: icol, irow
      !
      IF( dir == 'W' ) THEN
         !
         icol = rowid
         irow = colid
         !
         irdst = irow
         icdst = MOD( icol - irow-1 + np, np )
         !
         irow = rowid
         icol = MOD( colid + rowid+1 + np, np )
         !
         irsrc = icol
         icsrc = irow
         !
      ELSE IF( dir == 'N' ) THEN
         !
         icol = rowid
         irow = colid
         !
         icdst = icol
         irdst = MOD( irow - icol-1 + np, np )
         !
         irow = MOD( rowid + colid+1 + np, np )
         icol = colid
         !
         irsrc = icol
         icsrc = irow

      ELSE

         CALL errore( ' sqr_mm_cannon ', ' unknown shift_exch direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_PRECISION, &
           idest, tag, isour, tag, comm, istatus, ierr)
      !
#endif
      RETURN
   END SUBROUTINE shift_exch_block

END SUBROUTINE sqr_mm_cannon


!=----------------------------------------------------------------------------=!

SUBROUTINE sqr_zmm_cannon( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, desc )
   !
   !  Parallel square matrix multiplication with Cannon's algorithm
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , la_npc_ , &
                           la_comm_ , lambda_node_ , la_npr_ , la_npc_ , la_myr_ , la_myc_
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   COMPLEX(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   COMPLEX(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
   INTEGER, INTENT(IN) :: desc(*)
   !
   !  performs one of the matrix-matrix operations
   !
   !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
   !
   !  where  op( x ) is one of
   !
   !     OP( X ) = X   OR   OP( X ) = X',
   !
   !  alpha and beta are scalars, and a, b and c are square matrices
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   INTEGER :: ierr
   INTEGER :: np
   INTEGER :: i, j, nr, nc, nb, iter, rowid, colid
   LOGICAL :: ta, tb
   INTEGER :: comm
   !
   !
   COMPLEX(DP), ALLOCATABLE :: bblk(:,:), ablk(:,:)
   COMPLEX(DP) :: zone = ( 1.0_DP, 0.0_DP )
   COMPLEX(DP) :: zzero = ( 0.0_DP, 0.0_DP )
   !
#if defined (__MPI)
   !
   integer :: istatus( MPI_STATUS_SIZE )
   !
#endif
   !
   IF( desc( lambda_node_ ) < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF

   IF( desc( la_npr_ ) == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL zgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) &
      CALL errore( ' sqr_zmm_cannon ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' sqr_zmm_cannon ', ' n less or equal zero ', 1 )
   IF( n < desc( la_npr_ ) ) &
      CALL errore( ' sqr_zmm_cannon ', ' n less than the mesh size ', 1 )
   IF( n < desc( nlax_ ) ) &
      CALL errore( ' sqr_zmm_cannon ', ' n less than the block size ', 1 )

   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc( la_npr_ )
   comm  = desc( la_comm_ )
   rowid = desc( la_myr_  )
   colid = desc( la_myc_  )
   !
   !  Retrieve the size of the local block
   !
   nr    = desc( nlar_ ) 
   nc    = desc( nlac_ ) 
   nb    = desc( nlax_ )
   !
#if defined (__MPI)
   CALL MPI_BARRIER( comm, ierr )
#endif
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = zzero
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = zzero
      END DO
   END DO
   !
   !
   allocate( bblk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         bblk( i, j ) = b( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         bblk( i, j ) = zzero
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         bblk( i, j ) = zzero
      END DO
   END DO
   !
   !
   ta = ( TRANSA == 'C' .OR. TRANSA == 'c' )
   tb = ( TRANSB == 'C' .OR. TRANSB == 'c' )
   !
   !  Shift A rowid+1 places to the west
   ! 
   IF( ta ) THEN
      CALL shift_exch_block( ablk, 'W', 1 )
   ELSE
      CALL shift_block( ablk, 'W', rowid+1, 1 )
   END IF
   !
   !  Shift B colid+1 places to the north
   ! 
   IF( tb ) THEN
      CALL shift_exch_block( bblk, 'N', np+1 )
   ELSE
      CALL shift_block( bblk, 'N', colid+1, np+1 )
   END IF
   !
   !  Accumulate on C
   !
   CALL zgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, beta, c, ldc)
   !
   DO iter = 2, np
      !
      !  Shift A 1 places to the east
      ! 
      CALL shift_block( ablk, 'E', 1, iter )
      !
      !  Shift B 1 places to the south
      ! 
      CALL shift_block( bblk, 'S', 1, np+iter )
      !
      !  Accumulate on C
      !
      CALL zgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, zone, c, ldc)
      !
   END DO

   deallocate( ablk, bblk )
   
   RETURN

CONTAINS

   SUBROUTINE shift_block( blk, dir, ln, tag )
      !
      !   Block shift 
      !
      IMPLICIT NONE
      COMPLEX(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir      ! shift direction
      INTEGER,          INTENT(IN) :: ln       ! shift lenght
      INTEGER,          INTENT(IN) :: tag      ! communication tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      IF( dir == 'W' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid - ln + np, np )
         icsrc = MOD( colid + ln + np, np )
         !
      ELSE IF( dir == 'E' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid + ln + np, np )
         icsrc = MOD( colid - ln + np, np )
         !
      ELSE IF( dir == 'N' ) THEN

         irdst = MOD( rowid - ln + np, np )
         irsrc = MOD( rowid + ln + np, np )
         icdst = colid
         icsrc = colid

      ELSE IF( dir == 'S' ) THEN

         irdst = MOD( rowid + ln + np, np )
         irsrc = MOD( rowid - ln + np, np )
         icdst = colid
         icsrc = colid

      ELSE

         CALL errore( ' sqr_zmm_cannon ', ' unknown shift direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_COMPLEX, &
           idest, tag, isour, tag, comm, istatus, ierr)
      !
#endif
      RETURN
   END SUBROUTINE shift_block
   !
   SUBROUTINE shift_exch_block( blk, dir, tag )
      !
      !   Combined block shift and exchange
      !   only used for the first step
      !
      IMPLICIT NONE
      COMPLEX(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir
      INTEGER,          INTENT(IN) :: tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      INTEGER :: icol, irow
      !
      IF( dir == 'W' ) THEN
         !
         icol = rowid
         irow = colid
         !
         irdst = irow
         icdst = MOD( icol - irow-1 + np, np )
         !
         irow = rowid
         icol = MOD( colid + rowid+1 + np, np )
         !
         irsrc = icol
         icsrc = irow
         !
      ELSE IF( dir == 'N' ) THEN
         !
         icol = rowid
         irow = colid
         !
         icdst = icol
         irdst = MOD( irow - icol-1 + np, np )
         !
         irow = MOD( rowid + colid+1 + np, np )
         icol = colid
         !
         irsrc = icol
         icsrc = irow

      ELSE

         CALL errore( ' sqr_zmm_cannon ', ' unknown shift_exch direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_COMPLEX, &
           idest, tag, isour, tag, comm, istatus, ierr)
      !
#endif
      RETURN
   END SUBROUTINE shift_exch_block

END SUBROUTINE sqr_zmm_cannon

!
!
!
!

SUBROUTINE sqr_tr_cannon( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix transposition with Cannon's algorithm
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , la_npc_ , &
                           la_comm_ , lambda_node_ , la_npr_ , la_myr_ , la_myc_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP)            :: a(lda,*), b(ldb,*)
   INTEGER, INTENT(IN) :: desc(*)
   !
#if defined (__MPI)
   !
   INCLUDE 'mpif.h'
   !
#endif
   !
   INTEGER :: ierr
   INTEGER :: np, rowid, colid
   INTEGER :: i, j, nr, nc, nb
   INTEGER :: comm
   !
   REAL(DP), ALLOCATABLE :: ablk(:,:)
   !
#if defined (__MPI)
   !
   INTEGER :: istatus( MPI_STATUS_SIZE )
   !
#endif
   !
   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   IF( desc( la_npr_ ) == 1 ) THEN
      CALL mytranspose( a, lda, b, ldb, n, n )
      RETURN
   END IF

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) &
      CALL errore( ' sqr_tr_cannon ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' sqr_tr_cannon ', ' n less or equal zero ', 1 )
   IF( n < desc( la_npr_ ) ) &
      CALL errore( ' sqr_tr_cannon ', ' n less than the mesh size ', 1 )
   IF( n < desc( nlax_ ) ) &
      CALL errore( ' sqr_tr_cannon ', ' n less than the block size ', 1 )

   comm = desc( la_comm_ )


   rowid = desc( la_myr_ )
   colid = desc( la_myc_ )
   np    = desc( la_npr_ )
   !
   !  Compute the size of the local block
   !
   nr = desc( nlar_ ) 
   nc = desc( nlac_ ) 
   nb = desc( nlax_ )
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   CALL exchange_block( ablk )
   !
#if defined (__MPI)
   CALL MPI_BARRIER( comm, ierr )
#endif
   !
   DO j = 1, nr
      DO i = 1, nc
         b( j, i ) = ablk( i, j )
      END DO
   END DO
   !
   deallocate( ablk )
   
   RETURN

CONTAINS

   SUBROUTINE exchange_block( blk )
      !
      !   Block exchange ( transpose )
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      irdst = colid
      icdst = rowid
      irsrc = colid
      icsrc = rowid
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_PRECISION, &
           idest, np+np+1, isour, np+np+1, comm, istatus, ierr)
      !
#endif

      RETURN
   END SUBROUTINE


END SUBROUTINE

!
!
!

SUBROUTINE cyc2blk_redist( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (input) is cyclically distributed by rows across processors
   !  B (output) is distributed by block across 2D processors grid
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , lambda_node_ , la_npr_ , &
                           descla_siz_ , la_npc_ , la_n_ , la_me_ , la_comm_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP) :: a(lda,*), b(ldb,*)
   INTEGER  :: desc(*)
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr, itag
   integer :: np, ip, me, nproc, comm_a
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   !
   real(DP), allocatable :: rcvbuf(:,:,:)
   real(DP), allocatable :: sndbuf(:,:)
   integer, allocatable :: ip_desc(:,:)
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   np = desc( la_npr_ )  !  dimension of the processor mesh
   nb = desc( nlax_ )    !  leading dimension of the local matrix block
   me = desc( la_me_ )   !  my processor id (starting from 0)
   comm_a = desc( la_comm_ )
   nproc  = desc( la_npr_ ) * desc( la_npc_ )

   IF( np /= desc( la_npc_ ) ) &
      CALL errore( ' cyc2blk_redist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' cyc2blk_redist ', ' n less or equal zero ', 1 )
   IF( desc( la_n_ ) < nproc ) &
      CALL errore( ' cyc2blk_redist ', ' nb less than the number of proc ', 1 )

   ALLOCATE( ip_desc( descla_siz_ , nproc ) )

   CALL mpi_allgather( desc, descla_siz_ , mpi_integer, ip_desc, descla_siz_ , mpi_integer, comm_a, ierr )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   
   DO ip = 0, nproc - 1
      !
      IF( ip_desc( lambda_node_ , ip + 1 ) > 0 ) THEN

         ip_nr = ip_desc( nlar_ , ip + 1) 
         ip_nc = ip_desc( nlac_ , ip + 1) 
         ip_ir = ip_desc( ilar_ , ip + 1) 
         ip_ic = ip_desc( ilac_ , ip + 1) 
         !
         DO j = 1, ip_nc
            jj = j + ip_ic - 1
            il = 1
            DO i = 1, ip_nr
               ii = i + ip_ir - 1
               IF( MOD( ii - 1, nproc ) == me ) THEN
                  sndbuf( il, j ) = a( ( ii - 1 )/nproc + 1, jj )
                  il = il + 1
               END IF 
            END DO
         END DO
   
      END IF

      CALL mpi_gather( sndbuf, nbuf, mpi_double_precision, &
                       rcvbuf, nbuf, mpi_double_precision, ip, comm_a, ierr )
     
   END DO

   !
   nr  = desc( nlar_ ) 
   nc  = desc( nlac_ )
   ir  = desc( ilar_ )
   ic  = desc( ilac_ )
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               b( i, j ) = rcvbuf( il, j, ip+1 )
               il = il + 1
            END IF 
         END DO
      END DO
   END DO
   !
   !
   DEALLOCATE( ip_desc )
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   b( 1:n, 1:n ) = a( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE cyc2blk_redist


SUBROUTINE cyc2blk_zredist( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (input) is cyclically distributed by rows across processors
   !  B (output) is distributed by block across 2D processors grid
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , lambda_node_ , la_npr_ , &
                           descla_siz_ , la_npc_ , la_n_ , la_me_ , la_comm_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   COMPLEX(DP) :: a(lda,*), b(ldb,*)
   INTEGER  :: desc(*)
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr, itag
   integer :: np, ip, me, nproc, comm_a
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   !
   COMPLEX(DP), allocatable :: rcvbuf(:,:,:)
   COMPLEX(DP), allocatable :: sndbuf(:,:)
   integer, allocatable :: ip_desc(:,:)
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   np = desc( la_npr_ )  !  dimension of the processor mesh
   nb = desc( nlax_ )    !  leading dimension of the local matrix block
   me = desc( la_me_ )   !  my processor id (starting from 0)
   comm_a = desc( la_comm_ )
   nproc  = desc( la_npr_ ) * desc( la_npc_ )

   IF( np /= desc( la_npc_ ) ) &
      CALL errore( ' cyc2blk_zredist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' cyc2blk_zredist ', ' n less or equal zero ', 1 )
   IF( desc( la_n_ ) < nproc ) &
      CALL errore( ' cyc2blk_zredist ', ' nb less than the number of proc ', 1 )

   ALLOCATE( ip_desc( descla_siz_ , nproc ) )

   CALL mpi_allgather( desc, descla_siz_ , mpi_integer, ip_desc, descla_siz_ , mpi_integer, comm_a, ierr )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   
   DO ip = 0, nproc - 1
      !
      IF( ip_desc( lambda_node_ , ip + 1 ) > 0 ) THEN

         ip_nr = ip_desc( nlar_ , ip + 1) 
         ip_nc = ip_desc( nlac_ , ip + 1) 
         ip_ir = ip_desc( ilar_ , ip + 1) 
         ip_ic = ip_desc( ilac_ , ip + 1) 
         !
         DO j = 1, ip_nc
            jj = j + ip_ic - 1
            il = 1
            DO i = 1, ip_nr
               ii = i + ip_ir - 1
               IF( MOD( ii - 1, nproc ) == me ) THEN
                  sndbuf( il, j ) = a( ( ii - 1 )/nproc + 1, jj )
                  il = il + 1
               END IF 
            END DO
         END DO
   
      END IF

      CALL mpi_gather( sndbuf, nbuf, mpi_double_complex, &
                       rcvbuf, nbuf, mpi_double_complex, ip, comm_a, ierr )
     
   END DO

   !
   nr  = desc( nlar_ ) 
   nc  = desc( nlac_ )
   ir  = desc( ilar_ )
   ic  = desc( ilac_ )
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               b( i, j ) = rcvbuf( il, j, ip+1 )
               il = il + 1
            END IF 
         END DO
      END DO
   END DO
   !
   !
   DEALLOCATE( ip_desc )
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   b( 1:n, 1:n ) = a( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE cyc2blk_zredist


SUBROUTINE blk2cyc_redist( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (output) is cyclically distributed by rows across processors
   !  B (input) is distributed by block across 2D processors grid
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , lambda_node_ , la_npr_ , &
                           descla_siz_ , la_npc_ , la_n_ , la_me_ , la_comm_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP) :: a(lda,*), b(ldb,*)
   INTEGER  :: desc(*)
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr, itag
   integer :: np, ip, me, comm_a, nproc
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   !
   real(DP), allocatable :: rcvbuf(:,:,:)
   real(DP), allocatable :: sndbuf(:,:)
   integer, allocatable :: ip_desc(:,:)
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   np = desc( la_npr_ )  !  dimension of the processor mesh
   nb = desc( nlax_ )    !  leading dimension of the local matrix block
   me = desc( la_me_ )   !  my processor id (starting from 0)
   comm_a = desc( la_comm_ )
   nproc  = desc( la_npr_ ) * desc( la_npc_ )

   IF( np /= desc( la_npc_ ) ) &
      CALL errore( ' blk2cyc_redist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' blk2cyc_redist ', ' n less or equal zero ', 1 )
   IF( desc( la_n_ ) < nproc ) &
      CALL errore( ' blk2cyc_redist ', ' nb less than the number of proc ', 1 )

   ALLOCATE( ip_desc( descla_siz_ , nproc ) )

   CALL mpi_allgather( desc, descla_siz_ , mpi_integer, ip_desc, descla_siz_ , mpi_integer, comm_a, ierr )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   !
   nr  = desc( nlar_ ) 
   nc  = desc( nlac_ )
   ir  = desc( ilar_ )
   ic  = desc( ilac_ )
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               sndbuf( il, j ) = b( i, j )
               il = il + 1
            END IF 
         END DO
      END DO
      CALL mpi_gather( sndbuf, nbuf, mpi_double_precision, &
                       rcvbuf, nbuf, mpi_double_precision, ip, comm_a, ierr )
   END DO
   !
   
   DO ip = 0, nproc - 1
      !
      IF( ip_desc( lambda_node_ , ip + 1 ) > 0 ) THEN

         ip_nr = ip_desc( nlar_ , ip + 1) 
         ip_nc = ip_desc( nlac_ , ip + 1) 
         ip_ir = ip_desc( ilar_ , ip + 1) 
         ip_ic = ip_desc( ilac_ , ip + 1) 
         !
         DO j = 1, ip_nc
            jj = j + ip_ic - 1
            il = 1
            DO i = 1, ip_nr
               ii = i + ip_ir - 1
               IF( MOD( ii - 1, nproc ) == me ) THEN
                  a( ( ii - 1 )/nproc + 1, jj ) = rcvbuf( il, j, ip+1 )
                  il = il + 1
               END IF 
            END DO
         END DO
   
      END IF
     
   END DO
   !
   DEALLOCATE( ip_desc )
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   a( 1:n, 1:n ) = b( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE blk2cyc_redist


SUBROUTINE blk2cyc_zredist( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (output) is cyclically distributed by rows across processors
   !  B (input) is distributed by block across 2D processors grid
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , lambda_node_ , la_npr_ , &
                           descla_siz_ , la_npc_ , la_n_ , la_me_ , la_comm_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   COMPLEX(DP) :: a(lda,*), b(ldb,*)
   INTEGER  :: desc(*)
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr, itag
   integer :: np, ip, me, comm_a, nproc
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   !
   COMPLEX(DP), allocatable :: rcvbuf(:,:,:)
   COMPLEX(DP), allocatable :: sndbuf(:,:)
   integer, allocatable :: ip_desc(:,:)
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   np = desc( la_npr_ )  !  dimension of the processor mesh
   nb = desc( nlax_ )    !  leading dimension of the local matrix block
   me = desc( la_me_ )   !  my processor id (starting from 0)
   comm_a = desc( la_comm_ )
   nproc  = desc( la_npr_ ) * desc( la_npc_ )

   IF( np /= desc( la_npc_ ) ) &
      CALL errore( ' blk2cyc_zredist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' blk2cyc_zredist ', ' n less or equal zero ', 1 )
   IF( desc( la_n_ ) < nproc ) &
      CALL errore( ' blk2cyc_zredist ', ' nb less than the number of proc ', 1 )

   ALLOCATE( ip_desc( descla_siz_ , nproc ) )

   CALL mpi_allgather( desc, descla_siz_ , mpi_integer, ip_desc, descla_siz_ , mpi_integer, comm_a, ierr )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   !
   nr  = desc( nlar_ ) 
   nc  = desc( nlac_ )
   ir  = desc( ilar_ )
   ic  = desc( ilac_ )
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               sndbuf( il, j ) = b( i, j )
               il = il + 1
            END IF 
         END DO
      END DO
      CALL mpi_gather( sndbuf, nbuf, mpi_double_complex, &
                       rcvbuf, nbuf, mpi_double_complex, ip, comm_a, ierr )
   END DO
   !
   
   DO ip = 0, nproc - 1
      !
      IF( ip_desc( lambda_node_ , ip + 1 ) > 0 ) THEN

         ip_nr = ip_desc( nlar_ , ip + 1) 
         ip_nc = ip_desc( nlac_ , ip + 1) 
         ip_ir = ip_desc( ilar_ , ip + 1) 
         ip_ic = ip_desc( ilac_ , ip + 1) 
         !
         DO j = 1, ip_nc
            jj = j + ip_ic - 1
            il = 1
            DO i = 1, ip_nr
               ii = i + ip_ir - 1
               IF( MOD( ii - 1, nproc ) == me ) THEN
                  a( ( ii - 1 )/nproc + 1, jj ) = rcvbuf( il, j, ip+1 )
                  il = il + 1
               END IF 
            END DO
         END DO
   
      END IF
     
   END DO
   !
   DEALLOCATE( ip_desc )
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   a( 1:n, 1:n ) = b( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE blk2cyc_zredist
!
!
!
SUBROUTINE pzpotf( sll, ldx, n, desc )
   use descriptors, ONLY: descla_local_dims, descla_siz_ , la_myr_ , la_myc_ , la_me_ ,&
                          nlar_ , nlac_ , ilar_ , ilac_ , la_comm_ , la_nx_ , la_npr_ , la_npc_
   use parallel_include
   implicit none
   integer :: n, ldx
   integer :: desc( descla_siz_ )
   real*8  :: one, zero
   complex*16 :: sll( ldx, ldx ), cone, czero
   integer :: myrow, mycol, ierr
   integer :: jb, info, ib, kb
   integer :: jnr, jir, jic, jnc
   integer :: inr, iir, iic, inc
   integer :: knr, kir, kic, knc
   integer :: nr, nc
   integer :: rcomm, ccomm, color, key, myid, np
   complex*16, allocatable :: ssnd( :, : ), srcv( :, : )

   one   = 1.0d0
   cone  = 1.0d0
   zero  = 0.0d0
   czero = 0.0d0

#if defined __MPI

   myrow = desc( la_myr_ )
   mycol = desc( la_myc_ )
   myid  = desc( la_me_ )
   np    = desc( la_npr_ )

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) THEN
      CALL errore( ' pzpotf ', ' only square grid are allowed ', 1 ) 
   END IF

   nr = desc( nlar_ )
   nc = desc( nlac_ )

   ALLOCATE( ssnd( ldx, ldx ) )
   ALLOCATE( srcv( ldx, ldx ) )

   DO jb = 1, np
      !
      !    Update and factorize the current diagonal block and test
      !    for non-positive-definiteness.
      !
      CALL descla_local_dims( jir, jnr, n, desc( la_nx_ ), np, jb-1 )
      !
      !    since we loop on diagonal blocks/procs we have jnc == jnr
      !
      jnc = jnr   
      !
      !    prepare row and colum communicators
      IF( ( myrow >= ( jb-1 ) ) .AND. ( mycol <= ( jb-1 ) ) ) THEN 
          color = mycol
          key   = myrow
      ELSE
          color = np
          key   = myid
      END IF
      !
      CALL mpi_comm_split( desc( la_comm_ ) , color, key, ccomm, ierr )
      !  
      IF( myrow >= jb-1 .and. mycol <= jb-1 ) THEN
          color = myrow
          key   = mycol
      ELSE
          color = np
          key   = myid
      END IF
      !
      CALL mpi_comm_split( desc( la_comm_ ), color, key, rcomm, ierr )
      !
      !    here every process can work independently, then we need a reduce.
      !
      IF( jb > 1 ) THEN
         !
         DO ib = 1, jb - 1
            IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( ib - 1 ) ) ) THEN
               !
               !  remember: matrix ssnd is nr*nr, and procs on the diagonale have nr == nc 
               !
               CALL ZHERK( 'L', 'N', nr, nc, -ONE, sll, ldx, zero, ssnd, ldx )
               !
            END IF
         END DO
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            ssnd = sll
         END IF
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol <= ( jb - 1 ) ) ) THEN
            !
            !  accumulate on the diagonal block/proc
            !
            CALL MPI_REDUCE( ssnd, sll, ldx*ldx, MPI_DOUBLE_COMPLEX, MPI_SUM, jb-1, rcomm, ierr )
            !
         END IF
         !
      END IF
      !
      ! Only proj ( jb-1, jb-1 ) operates this
      !
      info = 0
      !
      IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
         CALL ZPOTF2( 'L', jnr, sll, ldx, INFO )
         IF( info /= 0 ) &
            CALL errore( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )
      END IF
      !
      IF( ( jb > 1 ) .AND. ( jb < np ) ) THEN
         !
         !           Compute the current block column.
         !
         ! processors ( 1 : jb - 1, jb ) should bcast their blocs
         ! along column to processor ( 1 : jb - 1, jb + 1 : nb )
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol < ( jb - 1 ) ) ) THEN
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol < ( jb - 1 ) ) ) THEN
            CALL mpi_bcast( srcv, ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
         END IF
         !
         DO ib = jb + 1, np
            CALL descla_local_dims( iir, inr, n, desc( la_nx_ ), np, ib-1 )
            DO kb = 1, jb - 1
               CALL descla_local_dims( kic, knc, n, desc( la_nx_ ), np, kb-1 )
               IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( kb - 1 ) ) ) THEN
                  CALL ZGEMM( 'N', 'C', inr, jnr, knc, -CONE, sll, ldx, srcv, ldx, czero, ssnd, ldx )
               END IF
            END DO
            IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
               ssnd = sll
            END IF
         END DO
         !
         ! processors ( jb, jb + 1 : nb ) should collect block along row,
         ! from processors  ( 1 : jb - 1, jb + 1 : nb )
         !
         DO kb = jb + 1, np
            IF( ( myrow == ( kb - 1 ) ) .AND. ( mycol <= ( jb - 1 ) ) ) THEN
               IF( ( jb == 1 ) ) THEN
                  IF( mycol == ( jb - 1 ) ) THEN
                     sll = ssnd
                  END IF
               ELSE
                  CALL MPI_REDUCE( ssnd, sll, ldx*ldx, MPI_DOUBLE_COMPLEX, MPI_SUM, jb-1, rcomm, ierr )
               END IF
            END IF
         END DO
         !
      END IF
      !
      IF( jb < np ) THEN
         !
         ! processor "jb,jb" should broadcast his block to procs ( jb+1 : nb, jb )
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            CALL mpi_bcast( srcv,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
         END IF
         !
         DO ib = jb + 1, np
            IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
               CALL ZTRSM( 'R', 'L', 'C', 'N', nr, nc, CONE, srcv, ldx, sll, ldx )
            END IF
         END DO
         !
      END IF
      !
      CALL mpi_comm_free( rcomm, ierr )
      CALL mpi_comm_free( ccomm, ierr )
      !
   END DO

   DEALLOCATE( srcv, ssnd )

#else

   CALL ZPOTRF( 'L', n, sll, ldx, info )

   IF( info /= 0 ) &
      CALL errore( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )

#endif

   return
END SUBROUTINE


  SUBROUTINE pztrtri ( sll, ldx, n, desc )
    
    ! pztrtri computes the parallel inversion of a lower triangular matrix 
    ! distribuited among the processes using a 2-D cyclic block partitioning. 
    ! The algorithm is based on the schema below and executes the model 
    ! recursively to each column C2 under the diagonal.     
    !
    !     |-------|-------|      |--------------------|--------------------|
    !     |   A1  |   0   |      |   C1 = trtri(A1)   |          0         |
    ! A = |-------|-------|  C = |--------------------|--------------------|
    !     |   A2  |   A3  |      | C2 = -C3 * A2 * C1 |   C3 = trtri(A3)   | 
    !     |-------|-------|      |--------------------|--------------------|
    !
    ! The recursive steps of multiplication (C2 = -C3 * A2 * C1) is based on the Cannon's algorithms 
    ! for parallel matrix multiplication and is done with BLACS(dgemm)
    !
    !
    ! Arguments
    ! ============
    !
    ! sll   = local block of data
    ! ldx    = maximum dimension of one block
    ! n     = size of the global array diributed among the blocks
    ! mpime = The coordinates of the process whose local array row or 
    !         column is to be determined.
    ! np    = size of the grid of processors
    ! comm  = global comunicator of the grid processors
    !
    !
    !

    USE kinds
    USE parallel_include
    USE descriptors, ONLY: descla_local_dims, descla_siz_ , la_myr_ , la_myc_ , la_me_ ,&
                          nlar_ , nlac_ , ilar_ , ilac_ , la_comm_ , la_nx_ , la_npr_ , la_npc_

    IMPLICIT NONE

    INTEGER, INTENT( IN ) :: n, ldx
    INTEGER, INTENT( IN ) :: desc( descla_siz_ )
    COMPLEX(DP), INTENT( INOUT ) :: sll( ldx, ldx )

    COMPLEX(DP), PARAMETER :: ONE = (1.0d0, 0.0d0)
    COMPLEX(DP), PARAMETER :: ZERO = (0.0d0, 0.0d0)

#if defined __MPI
    INTEGER :: status(MPI_STATUS_SIZE)
#endif
    INTEGER :: req(2), ierr, col_comm
    INTEGER :: send, recv, group_rank, group_size
    INTEGER :: myrow, mycol, np, myid, comm

    ! counters
    INTEGER :: k, i, j, count, step_count, shiftcount, cicle 
    INTEGER :: C3dim   ! Dimension of submatrix B
    INTEGER :: nc, nr ! Local dimension of block
    INTEGER :: info, sup_recv
    INTEGER :: idrowref, idcolref, idref, idrecv 

    ! B and BUF_RECV are used to overload the computation of matrix multiplication and the shift of the blocks
    COMPLEX(DP), ALLOCATABLE, DIMENSION( :, : ) :: B, C, BUF_RECV 
    COMPLEX(DP) :: first

    myrow = desc( la_myr_ )
    mycol = desc( la_myc_ )
    myid  = desc( la_me_ )
    np    = desc( la_npr_ )
    comm  = desc( la_comm_ )

    IF( desc( la_npr_ ) /= desc( la_npc_ ) ) THEN
       CALL errore( ' pztrtri ', ' only square grid are allowed ', 1 ) 
    END IF

    nr = desc( nlar_ )
    nc = desc( nlac_ )

    DO j = nc+1, ldx
       DO i = 1, ldx
          sll( i, j ) = ( 0.0_DP , 0.0_DP )
       END DO
    END DO
    DO j = 1, ldx
       DO i = nr+1, ldx
          sll( i, j ) = ( 0.0_DP , 0.0_DP )
       END DO
    END DO

    ! Compute the inverse of a lower triangular 
    ! along the diagonal of the global array with BLACS(dtrtri) 
    IF( mycol == myrow ) THEN
       DO j = 1, ldx
          DO i = 1, j-1
             sll ( i, j ) = ( 0.0_DP , 0.0_DP )
          END DO
       END DO
       CALL ztrtri( 'L', 'N', nr, sll, ldx, info )
       IF( info /= 0 ) THEN
          CALL errore( ' pztrtri ', ' problem in the local inversion ', info )
       END IF
    ELSE
       buf_recv = sll
    END IF

#if defined __MPI

    ALLOCATE( B( ldx, ldx ) )
    ALLOCATE( C( ldx, ldx ) )
    ALLOCATE( BUF_RECV ( ldx, ldx ) )

    ! Broadcast the blocks along the diagonal at the processors under the diagonal
    IF( myrow >= mycol ) THEN
       CALL MPI_Comm_split( comm, mycol, myrow, col_comm, ierr )
       CALL MPI_Comm_size( col_comm, group_size, ierr )
       CALL MPI_Comm_rank( col_comm, group_rank, ierr )
       CALL MPI_Bcast( sll, ldx*ldx, MPI_DOUBLE_COMPLEX, 0, col_comm, ierr )
    ELSE
       CALL MPI_Comm_split( comm, MPI_UNDEFINED, MPI_UNDEFINED, col_comm, ierr )
       sll = ( 0.0_DP , 0.0_DP )
    END IF
 
    ! Compute A2 * C1 and start the Cannon's algorithm shifting the blocks of column one place to the North
    IF( myrow > mycol ) THEN
       CALL zgemm( 'N', 'N', nr, nc, nc, ONE, buf_recv, ldx, sll, ldx, ZERO, c, ldx )
       send = shift ( 1, group_rank, 1, ( group_size - 1 ), 'N' )
       recv = shift(1, group_rank, 1, (group_size-1), 'S')
       CALL MPI_Sendrecv( c, ldx*ldx, MPI_DOUBLE_COMPLEX, send, 0, buf_recv, &
            ldx*ldx, MPI_DOUBLE_COMPLEX, recv, 0, col_comm, status, ierr )
    END IF

    ! Execute the Cannon's algorithm to compute ricorsively the multiplication of C2 = -C3 * A2 * C1
    DO count = ( np - 2 ), 0, -1
       C3dim = (np-1) - count ! Dimension of the submatrix C3
       first = ZERO
       cicle = 0
       IF( ( myrow > count ) .AND. ( mycol >= count ) ) THEN
          idcolref = count + 1
          idrowref = myrow
          CALL GRID2D_RANK( 'R', np, np, idrowref, idcolref, idref )
          idrecv = idref - 1
          ! Compute C2 = -C3 * A2 * C1
          DO shiftcount = count, np-2
             IF(mycol>count)THEN
                ! Execute the virtual shift of the matrix C3 along the row in order to know which processor 
                ! have to send the block to C2 
                IF( cicle == 0)THEN
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   send = shift(idref, myid, myrow-count, C3dim, 'W')
                ELSE
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   send = shift(idref, send, 1, C3dim, 'E')
                END IF
                IF(send==idref)THEN
                   CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_COMPLEX, idrecv, myid, comm, ierr)
                END IF
             ELSE
                IF( cicle == 0)THEN
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   sup_recv = shift(idref, myid+1, myrow-count, C3dim, 'E')
                ELSE
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   sup_recv = shift(idref, sup_recv, 1, C3dim, 'W')
                END IF
                CALL MPI_Recv(C, ldx*ldx, MPI_DOUBLE_COMPLEX, sup_recv, sup_recv, comm, status, ierr)
                send = shift(1, group_rank, 1, (group_size-1), 'S')
                recv = shift(1, group_rank, 1, (group_size-1), 'N')
                ! with the no-blocking communication the computation and the shift of the column block are overapped  
                IF( MOD( cicle, 2 ) == 0 ) THEN
                   CALL MPI_Isend(BUF_RECV, ldx*ldx, MPI_DOUBLE_COMPLEX, send, group_rank+cicle, col_comm, req(1), ierr)
                   CALL MPI_Irecv(B, ldx*ldx, MPI_DOUBLE_COMPLEX, recv, recv+cicle, col_comm, req(2), ierr)
                   CALL zgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, BUF_RECV, ldx, first, sll, ldx)
                ELSE
                   CALL MPI_Isend(B, ldx*ldx, MPI_DOUBLE_COMPLEX, send, group_rank+cicle, col_comm, req(1), ierr)
                   CALL MPI_Irecv(BUF_RECV, ldx*ldx, MPI_DOUBLE_COMPLEX, recv, recv+cicle, col_comm, req(2), ierr)
                   CALL zgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, B, ldx, ONE, sll, ldx)
                END IF
                CALL MPI_Wait(req(1), status, ierr)
                CALL MPI_Wait(req(2), status, ierr)
             END IF
             cicle = cicle + 1
             first = ONE 
          END DO
       END IF
    END DO

    IF( myrow >= mycol ) THEN
       IF( myrow == mycol ) THEN
          DO j = 1, ldx
             DO i = 1, j-1 
                sll ( i, j ) = ( 0.0_DP , 0.0_DP )
             END DO
          END DO
       END IF
       CALL mpi_comm_free( col_comm, ierr )
    END IF

    DEALLOCATE(B)
    DEALLOCATE(C)
    DEALLOCATE(BUF_RECV)

#endif

  CONTAINS

  INTEGER FUNCTION shift ( idref, id, pos, size, dir )

    IMPLICIT NONE

    INTEGER :: idref, id, pos, size
    CHARACTER ( LEN = 1 ) :: dir

    IF( ( dir == 'E' ) .OR. ( dir == 'S' ) ) THEN
       shift = idref + MOD ( ( id - idref ) + pos, size )
    ELSE IF( ( dir == 'W' ) .OR. ( dir == 'N' ) ) THEN
       shift = idref + MOD ( ( id - idref ) - pos + size, size )
    ELSE
       shift = -1
    END IF

    RETURN

  END FUNCTION shift

  END SUBROUTINE pztrtri
