!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE scala_cdiag (n, a, ilda, w, z, ildz)
  !----------------------------------------------------------------------------
  !
  !     Solve the standard diagonalization problem:Ax=lambdax for
  !     hermitian comples matrix using Blacs and
  !     Scalapack routines.
  !
  !
  !     n          integer (input)
  !                Size of the problem.
  !     a(ilda,n)  complex*16  (input)
  !                Matrix to be diagonalized ( identical on each processor
  !     ilda       integer (input)
  !                Leading dimension of matrix a.
  !     w(n)       real*8  (output)
  !                Eigenvalues.
  !     z(ildz,n)  complex*16  (output)
  !                Eigenvectors ( identical on each processor ).
  !     ildz       integer (input)
  !                Leading dimension of matrix z.
  !     ..
  !
  USE kinds, ONLY : DP
  USE io_global,  ONLY :  stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: n, ilda, ildz
    ! input: size of the matrix to diagonalize
  REAL (DP) :: w (n)
   ! output: eigenvalues
  COMPLEX (DP) :: a (ilda, n), z (ildz, n)
   ! input: matrix to diagonalize
   ! output: eigenvectors
  !
#if defined (__T3E)
  !  
  REAL (DP) :: zero, mone
    ! zero
    ! minus one

  PARAMETER (zero = 0.0d+0, mone = - 1.d0)

  INTEGER :: maxprocs, maxn
    ! maximum number of processors
    ! maximum dimension of the matrix

  PARAMETER (maxn = 2000, maxprocs = 256)


  INTEGER :: context, i, j, info, npcol, nprow, mycol, myrow, &
       rowsize, colsize, lda, nb, idum, iam_blacs, nprocs_blacs, lwork, &
       lrwork, liwork, m, nz, desca (50), descz (50), iclustr (maxprocs * &
       2), ifail (maxn)
  ! the context for the matrix division
  ! counters
  ! integer of error checking
  ! number of row and columns of the proc. gri
  ! row and column of this processor
  ! size of the local matrices
  ! leading dimension of the local matrix
  ! blocking factor for matrix division
  ! used for input of the diagonalization rout
  ! the processor me
  ! number of processors
  ! sizes of the work space
  ! number of computed eigenvalues and eigenve
  ! description of matrix divisions
  !
  ! clustering of eigenvalues
  ! on output tells which eigenvalues not conv

  REAL (DP) :: gap (maxprocs), abstol
  ! gap between eigenvalues
  ! absolute tolerance (see man pcheevx)
  ! local matrices of this proc
  COMPLEX (DP), ALLOCATABLE::at (:, :)

  COMPLEX (DP), ALLOCATABLE::zt (:, :)
  ! some variables needed to determine the workspace's size

  INTEGER :: nn, np, np0, npo, nq0, mq0
  COMPLEX (DP), ALLOCATABLE::work (:)
  ! work space
  REAL (DP), ALLOCATABLE::rwork (:)
  INTEGER, ALLOCATABLE::iwork (:)
  INTEGER :: pslamch, iceil, numroc

  EXTERNAL pslamch, iceil, numroc

  INTEGER :: iws

  IF (n.GT.maxn) CALL errore ('scala_cdiag', 'n is too large', n)
  !
  !   Ask for the number of processors
  !
  CALL blacs_pinfo (iam_blacs, nprocs_blacs)

  IF (nprocs_blacs.GT.32) CALL errore ('scala_cdiag', &
       'work space not  optimized ',  - nprocs_blacs)
  !
  !     Set the dimension of the 2D processors grid
  !

  CALL gridsetup_local (nprocs_blacs, nprow, npcol)
  !
  !     Initialize a single total BLACS context
  !
  CALL blacs_get (0, 0, context)
  CALL blacs_gridinit (context, 'r', nprow, npcol)

  CALL blacs_gridinfo (context, nprow, npcol, myrow, mycol)
  !
  !     Calculate the blocking factor for the matrix
  !
  CALL blockset_priv (nb, 32, n, nprow, npcol)
  !
  !     These are basic array descriptors
  !
  rowsize = n / nprow + nb + 1
  colsize = n / npcol + nb + 1
  lda = rowsize
  !
  !     allocate local array
  !
  ALLOCATE (at (rowsize, colsize) )
  ALLOCATE (zt (rowsize, colsize) )
  !
  !     initialize the data structure of description of the matrix partiti
  !
  CALL descinit (desca, n, n, nb, nb, 0, 0, context, lda, info)

  IF (info.NE.0) CALL errore ('scala_cdiag', &
       'something wrong with descinit1', info)
  CALL descinit (descz, n, n, nb, nb, 0, 0, context, lda, info)

  IF (info.NE.0) CALL errore ('scala_cdiag', &
       'something wrong with descinit2', info)
  !
  !   Now try to determine the size of the workspace
  !
  !
  !  NB: The values of lwork, lrwork, liwork are important for the
  !      performance of the routine,
  !      they seem appropriate for  150 < n < 2000. The routine should
  !      work also for n < 150, but it is slower than the scalar one
  !      The sizes could be adjusted if something better is found
  !      only 4 < nprocs_blacs =< 32 has been tested
  !
  iws = 0

  IF (iws.EQ.0) THEN
     lwork = 40 * n / (nprocs_blacs / 16.d0)
     IF (n.GT.1000) THEN
        lrwork = 300 * n / (nprocs_blacs / 16.d0)
        liwork = 8 * n
     ELSE
        liwork = 10 * n
        lrwork = 300 * n / (nprocs_blacs / 16.d0)
     ENDIF
     lwork = MAX (10000, lwork)
     lrwork = MAX (2000, lrwork)
     liwork = MAX (2000, liwork)
  ENDIF
  !
  !  following the notes included in the man page
  !
  IF (iws.EQ.1) THEN
     nn = MAX (n, nb, 2)
     npo = numroc (nn, nb, 0, 0, nprow)
     nq0 = MAX (numroc (n, nb, 0, 0, npcol), nb)

     mq0 = numroc (MAX (n, nb, 2), nb, 0, 0, npcol)


     lrwork = 4 * n + MAX (5 * nn, npo * mq0) + iceil (n, nprow * &
          npcol) * nn
     ! note: the following few lines from the man page are wrong: the right s
     ! the other way around !!!!!!!
     !
     !     lwork   Integer.  (locak input)
     !             Size of work array.  If only eigenvalues are requested, lw
     !             >= N + (NPO + MQP + NB) * NB.  If eigenvectors are request
     !             lwork >= N + MAX(NB*(NPO+1),3)
     !
     !

     lwork = n + nb * (npo + mq0 + nb)
     liwork = 6 * MAX (n, nprow * npcol + 1, 4)

  ENDIF

  IF (iws.EQ.2) THEN
     ! the first way we did ( just to compare)
     nb = desca (6)
     nn = MAX (n, nb, 2)

     np = numroc (n, nb, myrow, 0, nprow)

     liwork = 6 * MAX (n, nprow * npcol + 1, 4)
     npo = numroc (nn, nb, 0, 0, nprow)
     nq0 = MAX (numroc (n, nb, 0, 0, npcol), nb)

     mq0 = numroc (MAX (n, nb, 2), nb, 0, 0, npcol)
     lrwork = 4 * n + MAX (5 * nn, npo * mq0) + iceil (n, nprow * &
          npcol) * nn

     lwork = n + MAX (nb * (npo + 1), 3)

     lwork = 3 * (MAX (liwork, lrwork, 2 * lwork) + 3) + 200000
     lwork = lwork / 3
     lrwork = lwork

     liwork = lwork
  ENDIF
  !
  !  and allocate sufficient work space:
  !
  ALLOCATE (work (lwork) )
  ALLOCATE (rwork (lrwork) )

  ALLOCATE (iwork (liwork) )
  !
  !     copy the elements on the local matrices
  !
  DO j = 1, n
     DO i = 1, n
        CALL pcelset (at, i, j, desca, a (i, j) )
     ENDDO

  ENDDO

  abstol = pslamch (desca (2) , 's')
  CALL pcheevx ('v', 'a', 'u', n, at, 1, 1, desca, zero, zero, idum, &
       idum, abstol, m, nz, w, 1.0d-5, zt, 1, 1, descz, work, lwork, &
       rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
  !
  IF (ABS (info) .GT.2) THEN
     WRITE( stdout, * ) 'info ', info, m, nz
     CALL errore ('scala_cdiag', 'wrong info', 1)
  ENDIF
  !
  !
  !    compute the eigenvalues
  !

  CALL eigen (n, z, ildz, zt, descz, work)
  DEALLOCATE (at)
  DEALLOCATE (zt)
  DEALLOCATE (iwork)
  DEALLOCATE (rwork)
  DEALLOCATE (work)
  !
  CALL blacs_gridexit (context)
  !
#endif
  RETURN

END SUBROUTINE scala_cdiag

