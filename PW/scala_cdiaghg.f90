!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE scala_cdiaghg (n, a, ilda, b, ildb, w, z, ildz)
  !----------------------------------------------------------------------------
  !
  !     Solves the generalized problem Ax=lambdaBx for
  !     hermitian comples matrix using Blacs and
  !     Scalapack routines.
  !
  !     n          integer (input)
  !                Size of the problem.
  !     a(ilda,n)  complex*16  (input/ouput)
  !                Matrix to be diagonalized ( identical on each processor
  !     ilda       integer (input)
  !                Leading dimension of matrix a.
  !     b(ildb,n)  complex*16  (input/ouput)
  !                Overlap matrix  ( identical on each processor ).
  !     ildb       integer (input)
  !                Leading dimension of matrix b.
  !
  !     w(n)       real*8  (output)
  !                Eigenvalues.
  !     z(ildz,n)  complex*16 (output)
  !                Matrix containing eigenvectors (on each processor)
  !     iopt       integer (input)
  !                if iopt > 0  the routine calculate both eigenvalues and
  !                eigenvectors, if iopt <= 0 only the eigenvalues are
  !                computed.
  !
  ! This  driver work just if npool=1
  !
  USE kinds, ONLY : DP
  USE io_global,  ONLY :  stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: n, ilda, ildb, ildz
  ! input: size of the matrix to diagonalize
  ! input: see before

  REAL (kind=DP) :: w (n)
  ! output: eigenvalues

  COMPLEX (kind=DP) :: a (ilda, n), b (ildb, n), z (ildz, n)
  ! input: matrix to diagonalize
  ! input: overlap matrix
  ! output: eigenvectors
  !
#if defined (__T3E)
#  if defined (__AIX)
#define PCELSET pzelset
#define PCHEGVX pzhegvx
#endif
  !
  REAL (kind=DP) :: zero, mone
  ! zero
  ! minus one

  PARAMETER (zero = 0.0d+0, mone = - 1.d0)

  INTEGER :: maxprocs, maxn
  ! maximum number of processors
  ! maximum dimension of the matrix

  PARAMETER (maxn = 2000, maxprocs = 256)

  INTEGER :: context, i, j, info, npcol, nprow, mycol, myrow, &
       rowsize, colsize, lda, nb, idum, iam_blacs, nprocs_blacs, lwork, &
       lrwork, liwork, m, nz, desca (50), descb (50), descz (50), &
       iclustr (maxprocs * 2), ifail (maxn)
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
  !
  !  description of matrix divisions
  !
  ! clustering of eigenvalues
  ! on output tells which eigenvalues not conv

  REAL (kind=DP) :: gap (maxprocs)
  ! gap between eigenvalues
  COMPLEX (kind=DP), ALLOCATABLE::at (:, :)
  ! local matrices of this p
  COMPLEX (kind=DP), ALLOCATABLE::bt (:, :)

  COMPLEX (kind=DP), ALLOCATABLE::zt (:, :)
  ! some variables needed to determine the workspace's size

  INTEGER :: nn, np, npO, nq0, mq0
  COMPLEX (kind=DP), ALLOCATABLE::work (:)
  ! work space
  REAL (kind=DP), ALLOCATABLE::rwork (:)
  INTEGER, ALLOCATABLE::iwork (:)
  INTEGER :: PSLAMCH, iceil, numroc

  EXTERNAL PSLAMCH, iceil, numroc
  INTEGER :: iws, first

  REAL (kind=DP) :: time1, time2, tt
  IF (n.GT.maxn) CALL errore ('scala_cdiaghg', 'n is too large', n)
  !
  !   Ask for the number of processors
  !
  CALL blacs_pinfo (iam_blacs, nprocs_blacs)
  IF (nprocs_blacs.GT.32) CALL errore ('scala_cdiaghg', &
       'work space not optimized ',  - nprocs_blacs)
  !
  !     Set the dimension of the 2D processors grid
  !
  CALL gridsetup_local (nprocs_blacs, nprow, npcol)
  !
  !     Initialize a single total BLACS context
  !
  CALL blacs_get ( - 1, 0, context)
  CALL blacs_gridinit (context, 'R', nprow, npcol)
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
  ALLOCATE (bt (rowsize, colsize) )
  ALLOCATE (zt (rowsize, colsize) )
  !
  !     initialize the data structure of description of the matrix partiti
  !
  CALL descinit (desca, n, n, nb, nb, 0, 0, context, lda, info)
  IF (info.NE.0) CALL errore ('scala_cdiaghg', &
       'something wrong with descinit1', info)
  CALL descinit (descb, n, n, nb, nb, 0, 0, context, lda, info)
  IF (info.NE.0) CALL errore ('scala_cdiaghg', &
       'something wrong with descinit2', info)
  CALL descinit (descz, n, n, nb, nb, 0, 0, context, lda, info)
  IF (info.NE.0) CALL errore ('scala_cdiaghg', &
       'something wrong with descinit3', info)
  !
  !   allocate  workspace needed
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
  !   allocate sufficient work space:
  !
  ALLOCATE (work (lwork) )
  ALLOCATE (rwork (lrwork) )

  ALLOCATE (iwork (liwork) )
  !
  !     copy the elements on the local matrices
  !
  DO i = 1, n
     DO j = i, n
        CALL PCELSET (at, i, j, desca, a (i, j) )
        CALL PCELSET (bt, i, j, descb, b (i, j) )
     ENDDO
  ENDDO
  !
  !     compute work parameters and allocate workspace
  !

  CALL PCHEGVX (1, 'V', 'A', 'U', n, at, 1, 1, desca, bt, 1, 1, &
       descb, zero, zero, idum, idum, mone, m, nz, w, mone, zt, 1, 1, &
       descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, &
       gap, info)
  !
  !   check the info variable to detect problems
  !
#ifdef DEBUG
  IF (info.NE.0) THEN


     IF (iam_blacs.EQ.0) THEN
        IF (info.LT.0) THEN
           IF (info> - 100) THEN
              WRITE( stdout,  * ) 'scala_cdiaghg: Argument',  - info, &
                   'to PCHEGVX had an illegal value'

           ENDIF
           IF (info< - 100) THEN
              i = - info / 100
              j = MOD ( - info, 100)
              WRITE( stdout,  * ) 'scala_cdiagh: Element', j, 'of argument', i, &
                   'to PCHEGVX had an illegal value'
           ENDIF
        ENDIF
        WRITE( stdout, * ) 'given and requested lwork', lwork, work (1)
        WRITE( stdout, * ) 'given and requested lrwork', lrwork, rwork (1)


        WRITE( stdout, * ) 'given and requested liwork', liwork, iwork (1)
        IF (info.GT.0) THEN
           IF (MOD (info, 2) .NE.0) THEN
              WRITE( stdout,  * ) 'scala_cdiaghg: PCHEGVX: Calculation failed', &
                   ' to converge'

           ENDIF
           IF (MOD (info / 2, 2) .NE.0) THEN
              WRITE( stdout,  * ) 'scala_cdiaghg: PCHEGVX: Insufficient workspace', &
                   ' to orthogonalize eigenvectors'

           ENDIF
           IF (MOD (info / 4, 2) .NE.0) THEN
              WRITE( stdout,  * ) 'scala_cdiaghg: PCHEGVX: Insufficient workspace', &
                   ' to compute all eigenvectors'
           ENDIF
        ENDIF
     ENDIF
  ENDIF
#endif
  !
  !    compute the eigenvalues
  !

  CALL eigen (n, z, ildz, zt, descz, work)
  DEALLOCATE (at)
  DEALLOCATE (bt)
  DEALLOCATE (zt)
  DEALLOCATE (iwork)
  DEALLOCATE (rwork)
  DEALLOCATE (work)
  !
  CALL blacs_gridexit (context)
  !
#endif
  RETURN
END SUBROUTINE scala_cdiaghg

