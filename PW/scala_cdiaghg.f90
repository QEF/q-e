!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine scala_cdiaghg (n, a, ilda, b, ildb, w, z, ildz)
#ifdef T3D
#ifdef AIX
#define PCELSET pzelset
#define PCHEGVX pzhegvx
#endif
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
  use parameters, only : DP
  implicit none
  integer :: n, ilda, ildb, ildz
  ! input: size of the matrix to diagonalize
  ! input: see before

  real (kind=DP) :: w (n)
  ! output: eigenvalues

  complex (kind=DP) :: a (ilda, n), b (ildb, n), z (ildz, n)
  ! input: matrix to diagonalize
  ! input: overlap matrix
  ! output: eigenvectors
  real (kind=DP) :: zero, mone
  ! zero
  ! minus one

  parameter (zero = 0.0d+0, mone = - 1.d0)

  integer :: maxprocs, maxn
  ! maximum number of processors
  ! maximum dimension of the matrix

  parameter (maxn = 2000, maxprocs = 256)

  integer :: context, i, j, info, npcol, nprow, mycol, myrow, &
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

  real (kind=DP) :: gap (maxprocs)
  ! gap between eigenvalues
  complex (kind=DP), allocatable::at (:, :)
  ! local matrices of this p
  complex (kind=DP), allocatable::bt (:, :)

  complex (kind=DP), allocatable::zt (:, :)
  ! some variables needed to determine the workspace's size

  integer :: nn, np, npO, nq0, mq0
  complex (kind=DP), allocatable::work (:)
  ! work space
  real (kind=DP), allocatable::rwork (:)
  integer, allocatable::iwork (:)
  integer :: PSLAMCH, iceil, numroc

  external PSLAMCH, iceil, numroc
  integer :: iws, first

  real (kind=DP) :: time1, time2, tt
  if (n.gt.maxn) call error ('scala_cdiaghg', 'n is too large', n)
  !
  !   Ask for the number of processors
  !
  call blacs_pinfo (iam_blacs, nprocs_blacs)
  if (nprocs_blacs.gt.32) call error ('scala_cdiaghg', &
       'work space not optimized ',  - nprocs_blacs)
  !
  !     Set the dimension of the 2D processors grid
  !
  call gridsetup_local (nprocs_blacs, nprow, npcol)
  !
  !     Initialize a single total BLACS context
  !
  call blacs_get ( - 1, 0, context)
  call blacs_gridinit (context, 'R', nprow, npcol)
  call blacs_gridinfo (context, nprow, npcol, myrow, mycol)
  !
  !     Calculate the blocking factor for the matrix
  !
  call blockset_priv (nb, 32, n, nprow, npcol)
  !
  !     These are basic array descriptors
  !
  rowsize = n / nprow + nb + 1
  colsize = n / npcol + nb + 1
  lda = rowsize
  !
  !     allocate local array
  !
  allocate (at (rowsize, colsize) )
  allocate (bt (rowsize, colsize) )
  allocate (zt (rowsize, colsize) )
  !
  !     initialize the data structure of description of the matrix partiti
  !
  call descinit (desca, n, n, nb, nb, 0, 0, context, lda, info)
  if (info.ne.0) call error ('scala_cdiaghg', &
       'something wrong with descinit1', info)
  call descinit (descb, n, n, nb, nb, 0, 0, context, lda, info)
  if (info.ne.0) call error ('scala_cdiaghg', &
       'something wrong with descinit2', info)
  call descinit (descz, n, n, nb, nb, 0, 0, context, lda, info)
  if (info.ne.0) call error ('scala_cdiaghg', &
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

  if (iws.eq.0) then
     lwork = 40 * n / (nprocs_blacs / 16.d0)
     if (n.gt.1000) then
        lrwork = 300 * n / (nprocs_blacs / 16.d0)
        liwork = 8 * n
     else
        liwork = 10 * n
        lrwork = 300 * n / (nprocs_blacs / 16.d0)
     endif
     lwork = max (10000, lwork)
     lrwork = max (2000, lrwork)
     liwork = max (2000, liwork)
  endif
  !
  !  following the notes included in the man page
  !

  if (iws.eq.1) then
     nn = max (n, nb, 2)
     npo = numroc (nn, nb, 0, 0, nprow)
     nq0 = max (numroc (n, nb, 0, 0, npcol), nb)

     mq0 = numroc (max (n, nb, 2), nb, 0, 0, npcol)



     lrwork = 4 * n + max (5 * nn, npo * mq0) + iceil (n, nprow * &
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
     liwork = 6 * max (n, nprow * npcol + 1, 4)

  endif

  if (iws.eq.2) then
     ! the first way we did ( just to compare)
     nb = desca (6)
     nn = max (n, nb, 2)

     np = numroc (n, nb, myrow, 0, nprow)

     liwork = 6 * max (n, nprow * npcol + 1, 4)
     npo = numroc (nn, nb, 0, 0, nprow)
     nq0 = max (numroc (n, nb, 0, 0, npcol), nb)

     mq0 = numroc (max (n, nb, 2), nb, 0, 0, npcol)
     lrwork = 4 * n + max (5 * nn, npo * mq0) + iceil (n, nprow * &
          npcol) * nn

     lwork = n + max (nb * (npo + 1), 3)

     lwork = 3 * (max (liwork, lrwork, 2 * lwork) + 3) + 200000
     lwork = lwork / 3
     lrwork = lwork

     liwork = lwork


  endif
  !
  !   allocate sufficient work space:
  !
  allocate (work (lwork) )
  allocate (rwork (lrwork) )

  allocate (iwork (liwork) )
  !
  !     copy the elements on the local matrices
  !
  do i = 1, n
     do j = i, n
        call PCELSET (at, i, j, desca, a (i, j) )
        call PCELSET (bt, i, j, descb, b (i, j) )
     enddo
  enddo
  !
  !     compute work parameters and allocate workspace
  !

  call PCHEGVX (1, 'V', 'A', 'U', n, at, 1, 1, desca, bt, 1, 1, &
       descb, zero, zero, idum, idum, mone, m, nz, w, mone, zt, 1, 1, &
       descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, &
       gap, info)
  !
  !   check the info variable to detect problems
  !
#ifdef DEBUG
  if (info.ne.0) then


     if (iam_blacs.eq.0) then
        if (info.lt.0) then
           if (info> - 100) then
              write (6,  * ) 'scala_cdiaghg: Argument',  - info, &
                   'to PCHEGVX had an illegal value'

           endif
           if (info< - 100) then
              i = - info / 100
              j = mod ( - info, 100)
              write (6,  * ) 'scala_cdiagh: Element', j, 'of argument', i, &
                   'to PCHEGVX had an illegal value'
           endif
        endif
        write (6, * ) 'given and requested lwork', lwork, work (1)
        write (6, * ) 'given and requested lrwork', lrwork, rwork (1)


        write (6, * ) 'given and requested liwork', liwork, iwork (1)
        if (info.gt.0) then
           if (mod (info, 2) .ne.0) then
              write (6,  * ) 'scala_cdiaghg: PCHEGVX: Calculation failed', &
                   ' to converge'

           endif
           if (mod (info / 2, 2) .ne.0) then
              write (6,  * ) 'scala_cdiaghg: PCHEGVX: Insufficient workspace', &
                   ' to orthogonalize eigenvectors'

           endif
           if (mod (info / 4, 2) .ne.0) then
              write (6,  * ) 'scala_cdiaghg: PCHEGVX: Insufficient workspace', &
                   ' to compute all eigenvectors'
           endif
        endif
     endif
  endif
#endif
  !
  !    compute the eigenvalues
  !

  call eigen (n, z, ildz, zt, descz, work)
  deallocate (at)
  deallocate (bt)
  deallocate (zt)
  deallocate (iwork)
  deallocate (rwork)
  deallocate (work)
  !
  call blacs_gridexit (context)
  !
#endif
  return
end subroutine scala_cdiaghg

