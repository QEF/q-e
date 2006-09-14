      SUBROUTINE INFOG2L( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW,    &
     &                    MYCOL, LRINDX, LCINDX, RSRC, CSRC )
!
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            CSRC, GCINDX, GRINDX, LRINDX, LCINDX, MYCOL,   &
     &                   MYROW, NPCOL, NPROW, RSRC
!     ..
!     .. Array Arguments ..
      INTEGER            DESC( * )
!     ..
!
!  Purpose
!  =======
!
!  INFOG2L computes the starting local indexes LRINDX, LCINDX corres-
!  ponding to the distributed submatrix starting globally at the entry
!  pointed by GRINDX, GCINDX. This routine returns the coordinates in
!  the grid of the process owning the matrix entry of global indexes
!  GRINDX, GCINDX, namely RSRC and CSRC.
!
!  Notes
!  =====
!
!  Each global data object is described by an associated description
!  vector.  This vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  Let A be a generic term for any 2D block cyclicly distributed array.
!  Such a global array has an associated description vector DESCA.
!  In the following comments, the character _ should be read as
!  "of the global array".
!
!  NOTATION        STORED IN      EXPLANATION
!  --------------- -------------- --------------------------------------
!  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
!                                 DTYPE_A = 1.
!  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
!                                 the BLACS process grid A is distribu-
!                                 ted over. The context itself is glo-
!                                 bal, but the handle (the integer
!                                 value) may vary.
!  M_A    (global) DESCA( M_ )    The number of rows in the global
!                                 array A.
!  N_A    (global) DESCA( N_ )    The number of columns in the global
!                                 array A.
!  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
!                                 the rows of the array.
!  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
!                                 the columns of the array.
!  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
!                                 row of the array A is distributed.
!  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
!                                 first column of the array A is
!                                 distributed.
!  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
!                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
!
!  Let K be the number of rows or columns of a distributed matrix,
!  and assume that its process grid has dimension p x q.
!  LOCr( K ) denotes the number of elements of K that a process
!  would receive if K were distributed over the p processes of its
!  process column.
!  Similarly, LOCc( K ) denotes the number of elements of K that a
!  process would receive if K were distributed over the q processes of
!  its process row.
!  The values of LOCr() and LOCc() may be determined via a call to the
!  ScaLAPACK tool function, NUMROC:
!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
!  An upper bound for these quantities may be computed by:
!          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
!          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
!
!  Arguments
!  =========
!
!  GRINDX    (global input) INTEGER
!            The global row starting index of the submatrix.
!
!  GCINDX    (global input) INTEGER
!            The global column starting index of the submatrix.
!
!  DESC      (input) INTEGER array of dimension DLEN_.
!            The array descriptor for the underlying distributed matrix.
!
!  NPROW     (global input) INTEGER
!            The total number of process rows over which the distributed
!            matrix is distributed.
!
!  NPCOL     (global input) INTEGER
!            The total number of process columns over which the
!            distributed matrix is distributed.
!
!  MYROW     (local input) INTEGER
!            The row coordinate of the process calling this routine.
!
!  MYCOL     (local input) INTEGER
!            The column coordinate of the process calling this routine.
!
!  LRINDX    (local output) INTEGER
!            The local rows starting index of the submatrix.
!
!  LCINDX    (local output) INTEGER
!            The local columns starting index of the submatrix.
!
!  RSRC      (global output) INTEGER
!            The row coordinate of the process that possesses the first
!            row and column of the submatrix.
!
!  CSRC      (global output) INTEGER
!            The column coordinate of the process that possesses the
!            first row and column of the submatrix.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
     &                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
     &                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
     &                      RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!     ..
!     .. Local Scalars ..
      INTEGER            CBLK, GCCPY, GRCPY, RBLK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
      GRCPY = GRINDX-1
      GCCPY = GCINDX-1
!
      RBLK = GRCPY / DESC(MB_)
      CBLK = GCCPY / DESC(NB_)
      RSRC = MOD( RBLK + DESC(RSRC_), NPROW )
      CSRC = MOD( CBLK + DESC(CSRC_), NPCOL )
!
      LRINDX = ( RBLK / NPROW + 1 ) * DESC(MB_) + 1
      LCINDX = ( CBLK / NPCOL + 1 ) * DESC(NB_) + 1
!
      IF( MOD( MYROW+NPROW-DESC(RSRC_), NPROW ) .GE.                    &
     &    MOD( RBLK, NPROW ) ) THEN
         IF( MYROW.EQ.RSRC )                                            &
     &      LRINDX = LRINDX + MOD( GRCPY, DESC(MB_) )
         LRINDX = LRINDX - DESC(MB_)
      END IF
!
      IF( MOD( MYCOL+NPCOL-DESC(CSRC_), NPCOL ) .GE.                    &
     &    MOD( CBLK, NPCOL ) ) THEN
         IF( MYCOL.EQ.CSRC )                                            &
     &      LCINDX = LCINDX + MOD( GCCPY, DESC(NB_) )
         LCINDX = LCINDX - DESC(NB_)
      END IF
!
      RETURN
!
!     End of INFOG2L
!
      END
