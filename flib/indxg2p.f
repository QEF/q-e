      INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  INDXG2P computes the process coordinate which posseses the entry of a
!  distributed matrix specified by a global index INDXGLOB.
!
!  Arguments
!  =========
!
!  INDXGLOB  (global input) INTEGER
!            The global index of the element.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process that possesses the first
!            row/column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the matrix is
!            distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
      INDXG2P = MOD( ISRCPROC + (INDXGLOB - 1) / NB, NPROCS )
!
      RETURN
!
!     End of INDXG2P
!
      END
