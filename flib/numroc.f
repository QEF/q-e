      INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER              IPROC, ISRCPROC, N, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  NUMROC computes the NUMber of Rows Or Columns of a distributed
!  matrix owned by the process indicated by IPROC.
!
!  Arguments
!  =========
!
!  N         (global input) INTEGER
!            The number of rows/columns in distributed matrix.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local input) INTEGER
!            The coordinate of the process whose local array row or
!            column is to be determined.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process that possesses the first
!            row or column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the matrix is
!            distributed.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER              EXTRABLKS, MYDIST, NBLOCKS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC            MOD
!     ..
!     .. Executable Statements ..
!
!     Figure PROC's distance from source process
!
      MYDIST = MOD( NPROCS+IPROC-ISRCPROC, NPROCS )
!
!     Figure the total number of whole NB blocks N is split up into
!
      NBLOCKS = N / NB
!
!     Figure the minimum number of rows/cols a process can have
!
      NUMROC = (NBLOCKS/NPROCS) * NB
!
!     See if there are any extra blocks
!
      EXTRABLKS = MOD( NBLOCKS, NPROCS )
!
!     If I have an extra block
!
      IF( MYDIST.LT.EXTRABLKS ) THEN
          NUMROC = NUMROC + NB
!
!         If I have last block, it may be a partial block
!
      ELSE IF( MYDIST.EQ.EXTRABLKS ) THEN
          NUMROC = NUMROC + MOD( N, NB )
      END IF
!
      RETURN
!
!     End of NUMROC
!
      END
