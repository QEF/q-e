      INTEGER FUNCTION NPREROC( N, NB, IPROC, ISRCPROC, NPROCS )
*
*  -- ScaLAPACK tools routine (version 1.5) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER              IPROC, ISRCPROC, N, NB, NPROCS
*     ..
*
*  Purpose
*  =======
*
*  NPREROC computes the Number of PREceeding Rows Or Columns of a
*  distributed matrix that are possessed by processes closer to
*  ISRCPROC than IPROC.  Therefore, if ISRCPROC=0 and IPROC=4, then
*  NPREROC returns the number of distributed matrix rows or columns
*  owned by processes 0, 1, 2, and 3.
*
*  Arguments
*  =========
*
*  N         (global input) INTEGER
*            The number of rows or columns in the distributed matrix.
*
*  NB        (global input) INTEGER
*            Block size, size of the blocks the distributed matrix is
*            split into.
*
*  IPROC     (local intput) INTEGER
*            The coordinate of the process whose preceeding distributed
*            matrix rows or columns are to be determined.
*
*  ISRCPROC  (global input) INTEGER
*            The coordinate of the process that possesses the first
*            row or column of the distributed matrix.
*
*  NPROCS    (global input) INTEGER
*            The total number processes over which the matrix is
*            distributed.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER              EXTRABLKS, MYDIST, NBLOCKS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC            MOD
*     ..
*     .. Executable Statements ..
*
*     Figure PROC's distance from source process
*
      MYDIST = MOD( NPROCS+IPROC-ISRCPROC, NPROCS )
*
*     Figure the total number of whole NB blocks N is split up into
*
      NBLOCKS = N / NB
*
*     Figure the minimum number of rows/cols previous processes could have
*
      NPREROC = (NBLOCKS/NPROCS) * NB * MYDIST
*
*     See if there are any extra blocks
*
      EXTRABLKS = MOD( NBLOCKS, NPROCS )
*
*     If I have an extra block, all processes in front of me got one too
*
      IF( MYDIST.LE.EXTRABLKS ) THEN
         NPREROC = NPREROC + NB*MYDIST
*
*     If I have don't have an extra block, add in extra blocks of
*     preceeding processes and the partial block, if it exists
*
      ELSE
         NPREROC = NPREROC + EXTRABLKS*NB + MOD( N, NB )
      END IF
*
      RETURN
*
*     End of NPREROC
*
      END
