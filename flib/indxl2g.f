      INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            INDXLOC, IPROC, ISRCPROC, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  INDXL2G computes the global index of a distributed matrix entry
!  pointed to by the local index INDXLOC of the process indicated by
!  IPROC.
!
!  Arguments
!  =========
!
!  INDXLOC   (global input) INTEGER
!            The local index of the distributed matrix entry.
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
!            row/column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the distributed
!            matrix is distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC            MOD
!     ..
!     .. Executable Statements ..
!
      INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) +        &
     &          MOD(NPROCS+IPROC-ISRCPROC, NPROCS)*NB + 1
!
      RETURN
!
!     End of INDXL2G
!
      END
