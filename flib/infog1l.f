      SUBROUTINE INFOG1L( GINDX, NB, NPROCS, MYROC, ISRCPROC, LINDX,    &
     &                    ROCSRC )
!
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            GINDX, ISRCPROC, LINDX, MYROC, NB, NPROCS,     &
     &                   ROCSRC
!     ..
!
!  Purpose
!  =======
!
!  INFOG1L computes the starting local indexes LINDX corresponding to
!  the distributed submatrix starting globally at the entry pointed by
!  GINDX.  This routine returns the coordinates of the process in the
!  grid owning the submatrix entry of global index GINDX: ROCSRC.
!  INFOG1L is a 1-dimensional version of INFOG2L.
!
!  Arguments
!  =========
!
!  GINDX     (global input) INTEGER
!            The global starting index of the submatrix.
!
!  NB        (global input) INTEGER
!            The block size.
!
!  NPROCS    (global input) INTEGER
!            The total number of processes over which the distributed
!            submatrix is distributed.
!
!  MYROC     (local input) INTEGER
!            The coordinate of the process calling this routine.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process having the first entry of
!            the distributed submatrix.
!
!  LINDX     (local output) INTEGER
!            The local starting indexes of the distributed submatrix.
!
!  ROCSRC    (global output) INTEGER
!            The coordinate of the process that possesses the first
!            row and column of the submatrix.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            GCPY, IBLK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
      GCPY = GINDX-1
      IBLK = GCPY / NB
      ROCSRC = MOD( IBLK + ISRCPROC, NPROCS )
!
      LINDX = ( IBLK / NPROCS + 1 ) * NB + 1
!
      IF( MOD(MYROC+NPROCS-ISRCPROC,NPROCS).GE.MOD(IBLK, NPROCS) ) THEN
         IF( MYROC.EQ.ROCSRC ) THEN
            LINDX = LINDX + MOD( GCPY, NB )
         END IF
         LINDX = LINDX - NB
      END IF
!
      RETURN
!
!     End of INFOG1L
!
      END
