!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


       INTEGER FUNCTION ldim_cyclic(gdim, np, me)

!   gdim = global dimension of distributed array
!   np   = number of processors
!   me   = index of the calling processor (starting from 0)
!  
!   this function return the number of elements of the distributed array
!   stored in the local memory of the processor "me" for a cyclic 
!   data distribution.
!   Example of the cyclic distribution of a 10 elements array on 4 processors
!   array elements  |  PEs
!    a(1)           |   0
!    a(2)           |   1
!    a(3)           |   2
!    a(4)           |   3
!    a(5)           |   0
!    a(6)           |   1
!    a(7)           |   2
!    a(8)           |   3
!    a(9)           |   0
!    a(10)          |   1

       IMPLICIT NONE
       INTEGER :: gdim, np, me, r, q

       IF( me >= np .OR. me < 0 ) THEN
         WRITE(6,*) ' ** ldim_cyclic: arg no. 3 out of range '
         STOP
       END IF

       q = INT(gdim / np)
       r = MOD(gdim, np)

       IF( me .LT. r ) THEN
         ldim_cyclic = q+1
       ELSE
         ldim_cyclic = q
       END IF
 
       RETURN
       END FUNCTION ldim_cyclic

!=----------------------------------------------------------------------------=!

       INTEGER FUNCTION ldim_block(gdim, np, me)

!   gdim = global dimension of distributed array
!   np   = number of processors
!   me   = index of the calling processor (starting from 0)
!  
!   this function return the number of elements of the distributed array
!   stored in the local memory of the processor "me" for a balanced block 
!   data distribution, with the larger block on the lower index processors.
!   Example of the block distribution of 10 elements array a on 4 processors
!   array elements  |  PEs
!    a(1)           |   0
!    a(2)           |   0
!    a(3)           |   0
!    a(4)           |   1
!    a(5)           |   1
!    a(6)           |   1
!    a(7)           |   2
!    a(8)           |   2
!    a(9)           |   3
!    a(10)          |   3

       IMPLICIT NONE
       INTEGER :: gdim, np, me, r, q

       IF( me >= np .OR. me < 0 ) THEN
         WRITE(6,*) ' ** ldim_block: arg no. 3 out of range '
         STOP
       END IF

       q = INT(gdim / np)
       r = MOD(gdim, np)

       IF( me .LT. r ) THEN
! ...    if my index is less than the reminder I got an extra element
         ldim_block = q+1
       ELSE
         ldim_block = q
       END IF
 
       RETURN
       END FUNCTION ldim_block



      INTEGER FUNCTION ldim_block_cyclic( N, NB, NPROCS, IPROC )

!  -- Derived from:  NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER              IPROC, ISRCPROC, N, NB, NPROCS, NUMROC
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
      ISRCPROC = 0
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

      ldim_block_cyclic = numroc
      RETURN
!
!     End of NUMROC
!
      END FUNCTION ldim_block_cyclic
