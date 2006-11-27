!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


       INTEGER FUNCTION owner_block(ig, nx, np)
!
!   INPUT :
!           ig  global index of the x dimension of array element
!           nx  x dimension of the global array
!           np  number of processor in the x dimension of the processors grid
!
!   OUTPUT :
!   This function return the index (starting from 0) of the
!   processor owning the eleent ig for the balanced block
!   distribution
!

       IMPLICIT NONE

       INTEGER :: ig, nx, np, r, q

       q = INT(nx/np)
       r = MOD(nx,np)

       IF ( ig .le. ((q+1)*r) ) then
         owner_block = INT( (ig-1) / (q+1) )
       else
         owner_block = INT( (ig-1-r*(q+1)) / q ) + r
       end if

       RETURN
       END FUNCTION owner_block

!=----------------------------------------------------------------------------=!

       INTEGER FUNCTION owner_cyclic(ig, nx, np)
!
!   INPUT :
!           ig  global index of the x dimension of array element
!           nx  x dimension of the global array
!           np  number of processor in the x dimension of the processors grid
!
!   OUTPUT :
!   This function return the index (starting from 0) of the
!   processor owning the eleent ig for the cyclic distribution
!

       IMPLICIT NONE

       INTEGER :: ig, nx, np

       owner_cyclic = MOD( ig-1, np )

       RETURN
       END FUNCTION owner_cyclic

!=----------------------------------------------------------------------------=!

      INTEGER FUNCTION owner_block_cyclic( INDXGLOB, NB, NPROCS )

!     Derived from: INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER       INDXGLOB, IPROC, ISRCPROC, NB, NPROCS, INDXG2P
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
      ISRCPROC = 0
      IPROC    = 0
      INDXG2P = MOD( ISRCPROC + (INDXGLOB - 1) / NB, NPROCS )
!
      owner_block_cyclic = INDXG2P
      RETURN
!
!     End of INDXG2P
!
      END FUNCTION owner_block_cyclic
