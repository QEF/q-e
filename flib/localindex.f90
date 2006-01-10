!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


       INTEGER FUNCTION lind_block(ig, nx, np, me)
!
!   INPUT :
!      ig  global index of the x dimension of array element
!      nx  dimension of the global array
!      np  number of processor in the x dimension of the processors grid
!      me  index of the local processor in the processor grid
!                (starting from zero)
!
!   OUTPUT :
!
!      lind_block return the local index corresponding to the
!      global index "ig" for a balanced block distribution
!   

       IMPLICIT NONE

       INTEGER :: ig, nx, np, me, r, q

       q = INT(nx/np)
       r = MOD(nx,np) 

       IF( me < r ) THEN
         lind_block = ig - (q+1) * me
       ELSE
         lind_block = ig - (q+1) * r - q * (me - r)
       END IF

       RETURN 
       END FUNCTION lind_block


       INTEGER FUNCTION lind_cyclic(ig, nx, np, me)
!
!   INPUT :
!      ig  global index of the x dimension of array element
!      nx  dimension of the global array
!      np  number of processor in the x dimension of the processors grid
!      me  index of the local processor in the processor grid
!                (starting from zero)
!
!   OUTPUT :
!
!      lind_cyclic return the local index corresponding to the
!      global index "ig" for a cyclic distribution
!   

       IMPLICIT NONE

       INTEGER :: ig, nx, np, me

       lind_cyclic = (ig-1)/np + 1

       RETURN 
       END FUNCTION lind_cyclic



       INTEGER FUNCTION lind_block_cyclic( INDXGLOB, NB, NPROCS, IPROC )      

!     Derived from:  INDXG2L( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      IMPLICIT NONE
      INTEGER  INDXGLOB, IPROC, ISRCPROC, NB, NPROCS, INDXG2L
!     ..
!
!  Purpose
!  =======
!
!  INDXG2L computes the local index of a distributed matrix entry
!  pointed to by the global index INDXGLOB.
!
!  Arguments
!  =========
!
!  INDXGLOB  (global input) INTEGER
!            The global index of the distributed matrix entry.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  ISRCPROC  (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the distributed
!            matrix is distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
      ISRCPROC = 0      
      INDXG2L = NB*((INDXGLOB-1)/(NB*NPROCS))+MOD(INDXGLOB-1,NB)+1
      lind_block_cyclic = INDXG2L
!
      RETURN
!
!     End of INDXG2L
!
      END FUNCTION lind_block_cyclic


!=----------------------------------------------------------------------------=!


      INTEGER FUNCTION gind_cyclic( lind, n, np, me )

!  This function computes the global index of a distributed array entry
!  pointed to by the local index lind of the process indicated by me.
!  lind      local index of the distributed matrix entry.
!  N         is the size of the global array.
!  me        The coordinate of the process whose local array row or
!            column is to be determined.
!  np        The total number processes over which the distributed
!            matrix is distributed.
!

            INTEGER, INTENT(IN) :: lind, n, me, np
            INTEGER r, q

            gind_cyclic = (lind-1) * np + me + 1

            RETURN
      END FUNCTION gind_cyclic


!=----------------------------------------------------------------------------=!


      INTEGER FUNCTION gind_block( lind, n, np, me )

!  This function computes the global index of a distributed array entry
!  pointed to by the local index lind of the process indicated by me.
!  lind      local index of the distributed matrix entry.
!  N         is the size of the global array.
!  me        The coordinate of the process whose local array row or
!            column is to be determined.
!  np        The total number processes over which the distributed
!            matrix is distributed.


            INTEGER, INTENT(IN) :: lind, n, me, np
            INTEGER r, q

              q = INT(n/np)
              r = MOD(n,np)
              IF( me < r ) THEN
                gind_block = (Q+1)*me + lind
              ELSE
                gind_block = Q*me + R + lind
              END IF

         RETURN
      END FUNCTION gind_block

!=----------------------------------------------------------------------------=!

      INTEGER FUNCTION gind_block_cyclic( lind, n, nb, np, me )

!  This function computes the global index of a distributed array entry
!  pointed to by the local index lind of the process indicated by me.
!  lind      local index of the distributed matrix entry.
!  N         is the size of the global array.
!  NB        size of the blocks the distributed matrix is split into.
!  me        The coordinate of the process whose local array row or
!            column is to be determined.
!  np        The total number processes over which the distributed
!            matrix is distributed.


            INTEGER, INTENT(IN) :: lind, n, nb, me, np
            INTEGER r, q, isrc

            isrc = 0
            gind_block_cyclic = np*NB*((lind-1)/NB) + &
                MOD(lind-1,NB) + MOD(np+me-isrc, np)*NB + 1

         RETURN
      END FUNCTION gind_block_cyclic

