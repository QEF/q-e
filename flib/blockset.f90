!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!----------------------------------------------------------------------
!

      SUBROUTINE BLOCKSET( NB, NBUSER, N, NPROW, NPCOL)
!     ..
!     This subroutine try to choose an optimal block size
!     for the distributed matrix.
!
!     Written by Carlo Cavazzoni
!     ..
        IMPLICIT NONE
        INTEGER NB, N, NPROW, NPCOL, NBUSER
!     ..
        NB = MIN ( N/NPROW, N/NPCOL )
        IF(NBUSER.gt.0) then
          NB = MIN ( NB, NBUSER )
        ENDIF
        NB = MAX(NB,1)

      return
      end SUBROUTINE BLOCKSET

