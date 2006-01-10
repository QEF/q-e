!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
!-----------------------------------------------------------------------
!

SUBROUTINE GRID2D_DIMS( nproc, nprow, npcol )
   !
   ! This subroutine factorizes the number of processors (NPROC)
   ! into NPROW and NPCOL,  that are the sizes of the 2D processors mesh.
   !
   !    Written by Carlo Cavazzoni
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: nproc
   INTEGER, INTENT(OUT) :: nprow, npcol
   integer sqrtnp,i
   sqrtnp = int( sqrt( dble(nproc) ) + 1 )
   do i=1,sqrtnp
      if(mod(nproc,i).eq.0) nprow = i
   end do
   npcol = nproc/nprow
   RETURN
END SUBROUTINE

SUBROUTINE GRID2D_COORDS( rank, nprow, npcol, row, col )
   IMPLICIT NONE
   INTEGER, INTENT(IN)  ::  rank          ! process index starting from 0
   INTEGER, INTENT(IN)  ::  nprow, npcol  ! dimensions of the processor grid
   INTEGER, INTENT(OUT) ::  row, col
   row = MOD( rank, nprow )
   col = rank / nprow
   RETURN
END SUBROUTINE

SUBROUTINE GRID2D_RANK( nprow, npcol, row, col, rank )
   IMPLICIT NONE
   INTEGER, INTENT(OUT) ::  rank         ! process index starting from 0
   INTEGER, INTENT(IN)  ::  nprow, npcol  ! dimensions of the processor grid
   INTEGER, INTENT(IN)  ::  row, col
   rank = row + col * nprow
   RETURN
END SUBROUTINE
