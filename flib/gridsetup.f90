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

SUBROUTINE GRID2D_DIMS( grid_shape, nproc, nprow, npcol )
   !
   ! This subroutine factorizes the number of processors (NPROC)
   ! into NPROW and NPCOL according to the shape
   !
   !    Written by Carlo Cavazzoni
   !
   IMPLICIT NONE
   CHARACTER, INTENT(IN) :: grid_shape
   INTEGER, INTENT(IN)  :: nproc
   INTEGER, INTENT(OUT) :: nprow, npcol
   INTEGER :: sqrtnp, i
   !
   sqrtnp = INT( SQRT( REAL( nproc ) + 0.1 ) )
   !
   IF( grid_shape == 'S' ) THEN
      ! Square grid
      nprow = sqrtnp
      npcol = sqrtnp
   ELSE
      ! Rectangular grid
      DO i = 1, sqrtnp + 1
         IF( MOD( nproc, i ) == 0 ) nprow = i
      end do
      npcol = nproc / nprow
   END IF
   RETURN
END SUBROUTINE

SUBROUTINE GRID2D_COORDS( order, rank, nprow, npcol, row, col )
   !
   !  this subroutine compute the cartesian coordinetes "row" and "col"
   !  of the processor whose MPI task id is "rank". 
   !  Note that if the rank is larger that the grid size
   !  all processors whose MPI task id is greather or equal 
   !  than nprow * npcol are placed on the diagonal extension of the grid itself
   !
   IMPLICIT NONE
   CHARACTER, INTENT(IN) :: order
   INTEGER, INTENT(IN)  ::  rank          ! process index starting from 0
   INTEGER, INTENT(IN)  ::  nprow, npcol  ! dimensions of the processor grid
   INTEGER, INTENT(OUT) ::  row, col
   IF( rank >= 0 .AND. rank < nprow * npcol ) THEN
      IF( order == 'C' .OR. order == 'c' ) THEN
         !  grid in COLUMN MAJOR ORDER
         row = MOD( rank, nprow )
         col = rank / nprow
      ELSE
         !  grid in ROW MAJOR ORDER
         row = rank / npcol
         col = MOD( rank, npcol )
      END IF
   ELSE
      row = rank
      col = rank
   END IF
   RETURN
END SUBROUTINE

SUBROUTINE GRID2D_RANK( order, nprow, npcol, row, col, rank )
   !
   !  this subroutine compute the processor MPI task id "rank" of the processor  
   !  whose cartesian coordinate are "row" and "col".
   !  Note that the subroutine assume cyclic indexing ( row = nprow = 0 )
   !
   IMPLICIT NONE
   CHARACTER, INTENT(IN) :: order
   INTEGER, INTENT(OUT) ::  rank         ! process index starting from 0
   INTEGER, INTENT(IN)  ::  nprow, npcol  ! dimensions of the processor grid
   INTEGER, INTENT(IN)  ::  row, col
   
   IF( order == 'C' .OR. order == 'c' ) THEN
     !  grid in COLUMN MAJOR ORDER
     rank = MOD( row + nprow, nprow ) + MOD( col + npcol, npcol ) * nprow
   ELSE
     !  grid in ROW MAJOR ORDER
     rank = MOD( col + npcol, npcol ) + MOD( row + nprow, nprow ) * npcol
   END IF
   !
   RETURN
END SUBROUTINE
