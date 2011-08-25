!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
module betax
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE 
  SAVE
  !
  INTEGER              :: mmx = 5000
  REAL(DP)             :: refg
  REAL(DP),ALLOCATABLE :: betagx(:,:,:), dbetagx(:,:,:), &
                          qradx(:,:,:,:), dqradx(:,:,:,:)
  !
  CONTAINS
  !
  SUBROUTINE deallocate_betax()
    !
    IF ( ALLOCATED( betagx ) )  DEALLOCATE( betagx )
    IF ( ALLOCATED( dbetagx ) ) DEALLOCATE( dbetagx )
    IF ( ALLOCATED( qradx ) )   DEALLOCATE( qradx )
    IF ( ALLOCATED( dqradx ) )  DEALLOCATE( dqradx )
    !
  END SUBROUTINE deallocate_betax
  !
END MODULE betax
!
