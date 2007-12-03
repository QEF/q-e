!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE becmod
  !
  USE control_flags, ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, ONLY :  DP
  !
  SAVE
  ! ... *bec* contain <beta|psi> - used in h_psi, s_psi, many other places
  REAL(DP), ALLOCATABLE :: &
       rbecp(:,:) !   <beta|psi> for real (at Gamma) wavefunctions 
  COMPLEX(DP), ALLOCATABLE ::  &
       becp (:,:), &            !  as above for complex wavefunctions
       becp_nc(:,:,:)           !  as above for spinors
  !
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE allocate_bec ( nkb, nbnd )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nkb, nbnd
    !
    IF ( gamma_only ) THEN 
       !
       ALLOCATE( rbecp( nkb, nbnd ) )
       !
    ELSE IF ( noncolin) THEN
       !
       ALLOCATE( becp_nc( nkb, npol, nbnd ) )
       !
    ELSE
       !
       ALLOCATE( becp( nkb, nbnd ) )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE allocate_bec
  !
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_bec ()
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    IF ( gamma_only ) THEN 
       !
       DEALLOCATE( rbecp )
       !
    ELSE IF ( noncolin) THEN
       !
       DEALLOCATE( becp_nc )
       !
    ELSE
       !
       DEALLOCATE( becp )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE deallocate_bec
  !
END MODULE becmod
