!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
MODULE supercell
  !---------------------------------------------------------------------------
  !
  USE parameters,  ONLY : DP
  USE brilz,       ONLY : alat, at
  !
  IMPLICIT NONE
  !
  CONTAINS
     !
     PURE FUNCTION pbc( vect )
       !
       IMPLICIT NONE    
       !
       REAL (KIND=DP), DIMENSION(:), INTENT(IN) :: vect
       REAL (KIND=DP), DIMENSION(SIZE( vect ))  :: pbc
       REAL (KIND=DP), DIMENSION(3)             :: crystal
       INTEGER                                  :: i, j, index  
       !
       !
       pbc = 0.D0
       !
       DO i = 1, SIZE( vect ), 3
          !
          index = ( i + 3 - 1 )
          !
          DO j = 1, 3
             !
             ! ... conversion to crystal coordinates :
             ! ... for each atom the projection of the position vector on the 
             ! ... new basis vectors is taken. Moreover they are normalized
             !
             crystal(j) = DOT_PRODUCT( vect(i:index), at(:,j) ) / &
                          ( DOT_PRODUCT( at(:,j), at(:,j) ) * alat )
             !
             ! ... cubic pbc are applied to the position vector
             !
             crystal(j) = crystal(j) - ANINT( crystal(j) )
             !
             ! ... back conversion to cartesian coordinates
             !
             pbc(i:index) = pbc(i:index) + crystal(j) * at(:,j) * alat
             !
          END DO
          !
       END DO  
       !
       RETURN
       !
     END FUNCTION pbc
     !
END MODULE supercell 
