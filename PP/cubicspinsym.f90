!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE cubicspinsym(d_spin)
  !-----------------------------------------------------------------------
  !
  ! Provides symmetry operations in the spin space for all cubic and
  ! lower-symmetry (excepted Hexagonal and Trigonal) bravais lattices
  !
  USE kinds
  USE symm_base, ONLY : nsym, sname
  !
  IMPLICIT NONE
  !
  !     first the input/output variables
  !
  COMPLEX(DP) :: d_spin(2,2,48)
  ! output: the symmetries of the system in the j=1/2 subspace
  !
  !    here the local parameters
  !
  REAL(DP), PARAMETER :: cost1=0.707106781186547d0, cost2=0.5d0, &
                         mcost1=-0.707106781186547d0, mcost2=-0.5d0
  !
  !   and the local variables
  !
  COMPLEX(DP) :: id_spin(2,2,48)
  !
  ! symmetry matrices in the j=1/2 subspace
  !
  INTEGER :: i,j
  !
  ! counter over the rotations
  ! counter over the rotations
  !
  CHARACTER :: isname (48) * 45
  !
  ! full name of the rotational part of each symmetry operation
  !
  data isname/&
       &        'identity                                    ',&
       &        '180 deg rotation - cart. axis [0,0,1]       ',&
       &        '180 deg rotation - cart. axis [0,1,0]       ',&
       &        '180 deg rotation - cart. axis [1,0,0]       ',&
       &        '180 deg rotation - cart. axis [1,1,0]       ',&
       &        '180 deg rotation - cart. axis [1,-1,0]      ',&
       &        ' 90 deg rotation - cart. axis [0,0,-1]      ',&
       &        ' 90 deg rotation - cart. axis [0,0,1]       ',&
       &        '180 deg rotation - cart. axis [1,0,1]       ',&
       &        '180 deg rotation - cart. axis [-1,0,1]      ',&
       &        ' 90 deg rotation - cart. axis [0,1,0]       ',&
       &        ' 90 deg rotation - cart. axis [0,-1,0]      ',&
       &        '180 deg rotation - cart. axis [0,1,1]       ',&
       &        '180 deg rotation - cart. axis [0,1,-1]      ',&
       &        ' 90 deg rotation - cart. axis [-1,0,0]      ',&
       &        ' 90 deg rotation - cart. axis [1,0,0]       ',&
       &        '120 deg rotation - cart. axis [-1,-1,-1]    ',&
       &        '120 deg rotation - cart. axis [-1,1,1]      ',&
       &        '120 deg rotation - cart. axis [1,1,-1]      ',&
       &        '120 deg rotation - cart. axis [1,-1,1]      ',&
       &        '120 deg rotation - cart. axis [1,1,1]       ',&
       &        '120 deg rotation - cart. axis [-1,1,-1]     ',&
       &        '120 deg rotation - cart. axis [1,-1,-1]     ',&
       &        '120 deg rotation - cart. axis [-1,-1,1]     ',&
       &        'inversion                                    ',&
       &        'inv. 180 deg rotation - cart. axis [0,0,1]  ',&
       &        'inv. 180 deg rotation - cart. axis [0,1,0]  ',&
       &        'inv. 180 deg rotation - cart. axis [1,0,0]  ',&
       &        'inv. 180 deg rotation - cart. axis [1,1,0]  ',&
       &        'inv. 180 deg rotation - cart. axis [1,-1,0] ',&
       &        'inv.  90 deg rotation - cart. axis [0,0,-1] ',&
       &        'inv.  90 deg rotation - cart. axis [0,0,1]  ',&
       &        'inv. 180 deg rotation - cart. axis [1,0,1]  ',&
       &        'inv. 180 deg rotation - cart. axis [-1,0,1] ',&
       &        'inv.  90 deg rotation - cart. axis [0,1,0]  ',&
       &        'inv.  90 deg rotation - cart. axis [0,-1,0] ',&
       &        'inv. 180 deg rotation - cart. axis [0,1,1]  ',&
       &        'inv. 180 deg rotation - cart. axis [0,1,-1] ',&
       &        'inv.  90 deg rotation - cart. axis [-1,0,0] ',&
  &        'inv.  90 deg rotation - cart. axis [1,0,0]  ',&
  &        'inv. 120 deg rotation - cart. axis [-1,-1,-1]',&
  &        'inv. 120 deg rotation - cart. axis [-1,1,1] ',&
  &        'inv. 120 deg rotation - cart. axis [1,1,-1]' ,&
  &        'inv. 120 deg rotation - cart. axis [1,-1,1] ',&
  &        'inv. 120 deg rotation - cart. axis [1,1,1]  ',&
  &        'inv. 120 deg rotation - cart. axis [-1,1,-1] ',&
  &        'inv. 120 deg rotation - cart. axis [1,-1,-1]',&
  &        'inv. 120 deg rotation - cart. axis [-1,-1,1] ' /
  data id_spin / (-1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (-1.d0,0.d0), &
                (0.d0,-1.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,1.d0), &
                (0.d0,0.d0), (1.d0,0.d0), (-1.d0,0.d0), (0.d0,0.d0), &
                (0.d0,0.d0), (0.d0,-1.d0), (0.d0,-1.d0), (0.d0,0.d0), &
                (0.d0,0.d0), (mcost1,mcost1), (cost1,mcost1),(0.d0,0.d0), &
                (0.d0,0.d0), (cost1,mcost1), (mcost1,mcost1), (0.d0,0.d0), &
                (cost1,mcost1), (0.d0,0.d0), (0.d0,0.d0), (cost1,cost1), &
                (mcost1,mcost1),(0.d0,0.d0), (0.d0,0.d0), (mcost1,cost1), &
                (0.d0,mcost1), (0.d0,mcost1), (0.d0,mcost1), (0.d0,cost1), &
                (0.d0,mcost1), (0.d0,cost1), (0.d0,cost1), (0.d0,cost1), &
                (mcost1,0.d0), (mcost1,0.d0), (cost1,0.d0), (mcost1,0.d0), &
                (cost1,0.d0), (mcost1,0.d0), (cost1,0.d0), (cost1,0.d0), &
                (0.d0,mcost1), (mcost1,0.d0), (cost1,0.d0), (0.d0,cost1), &
                (0.d0,cost1), (mcost1,0.d0), (cost1,0.d0), (0.d0,mcost1), &
                (mcost1,0.d0), (0.d0,cost1), (0.d0,cost1), (mcost1,0.d0), &
                (mcost1,0.d0), (0.d0,mcost1), (0.d0,mcost1), (mcost1,0.d0), &
                (cost2,mcost2),(mcost2,mcost2),(cost2,mcost2),(cost2,cost2), &
                (mcost2,mcost2),(mcost2,cost2),(cost2,cost2),(mcost2,cost2), &
                (mcost2,cost2),(mcost2,mcost2),(cost2,mcost2),(mcost2,mcost2),&
                (mcost2,mcost2),(cost2,mcost2),(mcost2,mcost2),(mcost2,cost2),&
                (mcost2,mcost2),(mcost2,mcost2),(cost2,mcost2),(mcost2,cost2),&
                (mcost2,cost2),(mcost2,cost2),(cost2,cost2),(mcost2,mcost2),&
                (cost2,mcost2),(mcost2,cost2),(cost2,cost2),(cost2,cost2), &
                (mcost2,mcost2),(cost2,cost2),(mcost2,cost2),(mcost2,cost2), &
                (-1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (-1.d0,0.d0), &
                (0.d0,-1.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,1.d0), &
                (0.d0,0.d0), (1.d0,0.d0), (-1.d0,0.d0), (0.d0,0.d0), &
                (0.d0,0.d0), (0.d0,-1.d0), (0.d0,-1.d0), (0.d0,0.d0), &
                (0.d0,0.d0), (mcost1,mcost1), (cost1,mcost1),(0.d0,0.d0), &
                (0.d0,0.d0), (cost1,mcost1), (mcost1,mcost1), (0.d0,0.d0), &
                (cost1,mcost1), (0.d0,0.d0), (0.d0,0.d0), (cost1,cost1), &
                (mcost1,mcost1),(0.d0,0.d0), (0.d0,0.d0), (mcost1,cost1), &
                (0.d0,mcost1), (0.d0,mcost1), (0.d0,mcost1), (0.d0,cost1), &
                (0.d0,mcost1), (0.d0,cost1), (0.d0,cost1), (0.d0,cost1), &
                (mcost1,0.d0), (mcost1,0.d0), (cost1,0.d0), (mcost1,0.d0), &
                (cost1,0.d0), (mcost1,0.d0), (cost1,0.d0), (cost1,0.d0), &
                (0.d0,mcost1), (mcost1,0.d0), (cost1,0.d0), (0.d0,cost1), &
                (0.d0,cost1), (mcost1,0.d0), (cost1,0.d0), (0.d0,mcost1), &
                (mcost1,0.d0), (0.d0,cost1), (0.d0,cost1), (mcost1,0.d0), &
                (mcost1,0.d0), (0.d0,mcost1), (0.d0,mcost1), (mcost1,0.d0), &
                (cost2,mcost2),(mcost2,mcost2),(cost2,mcost2),(cost2,cost2), &
                (mcost2,mcost2),(mcost2,cost2),(cost2,cost2),(mcost2,cost2), &
                (mcost2,cost2),(mcost2,mcost2),(cost2,mcost2),(mcost2,mcost2),&
                (mcost2,mcost2),(cost2,mcost2),(mcost2,mcost2),(mcost2,cost2),&
                (mcost2,mcost2),(mcost2,mcost2),(cost2,mcost2),(mcost2,cost2),&
                (mcost2,cost2),(mcost2,cost2),(cost2,cost2),(mcost2,mcost2), &
                (cost2,mcost2),(mcost2,cost2),(cost2,cost2),(cost2,cost2), &
                (mcost2,mcost2),(cost2,cost2), (mcost2,cost2), (mcost2,cost2) /

  d_spin = (0.d0,0.d0)
  DO i=1,nsym
     DO j=1,48
        IF (sname(i)==isname(j)) THEN
           d_spin(1,1,i)=id_spin(1,1,j)
           d_spin(2,1,i)=id_spin(2,1,j)
           d_spin(1,2,i)=id_spin(1,2,j)
           d_spin(2,2,i)=id_spin(2,2,j)
        ENDIF
     ENDDO
  ENDDO
  RETURN
  !
  END SUBROUTINE cubicspinsym
