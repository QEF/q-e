!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hexspinsym(d_spin)
  !-----------------------------------------------------------------------
  !
  ! Provides symmetry operations in the spin space for Hexagonal and
  ! Trigonal lattices.
  ! The c axis is assumed to be along the z axis
  !
  !
  USE kinds
  USE symm_base, ONLY : nsym, sname
  !
  IMPLICIT NONE
  !
  !     first the input/output variables
  !
  COMPLEX(DP) :: d_spin(2,2,48)
  !
  ! output: the symmetries of the system in the j=1/2 subspace
  !
  !    here the local parameters
  !
  ! sin3 = sin(pi/3), cos3 = cos(pi/3), msin3 = -sin(pi/3), mcos3 = -sin(pi/3)
  !
  REAL(DP), PARAMETER :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
                         msin3 =-0.866025403784438597d0, mcos3 = -0.5d0
  !
  !   and the local variables
  !
  COMPLEX(DP) :: id_spin(2,2,24)
  !
  ! symmetry matrices in the j=1/2 subspace
  !
  INTEGER :: i,j
  !
  ! counter over the rotations
  ! counter over the rotations
  !
  CHARACTER :: isname (24) * 45
  !
  ! full name of the rotational part of each symmetry operation
  !
  data isname / 'identity                                     ',&
               '180 deg rotation - cryst. axis [0,0,1]       ',&
               '180 deg rotation - cryst. axis [1,2,0]       ',&
               '180 deg rotation - cryst. axis [1,0,0]       ',&
               ' 60 deg rotation - cryst. axis [0,0,1]       ',&
               ' 60 deg rotation - cryst. axis [0,0,-1]      ',&
               '120 deg rotation - cryst. axis [0,0,1]       ',&
               '120 deg rotation - cryst. axis [0,0,-1]      ',&
               '180 deg rotation - cryst. axis [1,-1,0]      ',&
               '180 deg rotation - cryst. axis [2,1,0]       ',&
               '180 deg rotation - cryst. axis [0,1,0]       ',&
               '180 deg rotation - cryst. axis [1,1,0]       ',&
               'inversion                                    ',&
               'inv. 180 deg rotation - cryst. axis [0,0,1]  ',&
               'inv. 180 deg rotation - cryst. axis [1,2,0]  ',&
               'inv. 180 deg rotation - cryst. axis [1,0,0]  ',&
               'inv.  60 deg rotation - cryst. axis [0,0,1]  ',&
               'inv.  60 deg rotation - cryst. axis [0,0,-1] ',&
               'inv. 120 deg rotation - cryst. axis [0,0,1]  ',&
               'inv. 120 deg rotation - cryst. axis [0,0,-1] ',&
               'inv. 180 deg rotation - cryst. axis [1,-1,0] ',&
               'inv. 180 deg rotation - cryst. axis [2,1,0]  ',&
               'inv. 180 deg rotation - cryst. axis [0,1,0]  ',&
               'inv. 180 deg rotation - cryst. axis [1,1,0]  ' /
  data id_spin / (-1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (-1.d0,0.d0), &
                (0.d0,-1.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,1.d0), &
                (0.d0,0.d0), (1.d0,0.d0), (-1.d0,0.d0), (0.d0,0.d0), &
                (0.d0,0.d0), (0.d0,-1.d0), (0.d0,-1.d0), (0.d0,0.d0), &
                (msin3,mcos3), (0.d0,0.d0), (0.d0,0.d0), (msin3,cos3), &
                (sin3,mcos3), (0.d0,0.d0), (0.d0,0.d0), (sin3,cos3), &
                (mcos3,msin3), (0.d0,0.d0), (0.d0,0.d0), (mcos3,sin3), &
                (cos3,msin3), (0.d0,0.d0), (0.d0,0.d0), (cos3,sin3), &
                (0.d0,0.d0), (cos3,msin3), (mcos3,msin3), (0.d0,0.d0), &
                (0.d0,0.d0), (mcos3,msin3), (cos3,msin3), (0.d0,0.d0), &
                (0.d0,0.d0), (sin3,mcos3), (msin3,mcos3), (0.d0,0.d0), &
                (0.d0,0.d0), (msin3,mcos3), (sin3,mcos3), (0.d0,0.d0), &
                (-1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (-1.d0,0.d0), &
                (0.d0,-1.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,1.d0), &
                (0.d0,0.d0), (1.d0,0.d0), (-1.d0,0.d0), (0.d0,0.d0), &
                (0.d0,0.d0), (0.d0,-1.d0), (0.d0,-1.d0), (0.d0,0.d0), &
                (msin3,mcos3), (0.d0,0.d0), (0.d0,0.d0), (msin3,cos3), &
                (sin3,mcos3), (0.d0,0.d0), (0.d0,0.d0), (sin3,cos3), &
                (mcos3,msin3), (0.d0,0.d0), (0.d0,0.d0), (mcos3,sin3), &
                (cos3,msin3), (0.d0,0.d0), (0.d0,0.d0), (cos3,sin3), &
                (0.d0,0.d0), (cos3,msin3), (mcos3,msin3), (0.d0,0.d0), &
                (0.d0,0.d0), (mcos3,msin3), (cos3,msin3), (0.d0,0.d0), &
                (0.d0,0.d0), (sin3,mcos3), (msin3,mcos3), (0.d0,0.d0), &
                (0.d0,0.d0), (msin3,mcos3), (sin3,mcos3), (0.d0,0.d0) /

  d_spin = (0.d0,0.d0)
  DO i= 1, nsym
     DO j=1,24
        IF (sname(i)==isname(j)) THEN
           d_spin(:,:,i)=id_spin(:,:,j)
        ENDIF
     ENDDO
  ENDDO
  RETURN
  !
  END SUBROUTINE hexspinsym
