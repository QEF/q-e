!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine hexsym (at, is, isname, nrot)  !-----------------------------------------------------------------------
  !
  ! Provides symmetry operations for Hexagonal and Trigonal lattices.
  ! The c axis is assumed to be along the z axis
  !
  use parameters
  implicit none  
  !
  !     first the input variables
  !
  real(kind=DP) :: at (3, 3)  
  ! input: the direct lattice vectors
  integer :: is (3, 3, 48), nrot  
  ! output: the symmetry matrices
  ! output: the number of symmetry matrice
  character :: isname (48) * 45  
  ! output: full name of the rotational part of each selected sym.op.
  !
  !    here the local parameters
  !
  real(kind=DP) :: sin3, cos3, msin3, mcos3  
  ! sin(pi/3), cos(pi/3), -sin(pi/3), -sin(pi/3)
  parameter (sin3 = 0.866025403784438597d0, cos3 = 0.5d0, msin3 = &
       - 0.866025403784438597d0, mcos3 = - 0.5d0)
  !
  !   and the local variables
  !
  real(kind=DP) :: s (3, 3, 12), overlap (3, 3), rat (3), rot (3, 3), &
       value
  ! the s matrices in real variables
  ! overlap matrix between direct lattice
  ! the rotated of a direct vector (cartesian)
  ! the rotated of a direct vector (crystal axis)
  integer :: jpol, kpol, mpol, irot
  ! counters over polarizations and rotations
  character :: sname (24) * 45  
  ! full name of the rotation part of each symmetry operation

  data s / 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
       - 1.d0, 0.d0, 0.d0, 0.d0, - 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, - 1.d0, &
       0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, - 1.d0, 1.d0, 0.d0, &
       0.d0, 0.d0, - 1.d0, 0.d0, 0.d0, 0.d0, - 1.d0, cos3, sin3, 0.d0, &
       msin3, cos3, 0.d0, 0.d0, 0.d0, 1.d0, cos3, msin3, 0.d0, sin3, &
       cos3, 0.d0, 0.d0, 0.d0, 1.d0, mcos3, sin3, 0.d0, msin3, mcos3, &
       0.d0, 0.d0, 0.d0, 1.d0, mcos3, msin3, 0.d0, sin3, mcos3, 0.d0, &
       0.d0, 0.d0, 1.d0, cos3, msin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, &
       0.d0, - 1.d0, cos3, sin3, 0.d0, sin3, mcos3, 0.d0, 0.d0, 0.d0, &
       - 1.d0, mcos3, msin3, 0.d0, msin3, cos3, 0.d0, 0.d0, 0.d0, &
       - 1.d0, mcos3, sin3, 0.d0, sin3, cos3, 0.d0, 0.d0, 0.d0, - 1.d0 /

  data sname / 'identity                                     ',&
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
  !
  !   first compute the overlap matrix between direct lattice vectors
  !
  do jpol = 1, 3  
     do kpol = 1, 3  
        overlap (kpol, jpol) = at (1, kpol) * at (1, jpol) + at (2, kpol) &
             * at (2, jpol) + at (3, kpol) * at (3, jpol)
     enddo
  enddo
  !
  !    then its inverse
  !
  call invmat (overlap, overlap, 3)  
  nrot = 1  
  do irot = 1, 12  
     !
     !   for each possible symmetry
     !
     do jpol = 1, 3  
        do mpol = 1, 3  
           !
           !   compute, in cartesian coordinates the rotated vector
           !
           rat (mpol) = s (mpol, 1, irot) * at (1, jpol) + s (mpol, 2, irot) &
                * at (2, jpol) + s (mpol, 3, irot) * at (3, jpol)
        enddo
        do kpol = 1, 3  
           !
           !   the rotated vector is projected on the direct lattice
           !
           rot (kpol, jpol) = at (1, kpol) * rat (1) + at (2, kpol) * rat (2) &
                + at (3, kpol) * rat (3)
        enddo
     enddo
     !
     !  and the inverse of the overlap matrix is applied
     !
     do jpol = 1, 3  
        do kpol = 1, 3  
           value = overlap (jpol, 1) * rot (1, kpol) + overlap (jpol, 2) &
                * rot (2, kpol) + overlap (jpol, 3) * rot (3, kpol)
           if (abs (float (nint (value) ) - value) .gt.1.0d-8) then  
              !
              ! if a noninteger is obtained, this implies that this operation
              ! is not a symmetry operation for the given lattice
              !
              goto 10  
           endif
           is (kpol, jpol, nrot) = nint (value)  
           isname (nrot) = sname (irot)  
        enddo
     enddo
     nrot = nrot + 1  
10   continue  
  enddo
  nrot = nrot - 1  
  !
  !   set the inversion symmetry ( Bravais lattices have always inversion)
  !
  do irot = 1, nrot  
     do kpol = 1, 3  
        do jpol = 1, 3  
           is (kpol, jpol, irot + nrot) = - is (kpol, jpol, irot)  
           isname (irot + nrot) = sname (irot + 12)  
        enddo
     enddo

  enddo

  nrot = 2 * nrot  
  return  
end subroutine hexsym
