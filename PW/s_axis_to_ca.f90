!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine s_axis_to_cart (s, sr, at, bg)  
  !----------------------------------------------------------------------
  !
  !     This routine transform a symmetry matrix expressed in the
  !     basis of the crystal axis in the cartesian basis.
  !
  !     last revised 3 may 1995 by A. Dal Corso
  !
  !
  use parameters
  implicit none  
  !
  !     first the input parameters
  !
  integer :: s (3, 3)  
  ! input: matrix in crystal axis
  real(kind=DP) :: sr (3, 3), at (3, 3), bg (3, 3)  
  ! output: matrix in cartesian axis
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors
  !
  !     here the local variable
  !

  integer :: apol, bpol, kpol, lpol  
  !\
  !
  !     counters on polarizations
  !/
  do apol = 1, 3  
     do bpol = 1, 3  
        sr (apol, bpol) = 0.d0  
        do kpol = 1, 3  
           do lpol = 1, 3  
              sr (apol, bpol) = sr (apol, bpol) + at (apol, kpol) * float (s ( &
                   lpol, kpol) ) * bg (bpol, lpol)
           enddo
        enddo
     enddo
  enddo
  return  
end subroutine s_axis_to_cart
