!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine trnvect (vect, at, bg, iflag)  
  !-----------------------------------------------------------------------
  !
  !  This routine transforms a vector (like forces which in the
  !  crystal axis is represented on the basis of the reciprocal lattice
  !  vectors) from crystal to cartesian axis (iflag.gt.0)
  !  and viceversa (iflag.le.0)
  !
  use parameters
  implicit none  
  integer :: iflag  
  ! input: gives the versus of the transformati

  real(kind=DP) :: vect (3), at (3, 3), bg (3, 3)  
  ! inp/out: the vector to transform
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors
  real(kind=DP) :: work (3)  
  ! a working array

  integer :: ipol, ialpha  
  ! counter on crystal coordinates
  ! counter on cartesian coordinates
  if (iflag.gt.0) then  
     !
     !     forward transformation, from crystal to cartesian axis
     !
     do ipol = 1, 3  
        work (ipol) = vect (ipol)  
     enddo
     do ialpha = 1, 3  
        vect (ialpha) = 0.d0  
        do ipol = 1, 3  
           vect (ialpha) = vect (ialpha) + work (ipol) * bg (ialpha, ipol)  
        enddo
     enddo
  else  
     !
     !    backward transformation, from cartesian to crystal axis
     !
     do ipol = 1, 3  
        work (ipol) = 0.d0  
        do ialpha = 1, 3  
           work (ipol) = work (ipol) + vect (ialpha) * at (ialpha, ipol)  
        enddo
     enddo
     do ipol = 1, 3  
        vect (ipol) = work (ipol)  
     enddo
  endif
  return  
end subroutine trnvect
