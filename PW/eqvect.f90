!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
logical function eqvect (x, y, f)  
  !-----------------------------------------------------------------------
  !
  !   This function test if the difference between two tridimensional
  !   vectors is an integer. The presence of a fractionary translation
  !   may be required. (f)
  !
  !   Last revision  June 1997 (PG+SdG)
  !
  use parameters
  implicit none  
  real(kind=DP) :: x (3), y (3), f (3)  
  ! input: input vector
  ! input: second input vector
  ! input: fractionary translation
  real(kind=DP) :: accep  
  ! acceptance parameter
  parameter (accep = 1.0d-5)  
  !
  eqvect = abs (x (1) - y (1) - f (1) - nint (x (1) - y (1) - f (1) &
       ) ) .lt.accep.and.abs (x (2) - y (2) - f (2) - nint (x (2) &
       - y (2) - f (2) ) ) .lt.accep.and.abs (x (3) - y (3) - f (3) &
       - nint (x (3) - y (3) - f (3) ) ) .lt.accep
  return  
end function eqvect
