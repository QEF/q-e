!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
integer function good_fft_dimension (n)
  !
  ! Determines the optimal maximum dimensions of fft arrays
  ! Useful on some machines to avoid memory conflicts
  !
#include "machine.h"
  use parameters
  implicit none
  integer :: n, nx
  ! this is the default: max dimension = fft dimension
  nx = n
#if defined(AIX) || defined(DXML)
  if ( n==8 .or. n==16 .or. n==32 .or. n==64 .or. n==128 .or. n==256) &
       nx = n + 1
#endif
#if defined(CRAYY) || defined(__SX4)
  if (mod (nr1, 2) ==0) nx = n + 1
#endif
  good_fft_dimension = nx
  return
end function good_fft_dimension

