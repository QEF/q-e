!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!----------------------------------------------------------------------
subroutine set_pencils (nks, xk, ngm, gg, nl, gcut, nrx1, nr1, &
     nr2, nr3)
  !----------------------------------------------------------------------
  !
  !  Calculates active pencils for use in "short" fft (for wavefunctions)
  !  Warning: common /pencils/ is used to transmit variables to cfts_3
  !
#include "machine.h"
  use parameters
  use pencils
  implicit none
  integer :: nks, ngm, nl (ngm), nr1, nr2, nr3, nrx1
  real(kind=DP) :: xk (3, nks), gg (ngm), gcut
  integer :: kpoint, ng, n2, n3, n321, n21

  real(kind=DP) :: qmax

  allocate (do_fft_x ( nr2 , nr3))    
  allocate (do_fft_y ( nr3))    
  qmax = 0.0
  do kpoint = 1, nks
     qmax = max (qmax, xk (1, kpoint) **2 + xk (2, kpoint) **2 + xk (3, &
          kpoint) **2)
  enddo

  do_fft_y = 0
  do_fft_x = 0
  

  do ng = 1, ngm
     if (sqrt (gg (ng) ) .gt.sqrt (gcut) + sqrt (qmax) ) return
     n321 = nl (ng) - 1
     n3 = n321 / (nrx1 * nr2) + 1
     if (n3.le.0.or.n3.gt.nr3) call errore ('set_pencils', 'something wrong', 3)
     n21 = mod (n321, nrx1 * nr2)
     n2 = n21 / nrx1 + 1

     if (n2.le.0.or.n2.gt.nr2) call errore ('set_pencils', 'something wrong', 2)
     do_fft_x (n2, n3) = 1

     do_fft_y (n3) = 1

  enddo
  return

end subroutine set_pencils

