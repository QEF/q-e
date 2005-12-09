!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine symrho (rho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
  !-----------------------------------------------------------------------
  !
  !     symmetrize the charge density.
  !
#include "f_defs.h"
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s (3, 3, 48), ftau (3, 48)
  !
  ! input:  dimensions of the FFT mesh
  ! input: the number of symmetries
  ! input: the symmetry matrices
  ! input: the fractionary translations
  !
  real(DP) :: rho (nrx1, nrx2, nrx3)
  ! inp/out: the charge density
  integer , allocatable :: symflag (:,:,:)
  integer :: ri (48), rj (48), rk (48), i, j, k, isym
  real(DP) :: sum

  if (nsym.eq.1) return

  allocate (symflag(nrx1, nrx2, nrx3))    
  do k = 1, nr3
     do j = 1, nr2
        do i = 1, nr1
           symflag (i, j, k) = 0
        enddo
     enddo
  enddo
  do k = 1, nr3
     do j = 1, nr2
        do i = 1, nr1
           if (symflag (i, j, k) .eq.0) then
              sum = 0.d0
              do isym = 1, nsym
                 call ruotaijk (s (1, 1, isym), ftau (1, isym), i, j, k, nr1, &
                      nr2, nr3, ri (isym), rj (isym), rk (isym) )
                 sum = sum + rho (ri (isym), rj (isym), rk (isym) )
              enddo
              sum = sum / nsym
              !
              !     sum contains the symmetrized charge density at point r.
              !     now fill the star of r with this sum.
              !
              do isym = 1, nsym
                 rho (ri (isym), rj (isym), rk (isym) ) = sum
                 symflag (ri (isym), rj (isym), rk (isym) ) = 1
              enddo
           endif
        enddo
     enddo

  enddo

  deallocate(symflag)
  return
end subroutine symrho

