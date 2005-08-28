!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine symrho_mag (rho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, &
     ftau, bg,at)
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
  integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s (3, 3, 48), &
       ftau (3, 48)
  !
  ! input:  dimensions of the FFT mesh
  ! input: the number of symmetries
  ! input: the symmetry matrices
  ! input: the fractionary translations
  !
  REAL(DP) :: bg(3,3), at(3,3)
  real(DP) :: rho (nrx1, nrx2, nrx3, 3)
  ! inp/out: the charge density
  integer , allocatable :: symflag (:,:,:)
  integer :: ri (48), rj (48), rk (48), kpol, i, j, k, isym
  real(DP) :: sumx, sumy, sumz, mag(3), magrot(3)
  ! auxiliary variables

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
              sumx = 0.d0
              sumy = 0.d0
              sumz = 0.d0
              do isym = 1, nsym
                 call ruotaijk (s (1, 1, isym), ftau (1, isym), i, j, k, nr1, &
                      nr2, nr3, ri (isym), rj (isym), rk (isym) )
! put the magnetic moment in crystal coordinates
                 do kpol = 1, 3
                    mag(kpol)=bg(1,kpol)* &
                    rho(ri(isym),rj(isym),rk(isym),1) + &
                    bg(2,kpol)*rho(ri(isym),rj(isym),rk(isym),2) + &
                    bg(3,kpol)*rho(ri(isym),rj(isym),rk(isym),3)
                 enddo
! rotate the magnetic moment
                do kpol = 1, 3
                    magrot(kpol) = s(1,kpol,isym)*mag(1) + &
                    s(2,kpol,isym)*mag(2) + &
                    s(3,kpol,isym)*mag(3)
                enddo
                sumx = sumx + magrot(1)
                sumy = sumy + magrot(2)
                sumz = sumz + magrot(3)
             enddo
             sumx = sumx/nsym
             sumy = sumy/nsym
             sumz = sumz/nsym
!
!     sum contains the symmetrised magnetisation at point r.
!     now fill the star of r with this (rotated) sum.
!
                  do isym = 1,nsym
                     mag(1) = sumx
                     mag(2) = sumy
                     mag(3) = sumz
! rotate the magnetic moment
                     do kpol = 1, 3
                        magrot(kpol) = s(1,kpol,isym)*mag(1) + &
                        s(2,kpol,isym)*mag(2) + &
                        s(3,kpol,isym)*mag(3)
                     enddo
! go back to carthesian coordinates
                     do kpol = 1, 3
                        mag(kpol)=at(kpol,1)*magrot(1) + &
                        at(kpol,2)*magrot(2) + &
                        at(kpol,3)*magrot(3)
                     enddo
                     rho(ri(isym),rj(isym),rk(isym),1) = mag(1)
                     rho(ri(isym),rj(isym),rk(isym),2) = mag(2)
                     rho(ri(isym),rj(isym),rk(isym),3) = mag(3)
                     symflag(ri(isym),rj(isym),rk(isym)) = 1
                  enddo
               endif
            enddo
         enddo
      enddo

      DEALLOCATE (symflag)

      RETURN
      end subroutine symrho_mag

