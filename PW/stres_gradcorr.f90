!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!
!--------------------------------------------------------------------
subroutine stres_gradcorr (rho, rho_core, nspin, nr1, nr2, nr3, &
     nrx1, nrx2, nrx3, nrxx, nl, ngm, g, alat, omega, sigmaxc)
  !--------------------------------------------------------------------
#include "f_defs.h"
  USE kinds
  use funct
  implicit none

  integer :: nspin, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, &
       nl (ngm)
  real(kind=DP) :: rho (nrxx, nspin), rho_core (nrxx), g (3, ngm), &
       alat, omega, sigmaxc (3, 3)
  integer :: k, l, m, ipol, is
  real(kind=DP) , allocatable :: grho (:,:,:)
  real(kind=DP), parameter :: epsr = 1.0d-6, epsg = 1.0d-10, e2 = 2.d0
  real(kind=DP) :: grh2, grho2 (2), sx, sc, v1x, v2x, v1c, v2c, fac, &
       v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw, zeta, rh, &
       sigma_gradcorr (3, 3)

  if (igcx.eq.0.and.igcc.eq.0) return

  sigma_gradcorr(:,:) = 0.d0

  allocate (grho( 3, nrxx, nspin))    
  fac = 1.d0 / dble (nspin)
  !
  !    calculate the gradient of rho+rhocore in real space
  !
  do is = 1, nspin
     call DAXPY (nrxx, fac, rho_core, 1, rho (1, is), 1)
     call gradient (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, rho (1, is), &
          ngm, g, nl, alat, grho (1, 1, is) )

  enddo
  if (nspin.eq.1) then
     !
     !    This is the LDA case
     !
     ! sigma_gradcor_{alpha,beta} ==
     !     omega^-1 \int (grad_alpha rho) ( D(rho*Exc)/D(grad_alpha rho) ) d3
     !
     do k = 1, nrxx
        grho2 (1) = grho (1, k, 1) **2 + grho (2, k, 1) **2 + grho (3, &
             k, 1) **2
        if (abs (rho (k, 1) ) .gt.epsr.and.grho2 (1) .gt.epsg) then
           call gcxc (rho (k, 1), grho2, sx, sc, v1x, v2x, v1c, v2c)
           do l = 1, 3
              do m = 1, l
                 sigma_gradcorr (l, m) = sigma_gradcorr (l, m) + grho (l, k, &
                      1) * grho (m, k, 1) * e2 * (v2x + v2c)
              enddo
           enddo
        endif

     enddo
  else
     !
     !    This is the LSDA case
     !
     do k = 1, nrxx
        grho2 (1) = grho (1, k, 1) **2 + grho (2, k, 1) **2 + grho (3, &
             k, 1) **2
        grho2 (2) = grho (1, k, 2) **2 + grho (2, k, 2) **2 + grho (3, &
             k, 2) **2

        if ( (abs (rho (k, 1) ) .gt.epsr.and.grho2 (1) .gt.epsg) .and. &
             (abs (rho (k, 2) ) .gt.epsr.and.grho2 (2) .gt.epsg) ) then
           call gcx_spin (rho (k, 1), rho (k, 2), grho2 (1), grho2 (2), &
                sx, v1xup, v1xdw, v2xup, v2xdw)
           rh = rho (k, 1) + rho (k, 2)
           if (rh.gt.epsr) then
              zeta = (rho (k, 1) - rho (k, 2) ) / rh
              grh2 = (grho (1, k, 1) + grho (1, k, 2) ) **2 + (grho (2, &
                   k, 1) + grho (2, k, 2) ) **2 + (grho (3, k, 1) + grho (3, &
                   k, 2) ) **2
              call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
           else
              sc = 0.d0
              v1cup = 0.d0
              v1cdw = 0.d0
              v2c = 0.d0

           endif
           do l = 1, 3
              do m = 1, l
                 !    exchange
                 sigma_gradcorr (l, m) = sigma_gradcorr (l, m) + grho (l, k, &
                      1) * grho (m, k, 1) * e2 * v2xup + grho (l, k, 2) * grho (m, &
                      k, 2) * e2 * v2xdw
                 !    correlation
                 sigma_gradcorr (l, m) = sigma_gradcorr (l, m) + (grho (l, k, &
                      1) + grho (l, k, 2) ) * (grho (m, k, 1) + grho (m, k, 2) ) &
                      * e2 * v2c
              enddo
           enddo
        endif
     enddo

  endif
  do l = 1, 3
     do m = 1, l - 1
        sigma_gradcorr (m, l) = sigma_gradcorr (l, m)
     enddo

  enddo
#ifdef __PARA
  call reduce (9, sigma_gradcorr)
#endif

  call DSCAL (9, 1.0d0 / (nr1 * nr2 * nr3), sigma_gradcorr, 1)

  call DAXPY (9, 1.d0, sigma_gradcorr, 1, sigmaxc, 1)
  do is = 1, nspin
     call DAXPY (nrxx, - fac, rho_core, 1, rho (1, is), 1)

  enddo
  deallocate(grho)
  return

end subroutine stres_gradcorr

