!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine gradcorr (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, &
     nrx3, nrxx, nl, ngm, g, alat, omega, nspin, etxc, vtxc, v)
  !     ===================
  !--------------------------------------------------------------------
#include "machine.h"
  USE kinds
  use funct
  implicit none
  !
  integer :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, nl (ngm), &
       nspin

  real(kind=DP) :: rho (nrxx, nspin), rho_core (nrxx), v (nrxx, nspin), &
       g (3, ngm), vtxc, etxc, alat, omega, zeta, rh, grh2
  integer :: k, ipol, is
  real(kind=DP), allocatable :: grho (:,:,:), h (:,:,:), dh (:)
  real(kind=DP) :: grho2 (2), sx, sc, v1x, v2x, v1c, v2c, v1xup, v1xdw, &
       v2xup, v2xdw, v1cup, v1cdw , etxcgc, vtxcgc, segno, arho, fac
  real(kind=DP), parameter :: e2 = 2.d0, epsr = 1.0d-6, epsg = 1.0d-10

  if (igcx == 0 .and. igcc == 0) return
  etxcgc = 0.d0
  vtxcgc = 0.d0

  allocate (h( 3, nrxx, nspin))    
  allocate (grho( 3, nrxx, nspin))    

  ! calculate the gradient of rho+rho_core in real space
  fac = 1.d0 / dble (nspin)
  do is = 1, nspin
     call DAXPY (nrxx, fac, rho_core, 1, rho (1, is), 1)
     call gradient (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, rho (1, is), &
          ngm, g, nl, alat, grho (1, 1, is) )
  enddo
  do k = 1, nrxx
     do is = 1, nspin
        grho2 (is) = grho(1, k, is)**2 + grho(2, k, is)**2 + grho(3, k, is)**2
     enddo
     if (nspin == 1) then
        !
        !    This is the spin-unpolarised case
        !
        arho = abs (rho (k, 1) )
        segno = sign (1.d0, rho (k, 1) )
        if (arho.gt.epsr.and.grho2 (1) .gt.epsg) then

           call gcxc (arho, grho2, sx, sc, v1x, v2x, v1c, v2c)
           !
           ! first term of the gradient correction : D(rho*Exc)/D(rho)

           v (k, 1) = v (k, 1) + e2 * (v1x + v1c)
           ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
           do ipol = 1, 3
              h (ipol, k, 1) = e2 * (v2x + v2c) * grho (ipol, k, 1)
           enddo
           vtxcgc = vtxcgc + e2 * (v1x + v1c) * (rho (k, 1) - rho_core(k) )
           etxcgc = etxcgc + e2 * (sx + sc) * segno
        else
           do ipol = 1, 3
              h (ipol, k, 1) = 0.d0
           enddo
        endif
     else
        !
        !    spin-polarised case
        !
        call gcx_spin (rho (k, 1), rho (k, 2), grho2 (1), grho2 (2), &
             sx, v1xup, v1xdw, v2xup, v2xdw)
        rh = rho (k, 1) + rho (k, 2)
        if (rh.gt.epsr) then
           zeta = (rho (k, 1) - rho (k, 2) ) / rh
           !
           ! ... ( - 1.0 + epsr )  <  zeta  <  ( 1.0 - epsr )
           !
           zeta = SIGN( MIN( ABS( zeta ), ( 1.D0 - epsr ) ) , zeta )
           !
           grh2 = (grho (1, k, 1) + grho (1, k, 2) ) **2 + &
                  (grho (2, k, 1) + grho (2, k, 2) ) **2 + &
                  (grho (3, k, 1) + grho (3, k, 2) ) **2
           call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
        else
           sc = 0.d0
           v1cup = 0.d0
           v1cdw = 0.d0
           v2c = 0.d0
        endif
        !
        ! first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        v (k, 1) = v (k, 1) + e2 * (v1xup + v1cup)
        v (k, 2) = v (k, 2) + e2 * (v1xdw + v1cdw)
        !
        ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        do ipol = 1, 3
           h (ipol, k, 1) = e2 * ( (v2xup + v2c) * grho (ipol, k, 1) &
                + v2c * grho (ipol, k, 2) )
           h (ipol, k, 2) = e2 * ( (v2xdw + v2c) * grho (ipol, k, 2) &
                + v2c * grho (ipol, k, 1) )
        enddo
        vtxcgc = vtxcgc + e2 * (v1xup + v1cup) * (rho (k, 1) - &
             rho_core (k) * fac)
        vtxcgc = vtxcgc + e2 * (v1xdw + v1cdw) * (rho (k, 2) - &
             rho_core (k) * fac)
        etxcgc = etxcgc + e2 * (sx + sc)
     endif
  enddo
  do is = 1, nspin
     call DAXPY (nrxx, - fac, rho_core, 1, rho (1, is), 1)
  enddo
  deallocate(grho)
  allocate (dh( nrxx))    
  !
  ! second term of the gradient correction :
  ! \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  do is = 1, nspin
     call grad_dot (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, h (1, 1, is), &
          ngm, g, nl, alat, dh)
     do k = 1, nrxx
        v (k, is) = v (k, is) - dh (k)
        vtxcgc = vtxcgc - dh (k) * rho (k, is)
     enddo
  enddo

  vtxc = vtxc + omega * vtxcgc / (nr1 * nr2 * nr3)
  etxc = etxc + omega * etxcgc / (nr1 * nr2 * nr3)

  deallocate (dh)
  deallocate (h)
  return

end subroutine gradcorr
!--------------------------------------------------------------------

subroutine gradient (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, a, &
     ngm, g, nl, alat, ga)
  !--------------------------------------------------------------------
  !
  ! Calculates ga = \grad a in R-space (a is also in R-space)
  USE kinds
  use gvect, only: nlm
  use wvfct, only: gamma_only
  implicit none

  integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, ngm, nl (ngm)
  real(kind=DP) :: a (nrxx), g (3, ngm), ga (3, nrxx), alat
  integer :: n, ipol
  real(kind=DP), allocatable :: aux (:,:), gaux (:,:)
  real(kind=DP) ::  tpi, tpiba
  parameter (tpi = 2.d0 * 3.14159265358979d0)

  allocate (aux( 2,nrxx))    
  allocate (gaux(2,nrxx))    

  tpiba = tpi / alat
  !
  ! copy a(r) to complex array...
  !
  aux(2,:) = 0.d0
  call DCOPY (nrxx, a, 1, aux, 2)
  !
  ! bring a(r) to G-space, a(G) ...
  !
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  ! multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  ga(:,:) = 0.d0

  do ipol = 1, 3
     gaux(:,:) = 0.d0
     do n = 1, ngm
        gaux (1, nl (n) ) = - g (ipol, n) * aux (2, nl (n) )
        gaux (2, nl (n) ) =   g (ipol, n) * aux (1, nl (n) )
     enddo
     if (gamma_only) then
        do n = 1, ngm
           gaux (1, nlm(n) ) =   gaux (1, nl(n) )
           gaux (2, nlm(n) ) = - gaux (2, nl(n) )
        enddo
     end if
     !
     ! bring back to R-space, (\grad_ipol a)(r) ...
     !
     call cft3 (gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !

     call DAXPY (nrxx, tpiba, gaux, 2, ga (ipol, 1), 3)

  enddo
  deallocate (gaux)
  deallocate (aux)
  return

end subroutine gradient
!--------------------------------------------------------------------

subroutine grad_dot (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, a, &
     ngm, g, nl, alat, da)
  !--------------------------------------------------------------------
  !
  ! Calculates da = \sum_i \grad_i a_i in R-space
  USE kinds
  use gvect, only: nlm
  use wvfct, only: gamma_only
  implicit none
  integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, ngm, nl (ngm)
  real(kind=DP) :: a (3, nrxx), g (3, ngm), da (nrxx), alat
  integer :: n, ipol
  real(kind=DP), allocatable :: aux (:,:), gaux (:,:)
  real(kind=DP) ::  tpi, tpiba
  parameter (tpi = 2.d0 * 3.14159265358979d0)

  allocate (aux( 2,nrxx))    
  allocate (gaux(2,nrxx))    

  gaux(:,:) = 0.d0
  do ipol = 1, 3
     !
     ! copy a(ipol,r) to a complex array...
     !
     aux(2,:) = 0.d0
     call DCOPY (nrxx, a (ipol, 1), 3, aux, 2)
     !
     ! bring a(ipol,r) to G-space, a(G) ...
     !
     call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     !
     ! multiply by (iG) to get (\grad_ipol a)(G) ...
     !
     do n = 1, ngm
        gaux (1, nl (n) ) = gaux (1, nl (n) ) - g (ipol, n) * aux (2,nl(n))
        gaux (2, nl (n) ) = gaux (2, nl (n) ) + g (ipol, n) * aux (1,nl(n))
     enddo
  enddo
  if (gamma_only) then
     do n = 1, ngm
        gaux (1, nlm(n) ) =   gaux (1, nl (n) )
        gaux (2, nlm(n) ) = - gaux (2, nl (n) )
     enddo
  end if
  !
  !  bring back to R-space, (\grad_ipol a)(r) ...
  !
  call cft3 (gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  ! ...add the factor 2\pi/a  missing in the definition of G and sum
  !
  tpiba = tpi / alat
  do n=1,nrxx
     da(n) = gaux(1,n)*tpiba
  end do
  !
  deallocate (gaux)
  deallocate (aux)
  return
end subroutine grad_dot

