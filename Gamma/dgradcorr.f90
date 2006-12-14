!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine dgradcor1 (rho, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, &
     drho, drhoc, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin, &
     nl, nlm, ngm, g, alat, omega, dvxc)
  !     ===================
  !--------------------------------------------------------------------
  !  ADD Gradient Correction contibution to screening potential
  !  phonon calculation, half G-vectors
#include "f_defs.h"
  USE kinds, only : DP
  implicit none

  !
  integer :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, nspin, &
       nl (ngm), nlm(ngm)

  real(DP) :: rho (nrxx, nspin), grho (3, nrxx, nspin), &
       dvxc_rr(nrxx, nspin, nspin), dvxc_sr (nrxx, nspin, nspin), &
       dvxc_ss (nrxx,nspin, nspin), dvxc_s (nrxx, nspin, nspin),&
       drho (nrxx,nspin), g (3, ngm), alat, omega

  complex(DP) :: drhoc(nrxx, nspin), dvxc (nrxx, nspin)
  integer :: k, ipol, is, js, ks, ls
  real(DP) :: epsr, epsg, grho2
  complex(DP) :: s1
  complex(DP) :: a (2, 2, 2), b (2, 2, 2, 2), c (2, 2, 2), &
                      ps (2, 2), ps1 (3, 2, 2), ps2 (3, 2, 2, 2)
  real(DP), allocatable  :: gdrho (:,:,:)
  complex(DP), allocatable :: h (:,:,:), dh (:)
  parameter (epsr = 1.0d-6, epsg = 1.0d-10)

  allocate (gdrho( 3, nrxx , nspin))    
  allocate (h(  3, nrxx , nspin))    
  allocate (dh( nrxx))    

  h (:,:,:) = (0.d0, 0.d0)
  do is = 1, nspin
     call gradient1 (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
         drhoc(1, is), ngm, g, nl, nlm, alat, gdrho (1, 1, is) )
  enddo
  do k = 1, nrxx
     grho2 = grho(1, k, 1)**2 + grho(2, k, 1)**2 + grho(3, k, 1)**2
     if (nspin.eq.1) then
        !
        !    LDA case
        !
        if (abs (rho (k, 1) ) .gt.epsr.and.grho2.gt.epsg) then
           s1 = grho (1, k, 1) * gdrho (1, k, 1) + &
                grho (2, k, 1) * gdrho (2, k, 1) + &
                grho (3, k, 1) * gdrho (3, k, 1)
           !
           ! linear variation of the first term
           !
           dvxc (k, 1) = dvxc (k, 1) + dvxc_rr (k, 1, 1) * drho (k, 1) &
                + dvxc_sr (k, 1, 1) * s1
           do ipol = 1, 3
              h (ipol, k, 1) = (dvxc_sr(k, 1, 1) * drho(k, 1) + &
                                dvxc_ss(k, 1, 1) * s1 )*grho(ipol, k, 1) + &
                                dvxc_s (k, 1, 1) * gdrho (ipol, k, 1)
           enddo
        else
           do ipol = 1, 3
              h (ipol, k, 1) = (0.d0, 0.d0)
           enddo
        endif
     else
        !
        !    LSDA case
        !
        ps (:,:) = (0.d0, 0.d0)
        do is = 1, nspin
           do js = 1, nspin
              do ipol = 1, 3
                 ps1(ipol, is, js) = drho (k, is) * grho (ipol, k, js)
                 ps(is, js) = ps(is, js) + grho(ipol,k,is)*gdrho(ipol,k,js)
              enddo
              do ks = 1, nspin
                 if (is.eq.js.and.js.eq.ks) then
                    a (is, js, ks) = dvxc_sr (k, is, is)
                    c (is, js, ks) = dvxc_sr (k, is, is)
                 else
                    if (is.eq.1) then
                       a (is, js, ks) = dvxc_sr (k, 1, 2)
                    else
                       a (is, js, ks) = dvxc_sr (k, 2, 1)
                    endif
                    if (js.eq.1) then
                       c (is, js, ks) = dvxc_sr (k, 1, 2)
                    else
                       c (is, js, ks) = dvxc_sr (k, 2, 1)
                    endif
                 endif
                 do ipol = 1, 3
                    ps2 (ipol, is, js, ks) = ps (is, js) * grho (ipol, k, ks)
                 enddo
                 do ls = 1, nspin
                    if (is.eq.js.and.js.eq.ks.and.ks.eq.ls) then
                       b (is, js, ks, ls) = dvxc_ss (k, is, is)
                    else
                       if (is.eq.1) then
                          b (is, js, ks, ls) = dvxc_ss (k, 1, 2)
                       else
                          b (is, js, ks, ls) = dvxc_ss (k, 2, 1)
                       endif
                    endif
                 enddo
              enddo
           enddo
        enddo
        do is = 1, nspin
           do js = 1, nspin
              dvxc (k, is) = dvxc (k, is) + dvxc_rr (k, is, js) * drho (k, &
                   js)
              do ipol = 1, 3
                 h (ipol, k, is) = h (ipol, k, is) + &
                      dvxc_s (k, is, js) * gdrho(ipol, k, js)
              enddo
              do ks = 1, nspin
                 dvxc (k, is) = dvxc (k, is) + a (is, js, ks) * ps (js, ks)
                 do ipol = 1, 3
                    h (ipol, k, is) = h (ipol, k, is) + &
                         c (is, js, ks) * ps1 (ipol, js, ks)
                 enddo
                 do ls = 1, nspin
                    do ipol = 1, 3
                       h (ipol, k, is) = h (ipol, k, is) + &
                            b (is, js, ks, ls) * ps2 (ipol, js, ks, ls)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     endif
  enddo
  ! linear variation of the second term
  do is = 1, nspin
     call grad_dot1 (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
          h (1, 1, is), ngm, g, nl, nlm, alat, dh)
     do k = 1, nrxx
        dvxc (k, is) = dvxc (k, is) - dh (k)
     enddo
  enddo
  deallocate (dh)
  deallocate (h)
  deallocate (gdrho)
  return
end subroutine dgradcor1
!
!--------------------------------------------------------------------
subroutine gradient1(nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
     a, ngm, g, nl, nlm, alat, ga)
  !--------------------------------------------------------------------
  ! Calculates ga = \grad a in R-space (a is G-space)
  USE kinds, only : DP
  USE constants, ONLY : tpi
  implicit none
  integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, ngm, nl (ngm), &
        nlm(ngm)
  complex(DP) :: a (nrxx)
  real(DP) :: ga (3, nrxx), g (3, ngm), alat
  integer :: n, ipol
  real(DP) :: tpiba
  complex(DP), allocatable :: gaux (:)

  allocate (gaux(  nrxx))    

  tpiba = tpi / alat

   ! a(G) multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
  !  do ipol = 1, 3
  ! x, y
     ipol=1
     do n = 1, nrxx
        gaux (n) = (0.d0, 0.d0)
     enddo
     do n = 1, ngm
        gaux(nl (n)) = CMPLX(0.d0, g(ipol  , n))* a (nl(n)) - &
                                   g(ipol+1, n) * a (nl(n))
        gaux(nlm(n)) = CMPLX(0.d0, - g(ipol  , n))* CONJG(a (nl(n))) + &
                                     g(ipol+1, n) * CONJG(a (nl(n)))
     enddo
     ! bring back to R-space, (\grad_ipol a)(r) ...

     call cft3 (gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     ! ...and add the factor 2\pi/a  missing in the definition of q+G
     do n = 1, nrxx
        ga (ipol  , n) =  DBLE(gaux (n)) * tpiba
        ga (ipol+1, n) = AIMAG(gaux (n)) * tpiba
     enddo
  ! z
     ipol=3
     do n = 1, nrxx
        gaux (n) = (0.d0, 0.d0)
     enddo
     do n = 1, ngm
        gaux(nl (n)) = CMPLX(0.d0, g(ipol, n)) * a (nl(n))
        gaux(nlm(n)) = CONJG(gaux(nl(n)))
     enddo
     ! bring back to R-space, (\grad_ipol a)(r) ...
     call cft3 (gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     ! ...and add the factor 2\pi/a  missing in the definition of q+G
     do n = 1, nrxx
        ga (ipol, n) =  DBLE(gaux (n)) * tpiba
     enddo
!  enddo
  deallocate (gaux)
  return

end subroutine gradient1
!--------------------------------------------------------------------
subroutine grad_dot1 (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
     a, ngm, g, nl, nlm, alat, da)
  !--------------------------------------------------------------------
  ! Calculates da = \sum_i \grad_i a_i in R-space
  USE kinds, only : DP
  USE constants, ONLY : tpi
  implicit none
  integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, ngm, nl (ngm), &
        nlm(ngm)
  complex(DP) :: a (3, nrxx), da (nrxx)

  real(DP) :: g (3, ngm), alat
  integer :: n, ipol
  real(DP) :: tpiba
  complex(DP), allocatable :: aux (:)
  complex(DP) :: fp, fm, aux1, aux2

  allocate (aux (  nrxx))    

  tpiba = tpi / alat
  do n = 1, nrxx
     da(n) = (0.d0, 0.d0)
  enddo
!!!  do ipol = 1, 3
     ! x, y
     ipol=1
     ! copy a(ipol,r) to a complex array...
     do n = 1, nrxx
        aux (n) = CMPLX( DBLE(a(ipol, n)), DBLE(a(ipol+1, n)))
     enddo
     ! bring a(ipol,r) to G-space, a(G) ...
     call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
     do n = 1, ngm
        fp = (aux(nl (n)) + aux (nlm(n)))*0.5d0
        fm = (aux(nl (n)) - aux (nlm(n)))*0.5d0
        aux1 = CMPLX( DBLE(fp), AIMAG(fm))
        aux2 = CMPLX(AIMAG(fp),- DBLE(fm))
        da (nl(n)) = da (nl(n)) + CMPLX(0.d0, g(ipol  , n)) * aux1 + &
                                  CMPLX(0.d0, g(ipol+1, n)) * aux2
     end do
     ! z
     ipol=3
     ! copy a(ipol,r) to a complex array...
     do n = 1, nrxx
        aux (n) = a(ipol, n)
     enddo
     ! bring a(ipol,r) to G-space, a(G) ...
     call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
     do n = 1, ngm
        da (nl(n)) = da (nl(n)) + CMPLX(0.d0, g(ipol, n)) * aux(nl(n))
     enddo
!!!  enddo
  do n = 1, ngm
     da(nlm(n)) = CONJG(da(nl(n)))
  enddo
  !  bring back to R-space, (\grad_ipol a)(r) ...
  call cft3 (da, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  ! ...add the factor 2\pi/a  missing in the definition of q+G and sum
  do n = 1, nrxx
     da (n) = da (n) * tpiba
  enddo
  deallocate (aux)

  return
end subroutine grad_dot1
