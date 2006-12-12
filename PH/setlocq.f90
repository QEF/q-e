!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine setlocq (xq, lloc, lmax, numeric, mesh, msh, rab, r, &
     vloc_at, cc, alpc, nlc, nnl, zp, aps, alps, tpiba2, ngm, g, omega, &
     vloc)
  !----------------------------------------------------------------------
  !
  !    This routine computes the Fourier transform of the local
  !    part of the pseudopotential in the q+G vectors.
  !    Two types of local potentials are allowed:
  !
  !    a) The pseudopotential is in analytic form and its fourier
  !       transform is computed analytically
  !    b) The pseudopotential is in numeric form and its fourier
  !       transform is computed numerically
  !
  !    The local pseudopotential of the US case is always in
  !    numerical form, expressed in Ry units.
  !
  !    Last revision 5 oct. 1995 by Andrea Dal Corso
  !
  !
#include "f_defs.h"
  USE kinds, only  : DP
  USE constants, ONLY : e2, fpi, pi
  !
  implicit none
  !
  !    first the dummy variables
  !
  integer :: nlc, nnl, ngm, lloc, lmax, mesh, msh
  ! input: analytic, number of erf functions
  ! input: analytic, number of gaussian funct
  ! input: the number of G vectors
  ! input: the non local part which becomes l
  ! input: the maximum non local angular mome
  ! input: numeric, the dimensions of the mes
  ! input: mesh points for radial integration

  real(DP) :: xq (3), cc (2), alpc (2), alps (3, 0:3), aps (6, 0:3), &
       zp, rab (mesh), r (mesh), vloc_at(mesh), tpiba2, omega, g (3, ngm), &
       vloc (ngm)
  ! input: the q point
  ! input: analytic, c of the erf functions
  ! input: analytic, alpha of the erf
  ! input: analytic, alpha of the gaussians
  ! input: analytic, a and b of the gaussians
  ! input: valence pseudocharge
  ! input: numeric, the derivative of mesh po
  ! input: numeric, the mesh points
  ! input: numeric, the pseudo on the radial
  ! input: 2 pi / alat
  ! input: the volume of the unit cell
  ! input: the g vectors coordinates
  ! output: the fourier transform of the pote
  logical :: numeric
  ! input: if true the pseudo is numeric
  !
  !
  !    and the local variables
  !
  real(DP), parameter :: eps = 1.d-8
  real(DP) :: vlcp, vloc0, fac, den1, den2, g2a, g2a1, aux (mesh), &
       aux1 (mesh), gx
  ! auxiliary variables
  ! gx = modulus of g vectors
  real(DP), external :: erf
  ! the erf function
  integer :: i, ig, l, ir
  ! counters
  !
  if (.not.numeric) then
     vloc (:) = 0.d0
     !
     !    In this case the potential is given in analytic form
     !
     do ig = 1, ngm
        g2a = (xq (1) + g (1, ig) ) **2 + (xq (2) + g (2, ig) ) **2 + &
             (xq (3) + g (3, ig) ) **2
        if (g2a < eps) then
           do i = 1, nlc
              vloc (ig) = vloc (ig) + cc (i) * tpiba2 * 0.25d0 / alpc (i)
           enddo
        else
           do i = 1, nlc
              den1 = 0.25d0 * tpiba2 / alpc (i)
              vlcp = - cc (i) * exp ( - g2a * den1) / g2a
              vloc (ig) = vloc (ig) + vlcp
           enddo
        endif
     enddo
     den1 = zp * e2 * fpi / tpiba2 / omega
     vloc (:) = vloc(:) * den1
     !
     ! Add the local part l=lloc term (only if l <= lmax)
     !
     l = lloc

     if (l <= lmax) then
        !
        !   and here all the other g components
        !
        do ig = 1, ngm
           g2a = ( (xq (1) + g (1, ig) ) **2 + (xq (2) + g (2, ig) ) ** &
                2 + (xq (3) + g (3, ig) ) **2)
           if (g2a < eps) then
              do i = 1, nnl
                 fac = (pi / alps (i, l) ) **1.5d0 * e2 / omega
                 den1 = aps (i + 3, l) / alps (i, l)
                 vloc (ig) = vloc (ig) + fac * (aps (i, l) + den1 * 1.5d0)
              enddo
           else
              do i = 1, nnl
                 den1 = aps (i + 3, l) / alps (i, l)
                 den2 = 0.25d0 * tpiba2 / alps (i, l)
                 fac = (pi / alps (i, l) ) **1.5d0 * e2 / omega
                 g2a1 = g2a * den2
                 vlcp = fac * exp ( - g2a1) * (aps (i, l) + den1 * &
                      (1.5d0 - g2a1) )
                 vloc (ig) = vloc (ig) + vlcp
              enddo
           endif
        enddo
     endif
  else
     !
     ! Pseudopotentials in numerical form (Vnl(lloc) contain the local part)
     ! in order to perform the Fourier transform, a term erf(r)/r is
     ! subtracted in real space and added again in G space
     !
     !
     ! first the G=0 term
     !
     do ir = 1, msh
        aux (ir) = r (ir) * (r (ir) * vloc_at (ir) + zp * e2)
     enddo
     call simpson (msh, aux, rab, vloc0)
     !
     !   here the G<>0 terms, we first compute the part of the integrand func
     !   indipendent of |G| in real space
     !
     do ir = 1, msh
        aux1 (ir) = r (ir) * vloc_at (ir) + zp * e2 * erf (r (ir) )
     enddo
     fac = zp * e2 / tpiba2
     !
     !    and here we perform the integral, after multiplying for the |G|
     !    dependent  part
     !
     do ig = 1, ngm
        g2a = (xq (1) + g (1, ig) ) **2 + (xq (2) + g (2, ig) ) **2 + &
             (xq (3) + g (3, ig) ) **2
        if (g2a < eps) then
           vloc (ig) = vloc0
        else
           gx = sqrt (g2a * tpiba2)
           do ir = 1, msh
              aux (ir) = aux1 (ir) * sin (gx * r (ir) ) / gx
           enddo
           call simpson (msh, aux, rab, vlcp)
           !
           !     here we add the analytic fourier transform of the erf function
           !
           vlcp = vlcp - fac * exp ( - g2a * tpiba2 * 0.25d0) / g2a
           vloc (ig) = vlcp
        endif
     enddo

     vloc(:) = vloc(:) * fpi / omega
  endif
  return
end subroutine setlocq
