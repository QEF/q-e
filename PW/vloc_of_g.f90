!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine vloc_of_g (lloc, lmax, numeric, mesh, msh, rab, r, vnl, &
     cc, alpc, nlc, nnl, zp, aps, alps, tpiba2, ngl, gl, omega, vloc)
  !----------------------------------------------------------------------
  !
  !    This routine computes the Fourier transform of the local
  !    part of the pseudopotential. Two types of local potentials
  !    are allowed:
  !
  !    a) The pseudopotential is in analytic form and its fourier
  !       transform is computed analytically
  !    b) The pseudopotential is in numeric form and its fourier
  !       transform is computed numerically
  !
  !    The local pseudopotential of the US case is always in
  !    numerical form, expressed in Ry units.
  !
#include "machine.h"
  use parameters
  implicit none
  !
  !    first the dummy variables
  !
  integer :: nlc, nnl, ngl, lloc, lmax, mesh, msh
  ! input: analytic, number of erf functions
  ! input: analytic, number of gaussian functions
  ! input: the number of shell of G vectors
  ! input: the l taken as local part
  ! input: the maximum non local angular momentum
  ! input: numeric, the dimensions of the mesh
  ! input: numeric, number of mesh points for radial integration

  real(kind=DP) :: cc (2), alpc (2), alps (3, 0:3), aps (6, 0:3), &
       zp, rab (mesh), r (mesh), vnl (mesh), tpiba2, omega, gl (ngl), &
       vloc (ngl)
  ! input: analytic, c of the erf functions
  ! input: analytic, alpha of the erf
  ! input: analytic, alpha of the gaussians
  ! input: analytic, a and b of the gaussians
  ! input: valence pseudocharge
  ! input: numeric, the derivative of mesh points
  ! input: numeric, the mesh points
  ! input: numeric, the pseudo on the radial mesh
  ! input: 2 pi / alat
  ! input: the volume of the unit cell
  ! input: the moduli of g vectors for each shekk
  ! output: the fourier transform of the potential
  logical :: numeric
  ! input: if true the pseudo is numeric
  !
  real(kind=DP), parameter :: pi = 3.14159265358979d0, fpi= 4.d0 * pi, &
                              e2 = 2.d0, eps= 1.d-8
  !    local variables
  !
  real(kind=DP) :: vlcp, fac, den1, den2, g2a, erf, gx
  real(kind=DP), allocatable :: aux (:), aux1 (:)
  !  auxiliary variables
  integer :: i, igl, igl0, l, ir
  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shells with g != 0
  ! the angular momentum
  ! counter on mesh points

  if (.not.numeric) then

     vloc(:) = 0.d0
     do i = 1, nlc
        if (gl (1) .lt.eps) then
           !
           !    This is the G=0 component of the local potential
           !    giving rise to the so-called "alpha*Z" term in the energy
           !
           vloc (1) = vloc (1) + cc (i) * tpiba2 * 0.25d0 / alpc (i)
           igl0 = 2
        else
           igl0 = 1
        endif
        !
        !    here the G<>0 terms
        !
        den1 = 0.25d0 * tpiba2 / alpc (i)
        do igl = igl0, ngl
           vlcp = - cc (i) * exp ( - gl (igl) * den1) / gl (igl)
           vloc (igl) = vloc (igl) + vlcp
        enddo
     enddo
     den1 = zp * e2 * fpi / tpiba2 / omega
     vloc(:) = vloc (:) * den1
     !
     ! Add the local part l=lloc term (only if l <= lmax)
     !
     l = lloc
     if (l.le.lmax) then
        do i = 1, nnl
           fac = (pi / alps (i, l) ) **1.5d0 * e2 / omega
           den1 = aps (i + 3, l) / alps (i, l)
           den2 = 0.25d0 * tpiba2 / alps (i, l)
           if (gl (1) .lt.eps) then
              !
              !      first the G=0 component
              !
              vloc (1) = vloc (1) + fac * (aps (i, l) + den1 * 1.5d0)
              igl0 = 2
           else
              igl0 = 1
           endif
           !
           !   and here all the other g components
           !
           do igl = igl0, ngl
              g2a = gl (igl) * den2
              vlcp = fac * exp ( - g2a) * (aps (i, l) + den1 * (1.5d0 - g2a) )
              vloc (igl) = vloc (igl) + vlcp
           enddo
        enddo
     endif
  else
     !
     ! Pseudopotentials in numerical form (Vnl(lloc) contain the local part)
     ! in order to perform the Fourier transform, a term erf(r)/r is
     ! subtracted in real space and added again in G space
     !
     allocate ( aux(mesh), aux1(mesh) )
     if (gl (1) .lt.eps) then
        !
        ! first the G=0 term
        !
        do ir = 1, msh
           aux (ir) = r (ir) * (r (ir) * vnl (ir) + zp * e2)
        enddo
        call simpson (msh, aux, rab, vlcp)
        vloc (1) = vlcp        
        igl0 = 2
     else
        igl0 = 1
     endif
     !
     !   here the G<>0 terms, we first compute the part of the integrand func
     !   indipendent of |G| in real space
     !
     do ir = 1, msh
        aux1 (ir) = r (ir) * vnl (ir) + zp * e2 * erf (r (ir) )
     enddo
     fac = zp * e2 / tpiba2
     !
     !    and here we perform the integral, after multiplying for the |G|
     !    dependent  part
     !
     do igl = igl0, ngl
        gx = sqrt (gl (igl) * tpiba2)
        do ir = 1, msh
           aux (ir) = aux1 (ir) * sin (gx * r (ir) ) / gx
        enddo
        call simpson (msh, aux, rab, vlcp)
        !
        !     here we add the analytic fourier transform of the erf function
        !
        vlcp = vlcp - fac * exp ( - gl (igl) * tpiba2 * 0.25d0) &
             / gl (igl)
        vloc (igl) = vlcp
     enddo
     vloc (:) = vloc(:) * fpi / omega
     deallocate (aux, aux1)
  endif
  return
end subroutine vloc_of_g

