!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvloc_of_g (mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, gl, &
     omega, dvloc)
  !----------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds
  USE constants , ONLY : pi, fpi, e2, eps8
  implicit none
  !
  !    first the dummy variables
  !
  integer :: ngl, mesh, msh
  ! input: the number of shell of G vectors
  ! input: the maximum non local angular mome
  ! input: numeric, the dimensions of the mes
  ! number of mesh points for radial integrat

  real(DP) :: zp, rab (mesh), r (mesh), vloc_at (mesh), tpiba2, omega,&
       gl (ngl), dvloc (ngl)
  ! input: valence pseudocharge
  ! input: the derivative of the radial grid
  ! input: the radial grid
  ! input: the pseudo on the radial grid
  ! input: 2 pi / alat
  ! input: the volume of the unit cell
  ! input: the moduli of g vectors for each s
  ! output: the fourier transform dVloc/dG
  !
  !    and the local variables
  !
  real(DP) :: vlcp, g2a, gx
  real(DP), allocatable ::  aux (:), aux1 (:)
  real(DP), external ::  erf

  integer :: i, igl, igl0, l
  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shells with g != 0
  ! the angular momentum

  ! the  G=0 component is not computed
  if (gl (1) < eps8) then
     dvloc (1) = 0.0d0
     igl0 = 2
  else
     igl0 = 1
  endif

  ! Pseudopotentials in numerical form (Vloc contains the local part)
  ! In order to perform the Fourier transform, a term erf(r)/r is
  ! subtracted in real space and added again in G space
  allocate (aux( mesh))    
  allocate (aux1( mesh))    
  !
  !   This is the part of the integrand function
  !   indipendent of |G| in real space
  !
  do i = 1, msh
     aux1 (i) = r (i) * vloc_at (i) + zp * e2 * erf (r (i) )
  enddo
  do igl = igl0, ngl
     gx = sqrt (gl (igl) * tpiba2)
     !
     !    and here we perform the integral, after multiplying for the |G|
     !    dependent  part
     !
     ! DV(g)/Dg = Integral of r (Dj_0(gr)/Dg) V(r) dr
     do i = 1, msh
        aux (i) = aux1 (i) * (r (i) * cos (gx * r (i) ) / gx - sin (gx &
             * r (i) ) / gx**2)
     enddo
     call simpson (msh, aux, rab, vlcp)
     ! DV(g^2)/Dg^2 = (DV(g)/Dg)/2g
     vlcp = fpi / omega / 2.0d0 / gx * vlcp
     ! subtract the long-range term
     g2a = gl (igl) * tpiba2 / 4.d0
     vlcp = vlcp + fpi / omega * zp * e2 * exp ( - g2a) * (g2a + &
          1.d0) / (gl (igl) * tpiba2) **2
     dvloc (igl) = vlcp
  enddo
  deallocate (aux1)
  deallocate (aux)

  return
end subroutine dvloc_of_g

