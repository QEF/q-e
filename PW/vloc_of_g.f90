!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine vloc_of_g (mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, &
     gl, omega, vloc)
  !----------------------------------------------------------------------
  !
  !    This routine computes the Fourier transform of the local
  !    part of the pseudopotential.
  !
  !    The local pseudopotential of the US case is always in
  !    numerical form, expressed in Ry units.
  !
#include "f_defs.h"
  USE kinds
  USE constants, ONLY : pi, fpi, e2, eps8
  implicit none
  !
  !    first the dummy variables
  !
  integer :: ngl, mesh, msh
  ! input: the number of shells of G vectors
  ! input: the dimensions of the mesh
  ! input: number of mesh points for radial integration

  real(DP) :: zp, rab (mesh), r (mesh), vloc_at (mesh), tpiba2, omega, &
       gl (ngl), vloc (ngl)
  ! input: valence pseudocharge
  ! input: the derivative of mesh points
  ! input: the mesh points
  ! input: the pseudo on the radial mesh
  ! input: 2 pi / alat
  ! input: the volume of the unit cell
  ! input: the moduli of g vectors for each shell
  ! output: the fourier transform of the potential
  !
  !    local variables
  !
  real(DP) :: vlcp, fac, den1, den2, g2a, gx
  real(DP), allocatable :: aux (:), aux1 (:)
  !  auxiliary variables
  real(DP), external :: erf
  integer :: i, igl, igl0, l, ir
  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shells with g != 0
  ! the angular momentum
  ! counter on mesh points

  ! Pseudopotentials in numerical form (Vloc_at) contain the local part)
  ! in order to perform the Fourier transform, a term erf(r)/r is
  ! subtracted in real space and added again in G space
  !
  allocate ( aux(mesh), aux1(mesh) )
  if (gl (1) < eps8) then
     !
     ! first the G=0 term
     !
     do ir = 1, msh
        aux (ir) = r (ir) * (r (ir) * vloc_at (ir) + zp * e2)
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
     aux1 (ir) = r (ir) * vloc_at (ir) + zp * e2 * erf (r (ir) )
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

return
end subroutine vloc_of_g

