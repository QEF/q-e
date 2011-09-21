!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine setlocq (xq, mesh, msh, rab, r, vloc_at, zp, tpiba2, ngm, &
     g, omega, vloc)
  !----------------------------------------------------------------------
  !
  !    This routine computes the Fourier transform of the local
  !    part of the pseudopotential in the q+G vectors.
  !
  !    The local pseudopotential of the US case is always in
  !    numerical form, expressed in Ry units.
  !
  USE kinds, only  : DP
  USE constants, ONLY : e2, fpi, pi
  !
  implicit none
  !
  !    first the dummy variables
  !
  integer :: ngm, mesh, msh
  ! input: the number of G vectors
  ! input: the dimensions of the mesh
  ! input: mesh points for radial integration

  real(DP) :: xq (3), zp, rab (mesh), r (mesh), vloc_at(mesh), tpiba2,&
       omega, g (3, ngm), vloc (ngm)
  ! input: the q point
  ! input: valence pseudocharge
  ! input: the derivative of mesh points
  ! input: the mesh points
  ! input: the pseudo on the radial
  ! input: 2 pi / alat
  ! input: the volume of the unit cell
  ! input: the g vectors coordinates
  ! output: the fourier transform of the potential
  !
  !    and the local variables
  !
  real(DP), parameter :: eps = 1.d-8
  real(DP) :: vlcp, vloc0, fac, g2a, aux (mesh), &
       aux1 (mesh), gx
  ! auxiliary variables
  ! gx = modulus of g vectors
  real(DP), external :: qe_erf
  ! the erf function
  integer :: ig, ir
  ! counters
  !
  ! Pseudopotentials in numerical form (Vnl(lloc) contain the local part)
  ! in order to perform the Fourier transform, a term erf(r)/r is
  ! subtracted in real space and added again in G space
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
     aux1 (ir) = r (ir) * vloc_at (ir) + zp * e2 * qe_erf (r (ir) )
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

  return
end subroutine setlocq

!----------------------------------------------------------------------
subroutine setlocq_coul (xq, zp, tpiba2, ngm, g, omega, vloc)
 !----------------------------------------------------------------------
 !
 !    Fourier transform of the Coulomb potential - For all-electron
 !    calculations, in specific cases only, for testing purposes
 !
 USE kinds, ONLY: DP
 USE constants, ONLY : fpi, e2, eps8
 implicit none
 !
 integer, intent(in) :: ngm
 real(DP) :: xq (3), zp, tpiba2, omega, g(3,ngm)
 real(DP), intent (out) :: vloc(ngm)
 !
 real(DP) :: g2a
 integer :: ig

 do ig = 1, ngm
  g2a = (xq (1) + g (1, ig) ) **2 + (xq (2) + g (2, ig) ) **2 + &
        (xq (3) + g (3, ig) ) **2
  if (g2a < eps8) then
       vloc (ig) = 0.d0
  else
       vloc (ig) = - fpi * zp *e2 / omega / tpiba2 / g2a
  endif
 enddo

end subroutine setlocq_coul
