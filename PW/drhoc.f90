!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine drhoc (ngl, gl, omega, tpiba2, numeric, a_nlcc, b_nlcc, &
     alpha_nlcc, mesh, r, rab, rhoc, rhocg)
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: ngl, mesh
  ! input: the number of g shell
  ! input: the number of radial mesh points

  real(DP) :: gl (ngl), r (mesh), rab (mesh), rhoc (mesh), omega, &
       tpiba2, a_nlcc, b_nlcc, alpha_nlcc, rhocg (ngl)
  ! input: the number of G shells
  ! input: the radial mesh
  ! input: the derivative of theradial mesh
  ! input: the radial core charge
  ! input: the volume of the unit cell
  ! input: 2 times pi / alat
  ! input: the a_c of the analitycal form
  ! input: the b_c of the analitical form
  ! input: the alpha of the analytical form
  ! output: the fourier transform of the co
  logical :: numeric
  ! input: if true the charge is in numeric
  !
  !     two parameters
  !
  real(DP), parameter :: pi = 3.14159265358979d0, fpi = 4.d0 * pi
  !
  !     here the local variables
  !
  real(DP) :: gx, g2a, rhocg1
  real(DP), allocatable ::  aux (:)
  ! the modulus of g for a given shell
  ! the argument of the exponential
  ! the fourier transform
  ! auxiliary memory for integration

  integer :: ir, igl, igl0
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl

  !
  !    here we compute the fourier transform is the charge in numeric form
  !
  if (numeric) then
     allocate (aux( mesh))     
     !
     ! G=0 term
     !
     if (gl (1) < 1.0e-8) then
        do ir = 1, mesh
           aux (ir) = r (ir) **2 * rhoc (ir)
        enddo
        call simpson (mesh, aux, rab, rhocg1)
        rhocg (1) = fpi * rhocg1 / omega
        igl0 = 2
     else
        igl0 = 1
     endif
     !
     ! G <> 0 term
     !
     do igl = igl0, ngl
        gx = sqrt (gl (igl) * tpiba2)
        call sph_bes (mesh, r, gx, 0, aux)
        do ir = 1, mesh
           aux (ir) = r (ir) **2 * rhoc (ir) * aux (ir)
        enddo
        call simpson (mesh, aux, rab, rhocg1)
        rhocg (igl) = fpi * rhocg1 / omega
     enddo
     deallocate(aux)
  else
     !
     !    here the case where the charge is in analytic form
     !
     do igl = 1, ngl
        g2a = gl (igl) * tpiba2 * 0.25d0 / alpha_nlcc
        rhocg (igl) = (pi / alpha_nlcc) **1.5d0 * exp ( - g2a) * &
             (a_nlcc + b_nlcc / alpha_nlcc * (1.5d0 - g2a) ) / omega
     enddo

  endif
  return
end subroutine drhoc

