!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine set_drhoc (q)
  !---------------------------------------------------------------------
  ! calculate the fourier trasform of the core charge for all pseudo
  ! without structure factor and put it in drc, at q point
  ! used to calculate derivatives of the core charge
  !
#include "f_defs.h"
  USE ions_base, ONLY : ntyp => nsp
  use pwcom
  USE radial_grids, ONLY: ndmx
  USE atom, ONLY : nlcc, msh, rgrid, rho_atc
  USE kinds, only : DP
  use phcom
  implicit none
  !
  !     here the local variables
  !
  real(DP) :: gq2, g2a, gx, rhocgnt, rhocgip, q (3)

  integer :: ir, ng, nt
  ! counter on radial mesh points
  ! counter on plane waves
  ! counter on atomic types
  real(DP), allocatable :: aux (:)

  call start_clock('set_drhoc')
  !
  drc (:,:) = (0.d0, 0.d0)
  if ( nlcc_any ) go to 10
  !
  allocate (aux ( ndmx))
  do ng = 1, ngm
     gq2 = (g (1, ng) + q (1) ) **2 + (g (2, ng) + q (2) ) **2 + &
          (g (3, ng) + q (3) ) **2
     gq2 = gq2 * tpiba2
     do nt = 1, ntyp
        rhocgnt = 0.d0
        if (nlcc (nt) ) then
           if (gq2 < 1.0d-8) then
              do ir = 1, msh (nt)
                 aux (ir) = rgrid(nt)%r(ir) **2 * rho_atc (ir, nt)
              enddo
              call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgip)
           else
              gx = sqrt (gq2)
              call sph_bes (msh (nt), rgrid(nt)%r, gx, 0, aux)
              do ir = 1, msh (nt)
                 aux (ir) = rgrid(nt)%r(ir) **2 * rho_atc (ir, nt) * aux (ir)
              enddo
              call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgip)
           endif
           rhocgnt = rhocgip * fpi
        endif
        drc (ng, nt) = rhocgnt / omega
     enddo

  enddo

  deallocate(aux)
10  call stop_clock('set_drhoc')
  return
end subroutine set_drhoc
