!
! Copyright (C) 2001 PWSCF group
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
  USE atom, ONLY : numeric, nlcc, msh, rgrid, rho_atc
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
  logical :: some_numeric
  ! .true. if some nlcc-pseudo is numeric
  ! in this case an auxiliary vector is needed
  real(DP), allocatable :: aux (:)

  call start_clock('set_drhoc')
  !
  ! check if some nlcc-pseudo is numeric and allocate memory if needed
  !
  some_numeric = .false.
  do nt = 1, ntyp
     some_numeric = some_numeric.or. (nlcc (nt) .and.numeric (nt) )
     if (nlcc (nt) .and. numeric (nt) .and. msh (nt) > ndmx) call errore &
          ('set_drhoc', 'too many mesh points', msh (nt) )
  enddo

  if (some_numeric) allocate (aux ( ndmx))
  drc (:,:) = (0.d0, 0.d0)
  do ng = 1, ngm
     gq2 = (g (1, ng) + q (1) ) **2 + (g (2, ng) + q (2) ) **2 + &
          (g (3, ng) + q (3) ) **2
     gq2 = gq2 * tpiba2
     do nt = 1, ntyp
        rhocgnt = 0.d0
        if (nlcc (nt) ) then
           if (numeric (nt) ) then
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
           else
              g2a = gq2 / 4.d0 / alpha_nlcc (nt)
              rhocgnt = (pi / alpha_nlcc (nt) ) **1.5d0 * exp ( - g2a) &
                   * (a_nlcc (nt) + b_nlcc (nt) / alpha_nlcc (nt) * (1.5d0 - &
                   g2a) )
           endif
        endif
        drc (ng, nt) = rhocgnt / omega
     enddo

  enddo

  if (some_numeric) deallocate(aux)
  call stop_clock('set_drhoc')
  return
end subroutine set_drhoc
