!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
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
  USE constants, ONLY : fpi
  USE cell_base, ONLY : omega, tpiba2
  USE gvect,     ONLY : g, ngm
  USE ions_base, ONLY : ntyp => nsp
  USE atom, ONLY : msh, rgrid
  USE uspp_param, ONLY : upf
  USE kinds, only : DP
  USE nlcc_ph,  ONLY : nlcc_any, drc
  implicit none
  !
  !     here the local variables
  !
  real(DP) :: gq2, gx, rhocgnt, rhocgip, q (3)

  integer :: ir, ng, nt, ndm
  ! counter on radial mesh points
  ! counter on plane waves
  ! counter on atomic types
  real(DP), allocatable :: aux (:)

  call start_clock('set_drhoc')
  !
  drc (:,:) = (0.d0, 0.d0)
  if ( .not. nlcc_any ) go to 10
  !
  ndm = MAXVAL (msh(1:ntyp))
  allocate (aux ( ndm))
  do ng = 1, ngm
     gq2 = (g (1, ng) + q (1) ) **2 + (g (2, ng) + q (2) ) **2 + &
          (g (3, ng) + q (3) ) **2
     gq2 = gq2 * tpiba2
     do nt = 1, ntyp
        rhocgnt = 0.d0
        if ( upf(nt)%nlcc ) then
           if (gq2 < 1.0d-8) then
              do ir = 1, msh (nt)
                 aux (ir) = rgrid(nt)%r(ir) **2 * upf(nt)%rho_atc (ir)
              enddo
              call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgip)
           else
              gx = sqrt (gq2)
              call sph_bes (msh (nt), rgrid(nt)%r, gx, 0, aux)
              do ir = 1, msh (nt)
                 aux (ir) = rgrid(nt)%r(ir) **2 * upf(nt)%rho_atc(ir) * aux(ir)
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
