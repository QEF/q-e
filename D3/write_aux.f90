!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine write_aux (isw)
  !-----------------------------------------------------------------------
  !
  ! Writes on files partial computation of d3dyn
  !
  USE ions_base,  ONLY : nat
  use pwcom
  use phcom
  use d3com
!
  implicit none
  integer :: isw
!
  if (isw.eq.1) then
     d3dyn_aux1 = (0.0_dp,0.0_dp)
     d3dyn_aux2 = (0.0_dp,0.0_dp)
     d3dyn_aux3 = (0.0_dp,0.0_dp)
     d3dyn_aux4 = (0.0_dp,0.0_dp)
     d3dyn_aux5 = (0.0_dp,0.0_dp)
     d3dyn_aux6 = (0.0_dp,0.0_dp)
     d3dyn_aux7 = (0.0_dp,0.0_dp)
     d3dyn_aux8 = (0.0_dp,0.0_dp)
     d3dyn_aux9 = (0.0_dp,0.0_dp)
  elseif (isw.eq.2) then
     call zcopy (27 * nat * nat * nat, d3dyn, 1, d3dyn_aux4, 1)
     call daxpy (2 * 27 * nat * nat * nat, - 1.d0, d3dyn_aux5, 1, &
          d3dyn_aux4, 1)
     call daxpy (2 * 27 * nat * nat * nat, - 1.d0, d3dyn_aux6, 1, &
          d3dyn_aux4, 1)
     call daxpy (2 * 27 * nat * nat * nat, - 1.d0, d3dyn_aux7, 1, &
          d3dyn_aux4, 1)
     call daxpy (2 * 27 * nat * nat * nat, - 1.d0, d3dyn_aux8, 1, &
          d3dyn_aux4, 1)
     call daxpy (2 * 27 * nat * nat * nat, - 1.d0, d3dyn_aux9, 1, &
          d3dyn_aux4, 1)
     call daxpy (2 * 27 * nat * nat * nat, - 1.d0, d3dyn_aux1, 1, &
          d3dyn_aux4, 1)
     call writed3dyn_5 (d3dyn_aux1, 'd3mat.1', - 1)
     call writed3dyn_5 (d3dyn_aux4, 'd3mat.4', - 1)
     call writed3dyn_5 (d3dyn_aux5, 'd3mat.5', - 1)
     call writed3dyn_5 (d3dyn_aux6, 'd3mat.6', - 1)
     call writed3dyn_5 (d3dyn_aux7, 'd3mat.7', - 1)
     call writed3dyn_5 (d3dyn_aux8, 'd3mat.8', - 1)
     call writed3dyn_5 (d3dyn_aux9, 'd3mat.9', - 1)
     call writed3dyn_5 (d3dyn, 'd3mat.ns', - 1)
  elseif (isw.eq.3) then
     call writed3dyn_5 (d3dyn, 'd3mat.sy', 1)

  endif
  return
end subroutine write_aux
