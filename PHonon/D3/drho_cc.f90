!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine drho_cc (iflag)
  !-----------------------------------------------------------------------
  !
  !   Used when non_linear_core_correction are present to change the files
  !   containing the variation of the charge
  !   iflag = +1 :
  !       adds the variation of the core charge to the variation of the
  !       valence charge ( both for xq.eq.0 and xq.ne.0 )
  !
  !   iflag = -1 :
  !       subtracts the variation of the core charge to the variation of
  !       the total charge --used to set drho and d0rho as they were
  !       before the first call of drho_cc--
  !
  USE kinds, only : DP
  use pwcom
  use phcom
  use d3com

  use qpoint,     ONLY : xq
  use control_lr, ONLY : lgamma
  use uspp,       ONLY : nlcc_any

  implicit none
  integer :: iflag
  real (DP) :: xq0 (3), scalef

  if (.not.nlcc_any) return
  if (iflag.eq. - 1) then
     iudrho = iudrho-1000
     iud0rho=iud0rho-1000
     RETURN
  else
     scalef = 1.d0
  end if
  xq0 = 0.d0
  call drho_drc (iud0rho, ug0, xq0, d0rc, scalef)
  if (.not.lgamma) call drho_drc (iudrho, u, xq, drc, scalef)

  return
end subroutine drho_cc
