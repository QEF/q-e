!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_ramtns (iudyn, ramtns)
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE constants, ONLY : fpi, BOHR_RADIUS_ANGS
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat
  USE control_ph, ONLY : xmldyn
  !
  implicit none
  integer, intent(in) :: iudyn  ! unit
  real(DP), intent(in) :: ramtns(3, 3, 3, nat) ! raman tensor

  ! local variables
  integer :: na, ic, jc, kc
  ! counters
  real (DP), parameter ::   convfact = BOHR_RADIUS_ANGS**2
  ! conversion factor from au^2 to A^2
  !
  ! write raman tensor (D chi/d tau in A^2) to iudyn
  !
  IF (xmldyn) RETURN
  write(iudyn,'(/5x,"Raman tensor (A^2)",/)')
  do na = 1, nat
     do kc = 1, 3
        write (iudyn,'(5x,"atom # ",i4,"    pol.",i3)') na, kc
        write (iudyn, '(3e24.12)') ( (ramtns(ic, jc, kc, na) * &
              omega/fpi*convfact, ic = 1, 3), jc = 1, 3)
     enddo
  enddo

  return
end subroutine write_ramtns
