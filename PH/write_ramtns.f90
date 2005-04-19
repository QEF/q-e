!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_ramtns (ramtns, nat, iudyn)
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE cell_base, ONLY : omega
  USE constants, ONLY : fpi
  implicit none
  ! input variables
  integer :: iudyn, nat
  ! unit number
  ! number of atom in the unit cell

  real(kind=DP) :: ramtns(3, 3, 3, nat)
  ! the raman tensor

  ! local variables
  integer :: na, ic, jc, kc
  ! counter on atoms
  ! cartesian coordinate counters
  real (kind=DP), parameter ::   convfact = 0.529177**2
  ! conversion factor from au^2 to A^2
  !
  ! write raman tensor (D chi/d tau in A^2) to iudyn
  !
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
