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
  use kinds, only : DP
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

  !
  ! write raman tensor on iudyn
  !
  write(iudyn,'(/5x,"Raman tensor",/)')
  do na = 1, nat
     write (iudyn,'(5x,"atom # ",i4)') na
     do kc = 1, 3
        write (iudyn, '(3e24.12)') ( (ramtns(ic, jc, kc, na), &
              ic = 1, 3), jc = 1, 3)
        write(iudyn, '(10x)')
     enddo
  enddo

  return
end subroutine write_ramtns
