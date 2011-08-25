!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_dyn_on_file (xq, phi, nat, iudyn)
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  implicit none
  ! input variables
  integer :: iudyn, nat
  ! unit number
  ! number of atom in the unit cell
  complex(DP) :: phi (3, 3, nat, nat)
  !  the dynamical matrix
  real(DP) :: xq (3)
  ! the q vector
  ! local variables

  integer :: na, nb, icar, jcar
  ! counters on atoms
  ! cartesian coordinate counters
  write (iudyn, 9000) (xq (icar), icar = 1, 3)
  do na = 1, nat
     do nb = 1, nat
        write (iudyn, '(2i5)') na, nb
        do icar = 1, 3
!           write (iudyn, '(3e24.12)') (phi(icar,jcar,na,nb), jcar=1,3)
           write (iudyn, '(3(2f12.8,2x))') (phi(icar,jcar,na,nb), jcar=1,3)
        enddo
     enddo
  enddo

  return
9000 format(/,5x,'Dynamical  Matrix in cartesian axes', &
       &       //,5x,'q = ( ',3f14.9,' ) ',/)
end subroutine write_dyn_on_file
