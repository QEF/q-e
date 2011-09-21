!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine write_d3dyn (xq, phi, nat, iudyn, wrmode)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  !
  ! input variables
  !
  integer :: iudyn, nat
  ! unit number
  ! number of atom in the unit cell
  complex (DP) :: phi (3, 3, 3, nat, nat, nat)
  !  derivative of the dynamical matrix
  real (DP) :: xq (3)
  ! the q vector
  logical :: wrmode (3 * nat)
  ! if .true. this mode is to be written
  !
  ! local variables
  !
  integer :: na, nb, nc, icar, jcar, kcar, i
  ! counters on atoms
  ! cartesian coordinate counters
  ! generic counter
  write (iudyn, 9000) (xq (icar), icar = 1, 3)
  do i = 1, 3 * nat
     if (wrmode (i) ) then
        write (iudyn, '(/,12x,"modo:",i5,/)') i
        nc = (i - 1) / 3 + 1
        kcar = i - 3 * (nc - 1)
        do na = 1, nat
           do nb = 1, nat
              write (iudyn, '(2i3)') na, nb
              do icar = 1, 3
                 write (iudyn, '(3e24.12)') (phi (kcar, icar, jcar, nc, na, nb) &
                      , jcar = 1, 3)
              enddo
           enddo
        enddo
     endif

  enddo

  return
9000 format(/,5x,'Third derivative in cartesian axes', &
       &       //,5x,'q = ( ',3f14.9,' ) ',/)
end subroutine write_d3dyn
