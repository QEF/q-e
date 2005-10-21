!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine inverse_s (nrot, s, table, invs)
  !---------------------------------------------------------------------
  implicit none
  integer :: nrot, s (3, 3, 48), table (48, 48), invs (3, 3, 48)
  ! input: number of symmetries of the original
  ! input: matrices of the symmetry operations
  ! input: multiplication table of the group
  ! output: contains the inverse of each rotati


  integer :: irot, jrot, ipol, jpol
  ! counter over the rotations
  ! counter over the rotations
  ! counter over the polarizations
  ! counter over the polarizations
  do irot = 1, nrot
     do jrot = 1, nrot
        if (table (irot, jrot) .eq.1) then
           do ipol = 1, 3
              do jpol = 1, 3
                 invs (ipol, jpol, irot) = s (ipol, jpol, jrot)
              enddo
           enddo
        endif
     enddo

  enddo
  return
end subroutine inverse_s

