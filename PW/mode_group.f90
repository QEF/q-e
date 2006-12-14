!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------

subroutine mode_group (modenum, xq, at, bg, nat, nrot, s, irt, &
     rtau, sym, minus_q)
  !-----------------------------------------------------------------------
  !
  ! This routine selects, among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave a given mode unchang
  ! For the moment it assume that the mode modenum displaces the atom
  ! modenum/3 in the direction mod(modenum,3)+1
  ! Also the minus_q operation is tested.
  !
  !  input-output variables
  !
  USE kinds
  USE constants, ONLY : tpi
  implicit none

  integer :: nat, s (3, 3, 48), irt (48, nat), nrot, modenum
  ! input: the number of atoms of the system
  ! input: the symmetry matrices
  ! input: the rotated atom
  ! input: number of symmetry operations
  ! input: the displacement pattern


  real(DP) :: xq (3), rtau (3, 48, nat), bg (3, 3), at (3, 3)
  ! input: the q point
  ! input: the translations of each atom
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  logical :: minus_q, sym (48)
  ! input: if true minus_q symmetry is used
  ! input-output: .true. if symm. op. do not change
  ! mode
  !
  !  local variables
  !

  integer :: isym, nas, ipols, na, sna, ipol, jpol
  ! counters
  ! counter on polarizations
  ! counter on polarizations

  real(DP) :: arg
  ! auxiliary

  complex(DP), allocatable :: u (:,:)
  ! the original pattern
  complex(DP)              :: fase, sum
  ! the phase of the mode
  ! check for orthogonality
  complex(DP), allocatable :: work_u (:,:), work_ru (:,:)
  ! the working pattern
  ! the rotated working pattern


  allocate(u(3, nat), work_u(3, nat), work_ru (3, nat))

  if (modenum.gt.3 * nat.or.modenum.lt.1) call errore ('mode_group', &
       'wrong modenum', 1)
  nas = (modenum - 1) / 3 + 1
  ipols = mod (modenum - 1, 3) + 1
  u (:,:) = (0.d0, 0.d0)
  u (ipols, nas) = (1.d0, 0.d0)
  do na = 1, nat
     call trnvecc (u (1, na), at, bg, - 1)
  enddo
  do isym = 1, nrot
     if (sym (isym) ) then
        do na = 1, nat
           do ipol = 1, 3
              work_u (ipol, na) = u (ipol, na)
           enddo
        enddo
        work_ru (:,:) = (0.d0, 0.d0)
        do na = 1, nat
           sna = irt (isym, na)
           arg = 0.d0
           do ipol = 1, 3
              arg = arg + xq (ipol) * rtau (ipol, isym, na)
           enddo
           arg = arg * tpi
           if (isym.eq.nrot.and.minus_q) then
              fase = CMPLX (cos (arg), sin (arg) )
           else
              fase = CMPLX (cos (arg), - sin (arg) )
           endif
           do ipol = 1, 3
              do jpol = 1, 3
                 work_ru (ipol, sna) = work_ru (ipol, sna) + s (jpol, ipol, &
                      isym) * work_u (jpol, na) * fase
              enddo
           enddo
        enddo
        !
        !    Transform back the rotated pattern
        !
        do na = 1, nat
           call trnvecc (work_ru (1, na), at, bg, 1)
           call trnvecc (work_u (1, na), at, bg, 1)
        enddo
        !
        !   only if the pattern remain the same ap to a phase we keep
        !   the symmetry
        !
        sum = (0.d0, 0.d0)
        do na = 1, nat
           do ipol = 1, 3
              sum = sum + CONJG(work_u (ipol, na) ) * work_ru (ipol, na)
           enddo
        enddo
        sum = abs (sum)
        if (abs (sum - 1.d0) .gt.1.d-7) sym (isym) = .false.
     endif

  enddo
  deallocate ( work_ru, work_u, u)
  return

end subroutine mode_group

