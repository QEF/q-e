!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine mode_group &
    (modenum, xq, at, bg, nat, nrot, s, irt, minus_q, rtau, sym)
  !-----------------------------------------------------------------------
  !
  ! This routine selects, among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave a given mode unchanged
  ! For the moment it assumes that the mode modenum displaces the atom
  ! modenum/3 in the direction mod(modenum,3)+1
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : tpi
  implicit none

  integer, intent(in) :: nat, s (3, 3, 48), irt (48, nat), nrot, modenum
  ! nat : the number of atoms of the system
  ! s   : the symmetry matrices
  ! irt : the rotated atom
  ! nrot: number of symmetry operations
  ! modenum: which displacement pattern

  real(DP), intent(in) :: xq (3), rtau (3, 48, nat), bg (3, 3), at (3, 3)
  ! xq  : the q point
  ! rtau: the translations of each atom
  ! bg  : the reciprocal lattice vectors
  ! at  : the direct lattice vectors
  logical, intent(in) :: minus_q
  ! if true Sq=>-q+G  symmetry is used
  logical, intent(inout) :: sym (48)
  ! on  input: .true. if symm. op. has to be tested
  ! on output: .true. if symm. op. does not change mode modenum
  !
  integer :: isym, nas, ipols, na, sna, ipol, jpol
  ! counters
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

  if (modenum > 3*nat .or. modenum < 1) call errore ('mode_group', &
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
              fase = CMPLX(cos (arg), sin (arg) ,kind=DP)
           else
              fase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
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
        !   only if the pattern remain the same up to a phase we keep
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

