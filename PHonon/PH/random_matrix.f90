!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine random_matrix_new (irt, nsymq, minus_q, irotmq, nat, &
     wdyn, lgamma)
  !----------------------------------------------------------------------
  !
  !   Create a random hermitian matrix with non zero elements similar to
  !   the dynamical matrix of the system
  !
  !
  USE kinds, only : DP
  USE random_numbers, ONLY : randy
  implicit none
  !
  !    The dummy variables
  !

  integer :: nat, irt (48, nat), nsymq, irotmq
  ! input: number of atoms
  ! input: index of the rotated atom
  ! input: the small group of q
  ! input: the order of the small group
  ! input: the rotation sending q -> -q

  complex(DP) :: wdyn (3, 3, nat, nat)
  ! output: random matrix
  logical :: lgamma, minus_q
  ! input: if true q=0
  ! input: if true there is a symmetry
  !
  !    The local variables
  !
  integer :: na, nb, ipol, jpol, isymq, irot, ira, iramq
  ! counters
  ! ira:   rotated atom
  ! iramq: rotated atom with the q->-q+G symmetry
  !
  !
  wdyn (:, :, :, :) = (0d0, 0d0)
  do na = 1, nat
     do ipol = 1, 3
        wdyn (ipol, ipol, na, na) = CMPLX(2.0_DP * randy () - 1.0_DP, 0.d0,kind=DP)
        do jpol = ipol + 1, 3
           if (lgamma) then
              wdyn (ipol, jpol, na, na) = CMPLX(2.0_DP * randy () - 1.0_DP, 0.d0,kind=DP)
           else
              wdyn (ipol, jpol, na, na) = &
                   CMPLX(2.0_DP * randy () - 1.0_DP, 2.0_DP * randy () - 1.0_DP,kind=DP)
           endif
           wdyn (jpol, ipol, na, na) = CONJG(wdyn (ipol, jpol, na, na) )
        enddo
        do nb = na + 1, nat
           do isymq = 1, nsymq
              irot = isymq
              ira = irt (irot, na)
              if (minus_q) then
                 iramq = irt (irotmq, na)
              else
                 iramq = 0
              endif
              if ( (nb == ira) .or. (nb == iramq) ) then
                 do jpol = 1, 3
                    if (lgamma) then
                       wdyn (ipol, jpol, na, nb) = CMPLX(2.0_DP*randy () - 1.0_DP, 0.d0,kind=DP)
                    else
                       wdyn (ipol, jpol, na, nb) = &
                            CMPLX(2.0_DP*randy()-1.0_DP, 2.0_DP*randy()-1.0_DP,kind=DP)
                    endif
                    wdyn(jpol, ipol, nb, na) = CONJG(wdyn(ipol, jpol, na, nb))
                 enddo
                 goto 10
              endif
           enddo
10         continue
        enddo
     enddo
  enddo
  return
end subroutine random_matrix_new
