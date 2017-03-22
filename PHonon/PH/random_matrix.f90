!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! The random matrix (see later) can be populated either with uniformly distributed
! random numbers or with normal-distributed randm numbers. The former has been the default
! until QE 6.1, however it sometimes produces accidentally degenerate eigenvalue, especially
! when dealing with large number of atoms.
! A matrix of normal-distributed numbers should have (on average) more evenly spaced
! eigenvalues, reducing the chance of collision.
!
! See <http://web.math.princeton.edu/mathlab/projects/ranmatrices/yl/randmtx.PDF>
! (If I understand it correctly)
! LP 2017
!
!!#define __UNIFORM_DISTRIB
!
!----------------------------------------------------------------------
subroutine random_matrix_new (irt, nsymq, minus_q, irotmq, nat, &
     wdyn, lgamma)
  !----------------------------------------------------------------------
  !
  !   Create a random hermitian matrix with non zero elements similar to
  !   the dynamical matrix of the system
  !
#if defined (__UNIFORM_DISTRIB)
#define __RANDOM_DBLE  CMPLX(2.0_DP*randy () - 1.0_DP, 0.d0,kind=DP)
#define __RANDOM_CMPLX CMPLX(2.0_DP * randy () - 1.0_DP, 2.0_DP * randy () - 1.0_DP,kind=DP)
#else
#define __RANDOM_DBLE  gauss_dist_scal(0._dp, 1._dp)
#define __RANDOM_CMPLX gauss_dist_cmplx(0._dp, 1._dp)
  USE random_numbers, ONLY : gauss_dist_scal, gauss_dist_cmplx
#endif
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
        wdyn (ipol, ipol, na, na) = 2*__RANDOM_DBLE
        do jpol = ipol + 1, 3
           if (lgamma) then
              wdyn (ipol, jpol, na, na) = __RANDOM_DBLE
           else
              wdyn (ipol, jpol, na, na) = __RANDOM_CMPLX
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
                       wdyn (ipol, jpol, na, nb) = __RANDOM_DBLE
                    else
                       wdyn (ipol, jpol, na, nb) = __RANDOM_CMPLX
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
