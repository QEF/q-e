!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine set_kup_and_kdw (xk, wk, isk, nks, npk)
  !-----------------------------------------------------------------------
  !     This routine sets the k vectors for the up and down spin wfc
  !
  !     on input: xk and wk contain k-points and corresponding weights
  !
  !     on output: the number of points is doubled and xk and wk in the
  !                first (nks/2) positions correspond to up spin
  !                those in the second (nks/2) one correspond to down spin
  !                the weights wk are halved !
  !
  use parameters, only : DP
  implicit none
  !
  ! I/O variables first
  !

  integer :: npk, isk (npk), nks
  ! input-output: maximum allowed number of k
  ! output: spin associated to a given k-point
  ! input-output: starting and ending number of
  real(kind=DP) :: xk (3, npk), wk (npk)
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  !
  !    And then the local variables
  !

  integer :: ik, j
  ! counter on k
  ! counter

  if (2*nks.gt.npk) call errore ('set_kup&kdw','too many k points',nks)
  do ik = 1, nks
     do j = 1, 3
        xk(j,ik+nks) = xk(j,ik)
     enddo
     wk (ik)     = 0.5 * wk(ik)
     wk (ik+nks) =       wk(ik)
     isk(ik)     = 1
     isk(ik+nks) = 2
  enddo

  nks = 2 * nks
  return

end subroutine set_kup_and_kdw
