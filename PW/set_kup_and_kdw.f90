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
  !                those in the second (nks/2) ones correspond to down spin
  !
  USE kinds, ONLY : DP
#if defined (EXX)
  USE exx, ONLY : exx_grid_check, lexx, xkq, index_xkq, index_xk, index_sym, &
                  nkqs, nqs
  USE klist, ONLY: nks_ => nks
#endif
  implicit none
  !
  ! I/O variables first
  !
  integer :: npk, isk (npk), nks
  ! input: maximum allowed number of k-points
  ! output: spin associated to a given k-point
  ! input-output: starting and ending number of k-points 
  real(kind=DP) :: xk (3, npk), wk (npk)
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  !
  integer :: ik, iq, ikq
  !
  !
  if (2*nks > npk) call errore ('set_kup&kdw','too many k points',nks)
  do ik = 1, nks
     xk(:,ik+nks)= xk(:,ik)
     wk (ik+nks) = wk(ik)
     isk(ik)     = 1
     isk(ik+nks) = 2
  enddo
#if defined (EXX)
  if (lexx) then
     if (nks /= nks_) call errore ('set_kup_and_kdw', 'wrong nks', 1)
     do ik =1, nks
        do iq =1, nqs
           index_xkq(ik + nks,iq) = index_xkq(ik,iq) + nkqs
        end do
     end do
     do ikq=1,nkqs
        xkq(:,ikq + nkqs)     = xkq(:,ikq)
        index_xk(ikq + nkqs)  = index_xk(ikq) + nks
        index_sym(ikq + nkqs) = index_sym(ikq)
     end do
  
     nkqs = 2 * nkqs
  end if
#endif
  nks = 2 * nks

#if defined (EXX)
  if (lexx) call exx_grid_check
#endif

  return

end subroutine set_kup_and_kdw
