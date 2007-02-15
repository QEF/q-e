!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine set_kup_and_kdw (xk, wk, isk, nkstot, npk)
  !-----------------------------------------------------------------------
  !     This routine sets the k vectors for the up and down spin wfc
  !
  !     on input: xk and wk contain k-points and corresponding weights
  !
  !     on output: the number of points is doubled and xk and wk in the
  !                first (nkstot/2) positions correspond to up spin
  !                those in the second (nkstot/2) ones correspond to down spin
  !
  USE kinds, ONLY : DP
#if defined (EXX)
  USE exx, ONLY : exx_grid_check, xkq, index_xkq, index_xk, index_sym, nkqs, nqs
  USE funct, ONLY: dft_is_hybrid
  USE klist, ONLY: nkstot_ => nkstot
#endif
  implicit none
  !
  ! I/O variables first
  !
  integer :: npk, isk (npk), nkstot
  ! input: maximum allowed number of k-points
  ! output: spin associated to a given k-point
  ! input-output: starting and ending number of k-points 
  real(DP) :: xk (3, npk), wk (npk)
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  !
  integer :: ik, iq, ikq
  !
  !
  if (2*nkstot > npk) call errore ('set_kup&kdw','too many k points',nkstot)
  do ik = 1, nkstot
     xk(:,ik+nkstot)= xk(:,ik)
     wk (ik+nkstot) = wk(ik)
     isk(ik)     = 1
     isk(ik+nkstot) = 2
  enddo
#if defined (EXX)
  if (dft_is_hybrid()) then
     if (nkstot /= nkstot_) call errore ('set_kup_and_kdw', 'wrong nkstot', 1)
     do ik =1, nkstot
        do iq =1, nqs
           index_xkq(ik + nkstot,iq) = index_xkq(ik,iq) + nkqs
        end do
     end do
     do ikq=1,nkqs
        xkq(:,ikq + nkqs)     = xkq(:,ikq)
        index_xk(ikq + nkqs)  = index_xk(ikq) + nkstot
        index_sym(ikq + nkqs) = index_sym(ikq)
     end do
  
     nkqs = 2 * nkqs
  end if
#endif
  nkstot = 2 * nkstot

#if defined (EXX)
  if (dft_is_hybrid()) call exx_grid_check
#endif

  return

end subroutine set_kup_and_kdw
