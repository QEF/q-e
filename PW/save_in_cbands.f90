!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine save_in_cbands (iter, ik_, dr2)  
  !-----------------------------------------------------------------------
  use pwcom  
  implicit none
  character :: where * 20  
  ! are we in the right place?
  integer :: ik, ibnd, ik_, iter  
  ! counters
  ! last completed kpoint
  ! last completed iteration
  logical :: exst  

  real(kind=DP) :: dr2  
  if (reduce_io) return  
  !
  ! open recover file
  !
  call seqopn (iunres, 'restart', 'unformatted', exst)  
  !
  ! save restart information
  !
  where = 'ELECTRONS'  
  write (iunres) where  
  write (iunres) ( (et (ibnd, ik), ibnd = 1, nbnd), ik = 1, nks)  

  write (iunres) iter, ik_, dr2, tr2, ethr  

  close (unit = iunres, status = 'keep')  
  !
  return  
end subroutine save_in_cbands
