!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine write_igk
  !
  USE klist,      ONLY: ngk, igk_k
  USE wvfct,      ONLY: npwx
  USE io_files,   ONLY: iunigk
  USE qpoint,     ONLY: nksq, ikks, ikqs
  USE control_lr, ONLY: lgamma

  implicit none
  integer :: ik, ikk, ikq
  
  rewind (unit = iunigk)
  do ik =1,nksq
     ikk=ikks(ik)
     write (iunigk) ngk(ikk), igk_k(:,ikk)
     if (.not.lgamma) then
        ikq=ikqs(ik)
        write (iunigk) ngk(ikq), igk_k(:,ikq)
     end if
  end do
  return
end subroutine write_igk
