!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine ortho
  !
#include "machine.h"
  !
  ! orthonormalize the old wavefunctions with the new S matrix, it must be
  ! after any move of atoms and only if the US problem is solved with
  ! overlap=.false.
  !
  use pwcom
  USE io_files, ONLY: iunwfc, nwordwfc, iunigk
  USE wavefunctions_module,    ONLY : evc
  use becmod
  implicit none
  integer :: ik

  complex(kind=DP), allocatable :: sevc (:,:), dummy (:,:)

  allocate (sevc ( npwx , nbnd))    
  allocate (dummy( npwx , nbnd))    
  sevc(:,:) = (0.d0,0.d0)
  dummy(:,:) = (0.d0,0.d0)
  if (nks.gt.1) rewind (iunigk)
  do ik = 1, nks
     if (nks.gt.1) read (iunigk) npw, igk
     !
     ! read the wavefunctions
     !
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     !
     ! calculate becp
     !
     call init_us_2 (npw, igk, xk (1, ik), vkb)

     call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     !
     !  find S|psi> with the new "moved" S
     !
     call s_psi (npwx, npw, nbnd, evc, sevc)
     !
     !  orthonormalize
     !
     call cgramg1 (npwx, nbndx, npw, 1, nbnd, evc, sevc, dummy)
     !
     !  now <psi_i|S|psi_j> = delta_ij; write the wavefunctions on file
     !
     call davcio (evc, nwordwfc, iunwfc, ik, 1)

  enddo
  deallocate (sevc)
  deallocate (dummy)
  return
end subroutine ortho

