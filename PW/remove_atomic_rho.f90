!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine remove_atomic_rho
  !-----------------------------------------------------------------------
#include "machine.h"
  USE io_global, ONLY: stdout
  USE io_files, ONLY: output_drho
  USE kinds, ONLY: DP
  USE gvect, ONLY: nrxx
  USE lsda_mod, ONLY: lsda, nspin
  USE scf, ONLY: rho
  implicit none
  integer :: ir
  ! do-loop variable on FFT grid

  real(kind=DP), allocatable :: work (:)
  real(kind=DP)          :: charge
  ! workspace, is the difference between t
  ! charge density and the atomic one at t
  ! charge

  allocate (work( nrxx))    
  work(:) = 0.d0
  !
  if (lsda) call errore ('rmv_at_rho', 'lsda not allowed', 1)

  WRITE( stdout, '(/5x,"remove atomic charge density from scf rho")')
  !
  !     subtract the old atomic charge density
  !
  call atomic_rho (work, nspin)
  call DSCAL (nrxx, - 1.0d0, work, 1)
  call DAXPY (nrxx, + 1.0d0, rho, 1, work, 1)

  call io_pot ( + 1, output_drho, work, nspin)

  deallocate(work)
  return

end subroutine remove_atomic_rho

