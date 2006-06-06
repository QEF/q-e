!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine raman
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  use kinds, only : DP
  use pwcom
  use phcom
  USE ramanm, ONLY: lraman, elop
  implicit none

  if (okvan) &
      call errore ('raman','Ultrasoft pseudopotentials not implemented',1)
  if (lsda) call errore ('raman',' spin-polarized case not implemented',1)
  if (degauss.ne.0.d0 .or..not.lgamma) &
      call errore ('raman','called in the wrong case',1)
  if (epsil.and..not.convt) &
      call errore ('raman','epsil calcul. not converged',1)
  !
  ! Computes Pc [DH,Drho] |psi>
  !
  IF (irr0 == -10) THEN
     ! restart from a previous calculation
     write (6,'(/,5x,''Skipping computation of Pc [DH,Drho] |psi> '')')
  ELSE
     write (6,'(/,5x,''Computing Pc [DH,Drho] |psi> '')')
     call dhdrhopsi
  END IF
  !
  ! Computes the electro-optic tensor
  !
  if (elop) call el_opt
  if (.not.lraman) return
  write (6,'(/,5x,''Computing Second order response '')')
  !
  ! Computes the potential that remains unchanged in the scf-cycle
  !
  call dvpsi_e2
  !
  ! Self-consistent cycle to compute the second derivative of the charge
  !
  call solve_e2
  !
  ! Computes and writes the Raman tensor
  !
  call raman_mat

  return
end subroutine raman
