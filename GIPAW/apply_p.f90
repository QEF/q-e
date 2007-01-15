! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE apply_p(psi, p_psi, ik, ipol, q)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the kinetic part of the velocity operator
  ! ... |p_psi> = (G+k+q/2)_{ipol} |psi>
  !  
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk  
  USE pwcom
  USE nmr_module

  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik               ! k-point
  INTEGER, INTENT(IN) :: ipol             ! cartesian direction (1..3)
  REAL(DP), INTENT(IN) :: q(3)
  COMPLEX(DP), INTENT(IN) :: psi(npwx,nbnd)
  COMPLEX(DP), INTENT(OUT) :: p_psi(npwx,nbnd)

  !-- local variables ----------------------------------------------------
  REAL(DP) :: gk
  INTEGER :: ig, ibnd

  call start_clock('apply_p')

  do ibnd = 1, nbnd_occ(ik)
    do ig = 1, npw
      gk = xk(ipol,ik) + g(ipol,igk(ig)) + q(ipol)
      p_psi(ig,ibnd) = gk * tpiba * psi(ig,ibnd)
    enddo
  enddo

  call stop_clock('apply_p')

END SUBROUTINE apply_p

 
