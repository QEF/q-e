!
! Copyright (C) 2023-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc( ik, wfcatom )
  !-----------------------------------------------------------------------
  !! Wrapper routine: calls atomic_wfc_acc to compute the (non-orthogonal)
  !! superposition of atomic wavefunctions at the ik-th k-point - output in 
  !! "wfcatom", on GPU if input is an ACC variable, copied to CPU otherwise
  !!
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, tau, nsp, ityp
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : xk, ngk, igk_k
  USE wvfct,            ONLY : npwx
  USE noncollin_module, ONLY : noncolin, domag, npol, angle1, angle2, &
                               starting_spin_angle
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx,npol,natomwfc)
  !! Superposition of atomic wavefunctions
  !
  ! ... local variables
  !
  CALL start_clock( 'atomic_wfc' )
  !
  !$acc data present_or_copyout(wfcatom)
  CALL atomic_wfc_acc( xk(1,ik), ngk(ik), igk_k(1,ik), nat, nsp, ityp, tau, &
       noncolin, domag, angle1, angle2, starting_spin_angle, &
       npwx, npol, natomwfc, wfcatom )
  !$acc end data
  !
  CALL stop_clock( 'atomic_wfc' )
  !
END SUBROUTINE atomic_wfc

