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
  USE noncollin_module, ONLY : noncolin, lspinorb, domag, npol, &
                               angle1, angle2, starting_spin_angle
  USE atomic_wfc_mod,   ONLY : atomic_wfc_acc
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
  LOGICAL :: starting_spin_angle_ , updown, spinorbit
  !
  CALL start_clock( 'atomic_wfc' )
  !
  starting_spin_angle_ = starting_spin_angle .OR. .NOT.domag
  updown = .FALSE.
  spinorbit = .FALSE.
  !$acc data present_or_copyout(wfcatom)
  CALL atomic_wfc_acc( xk(1,ik), ngk(ik), igk_k(1,ik), nat, nsp, ityp, tau, &
       noncolin, spinorbit, updown, angle1, angle2, starting_spin_angle_, &
       npwx, npol, natomwfc, wfcatom )
  !$acc end data
  !
  CALL stop_clock( 'atomic_wfc' )
  !
END SUBROUTINE atomic_wfc
!
!
!--------------------------------------------------------------------------
SUBROUTINE atomic_wfc_nc_updown( ik, wfcatom )
  !-----------------------------------------------------------------------
  !! Wrapper routine like "atomic_wfc", for the noncolinear case,
  !! when pure spin-up or spin-down atomic wavefunctions are required
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, tau, nsp, ityp
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : xk, ngk, igk_k
  USE wvfct,            ONLY : npwx
  USE noncollin_module, ONLY : noncolin, lspinorb, npol, angle1, angle2
  !! Note: lspinorb, angle1 and angle2 are not used
  USE atomic_wfc_mod,   ONLY : atomic_wfc_acc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx,npol,natomwfc)
  !! the superposition of atomic wavefunctions (up or down)
  !
  ! ... local variables
  !
  LOGICAL :: starting_spin_angle_ , updown, spinorbit
  !
  starting_spin_angle_ = .FALSE.
  updown = .TRUE.
  spinorbit = .FALSE.
  !! FIXME: the treatment of the spinorbit case without spinorbit potentials
  !! FIXME: is the only difference between this routine and atomic_wfc_nc_proj
  !! FIXME: Not sure whether this difference makes sense, and if so. why
  CALL start_clock( 'atomic_wfc' )
  !$acc data present_or_copyout(wfcatom)
  CALL atomic_wfc_acc( xk(1,ik), ngk(ik), igk_k(1,ik), nat, nsp, ityp, tau, &
       noncolin, spinorbit, updown, angle1, angle2, starting_spin_angle_, &
       npwx, npol, natomwfc, wfcatom )
  !$acc end data
  CALL stop_clock( 'atomic_wfc' )
  !
END SUBROUTINE atomic_wfc_nc_updown  
!
!--------------------------------------------------------------------------
SUBROUTINE atomic_wfc_nc_proj( ik, wfcatom )
  !-----------------------------------------------------------------------
  !! Wrapper routine like "atomic_wfc", for the noncolinear case.
  !! If lspinorb=.TRUE. it makes linear combinations of eigenstates of
  !! the atomic total angular momenta j and j_z; otherwise, of eigenstates of
  !! the orbital angular momenta l, l_z and of s_z (the z-component of spin).
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, tau, nsp, ityp
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : xk, ngk, igk_k
  USE wvfct,            ONLY : npwx
  USE noncollin_module, ONLY : noncolin, lspinorb, npol, angle1, angle2
  !! Note: angle1 and angle2 are not used
  USE atomic_wfc_mod,   ONLY : atomic_wfc_acc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx,npol,natomwfc)
  !! the superposition of atomic wavefunctions (up or down)
  !
  ! ... local variables
  !
  LOGICAL :: starting_spin_angle_ , updown, spinorbit
  !
  starting_spin_angle_ = .TRUE.
  updown = .TRUE.
  spinorbit = lspinorb
  !
  CALL start_clock( 'atomic_wfc' )
  !$acc data present_or_copyout(wfcatom)
  CALL atomic_wfc_acc( xk(1,ik), ngk(ik), igk_k(1,ik), nat, nsp, ityp, tau, &
       noncolin, spinorbit, updown, angle1, angle2, starting_spin_angle_, &
       npwx, npol, natomwfc, wfcatom )
  !$acc end data
  CALL stop_clock( 'atomic_wfc' )
  !
END SUBROUTINE atomic_wfc_nc_proj
