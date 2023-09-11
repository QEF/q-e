
!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc( ik, wfcatom )
   !-----------------------------------------------------------------------
   !! This routine computes the superposition of atomic wavefunctions
   !! for k-point "ik" - output in "wfcatom".
   !
   USE kinds,            ONLY : DP
   USE cell_base,        ONLY : omega, tpiba
   USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau
   USE basis,            ONLY : natomwfc
   USE gvect,            ONLY : mill, eigts1, eigts2, eigts3, g
   USE klist,            ONLY : xk, igk_k, ngk
   USE wvfct,            ONLY : npwx
   USE uspp_param,       ONLY : upf, nwfcm
   USE noncollin_module, ONLY : noncolin, domag, npol, angle1, angle2, &
                                starting_spin_angle
   USE upf_spinorb,      ONLY : rot_ylm
   USE atomic_wfc_init, only : atomic_wfc_nog
   !use atomic_wfc_noglob
   !
   implicit none
   INTEGER, INTENT(IN) :: ik
   !! k-point index
   COMPLEX(DP), INTENT(INOUT) :: wfcatom( npwx, npol, natomwfc )
   !! Superposition of atomic wavefunctions

   call atomic_wfc_nog(ik, wfcatom, omega, tpiba, nat, ntyp, ityp, tau, natomwfc, &
         mill, eigts1, eigts2, eigts3, g, xk, igk_k, ngk, npwx, &
         upf, nwfcm, noncolin, domag, npol, angle1, angle2, starting_spin_angle, &
         rot_ylm)

END SUBROUTINE atomic_wfc