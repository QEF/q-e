!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE H_h(npw,e,h,Ah)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY: DP
  USE wvfct, ONLY: nbnd, npwx, g2kin
  USE gvect, ONLY : gstart
  USE uspp,     ONLY : vkb, nkb
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs
  USE becmod, ONLY: bec_type, becp, calbec
  USE cgcom
  !
  IMPLICIT NONE
  !
  INTEGER :: npw
  real(DP):: e(nbnd)
  COMPLEX(DP):: h(npwx,nbnd), Ah(npwx,nbnd)
  !
  INTEGER:: j,ibnd
  !
  CALL start_clock('h_h')
  !
  ! [(k+G)^2 - e ]psi
  DO ibnd = 1,nbnd
     ! set to zero the imaginary part of h at G=0
     !  needed for numerical stability
     IF (gstart==2) h(1,ibnd) = cmplx( dble(h(1,ibnd)),0.d0,kind=DP)
     DO j = 1,npw
        ah(j,ibnd) = (g2kin(j)-e(ibnd)) * h(j,ibnd)
     ENDDO
  ENDDO
  ! V_Loc psi
  CALL vloc_psi_gamma(npwx, npw, nbnd, h, vrs(1,current_spin), ah)
  ! V_NL psi
   CALL calbec  ( npw, vkb, h, becp )
  IF (nkb > 0) CALL add_vuspsi (npwx, npw, nbnd, ah)
  ! set to zero the imaginary part of ah at G=0
  !  needed for numerical stability
  IF (gstart==2) THEN
     DO ibnd = 1, nbnd
        ah(1,ibnd) = cmplx( dble(ah(1,ibnd)),0.d0,kind=DP)
     ENDDO
  ENDIF
  !
  CALL stop_clock('h_h')
  !
  RETURN
END SUBROUTINE H_h
