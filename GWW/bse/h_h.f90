!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE H_h(e,h,Ah)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY: DP
  USE wvfct, ONLY: npwx, npw, g2kin
  USE gvect, ONLY : gstart
  USE uspp,     ONLY : vkb, nkb
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs
  USE becmod, ONLY: bec_type, becp, calbec
!  USE cgcom
  USE electrons_base, ONLY: nel
  use bse_wannier, ONLY:num_nbndv
  !
  IMPLICIT NONE
  !
  real(DP):: e(num_nbndv(1))
  COMPLEX(DP):: h(npwx,num_nbndv(1)), Ah(npwx,num_nbndv(1))
!  real(DP), allocatable     :: e(:)
!  COMPLEX(DP), allocatable  :: h(:,:), Ah(:,:)
  !
  INTEGER:: j,ibnd,nbnd
  !
  CALL start_clock('h_h')

  ! valid only for non-spin resolved calculations
  nbnd=num_nbndv(1)

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
   CALL calbec  ( npw, vkb, h, becp, nbnd )
  IF (nkb > 0) CALL add_vuspsi (npwx, npw, nbnd, ah)
  ! set to zero the imaginary part of ah at G=0
  !  needed for numerical stability
  IF (gstart==2) THEN
     DO ibnd = 1, nbnd
        ah(1,ibnd) = cmplx( dble(ah(1,ibnd)),0.d0,kind=DP)
     ENDDO
  ENDIF

!  DEALLOCATE(h)
!  DEALLOCATE(Ah)
!  DEALLOCATE(e)
  !
  CALL stop_clock('h_h')
  !
  RETURN
END SUBROUTINE H_h
