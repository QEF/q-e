!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE h_psiq_vdw (lda, n, m, psi, hpsi, spsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with a m  wavefunctions  contained
  !     in psi. It first computes the bec matrix of these
  !     wavefunctions and then with the routines hus_1psi and
  !     s_psi computes for each band the required products
  !
  USE fft_base, ONLY : dffts
  USE fft_interfaces, ONLY : fwfft, invfft
  USE pwcom
  USE scf, ONLY : vrs
  USE wavefunctions_module,  ONLY: psic
  USE becmod, ONLY: becp
  USE kinds, ONLY : DP
  USE phcom
  IMPLICIT NONE
  !
  !     Here the local variables
  !
  INTEGER :: ibnd
  ! counter on bands

  INTEGER :: lda, n, m
  ! input: the leading dimension of the array psi
  ! input: the real dimension of psi
  ! input: the number of psi to compute
  INTEGER :: j
  ! do loop index

  COMPLEX(kind=DP) :: psi (lda, m), hpsi (lda, m), spsi (lda, m)
  ! input: the functions where to apply H and S
  ! output: H times psi
  ! output: S times psi (Us PP's only)

  CALL start_clock ('h_psiq')
  CALL start_clock ('init')
!  call calbec ( n, vkb, psi, becp, m)   ! no need in TFvW
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  DO ibnd = 1, m
     DO j = 1, n
        hpsi (j, ibnd) = g2kin (j) * psi (j, ibnd)
     ENDDO
  ENDDO
  CALL stop_clock ('init')
  !
  ! the local potential V_Loc psi. First the psi in real space
  !

  DO ibnd = 1, m
     CALL start_clock ('firstfft')
     psic(:) = (0.d0, 0.d0)
     DO j = 1, n
        psic (nls(igkq(j))) = psi (j, ibnd)
     ENDDO
     CALL invfft ('Wave', psic, dffts)
     CALL stop_clock ('firstfft')
     !
     !   and then the product with the potential vrs = (vltot+vr) on the smoo
     !
     DO j = 1, dffts%nnr
        psic (j) = psic (j) * vrs (j, current_spin)
     ENDDO
     !
     !   back to reciprocal space
     !
     CALL start_clock ('secondfft')
     CALL fwfft ('Wave', psic, dffts)
     !
     !   addition to the total product
     !
     DO j = 1, n
        hpsi (j, ibnd) = hpsi (j, ibnd) + psic (nls(igkq(j)))
     ENDDO
     CALL stop_clock ('secondfft')
  ENDDO
  !
  !  Here the product with the non local potential V_NL psi
  !

!  call add_vuspsi (lda, n, m, hpsi)    ! no need in TFvW

  CALL s_psi (lda, n, m, psi, spsi)     ! no need in TFvW

  CALL stop_clock ('h_psiq')
  RETURN
END SUBROUTINE h_psiq_vdw
