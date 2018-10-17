!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE hs_1psi_gpu( lda, n, psi_d, hpsi_d, spsi_d )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the Hamiltonian and the S matrix
  ! ... to a vector psi and puts the result in hpsi and spsi
  ! ... Wrapper routine - calls h_psi and s_psi
  !
  ! ... No bgrp parallelization here !
  !
  USE kinds,  ONLY: DP
  USE control_flags, ONLY : gamma_only
  USE bp,     ONLY: lelfield
  USE noncollin_module, &
              ONLY: npol 
  USE realus, ONLY : real_space, &
                     invfft_orbital_gamma, fwfft_orbital_gamma, s_psir_gamma, &
                     invfft_orbital_k, fwfft_orbital_k, s_psir_k
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n
  COMPLEX (DP) :: psi_d(lda*npol,1), hpsi_d(n), spsi_d(n,1)
#if defined (__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, spsi_d
#endif
  COMPLEX (DP), ALLOCATABLE :: psi_h(:,:), spsi_h(:,:)
  !
  !
  CALL start_clock( 'hs_1psi' )
  ! 
  !OBM: I know this form is somewhat inelegant but, leaving the pre-real_space part intact
  !     makes it easier to debug probable errors, please do not "beautify" 
        if (real_space) then
           ALLOCATE(psi_h(lda*npol,1), spsi_h(n,1))
           psi_h(1:lda*npol,1) = psi_d(1:lda*npol,1)
           CALL h_psi_gpu( lda, n, 1, psi_d, hpsi_d )
           if (gamma_only) then
             call invfft_orbital_gamma(psi_h,1,1) !transform the orbital to real space
             call s_psir_gamma(1,1)
             call fwfft_orbital_gamma(spsi_h,1,1)
           else
             call invfft_orbital_k(psi_h,1,1) !transform the orbital to real space
             call s_psir_k(1,1)
             call fwfft_orbital_k(spsi_h,1,1)
           end if
           spsi_d(1:n,1) = spsi_h(1:n,1)
           DEALLOCATE(psi_h, spsi_h)
        else   
  CALL h_psi_gpu( lda, n, 1, psi_d, hpsi_d ) ! apply H to a single wfc (no bgrp parallelization here)
  CALL s_psi_gpu( lda, n, 1, psi_d, spsi_d ) ! apply S to a single wfc (no bgrp parallelization here)
       endif
  !
  CALL stop_clock( 'hs_1psi' )
  !
  RETURN
  !
END SUBROUTINE hs_1psi_gpu
