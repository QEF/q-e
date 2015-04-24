!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE s_1psi( npwx, n, psi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... spsi = S*psi for one wavefunction
  ! ... Wrapper routine - calls calbec and s_psi
  !
  USE kinds,  ONLY : DP
  USE uspp,   ONLY : vkb, nkb
  USE becmod, ONLY : bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only 
  USE noncollin_module, ONLY : noncolin, npol 
  USE realus,         ONLY : real_space, invfft_orbital_gamma, fwfft_orbital_gamma, &
                             calbec_rs_gamma, s_psir_gamma, initialisation_level
  USE wvfct,                ONLY: nbnd
  !
  IMPLICIT NONE
  !
  INTEGER          :: npwx, n, ibnd
  COMPLEX(DP) :: psi(npwx*npol,1), spsi(npwx*npol,1)
  !
  !
  CALL start_clock( 's_1psi' )
  !
  IF ( gamma_only  .and.  real_space) then
     do ibnd=1,nbnd,2
        ! transform the orbital to real space
        call invfft_orbital_gamma(psi,ibnd,nbnd) 
        ! global becp%r is updated
        call calbec_rs_gamma(ibnd,nbnd,becp%r) 
     enddo
     call s_psir_gamma(1,1)
     call fwfft_orbital_gamma(spsi,1,1)
     !
  ELSE
     !
     CALL calbec( n, vkb, psi, becp )
     !
  END IF
  !
  if (.not. real_space) CALL s_psi( npwx, n, 1, psi, spsi )
  !
  CALL stop_clock( 's_1psi' )
  !
  RETURN
  !
END SUBROUTINE s_1psi
