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
  USE realus,         ONLY : real_space, fft_orbital_gamma, bfft_orbital_gamma, &
                             calbec_rs_gamma, s_psir_gamma, initialisation_level,check_fft_orbital_gamma
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
  IF ( gamma_only ) THEN
     ! 
     if (real_space) then
            !if (.not. initialisation_level == 15) CALL errore ('s_1psi', 'improper initialisation of real space routines' , 4)
            !print *, "s1_psi rolling the real space!" 
            !do ibnd = 1 , nbnd , 2
             !call check_fft_orbital_gamma(psi,1,1)
             do ibnd=1,nbnd,2
              call fft_orbital_gamma(psi,ibnd,nbnd) !transform the orbital to real space
              call calbec_rs_gamma(ibnd,nbnd,becp%r) !global becp%r is updated
             enddo
             call s_psir_gamma(1,1)
             call bfft_orbital_gamma(spsi,1,1)
            !enddo
      !call calbec( n, vkb, psi, becp, -nbnd )
      !call s_psi( npwx, n, 1, psi, spsi )
     else
      CALL calbec( n, vkb, psi, becp )
     endif
     !
  ELSE IF ( noncolin ) THEN
     !
     CALL calbec( n, vkb, psi, becp )
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
