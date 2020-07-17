!
! Copyright (C) 2010- Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_nonloc_dft( rho, rho_core, nspin, sigma_nonloc_dft )

  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE funct,            ONLY : get_igcc, get_inlc 
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE vdW_DF,           ONLY : vdW_DF_stress
  USE rVV10,            ONLY : rVV10_stress
  USE io_global,        ONLY : stdout
  !
  IMPLICIT NONE
  !
  integer,  intent(in)     :: nspin
  real(DP), intent(in)     :: rho (dfftp%nnr), rho_core (dfftp%nnr)
  real(DP), intent(inout)  :: sigma_nonloc_dft (3, 3)

  integer :: l, m, inlc


  sigma_nonloc_dft(:,:) = 0.D0
  inlc = get_inlc()

  IF ( inlc > 0 .AND. inlc < 26 ) THEN
     CALL vdW_DF_stress (rho, rho_core, nspin, sigma_nonloc_dft)
  ELSEIF ( inlc == 26 ) THEN
     CALL rVV10_stress  (rho, rho_core, nspin, sigma_nonloc_dft)
  END IF

END SUBROUTINE stres_nonloc_dft
