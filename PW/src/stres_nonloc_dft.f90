!
! Copyright (C) 2010- Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
subroutine stres_nonloc_dft( rho, rho_core, nspin, sigma_nonloc_dft )

  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  use funct,            ONLY : get_igcc, get_inlc 
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE vdW_DF,           ONLY : stress_vdW_DF
  USE rVV10,            ONLY : stress_rVV10
  !
  IMPLICIT NONE
  !
  integer, intent(in) ::nspin
  real(DP), intent(in)     :: rho (dfftp%nnr, nspin), rho_core (dfftp%nnr)
  real(DP), intent(inout)  :: sigma_nonloc_dft (3, 3)

  integer :: l, m, inlc


  sigma_nonloc_dft(:,:) = 0.d0
  inlc = get_inlc()

  if ( inlc==1 .or. inlc==2 .or. inlc==4 .or. inlc==5 .or. inlc==6 ) then
     if ( nspin>1 ) call errore ('stres_vdW_DF',   'vdW stress not implemented for nspin > 1',1)
     CALL stress_vdW_DF(rho, rho_core, nspin, sigma_nonloc_dft)
  elseif ( inlc == 3 ) then
     if ( nspin>2 ) call errore ('stress_rVV10', 'rVV10 stress not implemented for nspin > 2',1)
     CALL stress_rVV10(rho, rho_core, nspin, sigma_nonloc_dft)
  end if

  return

end subroutine stres_nonloc_dft

