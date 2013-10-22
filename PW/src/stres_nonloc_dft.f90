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
  use funct,            ONLY : gcxc, gcx_spin, gcc_spin, gcc_spin_more, &
                               dft_is_gradient, get_igcc, get_inlc 
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE vdW_DF,           ONLY : stress_vdW_DF, print_sigma 
  !
  IMPLICIT NONE
  !
  real(DP), intent(in)     :: rho (dfftp%nnr, nspin), rho_core (dfftp%nnr)
  real(DP), intent(inout)  :: sigma_nonloc_dft (3, 3)
  integer ::nspin, inlc

  integer :: l, m


  sigma_nonloc_dft(:,:) = 0.d0
  inlc = get_inlc()

  if (inlc==1 .or. inlc==2) then
     if (nspin>1) call errore('stres_vdW_DF', &
                  'vdW+DF non implemented in spin polarized calculations',1)
     CALL stress_vdW_DF(rho, rho_core, sigma_nonloc_dft)

  end if

  return

end subroutine stres_nonloc_dft

