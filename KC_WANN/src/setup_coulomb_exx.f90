!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------
subroutine setup_coulomb_exx ( )
  !--------------------------------------------------------------------
  !
  !	
  !!  This routine initialize the coulomb kernel according to the 
  !!  Gygi-Balderschi scheme. This is compatible with peridic 
  !!  systems and finite sytems calculation. For finite system 
  !!  one can also use assume_isolated (see kc_readin) 
  !
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout
  USE funct,          ONLY : init_dft_exxrpa, dft_is_hybrid, &
                             exx_is_active, stop_exx
  USE exx_base,       ONLY : exx_grid_init, exx_div_check, &
                             exxdiv_treatment, x_gamma_extrapolation, &
                             exx_divergence, nq1, nq2, nq3, &
                             exx_grid_initialized,  exxdiv
  USE exx,            ONLY : exxinit, exxenergy, exxenergy2, deallocate_exx, ecutfock
  USE gvect,          ONLY : ecutrho
  USE control_kc_wann, ONLY : mp1, mp2, mp3, l_vcut
  !
  implicit none
  !
  !
  CALL start_clock( 'Coulomb setup' )
  !
  ecutfock = ecutrho
  WRITE (stdout,'(/, 5X, "INFO: Coulomb cut-off  ", f12.8)') ecutfock
  !
  call deallocate_exx()
  ! 
  nq1=mp1; nq2=mp2; nq3=mp3
  !
  x_gamma_extrapolation = .false.
! DEBUG >>>
!  x_gamma_extrapolation = .true. 
! DEBUG <<<
  ! The Screened contribution does not work with x_gamma_extrapolation
  !
  exxdiv_treatment = 'none'
  IF (l_vcut) exxdiv_treatment = 'gb'
  !
  exx_grid_initialized = .FALSE.  
  CALL exx_grid_init()
  CALL exx_div_check()
  exxdiv = exx_divergence ()
  ! Compute the divergence and set the q+G=0 term (exxdiv, to be used in
  ! g2_convolution called by bare_pot.f90
  !
  WRITE (stdout,'(/,5X, "INFO: Divergence         ", 3x, 1A8)')       exxdiv_treatment
  WRITE (stdout,'(  5X, "INFO: Gamma Extrapolation", 3x, 1L5 )')      x_gamma_extrapolation
  WRITE (stdout,'(  5X, "INFO: Coulomb G0         ", 3x, 1ES15.5 )')  exx_divergence ()
  !
  CALL stop_clock( 'Coulomb setup' )
  !
  return
  !
end subroutine setup_coulomb_exx
