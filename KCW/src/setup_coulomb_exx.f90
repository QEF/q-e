!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
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
  !!  one can also use assume_isolated (see kcw_readin) 
  !
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE exx_base,             ONLY : exx_grid_init, exx_div_check, &
                                   exxdiv_treatment, x_gamma_extrapolation, &
                                   exx_divergence, nq1, nq2, nq3, &
                                   exx_grid_initialized,  exxdiv
  USE exx,                  ONLY : exxinit, deallocate_exx
  USE control_kcw,          ONLY : mp1, mp2, mp3, l_vcut, eps_inf
  USE command_line_options, ONLY : nband_, ntg_
  USE mp_exx,               ONLY : mp_start_exx
  USE mp_pools,             ONLY : intra_pool_comm
  !
  implicit none
  !
  !
  CALL start_clock( 'Coulomb setup' )
  !
  CALL deallocate_exx()
  ! 
  nq1=mp1; nq2=mp2; nq3=mp3
  !
  x_gamma_extrapolation = .false.
  ! NOTABENE: The Screened contribution does not work with x_gamma_extrapolation
  !
  exxdiv_treatment = 'none'
  IF (l_vcut) exxdiv_treatment = 'gb'
  !
  CALL mp_start_exx ( nband_, ntg_, intra_pool_comm )
  exx_grid_initialized = .FALSE.  
  CALL exx_grid_init()
  CALL exx_div_check()
  exxdiv = exx_divergence ()
  ! Compute the divergence and set the q+G=0 term (exxdiv, to be used in
  ! g2_convolution called by bare_pot.f90
  !
  WRITE (stdout,'(/,5X, "INFO: Divergence         ", 3x, 1A8)')       exxdiv_treatment
  WRITE (stdout,'(  5X, "INFO: Gamma Extrapolation", 3x, 1L5 )')      x_gamma_extrapolation
  WRITE (stdout,'(  5X, "INFO: Bare Coulomb G0    ", 3x, 1ES15.5 )')  exx_divergence()
  WRITE (stdout,'(  5X, "INFO: Screened Coulomb G0", 3x, 1ES15.5 )')  exx_divergence()/eps_inf
  !
  CALL stop_clock( 'Coulomb setup' )
  !
  RETURN
  !
end subroutine setup_coulomb_exx
