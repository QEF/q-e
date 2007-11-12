!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock_pw()
   !---------------------------------------------------------------------------
   !
   ! ... this routine prints out the clocks at the end of the run
   ! ... it tries to construct the calling tree of the program.
   !
   USE io_global,          ONLY : stdout
   USE control_flags,      ONLY : isolve
   USE force_mod,          ONLY : lforce, lstres
   USE mp_global,          ONLY : mpime, root
   USE paw_variables,      ONLY : okpaw
   !
   IMPLICIT NONE
   !
   !
   IF ( mpime /= root ) &
      OPEN( UNIT = stdout, FILE = '/dev/null', STATUS = 'UNKNOWN' )
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'PWSCF' )
   CALL print_clock( 'init_run' )
   CALL print_clock( 'electrons' )
   CALL print_clock( 'update_pot' )
   !
   IF ( lforce ) CALL print_clock( 'forces' )
   IF ( lstres ) CALL print_clock( 'stress' )
   !
   WRITE( stdout, * )
   CALL print_clock( 'qpointlist' )
   CALL print_clock( 'realus:boxes' )
   CALL print_clock( 'realus:spher' )
   CALL print_clock( 'realus:qsave' )
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'electrons' )
   CALL print_clock( 'c_bands' )
   CALL print_clock( 'sum_band' )
   CALL print_clock( 'v_of_rho' )
   CALL print_clock( 'v_h' )
   CALL print_clock( 'v_xc' )
   CALL print_clock( 'newd' )
   !
#ifdef DEBUG_NEWD
   !
   CALL print_clock( 'newd:fftvg' )
   CALL print_clock( 'newd:qvan2' )
   CALL print_clock( 'newd:int1' )
   CALL print_clock( 'newd:int2' )
#endif
   !
   CALL print_clock( 'mix_rho' )
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'c_bands' )
   CALL print_clock( 'init_us_2' )
   CALL print_clock( 'cegterg' )
   CALL print_clock( 'ccgdiagg' )
   CALL print_clock( 'diis' )
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'sum_band' )
   CALL print_clock( 'becsum' )
   !
   CALL print_clock( 'addusdens' )
   !
#ifdef DEBUG_ADDUSDENS
   CALL print_clock( 'addus:qvan2' )
   CALL print_clock( 'addus:strf' )
   CALL print_clock( 'addus:aux2' )
   CALL print_clock( 'addus:aux' )
#endif
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'wfcrot' )
   CALL print_clock( 'wfcrot1' )
   CALL print_clock( 'cegterg' )
   CALL print_clock( 'ccdiagg' )
   CALL print_clock( 'cdiisg' )
   !
   IF ( isolve == 0 ) THEN
      !
      CALL print_clock( 'h_psi' )
      CALL print_clock( 'g_psi' )
      CALL print_clock( 'overlap' )
      CALL print_clock( 'diagh' )
      CALL print_clock( 'diaghg' )
      CALL print_clock( 'choldc' )
      CALL print_clock( 'inversion' )
      CALL print_clock( 'paragemm' )
      !
      CALL print_clock( 'update' )
      CALL print_clock( 'last' )
      !
      WRITE( stdout, * )
      !
      CALL print_clock( 'h_psi' )
      CALL print_clock( 'init' )
      CALL print_clock( 'firstfft' )
      CALL print_clock( 'secondfft' )
      CALL print_clock( 'add_vuspsi' )
      CALL print_clock( 'h_psi_meta' )
      CALL print_clock( 's_psi' )
      !
   ELSE IF ( isolve == 1 ) THEN
      !
      CALL print_clock( 'h_1psi' )
      CALL print_clock( 's_1psi' )
      CALL print_clock( 'cdiaghg' )
      !
      WRITE( stdout, * )
      !
      CALL print_clock( 'h_1psi' )
      CALL print_clock( 'init' )
      CALL print_clock( 'firstfft' )
      CALL print_clock( 'secondfft' )
      CALL print_clock( 'add_vuspsi' )
      CALL print_clock( 'h_psi_meta' )
      !
   ELSE
      !
      CALL print_clock( 'h_psi' )
      CALL print_clock( 's_psi' )
      CALL print_clock( 'g_psi' )
      CALL print_clock( 'cdiaghg' )
      CALL print_clock( 'cholortho' )
      !
      WRITE( stdout, * )
      !
      CALL print_clock( 'h_psi' )
      CALL print_clock( 'init' )
      CALL print_clock( 'firstfft' )
      CALL print_clock( 'secondfft' )
      CALL print_clock( 'add_vuspsi' )
      CALL print_clock( 'h_psi_meta' )
      !   
   END IF
   !
   WRITE( stdout, * )
   WRITE( stdout, '(5X,"General routines")' )
   !
   CALL print_clock( 'ccalbec' )
   CALL print_clock( 'cft3' )
   CALL print_clock( 'cft3s' )
   CALL print_clock( 'interpolate' )
   CALL print_clock( 'davcio' )
   !    
   WRITE( stdout, * )
   !
#if defined (__PARA)
   WRITE( stdout, '(5X,"Parallel routines")' )
   !
   CALL print_clock( 'reduce' )
   CALL print_clock( 'fft_scatter' )
   CALL print_clock( 'poolreduce' )
#endif
   !
#ifdef EXX
   WRITE( stdout, '(5X,"EXX routines")' )
   !
   CALL print_clock( 'exx_grid' )
   CALL print_clock( 'exxinit' )
   CALL print_clock( 'vexx' )
   CALL print_clock( 'exxenergy' )
   CALL print_clock( 'exxen2' )
   CALL print_clock ('cycleig')
#endif
   !
   IF ( okpaw ) THEN
      WRITE( stdout, * )
      WRITE( stdout, '(5X,"PAW routines")' )
      ! radial routines:
      CALL print_clock ('PAW_pot')
      CALL print_clock ('PAW_newd')
      CALL print_clock ('PAW_int')
      CALL print_clock ('PAW_ddot')
      CALL print_clock ('PAW_rad_init')
      CALL print_clock ('PAW_energy')
      CALL print_clock ('PAW_symme')
      ! second level routines:
      CALL print_clock ('PAW_rho_lm')
      CALL print_clock ('PAW_h_pot')
      CALL print_clock ('PAW_xc_pot')
      CALL print_clock ('PAW_lm2rad')
      CALL print_clock ('PAW_rad2lm')
      ! third level, or deeper:
      CALL print_clock ('PAW_gcxc_v')
      CALL print_clock ('PAW_div')
      CALL print_clock ('PAW_grad')
   END IF
   !
   RETURN
   !
END SUBROUTINE print_clock_pw
