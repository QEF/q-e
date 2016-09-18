!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock_lr()
   !---------------------------------------------------------------------------
   !
   ! This subroutine prints out the clocks at the end of the run.
   ! It constructs the calling tree of the program.
   !
   USE io_global,        ONLY : stdout
   USE mp_world,         ONLY : mpime, root
   USE realus,           ONLY : real_space,real_space_debug
   USE lr_variables,     ONLY : davidson, eels
   USE funct,            ONLY : dft_is_hybrid
#if defined(__ENVIRON)
   USE plugin_flags,     ONLY : use_environ
   USE environ_info,     ONLY : environ_clock
#endif
   !
   IMPLICIT NONE
   !
   WRITE( stdout, * )
   !
   IF (eels) THEN
      CALL print_clock( 'lr_eels_main' )
   ELSEIF (davidson) THEN
      CALL print_clock( 'lr_dav_main' )
   ELSE
      CALL print_clock( 'lr_main' )
   ENDIF   
   !
   IF (.NOT.eels) CALL print_clock( 'read_wf' )
   CALL print_clock( 'lr_solve_e' )
   !
   IF (davidson) THEN
     CALL print_clock( 'calc_residue' )
     CALL print_clock( 'expan_basis' )
     CALL print_clock( 'matrix')
     CALL print_clock( 'mGS_orth' )
     CALL print_clock( 'mGS_orth_pp' )
   ENDIF
   !
   CALL print_clock( 'one_step' )
   !
   WRITE( stdout, * )
   !
   CALL print_clock('lr_apply')
   CALL print_clock('lr_apply_int')
   CALL print_clock('lr_apply_no')
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'h_psi' )
   CALL print_clock( 'lr_calc_dens' )
   IF (eels) CALL print_clock( 'incdrhoscf' )
   CALL print_clock( 'lr_dvpsi_e' )
   CALL print_clock( 'lr_dv_setup' )
   CALL print_clock( 'dv_of_drho' )
   CALL print_clock( 'interaction' )
   CALL print_clock( 'lr_dot' )
   CALL print_clock( 'ortho' )
   IF (davidson) CALL print_clock( 'lr_ortho' )
   !
   WRITE( stdout, * ) 
   CALL print_clock( 'lr_exx_int')
   CALL print_clock( 'lr_exx_noint')
   !
   WRITE( stdout, * )
   WRITE( stdout, '(5X,"US routines")' )
   !
   CALL print_clock( 's_psi' )
   CALL print_clock( 'sd0psi' )
   CALL print_clock( 'lr_apply_s' )
   CALL print_clock( 'lr_dot_us' )
   IF (eels) THEN
    CALL print_clock( 'addusdbec' )
    CALL print_clock( 'addusdbec_nc' )
    CALL print_clock( 'lr_addusddens' )
    CALL print_clock( 'lr_addus_dvpsi' )
   ENDIF
   IF (eels) THEN
      CALL print_clock( 'lr_sm1_psiq' )
   ELSE
      CALL print_clock( 'lr_sm1_psi' )
   ENDIF
   !
   IF (real_space_debug>0) THEN
    WRITE( stdout, '(5X,"US routines, RS")' )
    CALL print_clock ( 'realus' )
    CALL print_clock ( 'betapointlist' )
    CALL print_clock ( 'calbec_rs' )
    CALL print_clock ( 's_psir' )
    CALL print_clock ( 'add_vuspsir' )
    CALL print_clock ( 'invfft_orbital' )
    CALL print_clock ( 'fwfft_orbital' )
    CALL print_clock ( 'v_loc_psir' )
   ENDIF
   !
   WRITE( stdout, * )
   WRITE( stdout, '(5X,"General routines")' )
   !
   CALL print_clock( 'calbec' )
   CALL print_clock( 'fft' )
   CALL print_clock( 'ffts' )
   CALL print_clock( 'fftc' )
   CALL print_clock( 'fftw' )
   CALL print_clock( 'fftcw' )
   CALL print_clock( 'interpolate' )
   CALL print_clock( 'davcio' )
   CALL print_clock( 'newq' )
   !
   WRITE( stdout, * )
   !
#if defined (__MPI)
   WRITE( stdout, '(5X,"Parallel routines")' )
   CALL print_clock( 'fft_scatter' )
   CALL print_clock ('mp_sum')
   WRITE( stdout, * )
#endif
   !
#if defined(__ENVIRON)
   IF ( use_environ ) CALL environ_clock( stdout )
#endif
   !
   IF (dft_is_hybrid()) THEN
    !
    WRITE( stdout, '(5X,"EXX routines")' )
    CALL print_clock( 'exx_grid' )
    CALL print_clock( 'exxinit' )
    CALL print_clock( 'vexx' )
    CALL print_clock( 'exxenergy' )
    CALL print_clock( 'exxen2' )
    CALL print_clock ('cycleig')
    WRITE( stdout, * )
    !
   ENDIF
   !
   IF (eels) THEN
      !
      WRITE( stdout, '(5X,"EELS routines")' )
      CALL print_clock( 'lr_run_nscf' )
      CALL print_clock( 'lr_setup_nscf' )
      CALL print_clock( 'lr_calc_dens' )
      CALL print_clock( 'lr_dvpsi_eels' )
      CALL print_clock( 'lr_sym_eels' )
      CALL print_clock( 'lr_psym_eels' )
      CALL print_clock( 'lr_smallgq' )
      CALL print_clock( 'lr_summary' )
      WRITE( stdout, * )
      !
   ENDIF
   !
   CALL print_clock( 'post-processing' )
   !
   RETURN
   !
END SUBROUTINE print_clock_lr
