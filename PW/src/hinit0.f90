!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hinit0()
  !-----------------------------------------------------------------------
  !! Hamiltonian initialization: atomic position independent initialization
  !! for nonlocal PP, structure factors, local potential, core charge.
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, nsp, ityp, tau
  USE basis,            ONLY : startingconfig
  USE cell_base,        ONLY : alat, at, bg, omega
  USE cellmd,           ONLY : omega_old, at_old, lmovecell, calc
  USE dynamics_module,  ONLY : verlet_read_tau_from_conf
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ecutrho, ngm, g, gg, eigts1, eigts2, eigts3
#if defined (__CUDA)
  USE gvect,            ONLY : eigts1_d, eigts2_d, eigts3_d
#endif
  USE gvecw,            ONLY : ecutwfc
  USE vlocal,           ONLY : strf
  USE realus,           ONLY : generate_qpointlist, betapointlist, &
                               init_realspace_vars, real_space
  USE ldaU,             ONLY : lda_plus_U, Hubbard_projectors
  USE control_flags,    ONLY : tqr, tq_smoothing, tbeta_smoothing, restart
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : report
  USE mp_bands,         ONLY : intra_bgrp_comm
  !
#if defined (__ENVIRON)
  USE plugin_flags,        ONLY : use_environ
  USE environ_base_module, ONLY : update_environ_ions, update_environ_cell
#endif
#if defined (__OSCDFT)
  USE plugin_flags,        ONLY : use_oscdft
  USE oscdft_base,         ONLY : oscdft_ctx
  USE oscdft_context,      ONLY : oscdft_init_context
#endif
  !
  IMPLICIT NONE
  REAL (dp) :: alat_old
  LOGICAL   :: is_tau_read = .FALSE.
  !
#if defined (__ENVIRON)
  REAL(DP) :: at_scaled(3, 3)
  REAL(DP) :: tau_scaled(3, nat)
#endif
  !
  CALL start_clock( 'hinit0' )
  !
  ! ... calculate the Fourier coefficients of the local part of the PP
  !
  CALL init_vloc()
  !
  ! ... k-point independent parameters of non-local pseudopotentials
  !
  IF (tbeta_smoothing) CALL init_us_b0(ecutwfc,intra_bgrp_comm)
  IF (tq_smoothing) CALL init_us_0(ecutrho,intra_bgrp_comm)
  CALL init_us_1(nat, ityp, omega, ngm, g, gg, intra_bgrp_comm)
  IF ( lda_plus_U .AND. ( Hubbard_projectors == 'pseudo' ) ) CALL init_q_aeps()
  CALL init_tab_atwfc (omega, intra_bgrp_comm)
  !
  IF ( restart .AND. startingconfig == 'file' ) THEN
     !
     IF ( lmovecell ) THEN
        !
        ! ... store initial values of at and omega, used for G-vectors etc.
        !
        at_old    = at
        omega_old = omega
        !
        CALL read_conf_from_file( lmovecell, nat, nsp, tau, alat, at, &
                                  is_tau_read )
        CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
        CALL volume (alat, at(:,1), at(:,2), at(:,3), omega)
        CALL scale_h( )
        !
     ELSE
        !
        CALL read_conf_from_file( lmovecell, nat, nsp, tau, alat_old, at_old, &
                                  is_tau_read )
        !
        IF (.NOT.is_tau_read .AND. calc == 'vd') THEN
           !
           ! Restart of Verlet MD is requested. Failed to read coordinates from
           ! XML. Try to restart from the .md file.
           !
           CALL verlet_read_tau_from_conf()
           !
        END IF
        !
     END IF
     !
  ENDIF
  !
  ! ... initialize the structure factor
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                   dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                   strf, eigts1, eigts2, eigts3 )
  ! sync duplicated version
#if defined(__CUDA)
  eigts1_d = eigts1
  eigts2_d = eigts2
  eigts3_d = eigts3
#endif
  !$acc update device(eigts1, eigts2, eigts3) 
  !
  ! these routines can be used to patch quantities that are dependent
  ! on the ions and cell parameters
  !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_init_ions(tau) 
  CALL plugin_init_cell() 
#endif 
#if defined (__ENVIRON)
  IF (use_environ) THEN
     at_scaled = at * alat
     tau_scaled = tau * alat
     CALL update_environ_ions(tau_scaled)
     CALL update_environ_cell(at_scaled)
  END IF
#endif
#if defined (__OSCDFT)
  IF (use_oscdft) THEN
     CALL oscdft_init_context(oscdft_ctx)
  END IF
#endif
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()
  !
  ! ... calculate the core charge (if any) for the nonlinear core correction
  !
  CALL set_rhoc()
  !
  ! ... more position-dependent initializations
  !
  IF ( tqr ) CALL generate_qpointlist()
  !
  IF (real_space ) THEN
     CALL betapointlist()
     CALL init_realspace_vars()
     WRITE(stdout,'(5X,"Real space initialisation completed")')    
  ENDIF
  !
  IF ( report /= 0 ) CALL make_pointlists( )
  !
  CALL stop_clock( 'hinit0' )
  !
  RETURN
  !
END SUBROUTINE hinit0

