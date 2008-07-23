!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM phonon
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the phonon program. It controls
  ! ... the initialization routines and the self-consistent cycle.
  ! ... At the end of the self-consistent run the dynamical matrix is
  ! ... computed. In the case q=0 the dielectric constant and the effective
  ! ... charges are computed.
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, ionode
  USE control_flags,   ONLY : conv_ions
  USE klist,           ONLY : xqq, lgauss, nks
  USE basis,           ONLY : startingwfc, startingpot, startingconfig
  USE force_mod,       ONLY : force
  USE io_files,        ONLY : prefix, tmp_dir, nd_nmbr
  USE input_parameters, ONLY: pseudo_dir
  USE ions_base,       ONLY : nat
  USE noncollin_module, ONLY : noncolin
  USE start_k,         ONLY : xk_start, wk_start, nks_start
  USE control_flags,   ONLY : restart, lphonon, tr2, ethr, &
                              mixing_beta, david, isolve
  USE qpoint,          ONLY : xq, nksq
  USE modes,           ONLY : nirr
  USE partial,         ONLY : done_irr, comp_irr
  USE disp,            ONLY : nqs, x_q, done_iq
  USE control_ph,      ONLY : ldisp, lnscf, lgamma, lgamma_gamma, convt, &
                              epsil, trans, elph, zue, recover, rec_code, &
                              lnoloc, lrpa, done_bands, xml_not_of_pw,   &
                              start_q,last_q,start_irr,last_irr,current_iq,&
                              reduce_io
  USE freq_ph
  USE output,          ONLY : fildyn, fildrho
  USE global_version,  ONLY : version_number
  USE ramanm,          ONLY : lraman, elop
  USE check_stop,      ONLY : check_stop_init
  USE ph_restart,      ONLY : ph_readfile, ph_writefile
  USE save_ph,         ONLY : save_ph_input_variables,  &
                              restore_ph_input_variables, clean_input_variables
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, iq_start, ierr, iu
  INTEGER :: irr
  LOGICAL :: exst, do_band
  CHARACTER (LEN=9)   :: code = 'PHONON'
  CHARACTER (LEN=256) :: auxdyn
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
#if defined __INTEL
  ! ... Intel compilers v .ge.8 allocate a lot of stack space
  ! ... Stack limit is often small, thus causing SIGSEGV and crash
  CALL remove_stack_limit ( )
#endif
  !
  CALL init_clocks( .TRUE. )
  !
  CALL start_clock( 'PHONON' )
  !
  CALL startup( nd_nmbr, code, version_number )
  !
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
  !
  ! ... and begin with the initialization part
  !
  CALL phq_readin()
  !
  CALL save_ph_input_variables()
  !
  CALL check_stop_init()
  !
  ! ... Checking the status of the calculation
  !
  IF (recover) THEN
     CALL ph_readfile('init',ierr)
     CALL check_restart_recover(iq_start,start_q,current_iq)
     IF ( .NOT.(ldisp .OR. lnscf )) THEN
        last_q=1
     ELSEIF (ierr == 0) THEN
        IF (last_q<1.OR.last_q>nqs) last_q=nqs
        IF (ldisp) auxdyn = fildyn
     ENDIF
     IF (ierr /= 0) recover=.FALSE.
  ELSE
     ierr=1
  ENDIF
  IF (ierr /= 0) THEN
     !
     ! recover file not found or not looked for
     !
     done_bands=.FALSE.
     xml_not_of_pw=.FALSE.
     iq_start=start_q
     IF (ldisp) THEN
        !
        ! ... Calculate the q-points for the dispersion
        !
        CALL q_points()
        !
        ! ... Store the name of the matdyn file in auxdyn
        !
        auxdyn = fildyn
        !
        ! ... do always a non-scf calculation
        !
        lnscf = .TRUE.
        !
        IF (last_q<1.or.last_q>nqs) last_q=nqs
        !
        ALLOCATE(done_iq(nqs))
        done_iq=0
        !
      ELSE IF ( lnscf ) THEN
        !
        ! ... xq  is the q-point for   phonon calculation (read from input)
        ! ... xqq is the q-point for the nscf calculation (read from data file)
        ! ... if the nscf calculation is to be performed, discard the latter
        !
        xqq = xq
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        x_q(:,1)=xqq(:)
        ALLOCATE(done_iq(1))
        done_iq=0
        !
     ELSE
        !
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        x_q(:,1)=xq(:)
        ALLOCATE(done_iq(1))
        done_iq=0
        !
     END IF
  END IF
  !
  IF (nks_start==0) CALL errore('phonon','wrong starting k',1)
  !
  IF ( lnscf ) CALL start_clock( 'PWSCF' )
  !
  DO iq = iq_start, last_q
     !
     IF (done_iq(iq)==1) CYCLE
     !
     current_iq=iq
     !
     IF ( ldisp ) THEN
        !
        ! ... set the name for the output file
        !
        fildyn = TRIM( auxdyn ) // TRIM( int_to_char( iq ) )
        !
        ! ... set the q point
        !
        xqq(1:3) = x_q(1:3,iq)
        xq(1:3)  = x_q(1:3,iq)
        !
        lgamma = ( xqq(1) == 0.D0 .AND. xqq(2) == 0.D0 .AND. xqq(3) == 0.D0 )
        !
        IF ( lgamma ) THEN
           !
           IF ( .NOT. lgauss ) THEN
              !
              ! ... in the case of an insulator at q=0 one has to calculate 
              ! ... the dielectric constant and the Born eff. charges
              !
              epsil = .TRUE.
              zue   = .TRUE.
              !
           ELSE
              !
              epsil = .FALSE.
              zue   = .FALSE.
              !
           END IF
           !
        ELSE
           !
           ! ... for q != 0 no calculation of the dielectric tensor 
           ! ...           and Born eff. charges
           !
           epsil = .FALSE.
           zue   = .FALSE.
           !
           ! ... non-scf calculation needed:
           ! ... reset the k-points to their starting values. Note that
           ! ... in LSDA case k-points are already doubled to account for
           ! ... spin polarization: restore the original number of k-points
           !
        END IF
     ENDIF
     !
     !  Save the current status of the run
     !
     CALL ph_writefile('init',0)
     !
     ! ... In the case of q != 0, we make first a non selfconsistent run
     !
     do_band=(start_irr /= 0).OR.(last_irr /= 0)
     IF ( lnscf .AND.(.NOT.lgamma.OR.xml_not_of_pw) &
                .AND..NOT. done_bands.and.do_band) THEN
        !
        !
        WRITE( stdout, '(/,5X,"Calculation of q = ",3F12.7)') xqq
        !
        CALL clean_pw( .FALSE. )
        !
        CALL close_files()
        !
        ! ... Setting the values for the nscf run
        !
        CALL set_defaults_pw()
        lphonon           = .TRUE.
        startingconfig    = 'input'
        startingpot       = 'file'
        startingwfc       = 'atomic'
        restart = recover
        pseudo_dir= TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
        CALL restart_from_file()
        conv_ions=.true.
        !
        ! ... the threshold for diagonalization ethr is calculated via
        ! ... the threshold on self-consistency tr2 - the value used
        ! ... here should be good enough for all cases
        !
        tr2 = 1.D-9
        ethr = 0.d0
        mixing_beta = 0.d0
        !
        ! ... Assume davidson diagonalization
        !
        isolve = 0
        david = 4
        !
        IF ( .NOT. ALLOCATED( force ) ) ALLOCATE( force( 3, nat ) )
        !
        CALL init_run()
        !
        CALL electrons()
        !
        IF (.NOT.reduce_io) THEN
           write(6,*) 'call punch'
           CALL punch( 'all' )
           write(6,*) 'done punch'
           done_bands=.TRUE.
           xml_not_of_pw=.TRUE.
        ENDIF
        !
        CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
        !
        CALL close_files()
        !
     END IF
     !
     ! ... Setting nksq
     !
     IF ( lgamma ) THEN
        !
        nksq = nks
        !
     ELSE
        !
        nksq = nks / 2
        !
     END IF
     !
     CALL ph_writefile('init',0)
     !
     ! ... Calculation of the dispersion: do all modes 
     !
     CALL allocate_phq()
     !
     !  read the displacement patterns if available in the recover file
     !
     rec_code=0
     IF (recover) CALL ph_readfile('data_u',ierr)
     CALL phq_setup()
     CALL phq_recover()
     CALL phq_summary()
     !
     CALL openfilq()
     !
     CALL phq_init()
     !
     CALL print_clock( 'PHONON' )
     !
     IF ( trans .AND. (done_irr(0)==0.AND.comp_irr(0)==1) ) CALL dynmat0()
     !
     IF ( epsil .AND. rec_code <=  0 ) THEN
        !
        IF (fpol) THEN    ! calculate freq. dependent polarizability
           !
           WRITE( stdout, '(/,5X,"Frequency Dependent Polarizability Calculation",/)' )
           !
           iu = nfs
           !
           freq_loop : DO WHILE ( iu .gt. 0)
              !
              CALL solve_e_fpol( fiu(iu) )
              IF ( convt ) CALL polariz ( fiu(iu) )
              iu = iu - 1
              !
           END DO freq_loop
           !
           WRITE( stdout, '(/,5X,"End of Frequency Dependent Polarizability Calculation")' )
           !
        ENDIF
        !
        WRITE( stdout, '(/,5X,"Electric Fields Calculation")' )
        !
        CALL solve_e()
        !
        WRITE( stdout, '(/,5X,"End of electric fields calculation")' )
        !
        IF ( convt ) THEN
           !
           ! ... calculate the dielectric tensor epsilon
           !
           CALL dielec()
           !
           ! ... calculate the effective charges Z(E,Us) (E=scf,Us=bare)
           !
           IF (.NOT.(lrpa.OR.lnoloc)) CALL zstar_eu()
           !
           IF ( fildrho /= ' ' ) CALL punch_plot_e()
           !
        ELSE
           !
           CALL stop_ph( .FALSE. )
           !
        END IF
        !
        IF (( lraman .OR. elop ).AND..NOT.noncolin) CALL raman()
        !
     END IF
     !
     IF ( trans ) THEN
        !
        CALL phqscf()
        CALL dynmatrix()
        !
        IF ( fildrho /= ' ' ) CALL punch_plot_ph()
        !
     END IF
     !
     IF ( elph ) THEN
        !
        IF (noncolin) CALL errore('phonon','e-ph and noncolin not programed',1)
        IF ( .NOT. trans ) THEN
           ! 
           CALL dvanqq()
           CALL elphon()
           !
        END IF
        !
        CALL elphsum()
        !
     END IF
     !
     ! ... cleanup of the variables
     !
     done_bands=.FALSE.
     done_iq(iq)=1
     DO irr=1,nirr
        IF (done_irr(irr)==0) done_iq(iq)=0
     ENDDO
     CALL clean_pw( .FALSE. )
     CALL deallocate_phq()
     !
     ! ... Close the files
     !
     CALL close_phq( .TRUE. )
     !
     CALL restore_ph_input_variables()
     !
  END DO

  CALL ph_writefile('init',0)
  CALL clean_input_variables()
  !
  IF ( ALLOCATED( xk_start ) ) DEALLOCATE( xk_start )
  IF ( ALLOCATED( wk_start ) ) DEALLOCATE( wk_start )
  !
  IF ( lnscf ) CALL print_clock_pw()
  !
  CALL stop_ph( .TRUE. )
  !
  STOP
  !
END PROGRAM phonon
