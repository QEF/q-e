!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM phonon
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the phonon code.
  ! ... It reads all the quantities calculated by pwscf, it
  ! ... checks if some recover file is present and determines
  ! ... which calculation needs to be done. Finally, it makes
  ! ... a loop over the q points. At a generic q, if necessary it
  ! ... recalculates the band structure calling pwscf again.
  ! ... Then it can calculate the response to an atomic displacement,
  ! ... the dynamical matrix at that q, and the electron-phonon 
  ! ... interaction at that q. At q=0 it can calculate the linear response
  ! ... to an electric field perturbation and hence the dielectric
  ! ... constant, the Born effective charges and the polarizability
  ! ... at imaginary frequencies. 
  ! ... At q=0, from the second order response to an electric field,
  ! ... it can calculate also the electro-optic and the raman tensors.
  ! ... Presently implemented: 
  ! ... dynamical matrix (q/=0)   NC [4], US [4], PAW [3]
  ! ... dynamical matrix (q=0)    NC [5], US [5], PAW [3]
  ! ... dielectric constant       NC [5], US [5], PAW [3] 
  ! ... born effective charges    NC [5], US [5], PAW [3]
  ! ... polarizability (iu)       NC [2], US [2]
  ! ... elctron-phonon            NC [3], US [3]
  ! ... electro-optic             NC [1]
  ! ... raman tensor              NC [1]
  !
  ! NC = norm conserving pseudopotentials
  ! US = ultrasoft pseudopotentials
  ! PAW = projector augmented-wave
  ! [1] LDA, [2] [1]+GGA, [3] [2]+LSDA/sGGA, [4] [3]+Spin-orbit/nonmagnetic,
  ! [5] [4]+Spin-orbit/magnetic
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, ionode
  USE control_flags,   ONLY : conv_ions, modenum, twfcollect
  USE klist,           ONLY : lgauss, nks
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE force_mod,       ONLY : force
  USE io_files,        ONLY : prefix, tmp_dir, nd_nmbr
  USE input_parameters,ONLY : pseudo_dir
  USE ions_base,       ONLY : nat
  USE start_k,         ONLY : xk_start, wk_start, nks_start
  USE noncollin_module,ONLY : noncolin
  USE control_flags,   ONLY : restart
  USE scf,             ONLY : rho
  USE lsda_mod,        ONLY : nspin
  USE io_rho_xml,      ONLY : write_rho
  USE qpoint,          ONLY : xq, nksq
  USE modes,           ONLY : nirr
  USE partial,         ONLY : done_irr
  USE disp,            ONLY : nqs, x_q, done_iq, rep_iq, done_rep_iq
  USE control_ph,      ONLY : ldisp, lgamma, lgamma_gamma, convt, &
                              epsil, trans, elph, zue, recover, rec_code, &
                              lnoloc, lrpa, done_bands,   &
                              start_q,last_q,start_irr,last_irr,current_iq,&
                              reduce_io, all_done, where_rec, tmp_dir_ph
  USE freq_ph
  USE output,          ONLY : fildyn, fildrho
  USE global_version,  ONLY : version_number
  USE ramanm,          ONLY : lraman, elop
  USE check_stop,      ONLY : check_stop_init
  USE ph_restart,      ONLY : ph_readfile, ph_writefile, check_status_run, &
                              init_status_run, destroy_status_run
  USE save_ph,         ONLY : save_ph_input_variables, tmp_dir_save, &
                              restore_ph_input_variables, clean_input_variables
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, iq_start, ierr, iu, ik
  INTEGER :: irr
  LOGICAL :: do_band, do_iq, setup_pw
  LOGICAL :: exst, exst_recover, exst_restart
  CHARACTER (LEN=9)   :: code = 'PHONON'
  CHARACTER (LEN=256) :: auxdyn
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  ! ... Intel compilers v .ge.8 allocate a lot of stack space
  ! ... Stack limit is often small, thus causing SIGSEGV and crash
  CALL remove_stack_limit ( )
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
  CALL check_stop_init()
  !
  ! ... Checking the status of the calculation
  !
  tmp_dir=tmp_dir_ph
  IF (recover) THEN
     ierr=0
     CALL check_restart_recover(exst_recover, exst_restart)
     IF (.NOT.exst_recover.AND..NOT.exst_restart) THEN
        iq_start=start_q
     ELSE
        iq_start=current_iq
     ENDIF
     IF (ierr == 0 )CALL check_status_run()
     IF ( .NOT.(ldisp)) THEN
        last_q=1
     ELSEIF (ierr == 0) THEN
        IF (last_q<1.OR.last_q>nqs) last_q=nqs
        IF (ldisp) auxdyn = fildyn
     ENDIF
     IF (ierr /= 0) THEN
        recover=.FALSE.
     ELSE
        WRITE(stdout, &
            '(5x,i4," /",i4," q-points for this run, from", i3,&
               & " to", i3,":")') last_q-iq_start+1, nqs, iq_start, last_q
        WRITE(stdout, '(5x,"  N       xq(1)       xq(2)       xq(3) " )')
        DO iq = 1, nqs
           WRITE(stdout, '(5x,i3, 3f12.7,l6)') iq, x_q(1,iq), x_q(2,iq), &
                                x_q(3,iq)
        END DO
        WRITE(stdout, *)
     ENDIF
  ELSE
     ierr=1
  ENDIF
  !
  ! We copy the charge density in the directory with the _ph prefix
  ! to calculate the bands
  !
  IF (ldisp.OR..NOT.lgamma.OR.modenum/=0) CALL write_rho( rho, nspin )
  !
  CALL save_ph_input_variables()
  !
  IF (ierr /= 0) THEN
     !
     ! recover file not found or not looked for
     !
     done_bands=.FALSE.
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
        IF (last_q<1.or.last_q>nqs) last_q=nqs
        !
        CALL init_status_run()
        !
     ELSE 
        !
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        x_q(:,1)=xq(:)
        CALL init_status_run()
        !
     END IF
  END IF
  !
  IF (nks_start==0) CALL errore('phonon','wrong starting k',1)
  !
  CALL start_clock( 'PWSCF' )
  !
  DO iq = iq_start, last_q
     !
     CALL prepare_q(auxdyn, do_band, do_iq, setup_pw, iq)
     !
     !  If this q is not done in this run, cycle
     !
     IF (.NOT.do_iq) CYCLE
     !
     IF (setup_pw) THEN
        !
        CALL clean_pw( .FALSE. )
        !
        CALL close_files()
        !
        ! ... Setting the values for the nscf run
        !
        tmp_dir=tmp_dir_ph
        startingconfig    = 'input'
        starting_pot      = 'file'
        starting_wfc      = 'atomic'
        restart = recover
        pseudo_dir= TRIM( tmp_dir_save ) // TRIM( prefix ) // '.save'
        CALL restart_from_file()
        conv_ions=.true.
        !
        IF ( .NOT. ALLOCATED( force ) ) ALLOCATE( force( 3, nat ) )
        !
        CALL setup_nscf (xq)
        CALL init_run()
        !
        IF (do_band) CALL electrons()
        !
        IF (.NOT.reduce_io.and.do_band) THEN
           twfcollect=.FALSE. 
           CALL punch( 'all' )
           done_bands=.TRUE.
        ENDIF
        !
        CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
        !
        CALL close_files()
        !
     END IF
     !
     !  Initialization of the phonon calculation
     !
     CALL initialize_ph()
     !
     IF ( trans .AND..NOT.all_done ) CALL dynmat0()
     !
     !  electric field perturbation
     !
     IF (epsil) CALL phescf()
     !
     !  phonon perturbation
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
     ! ... cleanup of the variables for the next q point
     !
     CALL clean_pw_ph(iq)
     !
  END DO

  CALL ph_writefile('init',0)
  CALL clean_input_variables()
  CALL destroy_status_run()
  !
  IF ( ALLOCATED( xk_start ) ) DEALLOCATE( xk_start )
  IF ( ALLOCATED( wk_start ) ) DEALLOCATE( wk_start )
  !
  CALL print_clock_pw()
  !
  CALL stop_ph( .TRUE. )
  !
  STOP
  !
END PROGRAM phonon
