!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Taken from Phonon code (phonon.f90)
! Modified by G. Stenuit
!
!-----------------------------------------------------------------------
PROGRAM head
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the phonon program. It controls
  ! ... the initialization routines and the self-consistent cycle.
  ! ... At the end of the self-consistent run the dynamical matrix is
  ! ... computed. In the case q=0 the dielectric constant and the effective
  ! ... charges are computed.
  !
  ! This phonon pgrm has been modified to compute the head
  ! of the dielectric constant at all imaginary time
  ! USED in the context of the GWW calculations
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, ionode
  USE control_flags,   ONLY : conv_ions, modenum, twfcollect, restart
  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE klist,           ONLY : lgauss, nks
  USE start_k,         ONLY : xk_start, wk_start, nks_start
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE force_mod,       ONLY : force
  USE io_files,        ONLY : prefix, tmp_dir
  USE input_parameters,ONLY : pseudo_dir
  USE ions_base,       ONLY : nat
  USE noncollin_module,ONLY : noncolin
  USE scf,             ONLY : rho
  USE lsda_mod,        ONLY : nspin
  USE io_rho_xml,      ONLY : write_rho
  USE qpoint,          ONLY : xq, nksq, ikks, ikqs
  USE modes,           ONLY : nirr
  USE partial,         ONLY : done_irr
  USE disp,            ONLY : nqs, x_q, done_iq, rep_iq, done_rep_iq
  USE control_ph,      ONLY : ldisp, lnscf, lgamma, lgamma_gamma, convt, &
                              epsil, trans, zue, recover, rec_code, &
                              lnoloc, lrpa, done_bands,   &
                              start_q,last_q,start_irr,last_irr,current_iq,&
                              reduce_io, all_done, where_rec, tmp_dir_ph
  USE freq_ph
  USE output,          ONLY : fildyn, fildrho
  USE ramanm,          ONLY : lraman, elop
  USE check_stop,      ONLY : check_stop_init
  USE ph_restart,      ONLY : ph_readfile, ph_writefile, check_status_run, &
                              init_status_run, destroy_status_run
  USE save_ph,         ONLY : save_ph_input_variables, tmp_dir_save, &
                              restore_ph_input_variables, clean_input_variables
  !
!#ifdef __GWW
  USE wannier_gw,      ONLY : l_head
!#endif __GWW
!
  IMPLICIT NONE
  !
  INTEGER :: iq, iq_start, ierr, iu, ik
  INTEGER :: irr
  LOGICAL :: exst, do_band, exst_recover, exst_restart
  CHARACTER (LEN=9)   :: code = 'PHONON'
  CHARACTER (LEN=256) :: auxdyn
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  ! Initialize MPI, clocks, print initial messages
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( code )
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
     IF ( .NOT.(ldisp .OR. lnscf )) THEN
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
        lnscf = .TRUE.
        !
        IF (last_q<1.or.last_q>nqs) last_q=nqs
        !
        CALL init_status_run()
        !
     ELSE IF ( lnscf ) THEN
        !
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        x_q(:,1)=xq(:)
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
  IF ( lnscf ) CALL start_clock( 'PWSCF' )
  !
  DO iq = iq_start, last_q
     !
     IF (done_iq(iq)==1) CYCLE
     IF (start_irr>rep_iq(iq)) THEN
        WRITE(6,'(5x,"Exiting... start_irr,",i4,&
               & " > number of representations,",i4   )') &
               start_irr, rep_iq(iq)
        CYCLE
     ENDIF
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
        xq(1:3)  = x_q(1:3,iq)
        !
        lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
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
     do_band=(.NOT.trans)
     IF (.NOT. do_band) THEN
        DO irr=start_irr, MIN(ABS(last_irr),rep_iq(iq))
           IF (done_rep_iq(irr,iq) /= 1) THEN
               do_band=.TRUE.
               EXIT
           ENDIF
        ENDDO
     ENDIF
     !
     IF ( lnscf .AND.(.NOT.lgamma.OR.modenum /= 0) &
                .AND..NOT. done_bands) THEN
        !
        WRITE( stdout, '(/,5X,"Calculation of q = ",3F12.7)') xq
        !
        CALL clean_pw( .FALSE. )
        !
        CALL close_files(.TRUE.)
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
        CALL close_files(.TRUE.)
        !
     END IF
     !
     ! ... nksq is the number of k-points, NOT including k+q points
     !
     IF ( lgamma ) THEN
        !
        nksq = nks
        ALLOCATE(ikks(nksq), ikqs(nksq))
        DO ik=1,nksq
           ikks(ik) = ik
           ikqs(ik) = ik
        ENDDO
        !
     ELSE
        !
        nksq = nks / 2
        ALLOCATE(ikks(nksq), ikqs(nksq))
        DO ik=1,nksq
           ikks(ik) = 2 * ik - 1
           ikqs(ik) = 2 * ik
        ENDDO
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
     IF ( trans .AND..NOT.all_done ) CALL dynmat0()
     !
     IF ( epsil .AND. rec_code <=  1 ) THEN
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
        where_rec='after_diel'
        rec_code=2
        CALL ph_writefile('data',0)
     END IF
     !
!#ifdef __GWW
     IF(l_head) then
         call solve_head
     endif
!#endif __GWW
!
     ! ... cleanup of the variables
     !
     done_bands=.FALSE.
     done_iq(iq)=1
     DO irr=1,nirr
        IF (done_irr(irr)==0) done_iq(iq)=0
     ENDDO
     twfcollect=.FALSE.
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
  CALL destroy_status_run()
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
END PROGRAM head
