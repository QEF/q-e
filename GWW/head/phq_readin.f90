!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Taken from Phonon code
! Modified by P. Umari and G. Stenuit
!
!----------------------------------------------------------------------------
SUBROUTINE phq_readin()
  !----------------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program phononq.
  !    from standard input (unit 5).
  !    A second routine readfile reads the variables saved on a file
  !    by the self-consistent program.
  !
  !
  USE kinds,         ONLY : DP
  USE parameters,    ONLY : nsx
  USE constants,     ONLY : amconv
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_bcast
  USE input_parameters, ONLY : max_seconds
  USE ions_base,     ONLY : amass, atm
  USE klist,         ONLY : xk, nks, nkstot, lgauss, two_fermi_energies
  USE control_flags, ONLY : gamma_only, tqr, restart, lkpoint_dir
  USE uspp,          ONLY : okvan
  USE fixed_occ,     ONLY : tfixed_occ
  USE lsda_mod,      ONLY : lsda, nspin
  USE run_info, ONLY : title
  USE control_ph,    ONLY : maxter, alpha_mix, lgamma, lgamma_gamma, epsil, &
                            zue, trans, reduce_io, nogg, &
                            elph, tr2_ph, niter_ph, nmix_ph, lnscf, &
                            ldisp, recover, lrpa, lnoloc, start_irr, &
                            last_irr, start_q, last_q, current_iq, tmp_dir_ph
  USE save_ph,       ONLY : tmp_dir_save
  USE ph_restart,    ONLY : ph_readfile
  USE gamma_gamma,   ONLY : asr
  USE qpoint,        ONLY : nksq, xq
  USE partial,       ONLY : atomo, nat_todo
  USE output,        ONLY : fildyn, fildvscf, fildrho
  USE disp,          ONLY : nq1, nq2, nq3
  USE io_files,      ONLY : tmp_dir, prefix
  USE noncollin_module, ONLY : i_cons, noncolin
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : iverbosity, modenum
  USE io_global,     ONLY : ionode
  USE mp_global,     ONLY : nproc, nproc_pool, nproc_file, nproc_pool_file
  USE control_flags, ONLY : twfcollect
  USE paw_variables, ONLY : okpaw
  USE ramanm,        ONLY : eth_rps, eth_ns, lraman, elop, dek
  USE freq_ph,       ONLY : fpol, fiu, nfs
  !
!#ifdef __GWW
!added the varibles for the calculation of the head and wings
  USE wannier_gw,    ONLY : l_head, omega_gauss, n_gauss, grid_type
!#endif
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios, ipol, iter, na, it
    ! integer variable for I/O control
    ! counter on polarizations
    ! counter on iterations
    ! counter on atoms
    ! counter on types
  REAL(DP) :: amass_input(nsx)
    ! save masses read from input here
  CHARACTER (LEN=256) :: outdir
  !
  CHARACTER(LEN=256)         :: input_line
  CHARACTER(LEN=80)          :: card
  CHARACTER(LEN=1), EXTERNAL :: capital
  LOGICAL                    :: tend
  LOGICAL                    :: end_of_file
  LOGICAL                    :: exst_restart, exst_recover
  INTEGER                    :: i, ierr
  INTEGER, EXTERNAL :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight

  !
  NAMELIST / INPUTPH / tr2_ph, amass, alpha_mix, niter_ph, nmix_ph,  &
                       nat_todo, iverbosity, outdir, epsil,  &
                       trans, elph, zue, max_seconds, reduce_io, &
                       modenum, prefix, fildyn, fildvscf, fildrho,   &
                       lnscf, ldisp, nq1, nq2, nq3,                 &
                       eth_rps, eth_ns, lraman, elop, dek, recover,  &
                       fpol, asr, lrpa, lnoloc, start_irr, last_irr, &
                       start_q, last_q, nogg, &
!#ifdef __GWW
                       l_head, omega_gauss, n_gauss, grid_type
!#endif

  ! tr2_ph       : convergence threshold
  ! amass        : atomic masses
  ! alpha_mix    : the mixing parameter
  ! niter_ph     : maximum number of iterations
  ! nmix_ph      : number of previous iterations used in mixing
  ! nat_todo     : number of atom to be displaced
  ! iverbosity   : verbosity control
  ! outdir       : directory where input, output, temporary files reside
  ! epsil        : if true calculate dielectric constant
  ! trans        : if true calculate phonon
  ! elph         : if true calculate electron-phonon coefficients
  ! zue          : if true calculate effective charges (alternate way)
  ! lraman       : if true calculate raman tensor
  ! elop         : if true calculate electro-optic tensor
  ! max_seconds  : maximum cputime for this run
  ! reduce_io    : reduce I/O to the strict minimum
  ! modenum      : single mode calculation
  ! prefix       : the prefix of files produced by pwscf
  ! fildyn       : output file for the dynamical matrix
  ! fildvscf     : output file containing deltavsc
  ! fildrho      : output file containing deltarho
  ! eth_rps      : threshold for calculation of  Pc R |psi> (Raman)
  ! eth_ns       : threshold for non-scf wavefunction calculation (Raman)
  ! dek          : delta_xk used for wavefunctions derivation (Raman)
  ! recover      : recover=.true. to restart from an interrupted run
  ! asr          : in the gamma_gamma case apply acoustic sum rule
  ! start_q      : in q list does the q points from start_q to last_q
  ! last_q       :
  ! start_irr    : does the irred. representation from start_irr to last_irr
  ! last_irr     :
  ! nogg         : if .true. gamma_gamma tricks are not used

  !
  IF ( .NOT. ionode ) GOTO 400
  !
  ! ... Input from file ?
  !
  CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
  READ( 5, '(A)', ERR = 100, IOSTAT = ios ) title
  !
100 CALL errore( 'phq_readin', 'reading title ', ABS( ios ) )
  !
  ! ... set default values for variables in namelist
  !
  tr2_ph       = 1.D-12
  eth_rps      = 1.D-9
  eth_ns       = 1.D-12
  amass(:)     = 0.D0
  alpha_mix(:) = 0.D0
  alpha_mix(1) = 0.7D0
  niter_ph     = maxter
  nmix_ph      = 4
  nat_todo     = 0
  modenum      = 0
  iverbosity   = 0
  trans        = .TRUE.
  lrpa         = .FALSE.
  lnoloc       = .FALSE.
  epsil        = .FALSE.
  zue          = .FALSE.
  fpol         = .FALSE.
  elph         = .FALSE.
  lraman       = .FALSE.
  elop         = .FALSE.
  max_seconds  =  1.E+7_DP
  reduce_io    = .FALSE.
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefix       = 'pwscf'
  fildyn       = 'matdyn'
  fildrho      = ' '
  fildvscf     = ' '
  lnscf        = .TRUE.
  ldisp        = .FALSE.
  nogg         = .FALSE.
  nq1          = 0
  nq2          = 0
  nq3          = 0
  dek          = 1.0d-3
  recover      = .FALSE.
  asr          = .FALSE.
  start_irr    = 0
  last_irr     = -1000
  start_q      = 1
  last_q       =-1000
  !
!#ifdef __GWW
  l_head       = .FALSE.
  n_gauss      = 79
  omega_gauss  = 20
  grid_type    = 2
!#endif
  !
  ! ...  reading the namelist inputph
  !
  READ( 5, INPUTPH, ERR = 200, IOSTAT = ios )
  !
200 CALL errore( 'phq_readin', 'reading inputph namelist', ABS( ios ) )
  !
  ! ... Check all namelist variables
  !
  IF (.NOT. lnscf) CALL errore &
                   (' phq_readin', ' lnscf=.false. no longer allowed', 1)
  IF (tr2_ph <= 0.D0) CALL errore (' phq_readin', ' Wrong tr2_ph ', 1)
  IF (eth_rps<= 0.D0) CALL errore ( 'phq_readin', ' Wrong eth_rps', 1)
  IF (eth_ns <= 0.D0) CALL errore ( 'phq_readin', ' Wrong eth_ns ', 1)

  DO iter = 1, maxter
     IF (alpha_mix (iter) .LT.0.D0.OR.alpha_mix (iter) .GT.1.D0) CALL &
          errore ('phq_readin', ' Wrong alpha_mix ', iter)
  ENDDO
  IF (niter_ph.LT.1.OR.niter_ph.GT.maxter) CALL errore ('phq_readin', &
       ' Wrong niter_ph ', 1)
  IF (nmix_ph.LT.1.OR.nmix_ph.GT.5) CALL errore ('phq_readin', ' Wrong &
       &nmix_ph ', 1)
  IF (iverbosity.NE.0.AND.iverbosity.NE.1) CALL errore ('phq_readin', &
       &' Wrong  iverbosity ', 1)
  IF (fildyn.EQ.' ') CALL errore ('phq_readin', ' Wrong fildyn ', 1)
  IF (max_seconds.LT.1.D0) CALL errore ('phq_readin', ' Wrong max_seconds', 1)

  IF (modenum < 0) CALL errore ('phq_readin', ' Wrong modenum ', 1)
  IF (modenum > 0 .AND. .NOT. lnscf) CALL errore ('phq_readin', &
       &' non-scf calculation must precede single-mode calculation', 1)
  IF (dek <= 0.d0) CALL errore ( 'phq_readin', ' Wrong dek ', 1)
  epsil = epsil .OR. lraman .OR. elop
  IF ( (lraman.OR.elop) .AND. fildrho == ' ') fildrho = 'drho'
  !
  !    reads the q point (just if ldisp = .false.)
  !
  IF (.NOT. ldisp) THEN
     READ (5, *, err = 300, iostat = ios) (xq (ipol), ipol = 1, 3)
  END IF
300 CALL errore ('phq_readin', 'reading xq', ABS (ios) )
  lgamma = xq (1) .EQ.0.D0.AND.xq (2) .EQ.0.D0.AND.xq (3) .EQ.0.D0
  IF ( (epsil.OR.zue) .AND..NOT.lgamma) CALL errore ('phq_readin', &
       'gamma is needed for elec.field', 1)
  IF (zue.AND..NOT.trans) CALL errore ('phq_readin', 'trans must be &
       &.t. for Zue calc.', 1)

  IF (trans.AND.(lrpa.OR.lnoloc)) CALL errore('phq_readin', &
                    'only dielectric constant with lrpa or lnoloc',1)
  !
  ! reads the frequencies ( just if fpol = .true. )
  !
  IF ( fpol ) THEN
     READ (5, *, err = 10, iostat = ios) card
     IF ( TRIM(card) == 'FREQUENCIES' ) THEN
        READ (5, *, err = 10, iostat = ios) nfs
        DO i = 1, nfs
           READ (5, *, err = 10, iostat = ios) fiu(i)
        END DO
     END IF
  END IF
10 CALL errore ('phq_readin', 'reading FREQUENCIES card', ABS(ios) )
  !
  tmp_dir = trimcheck (outdir)
  !
400 CONTINUE
  CALL bcast_ph_input ( )
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  !   amass will also be read from file:
  !   save its content in auxiliary variables
  !
  amass_input(:)= amass(:)
  !
  tmp_dir_ph= TRIM (tmp_dir) // '_ph'
  tmp_dir_save=tmp_dir
!
! The recover file is always written with the _ph prefix
!
  IF (recover) THEN
     CALL ph_readfile('init',ierr)
     IF (ierr /= 0 ) THEN
        recover=.FALSE.
        GOTO 1001
     ENDIF
     tmp_dir=tmp_dir_ph
     CALL check_restart_recover(exst_recover, exst_restart)
     tmp_dir=tmp_dir_save
     IF (ldisp) lgamma = (current_iq==1)
!
!  In this case ph.x has already saved its own data-file and we read
!  initial information from that file
!
     IF ((exst_recover.OR.exst_restart).AND..NOT.lgamma) tmp_dir=tmp_dir_ph
  ENDIF
1001 CONTINUE

  CALL read_file ( )
  !
  IF (modenum > 3*nat) CALL errore ('phq_readin', ' Wrong modenum ', 2)

  IF (gamma_only) CALL errore('phq_readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)

  IF (lda_plus_u) CALL errore('phq_readin',&
     'The phonon code with LDA+U is not yet available',1)

  IF (okpaw.and.(lraman.or.elop.or.elph)) CALL errore('phq_readin',&
     'The phonon code with paw and raman, elop or elph is not yet available',1)

  IF (nproc /= nproc_file .and. .not. twfcollect)  &
     CALL errore('phq_readin',&
     'pw.x run with a different number of processors. Use wf_collect=.true.',1)

  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('phq_readin',&
     'pw.x run with a different number of pools. Use wf_collect=.true.',1)

  IF (two_fermi_energies.or.i_cons /= 0) &
     CALL errore('phq_readin',&
     'The phonon code with constrained magnetization is not yet available',1)

  IF (tqr) CALL errore('phq_readin',&
     'The phonon code with Q in real space not available',1)

  IF (start_irr < 0 ) CALL errore('phq_readin', 'wrong start_irr',1)
  !
  IF (start_q <= 0 ) CALL errore('phq_readin', 'wrong start_q',1)
  !
  IF (noncolin.and.fpol) &
              CALL errore('phonon','noncolinear and fpol not programed',1)
  !
  ! If a band structure calculation needs to be done do not open a file
  ! for k point
  !
  IF (lnscf.or.ldisp) lkpoint_dir=.FALSE.
  restart = recover
  !
  !  set masses to values read from input, if available;
  !  leave values read from file otherwise
  !
  DO it = 1, ntyp
     IF (amass_input(it) < 0.0_DP) amass_input(it)= &
              atom_weight(atomic_number(TRIM(atm(it))))
     IF (amass_input(it) > 0.D0) amass(it) = amass_input(it)
     IF (amass(it) <= 0.D0) CALL errore ('phq_readin', 'Wrong masses', it)
  ENDDO
  lgamma_gamma=.FALSE.
  IF (nkstot==1.OR.(nkstot==2.AND.nspin==2)) THEN
     lgamma_gamma=(lgamma.AND.(ABS(xk(1,1))<1.D-12) &
                         .AND.(ABS(xk(2,1))<1.D-12) &
                         .AND.(ABS(xk(3,1))<1.D-12) )
  ENDIF
  IF (nogg) lgamma_gamma=.FALSE.
  !
  IF (lgamma) THEN
     nksq = nks
  ELSE
     nksq = nks / 2
  ENDIF
  !
  IF (tfixed_occ) &
     CALL errore('phq_readin','phonon with arbitrary occupations not tested',1)
  !
  IF (elph.AND..NOT.lgauss) CALL errore ('phq_readin', 'Electron-&
       &phonon only for metals', 1)
  IF (elph.AND.lsda) CALL errore ('phq_readin', 'El-ph and spin not &
       &implemented', 1)
!  IF (elph.AND.ldisp.AND..NOT.trans) CALL errore ('phq_readin', &
!       'El-ph on a grid requires phonon calculation', 1)
  IF (elph.AND.fildvscf.EQ.' ') CALL errore ('phq_readin', 'El-ph needs &
       &a DeltaVscf file', 1)
  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  CALL allocate_part ( nat )
  !
  IF ( .NOT. ionode ) GOTO 800

  IF (nat_todo.LT.0.OR.nat_todo.GT.nat) CALL errore ('phq_readin', &
       'nat_todo is wrong', 1)
  IF (nat_todo.NE.0) THEN
     READ (5, *, err = 600, iostat = ios) (atomo (na), na = 1, &
          nat_todo)
600  CALL errore ('phq_readin', 'reading atoms', ABS (ios) )
  ENDIF

800 CONTINUE

  CALL bcast_ph_input1 ( )

  IF (epsil.AND.lgauss) &
        CALL errore ('phq_readin', 'no elec. field with metals', 1)
  IF (MOD (nkstot, 2) .NE.0.AND..NOT.lgamma.and..not.lnscf) &
           CALL errore ('phq_readin', 'k-points are odd', nkstot)
  IF (modenum > 0) THEN
     IF ( ldisp ) &
          CALL errore('phq_readin','Dispersion calculation and &
          & single mode calculation not possibile !',1)
  ENDIF
  IF ((nat_todo /= 0 ) .and. lgamma_gamma) CALL errore( &
     'phq_readin', 'gamma_gamma tricks with nat_todo &
      & not available. Use nogg=.true.', 1)

  IF (modenum > 0 .OR. ldisp .OR. lraman ) lgamma_gamma=.FALSE.
  IF (.not.lgamma_gamma) asr=.FALSE.
  !
  IF (ldisp .AND. (nq1 .LE. 0 .OR. nq2 .LE. 0 .OR. nq3 .LE. 0)) &
       CALL errore('phq_readin','nq1, nq2, and nq3 must be greater than 0',1)
  !
  RETURN
  !
END SUBROUTINE phq_readin
