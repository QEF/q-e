!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm
  USE ions_base,     ONLY : amass, atm
  USE check_stop,    ONLY : max_seconds
  USE start_k,       ONLY : reset_grid
  USE klist,         ONLY : xk, nks, nkstot, lgauss, two_fermi_energies, ltetra
  USE control_flags, ONLY : gamma_only, tqr, restart, io_level, &
                            ts_vdw, ldftd3, lxdm, isolve
  USE funct,         ONLY : dft_is_meta, dft_is_hybrid
  USE uspp,          ONLY : okvan
  USE fixed_occ,     ONLY : tfixed_occ
  USE lsda_mod,      ONLY : lsda, nspin
  USE fft_base,      ONLY : dffts
  USE spin_orb,      ONLY : domag
  USE cellmd,        ONLY : lmovecell
  USE run_info,      ONLY : title
  USE control_ph,    ONLY : maxter, alpha_mix, lgamma_gamma, epsil, &
                            zue, zeu, xmldyn, newgrid,                      &
                            trans, reduce_io, tr2_ph, niter_ph,       &
                            nmix_ph, ldisp, recover, lnoloc, start_irr, &
                            last_irr, start_q, last_q, current_iq, tmp_dir_ph, &
                            ext_recover, ext_restart, u_from_file, ldiag, &
                            search_sym, lqdir, electron_phonon, tmp_dir_phq, &
                            rec_code_read, qplot, only_init, only_wfc, &
                            low_directory_check, nk1, nk2, nk3, k1, k2, k3
  USE save_ph,       ONLY : tmp_dir_save, save_ph_input_variables
  USE gamma_gamma,   ONLY : asr
  USE partial,       ONLY : atomo, nat_todo, nat_todo_input
  USE output,        ONLY : fildyn, fildvscf, fildrho
  USE disp,          ONLY : nq1, nq2, nq3, x_q, wq, nqs, lgamma_iq
  USE io_files,      ONLY : tmp_dir, prefix, postfix, create_directory, &
                            check_tempdir, xmlpun_schema
  USE noncollin_module, ONLY : i_cons, noncolin
  USE control_flags, ONLY : iverbosity, modenum
  USE io_global,     ONLY : meta_ionode, meta_ionode_id, ionode, ionode_id, stdout
  USE mp_images,     ONLY : nimage, my_image_id, intra_image_comm,   &
                            me_image, nproc_image
  USE mp_pools,      ONLY : npool
  USE paw_variables, ONLY : okpaw
  USE ramanm,        ONLY : eth_rps, eth_ns, lraman, elop, dek
  USE freq_ph,       ONLY : fpol, fiu, nfs
  USE cryst_ph,      ONLY : magnetic_sym
  USE ph_restart,    ONLY : ph_readfile
  USE el_phon,       ONLY : elph,elph_mat,elph_simple,elph_epa,elph_nbnd_min, elph_nbnd_max, &
                            el_ph_sigma, el_ph_nsigma, el_ph_ngauss,auxdvscf
  USE dfile_star,    ONLY : drho_star, dvscf_star

  USE qpoint,        ONLY : nksq, xq
  USE control_lr,    ONLY : lgamma, lrpa
  ! YAMBO >
  USE YAMBO,         ONLY : elph_yambo,dvscf_yambo
  ! YAMBO <
  USE elph_tetra_mod,ONLY : elph_tetra, lshift_q, in_alpha2f
  USE ktetra,        ONLY : tetra_type
  USE ldaU,          ONLY : lda_plus_u, U_projection, lda_plus_u_kind
  USE ldaU_ph,       ONLY : read_dns_bare, d2ns_type
  USE dvscf_interpolate, ONLY : ldvscf_interpolate, do_long_range, &
      do_charge_neutral, wpot_dir
  USE ahc,           ONLY : elph_ahc, ahc_dir, ahc_nbnd, ahc_nbndskip, &
      skip_upperfan
  USE read_namelists_module, ONLY : check_namelist_read
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios, ipol, iter, na, it, ierr, ierr1
    ! integer variable for I/O control
    ! counter on polarizations
    ! counter on iterations
    ! counter on atoms
    ! counter on types
  REAL(DP), ALLOCATABLE :: amass_input(:)
    ! save masses read from input here
  CHARACTER (LEN=256) :: outdir, filename
  CHARACTER (LEN=8)   :: verbosity
  CHARACTER(LEN=80)          :: card
  CHARACTER(LEN=6) :: int_to_char
  INTEGER                    :: i
  LOGICAL                    :: nogg
  LOGICAL      :: q2d, q_in_band_form
  INTEGER, EXTERNAL  :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  LOGICAL, EXTERNAL  :: imatches
  LOGICAL, EXTERNAL  :: has_xml
  LOGICAL :: exst, parallelfs
  REAL(DP), ALLOCATABLE :: xqaux(:,:)
  INTEGER, ALLOCATABLE :: wqaux(:)
  INTEGER :: nqaux, iq
  CHARACTER(len=80) :: diagonalization='david'
  !
  NAMELIST / INPUTPH / tr2_ph, amass, alpha_mix, niter_ph, nmix_ph,  &
                       nat_todo, verbosity, iverbosity, outdir, epsil,  &
                       trans,  zue, zeu, max_seconds, reduce_io, &
                       modenum, prefix, fildyn, fildvscf, fildrho, &
                       ldisp, nq1, nq2, nq3, &
                       eth_rps, eth_ns, lraman, elop, dek, recover,  &
                       fpol, asr, lrpa, lnoloc, start_irr, last_irr, &
                       start_q, last_q, nogg, ldiag, search_sym, lqdir, &
                       nk1, nk2, nk3, k1, k2, k3, &
                       drho_star, dvscf_star, only_init, only_wfc, &
                       elph_nbnd_min, elph_nbnd_max, el_ph_ngauss, &
                       el_ph_nsigma, el_ph_sigma,  electron_phonon, &
                       q_in_band_form, q2d, qplot, low_directory_check, &
                       lshift_q, read_dns_bare, d2ns_type, diagonalization, &
                       ldvscf_interpolate, do_long_range, do_charge_neutral, &
                       wpot_dir, ahc_dir, ahc_nbnd, ahc_nbndskip, &
                       skip_upperfan

  ! tr2_ph       : convergence threshold
  ! amass        : atomic masses
  ! alpha_mix    : the mixing parameter
  ! niter_ph     : maximum number of iterations
  ! nmix_ph      : number of previous iterations used in mixing
  ! nat_todo     : number of atom to be displaced
  ! verbosity    : verbosity control (iverbosity is obsolete)
  ! outdir       : directory where input, output, temporary files reside
  ! epsil        : if true calculate dielectric constant
  ! trans        : if true calculate phonon
  ! electron-phonon : select the kind of electron-phonon calculation
  ! elph         : if true calculate electron-phonon coefficients
  ! elph_mat     : if true eph coefficients for wannier
  ! zue          : if .true. calculate effective charges ( d force / dE )
  ! zeu          : if .true. calculate effective charges ( d P / du )
  ! lraman       : if true calculate raman tensor
  ! elop         : if true calculate electro-optic tensor
  ! max_seconds  : maximum cputime for this run
  ! reduce_io    : reduce I/O to the strict minimum
  ! modenum      : single mode calculation
  ! prefix       : the prefix of files produced by pwscf
  ! fildyn       : output file for the dynamical matrix
  ! fildvscf     : output file containing deltavsc
  ! fildrho      : output file containing deltarho
  ! fildrho_dir  : directory where fildrho files will be stored (default: outdir or ESPRESSO_FILDRHO_DIR variable)
  ! eth_rps      : threshold for calculation of  Pc R |psi> (Raman)
  ! eth_ns       : threshold for non-scf wavefunction calculation (Raman)
  ! dek          : delta_xk used for wavefunctions derivation (Raman)
  ! recover      : recover=.true. to restart from an interrupted run
  ! asr          : in the gamma_gamma case apply acoustic sum rule
  ! start_q      : in q list does the q points from start_q to last_q
  ! last_q       :
  ! start_irr    : does the irred. representation from start_irr to last_irr
  ! last_irr     :
  ! nogg         : if .true. lgamma_gamma tricks are not used
  ! ldiag        : if .true. force diagonalization of the dyn mat
  ! lqdir        : if .true. each q writes in its own directory
  ! search_sym   : if .true. analyze symmetry if possible
  ! nk1,nk2,nk3, k1,k2,k3 : 
  !              when specified in input, the phonon run uses a different
  !              k-point mesh from that used for the charge density.
  !
  ! dvscf_star%open : if .true. write in dvscf_star%dir the dvscf_q
  !                   'for all q' in the star of q with suffix dvscf_star%ext. 
  !                   The dvscf_q' is written in the basis dvscf_star%basis;
  !                   if dvscf_star%pat is .true. also save a pattern file.
  ! dvscf_star%dir, dvscf_star%ext, dvscf_star%basis : see dvscf_star%open
  ! drho_star%open  : like dvscf_star%open but for drho_q
  ! drho_star%dir, drho_star%ext, drho_star%basis : see drho_star%open
   !
  ! elph_nbnd_min,
  ! elph_nbnd_max: if (elph_mat=.true.) it dumps the eph matrix element from elph_nbnd_min
  !                  to elph_nbnd_max 
  ! el_ph_ngauss,
  ! el_ph_nsigma, 
  ! el_ph_sigma  :  if (elph_mat=.true.) it defines the kind and the val-ue of the
  !                 smearing to be used in the eph coupling calculation.
  ! qplot, : if true a list of q points is given in input
  ! q_in_band_form: if true the input list of q points defines paths 
  ! q2d, : if .true. the q list define a mesh in a square.
  ! low_directory_check : if .true. only the requested representations
  !                       are searched on file
  !
  ! read_dns_bare : If .true. the code tries to read three files in DFPT+U calculations:
  !                 dnsorth, dnsbare, d2nsbare 
  ! d2ns_type     : DFPT+U - the 2nd bare derivative of occupation matrices ns 
  !                 (d2ns_bare matrix). Experimental! This is why it is not documented in Doc.
  !                 d2ns_type='full': matrix calculated with no approximation. 
  !                 d2ns_type='fmmp': assume a m <=> m' symmetry. 
  !                 d2ns_type='diag': if okvan=.true. the matrix is calculated retaining only
  !                                     for <\beta_J|\phi_I> products where for J==I.   
  !                 d2ns_type='dmmp': same as 'diag', but also assuming a m <=> m'.
  ! diagonalization : diagonalization method used in the nscf calc
  ! ldvscf_interpolate: if .true., use Fourier interpolation of phonon potential
  ! do_long_range: if .true., add the long-range part of the potential to dvscf
  ! do_charge_neutral: if .true., impose the neutrality condition on Born effective charges
  ! wpot_dir: folder where w_pot binary files are located
  ! ahc_dir: Directory where the output binary files for AHC e-ph coupling are written
  ! ahc_nbnd: Number of bands where the electron self-energy is to be computed.
  ! ahc_nbndskip: Number of bands to exclude when computing the self-energy.
  ! skip_upperfan: If .true., skip the calculation of upper Fan self-energy.
  !
  ! Note: meta_ionode is a single processor that reads the input
  !       (ionode is also a single processor but per image)
  !       Data read from input is subsequently broadcast to all processors
  !       from ionode_id (using the default communicator world_comm)
  !
  IF (meta_ionode) THEN
  !
  ! ... Input from file ?
  !
     CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
     READ( 5, '(A)', IOSTAT = ios ) title
  !
  ENDIF
  !
  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore( 'phq_readin', 'reading title ', ABS( ios ) )
  CALL mp_bcast(title, meta_ionode_id, world_comm  )
  !
  ! Rewind the input if the title is actually the beginning of inputph namelist
  !
  IF( imatches("&inputph", title) ) THEN
    WRITE(stdout,'(6x,a)') "Title line not specified: using 'default'."
    title='default'
    IF (meta_ionode) REWIND(5, iostat=ios)
    CALL mp_bcast(ios, meta_ionode_id, world_comm  )
    CALL errore('phq_readin', 'Title line missing from input.', abs(ios))
  ENDIF
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
  iverbosity   = 1234567
  verbosity    = 'default'
  trans        = .TRUE.
  lrpa         = .FALSE.
  lnoloc       = .FALSE.
  epsil        = .FALSE.
  zeu          = .TRUE.
  zue          = .FALSE.
  fpol         = .FALSE.
  electron_phonon=' '
  elph_nbnd_min = 1
  elph_nbnd_max = 0
  el_ph_sigma  = 0.02
  el_ph_nsigma = 10
  el_ph_ngauss = 1
  lraman       = .FALSE.
  elop         = .FALSE.
  max_seconds  =  1.E+7_DP
  reduce_io    = .FALSE.
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefix       = 'pwscf'
  fildyn       = 'matdyn'
  fildrho      = ' '
  fildvscf     = ' '
  ldisp        = .FALSE.
  nq1          = 0
  nq2          = 0
  nq3          = 0
  dek          = 1.0d-3
  nogg         = .FALSE.
  recover      = .FALSE.
  asr          = .FALSE.
  start_irr    = 0
  last_irr     = -1000
  start_q      = 1
  last_q       =-1000
  ldiag        =.FALSE.
  lqdir        =.FALSE.
  qplot        =.FALSE.
  q_in_band_form=.FALSE.
  q2d         = .FALSE.
  only_init  = .FALSE.
  only_wfc    = .FALSE.
  low_directory_check=.FALSE.
  search_sym   =.TRUE.
  nk1       = 0
  nk2       = 0
  nk3       = 0
  k1       = 0
  k2       = 0
  k3       = 0
  !
  ! dvscf_interpolate
  ldvscf_interpolate = .FALSE.
  do_charge_neutral = .FALSE.
  do_long_range = .FALSE.
  wpot_dir = ' '
  !
  ! electron_phonon == 'ahc'
  ahc_dir = ' '
  ahc_nbnd = 0
  ahc_nbndskip = 0
  skip_upperfan = .FALSE.
  elph_ahc = .FALSE.
  !
  drho_star%open = .FALSE.
  drho_star%basis = 'modes'
  drho_star%pat  = .TRUE.
  drho_star%ext = 'drho'
  CALL get_environment_variable( 'ESPRESSO_FILDRHO_DIR', drho_star%dir)
  IF ( TRIM( drho_star%dir ) == ' ' ) &
      drho_star%dir = TRIM(outdir)//"/Rotated_DRHO/"
  !
  dvscf_star%open = .FALSE.
  dvscf_star%basis = 'modes'
  dvscf_star%pat  = .FALSE.
  dvscf_star%ext = 'dvscf'
  CALL get_environment_variable( 'ESPRESSO_FILDVSCF_DIR', dvscf_star%dir)
  IF ( TRIM( dvscf_star%dir ) == ' ' ) &
      dvscf_star%dir = TRIM(outdir)//"/Rotated_DVSCF/"
  !
  lshift_q = .false.
  read_dns_bare =.false.
  d2ns_type = 'full'
  !
  ! ...  reading the namelist inputph
  !
  IF (meta_ionode) THEN
     !READ( 5, INPUTPH, ERR=30, IOSTAT = ios )
     READ( 5, INPUTPH, IOSTAT = ios )
     !
     ! ...  iverbosity/verbosity hack
     !
     IF ( iverbosity == 1234567 ) THEN
        SELECT CASE( trim( verbosity ) )
           CASE( 'debug', 'high', 'medium' )
              iverbosity = 1
           CASE( 'low', 'default', 'minimal' )
              iverbosity = 0
           CASE DEFAULT
              iverbosity = 0
         END SELECT
     ELSE
        ios = 1234567
     END IF
  END IF
  CALL check_namelist_read(ios, 5, "inputph")
30 CONTINUE
  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  IF ( ios == 1234567 ) THEN
     CALL infomsg( 'phq_readin' , &
                 'iverbosity is obsolete, use "verbosity" instead' )
  ELSE IF ( ABS(ios) /= 0 ) THEN
     CALL errore( 'phq_readin', 'reading inputph namelist', ABS( ios ) )
  END IF
  ! diagonalization option
  SELECT CASE(TRIM(diagonalization))
  CASE ('david','davidson')
     isolve = 0
  CASE ('cg')
     isolve = 1
  CASE DEFAULT
     CALL errore('phq_readin','diagonalization '//trim(diagonalization)//' not implemented',1)
  END SELECT

  !
  ! ...  broadcast all input variables
  !
  tmp_dir = trimcheck (outdir)
  CALL bcast_ph_input ( )
  CALL mp_bcast(nogg, meta_ionode_id, world_comm  )
  CALL mp_bcast(q2d, meta_ionode_id, world_comm  )
  CALL mp_bcast(q_in_band_form, meta_ionode_id, world_comm  )
  !
  drho_star%dir=trimcheck(drho_star%dir)
  dvscf_star%dir=trimcheck(dvscf_star%dir)
  ! filename for the star must always be automatically generated:
  IF(drho_star%ext(1:5)/='auto:')  drho_star%ext  = 'auto:'//drho_star%ext
  IF(dvscf_star%ext(1:5)/='auto:') dvscf_star%ext = 'auto:'//dvscf_star%ext
  !
  ! ... Check all namelist variables
  !
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
  IF (max_seconds.LT.0.1D0) CALL errore ('phq_readin', ' Wrong max_seconds', 1)

  IF (only_init.AND.only_wfc) CALL errore('phq_readin', &
                        'only_init or only_wfc can be .true.', 1)

  IF (modenum < 0) CALL errore ('phq_readin', ' Wrong modenum ', 1)
  IF (dek <= 0.d0) CALL errore ( 'phq_readin', ' Wrong dek ', 1)
  !
  !
  elph_tetra = 0
  SELECT CASE( trim( electron_phonon ) )
  CASE( 'simple'  )
     elph=.true.
     elph_mat=.false.
     elph_simple=.true. 
     elph_epa=.false.
  CASE( 'epa' )
     elph=.true.
     elph_mat=.false.
     elph_simple=.false.
     elph_epa=.true.
  CASE( 'Wannier' )
     elph=.true.
     elph_mat=.true.
     elph_simple=.false.
     elph_epa=.false.
     auxdvscf=trim(fildvscf)
  CASE( 'interpolated' )
     elph=.true.
     elph_mat=.false.
     elph_simple=.false.
     elph_epa=.false.
  ! YAMBO >
  CASE( 'yambo' )
     elph=.true.
     elph_mat=.false.
     elph_simple=.false.
     elph_epa=.false.
     elph_yambo=.true.
     nogg=.true.
     auxdvscf=trim(fildvscf)
  CASE( 'dvscf' )
     elph=.false.
     elph_mat=.false.
     elph_simple=.false.
     elph_epa=.false.
     elph_yambo=.false.
     dvscf_yambo=.true.
     nogg=.true.
     auxdvscf=trim(fildvscf)
  ! YAMBO <
  CASE( 'lambda_tetra'  )
     elph=.true.
     elph_mat=.false.
     elph_simple=.false.
     trans = .false.
     elph_tetra = 1
  CASE( 'gamma_tetra'  )
     elph=.true.
     elph_mat=.false.
     elph_simple=.false.
     trans = .false.
     elph_tetra = 2
  CASE( 'scdft_input'  )
     elph=.true.
     elph_mat=.false.
     elph_simple=.false.
     trans = .false.
     elph_tetra = 3
  CASE( 'ahc' )
     elph = .true.
     elph_ahc = .true.
     elph_mat = .false.
     elph_simple = .false.
     elph_epa = .false.
     trans = .false.
     nogg = .true.
  CASE DEFAULT
     elph=.false.
     elph_mat=.false.
     elph_simple=.false.
     elph_epa=.false.
  END SELECT

  ! YAMBO >
  IF (.not.elph_yambo.and..not.dvscf_yambo) then
    ! YAMBO <
    IF (qplot.AND..NOT.ldisp) CALL errore('phq_readin','qplot requires ldisp=.true.',1)
    !
  ENDIF

  IF (ldisp.AND.only_init.AND.(.NOT.lqdir)) &
     CALL errore('phq_readin', &
                 'only_init=.TRUE. requires lqdir=.TRUE. or data are lost',1)

  epsil = epsil .OR. lraman .OR. elop

  IF (modenum /= 0) search_sym=.FALSE.

  IF (elph_mat .OR. elph_ahc) THEN
     trans = .FALSE.
  ELSEIF (.NOT. elph) THEN
     trans = trans .OR. ldisp
  ENDIF
  !
  ! Set default value for fildrho and fildvscf if they are required
  IF ( (lraman.OR.elop.OR.drho_star%open) .AND. fildrho == ' ') fildrho = 'drho'
  IF ( (elph_mat.OR.dvscf_star%open) .AND. fildvscf == ' ') fildvscf = 'dvscf'
  !
  !  We can calculate  dielectric, raman or elop tensors and no Born effective
  !  charges dF/dE, but we cannot calculate Born effective charges dF/dE
  !  without epsil.
  !
  IF (zeu) zeu = epsil
  !
  !    reads the q point (just if ldisp = .false.)
  !
  IF (meta_ionode) THEN
     ios = 0
     IF (qplot) THEN
        READ (5, *, iostat = ios) nqaux
     ELSE
        IF (.NOT. ldisp) READ (5, *, iostat = ios) (xq (ipol), ipol = 1, 3)
     ENDIF
  END IF
  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore ('phq_readin', 'reading xq', ABS (ios) )
  IF (qplot) THEN
     CALL mp_bcast(nqaux, meta_ionode_id, world_comm )
     ALLOCATE(xqaux(3,nqaux))
     ALLOCATE(wqaux(nqaux))
     IF (meta_ionode) THEN
        DO iq=1, nqaux
           READ (5, *, iostat = ios) (xqaux (ipol,iq), ipol = 1, 3), wqaux(iq)
        ENDDO
     ENDIF
     CALL mp_bcast(ios, meta_ionode_id, world_comm )
     CALL errore ('phq_readin', 'reading xq', ABS (ios) )
     CALL mp_bcast(xqaux, meta_ionode_id, world_comm )
     CALL mp_bcast(wqaux, meta_ionode_id, world_comm )
  ELSE
     CALL mp_bcast(xq, meta_ionode_id, world_comm  )
  ENDIF
  
  IF (.NOT.ldisp) THEN
     lgamma = xq (1) .EQ.0.D0.AND.xq (2) .EQ.0.D0.AND.xq (3) .EQ.0.D0
     IF ( (epsil.OR.zue.or.lraman.or.elop) .AND..NOT.lgamma) &
                CALL errore ('phq_readin', 'gamma is needed for elec.field', 1)
  ENDIF
  IF (zue.AND..NOT.trans) CALL errore ('phq_readin', 'trans must be &
       &.t. for Zue calc.', 1)

  IF (trans.AND.(lrpa.OR.lnoloc)) CALL errore('phq_readin', &
                    'only dielectric constant with lrpa or lnoloc',1)
  IF (lrpa.or.lnoloc) THEN
     zeu=.FALSE.
     lraman=.FALSE.
     elop = .FALSE.
  ENDIF
  !
  ! dvscf_interpolate
  !
  IF (ldvscf_interpolate) THEN
    !
    IF (wpot_dir == ' ') wpot_dir = TRIM(tmp_dir) // "w_pot/"
    wpot_dir = trimcheck(wpot_dir)
    !
    IF (do_charge_neutral .AND. (.NOT. do_long_range)) THEN
      WRITE(stdout, '(5x,a)') 'charge neutrality for dvscf_interpolate is &
        & meaningful only if do_long_range is true.'
      WRITE(stdout, '(5x,a)') 'Set do_charge_neutral = .false.'
      do_charge_neutral = .FALSE.
    ENDIF
    !
  ENDIF
  !
  IF (trans .AND. ldvscf_interpolate) CALL errore ('phq_readin', &
    'ldvscf_interpolate should be used only when trans = .false.', 1)
  IF (domag .AND. ldvscf_interpolate) CALL errore ('phq_readin', &
    'ldvscf_interpolate and magnetism not implemented', 1)
  IF (okpaw .AND. ldvscf_interpolate) CALL errore ('phq_readin', &
    'PAW and ldvscf_interpolate not tested.', 1)
  !
  ! AHC e-ph coupling
  !
  IF (elph_ahc) THEN
    !
    IF (ahc_nbnd <= 0) CALL errore('phq_readin', &
      'ahc_nbnd must be specified as a positive integer')
    IF (ahc_nbndskip < 0) CALL errore('phq_readin', &
      'ahc_nbndskip cannot be negative')
    !
    IF (ahc_dir == ' ') ahc_dir = TRIM(tmp_dir) // "ahc_dir/"
    ahc_dir = trimcheck(ahc_dir)
    !
  ENDIF ! elph_ahc
  !
  IF (elph_ahc .AND. trans) CALL errore ('phq_readin', &
    'elph_ahc can be used only when trans = .false.', 1)
  !
  ! reads the frequencies ( just if fpol = .true. )
  !
  IF ( fpol ) THEN
     IF ( .NOT. epsil) CALL errore ('phq_readin', &
                                    'fpol=.TRUE. needs epsil=.TRUE.', 1 )
     nfs=0
     IF (meta_ionode) THEN
        READ (5, *, iostat = ios) card
        IF ( TRIM(card)=='FREQUENCIES'.OR. &
             TRIM(card)=='frequencies'.OR. &
             TRIM(card)=='Frequencies') THEN
           READ (5, *, iostat = ios) nfs
        ENDIF
     ENDIF
     CALL mp_bcast(ios, meta_ionode_id, world_comm  )
     CALL errore ('phq_readin', 'reading number of FREQUENCIES', ABS(ios) )
     CALL mp_bcast(nfs, meta_ionode_id, world_comm  )
     if (nfs < 1) call errore('phq_readin','Too few frequencies',1)
     ALLOCATE(fiu(nfs))
     IF (meta_ionode) THEN
        IF ( TRIM(card) == 'FREQUENCIES' .OR. &
             TRIM(card) == 'frequencies' .OR. &
             TRIM(card) == 'Frequencies' ) THEN
           DO i = 1, nfs
              READ (5, *, iostat = ios) fiu(i)
           END DO
        END IF
     END IF
     CALL mp_bcast(ios, meta_ionode_id, world_comm )
     CALL errore ('phq_readin', 'reading FREQUENCIES card', ABS(ios) )
     CALL mp_bcast(fiu, meta_ionode_id, world_comm  )
  ELSE
     nfs=1
     ALLOCATE(fiu(1))
     fiu=0.0_DP
  END IF
  !
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  !   amass will be read again from file:
  !   save its content in auxiliary variables
  !
  ALLOCATE ( amass_input( SIZE(amass) ) )
  amass_input(:)= amass(:)
  !
  tmp_dir_save=tmp_dir
  tmp_dir_ph= trimcheck( TRIM (tmp_dir) // '_ph' // int_to_char(my_image_id) )
  CALL check_tempdir ( tmp_dir_ph, exst, parallelfs )
  tmp_dir_phq=tmp_dir_ph

  ext_restart=.FALSE.
  ext_recover=.FALSE.
  rec_code_read=-1000
  IF (recover) THEN
!
!    With a recover run we read here the mesh of q points, the current iq,
!    and the current frequency
!
     CALL ph_readfile('init', 0, 0, ierr)
     CALL ph_readfile('status_ph', 0, 0, ierr1)
!
!   If some error occured here, we cannot recover the run
!
     IF (ierr /= 0 .OR. ierr1 /= 0 ) THEN
        write(stdout,'(5x,"Run is not recoverable starting from scratch")')
        recover=.FALSE.
        goto 1001
     ENDIF
!
!   We check if the bands and the information on the pw run are in the directory
!   written by the phonon code for the current q point. If the file exists
!   we read from there, otherwise use the information in tmp_dir.
!
     IF (lqdir) THEN
        tmp_dir_phq= trimcheck ( TRIM(tmp_dir_ph) // TRIM(prefix) // &
                                & '.q_' // int_to_char(current_iq) )
        CALL check_restart_recover(ext_recover, ext_restart)
        IF (.NOT.ext_recover.AND..NOT.ext_restart) tmp_dir_phq=tmp_dir_ph
     ENDIF
     !
     filename=TRIM(tmp_dir_phq)//TRIM(prefix)//postfix//xmlpun_schema
     IF (ionode) inquire (file =TRIM(filename), exist = exst)
     !
     CALL mp_bcast( exst, ionode_id, intra_image_comm )
     !
     !  If this file exist we use it to recover the pw.x informations
     !
     IF (exst) tmp_dir=tmp_dir_phq
     u_from_file=.true.
  ENDIF
1001 CONTINUE

  IF (qplot.AND..NOT.recover) THEN
     IF (q2d) THEN
        nqs=wqaux(2)*wqaux(3)
        ALLOCATE(x_q(3,nqs))
        ALLOCATE(wq(nqs))
        CALL generate_k_in_plane(nqaux, xqaux, wqaux, x_q, wq, nqs)
     ELSEIF (q_in_band_form) THEN
        nqs=SUM(wqaux(1:nqaux-1))+1
        DO i=1,nqaux-1
           IF (wqaux(i)==0) nqs=nqs+1
        ENDDO
        ALLOCATE(x_q(3,nqs))
        ALLOCATE(wq(nqs))
        CALL generate_k_along_lines(nqaux, xqaux, wqaux, x_q, wq, nqs)
     ELSE
        nqs=nqaux
        ALLOCATE(x_q(3,nqs))
        ALLOCATE(wq(nqs))
        wq(:)=wqaux(:)
        x_q(:,1:nqs)=xqaux(:,1:nqs)  
     ENDIF
     DEALLOCATE(xqaux)
     DEALLOCATE(wqaux)
     ALLOCATE(lgamma_iq(nqs))
     DO iq=1, nqs
        lgamma_iq(iq)= ( ABS(x_q(1,iq)) .LT. 1.0e-10_dp ) .AND. &
                       ( ABS(x_q(2,iq)) .LT. 1.0e-10_dp ) .AND. &
                       ( ABS(x_q(3,iq)) .LT. 1.0e-10_dp )
     ENDDO
     WRITE(stdout, '(//5x,"Dynamical matrices for q-points given in input")') 
     WRITE(stdout, '(5x,"(",i4,"q-points):")') nqs
     WRITE(stdout, '(5x,"  N         xq(1)         xq(2)         xq(3) " )')
     DO iq = 1, nqs
        WRITE(stdout, '(5x,i3, 3f14.9)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq)
     END DO
  ENDIF
  !
  ! DFPT+U: the occupation matrix ns is read via read_file
  !
  CALL read_file ( )
  !
  magnetic_sym=noncolin .AND. domag
  !
  ! init_start_grid returns .true. if a new k-point grid is set from values
  ! read from input (this happens if nk1*nk2*nk3 > 0; otherwise reset_grid
  ! returns .false., leaves the current values, read in read_file, unchanged)
  !
  newgrid = reset_grid (nk1, nk2, nk3, k1, k2, k3) 
  !
  tmp_dir=tmp_dir_save
  !
  IF (modenum > 3*nat) CALL errore ('phq_readin', ' Wrong modenum ', 2)

  IF (gamma_only) CALL errore('phq_readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)

  IF (lda_plus_u) THEN
     ! 
     WRITE(stdout,'(/5x,a)') "Phonon calculation with DFPT+U; please cite"
     WRITE(stdout,'(5x,a)')  "A. Floris et al., Phys. Rev. B 84, 161102(R) (2011)"
     WRITE(stdout,'(5x,a)')  "A. Floris et al., Phys. Rev. B 101, 064305 (2020)"
     WRITE(stdout,'(5x,a)')  "in publications or presentations arising from this work."
     ! 
     IF (U_projection.NE."atomic") CALL errore("phq_readin", &
          " The phonon code for this U_projection_type is not implemented",1)
     IF (lda_plus_u_kind.NE.0) CALL errore("phq_readin", &
          " The phonon code for this lda_plus_u_kind is not implemented",1)
     IF (elph) CALL errore("phq_readin", &
          " Electron-phonon with Hubbard U is not supported",1)
     IF (lraman) CALL errore("phq_readin", &
          " The phonon code with Raman and Hubbard U is not implemented",1)
     !
  ENDIF
  ! checks
  IF (elph_ahc .AND. domag) CALL errore ('phq_readin', &
    'elph_ahc and magnetism not implemented', 1)
  IF (elph_ahc .AND. okpaw) CALL errore ('phq_readin', &
    'elph_ahc and PAW not tested.', 1)
  IF (elph_ahc .AND. okvan) CALL errore ('phq_readin', &
    'elph_ahc and PAW not tested.', 1)
  IF (elph_ahc .AND. lda_plus_u) CALL errore ('phq_readin', &
    'elph_ahc and lda_plus_u not tested.', 1)

  IF (ts_vdw) CALL errore('phq_readin',&
     'The phonon code with TS-VdW is not yet available',1)
  
  IF (ldftd3) CALL errore('phq_readin',&
     'The phonon code with Grimme''s DFT-D3 is not yet available',1)

  IF ( dft_is_meta() ) CALL errore('phq_readin',&
     'The phonon code with meta-GGA functionals is not yet available',1)

  IF ( dft_is_hybrid() ) CALL errore('phq_readin',&
     'The phonon code with hybrid functionals is not yet available',1)

  IF (okpaw.and.(lraman.or.elop)) CALL errore('phq_readin',&
     'The phonon code with paw and raman or elop is not yet available',1)

  IF (magnetic_sym) THEN 
     
     WRITE(stdout,'(/5x,a)') "Phonon calculation in the non-collinear magnetic case;"
     WRITE(stdout,'(5x,a)')  "please cite A. Urru and A. Dal Corso, Phys. Rev. B 100," 
     WRITE(stdout,'(5x,a)')  "045115 (2019) for the theoretical background."

     !IF (epsil) CALL errore('phq_readin',&
     !     'The calculation of Born effective charges in the non collinear &
     !      magnetic case does not work yet and is temporarily disabled',1)

     IF (okpaw) CALL errore('phq_readin',&
          'The phonon code with paw and domag is not available yet',1)
  ENDIF

  IF (okvan.and.(lraman.or.elop)) CALL errore('phq_readin',&
     'The phonon code with US-PP and raman or elop not yet available',1)

  IF (noncolin.and.(lraman.or.elop)) CALL errore('phq_readin', &
      'lraman, elop, and noncolin not programed',1)

  IF (lmovecell) CALL errore('phq_readin', &
      'The phonon code is not working after vc-relax',1)

  IF (reduce_io) io_level=1

  if(elph_mat.and.fildvscf.eq.' ') call errore('phq_readin',&
       'el-ph with wannier requires fildvscf',1)

  IF(elph_mat.and.npool.ne.1) call errore('phq_readin',&
       'el-ph with wannier : pools not implemented',1)

  IF(elph.and.nimage>1) call errore('phq_readin',&
       'el-ph with images not implemented',1)
  
  IF (elph.OR.fildvscf /= ' ') lqdir=.TRUE.

  IF(dvscf_star%open.and.nimage>1) CALL errore('phq_readin',&
       'dvscf_star with image parallelization is not yet available',1)
  IF(drho_star%open.and.nimage>1) CALL errore('phq_readin',&
       'drho_star with image parallelization is not yet available',1)

  IF (lda_plus_u .AND. read_dns_bare .AND. ldisp) lqdir=.TRUE.

  IF (.NOT.ldisp) lqdir=.FALSE.

  IF (i_cons /= 0) &
     CALL errore('phq_readin',&
     'The phonon code with constrained magnetization is not yet available',1)

  IF (two_fermi_energies .AND. (ltetra .OR. lgauss)) &
     CALL errore('phq_readin',&
     'The phonon code with two fermi energies is not available for metals',1)

  IF (tqr) CALL errore('phq_readin',&
     'The phonon code with Q in real space not available',1)

  IF (start_irr < 0 ) CALL errore('phq_readin', 'wrong start_irr',1)
  !
  IF (start_q <= 0 ) CALL errore('phq_readin', 'wrong start_q',1)
  !
  !  the dynamical matrix is written in xml format if fildyn ends in
  !  .xml or in the noncollinear case.
  !
  xmldyn=has_xml(fildyn)
  !IF (noncolin) xmldyn=.TRUE.
  !
  ! If a band structure calculation needs to be done do not open a file
  ! for k point
  !
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
  DEALLOCATE (amass_input)
  lgamma_gamma=.FALSE.
  IF (.NOT.ldisp) THEN
     IF (nkstot==1.OR.(nkstot==2.AND.nspin==2)) THEN
        lgamma_gamma=(lgamma.AND.(ABS(xk(1,1))<1.D-12) &
                            .AND.(ABS(xk(2,1))<1.D-12) &
                            .AND.(ABS(xk(3,1))<1.D-12) )
     ENDIF
     IF (nogg) lgamma_gamma=.FALSE.
     IF ((nat_todo /= 0) .and. lgamma_gamma) CALL errore( &
        'phq_readin', 'gamma_gamma tricks with nat_todo &
       & not available. Use nogg=.true.', 1)
     !
     IF (nimage > 1 .AND. lgamma_gamma) CALL errore( &
        'phq_readin','gamma_gamma tricks with images not implemented',1)
     IF (lgamma) THEN
        nksq = nks
     ELSE
        nksq = nks / 2
     ENDIF
  ENDIF
  IF (lgamma_gamma.AND.ldiag) CALL errore('phq_readin','incompatible flags',1)
  IF (lgamma_gamma.AND.elph ) CALL errore('phq_readin','lgamma_gamma and electron-phonon are incompatible',1)
  !
  IF (tfixed_occ) &
     CALL errore('phq_readin','phonon with arbitrary occupations not tested',1)
  !
  !YAMBO >
  IF (elph .AND. .NOT.(lgauss .OR. ltetra) &
      .AND. .NOT. (elph_yambo .OR. elph_ahc)) &
          CALL errore ('phq_readin', 'Electron-phonon only for metals', 1)
  !YAMBO <
  IF (elph .AND. fildvscf.EQ.' ' .AND. .NOT. ldvscf_interpolate) &
      CALL errore ('phq_readin', 'El-ph needs a DeltaVscf file', 1)
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  CALL allocate_part ( nat )
  !
  IF ( nat_todo < 0 .OR. nat_todo > nat ) &
     CALL errore ('phq_readin', 'nat_todo is wrong', 1)
  IF (nat_todo.NE.0) THEN
     IF (meta_ionode) READ (5, *, iostat = ios) (atomo (na), na = 1, nat_todo)
     CALL mp_bcast(ios, meta_ionode_id, world_comm  )
     CALL errore ('phq_readin', 'reading atoms', ABS (ios) )
     CALL mp_bcast(atomo, meta_ionode_id, world_comm  )
  ENDIF
  nat_todo_input=nat_todo

  IF (epsil.AND.(lgauss .OR. ltetra)) &
        CALL errore ('phq_readin', 'no elec. field with metals', 1)
  IF (modenum > 0) THEN
     IF ( ldisp ) &
          CALL errore('phq_readin','Dispersion calculation and &
          & single mode calculation not possibile !',1)
     nat_todo = 0
  ENDIF
  !
  ! Tetrahedron method for DFPT and El-Ph
  !
  IF(ltetra .AND. tetra_type == 0) &
  &  CALL errore ('phq_readin', 'DFPT with the Blochl correction is not implemented', 1)
  !
  IF(.NOT. ltetra .AND. elph_tetra /= 0) &
  &  CALL errore ('phq_readin', '"electron_phonon" and "occupation" are inconsistent.', 1)
  !
  IF (modenum > 0 .OR. lraman ) lgamma_gamma=.FALSE.
  IF (.NOT.lgamma_gamma) asr=.FALSE.
  !
  IF ((ldisp.AND..NOT.qplot) .AND. &
                  (nq1 .LE. 0 .OR. nq2 .LE. 0 .OR. nq3 .LE. 0)) &
       CALL errore('phq_readin','nq1, nq2, and nq3 must be greater than 0',1)

  CALL save_ph_input_variables()
  !
  RETURN
  !
END SUBROUTINE phq_readin
