!
! Copyright (C) 2002-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
#include "machine.h"
!
!-----------------------------------------------------------------------
SUBROUTINE iosys
  !-----------------------------------------------------------------------
  !   this subroutine reads input data from standard input (unit 5)
  !     ------------------------------------------------------------------
  !
  !    access the modules renaming the variables that have the same name 
  !    as the input parameters, this is required in order to use a code-
  !    independent input parser
  !
  !
  USE io_global,     ONLY : stdout
  USE constants,     ONLY : AU, eV_to_kelvin
  USE bp,            ONLY : nppstr_ => nppstr, &
                            gdir_   => gdir, &
                            lberry_ => lberry
  USE brilz,         ONLY : at, alat, omega, &
                            celldm_ => celldm, &
                            ibrav_  => ibrav
  USE basis,         ONLY : nat_  => nat, &
                            ntyp_ => ntyp, &
                            ityp, tau, atomic_positions, atm, &
                            startingwfc_ => startingwfc, &
                            startingpot_ => startingpot, &
                            startingconfig
  USE char,          ONLY : title_ => title, &
                            crystal
  USE cellmd,        ONLY : cmass, ttol, omega_old, at_old, ntcheck, &
                            cell_factor_ => cell_factor , &
                            press_       => press, &
                            calc, lmovecell
  USE constants,     ONLY : pi, rytoev, uakbar, amconv, bohr_radius_angs
  USE dynam,         ONLY : dt_ => dt, &
                            temperature, amass, delta_t, nraise
  USE extfield,      ONLY : tefield_  => tefield, &
                            dipfield_ => dipfield, &
                            edir_     => edir, &
                            emaxpos_  => emaxpos, &
                            eopreg_   => eopreg, &
                            eamp_     => eamp, &
                            forcefield
  USE io_files,      ONLY : input_drho, output_drho
  USE force_mod,     ONLY : lforce, lstres, force
  USE gvect,         ONLY : dual, &
                            nr1_     => nr1, &
                            nr2_     => nr2, &
                            nr3_     => nr3,  &
                            ecutwfc_ => ecutwfc, &
                            ecfixed_ => ecfixed, &
                            qcutz_   => qcutz, &
                            q2sigma_ => q2sigma
  USE gsmooth,       ONLY : nr1s_ => nr1s, &
                            nr2s_ => nr2s, &
                            nr3s_ => nr3s
  USE klist,         ONLY : xk, wk, nks, ngauss,&
                            xqq_     => xqq, &
                            degauss_ => degauss, &
                            nelec_   => nelec
  USE ktetra,        ONLY : nk1, nk2, nk3, k1, k2, k3, ltetra
  USE ldaU,          ONLY : Hubbard_U_     => hubbard_u, &
                            Hubbard_alpha_ => hubbard_alpha, &
                            niter_with_fixed_ns, &
                            lda_plus_u_ => lda_plus_u
  USE lsda_mod,      ONLY : nspin_                  => nspin, &
                            starting_magnetization_ => starting_magnetization, &
                            lsda
  USE io_files,      ONLY : tmp_dir, &
                            prefix_     => prefix, &
                            pseudo_dir_ => pseudo_dir, &
                            psfile
  USE relax,         ONLY : if_pos, epsf, starting_scf_threshold, &
                            restart_bfgs, epse
  USE varie,         ONLY : diis_ethr_cg, diis_ndim, diis_wfc_keep, isolve, &
                            max_cg_iter, diis_buff, david, imix, nmix, &
                            iverbosity, tr2, niter, order, iswitch, ntypx, &
                            upscale_     => upscale, &
                            mixing_beta_ => mixing_beta, &
                            nstep_       => nstep, &
                            iprint_      => iprint, &
                            nosym_       => nosym, &
                            modenum_     => modenum, &
                            reduce_io, ethr, lscf, noinv, time_max, restart
  USE wvfct,         ONLY : ibm_baco2, nbnd_ => nbnd
  USE fixed_occ,     ONLY : tfixed_occ
  USE control_flags, ONLY : twfcollect 
  USE noncollin_module,     ONLY : baco_ibm_xlf, &
                            noncolin_  => noncolin, &
                            lambda_    => lambda, &
                            i_cons_    => i_cons, &
                            mcons_     => mcons, &
                            angle1_    => angle1, &
                            angle2_    => angle2, &
                            report_    => report

  !
  ! CONTROL namelist
  !
  USE input_parameters, ONLY : title, calculation, verbosity, &
                               restart_mode, nstep, iprint, tstress, tprnfor, &
                               dt, outdir, prefix, max_seconds, &
                               etot_conv_thr, forc_conv_thr, pseudo_dir, &
                               disk_io, tefield, dipfield, lberry, gdir, &
                               nppstr, wf_collect
  !
  ! SYSTEM namelist
  !
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               nat, ntyp, nbnd, nelec, ecutwfc, ecutrho, &
                               nr1, nr2, nr3, nr1s, nr2s, nr3s, &
                               nosym, starting_magnetization, &
                               occupations, degauss, smearing, &
                               nspin, ecfixed, qcutz, q2sigma, &
                               lda_plus_U, Hubbard_U, Hubbard_alpha, &
                               edir, emaxpos, eopreg, eamp, &
                               noncolin, lambda, i_cons, mcons, angle1, &
                               angle2, report
  !
  ! ELECTRONS namelist
  !
  USE input_parameters, ONLY : electron_maxstep, electron_dynamics, &
                               mixing_mode, mixing_beta, mixing_ndim, &
                               mixing_fixed_ns, diago_cg_maxiter, &
                               diago_david_ndim, diago_diis_buff, &
                               diago_diis_ndim, diago_diis_keep, &
                               diago_diis_ethr, diagonalization, startingwfc, &
                               startingpot, conv_thr
  !
  ! IONS namelist
  !
  USE input_parameters, ONLY : ion_dynamics, ion_positions, ion_temperature, &
                               tempw, tolp, upscale, potential_extrapolation
  !
  ! CELL namelist
  !
  USE input_parameters, ONLY : cell_parameters, cell_dynamics, press, &
                               wmass, cell_temperature, cell_dofree, &
                               cell_factor
  !
  ! PHONON namelist
  !
  USE input_parameters, ONLY : phonon, modenum, xqq
  !
  USE read_namelists_module, ONLY: read_namelists
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: unit = 5, i, ia, ios, ierr, ilen, is
  !
  !
  CALL getenv( 'HOME', pseudo_dir )
  !
  pseudo_dir = TRIM( pseudo_dir )//'/pw/pseudo/'
  !
  CALL read_namelists( 'PW' )
  !
  nraise  = 100
  delta_t = 1.D0
  !
  ! ... translate from input to internals of PWscf, various checks
  !
  time_max = max_seconds
  !
  IF ( tefield .AND. ( .NOT. nosym ) ) THEN
     nosym = .TRUE.
     WRITE( stdout,'(5x,"Presently no symmetry can be used with electric field",/)')
  END IF
  IF ( tefield .AND. ( tstress ) ) THEN
     tstress = .FALSE.
     WRITE( stdout,'(5x,"Presently stress not available with electric field",/)')
  END IF
  IF ( tefield .AND. ( nspin == 2 ) ) THEN
     CALL errore('input','LSDA not available with electric field',1)
  END IF
  !
  twfcollect = wf_collect
  !
  ! ...   Set Values for electron and bands
  !
  tfixed_occ = .FALSE.
  !
  SELECT CASE ( TRIM( occupations ) )
     CASE ( 'fixed' )
        ngauss = 0
        ltetra = .FALSE.
        IF ( degauss /= 0.D0) THEN
           CALL errore(' iosys ',&
                       ' fixed occupations, gauss. broadening ignored',-1 )
           degauss = 0.D0
        END IF
     CASE ( 'smearing' )
        ltetra = .FALSE.
        IF ( degauss == 0.D0 ) THEN
           CALL errore(' iosys ',&
                       ' smearing requires gaussian broadening', 1 )
        END IF
        SELECT CASE ( TRIM( smearing ) )
        CASE ( 'gaussian' , 'gauss' )
           ngauss = 0
        CASE ( 'methfessel-paxton' , 'm-p' , 'mp' )
           ngauss = 1
        CASE ( 'marzari-vanderbilt' , 'cold' , 'm-v' , 'mv' )
           ngauss = -1
        CASE ( 'fermi-dirac' , 'f-d' , 'fd' )
           ngauss = -99
        END SELECT
     CASE ( 'tetrahedra' )
        ngauss = 0
        ltetra = .TRUE.
     CASE ( 'from_input' )
        ngauss     = 0
        ltetra     = .FALSE.
        tfixed_occ = .TRUE.
     CASE DEFAULT
        CALL errore(' iosys ',' occupations '//TRIM( occupations )// &
                    'not implemented', 1 )
  END SELECT
  !
  IF( nbnd < 1 ) THEN
     CALL errore(' iosys ',' nbnd less than 1 ', nbnd )
  END IF
  IF( nelec < 0 ) THEN
     CALL errore(' iosys ',' nelec less than 0 ', nelec )
  END IF
  lsda = ( nspin == 2 )
  !noncolin = ( nspin == 4 )
  WRITE( stdout, '(/,5X,"noncolin = ",L1)' ) noncolin
  !  starting_magnetization(ia) = -2.d0 means "not set" -- set it to 0
  DO ia = 1, ntyp
     IF (starting_magnetization(ia) == -2.d0) &
          starting_magnetization(ia) = 0.d0 
  END DO
  !
  IF ( ecutrho <= 0.D0 ) THEN
     dual = 4.D0
  ELSE
     dual = ecutrho / ecutwfc
     IF( dual <= 1.D0 ) THEN
        CALL errore(' iosys ',' invalid dual? ', 1)
     END IF
  END IF
  !
  SELECT CASE ( TRIM( restart_mode ) )
  CASE ( 'from_scratch' )
     restart        = .FALSE.
     restart_bfgs   = .FALSE.
     startingconfig = 'input'
  CASE ( 'restart' )
     restart      = .TRUE.
     restart_bfgs = .TRUE.
     startingpot  = 'file'
     startingwfc  = 'file'
     IF ( TRIM( ion_positions ) == 'from_input') THEN
        startingconfig = 'input'
     ELSE
        startingconfig = 'file'
     END IF
  CASE DEFAULT
     CALL errore(' iosys ',' unknown restart_mode '//TRIM( restart_mode ), 1 )
  END SELECT
  !
  SELECT CASE ( TRIM( disk_io ) )
  CASE ( 'high' )
     reduce_io = .FALSE.
  CASE DEFAULT
     reduce_io = .TRUE.
     restart   = .FALSE.
  END SELECT
  !
  Hubbard_U(:)    = Hubbard_U(:) / rytoev
  Hubbard_alpha(:)= Hubbard_alpha(:) / rytoev
  !
  ethr = 0.D0
  !
  SELECT CASE ( TRIM( calculation ) )
  CASE ( 'scf' )
     lscf      = .TRUE.
     iswitch   = 0
     lforce    = tprnfor
     lmovecell = .FALSE.
     nstep     = 1
 CASE ( 'nscf' )
     lscf      = .FALSE.
     iswitch   = -1
     lforce    = .FALSE.
     lmovecell = .FALSE.
     nstep     = 1
     ! ethr      = 1.D-6
     ! ... I think ethr should not be more strict than that in a simple band
     ! ... structure calculation but there is still something unsatisfactory 
     ! ... in the Davidson diagonalization convergence. SdG 20/03/2003
     !
  CASE ( 'relax' )
     lscf      = .TRUE.
     iswitch   = 1
     lforce    = .TRUE.
     lmovecell = .FALSE.
  CASE ( 'md' )
     lscf      = .TRUE.
     iswitch   = 3
     lforce    = .TRUE.
     lmovecell = .FALSE.
  CASE ( 'vc-relax' , 'vc-md' )
     lscf      = .TRUE.
     iswitch   = 3
     lmovecell = .TRUE.
     lforce    = .TRUE.
 CASE ( 'phonon' )
     lscf      = .FALSE.
     iswitch   = -2
     lforce    = .FALSE.
     lmovecell = .FALSE.
     nstep     = 1
  CASE DEFAULT
     CALL errore(' iosys ',' calculation '//&
                 TRIM( calculation )//' not implemented',1)
  END SELECT
  !
  IF ( modenum /= 0 ) THEN
     iswitch = -4
  END IF
  !
  IF ( startingpot /= 'atomic' .AND. startingpot /= 'file' ) THEN
     CALL errore(' iosys','wrong startingpot: use default',-1)
     IF (      lscf ) startingpot = 'atomic'
     IF ( .NOT.lscf ) startingpot = 'file'
  END IF
  !
  IF ( .NOT. lscf .AND. startingpot /= 'file' ) THEN
     CALL errore(' iosys','wrong startingpot: use default',-1)
     startingpot = 'file'
  END IF
  !
  IF ( startingwfc /= 'atomic' .AND. startingwfc /= 'random' .AND. &
       startingwfc /= 'file' ) THEN
     CALL errore(' iosys','wrong startingwfc: use default',-1)
     startingwfc = 'atomic'
  END IF
  !
  SELECT CASE ( TRIM( electron_dynamics ) )
  CASE ( 'none' )
     CONTINUE
  CASE DEFAULT
     CALL errore(' iosys ',' unknown electron_dynamics '// &
                 TRIM( electron_dynamics ),1)
  END SELECT
  !
  SELECT CASE ( TRIM( diagonalization ) )
  CASE ( 'cg' )
     isolve = 1
     max_cg_iter = diago_cg_maxiter
  CASE ( 'diis' )
     isolve = 2
     max_cg_iter   = diago_cg_maxiter
     diis_buff     = diago_diis_buff
     diis_ethr_cg  = diago_diis_ethr   ! ... SF
     diis_ndim     = diago_diis_ndim   ! ... SF
     diis_wfc_keep = diago_diis_keep
  CASE ( 'david' )
     isolve = 0
     david = diago_david_ndim
  CASE DEFAULT
     isolve = 0
     david = diago_david_ndim
  END SELECT
  !
  tr2   = conv_thr
  niter = electron_maxstep
  !
  SELECT CASE ( TRIM( potential_extrapolation ) )
  CASE ( 'none' )
     order = 0
  CASE ( 'atomic' )
     order = 1
  CASE ( 'wfc' )
     order = 2
  CASE ( 'wfc2' )
     order = 3
  CASE DEFAULT
     order = 1
  END SELECT
  !
  IF ( occupations == 'fixed' .AND. nspin == 2  .AND. lscf ) THEN
     CALL errore(' iosys ',&
                 ' fixed occupations and lsda not implemented ', 1 )
  END IF
  !
  calc = ' '
  !
  IF ( TRIM( calculation ) == 'relax' ) THEN
     SELECT CASE ( TRIM( ion_dynamics ) )
     CASE ( 'bfgs' )
        iswitch = 1
        epse = etot_conv_thr
        epsf = forc_conv_thr
     CASE ( 'constrained-bfgs' )
        iswitch = 2
        epse = etot_conv_thr
        epsf = forc_conv_thr
     CASE ( 'damp' )
        iswitch = 3
        calc = 'mm'
        epse = etot_conv_thr
        epsf = forc_conv_thr
        ntcheck=nstep+1
     CASE DEFAULT
        CALL errore(' iosys ','calculation='//TRIM( calculation )// &
                  ': ion_dynamics='//TRIM( ion_dynamics )//' not supported', 1 )
     END SELECT
  END IF
  !
  IF ( TRIM(calculation) == 'md' ) THEN
     SELECT CASE ( TRIM(ion_dynamics) )
     CASE ('verlet')
        iswitch = 3
     CASE ('constrained-verlet')
        iswitch = 4
     CASE ('beeman')
        iswitch = 3
        calc = 'md'
        ntcheck=nstep+1
     CASE DEFAULT
        CALL errore(' iosys ','calculation='//TRIM( calculation )// &
                 ': ion_dynamics='//TRIM( ion_dynamics )//' not supported', 1 )
     END SELECT
  END IF
  !
  IF ( TRIM( calculation ) == 'vc-relax' ) THEN
     SELECT CASE ( TRIM( cell_dynamics ) )
     CASE ( 'none' )
        epse = etot_conv_thr
        epsf = forc_conv_thr
        iswitch = 3
        calc = 'mm'
        ntcheck = nstep + 1
     CASE ( 'damp-pr' )
        epse = etot_conv_thr
        epsf = forc_conv_thr
        iswitch = 3
        calc = 'cm'
        ntcheck = nstep + 1
     CASE ( 'damp-w' )
        epse = etot_conv_thr
        epsf = forc_conv_thr
        iswitch = 3
        calc = 'nm'
        ntcheck = nstep + 1
     CASE DEFAULT
        CALL errore(' iosys ','calculation='//TRIM( calculation )// &
               ': cell_dynamics='//TRIM( cell_dynamics )//' not supported', 1 )
     END SELECT
     IF ( TRIM( ion_dynamics ) /= 'damp' ) THEN
        CALL errore(' iosys ','calculation='//TRIM( calculation )// &
                 ': ion_dynamics='//TRIM( ion_dynamics )//' not supported', 1 )
     END IF
  END IF
  !
  IF ( TRIM( calculation ) == 'vc-md' ) THEN
     SELECT CASE ( TRIM( cell_dynamics ) )
     CASE ( 'none')
        iswitch = 3
        calc = 'md'
        ntcheck = nstep + 1
     CASE ( 'pr' )
        iswitch = 3
        calc = 'cd'
        ntcheck = nstep + 1
     CASE ( 'w' )
        iswitch = 3
        calc = 'nd'
        ntcheck = nstep + 1
     CASE DEFAULT
        CALL errore(' iosys ','calculation='//TRIM( calculation )// &
                 ': ion_dynamics='//TRIM( ion_dynamics )//' not supported', 1 )
     END SELECT
     IF ( TRIM( ion_dynamics ) /= 'beeman' ) THEN
        CALL errore(' iosys ','calculation='//TRIM( calculation )// &
                 ': ion_dynamics='//TRIM( ion_dynamics )//' not supported', 1 )
     END IF
  END IF
  !
  SELECT CASE ( TRIM( ion_temperature ) )
  CASE ( 'not_controlled' )
     CONTINUE
  CASE ( 'rescaling' )
     temperature = tempw
     ttol        = tolp
  CASE DEFAULT
     CALL errore(' iosys ',' unknown ion_temperature '//&
          TRIM( ion_temperature ), 1 )
  END SELECT
  !
  SELECT CASE ( TRIM( cell_temperature ) )
  CASE ( 'not_controlled' )
     CONTINUE
  CASE DEFAULT
     CALL errore(' iosys ',' unknown cell_temperature '//&
          TRIM( cell_temperature ), 1 )
  END SELECT
  !
  SELECT CASE ( TRIM( cell_dofree ) )
  CASE ( 'all' )
     CONTINUE
  CASE DEFAULT
     CALL errore(' iosys ',' unknown cell_dofree '//TRIM( cell_dofree ), 1 )
  END SELECT
  !
  SELECT CASE ( TRIM( mixing_mode ) )
  CASE ( 'plain' )
     imix = 0
     starting_scf_threshold = tr2
  CASE ( 'TF' )
     imix = 1
     starting_scf_threshold = tr2
  CASE ( 'local-TF' )
     imix = 2
     starting_scf_threshold = tr2
  CASE ( 'potential' )
     imix = -1
     starting_scf_threshold = sqrt(tr2)
  CASE DEFAULT
     CALL errore(' iosys ',' unknown mixing '//TRIM( mixing_mode ), 1)
  END SELECT
  !
  nmix = mixing_ndim
  niter_with_fixed_ns = mixing_fixed_ns
  !
  SELECT CASE ( TRIM( verbosity ) )
  CASE ( 'high' )
     iverbosity = 1
  CASE DEFAULT
     iverbosity = 0
  END SELECT
  !
  IF ( ( lberry ) .AND. ( iswitch /= -1 ) ) CALL errore(" iosys", &
       "calculation='nscf' is mandatory if lberry=.true.",1)
  tmp_dir = TRIM( outdir )
  lstres = ( tstress .AND. lscf )
  !
  ! ... Copy values from input module to PW internals
  !
  nppstr_     = nppstr
  gdir_       = gdir
  lberry_     = lberry
  title_      = title
  dt_         = dt
  tefield_    = tefield
  dipfield_   = dipfield
  prefix_     = TRIM( prefix )
  pseudo_dir_ = TRIM( pseudo_dir )
  nstep_      = nstep
  iprint_     = iprint
  !
  celldm_  = celldm
  ibrav_   = ibrav
  nat_     = nat 
  ntyp_    = ntyp
  edir_    = edir
  emaxpos_ = emaxpos
  eopreg_  = eopreg
  eamp_    = eamp
  nr1_     = nr1
  nr2_     = nr2
  nr3_     = nr3
  ecutwfc_ = ecutwfc
  ecfixed_ = ecfixed
  qcutz_   = qcutz
  q2sigma_ = q2sigma
  nr1s_    = nr1s
  nr2s_    = nr2s
  nr3s_    = nr3s
  degauss_ = degauss
  nelec_   = nelec
  noncolin_ = noncolin
  angle1_  = angle1
  angle2_  = angle2
  report_  = report
  i_cons_  = i_cons
  mcons_   = mcons
  lambda_  = lambda
  !
  Hubbard_U_( 1 : ntyp )     = hubbard_u( 1 : ntyp )
  Hubbard_alpha_( 1 : ntyp ) = hubbard_alpha( 1 : ntyp )
  lda_plus_u_                = lda_plus_u
  nspin_                     = nspin
  starting_magnetization_    = starting_magnetization
  
  nosym_                     = nosym
  nbnd_                      = nbnd
  !
  startingwfc_ = startingwfc
  startingpot_ = startingpot
  mixing_beta_ = mixing_beta
  !
  upscale_     = upscale
  press_       = press
  cell_factor_ = cell_factor
  modenum_     = modenum
  xqq_         = xqq
  !
  ! ... read following cards
  !
  ALLOCATE( tau( 3, nat_ ) )
  ALLOCATE( ityp( nat_ ) )
  ALLOCATE( force( 3, nat_ ) )  ! ... compatibility with old readin
  ALLOCATE( if_pos( 3, nat_ ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat_ ) )
  !
  CALL read_cards( psfile, atomic_positions )
  !
  ! set up atomic positions and crystal lattice
  !
  IF ( celldm_(1) == 0.D0 .AND. a /= 0.D0 ) THEN
     IF (ibrav_ == 0) ibrav = 14 
     celldm_ (1) = a / bohr_radius_angs
     celldm_ (2) = b / a
     celldm_ (3) = c / a
     celldm_ (4) = cosab
     celldm_ (5) = cosac
     celldm_ (6) = cosbc 
  ELSE IF ( celldm_(1) /= 0.D0 .AND. a /= 0.D0 ) THEN
     CALL errore('input', ' do not specify both celldm and a,b,c!',1)
  END IF
  !
  IF ( ibrav_ == 0 .AND. celldm_(1) /= 0.D0 ) THEN
     ! input at are in units of alat
     alat = celldm_(1)
  ELSE IF ( ibrav_ == 0 .AND. celldm_(1) == 0.D0 ) THEN
     ! input at are in atomic units: define alat
     celldm_(1) = SQRT( at(1,1)**2 + at(1,2)**2 + at(1,3)**2 )
     alat = celldm_(1)
     ! bring at to alat units
     at(:,:) = at(:,:) / alat
  ELSE
     ! generate at (atomic units)
     CALL latgen( ibrav, celldm_, at(1,1), at(1,2), at(1,3), omega )
     alat = celldm_(1) 
     ! bring at to alat units
     at(:,:) = at(:,:) / alat
  END IF
  !
  CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
  !
  SELECT CASE ( atomic_positions )
     !
     !  convert input atomic positions to internally used format:
     !  tau in a0 units
     !
  CASE ( 'alat' )
     !
     !  input atomic positions are divided by a0: do nothing
     !
  CASE ( 'bohr' )
     !
     !  input atomic positions are in a.u.: divide by alat
     !
     tau = tau / alat
  CASE ( 'crystal' )
     !
     !  input atomic positions are in crystal axis
     !
     CALL cryst_to_cart( nat_, tau, at, 1)
  CASE ( 'angstrom' )
     !
     !  atomic positions in A: convert to a.u. and divide by alat
     !
     tau = tau / bohr_radius_angs / alat
  CASE DEFAULT
     CALL errore(' iosys ',' atomic_positions='//TRIM( atomic_positions )// &
          ' not implemented ', 1 )
  END SELECT
  !
  ! ... set default value of wmass
  !
  IF ( wmass == 0.D0 ) THEN
     IF ( calc == 'nd' .OR. calc == 'nm' ) THEN
        DO ia = 1, nat
           wmass = wmass + amass(ityp(ia))
        END DO
        wmass =  0.75D0 * wmass / pi / pi / omega**( 2.D0 / 3.D0 )
     END IF
     IF ( calc == 'cd' .OR. calc == 'cm' ) THEN
        DO ia = 1, nat
           wmass = wmass + amass(ityp(ia))
        END DO
        wmass =  0.75D0 * wmass / pi / pi
     END IF
  END IF
  !
  ! ... unit conversion for cell mass and pressure
  !
  cmass  = wmass * amconv
  press_ = press_ / uakbar
  !
  ! ... read pseudopotentials
  !
  CALL readpp
  !
  ! ... Renata's dynamics uses masses in atomic units
  !
  IF ( calc /= ' ' ) THEN
     amass = amass * amconv
  END IF
  !
  ! ... In the case of variable cell dynamics save old cell variables
  ! ... and initialize a few other variables
  !
  IF ( lmovecell ) THEN
     at_old    = at
     omega_old = omega
     lstres    = .TRUE.
     IF ( cell_factor_ <= 0.D0 ) cell_factor_ = 1.2D0
     IF ( cmass <= 0.D0 ) &
        CALL errore('readin',&
                    'vcsmd: a positive value for cell mass is required',1)
  ELSE
     cell_factor_ = 1.D0
  END IF
  !
  CALL verify_tmpdir
  !
!  WRITE( stdout,'(/5x,"current restart_mode = ",a)') TRIM( restart_mode )
!  WRITE( stdout,'( 5x,"current disk_io mode = ",a)') TRIM( disk_io )
  CALL restart_from_file
  !
  IF ( startingconfig == 'file' ) CALL read_config_from_file
  CALL write_config_to_file
  !
  ! ... Files
  !
  input_drho  = ' '
  output_drho = ' '
  crystal     = ' '
  !
  RETURN
  !
END SUBROUTINE iosys
!
!
!-----------------------------------------------------------------------
SUBROUTINE read_cards ( psfile, atomic_positions_ )
  !-----------------------------------------------------------------------
  !
  USE parser
  USE brilz,             ONLY : at, ibrav, symm_type
  USE basis,             ONLY : nat, ntyp, ityp, tau, atm
  USE klist,             ONLY : nks
  USE ktetra,            ONLY : nk1_   => nk1, &
                                nk2_   => nk2, &
                                nk3_   => nk3, &
                                k1_    => k1,  &
                                k2_    => k2,  &
                                k3_    => k3
  USE klist,             ONLY : lxkcry, &
                                xk_    => xk, &
                                wk_    => wk
  USE fixed_occ,         ONLY : tfixed_occ, &
                                f_inp_ => f_inp
  USE relax,             ONLY : fixatom, &
                                if_pos_ =>  if_pos
  USE dynam,             ONLY : amass
  USE input_parameters,  ONLY : atom_label, atom_pfile, atom_mass, atom_ptyp, &
                                taspc, tapos, rd_pos, atomic_positions, &
                                if_pos, sp_pos, k_points, xk, wk, &
                                nk1, nk2, nk3, k1, k2, k3, nkstot, &
                                cell_symmetry, rd_ht, trd_ht, f_inp
  USE read_cards_module, ONLY : read_cards_base => read_cards
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=80) :: psfile(ntyp)
  CHARACTER (LEN=30) :: atomic_positions_
  !
  LOGICAL :: tcell = .FALSE.
  INTEGER :: i, is, ns, ia, ik
  !
  !
  amass = 0
  !
  CALL read_cards_base( 'PW' )
  !
  IF ( .NOT. taspc ) &
       CALL errore(' cards ',' atomic species info missing', 1 )
  IF ( .NOT. tapos ) &
       CALL errore(' cards ',' atomic position info missing', 1 )
  !
  DO is = 1, ntyp
    amass(is)  = atom_mass(is)
    psfile(is) = atom_pfile(is)
    atm(is)    = atom_label(is)
    IF( amass(is) <= 0.D0 ) THEN
      CALL errore(' iosys ',' invalid  mass ', is)
    END IF
  END DO
  !
  DO ia = 1, nat
    tau(:,ia) = rd_pos(:,ia)
    ityp(ia)  = sp_pos(ia)
  END do
  !
  ! TEMP: calculate fixatom
  !
  fixatom = 0
  fix1: DO ia = nat, 1, -1
    IF ( if_pos(1,ia) /= 0 .or. &
         if_pos(2,ia) /= 0 .or. &
         if_pos(3,ia) /= 0 ) EXIT fix1
    fixatom = fixatom + 1
  END DO fix1
  !
  ! ... The constrain on fixed coordinates is implemented using the array 
  ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1 
  ! ... otherwise. fixatom is manteined for compatibility. ( C.S. 15/10/2003 )
  !
  if_pos_ = 1
  if_pos_(:,:) = if_pos(:,1:nat)
  !
  atomic_positions_ = TRIM( atomic_positions )
  !
  IF ( k_points == 'automatic' ) THEN
    !  automatic generation of k-points
    lxkcry = .FALSE.
    nks  = 0
    nk1_ = nk1
    nk2_ = nk2
    nk3_ = nk3
    k1_  = k1
    k2_  = k2
    k3_  = k3
  ELSE IF ( k_points == 'tpiba' ) THEN
    !  input k-points are in 2pi/a units
    lxkcry = .FALSE.
    nks  = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
  ELSE IF ( k_points == 'crystal' ) THEN
    !  input k-points are in crystal (reciprocal lattice) axis
    lxkcry = .TRUE.
    nks  = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
  ELSE IF ( k_points == 'gamma' ) THEN
    !  Only Gamma (k=0) is used
    lxkcry = .FALSE.
    nks = 1
    xk_(:,1) = 0.0
    wk_(1)   = 1.0
  ELSE
    !  default: input k-points are in 2pi/a units
    lxkcry = .FALSE.
    nks  = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
  END IF
  !
  IF ( tfixed_occ ) THEN
     IF ( nks > 1 .OR. ( nk1 * nk2 * nk3 ) > 1 ) &
        CALL errore('read_cards','only one k point with fixed occupations',1)
     f_inp_ = f_inp
  ENDIF
  !
  IF ( trd_ht ) THEN
    symm_type = cell_symmetry 
    at        = TRANSPOSE( rd_ht )
    tcell     = .TRUE.
  END IF
  !
  IF ( ibrav == 0 .AND. .NOT. tcell ) &
       CALL errore(' cards ',' ibrav=0: must read cell parameters', 1 )
  IF ( ibrav /= 0 .AND. tcell ) &
       CALL errore(' cards ',' redundant data for cell parameters', 2 )
  !
  RETURN
  !
END SUBROUTINE read_cards
!
!
!-----------------------------------------------------------------------
SUBROUTINE verify_tmpdir
  !-----------------------------------------------------------------------
  !
  USE io_files, ONLY: tmp_dir, nd_nmbr
  !
  IMPLICIT NONE
  !
  INTEGER :: l, ios
  !
  !
  ! ... verify if tmp_dir ends with /, add one if needed
  !
  l = LEN_TRIM( tmp_dir )
  IF ( tmp_dir(l:l) /= '/' ) THEN
     IF ( l > 0 .AND. l < LEN( tmp_dir ) ) THEN
        tmp_dir(l+1:l+1)='/'
     ELSE
        CALL errore('reading',tmp_dir//' truncated or empty',1)
     END IF
  END IF
  ios = 0
  OPEN( UNIT = 4, FILE = TRIM( tmp_dir )//'pwscf'//nd_nmbr, &
        STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
  IF ( ios /= 0 ) CALL errore('reading',TRIM( tmp_dir )// &
                              ' non existent or non writable',1)
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  RETURN
  !
END SUBROUTINE verify_tmpdir
