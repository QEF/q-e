!
! Copyright (C) 2002-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE iosys()
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine reads input data from standard input ( unit 5 )
  ! ...  ---------------------------------------------------------------
  !
  ! ...  access the modules renaming the variables that have the same name 
  ! ...  as the input parameters, this is required in order to use a code
  ! ...  independent input parser
  !
  !
  USE constants,     ONLY : AU, eV_to_kelvin
  !
  USE io_global,     ONLY : stdout
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
  USE control_flags, ONLY : diis_ethr_cg, diis_ndim, diis_wfc_keep, isolve, &
                            max_cg_iter, diis_buff, david, imix, nmix, &
                            iverbosity, tr2, niter, order, iswitch, &
                            upscale_     => upscale, &
                            mixing_beta_ => mixing_beta, &
                            nstep_       => nstep, &
                            iprint_      => iprint, &
                            nosym_       => nosym, &
                            modenum_     => modenum, &
                            reduce_io, ethr, lscf, lbfgs, lmd, lneb, lphonon, &
                            noinv, time_max, restart, loldbfgs
  USE wvfct,         ONLY : ibm_baco2, &
                            nbnd_ => nbnd
  USE fixed_occ,     ONLY : tfixed_occ
  USE control_flags, ONLY : twfcollect 
  USE neb_variables, ONLY : lsteep_des, lquick_min , ldamped_dyn, lmol_dyn, &
                            num_of_images_  => num_of_images, &
                            CI_scheme_      => CI_scheme, &
                            VEC_scheme_     => VEC_scheme, &
                            optimization_   => optimization, &
                            damp_           => damp, &
                            temp_req_       => temp_req, &
                            ds_             => ds, &
                            k_max_          => k_max, & 
                            k_min_          => k_min, &
                            neb_thr_        => neb_thr
  USE noncollin_module, ONLY : baco_ibm_xlf, &
                               noncolin_  => noncolin, &
                               lambda_    => lambda, &
                               i_cons_    => i_cons, &
                               mcons_     => mcons, &
                               angle1_    => angle1, &
                               angle2_    => angle2, &
                               report_    => report
  USE bfgs_module,   ONLY : bfgs_xlf_bug, &
                            lbfgs_ndim_       => lbfgs_ndim, &
                            trust_radius_max_ => trust_radius_max, &
                            trust_radius_min_ => trust_radius_min, &
                            trust_radius_ini_ => trust_radius_ini, &
                            trust_radius_end_ => trust_radius_end, &
                            w_1_              => w_1, & 
                            w_2_              => w_2 
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
                               startingpot, conv_thr, diago_thr_init
  !
  ! IONS namelist
  !
  USE input_parameters, ONLY : ion_dynamics, ion_positions, ion_temperature, &
                               tempw, tolp, upscale, potential_extrapolation, &
                               CI_scheme, VEC_scheme, minimization_scheme, &
                               num_of_images, optimization, damp, temp_req, &
                               ds, k_max, k_min, neb_thr, &
                               trust_radius_max, trust_radius_min, &
                               trust_radius_ini, trust_radius_end, &
                               w_1, w_2, lbfgs_ndim
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
  ! ... NEB specific
  !
  USE input_parameters, ONLY : pos 
  !
  USE read_namelists_module, ONLY: read_namelists
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER             :: unit = 5, &! standard input unit
                         i, iiarg, nargs, ia, ios, ierr, ilen, is, image
  INTEGER             :: iargc
  EXTERNAL               iargc                 
  CHARACTER (LEN=80)  :: input_file
  !
  !
  CALL getenv( 'HOME', pseudo_dir )
  !
  pseudo_dir = TRIM( pseudo_dir ) // '/pw/pseudo/'
  !
  ! ... Input from file ?
  !
  nargs = iargc() 
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL getarg( iiarg, input_file )  
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )  
        OPEN ( UNIT = unit, FILE = input_file, FORM = 'FORMATTED', &
               STATUS = 'OLD', IOSTAT = ierr )
        CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                   & ' not found' , ierr )
        !
     END IF
     !
  END DO
  !
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
     WRITE( stdout, &
            '(5x,"Presently no symmetry can be used with electric field",/)' )
  END IF
  IF ( tefield .AND. tstress ) THEN
     tstress = .FALSE.
     WRITE( stdout, &
            '(5x,"Presently stress not available with electric field",/)' )
  END IF
  IF ( tefield .AND. ( nspin == 2 ) ) THEN
     CALL errore( 'input', 'LSDA not available with electric field' , 1 )
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
        IF ( degauss /= 0.D0 ) THEN
           CALL errore( ' iosys ', &
                      & ' fixed occupations, gauss. broadening ignored', -1 )
           degauss = 0.D0
        END IF
     CASE ( 'smearing' )
        ltetra = .FALSE.
        IF ( degauss == 0.D0 ) THEN
           CALL errore( ' iosys ', &
                      & ' smearing requires gaussian broadening', 1 )
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
        CALL errore( ' iosys ',' occupations ' // TRIM( occupations ) // &
                   & 'not implemented', 1 )
  END SELECT
  !
  IF( nbnd < 1 ) THEN
     CALL errore( ' iosys ', ' nbnd less than 1 ', nbnd )
  END IF
  IF( nelec < 0 ) THEN
     CALL errore( ' iosys ', ' nelec less than 0 ', nelec )
  END IF
  !
  lsda = ( nspin == 2 )
  !
  ! ... starting_magnetization(ia) = -2.D0 means "not set" -- set it to 0
  !
  DO ia = 1, ntyp
     IF ( starting_magnetization(ia) == -2.D0 ) &
        starting_magnetization(ia) = 0.d0 
  END DO  
  !
  IF ( ecutrho <= 0.D0 ) THEN
     dual = 4.D0
  ELSE
     dual = ecutrho / ecutwfc
     IF( dual <= 1.D0 ) THEN
        CALL errore( ' iosys ', ' invalid dual? ', 1 )
     END IF
  END IF
  !
  SELECT CASE ( TRIM( restart_mode ) )
  CASE ( 'from_scratch' )
     restart        = .FALSE.
     restart_bfgs   = .FALSE.
     startingconfig = 'input'
  CASE ( 'restart' )
     !
     ! ... NEB specific
     !       
     IF ( calculation == 'neb' ) THEN
        restart      = .FALSE.
        restart_bfgs = .FALSE.     
     ELSE
        restart      = .TRUE.
        restart_bfgs = .TRUE.
        startingpot  = 'file'
        startingwfc  = 'file'
        IF ( TRIM( ion_positions ) == 'from_input' ) THEN
           startingconfig = 'input'
        ELSE
           startingconfig = 'file'
        END IF
     END IF
  CASE DEFAULT
     CALL errore( ' iosys ', ' unknown restart_mode ' // &
                & TRIM( restart_mode ), 1 )
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
  ethr = diago_thr_init
  !
  SELECT CASE ( TRIM( calculation ) )
  CASE ( 'scf' )
     lscf      = .TRUE.
     lbfgs     = .FALSE.
     lmd       = .FALSE.
     lneb      = .FALSE.     
     iswitch   = 0 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
     lforce    = tprnfor
     lmovecell = .FALSE.
     lphonon   = .FALSE.
     nstep     = 1
 CASE ( 'nscf' )
     lscf      = .FALSE.
     lbfgs     = .FALSE.
     lmd       = .FALSE.
     lneb      = .FALSE.     
     iswitch   = -1
     lforce    = .FALSE.
     lmovecell = .FALSE.
     lphonon   = .FALSE.
     nstep     = 1
  CASE ( 'relax' )
     lscf      = .TRUE.
     lbfgs     = .TRUE.
     lmd       = .FALSE.
     lneb      = .FALSE.          
     iswitch   = 1 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
     lforce    = .TRUE.
     lmovecell = .FALSE.
     lphonon   = .FALSE.
  CASE ( 'md' )
     lscf      = .TRUE.
     lbfgs     = .FALSE.
     lmd       = .TRUE.
     lneb      = .FALSE.               
     iswitch   = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
     lforce    = .TRUE.
     lmovecell = .FALSE.
     lphonon   = .FALSE.
  CASE ( 'vc-relax' , 'vc-md' )
     lscf      = .TRUE.
     lbfgs     = .FALSE.
     lmd       = .TRUE.
     lneb      = .FALSE.     
     iswitch   = 3
     lmovecell = .TRUE.
     lforce    = .TRUE.
     lphonon   = .FALSE.
 CASE ( 'phonon' )
     lscf      = .FALSE.
     lbfgs     = .FALSE.
     lmd       = .FALSE.
     lneb      = .FALSE.
     iswitch   = -2 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
     lforce    = .FALSE.
     lmovecell = .FALSE.
     lphonon   = .TRUE.
     nstep     = 1
  !
  ! ... NEB specific
  !   
  CASE ( 'neb' )
     lscf      = .TRUE.
     lbfgs     = .FALSE.
     lmd       = .FALSE.     
     lneb      = .TRUE.
     iswitch   = 1 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
     lforce    = tprnfor
     lmovecell = .FALSE.
     lphonon   = .FALSE.
  CASE DEFAULT
     CALL errore( ' iosys ', ' calculation ' // &
                & TRIM( calculation ) // ' not implemented', 1 )
  END SELECT
  !
  IF ( modenum /= 0 ) THEN
     iswitch = - 4
  END IF
  !
  IF ( startingpot /= 'atomic' .AND. startingpot /= 'file' ) THEN
     CALL errore( ' iosys', 'wrong startingpot: use default', -1 )
     IF (       lscf ) startingpot = 'atomic'
     IF ( .NOT. lscf ) startingpot = 'file'
  END IF
  !
  IF ( .NOT. lscf .AND. startingpot /= 'file' ) THEN
     CALL errore( ' iosys', 'wrong startingpot: use default', -1 )
     startingpot = 'file'
  END IF
  !
  IF ( startingwfc /= 'atomic' .AND. startingwfc /= 'random' .AND. &
       startingwfc /= 'file' ) THEN
     CALL errore( ' iosys', 'wrong startingwfc: use default', -1 )
     startingwfc = 'atomic'
  END IF
  !
  SELECT CASE ( TRIM( electron_dynamics ) )
  CASE ( 'none' )
     CONTINUE
  CASE DEFAULT
     CALL errore( ' iosys ',' unknown electron_dynamics '// &
                & TRIM( electron_dynamics ), 1 )
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
     CALL errore( ' iosys ', &
                & ' fixed occupations and lsda not implemented ', 1 )
  END IF
  !
  calc = ' '
  !
  IF ( TRIM( calculation ) == 'relax' ) THEN
     SELECT CASE ( TRIM( ion_dynamics ) )
     CASE ( 'bfgs' )
        iswitch = 1 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        lbfgs    = .TRUE.
        loldbfgs = .FALSE.
        epse    = etot_conv_thr
        epsf    = forc_conv_thr
        !
        IF ( epse <= 20.D0 * ( tr2 / upscale ) ) &
           CALL errore( ' iosys ', ' required etot_conv_thr is too small:' // &
                      & ' conv_thr must be reduced', 1 )   
        !
     CASE ( 'old-bfgs' )
        iswitch = 1 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        lbfgs    = .FALSE.
        loldbfgs = .TRUE.
        epse    = etot_conv_thr
        epsf    = forc_conv_thr
        !
        IF ( epse <= 20.D0 * ( tr2 / upscale ) ) &
           CALL errore( ' iosys ', ' required etot_conv_thr is too small:' // &
                      & ' conv_thr must be reduced', 1 )   
        !
     CASE ( 'constrained-bfgs' )
        iswitch = 2 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        lbfgs   = .TRUE.
        epse    = etot_conv_thr
        epsf    = forc_conv_thr
     CASE ( 'damp' )
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        lbfgs   = .FALSE.
        lmd     = .TRUE.
        calc    = 'mm'
        epse    = etot_conv_thr
        epsf    = forc_conv_thr
        ntcheck = nstep + 1
     CASE DEFAULT
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     END SELECT
  END IF
  !
  IF ( TRIM( calculation ) == 'md' ) THEN
     SELECT CASE ( TRIM( ion_dynamics ) )
     CASE ( 'verlet' )
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
     CASE ( 'constrained-verlet' )
        iswitch = 4 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
     CASE ( 'beeman' )
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        calc = 'md'
        ntcheck = nstep + 1
     CASE DEFAULT
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     END SELECT
  END IF
  !
  IF ( TRIM( calculation ) == 'vc-relax' ) THEN
     SELECT CASE ( TRIM( cell_dynamics ) )
     CASE ( 'none' )
        epse    = etot_conv_thr
        epsf    = forc_conv_thr
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        calc    = 'mm'
        ntcheck = nstep + 1
     CASE ( 'damp-pr' )
        epse    = etot_conv_thr
        epsf    = forc_conv_thr
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        calc    = 'cm'
        ntcheck = nstep + 1
     CASE ( 'damp-w' )
        epse    = etot_conv_thr
        epsf    = forc_conv_thr
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        calc    = 'nm'
        ntcheck = nstep + 1
     CASE DEFAULT
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': cell_dynamics=' // TRIM( cell_dynamics ) // &
                   & ' not supported', 1 )
     END SELECT
     IF ( TRIM( ion_dynamics ) /= 'damp' ) THEN
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     END IF
  END IF
  !
  IF ( TRIM( calculation ) == 'vc-md' ) THEN
     SELECT CASE ( TRIM( cell_dynamics ) )
     CASE ( 'none' )
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        calc    = 'md'
        ntcheck = nstep + 1
     CASE ( 'pr' )
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        calc    = 'cd'
        ntcheck = nstep + 1
     CASE ( 'w' )
        iswitch = 3 ! ... obsolescent: do not use in new code ( 29/10/2003 C.S.)
        calc    = 'nd'
        ntcheck = nstep + 1
     CASE DEFAULT
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     END SELECT
     IF ( TRIM( ion_dynamics ) /= 'beeman' ) THEN
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     END IF
  END IF
  !
  ! ... NEB specific
  !
  IF ( TRIM( calculation ) == 'neb' ) THEN
     !
     IF ( num_of_images < 2 ) THEN
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': num_of_images must be at least 2', 1 )
     END IF
     !
     IF ( ( CI_scheme /= "no-CI"      ) .AND. &
          ( CI_scheme /= "highest-TS" ) .AND. &
          ( CI_scheme /= "all-SP"     ) .AND. &
          ( CI_scheme /= "manual"     ) ) THEN
        !
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': unknown CI_scheme', 1 )  
        !   
     END IF
     IF ( ( VEC_scheme /= "energy-weighted" )   .AND. &
          ( VEC_scheme /= "gradient-weighted" ) ) THEN
        !
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': unknown VEC_scheme', 1 )
        !
     END IF
     !
     SELECT CASE ( minimization_scheme )
     !
     CASE ( "sd" )
        lsteep_des  = .TRUE.
        lquick_min  = .FALSE.
        ldamped_dyn = .FALSE.
        lmol_dyn    = .FALSE.
     CASE ( "quick-min" )
        lsteep_des  = .FALSE.
        lquick_min  = .TRUE.
        ldamped_dyn = .FALSE.
        lmol_dyn    = .FALSE.
     CASE ( "damped-dyn" )
        lsteep_des  = .FALSE.
        lquick_min  = .FALSE.
        ldamped_dyn = .TRUE.
        lmol_dyn    = .FALSE.
     CASE ( "mol-dyn" )
        lsteep_des  = .FALSE.
        lquick_min  = .FALSE.
        ldamped_dyn = .FALSE.
        lmol_dyn    = .TRUE.
        IF ( temp_req == 0 ) &
           WRITE( stdout,'(/,T2,"WARNING: tepm_req has not been set" )')
        !
        temp_req = temp_req / ( eV_to_kelvin * AU )
        !    
     CASE default
        !
        CALL errore( ' iosys ','calculation=' // TRIM( calculation ) // &
                   & ': unknown minimization_scheme', 1 )  
        !
     END SELECT             
     !
  END IF
  !
  SELECT CASE ( TRIM( ion_temperature ) )
  CASE ( 'not_controlled' )
     CONTINUE
  CASE ( 'rescaling' )
     temperature = tempw
     ttol        = tolp
  CASE DEFAULT
     CALL errore( ' iosys ', ' unknown ion_temperature ' // &
                & TRIM( ion_temperature ), 1 )
  END SELECT
  !
  SELECT CASE ( TRIM( cell_temperature ) )
  CASE ( 'not_controlled' )
     CONTINUE
  CASE DEFAULT
     CALL errore( ' iosys ', ' unknown cell_temperature ' // &
                & TRIM( cell_temperature ), 1 )
  END SELECT
  !
  SELECT CASE ( TRIM( cell_dofree ) )
  CASE ( 'all' )
     CONTINUE
  CASE DEFAULT
     CALL errore( ' iosys ', &
                & ' unknown cell_dofree ' // TRIM( cell_dofree ), 1 )
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
     starting_scf_threshold = SQRT( tr2 )
  CASE DEFAULT
     CALL errore( ' iosys ', ' unknown mixing ' // TRIM( mixing_mode ), 1 )
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
  !
  noncolin_ = noncolin
  angle1_   = angle1
  angle2_   = angle2
  report_   = report
  i_cons_   = i_cons
  mcons_    = mcons
  lambda_   = lambda
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
  ! ... NEB specific
  !
  num_of_images_ = num_of_images
  CI_scheme_     = CI_scheme
  VEC_scheme_    = VEC_scheme
  optimization_  = optimization
  damp_          = damp
  temp_req_      = temp_req
  ds_            = ds
  k_max_         = k_max 
  k_min_         = k_min
  neb_thr_       = neb_thr
  !
  ! ... new BFGS specific
  !
  lbfgs_ndim_       = lbfgs_ndim
  trust_radius_max_ = trust_radius_max
  trust_radius_min_ = trust_radius_min
  trust_radius_ini_ = trust_radius_ini
  trust_radius_end_ = trust_radius_end
  w_1_              = w_1
  w_2_              = w_2  
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
  ! ... set up atomic positions and crystal lattice
  !
  IF ( celldm_(1) == 0.D0 .AND. a /= 0.D0 ) THEN
     IF ( ibrav_ == 0 ) ibrav = 14 
     celldm_(1) = a / bohr_radius_angs
     celldm_(2) = b / a
     celldm_(3) = c / a
     celldm_(4) = cosab
     celldm_(5) = cosac
     celldm_(6) = cosbc 
  ELSE IF ( celldm_(1) /= 0.D0 .AND. a /= 0.D0 ) THEN
     CALL errore( 'input', ' do not specify both celldm and a,b,c!', 1 )
  END IF
  !
  IF ( ibrav_ == 0 .AND. celldm_(1) /= 0.D0 ) THEN
     !
     ! ... input at are in units of alat
     !
     alat = celldm_(1)
  ELSE IF ( ibrav_ == 0 .AND. celldm_(1) == 0.D0 ) THEN
     !
     ! ... input at are in atomic units: define alat
     !
     celldm_(1) = SQRT( at(1,1)**2 + at(1,2)**2 + at(1,3)**2 )
     alat = celldm_(1)
     !
     ! ... bring at to alat units
     !
     at(:,:) = at(:,:) / alat
  ELSE
     !
     ! ... generate at (atomic units)
     !
     CALL latgen( ibrav, celldm_, at(1,1), at(1,2), at(1,3), omega )
     alat = celldm_(1) 
     !
     ! ... bring at to alat units
     !
     at(:,:) = at(:,:) / alat
  END IF
  !
  CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
  !
  IF ( calculation == 'neb' ) THEN
     !
     ! ... NEB specific
     !
     DO image = 1, num_of_images_
        !
        tau = RESHAPE( SOURCE = pos(1:3*nat_,image), SHAPE = (/ 3 , nat_ /) )
        !
        SELECT CASE ( atomic_positions )
           !
           ! ... convert input atomic positions to internally used format:
           ! ... tau in a0 units
           !
        CASE ( 'alat' )
           !
           ! ... input atomic positions are divided by a0: do nothing
           !
        CASE ( 'bohr' )
           !
           ! ... input atomic positions are in a.u.: divide by alat
           !
           tau = tau / alat
        CASE ( 'crystal' )
           !
           ! ... input atomic positions are in crystal axis
           !
           CALL cryst_to_cart( nat_, tau, at, 1 )
        CASE ( 'angstrom' )
           !
           ! ... atomic positions in A: convert to a.u. and divide by alat
           !
           tau = tau / bohr_radius_angs / alat
        CASE DEFAULT
           CALL errore( ' iosys ',' atomic_positions=' // &
                      & TRIM( atomic_positions ) // ' not implemented ', 1 )
        END SELECT
        !
        pos(1:3*nat_,image) = RESHAPE( SOURCE = tau, SHAPE = (/ 3 * nat_ /) )
        !
     END DO 
     !
  ELSE
     !
     SELECT CASE ( atomic_positions )
        !
        ! ... convert input atomic positions to internally used format:
        ! ... tau in a0 units
        !
     CASE ( 'alat' )
        !
        ! ... input atomic positions are divided by a0: do nothing
        !
     CASE ( 'bohr' )
        !
        ! ... input atomic positions are in a.u.: divide by alat
        !
        tau = tau / alat
     CASE ( 'crystal' )
        !
        ! ... input atomic positions are in crystal axis
        !
        CALL cryst_to_cart( nat_, tau, at, 1 )
     CASE ( 'angstrom' )
        !
        ! ... atomic positions in A: convert to a.u. and divide by alat
        !
        tau = tau / bohr_radius_angs / alat
     CASE DEFAULT
        CALL errore( ' iosys ',' atomic_positions=' // &
                   & TRIM( atomic_positions ) // ' not implemented ', 1 )
     END SELECT
     !
  END IF        
  !
  ! ... set default value of wmass
  !
  IF ( wmass == 0.D0 ) THEN
     IF ( calc == 'nd' .OR. calc == 'nm' ) THEN
        DO ia = 1, nat
           wmass = wmass + amass(ityp(ia))
        END DO
        wmass = 0.75D0 * wmass / pi / pi / omega**( 2.D0 / 3.D0 )
     END IF
     IF ( calc == 'cd' .OR. calc == 'cm' ) THEN
        DO ia = 1, nat
           wmass = wmass + amass(ityp(ia))
        END DO
        wmass = 0.75D0 * wmass / pi / pi
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
  CALL readpp()
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
        CALL errore( 'readin', &
                   & 'vcsmd: a positive value for cell mass is required', 1 )
  ELSE
     cell_factor_ = 1.D0
  END IF
  !
  CALL verify_tmpdir()
  !
  CALL restart_from_file()
  !
  IF ( startingconfig == 'file' ) CALL read_config_from_file()
  !
  CALL write_config_to_file()
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
SUBROUTINE read_cards( psfile, atomic_positions_ )
  !-----------------------------------------------------------------------
  !
  USE parser
  USE wvfct,             ONLY : gamma_only
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
  USE input_parameters,  ONLY : atom_label, atom_pfile, atom_mass, &
                                atom_ptyp, taspc, tapos, rd_pos, &
                                atomic_positions, if_pos, sp_pos, &
                                k_points, xk, wk, nk1, nk2, nk3, &
                                k1, k2, k3, nkstot, cell_symmetry, rd_ht, &
                                trd_ht, f_inp, calculation
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
       CALL errore( ' cards ', ' atomic species info missing', 1 )
  IF ( .NOT. tapos ) &
       CALL errore( ' cards ', ' atomic position info missing', 1 )
  !
  DO is = 1, ntyp
    amass(is)  = atom_mass(is)
    psfile(is) = atom_pfile(is)
    atm(is)    = atom_label(is)
    IF( amass(is) <= 0.D0 ) THEN
      CALL errore( ' iosys ', ' invalid  mass ', is )
    END IF
  END DO
  !
  DO ia = 1, nat
     tau(:,ia) = rd_pos(:,ia)
     ityp(ia)  = sp_pos(ia)
  END DO
  !
  ! ... TEMP: calculate fixatom (to be removed)
  !
  fixatom = 0
  fix1: DO ia = nat, 1, -1
    IF ( if_pos(1,ia) /= 0 .OR. &
         if_pos(2,ia) /= 0 .OR. &
         if_pos(3,ia) /= 0 ) EXIT fix1
    fixatom = fixatom + 1
  END DO fix1
  !
  ! ... The constrain on fixed coordinates is implemented using the array 
  ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1 
  ! ... otherwise. fixatom is maintained for compatibility. ( C.S. 15/10/2003 )
  !
  if_pos_ = 1
  if_pos_(:,:) = if_pos(:,1:nat)
  !
  atomic_positions_ = TRIM( atomic_positions )
  !
  IF ( k_points == 'automatic' ) THEN
    !
    ! ... automatic generation of k-points
    !
    gamma_only = .FALSE.
    lxkcry     = .FALSE.
    nks        = 0
    ! nk1,nk2,nk3 and k1,k2,k3 are initialized even when not used
    nk1_       = nk1
    nk2_       = nk2
    nk3_       = nk3
    k1_        = k1
    k2_        = k2
    k3_        = k3
  ELSE IF ( k_points == 'tpiba' ) THEN
    !
    ! ... input k-points are in 2pi/a units
    !
    gamma_only   = .FALSE.
    lxkcry       = .FALSE.
    nks          = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    nk1_       = 0
    nk2_       = 0
    nk3_       = 0
    k1_        = 0
    k2_        = 0
    k3_        = 0
  ELSE IF ( k_points == 'crystal' ) THEN
    !
    ! ... input k-points are in crystal (reciprocal lattice) axis
    !
    gamma_only   = .FALSE.
    lxkcry       = .TRUE.
    nks          = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    nk1_       = 0
    nk2_       = 0
    nk3_       = 0
    k1_        = 0
    k2_        = 0
    k3_        = 0
  ELSE IF ( k_points == 'gamma' ) THEN
    !
    ! ... Only Gamma (k=0) is used
    ! ... specialized routines are used
    !
    gamma_only = .TRUE.
    lxkcry     = .FALSE.
    nks        = 1
    xk_(:,1)   = 0.0
    wk_(1)     = 1.0
    nk1_       = 0
    nk2_       = 0
    nk3_       = 0
    k1_        = 0
    k2_        = 0
    k3_        = 0
  ELSE
    !
    ! ... default: input k-points are in 2pi/a units
    !
    gamma_only   = .FALSE.
    lxkcry       = .FALSE.
    nks          = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    nk1_       = 0
    nk2_       = 0
    nk3_       = 0
    k1_        = 0
    k2_        = 0
    k3_        = 0
  END IF
  !
  IF ( tfixed_occ ) THEN
     IF ( nks > 1 .OR. ( nk1 * nk2 * nk3 ) > 1 ) &
        CALL errore( 'read_cards', &
                   & 'only one k point with fixed occupations', 1 )
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
     CALL errore( ' cards ', ' ibrav=0: must read cell parameters', 1 )
  IF ( ibrav /= 0 .AND. tcell ) &
     CALL errore( ' cards ', ' redundant data for cell parameters', 2 )
  !
  RETURN
  !
END SUBROUTINE read_cards
!
!
!-----------------------------------------------------------------------
SUBROUTINE verify_tmpdir()
  !-----------------------------------------------------------------------
  !
  USE input_parameters, ONLY : restart_mode
  USE control_flags,            ONLY : lneb
  USE io_files,         ONLY : prefix, tmp_dir, nd_nmbr
  USE neb_variables,    ONLY : num_of_images
#if defined (__PARA)
  USE para,             ONLY : me, mypool
  USE mp,               ONLY : mp_barrier
#endif  
  !
  USE parser,           ONLY : int_to_char
  !
  IMPLICIT NONE
  !
  INTEGER            :: l, ios, image
  CHARACTER (LEN=80) :: tmp_dir_saved
  INTEGER            :: c_mkdir
  EXTERNAL              c_mkdir
  !
  !
  ! ... verify if tmp_dir ends with /, add one if needed
  !
  l = LEN_TRIM( tmp_dir )
  IF ( tmp_dir(l:l) /= '/' ) THEN
     IF ( l > 0 .AND. l < LEN( tmp_dir ) ) THEN
        tmp_dir(l+1:l+1)='/'
     ELSE
        CALL errore( 'outdir: ', tmp_dir // ' truncated or empty', 1 )
     END IF
  END IF
  !
  ios = 0
  !
  OPEN( UNIT = 4, FILE = TRIM( tmp_dir ) // 'pwscf' // nd_nmbr, &
        STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  IF ( ios /= 0 ) CALL errore( 'outdir: ', TRIM( tmp_dir ) // &
                             & ' non existent or non writable', 1 )
  !
  ! ... if starting from scratch all temporary files are removed
  ! ... from tmp_dir ( only by the master node )  
  !
  IF ( restart_mode == 'from_scratch' ) THEN
     !
#if defined (__PARA)
     IF ( me == 1 .AND. mypool == 1 ) THEN
#endif        
        !
        ! ... BFGS rstart file is removed       
        !     
        OPEN( UNIT = 4, FILE = TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs', &
              STATUS = 'UNKNOWN' )
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
        !
#if defined (__PARA)
     END IF
#endif 
     !
  END IF    
  !
  ! ... NEB specific
  ! ... in the scratch directory the tree of subdirectories needed by neb are
  ! ... created
  !
  IF ( lneb ) THEN
     !
     tmp_dir_saved = tmp_dir
     ! 
     DO image = 1, num_of_images
        !
        ios = 0
        tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) //"_" // &
                  TRIM( int_to_char( image ) ) // '/'
        !
#if defined (__PARA)
        IF ( me == 1 .AND. mypool == 1 ) THEN
#endif   
           !
           ! ... a scratch directory for this image of the elastic band is
           ! ... created ( only by the master node )
           !
           ios = c_mkdir( TRIM( tmp_dir ), LEN_TRIM( tmp_dir ) )
           !
#if defined (__PARA)
        END IF
        !
        ! ... all jobs are syncronized
        !
        CALL mp_barrier()
#endif          
        !
        ! ... each job checks whether the scratch directory is accessible
        ! ... or not
        ! 
        OPEN( UNIT = 4, FILE = TRIM( tmp_dir ) // 'pwscf' // nd_nmbr, &
              STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
        CLOSE( UNIT = 4, STATUS = 'DELETE' ) 
        ! 
        IF ( ios /= 0 ) &
           CALL errore( 'outdir: ', TRIM( tmp_dir ) // &
                      & ' non existent or non writable', 1 )
        !
        ! ... if starting from scratch all temporary files are removed 
        ! ... from tmp_dir ( only by the master node )
        !
        IF ( restart_mode == 'from_scratch' ) THEN
           !
#if defined (__PARA)
           IF ( me == 1 .AND. mypool == 1 ) THEN
#endif                 
              !
              ! ... standard output of the self consistency is removed
              !      
              OPEN( UNIT = 4, FILE = TRIM( tmp_dir ) // 'PW.out', &
                    STATUS = 'UNKNOWN' )
              CLOSE( UNIT = 4, STATUS = 'DELETE' )
              !
#if defined (__PARA)
           END IF
#endif      
           !
        END IF  
        !
     END DO
     !
     tmp_dir = tmp_dir_saved 
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE verify_tmpdir
