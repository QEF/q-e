!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE iosys()
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine reads input data from standard input ( unit 5 )
  ! ...  Use "-input filename" to read input from file "filename":
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  ---------------------------------------------------------------
  !
  ! ...  access the modules renaming the variables that have the same name 
  ! ...  as the input parameters, this is required in order to use a code
  ! ...  independent input parser
  !
  !
  USE kinds,         ONLY : DP
  USE funct,         ONLY : enforce_input_dft
  USE constants,     ONLY : au, eV_to_kelvin, pi, rytoev, &
                            uakbar, amconv, bohr_radius_angs
  USE mp_global,     ONLY : npool, nproc_pool
  !
  USE io_global,     ONLY : stdout, ionode
  !
  USE bp,            ONLY : nppstr_    => nppstr, &
                            gdir_      => gdir, &
                            lberry_    => lberry, &
                            lelfield_  => lelfield, &
                            efield_    => efield, &
                            nberrycyc_ => nberrycyc

  !
  USE cell_base,     ONLY : at, bg, alat, omega, &
                            celldm_ => celldm, &
                            ibrav_  => ibrav
  !
  USE ions_base,     ONLY : if_pos, &
                            ntyp_ => nsp, &
                            nat_  => nat, &
                            ityp, tau, atm  
  !
  USE basis,         ONLY : atomic_positions, &
                            startingwfc_ => startingwfc, &
                            startingpot_ => startingpot, &
                            startingconfig
  !
  USE char,          ONLY : title_ => title, &
                            crystal
  !
  USE cellmd,        ONLY : cmass, omega_old, at_old, ntcheck, &
                            cell_factor_ => cell_factor , &
                            press_       => press, &
                            calc, lmovecell
  !
  USE dynam,         ONLY : temperature, amass, &
                            dt_      => dt, &
                            delta_t_ => delta_t, &
                            nraise_  => nraise
  !
  USE extfield,      ONLY : tefield_  => tefield, &
                            dipfield_ => dipfield, &
                            edir_     => edir, &
                            emaxpos_  => emaxpos, &
                            eopreg_   => eopreg, &
                            eamp_     => eamp, &
                            forcefield
  !
  USE io_files,      ONLY : input_drho, output_drho, &
                            psfile, tmp_dir, wfc_dir, &
                            prefix_     => prefix, &
                            pseudo_dir_ => pseudo_dir
  
  !
  USE force_mod,     ONLY : lforce, lstres, force
  !
  USE gvect,         ONLY : dual, &
                            nr1_     => nr1, &
                            nr2_     => nr2, &
                            nr3_     => nr3,  &
                            ecutwfc_ => ecutwfc, &
                            ecfixed_ => ecfixed, &
                            qcutz_   => qcutz, &
                            q2sigma_ => q2sigma
  !
  USE gsmooth,       ONLY : nr1s_ => nr1s, &
                            nr2s_ => nr2s, &
                            nr3s_ => nr3s
  !
  USE klist,         ONLY : xk, wk, nks, ngauss, two_fermi_energies, &
                            xqq_     => xqq, &
                            degauss_ => degauss, &
                            nelec_   => nelec, &
                            nelup_   => nelup, &
                            neldw_   => neldw, &
                            b_length_ => b_length, &
                            lcart_   => lcart, &
                            tot_charge_ => tot_charge, &
                            tot_magnetization_ => tot_magnetization, &
                            multiplicity_ => multiplicity
  !
  USE ktetra,        ONLY : nk1, nk2, nk3, k1, k2, k3, ltetra
  !
  USE ldaU,          ONLY : Hubbard_U_     => hubbard_u, &
                            Hubbard_alpha_ => hubbard_alpha, &
                            niter_with_fixed_ns, starting_ns, U_projection, &
                            lda_plus_u_ => lda_plus_u
  !
  USE a2F,           ONLY : la2F_ => la2F
  !
  USE exx,           ONLY : nqx1_   => nq1,  &
                            nqx2_   => nq2,  &
                            nqx3_   => nq3

  USE realus,        ONLY : tqr_    => tqr
  !
  USE lsda_mod,      ONLY : nspin_                  => nspin, &
                            starting_magnetization_ => starting_magnetization, &
                            lsda
  !
  USE relax,         ONLY : epsf, starting_scf_threshold, epse
  !
  USE control_flags, ONLY : diis_ndim, isolve, &
                            max_cg_iter, diis_buff, david, imix, nmix, &
                            iverbosity, tr2, niter, pot_order, wfc_order, &
                            tolp_        => tolp, &
                            upscale_     => upscale, &
                            mixing_beta_ => mixing_beta, &
                            nstep_       => nstep, &
                            iprint_      => iprint, &
                            nosym_       => nosym, &
                            modenum_     => modenum, &
                            reduce_io, ethr, lscf, lbfgs, lmd, lpath, lneb, &
                            lsmd, lphonon, ldamped, lraman, lrescale_t, &
                            lmetadyn, lconstrain, lcoarsegrained, restart, &
                            twfcollect
  !
  USE wvfct,         ONLY : nbnd_ => nbnd
  !
  USE fixed_occ,     ONLY : tfixed_occ
  !
  USE path_variables, ONLY : nstep_path, lsteep_des, lquick_min, lbroyden, &
                             ldamped_dyn, lmol_dyn, llangevin, &
                             ds_                 => ds, &
                             write_save_         => write_save, &
                             use_masses_         => use_masses, &
                             init_num_of_images_ => init_num_of_images, &
                             use_fourier_        => use_fourier, &
                             use_multistep_      => use_multistep, &
                             CI_scheme_          => CI_scheme, &
                             fixed_tan_          => fixed_tan, &
                             free_energy_        => free_energy, &
                             use_freezing_       => use_freezing, &
                             k_max_              => k_max, & 
                             k_min_              => k_min, &
                             num_of_images_      => num_of_images, &
                             first_last_opt_     => first_last_opt, &
                             damp_               => damp, &
                             temp_req_           => temp_req, &
                             path_thr_           => path_thr
  !
  USE noncollin_module, ONLY : i_cons, mcons, &
                               noncolin_  => noncolin, &
                               lambda_    => lambda, &
                               angle1_    => angle1, &
                               angle2_    => angle2, &
                               report_    => report
  !
  USE spin_orb, ONLY : lspinorb_ => lspinorb
  !
  USE bfgs_module,   ONLY : bfgs_ndim_        => bfgs_ndim, &
                            trust_radius_max_ => trust_radius_max, &
                            trust_radius_min_ => trust_radius_min, &
                            trust_radius_ini_ => trust_radius_ini, &
                            w_1_              => w_1, & 
                            w_2_              => w_2
  !
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : title, calculation, verbosity, restart_mode, &
                               nstep, iprint, tstress, tprnfor, dt, outdir, &
                               wfcdir, prefix, etot_conv_thr, forc_conv_thr, &
                               pseudo_dir, disk_io, tefield, dipfield, lberry, &
                               gdir, nppstr, wf_collect,lelfield, efield, tqr, &
                               nberrycyc
  !
  ! ... SYSTEM namelist
  !
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               nat, ntyp, nbnd, nelec, nelup, neldw, &
                               tot_charge, tot_magnetization, multiplicity, &
                               ecutwfc, ecutrho, &
                               nr1, nr2, nr3, nr1s, nr2s, nr3s, &
                               nosym, starting_magnetization, &
                               occupations, degauss, smearing, &
                               nspin, ecfixed, qcutz, q2sigma, &
                               lda_plus_U, Hubbard_U, Hubbard_alpha, &
                               starting_ns_eigenvalue, U_projection_type, &
                               input_dft, la2F, &
#if defined (EXX)                               
                               nqx1, nqx2, nqx3, &
#endif
                               edir, emaxpos, eopreg, eamp, &
                               noncolin, lambda, angle1, angle2, &
                               constrained_magnetization, B_field, &
                               fixed_magnetization, report, lspinorb
  !
  ! ... ELECTRONS namelist
  !
  USE input_parameters, ONLY : electron_maxstep, mixing_mode, mixing_beta, &
                               mixing_ndim, mixing_fixed_ns, conv_thr,     &
                               diago_thr_init, diago_cg_maxiter,           &
                               diago_david_ndim, diago_diis_ndim,          &
                               diagonalization, startingwfc, startingpot
  !
  ! ... IONS namelist
  !
  USE input_parameters, ONLY : phase_space, ion_dynamics, ion_positions, tolp, &
                               tempw, delta_t, nraise, ion_temperature,        &
                               upscale, pot_extrapolation,  wfc_extrapolation, &
                               num_of_images, path_thr, CI_scheme, opt_scheme, &
                               use_masses, use_multistep, first_last_opt,      &
                               damp, init_num_of_images, temp_req, k_max,      &
                               k_min, ds, use_fourier, use_freezing,           &
                               fixed_tan, free_energy, write_save, w_1, w_2,   &
                               trust_radius_max, trust_radius_min, bfgs_ndim,  &
                               trust_radius_ini
  !
  ! ... CELL namelist
  !
  USE input_parameters, ONLY : cell_parameters, cell_dynamics, press, &
                               wmass, cell_temperature, cell_factor
  !
  ! ... PHONON namelist
  !
  USE input_parameters, ONLY : phonon, modenum, xqq
  !
  ! ... RAMAN namelist
  !
  USE input_parameters, ONLY : b_length, lcart
  !
  ! ... "path" specific
  !
  USE input_parameters, ONLY : pos, full_phs_path_flag
  !
  !
  USE constraints_module,    ONLY : init_constraint
  USE metadyn_vars,          ONLY : init_metadyn_vars
  USE read_namelists_module, ONLY : read_namelists, sm_not_set
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER  :: i, ia, ios, is, image, nt
  LOGICAL  :: ltest
  REAL(DP) :: theta, phi
  !
  !
  CALL getenv( 'HOME', pseudo_dir )
  !
  pseudo_dir = TRIM( pseudo_dir ) // '/pw/pseudo/'
  !
  IF ( ionode ) CALL input_from_file()
  !
  ! ... all namelists are read
  !
  CALL read_namelists( 'PW' )
  !
  ! ... translate from input to internals of PWscf, various checks
  !
  if (input_dft /='none') call enforce_input_dft (input_dft)
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
  IF ( tefield .AND. ( nspin > 2 ) ) THEN
     CALL errore( 'iosys', 'LSDA not available with electric field' , 1 )
  END IF
  !
  twfcollect = wf_collect
  !
  ! ... Set Values for electron and bands
  !
  tfixed_occ = .FALSE.
  !
  SELECT CASE( TRIM( occupations ) )
  CASE( 'fixed' )
     !
     ngauss = 0
     ltetra = .FALSE.
     IF ( degauss /= 0.D0 ) THEN
        CALL errore( ' iosys ', &
                   & ' fixed occupations, gauss. broadening ignored', -1 )
        degauss = 0.D0
     END IF
     !
  CASE( 'smearing' )
     !
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
     !
  CASE( 'tetrahedra' )
     !
     ngauss = 0
     ltetra = .TRUE.
     !
  CASE( 'from_input' )
     !
     ngauss     = 0
     ltetra     = .FALSE.
     tfixed_occ = .TRUE.
     !
  CASE DEFAULT
     !
     CALL errore( ' iosys ',' occupations ' // TRIM( occupations ) // &
                & 'not implemented', 1 )
     !
  END SELECT
  !
  IF( nbnd < 1 ) &
     CALL errore( ' iosys ', ' nbnd less than 1 ', nbnd )
  !
  IF( nelec < 0 ) &
     CALL errore( ' iosys ', ' nelec less than 0 ', 1 )
  !
  IF ( nelup < 0 ) &
     CALL errore( ' iosys ', ' nelup less than 0 ', 1 )
  !
  IF ( neldw < 0 ) &
     CALL errore( ' iosys ', ' neldw less than 0 ', 1 )
  !
  SELECT CASE( nspin )
  CASE( 1 ) 
     !
     IF ( noncolin ) nspin = 4
     !
  CASE( 2 )
     !
     lsda = .TRUE.
     IF ( noncolin ) &
        CALL errore( ' iosys ', &
                     ' noncolin .and. nspin==2 are conflicting flags ', 1 )
     !
  CASE( 4 )
     !
     noncolin = .TRUE.
     !
  CASE DEFAULT 
     !
     CALL errore( ' iosys ', 'wrong input value for nspin ', 1 )
     !
  END SELECT
  !
  IF ( nelup == 0.D0 .AND. neldw == 0.D0 .AND. &
       tot_magnetization < 0 .AND. multiplicity == 0) THEN
     !
     two_fermi_energies = .FALSE.
     !
  ELSE
     !
     two_fermi_energies = .TRUE.
     !
     IF ( .NOT. lsda ) &
        CALL errore( ' iosys ', ' fixed nelup/neldw requires nspin=2 ', 1 )
     !
     IF ( ABS( nelup + neldw - nelec ) > 1.D-10 ) &
        CALL errore( ' iosys ', ' nelup + neldw must be equal to nelec ', 1 )
     !
     !
  END IF
  !
  SELECT CASE( TRIM( constrained_magnetization ) )
  CASE( 'none' )
     !
     i_cons = 0
     !
  CASE( 'total' )
     !
     IF ( nspin == 4 ) THEN
        !
        i_cons = 3
        !
        mcons(1,1) = fixed_magnetization(1)
        mcons(2,1) = fixed_magnetization(2)
        mcons(3,1) = fixed_magnetization(3)
        !
     ELSE IF ( nspin == 2 ) THEN
        !
        i_cons = 5
        !
        two_fermi_energies = .TRUE.
        !
        mcons(3,1) = fixed_magnetization(3)
        !
        IF ( fixed_magnetization(1) /= 0.D0 .OR. &
             fixed_magnetization(2) /= 0.D0 ) &
           CALL errore( 'iosys', 'only fixed_magnetization(3)' // &
                      & ' can be specified with nspin=2 ', 1 )
        !
     ELSE
        !
        CALL errore( 'iosys',' constrained total magnetization ' // &
                   & 'requires nspin=2 or 4 ', 1 )
        !
     END IF
     !
  CASE( 'atomic' )
     !
     IF ( nspin == 1 ) &
        CALL errore( 'iosys',' constrained atomic magnetizations ' // &
                   & 'require nspin=2 or 4 ', 1 )
     !
     i_cons = 1
     !
     DO nt = 1, ntyp
        !
        theta = pi * angle1(nt) / 180.D0
        phi   = pi * angle2(nt) / 180.D0
        !
        mcons(1,nt) = starting_magnetization(nt) * SIN( theta ) * COS( phi )
        mcons(2,nt) = starting_magnetization(nt) * SIN( theta ) * SIN( phi )
        mcons(3,nt) = starting_magnetization(nt) * COS( theta )
        !
     END DO
     !
  CASE( 'atomic direction' )
     !
     IF ( nspin == 1 ) &
        CALL errore( 'iosys',' constrained atomic magnetization ' // &
                   & 'directions require nspin=2 or 4 ', 1 )
     !
     i_cons = 2
     !
     DO nt = 1, ntyp
        !
        theta = pi * angle1(nt) / 180.D0
        !
        mcons(3,nt) = cos(theta)
        !
     END DO
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys',' constrained magnetization ' // &
                & TRIM( constrained_magnetization ) // 'not implemented', 1 )
     !
  END SELECT
  !
  IF ( B_field(1) /= 0.D0 .OR. &
       B_field(2) /= 0.D0 .OR. &
       B_field(3) /= 0.D0 ) THEN
     !
     IF ( nspin == 1 ) &
        CALL errore( 'iosys', &
                   & 'non-zero external B_field requires nspin=2 or 4', 1 )
     !
     IF ( TRIM( constrained_magnetization ) /= 'none' ) &
        CALL errore( 'iosys', 'constrained_magnetization and ' // &
                   & 'non-zero external B_field are conflicting flags', 1 )
     !
     IF ( nspin == 2 .AND. &
          ( B_field(1) /= 0.D0 .OR. B_field(2) /= 0.D0 ) ) &
        CALL errore( 'iosys', &
                   & 'only B_field(3) can be specified with nspin=2', 1 )
     !
  END IF
  !
  ! ... starting_magnetization(ia) = sm_not_set means "not set" -- set it to 0
  !
  DO ia = 1, ntyp
     !
     IF ( starting_magnetization(ia) == sm_not_set ) &
        starting_magnetization(ia) = 0.D0
     !
  END DO  
  !
  IF ( ecutrho <= 0.D0 ) THEN
     !
     dual = 4.D0
     !
  ELSE
     !
     dual = ecutrho / ecutwfc
     !
     IF ( dual <= 1.D0 ) &
        CALL errore( ' iosys ', ' invalid dual? ', 1 )
     !
  END IF
  !
  SELECT CASE( TRIM( restart_mode ) )
  CASE( 'from_scratch' )
     !
     restart        = .FALSE.
     startingconfig = 'input'
     !
  CASE( 'restart' )
     !
     IF ( calculation == 'neb' .OR. calculation == 'smd' ) THEN
        !
        ! ... "path" specific
        !
        restart = .FALSE.
        !
     ELSE
        !
        restart = .TRUE.
        !
        IF ( TRIM( ion_positions ) == 'from_input' ) THEN
           !
           startingconfig = 'input'
           !
        ELSE
           !
           startingconfig = 'file'
           !
        END IF
        !
     END IF
     !
  CASE DEFAULT
     !
     CALL errore( ' iosys ', &
                & ' unknown restart_mode ' // TRIM( restart_mode ), 1 )
     !
  END SELECT
  !
  SELECT CASE( TRIM( disk_io ) )
  CASE( 'high' )
     !
     reduce_io = .FALSE.
     !
  CASE DEFAULT
     !
     reduce_io = .TRUE.
     restart   = .FALSE.
     !
  END SELECT
  !
  Hubbard_U(:)    = Hubbard_U(:) / rytoev
  Hubbard_alpha(:)= Hubbard_alpha(:) / rytoev
  !
  ethr = diago_thr_init
  !
  SELECT CASE( TRIM( phase_space ) )
  CASE( 'full' )
     !
     lcoarsegrained  = .FALSE.
     !
  CASE ( 'coarse-grained' )
     !
     lcoarsegrained  = .TRUE.
     !
  END SELECT
  !
  ! ... various initializations of control variables
  !
  lscf      = .FALSE.
  lmd       = .FALSE.
  lmetadyn  = .FALSE.
  lpath     = .FALSE.
  lneb      = .FALSE.
  lsmd      = .FALSE.
  lmovecell = .FALSE.
  lphonon   = .FALSE.
  lraman    = .FALSE.
  lbfgs     = .FALSE.
  ldamped   = .FALSE.
  lconstrain= .FALSE.  
  lforce    = tprnfor
  calc      = ' '
  !
  SELECT CASE( TRIM( calculation ) )
  CASE( 'scf' )
     !
     lscf  = .TRUE.   
     nstep = 1
     !
  CASE( 'nscf' )
     !
     lforce = .FALSE.
     nstep  = 1
     !
  CASE( 'phonon' )
     !
     lforce  = .FALSE.
     lphonon = .TRUE.
     !
     nstep = 1
     !
  CASE( 'relax' )
     !
     lscf   = .TRUE.
     lforce = .TRUE.
     !
     epse = etot_conv_thr
     epsf = forc_conv_thr
     !
     SELECT CASE( TRIM( ion_dynamics ) )
     CASE( 'bfgs' )
        !
        lbfgs = .TRUE.
        !
        IF ( epse <= 20.D0 * ( tr2 / upscale ) ) &
           CALL errore( 'iosys ', 'required etot_conv_thr is too small:' // &
                      & ' conv_thr must be reduced', 1 )   
        !
     CASE( 'constrained-damp' )
        !
        lmd        = .TRUE.
        ldamped    = .TRUE.
        lconstrain = .TRUE.
        !
     CASE ( 'damp' )
        !
        lmd     = .TRUE.
        ldamped = .TRUE.
        !
        ntcheck = nstep + 1
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
  CASE( 'md' )
     !
     lscf   = .TRUE.
     lmd    = .TRUE.          
     lforce = .TRUE.
     !
     SELECT CASE( TRIM( ion_dynamics ) )
     CASE( 'verlet' )
        !
        CONTINUE
        !
     CASE( 'constrained-verlet' )
        !
        lconstrain = .TRUE.
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     END SELECT
     !
  CASE( 'vc-relax' )
     !
     lscf      = .TRUE.
     lmd       = .TRUE.
     lmovecell = .TRUE.
     lforce    = .TRUE.
     ldamped   = .TRUE.
     !
     epse = etot_conv_thr
     epsf = forc_conv_thr
     !
     SELECT CASE( TRIM( cell_dynamics ) )
     CASE( 'none' )
        !
        calc    = 'mm'
        ntcheck = nstep + 1
        !
     CASE( 'damp-pr' )
        !
        calc    = 'cm'
        ntcheck = nstep + 1
        !
     CASE( 'damp-w' )
        !
        calc    = 'nm'
        ntcheck = nstep + 1
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': cell_dynamics=' // TRIM( cell_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
     IF ( TRIM( ion_dynamics ) /= 'damp' ) &
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     !
  CASE( 'vc-md' )
     !
     lscf      = .TRUE.
     lmd       = .TRUE.
     lmovecell = .TRUE.
     lforce    = .TRUE.
     !
     ntcheck = nstep + 1
     !
     SELECT CASE( TRIM( cell_dynamics ) )
     CASE( 'none' )
        !
        calc = 'md'
        !
     CASE( 'pr' )
        !
        calc = 'cd'
        !
     CASE( 'w' )
        !
        calc = 'nd'
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
     IF ( TRIM( ion_dynamics ) /= 'beeman' ) &
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
     !
  CASE( 'raman' )
     !
     lraman    = .TRUE.
     nstep     = 1
     !
  CASE( 'neb' )
     !
     lscf  = .TRUE.
     lpath = .TRUE.
     lneb  = .TRUE.
     !
  CASE( 'smd' )
     !
     lscf  = .TRUE.
     lpath = .TRUE.
     lsmd  = .TRUE.
     !
  CASE( 'metadyn' )
     !
     lscf           = .TRUE.
     lmetadyn       = .TRUE.
     lcoarsegrained = .TRUE.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys ', ' calculation ' // &
                & TRIM( calculation ) // ' not implemented', 1 )
     !
  END SELECT
  !
  IF ( startingpot /= 'atomic' .AND. startingpot /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingpot: use default', -1 )
     !
     IF (       lscf ) startingpot = 'atomic'
     IF ( .NOT. lscf ) startingpot = 'file'
     !
  END IF
  !
  IF ( .NOT. lscf .AND. startingpot /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingpot: use default', -1 )
     !
     startingpot = 'file'
     !
  END IF
  !
  IF ( startingwfc /= 'atomic' .AND. &
       startingwfc /= 'random' .AND. &
       startingwfc /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingwfc: use default', -1 )
     !
     startingwfc = 'atomic'
     !
  END IF
  !
  SELECT CASE( TRIM( diagonalization ) )
  CASE ( 'cg' )
     !
     isolve = 1
     !
     max_cg_iter = diago_cg_maxiter
     !
  CASE ( 'diis' )
     !
     isolve = 2
     !
     max_cg_iter = diago_cg_maxiter
     diis_ndim   = diago_diis_ndim
     !
  CASE ( 'david' )
     !
     isolve = 0
     !
     david = diago_david_ndim
     !
  CASE DEFAULT
     !
     isolve = 0
     !
     david = diago_david_ndim
     !
  END SELECT
  !
  tr2   = conv_thr
  niter = electron_maxstep
  !
  SELECT CASE( TRIM( pot_extrapolation ) )
  CASE( 'none' )
     !
     pot_order = 0
     !
  CASE( 'atomic' )
     !
     pot_order = 1
     !
  CASE( 'first_order' )
     !
     pot_order = 2
     !
  CASE( 'second_order' )
     !
     pot_order = 3
     !
  CASE DEFAULT
     !
     pot_order = 1
     !
  END SELECT
  !
  SELECT CASE( TRIM( wfc_extrapolation ) )
  CASE( 'none' )
     !
     wfc_order = 0
     !
  CASE( 'first_order' )
     !
     wfc_order = 2
     !
  CASE( 'second_order' )
     !
     wfc_order = 3
     !
  CASE DEFAULT
     !
     wfc_order = 0
     !
  END SELECT
  !
  IF ( wfc_order > 0 .AND. noncolin ) THEN
     !
     CALL errore( 'iosys ', &
                & 'wfc extrapolation not implemented in the ' // &
                & 'noncollinear case', -1 )
     !
  END IF
  !
  IF ( occupations == 'fixed' .AND. nspin == 2  .AND. lscf ) THEN
     !
     IF ( two_fermi_energies ) THEN
        !
        IF ( ABS( NINT( nelup ) - nelup ) > 1.D-10 ) &
           CALL errore( 'iosys ', &
                      & 'fixed occupations requires integer nelup ', 1 )
        IF ( ABS( NINT( neldw ) - neldw ) > 1.D-10 ) &
           CALL errore( 'iosys ', &
                      & 'fixed occupations requires integer neldw ', 1 )
        !
     ELSE
        !
        CALL errore( 'iosys ', &
                   & 'fixed occupations and lsda need nelup and neldw ', 1 )
        !
     END IF
     !
  END IF
  !
  IF ( lcoarsegrained ) THEN
     !
     lmd        = .TRUE.
     lconstrain = .TRUE.
     !
     SELECT CASE( TRIM( ion_dynamics ) )
     CASE( 'constrained-verlet' )
        !
        CONTINUE
        !
     CASE( 'constrained-damp' )
        !
        ldamped = .TRUE.
        !
        epse = etot_conv_thr
        epsf = forc_conv_thr
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
  END IF
  !
  ! ... "path" specific initialization of control variables
  !
  IF ( lpath ) THEN
     !
     nstep_path = nstep
     !
     IF ( num_of_images < 2 ) &
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': num_of_images must be at least 2', 1 )
     !
     IF ( ( CI_scheme /= "no-CI"      ) .AND. &
          ( CI_scheme /= "highest-TS" ) .AND. &
          ( CI_scheme /= "all-SP"     ) .AND. &
          ( CI_scheme /= "manual"     ) ) THEN
        !
        CALL errore( 'iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': unknown CI_scheme', 1 )  
        !   
     END IF
     !
     ! ... initialization of logical variables
     !
     lsteep_des  = .FALSE.     
     lquick_min  = .FALSE.
     lbroyden    = .FALSE.
     ldamped_dyn = .FALSE.
     lmol_dyn    = .FALSE.     
     !
     SELECT CASE( opt_scheme )
     CASE( "sd" )
        !
        lsteep_des = .TRUE.
        !
     CASE( "quick-min" )
        !
        lquick_min = .TRUE.
        !
     CASE( "broyden" )
        !
        lbroyden     = .TRUE.
        !
     CASE( "damped-dyn" )
        !
        ldamped_dyn = .TRUE.
        !
     CASE( "mol-dyn" )
        !
        lmol_dyn = .TRUE.
        !
     CASE( "langevin" )
        !
        llangevin = .TRUE.
        !
        IF ( lneb ) &
           CALL errore( 'iosys ','calculation=' // TRIM( calculation ) // &
                      & ': langevin dynamics not implemented', 1 )
        !
        temp_req = temp_req / ( eV_to_kelvin * au )
        !
        IF ( temp_req <= 0.D0 ) &
           CALL errore( 'iosys ','calculation=' // TRIM( calculation ) // &
                      & ': tepm_req has not been set', 1 )
        !
        IF ( use_freezing ) &
           WRITE( UNIT = stdout, &
                  FMT = '(5X,"warning: freezing cannot be used in langevin")' )
        !
        IF ( use_multistep ) &
           WRITE( UNIT = stdout, &
                  FMT = '(5X,"warning: multistep cannot be used in langevin")' )
        !
        use_freezing  = .FALSE.
        use_multistep = .FALSE.
        !
     CASE default
        !
        CALL errore( 'iosys ','calculation=' // TRIM( calculation ) // &
                   & ': unknown opt_scheme', 1 )  
        !
     END SELECT             
     !
  END IF
  !
  SELECT CASE( TRIM( ion_temperature ) )
  CASE( 'not_controlled' )
     !
     lrescale_t = .FALSE.
     !
  CASE( 'rescaling' )
     !
     lrescale_t  = .TRUE.
     temperature = tempw
     tolp_       = tolp
     !
  CASE DEFAULT
     !
     CALL errore( ' iosys ', 'unknown ion_temperature ' // &
                & TRIM( ion_temperature ), 1 )
     !
  END SELECT
  !
  SELECT CASE( TRIM( mixing_mode ) )
  CASE( 'plain' )
     !
     imix = 0
     !
  CASE( 'TF' )
     !
     imix = 1
     !
  CASE( 'local-TF' )
     !
     imix = 2
     !
  CASE( 'potential' )
     !
     CALL errore( 'iosys ', 'potential mixing no longer implemented', 1 )
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys ', 'unknown mixing ' // TRIM( mixing_mode ), 1 )
     !
  END SELECT
  !
  IF ( imix > 0 .AND. noncolin )  CALL errore( ' iosys ', &
      ' only simple mixing is implemented in noncolinear calculations', 1 )
  !
  starting_scf_threshold = tr2
  nmix = mixing_ndim
  niter_with_fixed_ns = mixing_fixed_ns
  !
  SELECT CASE( TRIM( verbosity ) )
  CASE( 'high' )
     !
     iverbosity = 1
     !
  CASE DEFAULT
     !
     iverbosity = 0
     !
  END SELECT
  !
  tmp_dir = TRIM( outdir )
  lstres = ( tstress .AND. lscf )
  !
  IF ( lberry .AND. npool > 1 ) &
     CALL errore( 'iosys ', 'Berry Phase not implemented with pools ', 1 )
  !
  IF ( lberry .AND. nproc_pool > 1 .AND. gdir /= 3 ) &
     CALL errore( 'iosys ', 'Berry Phase in parallel only for gdir=3', 1 )
  !
  ! ... Copy values from input module to PW internals
  !
  nppstr_     = nppstr
  gdir_       = gdir
  lberry_     = lberry
  lelfield_   = lelfield
  efield_     = efield
  nberrycyc_  = nberrycyc
  tqr_        = tqr

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
  nelup_   = nelup
  neldw_   = neldw
  tot_charge_ = tot_charge
  tot_magnetization_ = tot_magnetization
  multiplicity_ = multiplicity
  !
  lspinorb_ = lspinorb
  noncolin_ = noncolin
  angle1_   = angle1
  angle2_   = angle2
  report_   = report
  lambda_   = lambda
  !
  Hubbard_U_(1:ntyp)      = hubbard_u(1:ntyp)
  Hubbard_alpha_(1:ntyp)  = hubbard_alpha(1:ntyp)
  lda_plus_u_             = lda_plus_u
  la2F_                   = la2F
  nspin_                  = nspin
  starting_magnetization_ = starting_magnetization
  starting_ns             = starting_ns_eigenvalue
  U_projection            = U_projection_type
  nosym_                  = nosym
  nbnd_                   = nbnd
  !
#if defined (EXX)
  !
  nqx1_   = nqx1
  nqx2_   = nqx2
  nqx3_   = nqx3
  !
#endif
  !
  startingwfc_ = startingwfc
  startingpot_ = startingpot
  mixing_beta_ = mixing_beta
  !
  upscale_     = upscale
  delta_t_     = delta_t
  nraise_      = nraise
  press_       = press
  cell_factor_ = cell_factor
  modenum_     = modenum
  xqq_         = xqq
  !
  b_length_    = b_length
  lcart_       = lcart
  !
  ! ... general "path" variables
  !
  ds_             = ds
  num_of_images_  = num_of_images
  first_last_opt_ = first_last_opt
  use_masses_     = use_masses
  write_save_     = write_save
  free_energy_    = free_energy
  use_fourier_    = use_fourier
  use_freezing_   = use_freezing
  damp_           = damp
  temp_req_       = temp_req
  path_thr_       = path_thr 
  !
  ! ... NEB specific
  !
  CI_scheme_ = CI_scheme
  k_max_     = k_max 
  k_min_     = k_min
  !
  ! ... SMD specific
  !
  init_num_of_images_ = init_num_of_images
  use_multistep_      = use_multistep
  fixed_tan_          = fixed_tan 
  !
  ! ... BFGS specific
  !
  bfgs_ndim_        = bfgs_ndim
  trust_radius_max_ = trust_radius_max
  trust_radius_min_ = trust_radius_min
  trust_radius_ini_ = trust_radius_ini
  w_1_              = w_1
  w_2_              = w_2
  !
  ! ... read following cards
  !
  ALLOCATE( ityp( nat_ ) )
  ALLOCATE( tau(    3, nat_ ) )
  ALLOCATE( force(  3, nat_ ) )  ! ... compatibility with old readin
  ALLOCATE( if_pos( 3, nat_ ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat_ ) )
  !
  CALL read_cards( psfile, atomic_positions )
  !
  ! ... set up atomic positions and crystal lattice
  !
  IF ( celldm_(1) == 0.D0 .AND. a /= 0.D0 ) THEN
     !
     IF ( ibrav_ == 0 ) ibrav = 14
     !
     celldm_(1) = a / bohr_radius_angs
     celldm_(2) = b / a
     celldm_(3) = c / a
     celldm_(4) = cosab
     celldm_(5) = cosac
     celldm_(6) = cosbc
     !
  ELSE IF ( celldm_(1) /= 0.D0 .AND. a /= 0.D0 ) THEN
     !
     CALL errore( 'input', ' do not specify both celldm and a,b,c!', 1 )
     !
  END IF
  !
  ! ... generate at (in atomic units) from ibrav and celldm
  !
  CALL latgen( ibrav_, celldm_, at(1,1), at(1,2), at(1,3), omega )
  !
  ! ... define alat
  !
  alat = celldm_(1)
  !
  ! ... convert at to unit of alat
  !
  at = at / alat
  !
  ! ... Generate the reciprocal lattice vectors
  !
  CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  !
  IF ( full_phs_path_flag ) THEN
     !
     ! ... "path" optimizations specific
     !
     DO image = 1, num_of_images_
        !
        tau = RESHAPE( pos(1:3*nat_,image), (/ 3 , nat_ /) )
        !
        ! ... convert input atomic positions to internally used format:
        ! ... tau in a0 units
        !
        SELECT CASE( atomic_positions )
        CASE( 'alat' )
           !
           ! ... input atomic positions are divided by a0: do nothing
           !
        CASE( 'bohr' )
           !
           ! ... input atomic positions are in a.u.: divide by alat
           !
           tau = tau / alat
           !
        CASE( 'crystal' )
           !
           ! ... input atomic positions are in crystal axis
           !
           CALL cryst_to_cart( nat_, tau, at, 1 )
           !
        CASE( 'angstrom' )
           !
           ! ... atomic positions in A: convert to a.u. and divide by alat
           !
           tau = tau / bohr_radius_angs / alat
           !
        CASE DEFAULT
           !
           CALL errore( 'iosys ','atomic_positions=' // &
                      & TRIM( atomic_positions ) // ' not implemented ', 1 )
           !
        END SELECT
        !
        pos(1:3*nat_,image) = RESHAPE( tau, (/ 3 * nat_ /) )
        !
     END DO 
     !
  ELSE
     !
     ! ... convert input atomic positions to internally used format:
     ! ... tau in a0 units
     !
     SELECT CASE( atomic_positions )
     CASE( 'alat' )
        !
        ! ... input atomic positions are divided by a0: do nothing
        !
     CASE( 'bohr' )
        !
        ! ... input atomic positions are in a.u.: divide by alat
        !
        tau = tau / alat
        !
     CASE( 'crystal' )
        !
        ! ... input atomic positions are in crystal axis
        !
        CALL cryst_to_cart( nat_, tau, at, 1 )
        !
     CASE( 'angstrom' )
        !
        ! ... atomic positions in A: convert to a.u. and divide by alat
        !
        tau = tau / bohr_radius_angs / alat
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys ','atomic_positions=' // &
                   & TRIM( atomic_positions ) // ' not implemented ', 1 )
        !
     END SELECT
     !
  END IF
  !
  ! ... Renata's dynamics uses masses in atomic units
  !
  IF ( calc /= ' ' ) amass = amass * amconv
  !
  ! ... set default value of wmass
  !
  IF ( wmass == 0.D0 ) THEN
     !
     IF ( calc == 'nd' .OR. calc == 'nm' ) THEN
        !
        wmass = SUM( amass(ityp(:)) )
        !
        wmass = 0.75D0 * wmass / pi / pi / omega**( 2.D0 / 3.D0 )
        !
     END IF
     !
     IF ( calc == 'cd' .OR. calc == 'cm' ) THEN
        !
        wmass = SUM( amass(ityp(:)) )
        !
        wmass = 0.75D0 * wmass / pi / pi
        !
     END IF
     !
     cmass  = wmass
     !
  ELSE
     !
     cmass  = wmass * amconv
     !
  END IF
  !
  ! ... unit conversion for pressure
  !
  press_ = press_ / uakbar
  !
  ! ... read pseudopotentials
  !
  CALL readpp()
  !
  ! ... In the case of variable cell dynamics save old cell variables
  ! ... and initialize a few other variables
  !
  IF ( lmovecell ) THEN
     !
     at_old    = at
     omega_old = omega
     lstres    = .TRUE.
     !
     IF ( cell_factor_ <= 0.D0 ) cell_factor_ = 1.2D0
     !
     IF ( cmass <= 0.D0 ) &
        CALL errore( 'readin', &
                   & 'vcsmd: a positive value for cell mass is required', 1 )
     !
  ELSE
     !
     cell_factor_ = 1.D0
     !
  END IF
  !
  ! ... set constraints
  !
  IF ( lconstrain ) CALL init_constraint( nat, tau, alat, ityp )
  !
  ! ... set variables for metadynamics
  !
  IF ( lcoarsegrained ) CALL init_metadyn_vars()
  !
  CALL verify_tmpdir( tmp_dir )
  !
  IF ( .NOT. TRIM( wfcdir ) == 'undefined' ) THEN
     !
     wfc_dir = TRIM( wfcdir )
     !
     CALL verify_tmpdir( wfc_dir )
     !
  ENDIF
  !
  CALL restart_from_file()
  !
  IF ( startingconfig == 'file' ) CALL read_config_from_file()
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
!----------------------------------------------------------------------------
SUBROUTINE read_cards( psfile, atomic_positions_ )
  !----------------------------------------------------------------------------
  !
  USE input_parameters,   ONLY : atom_label, atom_pfile, atom_mass, taspc, &
                                 tapos, rd_pos, atomic_positions, if_pos,  &
                                 sp_pos, k_points, xk, wk, nk1, nk2, nk3,  &
                                 k1, k2, k3, nkstot, cell_symmetry, rd_ht, &
                                 trd_ht, f_inp, calculation
  USE wvfct,              ONLY : gamma_only
  USE cell_base,          ONLY : at, ibrav, symm_type, celldm
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp, tau, atm
  USE klist,              ONLY : nks
  USE ktetra,             ONLY : nk1_   => nk1, &
                                 nk2_   => nk2, &
                                 nk3_   => nk3, &
                                 k1_    => k1,  &
                                 k2_    => k2,  &
                                 k3_    => k3
  USE klist,              ONLY : lxkcry, &
                                 xk_    => xk, &
                                 wk_    => wk
  USE fixed_occ,          ONLY : tfixed_occ, &
                                 f_inp_ => f_inp
  USE ions_base,          ONLY : fixatom, &
                                 if_pos_ =>  if_pos
  USE dynam,              ONLY : amass
  USE control_flags,      ONLY : lfixatom
  USE read_cards_module,  ONLY : read_cards_base => read_cards
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256) :: psfile(ntyp)
  CHARACTER (LEN=30)  :: atomic_positions_
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
       CALL errore( 'cards ', 'atomic species info missing', 1 )
  IF ( .NOT. tapos ) &
       CALL errore( 'cards ', 'atomic position info missing', 1 )
  !
  DO is = 1, ntyp
     !
     amass(is)  = atom_mass(is)
     psfile(is) = atom_pfile(is)
     atm(is)    = atom_label(is)
     !
     IF ( amass(is) <= 0.D0 ) &
        CALL errore( ' iosys ', ' invalid  mass ', is )
     !
  END DO
  !
  DO ia = 1, nat
     !
     tau(:,ia) = rd_pos(:,ia)
     ityp(ia)  = sp_pos(ia)
     !
  END DO
  !
  ! ... calculate fixatom
  !
  fixatom = 0
  !
  IF ( ANY( if_pos(:,1:nat) == 0 ) ) lfixatom = .TRUE.
  !
  DO ia = 1, nat
     !
     IF ( if_pos(1,ia) /= 0 .OR. &
          if_pos(2,ia) /= 0 .OR. &
          if_pos(3,ia) /= 0 ) CYCLE
     !
     fixatom = fixatom + 1
     !
  END DO
  !
  ! ... The constrain on fixed coordinates is implemented using the array 
  ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1 
  ! ... otherwise. fixatom is maintained for compatibility. ( C.S. 15/10/2003 )
  !
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
    !
    ! ... nk1,nk2,nk3 and k1,k2,k3 are initialized even when not used
    !
    nk1_ = nk1
    nk2_ = nk2
    nk3_ = nk3
    k1_  = k1
    k2_  = k2
    k3_  = k3
    !
  ELSE IF ( k_points == 'tpiba' ) THEN
    !
    ! ... input k-points are in 2pi/a units
    !
    gamma_only   = .FALSE.
    lxkcry       = .FALSE.
    nks          = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    !
    nk1_ = 0
    nk2_ = 0
    nk3_ = 0
    k1_  = 0
    k2_  = 0
    k3_  = 0
    !
  ELSE IF ( k_points == 'crystal' ) THEN
    !
    ! ... input k-points are in crystal (reciprocal lattice) axis
    !
    gamma_only   = .FALSE.
    lxkcry       = .TRUE.
    nks          = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    !
    nk1_ = 0
    nk2_ = 0
    nk3_ = 0
    k1_  = 0
    k2_  = 0
    k3_  = 0
    !
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
    !
    nk1_ = 0
    nk2_ = 0
    nk3_ = 0
    k1_  = 0
    k2_  = 0
    k3_  = 0
    !
  ELSE
    !
    ! ... default: input k-points are in 2pi/a units
    !
    gamma_only   = .FALSE.
    lxkcry       = .FALSE.
    nks          = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    !
    nk1_ = 0
    nk2_ = 0
    nk3_ = 0
    k1_  = 0
    k2_  = 0
    k3_  = 0
    !
  END IF
  !
  IF ( tfixed_occ ) THEN
     !
     IF ( nks > 1 .OR. ( nk1 * nk2 * nk3 ) > 1 ) &
        CALL errore( 'read_cards', &
                   & 'only one k point with fixed occupations', 1 )
     !
     f_inp_ = f_inp
     !
  END IF
  !
  IF ( trd_ht ) THEN
    !
    symm_type = cell_symmetry 
    at        = TRANSPOSE( rd_ht )
    tcell     = .TRUE.
    !
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
!-----------------------------------------------------------------------
SUBROUTINE verify_tmpdir( tmp_dir )
  !-----------------------------------------------------------------------
  !
  USE input_parameters, ONLY : restart_mode
  USE control_flags,    ONLY : lpath
  USE io_files,         ONLY : prefix, nd_nmbr
  USE path_variables,   ONLY : num_of_images
  USE mp_global,        ONLY : mpime, nproc
  USE io_global,        ONLY : ionode
  USE mp,               ONLY : mp_barrier
  USE parser,           ONLY : int_to_char, delete_if_present
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(INOUT) :: tmp_dir
  !
  INTEGER             :: l, ios, image, proc
  CHARACTER (LEN=256) :: file_path, tmp_dir_saved
  INTEGER, EXTERNAL   :: c_mkdir
  !
  !
  ! ... verify if tmp_dir ends with /, add one if needed
  !
  l = LEN_TRIM( tmp_dir )
  !
  IF ( tmp_dir(l:l) /= '/' ) THEN
     !
     IF ( l > 0 .AND. l < LEN( tmp_dir ) ) THEN
        !
        tmp_dir(l+1:l+1) = '/'
        !
     ELSE
        !
        CALL errore( 'outdir: ', tmp_dir // ' truncated or empty', 1 )
        !
     END IF
     !
  END IF
  !
  ios = 0
  !
  file_path = TRIM( tmp_dir ) // 'pwscf'
  !
  OPEN( UNIT = 4, FILE = TRIM( file_path ) // TRIM( int_to_char( mpime ) ), &
      & STATUS = 'UNKNOWN',  FORM = 'UNFORMATTED', IOSTAT = ios )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  IF ( ios /= 0 ) CALL errore( 'outdir: ', TRIM( tmp_dir ) // &
                             & ' non existent or non writable', 1 )
  !
  ! ... if starting from scratch all temporary files are removed
  ! ... from tmp_dir ( only by the master node )  
  !
  file_path = TRIM( tmp_dir ) // TRIM( prefix )  
  !
  IF ( restart_mode == 'from_scratch' ) THEN
     !
     IF ( ionode ) THEN
        !
        ! ... extrapolation file is removed
        !     
        CALL delete_if_present( TRIM( file_path ) // '.update' )
        !
        ! ... MD restart file is removed
        !
        CALL delete_if_present( TRIM( file_path ) // '.md' )
        !
        ! ... BFGS restart file is removed       
        !     
        CALL delete_if_present( TRIM( file_path ) // '.bfgs' )
        !
     END IF
     !
  END IF    
  !
  ! ... "path" optimisation specific :
  ! ... in the scratch directory the tree of subdirectories needed by "path"
  ! ... calculations are created
  !
  IF ( lpath ) THEN
     !
     IF ( ionode ) THEN
        !
        ! ... files needed by parallelization among images are removed
        !               
        CALL delete_if_present( TRIM( file_path ) // '.BLOCK' )
        CALL delete_if_present( TRIM( file_path ) // '.newimage' )
        !
        ! ... file containing the broyden's history
        !
        IF ( restart_mode == 'from_scratch' ) THEN
           !
           CALL delete_if_present( TRIM( tmp_dir ) // &
                                 & TRIM( prefix ) // '.broyden' )
           !
        END IF
        !
     END IF   
     !
     tmp_dir_saved = tmp_dir
     ! 
     DO image = 1, num_of_images
        !
        tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) //"_" // &
                  TRIM( int_to_char( image ) ) // '/'
        !
        DO proc = 0, nproc - 1
           !
           ! ... a scratch directory for this image is created sequentially
           ! ... by all the cpus
           !
           IF ( proc == mpime ) &
              ios = c_mkdir( TRIM( tmp_dir ), LEN_TRIM( tmp_dir ) )
           !
           CALL mp_barrier()
           !
        END DO
        !
        ! ... each job checks whether the scratch directory is accessible
        ! ... or not
        !
        OPEN( UNIT = 4, FILE = TRIM( tmp_dir ) // TRIM( prefix ) // &
            & TRIM( int_to_char( mpime ) ), &
            & STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
        CLOSE( UNIT = 4, STATUS = 'DELETE' ) 
        ! 
        IF ( ios /= 0 ) &
           CALL errore( 'outdir: ', TRIM( tmp_dir ) // &
                      & ' non existent or non writable', 1 )
        !
        ! ... if starting from scratch all temporary files are removed 
        ! ... from tmp_dir ( by all the cpus in sequence )
        !
        IF ( restart_mode == 'from_scratch' ) THEN
           !
           DO proc = 0, nproc - 1
              !
              IF ( proc == mpime ) THEN
                 !
                 ! ... BFGS restart file is removed       
                 !
                 CALL delete_if_present( TRIM( tmp_dir ) // &
                                       & TRIM( prefix ) // '.bfgs' )
                 !
                 ! ... extrapolation file is removed
                 !     
                 CALL delete_if_present( TRIM( tmp_dir ) // &
                                       & TRIM( prefix ) // '.update' )
                 !
                 ! ... standard output of the self-consistency is removed
                 !      
                 CALL delete_if_present( TRIM( tmp_dir ) // 'PW.out' )
                 !
              END IF
              !
              CALL mp_barrier()
              !
           END DO
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
