!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE input
   !---------------------------------------------------------------------------
   !
   USE kinds,      ONLY: DP
   USE read_input, ONLY : has_been_read
   !
   IMPLICIT NONE
   SAVE
   !
   PRIVATE                      !  Input Subroutines
                                !  should be called in the following order
                                !  a) read input file (module read_input)
   PUBLIC :: iosys_pseudo       !  b) then read pseudo files
   PUBLIC :: iosys              !  c) finally copy variables to modules
   PUBLIC :: modules_setup, set_control_flags
   !
   CHARACTER(LEN=256), EXTERNAL :: trimcheck
   !
   CONTAINS
   !
   !-------------------------------------------------------------------------
   SUBROUTINE iosys_pseudo()
     !-------------------------------------------------------------------------
     !
     USE input_parameters,        ONLY : atom_pfile, pseudo_dir, ntyp, nat, &
                                         prefix, outdir, input_dft
     USE read_pseudo_mod,         ONLY : readpp, check_order
     USE io_global,               ONLY : stdout
     USE io_files,                ONLY : psfile_     => psfile , &
                                         pseudo_dir_ => pseudo_dir, &
                                         prefix_     => prefix, &
                                         tmp_dir
     USE ions_base,               ONLY : nsp_ => nsp, nat_ => nat
     USE input_parameters,        ONLY : title
     USE run_info,                ONLY : title_ => title
     !
     !
     IF ( .NOT. has_been_read ) &
        CALL errore( 'iosys_pseudo ', 'input file has not been read yet!', 1 )
     !
     ! ... Set job title and print it on standard output
     !
     title_ = title
     WRITE( stdout, '(/,3X,"Job Title: ",A )' ) TRIM( title_ )
     !
     prefix_  = TRIM( prefix  )
     tmp_dir  = trimcheck( outdir )
     !
     ! ... Set internal variables for the number of species and number of atoms
     !
     nsp_ = ntyp
     nat_ = nat
     !
     psfile_         = ' '
     psfile_(1:nsp_) = atom_pfile(1:nsp_)
     pseudo_dir_     = trimcheck( pseudo_dir  )
     !
     ! ... read in pseudopotentials and wavefunctions files
     !
     CALL readpp( input_dft, .TRUE. )
     CALL check_order ( )
     !
     RETURN
     !
   END SUBROUTINE iosys_pseudo
   !
   !-------------------------------------------------------------------------
   SUBROUTINE iosys()
     !-------------------------------------------------------------------------
     !
     USE control_flags,      ONLY : fix_dependencies, lconstrain
     USE io_global,          ONLY : ionode, stdout
     USE ions_base,          ONLY : nat, tau, ityp
     USE constraints_module, ONLY : init_constraint
     !
     IMPLICIT NONE
     !
     !
     IF ( ionode ) THEN
        !
        WRITE( UNIT = stdout, &
               FMT = "(//,3X,'Main Simulation Parameters (from input)',/ &
                     &   ,3X,'---------------------------------------')" )
        !
     END IF
     !
     ! ... Set internal flags according to the input
     !
     CALL set_control_flags()
     !
     ! ... write to stdout basic simulation parameters
     !
     CALL input_info()
     !
     ! ... call the module specific setup routine
     !
     CALL modules_setup()
     !
     IF ( lconstrain ) CALL init_constraint( nat, tau, ityp, 1.D0 )
     !
     ! ... write to stdout input module information
     !
     CALL modules_info()
     !
     RETURN
     !
   END SUBROUTINE iosys
   !
   !-------------------------------------------------------------------------
   SUBROUTINE set_control_flags()
     !-------------------------------------------------------------------------
     !
     USE io_global,     ONLY : stdout
     USE autopilot,     ONLY : auto_check
     USE autopilot,     ONLY : restart_p
     USE control_flags, ONLY : ndw_        => ndw, &
                               ndr_        => ndr, &
                               iprint_     => iprint, &
                               isave_      => isave, &
                               tstress_    => tstress, &
                               tprnfor_    => tprnfor, &
                               ampre_      => ampre, &
                               trane_      => trane, &
                               nomore_     => nomore, &
                               memchk_     => memchk, &
                               tpre_       => tpre, &
                               timing_     => timing, &
                               iverbosity_ => iverbosity, &
                               taurdr_     => taurdr, &
                               nbeg_       => nbeg, &
                               gamma_only_ => gamma_only, &
                               tortho_     => tortho,   &
                               nstep_      => nstep
     USE control_flags, ONLY : tsde_          => tsde, &
                               tzeroe_        => tzeroe, &
                               trhor_         => trhor, &
                               trhow_         => trhow, &
                               tksw_          => tksw,  &
                               ortho_eps_     => ortho_eps, &
                               ortho_max_     => ortho_max, &
                               tnosee_        => tnosee
     USE control_flags, ONLY : tfor_      => tfor, &
                               tsdp_      => tsdp
     USE control_flags, ONLY : tnosep_ => tnosep, &
                               tcap_   => tcap, &
                               tcp_    => tcp, &
                               tolp_   => tolp, &
                               tzerop_ => tzerop, &
                               tv0rd_  => tv0rd, &
                               tranp_  => tranp, &
                               amprp_  => amprp, &
                               dt_old_ => dt_old
     USE control_flags, ONLY : tionstep_ => tionstep, &
                               nstepe_   => nstepe
     USE control_flags, ONLY : tzeroc_ => tzeroc, &
                               tnoseh_ => tnoseh, &
                               thdyn_  => thdyn, &
                               tsdc_   => tsdc, &
                               tbeg_   => tbeg
     USE control_flags, ONLY : ekin_conv_thr_ => ekin_conv_thr, &
                               etot_conv_thr_ => etot_conv_thr, &
                               forc_conv_thr_ => forc_conv_thr, &
                               ekin_maxiter_  => ekin_maxiter, &
                               etot_maxiter_  => etot_maxiter, &
                               forc_maxiter_  => forc_maxiter
     USE control_flags, ONLY : force_pairing_ => force_pairing
     USE control_flags, ONLY : remove_rigid_rot_ => remove_rigid_rot
     USE control_flags, ONLY : iesr_ => iesr
     USE control_flags, ONLY : textfor
     USE control_flags, ONLY : do_makov_payne, twfcollect
     USE control_flags, ONLY : lwf, lwfnscf, lwfpbe0nscf
     USE control_flags, ONLY : smallmem
     USE control_flags, ONLY : tconvthrs
     !
     ! ...  Other modules
     !
     USE cp_main_variables,        ONLY : nprint_nfi
     USE wave_base,          ONLY : frice_ => frice
     USE ions_base,          ONLY : fricp_ => fricp
     USE cell_base,          ONLY : frich_ => frich
     USE time_step,          ONLY : set_time_step
     USE cp_electronic_mass, ONLY : emass_ => emass, &
                                    emaec_ => emass_cutoff
     !
     USE efield_module,      ONLY : tefield_    => tefield,  &
                                    epol_       => epol,     &
                                    efield_     => efield,   &
                                    tefield2_    => tefield2,  &
                                    epol2_       => epol2,     &
                                    efield2_     => efield2
     !
     USE uspp_param,         ONLY : nvb
     !
     USE input_parameters,   ONLY: &
        electron_dynamics, electron_damping, electron_temperature,   &
        ion_dynamics, ekin_conv_thr, etot_conv_thr, forc_conv_thr, ion_maxstep,&
        electron_maxstep, ion_damping, ion_temperature, ion_velocities, tranp, &
        amprp, ion_nstepe, cell_nstepe, cell_dynamics, cell_damping, london,   &
        cell_parameters, cell_velocities, cell_temperature, force_pairing,     &
        tapos, tavel, ecutwfc, emass, emass_cutoff, taspc, trd_ht, ibrav,      &
        ortho_eps, ortho_max, ntyp, tolp, calculation, disk_io, dt,            &
        tcg, ndr, ndw, iprint, isave, tstress, k_points, tprnfor, verbosity,   &
        ampre, nstep, restart_mode, ion_positions, startingwfc,                &
        orthogonalization, electron_velocities, nat, if_pos,                   &
        tefield, epol, efield, tefield2, epol2, efield2, remove_rigid_rot,     &
        iesr, saverho, rd_for, assume_isolated, wf_collect,                    &
        memory, ref_cell, tcpbo
     USE funct,              ONLY : dft_is_hybrid
     !
     IMPLICIT NONE
     !
     !
     IF ( .NOT. has_been_read ) &
        CALL errore( 'iosys ', 'input file has not been read yet!', 1 )
     !
     ndr_           = ndr
     ndw_           = ndw
     iprint_        = iprint
     isave_         = isave
     tstress_       = tstress
     tpre_          = tstress
     gamma_only_    = ( TRIM( k_points ) == 'gamma' )
     tprnfor_       = tprnfor
     ekin_conv_thr_ = ekin_conv_thr
     etot_conv_thr_ = etot_conv_thr
     forc_conv_thr_ = forc_conv_thr
     ekin_maxiter_  = electron_maxstep
     iesr_          = iesr
     remove_rigid_rot_ = remove_rigid_rot
     !
     ! ... define memory- and disk-related internal switches
     !
     smallmem = ( TRIM( memory ) == 'small' )
     twfcollect = wf_collect
     !
     ! Options for isolated system
     SELECT CASE( TRIM( assume_isolated ) )
         !
       CASE( 'makov-payne', 'm-p', 'mp' )
         !
         do_makov_payne = .TRUE.
         !
       CASE( 'none' )
         !
         do_makov_payne = .FALSE.
         !
       CASE DEFAULT
         !
         do_makov_payne = .FALSE.
         !
     END SELECT
     !
     tefield_  = tefield
     epol_     = epol
     efield_   = efield
     tefield2_ = tefield2
     epol2_    = epol2
     efield2_  = efield2
     !
     ! ... Set internal time step variables ( delt, twodelt, dt2 ... )
     !
     CALL set_time_step( dt )
     !
     ! ... Set electronic fictitius mass and its cut-off for fourier
     ! ... acceleration
     !
     emass_ = emass
     emaec_ = emass_cutoff
!====================================================================
!exx_wf related
     lwf = ( TRIM( calculation ) == 'cp-wf'      .OR. &
             TRIM( calculation ) == 'vc-cp-wf'   .OR. &
             TRIM( calculation ) == 'cp-wf-nscf')
     lwfnscf     = ( TRIM( calculation ) == 'cp-wf-nscf' )
     lwfpbe0nscf = ( dft_is_hybrid() .AND. lwfnscf  )
!====================================================================

     !
     ! ... set the level of output, the code verbosity 
     !
     trhor_ = ( TRIM( calculation ) == 'nscf'   .OR. &
                TRIM( calculation ) == 'cp-wf-nscf')
     trhow_ = saverho
     tksw_  = ( TRIM( disk_io ) == 'high' )
     !
     iverbosity_ = 0
     timing_ = .FALSE.
          ! The code write to files fort.8 fort.41 fort.42 fort.43
          ! a detailed report of subroutines timing
     memchk_ = .FALSE.
          ! The code performs a memory check, write on standard
          ! output the allocated memory at each step.
          ! Architecture Dependent
     !
     SELECT CASE( TRIM( verbosity ) )
       CASE( 'minimal' )
         !
         iverbosity_ =-1
         !
       CASE( 'low', 'default' )
         !
         iverbosity_ = 0
         timing_ = .TRUE.
         !
       CASE( 'medium' )
         !
         iverbosity_ = 1
         timing_   = .TRUE.
         !
       CASE( 'high' )
         !
         iverbosity_ = 2
         memchk_   = .TRUE.
         timing_   = .TRUE.
         !
       CASE( 'debug' )
         !
         iverbosity_ = 3
         memchk_   = .TRUE.
         timing_   = .TRUE.
         !
       CASE DEFAULT
         !
         CALL errore( 'control_flags ', &
                      'unknown verbosity ' // TRIM( verbosity ), 1 )
         !
     END SELECT
     !
     ! ... set the restart flags
     !
     trane_  = .FALSE.
     ampre_  = ampre
     taurdr_ = .FALSE.
     !
     SELECT CASE ( TRIM( restart_mode ) )
       !
       CASE( 'from_scratch' )
          !
          nbeg_   = -1
          nomore_ = nstep
          nstep_  = nstep
          trane_  = ( startingwfc == 'random' )
          !
          IF ( ampre_ == 0.D0 ) ampre_ = 0.02D0
          !
       CASE( 'reset_counters' )
          !
          nbeg_   =  0
          nomore_ = nstep
          nstep_  = nstep
          !
       CASE( 'restart' )
          !
          nbeg_   =  1
          nomore_ = nstep
          nstep_  = nstep
          nprint_nfi = -2
          !
       CASE( 'auto' )
          !
          IF ( auto_check( ndr, ' ' ) ) THEN
             !
             WRITE( stdout, '("autopilot: Auto Check detects restart.xml")' )
             WRITE( stdout, '("           adjusting restart_mode to restart")' )
             !
             restart_mode = 'restart'
             !
             nbeg_ = 1
             !
             ! ... Also handle NSTEPS adjustment so that
             ! ... nomore does not include past nfi in cpr.f90
             !
             restart_p = .TRUE.
             nomore_   = nstep
             nstep_    = nstep
             nprint_nfi = -2
             !
             IF ( ion_positions == 'from_input' ) THEN
                !
                taurdr_ = .TRUE.
                nbeg_   = -1
                !
             END IF
             !
          ELSE
             !
             WRITE( stdout, &
                    '("autopilot: Auto Check did not detect restart.xml")' )
             !
             WRITE( stdout, &
                    '("           adjusting restart_mode to from_scratch")' ) 
             !
             restart_mode = 'from_scratch'
             !
             nbeg_ = -2
             !
             IF ( ion_positions == 'from_input' ) nbeg_ = -1
             !
             nomore_ = nstep
             nstep_  = nstep
             !
             trane_ = ( startingwfc == 'random' )
             !
             IF ( ampre_ == 0.d0 ) ampre_ = 0.02D0
             !
          END IF
          !
       CASE DEFAULT
          !
          CALL errore( 'iosys ', &
                       'unknown restart_mode ' // TRIM( restart_mode ), 1 )
          !
     END SELECT
     !
     ! ... Starting/Restarting Atomic positions
     !
     SELECT CASE ( TRIM(ion_positions) )
       CASE ( 'from_input' )
         taurdr_ = .TRUE.   ! Positions read from standard input
       CASE ( 'default' )
         taurdr_ = .FALSE.
       CASE DEFAULT
         CALL errore(' control_flags ',' unknown ion_positions '//TRIM(ion_positions), 1 )
     END SELECT

     ! ... Electronic randomization
        
     SELECT CASE ( TRIM(startingwfc) )
       CASE ('default','none')
         trane_ = .FALSE.
       CASE ('random')
         trane_ = .TRUE.
       CASE DEFAULT
         CALL errore(' control_flags ',' unimplemented startingwfc='//TRIM(startingwfc), 1 )
     END SELECT
     IF( ampre_ == 0 ) trane_ = .FALSE.

      ! ...   TORTHO

      SELECT CASE ( orthogonalization )
      CASE ('Gram-Schmidt')
         tortho_ = .FALSE.
      CASE ('ortho')
         tortho_ = .TRUE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown orthogonalization '//&
              TRIM(orthogonalization), 1 )
      END SELECT

      ortho_max_  = ortho_max
      ortho_eps_  = ortho_eps

      ! ... Electrons initial velocity

      SELECT CASE ( TRIM(electron_velocities) )
        CASE ('default')
          tzeroe_ = .FALSE.
        CASE ('zero')
          tzeroe_ = .TRUE.
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown electron_velocities '//TRIM(electron_velocities), 1 )
      END SELECT

      ! ... Electron dynamics

      frice_ = 0.d0
      SELECT CASE ( TRIM(electron_dynamics) )
        CASE ('sd', 'default')
          tsde_ = .TRUE.
        CASE ('verlet')
          tsde_ = .FALSE.
        CASE ('cg')
          tsde_      = .FALSE.
          tcg = .TRUE.
          tortho_ = .FALSE.
        CASE ('damp')
          tsde_   = .FALSE.
          frice_ = electron_damping
        CASE ('diis')
          CALL errore( "iosys ", " electron_dynamics keyword diis not yet implemented ", 1 )
        CASE ('none')
          tsde_ = .FALSE.
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown electron_dynamics '//TRIM(electron_dynamics), 1 )
      END SELECT

      ! ... Electronic Temperature

      tnosee_ = .FALSE.
      SELECT CASE ( TRIM(electron_temperature) )
        !         temperature control of electrons via Nose' thermostat
        CASE ('nose')
          tnosee_ = .TRUE.
        CASE ('not_controlled', 'default')
          tnosee_ = .FALSE.
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown electron_temperature '//TRIM(electron_temperature), 1 )
      END SELECT
      
      ! ... Ions dynamics

      tconvthrs%active = .FALSE.
      tconvthrs%nstep  = 1
      tconvthrs%ekin   = 0.0d0
      tconvthrs%derho  = 0.0d0
      tconvthrs%force  = 0.0d0
      SELECT CASE ( TRIM(ion_dynamics) )
        CASE ('sd')
          tsdp_  = .TRUE.
          tfor_  = .TRUE.
          fricp_ = 0.d0
          tconvthrs%ekin   = ekin_conv_thr
          tconvthrs%derho  = etot_conv_thr
          tconvthrs%force  = forc_conv_thr
          tconvthrs%active = .TRUE.
          tconvthrs%nstep  = 1
        CASE ('verlet')
          tsdp_  = .FALSE.
          tfor_  = .TRUE.
          fricp_ = 0.d0
        CASE ('cg')       ! Conjugate Gradient minimization for ions
          CALL errore( 'iosys ', ' ion_dynamics = '//TRIM(ion_dynamics)//' not yet implemented ', 1 )
        CASE ('damp')
          tsdp_      = .FALSE.
          tfor_      = .TRUE.
          fricp_     = ion_damping
          tconvthrs%ekin   = ekin_conv_thr
          tconvthrs%derho  = etot_conv_thr
          tconvthrs%force  = forc_conv_thr
          tconvthrs%active = .TRUE.
          tconvthrs%nstep  = 1
        CASE ('none', 'default')
          tsdp_  = .FALSE.
          tfor_  = .FALSE.
          fricp_ = 0.d0
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown ion_dynamics '//TRIM(ion_dynamics), 1 )
      END SELECT

      ! External Forces on Ions has been specified 
      !
      IF ( ANY( rd_for(:,1:nat) /= 0.0_DP ) ) textfor = .TRUE.

      ! ... Ionic Temperature

      tcp_      = .FALSE.
      tnosep_   = .FALSE.
      tolp_     = tolp
      SELECT CASE ( TRIM(ion_temperature) )
        !         temperature control of ions via Nose' thermostat
        CASE ('nose')
          tnosep_ = .TRUE.
          tcp_ = .FALSE.
        CASE ('not_controlled', 'default')
          tnosep_ = .FALSE.
          tcp_ = .FALSE.
        CASE ('rescaling' )
          tnosep_ = .FALSE.
          tcp_ = .TRUE.
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown ion_temperature '//TRIM(ion_temperature), 1 )
      END SELECT

      ! ... Starting/Restarting ionic velocities

      tzerop_= .FALSE.
      tv0rd_ = .FALSE.
      tcap_  = .FALSE.
      SELECT CASE ( TRIM(ion_velocities) )
        CASE ('default')
          CONTINUE
        CASE ('change_step')
          dt_old_ = tolp
        CASE ('zero')
          tzerop_= .TRUE.
        CASE ('from_input')
          tv0rd_  = .TRUE.
          IF( .NOT. tavel  ) CALL errore(' iosys ', &
                               ' ION_VELOCITIES not present in stdin ', 1 )
        CASE ('random')
          tcap_ = .TRUE.
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown ion_velocities '//TRIM(ion_velocities), 1 )
      END SELECT

      ! ... Ionic randomization

      tranp_ ( 1 : ntyp ) =  tranp ( 1 : ntyp )
      amprp_ ( 1 : ntyp ) =  amprp ( 1 : ntyp )

      ! ... Ionic/electronic step ratio

      tionstep_ = .FALSE.
      nstepe_   = 1
      IF( ( ion_nstepe > 1 ) .OR. ( cell_nstepe > 1 ) ) THEN
        !         This card is used to control the ionic step, when active ionic step are
        !         allowed only when the two criteria are met, i.e. the ions are allowed
        !         to move if MOD( NFI, NSTEP ) == 0 and EKIN < EKIN_THR .
        tionstep_ = .TRUE.
        nstepe_   = MAX( ion_nstepe, cell_nstepe )
        WRITE(stdout, * ) "  ion_nstepe or cell_nstepe have no effects "
      END IF

      !   Cell dynamics
         
      SELECT CASE ( TRIM(cell_dynamics) )
        CASE ('sd')
          tpre_ = .TRUE.
          thdyn_ = .TRUE.
          tsdc_ = .TRUE.
          frich_= 0.d0
        CASE ( 'damp', 'damp-pr' )
          thdyn_ = .TRUE.
          tsdc_ = .FALSE.
          frich_ = cell_damping
          tpre_  = .TRUE.
        CASE ('pr')
          thdyn_ = .TRUE.
          tsdc_ = .FALSE.
          tpre_ = .TRUE.
          frich_= 0.d0
        CASE ('none', 'default')
          thdyn_ = .FALSE.
          tsdc_ = .FALSE.
          frich_= 0.d0
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown cell_dynamics '//TRIM(cell_dynamics), 1 )
      END SELECT

      ! ... Starting/Restarting Cell parameters

      SELECT CASE ( TRIM(cell_parameters) )
        CASE ('default')
          tbeg_ = .FALSE.
        CASE ('from_input')
          tbeg_ = .TRUE.
          IF( force_pairing_) &
            WRITE(stdout) " cell_parameters have no effects "
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown cell_parameters '//TRIM(cell_parameters), 1 )
      END SELECT

      ! ... Cell initial velocities

      SELECT CASE ( TRIM(cell_velocities) )
        CASE ('default')
          tzeroc_ = .FALSE.
        CASE ('zero')
          tzeroc_ = .TRUE.
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown cell_velocities '//TRIM(cell_velocities), 1 )
      END SELECT

      ! ... Cell Temperature

      SELECT CASE ( TRIM(cell_temperature) )
!         cell temperature control of ions via Nose' thermostat
        CASE ('nose')
          tnoseh_ = .TRUE.
        CASE ('not_controlled', 'default')
          tnoseh_ = .FALSE.
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown cell_temperature '//TRIM(cell_temperature), 1 )
      END SELECT

      ! .. If only electron are allowed to move 
      ! .. check for SCF convergence on the ground state
     
      IF( ion_dynamics == 'none' .AND. cell_dynamics == 'none' ) THEN
        tconvthrs%ekin   = ekin_conv_thr
        tconvthrs%derho  = etot_conv_thr
        tconvthrs%force  = 1.D+10
        tconvthrs%active = .TRUE.
        tconvthrs%nstep  = 1
      END IF

      ! force pairing

      force_pairing_ = force_pairing

      ! ... having set all input keywords, read plugins' input file(s)

      CALL plugin_read_input()

      !
      ! ... the 'ATOMIC_SPECIES' card must be present, check it

      IF( .NOT. taspc ) &
        CALL errore(' iosys ',' ATOMIC_SPECIES not found in stdin ',1)

      ! ... the 'ATOMIC_POSITIONS' card must be present, check it

      IF( .NOT. tapos ) &
        CALL errore(' iosys ',' ATOMIC_POSITIONS not found in stdin ',1)

      IF( .NOT. trd_ht .AND. TRIM(cell_parameters)=='from_input' ) &
        CALL errore(' iosys ',' CELL_PARAMETERS not present in stdin ', 1 )

      IF( .NOT. trd_ht .AND. ibrav == 0 ) &
        CALL errore(' iosys ',' ibrav = 0 but CELL_PARAMETERS not present in stdin ', 1 )

      RETURN
   END SUBROUTINE set_control_flags
   !
   !-------------------------------------------------------------------------
   SUBROUTINE modules_setup()
     !-------------------------------------------------------------------------
     !
     USE input_parameters, ONLY: ibrav , celldm , trd_ht, dt,                 &
           rd_ht, a, b, c, cosab, cosac, cosbc, ntyp , nat ,                  &
           na_inp , sp_pos , rd_pos , rd_vel, atom_mass, atom_label, if_pos,  &
           atomic_positions, sic, sic_epsilon, ecutwfc,                       &
           ecutrho, ecfixed, qcutz, q2sigma, tk_inp, wmass,                   &
           ion_radius, emass, emass_cutoff, temph, fnoseh, nr1b, nr2b, nr3b,  &
           tempw, fnosep, nr1, nr2, nr3, nr1s, nr2s, nr3s, ekincw, fnosee,    &
           outdir, prefix, nkstot, xk, vdw_table_name,                        &
           occupations, n_inner, fermi_energy, rotmass, occmass,              &
           rotation_damping, occupation_damping, occupation_dynamics,         &
           rotation_dynamics, degauss, smearing, nhpcl, nhptyp, ndega,        &
           nhgrp, fnhscl, cell_units, restart_mode, sic_alpha ,               &
           niter_cold_restart, lambda_cold, rd_for, ref_cell, rd_ref_ht,      &
           ref_cell_units, ref_alat

     USE input_parameters, ONLY: nconstr_inp, iprnks, nprnks,                  &
           etot_conv_thr, ekin_conv_thr, nspin, f_inp, nbnd,                   &
           press, cell_damping, cell_dofree, tf_inp,                           &
           refg, greash, grease, greasp, epol, efield, tcg, maxiter, conv_thr, &
           passop, tot_charge, tot_magnetization, niter_cg_restart
     !
     USE input_parameters, ONLY : wf_efield, wf_switch, sw_len, efx0, efy0,    &
                                  efz0, efx1, efy1, efz1, wfsd, wfdt, maxwfdt, &
                                  wf_q, wf_friction, nit, nsd, nsteps, tolw,   &
                                  adapt, calwf, nwf, wffort, writev,           &
                                  wannier_index
!===============================================================
!exx_wf related
     USE input_parameters, ONLY : neigh=>exx_neigh, vnbsp,&  
                                  poisson_eps=>exx_poisson_eps,&
                                  dis_cutoff=>exx_dis_cutoff,&
                                  exx_ps_rcut_s=>exx_ps_rcut_self,&
                                  exx_me_rcut_s=>exx_me_rcut_self,&
                                  exx_ps_rcut_p=>exx_ps_rcut_pair,&
                                  exx_me_rcut_p=>exx_me_rcut_pair
!===============================================================
     !
     USE input_parameters, ONLY : abivol, abisur, pvar, fill_vac,     &
                                  scale_at, t_gauss, jellium, cntr,   &
                                  P_ext, P_in, P_fin, rho_thr,        &
                                  step_rad, Surf_t, dthr, R_j, h_j,   &
                                  delta_eps, delta_sigma, n_cntr,     &
                                  axis
     USE input_parameters, ONLY : lda_plus_u, Hubbard_U
     USE input_parameters, ONLY : step_pen, A_pen, alpha_pen, sigma_pen
     USE input_parameters, ONLY : vdw_corr, london, london_s6, london_rcut, &
                                  ts_vdw, ts_vdw_isolated, ts_vdw_econv_thr
     !
     USE constants,        ONLY : amu_au, pi
     USE control_flags,    ONLY : lconstrain, tpre, thdyn, tksw
     USE ions_base,        ONLY : zv
     USE cell_base,        ONLY : cell_base_init, cell_dyn_init, at, cell_alat
     USE cell_base,        ONLY : ref_cell_base_init
     USE cell_nose,        ONLY : cell_nose_init
     USE ions_base,        ONLY : ions_base_init, greasp_ => greasp
     USE sic_module,       ONLY : sic_initval
     USE ions_nose,        ONLY : ions_nose_init
     USE wave_base,        ONLY : grease_ => grease
     USE electrons_nose,   ONLY : electrons_nose_init
     USE printout_base,    ONLY : printout_base_init
     USE efield_module,    ONLY : efield_init
     USE cg_module,        ONLY : cg_init
     USE pres_ai_mod,      ONLY : pres_ai_init
     USE ldaU_cp,          ONLY : ldaU_init0
     USE step_penalty,     ONLY : ldaUpen_init
     USE fft_base,         ONLY : dfftp, dffts, dfftb
     USE kohn_sham_states, ONLY : ks_states_init
     USE electrons_module, ONLY : electrons_setup
     USE electrons_base,   ONLY : electrons_base_initval
     USE ensemble_dft,     ONLY : ensemble_initval,tens
     USE wannier_base,     ONLY : wannier_init
     USE efield_module,    ONLY : tefield
     USE funct,            ONLY : dft_is_nonlocc, get_inlc
     USE kernel_table,     ONLY : vdw_table_name_ => vdw_table_name, &
                                  initialize_kernel_table
     USE control_flags,    ONLY : llondon, ts_vdw_ => ts_vdw
     USE london_module,    ONLY : init_london, scal6, lon_rcut
     USE tsvdw_module,     ONLY : vdw_isolated, vdw_econv_thr
     !
     IMPLICIT NONE
     !
     REAL(DP) :: alat_ , massa_totale
     ! ...   DIIS
     INTEGER :: ia, iss, inlc
     LOGICAL :: ltest
     !
     !   Subroutine Body
     !
     IF( .NOT. has_been_read ) &
       CALL errore( ' modules_setup ', ' input file has not been read yet! ', 1 )
     !
     ! ...  Set cell base module
     !
     massa_totale = SUM( atom_mass(1:ntyp)*na_inp(1:ntyp) )
     !
     CALL cell_base_init( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                          trd_ht, rd_ht, cell_units )
     CALL cell_dyn_init ( trd_ht, rd_ht, wmass, massa_totale, press, &
                          cell_damping, greash, cell_dofree )
     !
     alat_ = cell_alat()

     ! ...  Set ions base module

     CALL ions_base_init( ntyp , nat , na_inp , sp_pos , rd_pos , rd_vel,  &
                          atom_mass, atom_label, if_pos, atomic_positions, &
                          alat_ , at, ion_radius, rd_for )

     ! ...   Set Values for the cutoff

     CALL ecutoffs_setup( ecutwfc, ecutrho, ecfixed, qcutz, q2sigma, refg )
     !
     if (.not. allocated(xk)) then
       allocate(xk(3,1))
       xk = 0.d0
     endif
     !
     IF ( ref_cell ) THEN
       CALL ref_cell_base_init( ref_alat, rd_ref_ht, ref_cell_units )
       CALL gcutoffs_setup( ref_alat , tk_inp, nkstot, xk )
     ELSE
       CALL gcutoffs_setup( alat_ , tk_inp, nkstot, xk )
     END IF

     ! ... 
     
     grease_ = grease
     greasp_ = greasp
     !
     ! ... set thermostat parameter for cell, ions and electrons
     !
     CALL cell_nose_init( temph, fnoseh )
     !
     CALL ions_nose_init( tempw, fnosep, nhpcl, nhptyp, ndega, nhgrp, fnhscl)
     !
     CALL electrons_nose_init( ekincw , fnosee )

     ! set box grid module variables

     dfftb%nr1 = nr1b  
     dfftb%nr2 = nr2b
     dfftb%nr3 = nr3b

     ! set size for potentials and charge density
     ! (re-calculated automatically)

     dfftp%nr1  = nr1
     dfftp%nr2  = nr2
     dfftp%nr3  = nr3

     ! set size for wavefunctions
     ! (re-calculated automatically)

     dffts%nr1 = nr1s
     dffts%nr2 = nr2s
     dffts%nr3 = nr3s

     CALL efield_init( epol, efield )

     CALL cg_init( tcg , maxiter , conv_thr , passop ,niter_cg_restart)

     !
     IF( ( TRIM( sic ) /= 'none' ) .and. ( tpre .or. thdyn ) ) &
        CALL errore( ' module setup ', ' Stress is not yet implemented with SIC ', 1 )
     !
     CALL sic_initval( nat, sic, sic_epsilon, sic_alpha )
     !
     CALL ks_states_init( nspin, nprnks, iprnks )
     !
     !  kohn-sham states implies disk-io = 'high' 
     !
     DO iss = 1, nspin
        tksw = tksw .OR. ( nprnks(iss) > 0 )
     END DO

     CALL electrons_base_initval( zv, na_inp, ntyp, nbnd, nspin, &
                                  occupations, f_inp, &
                                  tot_charge, tot_magnetization )

     CALL electrons_setup( emass, emass_cutoff )

     CALL ensemble_initval( occupations, n_inner, fermi_energy,&
                            niter_cold_restart, lambda_cold, rotmass, &
                            occmass, rotation_damping, occupation_damping, &
                            occupation_dynamics, rotation_dynamics, degauss, &
                            smearing )
     IF( .NOT.tcg .AND. tens ) &
          CALL errore(' modules_setup ', 'Ensemble DFT implemented only with CG   ', 1 )
     !
     ! ... variables for constrained dynamics are set here
     !
     lconstrain = ( nconstr_inp > 0 )
     !
!========================================================================
!exx_wf related
     CALL wannier_init( wf_efield, wf_switch, sw_len, efx0, efy0, efz0, &
                        efx1, efy1, efz1, wfsd, wfdt, neigh, poisson_eps,&
                        dis_cutoff, exx_ps_rcut_s, exx_me_rcut_s,&
                        exx_ps_rcut_p, exx_me_rcut_p, vnbsp,&
                        maxwfdt, wf_q, &
                        wf_friction, nit, nsd, nsteps, tolw, adapt,     &
                        calwf, nwf, wffort, writev, wannier_index,      &
                        restart_mode )
!========================================================================
     !
     ! ... initialize variables for clusters under pressure 
     !
     CALL pres_ai_init( abivol, abisur, pvar, fill_vac, scale_at,       &
                        t_gauss, jellium, cntr, P_ext, P_in, P_fin,     &
                        rho_thr, step_rad, Surf_t, dthr, R_j, h_j,      &
                        delta_eps, delta_sigma, n_cntr, axis )
     !
     ! ... initialize variables for lda+U calculations
     !
     CALL ldaU_init0 ( ntyp, lda_plus_u, Hubbard_U )
     CALL ldaUpen_init( SIZE(sigma_pen), step_pen, sigma_pen, alpha_pen, A_pen )
     !
     !  ... initialize variables for vdW (dispersions) corrections
     !
     SELECT CASE( TRIM( vdw_corr ) )
         !
       CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
         !
         llondon= .TRUE.
         ts_vdw_= .FALSE.
         !
       CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
         !
         llondon= .FALSE.
         ts_vdw_= .TRUE.
         !
       CASE DEFAULT
         !
         llondon= .FALSE.
         ts_vdw_= .FALSE.
         !
     END SELECT
     ! 
     IF ( ts_vdw ) THEN
        CALL infomsg("iosys","ts_vdw is obsolete, use ''vdw_corr='ts-vdw''' instead")
        ts_vdw_ = .TRUE.
     END IF
     IF ( london ) THEN
        CALL infomsg("iosys","london is obsolete, use ''vdw_corr='grimme-d2''' instead")
        llondon = .TRUE.
     END IF
     IF (ts_vdw_.AND.llondon) CALL errore("iosys", &
                                    "must choose a unique vdW correction!", 1)
     IF ( llondon) THEN
        lon_rcut    = london_rcut
        scal6       = london_s6
        CALL init_london ( )
     ELSE IF ( ts_vdw_ ) THEN
        vdw_isolated = ts_vdw_isolated
        vdw_econv_thr= ts_vdw_econv_thr
     END IF
     !
     ! ... initialize kernel table for nonlocal functionals
     !
     IF ( dft_is_nonlocc( ) ) THEN
        vdw_table_name_ = vdw_table_name
        inlc = get_inlc()
        call initialize_kernel_table(inlc)
     ENDIF
     !
     RETURN
     !
  END SUBROUTINE modules_setup
  !
  !     --------------------------------------------------------
  !
  !     print out heading
  !
  SUBROUTINE input_info()

    ! this subroutine print to standard output some parameters read from input
    ! ----------------------------------------------

    USE input_parameters,   ONLY: restart_mode
    USE control_flags,      ONLY: nbeg, iprint, ndr, ndw, nomore
    USE time_step,          ONLY: delt
    USE cp_electronic_mass, ONLY: emass, emass_cutoff
    USE io_global,          ONLY: ionode, stdout

    IMPLICIT NONE

    IF( .NOT. has_been_read ) &
      CALL errore( ' iosys ', ' input file has not been read yet! ', 1 )

    IF( ionode ) THEN
      WRITE( stdout, 500) nbeg, restart_mode, nomore, iprint, ndr, ndw
      WRITE( stdout, 505) delt
      WRITE( stdout, 510) emass
      WRITE( stdout, 511) emass_cutoff
    END IF

500 FORMAT(   3X,'Restart Mode       = ',I7, 3X, A15, /, &
              3X,'Number of MD Steps = ',I7,  /, &
              3X,'Print out every      ',I7, ' MD Steps',/  &
              3X,'Reads from unit    = ',I7,  /, &
              3X,'Writes to unit     = ',I7)
505 FORMAT(   3X,'MD Simulation time step            = ',F10.2)
510 FORMAT(   3X,'Electronic fictitious mass (emass) = ',F10.2)
511 FORMAT(   3X,'emass cut-off                      = ',F10.2)

    RETURN
  END SUBROUTINE input_info
  !
  ! ----------------------------------------------------------------
  !
  SUBROUTINE modules_info()

    USE input_parameters, ONLY: electron_dynamics, electron_temperature, &
      orthogonalization

    USE control_flags, ONLY:  tortho, tnosee, trane, ampre, &
                              trhor, tksw, tfor, tnosep, iverbosity, &
                              thdyn, tnoseh
    !
    USE electrons_nose,       ONLY: electrons_nose_info
    USE sic_module,           ONLY: sic_info
    USE wave_base,            ONLY: frice, grease
    USE ions_base,            ONLY: fricp
    USE ions_nose,            ONLY: ions_nose_info
    USE cell_nose,            ONLY: cell_nose_info
    USE cell_base,            ONLY: frich
    USE efield_module,        ONLY: tefield, efield_info, tefield2, efield_info2
    USE io_global,            ONLY: ionode, stdout
    USE time_step,            ONLY: delt
    !
    !
    IMPLICIT NONE

    INTEGER :: is

    IF( .NOT. has_been_read ) &
      CALL errore( ' iosys ', ' input file has not been read yet! ', 1 )

    IF( ionode ) THEN
      !
      CALL cutoffs_print_info( )
      !
      IF( tortho ) THEN
        CALL orthogonalize_info( )
      ELSE
        WRITE( stdout,512)
      END IF
      !
      IF( TRIM(electron_dynamics) == 'sd' ) THEN
          WRITE( stdout,513)
      ELSE IF( TRIM(electron_dynamics) == 'verlet' ) THEN
          WRITE( stdout,510)
          frice = 0.d0
      ELSE IF( TRIM(electron_dynamics) == 'damp' ) THEN
          tnosee = .FALSE.
          WRITE( stdout,509)
          WRITE( stdout,514) frice, grease
      ELSE IF( TRIM(electron_dynamics) == 'cg' ) THEN
          WRITE( stdout,511)
      ELSE
          CALL errore(' input_info ', ' unknown electron dynamics ', 1 )
      END IF
      !
      IF( tnosee ) THEN
        WRITE( stdout,590)
        CALL electrons_nose_info(delt)
      ELSE 
        WRITE( stdout,535)
      END IF
      !
      IF( trane ) THEN
         WRITE( stdout,515) ampre
      ENDIF
      !
      CALL electrons_print_info( )
      !
      CALL exch_corr_print_info( )

      IF ( trhor ) THEN
         WRITE( stdout,720)
      ENDIF
      IF( tksw )THEN
         WRITE( stdout,722)
      ENDIF
      !
      IF( tfor .AND. tnosep ) fricp = 0.0d0
      !
      CALL ions_print_info( )
      !
      IF( tfor .AND. tnosep ) CALL ions_nose_info(delt)
      !
      CALL constraint_info( )
      !
      IF( thdyn .AND. tnoseh ) frich = 0.0d0
      !
      CALL cell_print_info( )
      !
      IF( thdyn .AND. tnoseh ) CALL cell_nose_info (delt)
      !
      !   CALL sic_info()  ! maybe useful
      !
      CALL plugin_print_info( )
      !
      IF(tefield) call efield_info( ) 
      IF(tefield2) call efield_info2( )

      WRITE( stdout,700) iverbosity

    END IF
    !
    RETURN
    !
509 FORMAT( 3X,'verlet algorithm for electron dynamics')
510 FORMAT( 3X,'Electron dynamics with newton equations')
511 FORMAT( 3X,'Electron dynamics with conjugate gradient')
512 FORMAT( 3X,'Orthog. with Gram-Schmidt')
513 FORMAT( 3X,'Electron dynamics with steepest descent')
514 FORMAT( 3X,'with friction frice = ',f7.4,' , grease = ',f7.4)
515 FORMAT( 3X,'initial random displacement of el. coordinates with ',   &
               ' amplitude=',f10.6)
535 FORMAT( 3X,'Electron dynamics : the temperature is not controlled')
590 FORMAT( 3X,'Electron temperature control via nose thermostat')
    !
700 FORMAT( /,3X, 'Verbosity: iverbosity = ',i2,/)
720 FORMAT(   3X, 'charge density is read from file')
722 FORMAT(   3X, 'Wavefunctions will be written to file as Kohn-Sham states')
    !
  END SUBROUTINE modules_info
  !
END MODULE input
