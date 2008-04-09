!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE input
   !---------------------------------------------------------------------------
   !
   USE kinds,     ONLY: DP
   !
   IMPLICIT NONE
   SAVE
   !
   PRIVATE                      !  Input Subroutines
                                !  should be called in the following order
                                !
   PUBLIC :: read_input_file    !  a) This sub. should be called first
   PUBLIC :: iosys_pseudo       !  b) then read pseudo files
   PUBLIC :: iosys              !  c) finally copy variables to modules
   PUBLIC :: modules_setup, set_control_flags
   !
   LOGICAL :: has_been_read = .FALSE.
   !
   CONTAINS
   !
   !-------------------------------------------------------------------------
   SUBROUTINE read_input_file()
     !-------------------------------------------------------------------------
     !
     USE read_namelists_module, ONLY : read_namelists
     USE read_cards_module,     ONLY : read_cards
     USE input_parameters,      ONLY : calculation, title
     USE control_flags,         ONLY : lneb, lpath, lwf, lmetadyn, &
                                       program_name
     USE printout_base,         ONLY : title_ => title
     USE io_global,             ONLY : meta_ionode, stdout
     USE xml_input,             ONLY : xml_input_dump
     !
     IMPLICIT NONE
     !
     CHARACTER(LEN=2) :: prog
     !
     !
     prog = 'CP'
     !
     IF ( meta_ionode ) THEN
        CALL xml_input_dump()
        CALL input_from_file()
     END IF
     !
     ! ... Read NAMELISTS 
     !
     CALL read_namelists( prog )
     !
     ! ... Read CARDS 
     !
     CALL read_cards( prog )
     !
     IF( TRIM( calculation ) == 'fpmd' ) program_name = 'FPMD'
     !
     IF( TRIM( calculation ) == 'fpmd-neb' ) THEN
        !
        program_name = 'FPMD'
        lneb = .TRUE.
        !
     ELSE
        !
        lneb = ( TRIM( calculation ) == 'neb' )
        !
     END IF
     !
     lpath = lneb
     !
     lmetadyn = ( TRIM( calculation ) == 'metadyn' )
     !
     lwf = ( TRIM( calculation ) == 'cp-wf' )
     !
     ! ... Set job title and print it on standard output
     !
     title_ = title
     !
     WRITE( stdout, '(/,3X,"Job Title: ",A )' ) TRIM( title_ )
     !
     has_been_read = .TRUE.
     !
     RETURN
     !
   END SUBROUTINE read_input_file
   !
   !-------------------------------------------------------------------------
   SUBROUTINE iosys_pseudo()
     !-------------------------------------------------------------------------
     !
     USE input_parameters,        ONLY : atom_pfile, pseudo_dir, ntyp, nat, &
                                         prefix, outdir, xc_type
     USE control_flags,           ONLY : program_name
     USE parameters,              ONLY : nsx
     USE read_pseudo_module_fpmd, ONLY : readpp
     USE io_files,                ONLY : psfile_     => psfile , &
                                         pseudo_dir_ => pseudo_dir, &
                                         outdir_     => outdir, &
                                         prefix_     => prefix
     USE ions_base,               ONLY : nsp_ => nsp, nat_ => nat
     !
     IMPLICIT NONE
     !
     !
     IF ( .NOT. has_been_read ) &
        CALL errore( 'iosys_pseudo ', 'input file has not been read yet!', 1 )
     !
     prefix_  = TRIM( prefix  )
     outdir_  = TRIM( outdir )
     !
     ! ... Set internal variables for the number of species and number of atoms
     !
     nsp_ = ntyp
     nat_ = nat
     !
     psfile_         = ' '
     psfile_(1:nsp_) = atom_pfile(1:nsp_)
     pseudo_dir_     = TRIM( pseudo_dir  )
     !
     ! ... read in pseudopotentials and wavefunctions files
     !
     CALL readpp( xc_type )
     !
     RETURN
     !
   END SUBROUTINE iosys_pseudo
   !
   !-------------------------------------------------------------------------
   SUBROUTINE iosys()
     !-------------------------------------------------------------------------
     !
     USE control_flags,      ONLY : fix_dependencies, program_name, &
                                    lconstrain, lmetadyn
     USE io_global,          ONLY : meta_ionode, stdout
     USE ions_base,          ONLY : nat, tau, ityp
     USE constraints_module, ONLY : init_constraint
     USE metadyn_vars,       ONLY : init_metadyn_vars     
     !
     IMPLICIT NONE
     !
     !
     IF ( meta_ionode ) THEN
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
     IF ( lmetadyn ) CALL init_metadyn_vars()
     !
     ! ... fix values for dependencies
     !
     IF ( program_name == 'FPMD' ) THEN 
        !
        IF ( lconstrain .OR. lmetadyn ) THEN

           ! ...  Apply sort to constraints atomic index

           CALL new_atomind_constraints()

        END IF 
 
        CALL fix_dependencies()
        !
     END IF
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
     USE control_flags, ONLY : lcoarsegrained, ldamped, lmetadyn
     USE control_flags, ONLY : program_name
     USE control_flags, ONLY : ndw_        => ndw, &
                               ndr_        => ndr, &
                               iprint_     => iprint, &
                               isave_      => isave, &
                               tstress_    => tstress, &
                               tprnfor_    => tprnfor, &
                               tprnsfac_   => tprnsfac, &
                               ampre_      => ampre, &
                               trane_      => trane, &
                               newnfi_     => newnfi, &
                               tnewnfi_    => tnewnfi, &
                               tdipole_    => tdipole, &
                               nomore_     => nomore, &
                               memchk_     => memchk, &
                               tpre_       => tpre, &
                               timing_     => timing, &
                               iprsta_     => iprsta, &
                               taurdr_     => taurdr, &
                               nbeg_       => nbeg, &
                               gamma_only_ => gamma_only, &
                               tchi2_      => tchi2, &
                               tatomicwfc_ => tatomicwfc, &
                               printwfc_   => printwfc, &
                               tortho_     => tortho,   &
                               nstep_      => nstep
     USE control_flags, ONLY : tsde_          => tsde, &
                               tsteepdesc_    => tsteepdesc, &
                               tzeroe_        => tzeroe, &
                               tdamp_         => tdamp, &
                               trhor_         => trhor, &
                               trhow_         => trhow, &
                               tksw_          => tksw,  &
                               ortho_eps_     => ortho_eps, &
                               ortho_max_     => ortho_max, &
                               tnosee_        => tnosee
     USE control_flags, ONLY : tdampions_ => tdampions, &
                               tfor_      => tfor, &
                               tsdp_      => tsdp, &
                               lfixatom, tconvthrs
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
     USE control_flags, ONLY : iesr, tvhmean, vhrmin, vhrmax, vhasse
     USE control_flags, ONLY : tprojwfc
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
     USE cvan,               ONLY : nvb
     USE control_flags,      ONLY : ortho_para_ => ortho_para
     !
     USE input_parameters,   ONLY: &
        electron_dynamics, electron_damping, electron_temperature,   &
        ion_dynamics, ekin_conv_thr, etot_conv_thr, forc_conv_thr, ion_maxstep,&
        electron_maxstep, ion_damping, ion_temperature, ion_velocities, tranp, &
        amprp, ion_nstepe, cell_nstepe, cell_dynamics, cell_damping,           &
        cell_parameters, cell_velocities, cell_temperature, force_pairing,     &
        tapos, tavel, ecutwfc, emass, emass_cutoff, taspc, trd_ht, ibrav,      &
        ortho_eps, ortho_max, ntyp, tolp, tchi2_inp, calculation, disk_io, dt, &
        tcg, ndr, ndw, iprint, isave, tstress, k_points, tprnfor, verbosity,   &
        tdipole_card, tnewnfi_card, newnfi_card,                               &
        ampre, nstep, restart_mode, ion_positions, startingwfc, printwfc,      &
        orthogonalization, electron_velocities, nat, if_pos, phase_space,      &
        tefield, epol, efield, tefield2, epol2, efield2, remove_rigid_rot,     &
        iesr_inp, vhrmax_inp, vhrmin_inp, tvhmean_inp, vhasse_inp, saverho,    &
        ortho_para
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
     printwfc_      = printwfc
     tchi2_         = tchi2_inp
     ekin_conv_thr_ = ekin_conv_thr
     etot_conv_thr_ = etot_conv_thr
     forc_conv_thr_ = forc_conv_thr
     ekin_maxiter_  = electron_maxstep
     iesr           = iesr_inp
     !
     tvhmean = tvhmean_inp
     vhrmin  = vhrmin_inp
     vhrmax  = vhrmax_inp
     vhasse  = vhasse_inp
     !
     remove_rigid_rot_ = remove_rigid_rot
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
     !
     ! ... set the level of output, the code verbosity 
     !
     iprsta_ = 1
     timing_ = .FALSE.
          ! The code write to files fort.8 fort.41 fort.42 fort.43
          ! a detailed report of subroutines timing
     memchk_ = .FALSE.
          ! The code performs a memory check, write on standard
          ! output the allocated memory at each step.
          ! Architecture Dependent
     tprnsfac_  = .FALSE.
          ! Print on file STRUCTURE_FACTOR the structure factor
          ! gvectors and charge density, in reciprocal space.
     !
     trhor_ = ( TRIM( calculation ) == 'nscf' )
     trhow_ = saverho
     tksw_  = ( TRIM( disk_io ) == 'high' )
     !
     SELECT CASE( TRIM( verbosity ) )
       CASE( 'minimal' )
         !
         iprsta_ = 0
         !
       CASE( 'low', 'default' )
         !
         iprsta_ = 1
         timing_ = .TRUE.
         !
       CASE( 'default+projwfc' )
         !
         iprsta_  = 1
         timing_  = .TRUE.
         tprojwfc = .TRUE.
         !
       CASE( 'medium' )
         !
         iprsta_   = 2
         timing_   = .TRUE.
         tprnsfac_ = .TRUE.
         !
       CASE( 'high' )
         !
         iprsta_   = 3
         memchk_   = .TRUE.
         timing_   = .TRUE.
         tprnsfac_ = .TRUE.
         tprojwfc  = .TRUE.
         !
       CASE DEFAULT
         !
         CALL errore( 'control_flags ', &
                      'unknown verbosity ' // TRIM( verbosity ), 1 )
         !
     END SELECT
     !
     tdipole_  = tdipole_card
     newnfi_   = newnfi_card
     tnewnfi_  = tnewnfi_card
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
        
     tatomicwfc_ = .FALSE.
     SELECT CASE ( TRIM(startingwfc) )
       CASE ('default','none')
         trane_ = .FALSE.
       CASE ('random')
         trane_ = .TRUE.
       CASE ('atomic')
         tatomicwfc_ = .TRUE.
       CASE DEFAULT
         CALL errore(' control_flags ',' unknown startingwfc '//TRIM(startingwfc), 1 )
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
      ortho_para_ = ortho_para

      ! ... Electrons initial velocity

      SELECT CASE ( TRIM(electron_velocities) )
        CASE ('default')
          tzeroe_ = .FALSE.
        CASE ('zero')
          tzeroe_ = .TRUE.
          IF( program_name == 'CP90' ) &
            WRITE( stdout, &
                   '("Warning: electron_velocities keyword has no effect")' )
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown electron_velocities '//TRIM(electron_velocities), 1 )
      END SELECT

      ! ... Electron dynamics

      tdamp_          = .FALSE.
      tsteepdesc_     = .FALSE.
      frice_ = 0.d0
      SELECT CASE ( TRIM(electron_dynamics) )
        CASE ('sd', 'default')
          tsde_ = .TRUE.
        CASE ('verlet')
          tsde_ = .FALSE.
        CASE ('cg')
          tsde_      = .FALSE.
          IF( program_name == 'CP90' ) THEN
             tcg = .TRUE.
             tortho_ = .FALSE.
          ELSE
             CALL errore(' control_flags ',' conjugate gradient not yet implemented in FPMD ', 1 )
          ENDIF
        CASE ('damp')
          tsde_   = .FALSE.
          tdamp_  = .TRUE.
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
      IF ( lmetadyn ) lcoarsegrained  = .TRUE.

      ! ... Ions dynamics

      tdampions_       = .FALSE.
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
          CALL errore( "iosys ", " ion_dynamics = '//TRIM(ion_dynamics)//' not yet implemented ", 1 )
        CASE ('damp')
          ldamped    = .TRUE.
          tsdp_      = .FALSE.
          tfor_      = .TRUE.
          tdampions_ = .TRUE.
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

      

      IF ( ANY( if_pos(:,1:nat) == 0 ) ) lfixatom = .TRUE.

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

      tcap_         = .FALSE.
      SELECT CASE ( TRIM(ion_velocities) )
        CASE ('default')
          tzerop_ = .FALSE.
          tv0rd_ = .FALSE.
          tcap_ = .FALSE.
        CASE ('change_step')
          tzerop_ = .FALSE.
          tv0rd_ = .FALSE.
          tcap_ = .FALSE.
          dt_old_ = tolp
        CASE ('zero')
          tzerop_ = .TRUE.
          tv0rd_ = .FALSE.
        CASE ('from_input')
          tzerop_ = .TRUE.
          tv0rd_  = .TRUE.
        CASE ('random')
          tcap_ = .TRUE.
          IF( program_name == 'FPMD' ) &
            WRITE(stdout) " ion_velocities = '//TRIM(ion_velocities)//' has no effects "
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
        IF( program_name == 'CP90' ) &
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
          IF( program_name == 'CP90' .AND. force_pairing_) &
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
      !
      IF( ( nvb > 0 ) .and. ( program_name == 'FPMD' ) ) &
        CALL errore(' iosys ',' USPP not yet implemented in FPMD ',1)

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

      IF( .NOT. tavel .AND. TRIM(ion_velocities)=='from_input' ) &
        CALL errore(' iosys ',' ION_VELOCITIES not present in stdin ', 1 )

      IF( ( TRIM( calculation ) == 'smd' ) .AND. ( TRIM( cell_dynamics ) /= 'none' ) ) THEN
        CALL errore(' smiosys ',' cell_dynamics not implemented : '//TRIM(cell_dynamics), 1 )
      END IF

      RETURN
   END SUBROUTINE set_control_flags
   !
   !-------------------------------------------------------------------------
   SUBROUTINE modules_setup()
     !-------------------------------------------------------------------------
     !
     USE control_flags,    ONLY : program_name, lconstrain, lneb, lmetadyn, &
                                  tpre, thdyn, tksw

     USE constants,        ONLY : amu_au, pi
     !
     USE input_parameters, ONLY: ibrav , celldm , trd_ht, dt,    &
           cell_symmetry, rd_ht, a, b, c, cosab, cosac, cosbc, ntyp , nat ,   &
           na_inp , sp_pos , rd_pos , rd_vel, atom_mass, atom_label, if_pos,  &
           atomic_positions, id_loc, sic, sic_epsilon, sic_rloc, ecutwfc,     &
           ecutrho, ecfixed, qcutz, q2sigma, tk_inp, wmass,                   &
           ion_radius, emass, emass_cutoff, temph, fnoseh, nr1b, nr2b, nr3b,  &
           tempw, fnosep, nr1, nr2, nr3, nr1s, nr2s, nr3s, ekincw, fnosee,    &
           tturbo_inp, nturbo_inp, outdir, prefix,                            &
           k_points, nkstot, nk1, nk2, nk3, k1, k2, k3,                       &
           xk, wk, occupations, n_inner, fermi_energy, rotmass, occmass,      &
           rotation_damping, occupation_damping, occupation_dynamics,         &
           rotation_dynamics, degauss, smearing, nhpcl, nhptyp, ndega,        &
           nhgrp, fnhscl, cell_units, restart_mode, sic_alpha ,               &
           niter_cold_restart, lambda_cold

     USE input_parameters, ONLY: empty_states_maxstep,                         &
           empty_states_ethr, empty_states_nbnd,                               &
           iprnks_empty, nconstr_inp, iprnks, nprnks,                          &
           etot_conv_thr, ekin_conv_thr, nspin, f_inp, nelup, neldw, nbnd,     &
           nelec, press, cell_damping, cell_dofree, tf_inp, nprnks_empty,      &
           refg, greash, grease, greasp, epol, efield, tcg, maxiter, conv_thr, &
           passop, tot_charge, multiplicity, tot_magnetization, ncolvar_inp,   &
           niter_cg_restart
     !
     USE input_parameters, ONLY : wf_efield, wf_switch, sw_len, efx0, efy0,    &
                                  efz0, efx1, efy1, efz1, wfsd, wfdt, maxwfdt, &
                                  wf_q, wf_friction, nit, nsd, nsteps, tolw,   &
                                  adapt, calwf, nwf, wffort, writev,           &
                                  wannier_index
     !
     USE input_parameters, ONLY : abivol, abisur, pvar, fill_vac,     &
                                  scale_at, t_gauss, jellium, cntr,   &
                                  P_ext, P_in, P_fin, rho_thr,        &
                                  step_rad, Surf_t, dthr, R_j, h_j,   &
                                  delta_eps, delta_sigma, n_cntr,     &
                                  axis
     !
     USE ions_base,        ONLY : tau, ityp, zv
     USE cell_base,        ONLY : cell_base_init, a1, a2, a3, cell_alat
     USE cell_nose,        ONLY : cell_nose_init
     USE ions_base,        ONLY : ions_base_init, greasp_ => greasp
     USE sic_module,       ONLY : sic_initval
     USE ions_nose,        ONLY : ions_nose_init
     USE wave_base,        ONLY : grease_ => grease
     USE electrons_nose,   ONLY : electrons_nose_init
     USE printout_base,    ONLY : printout_base_init
     USE turbo,            ONLY : turbo_init
     USE efield_module,    ONLY : efield_init
     USE cg_module,        ONLY : cg_init
     USE pres_ai_mod,      ONLY : pres_ai_init
     !
     USE smallbox_grid_dimensions, ONLY: &
           nnrbx, &  !  variable is used to workaround internal compiler error (IBM xlf)
           nr1b_ => nr1b, &
           nr2b_ => nr2b, &
           nr3b_ => nr3b
     USE grid_dimensions,          ONLY: &
           nnrx, &  !  variable is used to workaround internal compiler error (IBM xlf)
           nr1_ => nr1, &
           nr2_ => nr2, &
           nr3_ => nr3
     USE smooth_grid_dimensions,   ONLY: &
           nnrsx, &  !  variable is used to workaround internal compiler error (IBM xlf)
           nr1s_ => nr1s, &
           nr2s_ => nr2s, &
           nr3s_ => nr3s
     USE charge_mix,         ONLY : charge_mix_setup
     USE kohn_sham_states,   ONLY : ks_states_init
     USE electrons_module,   ONLY : electrons_setup, empty_init
     USE electrons_base,     ONLY : electrons_base_initval
     USE ensemble_dft,       ONLY : ensemble_initval,tens
     USE wannier_base,       ONLY : wannier_init
     USE efield_module,      ONLY : tefield
     !
     IMPLICIT NONE
     !
     REAL(DP) :: alat_ , massa_totale
     REAL(DP) :: ethr_emp_inp
     ! ...   DIIS
     INTEGER :: ia, iss
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
     CALL cell_base_init( ibrav , celldm , trd_ht, cell_symmetry, rd_ht, &
                          cell_units, a, b, c, cosab, cosac, cosbc , wmass, &
                          massa_totale, press, cell_damping, greash, &
                          cell_dofree )
     !
     alat_ = cell_alat()

     ! ...  Set ions base module

     CALL ions_base_init( ntyp , nat , na_inp , sp_pos , rd_pos , rd_vel,  &
                          atom_mass, atom_label, if_pos, atomic_positions, &
                          alat_ , a1, a2, a3, ion_radius )

     ! ...   Set Values for the cutoff

     CALL ecutoffs_setup( ecutwfc, ecutrho, ecfixed, qcutz, q2sigma, refg )

     CALL gcutoffs_setup( alat_ , tk_inp, nkstot, xk )

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

     nr1b_ = nr1b  
     nr2b_ = nr2b
     nr3b_ = nr3b

     ! set size for potentials and charge density
     ! (re-calculated automatically)

     nr1_  = nr1
     nr2_  = nr2
     nr3_  = nr3

     ! set size for wavefunctions
     ! (re-calculated automatically)

     nr1s_ = nr1s
     nr2s_ = nr2s
     nr3s_ = nr3s

     CALL turbo_init( tturbo_inp, nturbo_inp )

     IF ( .NOT. lneb ) &
        CALL printout_base_init( outdir, prefix )

     CALL efield_init( epol, efield )

     CALL cg_init( tcg , maxiter , conv_thr , passop ,niter_cg_restart)

     !
     IF( ( TRIM( sic ) /= 'none' ) .and. ( tpre .or. thdyn ) ) &
        CALL errore( ' module setup ', ' Stress is not yet implemented with SIC ', 1 )
     !
     CALL sic_initval( nat, id_loc, sic, sic_epsilon, sic_alpha, sic_rloc  )

     !
     !  empty states
     !
     ethr_emp_inp  = ekin_conv_thr
     IF( empty_states_ethr > 0.d0 )  ethr_emp_inp  = empty_states_ethr
     CALL empty_init( empty_states_maxstep, ethr_emp_inp )

     !
     CALL ks_states_init( nspin, nprnks, iprnks, nprnks_empty, iprnks_empty )
     !
     !  kohn-sham states implies disk-io = 'high' 
     !
     DO iss = 1, nspin
        tksw = tksw .OR. ( nprnks(iss) > 0 )
        tksw = tksw .OR. ( nprnks_empty(iss) > 0 )
     END DO

     CALL electrons_base_initval( zv, na_inp, ntyp, nelec, nelup,         &
                                  neldw, nbnd, nspin, occupations, f_inp, &
                                  tot_charge, multiplicity, tot_magnetization )

     CALL electrons_setup( empty_states_nbnd, emass, emass_cutoff )

     CALL ensemble_initval( occupations, n_inner, fermi_energy,&
                            niter_cold_restart, lambda_cold, rotmass, &
                            occmass, rotation_damping, occupation_damping, &
                            occupation_dynamics, rotation_dynamics, degauss, &
                            smearing )
     IF( ( program_name == 'CP90' ) .AND. .NOT.tcg .AND. tens ) &
          CALL errore(' modules_setup ', 'Ensemble DFT implemented only with CG   ', 1 )
     !
     ! ... variables for constrained dynamics are set here
     !
     lconstrain = ( ncolvar_inp + nconstr_inp > 0 )
     !
     !
     CALL wannier_init( wf_efield, wf_switch, sw_len, efx0, efy0, efz0, &
                        efx1, efy1, efz1, wfsd, wfdt, maxwfdt, wf_q,    &
                        wf_friction, nit, nsd, nsteps, tolw, adapt,     &
                        calwf, nwf, wffort, writev, wannier_index,      &
                        restart_mode )
     !
     ! ... initialize variables for clusters under pressure 
     !
     CALL pres_ai_init( abivol, abisur, pvar, fill_vac, scale_at,       &
                        t_gauss, jellium, cntr, P_ext, P_in, P_fin,     &
                        rho_thr, step_rad, Surf_t, dthr, R_j, h_j,      &
                        delta_eps, delta_sigma, n_cntr, axis )
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
    USE io_global,          ONLY: meta_ionode, stdout

    IMPLICIT NONE

    IF( .NOT. has_been_read ) &
      CALL errore( ' iosys ', ' input file has not been read yet! ', 1 )

    IF( meta_ionode ) THEN
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
509 FORMAT(   3X,'Verlet algorithm for electron dynamics')
502 FORMAT(   3X,'An initial quench is performed')

    RETURN
  END SUBROUTINE input_info
  !
  ! ----------------------------------------------------------------
  !
  SUBROUTINE modules_info()

    USE input_parameters, ONLY: electron_dynamics, electron_temperature, &
      orthogonalization

    USE control_flags, ONLY:  program_name, tortho, tnosee, trane, ampre, &
                              trhor, tksw, tfor, tnosep, iprsta, &
                              thdyn, tnoseh
    !
    USE electrons_nose,       ONLY: electrons_nose_info
    USE electrons_module,     ONLY: empty_print_info
    USE sic_module,           ONLY: sic_info
    USE wave_base,            ONLY: frice, grease
    USE ions_base,            ONLY: fricp
    USE ions_nose,            ONLY: ions_nose_info
    USE cell_nose,            ONLY: cell_nose_info
    USE cell_base,            ONLY: frich
    USE efield_module,        ONLY: tefield, efield_info, tefield2, efield_info2
    USE io_global,            ONLY: meta_ionode, stdout
    !
    !
    IMPLICIT NONE

    INTEGER :: is

    IF( .NOT. has_been_read ) &
      CALL errore( ' iosys ', ' input file has not been read yet! ', 1 )

    IF( meta_ionode ) THEN
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
        CALL electrons_nose_info()
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
      IF( program_name == 'FPMD' ) THEN
        CALL empty_print_info( stdout )
      END IF
      !
      IF( tfor .AND. tnosep ) fricp = 0.0d0
      !
      CALL ions_print_info( )
      !
      IF( tfor .AND. tnosep ) CALL ions_nose_info()
      !
      CALL constraint_info( )
      !
      IF( thdyn .AND. tnoseh ) frich = 0.0d0
      !
      CALL cell_print_info( )
      !
      IF( thdyn .AND. tnoseh ) CALL cell_nose_info()
      !
      IF ( program_name == 'FPMD' ) THEN
         !
         CALL potential_print_info( stdout )
         CALL sic_info()
         !
      END IF
      !
      IF(tefield) call efield_info( ) 
      IF(tefield2) call efield_info2( )

      WRITE( stdout,700) iprsta

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
700 FORMAT( /,3X, 'Verbosity: iprsta = ',i2,/)
720 FORMAT(   3X, 'charge density is read from file')
722 FORMAT(   3X, 'Wavefunctions will be written to file as Kohn-Sham states')
    !
  END SUBROUTINE modules_info
  !
END MODULE input
