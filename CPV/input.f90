!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


! ----------------------------------------------------------------
  MODULE input
! ----------------------------------------------------------------

   USE kinds, ONLY: dbl
   USE io_global, ONLY: ionode, stdout

   IMPLICIT NONE
   SAVE

   PRIVATE

   PUBLIC :: read_input_file
   PUBLIC :: iosys_pseudo
   PUBLIC :: iosys

   LOGICAL :: has_been_read = .FALSE.

! ----------------------------------------------------------------
  CONTAINS
! ----------------------------------------------------------------


   SUBROUTINE read_input_file( lneb, lsmd, lwf )
        USE read_namelists_module, ONLY: read_namelists
        USE read_cards_module, ONLY: read_cards
        USE input_parameters, ONLY: calculation
        USE control_flags, ONLY: program_name
        IMPLICIT NONE

        LOGICAL, OPTIONAL, INTENT(OUT) :: lneb, lsmd, lwf
        CHARACTER(LEN=2) :: prog

        IF( program_name == 'FPMD' ) prog = 'FP'
        IF( program_name == 'CP90' ) prog = 'CP'

        ! . Read NAMELISTS ..................................................!

        CALL read_namelists( prog )

        ! . Read CARDS ......................................................!

        CALL read_cards( prog )

        IF( PRESENT(lneb) ) THEN
          lneb = ( TRIM( calculation ) == 'neb' )
        END IF
        IF( PRESENT(lsmd) ) THEN
          lsmd = ( TRIM( calculation ) == 'smd' )
        END IF
        IF( PRESENT(lwf) ) THEN
          lwf  = ( TRIM( calculation ) == 'cp-wf' )
        END IF

         has_been_read = .TRUE.

        RETURN
   END SUBROUTINE

   !  ----------------------------------------------

   SUBROUTINE iosys_pseudo( )

        use input_parameters, only:  atom_pfile, pseudo_dir, ntyp, prefix, scradir
        use control_flags, only:  program_name
        use parameters, only: nsx
        use read_pseudo_module_fpmd, only: readpp
        USE io_files, ONLY: &
              psfile_ => psfile , &
              pseudo_dir_ => pseudo_dir, &
              scradir_ => scradir, &
              prefix_ => prefix
        USE ions_base, ONLY: nsp_ => nsp

        implicit none

        IF( .NOT. has_been_read ) &
          CALL errore( ' iosys_pseudo ', ' input file has not been read yet! ', 1 )

        prefix_  = TRIM( prefix  )
        scradir_ = TRIM( scradir )

        nsp_   = ntyp
        psfile_= ' '
        psfile_ ( 1:nsp_ ) = atom_pfile( 1:nsp_ )
        pseudo_dir_        = TRIM( pseudo_dir  )
        !
        !  read in pseudopotentials and wavefunctions files
        !
        call readpp( )
        !
        return
   end subroutine

   !  ----------------------------------------------

   SUBROUTINE iosys

     USE control_flags, ONLY: fix_dependencies
     USE control_flags, only:  program_name

     IMPLICIT NONE

        ! . Set internal flags according to the input .......................!

        CALL set_control_flags( )
        CALL modules_setup()

        ! . CALL the Module specific setup routine ..........................!

        CALL modules_setup_cp()

        ! . Fix values for dependencies .....................................!

        IF( program_name == 'FPMD' ) THEN
          CALL fix_dependencies()
        END IF


        ! . Write to stdout input module information ........................!

        IF( program_name == 'FPMD' ) THEN
          CALL input_info()
          CALL modules_info()
        ELSE
          CALL input_info_cp()
        END IF

     RETURN
   END SUBROUTINE iosys

   !  ----------------------------------------------

   SUBROUTINE set_control_flags( )
     !
     USE kinds, ONLY: dbl
     USE control_flags, ONLY: program_name
     USE io_global, ONLY: stdout


     USE control_flags, ONLY: &
        ndw_     => ndw, &
        ndr_     => ndr, &
        iprint_  => iprint, &
        isave_   => isave, &
        tstress_ => tstress, &
        tprnfor_ => tprnfor, &
        tprnsfac_ => tprnsfac, &
        toptical_ => toptical, &
        ampre_   => ampre, &
        trane_   => trane, &
        newnfi_  => newnfi, &
        tnewnfi_ => tnewnfi, &
        rhoout_  => rhoout, &
        tdipole_ => tdipole, &
        nomore_  => nomore, &
        memchk_  => memchk, &
        tpre_    => tpre, &
        prn_     => prn, &
        timing_  => timing, &
        iprsta_  => iprsta, &
        taurdr_  => taurdr, &
        nbeg_    => nbeg, &
        gamma_only_ => gamma_only, &
        tchi2_   => tchi2, &
        tatomicwfc_ => tatomicwfc, &
        printwfc_ => printwfc, &
        tortho_  => tortho

     USE control_flags, ONLY: &
        t_diis_simple_ => t_diis_simple, &
        t_diis_ => t_diis, &
        tsde_ => tsde, &
        t_diis_rot_ => t_diis_rot, &
        tconjgrad_ => tconjgrad, &
        tsteepdesc_ => tsteepdesc, &
        tzeroe_ => tzeroe, &
        tdamp_ => tdamp, &
        trhor_ => trhor, &
        trhow_ => trhow, &
        tvlocw_ => tvlocw, &
        ortho_eps_ => ortho_eps, &
        ortho_max_ => ortho_max, &
        tnosee_ => tnosee

     USE control_flags, ONLY: &
        tdampions_ => tdampions, &
        tfor_ => tfor, &
        tsdp_ => tsdp, &
        tconvthrs, tconjgrad_ion

     USE control_flags, ONLY: &
        tnosep_ => tnosep, &
        tcap_ => tcap, &
        tcp_ => tcp, &
        tolp_ => tolp, &
        tzerop_ => tzerop, &
        tv0rd_ => tv0rd, &
        tranp_ => tranp, &
        amprp_ => amprp

      USE control_flags, ONLY: &
        tionstep_ => tionstep, &
        nstepe_ => nstepe

      USE control_flags, ONLY: &
        tzeroc_ => tzeroc, &
        tnoseh_ => tnoseh, &
        thdyn_ => thdyn, &
        tsdc_ => tsdc, &
        tbeg_ => tbeg

      USE control_flags, ONLY: &
          ekin_conv_thr_ => ekin_conv_thr, &
          etot_conv_thr_ => etot_conv_thr, &
          forc_conv_thr_ => forc_conv_thr, &
          ekin_maxiter_ => ekin_maxiter, &
          etot_maxiter_ => etot_maxiter, &
          forc_maxiter_ => forc_maxiter


      USE control_flags, ONLY: &
        force_pairing_ => force_pairing
     
     !   Other module

     USE wave_base, ONLY: frice_ => frice
     USE ions_base, ONLY: fricp_ => fricp
     USE cell_base, ONLY: frich_ => frich

     USE input_parameters, ONLY: &
        electron_dynamics, electron_damping, diis_rot, electron_temperature, &
        ion_dynamics, ekin_conv_thr, etot_conv_thr, forc_conv_thr, ion_maxstep, &
        electron_maxstep, ion_damping, ion_temperature, ion_velocities, tranp, &
        amprp, ion_nstepe, cell_nstepe, cell_dynamics, cell_damping,  &
        cell_parameters, cell_velocities, cell_temperature, force_pairing, &
        tapos, tavel, ecutwfc, emass_cutoff, taspc, trd_ht, ibrav, ortho_eps, &
        ortho_max, ntyp, tolp, tchi2_inp, calculation, disk_io

     USE input_parameters, ONLY: ndr, ndw, iprint, isave, tstress, k_points, &
        tprnfor, verbosity, tprnrho, tdipole_card, toptical_card, &
        tnewnfi_card, newnfi_card, ampre, nstep, restart_mode, ion_positions, &
        startingwfc, printwfc, orthogonalization, electron_velocities


     !
     IMPLICIT NONE
     !

     IF( .NOT. has_been_read ) &
       CALL errore( ' iosys ', ' input file has not been read yet! ', 1 )

     ndr_        = ndr
     ndw_        = ndw
     iprint_     = iprint
     isave_      = isave
     tstress_    = tstress
     tpre_       = tstress
     gamma_only_ = ( TRIM( k_points ) == 'gamma' )
     tprnfor_    = tprnfor
     printwfc_   = printwfc
     tchi2_      = tchi2_inp
     ekin_conv_thr_ = ekin_conv_thr
     etot_conv_thr_ = etot_conv_thr
     forc_conv_thr_ = forc_conv_thr
     ekin_maxiter_  = electron_maxstep
!     etot_maxiter_ = etot_maxiter, &
!     forc_maxiter_ = forc_maxiter

     !
     !  set the level of output, the code verbosity 
     !

     iprsta_ = 1
     prn_    = .FALSE.
     timing_ = .FALSE.
          ! The code write to files fort.8 fort.41 fort.42 fort.43
          ! a detailed report of subroutines timing
     rhoout_ = .false.
          ! save charge density to file  CHARGEDENSITY if nspin = 1, and
          ! CHARGEDENSITY.UP CHARGEDENSITY.DOWN if nspin = 2
     memchk_ = .FALSE.
          ! The code performs a memory check, write on standard
          ! output the allocated memory at each step.
          ! Architecture Dependent
     tprnsfac_   = .FALSE.
          ! Print on file STRUCTURE_FACTOR the structure factor
          ! gvectors and charge density, in reciprocal space.

     trhor_  = ( TRIM( calculation ) == 'nscf' )
     trhow_  = ( TRIM( disk_io ) == 'high' )
     tvlocw_ = .false. ! temporaneo

     SELECT CASE ( TRIM(verbosity) )
       CASE ('minimal')
         prn_ = .FALSE.
       CASE ('low', 'default')
         prn_ = .FALSE.
         timing_ = .TRUE.
       CASE ('medium')
         prn_ = .FALSE.
         timing_ = .TRUE.
         rhoout_ = .TRUE.
         tprnsfac_ = .TRUE.
       CASE ('high')
         iprsta_ = 3
         prn_ = .TRUE.
         memchk_ = .TRUE.
         timing_ = .TRUE.
         rhoout_ = .TRUE.
         tprnsfac_ = .TRUE.
       CASE DEFAULT
         CALL errore(' control_flags ',' unknown verbosity '//TRIM(verbosity), 1 )
     END SELECT

     ! ...   If explicitly requested force the charge density to be printed
     IF( tprnrho ) rhoout_ = .TRUE.

     tdipole_  = tdipole_card
     toptical_ = toptical_card
     newnfi_   = newnfi_card
     tnewnfi_  = tnewnfi_card

     !
     !   set the restart flags
     !

     trane_ = .FALSE.
     ampre_ = ampre
     SELECT CASE ( TRIM( restart_mode ) )
       CASE ('from_scratch')
         nbeg_ = -2
         if ( ion_positions == 'from_input' ) nbeg_ = -1
         nomore_ = nstep
         trane_  = ( startingwfc == 'random' )
         if ( ampre_ == 0.d0 ) ampre_ = 0.02
       CASE ('reset_counters')
         nbeg_ = 0
         nomore_ = nstep
       CASE ('restart')
         nbeg_ = 1
         nomore_ = nstep
         if ( ion_positions == 'from_input' ) then
           taurdr_ = .TRUE.
           nbeg_ = -1
         end if
       CASE DEFAULT
         CALL errore(' iosys ',' unknown restart_mode '//trim(restart_mode), 1 )
     END SELECT

     ! ... Starting/Restarting Atomic positions
     taurdr_ = .FALSE.
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
              trim(orthogonalization), 1 )
      END SELECT

      ortho_max_ = ortho_max
      ortho_eps_ = ortho_eps

      ! ... Electrons initial velocity

      SELECT CASE ( TRIM(electron_velocities) )
        CASE ('default')
          tzeroe_ = .FALSE.
        CASE ('zero')
          tzeroe_ = .TRUE.
          IF( program_name == 'CP90' ) &
            print '("Warning: electron_velocities keyword has no effect")'
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown electron_velocities '//TRIM(electron_velocities), 1 )
      END SELECT


      ! ... Electron dynamics

      tdamp_          = .FALSE.
      tconjgrad_      = .FALSE.
      tsteepdesc_     = .FALSE.
      t_diis_        = .FALSE.
      t_diis_simple_ = .FALSE.
      t_diis_rot_    = .FALSE.
      frice_ = 0.d0
      SELECT CASE ( TRIM(electron_dynamics) )
        CASE ('sd', 'default')
          tsde_ = .TRUE.
        CASE ('verlet')
          tsde_ = .FALSE.
        CASE ('cg')
          tsde_      = .FALSE.
          tconjgrad_ = .TRUE.
          IF( program_name == 'CP90' ) &
            CALL errore( "iosys ", " electron_dynamics keyword not yet implemented ", 1 )
        CASE ('damp')
          tsde_   = .FALSE.
          tdamp_  = .TRUE.
          frice_ = electron_damping
        CASE ('diis')
          IF( program_name == 'CP90' ) &
            CALL errore( "iosys ", " electron_dynamics keyword not yet implemented ", 1 )
          tsde_   = .FALSE.
          t_diis_ = .TRUE.
          IF( diis_rot ) THEN
            t_diis_rot_    = .TRUE.
          ELSE
            t_diis_simple_ = .TRUE.
          END IF
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

      tdampions_        = .FALSE.
      tconvthrs%active = .FALSE.
      tconvthrs%nstep  = 1
      tconvthrs%ekin   = 0.0d0
      tconvthrs%derho  = 0.0d0
      tconvthrs%force  = 0.0d0
      tconjgrad_ion%active = .FALSE.
      tconjgrad_ion%nstepix = 1
      tconjgrad_ion%nstepex = 1
      tconjgrad_ion%ionthr = 1.0d+10
      tconjgrad_ion%elethr = 1.0d+10
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
          tsdp_ = .FALSE.
          tfor_ = .TRUE.
          tconjgrad_ion%active  = .TRUE.
          tconjgrad_ion%nstepix = ion_maxstep    ! maximum number of iteration
          tconjgrad_ion%nstepex = electron_maxstep  ! maximum number of iteration for the electronic minimization
          tconjgrad_ion%ionthr  = etot_conv_thr ! energy threshold for convergence
          tconjgrad_ion%elethr  = ekin_conv_thr ! energy threshold for convergence in the electrons minimization
          tconvthrs%ekin   = ekin_conv_thr
          tconvthrs%derho  = etot_conv_thr
          tconvthrs%force  = forc_conv_thr
          tconvthrs%active = .TRUE.
          tconvthrs%nstep  = 1
          IF( program_name == 'CP90' ) &
            CALL errore( "iosys ", " ion_dynamics = '//TRIM(ion_dynamics)//' not yet implemented ", 1 )
        CASE ('damp')
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

      ! ... Ionic Temperature

      tcp_      = .FALSE.
      tnosep_   = .FALSE.
      tolp_     = tolp
      SELECT CASE ( TRIM(ion_temperature) )
        !         temperature control of ions via Nose' thermostat
        CASE ('nose')
          tnosep_ = .TRUE.
          tcp_ = .false.
        CASE ('not_controlled', 'default')
          tnosep_ = .FALSE.
          tcp_ = .false.
        CASE ('rescaling' )
          tnosep_ = .FALSE.
          tcp_ = .true.
        CASE DEFAULT
          CALL errore(' control_flags ',' unknown ion_temperature '//TRIM(ion_temperature), 1 )
      END SELECT


      ! ... Starting/Restarting ionic velocities

      tcap_         = .FALSE.
      SELECT CASE ( TRIM(ion_velocities) )
        CASE ('default')
          tzerop_ = .FALSE.
          tv0rd_ = .FALSE.
          tcap_ = .false.
        CASE ('zero')
          tzerop_ = .TRUE.
          tv0rd_ = .FALSE.
        CASE ('from_input')
          tzerop_ = .TRUE.
          tv0rd_  = .TRUE.
        CASE ('random')
          tcap_ = .true.
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
            WRITE(stdout) " ion_nstepe or cell_nstepe have no effects "
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
          tpre_ = .FALSE.
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
        tconvthrs%force  = 10d+10
        tconvthrs%active = .TRUE.
        tconvthrs%nstep  = 1
      END IF

      ! force pairing

      force_pairing_ = force_pairing
      IF( program_name == 'CP90' .AND. force_pairing_) &
            WRITE(stdout) " force_pairing have no effects "


      ! . Set internal flags according to the input .......................!

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

      IF( TRIM(electron_dynamics) /= 'diis' .AND. TRIM(orthogonalization) /= 'ortho' ) THEN
        IF( emass_cutoff < ecutwfc ) &
          CALL errore(' IOSYS ', ' FOURIER ACCELERATION WITHOUT ORTHO',0)
      END IF



      RETURN
   END SUBROUTINE


   !  ----------------------------------------------
   !
   !
   !  
   !  modules setup ( copy-in variables )
   !
   !
   !
   !  ----------------------------------------------

   SUBROUTINE modules_setup( )
     !
     USE kinds,            ONLY: dbl
     USE control_flags,    ONLY: lneb, program_name
     USE constants,        ONLY: UMA_AU, pi
     USE io_global,        ONLY: stdout
   
     USE input_parameters, ONLY: title, max_seconds, ibrav , celldm , trd_ht, &
           cell_symmetry, rd_ht, a, b, c, cosab, cosac, cosbc, ntyp , nat ,   &
           na_inp , sp_pos , rd_pos , rd_vel, atom_mass, atom_label, if_pos,  &
           atomic_positions, id_loc, sic, sic_epsilon, sic_rloc, ecutwfc,     &
           ecutrho, ecfixed, qcutz, q2sigma, tk_inp, nkstot, xk, dt, wmass,   &
           ion_radius, emass, emass_cutoff, temph, fnoseh, nr1b, nr2b, nr3b,  &
           tempw, fnosep, nr1, nr2, nr3, nr1s, nr2s, nr3s, ekincw, fnosee,    &
           tturbo_inp, nturbo_inp, outdir, prefix, xc_type, woptical,         &
           noptical, boptical, k_points, nkstot, nk1, nk2, nk3, k1, k2, k3,   &
           xk, wk

     USE input_parameters, ONLY: diis_achmix, diis_ethr, diis_wthr, diis_delt, &
           diis_nreset, diis_temp, diis_nrot, diis_maxstep, diis_fthr,         &
           diis_size, diis_hcut, diis_rothr, diis_chguess, diis_g0chmix,       &
           diis_nchmix, diis_g1chmix, empty_states_maxstep, empty_states_delt, &
           empty_states_emass, empty_states_ethr, empty_states_nbnd,           &
           tprnks_empty, vhrmax_inp, vhnr_inp, vhiunit_inp, vhrmin_inp,        &
           tvhmean_inp, vhasse_inp, ion_damping, anner_inp, constr_dist_inp,   &
           constr_inp, nconstr_inp, constr_tol_inp, constr_type_inp, iesr_inp, &
           etot_conv_thr, ekin_conv_thr, nspin, f_inp, nelup, neldw, nbnd,     &
           nelec, anne_inp, tprnks, ks_path, tneighbo, neighbo_radius, press,  &
           cell_damping, cell_dofree, tf_inp, tpstab_inp, pstab_size_inp,      &
           greash, grease, greasp, epol, efield, tcg, maxiter, etresh, passop

     !
     USE check_stop,       ONLY: check_stop_init
     !
     USE printout_base,    ONLY: title_ => title
     !
     USE cell_base,        ONLY: cell_base_init, a1, a2, a3
     USE cell_nose,        ONLY: cell_nose_init
     USE ions_base,        ONLY: ions_base_init, greasp_ => greasp
     USE ions_nose,        ONLY: ions_nose_init
     USE wave_base,        ONLY: grease_ => grease
     USE electrons_nose,   ONLY: electrons_nose_init
     USE printout_base,    ONLY: printout_base_init
     USE time_step,        ONLY: set_time_step
     USE turbo,            ONLY: turbo_init
     USE efield_module,    ONLY: efield_init
     USE cg_module,        ONLY: cg_init
     !
     USE reciprocal_space_mesh,    ONLY: recvecs_units
     !
     USE cp_electronic_mass,       only: &
           emass_ => emass, &
           emaec_ => emass_cutoff

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
     !
     USE exchange_correlation,     ONLY: exch_corr_init
     USE brillouin,                ONLY: kpoint_setup
     USE optical_properties,       ONLY: optical_setup
     USE pseudopotential,          ONLY: pseudopotential_setup
     USE cell_base,                ONLY: a1, a2, a3, cell_alat
     USE guess,                    ONLY: guess_setup
     USE ions_module,              ONLY: ions_setup
     USE empty_states,             ONLY: empty_init
     USE diis,                     ONLY: diis_setup
     USE charge_mix,               ONLY: charge_mix_setup
     USE potentials,               ONLY: potential_init
     USE kohn_sham_states,         ONLY: ks_states_init
     USE electrons_module,         ONLY: electrons_setup


     !
     IMPLICIT NONE

     REAL(dbl) :: alat_ , massa_totale
     REAL(dbl)  :: delt_emp_inp, emass_emp_inp, ethr_emp_inp
     ! ...   DIIS
     REAL(dbl) :: tol_diis_inp, delt_diis_inp, tolene_inp
     LOGICAL :: o_diis_inp, oqnr_diis_inp
     !
     !   Subroutine Body
     !
     IF( .NOT. has_been_read ) &
       CALL errore( ' modules_setup ', ' input file has not been read yet! ', 1 )
     !
     title_ = title
     !
     IF( .NOT. lneb ) THEN
        CALL check_stop_init( max_seconds )
     END IF

     ! ...  Set cell base module

     IF( .not. lneb ) THEN
        massa_totale = SUM( atom_mass(1:ntyp)*na_inp(1:ntyp) )
        CALL cell_base_init( ibrav , celldm , trd_ht, cell_symmetry, rd_ht,  &
               a, b, c, cosab, cosac, cosbc , wmass , massa_totale , press , &
               cell_damping, greash , cell_dofree, alat_ )
     END IF

     alat_ = cell_alat()


     ! ...  Set ions base module

     if( .not. lneb ) then
        IF( program_name == 'CP90' ) THEN
          CALL ions_base_init( ntyp , nat , na_inp , sp_pos , rd_pos , rd_vel, atom_mass, &
             atom_label, if_pos, atomic_positions , alat_ , a1, a2, a3, ion_radius  )
        ELSE
          CALL ions_base_init( ntyp , nat , na_inp , sp_pos , rd_pos , rd_vel, atom_mass, &
               atom_label, if_pos, atomic_positions, alat_ , a1 , a2 , a3 , ion_radius,   &
               id_loc, sic, sic_epsilon, sic_rloc  )
        END IF
     end if

     ! ...   Set units for Reciprocal vectors ( 2PI/alat by convention )

     CALL recvecs_units( alat_ )

     ! ...   Set Values for the cutoff

     CALL ecutoffs_setup( ecutwfc, ecutrho, ecfixed, qcutz, q2sigma )

     CALL gcutoffs_setup( alat_ , tk_inp, nkstot, xk )

     emass_ = emass
     emaec_ = emass_cutoff

     ! ...   Set internal time step variables ( delt, twodelt, dt2 ... )

     CALL set_time_step( dt )

     ! ... 

     grease_ = grease
     greasp_ = greasp

     !   set thermostat parameter for cell, ions and electrons

     CALL cell_nose_init( temph, fnoseh )

     CALL ions_nose_init( tempw, fnosep, nat )
  
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

     IF( program_name == 'FPMD' ) CALL exch_corr_init( xc_type )

     CALL turbo_init( tturbo_inp, nturbo_inp )

     CALL printout_base_init( outdir, prefix )

     IF( noptical > 0 ) then
        CALL optical_setup( woptical, noptical, boptical )
     END IF

     CALL kpoint_setup( k_points, nkstot, nk1, nk2, nk3, k1, k2, k3, xk, wk )

     CALL efield_init( epol, efield )

     CALL cg_init( tcg , maxiter , etresh , passop )

     !
     !  empty states
     !
     delt_emp_inp  = dt
     ethr_emp_inp  = ekin_conv_thr
     IF( empty_states_delt > 0.d0 )  delt_emp_inp  = empty_states_delt
     IF( empty_states_ethr > 0.d0 )  ethr_emp_inp  = empty_states_ethr
     CALL empty_init( empty_states_maxstep, delt_emp_inp, ethr_emp_inp )

     !
     CALL potential_init( tvhmean_inp,vhnr_inp, vhiunit_inp, &
               vhrmin_inp, vhrmax_inp, vhasse_inp, iesr_inp)

     CALL ks_states_init( nspin, tprnks, tprnks_empty )

    
     IF( program_name == 'FPMD' ) THEN

        CALL pseudopotential_setup( ntyp, tpstab_inp, pstab_size_inp, ion_radius )

        CALL electrons_setup( tf_inp, nbnd, nint(nelec), nint(nelup), &
          nint(neldw), nspin, empty_states_nbnd, emass, emass_cutoff, &
          f_inp, nkstot )

        CALL ions_setup( anne_inp, anner_inp,   &
             nconstr_inp, constr_tol_inp, constr_type_inp, constr_dist_inp, &
             constr_inp, ion_damping, tneighbo, neighbo_radius )
        !

        o_diis_inp        = .TRUE.
        oqnr_diis_inp     = .TRUE.
        tolene_inp        = etot_conv_thr
        tol_diis_inp      = ekin_conv_thr
        delt_diis_inp     = dt
        IF( diis_ethr > 0.0d0 ) tolene_inp    = diis_ethr
        IF( diis_wthr > 0.0d0 ) tol_diis_inp  = diis_wthr
        IF( diis_delt > 0.0d0 ) delt_diis_inp = diis_delt
        CALL diis_setup( diis_fthr, oqnr_diis_inp, o_diis_inp, &
          diis_size, diis_hcut, tol_diis_inp, diis_maxstep, diis_nreset, delt_diis_inp, &
          diis_temp, diis_nrot(1), diis_nrot(2), diis_nrot(3), &
          diis_rothr(1), diis_rothr(2), diis_rothr(3), tolene_inp)
        CALL guess_setup( diis_chguess )
        CALL charge_mix_setup(diis_achmix, diis_g0chmix, diis_nchmix, diis_g1chmix)


     END IF

     RETURN
   END SUBROUTINE modules_setup

   !
   !

!

    subroutine modules_setup_cp( )

!     this subroutine copies variables from input module to other modules
!     -------------------------------------------------------------------

      use input_parameters, only: &
           nr1, nr2, nr3, nr2s, nr3s, nr1s, &
           tempw, atomic_positions, nelec, &
           if_pos, rd_ht, trd_ht, a, b, c, cosab, cosac, cosbc, cell_symmetry, nelup, &
           neldw, occupations, f_inp, pos, pseudo_dir, &
           sp_pos, atom_mass, atom_pfile, &
           startingwfc, ion_dynamics, ion_damping, &
           cell_velocities, electron_dynamics, ion_velocities, &
           celldm, nbnd, nspin, calculation, ntyp, ibrav, restart_mode, ion_positions, &
           ecutwfc, ecutrho, ortho_eps, ortho_max, qcutz, q2sigma, &
           ecfixed, fnosep, nat, ion_temperature, &
           cell_temperature, cell_dynamics, cell_damping, electron_temperature, &
           dt, emass, emass_cutoff, ion_radius, &
           ekin_conv_thr, etot_conv_thr, na_inp, rd_pos, atom_label, rd_vel, &
           smd_polm, smd_kwnp, smd_linr, smd_stcd, smd_stcd1, smd_stcd2, smd_stcd3, smd_codf, &
           smd_forf, smd_smwf, smd_lmfreq, smd_tol, smd_maxlm, smd_smcp, smd_smopt, smd_smlm, &
           num_of_images, smd_ene_ini, smd_ene_fin, &
           n_inner, fermi_energy, rotmass, occmass, rotation_damping,                       &
           occupation_damping, occupation_dynamics, rotation_dynamics,                      &
           degauss, smearing, tk_inp, nkstot, xk


      use constants, only: pi, scmass, factem, eps8, uma_au, terahertz 

      use parameters, only: natx

      use io_global, only: ionode, stdout

      use control_flags, only:  &
            tconvthrs, lneb, lsmd, tzerop, tzeroe, tzeroc, nbeg,  &
            ndr_ => ndr, &
            ndw_ => ndw, &
            nomore_ => nomore, &
            iprint_ => iprint, &
            iprsta_ => iprsta, &
            ortho_eps_ => ortho_eps, &
            ortho_max_ => ortho_max, &
            trane_ => trane, &
            ampre_ => ampre, &
            tfor_ => tfor, &
            tsdp_ => tsdp, &
            tcp_ => tcp, &
            tcap_ => tcap, &
            tnosep_ => tnosep, &
            tnosee_ => tnosee, &
            tnoseh_ => tnoseh, &
            tpre_ => tpre, &
            thdyn_ => thdyn, &
            tsde_ => tsde

      use mp, only: mp_bcast
      !
      USE ions_base, ONLY: tau_srt, ind_srt, &
           rcmax_ => rcmax, &
           fricp_ => fricp
      !
      USE ions_positions, ONLY: &
           tau0_ => tau0
      ! 
      USE cell_base, ONLY: cell_alat, a1, a2, a3

      use gvecw, only: agg => ecutz, sgg => ecsig, e0gg => ecfix

      USE time_step, ONLY: set_time_step, &
           delt_ => delt

      USE cp_electronic_mass, only: &
           emass_ => emass, &
           emaec_ => emass_cutoff

      USE wave_base, ONLY: frice_ => frice

      USE ions_nose, ONLY: &
           qnp_ => qnp, &
           tempw_ => tempw
      USE electrons_base, ONLY: &
           nupdwn_ => nupdwn, &
           iupdwn_ => iupdwn, &
           nel_ => nel, &
           n_ => nbnd, &
           nx_ => nbndx, &
           f_ => f, &
           ispin_ => fspin, &
           nspin_ => nspin

      USE path_variables, ONLY: &
           sm_p_ => smd_p, &
           smcp_ => smd_cp, &
           smlm_ => smd_lm, &
           smopt_ => smd_opt, &
           linr_ => smd_linr, &
           polm_ => smd_polm, &
           kwnp_ => smd_kwnp, &
           codfreq_ => smd_codfreq, &
           forfreq_ => smd_forfreq, &
           smwfreq_ => smd_wfreq, &
           tol_ => smd_tol, &
           lmfreq_ => smd_lmfreq, &
           maxlm_ => smd_maxlm, &
           ene_ini_ => smd_ene_ini, &
           ene_fin_ => smd_ene_fin


      USE ensemble_dft, ONLY: &
           tens_ => tens, &
           tgrand_ => tgrand, &
           ninner_ => ninner, &
           ismear_ => ismear, &
           etemp_ => etemp, &
           ef_  => ef, &
           tdynz_ => tdynz, &
           tdynf_ => tdynf, &
           zmass_ => zmass, &
           fmass_ => fmass, &
           fricz_ => fricz, &
           fricf_ => fricf


      !
      implicit none
      !
      !
      ! local variables
      !

      real(kind=8) :: ocp, fsum
      integer :: i, ia, is, iss, in, isa
      real(kind=8) :: alat_

      !
      ! Subroutine body
      !

      alat_ = cell_alat()

      ! 
      !     translate from input to internals of SMCP, 
      !
      ! ... SM_P  
      !
      sm_p_ = num_of_images -1
      !
      ! ... what to do
      !
      smcp_   = smd_smcp
      smopt_  = smd_smopt
      smlm_   = smd_smlm
      !      
      ! ... initial path info
      !
      linr_ = smd_linr
      polm_ = smd_polm
      kwnp_ = smd_kwnp
      !
      ! ...  Frequencey of wiriting
      !
      codfreq_ = smd_codf
      forfreq_ = smd_forf
      smwfreq_ = smd_smwf
      !
      ! ... Lagrange multiplier info.
      !
      lmfreq_ = smd_lmfreq
      tol_    = smd_tol
      maxlm_  = smd_maxlm
      !
      ! ... if smlm
      !
      IF( smd_smlm .AND. ( smd_ene_ini >= 0.d0 .OR. smd_ene_fin >= 0.d0 ) ) THEN
         CALL errore(' start : ',' Check : ene_ini & ene_fin ', 1 )
      END IF
      !
      ene_ini_ = smd_ene_ini
      ene_fin_ = smd_ene_fin
      !
      IF( lsmd .AND. cell_dynamics /= 'none' ) THEN
        CALL errore(' smiosys ',' cell_dynamics not implemented : '//trim(cell_dynamics), 1 )
      END IF
      !
      IF( lsmd ) THEN

         !
         ! How to obtain the initial trial path.
         !

         IF(smd_smopt) THEN
   
          CALL init_path(sm_p_,kwnp_,smd_stcd,ntyp,nat,alat_,nbeg,1)

         ELSEIF(smd_linr) THEN

          CALL init_path(sm_p_,kwnp_,smd_stcd,ntyp,nat,alat_,nbeg,2)

         ELSEIF(smd_polm .AND. (smd_kwnp < num_of_images) ) THEN

          CALL init_path(sm_p_,kwnp_,smd_stcd,ntyp,nat,alat_,nbeg,3)

         ELSEIF(smd_kwnp == num_of_images ) THEN

          CALL init_path(sm_p_,kwnp_,smd_stcd,ntyp,nat,alat_,nbeg,4)

         ENDIF

      ELSE

         tau0_ = 0.0d0
         tau0_ ( 1:3 , 1:nat ) = tau_srt ( 1:3 , 1:nat )

      END IF


      !
      !
      !  set occupancies
      !

      ! ...   Set Values for bands and spin

      n_     = nbnd * nspin
      nspin_ = nspin
      
      IF( nelec < 1 ) THEN
         CALL errore(' iosys ',' nelec less than 1 ', int(nelec) )
      END IF
      IF( nint(nelec) - nelec > eps8 ) THEN
         CALL errore(' iosys ',' nelec must be integer', int(nelec) )
      END IF

      if( mod( n_ , 2 ) .ne. 0 ) then
         nx_ = n_ + 1
      else
         nx_= n_
      end if

      ALLOCATE( f_ ( nx_ ) )
      ALLOCATE( ispin_ ( nx_ ) )
      f_     = 0.0d0
      ispin_ = 0

      iupdwn_ ( 1 ) = 1
      nel_ = 0


      SELECT CASE ( TRIM(occupations) ) 
      CASE ('bogus')
         !
         ! empty-states calculation: occupancies have a (bogus) finite value
         !
         ! bogus to ensure \sum_i f_i = Nelec  (nelec is integer)
         !
         f_ ( : ) = nelec / n_         
         nel_ (1) = nint(nelec)
         nupdwn_ (1) = n_
         if ( nspin_ == 2 ) then
            !
            ! bogus to ensure Nelec = Nup + Ndw
            !
            nel_ (1) = ( nint(nelec) + 1 ) / 2
            nel_ (2) =   nint(nelec)       / 2
            nupdwn_ (1)=nbnd
            nupdwn_ (2)=nbnd
            iupdwn_ (2)=nbnd+1
         end if
      CASE ('from_input')
         !
         ! occupancies have been read from input
         !
         f_ ( 1:nbnd ) = f_inp( 1:nbnd, 1 )
         if( nspin_ == 2 ) f_ ( nbnd+1 : 2*nbnd ) = f_inp( 1:nbnd, 2 ) 
         if( nelec == 0.d0 ) nelec = SUM ( f_ ( 1:n_ ) )
         if( nspin_ == 2 .and. nelup == 0) nelup = SUM ( f_ ( 1:nbnd ) )
         if( nspin_ == 2 .and. neldw == 0) neldw = SUM ( f_ ( nbnd+1 : 2*nbnd ) )

         if( nspin_ == 1 ) then 
           nel_ (1) = nint(nelec)
           nupdwn_ (1) = n_
         else
           IF ( ABS (nelup + neldw - nelec) > eps8 ) THEN
              CALL errore(' iosys ',' wrong # of up and down spin', 1 )
           END IF
           nel_ (1) = nint(nelup)
           nel_ (2) = nint(neldw)
           nupdwn_ (1)=nbnd
           nupdwn_ (2)=nbnd
           iupdwn_ (2)=nbnd+1
         end if

      CASE ('fixed')

         if( nspin_ == 1 ) then
            nel_ (1) = nint(nelec)
            nupdwn_ (1) = n_
         else
            IF ( nelup + neldw /= nelec  ) THEN
               CALL errore(' iosys ',' wrong # of up and down spin', 1 )
            END IF
            nel_ (1) = nint(nelup)
            nel_ (2) = nint(neldw)
            nupdwn_ (1)=nbnd
            nupdwn_ (2)=nbnd
            iupdwn_ (2)=nbnd+1
         end if

         ! ocp = 2 for spinless systems, ocp = 1 for spin-polarized systems
         ocp = 2.d0 / nspin_
         ! default filling: attribute ocp electrons to each states
         !                  until the good number of electrons is reached
         do iss = 1, nspin_
            fsum = 0.0d0
            do in = iupdwn_ ( iss ), iupdwn_ ( iss ) - 1 + nupdwn_ ( iss )
               if ( fsum + ocp < nel_ ( iss ) + 0.0001 ) then
                  f_ (in) = ocp
               else
                  f_ (in) = max( nel_ ( iss ) - fsum, 0.d0 )
               end if
                fsum=fsum + f_(in)
            end do
         end do

      CASE ('grand-canonical','g-c','gc')
          tens_    =.true.
          tgrand_  =.true.
          CALL errore(' iosys ','grand-canonical not yet implemented ', 1 )

      CASE ('ensemble','ensemble-dft','edft')
          tens_    =.true.
          ninner_  = n_inner
          etemp_   = degauss
          ef_      = fermi_energy
          fricz_   = rotation_damping
          fricf_   = occupation_damping
          zmass_   = rotmass
          fmass_   = occmass

          SELECT CASE (rotation_dynamics)
            CASE ( 'line-minimization','l-m','lm' )
              tdynz_ = .FALSE.
              fricz_ = 0.0d0
              zmass_ = 0.0d0
            CASE DEFAULT
              CALL errore(' iosys ',' rotation_dynamics not implemented ', 1 )
          END SELECT

          SELECT CASE (occupation_dynamics)
            CASE ( 'line-minimization','l-m','lm' )
              tdynf_ = .FALSE.
              fricf_ = 0.0d0
              fmass_ = 0.0d0
            CASE DEFAULT
              CALL errore(' iosys ',' occupation_dynamics not implemented ', 1 )
          END SELECT

          if ( nspin_ == 1 ) then
            n_       = nbnd
            f_ ( : ) = nelec / n_
            nel_ (1) = nint(nelec)
            nupdwn_ (1) = n_
          else
            n_       = 2*nbnd
            if (nelup.ne.0) then
              if ((nelup+neldw).ne.nelec) then
                 CALL errore(' iosys ',' nelup+neldw .ne. nelec', 1 )
              end if
              nel_ (1) = nelup
              nel_ (2) = neldw
            else
              nel_ (1) = ( nint(nelec) + 1 ) / 2
              nel_ (2) =   nint(nelec)       / 2
            end if
            nupdwn_ (1) = nbnd
            nupdwn_ (2) = nbnd
            iupdwn_ (2) = nbnd+1
            do iss = 1, nspin_
             do i = iupdwn_ ( iss ), iupdwn_ ( iss ) - 1 + nupdwn_ ( iss )
                f_ (i) =  nel_ (iss) / real (nupdwn_ (iss))
             end do
            end do
          end if

          SELECT CASE (smearing)
            CASE ( 'gaussian','g' )
              ismear_ = 1
            CASE ( 'fermi-dirac','f-d', 'fd' )
              ismear_ = 2
            CASE ( 'hermite-delta','h-d','hd' )
              ismear_ = 3
            CASE ( 'gaussian-splines','g-s','gs' )
              ismear_ = 4
            CASE ( 'cold-smearing','c-s','cs','cs1' )
              ismear_ = 5
            CASE ( 'marzari-vanderbilt','m-v','mv','cs2' )
              ismear_ = 6
            CASE ( '0')
              ismear_ = 0
            CASE ( '-1')
              ismear_ = -1

            CASE DEFAULT
              CALL errore(' iosys ',' smearing not implemented', 1 )
          END SELECT

      CASE DEFAULT
         CALL errore(' iosys ',' occupation method not implemented', 1 )
      END SELECT

      do iss = 1, nspin_
         do in = iupdwn_(iss), iupdwn_(iss) - 1 + nupdwn_(iss)
            ispin_(in) = iss
         end do
      end do

      !

      RETURN

  END SUBROUTINE

!
!     --------------------------------------------------------
!     print out heading
!


  SUBROUTINE input_info_cp()

      use constants, only: pi, scmass, factem, eps8, uma_au, terahertz, gpa_au

      use parameters, only: natx

      use io_global, only: ionode, stdout

      use control_flags, only:  &
            tconvthrs, lneb, lsmd, tzerop, tzeroe, tzeroc, nbeg,  &
            ndr_ => ndr, &
            ndw_ => ndw, &
            nomore_ => nomore, &
            iprint_ => iprint, &
            iprsta_ => iprsta, &
            tortho_ => tortho, &
            ortho_eps_ => ortho_eps, &
            ortho_max_ => ortho_max, &
            trane_ => trane, &
            ampre_ => ampre, &
            tranp_ => tranp, &
            amprp_ => amprp, &
            tfor_ => tfor, &
            tsdp_ => tsdp, &
            tcp_ => tcp, &
            tcap_ => tcap, &
            tolp_ => tolp, &
            trhor_ => trhor, &
            trhow_ => trhow, &
            tvlocw_ => tvlocw, &
            tnosep_ => tnosep, &
            tnosee_ => tnosee, &
            tnoseh_ => tnoseh, &
            tpre_ => tpre, &
            thdyn_ => thdyn, &
            tsde_ => tsde

      use mp, only: mp_bcast
      !
      USE ions_base, ONLY: tau_srt, ind_srt, &
           rcmax_ => rcmax, &
           fricp_ => fricp, &
           greasp_ => greasp, &
           nsp
      !
      USE ions_positions, ONLY: &
           tau0_ => tau0
      ! 
      USE cell_base, ONLY: cell_alat, a1, a2, a3, &
           press_ => press, &
           frich_ => frich, &
           greash_ => greash, &
           thdiag_ => thdiag, &
           iforceh_ => iforceh
      USE cell_nose, ONLY: &
           qnh_ => qnh, &
           temph_ => temph

      use gvecw, only: agg => ecutz, sgg => ecsig, e0gg => ecfix

      USE time_step, ONLY: set_time_step, &
           delt_ => delt
      USE cp_electronic_mass, only: &
           emass_ => emass, &
           emaec_ => emass_cutoff
      USE wave_base, ONLY: &
           frice_ => frice, &
           grease_ => grease
      USE ions_nose, ONLY: &
           qnp_ => qnp, &
           tempw_ => tempw
      USE electrons_base, ONLY: &
           nupdwn_ => nupdwn, &
           iupdwn_ => iupdwn, &
           nel_ => nel, &
           n_ => nbnd, &
           nx_ => nbndx, &
           f_ => f, &
           ispin_ => fspin, &
           nspin_ => nspin
      USE electrons_nose, ONLY: &
           qne_ => qne, &
           ekincw_ => ekincw
      USE path_variables, ONLY: &
           sm_p_ => smd_p, &
           smcp_ => smd_cp, &
           smlm_ => smd_lm, &
           smopt_ => smd_opt, &
           linr_ => smd_linr, &
           polm_ => smd_polm, &
           kwnp_ => smd_kwnp, &
           codfreq_ => smd_codfreq, &
           forfreq_ => smd_forfreq, &
           smwfreq_ => smd_wfreq, &
           tol_ => smd_tol, &
           lmfreq_ => smd_lmfreq, &
           maxlm_ => smd_maxlm, &
           ene_ini_ => smd_ene_ini, &
           ene_fin_ => smd_ene_fin

      USE grid_dimensions, ONLY: &
           nnrx, &  !  variable is used to workaround internal compiler error (IBM xlf)
           nr1_ => nr1, &
           nr2_ => nr2, &
           nr3_ => nr3
      USE smooth_grid_dimensions, ONLY: &
           nnrsx, &  !  variable is used to workaround internal compiler error (IBM xlf)
           nr1s_ => nr1s, &
           nr2s_ => nr2s, &
           nr3s_ => nr3s
      !
      USE io_files, ONLY: &
           ntypx, &     !  variable is used to workaround internal compiler error (IBM xlf)
           pseudo_dir_ => pseudo_dir , &
           psfile_     => psfile

      USE ensemble_dft, ONLY: &
           tens_ => tens, &
           tgrand_ => tgrand, &
           ninner_ => ninner, &
           ismear_ => ismear, &
           etemp_ => etemp, &
           ef_  => ef, &
           tdynz_ => tdynz, &
           tdynf_ => tdynf, &
           zmass_ => zmass, &
           fmass_ => fmass, &
           fricz_ => fricz, &
           fricf_ => fricf
      USE cg_module, ONLY: &
           tcg_ => tcg, &
           maxiter_ => maxiter, &
           etresh_ => etresh, &
           passop_ => passop

      !

      IMPLICIT NONE
      !

      INTEGER :: is

      WRITE( stdout,500) nbeg , nomore_ , iprint_ , ndr_ , ndw_
      WRITE( stdout,505) delt_
      WRITE( stdout,510) emass_ , emaec_
!
      if( tortho_ ) then
         WRITE( stdout,511) ortho_eps_ , ortho_max_
      else
         WRITE( stdout,512)
      endif
!
      if( tsde_ ) then
         WRITE( stdout,513)
      else
         if ( tnosee_ ) frice_ = 0.
         WRITE( stdout,509)
         WRITE( stdout,514) frice_ , grease_
      endif
!
      if ( trhor_ ) then
         WRITE( stdout,720)
      endif
!
      if( .not. trhor_ .and. trhow_ )then
         WRITE( stdout,721)
      endif
!
      if( tvlocw_ )then
         WRITE( stdout,722)
      endif
!
      if( trane_ ) then
         WRITE( stdout,515) ampre_
      endif
      WRITE( stdout,516)
      do is =1, nsp
         if(tranp_(is)) WRITE( stdout,517) is, amprp_(is)
      end do
!
      if(tfor_) then
         if(tnosep_) fricp_ = 0.
         WRITE( stdout,520)
         if(tsdp_)then
            WRITE( stdout,521)
         else
            WRITE( stdout,522) fricp_ , greasp_
         endif
      else
         WRITE( stdout,518)
      endif
!
      if( tfor_ ) then
         if(( tcp_ .or. tcap_ .or. tnosep_ ) .and. tsdp_ ) then
            call errore(' main',' t contr. for ions when tsdp=.t.',0)
         endif
         if(.not. tcp_ .and. .not. tcap_ .and. .not. tnosep_ ) then
            WRITE( stdout,550)
         else if(tcp_ .and. tcap_ ) then
            call errore(' main',' tcp and tcap both true',0)
         else if(tcp_ .and. tnosep_ ) then
            call errore(' main',' tcp and tnosep both true',0)
         else if(tcap_ .and. tnosep_ ) then
            call errore(' main',' tcap and tnosep both true',0)
         else if(tcp_ ) then
            WRITE( stdout,555) tempw_ , tolp_
         else if(tcap_) then
            WRITE( stdout,560) tempw_ , tolp_
         else if(tnosep_ ) then
            WRITE( stdout,562) tempw_ , qnp_
         end if
         if(tnosee_) then
            WRITE( stdout,566) ekincw_ , qne_
         end if
      end if
!
      if(tpre_) then
         WRITE( stdout,600)
         if(thdyn_) then
            if(thdiag_) WRITE( stdout,608)
            if(tnoseh_) then
               frich_=0.
               WRITE( stdout,604) temph_,qnh_,press_ / gpa_au
            else
               WRITE( stdout,602) frich_,greash_,press_ / gpa_au
            endif
         else
            WRITE( stdout,606)
         endif
      endif
      if ( agg .ne. 0.d0) then
            WRITE( stdout,650) agg, sgg, e0gg
      end if
      WRITE( stdout,700) iprsta_

!     

 500  format(//                                                         &
     &       ' nbeg=',i3,' nomore=',i7,3x,' iprint=',i4,/               &
     &       ' reads from',i3,' writes on',i3)
 505  format(' time step = ',f9.4/)
 510  format(' parameters for electron dynamics:'/                      &
     &       ' emass= ',f10.2,2x,'emaec= ',f10.2,'ry')
 511  format(' orthog. with lagrange multipliers: eps=',e10.2,          &
     &         ' max=',i3)
 512  format(' orthog. with gram-schmidt')
 513  format(' electron dynamics with steepest descent')
 509  format(' verlet algorithm for electron dynamics')
 514  format(' with friction frice = ',f7.4,' , grease = ',f7.4)
 720  format(' charge density is read from unit 47',/)
 721  format(' charge density is written in unit 47',/)
 722  format(' local potential is written in unit 46',/)
 515  format(' initial random displacement of el. coordinates with ',   &
     &       ' amplitude=',f10.6,/                                      &
     &       ' trane not to be used with mass preconditioning')
 516  format(/)
 517  format(' initial random displacement of ionic coord. for species ',&
     &       i4,' : amplitude=',f10.6)
 518  format(' ions are not allowed to move'/)
 520  format(' ions are allowed to move')
 521  format(' ion dynamics with steepest descent')
 522  format(' ion dynamics with fricp = ',f7.4,' and greasp = ',f7.4)
 550  format(' ion dynamics: the temperature is not controlled'//)
 555  format(' ion dynamics with rescaling of velocities:'/             &
     &       ' temperature required=',f10.5,'(kelvin)',' tolerance=',   &
     &       f10.5//)
 560  format(' ion dynamics with canonical temp. control:'/             &
     &       ' temperature required=',f10.5,'(kelvin)',' tolerance=',   &
     &       f10.5//)
 562  format(' ion dynamics with nose` temp. control:'/                 &
     &       ' temperature required=',f10.5,'(kelvin)',' nose` mass = ',&
     &       f10.3//)
 566  format(' electronic dynamics with nose` temp. control:'/          &
     &       ' elec. kin. en. required=',f10.5,'(hartree)',             &
     &       ' nose` mass = ',f10.3//)
 600  format(' internal stress tensor calculated')
 602  format(' cell parameters dynamics with frich = ',f7.4,            &
     &       ' and greash = ',f7.4,/                                    &
     &       ' external pressure = ',f11.7,'(gpa)'//)
 604  format(' cell parameters dynamics with nose` temp. control:'/     &
     &       ' cell temperature required = ',f10.5,'(kelvin)',          &
     &       ' nose` mass = ',f10.3,/                                   &
     &       ' external pressure = ',f11.7,'(gpa)'//)
 606  format(' cell parameters are not allowed to move'//)
 608  format(' frozen off-diagonal cell parameters'//)
 650  format(' modified kinetic energy functional, with parameters:'/   &
           & ' agg = ',f8.4,'  sgg = ', f7.4,'  e0gg = ',f6.2)
 700  format(' iprsta = ',i2/)

      RETURN

    END SUBROUTINE
!


! ----------------------------------------------------------------

  SUBROUTINE input_info()

    ! this subroutine print to standard output some parameters read from input
    ! ----------------------------------------------

    USE input_parameters, ONLY: restart_mode, nstep, iprint, ndr, ndw, &
      celldm, dt, emass, title

    IMPLICIT NONE

    IF( ionode ) THEN
      WRITE( stdout, * )
      WRITE( stdout, 10)
      WRITE( stdout, 400) title
      WRITE( stdout, 500) restart_mode, nstep, iprint, ndr, ndw
      WRITE( stdout, 501) celldm(1), celldm
      WRITE( stdout, 505) dt
      WRITE( stdout, 510) emass
      WRITE( stdout, 509)
    END IF



    RETURN

 10   FORMAT(//,3X,'MD PARAMETERS READ FROM STANDARD INPUT',/ &
               ,3X,'-------------------------------------')
400   FORMAT(/  3X, 'Job Title: ', A )
500   FORMAT(   3X,'Restart Mode = ',A15,', Number of MD Steps = ',I7,/ &
               ,3X,'Print out every ',I4,' MD Steps',/  &
               ,3X,'Reads from unit = ',I3,', Writes to unit = ',I3)
501   FORMAT( 3X,'Alat   = ', F10.6,/  &
               ,3X,'Celldm = ',6F10.6)
502   FORMAT(   3X,'An initial quench is performed')
505   FORMAT(   3X,'MD Simulation time step    = ',F9.4)
509   FORMAT(   3X,'Verlet algorithm for electron dynamics')
510   FORMAT(   3X,'Electronic fictitious MASS = ',F10.2)



  END SUBROUTINE input_info

! ----------------------------------------------------------------
! ----------------------------------------------------------------

  SUBROUTINE modules_info()

    USE input_parameters, ONLY: electron_dynamics, electron_temperature, &
      orthogonalization

    USE cell_module, ONLY: metric_print_info
    USE cutoffs, ONLY: cutoffs_print_info
    USE ions_module, ONLY: print_scaled_positions, ions_print_info
    USE empty_states, ONLY: empty_print_info
    USE electrons_module, ONLY: electrons_print_info
    USE exchange_correlation, ONLY: exch_corr_print_info
    USE diis, ONLY:  diis_print_info
    USE potentials, ONLY:  potential_print_info
    USE orthogonalize, ONLY: orthogonalize_info
    USE brillouin, ONLY: kpoint_info
    USE runcg_module, ONLY: runcg_info

    IMPLICIT NONE

! ..  Write summary input infos on stdout
      IF( ionode ) THEN
        CALL cutoffs_print_info( stdout )
        CALL electrons_print_info( stdout )
        CALL empty_print_info( stdout )
        IF( TRIM(electron_dynamics) == 'diis' ) THEN
          CALL diis_print_info( stdout )
        ELSE IF( TRIM(electron_dynamics) == 'cg'  ) THEN
          CALL runcg_info( stdout )
        ELSE IF( TRIM(electron_dynamics) == 'sd' ) THEN
          WRITE( stdout,513)
        ELSE
          WRITE( stdout,514)
        END IF
        IF( TRIM(electron_temperature) /= 'nose') THEN
          WRITE( stdout,535)
        ELSE 
          WRITE( stdout,590)
        END IF
        IF( TRIM(orthogonalization) == 'ortho' ) THEN
          CALL orthogonalize_info( stdout )
        ELSE
          WRITE( stdout,512)
        END IF
        CALL kpoint_info( stdout )
        CALL exch_corr_print_info( stdout )
        CALL ions_print_info( stdout )
        CALL potential_print_info( stdout )
        CALL metric_print_info( stdout )
        CALL sic_info( stdout )
      END IF



    RETURN
  512   FORMAT(   3X,'Orthog. with Gram-Schmidt')
  513   FORMAT(   3X,'Electron dynamics with steepest descent')
  514   FORMAT(   3X,'Electron dynamics with newton equations')
  523   FORMAT(   3X,'Dynamic quenching for ions/cell')
  535   FORMAT(   3X,'Electron dynamics : the temperature is not controlled')
  540   FORMAT(   3X,'Electron dynamics with rescaling of velocities :',/ &
               ,3X,'Average kinetic energy required = ',F11.6,'(A.U.)' &
                  ,'Tolerance = ',F11.6)
  545   FORMAT(   3X,'Electron dynamics with canonical temp. control : ',/ &
               ,3X,'Average kinetic energy required = ',F11.6,'(A.U.)' &
                  ,'Tolerance = ',F11.6)
  580   FORMAT(   3X,'Nstepe = ',I3  &
                  ,' purely electronic steepest descent steps',/ &
               ,3X,'are performed for every ionic step in the program')
  590   FORMAT(   3X,'Electron temperature control via nose thermostat')
END SUBROUTINE modules_info



SUBROUTINE sic_info( stdout )
  USE ions_base, ONLY: self_interaction
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: stdout
  !
  ! prints the type of USIC we will do :
  !

  IF( self_interaction == 0 ) THEN
    RETURN
  END IF

        WRITE(stdout, 591)
        WRITE(stdout, 592) self_interaction
        WRITE(stdout, 593)
        select case (self_interaction)
        case (1)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_PZ = - U_hartree[rho_up-rhp_down] - E_excor[rho_up-rho_down,0] '
        case (2)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_MAC = - U_hartree[rho_up-rhp_down] - E_excor[rho_up,rho_down] + E[rho_down, rho_down] '
        case (3)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_basis = - U_hartree[rho_up-rhp_down] '
        case (-1)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_nohartree_PZ =  - E_excor[rho_up-rho_down,0] '
        case (-2)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_nohartree_MAC =  - E_excor[rho_up,rho_down] + E[rho_down, rho_down] '
        case default
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  No unpaired-electron self-interaction correction '
        end select
  591 FORMAT(   3X,'')
  592 FORMAT(   3X,'Introducing a Self_Interaction Correction case: ', I3)
  593 FORMAT(   3X,'----------------------------------------')

  RETURN
END SUBROUTINE

! ----------------------------------------------------------------
  END MODULE input
! ----------------------------------------------------------------
