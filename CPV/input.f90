!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


MODULE input

   IMPLICIT NONE
   SAVE

   PRIVATE

   PUBLIC :: read_input_file
   PUBLIC :: iosys_pseudo
   PUBLIC :: set_control_flags
   PUBLIC :: modules_setup

   LOGICAL :: has_been_read = .FALSE.

CONTAINS

   !  ----------------------------------------------

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

   subroutine iosys_pseudo( )
        use input_parameters, only:  atom_pfile, pseudo_dir, ntyp
        use control_flags, only:  program_name
        use parameters, only: nsx
        use read_pseudo_module_fpmd, only: readpp
        USE io_files, ONLY: psfile_ => psfile , pseudo_dir_ => pseudo_dir
        USE ions_base, ONLY: nsp_ => nsp

        implicit none

        IF( .NOT. has_been_read ) &
          CALL errore( ' iosys_pseudo ', ' input file has not been read yet! ', 1 )

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

   SUBROUTINE set_control_flags( )
     !
     USE kinds, ONLY: dbl
     USE input_parameters, ONLY: ndr, ndw, iprint, isave, tstress, k_points, &
        tprnfor, verbosity, tprnrho, tdipole_card, toptical_card, &
        tnewnfi_card, newnfi_card, ampre, nstep, restart_mode, ion_positions, &
        startingwfc, printwfc
     USE control_flags, ONLY: &
        ndr_     => ndr, &
        ndw_     => ndw, &
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
        tatomicwfc_ => tatomicwfc, &
        printwfc_ => printwfc
     !
     IMPLICIT NONE

     ndr_  = ndr
     ndw_  = ndw
     iprint_  = iprint
     isave_   = isave
     tstress_ = tstress
     if ( tstress ) tpre_ = .true.
     gamma_only_ = ( TRIM( k_points ) == 'gamma' )
     tprnfor_ = tprnfor
     printwfc_ = printwfc

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

     RETURN
   END SUBROUTINE

   !  ----------------------------------------------
   !  
   !  Initialise modules 

   SUBROUTINE modules_setup( )
     USE electrons_base, ONLY: electrons_base_init
     IMPLICIT NONE
     CALL electrons_base_init
     RETURN
   END SUBROUTINE

END MODULE input

! ----------------------------------------------------------------
!
!
!
!
! ----------------------------------------------------------------

MODULE input_cp

   IMPLICIT NONE
   SAVE

   PRIVATE

   PUBLIC :: iosys

CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------

    subroutine iosys( )

!     this subroutine copies variables from input module to other modules
!     -------------------------------------------------------------------

      use input_parameters, only: &
           nr1, nr2, nr3, greash, press, nr2s, nr3s, nr1s, tolp, temph, grease, &
           tempw, fnoseh, amprp, greasp, tranp, atomic_positions, nelec, &
           if_pos, rd_ht, trd_ht, a, b, c, cosab, cosac, cosbc, cell_symmetry, nelup, &
           neldw, occupations, f_inp, pos, nr3b, pseudo_dir, &
           nr1b, nr2b, sp_pos, atom_mass, atom_pfile, orthogonalization, &
           electron_velocities, startingwfc, ion_dynamics, ion_damping, &
           cell_velocities, electron_dynamics, electron_damping, ion_velocities, &
           celldm, nbnd, nspin, calculation, ntyp, ibrav, restart_mode, ion_positions, &
           ecutwfc, ecutrho, ortho_eps, ortho_max, wmass, qcutz, q2sigma, &
           ecfixed, ekincw, fnosep, nat, disk_io, fnosee, ion_temperature, &
           cell_temperature, cell_dofree, cell_dynamics, cell_damping, electron_temperature, &
           dt, emass, emass_cutoff, ion_radius, &
           ekin_conv_thr, etot_conv_thr, max_seconds, na_inp, rd_pos, atom_label, rd_vel, &
           smd_polm, smd_kwnp, smd_linr, smd_stcd, smd_stcd1, smd_stcd2, smd_stcd3, smd_codf, &
           smd_forf, smd_smwf, smd_lmfreq, smd_tol, smd_maxlm, smd_smcp, smd_smopt, smd_smlm, &
           num_of_images, smd_ene_ini, smd_ene_fin, title, &
           n_inner, fermi_energy, rotmass, occmass, rotation_damping,                       &
           occupation_damping, occupation_dynamics, rotation_dynamics,                      &
           degauss, smearing, tcg, maxiter, etresh, passop, epol, efield


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
      USE check_stop, ONLY: check_stop_init
      !
      USE ions_base, ONLY: ions_base_init, tau_srt, ind_srt, &
           rcmax_ => rcmax, &
           fricp_ => fricp, &
           greasp_ => greasp
      !
      USE ions_positions, ONLY: &
           tau0_ => tau0
      ! 
      USE cell_base, ONLY: cell_base_init, a1, a2, a3, &
           press_ => press, &
           frich_ => frich, &
           greash_ => greash, &
           wmass_ => wmass, &
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
      USE printout_base, ONLY: &
           title_ => title

      USE grid_dimensions, ONLY: &
           nnrx, &  !  variable is used to workaround internal compiler error (IBM xlf)
           nr1_ => nr1, &
           nr2_ => nr2, &
           nr3_ => nr3
      USE smallbox_grid_dimensions, ONLY: &
           nnrbx, &  !  variable is used to workaround internal compiler error (IBM xlf)
           nr1b_ => nr1b, &
           nr2b_ => nr2b, &
           nr3b_ => nr3b
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
      USE efield_module, ONLY: &
           epol_ => epol, &
           efield_ => efield


      USE input, ONLY: set_control_flags, modules_setup

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

      title_ = title   ! simulation title

      CALL set_control_flags( )

      IF( TRIM( calculation ) == 'nscf' ) trhor_ = .true.
     
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
      IF(smd_smlm) THEN
         IF(smd_ene_ini >= 0.d0 .OR. smd_ene_fin >= 0.d0) &
           & CALL errore(' start : ',' Check : ene_ini & ene_fin ', 1 )
      ENDIF
      !
      ene_ini_ = smd_ene_ini
      ene_fin_ = smd_ene_fin

!

      ! ...  Set cell base module

      if( .not. lneb ) then
        CALL cell_base_init( ibrav , celldm , trd_ht, cell_symmetry, rd_ht, &
               a, b, c, cosab, cosac, cosbc , alat_ )
      end if

      ! ...  Set ions base module

      if( .not. lneb ) then
        CALL ions_base_init( ntyp , nat , na_inp , sp_pos , rd_pos , rd_vel, atom_mass, &
             atom_label, if_pos, atomic_positions , alat_ , a1, a2, a3  )
      end if

      ! ...   Set Values for bands and spin

      n_     = nbnd * nspin
      nspin_ = nspin

      ! ...   Set Values for the cutoff

      CALL ecutoffs_setup( ecutwfc, ecutrho, ecfixed, qcutz, q2sigma )



      IF( .NOT. lneb ) THEN
        CALL check_stop_init( max_seconds )
      END IF

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


      SELECT CASE ( electron_velocities ) 
      CASE ('default')
         continue
      CASE ('zero')
         print '("Warning: electron_velocities keyword has no effect")'
      CASE DEFAULT
         CALL errore(' iosys ',' electron_velocities='// &
              trim(electron_velocities)//' not implemented', 1 )
      END SELECT

      ! ...   TSDE

      SELECT CASE ( electron_dynamics ) 
      CASE ('sd')
         tsde_  = .TRUE.
         frice_ = 0.d0
      CASE ('verlet')
         tsde_  = .FALSE.
         frice_ = 0.d0
      CASE ('damp')
         tsde_  = .FALSE.
         frice_ = electron_damping
      CASE ('none')
         tsde_  = .FALSE.
         frice_ = 0.d0
      CASE DEFAULT
         CALL errore(' iosys ',' unknown electron_dynamics '//&
              trim(electron_dynamics),1)
      END SELECT

      SELECT CASE ( electron_velocities )
      CASE ('zero') 
         tzeroe  = .TRUE.
      CASE ('default')
         tzeroe  = .FALSE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown electron_velocities '//&
              trim(electron_dynamics),1)
      END SELECT


      ! Ion velocities

      SELECT CASE ( ion_velocities ) 
      CASE ('default')
         tcap_ = .false.
      CASE ('random')
         tcap_ = .true.
      CASE ('zero')
         tzerop = .TRUE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown ion_velocities '//trim(ion_velocities),1)
      END SELECT

      ! ...   TFOR TSDP

      SELECT CASE ( ion_dynamics ) 
      CASE ('sd')
         tsdp_ = .TRUE.
         tfor_ = .TRUE.
         fricp_= 0.d0
      CASE ('verlet')
         tsdp_ = .FALSE.
         tfor_ = .TRUE.
         fricp_= 0.d0
      CASE ('damp')
         tsdp_ = .FALSE.
         tfor_ = .TRUE.
         fricp_= ion_damping
      CASE ('none')
         tsdp_ = .FALSE.
         tfor_ = .FALSE.
         fricp_= 0.d0
      CASE DEFAULT
         CALL errore(' iosys ',' unknown ion_dynamics '//trim(ion_dynamics), 1 )
      END SELECT

      !

      SELECT CASE ( cell_velocities ) 
      CASE ('default')
        tzeroc = .FALSE.
      CASE ('zero')
        tzeroc = .TRUE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown cell_velocities '//trim(cell_velocities),1)
      END SELECT

      !
     
      ! For SMD 

      IF( lsmd ) THEN

        SELECT CASE ( cell_dynamics ) 
        CASE ('none')
           tpre_ = .FALSE.
           thdyn_= .FALSE.
           frich_= 0.d0
        CASE DEFAULT
          CALL errore(' smiosys ',' cell_dynamics not implemented : '//trim(cell_dynamics), 1 )
        END SELECT

      ELSE

        SELECT CASE ( cell_dynamics )
        CASE ('sd')
           tpre_ = .TRUE.
           thdyn_= .TRUE.
           frich_= 0.d0
        CASE ('pr')
           tpre_ = .TRUE.
           thdyn_= .TRUE.
           frich_= 0.d0
        CASE ('damp-pr')
           tpre_  = .TRUE.
           thdyn_= .TRUE.
           frich_ = cell_damping
        CASE ('none')
           tpre_ = .FALSE.
           thdyn_= .FALSE.
           frich_= 0.d0
        CASE DEFAULT
           CALL errore(' iosys ',' unknown cell_dynamics '//trim(cell_dynamics), 1 )
        END SELECT

      END IF

      !

      SELECT CASE ( electron_temperature ) 
         !         temperature control of electrons via Nose' thermostat
         !         EKINW (REAL(DBL))  average kinetic energy (in atomic units)
         !         FNOSEE (REAL(DBL))  frequency (in terahertz)
      CASE ('nose')
         tnosee_ = .TRUE.
      CASE ('not_controlled')
         tnosee_ = .FALSE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown electron_temperature '//&
              trim(electron_temperature), 1 )
      END SELECT

      !

      SELECT CASE ( ion_temperature ) 
         !         temperature control of ions via Nose' thermostat
         !         TEMPW (REAL(DBL))  frequency (in which units?)
         !         FNOSEP (REAL(DBL))  temperature (in which units?)
      CASE ('nose')
         tnosep_ = .TRUE.
         tcp_ = .false.
      CASE ('not_controlled')
         tnosep_ = .FALSE.
         tcp_ = .false.
      CASE ('rescaling' )
         tnosep_ = .FALSE.
         tcp_ = .true.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown ion_temperature '//&
              trim(ion_temperature), 1 )
      END SELECT


      SELECT CASE ( cell_temperature ) 
         !         cell temperature control of ions via Nose' thermostat
         !         FNOSEH (REAL(DBL))  frequency (in which units?)
         !         TEMPH (REAL(DBL))  temperature (in which units?)
      CASE ('nose')
         tnoseh_ = .TRUE.
      CASE ('not_controlled')
         tnoseh_ = .FALSE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown cell_temperature '//&
              trim(cell_temperature), 1 )
      END SELECT


      SELECT CASE ( cell_dofree )
      CASE ('all')
         thdiag_ =.false.
      CASE ('xyz')
         thdiag_ =.true.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown cell_dofree '//trim(cell_dofree), 1 )
      END SELECT
      if(thdyn_) then
         if(thdiag_) then
            iforceh_=0
            do i=1,3
               iforceh_(i,i)=1
            enddo
         else
            iforceh_=1
         endif
      endif



      ! ...  radii, masses

      DO is = 1, ntyp
         rcmax_ (is) = ion_radius(is)
         IF( ion_radius(is) <= 0.d0 ) THEN
            CALL errore(' iosys ',' invalid  ion_radius ', is) 
         END IF
      END DO

      !
      ! compatibility between FPMD and CP90
      !


      tconvthrs%active = .FALSE.
      IF( ion_dynamics == 'none' .AND. cell_dynamics == 'none' ) THEN
        tconvthrs%ekin   = ekin_conv_thr
        tconvthrs%derho  = etot_conv_thr
        tconvthrs%force  = 10d+10
        tconvthrs%active = .TRUE.
        tconvthrs%nstep  = 1
      END IF

      CALL set_time_step( dt ) 
      emass_ = emass
      emaec_ = emass_cutoff
      ortho_eps_ = ortho_eps
      ortho_max_ = ortho_max


      trhow_ = ( trim( disk_io ) == 'high' )
      tvlocw_ = .false. ! temporaneo
      !
      qne_ = 0.0d0
      qnp_ = 0.0d0
      qnh_ = 0.0d0
      if( fnosee > 0.0d0 ) qne_ = 4.d0*ekincw/(fnosee*(2.d0*pi)*terahertz)**2
      if( fnosep > 0.0d0 ) qnp_ = 2.d0*(3*nat)*tempw/factem/(fnosep*(2.d0*pi)*terahertz)**2
      if( fnoseh > 0.0d0 ) qnh_ = 2.d0*(3*3  )*temph/factem/(fnoseh*(2.d0*pi)*terahertz)**2
      tempw_ = tempw
      temph_ = temph
      ekincw_ = ekincw

      tranp_ ( 1 : ntyp ) =  tranp ( 1 : ntyp )
      amprp_ ( 1 : ntyp ) =  amprp ( 1 : ntyp )
 
      grease_ = grease
      greasp_ = greasp
      tolp_ = tolp
      greash_ = greash
      press_ = press

      nr1_ = nr1
      nr2_ = nr2
      nr3_ = nr3
      nr1s_ = nr1s
      nr2s_ = nr2s
      nr3s_ = nr3s
      nr1b_ = nr1b
      nr2b_ = nr2b
      nr3b_ = nr3b


      !  set pseudopotentials file and directory
      !
      pseudo_dir_ = pseudo_dir
      psfile_= ' '
      psfile_ ( 1:ntyp ) = atom_pfile( 1:ntyp )


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

! ... Set the default value for the cell mass
      wmass_ = wmass
      IF( wmass_ == 0.d0 ) THEN
        wmass_ = 3.d0 / (4.d0 * pi**2 ) * SUM( atom_mass(1:ntyp)*na_inp(1:ntyp) )
        wmass_ = wmass_ * UMA_AU
        WRITE( stdout,999) wmass_
      ELSE
        WRITE( stdout,998) wmass_
      END IF
 998  format('   wmass (read from input) = ',f15.2,/)
 999  format('   wmass (calculated) = ',f15.2,/)


      !
      !
      !  set occupancies
      !
      
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

      !------------------Gradiente Coniugato------------------
      !
      tcg_=tcg
      maxiter_=maxiter
      etresh_=etresh
      passop_=passop
      !
      !---------------campo elettrico
      !
      epol_=epol
      efield_=efield

!
!     --------------------------------------------------------
!     print out heading
!
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
      do is =1, ntyp
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
               WRITE( stdout,604) temph_,qnh_,press_
            else
               WRITE( stdout,602) frich_,greash_,press_
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
      !  Modules setup     
      !     

      CALL modules_setup()
!     
!     
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
      return
      end subroutine
!

END MODULE input_cp


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni
!  ----------------------------------------------
!  BEGIN manual

MODULE input_fpmd

!  this module handles the reading of input data
!  ----------------------------------------------
!
!  END manual
! ----------------------------------------------------------------------

        USE kinds, ONLY: dbl
        USE constants, ONLY: UMA_AU, pi
        USE io_global, ONLY: ionode, stdout

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: iosys

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE iosys

        USE control_flags, ONLY: fix_dependencies
        USE input, ONLY: modules_setup

        IMPLICIT NONE

! . Set internal flags according to the input .......................!

        CALL set_internal_flags_fpmd()

! . Fix values for dependencies .....................................!

        CALL fix_dependencies()

! . Write to stdout input module information ........................!

        CALL input_info()

! . CALL the Module specific setup routine ..........................! 

        CALL modules_setup_fpmd()
        CALL modules_setup()

! . Write to stdout module specifice information  ...................!

        CALL modules_info()

        RETURN
      END SUBROUTINE iosys


! ----------------------------------------------------------------
! ----------------------------------------------------------------


      SUBROUTINE modules_setup_fpmd()

        USE input_parameters, ONLY: &
          tturbo_inp, nturbo_inp, diis_achmix, electron_damping, vhrmax_inp, &
          vhasse_inp, iesr_inp, diis_ethr, diis_wthr, diis_delt, &
          etot_conv_thr, diis_g0chmix, diis_nchmix, diis_g1chmix, xc_type, &
          empty_states_maxstep, empty_states_delt, &
          empty_states_emass, empty_states_ethr, vhnr_inp, vhiunit_inp, &
          vhrmin_inp, tvhmean_inp, tscra_inp, scradir, max_seconds, wk, &
          outdir, prefix, k3, nk3, k1, k2, woptical, noptical, boptical, &
          electron_maxstep, tprnks, tprnks_empty, ks_path, diis_nreset, &
          diis_temp, diis_nrot, diis_maxstep, diis_fthr, diis_size, diis_hcut, &
          k_points, nk1, nk2, diis_rothr, tchi2_inp, ortho_max, ortho_eps, &
          ekin_conv_thr, na_inp, ibrav, press, atom_mass, nkstot, xk, wmass, &
          nbnd, nelec, nelup, tf_inp, cell_dofree, t2dpegrid_inp, &
          diis_chguess, tpstab_inp, pstab_size_inp, ion_radius, atom_pfile, &
          dt, ntyp, pseudo_dir, qcutz, q2sigma, tk_inp, ecfixed, celldm, &
          ecutwfc, ecutrho, constr_dist_inp, constr_inp, ion_damping, &
          constr_type_inp, anner_inp, nconstr_inp, constr_tol_inp, tolp, &
          fnosee, ekincw, tempw, tneighbo, neighbo_radius, fnosep, &
          emass_cutoff, f_inp, sp_pos, emass, neldw, nspin, empty_states_nbnd, &
          atom_label, anne_inp, rd_vel, rd_pos, if_pos, sp_vel, rd_ht, &
          a, b, c, cosab, cosac, cosbc, cell_symmetry, trd_ht, nat, id_loc, &
          sic, sic_epsilon, sic_rloc, atomic_positions, cell_damping, greash, &
          fnoseh, temph, nr1b, nr2b, nr3b, title

        USE io_files, ONLY: &
          psfile, &
          pseudo_dir_ => pseudo_dir, &
          prefix_     => prefix, &
          scradir_   => scradir

        USE control_flags, ONLY: program_name, lneb, tnoseh, &
              tolp_      => tolp, &
              ortho_max_ => ortho_max, &
              ortho_eps_ => ortho_eps

        USE nose_ions, ONLY: nose_ions_setup
        USE nose_electrons, ONLY: nose_electrons_setup
        USE time_step, ONLY: set_time_step
        USE cell_module, ONLY: metric_setup
        USE ions_module, ONLY: ions_setup, tc_ions_setup
        USE guess, ONLY: guess_setup
        USE real_space_mesh, ONLY: real_space_mesh_setup
        USE bands_mesh, ONLY: bands_mesh_setup
        USE empty_states, ONLY: empty_setup
        USE electrons_module, ONLY: electrons_setup
        USE pseudopotential, ONLY: pseudopotential_setup
        USE exchange_correlation, ONLY: exch_corr_setup
        USE turbo, ONLY: turbo_setup
        USE diis, ONLY: diis_setup
        USE charge_mix, ONLY: charge_mix_setup
        USE chi2, ONLY: chi2_setup
        USE potentials, ONLY:  potential_setup
        USE wave_base, ONLY: wave_base_init
        USE orthogonalize, ONLY: orthogonalize_setup
        USE stress, ONLY: stress_setup
        USE brillouin, ONLY: kpoint_setup
        USE print_out_module, ONLY: printout_setup
        USE kohn_sham_states, ONLY: ks_states_setup
        USE runcg_module, ONLY: runcg_setup
        USE optical_properties, ONLY: optical_setup
        USE reciprocal_space_mesh, ONLY: recvecs_units
        USE cell_base, ONLY: cell_base_init, a1, a2, a3
        USE cell_nose, ONLY: cell_nose_init
        USE ions_base, ONLY: ions_base_init, self_interaction
        USE ions_nose, ONLY: ions_noseinit
        USE check_stop, ONLY: check_stop_init
        USE constants, ONLY: factem, terahertz, pi, eps8
        USE smallbox_grid_dimensions, ONLY: &
           nnrbx, &  !  variable is used to workaround internal compiler error (IBM xlf)
           nr1b_ => nr1b, nr2b_ => nr2b, nr3b_ => nr3b
        !
        USE printout_base, ONLY: &
              title_ => title


        IMPLICIT NONE

! ...   Declare Variables
        REAL(dbl)  :: delt_emp_inp, emass_emp_inp, ethr_emp_inp
! ...   DIIS
        REAL(dbl) :: tol_diis_inp, delt_diis_inp, tolene_inp
        REAL(dbl) :: alat_ 
        LOGICAL :: o_diis_inp, oqnr_diis_inp
        INTEGER :: icapstep
        LOGICAL :: timing = .TRUE.

! ...   end of declarations
!  ----------------------------------------------

        title_ = title     !  simulation title

        !  set module variables
  
        psfile( 1 : ntyp ) = atom_pfile( 1 : ntyp )
        pseudo_dir_        = pseudo_dir
        prefix_            = prefix

! ...   auxiliary input and/or setup routines

        IF( .not. lneb ) THEN
          CALL cell_base_init( ibrav , celldm , trd_ht, cell_symmetry, rd_ht, &
               a, b, c, cosab, cosac, cosbc , alat_ )
        END IF

        CALL set_time_step( dt )
        CALL pseudopotential_setup( ntyp, pseudo_dir, atom_pfile, tpstab_inp, pstab_size_inp, ion_radius )
        CALL recvecs_units( alat_ )

        ! cut offs

        CALL ecutoffs_setup( ecutwfc, ecutrho, ecfixed, qcutz, q2sigma )

        CALL gcutoffs_setup( alat_ , tk_inp, nkstot, xk)

! ... Set the default value for the cell mass
        IF( wmass == 0.d0 ) THEN  
          wmass = 3.d0 / (4.d0 * pi**2 ) * SUM( atom_mass(1:ntyp)*na_inp(1:ntyp) )
          wmass = wmass * UMA_AU
        END IF


        CALL metric_setup( wmass, press, cell_damping, greash, cell_dofree )

        ! set thermostat parameter for cell
        CALL cell_nose_init( temph, fnoseh )

        CALL real_space_mesh_setup( t2dpegrid_inp )

        nr1b_ = nr1b  ! set box grid module variables
        nr2b_ = nr2b
        nr3b_ = nr3b

        CALL guess_setup( diis_chguess )

        CALL bands_mesh_setup( )

        CALL electrons_setup( tf_inp, nbnd, nint(nelec), nint(nelup), &
          nint(neldw), nspin, empty_states_nbnd, emass, emass_cutoff, &
          f_inp, nkstot )

        ! call errore( ' qui ', ' qui ', 1 )

        IF( .not. lneb ) THEN
          CALL ions_base_init( ntyp , nat , na_inp , sp_pos , rd_pos , rd_vel, atom_mass, &
               atom_label, if_pos, atomic_positions, alat_ , a1 , a2 , a3 , &
               id_loc, sic, sic_epsilon, sic_rloc  )
        END IF

        CALL ions_setup( anne_inp, anner_inp,   &
             nconstr_inp, constr_tol_inp, constr_type_inp, constr_dist_inp, &
             constr_inp, ion_damping, tneighbo, neighbo_radius )

        IF( ABS(self_interaction) /= 0 ) CALL sic_info()

        CALL nose_ions_setup( fnosep, tempw )
        CALL ions_noseinit( tempw, fnosep, nat )

        icapstep = 1
        CALL tc_ions_setup(tempw, icapstep)
        tolp_ = tolp
        !
        CALL nose_electrons_setup( fnosee, ekincw )

        delt_emp_inp  = dt
        emass_emp_inp = emass
        ethr_emp_inp  = ekin_conv_thr
        IF( empty_states_delt > 0.d0 )  delt_emp_inp  = empty_states_delt
        IF( empty_states_emass > 0.d0 ) emass_emp_inp = empty_states_emass
        IF( empty_states_ethr > 0.d0 )  ethr_emp_inp  = empty_states_ethr
        CALL empty_setup( empty_states_maxstep, delt_emp_inp, ethr_emp_inp )

        CALL exch_corr_setup( xc_type )
        CALL check_stop_init( max_seconds )
        !
        scradir_ = TRIM( scradir )
        !
        CALL potential_setup(tvhmean_inp,vhnr_inp, vhiunit_inp, &
               vhrmin_inp, vhrmax_inp, vhasse_inp, timing, iesr_inp)
        CALL wave_base_init( electron_damping )
        CALL turbo_setup(tturbo_inp, nturbo_inp)
        CALL charge_mix_setup(diis_achmix, diis_g0chmix, diis_nchmix, diis_g1chmix)

        o_diis_inp        = .TRUE.
        oqnr_diis_inp     = .TRUE.
        tolene_inp = etot_conv_thr
        tol_diis_inp = ekin_conv_thr
        delt_diis_inp = dt
        IF( diis_ethr > 0.0d0 ) tolene_inp = diis_ethr
        IF( diis_wthr > 0.0d0 ) tol_diis_inp = diis_wthr
        IF( diis_delt > 0.0d0 ) delt_diis_inp = diis_delt
        CALL diis_setup( diis_fthr, oqnr_diis_inp, o_diis_inp, &
          diis_size, diis_hcut, tol_diis_inp, diis_maxstep, diis_nreset, delt_diis_inp, &
          diis_temp, diis_nrot(1), diis_nrot(2), diis_nrot(3), &
          diis_rothr(1), diis_rothr(2), diis_rothr(3), tolene_inp)

        CALL chi2_setup(tchi2_inp)
        !
        CALL orthogonalize_setup( timing )
        ortho_max_ = ortho_max
        ortho_eps_ = ortho_eps
        !
        CALL stress_setup(timing)
        CALL kpoint_setup(k_points, nkstot, nk1, nk2, nk3, k1, k2, k3, xk, wk)
        CALL printout_setup( outdir, prefix )
        CALL ks_states_setup(nspin, tprnks, tprnks_empty, ks_path)
        CALL runcg_setup( etot_conv_thr, electron_maxstep )
        CALL optical_setup( woptical, noptical, boptical )

        RETURN
      END SUBROUTINE modules_setup_fpmd



! ----------------------------------------------------------------
! ----------------------------------------------------------------

    SUBROUTINE set_internal_flags_fpmd()

      USE input_parameters, ONLY: &
        cell_velocities, &
        cell_dynamics, cell_parameters, cell_temperature, trd_ht, &
        tapos, tavel, ecutwfc, emass_cutoff, taspc,  &
        electron_temperature, &
        electron_velocities, diis_rot, startingwfc, electron_dynamics, &
        orthogonalization, restart_mode, tranp, &
        ion_velocities, amprp, cell_nstepe, ion_nstepe, ion_positions, &
        ekin_conv_thr, ion_dynamics, ion_maxstep, ion_temperature, &
        electron_maxstep, etot_conv_thr, forc_conv_thr, ibrav
      USE input_parameters, ONLY: force_pairing_ => force_pairing

      USE input, ONLY: set_control_flags

      USE control_flags, ONLY: &
        tortho_  => tortho

      USE control_flags, ONLY: &
        t_diis_simple_ => t_diis_simple, &
        t_diis_ => t_diis, &
        tsde_ => tsde, &
        t_diis_rot_ => t_diis_rot, &
        ekin_maxiter_ => ekin_maxiter, &
        tconjgrad_ => tconjgrad, &
        tdamp_ => tdamp, &
        ekin_conv_thr_ => ekin_conv_thr, &
        tsteepdesc_ => tsteepdesc, &
        tzeroe_ => tzeroe, &
        tnosee_ => tnosee

      USE control_flags, ONLY: &
        tzeroc_ => tzeroc, &
        tnoseh_ => tnoseh, &
        thdyn_ => thdyn, &
        tsdc_ => tsdc, &
        tbeg_ => tbeg

      USE control_flags, ONLY: &
        tionstep_ => tionstep, &
        nstepe_ => nstepe

      USE control_flags, ONLY: &
        tdampions_ => tdampions, &
        tfor_ => tfor, &
        tsdp_ => tsdp, &
        tconvthrs, tconjgrad_ion

      USE control_flags, ONLY: &
        tnosep_ => tnosep, &
        tcap_ => tcap, &
        tcp_ => tcp, &
        tzerop_ => tzerop, &
        tv0rd_ => tv0rd, &
        tranp_ => tranp, &
        amprp_ => amprp

      USE control_flags, ONLY: &
        force_pairing

      IMPLICIT NONE


! . Set internal flags according to the input .......................!

        CALL set_control_flags( )

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

        ! ... Electron dynamics

        tdamp_          = .FALSE.
        tconjgrad_      = .FALSE.
        tsteepdesc_     = .FALSE.
        ekin_conv_thr_  = 1.0d+10
        ekin_maxiter_ = 1
        t_diis_        = .FALSE.
        t_diis_simple_ = .FALSE.
        t_diis_rot_    = .FALSE.
        SELECT CASE ( TRIM(electron_dynamics) )
          CASE ('sd', 'default')
            tsde_ = .TRUE.
          CASE ('verlet')
            tsde_ = .FALSE.
          CASE ('cg')
            tsde_      = .FALSE.
            tconjgrad_ = .TRUE.
          CASE ('damp')
            tsde_   = .FALSE.
            tdamp_  = .TRUE.
          CASE ('diis')
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

        ! ... Electrons initial velocity

        SELECT CASE ( TRIM(electron_velocities) )
          CASE ('default')
            tzeroe_ = .FALSE.
          CASE ('zero')
            tzeroe_ = .TRUE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown electron_velocities '//TRIM(electron_velocities), 1 )
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
            tsdp_ = .TRUE.
            tfor_ = .TRUE.
            tconvthrs%ekin   = ekin_conv_thr
            tconvthrs%derho  = etot_conv_thr
            tconvthrs%force  = forc_conv_thr
            tconvthrs%active = .TRUE.
            tconvthrs%nstep  = 1
          CASE ('verlet')
            tsdp_ = .FALSE.
            tfor_ = .TRUE.
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
          CASE ('damp')
            tsdp_ = .FALSE.
            tfor_ = .TRUE.
            tdampions_ = .TRUE.
            tconvthrs%ekin   = ekin_conv_thr
            tconvthrs%derho  = etot_conv_thr
            tconvthrs%force  = forc_conv_thr
            tconvthrs%active = .TRUE.
            tconvthrs%nstep  = 1
          CASE ('none', 'default')
            tsdp_ = .FALSE.
            tfor_ = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown ion_dynamics '//TRIM(ion_dynamics), 1 )
        END SELECT

        ! ... Ionic Temperature

        tcap_         = .FALSE.
        tcp_          = .FALSE.
        tnosep_       = .FALSE.
        SELECT CASE ( TRIM(ion_temperature) )
          !         temperature control of ions via Nose' thermostat
          CASE ('nose')
            tnosep_ = .TRUE.
          CASE ('not_controlled', 'default')
            tnosep_ = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown ion_temperature '//TRIM(ion_temperature), 1 )
        END SELECT
        ! ... Starting/Restarting ionic velocities
        SELECT CASE ( TRIM(ion_velocities) )
          CASE ('default')
            tzerop_ = .FALSE.
            tv0rd_ = .FALSE.
          CASE ('zero')
            tzerop_ = .TRUE.
            tv0rd_ = .FALSE.
          CASE ('from_input')
            tzerop_ = .TRUE.
            tv0rd_  = .TRUE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown ion_velocities '//TRIM(ion_velocities), 1 )
        END SELECT
        ! ... Ionic randomization
        tranp_ = tranp
        amprp_ = amprp


        tionstep_ = .FALSE.
        nstepe_   = 1
        IF( ( ion_nstepe > 1 ) .OR. ( cell_nstepe > 1 ) ) THEN
!         This card is used to control the ionic step, when active ionic step are
!         allowed only when the two criteria are met, i.e. the ions are allowed
!         to move if MOD( NFI, NSTEP ) == 0 and EKIN < EKIN_THR .
          tionstep_ = .TRUE.
          nstepe_   = MAX( ion_nstepe, cell_nstepe )
        END IF
        ekin_conv_thr_   = ekin_conv_thr

           
        SELECT CASE ( TRIM(cell_dynamics) )
          CASE ('sd')
            thdyn_ = .TRUE.
            tsdc_ = .TRUE.
          CASE ('damp')
            thdyn_ = .TRUE.
            tsdc_ = .FALSE.
          CASE ('pr')
            thdyn_ = .TRUE.
            tsdc_ = .FALSE.
          CASE ('none', 'default')
            thdyn_ = .FALSE.
            tsdc_ = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown cell_dynamics '//TRIM(cell_dynamics), 1 )
        END SELECT
! ... Starting/Restarting Cell parameters
        SELECT CASE ( TRIM(cell_parameters) )
          CASE ('default')
            tbeg_ = .FALSE.
          CASE ('from_input')
            tbeg_ = .TRUE.
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
        force_pairing = force_pairing_

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
   END SUBROUTINE set_internal_flags_fpmd



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

    USE input_parameters, ONLY: electron_dynamics, &
      electron_velocities, electron_temperature, &
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


SUBROUTINE sic_info()
  USE io_global, ONLY: ionode, stdout
  USE ions_base, ONLY: self_interaction
  IMPLICIT NONE
! prints the type of USIC we will do :
      if (ionode) then
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
      endif
  591 FORMAT(   3X,'')
  592 FORMAT(   3X,'Introducing a Self_Interaction Correction case: ', I3)
  593 FORMAT(   3X,'----------------------------------------')
END SUBROUTINE

!=----------------------------------------------------------------------------=!
END MODULE input_fpmd
!=----------------------------------------------------------------------------=!
