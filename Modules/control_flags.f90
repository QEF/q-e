!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!=----------------------------------------------------------------------------=!
  MODULE control_flags
!=----------------------------------------------------------------------------=!

!
!  this module contains all basic variables that controls the istuctions
!  execution flow
!  ----------------------------------------------
!
        USE kinds
        USE parameters

        IMPLICIT NONE
        PRIVATE
        SAVE

        TYPE convergence_criteria
          LOGICAL :: active
          INTEGER :: nstep
          REAL(dbl) :: ekin
          REAL(dbl) :: derho
          REAL(dbl) :: force
        END TYPE

        TYPE do_ionic_step
          LOGICAL :: active
          INTEGER :: nstep
          REAL(dbl) :: ekin
        END TYPE

        TYPE electronic_steepest_descent
          LOGICAL :: active
          INTEGER :: nstep
          REAL(dbl) :: ekin
        END TYPE

        TYPE ionic_conjugate_gradient
          LOGICAL :: active
          INTEGER :: nstepix
          INTEGER :: nstepex
          REAL(dbl) :: ionthr
          REAL(dbl) :: elethr
        END TYPE



        PUBLIC :: tbeg, htbeg, nomore, &
                  nbeg, isave, iprint, tv0rd, nv0rd, tzeroc, tzerop, newnfi, tnewnfi, &
                  tfor, tpre, tzeroe, tsde, tsdp, tsdc, taurdr, ndr, &
                  ndw, tortho, tstress, tprnfor, prn, timing, &
                  memchk, tconjgrad, tprnsfac, toptical, &
                  tcarpar, rhoout, trande, arande, trandp, arandp, &
                  tdipole, t_diis, t_diis_simple, t_diis_rot, tnosee, tnosep, tnoseh, &
                  tcp, tcap, tnodump, tdamp, tdampions, tconvthrs, &
                  convergence_criteria, do_ionic_step, tionstep, tsteepdesc, &
                  electronic_steepest_descent, ionic_conjugate_gradient, tconjgrad_ion, &
                  tatomicwfc, tscreen, gamma_only


        PUBLIC ::  control_flags_setup, set_restart_mode, set_verbosity, set_ortho, &
          set_electron_flags, fix_dependencies, check_flags, set_cell_flags, &
          set_ele_steps, set_ion_dynamics, set_ion_flags

        PUBLIC :: tbuff, tvlocw, trhor, trhow, thdyn, iprsta
        


! ...   declare execution control variables

        LOGICAL :: tbuff      ! save wfc on unit 21  ( only cp, never used)
        LOGICAL :: tvlocw     ! write potential to unit 46 ( only cp, seldom used)
        LOGICAL :: trhor      ! read rho from unit 47  (only cp, seldom used)
        LOGICAL :: trhow      ! write rho to  unit 47 (only cp, seldom used)

        LOGICAL :: prn        ! verbosity
        LOGICAL :: tsde       ! electronic steepest descent
        LOGICAL :: tzeroe     ! set to zero the electronic velocities
        LOGICAL :: tfor       ! move the ions ( calculate forces )
        LOGICAL :: tsdp       ! ionic steepest descent
        LOGICAL :: tzerop     ! set to zero the ionic velocities
        LOGICAL :: tprnfor    ! print forces to standard output
        LOGICAL :: taurdr     ! read ionic position from standard input
        LOGICAL :: tv0rd      ! read ionic velocities from standard input
        LOGICAL :: tpre       ! calculate stress, and (in fpmd) variable cell dynamic
        LOGICAL :: thdyn      ! variable-cell dynamics (only cp)
        LOGICAL :: tsdc       ! cell geometry steepest descent
        LOGICAL :: tzeroc     ! set to zero the cell geometry velocities
        LOGICAL :: tstress    ! print stress to standard output
        LOGICAL :: tortho     ! use iterative orthogonalization 
        LOGICAL :: tconjgrad  ! use conjugate gradient electronic minimization
        LOGICAL :: timing     ! print out timing information
        LOGICAL :: memchk     ! check for memory leakage
        LOGICAL :: tprnsfac   ! print out structure factor 
        LOGICAL :: toptical   ! print out optical properties
        LOGICAL :: tcarpar    ! tcarpar is set TRUE for a "pure" Car Parrinello simulation
        LOGICAL :: rhoout     ! print out charge densities
        LOGICAL :: tnodump    ! if true avoid the dump of the wavefunctions in resetart
        LOGICAL :: tdamp      ! Use damped dinamics for electrons
        LOGICAL :: tdampions  ! Use damped dinamics for electrons
        LOGICAL :: tatomicwfc ! Use atomic wavefunctions as starting guess for ch. density
        LOGICAL :: tscreen    ! Use screened coulomb potentials for cluster calculations

        TYPE (convergence_criteria) :: tconvthrs
                              !  thresholds used to check GS convergence 

        TYPE (do_ionic_step) :: tionstep
                              !  parameters to control how many electronic steps 
                              !  between ions move

        TYPE (electronic_steepest_descent) :: tsteepdesc
                              !  parameters for electronic steepest desceent

        TYPE (ionic_conjugate_gradient) :: tconjgrad_ion
                              !  conjugate gradient for ionic minimization

        INTEGER :: nbeg   ! internal code for initialization ( -1, 0, 1, 2, .. )
        INTEGER :: ndw
        INTEGER :: ndr
        INTEGER :: nomore
        INTEGER :: iprint ! print output every iprint step
        INTEGER :: isave  ! write restart to ndr unit every isave step
        INTEGER :: nv0rd
        INTEGER :: iprsta ! output verbosity (increasing from 0 to infinity)

        LOGICAL :: gamma_only !  true if only gamma point is used


        LOGICAL :: tnewnfi = .FALSE.
        INTEGER :: newnfi  = 0

! ...   Wave function randomization
        LOGICAL :: trande
        REAL(dbl) :: arande

! ...   Ionic position randomization
        LOGICAL :: trandp(nsx)
        REAL(dbl) :: arandp(nsx)

! ...   Read the cell from standard input
        LOGICAL :: tbeg = .FALSE.
        REAL(dbl) :: htbeg(3,3)

! ...   This flags control the calculation of the Dipole Moments
        LOGICAL :: tdipole = .FALSE.

! ...   Flags that controls DIIS electronic minimization
        LOGICAL :: t_diis        = .FALSE.
        LOGICAL :: t_diis_simple = .FALSE.
        LOGICAL :: t_diis_rot    = .FALSE.

! ...   Flag controlling the Nose thermostat for electrons
        LOGICAL :: tnosee = .FALSE.

! ...   Flag controlling the Nose thermostat for the cell
        LOGICAL :: tnoseh = .FALSE.

! ...   Flag controlling the Nose thermostat for ions
        LOGICAL :: tnosep = .FALSE.
        LOGICAL :: tcap = .FALSE.
        LOGICAL :: tcp = .FALSE.

!  end of module-scope declarations
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

       SUBROUTINE control_flags_setup(ndr_inp, ndw_inp, iprint_inp, isave_inp, &
         tstress_inp, tprnfor_inp, tdipole_inp, toptical_inp, newnfi_inp,      &
         tnewnfi_inp, gamma_inp )
         IMPLICIT NONE
          INTEGER, INTENT(IN) :: ndr_inp
          INTEGER, INTENT(IN) :: ndw_inp
          INTEGER, INTENT(IN) :: iprint_inp
          INTEGER, INTENT(IN) :: isave_inp
          LOGICAL, INTENT(IN) :: tstress_inp
          LOGICAL, INTENT(IN) :: tprnfor_inp
          LOGICAL, INTENT(IN) :: tdipole_inp
          LOGICAL, INTENT(IN) :: toptical_inp
          INTEGER, INTENT(IN) :: newnfi_inp
          LOGICAL, INTENT(IN) :: tnewnfi_inp
          LOGICAL, INTENT(IN) :: gamma_inp
          ndr     = ndr_inp
          ndw     = ndw_inp
          iprint  = iprint_inp
          isave   = isave_inp
          tstress = tstress_inp
          tprnfor = tprnfor_inp
          tdipole = tdipole_inp
          toptical = toptical_inp
          newnfi = newnfi_inp
          tnewnfi = tnewnfi_inp
          gamma_only = gamma_inp
        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!
        
      SUBROUTINE set_restart_mode(nstep, restart_mode)
! ...   this subroutine set internal flags to restart the job properly
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nstep
        CHARACTER(LEN=80), INTENT(IN) :: restart_mode
        SELECT CASE ( TRIM(restart_mode) )
          CASE ('from_scratch')
            nbeg = -1
            nomore = nstep
          CASE ('reset_counters')
            nbeg = 0
            nomore = nstep
          CASE ('upto')
            nbeg = 1
            nomore = nstep
          CASE ('restart', 'default' )
            nbeg = 2
            nomore = nstep
          CASE DEFAULT
            CALL errore(' control_flag ',' unknown restart_mode '//TRIM(restart_mode), 1 )
        END SELECT
        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

      SUBROUTINE set_verbosity( verbosity, tprnrho )

! ...   this subroutine set the internal flags to the appropriate degree of 
! ...   verbosity

        IMPLICIT NONE
        CHARACTER(LEN=80), INTENT(IN) :: verbosity
        LOGICAL :: tprnrho

        prn     = .FALSE.
        tprnsfac   = .FALSE.
          ! Print on file STRUCTURE_FACTOR the structure factor
          ! gvectors and charge density, in reciprocal space.
        rhoout = .false.
          ! save charge density to file  CHARGEDENSITY if nspin = 1, and
          ! CHARGEDENSITY.UP CHARGEDENSITY.DOWN if nspin = 2
        memchk = .FALSE.
          ! The code performs a memory check, write on standard
          ! output the allocated memory at each step.
          ! Architecture Dependent
        timing = .FALSE.
          ! The code write to files fort.8 fort.41 fort.42 fort.43
          ! a detailed report of subroutines timing

        SELECT CASE ( TRIM(verbosity) )
          CASE ('minimal')
            prn = .FALSE.
          CASE ('low', 'default')
            prn = .FALSE.
            timing = .TRUE.
          CASE ('medium')
            prn = .FALSE.
            timing = .TRUE.
            rhoout = .TRUE.
            tprnsfac = .TRUE.
          CASE ('high')
            prn = .TRUE.
            memchk = .TRUE.
            timing = .TRUE.
            rhoout = .TRUE.
            tprnsfac = .TRUE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown verbosity '//TRIM(verbosity), 1 )
        END SELECT

! ...   If explicitly requested force the charge density to be printed 

        IF( tprnrho ) rhoout = .TRUE.

        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

      SUBROUTINE set_ortho(orthogonalization)
! ...   this subroutine set the orthogonalization method to be used
        IMPLICIT NONE
        CHARACTER(LEN=80), INTENT(IN) :: orthogonalization
        SELECT CASE ( TRIM(orthogonalization) )
          CASE ('Gram-Schmidt')
            tortho = .FALSE.
          CASE ('ortho', 'default')
            tortho = .TRUE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown orthogonalization '//TRIM(orthogonalization), 1 )
        END SELECT
        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

      SUBROUTINE set_electron_flags(electron_dynamics, electron_velocities, electron_temperature, &
        diis_rot, startingwfc, arande_inp)
        IMPLICIT NONE
        CHARACTER(LEN=80), INTENT(IN) :: electron_dynamics
        CHARACTER(LEN=80), INTENT(IN) :: electron_velocities
        CHARACTER(LEN=80), INTENT(IN) :: electron_temperature
        CHARACTER(LEN=80), INTENT(IN) :: startingwfc
        LOGICAL, INTENT(IN) :: diis_rot
        REAL(dbl), INTENT(IN) :: arande_inp
! ... Electronic randomization
        trande = .FALSE.
        tatomicwfc = .FALSE.
        SELECT CASE ( TRIM(startingwfc) )
          CASE ('default','none')
            trande = .FALSE.
          CASE ('random')
            trande = .TRUE.
          CASE ('atomic')
            tatomicwfc = .TRUE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown startingwfc '//TRIM(startingwfc), 1 )
        END SELECT
        arande = arande_inp
! ... Electron dynamics
        tdamp          = .FALSE.
        tconjgrad      = .FALSE.
        tsteepdesc%active = .FALSE.
        tsteepdesc%ekin  = 1.0d+10
        tsteepdesc%nstep = 1
        t_diis        = .FALSE.
        t_diis_simple = .FALSE.
        t_diis_rot    = .FALSE.
        SELECT CASE ( TRIM(electron_dynamics) )
          CASE ('sd', 'default')
            tsde = .TRUE.
          CASE ('verlet')
            tsde = .FALSE.
          CASE ('cg')
            tsde      = .FALSE.
            tconjgrad = .TRUE.
          CASE ('damp')
            tsde   = .FALSE.
            tdamp  = .TRUE.
          CASE ('diis')
            tsde   = .FALSE.
            t_diis = .TRUE.
            IF( diis_rot ) THEN
              t_diis_rot    = .TRUE.
            ELSE
              t_diis_simple = .TRUE.
            END IF
          CASE ('none')
            tsde = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown electron_dynamics '//TRIM(electron_dynamics), 1 )
        END SELECT
! ... Electrons initial velocity
        SELECT CASE ( TRIM(electron_velocities) )
          CASE ('default')
            tzeroe = .FALSE.
          CASE ('zero')
            tzeroe = .TRUE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown electron_velocities '//TRIM(electron_velocities), 1 )
        END SELECT
! ... Electronic Temperature
        tnosee = .FALSE.
        SELECT CASE ( TRIM(electron_temperature) )
!         temperature control of electrons via Nose' thermostat
          CASE ('nose')
            tnosee = .TRUE.
          CASE ('not_controlled', 'default')
            tnosee = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown electron_temperature '//TRIM(electron_temperature), 1 )
        END SELECT
        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

      SUBROUTINE set_cell_flags(cell_dynamics, cell_velocities, cell_parameters, cell_temperature, rd_ht )
        IMPLICIT NONE
        CHARACTER(LEN=80), INTENT(IN) :: cell_dynamics
        CHARACTER(LEN=80), INTENT(IN) :: cell_velocities
        CHARACTER(LEN=80), INTENT(IN) :: cell_parameters
        CHARACTER(LEN=80), INTENT(IN) :: cell_temperature
        REAL(dbl), INTENT(IN) :: rd_ht(3,3)
! ... Cell dynamics
        SELECT CASE ( TRIM(cell_dynamics) )
          CASE ('sd')
            tpre = .TRUE.
            tsdc = .TRUE.
          CASE ('pr')
            tpre = .TRUE.
            tsdc = .FALSE.
          CASE ('none', 'default')
            tpre = .FALSE.
            tsdc = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown cell_dynamics '//TRIM(cell_dynamics), 1 )
        END SELECT
! ... Starting/Restarting Cell parameters
        SELECT CASE ( TRIM(cell_parameters) )
          CASE ('default')
            tbeg = .FALSE.
          CASE ('from_input')
            tbeg = .TRUE.
            htbeg = rd_ht
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown cell_parameters '//TRIM(cell_parameters), 1 )
        END SELECT
! ... Cell initial velocities
        SELECT CASE ( TRIM(cell_velocities) )
          CASE ('default')
            tzeroc = .FALSE.
          CASE ('zero')
            tzeroc = .TRUE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown cell_velocities '//TRIM(cell_velocities), 1 )
        END SELECT
! ... Cell Temperature
        SELECT CASE ( TRIM(cell_temperature) )
!         cell temperature control of ions via Nose' thermostat
          CASE ('nose')
            tnoseh = .TRUE.
            CALL errore(' control_flags ', &
              ' cell_temperature = '//TRIM(cell_temperature)//' not yet implemented ', 1 )
          CASE ('not_controlled', 'default')
            tnoseh = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown cell_temperature '//TRIM(cell_temperature), 1 )
        END SELECT
        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

      SUBROUTINE set_ion_flags(ion_temperature, ion_positions, ion_velocities, &
        trandp_inp, arandp_inp)
        IMPLICIT NONE
        CHARACTER(LEN=80), INTENT(IN) :: ion_temperature
        CHARACTER(LEN=80), INTENT(IN) :: ion_positions
        CHARACTER(LEN=80), INTENT(IN) :: ion_velocities
        LOGICAL, INTENT(IN) :: trandp_inp(:)
        REAL(dbl), INTENT(IN) :: arandp_inp(:)
! ... Ionic Temperature
        tcap         = .FALSE.
        tcp          = .FALSE.
        tnosep       = .FALSE.
        SELECT CASE ( TRIM(ion_temperature) )
!         temperature control of ions via Nose' thermostat
          CASE ('nose')
            tnosep = .TRUE.
          CASE ('not_controlled', 'default')
            tnosep = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown ion_temperature '//TRIM(ion_temperature), 1 )
        END SELECT
! ... Starting/Restarting Atomic positions
        taurdr = .FALSE.
        SELECT CASE ( TRIM(ion_positions) )
          CASE ( 'from_input' )
            taurdr = .TRUE.   ! Positions read from standard input
          CASE ( 'default' )
            taurdr = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown ion_positions '//TRIM(ion_positions), 1 )
        END SELECT
! ... Starting/Restarting ionic velocities
        SELECT CASE ( TRIM(ion_velocities) )
          CASE ('default')
            tzerop = .FALSE.
            tv0rd = .FALSE.
          CASE ('zero')
            tzerop = .TRUE.
            tv0rd = .FALSE.
          CASE ('from_input')
            tzerop = .TRUE.
            tv0rd  = .TRUE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown ion_velocities '//TRIM(ion_velocities), 1 )
        END SELECT
! ... Ionic randomization
        trandp = trandp_inp
        arandp = arandp_inp
        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

! ... Ionic vs Electronic step frequency
! ... When "ion_nstep > 1" and "electron_dynamics = 'md' | 'sd' ", ions are 
! ... propagated every "ion_nstep" electronic step only if the electronic 
! ... "ekin" is lower than "ekin_conv_thr"
!
      SUBROUTINE set_ele_steps(ion_nstepe, cell_nstepe, ekin_conv_thr)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ion_nstepe, cell_nstepe
        REAL(dbl), INTENT(IN) :: ekin_conv_thr
        tionstep%active = .FALSE.
        tionstep%ekin  = 1.0d+10
        tionstep%nstep = 1
        IF( ( ion_nstepe > 1 ) .OR. ( cell_nstepe > 1 ) ) THEN
!         This card is used to control the ionic step, when active ionic step are
!         allowed only when the two criteria are met, i.e. the ions are allowed
!         to move if MOD( NFI, NSTEP ) == 0 and EKIN < EKIN_THR .
          tionstep%active = .TRUE.
          tionstep%nstep  = MAX( ion_nstepe, cell_nstepe )
          tionstep%ekin   = ekin_conv_thr
        END IF
        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

      SUBROUTINE set_ion_dynamics(ion_dynamics, ekin_conv_thr, etot_conv_thr, forc_conv_thr, ion_maxstep, electron_maxstep )
        IMPLICIT NONE
        CHARACTER(LEN=80), INTENT(IN) :: ion_dynamics
        REAL(dbl), INTENT(IN) :: ekin_conv_thr, etot_conv_thr, forc_conv_thr
        INTEGER, INTENT(IN) :: ion_maxstep, electron_maxstep 
! ... Ions dynamics
        tdampions        = .FALSE.
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
            tsdp = .TRUE.
            tfor = .TRUE.
            tconvthrs%ekin   = ekin_conv_thr
            tconvthrs%derho  = etot_conv_thr
            tconvthrs%force  = forc_conv_thr
            tconvthrs%active = .TRUE.
            tconvthrs%nstep  = 1
          CASE ('verlet')
            tsdp = .FALSE.
            tfor = .TRUE.
          CASE ('cg')       ! Conjugate Gradient minimization for ions
            tsdp = .FALSE.
            tfor = .TRUE.
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
            tsdp = .FALSE.
            tfor = .TRUE.
            tdampions = .TRUE.
            tconvthrs%ekin   = ekin_conv_thr
            tconvthrs%derho  = etot_conv_thr
            tconvthrs%force  = forc_conv_thr
            tconvthrs%active = .TRUE.
            tconvthrs%nstep  = 1
          CASE ('none', 'default')
            tsdp = .FALSE.
            tfor = .FALSE.
          CASE DEFAULT
            CALL errore(' control_flags ',' unknown ion_dynamics '//TRIM(ion_dynamics), 1 )
        END SELECT
        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

   SUBROUTINE fix_dependencies()
     IMPLICIT NONE
! ...   Restart file dump
        tnodump = .false.
          ! this flag used for benchmarking and debug
          ! if true Do not save wave functions and other system
          ! properties on the writefile subroutine
        IF( ndw < 0  ) THEN
          tnodump = .true.
        END IF

! ...   Car Parrinello simulation
        tcarpar     = .TRUE.
        IF( tconjgrad .OR. t_diis .OR. tsteepdesc%active ) THEN
          tcarpar = .FALSE.
        END IF

! ...   if TPRE = .FALSE. set TSDC and TZEROC to .FALSE. too.
        IF( .NOT. tpre ) THEN
          tsdc   = .FALSE.
          tzeroc = .FALSE.
        END IF

        IF( .NOT. tfor ) THEN
          tzerop     = .FALSE.
          tv0rd      = .FALSE.
          tsdp       = .FALSE.
          tcp        = .FALSE.
          tcap       = .FALSE.
          tnosep     = .FALSE.
        ELSE
          IF( tsdp ) THEN
            tcp    = .FALSE.
            tcap   = .FALSE.
            tnosep     = .FALSE.
            tv0rd      = .FALSE.
          END IF
          IF( tv0rd ) THEN
            tzerop = .TRUE.
          END IF
        END IF

        IF ( tsde ) THEN
          tnosee = .FALSE.
        END IF

        CALL check_flags()

     RETURN
   END SUBROUTINE fix_dependencies


!=----------------------------------------------------------------------------=!

   SUBROUTINE check_flags()
! ...   do some checks for consistency
        IF ( tnosee .AND. t_diis ) THEN
          CALL errore(' control_flags ','DIIS + ELECT. NOSE ? ',0)
        END IF
        IF ( tortho .AND. t_diis ) THEN
!          CALL errore(' control_flags ','DIIS, ORTHO NOT PERMITTED',0)
        END IF
        IF(tnosep .AND. tcp) THEN
          CALL errore(' control_flags ',' TCP AND TNOSEP BOTH TRUE',0)
        END IF
        IF(tnosep .AND. tcap) THEN
          CALL errore(' control_flags ',' TCAP AND TNOSEP BOTH TRUE',0)
        END IF
        IF(tcp .AND. tcap) THEN
          CALL errore(' control_flags ',' TCP AND TCAP BOTH TRUE',0)
        END IF
        IF( tdipole .AND. tpre ) THEN
          CALL errore('  control_flags  ',' DIPOLE WITH CELL DYNAMICS ',0)
        END IF
        IF( tv0rd .AND. tsdp ) THEN
          CALL errore(' control_flags ',' READING IONS VELOCITY WITH STEEPEST D.',0)
        END IF
     RETURN
   END SUBROUTINE check_flags



!=----------------------------------------------------------------------------=!
  END MODULE control_flags
!=----------------------------------------------------------------------------=!
