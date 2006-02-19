!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE control_flags
  !=--------------------------------------------------------------------------=!
  !
  ! ... this module contains all basic variables that controls the istuctions
  ! ... execution flow
  !----------------------------------------------
  !
  USE kinds
  USE parameters
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  TYPE convergence_criteria
     !
     LOGICAL  :: active
     INTEGER  :: nstep
     REAL(DP) :: ekin
     REAL(DP) :: derho
     REAL(DP) :: force
     !
  END TYPE convergence_criteria
  !
  TYPE ionic_conjugate_gradient
     !
     LOGICAL  :: active
     INTEGER  :: nstepix
     INTEGER  :: nstepex
     REAL(DP) :: ionthr
     REAL(DP) :: elethr
     !
  END TYPE ionic_conjugate_gradient
  !
  PUBLIC :: tbeg, nomore, nbeg, isave, iprint, tv0rd, nv0rd, tzeroc, tzerop, &
            newnfi, tnewnfi, tfor, tpre, tzeroe, tsde, tsdp, tsdc, taurdr,   &
            ndr, ndw, tortho, ortho_eps, ortho_max, tstress, tprnfor,        &
            timing, memchk, tconjgrad, tprnsfac, toptical, tcarpar,          &
            trane, ampre, tranp, amprp, tdipole, t_diis, t_diis_simple,      &
            t_diis_rot, tnosee, tnosep, tnoseh, tcp, tcap, tdamp, tdampions, &
            tconvthrs, tolp, convergence_criteria, tionstep, nstepe,         &
            tsteepdesc, ionic_conjugate_gradient, tconjgrad_ion, tatomicwfc, &
            tscreen, gamma_only, force_pairing, tchi2
  !
  PUBLIC :: fix_dependencies, check_flags
  PUBLIC :: tbuff, tvlocw, trhor, thdyn, iprsta
  PUBLIC :: twfcollect, printwfc
  PUBLIC :: tuspp
  PUBLIC :: program_name
  !
  ! ...   declare execution control variables
  !
  CHARACTER(LEN=4) :: program_name = ' '  !  used to control execution flow inside module
                                          !  main_program = 'PWSC'  pwscf
                                          !                 'FPMD'  fpmd
                                          !                 'CPVC'  cp
                                          !                 'SMCP'  smcp
  !
  LOGICAL :: tbuff     = .FALSE. ! save wfc on        unit 21 (only cp, never used)
  LOGICAL :: tvlocw    = .FALSE. ! write potential to unit 46 (only cp, seldom used)
  LOGICAL :: trhor     = .FALSE. ! read rho from      unit 47 (only cp, seldom used)
  !
  LOGICAL :: tsde          = .FALSE. ! electronic steepest descent
  LOGICAL :: tzeroe        = .FALSE. ! set to zero the electronic velocities
  LOGICAL :: tfor          = .FALSE. ! move the ions ( calculate forces )
  LOGICAL :: tsdp          = .FALSE. ! ionic steepest descent
  LOGICAL :: tzerop        = .FALSE. ! set to zero the ionic velocities
  LOGICAL :: tprnfor       = .FALSE. ! print forces to standard output
  LOGICAL :: taurdr        = .FALSE. ! read ionic position from standard input
  LOGICAL :: tv0rd         = .FALSE. ! read ionic velocities from standard input
  LOGICAL :: tpre          = .FALSE. ! calculate stress, and (in fpmd) variable cell dynamic
  LOGICAL :: thdyn         = .FALSE. ! variable-cell dynamics (only cp)
  LOGICAL :: tsdc          = .FALSE. ! cell geometry steepest descent
  LOGICAL :: tzeroc        = .FALSE. ! set to zero the cell geometry velocities
  LOGICAL :: tstress       = .FALSE. ! print stress to standard output
  LOGICAL :: tortho        = .FALSE. ! use iterative orthogonalization 
  LOGICAL :: tconjgrad     = .FALSE. ! use conjugate gradient electronic minimization
  LOGICAL :: timing        = .FALSE. ! print out timing information
  LOGICAL :: memchk        = .FALSE. ! check for memory leakage
  LOGICAL :: tprnsfac      = .FALSE. ! print out structure factor 
  LOGICAL :: toptical      = .FALSE. ! print out optical properties
  LOGICAL :: tcarpar       = .FALSE. ! tcarpar is set TRUE for a "pure" Car Parrinello simulation
  LOGICAL :: tdamp         = .FALSE. ! Use damped dinamics for electrons
  LOGICAL :: tdampions     = .FALSE. ! Use damped dinamics for electrons
  LOGICAL :: tatomicwfc    = .FALSE. ! Use atomic wavefunctions as starting guess for ch. density
  LOGICAL :: tscreen       = .FALSE. ! Use screened coulomb potentials for cluster calculations
  LOGICAL :: twfcollect    = .FALSE. ! Collect wave function in the restart file at the end of run.
  LOGICAL :: tuspp         = .FALSE. ! Ultra-soft pseudopotential are being used
  INTEGER :: printwfc      = -1      ! Print wave functions, temporarely used only by ensemble-dft
  LOGICAL :: force_pairing = .FALSE. ! ...   Force pairing
  LOGICAL :: tchi2         = .FALSE. ! Compute Chi^2
  !
  TYPE (convergence_criteria) :: tconvthrs
                              !  thresholds used to check GS convergence 
  !
  ! ... Ionic vs Electronic step frequency
  ! ... When "ion_nstep > 1" and "electron_dynamics = 'md' | 'sd' ", ions are 
  ! ... propagated every "ion_nstep" electronic step only if the electronic 
  ! ... "ekin" is lower than "ekin_conv_thr"
  !
  LOGICAL :: tionstep = .FALSE.
  INTEGER :: nstepe   = 1  
                            !  parameters to control how many electronic steps 
                            !  between ions move

  LOGICAL :: tsteepdesc = .FALSE.
                            !  parameters for electronic steepest desceent

  TYPE (ionic_conjugate_gradient) :: tconjgrad_ion
                            !  conjugate gradient for ionic minimization

  INTEGER :: nbeg   = 0 ! internal code for initialization ( -1, 0, 1, 2, .. )
  INTEGER :: ndw    = 0 !
  INTEGER :: ndr    = 0 !
  INTEGER :: nomore = 0 !
  INTEGER :: iprint = 0 ! print output every iprint step
  INTEGER :: isave  = 0 ! write restart to ndr unit every isave step
  INTEGER :: nv0rd  = 0 !
  INTEGER :: iprsta = 0 ! output verbosity (increasing from 0 to infinity)
  !
  ! ... .TRUE. if only gamma point is used
  !
  LOGICAL :: gamma_only = .TRUE.
  !
  LOGICAL :: tnewnfi = .FALSE.
  INTEGER :: newnfi  = 0
  !
  ! ... Wave function randomization
  !
  LOGICAL  :: trane = .FALSE.
  REAL(DP) :: ampre = 0.D0
  !
  ! ... Ionic position randomization
  !
  LOGICAL  :: tranp(nsx) = .FALSE.
  REAL(DP) :: amprp(nsx) = 0.D0
  !
  ! ... Read the cell from standard input
  !
  LOGICAL :: tbeg = .FALSE.
  !
  ! ... This flags control the calculation of the Dipole Moments
  !
  LOGICAL :: tdipole = .FALSE.
  !
  ! ... Flags that controls DIIS electronic minimization
  !
  LOGICAL :: t_diis        = .FALSE.
  LOGICAL :: t_diis_simple = .FALSE.
  LOGICAL :: t_diis_rot    = .FALSE.
  !
  ! ... Flag controlling the Nose thermostat for electrons
  !
  LOGICAL :: tnosee = .FALSE.
  !
  ! ... Flag controlling the Nose thermostat for the cell
  !
  LOGICAL :: tnoseh = .FALSE.
  !
  ! ... Flag controlling the Nose thermostat for ions
  !
  LOGICAL  :: tnosep = .FALSE.
  LOGICAL  :: tcap   = .FALSE.
  LOGICAL  :: tcp    = .FALSE.
  REAL(DP) :: tolp   = 0.D0   !  tolerance for temperature variation
  !
  REAL(DP), PUBLIC :: &
       ekin_conv_thr = 0.D0, &!  conv. threshold for fictitious e. kinetic energy
       etot_conv_thr = 0.D0, &!  conv. threshold for DFT energy
       forc_conv_thr = 0.D0   !  conv. threshold for atomic forces
  INTEGER, PUBLIC :: &
       ekin_maxiter = 100,   &!  max number of iter. for ekin convergence
       etot_maxiter = 100,   &!  max number of iter. for etot convergence
       forc_maxiter = 100     !  max number of iter. for atomic forces conv.
  !
  ! ... Several variables controlling the run ( used mainly in PW calculations )
  !
  ! ... logical flags controlling the execution
  !
  LOGICAL, PUBLIC :: &
    lfixatom,         &! if .TRUE. some atom is kept fixed
    lscf,             &! if .TRUE. the calc. is selfconsistent
    lbfgs,            &! if .TRUE. the calc. is a relaxation based on new BFGS scheme
    lmd,              &! if .TRUE. the calc. is a dynamics
    lmetadyn,         &! if .TRUE. the calc. is a meta-dynamics
    lpath,            &! if .TRUE. the calc. is a path optimizations
    lneb,             &! if .TRUE. the calc. is NEB dynamics
    lsmd,             &! if .TRUE. the calc. is string dynamics
    lwf,              &! if .TRUE. the calc. is with wannier functions
    lphonon,          &! if .TRUE. the calc. is phonon
    lraman,           &! if .TRUE. the calc. is raman
    lconstrain,       &! if .TRUE. the calc. is constraint
    ldamped,          &! if .TRUE. the calc. is a damped dynamics
    lrescale_t,       &! if .TRUE. the ionic temperature is rescaled
    lcoarsegrained,   &! if .TRUE. a coarse-grained phase-space is used
    restart            ! if .TRUE. restart from results of a preceding run
  !
  ! ... pw self-consistency
  !
  INTEGER, PUBLIC :: &
    ngm0,             &! used in mix_rho
    niter,            &! the maximum number of iteration
    nmix,             &! the number of iteration kept in the history
    imix               ! the type of mixing (0=plain,1=TF,2=local-TF)
  REAL(DP), PUBLIC  :: &
    mixing_beta,      &! the mixing parameter
    tr2                ! the convergence threshold for potential
  LOGICAL, PUBLIC :: &
    conv_elec          ! if .TRUE. electron convergence has been reached
  !
  ! ... pw diagonalization
  !
  REAL(DP), PUBLIC  :: &
    ethr               ! the convergence threshold for eigenvalues  
  INTEGER, PUBLIC :: &
    david,            &! used on Davidson diagonalization
    isolve,           &! Davidson or CG or DIIS diagonalization
    max_cg_iter,      &! maximum number of iterations in a CG di
    diis_buff,        &! dimension of the buffer in diis
    diis_ndim          ! dimension of reduced basis in DIIS
  !
  ! ... wfc and rho extrapolation
  !
  REAL(DP), PUBLIC  :: &
    alpha0,           &! the mixing parameters for the extrapolation
    beta0              ! of the starting potential
  INTEGER, PUBLIC :: &
    history,          &! number of old steps available for potential updating
    pot_order,        &! type of potential updating ( see update_pot )
    wfc_order          ! type of wavefunctions updating ( see update_pot )
  !
  ! ... ionic dynamics
  !
  INTEGER, PUBLIC :: &
    nstep,            &! number of ionic steps
    istep              ! current ionic step  
  LOGICAL, PUBLIC :: &
    conv_ions          ! if .TRUE. ionic convergence has been reached
  REAL(DP), PUBLIC  :: &
    upscale            ! maximum reduction of convergence threshold
  !
  ! ... system's symmetries
  !
  LOGICAL, PUBLIC :: &
    nosym,            &! if .TRUE. no symmetry is used
    noinv = .FALSE.    ! if .TRUE. eliminates inversion symmetry
  !
  ! ... phonon calculation
  !
  INTEGER, PUBLIC :: &
    modenum            ! for single mode phonon calculation
  !
  ! ... printout control
  !
  LOGICAL, PUBLIC :: &
    reduce_io          ! if .TRUE. reduce the I/O to the strict minimum
  INTEGER, PUBLIC :: &
    iverbosity         ! type of printing ( 0 few, 1 all )
  !
  ! ... miscellany
  !
  LOGICAL, PUBLIC :: &
    use_para_diago = .FALSE.  ! if .TRUE. a parallel Householder algorithm can used
  INTEGER, PUBLIC :: &
   para_diago_dim = 0         ! minimum matrix dimension above which a parallel algorithm is used
  LOGICAL, PUBLIC :: &
    wg_set = .FALSE.   ! if .TRUE. the occupations are known
  !
  INTEGER  :: ortho_max = 0    ! maximum number of iterations in routine ortho
  REAL(DP) :: ortho_eps = 0.D0 ! threshold for convergence in routine ortho
  !
  ! ...  end of module-scope declarations
  ! 
  !=--------------------------------------------------------------------------=!
  CONTAINS
  !=--------------------------------------------------------------------------=!
    !
    !------------------------------------------------------------------------
    SUBROUTINE fix_dependencies()
      !------------------------------------------------------------------------
      ! 
      IMPLICIT NONE
      !
      ! ... Car Parrinello simulation
      !
      tcarpar = .TRUE.
      !
      IF ( tconjgrad .OR. t_diis .OR. tsteepdesc ) THEN
         !
         tcarpar = .FALSE.
         !
      END IF
      !
      ! ... if thdyn = .FALSE. set TSDC and TZEROC to .FALSE. too.
      !    
      IF ( .NOT. thdyn ) THEN
         !
         tsdc   = .FALSE.
         tzeroc = .FALSE.
         !
      END IF
      !
      IF ( .NOT. tfor ) THEN
         !
         tzerop = .FALSE.
         tv0rd  = .FALSE.
         tsdp   = .FALSE.
         tcp    = .FALSE.
         tcap   = .FALSE.
         tnosep = .FALSE.
         !
      ELSE
         !
         IF ( tsdp ) THEN
            !
            tcp    = .FALSE.
            tcap   = .FALSE.
            tnosep = .FALSE.
            tv0rd  = .FALSE.
            !
         END IF
         !
         IF ( tv0rd ) tzerop = .TRUE.
         !
      END IF
      !
      IF ( tsde ) tnosee = .FALSE.
      !
      CALL check_flags()
      !
      RETURN
      !
    END SUBROUTINE fix_dependencies
    !
    !------------------------------------------------------------------------
    SUBROUTINE check_flags()
      !------------------------------------------------------------------------
      !
      ! ...  do some checks for consistency
      !
      IF ( tnosee .AND. t_diis ) &
         CALL errore( ' control_flags ', 'DIIS + ELECT. NOSE ? ', 0 )
      !
      !IF ( tortho .AND. t_diis ) &
      !   CALL errore(' control_flags ','DIIS, ORTHO NOT PERMITTED',0)
      !
      IF ( tnosep .AND. tcp ) &
         CALL errore( ' control_flags ', ' TCP AND TNOSEP BOTH TRUE', 0 )
      !
      IF ( tnosep .AND. tcap ) &
         CALL errore( ' control_flags ', ' TCAP AND TNOSEP BOTH TRUE', 0 )
      !
      IF ( tcp .AND. tcap ) &
         CALL errore( ' control_flags ', ' TCP AND TCAP BOTH TRUE', 0 )
      !
      IF ( tdipole .AND. thdyn ) &
         CALL errore( '  control_flags  ', ' DIPOLE WITH CELL DYNAMICS ', 0 )
      !
      IF ( tv0rd .AND. tsdp ) &
         CALL errore( ' control_flags ', &
                    & ' READING IONS VELOCITY WITH STEEPEST D.', 0 )
      !
      RETURN
      !
    END SUBROUTINE check_flags
    !
END MODULE control_flags

