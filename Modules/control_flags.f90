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
        SAVE

        PRIVATE

        TYPE convergence_criteria
          LOGICAL :: active
          INTEGER :: nstep
          REAL(dbl) :: ekin
          REAL(dbl) :: derho
          REAL(dbl) :: force
        END TYPE


        TYPE ionic_conjugate_gradient
          LOGICAL :: active
          INTEGER :: nstepix
          INTEGER :: nstepex
          REAL(dbl) :: ionthr
          REAL(dbl) :: elethr
        END TYPE


        PUBLIC :: tbeg, nomore, &
                  nbeg, isave, iprint, tv0rd, nv0rd, tzeroc, tzerop, newnfi, tnewnfi, &
                  tfor, tpre, tzeroe, tsde, tsdp, tsdc, taurdr, ndr, &
                  ndw, tortho, tstress, tprnfor, prn, timing, &
                  memchk, tconjgrad, tprnsfac, toptical, &
                  tcarpar, rhoout, trane, ampre, tranp, amprp, &
                  tdipole, t_diis, t_diis_simple, t_diis_rot, tnosee, tnosep, tnoseh, &
                  tcp, tcap, tnodump, tdamp, tdampions, tconvthrs, &
                  convergence_criteria, tionstep, nstepe, tsteepdesc, &
                  ionic_conjugate_gradient, tconjgrad_ion, &
                  tatomicwfc, tscreen, gamma_only, ekin_conv_thr, ekin_maxiter, force_pairing


        PUBLIC ::  fix_dependencies, check_flags

        PUBLIC :: tbuff, tvlocw, trhor, trhow, thdyn, iprsta

        PUBLIC :: twfcollect
        


! ...   declare execution control variables

        LOGICAL :: tbuff     = .FALSE.  ! save wfc on unit 21  ( only cp, never used)
        LOGICAL :: tvlocw    = .FALSE.  ! write potential to unit 46 ( only cp, seldom used)
        LOGICAL :: trhor     = .FALSE.  ! read rho from unit 47  (only cp, seldom used)
        LOGICAL :: trhow     = .FALSE.  ! write rho to  unit 47 (only cp, seldom used)

        LOGICAL :: prn       = .FALSE.  ! verbosity
        LOGICAL :: tsde      = .FALSE.  ! electronic steepest descent
        LOGICAL :: tzeroe    = .FALSE.  ! set to zero the electronic velocities
        LOGICAL :: tfor      = .FALSE.  ! move the ions ( calculate forces )
        LOGICAL :: tsdp      = .FALSE.  ! ionic steepest descent
        LOGICAL :: tzerop    = .FALSE.  ! set to zero the ionic velocities
        LOGICAL :: tprnfor   = .FALSE.  ! print forces to standard output
        LOGICAL :: taurdr    = .FALSE.  ! read ionic position from standard input
        LOGICAL :: tv0rd     = .FALSE.  ! read ionic velocities from standard input
        LOGICAL :: tpre      = .FALSE.  ! calculate stress, and (in fpmd) variable cell dynamic
        LOGICAL :: thdyn     = .FALSE.  ! variable-cell dynamics (only cp)
        LOGICAL :: tsdc      = .FALSE.  ! cell geometry steepest descent
        LOGICAL :: tzeroc    = .FALSE.  ! set to zero the cell geometry velocities
        LOGICAL :: tstress   = .FALSE.  ! print stress to standard output
        LOGICAL :: tortho    = .FALSE.  ! use iterative orthogonalization 
        LOGICAL :: tconjgrad = .FALSE.  ! use conjugate gradient electronic minimization
        LOGICAL :: timing    = .FALSE.  ! print out timing information
        LOGICAL :: memchk    = .FALSE.  ! check for memory leakage
        LOGICAL :: tprnsfac  = .FALSE.  ! print out structure factor 
        LOGICAL :: toptical  = .FALSE.  ! print out optical properties
        LOGICAL :: tcarpar   = .FALSE.  ! tcarpar is set TRUE for a "pure" Car Parrinello simulation
        LOGICAL :: rhoout    = .FALSE.  ! print out charge densities
        LOGICAL :: tnodump   = .FALSE.  ! if true avoid the dump of the wavefunctions in resetart
        LOGICAL :: tdamp     = .FALSE.  ! Use damped dinamics for electrons
        LOGICAL :: tdampions = .FALSE.  ! Use damped dinamics for electrons
        LOGICAL :: tatomicwfc= .FALSE.  ! Use atomic wavefunctions as starting guess for ch. density
        LOGICAL :: tscreen   = .FALSE.  ! Use screened coulomb potentials for cluster calculations
        LOGICAL :: twfcollect = .FALSE. ! Collect wave function in the restart file at the end of run.

! ...   Force pairing
        LOGICAL :: force_pairing


        TYPE (convergence_criteria) :: tconvthrs
                              !  thresholds used to check GS convergence 

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

        LOGICAL :: gamma_only = .TRUE. !  true if only gamma point is used


        LOGICAL :: tnewnfi = .FALSE.
        INTEGER :: newnfi  = 0

! ...   Wave function randomization
        LOGICAL   :: trane = .FALSE.
        REAL(dbl) :: ampre = 0.0d0

! ...   Ionic position randomization
        LOGICAL   :: tranp(nsx) = .FALSE.
        REAL(dbl) :: amprp(nsx) = 0.0d0

! ...   Read the cell from standard input
        LOGICAL   :: tbeg = .FALSE.

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

        REAL(dbl) :: ekin_conv_thr = 0.0d0
        INTEGER   :: ekin_maxiter = 100
        REAL(dbl) :: etot_conv_thr = 0.0d0
        INTEGER   :: etot_maxiter = 100
        REAL(dbl) :: forc_conv_thr = 0.0d0
        INTEGER   :: forc_maxiter = 100

  !
  ! ... Several variables controlling the run
  !
  ! ... Control variables use mainly in PW calculations
  !
  !
  REAL(KIND=DP), PUBLIC  :: &
       mixing_beta,      &! the mixing parameter
       tr2,              &! the convergence threshold for potential
       upscale,          &! maximum reduction of convergence threshold
       ethr,             &! the convergence threshold for eigenvalues
       alpha0,           &! the mixing parameters for the extrapolation
       beta0,            &! of the starting potential
       diis_ethr_cg       ! threshold in eigval for starting DIIS

  INTEGER, PUBLIC :: &
       ngm0,             &! used in mix_rho
       niter,            &! the maximum number of iteration
       nmix,             &! the number of iteration kept in the history
       imix,             &! the type of mixing (0=plain,1=TF,2=local-TF)
       ! iprint,           &! the interval between full writing of results
       iverbosity,       &! type of printing ( 0 few, 1 all )
       david,            &! used on Davidson diagonalization
       nstep,            &! number of minimization steps
       istep,            &! current minimization step
       isolve,           &! Davidson or CG diagonalization
       iswitch,          &! general switch for the calculation type
       modenum,          &! used with iswitch=-4
       max_cg_iter,      &! maximum number of iterations in a CG di
       diis_buff,        &! dimension of the buffer in diis
       diis_ndim,        &! dimension of reduced basis in DIIS
       history,          &! number of old steps available for potential updating
       order              ! type of potential updating ( see update_pot )
  !
  LOGICAL, PUBLIC :: &
       lscf,             &! if .TRUE. the calculation is selfconsistent
       lbfgs,            &! if .TRUE. the calculation is a relaxation based on new BFGS scheme
       loldbfgs,         &! if .TRUE. the calculation is a bfgs-type relaxation based on the old scheme
       lmd,              &! if .TRUE. the calculation is a dynamics
       lneb,             &! if .TRUE. the calculation is neb
       lphonon,          &! if .TRUE. the calculation is phonon
       lraman,           &! if .TRUE. the calculation is raman
       lconstrain,       &! if .TRUE. the calculation is constraint
       ldamped,          &! if .TRUE. the calculation is a damped dynamics
       conv_elec,        &! if .TRUE. electron convergence has been reached
       conv_ions,        &! if .TRUE. ionic convergence has been reached
       nosym,            &! if .TRUE. no symmetry is used
       noinv = .FALSE.,  &! if .TRUE. eliminates inversion symmetry
       diis_wfc_keep,    &! if .TRUE. keeps old wfc for starting
       restart,          &! if .TRUE. restart from results of a preceding run
       reduce_io          ! if .TRUE. reduce the I/O to the strict minimum
  !


!  end of module-scope declarations
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
  CONTAINS
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
        IF( tconjgrad .OR. t_diis .OR. tsteepdesc ) THEN
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
