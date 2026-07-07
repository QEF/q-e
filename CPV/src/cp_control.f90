!
! Copyright (C) 2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE cp_control
  !=--------------------------------------------------------------------------=!
  !! This module contains flags specific for the CP code
  !----------------------------------------------
  !
  USE kinds
  USE parameters, ONLY : nsx
  !
  IMPLICIT NONE
  !
  SAVE
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
  ! ... thresholds used to check GS convergence
  !
  TYPE (convergence_criteria) :: tconvthrs
  !
  ! ...   declare execution control variables
  !
  LOGICAL :: trhor     = .FALSE. ! read rho from unit 47 (only cp, seldom used)
  LOGICAL :: trhow     = .FALSE. ! code, write rho to restart dir
  LOGICAL :: tksw      = .FALSE. ! write Kohn-Sham states to restart dir
  LOGICAL :: tfirst    = .TRUE.  ! true if first iteration after restart
  LOGICAL :: tlast     = .FALSE. ! true if last iteration before ending
  LOGICAL :: tprint    = .FALSE. ! set to true when calculation of time
                                 ! derivatives of wave functions must be 
                                 ! computed via projection on occupied manifold 
  LOGICAL :: tbeg      = .FALSE. ! Read the cell from standard input
  LOGICAL :: tcap      = .FALSE. ! Unknown, please document
  LOGICAL :: tcp       = .FALSE. ! Unknown, please document
  LOGICAL :: taurdr    = .FALSE. ! read ionic position from standard input
  LOGICAL :: trescalee = .FALSE. ! rescale the electronics velocities
  LOGICAL :: tprnfor   = .FALSE. ! print forces
  LOGICAL :: tfor      = .FALSE. ! move the ions ( calculate forces )
  LOGICAL :: thdyn     = .FALSE. ! variable-cell dynamics (only cp)
  LOGICAL :: tpre      = .FALSE. ! calculate stress
  !
  LOGICAL :: tsde      = .FALSE. ! electronic steepest descent
  LOGICAL :: tzeroe    = .FALSE. ! set to zero the electronic velocities
  LOGICAL :: tsdp      = .FALSE. ! ionic steepest descent
  LOGICAL :: tzerop    = .FALSE. ! set to zero the ionic velocities
  LOGICAL :: tsdc      = .FALSE. ! cell geometry steepest descent
  LOGICAL :: tzeroc    = .FALSE. ! set to zero the cell geometry velocities
  LOGICAL :: tnosee    = .FALSE. ! Nosé for electrons
  LOGICAL :: force_pairing = .FALSE. ! Force pairing
  ! next variable never used
  ! LOGICAL :: tscreen       = .FALSE. ! Use screened coulomb potentials for cluster calculations
  INTEGER :: nomore    = 0       ! Unknown,please document
  INTEGER :: nbeg      = 0       ! internal code for initialization (-1,0,1,2,... )
  INTEGER :: ndw       = 0       ! output unit
  INTEGER :: ndr       = 0       ! input unit
  INTEGER :: isave     = 0       ! write restart to ndw unit every isave step
  !
  ! ... Ionic vs Electronic step frequency
  ! ... When "ion_nstep > 1" and "electron_dynamics = 'md' | 'sd' ", ions are
  ! ... propagated every "ion_nstep" electronic step only if the electronic
  ! ... "ekin" is lower than "ekin_conv_thr"
  !
  LOGICAL :: tionstep = .FALSE.
  INTEGER :: nstepe   = 1
  !
  ! ... Wave function randomization
  !
  LOGICAL  :: trane = .FALSE.
  REAL(DP) :: ampre = 0.0_DP
  !
  ! ... Ionic position randomization
  !
  LOGICAL  :: tranp(nsx) = .FALSE.
  REAL(DP) :: amprp(nsx) = 0.0_DP
  !
  ! ... Variables used whenever a timestep change is requested
  ! ... dt_xlm_old is necessary to mantain compatibility with the old way 
  ! ... of changing the molecular dynamics integration timestep.
  ! ... The code needs to check, in case the old method is used, that 
  ! ... the input old timestep and the xml old timestep are the same
  !
  REAL(DP) :: dt_old = -1.0_DP
  REAL(DP) :: dt_xml_old = -1.0_DP 
  !
  ! exx_wf related 
  LOGICAL :: lwf         =.FALSE. ! if .TRUE. the calc. is with wannier functions
  LOGICAL :: lwfnscf     = .FALSE.
  LOGICAL :: lwfpbe0nscf = .FALSE.
  !
  ! ... iterative orthonormalization
  !
  LOGICAL :: tortho    = .FALSE. ! use iterative orthogonalization
  INTEGER  :: ortho_max = 0      ! maximum number of iterations in routine ortho
  REAL(DP) :: ortho_eps = 0.0_DP ! threshold for convergence in routine ortho
  !
  ! ... Number of neighbouring cell to consider in ewald sum
  !
  INTEGER :: iesr = 1
  !
END MODULE cp_control

