!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=----------------------------------------------------------------------------=!
!
MODULE path_input_parameters_module
!
!=----------------------------------------------------------------------------=!
!
!  this module contains
!  1) the definition of all input parameters for NEB
!  2) the definition of namelist PATH
!  3) routines that allocate/deallocate data needed in input
!  Based upon original NEB implementation ( C.S. 17/10/2003 )
!
!=----------------------------------------------------------------------------=!
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
!=----------------------------------------------------------------------------=!
! BEGIN manual
!
!
! * DESCRIPTION OF THE INPUT FILE
!  (to be given as standard input)
!
!  The input file has the following layout:
!
!     &PATH
!       path_parameter_1,
!       path_parameter_2,
!       .......
!       path_parameter_Lastone
!     /
!     ATOMIC_SPECIES
!      slabel_1 mass_1 pseudo_file_1
!      slabel_2 mass_2 pseudo_file_2
!      .....
!     PATH_ATOMIC_POSITIONS
!      alabel_1  px_1 py_1 pz_1
!      alabel_2  px_2 py_2 pz_2
!      .....
!     CARD_3
!     ....
!     CARD_N
!
!  -- end of input file --
!
! ... variables added for "path" calculations
!
  INTEGER :: nat = 1
  REAL(DP) :: alat
  !
  CHARACTER(len=80) :: restart_mode
  ! specify how to start/restart the simulation
  CHARACTER(len=80) :: restart_mode_allowed(3)
  DATA restart_mode_allowed / 'from_scratch', 'restart', 'reset_counters' /
  !
  INTEGER :: nstep_path
  !
  CHARACTER(len=80) :: string_method = 'neb' 
  ! 'neb' traditional neb as described by Jonsson
  !  'sm' strings method
  CHARACTER(len=80) :: string_method_scheme_allowed(2)
  DATA string_method_scheme_allowed / 'neb', 'sm' /
  !
  INTEGER :: input_images = 0
  !
  INTEGER :: num_of_images = 0
  !
  CHARACTER(len=80) :: CI_scheme = 'no-CI'
  ! CI_scheme = 'no-CI' | 'auto' | 'manual'
  ! set the Climbing Image scheme
  ! 'no-CI'       Climbing Image is not used
  ! 'auto'        Standard Climbing Image
  ! 'manual'      the image is selected by hand
  !
  CHARACTER(len=80) :: CI_scheme_allowed(3)
  DATA CI_scheme_allowed / 'no-CI', 'auto', 'manual' /
  !
  LOGICAL :: first_last_opt = .false.
  LOGICAL :: minimum_image  = .false.
  LOGICAL :: use_masses     = .false.
  LOGICAL :: use_freezing   = .false.
  LOGICAL :: fixed_tan      = .false.
  !
  CHARACTER(len=80) :: opt_scheme = 'quick-min'
  ! minimization_scheme = 'quick-min' | 'damped-dyn' |
  !                       'mol-dyn'   | 'sd'
  ! set the minimization algorithm
  ! 'quick-min'   projected molecular dynamics
  ! 'sd'          steepest descent
  ! 'broyden'     broyden acceleration
  ! 'broyden2'    broyden acceleration - better ?
  ! 'langevin'    langevin dynamics
  !
  CHARACTER(len=80) :: opt_scheme_allowed(5)
  DATA opt_scheme_allowed / 'quick-min', 'broyden', 'broyden2', 'sd', 'langevin' /
  !
  REAL (DP)  :: temp_req = 0.0_DP
  ! meaningful only when minimization_scheme = 'sim-annealing'
  REAL (DP)  :: ds = 1.0_DP
  !
  REAL (DP)  :: k_max = 0.1_DP, k_min = 0.1_DP
  !
  REAL (DP)  :: path_thr = 0.05_DP
  !
  LOGICAL    :: lfcp    = .FALSE.
  REAL(DP)   :: fcp_mu  = 0.0_DP  ! in eV
  REAL(DP)   :: fcp_thr = 0.01_DP ! in eV
  !
  CHARACTER(len=80) :: fcp_scheme = 'lm'
  ! fcp_scheme = 'lm' | 'newton' | 'coupled'
  ! set the minimization algorithm for FCP
  ! 'lm'          line-minimization
  ! 'newton'      newton-raphson
  ! 'coupled'     coupled-method with ionic postions
  !
  CHARACTER(len=80) :: fcp_scheme_allowed(8)
  DATA fcp_scheme_allowed / 'lm', 'line-min', 'line-minimization', 'line-minimisation', &
                          & 'newton', 'couple', 'coupled', 'coupling' /
  !
  INTEGER    :: fcp_ndiis    = 4
  REAL(DP)   :: fcp_rdiis    = 1.0_DP
  REAL(DP)   :: fcp_max_volt = 1.0_DP ! in eV
  !
  !
  NAMELIST / PATH / &
                    restart_mode,                                     &
                    string_method, nstep_path, num_of_images,         &
                    CI_scheme, opt_scheme, use_masses,                &
                    first_last_opt, ds, k_max, k_min, temp_req,       &
                    path_thr, fixed_tan, use_freezing, minimum_image, &
                    lfcp, fcp_mu, fcp_thr, fcp_scheme,                &
                    fcp_ndiis, fcp_rdiis, fcp_max_volt
!
!    ATOMIC_POSITIONS
!
        REAL(DP), ALLOCATABLE :: pos(:,:)
        INTEGER, ALLOCATABLE :: typ(:)
        !
!
!   TOTAL_CHARGE
!
      REAL(DP), ALLOCATABLE :: tot_charge(:)
      !
!
!   CLIMBING_IMAGES
!
      LOGICAL, ALLOCATABLE :: climbing( : )
! ----------------------------------------------------------------------

CONTAINS

  SUBROUTINE allocate_path_input_ions( num_of_images )
    !
    INTEGER, INTENT(in) :: num_of_images
    !
    IF ( allocated( pos ) ) DEALLOCATE( pos )
    IF ( allocated( typ ) ) DEALLOCATE( typ )
    !
    IF ( allocated( tot_charge ) ) DEALLOCATE( tot_charge )
    !
    ALLOCATE( pos( 3*nat, num_of_images ) )
    ALLOCATE( typ( nat ) )
    !
    ALLOCATE( tot_charge( num_of_images ) )
    !
    pos(:,:) = 0.0_DP
    !
    tot_charge = 0.0_DP
    !
    RETURN
    !
  END SUBROUTINE allocate_path_input_ions
  !
  SUBROUTINE deallocate_path_input_ions()
    !
    IF ( allocated( pos ) ) DEALLOCATE( pos )
    IF ( allocated( typ ) ) DEALLOCATE( typ )
    !
    IF ( allocated( tot_charge ) ) DEALLOCATE( tot_charge )
    !
    IF ( allocated( climbing ) ) DEALLOCATE( climbing )
    !
    RETURN
    !
  END SUBROUTINE deallocate_path_input_ions
  !
!=----------------------------------------------------------------------------=!
!
END MODULE path_input_parameters_module
!
!=----------------------------------------------------------------------------=!
