!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE io_files
!=----------------------------------------------------------------------------=!
  !
  USE parameters, ONLY: ntypx
  !
  ! ... The name of the files
  !
  IMPLICIT NONE
  !
  SAVE
  !
  CHARACTER(len=256) :: tmp_dir = './'  ! directory for temporary files
  CHARACTER(len=256) :: wfc_dir = 'undefined'  ! directory for large files on each node, should be kept 'undefined' if not known 
  CHARACTER(len=256) :: prefix  = 'os'     ! prepended to file names
  CHARACTER(len=6)   :: nd_nmbr = '000000' ! node number (used only in parallel case)
  CHARACTER(len=256) :: pseudo_dir = './'
  CHARACTER(len=256) :: psfile( ntypx ) = 'UPF'
  CHARACTER(len=256) :: scradir = './'
  CHARACTER(len=256) :: outdir  = './'
  !
  CHARACTER(LEN=256) :: input_drho = ' '   ! name of the file with the input drho
  CHARACTER(LEN=256) :: output_drho = ' '  ! name of the file with the output drho
  !
  CHARACTER(LEN=19) :: band_file = ' '
  CHARACTER(LEN=19) :: tran_file = ' '
  CHARACTER(LEN=256) :: prefixt   = ' '
  CHARACTER(LEN=256) :: prefixl   = ' '
  CHARACTER(LEN=256) :: prefixs   = ' '
  CHARACTER(LEN=256) :: prefixr   = ' '
  CHARACTER(LEN=256) :: save_file = ' '
  CHARACTER(LEN=256) :: fil_loc = ' '      !  file with 2D eigenvectors and eigenvalues
  !
  CHARACTER(LEN=14), PARAMETER :: rho_name      = 'CHARGE_DENSITY'
  CHARACTER(LEN=17), PARAMETER :: rho_name_up   = 'CHARGE_DENSITY.UP'
  CHARACTER(LEN=19), PARAMETER :: rho_name_down = 'CHARGE_DENSITY.DOWN'
  CHARACTER(LEN=14), PARAMETER :: rho_name_avg  = 'CHARGE_AVERAGE'
  !
  CHARACTER(LEN=4 ), PARAMETER :: chifile       = 'CHI2'
  CHARACTER(LEN=7 ), PARAMETER :: dielecfile    = 'EPSILON'
  !
  CHARACTER(LEN=15), PARAMETER :: empty_file    = 'EMPTY_STATES.WF'
  CHARACTER(LEN=5 ), PARAMETER :: crash_file    = 'CRASH'
  CHARACTER(LEN=7 ), PARAMETER :: stop_file     = '.cpstop'
  CHARACTER(LEN=2 ), PARAMETER :: ks_file       = 'KS'
  CHARACTER(LEN=6 ), PARAMETER :: ks_emp_file   = 'KS_EMP'
  CHARACTER(LEN=16), PARAMETER :: sfac_file     = 'STRUCTURE_FACTOR'
  CHARACTER (LEN=256) :: &
    dat_file      = 'os.dat',    &! file containing the enegy profile
    int_file      = 'os.int',    &! file containing the interpolated energy profile
    path_file     = 'os.path',   &! file containing informations needed to restart a path simulation
    xyz_file      = 'os.xyz',    &! file containing coordinates of all images in xyz format
    axsf_file     = 'os.axsf',   &! file containing coordinates of all images in axsf format
    broy_file     = 'os.broyden'  ! file containing broyden's history
  CHARACTER (LEN=261) :: &
    exit_file = "os.EXIT"    ! file required for a soft exit  
  !
  CHARACTER (LEN=9),  PARAMETER :: xmlpun_base = 'data-file'
  CHARACTER (LEN=13), PARAMETER :: xmlpun      = xmlpun_base // '.xml'
  !
  ! ... The units where various variables are saved
  !
  INTEGER :: rhounit     = 17
  INTEGER :: emptyunit   = 19
  INTEGER :: crashunit   = 15
  INTEGER :: stopunit    = 7
  INTEGER :: ksunit      = 18
  INTEGER :: sfacunit    = 20
  INTEGER :: pseudounit  = 10
  INTEGER :: chiunit     = 20
  INTEGER :: dielecunit  = 20
  INTEGER :: opt_unit    = 20 ! optional unit 
  !
  ! ... units in pwscf
  !
  INTEGER :: iunres      =  1 ! unit for the restart of the run
  INTEGER :: iunpun      =  4 ! unit for saving the final results
  INTEGER :: iunwfc      = 10 ! unit with wavefunctions
  INTEGER :: iunoldwfc   = 11 ! unit with old wavefunctions
  INTEGER :: iunoldwfc2  = 12 ! as above at step -2
  INTEGER :: iunat       = 13 ! unit for saving (orthogonal) atomic wfcs 
  INTEGER :: iunsat      = 14 ! unit for saving (orthogonal) atomic wfcs * S
  INTEGER :: iunocc      = 15 ! unit for saving the atomic n_{ij}
  INTEGER :: iunigk      = 16 ! unit for saving indices
  INTEGER :: iunpaw      = 17 ! unit for saving paw becsum and D_Hxc
  !
  INTEGER :: iunexit     = 26 ! unit for a soft exit  
  INTEGER :: iunupdate   = 27 ! unit for saving old positions (extrapolation)
  INTEGER :: iunnewimage = 28 ! unit for parallelization among images
  INTEGER :: iunlock     = 29 ! as above (locking file)
  !
  INTEGER :: iunbfgs     = 30 ! unit for the bfgs restart file
  INTEGER :: iunatsicwfc = 31 ! unit for sic wfc
  !
  INTEGER :: nwordwfc    =  2 ! lenght of record in wavefunction file
  INTEGER :: nwordatwfc  =  2 ! lenght of record in atomic wfc file
  INTEGER :: nwordwann   =  2 ! lenght of record in sic wfc file
  !
  ! ... "path" specific
  !
  INTEGER :: iunpath     =  6 ! unit for string output ( stdout or what else )
  INTEGER :: iunrestart  = 21 ! unit for saving the restart file ( neb_file )
  INTEGER :: iundat      = 22 ! unit for saving the enegy profile
  INTEGER :: iunint      = 23 ! unit for saving the interpolated energy profile
  INTEGER :: iunxyz      = 24 ! unit for saving coordinates ( xyz format )
  INTEGER :: iunaxsf     = 25 ! unit for saving coordinates ( axsf format )
  INTEGER :: iunbroy     = 26 ! unit for saving broyden's history
  !
  ! ... meta-dynamics
  !
  INTEGER :: iunmeta     = 77 ! unit for saving meta-dynamics history
  !
  ! ... Y. Kanai combined smd/cp method
  !
  INTEGER :: smwout      = 20 ! base value to compute index for replica files
  !
  INTEGER :: vib_out     = 20 ! output of phrozen phonon vibrational calculation
  INTEGER :: vib_mass    = 21 ! isotope masses used for the dynamical matrix
  !
  !... finite electric field (Umari)
  !
  INTEGER :: iunefield   = 31 ! unit to store wavefunction for calculatin electric field operator
  !
  INTEGER :: iunefieldm  = 32 !unit to store projectors for hermitean electric field potential
  !
  INTEGER :: iunefieldp  = 33 !unit to store projectors for hermitean electric field potential
CONTAINS
  !
  !-----------------------------------------------------------------------
  FUNCTION trimcheck ( directory )
    !-----------------------------------------------------------------------
    !
    ! ... verify if directory ends with /, add one if needed; 
    ! ... trim white spaces and put the result in trimcheck
    !
    IMPLICIT NONE
    !
    CHARACTER (LEN=*), INTENT(IN) :: directory
    CHARACTER (LEN=256) :: trimcheck
    INTEGER  :: l
    !
    l = LEN_TRIM( directory )
    IF ( l == 0 ) CALL errore( 'trimcheck', ' input name empty', 1)
    !
    IF ( directory(l:l) == '/' ) THEN
       trimcheck = TRIM ( directory)
    ELSE
       IF ( l < LEN( trimcheck ) ) THEN
          trimcheck = TRIM ( directory ) // '/'
       ELSE
          CALL errore(  'trimcheck', ' input name too long', l )
       END IF
    END IF
    !
    RETURN
    !
  END FUNCTION trimcheck
  !
  !--------------------------------------------------------------------------
  FUNCTION find_free_unit()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: find_free_unit
    INTEGER :: iunit
    LOGICAL :: opnd
    !
    !
    unit_loop: DO iunit = 99, 1, -1
       !
       INQUIRE( UNIT = iunit, OPENED = opnd )
       !
       IF ( .NOT. opnd ) THEN
          !
          find_free_unit = iunit
          !
          RETURN
          !
       END IF
       !
    END DO unit_loop
    !
    CALL errore( 'find_free_unit()', 'free unit not found ?!?', 1 )
    !
    RETURN
    !
  END FUNCTION find_free_unit
  !
  !--------------------------------------------------------------------------
  SUBROUTINE delete_if_present( filename, in_warning )
    !--------------------------------------------------------------------------
    !
    USE io_global, ONLY : ionode, stdout
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),  INTENT(IN) :: filename
    LOGICAL, OPTIONAL, INTENT(IN) :: in_warning
    LOGICAL                       :: exst, warning
    INTEGER                       :: iunit
    !
    IF ( .NOT. ionode ) RETURN
    !
    INQUIRE( FILE = filename, EXIST = exst )
    !
    IF ( exst ) THEN
       !
       iunit = find_free_unit()
       !
       warning = .FALSE.
       !
       IF ( PRESENT( in_warning ) ) warning = in_warning
       !
       OPEN(  UNIT = iunit, FILE = filename , STATUS = 'OLD' )
       CLOSE( UNIT = iunit, STATUS = 'DELETE' )
       !
       IF ( warning ) &
          WRITE( UNIT = stdout, FMT = '(/,5X,"WARNING: ",A, &
               & " file was present; old file deleted")' ) filename
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE delete_if_present
  !
!=----------------------------------------------------------------------------=!
END MODULE io_files
!=----------------------------------------------------------------------------=!
