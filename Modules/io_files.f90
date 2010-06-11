!
! Copyright (C) 2002-2010 Quantum ESPRESSO group
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
  CHARACTER(len=256) :: tmp_dir = './'            ! directory for temporary files
  CHARACTER(len=256) :: wfc_dir = 'undefined'     ! directory for large files on each node, should be kept 'undefined' if not known 
  CHARACTER(len=256) :: prefix  = 'os'            ! prepended to file names
  CHARACTER(len=6)   :: nd_nmbr = '000000'        ! node number (used only in parallel case)
  CHARACTER(len=256) :: pseudo_dir = './'
  CHARACTER(len=256) :: psfile( ntypx ) = 'UPF'
  CHARACTER(len=256) :: outdir  = './'
  !
  CHARACTER(len=256) :: qexml_version = ' '       ! the format of the current qexml datafile 
  LOGICAL            :: qexml_version_init = .FALSE.  ! whether the fmt has been read or not
  !
  CHARACTER(LEN=256) :: input_drho = ' '          ! name of the file with the input drho
  CHARACTER(LEN=256) :: output_drho = ' '         ! name of the file with the output drho
  !
  CHARACTER(LEN=5 ), PARAMETER :: crash_file    = 'CRASH'
  CHARACTER (LEN=256) :: &
    dat_file      = 'os.dat',    &! file containing the enegy profile
    int_file      = 'os.int',    &! file containing the interpolated energy profile
    crd_file      = 'os.crd',    &! file containing path coordinates in pw.x input format
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
  INTEGER :: crashunit   = 15
  INTEGER :: pseudounit  = 10
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
  INTEGER :: iuntmp      = 90 ! temporary unit, when used must be closed ASAP
  !
  INTEGER :: nwordwfc    =  2 ! length of record in wavefunction file
  INTEGER :: nwordatwfc  =  2 ! length of record in atomic wfc file
  INTEGER :: nwordwann   =  2 ! length of record in sic wfc file
  !
  ! ... "path" specific
  !
  INTEGER :: iunpath     =  6 ! unit for string output ( stdout or what else )
  INTEGER :: iunrestart  = 2021 ! unit for saving the restart file ( neb_file )
  INTEGER :: iundat      = 2022 ! unit for saving the enegy profile
  INTEGER :: iunint      = 2023 ! unit for saving the interpolated energy profile
  INTEGER :: iunxyz      = 2024 ! unit for saving coordinates ( xyz format )
  INTEGER :: iunaxsf     = 2025 ! unit for saving coordinates ( axsf format )
  INTEGER :: iunbroy     = 2026 ! unit for saving broyden's history
  INTEGER :: iuncrd      = 2027 ! unit for saving coordinates in pw.x input format
  !
  !... finite electric field (Umari)
  !
  INTEGER :: iunefield   = 31 ! unit to store wavefunction for calculatin electric field operator
  !
  INTEGER :: iunefieldm  = 32 !unit to store projectors for hermitean electric field potential
  !
  INTEGER :: iunefieldp  = 33 !unit to store projectors for hermitean electric field potential
  !
  ! ... For Wannier Hamiltonian
  !
  INTEGER :: iunwpp   = 113
  INTEGER :: iunwf    = 114
  INTEGER :: nwordwpp = 2
  INTEGER :: nwordwf  = 2

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
  !--------------------------------------------------------------------------
  FUNCTION check_writable ( file_path, process_id ) RESULT ( ios )
    !--------------------------------------------------------------------------
    !
    ! ... if run by multiple processes, specific "process_id" to avoid
    ! ... opening, closing, deleting the same file from different processes
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),  INTENT(IN) :: file_path
    INTEGER, OPTIONAL, INTENT(IN) :: process_id
    !
    INTEGER :: ios
    !
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    ! ... check whether the scratch directory is writable
    ! ... note that file_path should end by a "/"
    !
    IF ( PRESENT (process_id ) ) THEN
       OPEN( UNIT = 4, FILE = TRIM(file_path) // 'test' // &
           & TRIM( int_to_char ( process_id ) ), &
           & STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
    ELSE
       OPEN( UNIT = 4, FILE = TRIM(file_path) // 'test', &
             STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
    END IF
    !
    CLOSE( UNIT = 4, STATUS = 'DELETE' )
    !
  END FUNCTION check_writable 
!=----------------------------------------------------------------------------=!
END MODULE io_files
!=----------------------------------------------------------------------------=!
