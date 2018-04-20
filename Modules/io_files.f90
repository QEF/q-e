!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
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
  USE kinds,      ONLY: dp
  USE io_global,  ONLY: ionode, ionode_id, stdout
  !
  ! ... I/O related variables: file names, units, utilities
  ! ... IMPORTANT: when directory names are set, they must always end with "/"
  !
  IMPLICIT NONE
  !
  SAVE
  PUBLIC :: create_directory, check_tempdir, clean_tempdir, check_file_exist, &
       delete_if_present, check_writable
  !
  ! ... directory for all temporary files
  CHARACTER(len=256) :: tmp_dir = './'
  ! ... directory for large files on each node. Default: same as tmp_dir
  CHARACTER(len=256) :: wfc_dir = 'undefined'
  ! ... prefix is prepended to all file (and directory) names 
  CHARACTER(len=256) :: prefix  = 'os'
  ! ... for parallel case and distributed I/O: node number
  CHARACTER(len=6)   :: nd_nmbr = '000000'
  ! ... directory where pseudopotential files are found
  CHARACTER(len=256) :: pseudo_dir = './'
  ! ... location of PP files after a restart from file
  CHARACTER(len=256) :: pseudo_dir_cur = ' '
  CHARACTER(len=256) :: psfile( ntypx ) = 'UPF'
  !
  CHARACTER(len=256) :: qexml_version = ' '       ! the format of the current qexml datafile 
  LOGICAL            :: qexml_version_init = .FALSE.  ! whether the fmt has been read or not
  !
!
  CHARACTER(LEN=256) :: qexsd_fmt = ' ', qexsd_version = ' '
  LOGICAL            :: qexsd_init = .FALSE. 
!
  CHARACTER(LEN=256) :: input_drho = ' '          ! name of the file with the input drho
  CHARACTER(LEN=256) :: output_drho = ' '         ! name of the file with the output drho
  !
  CHARACTER(LEN=5 ), PARAMETER :: crash_file  = 'CRASH'
  CHARACTER (LEN=261) :: exit_file = 'os.EXIT' ! file required for a soft exit  
  !
  CHARACTER (LEN=13), PARAMETER :: xmlpun      = 'data-file.xml'
  !
  CHARACTER (LEN=20), PARAMETER :: xmlpun_schema = 'data-file-schema.xml'
  !
  ! ... The units where various variables are saved
  ! ... Only units that are kept open during the run should be listed here
  !
  INTEGER :: iunres      =  1 ! unit for the restart of the run
  INTEGER :: iunpun      =  4 ! unit for saving the final results (data-file.xml)
  INTEGER :: iunwfc      = 10 ! unit with wavefunctions
  INTEGER :: iunoldwfc   = 11 ! unit with old wavefunctions
  INTEGER :: iunoldwfc2  = 12 ! as above at step -2
  INTEGER :: iunhub      = 13 ! unit for saving Hubbard-U atomic wfcs 
  INTEGER :: iunsat      = 14 ! unit for saving (orthogonal) atomic wfcs * S
  INTEGER :: iunmix      = 15 ! unit for saving mixing information
  INTEGER :: iunwfc_exx  = 16 ! unit with exx wavefunctions
  !
  INTEGER :: iunexit     = 26 ! unit for a soft exit  
  INTEGER :: iunupdate   = 27 ! unit for saving old positions (extrapolation)
  ! NEB
  INTEGER :: iunnewimage = 28 ! unit for parallelization among images
  INTEGER :: iunlock     = 29 ! as above (locking file)
  !
  INTEGER :: iunbfgs     = 30 ! unit for the bfgs restart file
  !
  INTEGER :: iuntmp      = 90 ! temporary unit, when used must be closed ASAP
  !
  INTEGER :: nwordwfc    =  2 ! length of record in wavefunction file
  INTEGER :: nwordatwfc  =  2 ! length of record in atomic wfc file
  INTEGER :: nwordwfcU   =  2 ! length of record in atomic hubbard wfc file
  INTEGER :: nwordwann   =  2 ! length of record in sic wfc file
  !
  !... finite electric field
  !
  INTEGER :: iunefield   = 31 ! unit to store wavefunction for calculating
                              ! electric field operator
  INTEGER :: iunefieldm  = 32 ! unit to store projectors for hermitean
                              ! electric field potential
  INTEGER :: iunefieldp  = 33 ! unit to store projectors for hermitean 
                              ! electric field potential
  !
  ! ... For Wannier Hamiltonian
  !
  INTEGER :: iunwpp   = 113
  INTEGER :: iunwf    = 114
  INTEGER :: nwordwpp = 2
  INTEGER :: nwordwf  = 2
  !
CONTAINS
  !
  !------------------------------------------------------------------------
  SUBROUTINE create_directory( dirname )
    !------------------------------------------------------------------------
    !
    USE wrappers,  ONLY : f_mkdir_safe
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_images, ONLY : me_image, intra_image_comm
    !
    CHARACTER(LEN=*), INTENT(IN) :: dirname
    !
    INTEGER                    :: ierr
    !
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    IF ( ionode ) ierr = f_mkdir_safe( TRIM( dirname ) )
    CALL mp_bcast ( ierr, ionode_id, intra_image_comm )
    !
    CALL errore( 'create_directory', &
         'unable to create directory ' // TRIM( dirname ), ierr )
    !
    ! ... syncronize all jobs (not sure it is really useful)
    !
    CALL mp_barrier( intra_image_comm )
    !
    ! ... check whether the scratch directory is writable
    !
    IF ( ionode ) ierr = check_writable ( dirname, me_image )
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    !
    CALL errore( 'create_directory:', &
         TRIM( dirname ) // ' non existent or non writable', ierr )
    !
    RETURN
    !
  END SUBROUTINE create_directory
  !
  !-----------------------------------------------------------------------
  SUBROUTINE check_tempdir ( tmp_dir, exst, pfs )
    !-----------------------------------------------------------------------
    !
    ! ... Verify if tmp_dir exists, creates it if not
    ! ... On output:
    ! ...    exst= .t. if tmp_dir exists
    ! ...    pfs = .t. if tmp_dir visible from all procs of an image
    !
    USE wrappers,      ONLY : f_mkdir_safe
    USE mp_images,     ONLY : intra_image_comm, nproc_image, me_image
    USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), INTENT(in) :: tmp_dir
    LOGICAL, INTENT(out)         :: exst, pfs
    !
    INTEGER             :: ios, image, proc, nofi
    CHARACTER (len=256) :: file_path, filename
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !
    ! ... create tmp_dir on ionode
    ! ... f_mkdir_safe returns -1 if tmp_dir already exists
    ! ...                       0 if         created
    ! ...                       1 if         cannot be created
    !
    IF ( ionode ) ios = f_mkdir_safe( TRIM(tmp_dir) )
    CALL mp_bcast ( ios, ionode_id, intra_image_comm )
    exst = ( ios == -1 )
    IF ( ios > 0 ) CALL errore ('check_tempdir','tmp_dir cannot be opened',1)
    !
    ! ... let us check now if tmp_dir is visible on all nodes
    ! ... if not, a local tmp_dir is created on each node
    !
    ios = f_mkdir_safe( TRIM(tmp_dir) )
    CALL mp_sum ( ios, intra_image_comm )
    pfs = ( ios == -nproc_image ) ! actually this is true only if .not.exst 
    !
    RETURN
    !
  END SUBROUTINE check_tempdir
  !
  !-----------------------------------------------------------------------
  SUBROUTINE clean_tempdir( tmp_dir )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), INTENT(in) :: tmp_dir
    !
    CHARACTER (len=256) :: file_path, filename
    !
    ! ... remove temporary files from tmp_dir ( only by the master node )
    !
    file_path = trim( tmp_dir ) // trim( prefix )
    IF ( ionode ) THEN
       CALL delete_if_present( trim( file_path ) // '.update' )
       CALL delete_if_present( trim( file_path ) // '.md' )
       CALL delete_if_present( trim( file_path ) // '.bfgs' )
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE clean_tempdir
  !
  !------------------------------------------------------------------------
  FUNCTION check_file_exist( filename )
    !------------------------------------------------------------------------
    !
    USE mp,        ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    !
    IMPLICIT NONE
    !
    LOGICAL          :: check_file_exist
    CHARACTER(LEN=*) :: filename
    !
    LOGICAL :: lexists
    !
    IF ( ionode ) THEN 
       !
       INQUIRE( FILE = TRIM( filename ), EXIST = lexists )
       !
    ENDIF
    !
    CALL mp_bcast ( lexists, ionode_id, intra_image_comm )
    !
    check_file_exist = lexists
    RETURN
    !
  END FUNCTION check_file_exist
  !
  !--------------------------------------------------------------------------
  SUBROUTINE delete_if_present( filename, in_warning )
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),  INTENT(IN) :: filename
    LOGICAL, OPTIONAL, INTENT(IN) :: in_warning
    LOGICAL                       :: exst, warning
    INTEGER                       :: iunit
    INTEGER, EXTERNAL :: find_free_unit
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
    !-----------------------------------------------------------------------
  END FUNCTION check_writable 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine diropn (unit, extension, recl, exst, tmp_dir_)
  !-----------------------------------------------------------------------
  !
  !     Opens a direct-access file named "prefix"."extension" in directory
  !     "tmp_dir_" if specified, in "tmp_dir" otherwise. 
  !     In parallel execution, the node number is added to the file name.
  !     The record length is "recl" double-precision numbers.
  !     On output, "exst" is .T. if opened file already exists
  !     If recl=-1, the file existence is checked, nothing else is done
  !
  implicit none
  !
  !    first the input variables
  !
  character(len=*) :: extension
  ! input: name of the file to open
  character(len=*), optional :: tmp_dir_
  ! optional variable, if present it is used as tmp_dir
  integer :: unit, recl
  ! input: unit of the file to open
  ! input: length of the records
  logical :: exst
  ! output: if true the file exists
  !
  !    local variables
  !
  character(len=256) :: tempfile, filename
  ! complete file name
  real(dp):: dummy
  integer*8 :: unf_recl
  ! double precision to prevent integer overflow
  integer :: ios, direct_io_factor
  logical :: opnd
  !
  !    initial checks
  !
  if (unit < 0) call errore ('diropn', 'wrong unit', 1)
  !
  !    ifirst we check that the file is not already openend
  !
  ios = 0
  inquire (unit = unit, opened = opnd)
  if (opnd) call errore ('diropn', "can't open a connected unit", abs(unit))
  !
  !    then we check the filename extension
  !
  if (extension == ' ') call errore ('diropn','filename extension not given',2)
  filename = trim(prefix) // "." // trim(extension)
  if (present(tmp_dir_)) then
     tempfile = trim(tmp_dir_) // trim(filename) //nd_nmbr
  else
     tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
  endif

  inquire (file = tempfile, exist = exst)
  if (recl == -1) RETURN
  !
  ! the  record length in direct-access I/O is given by the number of
  ! real*8 words times direct_io_factor (may depend on the compiler)
  !
  INQUIRE (IOLENGTH=direct_io_factor) dummy
  unf_recl = direct_io_factor * int(recl, kind=kind(unf_recl))
  if (unf_recl <= 0) call errore ('diropn', 'wrong record length', 3)
  !
  open (unit, file=trim(adjustl(tempfile)), iostat=ios, form='unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)

  if (ios /= 0) call errore ('diropn', 'error opening '//trim(tempfile), unit)
  return
  !-----------------------------------------------------------------------
end subroutine diropn
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine seqopn (unit, extension, formatt, exst, tmp_dir_)
  !-----------------------------------------------------------------------
  !
  !     this routine opens a file named "prefix"."extension"
  !     in tmp_dir for sequential I/O access
  !     If appropriate, the node number is added to the file name
  !
  implicit none
  !
  !    first the dummy variables
  !
  character(len=*) :: formatt, extension
  ! input: name of the file to connect
  ! input: 'formatted' or 'unformatted'
  character(len=*), optional :: tmp_dir_
  ! optional variable, if present it is used as tmp_dir
  integer :: unit
  ! input: unit to connect
  logical :: exst
  ! output: true if the file already exist
  !
  !    here the local variables
  !
  character(len=256) :: tempfile, filename
  ! complete file name
  integer :: ios
  ! integer variable to test I/O status
  logical :: opnd
  ! true if the file is already opened


  if (unit < 1) call errore ('seqopn', 'wrong unit', 1)
  !
  !    test if the file is already opened
  !
  ios = 0
  inquire (unit = unit, opened = opnd)
  if (opnd) call errore ('seqopn', "can't open a connected unit", &
       abs (unit) )
  !
  !      then we check the extension of the filename
  !
  if (extension.eq.' ') call errore ('seqopn','filename extension  not given',2)
  filename = trim(prefix) // "." // trim(extension)
  ! Use the tmp_dir from input, if available
  if ( present(tmp_dir_) ) then
    tempfile = trim(tmp_dir_) // trim(filename)
  else
    tempfile = trim(tmp_dir) // trim(filename)
  end if
  if ( trim(nd_nmbr) /= '1'     .and. trim(nd_nmbr) /= '01'   .and. &
       trim(nd_nmbr) /= '001'   .and. trim(nd_nmbr) /= '0001' .and. &
       trim(nd_nmbr) /= '00001' .and. trim(nd_nmbr) /= '000001' ) then
     !
     ! do not add processor number to files opened by processor 1
     ! in parallel execution: if only the first processor writes,
     ! we do not want the filename to be dependent on the number
     ! of processors
     !
     tempfile = trim(tempfile) // nd_nmbr
  end if
  inquire (file = tempfile, exist = exst)
  !
  !    Open the file
  !
  open (unit = unit, file = tempfile, form = formatt, status = &
       'unknown', iostat = ios)

  if (ios /= 0) call errore ('seqopn', 'error opening '//trim(tempfile), unit)
  return
  !-----------------------------------------------------------------------
end subroutine seqopn
!-----------------------------------------------------------------------
!
!=----------------------------------------------------------------------------=!
END MODULE io_files
!=----------------------------------------------------------------------------=!
!
!----------------------------------------------------------------------------
SUBROUTINE davcio( vect, nword, unit, nrec, io )
  !----------------------------------------------------------------------------
  !
  ! ... direct-access vector input/output
  ! ... read/write nword words starting from the address specified by vect
  !
  USE kinds ,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nword, unit, nrec, io
    ! input: the dimension of vect
    ! input: the unit where to read/write
    ! input: the record where to read/write
    ! input: flag if < 0 reading if > 0 writing
  REAL(DP), INTENT(INOUT) :: vect(nword)
   ! input/output: the vector to read/write
  !
  INTEGER :: ios
    ! integer variable for I/O control
  LOGICAL :: opnd
  CHARACTER*256 :: name
  !
  !
  CALL start_clock( 'davcio' )
  !
  IF ( unit  <= 0 ) CALL errore(  'davcio', 'wrong unit', 1 )
  IF ( nrec  <= 0 ) CALL errore(  'davcio', 'wrong record number', 2 )
  IF ( nword <= 0 ) CALL errore(  'davcio', 'wrong record length', 3 )
  IF ( io    == 0 ) CALL infomsg( 'davcio', 'nothing to do?' )
  !
  INQUIRE( UNIT = unit, OPENED = opnd, NAME = name )
  !
  IF ( .NOT. opnd ) &
     CALL errore(  'davcio', 'unit is not opened', unit )
  !
  ios = 0
  !
  IF ( io < 0 ) THEN
     !
     READ( UNIT = unit, REC = nrec, IOSTAT = ios ) vect
     IF ( ios /= 0 ) CALL errore( 'davcio', &
         & 'error while reading from file "' // TRIM(name) // '"', unit )
     !
  ELSE IF ( io > 0 ) THEN
     !
     WRITE( UNIT = unit, REC = nrec, IOSTAT = ios ) vect
     IF ( ios /= 0 ) CALL errore( 'davcio', &
         & 'error while writing from file "' // TRIM(name) // '"', unit )
     !
  END IF
  !
  CALL stop_clock( 'davcio' )
  !
  RETURN
  !
END SUBROUTINE davcio

