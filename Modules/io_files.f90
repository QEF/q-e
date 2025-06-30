!
! Copyright (C) 2002-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE io_files
!=----------------------------------------------------------------------------=!
  !! I/O related variables: file names, units, utilities.
  !
  ! ... IMPORTANT: when directory names are set, they must always end with "/"
  !
  USE parameters, ONLY: ntypx, nsolx
  USE kinds,      ONLY: dp
  USE io_global,  ONLY: ionode, ionode_id, stdout
  USE mp,         ONLY: mp_barrier, mp_bcast, mp_sum
  USE mp_images,  ONLY: me_image, intra_image_comm, nproc_image
  !
  IMPLICIT NONE
  !
  SAVE
  PUBLIC :: create_directory, check_tempdir, clean_tempdir, check_file_exist, &
       delete_if_present, check_writable, restart_dir, xmlfile
  !
  CHARACTER(len=256) :: tmp_dir = './'
  !! directory for all temporary files
  CHARACTER(len=256) :: wfc_dir = 'undefined'
  !! directory for large files on each node. Default: same as tmp_dir
  CHARACTER(len=256) :: prefix  = 'os'
  !! prefix is prepended to all file (and directory) names 
  !
#if defined (_WIN32)
#if defined (__PGI)
  CHARACTER(len=6) :: postfix  = '.save/'
  !! postfix is appended to directory names
#else
  CHARACTER(len=6) :: postfix  = '.save\'
#endif
#else
  CHARACTER(len=6) :: postfix  = '.save/'
#endif
  !
  CHARACTER(len=6)   :: nd_nmbr = '000000'
  !! for parallel case and distributed I/O: node number
  CHARACTER(len=256) :: pseudo_dir = './'
  !! directory where pseudopotential files are found
  CHARACTER(len=256) :: pseudo_dir_cur = ' '
  !! location of PP files after a restart from file
  CHARACTER(len=256) :: psfile( ntypx ) = 'UPF'
  !! default: UPF
  CHARACTER(len=256) :: molfile( nsolx ) = 'MOL'
  !
  CHARACTER(LEN=256) :: qexsd_fmt = ' ', qexsd_version = ' '
  LOGICAL            :: qexsd_init = .FALSE. 
  ! ... next two variables are no longer read from input but can be set
  ! ... by external codes using QE routines to perform an interpolation
  ! ... of valence electrons only, without the atomic-like part
  CHARACTER(LEN=256) :: input_drho = ' '
  CHARACTER(LEN=256) :: output_drho= ' '
  !
  CHARACTER(LEN=5 ), PARAMETER :: crash_file  = 'CRASH'
  CHARACTER (LEN=320) :: exit_file = 'os.EXIT'
  !! file required for a soft exit  
  !
  CHARACTER (LEN=20), PARAMETER :: xmlpun_schema = 'data-file-schema.xml'
  !
  ! ... The units where various variables are saved
  ! ... Only units that are kept open during the run should be listed here
  !
  INTEGER :: iunres      =  1
  !! unit for the restart of the run
  INTEGER :: iunpun      =  4
  !! unit for saving the final results (data-file.xml)
  INTEGER :: iunwfc      = 10
  !! unit with wavefunctions
  INTEGER :: iunoldwfc   = 11
  !! unit with old wavefunctions
  INTEGER :: iunoldwfc2  = 12
  !! as above at step -2
  INTEGER :: iunhub      = 13
  !! unit for saving Hubbard-U atomic wfcs * S 
  INTEGER :: iunsat      = 14
  !! unit for saving (orthogonal) atomic wfcs * S
  INTEGER :: iunmix      = 15
  !! unit for saving mixing information
  INTEGER :: iunwfc_exx  = 16
  !! unit with exx wavefunctions
  INTEGER :: iunhub_noS  = 17
  !! unit for saving Hubbard-U atomic wfcs
  !
  INTEGER :: iunexit     = 26
  !! unit for a soft exit  
  INTEGER :: iunupdate   = 27
  !! unit for saving old positions (extrapolation)
  !
  ! NEB
  INTEGER :: iunnewimage = 28
  !! unit for parallelization among images
  INTEGER :: iunlock     = 29
  !! as above (locking file)
  !
  INTEGER :: iuntmp      = 90
  !! temporary unit, when used must be closed ASAP
  !
  INTEGER :: nwordwfc    =  2
  !! length of record in wavefunction file
  INTEGER :: nwordatwfc  =  2
  !! length of record in atomic wfc file
  INTEGER :: nwordwfcU   =  2
  !! length of record in atomic hubbard wfc file
  INTEGER :: nwordwann   =  2
  !! length of record in sic wfc file
  !
  !... finite electric field
  !
  INTEGER :: iunefield   = 31
  !! unit to store wavefunction for calculating electric field operator
  INTEGER :: iunefieldm  = 32
  !! unit to store projectors for hermitean electric field potential
  INTEGER :: iunefieldp  = 33
  !! unit to store projectors for hermitean electric field potential
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
    USE clib_wrappers, ONLY : f_mkdir_safe
    !
    CHARACTER(LEN=*), INTENT(IN) :: dirname
    !
    INTEGER                    :: ierr, length
    !
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    length = LEN_TRIM(dirname)
#if defined (_WIN32)
    ! Windows returns error if tmp_dir ends with a backslash
#if defined (__PGI)
    IF ( dirname(length:length) == '\\' ) length=length-1
#else
    IF ( dirname(length:length) == '\' ) length=length-1
#endif
#endif
    IF ( ionode ) ierr = f_mkdir_safe( dirname(1:length ) )
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
    !! Verify if \(\text{tmp_dir}\) exists, creates it if not.
    !
    USE clib_wrappers, ONLY : f_mkdir_safe
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), INTENT(in) :: tmp_dir
    !! directory to check
    LOGICAL, INTENT(out)         :: exst
    !! TRUE if \(\text{tmp_dir}\) exists
    LOGICAL, INTENT(out)         :: pfs
    !! TRUE if tmp_dir visible from all procs of an image
    !
    INTEGER             :: ios, image, proc, nofi, length
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !
    ! ... create tmp_dir on ionode
    ! ... f_mkdir_safe returns -1 if tmp_dir already exists
    ! ...                       0 if         created
    ! ...                       1 if         cannot be created
    !
    length = LEN_TRIM(tmp_dir)
#if defined (_WIN32)
    ! Windows returns error if tmp_dir ends with a backslash
#if defined (__PGI)
    IF ( tmp_dir(length:length) == '\\' ) length=length-1
#else
    IF ( tmp_dir(length:length) == '\' ) length=length-1
#endif
#endif
    IF ( ionode ) ios = f_mkdir_safe( tmp_dir(1:length) )
    CALL mp_bcast ( ios, ionode_id, intra_image_comm )
    exst = ( ios == -1 )
    IF ( ios > 0 ) CALL errore ('check_tempdir', 'temporary directory ' &
            & // tmp_dir(1:length) // ' cannot be created or accessed',1)
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
    !! Remove temporary files from \(\text{tmp_dir}\) (only by the master node).
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), INTENT(in) :: tmp_dir
    !
    CHARACTER (len=256) :: file_path
    !
    file_path = trim( tmp_dir ) // trim( prefix )
    IF ( ionode ) THEN
       CALL delete_if_present( trim( file_path ) // '.update' )
       CALL delete_if_present( trim( file_path ) // '.md' )
       CALL delete_if_present( trim( file_path ) // '.bfgs' )
       CALL delete_if_present( trim( file_path ) // '.fire' )
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
  SUBROUTINE delete_if_present(filename, para)
  !--------------------------------------------------------------------------
  !!
  !! As the name says - if \(\text{para}\) is present and \(\text{para}\)=.TRUE., 
  !! \(\text{filename}\) is deleted by all cores; otherwise, on ionode only.
  !! (SP - Jan 2020).
  !!  
  !
  IMPLICIT NONE
  !
  CHARACTER(len = *), INTENT(in) :: filename
  !! Name of the file to be deleted
  LOGICAL, OPTIONAL, INTENT(in) :: para
  !! Optionally, the remove can be done by all the cores. 
  ! 
  ! Local variables
  LOGICAL :: exst
  !! Check if the file exist
  INTEGER :: iunit
  !! Unit of the file 
  ! 
  IF (PRESENT(para)) THEN
    IF (.NOT. para) THEN
      IF (.NOT. ionode) RETURN
    ENDIF
  ELSE ! Default if not present
    IF (.NOT. ionode) RETURN
  ENDIF
  !
  INQUIRE(FILE = filename, EXIST = exst)
  !
  IF (exst) THEN
    !
    OPEN(NEWUNIT = iunit, FILE = filename, STATUS = 'OLD')
    CLOSE(UNIT = iunit, STATUS = 'DELETE')
    !
    WRITE(UNIT = stdout, FMT = '(/,5X,"File ", A, " deleted, as requested")') TRIM(filename)
    !
  ENDIF
  !
  RETURN
  ! 
  !--------------------------------------------------------------------------
  END SUBROUTINE delete_if_present
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  FUNCTION check_writable ( file_path, process_id ) RESULT ( ios )
    !--------------------------------------------------------------------------
    !! If run by multiple processes, specific "process_id" to avoid
    !! opening, closing, deleting the same file from different processes.
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
  !------------------------------------------------------------------------
  FUNCTION restart_dir( runit )
    !------------------------------------------------------------------------
    !! Main restart directory (contains final / or Windows equivalent).
    !
    CHARACTER(LEN=256)  :: restart_dir
    INTEGER, INTENT(IN), OPTIONAL :: runit
    !
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    IF ( PRESENT (runit) ) THEN
       restart_dir = TRIM(tmp_dir) // TRIM(prefix) // '_' // &
               TRIM(int_to_char(runit)) // postfix
    ELSE
       restart_dir = TRIM(tmp_dir) // TRIM(prefix) // postfix
    END IF
    !
    RETURN
    !
  END FUNCTION restart_dir
  !
  !------------------------------------------------------------------------
  FUNCTION xmlfile ( runit )
    !------------------------------------------------------------------------
    !! \(\texttt{xml}\) file in main restart directory.
    !
    CHARACTER(LEN=320)  :: xmlfile
    INTEGER, INTENT(IN), OPTIONAL :: runit
    !
    xmlfile = TRIM( restart_dir(runit) ) // xmlpun_schema
    !
    RETURN
    !
  END FUNCTION xmlfile
!
!-----------------------------------------------------------------------
subroutine diropn (unit, extension, recl, exst, tmp_dir_)
  !-----------------------------------------------------------------------
  !! Opens a direct-access file named "prefix"."extension" in directory
  !! \(\text{tmp_dir_}\) if specified, in "tmp\_dir" otherwise.  
  !! In parallel execution, the node number is added to the file name.
  !! The record length is \(\text{recl}\) double-precision numbers.
  !! On output, \(\text{exst}\) is TRUE if opened file already exists.
  !! If \(\text{recl}=-1\), the file existence is checked, nothing else
  !! is done.
  !
  USE kinds, ONLY: i8b
  implicit none
  !
  !    first the input variables
  !
  character(len=*) :: extension
  !! input: name of the file to open
  character(len=*), optional :: tmp_dir_
  !! optional variable, if present it is used as tmp_dir
  integer :: unit
  !! input: unit of the file to open
  integer :: recl
  !! input: length of the records
  logical :: exst
  !! output: if true the file exists
  !
  !    local variables
  !
  character(len=320) :: tempfile
  ! complete file name
  real(dp):: dummy
  integer(kind=i8b) :: unf_recl
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
  if (present(tmp_dir_)) then
     tempfile = trim(tmp_dir_)// trim(prefix) //"."// trim(extension)//nd_nmbr
  else
     tempfile = trim(tmp_dir) // trim(prefix) //"."// trim(extension)//nd_nmbr
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
  !! This routine opens a file named "prefix"."extension"
  !! in \(\text{tmp_dir}\) for sequential I/O access.  
  !! If appropriate, the node number is added to the file name.
  !
  implicit none
  !
  !    first the dummy variables
  !
  character(len=*) :: extension
  !! input: name of the file to connect
  character(len=*) :: formatt
  !! input: 'formatted' or 'unformatted'
  character(len=*), optional :: tmp_dir_
  !! optional variable, if present it is used as tmp\_dir
  integer :: unit
  !! input: unit to connect
  logical :: exst
  !! output: true if the file already exist
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
  !! Direct-access vector input/output.  
  !! read/write \(\text{nword}\) words starting from the address specified by
  !! \(\text{vect}\).
  !
  USE kinds ,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nword
  !! the dimension of vect
  INTEGER, INTENT(IN) :: unit
  !! the unit where to read/write
  INTEGER, INTENT(IN) :: nrec
  !! the record where to read/write
  INTEGER, INTENT(IN) :: io
  !! flag if < 0 reading if > 0 writing
  REAL(DP), INTENT(INOUT) :: vect(nword)
  !! input/output: the vector to read/write
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
         & 'error reading file "' // TRIM(name) // '"', unit )
     !
  ELSE IF ( io > 0 ) THEN
     !
     WRITE( UNIT = unit, REC = nrec, IOSTAT = ios ) vect
     IF ( ios /= 0 ) CALL errore( 'davcio', &
         & 'error writing file "' // TRIM(name) // '"', unit )
     !
  END IF
  !
  CALL stop_clock( 'davcio' )
  !
  RETURN
  !
END SUBROUTINE davcio
