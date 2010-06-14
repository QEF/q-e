!
! Copyright (C) 2010 A. Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!========================
  PROGRAM qexml
  !========================
  !
  ! Simple example to show how to use the QEXML library
  ! to read data from the .save directory written by QE
  !
  ! General comments:
  !
  !   - first init the library
  !
  !   - for each data required, first get the
  !     dimensions, then allocate the target array
  !     and finally read the data with a second call
  !     to the proper qexml read routine
  !     (shown below)
  !
  !   - data that don't need any dynamical allocation
  !     (scalar or small arrays) can be read directly
  !
  !   - explicit error handling through the use of IERR arguments
  !     is required
  !
  USE qexml_module
  IMPLICIT NONE

  !
  ! parameters
  !
  INTEGER, PARAMETER :: iunit = 10
  INTEGER, PARAMETER :: stdin = 5
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: DP=kind(1.0d0)

  !
  ! input variables
  !
  CHARACTER(256)     :: prefix       ! used to locate the .save directory
  CHARACTER(256)     :: work_dir     !
  !
  INTEGER            :: ik           ! indexes of the kpt and band wfc
  INTEGER            :: ib           ! to be read

  NAMELIST /INPUT/ prefix, work_dir, ik, ib


  !
  ! local vars
  !
  CHARACTER(7) :: subname='example'
  !
  REAL(DP)     :: avec(3,3), alat
  INTEGER      :: npw, ngm
  !
  CHARACTER(256)           :: dirname, filename, str_units
  INTEGER,     ALLOCATABLE :: igv(:,:), igk(:)
  COMPLEX(DP), ALLOCATABLE :: wfc(:,:)
  !
  INTEGER      :: ierr, ig

!
!----------------------------------------
! main Body
!----------------------------------------
!

  WRITE(stdout, "(/,'< QEXML example> ')" )

  !
  ! read stdin
  !
  WRITE(stdout, "(/,'Reading input namelist...')" )
  !
  ik = 0
  ib = 0
  prefix = ' '
  work_dir = './'
  !
  READ( stdin, INPUT, IOSTAT=ierr)
  IF ( ierr/=0 ) CALL errore(subname,'reading INPUT namelist',abs(ierr))
  !

  !
  !==========================
  ! init the qexml library
  !==========================
  !
  WRITE(stdout, "(/, 'Init QEXML library...')" )
  !
  dirname = trim(work_dir) // '/' // trim(prefix) // '.save/'
  CALL qexml_init( iunit, DIR=dirname )

  filename = trim(dirname) // "data-file.xml"
  !
  CALL qexml_openfile( filename, "read", IERR=ierr )
  IF ( ierr/=0) CALL errore(subname,'opening dftdata file',abs(ierr))


  !
  !==========================
  ! read lattice data
  !==========================
  ! how to read data directly
  ! units can be read as well
  !
  WRITE(stdout, "(/, 'Read direct lattice data...')" )
  !
  CALL qexml_read_cell( ALAT=alat, &
                        A1=avec(:,1), A2=avec(:,2), A3=avec(:,3), &
                        A_UNITS=str_units, IERR=ierr)
  IF (ierr/=0) CALL errore(subname,'reading lattice',abs(ierr))

  !
  ! reports to stdout
  !
  WRITE(stdout, "(2x,' Direct lattice  [',a,']')") trim(str_units)
  WRITE(stdout, "(2x,' alat:  ',f15.9)") alat
  WRITE(stdout, "(2x,' a(1):  ',3f15.9)") avec(:,1)
  WRITE(stdout, "(2x,' a(2):  ',3f15.9)") avec(:,2)
  WRITE(stdout, "(2x,' a(3):  ',3f15.9)") avec(:,3)


  !
  !==========================
  ! read G grid for the density
  !==========================
  !
  ! First we read the whole G vectors grid
  ! For each kpt, a map of the corresponding G-vectors to
  ! those in the density G grid is then given
  !
  ! - read the dimensions
  ! - allocate
  ! - read the massive data
  !
  WRITE(stdout, "(/, 'Read main G grid...')" )
  !
  CALL qexml_read_planewaves( NGM=ngm, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading PW dims',abs(ierr))
  !
  ALLOCATE( igv(3,ngm), STAT=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'allocating igv',abs(ierr))
  !
  CALL qexml_read_planewaves( IGV=igv, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading main G grid',abs(ierr))
  !
  !
  WRITE(stdout, "(2x, ' Main grid dim: ',i5)" ) ngm
  WRITE(stdout, "(2x, ' Reporting the first 10 elements (check gvectors.dat)')" )
  DO ig = 1, 10
     WRITE(stdout, "(2x, ' ig(',i3, ') : ',3i5 )" ) ig, igv(:,ig)
  ENDDO


  !
  ! now read data specific to a given kpt
  !
  WRITE(stdout, "(/, 'Read ik-specific dims...')" )
  !
  CALL qexml_read_gk( ik, NPWK=npw, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ik dims',abs(ierr))
  !
  ALLOCATE( igk(npw), STAT=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'allocating igk',abs(ierr))
  !
  ! the second dimension is the # of bands to be read
  ALLOCATE( wfc(npw,1), STAT=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'allocating wfc',abs(ierr))
  !
  CALL qexml_read_gk( ik, index=igk, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading k-grid map',abs(ierr))
  !
  WRITE(stdout, "(2x, ' ik:',i3,'   dim: ',i5)" ) ik, npw

  !
  !==========================
  ! read wavefunctions
  !==========================
  !
  ! Here we do not need to get any dimension since
  ! we already have all we need
  !
  WRITE(stdout, "(/, 'Read a given wfc...')" )
  !
  CALL qexml_read_wfc( IBNDS=ib, IBNDE=ib, IK=ik, WF=wfc, IERR=ierr)
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',abs(ierr))
  !
  ! report to stdout
  !
  WRITE(stdout, "(2x, ' ik:',i3,'   ib: ',i3)" ) ik, ib
  WRITE(stdout, "(2x, ' Reporting the first 10 elements (check evc.dat)')" )
  DO ig = 1, 10
     WRITE(stdout, "(2x, ' ig(',i3, ') : ',2f15.9 )" ) ig, wfc(ig,1)
  ENDDO



  !
  !==========================
  ! finalize the qexml read
  !==========================
  !
  WRITE(stdout, "(/,'Finalize QEXML...')" )
  !
  CALL qexml_closefile ( "read", IERR=ierr )
  IF ( ierr/=0) CALL errore(subname,'closing dftdata file',abs(ierr))


  !
  ! local cleanup
  !
  DEALLOCATE( igv, igk, wfc)

CONTAINS

!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE errore( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an error message to output:
  ! ... if ierr <= 0 it does nothing,
  ! ... if ierr  > 0 it stops.
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: calling_routine, message
    ! the name of the calling calling_routinee
    ! the output messagee
  INTEGER,          INTENT(in) :: ierr
    ! the error flag
    !
  !
  !
  IF ( ierr <= 0 ) RETURN
  !
  ! ... the error message is written un the "*" unit
  !
  WRITE( UNIT = *, FMT = '(/,1X,78("%"))' )
  WRITE( UNIT = *, &
         FMT = '(5X,"from ",A," : error #",I10)' ) calling_routine, ierr
  WRITE( UNIT = *, FMT = '(5X,A)' ) message
  WRITE( UNIT = *, FMT = '(1X,78("%"),/)' )
  !
  WRITE( *, '("     stopping ...")' )
  !
  STOP 2
  !
  RETURN
  !
END SUBROUTINE errore

END PROGRAM qexml

