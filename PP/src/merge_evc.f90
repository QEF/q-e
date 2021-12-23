!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Riccardo De Gennaro, EPFL (Nov 2021).
!
!
!---------------------------------------------------------------------
PROGRAM merge_evc
    !-----------------------------------------------------------------
    !
    ! This program merges all the wavefunction files passed as
    ! standard input arguments. The program can be run as below: 
    ! 
    ! merge_evc.x -nr 4 -i path/to/file1 -i path/to/file2 -o output
    ! 
    ! The -nr argument is mandatory and it is followed by an integer
    ! specifying the number of R-vectors, i.e. the number of
    ! repetitions of the unit cell within the supercell
    !
    ! The -i argument is followed by the path to the input file(s)
    ! to be merged
    !
    ! The -o argument is optional and it can be used to specify
    ! the output file; if not present the default output name is
    ! evcw.dat
    !
    ! By adding -v or --verbose as a first argument, the standard
    ! output of a QE program is printed out.
    ! No output is printed by default.
    !
    ! For the moment the format of the files is that defined in
    ! cp_files.f90, used by the Koopmans-CP code.
    !
    ! NB: since this is just a read/write program there is no
    ! parallelization
    !
    !
    USE io_global,              ONLY : stdout, ionode
    USE mp_global,              ONLY : mp_startup
    USE mp_world,               ONLY : nproc
    USE environment,            ONLY : environment_start, environment_end
    USE kinds,                  ONLY : DP
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=20) :: arg
    CHARACTER(LEN=256), ALLOCATABLE :: io_files(:)
    CHARACTER(LEN=256) :: filename, ofilename, nfile
    LOGICAL :: verbosity = .false.
    LOGICAL :: output_exst
    INTEGER :: ifile, ios, n_files
    INTEGER :: iunit = 8
    INTEGER :: ounit = 24
    INTEGER :: iunit_
    INTEGER :: npw, npw_ref
    INTEGER :: ipw, ibnd, ibnd_, ir
    INTEGER :: nbnd, nrtot, istart
    INTEGER, ALLOCATABLE :: nbnd_i(:)
    COMPLEX(DP), ALLOCATABLE :: evc(:)
    !
    INTERFACE
      SUBROUTINE parse_args( istart, io_files, nrtot, output_exst )
        INTEGER, INTENT(IN) :: istart
        CHARACTER(LEN=256), ALLOCATABLE, INTENT(OUT) :: io_files(:)
        INTEGER, INTENT(OUT) :: nrtot
        LOGICAL, INTENT(OUT) :: output_exst
      END SUBROUTINE
    END INTERFACE
    !
    !
#if defined(__MPI)
    CALL mp_startup( )
#endif
    !
    CALL get_command_argument( 1, arg )
    IF ( trim(arg) == '-v' .or. trim(arg) == '--verbose' ) verbosity = .true.
    !
    IF ( verbosity ) THEN
      CALL environment_start( 'MERGE_EVC' )
      WRITE( stdout, * )
      istart = 2
    ELSE
      istart = 1
    ENDIF
    !
!    IF ( nproc .gt. 1 ) &
!      WRITE( stdout, fmt='(5X,A)' ) "WARNING: Parallelization not supported"
    !
    CALL parse_args( istart, io_files, nrtot, output_exst )
    !
    IF ( output_exst ) THEN
      ofilename = trim( io_files(size(io_files)) )
      n_files = size( io_files ) - 1
    ELSE
      ofilename = 'evcw.dat'
      n_files = size( io_files )
    ENDIF
    !
    IF ( ionode ) THEN
      !
      ! ... reading number of PWs and number of bands from each file
      !
      nbnd = 0
      ALLOCATE( nbnd_i(n_files) )
      !
      DO ifile = 1, n_files
        !
        iunit_ = iunit + ifile
        filename = trim(io_files(ifile))
        OPEN( UNIT=iunit_, FILE=trim(filename), STATUS='old', FORM='unformatted', IOSTAT=ios )
        !
        IF ( ios .ne. 0 ) THEN
          CLOSE( UNIT=ounit, STATUS='delete' )
          CALL errore( 'merge_evc', 'problem opening the input files', ifile )
        ENDIF
        !
        READ( iunit_ ) npw, nbnd_i(ifile)
        !
        IF ( mod( nbnd_i(ifile), nrtot ) .ne. 0 ) &
          CALL errore( 'merge_evc', 'incosistent number of R-vectors', mod( nbnd_i(ifile), nrtot ) )
        !
        IF ( ifile == 1 ) npw_ref = npw
        !
        IF ( npw .ne. npw_ref ) THEN
          CLOSE( UNIT=ounit, STATUS='delete' )
          CALL errore( 'merge_evc', 'unconsistent number of PW between files', ifile )
        ENDIF
        !
      ENDDO
      !
      nbnd = SUM( nbnd_i )
      OPEN( UNIT=ounit, FILE=trim(ofilename), STATUS='unknown', FORM='unformatted' )
      WRITE( ounit ) npw, nbnd
      ALLOCATE( evc(npw) )
      !
      ! ... reading wave functions and writing them to a single file
      !
      DO ir = 1, nrtot
        !
        DO ifile = 1, n_files
          !
          iunit_ = iunit + ifile
          !
          IF ( ir == 1 ) THEN
            WRITE( nfile, 100 ) ifile
            CALL start_clock( nfile )
          ENDIF
          !
          DO ibnd = 1, nbnd_i(ifile)/nrtot
            !
            evc(:) = ( 0.d0, 0.d0 )
            READ ( iunit_ ) ( evc(ipw), ipw = 1, npw )
            WRITE( ounit ) ( evc(ipw), ipw = 1, npw )
            !
          ENDDO
          !
          IF ( ir == nrtot ) THEN
            CLOSE( UNIT=iunit_ )
            WRITE( nfile, 100 ) ifile
            CALL stop_clock( nfile )
            IF ( verbosity ) CALL print_clock( nfile )
          ENDIF
          !
        ENDDO
        !
      ENDDO
      !
      CLOSE( UNIT=ounit )
      !
    ENDIF
    !
    !
    IF ( verbosity ) THEN
      WRITE( stdout, * )
      WRITE( stdout, fmt='(5X,A,A)' ) 'Merged wave functions written into : ', trim(ofilename)
      CALL environment_end( 'MERGE_EVC' )
      WRITE( stdout, * )
    ENDIF
    !
    !
100 FORMAT( 'r/w file', 1X, I2 )
    !
    !
END PROGRAM merge_evc
!
!
!---------------------------------------------------------------------
SUBROUTINE parse_args( istart, io_files, nrtot, output_exst )
    !-----------------------------------------------------------------
    !
    ! ... routine parsing arguments passed in input
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: istart
    CHARACTER(LEN=256), ALLOCATABLE, INTENT(OUT) :: io_files(:)
    INTEGER, INTENT(OUT) :: nrtot
    LOGICAL, INTENT(OUT) :: output_exst
    !
    INTEGER :: nargs, iarg, n_files
    INTEGER :: counter, ios
    CHARACTER(LEN=256) :: arg
    LOGICAL :: input_exst = .false.
    LOGICAL :: nrtot_exst = .false.
    !
    !
    nargs = command_argument_count()
    !   
    IF ( nargs - istart < 3 .or. mod( nargs-istart, 2 ) == 0 ) &
      CALL errore( 'parse_args', 'wrong number of positional arguments', 1 )
    !
    output_exst = .false.
    n_files = ( nargs - istart - 1 ) / 2
    ALLOCATE( io_files(n_files) )
    io_files(:) = ' '
    iarg = istart
    counter = 0
    !
    DO WHILE ( iarg < nargs )
      !
      CALL get_command_argument( iarg, arg )
      iarg = iarg + 1
      !
      SELECT CASE ( trim(arg) )
      !
      CASE ( '-nr' )
        CALL get_command_argument( iarg, arg )
        READ( arg, *, IOSTAT=ios )  nrtot
        IF ( ios .ne. 0 ) &
          CALL errore( 'parse_args', 'error while reading number of R-vectors', abs(ios) )
        IF ( nrtot_exst ) &
          CALL errore( 'parse_args', 'nrtot can be defined only once', 1 )
        nrtot_exst = .true.
      !
      CASE ( '-i', '-in', '-input' )
        CALL get_command_argument( iarg, arg )
        counter = counter + 1
        IF ( counter > n_files ) &
          CALL errore( 'parse_args', 'something wrong in input. Did you provide -nr ?', 1 )
        io_files(counter) = trim(arg)
        IF ( .not. input_exst ) input_exst = .true.
      !
      CASE ( '-o', '-out', '-output' )
        CALL get_command_argument( iarg, arg )
        IF ( output_exst ) &
          CALL errore( 'parse_args', 'only one output file can be defined', 1 )
        output_exst = .true.
        io_files(n_files) = trim(arg)
      !
      CASE DEFAULT
        CALL errore( 'parse_args', 'unrecognised argument option', abs(iarg-1) )
      !
      END SELECT
      !
      iarg = iarg + 1
      !
    ENDDO
    !
    IF ( .not. nrtot_exst ) CALL errore( 'parse_args', 'nrtot was not provided', 1)
    !
    IF ( .not. input_exst ) CALL errore( 'parse_args', 'no input files provided', 1 )
    IF ( nrtot .le. 0 ) CALL errore( 'parse_args', 'wrong number of R-vectors', abs(nrtot) )
    !
    !
END SUBROUTINE parse_args