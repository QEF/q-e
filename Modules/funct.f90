!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
MODULE funct
  !-------------------------------------------------------------------
  !! This module contains data defining the DFT functional in use
  !! and a number of functions and subroutines to manage them.
  !! All the data and routines related to LDA, GGA and MGGA 
  !! functionals have been moved into the library XClib.  
  !! Here the combinations with nonlocal functionals are still
  !! managed.
  !
  ! Data are PRIVATE and are accessed and set only by function calls.
  !
  !
  USE io_global,      ONLY: stdout, ionode
  USE kinds,          ONLY: DP
  USE beef_interface, ONLY: beef_set_type
  USE xc_interfaces
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  ! subroutines/functions managing dft name and indices
  PUBLIC :: set_dft_from_indices, set_dft_from_name
  PUBLIC :: enforce_input_dft, write_dft_name
  PUBLIC :: get_dft_name
  PUBLIC :: get_dft_short, get_dft_long
  PUBLIC :: get_nonlocc_name
  PUBLIC :: get_inlc
  PUBLIC :: dft_is_nonlocc
  ! driver subroutine computing XC non local
  PUBLIC  :: nlc
  ! XC non local index
  PRIVATE :: inlc
  !
  CHARACTER(LEN=25) :: dft = 'not set'
  !
  ! ------------------------------------------------------------------------
  !
  INTEGER, PARAMETER :: notset = -1
  !
  INTEGER :: inlc = notset
  !
  LOGICAL :: isnonlocc = .FALSE.
  !
  LOGICAL :: discard_input_dft = .FALSE.
  !
  INTEGER  :: beeftype = -1
  INTEGER  :: beefvdw = 0
  !
  INTEGER, PARAMETER :: ncnl = 26
  CHARACTER(LEN=4) :: nonlocc
  DIMENSION :: nonlocc(0:ncnl)
  !
  DATA nonlocc/ 'NONE', 'VDW1', 'VDW2', 'W31C', 'W32C', 'WC6', 20*'NONE', 'VV10' /
  !
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE set_dft_from_name( dft_ )
    !-----------------------------------------------------------------------
    !! It sets the dft functional IDs and parameters from the input name. It 
    !! directly calls the XClib routines and functions to set the LDA, GGA,
    !! MGGA terms (but not the non-local one).
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    !
    ! ... local variables
    !
    INTEGER :: len, l, i
    CHARACTER(len=150) :: dftout
    LOGICAL :: dft_defined
    LOGICAL :: check_libxc
    !
    CHARACTER(LEN=1), EXTERNAL :: capital
    CHARACTER(LEN=4) :: lda_exch, lda_corr, gga_exch, gga_corr
    !
    INTEGER :: save_inlc
    INTEGER :: iexch, icorr, igcx, igcc, imeta
    !
    ! Exit if set to discard further input dft
    !
    IF ( discard_input_dft ) RETURN
    !
    ! save current status of XC indices
    !
    dft_defined = .FALSE.
    !
    save_inlc  = inlc
    !
    ! convert to uppercase
    !
    len = LEN_TRIM(dft_)
    dftout = ' '
    !
    DO l = 1, len
       dftout(l:l) = capital( dft_(l:l) )
    ENDDO
    !
    !
    ! ----------------------------------------------
    ! NOW WE CHECK ALL THE SHORT NAMES
    ! Note: comparison is done via exact matching
    ! ----------------------------------------------
    !
    SELECT CASE( TRIM(dftout) )
    ! special case : case BEEF (default: BEEF-vdW-DF2)
    CASE('BEEF', 'BEEF-VDW')
       IF (LEN_TRIM(dftout) == 4) THEN
          beeftype = 0
       ELSE
          SELECT CASE(TRIM(dftout(5:)))
             CASE('-VDW')
                beeftype = 0
             CASE DEFAULT
                READ(dftout(5:), '(i1)', IOSTAT=i) beeftype
                IF (i /= 0) CALL errore('set_dft_from_name', &
                                      & 'unknown BEEF type', 1)
          END SELECT
       ENDIF
       IF (.NOT. beef_set_type(beeftype, ionode)) &
       & CALL errore('set_dft_from_name', 'unknown BEEF type number', 1)
       SELECT CASE(beeftype)
          CASE(0)
             ! turn on vdW-DF2 type interactions for BEEF-vdW
             beefvdw = 2
       END SELECT
       dft_defined = xclib_set_dft_IDs(1,4,43,14,0,0)
       inlc = beefvdw
    ! Special case vdW-DF
    CASE( 'VDW-DF' )
       dft_defined = xclib_set_dft_IDs(1,4,4,0,0,0)
       inlc = 1
    ! Special case vdW-DF2
    CASE( 'VDW-DF2' )
       dft_defined = xclib_set_dft_IDs(1,4,13,0,0,0)
       inlc = 2
    ! Special case vdW-DF3-opt1
    CASE( 'VDW-DF3-OPT1' )
       dft_defined = xclib_set_dft_IDs(1,4,45,0,0,0)
       inlc = 3
    ! Special case vdW-DF3-opt2
    CASE( 'VDW-DF3-OPT2' )
       dft_defined = xclib_set_dft_IDs(1,4,46,0,0,0)
       inlc = 4
    ! Special case vdW-DF-C6
    CASE( 'VDW-DF-C6' )
       dft_defined = xclib_set_dft_IDs(1,4,26,0,0,0)
       inlc = 5
    ! Special case vdW-DF with C09 exchange
    CASE( 'VDW-DF-C09' )
       dft_defined = xclib_set_dft_IDs(1,4,16,0,0,0)
       inlc = 1
    ! Special case vdW-DF2 with C09 exchange
    CASE( 'VDW-DF2-C09' )
       dft_defined = xclib_set_dft_IDs(1,4,16,0,0,0)
       inlc = 2
    ! Special case vdW-DF-obk8, or vdW-DF + optB88
    CASE( 'VDW-DF-OBK8' )
       dft_defined = xclib_set_dft_IDs(1,4,23,0,0,0)
       inlc = 1
    ! Special case vdW-DF-ob86, or vdW-DF + optB86
    CASE( 'VDW-DF-OB86' )
       dft_defined = xclib_set_dft_IDs(1,4,24,0,0,0)
       inlc = 1
    ! Special case vdW-DF2 with B86R
    CASE( 'VDW-DF2-B86R' )
       dft_defined = xclib_set_dft_IDs(1,4,26,0,0,0)
       inlc = 2
    ! Special case vdW-DF-CX
    CASE( 'VDW-DF-CX' )
       dft_defined = xclib_set_dft_IDs(1,4,27,0,0,0)
       inlc = 1
    ! Special case vdW-DF-CX0
    CASE( 'VDW-DF-CX0' )
       dft_defined = xclib_set_dft_IDs(6,4,29,0,0,0)
       inlc = 1
    ! Special case vdW-DF-CX0P
    CASE( 'VDW-DF-CX0P' )
       dft_defined = xclib_set_dft_IDs(6,4,31,0,0,0)
       inlc = 1
    ! Special case vdW-DF2-0
    CASE( 'VDW-DF2-0' )
       dft_defined = xclib_set_dft_IDs(6,4,30,0,0,0)
       inlc = 2
    ! Special case vdW-DF2-BR0
    CASE( 'VDW-DF2-BR0' )
       dft_defined = xclib_set_dft_IDs(6,4,38,0,0,0)
       inlc = 2
    ! Special case vdW-DF-C090
    CASE( 'VDW-DF-C090' )
       dft_defined = xclib_set_dft_IDs(6,4,40,0,0,0)
       inlc = 1
    ! Special case rVV10
    CASE( 'RVV10' )
       dft_defined = xclib_set_dft_IDs(1,4,13,4,0,0)
       inlc = 26
    ! Special case rVV10+scan
    CASE( 'RVV10-SCAN' )
       dft_defined = xclib_set_dft_IDs(0,0,0,0,5,0)
       inlc = 26
    !
    CASE( 'REV-VDW-DF2' )
       CALL errore( 'set_dft_from_name', 'obsolete XC label, use VDW-DF2-B86R', 1 )
    !
    CASE( 'VDW-DF3' )
       CALL errore( 'set_dft_from_name', 'obsolete XC label, use VDW-DF-OBK8', 1 )
    !
    CASE( 'VDW-DF4', 'OPTB86B-VDW' )
       CALL errore( 'set_dft_from_name', 'obsolete XC label, use VDW-DF-OB86', 1 )
    ! Special case vdW-DF-X
    CASE( 'VDW-DF-X' )
       CALL errore( 'set_dft_from_name', 'functional not yet implemented', 1 )
    ! Special case vdW-DF-Y
    CASE( 'VDW-DF-Y' )
       CALL errore( 'set_dft_from_name', 'functional not yet implemented', 1 )
    ! Special case vdW-DF-Z
    CASE( 'VDW-DF-Z' )
       CALL errore( 'set_dft_from_name', 'functional not yet implemented', 1 )
    ! Case for old RRKJ format, containing indices instead of label
    CASE DEFAULT
       !
       IF ('INDEX:' ==  dftout(1:6)) THEN
          READ( dftout(7:18), '(6i2)') iexch, icorr, igcx, igcc, inlc, imeta
          dft_defined = xclib_set_dft_IDs(iexch, icorr, igcx, igcc, imeta, 0)
          CALL xclib_get_name('LDA','EXCH', lda_exch)
          CALL xclib_get_name('LDA','CORR', lda_corr)
          CALL xclib_get_name('GGA','EXCH', gga_exch)
          CALL xclib_get_name('GGA','CORR', gga_corr)
          !
          dftout = TRIM(lda_exch) //'-'// &
                   TRIM(lda_corr) //'-'// &
                   TRIM(gga_exch) //'-'// &
                   TRIM(gga_corr) //'-'// nonlocc(inlc)
       ELSE
          CALL xclib_set_dft_from_name( TRIM(dftout) )
          inlc = matching( dftout, ncnl, nonlocc )
          dft_defined = .TRUE.
       ENDIF
       !
    END SELECT
    !
    !----------------------------------------------------------------
    ! Last check
    ! No more defaults, the code exits if the dft is not defined
    !----------------------------------------------------------------
    !
    iexch = xclib_get_id('LDA','EXCH')
    icorr = xclib_get_id('LDA','CORR')
    igcx  = xclib_get_id('GGA','EXCH')
    igcc  = xclib_get_id('GGA','CORR')
    imeta = xclib_get_id('MGGA','EXCH')
    !
    IF (igcx == 6) CALL infomsg( 'set_dft_from_name', 'OPTX untested! please test' )
    !
    ! check for unrecognized labels
    !
    IF ( iexch<=0 .AND. icorr<=0 .AND. igcx<=0 .AND. igcc<=0 .AND. imeta<=0 ) THEN
       IF ( inlc <= 0 .AND. TRIM(dftout) /= 'NOX-NOC') THEN
          CALL errore( 'set_dft_from_name', TRIM(dftout)//': unrecognized dft', 1 )
       ELSE
          ! if inlc is the only nonzero index the label is likely wrong
          CALL errore( 'set_dft_from_name', TRIM(dftout)//': strange dft, please check', inlc )
       ENDIF
    ENDIF
    !
    ! Fill variables and exit
    !
    dft = dftout
    !
    dft_defined = .TRUE.
    !
    isnonlocc = (inlc > 0)
    !
    CALL xclib_set_auxiliary_flags
    !
    ! check non-local term has not been previously set differently
    !
    IF (save_inlc /= notset  .AND. save_inlc /= inlc)   THEN
       WRITE (stdout,*) inlc, save_inlc
       CALL errore( 'set_dft_from_name', ' conflicting values for inlc', 1 )
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE set_dft_from_name
  !
  !
  !-----------------------------------------------------------------
  INTEGER FUNCTION matching( dft, n, name )
    !-----------------------------------------------------------------
    !! Looks for matches between the names of each single term of the 
    !! xc-functional and the input dft string.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN):: n
    CHARACTER(LEN=*), INTENT(IN):: name(0:n)
    CHARACTER(LEN=*), INTENT(IN):: dft
    INTEGER :: i
    LOGICAL, EXTERNAL :: matches
    !
    matching = notset
    !
    DO i = n, 0, -1
       IF ( matches(name(i), TRIM(dft)) ) THEN
          !
          IF ( matching == notset ) THEN
           !WRITE(*, '("matches",i2,2X,A,2X,A)') i, name(i), TRIM(dft)
           matching = i
          ELSE
             WRITE(*, '(2(2X,i2,2X,A))') i, TRIM(name(i)), &
                                  matching, TRIM(name(matching))
             CALL errore( 'set_dft', 'two conflicting matching values', 1 )
          ENDIF
       ENDIF
    ENDDO
    !
    IF (matching == notset) matching = 0
    !
  END FUNCTION matching
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE enforce_input_dft( dft_, nomsg )
    !---------------------------------------------------------------------
    !! Translates a string containing the exchange-correlation name
    !! into internal indices and force any subsequent call to 
    !! \(\textrm{set_dft_from_name}\) to return without changing them.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    LOGICAL, INTENT(IN), OPTIONAL :: nomsg
    !
    CALL set_dft_from_name( dft_ )
    IF (dft == 'not set') CALL errore( 'enforce_input_dft', 'cannot fix unset dft', 1 )
    discard_input_dft = .TRUE.
    !
    IF ( PRESENT(nomsg) ) RETURN
    !
    WRITE(stdout,'(/,5x,a)') "IMPORTANT: XC functional enforced from input :"
    CALL write_dft_name
    WRITE(stdout,'(5x,a)') "Any further DFT definition will be discarded"
    WRITE(stdout,'(5x,a/)') "Please, verify this is what you really want"
    !
    RETURN
    !
  END SUBROUTINE enforce_input_dft
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION get_inlc()
    !! Get dft index for non-local term.
    INTEGER :: get_inlc
    get_inlc = inlc
    RETURN
  END FUNCTION get_inlc
  !-----------------------------------------------------------------------
  FUNCTION get_nonlocc_name()
    !! Get dft name for non-local term.
    CHARACTER(10) get_nonlocc_name
    get_nonlocc_name = TRIM(nonlocc(inlc))
    RETURN
  END FUNCTION get_nonlocc_name
  !-----------------------------------------------------------------------
  FUNCTION dft_is_nonlocc()
    !! TRUE if dft is non-local.
    LOGICAL :: dft_is_nonlocc
    dft_is_nonlocc = isnonlocc
    RETURN
  END FUNCTION dft_is_nonlocc
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION get_dft_name()
    !! Get the string with the full dft name.
    CHARACTER(LEN=25) :: get_dft_name
    get_dft_name = dft
    RETURN
  END FUNCTION get_dft_name
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE set_dft_from_indices( iexch_, icorr_, igcx_, igcc_, imeta_, inlc_ )
     !--------------------------------------------------------------------
     !! Set dft functional from the IDs of each term.
     !
     IMPLICIT NONE
     !
     INTEGER :: iexch_, icorr_, igcx_, igcc_, imeta_, inlc_
     INTEGER :: iexch, icorr, igcx, igcc, imeta
     CHARACTER(LEN=4) :: lda_exch, lda_corr, gga_exch, gga_corr
     LOGICAL :: dft_defined
     !
     iexch  = xclib_get_id( 'LDA', 'EXCH' )
     icorr  = xclib_get_id( 'LDA', 'CORR' )
     igcx   = xclib_get_id( 'GGA', 'EXCH' )
     igcc   = xclib_get_id( 'GGA', 'CORR' )
     imeta  = xclib_get_id( 'MGGA','EXCH' )
     !
     IF ( discard_input_dft ) RETURN
     !
     IF (iexch == notset) iexch = iexch_
     IF (iexch /= iexch_) THEN
        write (stdout,*) iexch, iexch_
        CALL errore( 'set_dft', ' conflicting values for iexch', 1 )
     ENDIF
     IF (icorr == notset) icorr = icorr_
     IF (icorr /= icorr_) THEN
        write (stdout,*) icorr, icorr_
        CALL errore( 'set_dft', ' conflicting values for icorr', 1 )
     ENDIF
     IF (igcx  == notset) igcx = igcx_
     IF (igcx /= igcx_) THEN
        write (stdout,*) igcx, igcx_
        CALL errore( 'set_dft', ' conflicting values for igcx', 1 )
     ENDIF
     IF (igcc  == notset) igcc = igcc_
     IF (igcc /= igcc_) THEN
        write (stdout,*) igcc, igcc_
        CALL errore( 'set_dft', ' conflicting values for igcc', 1 )
     ENDIF
     IF (imeta  == notset) imeta = imeta_
     IF (imeta /= imeta_) THEN
        write (stdout,*) imeta, imeta_
        CALL errore( 'set_dft', ' conflicting values for imeta', 1 )
     ENDIF     
     IF (inlc  == notset) inlc = inlc_
     IF (inlc /= inlc_) THEN
        write (stdout,*) inlc, inlc_
        CALL errore( 'set_dft', ' conflicting values for inlc', 1 )
     ENDIF
     CALL xclib_get_name('LDA','EXCH', lda_exch)
     CALL xclib_get_name('LDA','CORR', lda_corr)
     CALL xclib_get_name('GGA','EXCH', gga_exch)
     CALL xclib_get_name('GGA','CORR', gga_corr)
     !
     dft = TRIM(lda_exch) //'-'// &
           TRIM(lda_corr) //'-'// &
           TRIM(gga_exch) //'-'// &
           TRIM(gga_corr) //'-'// nonlocc(inlc)
     !
     dft_defined = xclib_set_dft_IDs(iexch,icorr,igcx,igcc,imeta,0)
     !
     ! WRITE( stdout,'(a)') dft
     CALL xclib_set_auxiliary_flags
     RETURN
  END SUBROUTINE set_dft_from_indices
  !
  !
  !-------------------------------------------------------------------------
  FUNCTION get_dft_short()
    !---------------------------------------------------------------------
    !! It gets a short version (if exists) of the name of the dft in use.  
    !! If there is no non-local term directly calls the xclib analogous 
    !! routine (\(\texttt{xclib_get_dft_short}\)).
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=26) :: get_dft_short
    CHARACTER(LEN=26) :: shortname
    INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac
    !
    shortname = 'no shortname'
    !
    IF (inlc == 0) THEN
      shortname = xclib_get_dft_short()
    ELSE
      !
      iexch  = xclib_get_id( 'LDA', 'EXCH' )
      icorr  = xclib_get_id( 'LDA', 'CORR' )
      igcx   = xclib_get_id( 'GGA', 'EXCH' )
      igcc   = xclib_get_id( 'GGA', 'CORR' )
      !
      IF (inlc==1) THEN
        !
        IF (iexch==1 .AND. icorr==4 .AND. igcx==4 .AND. igcc==0) THEN
           shortname = 'VDW-DF'
        ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==27 .AND. igcc==0) THEN
           shortname = 'VDW-DF-CX'
        ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==29 .AND. igcc==0) THEN
           shortname = 'VDW-DF-CX0'
        ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==31 .AND. igcc==0) THEN
           shortname = 'VDW-DF-CX0P'
        ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==16 .AND. igcc==0) THEN
           shortname = 'VDW-DF-C09'
        ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==24 .AND. igcc==0) THEN
           shortname = 'VDW-DF-OB86'
        ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==23 .AND. igcc==0) THEN
           shortname = 'VDW-DF-OBK8'
        ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==40 .AND. igcc==0) THEN
           shortname = 'VDW-DF-C090'
        ENDIF
        !
      ELSEIF (inlc==2) THEN
        !
        IF (iexch==1 .AND. icorr==4  .AND. igcx==43 .AND. igcc==14) THEN
           shortname = 'BEEF'
        ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==13 .AND. igcc==0) THEN
           shortname = 'VDW-DF2'
        ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==16 .AND. igcc==0) THEN
           shortname = 'VDW-DF2-C09'
        ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==26 .AND. igcc==0) THEN
           shortname = 'VDW-DF2-B86R'
        ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==30 .AND. igcc==0) THEN
           shortname = 'VDW-DF2-0'
        ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==38 .AND. igcc==0) THEN
           shortname = 'VDW-DF2-BR0'
        ENDIF
        !
      ELSEIF (inlc==3) THEN
        !
        shortname = 'VDW-DF3-OPT1'
        !
      ELSEIF (inlc==4) THEN
        !
        shortname = 'VDW-DF3-OPT2'
        !
      ELSEIF (inlc==5) THEN
        !
        shortname = 'VDW-DF-C6'
        !
      ELSEIF (inlc==26) THEN
        !
        shortname = 'RVV10'
        !
      ENDIF
      !
    ENDIF
    !
    get_dft_short = shortname
    !
  END FUNCTION get_dft_short
  !
  !
  !---------------------------------------------------------------------
  FUNCTION get_dft_long()
    !---------------------------------------------------------------------
    !! Returns a string containing the name of each term of the dft functional.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=25) :: get_dft_long
    CHARACTER(LEN=25) :: longname
    !
    !WRITE(longname,'(4a5)') exc(iexch), corr(icorr), gradx(igcx), gradc(igcc)
    !
    longname = xclib_get_dft_long()
    !
    IF ( inlc > 0 ) longname = longname(1:20)//TRIM(nonlocc(inlc))
    !
    get_dft_long = longname
    !
  END FUNCTION get_dft_long
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE write_dft_name
    !-----------------------------------------------------------------------
    !! Print on output the name of each term of the dft functional.
    !
    IMPLICIT NONE
    !
    INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac
    !
    WRITE( stdout, '(5X,"Exchange-correlation= ",A)') TRIM( dft )
    iexch  = xclib_get_id( 'LDA', 'EXCH' )
    icorr  = xclib_get_id( 'LDA', 'CORR' )
    igcx   = xclib_get_id( 'GGA', 'EXCH' )
    igcc   = xclib_get_id( 'GGA', 'CORR' )
    imeta  = xclib_get_id( 'MGGA','EXCH' )
    imetac = xclib_get_id( 'MGGA','CORR' )
    !
    WRITE( stdout, '(27X,"(",I4,3I4,3I4,")")' ) iexch, icorr, igcx, igcc, inlc, &
                                                imeta, imetac
    IF ( xclib_get_exx_fraction() > 0.0_dp ) WRITE( stdout, &
         '(5X,"EXX-fraction              =",F12.2)') xclib_get_exx_fraction()
    RETURN
  END SUBROUTINE write_dft_name
  !
  !
  !-----------------------------------------------------------------------
  !------- NONLOCAL CORRECTIONS DRIVER ----------------------------------
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE nlc( rho_valence, rho_core, nspin, enl, vnl, v )
    !-----------------------------------------------------------------------
    !! Non-local contribution to the correlation energy.
    !
    !     input      :  rho_valence, rho_core
    !     definition :  E_nl = \int E_nl(rho',grho',rho'',grho'',|r'-r''|) dr
    !     output     :  enl = E^nl_c
    !                   vnl = D(E^nl_c)/D(rho)
    !                   v   = non-local contribution to the potential
    !
    !
    USE vdW_DF, ONLY: xc_vdW_DF, xc_vdW_DF_spin, inlc_ => inlc
    USE rVV10,  ONLY: xc_rVV10
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)    :: rho_valence(:,:), rho_core(:)
    INTEGER,  INTENT(IN)    :: nspin
    REAL(DP), INTENT(INOUT) :: v(:,:)
    REAL(DP), INTENT(INOUT) :: enl, vnl
    !
    IF ( inlc > 0 .AND. inlc < 26 ) THEN
      !
      inlc_ = inlc
      IF ( nspin == 1 ) THEN
         CALL xc_vdW_DF      (rho_valence, rho_core, enl, vnl, v)
      ELSE IF ( nspin == 2 ) THEN
         CALL xc_vdW_DF_spin (rho_valence, rho_core, enl, vnl, v)
      ELSE
         CALL errore ('nlc', 'vdW-DF not available for noncollinear spin case',1)
      END If
      !
    ELSE IF ( inlc == 26 ) THEN
      !
      IF ( xclib_get_id('MGGA','EXCH') == 0 ) THEN
        CALL xc_rVV10 (rho_valence(:,1), rho_core, nspin, enl, vnl, v)
      ELSE
        CALL xc_rVV10 (rho_valence(:,1), rho_core, nspin, enl, vnl, v, 15.7_dp)
      END IF
      !
    ELSE
      !
      CALL errore ('nlc', 'inlc choice for E^nl_c not implemented',1)
      !
    END IF
    !
    RETURN
    !
  END SUBROUTINE nlc
  !
  !
END MODULE funct
