!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
MODULE xmltools
  !--------------------------------------------------------
  !
  ! Poor-man tools for reading and writing xml files
  !
  USE upf_kinds, ONLY : dp
  IMPLICIT NONE
  !
  ! internal variables for reading and writing
  !
  INTEGER :: xmlunit
  character(len=1024) :: line
  integer :: eot
  integer :: nattr
  CHARACTER(LEN=:), ALLOCATABLE :: attrlist, attrvals
  !
  ! variables used keep track of open tags
  !
  INTEGER :: nlevel = -1
  INTEGER, PARAMETER :: maxlength=80, maxlevel=9
  CHARACTER(LEN=maxlength), DIMENSION(0:maxlevel) :: open_tags
  !
  PRIVATE
  PUBLIC :: xml_openfile, xml_closefile
  PUBLIC :: xmlw_writetag, xmlw_opentag, xmlw_closetag
  PUBLIC :: xml_protect
  PUBLIC :: i2c, l2c, r2c
  PUBLIC :: xmlr_readtag, xmlr_opentag, xmlr_closetag
  PUBLIC :: get_attr
  !
  INTERFACE xmlr_readtag
     MODULE PROCEDURE readtag_c, readtag_r, readtag_l, readtag_i, readtag_rv
  END INTERFACE xmlr_readtag
  !
  INTERFACE xmlw_writetag
     MODULE PROCEDURE writetag_c, writetag_r, writetag_l, writetag_i, writetag_rv
  END INTERFACE xmlw_writetag
  !
  INTERFACE get_attr
     MODULE PROCEDURE get_i_attr, get_l_attr, get_r_attr, get_c_attr
  END INTERFACE get_attr
  
CONTAINS

  SUBROUTINE get_i_attr ( attrname, attrval_i )
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    INTEGER, INTENT(OUT) :: attrval_i
    !
    CHARACTER(LEN=80) :: attrval_c
    !
    CALL get_c_attr ( attrname, attrval_c )
    if ( len_trim(attrval_c) > 0 ) then
       READ (attrval_c,*) attrval_i
    else
       attrval_i = 0
    end if
    !
  END SUBROUTINE get_i_attr
  !
  SUBROUTINE get_l_attr ( attrname, attrval_l )
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    LOGICAL, INTENT(OUT) :: attrval_l
    !
    CHARACTER(LEN=80) :: attrval_c
    !
    CALL get_c_attr ( attrname, attrval_c )
    if ( len_trim(attrval_c) > 0 ) then
       READ (attrval_c,*) attrval_l
    else
       attrval_l = .false.
    end if
    !
  END SUBROUTINE get_l_attr
  !
  SUBROUTINE get_r_attr ( attrname, attrval_r )
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    REAL(dp), INTENT(OUT) :: attrval_r
    !
    CHARACTER(LEN=80) :: attrval_c
    !
    CALL get_c_attr ( attrname, attrval_c )
    if ( len_trim(attrval_c) > 0 ) then
       READ (attrval_c,*) attrval_r
    else
       attrval_r = 0.0_dp
    end if
    !
  END SUBROUTINE get_r_attr
  !
  SUBROUTINE get_c_attr ( attrname, attrval_c )
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    CHARACTER(LEN=*), INTENT(OUT) :: attrval_c
    !
    INTEGER :: i0, i1, i, j0, j1
    !
    ! search for attrname in attrlist
    !
    j0 = 1
    do while ( j0 < len_trim(attrlist) )
       j1 = index ( attrlist(j0:), ',' )
       if ( j1 == 0 ) then
          ! no more commas: check if found
          if ( attrname == attrlist(j0:) ) exit
          ! here if not found
          attrval_c = ' '
          return
       else
          ! check if found: need exact match between commas
          if ( trim(attrname) == attrlist(j0:j0+j1-2) ) exit
       end if
       j0 = j0+j1
    end do
    !
    ! for every comma found in attrlist before attrname is found
    ! (at position j0) skip to the field after a comma in attrvals
    !
    i0 = 1
    DO i = 1, j0
       if ( attrlist(i:i) == ',' ) then
          i1 = index ( attrvals(i0:), ',' )
          i0 = i1 + i0
       end if
    END DO
    ! now find next comma
    i1 = index ( attrvals(i0:), ',' )
    IF ( i1 == 0 ) THEN
       ! no more commas?
       attrval_c = attrvals(i0:)
    ELSE IF ( i1 == 1 ) THEN
       ! two commas, one after the other (,,) ?
       attrval_c = ' '
    ELSE
       ! take field between two commas
       attrval_c = attrvals(i0:i0+i1-2)
    END IF
    !
  END SUBROUTINE get_c_attr
  !
  FUNCTION xml_openfile ( filexml ) RESULT (iun)
    !
    ! returns on output the opened unit number if opened successfully
    ! returns -1 otherwise
    !
    CHARACTER(LEN=*), INTENT(in) :: filexml
    INTEGER :: iun, ios
    !
    OPEN ( NEWUNIT=iun, FILE=filexml, FORM='formatted', STATUS='unknown', &
         IOSTAT=ios)
    IF ( ios /= 0 ) iun = -1
    xmlunit = iun
    nlevel = 0
    open_tags(nlevel) = 'root'
    !
  END FUNCTION xml_openfile
  !
  SUBROUTINE xml_closefile ( )
    !
    CLOSE ( UNIT=xmlunit, STATUS='keep' )
    xmlunit = -1
    IF ( nlevel > 0 ) THEN
       print '("severe error: file closed at level ",i1," with tag ",A," open")', &
            nlevel, trim(open_tags(nlevel))
    END IF
    nlevel = 0
    !
  END SUBROUTINE xml_closefile
  !
  SUBROUTINE xmlw_opentag (name, attrlist, attrvals, ierr )
    ! On input:
    ! name      required, character: tag name
    ! attrlist  optional, character: list of comma-separated attributes
    ! attrvals  optional, character: list of comma-separated attrbute values
    ! On output: the tag is left open, ready for addition of data -
    !            the tag must be subsequently closed with close_xml_tag
    ! If ierr is present, the following value is returned:
    ! ierr = 0     normal execution
    ! ierr = 1     cannot write to unit "xmlunit"
    ! ierr = 2     tag name too long
    ! ierr = 3     too many tag levels
    ! ierr =10     wrong number of values for attributes
    ! If absent, the above error messages are printed.
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrlist
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrvals
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    INTEGER :: ier_
    CHARACTER(LEN=1) :: tag_end='>'
    !
    ier_ = write_tag_and_attr (name, attrlist, attrvals )
    IF ( ier_ < 0 ) ier_ = 0
    ! complete tag, leaving it open for further data
    WRITE (xmlunit, "(A1)", ERR=100) tag_end
    ! exit here
100 IF ( present(ierr) ) THEN
       ierr = ier_
    ELSE IF ( ier_ > 0 ) THEN
       print '("Fatal error ",i2," in xmlw_opentag!")', ier_
    END IF
    !
  END SUBROUTINE xmlw_opentag

  SUBROUTINE writetag_c (name, cval, attrlist, attrvals, ierr )
    ! On input, same as xmlw_opentag, plus:
    ! cval   character, value of the tag.
    ! If cval=' ' write <name  attr1="val1" attr2="val2" ... /> 
    ! If cval='?' write <?name attr1="val1" attr2="val2" ...?> 
    ! otherwise,  write <name  attr1="val1" attr2="val2" ...>cval</name> 
    !             (su di una stessa riga)
    ! On output, same as xmlw_opentag
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(IN) :: cval
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrlist
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrvals
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    INTEGER :: ier_
    LOGICAL :: is_proc
    !
    is_proc = (LEN_TRIM(cval) == 1 .AND. cval(1:1) == '?')
    IF (is_proc) THEN
       ier_ = write_tag_and_attr ( '?'//name, attrlist, attrvals)
    ELSE
       ier_ = write_tag_and_attr ( name, attrlist, attrvals)
    END IF
    IF ( ier_ > 0 ) GO TO 10
    !
    ! all is well: write tag value if any, close otherwise
    !
    IF ( LEN_TRIM(cval) == 0 ) THEN
       ! empty tag value: close here the tag
       CALL xmlw_closetag ( '' )
    ELSE IF ( is_proc ) THEN
       ! close "process" tag (e.g. <?xml ... ?>
       CALL xmlw_closetag ( '?' )
    ELSE
       ! write value (character)
       WRITE (xmlunit, "('>',A)", ADVANCE='no') trim(cval)
       ! close here the tag
       CALL xmlw_closetag ( name )
    END IF
    ! in case of exit error close the tag anyway
10  IF ( ier_ /= 0 ) WRITE (xmlunit, "('>')", ERR=100)
100 IF ( present(ierr) ) THEN
       ierr = ier_
    ELSE IF ( ier_ > 0 ) THEN
       print '("Fatal error ",i2," in xmlw_writetag!")', ier_
       stop
    END IF
    !
  END SUBROUTINE writetag_c
  !
  SUBROUTINE writetag_i (name, ival, attrlist, attrvals, ierr )
    !
    ! As writetag_c, for integer value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN)          :: ival
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrlist
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrvals
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL writetag_c (name, i2c(ival), attrlist, attrvals, ierr )
    !
  END SUBROUTINE writetag_i
  !
  SUBROUTINE writetag_l (name, lval, attrlist, attrvals, ierr )
    !
    ! As writetag_c, for logical value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN)          :: lval
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrlist
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrvals
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL writetag_c (name, l2c(lval), attrlist, attrvals, ierr )
    !
  END SUBROUTINE writetag_l
  !
  SUBROUTINE writetag_r (name, rval, attrlist, attrvals, ierr )
    !
    ! As writetag_c, for real value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(dp), INTENT(IN)         :: rval
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrlist
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrvals
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL writetag_c (name, r2c(rval), attrlist, attrvals, ierr )
    !
  END SUBROUTINE writetag_r
  !
  SUBROUTINE writetag_rv (name, rval, attrlist, attrvals, ierr )
    !
    ! As writetag_c, for an array of real values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(dp), INTENT(IN)         :: rval(:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrlist
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrvals
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL xmlw_opentag (name, attrlist, attrvals, ierr )
    WRITE( xmlunit, *) rval
    CALL xmlw_closetag ( )
    !
  END SUBROUTINE writetag_rv
  
  FUNCTION write_tag_and_attr (name, attrlist, attrvals) RESULT (ierr)
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrlist
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: attrvals
    INTEGER :: ierr
    !
    LOGICAL :: have_list, have_vals
    INTEGER :: i, la, lv, n1a,n2a, n1v, n2v
    !
    IF ( LEN_TRIM(name) > maxlength ) THEN
       ierr = 2
       RETURN
    END IF
    !
    IF ( nlevel+1 > maxlevel ) THEN
       ierr = 3
       RETURN
    END IF
    nlevel = nlevel+1
    open_tags(nlevel) = TRIM(name)
    !
    ! pretty (?) printing
    !
    ierr = 1
    DO i=2,nlevel
       WRITE (xmlunit, "('  ')", ADVANCE="no", ERR=10)
    END DO
    WRITE (xmlunit, "('<',A)", ADVANCE="no", ERR=10) trim(name)
    print '("opened at level ",i1," tag ",A)', nlevel, trim(open_tags(nlevel))
    !
    ! attributes (if present)
    !
    ierr = 10
    have_list = PRESENT(attrlist)
    have_vals = PRESENT(attrvals)
    ! return with error code
    IF ( ( have_list .AND. .NOT. have_vals ) .OR. &
         ( have_vals .AND. .NOT. have_list ) ) RETURN
    IF ( have_list .AND. have_vals ) THEN
       !
       la=len_trim(attrlist)
       lv=len_trim(attrvals)
       ! skip initial white spaces. Alternatively:
       ! n1a=1; do while (attrlist(n1a:n1a) == ' '); n1a=n1a+1; end do
       do n1a = 1, la
          if ( attrlist(n1a:n1a) /= ' ' ) exit
       end do
       do n1v = 1, lv
          if ( attrvals(n1v:n1v) /= ' ' ) exit
       end do
       do while ( (n1a <= la) .AND. (n1v <= lv) )
          !
          ! comma is the separator
          !
          n2a = INDEX( attrlist(n1a:), ',' )
          n2v = INDEX( attrvals(n1v:), ',' )
          ! mismatch between the number of attributes and of values:
          ! return with error code
          IF ( ( n2a == 0 .and. n2v /= 0 ) .or. &
               ( n2a /= 0 .and. n2v == 0 ) ) RETURN
          !
          IF ( ( n2a == 0 .and. n2v == 0 ) ) THEN
             ! last attribute and respective value
             WRITE (xmlunit, "(' ',A,'=""',A,'""')", ADVANCE='no', ERR=10) &
                  attrlist(n1a:la), attrvals(n1v:lv)
             EXIT
          ELSE
             WRITE (xmlunit, "(' ',A,'=""',A,'""')", ADVANCE='no', ERR=10) &
                  attrlist(n1a:n1a+n2a-2), attrvals(n1v:n1v+n2v-2)
             n1a = n1a+n2a
             n1v = n1v+n2v
          END IF
          !
       END DO
       !
    END IF
    ! normal exit here
    ierr = 0
10  RETURN
    !
  END FUNCTION write_tag_and_attr
  !
  SUBROUTINE xmlw_closetag ( tag )
    ! tag   not present: close current open tag with </tag>
    ! empty tag present: close current open tag with />
    ! tag='?'   present: close current open tag with ?>
    ! otherwise,close specified tag with </tag>
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
    INTEGER :: i
    !
    IF ( nlevel < 0 ) &
         print '("severe error: closing tag that was never opened")'
    IF ( .NOT.PRESENT(tag) ) THEN
       DO i=2,nlevel
          WRITE (xmlunit, '("  ")', ADVANCE='NO')
       END DO
       WRITE (xmlunit, '("</",A,">")') trim(open_tags(nlevel))
    ELSE
       i = len_trim(tag)
       IF ( i == 0 ) THEN
          WRITE (xmlunit, '("/>")')
       ELSE IF ( i == 1 .AND. tag(1:1) == '?' ) THEN
          WRITE (xmlunit, '("?>")')
       ELSE
          WRITE (xmlunit, '("</",A,">")') trim(tag)
       END IF
    END IF
    print '("closed at level ",i1," tag ",A)', nlevel, trim(open_tags(nlevel))
    nlevel = nlevel-1
  END SUBROUTINE xmlw_closetag
  !
  !--------------------------------------------------------
  function xml_protect ( data_in ) result (data_out)
    !--------------------------------------------------------
    ! 
    ! poor-man escaping of a string so that it conforms to xml standard:
    ! replace & with @, < and > with *. data_out is also left-justified
    !
    character(len=*), intent(in) :: data_in
    character(len=:), allocatable :: data_out
    character(len=1) :: c
    integer:: n, i
    !
    n = len_trim(adjustl(data_in)) 
    ! Alternative version with CDATA:
    ! allocate(character(len=n+12):: data_out)
    ! data_out = '<![CDATA['//trim(adjustl(data_in))//']]>'
    data_out = trim(adjustl(data_in))
    do i=1,n
       if ( data_out(i:i) == '&' ) data_out(i:i) = '@'
       if ( data_out(i:i) == '<' .or. data_out(i:i) == '>') data_out(i:i) = '*'
    end do
    ! a more complete version should escape & as &amp;, < as &lt;
    ! (escaping > as &gt; , " as &quotes; , ' as &apo; is not strictly needed) 
    ! BUT taking care not to escape &amp; into &amp;amp;

  end function xml_protect
  
  ! Poor-man conversion utilities from integer, logical, real to character
  ! To be used in conjunction with routines in module xmlw to write xml
  !
  function i2c (i) result (c)
    integer, intent(in) :: i
    character(len=:), allocatable :: c
    character(len=11) :: caux
    !
    write(caux,'(i11)') i
    c = trim(adjustl(caux))
    !
  end function i2c
  
  function l2c (l) result (c)
    logical, intent(in) :: l
    character(len=:), allocatable :: c
    !
    if (l) then
       c='true'
    else
       c='false'
    endif
    !
  end function l2c
  
  function r2c (f) result (c)
    real(dp), intent(in) :: f
    character(len=:), allocatable :: c
    character(len=30) :: caux
    !
    integer :: n, m, i
    ! The format of real numbers can be vastly improved
    ! this is just the simplest solution
    write(caux,*) f
    c = trim(adjustl(caux))
    !
  end function r2c
  !
  SUBROUTINE readtag_i (name, ival, ierr )
    !
    ! As readtag_c, for integer value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT)         :: ival
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    CHARACTER(LEN=80) :: cval
    !
    CALL readtag_c (name, cval, ierr )
    if ( len_trim(cval) > 0 ) then
       READ (cval,*) ival
    else
       ival = 0
    end if
    !
  END SUBROUTINE readtag_i
  !
  SUBROUTINE readtag_l (name, lval, ierr )
    !
    ! As readtag_c, for logical value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(OUT)         :: lval
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    CHARACTER(LEN=80) :: cval
    !
    CALL readtag_c (name, cval, ierr )
    if ( len_trim(cval) > 0 ) then
       READ (cval,*) lval
    else
       lval = .false.
    end if
    !
  END SUBROUTINE readtag_l
  !
  SUBROUTINE readtag_r (name, rval, ierr )
    !
    ! As readtag_c, for real value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(8), INTENT(OUT)         :: rval
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    CHARACTER(LEN=80) :: cval
    !
    CALL readtag_c (name, cval, ierr )
    if ( len_trim(cval) > 0 ) then
       READ (cval,*) rval
    else
       rval = 0.0_dp
    end if
    !
  END SUBROUTINE readtag_r
  !
  SUBROUTINE readtag_rv (name, rval, ierr)
    !
    ! As readtag_c, for an array of real values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(8), INTENT(OUT)         :: rval(:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    INTEGER :: ier_
    CHARACTER(LEN=80) :: cval    
    !
    CALL xmlr_opentag (name, ier_)
    if ( ier_ == 0  ) then
       READ(xmlunit, *) rval
    else
       rval = 0.0_dp
    end if
    CALL xmlr_closetag ( )
    IF ( present (ierr) ) ierr = ier_
    !
  END SUBROUTINE readtag_rv
  !
  subroutine readtag_c ( tag, cval, ierr)
    !
    implicit none
    !
    character(len=*), intent(in) :: tag
    character(len=*), intent(out):: cval
    integer, intent(out), optional :: ierr
    ! 0: tag found and read
    !-1: tag not found
    ! 1: error parsing file
    ! 2: error in arguments
    !
    integer ::  i, j, lt
    character(len=1) :: endtag
    !
    call xmlr_opentag ( tag, ierr )
    !
    if ( eot > 0 ) then
       j = eot
       lt = len_trim(tag)
       ! beginning of val at line(j:j): search for end tag
       i = index ( line(j:), '</'//trim(tag) )
       if ( i < 1 ) then
          ! </tag> not found on this line
          ! print *, 'tag </',trim(tag),'> not found'
          ierr = 1
          return
       else
          ! maybe found end tag?
          endtag = adjustl( line(j+i+1+lt:) )
          if ( endtag /= '>' ) then
             ! print *, 'tag ',trim(tag),' not correctly closed'
             if (present(ierr)) ierr = 1
          else
             ! <tag ....>val</tag> found, exit
             cval = adjustl(trim(line(j:j+i-2)))
             ! print *, 'value=',cval
          end if
          nlevel = nlevel -1
          !
          return
          !
       endif
    else if ( eot == 0 ) then
       ! print *, 'end of file reached, tag not found'
       if ( present(ierr) ) ierr =-1
    else if ( eot < 0 ) then
       ! print *, 'tag found, no value to read on line'
       cval = ''
    end if
    !
  end subroutine readtag_c
  !
  subroutine xmlr_opentag ( tag, ierr)
    !
    implicit none
    !
    character(len=*), intent(in) :: tag
    integer, intent(out), optional :: ierr
    ! 0: tag found and read
    !-1: tag not found
    ! 1: error parsing file
    ! 2: line too long
    ! 3: too many levels of tags
    !
    integer :: stat, ll, lt, i, j, j0
    ! stat= 0: begin
    ! stat=-1: in comment
    ! stat=1 : tag found
    !
    character(len=1) :: quote
    !
    nattr =0
    if ( allocated (attrlist) ) deallocate(attrlist)
    if ( allocated (attrvals) ) deallocate(attrvals)
    attrlist=''
    attrvals=''
    !
    lt = len_trim(tag)
    stat=0
    eot =0
    do while (.true.)
       read(xmlunit,'(a)', end=10) line
       ll = len_trim(line)
       if ( ll == 1024 ) then
          print *, 'line too long'
          if (present(ierr)) ierr = 2
          return
       end if
       ! j is the current scan position
       j = 1
       parse: do while ( j <= ll )
          !
          if ( stat ==-1 ) then
             !
             ! scanning a comment
             i = index(line(j:),'-->')
             if ( i == 0 ) then
                ! no end of comment found on this line
                exit parse
             else
                ! end of comment found
                stat = 0
                j = j+i+3
             end if
             !
          else if ( stat == 0 ) then
             !
             ! searching for tag
             !
             i = index( line(j:),'<'//trim(tag) )
             if ( i == 0 ) then
                ! no tag found on this line
                exit parse
             else
                ! tag found? check what follows our would-be tag
                j = j+i+lt
                if ( j > ll ) then
                   print *, 'oops... opened tag not closed on same line'
                   exit parse
                else if ( line(j:j) == ' ' .or. line(j:j) == '>') then
                   ! print *, '<tag found'
                   stat = 1
                end if
             end if
             !
          else if ( stat == 1 ) then
             ! tag found, search for attributes if any or end of tag
             if (line(j:j) == ' ' ) then
                ! skip blanks: there is at least one if attributes are present
                j = j+1
                ! save value of j into j0: beginning of an attribute 
                j0= j
             else if ( line(j:j+1) == '/>' ) then
                ! <tag ... /> found : return
                if (present(ierr)) ierr = 0
                ! eot = -2: tag with no value found
                eot = -2
                !
                return
                !
             else if ( line(j:j) == '>' ) then
                ! <tag ... > found
                if ( j+1 > ll ) then
                   ! eot = -1: tag found, line ends
                   eot = -1
                else
                   ! eot points to the rest of the line
                   eot = j+1
                end if
                if (present(ierr)) ierr = 0
                nlevel = nlevel+1
                IF ( nlevel > maxlevel ) THEN
                   print *, ' too many levels'
                   if (present(ierr)) ierr = 3
                else
                   open_tags(nlevel) = trim(tag)
                end if
                !
                return
                !
             else if ( line(j:j) == '=' ) then
                ! end of attribute located: save attribute
                nattr=nattr+1
                ! print *, 'attr=',line(j0:j-1)
                if ( nattr == 1 ) then
                   attrlist = line(j0:j-1)
                else
                   attrlist = attrlist//','//line(j0:j-1)
                end if
                ! continue searching for attribute value
                j = j+1
             else if ( line(j:j) == '"' .or. line(j:j) =="'" ) then
                ! first occurrence of ' or " found, look for next
                quote = line(j:j)
                j=j+1
                i = index(line(j:),quote)
                if ( i < 1 ) then
                   print *, 'Error: matching quote not found'
                else
                   ! save attribute value and continue scanning
                   ! print *, 'attrval=',line(j:j+i-2)
                   if ( nattr == 1 ) then
                      attrvals = line(j:j+i-2)
                   else
                      attrvals = attrvals//','//line(j:j+i-2)
                   end if
                   j = j+i
                end if
             else
                ! continue scanning until end of attribute
                j = j+1
             endif
             !
          end if
       end do parse
       !
    end do
    !
10  if ( stat == 0 ) then
       if ( present(ierr) ) then
          ierr =-1
       else
          print *, 'end of file reached, tag '//trim(tag)//' not found'
       end if
    else
       print *, 'parsing error'
       if ( present(ierr) ) ierr = 1
    end if
    !
  end subroutine xmlr_opentag
  !
  subroutine xmlr_closetag ( tag, ierr)
    !
    implicit none
    !
    character(len=*), intent(in), optional :: tag
    integer, intent(out), optional :: ierr
    ! 0: </tag> found
    ! 1: </tag> not found
    ! 2: error parsing file
    !
    integer :: stat, ll, lt, i, j
    ! stat=-1: in comment
    ! stat= 0: begin
    ! stat= 1: end
    !
    IF ( nlevel < 0 ) &
         print '("severe error: closing tag that was never opened")'
    stat=0
    do while (.true.)
       read(xmlunit,'(a)', end=10) line
       ll = len_trim(line)
       if ( ll == 1024 ) then
          print *, 'line too long'
          if (present(ierr)) ierr = 1
          return
       end if
       ! j is the current scan position
       j = 1
       parse: do while ( j <= ll )
          !
          if ( stat ==-1 ) then
             !
             ! scanning a comment
             i = index(line(j:),'-->')
             if ( i == 0 ) then
                ! no end of comment found on this line
                exit parse
             else
                ! end of comment found
                stat = 0
                j = j+i+3
             end if
             !
          else if ( stat == 0 ) then
             !
             ! searching for closing tag
             !
             IF ( .NOT.PRESENT(tag) ) THEN
                i = index( line(j:),'</'//trim(open_tags(nlevel)) )
                lt= len_trim(open_tags(nlevel))
             ELSE
                i = index( line(j:),'</'//trim(tag) )
                lt= len_trim(tag)
             END IF
             if ( i == 0 ) then
                ! no tag found on this line
                exit parse
             else
                ! tag found? check what follows our would-be tag
                j = j+i+1+lt
                if ( j > ll ) then
                   print *, 'oops... opened tag not closed on same line'
                   exit parse
                else if ( line(j:j) == ' ' .or. line(j:j) == '>') then
                   ! print *, '</tag found'
                   stat = 1
                end if
             end if
             !
          else if ( stat == 1 ) then
             !
             ! </tag found, search for end of tag
             !
             if (line(j:j) == ' ' ) then
                ! skip blanks
                j = j+1
             else if ( line(j:j) == '>' ) then
                ! </tag ... > found
                ! print *, '</tag> found'
                if ( present(ierr) ) ierr = 0
                nlevel = nlevel - 1
                !
                return
                !
             endif
             !
          end if
       end do parse
       !
    end do
    !
10  print *, 'end of file reached, closing tag not found'
    if ( present(ierr) ) ierr = 1
    !
  end subroutine xmlr_closetag

END MODULE xmltools
