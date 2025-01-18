!
! Copyright (C) 2020-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
MODULE xmltools
  !--------------------------------------------------------
  !
  ! Poor-man set of tools for reading and writing xml files
  ! Similar to iotk but much simpler - Paolo Giannozzi, June 2020
  ! Limitations: too many to be listed in detail. Main ones:
  ! * works on a single opened file at the time. Exception:
  !   while a file is opened, one can open, R/W, close another file,
  !   but it is not possible to operate on both files at the same time
  ! * lines no more than 1024 characters long (see maxline parameter)
  ! * no more than 9 levels of tags (see maxlevel parameter)
  ! * length of tags no more than 80 characters (see maxlength parameter)
  ! * can read tags only in the correct order. If a tag is not found, the
  !   file is rewound. If "ierr" is present, a second attempt to find the
  !   tag is done starting from the top of the file - may work if the searched
  !   tag is found only above the current position, and nowhere else
  ! * only single values (e.g. no vectors) in attributes
  ! * attributes should not contain commas or strange characters
  ! * xml comments (<!-- ...  -->) or <![CDATA[ ... ]]> cannot be mixed
  !   with numerical fields
  !
#if defined(__XML_STANDALONE)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP_XML = selected_real_kind(14,200)
  PUBLIC :: DP_XML
#else
  USE upf_kinds, ONLY : DP_XML => dp
  IMPLICIT NONE
#endif
  !
#undef __debug
  !! define __debug to print information on opened and closed tags
  LOGICAL, PARAMETER :: one_line_tags=.true.
  !! if true, write tags with one value in a single line: 
  !! <tag attr1="val1" ... >value</tag>
  !! otherwise, as in iotk:
  !! <tag attr1="val1" ... >
  !! value
  !! </tag>
  !! Only for single values; arrays are always written as in iotk
  !
  ! internal variables for reading and writing
  !
  INTEGER :: xmlunit = -1
  INTEGER, PARAMETER :: maxline=1024, maxdim=maxline+16
  CHARACTER(LEN=maxdim) :: line
  INTEGER :: xmlsave = -1, nopen = 0
  INTEGER :: eot
  ! eot points to the end of tag in line just scanned
  INTEGER :: nattr
  CHARACTER(LEN=:), ALLOCATABLE :: attrlist
  !
  ! variables used to keep track of open tags
  !
  INTEGER :: nlevel = -1, level0 = 0
  INTEGER, PARAMETER :: maxlength=80, maxlevel=9
  CHARACTER(LEN=maxlength), DIMENSION(0:maxlevel) :: open_tags
  !
  PRIVATE
  ! general subroutines
  PUBLIC :: xml_open_file, xml_closefile
  ! subroutines for writing
  PUBLIC :: add_attr
  PUBLIC :: xmlw_writetag, xmlw_opentag, xmlw_closetag
  ! subroutines for reading
  PUBLIC :: xmlr_readtag, xmlr_opentag, xmlr_closetag
  PUBLIC :: get_attr
  ! utility functions
  PUBLIC :: xml_protect, i2c, l2c, r2c
  !
  ! Error codes returned by xmlr_opentag / xml_readtag:
  !  -11  as -10 + -1
  !  -10  tag found and read above the current position, after a file rewind:
  !       may or may not what you desire
  !  -1   tag with no value (e.g. <tag attr="val"/>) found (no error)
  !   0   tag found and read (no error)
  !   1   tag not found
  !   2   error parsing file
  !   3   line too long
  !   4   too many levels of tags
  ! 
  ! Error codes returned by xmlw_opentag / xml_writetag:
  !   0     tag open and/or written (no error)
  !   1     cannot write to unit "xmlunit"
  !   2     tag name too long
  !   3     wrong number of values for attributes
  !   4     too many levels of tag
  ! 
  INTERFACE xmlr_readtag
     MODULE PROCEDURE readtag_c, readtag_r, readtag_l, readtag_i, &
          readtag_iv, readtag_rv, readtag_rm, readtag_rt, &
          readtag_zv, readtag_zm, readtag_zt
  END INTERFACE xmlr_readtag
  !
  ! IMPORTANT NOTICE: complex numbers, z=a+ib, are written as two reals:
  ! "a b", not in fortran free format as "(a,b)". Reason:
  ! make the file readable by non-fortran tools, e.g. python
  !
  INTERFACE xmlw_writetag
     MODULE PROCEDURE writetag_c, writetag_r, writetag_l, writetag_i, &
          writetag_iv, writetag_rv, writetag_rm, writetag_rt, &
          writetag_zv, writetag_zm, writetag_zt
  END INTERFACE xmlw_writetag
  !
  INTERFACE get_attr
     MODULE PROCEDURE get_i_attr, get_l_attr, get_r_attr, get_c_attr
  END INTERFACE get_attr

  INTERFACE add_attr
     MODULE PROCEDURE add_i_attr, add_l_attr, add_r_attr, add_c_attr
  END INTERFACE add_attr
  
CONTAINS

  SUBROUTINE get_i_attr ( attrname, attrval_i )
    !
    ! returns attrval_i=0 if not found or not readable
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    INTEGER, INTENT(OUT) :: attrval_i
    !
    CHARACTER(LEN=80) :: attrval_c
    !
    CALL get_c_attr ( attrname, attrval_c )
    if ( len_trim(attrval_c) > 0 ) then
       READ (attrval_c,*, err=1) attrval_i
       return
1      print '("Error reading attribute ",a,": expected integer, found ",a)', &
            trim(attrname), trim(attrval_c)
    end if
    attrval_i = 0
    !
  END SUBROUTINE get_i_attr
  !
  SUBROUTINE get_l_attr ( attrname, attrval_l )
    !
    ! returns attrval_l=.false. if not found or not readable
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    LOGICAL, INTENT(OUT) :: attrval_l
    !
    CHARACTER(LEN=80) :: attrval_c
    !
    CALL get_c_attr ( attrname, attrval_c )
    if ( len_trim(attrval_c) > 0 ) then
       READ (attrval_c,*, err=1) attrval_l
       return
1      print '("Error reading attribute ",a,": expected logical, found ",a)', &
            trim(attrname), trim(attrval_c)
    end if
    attrval_l = .false.
    !
  END SUBROUTINE get_l_attr
  !
  SUBROUTINE get_r_attr ( attrname, attrval_r )
    !
    ! returns attrval_r=0 if not found or not readable
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    REAL(DP_XML), INTENT(OUT) :: attrval_r
    !
    CHARACTER(LEN=80) :: attrval_c
    !
    CALL get_c_attr ( attrname, attrval_c )
    if ( len_trim(attrval_c) > 0 ) then
       READ (attrval_c,*, err=1) attrval_r
       return
1      print '("Error reading attribute ",a,": expected real, found ",a)', &
            trim(attrname), trim(attrval_c)
    end if
    attrval_r = 0.0_DP_XML
    !
  END SUBROUTINE get_r_attr
  !
  SUBROUTINE get_c_attr ( attrname, attrval_c )
    !
    ! returns attrval_c='' if not found
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    CHARACTER(LEN=*), INTENT(OUT) :: attrval_c
    !
    CHARACTER(LEN=1) :: quote
    INTEGER :: j0, j1
    LOGICAL :: found
    !
    ! search for attribute name in attrlist: attr1="val1" attr2="val2" ...
    !
    attrval_c = ''
    if ( .not. allocated(attrlist) ) return
    if ( len_trim(attrlist) < 1 ) return
    !
    j0 = 1
    do while ( j0 < len_trim(attrlist) )
       ! locate = and first quote
       j1 = index ( attrlist(j0:), '=' )
       quote = attrlist(j0+j1:j0+j1)
       ! next line: something is not right
       if ( quote /= '"' .and. quote /= "'" ) return
       ! check if attribute found: need exact match 
       found =  ( trim(attrname) == adjustl(trim(attrlist(j0:j0+j1-2))) )
       ! locate next quote
       j0 = j0+j1+1
       j1 = index ( attrlist(j0:), quote )
       if ( found) then
          if ( j1 == 1 ) then
             ! two quotes, one after the other ("")
             attrval_c = ' '
          else
             ! get value between two quotes
             attrval_c = adjustl(trim(attrlist(j0:j0+j1-2)))
          end if
          return
       end if
       j0 = j0+j1
    end do
    !
  END SUBROUTINE get_c_attr
  !
  SUBROUTINE add_i_attr ( attrname, attrval_i )
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    INTEGER, INTENT(IN) :: attrval_i
    !
    CALL add_c_attr ( attrname, i2c(attrval_i) )
    !
  END SUBROUTINE add_i_attr
  !
  SUBROUTINE add_l_attr ( attrname, attrval_l )
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    LOGICAL, INTENT(IN) :: attrval_l
    !
    CALL add_c_attr ( attrname, l2c(attrval_l) )
    !
  END SUBROUTINE add_l_attr
  !
  SUBROUTINE add_r_attr ( attrname, attrval_r )
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    REAL(DP_XML), INTENT(IN) :: attrval_r
    !
    CALL add_c_attr ( attrname, r2c(attrval_r) )
    !
  END SUBROUTINE add_r_attr
  !
  SUBROUTINE add_c_attr ( attrname, attrval_c )
    !
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: attrname, attrval_c
    !
    IF ( .NOT. ALLOCATED(attrlist) ) THEN
       attrlist = ' '//TRIM(attrname)//'="'//TRIM(attrval_c)//'"'
    ELSE
       attrlist = attrlist // ' ' // TRIM(attrname)//'="'//TRIM(attrval_c)//'"'
    END IF
    !
  END SUBROUTINE add_c_attr
  !
  FUNCTION xml_open_file ( filexml ) RESULT (iun)
    !
    ! returns on output the opened unit number if opened successfully
    ! returns -1 otherwise
    !
    CHARACTER(LEN=*), INTENT(in) :: filexml
    INTEGER :: iun, ios
    !
    IF ( nopen > 1 ) THEN
       print "('cannot open file ',a,': two xml files already opened')",&
               trim(filexml)
       iun = -1
       RETURN
    END IF
    OPEN ( NEWUNIT=iun, FILE=filexml, FORM='formatted', STATUS='unknown', &
         IOSTAT=ios)
    IF ( ios /= 0 ) THEN
       CLOSE ( UNIT=iun, STATUS='keep')
       iun = -1
#if defined ( __debug )
       print "('file ',a,' not opened, ios = ',i5)",trim(filexml),ios
#endif
    ELSE
       nopen = nopen + 1
       IF ( nopen > 1 ) THEN
       ! a second file is opened: keep track of the status of the previous file
          xmlsave = xmlunit
          level0  = nlevel
       ELSE
          nlevel = 0
          open_tags(nlevel) = 'root'
       END IF
#if defined ( __debug )
       print "('file ',a,' opened with unit ',i5)",trim(filexml),iun
#endif
    END IF
    xmlunit = iun
    if ( allocated(attrlist) ) DEALLOCATE ( attrlist)
    !
  END FUNCTION xml_open_file
  !
  SUBROUTINE xml_closefile ( )
    !
    integer :: ios
    !
    IF ( xmlunit == -1 ) RETURN
    CLOSE ( UNIT=xmlunit, STATUS='keep' )
#if defined ( __debug )
    print "('unit ',i5,': file closed')", xmlunit
#endif
    xmlunit = xmlsave
    nopen = nopen - 1
    xmlsave = -1
    IF (nlevel > level0 ) print '("warning: file closed at level ",i1,&
            & " with tag ",A," open")', nlevel, trim(open_tags(nlevel))
    IF ( nopen == 1 ) THEN
       ! the second file is opened: return to the status of previous file
       nlevel  = level0
    ELSE
       level0 = 0
    END IF
    !
  END SUBROUTINE xml_closefile
  !
  SUBROUTINE xmlw_opentag (name, ierr, noadv )
    ! On input:
    ! name      required, character: tag name
    ! On output: the tag is left open, ready for addition of data -
    !            the tag must be subsequently closed with close_xml_tag
    ! If ierr is present, the error code set in write_tag_and_attr is returned
    ! If ierr is absent,  the above error code is reprinted on output
    ! If noadv is present and true, stay on the same line
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    LOGICAL, INTENT(IN), OPTIONAL :: noadv
    !
    INTEGER :: ier_
    CHARACTER(LEN=1) :: tag_end='>'
    LOGICAL :: noadv_
    !
    ier_ = write_tag_and_attr (name)
    IF ( ier_ < 0 ) ier_ = 0
    ! complete tag, leaving it open for further data
    noadv_ = present(noadv)
    IF ( noadv_ ) noadv_ = noadv
    IF ( noadv_ ) THEN
       WRITE (xmlunit, "(A1)", ADVANCE="no",ERR=100) tag_end
    ELSE
       WRITE (xmlunit, "(A1)", ERR=100) tag_end
    END IF
    ! exit here
100 IF ( present(ierr) ) THEN
       ierr = ier_
    ELSE IF ( ier_ > 0 ) THEN
       print '("Fatal error ",i2," in xmlw_opentag!")', ier_
    END IF
    !
  END SUBROUTINE xmlw_opentag

  SUBROUTINE writetag_c (name, cval, ierr )
    ! On input, same as xmlw_opentag, plus:
    ! cval   character, value of the tag.
    ! If cval=' ' write <name  attr1="val1" attr2="val2" ... /> 
    ! If cval='?' write <?name attr1="val1" attr2="val2" ...?> 
    ! otherwise,  write <name  attr1="val1" attr2="val2" ...>cval</name> 
    !             (on a same line if one_line_tags=.true.)
    ! On output, same as xmlw_opentag
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(IN) :: cval
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    INTEGER :: ier_
    LOGICAL :: is_proc
    !
    is_proc = (LEN_TRIM(cval) == 1)
    IF ( is_proc ) is_proc = is_proc .AND. ( cval(1:1) == '?')
    IF (is_proc) THEN
       ier_ = write_tag_and_attr ( '?'//name )
    ELSE
       ier_ = write_tag_and_attr ( name )
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
       IF (one_line_tags) THEN
          WRITE (xmlunit, "('>',A)", ADVANCE='no') trim(cval)
       ELSE
          WRITE (xmlunit, "('>',/,A)") trim(cval)
       ENDIF
       ! close here the tag
       CALL xmlw_closetag ( name )
    END IF
    ! in case of exit error close the tag anyway
10  IF ( ier_ /= 0 ) WRITE (xmlunit, "('>')", ERR=100)
100 IF ( present(ierr) ) THEN
       ierr = ier_
    ELSE IF ( ier_ > 0 ) THEN
       print '("Fatal error ",i2," in xmlw_writetag!")', ier_
    END IF
    !
  END SUBROUTINE writetag_c
  !
  SUBROUTINE writetag_i (name, ival, ierr )
    !
    ! As writetag_c, for integer value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN)          :: ival
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL writetag_c (name, i2c(ival), ierr )
    !
  END SUBROUTINE writetag_i
  !
  SUBROUTINE writetag_iv (name, ivec, ierr )
    !
    ! As writetag_c, for integer value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN)          :: ivec(:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL xmlw_opentag (name, ierr )
    WRITE( xmlunit, '(4I18)') ivec
    CALL xmlw_closetag ( )
    !
  END SUBROUTINE writetag_iv
  !
  SUBROUTINE writetag_l (name, lval, ierr )
    !
    ! As writetag_c, for logical value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN)          :: lval
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL writetag_c (name, l2c(lval), ierr )
    !
  END SUBROUTINE writetag_l
  !
  SUBROUTINE writetag_r (name, rval, ierr )
    !
    ! As writetag_c, for real value
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(DP_XML), INTENT(IN)         :: rval
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL writetag_c (name, r2c(rval), ierr )
    !
  END SUBROUTINE writetag_r
  !
  SUBROUTINE writetag_rv (name, rvec, ierr )
    !
    ! As writetag_c, for a vector of real values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(DP_XML), INTENT(IN)         :: rvec(:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL xmlw_opentag (name, ierr )
    WRITE( xmlunit, '(1p3es24.15)') rvec
    CALL xmlw_closetag ( )
    !
  END SUBROUTINE writetag_rv
  !
  SUBROUTINE writetag_rm (name, rmat, ierr )
    !
    ! As writetag_c, for a matrix of real values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(DP_XML), INTENT(IN)         :: rmat(:,:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL xmlw_opentag (name, ierr )
    WRITE( xmlunit, '(1p3es24.15)') rmat
    CALL xmlw_closetag ( )
    !
  END SUBROUTINE writetag_rm
  !
  SUBROUTINE writetag_rt (name, rtens, ierr )
    !
    ! As writetag_c, for a 3-dim tensor of real values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(DP_XML), INTENT(IN)         :: rtens(:,:,:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    CALL xmlw_opentag (name, ierr )
    WRITE( xmlunit, '(1p3es24.15)') rtens
    CALL xmlw_closetag ( )
    !
  END SUBROUTINE writetag_rt
  !
  SUBROUTINE writetag_zv (name, zvec, ierr )
    !
    ! As writetag_c, for a vector of complex values
    !
    USE iso_c_binding
    CHARACTER(LEN=*), INTENT(IN) :: name
    COMPLEX(DP_XML), INTENT(IN), TARGET:: zvec(:)
    INTEGER, INTENT(OUT), OPTIONAL :: ierr
    !
    ! Casts a real pointer (rvec) to a complex array (zvec) via C pointer (!)
    ! in order to write complexes as two reals. Some compilers require that
    ! the argument of c_loc (zvec) is a pointer or has the "target" attribute
    !
    TYPE (c_ptr) :: cp
    REAL(DP_XML), POINTER  :: rvec(:)
    INTEGER :: n, ndim
    !
    NULLIFY (rvec)
    cp = c_loc(zvec)
    CALL c_f_pointer (cp, rvec, shape(zvec)*[2])
    CALL xmlw_opentag (name, ierr )
    ndim = SIZE (zvec)
    DO n=1,2*ndim,2
       WRITE( xmlunit, *) rvec(n), rvec(n+1)
    END DO
    CALL xmlw_closetag ( )
    !
  END SUBROUTINE writetag_zv
  !
  SUBROUTINE writetag_zm (name, zmat, ierr )
    !
    ! As writetag_c for a matrix of complex values - see comments in writetag_zv
    !
    USE iso_c_binding
    CHARACTER(LEN=*), INTENT(IN) :: name
    COMPLEX(DP_XML), INTENT(IN), TARGET:: zmat(:,:)
    INTEGER, INTENT(OUT), OPTIONAL :: ierr
    !
    TYPE (c_ptr) :: cp
    REAL(DP_XML), POINTER  :: rmat(:,:)
    !
    NULLIFY (rmat)
    cp = c_loc(zmat)
    CALL c_f_pointer (cp, rmat, shape(zmat)*[2,1])
    !
    CALL xmlw_opentag (name, ierr )
    WRITE( xmlunit, '(1p2es24.15)') rmat
    CALL xmlw_closetag ( )
    !
  END SUBROUTINE writetag_zm
  !
  SUBROUTINE writetag_zt (name, ztens, ierr )
    !
    ! As writetag_c for a matrix of complex values - see comments in writetag_zv
    !
    USE iso_c_binding
    CHARACTER(LEN=*), INTENT(IN) :: name
    COMPLEX(DP_XML), INTENT(IN), TARGET:: ztens(:,:,:)
    INTEGER, INTENT(OUT), OPTIONAL :: ierr
    !
    TYPE (c_ptr) :: cp
    REAL(DP_XML), POINTER  :: rtens(:,:,:)
    !
    NULLIFY (rtens)
    cp = c_loc(ztens)
    CALL c_f_pointer (cp, rtens, shape(ztens)*[2,1,1])
    !
    CALL xmlw_opentag (name, ierr )
    WRITE( xmlunit, '(2es24.15)') rtens
    CALL xmlw_closetag ( )
    !
  END SUBROUTINE writetag_zt
  !
  FUNCTION write_tag_and_attr (name) RESULT (ierr)
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: ierr
    ! See list of error codes in the header of this file
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
       ierr = 4
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
#if defined ( __debug )
    print '("opened (write) level-",i1," tag ",A)', nlevel, trim(open_tags(nlevel))
#endif
    !
    ! attributes (if present)
    !
    ierr = 3
    if ( allocated (attrlist) ) then
       WRITE (xmlunit, "(A)", ADVANCE='no', ERR=10) attrlist
       deallocate (attrlist)
    end if
    ! normal exit here
    ierr = 0
10  RETURN
    !
  END FUNCTION write_tag_and_attr
  !
  SUBROUTINE xmlw_closetag ( tag, noind )
    ! tag   not present: close current open tag with </tag>
    ! empty tag present: close current open tag with />
    ! tag='?'   present: close current open tag with ?>
    ! otherwise,close specified tag with </tag>
    ! If noind is present and true, do not indent
    !
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
    LOGICAL, INTENT(IN), OPTIONAL :: noind
    INTEGER :: i
    LOGICAL :: indent
    !
    IF ( nlevel < 1 ) THEN
       IF ( nlevel < 0 ) print "('xmlw_closetag: severe error, closing tag that was never opened')"
      RETURN
    END IF
    IF ( .NOT.PRESENT(tag) ) THEN
       indent = .NOT. PRESENT(noind)
       IF ( .NOT. indent ) indent = .NOT.noind
       IF ( indent ) THEN
          DO i=2,nlevel
             WRITE (xmlunit, '("  ")', ADVANCE='NO')
          END DO
       END IF
       WRITE (xmlunit, '("</",A,">")') trim(open_tags(nlevel))
#if defined ( __debug )
    print '("closed (write) level-",i1," tag ",A)', nlevel, trim(open_tags(nlevel))
#endif
    ELSE
       i = len_trim(tag)
       IF ( i == 0 ) THEN
          WRITE (xmlunit, '("/>")')
#if defined ( __debug )
          print '("closed (write) level-",i1," tag ",A)', &
               nlevel, trim(open_tags(nlevel))
#endif
       ELSE IF ( i == 1 .AND. tag(1:1) == '?' ) THEN
          WRITE (xmlunit, '("?>")')
#if defined ( __debug )
          print '("closed (write) level-",i1," tag ",A)', nlevel, tag
#endif
       ELSE
          WRITE (xmlunit, '("</",A,">")') trim(tag)
#if defined ( __debug )
          print '("closed (write) level-",i1," tag ",A)', nlevel, tag
#endif
       END IF
    END IF
    nlevel = nlevel-1
    !
  END SUBROUTINE xmlw_closetag
  !
  !--------------------------------------------------------
  function xml_protect ( data_in ) result (data_out)
    !--------------------------------------------------------
    ! 
    ! poor-man escaping of a string so that it conforms to xml standard:
    ! replace & with @, < and > with *.
    ! To prevent problems with attributes, double quotes " are replaced
    ! with single quotes '. data_out is left-justified
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
       if ( data_out(i:i) == '"' ) data_out(i:i) = "'"
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
    real(DP_XML), intent(in) :: f
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
  SUBROUTINE readtag_iv (name, ivec, ierr)
    !
    ! As readtag_c, for a vector of integer values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT)         :: ivec(:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    INTEGER :: ier_
    !
    CALL xmlr_opentag (name, ier_)
    if ( ier_ == 0 .or. ier_ == -10 ) then
       READ(xmlunit, *) ivec
       CALL xmlr_closetag ( )
    else
       ivec = 0.0_DP_XML
    end if
    IF ( present (ierr) ) ierr = ier_
    !
  END SUBROUTINE readtag_iv
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
    REAL(DP_XML), INTENT(OUT)         :: rval
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    CHARACTER(LEN=80) :: cval
    !
    CALL readtag_c (name, cval, ierr )
    if ( len_trim(cval) > 0 ) then
       READ (cval,*) rval
    else
       rval = 0.0_DP_XML
    end if
    !
  END SUBROUTINE readtag_r
  !
  SUBROUTINE readtag_rv (name, rvec, ierr)
    !
    ! As readtag_c, for a vector of real values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(DP_XML), INTENT(OUT)         :: rvec(:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    INTEGER :: ier_
    CHARACTER(LEN=90) :: cval
    !
    if ( SIZE(rvec) <= 3 ) then
       ! allow for <name>r1 r2 r3</name> on a single line
       CALL readtag_c (name, cval, ier_ )
       if ( ier_ == 0 .and. len_trim(cval) > 0 ) then
          READ (cval, *, iostat=ier_) rvec
       else
          rvec = 0.0_DP_XML
       end if
    else
       CALL xmlr_opentag (name, ier_)
       if ( ier_ == 0 .or. ier_ == -10 ) then
          READ(xmlunit, *, iostat=ier_) rvec
          CALL xmlr_closetag ( )
       else
          rvec = 0.0_DP_XML
       end if
    end if
    IF ( present (ierr) ) ierr = ier_
    !
  END SUBROUTINE readtag_rv
  !
  SUBROUTINE readtag_rm (name, rmat, ierr)
    !
    ! As readtag_c, for a matrix of real values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(DP_XML), INTENT(OUT)         :: rmat(:,:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    INTEGER :: ier_
    !
    CALL xmlr_opentag (name, ier_)
    if ( ier_ == 0 .or. ier_ == -10 ) then
       READ(xmlunit, *) rmat
       CALL xmlr_closetag ( )
    else
       rmat = 0.0_DP_XML
    end if
    IF ( present (ierr) ) ierr = ier_
    !
  END SUBROUTINE readtag_rm
  !
  SUBROUTINE readtag_rt (name, rtens, ierr)
    !
    ! As readtag_c, for a 3-dim tensor of real values
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(DP_XML), INTENT(OUT)         :: rtens(:,:,:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    INTEGER :: ier_
    !
    CALL xmlr_opentag (name, ier_)
    if ( ier_ == 0 .or. ier_ == -10 ) then
       READ(xmlunit, *) rtens
       CALL xmlr_closetag ( )
    else
       rtens = 0.0_DP_XML
    end if
    IF ( present (ierr) ) ierr = ier_
    !
  END SUBROUTINE readtag_rt
  !
  SUBROUTINE readtag_zv (name, zvec, ierr)
    !
    ! As readtag_c, for a vector of complex values - see comments in writetag_zv
    !
    USE iso_c_binding
    CHARACTER(LEN=*), INTENT(IN) :: name
    COMPLEX(DP_XML), INTENT(OUT), target     :: zvec(:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    !
    TYPE (c_ptr) :: cp
    REAL(DP_XML), POINTER  :: rvec(:)
    INTEGER :: ier_
    !
    CALL xmlr_opentag (name, ier_)
    if ( ier_ == 0 .or. ier_ == -10 ) then
       NULLIFY (rvec)
       cp = c_loc(zvec)
       CALL c_f_pointer ( cp, rvec, shape(zvec)*[2])
       READ( xmlunit, *) rvec
       CALL xmlr_closetag ( )
    else
       zvec = 0.0_DP_XML
    end if
    IF ( present (ierr) ) ierr = ier_
    !
  END SUBROUTINE readtag_zv
  !
  SUBROUTINE readtag_zm (name, zmat, ierr)
    !
    ! As readtag_c, for a matrix of complex values - see comments in writetag_zv
    !
    USE iso_c_binding
    CHARACTER(LEN=*), INTENT(IN) :: name
    COMPLEX(DP_XML), INTENT(OUT), target     :: zmat(:,:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    TYPE (c_ptr) :: cp
    REAL(DP_XML), POINTER  :: rmat(:,:)
    INTEGER :: ier_
    !
    CALL xmlr_opentag (name, ier_)
    if ( ier_ == 0 .or. ier_ == -10 ) then
       NULLIFY (rmat)
       cp = c_loc(zmat)
       CALL c_f_pointer (cp, rmat, shape(zmat)*[2,1])
       READ(xmlunit, *) rmat
       CALL xmlr_closetag ( )
    else
       zmat = 0.0_DP_XML
    end if
    IF ( present (ierr) ) ierr = ier_
    !
  END SUBROUTINE readtag_zm
    !
  SUBROUTINE readtag_zt (name, ztens, ierr)
    !
    ! As readtag_c, for a matrix of complex values - see comments in writetag_zv
    !
    USE iso_c_binding
    CHARACTER(LEN=*), INTENT(IN) :: name
    COMPLEX(DP_XML), INTENT(OUT), target     :: ztens(:,:,:)
    INTEGER, INTENT(OUT),OPTIONAL :: ierr
    TYPE (c_ptr) :: cp
    REAL(DP_XML), POINTER  :: rtens(:,:,:)
    INTEGER :: ier_
    !
    CALL xmlr_opentag (name, ier_)
    if ( ier_ == 0 .or. ier_ == -10 ) then
       NULLIFY (rtens)
       cp = c_loc(ztens)
       CALL c_f_pointer (cp, rtens, shape(ztens)*[2,1,1])
       READ(xmlunit, *) rtens
       CALL xmlr_closetag ( )
    else
       ztens = 0.0_DP_XML
    end if
    IF ( present (ierr) ) ierr = ier_
    !
  END SUBROUTINE readtag_zt
    !
  subroutine readtag_c ( tag, cval, ierr)
    !
    implicit none
    !
    character(len=*), intent(in) :: tag
    character(len=*), intent(out):: cval
    integer, intent(out), optional :: ierr
    !
    integer ::  i, j, lt, ll
    character(len=1) :: endtag
    !
    call xmlr_opentag ( tag, ierr )
    !
    cval = ''
    if ( eot < 0 ) then
       if ( .not. present(ierr) ) then
          print *, 'end of file reached, tag not found'
       else
          ierr = 1
       end if
       return
    else if ( eot == 0 ) then
       ! print *, 'tag found, no value to read on line'
       return
    else
       ! scan current line if there is something after the end of tag 
       ! (variable "eot"); read a new line otherwise
       do while(.true.)
          if ( eot > len_trim(line) ) then
             read(xmlunit,'(a)', end=10) line
             j = 1
          else
             j = eot
          end if
          ! beginning of val at line(j:j): search for end tag
          i = index ( line(j:), '</'//trim(tag) )
          if ( i < 1 ) then
             ! </tag> not found on this line: read value and continue
             cval = trim(cval) // adjustl(trim(line(j:)))
             eot = MAXLINE+1
          else
             ! possible end tag found
             lt = len_trim(tag)
             endtag = adjustl( line(j+i+1+lt:) )
             if ( endtag /= '>' ) then
                if ( .not.present(ierr)) then
                   print *, 'tag ',trim(tag),' not correctly closed'
                else
                   ierr = 2
                endif
             else
                ! end of tag found, read value (if any) and exit
                if ( i > 1 ) cval = trim(cval) // adjustl(trim(line(j:j+i-2)))
                ! print *, 'value=',cval
             end if
#if defined ( __debug )
             print '("closed (read) level-",i1," tag ",A)', &
                     nlevel, trim(open_tags(nlevel))
#endif
             nlevel = nlevel -1
             !
             return
             !
          endif
          !
       end do
       !
    end if
    ! print *, 'tag </',trim(tag),'> not found'
10  if ( present(ierr) ) then
       ierr = 1
    else
       print *, 'end of file reached, tag </'//trim(tag)//'> not found'
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
    ! See list of error codes in the header of this file
    integer :: stat, ntry, ll, lt, i, j, j0
    ! stat=-1: in comment (not actually used)
    ! stat= 0: begin
    ! stat= 1: tag found
    !
    character(len=1) :: quote
    !
    nattr=0
    ntry =0
    if ( allocated(attrlist) ) deallocate (attrlist)
    lt = len_trim(tag)
    !
 1  ntry = ntry+1
    stat=0
    eot =-1
    do while (.true.)
       read(xmlunit,'(a)', end=10) line
       ll = len_trim(line)
       if ( ll > maxline ) then
          print *, 'xmlr_opentag: severe error, line too long'
          if (present(ierr)) ierr = 3
          return
       end if
       ! j is the current scan position
       j = 1
       ! j0 is the start of attributes and values
       j0 = 1
       parse: do while ( j <= ll )
          !
          ! following case is never set and unnecessary:
          !if ( stat ==-1 ) then
          !   ! scanning a comment
          !   i = index(line(j:),'-->')
          !   if ( i == 0 ) then
          !      ! no end of comment found on this line
          !      exit parse
          !   else
          !      ! end of comment found
          !      stat = 0
          !      j = j+i+3
          !   end if
          !else if ( stat == 0 ) then
          if ( stat == 0 ) then
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
                   stat = 1
                   ! <tag continues in next line
                   exit parse
                else if ( line(j:j) == ' ' .or. line(j:j) == '>' &
                                           .or. line(j:j+1)=='/>') then
                   ! <tag or <tag> or <tag/> found
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
                if (present(ierr) .and. ntry == 1) then
                   ierr =-1
                else if (present(ierr) .and. ntry == 2) then
                   ierr =-11
                end if
                ! eot = 0: tag with no value found
                eot = 0
                !
                return
                !
             else if ( line(j:j) == '>' ) then
                ! <tag ... > found
                ! eot points to the rest of the line
                eot = j+1
                if (present(ierr) .and. ntry == 1) then
                   ierr = 0
                else if (present(ierr) .and. ntry == 2) then
                   ierr =-10
                end if
                nlevel = nlevel+1
                IF ( nlevel > maxlevel ) THEN
                   print *, 'xmlr_opentag: severe error, too many levels'
                   if (present(ierr)) ierr = 4
                else
                   open_tags(nlevel) = trim(tag)
#if defined ( __debug )
                   print '("opened (read) level-",i1," tag ",A)',&
                       nlevel, trim(open_tags(nlevel))
#endif
                end if
                !
                return
                !
             else if ( line(j:j) == '=' ) then
                ! end of attribute located: save attribute (with final =)
                nattr=nattr+1
                ! print *, 'attr=',line(j0:j-1)
                if ( nattr == 1 ) then
                   attrlist = line(j0:j)
                else
                   attrlist = attrlist//' '//line(j0:j)
                end if
                ! continue searching for attribute value
                j = j+1
             else if ( line(j:j) == '"' .or. line(j:j) =="'" ) then
                ! first occurrence of ' or " found, look for next
                quote = line(j:j)
                i = index(line(j+1:),quote)
                if ( i < 1 ) then
                   ! print *, 'Error: matching quote not found'
                   go to 10
                else
                   ! save attribute value (with quotes) and continue scanning
                   ! print *, 'attrval=',line(j:j+i-2)
                   attrlist = attrlist//line(j:j+i)
                   j = j+i+1
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
          ierr = 1
          ! quick-and-dirty pseudo-fix to deal with tags not found:
          ! rewind and try again - will work if the desired tag is
          ! found above the current position (and nowhere else)
          rewind(xmlunit)
          if ( ntry == 1 ) go to 1
       else
          print *, 'end of file reached, tag '//trim(tag)//' not found'
       end if
    else
       print *, 'xmlr_opentag: severe parsing error'
       if ( present(ierr) ) ierr = 2
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
    ! stat=-1: in comment (not actually used)
    ! stat= 0: begin
    ! stat= 1: end
    !
    if ( nlevel < 0 ) &
         print '("xmlr_closetag: severe error, closing tag that was never opened")'
    stat=0
#if defined ( __debug )
    if ( .not. present(tag) ) then
       print '("closed (read) level-",i1," tag ",A)', &
            nlevel, trim(open_tags(nlevel))
    else
       print '("closed (read) level-",i1," tag ",A)', nlevel, tag
    end if
#endif
    do while (.true.)
       read(xmlunit,'(a)', end=10) line
       ll = len_trim(line)
       if ( ll > maxline ) then
          print *, 'Fatal error: line too long'
          if (present(ierr)) ierr = 2
          return
       end if
       ! j is the current scan position
       j = 1
       parse: do while ( j <= ll )
          !
          ! following case is never set and unnecessary:
          !if ( stat ==-1 ) then
          ! scanning a comment
          !   i = index(line(j:),'-->')
          !   if ( i == 0 ) then
          ! no end of comment found on this line
          !      exit parse
          !   else
          ! end of comment found
          !      stat = 0
          !      j = j+i+3
          !   end if
          !else if ( stat == 0 ) then
          if ( stat == 0 ) then
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
                   stat = 1
                   ! </tag continues in next line
                   exit parse
                else if ( line(j:j) == ' ' .or. line(j:j) == '>') then
                   ! </tag or </tag> found
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
                !print '("closed")'
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
