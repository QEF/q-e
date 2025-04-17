!
! Copyright (C) 2021-2022 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#undef __debug
  !! define __debug to print information on opened and closed tags
module dom
  !!
  !! Poor-man FoX_dom replacement - Paolo Giannozzi, 2022-2024
  !!
  implicit none
  !
  integer, parameter :: dp = selected_real_kind(14,200)
  ! Type 'nodelist' is a "linked list" of pointers to type 'node'
  type :: nodelist
     type(node), pointer :: node
     type(nodelist), pointer :: nextlist
  end type nodelist
  ! Type 'node' contains: tag name, attributes, data, pointer to parent tag,
  ! a linked list of pointers to tags contained in this node
  type :: node
     character(:), allocatable :: tag
     character(:), allocatable :: attr
     character(:), allocatable :: data
     type (node),     pointer  :: prev => null()
     type (nodelist), pointer  :: linklist => null()
  end type node
  ! Used for check: parsing and cleaning should end at level -1
  integer :: nlevel = -1
  ! The DOM is stored here
  type(node), target, save :: root
  ! Used to simulate FoX behavior
  type :: domexception
     integer :: code
  end type domexception
  ! The following machinery is used only to ensure that linked lists
  ! produced by "getelementsbytagname" are properly deallocated when
  ! "destroy" is called. Not sure this is the smartest way to do that.
  ! "meta-list" ml is a linked list to linked lists (!)
  type :: metalist
     type(nodelist), pointer :: linklist => null()
     type(metalist), pointer :: prevmeta => null()
  end type metalist
  type(metalist), pointer :: ml => null()
  !
  private
  ! Callable routines or interfaces
  public :: node, nodelist
  public :: parsefile, item, getelementsbytagname, getlength, destroy
  public :: gettagname, hasattribute, extractdataattribute, extractdatacontent
  public :: getfirstchild, domexception, getexceptioncode
  public :: parsestring
  !
  interface extractdataattribute
     module procedure extractdataattribute_c, &
          extractdataattribute_i, &
          extractdataattribute_iv,&
          extractdataattribute_l, &
          extractdataattribute_r
  end interface extractdataattribute
  !
  interface extractdatacontent
     module procedure extractdatacontent_c, &
          extractdatacontent_cv,&
          extractdatacontent_l, &
          extractdatacontent_i, &
          extractdatacontent_iv,&
          extractdatacontent_r, &
          extractdatacontent_rv,&
          extractdatacontent_rm
  end interface extractdatacontent
  !
CONTAINS
  !
  function parsestring( string, ex )
    !
    character(len=*), intent(in) :: string
    type(domexception), intent (out), optional :: ex
    type(node), pointer :: parsestring
    !
    parsestring => parse( strbuf=string, ex=ex )
    !
  end function parsestring
  !
  function parsefile ( filename, ex )
    !
    character(len=*), intent (in) :: filename
    type(domexception), intent (out), optional :: ex
    type(node), pointer :: parsefile
    integer :: iunit, ierr
    !
    open(newunit=iunit, file=filename, form='formatted', status='old', &
         iostat=ierr)
    !
    if ( ierr /= 0 ) then
       if ( present(ex) ) then
          ex%code = ierr
       else
          print *,'error opening file: ierr=',ierr
          stop
       end if
    else
       !
       parsefile => parse( iun=iunit, ex=ex )
       !
    end if
    !
  end function parsefile
  !
  function parse ( iun, strbuf, ex )
    ! This is where the action is: parse either from unit "iun"
    ! or from a buffer "strbuf", returns a pointer to the node "root"
    ! and optionally an error code in "ex"
    character(len=*), intent (in), optional :: strbuf
    integer, intent(in), optional :: iun
    type(domexception), intent (inout), optional :: ex
    type(node), pointer :: parse
    !
    integer, parameter :: maxline=1024, maxdim=maxline+16
    character(len=maxdim) :: line
    !
    integer :: ierr, i0, i1
    logical :: firstline
    !
    nlevel = -1
    ierr = 0
    i0 = 1
    firstline = .true.
    !
    readline: do
       !
       if ( present(iun) .and. .not. present(strbuf) ) then
          ! read from file
          read(iun,'(a)',end=10) line
       else if ( present(strbuf) .and. .not. present(iun) ) then
          ! read from buffer
          if ( i0 > len(strbuf) ) go to 10
          ! locate newline (ascii n.10)
          i1= index( strbuf(i0:), char(10))
          if ( i1 > 1 ) then
             !  skip LF and go to next line
             line = strbuf(i0:i0+i1-2)
             i0=i0+i1
          else if ( i1 == 1 ) then
             ! empty line
             line = strbuf(i0:i0)
             i0=i0+i1
          else
             ! i1=0: last line (if no LF at the end)
             line = strbuf(i0:)
             i0 = len(strbuf)+1
          end if
       else
          if ( .not.present(ex) ) then
             print *, 'error: both unit and string, or none, in input'
          else
             ierr = 1001
          end if
          exit readline
       end if
  
       if ( len_trim(line) > maxline ) then
          if ( .not.present(ex) ) &
               print *, 'error: line exceeds ', maxline, ' characters'
          ierr = 1
          exit readline
       end if
       ! print *,'debug:',trim(line)
       !
       ierr = parseline ( firstline, line, ex ) 
       if ( ierr /= 0 ) exit readline
       !
    end do readline
    !
10  continue
    if (present(iun) ) close(iun)
    if ( present(ex) ) ex%code = ierr
    !
    if ( ierr == 0 .and. nlevel /= -1) then
       if ( present(ex) ) then
          ex%code = nlevel
       else
          print *, 'error: parsing ended with ',nlevel+1,' level(s) open'
       end if
    else if ( ierr > 0 .and. .not. present(ex) ) then
       print *,'error in parsing: ierr=',ierr
       stop
    end if
    parse => root
    !
  end function parse
  !
  function parseline ( firstline, line, ex )
    ! result is returned in pointer "root"
    ! input line to be parsed
    logical, intent(inout) :: firstline
    character(len=*), intent(in) :: line
    ! optional exception code - if present, error messages are not printed
    type(domexception), intent (in), optional :: ex
    ! error code: 0 = success, otherwise parsing failed
    integer :: parseline
    !
    integer, parameter :: maxlength=80
    ! for checks, not really needed:
    ! integer, parameter :: maxlevel=9
    ! character(len=maxlength), dimension(0:maxlevel) :: open_tags
    !
    ! variables keeping track of parsing status, must be conserved between calls
    ! (tag is actually a local variable, except maybe when used for debugging)
    logical, save :: in_comment
    logical, save :: in_attribute
    logical, save :: in_data
    character(len=maxlength), save :: tag
    ! pointer to current node
    type(node), pointer, save :: curr
    ! next node: the "root" pointer points to the memory allocated in "next"
    type(node), pointer :: next
    ! local variables
    type(node), pointer :: prev
    integer :: n, nl, n1, n2, n3, m
    logical :: is_found
#if defined(_AOCC)
    character(:), allocatable:: temp
#endif
    !
    ! Initialization
    if ( firstline ) then
       firstline = .false.
       curr => null()
       in_comment = .false.
       in_attribute = .false.
       in_data = .false.
       tag = ' '
    end if
    n = 1
    nl= len_trim(line)
    !
    parseline = 0 
    scanline: do while ( n < nl+1 )
       !print *, 'debug: n=',n, line(n:n)
       if ( in_comment ) then
          ! trick to avoid trespassing the EOL
          n1 = min (n+1,nl)
          n2 = min (n+2,nl)
          if ( line(n:n2) == ']]>' ) then
             in_comment = .false.
             n = n+3
#if defined ( __debug )
             print *, 'debug: cdata ends'
#endif
          else if ( line(n:n2) == '-->' ) then
             in_comment = .false.
             n = n+3
#if defined ( __debug )
             print *, 'debug: comment ends'
#endif
          else if ( line(n:n1) == '?>' ) then
             in_comment = .false.
             n = n+2
#if defined ( __debug )
             print *, 'debug: process ends'
#endif
          else
             n = n+1
          end if
       else
          if ( line(n:n) == '<' ) then
             ! trick to avoid trespassing the EOL
             n1 = min (n+1,nl)
             n2 = min (n+3,nl)
             n3 = min (n+8,nl)
             if ( line(n1:n2) == '!--' ) then
                n = n+4
                in_comment = .true.
#if defined ( __debug )
                print *, 'debug: comment begins'
#endif
             else if ( line(n1:n3) == '![CDATA[' ) then
                n = n+9
                in_comment = .true.
#if defined ( __debug )
                print *, 'debug: cdata begins'
#endif
             else if ( line(n1:n1) == '?' ) then
                n = n+2
                in_comment = .true.
#if defined ( __debug )
                print *, 'debug: process begins'
#endif
             else if ( line(n1:n1) == '/' ) then
                ! tag = trim( open_tags(nlevel) )
                tag = curr%tag
                n = n+2
                m = min(n+len_trim(tag)+1,nl)
                if ( line(n:m) == trim(tag)//'>' ) then
#if defined ( __debug )
                   print *, 'debug: closing tag </',trim(tag),'> found'
#endif
                   prev => curr%prev
                   curr => prev
                   in_data = .false.
                else
                   if ( .not.present(ex) ) then
                      print *, n,m,nlevel,tag
                      print *, 'error: unexpected closing tag </',line(n:nl),'> found'
                   end if
                   parseline = 2
                   return
                end if
                nlevel = nlevel - 1
                n = m+1
             else
                scantag: do m = n+1, nl
                   is_found = .false.
                   if ( line(m:m) == '>' ) then
                      if ( m == n+1 ) then
                         if ( .not.present(ex) ) &
                              print *, 'error: empty tag <>'
                          parseline = 3
                         return
                      end if
                      is_found = .true.
                      in_data  = .true.
                      in_attribute = .false.
                   else if ( line(m:m) == ' ' .or. line(m:m) == '/' &
                                              .or. m == nl ) then
                      ! case '/' may occur for empty tags like "<tag/>"
                      if ( m == n+1 ) then
                         if ( .not.present(ex) ) &
                              print *, 'error: space after <'
                         parseline = 4
                         return
                      end if
                      is_found = .true.
                      in_data  = .false.
                      in_attribute = .true.
                   end if
                   if ( is_found ) then
                      tag = line(n+1:m-1)
#if defined ( __debug )
                      if ( in_attribute ) then
                         print *, 'debug: tag with attributes ',trim(tag),'...'
                      else
                         print *, 'debug: tag <',trim(tag),'> found'
                      endif
#endif
                      nlevel = nlevel + 1
                      ! open_tags(nlevel) = trim(tag)
                      allocate(next)
                      next%tag  = trim(tag)
                      if ( in_attribute ) next%attr=' '
                      if ( in_data      ) next%data=' '
                      if ( associated(curr) ) then
                         next%prev => curr
                         call add_to_list(curr%linklist,next)
                         curr => next
                      else
                         if ( allocated(root%tag) ) then
                            if ( .not.present(ex) ) &
                                 print *, 'error: more than one root tag'
                            parseline = 5
                            return
                         end if
                         curr => root
                         root = next
                         deallocate(next)
                      end if
                      !
                      n = m+1
                      exit scantag
                   end if
                end do scantag
                if ( m > nl) then
                   tag = ' '
                   n = nl+1
                   exit scanline
                end if
             end if
          else if ( line(n:n) == '>' ) then
             if ( in_attribute ) then
                if ( line(n-1:n-1) == '/' ) then
#if defined ( __debug )
                   print *, 'info short tag ',trim(tag),' found'
#endif
                   ! remove slash from attribute
                   curr%attr(len(curr%attr):len(curr%attr)) = ' '
                   prev => curr%prev
                   curr => prev
                   nlevel = nlevel - 1
                else
#if defined ( __debug )
                   print *, 'debug: tag with attributes ',trim(tag),' found'
#endif
                   in_data = .true.
                end if
                in_attribute = .false.
             else
                if ( .not.present(ex) ) &
                     print *, 'error: closed tag that was not open'
                parseline = 6
                return
             end if
             n = n+1
          else
             if ( in_attribute ) then
                if ( .not. allocated(curr%attr) ) curr%attr = ' '
#if defined(_AOCC)
                temp = curr%attr // line(n:n) 
                curr%attr = temp
#else
                curr%attr = curr%attr // line(n:n) 
#endif
             end if
             if ( in_data      ) then
                if ( .not. allocated(curr%data) ) curr%data = ' '
#if defined(_AOCC)
                temp = curr%data // line(n:n) 
                curr%data = temp
#else
                curr%data = curr%data // line(n:n) 
#endif
             end if
             n = n+1
          end if
       end if
    end do scanline
    ! if data extends over more than one line, add space between lines
    if ( in_data .and. associated(curr) ) then
       if ( allocated(curr%data) ) then
#if defined(_AOCC)
          temp = curr%data // ' ' 
          curr%data = temp
#else
          curr%data = curr%data // ' '
#endif
          end if
    end if

  end function parseline
  
  integer function getexceptioncode(ex)
     type(domexception), intent(in):: ex
     getexceptioncode = ex%code
  end function getexceptioncode
  !
  subroutine add_to_list(linklist, next)
    type(node), pointer :: next
    type(nodelist), pointer :: linklist
    type(nodelist), pointer :: nextlist
    type(nodelist), pointer :: currlist
    !
    if ( .not. associated(linklist) ) then
       allocate(linklist)
       linklist%node => next
       linklist%nextlist => null()
    else
       currlist => linklist
       do while ( associated(currlist%nextlist) )
          currlist => currlist%nextlist
       end do
       allocate(nextlist)
       nextlist%node => next
       nextlist%nextlist => null()
       currlist%nextlist => nextlist
    end if
    !
  end subroutine add_to_list
  !
  recursive subroutine destroy ( curr, iun )
    !
    ! This (obscure) code goes down recursively into the "curr" tree,
    ! then deallocates (hopefully) everything. If "iun" is present,
    ! the tree is reprinted to unit "iun". Useful for debugging: the
    ! reprinted tree should have the same structure as the original file
    !
    type(node), pointer :: curr, next
    type(nodelist), pointer :: linklist, nextlist
    integer, intent(in), optional :: iun
    !
    nlevel = nlevel + 1
    ! print *, nlevel, '<', curr%tag,'>, ',curr%attr
    ! print *, curr%data(1:min(80,len(curr%data)))
    if ( present(iun ) ) then
       if ( allocated(curr%attr) ) then
          write(iun,'("<",A," ",A,">")') trim(curr%tag),trim(curr%attr)
       else
          write(iun,'("<",A,">")') trim(curr%tag)
       end if
       if ( allocated(curr%data) ) write(iun,'(A)') trim(curr%data)
    end if
    ! Go down recursively on the tree (note the call to itself below)
    linklist => curr%linklist
    do while ( associated(linklist) )
       call destroy(linklist%node, iun)
       nextlist => linklist%nextlist
       deallocate (linklist)
       ! The linked list must be explicitly deallocated to avoid memory leaks
       linklist => nextlist
    end do
    !
    if ( present(iun ) ) write(iun,'("</",A,">")') trim(curr%tag)
    nlevel = nlevel - 1
    ! now deallocate all memory
    if ( allocated(curr%tag) ) deallocate (curr%tag)
    if ( allocated(curr%data) ) deallocate (curr%data)
    if ( allocated(curr%attr) ) deallocate (curr%attr)
    !
    if ( associated(curr%prev) ) then
       ! go down one level and deallocate
       next => curr%prev
       deallocate(curr)
       curr => next
    else
       call destroyml ( )
       if ( nlevel /= -1 ) print *, 'destroy: did not reach root level?'
    end if
    !
  end subroutine destroy
  !
  subroutine destroyll (linklist)
    type(nodelist), pointer :: linklist
    type(nodelist), pointer :: nextlist
    do while ( associated(linklist) )
       nextlist => linklist%nextlist
       deallocate(linklist)
       ! if ( .not.associated(nextlist) ) exit
       linklist =>nextlist
    end do
  end subroutine destroyll
  !
  subroutine destroyml ( )
    type(metalist), pointer :: prevml
    ! deallocate all linked lists by going back into "ml"
    do while ( associated(ml) )
       call destroyll(ml%linklist)
       prevml => ml%prevmeta
       deallocate (ml)
       ml => prevml
    end do
  end subroutine destroyml
  !
  function getelementsbytagname(root,tag)
    !
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: tag
    type(nodelist), pointer :: getelementsbytagname
    !
    type(nodelist), pointer :: linklist, outlist, newlist
    type(metalist), pointer :: nextml
    integer :: n
    !
    n = -1
    getelementsbytagname => null()
    if ( associated( root%linklist ) ) then
       linklist => root%linklist
       lista: do
          if ( trim(adjustl(tag)) == linklist%node%tag ) then
             n = n+1
             ! print *, 'info: tag: ',tag,' found: n=',n
             allocate(newlist)
             newlist%node => linklist%node
             newlist%nextlist => null()
             if ( n == 0 ) then
                getelementsbytagname => newlist
                outlist => getelementsbytagname
             else
                outlist%nextlist => newlist
                outlist => newlist
             end if
          end if
          if ( .not. associated( linklist%nextlist ) ) exit lista
          linklist => linklist%nextlist
       end do lista
       !
       ! Store linked list at the end of "meta-list" ml,
       ! keeping track of the previous one, for later deallocation
       !
       if ( .not.associated(ml) ) then
          allocate(ml)
       else
          allocate(nextml)
          nextml%prevmeta => ml
          ml => nextml
       end if
       ml%linklist => getelementsbytagname
       !
    end if
    ! if ( n < 0 ) print *, ' tag: ',tag,' not found'
    !
  end function getelementsbytagname
  !
  function getfirstchild(root, ex)
    !
    type(node), pointer, intent(in) :: root
    type(domexception), intent (out), optional :: ex
    type(node), pointer :: getfirstchild
    !
    if ( associated( root ) ) then
       getfirstchild => root
       if ( present(ex) ) ex%code = 0
    else
       getfirstchild => null()
       if ( present(ex) ) ex%code = 1
    endif
    !
  end function getfirstchild
  !
  function gettagname(root, ex)
    !
    type(node), pointer, intent(in) :: root
    type(domexception), intent (out), optional :: ex
    character(len=:), allocatable   :: gettagname
    !
    gettagname = root%tag
    ! ignored
    if ( present(ex) ) ex%code = 0
    !
  end function gettagname
  !
  integer function getlength (llist) result(n)
    !
    type(nodelist), pointer, intent(in) :: llist
    type(nodelist), pointer :: mylist
    !
    n = 0
    if ( .not.associated(llist) ) return
    n = 1
    mylist => llist
    lista: do while( associated(mylist%nextlist) )
       n = n + 1
       mylist => mylist%nextlist
    end do lista
    !
  end function getlength
  !
  function item (llist,n)
    !
    type(nodelist), pointer, intent(in) :: llist
    integer, intent(in) :: n
    type(node), pointer :: item
    !
    type(nodelist), pointer :: mylist
    integer :: i
    !
    item => null()
    if ( .not.associated(llist) ) return
    mylist => llist
    lista: do i=0,n-1
       mylist => mylist%nextlist
    end do lista
    item => mylist%node
    !
  end function item
  !
  logical function hasattribute(root,attr,val) result(found)
    !
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    character(len=*), intent(out), optional :: val
    !
    integer :: la, l0, i1, i2, i
    logical :: in_attrval
    character(len=1) :: delimiter
    !
    if(allocated(root%attr)) then
       la = len_trim(root%attr)
       l0 = len_trim(adjustl(root%attr))
    else
       la = 0
       l0 = 0
    endif
    in_attrval=.false.
    found = .false.
    i1 = 0
scan: do i=la-l0+1,la
       if ( .not. in_attrval) then
          if (root%attr(i:i) == '"' .or. root%attr(i:i) == "'") then
             in_attrval=.true.
             delimiter = root%attr(i:i)
             ! write(*,'("attr:",a,", ")',advance='no') root%attr(i1:i2)
             found = ( attr == root%attr(i1:i2) )
             i1 = i+1
             i2 = 0
          else if (i1 == 0 .and. root%attr(i:i) /= ' ') then
             i1 = i
          else if (i1  > 0 .and. root%attr(i:i) /= ' ' .and. root%attr(i:i) /= '=' ) then
             i2 = i
          end if
       else
          if (root%attr(i:i) == delimiter ) then
             in_attrval=.false.
             i2 = i-1
             ! write(*,'("value:",a,".")') root%attr(i1:i2)
             if ( present(val) ) val = root%attr(i1:i2)
             if ( found ) return
             i1 = 0
          end if
       end if
    end do scan
    !
  end function hasattribute
  !
  subroutine extractdataattribute_c(root, attr, cval, iostat)
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    character(len=*), intent(out) :: cval
    integer, intent(out), optional:: iostat
    !
    if (present(iostat)) iostat=0
    if ( hasattribute(root, attr, cval) ) return
    if (present(iostat)) iostat=1
    cval = ' '
    !
  end subroutine extractdataattribute_c
  !
  subroutine extractdataattribute_l(root, attr, lval, iostat)
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    logical, intent(out) :: lval
    character(len=80) :: val
    integer, intent(out), optional:: iostat
    !
    if (present(iostat)) iostat=0
    if ( hasattribute(root, attr, val) ) then
       read(val,*, end=10,err=10) lval
       return
    end if
    ! not found or not readable
10  lval = .false.
    if (present(iostat)) iostat=1
    !
  end subroutine extractdataattribute_l
  !
  subroutine extractdataattribute_i(root, attr, ival, iostat)
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    integer, intent(out) :: ival
    character(len=80) :: val
    integer, intent(out), optional:: iostat
    !
    if (present(iostat)) iostat=0
    if ( hasattribute(root, attr, val) ) then
       read(val,*, end=10,err=10) ival
       return
    end if
    ! not found or not readable
10  ival = 0
    if (present(iostat)) iostat=1
    !
  end subroutine extractdataattribute_i
  !
  subroutine extractdataattribute_iv(root, attr, ivec, iostat)
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    integer, intent(out) :: ivec(:)
    character(len=80) :: val
    integer, intent(out), optional:: iostat
    !
    if (present(iostat)) iostat=0
    if ( hasattribute(root, attr, val) ) then
       read(val,*, end=10,err=10) ivec
       return
    end if
    ! not found or not readable
10  ivec = 0
    if (present(iostat)) iostat=1
    !
  end subroutine extractdataattribute_iv
  !
  subroutine extractdataattribute_r(root, attr, rval, iostat)
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    real(dp), intent(out) :: rval
    character(len=80) :: val
    integer, intent(out), optional:: iostat
    !
    if (present(iostat)) iostat=0
    if ( hasattribute(root, attr, val) ) then
       read(val,*, end=10,err=10) rval
       return
    end if
    ! not found or not readable
10  rval = 0
    if (present(iostat)) iostat=1
    !
  end subroutine extractdataattribute_r
  !
  subroutine extractdatacontent_c(root, cval, iostat)
    type(node), pointer, intent(in) :: root
    character(len=*), intent(out) :: cval
    integer, intent(out), optional :: iostat
    integer :: ios
    !
    cval=' '
    ios = 0
    if ( allocated(root%data) ) then
       if ( len_trim(root%data) > 0 ) read(root%data,*,iostat=ios) cval
    else
       ios = 1
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_c
  !
  subroutine extractdatacontent_cv(root, cvec, iostat)
    type(node), pointer, intent(in) :: root
    character(len=*), intent(inout), pointer :: cvec(:)
    integer, intent(out), optional :: iostat
    integer :: ibeg, iend, n, ios
    !
    cvec(:) = ' '
    ios = 0
    if ( allocated(root%data) ) then
       if ( len_trim(root%data) > 0 ) then
          iend = 0
          do n=1,size(cvec)
             ios = find_token( root%data, ibeg, iend)
             if ( ios == 0 ) then
                cvec(n) = root%data(ibeg:iend)
             else
                cvec(n) = ' '
             end if
          end do
       end if
    else
       ios = 1
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_cv
  !
  subroutine extractdatacontent_l(root, lval, iostat)
    type(node), pointer, intent(in) :: root
    logical, intent(out) :: lval
    integer, intent(out), optional :: iostat
    integer :: ios
    if ( allocated(root%data) ) then
       read(root%data,*,iostat=ios) lval
    else
       lval=.false.
       ios = 1
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_l
  !
  subroutine extractdatacontent_i(root, ival, iostat)
    type(node), pointer, intent(in) :: root
    integer, intent(out) :: ival
    integer, intent(out), optional :: iostat
    integer :: ios
    !
    if ( allocated(root%data) ) then
       read(root%data,*,iostat=ios) ival
    else
       ival=0
       ios = 1
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_i
  !
  subroutine extractdatacontent_iv(root, ivec, iostat)
    type(node), pointer, intent(in) :: root
    integer, intent(out) :: ivec(:)
    integer, intent(out), optional :: iostat
    integer :: ios, n, iend, ibeg
    !
    ios = 1
    if ( allocated(root%data) ) then
       ! the simple solution fails if root%data > 1024 characters:
       ! read(root%data,*,iostat=ios) ivec
       iend = 0
       do n=1,size(ivec)
          ios = find_token( root%data, ibeg, iend)
          if ( ios == 0 ) then
             read(root%data(ibeg:iend),*,iostat=ios) ivec(n)
          else
             ivec(n) = 0
          end if
       end do
    else
       ivec(:) = 0
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_iv
  !
  subroutine extractdatacontent_r(root, rval, iostat)
    type(node), pointer, intent(in) :: root
    real(dp), intent(out) :: rval
    integer, intent(out), optional :: iostat
    integer :: ios
    !
    if ( allocated(root%data) ) then
       read(root%data,*,iostat=ios) rval
    else
       rval=0.0_dp
       ios = 1
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_r
  !
  subroutine extractdatacontent_rv(root, rvec, iostat)
    type(node), pointer, intent(in) :: root
    real(dp), intent(out) :: rvec(:)
    integer, intent(out), optional :: iostat
    integer :: ios, n, iend, ibeg
    !
    ios = 1 
    if ( allocated(root%data) ) then
       ! the simple solution fails if root%data > 1024 characters:
       ! read(root%data,*,iostat=ios) rvec
       iend = 0
       do n=1,size(rvec)
          ios = find_token( root%data, ibeg, iend)
          if ( ios == 0 ) then
             read(root%data(ibeg:iend),*,iostat=ios) rvec(n)
          else
             rvec(n) = 0.0_dp
          end if
       end do
    else
       rvec(:) = 0.0_dp
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_rv
  !
  subroutine extractdatacontent_rm(root, rmat, iostat)
    type(node), pointer, intent(in) :: root
    real(dp), intent(out) :: rmat(:,:)
    integer, intent(out), optional :: iostat
    integer :: ios, n, m, iend, ibeg
    !
    ios = 1
    if ( allocated(root%data) ) then
       ! the simple solution fails if root%data > 1024 characters:
       ! read(root%data,*,iostat=ios) rmat
       iend = 0
       do m=1,size(rmat,2)
          do n=1,size(rmat,1)
             ios = find_token( root%data, ibeg, iend)
             if ( ios == 0 ) then
                read(root%data(ibeg:iend),*,iostat=ios) rmat(n,m)
             else
                rmat(n,m) = 0.0_dp
             end if
          end do
       end do
    else
       rmat = 0.0_dp
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_rm
  !
  integer function find_token ( data, ibeg, iend )
    !
    ! Locate tokens (numbers, fields) in a string
    ! Tokens are assumed to be separated by space or commas
    !
    ! on input:
    !    data    string containing tokens
    !    iend    0 on first run, end position of previous token otherwise
    ! On output:
    !    find_token   0 if token found, 1 otherwise
    !    ibeg, iend   if find_token, data(ibeg:iend) contains a token
    !                 if not, ibeg=iend, iend unchanged
    !
    ! Beware: will not work if empty tokens and multiple commas and present,
    !         e.g.: "field1, ,field3" where field2 is empty
    !
    character(len=:), allocatable, intent(in) :: data
    integer, intent(out)  :: ibeg
    integer, intent(inout):: iend
    integer:: lt
    !
    lt = len_trim(data) 
    find_token = 1
    do ibeg = iend+1, lt
       if ( data(ibeg:ibeg) == ' ' .or.  data(ibeg:ibeg) == ','  ) then
          cycle
       else
          find_token = 0
          exit
       end if
    end do
    if ( find_token == 0 ) then
       do iend = ibeg, lt
          if ( data(iend:iend) /= ' ' .and.  data(iend:iend) /= ','  ) then
             cycle
          else
             exit
          end if
       end do
       iend = min(iend, lt)
    else
       ibeg = iend
    end if
    !
  end function find_token
  !
end module dom
