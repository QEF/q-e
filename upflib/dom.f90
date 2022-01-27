!
! Copyright (C) 2021-2022 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module dom
  !
  ! Poor-man FoX_dom replacement - Paolo Giannozzi, 2022
  !
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
     type (node), pointer  :: prev => null()
     type (nodelist), allocatable :: linklist
  end type node
  ! Used for check: parsing and cleaning should end at level -1
  integer :: nlevel = -1
  ! The DOM is stored here
  type(node), target, save :: root
  !
  private
  public :: node, nodelist
  public :: parsefile, item, getelementsbytagname, getlength, destroy
  public :: gettagname, hasattribute, extractdataattribute, extractdatacontent
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
  function parsefile ( filename )
    !
    implicit none
    character(len=*), intent (in) :: filename
    type(node), pointer :: parsefile
    integer :: iun, ierr
    !
    type(node), pointer :: curr, next, prev
    integer, parameter :: maxline=1024, maxdim=maxline+16
    character(len=maxdim) :: line
    integer, parameter :: maxlength=80
    character(len=maxlength) :: tag
    ! for checks, not really needed:
    ! integer, parameter :: maxlevel=9
    ! character(len=maxlength), dimension(0:maxlevel) :: open_tags
    !
    logical :: in_comment
    logical :: is_found
    logical :: in_attribute
    logical :: in_data
    integer:: nl, n, n1, n2, m, j
    !
    curr => null()
    in_comment = .false.
    in_attribute = .false.
    in_data = .false.
    nlevel = -1
    ierr = 0 
    !
    open(newunit=iun, file=filename, form='formatted', status='old')
    !
    readline: do
       read(iun,'(a)',end=10) line
       n = 1
       nl = len_trim(line)
       if ( nl > maxline ) then
          ! print *, 'error: line exceeds ', maxline, ' characters'
          ierr = 1
          exit readline
       end if
       tag = ' '
       ! print *,'debug:',trim(line)
       scanline: do while ( n < nl+1 )
          !print *, 'debug: n=',n, line(n:n)
          if ( in_comment ) then
             ! trick to avoid trespassing the EOL
             n1 = min (n+1,nl)
             n2 = min (n+2,nl)
             if ( line(n:n2) == ']]>' ) then
                in_comment = .false.
                n = n+3
                ! print *, 'debug: cdata ended'
             else if ( line(n:n2) == '-->' ) then
                in_comment = .false.
                n = n+3
                ! print *, 'debug: comment ended'
             else if ( line(n:n1) == '?>' ) then
                in_comment = .false.
                n = n+2
                ! print *, 'debug: process ended'
             else
                n = n+1
             end if
          else
             if ( line(n:n) == '<' ) then
                ! trick to avoid trespassing the EOL
                n1 = min (n+1,nl)
                n2 = min (n+3,nl)
                if ( line(n1:n2) == '!--' ) then
                   n = n+4
                   in_comment = .true.
                   ! print *, 'debug: comment begin'
                else if ( line(n1:n1) == '?' ) then
                   n = n+2
                   in_comment = .true.
                   ! print *, 'debug: process begin'
                else if ( line(n1:n1) == '/' ) then
                   ! tag = trim( open_tags(nlevel) )
                   tag = curr%tag
                   n = n+2
                   m = min(n+len_trim(tag)+1,nl)
                   if ( line(n:m) == trim(tag)//'>' ) then
                      ! print *, 'debug: closing tag </',trim(tag),'> found'
                      prev => curr%prev
                      curr => prev
                      in_data = .false.
                   else
                      ! print *, n,m,nlevel,tag
                      ! print *, 'error: unexpected closing tag </',line(n:nl)'> found'
                      ierr = 2
                      exit readline
                   end if
                   nlevel = nlevel - 1
                   n = m+1
                else
                   scantag: do m = n+1, nl
                      is_found = .false.
                      if ( line(m:m) == '>' ) then
                         if ( m == n+1 ) then
                            ! print *, 'error: empty tag <>'
                            ierr = 3
                            exit readline
                         end if
                         is_found = .true.
                         in_data  = .true.
                         in_attribute = .false.
                      else if ( line(m:m) == ' ' .or. m == nl ) then
                         if ( m == n+1 ) then
                            ! print *, 'error: space after <'
                            ierr = 4
                            exit readline
                         end if
                         is_found = .true.
                         in_data  = .false.
                         in_attribute = .true.
                      end if
                      if ( is_found ) then
                         tag = line(n+1:m-1)
                         !if ( in_attribute ) then
                         !   print *, 'debug: tag with attributes ',trim(tag),'...'
                         !else
                         !   print *, 'debug: tag <',trim(tag),'> found'
                         !endif
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
                               ! print *, 'error: more than one root tag'
                               ierr = 5
                               exit readline
                            end if
                            curr => root
                            root = next
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
                      ! print *, 'info short tag ',trim(tag),' found'
                      ! remove slash from attribute
                      curr%attr(len(curr%attr):len(curr%attr)) = ' '
                      prev => curr%prev
                      curr => prev
                      nlevel = nlevel - 1
                   else
                      ! print *, 'debug: tag with attributes ',trim(tag),' found'
                      in_data = .true.
                   end if
                   in_attribute = .false.
                else
                   ! print *, 'error: closed tag that was not open'
                   ierr = 6
                   exit readline
                end if
                n = n+1
             else
                if ( in_attribute ) curr%attr = curr%attr // line(n:n) 
                if ( in_data      ) curr%data = curr%data // line(n:n) 
                n = n+1
             end if
          end if
       end do scanline
       ! if data extends over more than one line, add space between lines
       if ( in_data ) curr%data = curr%data // ' '
    end do readline
10  continue
    if ( ierr == 0 .and. nlevel /= -1) &
         print *, 'error: parsing ended with ',nlevel+1,' level(s) open'
    !
    close(iun)
    if ( ierr > 0 ) then
       print *,'error in parsing: ierr=',ierr
       stop
    end if
    parsefile => root
    !
  end function parsefile
  !
  subroutine add_to_list(linklist, next)
    type(node), pointer :: next
    type(nodelist), allocatable, target :: linklist
    type(nodelist), pointer :: nextlist
    type(nodelist), pointer :: currlist
    !
    if ( .not. allocated(linklist) ) then
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
  recursive subroutine destroy ( curr )
    !
    implicit none
    type(node), pointer :: curr, next
    type(nodelist), pointer :: linklist
    !
    nlevel = nlevel + 1
    ! print *, nlevel, '<', curr%tag,'>, ',curr%attr
    ! print *, curr%data(1:min(80,len(curr%data)))
    if ( allocated( curr%linklist ) ) then
       linklist => curr%linklist
       next  => linklist%node
       lista: do
          call destroy ( next )
          if ( .not. associated( linklist%nextlist ) ) exit lista
          linklist => linklist%nextlist
          next  =>  linklist%node
       end do lista
    end if
    !
    nlevel = nlevel - 1
    if ( associated(curr%prev) ) then
       next => curr%prev
       deallocate(curr)
       curr => next
    else
       ! if ( nlevel /= -1 ) print *, 'destroy: something not right'
       if ( allocated(curr%tag ) ) deallocate (curr%tag)
       if ( allocated(curr%data) ) deallocate (curr%data)
       if ( allocated(curr%attr) ) deallocate (curr%attr)
       if ( allocated(curr%linklist) ) deallocate (curr%linklist)
    end if
    !
  end subroutine destroy
  !
  function getelementsbytagname(root,tag)
    !
    implicit none
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: tag
    type(nodelist), pointer :: getelementsbytagname
    !
    type(nodelist), pointer :: linklist, outlist, newlist
    integer :: n
    !
    n = -1
    getelementsbytagname => null()
    if ( allocated( root%linklist ) ) then
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
    end if
    ! if ( n < 0 ) print *, ' tag: ',tag,' not found'
    !
  end function getelementsbytagname
  !
  function gettagname(root)
    !
    type(node), pointer, intent(in) :: root
    character(len=:), allocatable   :: gettagname
    !
    gettagname = root%tag
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
    la = len_trim(root%attr)
    l0 = len_trim(adjustl(root%attr))
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
  subroutine extractdataattribute_c(root, attr, cval)
    implicit none
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    character(len=*), intent(out) :: cval
    if ( hasattribute(root, attr, cval) ) return
    cval = ''
  end subroutine extractdataattribute_c
  !
  subroutine extractdataattribute_l(root, attr, lval)
    implicit none
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    logical, intent(out) :: lval
    character(len=80) :: val
    if ( hasattribute(root, attr, val) ) then
       read(val,*, end=10,err=10) lval
       return
    end if
    ! not found or not readable
10  lval = .false.
  end subroutine extractdataattribute_l
  !
  subroutine extractdataattribute_i(root, attr, ival)
    implicit none
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    integer, intent(out) :: ival
    character(len=80) :: val
    if ( hasattribute(root, attr, val) ) then
       read(val,*, end=10,err=10) ival
       return
    end if
    ! not found or not readable
10  ival = 0
  end subroutine extractdataattribute_i
  !
  subroutine extractdataattribute_iv(root, attr, ivec)
    implicit none
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    integer, intent(out) :: ivec(:)
    character(len=80) :: val
    if ( hasattribute(root, attr, val) ) then
       read(val,*, end=10,err=10) ivec
       return
    end if
    ! not found or not readable
10  ivec = 0
  end subroutine extractdataattribute_iv
  !
  subroutine extractdataattribute_r(root, attr, rval)
    implicit none
    type(node), pointer, intent(in) :: root
    character(len=*), intent(in) :: attr
    real(dp), intent(out) :: rval
    character(len=80) :: val
    if ( hasattribute(root, attr, val) ) then
       read(val,*, end=10,err=10) rval
       return
    end if
    ! not found or not readable
10  rval = 0
  end subroutine extractdataattribute_r
  !
  subroutine extractdatacontent_c(root, cval, iostat)
    implicit none
    type(node), pointer, intent(in) :: root
    character(len=*), intent(out) :: cval
    integer, intent(out), optional :: iostat
    integer :: ios
    if ( len_trim(root%data) > 0 ) then
       read(root%data,*,iostat=ios) cval
    else
       cval=''
       ios = 0
    end if
    if ( present(iostat) ) iostat=ios
  end subroutine extractdatacontent_c
  !
  subroutine extractdatacontent_l(root, lval, iostat)
    implicit none
    type(node), pointer, intent(in) :: root
    logical, intent(out) :: lval
    integer, intent(out), optional :: iostat
    integer :: ios
    if ( allocated(root%data) ) then
       read(root%data,*,iostat=ios) lval
    else
       lval=.false.
       ios = -1
    end if
    if ( present(iostat) ) iostat=ios
    !
  end subroutine extractdatacontent_l
  !
  subroutine extractdatacontent_i(root, ival, iostat)
    implicit none
    type(node), pointer, intent(in) :: root
    integer, intent(out) :: ival
    integer, intent(out), optional :: iostat
    integer :: ios
    if ( allocated(root%data) ) then
       read(root%data,*,iostat=ios) ival
    else
       ival=0
       ios = -1
    end if
    if ( present(iostat) ) iostat=ios
  end subroutine extractdatacontent_i
  !
  subroutine extractdatacontent_iv(root, ivec, iostat)
    implicit none
    type(node), pointer, intent(in) :: root
    integer, intent(out) :: ivec(:)
    integer, intent(out), optional :: iostat
    integer :: ios, n, iend, ibeg
    !
    if ( allocated(root%data) ) then
       ! the simple solution fails if root%data > 1024 characters:
       ! read(root%data,*,iostat=ios) ivec
       iend = 1
       do n=1,size(ivec)
          call find_token( root%data, ibeg, iend)
          read(root%data(ibeg:iend),*,iostat=ios) ivec(n)
       end do
    else
       ivec= 0
       ios = -1
    end if
    if ( present(iostat) ) iostat=ios
  end subroutine extractdatacontent_iv
  !
  subroutine extractdatacontent_r(root, rval, iostat)
    implicit none
    type(node), pointer, intent(in) :: root
    real(dp), intent(out) :: rval
    integer, intent(out), optional :: iostat
    integer :: ios
    if ( allocated(root%data) ) then
       read(root%data,*,iostat=ios) rval
    else
       rval=0.0_dp
       ios = -1
    end if
    if ( present(iostat) ) iostat=ios
  end subroutine extractdatacontent_r
  !
  subroutine extractdatacontent_rv(root, rvec, iostat)
    implicit none
    type(node), pointer, intent(in) :: root
    real(dp), intent(out) :: rvec(:)
    integer, intent(out), optional :: iostat
    integer :: ios, n, iend, ibeg
    !
    if ( allocated(root%data) ) then
       ! the simple solution fails if root%data > 1024 characters:
       ! read(root%data,*,iostat=ios) rvec
       iend = 1
       do n=1,size(rvec)
          call find_token( root%data, ibeg, iend)
          read(root%data(ibeg:iend),*,iostat=ios) rvec(n)
       end do
    else
       rvec= 0.0_dp
       ios = -1
    end if
    if ( present(iostat) ) iostat=ios
  end subroutine extractdatacontent_rv
  !
  subroutine extractdatacontent_rm(root, rmat, iostat)
    implicit none
    type(node), pointer, intent(in) :: root
    real(dp), intent(out) :: rmat(:,:)
    integer, intent(out), optional :: iostat
    integer :: ios, n, m, iend, ibeg
    !
    if ( allocated(root%data) ) then
       ! the simple solution fails if root%data > 1024 characters:
       ! read(root%data,*,iostat=ios) rmat
       iend = 1
       do m=1,size(rmat,2)
          do n=1,size(rmat,1)
             call find_token( root%data, ibeg, iend)
             read(root%data(ibeg:iend),*,iostat=ios) rmat(n,m)
          end do
       end do
    else
       rmat= 0.0_dp
       ios = -1
    end if
    if ( present(iostat) ) iostat=ios
  end subroutine extractdatacontent_rm
  !
  subroutine find_token ( data, ibeg, iend)
    ! on input:
    !    data    data to be read
    !    iend    1 on first run, end position of previous token otherwise
    ! On output:
    !    data(ibeg:iend) containing a token (a number)
    ! Tokens are assumed to be separated by space or commas
    ! Beware: will not work if empty tokens and multiple commas and present
    !
    character(len=:), allocatable, intent(in) :: data
    integer, intent(out)  :: ibeg
    integer, intent(inout):: iend
    integer:: iscan
    !
    do ibeg = iend, len_trim(data) 
       if ( data(ibeg:ibeg) == ' ' .or.  data(ibeg:ibeg) == ','  ) then
          cycle
       else
          exit
       end if
    end do
    ibeg = min(ibeg, len_trim(data))
    do iend = ibeg, len_trim(data)
       if ( data(iend:iend) /= ' ' .and.  data(iend:iend) /= ','  ) then
          cycle
       else
          exit
       end if
    end do
    iend = min(iend, len_trim(data))
    !
  end subroutine find_token
  !
end module dom
