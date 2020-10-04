! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004-2006 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

#include "iotk_auxmacros.h"


!-<
! ... doesn't work with qe_syntax
!->
subroutine iotk_copy_tag_x(source,dest,dummy,maxsize,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_scan_interf
  use iotk_write_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,           intent(in)  :: source
  integer,           intent(in)  :: dest
  type(iotk_dummytype),optional  :: dummy
  integer, optional, intent(in)  :: maxsize
  integer, optional, intent(out) :: ierr
  logical :: source_binary,dest_binary,source_stream,dest_stream
  integer :: ierrl,control,maxsizel
  character(iotk_taglenx) :: tag
  character(iotk_namlenx) :: name
  character(iotk_attlenx) :: attr
  character(iotk_vallenx) :: type
  type(iotk_unit), pointer :: this
  ierrl = 0         
  maxsizel = -1     
  if(present(maxsize)) maxsizel = maxsize 
  call iotk_inquire(source,source_binary,source_stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
    goto 1
  end if
  call iotk_inquire(dest,dest_binary,dest_stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")            
    goto 1
  end if
  call iotk_unit_get(source,pointer=this)
  if(.not.associated(this)) then
    call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'unit')  
    goto 1
  end if
  do
    call iotk_scan_tag(source,+1,control,tag,source_binary,source_stream,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
      goto 1
    end if
    if(control/=4) then ! SKIP FOR COMMENTS
      call iotk_tag_parse(tag,name,attr,ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
        goto 1
      end if
    end if
    if(iotk_strcomp(name,this%root)) then
      call iotk_scan_tag(source,-1,control,tag,source_binary,source_stream,ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
        goto 1
      end if
      return
    end if
    select case(control)
    case(1)
      call iotk_scan_attr(attr,"type",type,ierr=ierrl,eos=.true.,default=" ")
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
        goto 1
      end if
      if((iotk_strcomp(type,"real") .or. iotk_strcomp(type,"integer") .or. iotk_strcomp(type,"logical") &
                                    .or. iotk_strcomp(type,"character") .or. iotk_strcomp(type,"complex")) .and. control==1) then
        call iotk_copy_dat(source,dest,source_binary,dest_binary,name,attr,maxsize=maxsizel,ierr=ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
          goto 1
        end if
        call iotk_scan_tag(source,+1,control,tag,source_binary,source_stream,ierrl)
      else
        call iotk_write_begin(dest,name,attr,ierr=ierrl)
      end if
    case(2)
      call iotk_write_end(dest,name,ierr=ierrl)
    case(3)
      call iotk_write_empty(dest,name,attr,ierr=ierrl)
    case(4)
      call iotk_write_comment(dest,tag,ierr=ierrl)
    case(5)
      call iotk_write_pi(dest,name,attr,ierr=ierrl)
    end select
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_copy_tag",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
      goto 1
    end if
  end do
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_copy_tag_x

subroutine iotk_parse_dat_x(attr,type,ikind,isize,ilen,fmt,columns,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in)  :: attr
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: type
#else
  character(len=*), intent(out) :: type
#endif
  integer,          intent(out) :: ikind
  integer,          intent(out) :: isize
  integer,          intent(out) :: ilen
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: fmt
#else
  character(len=*), intent(out) :: fmt
#endif
  integer,          intent(out) :: columns
  integer,          intent(out) :: ierr
  character(iotk_vallenx) :: typename
  ierr = 0
  call iotk_scan_attr(attr,"type",typename,ierr=ierr,eos=.true.,default=iotk_eos)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_parse_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  type = iotk_toupper(typename)
  call iotk_scan_attr(attr,"kind",ikind,ierr=ierr,default=-1)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_parse_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  call iotk_scan_attr(attr,"size",isize,ierr=ierr,default=-1)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_parse_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  call iotk_scan_attr(attr,"len", ilen, ierr=ierr,default=-1)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_parse_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  call iotk_scan_attr(attr,"fmt", fmt, ierr=ierr,eos=.true.,default="!"//iotk_eos)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_parse_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  call iotk_scan_attr(attr,"columns",columns,ierr=ierr,default=1)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_parse_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
end subroutine iotk_parse_dat_x

subroutine iotk_set_x(dummy,unitmin,unitmax,getline_buffer,error_warn_overflow, &
                      linlen,indent,maxindent,error_unit,output_unit,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_dummytype),optional  :: dummy
  integer, optional, intent(in)  :: unitmin
  integer, optional, intent(in)  :: unitmax
  integer, optional, intent(in)  :: getline_buffer
  logical, optional, intent(in)  :: error_warn_overflow
  integer, optional, intent(in)  :: linlen
  integer, optional, intent(in)  :: indent
  integer, optional, intent(in)  :: maxindent
  integer, optional, intent(in)  :: error_unit
  integer, optional, intent(in)  :: output_unit
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  ierrl = 0
  if(present(error_warn_overflow)) then
    iotk_error_warn_overflow = error_warn_overflow
  end if
  if(present(unitmin)) then
    if(unitmin<0) then
      call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Wrong value for unitmin')
call iotk_error_write(ierrl,"unitmin",unitmin)
      goto 1 
    end if
    iotk_unitmin = unitmin 
  end if
  if(present(unitmax)) then
    if(unitmax<0) then
      call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Wrong value for unitmax')
call iotk_error_write(ierrl,"unitmax",unitmax)
      goto 1
    end if
    iotk_unitmax = unitmax
  end if
  if(iotk_unitmin>iotk_unitmax) then
    call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Inconsistency: unitmin should be less then unitmax')
call iotk_error_write(ierrl,"iotk_unitmin",iotk_unitmin)
call iotk_error_write(ierrl,"iotk_unitmax",iotk_unitmax)
    goto 1
  end if
  if(present(getline_buffer)) then
    if(getline_buffer<1) then
      call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Wrong value for getline_buffer')
call iotk_error_write(ierrl,"getline_buffer",getline_buffer)
      goto 1
    end if
    iotk_getline_buffer = getline_buffer
  end if
  if(present(linlen)) then
    if(linlen<1) then
      call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Wrong value for linlen')
call iotk_error_write(ierrl,"linlen",linlen)
      goto 1
    end if
    iotk_linlen = linlen
  end if
  if(present(indent)) then
    if(indent<0) then
      call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Wrong value for indent')
call iotk_error_write(ierrl,"indent",indent)
      goto 1
    end if
    iotk_indent = indent
  end if
  if(present(maxindent)) then
    if(maxindent<0 .or. maxindent>iotk_linlenx) then
      call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Wrong value for maxindent, should be between 0 and iotk_linlenx')
call iotk_error_write(ierrl,"maxindent",maxindent)
call iotk_error_write(ierrl,"iotk_linlenx",iotk_linlenx)
      goto 1
    end if
    iotk_maxindent = maxindent
  end if
  if(present(error_unit)) then
    if(error_unit<0) then
      call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Wrong value for error_unit')
call iotk_error_write(ierrl,"error_unit",error_unit)
      goto 1
    end if
    iotk_error_unit = error_unit
  end if
  if(present(output_unit)) then
    if(output_unit<0) then
      call iotk_error_issue(ierrl,"iotk_set",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
call iotk_error_msg(ierrl,'Wrong value for output_unit')
call iotk_error_write(ierrl,"output_unit",output_unit)
      goto 1
    end if
    iotk_output_unit = output_unit
  end if
1 continue 
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_set_x

subroutine iotk_get_x(dummy,unitmin,unitmax,getline_buffer,error_warn_overflow, &
                              linlen,indent,maxindent,error_unit,output_unit)
  use iotk_base
  use iotk_misc_interf
  implicit none
  type(iotk_dummytype),optional  :: dummy
  integer, optional, intent(out) :: unitmin
  integer, optional, intent(out) :: unitmax
  integer, optional, intent(out) :: getline_buffer
  logical, optional, intent(out) :: error_warn_overflow
  integer, optional, intent(out) :: linlen
  integer, optional, intent(out) :: indent
  integer, optional, intent(out) :: maxindent
  integer, optional, intent(out) :: error_unit
  integer, optional, intent(out) :: output_unit
  if(present(unitmin)) unitmin = iotk_unitmin
  if(present(unitmax)) unitmax = iotk_unitmax
  if(present(getline_buffer)) getline_buffer = iotk_getline_buffer
  if(present(error_warn_overflow)) error_warn_overflow = iotk_error_warn_overflow
  if(present(linlen)) linlen = iotk_linlen
  if(present(indent)) indent = iotk_indent
  if(present(maxindent)) maxindent = iotk_maxindent
  if(present(error_unit)) error_unit = iotk_error_unit
  if(present(output_unit)) output_unit = iotk_output_unit
end subroutine iotk_get_x

subroutine iotk_print_kinds_x
  use iotk_base
  use iotk_misc_interf
  use iotk_xtox_interf
  implicit none
  character(100) :: string
  write(*,"(a,i5)") "Maximum rank            : ", iotk_maxrank
  write(*,"(a,i5)") "Maximum rank hard limit : ", iotk_maxrank
#ifdef __IOTK_LOGICAL1
  string = "logical(kind="//trim(iotk_itoa(iotk_LOGICAL1))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_LOGICAL2
  string = "logical(kind="//trim(iotk_itoa(iotk_LOGICAL2))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_LOGICAL3
  string = "logical(kind="//trim(iotk_itoa(iotk_LOGICAL3))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_LOGICAL4
  string = "logical(kind="//trim(iotk_itoa(iotk_LOGICAL4))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_INTEGER1
  string = "integer(kind="//trim(iotk_itoa(iotk_INTEGER1))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_INTEGER2
  string = "integer(kind="//trim(iotk_itoa(iotk_INTEGER2))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_INTEGER3
  string = "integer(kind="//trim(iotk_itoa(iotk_INTEGER3))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_INTEGER4
  string = "integer(kind="//trim(iotk_itoa(iotk_INTEGER4))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_REAL1
  string = "real(kind="//trim(iotk_itoa(iotk_REAL1))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_REAL2
  string = "real(kind="//trim(iotk_itoa(iotk_REAL2))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_REAL3
  string = "real(kind="//trim(iotk_itoa(iotk_REAL3))//")"
  write(*,"(a)") trim(string)
#endif
#ifdef __IOTK_REAL4
  string = "real(kind="//trim(iotk_itoa(iotk_REAL4))//")"
  write(*,"(a)") trim(string)
#endif
  string = "character(kind="//trim(iotk_itoa(iotk_CHARACTER1))//")"
  write(*,"(a)") trim(string)
end subroutine iotk_print_kinds_x

!-<
! ... doesn't work with qe_syntax
!->
subroutine iotk_copy_dat_aux_x(source,dest,source_binary,dest_binary,name,type,ikind,isize, &
                               ilen,fmt,columns,attr,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,      intent(in)  :: source
  integer,      intent(in)  :: dest
  logical,      intent(in)  :: source_binary
  logical,      intent(in)  :: dest_binary
  character(*), intent(in)  :: name
  character(*), intent(in)  :: type
  integer,      intent(in)  :: ikind
  integer,      intent(in)  :: isize
  integer,      intent(in)  :: ilen
  character(*), intent(in)  :: fmt
  integer,      intent(in)  :: columns
  character(*), intent(in)  :: attr
  integer,      intent(out) :: ierr
  
  integer :: tmpkind
#ifdef __IOTK_LOGICAL1
  LOGICAL (kind=iotk_LOGICAL1), allocatable :: dat_LOGICAL1 (:)
#endif
#ifdef __IOTK_LOGICAL2
  LOGICAL (kind=iotk_LOGICAL2), allocatable :: dat_LOGICAL2 (:)
#endif
#ifdef __IOTK_LOGICAL3
  LOGICAL (kind=iotk_LOGICAL3), allocatable :: dat_LOGICAL3 (:)
#endif
#ifdef __IOTK_LOGICAL4
  LOGICAL (kind=iotk_LOGICAL4), allocatable :: dat_LOGICAL4 (:)
#endif
#ifdef __IOTK_INTEGER1
  INTEGER (kind=iotk_INTEGER1), allocatable :: dat_INTEGER1 (:)
#endif
#ifdef __IOTK_INTEGER2
  INTEGER (kind=iotk_INTEGER2), allocatable :: dat_INTEGER2 (:)
#endif
#ifdef __IOTK_INTEGER3
  INTEGER (kind=iotk_INTEGER3), allocatable :: dat_INTEGER3 (:)
#endif
#ifdef __IOTK_INTEGER4
  INTEGER (kind=iotk_INTEGER4), allocatable :: dat_INTEGER4 (:)
#endif
#ifdef __IOTK_REAL1
  REAL (kind=iotk_REAL1), allocatable :: dat_REAL1 (:)
#endif
#ifdef __IOTK_REAL2
  REAL (kind=iotk_REAL2), allocatable :: dat_REAL2 (:)
#endif
#ifdef __IOTK_REAL3
  REAL (kind=iotk_REAL3), allocatable :: dat_REAL3 (:)
#endif
#ifdef __IOTK_REAL4
  REAL (kind=iotk_REAL4), allocatable :: dat_REAL4 (:)
#endif
#ifdef __IOTK_COMPLEX1
  COMPLEX (kind=iotk_COMPLEX1), allocatable :: dat_COMPLEX1 (:)
#endif
#ifdef __IOTK_COMPLEX2
  COMPLEX (kind=iotk_COMPLEX2), allocatable :: dat_COMPLEX2 (:)
#endif
#ifdef __IOTK_COMPLEX3
  COMPLEX (kind=iotk_COMPLEX3), allocatable :: dat_COMPLEX3 (:)
#endif
#ifdef __IOTK_COMPLEX4
  COMPLEX (kind=iotk_COMPLEX4), allocatable :: dat_COMPLEX4 (:)
#endif
#ifdef __IOTK_CHARACTER1
  CHARACTER (kind=iotk_CHARACTER1,len=ilen), allocatable :: dat_CHARACTER1 (:)
#endif

! Here is the rule:
! IF SOURCE IS BINARY: use the kind of source
! IF SOURCE IS TEXTUAL: use the default kind if available
!                       otherwise use the first kind found
! anyway, kind is computed runtime and in a future implementation
! it might be also asked to the user
!
  ierr=0
  select case(type(1:iotk_strlen(type)))
  case("LOGICAL")
    if(source_binary) then
      tmpkind=ikind
    else
      tmpkind=0
#ifdef __IOTK_LOGICAL1
      if(tmpkind==0) tmpkind=iotk_LOGICAL1
      if(iotk_LOGICAL1 == iotk_LOGICAL_defkind) then
        tmpkind=iotk_LOGICAL_defkind
      end if
#endif
#ifdef __IOTK_LOGICAL2
      if(tmpkind==0) tmpkind=iotk_LOGICAL2
      if(iotk_LOGICAL2 == iotk_LOGICAL_defkind) then
        tmpkind=iotk_LOGICAL_defkind
      end if
#endif
#ifdef __IOTK_LOGICAL3
      if(tmpkind==0) tmpkind=iotk_LOGICAL3
      if(iotk_LOGICAL3 == iotk_LOGICAL_defkind) then
        tmpkind=iotk_LOGICAL_defkind
      end if
#endif
#ifdef __IOTK_LOGICAL4
      if(tmpkind==0) tmpkind=iotk_LOGICAL4
      if(iotk_LOGICAL4 == iotk_LOGICAL_defkind) then
        tmpkind=iotk_LOGICAL_defkind
      end if
#endif
    end if
    select case(tmpkind)
#ifdef __IOTK_LOGICAL1
    case(iotk_LOGICAL1)
      allocate(dat_LOGICAL1(isize))
      call iotk_scan_dat_aux(source,dat_LOGICAL1,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_LOGICAL1,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_LOGICAL1)
#endif
#ifdef __IOTK_LOGICAL2
    case(iotk_LOGICAL2)
      allocate(dat_LOGICAL2(isize))
      call iotk_scan_dat_aux(source,dat_LOGICAL2,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_LOGICAL2,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_LOGICAL2)
#endif
#ifdef __IOTK_LOGICAL3
    case(iotk_LOGICAL3)
      allocate(dat_LOGICAL3(isize))
      call iotk_scan_dat_aux(source,dat_LOGICAL3,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_LOGICAL3,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_LOGICAL3)
#endif
#ifdef __IOTK_LOGICAL4
    case(iotk_LOGICAL4)
      allocate(dat_LOGICAL4(isize))
      call iotk_scan_dat_aux(source,dat_LOGICAL4,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_LOGICAL4,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_LOGICAL4)
#endif
    case default
      call iotk_error_issue(ierr,"iotk_copy_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'internal error')
    end select
  case("INTEGER")
    if(source_binary) then
      tmpkind=ikind
    else
      tmpkind=0
#ifdef __IOTK_INTEGER1
      if(tmpkind==0) tmpkind=iotk_INTEGER1
      if(iotk_INTEGER1 == iotk_INTEGER_defkind) then
        tmpkind=iotk_INTEGER_defkind
      end if
#endif
#ifdef __IOTK_INTEGER2
      if(tmpkind==0) tmpkind=iotk_INTEGER2
      if(iotk_INTEGER2 == iotk_INTEGER_defkind) then
        tmpkind=iotk_INTEGER_defkind
      end if
#endif
#ifdef __IOTK_INTEGER3
      if(tmpkind==0) tmpkind=iotk_INTEGER3
      if(iotk_INTEGER3 == iotk_INTEGER_defkind) then
        tmpkind=iotk_INTEGER_defkind
      end if
#endif
#ifdef __IOTK_INTEGER4
      if(tmpkind==0) tmpkind=iotk_INTEGER4
      if(iotk_INTEGER4 == iotk_INTEGER_defkind) then
        tmpkind=iotk_INTEGER_defkind
      end if
#endif
    end if
    select case(tmpkind)
#ifdef __IOTK_INTEGER1
    case(iotk_INTEGER1)
      allocate(dat_INTEGER1(isize))
      call iotk_scan_dat_aux(source,dat_INTEGER1,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_INTEGER1,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_INTEGER1)
#endif
#ifdef __IOTK_INTEGER2
    case(iotk_INTEGER2)
      allocate(dat_INTEGER2(isize))
      call iotk_scan_dat_aux(source,dat_INTEGER2,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_INTEGER2,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_INTEGER2)
#endif
#ifdef __IOTK_INTEGER3
    case(iotk_INTEGER3)
      allocate(dat_INTEGER3(isize))
      call iotk_scan_dat_aux(source,dat_INTEGER3,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_INTEGER3,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_INTEGER3)
#endif
#ifdef __IOTK_INTEGER4
    case(iotk_INTEGER4)
      allocate(dat_INTEGER4(isize))
      call iotk_scan_dat_aux(source,dat_INTEGER4,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_INTEGER4,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_INTEGER4)
#endif
    case default
      call iotk_error_issue(ierr,"iotk_copy_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'internal error')
    end select
  case("REAL")
    if(source_binary) then
      tmpkind=ikind
    else
      tmpkind=0
#ifdef __IOTK_REAL1
      if(tmpkind==0) tmpkind=iotk_REAL1
      if(iotk_REAL1 == iotk_REAL_defkind) then
        tmpkind=iotk_REAL_defkind
      end if
#endif
#ifdef __IOTK_REAL2
      if(tmpkind==0) tmpkind=iotk_REAL2
      if(iotk_REAL2 == iotk_REAL_defkind) then
        tmpkind=iotk_REAL_defkind
      end if
#endif
#ifdef __IOTK_REAL3
      if(tmpkind==0) tmpkind=iotk_REAL3
      if(iotk_REAL3 == iotk_REAL_defkind) then
        tmpkind=iotk_REAL_defkind
      end if
#endif
#ifdef __IOTK_REAL4
      if(tmpkind==0) tmpkind=iotk_REAL4
      if(iotk_REAL4 == iotk_REAL_defkind) then
        tmpkind=iotk_REAL_defkind
      end if
#endif
    end if
    select case(tmpkind)
#ifdef __IOTK_REAL1
    case(iotk_REAL1)
      allocate(dat_REAL1(isize))
      call iotk_scan_dat_aux(source,dat_REAL1,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_REAL1,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_REAL1)
#endif
#ifdef __IOTK_REAL2
    case(iotk_REAL2)
      allocate(dat_REAL2(isize))
      call iotk_scan_dat_aux(source,dat_REAL2,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_REAL2,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_REAL2)
#endif
#ifdef __IOTK_REAL3
    case(iotk_REAL3)
      allocate(dat_REAL3(isize))
      call iotk_scan_dat_aux(source,dat_REAL3,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_REAL3,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_REAL3)
#endif
#ifdef __IOTK_REAL4
    case(iotk_REAL4)
      allocate(dat_REAL4(isize))
      call iotk_scan_dat_aux(source,dat_REAL4,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_REAL4,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_REAL4)
#endif
    case default
      call iotk_error_issue(ierr,"iotk_copy_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'internal error')
    end select
  case("COMPLEX")
    if(source_binary) then
      tmpkind=ikind
    else
      tmpkind=0
#ifdef __IOTK_COMPLEX1
      if(tmpkind==0) tmpkind=iotk_COMPLEX1
      if(iotk_COMPLEX1 == iotk_COMPLEX_defkind) then
        tmpkind=iotk_COMPLEX_defkind
      end if
#endif
#ifdef __IOTK_COMPLEX2
      if(tmpkind==0) tmpkind=iotk_COMPLEX2
      if(iotk_COMPLEX2 == iotk_COMPLEX_defkind) then
        tmpkind=iotk_COMPLEX_defkind
      end if
#endif
#ifdef __IOTK_COMPLEX3
      if(tmpkind==0) tmpkind=iotk_COMPLEX3
      if(iotk_COMPLEX3 == iotk_COMPLEX_defkind) then
        tmpkind=iotk_COMPLEX_defkind
      end if
#endif
#ifdef __IOTK_COMPLEX4
      if(tmpkind==0) tmpkind=iotk_COMPLEX4
      if(iotk_COMPLEX4 == iotk_COMPLEX_defkind) then
        tmpkind=iotk_COMPLEX_defkind
      end if
#endif
    end if
    select case(tmpkind)
#ifdef __IOTK_COMPLEX1
    case(iotk_COMPLEX1)
      allocate(dat_COMPLEX1(isize))
      call iotk_scan_dat_aux(source,dat_COMPLEX1,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_COMPLEX1,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_COMPLEX1)
#endif
#ifdef __IOTK_COMPLEX2
    case(iotk_COMPLEX2)
      allocate(dat_COMPLEX2(isize))
      call iotk_scan_dat_aux(source,dat_COMPLEX2,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_COMPLEX2,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_COMPLEX2)
#endif
#ifdef __IOTK_COMPLEX3
    case(iotk_COMPLEX3)
      allocate(dat_COMPLEX3(isize))
      call iotk_scan_dat_aux(source,dat_COMPLEX3,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_COMPLEX3,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_COMPLEX3)
#endif
#ifdef __IOTK_COMPLEX4
    case(iotk_COMPLEX4)
      allocate(dat_COMPLEX4(isize))
      call iotk_scan_dat_aux(source,dat_COMPLEX4,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_COMPLEX4,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_COMPLEX4)
#endif
    case default
      call iotk_error_issue(ierr,"iotk_copy_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'internal error')
    end select
  case("CHARACTER")
    tmpkind=iotk_CHARACTER_defkind
    select case(tmpkind)
#ifdef __IOTK_CHARACTER1
    case(iotk_CHARACTER1)
      allocate(dat_CHARACTER1(isize))
      call iotk_scan_dat_aux(source,dat_CHARACTER1,ikind,ilen,fmt,ierr)
      if(ierr==0) call iotk_write_dat(dest,name,dat_CHARACTER1,attr=attr,ierr=ierr,fmt=fmt,columns=columns)
      deallocate(dat_CHARACTER1)
#endif
    case default
      call iotk_error_issue(ierr,"iotk_copy_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'internal error')
    end select
  case default
    call iotk_error_issue(ierr,"iotk_copy_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'internal error')
  end select
  
end subroutine iotk_copy_dat_aux_x


subroutine iotk_copy_dat_x(source,dest,source_binary,dest_binary,name,attr,maxsize,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_write_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,      intent(in)  :: source
  integer,      intent(in)  :: dest
  logical,      intent(in)  :: source_binary
  logical,      intent(in)  :: dest_binary
  character(*), intent(in)  :: name
  character(*), intent(in)  :: attr
  integer,      intent(in)  :: maxsize
  integer,      intent(out) :: ierr
  character(9) :: type
  integer :: ikind,isize,ilen,columns
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr1
  ierr = 0
  call iotk_parse_dat(attr,type,ikind,isize,ilen,fmt,columns,ierr)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_copy_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  if(iotk_strcomp(type,iotk_eos)) then
    call iotk_error_issue(ierr,"iotk_copy_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  if(isize==-1) then
    call iotk_error_issue(ierr,"iotk_copy_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  if(ilen==-1 .and. iotk_strcomp(type,"CHARACTER")) then
    call iotk_error_issue(ierr,"iotk_copy_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  if(isize<=maxsize .or. maxsize==-1 .or. dest_binary) then
    call iotk_copy_dat_aux(source,dest,source_binary,dest_binary,name,type,ikind,isize, &
                           ilen,fmt,columns,attr,ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_copy_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
      return
    end if  
  else    
    call iotk_strcpy(attr1,attr,ierr=ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_copy_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
      return
    end if
    call iotk_write_attr (attr1,"trunc",.true.,ierr=ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_copy_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
      return
    end if
    call iotk_write_empty(dest,name,attr=attr1,ierr=ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_copy_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
      return
    end if
  end if
end subroutine iotk_copy_dat_x

subroutine iotk_check_iotk_attr_x(unit,attr,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_xtox_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                 intent(in)  :: unit
  character(iotk_attlenx), intent(in)  :: attr
  integer,                 intent(out) :: ierr
  character(iotk_vallenx) :: file_version,extensions
  logical :: binary,rbinary,check,found,stream
  integer :: pos1,pos2,attlen,itmp_major,itmp_minor
  ierr = 0
  call iotk_scan_attr(attr,"file_version",file_version,eos=.true.,ierr=ierr,found=found)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  if(found) then
    attlen = iotk_strlen(file_version)
    pos1   = iotk_strscan(file_version,".")
    if(pos1<=1 .or. pos1>=attlen) then
      call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'Problems reading file version')
call iotk_error_write(ierr,"file_version",file_version)
call iotk_error_write(ierr,"attlen",attlen)
call iotk_error_write(ierr,"pos1",pos1)
      return
    end if
    pos2   = pos1 + verify(file_version(pos1+1:attlen),numbers)
    if(pos2==pos1+1) then
      call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'Problems reading file version')
call iotk_error_write(ierr,"file_version",file_version)
call iotk_error_write(ierr,"attlen",attlen)
call iotk_error_write(ierr,"pos1",pos1)
call iotk_error_write(ierr,"pos2",pos2)
      return
    end if
    if(pos2==pos1) pos2 = attlen+1
    call iotk_atoi(itmp_major,file_version(1:pos1-1),check)
    if(.not.check) then
      call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'Problems reading file version')
call iotk_error_write(ierr,"file_version",file_version)
      return
    end if
    call iotk_atoi(itmp_minor,file_version(pos1+1:pos2-1),check)
    if(.not.check) then
      call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'Problems reading file version')
call iotk_error_write(ierr,"file_version",file_version)
      return
    end if
    if(itmp_major > iotk_file_version_major .or. &
      (itmp_major==iotk_file_version_major .and. itmp_minor > iotk_file_version_minor) ) then
      call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'File version is newer than internal version')
call iotk_error_write(ierr,"file_version",file_version)
call iotk_error_write(ierr,"internal_version",iotk_file_version)
      return
    end if
  end if
  call iotk_scan_attr(attr,"binary",rbinary,ierr=ierr,found=found)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  if(found) then
    call iotk_inquire(unit,binary,stream,ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
      return
    end if
    if(rbinary .neqv. binary) then
      call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
      return
    end if
  end if
  call iotk_scan_attr(attr,"extensions",extensions,ierr=ierr,found=found,eos=.true.)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
    return
  end if
  if(found) then
    if(iotk_strlen(extensions) > 0) then
      call iotk_error_issue(ierr,"iotk_check_iotk_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'Extensions are not supported in this version')
call iotk_error_write(ierr,"extensions",extensions)
      return
    end if
  end if
end subroutine iotk_check_iotk_attr_x

function iotk_index_scal(index)
  use iotk_base
  use iotk_xtox_interf
  use iotk_misc_interf
  integer,           intent(in) :: index
  character(len=range(index)+3) :: iotk_index_scal
  iotk_index_scal="."//iotk_itoa(index)
end function iotk_index_scal
  
function iotk_index_vec(index)
  use iotk_base
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  integer,                         intent(in) :: index(:)
  character(len=(range(index)+3)*size(index)) :: iotk_index_vec
  integer :: length,i
  length = 0
  iotk_index_vec = " "
  do i = 1,size(index)
    iotk_index_vec(length+1:length+1+(range(index)+3)) = "."//iotk_itoa(index(i))
    length = len_trim(iotk_index_vec)
  end do
end function iotk_index_vec


subroutine iotk_tag_parse_x(tag,name,attr,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  character(iotk_taglenx), intent(in)  :: tag
  character(iotk_namlenx), intent(out) :: name
  character(iotk_attlenx), intent(out) :: attr
  integer,                 intent(out) :: ierr
  integer :: pos,lenatt,lentag
  ierr = 0
  lentag=iotk_strlen(tag)
  if(verify(tag(1:1),iotk_namcharfirst)/=0) then
    call iotk_error_issue(ierr,"iotk_tag_parse",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'Wrong syntax in tag')
    call iotk_error_write(ierr,"tag",tag(1:lentag))
    return
  end if
  pos = scan(tag(1:lentag)," ")
  if(pos==0) pos=lentag+1
  if(pos>len(name)+1) then
    call iotk_error_issue(ierr,"iotk_tag_parse",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'Tag name too long')
    return
   end if
  name = tag(1:pos-1)
  if(pos<=len(name)) name(pos:pos) = iotk_eos
  lenatt = len_trim(tag(pos:lentag))
  if(lenatt>iotk_attlenx) then
    call iotk_error_issue(ierr,"iotk_tag_parse",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'Attribute string too long')
    return
  end if
  if(lenatt>0) then
    attr(1:lenatt) = tag(pos:pos+lenatt-1)
    if(lenatt+1<=len(attr)) attr(lenatt+1:lenatt+1)=iotk_eos
  else
    attr(1:1)=iotk_eos
  end if
end subroutine iotk_tag_parse_x

function iotk_complete_filepath_x(newfile,oldfile)
  use iotk_base
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in) :: newfile
  character(len=*), intent(in) :: oldfile
  character(len=len(newfile)+len(oldfile)) :: iotk_complete_filepath_x
  character(len=len(oldfile)) :: prefix
  integer :: pos
  if(newfile(1:1)=="/") then
    iotk_complete_filepath_x = newfile
  else
    pos = scan(oldfile,"/",back=.true.)
    prefix = " "
    if(pos>0) prefix = oldfile(1:pos)
    iotk_complete_filepath_x = trim(prefix)//trim(newfile)
  end if
end function iotk_complete_filepath_x

function iotk_check_name_x(name)
  use iotk_base
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  character(len=*), intent(in) :: name
  logical                      :: iotk_check_name_x
! Checks a single name
  integer :: len_name
  iotk_check_name_x = .true.
  len_name = iotk_strlen_trim(name)
  if(len_name>iotk_namlenx) iotk_check_name_x = .false.
  if(verify(name(1:1),iotk_namcharfirst)/=0) iotk_check_name_x = .false.
  if(len_name>1) then
    if(verify(name(2:len_name),iotk_namchar)/=0) iotk_check_name_x = .false.
  end if
end function iotk_check_name_x

subroutine iotk_delete_attr_x(attr,name,ierr)
  use iotk_base
  use iotk_str_interf
  use iotk_error_interf
  implicit none
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  integer,          intent(out)   :: ierr
  integer :: attlen,pos,equal,begin
  logical :: foundl
  character :: delim
  ierr = 0
  attlen=iotk_strlen(attr)
  foundl = .false.
  equal = 0
  begin = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierr,"iotk_delete_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'')
call iotk_error_write(ierr,"attr",attr)
call iotk_error_write(ierr,"equal",equal)
      return
    end if
    if(trim(attr(equal:equal+pos-1))==trim(name)) foundl = .true.
    begin = equal
    equal = equal + pos
    pos   = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierr,"iotk_delete_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
      return
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierr,"iotk_delete_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
call iotk_error_msg(ierr,'delim')
      return
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierr,"iotk_delete_attr",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.39 ")
      return
    end if
    equal = equal + pos
    if(foundl) exit
  end do
  if(foundl) then
    if(equal<attlen) then
      pos = verify(attr(equal+1:attlen)," ")
      if(pos==0) then
        equal=attlen
      else
        equal=equal+pos-1
      end if
    end if
    if(equal<attlen) then
      attr(begin:begin+attlen-equal) = attr(equal+1:attlen)//iotk_eos
    else
      attr(begin:begin)=iotk_eos
    end if
  end if
end subroutine iotk_delete_attr_x

subroutine iotk_readcmdline_x(args,nargs,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  implicit none
  character(len=*),  intent(out) :: args(:)
  integer,           intent(out) :: nargs
  logical, optional, intent(in)  :: eos
  integer, optional, intent(out) :: ierr
  integer :: iarg,iostat,pos,ierrl
  logical :: leos

  ierrl = 0
  iarg = 0
  leos=.false.
  nargs = 0
  if(present(eos))leos=eos
  do
    iarg = iarg + 1
    if(iarg>size(args)) then
      call iotk_error_issue(ierrl,"iotk_delete_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
      goto 1
    end if
    read(*,"(a)",iostat=iostat) args(iarg)
    if(iostat<0) exit
    if(iostat>0) then
      call iotk_error_issue(ierrl,"iotk_delete_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.39 ")
      goto 1
    end if
    pos = scan(args(iarg),"|",back=.true.)
    if(pos>0) then
      if(args(iarg)(pos:)=="|") then
        args(iarg)(pos:pos) = iotk_eos
      end if
    else
      pos = len_trim(args(iarg)) + 1
      if(pos<=len(args)) args(iarg)(pos:pos) = iotk_eos
    end if
    if(.not. leos) then
      pos = iotk_strlen(args(iarg))
      if(pos<len(args(iarg))) args(iarg)(pos+1:) = " "
    end if
  end do
  iarg = iarg -1
  nargs = iarg
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_readcmdline_x

subroutine iotk_init_static_vars_x()
  use iotk_base
  implicit none
  logical, save :: first=.true.
  if(first) then
    iotk_error_pool_used     = .false.
    iotk_error_pool_order    = 0
    iotk_unitmin             = iotk_unitmin_def
    iotk_unitmax             = iotk_unitmax_def
    iotk_error_unit          = iotk_error_unit_def
    iotk_output_unit         = iotk_output_unit_def
    iotk_getline_buffer      = iotk_getline_buffer_def
    iotk_linlen              = iotk_linlen_def
    iotk_indent              = iotk_indent_def
    iotk_maxindent           = iotk_maxindent_def
    iotk_error_warn_overflow = iotk_error_warn_overflow_def
  end if
  first=.false.
end subroutine iotk_init_static_vars_x


