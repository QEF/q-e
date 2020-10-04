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


subroutine iotk_copyfile_x(dummy,source,dest,source_unit,dest_unit,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  type(iotk_dummytype), optional         :: dummy
  character(len=*), optional, intent(in) :: source
  character(len=*), optional, intent(in) :: dest
  integer,          optional, intent(in) :: source_unit
  integer,          optional, intent(in) :: dest_unit
  integer,          optional, intent(out):: ierr
  integer :: ierrl,unit1,unit2
  integer :: iostat,length
  character(len=iotk_linlenx) :: line
  iostat = 0
  ierrl  = 0
  if(present(source) .eqv. present(source_unit)) then
    call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Use exactly one between source and source_unit')
    goto 1
  end if
  if(present(dest)   .eqv. present(dest_unit)) then
    call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Use exactly one between dest and dest_unit')
    goto 1
  end if
  if(present(source)) then
    call iotk_free_unit(unit1,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error searching for a free unit')
      goto 1
    end if
    open(unit1,file=trim(iotk_strpad(source)),iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'messaggio')
call iotk_error_write(ierrl,"sourcefile",trim(iotk_strpad(source)))
call iotk_error_write(ierrl,"sourceunit",unit1)
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  else
    unit1=source_unit
  end if
  if(present(dest)) then
    call iotk_free_unit(unit2,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
      goto 1
    end if
    open(unit2,file=trim(iotk_strpad(dest)),iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error opening destination file')
call iotk_error_write(ierrl,"destfile",trim(iotk_strpad(dest)))
call iotk_error_write(ierrl,"destunit",unit2)
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  else
    unit2=dest_unit
  end if
  do
    call iotk_getline(unit1,line,length,ierrl)
    if(ierrl/=0) then
      call iotk_error_scan(ierrl,"iostat",iostat)
      if(iostat<0) then
        call iotk_error_clear(ierrl)
        exit
      end if
      call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error reading source file')
call iotk_error_write(ierrl,"sourceunit",unit1)
      goto 1
    end if
    write(unit2,"(a)",iostat=iostat) line(1:length)
    if(iostat/=0) then
       call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing destination file')
call iotk_error_write(ierrl,"destunit",unit2)
call iotk_error_write(ierrl,"iostat",iostat)
       goto 1
    end if 
  end do
  iostat=0
  if(present(source)) then
    close(unit1,iostat=iostat)
    if(iostat/=0) then
       call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error closing source file')
call iotk_error_write(ierrl,"sourcefile",trim(iotk_strpad(source)))
call iotk_error_write(ierrl,"sourceunit",unit1)
call iotk_error_write(ierrl,"iostat",iostat)
       goto 1
    end if 
  end if
  if(present(dest)) then
    close(unit2,iostat=iostat)
    if(iostat/=0) then
       call iotk_error_issue(ierrl,"iotk_copyfile_x",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error closing destination file')
call iotk_error_write(ierrl,"destfile",trim(iotk_strpad(dest)))
call iotk_error_write(ierrl,"destunit",unit2)
call iotk_error_write(ierrl,"iostat",iostat)
       goto 1
    end if
  end if
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_copyfile_x

subroutine iotk_link_x(unit,name,file,dummy,binary,raw,create,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_files_interf
  use iotk_str_interf
  use iotk_write_interf
  use iotk_misc_interf
  use iotk_unit_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*),           intent(in)  :: file
  type(iotk_dummytype), optional      :: dummy
  logical,      optional, intent(in)  :: binary
  logical,      optional, intent(in)  :: raw
  logical,      optional, intent(in)  :: create
  integer,      optional, intent(out) :: ierr
  logical :: lbinary,lraw,lcreate
  integer :: ierrl,iostat
  integer :: lunit,link_unit
  type(iotk_unit), pointer :: this_unit
  character(iotk_attlenx) :: attr
  character(iotk_fillenx) :: oldfile
  ierrl  = 0
  iostat = 0
  lbinary=.false.
  lraw   =.false.
  lcreate=.false.
  if(present(binary)) lbinary = binary
  if(present(raw))    lraw = raw
  if(present(create)) lcreate = create
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this_unit)
  if(.not.associated(this_unit)) then
    call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Links do not apply to units which are not explicitly connected')
    goto 1
  end if
  call iotk_write_attr(attr,"iotk_link",iotk_strtrim(file),ierr=ierrl,first=.true.)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  if(lraw) then
    if(lbinary) then
      call iotk_write_attr(attr,"iotk_binary",lbinary,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
        goto 1
      end if
    end if
    call iotk_write_attr(attr,"iotk_raw",lraw,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
      goto 1
    end if
  end if
  call iotk_write_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  call iotk_write_comment(unit,"This is a link to the file indicated in the iotk_link attribute",ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  call iotk_write_end  (unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  if(lcreate) then
    call iotk_free_unit(link_unit,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
      goto 1
    end if
    inquire(unit=lunit,name=oldfile,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error inquiring')
call iotk_error_write(ierrl,"unit",lunit)
call iotk_error_write(ierrl,"file",trim(oldfile))
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
    call iotk_open_write(link_unit,file=iotk_complete_filepath(file,trim(oldfile)), &
!-<
! ... inside link we use the same syntax (qe or not)
                                 binary=lbinary,raw=lraw,qe_syntax=this_unit%qe_syntax,skip_root=.true.,ierr=ierrl)
!->
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
      goto 1
    end if
    call iotk_unit_parent(parent=lunit,son=link_unit,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_link",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
      goto 1
    end if
  end if
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_link_x

!-<
subroutine iotk_open_write_x(unit,file,dummy,attr,binary,new,raw,root,qe_syntax,skip_root,skip_head,ierr)
!->
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_write_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*), optional, intent(in)  :: file
  type(iotk_dummytype), optional      :: dummy
  character(*), optional, intent(in)  :: attr
  logical,      optional, intent(in)  :: binary
  logical,      optional, intent(in)  :: new
  logical,      optional, intent(in)  :: raw
  character(*), optional, intent(in)  :: root
!-<
  logical,      optional, intent(in)  :: qe_syntax
!->
  logical,      optional, intent(in)  :: skip_root
  logical,      optional, intent(in)  :: skip_head
  integer,      optional, intent(out) :: ierr
! Opens a file properly
  integer :: iostat
  character(50) :: status,form
  character(iotk_namlenx) :: lroot
  character(iotk_attlenx) :: lattr
  integer :: ierrl
  logical :: lbinary,lraw,lnew,lskip_root,lskip_head,lstream
!-<
  logical :: lqe_syntax
!->
  type (iotk_unit), pointer :: this
  ierrl = 0
  iostat = 0
  lroot = "Root"
  lraw = .false.
  lnew = .false.
!-<
  lqe_syntax = .false.
!->
  lbinary = .false.
  lskip_root = .false.
  lskip_head = .false.
  if(present(root)) lroot = root
  if(present(raw)) lraw=raw
  if(present(binary)) lbinary = binary
  if(present(new)) lnew = new
!-<
  if(present(qe_syntax)) lqe_syntax = qe_syntax
!->
  if(present(skip_root)) lskip_root = skip_root
  if(lskip_root) lroot=""
  if(present(skip_head)) lskip_head = skip_head
  if(present(file)) then
    form = "formatted"
    if(lbinary) form = "unformatted"
    status = "unknown"
    if(lnew) status = "new"
    open(unit=unit,file=file,status=status,form=form,position="rewind",iostat=iostat,action="write")
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error opening file')
call iotk_error_write(ierrl,"unit",unit)
call iotk_error_write(ierrl,"file",file)
call iotk_error_write(ierrl,"binary",lbinary)
call iotk_error_write(ierrl,"new",lnew)
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  else
    call iotk_inquire(unit,lbinary,lstream,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
      goto 1
    end if
  end if
  if(.not.lraw) then
    if(.not.lskip_head) then
      if(.not. lbinary) then
        write(unit,"(a)",iostat=iostat) '<?xml version="1.0"?>'
        if(iostat/=0) then
          call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing XML tag')
call iotk_error_write(ierrl,"unit",unit)
call iotk_error_write(ierrl,"iostat",iostat)
          goto 1
        end if
      end if
      call iotk_write_attr(lattr,"version",trim(iotk_version),first=.true.,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing version attribute')
        goto 1
      end if
      call iotk_write_pi(unit,"iotk",lattr,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing version tag')
        goto 1
      end if
      call iotk_write_attr(lattr,"file_version",trim(iotk_file_version),first=.true.,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing file_version attribute')
        goto 1
      end if
      call iotk_write_pi(unit,"iotk",lattr,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing version tag')
        goto 1
      end if
      call iotk_write_attr(lattr,"binary",lbinary,first=.true.,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing binary attribute')
        goto 1
      end if
      call iotk_write_pi(unit,"iotk",lattr,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing binary tag')
        goto 1
      end if
!-<
      call iotk_write_attr(lattr,"qe_syntax",lqe_syntax,first=.true.,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing qe_syntax attribute')
        goto 1
      end if
      call iotk_write_pi(unit,"iotk",lattr,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing qe_syntax tag')
        goto 1
      end if
!->
    end if
    if(.not.lskip_root) then
      lattr(1:1) = iotk_eos
      if(present(attr)) then
        call iotk_strcpy(lattr,attr,ierr=ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing attributes from the root tag')
          goto 1
        end if
      end if
      call iotk_write_begin(unit,lroot,attr=lattr,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error writing the root tag')
        goto 1
      end if
    end if
  end if
  call iotk_unit_add(unit,this,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_open_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'Error adding the unit to the list')
    goto 1
  end if
  this%root=lroot
  this%raw=lraw
!-<
  this%qe_syntax = lqe_syntax
!->
  this%close_at_end=present(file)
  this%skip_root=lskip_root
  if(lskip_root) this%level = -1
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_open_write_x

recursive subroutine iotk_close_write_x(unit,dummy,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_write_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  type(iotk_dummytype), optional      :: dummy
  integer,      optional, intent(out) :: ierr
! Closes a file properly
  logical :: binary,stream
  integer :: ierrl,iostat
  type(iotk_unit), pointer :: this
  nullify(this)
  ierrl = 0
  iostat = 0
  call iotk_unit_get(unit,pointer=this)
  if(.not.associated(this)) then
    call iotk_error_issue(ierrl,"iotk_close_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  call iotk_inquire(unit,binary,stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_close_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  if(.not.this%raw) then
    if(.not.this%skip_root) then
      call iotk_write_end(unit,this%root,ierr=ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_close_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
        goto 1
      end if
    end if
  end if
  if(this%close_at_end) then
    if(.not.binary) then
      write(unit,*,iostat=iostat)
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_close_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'unit')
call iotk_error_write(ierrl,"iostat",iostat)
        goto 1
      end if
    end if
    close(unit,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_close_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'unit')
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
  call iotk_unit_del(unit,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_close_write",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_close_write_x


!-<
subroutine iotk_open_read_x(unit,file,dummy,attr,binary,stream,raw,qe_syntax,root,ierr)
!->
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_attr_interf
  use iotk_scan_interf
  use iotk_unit_interf
  use iotk_misc_interf
  use iotk_files_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*), optional, intent(in)  :: file
  type(iotk_dummytype), optional      :: dummy
  logical,      optional, intent(in)  :: binary
!-<
  logical,      optional, intent(in)  :: qe_syntax
!->
  logical,      optional, intent(in)  :: stream
  logical,      optional, intent(in)  :: raw
#ifdef __IOTK_WORKAROUND6
  character(len=*), optional              :: attr
  character(len=*), optional              :: root
#else
  character(len=*), optional, intent(out) :: attr
  character(len=*), optional, intent(out) :: root
#endif
  integer,      optional, intent(out) :: ierr
  character(50)           :: form,access
  character(iotk_attlenx) :: lattr
  character(iotk_taglenx) :: tag
  character(iotk_namlenx) :: lroot
  type(iotk_unit),pointer :: this
  integer                 :: ierrl,control,iostat
  logical                 :: lbinary,lraw,lstream
!-<
  logical                 :: lqe_syntax,found
!->
  ierrl = 0
  iostat = 0
  lbinary=.false.
!-<
  lqe_syntax=.false.
!->
  lstream=.false.
  lraw=.false.
  lroot = " "
  lattr(1:1) = iotk_eos
  if(present(raw)) lraw=raw
  if(present(file)) then
    if(present(binary)) lbinary = binary
    if(present(stream)) lstream = stream
    if(.not.lbinary .and. .not. lraw) call iotk_magic(file,lbinary)
    form = "formatted"
    if(lbinary) form = "unformatted"
    access="sequential"
#ifdef __IOTK_STREAMS
    if(lstream .and. lbinary) access= "stream"
#endif
    open(unit=unit,file=trim(file(1:iotk_strlen(file))),status="old",form=form,position="rewind",iostat=iostat,action="read", &
         access=access)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'unit')
call iotk_error_write(ierrl,"file",file)
call iotk_error_write(ierrl,"binary",lbinary)
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  else
    call iotk_inquire(unit,lbinary,lstream,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
      goto 1
    end if
  end if
  if(.not.lraw) then
    do
      call iotk_scan_tag(unit,+1,control,tag,lbinary,lstream,ierrl)
      if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
        goto 1
      end if
      select case(control)
      case(1)
        call iotk_tag_parse(tag,lroot,lattr,ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
          goto 1
        end if
        exit
      case(2:3)
        call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'End or empty tag at the beginning of a file')
call iotk_error_write(ierrl,"unit",unit)
call iotk_error_write(ierrl,"file",trim(file(1:iotk_strlen(file))))
call iotk_error_write(ierrl,"binary",lbinary)
call iotk_error_write(ierrl,"iostat",iostat)
        goto 1
      case(5)
        call iotk_tag_parse(tag,lroot,lattr,ierrl)
        if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
          goto 1
        end if
        if(iotk_strcomp(lroot,"iotk")) then
          call iotk_check_iotk_attr(unit,lattr,ierrl)
          if(ierrl/=0) then
            call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
            goto 1
          end if
!-<
  call iotk_scan_attr(lattr,"qe_syntax",lqe_syntax,default=.false.,ierr=ierrl,found=found)
          if(ierrl/=0) then
            call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
            goto 1
          end if
!->
        end if
      end select
    end do
  end if
!-<
  if(present(qe_syntax)) lqe_syntax = qe_syntax
!->
  if(present(root)) root = lroot
  if(present(attr)) call iotk_strcpy(attr,lattr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  call iotk_unit_add(unit,this,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_open_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  this%root=lroot
  this%raw=lraw
  this%close_at_end=present(file)
  this%skip_root=.false.
!-<
  this%qe_syntax = lqe_syntax
!->
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_open_read_x

subroutine iotk_close_read_x(unit,dummy,ierr) 
  use iotk_base
  use iotk_error_interf
  use iotk_scan_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  type(iotk_dummytype), optional      :: dummy
  integer,      optional, intent(out) :: ierr
  integer                 :: ierrl
  integer :: iostat
  type(iotk_unit), pointer :: this
  character(iotk_namlenx) :: root
  logical :: raw
  logical :: close_at_end
  ierrl = 0
  iostat = 0
  call iotk_unit_get(unit,pointer=this)
  if(.not.associated(this)) then
    call iotk_error_issue(ierrl,"iotk_close_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  root = this%root
  close_at_end = this%close_at_end
  raw = this%raw
  call iotk_unit_del(unit,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_close_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
    goto 1
  end if
  if(.not.raw) then      
    call iotk_scan_end(unit,root,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_close_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
      goto 1
    end if
  end if
  if(close_at_end) then
    close(unit,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_close_read",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.20 ")
call iotk_error_msg(ierrl,'unit')
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
  end if
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_close_read_x

subroutine iotk_magic_x(file,binary)
  use iotk_base
  use iotk_str_interf
  use iotk_error_interf
  use iotk_scan_interf
  use iotk_misc_interf
  use iotk_unit_interf
  use iotk_attr_interf
  character(len=*), intent(in) :: file
  logical,          intent(out):: binary
  integer :: iostat,unit,control,ierrl
  logical :: found,opened
  character(len=iotk_taglenx) :: tag
  character(len=iotk_namlenx) :: name
  character(len=iotk_attlenx) :: attr
  ierrl=0
  binary=.false.
  call iotk_free_unit(unit)
  open(unit=unit,file=trim(file(1:iotk_strlen(file))),status="old",form="unformatted", &
       position="rewind",iostat=iostat,action="read")
  if(iostat/=0) goto 1
  do
    call iotk_scan_tag(unit,+1,control,tag,.true.,.false.,ierrl)
    if(ierrl/=0) goto 1
    if(control==1) then
      exit
    else if(control==5) then
      call iotk_tag_parse(tag,name,attr,ierrl)
      if(iotk_strcomp(name,"iotk")) then
        call iotk_scan_attr(attr,"binary",binary,found=found,ierr=ierrl)
        if(ierrl/=0) goto 1
        if(found) goto 1
      end if
    end if
  end do
1 continue
  if(ierrl/=0) call iotk_error_clear(ierrl)
  inquire(unit=unit,opened=opened)
  if(opened) close(unit,iostat=iostat)
end subroutine iotk_magic_x
