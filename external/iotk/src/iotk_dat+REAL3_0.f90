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


#ifdef __IOTK_REAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                    intent(in)  :: unit
  character(len=*),                           intent(in)  :: name
  REAL (kind=this_kind),           intent(in)  :: dat  
  type(iotk_dummytype),             optional              :: dummy
  character(len=*),                 optional, intent(in)  :: attr
  integer,                          optional, intent(in)  :: columns
  character(len=*),                 optional, intent(in)  :: sep
  character(len=*),                 optional, intent(in)  :: fmt
  integer,                          optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw,stream
  integer :: lcolumns
!-<
  integer :: dat_rank,i
  integer, dimension(:), allocatable :: dat_shape
  logical :: qe_syntax
!->
  integer(iotk_header_kind), parameter :: idummy=0
  character(100) :: lsep
  character(300) :: usefmt
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
!-<
  character(len=10) :: tmpstr
!->
  type (iotk_unit), pointer :: this
  REAL (kind=this_kind),allocatable :: dattmp(:)
  integer :: itmp
  ierrl = 0
  iostat = 0
  lcolumns = 1
  lsep(1:2) = " "//iotk_eos
  if(present(columns)) lcolumns = columns
  if(present(sep)) then
    call iotk_strcpy(lsep,sep,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
!-<  
  qe_syntax = .false.
  if (associated(this)) then
     qe_syntax = this%qe_syntax
  end if
!->
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary,stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
call iotk_error_write(ierrl,"unit",unit)
call iotk_error_write(ierrl,"name",trim(name))
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  call iotk_write_attr(lattr,"type",iotk_tolower("REAL"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",1,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  end if
!->
!-<
  if (qe_syntax) then
    lattr(1:1) = iotk_eos
    dat_rank = size(shape(dat))
    if (dat_rank>0) then
       allocate(dat_shape(dat_rank))
       dat_shape = shape(dat)
       call iotk_write_attr(lattr,"rank",dat_rank,ierr=ierrl,first=.true.)
       if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
          goto 1
       end if
       do i=1,dat_rank
       write(tmpstr,'(i10)') i
  tmpstr = adjustl(tmpstr)
          call iotk_write_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),dat_shape(i),ierr=ierrl)
          if(ierrl/=0) then
             call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
             goto 1
          end if
       end do
       deallocate(dat_shape)
    end if
  end if
!->
  if(binary) then
    call iotk_write_attr(lattr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
  if(.not.iotk_strcomp(usefmt,"!") .and. .not.qe_syntax) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(lcolumns/=1) call iotk_write_attr(lattr,"columns",lcolumns,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"columns",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  else
    if(present(attr)) then
      call iotk_write_begin(unit,name,attr,ierr=ierrl)
    else
      call iotk_write_begin(unit,name,ierr=ierrl)
    end if
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_write_begin(unit,iotk_tolower("REAL"),lattr,ierr=ierrl,new_line=.false.)
  end if
!->
  allocate(dattmp(1))
     dattmp(1) = dat
  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
  else
    if(raw) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
!-<
    if ((.not. qe_syntax) .or. (size(dattmp)>1)) then
     write(lunit,fmt=trim(iotk_wfmt("REAL",kind(dattmp),lcolumns,-1,lsep)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
     end if
    else
     write(lunit,fmt=trim(iotk_wfmt("REAL",kind(dattmp),lcolumns,-1,lsep)),advance='no',iostat=iostat) dattmp(1)
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
     end if
    end if
!->
    else
!-<
    if ((.not. qe_syntax) .or. (size(dattmp)>1)) then
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),advance='no',iostat=iostat) dattmp(1)
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
!->
    end if
  end if
!-<
  if (qe_syntax) then
    call iotk_write_end(unit,iotk_tolower("REAL"),indentation=.false.,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
!->
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_REAL3_0



subroutine iotk_scan_dat_REAL3_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
!-<
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
!->
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
  REAL (kind=this_kind),              allocatable :: tmpdat(:)
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
!-<
! ... necessary because i need the syntax to use
  type (iotk_unit), pointer :: this
! ... necessary to read the tag describing the type
  character(iotk_taglenx) :: ltag
  character(iotk_namlenx) :: r_name
  character(iotk_attlenx) :: lattr2
  character(len=20) :: tmpstr
! ... necessary for scan_tag
  logical :: binary,stream,qe_syntax
  integer :: r_control,rrank
  integer :: rshape,i
!->
  integer :: columns
  logical :: inside,foundl
!-<
  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)
  qe_syntax = .false.
  if (associated(this)) THEN
     qe_syntax = this%qe_syntax
  end if    
!->
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,lattr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  if(present(attr)) call iotk_strcpy(attr,lattr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!-<
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,columns,ierrl)
! Note that columns is not effectively used
  if(ierrl/=0) goto 1
!-<
  else
     call iotk_inquire(iotk_phys_unit(unit),binary,stream,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
     end if
     do
        call iotk_scan_tag(unit,+1,r_control,ltag,binary,stream,ierrl)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           goto 1
        end if
       if (r_control==4) then
   cycle
else
   exit
end if
     end do

     call iotk_tag_parse(ltag,r_name,lattr2,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
goto 1
     end if
     call iotk_strcpy(rtype,r_name,ierrl)
     rtype = iotk_toupper(rtype)

      rlen = -1

     if (rtype(1:iotk_strlen(rtype))/="STRING") then
        call iotk_scan_attr(lattr2,"rank",rrank,ierr=ierrl,default=0)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           return
        end if
        rsize = 1
        if (rrank>0) then
           do i=1,rrank
           write(tmpstr,'(i20)') i
      tmpstr = adjustl(tmpstr)
      call iotk_scan_attr(lattr2,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
              if(ierrl/=0) then
                 call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
                 return
              end if
      rsize = rsize*rshape
   end do
end if
     else
        call iotk_strcpy(rtype,"CHARACTER",ierrl)
        rsize = -1        
     end if

     call iotk_scan_attr(lattr2,"kind",rkind,ierr=ierrl,default=-1)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
     if(ierrl/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"columns",columns,ierr=ierrl,default=1)
     if(ierrl/=0) then
       call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
       return
     end if
  end if
!->
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"REAL") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,' ')
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
call iotk_error_write(ierrl,"type","REAL")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==1) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
  allocate(tmpdat(1))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Error reading data')
call iotk_error_write(ierrl,"name",name)
call iotk_error_write(ierrl,"rkind",rkind)
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if
     dat = tmpdat(1)
1 continue
!-<
  if ( allocated(tmpdat) ) deallocate(tmpdat)
!->
  if(inside) then
!-<
    if (qe_syntax) then 
      call iotk_scan_end(unit,iotk_tolower("REAL"),ierr=ierrl2) 
       if(ierrl2/=0) then
          call iotk_error_clear(ierrl)
          ierrl=ierrl2
       end if
    end if
!->
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Dat not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_REAL3_0

!-<
! ... new procedure: scan dat when you are already inside the data tag
subroutine iotk_scan_dat_inside_REAL3_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif

  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
  REAL (kind=this_kind),              allocatable :: tmpdat(:)
  character(len=20) :: tmpstr
  integer :: rrank
  integer :: rshape,i
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  character(iotk_namlenx) :: rname
  type (iotk_unit), pointer :: this
  integer :: columns
  logical :: inside,foundl,qe_syntax


  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)

  qe_syntax = .true.
  IF (associated(this)) then
     qe_syntax = this%qe_syntax
  END IF

  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.

  if (.not.qe_syntax) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if

  rname = iotk_tolower("REAL")
  call iotk_scan_begin(unit,iotk_tolower("REAL"),lattr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1


  foundl = .true.
  inside = .true.
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  
  rlen = -1


     call iotk_scan_attr(lattr,"rank",rrank,ierr=ierrl,default=0)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     rsize = 1
     if (rrank>0) then
        do i=1,rrank
   write(tmpstr,'(i20)') i
   tmpstr = adjustl(tmpstr)
           call iotk_scan_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
           if(ierrl/=0) then
              call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
              return
           end if
           rsize = rsize*rshape
        end do
     end if


  call iotk_scan_attr(lattr,"kind",rkind,ierr=ierrl,default=-1)
  if(ierrl/=0) then
     call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
  if(ierrl/=0) then
     call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"columns",columns,ierr=ierrl,default=1)
  if(ierrl/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
    return
  end if

  if(.not. (rsize==-1 .or. rsize==1) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)

  allocate(tmpdat(1))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Error reading data')
call iotk_error_write(ierrl,"rname",rname)
call iotk_error_write(ierrl,"rkind",rkind)
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if

     dat = tmpdat(1)
  deallocate(tmpdat)
1 continue
  if(inside) then
     call iotk_scan_end(unit,iotk_tolower("REAL"),ierr=ierrl2) 
     if(ierrl2/=0) then
        call iotk_error_clear(ierrl)
        ierrl=ierrl2
     end if
  end if

  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Dat not found')
call iotk_error_write(ierrl,"rname",rname)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if

end subroutine iotk_scan_dat_inside_REAL3_0
!->

#endif
#endif

subroutine iotk_dat_dummy_REAL3_0
  write(0,*)
end subroutine iotk_dat_dummy_REAL3_0





!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

#include "iotk_auxmacros.h"


#ifdef __IOTK_REAL3
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                    intent(in)  :: unit
  character(len=*),                           intent(in)  :: name
  REAL (kind=this_kind),           intent(in)  :: dat (:) 
  type(iotk_dummytype),             optional              :: dummy
  character(len=*),                 optional, intent(in)  :: attr
  integer,                          optional, intent(in)  :: columns
  character(len=*),                 optional, intent(in)  :: sep
  character(len=*),                 optional, intent(in)  :: fmt
  integer,                          optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw,stream
  integer :: lcolumns
!-<
  integer :: dat_rank,i
  integer, dimension(:), allocatable :: dat_shape
  logical :: qe_syntax
!->
  integer(iotk_header_kind), parameter :: idummy=0
  character(100) :: lsep
  character(300) :: usefmt
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
!-<
  character(len=10) :: tmpstr
!->
  type (iotk_unit), pointer :: this
  REAL (kind=this_kind),allocatable :: dattmp(:)
  integer :: itmp
  ierrl = 0
  iostat = 0
  lcolumns = 1
  lsep(1:2) = " "//iotk_eos
  if(present(columns)) lcolumns = columns
  if(present(sep)) then
    call iotk_strcpy(lsep,sep,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
!-<  
  qe_syntax = .false.
  if (associated(this)) then
     qe_syntax = this%qe_syntax
  end if
!->
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary,stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
call iotk_error_write(ierrl,"unit",unit)
call iotk_error_write(ierrl,"name",trim(name))
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  call iotk_write_attr(lattr,"type",iotk_tolower("REAL"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  end if
!->
!-<
  if (qe_syntax) then
    lattr(1:1) = iotk_eos
    dat_rank = size(shape(dat))
    if (dat_rank>0) then
       allocate(dat_shape(dat_rank))
       dat_shape = shape(dat)
       call iotk_write_attr(lattr,"rank",dat_rank,ierr=ierrl,first=.true.)
       if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
          goto 1
       end if
       do i=1,dat_rank
       write(tmpstr,'(i10)') i
  tmpstr = adjustl(tmpstr)
          call iotk_write_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),dat_shape(i),ierr=ierrl)
          if(ierrl/=0) then
             call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
             goto 1
          end if
       end do
       deallocate(dat_shape)
    end if
  end if
!->
  if(binary) then
    call iotk_write_attr(lattr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
  if(.not.iotk_strcomp(usefmt,"!") .and. .not.qe_syntax) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(lcolumns/=1) call iotk_write_attr(lattr,"columns",lcolumns,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"columns",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  else
    if(present(attr)) then
      call iotk_write_begin(unit,name,attr,ierr=ierrl)
    else
      call iotk_write_begin(unit,name,ierr=ierrl)
    end if
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_write_begin(unit,iotk_tolower("REAL"),lattr,ierr=ierrl,new_line=.true.)
  end if
!->
  allocate(dattmp(size(dat)))
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
     call iotk_private_pack_REAL3(dattmp,dat,size(dattmp),1)
#else
     dattmp = pack(dat,mask=.true.)
#endif
  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
  else
    if(raw) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
!-<
    if ((.not. qe_syntax) .or. (size(dattmp)>1)) then
     write(lunit,fmt=trim(iotk_wfmt("REAL",kind(dattmp),lcolumns,-1,lsep)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
     end if
    else
     write(lunit,fmt=trim(iotk_wfmt("REAL",kind(dattmp),lcolumns,-1,lsep)),advance='no',iostat=iostat) dattmp(1)
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
     end if
    end if
!->
    else
!-<
    if ((.not. qe_syntax) .or. (size(dattmp)>1)) then
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),advance='no',iostat=iostat) dattmp(1)
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
!->
    end if
  end if
!-<
  if (qe_syntax) then
    call iotk_write_end(unit,iotk_tolower("REAL"),indentation=.true.,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
!->
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_REAL3_1


recursive subroutine iotk_scan_dat_aux_REAL3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only: iotk_read
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  use iotk_stream_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)              :: dat (:)
#else
  REAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
  integer(iotk_header_kind) :: idummy
  logical :: raw,binary,stream
  integer :: lunit
  integer :: i
#ifdef __IOTK_WORKAROUND3
  integer :: j
#endif
  integer :: index,length,nexttag,iostat,altlength
  type(iotk_unit), pointer :: this
  character(len=iotk_linlenx) :: line,altline
#ifdef __IOTK_REAL1
  REAL (kind=iotk_REAL1), allocatable :: dat1 (:)
#endif
#ifdef __IOTK_REAL2
  REAL (kind=iotk_REAL2), allocatable :: dat2 (:)
#endif
#ifdef __IOTK_REAL4
  REAL (kind=iotk_REAL4), allocatable :: dat4 (:)
#endif
  lunit = iotk_phys_unit(unit)
  ierr = 0
  iostat = 0
  idummy = 0
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary,stream,ierr)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
    return
  end if
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
#ifdef __IOTK_WORKAROUND3
        read(lunit,iostat=iostat) ( dat(j), j=1,ubound(dat,1) )
#else
        read(lunit,iostat=iostat) dat
#endif
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        if(stream) then
          call iotk_stream_read(lunit,idummy,dat,ierr=ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
            return
          end if
        else
#ifdef __IOTK_WORKAROUND3
          read(lunit,iostat=iostat) idummy, ( dat(j), j=1,ubound(dat,1) )
#else
          read(lunit,iostat=iostat) idummy, dat
#endif
        end if
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
#ifdef __IOTK_REAL1
    case(kind(dat1))
      ! for the sake of completeness: if the file is raw, there are no
      ! information about kind and this line cannot be reached
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
        return
      end if
      allocate(dat1(ubound(dat,1)))
      if(stream) then
        call iotk_stream_read(lunit,idummy,dat1,ierr=ierr)
        if(ierr/=0) then
          deallocate(dat1)
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,( dat1(i), i=1,ubound(dat1,1) )
        if(iostat/=0) then
          deallocate(dat1)
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
      dat = real(dat1,kind=kind(dat))
      deallocate(dat1)
#endif
#ifdef __IOTK_REAL2
    case(kind(dat2))
      ! for the sake of completeness: if the file is raw, there are no
      ! information about kind and this line cannot be reached
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
        return
      end if
      allocate(dat2(ubound(dat,1)))
      if(stream) then
        call iotk_stream_read(lunit,idummy,dat2,ierr=ierr)
        if(ierr/=0) then
          deallocate(dat2)
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,( dat2(i), i=1,ubound(dat2,1) )
        if(iostat/=0) then
          deallocate(dat2)
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
      dat = real(dat2,kind=kind(dat))
      deallocate(dat2)
#endif
#ifdef __IOTK_REAL4
    case(kind(dat4))
      ! for the sake of completeness: if the file is raw, there are no
      ! information about kind and this line cannot be reached
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
        return
      end if
      allocate(dat4(ubound(dat,1)))
      if(stream) then
        call iotk_stream_read(lunit,idummy,dat4,ierr=ierr)
        if(ierr/=0) then
          deallocate(dat4)
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,( dat4(i), i=1,ubound(dat4,1) )
        if(iostat/=0) then
          deallocate(dat4)
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
      dat = real(dat4,kind=kind(dat))
      deallocate(dat4)
#endif
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,'Kind incompatibility')
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(raw) then
#ifdef __IOTK_WORKAROUND3
      read(lunit,fmt=*,iostat=iostat) ( dat(j), j=1,ubound(dat,1) )
#else
      read(lunit,fmt=*,iostat=iostat) dat
#endif
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"*")) then
#ifdef __IOTK_WORKAROUND3
      read(lunit,fmt=*,iostat=iostat) ( dat(j), j=1,ubound(dat,1) )
#else
      read(lunit,fmt=*,iostat=iostat) dat
#endif
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
          return
        end if
        nexttag = scan(line(1:length),"<")
        if(nexttag==0) then
          nexttag = length + 1
        else
! adjust the positioning if there is a tag on this line
! implementation to be improved
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_str_clean(line(1:nexttag - 1))
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,'Error reading REAL data')
          return
        end if
        if(index == size(dat)) exit
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
          return
        end if
      end do
    else
#ifdef __IOTK_WORKAROUND3
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) ( dat(j), j=1,ubound(dat,1) )
#else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
#endif
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
    return
  end if
end subroutine iotk_scan_dat_aux_REAL3

subroutine iotk_scan_dat_REAL3_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
!-<
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
!->
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
  REAL (kind=this_kind),              allocatable :: tmpdat(:)
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
!-<
! ... necessary because i need the syntax to use
  type (iotk_unit), pointer :: this
! ... necessary to read the tag describing the type
  character(iotk_taglenx) :: ltag
  character(iotk_namlenx) :: r_name
  character(iotk_attlenx) :: lattr2
  character(len=20) :: tmpstr
! ... necessary for scan_tag
  logical :: binary,stream,qe_syntax
  integer :: r_control,rrank
  integer :: rshape,i
!->
  integer :: columns
  logical :: inside,foundl
!-<
  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)
  qe_syntax = .false.
  if (associated(this)) THEN
     qe_syntax = this%qe_syntax
  end if    
!->
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,lattr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  if(present(attr)) call iotk_strcpy(attr,lattr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!-<
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,columns,ierrl)
! Note that columns is not effectively used
  if(ierrl/=0) goto 1
!-<
  else
     call iotk_inquire(iotk_phys_unit(unit),binary,stream,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
     end if
     do
        call iotk_scan_tag(unit,+1,r_control,ltag,binary,stream,ierrl)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           goto 1
        end if
       if (r_control==4) then
   cycle
else
   exit
end if
     end do

     call iotk_tag_parse(ltag,r_name,lattr2,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
goto 1
     end if
     call iotk_strcpy(rtype,r_name,ierrl)
     rtype = iotk_toupper(rtype)

      rlen = -1

     if (rtype(1:iotk_strlen(rtype))/="STRING") then
        call iotk_scan_attr(lattr2,"rank",rrank,ierr=ierrl,default=0)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           return
        end if
        rsize = 1
        if (rrank>0) then
           do i=1,rrank
           write(tmpstr,'(i20)') i
      tmpstr = adjustl(tmpstr)
      call iotk_scan_attr(lattr2,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
              if(ierrl/=0) then
                 call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
                 return
              end if
      rsize = rsize*rshape
   end do
end if
     else
        call iotk_strcpy(rtype,"CHARACTER",ierrl)
        rsize = -1        
     end if

     call iotk_scan_attr(lattr2,"kind",rkind,ierr=ierrl,default=-1)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
     if(ierrl/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"columns",columns,ierr=ierrl,default=1)
     if(ierrl/=0) then
       call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
       return
     end if
  end if
!->
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"REAL") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,' ')
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
call iotk_error_write(ierrl,"type","REAL")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Error reading data')
call iotk_error_write(ierrl,"name",name)
call iotk_error_write(ierrl,"rkind",rkind)
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
     call iotk_private_pack_REAL3(dat,tmpdat,size(tmpdat),1)
#else
     dat = reshape(tmpdat,shape(dat))
#endif
1 continue
!-<
  if ( allocated(tmpdat) ) deallocate(tmpdat)
!->
  if(inside) then
!-<
    if (qe_syntax) then 
      call iotk_scan_end(unit,iotk_tolower("REAL"),ierr=ierrl2) 
       if(ierrl2/=0) then
          call iotk_error_clear(ierrl)
          ierrl=ierrl2
       end if
    end if
!->
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Dat not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_REAL3_1

!-<
! ... new procedure: scan dat when you are already inside the data tag
subroutine iotk_scan_dat_inside_REAL3_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif

  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
  REAL (kind=this_kind),              allocatable :: tmpdat(:)
  character(len=20) :: tmpstr
  integer :: rrank
  integer :: rshape,i
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  character(iotk_namlenx) :: rname
  type (iotk_unit), pointer :: this
  integer :: columns
  logical :: inside,foundl,qe_syntax


  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)

  qe_syntax = .true.
  IF (associated(this)) then
     qe_syntax = this%qe_syntax
  END IF

  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.

  if (.not.qe_syntax) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if

  rname = iotk_tolower("REAL")
  call iotk_scan_begin(unit,iotk_tolower("REAL"),lattr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1


  foundl = .true.
  inside = .true.
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  
  rlen = -1


     call iotk_scan_attr(lattr,"rank",rrank,ierr=ierrl,default=0)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     rsize = 1
     if (rrank>0) then
        do i=1,rrank
   write(tmpstr,'(i20)') i
   tmpstr = adjustl(tmpstr)
           call iotk_scan_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
           if(ierrl/=0) then
              call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
              return
           end if
           rsize = rsize*rshape
        end do
     end if


  call iotk_scan_attr(lattr,"kind",rkind,ierr=ierrl,default=-1)
  if(ierrl/=0) then
     call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
  if(ierrl/=0) then
     call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"columns",columns,ierr=ierrl,default=1)
  if(ierrl/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
    return
  end if

  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Error reading data')
call iotk_error_write(ierrl,"rname",rname)
call iotk_error_write(ierrl,"rkind",rkind)
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if

#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
     call iotk_private_pack_REAL3(dat,tmpdat,size(tmpdat),1)
#else
     dat = reshape(tmpdat,shape(dat))
#endif
  deallocate(tmpdat)
1 continue
  if(inside) then
     call iotk_scan_end(unit,iotk_tolower("REAL"),ierr=ierrl2) 
     if(ierrl2/=0) then
        call iotk_error_clear(ierrl)
        ierrl=ierrl2
     end if
  end if

  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Dat not found')
call iotk_error_write(ierrl,"rname",rname)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if

end subroutine iotk_scan_dat_inside_REAL3_1
!->

#endif
#endif

subroutine iotk_dat_dummy_REAL3_1
  write(0,*)
end subroutine iotk_dat_dummy_REAL3_1





!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

#include "iotk_auxmacros.h"


#ifdef __IOTK_REAL3
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                    intent(in)  :: unit
  character(len=*),                           intent(in)  :: name
  REAL (kind=this_kind),           intent(in)  :: dat (:,:) 
  type(iotk_dummytype),             optional              :: dummy
  character(len=*),                 optional, intent(in)  :: attr
  integer,                          optional, intent(in)  :: columns
  character(len=*),                 optional, intent(in)  :: sep
  character(len=*),                 optional, intent(in)  :: fmt
  integer,                          optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw,stream
  integer :: lcolumns
!-<
  integer :: dat_rank,i
  integer, dimension(:), allocatable :: dat_shape
  logical :: qe_syntax
!->
  integer(iotk_header_kind), parameter :: idummy=0
  character(100) :: lsep
  character(300) :: usefmt
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
!-<
  character(len=10) :: tmpstr
!->
  type (iotk_unit), pointer :: this
  REAL (kind=this_kind),allocatable :: dattmp(:)
  integer :: itmp
  ierrl = 0
  iostat = 0
  lcolumns = 1
  lsep(1:2) = " "//iotk_eos
  if(present(columns)) lcolumns = columns
  if(present(sep)) then
    call iotk_strcpy(lsep,sep,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
!-<  
  qe_syntax = .false.
  if (associated(this)) then
     qe_syntax = this%qe_syntax
  end if
!->
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary,stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
call iotk_error_write(ierrl,"unit",unit)
call iotk_error_write(ierrl,"name",trim(name))
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  call iotk_write_attr(lattr,"type",iotk_tolower("REAL"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  end if
!->
!-<
  if (qe_syntax) then
    lattr(1:1) = iotk_eos
    dat_rank = size(shape(dat))
    if (dat_rank>0) then
       allocate(dat_shape(dat_rank))
       dat_shape = shape(dat)
       call iotk_write_attr(lattr,"rank",dat_rank,ierr=ierrl,first=.true.)
       if(ierrl/=0) then
          call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
          goto 1
       end if
       do i=1,dat_rank
       write(tmpstr,'(i10)') i
  tmpstr = adjustl(tmpstr)
          call iotk_write_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),dat_shape(i),ierr=ierrl)
          if(ierrl/=0) then
             call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
             goto 1
          end if
       end do
       deallocate(dat_shape)
    end if
  end if
!->
  if(binary) then
    call iotk_write_attr(lattr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
  if(.not.iotk_strcomp(usefmt,"!") .and. .not.qe_syntax) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(lcolumns/=1) call iotk_write_attr(lattr,"columns",lcolumns,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!->
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"columns",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  else
    if(present(attr)) then
      call iotk_write_begin(unit,name,attr,ierr=ierrl)
    else
      call iotk_write_begin(unit,name,ierr=ierrl)
    end if
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
    call iotk_write_begin(unit,iotk_tolower("REAL"),lattr,ierr=ierrl,new_line=.true.)
  end if
!->
  allocate(dattmp(size(dat)))
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
     call iotk_private_pack_REAL3(dattmp,dat,size(dattmp),1)
#else
     dattmp = pack(dat,mask=.true.)
#endif
  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
  else
    if(raw) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
!-<
    if ((.not. qe_syntax) .or. (size(dattmp)>1)) then
     write(lunit,fmt=trim(iotk_wfmt("REAL",kind(dattmp),lcolumns,-1,lsep)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
     end if
    else
     write(lunit,fmt=trim(iotk_wfmt("REAL",kind(dattmp),lcolumns,-1,lsep)),advance='no',iostat=iostat) dattmp(1)
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
     end if
    end if
!->
    else
!-<
    if ((.not. qe_syntax) .or. (size(dattmp)>1)) then
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),advance='no',iostat=iostat) dattmp(1)
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
      end if
    end if
!->
    end if
  end if
!-<
  if (qe_syntax) then
    call iotk_write_end(unit,iotk_tolower("REAL"),indentation=.true.,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
      goto 1
    end if
  end if
!->
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_REAL3_2



subroutine iotk_scan_dat_REAL3_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
!-<
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
!->
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
  REAL (kind=this_kind),              allocatable :: tmpdat(:)
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
!-<
! ... necessary because i need the syntax to use
  type (iotk_unit), pointer :: this
! ... necessary to read the tag describing the type
  character(iotk_taglenx) :: ltag
  character(iotk_namlenx) :: r_name
  character(iotk_attlenx) :: lattr2
  character(len=20) :: tmpstr
! ... necessary for scan_tag
  logical :: binary,stream,qe_syntax
  integer :: r_control,rrank
  integer :: rshape,i
!->
  integer :: columns
  logical :: inside,foundl
!-<
  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)
  qe_syntax = .false.
  if (associated(this)) THEN
     qe_syntax = this%qe_syntax
  end if    
!->
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,lattr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  if(present(attr)) call iotk_strcpy(attr,lattr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
!-<
  if (.not.qe_syntax) then
!-<
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,columns,ierrl)
! Note that columns is not effectively used
  if(ierrl/=0) goto 1
!-<
  else
     call iotk_inquire(iotk_phys_unit(unit),binary,stream,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        goto 1
     end if
     do
        call iotk_scan_tag(unit,+1,r_control,ltag,binary,stream,ierrl)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           goto 1
        end if
       if (r_control==4) then
   cycle
else
   exit
end if
     end do

     call iotk_tag_parse(ltag,r_name,lattr2,ierrl)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
goto 1
     end if
     call iotk_strcpy(rtype,r_name,ierrl)
     rtype = iotk_toupper(rtype)

      rlen = -1

     if (rtype(1:iotk_strlen(rtype))/="STRING") then
        call iotk_scan_attr(lattr2,"rank",rrank,ierr=ierrl,default=0)
        if(ierrl/=0) then
           call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
           return
        end if
        rsize = 1
        if (rrank>0) then
           do i=1,rrank
           write(tmpstr,'(i20)') i
      tmpstr = adjustl(tmpstr)
      call iotk_scan_attr(lattr2,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
              if(ierrl/=0) then
                 call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
                 return
              end if
      rsize = rsize*rshape
   end do
end if
     else
        call iotk_strcpy(rtype,"CHARACTER",ierrl)
        rsize = -1        
     end if

     call iotk_scan_attr(lattr2,"kind",rkind,ierr=ierrl,default=-1)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
     if(ierrl/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
        return
     end if
     call iotk_scan_attr(lattr2,"columns",columns,ierr=ierrl,default=1)
     if(ierrl/=0) then
       call iotk_error_issue(ierr,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
       return
     end if
  end if
!->
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"REAL") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,' ')
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
call iotk_error_write(ierrl,"type","REAL")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Error reading data')
call iotk_error_write(ierrl,"name",name)
call iotk_error_write(ierrl,"rkind",rkind)
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
     call iotk_private_pack_REAL3(dat,tmpdat,size(tmpdat),1)
#else
     dat = reshape(tmpdat,shape(dat))
#endif
1 continue
!-<
  if ( allocated(tmpdat) ) deallocate(tmpdat)
!->
  if(inside) then
!-<
    if (qe_syntax) then 
      call iotk_scan_end(unit,iotk_tolower("REAL"),ierr=ierrl2) 
       if(ierrl2/=0) then
          call iotk_error_clear(ierrl)
          ierrl=ierrl2
       end if
    end if
!->
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Dat not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_REAL3_2

!-<
! ... new procedure: scan dat when you are already inside the data tag
subroutine iotk_scan_dat_inside_REAL3_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  use iotk_unit_interf
  use iotk_attr_interf, only : iotk_scan_attr
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif

  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
  REAL (kind=this_kind),              allocatable :: tmpdat(:)
  character(len=20) :: tmpstr
  integer :: rrank
  integer :: rshape,i
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  character(iotk_namlenx) :: rname
  type (iotk_unit), pointer :: this
  integer :: columns
  logical :: inside,foundl,qe_syntax


  call iotk_unit_get(iotk_phys_unit(unit),pointer=this)

  qe_syntax = .true.
  IF (associated(this)) then
     qe_syntax = this%qe_syntax
  END IF

  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.

  if (.not.qe_syntax) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if

  rname = iotk_tolower("REAL")
  call iotk_scan_begin(unit,iotk_tolower("REAL"),lattr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1


  foundl = .true.
  inside = .true.
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  
  rlen = -1


     call iotk_scan_attr(lattr,"rank",rrank,ierr=ierrl,default=0)
     if(ierrl/=0) then
        call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
        return
     end if
     rsize = 1
     if (rrank>0) then
        do i=1,rrank
   write(tmpstr,'(i20)') i
   tmpstr = adjustl(tmpstr)
           call iotk_scan_attr(lattr,"n"//tmpstr(1:len_trim(tmpstr)),rshape,ierr=ierrl,default=0)
           if(ierrl/=0) then
              call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
              return
           end if
           rsize = rsize*rshape
        end do
     end if


  call iotk_scan_attr(lattr,"kind",rkind,ierr=ierrl,default=-1)
  if(ierrl/=0) then
     call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"fmt", fmt, ierr=ierrl,eos=.true.,default="!"//iotk_eos)
  if(ierrl/=0) then
     call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
     return
  end if
  call iotk_scan_attr(lattr,"columns",columns,ierr=ierrl,default=1)
  if(ierrl/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.27 ")
    return
  end if

  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Error reading data')
call iotk_error_write(ierrl,"rname",rname)
call iotk_error_write(ierrl,"rkind",rkind)
call iotk_error_write(ierrl,"rlen",rlen)
    goto 1
  end if

#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
     call iotk_private_pack_REAL3(dat,tmpdat,size(tmpdat),1)
#else
     dat = reshape(tmpdat,shape(dat))
#endif
  deallocate(tmpdat)
1 continue
  if(inside) then
     call iotk_scan_end(unit,iotk_tolower("REAL"),ierr=ierrl2) 
     if(ierrl2/=0) then
        call iotk_error_clear(ierrl)
        ierrl=ierrl2
     end if
  end if

  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_dat_inside",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.27 ")
call iotk_error_msg(ierrl,'Dat not found')
call iotk_error_write(ierrl,"rname",rname)
    ierrl = - ierrl
  end if 
  if(present(default) .and. .not. foundl) then
    dat=default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if

end subroutine iotk_scan_dat_inside_REAL3_2
!->

#endif
#endif

subroutine iotk_dat_dummy_REAL3_2
  write(0,*)
end subroutine iotk_dat_dummy_REAL3_2



