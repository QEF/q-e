# 47 "iotk_attr.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_attr.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_CHARACTER1
#if 0 <= __IOTK_MAXRANK

# 64 "iotk_attr.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_CHARACTER1(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 70 "iotk_attr.spp"
    CHARACTER (kind=__IOTK_CHARACTER1,len=l), intent(out) :: out(n)
    CHARACTER (kind=__IOTK_CHARACTER1,len=l), intent(in)  :: in(n)
# 76 "iotk_attr.spp"
    out = in
end subroutine iotk_private_pack_CHARACTER1

# 137 "iotk_attr.spp"

# 210 "iotk_attr.spp"

# 213 "iotk_attr.spp"
subroutine iotk_write_attr_CHARACTER1_0(attr,name,val,dummy,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  CHARACTER(kind=__IOTK_CHARACTER1,len=*), intent(in)  :: val 
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 233 "iotk_attr.spp"
  logical :: lquot,lapos
# 235 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 242 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 242 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 242 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 249 "iotk_attr.spp"
  lquot=iotk_strscan(val,'"')>0
  lapos=iotk_strscan(val,"'")>0
  if(.not.lquot) then
    delim='"'
    call iotk_deescape(tmpval,val)
  else if(.not.lapos) then
    delim="'"
    call iotk_deescape(tmpval,val)
  else
    delim='"'
    call iotk_deescape(tmpval,val,quot=.true.,apos=.true.)
  end if
# 273 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 275 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 275 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  attr(attlen+1:attlen+vallen+namlen+5) = " "//trim(name)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_CHARACTER1_0

# 288 "iotk_attr.spp"
subroutine iotk_scan_attr_CHARACTER1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=__IOTK_CHARACTER1,len=*)                        :: val 
#else
  CHARACTER(kind=__IOTK_CHARACTER1,len=*), intent(out)           :: val 
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  CHARACTER(kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 313 "iotk_attr.spp"
  character(iotk_vallenx) :: valctmp
  integer :: vallen,defaultlen
  logical :: leos
# 320 "iotk_attr.spp"
  ierrl = 0
  attlen=iotk_strlen(attr)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 330 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 330 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 330 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 337 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 343 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 348 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 357 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
      goto 1
    end if
  else
    goto 1
  end if
# 364 "iotk_attr.spp"
  call iotk_escape(valctmp,valc)
  vallen = iotk_strlen(valctmp)
  if(len(val) < vallen) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 367 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
    goto 1
  end if
  leos=.false.
  if(present(eos)) leos=eos
  val(1:vallen) = valctmp(1:vallen)
  if(len(val) > vallen) then
    val(vallen+1:vallen+1) = iotk_eos
    if(.not.leos) then
      val(vallen+1:)=" "
    end if
  end if
# 402 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 406 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.7 ")
# 406 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 406 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 411 "iotk_attr.spp"
    if(leos) then
      defaultlen = min(iotk_strlen(default),len(val))
      val(1:defaultlen) = default(1:defaultlen)
      if(defaultlen<len(val)) val(defaultlen+1:defaultlen+1)=iotk_eos
    else
      val = default
    end if
# 421 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_CHARACTER1_0
# 429 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_CHARACTER1_0
  write(0,*)
end subroutine iotk_attr_dummy_CHARACTER1_0

# 47 "iotk_dat.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_CHARACTER1
#if 0 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_CHARACTER1_0(unit,name,dat,dummy,attr,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat  
  type(iotk_dummytype), optional      :: dummy
  character(len=*), optional, intent(in)  :: attr
  character(len=*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 93 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 104 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 109 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 114 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(lattr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 123 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",1,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 128 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 132 "iotk_dat.spp"
  call iotk_write_attr(lattr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 146 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 148 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 155 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 160 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 165 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 170 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 180 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 186 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if

  allocate(dattmp(1))
# 197 "iotk_dat.spp"
     dattmp(1) = dat
# 209 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 227 "iotk_dat.spp"
      write(lunit,"(a)",iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 231 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 238 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 243 "iotk_dat.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
        end if
     end do
# 258 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 261 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 268 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_CHARACTER1_0


# 545 "iotk_dat.spp"

# 547 "iotk_dat.spp"
subroutine iotk_scan_dat_CHARACTER1_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat 
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat 
#endif
  type(iotk_dummytype), optional         :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),optional              :: attr
#else
  character(len=*),optional, intent(out) :: attr
#endif
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 572 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 576 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  logical :: inside,foundl
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
# 592 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==1) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 602 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 607 "iotk_dat.spp"
  if(rlen ==-1) rlen  = len(dat)
# 609 "iotk_dat.spp"

  allocate(tmpdat(1))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 613 "iotk_dat.spp"
        dat = tmpdat(1)
# 617 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
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
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 629 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_CHARACTER1_0


#endif
#endif

subroutine iotk_dat_dummy_CHARACTER1_0
  write(0,*)
end subroutine iotk_dat_dummy_CHARACTER1_0

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_CHARACTER1
#if 1 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_CHARACTER1_1(unit,name,dat,dummy,attr,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:) 
  type(iotk_dummytype), optional      :: dummy
  character(len=*), optional, intent(in)  :: attr
  character(len=*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 93 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 104 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 109 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 114 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(lattr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 123 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 128 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 132 "iotk_dat.spp"
  call iotk_write_attr(lattr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 146 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 148 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 155 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 160 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 165 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 170 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 180 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 186 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 199 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 201 "iotk_dat.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 205 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 209 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 227 "iotk_dat.spp"
      write(lunit,"(a)",iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 231 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 238 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 243 "iotk_dat.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
        end if
     end do
# 258 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 261 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 268 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_CHARACTER1_1


# 283 "iotk_dat.spp"
recursive subroutine iotk_scan_dat_aux_CHARACTER1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only: iotk_read
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat (:)
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:)
#endif
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
  integer(iotk_header_kind) :: idummy
  logical :: raw,binary
  integer :: lunit
  integer :: index,length,nexttag,iostat,altlength
  type(iotk_unit), pointer :: this
  character(len=iotk_linlenx) :: line,altline
# 309 "iotk_dat.spp"
  CHARACTER (kind=kind(dat), len=rlen) :: dattmp(ubound(dat,1))
# 319 "iotk_dat.spp"
  lunit = iotk_phys_unit(unit)
  ierr = 0
  iostat = 0
  idummy = 0
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(unit=lunit,binary=binary,ierr=ierr)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 330 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
    return
  end if
# 334 "iotk_dat.spp"
  if(binary) then
    if(raw) then
      read(lunit,iostat=iostat) dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 338 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 338 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 338 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else
      read(lunit,iostat=iostat) idummy,dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 344 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 344 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 344 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  else
    if(raw) then
      read(lunit,"(a)",iostat=iostat) dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 352 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 352 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 352 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,  iostat=iostat) dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 358 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 358 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 358 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      iostat = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 367 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
          return
        end if
        nexttag = scan(line(1:length),"<")
        if(nexttag==0) then
          nexttag = length + 1
        else
! AGGIUSTA LA POSIZIONE SE C'E' UNA TAG SU QUESTA LINEA
! E' UN PO' CASERECCIO MA FUNZIONA
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 378 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 378 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 378 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 383 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 388 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 388 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 388 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 393 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 393 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 393 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        index = index + 1
        call iotk_escape(to=dattmp(index),from=line(1:nexttag - 1))
        if(iotk_strlen(dattmp(index)) < len(dattmp)) dattmp(index)(iotk_strlen(dattmp(index))+1:) = " "
        if(index == size(dat)) exit
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 402 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 402 "iotk_dat.spp"
call iotk_error_msg(ierr,'Missing dat')
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 409 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
# 409 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 409 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
  if(len(dattmp) <= len(dat)) then
    dat (:) = dattmp (:)
  else
    dat (:) = dattmp (:) (1:len(dat))
  end if
# 539 "iotk_dat.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 540 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.12 ")
    return
  end if
end subroutine iotk_scan_dat_aux_CHARACTER1
# 545 "iotk_dat.spp"

# 547 "iotk_dat.spp"
subroutine iotk_scan_dat_CHARACTER1_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat (:)
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype), optional         :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),optional              :: attr
#else
  character(len=*),optional, intent(out) :: attr
#endif
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 572 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 576 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  logical :: inside,foundl
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
# 592 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 602 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 607 "iotk_dat.spp"
  if(rlen ==-1) rlen  = len(dat)
# 609 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 615 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 617 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
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
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 629 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_CHARACTER1_1


#endif
#endif

subroutine iotk_dat_dummy_CHARACTER1_1
  write(0,*)
end subroutine iotk_dat_dummy_CHARACTER1_1

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_CHARACTER1
#if 2 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_CHARACTER1_2(unit,name,dat,dummy,attr,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:) 
  type(iotk_dummytype), optional      :: dummy
  character(len=*), optional, intent(in)  :: attr
  character(len=*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 93 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 104 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 109 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 114 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(lattr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 123 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 128 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 132 "iotk_dat.spp"
  call iotk_write_attr(lattr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 146 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 148 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 155 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 160 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 165 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 170 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 180 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 186 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 199 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 201 "iotk_dat.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 205 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 209 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 227 "iotk_dat.spp"
      write(lunit,"(a)",iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 231 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 238 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 243 "iotk_dat.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
        end if
     end do
# 258 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 261 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 268 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_CHARACTER1_2


# 545 "iotk_dat.spp"

# 547 "iotk_dat.spp"
subroutine iotk_scan_dat_CHARACTER1_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat (:,:)
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype), optional         :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),optional              :: attr
#else
  character(len=*),optional, intent(out) :: attr
#endif
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 572 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 576 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  logical :: inside,foundl
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
# 592 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 602 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 607 "iotk_dat.spp"
  if(rlen ==-1) rlen  = len(dat)
# 609 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 615 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 617 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
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
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 629 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_CHARACTER1_2


#endif
#endif

subroutine iotk_dat_dummy_CHARACTER1_2
  write(0,*)
end subroutine iotk_dat_dummy_CHARACTER1_2

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_CHARACTER1
#if 3 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_CHARACTER1_3(unit,name,dat,dummy,attr,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:) 
  type(iotk_dummytype), optional      :: dummy
  character(len=*), optional, intent(in)  :: attr
  character(len=*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 93 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 104 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 109 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 114 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(lattr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 123 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 128 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 132 "iotk_dat.spp"
  call iotk_write_attr(lattr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 146 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 148 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 155 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 160 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 165 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 170 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 180 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 186 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 199 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 201 "iotk_dat.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 205 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 209 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 227 "iotk_dat.spp"
      write(lunit,"(a)",iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 231 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 238 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 243 "iotk_dat.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
        end if
     end do
# 258 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 261 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 268 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_CHARACTER1_3


# 545 "iotk_dat.spp"

# 547 "iotk_dat.spp"
subroutine iotk_scan_dat_CHARACTER1_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat (:,:,:)
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype), optional         :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),optional              :: attr
#else
  character(len=*),optional, intent(out) :: attr
#endif
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 572 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 576 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  logical :: inside,foundl
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
# 592 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 602 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 607 "iotk_dat.spp"
  if(rlen ==-1) rlen  = len(dat)
# 609 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 615 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 617 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
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
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 629 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_CHARACTER1_3


#endif
#endif

subroutine iotk_dat_dummy_CHARACTER1_3
  write(0,*)
end subroutine iotk_dat_dummy_CHARACTER1_3

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_CHARACTER1
#if 4 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_CHARACTER1_4(unit,name,dat,dummy,attr,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:,:) 
  type(iotk_dummytype), optional      :: dummy
  character(len=*), optional, intent(in)  :: attr
  character(len=*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 93 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 104 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 109 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 114 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(lattr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 123 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 128 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 132 "iotk_dat.spp"
  call iotk_write_attr(lattr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 146 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 148 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 155 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 160 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 165 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 170 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 180 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 186 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 199 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 201 "iotk_dat.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 205 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 209 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 227 "iotk_dat.spp"
      write(lunit,"(a)",iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 231 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 238 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 243 "iotk_dat.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
        end if
     end do
# 258 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 261 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 268 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_CHARACTER1_4


# 545 "iotk_dat.spp"

# 547 "iotk_dat.spp"
subroutine iotk_scan_dat_CHARACTER1_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat (:,:,:,:)
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype), optional         :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),optional              :: attr
#else
  character(len=*),optional, intent(out) :: attr
#endif
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 572 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 576 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  logical :: inside,foundl
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
# 592 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 602 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 607 "iotk_dat.spp"
  if(rlen ==-1) rlen  = len(dat)
# 609 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 615 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 617 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
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
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 629 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_CHARACTER1_4


#endif
#endif

subroutine iotk_dat_dummy_CHARACTER1_4
  write(0,*)
end subroutine iotk_dat_dummy_CHARACTER1_4

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_CHARACTER1
#if 5 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_CHARACTER1_5(unit,name,dat,dummy,attr,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:,:,:) 
  type(iotk_dummytype), optional      :: dummy
  character(len=*), optional, intent(in)  :: attr
  character(len=*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 93 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 104 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 109 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 114 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(lattr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 123 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 128 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 132 "iotk_dat.spp"
  call iotk_write_attr(lattr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 146 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 148 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 155 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 160 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 165 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 170 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 180 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 186 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 199 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 201 "iotk_dat.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 205 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 209 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 227 "iotk_dat.spp"
      write(lunit,"(a)",iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 231 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 238 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 243 "iotk_dat.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
        end if
     end do
# 258 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 261 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 268 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_CHARACTER1_5


# 545 "iotk_dat.spp"

# 547 "iotk_dat.spp"
subroutine iotk_scan_dat_CHARACTER1_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat (:,:,:,:,:)
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype), optional         :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),optional              :: attr
#else
  character(len=*),optional, intent(out) :: attr
#endif
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 572 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 576 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  logical :: inside,foundl
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
# 592 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 602 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 607 "iotk_dat.spp"
  if(rlen ==-1) rlen  = len(dat)
# 609 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 615 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 617 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
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
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 629 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_CHARACTER1_5


#endif
#endif

subroutine iotk_dat_dummy_CHARACTER1_5
  write(0,*)
end subroutine iotk_dat_dummy_CHARACTER1_5

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_CHARACTER1
#if 6 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_CHARACTER1_6(unit,name,dat,dummy,attr,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:,:,:,:) 
  type(iotk_dummytype), optional      :: dummy
  character(len=*), optional, intent(in)  :: attr
  character(len=*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 93 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 104 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 109 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 114 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(lattr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 123 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 128 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 132 "iotk_dat.spp"
  call iotk_write_attr(lattr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 146 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 148 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 155 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 160 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 165 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 170 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 180 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 186 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 199 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 201 "iotk_dat.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 205 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 209 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 227 "iotk_dat.spp"
      write(lunit,"(a)",iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 231 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 238 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 243 "iotk_dat.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
        end if
     end do
# 258 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 261 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 268 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_CHARACTER1_6


# 545 "iotk_dat.spp"

# 547 "iotk_dat.spp"
subroutine iotk_scan_dat_CHARACTER1_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat (:,:,:,:,:,:)
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype), optional         :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),optional              :: attr
#else
  character(len=*),optional, intent(out) :: attr
#endif
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 572 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 576 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  logical :: inside,foundl
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
# 592 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 602 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 607 "iotk_dat.spp"
  if(rlen ==-1) rlen  = len(dat)
# 609 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 615 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 617 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
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
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 629 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_CHARACTER1_6


#endif
#endif

subroutine iotk_dat_dummy_CHARACTER1_6
  write(0,*)
end subroutine iotk_dat_dummy_CHARACTER1_6

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 55 "iotk_dat.spp"
#include "iotk_auxmacros.h"
# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_CHARACTER1
#if 7 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_CHARACTER1_7(unit,name,dat,dummy,attr,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_write_attr
  use iotk_write_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:,:,:,:,:) 
  type(iotk_dummytype), optional      :: dummy
  character(len=*), optional, intent(in)  :: attr
  character(len=*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: lattr
  character(iotk_attlenx) :: attr_tmp
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 93 "iotk_dat.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,pointer=this)
  raw = .false.
  if(associated(this)) then
    raw = this%raw
  end if
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 104 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 109 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 114 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 118 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 118 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(lattr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 123 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_attr(lattr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 128 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 132 "iotk_dat.spp"
  call iotk_write_attr(lattr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 134 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
# 146 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(lattr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 148 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(present(attr)) then
    attr_tmp(1:1)=iotk_eos
    call iotk_strcpy(attr_tmp,attr,ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 155 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"type",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 160 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"kind",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 165 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"size",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 170 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"fmt",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    call iotk_delete_attr(attr_tmp,"len",ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 180 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
      goto 1
    end if
    if(iotk_strlen_trim(attr_tmp)>0) call iotk_strcat(lattr,iotk_strtrim(attr_tmp),ierr=ierrl)
  end if
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 186 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,lattr,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 199 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 201 "iotk_dat.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 205 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 209 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 214 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 227 "iotk_dat.spp"
      write(lunit,"(a)",iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 231 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 232 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 238 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 243 "iotk_dat.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 247 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
        end if
     end do
# 258 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 261 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 268 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_CHARACTER1_7


# 545 "iotk_dat.spp"

# 547 "iotk_dat.spp"
subroutine iotk_scan_dat_CHARACTER1_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER (kind=__IOTK_CHARACTER1,len=*)                        :: dat (:,:,:,:,:,:,:)
#else
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype), optional         :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),optional              :: attr
#else
  character(len=*),optional, intent(out) :: attr
#endif
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 572 "iotk_dat.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 576 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: lattr
  logical :: inside,foundl
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
# 592 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  call iotk_parse_dat(lattr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 598 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 598 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 602 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 607 "iotk_dat.spp"
  if(rlen ==-1) rlen  = len(dat)
# 609 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 615 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 617 "iotk_dat.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
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
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.12 ")
# 629 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 629 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_CHARACTER1_7


#endif
#endif

subroutine iotk_dat_dummy_CHARACTER1_7
  write(0,*)
end subroutine iotk_dat_dummy_CHARACTER1_7

