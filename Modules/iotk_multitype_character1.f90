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

