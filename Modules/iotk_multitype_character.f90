# 46 "iotk_multitype.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 56 "iotk_multitype.spp"

# 58 "iotk_multitype.spp"

#ifdef __IOTK_CHARACTER1
#if 0 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_CHARACTER1_0(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 77 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 82 "iotk_multitype.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 89 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 94 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 99 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"unit",unit)
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 108 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_attr(attr,"size",iotk_size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 113 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 117 "iotk_multitype.spp"
  call iotk_write_attr(attr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 119 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 131 "iotk_multitype.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 133 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if

  allocate(dattmp(iotk_size(dat)))
# 144 "iotk_multitype.spp"
     dattmp(1) = dat
# 156 "iotk_multitype.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 161 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 167 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  else
    if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 180 "iotk_multitype.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 184 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
        end if
     end do
# 195 "iotk_multitype.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 198 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 205 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_CHARACTER1_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 483 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 487 "iotk_multitype.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  call iotk_scan_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) goto 1
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 511 "iotk_multitype.spp"
  if(rlen ==-1) rlen  = len(dat)
# 513 "iotk_multitype.spp"

  allocate(tmpdat(iotk_size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 517 "iotk_multitype.spp"
        dat = tmpdat(1)
# 521 "iotk_multitype.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
    if(ierrl/=0) dat = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_CHARACTER1_0

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_CHARACTER1_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  CHARACTER(kind=__IOTK_CHARACTER1,len=*), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 685 "iotk_multitype.spp"
  logical :: lquot,lapos
# 687 "iotk_multitype.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 694 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 694 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 694 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 701 "iotk_multitype.spp"
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
# 725 "iotk_multitype.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 727 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 727 "iotk_multitype.spp"
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

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_CHARACTER1_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  CHARACTER(kind=__IOTK_CHARACTER1,len=*),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  CHARACTER(kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 757 "iotk_multitype.spp"
  character(iotk_vallenx) :: valctmp
  integer :: vallen,defaultlen
  logical :: leos
# 764 "iotk_multitype.spp"
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
# 774 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 781 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 787 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 792 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 801 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
  else
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 805 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    ierrl = - ierrl
    goto 1
  end if
# 810 "iotk_multitype.spp"
  call iotk_escape(valctmp,valc)
  vallen = iotk_strlen(valctmp)
  if(len(val) < vallen) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 813 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 855 "iotk_multitype.spp"
    if(leos) then
      if(ierrl/=0) then
        defaultlen = min(iotk_strlen(default),len(val))
        val(1:defaultlen) = default(1:defaultlen)
        if(defaultlen<len(val)) val(defaultlen+1:defaultlen+1)=iotk_eos
      end if
    else
      if(ierrl/=0) val = default
    end if
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_CHARACTER1_0
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_CHARACTER1_0
  write(0,*)
end subroutine iotk_dummy_CHARACTER1_0

# 44 "iotk_multitype.spp"

# 46 "iotk_multitype.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 56 "iotk_multitype.spp"

# 58 "iotk_multitype.spp"

#ifdef __IOTK_CHARACTER1
#if 1 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_CHARACTER1_1(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 77 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 82 "iotk_multitype.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 89 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 94 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 99 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"unit",unit)
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 108 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_attr(attr,"size",iotk_size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 113 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 117 "iotk_multitype.spp"
  call iotk_write_attr(attr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 119 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 131 "iotk_multitype.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 133 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if

  allocate(dattmp(iotk_size(dat)))
# 146 "iotk_multitype.spp"
#ifdef __IOTK_WORKAROUND3
# 148 "iotk_multitype.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 152 "iotk_multitype.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 156 "iotk_multitype.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 161 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 167 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  else
    if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 180 "iotk_multitype.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 184 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
        end if
     end do
# 195 "iotk_multitype.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 198 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 205 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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

# 218 "iotk_multitype.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_CHARACTER1(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 224 "iotk_multitype.spp"
    CHARACTER (kind=__IOTK_CHARACTER1,len=l), intent(out) :: out(n)
    CHARACTER (kind=__IOTK_CHARACTER1,len=l), intent(in)  :: in(n)
# 230 "iotk_multitype.spp"
    out = in
end subroutine iotk_private_pack_CHARACTER1

# 234 "iotk_multitype.spp"
recursive subroutine iotk_scan_dat_aux_CHARACTER1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,         intent(in)  :: unit
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(out) :: dat (:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
  integer(iotk_header_kind) :: idummy
  logical :: raw,binary
  integer :: lunit
  integer :: index,length,nexttag,iostat,altlength
  character(len=iotk_linlenx) :: line,altline
# 250 "iotk_multitype.spp"
  CHARACTER (kind=kind(dat), len=rlen) :: dattmp(ubound(dat,1))
# 260 "iotk_multitype.spp"
  lunit = iotk_phys_unit(unit)
  ierr = 0
  iostat = 0
  idummy = 0
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(unit=lunit,binary=binary,ierr=ierr)
  if(ierr/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 267 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
# 271 "iotk_multitype.spp"
  if(binary) then
    if(raw) then
      read(lunit,iostat=iostat) dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 275 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 275 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 275 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else
      read(lunit,iostat=iostat) idummy,dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 281 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 281 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 281 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  else
    if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,  iostat=iostat) dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 289 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 289 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 289 "iotk_multitype.spp"
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
# 298 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
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
# 309 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 309 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 309 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 314 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 319 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 319 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 319 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 324 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 324 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 324 "iotk_multitype.spp"
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
# 333 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 333 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Missing dat')
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dattmp
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 340 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 340 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 340 "iotk_multitype.spp"
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
# 464 "iotk_multitype.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 465 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
end subroutine iotk_scan_dat_aux_CHARACTER1
# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_CHARACTER1_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 483 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 487 "iotk_multitype.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  call iotk_scan_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) goto 1
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 511 "iotk_multitype.spp"
  if(rlen ==-1) rlen  = len(dat)
# 513 "iotk_multitype.spp"

  allocate(tmpdat(iotk_size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 519 "iotk_multitype.spp"
        dat = reshape(tmpdat,shape(dat))
# 521 "iotk_multitype.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
    if(ierrl/=0) dat = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_CHARACTER1_1

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_CHARACTER1_1
  write(0,*)
end subroutine iotk_dummy_CHARACTER1_1

# 44 "iotk_multitype.spp"

# 46 "iotk_multitype.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 56 "iotk_multitype.spp"

# 58 "iotk_multitype.spp"

#ifdef __IOTK_CHARACTER1
#if 2 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_CHARACTER1_2(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 77 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 82 "iotk_multitype.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 89 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 94 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 99 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"unit",unit)
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 108 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_attr(attr,"size",iotk_size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 113 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 117 "iotk_multitype.spp"
  call iotk_write_attr(attr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 119 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 131 "iotk_multitype.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 133 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if

  allocate(dattmp(iotk_size(dat)))
# 146 "iotk_multitype.spp"
#ifdef __IOTK_WORKAROUND3
# 148 "iotk_multitype.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 152 "iotk_multitype.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 156 "iotk_multitype.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 161 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 167 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  else
    if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 180 "iotk_multitype.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 184 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
        end if
     end do
# 195 "iotk_multitype.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 198 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 205 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_CHARACTER1_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 483 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 487 "iotk_multitype.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  call iotk_scan_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) goto 1
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 511 "iotk_multitype.spp"
  if(rlen ==-1) rlen  = len(dat)
# 513 "iotk_multitype.spp"

  allocate(tmpdat(iotk_size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 519 "iotk_multitype.spp"
        dat = reshape(tmpdat,shape(dat))
# 521 "iotk_multitype.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
    if(ierrl/=0) dat = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_CHARACTER1_2

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_CHARACTER1_2
  write(0,*)
end subroutine iotk_dummy_CHARACTER1_2

# 44 "iotk_multitype.spp"

# 46 "iotk_multitype.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 56 "iotk_multitype.spp"

# 58 "iotk_multitype.spp"

#ifdef __IOTK_CHARACTER1
#if 3 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_CHARACTER1_3(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 77 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 82 "iotk_multitype.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 89 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 94 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 99 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"unit",unit)
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 108 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_attr(attr,"size",iotk_size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 113 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 117 "iotk_multitype.spp"
  call iotk_write_attr(attr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 119 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 131 "iotk_multitype.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 133 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if

  allocate(dattmp(iotk_size(dat)))
# 146 "iotk_multitype.spp"
#ifdef __IOTK_WORKAROUND3
# 148 "iotk_multitype.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 152 "iotk_multitype.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 156 "iotk_multitype.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 161 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 167 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  else
    if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 180 "iotk_multitype.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 184 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
        end if
     end do
# 195 "iotk_multitype.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 198 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 205 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_CHARACTER1_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 483 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 487 "iotk_multitype.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  call iotk_scan_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) goto 1
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 511 "iotk_multitype.spp"
  if(rlen ==-1) rlen  = len(dat)
# 513 "iotk_multitype.spp"

  allocate(tmpdat(iotk_size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 519 "iotk_multitype.spp"
        dat = reshape(tmpdat,shape(dat))
# 521 "iotk_multitype.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
    if(ierrl/=0) dat = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_CHARACTER1_3

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_CHARACTER1_3
  write(0,*)
end subroutine iotk_dummy_CHARACTER1_3

# 44 "iotk_multitype.spp"

# 46 "iotk_multitype.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 56 "iotk_multitype.spp"

# 58 "iotk_multitype.spp"

#ifdef __IOTK_CHARACTER1
#if 4 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_CHARACTER1_4(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 77 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 82 "iotk_multitype.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 89 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 94 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 99 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"unit",unit)
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 108 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_attr(attr,"size",iotk_size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 113 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 117 "iotk_multitype.spp"
  call iotk_write_attr(attr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 119 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 131 "iotk_multitype.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 133 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if

  allocate(dattmp(iotk_size(dat)))
# 146 "iotk_multitype.spp"
#ifdef __IOTK_WORKAROUND3
# 148 "iotk_multitype.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 152 "iotk_multitype.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 156 "iotk_multitype.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 161 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 167 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  else
    if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 180 "iotk_multitype.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 184 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
        end if
     end do
# 195 "iotk_multitype.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 198 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 205 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_CHARACTER1_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 483 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 487 "iotk_multitype.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  call iotk_scan_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) goto 1
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 511 "iotk_multitype.spp"
  if(rlen ==-1) rlen  = len(dat)
# 513 "iotk_multitype.spp"

  allocate(tmpdat(iotk_size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 519 "iotk_multitype.spp"
        dat = reshape(tmpdat,shape(dat))
# 521 "iotk_multitype.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
    if(ierrl/=0) dat = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_CHARACTER1_4

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_CHARACTER1_4
  write(0,*)
end subroutine iotk_dummy_CHARACTER1_4

# 44 "iotk_multitype.spp"

# 46 "iotk_multitype.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 56 "iotk_multitype.spp"

# 58 "iotk_multitype.spp"

#ifdef __IOTK_CHARACTER1
#if 5 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_CHARACTER1_5(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 77 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 82 "iotk_multitype.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 89 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 94 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 99 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"unit",unit)
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 108 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_attr(attr,"size",iotk_size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 113 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 117 "iotk_multitype.spp"
  call iotk_write_attr(attr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 119 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 131 "iotk_multitype.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 133 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if

  allocate(dattmp(iotk_size(dat)))
# 146 "iotk_multitype.spp"
#ifdef __IOTK_WORKAROUND3
# 148 "iotk_multitype.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 152 "iotk_multitype.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 156 "iotk_multitype.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 161 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 167 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  else
    if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 180 "iotk_multitype.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 184 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
        end if
     end do
# 195 "iotk_multitype.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 198 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 205 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_CHARACTER1_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 483 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 487 "iotk_multitype.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  call iotk_scan_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) goto 1
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 511 "iotk_multitype.spp"
  if(rlen ==-1) rlen  = len(dat)
# 513 "iotk_multitype.spp"

  allocate(tmpdat(iotk_size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 519 "iotk_multitype.spp"
        dat = reshape(tmpdat,shape(dat))
# 521 "iotk_multitype.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
    if(ierrl/=0) dat = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_CHARACTER1_5

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_CHARACTER1_5
  write(0,*)
end subroutine iotk_dummy_CHARACTER1_5

# 44 "iotk_multitype.spp"

# 46 "iotk_multitype.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 56 "iotk_multitype.spp"

# 58 "iotk_multitype.spp"

#ifdef __IOTK_CHARACTER1
#if 6 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_CHARACTER1_6(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 77 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 82 "iotk_multitype.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 89 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 94 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 99 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"unit",unit)
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 108 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_attr(attr,"size",iotk_size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 113 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 117 "iotk_multitype.spp"
  call iotk_write_attr(attr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 119 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 131 "iotk_multitype.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 133 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if

  allocate(dattmp(iotk_size(dat)))
# 146 "iotk_multitype.spp"
#ifdef __IOTK_WORKAROUND3
# 148 "iotk_multitype.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 152 "iotk_multitype.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 156 "iotk_multitype.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 161 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 167 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  else
    if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 180 "iotk_multitype.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 184 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
        end if
     end do
# 195 "iotk_multitype.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 198 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 205 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_CHARACTER1_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 483 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 487 "iotk_multitype.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  call iotk_scan_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) goto 1
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 511 "iotk_multitype.spp"
  if(rlen ==-1) rlen  = len(dat)
# 513 "iotk_multitype.spp"

  allocate(tmpdat(iotk_size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 519 "iotk_multitype.spp"
        dat = reshape(tmpdat,shape(dat))
# 521 "iotk_multitype.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
    if(ierrl/=0) dat = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_CHARACTER1_6

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_CHARACTER1_6
  write(0,*)
end subroutine iotk_dummy_CHARACTER1_6

# 44 "iotk_multitype.spp"

# 46 "iotk_multitype.spp"

!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 56 "iotk_multitype.spp"

# 58 "iotk_multitype.spp"

#ifdef __IOTK_CHARACTER1
#if 7 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_CHARACTER1_7(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 77 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)),allocatable :: dattmp(:)
  character(len=iotk_linlenx) :: linetmp
# 82 "iotk_multitype.spp"
  integer :: itmp
  ierrl = 0
  iostat = 0
  lunit = iotk_phys_unit(unit)
  call iotk_unit_get(lunit,raw=raw)
  call iotk_inquire(lunit,binary=binary,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 89 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 94 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 99 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 103 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"unit",unit)
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 103 "iotk_multitype.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("CHARACTER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 108 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_attr(attr,"size",iotk_size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 113 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 117 "iotk_multitype.spp"
  call iotk_write_attr(attr,"len",len(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 119 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 131 "iotk_multitype.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 133 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if

  allocate(dattmp(iotk_size(dat)))
# 146 "iotk_multitype.spp"
#ifdef __IOTK_WORKAROUND3
# 148 "iotk_multitype.spp"
     call iotk_private_pack_CHARACTER1(dattmp,dat,size(dattmp),len(dattmp))
# 152 "iotk_multitype.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 156 "iotk_multitype.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 161 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 167 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  else
    if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 175 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 180 "iotk_multitype.spp"
     do itmp = 1 , size(dattmp)
       call iotk_deescape(linetmp,dattmp(itmp))
       write(lunit,"(a)",iostat=iostat) linetmp(1:iotk_strlen(linetmp))
       if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 184 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
        end if
     end do
# 195 "iotk_multitype.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 198 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 205 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
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

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_CHARACTER1_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 483 "iotk_multitype.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=len(dat)), allocatable :: tmpdat(:)
# 487 "iotk_multitype.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  call iotk_scan_begin(unit,name,attr,ierr=ierrl)
  if(ierrl/=0) goto 1
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"CHARACTER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","CHARACTER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 511 "iotk_multitype.spp"
  if(rlen ==-1) rlen  = len(dat)
# 513 "iotk_multitype.spp"

  allocate(tmpdat(iotk_size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 519 "iotk_multitype.spp"
        dat = reshape(tmpdat,shape(dat))
# 521 "iotk_multitype.spp"
  deallocate(tmpdat)
1 continue
  if(inside) then
    call iotk_scan_end(unit,name,ierr=ierrl2)
    if(ierrl2/=0) then
      call iotk_error_clear(ierrl)
      ierrl=ierrl2
    end if
  end if
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
    if(ierrl/=0) dat = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. (.not.present(found) .and. .not.present(default))) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_dat_CHARACTER1_7

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_CHARACTER1_7
  write(0,*)
end subroutine iotk_dummy_CHARACTER1_7

