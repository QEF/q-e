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

#ifdef __IOTK_COMPLEX1
#if 0 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX1_0(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX1_0

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX1_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX1_0

# 547 "iotk_multitype.spp"
subroutine iotk_write_COMPLEX1(val,string,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in) :: val(:)
  character(len=*), intent(out) :: string
  integer, intent(out) :: ierr
  character(len=100) :: tmpval
  integer :: index,iostat
  ierr = 0
  iostat = 0 
  string(1:1) = iotk_eos
  if(size(val)==0) return
  if(len(string)==0) then
    call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 561 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
  do index=1,size(val)
# 578 "iotk_multitype.spp"
    write(tmpval,trim(iotk_wfmt("COMPLEX",kind(val),size(val),-1)),iostat=iostat) val(index)
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 580 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 580 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 580 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    call iotk_strcat(string,trim(adjustl(tmpval))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 585 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
      return
    end if
# 589 "iotk_multitype.spp"
  end do
! taglio l'ultimo spazio
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_COMPLEX1
# 595 "iotk_multitype.spp"

# 599 "iotk_multitype.spp"
subroutine iotk_read_COMPLEX1(val,string,index,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX1), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  logical :: check
  integer :: pos,pos1,iostat
  integer :: maxindex
# 611 "iotk_multitype.spp"
  real(kind=__IOTK_COMPLEX1) :: tmpreal
  complex(kind=__IOTK_COMPLEX1) :: tmpcomplex
# 614 "iotk_multitype.spp"
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 619 "iotk_multitype.spp"
   maxindex = 2 * size(val)
# 623 "iotk_multitype.spp"
! PER ORA CONSIDERA LE VIRGOLE COME SPAZII
  do
    pos = verify(string(pos1+1:)," ,")+pos1
    if(pos==pos1) exit
    pos = pos - 1
    pos1 = scan(string(pos+1:)," ,")+pos
    if(pos1==pos) pos1 = len(string) + 1
!LEGGI string(pos+1:pos1-1)
    index = index+1
    if(index>maxindex) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 633 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 633 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 642 "iotk_multitype.spp"
    read(string(pos+1:pos1-1),"(G100.95)",iostat=iostat) tmpreal
    if(modulo(index,2)==1) then
      tmpcomplex = cmplx(tmpreal,aimag((val((index+1)/2))))
    else
      tmpcomplex = cmplx(real(val((index+1)/2)),tmpreal)
    end if
    val((index+1)/2) = tmpcomplex
# 656 "iotk_multitype.spp"
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 657 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 657 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Error reading from string')
# 657 "iotk_multitype.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
# 657 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
# 661 "iotk_multitype.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_COMPLEX1
# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX1_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 716 "iotk_multitype.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX1_0

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX1_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX1), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 842 "iotk_multitype.spp"
  val = tmpval(1)
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_0
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX1_0
  write(0,*)
end subroutine iotk_dummy_COMPLEX1_0

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

#ifdef __IOTK_COMPLEX1
#if 1 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX1_1(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX1(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX1_1

# 218 "iotk_multitype.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_COMPLEX1(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 227 "iotk_multitype.spp"
    COMPLEX (kind=__IOTK_COMPLEX1), intent(out) :: out(n)
    COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: in(n)
# 230 "iotk_multitype.spp"
    out = in
end subroutine iotk_private_pack_COMPLEX1

# 234 "iotk_multitype.spp"
recursive subroutine iotk_scan_dat_aux_COMPLEX1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,         intent(in)  :: unit
  COMPLEX (kind=__IOTK_COMPLEX1), intent(out) :: dat (:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
  integer(iotk_header_kind) :: idummy
  logical :: raw,binary
  integer :: lunit
  integer :: index,length,nexttag,iostat,altlength
  character(len=iotk_linlenx) :: line,altline
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX2
  COMPLEX (__IOTK_COMPLEX2), allocatable :: dat2 (:)
#endif
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX3
  COMPLEX (__IOTK_COMPLEX3), allocatable :: dat3 (:)
#endif
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX4
  COMPLEX (__IOTK_COMPLEX4), allocatable :: dat4 (:)
#endif
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
# 351 "iotk_multitype.spp"
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
        read(lunit,iostat=iostat) dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 357 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 357 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 357 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 363 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 363 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 363 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX2
    case(kind(dat2))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat2(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat2
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat2
# 391 "iotk_multitype.spp"
      deallocate(dat2)
#endif
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX3
    case(kind(dat3))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat3(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat3
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat3
# 391 "iotk_multitype.spp"
      deallocate(dat3)
#endif
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX4
    case(kind(dat4))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat4(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat4
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat4
# 391 "iotk_multitype.spp"
      deallocate(dat4)
#endif
# 395 "iotk_multitype.spp"
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 396 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 396 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Kind incompatibility')
# 396 "iotk_multitype.spp"
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 402 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 402 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 402 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 410 "iotk_multitype.spp"
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
# 421 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 421 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 421 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 426 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 431 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 431 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 431 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 436 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 436 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 436 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 442 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
          return
        end if
# 446 "iotk_multitype.spp"
        if(index == 2 * size(dat)) exit
# 450 "iotk_multitype.spp"
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 451 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 458 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 458 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 458 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
# 464 "iotk_multitype.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 465 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
end subroutine iotk_scan_dat_aux_COMPLEX1
# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX1_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX1_1

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX1_1(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX1_1

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX1_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX1), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_1
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX1_1
  write(0,*)
end subroutine iotk_dummy_COMPLEX1_1

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

#ifdef __IOTK_COMPLEX1
#if 2 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX1_2(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX1(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX1_2

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX1_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX1_2

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX1_2(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX1_2

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX1_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX1), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_2
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX1_2
  write(0,*)
end subroutine iotk_dummy_COMPLEX1_2

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

#ifdef __IOTK_COMPLEX1
#if 3 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX1_3(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX1(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX1_3

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX1_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX1_3

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX1_3(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX1_3

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX1_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX1), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_3
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX1_3
  write(0,*)
end subroutine iotk_dummy_COMPLEX1_3

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

#ifdef __IOTK_COMPLEX1
#if 4 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX1_4(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX1(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX1_4

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX1_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX1_4

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX1_4(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX1_4

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX1_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX1), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_4
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX1_4
  write(0,*)
end subroutine iotk_dummy_COMPLEX1_4

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

#ifdef __IOTK_COMPLEX1
#if 5 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX1_5(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX1(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX1_5

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX1_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX1_5

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX1_5(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX1_5

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX1_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX1), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_5
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX1_5
  write(0,*)
end subroutine iotk_dummy_COMPLEX1_5

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

#ifdef __IOTK_COMPLEX1
#if 6 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX1_6(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX1(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX1_6

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX1_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX1_6

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX1_6(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX1_6

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX1_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX1), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_6
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX1_6
  write(0,*)
end subroutine iotk_dummy_COMPLEX1_6

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

#ifdef __IOTK_COMPLEX1
#if 7 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX1_7(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX1(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX1_7

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX1_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX1_7

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX1_7(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX1_7

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX1_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX1), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX1_7
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX1_7
  write(0,*)
end subroutine iotk_dummy_COMPLEX1_7

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

#ifdef __IOTK_COMPLEX2
#if 0 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX2_0(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX2_0

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX2_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX2_0

# 547 "iotk_multitype.spp"
subroutine iotk_write_COMPLEX2(val,string,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in) :: val(:)
  character(len=*), intent(out) :: string
  integer, intent(out) :: ierr
  character(len=100) :: tmpval
  integer :: index,iostat
  ierr = 0
  iostat = 0 
  string(1:1) = iotk_eos
  if(size(val)==0) return
  if(len(string)==0) then
    call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 561 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
  do index=1,size(val)
# 578 "iotk_multitype.spp"
    write(tmpval,trim(iotk_wfmt("COMPLEX",kind(val),size(val),-1)),iostat=iostat) val(index)
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 580 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 580 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 580 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    call iotk_strcat(string,trim(adjustl(tmpval))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 585 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
      return
    end if
# 589 "iotk_multitype.spp"
  end do
! taglio l'ultimo spazio
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_COMPLEX2
# 595 "iotk_multitype.spp"

# 599 "iotk_multitype.spp"
subroutine iotk_read_COMPLEX2(val,string,index,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX2), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  logical :: check
  integer :: pos,pos1,iostat
  integer :: maxindex
# 611 "iotk_multitype.spp"
  real(kind=__IOTK_COMPLEX2) :: tmpreal
  complex(kind=__IOTK_COMPLEX2) :: tmpcomplex
# 614 "iotk_multitype.spp"
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 619 "iotk_multitype.spp"
   maxindex = 2 * size(val)
# 623 "iotk_multitype.spp"
! PER ORA CONSIDERA LE VIRGOLE COME SPAZII
  do
    pos = verify(string(pos1+1:)," ,")+pos1
    if(pos==pos1) exit
    pos = pos - 1
    pos1 = scan(string(pos+1:)," ,")+pos
    if(pos1==pos) pos1 = len(string) + 1
!LEGGI string(pos+1:pos1-1)
    index = index+1
    if(index>maxindex) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 633 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 633 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 642 "iotk_multitype.spp"
    read(string(pos+1:pos1-1),"(G100.95)",iostat=iostat) tmpreal
    if(modulo(index,2)==1) then
      tmpcomplex = cmplx(tmpreal,aimag((val((index+1)/2))))
    else
      tmpcomplex = cmplx(real(val((index+1)/2)),tmpreal)
    end if
    val((index+1)/2) = tmpcomplex
# 656 "iotk_multitype.spp"
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 657 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 657 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Error reading from string')
# 657 "iotk_multitype.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
# 657 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
# 661 "iotk_multitype.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_COMPLEX2
# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX2_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 716 "iotk_multitype.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX2_0

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX2_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX2), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 842 "iotk_multitype.spp"
  val = tmpval(1)
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_0
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX2_0
  write(0,*)
end subroutine iotk_dummy_COMPLEX2_0

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

#ifdef __IOTK_COMPLEX2
#if 1 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX2_1(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX2(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX2_1

# 218 "iotk_multitype.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_COMPLEX2(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 227 "iotk_multitype.spp"
    COMPLEX (kind=__IOTK_COMPLEX2), intent(out) :: out(n)
    COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: in(n)
# 230 "iotk_multitype.spp"
    out = in
end subroutine iotk_private_pack_COMPLEX2

# 234 "iotk_multitype.spp"
recursive subroutine iotk_scan_dat_aux_COMPLEX2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,         intent(in)  :: unit
  COMPLEX (kind=__IOTK_COMPLEX2), intent(out) :: dat (:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
  integer(iotk_header_kind) :: idummy
  logical :: raw,binary
  integer :: lunit
  integer :: index,length,nexttag,iostat,altlength
  character(len=iotk_linlenx) :: line,altline
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX1
  COMPLEX (__IOTK_COMPLEX1), allocatable :: dat1 (:)
#endif
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX3
  COMPLEX (__IOTK_COMPLEX3), allocatable :: dat3 (:)
#endif
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX4
  COMPLEX (__IOTK_COMPLEX4), allocatable :: dat4 (:)
#endif
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
# 351 "iotk_multitype.spp"
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
        read(lunit,iostat=iostat) dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 357 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 357 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 357 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 363 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 363 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 363 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX1
    case(kind(dat1))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat1(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat1
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat1
# 391 "iotk_multitype.spp"
      deallocate(dat1)
#endif
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX3
    case(kind(dat3))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat3(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat3
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat3
# 391 "iotk_multitype.spp"
      deallocate(dat3)
#endif
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX4
    case(kind(dat4))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat4(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat4
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat4
# 391 "iotk_multitype.spp"
      deallocate(dat4)
#endif
# 395 "iotk_multitype.spp"
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 396 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 396 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Kind incompatibility')
# 396 "iotk_multitype.spp"
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 402 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 402 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 402 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 410 "iotk_multitype.spp"
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
# 421 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 421 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 421 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 426 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 431 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 431 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 431 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 436 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 436 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 436 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 442 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
          return
        end if
# 446 "iotk_multitype.spp"
        if(index == 2 * size(dat)) exit
# 450 "iotk_multitype.spp"
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 451 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 458 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 458 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 458 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
# 464 "iotk_multitype.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 465 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
end subroutine iotk_scan_dat_aux_COMPLEX2
# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX2_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX2_1

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX2_1(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX2_1

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX2_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX2), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_1
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX2_1
  write(0,*)
end subroutine iotk_dummy_COMPLEX2_1

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

#ifdef __IOTK_COMPLEX2
#if 2 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX2_2(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX2(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX2_2

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX2_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX2_2

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX2_2(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX2_2

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX2_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX2), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_2
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX2_2
  write(0,*)
end subroutine iotk_dummy_COMPLEX2_2

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

#ifdef __IOTK_COMPLEX2
#if 3 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX2_3(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX2(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX2_3

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX2_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX2_3

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX2_3(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX2_3

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX2_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX2), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_3
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX2_3
  write(0,*)
end subroutine iotk_dummy_COMPLEX2_3

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

#ifdef __IOTK_COMPLEX2
#if 4 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX2_4(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX2(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX2_4

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX2_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX2_4

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX2_4(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX2_4

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX2_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX2), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_4
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX2_4
  write(0,*)
end subroutine iotk_dummy_COMPLEX2_4

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

#ifdef __IOTK_COMPLEX2
#if 5 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX2_5(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX2(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX2_5

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX2_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX2_5

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX2_5(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX2_5

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX2_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX2), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_5
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX2_5
  write(0,*)
end subroutine iotk_dummy_COMPLEX2_5

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

#ifdef __IOTK_COMPLEX2
#if 6 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX2_6(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX2(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX2_6

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX2_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX2_6

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX2_6(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX2_6

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX2_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX2), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_6
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX2_6
  write(0,*)
end subroutine iotk_dummy_COMPLEX2_6

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

#ifdef __IOTK_COMPLEX2
#if 7 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX2_7(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX2(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX2_7

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX2_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX2_7

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX2_7(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX2_7

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX2_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX2), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_7
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX2_7
  write(0,*)
end subroutine iotk_dummy_COMPLEX2_7

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

#ifdef __IOTK_COMPLEX3
#if 0 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX3_0(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX3_0

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX3_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX3_0

# 547 "iotk_multitype.spp"
subroutine iotk_write_COMPLEX3(val,string,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in) :: val(:)
  character(len=*), intent(out) :: string
  integer, intent(out) :: ierr
  character(len=100) :: tmpval
  integer :: index,iostat
  ierr = 0
  iostat = 0 
  string(1:1) = iotk_eos
  if(size(val)==0) return
  if(len(string)==0) then
    call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 561 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
  do index=1,size(val)
# 578 "iotk_multitype.spp"
    write(tmpval,trim(iotk_wfmt("COMPLEX",kind(val),size(val),-1)),iostat=iostat) val(index)
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 580 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 580 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 580 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    call iotk_strcat(string,trim(adjustl(tmpval))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 585 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
      return
    end if
# 589 "iotk_multitype.spp"
  end do
! taglio l'ultimo spazio
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_COMPLEX3
# 595 "iotk_multitype.spp"

# 599 "iotk_multitype.spp"
subroutine iotk_read_COMPLEX3(val,string,index,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX3), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  logical :: check
  integer :: pos,pos1,iostat
  integer :: maxindex
# 611 "iotk_multitype.spp"
  real(kind=__IOTK_COMPLEX3) :: tmpreal
  complex(kind=__IOTK_COMPLEX3) :: tmpcomplex
# 614 "iotk_multitype.spp"
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 619 "iotk_multitype.spp"
   maxindex = 2 * size(val)
# 623 "iotk_multitype.spp"
! PER ORA CONSIDERA LE VIRGOLE COME SPAZII
  do
    pos = verify(string(pos1+1:)," ,")+pos1
    if(pos==pos1) exit
    pos = pos - 1
    pos1 = scan(string(pos+1:)," ,")+pos
    if(pos1==pos) pos1 = len(string) + 1
!LEGGI string(pos+1:pos1-1)
    index = index+1
    if(index>maxindex) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 633 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 633 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 642 "iotk_multitype.spp"
    read(string(pos+1:pos1-1),"(G100.95)",iostat=iostat) tmpreal
    if(modulo(index,2)==1) then
      tmpcomplex = cmplx(tmpreal,aimag((val((index+1)/2))))
    else
      tmpcomplex = cmplx(real(val((index+1)/2)),tmpreal)
    end if
    val((index+1)/2) = tmpcomplex
# 656 "iotk_multitype.spp"
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 657 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 657 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Error reading from string')
# 657 "iotk_multitype.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
# 657 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
# 661 "iotk_multitype.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_COMPLEX3
# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX3_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 716 "iotk_multitype.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX3_0

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX3_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX3), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 842 "iotk_multitype.spp"
  val = tmpval(1)
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX3_0
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX3_0
  write(0,*)
end subroutine iotk_dummy_COMPLEX3_0

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

#ifdef __IOTK_COMPLEX3
#if 1 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX3_1(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX3(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX3_1

# 218 "iotk_multitype.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_COMPLEX3(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 227 "iotk_multitype.spp"
    COMPLEX (kind=__IOTK_COMPLEX3), intent(out) :: out(n)
    COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: in(n)
# 230 "iotk_multitype.spp"
    out = in
end subroutine iotk_private_pack_COMPLEX3

# 234 "iotk_multitype.spp"
recursive subroutine iotk_scan_dat_aux_COMPLEX3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,         intent(in)  :: unit
  COMPLEX (kind=__IOTK_COMPLEX3), intent(out) :: dat (:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
  integer(iotk_header_kind) :: idummy
  logical :: raw,binary
  integer :: lunit
  integer :: index,length,nexttag,iostat,altlength
  character(len=iotk_linlenx) :: line,altline
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX1
  COMPLEX (__IOTK_COMPLEX1), allocatable :: dat1 (:)
#endif
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX2
  COMPLEX (__IOTK_COMPLEX2), allocatable :: dat2 (:)
#endif
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX4
  COMPLEX (__IOTK_COMPLEX4), allocatable :: dat4 (:)
#endif
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
# 351 "iotk_multitype.spp"
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
        read(lunit,iostat=iostat) dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 357 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 357 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 357 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 363 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 363 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 363 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX1
    case(kind(dat1))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat1(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat1
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat1
# 391 "iotk_multitype.spp"
      deallocate(dat1)
#endif
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX2
    case(kind(dat2))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat2(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat2
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat2
# 391 "iotk_multitype.spp"
      deallocate(dat2)
#endif
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX4
    case(kind(dat4))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat4(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat4
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat4
# 391 "iotk_multitype.spp"
      deallocate(dat4)
#endif
# 395 "iotk_multitype.spp"
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 396 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 396 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Kind incompatibility')
# 396 "iotk_multitype.spp"
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 402 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 402 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 402 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 410 "iotk_multitype.spp"
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
# 421 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 421 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 421 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 426 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 431 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 431 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 431 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 436 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 436 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 436 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 442 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
          return
        end if
# 446 "iotk_multitype.spp"
        if(index == 2 * size(dat)) exit
# 450 "iotk_multitype.spp"
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 451 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 458 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 458 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 458 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
# 464 "iotk_multitype.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 465 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
end subroutine iotk_scan_dat_aux_COMPLEX3
# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX3_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX3_1

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX3_1(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX3_1

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX3_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX3), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX3_1
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX3_1
  write(0,*)
end subroutine iotk_dummy_COMPLEX3_1

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

#ifdef __IOTK_COMPLEX3
#if 2 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX3_2(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX3(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX3_2

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX3_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX3_2

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX3_2(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX3_2

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX3_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX3), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX3_2
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX3_2
  write(0,*)
end subroutine iotk_dummy_COMPLEX3_2

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

#ifdef __IOTK_COMPLEX3
#if 3 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX3_3(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX3(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX3_3

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX3_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX3_3

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX3_3(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX3_3

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX3_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX3), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX3_3
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX3_3
  write(0,*)
end subroutine iotk_dummy_COMPLEX3_3

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

#ifdef __IOTK_COMPLEX3
#if 4 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX3_4(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX3(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX3_4

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX3_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX3_4

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX3_4(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX3_4

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX3_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX3), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX3_4
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX3_4
  write(0,*)
end subroutine iotk_dummy_COMPLEX3_4

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

#ifdef __IOTK_COMPLEX3
#if 5 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX3_5(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX3(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX3_5

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX3_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX3_5

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX3_5(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX3_5

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX3_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX3), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX3_5
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX3_5
  write(0,*)
end subroutine iotk_dummy_COMPLEX3_5

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

#ifdef __IOTK_COMPLEX3
#if 6 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX3_6(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX3(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX3_6

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX3_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX3_6

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX3_6(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX3_6

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX3_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX3), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX3_6
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX3_6
  write(0,*)
end subroutine iotk_dummy_COMPLEX3_6

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

#ifdef __IOTK_COMPLEX3
#if 7 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX3_7(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX3(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX3_7

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX3_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX3_7

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX3_7(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX3_7

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX3_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX3), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX3_7
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX3_7
  write(0,*)
end subroutine iotk_dummy_COMPLEX3_7

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

#ifdef __IOTK_COMPLEX4
#if 0 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX4_0(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX4_0

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX4_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX4_0

# 547 "iotk_multitype.spp"
subroutine iotk_write_COMPLEX4(val,string,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in) :: val(:)
  character(len=*), intent(out) :: string
  integer, intent(out) :: ierr
  character(len=100) :: tmpval
  integer :: index,iostat
  ierr = 0
  iostat = 0 
  string(1:1) = iotk_eos
  if(size(val)==0) return
  if(len(string)==0) then
    call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 561 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
  do index=1,size(val)
# 578 "iotk_multitype.spp"
    write(tmpval,trim(iotk_wfmt("COMPLEX",kind(val),size(val),-1)),iostat=iostat) val(index)
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 580 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 580 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 580 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    call iotk_strcat(string,trim(adjustl(tmpval))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 585 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
      return
    end if
# 589 "iotk_multitype.spp"
  end do
! taglio l'ultimo spazio
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_COMPLEX4
# 595 "iotk_multitype.spp"

# 599 "iotk_multitype.spp"
subroutine iotk_read_COMPLEX4(val,string,index,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX4), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  logical :: check
  integer :: pos,pos1,iostat
  integer :: maxindex
# 611 "iotk_multitype.spp"
  real(kind=__IOTK_COMPLEX4) :: tmpreal
  complex(kind=__IOTK_COMPLEX4) :: tmpcomplex
# 614 "iotk_multitype.spp"
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 619 "iotk_multitype.spp"
   maxindex = 2 * size(val)
# 623 "iotk_multitype.spp"
! PER ORA CONSIDERA LE VIRGOLE COME SPAZII
  do
    pos = verify(string(pos1+1:)," ,")+pos1
    if(pos==pos1) exit
    pos = pos - 1
    pos1 = scan(string(pos+1:)," ,")+pos
    if(pos1==pos) pos1 = len(string) + 1
!LEGGI string(pos+1:pos1-1)
    index = index+1
    if(index>maxindex) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 633 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 633 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 642 "iotk_multitype.spp"
    read(string(pos+1:pos1-1),"(G100.95)",iostat=iostat) tmpreal
    if(modulo(index,2)==1) then
      tmpcomplex = cmplx(tmpreal,aimag((val((index+1)/2))))
    else
      tmpcomplex = cmplx(real(val((index+1)/2)),tmpreal)
    end if
    val((index+1)/2) = tmpcomplex
# 656 "iotk_multitype.spp"
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 657 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 657 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Error reading from string')
# 657 "iotk_multitype.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
# 657 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
# 661 "iotk_multitype.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_COMPLEX4
# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX4_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 716 "iotk_multitype.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX4_0

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX4_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX4), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 842 "iotk_multitype.spp"
  val = tmpval(1)
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_0
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX4_0
  write(0,*)
end subroutine iotk_dummy_COMPLEX4_0

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

#ifdef __IOTK_COMPLEX4
#if 1 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX4_1(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX4(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX4_1

# 218 "iotk_multitype.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_COMPLEX4(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 227 "iotk_multitype.spp"
    COMPLEX (kind=__IOTK_COMPLEX4), intent(out) :: out(n)
    COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: in(n)
# 230 "iotk_multitype.spp"
    out = in
end subroutine iotk_private_pack_COMPLEX4

# 234 "iotk_multitype.spp"
recursive subroutine iotk_scan_dat_aux_COMPLEX4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,         intent(in)  :: unit
  COMPLEX (kind=__IOTK_COMPLEX4), intent(out) :: dat (:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
  integer(iotk_header_kind) :: idummy
  logical :: raw,binary
  integer :: lunit
  integer :: index,length,nexttag,iostat,altlength
  character(len=iotk_linlenx) :: line,altline
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX1
  COMPLEX (__IOTK_COMPLEX1), allocatable :: dat1 (:)
#endif
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX2
  COMPLEX (__IOTK_COMPLEX2), allocatable :: dat2 (:)
#endif
# 254 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX3
  COMPLEX (__IOTK_COMPLEX3), allocatable :: dat3 (:)
#endif
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
# 351 "iotk_multitype.spp"
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
        read(lunit,iostat=iostat) dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 357 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 357 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 357 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 363 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 363 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 363 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX1
    case(kind(dat1))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat1(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat1
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat1
# 391 "iotk_multitype.spp"
      deallocate(dat1)
#endif
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX2
    case(kind(dat2))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat2(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat2
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat2
# 391 "iotk_multitype.spp"
      deallocate(dat2)
#endif
# 369 "iotk_multitype.spp"
#ifdef __IOTK_COMPLEX3
    case(kind(dat3))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 373 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
        return
      end if
      allocate(dat3(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat3
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 379 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 379 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 389 "iotk_multitype.spp"
      dat = dat3
# 391 "iotk_multitype.spp"
      deallocate(dat3)
#endif
# 395 "iotk_multitype.spp"
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 396 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 396 "iotk_multitype.spp"
call iotk_error_msg(ierr,'Kind incompatibility')
# 396 "iotk_multitype.spp"
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 402 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 402 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 402 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 410 "iotk_multitype.spp"
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
# 421 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 421 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 421 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 426 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 431 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 431 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 431 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 436 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 436 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 436 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 442 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
          return
        end if
# 446 "iotk_multitype.spp"
        if(index == 2 * size(dat)) exit
# 450 "iotk_multitype.spp"
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 451 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 458 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
# 458 "iotk_multitype.spp"
call iotk_error_msg(ierr,' ')
# 458 "iotk_multitype.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
# 464 "iotk_multitype.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 465 "iotk_multitype.spp"
call iotk_error_msg(ierr,"CVS $Revision: 1.1 $")
    return
  end if
end subroutine iotk_scan_dat_aux_COMPLEX4
# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX4_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX4_1

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX4_1(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX4_1

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX4_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX4), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_1
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX4_1
  write(0,*)
end subroutine iotk_dummy_COMPLEX4_1

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

#ifdef __IOTK_COMPLEX4
#if 2 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX4_2(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX4(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX4_2

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX4_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX4_2

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX4_2(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX4_2

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX4_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX4), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_2
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX4_2
  write(0,*)
end subroutine iotk_dummy_COMPLEX4_2

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

#ifdef __IOTK_COMPLEX4
#if 3 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX4_3(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX4(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX4_3

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX4_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX4_3

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX4_3(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX4_3

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX4_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX4), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_3
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX4_3
  write(0,*)
end subroutine iotk_dummy_COMPLEX4_3

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

#ifdef __IOTK_COMPLEX4
#if 4 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX4_4(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX4(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX4_4

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX4_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX4_4

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX4_4(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX4_4

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX4_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX4), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_4
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX4_4
  write(0,*)
end subroutine iotk_dummy_COMPLEX4_4

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

#ifdef __IOTK_COMPLEX4
#if 5 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX4_5(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX4(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX4_5

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX4_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX4_5

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX4_5(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX4_5

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX4_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX4), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_5
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX4_5
  write(0,*)
end subroutine iotk_dummy_COMPLEX4_5

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

#ifdef __IOTK_COMPLEX4
#if 6 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX4_6(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX4(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX4_6

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX4_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX4_6

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX4_6(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX4_6

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX4_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX4), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_6
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX4_6
  write(0,*)
end subroutine iotk_dummy_COMPLEX4_6

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

#ifdef __IOTK_COMPLEX4
#if 7 <= __IOTK_MAXRANK
# 62 "iotk_multitype.spp"
subroutine iotk_write_dat_COMPLEX4_7(unit,name,dat,fmt,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
# 80 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),allocatable :: dattmp(:)
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
  call iotk_write_attr(attr,"type",iotk_tolower("COMPLEX"),first=.true.,ierr=ierrl)
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
# 123 "iotk_multitype.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 126 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
    end if
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
# 150 "iotk_multitype.spp"
     call iotk_private_pack_COMPLEX4(dattmp,dat,size(dattmp),1)
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
# 189 "iotk_multitype.spp"
     write(lunit,fmt=iotk_wfmt("COMPLEX",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
      goto 1
     end if
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
end subroutine iotk_write_dat_COMPLEX4_7

# 470 "iotk_multitype.spp"

# 472 "iotk_multitype.spp"
subroutine iotk_scan_dat_COMPLEX4_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 485 "iotk_multitype.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),              allocatable :: tmpdat(:)
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
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"COMPLEX") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 502 "iotk_multitype.spp"
call iotk_error_msg(ierrl,' ')
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 502 "iotk_multitype.spp"
call iotk_error_write(ierrl,"type","COMPLEX")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==iotk_size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 506 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
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
end subroutine iotk_scan_dat_COMPLEX4_7

# 595 "iotk_multitype.spp"

# 666 "iotk_multitype.spp"

# 669 "iotk_multitype.spp"
subroutine iotk_write_attr_COMPLEX4_7(attr,name,val,first,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
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
# 714 "iotk_multitype.spp"
  delim = '"'
# 718 "iotk_multitype.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 720 "iotk_multitype.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 721 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
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
end subroutine iotk_write_attr_COMPLEX4_7

# 740 "iotk_multitype.spp"
subroutine iotk_scan_attr_COMPLEX4_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_interface
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 761 "iotk_multitype.spp"
  integer :: index
  COMPLEX(kind=__IOTK_COMPLEX4), allocatable :: tmpval (:)
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
# 826 "iotk_multitype.spp"
  allocate(tmpval(iotk_size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 830 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
    goto 1
  end if
# 834 "iotk_multitype.spp"
  if(index/=2*iotk_size(val)) then
# 838 "iotk_multitype.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,"CVS $Revision: 1.1 $")
# 838 "iotk_multitype.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 838 "iotk_multitype.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 844 "iotk_multitype.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 846 "iotk_multitype.spp"
  deallocate(tmpval)
# 848 "iotk_multitype.spp"
1 continue
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(default)) then
# 865 "iotk_multitype.spp"
    if(ierrl/=0) val = default
# 867 "iotk_multitype.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else if((present(found) .or. present(default)) .and. ierrl<0) then
    call iotk_error_clear(ierrl)
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_7
# 877 "iotk_multitype.spp"

#endif
#endif

subroutine iotk_dummy_COMPLEX4_7
  write(0,*)
end subroutine iotk_dummy_COMPLEX4_7

# 44 "iotk_multitype.spp"

