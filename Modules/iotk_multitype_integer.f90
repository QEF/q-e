# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER1
#if 0 <= __IOTK_MAXRANK

# 64 "iotk_attr.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_INTEGER1(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 73 "iotk_attr.spp"
    INTEGER (kind=__IOTK_INTEGER1), intent(out) :: out(n)
    INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: in(n)
# 76 "iotk_attr.spp"
    out = in
end subroutine iotk_private_pack_INTEGER1

# 81 "iotk_attr.spp"
subroutine iotk_write_INTEGER1(val,string,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  INTEGER(kind=__IOTK_INTEGER1), intent(in) :: val(:)
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
# 99 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
    return
  end if
  do index=1,size(val)
# 110 "iotk_attr.spp"
    call iotk_strcat(string,trim(iotk_itoa(val(index)))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 112 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
      return
    end if
# 127 "iotk_attr.spp"
  end do
! taglio l'ultimo spazio
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_INTEGER1
# 133 "iotk_attr.spp"

# 137 "iotk_attr.spp"
subroutine iotk_read_INTEGER1(val,string,index,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  INTEGER(kind=__IOTK_INTEGER1), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  logical :: check
  integer :: pos,pos1,iostat
  integer :: maxindex
# 154 "iotk_attr.spp"
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 161 "iotk_attr.spp"
    maxindex = size(val)
# 163 "iotk_attr.spp"
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
# 173 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
# 173 "iotk_attr.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 176 "iotk_attr.spp"
    call iotk_atoi(val(index),string(pos+1:pos1-1),check=check)
# 191 "iotk_attr.spp"
    if(.not.check) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 192 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
# 192 "iotk_attr.spp"
call iotk_error_msg(ierr,'Wrong string')
# 192 "iotk_attr.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
      return
    end if
# 201 "iotk_attr.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_INTEGER1
# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER1_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 259 "iotk_attr.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER1_0

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER1_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER1), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(1))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=1) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 386 "iotk_attr.spp"
  val = tmpval(1)
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER1_0
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER1_0
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER1_0

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER1
#if 1 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER1_1(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER1_1

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER1_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER1), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER1_1
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER1_1
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER1_1

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER1
#if 2 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER1_2(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER1_2

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER1_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER1), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER1_2
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER1_2
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER1_2

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER1
#if 3 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER1_3(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER1_3

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER1_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER1), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER1_3
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER1_3
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER1_3

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER1
#if 4 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER1_4(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER1_4

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER1_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER1), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER1_4
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER1_4
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER1_4

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER1
#if 5 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER1_5(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER1_5

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER1_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER1), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER1_5
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER1_5
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER1_5

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER1
#if 6 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER1_6(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER1_6

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER1_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER1), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER1_6
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER1_6
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER1_6

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER1
#if 7 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER1_7(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER1_7

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER1_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER1), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER1_7
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER1_7
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER1_7

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER2
#if 0 <= __IOTK_MAXRANK

# 64 "iotk_attr.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_INTEGER2(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 73 "iotk_attr.spp"
    INTEGER (kind=__IOTK_INTEGER2), intent(out) :: out(n)
    INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: in(n)
# 76 "iotk_attr.spp"
    out = in
end subroutine iotk_private_pack_INTEGER2

# 81 "iotk_attr.spp"
subroutine iotk_write_INTEGER2(val,string,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  INTEGER(kind=__IOTK_INTEGER2), intent(in) :: val(:)
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
# 99 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
    return
  end if
  do index=1,size(val)
# 110 "iotk_attr.spp"
    call iotk_strcat(string,trim(iotk_itoa(val(index)))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 112 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
      return
    end if
# 127 "iotk_attr.spp"
  end do
! taglio l'ultimo spazio
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_INTEGER2
# 133 "iotk_attr.spp"

# 137 "iotk_attr.spp"
subroutine iotk_read_INTEGER2(val,string,index,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  INTEGER(kind=__IOTK_INTEGER2), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  logical :: check
  integer :: pos,pos1,iostat
  integer :: maxindex
# 154 "iotk_attr.spp"
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 161 "iotk_attr.spp"
    maxindex = size(val)
# 163 "iotk_attr.spp"
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
# 173 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
# 173 "iotk_attr.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 176 "iotk_attr.spp"
    call iotk_atoi(val(index),string(pos+1:pos1-1),check=check)
# 191 "iotk_attr.spp"
    if(.not.check) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 192 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
# 192 "iotk_attr.spp"
call iotk_error_msg(ierr,'Wrong string')
# 192 "iotk_attr.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
      return
    end if
# 201 "iotk_attr.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_INTEGER2
# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER2_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 259 "iotk_attr.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER2_0

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER2_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER2), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(1))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=1) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 386 "iotk_attr.spp"
  val = tmpval(1)
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER2_0
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER2_0
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER2_0

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER2
#if 1 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER2_1(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER2_1

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER2_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER2), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER2_1
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER2_1
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER2_1

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER2
#if 2 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER2_2(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER2_2

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER2_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER2), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER2_2
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER2_2
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER2_2

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER2
#if 3 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER2_3(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER2_3

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER2_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER2), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER2_3
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER2_3
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER2_3

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER2
#if 4 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER2_4(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER2_4

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER2_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER2), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER2_4
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER2_4
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER2_4

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER2
#if 5 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER2_5(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER2_5

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER2_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER2), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER2_5
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER2_5
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER2_5

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER2
#if 6 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER2_6(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER2_6

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER2_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER2), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER2_6
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER2_6
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER2_6

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER2
#if 7 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER2_7(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER2_7

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER2_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER2), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER2_7
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER2_7
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER2_7

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER3
#if 0 <= __IOTK_MAXRANK

# 64 "iotk_attr.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_INTEGER3(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 73 "iotk_attr.spp"
    INTEGER (kind=__IOTK_INTEGER3), intent(out) :: out(n)
    INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: in(n)
# 76 "iotk_attr.spp"
    out = in
end subroutine iotk_private_pack_INTEGER3

# 81 "iotk_attr.spp"
subroutine iotk_write_INTEGER3(val,string,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  INTEGER(kind=__IOTK_INTEGER3), intent(in) :: val(:)
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
# 99 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
    return
  end if
  do index=1,size(val)
# 110 "iotk_attr.spp"
    call iotk_strcat(string,trim(iotk_itoa(val(index)))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 112 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
      return
    end if
# 127 "iotk_attr.spp"
  end do
! taglio l'ultimo spazio
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_INTEGER3
# 133 "iotk_attr.spp"

# 137 "iotk_attr.spp"
subroutine iotk_read_INTEGER3(val,string,index,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  INTEGER(kind=__IOTK_INTEGER3), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  logical :: check
  integer :: pos,pos1,iostat
  integer :: maxindex
# 154 "iotk_attr.spp"
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 161 "iotk_attr.spp"
    maxindex = size(val)
# 163 "iotk_attr.spp"
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
# 173 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
# 173 "iotk_attr.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 176 "iotk_attr.spp"
    call iotk_atoi(val(index),string(pos+1:pos1-1),check=check)
# 191 "iotk_attr.spp"
    if(.not.check) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 192 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
# 192 "iotk_attr.spp"
call iotk_error_msg(ierr,'Wrong string')
# 192 "iotk_attr.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
      return
    end if
# 201 "iotk_attr.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_INTEGER3
# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER3_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 259 "iotk_attr.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER3_0

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER3_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER3), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(1))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=1) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 386 "iotk_attr.spp"
  val = tmpval(1)
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER3_0
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER3_0
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER3_0

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER3
#if 1 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER3_1(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER3_1

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER3_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER3), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER3_1
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER3_1
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER3_1

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER3
#if 2 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER3_2(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER3_2

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER3_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER3), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER3_2
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER3_2
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER3_2

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER3
#if 3 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER3_3(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER3_3

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER3_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER3), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER3_3
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER3_3
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER3_3

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER3
#if 4 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER3_4(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER3_4

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER3_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER3), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER3_4
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER3_4
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER3_4

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER3
#if 5 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER3_5(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER3_5

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER3_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER3), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER3_5
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER3_5
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER3_5

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER3
#if 6 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER3_6(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER3_6

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER3_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER3), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER3_6
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER3_6
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER3_6

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER3
#if 7 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER3_7(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER3_7

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER3_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER3), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER3_7
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER3_7
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER3_7

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER4
#if 0 <= __IOTK_MAXRANK

# 64 "iotk_attr.spp"
! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_INTEGER4(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
# 73 "iotk_attr.spp"
    INTEGER (kind=__IOTK_INTEGER4), intent(out) :: out(n)
    INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: in(n)
# 76 "iotk_attr.spp"
    out = in
end subroutine iotk_private_pack_INTEGER4

# 81 "iotk_attr.spp"
subroutine iotk_write_INTEGER4(val,string,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  INTEGER(kind=__IOTK_INTEGER4), intent(in) :: val(:)
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
# 99 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
    return
  end if
  do index=1,size(val)
# 110 "iotk_attr.spp"
    call iotk_strcat(string,trim(iotk_itoa(val(index)))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
# 112 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
      return
    end if
# 127 "iotk_attr.spp"
  end do
! taglio l'ultimo spazio
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_INTEGER4
# 133 "iotk_attr.spp"

# 137 "iotk_attr.spp"
subroutine iotk_read_INTEGER4(val,string,index,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  INTEGER(kind=__IOTK_INTEGER4), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  logical :: check
  integer :: pos,pos1,iostat
  integer :: maxindex
# 154 "iotk_attr.spp"
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
# 161 "iotk_attr.spp"
    maxindex = size(val)
# 163 "iotk_attr.spp"
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
# 173 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
# 173 "iotk_attr.spp"
call iotk_error_msg(ierr,'Too many data')
    end if
# 176 "iotk_attr.spp"
    call iotk_atoi(val(index),string(pos+1:pos1-1),check=check)
# 191 "iotk_attr.spp"
    if(.not.check) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
# 192 "iotk_attr.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.1 ")
# 192 "iotk_attr.spp"
call iotk_error_msg(ierr,'Wrong string')
# 192 "iotk_attr.spp"
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
      return
    end if
# 201 "iotk_attr.spp"
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_INTEGER4
# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER4_0(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 259 "iotk_attr.spp"
  call iotk_write((/val/),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER4_0

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER4_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER4), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(1))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=1) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 386 "iotk_attr.spp"
  val = tmpval(1)
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER4_0
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER4_0
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER4_0

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER4
#if 1 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER4_1(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER4_1

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER4_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER4), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER4_1
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER4_1
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER4_1

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER4
#if 2 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER4_2(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER4_2

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER4_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER4), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER4_2
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER4_2
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER4_2

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER4
#if 3 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER4_3(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER4_3

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER4_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER4), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER4_3
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER4_3
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER4_3

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER4
#if 4 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER4_4(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER4_4

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER4_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER4), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER4_4
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER4_4
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER4_4

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER4
#if 5 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER4_5(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER4_5

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER4_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER4), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER4_5
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER4_5
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER4_5

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER4
#if 6 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER4_6(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER4_6

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER4_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER4), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER4_6
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER4_6
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER4_6

# 45 "iotk_attr.spp"

# 47 "iotk_attr.spp"

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

# 57 "iotk_attr.spp"

# 59 "iotk_attr.spp"

#ifdef __IOTK_INTEGER4
#if 7 <= __IOTK_MAXRANK

# 133 "iotk_attr.spp"

# 206 "iotk_attr.spp"

# 209 "iotk_attr.spp"
subroutine iotk_write_attr_INTEGER4_7(attr,name,val,first,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  integer :: iostat
  character :: delim
# 230 "iotk_attr.spp"
  character(iotk_vallenx) :: tmpval
  ierrl = 0
  iostat = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 237 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Wrong tag name')
# 237 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",trim(name))
    goto 1
  end if
  attlen = iotk_strlen(attr)
  if(attlen==len(attr)) attlen = len_trim(attr)
  namlen = len_trim(name)
# 257 "iotk_attr.spp"
  delim = '"'
# 261 "iotk_attr.spp"
  call iotk_write(pack(val,mask=.true.),tmpval,ierrl)
# 263 "iotk_attr.spp"
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 264 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 268 "iotk_attr.spp"
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
# 270 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 270 "iotk_attr.spp"
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
end subroutine iotk_write_attr_INTEGER4_7

# 283 "iotk_attr.spp"
subroutine iotk_scan_attr_INTEGER4_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
# 307 "iotk_attr.spp"
  integer :: index
  INTEGER(kind=__IOTK_INTEGER4), allocatable :: tmpval (:)
# 310 "iotk_attr.spp"
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
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 320 "iotk_attr.spp"
call iotk_error_msg(ierrl,'')
# 320 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==trim(name)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 327 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 333 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 338 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 347 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
      goto 1
    end if
  else
    goto 1
  end if
# 370 "iotk_attr.spp"
  allocate(tmpval(size(val)))
  index = 0
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 374 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
    goto 1
  end if
# 380 "iotk_attr.spp"
  if(index/=size(val)) then
# 382 "iotk_attr.spp"
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 382 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute size does not match')
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"attr",valc(1:iotk_strlen(valc)))
# 382 "iotk_attr.spp"
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
# 388 "iotk_attr.spp"
  val = reshape (source=tmpval,shape=shape(val))
# 390 "iotk_attr.spp"
  deallocate(tmpval)
# 392 "iotk_attr.spp"
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.1 ")
# 396 "iotk_attr.spp"
call iotk_error_msg(ierrl,'Attribute not found')
# 396 "iotk_attr.spp"
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
# 409 "iotk_attr.spp"
    val = default
# 411 "iotk_attr.spp"
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_INTEGER4_7
# 419 "iotk_attr.spp"

#endif
#endif

subroutine iotk_attr_dummy_INTEGER4_7
  write(0,*)
end subroutine iotk_attr_dummy_INTEGER4_7

# 45 "iotk_attr.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 0 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_0(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",1,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(1))
# 156 "iotk_dat.spp"
     dattmp(1) = dat
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_0


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==1) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(1))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 553 "iotk_dat.spp"
        dat = tmpdat(1)
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER1_0


#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_0
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_0

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 1 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_1(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_1


# 242 "iotk_dat.spp"
recursive subroutine iotk_scan_dat_aux_INTEGER1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only: iotk_read
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,         intent(in)  :: unit
  INTEGER (kind=__IOTK_INTEGER1), intent(out) :: dat (:)
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
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER2
  INTEGER (__IOTK_INTEGER2), allocatable :: dat2 (:)
#endif
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER3
  INTEGER (__IOTK_INTEGER3), allocatable :: dat3 (:)
#endif
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER4
  INTEGER (__IOTK_INTEGER4), allocatable :: dat4 (:)
#endif
# 274 "iotk_dat.spp"
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
# 285 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
    return
  end if
# 375 "iotk_dat.spp"
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
        read(lunit,iostat=iostat) dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 381 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 381 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 381 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 387 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 387 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 387 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER2
    case(kind(dat2))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat2(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat2
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat2
# 415 "iotk_dat.spp"
      deallocate(dat2)
#endif
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER3
    case(kind(dat3))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat3(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat3
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat3
# 415 "iotk_dat.spp"
      deallocate(dat3)
#endif
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER4
    case(kind(dat4))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat4(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat4
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat4
# 415 "iotk_dat.spp"
      deallocate(dat4)
#endif
# 419 "iotk_dat.spp"
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 420 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 420 "iotk_dat.spp"
call iotk_error_msg(ierr,'Kind incompatibility')
# 420 "iotk_dat.spp"
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(raw) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 426 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 426 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 426 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 432 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 432 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 432 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 440 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
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
# 451 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 451 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 451 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 456 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 461 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 461 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 461 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 466 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 466 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 466 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 472 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
          return
        end if
# 478 "iotk_dat.spp"
        if(index == size(dat)) exit
# 480 "iotk_dat.spp"
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 481 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 488 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 488 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 488 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
# 494 "iotk_dat.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 495 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
    return
  end if
end subroutine iotk_scan_dat_aux_INTEGER1
# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER1_1


#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_1
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_1

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 2 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_2(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_2


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER1_2


#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_2
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_2

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 3 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_3(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_3


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER1_3


#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_3
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_3

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 4 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_4(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_4


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER1_4


#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_4
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_4

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 5 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_5(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_5


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER1_5


#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_5
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_5

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 6 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_6(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_6


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER1_6


#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_6
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_6

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER1
#if 7 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER1_7(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER1), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER1(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER1_7


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER1_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER1),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER1_7


#endif
#endif

subroutine iotk_dat_dummy_INTEGER1_7
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER1_7

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER2
#if 0 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER2_0(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",1,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(1))
# 156 "iotk_dat.spp"
     dattmp(1) = dat
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER2_0


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER2_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==1) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(1))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 553 "iotk_dat.spp"
        dat = tmpdat(1)
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER2_0


#endif
#endif

subroutine iotk_dat_dummy_INTEGER2_0
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER2_0

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER2
#if 1 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER2_1(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER2(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER2_1


# 242 "iotk_dat.spp"
recursive subroutine iotk_scan_dat_aux_INTEGER2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only: iotk_read
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,         intent(in)  :: unit
  INTEGER (kind=__IOTK_INTEGER2), intent(out) :: dat (:)
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
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER1
  INTEGER (__IOTK_INTEGER1), allocatable :: dat1 (:)
#endif
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER3
  INTEGER (__IOTK_INTEGER3), allocatable :: dat3 (:)
#endif
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER4
  INTEGER (__IOTK_INTEGER4), allocatable :: dat4 (:)
#endif
# 274 "iotk_dat.spp"
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
# 285 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
    return
  end if
# 375 "iotk_dat.spp"
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
        read(lunit,iostat=iostat) dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 381 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 381 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 381 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 387 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 387 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 387 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER1
    case(kind(dat1))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat1(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat1
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat1
# 415 "iotk_dat.spp"
      deallocate(dat1)
#endif
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER3
    case(kind(dat3))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat3(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat3
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat3
# 415 "iotk_dat.spp"
      deallocate(dat3)
#endif
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER4
    case(kind(dat4))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat4(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat4
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat4
# 415 "iotk_dat.spp"
      deallocate(dat4)
#endif
# 419 "iotk_dat.spp"
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 420 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 420 "iotk_dat.spp"
call iotk_error_msg(ierr,'Kind incompatibility')
# 420 "iotk_dat.spp"
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(raw) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 426 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 426 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 426 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 432 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 432 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 432 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 440 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
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
# 451 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 451 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 451 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 456 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 461 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 461 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 461 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 466 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 466 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 466 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 472 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
          return
        end if
# 478 "iotk_dat.spp"
        if(index == size(dat)) exit
# 480 "iotk_dat.spp"
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 481 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 488 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 488 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 488 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
# 494 "iotk_dat.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 495 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
    return
  end if
end subroutine iotk_scan_dat_aux_INTEGER2
# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER2_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER2_1


#endif
#endif

subroutine iotk_dat_dummy_INTEGER2_1
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER2_1

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER2
#if 2 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER2_2(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER2(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER2_2


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER2_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER2_2


#endif
#endif

subroutine iotk_dat_dummy_INTEGER2_2
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER2_2

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER2
#if 3 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER2_3(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER2(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER2_3


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER2_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER2_3


#endif
#endif

subroutine iotk_dat_dummy_INTEGER2_3
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER2_3

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER2
#if 4 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER2_4(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER2(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER2_4


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER2_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER2_4


#endif
#endif

subroutine iotk_dat_dummy_INTEGER2_4
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER2_4

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER2
#if 5 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER2_5(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER2(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER2_5


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER2_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER2_5


#endif
#endif

subroutine iotk_dat_dummy_INTEGER2_5
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER2_5

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER2
#if 6 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER2_6(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER2(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER2_6


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER2_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER2_6


#endif
#endif

subroutine iotk_dat_dummy_INTEGER2_6
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER2_6

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER2
#if 7 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER2_7(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER2), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER2(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER2_7


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER2_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER2),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER2_7


#endif
#endif

subroutine iotk_dat_dummy_INTEGER2_7
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER2_7

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER3
#if 0 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER3_0(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",1,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(1))
# 156 "iotk_dat.spp"
     dattmp(1) = dat
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER3_0


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER3_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==1) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(1))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 553 "iotk_dat.spp"
        dat = tmpdat(1)
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER3_0


#endif
#endif

subroutine iotk_dat_dummy_INTEGER3_0
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER3_0

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER3
#if 1 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER3_1(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER3(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER3_1


# 242 "iotk_dat.spp"
recursive subroutine iotk_scan_dat_aux_INTEGER3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only: iotk_read
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,         intent(in)  :: unit
  INTEGER (kind=__IOTK_INTEGER3), intent(out) :: dat (:)
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
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER1
  INTEGER (__IOTK_INTEGER1), allocatable :: dat1 (:)
#endif
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER2
  INTEGER (__IOTK_INTEGER2), allocatable :: dat2 (:)
#endif
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER4
  INTEGER (__IOTK_INTEGER4), allocatable :: dat4 (:)
#endif
# 274 "iotk_dat.spp"
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
# 285 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
    return
  end if
# 375 "iotk_dat.spp"
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
        read(lunit,iostat=iostat) dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 381 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 381 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 381 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 387 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 387 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 387 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER1
    case(kind(dat1))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat1(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat1
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat1
# 415 "iotk_dat.spp"
      deallocate(dat1)
#endif
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER2
    case(kind(dat2))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat2(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat2
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat2
# 415 "iotk_dat.spp"
      deallocate(dat2)
#endif
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER4
    case(kind(dat4))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat4(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat4
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat4
# 415 "iotk_dat.spp"
      deallocate(dat4)
#endif
# 419 "iotk_dat.spp"
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 420 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 420 "iotk_dat.spp"
call iotk_error_msg(ierr,'Kind incompatibility')
# 420 "iotk_dat.spp"
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(raw) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 426 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 426 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 426 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 432 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 432 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 432 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 440 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
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
# 451 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 451 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 451 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 456 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 461 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 461 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 461 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 466 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 466 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 466 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 472 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
          return
        end if
# 478 "iotk_dat.spp"
        if(index == size(dat)) exit
# 480 "iotk_dat.spp"
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 481 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 488 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 488 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 488 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
# 494 "iotk_dat.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 495 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
    return
  end if
end subroutine iotk_scan_dat_aux_INTEGER3
# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER3_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER3_1


#endif
#endif

subroutine iotk_dat_dummy_INTEGER3_1
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER3_1

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER3
#if 2 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER3_2(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER3(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER3_2


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER3_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER3_2


#endif
#endif

subroutine iotk_dat_dummy_INTEGER3_2
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER3_2

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER3
#if 3 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER3_3(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER3(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER3_3


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER3_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER3_3


#endif
#endif

subroutine iotk_dat_dummy_INTEGER3_3
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER3_3

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER3
#if 4 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER3_4(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER3(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER3_4


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER3_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER3_4


#endif
#endif

subroutine iotk_dat_dummy_INTEGER3_4
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER3_4

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER3
#if 5 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER3_5(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER3(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER3_5


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER3_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER3_5


#endif
#endif

subroutine iotk_dat_dummy_INTEGER3_5
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER3_5

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER3
#if 6 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER3_6(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER3(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER3_6


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER3_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER3_6


#endif
#endif

subroutine iotk_dat_dummy_INTEGER3_6
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER3_6

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER3
#if 7 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER3_7(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER3), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER3(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER3_7


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER3_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER3),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER3_7


#endif
#endif

subroutine iotk_dat_dummy_INTEGER3_7
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER3_7

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER4
#if 0 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER4_0(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: dat  
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",1,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(1))
# 156 "iotk_dat.spp"
     dattmp(1) = dat
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER4_0


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER4_0(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==1) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(1))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 553 "iotk_dat.spp"
        dat = tmpdat(1)
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER4_0


#endif
#endif

subroutine iotk_dat_dummy_INTEGER4_0
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER4_0

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER4
#if 1 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER4_1(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: dat (:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER4(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER4_1


# 242 "iotk_dat.spp"
recursive subroutine iotk_scan_dat_aux_INTEGER4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only: iotk_read
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,         intent(in)  :: unit
  INTEGER (kind=__IOTK_INTEGER4), intent(out) :: dat (:)
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
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER1
  INTEGER (__IOTK_INTEGER1), allocatable :: dat1 (:)
#endif
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER2
  INTEGER (__IOTK_INTEGER2), allocatable :: dat2 (:)
#endif
# 268 "iotk_dat.spp"
#ifdef __IOTK_INTEGER3
  INTEGER (__IOTK_INTEGER3), allocatable :: dat3 (:)
#endif
# 274 "iotk_dat.spp"
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
# 285 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
    return
  end if
# 375 "iotk_dat.spp"
  if(binary) then
    select case(rkind)
    case(kind(dat))
      if(raw) then
        read(lunit,iostat=iostat) dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 381 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 381 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 381 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      else
        read(lunit,iostat=iostat) idummy,dat
        if(iostat/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 387 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 387 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 387 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
          return
        end if
      end if
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER1
    case(kind(dat1))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat1(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat1
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat1
# 415 "iotk_dat.spp"
      deallocate(dat1)
#endif
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER2
    case(kind(dat2))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat2(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat2
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat2
# 415 "iotk_dat.spp"
      deallocate(dat2)
#endif
# 393 "iotk_dat.spp"
#ifdef __IOTK_INTEGER3
    case(kind(dat3))
      ! Giusto per scrupolo. Se e' raw non ci sono info sul kind, quindi questa linea e' irraggiungibile
      if(raw) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 397 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
        return
      end if
      allocate(dat3(ubound(dat,1)))
      read(lunit,iostat=iostat) idummy,dat3
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 403 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 403 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
# 413 "iotk_dat.spp"
      dat = dat3
# 415 "iotk_dat.spp"
      deallocate(dat3)
#endif
# 419 "iotk_dat.spp"
    case default
      call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 420 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 420 "iotk_dat.spp"
call iotk_error_msg(ierr,'Kind incompatibility')
# 420 "iotk_dat.spp"
call iotk_error_write(ierr,"kind",rkind)
    end select
  else
    if(raw) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 426 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 426 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 426 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"*")) then
      read(lunit,fmt=*,iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 432 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 432 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 432 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    else if(iotk_strcomp(fmt,"!")) then
      index = 0
      do
        call iotk_getline(lunit,line,length,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 440 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
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
# 451 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 451 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 451 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          call iotk_getline(lunit,altline,altlength,ierr)
          if(ierr/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 456 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
            return
          end if
          backspace(lunit,iostat=iostat)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 461 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 461 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 461 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
          read(lunit,"(a)",advance="no",iostat=iostat) altline(1:nexttag-1 + altlength - length)
          if(iostat/=0) then
            call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 466 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 466 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 466 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
            return
          end if
        end if
        call iotk_read(dat,line(1:nexttag - 1),index,ierr)
        if(ierr/=0) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 472 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
          return
        end if
# 478 "iotk_dat.spp"
        if(index == size(dat)) exit
# 480 "iotk_dat.spp"
        if(nexttag/=length + 1) then
          call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 481 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
          return
        end if
      end do
    else
      read(lunit,fmt=fmt(1:iotk_strlen(fmt)),iostat=iostat) dat
      if(iostat/=0) then
        call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 488 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
# 488 "iotk_dat.spp"
call iotk_error_msg(ierr,' ')
# 488 "iotk_dat.spp"
call iotk_error_write(ierr,"iostat",iostat)
        return
      end if
    end if
  end if
# 494 "iotk_dat.spp"
  if(idummy/=0) then
    call iotk_error_issue(ierr,"iotk_scan_dat_aux",__FILE__,__LINE__)
# 495 "iotk_dat.spp"
call iotk_error_msg(ierr,"CVS Revision: 1.2 ")
    return
  end if
end subroutine iotk_scan_dat_aux_INTEGER4
# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER4_1(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER4_1


#endif
#endif

subroutine iotk_dat_dummy_INTEGER4_1
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER4_1

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER4
#if 2 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER4_2(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: dat (:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER4(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER4_2


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER4_2(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER4_2


#endif
#endif

subroutine iotk_dat_dummy_INTEGER4_2
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER4_2

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER4
#if 3 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER4_3(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: dat (:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER4(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER4_3


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER4_3(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER4_3


#endif
#endif

subroutine iotk_dat_dummy_INTEGER4_3
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER4_3

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER4
#if 4 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER4_4(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: dat (:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER4(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER4_4


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER4_4(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER4_4


#endif
#endif

subroutine iotk_dat_dummy_INTEGER4_4
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER4_4

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER4
#if 5 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER4_5(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: dat (:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER4(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER4_5


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER4_5(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER4_5


#endif
#endif

subroutine iotk_dat_dummy_INTEGER4_5
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER4_5

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER4
#if 6 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER4_6(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: dat (:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER4(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER4_6


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER4_6(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER4_6


#endif
#endif

subroutine iotk_dat_dummy_INTEGER4_6
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER4_6

# 45 "iotk_dat.spp"

# 47 "iotk_dat.spp"

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

# 57 "iotk_dat.spp"

# 59 "iotk_dat.spp"

#ifdef __IOTK_INTEGER4
#if 7 <= __IOTK_MAXRANK
# 63 "iotk_dat.spp"
subroutine iotk_write_dat_INTEGER4_7(unit,name,dat,fmt,ierr)
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
  INTEGER (kind=__IOTK_INTEGER4), intent(in)  :: dat (:,:,:,:,:,:,:) 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
  integer :: ierrl,lunit,iostat
  logical :: binary,raw
  integer(iotk_header_kind), parameter :: idummy=0
  character(300) :: usefmt,usefmt1
  character(iotk_attlenx) :: attr
  type (iotk_unit), pointer :: this
# 88 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),allocatable :: dattmp(:)
# 90 "iotk_dat.spp"
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
# 101 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_strcpy(usefmt,"!",ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 106 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(present(fmt) .and. .not. raw) call iotk_strcpy(usefmt,iotk_strtrim(fmt),ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 111 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(iotk_strscan(usefmt,"<>&")/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 115 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Special characters (<>&) found in fmt string')
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"unit",unit)
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"name",trim(name))
# 115 "iotk_dat.spp"
call iotk_error_write(ierrl,"fmt",trim(fmt))
    goto 1
  end if
  call iotk_write_attr(attr,"type",iotk_tolower("INTEGER"),first=.true.,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 120 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_attr(attr,"size",size(dat),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 125 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
# 135 "iotk_dat.spp"
  if(binary) then
    call iotk_write_attr(attr,"kind",kind(dat),ierr=ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 138 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
    end if
  end if
# 143 "iotk_dat.spp"
  if(.not.iotk_strcomp(usefmt,"!")) call iotk_write_attr(attr,"fmt",iotk_strtrim(usefmt),ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 145 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  call iotk_write_begin(unit,name,attr,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 150 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if

  allocate(dattmp(size(dat)))
# 158 "iotk_dat.spp"
#if defined(__IOTK_WORKAROUND3) || defined(__IOTK_WORKAROUND4)
# 162 "iotk_dat.spp"
     call iotk_private_pack_INTEGER4(dattmp,dat,size(dattmp),1)
# 164 "iotk_dat.spp"
#else
     dattmp = pack(dat,mask=.true.)
#endif
# 168 "iotk_dat.spp"

  if(binary) then
    if(raw) then
      write(lunit,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 173 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else
      write(lunit,iostat=iostat) idummy,(dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 179 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  else
    if(raw) then
# 188 "iotk_dat.spp"
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
# 190 "iotk_dat.spp"
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 191 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"*")) then
      write(lunit,*,iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 197 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    else if(iotk_strcomp(usefmt,"!")) then
# 211 "iotk_dat.spp"
     write(lunit,fmt=iotk_wfmt("INTEGER",kind(dattmp),1,-1),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
     if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 213 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
      goto 1
     end if
# 217 "iotk_dat.spp"
    else
      write(lunit,fmt=usefmt(1:iotk_strlen(usefmt)),iostat=iostat) (dattmp(itmp),itmp=1,size(dattmp))
      if(iostat/=0) then
        call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 220 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
        goto 1
      end if
    end if
  end if
  call iotk_write_end(unit,name,ierr=ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_dat",__FILE__,__LINE__)
# 227 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
1 continue
  if(allocated(dattmp)) deallocate(dattmp)
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_dat_INTEGER4_7


# 500 "iotk_dat.spp"

# 502 "iotk_dat.spp"
subroutine iotk_scan_dat_INTEGER4_7(unit,name,dat,found,default,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_dat_interf, only: iotk_scan_dat_aux
  use iotk_scan_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
# 519 "iotk_dat.spp"
  INTEGER (kind=__IOTK_INTEGER4),              allocatable :: tmpdat(:)
# 521 "iotk_dat.spp"
  integer :: ierrl,ierrl2
  integer :: rkind,rsize,rlen
  character(iotk_vallenx) :: rtype
  character(iotk_vallenx) :: fmt
  character(iotk_attlenx) :: attr
  logical :: inside,foundl
  inside = .false.
  ierrl = 0
  ierrl2 = 0
  foundl=.false.
  call iotk_scan_begin(unit,name,attr,found=foundl,ierr=ierrl)
  if(.not. foundl) goto 1
  foundl = .true.
  inside = .true.
  call iotk_parse_dat(attr,rtype,rkind,rsize,rlen,fmt,ierrl)
  if(ierrl/=0) goto 1
  if(.not. (iotk_strcomp(rtype,iotk_eos) .or. iotk_strcomp(rtype,"INTEGER") ) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 538 "iotk_dat.spp"
call iotk_error_msg(ierrl,' ')
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"rtype",rtype(1:iotk_strlen(rtype)))
# 538 "iotk_dat.spp"
call iotk_error_write(ierrl,"type","INTEGER")
    goto 1
  end if
  if(.not. (rsize==-1 .or. rsize==size(dat)) ) then
    call iotk_error_issue(ierrl,"iotk_scan_dat",__FILE__,__LINE__)
# 542 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
    goto 1
  end if
  if(rkind==-1) rkind = kind(dat)
# 549 "iotk_dat.spp"

  allocate(tmpdat(size(dat)))
  call iotk_scan_dat_aux(unit,tmpdat,rkind,rlen,fmt(1:iotk_strlen(fmt)),ierrl)
# 555 "iotk_dat.spp"
        dat = reshape(tmpdat,shape(dat))
# 557 "iotk_dat.spp"
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
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,"CVS Revision: 1.2 ")
# 569 "iotk_dat.spp"
call iotk_error_msg(ierrl,'Dat not found')
# 569 "iotk_dat.spp"
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
end subroutine iotk_scan_dat_INTEGER4_7


#endif
#endif

subroutine iotk_dat_dummy_INTEGER4_7
  write(0,*)
end subroutine iotk_dat_dummy_INTEGER4_7

# 45 "iotk_dat.spp"

