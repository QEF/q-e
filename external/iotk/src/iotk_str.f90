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


function iotk_toupper_x(str)
  use iotk_base
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: iotk_toupper_x
  integer :: i,pos
  do i = 1,len(str)
    if(str(i:i)==iotk_eos) exit
    pos=scan(lowalphabet,str(i:i))
    if(pos==0) then
      iotk_toupper_x(i:i) = str(i:i)
    else
      iotk_toupper_x(i:i) = upalphabet(pos:pos)
    end if
  end do
  if(i<=len(iotk_toupper_x)) iotk_toupper_x(i:i) = iotk_eos
end function iotk_toupper_x

function iotk_tolower_x(str)
  use iotk_base
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: iotk_tolower_x
  integer :: i,pos
  do i = 1,len(str)
    if(str(i:i)==iotk_eos) exit
    pos=scan(upalphabet,str(i:i))
    if(pos==0) then
      iotk_tolower_x(i:i) = str(i:i)
    else
      iotk_tolower_x(i:i) = lowalphabet(pos:pos)
    end if
  end do
  if(i<=len(iotk_tolower_x)) iotk_tolower_x(i:i) = iotk_eos
end function iotk_tolower_x

subroutine iotk_escape_x(to,from)
  use iotk_base
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  character(len=*), intent(in)  :: from
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: to
#else
  character(len=*), intent(out) :: to
#endif
  integer :: pos,pos1,semic,fromlen
  pos = 1
  pos1 = 1
  fromlen = iotk_strlen(from)
  do  
    if(pos>fromlen) exit
    if(from(pos:pos)=="&" .and. pos/=fromlen) then
      semic = scan(from(pos+1:fromlen),";")
      if(semic<=1) to(pos1:pos1)="&"
      select case(from(pos+1:pos+semic-1))
      case("amp")
        to(pos1:pos1)="&"
      case("lt")
        to(pos1:pos1)="<"
      case("gt")
        to(pos1:pos1)=">"
      case("quot")
        to(pos1:pos1)='"'
      case("apos")
        to(pos1:pos1)="'"
      case default
        to(pos1:pos1+semic) = from(pos:pos+semic)
        pos1 = pos1 + semic
      end select
      pos = pos + semic
    else
      to(pos1:pos1)=from(pos:pos)
    end if
    pos = pos + 1
    pos1 = pos1 + 1
    if(pos1>len(to)) exit
  end do  
  if(pos1<=len(to)) to(pos1:pos1)=iotk_eos
end subroutine iotk_escape_x

subroutine iotk_deescape_x(to,from,quot,apos)
  use iotk_base
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  character(len=*), intent(in)  :: from
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: to
#else
  character(len=*), intent(out) :: to
#endif
  logical, optional, intent(in) :: quot,apos
  logical :: lquot,lapos
  integer :: pos,pos1
  lquot=.false.
  lapos=.false.
  if(present(quot)) lquot = quot
  if(present(apos)) lapos = apos
  pos = 1
  pos1 = 1
  do
    if(pos>len(from) .or. pos1>len(to)) exit ! The two checks must be separated
    if(from(pos:pos)==iotk_eos) exit
    select case(from(pos:pos))
    case("&")
      if(pos1+4<=len(to)) to(pos1:pos1+4)="&amp;"
      pos1=pos1+4
    case("<")
      if(pos1+3<=len(to)) to(pos1:pos1+3)="&lt;"
      pos1=pos1+3
    case(">")
      if(pos1+3<=len(to)) to(pos1:pos1+3)="&gt;"
      pos1=pos1+3
    case('"')
      if(lquot) then
        if(pos1+5<=len(to)) to(pos1:pos1+5)="&quot;"
        pos1=pos1+5
      else
        to(pos1:pos1) = from(pos:pos)
      end if
    case("'")
      if(lapos) then
        if(pos1+5<=len(to)) to(pos1:pos1+5)="&apos;"
        pos1=pos1+5
      else
        to(pos1:pos1) = from(pos:pos)
      end if
    case default
      to(pos1:pos1) = from(pos:pos)
    end select
    pos = pos + 1
    pos1 = pos1 + 1
  end do
  if(pos1<=len(to)) to(pos1:pos1)=iotk_eos
end subroutine iotk_deescape_x

function iotk_strtrim_x(str)
  use iotk_base
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: iotk_strtrim_x
  integer :: lentrim
  lentrim = len_trim(str(1:iotk_strlen(str)))
  iotk_strtrim_x(1:lentrim) = str(1:lentrim)
  if(lentrim<len(iotk_strtrim_x)) iotk_strtrim_x(lentrim+1:lentrim+1) = iotk_eos
end function iotk_strtrim_x

function iotk_strlen_trim_x(str)
  use iotk_base
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in) :: str
  integer                      :: iotk_strlen_trim_x
  iotk_strlen_trim_x = len_trim(str(1:iotk_strlen(str)))
end function iotk_strlen_trim_x

function iotk_strscan_x(string,set,back)
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  character(len=*),  intent(in) :: string
  character(len=*),  intent(in) :: set
  logical, optional, intent(in) :: back
  integer                       :: iotk_strscan_x
  logical :: backl
  backl = .false.
  if(present(back)) backl=back
  iotk_strscan_x = scan(string(1:iotk_strlen(string)),set(1:iotk_strlen(set)),backl)
end function iotk_strscan_x

function iotk_strlen_x(str)
  use iotk_base
  implicit none
  character(len=*), intent(in) :: str
  integer :: iotk_strlen_x
  integer :: pos
  pos = scan(str,iotk_eos) - 1
  if(pos>=0) then
    iotk_strlen_x = pos
  else
    iotk_strlen_x = len(str)
  end if
end function iotk_strlen_x

function iotk_strpad_x(str)
  use iotk_base
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: iotk_strpad_x
  integer :: strlen
  strlen = iotk_strlen(str)
  iotk_strpad_x(1:strlen) = str(1:strlen)
  if(strlen<len(iotk_strpad_x)) iotk_strpad_x(strlen+1:) = " "
end function iotk_strpad_x

subroutine iotk_strcpy_x(to,from,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: to
#else
  character(len=*), intent(out) :: to
#endif
  character(len=*), intent(in)  :: from
  integer,          intent(out) :: ierr
  integer :: i
  ierr = 0
  do i=1,min(len(from),len(to))
    if(from(i:i)==iotk_eos) exit
    to(i:i)=from(i:i)
  end do
  if(i>len(to) .and. i<=len(from)) then
    if(from(i:i)/=iotk_eos) then
      call iotk_error_issue(ierr,"iotk_strcpy",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.14 ")
      return
    end if
  end if
  if(i<=len(to)) to(i:i) = iotk_eos
end subroutine iotk_strcpy_x

subroutine iotk_strcat_x(to,from,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(len=*), intent(inout):: to
  character(len=*), intent(in) :: from
  integer,          intent(out):: ierr
  integer :: tolen,fromlen
  ierr = 0
  tolen = iotk_strlen(to)
  fromlen = iotk_strlen(from)
  if(tolen+fromlen>len(to)) then
    call iotk_error_issue(ierr,"iotk_strcat",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.14 ")
  end if
  if(ierr/=0) return
  to(tolen+1:tolen+fromlen) = from(1:fromlen)
  if(tolen+fromlen+1<=len(to)) to(tolen+fromlen+1:tolen+fromlen+1)=iotk_eos
end subroutine iotk_strcat_x

function iotk_strcomp_x(str1,str2)
  use iotk_base
  implicit none
  logical :: iotk_strcomp_x
  character(len=*), intent(in) :: str1,str2
  integer :: i
  iotk_strcomp_x = .false.
  do i=1,min(len(str1),len(str2))
    if(str1(i:i)/=str2(i:i)) return
    if(str1(i:i)==iotk_eos) exit
  end do
  if(i>len(str1)) then
    if(i<=len(str2)) then
      if(str2(i:i)/=iotk_eos) return
    end if
  else if(i>len(str2)) then
    if(i<=len(str1)) then
      if(str1(i:i)/=iotk_eos) return
    end if
  end if
  iotk_strcomp_x = .true.
end function iotk_strcomp_x

subroutine iotk_str_clean_x(str)
! transforms all characters which are separators in blanks
  use iotk_base
  implicit none
  character(len=*), intent(inout) :: str
  integer :: i
  do i = 1 , len(str)
    if(str(i:i)==iotk_eos) exit
    if(scan(not_separator,str(i:i))==0) str(i:i)=" "
  end do
end subroutine iotk_str_clean_x
