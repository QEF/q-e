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


function iotk_atol_x(a,check)
  use iotk_base
  use iotk_misc_interf
  implicit none
  character(len=*),           intent(in)  :: a
  logical,          optional, intent(out) :: check
  logical :: iotk_atol_x
  integer :: i
  iotk_atol_x = .false.
  if(present(check)) check = .false.
  if(len(a)==0) return
  do i = 1 , len(a)
    if(a(i:i)/=" " .and. a(i:i)/=".") exit
  end do
  if(i>len(a)) return
  if(present(check)) check = .true.
  if(a(i:i)=="T" .or. a(i:i)=="t") then
    iotk_atol_x = .true.
  else if(a(i:i)=="F" .or. a(i:i)=="f") then
    iotk_atol_x = .false.
  else
    if(present(check)) check = .false.
  end if
end function iotk_atol_x

#ifdef __IOTK_INTEGER1
subroutine iotk_atoi1(i,a,check)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer1
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
  logical :: minus
  integer :: pos,ii
  integer(kind=this_kind) :: j
  integer :: index
#ifdef __IOTK_WORKAROUND5
  integer(kind=this_kind) :: limit(0:9)
  integer(kind=this_kind) :: hug
  hug = huge(j)
  limit(0:9) = (/ ((hug-j)/10,j=0,9) /)
#else
  integer(kind=this_kind), parameter :: limit(0:9) = (/ ((huge(j)-j)/10_this_kind,j=0,9) /) 
#endif
  minus = .false.
  i = 0
  if(present(check)) check = .false.
  if(len(a)==0) return
  do ii = 1 , len(a)
    if(a(ii:ii)/=" ") exit
  end do
  if(ii>len(a)) return
  if(a(ii:ii)=="-") then
    minus = .true.
    ii = ii + 1
  else if(a(ii:ii)=="+") then
    ii = ii + 1
  end if
  if(ii>len(a)) return
  pos = ii
  do ii=pos,len(a)
    index = int( iachar(a(ii:ii)) - iachar("0") )
    if(index<0 .or. index>9) exit
    if(i>limit(index)) exit ! Check sull'overflow
    i = i*10_this_kind + int(index,kind=this_kind)
  end do
  if(minus) i = - i
  if(present(check)) then
    pos = ii
    do ii=pos,len(a)
      if(a(ii:ii)/=" ") return
    end do
    check = .true.
  end if
end subroutine iotk_atoi1
#endif
#ifdef __IOTK_INTEGER2
subroutine iotk_atoi2(i,a,check)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer2
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
  logical :: minus
  integer :: pos,ii
  integer(kind=this_kind) :: j
  integer :: index
#ifdef __IOTK_WORKAROUND5
  integer(kind=this_kind) :: limit(0:9)
  integer(kind=this_kind) :: hug
  hug = huge(j)
  limit(0:9) = (/ ((hug-j)/10,j=0,9) /)
#else
  integer(kind=this_kind), parameter :: limit(0:9) = (/ ((huge(j)-j)/10_this_kind,j=0,9) /) 
#endif
  minus = .false.
  i = 0
  if(present(check)) check = .false.
  if(len(a)==0) return
  do ii = 1 , len(a)
    if(a(ii:ii)/=" ") exit
  end do
  if(ii>len(a)) return
  if(a(ii:ii)=="-") then
    minus = .true.
    ii = ii + 1
  else if(a(ii:ii)=="+") then
    ii = ii + 1
  end if
  if(ii>len(a)) return
  pos = ii
  do ii=pos,len(a)
    index = int( iachar(a(ii:ii)) - iachar("0") )
    if(index<0 .or. index>9) exit
    if(i>limit(index)) exit ! Check sull'overflow
    i = i*10_this_kind + int(index,kind=this_kind)
  end do
  if(minus) i = - i
  if(present(check)) then
    pos = ii
    do ii=pos,len(a)
      if(a(ii:ii)/=" ") return
    end do
    check = .true.
  end if
end subroutine iotk_atoi2
#endif
#ifdef __IOTK_INTEGER3
subroutine iotk_atoi3(i,a,check)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer3
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
  logical :: minus
  integer :: pos,ii
  integer(kind=this_kind) :: j
  integer :: index
#ifdef __IOTK_WORKAROUND5
  integer(kind=this_kind) :: limit(0:9)
  integer(kind=this_kind) :: hug
  hug = huge(j)
  limit(0:9) = (/ ((hug-j)/10,j=0,9) /)
#else
  integer(kind=this_kind), parameter :: limit(0:9) = (/ ((huge(j)-j)/10_this_kind,j=0,9) /) 
#endif
  minus = .false.
  i = 0
  if(present(check)) check = .false.
  if(len(a)==0) return
  do ii = 1 , len(a)
    if(a(ii:ii)/=" ") exit
  end do
  if(ii>len(a)) return
  if(a(ii:ii)=="-") then
    minus = .true.
    ii = ii + 1
  else if(a(ii:ii)=="+") then
    ii = ii + 1
  end if
  if(ii>len(a)) return
  pos = ii
  do ii=pos,len(a)
    index = int( iachar(a(ii:ii)) - iachar("0") )
    if(index<0 .or. index>9) exit
    if(i>limit(index)) exit ! Check sull'overflow
    i = i*10_this_kind + int(index,kind=this_kind)
  end do
  if(minus) i = - i
  if(present(check)) then
    pos = ii
    do ii=pos,len(a)
      if(a(ii:ii)/=" ") return
    end do
    check = .true.
  end if
end subroutine iotk_atoi3
#endif
#ifdef __IOTK_INTEGER4
subroutine iotk_atoi4(i,a,check)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer4
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
  logical :: minus
  integer :: pos,ii
  integer(kind=this_kind) :: j
  integer :: index
#ifdef __IOTK_WORKAROUND5
  integer(kind=this_kind) :: limit(0:9)
  integer(kind=this_kind) :: hug
  hug = huge(j)
  limit(0:9) = (/ ((hug-j)/10,j=0,9) /)
#else
  integer(kind=this_kind), parameter :: limit(0:9) = (/ ((huge(j)-j)/10_this_kind,j=0,9) /) 
#endif
  minus = .false.
  i = 0
  if(present(check)) check = .false.
  if(len(a)==0) return
  do ii = 1 , len(a)
    if(a(ii:ii)/=" ") exit
  end do
  if(ii>len(a)) return
  if(a(ii:ii)=="-") then
    minus = .true.
    ii = ii + 1
  else if(a(ii:ii)=="+") then
    ii = ii + 1
  end if
  if(ii>len(a)) return
  pos = ii
  do ii=pos,len(a)
    index = int( iachar(a(ii:ii)) - iachar("0") )
    if(index<0 .or. index>9) exit
    if(i>limit(index)) exit ! Check sull'overflow
    i = i*10_this_kind + int(index,kind=this_kind)
  end do
  if(minus) i = - i
  if(present(check)) then
    pos = ii
    do ii=pos,len(a)
      if(a(ii:ii)/=" ") return
    end do
    check = .true.
  end if
end subroutine iotk_atoi4
#endif

#ifdef __IOTK_INTEGER1
function iotk_itoa1(i,length)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer1
  integer(kind=this_kind),           intent(in)  :: i
  integer,                 optional, intent(out) :: length
  character(len=range(i)+2)                      :: iotk_itoa1
  integer(kind=this_kind) :: itmp
  integer :: pos,pos1
  character(len=range(i)+2) :: tmp
  itmp = abs(i)
  do pos = 1 , len(tmp)
    tmp(pos:pos) = achar( modulo(itmp,int(10,kind(itmp))) + iachar("0") )
    itmp = itmp/10_this_kind
    if(itmp==0) exit
    if(pos==len(tmp)) exit
  end do
  if(i<0) then
    tmp(pos+1:pos+1)="-"
    pos = pos + 1
  end if
  do pos1=1,pos
    iotk_itoa1(pos1:pos1) = tmp(pos-pos1+1:pos-pos1+1)
  end do
  if(present(length)) length = pos
  do pos1=pos+1,len(iotk_itoa1)
    iotk_itoa1(pos1:pos1) = " "
  end do
end function iotk_itoa1
#endif
#ifdef __IOTK_INTEGER2
function iotk_itoa2(i,length)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer2
  integer(kind=this_kind),           intent(in)  :: i
  integer,                 optional, intent(out) :: length
  character(len=range(i)+2)                      :: iotk_itoa2
  integer(kind=this_kind) :: itmp
  integer :: pos,pos1
  character(len=range(i)+2) :: tmp
  itmp = abs(i)
  do pos = 1 , len(tmp)
    tmp(pos:pos) = achar( modulo(itmp,int(10,kind(itmp))) + iachar("0") )
    itmp = itmp/10_this_kind
    if(itmp==0) exit
    if(pos==len(tmp)) exit
  end do
  if(i<0) then
    tmp(pos+1:pos+1)="-"
    pos = pos + 1
  end if
  do pos1=1,pos
    iotk_itoa2(pos1:pos1) = tmp(pos-pos1+1:pos-pos1+1)
  end do
  if(present(length)) length = pos
  do pos1=pos+1,len(iotk_itoa2)
    iotk_itoa2(pos1:pos1) = " "
  end do
end function iotk_itoa2
#endif
#ifdef __IOTK_INTEGER3
function iotk_itoa3(i,length)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer3
  integer(kind=this_kind),           intent(in)  :: i
  integer,                 optional, intent(out) :: length
  character(len=range(i)+2)                      :: iotk_itoa3
  integer(kind=this_kind) :: itmp
  integer :: pos,pos1
  character(len=range(i)+2) :: tmp
  itmp = abs(i)
  do pos = 1 , len(tmp)
    tmp(pos:pos) = achar( modulo(itmp,int(10,kind(itmp))) + iachar("0") )
    itmp = itmp/10_this_kind
    if(itmp==0) exit
    if(pos==len(tmp)) exit
  end do
  if(i<0) then
    tmp(pos+1:pos+1)="-"
    pos = pos + 1
  end if
  do pos1=1,pos
    iotk_itoa3(pos1:pos1) = tmp(pos-pos1+1:pos-pos1+1)
  end do
  if(present(length)) length = pos
  do pos1=pos+1,len(iotk_itoa3)
    iotk_itoa3(pos1:pos1) = " "
  end do
end function iotk_itoa3
#endif
#ifdef __IOTK_INTEGER4
function iotk_itoa4(i,length)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer4
  integer(kind=this_kind),           intent(in)  :: i
  integer,                 optional, intent(out) :: length
  character(len=range(i)+2)                      :: iotk_itoa4
  integer(kind=this_kind) :: itmp
  integer :: pos,pos1
  character(len=range(i)+2) :: tmp
  itmp = abs(i)
  do pos = 1 , len(tmp)
    tmp(pos:pos) = achar( modulo(itmp,int(10,kind(itmp))) + iachar("0") )
    itmp = itmp/10_this_kind
    if(itmp==0) exit
    if(pos==len(tmp)) exit
  end do
  if(i<0) then
    tmp(pos+1:pos+1)="-"
    pos = pos + 1
  end if
  do pos1=1,pos
    iotk_itoa4(pos1:pos1) = tmp(pos-pos1+1:pos-pos1+1)
  end do
  if(present(length)) length = pos
  do pos1=pos+1,len(iotk_itoa4)
    iotk_itoa4(pos1:pos1) = " "
  end do
end function iotk_itoa4
#endif

#ifdef __IOTK_LOGICAL1
function iotk_ltoa1(l)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer1
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa1
  if(l) then
    iotk_ltoa1 = "T"
  else
    iotk_ltoa1 = "F"
  end if
end function iotk_ltoa1
#endif
#ifdef __IOTK_LOGICAL2
function iotk_ltoa2(l)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer2
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa2
  if(l) then
    iotk_ltoa2 = "T"
  else
    iotk_ltoa2 = "F"
  end if
end function iotk_ltoa2
#endif
#ifdef __IOTK_LOGICAL3
function iotk_ltoa3(l)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer3
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa3
  if(l) then
    iotk_ltoa3 = "T"
  else
    iotk_ltoa3 = "F"
  end if
end function iotk_ltoa3
#endif
#ifdef __IOTK_LOGICAL4
function iotk_ltoa4(l)
  use iotk_base
  use iotk_misc_interf
  implicit none
  integer, parameter :: this_kind = iotk_integer4
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa4
  if(l) then
    iotk_ltoa4 = "T"
  else
    iotk_ltoa4 = "F"
  end if
end function iotk_ltoa4
#endif
