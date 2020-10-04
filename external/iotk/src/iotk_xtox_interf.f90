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

module iotk_xtox_interf
use iotk_base
implicit none
private

public :: iotk_atol
public :: iotk_ltoa
public :: iotk_atoi
public :: iotk_itoa

interface iotk_atol
function iotk_atol_x(a,check)
  character(len=*),           intent(in)  :: a
  logical,          optional, intent(out) :: check
  logical                                 :: iotk_atol_x
end function iotk_atol_x
end interface

interface iotk_atoi
#ifdef __IOTK_INTEGER1
subroutine iotk_atoi1(i,a,check)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer1
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
end subroutine iotk_atoi1
#endif
#ifdef __IOTK_INTEGER2
subroutine iotk_atoi2(i,a,check)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer2
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
end subroutine iotk_atoi2
#endif
#ifdef __IOTK_INTEGER3
subroutine iotk_atoi3(i,a,check)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer3
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
end subroutine iotk_atoi3
#endif
#ifdef __IOTK_INTEGER4
subroutine iotk_atoi4(i,a,check)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer4
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
end subroutine iotk_atoi4
#endif
end interface

interface iotk_itoa
#ifdef __IOTK_INTEGER1
function iotk_itoa1(i,length)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer1
  integer(kind=this_kind),             intent(in)  :: i
  integer,                   optional, intent(out) :: length
  character(len=range(i)+2)                        :: iotk_itoa1
end function iotk_itoa1
#endif
#ifdef __IOTK_INTEGER2
function iotk_itoa2(i,length)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer2
  integer(kind=this_kind),             intent(in)  :: i
  integer,                   optional, intent(out) :: length
  character(len=range(i)+2)                        :: iotk_itoa2
end function iotk_itoa2
#endif
#ifdef __IOTK_INTEGER3
function iotk_itoa3(i,length)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer3
  integer(kind=this_kind),             intent(in)  :: i
  integer,                   optional, intent(out) :: length
  character(len=range(i)+2)                        :: iotk_itoa3
end function iotk_itoa3
#endif
#ifdef __IOTK_INTEGER4
function iotk_itoa4(i,length)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer4
  integer(kind=this_kind),             intent(in)  :: i
  integer,                   optional, intent(out) :: length
  character(len=range(i)+2)                        :: iotk_itoa4
end function iotk_itoa4
#endif
end interface

interface iotk_ltoa
#ifdef __IOTK_LOGICAL1
function iotk_ltoa1(l)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_logical1
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa1
end function iotk_ltoa1
#endif
#ifdef __IOTK_LOGICAL2
function iotk_ltoa2(l)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_logical2
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa2
end function iotk_ltoa2
#endif
#ifdef __IOTK_LOGICAL3
function iotk_ltoa3(l)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_logical3
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa3
end function iotk_ltoa3
#endif
#ifdef __IOTK_LOGICAL4
function iotk_ltoa4(l)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_logical4
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa4
end function iotk_ltoa4
#endif
end interface

end module iotk_xtox_interf
