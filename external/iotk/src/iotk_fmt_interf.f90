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

module iotk_fmt_interf
use iotk_base
implicit none
private

public :: iotk_basefmt
public :: iotk_wfmt

interface iotk_basefmt
function iotk_basefmt_x(type,ikind,ilen)
  implicit none
  character(len=*), intent(in) :: type
  integer,          intent(in) :: ikind
  integer,          intent(in) :: ilen
  character(100)               :: iotk_basefmt_x
end function iotk_basefmt_x
end interface

interface iotk_wfmt
function iotk_wfmt_x(type,ikind,isize,ilen,sep)
  implicit none
  character(len=*), intent(in)  :: type
  integer,          intent(in)  :: ikind
  integer,          intent(in)  :: isize
  integer,          intent(in)  :: ilen
  character(len=*), intent(in)  :: sep
  character(150)                :: iotk_wfmt_x
end function iotk_wfmt_x
end interface

end module iotk_fmt_interf
