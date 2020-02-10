! Input/Output Tool Kit (IOTK)
! Copyright (C) 2006 Giovanni Bussi
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

module iotk_tool_interf
private
public :: iotk_tool
public :: iotk_tool_convert
public :: iotk_tool_dump
public :: iotk_tool_info
public :: iotk_tool_man

interface iotk_tool
subroutine iotk_tool_x(args)
  implicit none
  character(len=*), intent(in) :: args(:)
end subroutine iotk_tool_x
end interface

interface iotk_tool_convert
subroutine iotk_tool_convert_x(args,ierr)
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
end subroutine iotk_tool_convert_x
end interface

interface iotk_tool_dump
subroutine iotk_tool_dump_x(args,ierr)
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
end subroutine iotk_tool_dump_x
end interface

interface iotk_tool_info
subroutine iotk_tool_info_x(args,ierr)
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
end subroutine iotk_tool_info_x
end interface

interface iotk_tool_man
subroutine iotk_tool_man_x(args,ierr)
  implicit none
  character(len=*),           intent(in)  :: args(:)
  integer,          optional, intent(out) :: ierr
end subroutine iotk_tool_man_x
end interface

end module iotk_tool_interf
