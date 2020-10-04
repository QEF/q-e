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

module iotk_error_interf
implicit none
private

public :: iotk_error_init
public :: iotk_error_clear
public :: iotk_error_append
public :: iotk_error_add
public :: iotk_error_print
public :: iotk_error_issue
public :: iotk_error_check
public :: iotk_error_msg
public :: iotk_error_write
public :: iotk_error_scan
public :: iotk_error_handler
public :: iotk_error_pool_pending

interface iotk_error_init
subroutine iotk_error_init_e(error)
  use iotk_base
  implicit none
  type(iotk_error), intent(out) :: error
end subroutine iotk_error_init_e
subroutine iotk_error_init_i(ierr)
  use iotk_base
  implicit none
  integer, intent(out) :: ierr
end subroutine iotk_error_init_i
end interface

interface iotk_error_clear
subroutine iotk_error_clear_e(error)
  use iotk_base
  implicit none
  type(iotk_error), intent(inout) :: error
end subroutine iotk_error_clear_e
subroutine iotk_error_clear_i(ierr)
  use iotk_base
  implicit none
  integer, intent(inout) :: ierr
end subroutine iotk_error_clear_i
end interface

interface iotk_error_add
function iotk_error_add_x()
  use iotk_base
  implicit none
  integer :: iotk_error_add_x
end function iotk_error_add_x
end interface

interface iotk_error_append
subroutine iotk_error_append_e(error,str)
  use iotk_base
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: str
end subroutine iotk_error_append_e
subroutine iotk_error_append_i(ierr,str)
  use iotk_base
  implicit none
  integer, intent(inout) :: ierr
  character(len=*), intent(in)    :: str
end subroutine iotk_error_append_i
end interface

interface iotk_error_print
subroutine iotk_error_print_e(error,unit)
  use iotk_base
  implicit none
  type(iotk_error), intent(in) :: error
  integer,          intent(in) :: unit
end subroutine iotk_error_print_e
subroutine iotk_error_print_i(ierr,unit)
  use iotk_base
  implicit none
  integer, intent(in) :: ierr
  integer, intent(in) :: unit
end subroutine iotk_error_print_i
end interface

interface iotk_error_issue
subroutine iotk_error_issue_e(error,sub,file,line)
  use iotk_base
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: sub
  character(len=*), intent(in)    :: file
  integer,          intent(in)    :: line
end subroutine iotk_error_issue_e
subroutine iotk_error_issue_i(ierr,sub,file,line)
  use iotk_base
  implicit none
  integer,          intent(inout) :: ierr
  character(len=*), intent(in)    :: sub
  character(len=*), intent(in)    :: file
  integer,          intent(in)    :: line
end subroutine iotk_error_issue_i
end interface

interface iotk_error_msg
subroutine iotk_error_msg_e(error,msg)
  use iotk_base
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: msg
end subroutine iotk_error_msg_e
subroutine iotk_error_msg_i(ierr,msg)
  use iotk_base
  implicit none
  integer,          intent(inout) :: ierr
  character(len=*), intent(in)    :: msg
end subroutine iotk_error_msg_i
end interface

interface iotk_error_check
function iotk_error_check_e(error)
  use iotk_base
  implicit none
  type(iotk_error), intent(in) :: error
  logical :: iotk_error_check_e
end function iotk_error_check_e
function iotk_error_check_i(ierr)
  use iotk_base
  implicit none
  integer, intent(in) :: ierr
  logical :: iotk_error_check_i
end function iotk_error_check_i
end interface

interface iotk_error_write
subroutine iotk_error_write_character_e(error,name,val)
  use iotk_base
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: name
  character(len=*), intent(in)    :: val
end subroutine iotk_error_write_character_e
subroutine iotk_error_write_character_i(ierr,name,val)
  use iotk_base
  implicit none
  integer, intent(inout) :: ierr
  character(len=*), intent(in)    :: name
  character(len=*), intent(in)    :: val
end subroutine iotk_error_write_character_i
subroutine iotk_error_write_logical_e(error,name,val)
  use iotk_base
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: name
  logical,          intent(in)    :: val
end subroutine iotk_error_write_logical_e
subroutine iotk_error_write_logical_i(ierr,name,val)
  use iotk_base
  implicit none
  integer, intent(inout) :: ierr
  character(len=*), intent(in)    :: name
  logical,          intent(in)    :: val
end subroutine iotk_error_write_logical_i 
subroutine iotk_error_write_integer_e(error,name,val)
  use iotk_base
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: name
  integer,          intent(in)    :: val
end subroutine iotk_error_write_integer_e
subroutine iotk_error_write_integer_i(ierr,name,val)
  use iotk_base
  implicit none
  integer, intent(inout) :: ierr
  character(len=*), intent(in)    :: name
  integer,          intent(in)    :: val
end subroutine iotk_error_write_integer_i
end interface

interface iotk_error_scan
subroutine iotk_error_scan_character_e(error,name,val)
  use iotk_base
  implicit none
  type(iotk_error), intent(in) :: error
  character(len=*), intent(in) :: name
#ifdef __IOTK_WORKAROUND6
  character(len=*)             :: val
#else
  character(len=*), intent(out):: val
#endif
end subroutine iotk_error_scan_character_e
subroutine iotk_error_scan_character_i(ierr,name,val)
  use iotk_base
  implicit none
  integer, intent(in) :: ierr
  character(len=*), intent(in) :: name
#ifdef __IOTK_WORKAROUND6
  character(len=*)             :: val
#else
  character(len=*), intent(out):: val
#endif
end subroutine iotk_error_scan_character_i
subroutine iotk_error_scan_logical_e(error,name,val)
  use iotk_base
  implicit none
  type(iotk_error), intent(in) :: error
  character(len=*), intent(in) :: name
  logical,          intent(out):: val
end subroutine iotk_error_scan_logical_e
subroutine iotk_error_scan_logical_i(ierr,name,val)
  use iotk_base
  implicit none
  integer, intent(in) :: ierr
  character(len=*), intent(in) :: name
  logical,          intent(out):: val
end subroutine iotk_error_scan_logical_i
subroutine iotk_error_scan_integer_e(error,name,val)
  use iotk_base
  implicit none
  type(iotk_error), intent(in) :: error
  character(len=*), intent(in) :: name
  integer,          intent(out):: val
end subroutine iotk_error_scan_integer_e
subroutine iotk_error_scan_integer_i(ierr,name,val)
  use iotk_base
  implicit none
  integer, intent(in) :: ierr
  character(len=*), intent(in) :: name
  integer,          intent(out):: val
end subroutine iotk_error_scan_integer_i
end interface

interface iotk_error_pool_pending
function iotk_error_pool_pending_x()
  use iotk_base
  implicit none
  integer :: iotk_error_pool_pending_x
end function iotk_error_pool_pending_x
end interface

interface iotk_error_handler
subroutine iotk_error_handler_x(ierr)
  use iotk_base
  implicit none
  integer, intent(in) :: ierr
end subroutine iotk_error_handler_x
end interface

end module iotk_error_interf
