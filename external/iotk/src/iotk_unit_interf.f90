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

module iotk_unit_interf
use iotk_base
implicit none
private

public :: iotk_free_unit
public :: iotk_phys_unit
public :: iotk_unit_print
public :: iotk_unit_add
public :: iotk_inquire
public :: iotk_unit_del
public :: iotk_unit_parent
public :: iotk_unit_get
!-<
public :: iotk_rewind
!->

! This module contains the interfaces to all iotk routines

interface iotk_free_unit
subroutine iotk_free_unit_x(unit,ierr)
  implicit none
  integer,           intent(out) :: unit
  integer, optional, intent(out) :: ierr  
end subroutine iotk_free_unit_x
end interface

interface iotk_phys_unit
function iotk_phys_unit_x(unit)
  implicit none
  integer, intent(in) :: unit
  integer :: iotk_phys_unit_x
end function iotk_phys_unit_x
end interface

interface iotk_unit_print
subroutine iotk_unit_print_x(unit)
  implicit none
  integer, intent(in) :: unit
end subroutine iotk_unit_print_x
end interface

interface iotk_unit_add
subroutine iotk_unit_add_x(unit,this,ierr)
  use iotk_base
  implicit none
  integer,           intent(in)  :: unit
  type(iotk_unit),   pointer     :: this
  integer,           intent(out) :: ierr
end subroutine iotk_unit_add_x
end interface

interface iotk_inquire
subroutine iotk_inquire_x(unit,binary,stream,ierr)
  implicit none
  integer, intent(in)  :: unit
  logical, intent(out) :: binary
  logical, intent(out) :: stream
  integer, intent(out) :: ierr
end subroutine iotk_inquire_x
end interface

interface iotk_unit_del
subroutine iotk_unit_del_x(unit,ierr)
  implicit none
  integer, intent(in)  :: unit
  integer, intent(out) :: ierr
end subroutine iotk_unit_del_x
end interface

interface iotk_unit_parent
subroutine iotk_unit_parent_x(parent,son,ierr)
  implicit none
  integer, intent(in) :: parent,son
  integer, intent(out) :: ierr
end subroutine iotk_unit_parent_x
end interface

interface iotk_unit_get
subroutine iotk_unit_get_x(unit,pointer)
  use iotk_base
  implicit none
  integer,     intent(in)  :: unit
  type(iotk_unit), pointer :: pointer
end subroutine iotk_unit_get_x
end interface

!-<
interface iotk_rewind
subroutine iotk_rewind_x(unit ,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  integer, optional,               intent(out) :: ierr
end subroutine iotk_rewind_x
end interface
!->

end module iotk_unit_interf
