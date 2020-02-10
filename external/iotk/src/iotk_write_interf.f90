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

module iotk_write_interf
implicit none
private

public :: iotk_write_begin
public :: iotk_write_end
public :: iotk_write_pi
public :: iotk_write_comment
public :: iotk_write_empty
public :: iotk_write_tag

interface iotk_write_begin
!-<
subroutine iotk_write_begin_x(unit,name,attr,dummy,new_line,ierr)
!->
  use iotk_base
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: name
  character(len=*), optional, intent(in)  :: attr
  type(iotk_dummytype), optional          :: dummy
  integer,      optional, intent(out) :: ierr
!-<
  logical, optional, intent(in)           :: new_line
!->
end subroutine iotk_write_begin_x
end interface

interface iotk_write_end
!-<
subroutine iotk_write_end_x(unit,name,dummy,indentation,ierr)
!->
  use iotk_base
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: name
  type(iotk_dummytype), optional          :: dummy
!-<
  logical, optional,          intent(in)  :: indentation
!->
  integer,          optional, intent(out) :: ierr
end subroutine iotk_write_end_x
end interface

interface iotk_write_pi
subroutine iotk_write_pi_x(unit,name,attr,dummy,ierr)
  use iotk_base
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: name
  character(len=*), optional, intent(in)  :: attr
  type(iotk_dummytype), optional          :: dummy
  integer,          optional, intent(out) :: ierr
end subroutine iotk_write_pi_x
end interface

interface iotk_write_comment
subroutine iotk_write_comment_x(unit,text,dummy,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: text
  type(iotk_dummytype), optional      :: dummy
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_comment_x
end interface

interface iotk_write_empty
subroutine iotk_write_empty_x(unit,name,attr,dummy,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(in)  :: attr
  type(iotk_dummytype), optional      :: dummy
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_empty_x
end interface

interface iotk_write_tag
!-<
subroutine iotk_write_tag_x(unit,control,tag,binary,indent,ierr,new_line)
!->
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  integer,                   intent(in)  :: control
  character(iotk_taglenx),   intent(in)  :: tag
  logical,                   intent(in)  :: binary
  integer,                   intent(in)  :: indent
  integer,                   intent(out) :: ierr
!-<
  logical,optional,          intent(in)  :: new_line
!->
end subroutine iotk_write_tag_x
end interface

end module iotk_write_interf
