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

module iotk_scan_interf
implicit none
private

public :: iotk_scan_begin
public :: iotk_scan_end
public :: iotk_scan_pi
public :: iotk_scan_empty
public :: iotk_scan_tag
public :: iotk_scan
public :: iotk_getline

interface iotk_scan_begin
!-<
subroutine iotk_scan_begin_x(unit,name,attr,dummy,found,direction,ierr)
!->
  use iotk_base
  implicit none
  integer,                        intent(in)  :: unit
  character(len=*),               intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  character(len=*),     optional              :: attr
#else
  character(len=*),     optional, intent(out) :: attr
#endif
  type(iotk_dummytype), optional              :: dummy
  logical,              optional, intent(out) :: found
!-<
  integer,              optional, intent(out) :: direction
!->
  integer,              optional, intent(out) :: ierr
end subroutine iotk_scan_begin_x
end interface

interface iotk_scan_end
subroutine iotk_scan_end_x(unit,name,dummy,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  type(iotk_dummytype), optional      :: dummy
  integer,      optional, intent(out) :: ierr
end subroutine iotk_scan_end_x
end interface

interface iotk_scan_pi
subroutine iotk_scan_pi_x(unit,name,attr,dummy,found,ierr)
  use iotk_base
  implicit none
  integer,                        intent(in)  :: unit
  character(*),                   intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  character(len=*),     optional              :: attr
#else
  character(len=*),     optional, intent(out) :: attr
#endif
  type(iotk_dummytype), optional              :: dummy
  logical,              optional, intent(out) :: found
  integer,              optional, intent(out) :: ierr
end subroutine iotk_scan_pi_x
end interface

interface iotk_scan_empty
subroutine iotk_scan_empty_x(unit,name,attr,dummy,found,ierr)
  use iotk_base
  implicit none
  integer,                        intent(in)  :: unit
  character(len=*),               intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  character(len=*),     optional              :: attr
#else
  character(len=*),     optional, intent(out) :: attr
#endif
  type(iotk_dummytype), optional              :: dummy
  logical,              optional, intent(out) :: found
  integer,              optional, intent(out) :: ierr
end subroutine iotk_scan_empty_x
end interface

interface iotk_scan_tag
subroutine iotk_scan_tag_x(unit,direction,control,tag,binary,stream,ierr)
  use iotk_base
  implicit none
  integer,                 intent(in)  :: unit
  integer,                 intent(in)  :: direction
  integer,                 intent(out) :: control
  character(iotk_taglenx), intent(out) :: tag
  logical,                 intent(in)  :: binary
  logical,                 intent(in)  :: stream
  integer,                 intent(out) :: ierr
end subroutine iotk_scan_tag_x
end interface

interface iotk_scan
subroutine iotk_scan_x(unit,direction,control,name,attr,binary,stream,found,ierr)
  use iotk_base
  implicit none
  integer,                 intent(in)  :: unit
  integer,                 intent(in)  :: direction
  integer,                 intent(in)  :: control
  character(iotk_namlenx), intent(in)  :: name
  character(iotk_attlenx), intent(out) :: attr
  logical,                 intent(in)  :: binary
  logical,                 intent(in)  :: stream
  logical,                 intent(out) :: found
  integer,                 intent(out) :: ierr
end subroutine iotk_scan_x
end interface

interface iotk_getline
subroutine iotk_getline_x(unit,line,length,ierr)
  implicit none
  integer,           intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  character(len=*)               :: line
#else
  character(len=*),  intent(out) :: line
#endif
  integer, optional, intent(out) :: length
  integer, optional, intent(out) :: ierr
end subroutine iotk_getline_x
end interface

end module iotk_scan_interf
