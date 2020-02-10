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

module iotk_files_interf
use iotk_base
implicit none
private

public :: iotk_copyfile
public :: iotk_link
public :: iotk_open_write
public :: iotk_close_write
public :: iotk_open_read
public :: iotk_close_read
public :: iotk_magic

interface iotk_copyfile
subroutine iotk_copyfile_x(dummy,source,dest,source_unit,dest_unit,ierr)
  use iotk_base
  implicit none
  type(iotk_dummytype), optional         :: dummy
  character(len=*), optional, intent(in) :: source
  character(len=*), optional, intent(in) :: dest
  integer,          optional, intent(in) :: source_unit
  integer,          optional, intent(in) :: dest_unit
  integer,          optional, intent(out):: ierr
end subroutine iotk_copyfile_x
end interface

interface iotk_link
subroutine iotk_link_x(unit,name,file,dummy,binary,raw,create,ierr)
  use iotk_base
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: name
  character(len=*),           intent(in)  :: file
  type(iotk_dummytype), optional          :: dummy
  logical,          optional, intent(in)  :: binary
  logical,          optional, intent(in)  :: raw
  logical,          optional, intent(in)  :: create
  integer,          optional, intent(out) :: ierr
end subroutine iotk_link_x
end interface

interface iotk_open_write
!-<
subroutine iotk_open_write_x(unit,file,dummy,attr,binary,new,raw,root,qe_syntax,skip_root,skip_head,ierr)
!->
  use iotk_base
  implicit none
  integer,                    intent(in)  :: unit
  type(iotk_dummytype), optional          :: dummy
  character(len=*), optional, intent(in)  :: file
  character(len=*), optional, intent(in)  :: attr
  logical,          optional, intent(in)  :: binary
  logical,          optional, intent(in)  :: new
  logical,          optional, intent(in)  :: raw
  character(len=*), optional, intent(in)  :: root
!-<
  logical,      optional, intent(in)  :: qe_syntax
!->
  logical,          optional, intent(in)  :: skip_root
  logical,          optional, intent(in)  :: skip_head
  integer,          optional, intent(out) :: ierr
end subroutine iotk_open_write_x
end interface

interface iotk_close_write
subroutine iotk_close_write_x(unit,dummy,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  type(iotk_dummytype), optional      :: dummy
  integer,      optional, intent(out) :: ierr
end subroutine iotk_close_write_x
end interface

interface iotk_open_read
!-<
subroutine iotk_open_read_x(unit,file,dummy,attr,binary,stream,raw,qe_syntax,root,ierr)
!->
  use iotk_base
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*), optional, intent(in)  :: file
  type(iotk_dummytype), optional          :: dummy
  logical,          optional, intent(in)  :: binary 
  logical,          optional, intent(in)  :: stream
  logical,          optional, intent(in)  :: raw
!-<
  logical,      optional, intent(in)      :: qe_syntax
!->
#ifdef __IOTK_WORKAROUND6
  character(len=*), optional              :: attr
  character(len=*), optional              :: root
#else
  character(len=*), optional, intent(out) :: attr
  character(len=*), optional, intent(out) :: root
#endif
  integer,          optional, intent(out) :: ierr
end subroutine iotk_open_read_x
end interface

interface iotk_close_read
subroutine iotk_close_read_x(unit,dummy,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  type(iotk_dummytype), optional      :: dummy
  integer,      optional, intent(out) :: ierr
end subroutine iotk_close_read_x
end interface

interface iotk_magic
subroutine iotk_magic_x(file,binary)
  implicit none
  character(len=*), intent(in)  :: file
  logical,          intent(out) :: binary
end subroutine iotk_magic_x
end interface

end module iotk_files_interf
