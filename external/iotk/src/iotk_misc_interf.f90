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

module iotk_misc_interf
implicit none
private

public :: iotk_copy_tag
public :: iotk_parse_dat
public :: iotk_set
public :: iotk_get
public :: iotk_copy_dat_aux
public :: iotk_copy_dat
public :: iotk_print_kinds
public :: iotk_check_iotk_attr
public :: iotk_index
public :: iotk_check_name
public :: iotk_tag_parse
public :: iotk_complete_filepath
public :: iotk_delete_attr
public :: iotk_readcmdline
public :: iotk_init_static_vars

! This module contains the interfaces to all iotk routines

interface iotk_copy_tag
subroutine iotk_copy_tag_x(source,dest,dummy,maxsize,ierr)
  use iotk_base
  implicit none
  integer,           intent(in)  :: source
  integer,           intent(in)  :: dest
  type(iotk_dummytype), optional :: dummy
  integer, optional, intent(in)  :: maxsize
  integer, optional, intent(out) :: ierr
end subroutine iotk_copy_tag_x
end interface

interface iotk_parse_dat
subroutine iotk_parse_dat_x(attr,type,kind,isize,len,fmt,columns,ierr)
  implicit none
  character(len=*), intent(in)  :: attr
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: type
#else
  character(len=*), intent(out) :: type
#endif
  integer,          intent(out) :: kind
  integer,          intent(out) :: isize
  integer,          intent(out) :: len
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: fmt
#else
  character(len=*), intent(out) :: fmt
#endif
  integer,          intent(out) :: columns
  integer,          intent(out) :: ierr
end subroutine iotk_parse_dat_x
end interface

interface iotk_set
subroutine iotk_set_x(dummy,unitmin,unitmax,getline_buffer,error_warn_overflow, &
                              linlen,indent,maxindent,error_unit,output_unit,ierr)
  use iotk_base
  implicit none
  type(iotk_dummytype),optional  :: dummy
  integer, optional, intent(in)  :: unitmin
  integer, optional, intent(in)  :: unitmax
  integer, optional, intent(in)  :: getline_buffer
  logical, optional, intent(in)  :: error_warn_overflow
  integer, optional, intent(in)  :: linlen
  integer, optional, intent(in)  :: indent
  integer, optional, intent(in)  :: maxindent
  integer, optional, intent(in)  :: error_unit
  integer, optional, intent(in)  :: output_unit
  integer, optional, intent(out) :: ierr
end subroutine iotk_set_x
end interface

interface iotk_get
subroutine iotk_get_x(dummy,unitmin,unitmax,getline_buffer,error_warn_overflow, &
                              linlen,indent,maxindent,error_unit,output_unit)
  use iotk_base
  implicit none
  type(iotk_dummytype),optional  :: dummy
  integer, optional, intent(out) :: unitmin
  integer, optional, intent(out) :: unitmax
  integer, optional, intent(out) :: getline_buffer
  logical, optional, intent(out) :: error_warn_overflow
  integer, optional, intent(out) :: linlen
  integer, optional, intent(out) :: indent
  integer, optional, intent(out) :: maxindent
  integer, optional, intent(out) :: error_unit
  integer, optional, intent(out) :: output_unit
end subroutine iotk_get_x
end interface

interface iotk_copy_dat_aux
subroutine iotk_copy_dat_aux_x(source,dest,source_binary,dest_binary,name,type,ikind,isize, &
                               len,fmt,columns,attr,ierr)
  implicit none
  integer,      intent(in)  :: source
  integer,      intent(in)  :: dest
  logical,      intent(in)  :: source_binary
  logical,      intent(in)  :: dest_binary
  character(*), intent(in)  :: name
  character(*), intent(in)  :: type
  integer,      intent(in)  :: ikind
  integer,      intent(in)  :: isize
  integer,      intent(in)  :: len
  character(*), intent(in)  :: fmt
  integer,      intent(in)  :: columns
  character(*), intent(in)  :: attr
  integer,      intent(out) :: ierr
end subroutine iotk_copy_dat_aux_x
end interface

interface iotk_copy_dat
subroutine iotk_copy_dat_x(source,dest,source_binary,dest_binary,name,attr,maxsize,ierr)
  implicit none
  integer,      intent(in)  :: source
  integer,      intent(in)  :: dest
  logical,      intent(in)  :: source_binary
  logical,      intent(in)  :: dest_binary
  character(*), intent(in)  :: name
  character(*), intent(in)  :: attr
  integer,      intent(in)  :: maxsize
  integer,      intent(out) :: ierr
end subroutine iotk_copy_dat_x
end interface

interface iotk_print_kinds
subroutine iotk_print_kinds_x
end subroutine iotk_print_kinds_x
end interface

interface iotk_check_iotk_attr
subroutine iotk_check_iotk_attr_x(unit,attr,ierr)
  use iotk_base
  implicit none
  integer,                 intent(in)  :: unit
  character(iotk_attlenx), intent(in)  :: attr
  integer,                 intent(out) :: ierr
end subroutine iotk_check_iotk_attr_x
end interface


interface iotk_index
function iotk_index_scal(index)
  implicit none
  integer,           intent(in) :: index
  character(len=range(index)+3) :: iotk_index_scal
end function iotk_index_scal
function iotk_index_vec(index)
  implicit none
  integer,                         intent(in) :: index(:)
  character(len=(range(index)+3)*size(index)) :: iotk_index_vec
end function iotk_index_vec
end interface

interface iotk_tag_parse
subroutine iotk_tag_parse_x(tag,name,attr,ierr)
  use iotk_base
  implicit none
  character(iotk_taglenx), intent(in)  :: tag
  character(iotk_namlenx), intent(out) :: name
  character(iotk_attlenx), intent(out) :: attr
  integer,                 intent(out) :: ierr
end subroutine iotk_tag_parse_x
end interface

interface iotk_complete_filepath
function iotk_complete_filepath_x(newfile,oldfile)
  implicit none
  character(len=*), intent(in) :: newfile
  character(len=*), intent(in) :: oldfile
  character(len=len(newfile)+len(oldfile)) :: iotk_complete_filepath_x
end function iotk_complete_filepath_x
end interface

interface iotk_check_name
function iotk_check_name_x(name)
  implicit none
  character(len=*), intent(in) :: name
  logical                      :: iotk_check_name_x
end function iotk_check_name_x
end interface

interface iotk_delete_attr
subroutine iotk_delete_attr_x(attr,name,ierr)
  implicit none
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  integer,          intent(out)   :: ierr
end subroutine iotk_delete_attr_x
end interface

interface iotk_readcmdline
subroutine iotk_readcmdline_x(args,nargs,eos,ierr)
  implicit none
  character(len=*),  intent(out) :: args(:)
  integer,           intent(out) :: nargs
  logical, optional, intent(in)  :: eos
  integer, optional, intent(out) :: ierr
end subroutine iotk_readcmdline_x
end interface

interface iotk_init_static_vars
subroutine iotk_init_static_vars_x()
end subroutine iotk_init_static_vars_x
end interface

end module iotk_misc_interf
