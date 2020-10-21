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

module iotk_attr_interf
implicit none
private

public :: iotk_read
public :: iotk_write
public :: iotk_write_attr
public :: iotk_scan_attr


interface iotk_read
#ifdef __IOTK_LOGICAL1
subroutine iotk_read_LOGICAL1(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  LOGICAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_LOGICAL1
#endif
#ifdef __IOTK_LOGICAL2
subroutine iotk_read_LOGICAL2(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  LOGICAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_LOGICAL2
#endif
#ifdef __IOTK_LOGICAL3
subroutine iotk_read_LOGICAL3(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  LOGICAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_LOGICAL3
#endif
#ifdef __IOTK_LOGICAL4
subroutine iotk_read_LOGICAL4(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  LOGICAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_LOGICAL4
#endif
#ifdef __IOTK_INTEGER1
subroutine iotk_read_INTEGER1(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  INTEGER(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_INTEGER1
#endif
#ifdef __IOTK_INTEGER2
subroutine iotk_read_INTEGER2(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  INTEGER(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_INTEGER2
#endif
#ifdef __IOTK_INTEGER3
subroutine iotk_read_INTEGER3(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  INTEGER(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_INTEGER3
#endif
#ifdef __IOTK_INTEGER4
subroutine iotk_read_INTEGER4(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  INTEGER(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_INTEGER4
#endif
#ifdef __IOTK_REAL1
subroutine iotk_read_REAL1(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  REAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_REAL1
#endif
#ifdef __IOTK_REAL2
subroutine iotk_read_REAL2(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  REAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_REAL2
#endif
#ifdef __IOTK_REAL3
subroutine iotk_read_REAL3(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  REAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_REAL3
#endif
#ifdef __IOTK_REAL4
subroutine iotk_read_REAL4(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  REAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_REAL4
#endif
#ifdef __IOTK_COMPLEX1
subroutine iotk_read_COMPLEX1(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  COMPLEX(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_COMPLEX1
#endif
#ifdef __IOTK_COMPLEX2
subroutine iotk_read_COMPLEX2(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  COMPLEX(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_COMPLEX2
#endif
#ifdef __IOTK_COMPLEX3
subroutine iotk_read_COMPLEX3(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  COMPLEX(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_COMPLEX3
#endif
#ifdef __IOTK_COMPLEX4
subroutine iotk_read_COMPLEX4(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  COMPLEX(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_COMPLEX4
#endif
end interface

interface iotk_write
#ifdef __IOTK_LOGICAL1
subroutine iotk_write_LOGICAL1(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  LOGICAL(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_LOGICAL1
#endif
#ifdef __IOTK_LOGICAL2
subroutine iotk_write_LOGICAL2(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  LOGICAL(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_LOGICAL2
#endif
#ifdef __IOTK_LOGICAL3
subroutine iotk_write_LOGICAL3(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  LOGICAL(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_LOGICAL3
#endif
#ifdef __IOTK_LOGICAL4
subroutine iotk_write_LOGICAL4(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  LOGICAL(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_LOGICAL4
#endif
#ifdef __IOTK_INTEGER1
subroutine iotk_write_INTEGER1(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  INTEGER(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_INTEGER1
#endif
#ifdef __IOTK_INTEGER2
subroutine iotk_write_INTEGER2(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  INTEGER(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_INTEGER2
#endif
#ifdef __IOTK_INTEGER3
subroutine iotk_write_INTEGER3(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  INTEGER(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_INTEGER3
#endif
#ifdef __IOTK_INTEGER4
subroutine iotk_write_INTEGER4(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  INTEGER(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_INTEGER4
#endif
#ifdef __IOTK_REAL1
subroutine iotk_write_REAL1(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  REAL(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_REAL1
#endif
#ifdef __IOTK_REAL2
subroutine iotk_write_REAL2(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  REAL(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_REAL2
#endif
#ifdef __IOTK_REAL3
subroutine iotk_write_REAL3(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  REAL(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_REAL3
#endif
#ifdef __IOTK_REAL4
subroutine iotk_write_REAL4(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  REAL(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_REAL4
#endif
#ifdef __IOTK_COMPLEX1
subroutine iotk_write_COMPLEX1(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  COMPLEX(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_COMPLEX1
#endif
#ifdef __IOTK_COMPLEX2
subroutine iotk_write_COMPLEX2(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  COMPLEX(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_COMPLEX2
#endif
#ifdef __IOTK_COMPLEX3
subroutine iotk_write_COMPLEX3(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  COMPLEX(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_COMPLEX3
#endif
#ifdef __IOTK_COMPLEX4
subroutine iotk_write_COMPLEX4(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  COMPLEX(kind=this_kind), intent(in)  :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)                     :: string
#else
  character(len=*),        intent(out) :: string
#endif
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_COMPLEX4
#endif
end interface

interface iotk_write_attr
#ifdef __IOTK_LOGICAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_7
#endif
#endif
#ifdef __IOTK_LOGICAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_7
#endif
#endif
#ifdef __IOTK_LOGICAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_7
#endif
#endif
#ifdef __IOTK_LOGICAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_7
#endif
#endif
#ifdef __IOTK_INTEGER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_7
#endif
#endif
#ifdef __IOTK_INTEGER2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_7
#endif
#endif
#ifdef __IOTK_INTEGER3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_7
#endif
#endif
#ifdef __IOTK_INTEGER4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_7
#endif
#endif
#ifdef __IOTK_REAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_7
#endif
#endif
#ifdef __IOTK_REAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_7
#endif
#endif
#ifdef __IOTK_REAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_7
#endif
#endif
#ifdef __IOTK_REAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_7
#endif
#endif
#ifdef __IOTK_COMPLEX1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_7
#endif
#endif
#ifdef __IOTK_COMPLEX2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_7
#endif
#endif
#ifdef __IOTK_COMPLEX3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_7
#endif
#endif
#ifdef __IOTK_COMPLEX4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_5(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_6(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_7(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_7
#endif
#endif
#ifdef __IOTK_CHARACTER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_CHARACTER1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_CHARACTER1_0
#endif
#endif
end interface

interface iotk_scan_attr
#ifdef __IOTK_LOGICAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val 
#else
  LOGICAL(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_7
#endif
#endif
#ifdef __IOTK_LOGICAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val 
#else
  LOGICAL(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_7
#endif
#endif
#ifdef __IOTK_LOGICAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val 
#else
  LOGICAL(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_7
#endif
#endif
#ifdef __IOTK_LOGICAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val 
#else
  LOGICAL(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_7
#endif
#endif
#ifdef __IOTK_INTEGER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val 
#else
  INTEGER(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_7
#endif
#endif
#ifdef __IOTK_INTEGER2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val 
#else
  INTEGER(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_7
#endif
#endif
#ifdef __IOTK_INTEGER3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val 
#else
  INTEGER(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_7
#endif
#endif
#ifdef __IOTK_INTEGER4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val 
#else
  INTEGER(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_7
#endif
#endif
#ifdef __IOTK_REAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val 
#else
  REAL(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_7
#endif
#endif
#ifdef __IOTK_REAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val 
#else
  REAL(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_7
#endif
#endif
#ifdef __IOTK_REAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val 
#else
  REAL(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_7
#endif
#endif
#ifdef __IOTK_REAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val 
#else
  REAL(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_7
#endif
#endif
#ifdef __IOTK_COMPLEX1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val 
#else
  COMPLEX(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_7
#endif
#endif
#ifdef __IOTK_COMPLEX2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val 
#else
  COMPLEX(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_7
#endif
#endif
#ifdef __IOTK_COMPLEX3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val 
#else
  COMPLEX(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_7
#endif
#endif
#ifdef __IOTK_COMPLEX4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val 
#else
  COMPLEX(kind=this_kind),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_5(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_6(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_7(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: val (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_7
#endif
#endif
#ifdef __IOTK_CHARACTER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_CHARACTER1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: val 
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: val 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_CHARACTER1_0
#endif
#endif
end interface

end module iotk_attr_interf
