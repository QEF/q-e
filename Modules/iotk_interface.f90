! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_attr_interf.spp"

module iotk_attr_interf
implicit none
private

public :: iotk_read
public :: iotk_write
public :: iotk_write_attr
public :: iotk_scan_attr


interface iotk_read
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL1
subroutine iotk_read_LOGICAL1(val,string,index,ierr)
  use iotk_base
  implicit none
  LOGICAL(kind=__IOTK_LOGICAL1), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_LOGICAL1
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL2
subroutine iotk_read_LOGICAL2(val,string,index,ierr)
  use iotk_base
  implicit none
  LOGICAL(kind=__IOTK_LOGICAL2), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_LOGICAL2
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL3
subroutine iotk_read_LOGICAL3(val,string,index,ierr)
  use iotk_base
  implicit none
  LOGICAL(kind=__IOTK_LOGICAL3), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_LOGICAL3
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL4
subroutine iotk_read_LOGICAL4(val,string,index,ierr)
  use iotk_base
  implicit none
  LOGICAL(kind=__IOTK_LOGICAL4), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_LOGICAL4
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER1
subroutine iotk_read_INTEGER1(val,string,index,ierr)
  use iotk_base
  implicit none
  INTEGER(kind=__IOTK_INTEGER1), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_INTEGER1
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER2
subroutine iotk_read_INTEGER2(val,string,index,ierr)
  use iotk_base
  implicit none
  INTEGER(kind=__IOTK_INTEGER2), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_INTEGER2
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER3
subroutine iotk_read_INTEGER3(val,string,index,ierr)
  use iotk_base
  implicit none
  INTEGER(kind=__IOTK_INTEGER3), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_INTEGER3
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER4
subroutine iotk_read_INTEGER4(val,string,index,ierr)
  use iotk_base
  implicit none
  INTEGER(kind=__IOTK_INTEGER4), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_INTEGER4
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL1
subroutine iotk_read_REAL1(val,string,index,ierr)
  use iotk_base
  implicit none
  REAL(kind=__IOTK_REAL1), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_REAL1
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL2
subroutine iotk_read_REAL2(val,string,index,ierr)
  use iotk_base
  implicit none
  REAL(kind=__IOTK_REAL2), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_REAL2
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL3
subroutine iotk_read_REAL3(val,string,index,ierr)
  use iotk_base
  implicit none
  REAL(kind=__IOTK_REAL3), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_REAL3
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL4
subroutine iotk_read_REAL4(val,string,index,ierr)
  use iotk_base
  implicit none
  REAL(kind=__IOTK_REAL4), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_REAL4
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX1
subroutine iotk_read_COMPLEX1(val,string,index,ierr)
  use iotk_base
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX1), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_COMPLEX1
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX2
subroutine iotk_read_COMPLEX2(val,string,index,ierr)
  use iotk_base
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX2), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_COMPLEX2
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX3
subroutine iotk_read_COMPLEX3(val,string,index,ierr)
  use iotk_base
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX3), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_COMPLEX3
#endif
# 45 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX4
subroutine iotk_read_COMPLEX4(val,string,index,ierr)
  use iotk_base
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX4), intent(inout) :: val(:)
  character(len=*),                    intent(in)    :: string
  integer,                             intent(inout) :: index
  integer,                             intent(out) :: ierr
end subroutine iotk_read_COMPLEX4
#endif
# 58 "iotk_attr_interf.spp"
end interface

interface iotk_write
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL1
subroutine iotk_write_LOGICAL1(val,string,ierr)
  use iotk_base
  implicit none
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_LOGICAL1
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL2
subroutine iotk_write_LOGICAL2(val,string,ierr)
  use iotk_base
  implicit none
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_LOGICAL2
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL3
subroutine iotk_write_LOGICAL3(val,string,ierr)
  use iotk_base
  implicit none
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_LOGICAL3
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL4
subroutine iotk_write_LOGICAL4(val,string,ierr)
  use iotk_base
  implicit none
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_LOGICAL4
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER1
subroutine iotk_write_INTEGER1(val,string,ierr)
  use iotk_base
  implicit none
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_INTEGER1
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER2
subroutine iotk_write_INTEGER2(val,string,ierr)
  use iotk_base
  implicit none
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_INTEGER2
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER3
subroutine iotk_write_INTEGER3(val,string,ierr)
  use iotk_base
  implicit none
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_INTEGER3
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER4
subroutine iotk_write_INTEGER4(val,string,ierr)
  use iotk_base
  implicit none
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_INTEGER4
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL1
subroutine iotk_write_REAL1(val,string,ierr)
  use iotk_base
  implicit none
  REAL(kind=__IOTK_REAL1), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_REAL1
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL2
subroutine iotk_write_REAL2(val,string,ierr)
  use iotk_base
  implicit none
  REAL(kind=__IOTK_REAL2), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_REAL2
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL3
subroutine iotk_write_REAL3(val,string,ierr)
  use iotk_base
  implicit none
  REAL(kind=__IOTK_REAL3), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_REAL3
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL4
subroutine iotk_write_REAL4(val,string,ierr)
  use iotk_base
  implicit none
  REAL(kind=__IOTK_REAL4), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_REAL4
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX1
subroutine iotk_write_COMPLEX1(val,string,ierr)
  use iotk_base
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_COMPLEX1
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX2
subroutine iotk_write_COMPLEX2(val,string,ierr)
  use iotk_base
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_COMPLEX2
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX3
subroutine iotk_write_COMPLEX3(val,string,ierr)
  use iotk_base
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_COMPLEX3
#endif
# 64 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX4
subroutine iotk_write_COMPLEX4(val,string,ierr)
  use iotk_base
  implicit none
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val(:)
  character(len=*),                    intent(out) :: string
  integer,                             intent(out) :: ierr
end subroutine iotk_write_COMPLEX4
#endif
# 76 "iotk_attr_interf.spp"
end interface

interface iotk_write_attr
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL1
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL1_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL2
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL2_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL2_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL3
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL3_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL3_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL4
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_LOGICAL4_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL4_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER1
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER1_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER2
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER2_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER2_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER3
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER3_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER3_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER4
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_INTEGER4_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER4_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL1
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL1_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL2
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL2_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL3
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL3_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL3_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL4
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_REAL4_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL4_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX1
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX1_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX2
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX2_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX3
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX3_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX3_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX4
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_0
#endif
# 85 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_1(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_1
#endif
# 85 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_2(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_2
#endif
# 85 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_3(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_3
#endif
# 85 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_4(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_4
#endif
# 85 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_5(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_5
#endif
# 85 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_6(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_6
#endif
# 85 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_attr_COMPLEX4_7(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4), intent(in)  :: val (:,:,:,:,:,:,:)
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX4_7
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 82 "iotk_attr_interf.spp"
#ifdef __IOTK_CHARACTER1
# 85 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_attr_CHARACTER1_0(attr,name,val,first,ierr)
  use iotk_base
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
# 92 "iotk_attr_interf.spp"
  CHARACTER(kind=__IOTK_CHARACTER1,len=*), intent(in)  :: val 
  logical, optional, intent(in)  :: first
  integer, optional, intent(out) :: ierr
end subroutine iotk_write_attr_CHARACTER1_0
#endif
# 99 "iotk_attr_interf.spp"
#endif
# 103 "iotk_attr_interf.spp"
end interface

interface iotk_scan_attr
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL1
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL1), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL1_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL1),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL2
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL2), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL2_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL2),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL2_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL3
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL3), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL3_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL3),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL3_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_LOGICAL4
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL4), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_LOGICAL4_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  LOGICAL(kind=__IOTK_LOGICAL4),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  LOGICAL(kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL4_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER1
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER1_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER1),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER2
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER2_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER2),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER2_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER3
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER3_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER3),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER3_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_INTEGER4
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_INTEGER4_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  INTEGER(kind=__IOTK_INTEGER4),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  INTEGER(kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER4_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL1
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL1), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL1), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL1_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL1),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL2
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL2), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL2), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL2_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL2),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL3
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL3), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL3), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL3_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL3),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL3_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_REAL4
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL4), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL4), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_REAL4_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  REAL(kind=__IOTK_REAL4),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  REAL(kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL4_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX1
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX1_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX1),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX2
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX2_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX2),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX3
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX3_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX3),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX3_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_COMPLEX4
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_0
#endif
# 112 "iotk_attr_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_1(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_1
#endif
# 112 "iotk_attr_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_2(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_2
#endif
# 112 "iotk_attr_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_3(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_3
#endif
# 112 "iotk_attr_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_4(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_4
#endif
# 112 "iotk_attr_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_5(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_5
#endif
# 112 "iotk_attr_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_6(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_6
#endif
# 112 "iotk_attr_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_COMPLEX4_7(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  COMPLEX(kind=__IOTK_COMPLEX4),           intent(out) :: val (:,:,:,:,:,:,:)
  logical,        optional, intent(out) :: found
  COMPLEX(kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX4_7
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 109 "iotk_attr_interf.spp"
#ifdef __IOTK_CHARACTER1
# 112 "iotk_attr_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_attr_CHARACTER1_0(attr,name,val,found,default,eos,ierr)
  use iotk_base
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
# 119 "iotk_attr_interf.spp"
  CHARACTER(kind=__IOTK_CHARACTER1,len=*),           intent(out) :: val 
  logical,        optional, intent(out) :: found
  CHARACTER(kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
end subroutine iotk_scan_attr_CHARACTER1_0
#endif
# 128 "iotk_attr_interf.spp"
#endif
# 132 "iotk_attr_interf.spp"
end interface

end module iotk_attr_interf
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_dat_interf.spp"

module iotk_dat_interf
implicit none
private

public :: iotk_write_dat
public :: iotk_scan_dat
public :: iotk_scan_dat_aux


interface iotk_write_dat
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL1
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL2
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL3
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL4
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER1
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER2
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER3
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER4
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL1
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL2
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL3
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL4
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX1
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX2
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX3
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX4
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 44 "iotk_dat_interf.spp"
#ifdef __IOTK_CHARACTER1
# 46 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_0(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),        intent(in)  :: dat 
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_0
#endif
# 46 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_1(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),        intent(in)  :: dat (:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_1
#endif
# 46 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_2(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),        intent(in)  :: dat (:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_2
#endif
# 46 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_3(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),        intent(in)  :: dat (:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_3
#endif
# 46 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_4(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),        intent(in)  :: dat (:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_4
#endif
# 46 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_5(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),        intent(in)  :: dat (:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_5
#endif
# 46 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_6(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),        intent(in)  :: dat (:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_6
#endif
# 46 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_7(unit,name,dat,fmt,ierr)
  use iotk_base
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
# 53 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),        intent(in)  :: dat (:,:,:,:,:,:,:)
  character(*), optional, intent(in)  :: fmt
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_7
#endif
# 59 "iotk_dat_interf.spp"
#endif
# 63 "iotk_dat_interf.spp"
end interface

interface iotk_scan_dat
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL1
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL1), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL2
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL2), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL3
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL3), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL4
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL4), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  LOGICAL (kind=__IOTK_LOGICAL4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER1
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER2
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER3
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER4
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  INTEGER (kind=__IOTK_INTEGER4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL1
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL1), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL1), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL2
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL2), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL2), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL3
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL3), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL3), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL4
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL4), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL4), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  REAL (kind=__IOTK_REAL4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX1
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX1), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX2
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX2), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX3
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX3), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX4
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  COMPLEX (kind=__IOTK_COMPLEX4), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 69 "iotk_dat_interf.spp"
#ifdef __IOTK_CHARACTER1
# 71 "iotk_dat_interf.spp"
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_0(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat 
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default 
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_0
#endif
# 71 "iotk_dat_interf.spp"
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_1(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_1
#endif
# 71 "iotk_dat_interf.spp"
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_2(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_2
#endif
# 71 "iotk_dat_interf.spp"
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_3(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_3
#endif
# 71 "iotk_dat_interf.spp"
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_4(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_4
#endif
# 71 "iotk_dat_interf.spp"
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_5(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_5
#endif
# 71 "iotk_dat_interf.spp"
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_6(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_6
#endif
# 71 "iotk_dat_interf.spp"
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_7(unit,name,dat,found,default,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  character(*),              intent(in)  :: name
# 78 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*),           intent(out) :: dat (:,:,:,:,:,:,:)
  logical,         optional, intent(out) :: found
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_7
#endif
# 85 "iotk_dat_interf.spp"
#endif
# 89 "iotk_dat_interf.spp"
end interface

interface iotk_scan_dat_aux
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL1
subroutine iotk_scan_dat_aux_LOGICAL1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL1), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_LOGICAL1
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL2
subroutine iotk_scan_dat_aux_LOGICAL2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL2), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_LOGICAL2
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL3
subroutine iotk_scan_dat_aux_LOGICAL3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL3), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_LOGICAL3
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_LOGICAL4
subroutine iotk_scan_dat_aux_LOGICAL4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  LOGICAL (kind=__IOTK_LOGICAL4), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_LOGICAL4
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER1
subroutine iotk_scan_dat_aux_INTEGER1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER1), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_INTEGER1
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER2
subroutine iotk_scan_dat_aux_INTEGER2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER2), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_INTEGER2
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER3
subroutine iotk_scan_dat_aux_INTEGER3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER3), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_INTEGER3
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_INTEGER4
subroutine iotk_scan_dat_aux_INTEGER4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  INTEGER (kind=__IOTK_INTEGER4), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_INTEGER4
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL1
subroutine iotk_scan_dat_aux_REAL1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL1), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_REAL1
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL2
subroutine iotk_scan_dat_aux_REAL2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL2), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_REAL2
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL3
subroutine iotk_scan_dat_aux_REAL3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL3), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_REAL3
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_REAL4
subroutine iotk_scan_dat_aux_REAL4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  REAL (kind=__IOTK_REAL4), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_REAL4
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX1
subroutine iotk_scan_dat_aux_COMPLEX1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX1), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_COMPLEX1
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX2
subroutine iotk_scan_dat_aux_COMPLEX2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX2), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_COMPLEX2
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX3
subroutine iotk_scan_dat_aux_COMPLEX3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX3), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_COMPLEX3
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_COMPLEX4
subroutine iotk_scan_dat_aux_COMPLEX4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  COMPLEX (kind=__IOTK_COMPLEX4), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_COMPLEX4
#endif
# 95 "iotk_dat_interf.spp"
#ifdef __IOTK_CHARACTER1
subroutine iotk_scan_dat_aux_CHARACTER1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer,         intent(in)  :: unit
# 101 "iotk_dat_interf.spp"
  CHARACTER (kind=__IOTK_CHARACTER1,len=*), intent(out) :: dat(:)
  integer,         intent(in)  :: rkind
  integer,         intent(in)  :: rlen
  character(*),    intent(in)  :: fmt
  integer,         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_CHARACTER1
#endif
# 111 "iotk_dat_interf.spp"
end interface

end module iotk_dat_interf
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_error_interf.spp"

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
  character(len=*), intent(out):: val
end subroutine iotk_error_scan_character_e
subroutine iotk_error_scan_character_i(ierr,name,val)
  use iotk_base
  implicit none
  integer, intent(in) :: ierr
  character(len=*), intent(in) :: name
  character(len=*), intent(out):: val
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
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_files_interf.spp"

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
subroutine iotk_copyfile_x(source,dest,source_unit,dest_unit,ierr)
  implicit none
  character(len=*), optional, intent(in) :: source
  character(len=*), optional, intent(in) :: dest
  integer,          optional, intent(in) :: source_unit
  integer,          optional, intent(in) :: dest_unit
  integer,          optional, intent(out):: ierr
end subroutine iotk_copyfile_x
end interface

interface iotk_link
subroutine iotk_link_x(unit,name,file,binary,raw,create,ierr)
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: name
  character(len=*),           intent(in)  :: file
  logical,          optional, intent(in)  :: binary
  logical,          optional, intent(in)  :: raw
  logical,          optional, intent(in)  :: create
  integer,          optional, intent(out) :: ierr
end subroutine iotk_link_x
end interface

interface iotk_open_write
subroutine iotk_open_write_x(unit,file,attr,binary,new,raw,root,skip_root,skip_head,ierr)
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*), optional, intent(in)  :: file
  character(len=*), optional, intent(in)  :: attr
  logical,          optional, intent(in)  :: binary
  logical,          optional, intent(in)  :: new
  logical,          optional, intent(in)  :: raw
  character(len=*), optional, intent(in)  :: root
  logical,          optional, intent(in)  :: skip_root
  logical,          optional, intent(in)  :: skip_head
  integer,          optional, intent(out) :: ierr
end subroutine iotk_open_write_x
end interface

interface iotk_close_write
subroutine iotk_close_write_x(unit,ierr)
  implicit none
  integer,                intent(in)  :: unit
  integer,      optional, intent(out) :: ierr
end subroutine iotk_close_write_x
end interface

interface iotk_open_read
subroutine iotk_open_read_x(unit,file,attr,binary,raw,root,ierr)
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*), optional, intent(in)  :: file
  logical,          optional, intent(in)  :: binary 
  logical,          optional, intent(in)  :: raw
  character(len=*), optional, intent(out) :: attr
  character(len=*), optional, intent(out) :: root
  integer,          optional, intent(out) :: ierr
end subroutine iotk_open_read_x
end interface

interface iotk_close_read
subroutine iotk_close_read_x(unit,ierr)
  implicit none
  integer,                intent(in)  :: unit
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
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_fmt_interf.spp"

module iotk_fmt_interf
use iotk_base
implicit none
private

public :: iotk_basefmt
public :: iotk_wfmt

interface iotk_basefmt
function iotk_basefmt_x(type,ikind,ilen)
  implicit none
  character(100)           :: iotk_basefmt_x
  integer,      intent(in) :: ikind,ilen
  character(*), intent(in) :: type
end function iotk_basefmt_x
end interface

interface iotk_wfmt
function iotk_wfmt_x(type,ikind,isize,ilen)
  implicit none
  integer,       intent(in)  :: ikind
  character(*),  intent(in)  :: type
  integer,       intent(in)  :: isize
  integer,       intent(in)  :: ilen
  character(150)             :: iotk_wfmt_x
end function iotk_wfmt_x
end interface

end module iotk_fmt_interf
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_misc_interf.spp"

module iotk_misc_interf
use iotk_base
implicit none
private

public :: iotk_copy_tag
public :: iotk_parse_dat
public :: iotk_set_options
public :: iotk_get_options
public :: iotk_copy_dat_aux
public :: iotk_copy_dat
public :: iotk_print_kinds
public :: iotk_check_iotk_attr
public :: iotk_index
public :: iotk_check_name
public :: iotk_tag_parse
public :: iotk_complete_filepath



! This module contains the interfaces to all iotk routines

interface iotk_copy_tag
subroutine iotk_copy_tag_x(source,dest,maxsize,ierr)
  implicit none
  integer,           intent(in)  :: source
  integer,           intent(in)  :: dest
  integer, optional, intent(in)  :: maxsize
  integer, optional, intent(out) :: ierr
end subroutine iotk_copy_tag_x
end interface

interface iotk_parse_dat
subroutine iotk_parse_dat_x(attr,type,kind,isize,len,fmt,ierr)
  implicit none
  character(*), intent(in)  :: attr
  character(*), intent(out) :: type
  integer,      intent(out) :: kind
  integer,      intent(out) :: isize
  integer,      intent(out) :: len
  character(*), intent(out) :: fmt
  integer,      intent(out) :: ierr
end subroutine iotk_parse_dat_x
end interface

interface iotk_set_options
subroutine iotk_set_options_x(unitmin,unitmax,getline_buffer,error_warn_overflow,ierr)
  implicit none
  integer, optional, intent(in) :: unitmin
  integer, optional, intent(in) :: unitmax
  integer, optional, intent(in) :: getline_buffer
  logical, optional, intent(in) :: error_warn_overflow
  integer, optional, intent(out):: ierr
end subroutine iotk_set_options_x
end interface

interface iotk_get_options
subroutine iotk_get_options_x(unitmin,unitmax,getline_buffer,error_warn_overflow)
  implicit none
  integer, optional, intent(out):: unitmin
  integer, optional, intent(out):: unitmax
  logical, optional, intent(out):: error_warn_overflow
  integer, optional, intent(out):: getline_buffer
end subroutine iotk_get_options_x
end interface

interface iotk_copy_dat_aux
subroutine iotk_copy_dat_aux_x(source,dest,source_binary,dest_binary,name,type,ikind,isize,len,fmt,ierr)
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
  character(iotk_vallenx) :: version
end subroutine iotk_check_iotk_attr_x
end interface


interface iotk_index
function iotk_index_scal(index)
  integer,           intent(in) :: index
  character(len=range(index)+3) :: iotk_index_scal
end function iotk_index_scal
function iotk_index_vec(index)
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
  character(len=*), intent(in) :: newfile
  character(len=*), intent(in) :: oldfile
  character(len=len(newfile)+len(oldfile)) :: iotk_complete_filepath_x
  character(len=len(oldfile)) :: prefix
end function iotk_complete_filepath_x
end interface

interface iotk_check_name
function iotk_check_name_x(name)
  character(len=*), intent(in) :: name
  logical                      :: iotk_check_name_x
end function iotk_check_name_x
end interface

end module iotk_misc_interf
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_str_interf.spp"

module iotk_str_interf
use iotk_base
implicit none
private

public :: iotk_escape
public :: iotk_deescape
public :: iotk_strcpy
public :: iotk_strcat
public :: iotk_strlen
public :: iotk_strpad
public :: iotk_strscan
public :: iotk_strcomp
public :: iotk_strtrim
public :: iotk_strlen_trim
public :: iotk_toupper
public :: iotk_tolower

interface iotk_toupper
function iotk_toupper_x(str)
  implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: iotk_toupper_x
end function iotk_toupper_x
end interface

interface iotk_tolower
function iotk_tolower_x(str)
  implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: iotk_tolower_x
end function
end interface

interface iotk_escape
subroutine iotk_escape_x(to,from)
  implicit none
  character(len=*), intent(in) :: from
  character(len=*), intent(out):: to
end subroutine
end interface

interface iotk_deescape
subroutine iotk_deescape_x(to,from,quot,apos)
  implicit none
  character(len=*), intent(in)  :: from
  character(len=*), intent(out) :: to
  logical, optional, intent(in) :: quot,apos
end subroutine iotk_deescape_x
end interface

interface iotk_strscan
function iotk_strscan_x(string,set,back)
  implicit none
  character(len=*),  intent(in) :: string
  character(len=*),  intent(in) :: set
  logical, optional, intent(in) :: back
  integer                       :: iotk_strscan_x
end function iotk_strscan_x
end interface


interface iotk_strtrim
function iotk_strtrim_x(str)
  implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: iotk_strtrim_x
end function iotk_strtrim_x 
end interface

interface iotk_strlen_trim
function iotk_strlen_trim_x(str)
  implicit none
  character(len=*), intent(in) :: str
  integer                      :: iotk_strlen_trim_x
end function iotk_strlen_trim_x
end interface

interface iotk_strlen
function iotk_strlen_x(str)
  implicit none
  character(len=*), intent(in) :: str
  integer :: iotk_strlen_x
end function iotk_strlen_x
end interface

interface iotk_strpad
function iotk_strpad_x(str)
  implicit none
  character(len=*), intent(in) :: str
  character(len=len(str))      :: iotk_strpad_x
end function iotk_strpad_x
end interface

interface iotk_strcpy
subroutine iotk_strcpy_x(to,from,ierr)
  implicit none
  character(len=*), intent(out):: to
  character(len=*), intent(in) :: from
  integer,          intent(out):: ierr
end subroutine iotk_strcpy_x
end interface

interface iotk_strcat
subroutine iotk_strcat_x(to,from,ierr)
  implicit none
  character(len=*), intent(inout):: to
  character(len=*), intent(in) :: from
  integer,          intent(out):: ierr
end subroutine iotk_strcat_x
end interface

interface iotk_strcomp
function iotk_strcomp_x(str1,str2)
  implicit none
  logical :: iotk_strcomp_x
  character(len=*), intent(in) :: str1,str2
end function iotk_strcomp_x
end interface

end module iotk_str_interf
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_unit_interf.spp"

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
subroutine iotk_inquire_x(unit,binary,ierr)
  implicit none
  integer, intent(in)  :: unit
  logical, intent(out) :: binary
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
  integer,                intent(in)  :: unit
  type(iotk_unit), optional, pointer :: pointer
end subroutine iotk_unit_get_x
end interface

end module iotk_unit_interf
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_scan_interf.spp"

module iotk_scan_interf
use iotk_base
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
subroutine iotk_scan_begin_x(unit,name,attr,found,ierr)
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(out) :: attr
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
end subroutine iotk_scan_begin_x
end interface

interface iotk_scan_end
subroutine iotk_scan_end_x(unit,name,ierr)
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  integer,      optional, intent(out) :: ierr
end subroutine iotk_scan_end_x
end interface

interface iotk_scan_pi
subroutine iotk_scan_pi_x(unit,name,attr,found,ierr)
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(out) :: attr
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
end subroutine iotk_scan_pi_x
end interface

interface iotk_scan_empty
subroutine iotk_scan_empty_x(unit,name,attr,found,ierr)
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(out) :: attr
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
end subroutine iotk_scan_empty_x
end interface

interface iotk_scan_tag
subroutine iotk_scan_tag_x(unit,direction,control,tag,binary,ierr)
  use iotk_base
  implicit none
  integer,                 intent(in)  :: unit
  integer,                 intent(in)  :: direction
  integer,                 intent(out) :: control
  character(iotk_taglenx), intent(out) :: tag
  logical,                 intent(in)  :: binary
  integer,                 intent(out) :: ierr
end subroutine iotk_scan_tag_x
end interface

interface iotk_scan
subroutine iotk_scan_x(unit,direction,control,name,attr,binary,found,ierr)
  use iotk_base
  implicit none
  integer,                 intent(in)  :: unit
  integer,                 intent(in)  :: direction
  integer,                 intent(in)  :: control
  character(iotk_namlenx), intent(in)  :: name
  character(iotk_attlenx), intent(out) :: attr
  logical,                 intent(in)  :: binary
  logical,                 intent(out) :: found
  integer,                 intent(out) :: ierr
end subroutine iotk_scan_x
end interface

interface iotk_getline
subroutine iotk_getline_x(unit,line,length,ierr)
  implicit none
  integer,            intent(in)  :: unit
  character(len=*),   intent(out) :: line
  integer, optional,  intent(out) :: length
  integer, optional,  intent(out) :: ierr
end subroutine iotk_getline_x
end interface

end module iotk_scan_interf
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_write_interf.spp"

module iotk_write_interf
use iotk_base
implicit none
private

public :: iotk_write_begin
public :: iotk_write_end
public :: iotk_write_pi
public :: iotk_write_comment
public :: iotk_write_empty
public :: iotk_write_tag

interface iotk_write_begin
subroutine iotk_write_begin_x(unit,name,attr,ierr)
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: name
  character(len=*), optional, intent(in)  :: attr
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_begin_x
end interface

interface iotk_write_end
subroutine iotk_write_end_x(unit,name,ierr)
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: name
  integer,          optional, intent(out) :: ierr
end subroutine iotk_write_end_x
end interface

interface iotk_write_pi
subroutine iotk_write_pi_x(unit,name,attr,ierr)
  implicit none
  integer,                    intent(in)  :: unit
  character(len=*),           intent(in)  :: name
  character(len=*), optional, intent(in)  :: attr
  integer,          optional, intent(out) :: ierr
end subroutine iotk_write_pi_x
end interface

interface iotk_write_comment
subroutine iotk_write_comment_x(unit,text,ierr)
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: text
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_comment_x
end interface

interface iotk_write_empty
subroutine iotk_write_empty_x(unit,name,attr,ierr)
  implicit none
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(in)  :: attr
  integer,      optional, intent(out) :: ierr
end subroutine iotk_write_empty_x
end interface

interface iotk_write_tag
subroutine iotk_write_tag_x(unit,control,tag,binary,indent,ierr)
  use iotk_base
  implicit none
  integer,                   intent(in)  :: unit
  integer,                   intent(in)  :: control
  character(iotk_taglenx),   intent(in)  :: tag
  logical,                   intent(in)  :: binary
  integer,                   intent(in)  :: indent
  integer,                   intent(out) :: ierr
end subroutine iotk_write_tag_x
end interface

end module iotk_write_interf
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004,2005 Giovanni Bussi
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

# 2 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_AUXMACROS
#define __IOTK_AUXMACROS

! The macros are defined with -D option or inside iotk_config.h
! The default values are set here
! Maximum rank of an array
#ifndef __IOTK_MAXRANK
#  define __IOTK_MAXRANK 7
#endif
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
! Kind for header in binary files
#ifndef __IOTK_HEADER_KIND
#  define __IOTK_HEADER_KIND selected_int_kind(8)
#endif
! Character (or eventually string) for newline
! It may be adjusted for particular systems
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)
#ifndef __IOTK_NEWLINE
#  define __IOTK_NEWLINE achar(10)
#endif
! Character for EOS
#ifndef __IOTK_EOS
#  define __IOTK_EOS achar(0)
#endif
! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
#define __IOTK_CHARACTER1 iotk_defkind_character
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_LOGICAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_LOGICAL1 iotk_defkind_LOGICAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_INTEGER4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_INTEGER1 iotk_defkind_INTEGER
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL1
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL2
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL3
# 46 "../include/iotk_auxmacros.spp"
#ifndef __IOTK_REAL4
# 48 "../include/iotk_auxmacros.spp"
#define __IOTK_REAL1 iotk_defkind_REAL
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 51 "../include/iotk_auxmacros.spp"
#endif
# 54 "../include/iotk_auxmacros.spp"
! Complex are treated indentically to reals
! These lines map the definitions.
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL1
#define __IOTK_COMPLEX1 __IOTK_REAL1
#else
#undef __IOTK_COMPLEX1
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL2
#define __IOTK_COMPLEX2 __IOTK_REAL2
#else
#undef __IOTK_COMPLEX2
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL3
#define __IOTK_COMPLEX3 __IOTK_REAL3
#else
#undef __IOTK_COMPLEX3
#endif
# 57 "../include/iotk_auxmacros.spp"
#ifdef __IOTK_REAL4
#define __IOTK_COMPLEX4 __IOTK_REAL4
#else
#undef __IOTK_COMPLEX4
#endif
# 63 "../include/iotk_auxmacros.spp"
! If the binary format is not defined, use *
#ifndef __IOTK_BINARY_FORMAT
#define __IOTK_BINARY_FORMAT "*"
#endif

! Some check 
#if __IOTK_MAXRANK > 7
#  error
#endif
#if __IOTK_MAXRANK < 1
#  error
#endif

#endif

# 30 "iotk_xtox_interf.spp"

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
  character(len=*),  intent(in)  :: a
  logical, optional, intent(out) :: check
  logical                        :: iotk_atol_x
end function iotk_atol_x
end interface

interface iotk_atoi
# 51 "iotk_xtox_interf.spp"
#ifdef __IOTK_INTEGER1
subroutine iotk_atoi1(i,a,check)
  use iotk_base
  implicit none
  integer(kind=__IOTK_INTEGER1), intent(out) :: i
  character(len=*),                    intent(in)  :: a
  logical, optional,                   intent(out) :: check
end subroutine iotk_atoi1
#endif
# 51 "iotk_xtox_interf.spp"
#ifdef __IOTK_INTEGER2
subroutine iotk_atoi2(i,a,check)
  use iotk_base
  implicit none
  integer(kind=__IOTK_INTEGER2), intent(out) :: i
  character(len=*),                    intent(in)  :: a
  logical, optional,                   intent(out) :: check
end subroutine iotk_atoi2
#endif
# 51 "iotk_xtox_interf.spp"
#ifdef __IOTK_INTEGER3
subroutine iotk_atoi3(i,a,check)
  use iotk_base
  implicit none
  integer(kind=__IOTK_INTEGER3), intent(out) :: i
  character(len=*),                    intent(in)  :: a
  logical, optional,                   intent(out) :: check
end subroutine iotk_atoi3
#endif
# 51 "iotk_xtox_interf.spp"
#ifdef __IOTK_INTEGER4
subroutine iotk_atoi4(i,a,check)
  use iotk_base
  implicit none
  integer(kind=__IOTK_INTEGER4), intent(out) :: i
  character(len=*),                    intent(in)  :: a
  logical, optional,                   intent(out) :: check
end subroutine iotk_atoi4
#endif
# 61 "iotk_xtox_interf.spp"
end interface

interface iotk_itoa
# 65 "iotk_xtox_interf.spp"
#ifdef __IOTK_INTEGER1
function iotk_itoa1(i,length)
  use iotk_base
  implicit none
  integer(kind=__IOTK_INTEGER1), intent(in)  :: i
  integer, optional, intent(out)                 :: length
  character(len=range(i)+2)                      :: iotk_itoa1
end function iotk_itoa1
#endif
# 65 "iotk_xtox_interf.spp"
#ifdef __IOTK_INTEGER2
function iotk_itoa2(i,length)
  use iotk_base
  implicit none
  integer(kind=__IOTK_INTEGER2), intent(in)  :: i
  integer, optional, intent(out)                 :: length
  character(len=range(i)+2)                      :: iotk_itoa2
end function iotk_itoa2
#endif
# 65 "iotk_xtox_interf.spp"
#ifdef __IOTK_INTEGER3
function iotk_itoa3(i,length)
  use iotk_base
  implicit none
  integer(kind=__IOTK_INTEGER3), intent(in)  :: i
  integer, optional, intent(out)                 :: length
  character(len=range(i)+2)                      :: iotk_itoa3
end function iotk_itoa3
#endif
# 65 "iotk_xtox_interf.spp"
#ifdef __IOTK_INTEGER4
function iotk_itoa4(i,length)
  use iotk_base
  implicit none
  integer(kind=__IOTK_INTEGER4), intent(in)  :: i
  integer, optional, intent(out)                 :: length
  character(len=range(i)+2)                      :: iotk_itoa4
end function iotk_itoa4
#endif
# 75 "iotk_xtox_interf.spp"
end interface

interface iotk_ltoa
# 79 "iotk_xtox_interf.spp"
#ifdef __IOTK_LOGICAL1
function iotk_ltoa1(l)
  use iotk_base
  implicit none
  logical(kind=__IOTK_LOGICAL1), intent(in) :: l
  character                                     :: iotk_ltoa1
end function iotk_ltoa1
#endif
# 79 "iotk_xtox_interf.spp"
#ifdef __IOTK_LOGICAL2
function iotk_ltoa2(l)
  use iotk_base
  implicit none
  logical(kind=__IOTK_LOGICAL2), intent(in) :: l
  character                                     :: iotk_ltoa2
end function iotk_ltoa2
#endif
# 79 "iotk_xtox_interf.spp"
#ifdef __IOTK_LOGICAL3
function iotk_ltoa3(l)
  use iotk_base
  implicit none
  logical(kind=__IOTK_LOGICAL3), intent(in) :: l
  character                                     :: iotk_ltoa3
end function iotk_ltoa3
#endif
# 79 "iotk_xtox_interf.spp"
#ifdef __IOTK_LOGICAL4
function iotk_ltoa4(l)
  use iotk_base
  implicit none
  logical(kind=__IOTK_LOGICAL4), intent(in) :: l
  character                                     :: iotk_ltoa4
end function iotk_ltoa4
#endif
# 88 "iotk_xtox_interf.spp"
end interface

end module iotk_xtox_interf
