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

module iotk_dat_interf
implicit none
private

public :: iotk_write_dat
public :: iotk_scan_dat
!-<
public :: iotk_scan_dat_inside
!->
public :: iotk_scan_dat_aux


interface iotk_write_dat
#ifdef __IOTK_LOGICAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL1_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL1_7
#endif
#endif
#ifdef __IOTK_LOGICAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL2_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL2_7
#endif
#endif
#ifdef __IOTK_LOGICAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL3_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL3_7
#endif
#endif
#ifdef __IOTK_LOGICAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_LOGICAL4_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_LOGICAL4_7
#endif
#endif
#ifdef __IOTK_INTEGER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER1_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER1_7
#endif
#endif
#ifdef __IOTK_INTEGER2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER2_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER2_7
#endif
#endif
#ifdef __IOTK_INTEGER3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER3_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER3_7
#endif
#endif
#ifdef __IOTK_INTEGER4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_INTEGER4_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_INTEGER4_7
#endif
#endif
#ifdef __IOTK_REAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL1_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL1_7
#endif
#endif
#ifdef __IOTK_REAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL2_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL2_7
#endif
#endif
#ifdef __IOTK_REAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL3_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL3_7
#endif
#endif
#ifdef __IOTK_REAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_REAL4_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_REAL4_7
#endif
#endif
#ifdef __IOTK_COMPLEX1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX1_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX1_7
#endif
#endif
#ifdef __IOTK_COMPLEX2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX2_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX2_7
#endif
#endif
#ifdef __IOTK_COMPLEX3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX3_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX3_7
#endif
#endif
#ifdef __IOTK_COMPLEX4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_COMPLEX4_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_COMPLEX4_7
#endif
#endif
#ifdef __IOTK_CHARACTER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_0(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: dat 
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_1(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: dat (:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_2(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: dat (:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_3(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: dat (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_4(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: dat (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_5(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: dat (:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_6(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: dat (:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_write_dat_CHARACTER1_7(unit,name,dat,dummy,attr,columns,sep,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: dat (:,:,:,:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  character(len=*),                optional, intent(in)  :: sep
  character(len=*),                optional, intent(in)  :: fmt
  character(len=*),                optional, intent(in)  :: attr
  integer,                         optional, intent(in)  :: columns
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_dat_CHARACTER1_7
#endif
#endif
end interface

interface iotk_scan_dat
#ifdef __IOTK_LOGICAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat 
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL1_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL1_7
#endif
#endif
#ifdef __IOTK_LOGICAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat 
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL2_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL2_7
#endif
#endif
#ifdef __IOTK_LOGICAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat 
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL3_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL3_7
#endif
#endif
#ifdef __IOTK_LOGICAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat 
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_LOGICAL4_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_LOGICAL4_7
#endif
#endif
#ifdef __IOTK_INTEGER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat 
#else
  INTEGER(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER1_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER1_7
#endif
#endif
#ifdef __IOTK_INTEGER2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat 
#else
  INTEGER(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER2_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER2_7
#endif
#endif
#ifdef __IOTK_INTEGER3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat 
#else
  INTEGER(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER3_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER3_7
#endif
#endif
#ifdef __IOTK_INTEGER4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat 
#else
  INTEGER(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_INTEGER4_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_INTEGER4_7
#endif
#endif
#ifdef __IOTK_REAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL1_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL1_7
#endif
#endif
#ifdef __IOTK_REAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL2_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL2_7
#endif
#endif
#ifdef __IOTK_REAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL3_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL3_7
#endif
#endif
#ifdef __IOTK_REAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_REAL4_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_REAL4_7
#endif
#endif
#ifdef __IOTK_COMPLEX1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat 
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX1_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX1_7
#endif
#endif
#ifdef __IOTK_COMPLEX2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat 
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX2_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX2_7
#endif
#endif
#ifdef __IOTK_COMPLEX3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat 
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX3_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX3_7
#endif
#endif
#ifdef __IOTK_COMPLEX4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat 
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_COMPLEX4_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_COMPLEX4_7
#endif
#endif
#ifdef __IOTK_CHARACTER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_0(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat 
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_1(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_2(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_3(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_4(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_5(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_6(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:,:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_CHARACTER1_7(unit,name,dat,dummy,attr,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
  character(len=*),                          intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:,:,:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
#ifdef __IOTK_WORKAROUND6
  character(len=*),                optional              :: attr
#else
  character(len=*),                optional, intent(out) :: attr
#endif
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_CHARACTER1_7
#endif
#endif
end interface

!-<
interface iotk_scan_dat_inside
#ifdef __IOTK_LOGICAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL1_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat 
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL1_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL1_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL1_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL1_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL1_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL1_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL1_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL1_7
#endif
#endif
#ifdef __IOTK_LOGICAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL2_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat 
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL2_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL2_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL2_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL2_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL2_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL2_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL2_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL2_7
#endif
#endif
#ifdef __IOTK_LOGICAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL3_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat 
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL3_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL3_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL3_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL3_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL3_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL3_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL3_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL3_7
#endif
#endif
#ifdef __IOTK_LOGICAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL4_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat 
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL4_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL4_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL4_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL4_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL4_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL4_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_LOGICAL4_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  LOGICAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_LOGICAL4_7
#endif
#endif
#ifdef __IOTK_INTEGER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER1_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat 
#else
  INTEGER(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER1_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER1_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER1_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER1_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER1_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER1_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER1_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER1_7
#endif
#endif
#ifdef __IOTK_INTEGER2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER2_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat 
#else
  INTEGER(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER2_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER2_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER2_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER2_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER2_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER2_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER2_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER2_7
#endif
#endif
#ifdef __IOTK_INTEGER3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER3_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat 
#else
  INTEGER(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER3_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER3_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER3_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER3_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER3_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER3_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER3_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER3_7
#endif
#endif
#ifdef __IOTK_INTEGER4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER4_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat 
#else
  INTEGER(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER4_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER4_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER4_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER4_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER4_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER4_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_INTEGER4_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  INTEGER(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_INTEGER4_7
#endif
#endif
#ifdef __IOTK_REAL1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL1_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL1_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL1_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL1_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL1_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL1_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL1_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL1_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL1_7
#endif
#endif
#ifdef __IOTK_REAL2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL2_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL2_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL2_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL2_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL2_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL2_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL2_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL2_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL2_7
#endif
#endif
#ifdef __IOTK_REAL3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL3_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL3_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL3_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL3_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL3_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL3_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL3_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL3_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL3_7
#endif
#endif
#ifdef __IOTK_REAL4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL4_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat 
#else
  REAL(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL4_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL4_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL4_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL4_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL4_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL4_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_REAL4_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  REAL(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_REAL4_7
#endif
#endif
#ifdef __IOTK_COMPLEX1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX1_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat 
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX1_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX1_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX1_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX1_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX1_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX1_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX1_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX1_7
#endif
#endif
#ifdef __IOTK_COMPLEX2
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX2_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat 
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX2_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX2_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX2_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX2_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX2_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX2_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX2_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX2_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX2_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX2_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX2_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX2_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX2_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX2_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX2_7
#endif
#endif
#ifdef __IOTK_COMPLEX3
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX3_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat 
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX3_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX3_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX3_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX3_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX3_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX3_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX3_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX3_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX3_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX3_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX3_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX3_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX3_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX3_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX3_7
#endif
#endif
#ifdef __IOTK_COMPLEX4
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX4_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat 
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX4_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX4_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX4_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX4_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX4_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX4_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX4_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX4_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX4_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX4_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX4_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX4_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX4_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_COMPLEX4_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)                        :: dat (:,:,:,:,:,:,:)
#else
  COMPLEX(kind=this_kind),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_COMPLEX4_7
#endif
#endif
#ifdef __IOTK_CHARACTER1
#if 0 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_CHARACTER1_0(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat 
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat 
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default 
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_CHARACTER1_0
#endif
#if 1 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_CHARACTER1_1(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_CHARACTER1_1
#endif
#if 2 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_CHARACTER1_2(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_CHARACTER1_2
#endif
#if 3 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_CHARACTER1_3(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_CHARACTER1_3
#endif
#if 4 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_CHARACTER1_4(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_CHARACTER1_4
#endif
#if 5 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_CHARACTER1_5(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_CHARACTER1_5
#endif
#if 6 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_CHARACTER1_6(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:,:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_CHARACTER1_6
#endif
#if 7 <= __IOTK_MAXRANK
subroutine iotk_scan_dat_inside_CHARACTER1_7(unit,dat,dummy,found,default,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                                   intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)                        :: dat (:,:,:,:,:,:,:)
#else
  CHARACTER(kind=this_kind,len=*),           intent(out) :: dat (:,:,:,:,:,:,:)
#endif
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default (:,:,:,:,:,:,:)
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_dat_inside_CHARACTER1_7
#endif
#endif
end interface
!->

interface iotk_scan_dat_aux
#ifdef __IOTK_LOGICAL1
subroutine iotk_scan_dat_aux_LOGICAL1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)              :: dat (:)
#else
  LOGICAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_LOGICAL1
#endif
#ifdef __IOTK_LOGICAL2
subroutine iotk_scan_dat_aux_LOGICAL2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL2
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)              :: dat (:)
#else
  LOGICAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_LOGICAL2
#endif
#ifdef __IOTK_LOGICAL3
subroutine iotk_scan_dat_aux_LOGICAL3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL3
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)              :: dat (:)
#else
  LOGICAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_LOGICAL3
#endif
#ifdef __IOTK_LOGICAL4
subroutine iotk_scan_dat_aux_LOGICAL4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL4
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  LOGICAL(kind=this_kind)              :: dat (:)
#else
  LOGICAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_LOGICAL4
#endif
#ifdef __IOTK_INTEGER1
subroutine iotk_scan_dat_aux_INTEGER1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)              :: dat (:)
#else
  INTEGER(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_INTEGER1
#endif
#ifdef __IOTK_INTEGER2
subroutine iotk_scan_dat_aux_INTEGER2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER2
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)              :: dat (:)
#else
  INTEGER(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_INTEGER2
#endif
#ifdef __IOTK_INTEGER3
subroutine iotk_scan_dat_aux_INTEGER3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER3
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)              :: dat (:)
#else
  INTEGER(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_INTEGER3
#endif
#ifdef __IOTK_INTEGER4
subroutine iotk_scan_dat_aux_INTEGER4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER4
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  INTEGER(kind=this_kind)              :: dat (:)
#else
  INTEGER(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_INTEGER4
#endif
#ifdef __IOTK_REAL1
subroutine iotk_scan_dat_aux_REAL1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)              :: dat (:)
#else
  REAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_REAL1
#endif
#ifdef __IOTK_REAL2
subroutine iotk_scan_dat_aux_REAL2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)              :: dat (:)
#else
  REAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_REAL2
#endif
#ifdef __IOTK_REAL3
subroutine iotk_scan_dat_aux_REAL3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL3
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)              :: dat (:)
#else
  REAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_REAL3
#endif
#ifdef __IOTK_REAL4
subroutine iotk_scan_dat_aux_REAL4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL4
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  REAL(kind=this_kind)              :: dat (:)
#else
  REAL(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_REAL4
#endif
#ifdef __IOTK_COMPLEX1
subroutine iotk_scan_dat_aux_COMPLEX1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)              :: dat (:)
#else
  COMPLEX(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_COMPLEX1
#endif
#ifdef __IOTK_COMPLEX2
subroutine iotk_scan_dat_aux_COMPLEX2(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)              :: dat (:)
#else
  COMPLEX(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_COMPLEX2
#endif
#ifdef __IOTK_COMPLEX3
subroutine iotk_scan_dat_aux_COMPLEX3(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX3
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)              :: dat (:)
#else
  COMPLEX(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_COMPLEX3
#endif
#ifdef __IOTK_COMPLEX4
subroutine iotk_scan_dat_aux_COMPLEX4(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX4
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=this_kind)              :: dat (:)
#else
  COMPLEX(kind=this_kind), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_COMPLEX4
#endif
#ifdef __IOTK_CHARACTER1
subroutine iotk_scan_dat_aux_CHARACTER1(unit,dat,rkind,rlen,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  integer,                         intent(in)  :: unit
#ifdef __IOTK_WORKAROUND6
  CHARACTER(kind=this_kind,len=*)              :: dat (:)
#else
  CHARACTER(kind=this_kind,len=*), intent(out) :: dat (:)
#endif
  integer,                         intent(in)  :: rkind
  integer,                         intent(in)  :: rlen
  character(len=*),                intent(in)  :: fmt
  integer,                         intent(out) :: ierr
end subroutine iotk_scan_dat_aux_CHARACTER1
#endif
end interface

end module iotk_dat_interf
