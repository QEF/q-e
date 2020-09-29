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
public :: iotk_str_clean

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
  character(len=*), intent(in)  :: from
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: to
#else
  character(len=*), intent(out) :: to
#endif
end subroutine
end interface

interface iotk_deescape
subroutine iotk_deescape_x(to,from,quot,apos)
  implicit none
  character(len=*), intent(in)  :: from
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: to
#else
  character(len=*), intent(out) :: to
#endif
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
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: to
#else
  character(len=*), intent(out) :: to
#endif
  character(len=*), intent(in)  :: from
  integer,          intent(out) :: ierr
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

interface iotk_str_clean
subroutine iotk_str_clean_x(str)
! transforms all characters which are separators in blanks
  implicit none
  character(len=*), intent(inout) :: str
end subroutine iotk_str_clean_x
end interface

end module iotk_str_interf
