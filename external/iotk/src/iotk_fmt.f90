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


function iotk_basefmt_x(type,ikind,ilen)
  use iotk_base
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  character(len=*), intent(in) :: type
  integer,          intent(in) :: ikind
  integer,          intent(in) :: ilen
  character(100)               :: iotk_basefmt_x
  integer :: nexp,exp,ndig,baselen
  logical, save :: first_call = .true.
#ifdef __IOTK_INTEGER1
  integer (kind=iotk_integer1) :: example_integer1 = 0
  character(46), save :: save_basefmt_integer1 = ""
#endif
#ifdef __IOTK_REAL1
  real (kind=iotk_real1) :: example_real1 = 0.0
  character(46), save :: save_basefmt_real1 = ""
#endif
#ifdef __IOTK_INTEGER2
  integer (kind=iotk_integer2) :: example_integer2 = 0
  character(46), save :: save_basefmt_integer2 = ""
#endif
#ifdef __IOTK_REAL2
  real (kind=iotk_real2) :: example_real2 = 0.0
  character(46), save :: save_basefmt_real2 = ""
#endif
#ifdef __IOTK_INTEGER3
  integer (kind=iotk_integer3) :: example_integer3 = 0
  character(46), save :: save_basefmt_integer3 = ""
#endif
#ifdef __IOTK_REAL3
  real (kind=iotk_real3) :: example_real3 = 0.0
  character(46), save :: save_basefmt_real3 = ""
#endif
#ifdef __IOTK_INTEGER4
  integer (kind=iotk_integer4) :: example_integer4 = 0
  character(46), save :: save_basefmt_integer4 = ""
#endif
#ifdef __IOTK_REAL4
  real (kind=iotk_real4) :: example_real4 = 0.0
  character(46), save :: save_basefmt_real4 = ""
#endif
  if(first_call) then
#ifdef __IOTK_INTEGER1
    baselen = range(example_integer1) + 1
    save_basefmt_integer1 = "(i"//trim(iotk_itoa(baselen))//")"
#endif
#ifdef __IOTK_REAL1
    ndig = precision(example_real1)+1
    exp = range(example_real1)+1
    nexp = 1
    do
      if(exp < 10) exit
      exp = exp / 10
      nexp = nexp + 1
    end do
    baselen = nexp+ndig-1+5
    save_basefmt_real1 = "(ES"//trim(iotk_itoa(baselen))//"." &
                //trim(iotk_itoa(ndig-1))//"E"//trim(iotk_itoa(nexp))//")"

#endif
#ifdef __IOTK_INTEGER2
    baselen = range(example_integer2) + 1
    save_basefmt_integer2 = "(i"//trim(iotk_itoa(baselen))//")"
#endif
#ifdef __IOTK_REAL2
    ndig = precision(example_real2)+1
    exp = range(example_real2)+1
    nexp = 1
    do
      if(exp < 10) exit
      exp = exp / 10
      nexp = nexp + 1
    end do
    baselen = nexp+ndig-1+5
    save_basefmt_real2 = "(ES"//trim(iotk_itoa(baselen))//"." &
                //trim(iotk_itoa(ndig-1))//"E"//trim(iotk_itoa(nexp))//")"

#endif
#ifdef __IOTK_INTEGER3
    baselen = range(example_integer3) + 1
    save_basefmt_integer3 = "(i"//trim(iotk_itoa(baselen))//")"
#endif
#ifdef __IOTK_REAL3
    ndig = precision(example_real3)+1
    exp = range(example_real3)+1
    nexp = 1
    do
      if(exp < 10) exit
      exp = exp / 10
      nexp = nexp + 1
    end do
    baselen = nexp+ndig-1+5
    save_basefmt_real3 = "(ES"//trim(iotk_itoa(baselen))//"." &
                //trim(iotk_itoa(ndig-1))//"E"//trim(iotk_itoa(nexp))//")"

#endif
#ifdef __IOTK_INTEGER4
    baselen = range(example_integer4) + 1
    save_basefmt_integer4 = "(i"//trim(iotk_itoa(baselen))//")"
#endif
#ifdef __IOTK_REAL4
    ndig = precision(example_real4)+1
    exp = range(example_real4)+1
    nexp = 1
    do
      if(exp < 10) exit
      exp = exp / 10
      nexp = nexp + 1
    end do
    baselen = nexp+ndig-1+5
    save_basefmt_real4 = "(ES"//trim(iotk_itoa(baselen))//"." &
                //trim(iotk_itoa(ndig-1))//"E"//trim(iotk_itoa(nexp))//")"

#endif
    first_call = .false.
  end if
  select case(type)
  case("LOGICAL")
    iotk_basefmt_x = "(l1)"
  case("INTEGER")
    select case(ikind)
#ifdef __IOTK_INTEGER1
    case(iotk_integer1)
      iotk_basefmt_x = save_basefmt_integer1
#endif
#ifdef __IOTK_INTEGER2
    case(iotk_integer2)
      iotk_basefmt_x = save_basefmt_integer2
#endif
#ifdef __IOTK_INTEGER3
    case(iotk_integer3)
      iotk_basefmt_x = save_basefmt_integer3
#endif
#ifdef __IOTK_INTEGER4
    case(iotk_integer4)
      iotk_basefmt_x = save_basefmt_integer4
#endif
    end select
  case("REAL")
    select case(ikind)
#ifdef __IOTK_REAL1
    case(iotk_real1)
      iotk_basefmt_x = save_basefmt_real1
#endif
#ifdef __IOTK_REAL2
    case(iotk_real2)
      iotk_basefmt_x = save_basefmt_real2
#endif
#ifdef __IOTK_REAL3
    case(iotk_real3)
      iotk_basefmt_x = save_basefmt_real3
#endif
#ifdef __IOTK_REAL4
    case(iotk_real4)
      iotk_basefmt_x = save_basefmt_real4
#endif
    end select
  case("COMPLEX")
    select case(ikind)
#ifdef __IOTK_REAL1
    case(iotk_real1)
      iotk_basefmt_x = "("//trim(save_basefmt_real1)//",',',"//trim(save_basefmt_real1)//")"
#endif
#ifdef __IOTK_REAL2
    case(iotk_real2)
      iotk_basefmt_x = "("//trim(save_basefmt_real2)//",',',"//trim(save_basefmt_real2)//")"
#endif
#ifdef __IOTK_REAL3
    case(iotk_real3)
      iotk_basefmt_x = "("//trim(save_basefmt_real3)//",',',"//trim(save_basefmt_real3)//")"
#endif
#ifdef __IOTK_REAL4
    case(iotk_real4)
      iotk_basefmt_x = "("//trim(save_basefmt_real4)//",',',"//trim(save_basefmt_real4)//")"
#endif
    end select
  case("CHARACTER")
    if(ilen>=0) then
      iotk_basefmt_x = "(a"//trim(iotk_itoa(ilen))//")"
    else
      iotk_basefmt_x = "(a)"
    end if
  end select
end function iotk_basefmt_x

function iotk_wfmt_x(type,ikind,isize,ilen,sep)
  use iotk_base
  use iotk_xtox_interf
  use iotk_fmt_interf
  use iotk_misc_interf
  use iotk_str_interf
  implicit none
  character(len=*), intent(in)  :: type
  integer,          intent(in)  :: ikind
  integer,          intent(in)  :: isize
  integer,          intent(in)  :: ilen
  character(len=*), intent(in)  :: sep
  character(150)                :: iotk_wfmt_x
  if(isize==1) then
    iotk_wfmt_x = "("//trim(iotk_basefmt(type,ikind,ilen))//")"
  else
    iotk_wfmt_x = "("//trim(iotk_itoa(isize))//"("//trim(iotk_basefmt(type,ikind,ilen)) &
                //",:,'"//sep(1:iotk_strlen(sep))//"'))"
  end if
end function iotk_wfmt_x
