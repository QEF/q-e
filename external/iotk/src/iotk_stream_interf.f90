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



module iotk_stream_interf
implicit none
private

public :: iotk_stream_read
public :: iotk_stream_backspace

interface iotk_stream_read
  subroutine iotk_stream_read_x(unit,header,setpos,getpos,ierr)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: unit
    integer(iotk_header_kind),                  intent(out) :: header
    integer,                          optional, intent(in)  :: setpos
    integer,                          optional, intent(out) :: getpos
    integer,                          optional, intent(out) :: ierr
  end subroutine iotk_stream_read_x
#ifdef __IOTK_LOGICAL1
  subroutine iotk_stream_read_LOGICAL1(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    LOGICAL(kind=iotk_LOGICAL1),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_LOGICAL1
#endif
#ifdef __IOTK_LOGICAL2
  subroutine iotk_stream_read_LOGICAL2(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    LOGICAL(kind=iotk_LOGICAL2),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_LOGICAL2
#endif
#ifdef __IOTK_LOGICAL3
  subroutine iotk_stream_read_LOGICAL3(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    LOGICAL(kind=iotk_LOGICAL3),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_LOGICAL3
#endif
#ifdef __IOTK_LOGICAL4
  subroutine iotk_stream_read_LOGICAL4(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    LOGICAL(kind=iotk_LOGICAL4),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_LOGICAL4
#endif
#ifdef __IOTK_INTEGER1
  subroutine iotk_stream_read_INTEGER1(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    INTEGER(kind=iotk_INTEGER1),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_INTEGER1
#endif
#ifdef __IOTK_INTEGER2
  subroutine iotk_stream_read_INTEGER2(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    INTEGER(kind=iotk_INTEGER2),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_INTEGER2
#endif
#ifdef __IOTK_INTEGER3
  subroutine iotk_stream_read_INTEGER3(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    INTEGER(kind=iotk_INTEGER3),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_INTEGER3
#endif
#ifdef __IOTK_INTEGER4
  subroutine iotk_stream_read_INTEGER4(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    INTEGER(kind=iotk_INTEGER4),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_INTEGER4
#endif
#ifdef __IOTK_REAL1
  subroutine iotk_stream_read_REAL1(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    REAL(kind=iotk_REAL1),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_REAL1
#endif
#ifdef __IOTK_REAL2
  subroutine iotk_stream_read_REAL2(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    REAL(kind=iotk_REAL2),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_REAL2
#endif
#ifdef __IOTK_REAL3
  subroutine iotk_stream_read_REAL3(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    REAL(kind=iotk_REAL3),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_REAL3
#endif
#ifdef __IOTK_REAL4
  subroutine iotk_stream_read_REAL4(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    REAL(kind=iotk_REAL4),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_REAL4
#endif
#ifdef __IOTK_COMPLEX1
  subroutine iotk_stream_read_COMPLEX1(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    COMPLEX(kind=iotk_COMPLEX1),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_COMPLEX1
#endif
#ifdef __IOTK_COMPLEX2
  subroutine iotk_stream_read_COMPLEX2(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    COMPLEX(kind=iotk_COMPLEX2),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_COMPLEX2
#endif
#ifdef __IOTK_COMPLEX3
  subroutine iotk_stream_read_COMPLEX3(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    COMPLEX(kind=iotk_COMPLEX3),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_COMPLEX3
#endif
#ifdef __IOTK_COMPLEX4
  subroutine iotk_stream_read_COMPLEX4(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    COMPLEX(kind=iotk_COMPLEX4),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_COMPLEX4
#endif
#ifdef __IOTK_CHARACTER1
  subroutine iotk_stream_read_CHARACTER1(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    CHARACTER(kind=iotk_CHARACTER1,len=*),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_CHARACTER1
#endif
#ifdef __IOTK_CHARACTER2
  subroutine iotk_stream_read_CHARACTER2(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    CHARACTER(kind=iotk_CHARACTER2,len=*),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_CHARACTER2
#endif
#ifdef __IOTK_CHARACTER3
  subroutine iotk_stream_read_CHARACTER3(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    CHARACTER(kind=iotk_CHARACTER3,len=*),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_CHARACTER3
#endif
#ifdef __IOTK_CHARACTER4
  subroutine iotk_stream_read_CHARACTER4(unit,header,val,setpos,getpos,noval,ierr)
    use iotk_base
    implicit none
    integer,                                      intent(in)  :: unit
    integer(iotk_header_kind),                    intent(out) :: header
    CHARACTER(kind=iotk_CHARACTER4,len=*),      intent(out) :: val(:)
    integer,                            optional, intent(in)  :: setpos
    integer,                            optional, intent(out) :: getpos
    logical,                            optional, intent(in)  :: noval
    integer,                            optional, intent(out) :: ierr
  end subroutine iotk_stream_read_CHARACTER4
#endif
end interface

interface iotk_stream_backspace
  subroutine iotk_stream_backspace_x(unit,ierr)
    implicit none
    integer,           intent(in)  :: unit
    integer, optional, intent(out) :: ierr
  end subroutine iotk_stream_backspace_x
end interface

end module iotk_stream_interf
