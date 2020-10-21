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


! ERROR ROUTINES
subroutine iotk_error_init_e(error)
  use iotk_base
  implicit none
  type(iotk_error), intent(out) :: error
  nullify(error%str)
end subroutine iotk_error_init_e

subroutine iotk_error_init_i(ierr)
  implicit none
  integer, intent(out) :: ierr
  ierr = 0
end subroutine iotk_error_init_i

subroutine iotk_error_clear_e(error)
  use iotk_base
  implicit none
  type(iotk_error), intent(inout) :: error
  if(associated(error%str)) deallocate(error%str)
end subroutine iotk_error_clear_e

subroutine iotk_error_clear_i(ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(inout) :: ierr
  if(abs(ierr)>0 .and. abs(ierr)<=iotk_error_pool_size) then
    if(iotk_error_pool_used(abs(ierr))) then
      call iotk_error_clear(iotk_error_pool(abs(ierr)))
      iotk_error_pool_used(abs(ierr)) = .false.
      iotk_error_pool_order(abs(ierr)) = 0
    end if
  end if
  ierr = 0
end subroutine iotk_error_clear_i

function iotk_error_add_x()
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer :: i,ii(1),order
  integer :: iotk_error_add_x
  do i = 1 , iotk_error_pool_size
    if(.not. iotk_error_pool_used(i)) exit
  end do
  if(i>iotk_error_pool_size) then
    order=0
    do order=1,iotk_error_pool_size
      ii = minloc(iotk_error_pool_order,iotk_error_pool_order>=order)
      iotk_error_pool_order(ii(1)) = order
    end do
    if(iotk_error_warn_overflow) then
      write(iotk_error_unit,*) "Warning: ERROR OVERFLOW"
      call iotk_error_print(iotk_error_pool(iotk_error_pool_size),iotk_error_unit)
    end if
    ii = minloc(iotk_error_pool_order)
    i = ii(1)
    call iotk_error_clear(iotk_error_pool(i))
  end if
  iotk_error_pool_order(i) = maxval(iotk_error_pool_order)+1
  iotk_error_pool_used(i) = .true.
  call iotk_error_init(iotk_error_pool(i))
  iotk_error_add_x=i
end function iotk_error_add_x


subroutine iotk_error_append_e(error,str)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: str
  character, pointer :: tmp(:)
  integer :: i,strlen
  strlen = min(len(str),iotk_error_linelength)
  if(.not.associated(error%str)) then
    allocate(error%str(strlen+1))
    do i = 1 , strlen
      error%str(i) = str(i:i)
    end do
    error%str(strlen+1) = iotk_eos
  else
    tmp => error%str
    allocate(error%str(size(tmp)+strlen+1))
    error%str (1:size(tmp)) = tmp
    do i = 1 , strlen
      error%str (size(tmp)+i) = str(i:i)
    end do
    error%str(size(tmp)+strlen+1) = iotk_eos
    deallocate(tmp)
  end if
end subroutine iotk_error_append_e

subroutine iotk_error_append_i(ierr,str)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(inout) :: ierr
  character(len=*), intent(in)    :: str
  if(ierr==0) ierr = iotk_error_add()
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_append(iotk_error_pool(abs(ierr)),str)
end subroutine iotk_error_append_i

subroutine iotk_error_print_e(error,unit)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(in) :: error
  integer,          intent(in) :: unit
  integer :: i
  if(.not.associated(error%str)) return
  do i=1,size(error%str)
    if(error%str(i)==iotk_eos) then
      write(unit,"(a)")
    else
      write(unit,"(a)",advance='no') error%str(i)
    end if
  end do
end subroutine iotk_error_print_e

subroutine iotk_error_print_i(ierr,unit)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(in) :: ierr
  integer, intent(in) :: unit
  if(ierr==0) return
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_print(iotk_error_pool(abs(ierr)),unit)
end subroutine iotk_error_print_i

subroutine iotk_error_issue_e(error,sub,file,line)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: sub
  character(len=*), intent(in)    :: file
  integer,          intent(in)    :: line
  call iotk_error_append(error,"# ERROR IN: "//trim(sub)//" ("//trim(file)//":"//trim(iotk_itoa(line))//")")
end subroutine iotk_error_issue_e

subroutine iotk_error_issue_i(ierr,sub,file,line)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer,          intent(inout) :: ierr
  character(len=*), intent(in)    :: sub
  character(len=*), intent(in)    :: file
  integer,          intent(in)    :: line
  if(ierr==0) ierr = iotk_error_add()
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_issue(iotk_error_pool(abs(ierr)),sub,file,line)
end subroutine iotk_error_issue_i

subroutine iotk_error_msg_e(error,msg)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: msg
  call iotk_error_append(error,"# "//msg)
end subroutine iotk_error_msg_e

subroutine iotk_error_msg_i(ierr,msg)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer,          intent(inout) :: ierr
  character(len=*), intent(in)    :: msg
  if(ierr==0) ierr = iotk_error_add()
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_msg(iotk_error_pool(abs(ierr)),msg)
end subroutine iotk_error_msg_i

function iotk_error_check_e(error)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(in) :: error
  logical :: iotk_error_check_e
  iotk_error_check_e = .false.
  if(associated(error%str)) iotk_error_check_e = .true.
end function iotk_error_check_e

function iotk_error_check_i(ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(in) :: ierr
  logical :: iotk_error_check_i
  iotk_error_check_i = .false.
  if(ierr==0) return
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  iotk_error_check_i = .true.
end function iotk_error_check_i

subroutine iotk_error_write_character_e(error,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: name
  character(len=*), intent(in)    :: val
  integer :: namelen,vallen
  namelen=verify(name,alphabet_//numbers//".()%")-1
  if(namelen<0) namelen=len(name)
  vallen =scan  (val,iotk_eos)-1
  if(vallen<0) vallen=len(val)
  call iotk_error_append(error,name(1:namelen)//"="//val(1:vallen))
end subroutine iotk_error_write_character_e

subroutine iotk_error_write_character_i(ierr,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(inout) :: ierr
  character(len=*), intent(in)    :: name
  character(len=*), intent(in)    :: val
  if(ierr==0) ierr = iotk_error_add()
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_write(iotk_error_pool(abs(ierr)),name,val)
end subroutine iotk_error_write_character_i

subroutine iotk_error_write_logical_e(error,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: name
  logical,          intent(in)    :: val
  integer :: namelen
  character :: valc
  namelen=verify(name,alphabet_//numbers//".()%")-1
  if(namelen<0) namelen=len(name)
  valc="F"
  if(val) valc="T"
  call iotk_error_append(error,name(1:namelen)//"="//valc)
end subroutine iotk_error_write_logical_e

subroutine iotk_error_write_logical_i(ierr,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(inout) :: ierr
  character(len=*), intent(in)    :: name
  logical,          intent(in)    :: val
  if(ierr==0) ierr = iotk_error_add()
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_write(iotk_error_pool(abs(ierr)),name,val)
end subroutine iotk_error_write_logical_i 

subroutine iotk_error_write_integer_e(error,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(inout) :: error
  character(len=*), intent(in)    :: name
  integer,          intent(in)    :: val
  integer :: namelen
  namelen=verify(name,alphabet_//numbers//".()%")-1
  if(namelen<0) namelen=len(name)
  call iotk_error_append(error,name(1:namelen)//"="//trim(iotk_itoa(val)))
end subroutine iotk_error_write_integer_e

subroutine iotk_error_write_integer_i(ierr,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(inout) :: ierr
  character(len=*), intent(in)    :: name
  integer,          intent(in)    :: val
  if(ierr==0) ierr = iotk_error_add()
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_write(iotk_error_pool(abs(ierr)),name,val)
end subroutine iotk_error_write_integer_i

subroutine iotk_error_scan_character_e(error,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(in) :: error
  character(len=*), intent(in) :: name
#ifdef __IOTK_WORKAROUND6
  character(len=*)             :: val
#else
  character(len=*), intent(out):: val
#endif
  integer :: i1,i2,i3
  logical :: eos,found
  i2 = 0
  val=""
  found = .false.
  if(.not.associated(error%str)) return
  do i1 = size(error%str) , 0 , -1
    eos = .false.
    if(i1==0) eos = .true.
    if(.not.eos) then
      if(error%str(i1)==iotk_eos) eos = .true.
    end if
    if(eos) then
      do i2=1,len(name)
        if(i1+i2 > size(error%str)) goto 1
        if(error%str(i1+i2)/=name(i2:i2)) goto 1
      end do
      if(i1+i2 > size(error%str)) goto 1
      if(error%str(i1+i2)/="=") goto 1 
      found=.true.
      exit
    end if
1 continue
  end do
  val=""
  if(found) then
    do i3=1,len(val)
      if(i1+i2+i3>size(error%str)) exit
      if(error%str(i1+i2+i3)==iotk_eos) exit
      val(i3:i3)=error%str(i1+i2+i3)
    end do
  end if
end subroutine iotk_error_scan_character_e

subroutine iotk_error_scan_character_i(ierr,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(in) :: ierr
  character(len=*), intent(in) :: name
#ifdef __IOTK_WORKAROUND6
  character(len=*)             :: val
#else
  character(len=*), intent(out):: val
#endif
  val = ""
  if(ierr==0) return
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_scan(iotk_error_pool(abs(ierr)),name,val)
end subroutine iotk_error_scan_character_i

subroutine iotk_error_scan_logical_e(error,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(in) :: error
  character(len=*), intent(in) :: name
  logical,          intent(out):: val
  character :: valc
  val = .false.
  call iotk_error_scan(error,name,valc)
  if(valc=="T" .or. valc=="t") val=.true.
end subroutine iotk_error_scan_logical_e

subroutine iotk_error_scan_logical_i(ierr,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(in) :: ierr
  character(len=*), intent(in) :: name
  logical,          intent(out):: val
  val = .false.
  if(ierr==0) return
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_scan(iotk_error_pool(abs(ierr)),name,val)
end subroutine iotk_error_scan_logical_i

subroutine iotk_error_scan_integer_e(error,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  type(iotk_error), intent(in) :: error
  character(len=*), intent(in) :: name
  integer,          intent(out):: val
  character(range(val)+2) :: valc
  call iotk_error_scan(error,name,valc)
  call iotk_atoi(val,valc)
end subroutine iotk_error_scan_integer_e

subroutine iotk_error_scan_integer_i(ierr,name,val)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer, intent(in) :: ierr
  character(len=*), intent(in) :: name
  integer,          intent(out):: val
  val = 0
  if(ierr==0) return
  if(abs(ierr)>iotk_error_pool_size) return
  if(.not. iotk_error_pool_used(abs(ierr))) return
  call iotk_error_scan(iotk_error_pool(abs(ierr)),name,val)
end subroutine iotk_error_scan_integer_i

function iotk_error_pool_pending_x()
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  implicit none
  integer :: iotk_error_pool_pending_x
  iotk_error_pool_pending_x = count (iotk_error_pool_used)
end function iotk_error_pool_pending_x

subroutine iotk_error_handler_x(ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  integer, intent(in) :: ierr
  integer :: pending,i
#ifdef __IOTK_MPI_ABORT
  include 'mpif.h'
  integer :: ierrx
#endif
  if(ierr==0) return
  do i = 1 , iotk_error_linelength
    write(iotk_error_unit,"(a)",advance='no') "#"
  end do
  write(iotk_error_unit,*)
  pending = iotk_error_pool_pending()
  if(pending>1) then
    write(iotk_error_unit,"(a)") "# WARNING: there are pending errors"
    do i = 1 , iotk_error_pool_size
      if(iotk_error_pool_used(i) .and. i/=abs(ierr)) then
        write(iotk_error_unit,"(a)") "# PENDING ERROR (ierr="//trim(iotk_itoa(i))//")"
        call iotk_error_print(i,iotk_error_unit)
      end if
    end do
  end if
  write(iotk_error_unit,"(a)") "# FROM IOTK LIBRARY, VERSION "//trim(iotk_version)
  write(iotk_error_unit,"(a)") "# UNRECOVERABLE ERROR (ierr="//trim(iotk_itoa(ierr))//")"
  call iotk_error_print(ierr,0)
  do i = 1 , iotk_error_linelength
    write(iotk_error_unit,"(a)",advance='no') "#"
  end do
  write(iotk_error_unit,*)
#ifdef __IOTK_MPI_ABORT
  call MPI_Abort(MPI_COMM_WORLD,1,ierrx)
#else
  stop
#endif
end subroutine iotk_error_handler_x
