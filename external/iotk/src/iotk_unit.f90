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


subroutine iotk_free_unit_x(unit,ierr)
  use iotk_base
  use iotk_error_interf
  implicit none
! This subroutine sets 'unit' to the number of
! an I/O unit which is free (i.e. not already opened).
! The search is carried out starting from unit
! 'unitmin' in a range of 'nsearch' units.
! The starting unit for the search is increased at each
! call, so that a number of subsequent ask can be done
! obtaining different units.
  integer,           intent(out) :: unit
  integer, optional, intent(out) :: ierr
  integer, save :: offset = 0
  logical       :: opened,exist
  integer       :: isearch,nsearch,unitmin
  integer       :: ierrl
  integer       :: iostat
  iostat = 0
  unitmin = iotk_unitmin
  nsearch = iotk_unitmax - iotk_unitmin + 1
  ierrl = 0 
  do isearch=0,nsearch-1
    unit = modulo(isearch+offset,nsearch) + unitmin
    inquire(unit=unit,opened=opened,exist=exist,iostat=iostat)
    if(iostat/=0) then
      call iotk_error_issue(ierrl,"iotk_free_unit",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.16 ")
call iotk_error_msg(ierrl,'Error inquiring')
call iotk_error_write(ierrl,"unit",unit)
call iotk_error_write(ierrl,"iostat",iostat)
      goto 1
    end if
    if((.not.opened .and. exist) .or. iostat/=0) exit
  end do
  if(isearch>=nsearch) then
    call iotk_error_issue(ierrl,"iotk_free_unit",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.16 ")
call iotk_error_msg(ierrl,'There are no units left')
call iotk_error_write(ierrl,"iotk_unitmin",iotk_unitmin)
call iotk_error_write(ierrl,"iotk_unitmax",iotk_unitmax)
call iotk_error_write(ierrl,"offset",offset)
    goto 1
  end if 
  offset = modulo(unit - unitmin + 1,nsearch)
1 continue
  if(present(ierr)) then 
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_free_unit_x

function iotk_phys_unit_x(unit) result(result)
  use iotk_base
  use iotk_unit_interf
  use iotk_unit_list_module
  implicit none
  integer, intent(in) :: unit
  integer :: result
  type(iotk_unit), pointer :: this
  result = unit
  if(.not. iotk_units_init) then
    iotk_units_init = .true.
    call iotk_unit_list_init(iotk_units)
  end if
  call iotk_unit_get(unit,pointer=this)
  if(.not.associated(this)) return
  do
    if(.not. associated(this%son)) exit
    this => this%son
  end do
  result = this%unit
end function iotk_phys_unit_x

subroutine iotk_unit_print_x(unit)
  use iotk_base
  use iotk_str_interf
  implicit none 
  integer, intent(in) :: unit
!  type (iotk_unit), pointer :: this
!  stop 
!  this => iotk_units
!  write(unit,"(a)") "IOTK units"
!  do
!    if(.not. associated(this)) exit
!    write(unit,"(a,i8)") "Unit :",this%unit
!    write(unit,"(a,a,a,i8)") "Root :",this%root(1:iotk_strlen_trim(this%root)),"Level:",this%level
!    write(unit,"(a,l8)") "Raw  :",this%raw
!    if(associated(this%son)) then
!      write(unit,"(a,i8)") "Son :",this%son%unit
!    end if
!    if(associated(this%parent)) then
!      write(unit,"(a,i8)") "Parent :",this%parent%unit
!    end if
!!    this => this%next
!  end do
!  write(unit,"(a)") "end IOTK units"
end subroutine iotk_unit_print_x


subroutine iotk_unit_add_x(unit,this,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_unit_list_module
  implicit none
  integer,      intent(in)  :: unit
  type (iotk_unit), pointer :: this
  integer,      intent(out) :: ierr
  ierr = 0
  if(.not. iotk_units_init) then
    iotk_units_init = .true.
    call iotk_unit_list_init(iotk_units)
  end if
  call iotk_unit_list_search(iotk_units,this,unit=unit)
  if(associated(this)) then
    call iotk_error_issue(ierr,"iotk_unit_add",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
call iotk_error_msg(ierr,'unit')
    return
  end if
  call iotk_unit_list_add(iotk_units,this)
  this%unit         = unit
  this%root         = ""
  this%skip_root    = .false.
  this%raw          = .false.
!-<
  this%qe_syntax   = .false.
!->
  this%level        = 0
  this%close_at_end = .false.
  nullify(this%son)
  nullify(this%parent)
end subroutine iotk_unit_add_x

subroutine iotk_inquire_x(unit,binary,stream,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  integer, intent(in)  :: unit
  logical, intent(out) :: binary
  logical, intent(out) :: stream
  integer, intent(out) :: ierr
  character(50) :: form,access,pad,blank
  logical :: opened
  integer :: iostat
  iostat = 0
  ierr = 0
  inquire(unit=unit,form=form,iostat=iostat,access=access,pad=pad,blank=blank,opened=opened)
  if(iostat/=0) then
    call iotk_error_issue(ierr,"iotk_inquire",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
call iotk_error_msg(ierr,'Error inquiring')
    return
  end if
  if(opened .and. iotk_toupper(form)=="UNFORMATTED") then
    binary = .true.
  else
    binary = .false.
  end if
  stream = .false.
  if(opened) then
    select case(iotk_toupper(access))
    case("SEQUENTIAL")
    case("STREAM")
      stream = .true.
    case default
      call iotk_error_issue(ierr,"iotk_inquire",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
call iotk_error_msg(ierr,'Direct-access files are not allowed')
      return
    end select
  end if
  if(.not. binary) then
    if(opened .and. iotk_toupper(blank)/="NULL") then
      call iotk_error_issue(ierr,"iotk_inquire",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
      return
    end if
    if(opened .and. iotk_toupper(pad)  /="YES") then
      call iotk_error_issue(ierr,"iotk_inquire",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
      return
    end if
  end if
end subroutine iotk_inquire_x

subroutine iotk_unit_del_x(unit,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  use iotk_unit_list_module
  implicit none
  integer, intent(in)  :: unit
  integer, intent(out) :: ierr
  type (iotk_unit), pointer :: this
  ierr = 0
  if(.not. iotk_units_init) then
    iotk_units_init = .true.
    call iotk_unit_list_init(iotk_units)
  end if
  call iotk_unit_list_search(iotk_units,unit=unit,ptr=this)
  if(.not.associated(this)) then
    call iotk_error_issue(ierr,"iotk_unit_del",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
    return
  end if
  if(associated(this%parent)) nullify(this%parent%son)
  call iotk_unit_list_del(iotk_units,ptr=this)
end subroutine iotk_unit_del_x

subroutine iotk_unit_parent_x(parent,son,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_misc_interf
  use iotk_unit_interf
  implicit none
  integer, intent(in) :: parent,son
  integer, intent(out) :: ierr
  type(iotk_unit), pointer :: this_parent,this_son
  ierr = 0
  call iotk_unit_get(parent,pointer=this_parent)
  if(.not.associated(this_parent)) then
    call iotk_error_issue(ierr,"iotk_unit_parent",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
    return
  end if
  call iotk_unit_get(son,pointer=this_son)
  if(.not.associated(this_son)) then
    call iotk_error_issue(ierr,"iotk_unit_parent",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
    return
  end if
  if(associated(this_parent%son)) then
    call iotk_error_issue(ierr,"iotk_unit_parent",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
    return
  end if
  if(associated(this_son%parent)) then
    call iotk_error_issue(ierr,"iotk_unit_parent",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.16 ")
    return
  end if
  this_parent%son => this_son
  this_son%parent => this_parent
end subroutine iotk_unit_parent_x

subroutine iotk_unit_get_x(unit,pointer)
  use iotk_base
  use iotk_misc_interf
  use iotk_unit_list_module
  implicit none
  integer,                intent(in)  :: unit
  type(iotk_unit), pointer :: pointer
  if(.not. iotk_units_init) then
    iotk_units_init = .true.
    call iotk_unit_list_init(iotk_units)
  end if
  call iotk_unit_list_search(iotk_units,unit=unit,ptr=pointer)
end subroutine iotk_unit_get_x


!-<
! ... procedure that move the reading cursor at the beginning of the root node
subroutine iotk_rewind_x(unit,ierr) 
  use iotk_base
  use iotk_error_interf
  use iotk_scan_interf
  use iotk_unit_interf
  use iotk_misc_interf
  implicit none
  integer,                intent(in)  :: unit
  integer, optional,                intent(out) :: ierr
  integer                             :: control
  logical                             :: stream, binary
  character(iotk_taglenx) :: tag

  integer                 :: ierrl
  type(iotk_unit), pointer :: this

  ierrl = 0

  call iotk_unit_get(unit,pointer=this)
  call iotk_inquire(unit,binary,stream,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_rewind",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.16 ")
    goto 1
  end if

  if(.not.associated(this)) then
    call iotk_error_issue(ierrl,"iotk_rewind",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.16 ")
    goto 1
  end if

  rewind this%unit

  do
    call iotk_scan_tag(unit,+1,control,tag, binary,stream,ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_rewind",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.16 ")
      goto 1
    end if
    IF (control==1) EXIT
  end do
  this%level = 0

1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if

end subroutine iotk_rewind_x
!->
