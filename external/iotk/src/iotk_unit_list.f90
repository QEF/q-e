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

module iotk_unit_list_module
use iotk_base
implicit none

type iotk_unit_list
  type (iotk_unit_list), pointer :: next
  type (iotk_unit),            pointer :: ptr
end type iotk_unit_list

logical,              save :: iotk_units_init = .false.
type(iotk_unit_list), save :: iotk_units

contains


subroutine iotk_unit_list_init(list)
  type (iotk_unit_list), intent(out) :: list
  nullify(list%ptr)
  nullify(list%next)
end subroutine iotk_unit_list_init

subroutine iotk_unit_list_destroy(list)
  type (iotk_unit_list), intent(inout) :: list
  type (iotk_unit_list), pointer       :: this,next
  if(.not.associated(list%next)) return
  this=>list%next
  do
    if(associated(this%ptr))deallocate(this%ptr)
    next=>this%next
    deallocate(this)
    if(.not.associated(next)) exit
    this=>next
  end do
end subroutine iotk_unit_list_destroy

subroutine iotk_unit_list_add(list,ptr)
  type (iotk_unit_list), intent(inout) :: list
  type (iotk_unit),            pointer       :: ptr
  type (iotk_unit_list), pointer       :: this
  allocate(this)
  this%next => list%next
  list%next => this
  allocate(this%ptr)
  ptr => this%ptr
end subroutine iotk_unit_list_add

subroutine iotk_unit_list_del(list,ptr)
  type (iotk_unit_list), intent(inout) :: list
  type (iotk_unit),            pointer       :: ptr
  type (iotk_unit_list), pointer       :: this,next_save
  if(.not.associated(list%next)) return
  if(associated(list%next%ptr,ptr)) then
    deallocate(list%next%ptr)
    next_save => list%next%next
    deallocate(list%next)
    list%next => next_save
    nullify(ptr)
    return
  end if
  this => list%next
  do
    if(.not.associated(this%next)) return
    if(associated(this%next%ptr,ptr)) exit
    this => this%next
  end do
  deallocate(this%next%ptr)
  next_save => this%next%next
  deallocate(this%next)
  this%next => next_save
  nullify(ptr)
end subroutine iotk_unit_list_del


 subroutine iotk_unit_list_search(list,ptr,unit)
   type (iotk_unit_list), intent(in) :: list
   type (iotk_unit),            pointer    :: ptr
 
  integer, optional,intent(in) :: unit

   type (iotk_unit_list), pointer    :: this
   nullify(ptr)
   this => list%next
   if(.not.associated(this)) return
   do
     if(.not.associated(this%ptr)) goto 1000


     if(present(unit)) then
       if(this%ptr%unit /= unit) goto 1000
     end if


     ptr => this%ptr
     exit
1000 continue
     if(.not.associated(this%next)) exit
     this => this%next
   end do
 end subroutine iotk_unit_list_search
 

end module iotk_unit_list_module
