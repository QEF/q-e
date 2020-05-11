!
! Copyright (C) 2020-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! This soubroutine has been adapted from iotk_free_unit 
! originally included in the IOTK library.
!
!
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
!
! Minimum value used in iotk_free_unit
#ifndef __IOTK_UNITMIN
#  define __IOTK_UNITMIN 90000
#endif
!
! Maximum value used in iotk_free_unit
#ifndef __IOTK_UNITMAX
#  define __IOTK_UNITMAX 99999
#endif
!
subroutine upf_free_unit(unit,ierr)
  implicit none
! This subroutine sets 'unit' to the number of
! an I/O unit which is free (i.e. not already opened).
! The search is carried out starting from unit
! 'unitmin' in a range of 'nsearch' units.
! The starting unit for the search is increased at each
! call, so that a number of subsequent ask can be done
! obtaining different units.
  integer, parameter :: iotk_unitmin = __IOTK_UNITMIN 
  integer, parameter :: iotk_unitmax = __IOTK_UNITMAX
  !
  integer,           intent(out) :: unit
  integer, optional, intent(out) :: ierr
  integer, save :: offset = 0
  logical       :: opened,exist
  integer       :: isearch,nsearch,unitmin
  integer       :: ierrl
  integer       :: iostat
  character(256):: msg
  !
  iostat = 0
  unitmin = iotk_unitmin
  nsearch = iotk_unitmax - iotk_unitmin + 1
  ierrl = 0 
  do isearch=0,nsearch-1
    unit = modulo(isearch+offset,nsearch) + unitmin
    inquire(unit=unit,opened=opened,exist=exist,iostat=iostat)
    if(iostat/=0) then
      ierrl=10
      msg='Error inquiring'
      goto 1
    endif
    if((.not.opened .and. exist) .or. iostat/=0) exit
  end do
  if(isearch>=nsearch) then
    ierrl=20
    msg='There are no units left'
    goto 1
  end if 
  offset = modulo(unit - unitmin + 1,nsearch)
1 continue
  if(present(ierr)) then 
    ierr = ierrl
  else
    if(ierrl/=0) call upf_error("upf_free_unit",msg,ierrl)
  end if
end subroutine upf_free_unit

