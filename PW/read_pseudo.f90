!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine read_pseudo (is, iunps, ierr)
  !---------------------------------------------------------------------
  !
  !   read "is"-th pseudopotential in the Unified Pseudopotential Format
  !   from unit "iunps" - convert and copy to internal PWscf variables
  !   return error code in "ierr" (success: ierr=0)
  !
  use pseudo_types
  use read_pseudo_module
  use upf_to_internal
  !
  implicit none
  !
  integer :: is, iunps, ierr
  !
  TYPE (pseudo_upf) :: upf
  !
  !
  call read_pseudo_upf(iunps, upf, ierr)
  !
  if (ierr == 0) then
     call set_pseudo_upf (is, upf)
  end if
  !
  CALL deallocate_pseudo_upf( upf )
  !
end subroutine read_pseudo
