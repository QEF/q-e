!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_config_to_file
  !-----------------------------------------------------------------------
#if defined (__OLDPUNCH)
  USE control_flags, ONLY : lscf
  USE kinds, ONLY : DP
  USE io_files, ONLY : prefix, iunres
  USE restart_module, ONLY : writefile_new
  !
  implicit none
  !
  logical :: exst
  integer :: iunit
  integer :: kunittmp
  real(DP) :: et_g(1,1), wg_g(1,1)

  !
  ! do not modify the file if in a non-scf run..
  ! probably not needed precaution
  !
  if ( lscf ) then
    !
    ! write configuration information directly to the general pourpose
    ! ".save" file
    !
    kunittmp = 1
    call writefile_new( 'config', iunres, et_g, wg_g, kunittmp )
  end if
  !
#endif
  return
end subroutine write_config_to_file
