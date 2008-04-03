!
! Copyright (C) 2004-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module is USEd, for the time being, as an interface
! between the UPF pseudo type and the pseudo variables internal representation

!=----------------------------------------------------------------------------=!
  MODULE upf_to_internal
!=----------------------------------------------------------------------------=!

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: set_pseudo_upf
  SAVE

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!
!
!---------------------------------------------------------------------
subroutine set_pseudo_upf (is, upf, grid)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal variables
  !
  USE funct, ONLY: set_dft_from_name, set_dft_from_indices, dft_is_meta
  !
  USE pseudo_types
  USE radial_grids, ONLY: radial_grid_type, allocate_radial_grid
  !
  implicit none
  !
  integer :: is
  !
  !     Local variables
  !
  integer :: nb, mb, ijv
  integer :: iexch,icorr,igcx,igcc
  TYPE(radial_grid_type),target,optional :: grid ! if present reconstruct radial grid.
                              ! (only for old format pseudos)
  TYPE (pseudo_upf) :: upf
  !
  !
  ! workaround for rrkj format - it contains the indices, not the name
  if ( upf%dft(1:6)=='INDEX:') then
     read( upf%dft(7:10), '(4i1)') iexch,icorr,igcx,igcc
     call set_dft_from_indices(iexch,icorr,igcx,igcc)
  else
     call set_dft_from_name( upf%dft )
  end if
  !
  if(present(grid)) then
    call allocate_radial_grid(grid,upf%mesh)
    grid%dx   = upf%dx
    grid%xmin = upf%xmin
    grid%zmesh= upf%zmesh
    grid%mesh = upf%mesh
    !
    grid%r  (1:upf%mesh) = upf%r  (1:upf%mesh)
    grid%rab(1:upf%mesh) = upf%rab(1:upf%mesh)
    upf%grid => grid
  endif
  !
end subroutine set_pseudo_upf


!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
