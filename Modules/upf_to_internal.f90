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
subroutine set_pseudo_upf (is, upf)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal variables
  !
  USE atom,  ONLY: rgrid, rho_at, rho_atc, nlcc
  USE uspp_param, ONLY: tvanp
  USE funct, ONLY: set_dft_from_name, set_dft_from_indices, dft_is_meta
  !
  USE pseudo_types
  !
  implicit none
  !
  integer :: is
  !
  !     Local variables
  !
  integer :: nb, mb, ijv
  integer :: iexch,icorr,igcx,igcc
  TYPE (pseudo_upf) :: upf
  !
  !
  tvanp(is)=upf%tvanp
  nlcc(is) = upf%nlcc
  ! workaround for rrkj format - it contains the indices, not the name
  if ( upf%dft(1:6)=='INDEX:') then
     read( upf%dft(7:10), '(4i1)') iexch,icorr,igcx,igcc
     call set_dft_from_indices(iexch,icorr,igcx,igcc)
  else
     call set_dft_from_name( upf%dft )
  end if
  !
  rgrid(is)%dx   = upf%dx
  rgrid(is)%xmin = upf%xmin
  rgrid(is)%zmesh= upf%zmesh
  rgrid(is)%mesh = upf%mesh
  IF ( rgrid(is)%mesh > SIZE (rgrid(is)%r) ) &
     CALL errore('upf_to_internals', 'too many grid points', 1)
  !
  rgrid(is)%r  (1:upf%mesh) = upf%r  (1:upf%mesh)
  rgrid(is)%rab(1:upf%mesh) = upf%rab(1:upf%mesh)
  !
  if ( upf%nlcc) then
     rho_atc(1:upf%mesh, is) = upf%rho_atc(1:upf%mesh)
  else
     rho_atc(:,is) = 0.0_DP
  end if
  rho_at (1:upf%mesh, is) = upf%rho_at (1:upf%mesh)

end subroutine set_pseudo_upf


!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
