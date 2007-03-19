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
  USE parameters, ONLY: ndmx
  USE atom,  ONLY: mesh, r, rab, chi, oc, nchi, lchi, jchi, rho_at, &
                   rho_atc, nlcc
  USE uspp_param, ONLY: zp, vloc_at, dion, betar, qqq, qfcoef, qfunc, nqf, &
                        nqlc, rinner, nbeta, kkbeta, lll, jjj, psd, tvanp
  USE funct, ONLY: set_dft_from_name, dft_is_meta
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
  TYPE (pseudo_upf) :: upf
  !
  !
  zp(is)  = upf%zp
  psd (is)= upf%psd
  tvanp(is)=upf%tvanp
  nlcc(is) = upf%nlcc
  call set_dft_from_name( upf%dft )
  !
  mesh(is) = upf%mesh
  IF ( mesh(is) > ndmx ) &
     CALL errore('upf_to_internals', 'too many grid points', 1)
  !
  nchi(is) = upf%nwfc
  lchi(1:upf%nwfc, is) = upf%lchi(1:upf%nwfc)
  oc(1:upf%nwfc, is) = upf%oc(1:upf%nwfc)
  chi(1:upf%mesh, 1:upf%nwfc, is) = upf%chi(1:upf%mesh, 1:upf%nwfc)
  !
  nbeta(is)= upf%nbeta
  kkbeta(is)=0
  do nb=1,upf%nbeta
     kkbeta(is)=max(upf%kkbeta(nb),kkbeta(is))
  end do
  betar(1:upf%mesh, 1:upf%nbeta, is) = upf%beta(1:upf%mesh, 1:upf%nbeta)
  dion(1:upf%nbeta, 1:upf%nbeta, is) = upf%dion(1:upf%nbeta, 1:upf%nbeta)
  !
  nqlc(is) = upf%nqlc
  nqf (is) = upf%nqf
  lll(1:upf%nbeta,is) = upf%lll(1:upf%nbeta)
  rinner(1:upf%nqlc,is) = upf%rinner(1:upf%nqlc)
  qqq(1:upf%nbeta,1:upf%nbeta,is) = upf%qqq(1:upf%nbeta,1:upf%nbeta)
  do nb = 1, upf%nbeta
     do mb = nb, upf%nbeta
        ijv = mb * (mb-1) / 2 + nb
        qfunc (1:upf%mesh, ijv, is) = upf%qfunc(1:upf%mesh, nb, mb)
     end do
  end do
  qfcoef(1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta, is ) = &
       upf%qfcoef( 1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta )
  !
  r  (1:upf%mesh, is) = upf%r  (1:upf%mesh)
  rab(1:upf%mesh, is) = upf%rab(1:upf%mesh)

  if (upf%has_so) then
     jchi(1:upf%nwfc, is) = upf%jchi(1:upf%nwfc)
     jjj(1:upf%nbeta, is) = upf%jjj(1:upf%nbeta)
  else
     jchi(1:upf%nwfc, is) = 0.d0
     jjj(1:upf%nbeta, is) = 0.d0
  endif
  !
  if ( upf%nlcc) then
     rho_atc(1:upf%mesh, is) = upf%rho_atc(1:upf%mesh)
  else
     rho_atc(:,is) = 0.d0
  end if
  rho_at (1:upf%mesh, is) = upf%rho_at (1:upf%mesh)
  vloc_at(1:upf%mesh,is) = upf%vloc(1:upf%mesh)

end subroutine set_pseudo_upf

!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
