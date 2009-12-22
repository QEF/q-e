!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine symrho (rho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
  !-----------------------------------------------------------------------
  !
  !     Symmetrization of the charge density - Wrapper routine,
  !     calls the new reciprocal-space algorithm
  !
  ! nrx*, nr*: physical and actual dimensions of the FFT mesh
  ! nsym     : number of symmetry operations
  ! s, ftau  : symmetry operations (rotation + fractionary translation)
  !
  USE kinds
  USE gvect, ONLY : ngm, g, ngl, igtongl, nl, nrxx
  USE wavefunctions_module, ONLY : psic
  USE symme, ONLY : sym_rho
  !
  IMPLICIT NONE
  !
  integer, intent (IN) :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, &
                          s (3, 3, 48), ftau (3, 48)
  real(DP), intent (INOUT)  :: rho (nrxx)
  !
  complex(DP), allocatable :: rhog(:)
  !
  !
  if (nsym == 1) return
  psic(:) = rho(:)
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
  ALLOCATE ( rhog (ngm) )
  rhog(:) = psic(nl(:))
  !
  CALL sym_rho ( 1, rhog )
  !
  psic(:) = (0.0_dp, 0.0_dp)
  psic(nl(:)) = rhog(:)
  DEALLOCATE (rhog)
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  rho(:) = psic(:)
  !
  RETURN
END SUBROUTINE symrho

