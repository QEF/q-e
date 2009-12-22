!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine psymrho_mag (rho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, &
     ftau, bg, at)
  !-----------------------------------------------------------------------
  !
  !     Symmetrization of the magnetization density - Wrapper routine,
  !     calls the new reciprocal-space parallel algorithm
  !
  ! nrx*, nr*: physical and actual dimensions of the FFT mesh
  ! nsym     : number of symmetry operations
  ! s, ftau  : symmetry operations (rotation + fractionary translation)
  !
  USE kinds
  USE symme, ONLY : sym_rho
  USE gvect, ONLY : ngm, nl, nrxx
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE
  !
  integer, intent (IN) :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, &
                          s (3, 3, 48), ftau (3, 48)
  REAL(DP), INTENT (IN) :: bg(3,3), at(3,3)
  real(DP), intent (INOUT)  :: rho (nrxx,3)
  !
  integer :: is, nspin=4
  complex(DP), allocatable :: rhog(:,:)
  !
  !
  if (nsym == 1) return

  ALLOCATE ( rhog (ngm,nspin) )
  rhog(:,1) = (0.0_dp, 0.0_dp)
  DO is=2,nspin
     psic(:) = rho(:, is-1)
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
     rhog(:,is) = psic(nl(:))
  END DO
  !
  CALL sym_rho ( nspin, rhog )
  !
  DO is=2,nspin
     psic(:) = (0.0_dp, 0.0_dp)
     psic(nl(:)) = rhog(:,is)
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
     rho(:,is-1) = psic(:)
  END DO
  DEALLOCATE (rhog)
  RETURN
end subroutine psymrho_mag

