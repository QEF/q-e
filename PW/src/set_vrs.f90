!
! Copyright (C) 2001-2024 Quantum ESPESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
SUBROUTINE set_vrs( vrs, vltot, vr, kedtau, kedtaur, nrxx, nspin, doublegrid )
  !-----------------------------------------------------------------------------
  !! Set the total local potential vrs on the smooth mesh to be used in 
  !! \(\texttt{h_psi}\), adding the (spin dependent) scf (H+xc) part and
  !! the sum of all the local pseudopotential contributions.
  !
  USE kinds
  USE xc_lib,          ONLY : xclib_dft_is
  USE fft_base,        ONLY : dffts, dfftp
  USE fft_interfaces,  ONLY : fft_interpolate
  !
  IMPLICIT NONE
  !
  INTEGER :: nspin
  !! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  INTEGER :: nrxx
  !! input: the fft grid dimension
  REAL(DP) :: vrs(nrxx,nspin)
  !! output: total local potential on the smooth grid vrs=vltot+vr
  REAL(DP) :: vltot(nrxx)
  !! input: the total local pseudopotential
  REAL(DP) :: vr(nrxx,nspin)
  !! input: the scf(H+xc) part of the local potential
  REAL(DP) :: kedtau(dffts%nnr,nspin)
  !! position dependent kinetic energy enhancement factor
  REAL(DP) :: kedtaur(nrxx,nspin)
  !! the kinetic energy density in R-space
  LOGICAL :: doublegrid
  ! input: true if a doublegrid is used
  INTEGER :: is
  !
  !$acc data present(vrs)
  DO is = 1, nspin
     !
     ! define the total local potential (external + scf) for each spin ...
     !
     IF (is > 1 .AND. nspin == 4) THEN
        !
        ! noncolinear case: only the first component contains vltot
        !
        vrs(:,is) = vr(:,is)
     ELSE
        vrs(:,is) = vltot(:) + vr(:,is)
     ENDIF
     !
  ENDDO
  !
  ! ... interpolate it on the smooth mesh if necessary
  !
  DO is = 1, nspin
     IF (doublegrid) CALL fft_interpolate( dfftp, vrs(:, is), dffts, vrs(:, is) )
     IF (xclib_dft_is('meta')) CALL fft_interpolate( dfftp, kedtaur(:,is), dffts, kedtau(:,is) )
  ENDDO
  !
  !$acc update device(vrs)
  !
  !$acc end data
  ! 
  RETURN
  !
END SUBROUTINE set_vrs
