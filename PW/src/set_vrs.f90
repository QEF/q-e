!
! Copyright (C) 2001 PWSCF group
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
  USE funct,    ONLY : dft_is_meta
  USE fft_base, ONLY : dffts 
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
  !
  CALL sum_vrs( nrxx, nspin, vltot, vr, vrs )
  !
  CALL interpolate_vrs( nrxx, nspin, doublegrid, kedtau, kedtaur, vrs )
  ! 
  RETURN
  !
END SUBROUTINE set_vrs
!
!
!--------------------------------------------------------------------
SUBROUTINE sum_vrs( nrxx, nspin, vltot, vr, vrs )
  !--------------------------------------------------------------------
  !! Accumulates local potential contributions into vrs (the total local potential).
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER :: nspin
  !! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  INTEGER :: nrxx
  !! input: the fft grid dimension
  REAL(DP) :: vrs(nrxx,nspin)
  !! output: total local potential on the smooth grid:  
  !! \(\text{vrs\}=\text{vltot}+\text{vr}\)
  REAL(DP) :: vltot(nrxx)
  !! input: the total local pseudopotential
  REAL(DP) :: vr(nrxx,nspin)
  !! input: the scf(H+xc) part of the local potential
  !
  ! ... local variable
  !
  INTEGER :: is
  !
  !
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
  RETURN
  !
END SUBROUTINE sum_vrs
!
!--------------------------------------------------------------------------
SUBROUTINE interpolate_vrs( nrxx, nspin, doublegrid, kedtau, kedtaur, vrs )
  !--------------------------------------------------------------------------
  !! Interpolates local potential on the smooth mesh if necessary.
  !
  USE kinds
  USE funct,           ONLY : dft_is_meta
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
  !! output: total local potential interpolated on the smooth grid
  REAL(DP) :: kedtau(dffts%nnr,nspin)
  !! position dependent kinetic energy enhancement factor
  REAL(DP) :: kedtaur(nrxx,nspin)
  !! the kinetic energy density in R-space
  LOGICAL :: doublegrid
  !! input: true if a doublegrid is used
  !
  ! ... local variable
  !
  INTEGER :: is
  !
  ! ... interpolate it on the smooth mesh if necessary
  !
  DO is = 1, nspin
     IF (doublegrid) CALL fft_interpolate( dfftp, vrs(:, is), dffts, vrs(:, is) )
     IF (dft_is_meta()) CALL fft_interpolate( dfftp, kedtaur(:,is), dffts, kedtau(:,is) )
  ENDDO
  !
  RETURN
  !
END SUBROUTINE interpolate_vrs
