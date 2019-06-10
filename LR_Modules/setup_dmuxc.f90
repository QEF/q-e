!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE setup_dmuxc
  !-----------------------------------------------------------------------
  !! This subroutine computes dmuxc (derivative of the XC potential).
  !
  USE kinds,            ONLY : DP
  USE eqv,              ONLY : dmuxc
  USE lsda_mod,         ONLY : lsda
  USE fft_base,         ONLY : dfftp
  USE scf,              ONLY : rho, rho_core
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE spin_orb,         ONLY : domag
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: rho_aux
  ! auxiliary array for density
  INTEGER  :: ir, is, js, ns
  !
  CALL start_clock ('setup_dmuxc')
  !
  ns = 1
  IF ( lsda ) ns = 2
  IF ( (.NOT. lsda) .AND. noncolin .AND. domag ) ns = 4
  !
  ALLOCATE( rho_aux(dfftp%nnr,ns) )
  !
  dmuxc(:,:,:) = 0.d0
  !
  IF ( lsda ) THEN
     !
     rho_aux(:,1) = ( rho%of_r(:,1) + rho%of_r(:,2) + rho_core(:) )*0.5_DP
     rho_aux(:,2) = ( rho%of_r(:,1) - rho%of_r(:,2) + rho_core(:) )*0.5_DP
     !
     CALL dmxc( dfftp%nnr, 2, rho_aux, dmuxc )
     !
  ELSE
     !
     IF ( noncolin .AND. domag ) THEN
        !
        rho_aux(:,1) = rho%of_r(:,1) + rho_core(:)
        rho_aux(:,2:4) = rho%of_r(:,2:4)
        CALL dmxc( dfftp%nnr, 4, rho_aux, dmuxc )
        !
     ELSE
        !
        rho_aux(:,1) = rho%of_r(:,1) + rho_core(:)
        CALL dmxc( dfftp%nnr, 1, rho_aux, dmuxc )
        !
     ENDIF
     !
  ENDIF
  !
  DEALLOCATE( rho_aux )
  !
  CALL stop_clock ('setup_dmuxc')
  !
  RETURN
  !
END SUBROUTINE setup_dmuxc
