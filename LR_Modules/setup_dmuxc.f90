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
  !! This subroutine computes dmuxc (derivative of the XC potential)
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
  !! auxiliary array for density
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: sign_r
  !! auxiliary array to put the correct sign on dmuxc when not lsda
  INTEGER  :: ir, is, js, ns
  !
  CALL start_clock ('setup_dmuxc')
  !
  ns = 1
  IF ( lsda ) ns = 2
  ALLOCATE(rho_aux(dfftp%nnr,ns))
  !
  dmuxc(:,:,:) = 0.d0
  !
  IF ( lsda ) THEN
      !
      rho_aux(:,1) = ( rho%of_r(:,1) + rho%of_r(:,2) + rho_core(:) ) * 0.5d0
      rho_aux(:,2) = ( rho%of_r(:,1) - rho%of_r(:,2) + rho_core(:) ) * 0.5d0
      CALL dmxc_spin(dfftp%nnr, rho_aux, dmuxc )
      !
  ELSE
      !
      IF ( noncolin .AND. domag ) THEN
         !
         rho_aux(:,1) = rho%of_r (:,1) + rho_core(:)
         CALL dmxc_nc( dfftp%nnr, rho_aux(:,1), rho%of_r(:,2:4), dmuxc )
         !
      ELSE
         !
         ALLOCATE(sign_r(dfftp%nnr))
         rho_aux(:,1) = rho%of_r(:,1) + rho_core(:)
         sign_r = 1.0_DP
         DO ir = 1, dfftp%nnr
            IF ( rho_aux(ir,1) < -1.d-30 ) THEN
               sign_r(ir) = -1.0_DP
               rho_aux(ir,1) = -rho_aux(ir,1)
            ELSEIF ( rho_aux(ir,1) < 1.d-30 .AND.rho_aux(ir,1) > -1.d-30 ) THEN
               sign_r(ir) = 0.0_DP
               rho_aux(ir,1) = 0.5_DP
            ENDIF
         ENDDO
         !
         CALL dmxc( dfftp%nnr, rho_aux(:,1), dmuxc(:,1,1) )
         dmuxc(:,1,1) = dmuxc(:,1,1)*sign_r(:)
         !
         DEALLOCATE(sign_r)
         !
      ENDIF
      !
  ENDIF
  !
  DEALLOCATE(rho_aux)
  !
  CALL stop_clock ('setup_dmuxc')
  !
  RETURN
  !
END SUBROUTINE setup_dmuxc
