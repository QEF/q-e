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
  !
  ! This subroutine computes dmuxc (derivative of the XC potential)
  !
  USE kinds,            ONLY : DP
  USE eqv,              ONLY : dmuxc
  USE lsda_mod,         ONLY : lsda
  USE fft_base,         ONLY : dfftp
  USE scf,              ONLY : rho, rho_core
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : dmxc, dmxc_spin, dmxc_nc
  !
  IMPLICIT NONE
  !
  REAL(DP) :: rhotot, rhoup, rhodw
  ! total charge
  ! total up charge
  ! total down charge
  REAL(DP) :: auxdmuxc(4,4)
  INTEGER  :: ir, is, js
  !
  CALL start_clock ('setup_dmuxc')
  !
  dmuxc(:,:,:) = 0.d0
  !
  IF (lsda) THEN
     DO ir = 1, dfftp%nnr
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        CALL dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     ENDDO
  ELSE
     IF (noncolin.and.domag) THEN
        DO ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           CALL dmxc_nc (rhotot, rho%of_r(ir,2), rho%of_r(ir,3), rho%of_r(ir,4), auxdmuxc)
           DO is=1,nspin_mag
              DO js=1,nspin_mag
                 dmuxc(ir,is,js)=auxdmuxc(is,js)
              ENDDO
           ENDDO
        ENDDO
     ELSE
        DO ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           IF (rhotot.GT.1.d-30)    dmuxc (ir, 1, 1) =   dmxc (rhotot)
           IF (rhotot.LT. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
        ENDDO
     ENDIF
  ENDIF
  !
  CALL stop_clock ('setup_dmuxc')
  !
  RETURN
  !
END SUBROUTINE setup_dmuxc
