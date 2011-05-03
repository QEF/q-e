!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! Adapted to TDDFPT by Osman Baris Malcioglu (2009)
!-----------------------------------------------------------------------
SUBROUTINE lr_dv_setup
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares some variable which is needed for derivatives
  !  1) Set non linear core correction stuff
  !  2) computes dmuxc 3) with GC if needed
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : ntyp => nsp
  USE lsda_mod,      ONLY : nspin, lsda
  USE scf,           ONLY : rho, rho_core
  USE gvect,         ONLY : ngm
  USE grid_dimensions,ONLY: nrxx
  ! USE atom,          ONLY : nlcc
  USE uspp_param,    ONLY : upf
  !USE lr_variables,  ONLY : dmuxc, nlcc_any
  USE nlcc_ph,       ONLY : nlcc_any
  USE eqv,           ONLY : dmuxc
  USE funct,         ONLY : dmxc, dmxc_spin
  USE lr_variables,  ONLY : lr_verbosity
  USE io_global,     ONLY : stdout
 IMPLICIT NONE
  !
  real(DP) :: rhotot, rhoup, rhodw
  ! total charge
  ! total up charge
  ! total down charge
  !
  INTEGER :: nt, ir
  ! counter on mesh points
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_dv_setup>")')
  ENDIF
  CALL start_clock ('lr_dv_setup')
  !
  ! 1) Set non linear core correction stuff
  !
  nlcc_any = any ( upf(1:ntyp)%nlcc )
  !do nt = 1, ntyp
  !   nlcc_any = nlcc_any.or.nlcc (nt)
  !enddo
  !
  ! 2) Computes the derivative of the xc potential
  !
  dmuxc(:,:,:) = 0.0D0
  IF (lsda) THEN
     DO ir = 1, nrxx
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        CALL dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     ENDDO
  ELSE
     DO ir = 1, nrxx
        rhotot = rho%of_r (ir, 1) + rho_core (ir)
        IF (rhotot>1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
        IF (rhotot<-1.d-30) dmuxc (ir, 1, 1) = -dmxc ( -rhotot)
     ENDDO
  ENDIF
  DEALLOCATE(rho_core)
  !
  ! 3) Setup all gradient correction stuff
  !
  CALL lr_setup_dgc()
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<end of lr_dv_setup>")')
  CALL stop_clock ('lr_dv_setup')
  !
  RETURN
END SUBROUTINE lr_dv_setup
