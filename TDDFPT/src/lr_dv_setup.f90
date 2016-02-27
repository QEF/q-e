!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_dv_setup
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares some variables which are needed for derivatives
  !  1) Set the nonlinear core correction (nlcc) stuff 
  !  2) Computes dmuxc (derivative of the XC potential)
  !  3) Gradient correction staff (GGA)
  !
  ! Adapted to TDDFPT by Osman Baris Malcioglu (2009)
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : ntyp => nsp
  USE lsda_mod,              ONLY : nspin, lsda
  USE scf,                   ONLY : rho, rho_core
  USE fft_base,              ONLY : dfftp
  USE gvect,                 ONLY : ngm
  USE uspp_param,            ONLY : upf
  USE spin_orb,              ONLY : domag
  USE nlcc_ph,               ONLY : drc
  USE uspp,                  ONLY : nlcc_any
  USE noncollin_module,      ONLY : noncolin, nspin_mag
  USE eqv,                   ONLY : dmuxc
  USE funct,                 ONLY : dmxc, dmxc_spin, dmxc_nc
  USE lr_variables,          ONLY : lr_verbosity, lr_exx
  USE io_global,             ONLY : stdout
  USE funct,                 ONLY : dft_is_gradient, exx_is_active
  USE wavefunctions_module,  ONLY : psic
  !
  IMPLICIT NONE
  !
  REAL(DP) :: rhotot, rhoup, rhodw
  ! total charge
  ! total up charge
  ! total down charge
  REAL(DP) :: auxdmuxc(4,4)
  INTEGER :: nt, ir, is, js
  ! counter on mesh points
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_dv_setup>")')
  ENDIF
  !
  CALL start_clock ('lr_dv_setup')
  !
  ! 1) Set the nonlinear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  IF (nlcc_any) ALLOCATE(drc( ngm, ntyp))
  !
  ! 2) Computes the derivative of the XC potential
  !
  IF ( ( .not. exx_is_active() ) .AND. lr_exx ) THEN
     !
     ! If exx has been enforced from input but we are not using 
     ! any exact exchange (i.e. we are using the scissor approximation)
     ! the derivative of the XC potential is not used.
     !
     dmuxc(:,:,:) = 0.0_DP
     !
  ELSE
     !
     dmuxc(:,:,:) = 0.0_DP
     !
     IF (lsda) THEN
        !
        DO ir = 1, dfftp%nnr
           rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
           rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
           CALL dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                         dmuxc(ir,1,2), dmuxc(ir,2,2) )
        ENDDO
        !
     ELSE
        !
        IF (noncolin.and.domag) THEN
           DO ir = 1, dfftp%nnr
              rhotot = rho%of_r (ir, 1) + rho_core (ir)
              CALL dmxc_nc (rhotot, rho%of_r(ir,2), rho%of_r(ir,3), rho%of_r(ir,4), auxdmuxc)
              DO is=1,nspin_mag
                 DO js=1,nspin_mag
                    dmuxc(ir,is,js) = auxdmuxc(is,js)
                 ENDDO
              ENDDO
           ENDDO
        ELSE
           DO ir = 1, dfftp%nnr
              rhotot = rho%of_r (ir, 1) + rho_core (ir)
              IF (rhotot.GT.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
              IF (rhotot.LT. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
           ENDDO
        ENDIF
        !
     ENDIF
     !
  ENDIF
  !
  ! 3) Setup all gradient correction stuff
  !
  IF (dft_is_gradient()) THEN
     !
     IF (noncolin .AND. domag) THEN
        IF (.NOT.ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
        psic(:) = (0.0_dp, 0.0_dp)
     ENDIF
     !
     CALL setup_dgc()
     !
     IF (ALLOCATED(psic)) DEALLOCATE(psic)
     !
  ENDIF
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<end of lr_dv_setup>")')
  !
  CALL stop_clock ('lr_dv_setup')
  !
  RETURN
  !
END SUBROUTINE lr_dv_setup
