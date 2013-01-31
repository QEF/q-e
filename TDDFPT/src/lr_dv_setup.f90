!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
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
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : ntyp => nsp
  USE lsda_mod,      ONLY : nspin, lsda
  USE scf,           ONLY : rho, rho_core
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm
  ! USE atom,          ONLY : nlcc
  USE uspp_param,    ONLY : upf
  USE spin_orb,      ONLY : domag
  !USE lr_variables,  ONLY : dmuxc, nlcc_any
  USE nlcc_ph,       ONLY : drc,nlcc_any
  USE noncollin_module, ONLY : noncolin,nspin_mag
  USE eqv,           ONLY : dmuxc
  USE funct,         ONLY : dmxc, dmxc_spin
  USE lr_variables,  ONLY : lr_verbosity, lr_exx
  USE io_global,     ONLY : stdout
  USE funct,         ONLY : exx_is_active

 IMPLICIT NONE
  !
  real(DP) :: rhotot, rhoup, rhodw
  ! total charge
  ! total up charge
  ! total down charge
  !
  real(DP) :: auxdmuxc(4,4)

  !
  INTEGER :: nt, ir,is,js
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
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  if (nlcc_any) allocate (drc( ngm, ntyp))

  !
  ! 2) Computes the derivative of the xc potential
  !
  !
  IF ( ( .not. exx_is_active() ) .AND. lr_exx ) THEN
     ! If exx has been enforced from input but we are not using 
     ! any exact exchange (ie we are using the scissor approximation)
     ! the derivative of the xc potential is not used.
     dmuxc(:,:,:) = 0.0_DP
  ELSE
     dmuxc(:,:,:) = 0.0_DP
     IF (lsda) THEN
        DO ir = 1, dfftp%nnr
           rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
           rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
           CALL dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                dmuxc(ir,1,2), dmuxc(ir,2,2) )
        ENDDO
     ELSE
!    IF (noncolin.and.domag) THEN
!       do ir = 1, dfftp%nnr
!          rhotot = rho%of_r (ir, 1) + rho_core (ir)
!          call dmxc_nc (rhotot, rho%of_r(ir,2), rho%of_r(ir,3), rho%of_r(ir,4), auxdmuxc)
!          DO is=1,nspin_mag
!             DO js=1,nspin_mag
!                dmuxc(ir,is,js)=auxdmuxc(is,js)
!             END DO
!          END DO
!       enddo
!    ELSE
        DO ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           IF (rhotot.GT.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
           IF (rhotot.LT. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
        ENDDO
!     END IF
     ENDIF
  ENDIF
  !
  ! 3.1) Setup all gradient correction stuff
  !
  call setup_dgc
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<end of lr_dv_setup>")')
  CALL stop_clock ('lr_dv_setup')
  !
  RETURN
END SUBROUTINE lr_dv_setup
