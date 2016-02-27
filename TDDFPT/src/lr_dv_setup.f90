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
  !  1) Set the nonlinear core correction 
  !  2) Computes dmuxc (derivative of the XC potential)
  !  3) Set gradient correction (GGA) if needed
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : ntyp => nsp
  USE fft_base,              ONLY : dfftp
  USE uspp_param,            ONLY : upf
  USE spin_orb,              ONLY : domag
  USE uspp,                  ONLY : nlcc_any
  USE noncollin_module,      ONLY : noncolin
  USE eqv,                   ONLY : dmuxc
  USE lr_variables,          ONLY : lr_exx
  USE funct,                 ONLY : dft_is_gradient, exx_is_active
  USE wavefunctions_module,  ONLY : psic
  !
  IMPLICIT NONE
  !
  CALL start_clock ('lr_dv_setup')
  !
  ! 1) Set the nonlinear core correction
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  ! 2) Compute the derivative of the XC potential
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
     CALL setup_dmuxc()
     !
  ENDIF
  !
  ! 3) Setup gradient correction
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
  CALL stop_clock ('lr_dv_setup')
  !
  RETURN
  !
END SUBROUTINE lr_dv_setup
