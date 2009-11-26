!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! OBM
!  150608 (location of nlcc changed, new method of setting nlcc_any)
!         rho -> rho%of_r

!-----------------------------------------------------------------------
subroutine lr_dv_setup
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
  USE gvect,         ONLY : nrxx, ngm
  ! USE atom,          ONLY : nlcc
  USE uspp_param,           ONLY : upf
  USE lr_variables,  ONLY : dmuxc, nlcc_any
  USE funct,         ONLY : dmxc, dmxc_spin  
  USE lr_variables,   ONLY : lr_verbosity
   USE io_global,      ONLY : stdout
 implicit none
  !
  real(DP) :: rhotot, rhoup, rhodw
  ! total charge
  ! total up charge
  ! total down charge
  !
  integer :: nt, ir
  ! counter on mesh points
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_dv_setup>")')
  endif
  call start_clock ('lr_dv_setup')
  !
  ! 1) Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !do nt = 1, ntyp
  !   nlcc_any = nlcc_any.or.nlcc (nt)
  !enddo
  !
  ! 2) Computes the derivative of the xc potential
  !
  dmuxc(:,:,:) = 0.0D0
  if (lsda) then
     do ir = 1, nrxx
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        call dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     enddo
  else
     do ir = 1, nrxx
        rhotot = rho%of_r (ir, 1) + rho_core (ir)
        if (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
        if (rhotot.lt.-1.d-30) dmuxc (ir, 1, 1) = -dmxc ( -rhotot)
     enddo
  endif
  !
  ! 3) Setup all gradient correction stuff
  !
  call lr_setup_dgc()
  !
  if (lr_verbosity > 5) WRITE(stdout,'("<end of lr_dv_setup>")') 
  call stop_clock ('lr_dv_setup')
  !
  return
end subroutine lr_dv_setup
