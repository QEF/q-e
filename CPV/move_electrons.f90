!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE move_electrons( nfi, tfirst, tlast, b1, b2, b3, fion, &
                           enthal, enb, enbi, fccc, ccc, dt2bye )
  !----------------------------------------------------------------------------
  !
  ! ... this routine updates the electronic degrees of freedom
  !
  USE kinds,                ONLY : DP
  USE parameters,           ONLY : natx
  USE control_flags,        ONLY : lwf, trhow, tfor, tprnfor, thdyn
  USE cg_module,            ONLY : tcg
  USE cp_main_variables,    ONLY : eigr, bec, irb, eigrb, rhog, rhos, rhor, &
                                   rhopr, ei1, ei2, ei3, sfac, ema0bg, becdr, &
                                   taub, lambda, lambdam, lambdap
  USE wavefunctions_module, ONLY : c0, cm, phi => cp
  USE cell_base,            ONLY : omega, ibrav, h, press
  USE uspp,                 ONLY : becsum, vkb, nkb
  USE energies,             ONLY : ekin, enl, entropy, etot
  USE grid_dimensions,      ONLY : nnrx
  USE electrons_base,       ONLY : nbsp, nspin, f
  USE core,                 ONLY : nlcc_any, rhoc
  USE ions_positions,       ONLY : tau0
  USE stre,                 ONLY : stress
  USE dener,                ONLY : detot
  USE efield_module,        ONLY : tefield, ipolp, qmat, gqq, evalue
  !
  USE wannier_subroutines,  ONLY : get_wannier_center, wf_options, &
                                   write_charge_and_exit, ef_tune
  USE cpr_subroutines,      ONLY : compute_stress
  USE ensemble_dft,         ONLY : compute_entropy2
  USE efield_module,        ONLY : berry_energy
  USE runcp_module,         ONLY : runcp_uspp
  USE wave_constrains,      ONLY : interpolate_lambda
  USE gvecw,                ONLY : ngw
  USE orthogonalize_base,   ONLY : calphi
  !
  IMPLICIT NONE
  !
  INTEGER,        INTENT(IN)    :: nfi
  LOGICAL,        INTENT(IN)    :: tfirst, tlast
  REAL(DP), INTENT(IN)    :: b1(3), b2(3), b3(3)
  REAL(DP), INTENT(IN)    :: fion(3,natx)
  REAL(DP), INTENT(IN)    :: dt2bye
  REAL(DP), INTENT(IN)    :: fccc, ccc
  REAL(DP), INTENT(INOUT) :: enb, enbi
  REAL(DP), INTENT(INOUT) :: enthal
  !
  INTEGER :: i, is
  !
  !
  electron_dynamic: IF ( tcg ) THEN
     !
     CALL runcg_uspp( nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, &
                      fion, ema0bg, becdr, lambdap, lambda  )
     !
  ELSE
     !
     IF ( lwf ) &
          CALL get_wannier_center( tfirst, cm, bec, becdr, eigr, &
                                   eigrb, taub, irb, ibrav, b1, b2, b3 )
     !
     CALL rhoofr( nfi, c0, irb, eigrb, bec, &
                     becsum, rhor, rhog, rhos, enl, ekin )
        !
     !
     IF ( tfirst .OR. tlast ) rhopr = rhor
     !
     ! ... put core charge (if present) in rhoc(r)
     !
     IF ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc )
     !
     IF ( lwf ) THEN
        !
        CALL write_charge_and_exit( rhog )
        CALL ef_tune( rhog, tau0 )
        !
     END IF
     !
     CALL vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )
     !
     IF ( lwf ) CALL wf_options( tfirst, nfi, cm, becsum, bec, becdr, &
                                 eigr, eigrb, taub, irb, ibrav, b1,   &
                                 b2, b3, rhor, rhog, rhos, enl, ekin  )
     !
     CALL compute_stress( stress, detot, h, omega )
     !
     enthal = etot + press * omega
     !
     IF( tefield )  THEN
        !
        CALL berry_energy( enb, enbi, bec, c0(:,:,1,1), fion )
        !
        etot = etot + enb + enbi
        !
     END IF
     !
     !=======================================================================
     !
     !              verlet algorithm
     !
     !     loop which updates electronic degrees of freedom
     !     cm=c(t+dt) is obtained from cm=c(t-dt) and c0=c(t)
     !     the electron mass rises with g**2
     !
     !=======================================================================
     !
     CALL newd( rhor, irb, eigrb, becsum, fion )
     !
     CALL prefor( eigr, vkb )
     !
     CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, &
                      rhos, bec, c0(:,:,1,1), cm(:,:,1,1) )
     !
     !----------------------------------------------------------------------
     !                 contribution to fion due to lambda
     !----------------------------------------------------------------------
     !
     ! ... nlfq needs deeq bec
     !
     IF ( tfor .OR. tprnfor ) CALL nlfq( c0, eigr, bec, becdr, fion )
     !
     IF ( (tfor.or.tprnfor) .AND. tefield ) &
        CALL bforceion( fion, .TRUE. , ipolp, qmat, bec, becdr, gqq, evalue )
     !
     IF ( tfor .OR. thdyn ) then
        CALL interpolate_lambda( lambdap, lambda, lambdam )
     ELSE
        ! take care of the otherwise uninitialized lambdam
        lambdam(:,:)=lambda(:,:)
     END IF
     !
     ! ... calphi calculates phi
     ! ... the electron mass rises with g**2
     !
     CALL calphi( c0, ngw, bec, nkb, vkb, phi, nbsp, ema0bg )
     !
     ! ... begin try and error loop (only one step!)
     !
     ! ... nlfl and nlfh need: lambda (guessed) becdr
     !
     IF ( tfor .OR. tprnfor ) CALL nlfl( bec, becdr, lambda, fion )
     !
  END IF electron_dynamic
  !
  RETURN
  !
END SUBROUTINE move_electrons
