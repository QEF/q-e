!
! Copyright (C) 2002-2005 FPMD-CPV groups
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
  USE kinds,                ONLY : dbl
  USE parameters,           ONLY : natx
  USE control_flags,        ONLY : lwf, trhow, tfor, tprnfor, thdyn
  USE cg_module,            ONLY : tcg
  USE cp_main_variables,    ONLY : eigr, bec, irb, eigrb, rhog, rhos, rhor,  &
                                   ei1, ei2, ei3, sfac, ema0bg, becdr, &
                                   taub, lambda, lambdam, lambdap
  USE wavefunctions_module, ONLY : c0, cm, phi => cp
  USE ensemble_dft,         ONLY : tens, z0, c0diag, becdiag, bec0, v0s, &
                                   vhxcs, becdrdiag
  USE cell_base,            ONLY : omega, ibrav, h, press
  USE uspp,                 ONLY : becsum, vkb
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
  USE para_mod,             ONLY : 
  !
  IMPLICIT NONE
  !
  INTEGER,        INTENT(IN)    :: nfi
  LOGICAL,        INTENT(IN)    :: tfirst, tlast
  REAL(KIND=dbl), INTENT(IN)    :: b1(3), b2(3), b3(3)
  REAL(KIND=dbl), INTENT(IN)    :: fion(3,natx)
  REAL(KIND=dbl), INTENT(IN)    :: dt2bye
  REAL(KIND=dbl), INTENT(IN)    :: fccc, ccc
  REAL(KIND=dbl), INTENT(INOUT) :: enb, enbi
  REAL(KIND=dbl), INTENT(INOUT) :: enthal
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
     IF ( .NOT. tens ) THEN
        !
        CALL rhoofr( nfi, c0, irb, eigrb, bec, &
                     becsum, rhor, rhog, rhos, enl, ekin )
        !
     ELSE
        !
        ! ... calculation of the rotated quantities
        !
        CALL rotate( z0, c0(:,:,1,1), bec, c0diag, becdiag )
        !
        ! ... calculation of rho corresponding to the rotated wavefunctions
        !
        CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, &
                     becsum, rhor, rhog, rhos, enl, ekin )
        !
        bec0(:,:) = bec(:,:)
        !
     END IF
     !
#if defined (__PARA)
     IF ( trhow .AND. tlast ) CALL write_rho( 47, nspin, rhor )
#else
     IF ( trhow .AND. tlast ) &
          WRITE(47) ( ( rhor(i,is), i = 1, nnrx ), is = 1, nspin )
#endif
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
     IF ( .NOT. tens ) THEN
        !
        CALL vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, tlast, &
                     ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )
        !
     ELSE
        !
        IF ( tfirst ) CALL compute_entropy2( entropy, f, nbsp, nspin )
        !
        CALL vofrho2( nfi, rhor, rhog, rhos, rhoc, tfirst, tlast, ei1, &
                      ei2, ei3, irb, eigrb, sfac, tau0, fion, v0s, vhxcs )
     END IF
     !
     IF ( lwf ) CALL wf_options( tfirst, nfi, cm, becsum, bec, becdr, &
                                 eigr, eigrb, taub, irb, ibrav, b1,   &
                                 b2, b3, rhor, rhog, rhos, enl, ekin  )
     !
     CALL compute_stress( stress, detot, h, omega )
     !
     enthal = etot + press * omega
     !
     IF( tefield .AND. ( ( MOD( nfi-1, 10 ) == 0 ) .OR. tfirst ) ) THEN
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
     IF ( .NOT. tens ) THEN
        !
        IF ( tfor .OR. tprnfor ) CALL nlfq( c0, eigr, bec, becdr, fion )
        !
     ELSE
        !
        IF ( tfor .OR. tprnfor ) &
             CALL nlfq( c0diag, eigr, becdiag, becdrdiag, fion )
        !
     END IF
     !
     IF( tfor .AND. tefield ) &
          CALL bforceion( fion, tfor, ipolp, qmat, bec, becdr, gqq, evalue )
     !
     IF( tfor .OR. thdyn ) &
          CALL interpolate_lambda( lambdap, lambda, lambdam )
     !
     ! ... calphi calculates phi
     ! ... the electron mass rises with g**2
     !
     CALL calphi( c0, ema0bg, bec, vkb, phi )
     !
     ! ... begin try and error loop (only one step!)
     !
     ! ...   nlfl and nlfh need: lambda (guessed) becdr
     !
     IF ( .NOT. tens ) THEN
        !
        IF ( tfor .OR. tprnfor ) CALL nlfl( bec, becdr, lambda, fion )
        !
     ELSE
        !
        IF ( tfor .OR. tprnfor ) THEN
           !
           CALL nlsm2( eigr, c0(:,:,1,1), becdr )
           !
           CALL nlfl( bec, becdr, lambda, fion )
           !
        END IF
        !
     END IF
     !
  END IF electron_dynamic
  !
  RETURN
  !
END SUBROUTINE move_electrons
