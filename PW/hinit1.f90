!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE hinit1()
  !----------------------------------------------------------------------------
  !
  ! ... Atomic configuration dependent hamiltonian initialization
  !
  USE ions_base,     ONLY : nat, nsp, ityp, tau
  USE cell_base,     ONLY : at, bg, omega, tpiba2
  USE cellmd,        ONLY : lmovecell 
  USE gvect,         ONLY : nr1, nr2, nr3, nrxx, ngm, g, eigts1, eigts2, eigts3
  USE gsmooth,       ONLY : doublegrid
  USE ldaU,          ONLY : lda_plus_u
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : vrs, vltot, v, kedtau
  USE vlocal,        ONLY : strf
  USE control_flags, ONLY : pot_order, tqr
  USE realus,        ONLY : qpointlist
  !
  IMPLICIT NONE
  !
  !
  ! ... update the potential
  !
  CALL update_pot()
  !
  ! ... initialize structure factor array if it has not already been
  ! ... calculated in update_pot ( i.e. when pot_order > 0 )
  !
  IF ( pot_order == 0 ) THEN
     !
     IF ( lmovecell ) CALL scale_h()
     !
     CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                      nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
     !
     ! ... calculate the core charge (if any) for the nonlinear 
     ! ... core correction
     !
     CALL set_rhoc()
     !
  END IF
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()
  !
  ! ... define the total local potential (external+scf)
  !
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid )
  !
  IF ( tqr ) CALL qpointlist()
  !
  ! ... update the D matrix
  !
  CALL newd()
  !
  ! ... and recalculate the products of the S with the atomic wfcs used 
  ! ... in LDA+U calculations
  !
  IF ( lda_plus_u ) CALL orthoatwfc()
  !
  RETURN
  !
END SUBROUTINE hinit1

