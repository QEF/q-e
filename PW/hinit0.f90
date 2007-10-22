!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hinit0()
  !-----------------------------------------------------------------------
  !
  ! ... configuration-independent hamiltonian initialization
  !
  USE ions_base, ONLY : nat, nsp, ityp, tau
  USE basis,     ONLY : startingconfig
  USE cell_base, ONLY : at, bg, omega, tpiba2
  USE cellmd,    ONLY : omega_old, at_old, lmovecell
  USE klist,     ONLY : nks, xk
  USE gvect,     ONLY : nr1, nr2, nr3, ngm, ecutwfc, ig_l2g, &
                        g, eigts1, eigts2, eigts3
  USE vlocal,    ONLY : strf
  USE wvfct,     ONLY : npw, g2kin, igk
  USE io_files,  ONLY : iunigk
  USE realus,    ONLY : tqr, qpointlist
#ifdef __GRID_PAW
  USE grid_paw_routines,  ONLY : init_prad, set_paw_rhoc, init_paw_vloc, paw_grid_setlocal
#endif
  USE grid_paw_variables, ONLY : okpaw
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  ! counter on k points
  !
  ! ... calculate the local part of the pseudopotentials
  !
  CALL init_vloc()
#ifdef __GRID_PAW
  IF (okpaw) CALL init_paw_vloc() !!PAW!!
#endif
  !
  ! ... k-point independent parameters of non-local pseudopotentials
  !
  CALL init_us_1()
#ifdef __GRID_PAW
  IF (okpaw) CALL init_prad() !!PAW!!
#endif
  CALL init_at_1()
  !
  REWIND( iunigk )
  !
  DO ik = 1, nks
     !
     ! ... g2kin is used here as work space
     !
     CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
     !
     ! ... if there is only one k-point npw and igk stay in memory
     !
     IF ( nks > 1 ) WRITE( iunigk ) igk
     !
  END DO
  !
  IF ( lmovecell .AND. startingconfig == 'file' ) THEN
     !
     ! ... If lmovecell and restart are both true the cell shape is read from
     ! ... the restart file and stored. The xxx_old variables are used instead 
     ! ... of the current (read from input) ones.
     ! ... xxx and xxx_old are swapped, the atomic positions rescaled and 
     ! ... the hamiltonian scaled.
     !
     CALL cryst_to_cart( nat, tau, bg, - 1 )
     !
     CALL dswap( 9, at, 1, at_old, 1 )
     CALL dswap( 1, omega, 1, omega_old, 1 )
     !
     CALL cryst_to_cart( nat, tau, at, + 1 )
     !
     CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     CALL scale_h()
     !
  END IF
  !
  ! ... initialize the structure factor
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                   nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()
#ifdef __GRID_PAW
  IF (okpaw) CALL paw_grid_setlocal() !!PAW!!
#endif
  !
  ! ... calculate the core charge (if any) for the nonlinear core correction
  !
  CALL set_rhoc()
#ifdef __GRID_PAW
  IF (okpaw) CALL set_paw_rhoc() !!PAW!!
#endif
  !
  IF ( tqr ) CALL qpointlist()
  !
  RETURN
  !
END SUBROUTINE hinit0

