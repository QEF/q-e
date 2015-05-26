!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hinit0()
  !-----------------------------------------------------------------------
  !
  ! ... hamiltonian initialization: 
  ! ... atomic position independent initialization for nonlocal PP,
  ! ... structure factors, local potential, core charge
  !
  USE ions_base,    ONLY : nat, nsp, ityp, tau
  USE basis,        ONLY : startingconfig
  USE cell_base,    ONLY : at, bg, omega, tpiba2
  USE cellmd,       ONLY : omega_old, at_old, lmovecell
  USE klist,        ONLY : nks, xk
  USE fft_base,     ONLY : dfftp
  USE gvect,        ONLY : ngm, ig_l2g, g, eigts1, eigts2, eigts3
  USE vlocal,       ONLY : strf
  USE wvfct,        ONLY : npw, g2kin, igk, ecutwfc
  USE io_files,     ONLY : iunigk
  USE realus,       ONLY : generate_qpointlist,betapointlist,init_realspace_vars,real_space
  use ldaU,         ONLY : lda_plus_U, U_projection
  USE control_flags,ONLY : tqr 
  USE io_global,    ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  ! counter on k points
  !
  ! ... calculate the Fourier coefficients of the local part of the PP
  !
  CALL init_vloc()
  !
  ! ... k-point independent parameters of non-local pseudopotentials
  !
  CALL init_us_1()
  IF ( lda_plus_U .AND. ( U_projection == 'pseudo' ) ) CALL init_q_aeps()
  CALL init_at_1()
  !
  REWIND( iunigk )
  !
  ! ... The following loop must NOT be called more than once in a run
  ! ... or else there will be problems with variable-cell calculations
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
                   dfftp%nr1, dfftp%nr2, dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  !
  ! these routines can be used to patch quantities that are dependent
  ! on the ions and cell parameters
  !
  CALL plugin_init_ions()
  CALL plugin_init_cell()
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()
  !
  ! ... calculate the core charge (if any) for the nonlinear core correction
  !
  CALL set_rhoc()
  !
  IF ( tqr ) CALL generate_qpointlist()
  IF (real_space ) then
   call betapointlist()
   call init_realspace_vars()
   write(stdout,'(5X,"Real space initialisation completed")')    
  endif
  !
  RETURN
  !
END SUBROUTINE hinit0

