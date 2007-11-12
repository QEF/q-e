!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE init_run()
  !----------------------------------------------------------------------------
  !
  USE klist,              ONLY : nkstot
  USE wvfct,              ONLY : nbnd, et, wg, btype
  USE control_flags,      ONLY : lmd
  USE dynamics_module,    ONLY : allocate_dyn_vars
  USE paw_variables,      ONLY : okpaw
  USE paw_init,           ONLY : paw_init_onecenter, allocate_paw_internals
  USE bp,                 ONLY : lberry, lelfield
  !
  IMPLICIT NONE
  !
  !
  CALL start_clock( 'init_run' )
  !
  CALL setup()
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft()
  !
  ! ... generate reciprocal-lattice vectors and fft indices
  !
  CALL ggen()
  !
  CALL summary()
  !
  ! ... allocate memory for all other arrays (potentials, wavefunctions etc)
  !
  CALL allocate_nlpot()
  if (okpaw)  CALL allocate_paw_internals()
  if (okpaw)  CALL paw_init_onecenter()
  CALL allocate_locpot()
  CALL allocate_wfc()
  CALL allocate_bp_efield()
  IF( lberry .or. lelfield) call bp_global_map()
  CALL memory_report()
  !
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ), btype( nbnd, nkstot ) )
  !
  et(:,:) = 0.D0
  wg(:,:) = 0.D0
  !
  btype(:,:) = 1
  !
  CALL openfil()
  !
  CALL init_h()
  !
  IF ( lmd ) CALL allocate_dyn_vars()
  !
  CALL stop_clock( 'init_run' )
  !
  RETURN
  !
END SUBROUTINE init_run
!
!----------------------------------------------------------------------------
SUBROUTINE init_h()
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CALL hinit0()
  CALL potinit()
  !
  CALL newd()
  !
  CALL wfcinit()
  !
  RETURN
  !
END SUBROUTINE init_h
