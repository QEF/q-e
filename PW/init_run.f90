!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE init_run
  !-----------------------------------------------------------------------
  !
  USE parameters,  ONLY : ntypx, npk, lmaxx, nchix, ndm, nbrx, nqfm
  USE wvfct,       ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  !
  CALL start_clock( 'init_run' )
  !
  IF ( gamma_only ) THEN
     WRITE(6, '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials, ", &
                 &  "Gamma point")')
  ELSE
     WRITE(6, '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials")')
  END IF
  !
  WRITE(6, 9010) ntypx, npk, lmaxx, nchix, ndm, nbrx, nqfm
  !
  CALL iosys
  CALL setup
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft
  !
  ! ... generate reciprocal-lattice vectors and fft indices
  !
  CALL ggen
  CALL summary
  CALL allocate_nlpot
  CALL allocate_locpot
  CALL allocate_wfc
  !
  CALL openfil
  !
  CALL hinit0
  CALL potinit
  !
  CALL newd
  !
  CALL wfcinit
  !
  CALL stop_clock ( 'init_run' )
  WRITE(6, * )
  !
  CALL show_memory()
  !
  RETURN
  !
9010 FORMAT( /5X,'Current dimensions of program pwscf are:' &
             /5X,'ntypx =',I2,'   npk =',I5,'  lmax =',I2   &
             /5X,'nchix =',I2,'  ndim =',I5,'  nbrx =',I2,' nqfm =',I2 )
  !
END SUBROUTINE init_run

