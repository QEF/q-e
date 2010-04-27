!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE exx_loop( )
  !----------------------------------------------------------------------------
  !
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : conv_elec, conv_ions
  USE io_files,         ONLY : prefix, tmp_dir
  USE io_global,        ONLY : stdout, meta_ionode, ionode, ionode_id
  USE mp_global,        ONLY : my_image_id, nimage,  me_image, root_image
  USE mp,               ONLY : mp_barrier
  !
  IMPLICIT NONE
  !
  REAL(DP), EXTERNAL :: get_clock
  !
  CHARACTER (LEN=256)   :: tmp_dir_saved
  LOGICAL               :: opnd, exst
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  !
  CALL mp_barrier()
  !
  IF ( nimage > 1 ) THEN
     !
     ! ... Image parallelization of exchange: reset I/O nodes, one per image
     !
     ionode = ( me_image == root_image )
     ionode_id = root_image
     !
     ! ... one scratch directory per image
     !
     tmp_dir_saved = tmp_dir
     tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
          TRIM( int_to_char( my_image_id + 1 ) ) // "/"
     !
  END IF
  !
  CALL setup()
  CALL init_run()
  !
  main_loop: DO
     !
     ! ... electronic self-consistentcy
     !
     CALL electrons()
     !
     IF ( .NOT. conv_elec ) CALL stop_run( conv_elec )
     !
     ! ... if requested ions are moved
     !
     CALL ions()
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( conv_ions ) EXIT main_loop
     !
     ! ... the ionic part of the hamiltonian is reinitialized
     !
     CALL hinit1()
     !
  END DO main_loop
  !
  CALL punch( 'all' )
  !
1 CALL mp_barrier()
  !
  CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  ! ... all other files must be reopened and removed
  !
  CALL seqopn( 4, 'update', 'FORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  CALL seqopn( 4, 'para', 'FORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  IF ( nimage > 1 ) tmp_dir = tmp_dir_saved
  !
  RETURN
  !
END SUBROUTINE 
