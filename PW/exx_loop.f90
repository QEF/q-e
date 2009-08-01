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
  USE io_files,         ONLY : prefix, tmp_dir, iunpath
  USE path_formats,     ONLY : scf_fmt, scf_fmt_para
  USE io_global,        ONLY : stdout, ionode
  USE mp_global,        ONLY : my_image_id, nimage
  USE mp,               ONLY : mp_barrier
  !
  IMPLICIT NONE
  !
  REAL(DP), EXTERNAL :: get_clock
  !
  REAL(DP)              :: tcpu
  CHARACTER (LEN=256)   :: tmp_dir_saved
  LOGICAL               :: opnd
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  CALL flush_unit( iunpath )
  !
  tmp_dir_saved = tmp_dir
  !
  CALL mp_barrier()
  !
  CALL clean_pw( .FALSE. )
  !
  tcpu = get_clock( 'PWSCF' )
  !
  IF ( nimage > 1 ) THEN
         !
         WRITE( UNIT = iunpath, FMT = scf_fmt_para ) my_image_id, tcpu, 0
         tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
              TRIM( int_to_char( my_image_id + 1 ) ) // "/"
         !
  ELSE
         !
         WRITE( UNIT = iunpath, FMT = scf_fmt ) tcpu, 0
         !
  END IF
  !

  !
  ! ... unit stdout is connected to the appropriate file
  !
  IF ( ionode .and. nimage>1) THEN
         !
         INQUIRE( UNIT = stdout, OPENED = opnd )
         IF ( opnd ) CLOSE( UNIT = stdout )
         OPEN( UNIT = stdout, FILE = TRIM( tmp_dir ) // 'PW.out', &
               STATUS = 'UNKNOWN', POSITION = 'APPEND' )
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
 ! IF (flag) THEN
     !write(*,*) "sono nello stop_run"  
     CALL seqopn( 4, 'restart', 'UNFORMATTED', conv_ions )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !ENDIF

  !IF ( flag .AND. ionode ) THEN
     !
     ! ... all other files must be reopened and removed
     !
     CALL seqopn( 4, 'update', 'FORMATTED', conv_ions )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     CALL seqopn( 4, 'para', 'FORMATTED', conv_ions )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     CALL seqopn( 4, 'BLOCK', 'FORMATTED', conv_ions )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
  !END IF

  tmp_dir = tmp_dir_saved
  !
  RETURN
  !
END SUBROUTINE 
