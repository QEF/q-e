!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... This file contains "sub-main" subroutines that drive the different
! ... kinds of "meta"/"non-meta" dynamics
!
! ... set ion_positions = 'from_input' and rd_pos = +your_positions+
! ... to force cprmain to compute forces for +your_position+ configuration
!
!----------------------------------------------------------------------------
SUBROUTINE neb_loop( )
  !----------------------------------------------------------------------------
  !
  USE path_base,        ONLY : initialize_path, search_mep
  USE path_routines,    ONLY : iosys_path
  USE path_io_routines, ONLY : io_path_start, io_path_stop, path_summary
  !
  IMPLICIT NONE
  !
  CALL iosys_path()
  !
  CALL io_path_start()
  !
  CALL initialize_path()
  !
  CALL path_summary()
  !
  CALL search_mep()
  !
  CALL io_path_stop()
  !
  RETURN
  !
END SUBROUTINE neb_loop
!
!----------------------------------------------------------------------------
SUBROUTINE cpr_loop( nloop )
  !----------------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE ions_base,       ONLY : nat
  USE control_flags,   ONLY : lmetadyn, program_name
  USE metadyn_base,    ONLY : metadyn_init
  USE cp_interfaces,   ONLY : main_fpmd
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nloop
  !
  INTEGER               :: iloop
  REAL(DP), ALLOCATABLE :: tau(:,:)
  REAL(DP), ALLOCATABLE :: fion(:,:)
  REAL(DP)              :: etot
  !
  !
  IF ( nat > 0 ) THEN
     !
     ALLOCATE( tau(  3, nat ) )
     ALLOCATE( fion( 3, nat ) )
     !
  ELSE
     !
     CALL errore( ' cpr_loop ', ' nat less or equal 0 ', 1 )
     !
  END IF
  !
  IF ( lmetadyn ) THEN
     !
     CALL metadyn_init( 'CP', tau )
     !
     CALL metadyn()
     !
  ELSE
     !
     CALL init_run()
     !
     DO iloop = 1, nloop
        !
        IF( program_name == 'CP90' ) THEN
           !
           CALL cprmain( tau(1,1), fion(1,1), etot )
           !
        ELSE
           !
           CALL main_fpmd( tau, fion, etot )
           !
        END IF
        !
     END DO
     !
  END IF
  !
  CALL terminate_run()
  !
  DEALLOCATE( tau, fion )
  !
  RETURN
  !
END SUBROUTINE cpr_loop
