!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
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
  USE path_io_routines, ONLY : path_summary
  USE image_io_routines, ONLY : io_image_start, io_image_stop
  !
  IMPLICIT NONE
  !
  CALL iosys_path()
  !
  CALL io_image_start()
  !
  CALL initialize_path()
  !
  CALL path_summary()
  !
  CALL search_mep()
  !
  CALL io_image_stop()
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
  CALL init_run()
  !
  DO iloop = 1, nloop
     !
     CALL cprmain( tau(1,1), fion(1,1), etot )
     !
  END DO
  !
  CALL terminate_run()
  !
  DEALLOCATE( tau, fion )
  !
  RETURN
  !
END SUBROUTINE cpr_loop
