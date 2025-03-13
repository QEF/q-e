!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE cpr_loop( nloop )
  !----------------------------------------------------------------------------
  !! Set \(\text{ion_positions} = \text{from_input}\) and \(\text{rd_pos} = 
  !! \text{your_positions}\) to force \(\texttt{cprmain}\) to compute forces
  !! for \(\text{your_position}\) configuration.
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
     CALL cprmain( tau, fion, etot )
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
