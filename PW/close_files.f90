!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE close_files()
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes for a new scf calculation.
  !
  USE varie,    ONLY :  order
  USE io_files, ONLY :  prefix, iunwfc, iunoldwfc, iunoldwfc2, iunigk
#ifdef __PARA
  USE mp,       ONLY :  mp_barrier
#endif
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  LOGICAL :: exst
  !
  !  ... iunwfc contains wavefunctions and is kept open during
  !  ... the execution - close and save the file
  !
  CLOSE( UNIT = iunwfc, STATUS = 'KEEP' )
  !
  IF ( order > 1 ) &
     CLOSE( UNIT = iunoldwfc, STATUS = 'KEEP' ) 
  !
  IF ( order > 2 ) &
     CLOSE( UNIT = iunoldwfc2, STATUS = 'KEEP' )   
  !
  ! ... iunigk is kept open during the execution - close and remove
  !
  CLOSE( UNIT = iunigk, STATUS = 'DELETE' )
  !
#ifdef __PARA
  CALL mp_barrier()  
#endif
  !
#ifdef __T3E
  !
  ! ... set streambuffers off
  !
  CALL set_d_stream( 0 )
#endif
  !
  RETURN
  !
END SUBROUTINE close_files
