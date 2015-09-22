!
! Copyright (C) 2005 PWSCF-FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!----------------------------------------------------------------------------
SUBROUTINE flush_unit( unit_tobeflushed )
  !----------------------------------------------------------------------------
  !
  ! ... this is a wrapper to the standard flush routine
  !
  INTEGER, INTENT(IN) :: unit_tobeflushed
  LOGICAL             :: opnd
  !
  !
  INQUIRE( UNIT = unit_tobeflushed, OPENED = opnd )
  !
#if defined(__XLF)
  IF ( opnd ) CALL flush_( unit_tobeflushed )
#else
#if defined(__NAG)
  IF ( opnd ) FLUSH( unit_tobeflushed )
#else
  IF ( opnd ) CALL flush( unit_tobeflushed )
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE
