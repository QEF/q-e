!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE metadyn_io
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the I/O methods used by meta-dynamics
  !
  ! ... code written by Carlo Sbraccia (2005)
  !
  USE kinds, ONLY : DP
  !
  USE iotk_module
  !
  IMPLICIT NONE
  !
  CHARACTER(iotk_attlenx) :: attr
  !
  PRIVATE
  !
  PUBLIC :: write_axsf_file
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_axsf_file( image, tau, tau_units )
      !------------------------------------------------------------------------
      !
      USE io_files,         ONLY : iunaxsf
      USE constants,        ONLY : bohr_radius_angs
      USE ions_base,        ONLY : nat, ityp, atm
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN) :: image
      REAL(DP), INTENT(IN) :: tau(:,:)
      REAL(DP), INTENT(IN) :: tau_units
      !
      INTEGER :: ia
      LOGICAL :: opnd
      !
      !
      INQUIRE( UNIT = iunaxsf, OPENED = opnd )
      IF ( .NOT.opnd ) &
         CALL errore( 'write_axsf_file', &
                      'unit to write the axsf file is closed', 1 )
      !
      WRITE( UNIT = iunaxsf, FMT = '(" PRIMCOORD ",I5)' ) image
      WRITE( UNIT = iunaxsf, FMT = '(I5,"  1")' ) nat
      !
      DO ia = 1, nat
         !
         WRITE( UNIT = iunaxsf, FMT = '(A2,3(2X,F18.10))' ) &
                TRIM( atm(ityp(ia)) ), &
             tau(1,ia)*tau_units*bohr_radius_angs, &
             tau(2,ia)*tau_units*bohr_radius_angs, &
             tau(3,ia)*tau_units*bohr_radius_angs
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE write_axsf_file
    !
END MODULE metadyn_io
