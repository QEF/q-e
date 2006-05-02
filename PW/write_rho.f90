!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE write_rho( rho )
  !----------------------------------------------------------------------------
  !
  ! ... this routine writes the charge-density in xml format into the
  ! ... '.save' directory. The '.save' directory is created if not already
  ! ... present
  !
  USE kinds,       ONLY : DP
  USE io_files,    ONLY : tmp_dir, prefix
  USE sticks,      ONLY : dfftp
  USE gvect,       ONLY : nr1, nr2, nr3, nrx1, nrx2, nrxx
  USE lsda_mod,    ONLY : nspin
  USE xml_io_base, ONLY : create_directory, write_rho_xml
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho(nrxx,nspin)
  !
  CHARACTER(LEN=256)    :: dirname, rho_file_base
  REAL(DP), ALLOCATABLE :: rhosum(:)
  !
  !
  dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  CALL create_directory( dirname )
  !
  rho_file_base = TRIM( dirname ) // '/charge-density'
  !
  IF ( nspin == 1 ) THEN
     !
     CALL write_rho_xml( rho_file_base, rho(:,1), &
                         nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     ALLOCATE( rhosum( SIZE( rho, 1 ) ) )
     !
     rhosum = rho(:,1) + rho(:,2) 
     !
     CALL write_rho_xml( rho_file_base, rhosum, &
                         nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
     !
     DEALLOCATE( rhosum )
     !
     rho_file_base = TRIM( dirname ) // '/charge-density-up'
     !
     CALL write_rho_xml( rho_file_base, rho(:,1), &
                         nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
     !
     rho_file_base = TRIM( dirname ) // '/charge-density-dw'
     !
     CALL write_rho_xml( rho_file_base, rho(:,2), &
                         nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     CALL write_rho_xml( rho_file_base, rho(:,1), &
                         nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
     !
     rho_file_base = TRIM( dirname ) // '/magnetization.x'
     !
     CALL write_rho_xml( rho_file_base, rho(:,2), &
                         nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
     !
     rho_file_base = TRIM( dirname ) // '/magnetization.y'
     !
     CALL write_rho_xml( rho_file_base, rho(:,3), &
                          nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
     !
     rho_file_base = TRIM( dirname ) // '/magnetization.z'
     !
     CALL write_rho_xml( rho_file_base, rho(:,4), &
                         nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
     !
  END IF  
  !
  RETURN
  !
END SUBROUTINE write_rho
