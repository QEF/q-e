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
MODULE io_rho_xml
  !----------------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE xml_io_base, ONLY : create_directory, write_rho_xml, read_rho_xml
  !
  PRIVATE
  !
  PUBLIC :: write_rho, read_rho
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rho( rho, nspin, extension )
      !------------------------------------------------------------------------
      !
      ! ... this routine writes the charge-density in xml format into the
      ! ... '.save' directory
      ! ... the '.save' directory is created if not already present
      !
      USE io_files, ONLY : tmp_dir, prefix
      USE fft_base, ONLY : dfftp
      USE gvect,    ONLY : nr1, nr2, nr3, nrx1, nrx2, nrxx
      USE spin_orb, ONLY : domag
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN)           :: nspin
      REAL(DP),         INTENT(IN)           :: rho(nrxx,nspin)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      !
      CHARACTER(LEN=256)    :: dirname, file_base
      CHARACTER(LEN=256)    :: ext
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      !
      !
      ext = ' '
      !
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      CALL create_directory( dirname )
      !
      IF ( PRESENT( extension ) ) ext = '.' // TRIM( extension )
      !
      file_base = TRIM( dirname ) // '/charge-density' // TRIM( ext )
      !
      IF ( nspin == 1 ) THEN
         !
         CALL write_rho_xml( file_base, rho(:,1), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
      ELSE IF ( nspin == 2 ) THEN
         !
         ALLOCATE( rhoaux( nrxx ) )
         !
         rhoaux(:) = rho(:,1) + rho(:,2)
         !
         CALL write_rho_xml( file_base, rhoaux, &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         file_base = TRIM( dirname ) // '/spin-polarization' // TRIM( ext )
         !
         rhoaux(:) = rho(:,1) - rho(:,2)
         !
         CALL write_rho_xml( file_base, rhoaux, &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         DEALLOCATE( rhoaux )
         !
      ELSE IF ( nspin == 4 ) THEN
         !
         CALL write_rho_xml( file_base, rho(:,1), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         IF (domag) THEN
            file_base = TRIM( dirname ) // '/magnetization.x' // TRIM( ext )
            !
            CALL write_rho_xml( file_base, rho(:,2), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
            !
            file_base = TRIM( dirname ) // '/magnetization.y' // TRIM( ext )
            !
            CALL write_rho_xml( file_base, rho(:,3), &
                              nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
            !
            file_base = TRIM( dirname ) // '/magnetization.z' // TRIM( ext )
            !
            CALL write_rho_xml( file_base, rho(:,4), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         END IF
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_rho
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rho( rho, nspin, extension )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the charge-density in xml format from the
      ! ... files saved into the '.save' directory
      !
      USE io_files, ONLY : tmp_dir, prefix
      USE fft_base, ONLY : dfftp
      USE gvect,    ONLY : nr1, nr2, nr3, nrx1, nrx2, nrxx
      USE spin_orb, ONLY : domag
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN)           :: nspin
      REAL(DP),         INTENT(OUT)          :: rho(nrxx,nspin)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      !
      CHARACTER(LEN=256)    :: dirname, file_base
      CHARACTER(LEN=256)    :: ext
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      !
      !
      ext = ' '
      !
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      IF ( PRESENT( extension ) ) ext = '.' // TRIM( extension )
      !
      file_base = TRIM( dirname ) // '/charge-density' // TRIM( ext )
      CALL para_inquire( file_base )
      !
      IF ( nspin == 1 ) THEN
         !
         CALL read_rho_xml( file_base, rho(:,1), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
      ELSE IF ( nspin == 2 ) THEN
         !
         ALLOCATE( rhoaux( nrxx ) )
         !
         CALL read_rho_xml( file_base, rhoaux, &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         rho(:,1) = rhoaux(:)
         rho(:,2) = rhoaux(:)
         !
         file_base = TRIM( dirname ) // '/spin-polarization' // TRIM( ext )
         CALL para_inquire( file_base )
         !
         CALL read_rho_xml( file_base, rhoaux, &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         rho(:,1) = 0.5D0*( rho(:,1) + rhoaux(:) )
         rho(:,2) = 0.5D0*( rho(:,2) - rhoaux(:) )
         !
         DEALLOCATE( rhoaux )
         !
      ELSE IF ( nspin == 4 ) THEN
         !
         CALL read_rho_xml( file_base, rho(:,1), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         IF ( domag ) THEN
            !
            file_base = TRIM( dirname ) // '/magnetization.x' // TRIM( ext )
            !
            CALL para_inquire( file_base )
            !
            CALL read_rho_xml( file_base, rho(:,2), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
            !
            file_base = TRIM( dirname ) // '/magnetization.y' // TRIM( ext )
            !
            CALL para_inquire( file_base )
            !
            CALL read_rho_xml( file_base, rho(:,3), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
            !
            file_base = TRIM( dirname ) // '/magnetization.z' // TRIM( ext )
            !
            CALL para_inquire( file_base )
            !
            CALL read_rho_xml( file_base, rho(:,4), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
            !
         ELSE
            !
            rho(:,2:4) = 0.D0
            !
         END IF
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_rho
    !
    !------------------------------------------------------------------------
    SUBROUTINE para_inquire( file_base )
      !------------------------------------------------------------------------
      !
      ! ... same as fortran function 'inquire' but the check is performed
      ! ... on a single processor whether 'file_base'.xml exists
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      USE mp,        ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      LOGICAL          :: lexists
      CHARACTER(LEN=*) :: file_base
      !
      IF ( ionode ) &
           INQUIRE( FILE = TRIM( file_base ) // '.xml', EXIST = lexists )
      !
      CALL mp_bcast ( lexists, ionode_id, intra_image_comm )
      !
      IF ( .NOT. lexists ) &
         CALL errore( 'read_rho', 'file ' // &
                    & TRIM( file_base ) // '.xml nonexistent', 1 )
      !
      RETURN
      !
    END SUBROUTINE para_inquire
    !
END MODULE io_rho_xml
