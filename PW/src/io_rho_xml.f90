!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
  ! {read|write}_rho_only:    read or write the real space charge density
  ! {read|write}_rho_general: as above, plus read or write ldaU ns coeffs
  !                           and PAW becsum coeffs.

  INTERFACE write_rho
        MODULE PROCEDURE write_rho_only, write_rho_general
  END INTERFACE

  INTERFACE read_rho
        MODULE PROCEDURE read_rho_only, read_rho_general
  END INTERFACE

  CONTAINS

    SUBROUTINE write_rho_general( rho, nspin, extension )
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u
      USE funct,            ONLY : dft_is_meta
      USE noncollin_module, ONLY : noncolin
      USE io_files,         ONLY : seqopn
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE scf,              ONLY : scf_type
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast

      !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(IN)           :: rho
      INTEGER,          INTENT(IN)           :: nspin
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      LOGICAL :: lexist
      INTEGER :: iunocc, iunpaw, ierr
      INTEGER, EXTERNAL :: find_free_unit

      ! Use the equivalent routine to write real space density
      CALL write_rho_only( rho%of_r, nspin, extension )

      ! Then write the other terms to separate files

      IF ( lda_plus_u ) THEN
         !
         iunocc = find_free_unit ()
         IF ( ionode ) THEN
            CALL seqopn( iunocc, 'save/occup.txt', 'FORMATTED', lexist )
            if (noncolin) then
              WRITE( iunocc, * , iostat = ierr) rho%ns_nc
            else
              WRITE( iunocc, * , iostat = ierr) rho%ns
            endif
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('write_rho_general', 'Writing ldaU ns', 1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunocc, STATUS = 'KEEP' )
         ENDIF
         !
      END IF
      !
      IF ( okpaw ) THEN
         !
         iunpaw = find_free_unit ()
         IF ( ionode ) THEN
            CALL seqopn( iunpaw, 'save/paw.txt', 'FORMATTED', lexist )
            WRITE( iunpaw, * , iostat = ierr) rho%bec
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('write_rho_general', 'Writing PAW becsum',1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunpaw, STATUS = 'KEEP' )
         ENDIF
         !
      END IF
      !
      IF ( dft_is_meta() ) THEN
         WRITE(stdout,'(5x,"Writing meta-gga kinetic term")')
          CALL write_rho_only( rho%kin_r, nspin, 'kin' )
      ENDIF

      RETURN
    END SUBROUTINE write_rho_general

    SUBROUTINE read_rho_general( rho, nspin, extension )
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u, starting_ns
      USE noncollin_module, ONLY : noncolin
      USE funct,            ONLY : dft_is_meta
      USE io_files,         ONLY : seqopn
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE scf,              ONLY : scf_type
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast, mp_sum
      !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(INOUT)        :: rho
      INTEGER,          INTENT(IN)           :: nspin
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      LOGICAL :: lexist
      INTEGER :: iunocc, iunpaw, ierr
      INTEGER, EXTERNAL :: find_free_unit

      ! Use the equivalent routine to read real space density
      CALL read_rho_only( rho%of_r, nspin, extension )
      !
      IF ( lda_plus_u ) THEN
         !
         ! The occupations ns also need to be read in order to build up
         ! the potential
         !
         iunocc = find_free_unit ()
         IF ( ionode ) THEN
            CALL seqopn( iunocc, 'save/occup.txt', 'FORMATTED', lexist )
            if (noncolin) then
              READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns_nc
            else
              READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns
            endif
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('read_rho_general', 'Reading ldaU ns', 1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunocc, STATUS = 'KEEP')
         ELSE
            if (noncolin) then
              rho%ns_nc(:,:,:,:) = 0.D0
            else
              rho%ns(:,:,:,:) = 0.D0
            endif 
         END IF
         if (noncolin) then
           CALL mp_sum(rho%ns_nc, intra_image_comm)
         else
           CALL mp_sum(rho%ns, intra_image_comm)
         endif
         ! If projections on Hubbard manifold are read from file, there is no
         ! need to set starting values: reset them to prevent any problem
         starting_ns = -1.0_dp
      END IF
      !
      IF ( okpaw ) THEN
         !
         ! Also the PAW coefficients are needed:
         !
         iunpaw = find_free_unit ()
         IF ( ionode ) THEN
            CALL seqopn( iunpaw, 'save/paw.txt', 'FORMATTED', lexist )
            READ( UNIT = iunpaw, FMT = *, iostat=ierr ) rho%bec
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('read_rho_general', 'Reading PAW becsum',1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunpaw, STATUS = 'KEEP')
         ELSE
            rho%bec(:,:,:) = 0.D0
         END IF
         CALL mp_sum(rho%bec, intra_image_comm)
         !
      END IF
      !
      IF ( dft_is_meta() ) THEN
         WRITE(stdout,'(5x,"Reading meta-gga kinetic term")')
         CALL read_rho_only( rho%kin_r, nspin, 'kin' )
      END IF

      RETURN
    END SUBROUTINE read_rho_general
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rho_only( rho, nspin, extension )
      !------------------------------------------------------------------------
      !
      ! ... this routine writes the charge-density in xml format into the
      ! ... '.save' directory
      ! ... the '.save' directory is created if not already present
      !
      USE io_files, ONLY : tmp_dir, prefix
      USE fft_base, ONLY : dfftp
      USE spin_orb, ONLY : domag
      USE io_global,ONLY : ionode
      USE mp_bands, ONLY : intra_bgrp_comm, inter_bgrp_comm
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN)           :: nspin
      REAL(DP),         INTENT(IN)           :: rho(dfftp%nnr,nspin)
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
         CALL write_rho_xml( file_base, rho(:,1), dfftp%nr1, dfftp%nr2, &
                  dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                  ionode, intra_bgrp_comm, inter_bgrp_comm )
         !
      ELSE IF ( nspin == 2 ) THEN
         !
         ALLOCATE( rhoaux( dfftp%nnr ) )
         !
         rhoaux(:) = rho(:,1) + rho(:,2)
         !
         CALL write_rho_xml( file_base, rhoaux, dfftp%nr1, dfftp%nr2, &
                  dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                  ionode, intra_bgrp_comm, inter_bgrp_comm )
         !
         file_base = TRIM( dirname ) // '/spin-polarization' // TRIM( ext )
         !
         rhoaux(:) = rho(:,1) - rho(:,2)
         !
         CALL write_rho_xml( file_base, rhoaux,  dfftp%nr1, dfftp%nr2, &
                  dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                  ionode, intra_bgrp_comm, inter_bgrp_comm )
         !
         DEALLOCATE( rhoaux )
         !
      ELSE IF ( nspin == 4 ) THEN
         !
         CALL write_rho_xml( file_base, rho(:,1), dfftp%nr1, dfftp%nr2, &
                  dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                  ionode, intra_bgrp_comm, inter_bgrp_comm )
         !
         IF (domag) THEN
            file_base = TRIM( dirname ) // '/magnetization.x' // TRIM( ext )
            !
            CALL write_rho_xml( file_base, rho(:,2), dfftp%nr1, dfftp%nr2, &
                  dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                  ionode, intra_bgrp_comm, inter_bgrp_comm )
            !
            file_base = TRIM( dirname ) // '/magnetization.y' // TRIM( ext )
            !
            CALL write_rho_xml( file_base, rho(:,3), dfftp%nr1, dfftp%nr2, &
                  dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                  ionode, intra_bgrp_comm, inter_bgrp_comm )
            !
            file_base = TRIM( dirname ) // '/magnetization.z' // TRIM( ext )
            !
            CALL write_rho_xml( file_base, rho(:,4), dfftp%nr1, dfftp%nr2, &
                  dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                  ionode, intra_bgrp_comm, inter_bgrp_comm )
         END IF
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_rho_only
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rho_only( rho, nspin, extension )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the charge-density in xml format from the
      ! ... files saved into the '.save' directory
      !
      USE io_files,  ONLY : tmp_dir, prefix
      USE fft_base,  ONLY : dfftp
      USE spin_orb,  ONLY : domag
      USE io_global, ONLY : ionode
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN)           :: nspin
      REAL(DP),         INTENT(OUT)          :: rho(dfftp%nnr,nspin)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      !
      CHARACTER(LEN=256)    :: dirname, file_base
      CHARACTER(LEN=256)    :: ext
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      !
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      ext = ' '
      IF ( PRESENT( extension ) ) ext = '.' // TRIM( extension )
      file_base = TRIM( dirname ) // '/charge-density' // TRIM( ext )
      !
      CALL read_rho_xml ( file_base, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                 dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, rho(:,1) ) 
      !
      IF ( nspin == 2 ) THEN
         !
         rho(:,2) = rho(:,1)
         !
         ALLOCATE( rhoaux( dfftp%nnr ) )
         !
         file_base = TRIM( dirname ) // '/spin-polarization' // TRIM( ext )
         CALL read_rho_xml ( file_base, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                    dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, rhoaux ) 
         !
         rho(:,1) = 0.5D0*( rho(:,1) + rhoaux(:) )
         rho(:,2) = 0.5D0*( rho(:,2) - rhoaux(:) )
         !
         DEALLOCATE( rhoaux )
         !
      ELSE IF ( nspin == 4 ) THEN
         !
         IF ( domag ) THEN
            !
            file_base = TRIM( dirname ) // '/magnetization.x' // TRIM( ext )
            CALL read_rho_xml ( file_base, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                 dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, rho(:,2) ) 
            !
            file_base = TRIM( dirname ) // '/magnetization.y' // TRIM( ext )
            CALL read_rho_xml ( file_base, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                 dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, rho(:,3) ) 
            !
            file_base = TRIM( dirname ) // '/magnetization.z' // TRIM( ext )
            CALL read_rho_xml ( file_base, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                 dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, rho(:,4) ) 
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
    END SUBROUTINE read_rho_only
    !
END MODULE io_rho_xml
