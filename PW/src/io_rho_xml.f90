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
  USE xml_io_base, ONLY : create_directory, write_rho, read_rho
  !
  PRIVATE
  !
  PUBLIC :: write_scf, read_scf
  !
  ! {read|write}_rho: read or write the charge density
  ! {read|write}_scf: as above, plus ldaU ns, PAW becsum, meta-GGA

  CONTAINS

    SUBROUTINE write_scf ( rho, nspin )
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u
      USE funct,            ONLY : dft_is_meta
      USE noncollin_module, ONLY : noncolin
      USE spin_orb,         ONLY : domag
      USE io_files,         ONLY : seqopn
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE scf,              ONLY : scf_type
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast

      !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(IN)           :: rho
      INTEGER,          INTENT(IN)           :: nspin
      LOGICAL :: lexist
      INTEGER :: iunocc, iunpaw, ierr
      INTEGER, EXTERNAL :: find_free_unit

      ! Use the equivalent routine to write real space density

      ! in the following case do not read or write polarization
      IF ( noncolin .AND. .NOT.domag ) THEN
         CALL write_rho ( rho%of_r, 1 )
      ELSE
         CALL write_rho ( rho%of_r, nspin )
      END IF

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
         IF ( ierr/=0 ) CALL errore('write_scf', 'Writing ldaU ns', 1)
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
         IF ( ierr/=0 ) CALL errore('write_scf', 'Writing PAW becsum',1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunpaw, STATUS = 'KEEP' )
         ENDIF
         !
      END IF
      !
      IF ( dft_is_meta() ) THEN
         WRITE(stdout,'(5x,"Writing meta-gga kinetic term")')
          CALL write_rho ( rho%kin_r, nspin, 'kin' )
      ENDIF

      RETURN
    END SUBROUTINE write_scf

    SUBROUTINE read_scf ( rho, nspin )
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u, starting_ns
      USE noncollin_module, ONLY : noncolin
      USE spin_orb,         ONLY : domag
      USE funct,            ONLY : dft_is_meta
      USE io_files,         ONLY : seqopn, prefix, tmp_dir
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE scf,              ONLY : scf_type
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast, mp_sum
      !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(INOUT)        :: rho
      INTEGER,          INTENT(IN)           :: nspin
      CHARACTER(LEN=256) :: dirname
      LOGICAL :: lexist
      INTEGER :: iunocc, iunpaw, ierr
      INTEGER, EXTERNAL :: find_free_unit

      dirname = TRIM(tmp_dir) // TRIM(prefix) // '.save/'
      ! in the following case do not read or write polarization
      IF ( noncolin .AND. .NOT.domag ) THEN
         CALL read_rho ( dirname, rho%of_r, 1 )
         rho%of_r(:,2:4) = 0.0_dp
      ELSE
         CALL read_rho ( dirname, rho%of_r, nspin )
      END IF
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
         IF ( ierr/=0 ) CALL errore('read_scf', 'Reading ldaU ns', 1)
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
         IF ( ierr/=0 ) CALL errore('read_scf', 'Reading PAW becsum',1)
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
         CALL read_rho( dirname, rho%kin_r, nspin, 'kin' )
      END IF

      RETURN
    END SUBROUTINE read_scf
    !
END MODULE io_rho_xml
