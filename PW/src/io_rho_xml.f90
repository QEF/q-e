!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
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
  USE io_files,    ONLY : create_directory
  USE io_base,     ONLY : write_rhog, read_rhog
  !
  PRIVATE
  PUBLIC :: write_scf, read_scf
  !
  ! {read|write}_rho: read or write the charge density
  ! {read|write}_scf: as above, plus ldaU ns, PAW becsum, meta-GGA
  !
  CONTAINS

    SUBROUTINE write_scf ( rho, nspin )
      !
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u, hub_back, lda_plus_u_kind, nsg
      USE funct,            ONLY : dft_is_meta
      USE noncollin_module, ONLY : noncolin
      USE spin_orb,         ONLY : domag
      USE scf,              ONLY : scf_type
      !
      USE cell_base,        ONLY : bg, tpiba
      USE gvect,            ONLY : ig_l2g, mill
      USE control_flags,    ONLY : gamma_only
      USE io_files,         ONLY : restart_dir
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE mp_pools,         ONLY : my_pool_id
      USE mp_bands,         ONLY : my_bgrp_id, root_bgrp_id, &
                                   root_bgrp, intra_bgrp_comm
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast

      !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(IN)           :: rho
      INTEGER,          INTENT(IN)           :: nspin
      !
      CHARACTER (LEN=256) :: dirname
      INTEGER :: nspin_, iunocc, iunpaw, ierr
      INTEGER, EXTERNAL :: find_free_unit

      dirname = restart_dir ( )
      CALL create_directory( dirname )
      ! in the following case do not read or write polarization
      IF ( noncolin .AND. .NOT.domag ) THEN
         nspin_ = 1
      ELSE
         nspin_ = nspin
      ENDIF
      ! Write G-space density
      IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
           CALL write_rhog( TRIM(dirname) // "charge-density", &
           root_bgrp, intra_bgrp_comm, &
           bg(:,1)*tpiba, bg(:,2)*tpiba, bg(:,3)*tpiba, &
           gamma_only, mill, ig_l2g, rho%of_g(:,1:nspin_) )
      !
      ! Write kinetic energy density density
      IF ( dft_is_meta() ) THEN 
         IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
              CALL write_rhog( TRIM(dirname) // "ekin-density", &
              root_bgrp, intra_bgrp_comm, &
              bg(:,1)*tpiba, bg(:,2)*tpiba, bg(:,3)*tpiba, &
              gamma_only, mill, ig_l2g, rho%kin_g(:,1:nspin_) )
         WRITE(stdout,'(5x,"Writing meta-gga kinetic term")')
      ENDIF

      ! Then write the other terms to separate files

      IF ( lda_plus_u ) THEN
         !
         iunocc = find_free_unit ()
         IF ( ionode ) THEN
            OPEN ( UNIT=iunocc, FILE = TRIM(dirname) // 'occup.txt', &
                 FORM='formatted', STATUS='unknown' )
            IF (lda_plus_u_kind.EQ.0) THEN
               WRITE( iunocc, * , iostat = ierr) rho%ns
               IF (hub_back) WRITE( iunocc, * , iostat = ierr) rho%nsb
            ELSEIF (lda_plus_u_kind.EQ.1) THEN
               IF (noncolin) THEN
                  WRITE( iunocc, * , iostat = ierr) rho%ns_nc
               ELSE
                  WRITE( iunocc, * , iostat = ierr) rho%ns
               ENDIF
            ELSEIF (lda_plus_u_kind.EQ.2) THEN
               WRITE( iunocc, * , iostat = ierr) nsg
               ! Write Hubbard_V to file
               CALL write_V  
            ENDIF
         ENDIF
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
            OPEN ( UNIT=iunpaw, FILE = TRIM(dirname) // 'paw.txt', &
                 FORM='formatted', STATUS='unknown' )
            WRITE( iunpaw, * , iostat = ierr) rho%bec
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('write_scf', 'Writing PAW becsum',1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunpaw, STATUS = 'KEEP' )
         ENDIF
         !
      END IF

      RETURN
    END SUBROUTINE write_scf

    SUBROUTINE read_scf ( rho, nspin, gamma_only )
      !
      USE scf,              ONLY : scf_type
      USE paw_variables,    ONLY : okpaw
      USE ldaU,             ONLY : lda_plus_u, starting_ns, hub_back, &
                                   lda_plus_u_kind, nsg
      USE noncollin_module, ONLY : noncolin
      USE spin_orb,         ONLY : domag
      USE gvect,            ONLY : ig_l2g
      USE funct,            ONLY : dft_is_meta
      USE io_files,         ONLY : restart_dir
      USE io_global,        ONLY : ionode, ionode_id, stdout
      USE mp_bands,         ONLY : root_bgrp, intra_bgrp_comm
      USE mp_images,        ONLY : intra_image_comm
      USE mp,               ONLY : mp_bcast, mp_sum
      !
      IMPLICIT NONE
      TYPE(scf_type),   INTENT(INOUT)        :: rho
      INTEGER,          INTENT(IN)           :: nspin
      LOGICAL, OPTIONAL,INTENT(IN)           :: gamma_only
      !
      CHARACTER(LEN=256) :: dirname
      LOGICAL :: lexist
      INTEGER :: nspin_, iunocc, iunpaw, ierr
      INTEGER, EXTERNAL :: find_free_unit

      dirname = restart_dir ( )
      ! in the following case do not read or write polarization
      IF ( noncolin .AND. .NOT.domag ) THEN
         nspin_=1
      ELSE
         nspin_=nspin
      ENDIF
      ! read charge density
      CALL read_rhog( TRIM(dirname) // "charge-density", &
           root_bgrp, intra_bgrp_comm, &
           ig_l2g, nspin_, rho%of_g, gamma_only )
      IF ( nspin > nspin_) rho%of_g(:,nspin_+1:nspin) = (0.0_dp, 0.0_dp)
      !
      ! read kinetic energy density
      IF ( dft_is_meta() ) THEN
         CALL read_rhog( TRIM(dirname) // "ekin-density", &
              root_bgrp, intra_bgrp_comm, &
              ig_l2g, nspin_, rho%kin_g, gamma_only, ierr )
         IF ( ierr == 0 ) THEN
            WRITE(stdout,'(5x,"Reading meta-gga kinetic term")')
         ELSE
            rho%kin_g(:,:) = (0.0_dp, 0.0_dp)
            WRITE(stdout,'(5x,"BEWARE: kinetic-energy density file not found,",&
                    & " Kinetic-energy density set to 0")')
         ENDIF
      END IF

      IF ( lda_plus_u ) THEN
         !
         ! The occupations ns also need to be read in order to build up
         ! the potential
         !
         iunocc = find_free_unit ()
         IF ( ionode ) THEN
            OPEN ( UNIT=iunocc, FILE = TRIM(dirname) // 'occup.txt', &
                 FORM='formatted', STATUS='old', IOSTAT=ierr )
            IF (lda_plus_u_kind.EQ.0) THEN
               READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns
               IF (hub_back) READ( UNIT = iunocc, FMT = * , iostat = ierr) rho%nsb
            ELSEIF (lda_plus_u_kind.EQ.1) THEN
               IF (noncolin) THEN
                  READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns_nc
               ELSE
                  READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns
               ENDIF
            ELSEIF (lda_plus_u_kind.EQ.2) THEN
               READ( UNIT = iunocc, FMT = * , iostat = ierr) nsg 
            ENDIF
         ENDIF
         !
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('read_scf', 'Reading ldaU ns', 1)
         !
         IF ( ionode ) THEN
            CLOSE( UNIT = iunocc, STATUS = 'KEEP')
         ELSE
            IF (lda_plus_u_kind.EQ.0) THEN
               rho%ns(:,:,:,:) = 0.D0
               IF (hub_back) rho%nsb(:,:,:,:) = 0.D0
            ELSEIF (lda_plus_u_kind.EQ.1) THEN
               IF (noncolin) THEN
                  rho%ns_nc(:,:,:,:) = 0.D0
               ELSE
                  rho%ns(:,:,:,:) = 0.D0
               ENDIF 
            ELSEIF (lda_plus_u_kind.EQ.2) THEN
               nsg(:,:,:,:,:) = (0.d0, 0.d0) 
            ENDIF
         ENDIF
         !
         IF (lda_plus_u_kind.EQ.0) THEN
            CALL mp_sum(rho%ns, intra_image_comm) 
            IF (hub_back) CALL mp_sum(rho%nsb, intra_image_comm)   
         ELSEIF (lda_plus_u_kind.EQ.1) THEN
            IF (noncolin) THEN
               CALL mp_sum(rho%ns_nc, intra_image_comm)
            ELSE
               CALL mp_sum(rho%ns, intra_image_comm)
            ENDIF
         ELSEIF (lda_plus_u_kind.EQ.2) THEN
            CALL mp_sum(nsg, intra_image_comm)
         ENDIF
         !
         ! If projections on Hubbard manifold are read from file, there is no
         ! need to set starting values: reset them to prevent any problem
         starting_ns = -1.0_dp
         !
      ENDIF
      !
      IF ( okpaw ) THEN
         !
         ! Also the PAW coefficients are needed:
         !
         iunpaw = find_free_unit ()
         IF ( ionode ) THEN
            OPEN ( UNIT=iunpaw, FILE = TRIM(dirname) // 'paw.txt', &
                 FORM='formatted', STATUS='old', IOSTAT=ierr )
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
      RETURN
    END SUBROUTINE read_scf
    !
END MODULE io_rho_xml
