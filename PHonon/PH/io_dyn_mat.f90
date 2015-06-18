!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE io_dyn_mat
  !----------------------------------------------------------------------------
  !
  ! ... this module contains methods to read and write the dynamical
  !     matrix and the interatomic force constants files in xml format.
  !
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : ionode, ionode_id
  USE mp_images, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: write_dyn_mat_header, write_dyn_mat, write_dyn_mat_tail, &
            write_ifc, read_dyn_mat_param, read_dyn_mat_header, read_dyn_mat, &
            read_dyn_mat_tail, read_ifc, read_ifc_param
  !
  INTEGER, PRIVATE :: iunout
  !
  CHARACTER(iotk_attlenx)  :: attr
  !
  CONTAINS
    !
    SUBROUTINE write_dyn_mat_header( fildyn, ntyp, nat, ibrav, nspin_mag,  &
               celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
               nqs, epsil, zstareu, lraman, ramtns)
    !
  USE constants,     ONLY : FPI, BOHR_RADIUS_ANGS

    INTEGER, INTENT(IN) :: ntyp, nat, ibrav, nspin_mag, nqs
    CHARACTER(LEN=256), INTENT(IN) :: fildyn
    CHARACTER(LEN=3), INTENT(IN) :: atm(ntyp)
    REAL(DP), INTENT(IN) :: celldm(6)
    REAL(DP), INTENT(IN) :: at(3,3)
    REAL(DP), INTENT(IN) :: bg(3,3)
    REAL(DP), INTENT(IN) :: omega
    REAL(DP), INTENT(IN) :: amass(ntyp)
    REAL(DP), INTENT(IN) :: tau(3,nat)
    REAL(DP), INTENT(IN) :: m_loc(3,nat)
    REAL(DP), INTENT(IN), OPTIONAL :: epsil(3,3)
    REAL(DP), INTENT(IN), OPTIONAL :: zstareu(3,3,nat)
    LOGICAL, INTENT(IN), OPTIONAL  :: lraman
    REAL(DP), INTENT(IN), OPTIONAL :: ramtns(3,3,3,nat)

    INTEGER, INTENT(IN) :: ityp(nat)

    INTEGER :: ierr, na, nt, kc
    REAL(DP) :: aux(3,3)
    REAL (DP), PARAMETER ::   convfact = BOHR_RADIUS_ANGS**2


    IF ( ionode ) THEN
       !
       CALL iotk_free_unit( iunout, ierr )
       !
    END IF
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    !
    CALL errore( 'write_dyn_mat_header', 'no free units to write ', ierr )
    IF ( ionode ) THEN
       !
       ! ... open XML descriptor
       !
       ierr=0
       CALL iotk_open_write( iunout, FILE = TRIM( fildyn ) // '.xml', &
                          BINARY = .FALSE., IERR = ierr )
    ENDIF
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    !
    CALL errore( 'write_dyn_mat_header', 'error opening the dyn mat file ', ierr )
    !
    IF (ionode) THEN
       CALL iotk_write_begin(iunout, "GEOMETRY_INFO" )
       !
       CALL iotk_write_dat(iunout, "NUMBER_OF_TYPES", ntyp )
       CALL iotk_write_dat(iunout, "NUMBER_OF_ATOMS", nat )
       CALL iotk_write_dat(iunout, "BRAVAIS_LATTICE_INDEX", ibrav )
       CALL iotk_write_dat(iunout, "SPIN_COMPONENTS", nspin_mag )
       CALL iotk_write_dat(iunout, "CELL_DIMENSIONS", celldm )
       CALL iotk_write_dat(iunout, "AT", at, COLUMNS=3 )
       CALL iotk_write_dat(iunout, "BG", bg, COLUMNS=3 )
       CALL iotk_write_dat(iunout, "UNIT_CELL_VOLUME_AU", omega )
       DO nt=1, ntyp
          CALL iotk_write_dat(iunout,"TYPE_NAME"//TRIM(iotk_index(nt)),atm(nt))
          CALL iotk_write_dat(iunout,"MASS" // TRIM(iotk_index(nt)),amass(nt))
       ENDDO
       DO na=1,nat
          CALL iotk_write_attr( attr, "SPECIES", &
                           & atm( ityp(na) ), FIRST = .TRUE. )
          CALL iotk_write_attr( attr, "INDEX", ityp(na) )
          CALL iotk_write_attr( attr, "TAU", tau(:,na) )
          CALL iotk_write_empty( iunout, &
                           & "ATOM" // TRIM(iotk_index(na)), attr )
          IF (nspin_mag==4) &
             CALL iotk_write_dat(iunout,"STARTING_MAG_"//TRIM(iotk_index(na)),&
                                 m_loc(:,na),COLUMNS=3)
       END DO
       CALL iotk_write_dat(iunout,"NUMBER_OF_Q",nqs)

       CALL iotk_write_end(iunout, "GEOMETRY_INFO" )
       IF (present(epsil)) THEN
          CALL iotk_write_begin(iunout, "DIELECTRIC_PROPERTIES" )
          CALL iotk_write_dat(iunout,"EPSILON",epsil,COLUMNS=3)
          IF (present(zstareu)) THEN
             CALL iotk_write_begin(iunout, "ZSTAR" )
             DO na=1, nat
                CALL iotk_write_dat(iunout,"Z_AT_"//TRIM(iotk_index(na)),&
                                         zstareu(:,:,na),COLUMNS=3)
             ENDDO
             CALL iotk_write_end(iunout, "ZSTAR" )
          ENDIF
          IF (PRESENT(lraman)) THEN
             IF (lraman) THEN
                CALL iotk_write_begin(iunout,"RAMAN_TENSOR_A2")
                DO na = 1, nat
                   DO kc = 1, 3
                      aux(:,:) = ramtns(:, :, kc, na)*omega/fpi*convfact
                      CALL iotk_write_dat(iunout, &
                             "RAMAN_S_ALPHA"//TRIM(iotk_index(na)) &
                              // TRIM(iotk_index(kc)),aux, COLUMNS=3)
                   ENDDO
                ENDDO
                CALL iotk_write_END(iunout,"RAMAN_TENSOR_A2")
             ENDIF
          ENDIF
          CALL iotk_write_end(iunout, "DIELECTRIC_PROPERTIES" )
       ENDIF
    ENDIF
    !
    RETURN
    END SUBROUTINE write_dyn_mat_header

    SUBROUTINE write_dyn_mat(nat,iq,xq,phi)

    INTEGER, INTENT(IN) :: nat, iq
    REAL(DP), INTENT(IN) :: xq(3)
    COMPLEX(DP), INTENT(IN) :: phi(3,3,nat,nat)

    INTEGER :: na, nb

    IF (.NOT.ionode) RETURN

    CALL iotk_write_begin(iunout, "DYNAMICAL_MAT_"//TRIM(iotk_index(iq)) )

    CALL iotk_write_dat(iunout,"Q_POINT",xq,COLUMNS=3)

    DO na=1, nat
       DO nb=1,nat
          CALL iotk_write_dat(iunout,"PHI"//TRIM(iotk_index(na))&
                   &//TRIM(iotk_index(nb)),phi(:,:,na,nb),COLUMNS=1)
       ENDDO
    ENDDO

    CALL iotk_write_end(iunout, "DYNAMICAL_MAT_"//TRIM(iotk_index(iq)) )

    RETURN
    END SUBROUTINE write_dyn_mat

    SUBROUTINE write_dyn_mat_tail(nat,omega2,u)

    USE constants, ONLY : RY_TO_THZ, RY_TO_CMM1

    INTEGER, INTENT(IN) :: nat
    REAL(DP), INTENT(IN) :: omega2(3*nat)
    COMPLEX(DP), INTENT(IN) :: u(3*nat,3*nat)

    REAL(DP) :: omega(2), om
    INTEGER :: mu

    IF (.NOT. ionode) RETURN

    CALL iotk_write_begin( iunout, "FREQUENCIES_THZ_CMM1" )

    DO mu=1,3*nat
       om = SIGN( SQRT( ABS(omega2(mu)) ),  omega2(mu) )
       omega(1) = om * RY_TO_THZ
       omega(2) = om * RY_TO_CMM1
       CALL iotk_write_dat(iunout,"OMEGA"//TRIM(iotk_index(mu)),&
                                                omega, COLUMNS=2)
       CALL iotk_write_dat(iunout,"DISPLACEMENT"//TRIM(iotk_index(mu)),&
                                                u(:,mu), COLUMNS=1)
    END DO

    CALL iotk_write_end( iunout, "FREQUENCIES_THZ_CMM1" )

    CALL iotk_close_write( iunout )
    RETURN
    END SUBROUTINE write_dyn_mat_tail

    SUBROUTINE write_ifc( nr1, nr2, nr3, nat, phid)

    INTEGER, INTENT(IN) :: nr1, nr2, nr3, nat
    COMPLEX(DP), INTENT(IN) :: phid(nr1*nr2*nr3,3,3,nat,nat)

    INTEGER :: meshfft(3)
    INTEGER :: na, nb, nn, m1, m2, m3
    REAL(DP) :: aux(3,3)

    IF (.NOT.ionode) RETURN

    meshfft(1)=nr1
    meshfft(2)=nr2
    meshfft(3)=nr3
    CALL iotk_write_begin( iunout, "INTERATOMIC_FORCE_CONSTANTS" )
    CALL iotk_write_dat( iunout, "MESH_NQ1_NQ2_NQ3", meshfft, COLUMNS=3 )

    DO na=1,nat
       DO nb=1,nat
          nn=0
          DO m3=1,nr3
             DO m2=1,nr2
                DO m1=1,nr1
                   nn=nn+1
                   CALL iotk_write_begin( iunout, "s_s1_m1_m2_m3" //     &
                        TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) //  &
                        TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) //  &
                        TRIM(iotk_index(m3)) )
                   aux(:,:)=DBLE(phid(nn,:,:,na,nb))
                   CALL iotk_write_dat( iunout, 'IFC', aux, COLUMNS=3 )
                   CALL iotk_write_end( iunout, "s_s1_m1_m2_m3" //     &
                        TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) //  &
                        TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) //  &
                        TRIM(iotk_index(m3)) )
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    CALL iotk_write_end( iunout, "INTERATOMIC_FORCE_CONSTANTS" )
    CALL iotk_close_write( iunout )
    RETURN
    END SUBROUTINE write_ifc

    SUBROUTINE read_dyn_mat_param(fildyn, ntyp, nat )

    CHARACTER(LEN=256), INTENT(IN) :: fildyn
    INTEGER, INTENT(OUT) :: ntyp, nat
    INTEGER :: ierr

    IF ( ionode ) THEN
       !
       CALL iotk_free_unit( iunout, ierr )
       !
    END IF
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    !
    CALL errore( 'read_dyn_mat_param', 'no free units to write ', ierr )
    IF ( ionode ) THEN
       !
       ! ... open XML descriptor
       !
       ierr=0
       CALL iotk_open_read( iunout, FILE = TRIM( fildyn ) // '.xml', &
                          BINARY = .FALSE., IERR = ierr )
    ENDIF
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    !
    CALL errore( 'read_dyn_mat_param', 'error opening the dyn mat file ', ierr )
    !
    IF (ionode) THEN
       CALL iotk_scan_begin(iunout, "GEOMETRY_INFO" )
       !
       CALL iotk_scan_dat(iunout,"NUMBER_OF_TYPES",ntyp)
       CALL iotk_scan_dat(iunout,"NUMBER_OF_ATOMS",nat)

       CALL iotk_scan_end(iunout, "GEOMETRY_INFO" )
    ENDIF

    CALL mp_bcast(ntyp, ionode_id, intra_image_comm)
    CALL mp_bcast(nat, ionode_id, intra_image_comm)

    RETURN
    END SUBROUTINE read_dyn_mat_param

    SUBROUTINE read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag,  &
               celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
               nqs, lrigid, epsil, zstareu, lraman, ramtns)

    INTEGER, INTENT(IN) :: ntyp, nat
    INTEGER, INTENT(OUT) :: ibrav, nspin_mag, nqs
    CHARACTER(LEN=3), INTENT(OUT) :: atm(ntyp)
    REAL(DP), INTENT(OUT) :: celldm(6)
    REAL(DP), INTENT(OUT) :: at(3,3)
    REAL(DP), INTENT(OUT) :: bg(3,3)
    REAL(DP), INTENT(OUT) :: omega
    REAL(DP), INTENT(OUT) :: amass(ntyp)
    REAL(DP), INTENT(OUT) :: tau(3,nat)
    REAL(DP), INTENT(OUT) :: m_loc(3,nat)
    INTEGER,  INTENT(OUT) :: ityp(nat)
    REAL(DP), INTENT(OUT), OPTIONAL :: epsil(3,3)
    REAL(DP), INTENT(OUT), OPTIONAL :: zstareu(3,3,nat)
    LOGICAL, INTENT(OUT), OPTIONAL :: lrigid
    LOGICAL, INTENT(OUT), OPTIONAL :: lraman
    REAL(DP), INTENT(OUT), OPTIONAL :: ramtns(3,3,3,nat)

    REAL(DP) :: aux(3,3)
    INTEGER :: nt, na, kc
    LOGICAL :: found_z, lrigid_
    !
    IF (ionode) THEN
       CALL iotk_scan_begin( iunout, "GEOMETRY_INFO" )
       !
       CALL iotk_scan_dat( iunout, "BRAVAIS_LATTICE_INDEX", ibrav )
       CALL iotk_scan_dat( iunout, "SPIN_COMPONENTS", nspin_mag )
       CALL iotk_scan_dat( iunout, "CELL_DIMENSIONS", celldm )
       CALL iotk_scan_dat( iunout, "AT", at )
       CALL iotk_scan_dat( iunout, "BG", bg )
       CALL iotk_scan_dat( iunout, "UNIT_CELL_VOLUME_AU", omega )

       DO nt=1, ntyp
          CALL iotk_scan_dat(iunout,"TYPE_NAME"//TRIM(iotk_index(nt)),atm(nt))
          CALL iotk_scan_dat(iunout,"MASS" // TRIM(iotk_index(nt)),amass(nt))
       ENDDO
       DO na=1,nat
          CALL iotk_scan_empty( iunout,"ATOM" // TRIM( iotk_index(na) ), attr )
          CALL iotk_scan_attr( attr, "INDEX",  ityp(na) )
          CALL iotk_scan_attr( attr, "TAU", tau(:,na) )
          IF (nspin_mag==4) &
             CALL iotk_scan_dat(iunout,"STARTING_MAG_"//TRIM(iotk_index(na)),&
                           m_loc(:,na))

       ENDDO
       CALL iotk_scan_dat(iunout,"NUMBER_OF_Q",nqs)

       CALL iotk_scan_end(iunout, "GEOMETRY_INFO" )
       IF (PRESENT(lrigid)) lrigid=.FALSE.
       IF (PRESENT(epsil)) THEN
          CALL iotk_scan_begin(iunout, "DIELECTRIC_PROPERTIES", FOUND=lrigid_)
          IF (PRESENT(lrigid)) lrigid=lrigid_
          IF (lrigid_) THEN
             CALL iotk_scan_dat(iunout,"EPSILON",epsil)
             CALL iotk_scan_begin(iunout, "ZSTAR", FOUND=found_z )
             IF (found_z) THEN
                DO na=1, nat
                   CALL iotk_scan_dat(iunout,"Z_AT_"//TRIM(iotk_index(na)),&
                                        aux(:,:))
                   IF (PRESENT(zstareu)) zstareu(:,:,na)=aux
                ENDDO
                CALL iotk_scan_end(iunout, "ZSTAR" )
             ELSE
                IF (PRESENT(zstareu)) zstareu=0.0_DP
             ENDIF
             IF (PRESENT(lraman)) THEN
                CALL iotk_scan_begin(iunout,"RAMAN_TENSOR_A2",found=lraman)
                IF (lraman) THEN
                   DO na = 1, nat
                      DO kc = 1, 3
                         CALL iotk_scan_dat(iunout, &
                             "RAMAN_S_ALPHA"//TRIM(iotk_index(na)) &
                              // TRIM(iotk_index(kc)),aux)
                         IF (PRESENT(ramtns)) ramtns(:, :, kc, na) = aux(:,:)
                      ENDDO
                   ENDDO
                   CALL iotk_scan_END(iunout,"RAMAN_TENSOR_A2")
                ELSE
                   IF (PRESENT(ramtns)) ramtns=0.0_DP
                ENDIF
             ENDIF
             CALL iotk_scan_end(iunout, "DIELECTRIC_PROPERTIES" )
          ELSE
             IF (PRESENT(epsil)) epsil=0.0_DP
             IF (PRESENT(zstareu)) zstareu=0.0_DP
             IF (PRESENT(ramtns))  ramtns=0.0_DP
          ENDIF
       ENDIF
    ENDIF
    CALL mp_bcast(ibrav,ionode_id, intra_image_comm)
    CALL mp_bcast(nspin_mag,ionode_id, intra_image_comm)
    CALL mp_bcast(celldm,ionode_id, intra_image_comm)
    CALL mp_bcast(at,ionode_id, intra_image_comm)
    CALL mp_bcast(bg,ionode_id, intra_image_comm)
    CALL mp_bcast(omega,ionode_id, intra_image_comm)
    CALL mp_bcast(atm,ionode_id, intra_image_comm)
    CALL mp_bcast(amass,ionode_id, intra_image_comm)
    CALL mp_bcast(ityp,ionode_id, intra_image_comm)
    CALL mp_bcast(tau,ionode_id, intra_image_comm)
    CALL mp_bcast(m_loc,ionode_id, intra_image_comm)
    CALL mp_bcast(nqs,ionode_id, intra_image_comm)
    IF (PRESENT(lrigid)) CALL mp_bcast(lrigid,ionode_id, intra_image_comm)
    IF (PRESENT(epsil))  CALL mp_bcast(epsil,ionode_id, intra_image_comm)
    IF (PRESENT(zstareu)) CALL mp_bcast(zstareu,ionode_id, intra_image_comm)
    IF (PRESENT(lraman)) CALL mp_bcast(lraman,ionode_id, intra_image_comm)
    IF (PRESENT(ramtns)) CALL mp_bcast(ramtns,ionode_id, intra_image_comm)

    RETURN
    END SUBROUTINE read_dyn_mat_header

    SUBROUTINE read_dyn_mat(nat,iq,xq,dyn)
!
!   This routine reads the dynamical matrix file. The file is assumed to
!   be already opened. iq is the number of the dynamical matrix to read.
!
    INTEGER, INTENT(IN) :: nat, iq
    REAL(DP), INTENT(OUT) :: xq(3)
    COMPLEX(DP), INTENT(OUT) :: dyn(3,3,nat,nat)

    INTEGER :: na, nb

    IF (ionode) THEN
       CALL iotk_scan_begin(iunout, "DYNAMICAL_MAT_"//TRIM(iotk_index(iq)) )

       CALL iotk_scan_dat(iunout,"Q_POINT",xq)

       DO na=1, nat
          DO nb=1,nat
             CALL iotk_scan_dat(iunout,"PHI"//TRIM(iotk_index(na))&
                   &//TRIM(iotk_index(nb)),dyn(:,:,na,nb))
          ENDDO
       ENDDO

       CALL iotk_scan_end(iunout, "DYNAMICAL_MAT_"//TRIM(iotk_index(iq)) )
    ENDIF
    CALL mp_bcast(xq, ionode_id, intra_image_comm)
    CALL mp_bcast(dyn, ionode_id, intra_image_comm)

    RETURN
    END SUBROUTINE read_dyn_mat

    SUBROUTINE read_dyn_mat_tail(nat,omega,u)
!
!   The output of the routine in a.u.
!
    USE constants, ONLY : RY_TO_THZ

    INTEGER, INTENT(IN) :: nat
    REAL(DP), INTENT(OUT), OPTIONAL :: omega(3*nat)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: u(3*nat,3*nat)

    REAL(DP) :: omega_(2)
    INTEGER :: mu

    IF (PRESENT(u).AND..NOT.PRESENT(omega)) &
       CALL errore('read_dyn_mat_tail','omega must be present to read u',1)

    IF (ionode) THEN
       IF (PRESENT(omega)) THEN
          CALL iotk_scan_begin( iunout, "FREQUENCIES_THZ_CMM1" )
          DO mu=1,3*nat
             CALL iotk_scan_dat(iunout,"OMEGA"//TRIM(iotk_index(mu)), omega_)
             omega(mu)=omega_(1) / RY_TO_THZ
             IF (PRESENT(u)) CALL iotk_scan_dat(iunout, &
                             "DISPLACEMENT"//TRIM(iotk_index(mu)),u(:,mu))
          END DO
          CALL iotk_scan_end( iunout, "FREQUENCIES_THZ_CMM1" )
       ENDIF

       CALL iotk_close_read( iunout )
    END IF
    IF (PRESENT(omega)) CALL mp_bcast(omega, ionode_id, intra_image_comm)
    IF (PRESENT(u)) CALL mp_bcast(u, ionode_id, intra_image_comm)

    RETURN
    END SUBROUTINE read_dyn_mat_tail

    SUBROUTINE read_ifc_param( nr1, nr2, nr3 )
!
!   To read the interatomic force constant the following sequence should
!   be used:
!   read_dyn_mat_param
!   read_dyn_mat_header
!   read_ifc_param
!   read_ifc
!
    INTEGER, INTENT(OUT) :: nr1, nr2, nr3
    INTEGER :: meshfft(3)

    IF (ionode) THEN
       CALL iotk_scan_begin( iunout, "INTERATOMIC_FORCE_CONSTANTS" )
       CALL iotk_scan_dat( iunout, "MESH_NQ1_NQ2_NQ3", meshfft )
       nr1 = meshfft(1)
       nr2 = meshfft(2)
       nr3 = meshfft(3)

       CALL iotk_scan_end( iunout, "INTERATOMIC_FORCE_CONSTANTS" )
    ENDIF
    CALL mp_bcast(nr1, ionode_id, intra_image_comm)
    CALL mp_bcast(nr2, ionode_id, intra_image_comm)
    CALL mp_bcast(nr3, ionode_id, intra_image_comm)
    RETURN
    END SUBROUTINE read_ifc_param

    SUBROUTINE read_ifc( nr1, nr2, nr3, nat, phid)

    INTEGER, INTENT(IN) :: nr1, nr2, nr3, nat
    REAL(DP), INTENT(OUT) :: phid(nr1*nr2*nr3,3,3,nat,nat)
    INTEGER :: na, nb, nn, m1, m2, m3
    REAL(DP) :: aux(3,3)

    IF (ionode) THEN
       CALL iotk_scan_begin( iunout, "INTERATOMIC_FORCE_CONSTANTS" )

       DO na=1,nat
          DO nb=1,nat
             nn=0
             DO m3=1,nr3
                DO m2=1,nr2
                   DO m1=1,nr1
                      nn=nn+1
                      CALL iotk_scan_begin( iunout, "s_s1_m1_m2_m3" //     &
                          TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) //  &
                          TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) //  &
                          TRIM(iotk_index(m3)) )
                      CALL iotk_scan_dat( iunout, 'IFC', aux )
                      phid(nn,:,:,na,nb) = aux(:,:)
                      CALL iotk_scan_end( iunout, "s_s1_m1_m2_m3" //     &
                           TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) //  &
                           TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) //  &
                           TRIM(iotk_index(m3)) )
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       CALL iotk_scan_end( iunout, "INTERATOMIC_FORCE_CONSTANTS" )
       CALL iotk_close_read( iunout )
    ENDIF
    CALL mp_bcast(phid,ionode_id, intra_image_comm)

    RETURN
    END SUBROUTINE read_ifc

END MODULE io_dyn_mat
