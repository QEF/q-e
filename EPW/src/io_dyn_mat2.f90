  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2010 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !----------------------------------------------------------------------------
  MODULE io_dyn_mat2
  !----------------------------------------------------------------------------
  !!
  !!     This module contains methods to read and write the dynamical
  !!     matrix and the interatomic force constants files in xml format.
  !!
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  USE mp_images, ONLY : intra_image_comm
  USE io_global, ONLY : meta_ionode
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: read_dyn_mat_param, read_dyn_mat_header, read_dyn_mat, &
            read_ifc_xml, read_ifc_param
  !
  INTEGER, PRIVATE :: iunout
  !
  CHARACTER(iotk_attlenx)  :: attr
  !
  CONTAINS
    !
    SUBROUTINE read_dyn_mat_param(fildyn, ntyp, nat )
    !! 
    !! Read paramters from the dynamical matrix
    !! 
    CHARACTER(LEN=256), INTENT(IN) :: fildyn
    !! Name of the file to read
    INTEGER, INTENT(OUT) :: ntyp
    !! Number of type of atoms
    INTEGER, INTENT(OUT) :: nat
    !! Number of atoms
    !
    ! Local variables
    INTEGER :: ierr

    IF ( meta_ionode ) THEN
       !
       CALL iotk_free_unit( iunout, ierr )
       !
    END IF
    !
    CALL errore( 'read_dyn_mat_param', 'no free units to write ', ierr )
    IF ( meta_ionode ) THEN
       !
       ! ... open XML descriptor
       !
       ierr=0
       CALL iotk_open_read( iunout, FILE = TRIM( fildyn ) // '.xml', &
                          BINARY = .FALSE., IERR = ierr )
    ENDIF
    !
    CALL errore( 'read_dyn_mat_param', 'error opening the dyn mat file ', ierr )
    !
    IF (meta_ionode) THEN
       CALL iotk_scan_begin(iunout, "GEOMETRY_INFO" )
       !
       CALL iotk_scan_dat(iunout,"NUMBER_OF_TYPES",ntyp)
       CALL iotk_scan_dat(iunout,"NUMBER_OF_ATOMS",nat)

       CALL iotk_scan_end(iunout, "GEOMETRY_INFO" )
    ENDIF

    RETURN
    END SUBROUTINE read_dyn_mat_param

    SUBROUTINE read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag,  &
               celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
               nqs, lrigid, epsil, zstareu, lraman, ramtns)
    !!   
    !! Read the dynamical matrix
    !!
    INTEGER, INTENT(IN) :: ntyp
    !! Number of type of atoms
    INTEGER, INTENT(IN) :: nat
    !! Number of atoms
    INTEGER, INTENT(OUT) :: ibrav, nspin_mag, nqs
    CHARACTER(LEN=3), INTENT(OUT) :: atm(ntyp)
    REAL(DP), INTENT(OUT) :: celldm(6)
    REAL(DP), INTENT(OUT) :: at(3,3)
    !! Real-space lattice
    REAL(DP), INTENT(OUT) :: bg(3,3)
    !! Reciprocal-space latrice
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
    IF (meta_ionode) THEN
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
    RETURN
    END SUBROUTINE read_dyn_mat_header

    SUBROUTINE read_dyn_mat(nat,iq,xq,dyn)
    !!
    !!   This routine reads the dynamical matrix file. The file is assumed to
    !!   be already opened. iq is the number of the dynamical matrix to read.
    !!
    INTEGER, INTENT(IN) :: nat, iq
    REAL(DP), INTENT(OUT) :: xq(3)
    COMPLEX(DP), INTENT(OUT) :: dyn(3,3,nat,nat)

    INTEGER :: na, nb

    IF (meta_ionode) THEN
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

    RETURN
    END SUBROUTINE read_dyn_mat

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

    IF (meta_ionode) THEN
       CALL iotk_scan_begin( iunout, "INTERATOMIC_FORCE_CONSTANTS" )
       CALL iotk_scan_dat( iunout, "MESH_NQ1_NQ2_NQ3", meshfft )
       nr1 = meshfft(1)
       nr2 = meshfft(2)
       nr3 = meshfft(3)

       CALL iotk_scan_end( iunout, "INTERATOMIC_FORCE_CONSTANTS" )
    ENDIF
    !CALL mp_bcast(nr1, ionode_id, intra_image_comm)
    !CALL mp_bcast(nr2, ionode_id, intra_image_comm)
    !CALL mp_bcast(nr3, ionode_id, intra_image_comm)
    RETURN
    END SUBROUTINE read_ifc_param


    SUBROUTINE read_ifc_xml( nr1, nr2, nr3, nat, phid)

    INTEGER, INTENT(IN) :: nr1, nr2, nr3, nat
    REAL(DP), INTENT(OUT) :: phid(nr1*nr2*nr3,3,3,nat,nat)
    INTEGER :: na, nb, nn, m1, m2, m3
    REAL(DP) :: aux(3,3)

    IF (meta_ionode) THEN
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
    !CALL mp_bcast(phid,ionode_id, intra_image_comm)

    RETURN
    END SUBROUTINE read_ifc_xml



END MODULE io_dyn_mat2
