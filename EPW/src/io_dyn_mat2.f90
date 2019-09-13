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
  !! This module contains methods to read and write the dynamical
  !! matrix and the interatomic force constants files in xml format.
  !!
  USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                          iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                          iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                          iotk_scan_empty
  USE kinds,       ONLY : DP
  USE mp_images,   ONLY : intra_image_comm
  USE io_global,   ONLY : meta_ionode
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
  CHARACTER(iotk_attlenx) :: attr
  !
  CONTAINS
    !
    !----------------------------------------------------------------------------
    SUBROUTINE read_dyn_mat_param(fildyn, ntyp, nat)
    !----------------------------------------------------------------------------
    !! 
    !! Read paramters from the dynamical matrix
    !! 
    CHARACTER(LEN = 256), INTENT(in) :: fildyn
    !! Name of the file to read
    INTEGER, INTENT(out) :: ntyp
    !! Number of type of atoms
    INTEGER, INTENT(out) :: nat
    !! Number of atoms
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    ! 
    IF (meta_ionode) THEN
      !
      CALL iotk_free_unit(iunout, ierr)
      !
    ENDIF
    !
    CALL errore('read_dyn_mat_param', 'no free units to write ', ierr)
    IF (meta_ionode) THEN
      !
      ! Open XML descriptor
      ierr = 0
      CALL iotk_open_read(iunout, FILE = TRIM(fildyn) // '.xml', BINARY = .FALSE., IERR = ierr)
    ENDIF
    !
    CALL errore('read_dyn_mat_param', 'error opening the dyn mat file ', ierr)
    !
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iunout, "GEOMETRY_INFO")
      CALL iotk_scan_dat(iunout, "NUMBER_OF_TYPES", ntyp)
      CALL iotk_scan_dat(iunout, "NUMBER_OF_ATOMS", nat)
      CALL iotk_scan_end(iunout, "GEOMETRY_INFO")
    ENDIF
    !  
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_dyn_mat_param
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag,     &
               celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
               nqs, lrigid, epsil, zstareu, lraman, ramtns)
    !----------------------------------------------------------------------------
    !!   
    !! Read the dynamical matrix
    !!
    CHARACTER(LEN = 3), INTENT(out) :: atm(ntyp)
    !! Atom 
    LOGICAL, INTENT(out), OPTIONAL :: lrigid
    !! 
    LOGICAL, INTENT(out), OPTIONAL :: lraman
    !! Raman 
    INTEGER, INTENT(in) :: ntyp
    !! Number of type of atoms
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    INTEGER, INTENT(out) :: ibrav
    !! Bravais lattice
    INTEGER, INTENT(out) :: nspin_mag
    !! 
    INTEGER, INTENT(out) :: nqs 
    !! 
    INTEGER,  INTENT(out) :: ityp(nat)
    !! Atom type
    REAL(KIND = DP), INTENT(out) :: celldm(6)
    !! Celldm
    REAL(KIND = DP), INTENT(out) :: at(3, 3)
    !! Real-space lattice
    REAL(KIND = DP), INTENT(out) :: bg(3, 3)
    !! Reciprocal-space latrice
    REAL(KIND = DP), INTENT(out) :: omega
    !! Volume of primitive cell
    REAL(KIND = DP), INTENT(out) :: amass(ntyp)
    !! Atom mass
    REAL(KIND = DP), INTENT(out) :: tau(3, nat)
    !! Atom position
    REAL(KIND = DP), INTENT(out) :: m_loc(3, nat)
    !! 
    REAL(KIND = DP), INTENT(out), OPTIONAL :: epsil(3, 3)
    !! Dielectric cst
    REAL(KIND = DP), INTENT(out), OPTIONAL :: zstareu(3, 3, nat)
    !! 
    REAL(KIND = DP), INTENT(out), OPTIONAL :: ramtns(3, 3, 3, nat)
    !! 
    ! 
    ! Local work
    LOGICAL :: found_z
    !! 
    LOGICAL :: lrigid_
    !!  
    INTEGER :: nt
    !! Type of atoms
    INTEGER :: na
    !! Number of atoms
    INTEGER :: kc
    !! Cartesian direction
    REAL(KIND = DP) :: aux(3, 3)
    !! Auxillary
    !
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iunout, "GEOMETRY_INFO")
      CALL iotk_scan_dat(iunout, "BRAVAIS_LATTICE_INDEX", ibrav)
      CALL iotk_scan_dat(iunout, "SPIN_COMPONENTS", nspin_mag)
      CALL iotk_scan_dat(iunout, "CELL_DIMENSIONS", celldm)
      CALL iotk_scan_dat(iunout, "AT", at)
      CALL iotk_scan_dat(iunout, "BG", bg)
      CALL iotk_scan_dat(iunout, "UNIT_CELL_VOLUME_AU", omega)
      ! 
      DO nt = 1, ntyp
        CALL iotk_scan_dat(iunout, "TYPE_NAME"//TRIM(iotk_index(nt)), atm(nt))
        CALL iotk_scan_dat(iunout, "MASS" // TRIM(iotk_index(nt)), amass(nt))
      ENDDO
      DO na = 1, nat
        CALL iotk_scan_empty(iunout,"ATOM" // TRIM(iotk_index(na)), attr)
        CALL iotk_scan_attr(attr, "INDEX",  ityp(na))
        CALL iotk_scan_attr(attr, "TAU", tau(:, na))
        IF (nspin_mag == 4) THEN
          CALL iotk_scan_dat(iunout, "STARTING_MAG_" // TRIM(iotk_index(na)), m_loc(:, na))
        ENDIF 
      ENDDO
      CALL iotk_scan_dat(iunout, "NUMBER_OF_Q", nqs)

      CALL iotk_scan_end(iunout, "GEOMETRY_INFO")
      IF (PRESENT(lrigid)) lrigid = .FALSE.
      IF (PRESENT(epsil)) THEN
        CALL iotk_scan_begin(iunout, "DIELECTRIC_PROPERTIES", FOUND = lrigid_)
        IF (PRESENT(lrigid)) lrigid = lrigid_
        IF (lrigid_) THEN
          CALL iotk_scan_dat(iunout, "EPSILON", epsil)
          CALL iotk_scan_begin(iunout, "ZSTAR", FOUND = found_z)
          IF (found_z) THEN
            DO na = 1, nat
              CALL iotk_scan_dat(iunout, "Z_AT_" // TRIM(iotk_index(na)), aux(:, :))
              IF (PRESENT(zstareu)) zstareu(:, :, na) = aux
            ENDDO
            CALL iotk_scan_end(iunout, "ZSTAR")
          ELSE
            IF (PRESENT(zstareu)) zstareu = 0.0_DP
          ENDIF
          IF (PRESENT(lraman)) THEN
            CALL iotk_scan_begin(iunout, "RAMAN_TENSOR_A2", found = lraman)
            IF (lraman) THEN
              DO na = 1, nat
                DO kc = 1, 3
                  CALL iotk_scan_dat(iunout, "RAMAN_S_ALPHA" // TRIM(iotk_index(na)) &
                      // TRIM(iotk_index(kc)), aux)
                  IF (PRESENT(ramtns)) ramtns(:, :, kc, na) = aux(:, :)
                ENDDO
              ENDDO
              CALL iotk_scan_END(iunout, "RAMAN_TENSOR_A2")
            ELSE
              IF (PRESENT(ramtns)) ramtns = 0.0_DP
            ENDIF
          ENDIF
          CALL iotk_scan_end(iunout, "DIELECTRIC_PROPERTIES")
        ELSE
          IF (PRESENT(epsil)) epsil = 0.0_DP
          IF (PRESENT(zstareu)) zstareu = 0.0_DP
          IF (PRESENT(ramtns)) ramtns = 0.0_DP
        ENDIF
      ENDIF
    ENDIF
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_dyn_mat_header
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_dyn_mat(nat, iq, xq, dyn)
    !----------------------------------------------------------------------------
    !!
    !! This routine reads the dynamical matrix file. The file is assumed to
    !! be already opened. iq is the number of the dynamical matrix to read.
    !!
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    REAL(KIND = DP), INTENT(out) :: xq(3)
    !! Q-point value
    COMPLEX(KIND = DP), INTENT(out) :: dyn(3, 3, nat, nat)
    !! Dynamical matrix
    ! 
    ! Local variables
    INTEGER :: na, nb
    !! Number of atoms
    ! 
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iunout, "DYNAMICAL_MAT_" // TRIM(iotk_index(iq)))
      CALL iotk_scan_dat(iunout, "Q_POINT", xq)
      ! 
      DO na = 1, nat
        DO nb = 1, nat
          CALL iotk_scan_dat(iunout,"PHI"//TRIM(iotk_index(na)) // TRIM(iotk_index(nb)), dyn(:, :, na, nb))
        ENDDO
      ENDDO
      !  
      CALL iotk_scan_end(iunout, "DYNAMICAL_MAT_" // TRIM(iotk_index(iq)))
    ENDIF
    !  
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_dyn_mat
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_ifc_param(nr1, nr2, nr3)
    !----------------------------------------------------------------------------
    !! 
    !! Read IFC parameters 
    !! 
    INTEGER, INTENT(out) :: nr1, nr2, nr3
    !! Grid size
    ! Local varialbes
    INTEGER :: meshfft(3)
    !! Mesh
    ! 
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iunout, "INTERATOMIC_FORCE_CONSTANTS")
      CALL iotk_scan_dat(iunout, "MESH_NQ1_NQ2_NQ3", meshfft)
      nr1 = meshfft(1)
      nr2 = meshfft(2)
      nr3 = meshfft(3)
      CALL iotk_scan_end(iunout, "INTERATOMIC_FORCE_CONSTANTS")
    ENDIF
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_ifc_param
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_ifc_xml(nr1, nr2, nr3, nat, phid)
    !----------------------------------------------------------------------------
    !! 
    !! Read IFC in XML format
    !!  
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Grid size
    INTEGER, INTENT(in) :: nat
    !! Number of atoms  
    REAL(KIND = DP), INTENT(out) :: phid(nr1 * nr2 * nr3, 3, 3, nat, nat)
    !! 
    ! Local variables
    INTEGER :: na, nb
    !! Atoms
    INTEGER :: nn
    !! 
    INTEGER :: m1, m2, m3
    !! nr dimension 
    REAL(KIND = DP) :: aux(3, 3)
    !! Auxillary
    ! 
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iunout, "INTERATOMIC_FORCE_CONSTANTS")
      DO na = 1, nat
        DO nb = 1, nat
          nn = 0
          DO m3 = 1, nr3
            DO m2 = 1, nr2
              DO m1 = 1, nr1
                nn = nn + 1
                CALL iotk_scan_begin(iunout, "s_s1_m1_m2_m3" //     &
                    TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) // &
                    TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) // &
                    TRIM(iotk_index(m3)))
                CALL iotk_scan_dat(iunout, 'IFC', aux)
                phid(nn, :, :, na, nb) = aux(:, :)
                CALL iotk_scan_end(iunout, "s_s1_m1_m2_m3" //        &
                     TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) // &
                     TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) // &
                     TRIM(iotk_index(m3)))
              ENDDO ! m1
            ENDDO ! m2
          ENDDO ! m3
        ENDDO ! nb
      ENDDO ! na
      CALL iotk_scan_end(iunout, "INTERATOMIC_FORCE_CONSTANTS")
      CALL iotk_close_read(iunout)
    ENDIF ! meta_ionode
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_ifc_xml
    !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  END MODULE io_dyn_mat2
  !----------------------------------------------------------------------------
