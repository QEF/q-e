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
  USE kinds,     ONLY : DP
  USE io_global, ONLY : ionode, ionode_id
  USE mp_images, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast
  USE xmltools
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

    LOGICAL :: epsil_,raman_, zstareu_

    INTEGER :: na, nt, kc
    REAL(DP) :: aux(3,3)
    REAL (DP), PARAMETER ::   convfact = BOHR_RADIUS_ANGS**2

    IF ( ionode ) THEN
       !
       ! ... open XML descriptor
       !
       iunout = xml_openfile (TRIM( fildyn ) // '.xml' )
       !
    ENDIF
    CALL mp_bcast( iunout, ionode_id, intra_image_comm )
    !
    IF ( iunout == -1 ) CALL errore( 'write_dyn_mat_header', &
         'error opening the dyn mat file ', 1 )
    !
    IF (ionode) THEN
       !
       call add_attr( 'version','1.0')
       call add_attr( 'encoding','UTF-8')
       CALL xmlw_writetag ( 'xml', '?' )
       CALL xmlw_opentag ( 'Root' )
       !
       CALL xmlw_opentag("GEOMETRY_INFO" )
       !
       CALL xmlw_writetag ( "NUMBER_OF_TYPES", ntyp )
       CALL xmlw_writetag( "NUMBER_OF_ATOMS", nat )
       CALL xmlw_writetag( "BRAVAIS_LATTICE_INDEX", ibrav )
       CALL xmlw_writetag( "SPIN_COMPONENTS", nspin_mag )
       CALL xmlw_writetag( "CELL_DIMENSIONS", celldm )
       CALL xmlw_writetag( "AT", at )
       CALL xmlw_writetag( "BG", bg )
       CALL xmlw_writetag( "UNIT_CELL_VOLUME_AU", omega )
       DO nt=1, ntyp
          CALL xmlw_writetag( "TYPE_NAME."//i2c(nt),atm(nt))
          CALL xmlw_writetag( "MASS." // i2c(nt),amass(nt))
       ENDDO
       DO na=1,nat
          CALL add_attr( "SPECIES", atm(ityp(na)) )
          CALL add_attr( "INDEX", ityp(na) )
          CALL add_attr( "TAU", &
               r2c(tau(1,na)) //' '// r2c(tau(2,na)) //' '// r2c(tau(3,na)) )
          CALL xmlw_writetag( "ATOM." // i2c(na), '' )
          IF (nspin_mag==4) &
             CALL xmlw_writetag( "STARTING_MAG_."//i2c(na), m_loc(:,na) )
       END DO
       CALL xmlw_writetag( "NUMBER_OF_Q",nqs )
       CALL xmlw_closetag( )
       !
       epsil_=.false.
       IF (present(epsil)) epsil_=.true.
       zstareu_=.false.
       IF (present(zstareu)) zstareu_=.true.
       raman_=.false.
       IF (PRESENT(lraman)) raman_=.true.
       !
       CALL add_attr ( "epsil", epsil_)
       CALL add_attr ( "zstar", zstareu_)
       CALL add_attr ( "raman", raman_)
       CALL xmlw_opentag( "DIELECTRIC_PROPERTIES" )
       IF ( epsil_ ) THEN
          CALL xmlw_writetag( "EPSILON",epsil)
          IF ( zstareu_ ) THEN
             CALL xmlw_opentag( "ZSTAR" )
             DO na=1, nat
                CALL xmlw_writetag( "Z_AT_."//i2c(na), zstareu(:,:,na) )
             ENDDO
             CALL xmlw_closetag( )
          ENDIF
          IF ( raman_) THEN
             CALL xmlw_opentag( "RAMAN_TENSOR_A2")
             DO na = 1, nat
                DO kc = 1, 3
                   aux(:,:) = ramtns(:, :, kc, na)*omega/fpi*convfact
                   CALL xmlw_writetag( "RAMAN_S_ALPHA."//i2c(na)//'.'//i2c(kc),&
                           aux )
                ENDDO
             ENDDO
             CALL xmlw_closetag( )
          ENDIF
       ENDIF
       CALL xmlw_closetag( )
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

    CALL xmlw_opentag( "DYNAMICAL_MAT_."//i2c(iq) )

    CALL xmlw_writetag( "Q_POINT", xq )

    DO na=1, nat
       DO nb=1,nat
          CALL xmlw_writetag( "PHI."//i2c(na)//'.'//i2c(nb),phi(:,:,na,nb) )
       ENDDO
    ENDDO

    CALL xmlw_closetag( )

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

    CALL xmlw_opentag( "FREQUENCIES_THZ_CMM1" )

    DO mu=1,3*nat
       om = SIGN( SQRT( ABS(omega2(mu)) ),  omega2(mu) )
       omega(1) = om * RY_TO_THZ
       omega(2) = om * RY_TO_CMM1
       CALL xmlw_writetag( "OMEGA."//i2c(mu), omega )
       CALL xmlw_writetag( "DISPLACEMENT."//i2c(mu), u(:,mu) )
    END DO
    CALL xmlw_closetag( )
    !
    CALL xmlw_closetag( ) ! Root
    CALL xml_closefile( )
    !
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
    CALL xmlw_opentag( "INTERATOMIC_FORCE_CONSTANTS" )
    CALL xmlw_writetag( "MESH_NQ1_NQ2_NQ3", meshfft )

    DO na=1,nat
       DO nb=1,nat
          nn=0
          DO m3=1,nr3
             DO m2=1,nr2
                DO m1=1,nr1
                   nn=nn+1
                   CALL xmlw_opentag( "s_s1_m1_m2_m3." // i2c(na) //'.'// &
                      i2c(nb) //'.'// i2c(m1) //'.'// i2c(m2) //'.'// i2c(m3))
                   aux(:,:)=DBLE(phid(nn,:,:,na,nb))
                   CALL xmlw_writetag( 'IFC', aux )
                   CALL xmlw_closetag( )
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL xmlw_closetag()
    !
    CALL xmlw_closetag() ! Root
    CALL xml_closefile()
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE write_ifc
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_dyn_mat_param(fildyn, ntyp, nat)
    !----------------------------------------------------------------------------
    !!
    !! Read paramters from the dynamical matrix
    !!
    USE io_global,   ONLY : ionode
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256), INTENT(in) :: fildyn
    !! Name of the file to read
    INTEGER, INTENT(out) :: ntyp
    !! Number of type of atoms
    INTEGER, INTENT(out) :: nat
    !! Number of atoms
    !
    ! Open XML descriptor
    !
    IF (ionode) iunout = xml_openfile( TRIM(fildyn) // '.xml')
    !
    CALL mp_bcast(iunout, ionode_id, intra_image_comm)
    IF ( iunout == -1 ) &
         CALL errore('read_dyn_mat_param', 'error opening the dyn mat file ',1)
    !
    IF (ionode) THEN
      CALL xmlr_opentag( "GEOMETRY_INFO")
      CALL xmlr_readtag( "NUMBER_OF_TYPES", ntyp)
      CALL xmlr_readtag( "NUMBER_OF_ATOMS", nat)
      CALL xmlr_closetag() ! GEOMETRY_INFO
      REWIND(iunout)
    ENDIF
    ! 
    CALL mp_bcast(ntyp, ionode_id, intra_image_comm)
    CALL mp_bcast(nat, ionode_id, intra_image_comm)
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
    USE kinds,       ONLY : DP
    USE io_global,   ONLY : ionode
    USE xmltools
    !
    IMPLICIT NONE
    !
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
    CHARACTER(LEN = 80) :: dummy
    !! 
    LOGICAL :: found_z
    !!
    LOGICAL :: lrigid_
    !!
    LOGICAL :: raman_
    !! Is Raman present
    INTEGER :: nt
    !! Type of atoms
    INTEGER :: na
    !! Number of atoms
    INTEGER :: kc
    !! Cartesian direction
    INTEGER :: ierr
    REAL(KIND = DP) :: aux(3, 3)
    !! Auxiliary
    !
    IF (ionode) THEN
      CALL xmlr_opentag("GEOMETRY_INFO")
      CALL xmlr_readtag("BRAVAIS_LATTICE_INDEX", ibrav)
      CALL xmlr_readtag("SPIN_COMPONENTS", nspin_mag)
      CALL xmlr_readtag("CELL_DIMENSIONS", celldm)
      CALL xmlr_readtag("AT", at)
      CALL xmlr_readtag("BG", bg)
      CALL xmlr_readtag("UNIT_CELL_VOLUME_AU", omega)
      DO nt = 1, ntyp
        CALL xmlr_readtag("TYPE_NAME."//i2c(nt), atm(nt))
        CALL xmlr_readtag("MASS." // i2c(nt), amass(nt))
      ENDDO
      DO na = 1, nat
        CALL xmlr_readtag("ATOM." // i2c(na), dummy)
        CALL get_attr("INDEX",  ityp(na))
        CALL get_attr("TAU", dummy )
        READ(dummy,*) tau(1, na), tau(2, na), tau(3, na)
        IF (nspin_mag == 4) THEN
          CALL xmlr_readtag("STARTING_MAG_."//i2c(na), m_loc(:, na))
        ENDIF
      ENDDO
      CALL xmlr_readtag("NUMBER_OF_Q", nqs)
      CALL xmlr_closetag() ! GEOMETRY_INFO
      !
      IF (PRESENT(epsil)) THEN
        CALL xmlr_opentag("DIELECTRIC_PROPERTIES", ierr)
        IF (ierr == -1) THEN
          IF (PRESENT(lrigid))  lrigid = .false.
          IF (PRESENT(lraman))  lraman = .false.
          epsil = 0.0_dp
          IF (PRESENT(zstareu)) zstareu = 0.0_DP
          IF (PRESENT(ramtns))  ramtns = 0.0_DP
          GOTO 10
        ENDIF
        CALL get_attr("epsil", lrigid_)
        IF (PRESENT(lrigid)) lrigid = lrigid_
        CALL get_attr("zstar", found_z)
        CALL get_attr("raman", raman_)
        IF (PRESENT(lraman)) lraman = raman_
        IF (lrigid_) THEN
          CALL xmlr_readtag( "EPSILON", epsil)
          IF (found_z) THEN
            CALL xmlr_opentag( "ZSTAR" )
            DO na = 1, nat
              CALL xmlr_readtag( "Z_AT_."//i2c(na), aux(:, :))
              IF (PRESENT(zstareu)) zstareu(:, :, na) = aux
            ENDDO
            CALL xmlr_closetag() ! ZSTAR
          ELSE
            IF (PRESENT(zstareu)) zstareu = 0.0_DP
          ENDIF
          IF (raman_) THEN
            CALL xmlr_opentag("RAMAN_TENSOR_A2" )
            IF (PRESENT(ramtns)) THEN
              DO na = 1, nat
                DO kc = 1, 3
                  CALL xmlr_readtag("RAMAN_S_ALPHA."//i2c(na)//'.'//i2c(kc), aux)
                  IF (PRESENT(ramtns)) ramtns(:, :, kc, na) = aux(:, :)
                ENDDO
              ENDDO
            ELSE
              IF (PRESENT(ramtns)) ramtns = 0.0_DP
            ENDIF
            CALL xmlr_closetag() ! RAMAN_TENSOR_A2
          ENDIF
        ELSE
           IF (PRESENT(epsil)) epsil = 0.0_DP
           IF (PRESENT(zstareu)) zstareu = 0.0_DP
           IF (PRESENT(ramtns)) ramtns = 0.0_DP
        ENDIF ! lrigid
        CALL xmlr_closetag() ! DIELECTRIC_PROPERTIES
      ENDIF ! epsil 
      10 CONTINUE
    ENDIF ! ionode
    CALL mp_bcast(ibrav, ionode_id, intra_image_comm)
    CALL mp_bcast(nspin_mag, ionode_id, intra_image_comm)
    CALL mp_bcast(celldm, ionode_id, intra_image_comm)
    CALL mp_bcast(at, ionode_id, intra_image_comm)
    CALL mp_bcast(bg, ionode_id, intra_image_comm)
    CALL mp_bcast(omega, ionode_id, intra_image_comm)
    CALL mp_bcast(atm, ionode_id, intra_image_comm)
    CALL mp_bcast(amass, ionode_id, intra_image_comm)
    CALL mp_bcast(ityp, ionode_id, intra_image_comm)
    CALL mp_bcast(tau, ionode_id, intra_image_comm)
    CALL mp_bcast(m_loc, ionode_id, intra_image_comm)
    CALL mp_bcast(nqs, ionode_id, intra_image_comm)
    IF (PRESENT(lrigid)) CALL mp_bcast(lrigid, ionode_id, intra_image_comm)
    IF (PRESENT(epsil)) CALL mp_bcast(epsil, ionode_id, intra_image_comm)
    IF (PRESENT(zstareu)) CALL mp_bcast(zstareu, ionode_id, intra_image_comm)
    IF (PRESENT(lraman)) CALL mp_bcast(lraman, ionode_id, intra_image_comm)
    IF (PRESENT(ramtns)) CALL mp_bcast(ramtns, ionode_id, intra_image_comm)
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
    USE kinds,       ONLY : DP
    USE io_global,   ONLY : ionode
    !
    IMPLICIT NONE
    !
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
    IF (ionode) THEN
      CALL xmlr_opentag("DYNAMICAL_MAT_."//i2c(iq))
      CALL xmlr_readtag("Q_POINT", xq)
      DO na = 1, nat
        DO nb = 1,nat
          CALL xmlr_readtag( "PHI."//i2c(na)//'.'//i2c(nb), dyn(:, :, na, nb))
        ENDDO
      ENDDO
      CALL xmlr_closetag() ! DYNAMICAL_MAT_.
    ENDIF
    CALL mp_bcast(xq, ionode_id, intra_image_comm)
    CALL mp_bcast(dyn, ionode_id, intra_image_comm)
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_dyn_mat
    !----------------------------------------------------------------------------    
    ! 
    !----------------------------------------------------------------------------    
    SUBROUTINE read_dyn_mat_tail(nat, omega, u)
    !----------------------------------------------------------------------------    
    !!
    !! The output of the routine in a.u.
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : RY_TO_THZ
    ! 
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    REAL(KIND = DP), INTENT(out), OPTIONAL :: omega(3 * nat)
    !! Phonon freq.
    COMPLEX(KIND = DP), INTENT(out), OPTIONAL :: u(3 * nat, 3 * nat)
    !! Eigen displacement vectors
    ! 
    ! Local variables
    REAL(KIND = DP) :: omega_(2)
    !! Phonon freq
    INTEGER :: mu
    !! 
    ! 
    IF (PRESENT(u) .AND. .NOT. PRESENT(omega)) &
       CALL errore('read_dyn_mat_tail','omega must be present to read u',1)

    IF (ionode) THEN
      IF (PRESENT(omega)) THEN
        CALL xmlr_opentag("FREQUENCIES_THZ_CMM1")
        DO mu = 1, 3 * nat
          CALL xmlr_readtag("OMEGA."//i2c(mu), omega_)
          omega(mu) = omega_(1) / RY_TO_THZ
          IF (PRESENT(u)) CALL xmlr_readtag("DISPLACEMENT."//i2c(mu),u(:,mu))
        END DO
        CALL xmlr_closetag() ! FREQUENCIES_THZ_CMM1
      ENDIF
      CALL xml_closefile()
    ENDIF
    IF (PRESENT(omega)) CALL mp_bcast(omega, ionode_id, intra_image_comm)
    IF (PRESENT(u)) CALL mp_bcast(u, ionode_id, intra_image_comm)
    ! 
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_dyn_mat_tail
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE read_ifc_param(nr1, nr2, nr3)
    !----------------------------------------------------------------------------
    !!
    !! Read IFC parameters
    !!
    !! The following sequence should be used:
    !! read_dyn_mat_param
    !! read_dyn_mat_header
    !! read_ifc_param
    !! read_ifc
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: nr1, nr2, nr3
    !! Grid size
    ! Local varialbes
    INTEGER :: meshfft(3)
    !! Mesh
    IF (ionode) THEN
      CALL xmlr_opentag( "INTERATOMIC_FORCE_CONSTANTS")
      CALL xmlr_readtag( "MESH_NQ1_NQ2_NQ3", meshfft)
      nr1 = meshfft(1)
      nr2 = meshfft(2)
      nr3 = meshfft(3)
      CALL xmlr_closetag( )
    ENDIF
    CALL mp_bcast(nr1, ionode_id, intra_image_comm)
    CALL mp_bcast(nr2, ionode_id, intra_image_comm)
    CALL mp_bcast(nr3, ionode_id, intra_image_comm)
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_ifc_param
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE read_ifc(nr1, nr2, nr3, nat, phid)
    !----------------------------------------------------------------------------
    !!
    !! Read IFC in XML format
    !!
    USE kinds,       ONLY : DP
    USE io_global,   ONLY : ionode
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Grid size
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    REAL(KIND = DP), INTENT(out) :: phid(nr1 * nr2 * nr3, 3, 3, nat, nat)
    !!
    ! Local variables
    INTEGER :: na, nb
    !! Atoms
    INTEGER :: nn, ierr
    !!
    INTEGER :: m1, m2, m3
    !! nr dimension
    REAL(KIND = DP) :: aux(3, 3)
    !! Auxiliary
    ! 
    IF (ionode) THEN
      CALL xmlr_opentag( "INTERATOMIC_FORCE_CONSTANTS", ierr)
      DO na = 1, nat
        DO nb = 1, nat
          nn = 0
          DO m3 = 1, nr3
            DO m2 = 1, nr2
              DO m1 = 1, nr1
                nn = nn + 1
                CALL xmlr_opentag( "s_s1_m1_m2_m3." // i2c(na) //'.'// &
                      i2c(nb) //'.'// i2c(m1) //'.'// i2c(m2) //'.'// i2c(m3))
                CALL xmlr_readtag( 'IFC', aux)
                phid(nn, :, :, na, nb) = aux(:, :)
                CALL xmlr_closetag( )
              ENDDO ! m1
            ENDDO ! m2
          ENDDO ! m3
        ENDDO ! nb
      ENDDO ! na
      CALL xmlr_closetag( )
      CALL xmlr_closetag( ) ! Root
      CALL xml_closefile( )
    ENDIF
    CALL mp_bcast(phid, ionode_id, intra_image_comm)
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_ifc
    !----------------------------------------------------------------------------    
  !----------------------------------------------------------------------------    
  END MODULE io_dyn_mat
  !----------------------------------------------------------------------------    
