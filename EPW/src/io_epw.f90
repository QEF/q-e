  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !                                                                            
  !----------------------------------------------------------------------
  MODULE io_epw
  !----------------------------------------------------------------------
  !! 
  !! This module contains various writing or reading routines to files. 
  !! Most of them are for restart purposes. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------
    SUBROUTINE rwepmatw(epmatw, nbnd, np, nmodes, nrec, iun, iop)
    !-----------------------------------------------------------------
    !!
    !! A simple wrapper to the davcio routine to read/write arrays
    !! instead of vectors 
    !!
    !-----------------------------------------------------------------
    USE kinds, ONLY : DP
    USE mp,    ONLY : mp_barrier
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nbnd
    !! Total number of bands
    INTEGER, INTENT(in) :: np
    !! np is either nrr_k or nq (epmatwe and epmatwp have the same structure)
    INTEGER, INTENT(in) :: nmodes
    !! Number of modes
    INTEGER, INTENT(in) :: nrec
    !! Place where to start reading/writing
    INTEGER, INTENT(in) :: iun
    !! Record number
    INTEGER, INTENT(in) :: iop
    !! If -1, read and if +1 write the matrix
    COMPLEX(KIND = DP), INTENT(inout) :: epmatw(nbnd, nbnd, np, nmodes)
    !! El-ph matrix to read or write
    !
    ! Local variables
    INTEGER :: lrec
    !! Record length
    INTEGER :: i
    !! Index number 
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: imode
    !! Mode index
    INTEGER :: ip
    !! REal space index (either nrr_k or nq)
    COMPLEX(KIND = DP):: aux(nbnd * nbnd * np * nmodes)
    !! 1-D vector to store the matrix elements. 
    !
    lrec = 2 * nbnd * nbnd * np * nmodes
    !
    IF (iop == -1) then
      !
      !  read matrix
      !
      CALL davcio(aux, lrec, iun, nrec, -1)
      !
      i = 0
      DO imode = 1, nmodes
        DO ip = 1, np
          DO jbnd = 1, nbnd
            DO ibnd = 1, nbnd
              i = i + 1
              epmatw(ibnd, jbnd, ip, imode) = aux(i)
              ! 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
    ELSEIF (iop == 1) THEN 
      !
      !  write matrix
      !
      i = 0
      DO imode = 1, nmodes
        DO ip = 1, np
          DO jbnd = 1, nbnd
            DO ibnd = 1, nbnd
              i = i + 1
              aux(i) = epmatw(ibnd, jbnd, ip, imode)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      CALL davcio(aux, lrec, iun, nrec, +1)
      !
    ELSE
      !
      CALL errore('rwepmatw', 'iop not permitted', 1)
      !
    ENDIF
    !
    !----------------------------------------------------------------------
    END SUBROUTINE rwepmatw
    !----------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------
    SUBROUTINE epw_write(nrr_k, nrr_q, nrr_g, w_centers)
    !----------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : nbndsub, vme, eig_read, etf_mem
    USE pwcom,     ONLY : ef, nelec
    USE elph2,     ONLY : chw, rdw, cdmew, cvmew, chw_ks, &
                          zstar, epsi, epmatwp
    USE ions_base, ONLY : amass, ityp, nat, tau
    USE cell_base, ONLY : at, bg, omega, alat
    USE phcom,     ONLY : nmodes  
    USE io_var,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp, &
                          crystal
    USE noncollin_module, ONLY : noncolin              
    USE io_files,  ONLY : prefix, diropn
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id, stdout
    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS vectors for the electrons
    INTEGER, INTENT(in) :: nrr_q
    !! Number of WS vectors for the phonons
    INTEGER, INTENT(in) :: nrr_g
    !! Number of WS vectors for the electron-phonons
    REAL(KIND = DP), INTENT(in) :: w_centers(3, nbndsub)
    !! Wannier center
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: filint
    !! Name of the file
    LOGICAL             :: exst
    !! The file exists
    INTEGER :: ibnd, jbnd
    !! Band index
    INTEGER :: jmode, imode
    !! Mode index        
    INTEGER :: irk, irq, irg
    !! WS vector looping index on electron, phonons and el-ph
    INTEGER :: ipol
    !! Cartesian direction (polarison direction)
    INTEGER :: lrepmatw
    !! Record length
    !
    WRITE(stdout,'(/5x,"Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep to file"/)')
    !
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = epwdata, FILE = 'epwdata.fmt')
      OPEN(UNIT = crystal, FILE = 'crystal.fmt')
      IF (vme) THEN 
        OPEN(UNIT = iunvmedata, FILE = 'vmedata.fmt')
      ELSE
        OPEN(UNIT = iundmedata, FILE = 'dmedata.fmt')
      ENDIF
      IF (eig_read) OPEN(UNIT = iunksdata, FILE = 'ksdata.fmt')
      WRITE(crystal,*) nat
      WRITE(crystal,*) nmodes
      WRITE(crystal,*) nelec
      WRITE(crystal,*) at
      WRITE(crystal,*) bg
      WRITE(crystal,*) omega
      WRITE(crystal,*) alat
      WRITE(crystal,*) tau
      WRITE(crystal,*) amass
      WRITE(crystal,*) ityp
      WRITE(crystal,*) noncolin
      WRITE(crystal,*) w_centers
      !
      WRITE(epwdata,*) ef
      WRITE(epwdata,*) nbndsub, nrr_k, nmodes, nrr_q, nrr_g
      WRITE(epwdata,*) zstar, epsi
      !
      DO ibnd = 1, nbndsub
        DO jbnd = 1, nbndsub
          DO irk = 1, nrr_k
            WRITE (epwdata,*) chw(ibnd, jbnd, irk)
            IF (eig_read) WRITE (iunksdata,*) chw_ks(ibnd, jbnd, irk)
            DO ipol = 1, 3
              IF (vme) THEN 
                WRITE(iunvmedata,*) cvmew(ipol, ibnd, jbnd, irk)
              ELSE
                WRITE(iundmedata,*) cdmew(ipol, ibnd, jbnd, irk)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      DO imode = 1, nmodes
        DO jmode = 1, nmodes
          DO irq = 1, nrr_q
            WRITE(epwdata,*) rdw(imode, jmode, irq) 
          ENDDO
        ENDDO
      ENDDO
      !
      IF (etf_mem == 0) THEN
        ! SP: The call to epmatwp is now inside the loop
        !     This is important as otherwise the lrepmatw INTEGER 
        !     could become too large for integer(kind=4).
        !     Note that in Fortran the record length has to be a integer
        !     of kind 4. 
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
        filint   = TRIM(prefix)//'.epmatwp'
        CALL diropn(iunepmatwp, 'epmatwp', lrepmatw, exst)
        DO irg = 1, nrr_g
          CALL davcio(epmatwp(:, :, :, :, irg), lrepmatw, iunepmatwp, irg, +1)
        ENDDO
        ! 
        CLOSE(iunepmatwp)
      ENDIF 
      !
      CLOSE(epwdata)
      CLOSE(crystal)
      IF (vme) THEN 
        CLOSE(iunvmedata)
      ELSE
        CLOSE(iundmedata)
      ENDIF
      IF (eig_read) CLOSE(iunksdata)
      !
    ENDIF
    !--------------------------------------------------------------------------------
    END SUBROUTINE epw_write
    !--------------------------------------------------------------------------------
    ! 
    !--------------------------------------------------------------------------------
    SUBROUTINE epw_read(nrr_k, nrr_q, nrr_g)
    !--------------------------------------------------------------------------------
    USE epwcom,    ONLY : nbndsub, vme, eig_read, etf_mem, lifc 
    USE pwcom,     ONLY : ef
    USE elph2,     ONLY : chw, rdw, epmatwp, cdmew, cvmew, chw_ks, zstar, epsi
    USE ions_base, ONLY : nat
    USE phcom,     ONLY : nmodes  
    USE io_global, ONLY : stdout
    USE io_files,  ONLY : prefix, diropn
    USE io_var,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp
    USE constants_epw, ONLY : czero, zero
#if defined(__NAG)
    USE f90_unix_io,ONLY : flush
#endif
    USE io_global, ONLY : ionode_id
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : world_comm
    USE mp_world,  ONLY : mpime
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: nrr_k
    !! Number of WS vectors for the electrons 
    INTEGER, INTENT(out) :: nrr_q
    !! Number of WS vectors for the phonons
    INTEGER, INTENT(out) :: nrr_g
    !! Number of WS vectors for the electron-phonons
    ! 
    ! Local variables
    ! 
    CHARACTER(LEN = 256) :: filint                                                                                   
    !! Name of the file
    LOGICAL :: exst                           
    !! The file exists
    INTEGER :: ibnd, jbnd
    !! Band index
    INTEGER :: jmode, imode
    !! Mode index        
    INTEGER :: irk, irq, irg
    !! WS vector looping index on electron, phonons and el-ph
    INTEGER :: ipol
    !! Cartesian direction (polarison direction)
    INTEGER :: lrepmatw
    !! Record length
    INTEGER :: ios 
    !! Status of files
    INTEGER :: ierr
    !! Error status
    ! 
    WRITE(stdout,'(/5x,"Reading Hamiltonian, Dynamical matrix and EP vertex in Wann rep from file"/)')
    FLUSH(stdout)
    ! 
    ! This is important in restart mode as zstar etc has not been allocated
    ALLOCATE(zstar(3, 3, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating zstar', 1)
    ALLOCATE(epsi(3, 3), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating epsi', 1)
    ! 
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = epwdata, FILE = 'epwdata.fmt', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore ('epw_read', 'error opening epwdata.fmt', epwdata)
      IF (eig_read) OPEN(UNIT = iunksdata, FILE = 'ksdata.fmt', STATUS = 'old', IOSTAT = ios)
      IF (eig_read .AND. ios /= 0) CALL errore ('epw_read', 'error opening ksdata.fmt', iunksdata)
      IF (vme) THEN 
        OPEN(UNIT = iunvmedata, FILE = 'vmedata.fmt', STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) CALL errore ('epw_read', 'error opening vmedata.fmt', iunvmedata)
      ELSE
        OPEN(UNIT = iundmedata, FILE = 'dmedata.fmt', STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) CALL errore ('epw_read', 'error opening dmedata.fmt', iundmedata)
      ENDIF
      READ(epwdata,*) ef
      READ(epwdata,*) nbndsub, nrr_k, nmodes, nrr_q, nrr_g
      READ(epwdata,*) zstar, epsi
      ! 
    ENDIF
    CALL mp_bcast(ef,      ionode_id, world_comm)
    CALL mp_bcast(nbndsub, ionode_id, world_comm)
    CALL mp_bcast(nrr_k,   ionode_id, world_comm)
    CALL mp_bcast(nmodes,  ionode_id, world_comm)
    CALL mp_bcast(nrr_q,   ionode_id, world_comm)
    CALL mp_bcast(nrr_g,   ionode_id, world_comm)
    CALL mp_bcast(zstar,   ionode_id, world_comm)
    CALL mp_bcast(epsi,    ionode_id, world_comm)
    !
    ALLOCATE(chw(nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating chw', 1)
    ALLOCATE(chw_ks(nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating chw_ks', 1)
    ALLOCATE(rdw(nmodes, nmodes,  nrr_q), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating rdw', 1)
    IF (vme) THEN 
      ALLOCATE(cvmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error allocating cvmew', 1)
    ELSE
      ALLOCATE(cdmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error allocating cdmew', 1)
    ENDIF
    !
    IF (mpime == ionode_id) THEN
     !
     DO ibnd = 1, nbndsub
       DO jbnd = 1, nbndsub
         DO irk = 1, nrr_k
           READ(epwdata,*) chw(ibnd, jbnd, irk)
           IF (eig_read) READ(iunksdata,*) chw_ks(ibnd, jbnd, irk)
           DO ipol = 1,3
             IF (vme) THEN 
               READ(iunvmedata,*) cvmew(ipol, ibnd, jbnd, irk)
             ELSE
               READ(iundmedata,*) cdmew(ipol, ibnd, jbnd, irk)
             ENDIF
           ENDDO
         ENDDO
       ENDDO
     ENDDO
     !
     IF (.NOT. lifc) THEN
       DO imode = 1, nmodes
         DO jmode = 1, nmodes
           DO irq = 1, nrr_q
             READ(epwdata,*) rdw(imode, jmode, irq)
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     !
    ENDIF
    !
    CALL mp_bcast(chw, ionode_id, world_comm)
    !
    IF (eig_read) CALL mp_bcast(chw_ks, ionode_id, world_comm)
    IF (.NOT. lifc) CALL mp_bcast(rdw, ionode_id, world_comm)
    !
    IF (vme) THEN 
      CALL mp_bcast(cvmew, ionode_id, world_comm)
    ELSE
      CALL mp_bcast(cdmew, ionode_id, world_comm)
    ENDIF
    !
    IF (lifc) THEN
      CALL read_ifc
    ENDIF
    !
    IF (etf_mem == 0) THEN
      ALLOCATE(epmatwp(nbndsub, nbndsub, nrr_k, nmodes, nrr_g), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error allocating epmatwp', 1)
      epmatwp = czero
      IF (mpime == ionode_id) THEN
        ! SP: The call to epmatwp is now inside the loop
        !     This is important as otherwise the lrepmatw INTEGER 
        !     could become too large for integer(kind=4).
        !     Note that in Fortran the record length has to be a integer
        !     of kind 4.      
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
        filint   = TRIM(prefix)//'.epmatwp'
        CALL diropn(iunepmatwp, 'epmatwp', lrepmatw, exst)
        DO irg = 1, nrr_g
          CALL davcio(epmatwp(:, :, :, :, irg), lrepmatw, iunepmatwp, irg, -1)
        ENDDO
        !  
        CLOSE(iunepmatwp)
      ENDIF
      !
      CALL mp_bcast(epmatwp, ionode_id, world_comm)
      !
    ENDIF
    !
    !CALL mp_barrier(inter_pool_comm)
    IF (mpime == ionode_id) THEN
      CLOSE(epwdata)
      IF (vme) THEN 
        CLOSE(iunvmedata)
      ELSE
        CLOSE(iundmedata)
      ENDIF
    ENDIF
    !
    WRITE(stdout, '(/5x,"Finished reading Wann rep data from file"/)')
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE epw_read
    !------------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE read_dyn_mat_param(fildyn, ntyp, nat)
    !----------------------------------------------------------------------------
    !! 
    !! Read paramters from the dynamical matrix
    !! 
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                            iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                            iotk_scan_empty
    USE kinds,       ONLY : DP
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn
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
    ! Local variables
    INTEGER :: ierr
    !! Error status
    ! 
    IF (meta_ionode) THEN
      !
      CALL iotk_free_unit(iudyn, ierr)
      !
    ENDIF
    !
    CALL errore('read_dyn_mat_param', 'no free units to write ', ierr)
    IF (meta_ionode) THEN
      !
      ! Open XML descriptor
      ierr = 0
      CALL iotk_open_read(iudyn, FILE = TRIM(fildyn) // '.xml', BINARY = .FALSE., IERR = ierr)
    ENDIF
    !
    CALL errore('read_dyn_mat_param', 'error opening the dyn mat file ', ierr)
    !
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "GEOMETRY_INFO")
      CALL iotk_scan_dat(iudyn, "NUMBER_OF_TYPES", ntyp)
      CALL iotk_scan_dat(iudyn, "NUMBER_OF_ATOMS", nat)
      CALL iotk_scan_end(iudyn, "GEOMETRY_INFO")
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
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                            iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                            iotk_scan_empty
    USE kinds,       ONLY : DP
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn
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
    ! 
    ! Local work
    CHARACTER(iotk_attlenx) :: attr
    !! Attribute
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
      CALL iotk_scan_begin(iudyn, "GEOMETRY_INFO")
      CALL iotk_scan_dat(iudyn, "BRAVAIS_LATTICE_INDEX", ibrav)
      CALL iotk_scan_dat(iudyn, "SPIN_COMPONENTS", nspin_mag)
      CALL iotk_scan_dat(iudyn, "CELL_DIMENSIONS", celldm)
      CALL iotk_scan_dat(iudyn, "AT", at)
      CALL iotk_scan_dat(iudyn, "BG", bg)
      CALL iotk_scan_dat(iudyn, "UNIT_CELL_VOLUME_AU", omega)
      ! 
      DO nt = 1, ntyp
        CALL iotk_scan_dat(iudyn, "TYPE_NAME"//TRIM(iotk_index(nt)), atm(nt))
        CALL iotk_scan_dat(iudyn, "MASS" // TRIM(iotk_index(nt)), amass(nt))
      ENDDO
      DO na = 1, nat
        CALL iotk_scan_empty(iudyn,"ATOM" // TRIM(iotk_index(na)), attr)
        CALL iotk_scan_attr(attr, "INDEX",  ityp(na))
        CALL iotk_scan_attr(attr, "TAU", tau(:, na))
        IF (nspin_mag == 4) THEN
          CALL iotk_scan_dat(iudyn, "STARTING_MAG_" // TRIM(iotk_index(na)), m_loc(:, na))
        ENDIF 
      ENDDO
      CALL iotk_scan_dat(iudyn, "NUMBER_OF_Q", nqs)

      CALL iotk_scan_end(iudyn, "GEOMETRY_INFO")
      IF (PRESENT(lrigid)) lrigid = .FALSE.
      IF (PRESENT(epsil)) THEN
        CALL iotk_scan_begin(iudyn, "DIELECTRIC_PROPERTIES", FOUND = lrigid_)
        IF (PRESENT(lrigid)) lrigid = lrigid_
        IF (lrigid_) THEN
          CALL iotk_scan_dat(iudyn, "EPSILON", epsil)
          CALL iotk_scan_begin(iudyn, "ZSTAR", FOUND = found_z)
          IF (found_z) THEN
            DO na = 1, nat
              CALL iotk_scan_dat(iudyn, "Z_AT_" // TRIM(iotk_index(na)), aux(:, :))
              IF (PRESENT(zstareu)) zstareu(:, :, na) = aux
            ENDDO
            CALL iotk_scan_end(iudyn, "ZSTAR")
          ELSE
            IF (PRESENT(zstareu)) zstareu = 0.0_DP
          ENDIF
          IF (PRESENT(lraman)) THEN
            CALL iotk_scan_begin(iudyn, "RAMAN_TENSOR_A2", found = lraman)
            IF (lraman) THEN
              DO na = 1, nat
                DO kc = 1, 3
                  CALL iotk_scan_dat(iudyn, "RAMAN_S_ALPHA" // TRIM(iotk_index(na)) &
                      // TRIM(iotk_index(kc)), aux)
                  IF (PRESENT(ramtns)) ramtns(:, :, kc, na) = aux(:, :)
                ENDDO
              ENDDO
              CALL iotk_scan_END(iudyn, "RAMAN_TENSOR_A2")
            ELSE
              IF (PRESENT(ramtns)) ramtns = 0.0_DP
            ENDIF
          ENDIF
          CALL iotk_scan_end(iudyn, "DIELECTRIC_PROPERTIES")
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
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                            iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                            iotk_scan_empty
    USE kinds,       ONLY : DP
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn
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
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "DYNAMICAL_MAT_" // TRIM(iotk_index(iq)))
      CALL iotk_scan_dat(iudyn, "Q_POINT", xq)
      ! 
      DO na = 1, nat
        DO nb = 1, nat
          CALL iotk_scan_dat(iudyn, "PHI"//TRIM(iotk_index(na)) // TRIM(iotk_index(nb)), dyn(:, :, na, nb))
        ENDDO
      ENDDO
      !  
      CALL iotk_scan_end(iudyn, "DYNAMICAL_MAT_" // TRIM(iotk_index(iq)))
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
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end
    USE kinds,       ONLY : DP
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: nr1, nr2, nr3
    !! Grid size
    ! Local varialbes
    INTEGER :: meshfft(3)
    !! Mesh
    ! 
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "INTERATOMIC_FORCE_CONSTANTS")
      CALL iotk_scan_dat(iudyn, "MESH_NQ1_NQ2_NQ3", meshfft)
      nr1 = meshfft(1)
      nr2 = meshfft(2)
      nr3 = meshfft(3)
      CALL iotk_scan_end(iudyn, "INTERATOMIC_FORCE_CONSTANTS")
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
    USE iotk_module, ONLY : iotk_index, iotk_scan_begin, iotk_open_read,     &
                            iotk_attlenx, iotk_scan_dat, iotk_scan_end,      &
                            iotk_scan_attr, iotk_free_unit, iotk_close_read, &
                            iotk_scan_empty
    USE kinds,       ONLY : DP
    USE io_global,   ONLY : meta_ionode
    USE io_var,      ONLY : iudyn 
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
    INTEGER :: nn
    !! 
    INTEGER :: m1, m2, m3
    !! nr dimension 
    REAL(KIND = DP) :: aux(3, 3)
    !! Auxillary
    ! 
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iudyn, "INTERATOMIC_FORCE_CONSTANTS")
      DO na = 1, nat
        DO nb = 1, nat
          nn = 0
          DO m3 = 1, nr3
            DO m2 = 1, nr2
              DO m1 = 1, nr1
                nn = nn + 1
                CALL iotk_scan_begin(iudyn, "s_s1_m1_m2_m3" //     &
                    TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) // &
                    TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) // &
                    TRIM(iotk_index(m3)))
                CALL iotk_scan_dat(iudyn, 'IFC', aux)
                phid(nn, :, :, na, nb) = aux(:, :)
                CALL iotk_scan_end(iudyn, "s_s1_m1_m2_m3" //        &
                     TRIM(iotk_index(na)) // TRIM(iotk_index(nb)) // &
                     TRIM(iotk_index(m1)) // TRIM(iotk_index(m2)) // &
                     TRIM(iotk_index(m3)))
              ENDDO ! m1
            ENDDO ! m2
          ENDDO ! m3
        ENDDO ! nb
      ENDDO ! na
      CALL iotk_scan_end(iudyn, "INTERATOMIC_FORCE_CONSTANTS")
      CALL iotk_close_read(iudyn)
    ENDIF ! meta_ionode
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE read_ifc_xml
    !----------------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------------------
    SUBROUTINE read_ifc()
    !---------------------------------------------------------------------------------
    !!
    !! Read IFC in real space from the file generated by q2r. 
    !! Adapted from PH/matdyn.x by C. Verdi and S. Ponce
    !! 
    !
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : ifc, zstar, epsi
    USE epwcom,    ONLY : asr_typ, dvscf_dir, nqc1, nqc2, nqc3
    USE ions_base, ONLY : nat
    USE cell_base, ONLY : ibrav, omega, at, bg, celldm, alat
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iunifc
    USE noncollin_module, ONLY : nspin_mag
    USE io_global, ONLY : ionode_id
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : intra_pool_comm, inter_pool_comm, root_pool
    USE mp_world,  ONLY : mpime, world_comm
#if defined(__NAG)
    USE f90_unix_io, ONLY : flush
#endif
    !
    IMPLICIT NONE
    !
    ! Local variables 
    LOGICAL :: lpolar_
    !! Polar flag
    LOGICAL :: has_zstar
    !! Does it has Born effective charges
    LOGICAL :: is_xml_file
    !! Is the file XML
    CHARACTER(LEN = 80) :: line
    !! 
    CHARACTER(LEN = 256) :: tempfile
    !! 
    CHARACTER(LEN = 3), ALLOCATABLE :: atm(:)
    !! 
    INTEGER :: ios
    !! 
    INTEGER :: i, j
    !! 
    INTEGER :: m1, m2, m3
    !! 
    INTEGER :: na, nb
    !! 
    INTEGER :: idum
    !! 
    INTEGER :: ibid, jbid
    !! 
    INTEGER :: nabid
    !! 
    INTEGER :: nbbid
    !! 
    INTEGER :: m1bid, m2bid, m3bid
    !! 
    INTEGER :: ntyp_
    !! 
    INTEGER :: nat_
    !! 
    INTEGER :: ibrav_
    !! 
    INTEGER :: ityp_(nat)
    !! 
    INTEGER :: nqs
    !! 
    INTEGER :: ierr
    !! Error status
    INTEGER, PARAMETER :: ntypx = 10
    !! 
    REAL(KIND = DP):: tau_(3, nat)
    !! 
    REAL(KIND = DP):: amass2(ntypx)
    !! 
    REAL(KIND = DP), ALLOCATABLE :: m_loc(:, :)
    !! 
    !
    WRITE(stdout, '(/5x,"Reading interatomic force constants"/)')
    FLUSH(stdout)
    ! 
    ! Generic name for the ifc.q2r file. If it is xml, the file will be named ifc.q2r.xml instead
    tempfile = TRIM(dvscf_dir) // 'ifc.q2r'
    ! The following function will check if the file exists in xml format
    CALL check_is_xml_file(tempfile, is_xml_file)
    ! 
    IF (mpime == ionode_id) THEN
      ! 
      IF (is_xml_file) THEN
        ! pass the 'tempfile' as the '.xml' extension is added in the next routine
        CALL read_dyn_mat_param(tempfile, ntyp_, nat_)
        ALLOCATE(m_loc(3, nat_), STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc', 'Error allocating m_loc', 1)
        ALLOCATE(atm(ntyp_), STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc', 'Error allocating atm', 1)
        CALL read_dyn_mat_header(ntyp_, nat_, ibrav, nspin_mag, &
                 celldm, at, bg, omega, atm, amass2, &
                 tau_, ityp_,  m_loc, nqs, has_zstar, epsi, zstar)
        call volume(alat, at(1, 1), at(1, 2), at(1, 3), omega)
        CALL read_ifc_param(nqc1, nqc2, nqc3)
        CALL read_ifc_xml(nqc1, nqc2, nqc3, nat_, ifc)
        DEALLOCATE(m_loc, STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc', 'Error deallocating m_loc', 1)
        DEALLOCATE(atm, STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc', 'Error deallocating atm', 1)
        ! 
      ELSE
        !
        OPEN(UNIT = iunifc, FILE = tempfile, STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) call errore ('read_ifc', 'error opening ifc.q2r', iunifc)
        !
        !  read real-space interatomic force constants
        !
        READ(iunifc,'(3i4)') ntyp_ , nat_ , ibrav_
        IF (ibrav_ == 0) THEN
          DO i = 1, 3
            READ(iunifc, *) line
          ENDDO
        ENDIF
        DO i = 1, ntyp_
          READ(iunifc, '(a)') line
        ENDDO
        DO na = 1, nat
          READ(iunifc, *) idum, idum, (tau_(j, na), j = 1, 3)
        ENDDO
        READ(iunifc, *) lpolar_
        !
        IF (lpolar_) THEN
          READ (iunifc,*) ((epsi(i, j), j = 1, 3), i = 1, 3)
          DO na = 1, nat
             READ(iunifc, *) idum
             READ(iunifc, *) ((zstar(i, j, na), j = 1, 3), i = 1, 3)
          ENDDO
          WRITE(stdout, '(5x,a)') "Read Z* and epsilon"
        ENDIF
        !
        READ (iunifc,*) idum
        !
        ifc = 0.d0
        DO i = 1, 3
          DO j = 1, 3
            DO na = 1, nat
              DO nb = 1, nat
                READ(iunifc, *) ibid, jbid, nabid, nbbid
                IF (i /= ibid .OR. j /= jbid .OR. na /= nabid .OR. nb /= nbbid)  &
                  CALL errore('read_epw', 'error in reading ifc', 1)
                READ(iunifc, *) (((m1bid, m2bid, m3bid, ifc(m1, m2, m3, i, j, na, nb), &
                           m1 = 1, nqc1), m2 = 1, nqc2), m3 = 1, nqc3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF ! noncol
    ENDIF
    !
    ! It has to be casted like this because mpi cannot cast 7 indices
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            CALL mp_bcast(ifc(:, :, :, i, j, na, nb), ionode_id, inter_pool_comm)
            CALL mp_bcast(ifc(:, :, :, i, j, na, nb), root_pool, intra_pool_comm)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_bcast(zstar, ionode_id, world_comm)
    CALL mp_bcast(epsi, ionode_id, world_comm)
    CALL mp_bcast(tau_, ionode_id, world_comm)
    CALL mp_bcast(ibrav_, ionode_id, world_comm)
    !
    WRITE(stdout,'(5x,"IFC last ", 1f12.7)') ifc(nqc1, nqc2, nqc3, 3, 3, nat, nat)
    !
    CALL set_asr2(asr_typ, nqc1, nqc2, nqc3, ifc, zstar, nat, ibrav_, tau_)
    !
    !CALL mp_barrier(inter_pool_comm)
    IF (mpime == ionode_id) THEN
      CLOSE(iunifc)
    ENDIF
    !
    WRITE(stdout, '(/5x,"Finished reading ifcs"/)')
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE read_ifc
    !-------------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE set_asr2(asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
    !-----------------------------------------------------------------------
    !!
    !! Set the acoustic sum rule. 
    !! Taken directly from PHonon/PH/q2trans.f90
    !! It would be better to take the set_asr for /Modules/. 
    !! However they are different (frc) and to be consitent with q2r.x we take this one.
    !
    USE kinds,      ONLY : DP
    USE io_global,  ONLY : stdout
    !
    IMPLICIT NONE
    ! 
    CHARACTER(LEN = 10), INTENT(in) :: asr
    !! Acoustic sum rule
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! 
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    INTEGER, INTENT(in) :: ibrav
    !! Bravais lattice
    REAL(KIND = DP), INTENT(in) :: tau(3, nat)
    !! Atomic position
    REAL(KIND = DP), INTENT(inout) :: frc(nr1, nr2, nr3, 3, 3, nat, nat)
    !! Real-space IFC
    REAL(KIND = DP), INTENT(inout) :: zeu(3, 3, nat)
    !! Z*
    !
    ! Local variables
    TYPE vector
      REAL(KIND = DP), pointer :: vec(:, :, :, :, :, :, :)
    END TYPE vector
    !
    TYPE (vector) u(6 * 3 * nat)
    ! These are the "vectors" associated with the sum rules on force-constants
    INTEGER :: axis
    !! 
    INTEGER :: n
    !! 
    INTEGER :: i, j
    !! 
    INTEGER :: na, nb
    !! 
    INTEGER :: n1, n2, n3
    !! 
    INTEGER :: m, p, k, l, q, r
    !! 
    INTEGER :: i1, j1 
    !! 
    INTEGER :: na1
    !! 
    INTEGER :: ierr
    !! Error status
    INTEGER :: u_less(6 * 3 * nat)
    !! indices of the vectors u that are not independent to the preceding ones
    INTEGER :: n_less
    !! Number of index that are no independent
    INTEGER :: i_less
    !! temporary parameter
    INTEGER :: zeu_less(6 * 3)
    !! indices of the vectors zeu_u that are not independent to the preceding ones,
    INTEGER :: nzeu_less
    !! number of such vectors
    INTEGER :: izeu_less 
    !!  temporary parameter
    INTEGER, ALLOCATABLE :: ind_v(:, :, :)
    !! 
    REAL(KIND = DP) :: zeu_new(3, 3, nat)
    !! 
    REAL(KIND = DP) :: scal
    !! 
    REAL(KIND = DP) :: norm2
    !! 
    REAL(KIND = DP) :: sum
    !! 
    REAL(KIND = DP) :: zeu_u(6 * 3, 3, 3, nat)
    !! These are the "vectors" associated with the sum rules on effective charges
    REAL(KIND = DP) :: zeu_w(3, 3, nat)
    !! 
    REAL(KIND = DP) :: zeu_x(3, 3, nat)
    !! 
    REAL(KIND = DP), ALLOCATABLE :: w(:, :, :, :, :, :, :)
    !! temporary vectors and parameters
    REAL(KIND = DP), ALLOCATABLE :: x(:, :, :, :, :, :, :)
    !! temporary vectors and parameters
    REAL(KIND = DP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
    !
    REAL(KIND = DP), ALLOCATABLE :: v(:, :)
    !! These are the "vectors" associated with symmetry conditions, coded by
    !! indicating the positions (i.e. the seven indices) of the non-zero elements (there
    !! should be only 2 of them) and the value of that element. We do so in order
    !! to limit the amount of memory used.
    !
    ! Initialization. n is the number of sum rules to be considered (if
    ! asr/='simple')
    ! and 'axis' is the rotation axis in the case of a 1D system
    ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if
    ! axis='3')
    !
    IF ((asr /= 'simple') .AND. (asr /= 'crystal') .AND. (asr /= 'one-dim') .AND. (asr /= 'zero-dim')) THEN
      CALL errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
    ENDIF
    !
    IF (asr == 'simple') THEN
      !
      ! Simple Acoustic Sum Rule on effective charges
      !
      DO i = 1, 3
        DO j = 1, 3
          sum = 0.0d0
          DO na = 1, nat
            sum = sum + zeu(i, j, na)
          ENDDO
          DO na = 1, nat
            zeu(i, j, na) = zeu(i, j, na) - sum / nat
          ENDDO
        ENDDO
      ENDDO
      !
      ! Simple Acoustic Sum Rule on force constants in real space
      !
      DO i = 1, 3
        DO j = 1, 3
          DO na = 1, nat
            sum = 0.0d0
            DO nb = 1, nat
              DO n1 = 1, nr1
                DO n2 = 1, nr2
                  DO n3 = 1, nr3
                    sum = sum + frc(n1, n2, n3, i, j, na, nb)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            frc(1, 1, 1, i, j, na, na) = frc(1, 1, 1, i, j, na, na) - sum
          ENDDO
        ENDDO
      ENDDO
      WRITE (stdout,'(5x,a)') " Imposed simple ASR"
      RETURN
      !
    ENDIF ! simple
    ! 
    IF (asr == 'crystal') n = 3
    IF (asr == 'one-dim') THEN
      ! the direction of periodicity is the rotation axis
      ! It will work only if the crystal axis considered is one of
      ! the cartesian axis (typically, ibrav = 1, 6 or 8, or 4 along the z-direction)
      IF (nr1 * nr2 * nr3 == 1) axis = 3
      IF ((nr1 /= 1) .AND. (nr2 * nr3 == 1)) axis = 1
      IF ((nr2 /= 1) .AND. (nr1 * nr3 == 1)) axis = 2
      IF ((nr3 /= 1) .AND. (nr1 * nr2 == 1)) axis = 3
      IF (((nr1 /= 1) .AND. (nr2 /=1 )) .OR. ((nr2 /= 1) .AND. (nr3 /= 1)) .OR. ((nr1 /= 1) .AND. (nr3 /=1 ))) THEN
        CALL errore('set_asr', 'too many directions of periodicity in 1D system', axis)
      ENDIF
      IF ((ibrav /= 1) .AND. (ibrav /= 6) .AND. (ibrav /= 8) .AND. ((ibrav /= 4) .OR. (axis /= 3))) THEN
        WRITE(stdout, *) 'asr: rotational axis may be wrong'
      ENDIF
      WRITE(stdout, '("asr rotation axis in 1D system= ", I4)') axis
      n = 4
    ENDIF
    IF(asr == 'zero-dim') n = 6
    !
    ! Acoustic Sum Rule on effective charges
    !
    ! generating the vectors of the orthogonal of the subspace to project the effective charges matrix on
    !
    zeu_u(:, :, :, :) = 0.0d0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          zeu_new(i, j, na) = zeu(i, j, na)
        ENDDO
      ENDDO
    ENDDO
    !
    p = 0
    DO i = 1, 3
      DO j = 1, 3
         ! These are the 3*3 vectors associated with the
         ! translational acoustic sum rules
         p = p + 1
         zeu_u(p, i, j, :) = 1.0d0
         !
      ENDDO
    ENDDO
    !
    IF (n == 4) THEN
      DO i = 1, 3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p = p + 1
        DO na = 1, nat
          zeu_u(p, i, MOD(axis, 3) + 1, na) = -tau(MOD(axis + 1, 3) + 1,na)
          zeu_u(p, i, MOD(axis + 1, 3) + 1, na) = tau(MOD(axis, 3) + 1, na)
        ENDDO
      ENDDO
    ENDIF
    !
    IF (n == 6) THEN
      DO i = 1, 3
        DO j = 1, 3
          ! These are the 3*3 vectors associated with the
          ! three rotational sum rules (0D system - typ. molecule)
          p = p + 1
          DO na = 1, nat
            zeu_u(p, i, MOD(j, 3) + 1, na) = -tau(MOD(j + 1, 3) + 1, na)
            zeu_u(p, i, MOD(j + 1, 3) + 1, na) = tau(MOD(j, 3) + 1, na)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ! Gram-Schmidt orthonormalization of the set of vectors created.
    !
    nzeu_less = 0
    DO k = 1, p
      zeu_w(:, :, :) = zeu_u(k, :, :, :)
      zeu_x(:, :, :) = zeu_u(k, :, :, :)
      DO q = 1, k - 1
        r = 1
        DO izeu_less = 1, nzeu_less
          IF (zeu_less(izeu_less) == q) r = 0
        ENDDO
        IF (r /= 0) THEN
          CALL sp_zeu(zeu_x, zeu_u(q, :, :, :), nat,scal)
          zeu_w(:, :, :) = zeu_w(:, :, :) - scal * zeu_u(q, :, :, :)
        ENDIF
      ENDDO
      CALL sp_zeu(zeu_w, zeu_w, nat, norm2)
      IF (norm2 > 1.0d-16) THEN
        zeu_u(k,:,:,:) = zeu_w(:, :, :) / DSQRT(norm2)
      ELSE
        nzeu_less = nzeu_less + 1
        zeu_less(nzeu_less) = k
      ENDIF
    ENDDO
    !
    ! Projection of the effective charge "vector" on the orthogonal of the
    ! subspace of the vectors verifying the sum rules
    !
    zeu_w(:, :, :) = 0.0d0
    DO k = 1, p
      r = 1
      DO izeu_less = 1,nzeu_less
        IF (zeu_less(izeu_less) == k) r = 0
      ENDDO
      IF (r /= 0) THEN
        zeu_x(:, :, :) = zeu_u(k, :, :, :)
        CALL sp_zeu(zeu_x, zeu_new, nat, scal)
        zeu_w(:, :, :) = zeu_w(:, :, :) + scal * zeu_u(k, :, :, :)
      ENDIF
    ENDDO
    !
    ! Final substraction of the former projection to the initial zeu, to get
    ! the new "projected" zeu
    !
    zeu_new(:, :, :) = zeu_new(:, :, :) - zeu_w(:, :, :)
    CALL sp_zeu(zeu_w, zeu_w, nat, norm2)
    WRITE(stdout, '(5x,"Norm of the difference between old and new effective charges: ", 1f12.7)') DSQRT(norm2)
    !
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          zeu(i, j, na) = zeu_new(i, j, na)
        ENDDO
      ENDDO
    ENDDO
    !
    ! Acoustic Sum Rule on force constants
    !
    ! generating the vectors of the orthogonal of the subspace to project
    ! the force-constants matrix on
    !
    DO k = 1, 18 * nat
      ALLOCATE(u(k) % vec(nr1, nr2, nr3, 3, 3, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating u(k) % vec', 1)
      u(k) % vec (:, :, :, :, :, :, :) = 0.0d0
    ENDDO
    ALLOCATE(frc_new(nr1, nr2, nr3, 3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating frc_new', 1)
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = 1, nr1
              DO n2 = 1, nr2
                DO n3 = 1, nr3
                  frc_new(n1, n2, n3, i, j, na, nb) = frc(n1, n2, n3, i, j, na, nb)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    p = 0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          ! These are the 3*3*nat vectors associated with the
          ! translational acoustic sum rules
          p = p + 1
          u(p) % vec(:, :, :, i, j, na, :) = 1.0d0
          !
        ENDDO
      ENDDO
    ENDDO
    !
    IF (n == 4) THEN
      DO i = 1, 3
        DO na = 1, nat
          ! These are the 3*nat vectors associated with the
          ! single rotational sum rule (1D system)
          p = p + 1
          DO nb = 1, nat
            u(p) % vec(: ,:, :, i, MOD(axis, 3) + 1, na, nb) = -tau(MOD(axis + 1, 3) + 1, nb)
            u(p) % vec(:, :, :, i, MOD(axis + 1, 3) + 1, na, nb) = tau(MOD(axis, 3) + 1, nb)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    IF (n == 6) THEN
      DO i = 1, 3
        DO j = 1, 3
          DO na = 1, nat
            ! These are the 3*3*nat vectors associated with the
            ! three rotational sum rules (0D system - typ. molecule)
            p = p + 1
            DO nb = 1, nat
              u(p) % vec(:, :, :, i, MOD(j, 3) + 1, na, nb) = -tau(MOD(j + 1, 3) + 1, nb)
              u(p) % vec(:, :, :, i, MOD(j + 1, 3) + 1, na, nb) = tau(MOD(j, 3) + 1, nb)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ALLOCATE(ind_v(9 * nat * nat * nr1 * nr2 * nr3, 2, 7), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating ind_v', 1)
    ALLOCATE(v(9 * nat * nat * nr1 * nr2 * nr3, 2), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating v', 1)
    m = 0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = 1, nr1
              DO n2 = 1, nr2
                DO n3 = 1, nr3
                  ! These are the vectors associated with the symmetry
                  ! constraints
                  q = 1
                  l = 1
                  DO WHILE((l <= m) .AND. (q /= 0))
                    IF ((ind_v(l, 1, 1) == n1) .AND. (ind_v(l, 1, 2) == n2) .AND. &
                        (ind_v(l, 1, 3) == n3) .AND. (ind_v(l, 1, 4) == i) .AND. &
                        (ind_v(l, 1, 5) == j) .AND. (ind_v(l, 1, 6) == na) .AND. &
                        (ind_v(l, 1, 7) == nb)) q = 0
                    IF ((ind_v(l, 2, 1) == n1) .AND. (ind_v(l, 2, 2) == n2) .AND. &
                        (ind_v(l, 2, 3) == n3) .AND. (ind_v(l, 2, 4) == i) .AND. &
                        (ind_v(l, 2, 5) == j) .AND. (ind_v(l, 2, 6) == na) .AND. &
                        (ind_v(l, 2, 7) == nb)) q = 0
                    l = l + 1
                  ENDDO
                  IF ((n1 == MOD(nr1 + 1 - n1, nr1) + 1) .AND. (n2 == MOD(nr2 + 1 - n2, nr2) + 1) &
                       .AND. (n3 == MOD(nr3 + 1 - n3, nr3) + 1) .AND. (i == j) .AND. (na == nb)) q = 0
                  IF (q /= 0) THEN
                    m = m + 1
                    ind_v(m, 1, 1) = n1
                    ind_v(m, 1, 2) = n2
                    ind_v(m, 1, 3) = n3
                    ind_v(m, 1, 4) = i
                    ind_v(m, 1, 5) = j
                    ind_v(m, 1, 6) = na
                    ind_v(m, 1, 7) = nb
                    v(m, 1) = 1.0d0 / DSQRT(2.0d0)
                    ind_v(m, 2, 1) = MOD(nr1 + 1 - n1, nr1) + 1
                    ind_v(m, 2, 2) = MOD(nr2 + 1 - n2, nr2) + 1
                    ind_v(m, 2, 3) = MOD(nr3 + 1 - n3, nr3) + 1
                    ind_v(m, 2, 4) = j
                    ind_v(m, 2, 5) = i
                    ind_v(m, 2, 6) = nb
                    ind_v(m, 2, 7) = na
                    v(m, 2) = -1.0d0 / DSQRT(2.0d0)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    ! Gram-Schmidt orthonormalization of the set of vectors created.
    ! Note that the vectors corresponding to symmetry constraints are already
    ! orthonormalized by construction.
    !
    n_less = 0
    ALLOCATE(w(nr1, nr2, nr3, 3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating w', 1)
    ALLOCATE(x(nr1, nr2, nr3, 3, 3, nat, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error allocating x', 1)
    DO k = 1, p
      w(:, :, :, :, :, :, :) = u(k) % vec(:, :, :, :, :, :, :)
      x(:, :, :, :, :, :, :) = u(k) % vec(:, :, :, :, :, :, :)
      DO l = 1, m
        !
        CALL sp2(x, v(l, :), ind_v(l, :, :), nr1, nr2, nr3, nat, scal)
        DO r = 1, 2
          n1 = ind_v(l, r, 1)
          n2 = ind_v(l, r, 2)
          n3 = ind_v(l, r, 3)
          i = ind_v(l, r, 4)
          j = ind_v(l, r, 5)
          na = ind_v(l, r, 6)
          nb = ind_v(l, r, 7)
          w(n1, n2, n3, i, j, na, nb) = w(n1, n2, n3, i, j, na, nb) - scal * v(l, r)
        ENDDO
      ENDDO
      IF (k <= (9 * nat)) THEN
        na1 = MOD(k, nat)
        IF (na1 == 0) na1 = nat
        j1 = MOD((k - na1) / nat, 3) + 1
        i1 = MOD((((k - na1) / nat) - j1 + 1) / 3, 3) + 1
      ELSE
        q = k - 9 * nat
        IF (n == 4) THEN
          na1 = MOD(q, nat)
          IF (na1 == 0) na1 = nat
          i1 = MOD((q - na1) / nat, 3) + 1
        ELSE
          na1 = MOD(q, nat)
          IF (na1 == 0) na1 = nat
          j1 = MOD((q - na1) / nat, 3) + 1
          i1 = MOD((((q - na1) / nat) - j1 + 1) / 3, 3) + 1
        ENDIF
      ENDIF
      DO q = 1, k - 1
        r = 1
        DO i_less = 1, n_less
          IF (u_less(i_less) == q) r = 0
        ENDDO
        IF (r /= 0) THEN
          CALL sp3(x, u(q) % vec(:, :, :, :, :, :, :), i1, na1, nr1, nr2, nr3, nat, scal)
          w(:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) - scal * u(q) % vec(:, :, :, :, :, :, :)
        ENDIF
      ENDDO
      CALL sp1(w, w, nr1, nr2, nr3, nat, norm2)
      IF (norm2 > 1.0d-16) THEN
        u(k) % vec(:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) / DSQRT(norm2)
      ELSE
        n_less = n_less + 1
        u_less(n_less) = k
      ENDIF
    ENDDO
    !
    ! Projection of the force-constants "vector" on the orthogonal of the
    ! subspace of the vectors verifying the sum rules and symmetry contraints
    !
    w(:, :, :, :, :, :, :) = 0.0d0
    DO l = 1, m
      CALL sp2(frc_new, v(l, :), ind_v(l, :, :), nr1, nr2, nr3, nat, scal)
      DO r = 1, 2
        n1 = ind_v(l, r, 1)
        n2 = ind_v(l, r, 2)
        n3 = ind_v(l, r, 3)
        i = ind_v(l, r, 4)
        j = ind_v(l, r, 5)
        na = ind_v(l, r, 6)
        nb = ind_v(l, r, 7)
        w(n1, n2, n3, i, j, na, nb) = w(n1, n2, n3, i, j, na, nb) + scal * v(l, r)
      ENDDO
    ENDDO
    DO k = 1, p
      r = 1
      DO i_less = 1,n_less
        IF (u_less(i_less) == k) r = 0
      ENDDO
      IF (r /= 0) THEN
        x(:, :, :, :, :, :, :) = u(k) % vec(:, :, :, :, :, :, :)
        call sp1(x, frc_new, nr1, nr2, nr3, nat, scal)
        w(:, :, :, :, :, :, :) = w(:, :, :, :, :, :, :) + scal * u(k)%vec(:, :, :, :, :, :, :)
      ENDIF
      DEALLOCATE(u(k)%vec, STAT = ierr)
      IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating u(k)%vec', 1)
    ENDDO
    !
    ! Final substraction of the former projection to the initial frc, to get
    ! the new "projected" frc
    !
    frc_new(:, :, :, :, :, :, :) = frc_new(:, :, :, :, :, :, :) - w(:, :, :, :, :, :, :)
    CALL sp1(w, w, nr1, nr2, nr3, nat, norm2)
    ! 
    WRITE(stdout, '(5x,"Norm of the difference between old and new force-constants: ", 1f12.7)') DSQRT(norm2)
    !
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = 1, nr1
              DO n2 = 1, nr2
                DO n3 = 1, nr3
                  frc(n1, n2, n3, i, j, na, nb) = frc_new(n1, n2, n3, i, j, na, nb)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(x, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating x', 1)
    DEALLOCATE(w, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating w', 1)
    DEALLOCATE(v, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating v', 1)
    DEALLOCATE(ind_v, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating ind_v', 1)
    DEALLOCATE(frc_new, STAT = ierr)
    IF (ierr /= 0) CALL errore('set_asr2', 'Error deallocating frc_new', 1)
    WRITE(stdout, '(5x,a)') "Imposed crystal ASR"
    !
    RETURN
    !----------------------------------------------------------------------
    END SUBROUTINE set_asr2
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE sp_zeu(zeu_u, zeu_v, nat, scal)
    !-----------------------------------------------------------------------
    !!
    !! Does the scalar product of two effective charges matrices zeu_u and zeu_v
    !! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
    !!
    USE kinds, ONLY: DP
    ! 
    IMPLICIT NONE
    !  
    INTEGER, INTENT(in) :: nat
    !! 
    REAL(KIND = DP), INTENT(in) :: zeu_u(3, 3, nat)
    !! 
    REAL(KIND = DP), INTENT(in) :: zeu_v(3, 3, nat)
    !! 
    REAL(KIND = DP), INTENT(inout) :: scal
    !! 
    ! Local variables
    INTEGER :: i
    !! 
    INTEGER :: j
    !! 
    INTEGER :: na
    !! 
    !
    scal = 0.0d0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          scal = scal + zeu_u(i, j, na) * zeu_v(i, j, na)
        ENDDO
      ENDDO
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sp_zeu
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE sp1(u, v, nr1, nr2, nr3, nat, scal)
    !-----------------------------------------------------------------------
    !!
    !! Does the scalar product of two force-constants matrices u and v (considered as
    !! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
    !!
    USE kinds, ONLY: DP
    ! 
    IMPLICIT NONE 
    ! 
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Supercell dims
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    REAL(KIND = DP), INTENT(in) :: u(nr1, nr2, nr3, 3, 3, nat, nat)
    !! First force-constent matrix 
    REAL(KIND = DP), INTENT(in) :: v(nr1, nr2, nr3, 3, 3, nat, nat)
    !! Second force-constent matrix
    REAL(KIND = DP), INTENT(out) :: scal
    !! Scalar product
    ! 
    ! Local variables
    INTEGER ::  i, j 
    !! Cartesian direction
    INTEGER :: na, nb
    !! Atoms index
    INTEGER :: n1, n2, n3
    !! Supercell index
    !
    scal = 0.0d0
    DO i = 1, 3
      DO j = 1, 3
        DO na = 1, nat
          DO nb = 1, nat
            DO n1 = 1, nr1
              DO n2 = 1, nr2
                DO n3 = 1, nr3
                  scal = scal + u(n1, n2, n3, i, j, na, nb) * v(n1, n2, n3, i, j, na, nb)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    RETURN
    !
    !----------------------------------------------------------------------
    END SUBROUTINE sp1
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE sp2(u, v, ind_v, nr1, nr2, nr3, nat, scal)
    !-----------------------------------------------------------------------
    !!
    !! Does the scalar product of two force-constants matrices u and v (considered as
    !! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
    !! but v is coded as explained when defining the vectors corresponding to the symmetry constraints
    !!
    USE kinds, ONLY: DP
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Supercell dims
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    INTEGER, INTENT(in) :: ind_v(2, 7)
    !! Index vector
    REAL(KIND = DP), INTENT(in) :: u(nr1, nr2, nr3, 3, 3, nat, nat)
    !! Input vector
    REAL(KIND = DP), INTENT(in) :: v(2)
    !! input vector
    REAL(KIND = DP), INTENT(out) :: scal
    !! Output vector
    ! 
    ! Local variables
    INTEGER ::  i
    !! Index
    !
    scal = 0.0d0
    DO i = 1, 2
      scal = scal + u(ind_v(i, 1), ind_v(i, 2), ind_v(i, 3), ind_v(i, 4), ind_v(i, 5), ind_v(i, 6), &
           ind_v(i, 7)) * v(i)
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sp2
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE sp3(u, v, i, na, nr1, nr2, nr3, nat, scal)
    !-----------------------------------------------------------------------
    !
    ! like sp1, but in the particular case when u is one of the u(k)%vec
    ! defined in set_asr (before orthonormalization). In this case most of the
    ! terms are zero (the ones that are not are characterized by i and na), so
    ! that a lot of computer time can be saved (during Gram-Schmidt).
    !
    USE kinds, ONLY: DP
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nr1, nr2, nr3
    !! Supercell dims
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    REAL(KIND = DP), INTENT(in) :: u(nr1, nr2, nr3, 3, 3, nat, nat)
    !! First force-constent matrix 
    REAL(KIND = DP), INTENT(in) :: v(nr1, nr2, nr3, 3, 3, nat, nat)
    !! Second force-constent matrix
    REAL(KIND = DP), INTENT(out) :: scal
    !! Scalar product
    !
    ! Local variables
    INTEGER ::  i, j
    !! Cartesian direction
    INTEGER :: na, nb
    !! Atoms index
    INTEGER :: n1, n2, n3
    !! Supercell index
    !
    scal = 0.0d0
    DO j = 1,3
      DO nb = 1, nat
        DO n1 = 1, nr1
          DO n2 = 1, nr2
            DO n3 = 1, nr3
              scal = scal + u(n1, n2, n3, i, j, na, nb) * v(n1, n2, n3, i, j, na, nb)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    RETURN
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE sp3
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE check_is_xml_file(filename, is_xml_file)
    !-------------------------------------------------------------------------------
    !!
    !! This routine checks if a file is formatted in XML. It does so by
    !! checking if the file exists and if the file + '.xml' in its name exists.
    !! If both of them or none of them exists, an error is raised. If only one of
    !! them exists, it sets the 'is_xml_file' to .TRUE. of .FALSE. depending of
    !! the file found.
    !!
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256), INTENT(in) :: filename
    !! The name of the file to check if formatted in XML format
    !! This string is assumed to be trimmed
    LOGICAL, INTENT(out) :: is_xml_file
    !! Is .TRUE. if the file is in xml format. .FALSE. otherwise.
    !
    ! Local variables
    CHARACTER(LEN = 256) :: filename_xml
    !! File name
    CHARACTER(LEN = 1024) :: errmsg
    !! Error message 
    LOGICAL :: is_plain_text_file
    !! Plain tex 
    ! 
    filename_xml = TRIM(filename) // '.xml'
    filename_xml = TRIM(filename_xml)
    INQUIRE(FILE = filename, EXIST = is_plain_text_file)
    INQUIRE(FILE = filename_xml, EXIST = is_xml_file)
    ! Tell user if any inconsistencies
    IF (is_xml_file .AND. is_plain_text_file) THEN
      ! 2 different type of files exist => warn user
      errmsg = "Detected both: '" // filename // "' and '" // filename_xml // &
              &"' which one to choose?"
      CALL errore('check_is_xml_file', errmsg, 1)
    ELSEIF (.NOT. is_xml_file .AND. .NOT. is_plain_text_file) THEN
      errmsg = "Expected a file named either '" // filename //"' or '"&
              &// filename_xml // "' but none was found."
      CALL errore('check_is_xml_file', errmsg, 1)
    ENDIF
    ! else one of the file in exists
    !------------------------------------------------------------------------------
    END SUBROUTINE check_is_xml_file
    !------------------------------------------------------------------------------
    ! 
    !-------------------------------------------------------------
    SUBROUTINE readdvscf(dvscf, recn, iq, nqc)
    !-------------------------------------------------------------
    !!
    !! Open dvscf files as direct access, read, and close again
    !!
    !! SP - Nov 2017
    !! Replaced fstat by Fortran instric function inquire. 
    !! 
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !-------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_files,  ONLY : prefix
    USE units_ph,  ONLY : lrdrho
    USE fft_base,  ONLY : dfftp
    USE epwcom,    ONLY : dvscf_dir
    USE io_var,    ONLY : iudvscf
    USE low_lvl,   ONLY : set_ndnmbr 
    USE noncollin_module, ONLY : nspin_mag
    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: recn
    !! perturbation number
    INTEGER, INTENT(in) :: iq
    !! the current q-point
    INTEGER, INTENT(in) :: nqc
    !! the total number of q-points in the list
    COMPLEX(KIND = DP), INTENT(inout) :: dvscf(dfftp%nnr, nspin_mag) 
    !! dVscf potential is read from file
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: tempfile
    !! Temp file 
    CHARACTER(LEN = 3) :: filelab
    !! File number
    INTEGER :: unf_recl
    !! Rcl unit
    INTEGER :: ios
    !! Error number
    INTEGER(KIND = 8) :: mult_unit
    !! Record length
    INTEGER(KIND = 8) :: file_size
    !! File size
    REAL(KIND = DP) :: dummy
    !! Dummy variable 
    !
    !  the call to set_ndnmbr is just a trick to get quickly
    !  a file label by exploiting an existing subroutine
    !  (if you look at the sub you will find that the original 
    !  purpose was for pools and nodes)
    !   
    CALL set_ndnmbr(0, iq, 1, nqc, filelab)
    tempfile = TRIM(dvscf_dir) // TRIM(prefix) // '.dvscf_q' // filelab
    INQUIRE(IOLENGTH = unf_recl) dummy 
    unf_recl = unf_recl  * lrdrho
    mult_unit = unf_recl
    mult_unit = recn * mult_unit
    !
    !  open the dvscf file, read and close
    !
    OPEN(iudvscf, FILE = tempfile, FORM = 'unformatted', &
         ACCESS = 'direct', IOSTAT = ios, RECL = unf_recl, STATUS = 'old')
    IF (ios /= 0) CALL errore('readdvscf', 'error opening ' // tempfile, iudvscf)
    !
    ! check that the binary file is long enough
    INQUIRE(FILE = tempfile, SIZE = file_size)
    IF (mult_unit > file_size) CALL errore('readdvscf', &
         TRIM(tempfile) //' too short, check ecut', iudvscf)
    !
    READ(iudvscf, REC = recn) dvscf
    CLOSE(iudvscf, STATUS = 'keep')
    !
    RETURN
    !
    !-------------------------------------------------------------
    END SUBROUTINE readdvscf
    !-------------------------------------------------------------
    ! 
    !------------------------------------------------------------
    SUBROUTINE readwfc(ipool, recn, evc0)
    !------------------------------------------------------------
    !!
    !! Open wfc files as direct access, read, and close again
    !!
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !
    USE kinds,    ONLY : DP
    USE io_files, ONLY : prefix, tmp_dir
    USE units_lr, ONLY : lrwfc, iuwfc
    USE wvfct,    ONLY : npwx
    USE pwcom,    ONLY : nbnd
    USE low_lvl,  ONLY : set_ndnmbr
    USE noncollin_module, ONLY : npol
    USE mp_global,        ONLY : nproc_pool, me_pool, npool
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: recn
    !! kpoint number
    INTEGER, INTENT(in) :: ipool
    !! poolfile number to be read (not used in serial case)
    COMPLEX(KIND = DP), INTENT(out) :: evc0(npwx * npol, nbnd)
    !! wavefunction is read from file
    !
    ! Local variables
    CHARACTER(LEN = 256) :: tempfile
    !! Temp file 
    CHARACTER(LEN = 3) :: nd_nmbr0
    !! File number
    INTEGER :: unf_recl
    !! Rcl unit
    INTEGER :: ios
    !! Error number
    REAL(KIND = DP) :: dummy
    !! Dummy variable 
    !
    ! Open the wfc file, read and close
    CALL set_ndnmbr(ipool, me_pool, nproc_pool, npool, nd_nmbr0)
    !print*,'nd_nmbr0 ',nd_nmbr0
    !
#if defined(__MPI)
    tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc' // nd_nmbr0
# else
    tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc'
#endif
    INQUIRE(IOLENGTH = unf_recl) dummy
    unf_recl = unf_recl * lrwfc
    !
    OPEN(iuwfc, FILE = tempfile, FORM = 'unformatted', ACCESS = 'direct', IOSTAT = ios, RECL = unf_recl)
    IF (ios /= 0) CALL errore('readwfc', 'error opening wfc file', iuwfc)
    READ(iuwfc, REC = recn) evc0
    CLOSE(iuwfc, STATUS = 'keep')
    !
    RETURN
    !
    !------------------------------------------------------------
    END SUBROUTINE readwfc
    !------------------------------------------------------------
    ! 
    !--------------------------------------------------------------
    SUBROUTINE readgmap(nkstot, ngxx, ng0vec, g0vec_all_r, lower_bnd) 
    !--------------------------------------------------------------
    !!
    !! read map of G vectors G -> G-G_0 for a given q point
    !! (this is used for the folding of k+q into the first BZ) 
    !!
    ! 
    USE kinds,    ONLY : DP
    USE mp_global,ONLY : inter_pool_comm, world_comm
    USE mp,       ONLY : mp_bcast, mp_max
    use io_global,ONLY : meta_ionode, meta_ionode_id
    use io_var,   ONLY : iukgmap, iukmap
    use pwcom,    ONLY : nks
    use elph2,    ONLY : shift, gmap, igk_k_all, ngk_all
    USE io_files, ONLY : prefix
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nkstot
    !! Total number of k-points
    INTEGER, INTENT(out) :: ngxx
    !! Maximum number of G-vectors over all pools
    INTEGER, INTENT(out) :: ng0vec
    !! Number of G_0 vectors
    INTEGER, INTENT(in) :: lower_bnd
    !! Lower bound for the k-parallellization
    REAL(KIND = DP), INTENT(out) :: g0vec_all_r(3, 125)
    !! G_0 vectors needed to fold the k+q grid into the k grid, cartesian coord.
    !
    ! Lork variables
    INTEGER :: ik
    !! Counter on k-points 
    INTEGER :: ik1, itmp
    !! Temporary indeces when reading kmap and kgmap files
    INTEGER :: ig0
    !! Counter on G_0 vectors
    INTEGER :: ishift
    !! Counter on G_0 vectors
    INTEGER :: ig
    !! Counter on G vectors
    INTEGER :: ios
    !! Integer variable for I/O control
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: tmp
    !! Temporary variable
    !
    ! OBSOLETE: now we read directly the igkq to get the proper ngxx
    ! read only a piece of the map to save time 
    ! the proper allocation bound would be ngxx = max(max(igkq))
    ! where the max is taken over the ig and the ik
    ! Here I use a simpler estimate: take the sphere npwx + two
    ! extra shells. This may not work for strange shapes of the
    ! reciproc latt. In this case just set ngxx = ngm_g
    !
    ! ngxx = NINT(4./3.*3.14*(2+(3.0/4.0/3.14*DBLE(npwx))**(1./3.))**3.)
    !
    ! Note that the k+q point below does not correspond to the actual (true) 
    ! k+q, but since we only need to take the max over k and k+q this
    ! does not matter
    !
    ngxx = 0
    DO ik = 1, nks
      !
      IF (MAXVAL(igk_k_all(1:ngk_all(ik + lower_bnd - 1), ik + lower_bnd - 1)) > ngxx) THEN
        ngxx = MAXVAL(igk_k_all(1:ngk_all(ik + lower_bnd - 1), ik + lower_bnd - 1))
      ENDIF
      !
    ENDDO
    !
#if defined(__MPI)
    tmp = DBLE(ngxx)
    CALL mp_max(tmp, inter_pool_comm)  
    ngxx = NINT(tmp)
#endif
    !
    IF (meta_ionode) THEN
      !
      OPEN(iukgmap, FILE = TRIM(prefix)//'.kgmap', FORM = 'formatted', STATUS = 'old', IOSTAT = ios)
      IF (ios /=0) CALL errore('readgmap', 'error opening kgmap file', iukgmap)
      !
      DO ik = 1, nkstot
        READ(iukgmap, *) ik1, shift(ik1)
      ENDDO
      READ(iukgmap, *) ng0vec
      !
      !  the following seems crazy but I make it for compatibility
      !  with versions up to 2.1.5:
      !
      !  iukgmap has been created by ../PW/set_kplusq.f90 and has
      !  the correct gmap(), but the wrong shift() (actually the
      !  shift for a specific q-point)
      !
      !  since createkmap.f90 has regenerated the shifts for the
      !  present k-point I read them again in kmap.dat. The above 
      !  'fake' reading is because the gmap appears *after* the
      !  wrong kmap.
      !
      OPEN(iukmap, FILE = TRIM(prefix)//'.kmap', FORM = 'formatted', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore('readgmap', 'error opening kmap file', iukmap)
      DO ik = 1, nkstot
        READ(iukmap,*) ik1, itmp, shift(ik1)
      ENDDO
      CLOSE(iukmap) 
      !
    ENDIF
    !
    ! first node broadcasts ng0vec to all nodes for allocation of gmap
    !
    CALL mp_bcast(ng0vec, meta_ionode_id, world_comm)
    !
    ALLOCATE(gmap(ngxx * ng0vec), STAT = ierr)
    IF (ierr /= 0) CALL errore('readgmap', 'Error allocating gmap', 1)
    !
    IF (meta_ionode) THEN
       !
      DO ig0 = 1, ng0vec
        READ(iukgmap,*) g0vec_all_r(:,ig0)
      ENDDO
      DO ig = 1, ngxx
        ! 
        ! at variance with the nscf calculation, here gmap is read as a vector,
        ! 
        READ(iukgmap,*) (gmap(ng0vec * ( ig - 1 ) + ishift), ishift = 1, ng0vec)
      ENDDO
      !
      CLOSE(iukgmap)
      !
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    !
    CALL mp_bcast(g0vec_all_r, meta_ionode_id, world_comm)
    CALL mp_bcast(shift, meta_ionode_id, world_comm)
    CALL mp_bcast(gmap, meta_ionode_id, world_comm)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE readgmap
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE openfilepw()
    !-----------------------------------------------------------------------
    !!
    !! Adapted from the code PH/openfilq - Quantum-ESPRESSO group                
    !! This routine opens the WF files necessary for the EPW
    !! calculation.
    !!
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !-----------------------------------------------------------------------
    USE io_files,         ONLY : prefix, diropn, seqopn
    USE units_lr,         ONLY : iuwfc, lrwfc
    USE wvfct,            ONLY : nbnd, npwx
    USE noncollin_module, ONLY : npol, nspin_mag
    USE units_ph,         ONLY : lrdrho
    USE fft_base,         ONLY : dfftp
    !
    IMPLICIT NONE
    !
    ! Local variables
    LOGICAL :: exst
    !! logical variable to check file existe
    !
    IF (len_TRIM(prefix) == 0) CALL errore('openfilepw', 'wrong prefix', 1)
    !
    ! The file with the wavefunctions
    !
    iuwfc = 20
    lrwfc = 2 * nbnd * npwx * npol
    CALL diropn(iuwfc, 'wfc', lrwfc, exst)
    IF (.NOT. exst) CALL errore ('openfilepw', 'file ' // TRIM(prefix) // '.wfc' // ' not found', 1)
    !
    ! file for setting unitary gauges of eigenstates
    !
    lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
    !
    RETURN
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE openfilepw
    !----------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE io_epw
  !------------------------------------------------------------------------------
