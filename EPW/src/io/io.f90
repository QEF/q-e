  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !
  !----------------------------------------------------------------------
  MODULE io
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
    INTEGER :: ierr
    !! Error index
    COMPLEX(KIND = DP), ALLOCATABLE :: aux(:)
    !! 1-D vector to store the matrix elements.
    !
    ALLOCATE(aux(nbnd * nbnd * np * nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('rwepmatw', 'Error allocating aux', 1)
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
    DEALLOCATE(aux, STAT = ierr)
    IF (ierr /= 0) CALL errore('rwepmatw', 'Error deallocating aux', 1)
    !
    !----------------------------------------------------------------------
    END SUBROUTINE rwepmatw
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE epw_write(nrr_k, nrr_q, nrr_g, w_centers)
    !----------------------------------------------------------------------
    !!
    !! Routine to write files on real-space grid for fine grid interpolation
    !!
    USE kinds,     ONLY : DP
    USE input,     ONLY : nbndsub, eig_read, etf_mem, lifc, lwfpt
    USE pwcom,     ONLY : ef, nelec
    USE klist,     ONLY : degauss, ngauss
    USE global_var,ONLY : chw, rdw, cdmew, cvmew, chw_ks, zstar, epsi, &
                          epmatwp, crrw, L, do_cutoff_2D_epw, dwmatwe, cpmew, &
                          nbndskip
    USE ions_base, ONLY : amass, ityp, nat, tau
    USE cell_base, ONLY : at, bg, omega, alat
    USE modes,     ONLY : nmodes
    USE io_var,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp, &
                          crystal, iudwwe, iupmwe
    USE noncollin_module, ONLY : noncolin
    USE io_files,  ONLY : prefix, diropn, tmp_dir
    USE ep_constants,      ONLY : czero
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
    INTEGER*8 :: unf_recl
    !! Record length
    INTEGER :: direct_io_factor
    !! Type of IOlength
    INTEGER :: ierr
    !! Error index
    REAL(KIND = DP) :: dummy
    !! Dummy variable
    !
    WRITE(stdout,'(/5x,"Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep to file"/)')
    !
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = epwdata, FILE = 'epwdata.fmt')
      OPEN(UNIT = crystal, FILE = 'crystal.fmt')
      OPEN(UNIT = iunvmedata, FILE = 'vmedata.fmt')
      OPEN(UNIT = iundmedata, FILE = 'dmedata.fmt')
      IF (eig_read) OPEN(UNIT = iunksdata, FILE = 'ksdata.fmt')
      WRITE(crystal,*) nat
      WRITE(crystal,*) nmodes
      WRITE(crystal,*) nelec, nbndskip
      WRITE(crystal,*) at
      WRITE(crystal,*) bg
      WRITE(crystal,*) omega
      WRITE(crystal,*) alat
      WRITE(crystal,*) tau
      WRITE(crystal,*) amass
      WRITE(crystal,*) ityp
      WRITE(crystal,*) noncolin
      WRITE(crystal,*) do_cutoff_2D_epw
      WRITE(crystal,*) w_centers
      WRITE(crystal,*) L
      WRITE(crystal,*) degauss
      WRITE(crystal,*) ngauss
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
              WRITE(iunvmedata,*) cvmew(ipol, ibnd, jbnd, irk)
              WRITE(iunvmedata,*) crrw(ipol, ibnd, jbnd, irk)
              !
              WRITE(iundmedata,*) cdmew(ipol, ibnd, jbnd, irk)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      IF (.NOT. lifc) THEN
        DO imode = 1, nmodes
          DO jmode = 1, nmodes
            DO irq = 1, nrr_q
              WRITE(epwdata,*) rdw(imode, jmode, irq)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        ! output a dummy value if lifc=true.
        DO imode = 1, nmodes
          DO jmode = 1, nmodes
            DO irq = 1, nrr_q
              WRITE(epwdata,*) czero
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      !
      IF (etf_mem == 0) THEN
        ! SP: The call to epmatwp is now inside the loop
        !     This is important as otherwise the lrepmatw INTEGER
        !     could become too large for integer(kind=4).
        !     Note that in Fortran the record length has to be a integer
        !     of kind 4.
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
        filint   = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwp'
        INQUIRE(IOLENGTH = direct_io_factor) dummy
        unf_recl = direct_io_factor * INT(lrepmatw, KIND = KIND(unf_recl))
        IF (unf_recl <= 0) CALL errore('epw_write', 'wrong record length', 3)
        OPEN(iunepmatwp, FILE = TRIM(ADJUSTL(filint)), IOSTAT = ierr, form='unformatted', &
             STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
        IF (ierr /= 0) CALL errore('epw_write', 'error opening ' // TRIM(filint), 1)
        !
        !CALL diropn(iunepmatwp, 'epmatwp', lrepmatw, exst)
        DO irg = 1, nrr_g
          CALL davcio(epmatwp(:, :, :, :, irg), lrepmatw, iunepmatwp, irg, +1)
        ENDDO
        !
        CLOSE(iunepmatwp)
      ENDIF
      !
      ! Sternheimer and hopping correction matrices are written inside the loop
      ! over q-points in ephwann_shuffle
      !
      ! Write Debye-Waller and momentum matrices
      !
      IF (lwfpt) THEN
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * 3
        CALL diropn(iupmwe, 'cpmew', lrepmatw, exst)
        CALL davcio(cpmew, lrepmatw, iupmwe, 1, +1)
        CLOSE(iupmwe)
        !
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * 3 * nmodes
        CALL diropn(iudwwe, 'dwmatwe', lrepmatw, exst)
        CALL davcio(dwmatwe, lrepmatw, iudwwe, 1, +1)
        CLOSE(iudwwe)
      ENDIF ! lwfpt
      !
      CLOSE(epwdata)
      CLOSE(crystal)
      CLOSE(iunvmedata)
      CLOSE(iundmedata)
      IF (eig_read) CLOSE(iunksdata)
      !
    ENDIF ! ionode
    !--------------------------------------------------------------------------------
    END SUBROUTINE epw_write
    !--------------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------------
    SUBROUTINE epw_read(nrr_k, nrr_q, nrr_g)
    !--------------------------------------------------------------------------------
    !!
    !! Routine to read the real space quantities for fine grid interpolation
    !!
    USE kinds,     ONLY : DP
    USE input,     ONLY : nbndsub, eig_read, etf_mem, lifc, lwfpt, epw_memdist
    USE pwcom,     ONLY : ef
    USE global_var,ONLY : chw, rdw, epmatwp, cdmew, cvmew, chw_ks, zstar, &
                          epsi, crrw, dwmatwe, cpmew, epmatwp_dist
    USE ions_base, ONLY : nat
    USE modes,     ONLY : nmodes
    USE io_global, ONLY : stdout, ionode_id
    USE io_files,  ONLY : prefix, diropn, tmp_dir
    USE io_var,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, &
                          iunepmatwp, iudwwe, iupmwe
    USE ep_constants,     ONLY : czero, zero, eps10
    USE parallelism,      ONLY : para_bounds
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_MODE_RDONLY, MPI_DOUBLE_PRECISION, &
                                 MPI_STATUS_IGNORE, MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_INFO_NULL, &
                                 MPI_DOUBLE_COMPLEX
#endif
#if defined(__NAG)
    USE f90_unix_io,ONLY : flush
#endif
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : world_comm, inter_pool_comm
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
    !! WS-vector index (electrons, phonons, and electron-phonon)
    INTEGER :: ir_start, ir_stop
    !! locally start and end points of WS-vector index (electrons, phonons, and electron-phonon)
    INTEGER :: ipol
    !! Cartesian direction (polarison direction)
    INTEGER :: lrepmatw
    !! Record length
    INTEGER :: ios
    !! Status of files
    INTEGER :: ierr
    !! Error status
    INTEGER*8 :: unf_recl
    !! Record length
    INTEGER :: direct_io_factor
    !! Type of IOlength
    REAL(KIND = DP) :: dummy
    !! Dummy variable
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND) :: epmatwp_recl
    !! Record length for parallel reading epmatwp
    INTEGER(KIND = MPI_OFFSET_KIND) :: epmatwp_offset
    !! Offset for parallel reading epmatwp 
    INTEGER :: epmatwp_block_dtype
    !! MPI block datatype for parallel reading epmatwp
#endif
    INTEGER :: irn
    !! Mode*WS-vector Index
    INTEGER :: irn_loc
    !! Mode*WS-vector Index in local
    INTEGER :: nirn_loc
    !! Number of Mode*WS-vector Index in the local


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
      OPEN(UNIT = iunvmedata, FILE = 'vmedata.fmt', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore ('epw_read', 'error opening vmedata.fmt', iunvmedata)
      OPEN(UNIT = iundmedata, FILE = 'dmedata.fmt', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore ('epw_read', 'error opening dmedata.fmt', iundmedata)
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
    ALLOCATE(cvmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating cvmew', 1)
    ALLOCATE(crrw(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating crrw', 1)
    ALLOCATE(cdmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read', 'Error allocating cdmew', 1)
    chw(:, :, :) = czero
    chw_ks(:, :, :) = czero
    rdw(:, :, :) = czero
    cvmew(:, :, :, :) = czero
    crrw(:, :, :, :) = czero
    cdmew(:, :, :, :) = czero
    !
    IF (mpime == ionode_id) THEN
      !
      DO ibnd = 1, nbndsub
        DO jbnd = 1, nbndsub
          DO irk = 1, nrr_k
            READ(epwdata,*) chw(ibnd, jbnd, irk)
            IF (eig_read) READ(iunksdata,*) chw_ks(ibnd, jbnd, irk)
            DO ipol = 1,3
              READ(iunvmedata,*) cvmew(ipol, ibnd, jbnd, irk)
              READ(iunvmedata,*) crrw(ipol, ibnd, jbnd, irk)
              READ(iundmedata,*) cdmew(ipol, ibnd, jbnd, irk)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      DO imode = 1, nmodes
        DO jmode = 1, nmodes
          DO irq = 1, nrr_q
            READ(epwdata,*) rdw(imode, jmode, irq)
          ENDDO
        ENDDO
      ENDDO
      !
    ENDIF
    !
    CALL mp_bcast(chw, ionode_id, world_comm)
    !
    IF (eig_read) CALL mp_bcast(chw_ks, ionode_id, world_comm)
    CALL mp_bcast(rdw, ionode_id, world_comm)
    IF (lifc) THEN
      ! HM: If lifc=true was used in the previous calculation,
      !     the rdw should be zero.
      IF (MAXVAL(ABS(rdw(:, :, :))) > eps10) THEN
        CALL errore('epw_read', "Non-zero rdw is found in epwdata.fmt. "//&
                                "lifc may have changed since the last calculation.", 1)
      ENDIF
    ELSE
      ! HM: If lifc=false was used in the previous calculation,
      !     the maximum of rdw should have a non-zero value.
      IF (MAXVAL(ABS(rdw(:, :, :))) < eps10) THEN
        CALL errore('epw_read', "rdw is not found in epwdata.fmt. "//&
                                "lifc may have changed since the last calculation.", 1)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(cvmew, ionode_id, world_comm)
    CALL mp_bcast(crrw, ionode_id, world_comm)
    CALL mp_bcast(cdmew, ionode_id, world_comm)
    !
    IF (lifc) THEN
      CALL read_ifc_epw
    ENDIF
    !
    IF (etf_mem == 0) THEN
      !
#if defined(__MPI)
      CALL MPI_TYPE_CONTIGUOUS((nbndsub**2) * nrr_k, MPI_DOUBLE_COMPLEX, epmatwp_block_dtype, ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error creating epmatwp_block_dtype', 1)
      !
      CALL MPI_TYPE_COMMIT(epmatwp_block_dtype, ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error commiting epmatwp_block_dtype', 1)
#endif
      !
      IF (.NOT. epw_memdist) THEN
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
          filint   = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwp'
          !
          INQUIRE(IOLENGTH = direct_io_factor) dummy
          unf_recl = direct_io_factor * INT(lrepmatw, KIND = KIND(unf_recl))
          IF (unf_recl <= 0) CALL errore('epw_read', 'wrong record length', 3)
          OPEN(iunepmatwp, FILE = TRIM(ADJUSTL(filint)), IOSTAT = ierr, FORM = 'unformatted', &
              STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
          IF (ierr /= 0) CALL errore('epw_read', 'error opening ' // TRIM(filint), 1)
          !
          DO irg = 1, nrr_g
            CALL davcio(epmatwp(:, :, :, :, irg), lrepmatw, iunepmatwp, irg, -1)
          ENDDO
          !
          CLOSE(iunepmatwp)
        ENDIF
        !
#if defined(__MPI)
        CALL MPI_BCAST(epmatwp, nmodes * nrr_g, epmatwp_block_dtype, ionode_id, world_comm, ierr)
        IF (ierr /= 0) CALL errore('epw_read', 'Error broadcasting epmatwp', ierr)
#endif
        !
      ELSE ! epw_memdist
#if defined(__MPI)
        !
        WRITE(stdout,'(/5x,"EP vertexes are distributed into MPI processes"/)')
        FLUSH(stdout)
        !
        CALL para_bounds(ir_start, ir_stop, nrr_g * nmodes)
        nirn_loc = ir_stop - ir_start + 1
        !
        ALLOCATE(epmatwp_dist(nbndsub, nbndsub, nrr_k, nirn_loc), STAT = ierr)
        IF (ierr /= 0) CALL errore('epw_read', 'Error allocating epmatwp_dist', 1)
        epmatwp_dist = czero
        !
        filint = TRIM(tmp_dir) // TRIM(prefix) // '.epmatwp'
        CALL MPI_FILE_OPEN(inter_pool_comm, filint, MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmatwp, ierr)
        IF (ierr /= 0) CALL errore('epw_read', 'error in MPI_FILE_OPEN', 1)
        !
        epmatwp_offset = 2_MPI_OFFSET_KIND * 8_MPI_OFFSET_KIND * &
                                    INT(nbndsub , KIND = MPI_OFFSET_KIND) * &
                                    INT(nbndsub , KIND = MPI_OFFSET_KIND) * &
                                    INT(nrr_k, KIND = MPI_OFFSET_KIND) * &
                                    INT(ir_start - 1, KIND = MPI_OFFSET_KIND)
        !                                    
        CALL MPI_FILE_READ_AT_ALL(iunepmatwp, epmatwp_offset, epmatwp_dist, &
          nirn_loc, epmatwp_block_dtype, MPI_STATUS_IGNORE, ierr)
        IF (ierr /= 0) CALL errore('epw_read', 'Error reading epmatwp', 1)
        !
        CALL MPI_FILE_CLOSE(iunepmatwp, ierr)
        IF (ierr /= 0) CALL errore('epw_read', 'Error closing epmatwp', 1)
#else
        CALL errore('epw_read', 'Distributed storage of epmatwp can only be used with MPI', 1)
#endif
      ENDIF ! epw_memdist
      !
#if defined(__MPI)
      CALL MPI_TYPE_FREE(epmatwp_block_dtype, ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error freeing epmatwp_block_dtype', 1)
#endif
      !
    ENDIF
    !
    ! Sternheimer and hopping correction matrices are read inside the loop over
    ! q-points in ephwann_shuffle
    !
    ! Read Debye-Waller and momentum matrices
    !
    IF (lwfpt) THEN
      ALLOCATE(cpmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw_read', 'Error allocating cdmew', 1)
      cpmew = czero
      !
      ALLOCATE(dwmatwe(nbndsub, nbndsub, nrr_k, 3, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error allocating dwmatwe', 1)
      dwmatwe = czero
      !
      IF (mpime == ionode_id) THEN
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * 3
        CALL diropn(iupmwe, 'cpmew', lrepmatw, exst)
        CALL davcio(cpmew, lrepmatw, iupmwe, 1, -1)
        CLOSE(iupmwe)
        !
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * 3 * nmodes
        CALL diropn(iudwwe, 'dwmatwe', lrepmatw, exst)
        CALL davcio(dwmatwe, lrepmatw, iudwwe, 1, -1)
        CLOSE(iudwwe)
      ENDIF ! ionode
      !
      CALL mp_bcast(cpmew, ionode_id, world_comm)
      CALL mp_bcast(dwmatwe, ionode_id, world_comm)
      !
    ENDIF ! lwfpt
    !
    IF (mpime == ionode_id) THEN
      CLOSE(epwdata)
      CLOSE(iunvmedata)
      CLOSE(iundmedata)
    ENDIF
    !
    WRITE(stdout, '(/5x,"Finished reading Wann rep data from file"/)')
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE epw_read
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE epw_write_ws_data(dims, dims2, nrr_k, irvec_k, ndegen_k, wslen_k, &
      nrr_q, irvec_q, ndegen_q, wslen_q, nrr_g, irvec_g, ndegen_g, wslen_g)
    !------------------------------------------------------------------------------
    !!
    !! Write Wigner-Seitz data for restarts and third-party codes.
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : ionode
    USE io_var,    ONLY : iuwigner
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: dims, dims2
    !! Number of Wannier functions and atoms if use_ws, 1 otherwise
    INTEGER, INTENT(in) :: nrr_k, nrr_q, nrr_g
    !! Number of WS vectors (electrons, phonons, and electron-phonon)
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k), irvec_q(3, nrr_q), irvec_g(3, nrr_g)
    !! WS vectors (electrons, phonons, and electron-phonon)
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! WS degeneracies (electrons; indices: irvec_k, iw at 0, iw2 at irvec_k)
    INTEGER, INTENT(in) :: ndegen_q(nrr_q, dims2, dims2)
    !! WS degeneracies (phonons; indices: irvec_q, na at 0, nb at irvec_q)
    INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, dims2)
    !! WS degeneracies (electron-phonon; indices: iw at 0, irvec_g, nb at irvec_q)
    REAL(KIND = DP), INTENT(in) :: wslen_k(nrr_k), wslen_q(nrr_q), wslen_g(nrr_g)
    !! WS-vector lengths (electrons, phonons, and electron-phonon)
    !
    INTEGER :: ir, iw, na
    !! WS-vector, Wannier-function, and atom indices
    !
    IF (.NOT. ionode) RETURN
    !
    OPEN(UNIT=iuwigner, FILE='wigner.fmt', ACTION='write', STATUS='replace')
    !
    WRITE(iuwigner, 1) nrr_k, nrr_q, nrr_g, dims, dims2
    !
    DO ir = 1, nrr_k
      WRITE(iuwigner, 2) irvec_k(:, ir), wslen_k(ir)
      DO iw = 1, dims
        WRITE(iuwigner, 1) ndegen_k(ir, iw, :)
      ENDDO
    ENDDO
    !
    DO ir = 1, nrr_q
      WRITE(iuwigner, 2) irvec_q(:, ir), wslen_q(ir)
      DO na = 1, dims2
        WRITE(iuwigner, 1) ndegen_q(ir, na, :)
      ENDDO
    ENDDO
    !
    DO ir = 1, nrr_g
      WRITE(iuwigner, 2) irvec_g(:, ir), wslen_g(ir)
      DO iw = 1, dims
        WRITE(iuwigner, 1) ndegen_g(iw, ir, :)
      ENDDO
    ENDDO
    !
    CLOSE(iuwigner)
    !
    1 FORMAT(1000(I0, :, 1X)) ! 1000 instead of * for Fortran < 2008
    2 FORMAT(3I6, ES26.17E3) ! should preserve double precision
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE epw_write_ws_data
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE epw_read_ws_data(dims, dims2, nrr_k, irvec_k, ndegen_k, wslen_k, &
      nrr_q, irvec_q, ndegen_q, wslen_q, nrr_g, irvec_g, ndegen_g, wslen_g)
    !------------------------------------------------------------------------------
    !!
    !! Read Wigner-Seitz data during restarts.
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : ionode, ionode_id
    USE mp,        ONLY : mp_bcast
    USE mp_global, ONLY : world_comm
    USE io_var,    ONLY : iuwigner
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: dims, dims2
    !! Number of Wannier functions and atoms if use_ws, 1 otherwise
    INTEGER, INTENT(out) :: nrr_k, nrr_q, nrr_g
    !! Number of WS vectors (electrons, phonons, and electron-phonon)
    INTEGER, ALLOCATABLE, INTENT(out) :: irvec_k(:, :), irvec_q(:, :), irvec_g(:, :)
    !! WS vectors (electrons, phonons, and electron-phonon)
    INTEGER, ALLOCATABLE, INTENT(out) :: ndegen_k(:, :, :)
    !! WS degeneracies (electrons; indices: irvec_k, iw at 0, iw2 at irvec_k)
    INTEGER, ALLOCATABLE, INTENT(out) :: ndegen_q(:, :, :)
    !! WS degeneracies (phonons; indices: irvec_q, na at 0, nb at irvec_q)
    INTEGER, ALLOCATABLE, INTENT(out) :: ndegen_g(:, :, :)
    !! WS degeneracies (electron-phonon; indices: iw at 0, irvec_g, nb at irvec_q)
    REAL(KIND = DP), ALLOCATABLE, INTENT(out) :: wslen_k(:), wslen_q(:), wslen_g(:)
    !! WS-vector lengths (electrons, phonons, and electron-phonon)
    !
    INTEGER :: ir, iw, na
    !! WS-vector, Wannier-function, and atom indices
    INTEGER :: ierr
    !! Error status
    !
    IF (ionode) THEN
      OPEN(UNIT=iuwigner, FILE='wigner.fmt', ACTION='read', STATUS='old')
      !
      READ(iuwigner, *) nrr_k, nrr_q, nrr_g, dims, dims2
    ENDIF
    !
    CALL mp_bcast(nrr_k, ionode_id, world_comm)
    CALL mp_bcast(nrr_q, ionode_id, world_comm)
    CALL mp_bcast(nrr_g, ionode_id, world_comm)
    CALL mp_bcast(dims, ionode_id, world_comm)
    CALL mp_bcast(dims2, ionode_id, world_comm)
    !
    ALLOCATE(irvec_k(3, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating irvec_k', 1)
    ALLOCATE(irvec_q(3, nrr_q), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating irvec_q', 1)
    ALLOCATE(irvec_g(3, nrr_g), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating irvec_g', 1)
    ALLOCATE(ndegen_k(nrr_k, dims, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating ndegen_k', 1)
    ALLOCATE(ndegen_q(nrr_q, dims2, dims2), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating ndegen_q', 1)
    ALLOCATE(ndegen_g(dims, nrr_g, dims2), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating ndegen_g', 1)
    ALLOCATE(wslen_k(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating wslen_k', 1)
    ALLOCATE(wslen_q(nrr_q), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating wslen_q', 1)
    ALLOCATE(wslen_g(nrr_g), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_read_ws_data', 'Error allocating wslen_g', 1)
    !
    IF (ionode) THEN
      DO ir = 1, nrr_k
        READ(iuwigner, *) irvec_k(:, ir), wslen_k(ir)
        DO iw = 1, dims
          READ(iuwigner, *) ndegen_k(ir, iw, :)
        ENDDO
      ENDDO
      !
      DO ir = 1, nrr_q
        READ(iuwigner, *) irvec_q(:, ir), wslen_q(ir)
        DO na = 1, dims2
          READ(iuwigner, *) ndegen_q(ir, na, :)
        ENDDO
      ENDDO
      !
      DO ir = 1, nrr_g
        READ(iuwigner, *) irvec_g(:, ir), wslen_g(ir)
        DO iw = 1, dims
          READ(iuwigner, *) ndegen_g(iw, ir, :)
        ENDDO
      ENDDO
      !
      CLOSE(iuwigner)
    ENDIF
    !
    CALL mp_bcast(irvec_k, ionode_id, world_comm)
    CALL mp_bcast(irvec_q, ionode_id, world_comm)
    CALL mp_bcast(irvec_g, ionode_id, world_comm)
    CALL mp_bcast(ndegen_k, ionode_id, world_comm)
    CALL mp_bcast(ndegen_q, ionode_id, world_comm)
    CALL mp_bcast(ndegen_g, ionode_id, world_comm)
    CALL mp_bcast(wslen_k, ionode_id, world_comm)
    CALL mp_bcast(wslen_q, ionode_id, world_comm)
    CALL mp_bcast(wslen_g, ionode_id, world_comm)
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE epw_read_ws_data
    !------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------
    SUBROUTINE read_ifc_epw()
    !---------------------------------------------------------------------------------
    !!
    !! Read IFC in real space from the file generated by q2r.
    !! Adapted from PH/matdyn.x by C. Verdi and S. Ponce
    !!
    !
    USE kinds,     ONLY : DP
    USE global_var,ONLY : ifc, zstar, epsi, alph
    USE input,     ONLY : asr_typ, dvscf_dir, nqc1, nqc2, nqc3
    USE ions_base, ONLY : nat
    USE cell_base, ONLY : ibrav, omega, at, bg, celldm, alat
    USE io_global, ONLY : stdout, ionode_id
    USE io_var,    ONLY : iunifc
    USE noncollin_module, ONLY : nspin_mag
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : intra_pool_comm, inter_pool_comm, root_pool
    USE mp_world,  ONLY : mpime, world_comm
    USE io_dyn_mat,ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                          read_ifc_param, read_ifc
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
    IF (is_xml_file) THEN
      !IF (mpime == 0) THEN
      !  ionode = .TRUE.
      !ELSE
      !  ionode = .FALSE.
      !ENDIF
      !
      ! pass the 'tempfile' as the '.xml' extension is added in the next routine
      CALL read_dyn_mat_param(tempfile, ntyp_, nat_)
      ALLOCATE(m_loc(3, nat_), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_ifc_epw', 'Error allocating m_loc', 1)
      ALLOCATE(atm(ntyp_), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_ifc_epw', 'Error allocating atm', 1)
      CALL read_dyn_mat_header(ntyp_, nat_, ibrav, nspin_mag, &
               celldm, at, bg, omega, atm, amass2, &
               tau_, ityp_,  m_loc, nqs, has_zstar, epsi, zstar)
      CALL volume(alat, at(1, 1), at(1, 2), at(1, 3), omega)
      CALL read_ifc_param(nqc1, nqc2, nqc3)
      CALL read_ifc(alph, nqc1, nqc2, nqc3, nat_, ifc)
      DEALLOCATE(m_loc, STAT = ierr)
      IF (ierr /= 0) CALL errore('read_ifc_epw', 'Error deallocating m_loc', 1)
      DEALLOCATE(atm, STAT = ierr)
      IF (ierr /= 0) CALL errore('read_ifc_epw', 'Error deallocating atm', 1)
      !
    ELSE ! is_xml_file
      IF (mpime == ionode_id) THEN
        !
        OPEN(UNIT = iunifc, FILE = tempfile, STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) call errore ('read_ifc_epw', 'error opening ifc.q2r', iunifc)
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
        CLOSE(iunifc)
      ENDIF ! mpime
      ! It has to be casted like this because mpi cannot cast 7 indices
      DO i = 1, 3
        DO j = 1, 3
          DO na = 1, nat
            DO nb = 1, nat
              CALL mp_bcast(ifc(:, :, :, i, j, na, nb), ionode_id, world_comm)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL mp_bcast(zstar, ionode_id, world_comm)
      CALL mp_bcast(epsi, ionode_id, world_comm)
      CALL mp_bcast(tau_, ionode_id, world_comm)
      CALL mp_bcast(ibrav_, ionode_id, world_comm)
    ENDIF ! has_xml
    !
    WRITE(stdout,'(5x,"IFC last ", 1f12.7)') ifc(nqc1, nqc2, nqc3, 3, 3, nat, nat)
    !
    CALL set_asr2(asr_typ, nqc1, nqc2, nqc3, ifc, zstar, nat, ibrav_, tau_)
    !
    WRITE(stdout, '(/5x,"Finished reading ifcs"/)')
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE read_ifc_epw
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
    USE input,     ONLY : dvscf_dir
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
    CHARACTER(LEN = 4) :: filelab
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
    !-------------------------------------------------------------
    SUBROUTINE readint3paw(int3paw, recn, iq, nqc)
    !-------------------------------------------------------------
    !!
    !! Open int3paw files as direct access, read, and close again
    !!
    !! HL - Mar 2020 based on the subroutine of readdvscf
    !!
    !-------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_files,         ONLY : prefix
    USE units_ph,         ONLY : lint3paw
    USE input,            ONLY : dvscf_dir
    USE io_var,           ONLY : iuint3paw
    USE low_lvl,          ONLY : set_ndnmbr
    USE uspp_param,       ONLY : nhm
    USE ions_base,        ONLY : nat
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
    COMPLEX(KIND = DP), INTENT(inout) :: int3paw(nhm, nhm, nat, nspin_mag)
    !! int3paw is read from file
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: tempfile
    !! Temp file
    CHARACTER(LEN = 4) :: filelab
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
    tempfile = TRIM(dvscf_dir) // TRIM(prefix) // '.dvscf_paw_q' // filelab
    INQUIRE(IOLENGTH = unf_recl) dummy
    unf_recl = unf_recl  * lint3paw
    mult_unit = unf_recl
    mult_unit = recn * mult_unit
    !
    !  open the dvscf_paw file, read and close
    !
    OPEN(iuint3paw, FILE = tempfile, FORM = 'unformatted', &
         ACCESS = 'direct', IOSTAT = ios, RECL = unf_recl, STATUS = 'old')
    IF (ios /= 0) CALL errore('readint3paw', 'error opening ' // tempfile, iuint3paw)
    !
    ! check that the binary file is long enough
    INQUIRE(FILE = tempfile, SIZE = file_size)
    IF (mult_unit > file_size) CALL errore('readint3paw', &
         TRIM(tempfile) //' too short', iuint3paw)
    !
    READ(iuint3paw, REC = recn) int3paw
    CLOSE(iuint3paw, STATUS = 'keep')
    !
    RETURN
    !
    !-------------------------------------------------------------
    END SUBROUTINE readint3paw
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
    CHARACTER(LEN = 4) :: nd_nmbr0
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
    SUBROUTINE readgmap(nkstot)
    !--------------------------------------------------------------
    !!
    !! read the map of G vectors G -> G-G_0 for a given q point
    !! (this is used for the folding of k+q into the first BZ)
    !!
    !
    USE kinds,      ONLY : DP
    USE mp_global,  ONLY : world_comm
    USE mp,         ONLY : mp_bcast, mp_max
    USE io_global,  ONLY : meta_ionode, meta_ionode_id
    USE io_var,     ONLY : iukgmap
    USE global_var, ONLY : ngxxf, ngxx, ng0vec, shift, gmap, g0vec_all_r
    USE io_files,   ONLY : prefix
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nkstot
    !! Total number of k-points
    !
    ! Lork variables
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ik1
    !! Temporary indices when reading kgmap files
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
    !
    IF (meta_ionode) THEN
      !
      OPEN(iukgmap, FILE = TRIM(prefix)//'.kgmap', FORM = 'formatted', STATUS = 'old', IOSTAT = ios)
      IF (ios /=0) CALL errore('readgmap', 'error opening kgmap file', iukgmap)
      !
      READ(iukgmap, *) ngxxf
      !
      !! HL: The part below should be removed later since it is useless.
      !
      DO ik = 1, nkstot
        READ(iukgmap, *) ik1, shift(ik1)
      ENDDO
      !
      READ(iukgmap, *) ng0vec
      !
    ENDIF
    !
    ! first node broadcasts ng0vec to all nodes for allocation of gmap
    !
    CALL mp_bcast(ngxxf, meta_ionode_id, world_comm)
    CALL mp_bcast(ng0vec, meta_ionode_id, world_comm)
    !
    ALLOCATE(gmap(ngxx * ng0vec), STAT = ierr)
    IF (ierr /= 0) CALL errore('readgmap', 'Error allocating gmap', 1)
    gmap(:) = 0
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
        READ(iukgmap,*) (gmap(ng0vec * (ig - 1) + ishift), ishift = 1, ng0vec)
      ENDDO
      !
      CLOSE(iukgmap)
      !
    ENDIF
    !
    ! first node broadcasts arrays to all nodes
    !
    CALL mp_bcast(g0vec_all_r, meta_ionode_id, world_comm)
    CALL mp_bcast(gmap, meta_ionode_id, world_comm)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE readgmap
    !-----------------------------------------------------------------------
    !
    !--------------------------------------------------------------
    SUBROUTINE readkmap(nkstot)
    !--------------------------------------------------------------
    !!
    !! Read the index of G_0 such that k+q+G_0 belongs to the 1st BZ
    !!
    !
    USE kinds,      ONLY : DP
    USE mp_global,  ONLY : world_comm
    USE mp,         ONLY : mp_bcast
    USE io_global,  ONLY : meta_ionode, meta_ionode_id
    USE io_var,     ONLY : iukmap
    USE io_files,   ONLY : prefix
    USE global_var, ONLY : shift
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nkstot
    !! Total number of k-points
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ik1
    !! Temporary indices when reading kmap files
    INTEGER :: itmp
    !! Temporary indices when reading kmap files
    INTEGER :: ios
    !! Integer variable for I/O control
    !
    IF (meta_ionode) THEN
      !
      OPEN(iukmap, FILE = TRIM(prefix)//'.kmap', FORM = 'formatted', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore('readkmap', 'error opening kmap file', iukmap)
      DO ik = 1, nkstot
        READ(iukmap,*) ik1, itmp, shift(ik1)
      ENDDO
      CLOSE(iukmap)
      !
    ENDIF
    !
    ! first node broadcasts shift to all nodes
    !
    CALL mp_bcast(shift, meta_ionode_id, world_comm)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE readkmap
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
    USE units_ph,         ONLY : lrdrho, lint3paw
    USE fft_base,         ONLY : dfftp
    USE uspp_param,       ONLY : nhm
    USE ions_base,        ONLY : nat
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
    lint3paw = 2 * nhm * nhm * nat * nspin_mag
    !
    RETURN
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE openfilepw
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE param_get_range_vector(field, length, lcount, i_value)
    !----------------------------------------------------------------------------
    !!
    !!   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100
    !!   if(lcount) we return the number of states in length
    !!
    !!   HL - April 2020
    !!   Imported and adapted from the same name of subroutine in parameters.F90
    !!   in the directory of src in Wannier90
    !!
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(in)         :: field
    !! Field read for parsing
    INTEGER, INTENT(inout)           :: length
    !! Number of states
    LOGICAL, INTENT(in)              :: lcount
    !! If T only count states
    INTEGER, OPTIONAL, INTENT(out)   :: i_value(length)
    !! States specified in range vector
    !
    INTEGER :: loop
    !! Loop index
    INTEGER :: num1
    !! Integer number read
    INTEGER :: num2
    !! Integer number read
    INTEGER :: i_punc
    !! Position returned after scanning punctuation marks
    INTEGER :: counter
    !! Counter index
    INTEGER :: i_digit
    !! Position returned after scanning numbers
    INTEGER :: loop_r
    !! Loop index
    INTEGER :: range_size
    !! Size of range
    CHARACTER(LEN = 255) :: dummy
    !! Copy of field read for parsing
    CHARACTER(LEN = 10), PARAMETER :: c_digit = "0123456789"
    CHARACTER(LEN = 2), PARAMETER :: c_range = "-:"
    CHARACTER(LEN = 3), PARAMETER :: c_sep = " ,;"
    CHARACTER(LEN = 5), PARAMETER :: c_punc = " ,;-:"
    CHARACTER(LEN = 5) :: c_num1
    !! Number read
    CHARACTER(LEN = 5) :: c_num2
    !! Number read
    !
    IF (lcount .AND. PRESENT(i_value)) &
      CALL errore('param_get_range_vector', 'incorrect call', 1)
    !
    dummy = field
    dummy = ADJUSTL(dummy)
    !
    counter = 0
    IF (LEN_TRIM(dummy) == 0) THEN
      length = counter
      RETURN
    ENDIF
    !
    DO
      i_punc = SCAN(dummy, c_punc)
      IF (i_punc == 0) &
        CALL errore('param_get_range_vector', 'Error parsing field', 1)
      c_num1 = dummy(1:i_punc - 1)
      READ(c_num1, *, ERR = 101, END = 101) num1
      dummy = ADJUSTL(dummy(i_punc:))
      !look for range
      IF (SCAN(dummy, c_range) == 1) THEN
        i_digit = SCAN(dummy, c_digit)
        dummy = ADJUSTL(dummy(i_digit:))
        i_punc = SCAN(dummy, c_punc)
        c_num2 = dummy(1:i_punc - 1)
        READ(c_num2, *, ERR = 101, END = 101) num2
        dummy = ADJUSTL(dummy(i_punc:))
        range_size = ABS(num2 - num1) + 1
        DO loop_r = 1, range_size
          counter = counter + 1
          IF (.NOT. lcount) i_value(counter) = MIN(num1, num2) + loop_r - 1
        ENDDO
      ELSE
        counter = counter + 1
        IF (.NOT. lcount) i_value(counter) = num1
      ENDIF
      IF (SCAN(dummy, c_sep) == 1) dummy = ADJUSTL(dummy(2:))
      IF (SCAN(dummy, c_range) == 1) &
        CALL errore('param_get_range_vector', 'Error parsing field: incorrect range', 1)
      IF (INDEX(dummy, ' ') == 1) EXIT
    ENDDO
    !
    IF (lcount) length = counter
    IF (.NOT. lcount) THEN
      DO loop = 1, counter - 1
        DO loop_r = loop + 1, counter
          IF (i_value(loop) == i_value(loop_r)) &
            CALL errore('param_get_range_vector', 'Error parsing field: duplicate values', 1)
        ENDDO
      ENDDO
    ENDIF
    !
    RETURN
    !
101 CALL errore('param_get_range_vector', 'Error parsing field', 1)
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE param_get_range_vector
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE dynmat_asr(iq_irr, nqc_irr, nq, iq_first, sxq, imq, isq, &
                                invs, s, irt, rtau, sumr)
    !-----------------------------------------------------------------------
    !!
    !! read dynamical matrix for the q points, either in plain text or xml.
    !! iq_first, iq_first+1, ... iq_first+nq-1
    !!
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    use io_files,         ONLY : prefix
    USE cell_base,        ONLY : ibrav, celldm, omega, at, bg
    USE ions_base,        ONLY : amass, tau, nat, ntyp => nsp, ityp
    USE global_var,       ONLY : dynq, zstar, epsi, wf_temp
    USE symm_base,        ONLY : nsym
    USE input,            ONLY : dvscf_dir, lpolar, nqc1, nqc2, nqc3, exciton
    USE modes,            ONLY : nmodes
    USE control_flags,    ONLY : iverbosity
    USE noncollin_module, ONLY : nspin_mag
    USE io_global,        ONLY : ionode, stdout, ionode_id, meta_ionode_id
    USE ep_constants,     ONLY : cone, czero, twopi, rydcm1, eps6, zero
    USE io_var,           ONLY : iudyn
    USE low_lvl,          ONLY : set_ndnmbr, eqvect_strict
    USE io_dyn_mat,       ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                                 read_dyn_mat, read_dyn_mat_tail
    USE mp_world,         ONLY : mpime
    USE mp_global,        ONLY : world_comm
    USE mp,               ONLY : mp_bcast
    USE rigid,            ONLY : cdiagh2
    !
    IMPLICIT NONE
    !
    ! Input
    INTEGER, INTENT(in) :: iq_irr
    !! The index of the irreducible q point
    INTEGER, INTENT(in) :: nqc_irr
    !! The number of irreducible qpoints
    INTEGER, INTENT(in) :: nq
    !! The number of q points in the star of q
    INTEGER, INTENT(in) :: iq_first
    !! The index of the first qpoint to be read in the uniform q-grid
    INTEGER, INTENT(in) :: imq
    !! Flag which tells whether we have to consider the -q vectors
    INTEGER, INTENT(in) :: isq(48)
    !! Index of q in the star for a given symmetry.
    INTEGER, INTENT(in) :: invs(48)
    !! Matrix of inversed symmetry
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! Set of symmetry operations
    INTEGER, INTENT(in) :: irt(48, nat)
    !! For each atoms give the rotated atoms
    REAL(KIND = DP), INTENT(in) :: sxq(3, 48)
    !! Symmetry matrix
    REAL(KIND = DP), INTENT(inout) :: rtau(3, 48, nat)
    !! the relative position of the rotated atom to the original one
    REAL(KIND = DP), INTENT(inout) :: sumr(2, 3, nat, 3)
    !! Sum to impose the ASR
    !
    ! Local variables
    LOGICAL :: found
    !! Is it found
    LOGICAL :: lrigid
    !! ifc
    LOGICAL :: lraman
    !! is raman present
    LOGICAL :: nog
    !!
    LOGICAL :: is_xml_file
    !! Is the file XML
    CHARACTER(LEN = 3), ALLOCATABLE :: atm(:)
    !! dummy variable for atom types when xml is used
    CHARACTER(LEN = 3) :: atm_
    !! dummy variable for atom types when xml is not used
    CHARACTER(LEN = 4) :: filelab
    !!
    CHARACTER(LEN = 80) :: line
    !!
    CHARACTER(LEN = 256) :: tempfile
    !!
    INTEGER :: info
    !! "0" successful exit, "<0" i-th argument had an illegal value, ">0" i eigenvectors failed to converge.
    INTEGER :: isym
    !! Symmetry number
    INTEGER :: m1, m2, m3
    !! Supercell index
    INTEGER :: sna
    !! Rotation + translation
    INTEGER :: jsym
    !! Symmetry index
    INTEGER :: nsq
    !! nsq is the degeneracy of the small group for this iq in the star
    INTEGER :: sym_sgq(48)
    !! Degenerate sym
    INTEGER :: current_iq
    !! the q index in the dynq matrix
    INTEGER :: mq
    !! Size of q
    INTEGER :: ism1
    !! Inversion symmetry
    INTEGER :: imode
    !! i-mode
    INTEGER :: jmode
    !! j-mode
    INTEGER :: ntyp_
    !!
    INTEGER :: nat_
    !!
    INTEGER :: ibrav_
    !!
    INTEGER :: ityp_
    !!
    INTEGER :: ios
    !!
    INTEGER :: iq, jq
    !!
    INTEGER :: nt
    !! Number of type of atoms
    INTEGER :: na, nb
    !!
    INTEGER :: naa, nbb
    !!
    INTEGER :: nu, mu
    !!
    INTEGER :: i, j
    !!
    INTEGER :: ipol, jpol
    !!
    INTEGER :: nqs
    !!
    !INTEGER :: nrws
    !!
    INTEGER :: ierr
    !! Error status
    INTEGER, parameter :: ntypx = 10
    !!
    INTEGER, PARAMETER :: nrwsx = 200
    !!
    REAL(KIND = DP) :: arg
    !!
    REAL(KIND = DP) :: aq(3)
    !!
    REAL(KIND = DP) :: saq(3)
    !!
    REAL(KIND = DP) :: raq(3)
    !!
    REAL(KIND = DP) :: xq(3)
    !! The rotated q-vector
    REAL(KIND = DP) :: scart(3, 3)
    !!
    REAL(KIND = DP) :: massfac
    !!
    REAL(KIND = DP) :: w1(nmodes)
    !!
    REAL(KIND = DP) :: wtmp(nmodes)
    !!
    REAL(KIND = DP) :: celldm_(6)
    !!
    REAL(KIND = DP) :: amass_
    !!
    REAL(KIND = DP) :: tau_(3)
    !!
    REAL(KIND = DP) :: q(3, nqc1 * nqc2 * nqc3)
    !!
    REAL(KIND = DP) :: dynr(2, 3, nat, 3, nat)
    !!
    REAL(KIND = DP) :: sumz
    !!
    REAL(KIND = DP) :: qout(3)
    !!
    REAL(KIND = DP) :: amass2(ntypx)
    !!
    !REAL(KIND = DP) :: rws(0:3, nrwsx)
    !!
    !REAL(KIND = DP) :: atws(3, 3)
    !! WS cell dist.
    REAL(KIND = DP) :: zstar_(3, 3, nat)
    !!  Temporary array to store Z*
    REAL(KIND = DP) :: epsi_(3, 3)
    !! Temporary array to store dielectric function
    REAL(KIND = DP), ALLOCATABLE :: m_loc(:, :)
    !!
    REAL(KIND = DP), ALLOCATABLE :: dchi_dtau(:, :, :, :)
    !!
    COMPLEX(KIND = DP) :: cfac
    !!
    COMPLEX(KIND = DP) :: gamma(nmodes, nmodes)
    !! the Gamma matrix for the symmetry operation on the dyn mat
    COMPLEX(KIND = DP) :: dynq_tmp(nmodes, nmodes)
    !!
    COMPLEX(KIND = DP) :: cz1(nmodes, nmodes)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: dyn(:, :, :, :) ! 3,3,nat,nat
    !! Dynamical matrix
    !
    q(:, :) = zero
    !
    ! the call to set_ndnmbr is just a trick to get quickly
    ! a file label by exploiting an existing subroutine
    ! (if you look at the sub you will find that the original
    ! purpose was for pools and nodes)
    !
    CALL set_ndnmbr(0, iq_irr, 1, nqc_irr, filelab)
    tempfile = TRIM(dvscf_dir) // TRIM(prefix) // '.dyn_q' // TRIM(filelab)
    ! The following function will check either or not the file is formatted in
    ! xml. If no file is found, an error is raised
    CALL check_is_xml_file(tempfile, is_xml_file)
    !
    IF (is_xml_file) THEN
      CALL read_dyn_mat_param(tempfile, ntyp, nat)
      ALLOCATE(m_loc(3, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynmat', 'Error allocating m_loc', 1)
      ALLOCATE(atm(ntyp), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynmat', 'Error allocating atm', 1)
      ALLOCATE(dchi_dtau(3, 3, 3, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynmat', 'Error allocating dchi_dtau', 1)
      CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
              celldm, at, bg, omega, atm, amass2, tau, ityp, &
              m_loc, nqs, lrigid, epsi_, zstar_, lraman, dchi_dtau)
      DEALLOCATE(m_loc, STAT = ierr)
      IF (ierr /= 0) CALL errore('dynmat', 'Error deallocating m_loc', 1)
      DEALLOCATE(atm, STAT = ierr)
      IF (ierr /= 0) CALL errore('dynmat', 'Error deallocating atm', 1)
      DEALLOCATE(dchi_dtau, STAT = ierr)
      IF (ierr /= 0) CALL errore('dynmat', 'Error deallocating dchi_dtau', 1)
      !
      ALLOCATE(dyn(3, 3, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('dynmat', 'Error allocating dyn', 1)
      IF (ionode) THEN
        DO nt = 1, ntyp
          IF (amass(nt) <= 0.0d0) amass(nt) = amass2(nt)
        ENDDO
      ENDIF
      !
      IF (lrigid) THEN
        WRITE(stdout, '(5x,a)') 'Read dielectric tensor and effective charges'
        zstar = zstar_
        epsi = epsi_
        !ASR on effective charges
        DO i = 1, 3
          DO j = 1, 3
            sumz = 0.0d0
            DO na = 1,nat
              sumz = sumz + zstar(i, j, na)
            ENDDO
            DO na = 1,nat
              zstar(i, j, na) = zstar(i, j, na) - sumz / nat
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      !
      ! If time-reversal is not included in the star of q, then double the nq to
      ! search from.
      IF (imq == 0) THEN
        mq = 2 * nq
      ELSE
        mq = nq
      ENDIF
      !
      ! First read the data and then sort them. This is because the order of
      ! q in the star has changed (SP).
      !
      DO iq = 1, mq
        q(:, iq) = 0.0d0
        CALL read_dyn_mat(nat, iq, qout, dyn(:, :, :, :))
        q(:, iq) = qout(:)
        !
        DO na = 1,nat
          DO ipol = 1, 3
            DO jpol = 1, 3
              dynr(1, ipol, na, jpol, :) = REAL(dyn(ipol, jpol, na, :))
              dynr(2, ipol, na, jpol, :) = REAL(AIMAG(dyn(ipol, jpol, na, :)))
            ENDDO
          ENDDO
        ENDDO
        !
        ! Impose the acoustic sum rule (q=0 needs to be the first q point in the coarse grid)
        ! [Gonze and Lee, PRB 55, 10361 (1998), Eq. (45) and (81)]
        !
        IF (ABS(q(1, iq)) < eps6 .AND. ABS(q(2, iq)) < eps6 .AND. ABS(q(3, iq)) < eps6) THEN
          WRITE(stdout, '(5x,a)') 'Imposing acoustic sum rule on the dynamical matrix'
          IF (lpolar .AND. .NOT. lrigid) CALL errore('dynmat', &
            &'You set lpolar = .TRUE. but did not put epsil = true in the PH calculation at Gamma. ',1)
        ENDIF
        DO na = 1, nat
          DO ipol = 1, 3
            DO jpol = ipol, 3
              !
              IF (ABS(q(1, iq)) < eps6 .AND. ABS(q(2, iq)) < eps6 .AND. ABS(q(3, iq)) < eps6 ) THEN
                sumr(1, ipol, na, jpol) = SUM(dynr(1, ipol, na, jpol, :))
                sumr(2, ipol, na, jpol) = SUM(dynr(2, ipol, na, jpol, :))
              ENDIF
              !
              dynr(:, ipol, na, jpol, na) = dynr(:, ipol, na, jpol, na) - sumr(:, ipol, na, jpol)
              !
            ENDDO
          ENDDO
        ENDDO
        !
        ! Fill the two-indices dynamical matrix in cartesian coordinates
        ! the proper index in the complete list is iq_first+iq-1
        !
        DO na = 1, nat
          DO nb = 1, nat
            DO ipol = 1, 3
              DO jpol = 1, 3
                !
                mu = (na - 1) * 3 + ipol
                nu = (nb - 1) * 3 + jpol
                ! Only store the dyn of the first q in the star.
                IF (iq == 1) THEN
                  dynq(mu, nu, iq_first) = DCMPLX(dynr(1, ipol, na, jpol, nb), dynr(2, ipol, na, jpol, nb))
                ENDIF
                !
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !
      ENDDO !  iq = 1, mq
      !
      ! Close the dyn file
      CALL read_dyn_mat_tail(nat)
      !
    ELSE ! not a xml file
      IF (mpime == ionode_id) THEN
        OPEN(UNIT = iudyn, FILE = tempfile, STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) CALL errore('dynmat_asr', 'opening file' // tempfile, ABS(ios))
        !
        !  read header and run some checks
        !
        READ(iudyn, '(a)') line
        READ(iudyn, '(a)') line
        READ(iudyn, *) ntyp_, nat_, ibrav_, celldm_
        !
        ! We stop testing celldm as it can be different between scf and nscf
        !IF (ntyp/=ntyp_.OR.nat/=nat_.OR.ibrav_/=ibrav.OR.abs ( &
        !   celldm_ (1) - celldm (1) )  > 1.0d-5) call errore ('dynmat', &
        !   'inconsistent data', 1)
        IF (ntyp /= ntyp_ .OR. nat /= nat_ .OR. ibrav_ /= ibrav) CALL errore ('dynmat_asr', &
           'inconsistent data', 1)
        !
        !  skip reading of cell parameters here
        !
        IF (ibrav_ == 0) THEN
          DO i = 1, 4
            READ(iudyn, *) line
          ENDDO
        ENDIF
        DO nt = 1, ntyp
          READ(iudyn, *) i, atm_, amass_
          IF (nt /=i .OR. ABS(amass_ - amass(nt))  > 1.0d-2) THEN
            WRITE(stdout,*) amass_, amass(nt)
            CALL errore('dynmat_asr', 'inconsistent data', 1)
          ENDIF
        ENDDO
        DO na = 1, nat
          READ(iudyn, * ) i, ityp_, tau_
          IF (na /= i .OR. ityp_ /= ityp(na)) CALL errore('dynmat_asr', &
              'inconsistent data (names)', 10 + na)
        ENDDO
        !
        ! Read dyn mat only for the first q in the star and then reconstruct the
        ! other using symmetry.
        !
        ! If time-reversal is not included in the star of q, then double the nq to
        ! search from.
        IF (imq == 0) THEN
          mq = 2 * nq
        ELSE
          mq = nq
        ENDIF
        !
        ! Read the dyn for the first q in the star.
        ! For other q in the star, only read the value of the q for checking purposes.
        !
        DO iq = 1, mq
          !
          READ(iudyn, '(///a)') line
          READ(line(11:80), *) (q(i, iq), i = 1, 3)
          READ(iudyn, '(a)') line
          !
          DO na = 1, nat
            DO nb = 1, nat
              READ(iudyn, *) naa, nbb
              IF (na /= naa .OR. nb /= nbb) CALL errore('dynmat_asr', 'error reading file', nb)
              READ(iudyn, *)((dynr(1, i, na, j, nb), dynr(2, i, na, j, nb), j = 1, 3), i = 1, 3)
            ENDDO
          ENDDO
          !
          ! impose the acoustic sum rule (q=0 needs to be the first q point in the coarse grid)
          ! [Gonze and Lee, PRB 55, 10361 (1998), Eq. (45) and (81)]
          !
          IF (ABS(q(1, iq)) < eps6 .AND. ABS(q(2, iq)) < eps6 .AND. ABS(q(3, iq)) < eps6) THEN
            WRITE(stdout, '(5x,a)') 'Imposing acoustic sum rule on the dynamical matrix'
          ENDIF
          DO na = 1, nat
            DO ipol = 1,3
              DO jpol = ipol, 3
                !
                IF (ABS(q(1, iq)) < eps6 .AND. ABS(q(2, iq)) < eps6 .AND. ABS(q(3, iq)) < eps6) THEN
                  sumr(1, ipol, na, jpol) = SUM(dynr(1, ipol, na, jpol, :))
                  sumr(2, ipol, na, jpol) = SUM(dynr(2, ipol, na, jpol, :))
                ENDIF
                !
                dynr(:, ipol, na, jpol, na) = dynr(:, ipol, na, jpol, na) - sumr(:, ipol, na, jpol)
                !
              ENDDO
            ENDDO
          ENDDO
          !
          !  fill the two-indices dynamical matrix in cartesian coordinates
          !  the proper index in the complete list is iq_first+iq-1
          !
          DO na = 1, nat
            DO nb = 1, nat
              DO ipol = 1, 3
                DO jpol = 1, 3
                  !
                  mu = (na - 1) * 3 + ipol
                  nu = (nb - 1) * 3 + jpol
                  ! Only store the dyn of the first q in the star.
                  IF (iq == 1) THEN
                    dynq(mu, nu, iq_first) = DCMPLX(dynr(1, ipol, na, jpol, nb), dynr(2, ipol, na, jpol, nb))
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !
        IF (ABS(q(1, 1)) < eps6 .AND. ABS(q(2, 1)) < eps6 .AND. ABS(q(3, 1)) < eps6) THEN
          !  read dielectric tensor and effective charges if present
          !  SP: Warning zstar is not properly bcast at the moment
          READ(iudyn, '(a)') line
          READ(iudyn, '(a)') line
          !
          IF (line(6:15) .EQ. 'Dielectric') THEN
            READ(iudyn, '(a)') line
            READ(iudyn, *) ((epsi(i, j), j = 1, 3), i = 1, 3)
            READ(iudyn, '(a)') line
            READ(iudyn, '(a)') line
            READ(iudyn, '(a)') line
            DO na = 1, nat
              READ(iudyn, '(a)') line
              READ(iudyn, *) ((zstar(i, j, na), j = 1, 3), i = 1, 3)
            ENDDO
            WRITE(stdout,'(5x,a)') 'Read dielectric tensor and effective charges'
            !
            !ASR on effective charges
            DO i = 1, 3
              DO j = 1, 3
                sumz = 0.0d0
                DO na = 1, nat
                  sumz = sumz + zstar(i, j, na)
                ENDDO
                DO na = 1, nat
                  zstar(i, j, na) = zstar(i, j, na) - sumz / nat
                ENDDO
              ENDDO
            ENDDO
            !
          ELSE
            IF (lpolar) CALL errore('dynmat_asr', &
               'You set lpolar = .TRUE. but did not put epsil = true in the PH calculation at Gamma. ', 1)
          ENDIF
        ENDIF
        CLOSE(iudyn)
      ENDIF ! mpime
      CALL mp_bcast(zstar, meta_ionode_id, world_comm)
      CALL mp_bcast(epsi , meta_ionode_id, world_comm)
      CALL mp_bcast(dynq , meta_ionode_id, world_comm)
      CALL mp_bcast(q    , meta_ionode_id, world_comm)
      CALL mp_bcast(mq   , meta_ionode_id, world_comm)
    ENDIF ! not xml
    !
    !
    ! Now check that the dyn file is consistent with the current EPW run (SP)
    ! SP: Be careful, here time-reversal is not actual time reversal but is due to
    !     change in order and values of the q in the star between QE 4 and 5.
    !
    CALL cryst_to_cart(nqc1 * nqc2 * nqc3, q, at, -1)
    CALL cryst_to_cart(48, sxq, at, -1)
    !
    current_iq = iq_first
    !
    DO iq = 1, nq
      !
      found = .FALSE.
      DO jq = 1, mq
        DO m1 = -2, 2
          DO m2 = -2, 2
            DO m3 = -2, 2
              IF ((ABS(q(1, jq) - (sxq(1, iq) + m1)) < eps6 .AND. &
                   ABS(q(2, jq) - (sxq(2, iq) + m2)) < eps6 .AND. &
                   ABS(q(3, jq) - (sxq(3, iq) + m3)) < eps6 )) THEN
                found = .TRUE.
                EXIT ! exit loop
              ENDIF
            ENDDO
            IF (found) EXIT
          ENDDO
          IF (found) EXIT
        ENDDO
        IF (found) EXIT
      ENDDO
      current_iq = current_iq + 1
      IF (found .EQV. .FALSE.) THEN
        CALL errore('dynmat_asr', 'wrong qpoint', 1)
      ENDIF
    ENDDO
    ! Transform back the sxq in Cartesian
    CALL cryst_to_cart(48, sxq, bg, 1)
    !
    ! In case of reading of the IFC to impose the asr in real space
    ! We still call the above just to make the checks. The content of dynq
    ! will be re-written just below and NOT read from the dyn from the /save folder
    !IF (lifc) THEN
    !  !
    !  ! build the WS cell corresponding to the force constant grid
    !  atws(:, 1) = at(:, 1) * DBLE(nqc1)
    !  atws(:, 2) = at(:, 2) * DBLE(nqc2)
    !  atws(:, 3) = at(:, 3) * DBLE(nqc3)
    !  ! initialize WS r-vectors
    !  CALL wsinit(rws, nrwsx, nrws, atws)
    !  ! dynifc2blochc requires ifc
    !  CALL dynifc2blochc(nmodes, rws, nrws, q(:, 1), dynq_tmp)
    !  dynq(:, :, iq_first) = dynq_tmp
    !  WRITE(stdout, '(5x,a)') "Dyn mat calculated from ifcs"
    !  !
    !ENDIF
    !
    ! Now construct the other dyn matrix for the q in the star using sym.
    ! For this we use the gamma matrix.
    !
    current_iq = iq_first
    !
    DO iq = 1, nq
      !
      xq = sxq(:, iq)
      nsq = 0 ! nsq is the degeneracy of the small group for this iq in the star
      sym_sgq(:) = 0
      DO jsym = 1, nsym
        IF (isq(jsym) == iq ) then
          nsq = nsq + 1
          sym_sgq(nsq) = jsym
        ENDIF
      ENDDO
      !
      ! SP: We now need to select one symmetry among the small group of q (i.e. that respect
      !     Sq0+G=q ) that has G=0. There should always be such symmetry.
      !     We enforce this for later easiness.
      !
      aq = sxq(:, 1) ! This is xq0
      saq = xq
      call cryst_to_cart(1, aq, at, - 1)
      CALL cryst_to_cart(1, saq, at, -1)
      ! Initialize isym
      isym = 1
      DO jsym = 1, nsq
        ism1 = invs(sym_sgq(jsym))
        raq = 0.d0
        DO ipol = 1, 3
          DO jpol = 1, 3
            raq(ipol) = raq(ipol) + s(ipol, jpol, ism1) * aq(jpol)
          ENDDO
        ENDDO
        nog = eqvect_strict(raq, saq, eps6)
        IF (nog) THEN ! This is the symmetry such that Sq=q
          isym = sym_sgq(jsym)
          EXIT
        ENDIF
        ! If we enter into that loop it means that we have not found
        ! such symmetry within the small group of Q.
        IF (jsym == nsq) THEN
          CALL errore('dynmat ', 'No sym. such that Sxq0=iq was found in the sgq !', 1)
        ENDIF
      ENDDO
      !
      !  -----------------------------------------------------------------------
      !  the matrix gamma (Maradudin & Vosko, RMP, eq. 2.37)
      !  -----------------------------------------------------------------------
      !
      ism1 = invs(isym)
      !
      !  the symmetry matrix in cartesian coordinates
      !  (so that we avoid going back and forth with the dynmat)
      !  note the presence of both at and bg in the transform!
      !
      scart = DBLE(s(:, :, ism1))
      scart = MATMUL(MATMUL(bg, scart), TRANSPOSE(at))
      !
      gamma = czero
      DO na = 1, nat
        !
        ! the corresponding of na in {S|v}
        sna = irt(isym, na)
        !
        ! cfac = exp[iSq*(tau_K - {S|v} tau_k)]   (Maradudin&Vosko RMP Eq. 2.33)
        ! [v can be ignored since it cancels out, see endnotes. xq is really Sq]
        ! rtau(:,isym,na) = s(:,:,invs(isym)) * tau(:, na) - tau(:,irt(isym,na))) (cartesian)
        !
        arg = twopi * DOT_PRODUCT(xq, rtau (:, isym, na))
        cfac = DCMPLX(COS(arg),-SIN(arg))
        !
        !  the submatrix (sna,na) contains the rotation scart
        !
        gamma(3 * (sna - 1) + 1:3 * sna, 3 * (na - 1) + 1:3 * na) = cfac * scart
        !
      ENDDO
      !
      !  D_{Sq} = gamma * D_q * gamma^\dagger (Maradudin & Vosko, RMP, eq. 3.5)
      !
      CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma, &
             nmodes, dynq(:, :, iq_first), nmodes, czero, dynq_tmp, nmodes)
      CALL zgemm ('n', 'c', nmodes, nmodes, nmodes, cone, dynq_tmp, &
             nmodes, gamma, nmodes, czero, dynq(:, :,current_iq), nmodes)
      !
      DO nu = 1, nmodes
        DO mu = 1, nmodes
          IF (mu /= nu .AND. ABS(dynq(mu, nu, current_iq)) > eps6 ) CALL errore &
            ('rotate_eigenm', 'problem with rotated eigenmodes', 0)
          ENDDO
      ENDDO
      !
      ! DBSP-----------------------------------------------
      !  a simple check on the frequencies
      !
      IF (iverbosity == 1) THEN
        DO na = 1, nat
          DO nb = 1, nat
            massfac = 1.d0 / DSQRT(amass(ityp(na)) * amass(ityp(nb)))
            dynq_tmp(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
            dynq(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, current_iq) * massfac
          END DO
        END DO
        !
        DO jmode = 1, nmodes
          DO imode = 1, jmode
            dynq_tmp(imode, jmode) = &
              (dynq_tmp(imode, jmode) + CONJG(dynq_tmp(jmode, imode))) / 2.d0
          ENDDO
        ENDDO
        !
        CALL cdiagh2(nmodes, dynq_tmp, nmodes, w1, cz1)
        !
        DO nu = 1, nmodes
          IF (w1(nu) > 0.d0) THEN
            wtmp(nu) =  DSQRT(ABS(w1(nu)))
          ELSE
            wtmp(nu) = -DSQRT(ABS(w1(nu)))
          ENDIF
        ENDDO
        WRITE(stdout, '(5x,"Frequencies of the matrix for the current q in the star (cm^-1)")')
        WRITE(stdout, '(6(2x,f10.5))' ) (wtmp(nu) * rydcm1, nu = 1, nmodes)
      ENDIF
    ! ZD: Save eigenfrequencies for ex-ph calculations
      IF (exciton) THEN
        DO na = 1, nat
          DO nb = 1, nat
            massfac = 1.d0 / DSQRT(amass(ityp(na)) * amass(ityp(nb)))
            dynq_tmp(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
            dynq(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, current_iq) * massfac
          END DO
        END DO
        !
        DO jmode = 1, nmodes
          DO imode = 1, jmode
            dynq_tmp(imode, jmode) = &
              (dynq_tmp(imode, jmode) + CONJG(dynq_tmp(jmode, imode))) / 2.d0
          ENDDO
        ENDDO
        !
        CALL cdiagh2(nmodes, dynq_tmp, nmodes, w1, cz1)
        !
        DO nu = 1, nmodes
          IF (w1(nu) > 0.d0) THEN
            wf_temp(nu,iq) =  DSQRT(ABS(w1(nu)))
          ELSE
            wf_temp(nu,iq) = -DSQRT(ABS(w1(nu)))
          ENDIF
        ENDDO
        ! WRITE(stdout, '(5x,"Frequencies of the matrix for the current q in the star (cm^-1)")')
        ! WRITE(stdout, '(6(2x,f10.5))' ) (wtmp(nu) * rydcm1, nu = 1, nmodes)
        !
      ENDIF
      !
      !END --------------------------------------------------
      current_iq = current_iq + 1
      !
      ! SP Repeat the same but for minus_q one
      IF (imq == 0) then
        !
        xq = -sxq(:, iq)
        ism1 = invs(isym)
        scart = DBLE(s(:, :, ism1))
        scart = MATMUL(MATMUL(bg, scart), TRANSPOSE(at))
        !
        gamma = czero
        DO na = 1, nat
          !
          sna = irt(isym, na)
          arg = twopi * DOT_PRODUCT(xq, rtau(:, isym, na))
          cfac = dcmplx(COS(arg), -SIN(arg))
          gamma(3 * (sna - 1) + 1:3 * sna, 3 * (na - 1) + 1:3 * na) = cfac * scart
          !
        ENDDO
        !
        CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma  , &
               nmodes, CONJG(dynq(:, :, iq_first)) , nmodes, czero , dynq_tmp, nmodes)
        CALL ZGEMM('n', 'c', nmodes, nmodes, nmodes, cone, dynq_tmp, &
               nmodes, gamma, nmodes, czero , dynq(:, :, current_iq), nmodes)
        !
        DO nu = 1, nmodes
          DO mu = 1, nmodes
            IF (mu /= nu .AND. ABS(dynq(mu, nu, current_iq)) > eps6 ) CALL errore &
              ('rotate_eigenm', 'problem with rotated eigenmodes', 0)
          ENDDO
        ENDDO
        current_iq = current_iq + 1
      ENDIF
    ENDDO ! iq
    !
    !---------------------------------------------------------------------------------
    END SUBROUTINE dynmat_asr
    !---------------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE loadbm()
    !----------------------------------------------------------------------------
    !!
    !! Load the information on the band manifold determined in Wannierization step
    !!
    !-----------------------------------------------------------------------
    USE input,         ONLY : filukk
    USE io_var,        ONLY : iunukk
    USE global_var,    ONLY : nbndep, nbndskip, ibndkept
    USE io_global,     ONLY : meta_ionode_id, meta_ionode
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : world_comm
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ierr
    !! Error index when reading/writing a file
    !
    IF (meta_ionode) THEN
      !
      OPEN(iunukk, FILE = filukk, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
      IF (ios /=0) CALL errore('loadbm', 'error opening ukk file', iunukk)
      !
      READ(iunukk, *) nbndep, nbndskip
      !
      ALLOCATE(ibndkept(nbndep), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadbm', 'Error allocating ibndkept', 1)
      !
      DO ibnd = 1, nbndep
        READ(iunukk, *) ibndkept(ibnd)
      ENDDO
      !
      CLOSE(iunukk)
    ENDIF ! meta_ionode
    !
    CALL mp_bcast(nbndep, meta_ionode_id, world_comm)
    !
    IF (.NOT. meta_ionode) THEN
      ALLOCATE(ibndkept(nbndep), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadbm', 'Error allocating ibndkept', 1)
    ENDIF
    !
    CALL mp_bcast(nbndskip, meta_ionode_id, world_comm)
    CALL mp_bcast(ibndkept, meta_ionode_id, world_comm)
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE loadbm
    !-----------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE loadumat(nbndep, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, &
                        exband, w_centers)
    !----------------------------------------------------------------------------
    !!
    !! Wannier interpolation of e-p vertex:
    !! load rotation matrix on coarse mesh and distribute
    !!
    !! 07/2010 JN changed the way this was done.  Previously
    !! a pool scatter call that didn't work on hbar was performed
    !! and some crazy packing scheme was needed to use poolscatter
    !! subsequently it is just a bcast followed by an appropriate assignment
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE input,         ONLY : kmap
    USE input,         ONLY : filukk
    USE ep_constants,  ONLY : czero, zero
    USE io_var,        ONLY : iunukk
    USE io_global,     ONLY : meta_ionode_id, meta_ionode
    USE mp,            ONLY : mp_sum, mp_barrier, mp_bcast
    USE parallelism,   ONLY : fkbounds
    USE global_var,    ONLY : xkq
    USE kfold,         ONLY : createkmap2
    USE pwcom,         ONLY : nbnd
    USE mp_world,      ONLY : world_comm
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(out) :: lwin(nbndep, nks)
    !! Band windows at k
    LOGICAL, INTENT(out) :: lwinq(nbndep, nks)
    !! Band windows at k+q
    LOGICAL, INTENT(out) :: exband(nbnd)
    !! Band excluded
    INTEGER, INTENT(in) :: nbndep
    !! Number of bands
    INTEGER, INTENT(in) :: nbndsub
    !! number of bands in the optimal subspace
    INTEGER, INTENT(in) :: nks
    !! number of kpoints
    INTEGER, INTENT(in) :: nkstot
    !! total number of kpoints across pools
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! the qpoint for folding of U
    REAL(KIND = DP), INTENT(inout) :: w_centers(3, nbndsub)
    !! Wannier centers
    COMPLEX(KIND = DP), INTENT(out) :: cu(nbndep, nbndsub, nks)
    !! U(k) matrix for k-points in the pool
    COMPLEX(KIND = DP), INTENT(out) :: cuq(nbndep, nbndsub, nks)
    !! U(k+q) matrix for k+q-points in the pool
    !
    ! Local variables
    LOGICAL :: lwin_big(nbndep, nkstot)
    !! .TRUE. if the band ibnd lies within the outer window at k-point ik
    LOGICAL :: lwinq_big(nbndep, nkstot)
    !! .TRUE. if the band ibnd lies within the outer window at k+qpoint ikq
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: iw
    !! Counter on Wannier centers
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on wannierized bands
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik_start
    !! Index of first k-point in the pool
    INTEGER :: ik_stop
    !! Index of last k-point in the pool
    COMPLEX(KIND = DP), ALLOCATABLE :: cu_big(:, :, :)
    !! U(k) matrix for all k-points
    COMPLEX(KIND = DP), ALLOCATABLE :: cuq_big(:, :, :)
    !! U(k+q) matrix for all k+q-points
    !
    ALLOCATE(cu_big(nbndep, nbndsub, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadumat', 'Error allocating cu_big', 1)
    ALLOCATE(cuq_big(nbndep, nbndsub, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadumat', 'Error allocating cuq_big', 1)
    !
    cu_big = czero
    cuq_big = czero
    IF (meta_ionode) THEN
      !
      ! First proc read rotation matrix (coarse mesh) from file
      !
      OPEN(iunukk, FILE = filukk, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
      IF (ios /=0) CALL errore('loadumat', 'error opening ukk file', iunukk)
      !
      ! dummy operation for skipping unnecessary data (nbndep, nbndskip & ibndkept) here
      !
      READ(iunukk, *) ibnd, jbnd
      DO ibnd = 1, nbndep
        READ(iunukk, *) jbnd
      ENDDO
      !
      DO ik = 1, nkstot
        DO ibnd = 1, nbndep
          DO jbnd = 1, nbndsub
            READ(iunukk, *) cu_big(ibnd, jbnd, ik)
          ENDDO
        ENDDO
      ENDDO
      DO ik = 1, nkstot
        DO ibnd = 1, nbndep
          READ(iunukk, *) lwin_big(ibnd, ik)
        ENDDO
      ENDDO
      DO ibnd = 1, nbnd
        READ(iunukk, *) exband(ibnd)
      ENDDO
      ! Read the Wannier centers
      DO iw = 1, nbndsub
        READ(iunukk, *) w_centers(:, iw)
      ENDDO
      !
      CLOSE(iunukk)
      !
      ! Generate U(k+q) through the map
      ! here we create the map k+q --> k
      ! (first proc has a copy of all kpoints)
      ! Generates kmap(ik) for this xxq
      !
      xkq(:, :) = zero
      CALL createkmap2(xxq)
      !
      ! and we generate the matrix for the q-displaced mesh
      !
      DO ik = 1, nkstot
        cuq_big(:, :, ik) = cu_big(:, :, kmap(ik))
        lwinq_big(:, ik) = lwin_big(:, kmap(ik))
      ENDDO
    ENDIF ! meta_ionode
    !
    CALL mp_bcast(cu_big, meta_ionode_id, world_comm)
    CALL mp_bcast(cuq_big, meta_ionode_id, world_comm)
    CALL mp_bcast(lwin_big, meta_ionode_id, world_comm)
    CALL mp_bcast(lwinq_big, meta_ionode_id, world_comm)
    CALL mp_bcast(exband, meta_ionode_id, world_comm)
    CALL mp_bcast(w_centers, meta_ionode_id, world_comm)
    !
    CALL fkbounds(nkstot, ik_start, ik_stop)
    !
    IF ((ik_stop - ik_start + 1) /= nks) CALL errore('loadumat', "Improper parallel ukk load", 1)
    !
    cu = cu_big(:, :, ik_start:ik_stop)
    cuq = cuq_big(:, :, ik_start:ik_stop)
    lwin = lwin_big(:, ik_start:ik_stop)
    lwinq = lwinq_big(:, ik_start:ik_stop)
    !
    DEALLOCATE(cu_big, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadumat', 'Error deallocating cu_big', 1)
    DEALLOCATE(cuq_big, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadumat', 'Error deallocating cuq_big', 1)
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE loadumat
    !-----------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE io
  !------------------------------------------------------------------------------
