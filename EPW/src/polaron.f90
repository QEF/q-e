  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !-----------------------------------------------------------------------
  MODULE polaron
  !-----------------------------------------------------------------------
  !!
  !! This module contains variables and subroutines of polaron
  !! Authored by Chao Lian, Weng Hong (Denny) Sio, and Jon Lafuente-Bartolome
  !! Partial cleaning by SP (Nov 2023)
  !! Partial cleaning by STiwari (Nov 2023)
  !! Partial cleaning by JLB (Aug 2024)
  !! Adding variables by KL (Oct 2024)
  !! Optimization by DK, TYK, JLB (May 2025)
  !!
  USE kinds,     ONLY : DP
  USE buffers,   ONLY : open_buffer, get_buffer, save_buffer, close_buffer
  USE io_var,    ONLY : iepfall, ihamil, iMmn, irho, iUmn, iekanu
  !
  IMPLICIT NONE
  !
  LOGICAL :: test_tags_plrn(20) = .FALSE.
  !! The B matrix Bqu
  LOGICAL :: mem_save_h = .FALSE.
  !! The B matrix Bqu
  LOGICAL, ALLOCATABLE :: is_mirror_k(:)
  !! .true. if k is a mirror point, used for time-reversal symmetry
  LOGICAL, ALLOCATABLE :: is_mirror_q(:)
  !! .true. if q is a mirror point, used for time-reversal symmertry
  LOGICAL, ALLOCATABLE :: is_tri_k(:)
  !! .true. if k is a time-reversal invariant point, used for time-reversal symmetry
  LOGICAL, ALLOCATABLE :: is_tri_q(:)
  !! .true. if q is a time-reversal invariant point, used for time-reversal symmetry
  INTEGER :: nbnd_plrn
  !! Number of bands used in polaron calculations
  INTEGER :: nbnd_g_plrn
  !! Number of bands in which g is to be interpolated
  INTEGER :: lword_h
  !! Hamiltonian record size for I/O
  INTEGER :: lword_g
  !! el-ph matrix element record size for I/O
  INTEGER :: lword_m
  !! FIXME
  INTEGER :: io_level_g_plrn
  !! Write el-ph matrix elements to disk or store in memory
  INTEGER :: io_level_h_plrn
  !! Write Hamiltonian to disk or store in memory
  INTEGER :: hblocksize
  !! FIXME
  INTEGER :: band_pos
  !! Band in which CBM or VBM is located
  INTEGER :: ik_edge
  !! k-point in which CBM or VBM is located
  INTEGER :: nRp
  !! Number of unit cells on supercell for non-diagonal supercells
  INTEGER, ALLOCATABLE :: Rp(:,:)
  !! List of unit cell vectors within supercell
  INTEGER, ALLOCATABLE :: select_bands_plrn(:)
  !! Map from {start_band_plrn, end_band_plrn} to {1, end_band_plrn - start_band_plrn}
  INTEGER, ALLOCATABLE :: kpg_map(:)
  !! Map from a given k1 to its mirror point k2 = -k1 + G
  REAL(KIND = DP) :: wq_model
  !! Phonon freq in Frohlich model
  REAL(KIND = DP), ALLOCATABLE :: etf_model(:)
  !! Band structure in Frohlich model
  REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
  !! Gathered KS eigenvalues over the pools
  REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
  !! Gathered k-point coordinates over the pools
  COMPLEX(KIND = DP), ALLOCATABLE :: Hamil(:, :)
  !! Effective polaron Hamiltonian
  COMPLEX(KIND = DP), ALLOCATABLE :: eigvec(:, :)
  !! Polaron wave function coefficients in Bloch basis, Ank
  COMPLEX(KIND = DP), ALLOCATABLE :: Bmat(:,:)
  !! Polaron displacement coefficients in phono basis, Bqv
  COMPLEX(KIND = DP), ALLOCATABLE :: gq_model(:)
  !! el-ph matrix element in a simplified Fr\"ohlich model
  COMPLEX(KIND = DP), ALLOCATABLE :: epf(:, :, :, :)
  !! el-ph matrix element in for a given q
  COMPLEX(KIND = DP), ALLOCATABLE :: epfall(:, :, :, :, :)
  !! el-ph matrix element for all q
  COMPLEX(KIND = DP) :: berry_phase(1:3)
  !! FIXME
  PUBLIC  :: plrn_prepare
  !! Subroutine preparing variables for polaron calculations
  PUBLIC  :: plrn_flow_select
  !! Subroutine selecting polaron scf or post-processing
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE plrn_prepare(totq, iq_restart)
    !-----------------------------------------------------------------------
    !!
    !! Routine to prepare quantities for polaron calculation
    !!
    !-----------------------------------------------------------------------
    !
    USE input,         ONLY : start_band_plrn, end_band_plrn, nbndsub, nstate_plrn, debug_plrn, &
                              cal_psir_plrn, restart_plrn,  interp_Ank_plrn, interp_Bqu_plrn,   &
                              model_vertex_plrn, nhblock_plrn, g_start_band_plrn,               &
                              g_end_band_plrn, lrot, lphase, type_plrn, g_start_energy_plrn,    &
                              g_end_energy_plrn, model_enband_plrn, model_vertex_plrn,          &
                              g_power_order_plrn, io_lvl_plrn, nkf1, nkf2, nkf3, m_eff_plrn,    &
                              kappa_plrn, omega_LO_plrn, lfast_kmesh, scell_mat_plrn
    USE global_var,    ONLY : nkqf, nkf, nqf, nqtotf, nktotf, etf, xkf, xqf
    USE modes,         ONLY : nmodes
    USE cell_base,     ONLY : bg, omega, alat
    USE ep_constants,  ONLY : czero, cone, pi, ci, twopi, fpi, eps6, eps8, eps5, zero, ryd2ev
    USE parallelism,   ONLY : poolgather2
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE io_files,      ONLY : check_tempdir
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out)  :: iq_restart
    !! Restart q-point indexx
    INTEGER, INTENT(in)   :: totq
    !! Total q-point in the fsthick
    !
    ! Local variables
    LOGICAL :: debug
    !! Debug flag
    LOGICAL :: plrn_scf
    !! .true. if self-consistent polaron calculation is to be performed
    LOGICAL :: exst
    !! Checking on the file presence
    LOGICAL :: pfs
    !! FIXME
    INTEGER :: iq
    !! q-point index
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: ikq
    !! k+q point index
    INTEGER :: ik_global
    !! k-point index in global list
    INTEGER :: ierr
    !! Error status
    INTEGER :: ikGamma
    !! Index of \Gamma point in k-point list
    INTEGER :: iqGamma
    !! Index of \Gamma point in q-point list
    INTEGER :: minNBlock
    !! FIXME
    INTEGER :: ishift
    !! Index of neighbor G-vectors
    INTEGER :: lword_h_tmp
    !! Temporary Hamiltonian record size for I/O
    INTEGER :: lword_g_tmp
    !! Temporary el-ph matrix element record size for I/O
    INTEGER, PARAMETER :: maxword = HUGE(1)
    !! FIXME
    REAL(KIND = DP) :: xxk(3)
    !! k-point coordinates
    REAL(KIND = DP) :: efermi
    !! Fermi level (VBM or CBM)
    REAL(KIND = DP) :: klen
    !! Distance to the nearest \Gamma point
    REAL(KIND = DP) :: shift(3)
    !! Shift G-vector to find nearest \Gamma point
    REAL(KIND = DP) :: rfac
    !! Numerator of matrix element in Frohlich model
    REAL(KIND = DP), ALLOCATABLE :: rtmp2(:,:)
    !! Temporary array for gathering eigenvalues across pools
    !
    CALL start_clock('plrn_prepare')
    !
    IF(lfast_kmesh) THEN
      CALL errore('polaron_prepare', 'Polaron module not working with lfast_kmesh', 1)
    ENDIF
    !
    WRITE(stdout, '(5x,"fsthick not working in polaron module, selecting all the k/q points.")')
    !! type_plrn denotes whether electron polaron (-1) or hole polaron (+1)
    !! Legalize the type_plrn input,     in case that the user use an arbitrary number
    IF(type_plrn < 0) THEN
      type_plrn = -1
      WRITE(stdout, '(5x, "The electron polaron is calculated.")')
    ELSE
      type_plrn = 1
      WRITE(stdout, '(5x, "The hole polaron is calculated.")')
    ENDIF
    !
    lrot = .TRUE.
    lphase = .TRUE.
    !
    debug = debug_plrn
    IF (debug_plrn) CALL check_tempdir('test_out', exst, pfs)
    !
    WRITE(stdout,'(5x,a)') REPEAT('=',67)
    !
    IF (g_start_band_plrn == 0) g_start_band_plrn = 1
    IF (g_end_band_plrn == 0) g_end_band_plrn = nbndsub
    nbnd_g_plrn = g_end_band_plrn - g_start_band_plrn + 1
    !
    IF (start_band_plrn == 0) start_band_plrn = g_start_band_plrn
    IF (end_band_plrn == 0) end_band_plrn = g_end_band_plrn
    nbnd_plrn = end_band_plrn - start_band_plrn + 1
    !
    IF(g_start_band_plrn > start_band_plrn .OR. g_end_band_plrn < end_band_plrn) THEN
      CALL errore('polaron_prepare', 'Selecting more bands in polaron than saving g matrix', 1)
    ENDIF
    !
    ALLOCATE(select_bands_plrn(nbnd_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_prepare', 'Error allocating select_bands_plrn', 1)
    !
    select_bands_plrn = 0
    DO ibnd = 1, nbnd_plrn
      select_bands_plrn(ibnd) = start_band_plrn + ibnd - 1
    ENDDO
    !
    !! copy q(x,y,z) to xkf_all, save the copy of all kpoints
    !! Note that poolgather2 has the dimension of nktotf*2,
    !! which has k at ik and k+q at ik+1
    !! This is because that xkf has the dimension of 3, nkf*2
    !! where the ik is k and ik+1 is k+q
    ALLOCATE(xkf_all(3, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating xkf_all', 1)
    ALLOCATE(rtmp2(3, nktotf * 2), STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating rtmp2', 1)
    ALLOCATE(epf(nbnd_g_plrn, nbnd_g_plrn, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating epf', 1)
    epf = czero
    !
    xkf_all = zero
    rtmp2 = zero
    CALL poolgather2(3, nktotf * 2, nkqf, xkf, rtmp2)
    xkf_all(1:3, 1:nktotf) = rtmp2(1:3, 1:nktotf * 2:2)
    !
    DEALLOCATE(rtmp2, STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error deallocating rtmp2', 1)
    !
    WRITE(stdout, "(5x, 'Use the band from ',i0, ' to ', i0, ' total ', i0)") start_band_plrn, end_band_plrn, nbnd_plrn
    WRITE(stdout, "(5x, 'Including bands: ', 10i3)") select_bands_plrn
    WRITE(stdout, "(5x, 'Use the band from ',i0, ' to ', i0, ' total ', i0, ' in saving g')") &
          g_start_band_plrn, g_end_band_plrn, nbnd_g_plrn
    !
    WRITE(stdout, "(5x, 'Gathering eigenvalues of ', i0, ' bands and ', i0, ' k points')") nbndsub, nktotf
    !
    ALLOCATE(etf_all(nbndsub, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating etf_all', 1)
    !
    IF(model_enband_plrn) THEN
      etf = zero
      ! Find the distance to the nearest Gamma point in crystal coordinates
      DO ik = 1, 2 * nkf
        klen = 1E3
        xxk = xkf(:, ik)
        CALL cryst_to_cart(1, xxk, bg, 1)
        DO ishift = 1, 27
          shift(1:3) = REAL(index_shift(ishift), KIND = DP)
          xxk = xkf(:, ik) + shift
          CALL cryst_to_cart(1, xxk, bg, 1)
          klen = MIN(klen, NORM2(xxk))
        ENDDO
        etf(1, ik) = 0.5 / m_eff_plrn * (klen * twopi / alat)**2
      ENDDO
    ENDIF
    !
    etf_all = zero
    CALL gather_band_eigenvalues(etf, etf_all)
    !
    IF(model_vertex_plrn) THEN
      WRITE(stdout, '(5x, a, f8.3)') "Using model g vertex, with order ", g_power_order_plrn
      ALLOCATE(gq_model(nqf), STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating gq_model', 1)
      gq_model = zero
      !
      rfac = SQRT(fpi / omega * omega_LO_plrn / kappa_plrn)
      DO iq = 1, nqf
        klen = 1E3
        DO ishift = 1, 27
          shift(1:3) = REAL(index_shift(ishift), KIND = DP)
          xxk = xqf(:, iq) + shift
          CALL cryst_to_cart(1, xxk, bg, 1)
          klen = MIN(klen, NORM2(xxk))
        ENDDO
        IF(klen > eps8) THEN
          gq_model(iq) = rfac / ((klen * twopi / alat) ** g_power_order_plrn)
        ENDIF
      ENDDO
    ENDIF

    ! change unit from eV to Rydberg
    g_start_energy_plrn = g_start_energy_plrn / ryd2ev
    g_end_energy_plrn = g_end_energy_plrn / ryd2ev
    !
    CALL start_clock('find_EVBM')
    CALL find_band_extreme(type_plrn, etf_all, ik_edge, band_pos, efermi)
    !
    ! Determine the Fermi energy, read from the input or calculated from band structure
    WRITE(stdout, '(5x, "Fermi Energy is", f16.7, &
    &" (eV) located at kpoint ", i6, 3f8.3, " band ", i3)') efermi * ryd2ev, ik_edge, xkf_all(1:3, ik_edge), band_pos
    ! Shift the eigenvalues to make VBM/CBM zero
    etf_all(1:nbndsub, 1:nktotf) = etf_all(1:nbndsub, 1:nktotf) - efermi
    !
    CALL stop_clock('find_EVBM')
    !
    WRITE(stdout, "(5x, 'Allocating arrays and open files.')")
    !
    IF(interp_Ank_plrn .OR. interp_Bqu_plrn .OR. cal_psir_plrn) THEN
      plrn_scf = .FALSE.
      restart_plrn = .TRUE.
    ELSE
      plrn_scf = .TRUE.
    ENDIF
    !
    IF(restart_plrn) THEN
      iq_restart = totq + 1
    ELSE
      iq_restart = 1
    ENDIF
    !
    IF(interp_Ank_plrn) THEN
      ALLOCATE(eigvec(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating eigvec', 1)
      eigvec = czero
    ELSE IF(plrn_scf) THEN
       CALL check_tempdir('plrn_tmp', exst, pfs)
       !
       io_level_g_plrn = 1
       lword_g_tmp = nbnd_g_plrn * nbnd_g_plrn * nmodes * nkf
       IF(lword_g_tmp > maxword) THEN
         CALL errore('plrn_prepare', 'Record size larger than maximum, use more cores!', 1)
       ELSE
         lword_g = INT(lword_g_tmp)
       ENDIF
       !
       IF (io_lvl_plrn == 0) THEN
         ALLOCATE(epfall(nbnd_g_plrn, nbnd_g_plrn, nmodes, nkf, nqtotf), STAT = ierr)
         IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating epfall', 1)
         epfall = czero
       ELSE IF (io_lvl_plrn == 1) THEN
         CALL open_buffer( iepfall , 'ephf' , lword_g, io_level_g_plrn, exst, direc = 'plrn_tmp/')
       ENDIF
       !
       io_level_h_plrn = 1
       IF(nhblock_plrn < 1 .OR. nhblock_plrn > nkf * nbnd_plrn) THEN
         CALL errore('plrn_prepare','Illegal nhblock_plrn, should between 1 and nkf * nbnd_plrn', 1)
       ENDIF
       !
       minNBlock = CEILING(REAL(nkf * nbnd_plrn * nktotf * nbnd_plrn, dp) / maxword)
       !
       IF(minNBlock >  nhblock_plrn .AND. nhblock_plrn /= 1) THEN
         CALL errore('plrn_prepare', 'Record size larger than maximum, use more cores!', 1)
       ENDIF
       !
       hblocksize = CEILING(REAL(nkf * nbnd_plrn, dp) / nhblock_plrn)
       !
       lword_h_tmp = nktotf * nbnd_plrn * hblocksize
       !
       IF(nhblock_plrn /= 1) THEN
         IF(lword_h_tmp > maxword) THEN
           CALL errore('plrn_prepare', 'Record size larger than maximum, use more cores or larger nhblock_plrn!', 1)
         ELSE
           lword_h = lword_h_tmp
         ENDIF
       ENDIF
       !
       lword_m = nbnd_plrn * nbnd_plrn * nktotf * 3
       !
       ALLOCATE(Hamil(nktotf * nbnd_plrn, hblocksize), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating Hamil', 1)
       !
       ! Allocate and initialize the variables
       ALLOCATE(eigvec(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating eigvec', 1)
       eigvec = czero
       !
       ! Check whether the input is legal, otherwise print warning and stop
       ! Check the input now, because if inputs are illegal, we can stop the calculation
       ! before the heavy el-ph interpolation begins.
       IF(nktotf /= nqtotf .OR. nktotf < 1) CALL errore('plrn_prepare','Not identical k and q grid. Do use same nkf and nqf!', 1)
       !
       IF( (.NOT. scell_mat_plrn) .AND. (nkf1 == 0 .or. nkf2 == 0 .or. nkf3 == 0) ) THEN
          CALL errore('plrn_prepare','Try to use nkf and nqf to generate k and q grid, &
             &IF you are using a manual grid, also provide this information.', 1)
       ENDIF
       !
       IF(nkf < 1) CALL errore('plrn_prepare','Some node has no k points!', 1)
       IF(nqtotf /= nqf) CALL errore('plrn_prepare','Parallel over q is not available for polaron calculations.', 1)
       !
       WRITE(stdout, '(5x, "Polaron wavefunction calculation starts with k points ",&
          &i0, ", q points ", i0, " and KS band ", i0)') nktotf,  nqtotf,  nbnd_plrn
       !
       ! check whether the k and q mesh are identical
       ! This may not be theoretically necessary,
       ! but necessary in this implementation
       DO ik = 1, nkf
         ik_global = ikqLocal2Global(ik, nktotf)
         IF (ANY(ABS(xkf_all(1:3, ik_global) - xqf(1:3, ik_global)) > eps6)) THEN
           CALL errore('plrn_prepare', 'The k and q meshes must be exactly the same!', 1)
         ENDIF
       ENDDO
       !
       ! map iq to G-iq, ik to G-ik
       ! find the position of Gamma point in k and q grid
       ! Note that xqf is not a MPI-local variable, xqf = xqtotf otherwise
       ! the program gives wrong results
       ikGamma = indexGamma(xkf_all)
       iqGamma = indexGamma(xqf)
       IF (ikGamma == 0) CALL errore('plrn_prepare','k = 0 not included in k grid!', 1)
       IF (iqGamma == 0) CALL errore('plrn_prepare','q = 0 not included in q grid!', 1)
       WRITE(stdout, '(5x, "The index of Gamma point in k grid is ", i0, " &
          &and q grid IS ", i0)') ikGamma, iqGamma
       !
       ! Given k, find the index of -k+G, to impose the relation
       ! A_{n,-k+G} = A^*_{n,k} and B_{-q+G, \nu} = B^*_{q, \nu}
       ! the relation k1 + k2 = G is unique given k1
       ! Since both A and B have the dimension of nktotf/nqtotf,
       ! kpg_map should map all the k from 1 to nktotf (global)
       ALLOCATE(kpg_map(nqtotf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating kpg_map', 1)
       kpg_map = 0
       !
       ALLOCATE(is_mirror_k(nkf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_k', 1)
       is_mirror_k = .FALSE.
       ALLOCATE(is_mirror_q(nqf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       is_mirror_q = .FALSE.
       ALLOCATE(is_tri_k(nkf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       is_tri_k = .FALSE.
       ALLOCATE(is_tri_q(nqf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       is_tri_q = .FALSE.
       !
       WRITE(stdout, '(5x, a)') "Finding the index of -k for each k point."
       ! For two k points k1 and k2, if k1 = G - k2, G is any reciprocal vector
       ! then k2 is the mirror point of k1 if the index of k2 is larger than k1
       ! Same rule for q, while is_mirror_q(nqf) is global but is_mirror_k(nkf) is local
       DO ik = 1, nkf
         ik_global = ikqLocal2Global(ik, nktotf)
         DO ikq = 1, nktotf
           ! -k+G = k', i.e. k' + k = G, G may be (0, 0, 0)
           xxk = xkf_all(1:3, ik_global) + xkf_all(1:3, ikq)
           IF (isGVec(xxk)) kpg_map(ik_global) = ikq
         ENDDO
         IF (kpg_map(ik_global) == 0) CALL errore('plrn_prepare', 'Not legal k/q grid!', 1)
         !
         ikq = kpg_map(ik_global)
         IF (ik_global > ikq) THEN
           is_mirror_k(ik) = .TRUE.
         ELSE IF (ik == ikq) THEN
           is_tri_k(ik) = .TRUE.
         ENDIF
       ENDDO
       !
       CALL mp_sum(kpg_map, inter_pool_comm)
       !
       DO iq = 1, nqf
         ikq = kpg_map(iq)
         IF (iq > ikq) THEN
           is_mirror_q(iq) = .TRUE.
         ELSE IF (iq == ikq) THEN
           is_tri_q(iq) = .TRUE.
         ENDIF
       ENDDO
       !
       WRITE(stdout, '(5x, a)') "Checking the k + q is included in the mesh grid for each k and q."
       ! find the global index of ik_global, ikq with vector k and k+q.
       ! Different from kpg_map, it is used in constructing Hamiltonian or hpsi
       ! ik goes over all the local k points, to parallel the program
       DO ik = 1, nkf
         DO iq = 1, nqtotf
           ik_global = ikqLocal2Global(ik, nktotf)
           xxk = xkf_all(1:3, iq) + xkf_all(1:3, ik_global)
           IF (ikq_all(ik, iq) == 0) THEN
             CALL errore('plrn_prepare','Not commensurate k and q grid!', 1)
           ENDIF
         ENDDO
       ENDDO
    ENDIF
    !
    WRITE(stdout, "(5x, 'End of plrn_prepare')")
    !
    CALL stop_clock('plrn_prepare')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE plrn_prepare
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE plrn_save_g_to_file(iq, epfg, wfreq)
    !-----------------------------------------------------------------------
    !!
    !! Save el-ph matrix element to file
    !!
    USE modes,         ONLY : nmodes
    USE input,         ONLY : g_start_band_plrn, g_end_band_plrn, &
                              io_lvl_plrn, eps_acoustic, g_scale_plrn
    USE global_var,    ONLY : nkf
    USE ep_constants,  ONLY : eps8, two, czero
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! q-point index
    REAL(KIND = DP), INTENT(in) :: wfreq(:, :)
    !! Phonon frequency
    COMPLEX(KIND = DP), INTENT(inout) :: epfg(:, :, :, :)
    !! el-ph matrix element
    !
    ! Local variables
    INTEGER :: ik
    !! k-point index
    INTEGER :: ikq
    !! k+q point index
    INTEGER :: imode
    !! mode index
    !
    ! In polaron equations, g is not epf but epf/omega
    ! To ensure a Hermitian Hamiltonian, g_{mnu}(k, -q) is calculated as g*_{nmu}(k-q, q)
    DO ik = 1, nkf
      ikq = ikq_all(ik, iq)
      DO imode = 1, nmodes
        IF (wfreq(imode, iq) > eps_acoustic) THEN
          !epf(:, :, imode, ik) = &
          epf(:, :, imode, ik) = g_scale_plrn * &
          epfg(g_start_band_plrn:g_end_band_plrn, g_start_band_plrn:g_end_band_plrn, imode, ik) / DSQRT(two * wfreq(imode, iq))
        ENDIF
      ENDDO ! imode
    ENDDO ! ik
    !
    IF(io_lvl_plrn == 0) THEN
      epfall(:, :, : ,: ,iq) = epf(:, :, :, :)
    ELSE IF (io_lvl_plrn == 1) THEN
      CALL save_buffer(epf(:, :, :, :), nbnd_g_plrn * nbnd_g_plrn * nmodes * nkf, iepfall, iq)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE plrn_save_g_to_file
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION ikq_all(ik, iq)
    !-----------------------------------------------------------------------
    !!
    !! find the global index of k+q for the local ik and global iq
    !!
    USE global_var, ONLY : nktotf
    USE input,      ONLY : nkf1, nkf2, nkf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! k-point index
    INTEGER, INTENT(in) :: iq
    !! q-point index
    !
    ! Local variables
    INTEGER :: ikq
    ! k+q index
    INTEGER :: ik_global
    !! k-point index in global list
    INTEGER :: ikq_all
    !! k+q point index in global list
    INTEGER :: index_target(1:3)
    !! Auxiliary index
    INTEGER :: index_kq
    !! Auxiliary index of k+q point
    INTEGER :: ikq_loop
    !! Loop counter
    REAL(KIND = DP) :: xxk(1:3)
    !! k+q point coordinates
    REAL(KIND = DP) :: xxk_target(1:3)
    !! k+q point coordinates in 1BZ
    !
    ikq_all = 0
    !
    ik_global = ikqLocal2Global(ik, nktotf)
    xxk = xkf_all(1:3, iq) + xkf_all(1:3, ik_global)
    !
    xxk_target(1:3) = xxk(1:3) - INT(xxk(1:3))
    index_target(1:3) = NINT(xxk_target(1:3) * (/nkf1, nkf2, nkf3/))
    !
    index_kq = index_target(1) * nkf1 * nkf2 + index_target(2) * nkf2 + index_target(3) + 1
    !
    DO ikq_loop = index_kq - 1, nktotf + index_kq
      ! ik (local) + iq (global) = ikq (global)
      ! get ikq to locate the column of the Hamiltonian
      ikq = MOD(ikq_loop, nktotf) + 1
      IF (isGVec(xxk - xkf_all(1:3, ikq))) THEN
        ikq_all = ikq
        EXIT
      ENDIF
    ENDDO
    !
    IF (ikq_all == 0) CALL errore('ikq_all','k + q not found', 1)
    !
    !-----------------------------------------------------------------------
    END FUNCTION ikq_all
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION find_ik(xxk, xkf_global)
    !-----------------------------------------------------------------------
    !!
    !! Find k-point index
    !!
    USE global_var, ONLY : nktotf
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xxk(1:3)
    !! K-point position per cpu
    REAL(KIND = DP), INTENT(in) :: xkf_global(1:3, 1:nktotf)
    !! global k-points
    !
    ! Local variables
    INTEGER :: ik
    !! k-point index
    INTEGER :: find_ik
    !! k-point index of found k-point
    REAL(KIND = DP) :: xkq(1:3)
    !! k-point coordinates
    !
    CALL start_clock('find_k')
    !
    find_ik = 0
    DO ik = 1, nktotf
      xkq(1:3) = xkf_global(1:3, ik) - xxk(1:3)
      IF(isGVec(xkq)) THEN
        find_ik = ik
        EXIT
      ENDIF
    ENDDO
    !
    IF (find_ik == 0) CALL errore('find_ik','k not found', 1)
    !
    CALL stop_clock('find_k')
    !
    !-----------------------------------------------------------------------
    END FUNCTION find_ik
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE plrn_flow_select(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
    !-----------------------------------------------------------------------
    !!
    !! Driver which selects whether a self-consistent polaron
    !! or a post-processing calculation is to be performed.
    !!
    !-----------------------------------------------------------------------
    USE input,         ONLY : cal_psir_plrn,  interp_Ank_plrn, interp_Bqu_plrn, &
                              io_lvl_plrn, scell_mat_plrn
    USE io_global,     ONLY : stdout, ionode
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT (in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT (in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT (in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT (in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT (in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT (in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT (in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), INTENT (in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    !
    ! Local variable
    !
    ! Bqu Ank interpolation is not compatible with self-consistency process
    ! Added by Chao Lian for polaron calculations flow select
    ! If postprocess is ON, i.e. Bqu interpolation with saved dtau,
    ! Ank interpolation with saved Amp, and polaron visualization with saved Wannier function cube files,
    ! then self-consistent process is skipped.
    IF (.NOT. (interp_Bqu_plrn .OR. interp_Ank_plrn .OR. cal_psir_plrn)) THEN
      CALL polaron_scf(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
      IF (io_lvl_plrn == 1) CALL close_buffer(iepfall, 'KEEP')
      CALL close_buffer(ihamil, 'delete')
    ENDIF
    !
    IF (interp_Ank_plrn) THEN
      IF(ionode) WRITE(stdout, "(5x, 'Interpolating the Ank at given k-point set....')")
      CALL interp_plrn_wf(nrr_k, ndegen_k, irvec_r, dims)
    ENDIF
    !
    IF (interp_Bqu_plrn) THEN
      IF(ionode) THEN
        WRITE(stdout, "(5x, 'Interpolating the Bqu at given q-point set....')")
      ENDIF
      CALL interp_plrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws)
    ENDIF
    !
    IF (cal_psir_plrn) THEN
      IF(ionode) WRITE(stdout, "(5x, 'Calculating the real-space distribution of polaron wavefunction....')")
      IF (scell_mat_plrn) THEN
        CALL scell_write_real_space_wavefunction()
      ELSE
        CALL write_real_space_wavefunction()
      ENDIF
    ENDIF
    !
    ! Deallocate allocated arrays and close open files
    CALL plrn_close()
    !
    WRITE(stdout, '(/5x, "======================== Polaron Timers ===========================")')
    CALL print_clock('main_prln')
    CALL print_clock('plrn_prepare')
    CALL print_clock('write_files')
    CALL print_clock('Bqu_tran')
    CALL print_clock('Ank_trans')
    CALL print_clock('cal_E_Form')
    CALL print_clock('DiagonH')
    CALL print_clock('Setup_H')
    CALL print_clock('H_alloc')
    CALL print_clock('read_gmat')
    CALL print_clock('read_Hmat')
    CALL print_clock('Write_Hmat')
    CALL print_clock( 'cegterg' )
    CALL print_clock( 'cegterg:init' )
    CALL print_clock( 'cegterg:diag' )
    CALL print_clock( 'cegterg:update' )
    CALL print_clock( 'cegterg:overlap' )
    CALL print_clock( 'cegterg:last' )
    CALL print_clock('cal_bqu')
    CALL print_clock('init_Ank')
    CALL print_clock('find_EVBM')
    CALL print_clock('re_omega')
    CALL print_clock('cal_hpsi')
    WRITE(stdout, '(5x, "===================================================================")')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE plrn_flow_select
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE polaron_scf (nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
    !-----------------------------------------------------------------------
    !!
    !! Self consistency calculation of polaron wavefunction.
    !! Rewritten by Chao Lian based on the implementation by Denny Sio.
    !! SP: cleaning (Nov 2023)
    !!
    !
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : ryd2mev, one, ryd2ev, two, zero, twopi,           &
                              czero, cone, pi, ci, twopi, eps6, eps8, eps5
    USE input,         ONLY : type_plrn, full_diagon_plrn, debug_plrn,          &
                              init_sigma_plrn, init_k0_plrn, nstate_plrn,       &
                              conv_thr_plrn, init_plrn, niter_plrn,             &
                              nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, r0_plrn,      &
                              init_ntau_plrn, nbndsub, as, time_rev_A_plrn,     &
                              model_phfreq_plrn, omega_LO_plrn, scell_mat_plrn, &
                              acoustic_plrn, cal_acous_plrn, dtau_max_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE global_var,    ONLY : nqtotf, nktotf, wf
    USE cell_base,     ONLY : alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE parallelism,   ONLY : poolgather2
    USE mp_world,      ONLY : world_comm
    USE input,         ONLY : ethrdg_plrn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT(in) ::  dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    !
    ! local variables
    CHARACTER(LEN = 256) :: filename
    !! Output file name
    CHARACTER(LEN = 256) :: tmpch
    !! Temporary character to assign name to different displacement files
    LOGICAL :: debug
    !! .true. if extra output is to be printed for debugging
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: iq
    !! q-point counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ibnd
    !! band counter
    INTEGER :: ierr
    !! Error status
    INTEGER :: itau
    !! Atom index
    INTEGER :: iter
    !! Iteration counter
    INTEGER :: indexkn1
    !! Combined k-point and band index
    INTEGER :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER :: nkf3_p
    !! Fine k-point grid along b3
    INTEGER :: nbndsub_p
    !! Number of bands in polaron calculation
    INTEGER :: nktotf_p
    !! Number of k-points in the fine grid
    REAL(KIND = DP) :: estmteRt(nstate_plrn)
    !! Polaron eigenvalue in diagonalization
    REAL(KIND = DP) :: eigval(nstate_plrn)
    !! Polaron eigenvalue
    REAL(KIND = DP) :: esterr
    !! Difference of displacements between scf loops
    REAL(KIND = DP) :: eplrnelec
    !! Electron part of polaron formation energy
    REAL(KIND = DP) :: eplrnphon
    !! Phonon part of polaron formation energy
    REAL(KIND = DP) :: r_cry(1:3)
    !! Polaron center
    REAL(KIND = DP) :: dtau_diff
    !! Max difference between displacements between scf loops
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! Polaron wave function coefficients in Wannier basis, Amp
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! Polaron displacements
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau_acoustic(:, :)
    !! Polaron displacements by acoustic phonon modes
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau_save(:, :)
    !! Polaron displacements from previous iteration to compare
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau_list(:, :, :)
    !! List of displacements from which scf is to be initiated
    !
    ALLOCATE(dtau(nktotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau', 1)
    ALLOCATE(dtau_save(nktotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau_save', 1)
    dtau = czero
    dtau_save = czero
    debug = debug_plrn
    IF (cal_acous_plrn) THEN
      ALLOCATE(dtau_acoustic(nktotf, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau_acoustic', 1)
      dtau_acoustic = czero
    ENDIF
    !
    CALL start_clock('main_prln')
    ! Gather all the eigenvalues to determine the EBM/VBM,
    CALL start_clock('re_omega')
    ! Recalculate the frequency, when restart from save g
    CALL cal_phonon_eigenfreq(nrr_q, irvec_q, ndegen_q, rws, nrws, wf)
    !
    IF(model_phfreq_plrn) THEN
      wf = zero
      wf(nmodes, :) = omega_LO_plrn
    ENDIF
    !
    CALL stop_clock('re_omega')
    !
    !! Initialize Ac(k) based on profile
    !! TODO: ik_bm should be user-adjustable
    CALL start_clock('init_Ank')
    eigvec = czero
    SELECT CASE (init_plrn)
      CASE (1)
        ! If k0 has not been set on input, center gaussian at band edge
        IF (ALL(init_k0_plrn(:) == 1000.d0)) init_k0_plrn = xkf_all(1:3, ik_edge)
        !
        WRITE(stdout, '(5x, "Initializing polaron wavefunction using Gaussian wave &
           &packet with a width of", ES14.6)') init_sigma_plrn
        WRITE(stdout, '(5x, "centered at k=", 3f14.6)') init_k0_plrn !xkf_all(1:3, ik_edge)
        CALL init_plrn_gaussian((/zero, zero, zero/), xkf_all, init_k0_plrn, eigvec)
      CASE (3)
        ALLOCATE(eigvec_wan(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
        IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
        WRITE(stdout, '(5x, a)') "Initializing the polaron wavefunction with previously saved Amp.plrn file"
        CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nstate_plrn, 'Amp.plrn')
        CALL plrn_eigvec_tran('Wan2Bloch', time_rev_A_plrn, eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
           nrr_k, ndegen_k, irvec_r, dims, eigvec)
        DEALLOCATE(eigvec_wan, STAT = ierr)
        IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating eigvec_wan', 1)
      CASE (6)
        WRITE(stdout, '(5x, a, I6)') "Starting from displacements read from file; number of displacements:", init_ntau_plrn
        ALLOCATE(dtau_list(init_ntau_plrn, nktotf, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau_list', 1)
        dtau_list = CMPLX(0.d0, 0.d0)
        !
        IF (init_ntau_plrn == 1) THEN
          filename = 'dtau_disp.plrn'
          CALL read_plrn_dtau(dtau, nqtotf, nmodes, filename, scell_mat_plrn)
          dtau_list(1, :, :) = dtau(:, :)
        ELSE
          DO itau = 1, init_ntau_plrn
            WRITE(tmpch,'(I4)') itau
            filename = TRIM('dtau_disp.plrn_'//ADJUSTL(tmpch))
            CALL read_plrn_dtau(dtau, nqtotf, nmodes, filename, scell_mat_plrn)
            dtau_list(itau, :, :) = dtau(:, :)
          ENDDO
        ENDIF
        !
        CALL mp_bcast(dtau_list, meta_ionode_id, world_comm)
        ! Initialize Ank wavefunction for iterative diagonalization Gaussian
        ! If k0 has not been set on input,     center gaussian at band edge
        IF (ALL(init_k0_plrn(:) == 1000.d0)) init_k0_plrn = xkf_all(1:3, ik_edge)
        WRITE(stdout, '(5x, "Initializing polaron wavefunction using Gaussian wave &
        &packet with the width of", f15.7)') init_sigma_plrn
        WRITE(stdout, '(5x, "centered at k=", 3f15.7)') init_k0_plrn !xkf_all(1:3, ik_edge)
        CALL init_plrn_gaussian((/zero, zero, zero/), xkf_all, init_k0_plrn, eigvec)
        CALL norm_plrn_wf(eigvec, REAL(nktotf, DP))
      CASE DEFAULT
        CALL errore('polaron_scf','init_plrn not implemented!', 1)
    END SELECT
    !
    ! Only keep the coefficients in lowest/highest band,
    ! since the electron/hole localized at this band will be more stable.
    IF (init_plrn <= 2) THEN
      DO ik = 1, nktotf
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          IF (select_bands_plrn(ibnd) /= band_pos) eigvec(indexkn1, 1:nstate_plrn) = czero
        ENDDO
      ENDDO
      CALL norm_plrn_wf(eigvec, REAL(nktotf, DP))
    ENDIF
    CALL stop_clock('init_Ank')
    !
    WRITE(stdout, '(5x, "Starting the SCF cycles")')
    IF (full_diagon_plrn) THEN
      WRITE(stdout, '(5x, a)') "Using serial direct diagonalization"
    ELSE
      WRITE(stdout, '(5x, a)') "Using parallel iterative diagonalization"
      WRITE(stdout, '(5x, "Diagonalizing polaron Hamiltonian with a threshold of ", ES18.6)') ethrdg_plrn
      WRITE(stdout, '(5x, "Please check the results are convergent with this value")')
    ENDIF
    !
    WRITE(stdout, '(/5x, a)') "Starting the self-consistent process"
    WRITE(stdout, '( 5x, a)') REPEAT('-',80)
    WRITE(stdout, '(5x, " iter", 60a15)') "  Eigval/eV", "Phonon/eV", "Electron/eV", &
                                          "Formation/eV", "Error/eV"
    ALLOCATE(Bmat(nqtotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating Bmat', 1)
    !
    IF (scell_mat_plrn) THEN
      CALL read_Rp_in_S()
    ENDIF
    !
    !JLB: possibility of multiple displacements read from file, to calculate polaron energy landscape.
    !     Calculate and print the energies; .plrn files written to disk for last calculation only.
    DO itau = 1, init_ntau_plrn ! ntau_plrn=1 by default
      !
      IF (init_plrn == 6) dtau(:, :) = dtau_list(itau, :, :)
      !
      estmteRt = 1E3
      esterr = 1E5
      DO iter = 1, niter_plrn
        ! Enforce the relation A_k = A*_{G-k} and normalize |A| = 1
        ! Calculating $$ B_{\bq\nu} = \frac{1}{\omega_{\bq,\nu} N_p} \sum_\bk A^\dagger_{\bk+\bq} g_\nu(\bk,\bq) A_\bk $$
        CALL start_clock('cal_bqu')
        IF (init_plrn == 6 .AND. iter==1) THEN
          Bmat = czero
          IF (scell_mat_plrn) THEN
            CALL scell_plrn_bmat_tran('Dtau2Bmat', .true., dtau, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat)
          ELSE
            CALL plrn_bmat_tran('Dtau2Bmat', .true., dtau, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat)
          ENDIF
        ELSE
          CALL build_plrn_bmat(Bmat)
          !
          IF (scell_mat_plrn) THEN
            CALL scell_plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
          ELSE
            CALL plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
          ENDIF
          !
        ENDIF
        !
        dtau_diff = MAXVAL(ABS(REAL(dtau - dtau_save)))
        esterr = dtau_diff
        IF (dtau_diff < conv_thr_plrn .AND. iter > 1) THEN
          !! IF(MAXVAL(ABS(REAL(dtau))) > alat / 2.d0) THEN
          !! KL: dtau_max_plrn is a criteria of maximual polaron displacements
          !! of a physical solution, which is 0.5d0 by default. 
          !! When acoustic phonon modes are found to be important, one can try to increase. 
          IF(MAXVAL(ABS(REAL(dtau))) > alat * dtau_max_plrn) THEN
            CALL errore("polaron_scf","Non-physical solution, check initial guess and convergence.", 1)
          ENDIF
          ! converged, write the final value of eigenvalue
          WRITE(stdout,'(5x,a)') REPEAT('-',80)
          WRITE(stdout, '(5x,a,f10.6,a)' )  'End of self-consistent cycle'
          EXIT
        ELSE
          dtau_save = dtau
        ENDIF
        !
        CALL stop_clock('cal_bqu')
        !
        CALL start_clock('Setup_H')
        !
        CALL build_plrn_hamil(Bmat)
        CALL stop_clock('Setup_H')
        CALL start_clock('DiagonH')
        ! For hole polaron (type_plrn = 1),
        ! we need the highest eigenvalues instead of the lowest eigenvalues
        ! To use KS_solver, which only gives the lowest eigenvalues,
        ! we multiply -1 to the Hamiltonian to get the lowest eigenvalues
        IF (full_diagon_plrn) THEN
          ! Diagonalize Hamiltonian with Serial LAPACK subroutine
          ! Used for testing or robust benchmark
          CALL diag_serial(estmteRt, eigvec)
          !
          CALL mp_bcast(estmteRt, meta_ionode_id, world_comm)
          CALL mp_bcast(eigvec, meta_ionode_id, world_comm)
        ELSE
          ! Diagonalize Hamiltonian with Davidson Solver
          CALL diag_parallel(estmteRt, eigvec)
        ENDIF
        CALL stop_clock('DiagonH')
        !
        ! Reverse the eigenvalues if it is the hole polaron
        estmteRt(1:nstate_plrn) = (-type_plrn) * estmteRt(1:nstate_plrn)

        ! impose the time-reversal symmetry: A^T_k = A_k + A^*_{-k}
        IF (time_rev_A_plrn) CALL check_time_rev_sym(eigvec)
        CALL norm_plrn_wf(eigvec, REAL(nktotf, dp))
        !
        eigvec = (- type_plrn) * eigvec
        !
        CALL start_clock('cal_E_Form')
        CALL calc_form_energy(eplrnphon, eplrnelec)
        CALL stop_clock('cal_E_Form')
        !
        ! TODO : use exact number instead of 20 in 20e15.7
        r_cry(1:3) = IMAG(LOG(berry_phase(1:3) * EXP(- twopi * ci * r0_plrn(1:3)))) / twopi
        r_cry(1:3) = r_cry(1:3) - NINT(r_cry(1:3))
        WRITE(stdout, '(5x, i5, 60e15.4)') iter, estmteRt(1:nstate_plrn) * ryd2ev, eplrnphon * ryd2ev, &
           - eplrnelec * ryd2ev, (eplrnelec + eplrnphon) * ryd2ev, esterr
        eigval = estmteRt
      ENDDO
      !
      ! Calculate and write the energies
      WRITE(stdout, '(5x, a, 50f16.7)') '      Eigenvalue (eV): ', eigval * ryd2ev
      WRITE(stdout, '(5x, a, f16.7)')   '     Phonon part (eV): ', eplrnphon * ryd2ev
      WRITE(stdout, '(5x, a, f16.7)')   '   Electron part (eV): ', eplrnelec * ryd2ev
      IF (init_plrn == 6) THEN
        WRITE(stdout, '(5x, a, f16.7)') 'Formation Energy at this \dtau (eV): ', ((-type_plrn) * eigval - eplrnphon) * ryd2ev
      ELSE
        WRITE(stdout, '(5x, a, f16.7)')   'Formation Energy (eV): ', (eplrnelec + eplrnphon) * ryd2ev
      ENDIF
    ENDDO ! init_ntau_plrn
    !
    ! Calculate and write Density of State of Bqnu and Ank
    WRITE(stdout, '(5x, a)') "Calculating density of states to save in dos.plrn"
    CALL calc_den_of_state(eigvec, Bmat)
    !
    ! Do Bloch to Wannier transform, with U matrix
    CALL start_clock('Ank_trans')
    WRITE(stdout, '(5x, a)') "Generating the polaron wavefunction in Wannier basis to save in Amp.plrn"
    !
    ALLOCATE(eigvec_wan(nbndsub * nktotf, nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
    eigvec_wan = czero
    IF (scell_mat_plrn) THEN
      ! JLB: t_rev set to .true. before
      CALL scell_plrn_eigvec_tran('Bloch2Wan', time_rev_A_plrn, eigvec, nktotf, nRp, Rp, nbndsub, nrr_k, &
              ndegen_k, irvec_r, dims, eigvec_wan)
    ELSE
      ! JLB: t_rev set to .true. before
      CALL plrn_eigvec_tran('Bloch2Wan', time_rev_A_plrn, eigvec, nkf1, nkf2, nkf3, nbndsub, nrr_k, &
         ndegen_k, irvec_r, dims, eigvec_wan)
    ENDIF
    CALL stop_clock('Ank_trans')
    !
    ! Calculate displacements of ions dtau, which is B matrix in Wannier basis
    CALL start_clock('Bqu_tran')
    dtau = czero
    WRITE(stdout, '(5x, a)') "Generating the ionic displacements to save in dtau.plrn and dtau.plrn.xsf"
    IF (scell_mat_plrn) THEN
      CALL scell_plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
    ELSE
      CALL plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
      IF (cal_acous_plrn) THEN
        CALL plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, &
                            dtau_acoustic, acoustic_plrn=acoustic_plrn)
      ENDIF
    ENDIF
    CALL stop_clock('Bqu_tran')
    CALL start_clock('write_files')
    !
    IF (ionode) THEN
      ! Write Amp in Wannier basis
      CALL write_plrn_wf(eigvec_wan, 'Amp.plrn')
      ! Write Ank in Bloch basis
      CALL write_plrn_wf(eigvec, 'Ank.plrn',  etf_all)
      ! Write Bqnu
      CALL write_plrn_bmat(Bmat, 'Bmat.plrn', wf)
      ! Write dtau
      CALL write_plrn_bmat(dtau, 'dtau.plrn')
      IF (cal_acous_plrn) THEN
        CALL write_plrn_bmat(dtau_acoustic, 'dtau.acoustic.plrn')
      ENDIF
      ! Write dtau in a user-friendly format for visulization
      IF (scell_mat_plrn) THEN
        CALL scell_write_plrn_dtau_xsf(dtau, nqtotf, nRp, Rp, as, 'dtau.plrn.xsf')
      ELSE
        CALL write_plrn_dtau_xsf(dtau, nqf1, nqf2, nqf3, 'dtau.plrn.xsf')
      IF (cal_acous_plrn) THEN
        CALL write_plrn_dtau_xsf(dtau_acoustic, nqf1, nqf2, nqf3, 'dtau.acoustic.plrn.xsf')
      ENDIF
      ENDIF
    ENDIF
    CALL stop_clock('write_files')
    !
    DEALLOCATE(dtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating dtau', 1)
    DEALLOCATE(dtau_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating dtau', 1)
    DEALLOCATE(eigvec_wan, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating eigvec_wan', 1)
    DEALLOCATE(Bmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating Bmat', 1)
    IF (cal_acous_plrn) THEN
      DEALLOCATE(dtau_acoustic, STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating dtau_acoustic', 1)
    ENDIF
    !
    CALL stop_clock('main_prln')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE polaron_scf
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE calc_form_energy(eplrnphon, eplrnelec)
    !-----------------------------------------------------------------------
    !!
    !! Computes the polaron formation energy
    !! Note: Require etf_all to be properly initialized
    !!
    !-----------------------------------------------------------------------
    USE ep_constants,    ONLY : zero, czero, twopi, ci, two
    USE modes,           ONLY : nmodes
    USE global_var,      ONLY : xqf, wf, nqtotf, nktotf, nkf
    USE input,           ONLY : type_plrn, nstate_plrn
    USE mp,              ONLY : mp_sum
    USE mp_global,       ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: eplrnphon
    !! Phonon part of polaron formation energy
    REAL(KIND = DP), INTENT(out) :: eplrnelec
    !! Electron part of polaron formation energy
    !
    ! Local variables
    INTEGER :: iq
    !! q point index
    INTEGER :: ik
    !! k point index
    INTEGER :: ik_global
    !! k point global index
    INTEGER :: iq_global
    !! q point global index
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: inu
    !! Phonon mode index
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: iplrn
    !! Polaron state index
    INTEGER :: indexkn1
    !! Combined k-point and band index
    !
    ! Based on Eq. 41 of Ref. 2:
    ! E_{f,ph} = 1/N_p \sum_{q\nu}|B_{q\nu}|^2\hbar\omega_{q\nu}
    ! iq -> q, nqtotf -> N_p
    ! Bmat(iq, inu) -> B_{q\nu}
    ! wf(inu, iq) -> \hbar\omega_{q\nu}
    eplrnphon = zero
    DO iq = 1, nkf
      iq_global = ikqLocal2Global(iq, nqtotf)
      ! JLB - Swapped indices!
      ! I think it would be better to discard modes by looking at the frequencies, i.e. discard negative or zero frequency modes.
      IF(isGVec(xqf(1:3, iq_global))) THEN
        start_mode = 4
      ELSE
        start_mode = 1
      ENDIF
      DO inu = start_mode, nmodes
        eplrnphon = eplrnphon - ABS(Bmat(iq_global, inu))**2 * (wf(inu, iq_global) / nqtotf)
      ENDDO
    ENDDO
    CALL mp_sum(eplrnphon, inter_pool_comm)
    !
    ! E_{f,el} = 1/N_p \sum_{nk}|A_{nk}|^2(\epsilon_{nk}-\epsilon_{F})
    ! indexkn1 -> nk, nktotf -> N_p
    ! eigvec(indexkn1, iplrn) -> A_{nk}
    ! etf_all(select_bands_plrn(ibnd), ik) - ef -> \epsilon_{nk}-\epsilon_{F}
    eplrnelec = zero
    ! TODO: what should we do in iplrn
    DO iplrn = 1, nstate_plrn
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nqtotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik_global - 1) * nbnd_plrn + ibnd
          eplrnelec = eplrnelec - type_plrn * ABS(eigvec(indexkn1, iplrn))**2 / nktotf *&
             etf_all(select_bands_plrn(ibnd), ik_global)
        ENDDO
      ENDDO
    ENDDO
    CALL mp_sum(eplrnelec, inter_pool_comm)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE calc_form_energy
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE find_band_extreme(type_plrn, enk_all, ik_bm, band_loc, efermi)
    !-----------------------------------------------------------------------
    !!
    !! Determine the Fermi energy, read from the input or calculated from band structure
    !!
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, ryd2ev
    USE input,         ONLY : efermi_read, fermi_energy
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : nkf, nktotf
    USE mp,            ONLY : mp_max, mp_min, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: type_plrn
    !! Whether electron polaron (-1) or hole polaron (+1) is to be calculated
    INTEGER, INTENT(out) :: ik_bm
    !! k-point index for CBM or VBM
    INTEGER, INTENT(out) :: band_loc
    !! band index for CBM or VBM
    REAL(KIND = DP), INTENT(in) :: enk_all(:,:)
    !! KS eigenvalues
    REAL(KIND = DP), INTENT(out) :: efermi
    !! Fermi energy
    !
    ! Local variable
    INTEGER :: ik
    !! k-point index
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: k_extreme_local(npool)
    !! Auxiliary index of CBM or VBM k-point at different pools
    INTEGER :: ipool(1)
    !! Index to locate CBM or VBM
    REAL(KIND = DP) :: band_edge
    !! Temporary KS eigenvalue
    REAL(KIND = DP) :: extreme_local(npool)
    !! Auxiliary CBM or VBM energy at different pools
    !
    IF (type_plrn == -1 ) THEN
      band_loc = select_bands_plrn(1)
    ELSE IF ( type_plrn == 1 ) THEN
      band_loc = select_bands_plrn(nbnd_plrn)
    ENDIF
    !
    WRITE(stdout, '(5x, "The band extremes are at band ",  i0)') band_loc
    !
    ! Determine the Fermi energy, read from the input or calculated from band structure
    ! = 1E4*(-type_plrn)
    ik_bm = 0
    k_extreme_local = 0
    extreme_local = zero
    IF (efermi_read) THEN
      efermi = fermi_energy
      WRITE(stdout, '(5x, "Polaron Reference energy (VBM or CBM) is read from the input file: ",&
         &f16.6, " eV.")') efermi * ryd2ev
    ELSE
      IF(type_plrn == 1) THEN
        efermi = -1E5
      ELSE IF (type_plrn == -1) THEN
        efermi =  1E5
      ELSE
        CALL errore('','Wrong type_plrn, should be 1 or -1', 1)
      ENDIF
      !
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        band_edge = enk_all(band_loc, ik_global)
        !
        IF (type_plrn == 1) THEN
          ! For hole polaron (type_plrn = 1), find the highest eigenvalue
          IF (band_edge > efermi) THEN
            efermi = band_edge
            ik_bm = ik_global
          ENDIF
        ELSE IF (type_plrn == -1) THEN
          ! For electron polaron (type_plrn = -1), find the lowest eigenvalue
          IF (band_edge < efermi) THEN
            efermi = band_edge
            ik_bm = ik_global
          ENDIF
        ENDIF
      ENDDO
      !
      k_extreme_local(my_pool_id + 1) = ik_bm
      extreme_local(my_pool_id + 1) = efermi
      CALL mp_sum(k_extreme_local, inter_pool_comm)
      CALL mp_sum(extreme_local, inter_pool_comm)
      !
      IF (type_plrn == 1) THEN
        ipool = MAXLOC(extreme_local)
        ik_bm = k_extreme_local(ipool(1))
        efermi = MAXVAL(extreme_local)
      ELSE IF (type_plrn == -1) THEN
        ipool = MINLOC(extreme_local)
        ik_bm = k_extreme_local(ipool(1))
        efermi = MINVAL(extreme_local)
      ENDIF
    ENDIF
    !-----------------------------------------------------------------------
    END SUBROUTINE find_band_extreme
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE gather_band_eigenvalues(etf, enk_all)
    !-----------------------------------------------------------------------
    !! Gather all the eigenvalues to determine the EBM/VBM,
    !! and calculate the density state of Ank and Bqnu
    !-----------------------------------------------------------------------
    USE input,         ONLY : nbndsub
    USE global_var,    ONLY : nkqf, nktotf
    USE ep_constants,  ONLY : zero
    USE parallelism,  ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: etf(:, :)
    !! Eigenvalues per cpu
    REAL(KIND = DP), INTENT(out) :: enk_all(:, :)
    !! Eigenvalues (total)
    !
    ! Local variables
    INTEGER :: ierr
    !! Error index
    REAL(KIND = DP), ALLOCATABLE :: rtmp2(:, :)
    !! Temporary variable to gather eigenvalues across pools
    !
    ALLOCATE(rtmp2(nbndsub, nktotf*2), STAT = ierr)
    IF (ierr /= 0) CALL errore('gather_band_eigenvalues', 'Error allocating rtmp2', 1)
    rtmp2 = zero
    !
    CALL poolgather2 ( nbndsub, nktotf*2, nkqf, etf, rtmp2  )
    enk_all(1:nbndsub, 1:nktotf) = rtmp2(1:nbndsub, 1:nktotf*2:2)
    !
    DEALLOCATE(rtmp2, STAT = ierr)
    IF (ierr /= 0) CALL errore('gather_band_eigenvalues', 'Error deallocating rtmp2', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gather_band_eigenvalues
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE cal_phonon_eigenfreq(nrr_q, irvec_q, ndegen_q, rws, nrws, wfreq)
    !-----------------------------------------------------------------------
    !! Calculate the phonon eigen frequencies. This is needed when restarting the polaron
    !! calculation with recalculating el-ph vertex
    !-----------------------------------------------------------------------
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : xqf, nqtotf
    USE ep_constants,  ONLY : zero, eps8, czero
    USE wannier2bloch, ONLY : dynwan2bloch, dynifc2blochf
    USE input,         ONLY : lifc
    USE io_global,     ONLY : ionode, stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    REAL(KIND = DP), INTENT(out) :: wfreq(:, :)
    !! Phonon frequencies
    !
    ! Local variables
    INTEGER  :: inu
    !! Phonon mode counter
    INTEGER  :: iq
    !! q-point counter
    REAL(KIND = DP) :: w2(nmodes)
    !! Phonon frequency squared
    REAL(KIND = DP) :: xxq(3)
    !! q-point coordinates
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! Phonon eigenvectors
    !
    uf = czero
    w2 = zero
    !
    !TODO: make this part parallel over q
    wfreq = zero
    DO iq = 1, nqtotf
      ! iq -> q
      xxq = xqf(1:3, iq)
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf, w2, is_mirror_q(iq))
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq, uf, w2, is_mirror_q(iq))
      ENDIF
      DO inu = 1, nmodes
        IF (w2(inu) > -eps8) THEN
          wfreq(inu, iq) =  DSQRT(ABS(w2(inu)))
        ELSE
          IF (ionode) THEN
            WRITE(stdout, '(5x, "WARNING: Imaginary frequency mode ",&
            &I6, " at iq=", I6)') inu, iq
          ENDIF
          wfreq(inu, iq) = 0.d0
        ENDIF
      ENDDO ! inu
    ENDDO ! iq
    !-----------------------------------------------------------------------
    END SUBROUTINE cal_phonon_eigenfreq
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE init_plrn_random(eigvec_init)
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nktotf
    USE ep_constants,  ONLY : ci, cone
    USE input,         ONLY : nstate_plrn
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_init(:, :)
    !! Polaron wf coefficients, Ank, upon initialization
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP), ALLOCATABLE :: rmat_tmp(:, :)
    !! Temporary variable for random number matrix
    !
    CALL RANDOM_SEED()
    ALLOCATE(rmat_tmp(1:nktotf*nbnd_plrn, 1:nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('init_plrn_random', 'Error allocating rmat_tmp', 1)
    CALL RANDOM_NUMBER(rmat_tmp)
    eigvec_init(1:nktotf * nbnd_plrn, 1:nstate_plrn) = cone * rmat_tmp(1:nktotf * nbnd_plrn, 1:nstate_plrn)
    CALL RANDOM_NUMBER(rmat_tmp)
    eigvec_init(1:nktotf * nbnd_plrn, 1:nstate_plrn) = eigvec_init(1:nktotf*nbnd_plrn, 1:nstate_plrn) + &
                                                ci * rmat_tmp(1:nktotf * nbnd_plrn, 1:nstate_plrn)
    DEALLOCATE(rmat_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('init_plrn_random', 'Error deallocating rmat_tmp', 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE init_plrn_random
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE init_plrn_gaussian(r0, k_all, k0, eigvec_init)
    !-----------------------------------------------------------------------
    !! Initialize Ank coefficients with a Gaussian lineshape
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, cone, ci, twopi, one, zero
    USE input,         ONLY : init_sigma_plrn
    USE global_var,    ONLY : nktotf
    USE cell_base,     ONLY : bg, alat
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: r0(3)
    !! Center of Gaussian in real space
    REAL(KIND = DP), INTENT(in) :: k0(3)
    !! Center of Gaussian in k-space
    REAL(KIND = DP), INTENT(in) :: k_all(:, :)
    !! List with all k-point coordinates
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_init(:, :)
    !! Polaron wf coefficients, Ank, upon initialization
    !
    ! Local variables
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: indexkn1
    !! Combined k-point and band counter
    INTEGER :: ishift
    !! Shift counter
    REAL(KIND = DP) :: qcart(3)
    !! q-point coordinates in cartesian
    REAL(KIND = DP) :: xxq(3)
    !! q-point coordinates
    REAL(KIND = DP) :: shift(3)
    !! shift coordinates
    REAL(KIND = DP) :: disK
    !! Value of Gaussian in neighbor BZs
    COMPLEX(KIND = DP) :: ctemp
    !! Exponential prefactor from Fourier transform of Gaussian
    !
    DO ik = 1, nktotf
      xxq = k_all(1:3, ik) - (k0(:) - INT(k0(:))) ! shift k0 to 1BZ
      CALL dgemv('n', 3, 3, one, bg, 3, xxq, 1, zero, qcart, 1)
      ctemp = EXP(-ci * twopi * DOT_PRODUCT( qcart, r0 ))
      disK = -1
      ! Ensure periodicity of Ank checking distance to other equivalent BZ points
      DO ishift = 1, 27
        shift(1:3) = REAL(index_shift(ishift), KIND=DP)
        CALL dgemv('n', 3, 3, one, bg, 3, xxq+shift, 1, zero, qcart, 1)
        disK = MAX(disK, EXP(-init_sigma_plrn * NORM2(qcart) * twopi / alat )) ! for sigma to be in bohr
      ENDDO
      !
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        eigvec_init(indexkn1, :) = CONE * disK * ctemp
      ENDDO
    ENDDO
    !-----------------------------------------------------------------------
    END SUBROUTINE init_plrn_gaussian
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE build_plrn_bmat(bqv)
    !-----------------------------------------------------------------------
    !! Create the Bmat
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nkf, nqtotf, xqf, nktotf
    USE input,         ONLY : model_vertex_plrn, io_lvl_plrn,               &
                              g_start_energy_plrn, g_end_energy_plrn,       &
                              g_start_band_plrn, model_vertex_plrn
    USE ep_constants,  ONLY : czero, one, two, zero, cone, eps2, eps8
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(out) :: bqv(:, :)
    !! Polaron displacement coefficients in phonon basis, Bqv
    !
    ! Local variables
    INTEGER :: iq
    !! q-point counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: iqpg
    !! Mirror q-point index
    REAL(KIND = DP) :: eig
    !! KS eigenvalue
    !
    bqv = czero
    DO iq = 1, nqtotf
      IF (model_vertex_plrn) THEN
        epf = czero
        epf(1, 1, nmodes, 1:nkf) = gq_model(iq)
      ELSE
        IF (io_lvl_plrn == 0) THEN
          epf(:, :, :, :) = epfall(:, :, :, :, iq)
        ELSE IF (io_lvl_plrn == 1) THEN
          CALL get_buffer(epf, lword_g, iepfall, iq)
        ENDIF
      ENDIF
      IF (test_tags_plrn(1)) THEN
        epf = 1E-6
        epf(:, :, nmodes, :) = 2E-2
      ELSE IF (test_tags_plrn(2)) THEN
        epf = ABS(epf)
      ELSE IF (test_tags_plrn(3)) THEN
        IF (NORM2(xqf(:,iq)) > 1E-5) epf(:, :, nmodes, :) = 0.005 / NORM2(xqf(:, iq))
      ENDIF
      !
      ! energy cutoff for g
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_g_plrn
          eig = etf_all(ibnd + g_start_band_plrn - 1, ik_global)
          IF (eig < g_start_energy_plrn .OR. eig > g_end_energy_plrn) THEN
            epf(ibnd, :, :, ik) = czero
            epf(:, ibnd, :, ik) = czero
          ENDIF
        ENDDO
      ENDDO
      !
      iqpg = kpg_map(iq)
      ! if iq is the gamma point, the first three modes should be
      ! dropped because wf(q=0, 1:3) will be zero
      IF (isGVec(xqf(1:3, iq))) THEN
        start_mode = 4
      ELSE
        start_mode = 1
      ENDIF
      !
      IF ( iq > iqpg ) THEN
        DO inu = start_mode, nmodes!
          ! Enforce the relation B_q = B*_{G-q}
          bqv(iq, inu) = CONJG(bqv(iqpg, inu))
        ENDDO
      ELSE
        DO inu = start_mode, nmodes
          bqv(iq, inu)   = cal_Bmat(iq, inu)
        ENDDO
      ENDIF
      !TODO: to be consistent with Denny, whether this is correct?
      ! IF ( ABS(wf(1, iq)) < eps2 ) bqv(iq, 1) = czero
    ENDDO
    ! cal_bqv only sum over local k, so we have to do mp_sum
    CALL mp_sum(bqv, inter_pool_comm )
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE build_plrn_bmat
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE build_plrn_hamil(bqv)
    !-----------------------------------------------------------------------
    !! Build the effective polaron Hamiltonian,
    !! Eq.(61) of PRB 99, 235139 (2019)
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nkf, nqtotf, xqf, nktotf
    USE input,         ONLY : model_vertex_plrn, io_lvl_plrn,  type_plrn,           &
                              nhblock_plrn, g_start_energy_plrn, g_end_energy_plrn, &
                              g_start_band_plrn
    USE ep_constants,  ONLY : czero, one, two, zero, cone, eps2, eps8, twopi, ci
    USE mp,            ONLY : mp_sum
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(in) :: bqv(:,:)
    !! FIXME
    !
    ! Local variables
    INTEGER :: iq
    !! q-point counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ikq
    !! k+q point counter
    INTEGER :: ik_global
    !! Globar k-point index
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: jbnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: index_blk
    !! FIXME
    INTEGER :: index_loc
    !! FIXME
    REAL(KIND = DP) :: eig
    !! KS eigenvalue
    COMPLEX(KIND = DP) :: ctemp
    !! Prefactor
    !
    test_tags_plrn(1) = .FALSE.
    test_tags_plrn(2) = .FALSE.
    test_tags_plrn(3) = .FALSE.
    !
    ! Calculate the Hamiltonian with Bq $$H_{n\bk,n'\bk'} = \delta_{n\bk,n'\bk'}\varepsilon_{n\bk} -\frac{2}{N_p} \sum_{\nu} B^*_{\bq,\nu}g_{nn'\nu}(\bk',\bq)$$
    ! H_{n\bk,n'\bk'} -> Hamil(ik, ibnd, ikq, jbnd)
    ! B^*_{\bq,\nu} -> conj(bqv(iq, inu))
    ! g_{nn'\nu}(\bk',\bq) -> epf(ibnd, jbnd, inu, ikq, iq)
    ! if q == 0, \delta_{n\bk,n'\bk'}\varepsilon_{n\bk} is diagonal matrix with \varepsilon_{n\bk}
    !
    ! G == (0,0,0) means this is the diagonal term with k=k'
    ! ikq is the global index, i.e. the second index
    ! ik is the local index, i.e. the first index
    Hamil = czero
    DO iq = 1, nqtotf
      IF (model_vertex_plrn) THEN
        epf = czero
        epf(1, 1, nmodes, 1:nkf) = gq_model(iq)
      ELSE
        CALL start_clock('read_gmat') !nbndsub*nbndsub*nmodes*nkf
        IF(io_lvl_plrn == 0) THEN
          epf(:, :, :, :) = epfall(:, :, :, :, iq)
        ELSE IF (io_lvl_plrn == 1) THEN
          CALL get_buffer(epf, lword_g, iepfall, iq)
        ENDIF
        CALL stop_clock('read_gmat')
      ENDIF
      ! energy cutoff of g
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_g_plrn
          eig = etf_all(ibnd + g_start_band_plrn - 1, ik_global)
          IF (eig < g_start_energy_plrn .OR. eig > g_end_energy_plrn) THEN
            epf(ibnd, :, :, ik) = czero
            epf(:, ibnd, :, ik) = czero
          ENDIF
        ENDDO
      ENDDO
      !
      DO ik = 1, nkf
        ikq = ikq_all(ik, iq)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          !
          IF (nhblock_plrn == 1) THEN
            index_loc = indexkn1
            index_blk = 1
          ELSE
            index_loc = MOD(indexkn1 - 1, hblocksize) + 1
            index_blk = INT((indexkn1 - 1) / hblocksize) + 1
          ENDIF
          !
          IF (iq /= 1 .AND. index_loc == 1 .AND. nhblock_plrn /= 1) THEN
            CALL start_clock('read_Hmat')
            CALL get_buffer(Hamil, lword_h, ihamil, index_blk)
            CALL stop_clock('read_Hmat')
          ENDIF
          !
          IF (isGVec(xqf(1:3, iq))) THEN
            ! Note that, ik is local index while ikq is global index,
            ! so even when q=0, ik \= ikq, but ik_global == ikq
            ! delta_{nn' kk'} epsilon_{nk}
            ctemp = etf_all(select_bands_plrn(ibnd), ikq)
            indexkn2 = (ikq - 1) * nbnd_plrn + ibnd
            Hamil(indexkn2, index_loc) = Hamil(indexkn2, index_loc) + ctemp
          ENDIF
          !
          DO jbnd = 1, nbnd_plrn
            indexkn2 = (ikq - 1) * nbnd_plrn + jbnd
            DO inu = 1, nmodes
              ctemp = type_plrn * two / REAL(nqtotf, KIND = DP) * (bqv(iq, inu)) * &
                 CONJG(epf(select_bands_plrn(jbnd) - g_start_band_plrn + 1, &
                       select_bands_plrn(ibnd)- g_start_band_plrn + 1, inu, ik))
              Hamil(indexkn2, index_loc) = Hamil(indexkn2, index_loc) + ctemp
            ENDDO
          ENDDO
          IF (nhblock_plrn /= 1) THEN
            IF ((index_loc == hblocksize .OR. indexkn1 == nkf * nbnd_plrn)) THEN
              CALL start_clock('Write_Hmat')
              CALL save_buffer(Hamil, lword_h, ihamil, index_blk)
              Hamil = czero
              CALL stop_clock('Write_Hmat')
            ENDIF
          ENDIF
        ENDDO !ibnd
      ENDDO !ik
    ENDDO ! iq
    !-----------------------------------------------------------------------
    END SUBROUTINE build_plrn_hamil
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION isGVec(xxk)
    !-----------------------------------------------------------------------
    !! Return true if xxk integer times of the reciprocal vector
    !! if xxk is the difference of two vectors, then return true if these
    !! two vector are the same
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : eps6
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xxk(3)
    !! k-point coordinate
    !
    ! Local variable
    LOGICAL :: isGVec
    !! .true. if k-point is a G-vector
    !
    isGVec = &
      ABS(xxk(1) - NINT(xxk(1))) < eps6 .AND. &
      ABS(xxk(2) - NINT(xxk(2))) < eps6 .AND. &
      ABS(xxk(3) - NINT(xxk(3))) < eps6
    !-----------------------------------------------------------------------
    END FUNCTION isGVec
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION ikqLocal2Global(ikq, nkqtotf)
    !-----------------------------------------------------------------------
    !! Return the global index of the local k point ik
    !-----------------------------------------------------------------------
    USE parallelism, ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ikq
    !! k+q point counter
    INTEGER, INTENT(in) :: nkqtotf
    !! number of k-points in the fine grid
    !
    ! Local variable
    INTEGER :: ikqLocal2Global
    !! Index of k+q point in global list
    INTEGER :: startn
    !! Lower bound for k-points in pools
    INTEGER :: lastn
    !! Upper bound for k-points in pools
    !
    CALL fkbounds(nkqtotf, startn, lastn)
    !
    ikqLocal2Global = startn + ikq - 1
    IF (ikqLocal2Global > lastn) THEN
      CALL errore('ikqLocal2Global', 'Index of k/q is beyond this pool.', 1)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END FUNCTION ikqLocal2Global
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION ikGlobal2Local(ik_g, nktotf)
    !-----------------------------------------------------------------------
    !! Return the global index of the local k point ik
    !-----------------------------------------------------------------------
    USE parallelism, ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik_g
    !! Global k-point index
    INTEGER, INTENT(in) :: nktotf
    !! Number of k-points in global fine grid
    !
    ! Local variable
    INTEGER :: ikGlobal2Local
    !! Index of k-point in local pool list
    INTEGER :: startn
    !! Lower bound for k-points in pools
    INTEGER :: lastn
    !! Upper bound for k-points in pools
    !
    CALL fkbounds(nktotf, startn, lastn)
    !
    ikGlobal2Local = ik_g - startn + 1
    !
    IF(ikGlobal2Local <= 0) THEN
      ikGlobal2Local = 0
    ENDIF
    !-----------------------------------------------------------------------
    END FUNCTION ikGlobal2Local
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE cal_f_delta(energy, sigma, f_delta)
    !-----------------------------------------------------------------------
    !! Return EXP(-(energy/sigma)**2)
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : ryd2mev, one, ryd2ev, two, zero
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: energy(:)
    !! Energy variable
    REAL(KIND = DP), INTENT(in) :: sigma
    !! WIdth of Gaussian
    REAL(KIND = DP), INTENT(out) :: f_delta(:)
    !! Gaussian function
    !
    f_delta = EXP(-(energy / sigma)**2)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE cal_f_delta
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE h_psi_plrn(lda, n, m, psi, hpsi)
    !-----------------------------------------------------------------------
    ! Calculate Hpsi with psi as input to use the diagon solver in KS_solver cegterg
    ! cegterg take two external subroutine to calculate Hpsi and Spsi to calculate
    ! ( H - e S ) * evc = 0, since H and S is not saved due to their sizes
    ! Hamil need to be passed to h_psi because the parameter space is fixed
    ! to meet the requirement of Davidson diagonalization.
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nkf, nktotf
    USE input,         ONLY : type_plrn, nhblock_plrn
    USE ep_constants,  ONLY : czero, one, two, zero, cone, eps2, ci
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE ep_constants,  ONLY : czero
    USE global_var,    ONLY : nkf, nktotf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: lda
    !! leading dimension of arrays psi, spsi, hpsi, which is nkf * nbnd_plrn
    INTEGER, INTENT(in) :: n
    !! true dimension of psi, spsi, hpsi
    INTEGER, INTENT(in) :: m
    !! number of states psi
    COMPLEX(KIND = DP), INTENT(inout) :: psi(lda, m)
    !! the wavefunction
    COMPLEX(KIND = DP), INTENT(out) :: hpsi(lda, m)
    !! FIXME
    !
    ! Local variables
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ikq
    !! k+q point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: jbnd
    !! Electron band index
    INTEGER :: indexkn1
    !! Combined electron and k-point index
    INTEGER :: indexkn2
    !! Combined electron and k-point index
    INTEGER :: index_loc
    !! FIXME
    INTEGER :: index_blk
    !! FIXME
    !
    ! Gather psi (dimension nkf) to form eigvec (dimension nktotf)
    CALL start_clock('cal_hpsi')
    IF (lda < nkf * nbnd_plrn) CALL errore('h_psi_plrn', 'leading dimension of arrays psi is not correct', 1)
    eigvec = czero
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
        eigvec(indexkn2, 1:m) = psi(indexkn1, 1:m)
      ENDDO
    ENDDO
    CALL mp_sum(eigvec, inter_pool_comm)
    !
    ! Iterative diagonalization only get the lowest eigenvalues,
    ! however, we will need the highest eigenvalues if we are calculating hole polaron
    ! so, we multiply hpsi by -1, and get eigenvalues in diagonalization
    ! and then multiply eigenvalues by -1.
    hpsi(1:lda, 1:m) = czero
    DO ik = 1, nkf
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        IF(nhblock_plrn == 1) THEN
          index_loc = indexkn1
          index_blk = 1
        ELSE
          index_loc = MOD(indexkn1 - 1, hblocksize) + 1
          index_blk = INT((indexkn1 - 1) / hblocksize) + 1
        ENDIF
        IF (index_loc == 1 .AND. nhblock_plrn /= 1) CALL get_buffer(Hamil, lword_h, ihamil, index_blk)
        DO ikq = 1, nktotf
          DO jbnd = 1, nbnd_plrn
            indexkn2 = (ikq - 1) * nbnd_plrn + jbnd
            hpsi(indexkn1, 1:m) = hpsi(indexkn1, 1:m) - &
               type_plrn * Hamil(indexkn2, index_loc) * eigvec(indexkn2, 1:m)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL stop_clock('cal_hpsi')
    !-----------------------------------------------------------------------
    END SUBROUTINE h_psi_plrn
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE s_psi_plrn(lda, n, m, psi, spsi)
    !-----------------------------------------------------------------------
    !! FIXME ???
    !! SP - what is the point of this subroutine ? Can it be removed ?
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: lda
    !! FIXME
    INTEGER, INTENT(in) :: n
    !! FIXME
    INTEGER, INTENT(in) :: m
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: psi(lda,m)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: spsi(lda,m)
    !! FIXME
    CALL errore('s_psi_plrn', "WARNING: This function should not be called at all!", 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE s_psi_plrn
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE g_psi_plrn(lda, n, m, npol, psi, e)
    !-----------------------------------------------------------------------
    !! This routine computes an estimate of the inverse Hamiltonian
    !! and applies it to m wavefunctions.
    !! SP: This routines does nothing !!
    !! FIXME
    !
    IMPLICIT NONE
    !
    INTEGER     :: lda, n, m, npol
    COMPLEX(KIND = DP) :: psi(lda, npol, m)
    REAL(KIND = DP)    :: e(m)
    !-----------------------------------------------------------------------
    END SUBROUTINE g_psi_plrn
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE get_cfac(xk, nrr_k, irvec_r, cfac)
    !-----------------------------------------------------------------------
    !! Compute the exponential factor.
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : twopi, ci, czero
    USE kinds,         ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    REAL(KIND = DP), INTENT(in) :: xk(3) 
    !! k-point coordinates
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(out) :: cfac(nrr_k)
    !! Exponential prefactor 
    !
    ! Local Variables
    REAL(KIND = DP) :: rdotk(nrr_k)
    !! Dot product between k-point and R WS vector
    !
    cfac = czero
    rdotk = czero
    !
    CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xk, 1, 0.0_dp, rdotk, 1 )
    cfac(:) = EXP(ci * rdotk(:))
    !-----------------------------------------------------------------------
    END SUBROUTINE
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    FUNCTION cal_Bmat(iq, inu)
    !-----------------------------------------------------------------------
    !!
    !! This function calculates the Bq matrix:
    !! B_{qu} = 1/N_p \sum_{mnk} A^*_{mk+q}A_{nk} [g_{mnu}(k,q)/\hbar\omega_{qu}]
    !! Eq.(38) of PRB 99, 235139 (2019)
    !!
    !-----------------------------------------------------------------------
    USE global_var,    ONLY : nkf, nktotf, wf, nqtotf
    USE input,         ONLY : eps_acoustic, &
                              g_start_band_plrn
    USE ep_constants,  ONLY : czero, one, eps2, cone, eps8
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! q-point counter
    INTEGER, INTENT(in) :: inu
    !! Phonon mode counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ikq
    !! k+q point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: jbnd
    !! Electron-band index
    INTEGER :: iplrn
    !! Polaron state index
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    COMPLEX(KIND = DP) :: cal_Bmat
    !! Polaron displacement coefficients in phonon basis, Bqv
    COMPLEX(KIND = DP) :: prefac
    !! Prefactor variable
    !
    ! sum k = sum 1 to nkf + mp_sum (inter_pool)
    ! mp_sum is in polaron_scf
    cal_Bmat = czero
    DO ik = 1, nkf
      ikq = ikq_all(ik, iq)
      ik_global = ikqLocal2Global(ik, nktotf)
      ! TODO : what should do for iplrn?
      DO iplrn = 1, 1
        DO ibnd = 1, nbnd_plrn
          DO jbnd = 1, nbnd_plrn
            indexkn1 = (ikq - 1) * nbnd_plrn + ibnd
            indexkn2 = (ik_global - 1) * nbnd_plrn + jbnd
            IF (wf(inu, iq) > eps_acoustic ) THEN
              prefac = cone / (wf(inu, iq) * REAL(nqtotf, DP))
            ELSE
              prefac = czero
            ENDIF
            ! B_{q\nu} = \frac{1}{N_p}\sum_{nn'k}A^*_{n'k+q}\frac{g_{n'n\nu}(k, q)}{\hbar \omega_{q\nu}} A_{nk}
            cal_Bmat = cal_Bmat + prefac * (eigvec(indexkn2, iplrn)) * CONJG(eigvec(indexkn1, iplrn)) * &
               (epf(select_bands_plrn(ibnd) - g_start_band_plrn + 1, &
               select_bands_plrn(jbnd) - g_start_band_plrn + 1, inu, ik)) !conjg
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! JLB - discard zero or imaginary frequency modes
    IF (wf(inu, iq) < eps_acoustic) THEN
      cal_Bmat = czero
    ENDIF
    !-----------------------------------------------------------------------
    END FUNCTION cal_Bmat
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION indexGamma(k_all)
    !-----------------------------------------------------------------------
    !! Find the index of Gamma point i.e. (0, 0, 0) in xkf_all
    !! which contains all the crystal coordinates of the k/q points
    !! if Gamma point is not included, return 0
    !
    !-----------------------------------------------------------------------
    USE global_var,   ONLY : nkf, nktotf
    USE mp,           ONLY : mp_sum
    USE mp_global,    ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: k_all(:, :)
    !! crystal coordinates of k/q points. 
    !! Renamed from xkf_all to avoid variable shadowing.
    !
    ! Local variable
    INTEGER :: indexGamma
    !! Index of \Gamma point in global k-point list
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    !
    indexGamma = 0
    !
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      IF(isGVec(k_all(1:3, ik_global))) THEN
        indexGamma = ik_global
      ENDIF
    ENDDO
    CALL mp_sum(indexGamma, inter_pool_comm)
    !
    IF (.NOT. isGVec(k_all(1:3, indexGamma))) THEN
      CALL errore('indexGamma','The index of Gamma point is wrong!', 1)
    ENDIF
    !----------------------------------------------------------------------
    END FUNCTION indexGamma
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    SUBROUTINE norm_plrn_wf(eigvec_coef, norm_new)
    !----------------------------------------------------------------------
    !! Computes the norm of the polaron wavefunction
    !----------------------------------------------------------------------
    USE global_var,    ONLY : nktotf
    USE input,         ONLY : nstate_plrn
    USE ep_constants,  ONLY : czero, one, two, cone
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: norm_new
    !! Wave function normalization
    COMPLEX(KIND = DP), INTENT(inout) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in Bloch basis, Ank
    !
    ! Local variable
    INTEGER :: iplrn
    !! Polaron state counter
    REAL(KIND = DP) :: norm
    !! Normalization
    !
    DO iplrn = 1, nstate_plrn
      norm = REAL(DOT_PRODUCT(eigvec_coef(1:nbnd_plrn * nktotf, iplrn), eigvec_coef(1:nbnd_plrn * nktotf, iplrn)))
      eigvec_coef(:, iplrn) = eigvec_coef(:, iplrn) / DSQRT(norm) * SQRT(norm_new)
    ENDDO
    !-----------------------------------------------------------------------
    END SUBROUTINE norm_plrn_wf
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE check_time_rev_sym(eigvec_coef)
    !-----------------------------------------------------------------------
    !! Enforces TR symmetry on polaron wave function coefficients.
    !! Not used by default, only for testing purposes.
    !-----------------------------------------------------------------------
    USE global_var,     ONLY : nkf, nktotf
    USE ep_constants,   ONLY : czero, one, two, cone
    USE mp,             ONLY : mp_sum
    USE mp_global,      ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(inout) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in Bloch basis, Ank
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: ikpg
    !! Index of mirror k-point
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    INTEGER :: nPlrn_l
    !! Number of polaron states
    REAL(KIND = DP) :: norm
    !! Norm of polaron wave function
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_save(:, :)
    !! Auxiliary array for polaron wave function coefficients
    !
    ! nstate_plrn
    nPlrn_l = 1
    !
    ALLOCATE(eigvec_save(nktotf * nbnd_plrn, nPlrn_l), STAT = ierr)
    IF (ierr /= 0) CALL errore('check_time_rev_sym', 'Error allocating eigvec_save', 1)
    eigvec_save = czero
    !
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      ikpg = kpg_map(ik_global)
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ikpg - 1) * nbnd_plrn + ibnd
        indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
        eigvec_save(indexkn1, 1:nPlrn_l)  = CONJG(eigvec(indexkn2, 1:nPlrn_l))
      ENDDO
    ENDDO
    CALL mp_sum(eigvec_save, inter_pool_comm)
    eigvec_coef(:, 1:nPlrn_l) = (eigvec_coef(:, 1:nPlrn_l) + eigvec_save(:, 1:nPlrn_l))
    !
    DO iplrn = 1, nPlrn_l
      norm = REAL(DOT_PRODUCT(eigvec_coef(1:nbnd_plrn*nktotf, iplrn), eigvec_coef(1:nbnd_plrn*nktotf, iplrn)))!nktotf*nbnd_plrn*
      eigvec_coef(:, iplrn) = eigvec_coef(:, iplrn)/DSQRT(norm)
    ENDDO
    !
    DEALLOCATE(eigvec_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('check_time_rev_sym', 'Error deallocating eigvec_save', 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE check_time_rev_sym
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE diag_serial(estmteRt, eigvec_coef)
    !-----------------------------------------------------------------------
    !! Serial diagonalization using LAPACK library
    !-----------------------------------------------------------------------
    USE ep_constants,        ONLY : czero, twopi, ci, cone, zero
    USE global_var,          ONLY : nkf, nktotf
    USE input,               ONLY : nstate_plrn, &
                                    type_plrn, nhblock_plrn
    USE io_global,           ONLY : ionode, meta_ionode_id
    USE mp_world,            ONLY : world_comm
    USE mp_global,           ONLY : inter_pool_comm
    USE mp,                  ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: estmteRt(:)
    !! Polaron eigenvalue
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_coef(:, :)
    !! Polaron eigenvector coefficients
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    INTEGER :: lwork
    !! FIXME
    INTEGER :: info
    !! FIXME
    INTEGER :: mm
    !! FIXME
    INTEGER :: index_loc
    !! FIXME
    INTEGER :: index_blk
    !! FIXME
    INTEGER, ALLOCATABLE :: iwork(:)
    !! FIXME
    INTEGER, ALLOCATABLE :: ifail(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !! FIXME
    COMPLEX(KIND = DP),  ALLOCATABLE :: work(:)
    !! FIXME
    COMPLEX(KIND = DP),  ALLOCATABLE :: Hamil_save(:,:)
    !! Auxiliary array for storing Hamiltonian for all k-points
    COMPLEX(KIND = DP),  ALLOCATABLE :: Identity(:,:)
    !! FIXME
    !
    ALLOCATE(Hamil_save(nktotf * nbnd_plrn, nktotf * nbnd_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating Hamil_save', 1)
    !
    Hamil_save = czero
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        !
        index_loc = MOD(indexkn1 - 1, hblocksize) + 1
        index_blk = INT((indexkn1 - 1) / hblocksize) + 1
        IF (index_loc == 1 .AND. nhblock_plrn /= 1) CALL get_buffer(Hamil, lword_h, ihamil, index_blk)
        !
        indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
        Hamil_save(indexkn2, 1:nktotf * nbnd_plrn) = - type_plrn * Hamil(1:nktotf * nbnd_plrn, index_loc)
      ENDDO
    ENDDO
    CALL mp_sum(Hamil_save, inter_pool_comm)
    IF (ionode) THEN
      ALLOCATE(Identity(nktotf * nbnd_plrn, nktotf * nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating Identity', 1)
      lwork = 5 * nktotf * nbnd_plrn
      ALLOCATE(rwork(7 * nktotf * nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating rwork', 1)
      ALLOCATE(iwork(5 * nktotf * nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating iwork', 1)
      ALLOCATE(ifail(nktotf * nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating ifail', 1)
      ALLOCATE(work(lwork), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error allocating work', 1)
      Identity = czero
      !
      DO ibnd = 1, nbnd_plrn * nktotf
        Identity(ibnd, ibnd) = cone
      ENDDO
      !
      eigvec_coef = czero
      estmteRt = zero
      !
      CALL ZHEGVX( 1, 'V', 'I', 'U', nktotf * nbnd_plrn, Hamil_save, nktotf * nbnd_plrn, Identity,&
         nktotf * nbnd_plrn, zero, zero, 1, nstate_plrn, zero, mm, estmteRt(1:nstate_plrn), &
         eigvec_coef, nktotf * nbnd_plrn, work, lwork, rwork, iwork, ifail, info)
      !
      IF (info /= 0) CALL errore('diag_serial','Polaron: diagonal error.', 1)
      DEALLOCATE(rwork, iwork, ifail, work, Identity, STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error deallocating rwork,', 1)
      !
    ENDIF
    !
    DEALLOCATE(Hamil_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_serial', 'Error deallocating Hamil_save,', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE diag_serial
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE diag_parallel(estmteRt, eigvec_coef)
    !-----------------------------------------------------------------------
    !! Parallel diagonalization using Davidson library from QE
    !-----------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, twopi, ci, eps5, eps6, eps4, eps2, eps8, eps10
    USE global_var,    ONLY : nkf, nktotf
    USE input,         ONLY : nstate_plrn, ethrdg_plrn, &
                              adapt_ethrdg_plrn, init_ethrdg_plrn, nethrdg_plrn, &
                              david_ndim_plrn
    USE io_global,     ONLY : stdout, ionode
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum, mp_bcast, mp_size, mp_max
    USE mp_bands,      ONLY : inter_bgrp_comm, mp_start_bands
    USE mp_bands_util, ONLY : intra_bgrp_comm_ => intra_bgrp_comm, &
                              inter_bgrp_comm_ => inter_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: estmteRt(:)
    !! Polaron eigenvalue
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_coef(:, :)
    !! Polaron eigenvector coefficients
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER  :: ik
    !! k-point counter
    INTEGER  :: ik_global
    !! Global k-point index
    INTEGER  :: ibnd
    !! Electron band counter
    INTEGER  :: itemp
    !! FIXME
    INTEGER  :: jtemp
    !! FIXME
    INTEGER  :: indexkn1
    !! Combined band and k-point index
    INTEGER  :: indexkn2
    !! Combined band and k-point index
    INTEGER  :: ithr
    !! Counter for incremental threshold
    INTEGER  :: nthr
    !! Number of incremental threshold steps
    INTEGER  :: npw
    !! FIXME
    INTEGER  :: npwx
    !! FIXME
    INTEGER  :: dav_iter
    !! FIXME
    INTEGER  :: notcnv
    !! Number of non-converged eigenvalues
    INTEGER  :: btype(nstate_plrn)
    !! FIXME
    INTEGER  :: nhpsi
    !! FIXME
    REAL(KIND = DP) :: ethrdg_init
    !! Initial incremental threshold
    REAL(KIND = DP) :: ethrdg
    !! Threshold for convergence in diagonalization
    COMPLEX(KIND = DP), ALLOCATABLE :: psi(:, :)
    !! Eigenvector
    !
    npw = nkf * nbnd_plrn
    npwx = npw
    CALL mp_max(npwx, inter_pool_comm)
    !
    ALLOCATE(psi(1:npwx, 1:nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_parallel', 'Error allocating psi', 1)
    !
    ! JLB: Option for adaptive threshold
    IF (adapt_ethrdg_plrn) THEN
      ethrdg_init = init_ethrdg_plrn
      nthr = nethrdg_plrn
      IF(ionode) THEN
        WRITE(stdout, "(a)") "     Adaptive threshold on iterative diagonalization activated:"
        WRITE(stdout, "(a)") "     threshold, # iterations, eigenvalue(Ry)"
      ENDIF
    ELSE
      nthr = 1
    ENDIF
    !
    DO ithr = 1, nthr
      psi = czero
      btype(1:nstate_plrn) = 1
      !
      IF (adapt_ethrdg_plrn) THEN
        ethrdg = 10**(LOG10(ethrdg_init) + (ithr - 1) * (LOG10(ethrdg_plrn) - LOG10(ethrdg_init)) /(nthr - 1))
      ELSE
        ethrdg = ethrdg_plrn
      ENDIF
      !
      ! split eigvector (nqtotf) into parallel pieces psi (nkf), contains corresponding part with Hpsi
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
          psi(indexkn1, 1:nstate_plrn) = eigvec_coef(indexkn2, 1:nstate_plrn)
        ENDDO
      ENDDO
      ! inter_bgrp_comm should be some non-existing number,
      ! to make the nodes in bgrp equal to 1
      ! intra_bgrp_comm is parallel PW in pwscf
      ! but here it should be parallel K.
      ! Save them before change them
      itemp = intra_bgrp_comm_
      jtemp = inter_bgrp_comm_
      !
      intra_bgrp_comm_ = inter_pool_comm
      inter_bgrp_comm_ = inter_bgrp_comm
      !
      CALL start_clock('cegterg_prln')
      CALL cegterg( h_psi_plrn, s_psi_plrn, .FALSE., g_psi_plrn, &
        npw, npwx, nstate_plrn, nstate_plrn * david_ndim_plrn, 1, psi, ethrdg, &
        estmteRt, btype, notcnv, .FALSE., dav_iter, nhpsi)
      CALL start_clock('cegterg_prln')
      IF(adapt_ethrdg_plrn .AND. ionode) WRITE(stdout, "(a, E14.6, I6, E14.6)") "   ", ethrdg, dav_iter, estmteRt
      IF(notcnv > 0 .AND. ionode) WRITE(stdout, "(a)") "   WARNING: Some eigenvalues not converged, &
      &check initialization, ethrdg_plrn or try adapt_ethrdg_plrn"
      !
      intra_bgrp_comm_ = itemp
      inter_bgrp_comm_ = jtemp
      !
      eigvec_coef = czero
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
          eigvec_coef(indexkn2, 1:nstate_plrn) = psi(indexkn1, 1:nstate_plrn)
        ENDDO
      ENDDO
      CALL mp_sum(eigvec_coef, inter_pool_comm)
      !
    ENDDO
    !
    DEALLOCATE(psi, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_parallel', 'Error deallocating psi', 1)
    !------------------------------------------------------------------------
    END SUBROUTINE diag_parallel
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    SUBROUTINE write_plrn_dtau_xsf(dtau, nqf1, nqf2, nqf3, filename, species)
    !------------------------------------------------------------------------
    !! Write ionic positions and displacements in XSF format
    !------------------------------------------------------------------------
    USE ep_constants,   ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
    USE mp,             ONLY : mp_sum, mp_bcast
    USE ions_base,      ONLY : nat, ityp, tau, ntypx
    USE cell_base,      ONLY : at, alat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqf1
    !! Fine q-point grid along b1
    INTEGER, INTENT(in) :: nqf2
    !! Fine q-point grid along b2
    INTEGER, INTENT(in) :: nqf3
    !! Fine q-point grid along b3
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    COMPLEX(KIND = DP), INTENT(in) :: dtau(:, :)
    !! Polaron displacements in real space
    INTEGER, INTENT(in), OPTIONAL :: species(50)
    !! Atomic species in unit cell
    !
    ! Local variables
    INTEGER :: ierr
    !! Error index
    INTEGER :: nat_all
    !! Number of atoms in supercell
    INTEGER :: nptotf
    !! Number of unit cells in supercell
    INTEGER :: nqf_s(1:3)
    !! q-point grid
    INTEGER :: iRp
    !! Counter for unit cell vectors in supercell
    INTEGER :: iatm
    !! Counter for atoms in unit cell
    INTEGER :: iatm_all
    !! Counter for atoms in supercell
    INTEGER :: ika
    !! Combined atom and cartesian direction index
    INTEGER :: isp
    !! Atomic species counter
    INTEGER :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER, ALLOCATABLE :: elements(:)
    !! Atomic species array
    REAL(KIND = DP) :: cell(3, 3)
    !! Supercell coordinates
    REAL(KIND = DP) :: shift(1:3)
    !! Shift vector coordinates
    REAL(KIND = DP), ALLOCATABLE :: atoms(:,:)
    !! Atomic position coordinates in supercell
    REAL(KIND = DP), ALLOCATABLE :: displacements(:,:)
    !! Displacement coordinates in supercell
    !
    ! total number of atoms is
    ! (number of atoms in unit cell) x (number of cells in the supercell)
    nptotf =  nqf1 * nqf2 * nqf3
    nqf_s  = (/nqf1, nqf2, nqf3/)
    nat_all = nat * nptotf
    !
    ALLOCATE(atoms(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating atoms', 1)
    ALLOCATE(elements(nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating elements', 1)
    ALLOCATE(displacements(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error allocating displacements', 1)
    !
    atoms = zero
    elements = 0
    displacements = zero
    !
    cell(1:3, 1) = at(1:3, 1) * nqf1
    cell(1:3, 2) = at(1:3, 2) * nqf2
    cell(1:3, 3) = at(1:3, 3) * nqf3
    !
    iatm_all = 0
    DO isp = 1, ntypx
      DO iRp = 1, nptotf
        Rp_vec(1:3) = index_Rp(iRp, nqf_s)
        DO iatm = 1, nat
          IF(ityp(iatm) == isp) THEN
            iatm_all = iatm_all + 1
            ika = (iatm - 1) * 3 + 1
            !Rp(1:3) = (ix - nqf1/2) * at(1:3, 1) + (iy - nqf2/2) * at(1:3, 2) + (iz - nqf3/2) * at(1:3, 3)
            shift(1:3) = Rp_vec(1) * at(1:3, 1) + Rp_vec(2) * at(1:3, 2) + Rp_vec(3) * at(1:3, 3)
            IF(PRESENT(species)) THEN
              elements(iatm_all) = species(ityp(iatm))
            ELSE
              elements(iatm_all) = ityp(iatm)
            ENDIF
            atoms(1:3, iatm_all) = tau(1:3, iatm) + shift(1:3)
            displacements(1:3, iatm_all) = REAL(dtau(iRp, ika:ika + 2))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    cell = cell * alat
    atoms = atoms * alat
    !
    CALL write_xsf_file(filename, cell * bohr2ang, elements, atoms * bohr2ang, displacements * bohr2ang)
    !
    DEALLOCATE(atoms, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating atoms', 1)
    DEALLOCATE(elements, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating elements', 1)
    DEALLOCATE(displacements, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_plrn_dtau_xsf', 'Error deallocating displacements', 1)
    !----------------------------------------------------------------------------------------
    END SUBROUTINE write_plrn_dtau_xsf
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    SUBROUTINE scell_write_plrn_dtau_xsf(dtau, nqtotf_p, nRp_p, Rp_p, as_p, filename, species)
    !----------------------------------------------------------------------------------------
    !! JLB: Write ionic positions and displacements for transformed supercell
    !----------------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
    USE mp,            ONLY : mp_sum, mp_bcast
    USE ions_base,     ONLY : nat, ityp, tau, ntypx
    USE cell_base,     ONLY : at, alat
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !! Output file name
    INTEGER, INTENT(in) :: nqtotf_p
    !! Number of q-points in fine grid
    INTEGER, INTENT(in) :: nRp_p
    !! Number of unit cells within supercell
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! Coordinates of lattice vectors in supercell
    REAL(KIND = DP), INTENT(in) :: as_p(3,3)
    !! Supercell lattice vectors
    COMPLEX(KIND = DP), INTENT(in) :: dtau(:, :)
    !! Polaron displacement coordinates
    INTEGER, INTENT(in), OPTIONAL :: species(50)
    !! Atomic species in unit cell
    !
    ! Local variable
    INTEGER :: ierr
    !! Error index
    INTEGER :: nat_all
    !! Total number of atoms in supercell
    INTEGER :: iRp
    !! Lattice vector counter in supercell
    INTEGER :: iatm
    !! Atom counter in unit cell
    INTEGER :: iatm_all
    !! Atom counter in supercell
    INTEGER :: ika
    !! Combined atom and cartesian direction index
    INTEGER :: isp
    !! Atomic species counter
    INTEGER :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER, ALLOCATABLE :: elements(:)
    !! Atomic species
    REAL(KIND = DP) :: cell(3, 3)
    !! Supercell lattice vectors
    REAL(KIND = DP) :: shift(1:3)
    !! Shift vector coordinates
    REAL(KIND = DP), ALLOCATABLE :: atoms(:,:)
    !! Atomic position coordinates in supercell
    REAL(KIND = DP), ALLOCATABLE :: displacements(:,:)
    !! Atomic displacement coordinates in supercell
    !
    ! total number of atoms is
    ! (number of atoms in unit cell) x (number of cells in the supercell)
    nat_all = nat * nRp_p
    !
    ALLOCATE(atoms(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error allocating atoms', 1)
    ALLOCATE(elements(nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error allocating elements', 1)
    ALLOCATE(displacements(3, nat_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error allocating displacements', 1)
    !
    atoms = zero
    elements = 0
    displacements = zero
    !
    cell(1:3, 1) = as_p(1, 1:3)
    cell(1:3, 2) = as_p(2, 1:3)
    cell(1:3, 3) = as_p(3, 1:3)
    !
    iatm_all = 0
    DO isp = 1, ntypx
      DO iRp = 1, nRp_p
        Rp_vec(1:3) = Rp_p(1:3, iRp)
        DO iatm = 1, nat
          IF(ityp(iatm) == isp) THEN
            iatm_all = iatm_all + 1
            ika = (iatm - 1) * 3 + 1
            !Rp(1:3) = (ix - nqf1/2) * at(1:3, 1) + (iy - nqf2/2) * at(1:3, 2) + (iz - nqf3/2) * at(1:3, 3)
            shift(1:3) = Rp_vec(1) * at(1:3, 1) + Rp_vec(2) * at(1:3, 2) + Rp_vec(3) * at(1:3, 3)
            IF(PRESENT(species)) THEN
              elements(iatm_all) = species(ityp(iatm))
            ELSE
              elements(iatm_all) = ityp(iatm)
            ENDIF
            atoms(1:3, iatm_all) = tau(1:3, iatm) + shift(1:3)
            displacements(1:3, iatm_all) = REAL(dtau(iRp, ika:ika + 2))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    cell = cell * alat
    atoms = atoms * alat
    !
    CALL write_xsf_file(filename, cell * bohr2ang, elements, atoms * bohr2ang, displacements * bohr2ang)
    !
    DEALLOCATE(atoms, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error deallocating atoms', 1)
    DEALLOCATE(elements, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error deallocating elements', 1)
    DEALLOCATE(displacements, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_plrn_dtau_xsf', 'Error deallocating displacements', 1)
    !---------------------------------------------------------------------------
    END SUBROUTINE scell_write_plrn_dtau_xsf
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    SUBROUTINE write_xsf_file(filename, cell, elements, atoms, forces, data_cube)
    !---------------------------------------------------------------------------
    !! Write xsf to file
    !---------------------------------------------------------------------------
    USE io_var,            ONLY : ixsfplrn
    !  
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in):: filename
    !! Output file name
    INTEGER, INTENT(in) :: elements(:)
    !! Atomic species
    REAL(KIND = DP), INTENT(in) :: cell(3, 3)
    !! Supercell lattice vectors
    REAL(KIND = DP), INTENT(in) :: atoms(:, :)
    !! Atomic position coordinates
    REAL(KIND = DP), INTENT(in), OPTIONAL  :: forces(:, :)
    !! Atomic displacement coordinates
    REAL(KIND = DP), INTENT(in), OPTIONAL  :: data_cube(:, :, :)
    !! Data for isosurface
    !
    ! Local variables
    INTEGER :: ix
    !! x-coordinate counter
    INTEGER :: iy
    !! y-coordinate counter
    INTEGER :: iz
    !! z-coordinate counter
    INTEGER :: iatm
    !! Atom counter
    INTEGER :: natm
    !! Number of atoms in supercell
    INTEGER :: shapeTemp(3)
    !! Shape of isosurface data
    !
    natm = UBOUND(elements, DIM = 1)
    !
    OPEN(UNIT = ixsfplrn, FILE = TRIM(filename), FORM = 'formatted', STATUS = 'unknown')
    !
    WRITE(ixsfplrn, '(a)') '#'
    WRITE(ixsfplrn, '(a)') '# Generated by the EPW polaron code'
    WRITE(ixsfplrn, '(a)') '#'
    WRITE(ixsfplrn, '(a)') '#'
    WRITE(ixsfplrn, '(a)') 'CRYSTAL'
    WRITE(ixsfplrn, '(a)') 'PRIMVEC'
    WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 1)
    WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 2)
    WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 3)
    WRITE(ixsfplrn, '(a)') 'PRIMCOORD'
    ! The second number is always 1 for PRIMCOORD coordinates,
    ! according to http://www.xcrysden.org/doc/XSF.html
    WRITE(ixsfplrn, '(2i6)')  natm, 1
    !
    DO iatm = 1, natm
      IF (PRESENT(forces)) THEN
        WRITE(ixsfplrn,'(I3, 3x, 3f15.9, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm), forces(1:3, iatm)
      ELSE
        WRITE(ixsfplrn,'(I3, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm)
      ENDIF
    ENDDO
    !
    IF(PRESENT(data_cube)) THEN
      shapeTemp = SHAPE(data_cube)
      WRITE(ixsfplrn, '(/)')
      WRITE(ixsfplrn, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
      WRITE(ixsfplrn, '(3i6)') SHAPE(data_cube)
      WRITE(ixsfplrn, '(3f12.6)') 0.0, 0.0, 0.0
      WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 1)
      WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 2)
      WRITE(ixsfplrn, '(3f12.7)') cell(1:3, 3)
      ! TODO: data cube is probably to large to take in the same way of lattice information
      ! May be usefull and implemented in the furture
      WRITE(ixsfplrn, *) (((data_cube(ix, iy, iz), ix = 1, shapeTemp(1)), &
        iy = 1, shapeTemp(2)), iz = 1, shapeTemp(3))
      WRITE(ixsfplrn, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
    ENDIF
    CLOSE(ixsfplrn)
    !----------------------------------------------------------------------------------------
    END SUBROUTINE write_xsf_file
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    SUBROUTINE write_plrn_wf(eigvec_coef, filename, enk_all)
    !----------------------------------------------------------------------------------------
    !! Write polaron wavefunction coefficients
    !----------------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, ryd2ev, ryd2mev
    USE global_var,    ONLY : nktotf
    USE io_var,        ONLY : iwfplrn
    USE input,         ONLY : nstate_plrn, nkf1, nkf2, nkf3, nbndsub, scell_mat_plrn
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    REAL(KIND = DP), INTENT(in), OPTIONAL :: enk_all(:, :)
    !! KS eigenvalues on global fine grid
    COMPLEX(KIND = DP), INTENT(in) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in Bloch (Ank) or Wannier (Amp) basis
    !
    ! Local variables
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: nbnd_out
    !! Number of bands in polaron wave function expansion
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: iplrn
    !! Polaron state counter
    !
    IF(PRESENT(enk_all)) THEN
      nbnd_out = nbnd_plrn
    ELSE
      nbnd_out = nbndsub
    ENDIF
    !
    OPEN(UNIT = iwfplrn, FILE = TRIM(filename))
    !
    IF (scell_mat_plrn) THEN
      WRITE(iwfplrn, '(a, 3I10)') 'Scell', nktotf, nbndsub, nstate_plrn
    ELSE
      WRITE(iwfplrn, '(6I10)') nkf1, nkf2, nkf3, nktotf, nbndsub, nstate_plrn
    ENDIF
    !
    DO ik = 1, nktotf
      DO ibnd = 1, nbnd_out
        DO iplrn = 1, nstate_plrn
          indexkn1 = (ik - 1) * nbnd_out + ibnd
          IF (PRESENT(enk_all)) THEN
            WRITE(iwfplrn, '(2I5, 4f15.7)') ik, ibnd, enk_all(select_bands_plrn(ibnd), ik) * ryd2ev, &
               eigvec_coef(indexkn1, iplrn), ABS(eigvec_coef(indexkn1, iplrn))
          ELSE
            WRITE(iwfplrn, '(2f15.7)') eigvec_coef(indexkn1, iplrn)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    CLOSE(iwfplrn)
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE write_plrn_wf
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE read_plrn_wf_grid(nkf1_p, nkf2_p, nkf3_p, nktotf_p, &
               nbndsub_p, nplrn_p, filename, scell)
    !----------------------------------------------------------------------------
    !! Read k-point grid in which polaron wave function has been written to file.
    !! Needed for correct allocation when interp_plrn_wf.
    !----------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : iwfplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    LOGICAL, INTENT(in), OPTIONAL  :: scell
    !! .true. for non-diagonal supercell calculation
    INTEGER, INTENT(out) :: nkf1_p
    !! Number of k-points in fine grid along b1
    INTEGER, INTENT(out) :: nkf2_p
    !! Number of k-points in fine grid along b2
    INTEGER, INTENT(out) :: nkf3_p
    !! Number of k-points in fine grid along b3
    INTEGER, INTENT(out) :: nktotf_p
    !! Total number of k-points in fine grid
    INTEGER, INTENT(out) :: nbndsub_p
    !! Number of bands in polaron wave function expansion
    INTEGER, INTENT(out) :: nplrn_p
    !! Number of polaron states    
    !
    ! Local variables
    LOGICAL :: scell_
    !! Dummy variable to set a default value (.FALSE.) for scell
    CHARACTER(LEN = 5) :: dmmy
    !! Dummy variable to read from scell wf file
    !
    !
    OPEN(UNIT = iwfplrn, FILE = TRIM(filename))
    !
    scell_ = .FALSE.
    IF (PRESENT(scell)) scell_ = scell
    !
    IF (scell_) THEN
      READ(iwfplrn, '(a, 3I10)') dmmy, nktotf_p, nbndsub_p, nPlrn_p
      ! nkf1_p, nkf2_p, nkf3_p should never be called if scell=.true.
      ! Just assigning an arbitrary value
      nkf1_p = 0
      nkf2_p = 0
      nkf3_p = 0
    ELSE
      READ(iwfplrn, '(6I10)') nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nPlrn_p
      IF(nkf1_p * nkf2_p * nkf3_p /= nktotf_p) THEN
        CALL errore("read_plrn_wf_grid", 'Amp.plrn'//'Not generated from the uniform grid!', 1)
      ENDIF
    ENDIF
    !
    CLOSE(iwfplrn)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_wf_grid
    !-----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE read_plrn_wf(eigvec_coef, nkf1_p, nkf2_p, nkf3_p, nktotf_p, &
               nbndsub_p, nplrn_p, filename, enk_all)
    !----------------------------------------------------------------------------
    !! Read polaron wavefunction coefficients.
    !----------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : iwfplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    COMPLEX(KIND = DP), INTENT(out) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in Bloch basis, Ank
    INTEGER, INTENT(in) :: nkf1_p
    !! Number of k-points in fine grid along b1
    INTEGER, INTENT(in) :: nkf2_p
    !! Number of k-points in fine grid along b2
    INTEGER, INTENT(in) :: nkf3_p
    !! Number of k-points in fine grid along b3
    INTEGER, INTENT(in) :: nktotf_p
    !! Total number of k-points in fine grid
    INTEGER, INTENT(in) :: nbndsub_p
    !! Number of bands in polaron wave function expansion
    INTEGER, INTENT(in) :: nplrn_p
    !! Number of polaron states    
    REAL(KIND = DP), INTENT(in), OPTIONAL :: enk_all(:, :)
    !! KS eigenvalues in global fine k-point grid
    !
    ! Local variables
    INTEGER:: ierr
    !! Error status
    INTEGER :: indexkn1
    !! Combined band and k-poin index
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: i1
    !! Dummy integer to read from file
    INTEGER :: i2
    !! Dummy integer to read from file
    REAL(KIND = DP) :: r1
    !! Dummy real to read from file
    !
    OPEN(UNIT = iwfplrn, FILE = TRIM(filename))
    !
    ! First line should have been read already,
    ! so that eigvec has been properly allocated
    READ(iwfplrn, *)
    !
    eigvec_coef = czero
    DO ik = 1, nktotf_p
      DO ibnd = 1, nbndsub_p
        DO iplrn = 1, nplrn_p
          indexkn1 = (ik - 1) * nbndsub_p + ibnd
          IF(PRESENT(enk_all)) THEN ! Ank.plrn is read
            READ(iwfplrn, '(2I5, 3f15.7)') i1, i2, r1, eigvec_coef(indexkn1, iplrn)
          ELSE ! Amp.plrn is read
            READ(iwfplrn, '(2f15.7)') eigvec_coef(indexkn1, iplrn)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    CLOSE(iwfplrn)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_wf
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE write_plrn_bmat(bqv_coef, filename, enk_all)
    !-----------------------------------------------------------------------------------
    !!
    !! Write Bqv coeffients (or dtau displacements) and phonon frequencies to filename
    !!
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, ryd2ev, ryd2mev
    USE global_var,    ONLY : nqtotf
    USE io_var,        ONLY : idtauplrn
    USE input,         ONLY : nqf1, nqf2, nqf3, scell_mat_plrn
    USE mp,            ONLY : mp_sum, mp_bcast
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    COMPLEX(KIND = DP), INTENT(in) :: bqv_coef(:, :)
    !! Polaron displacement coefficients
    REAL(KIND = DP), INTENT(in), OPTIONAL :: enk_all(:, :)
    !! KS eigenvalues in fine grid
    !
    ! Local variables
    INTEGER :: iq
    !! q-point counter
    INTEGER :: imode
    !! Phonon mode counter
    !
    OPEN(UNIT = idtauplrn, FILE = TRIM(filename))
    IF (scell_mat_plrn) THEN
      WRITE(idtauplrn, '(a, 2I10)') 'Scell', nqtotf, nmodes
    ELSE
      WRITE(idtauplrn, '(5I10)') nqf1, nqf2, nqf3, nqtotf, nmodes
    ENDIF
    !
    DO iq = 1, nqtotf
      DO imode = 1, nmodes ! p
        IF (PRESENT(enk_all)) THEN ! write Bqv
          !JLB: Changed format for improved accuracy
          WRITE(idtauplrn, '(2I5, 4ES18.10)') iq, imode, enk_all(imode, iq) * ryd2mev, &
                                                  bqv_coef(iq, imode), ABS(bqv_coef(iq, imode))
        ELSE ! write \dtau
          !JLB: Changed format for improved accuracy
          WRITE(idtauplrn, '(2ES18.10)') bqv_coef(iq, imode)
        ENDIF
      ENDDO
    ENDDO
    CLOSE(idtauplrn)
    !---------------------------------------------------------------------------------
    END SUBROUTINE write_plrn_bmat
    !---------------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE read_plrn_dtau_grid(nqf1_p, nqf2_p, nqf3_p, nqtotf_p, &
               nmodes_p, filename, scell)
    !----------------------------------------------------------------------------
    !! Read q-point grid in which polaron displacements have been written to file.
    !! Needed for correct allocation when interp_plrn_bq.
    !----------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : idtauplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    LOGICAL, INTENT(in), OPTIONAL  :: scell
    !! .true. for non-diagonal supercell calculation
    INTEGER, INTENT(out) :: nqf1_p
    !! Number of q-points in fine grid along b1
    INTEGER, INTENT(out) :: nqf2_p
    !! Number of q-points in fine grid along b2
    INTEGER, INTENT(out) :: nqf3_p
    !! Number of q-points in fine grid along b3
    INTEGER, INTENT(out) :: nqtotf_p
    !! Total number of q-points in fine grid
    INTEGER, INTENT(out) :: nmodes_p
    !! Number of phonon modes 
    !
    ! Local variables
    LOGICAL :: scell_
    !! Dummy variable to set a default value (.FALSE.) for scell
    CHARACTER(LEN = 5) :: dmmy
    !! Dummy variable to read from scell wf file
    !
    !
    OPEN(UNIT = idtauplrn, FILE = 'dtau.plrn')
    !
    scell_ = .FALSE.
    IF (PRESENT(scell)) scell_ = scell
    !
    IF (scell_) THEN
      READ(idtauplrn, '(a, 3I10)') dmmy, nqtotf_p, nmodes_p
      ! nkf1_p, nkf2_p, nkf3_p should never be called if scell=.true.
      ! Just assigning an arbitrary value
      nqf1_p = 0
      nqf2_p = 0
      nqf3_p = 0
    ELSE
      READ(idtauplrn, '(6I10)') nqf1_p, nqf2_p, nqf3_p, nqtotf_p, nmodes_p
      IF(nqf1_p * nqf2_p * nqf3_p /= nqtotf_p) THEN
        CALL errore("read_plrn_dtau_grid", 'dtau.plrn'//'Not generated from the uniform grid!', 1)
      ENDIF
    ENDIF
    !
    CLOSE(idtauplrn)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_dtau_grid
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE read_plrn_dtau(dtau_read, nqtotf_p, nmodes_p,&
                  filename, scell, wfreq)
    !-----------------------------------------------------------------------------------
    !! Read displacement coefficients in phonon (Bqv) or real-space (dtau) basis from file
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : idtauplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! Output file name
    LOGICAL, INTENT(in), OPTIONAL :: scell
    !! .true. if non-diagonal supercell has been used in polaron calculation
    INTEGER, INTENT(in) :: nqtotf_p
    !! Total number of q-points in fine grid
    INTEGER, INTENT(in) :: nmodes_p
    !! Number of phonon modes
    REAL(KIND = DP), INTENT(in), OPTIONAL :: wfreq(:, :)
    !! Phonon frequencies
    COMPLEX(KIND = DP), INTENT(out) :: dtau_read(:, :)
    !! Polaron displacement coefficients
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    INTEGER :: iatm
    !! Atom counter
    INTEGER :: iq
    !! q-point counter
    INTEGER :: i1
    !! Dummy integer to read from file
    INTEGER :: i2
    !! Dummy integer to read from file
    REAL(KIND = DP) :: r1
    !! Dummy real to read from file
    !
    OPEN(UNIT = idtauplrn, FILE = TRIM(filename))
    !
    ! JLB:
    ! First line should have been read already,
    ! so that dtau has been properly allocated
    READ(idtauplrn, *)
    !
    dtau_read = czero
    DO iq = 1, nqtotf_p
      DO iatm = 1, nmodes_p
        IF(PRESENT(wfreq)) THEN ! read Bqv
          !JLB: Changed format for improved accuracy
          READ(idtauplrn, '(2I5, 3ES18.10)') i1, i2, r1, dtau_read(iq, iatm)
        ELSE ! read \dtau
          !JLB: Changed format for improved accuracy
          READ(idtauplrn, '(2ES18.10)') dtau_read(iq, iatm)
        ENDIF
      ENDDO
    ENDDO
    CLOSE(idtauplrn)
    !--------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_dtau
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    SUBROUTINE plrn_eigvec_tran(ttype, t_rev, eigvecin, nkf1_p, nkf2_p, nkf3_p, &
               nbndsub_p, nrr_k, ndegen_k, irvec_r, dims, eigvecout, ip_center)
    !--------------------------------------------------------------------------------
    !! Fourier transform from eigvecin to eigvecout
    !! ttype is 'Bloch2Wan' or 'Wan2Bloch'
    !! Parallel version, each pool calculates its own k point set (nkf),
    !! then the mp_sum is used to sum over different pools.
    !! require the correct initialization of Rp_array
    !--------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, twopi, ci, cone, two
    USE global_var,    ONLY : nkf, xkf, chw, nktotf
    USE input,         ONLY : nstate_plrn, nbndsub
    USE wannier2bloch, ONLY : hamwan2bloch !!=> hamwan2bloch_old
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! Transformation direction, 'Bloch2Wan' or 'Wan2Bloch'
    LOGICAL, INTENT(in) :: t_rev
    !! .true. if time reversal symmetry is to be imposed in the Ank coefficients
    INTEGER, INTENT(in) :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER, INTENT(in) :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER, INTENT(in) :: nkf3_p
    !! Fine k-point grid along b2
    INTEGER, INTENT(in) :: nbndsub_p
    !! Number of bands
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(in), OPTIONAL :: ip_center(1:3)
    !! Center of polaron wave function, to shift supercell accordingly
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(out) :: eigvecout(:, :)
    !! Output wave function coefficients
    COMPLEX(KIND = DP), INTENT(in) :: eigvecin(:, :)
    !! Input wave function coefficients
    !
    ! Local variables
    LOGICAL :: is_mirror
    !! .true. if k-point is a time-reversal mirror point
    INTEGER :: itype
    !! Transformation direction
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: ikpg
    !! Index of mirror k-point
    INTEGER :: ikglob
    !! Global inner k-point counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: jbnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    INTEGER :: i_vec(3)
    !! Shifted lattice vector coordinates
    INTEGER :: center_shift(1:3)
    !! Shift coordinates
    INTEGER :: nkf_p(3)
    !! k-point coordinates for shift
    REAL(KIND = DP) :: xxk(3)
    !! k-point coordinates
    REAL(KIND = DP) :: etf_tmp(nbndsub)
    !! Eigenvalues after interpolated KS Hamiltonian diagonalization
    COMPLEX(KIND = DP) :: ctemp
    !! Exponential prefactor
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! U_{mn} matrices after interpolated KS Hamiltonian diagonalization
    COMPLEX(KIND = DP) :: cfac(nrr_k)
    !!! Exponential factor for Hamiltonian transformation
    !
    nkf_p(1:3) = (/nkf1_p, nkf2_p, nkf3_p/)
    IF (nbndsub_p /= nbndsub) CALL errore('plrnwfwan2bloch','Different bands included in last calculation!', 1)
    IF (ttype == 'Bloch2Wan') THEN
      itype =  1
    ELSE IF (ttype == 'Wan2Bloch') THEN
      itype = -1
    ELSE
      CALL errore('plrn_eigvec_tran', 'Illegal translate form; should be Bloch2Wan or Wan2Bloch!', 1)
    ENDIF
    !
    IF(PRESENT(ip_center)) THEN
      center_shift(1:3) = nkf_p / 2 - ip_center
    ELSE
      center_shift(1:3) = 0
    ENDIF
    !! itype =  1 : Bloch2Wan: A_{mp} =  \frac{1}{N_p} \sum_{nk}A_{nk} \exp\left(ik\cdot R_p\right)U^\dagger_{mnk}
    !! itype = -1 : Wan2Bloch: A_{nk} = \sum_{mp}A_{mp}\exp(-ik\cdot R_p) U_{mnk}
    !! ibnd -> m, jbnd -> n
    !! R_p from 1 to nkf1/2/3_p, note that loop in the sequence of ix, iy, and iz,
    !! This sequence need to be consistent every time transpose between eigvec_wann and eigvec
    eigvecout = czero
    DO ik = 1, nkf
      xxk = xkf(1:3, 2 * ik - 1)
      ik_global = ikqLocal2Global(ik, nktotf)
      !
      IF (t_rev) THEN
        ikpg = kpg_map(ik_global)
        is_mirror = (ik_global > ikpg)
      END IF
      !
      CALL get_cfac(xxk, nrr_k, irvec_r, cfac)
      !
      CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
         etf_tmp, chw, cfac, is_mirror)
      !
      IF(itype == 1) cufkk(1:nbndsub, 1:nbndsub) = CONJG(TRANSPOSE(cufkk(1:nbndsub, 1:nbndsub)))
      DO iplrn = 1, nstate_plrn
        !ikglob = 0
        !! loop over all Wannier position p
        IF (nkf1_p == 0 .OR. nkf2_p == 0 .OR. nkf3_p == 0) THEN
          CALL errore('plrn_eigvec_tran','Wrong k grid, use nkf1/2/3 to give k grid!', 1)
        ENDIF
        DO ikglob = 1, nkf1_p * nkf2_p * nkf3_p
          i_vec(1:3) = MODULO(index_Rp(ikglob, nkf_p) + center_shift, nkf_p)
          ctemp = EXP(CMPLX(0.0_DP, twopi * DOT_PRODUCT(xxk, i_vec), KIND = DP))
          DO ibnd = 1, nbndsub_p ! loop over all Wannier state m
            DO jbnd = 1, nbnd_plrn ! loop over all Bloch state n
              indexkn1 = (ikglob - 1) * nbndsub + ibnd !mp
              indexkn2 = (ik_global - 1) * nbnd_plrn + jbnd !nk
              SELECT CASE(itype)
                CASE(1)  ! Bloch2Wan !
                   eigvecout(indexkn1, iplrn) = eigvecout(indexkn1, iplrn) + &
                      eigvecin(indexkn2, iplrn) * ctemp / nktotf * cufkk(ibnd, select_bands_plrn(jbnd)) !JLB: Conjugate transpose taken above!
                CASE(-1) ! Wan2Bloch !
                   eigvecout(indexkn2, iplrn) = eigvecout(indexkn2, iplrn) + &
                      eigvecin(indexkn1, iplrn) * CONJG(ctemp) * cufkk(select_bands_plrn(jbnd), ibnd) !JLB
              END SELECT
            ENDDO ! jbnd
          ENDDO ! ibnd
        ENDDO ! ikglob
      ENDDO ! iplrn
    ENDDO ! ik
    ! MPI sum due to the loop ik is within local k set
    CALL mp_sum(eigvecout, inter_pool_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE plrn_eigvec_tran
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE scell_plrn_eigvec_tran(ttype, t_rev, eigvecin, nktotf_p, nRp_p, Rp_p, &
                                nbndsub_p, nrr_k, ndegen_k, irvec_r, dims, eigvecout)
    !-----------------------------------------------------------------------------------
    !! JLB: Fourier transform for non-diagonal supercells
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : czero, twopi, ci, cone, two
    USE global_var,    ONLY : nkf, xkf, chw, nktotf
    USE input,         ONLY : nstate_plrn, nbndsub
    USE wannier2bloch, ONLY : hamwan2bloch !!=> hamwan2bloch_old
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! Transformation direction, 'Bloch2Wan' or 'Wan2Bloch'
    LOGICAL, INTENT(in) :: t_rev
    !! .true. if time-reversal symmetry is to be imposed
    INTEGER, INTENT(in) :: nktotf_p
    !! Number of k-points in fine grid
    INTEGER, INTENT(in) :: nRp_p
    !! Number of unit cells within supercell
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! Lattice vector coefficients in supercell
    INTEGER, INTENT(in) :: nbndsub_p
    !! Number of bands in polaron expansion
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(out) :: eigvecout(:, :)
    !! Output wave function coefficients
    COMPLEX(KIND = DP), INTENT(in) :: eigvecin(:, :)
    !! Input wave function coefficients
    !
    ! Local Variables
    REAL(KIND = DP) :: xxk(3)
    !! k-point coordinates
    REAL(KIND = DP) :: etf_tmp(nbndsub)
    !! Eigenvalues after interpolated KS Hamiltonian diagonalization
    COMPLEX(KIND = DP) :: ctemp
    !! Exponential prefactor
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! U_{mn} matrices after interpolated KS Hamiltonian diagonalization
    INTEGER :: itype
    !! Transformation direction
    INTEGER :: ik
    !! k-point counter
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: ikpg
    !! Index of mirror k-point
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: jbnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: indexkn2
    !! Combined band and k-point index
    LOGICAL :: is_mirror
    !! .true. if k-point is time-reversal mirror point
    INTEGER :: ierr
    !! Error status
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:, :, :)
    !! Exponential prefactor
    !
    IF (nbndsub_p /= nbndsub) CALL errore('scell_plrn_eigvec_tran','Different bands included in last calculation!',1)
    IF (ttype == 'Bloch2Wan') THEN
      itype =  1
    ELSEIF (ttype == 'Wan2Bloch') THEN
      itype = -1
    ELSE
      CALL errore('scell_plrn_eigvec_tran', 'Illegal translate form; should be Bloch2Wan or Wan2Bloch!', 1)
    ENDIF
    !
    ALLOCATE(cfac(nrr_k, dims, dims), STAT = ierr)
    IF(ierr /= 0) CALL errore('scell_plrn_eigvec_tran', 'Error allocating cfac', 1)
    !
    !! itype =  1 : Bloch2Wan: A_{mp} =  \frac{1}{N_p} \sum_{nk}A_{nk} \exp\left(ik\cdot R_p\right)U^\dagger_{mnk}
    !! itype = -1 : Wan2Bloch: A_{nk} = \sum_{mp}A_{mp}\exp(-ik\cdot R_p) U_{mnk}
    !! ibnd -> m, jbnd -> n
    !! R_p from 1 to nktotf_p
    !! This sequence need to be consistent every time transpose between eigvec_wann and eigvec
    eigvecout = czero
    DO ik = 1, nkf
      xxk = xkf(1:3, 2 * ik - 1)
      ik_global = ikqLocal2Global(ik, nktotf)
      is_mirror = (t_rev .AND. (ik_global > ikpg))
      !
      CALL get_cfac(xxk, nrr_k, irvec_r, cfac)
      CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
         etf_tmp, chw, cfac, is_mirror)
      IF(itype == 1) cufkk(1:nbndsub, 1:nbndsub) = CONJG(TRANSPOSE(cufkk(1:nbndsub, 1:nbndsub)))
      !
      DO iplrn = 1, nstate_plrn
        !icount = 0
        !! loop over all Wannier position p
        DO iRp = 1, nRp_p
          ctemp = EXP(twopi * ci * DOT_PRODUCT(xxk, Rp_p(1:3, iRp)))
          DO ibnd = 1, nbndsub_p ! loop over all Wannier state m
            DO jbnd = 1, nbnd_plrn ! loop over all Bloch state n
              indexkn1 = (iRp - 1) * nbndsub + ibnd !mp
              indexkn2 = (ik_global - 1) * nbnd_plrn + jbnd !nk
              SELECT CASE(itype)
                CASE(1)  ! Bloch2Wan !
                   eigvecout(indexkn1, iplrn) = eigvecout(indexkn1, iplrn) + &
                      eigvecin(indexkn2, iplrn) * ctemp / nktotf * cufkk(ibnd, select_bands_plrn(jbnd))
                CASE(-1) ! Wan2Bloch !
                   eigvecout(indexkn2, iplrn) = eigvecout(indexkn2, iplrn) + &
                      eigvecin(indexkn1, iplrn) * CONJG(ctemp) * cufkk(select_bands_plrn(jbnd), ibnd)
              END SELECT
            ENDDO ! jbnd
          ENDDO ! ibnd
        ENDDO
      ENDDO !iplrn
    ENDDO ! ik
    ! MPI sum due to the loop ik is within local k set
    CALL mp_sum(eigvecout, inter_pool_comm)
    !
    DEALLOCATE(cfac, STAT = ierr)
    IF(ierr /= 0) CALL errore('scell_plrn_eigvec_tran', 'Error deallocating cfac', 1)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE scell_plrn_eigvec_tran
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE interp_plrn_wf(nrr_k, ndegen_k, irvec_r, dims, scell)
    !-----------------------------------------------------------------------------------
    !! Interpolate polaron wave function coeffcients (Ank) and write to Ank.band.plrn.
    !! Mostly used to visualize contributions from different bands and k-points.
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, ryd2ev, czero
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE io_var,        ONLY : iwfplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! Number of electronic WS points
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    LOGICAL, INTENT(in), OPTIONAL  :: scell
    !! .true. for non-diagonal supercell calculation
    !
    ! Local variables
    INTEGER :: ierr
    !! Error code when reading file
    INTEGER :: i_center(2)
    !! Index of polaron center lattice vector
    INTEGER :: ip_center(3)
    !! Coordinates of polaron center
    INTEGER :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER :: nkf3_p
    !! Fine k-point grid along b3
    INTEGER :: nktotf_p
    !! Number of k-points in fine grid
    INTEGER :: nbndsub_p
    !! Number of bands in polaron expansion
    INTEGER :: nplrn_p
    !! Number of polaron states
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! Polaron wave function coefficients in Wannier basis, Amp
    !
    IF (ionode) WRITE(stdout, "(5x, a)") "Start of interpolation of electronic band structure."
    !
    ! read Amp.plrn, save eigvec_wan for the latter use
    IF(ionode) THEN
      CALL read_plrn_wf_grid(nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(nkf1_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf2_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf3_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nktotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nbndsub_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nplrn_p,  meta_ionode_id, world_comm)
    !
    ALLOCATE(eigvec_wan(nktotf_p * nbndsub_p, nplrn_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_wf', 'Error allocating eigvec_wan', 1)    
    !
    IF (ionode) THEN
      CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(eigvec_wan, meta_ionode_id, world_comm)
    !
    i_center = MAXLOC(ABS(eigvec_wan))
    !
    ip_center = index_Rp(i_center(1) / nbndsub_p + 1, (/nkf1_p, nkf2_p, nkf3_p/))
    WRITE(stdout, '(5x, a, i8, 3i5)') "The largest Amp ", i_center(1), ip_center
    !
    ! JLB: kpg_map cannot be generally defined in interpolation k-paths,
    !      thus t_rev set to .false.
    CALL plrn_eigvec_tran('Wan2Bloch', .FALSE., eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
       nrr_k, ndegen_k, irvec_r, dims, eigvec, ip_center)
    !
    CALL write_plrn_wf(eigvec, 'Ank.band.plrn', etf_all)
    !
    DEALLOCATE(eigvec_wan, STAT = ierr)
    IF(ierr /= 0) CALL errore('interp_plrn_wf', 'Error deallocating eigvec_wan', 1)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE interp_plrn_wf
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE interp_plrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws, scell)
    !-----------------------------------------------------------------------------------
    !! Interpolate polaron displacements coefficients (Bqv) and write to Bmat.band.plrn.
    !! Mostly used to visualize contributions from each phonon mode and q-point.
    !-----------------------------------------------------------------------------------
    USE global_var,    ONLY : wf, nqtotf
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : czero
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE io_var,        ONLY : idtauplrn
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER,  INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    LOGICAL, INTENT(in), OPTIONAL  :: scell
    !! .true. for non-diagonal supercell calculation
    !
    ! Local variables
    CHARACTER(LEN = 5) :: dmmy
    !! Dummy variables read from file
    INTEGER :: nqf1_p
    !! Fine q-point grid along b1
    INTEGER :: nqf2_p
    !! Fine q-point grid along b2
    INTEGER :: nqf3_p
    !! Fine q-point grid along b3
    INTEGER :: nqtotf_p
    !! Number of q-points in fine grid
    INTEGER :: nmodes_p
    !! Number of phonon modes
    INTEGER :: ierr
    !! Error status
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: ina
    !! Atom counter
    INTEGER :: i_center(2)
    !! Index of polaron center
    INTEGER :: ip_center(3)
    !! Coordinates of polaron center
    COMPLEX(KIND = DP), ALLOCATABLE :: bqv_coef(:,:)
    !! Polaron displacement coefficients in phonon basis, Bqv
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! Polaron displacements in real space
    REAL(KIND = DP),    ALLOCATABLE :: dtau_r(:, :)
    !! Auxiliary polaron displacements to find polaron center
    !
    IF(ionode) THEN
      CALL read_plrn_dtau_grid(nqf1_p, nqf2_p, nqf3_p, nqtotf_p, nmodes_p, 'dtau.plrn')
    ENDIF
    CALL mp_bcast(nqf1_p,   meta_ionode_id, world_comm)
    CALL mp_bcast(nqf2_p,   meta_ionode_id, world_comm)
    CALL mp_bcast(nqf3_p,   meta_ionode_id, world_comm)
    CALL mp_bcast(nqtotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nmodes_p, meta_ionode_id, world_comm)
    !
    ALLOCATE(dtau(nqtotf_p, nmodes_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating dtau', 1)    
    !
    IF (ionode) THEN
      CALL read_plrn_dtau(dtau, nqtotf_p, nmodes_p, 'dtau.plrn')
    END IF
    CALL mp_bcast(dtau, meta_ionode_id, world_comm)
    !
    ! Locate max displacement to center supercell
    ALLOCATE(dtau_r(nqtotf_p, nmodes/3), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating dtau_e', 1)
    dtau_r = czero
    DO iRp = 1, nqtotf_p
      DO ina = 1, nmodes / 3 ! ika -> kappa alpha
        dtau_r(iRp, ina) = NORM2(REAL(dtau(iRp, (ina - 1) * 3 + 1:ina * 3)))
      ENDDO
    ENDDO
    i_center = MAXLOC(ABS(dtau_r))
    ip_center = index_Rp(i_center(1), (/nqf1_p, nqf2_p, nqf3_p/))
    DEALLOCATE(dtau_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating dtau_r', 1)
    !
    ALLOCATE(bqv_coef(nqtotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating Bmat', 1)
    bqv_coef = czero
    !
    CALL plrn_bmat_tran('Dtau2Bmat', .false., dtau, nqf1_p, nqf2_p, nqf3_p, &
       nrr_q, ndegen_q, irvec_q, rws, nrws, bqv_coef, ip_center)
    !
    IF (ionode) CALL write_plrn_bmat(bqv_coef, 'Bmat.band.plrn', wf)
    !
    DEALLOCATE(dtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating dtau', 1)
    DEALLOCATE(bqv_coef, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating Bmat', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE interp_plrn_bq
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE plrn_bmat_tran(ttype, t_rev, mat_in, nqf1_p, nqf2_p, nqf3_p, &
          nrr_q, ndegen_q, irvec_q, rws, nrws, mat_out, ip_center, acoustic_plrn)
    !-----------------------------------------------------------------------------------
    !! Fourier transform between Bmat and dtau,
    !! Eq.(39) of PRB 99, 235139 (2019).
    !! Dtau2Bmat : B_{q\nu} = -1/N_p\sum_{\kappa\alpha p}C_{q\kappa \nu}\Delta\tau_{\kappa\alpha p}  e_{\kappa\alpha\nu}(q)\exp(iqR_p)
    !! Bmat2Dtau : \Delta \tau_{\kappa\alpha p} = -\sum_{q\nu} 1/(C_{q\kappa \nu}) B^*_{q\nu} e_{\kappa\alpha,\nu}(q) \exp(iqR_p)
    !! C_{q\kappa \nu} = N_p\left(\frac{M_k\omega_{q\nu}}{2\hbar}\right)^{\frac{1}{2}} = N_p(M_k)^{\frac{1}{2}}D_{q\nu}
    !! D_{q \nu} = \left(\frac{\omega_{q\nu}}{2\hbar}\right)^{\frac{1}{2}}
    !-----------------------------------------------------------------------------------
    USE global_var,    ONLY : xqf, nqtotf, wf
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : eps8, czero, one, two, twopi, zero, ci, cone
    USE ions_base,     ONLY : amass, ityp
    USE wannier2bloch, ONLY : dynwan2bloch, dynifc2blochf
    USE input,         ONLY : lifc, type_plrn, eps_acoustic
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE parallelism,   ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! Transformation direction, 'Bloch2Wan' or 'Wan2Bloch'
    LOGICAL, INTENT(in) :: t_rev
    !! .true. if time-reversal symmetry is to be imposed in the Bqv coefficients
    INTEGER, INTENT(in) :: nqf1_p
    !! Fine q-point grid along b1
    INTEGER, INTENT(in) :: nqf2_p
    !! Fine q-point grid along b2
    INTEGER, INTENT(in) :: nqf3_p
    !! Fine q-point grid along b3
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    INTEGER, INTENT(in), OPTIONAL :: ip_center(1:3)
    !! Coordinates of polaron center, for shifting supercell
    REAL(KIND = DP), INTENT(in), OPTIONAL :: acoustic_plrn
    !! the cutoff frequency of acoustic phonon modes in dtau.acoustic.plrn.xsf
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    COMPLEX(KIND = DP), INTENT(in) :: mat_in(:, :)
    !! Input matrix with Bqv/dtau coefficients
    COMPLEX(KIND = DP), INTENT(out) :: mat_out(:, :)
    !! Output matrix with dtau/Bqv coefficients
    !
    ! Local variables
    LOGICAL :: mirror_q
    !! .true. if q1 is TR mirror of another q2 point
    INTEGER :: iq
    !! q-point counter
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: itype
    !! Transformation direction
    INTEGER :: ika
    !! Combined atom and cartesian direction counter
    INTEGER :: ip_start
    !! Initial lattice vector in this pool
    INTEGER :: ip_end
    !! Final lattice vector in this pool
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: nqf_p(1:3)
    !! Fine q-point grid
    INTEGER :: ina
    !! Atom counter
    INTEGER :: nptotf
    !! Lattice vector counter
    INTEGER :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER :: center_shift(1:3)
    !! Coordinates of shift to center supercell around polaron
    REAL(KIND = DP) :: xxq(3)
    !! q-point coordinate
    REAL(KIND = DP) :: xxq_r(3)
    !! q-point coordinate in case time-reversal has to be taken
    REAL(KIND = DP) :: ctemp
    !! Prefactor
    REAL(KIND = DP) :: w2(nmodes)
    !! Phonon frequency squared
    COMPLEX(KIND = DP) :: dtemp
    !! Prefactor
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! Phonon eigenvectors
    !
    nptotf = nqf1_p * nqf2_p * nqf3_p
    nqf_p(1:3) = (/nqf1_p, nqf2_p, nqf3_p/)
    !
    IF (nptotf <= 0) CALL errore('plrn_eigvec_tran', 'Use correct .plrn file with nqf1_p \= 0!', 1)
    IF (ttype == 'Bmat2Dtau') THEN
      itype =  1
    ELSE IF (ttype == 'Dtau2Bmat') THEN
      itype = -1
    ELSE
      CALL errore('plrn_eigvec_tran', 'Illegal translation form; should be Bmat2Dtau or Dtau2Bmat!', 1)
    ENDIF
    !
    uf = czero
    w2 = zero
    wf = zero
    !
    mat_out = czero
    !
    CALL fkbounds(nptotf, ip_start, ip_end)
    !
    DO iq = 1, nqtotf ! iq -> q
      xxq = xqf(1:3, iq)
      xxq_r = xxq(1:3)
      mirror_q = .false.
      ! if we need to force the time-rev symmetry, we have to ensure that the phase of uf is fixed
      ! i.e. uf = uf*(-q)
      IF (t_rev) THEN
        IF (is_mirror_q (iq)) THEN
          xxq_r = xqf(1:3, kpg_map(iq))
          mirror_q = .true.
        ENDIF
      ENDIF
      !
      ! Get phonon eigenmode and eigenfrequencies
      IF (.NOT. lifc) THEN
        ! Incompatible bugs found 9/4/2020 originated from the latest EPW changes.
        ! parallel q is not working any more due to mp_sum in rgd_blk
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, mirror_q)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, mirror_q)
      ENDIF
      !
      DO inu = 1, nmodes
        IF (w2(inu) > -eps8) THEN
          wf(inu, iq) =  DSQRT(ABS(w2(inu)))
        ELSE
          wf(inu, iq) = 0.d0
        ENDIF
      ENDDO
      !
      IF (PRESENT(ip_center)) THEN
        center_shift(1:3) = nqf_p / 2 - ip_center
      ELSE
        center_shift(1:3) = 0
      ENDIF
      ! For mirror q, calculate the time-symmetric q' and get uf from q'
      ! e_{\kappa\alpha\nu}(-q)= e^*_{\kappa\alpha\nu}(q)
      !!IF(t_rev .and. iq > iqpg) uf = CONJG(uf) !transpose
      DO inu = 1, nmodes ! inu -> nu
        IF (wf(inu, iq) < eps_acoustic) CYCLE !JLB - cycle zero and imaginary frequency modes
        IF (PRESENT(acoustic_plrn)) THEN
          IF (wf(inu, iq) > acoustic_plrn) CYCLE 
          ! KL - include only modes with frequency lower than acoustic_plrn
        ENDIF
        DO ika = 1, nmodes ! ika -> kappa alpha
          ina = (ika - 1) / 3 + 1
          ctemp = DSQRT(two / (wf(inu, iq) * amass(ityp(ina))))
          ! Parallel run, only calculate the local cell ip
          ! Note that, ip_end obtained from fkbounds should be included
          ! If you have 19 kpts and 2 pool,
          ! lower_bnd= 1 and upper_bnd=10 for the first pool
          ! lower_bnd= 1 and upper_bnd=9 for the second pool
          DO iRp = ip_start, ip_end !, (nqf1_p + 1)/2
            Rp_vec(1:3) = MODULO(index_Rp(iRp, nqf_p) + center_shift, nqf_p)
            ! D_{\kappa\alpha\nu,p}(q) = e_{\kappa\alpha,\nu}(q) \exp(iq\cdot R_p)
            dtemp = uf(ika, inu) * EXP(CMPLX(0.0_DP, twopi * DOT_PRODUCT(xxq, Rp_vec), KIND = DP))
            IF (itype == 1) THEN ! Bqv -> dtau
              ! \Delta \tau_{\kappa\alpha p} = -\frac{1}{N_p} \sum_{q\nu} C_{\kappa\nu q} D_{\kappa\alpha\nu q}  B^*_{q\nu}
              ! Dtau(iRp, ika) = Dtau(iRp, ika) + conjg(B(iq, inu)) * ctemp * dtemp
              mat_out(iRp, ika) = mat_out(iRp, ika) -  cone / REAL(nptotf, dp) * dtemp * ctemp &
                 * (-type_plrn) * CONJG(mat_in(iq, inu))
            ELSE IF(itype == -1) THEN
              !  B_{q\nu} = \frac{1}{N_p} \sum_{\kappa\alpha p} D_{\kappa \alpha\nu, p}(q) C_{q}\nu \Delta\tau_{\kappa\alpha p}
              mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp / ctemp * CONJG(mat_in(iRp, ika)) !JLB: dtau should be real but just in case
              !mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * mat_in(iRp, ika)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! sum all the cell index ip
    CALL mp_sum(mat_out, inter_pool_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE plrn_bmat_tran
    !-----------------------------------------------------------------------------------
    SUBROUTINE scell_plrn_bmat_tran(ttype, t_rev, mat_in, nqtotf_p, nRp_p, Rp_p, &
          nrr_q, ndegen_q, irvec_q, rws, nrws, mat_out)
    !-----------------------------------------------------------------------------------
    !! JLB: Fourier transform between Bmat and dtau for non-diagonal supercells
    !-----------------------------------------------------------------------------------
    USE global_var,    ONLY : xqf, wf
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : eps8, czero, one, two, twopi, zero, ci, cone
    USE ions_base,     ONLY : amass, ityp
    USE wannier2bloch, ONLY : dynwan2bloch, dynifc2blochf
    USE input,         ONLY : lifc, type_plrn, eps_acoustic
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE parallelism,   ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! Transformation direction, 'Bloch2Wan' or 'Wan2Bloch'
    LOGICAL, INTENT(in) :: t_rev
    !! .true. if time-reversal symmetry is to be imposed in Bqv coefficients
    INTEGER, INTENT(in) :: nqtotf_p
    !! Number of q-points in fine grid
    INTEGER, INTENT(in) :: nRp_p
    !! Number of unit cells within supercell
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! Lattice vector coordinates within supercell
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: nrws
    !! Number of real-space Wigner-Seitz
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! Real-space wigner-Seitz vectors
    COMPLEX(KIND = DP), INTENT(in) :: mat_in(:, :)
    !! Input matrix with Bqv/dtau coefficients
    COMPLEX(KIND = DP), INTENT(out) :: mat_out(:, :)
    !! Output matrix with dtau/Bqv coefficients
    !
    ! Local variables
    LOGICAL :: mirror_q
    !! .true. if q1 is a TR mirror point of another q2 point
    INTEGER :: iq
    !! q-point counter
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: itype
    !! Atom type counter
    INTEGER :: ika
    !! Combined atom and cartesian direction counter
    INTEGER :: ip_start
    !! Initial lattice vector within this pool
    INTEGER :: ip_end
    !! Final lattice vector within this pool
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: ina
    !! Atom counter
    REAL(KIND = DP) :: xxq(3)
    !! q-point coordinate
    REAL(KIND = DP) :: xxq_r(3)
    !! auxiliary q-point coordinate in case TR is to be imposed
    REAL(KIND = DP) :: ctemp
    !! Prefactor
    REAL(KIND = DP) :: w2(nmodes)
    !! Phonon frequency squared
    COMPLEX(KIND = DP) :: dtemp
    !! Prefactor
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! Phonon eigenvectors
    !
    IF (ttype == 'Bmat2Dtau') THEN
      itype =  1
    ELSE IF (ttype == 'Dtau2Bmat') THEN
      itype = -1
    ELSE
      CALL errore('scell_plrn_bmat_tran', 'Illegal translation form; should be Bmat2Dtau or Dtau2Bmat!', 1)
    ENDIF
    !
    uf = czero
    w2 = zero
    wf = zero
    !
    mat_out = czero
    !
    CALL fkbounds(nRp_p, ip_start, ip_end)
    !
    DO iq = 1, nqtotf_p ! iq -> q
      xxq = xqf(1:3, iq)
      xxq_r = xxq(1:3)
      mirror_q = .false.
      ! if we need to force the time-rev symmetry, we have to ensure that the phase of uf is fixed
      ! i.e. uf = uf*(-q)
      IF (t_rev) THEN
        IF (is_mirror_q (iq)) THEN
          xxq_r = xqf(1:3, kpg_map(iq))
          mirror_q = .TRUE.
        ENDIF
      ENDIF
      !
      ! Get phonon eigenmode and eigenfrequencies
      IF (.NOT. lifc) THEN
        ! Incompatible bugs found 9/4/2020 originated from the latest EPW changes.
        ! parallel q is not working any more due to mp_sum in rgd_blk
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, mirror_q)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, mirror_q)
      ENDIF
      !
      DO inu = 1, nmodes
        IF (w2(inu) > -eps8) THEN
          wf(inu, iq) =  DSQRT(ABS(w2(inu)))
        ELSE
          wf(inu, iq) = 0.d0
        ENDIF
      ENDDO
      ! For mirror q, calculate the time-symmetric q' and get uf from q'
      ! e_{\kappa\alpha\nu}(-q)= e^*_{\kappa\alpha\nu}(q)
      !IF(t_rev .and. iq > iqpg) uf = CONJG(uf) !transpose
      DO inu = 1, nmodes ! inu -> nu
        IF (wf(inu, iq) < eps_acoustic) CYCLE !JLB - cycle zero and imaginary frequency modes
        DO ika = 1, nmodes ! ika -> kappa alpha
          ina = (ika - 1) / 3 + 1
          ctemp = DSQRT(two / (wf(inu, iq) * amass(ityp(ina))))
          !
          DO iRp = ip_start, ip_end
            ! D_{\kappa\alpha\nu,p}(q) = e_{\kappa\alpha,\nu}(q) \exp(iq\cdot R_p)
            dtemp = uf(ika, inu) * EXP( twopi * ci * DOT_PRODUCT(xxq(1:3), Rp_p(1:3, iRp)) )
            IF (itype == 1) THEN ! Bqv -> dtau
              ! \Delta \tau_{\kappa\alpha p} = -\frac{1}{N_p} \sum_{q\nu} C_{\kappa\nu q} D_{\kappa\alpha\nu q}  B^*_{q\nu}
              ! Dtau(iRp, ika) = Dtau(iRp, ika) + conjg(B(iq, inu)) * ctemp * dtemp
              mat_out(iRp, ika) = mat_out(iRp, ika) -  cone / REAL(nRp, dp) * dtemp * ctemp &
                 * (-type_plrn) * CONJG(mat_in(iq, inu))
            ELSE IF(itype == -1) THEN
              !  B_{q\nu} = \frac{1}{N_p} \sum_{\kappa\alpha p} D_{\kappa \alpha\nu, p}(q) C_{q}\nu \Delta\tau_{\kappa\alpha p}
              mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp / ctemp * CONJG(mat_in(iRp, ika)) !JLB: dtau should be real but just in case
              !mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * mat_in(iRp, ika)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! sum all the cell index ip
    CALL mp_sum(mat_out, inter_pool_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE scell_plrn_bmat_tran
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE calc_den_of_state(eigvec_coef, bqv_coef)
    !-----------------------------------------------------------------------------------
    !! Compute the DOS for Ank and Bqv coefficients
    !-----------------------------------------------------------------------------------
    USE input,         ONLY : nDOS_plrn, edos_max_plrn, edos_min_plrn, edos_sigma_plrn,   &
                              pdos_max_plrn, pdos_min_plrn, pdos_sigma_plrn
    USE global_var,    ONLY : nqtotf, nktotf, wf
    USE io_var,        ONLY : idosplrn
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : ryd2mev, czero, one, ryd2ev, two, zero, cone, pi, ci, twopi,&
                              eps6, eps8, eps5
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(in) :: eigvec_coef(:, :)
    !! Polaron wave function coefficients in the Bloch basis, Ank
    COMPLEX(KIND = DP), INTENT(in) :: bqv_coef(:, :)
    !! Polaron displacement coefficients in the phonon basis, Bqv
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    INTEGER :: idos
    !! Counter for DOS
    INTEGER :: iq
    !! q-point counter
    INTEGER :: inu
    !! Phonon mode counter
    INTEGER :: ik
    !! k-point counter
    INTEGER :: iplrn
    !! Polaron state counter
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    REAL(KIND = DP) :: temp
    !! Temporary max or min eigenvalue
    REAL(KIND = DP), ALLOCATABLE :: f_tmp(:)
    !! Temporary Gaussian function
    REAL(KIND = DP), ALLOCATABLE :: edos(:)
    !! Ank DOS
    REAL(KIND = DP), ALLOCATABLE :: pdos(:)
    !! Bqv DOS
    REAL(KIND = DP), ALLOCATABLE :: edos_all(:)
    !! Electron DOS
    REAL(KIND = DP), ALLOCATABLE :: pdos_all(:)
    !! Phonon DOS
    REAL(KIND = DP), ALLOCATABLE :: e_grid(:)
    !! Grid of energy points to compute Ank DOS
    REAL(KIND = DP), ALLOCATABLE :: p_grid(:)
    !! Grid of energy points to compute Bqv DOS
    !
    !Calculating DOS
    ALLOCATE(f_tmp(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating f_tmp', 1)
    ALLOCATE(e_grid(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating e_grid', 1)
    ALLOCATE(edos(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating edos', 1)
    ALLOCATE(edos_all(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating edos_all', 1)
    temp = MAXVAL(etf_all) * ryd2ev + one
    IF (edos_max_plrn < temp) edos_max_plrn = temp
    temp = MINVAL(etf_all) * ryd2ev - one
    IF (edos_min_plrn > temp) edos_min_plrn = temp
    !
    e_grid = zero
    DO idos = 1, nDOS_plrn
      e_grid(idos) = edos_min_plrn + idos * (edos_max_plrn - edos_min_plrn) / (nDOS_plrn)
    ENDDO
    !
    edos = zero
    edos_all = zero
    DO ik = 1, nktotf
      DO ibnd = 1, nbnd_plrn
        ! TODO : iplrn
        DO iplrn = 1, 1
          CALL cal_f_delta(e_grid - (etf_all(select_bands_plrn(ibnd), ik) * ryd2ev), &
             edos_sigma_plrn, f_tmp)
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          edos = edos + (ABS(eigvec_coef(indexkn1, iplrn))**2) * f_tmp
          edos_all = edos_all + f_tmp
        ENDDO
      ENDDO
    ENDDO
    !
    ALLOCATE(p_grid(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating p_grid', 1)
    ALLOCATE(pdos(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos', 1)
    ALLOCATE(pdos_all(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos_all', 1)
    !
    temp = MAXVAL(wf) * ryd2mev + 10.0_dp
    IF (pdos_max_plrn < temp) pdos_max_plrn = temp
    temp = MINVAL(wf) * ryd2mev - 10.0_dp
    IF (pdos_min_plrn > temp) pdos_min_plrn = temp
    !
    p_grid = zero
    DO idos = 1, nDOS_plrn
      p_grid(idos) = pdos_min_plrn + idos * (pdos_max_plrn - pdos_min_plrn) / (nDOS_plrn)
    ENDDO
    !
    pdos = zero
    pdos_all = zero
    DO iq = 1, nqtotf
      DO inu = 1, nmodes
        CALL cal_f_delta(p_grid - wf(inu, iq) * ryd2mev, pdos_sigma_plrn, f_tmp)
        pdos = pdos + (ABS(bqv_coef(iq, inu))**2) * f_tmp
        pdos_all = pdos_all + f_tmp
      ENDDO
    ENDDO
    !
    OPEN(UNIT = idosplrn, FILE = 'dos.plrn')
    WRITE(idosplrn, '(/2x, a/)') '#energy(ev)  A^2   edos  energy(mev)  B^2  pdos'
    DO idos = 1, nDOS_plrn
      WRITE(idosplrn, '(6f15.7)') e_grid(idos), edos(idos), &
         edos_all(idos), p_grid(idos), pdos(idos), pdos_all(idos)
    ENDDO
    CLOSE(idosplrn)
    !
    DEALLOCATE(pdos_all)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos_all', 1)
    DEALLOCATE(f_tmp)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating f_tmp', 1)
    DEALLOCATE(e_grid)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating e_grid', 1)
    DEALLOCATE(edos)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating edos', 1)
    DEALLOCATE(edos_all)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating edos_all', 1)
    DEALLOCATE(p_grid)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating p_grid', 1)
    DEALLOCATE(pdos)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE calc_den_of_state
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE write_real_space_wavefunction()
    !-----------------------------------------------------------------------------------
    !! Write polaron wave function in real space to file.
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, czero, cone, twopi, ci, bohr2ang
    USE input,         ONLY : nbndsub, step_wf_grid_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE io_var,        ONLY : ipsirplrn
    USE mp_world,      ONLY : world_comm
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE parallelism,   ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 60) :: plrn_file
    !! Name of output file
    INTEGER :: ierr
    !! Error status
    INTEGER :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER :: nkf3_p
    !! Fine k-point grid along b3
    INTEGER :: nktotf_p
    !! Number of k-points in fine grid
    INTEGER :: nbndsub_p
    !! Number of bands in polaron wave function expansion
    INTEGER :: nplrn_p
    !! Number of polaron states
    INTEGER :: nqf_p(3)
    !! Fine q-point grid
    INTEGER :: nqtotf_p
    !! Number of k-points in fine grid
    INTEGER :: nmodes_p
    !! Number of phonon modes
    INTEGER :: ibnd
    !! Electron band index
    INTEGER :: idir
    !! Cartesian direction counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: nxx, nyy, nzz
    !! Number of grid points in real space supercell along cartesian directions
    INTEGER :: ip_min
    !! Initial lattice vector within this pool
    INTEGER :: ip_max
    !! Final lattice vector within this pool
    INTEGER :: ig_vec(1:3)
    !! Supercell latice vector coordinates
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: n_grid(3)
    !! Number of grid points in cell where Wannier functions are written
    INTEGER :: grid_start(3)
    !! Initial grid point within this pool
    INTEGER :: grid_end(3)
    !! Final grid point within this pool
    INTEGER :: n_grid_super(3)
    !! Number of grid points in supercellcell where polaron wave function is written
    INTEGER :: r_grid_vec(3)
    !! Grid point coordinates
    INTEGER :: rpc(1:3)
    !! Grid point coordinates for lattice vectors
    INTEGER :: species(50)
    !! Atomic species counter
    INTEGER :: Rp_vec(1:3)
    !! Lattice vector coordinates
    INTEGER :: shift(1:3)
    !! Shift vector coordinates
    INTEGER :: ishift
    !! Shift counter
    REAL(KIND = DP) :: orig(3)
    !! Supercell origin coordinates
    REAL(KIND = DP) :: cell(3, 3)
    !! Supercell lattice vectors
    REAL(KIND = DP) :: r_cry(3)
    !! Polaron center in crystal coords
    REAL(KIND = DP) :: r_cart(3)
    !! Polaron center in cart coord
    REAL(KIND = DP), ALLOCATABLE :: wann_func(:, :, :, :)
    !! Wannier function in real space
    COMPLEX(KIND = DP) :: ctemp(1:3)
    !! Prefactor for polaron center calculation
    COMPLEX(KIND = DP) :: b_vec(1:3)
    !! Prefactor for polaron center calculation
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! Polaron wave function coefficients in Wannier basis, Amp
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! Polaron displacements in real space
    COMPLEX(KIND = DP), ALLOCATABLE :: cvec(:)
    !! Polaron wave function magnitude squared at each grid point
    !
    ! read Amp.plrn, save eigvec_wan for the latter use
    IF (ionode) THEN
      CALL read_plrn_wf_grid(nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(nkf1_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf2_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf3_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nktotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nbndsub_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nplrn_p,  meta_ionode_id, world_comm)
    !
    ALLOCATE(eigvec_wan(nktotf_p * nbndsub_p, nplrn_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_wf', 'Error allocating eigvec_wan', 1)  
    !
    IF (ionode) THEN
      CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(eigvec_wan, meta_ionode_id, world_comm)
    !
    ! read dtau.plrn, get the displacement.
    IF (ionode) THEN
      CALL read_plrn_dtau_grid(nqf_p(1), nqf_p(2), nqf_p(3), nqtotf_p, nmodes_p, 'dtau.plrn')
    END IF
    CALL mp_bcast(nqf_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nqtotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nmodes_p, meta_ionode_id, world_comm)
    !
    ALLOCATE(dtau(nqtotf_p, nmodes_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating dtau', 1)     
    !
    IF (ionode) THEN
      CALL read_plrn_dtau(dtau, nqtotf_p, nmodes_p, 'dtau.plrn')
    END IF
    CALL mp_bcast(dtau, meta_ionode_id, world_comm)
    !
    ! read cube files for the real-space Wannier function Wm(r)
    CALL read_wannier_cube(select_bands_plrn, wann_func, species, &
       n_grid, grid_start, grid_end)
    !
    cell(1:3, 1) = at(1:3, 1) * nqf_p(1) * alat
    cell(1:3, 2) = at(1:3, 2) * nqf_p(2) * alat
    cell(1:3, 3) = at(1:3, 3) * nqf_p(3) * alat
    !
    plrn_file = 'psir_plrn.xsf'
    ! Write the file head including information of structures
    IF (ionode) THEN
      CALL write_plrn_dtau_xsf(dtau, nqf_p(1), nqf_p(2), nqf_p(3), plrn_file, species)
    END IF
    !
    orig(1:3) = zero
    n_grid_super(1:3) = nqf_p(1:3) * n_grid(1:3)
    !
    IF (ionode) THEN
      OPEN(UNIT = ipsirplrn, FILE = TRIM(plrn_file), POSITION='APPEND')
      WRITE(ipsirplrn, '(/)')
      WRITE(ipsirplrn, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
      WRITE(ipsirplrn, '(3i6)')  n_grid_super / step_wf_grid_plrn
      WRITE(ipsirplrn, '(3f12.6)') zero, zero, zero
      WRITE(ipsirplrn, '(3f12.7)') cell(1:3, 1) * bohr2ang
      WRITE(ipsirplrn, '(3f12.7)') cell(1:3, 2) * bohr2ang
      WRITE(ipsirplrn, '(3f12.7)') cell(1:3, 3) * bohr2ang
    ENDIF
    !
    b_vec(1:3) = twopi * ci / REAL(n_grid_super(1:3))
    !
    ALLOCATE(cvec(1:n_grid_super(1)), STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating cvec', 1)
    !
    CALL fkbounds(nqtotf_p, ip_min, ip_max)
    !
    ctemp = czero
    DO nzz = 1, n_grid_super(3), step_wf_grid_plrn
      DO nyy = 1, n_grid_super(2), step_wf_grid_plrn
        cvec = czero
        DO nxx = 1, n_grid_super(1), step_wf_grid_plrn
          DO iRp = ip_min, ip_max !-nqf_p(3)/2, (nqf_p(3)+1)/2
            Rp_vec(1:3) = index_Rp(iRp, nqf_p)
            rpc(1:3) = Rp_vec(1:3) * n_grid(1:3) !- (n_grid_super(1:3))/2
            ! To make sure that all the nonzero points are included,
            ! we need to try from -1 to 1 neighbor supercells
            DO ishift = 1, 27
              shift(1:3) = index_shift(ishift)
              ig_vec(1:3) = (/nxx, nyy, nzz/) - rpc(1:3) + shift(1:3) * n_grid_super(1:3)
              IF (ALL(ig_vec(1:3) <= grid_end(1:3)) .AND. &
                ALL(ig_vec(1:3) >= grid_start(1:3))) THEN
                DO ibnd = 1, nbndsub !TODO change to nbndsub
                  indexkn1 = (iRp - 1) * nbndsub + ibnd
                  !TODO eigvec_wan(indexkn1, 1) should be eigvec_wan(indexkn1, iplrn)
                  cvec(nxx) = cvec(nxx)  +  eigvec_wan(indexkn1, 1) * wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd)
                ENDDO !ibnd
              ENDIF
            ENDDO ! iscx
          ENDDO ! ipx
        ENDDO ! nxx
        CALL mp_sum(cvec, inter_pool_comm)
        IF(ionode) THEN
          !JLB: Changed to |\Psi(r)|^{2}, I think it's physically more meaningful
          WRITE (ipsirplrn, '(5e13.5)', ADVANCE='yes') ABS(cvec(::step_wf_grid_plrn))**2
        ENDIF
        ! Calculate the center of polaron
        ! TODO: not parallel, all the processors are doing the same calculations
        DO nxx = 1, n_grid_super(1)
          r_grid_vec(1:3) = (/nxx - 1, nyy - 1, nzz - 1/)
          ctemp(1:3) = ctemp(1:3) + EXP(b_vec(1:3) * r_grid_vec(1:3)) * (cvec(nxx)**2)
        ENDDO
        ! End calculating the center of polaron
      ENDDO
    ENDDO
    !
    r_cry(1:3) = IMAG(LOG(ctemp(1:3)))/twopi
    ! make crystal coordinates with 0 to 1
    r_cry(1:3) = r_cry - FLOOR(r_cry)
    r_cart(1:3) = zero
    DO idir = 1, 3
      r_cart(1:3) = r_cart(1:3) + r_cry(idir) * cell(1:3, idir)
    ENDDO
    !
    IF (ionode) THEN
      WRITE(stdout, "(5x, 'The position of polaron:')")
      WRITE(stdout, "(5x, 3f9.4, ' in crystal coordinates')") r_cry(1:3)
      WRITE(stdout, "(5x, 3f9.4, ' in Cartesian coordinates (Angstrom)')") r_cart(1:3)
      WRITE(ipsirplrn, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
      CLOSE(ipsirplrn)
      WRITE(stdout, "(5x, '|\Psi(r)|^2 written to file.')")
    ENDIF
    DEALLOCATE(dtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating dtau', 1)
    DEALLOCATE(eigvec_wan , STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(wann_func, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(cvec , STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE write_real_space_wavefunction
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE scell_write_real_space_wavefunction()
    !-----------------------------------------------------------------------------------
    !!JLB: write psir in transformed supercell
    !!     xsf format no longer compatible,
    !!     as .cube files are written in primitive coords.
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, czero, cone, twopi, ci, bohr2ang
    USE input,         ONLY : nbndsub, step_wf_grid_plrn, as
    USE io_var,        ONLY : iunpsirscell
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE parallelism,   ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm
    USE low_lvl,       ONLY : matinv3
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 60) :: plrn_file
    !! FIXME
    INTEGER :: ierr
    !! Error status
    INTEGER :: nktotf_p
    !! Number of k-points in fine grid
    INTEGER :: nkf1_p
    !! Fine k-point grid along b1
    INTEGER :: nkf2_p
    !! Fine k-point grid along b2
    INTEGER :: nkf3_p
    !! Fine k-point grid along b3
    INTEGER :: nbndsub_p
    !! Number of bands in polaron wave function expansion
    INTEGER :: nplrn_p
    !! Number of polaron states
    INTEGER :: species(50)
    !! Atomic species in unit cell
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: indexkn1
    !! Combined band and k-point index
    INTEGER :: ig_vec(1:3)
    !! Supercell latice vector coordinates
    INTEGER :: ir1, ir2, ir3
    !! Real space grid counters
    INTEGER :: iRp1
    !! Lattice vector counter
    INTEGER :: iRp2
    !! Lattice vector counter
    INTEGER :: n_grid_total
    !! Number of points in real space grid within supercell
    INTEGER :: ip_min
    !! First lattice vector within this pool 
    INTEGER :: ip_max
    !! Last lattice vector within this pool
    INTEGER :: n_grid(3)
    !! Number of grid points in cell where Wannier functions are written
    INTEGER :: grid_start(3)
    !! Initial grid point within this pool
    INTEGER :: grid_end(3)
    !! Final grid point within this pool
    INTEGER :: r_in_crys_p_sup(3)
    !! Real space grid point indices in supercell
    INTEGER :: ishift
    !! Shift counter
    INTEGER :: shift(3)
    !! Shift vector coordinates
    REAL(KIND = DP) :: r_in_crys_p(3)
    !! Unit cell grid points in crystal coordinates
    REAL(KIND = DP) :: r_in_crys_s(3)
    !! Supercell grid points in crystal coordinates
    REAL(KIND = DP) :: r_in_cart(3)
    !! Supercell grid points in cartesian coordinates
    REAL(KIND = DP) :: p2s(3, 3)
    !! Primitive to supercell coordinate transformation
    REAL(KIND = DP) :: s2p(3,3)
    !! Supercell to primitive coordinate transformation
    REAL(KIND = DP), ALLOCATABLE :: wann_func(:, :, :, :)
    !! Wannier functions in real space grid
    COMPLEX(KIND = DP) :: cvec
    !! Polaron wave function modulus squared
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! Polaron wave function coefficient in Wannier basis, Amp
    !
    ! Broadcast supercell lattice vectors
    CALL mp_bcast(as, meta_ionode_id, world_comm)
    !
    ! read Amp.plrn, save eigvec_wan
    IF (ionode) THEN
      CALL read_plrn_wf_grid(nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(nkf1_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf2_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf3_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nktotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nbndsub_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nplrn_p,  meta_ionode_id, world_comm)
    !
    ALLOCATE(eigvec_wan(nktotf_p * nbndsub_p, nplrn_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_wf', 'Error allocating eigvec_wan', 1)  
    !
    IF (ionode) THEN
      CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nplrn_p, 'Amp.plrn')
    END IF
    CALL mp_bcast(eigvec_wan, meta_ionode_id, world_comm)
    !
    ! read cube files for the real-space Wannier function Wm(r)
    CALL read_wannier_cube(select_bands_plrn, wann_func, species, &
       n_grid, grid_start, grid_end)
    !
    ! Read list of Rp-s within supercell
    CALL read_Rp_in_S()
    WRITE(stdout, '(a, i12)') "     Number of unit cells within supercell:", nRp
    !
    ! Open file
    plrn_file = 'psir_plrn.scell.csv'
    IF (ionode) THEN
      OPEN(UNIT = iunpsirscell, FILE = TRIM(plrn_file), FORM = 'formatted', STATUS = 'unknown')
      WRITE(iunpsirscell, '(a)') "x , y , z, |\psi(r)|^2"
    ENDIF
    !
    ! Total number of grid points
    n_grid_total = nRp * n_grid(1) * n_grid(2) * n_grid(3)
    WRITE(stdout, '(a, i12)') "     Total grid points:", n_grid_total
    WRITE(stdout, '(a, i12)') "     Step:", step_wf_grid_plrn
    !
    ! Parallelize iRp
    CALL fkbounds(nRp, ip_min, ip_max)
    !
    ! Matrix to transform from primitive to supercell crystal coordinates
    p2s = matinv3(TRANSPOSE(as))
    p2s = MATMUL(p2s, at)
    ! Supercell to primitive coordinates
    s2p = matinv3(p2s)
    !
    ! Loop over all the grid points
    DO iRp1 = 1, nRp
      !
      DO ir1 = 1, n_grid(1), step_wf_grid_plrn
        DO ir2 = 1, n_grid(2), step_wf_grid_plrn
          DO ir3 = 1, n_grid(3), step_wf_grid_plrn
            !
            r_in_crys_p(1:3) = (/REAL(ir1 - 1, DP) / n_grid(1), REAL(ir2 - 1, DP) / n_grid(2), REAL(ir3 - 1, DP) / n_grid(3)/) &
                               + REAL(Rp(1:3, iRp1), DP)
            !
            ! Wannier functions stored in (1:ngrid*iRp) list
            r_in_crys_p_sup(1:3) = (/ir1, ir2, ir3/) +  Rp(1:3, iRp1) * n_grid(1:3)
            !
            ! Move the r-point to the first supercell and store in cartesian coordinates for plotting
            r_in_crys_s = MATMUL(p2s, r_in_crys_p)
            r_in_crys_s = MODULO(r_in_crys_s, (/1.d0, 1.d0, 1.d0/))
            r_in_cart   = MATMUL(TRANSPOSE(as), r_in_crys_s) * alat * bohr2ang
            !
            ! Sum over p, PRB 99, 235139 Eq.(47)
            cvec = czero
            DO iRp2 = ip_min, ip_max !1, nRp
              !
              DO ishift = 1, 27
                !
                shift(1:3) = index_shift(ishift)
                ig_vec(1:3) = r_in_crys_p_sup(1:3) - Rp(1:3, iRp2) * n_grid(1:3) + MATMUL(s2p, shift(1:3)) * n_grid(1:3)
                !
                IF (ALL(ig_vec(1:3) <= grid_end(1:3)) .AND. &
                  ALL(ig_vec(1:3) >= grid_start(1:3))) THEN
                  !
                  DO ibnd = 1, nbndsub ! sum over m
                    !
                    indexkn1 = (iRp2 - 1) * nbndsub + ibnd
                    cvec = cvec + eigvec_wan(indexkn1, 1) * wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd)
                    !
                  ENDDO !ibnd
                  !
                ENDIF
                !
              ENDDO ! ishift
              !
            ENDDO ! iRp2
            CALL mp_sum(cvec, inter_pool_comm)
            !
            ! Write |\psi(r)|^2 data point to file
            IF (ionode) THEN
              WRITE(iunpsirscell, '(f12.6,", ", f12.6,", ", f12.6,", ", E13.5)') r_in_cart(1:3), ABS(cvec)**2
            ENDIF
            !
          ENDDO ! ir3
        ENDDO ! ir2
      ENDDO ! ir1
    ENDDO ! iRp1
    !
    IF (ionode) THEN
      CLOSE(iunpsirscell)
      WRITE(stdout, "(5x, '|\Psi(r)|^2 written to file.')")
    ENDIF
    DEALLOCATE(wann_func, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(eigvec_wan, STAT = ierr)
    IF (ierr /= 0) CALL errore('scell_write_real_space_wavefunction', 'Error allocating eigvec_wan', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE scell_write_real_space_wavefunction
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE read_wannier_cube(select_bands, wann_func, species, n_grid, &
               grid_start_min, grid_end_max)
    !-----------------------------------------------------------------------------------
    !! Read the nth Wannier function from prefix_0000n.cube file
    !-----------------------------------------------------------------------------------
    USE ep_constants,  ONLY : zero, czero, cone
    USE io_var,        ONLY : iun_plot
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : ionode, meta_ionode_id
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE parallelism,   ONLY : fkbounds
    USE input,         ONLY : nbndsub
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: select_bands(:)
    !! Wannier functions in which polaron wave function has been expanded
    INTEGER, INTENT(out) :: species(50)
    !! Atomic species
    INTEGER, INTENT(out) :: n_grid(3)
    !! Number of points in real space grid where Wannier functions are written
    INTEGER, INTENT(out) :: grid_start_min(3)
    !! Initial grid point within this pool
    INTEGER, INTENT(out) :: grid_end_max(3)
    !! Final grid point within this pool
    REAL(KIND = DP), ALLOCATABLE, INTENT(out) :: wann_func(:, :, :, :)
    !! Wannier function in real space
    !
    ! Local variables
    CHARACTER(LEN = 60) :: wancube
    !! Name of file containing Wannier function
    CHARACTER(LEN = 60) :: temp_str
    !! Temporary string
    INTEGER :: ierr
    !! Error status
    INTEGER :: ibnd
    !! Electron band counter
    INTEGER :: ie
    !! Atomic species counter
    INTEGER :: idir
    !! Cartesian direction counter
    INTEGER :: i_species
    !! Atomic species index
    INTEGER :: nbnd
    !! Number of Wannier funcions in which polaron wave function is expanded
    INTEGER :: iline
    !! Counter along line
    INTEGER :: nAtoms
    !! Total number of atoms
    INTEGER :: nxx, nyy, nzz
    !! Number of grid points in Cartesian directions
    INTEGER :: n_len_z
    !! Number of grid points within this pool
    INTEGER :: grid_start(3)
    !! Initial grid point within this loop
    INTEGER :: grid_end(3)
    !! Final grid point within this loop
    REAL(KIND = DP) :: rtempvec(4)
    !! Temporary vector to be read from .cube file
    REAL(KIND = DP) :: norm
    !! Wannier function normalization
    !
    nbnd = SIZE(select_bands)
    ! find the max and min of real space grid of Wannier functions of all Wannier orbitals
    IF (ionode) THEN
      grid_start_min(:) = 100000
      grid_end_max(:) = -100000
      DO ibnd = 1, nbndsub
        WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
        OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
        READ(iun_plot, *) temp_str !, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str
        READ(iun_plot, *) n_grid, grid_start, grid_end
        DO idir = 1, 3
          IF (grid_start_min(idir) >= grid_start(idir)) grid_start_min(idir) = grid_start(idir)
          IF (grid_end_max(idir) <= grid_end(idir))   grid_end_max(idir) = grid_end(idir)
        ENDDO
        CLOSE(iun_plot)
      ENDDO
    ENDIF
    !
    CALL mp_bcast(n_grid,         meta_ionode_id, world_comm)
    CALL mp_bcast(grid_start_min, meta_ionode_id, world_comm)
    CALL mp_bcast(grid_end_max,   meta_ionode_id, world_comm)
    !
    ! Read the xth Wannier functions from prefix_0000x.cube in ionode
    ! and broadcast to all nodes
    ALLOCATE(wann_func(grid_start_min(1):grid_end_max(1), &
       grid_start_min(2):grid_end_max(2), &
       grid_start_min(3):grid_end_max(3), nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_wannier_cube', 'Error allocating wann_func', 1)
    wann_func = zero
    species = 0
    IF (ionode) THEN
      DO ibnd = 1, nbndsub
        WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
        OPEN(UNIT = iun_plot, FILE = TRIM(wancube), FORM = 'formatted', STATUS = 'unknown')
        READ(iun_plot, *) temp_str
        READ(iun_plot, *) n_grid, grid_start, grid_end
        READ(iun_plot, *) nAtoms, rtempvec(1:3)
        !
        DO iline = 1, 3
          READ(iun_plot, '(8A)') temp_str
        ENDDO
        ie = 1
        DO iline = 1, nAtoms
          READ(iun_plot, '(i4, 4f13.5)') i_species, rtempvec
          IF (iline == 1 ) THEN
            species(ie) = i_species
            ie = ie + 1
          ELSE IF (species(ie - 1) /= i_species) THEN
            species(ie) = i_species
            ie = ie + 1
          ENDIF
        ENDDO
        n_len_z = grid_end(3) - grid_start(3) + 1
        !
        DO nxx = grid_start(1), grid_end(1)
          DO nyy = grid_start(2), grid_end(2)
            DO nzz = grid_start(3), grid_end(3), 6
              IF (grid_end(3) - nzz < 6) THEN
                READ(iun_plot, *) wann_func(nxx, nyy, nzz:grid_end(3) - 1, ibnd)
              ELSE
                READ(iun_plot, '(6E13.5)') wann_func(nxx, nyy, nzz:nzz + 5, ibnd)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iun_plot)
        ! Wannier function is not well normalized
        ! Normalize here will make the calculations with Wannier functions easier
        norm = SUM(wann_func(:, :, :, ibnd) * wann_func(:, :, :, ibnd))
        wann_func(:, :, :, ibnd) = wann_func(:, :, :, ibnd) / SQRT(norm)
      ENDDO
    ENDIF
    CALL mp_bcast(wann_func, meta_ionode_id, world_comm)
    CALL mp_bcast(species, meta_ionode_id, world_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_wannier_cube
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE read_Rp_in_S()
    !-----------------------------------------------------------------------------------
    ! JLB
    !! Allocate and read list of Rp unit cell vectors contained on transformed supercell
    !-----------------------------------------------------------------------------------
    USE io_var,    ONLY : iunRpscell
    USE io_global, ONLY : ionode, meta_ionode_id
    USE mp,        ONLY : mp_bcast
    USE mp_world,  ONLY : world_comm
    USE global_var,ONLY : nqtotf
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: iRp
    !! Lattice vector counter
    INTEGER :: ierr
    !! Error status
    INTEGER :: nRp2
    !! Number of lattice vectors within supercell
    !
    IF (ionode) THEN
      OPEN(UNIT = iunRpscell, FILE = 'Rp.scell.plrn', FORM = 'formatted', STATUS = 'unknown')
      READ(iunRpscell, *) nRp
      IF (nRp /= nqtotf) CALL errore('read_Rp_in_S', 'nRp and nqtotf are not the same!',1)
      CLOSE(UNIT = iunRpscell)
    ENDIF
    CALL mp_bcast(nRp, meta_ionode_id, world_comm)
    ALLOCATE(Rp(3, nRp), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_Rp_in_S', 'Error allocating Rp', 1)
    Rp = 0
    IF (ionode) THEN
      OPEN(UNIT = iunRpscell, FILE = 'Rp.scell.plrn', FORM = 'formatted', STATUS = 'unknown')
      READ(iunRpscell, *) nRp2
      DO iRp = 1, nRp
        READ(iunRpscell, *) Rp(1:3, iRp)
      ENDDO
      CLOSE(UNIT = iunRpscell)
    ENDIF
    CALL mp_bcast(Rp, meta_ionode_id, world_comm)
    !
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_Rp_in_S
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    FUNCTION index_Rp(iRp, nqfs)
    !-----------------------------------------------------------------------------------
    !! Index
    !-----------------------------------------------------------------------------------
    USE input,      ONLY: nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iRp
    !! Lattice vector index
    INTEGER, INTENT(in), OPTIONAL  :: nqfs(1:3)
    !! k/q points along each direcion in grid
    !
    ! Local variable
    INTEGER  :: index_Rp(1:3)
    !! Lattice vector in crystal coords
    INTEGER  :: nqf_c(1:3)
    !! k/q points along each direcion in grid
    !
    IF (PRESENT(nqfs)) THEN
      nqf_c(1:3) = nqfs(1:3)
    ELSE
      nqf_c(1:3) = (/nqf1, nqf2, nqf3/)
    ENDIF
    !
    index_Rp(1) = (iRp - 1)/(nqf_c(2) * nqf_c(3))
    index_Rp(2) = MOD(iRp - 1, nqf_c(2) * nqf_c(3))/nqf_c(3)
    index_Rp(3) = MOD(iRp - 1, nqf_c(3))
    !
    IF (ANY(index_Rp < 0) .OR. ANY(index_Rp >= nqf_c)) THEN
      CALL errore('index_Rp','index_Rp not correct!',1)
    ENDIF
    !-----------------------------------------------------------------------------------
    END FUNCTION index_Rp
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    FUNCTION index_shift(ishift)
    !-----------------------------------------------------------------------------------
    !! Find supercell lattice vector in crystal coords for shift loop around neighbors
    !-----------------------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ishift
    !! Shift vector index
    !
    ! Local variable
    INTEGER  :: index_shift(1:3)
    !! Shift vector crystal coords.
    !
    index_shift(1) = (ishift - 1) / 9 - 1
    index_shift(2) = MOD(ishift - 1, 9) / 3 - 1
    index_shift(3) = MOD(ishift - 1, 3) - 1
    !
    IF (ANY(index_shift < -1) .OR. ANY(index_shift > 1)) THEN
      CALL errore('index_shift', 'index_shift not correct!', 1)
    ENDIF
    !-----------------------------------------------------------------------------------
    END FUNCTION index_shift
    !-----------------------------------------------------------------------------------
    !
    SUBROUTINE plrn_close()
    !-----------------------------------------------------------------------------------
    !! Deallocate arrays and close polaron calculations
    !-----------------------------------------------------------------------------------
    USE input,         ONLY : model_vertex_plrn, interp_Ank_plrn, &
                              interp_Bqu_plrn, cal_psir_plrn, scell_mat_plrn
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    IF (.NOT. (interp_Ank_plrn .OR. interp_Bqu_plrn .OR. cal_psir_plrn)) THEN
      DEALLOCATE(is_mirror_k, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating is_mirror_k', 1)
      DEALLOCATE(is_mirror_q, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating is_mirror_q', 1)
      DEALLOCATE(is_tri_k, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating is_tri_k', 1)
      DEALLOCATE(is_tri_q, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating is_tri_q', 1)
      DEALLOCATE(kpg_map, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating kpg_map', 1)
      DEALLOCATE(Hamil, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating Hamil', 1)
    ENDIF
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating Hamil', 1)
    DEALLOCATE(xkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating xkf_all', 1)
    DEALLOCATE(select_bands_plrn, STAT = ierr)
    IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating select_bands_plrn', 1)
    IF (interp_Ank_plrn .OR. (.NOT. (interp_Bqu_plrn .OR. cal_psir_plrn))) THEN
      DEALLOCATE(eigvec, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating eigvec', 1)
    END IF
    IF (scell_mat_plrn) THEN
      DEALLOCATE(Rp, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating Rp', 1)
    END IF
    IF (model_vertex_plrn) THEN
      DEALLOCATE(gq_model, STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_close', 'Error deallocating gq_model', 1)
    END IF
    !
    !-----------------------------------------------------------------------------------    
    END SUBROUTINE plrn_close
    !-----------------------------------------------------------------------------------
  !  
  !-----------------------------------------------------------------------------------
  END MODULE polaron
  !-----------------------------------------------------------------------------------
