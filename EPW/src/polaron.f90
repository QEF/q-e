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
  !! is this local k and global q point is a mirror point, used for time-reversal symmertry
  LOGICAL, ALLOCATABLE :: is_mirror_q(:)
  !! is this local k and global q point is a mirror point, used for time-reversal symmertry
  LOGICAL, ALLOCATABLE :: is_tri_q(:)
  !! is this local k and global q point is a mirror point, used for time-reversal symmertry
  LOGICAL, ALLOCATABLE :: is_tri_k(:)
  !! is this local k and global q point is a mirror point, used for time-reversal symmertry
  INTEGER :: nbnd_plrn
  !! FIXME
  INTEGER :: nbnd_g_plrn
  !! FIXME
  INTEGER :: lword_h
  !! FIXME
  INTEGER :: lword_g
  !! FIXME
  INTEGER :: lword_m
  !! FIXME
  INTEGER :: lword_umn
  !! FIXME
  INTEGER :: lword_ekanu
  !! FIXME
  INTEGER :: io_level_g_plrn
  !! FIXME
  INTEGER :: io_level_h_plrn
  !! FIXME
  INTEGER :: io_level_umn_plrn
  !! FIXME
  INTEGER :: io_level_ekanu_plrn
  !! FIXME
  INTEGER :: hblocksize
  !! FIXME
  INTEGER :: band_pos
  !! FIXME
  INTEGER :: ibvec(-3:3)
  !! FIXME
  INTEGER :: ik_edge
  !! FIXME
  INTEGER :: nRp
  !! Number of unit cells on supercell for non-diagonal supercells
  INTEGER, ALLOCATABLE :: Rp(:,:)
  !! List of unit cell vectors within supercell
  INTEGER, ALLOCATABLE :: select_bands_plrn(:)
  !! FIXME
  INTEGER, ALLOCATABLE :: kpg_map(:)
  !! Number of bands in subgroup used in polaron calculations
  REAL(KIND = DP) :: wq_model
  !! phonon freq in vertex model
  REAL(KIND = DP), ALLOCATABLE :: etf_model(:)
  !! band structure and in vertex model
  REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
  !! Gathered k points coordinates over the pools
  REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
  !! Generate the maps of k->k+q and k->G-k
  COMPLEX(KIND = DP), ALLOCATABLE :: Hamil(:, :)
  !! polaron eigenvector
  COMPLEX(KIND = DP), ALLOCATABLE :: eigVec(:, :)
  !! Gathered eigenvalues over the pools
  COMPLEX(KIND = DP), ALLOCATABLE :: hEigVec(:, :)
  !! Gathered eigenvalues over the pools
  COMPLEX(KIND = DP), ALLOCATABLE :: gq_model(:)
  !! el-ph matrix element in vertex model
  COMPLEX(KIND = DP), ALLOCATABLE :: M_mat(:, :, :, :)
  !! el-ph matrix element in vertex model
  COMPLEX(KIND = DP), ALLOCATABLE :: epf(:, :, :, :)
  !! el-ph matrix element in vertex model
  COMPLEX(KIND = DP), ALLOCATABLE :: epfall(:, :, :, :, :)
  !! el-ph matrix element in vertex model
  COMPLEX(KIND = DP), ALLOCATABLE :: rho_mat(:, :, :, :)
  !! FIXME
  COMPLEX(KIND = DP), ALLOCATABLE :: Bmat(:,:)
  !! FIXME
  COMPLEX(KIND = DP) :: berry_phase(1:3)
  !! FIXME
  PUBLIC  :: plrn_prepare
  !! FIXME
  PUBLIC  :: plrn_flow_select
  !! FIXME
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
    USE epwcom,        ONLY : start_band_plrn, end_band_plrn, nbndsub, nstate_plrn, debug_plrn, &
                              cal_psir_plrn, restart_plrn,  interp_Ank_plrn, interp_Bqu_plrn,   &
                              model_vertex_plrn, full_diagon_plrn, nhblock_plrn, Mmn_plrn, lifc,&
                              g_start_band_plrn, g_end_band_plrn, lrot, lphase, type_plrn,      &
                              g_start_energy_plrn, g_end_energy_plrn, nqf1, nqf2, nqf3,         &
                              model_enband_plrn, model_phfreq_plrn, model_vertex_plrn,          &
                              g_power_order_plrn, io_lvl_plrn, nkf1, nkf2, nkf3, seed_plrn,     &
                              m_eff_plrn, kappa_plrn, omega_LO_plrn, fsthick, etf_mem,          &
                              scell_mat_plrn
    USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf, etf, xkf, xqf, wf, xkq, chw
    USE modes,         ONLY : nmodes
    USE mp_world,      ONLY : mpime
    USE cell_base,     ONLY : bg, at, omega, alat
    USE constants_epw, ONLY : czero, cone, pi, ci, twopi, fpi, eps6, eps8, eps5, zero, ryd2ev
    USE poolgathering, ONLY : poolgather2
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE io_files,      ONLY : check_tempdir
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out)  :: iq_restart
    !! Restart q-point indexx
    INTEGER, INTENT(in)   :: totq
    !! Total q-point in the fsthick
    !
    ! Local variables
    CHARACTER(LEN = 100) :: fmt_mode
    !! printing
    LOGICAL :: debug
    !! Debug flag
    LOGICAL :: plrn_scf
    !! Self-consistent
    LOGICAL :: exst
    !! Checking on the file presence
    LOGICAL :: pfs
    !! FIXME
    INTEGER :: inu
    !! mode index
    INTEGER :: iq
    !! q-point index
    INTEGER :: ik
    !! k-point index
    INTEGER :: ikk
    !! k+q point index
    INTEGER :: jk
    !! FIXME
    INTEGER :: iibnd
    !! band index
    INTEGER :: jjbnd
    !! band index
    INTEGER :: ibnd
    !! band index
    INTEGER :: jbnd
    !! band index
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ierr
    !! Error status
    INTEGER :: icount
    !! FIXME
    INTEGER :: iter
    !! FIXME
    INTEGER :: ix, iy, iz
    !! FIXME
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: idos
    !! FIXME
    INTEGER :: iatm
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: ikGamma
    !! FIXME
    INTEGER :: iqGamma
    !! FIXME
    INTEGER :: io_level
    !! FIXME
    INTEGER :: minNBlock
    !! FIXME
    INTEGER :: ishift
    !! FIXME
    INTEGER :: lword_h_tmp
    !! FIXME
    INTEGER :: lword_g_tmp
    !! FIXME
    INTEGER :: lword_umn_tmp
    !! FIXME
    INTEGER :: lword_ekanu_tmp
    !! FIXME
    INTEGER, PARAMETER :: maxword = HUGE(1)
    !! FIXME
    REAL(KIND = DP) :: xxk(3)
    !! FIXME
    REAL(KIND = DP) :: xkf_cart(3, nkqf)
    !! FIXME
    REAL(KIND = DP) :: xqf_cart(3, nqf)
    !! FIXME
    REAL(KIND = DP) :: efermi
    !! FIXME
    REAL(KIND = DP) :: xkf_len(nkqf)
    !! FIXME
    REAL(KIND = DP) :: klen
    !! FIXME
    REAL(KIND = DP) :: shift(3)
    !! FIXME
    REAL(KIND = DP) :: rfac
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rtmp2(:,:)
    !! FIXME
    !
    CALL start_clock('plrn_prepare')
    !
    IF(etf_mem == 3) THEN
      CALL errore('polaron_prepare', 'Polaron module not working with etf_mem = 3', 1)
    ENDIF
    !
    WRITE(stdout, '(5x,"fsthick not working in polaron module, selecting all the k/q points.")')
    !! type_plrn denotes whether electron polaron (-1) or hole polaron (+1)
    !! Legalize the type_plrn input, in case that the user use an arbitrary number
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
    IF(debug_plrn) THEN
      DO ik = 1, nktotf
        WRITE(stdout, '(5x, 3f15.6)') xkf_all(1:3, ik)
      ENDDO
    ENDIF
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
      ALLOCATE(eigVec(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating eigVec', 1)
      eigVec = czero
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
       ALLOCATE(eigVec(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating eigVec', 1)
       eigVec = czero
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
       ALLOCATE(is_mirror_q(nqf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       ALLOCATE(is_tri_k(nkf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       is_tri_k = .FALSE.
       ALLOCATE(is_tri_q(nqf), STAT = ierr)
       IF (ierr /= 0) CALL errore('plrn_prepare', 'Error allocating is_mirror_q', 1)
       is_tri_q = .FALSE.
       is_mirror_q = .false.
       is_mirror_k = .false.
       !
       WRITE(stdout, '(5x, a)') "Finding the index of -k for each k point."
       ! For two k points k1 and k2, if k1 = G - k2, G is any reciprocal vector
       ! then k2 is the mirror point if the index of k2 is larger than k1
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
    SUBROUTINE plrn_save_g_to_file(iq, epfg, wf)
    !-----------------------------------------------------------------------
    !!
    !! Save el-ph matrix element to file
    !!
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : g_start_band_plrn, g_end_band_plrn,    &
                              nbndsub, g_tol_plrn, io_lvl_plrn,      &
                              g_start_energy_plrn, g_end_energy_plrn
    USE elph2,         ONLY : nkf, nktotf, nqtotf
    USE constants_epw, ONLY : eps8, two, czero
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! q-point index
    REAL(KIND = DP), INTENT(in) :: wf(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(inout) :: epfg(:, :, :, :)
    !! el-ph matrix element
    !
    ! Local variables
    INTEGER :: ik
    !! k-point index
    INTEGER :: ikq
    !! k+q point
    INTEGER :: imode
    !! mode index
    INTEGER :: ibnd
    !! band index
    INTEGER :: jbnd
    !! band index
    INTEGER :: ik_global
    !! global k-point
    INTEGER, SAVE :: g_count_sum
    !! FIXME
    REAL(KIND = DP) :: eig
    !! Eigenvalues
    COMPLEX(KIND = DP) :: ctemp
    !! Temporary variable
    !
    ! In polaron equations, g is not epf but epf/omega
    ! To ensure a Hermitian Hamiltonian, g_{mnu}(k, -q) is calculated as g*_{nmu}(k-q, q)
    DO ik = 1, nkf
      ikq = ikq_all(ik, iq)
      DO imode = 1, nmodes
        IF (wf(imode, iq) > eps8) THEN
          epf(:, :, imode, ik) = &
          epfg(g_start_band_plrn:g_end_band_plrn, g_start_band_plrn:g_end_band_plrn, imode, ik) / DSQRT(two * wf(imode, iq))
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
    USE elph2, ONLY : nktotf
    USE epwcom, ONLY : nkf1, nkf2, nkf3
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
    !! FIXME
    INTEGER :: ikq_all
    !! FIXME
    INTEGER :: index_target(1:3)
    !! FIXME
    INTEGER :: index_kq
    !! FIXME
    INTEGER :: ikq_loop
    !! FIXME
    REAL(KIND = DP) :: xxk(1:3)
    !! FIXME
    REAL(KIND = DP) :: xxk_target(1:3)
    !! FIXME
    !
    CALL start_clock('find_k+q')
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
    CALL stop_clock('find_k+q')
    !
    IF (ikq_all == 0) CALL errore('ikq_all','k + q not found', 1)
    !
    !-----------------------------------------------------------------------
    END FUNCTION ikq_all
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION find_ik(xxk, xkf_all)
    !-----------------------------------------------------------------------
    !!
    !! Find k-point index
    !!
    USE elph2, ONLY : nktotf
    USE epwcom, ONLY : nkf1, nkf2, nkf3
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xxk(1:3)
    !! K-point position per cpu
    REAL(KIND = DP), INTENT(in) :: xkf_all(1:3, 1:nktotf)
    !! global k-points
    !
    ! Local variables
    INTEGER :: ik
    !! k-point index
    INTEGER :: find_ik
    !! k-point index of found k-point
    REAL(KIND = DP) :: xkq(1:3)
    !! k+q point index
    !
    CALL start_clock('find_k')
    !
    find_ik = 0
    DO ik = 1, nktotf
      xkq(1:3) = xkf_all(1:3, ik) - xxk(1:3)
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
    !! FIXME
    !!
    USE epwcom,        ONLY : cal_psir_plrn, restart_plrn,  interp_Ank_plrn, interp_Bqu_plrn, &
                              io_lvl_plrn, scell_mat_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: nrr_k
    !! FIXME
    INTEGER, INTENT (in) :: dims
    !! FIXME
    INTEGER, INTENT (in) :: ndegen_k(:,:,:)
    !! FIXME
    INTEGER, INTENT (in) :: nrr_q
    !! FIXME
    INTEGER, INTENT (in) :: ndegen_q(:,:,:)
    !! FIXME
    INTEGER, INTENT (in) :: irvec_q(3, nrr_q)
    !! FIXME
    INTEGER, INTENT (in) :: nrws
    !! FIXME
    REAL(KIND = DP), INTENT (in) :: irvec_r(3, nrr_k)
    !! FIXME
    REAL(KIND = DP), INTENT (in) :: rws(:, :)
    !! FIXME
    !
    ! Local variable
    LOGICAL :: itsopen
    !! FIXME
    !
    ! Bqu Ank interpolation is not compatible with self-consistency process
    ! Added by Chao Lian for polaron calculations flow select
    ! If postprocess is ON, i.e. Bqu interpolation with saved dtau,
    ! Ank interpolation with saved Amp, and polaron visulation with saved Wannier function cube files,
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
    ! Clean up the allocated arrays, close open files inquire(unit=iepfall, opened=itsopen)
    ! if (itsopen) CLOSE(iepfall)
    !
    ! -- FIXME FIXME FIXME FIXME FIXME ------------------------------------------------------------
    ! SP - IF ALLOCATED is not permitted. One must correctly allocate and deallocate when needed
    IF (ALLOCATED(is_mirror_k))      DEALLOCATE(is_mirror_k)
    IF (ALLOCATED(is_mirror_q))      DEALLOCATE(is_mirror_q)
    IF (ALLOCATED(is_tri_q))         DEALLOCATE(is_tri_q) !JLB
    IF (ALLOCATED(is_tri_k))         DEALLOCATE(is_tri_k) !JLB
    IF (ALLOCATED(Hamil))            DEALLOCATE(Hamil)
    IF (ALLOCATED(eigVec))           DEALLOCATE(eigVec)
    IF (ALLOCATED(hEigVec))          DEALLOCATE(hEigVec)
    IF (ALLOCATED(kpg_map))          DEALLOCATE(kpg_map)
    IF (ALLOCATED(etf_all))          DEALLOCATE(etf_all)
    IF (ALLOCATED(xkf_all))          DEALLOCATE(xkf_all)
    IF (ALLOCATED(select_bands_plrn))DEALLOCATE(select_bands_plrn)
    IF (ALLOCATED(gq_model))         DEALLOCATE(gq_model)
    ! -- FIXME FIXME FIXME FIXME FIXME ------------------------------------------------------------
    !
    WRITE(stdout, '(/5x, "======================== Polaron Timers ===========================")')
    CALL print_clock('main_prln')
    CALL print_clock('find_k+q')
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
    CALL print_clock('HOffDiagTerm')
    CALL print_clock('HdiagTerm')
    CALL print_clock( 'cegterg' )
    CALL print_clock( 'cegterg:init' )
    CALL print_clock( 'cegterg:diag' )
    CALL print_clock( 'cegterg:update' )
    CALL print_clock( 'cegterg:overlap' )
    CALL print_clock( 'cegterg:last' )
    CALL print_clock('ik_l2g')
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
    USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, twopi,           &
                              czero, cone, pi, ci, twopi, eps6, eps8, eps5
    USE epwcom,        ONLY : type_plrn, full_diagon_plrn, lifc, debug_plrn,    &
                              init_sigma_plrn, init_k0_plrn, nstate_plrn,       &
                              conv_thr_plrn, mixing_Plrn, init_plrn, niter_plrn,&
                              restart_plrn, nkf1, nkf2, nkf3, seed_plrn,        &
                              nqf1, nqf2, nqf3, r0_plrn, init_ntau_plrn,        &
                              efermi_read, fermi_energy, nbndsub, as,           &
                              model_vertex_plrn, time_rev_A_plrn, beta_plrn,    &
                              Mmn_plrn, recal_Mmn_plrn, model_vertex_plrn,      &
                              model_enband_plrn, model_phfreq_plrn,             &
                              omega_LO_plrn, kappa_plrn, m_eff_plrn,            &
                              scell_mat_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, nkf, nqf, nqtotf,    &
                              nktotf, xkf, wf, xkq, chw, cu, cuq
    USE mp_global,     ONLY : inter_pool_comm
    USE cell_base,     ONLY : bg, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE poolgathering, ONLY : poolgather2
    USE ions_base,     ONLY : nat
    USE mp_world,      ONLY : mpime, world_comm
    USE epwcom,        ONLY : ethrdg_plrn
    USE io_var,        ONLY : iunRpscell
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! FIXME
    INTEGER, INTENT(in) ::  dims
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! FIXME
    INTEGER, INTENT(in) :: nrr_q
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! FIXME
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! FIXME
    INTEGER, INTENT(in) :: nrws
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! FIXME
    !
    ! local variables
    CHARACTER(LEN = 256) :: filename
    !! FIXME
    CHARACTER(LEN = 256) :: tmpch
    !! FIXME
    LOGICAL :: debug
    !! FIXME
    LOGICAL :: file_exist
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikk
    !! FIXME
    INTEGER :: jk
    !! FIXME
    INTEGER :: iibnd
    !! FIXME
    INTEGER :: jjbnd
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ierr
    !! Error status
    INTEGER :: iRp
    !! FIXME
    INTEGER :: itau
    !! Atom index
    INTEGER :: iter
    !! FIXME
    INTEGER :: icount
    !! FIXME
    INTEGER :: ix, iy, iz
    !! FIXME
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: idos
    !! FIXME
    INTEGER :: iatm
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: nkf1_p
    !! FIXME
    INTEGER :: nkf2_p
    !! FIXME
    INTEGER :: nkf3_p
    !! FIXME
    INTEGER :: nbnd_plrn_p
    !! FIXME
    INTEGER :: nbndsub_p
    !! FIXME
    INTEGER :: nPlrn_p
    !! FIXME
    INTEGER :: nktotf_p
    !! FIXME
    INTEGER :: iqpg
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: dos_file
    !! FIXME
    INTEGER :: wan_func_file
    !! FIXME
    INTEGER :: bloch_func_file
    !! FIXME
    INTEGER :: bmat_file
    !! FIXME
    INTEGER :: dtau_file
    !! FIXME
    INTEGER :: itemp
    !! FIXME
    INTEGER :: jtemp
    !! FIXME
    INTEGER :: ngrid(1:3)
    !! FIXME
    REAL(KIND = DP) :: estmteRt(nstate_plrn)
    !! FIXME
    REAL(KIND = DP) :: eigVal(nstate_plrn)
    !! FIXME
    REAL(KIND = DP) :: esterr
    !! FIXME
    REAL(KIND = DP) :: eb
    !! FIXME
    REAL(KIND = DP) :: EPlrnTot
    !! FIXME
    REAL(KIND = DP) :: EPlrnElec
    !! FIXME
    REAL(KIND = DP) :: EPlrnPhon
    !! FIXME
    REAL(KIND = DP) :: EPlrnBeta
    !! FIXME
    REAL(KIND = DP) :: EPlrnDisp
    !! FIXME
    REAL(KIND = DP) :: xxk(3)
    !! FIXME
    REAL(KIND = DP) :: xxq(3)
    !! FIXME
    REAL(KIND = DP) :: shift(3)
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    REAL(KIND = DP) :: disK
    !! FIXME
    REAL(KIND = DP) :: disK_t
    !! FIXME
    REAL(KIND = DP) :: prefac
    !! FIXME
    REAL(KIND = DP) :: norm
    !! FIXME
    REAL(KIND = DP) :: r_cry(1:3)
    !! FIXME
    REAL(KIND = DP) :: totVal_save
    !! FIXME
    REAL(KIND = DP) :: b_vec(1:3)
    !! FIXME
    REAL(KIND = DP) :: dtau_diff
    !! FIXME
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! FIXME
    COMPLEX(KIND = DP) :: cfac(nrr_k, dims, dims)
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: Bmat_save(:,:)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan_save(:, :)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau_save(:, :)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau_list(:, :, :)
    !! FIXME
    !
    ALLOCATE(dtau(nktotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau', 1)
    ALLOCATE(dtau_save(nktotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau_save', 1)
    dtau = czero
    dtau_save = czero
    b_vec(1:3) = (/one / nqf1, one / nqf2, one / nqf3/)
    debug = debug_plrn
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
    eigVec = czero
    SELECT CASE (init_plrn)
      CASE (1)
        ! If k0 has not been set on input, center gaussian at band edge
        IF (ALL(init_k0_plrn(:) == 1000.d0)) init_k0_plrn = xkf_all(1:3, ik_edge)
        !
        WRITE(stdout, '(5x, "Initializing polaron wavefunction using Gaussian wave &
           &packet with a width of", ES14.6)') init_sigma_plrn
        WRITE(stdout, '(5x, "centered at k=", 3f14.6)') init_k0_plrn !xkf_all(1:3, ik_edge)
        CALL init_plrn_gaussian((/zero, zero, zero/), xkf_all, init_k0_plrn, eigVec)
      CASE (3)
        ALLOCATE(eigvec_wan(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
        IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
        WRITE(stdout, '(5x, a)') "Initializing the polaron wavefunction with previously saved Amp.plrn file"
        CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, 'Amp.plrn', scell_mat_plrn)
        CALL plrn_eigvec_tran('Wan2Bloch', time_rev_A_plrn, eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
           nrr_k, ndegen_k, irvec_r, dims, eigVec)
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
          CALL read_plrn_dtau(dtau, nqf1, nqf2, nqf3, nqtotf, nmodes, filename, scell_mat_plrn)
          dtau_list(1, :, :) = dtau(:, :)
        ELSE
          DO itau = 1, init_ntau_plrn
            WRITE(tmpch,'(I4)') itau
            filename = TRIM('dtau_disp.plrn_'//ADJUSTL(tmpch))
            CALL read_plrn_dtau(dtau, nqf1, nqf2, nqf3, nqtotf, nmodes, filename, scell_mat_plrn)
            dtau_list(itau, :, :) = dtau(:, :)
          ENDDO
        ENDIF
        !
        CALL mp_bcast(dtau_list, meta_ionode_id, world_comm)
        ! Initialize Ank wavefunction for iterative diagonalization Gaussian
        ! If k0 has not been set on input, center gaussian at band edge
        IF (ALL(init_k0_plrn(:) == 1000.d0)) init_k0_plrn = xkf_all(1:3, ik_edge)
        WRITE(stdout, '(5x, "Initializing polaron wavefunction using Gaussian wave &
        &packet with the width of", f15.7)') init_sigma_plrn
        WRITE(stdout, '(5x, "centered at k=", 3f15.7)') init_k0_plrn !xkf_all(1:3, ik_edge)
        CALL init_plrn_gaussian((/zero, zero, zero/), xkf_all, init_k0_plrn, eigVec)
        CALL norm_plrn_wf(eigVec, REAL(nktotf, DP))
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
          IF (select_bands_plrn(ibnd) /= band_pos) eigVec(indexkn1, 1:nstate_plrn) = czero
        ENDDO
      ENDDO
      CALL norm_plrn_wf(eigVec, REAL(nktotf, DP))
    ENDIF
    CALL stop_clock('init_Ank')
    !
    IF (debug_plrn) THEN
      ! SP: IF allocated is not recommanded
      ! FIXME
      IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
      ALLOCATE(eigvec_wan(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
      DO ik = 1, nktotf
        DO ibnd = 1, nbnd_plrn
          eigvec_wan(ik + (ibnd-1) * nktotf, 1:nstate_plrn) &
             = eigVec((ik - 1) * nbnd_plrn + ibnd, 1:nstate_plrn)
        ENDDO
      ENDDO
      DEALLOCATE(eigvec_wan, STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating eigvec_wan', 1)
    ENDIF
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
    ALLOCATE(eigvec_wan(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
    ALLOCATE(eigvec_wan_save(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan_save', 1)
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
      eigvec_wan = czero
      eigvec_wan_save = czero
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
          CALL build_plrn_bmat(Bmat, iter == 1)
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
          IF(MAXVAL(ABS(REAL(dtau))) > alat / 2.d0) THEN
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
        IF (debug_plrn) THEN
          ! SP: If allocated not recommanded
          ! FIXME
          IF (ALLOCATED(Bmat_save)) DEALLOCATE(Bmat_save)
          ALLOCATE(Bmat_save(nmodes, nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating Bmat_save', 1)
          Bmat_save = czero
          DO inu = 1, nmodes
            DO iq = 1, nqtotf
              Bmat_save(inu, iq) = Bmat(iq, inu) * wf(inu, iq)
            ENDDO
          ENDDO
          Bmat_save = czero
        ENDIF
        !
        CALL start_clock('Setup_H')
        !
        ! S Tiwari: Quick Fix for GNU compiler, should be replaced in future
        ! 
        IF (ALLOCATED(Bmat_save)) DEALLOCATE(Bmat_save)
        ALLOCATE(Bmat_save(nmodes, nqtotf), STAT = ierr)
        ! FIXME 
        CALL build_plrn_hamil(Bmat, Bmat_save, iter)
        CALL stop_clock('Setup_H')
        CALL start_clock('DiagonH')
        ! For hole polaron (type_plrn = 1),
        ! we need the highest eigenvalues instead of the lowest eigenvalues
        ! To use KS_solver, which only gives the lowest eigenvalues,
        ! we multiply -1 to the Hamiltonian to get the lowest eigenvalues
        IF (full_diagon_plrn) THEN
          ! Diagonalize Hamiltonian with Serial LAPACK subroutine
          ! Used for testing or robust benchmark
          CALL diag_serial(estmteRt, eigVec)
        ELSE
          ! Diagonalize Hamiltonian with Davidson Solver
          CALL diag_parallel(estmteRt, eigVec)
        ENDIF
        CALL stop_clock('DiagonH')
        !
        ! Reverse the eigenvalues if it is the hole polaron
        estmteRt(1:nstate_plrn) = (-type_plrn) * estmteRt(1:nstate_plrn)

        ! enforce the time-reversal symmetry: A^T_k = A_k + A^*_{-k}
        IF (time_rev_A_plrn) CALL check_time_rev_sym(eigVec)
        CALL norm_plrn_wf(eigVec, REAL(nktotf, dp))
        !
        eigVec = (- type_plrn) * eigVec
        !
        IF (debug_plrn) THEN
          ! SP: If allocated not recommanded
          ! FIXME
          IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
          ALLOCATE(eigvec_wan(nktotf * nbnd_plrn, nstate_plrn), STAT = ierr)
          IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
          DO ik = 1, nktotf
            DO ibnd = 1, nbnd_plrn
              eigvec_wan(ik + (ibnd-1) * nktotf, 1:nstate_plrn) &
                 = eigVec((ik - 1) * nbnd_plrn + ibnd, 1:nstate_plrn)
            ENDDO
          ENDDO
          DEALLOCATE(eigvec_wan, STAT = ierr)
          IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating eigvec_wan,', 1)
        ENDIF
        !
        CALL start_clock('cal_E_Form')
        CALL calc_form_energy(EPlrnPhon, EPlrnElec, EPlrnBeta)
        CALL stop_clock('cal_E_Form')
        !
        ! TODO : use exact number instead of 20 in 20e15.7
        r_cry(1:3) = IMAG(LOG(berry_phase(1:3) * EXP(- twopi * ci * r0_plrn(1:3)))) / twopi
        r_cry(1:3) = r_cry(1:3) - NINT(r_cry(1:3))
        WRITE(stdout, '(5x, i5, 60e15.4)') iter, estmteRt(1:nstate_plrn) * ryd2ev, EPlrnPhon * ryd2ev, &
           - EPlrnElec * ryd2ev, (EPlrnElec + EPlrnPhon) * ryd2ev, esterr
        eigVal = estmteRt
        totVal_save = EPlrnElec + EPlrnPhon
      ENDDO
      !
      ! Calculate and write the energies
      WRITE(stdout, '(5x, a, 50f16.7)') '      Eigenvalue (eV): ', eigVal * ryd2ev
      WRITE(stdout, '(5x, a, f16.7)')   '     Phonon part (eV): ', EPlrnPhon * ryd2ev
      WRITE(stdout, '(5x, a, f16.7)')   '   Electron part (eV): ', EPlrnElec * ryd2ev
      IF (init_plrn == 6) THEN
        WRITE(stdout, '(5x, a, f16.7)') 'Formation Energy at this \dtau (eV): ', ((-type_plrn) * eigval - EPlrnPhon) * ryd2ev
      ELSE
        WRITE(stdout, '(5x, a, f16.7)')   'Formation Energy (eV): ', (EPlrnElec + EPlrnPhon) * ryd2ev
      ENDIF
    ENDDO ! init_ntau_plrn
    !
    ! Calculate and write Density of State of Bqnu and Ank
    WRITE(stdout, '(5x, a)') "Calculating density of states to save in dos.plrn"
    CALL calc_den_of_state(eigVec, Bmat)
    !
    ! Do Bloch to Wannier transform, with U matrix
    CALL start_clock('Ank_trans')
    WRITE(stdout, '(5x, a)') "Generating the polaron wavefunction in Wannier basis to save in Amp.plrn"
    !
    ! SP: If allocated is not recommanded
    ! FIXME
    IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
    !
    ALLOCATE(eigvec_wan(nbndsub * nktotf, nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
    eigvec_wan = czero
    IF (scell_mat_plrn) THEN
      CALL scell_plrn_eigvec_tran('Bloch2Wan',.TRUE., eigVec, nktotf, nRp, Rp, nbndsub, nrr_k, &
              ndegen_k, irvec_r, dims, eigvec_wan)
    ELSE
      CALL plrn_eigvec_tran('Bloch2Wan',.TRUE., eigVec, nkf1, nkf2, nkf3, nbndsub, nrr_k, &
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
      ! Write dtau in a user-friendly format for visulization
      IF (scell_mat_plrn) THEN
        CALL scell_write_plrn_dtau_xsf(dtau, nqtotf, nRp, Rp, as, 'dtau.plrn.xsf')
      ELSE
        CALL write_plrn_dtau_xsf(dtau, nqf1, nqf2, nqf3, 'dtau.plrn.xsf')
      ENDIF
    ENDIF
    CALL stop_clock('write_files')
    !
    DEALLOCATE(dtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating dtau', 1)
    DEALLOCATE(eigvec_wan, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating eigvec_wan', 1)
    DEALLOCATE(Bmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('polaron_scf', 'Error deallocating Bmat', 1)
    !
    ! SP: If allocated is not recommanded
    ! FIXME
    IF(ALLOCATED(Rp)) DEALLOCATE(Rp)
    !
    CALL stop_clock('main_prln')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE polaron_scf
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE calc_form_energy(EPlrnPhon, EPlrnElec, EPlrnBeta)
    !-----------------------------------------------------------------------
    !!
    !! Computes the polaron formation energy
    !! Note: Require etf_all to be properly initialized
    !!
    !!
    USE constants_epw,   ONLY : zero, czero, twopi, ci, two
    USE modes,           ONLY : nmodes
    USE elph2,           ONLY : xqf, wf, nqtotf, nktotf, nkf
    USE epwcom,          ONLY : type_plrn, nstate_plrn, beta_plrn
    USE mp,              ONLY : mp_sum
    USE mp_global,       ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: EPlrnPhon
    !! FIXME
    REAL(KIND = DP), INTENT(out) :: EPlrnElec
    !! FIXME
    REAL(KIND = DP), INTENT(out) :: EPlrnBeta
    !! FIXME
    !
    ! Local variable
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
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: indexkn3
    !! FIXME
    INTEGER :: ikmbi
    !! FIXME
    INTEGER :: ikpbi
    !! FIXME
    REAL(KIND = DP):: prefix
    !! FIXME
    COMPLEX(KIND = DP) :: Q_i
    !! FIXME
    COMPLEX(KIND = DP) :: Mmn(2)
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp(2)
    !! FIXME
    !
    ! Based on Eq. 41 of Ref. 2:
    ! E_{f,ph} = 1/N_p \sum_{q\nu}|B_{q\nu}|^2\hbar\omega_{q\nu}
    ! iq -> q, nqtotf -> N_p
    ! Bmat(iq, inu) -> B_{q\nu}
    ! wf(inu, iq) -> \hbar\omega_{q\nu}
    EPlrnPhon = zero
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
        EPlrnPhon = EPlrnPhon - ABS(Bmat(iq_global, inu))**2 * (wf(inu, iq_global) / nqtotf)
      ENDDO
    ENDDO
    CALL mp_sum(EPlrnPhon, inter_pool_comm)
    !
    ! E_{f,el} = 1/N_p \sum_{nk}|A_{nk}|^2(\epsilon_{nk}-\epsilon_{F})
    ! indexkn1 -> nk, nktotf -> N_p
    ! eigVec(indexkn1, iplrn) -> A_{nk}
    ! etf_all(select_bands_plrn(ibnd), ik) - ef -> \epsilon_{nk}-\epsilon_{F}
    EPlrnElec = zero
    ! TODO: what should we do in iplrn
    DO iplrn = 1, nstate_plrn
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nqtotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik_global - 1) * nbnd_plrn + ibnd
          EPlrnElec = EPlrnElec - type_plrn * ABS(eigVec(indexkn1, iplrn))**2 / nktotf *&
             etf_all(select_bands_plrn(ibnd), ik_global)
        ENDDO
      ENDDO
    ENDDO
    CALL mp_sum(EPlrnElec, inter_pool_comm)
    !
    EPlrnBeta = zero
    !-----------------------------------------------------------------------
    END SUBROUTINE calc_form_energy
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE find_band_extreme(type_plrn, etf_all, ik_bm, band_pos, efermi)
    !-----------------------------------------------------------------------
    !!
    !! type_plrn denotes whether electron polaron (-1) or hole polaron (+1)
    !! Determine the Fermi energy, read from the input or calculated from band structure
    !!
    !-----------------------------------------------------------------------
    USE constants_epw, ONLY : zero, ryd2ev
    USE epwcom,        ONLY : efermi_read, fermi_energy
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf
    USE mp,            ONLY : mp_max, mp_min, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: type_plrn
    !! FIXME
    INTEGER, INTENT(out) :: ik_bm
    !! FIXME
    INTEGER, INTENT(out) :: band_pos
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: etf_all(:,:)
    !!
    REAL(KIND = DP), INTENT(out) :: efermi
    !! Fermi energy
    !
    ! Local variable
    INTEGER :: ik
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: k_extreme_local(npool)
    !! FIXME
    INTEGER :: ipool(1)
    !! FIXME
    REAL(KIND = DP) :: band_edge
    !! FIXME
    REAL(KIND = DP) :: extreme_local(npool)
    !! FIXME
    ! type_plrn denotes whether electron polaron (-1) or hole polaron (+1)
    IF (type_plrn == -1 ) THEN
      band_pos = select_bands_plrn(1)
    ELSE IF ( type_plrn == 1 ) THEN
      band_pos = select_bands_plrn(nbnd_plrn)
    ENDIF
    !
    WRITE(stdout, '(5x, "The band extremes are at band ",  i0)') band_pos
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
        band_edge = etf_all(band_pos, ik_global)
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
    SUBROUTINE gather_band_eigenvalues(etf, etf_all)
    !-----------------------------------------------------------------------
    !! Gather all the eigenvalues to determine the EBM/VBM,
    !! and calculate the density state of Ank and Bqnu
    !-----------------------------------------------------------------------
    USE epwcom,        ONLY : nbndsub
    USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf, xkf, xqf, wf, xkq, chw
    USE constants_epw, ONLY : zero
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: etf(:, :)
    !! Eigenvalues per cpu
    REAL(KIND = DP), INTENT(out) :: etf_all(:, :)
    !! Eigenvalues (total)
    !
    ! Local variables
    INTEGER :: ierr
    !! Error index
    REAL(KIND = DP), ALLOCATABLE :: rtmp2(:, :)
    !! FIXME
    !
    ALLOCATE(rtmp2(nbndsub, nktotf*2), STAT = ierr)
    IF (ierr /= 0) CALL errore('gather_band_eigenvalues', 'Error allocating rtmp2', 1)
    rtmp2 = zero
    !
    CALL poolgather2 ( nbndsub, nktotf*2, nkqf, etf, rtmp2  )
    etf_all(1:nbndsub, 1:nktotf) = rtmp2(1:nbndsub, 1:nktotf*2:2)
    !
    DEALLOCATE(rtmp2, STAT = ierr)
    IF (ierr /= 0) CALL errore('gather_band_eigenvalues', 'Error deallocating rtmp2', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gather_band_eigenvalues
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE cal_phonon_eigenfreq(nrr_q, irvec_q, ndegen_q, rws, nrws, wf)
    !-----------------------------------------------------------------------
    !! Calculate the phonon eigen frequencies. This is needed when restarting the polaron
    !! calculation with recalculating el-ph vertex
    !-----------------------------------------------------------------------
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : xqf, xkq, chw, nkqf, nkf, nqf, nqtotf, nktotf
    USE constants_epw, ONLY : zero, eps8, czero
    USE wan2bloch,     ONLY : dynwan2bloch, dynifc2blochf
    USE epwcom,        ONLY : type_plrn, full_diagon_plrn, lifc
    USE io_global,     ONLY : ionode, stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_q
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! FIXME
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! FIXME
    INTEGER, INTENT(in) :: nrws
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! FIXME
    REAL(KIND = DP), INTENT(out) :: wf(:, :)
    !! FIXME
    !
    ! Local variables
    LOGICAL  :: mirror_q
    !! FIXME
    INTEGER  :: inu
    !! FIXME
    INTEGER  :: ierr
    !! FIXME
    INTEGER  :: iq
    !! FIXME
    REAL(KIND = DP) :: w2(nmodes)
    !! FIXME
    REAL(KIND = DP) :: xxq(3)
    !! FIXME
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! FIXME
    !
    uf = czero
    w2 = zero
    !
    !TODO: make this part parallel over q
    wf = zero
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
          wf(inu, iq) =  DSQRT(ABS(w2(inu)))
        ELSE
          IF (ionode) THEN
            WRITE(stdout, '(5x, "WARNING: Imaginary frequency mode ",&
            &I6, " at iq=", I6)') inu, iq
          ENDIF
          wf(inu, iq) = 0.d0
        ENDIF
      ENDDO ! inu
    ENDDO ! iq
    !-----------------------------------------------------------------------
    END SUBROUTINE cal_phonon_eigenfreq
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE init_plrn_random(eigVec)
    !-----------------------------------------------------------------------
    USE elph2,         ONLY : nktotf
    USE constants_epw, ONLY : ci, cone
    USE epwcom,        ONLY : nstate_plrn
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(out) :: eigVec(:, :)
    !! FIXME
    !
    ! Local variables
    INTEGER :: ierr
    !! Error variable
    REAL(KIND = DP), ALLOCATABLE :: rmat_tmp(:, :)
    !! FIXME
    !
    CALL RANDOM_SEED()
    ALLOCATE(rmat_tmp(1:nktotf*nbnd_plrn, 1:nstate_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('init_plrn_random', 'Error allocating rmat_tmp', 1)
    CALL RANDOM_NUMBER(rmat_tmp)
    eigVec(1:nktotf * nbnd_plrn, 1:nstate_plrn) = cone * rmat_tmp(1:nktotf * nbnd_plrn, 1:nstate_plrn)
    CALL RANDOM_NUMBER(rmat_tmp)
    eigVec(1:nktotf * nbnd_plrn, 1:nstate_plrn) = eigVec(1:nktotf*nbnd_plrn, 1:nstate_plrn) + &
                                                ci * rmat_tmp(1:nktotf * nbnd_plrn, 1:nstate_plrn)
    DEALLOCATE(rmat_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('init_plrn_random', 'Error deallocating rmat_tmp', 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE init_plrn_random
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE init_plrn_gaussian(r0, xkf_all, k0, eigVec)
    !-----------------------------------------------------------------------
    !! Initialize Ank coefficients with a Gaussian lineshape
    !-----------------------------------------------------------------------
    USE constants_epw, ONLY : czero, cone, ci, twopi, one, zero
    USE epwcom,        ONLY : nstate_plrn, init_sigma_plrn
    USE elph2,         ONLY : nktotf
    USE cell_base,     ONLY : bg, alat
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: r0(3)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: k0(3)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: xkf_all(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: eigVec(:, :)
    !! FIXME
    !
    ! Local variable
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ix, iy, iz
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: ishift
    !! FIXME
    REAL(KIND = DP) :: qcart(3)
    !! FIXME
    REAL(KIND = DP) :: xxq(3)
    !! FIXME
    REAL(KIND = DP) :: shift(3)
    !! FIXME
    REAL(KIND = DP) :: disK
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    !
    ! Calculating $$ B_{\bq\nu} = \frac{1}{\omega_{\bq,\nu} N_p} \sum_\bk A^\dagger_{\bk+\bq} g_\nu(\bk,\bq) A_\bk $$
    ! \sum_\bk in local 1 to nkf first, then a inter pool sum
    ! eq to code: k -> ik, q -> iq, \nu -> inu
    ! g_\nu(\bk,\bq) -> epfall(:,:, inu, ik, iq)
    ! \omega_{\bq,\nu} -> wf(inu, iq)
    ! A_\bk -> eigVec(ik,:), A_{\bk+\bq} -> eigVec(ikq,:)
    ! B_{\bq\nu} -> Bmat(iq, inu)
    ! Whole equation translate to: Bmat(iq, inu) = one/(wf(inu, iq)*nqtotf) \sum_\bk conj(eigVec(ikq,:)) * epfall(:,:, inu, ik, iq) * eigVec(ik,:)
    ! call cal_Bmat(eigVec, wf, kpg_map, ikq_all, epfall, Bmat)
    DO ik = 1, nktotf
      xxq = xkf_all(1:3, ik) - (k0(:) - INT(k0(:))) ! shift k0 to 1BZ
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
        eigVec(indexkn1, :) = CONE * disK * ctemp
      ENDDO
    ENDDO
    !-----------------------------------------------------------------------
    END SUBROUTINE init_plrn_gaussian
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE build_plrn_bmat(Bmat, first)
    !-----------------------------------------------------------------------
    !! Create the Bmat
    !-----------------------------------------------------------------------
    USE elph2,         ONLY : nkf, nqtotf, wf, xqf, nktotf, etf
    USE epwcom,        ONLY : model_vertex_plrn, nbndsub, debug_plrn, nstate_plrn,      &
                              mixing_Plrn, type_plrn, io_lvl_plrn, m_eff_plrn,          &
                              g_start_energy_plrn, g_end_energy_plrn, g_start_band_plrn,&
                              model_enband_plrn, model_phfreq_plrn, model_vertex_plrn,  &
                              omega_LO_plrn, kappa_plrn
    USE constants_epw, ONLY : czero, one, two, zero, cone, eps2, eps8
    USE mp_world,      ONLY : mpime, world_comm
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: first
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: Bmat(:, :)
    !! FIXME
    !
    ! Local variables
    INTEGER :: iq
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: indexkn3
    !! FIXME
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: iqpg
    !! FIXME
    INTEGER :: jnu
    !! FIXME
    INTEGER :: ndegen(nmodes)
    !! FIXME
    REAL(KIND = DP) :: eig
    !! FIXME
    COMPLEX(KIND = DP) :: prefac
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    COMPLEX(KIND = DP) :: Bmat_tmp(nmodes)
    !! FIXME
    !
    Bmat = czero
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
          Bmat(iq, inu) = CONJG(Bmat(iqpg, inu))
        ENDDO
      ELSE
        DO inu = start_mode, nmodes
          Bmat(iq, inu)   = cal_Bmat(iq, inu)
        ENDDO
      ENDIF
      !TODO: to be consistent with Denny, whether this is correct?
      ! IF ( ABS(wf(1, iq)) < eps2 ) Bmat(iq, 1) = czero
    ENDDO
    ! cal_Bmat only sum over local k, so we have to do mp_sum
    CALL mp_sum(Bmat, inter_pool_comm )
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE build_plrn_bmat
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE build_plrn_hamil(Bmat, Bmat_save, iter)
    !-----------------------------------------------------------------------
    !! Create the Hamiltonian
    !-----------------------------------------------------------------------
    USE elph2,         ONLY : nkf, nqtotf, wf, xqf, nktotf, etf
    USE epwcom,        ONLY : model_vertex_plrn, nbndsub, io_lvl_plrn, nqf1, nqf2, nqf3, &
                              nstate_plrn, mixing_Plrn, type_plrn, nhblock_plrn, r0_plrn,&
                              beta_plrn, g_start_energy_plrn, g_end_energy_plrn,         &
                              g_start_band_plrn
    USE constants_epw, ONLY : czero, one, two, zero, cone, eps2, eps8, twopi, ci
    USE mp_world,      ONLY : mpime, world_comm
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE modes,         ONLY : nmodes
    USE cell_base,     ONLY : at, alat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iter
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: Bmat_save(:,:)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: Bmat(:,:)
    !! FIXME
    !
    ! Local variables
    LOGICAL, ALLOCATABLE :: saved(:)
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: iqpg
    !! FIXME
    INTEGER :: index_blk
    !! FIXME
    INTEGER :: index_loc
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: jdir
    !! FIXME
    INTEGER :: ialpha
    !! FIXME
    INTEGER :: ibpi
    !! FIXME
    INTEGER :: ibmi
    !! FIXME
    INTEGER :: ibpj
    !! FIXME
    INTEGER :: ibmj
    !! FIXME
    INTEGER :: ivec
    !! FIXME
    INTEGER :: ibm
    !! FIXME
    INTEGER :: ikpbi
    !! FIXME
    INTEGER :: ikmbi
    !! FIXME
    INTEGER :: ikpbj
    !! FIXME
    INTEGER :: ikmbj
    !! FIXME
    REAL(KIND = DP) :: F_mat(1:3,1:3)
    !! FIXME
    REAL(KIND = DP) :: b_vec(1:3)
    !! FIXME
    REAL(KIND = DP) :: eta(1:3)
    !! FIXME
    REAL(KIND = DP) :: a_i
    !! FIXME
    REAL(KIND = DP) :: a_j
    !! FIXME
    REAL(KIND = DP) :: eig
    !! FIXME
    COMPLEX(KIND = DP) :: prefac
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    COMPLEX(KIND = DP) :: Q_i
    !! FIXME
    COMPLEX(KIND = DP) :: Q_j
    !! FIXME
    COMPLEX(KIND = DP) :: Mmn(1:4)
    !! FIXME
!    COMPLEX(KIND = DP) ::
    !! FIXME
!    COMPLEX(KIND = DP) ::
    !! FIXME
!    COMPLEX(KIND = DP) ::
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: Bmat_comp(:,:)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: Hamil_tmp(:, :)
    !! FIXME
    !
    test_tags_plrn(1) = .FALSE.
    test_tags_plrn(2) = .FALSE.
    test_tags_plrn(3) = .FALSE.
    !
    ! Calculate the Hamiltonian with Bq $$H_{n\bk,n'\bk'} = \delta_{n\bk,n'\bk'}\varepsilon_{n\bk} -\frac{2}{N_p} \sum_{\nu} B^*_{\bq,\nu}g_{nn'\nu}(\bk',\bq)$$
    ! H_{n\bk,n'\bk'} -> Hamil(ik, ibnd, ikq, jbnd)
    ! B^*_{\bq,\nu} -> conj(Bmat(iq, inu))
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
          CALL start_clock('HdiagTerm')
          IF (isGVec(xqf(1:3, iq))) THEN
            ! Note that, ik is local index while ikq is global index,
            ! so even when q=0, ik \= ikq, but ik_global == ikq
            ! delta_{nn' kk'} epsilon_{nk}
            ctemp = etf_all(select_bands_plrn(ibnd), ikq)
            indexkn2 = (ikq - 1) * nbnd_plrn + ibnd
            Hamil(indexkn2, index_loc) = Hamil(indexkn2, index_loc) + ctemp
          ENDIF
          CALL stop_clock('HdiagTerm')
          !
          CALL start_clock('HOffDiagTerm')
          DO jbnd = 1, nbnd_plrn
            indexkn2 = (ikq - 1) * nbnd_plrn + jbnd
            DO inu = 1, nmodes
              ctemp = type_plrn * two / REAL(nqtotf, KIND = DP) * (Bmat(iq, inu)) * &
                 CONJG(epf(select_bands_plrn(jbnd) - g_start_band_plrn + 1, &
                       select_bands_plrn(ibnd)- g_start_band_plrn + 1, inu, ik))
              Hamil(indexkn2, index_loc) = Hamil(indexkn2, index_loc) + ctemp
            ENDDO
          ENDDO
          CALL stop_clock('HOffDiagTerm')
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
    USE constants_epw, ONLY : eps6
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xxk(3)
    !! FIXME
    !
    ! Local variable
    LOGICAL :: isGVec
    !! FIXME
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
    USE division, ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ikq
    !! FIXME
    INTEGER, INTENT(in) :: nkqtotf
    !! FIXME
    !
    ! Local variable
    INTEGER :: ikqLocal2Global
    !! FIXME
    INTEGER :: startn
    !! FIXME
    INTEGER :: lastn
    !! FIXME
    !
    CALL start_clock('ik_l2g')
    CALL fkbounds(nkqtotf, startn, lastn)
    !
    ikqLocal2Global = startn + ikq - 1
    IF (ikqLocal2Global > lastn) THEN
      CALL errore('ikqLocal2Global', 'Index of k/q is beyond this pool.', 1)
    ENDIF
    CALL stop_clock('ik_l2g')
    !
    !-----------------------------------------------------------------------
    END FUNCTION ikqLocal2Global
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION ikGlobal2Local(ik_g, nktotf)
    !-----------------------------------------------------------------------
    !! Return the global index of the local k point ik
    !-----------------------------------------------------------------------
    USE division, ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik_g
    !! FIXME
    INTEGER, INTENT(in) :: nktotf
    !! FIXME
    !
    ! Local variable
    INTEGER :: ikGlobal2Local
    !! FIXME
    INTEGER :: startn
    !! FIXME
    INTEGER :: lastn
    !! FIXME
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
    USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: energy(:)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: sigma
    !! FIXME
    REAL(KIND = DP), INTENT(out) :: f_delta(:)
    !! FIXME
    !
    f_delta = EXP(-(energy / sigma)**2)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE cal_f_delta
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE h_psi_plrn(lda, n, m, psi, hpsi)
    !-----------------------------------------------------------------------
    ! Calculate Hpsi with psi as input to use the diagon sovler in KS_solver cegterg
    ! cegterg take two external subroutine to calculate Hpsi and Spsi to calculate
    ! ( H - e S ) * evc = 0, since H and S is not saved due to their sizes
    ! Hamil need to be passed to h_psi because the parameter space is fixed
    ! to meet the requirement of Davidson diagonalization.
    !-----------------------------------------------------------------------
    USE elph2,         ONLY : nkf, nqtotf, wf, xqf, nktotf, etf
    USE epwcom,        ONLY : model_vertex_plrn, time_rev_A_plrn, nbndsub,    &
                              mixing_Plrn, type_plrn, nhblock_plrn, beta_plrn
    USE constants_epw, ONLY : czero, one, two, zero, cone, eps2, ci
    USE mp_world,      ONLY : mpime, world_comm
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE modes,         ONLY : nmodes
    USE constants_epw, ONLY : czero
    USE elph2,         ONLY : nkf, nktotf
    USE io_global,     ONLY : stdout
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
    !! k-point index
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: start_mode
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: iqpg
    !! FIXME
    INTEGER :: startn
    !! FIXME
    INTEGER :: lastn
    !! FIXME
    INTEGER :: index_loc
    !! FIXME
    INTEGER :: index_blk
    !! FIXME
    INTEGER :: ibp
    !! FIXME
    INTEGER :: ibm
    !! FIXME
    INTEGER :: ikpb
    !! FIXME
    INTEGER :: ikmb
    !! FIXME
    INTEGER :: ivec
    !! FIXME
    COMPLEX(KIND = DP) :: prefac
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp2
    !! FIXME
    COMPLEX(KIND = DP) :: hamil_kq
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: hamiltonian(:)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: hpsi_global(:, :)
    !! FIXME
    ! Gather psi (dimension nkf) to form eigVec (dimension nktotf)
    CALL start_clock('cal_hpsi')
    IF (lda < nkf * nbnd_plrn) CALL errore('h_psi_plrn', 'leading dimension of arrays psi is not correct', 1)
    eigVec = czero
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ik - 1) * nbnd_plrn + ibnd
        indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
        eigVec(indexkn2, 1:m) = psi(indexkn1, 1:m)
      ENDDO
    ENDDO
    CALL mp_sum(eigVec, inter_pool_comm)
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
               type_plrn * Hamil(indexkn2, index_loc) * eigVec(indexkn2, 1:m)
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
    SUBROUTINE get_cfac(xk, nrr_k, ndegen_k, irvec_r, dims, cfac)
    !-----------------------------------------------------------------------
    !! Compute the exponential factor.
    !-----------------------------------------------------------------------
    USE epwcom,        ONLY : use_ws
    USE constants_epw, ONLY : twopi, ci, czero
    USE kinds,         ONLY : DP, i4b
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! FIXME
    INTEGER, INTENT(in) :: dims
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: xk(3), irvec_r(3, nrr_k)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: cfac(nrr_k, dims, dims)
    !! FIXME
    !
    ! Local Variables
    INTEGER :: ikk
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: iw
    !! FIXME
    INTEGER :: iw2
    !! FIXME
    INTEGER :: ir
    !! FIXME
    REAL(KIND = DP) :: rdotk(nrr_k)
    !! FIXME
    !
    cfac = czero
    rdotk = czero
    !
    CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xk, 1, 0.0_dp, rdotk, 1 )
    !
    IF (use_ws) THEN
      DO iw = 1, dims
        DO iw2 = 1, dims
          DO ir = 1, nrr_k
            IF (ndegen_k(ir, iw2, iw) > 0) THEN
              cfac(ir, iw2, iw) = EXP(ci * rdotk(ir)) / ndegen_k(ir, iw2, iw)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSE
      cfac(:, 1, 1) = EXP(ci * rdotk(:)) / ndegen_k(:, 1, 1)
    ENDIF
    !-----------------------------------------------------------------------
    END SUBROUTINE
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    FUNCTION cal_Bmat(iq, inu)
    !-----------------------------------------------------------------------
    !! This function calculates the Bq matrix
    !! B_{qu} = 1/N_p \sum_{mnk} A^*_{mk+q}A_{nk} [g_{mnu}(k,q)/\hbar\omega_{qu}]
    !
    !-----------------------------------------------------------------------
    USE elph2,         ONLY : nkf, nktotf, wf, nqtotf
    USE epwcom,        ONLY : nstate_plrn, model_vertex_plrn
    USE epwcom,        ONLY : g_start_band_plrn
    USE constants_epw, ONLY : czero, one, eps2, cone, eps8
    USE mp_world,      ONLY : mpime, world_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! q index
    INTEGER, INTENT(in) :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    COMPLEX(KIND = DP) :: cal_Bmat
    !! B_{qu}
    COMPLEX(KIND = DP) :: prefac
    !! FIXME
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
            IF (wf(inu, iq) > eps8 ) THEN
              prefac = cone / (wf(inu, iq) * REAL(nqtotf, DP))
            ELSE
              prefac = czero
            ENDIF
            ! B_{q\nu} = \frac{1}{N_p}\sum_{nn'k}A^*_{n'k+q}\frac{g_{n'n\nu}(k, q)}{\hbar \omega_{q\nu}} A_{nk}
            cal_Bmat = cal_Bmat + prefac * (eigVec(indexkn2, iplrn)) * CONJG(eigVec(indexkn1, iplrn)) * &
               (epf(select_bands_plrn(ibnd) - g_start_band_plrn + 1, &
               select_bands_plrn(jbnd) - g_start_band_plrn + 1, inu, ik)) !conjg
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !JLB - discard zero or imaginary frequency modes
    IF (wf(inu, iq) < eps8) THEN
      cal_Bmat = czero
    ENDIF
    !-----------------------------------------------------------------------
    END FUNCTION cal_Bmat
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION indexGamma(xkf_all)
    !-----------------------------------------------------------------------
    !! Find the index of Gamma point i.e. (0, 0, 0) in xkf_all
    !! which contains all the crystal coordinates of the k/q points
    !! if Gamma point is not included, return 0
    !
    !-----------------------------------------------------------------------
    USE elph2,        ONLY : nkf, nktotf
    USE mp,           ONLY : mp_sum
    USE mp_global,    ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: xkf_all(:, :)
    !! crystal co-ordinates of k/q points
    !
    ! Local variable
    INTEGER :: indexGamma
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    !
    indexGamma = 0
    !
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      IF(isGVec(xkf_all(1:3, ik_global))) THEN
        indexGamma = ik_global
      ENDIF
    ENDDO
    CALL mp_sum(indexGamma, inter_pool_comm)
    !
    IF (.NOT. isGVec(xkf_all(1:3, indexGamma))) THEN
      CALL errore('indexGamma','The index of Gamma point is wrong!', 1)
    ENDIF
    !----------------------------------------------------------------------
    END FUNCTION indexGamma
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    SUBROUTINE norm_plrn_wf(eigVec, norm_new)
    !----------------------------------------------------------------------
    !! Computes the norm of the polaron wavefunction
    !----------------------------------------------------------------------
    USE elph2,         ONLY : nkf, nqtotf, wf, nktotf
    USE epwcom,        ONLY : nstate_plrn, time_rev_A_plrn
    USE constants_epw, ONLY : czero, one, two, cone
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: norm_new
    !! FIXME
    COMPLEX(KIND = DP), INTENT(inout) :: eigVec(:, :)
    !! FIXME
    !
    ! Local variable
    INTEGER :: iplrn
    !! FIXME
    REAL(KIND = DP) :: norm
    !! FIXME
    !
    DO iplrn = 1, nstate_plrn
      norm = REAL(DOT_PRODUCT(eigVec(1:nbnd_plrn * nktotf, iplrn), eigVec(1:nbnd_plrn * nktotf, iplrn)))
      eigVec(:, iplrn) = eigVec(:, iplrn) / DSQRT(norm) * SQRT(norm_new)
    ENDDO
    !-----------------------------------------------------------------------
    END SUBROUTINE norm_plrn_wf
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE check_time_rev_sym(eigVec)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    USE elph2,          ONLY : nkf, nqtotf, wf, nktotf
    USE epwcom,         ONLY : nstate_plrn, time_rev_A_plrn
    USE constants_epw,  ONLY : czero, one, two, cone
    USE mp,             ONLY : mp_sum
    USE mp_global,      ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(inout) :: eigVec(:, :)
    !! FIXME
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: nPlrn_l
    !! FIXME
    REAL(KIND = DP) :: norm
    !! FIXME
    COMPLEX(KIND = DP) :: temp
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: eigVec_save(:, :)
    !! FIXME
    !
    ! nstate_plrn
    nPlrn_l = 1
    !
    ALLOCATE(eigVec_save(nktotf * nbnd_plrn, nPlrn_l), STAT = ierr)
    IF (ierr /= 0) CALL errore('check_time_rev_sym', 'Error allocating eigVec_save', 1)
    eigVec_save = czero
    !
    DO ik = 1, nkf
      ik_global = ikqLocal2Global(ik, nktotf)
      ikpg = kpg_map(ik_global)
      DO ibnd = 1, nbnd_plrn
        indexkn1 = (ikpg - 1) * nbnd_plrn + ibnd
        indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
        eigVec_save(indexkn1, 1:nPlrn_l)  = CONJG(eigVec(indexkn2, 1:nPlrn_l))
      ENDDO
    ENDDO
    CALL mp_sum(eigVec_save, inter_pool_comm)
    eigVec(:, 1:nPlrn_l) = (eigVec(:, 1:nPlrn_l) + eigVec_save(:, 1:nPlrn_l))
    !
    DO iplrn = 1, nPlrn_l
      norm = REAL(DOT_PRODUCT(eigVec(1:nbnd_plrn*nktotf, iplrn), eigVec(1:nbnd_plrn*nktotf, iplrn)))!nktotf*nbnd_plrn*
      eigVec(:, iplrn) = eigVec(:, iplrn)/DSQRT(norm)
    ENDDO
    !
    DEALLOCATE(eigVec_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('check_time_rev_sym', 'Error deallocating eigVec_save', 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE check_time_rev_sym
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE diag_serial(estmteRt, eigVec)
    !-----------------------------------------------------------------------
    !! Diagonalization
    !-----------------------------------------------------------------------
    USE constants_epw,       ONLY : czero, twopi, ci, cone, zero
    USE elph2,               ONLY : nkf, nqtotf, nktotf, xkf, etf, chw
    USE epwcom,              ONLY : nstate_plrn, nkf1, nkf2, nkf3, &
                                    type_plrn, nhblock_plrn
    USE io_global,           ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,            ONLY : world_comm
    USE mp_global,           ONLY : inter_pool_comm
    USE mp,                  ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: estmteRt(:)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: eigVec(:, :)
    !! FIXME
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikk
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: icount
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
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
    REAL(KIND = DP) :: rtemp
    !! FIXME
    REAL(KIND = DP) :: xxk(3)
    !! FIXME
    REAL(KIND = DP) :: shift(3)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    COMPLEX(KIND = DP),  ALLOCATABLE :: work(:)
    !! FIXME
    COMPLEX(KIND = DP),  ALLOCATABLE :: Hamil_save(:,:)
    !! FIXME
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
      eigVec = czero
      estmteRt = zero
      !
      CALL ZHEGVX( 1, 'V', 'I', 'U', nktotf * nbnd_plrn, Hamil_save, nktotf * nbnd_plrn, Identity,&
         nktotf * nbnd_plrn, zero, zero, 1, nstate_plrn, zero, mm, estmteRt(1:nstate_plrn), &
         eigVec, nktotf * nbnd_plrn, work, lwork, rwork, iwork, ifail, info)
      !
      IF (info /= 0) CALL errore('diag_serial','Polaron: diagonal error.', 1)
      DEALLOCATE(rwork, iwork, ifail, work, Identity, STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_serial', 'Error deallocating rwork,', 1)
      ! SP: the other array above should also be deallocated ?
      ! FIXME
    ENDIF
    DEALLOCATE(Hamil_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_serial', 'Error deallocating Hamil_save,', 1)
    CALL mp_bcast(estmteRt, meta_ionode_id, world_comm)
    CALL mp_bcast(eigVec, meta_ionode_id, world_comm)
    !-----------------------------------------------------------------------
    END SUBROUTINE diag_serial
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE diag_parallel(estmteRt, eigVec)
    !-----------------------------------------------------------------------
    !! Diagonalization routine
    !-----------------------------------------------------------------------
    USE constants_epw, ONLY : czero, twopi, ci, eps5, eps6, eps4, eps2, eps8, eps10
    USE elph2,         ONLY : nkf, nqtotf, nktotf, xkf, etf, chw
    USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3, ethrdg_plrn, &
                              adapt_ethrdg_plrn, init_ethrdg_plrn, nethrdg_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,      ONLY : world_comm
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum, mp_bcast, mp_size, mp_max
    USE mp_bands,      ONLY : intra_bgrp_comm, inter_bgrp_comm, mp_start_bands
    USE mp_bands_util, ONLY : intra_bgrp_comm_ => intra_bgrp_comm, &
                              inter_bgrp_comm_ => inter_bgrp_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: estmteRt(:)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: eigVec(:, :)
    !! FIXME
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    INTEGER  :: iq
    !! FIXME
    INTEGER  :: inu
    !! FIXME
    INTEGER  :: ik
    !! FIXME
    INTEGER  :: ikk
    !! FIXME
    INTEGER  :: ikq
    !! FIXME
    INTEGER  :: ik_global
    !! FIXME
    INTEGER  :: iplrn
    !! FIXME
    INTEGER  :: ikpg
    !! FIXME
    INTEGER  :: icount
    !! FIXME
    INTEGER  :: ibnd
    !! FIXME
    INTEGER  :: jbnd
    !! FIXME
    INTEGER  :: itemp
    !! FIXME
    INTEGER  :: jtemp
    !! FIXME
    INTEGER  :: indexkn1
    !! FIXME
    INTEGER  :: indexkn2
    !! FIXME
    INTEGER  :: ithr
    !! FIXME
    INTEGER  :: nthr
    !! FiXME
    INTEGER  :: npw
    !! FIXME
    INTEGER  :: npwx
    !! FIXME
    INTEGER  :: dav_iter
    !! FIXME
    INTEGER  :: notcnv
    !! FIXME
    INTEGER  :: btype(nstate_plrn)
    !! FIXME
    INTEGER  :: nhpsi
    !! FIXME
    INTEGER, ALLOCATABLE :: iwork(:)
    !! FiXME
    INTEGER, ALLOCATABLE :: ifail(:)
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    REAL(KIND = DP) :: xxk(3)
    !! FIXME
    REAL(KIND = DP) :: shift(3)
    !! FIXME
    REAL(KIND = DP) :: ethrdg_init
    !! FIXME
    REAL(KIND = DP) :: ethrdg
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rwork(:)
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: work(:)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: psi(:, :)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: Identity(:, :)
    !! FIXME
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
      ! split eigVector (nqtotf) into parallel pieces psi (nkf), contains corresponding part with Hpsi
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
          psi(indexkn1, 1:nstate_plrn) = eigVec(indexkn2, 1:nstate_plrn)
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
        npw, npwx, nstate_plrn, nstate_plrn * 10, 1, psi, ethrdg, &
        estmteRt, btype, notcnv, .FALSE., dav_iter, nhpsi)
      CALL start_clock('cegterg_prln')
      IF(adapt_ethrdg_plrn .AND. ionode) WRITE(stdout, "(a, E14.6, I6, E14.6)") "   ", ethrdg, dav_iter, estmteRt
      IF(notcnv > 0 .AND. ionode) WRITE(stdout, "(a)") "   WARNING: Some eigenvalues not converged, &
      &check initialization, ethrdg_plrn or try adapt_ethrdg_plrn"
      !
      intra_bgrp_comm_ = itemp
      inter_bgrp_comm_ = jtemp
      !
      eigVec = czero
      DO ik = 1, nkf
        ik_global = ikqLocal2Global(ik, nktotf)
        DO ibnd = 1, nbnd_plrn
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          indexkn2 = (ik_global - 1) * nbnd_plrn + ibnd
          eigVec(indexkn2, 1:nstate_plrn) = psi(indexkn1, 1:nstate_plrn)
        ENDDO
      ENDDO
      CALL mp_sum(eigVec, inter_pool_comm)
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
    !! Write ionic positions and displacements
    !------------------------------------------------------------------------
    USE constants_epw,  ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
    USE epwcom,         ONLY : nstate_plrn
    USE io_global,      ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,       ONLY : world_comm
    USE mp,             ONLY : mp_sum, mp_bcast
    USE modes,          ONLY : nmodes
    USE ions_base,      ONLY : nat, amass, ityp, tau, atm, nsp, na, ntypx
    USE cell_base,      ONLY : at, alat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqf1
    !! FIXME
    INTEGER, INTENT(in) :: nqf2
    !! FIXME
    INTEGER, INTENT(in) :: nqf3
    !! FIXME
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: dtau(:, :)
    !! FIXME
    INTEGER, INTENT(in), OPTIONAL :: species(50)
    !! FIXME
    !
    ! Local variables
    INTEGER :: ierr
    !! Error index
    INTEGER :: wan_func_file
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: nbnd_out
    !! FIXME
    INTEGER :: nat_all
    !! FIXME
    INTEGER :: nptotf
    !! FIXME
    INTEGER :: nqf_s(1:3)
    !! FIXME
    INTEGER :: ix
    !! FIXME
    INTEGER :: iy
    !! FIXME
    INTEGER :: iz
    !! FIXME
    INTEGER :: iRp_local
    !! FIXME
    INTEGER :: iRp
    !! FIXME
    INTEGER :: iatm
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: iatm_all
    !! FIXME
    INTEGER :: iatm_sp
    !! FIXME
    INTEGER :: ika
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: isp
    !! FIXME
    INTEGER :: Rp_vec(1:3)
    !! FIXME
    INTEGER, ALLOCATABLE :: elements(:)
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    REAL(KIND = DP) :: cell(3, 3)
    !! FIXME
    REAL(KIND = DP) :: shift(1:3)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: atoms(:,:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: displacements(:,:)
    !! FIXME
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
    USE constants_epw, ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
    USE epwcom,        ONLY : nstate_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    USE modes,         ONLY : nmodes
    USE ions_base,     ONLY : nat, amass, ityp, tau, atm, nsp, na, ntypx
    USE cell_base,     ONLY : at, alat
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !! FIXME
    INTEGER, INTENT(in) :: nqtotf_p
    !! FIXME
    INTEGER, INTENT(in) :: nRp_p
    !! FIXME
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: as_p(3,3)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: dtau(:, :)
    !! FIXME
    INTEGER, INTENT(in), OPTIONAL :: species(50)
    !! FIXME
    !
    ! Local variable
    INTEGER :: ierr
    !! Error index
    INTEGER :: wan_func_file
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: nbnd_out
    !! FIXME
    INTEGER :: nat_all
    !! FIXME
    INTEGER :: nqf_s(1:3)
    !! FIXME
    INTEGER :: ix
    !! FIXME
    INTEGER :: iy
    !! FIXME
    INTEGER :: iz
    !! FIXME
    INTEGER :: iRp_local
    !! FIXME
    INTEGER :: iRp
    !! FIXME
    INTEGER :: iatm
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: iatm_all
    !! FIXME
    INTEGER :: iatm_sp
    !! FIXME
    INTEGER :: ika
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: isp
    !! FIXME
    INTEGER :: Rp_vec(1:3)
    !! FIXME
    INTEGER, ALLOCATABLE :: elements(:)
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    REAL(KIND = DP) :: cell(3, 3)
    !! FIXME
    REAL(KIND = DP) :: shift(1:3)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: atoms(:,:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: displacements(:,:)
    !! FIXME
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
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in):: filename
    !! FIXME
    INTEGER, INTENT(in) :: elements(:)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: cell(3, 3)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: atoms(:, :)
    !! FIXME
    REAL(KIND = DP), INTENT(in), OPTIONAL  :: forces(:, :)
    !! FIXME
    REAL(KIND = DP), INTENT(in), OPTIONAL  :: data_cube(:, :, :)
    !! FIXME
    !
    ! Local variables
    REAL(KIND = DP) :: rtemp
    !! FIXME
    INTEGER :: file_unit
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: nbnd_out
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: ix
    !! FIXME
    INTEGER :: iy
    !! FIXME
    INTEGER :: iz
    !! FIXME
    INTEGER :: nx
    !! FIXME
    INTEGER :: ny
    !! FIXME
    INTEGER :: nz
    !! FIXME
    INTEGER :: iatm
    !! FIXME
    INTEGER :: natm
    !! FIXME
    INTEGER :: shapeTemp(3)
    !! FIXME
    !
    ! SP: File unit cannot be define here.
    ! FIXME
    file_unit = 602
    natm = UBOUND(elements, DIM = 1)
    !
    OPEN(UNIT = file_unit, FILE = TRIM(filename), FORM = 'formatted', STATUS = 'unknown')
    !
    WRITE(file_unit, '(a)') '#'
    WRITE(file_unit, '(a)') '# Generated by the EPW polaron code'
    WRITE(file_unit, '(a)') '#'
    WRITE(file_unit, '(a)') '#'
    WRITE(file_unit, '(a)') 'CRYSTAL'
    WRITE(file_unit, '(a)') 'PRIMVEC'
    WRITE(file_unit, '(3f12.7)') cell(1:3, 1)
    WRITE(file_unit, '(3f12.7)') cell(1:3, 2)
    WRITE(file_unit, '(3f12.7)') cell(1:3, 3)
    WRITE (file_unit, '(a)') 'PRIMCOORD'
    ! The second number is always 1 for PRIMCOORD coordinates,
    ! according to http://www.xcrysden.org/doc/XSF.html
    WRITE (file_unit, '(2i6)')  natm, 1
    !
    DO iatm = 1, natm
      IF (PRESENT(forces)) THEN
        WRITE(file_unit,'(I3, 3x, 3f15.9, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm), forces(1:3, iatm)
      ELSE
        WRITE(file_unit,'(I3, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm)
      ENDIF
    ENDDO
    !
    IF(PRESENT(data_cube)) THEN
      shapeTemp = SHAPE(data_cube)
      WRITE(file_unit, '(/)')
      WRITE(file_unit, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
      WRITE(file_unit, '(3i6)') SHAPE(data_cube)
      WRITE(file_unit, '(3f12.6)') 0.0, 0.0, 0.0
      WRITE(file_unit, '(3f12.7)') cell(1:3, 1)
      WRITE(file_unit, '(3f12.7)') cell(1:3, 2)
      WRITE(file_unit, '(3f12.7)') cell(1:3, 3)
      ! TODO: data cube is probably to large to take in the same way of lattice information
      ! May be usefull and implemented in the furture
      WRITE(file_unit, *) (((data_cube(ix, iy, iz), ix = 1, shapeTemp(1)), &
        iy = 1, shapeTemp(2)), iz = 1, shapeTemp(3))
      WRITE (file_unit, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
    ENDIF
    CLOSE(file_unit)
    !----------------------------------------------------------------------------------------
    END SUBROUTINE write_xsf_file
    !----------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------
    SUBROUTINE write_plrn_wf(eigvec_wan, filename, etf_all)
    !----------------------------------------------------------------------------------------
    !! Write Polaron wavefunction
    !----------------------------------------------------------------------------------------
    USE constants_epw, ONLY : czero, ryd2ev, ryd2mev
    USE elph2,         ONLY : nkf, nqtotf, nktotf
    USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3, nbndsub, scell_mat_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! FIXME
    REAL(KIND = DP), INTENT(in), OPTIONAL :: etf_all(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: eigvec_wan(:, :)
    !! FIXME
    !
    ! Local variables
    INTEGER :: wan_func_file
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: nbnd_out
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    !
    IF(PRESENT(etf_all)) THEN
      nbnd_out = nbnd_plrn
    ELSE
      nbnd_out = nbndsub
    ENDIF
    ! SP - file number should not be defined here
    ! FIXME
    wan_func_file = 602
    !
    OPEN(UNIT = wan_func_file, FILE = TRIM(filename))
    IF (scell_mat_plrn) THEN
      WRITE(wan_func_file, '(a, 3I10)') 'Scell', nktotf, nbndsub, nstate_plrn
    ELSE
      WRITE(wan_func_file, '(6I10)') nkf1, nkf2, nkf3, nktotf, nbndsub, nstate_plrn
    ENDIF
    !
    DO ik = 1, nktotf
      DO ibnd = 1, nbnd_out
        DO iplrn = 1, nstate_plrn
          indexkn1 = (ik - 1) * nbnd_out + ibnd
          IF (PRESENT(etf_all)) THEN
            WRITE(wan_func_file, '(2I5, 4f15.7)') ik, ibnd, etf_all(select_bands_plrn(ibnd), ik) * ryd2ev, &
               eigvec_wan(indexkn1, iplrn), ABS(eigvec_wan(indexkn1, iplrn))
          ELSE
            WRITE(wan_func_file, '(2f15.7)') eigvec_wan(indexkn1, iplrn)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    CLOSE(wan_func_file)
    !----------------------------------------------------------------------------
    END SUBROUTINE write_plrn_wf
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    SUBROUTINE read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, &
               nbndsub_p, filename, scell, etf_all)
    !----------------------------------------------------------------------------
    !! Read polaron wavefunction.
    !----------------------------------------------------------------------------
    USE constants_epw, ONLY : czero
    USE elph2,         ONLY : nkf, nqtotf, nktotf
    USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! FIXME
    LOGICAL, INTENT(in), OPTIONAL  :: scell
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: eigvec_wan(:, :)
    !! FIXME
    INTEGER, INTENT(out) :: nkf1_p
    !! FIXME
    INTEGER, INTENT(out) :: nkf2_p
    !! FIXME
    INTEGER, INTENT(out) :: nkf3_p
    !! FIXME
    INTEGER, INTENT(out) :: nktotf_p
    !! FIXME
    INTEGER, INTENT(out) :: nbndsub_p
    !! FIXME
    REAL(KIND = DP), INTENT(in), OPTIONAL :: etf_all(:, :)
    !! FIXME
    !
    ! Local variables
    CHARACTER(LEN = 5) :: dmmy
    !! Dummy variables read from file
    INTEGER:: ierr
    !! Error status
    INTEGER :: wan_func_file
    !! FIXME
    INTEGER :: nPlrn_p
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: i1
    !! FIXME
    INTEGER :: i2
    !! FIXME
    REAL(KIND = DP) :: r1
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    !
    ! SP - update IF allocated
    ! FIXME
    IF (ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
    !
    IF(ionode) THEN
      ! SP - File number should be placed in io_var.f90
      ! FIXME
      wan_func_file = 602
      OPEN(UNIT = wan_func_file, FILE = TRIM(filename))
      !
      IF (PRESENT(scell) .AND. scell) THEN
        READ(wan_func_file, '(a, 3I10)') dmmy, nktotf_p, nbndsub_p, nPlrn_p
        ! nkf1_p, nkf2_p, nkf3_p should never be called if scell=.true.
        ! Just assigning an arbitrary value
        nkf1_p = 0
        nkf2_p = 0
        nkf3_p = 0
      ELSE
        READ(wan_func_file, '(6I10)') nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nPlrn_p
        IF(nkf1_p * nkf2_p * nkf3_p /= nktotf_p) THEN
          CALL errore("read_plrn_wf", filename//'Not generated from the uniform grid!', 1)
        ENDIF
      ENDIF
      !
      ALLOCATE(eigvec_wan(nbndsub_p * nktotf_p, nPlrn_p), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_plrn_wf', 'Error allocating eigvec_wan', 1)
      !
      eigvec_wan = czero
      DO ik = 1, nktotf_p
        DO ibnd = 1, nbndsub_p
          DO iplrn = 1, nPlrn_p
            indexkn1 = (ik - 1) * nbndsub_p + ibnd
            IF(PRESENT(etf_all)) THEN
              READ(wan_func_file, '(2I5, 3f15.7)') i1, i2, r1, eigvec_wan(indexkn1, iplrn)
            ELSE
              READ(wan_func_file, '(2f15.7)') eigvec_wan(indexkn1, iplrn)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CLOSE(wan_func_file)
    ENDIF
    CALL mp_bcast(nkf1_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf2_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf3_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nktotf_p, meta_ionode_id, world_comm)
    CALL mp_bcast(nPlrn_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nbndsub_p, meta_ionode_id, world_comm)
    ! SP - IF ALLOCATED should be avoided
    ! FIXME
    IF (.NOT. ALLOCATED(eigvec_wan)) THEN
      ALLOCATE(eigvec_wan(nbndsub_p * nktotf_p, nPlrn_p))
      eigvec_wan = czero
    ENDIF
    CALL mp_bcast (eigvec_wan, meta_ionode_id, world_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_wf
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE write_plrn_bmat(eigvec_wan, filename, etf_all)
    !-----------------------------------------------------------------------------------
    !! Write Bmat and phonon frequency to filename
    !!
    !-----------------------------------------------------------------------------------
    USE constants_epw, ONLY : czero, ryd2ev, ryd2mev
    USE elph2,         ONLY : nkf, nqtotf
    USE epwcom,        ONLY : nstate_plrn, nqf1, nqf2, nqf3, scell_mat_plrn
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    USE modes,         ONLY : nmodes
    USE ions_base,     ONLY : amass, ityp
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: eigvec_wan(:, :)
    !! FIXME
    REAL(KIND = DP), INTENT(in), OPTIONAL :: etf_all(:, :)
    !! FIXME
    !
    ! Local variables
    INTEGER :: wan_func_file
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: nbnd_out
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    !
    nbnd_out = nmodes
    ! SP - file index should be placed in io_var
    ! FIXME
    wan_func_file = 602
    !
    OPEN(UNIT = wan_func_file, FILE = TRIM(filename))
    IF (scell_mat_plrn) THEN
      WRITE(wan_func_file, '(a, 2I10)') 'Scell', nqtotf, nmodes
    ELSE
      WRITE(wan_func_file, '(5I10)') nqf1, nqf2, nqf3, nqtotf, nmodes
    ENDIF
    !
    DO ik = 1, nqtotf
      DO ibnd = 1, nmodes ! p
        IF (PRESENT(etf_all)) THEN ! \kappa, \alpha
          !!JLB: Changed format for improved accuracy
          WRITE(wan_func_file, '(2I5, 4ES18.10)') ik, ibnd, etf_all(ibnd, ik) * ryd2mev, &
             eigvec_wan(ik, ibnd), ABS(eigvec_wan(ik, ibnd))
        ELSE
          !JLB: Changed format for improved accuracy
          WRITE(wan_func_file, '(2ES18.10)') eigvec_wan(ik, ibnd)
        ENDIF
      ENDDO
    ENDDO
    CLOSE(wan_func_file)
    !---------------------------------------------------------------------------------
    END SUBROUTINE write_plrn_bmat
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    SUBROUTINE read_plrn_dtau(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nmodes_p,&
                  filename, scell, etf_all)
    !---------------------------------------------------------------------------------
    !! Read dtau from filename
    !---------------------------------------------------------------------------------
    USE constants_epw, ONLY : czero
    USE elph2,         ONLY : nkf, nqtotf, nktotf
    USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_sum, mp_bcast
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = *), INTENT(in) :: filename
    !! FIXME
    LOGICAL, INTENT(in), OPTIONAL :: scell
    !! FIXME
    INTEGER, INTENT(out) :: nkf1_p
    !! FIXME
    INTEGER, INTENT(out) :: nkf2_p
    !! FIXME
    INTEGER, INTENT(out) :: nkf3_p
    !! FIXME
    INTEGER, INTENT(out) :: nktotf_p
    !! FIXME
    INTEGER, INTENT(out) :: nmodes_p
    !! FIXME
    REAL(KIND = DP), INTENT(in), OPTIONAL :: etf_all(:, :) !JLB
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: eigvec_wan(:, :)
    !! FIXME
    !
    ! Local variables
    CHARACTER(LEN = 5) :: dmmy
    !! FIXME
    INTEGER :: ierr
    !! Error status
    INTEGER :: wan_func_file
    !! FIXME
    INTEGER :: nPlrn_p
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: iatm
    !! FIXME
    INTEGER :: icount
    !! FIXME
    INTEGER :: i1
    !! FIXME
    INTEGER :: i2
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    REAL(KIND = DP) :: r1
    !! FIXME
    !
    ! SP - IF ALLOCATED should be avoided
    ! FIXME
    IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
    !
    IF(ionode) THEN
      ! SP - file index should be placed in io_var
      ! FIXME
      wan_func_file = 602
      OPEN(UNIT = wan_func_file, FILE = TRIM(filename))
      !
      IF (PRESENT(scell) .AND. scell) THEN
        READ(wan_func_file, '(a, 2I10)') dmmy, nktotf_p, nmodes_p
        ! nkf1_p, nkf2_p, nkf3_p should never be called if scell=.true.
        ! Just assigning an arbitrary value
        nkf1_p = 0
        nkf2_p = 0
        nkf3_p = 0
      ELSE
        READ(wan_func_file, '(5I10)') nkf1_p, nkf2_p, nkf3_p, nktotf_p, nmodes_p
        IF (nkf1_p * nkf2_p * nkf3_p /= nktotf_p) THEN
          CALL errore('read_plrn_dtau', filename//'Not generated from the uniform grid!', 1)
        ENDIF
      ENDIF
      !
      IF(nmodes /= nmodes_p) THEN
        CALL errore('read_plrn_dtau', "Number of phonon modes are different with last run", 1)
      ENDIF
      !
      ALLOCATE(eigvec_wan(nktotf_p, nmodes_p), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_plrn_dtau', 'Error allocating eigvec_wan', 1)
      !
      eigvec_wan = czero
      DO icount = 1, nktotf_p
        DO iatm = 1, nmodes_p
          IF(PRESENT(etf_all)) THEN
            !JLB: Changed format for improved accuracy
            READ(wan_func_file, '(2I5, 3ES18.10)') i1, i2, r1, eigvec_wan(icount, iatm)
          ELSE
            !JLB: Changed format for improved accuracy
            READ(wan_func_file, '(2ES18.10)') eigvec_wan(icount, iatm)
          ENDIF
        ENDDO
      ENDDO
      CLOSE(wan_func_file)
    ENDIF
    CALL mp_bcast(nkf1_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf2_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nkf3_p,  meta_ionode_id, world_comm)
    CALL mp_bcast(nktotf_p,meta_ionode_id, world_comm)
    CALL mp_bcast(nmodes_p, meta_ionode_id, world_comm)
    ! SP - If allocated should be avoided
    ! FIXME
    IF(.NOT. ALLOCATED(eigvec_wan)) THEN
      ALLOCATE(eigvec_wan(nktotf_p, nmodes_p))
      eigvec_wan = czero
    ENDIF
    CALL mp_bcast(eigvec_wan, meta_ionode_id, world_comm)
    !--------------------------------------------------------------------------------
    END SUBROUTINE read_plrn_dtau
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    SUBROUTINE plrn_eigvec_tran(ttype, t_rev, eigVecIn, nkf1_p, nkf2_p, nkf3_p, &
               nbndsub_p, nrr_k, ndegen_k, irvec_r, dims, eigVecOut, ip_center)
    !--------------------------------------------------------------------------------
    !! Fourier transform from eigVecIn to eigVecOut
    !! ttype is 'Bloch2Wan' or 'Wan2Bloch'
    !! Parallel version, each pool calculates its own k point set (nkf),
    !! then the mp_sum is used to sum over different pools.
    !! require the correct initialization of Rp_array
    !--------------------------------------------------------------------------------
    USE constants_epw, ONLY : czero, twopi, ci, cone, two
    USE elph2,         ONLY : nkf, xkf, etf, chw, nktotf
    USE epwcom,        ONLY : nstate_plrn, nbndsub, time_rev_A_plrn, nkf1, nkf2, nkf3
    USE wan2bloch,     ONLY : hamwan2bloch !!=> hamwan2bloch_old
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! FIXME
    LOGICAL, INTENT(in) :: t_rev
    !! FIXME
    INTEGER, INTENT(in) :: nkf1_p
    !! FIXME
    INTEGER, INTENT(in) :: nkf2_p
    !! FIXME
    INTEGER, INTENT(in) :: nkf3_p
    !! FIXME
    INTEGER, INTENT(in) :: nbndsub_p
    !! FIXME
    INTEGER, INTENT(in) :: nrr_k
    !! FIXME
    INTEGER, INTENT(in) :: dims
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! FIXME
    INTEGER, INTENT(in), OPTIONAL :: ip_center(1:3)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: eigVecOut(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: eigvecIn(:, :)
    !! FIXME
    !
    ! Local variables
    LOGICAL :: is_mirror
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: itype
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikk
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: icount
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: ix
    !! FIXME
    INTEGER :: iy
    !! FIXME
    INTEGER :: iz
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    INTEGER :: i_vec(3)
    !! FIXME
    INTEGER :: center_shift(1:3)
    !! FIXME
    INTEGER :: nkf_p(3)
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    REAL(KIND = DP) :: xxk(3)
    !! FIXME
    REAL(KIND = DP) :: shift(3)
    !! FIXME
    REAL(KIND = DP) :: etf_tmp(nbndsub)
    !! FIXME
    REAL(KIND = DP) :: phi
    !! FIXME
    REAL(KIND = DP) :: maxreal !JLB
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! FIXME
    COMPLEX(KIND = DP) :: cfac(nrr_k, dims, dims)
    !! FIXME
    COMPLEX(KIND = DP) :: cufkk_k(nbndsub, nbndsub, nktotf)
    !! FIXME
    COMPLEX(KIND = DP) :: phase !JLB
    !! FIXME
    COMPLEX(KIND = DP) :: expTable(3)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkkg ( :, :, :)
    !! FIXME
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
    eigVecOut = czero
    ! S Tiwari: quickfix for ikpg
    ikpg=0
    ! FIXME
    DO ik = 1, nkf
      xxk = xkf(1:3, 2 * ik - 1)
      expTable(1:3) = EXP( twopi * ci * xxk(1:3) )
      ik_global = ikqLocal2Global(ik, nktotf)
      is_mirror = (t_rev .AND. (ik_global > ikpg))
      !
      CALL get_cfac(xxk, nrr_k, ndegen_k, irvec_r, dims, cfac)
      !
      CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
         etf_tmp, chw, cfac, dims, is_mirror)
      !
      IF(itype == 1) cufkk(1:nbndsub, 1:nbndsub) = CONJG(TRANSPOSE(cufkk(1:nbndsub, 1:nbndsub)))
      DO iplrn = 1, nstate_plrn
        !icount = 0
        !! loop over all Wannier position p
        IF (nkf1_p == 0 .OR. nkf2_p == 0 .OR. nkf3_p == 0) THEN
          CALL errore('plrn_eigvec_tran','Wrong k grid, use nkf1/2/3 to give k grid!', 1)
        ENDIF
        DO icount = 1, nkf1_p * nkf2_p * nkf3_p
          i_vec(1:3) = MODULO(index_Rp(icount, nkf_p) + center_shift, nkf_p)
          ! Same as EXP(twopi * ci * DOT_PRODUCT((/ix, iy, iz/), xxk)), to save time
          ctemp = PRODUCT(expTable(1:3)**i_vec(1:3))
          DO ibnd = 1, nbndsub_p ! loop over all Wannier state m
            DO jbnd = 1, nbnd_plrn ! loop over all Bloch state n
              indexkn1 = (icount - 1) * nbndsub + ibnd !mp
              indexkn2 = (ik_global - 1) * nbnd_plrn + jbnd !nk
              SELECT CASE(itype)
                CASE(1)  ! Bloch2Wan !
                   eigVecOut(indexkn1, iplrn) = eigVecOut(indexkn1, iplrn) + &
                      eigVecIn(indexkn2, iplrn) * ctemp / nktotf * cufkk(ibnd, select_bands_plrn(jbnd)) !JLB: Conjugate transpose taken above!
                CASE(-1) ! Wan2Bloch !
                   eigVecOut(indexkn2, iplrn) = eigVecOut(indexkn2, iplrn) + &
                      eigVecIn(indexkn1, iplrn) * CONJG(ctemp) * cufkk(select_bands_plrn(jbnd), ibnd) !JLB
              END SELECT
            ENDDO ! jbnd
          ENDDO ! ibnd
        ENDDO
      ENDDO !iplrn
    ENDDO ! ik
    ! MPI sum due to the loop ik is within local k set
    CALL mp_sum(eigVecOut, inter_pool_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE plrn_eigvec_tran
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE scell_plrn_eigvec_tran(ttype, t_rev, eigVecIn, nktotf_p, nRp_p, Rp_p, &
                                nbndsub_p, nrr_k, ndegen_k, irvec_r, dims, eigVecOut)
    !-----------------------------------------------------------------------------------
    !! JLB: Fourier transform for non-diagonal supercells
    !-----------------------------------------------------------------------------------
    USE constants_epw, ONLY : czero, twopi, ci, cone, two
    USE elph2,         ONLY : nkf, xkf, etf, chw, nktotf
    USE epwcom,        ONLY : nstate_plrn, nbndsub, time_rev_A_plrn, scell_mat_plrn
    USE wan2bloch,     ONLY : hamwan2bloch !!=> hamwan2bloch_old
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : stdout, ionode
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! FIXME
    LOGICAL, INTENT(in) :: t_rev
    !! FIXME
    INTEGER, INTENT(in) :: nktotf_p
    !! FIXME
    INTEGER, INTENT(in) :: nRp_p
    !! FIXME
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! FIXME
    INTEGER, INTENT(in) :: nbndsub_p
    !! FIXME
    INTEGER, INTENT(in) :: nrr_k
    !! FIXME
    INTEGER, INTENT(in) :: dims
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: eigVecOut(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: eigvecIn(:, :)
    !! FIXME
    ! Local Variables
    COMPLEX(KIND=dp) :: expTable(3)
    !! FIXME
    REAL(KIND = DP) :: rtemp
    !! FIXME
    REAL(KIND = DP) :: xxk(3)
    !! FIXME
    REAL(KIND = DP) :: shift(3)
    !! FIXME
    REAL(KIND = DP) :: etf_tmp(nbndsub)
    !! FIXME
    REAL(KIND = DP) :: phi
    !! FIXME
    REAL(KIND = DP) :: maxreal !JLB
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp
    !! FIXME
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! FIXME
    COMPLEX(KIND = DP) :: cfac(nrr_k, dims, dims)
    !! FIXME
    COMPLEX(KIND = DP) :: cufkk_k(nbndsub, nbndsub, nktotf)
    !! FIXME
    COMPLEX(KIND = DP) :: phase !JLB
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkkg ( :, :, :)
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: itype
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikk
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: iRp
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: ix
    !! FIXME
    INTEGER :: iy
    !! FIXME
    INTEGER :: iz
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: indexkn2
    !! FIXME
    LOGICAL :: is_mirror
    !! FIXME
    !
    IF (nbndsub_p /= nbndsub) CALL errore('scell_plrn_eigvec_tran','Different bands included in last calculation!',1)
    IF (ttype == 'Bloch2Wan') THEN
      itype =  1
    ELSEIF (ttype == 'Wan2Bloch') THEN
      itype = -1
    ELSE
      CALL errore('scell_plrn_eigvec_tran', 'Illegal translate form; should be Bloch2Wan or Wan2Bloch!', 1)
    ENDIF
    !! itype =  1 : Bloch2Wan: A_{mp} =  \frac{1}{N_p} \sum_{nk}A_{nk} \exp\left(ik\cdot R_p\right)U^\dagger_{mnk}
    !! itype = -1 : Wan2Bloch: A_{nk} = \sum_{mp}A_{mp}\exp(-ik\cdot R_p) U_{mnk}
    !! ibnd -> m, jbnd -> n
    !! R_p from 1 to nktotf_p
    !! This sequence need to be consistent every time transpose between eigvec_wann and eigvec
    eigVecOut = czero
    DO ik = 1, nkf
      xxk = xkf(1:3, 2 * ik - 1)
      ik_global = ikqLocal2Global(ik, nktotf)
      is_mirror = (t_rev .AND. (ik_global > ikpg))
      !
      CALL get_cfac(xxk, nrr_k, ndegen_k, irvec_r, dims, cfac)
      CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
         etf_tmp, chw, cfac, dims, is_mirror)
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
                   eigVecOut(indexkn1, iplrn) = eigVecOut(indexkn1, iplrn) + &
                      eigVecIn(indexkn2, iplrn) * ctemp / nktotf * cufkk(ibnd, select_bands_plrn(jbnd))
                CASE(-1) ! Wan2Bloch !
                   eigVecOut(indexkn2, iplrn) = eigVecOut(indexkn2, iplrn) + &
                      eigVecIn(indexkn1, iplrn) * CONJG(ctemp) * cufkk(select_bands_plrn(jbnd), ibnd)
              END SELECT
            ENDDO ! jbnd
          ENDDO ! ibnd
        ENDDO
      ENDDO !iplrn
    ENDDO ! ik
    ! MPI sum due to the loop ik is within local k set
    CALL mp_sum(eigVecOut, inter_pool_comm)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE scell_plrn_eigvec_tran
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE interp_plrn_wf(nrr_k, ndegen_k, irvec_r, dims)
    !-----------------------------------------------------------------------------------
    !! Interpolate Ank and write to Ank.band.plrn,
    !! especially used to visualize phonon contribution to polaron in band-mode
    !-----------------------------------------------------------------------------------
    USE constants_epw, ONLY : zero, ryd2ev, czero
    USE io_global,     ONLY : stdout, ionode
    USE epwcom,        ONLY : type_plrn, nbndsub, nstate_plrn
    USE elph2,         ONLY : nktotf, etf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! FIXME
    INTEGER, INTENT(in) :: dims
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_k(:,:,:)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! FIXME
    !
    ! Local variables
    INTEGER :: iRp
    !! FIXME
    INTEGER :: Rp_vec(3)
    !! FIXME
    INTEGER :: i_center(2)
    !! FIXME
    INTEGER :: ip_center(3)
    !! FIXME
    INTEGER :: nkf1_p
    !! FIXME
    INTEGER :: nkf2_p
    !! FIXME
    INTEGER :: nkf3_p
    !! FIXME
    INTEGER :: nktotf_p
    !! FIXME
    INTEGER :: nbndsub_p
    !! FIXME
    INTEGER :: ik_bm
    !! FIXME
    INTEGER :: band_pos
    !! FIXME
    INTEGER :: ierr
    !! FIXME
    REAL(KIND = DP) :: efermi
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: dtau_r(:, :)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! FIXME
    !
    IF (ionode) WRITE(stdout, "(5x, a)") "Start of interpolation of electronic band structure."
    IF (.NOT. ALLOCATED(etf_all)) THEN
      CALL errore('interp_plrn_wf','etf_all should be correctly prepared before calling interp_plrn_wf', 1)
    ENDIF
    !
    !!FIXME: When selecting band in solving polaron, nbndsub_p should be changed when output
    CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, 'Amp.plrn')
    !
    i_center = MAXLOC(ABS(eigvec_wan))
    !
    ip_center = index_Rp(i_center(1) / nbndsub_p + 1, (/nkf1_p, nkf2_p, nkf3_p/))
    WRITE(stdout, '(5x, a, i8, 3i5)') "The largest Amp ", i_center(1), ip_center
    !
    CALL plrn_eigvec_tran('Wan2Bloch', .false., eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
       nrr_k, ndegen_k, irvec_r, dims, eigVec, ip_center)
    !
    CALL write_plrn_wf(eigVec, 'Ank.band.plrn', etf_all)
    !
    ! SP - If allocated should be avoided
    ! FIXME
    IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE interp_plrn_wf
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE interp_plrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws)
    !-----------------------------------------------------------------------------------
    !! Interpolate bmat and write to Bmat.band.plrn,
    !! especially used to visualize phonon contribution to polaron in band-mode
    !-----------------------------------------------------------------------------------
    USE elph2,         ONLY : xqf, wf, nqtotf
    USE modes,         ONLY : nmodes
    USE constants_epw, ONLY : czero
    USE io_global,     ONLY : stdout, ionode
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_q
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! FIXME
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! FIXME
    INTEGER,  INTENT(in) :: nrws
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! FIXME
    !
    ! Local variables
    INTEGER :: nqf1_p
    !! FIXME
    INTEGER :: nqf2_p
    !! FIXME
    INTEGER :: nqf3_p
    !! FIXME
    INTEGER :: nqtotf_p
    !! FIXME
    INTEGER :: nmodes_p
    !! FIXME
    INTEGER :: ierr
    !! FIXME
    INTEGER :: iRp
    !! FIXME
    INTEGER :: ina
    !! FIXME
    INTEGER :: Rp_vec(3)
    !! FIXME
    INTEGER :: i_center(2)
    !! FIXME
    INTEGER :: ip_center(3)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: Bmat(:,:)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! FIXME
    REAL(KIND = DP),    ALLOCATABLE :: dtau_r(:, :)
    !! FIXME
    !
    CALL read_plrn_dtau(dtau, nqf1_p, nqf2_p, nqf3_p, nqtotf_p, nmodes_p, 'dtau.plrn')
    !
    ALLOCATE(dtau_r(nqtotf_p, nmodes/3), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating Bmat', 1)
    dtau_r = czero
    DO iRp = 1, nqtotf_p
      DO ina = 1, nmodes / 3 ! ika -> kappa alpha
        dtau_r(iRp, ina) = NORM2(REAL(dtau(iRp, (ina - 1) * 3 + 1:ina * 3)))
      ENDDO
    ENDDO
    i_center = MAXLOC(ABS(dtau_r))
    ip_center = index_Rp(i_center(1), (/nqf1_p, nqf2_p, nqf3_p/))
    !
    ALLOCATE(Bmat(nqtotf, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating Bmat', 1)
    Bmat = czero
    !
    CALL plrn_bmat_tran('Dtau2Bmat', .false., dtau, nqf1_p, nqf2_p, nqf3_p, &
       nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat, ip_center)
    !
    IF (ionode) CALL write_plrn_bmat(Bmat, 'Bmat.band.plrn', wf)
    !
    DEALLOCATE(dtau, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating dtau', 1)
    DEALLOCATE(Bmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating Bmat', 1)
    DEALLOCATE(dtau_r, STAT = ierr)
    IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error deallocating dtau_r', 1)
    !-----------------------------------------------------------------------------------
    END SUBROUTINE interp_plrn_bq
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    SUBROUTINE plrn_bmat_tran(ttype, t_rev, mat_in, nqf1_p, nqf2_p, nqf3_p, &
          nrr_q, ndegen_q, irvec_q, rws, nrws, mat_out, ip_center)
    !-----------------------------------------------------------------------------------
    !! Fourier transform between Bmat and dtau
    !! Dtau2Bmat : B_{q\nu} = -1/N_p\sum_{\kappa\alpha p}C_{q\kappa \nu}\Delta\tau_{\kappa\alpha p}  e_{\kappa\alpha\nu}(q)\exp(iqR_p)
    !! Bmat2Dtau : \Delta \tau_{\kappa\alpha p} = -\sum_{q\nu} 1/(C_{q\kappa \nu}) B^*_{q\nu} e_{\kappa\alpha,\nu}(q) \exp(iqR_p)
    !! C_{q\kappa \nu} = N_p\left(\frac{M_k\omega_{q\nu}}{2\hbar}\right)^{\frac{1}{2}} = N_p(M_k)^{\frac{1}{2}}D_{q\nu}
    !! D_{q \nu} = \left(\frac{\omega_{q\nu}}{2\hbar}\right)^{\frac{1}{2}}
    !-----------------------------------------------------------------------------------
    USE elph2,         ONLY : xqf, nqtotf, nkf, wf
    USE modes,         ONLY : nmodes
    USE constants_epw, ONLY : eps8, czero, one, two, twopi, zero, ci, cone
    USE ions_base,     ONLY : amass, ityp
    USE wan2bloch,     ONLY : dynwan2bloch, dynifc2blochf
    USE epwcom,        ONLY : lifc, type_plrn
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE division,      ONLY : fkbounds
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! FIXME
    LOGICAL, INTENT(in) :: t_rev
    !! FIXME
    INTEGER, INTENT(in) :: nqf1_p
    !! FIXME
    INTEGER, INTENT(in) :: nqf2_p
    !! FIXME
    INTEGER, INTENT(in) :: nqf3_p
    !! FIXME
    INTEGER, INTENT(in) :: nrr_q
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! FIXME
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! FIXME
    INTEGER, INTENT(in) :: nrws
    !! FIXME
    INTEGER, INTENT(in), OPTIONAL :: ip_center(1:3)
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: mat_in(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: mat_out(:, :)
    !! FIXME
    !
    ! Local variables
    LOGICAL :: mirror_q
    !! FIXME
    INTEGER :: jnu
    !! FIXME
    INTEGER :: ndegen(nmodes)
    !! FIXME
    INTEGER :: imode
    !! FIXME
    INTEGER :: jmode
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ierr
    !! FIXME
    INTEGER :: imu
    !! FIXME
    INTEGER :: iatm
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: itype
    !! FIXME
    INTEGER :: ika
    !! FIXME
    INTEGER :: ip_start
    !! FIXME
    INTEGER :: ip_end
    !! FIXME
    INTEGER :: iRp
    !! FIXME
    INTEGER :: nqf_p(1:3)
    !! FIXME
    INTEGER :: ix, iy, iz
    !! FIXME
    INTEGER :: ina
    !! FIXME
    INTEGER :: nqtotf_p
    !! FIXME
    INTEGER :: iqpg
    !! FIXME
    INTEGER :: nqf
    !! FIXME
    INTEGER :: nptotf
    !! FIXME
    INTEGER :: start_modes
    !! FIXME
    INTEGER :: Rp_vec(1:3)
    !! FIXME
    INTEGER :: center_shift(1:3)
    !! FIXME
    REAL(KIND = DP) :: xxq(3)
    !! FIXME
    REAL(KIND = DP) :: xxq_r(3)
    !! FIXME
    REAL(KIND = DP) :: ctemp
    !! FIXME
    REAL(KIND = DP) :: w2(nmodes)
    !! FIXME
    COMPLEX(KIND = DP) :: dtemp
    !! FIXME
    COMPLEX(KIND = DP) :: shift(3)
    !! FIXME
    COMPLEX(KIND = DP) :: expTable(3)
    !! FIXME
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! FIXME
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
      expTable(1:3) = EXP(twopi * ci * xxq(1:3))
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
      start_modes = 1
      DO inu = start_modes, nmodes ! inu -> nu
        IF (wf(inu, iq) < eps8) CYCLE !JLB - cycle zero and imaginary frequency modes
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
            !dtemp = uf_q(ika, inu, iq) * PRODUCT(expTable(1:3)**Rp_vec(1:3))
            dtemp = uf(ika, inu) * PRODUCT(expTable(1:3)**Rp_vec(1:3))
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
    USE elph2,         ONLY : xqf, nqtotf, nkf, wf
    USE modes,         ONLY : nmodes
    USE constants_epw, ONLY : eps8, czero, one, two, twopi, zero, ci, cone
    USE ions_base,     ONLY : amass, ityp
    USE wan2bloch,     ONLY : dynwan2bloch, dynifc2blochf
    USE epwcom,        ONLY : lifc, type_plrn
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE division,      ONLY : fkbounds
    USE io_global,     ONLY : stdout, ionode
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 9), INTENT(in) :: ttype
    !! FIXME
    LOGICAL, INTENT(in) :: t_rev
    !! FIXME
    INTEGER, INTENT(in) :: nqtotf_p
    !! FIXME
    INTEGER, INTENT(in) :: nRp_p
    !! FIXME
    INTEGER, INTENT(in) :: Rp_p(:,:)
    !! FIXME
    INTEGER, INTENT(in) :: nrr_q
    !! FIXME
    INTEGER, INTENT(in) :: ndegen_q(:,:,:)
    !! FIXME
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! FIXME
    INTEGER, INTENT(in) :: nrws
    !! FIXME
    REAL(KIND = DP), INTENT(in) :: rws(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: mat_in(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: mat_out(:, :)
    !! FIXME
    !
    ! Local variables
    LOGICAL :: mirror_q
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ierr
    !! FIXME
    INTEGER :: imu
    !! FIXME
    INTEGER :: iatm
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: itype
    !! FIXME
    INTEGER :: ika
    !! FIXME
    INTEGER :: ip_start
    !! FIXME
    INTEGER :: ip_end
    !! FIXME
    INTEGER :: iRp
    !! FIXME
    INTEGER :: ix, iy, iz
    !! FIXME
    INTEGER :: ina
    !! FIXME
    INTEGER :: iqpg
    !! FIXME
    INTEGER :: nqf
    !! FIXME
    INTEGER :: start_modes
    !! FIXME
    INTEGER :: jnu
    !! FIXME
    INTEGER :: ndegen(nmodes)
    !! FIXME
    INTEGER :: imode
    !! FIXME
    INTEGER :: jmode
    !! FIXME
    REAL(KIND = DP) :: xxq(3)
    !! FIXME
    REAL(KIND = DP) :: xxq_r(3)
    !! FIXME
    REAL(KIND = DP) :: ctemp
    !! FIXME
    REAL(KIND = DP) :: w2(nmodes)
    !! FIXME
    COMPLEX(KIND = DP) :: dtemp
    !! FIXME
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! FIXME
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
        IF (wf(inu, iq) < eps8) CYCLE !JLB - cycle zero and imaginary frequency modes
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
    SUBROUTINE calc_den_of_state(eigVec, Bmat)
    !-----------------------------------------------------------------------------------
    !! Compute the DOS.
    !-----------------------------------------------------------------------------------
    USE epwcom,        ONLY : nDOS_plrn, edos_max_plrn, edos_min_plrn, edos_sigma_plrn,   &
                              pdos_max_plrn, pdos_min_plrn, pdos_sigma_plrn, nstate_plrn, &
                              nkf1, nkf2, nkf3
    USE elph2,         ONLY : nkf, nqtotf, nktotf, xkf, ibndmin, ibndmax, wf
    USE modes,         ONLY : nmodes
    USE constants_epw, ONLY : ryd2mev, czero, one, ryd2ev, two, zero, cone, pi, ci, twopi,&
                              eps6, eps8, eps5
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(in) :: eigVec(:, :)
    !! FIXME
    COMPLEX(KIND = DP), INTENT(in) :: Bmat(:, :)
    !! FIXME
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    INTEGER :: idos
    !! FIXME
    INTEGER :: iq
    !! FIXME
    INTEGER :: inu
    !! FIXME
    INTEGER :: ik
    !! FIXME
    INTEGER :: ikk
    !! FIXME
    INTEGER :: ikq
    !! FIXME
    INTEGER :: ik_global
    !! FIXME
    INTEGER :: iplrn
    !! FIXME
    INTEGER :: ikpg
    !! FIXME
    INTEGER :: icount
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: ix, iy, iz
    !! FIXME
    INTEGER :: dos_file
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    REAL(KIND = DP) :: temp
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rmat_tmp(:, :)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rvec_tmp(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: edos(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: pdos(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: edos_all(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: pdos_all(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: e_grid(:)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: p_grid(:)
    !! FIXME
    !
    !Calculating DOS
    ALLOCATE(rvec_tmp(nDOS_plrn), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating rvec_tmp', 1)
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
             edos_sigma_plrn, rvec_tmp)
          indexkn1 = (ik - 1) * nbnd_plrn + ibnd
          edos = edos + (ABS(eigVec(indexkn1, iplrn))**2) * rvec_tmp
          edos_all = edos_all + rvec_tmp
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
        CALL cal_f_delta(p_grid - wf(inu, iq) * ryd2mev, pdos_sigma_plrn, rvec_tmp)
        pdos = pdos + (ABS(Bmat(iq, inu))**2) * rvec_tmp
        pdos_all = pdos_all + rvec_tmp
      ENDDO
    ENDDO
    !
    ! SP - file accession should be done in io_var
    ! FIXME
    dos_file = 601
    OPEN(UNIT = dos_file, FILE = 'dos.plrn')
    WRITE(dos_file, '(/2x, a/)') '#energy(ev)  A^2   edos  energy(mev)  B^2  pdos'
    DO idos = 1, nDOS_plrn
      WRITE(dos_file, '(6f15.7)') e_grid(idos), edos(idos), &
         edos_all(idos), p_grid(idos), pdos(idos), pdos_all(idos)
    ENDDO
    CLOSE(dos_file)
    !
    DEALLOCATE(pdos_all)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating pdos_all', 1)
    DEALLOCATE(rvec_tmp)
    IF (ierr /= 0) CALL errore('calc_den_of_state', 'Error allocating rvec_tmp', 1)
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
    !! Write real space wavefunction to file.
    !-----------------------------------------------------------------------------------
    USE constants_epw, ONLY : zero, czero, cone, twopi, ci, bohr2ang
    USE epwcom,        ONLY : nbndsub, step_wf_grid_plrn
    USE io_var,        ONLY : iun_plot
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE division,      ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 60) :: plrn_file
    !! FIXME
    INTEGER :: ierr
    !! Error status
    INTEGER :: nkf1_p
    !! FIXME
    INTEGER :: nkf2_p
    !! FIXME
    INTEGER :: nkf3_p
    !! FIXME
    INTEGER :: nktotf_p
    !! FIXME
    INTEGER :: nbnd_plrn_p
    !! FIXME
    INTEGER :: nqf_p(3)
    !! FIXME
    INTEGER :: nqtotf_p
    !! FIXME
    INTEGER :: nmodes_p
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: jbnd
    !! FIXME
    INTEGER :: iline
    !! FIXME
    INTEGER :: nAtoms
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: file_unit
    !! FIXME
    INTEGER :: igrid
    !! FIXME
    INTEGER :: itemp
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: nxx, nyy, nzz
    !! FIXME
    INTEGER :: ipx, ipy, ipz
    !! FIXME
    INTEGER :: ip_min
    !! FIXME
    INTEGER :: ip_max
    !! FIXME
    INTEGER :: ig_vec(1:3)
    !! FIXME
    INTEGER :: iscx, iscy, iscz
    !! FIXME
    INTEGER :: ie
    !! FIXME
    INTEGER :: i_species
    !! FIXME
    INTEGER :: ivec
    !! FIXME
    INTEGER :: iRp
    !! FIXME
    INTEGER :: n_grid(3)
    !! FIXME
    INTEGER :: grid_start(3)
    !! FIXME
    INTEGER :: grid_end(3)
    !! FIXME
    INTEGER :: n_grid_super(3)
    !! FIXME
    INTEGER :: r_grid_vec(3)
    !! FIXME
    INTEGER :: rpc(1:3)
    !! FIXME
    INTEGER :: ipc(1:3)
    !! FIXME
    INTEGER :: shift_start(1:3)
    !! FIXME
    INTEGER :: species(50)
    !! FIXME
    INTEGER :: Rp_vec(1:3)
    !! FIXME
    INTEGER :: shift(1:3)
    !! FIXME
    INTEGER :: ishift
    !! FIXME
    REAL(KIND = DP) :: orig(3)
    !! FIXME
    REAL(KIND = DP) :: rtempvec(4)
    !! FIXME
    REAL(KIND = DP) :: cell(3, 3)
    !! FIXME
    REAL(KIND = DP) :: rtemp(5)
    !! FIXME
    REAL(KIND = DP) :: tempDen(5, 5, 5)
    !! FIXME
    REAL(KIND = DP) :: r_cry(3)
    !! FIXME
    REAL(KIND = DP) :: r_cart(3)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: wann_func(:, :, :, :)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rvec(:)
    !! FIXME
    COMPLEX(KIND = DP) :: Amp
    !! FIXME
    COMPLEX(KIND = DP) :: ctemp(1:3)
    !! FIXME
    COMPLEX(KIND = DP) :: b_vec(1:3)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: dtau(:, :)
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: cvec(:)
    !! FIXME
    !
    ! SP - File number should be allocated in io_var
    ! FIXME
    file_unit = 60512
    ! read Amp.plrn, save eigvec_wan for the latter use
    CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbnd_plrn_p, 'Amp.plrn')
    ! read dtau.plrn, get the displacement.
    CALL read_plrn_dtau(dtau, nqf_p(1), nqf_p(2), nqf_p(3), nqtotf_p, nmodes_p, 'dtau.plrn')
    ! read cube files for the real-space Wannier function Wm(r)
    CALL read_wannier_cube(select_bands_plrn, wann_func, species, &
       n_grid, grid_start, grid_end)
    !
    cell(1:3, 1) = at(1:3, 1) * nqf_p(1) * alat
    cell(1:3, 2) = at(1:3, 2) * nqf_p(2) * alat
    cell(1:3, 3) = at(1:3, 3) * nqf_p(3) * alat
    !
    plrn_file = 'psir_plrn.xsf'
    ! Write the file head including information of structures,
    ! using the same format of
    CALL write_plrn_dtau_xsf(dtau, nqf_p(1), nqf_p(2), nqf_p(3), plrn_file, species)
    !
    orig(1:3) = zero
    n_grid_super(1:3) = nqf_p(1:3) * n_grid(1:3)
    !
    IF (ionode) THEN
      OPEN(UNIT = file_unit, FILE = TRIM(plrn_file), POSITION='APPEND')
      WRITE(file_unit, '(/)')
      WRITE(file_unit, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
      WRITE(file_unit, '(3i6)')  n_grid_super / step_wf_grid_plrn
      WRITE(file_unit, '(3f12.6)') zero, zero, zero
      WRITE(file_unit, '(3f12.7)') cell(1:3, 1) * bohr2ang
      WRITE(file_unit, '(3f12.7)') cell(1:3, 2) * bohr2ang
      WRITE(file_unit, '(3f12.7)') cell(1:3, 3) * bohr2ang
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
    rtemp = czero
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
          WRITE (file_unit, '(5e13.5)', ADVANCE='yes') ABS(cvec(::step_wf_grid_plrn))**2
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
      WRITE (file_unit, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
      CLOSE(file_unit)
      WRITE(stdout, "(5x, '|\Psi(r)|^2 written to file.')")
    ENDIF
    DEALLOCATE(wann_func, STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(cvec , STAT = ierr)
    IF (ierr /= 0) CALL errore('write_real_space_wavefunction', 'Error allocating wann_func', 1)
    DEALLOCATE(eigvec_wan , STAT = ierr)
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
    USE constants_epw, ONLY : zero, czero, cone, twopi, ci, bohr2ang
    USE epwcom,        ONLY : nbndsub, step_wf_grid_plrn, as, scell_mat
    USE io_var,        ONLY : iun_plot, iunRpscell, iunpsirscell
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE division,      ONLY : fkbounds
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
    !! FIXME
    INTEGER :: nkf1_p
    !! FIXME
    INTEGER :: nkf2_p
    !! FIXME
    INTEGER :: nkf3_p
    !! FIXME
    INTEGER :: nbnd_plrn_p
    !! FIXME
    INTEGER :: species(50)
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: indexkn1
    !! FIXME
    INTEGER :: ig_vec(1:3)
    !! FIXME
    INTEGER :: ir1, ir2, ir3
    !! FIXME
    INTEGER :: iRp1
    !! FIXME
    INTEGER :: iRp2
    !! FIXME
    INTEGER :: n_grid_total
    !! FIXME
    INTEGER :: ip_min
    !! FIXME
    INTEGER :: ip_max
    !! FIXME
    INTEGER :: n_grid(3)
    !! FIXME
    INTEGER :: grid_start(3)
    !! FIXME
    INTEGER :: grid_end(3)
    !! FIXME
    INTEGER :: n_grid_super(3)
    !! FIXME
    INTEGER :: r_in_crys_p_sup(3)
    !! FIXME
    INTEGER :: ishift
    !! FIXME
    INTEGER :: shift(3)
    !! FIXME
    REAL(KIND = DP) :: r_in_crys_p(3)
    !! FIXME
    REAL(KIND = DP) :: r_in_crys_s(3)
    !! FIXME
    REAL(KIND = DP) :: r_in_cart(3)
    !! FIXME
    REAL(KIND = DP) :: Rp_in_cart(3)
    !! FIXME
    REAL(KIND = DP) :: p2s(3, 3)
    !! FIXME
    REAL(KIND = DP) :: s2p(3,3)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: wann_func(:, :, :, :)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE :: rvec(:)
    !! FIXME
    COMPLEX(KIND = DP) :: cvec
    !! FIXME
    COMPLEX(KIND = DP), ALLOCATABLE :: eigvec_wan(:, :)
    !! FIXME
    !
    ! Broadcast supercell lattice vectors
    CALL mp_bcast(as, meta_ionode_id, world_comm)
    !
    ! read Amp.plrn, save eigvec_wan
    CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbnd_plrn_p, 'Amp.plrn', .true.)
    ! read cube files for the real-space Wannier function Wm(r)
    CALL read_wannier_cube(select_bands_plrn, wann_func, species, &
       n_grid, grid_start, grid_end)
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
    USE constants_epw, ONLY : zero, czero, cone
    USE io_var,        ONLY : iun_plot
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : stdout, ionode, meta_ionode_id
    USE cell_base,     ONLY : at, alat
    USE mp,            ONLY : mp_sum, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE division,      ONLY : fkbounds
    USE mp_global,     ONLY : inter_pool_comm
    USE epwcom,        ONLY : nbndsub
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: select_bands(:)
    !! FIXME
    INTEGER, INTENT(out) :: species(50)
    !! FIXME
    INTEGER, INTENT(out) :: n_grid(3)
    !! FIXME
    INTEGER, INTENT(out) :: grid_start_min(3)
    !! FIXME
    INTEGER, INTENT(out) :: grid_end_max(3)
    !! FIXME
    REAL(KIND = DP), ALLOCATABLE, INTENT(out) :: wann_func(:, :, :, :)
    !! FIXME
    !
    ! Local variables
    CHARACTER(LEN = 60) :: wancube
    !! FIXME
    CHARACTER(LEN = 60) :: temp_str
    !! FIXME
    INTEGER :: ierr
    !! Error status
    INTEGER :: ibnd_index
    !! FIXME
    INTEGER :: ibnd
    !! FIXME
    INTEGER :: ie
    !! FIXME
    INTEGER :: idir
    !! FIXME
    INTEGER :: i_species
    !! FIXME
    INTEGER :: nbnd
    !! FIXME
    INTEGER :: iline
    !! FIXME
    INTEGER :: nAtoms
    !! FIXME
    INTEGER :: nxx, nyy, nzz
    !! FIXME
    INTEGER :: n_len_z
    !! FIXME
    INTEGER :: grid_start(3)
    !! FIXME
    INTEGER :: grid_end(3)
    !! FIXME
    REAL(KIND = DP) :: rtempvec(4)
    !! FIXME
    REAL(KIND = DP) :: norm
    !! FIXME
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
    !! Read list of Rp unit cell vectors contained on transformed supercell
    !-----------------------------------------------------------------------------------
    USE io_var,    ONLY : iunRpscell
    USE io_global, ONLY : stdout, ionode, meta_ionode_id
    USE mp,        ONLY : mp_bcast
    USE mp_world,  ONLY : world_comm
    USE elph2,     ONLY : nqtotf
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: iRp
    !! FIXME
    INTEGER :: ierr
    !! FIXME
    INTEGER :: nRp2
    !! FIXME
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
    USE epwcom, ONLY: nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iRp
    !! FIXME
    INTEGER, INTENT(in), OPTIONAL  :: nqfs(1:3)
    !! FIXME
    !
    ! Local variable
    INTEGER  :: index_Rp(1:3)
    !! FIXME
    INTEGER  :: nqf_c(1:3)
    !! FIXME
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
    FUNCTION index_xp(delta_p)
    !-----------------------------------------------------------------------------------
    !! FIXME
    !-----------------------------------------------------------------------------------
    USE elph2,  ONLY : nqtotf
    USE epwcom, ONLY : nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: delta_p(1:3)
    !! FIXME
    !
    ! Local variable
    INTEGER :: index_xp
    !! FIXME
    !
    index_xp = delta_p(1) * nqf2 * nqf3 + delta_p(2) * nqf3 + delta_p(3) + 1
    !
    IF (.NOT. ALL(index_Rp(index_xp) == delta_p(1:3))) THEN
       CALL errore('index_xp', 'index_Rp not correct!', 1)
    ENDIF
    !-----------------------------------------------------------------------------------
    END FUNCTION index_xp
    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    FUNCTION index_shift(ishift)
    !-----------------------------------------------------------------------------------
    !! FIXME
    !-----------------------------------------------------------------------------------
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ishift
    !! FIXME
    !
    ! Local variable
    INTEGER  :: index_shift(1:3)
    !! FIXME
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
  !-----------------------------------------------------------------------------------
  END MODULE polaron
  !-----------------------------------------------------------------------------------
