!
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!-----------------------------------------------------------------------
MODULE polaron
   !!
   !! This module contains variables and subroutines of polaron
   !! Authored by Chao Lian, Weng Hong (Denny) Sio, and Jon Lafuente-Bartolome
   !!
   USE kinds,  ONLY : dp
   USE buffers,       ONLY : open_buffer, get_buffer, save_buffer, close_buffer

   IMPLICIT NONE
   !! pool-dependent Hamiltonian
   COMPLEX(DP), ALLOCATABLE :: Hamil(:, :)
   !! polaron eigenvector
   COMPLEX(DP), ALLOCATABLE :: eigVec(:, :), hEigVec(:, :)
   !! Gathered eigenvalues over the pools
   COMPLEX(DP), ALLOCATABLE :: gq_model(:)
   !! el-ph matrix element in vertex model
   COMPLEX(DP), ALLOCATABLE :: M_mat(:, :, :, :)
   !! el-ph matrix element in vertex model
   COMPLEX(DP), ALLOCATABLE :: epf(:, :, :, :), epfall(:, :, :, :, :)
   !! el-ph matrix element in vertex model
   REAL(DP),    ALLOCATABLE :: etf_model(:), wq_model
   !! band structure and phonon freq in vertex model
   REAL(DP),    ALLOCATABLE :: etf_all(:, :)
   !! Gathered k points coordinates over the pools
   REAL(DP),    ALLOCATABLE :: xkf_all(:,:) !, xkf_save(:,:)
   !! Generate the maps of k->k+q and k->G-k
   !INTEGER,     ALLOCATABLE :: Rp_array(:, :, :)
   !! Generate the maps of k->k+q and k->G-k
   INTEGER,     ALLOCATABLE :: kpg_map(:) ! remove ikq_all(:, :), too large (nkf, nktotf), calculate when in use
   !! Number of bands in subgroup used in polaron calculations
   LOGICAL,     ALLOCATABLE :: is_mirror_k(:), is_mirror_q(:), is_tri_q(:), is_tri_k(:) !JLB
   !! is this local k and global q point is a mirror point, used for time-reversal symmertry
   COMPLEX(dp), ALLOCATABLE :: rho_mat(:, :, :, :)
   !!
   INTEGER                  :: nbnd_plrn, nbnd_g_plrn!, iq_save, ik_save
   !!
   INTEGER                  :: lword_h, lword_g, lword_m, lword_umn, lword_ekanu!, iq_save, ik_save
   !!
   INTEGER                  :: io_level_g_plrn, io_level_h_plrn, io_level_umn_plrn, io_level_ekanu_plrn!, iq_save, ik_save
   !!
   INTEGER                  :: hblocksize
   !!
   INTEGER                  :: band_pos
   !!
   INTEGER                  :: ibvec(-3:3)
   !!
   INTEGER                  :: ik_edge
   !!
   INTEGER,     ALLOCATABLE :: select_bands_plrn(:)
   !!
   INTEGER,       PARAMETER :: iepfall = 12315, ihamil = 12312, iMmn = 12344, irho = 12555, iUmn = 12325, iekanu = 12326
   !! The file unit to read el-ph matrix element
   LOGICAL                  :: test_tags_plrn(20) = .false., mem_save_h = .false.
   !! The B matrix Bqu
   COMPLEX(DP),  ALLOCATABLE :: Bmat(:,:)
   !!
   COMPLEX(DP) :: berry_phase(1:3)
   !!
   ! JLB: for non-diagonal supercells
   INTEGER :: nRp
   !! Number of unit cells on supercell
   INTEGER, ALLOCATABLE :: Rp(:,:)
   !! List of unit cell vectors within supercell
   !
   PUBLIC  :: plrn_prepare, plrn_flow_select
CONTAINS
   !
   !-----------------------------------------------------------------------
   SUBROUTINE plrn_prepare(totq, iq_restart)
      USE epwcom,        ONLY : start_band_plrn, end_band_plrn, nbndsub, nstate_plrn, debug_plrn
      USE epwcom,        ONLY : cal_psir_plrn, restart_plrn,  interp_Ank_plrn, interp_Bqu_plrn
      USE epwcom,        ONLY : model_vertex_plrn, full_diagon_plrn, nhblock_plrn, Mmn_plrn, lifc
      USE epwcom,        ONLY : g_start_band_plrn, g_end_band_plrn
      USE epwcom,        ONLY : g_start_energy_plrn, g_end_energy_plrn
      USE epwcom,        ONLY : model_enband_plrn, model_phfreq_plrn, model_vertex_plrn
      USE epwcom,        ONLY : g_power_order_plrn, io_lvl_plrn
      USE epwcom,        ONLY : m_eff_plrn, kappa_plrn, omega_LO_plrn, fsthick, etf_mem
      USE epwcom,        ONLY : lrot, lphase
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf, etf
      USE modes,         ONLY : nmodes
      USE mp_world,      ONLY : mpime
      USE cell_base,     ONLY : bg, at, omega, alat
      USE constants_epw, ONLY : czero, zero, ryd2ev
      USE elph2,         ONLY : xkf, xqf, wf, xkq, chw
      USE constants_epw, ONLY : czero, cone, pi, ci, twopi, fpi, eps6, eps8, eps5
      USE epwcom,        ONLY : nkf1, nkf2, nkf3, seed_plrn
      USE epwcom,        ONLY : nqf1, nqf2, nqf3, type_plrn
      USE epwcom,        ONLY : scell_mat_plrn
      USE poolgathering, ONLY : poolgather2
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum
      USE io_files,      ONLY : check_tempdir

      IMPLICIT NONE
      INTEGER, INTENT(OUT)  :: iq_restart
      INTEGER, INTENT(IN)   :: totq
      LOGICAL               :: debug, plrn_scf, exst, pfs
      INTEGER :: inu, iq, ik, ikk, jk, iibnd, jjbnd, ibnd, jbnd, ikq, ik_global, iplrn, ierr
      INTEGER :: iter, icount, ix, iy, iz, start_mode, idos, iatm, indexkn1, indexkn2
      INTEGER :: ikGamma, iqGamma, io_level, minNBlock, ishift
      !INTEGER(kind=8), parameter :: maxword = 2**31 - 1
      INTEGER, PARAMETER :: maxword = HUGE(1)
      INTEGER :: lword_h_tmp, lword_g_tmp, lword_umn_tmp, lword_ekanu_tmp
      REAL(DP), ALLOCATABLE :: rtmp2(:,:)
      REAL(DP) :: xxk(3), xkf_cart(3, nkqf), xqf_cart(3, nqf), efermi, xkf_len(nkqf), klen, shift(3), rfac
      CHARACTER(LEN=100) :: fmt_mode

      CALL start_clock('plrn_prepare')
      IF(etf_mem == 3) THEN
         CALL errore('polaron_prepare', 'Polaron module not working with etf_mem = 3', 1)
      END IF

      WRITE(stdout, '(5x,"fsthick not working in polaron module, selecting all the k/q points.")')
      !! type_plrn denotes whether electron polaron (-1) or hole polaron (+1)
      !! Legalize the type_plrn input, in case that the user use an arbitrary number
      IF(type_plrn < 0) THEN
         type_plrn = -1
         WRITE(stdout, '(5x,"The electron polaron is calculated.")')
      ELSE
         type_plrn = 1
         WRITE(stdout, '(5x,"The hole polaron is calculated.")')
      END IF

      lrot = .true.
      lphase = .true.

      debug = debug_plrn
      IF (debug_plrn) CALL check_tempdir('test_out', exst, pfs)

      WRITE(stdout,'(5x,a)') REPEAT('=',67)



      IF(g_start_band_plrn == 0) g_start_band_plrn = 1
      IF(g_end_band_plrn == 0) g_end_band_plrn = nbndsub
      nbnd_g_plrn = g_end_band_plrn - g_start_band_plrn + 1

      IF(start_band_plrn == 0) start_band_plrn = g_start_band_plrn
      IF(end_band_plrn == 0) end_band_plrn = g_end_band_plrn
      nbnd_plrn = end_band_plrn - start_band_plrn + 1

      IF(g_start_band_plrn > start_band_plrn .or. g_end_band_plrn < end_band_plrn) THEN
         CALL errore('polaron_prepare', 'Selecting more bands in polaron than saving g matrix', 1)
      END IF

      ALLOCATE(select_bands_plrn(nbnd_plrn), STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_prepare', 'Error allocating select_bands_plrn', 1)

      select_bands_plrn = 0
      DO ibnd = 1, nbnd_plrn
         select_bands_plrn(ibnd) = start_band_plrn + ibnd - 1
      END DO

      !! copy q(x,y,z) to xkf_all, save the copy of all kpoints
      !! Note that poolgather2 has the dimension of nktotf*2,
      !! which has k at ik and k+q at ik+1
      !! This is because that xkf has the dimension of 3, nkf*2
      !! where the ik is k and ik+1 is k+q
      ALLOCATE(xkf_all(3, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating xkf_all', 1)
      ALLOCATE(rtmp2(3, nktotf*2), STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating rtmp2', 1)
      ALLOCATE(epf(nbnd_g_plrn, nbnd_g_plrn, nmodes, nkf), STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating epf', 1)
      epf = czero
      !
      xkf_all = zero
      rtmp2 = zero
      CALL poolgather2 ( 3, nktotf*2, nkqf, xkf, rtmp2)
      xkf_all(1:3, 1:nktotf) = rtmp2(1:3, 1:nktotf*2:2)

      IF(debug_plrn) THEN
         DO ik = 1, nktotf
            WRITE(stdout, '(5x, 3f15.6)') xkf_all(1:3, ik)
         END DO
      END IF

      DEALLOCATE(rtmp2)

      WRITE(stdout, "(5x, 'Use the band from ',i0, ' to ', i0, ' total ', i0)") start_band_plrn, end_band_plrn, nbnd_plrn
      WRITE(stdout, "(5x, 'Including bands: ', 10i3)") select_bands_plrn
      WRITE(stdout, "(5x, 'Use the band from ',i0, ' to ', i0, ' total ', i0, ' in saving g')") &
            g_start_band_plrn, g_end_band_plrn, nbnd_g_plrn


      WRITE(stdout, "(5x, 'Gathering eigenvalues of ', i0, ' bands and ', i0, ' k points')") nbndsub, nktotf
      ALLOCATE(etf_all(nbndsub, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating etf_all', 1)

      IF(model_enband_plrn) THEN
         etf = zero
         ! Find the distance to the nearest Gamma point in crystal coordinates
         DO ik = 1, 2*nkf
            klen = 1E3
            xxk = xkf(:, ik)
            CALL cryst_to_cart(1, xxk, bg, 1)
            !WRITE(stdout, '(3f7.4)') xxk * twopi/alat
            DO ishift = 1, 27
               shift(1:3) = REAL(index_shift(ishift), KIND=DP)
               xxk = xkf(:, ik) + shift
               CALL cryst_to_cart(1, xxk, bg, 1)
               klen = MIN(klen, NORM2(xxk))
            END DO
            etf(1, ik) = 0.5/m_eff_plrn * (klen * twopi/alat)**2
         END DO
      END IF

      etf_all = zero
      CALL gather_band_eigenvalues(etf, etf_all)

      IF(model_vertex_plrn) THEN
         WRITE(stdout, '(5x, a, f8.3)') "Using model g vertex, with order ", g_power_order_plrn
         ALLOCATE(gq_model(nqf), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating gq_model', 1)
         gq_model = zero

         rfac = SQRT(fpi/omega*omega_LO_plrn/kappa_plrn)
         DO iq = 1, nqf
            klen = 1E3
            DO ishift = 1, 27
               shift(1:3) = REAL(index_shift(ishift), KIND=DP)
               xxk = xqf(:, iq) + shift
               CALL cryst_to_cart(1, xxk, bg, 1)
               klen = MIN(klen, NORM2(xxk))
            END DO
            IF(klen > eps8) THEN
               gq_model(iq) = rfac/((klen * twopi/alat) ** g_power_order_plrn)
            END IF
         END DO
      END IF

      ! change unit from eV to Rydberg
      g_start_energy_plrn = g_start_energy_plrn / ryd2ev
      g_end_energy_plrn = g_end_energy_plrn / ryd2ev

      CALL start_clock('find_EVBM')
      CALL find_band_extreme(type_plrn, etf_all, ik_edge, band_pos, efermi)

      !! Determine the Fermi energy, read from the input or calculated from band structure
      WRITE(stdout, '(5x, "Fermi Energy is", f16.7, &
      &" (eV) located at kpoint ", i6, 3f8.3, " band ", i3)') efermi*ryd2ev, ik_edge, xkf_all(1:3, ik_edge), band_pos
      ! Shift the eigenvalues to make VBM/CBM zero
      etf_all(1:nbndsub, 1:nktotf) = etf_all(1:nbndsub, 1:nktotf) - efermi

      CALL stop_clock('find_EVBM')

      WRITE(stdout, "(5x, 'Allocating arrays and open files.')")

      IF(interp_Ank_plrn .or. interp_Bqu_plrn .or. cal_psir_plrn) THEN
         plrn_scf = .false.
         restart_plrn = .true.
      ELSE
         plrn_scf = .true.
      END IF
      IF(restart_plrn) THEN
         iq_restart = totq + 1
      ELSE
         iq_restart = 1
      END IF
      IF(interp_Ank_plrn) THEN
         ALLOCATE(eigVec(nktotf*nbnd_plrn, nstate_plrn), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_prepare', 'Error allocating eigVec', 1)
         eigVec = czero
      ELSE IF(plrn_scf) THEN
         CALL check_tempdir('plrn_tmp', exst, pfs)

         io_level_g_plrn = 1
         lword_g_tmp = nbnd_g_plrn * nbnd_g_plrn * nmodes * nkf
         IF(lword_g_tmp > maxword) THEN
            CALL errore('polaron_prepare', 'Record size larger than maximum, use more cores!', 1)
         ELSE
            lword_g = INT(lword_g_tmp)
         END IF
         IF (io_lvl_plrn == 0) THEN
            ALLOCATE(epfall(nbnd_g_plrn, nbnd_g_plrn, nmodes, nkf, nqtotf), STAT = ierr)
            IF (ierr /= 0) CALL errore('polaron_prepare', 'Error allocating epfall', 1)
            epfall = czero
         ELSE IF (io_lvl_plrn == 1) THEN
            CALL open_buffer( iepfall , 'ephf' , lword_g, io_level_g_plrn, exst, direc='plrn_tmp/') !,
         END IF

         !!JLB
         !! Open unit to write umn matrices
         !OPEN(UNIT=iUmn, ACTION='write', FILE='plrn_tmp/umn.dat')
         !! Open unit to write ekanu matrices
         !OPEN(UNIT=iekanu, ACTION='write', FILE='plrn_tmp/ekanu.dat')
         !!JLB


         io_level_h_plrn = 1
         IF(nhblock_plrn < 1 .or. nhblock_plrn > nkf * nbnd_plrn ) THEN
            CALL errore('polaron_prepare','Illegal nhblock_plrn, should between 1 and nkf * nbnd_plrn', 1)
         END IF
         minNBlock = CEILING(REAL(nkf * nbnd_plrn * nktotf * nbnd_plrn, dp) / maxword)

         IF(minNBlock >  nhblock_plrn .and. nhblock_plrn /= 1) THEN
            CALL errore('polaron_prepare', 'Record size larger than maximum, use more cores!', 1)
         END IF

         hblocksize = CEILING(REAL(nkf * nbnd_plrn, dp) / nhblock_plrn)

         lword_h_tmp = nktotf * nbnd_plrn * hblocksize

         IF(nhblock_plrn /= 1) THEN
            IF(lword_h_tmp > maxword) THEN
               CALL errore('polaron_prepare', 'Record size larger than maximum, use more cores or larger nhblock_plrn!', 1)
            ELSE
               lword_h = lword_h_tmp
               !IF(nhblock_plrn /= 1) CALL open_buffer( ihamil , 'ham' , lword_h, io_level_h_plrn, exst, direc='plrn_tmp/')
            END IF
         END IF

         lword_m = nbnd_plrn * nbnd_plrn * nktotf * 3
         !CALL open_buffer( iMmn , 'mmn' , lword_m, io_level_g_plrn, exst, direc='plrn_tmp/')
         !CALL open_buffer( irho , 'rho' , lword_m, io_level_g_plrn, exst, direc='plrn_tmp/')

         ALLOCATE(Hamil(nktotf*nbnd_plrn, hblocksize), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_prepare', 'Error allocating Hamil', 1)

         ! Allocate and initialize the variables
         ALLOCATE(eigVec(nktotf*nbnd_plrn, nstate_plrn), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_prepare', 'Error allocating eigVec', 1)
         eigVec = czero

         !! Check whether the input is legal, otherwise print warning and stop
         !! Check the input now, because if inputs are illegal, we can stop the calculation
         !! before the heavy el-ph interpolation begins.
         IF(nktotf /= nqtotf .or. nktotf < 1) CALL errore('polaron_scf','Not identical k and q grid. Do use same nkf and nqf!', 1)

         IF( (.NOT. scell_mat_plrn) .AND. (nkf1 == 0 .or. nkf2 == 0 .or. nkf3 == 0) ) THEN
            CALL errore('polaron_scf','Try to use nkf and nqf to generate k and q grid, &
               IF you are using a manual grid, also provide this information.', 1)
         END IF
         IF(nkf < 1) CALL errore('polaron_scf','Some node has no k points!', 1)
         IF(nqtotf /= nqf) CALL errore('polaron_scf','Parallel over q is not available for polaron calculations.', 1)

         WRITE(stdout, '(5x, "Polaron wavefunction calculation starts with k points ",&
            &i0, ", q points ", i0, " and KS band ", i0)') nktotf,  nqtotf,  nbnd_plrn

         !! check whether the k and q mesh are identical
         !! This may not be theoretically necessary,
         !! but necessary in this implementation
         DO ik = 1, nkf
            ik_global = ikqLocal2Global(ik, nktotf)
            IF (ANY(ABS(xkf_all(1:3, ik_global) - xqf(1:3, ik_global)) > eps6)) THEN
               CALL errore('polaron_scf', 'The k and q meshes must be exactly the same!', 1)
            END IF
         END DO

         !! map iq to G-iq, ik to G-ik
         !! find the position of Gamma point in k and q grid
         !! Note that xqf is not a MPI-local variable, xqf = xqtotf otherwise
         !! the program gives wrong results
         ikGamma = indexGamma(xkf_all)
         iqGamma = indexGamma(xqf)
         IF (ikGamma == 0) CALL errore('polaron_scf','k=0 not included in k grid!', 1)
         IF (iqGamma == 0) CALL errore('polaron_scf','q=0 not included in q grid!', 1)
         WRITE(stdout, '(5x, "The index of Gamma point in k grid is ", i0, " &
            &and q grid IS ", i0)') ikGamma, iqGamma

         !! Given k, find the index of -k+G, to impose the relation
         !! A_{n,-k+G} = A^*_{n,k} and B_{-q+G, \nu} = B^*_{q, \nu}
         !! the relation k1 + k2 = G is unique given k1
         !! Since both A and B have the dimension of nktotf/nqtotf,
         !! kpg_map should map all the k from 1 to nktotf (global)
         ALLOCATE(kpg_map(nqtotf), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating kpg_map', 1)
         kpg_map = 0

         ALLOCATE(is_mirror_k(nkf), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating is_mirror_k', 1)
         ALLOCATE(is_mirror_q(nqf), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating is_mirror_q', 1)
         !JLB
         ALLOCATE(is_tri_k(nkf), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating is_mirror_q', 1)
         is_tri_k = .false.
         ALLOCATE(is_tri_q(nqf), STAT = ierr)
         IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating is_mirror_q', 1)
         is_tri_q = .false.

         is_mirror_q = .false.
         is_mirror_k = .false.


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
            END DO
            IF (kpg_map(ik_global) == 0) CALL errore('polaron_scf', 'Not legal k/q grid!', 1)

            ikq = kpg_map(ik_global)
            IF (ik_global > ikq) THEN
               is_mirror_k(ik) = .true.
            !JLB
            ELSE IF (ik == ikq) THEN
               is_tri_k(ik) = .true.
            !JLB
            END IF
         END DO
         CALL mp_sum(kpg_map, inter_pool_comm)

         DO iq = 1, nqf
            ikq = kpg_map(iq)
            IF (iq > ikq) THEN
               is_mirror_q(iq) = .true.
            !JLB
            ELSE IF (iq == ikq) THEN
               is_tri_q(iq) = .true.
            !JLB
            END IF
         END DO


         WRITE(stdout, '(5x, a)') "Checking the k + q is included in the mesh grid for each k and q."
         !! find the global index of ik_global, ikq with vector k and k+q.
         !! Different from kpg_map, it is used in constructing Hamiltonian or hpsi
         !! ik goes over all the local k points, to parallel the program
         DO ik = 1, nkf
            DO iq = 1, nqtotf
               ik_global = ikqLocal2Global(ik, nktotf)
               xxk = xkf_all(1:3, iq) + xkf_all(1:3, ik_global)
               IF (ikq_all(ik, iq) == 0) THEN
                  CALL errore('polaron_scf','Not commensurate k and q grid!', 1)
               END IF
            END DO
         END DO



      END IF

      WRITE(stdout, "(5x, 'End of plrn_prepare')")
      CALL stop_clock('plrn_prepare')
   END SUBROUTINE
   !-----------------------------------------------------------------------
   SUBROUTINE plrn_save_g_to_file(iq, epf17, wf)
      USE modes,         ONLY : nmodes
      USE epwcom,        ONLY : g_start_band_plrn, g_end_band_plrn
      USE epwcom,        ONLY : nbndsub, g_tol_plrn, io_lvl_plrn
      USE epwcom,        ONLY : g_start_energy_plrn, g_end_energy_plrn
      USE elph2,         ONLY : nkf, nktotf, nqtotf
      USE constants_epw, ONLY : eps8, two, czero
      USE mp,            ONLY : mp_sum, mp_bcast
      USE mp_global,     ONLY : inter_pool_comm
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id

      IMPLICIT NONE
      INTEGER :: ik, ikq, imode, ibnd, jbnd, ik_global
      INTEGER, SAVE :: g_count_sum
      COMPLEX(KIND = DP), INTENT(INOUT) :: epf17(:, :, :, :)
      REAL(KIND = DP), INTENT(IN) :: wf(:, :)
      REAL(KIND = DP) :: eig
      INTEGER, INTENT(IN) :: iq
      COMPLEX(DP) :: ctemp
      ! In polaron equations, g is not epf17 but epf17/omega
      ! To ensure a Hermitian Hamiltonian, g_{mnu}(k, -q) is calculated as g*_{nmu}(k-q, q)
      DO ik = 1, nkf
         ikq = ikq_all(ik, iq)
         DO imode = 1, nmodes
            IF (wf(imode, iq) > eps8) THEN
               epf(:, :, imode, ik) = &
               epf17(g_start_band_plrn:g_end_band_plrn, g_start_band_plrn:g_end_band_plrn, imode, ik)/DSQRT( two * wf(imode, iq))
            END IF
         END DO
      END DO

      IF(io_lvl_plrn == 0) THEN
         epfall(:, :, : ,: ,iq) = epf(:, :, :, :)
      ELSE IF (io_lvl_plrn == 1) THEN
         CALL save_buffer(epf(:, :, :, :), nbnd_g_plrn*nbnd_g_plrn*nmodes*nkf, iepfall, iq)
      END IF

   END SUBROUTINE
   !-----------------------------------------------------------------------
   FUNCTION ikq_all(ik, iq)
      ! find the global index of k+q for the local ik and global iq
      USE elph2, ONLY : nktotf
      USE epwcom, ONLY : nkf1, nkf2, nkf3

      IMPLICIT NONE

      INTEGER :: ikq, ik_global, ikq_all, index_target(1:3), index_kq, ikq_loop
      INTEGER, INTENT(IN) :: ik, iq ! ik is local and iq is global

      REAL(dp) :: xxk(1:3), xxk_target(1:3)

      CALL start_clock('find_k+q')
      ikq_all = 0

      ik_global = ikqLocal2Global(ik, nktotf)
      xxk = xkf_all(1:3, iq) + xkf_all(1:3, ik_global)

      xxk_target(1:3) = xxk(1:3) - INT(xxk(1:3))
      index_target(1:3) = NINT(xxk_target(1:3) * (/nkf1, nkf2, nkf3/))

      index_kq = index_target(1) * nkf1 * nkf2 + index_target(2) * nkf2 + index_target(3) + 1

      !print *, "is this correct ", xxk(1:3), xkf_all(1:3, index_kq), index_target(1:3), isGVec(xxk - xkf_all(1:3, index_kq))
      DO ikq_loop = index_kq - 1, nktotf + index_kq
         ! ik (local) + iq (global) = ikq (global)
         ! get ikq to locate the column of the Hamiltonian
         ikq = MOD(ikq_loop, nktotf) + 1
         IF (isGVec(xxk - xkf_all(1:3, ikq))) THEN
            ikq_all = ikq
            !print *, "found ikq at ", ikq_loop - index_kq, " cycles"
            EXIT
         END IF
      END DO
      CALL stop_clock('find_k+q')

      IF (ikq_all == 0) CALL errore('ikq_all','k + q not found', 1)

   END FUNCTION
   !-----------------------------------------------------------------------
   FUNCTION find_ik(xxk, xkf_all)
      USE elph2, ONLY : nktotf
      USE epwcom, ONLY : nkf1, nkf2, nkf3

      IMPLICIT NONE

      INTEGER :: ik, find_ik

      REAL(dp), INTENT(IN) :: xxk(1:3), xkf_all(1:3, 1:nktotf)
      REAL(dp) :: xkq(1:3)
      !xxk_target(1:3)

      CALL start_clock('find_k')
      find_ik = 0
      DO ik = 1, nktotf
         xkq(1:3) = xkf_all(1:3, ik) - xxk(1:3)
         IF(isGVec(xkq)) THEN
            find_ik = ik
            EXIT
         END IF
      END DO

      IF (find_ik == 0) CALL errore('find_ik','k not found', 1)
      CALL stop_clock('find_k')
   END FUNCTION
   !-----------------------------------------------------------------------
   SUBROUTINE plrn_flow_select(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
      USE epwcom,        ONLY : cal_psir_plrn, restart_plrn,  interp_Ank_plrn, interp_Bqu_plrn
      USE epwcom,        ONLY : io_lvl_plrn, scell_mat_plrn
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id


      IMPLICIT NONE
      INTEGER,  INTENT (IN) :: nrr_k, dims, ndegen_k(:,:,:)
      REAL(DP), INTENT (IN) :: irvec_r(3, nrr_k)

      INTEGER,  INTENT (IN) :: nrr_q, ndegen_q(:,:,:)
      INTEGER,  INTENT (IN) :: irvec_q(3, nrr_q)

      INTEGER,  INTENT (IN) :: nrws
      REAL(DP), INTENT (IN) :: rws(:, :)

      LOGICAL :: itsopen

      ! Bqu Ank interpolation is not compatible with self-consistency process
      ! Added by Chao Lian for polaron calculations flow select
      ! If postprocess is ON, i.e. Bqu interpolation with saved dtau,
      ! Ank interpolation with saved Amp, and polaron visulation with saved Wannier function cube files,
      ! then self-consistent process is skipped.
      IF(.NOT. (interp_Bqu_plrn .or. interp_Ank_plrn .or. cal_psir_plrn)) THEN
         CALL polaron_scf(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
         IF (io_lvl_plrn == 1) CALL close_buffer(iepfall, 'KEEP')
         CALL close_buffer(ihamil, 'delete')
         !!JLB
         !CLOSE(UNIT=iUmn)
         !CLOSE(UNIT=iekanu)
         !!JLB
      END IF

      IF(interp_Ank_plrn) THEN
         IF(ionode) WRITE(stdout, "(5x, 'Interpolating the Ank at given k-point set....')")
         CALL interp_plrn_wf(nrr_k, ndegen_k, irvec_r, dims)
      END IF

      IF(interp_Bqu_plrn) THEN
         IF(ionode) THEN
            WRITE(stdout, "(5x, 'Interpolating the Bqu at given q-point set....')")
         END IF
         CALL interp_plrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws)
      END IF

      IF(cal_psir_plrn) THEN
         IF(ionode) WRITE(stdout, "(5x, 'Calculating the real-space distribution of polaron wavefunction....')")
         IF (scell_mat_plrn) THEN
            CALL scell_write_real_space_wavefunction()
         ELSE
            CALL write_real_space_wavefunction()
         END IF
      END IF


      !Clean up the allocated arrays, close open files
      !inquire(unit=iepfall, opened=itsopen)
      !if (itsopen) CLOSE(iepfall)

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
      !IF (ALLOCATED(Rp_array))         DEALLOCATE(Rp_array)

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
   END SUBROUTINE
   !
   SUBROUTINE polaron_scf (nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
      !
      ! Self consistency calculation of polaron wavefunction.
      ! Rewritten by Chao Lian based on the implementation by Denny Sio.
      !
      USE modes,         ONLY : nmodes
      USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, twopi
      USE constants_epw, ONLY : czero, cone, pi, ci, twopi, eps6, eps8, eps5
      USE epwcom,        ONLY : type_plrn, full_diagon_plrn, lifc, debug_plrn
      USE epwcom,        ONLY : init_sigma_plrn, init_k0_plrn
      USE epwcom,        ONLY : nstate_plrn, conv_thr_plrn
      USE epwcom,        ONLY : mixing_Plrn, init_plrn, niter_plrn, restart_plrn
      USE epwcom,        ONLY : nkf1, nkf2, nkf3, seed_plrn
      USE epwcom,        ONLY : nqf1, nqf2, nqf3, r0_plrn
      USE epwcom,        ONLY : efermi_read, fermi_energy, nbndsub
      USE epwcom,        ONLY : model_vertex_plrn, time_rev_A_plrn
      USE epwcom,        ONLY : beta_plrn, Mmn_plrn, recal_Mmn_plrn
      USE epwcom,        ONLY : model_enband_plrn, model_phfreq_plrn, model_vertex_plrn
      USE epwcom,        ONLY : omega_LO_plrn, kappa_plrn, m_eff_plrn
      USE epwcom,        ONLY : scell_mat_plrn, as
      USE epwcom,        ONLY : init_ntau_plrn
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE elph2,         ONLY : etf, ibndmin, ibndmax
      USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf
      USE elph2,         ONLY : xkf, wf, xkq, chw
      USE elph2,         ONLY : cu, cuq
      USE mp_global,     ONLY : inter_pool_comm
      USE cell_base,     ONLY : bg, alat
      USE mp,            ONLY : mp_sum, mp_bcast
      USE poolgathering, ONLY : poolgather2
      USE ions_base,     ONLY : nat
      USE mp_world,      ONLY : mpime, world_comm
      USE epwcom,        ONLY : ethrdg_plrn
      USE io_var,        ONLY : iunRpscell


      IMPLICIT NONE

      INTEGER,  INTENT (IN) :: nrr_k, dims, ndegen_k(:,:,:)
      REAL(DP), INTENT (IN) :: irvec_r(3, nrr_k)

      INTEGER,  INTENT (IN) :: nrr_q, ndegen_q(:,:,:)
      INTEGER,  INTENT (IN) :: irvec_q(3, nrr_q)

      INTEGER,  INTENT (IN) :: nrws
      REAL(DP), INTENT (IN) :: rws(:, :)
      ! local variables
      CHARACTER(LEN=256) :: filename, tmpch

      LOGICAL  :: debug, file_exist

      INTEGER  :: inu, iq, ik, ikk, jk, iibnd, jjbnd, ibnd, jbnd, ikq, ik_global, iplrn, ierr, iRp, itau
      INTEGER  :: iter, icount, ix, iy, iz, start_mode, idos, iatm, indexkn1, indexkn2
      INTEGER  :: nkf1_p, nkf2_p, nkf3_p, nbnd_plrn_p, nbndsub_p, nPlrn_p, nktotf_p, iqpg, ikpg
      INTEGER  :: dos_file, wan_func_file, bloch_func_file, bmat_file, dtau_file, itemp, jtemp
      INTEGER  :: ngrid(1:3)

      COMPLEX(DP),  ALLOCATABLE :: Bmat_save(:,:)
      COMPLEX(DP),  ALLOCATABLE :: eigvec_wan(:, :), eigvec_wan_save(:, :)
      COMPLEX(DP),  ALLOCATABLE :: dtau(:, :), dtau_save(:, :), dtau_list(:, :, :)



      COMPLEX(KIND=dp) :: cufkk ( nbndsub, nbndsub ), cfac(nrr_k, dims, dims), ctemp

      REAL(dp) :: estmteRt(nstate_plrn),  eigVal(nstate_plrn), esterr, eb, EPlrnTot, EPlrnElec, EPlrnPhon, EPlrnBeta, EPlrnDisp
      REAL(dp) :: xxk(3), xxq(3), shift(3), rtemp, disK, disK_t, prefac, norm, r_cry(1:3)
      REAL(dp) :: totVal_save, b_vec(1:3), dtau_diff

      ALLOCATE(dtau(nktotf, nmodes), STAT = ierr)
      ALLOCATE(dtau_save(nktotf, nmodes), STAT = ierr)
      dtau = czero
      dtau_save = czero
      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau', 1)
      b_vec(1:3) = (/one/nqf1, one/nqf2, one/nqf3/)
      debug = debug_plrn

      CALL start_clock('main_prln')
      !! Gather all the eigenvalues to determine the EBM/VBM,

      CALL start_clock('re_omega')
      !! Recalculate the frequency, when restart from save g
      CALL cal_phonon_eigenfreq(nrr_q, irvec_q, ndegen_q, rws, nrws, wf)
      IF(model_phfreq_plrn) THEN
         wf = zero
         wf(nmodes, :) = omega_LO_plrn
      END IF


      CALL stop_clock('re_omega')


      !! Initialize Ac(k) based on profile
      !! TODO: ik_bm should be user-adjustable
      CALL start_clock('init_Ank')
      eigVec = czero
      SELECT CASE (init_plrn)
         CASE (1)
            ! If k0 has not been set on input, center gaussian at band edge
            IF (ALL(init_k0_plrn(:)==1000.d0)) init_k0_plrn = xkf_all(1:3, ik_edge)
            !
            WRITE(stdout, '(5x, "Initializing polaron wavefunction using Gaussian wave &
               &packet with a width of", ES14.6)') init_sigma_plrn
            WRITE(stdout, '(5x, "centered at k=", 3f14.6)') init_k0_plrn !xkf_all(1:3, ik_edge)
            CALL init_plrn_gaussian((/zero, zero, zero/), xkf_all, init_k0_plrn, eigVec)
         CASE (3)
            ALLOCATE(eigvec_wan(nktotf*nbnd_plrn, nstate_plrn))
            WRITE(stdout, '(5x, a)') "Initializing the polaron wavefunction with previously saved Amp.plrn file"
            CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, 'Amp.plrn', scell_mat_plrn)
            CALL plrn_eigvec_tran('Wan2Bloch', time_rev_A_plrn, eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
               nrr_k, ndegen_k, irvec_r, dims, eigVec)
            DEALLOCATE(eigvec_wan)
         !JLB
         CASE (6)
            WRITE(stdout, '(5x, a, I6)') "Starting from displacements read from file; number of displacements:", init_ntau_plrn
            !
            ALLOCATE(dtau_list(init_ntau_plrn, nktotf, nmodes), STAT = ierr)
            IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating dtau_list', 1)
            dtau_list = CMPLX(0.d0, 0.d0)
            !
            IF (init_ntau_plrn==1) THEN
               !write (filename, '(A14)') 'dtau_disp.plrn'
               filename = 'dtau_disp.plrn'
               CALL read_plrn_dtau(dtau, nqf1, nqf2, nqf3, nqtotf, nmodes, filename, scell_mat_plrn)
               dtau_list(1, :, :) = dtau(:, :)
            ELSE
               DO itau = 1, init_ntau_plrn
                  write (tmpch,'(I4)') itau
                  filename=TRIM('dtau_disp.plrn_'//ADJUSTL(tmpch))
                  CALL read_plrn_dtau(dtau, nqf1, nqf2, nqf3, nqtotf, nmodes, filename, scell_mat_plrn)
                  dtau_list(itau, :, :) = dtau(:, :)
               END DO
            END IF
            !
            CALL mp_bcast(dtau_list, meta_ionode_id, world_comm)
            ! Initialize Ank wavefunction for iterative diagonalization
            ! Gaussian
            ! If k0 has not been set on input, center gaussian at band edge
            IF (ALL(init_k0_plrn(:)==1000.d0)) init_k0_plrn = xkf_all(1:3, ik_edge)
            WRITE(stdout, '(5x, "Initializing polaron wavefunction using Gaussian wave &
            &packet with the width of", f15.7)') init_sigma_plrn
            WRITE(stdout, '(5x, "centered at k=", 3f15.7)') init_k0_plrn !xkf_all(1:3, ik_edge)
            CALL init_plrn_gaussian((/zero, zero, zero/), xkf_all, init_k0_plrn, eigVec)
            CALL norm_plrn_wf(eigVec, REAL(nktotf, DP))
         CASE DEFAULT
            CALL errore('polaron_scf','init_plrn not implemented!', 1)
      END SELECT

      !! Only keep the coefficients in lowest/highest band,
      !! since the electron/hole localized at this band will be more stable.
      IF(init_plrn <= 2) THEN
         DO ik = 1, nktotf
            DO ibnd = 1, nbnd_plrn
               indexkn1 = (ik-1)*nbnd_plrn + ibnd
               IF(select_bands_plrn(ibnd) /= band_pos)  &
                  eigVec(indexkn1, 1:nstate_plrn) = czero
            END DO
         END DO
         CALL norm_plrn_wf(eigVec, REAL(nktotf, DP))
      END IF
      CALL stop_clock('init_Ank')

      IF(debug_plrn) THEN
         IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
         ALLOCATE(eigvec_wan(nktotf*nbnd_plrn, nstate_plrn))
         DO ik = 1, nktotf
            DO ibnd = 1, nbnd_plrn
               eigvec_wan(ik + (ibnd-1)*nktotf, 1:nstate_plrn) &
                  =  eigVec((ik - 1)*nbnd_plrn + ibnd, 1:nstate_plrn)
            END DO
         END DO
         DEALLOCATE(eigvec_wan)
      END IF


      WRITE(stdout, '(5x, "Starting the SCF cycles")')
      IF (full_diagon_plrn) THEN
         WRITE(stdout, '(5x, a)') "Using serial direct diagonalization"
      ELSE
         WRITE(stdout, '(5x, a)') "Using parallel iterative diagonalization"
         WRITE(stdout, '(5x, "Diagonalizing polaron Hamiltonian with a threshold of ", ES18.6)') ethrdg_plrn
         WRITE(stdout, '(5x, "Please check the results are convergent with this value")')
      END IF

      WRITE(stdout, '(/5x, a)') "Starting the self-consistent process"
      WRITE(stdout, '( 5x, a)') REPEAT('-',80)
!      WRITE(stdout, '(14x, 60a15)') " phonon/eV    electron/eV     Formation/eV    esterr      Eigenvalues/eV", &
!                                   "         i           j             k"
      WRITE(stdout, '(5x, " iter", 60a15)') "  Eigval/eV", "Phonon/eV", "Electron/eV", & 
                                            "Formation/eV", "Error/eV"

      ALLOCATE(Bmat(nqtotf, nmodes))
      ALLOCATE(eigvec_wan(nktotf*nbnd_plrn, nstate_plrn))
      ALLOCATE(eigvec_wan_save(nktotf*nbnd_plrn, nstate_plrn))

      !JLB
      IF (scell_mat_plrn) THEN
         CALL read_Rp_in_S()
      END IF
      !JLB

      !JLB: possibility of multiple displacements read from file, to calculate polaron energy landscape.
      !     Calculate and print the energies; .plrn files written to disk for last calculation only.
      DO itau = 1, init_ntau_plrn ! ntau_plrn=1 by default

         IF (init_plrn==6) dtau(:, :) = dtau_list(itau, :, :)

         eigvec_wan = czero
         eigvec_wan_save = czero

         estmteRt = 1E3
         esterr = 1E5
         DO iter = 1, niter_plrn
            ! Enforce the relation A_k = A*_{G-k} and normalize |A| = 1
            ! Calculating $$ B_{\bq\nu} = \frac{1}{\omega_{\bq,\nu} N_p} \sum_\bk A^\dagger_{\bk+\bq} g_\nu(\bk,\bq) A_\bk $$
            CALL start_clock('cal_bqu')
            ! JLB
            IF(init_plrn == 6 .AND. iter==1) THEN
               Bmat = czero
               IF (scell_mat_plrn) THEN
                  CALL scell_plrn_bmat_tran('Dtau2Bmat', .true., dtau, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat)
               ELSE
                  CALL plrn_bmat_tran('Dtau2Bmat', .true., dtau, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat)
               END IF
            ELSE
               CALL build_plrn_bmat(Bmat, iter==1)
               !
               IF (scell_mat_plrn) THEN
                  CALL scell_plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
               ELSE
                  CALL plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
               END IF
               !
            END IF

            dtau_diff = MAXVAL(ABS(REAL(dtau - dtau_save)))
            esterr = dtau_diff
            IF(dtau_diff < conv_thr_plrn .and. iter > 1) THEN
               IF(MAXVAL(ABS(REAL(dtau))) > alat/2.d0) THEN
                    CALL errore("polaron_scf","Non-physical solution, check initial guess and convergence.", 1) 
               END IF
               ! converged, write the final value of eigenvalue
               WRITE(stdout,'(5x,a)') REPEAT('-',80)
               WRITE(stdout, '(5x,a,f10.6,a)' )  'End of self-consistent cycle'
               EXIT
            ELSE
               dtau_save = dtau
               !CALL plrn_bmat_tran('Dtau2Bmat', .true., dtau, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat)
            END IF
            !
            CALL stop_clock('cal_bqu')

            IF(debug_plrn) THEN
               IF(ALLOCATED(Bmat_save)) DEALLOCATE(Bmat_save)
               ALLOCATE(Bmat_save(nmodes, nqtotf))
               Bmat_save = czero
               DO inu = 1, nmodes
                  DO iq = 1, nqtotf
                     Bmat_save(inu, iq) = Bmat(iq, inu)*wf(inu, iq)
                  END DO
               END DO
               Bmat_save = czero
            END IF
            !
            CALL start_clock('Setup_H')
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
            END IF
            CALL stop_clock('DiagonH')
            !
            ! Reverse the eigenvalues if it is the hole polaron
            estmteRt(1:nstate_plrn) = (-type_plrn) * estmteRt(1:nstate_plrn)

            ! enforce the time-reversal symmetry: A^T_k = A_k + A^*_{-k}
            IF(time_rev_A_plrn) CALL check_time_rev_sym(eigVec)
            CALL norm_plrn_wf(eigVec, REAL(nktotf, dp))

            eigVec = (- type_plrn) * eigVec
            !esterr = MAXVAL(ABS(estmteRt(1:nstate_plrn) - eigVal(1:nstate_plrn)))

            IF(debug_plrn) THEN
               IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
               ALLOCATE(eigvec_wan(nktotf*nbnd_plrn, nstate_plrn))
               DO ik = 1, nktotf
                  DO ibnd = 1, nbnd_plrn
                     eigvec_wan(ik + (ibnd-1)*nktotf, 1:nstate_plrn) &
                        =  eigVec((ik - 1)*nbnd_plrn + ibnd, 1:nstate_plrn)
                  END DO
               END DO
               DEALLOCATE(eigvec_wan)
            END IF

            CALL start_clock('cal_E_Form')
            CALL calc_form_energy(EPlrnPhon, EPlrnElec, EPlrnBeta)
            CALL stop_clock('cal_E_Form')

            ! TODO : use exact number instead of 20 in 20e15.7
            r_cry(1:3) = IMAG(LOG(berry_phase(1:3) * EXP(- twopi * ci * r0_plrn(1:3))))/twopi
            r_cry(1:3) = r_cry(1:3) - NINT(r_cry(1:3))
            WRITE(stdout, '(5x, i5, 60e15.4)') iter, &
               estmteRt(1:nstate_plrn)*ryd2ev, EPlrnPhon*ryd2ev, &
               - EPlrnElec*ryd2ev, (EPlrnElec + EPlrnPhon)*ryd2ev, &
               esterr
            eigVal = estmteRt
            totVal_save = EPlrnElec + EPlrnPhon
         END DO

         !! Calculate and write the energies
         WRITE(stdout, '(5x, a, 50f16.7)') '      Eigenvalue (eV): ', eigVal*ryd2ev
         WRITE(stdout, '(5x, a, f16.7)')   '     Phonon part (eV): ', EPlrnPhon*ryd2ev
         WRITE(stdout, '(5x, a, f16.7)')   '   Electron part (eV): ', EPlrnElec*ryd2ev
         IF (init_plrn==6) THEN
            WRITE(stdout, '(5x, a, f16.7)') 'Formation Energy at this \dtau (eV): ', ((- type_plrn)*eigval - EPlrnPhon)*ryd2ev
         ELSE 
            WRITE(stdout, '(5x, a, f16.7)')   'Formation Energy (eV): ', (EPlrnElec + EPlrnPhon)*ryd2ev
         END IF

      END DO ! init_ntau_plrn


      !! Calculate and write Density of State of Bqnu and Ank
      WRITE(stdout, '(5x, a)') "Calculating density of states to save in dos.plrn"
      CALL calc_den_of_state(eigVec, Bmat)
      !! Do Bloch to Wannier transform, with U matrix
      CALL start_clock('Ank_trans')
      WRITE(stdout, '(5x, a)') "Generating the polaron wavefunction in Wannier basis to save in Amp.plrn"

      IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
      ALLOCATE(eigvec_wan(nbndsub * nktotf, nstate_plrn), STAT = ierr)

      IF (ierr /= 0) CALL errore('polaron_scf', 'Error allocating eigvec_wan', 1)
      eigvec_wan = czero
      ! .true.
      IF (scell_mat_plrn) THEN
         CALL scell_plrn_eigvec_tran('Bloch2Wan',.true., eigVec, nktotf, nRp, Rp, nbndsub, nrr_k, &
                 ndegen_k, irvec_r, dims, eigvec_wan)
      ELSE
         CALL plrn_eigvec_tran('Bloch2Wan',.true., eigVec, nkf1, nkf2, nkf3, nbndsub, nrr_k, &
            ndegen_k, irvec_r, dims, eigvec_wan)
      END IF
      CALL stop_clock('Ank_trans')

      !! Calculate displacements of ions dtau, which is B matrix in Wannier basis
      CALL start_clock('Bqu_tran')
      dtau = czero
      WRITE(stdout, '(5x, a)') "Generating the ionic displacements to save in dtau.plrn and dtau.plrn.xsf"
      IF (scell_mat_plrn) THEN
         CALL scell_plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqtotf, nRp, Rp, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
      ELSE
         CALL plrn_bmat_tran('Bmat2Dtau', .true., Bmat, nqf1, nqf2, nqf3, nrr_q, ndegen_q, irvec_q, rws, nrws, dtau)
      END IF
      CALL stop_clock('Bqu_tran')
      CALL start_clock('write_files')

      IF(ionode) THEN
         !! Write Amp in Wannier basis
         CALL write_plrn_wf(eigvec_wan, 'Amp.plrn')

         !! Write Ank in Bloch basis
         CALL write_plrn_wf(eigvec,     'Ank.plrn',  etf_all)

         !! Write Bqnu
         CALL write_plrn_bmat(Bmat, 'Bmat.plrn', wf)

         !! Write dtau
         CALL write_plrn_bmat(dtau, 'dtau.plrn')

         !! Write dtau in a user-friendly format for visulization
         IF (scell_mat_plrn) THEN
            CALL scell_write_plrn_dtau_xsf(dtau, nqtotf, nRp, Rp, as, 'dtau.plrn.xsf')
         ELSE
            CALL write_plrn_dtau_xsf(dtau, nqf1, nqf2, nqf3, 'dtau.plrn.xsf')
         END IF
      END IF
      CALL stop_clock('write_files')

      ! clean up
      DEALLOCATE(dtau)
      DEALLOCATE(eigvec_wan)
      DEALLOCATE(Bmat)
      IF(ALLOCATED(Rp)) DEALLOCATE(Rp)

      CALL stop_clock('main_prln')
   END SUBROUTINE

   !! Require etf_all to be properly initialized
   SUBROUTINE calc_form_energy(EPlrnPhon, EPlrnElec, EPlrnBeta)
      USE constants_epw,   ONLY : zero, czero, twopi, ci, two
      USE modes,           ONLY : nmodes
      USE elph2,           ONLY : xqf, wf, nqtotf, nktotf, nkf
      USE epwcom,          ONLY : type_plrn, nstate_plrn, beta_plrn
      USE mp,              ONLY : mp_sum
      USE mp_global,       ONLY : inter_pool_comm

      IMPLICIT NONE

      REAL(dp), INTENT(OUT) :: EPlrnPhon, EPlrnElec, EPlrnBeta
      REAL(dp):: prefix
      INTEGER :: start_mode, inu, ibnd, jbnd, iplrn
      INTEGER :: indexkn1, indexkn2, indexkn3
      INTEGER :: ik, iq, ik_global, iq_global, idir
      INTEGER :: ikmbi, ikpbi
      COMPLEX(dp) :: Q_i, Mmn(2), ctemp(2)

      !! Based on Eq. 41 of Ref. 2:
      !! E_{f,ph} = 1/N_p \sum_{q\nu}|B_{q\nu}|^2\hbar\omega_{q\nu}
      !! iq -> q, nqtotf -> N_p
      !! Bmat(iq, inu) -> B_{q\nu}
      !! wf(inu, iq) -> \hbar\omega_{q\nu}
      EPlrnPhon = zero
      DO iq = 1, nkf
         iq_global = ikqLocal2Global(iq, nqtotf)
         !JLB - Swapped indices!
         !JLB: I think it would be better to discard modes by looking at the frequencies, i.e. discard negative or zero frequency modes.
         IF(isGVec(xqf(1:3, iq_global))) THEN
         !IF(isGVec(xqf(iq_global, 1:3))) THEN
         !JLB
            start_mode = 4
         ELSE
            start_mode = 1
         END IF
         DO inu = start_mode, nmodes
            EPlrnPhon = EPlrnPhon - ABS(Bmat(iq_global, inu))**2 * (wf(inu, iq_global)/nqtotf)
         END DO
      END DO
      CALL mp_sum(EPlrnPhon, inter_pool_comm)

      !! E_{f,el} = 1/N_p \sum_{nk}|A_{nk}|^2(\epsilon_{nk}-\epsilon_{F})
      !! indexkn1 -> nk, nktotf -> N_p
      !! eigVec(indexkn1, iplrn) -> A_{nk}
      !! etf_all(select_bands_plrn(ibnd), ik) - ef -> \epsilon_{nk}-\epsilon_{F}
      EPlrnElec = zero
      ! TODO: what should we do in iplrn
      DO iplrn = 1, nstate_plrn
         DO ik = 1, nkf
            ik_global = ikqLocal2Global(ik, nqtotf)
            DO ibnd = 1, nbnd_plrn
               indexkn1 = (ik_global - 1) * nbnd_plrn + ibnd
               EPlrnElec = EPlrnElec - type_plrn * ABS(eigVec(indexkn1, iplrn))**2/nktotf*&
                  etf_all(select_bands_plrn(ibnd), ik_global)
            END DO
         END DO
      END DO
      CALL mp_sum(EPlrnElec, inter_pool_comm)

      EPlrnBeta = zero
   END SUBROUTINE
   !
   !! type_plrn denotes whether electron polaron (-1) or hole polaron (+1)
   !! Determine the Fermi energy, read from the input or calculated from band structure
   SUBROUTINE find_band_extreme(type_plrn, etf_all, ik_bm, band_pos, efermi)
      USE constants_epw, ONLY : zero, ryd2ev
      USE epwcom,        ONLY : efermi_read, fermi_energy
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf
      USE mp,            ONLY : mp_max, mp_min, mp_sum
      USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id

      IMPLICIT NONE

      INTEGER,  INTENT(IN)  :: type_plrn
      REAL(DP), INTENT(IN)  :: etf_all(:,:)

      INTEGER,  INTENT(OUT) :: ik_bm, band_pos
      REAL(DP), INTENT(OUT) :: efermi

      REAL(DP) :: band_edge, extreme_local(npool)
      INTEGER :: ik, ik_global, k_extreme_local(npool), ipool(1)


      !! type_plrn denotes whether electron polaron (-1) or hole polaron (+1)
      IF ( type_plrn .eq. -1 ) THEN
         band_pos = select_bands_plrn(1)
      ELSE IF ( type_plrn .eq. 1 ) THEN
         band_pos = select_bands_plrn(nbnd_plrn)
      END IF
      WRITE(stdout, '(5x, "The band extremes are at band ",  i0)') band_pos

      !! Determine the Fermi energy, read from the input or calculated from band structure
      ! = 1E4*(-type_plrn)
      ik_bm = 0
      k_extreme_local = 0
      extreme_local = zero
      IF(efermi_read) THEN
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
         END IF

         DO ik = 1, nkf
            ik_global = ikqLocal2Global(ik, nktotf)
            band_edge = etf_all(band_pos, ik_global)

            IF (type_plrn == 1) THEN
               !! For hole polaron (type_plrn = 1), find the highest eigenvalue
               IF (band_edge > efermi) THEN
                  efermi = band_edge
                  ik_bm = ik_global
               END IF
            ELSE IF (type_plrn == -1) THEN
               !! For electron polaron (type_plrn = -1), find the lowest eigenvalue
               IF (band_edge < efermi) THEN
                  efermi = band_edge
                  ik_bm = ik_global
               END IF
            END IF
         END DO

         k_extreme_local(my_pool_id + 1) = ik_bm
         extreme_local(my_pool_id + 1) = efermi
         CALL mp_sum(k_extreme_local, inter_pool_comm)
         CALL mp_sum(extreme_local, inter_pool_comm)

         IF(type_plrn == 1) THEN
            ipool = MAXLOC(extreme_local)
            ik_bm = k_extreme_local(ipool(1))
            efermi = MAXVAL(extreme_local)
         ELSE IF (type_plrn == -1) THEN
            ipool = MINLOC(extreme_local)
            ik_bm = k_extreme_local(ipool(1))
            efermi = MINVAL(extreme_local)
         END IF
      END IF
   END SUBROUTINE
   !
   !! Gather all the eigenvalues to determine the EBM/VBM,
   !! and calculate the density state of Ank and Bqnu
   SUBROUTINE gather_band_eigenvalues(etf, etf_all)
      USE epwcom,        ONLY : nbndsub
      USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf
      USE elph2,         ONLY : xkf, xqf, wf, xkq, chw
      USE constants_epw, ONLY : zero
      USE poolgathering, ONLY : poolgather2

      IMPLICIT NONE

      REAL(DP),     INTENT(IN)  :: etf(:,:)
      REAL(DP),     INTENT(OUT) :: etf_all(:,:)

      INTEGER :: ierr
      REAL(DP),     ALLOCATABLE :: rtmp2(:,:)

      ALLOCATE(rtmp2(nbndsub, nktotf*2), STAT = ierr)
      IF (ierr /= 0) CALL errore('gather_band_eigenvalues', 'Error allocating rtmp2', 1)
      rtmp2 = zero

      CALL poolgather2 ( nbndsub, nktotf*2, nkqf, etf, rtmp2  )
      etf_all(1:nbndsub, 1:nktotf) = rtmp2(1:nbndsub, 1:nktotf*2:2)

      DEALLOCATE(rtmp2)
   END SUBROUTINE
   !
   ! Calculate the phonon eigen frequencies. This is needed when restarting the polaron
   ! calculation with recalculating el-ph vertex
   SUBROUTINE cal_phonon_eigenfreq(nrr_q, irvec_q, ndegen_q, rws, nrws, wf)
      USE modes,         ONLY : nmodes
      USE elph2,         ONLY : xqf, xkq, chw
      USE constants_epw, ONLY : zero, eps8, czero
      USE wan2bloch,     ONLY : dynwan2bloch, dynifc2blochf
      USE elph2,         ONLY : nkqf, nkf, nqf, nqtotf, nktotf
      USE epwcom,        ONLY : type_plrn, full_diagon_plrn, lifc
      USE io_global,     ONLY : ionode, stdout

      IMPLICIT NONE

      INTEGER,  INTENT (IN) :: nrr_q, ndegen_q(:,:,:)
      INTEGER,  INTENT (IN) :: irvec_q(3, nrr_q)
      INTEGER,  INTENT (IN) :: nrws
      REAL(DP), INTENT (IN) :: rws(:, :)

      REAL(DP), INTENT (OUT) :: wf(:, :)

      COMPLEX(DP) :: uf(nmodes, nmodes)
      REAL(DP)    :: w2(nmodes)

      INTEGER  :: inu, ierr, iq
      REAL(DP) :: xxq(3)
      LOGICAL  :: mirror_q

      uf = czero
      w2 = zero

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
               !JLB
               IF (ionode) THEN
                  WRITE(stdout, '(5x, "WARNING: Imaginary frequency mode ",&
                  &I6, " at iq=", I6)') inu, iq
               END IF
               wf(inu, iq) = 0.d0
               !wf(inu, iq) = -DSQRT(ABS(w2(inu)))
               !JLB
            ENDIF
         END DO
      END DO
   END SUBROUTINE
   !
   SUBROUTINE init_plrn_random(eigVec)
      USE elph2, ONLY : nktotf
      USE constants_epw, ONLY : ci, cone
      USE epwcom, ONLY : nstate_plrn

      IMPLICIT NONE
      COMPLEX(dp), INTENT(OUT) :: eigVec(:, :)
      REAL(DP),     ALLOCATABLE :: rmat_tmp(:, :)

      CALL RANDOM_SEED()
      ALLOCATE(rmat_tmp(1:nktotf*nbnd_plrn, 1:nstate_plrn))
      CALL RANDOM_NUMBER(rmat_tmp)
      eigVec(1:nktotf*nbnd_plrn, 1:nstate_plrn) = cone*rmat_tmp(1:nktotf*nbnd_plrn, 1:nstate_plrn)
      CALL RANDOM_NUMBER(rmat_tmp)
      eigVec(1:nktotf*nbnd_plrn, 1:nstate_plrn) = eigVec(1:nktotf*nbnd_plrn, 1:nstate_plrn) + &
                                                 &ci*rmat_tmp(1:nktotf*nbnd_plrn, 1:nstate_plrn)
      DEALLOCATE(rmat_tmp)
   END SUBROUTINE
   !
   SUBROUTINE init_plrn_gaussian(r0, xkf_all, k0, eigVec)
      !!
      !! Initialize Ank coefficients with a Gaussian lineshape
      !!
      USE constants_epw, ONLY : czero, cone, ci, twopi, one, zero
      USE epwcom, ONLY : nstate_plrn, init_sigma_plrn
      USE elph2, ONLY : nktotf
      USE cell_base, ONLY : bg, alat

      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: r0(3), k0(3), xkf_all(:, :)
      COMPLEX(dp), INTENT(OUT) :: eigVec(:, :)


      INTEGER :: ibnd, ik, iplrn, ix, iy, iz, indexkn1, indexkn2, ishift
      REAL(dp) :: qcart(3), xxq(3), shift(3), disK
      COMPLEX(dp) :: ctemp

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
         xxq = xkf_all(1:3, ik) - (k0(:)-INT(k0(:))) ! shift k0 to 1BZ
         CALL dgemv('n', 3, 3, one, bg, 3, xxq, 1, zero, qcart, 1)
         ctemp = EXP( -ci* twopi * DOT_PRODUCT( qcart, r0 ))
         disK = -1
         ! Ensure periodicity of Ank checking distance to other equivalent BZ points
         DO ishift = 1, 27
            shift(1:3) = REAL(index_shift(ishift), KIND=DP)
            CALL dgemv('n', 3, 3, one, bg, 3, xxq+shift, 1, zero, qcart, 1)
            disK = MAX(disK, EXP(-init_sigma_plrn * NORM2(qcart)*twopi/alat )) ! for sigma to be in bohr
            !disK = MAX(disK, EXP(-init_sigma_plrn * NORM2(xxq + shift)))
         END DO
         !
         DO ibnd = 1, nbnd_plrn
            indexkn1 = (ik - 1)*nbnd_plrn + ibnd
            eigVec(indexkn1, :) = CONE * disK * ctemp
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE build_plrn_bmat(Bmat, first)
      USE elph2,         ONLY : nkf, nqtotf, wf, xqf, nktotf, etf
      USE epwcom,        ONLY : model_vertex_plrn, nbndsub, debug_plrn
      USE epwcom,        ONLY : nstate_plrn, mixing_Plrn, type_plrn
      USE epwcom,        ONLY : g_start_energy_plrn, g_end_energy_plrn, g_start_band_plrn
      USE epwcom,        ONLY : model_enband_plrn, model_phfreq_plrn, model_vertex_plrn
      USE epwcom,        ONLY : omega_LO_plrn, kappa_plrn, m_eff_plrn
      USE epwcom,        ONLY : io_lvl_plrn
      USE constants_epw, ONLY : czero, one, two, zero, cone, eps2, eps8 !JLB
      USE mp_world,      ONLY : mpime, world_comm
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum
      USE modes,         ONLY : nmodes

      IMPLICIT NONE

      COMPLEX(DP), INTENT(OUT) :: Bmat(:, :)
      LOGICAL, INTENT(IN) :: first
      INTEGER :: iq, ik, ikq, ik_global, ibnd, jbnd, iplrn, indexkn1, indexkn2, indexkn3
      INTEGER :: start_mode, inu, iqpg
      COMPLEX(DP) :: prefac, ctemp
      !JLB
      INTEGER :: jnu, ndegen(nmodes)
      COMPLEX(DP) :: Bmat_tmp(nmodes)
      REAL(DP) :: eig

      Bmat = czero
      !Bmat_comp = czero
      !REWIND(iepfall)
      DO iq = 1, nqtotf
         IF(model_vertex_plrn) THEN
            epf = czero
            epf(1, 1, nmodes, 1:nkf) = gq_model(iq)
         ELSE
            IF(io_lvl_plrn == 0) THEN
               epf(:, :, :, :) = epfall(:, :, :, :, iq)
            ELSE IF (io_lvl_plrn == 1) THEN
               CALL get_buffer(epf, lword_g, iepfall, iq)
            END IF
            !READ(iepfall) epf(1:nbnd_plrn, 1:nbnd_plrn, 1:nmodes, 1:nkf)
         END IF
         IF(test_tags_plrn(1)) THEN
            epf = 1E-6
            epf(:, :, nmodes, :) = 2E-2
         ELSE IF (test_tags_plrn(2)) THEN
            epf = ABS(epf)
         ELSE IF (test_tags_plrn(3)) THEN
             IF(NORM2(xqf(:,iq)) > 1E-5) epf(:, :, nmodes, :) = 0.005/NORM2(xqf(:,iq))
         END IF
         !
         ! energy cutoff for g
         DO ik = 1, nkf
            ik_global = ikqLocal2Global(ik, nktotf)
            DO ibnd = 1, nbnd_g_plrn
               eig = etf_all(ibnd + g_start_band_plrn - 1, ik_global)
               IF (eig < g_start_energy_plrn .OR. eig > g_end_energy_plrn) THEN
                  !print *, "triggered"
                   epf(ibnd, :, :, ik) = czero
                   epf(:, ibnd, :, ik) = czero
               END IF
            END DO
         END DO
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
         IF( iq > iqpg ) THEN
            DO inu = start_mode, nmodes!
               ! Enforce the relation B_q = B*_{G-q}
               !JLB
               Bmat(iq, inu) = CONJG(Bmat(iqpg, inu))
               !JLB
               !Bmat(iq, inu)   = cal_Bmat(iq, inu)
               !Bmat(iqpg, inu) = CONJG(Bmat(iq, inu))
            END DO
         !ELSE IF (iq == iqpg) THEN
         !   DO inu = start_mode, nmodes!
         !      ! Enforce the relation B_q = B*_{G-q}
         !      Bmat(iq, inu)   = REAL(cal_Bmat(iq, inu))
         !   END DO
         ELSE
            DO inu = start_mode, nmodes
               Bmat(iq, inu)   = cal_Bmat(iq, inu)
            END DO
         !JLB
         END IF
         !TODO: to be consistent with Denny, whether this is correct?
         ! IF ( ABS(wf(1, iq)) < eps2 ) Bmat(iq, 1) = czero
      END DO
      ! cal_Bmat only sum over local k, so we have to do mp_sum
      CALL mp_sum( Bmat, inter_pool_comm )
   END SUBROUTINE

   SUBROUTINE build_plrn_hamil(Bmat, Bmat_save, iter)
      USE elph2,         ONLY : nkf, nqtotf, wf, xqf, nktotf, etf
      USE epwcom,        ONLY : model_vertex_plrn, nbndsub
      USE epwcom,        ONLY : nstate_plrn, mixing_Plrn, type_plrn
      USE epwcom,        ONLY : nhblock_plrn, r0_plrn, beta_plrn
      USE epwcom,        ONLY : nqf1, nqf2, nqf3
      USE epwcom,        ONLY : g_start_energy_plrn, g_end_energy_plrn, g_start_band_plrn
      USE epwcom,        ONLY : io_lvl_plrn
      USE constants_epw, ONLY : czero, one, two, zero, cone, eps2, eps8, twopi, ci
      USE mp_world,      ONLY : mpime, world_comm
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum
      USE modes,         ONLY : nmodes
      USE cell_base,     ONLY : at, alat

      IMPLICIT NONE

      !COMPLEX(DP), INTENT(OUT) :: Hamiltonian(:)
      COMPLEX(DP), INTENT(IN) :: Bmat_save(:,:), Bmat(:,:)
      INTEGER, INTENT(IN) :: iter
      INTEGER :: iq, ik, ikq, ik_global, ibnd, jbnd, iplrn
      INTEGER :: indexkn1, indexkn2, indexkn3, indexkn4, indexkn5
      INTEGER :: start_mode, inu, iqpg, index_blk, index_loc, idir, jdir, ialpha
      INTEGER :: ibpi, ibmi, ibpj, ibmj, ivec, ibm, ikpbi, ikmbi, ikpbj, ikmbj
      COMPLEX(DP) :: prefac, ctemp, Q_i, Q_j, Mmn(1:4)
      COMPLEX(DP), ALLOCATABLE :: Bmat_comp(:,:), Hamil_tmp(:, :)
      LOGICAL, ALLOCATABLE :: saved(:)

      REAL(DP) :: F_mat(1:3,1:3), b_vec(1:3), eta(1:3), a_i, a_j, eig

      test_tags_plrn(1) = .false.
      test_tags_plrn(2) = .false.
      test_tags_plrn(3) = .false.



      ! To avoid the zero phonon frequency
      ! IF(wf(1, iq) < eps2) Bmat(iq, 1) = czero

      !IF(.NOT. first) Bmat = mixing_Plrn*Bmat + (1 - mixing_Plrn) * Bmat_save

      ! allocate(saved(1:nkf*nbnd_plrn))
      ! saved(1:nkf*nbnd_plrn) = .false.
      ! call start_clock('H_alloc')


      ! if(.not. mem_save_h) then
      ! allocate(Hamil_tmp(1:lword_h, blocksize))
      ! Hamil_tmp = czero
      ! end if
      ! call stop_clock('H_alloc')

      !
      ! Calculate the Hamiltonian with Bq $$H_{n\bk,n'\bk'} = \delta_{n\bk,n'\bk'}\varepsilon_{n\bk} -\frac{2}{N_p} \sum_{\nu} B^*_{\bq,\nu}g_{nn'\nu}(\bk',\bq)$$
      ! H_{n\bk,n'\bk'} -> Hamil(ik, ibnd, ikq, jbnd)
      ! B^*_{\bq,\nu} -> conj(Bmat(iq, inu))
      ! g_{nn'\nu}(\bk',\bq) -> epf(ibnd, jbnd, inu, ikq, iq)
      ! if q == 0, \delta_{n\bk,n'\bk'}\varepsilon_{n\bk} is diagonal matrix with \varepsilon_{n\bk}

      ! G == (0,0,0) means this is the diagonal term with k=k'
      ! ikq is the global index, i.e. the second index
      ! ik is the local index, i.e. the first index
      Hamil = czero
      !print *, "hblocksize, nhblock_plrn, nkf * nbnd_plrn:", hblocksize, nhblock_plrn, nkf * nbnd_plrn



      DO iq = 1, nqtotf
         IF(model_vertex_plrn) THEN
            epf = czero
            epf(1, 1, nmodes, 1:nkf) = gq_model(iq)
         ELSE
            CALL start_clock('read_gmat') !nbndsub*nbndsub*nmodes*nkf
            IF(io_lvl_plrn == 0) THEN
               epf(:, :, :, :) = epfall(:, :, :, :, iq)
            ELSE IF (io_lvl_plrn == 1) THEN
               CALL get_buffer(epf, lword_g, iepfall, iq)
            END IF
            CALL stop_clock('read_gmat')
         ENDIF
         ! energy cutoff of g
         DO ik = 1, nkf
            ik_global = ikqLocal2Global(ik, nktotf)
            DO ibnd = 1, nbnd_g_plrn
               eig = etf_all(ibnd + g_start_band_plrn - 1, ik_global)
               IF (eig < g_start_energy_plrn .OR. eig > g_end_energy_plrn) THEN
                  !print *, "triggered"
                   epf(ibnd, :, :, ik) = czero
                   epf(:, ibnd, :, ik) = czero
               END IF
            END DO
         END DO

         DO ik = 1, nkf
            ikq = ikq_all(ik, iq)

            DO ibnd = 1, nbnd_plrn
               indexkn1 = (ik - 1)*nbnd_plrn + ibnd

               IF(nhblock_plrn == 1) THEN
                  index_loc = indexkn1
                  index_blk = 1
               ELSE
                  index_loc = MOD(indexkn1 - 1, hblocksize) + 1
                  index_blk = INT((indexkn1 - 1) / hblocksize) + 1
               END IF

               IF(iq /= 1 .and. index_loc == 1 .and. nhblock_plrn /= 1) THEN
                  CALL start_clock('read_Hmat')
                  CALL get_buffer(Hamil, lword_h, ihamil, index_blk)
                  CALL stop_clock('read_Hmat')
               END IF

               CALL start_clock('HdiagTerm')
               IF (isGVec(xqf(1:3, iq))) THEN
                  ! Note that, ik is local index while ikq is global index,
                  ! so even when q=0, ik \= ikq, but ik_global == ikq
                  ! delta_{nn' kk'} epsilon_{nk}
                  ctemp = etf_all(select_bands_plrn(ibnd), ikq)
                  indexkn2 = (ikq - 1)*nbnd_plrn + ibnd
                  Hamil(indexkn2, index_loc) = Hamil(indexkn2, index_loc) + ctemp
               END IF
               CALL stop_clock('HdiagTerm')

               CALL start_clock('HOffDiagTerm')
               DO jbnd = 1, nbnd_plrn
                  indexkn2 = (ikq - 1) * nbnd_plrn + jbnd
                  DO inu = 1, nmodes
                     ctemp = type_plrn * two/REAL(nqtotf,KIND=dp)*(Bmat(iq, inu))*&
                        CONJG(epf(select_bands_plrn(jbnd) - g_start_band_plrn + 1, &
                              select_bands_plrn(ibnd)- g_start_band_plrn + 1, inu, ik))
                     Hamil(indexkn2, index_loc) = Hamil(indexkn2, index_loc) + ctemp
                  END DO
               END DO
               CALL stop_clock('HOffDiagTerm')
               IF(nhblock_plrn /= 1) THEN
                  IF((index_loc == hblocksize .or. indexkn1 == nkf*nbnd_plrn)) THEN
                     CALL start_clock('Write_Hmat')
                     CALL save_buffer(Hamil, lword_h, ihamil, index_blk)
                     Hamil = czero
                     CALL stop_clock('Write_Hmat')
                  END IF
               END IF
            END DO !ibnd
         END DO !ik
      END DO ! iq
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   ! Return true if xxk integer times of the reciprocal vector
   ! if xxk is the difference of two vectors, then return true if these
   ! two vector are the same
   FUNCTION isGVec(xxk)
      USE constants_epw, ONLY : eps6

      IMPLICIT NONE
      LOGICAL :: isGVec
      REAL(dp), INTENT(IN) :: xxk(3)
      isGVec = &
         ABS(xxk(1) - NINT(xxk(1))) < eps6 .and. &
         ABS(xxk(2) - NINT(xxk(2))) < eps6 .and. &
         ABS(xxk(3) - NINT(xxk(3))) < eps6
   END FUNCTION
   !
   !-----------------------------------------------------------------------
   ! Return the global index of the local k point ik
   FUNCTION ikqLocal2Global(ikq, nkqtotf)
      USE division,       ONLY : fkbounds

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ikq, nkqtotf ! ik or iq
      INTEGER :: ikqLocal2Global
      INTEGER :: startn, lastn

      CALL start_clock('ik_l2g')
      CALL fkbounds(nkqtotf, startn, lastn)

      ikqLocal2Global = startn + ikq - 1
      IF (ikqLocal2Global > lastn) THEN
         CALL errore('ikqLocal2Global', 'Index of k/q is beyond this pool.', 1)
      END IF
      CALL stop_clock('ik_l2g')

   END FUNCTION
   !
   !-----------------------------------------------------------------------
   ! Return the global index of the local k point ik
   FUNCTION ikGlobal2Local(ik_g, nktotf)
      USE division,       ONLY : fkbounds

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ik_g, nktotf ! ik or iq
      INTEGER :: ikGlobal2Local
      INTEGER :: startn, lastn

      CALL fkbounds(nktotf, startn, lastn)

      ikGlobal2Local = ik_g - startn + 1

      IF(ikGlobal2Local <= 0) THEN
         ikGlobal2Local = 0
      END IF
   END FUNCTION
   !
   !-----------------------------------------------------------------------
   ! Return EXP(-(energy/sigma)**2)
   SUBROUTINE cal_f_delta(energy, sigma, f_delta)
      USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero

      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: energy(:), sigma
      REAL(dp), INTENT(OUT) :: f_delta(:)

      f_delta = EXP(-(energy/sigma)**2)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   ! Calculate Hpsi with psi as input to use the diagon sovler in KS_solver cegterg
   ! cegterg take two external subroutine to calculate Hpsi and Spsi to calculate
   ! ( H - e S ) * evc = 0, since H and S is not saved due to their sizes
   ! Hamil need to be passed to h_psi because the parameter space is fixed
   ! to meet the requirement of Davidson diagonalization.
   SUBROUTINE h_psi_plrn(lda, n, m, psi, hpsi)
      USE elph2,         ONLY : nkf, nqtotf, wf, xqf, nktotf, etf
      USE epwcom,        ONLY : model_vertex_plrn, time_rev_A_plrn, nbndsub
      USE epwcom,        ONLY : mixing_Plrn, type_plrn, nhblock_plrn
      USE epwcom,        ONLY : beta_plrn
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
      INTEGER, INTENT(IN) :: lda
      !! leading dimension of arrays psi, spsi, hpsi, which is nkf * nbnd_plrn
      INTEGER, INTENT(IN) :: n
      !! true dimension of psi, spsi, hpsi
      INTEGER, INTENT(IN) :: m
      !! number of states psi
      COMPLEX(DP), INTENT(INOUT) :: psi(lda,m)
      !! the wavefunction
      COMPLEX(DP), INTENT(OUT) :: hpsi(lda,m)

      INTEGER :: iq, ik, ikq, ikpg, ik_global, ibnd, jbnd, iplrn, indexkn1, indexkn2, indexkn3
      INTEGER :: start_mode, inu, iqpg, startn, lastn, index_loc, index_blk
      INTEGER :: ibp, ibm, ikpb, ikmb, ivec
      COMPLEX(DP) :: prefac, ctemp, ctemp2, hamil_kq
      COMPLEX(DP), ALLOCATABLE :: hamiltonian(:), hpsi_global(:, :) ! One row of Hamiltonian

      ! Gather psi (dimension nkf) to form eigVec (dimension nktotf)
      CALL start_clock('cal_hpsi')
      IF(lda < nkf * nbnd_plrn) CALL errore('h_psi_plrn', 'leading dimension of arrays psi is not correct', 1)
      eigVec = czero
      DO ik = 1, nkf
         ik_global = ikqLocal2Global(ik, nktotf)
         DO ibnd = 1, nbnd_plrn
            indexkn1 = (ik-1)*nbnd_plrn + ibnd
            indexkn2 = (ik_global-1)*nbnd_plrn + ibnd
            eigVec(indexkn2, 1:m) = psi(indexkn1, 1:m)
         END DO
      END DO
      CALL mp_sum(eigVec, inter_pool_comm)
      !

      ! Iterative diagonalization only get the lowest eigenvalues,
      ! however, we will need the highest eigenvalues if we are calculating hole polaron
      ! so, we multiply hpsi by -1, and get eigenvalues in diagonalization
      ! and then multiply eigenvalues by -1.
      hpsi(1:lda,1:m) = czero
      DO ik = 1, nkf
         DO ibnd = 1, nbnd_plrn
            indexkn1 = (ik - 1) * nbnd_plrn + ibnd
            IF(nhblock_plrn == 1) THEN
               index_loc = indexkn1
               index_blk = 1
            ELSE
               index_loc = MOD(indexkn1 - 1, hblocksize) + 1
               index_blk = INT((indexkn1 - 1) / hblocksize) + 1
            END IF
            IF (index_loc == 1 .and. nhblock_plrn /= 1) CALL get_buffer(Hamil, lword_h, ihamil, index_blk)

            DO ikq = 1, nktotf
               DO jbnd = 1, nbnd_plrn
                  indexkn2 = (ikq - 1)*nbnd_plrn + jbnd
                  hpsi(indexkn1, 1:m) = hpsi(indexkn1, 1:m) - &
                     type_plrn * Hamil(indexkn2, index_loc) * eigVec(indexkn2, 1:m)
               END DO
            END DO
         END DO
      END DO

      CALL stop_clock('cal_hpsi')
   END SUBROUTINE
   !-----------------------------------------------------------------------
   SUBROUTINE s_psi_plrn(lda, n, m, psi, spsi)
      IMPLICIT NONE

      INTEGER,     INTENT(IN) :: lda, n, m
      COMPLEX(DP), INTENT(IN) :: psi(lda,m), spsi(lda,m)

      CALL errore('s_psi_plrn',"WARNING: This function should not be called at all!", 1)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE g_psi_plrn( lda, n, m, npol, psi, e )
      !! This routine computes an estimate of the inverse Hamiltonian
      !! and applies it to m wavefunctions.
      IMPLICIT NONE

      INTEGER     :: lda, n, m, npol
      COMPLEX(DP) :: psi(lda, npol, m)
      REAL(DP)    :: e(m)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE get_cfac(xk, nrr_k, ndegen_k, irvec_r, dims, cfac)
      USE epwcom, ONLY : use_ws
      USE constants_epw, ONLY : twopi, ci, czero
      USE kinds,         ONLY : dp, i4b

      IMPLICIT NONE


      INTEGER, INTENT(IN):: nrr_k, dims
      INTEGER, INTENT(IN):: ndegen_k(nrr_k, dims, dims)
      REAL(KIND=dp), INTENT (IN) :: xk(3), irvec_r(3, nrr_k)
      COMPLEX(KIND=dp), INTENT(OUT) :: cfac(nrr_k, dims, dims)
      ! Local Variables
      REAL(KIND=dp) :: rdotk(nrr_k)
      INTEGER:: ikk, ikq, iw, iw2, ir

      cfac = czero
      rdotk = czero

      CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xk, 1, 0.0_dp, rdotk, 1 )
      !
      IF (use_ws) THEN
         DO iw=1, dims
            DO iw2=1, dims
               DO ir = 1, nrr_k
                  IF (ndegen_k(ir,iw2,iw) > 0) THEN
                     cfac(ir,iw2,iw)  = EXP( ci*rdotk(ir) ) / ndegen_k(ir,iw2,iw)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSE
         cfac(:,1,1)   = EXP( ci*rdotk(:) ) / ndegen_k(:,1,1)
      ENDIF
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   ! B_{qu} = 1/N_p \sum_{mnk} A^*_{mk+q}A_{nk} [g_{mnu}(k,q)/\hbar\omega_{qu}]
   FUNCTION cal_Bmat(iq, inu)
      USE elph2,  ONLY : nkf, nktotf, wf, nqtotf
      USE epwcom, ONLY : nstate_plrn, model_vertex_plrn
      USE epwcom, ONLY : g_start_band_plrn
      USE constants_epw, ONLY : czero, one, eps2, cone, eps8
      USE mp_world,      ONLY : mpime, world_comm

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iq, inu
      INTEGER :: ik, ikq, ik_global, ibnd, jbnd, iplrn, indexkn1, indexkn2
      COMPLEX(DP) :: cal_Bmat
      COMPLEX(DP) :: prefac

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
                  indexkn1 = (ikq - 1)*nbnd_plrn + ibnd
                  indexkn2 = (ik_global - 1) * nbnd_plrn + jbnd
                  IF (wf(inu, iq) > eps8 ) THEN
                     prefac = cone/(wf(inu, iq) * REAL(nqtotf, dp))
                  ELSE
                     prefac = czero
                  END IF
                  ! B_{q\nu} = \frac{1}{N_p}\sum_{nn'k}A^*_{n'k+q}\frac{g_{n'n\nu}(k, q)}{\hbar \omega_{q\nu}} A_{nk}
                  cal_Bmat = cal_Bmat + prefac * (eigVec(indexkn2, iplrn)) * CONJG(eigVec(indexkn1, iplrn)) * &
                     (epf(select_bands_plrn(ibnd) - g_start_band_plrn + 1, &
                     select_bands_plrn(jbnd) - g_start_band_plrn + 1, inu, ik)) !conjg
               END DO
            END DO
         END DO
      END DO

      !JLB - discard zero or imaginary frequency modes
      IF (wf(inu, iq) < eps8) THEN
         cal_Bmat = czero
      END IF
      !JLB

   END FUNCTION
   !
   !-----------------------------------------------------------------------
   ! Find the index of Gamma point i.e. (0, 0, 0) in xkf_all
   ! which contains all the crystal coordinates of the k/q points
   ! if Gamma point is not included, return 0
   FUNCTION indexGamma(xkf_all)
      USE elph2,        ONLY : nkf, nktotf
      USE mp,           ONLY : mp_sum
      USE mp_global,    ONLY : inter_pool_comm

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: xkf_all(:,:)
      INTEGER :: indexGamma, ik, ik_global

      indexGamma = 0

      DO ik = 1, nkf
         ik_global = ikqLocal2Global(ik, nktotf)
         IF(isGVec(xkf_all(1:3, ik_global))) THEN
            indexGamma = ik_global
         END IF
      END DO
      CALL mp_sum(indexGamma, inter_pool_comm)

      IF(.NOT. isGVec(xkf_all(1:3, indexGamma))) THEN
         CALL errore('indexGamma','The index of Gamma point is wrong!', 1)
      END IF

   END FUNCTION
   !
   SUBROUTINE norm_plrn_wf(eigVec, norm_new)
      USE elph2,  ONLY : nkf, nqtotf, wf, nktotf
      USE epwcom, ONLY : nstate_plrn, time_rev_A_plrn
      USE constants_epw, ONLY : czero, one, two, cone
      USE mp,           ONLY : mp_sum
      USE mp_global,    ONLY : inter_pool_comm

      IMPLICIT NONE

      COMPLEX(DP), INTENT(INOUT) :: eigVec(:, :)
      REAL(DP), INTENT(IN) :: norm_new
      REAL(DP) :: norm
      INTEGER :: iplrn


      DO iplrn = 1, nstate_plrn
         norm = REAL(DOT_PRODUCT(eigVec(1:nbnd_plrn*nktotf, iplrn), eigVec(1:nbnd_plrn*nktotf, iplrn)))
         eigVec(:, iplrn) = eigVec(:, iplrn)/DSQRT(norm) * SQRT(norm_new)
      END DO

   END SUBROUTINE
   !-----------------------------------------------------------------------
   SUBROUTINE check_time_rev_sym(eigVec)
      USE elph2,  ONLY : nkf, nqtotf, wf, nktotf
      USE epwcom, ONLY : nstate_plrn, time_rev_A_plrn
      USE constants_epw, ONLY : czero, one, two, cone
      USE mp,           ONLY : mp_sum
      USE mp_global,    ONLY : inter_pool_comm

      IMPLICIT NONE

      INTEGER :: iq, inu, ik, ikq, ik_global, ibnd, iplrn, ikpg, indexkn1, indexkn2
      INTEGER :: nPlrn_l
      COMPLEX(DP) :: temp
      COMPLEX(DP), INTENT(INOUT) :: eigVec(:, :)
      COMPLEX(DP), ALLOCATABLE :: eigVec_save(:, :)
      REAL(DP) :: norm


      ! nstate_plrn
      nPlrn_l = 1
      ALLOCATE(eigVec_save(nktotf*nbnd_plrn, nPlrn_l))
      eigVec_save = czero

      DO ik = 1, nkf
         ik_global = ikqLocal2Global(ik, nktotf)
         ikpg = kpg_map(ik_global)
         !if(is_mirror_k(ik)) then
         DO ibnd = 1, nbnd_plrn
            indexkn1 = (ikpg - 1)*nbnd_plrn + ibnd
            indexkn2 = (ik_global - 1)*nbnd_plrn + ibnd
            eigVec_save(indexkn1, 1:nPlrn_l)  = CONJG(eigVec(indexkn2, 1:nPlrn_l))
         END DO
         !end if
      END DO
      CALL mp_sum(eigVec_save, inter_pool_comm)
      eigVec(:, 1:nPlrn_l) = (eigVec(:, 1:nPlrn_l) + eigVec_save(:, 1:nPlrn_l))

      DO iplrn = 1, nPlrn_l
         norm = REAL(DOT_PRODUCT(eigVec(1:nbnd_plrn*nktotf, iplrn), eigVec(1:nbnd_plrn*nktotf, iplrn)))!nktotf*nbnd_plrn*
         eigVec(:, iplrn) = eigVec(:, iplrn)/DSQRT(norm)
      END DO

      DEALLOCATE(eigVec_save)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE diag_serial(estmteRt, eigVec)
      USE constants_epw, ONLY : czero, twopi, ci, cone, zero
      USE elph2,         ONLY : nkf, nqtotf, nktotf, xkf
      USE elph2,         ONLY : etf, chw
      USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3, type_plrn, nhblock_plrn
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp_world,      ONLY : world_comm
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum, mp_bcast

      IMPLICIT NONE

      REAL(DP),    INTENT(OUT) :: estmteRt(:)
      COMPLEX(DP), INTENT(OUT) :: eigVec(:, :)

      REAL(DP)    :: rtemp, xxk(3), shift(3)
      COMPLEX(DP) :: ctemp
      INTEGER     :: iq, inu, ik, ikk, ikq, ik_global, iplrn, ikpg, icount
      INTEGER     :: ibnd, jbnd, indexkn1, indexkn2
      INTEGER     :: lwork, info, mm, index_loc, index_blk

      INTEGER,      ALLOCATABLE :: iwork(:), ifail(:)
      REAL(DP),     ALLOCATABLE :: rwork(:)
      COMPLEX(DP),  ALLOCATABLE :: work(:)

      COMPLEX(DP),  ALLOCATABLE :: Hamil_save(:,:)
      COMPLEX(DP),  ALLOCATABLE :: Identity(:,:)


      ALLOCATE(Hamil_save(nktotf*nbnd_plrn, nktotf*nbnd_plrn))
      Hamil_save = czero
      DO ik = 1, nkf
         ik_global = ikqLocal2Global(ik, nktotf)
         DO ibnd = 1, nbnd_plrn
            indexkn1 = (ik-1)*nbnd_plrn + ibnd

            index_loc = MOD(indexkn1 - 1, hblocksize) + 1
            index_blk = INT((indexkn1 - 1)/hblocksize) + 1
            IF (index_loc == 1 .and. nhblock_plrn /= 1) CALL get_buffer(Hamil, lword_h, ihamil, index_blk)

            indexkn2 = (ik_global-1)*nbnd_plrn + ibnd
            Hamil_save(indexkn2, 1:nktotf*nbnd_plrn) = - type_plrn * Hamil(1:nktotf*nbnd_plrn, index_loc)
         END DO
      END DO
      CALL mp_sum(Hamil_save, inter_pool_comm)
      IF (ionode) THEN
         ALLOCATE(Identity(nktotf*nbnd_plrn, nktotf*nbnd_plrn))
         lwork = 5*nktotf*nbnd_plrn
         ALLOCATE( rwork( 7*nktotf*nbnd_plrn ) )
         ALLOCATE( iwork( 5*nktotf*nbnd_plrn ) )
         ALLOCATE( ifail( nktotf*nbnd_plrn ) )
         ALLOCATE( work( lwork ) )
         Identity = czero

         DO ibnd = 1, nbnd_plrn*nktotf
            Identity(ibnd, ibnd) = cone
         END DO

         eigVec = czero
         estmteRt = zero

         CALL ZHEGVX( 1, 'V', 'I', 'U', nktotf*nbnd_plrn, Hamil_save, nktotf*nbnd_plrn, Identity, nktotf*nbnd_plrn, &
            zero, zero, 1, nstate_plrn, zero, mm, estmteRt(1:nstate_plrn), eigVec, nktotf*nbnd_plrn, &
            work, lwork, rwork, iwork, ifail, info )

         IF (info /= 0) CALL errore('diag_serial','Polaron: diagonal error.', 1)
         DEALLOCATE(rwork, iwork, ifail, work, Identity)
      END IF
      DEALLOCATE(Hamil_save)
      CALL mp_bcast( estmteRt, meta_ionode_id, world_comm )
      CALL mp_bcast( eigVec, meta_ionode_id, world_comm )
   END SUBROUTINE
   !-----------------------------------------------------------------------
   !
   SUBROUTINE diag_parallel(estmteRt, eigVec)
      USE constants_epw, ONLY : czero, twopi, ci
      USE constants_epw, ONLY : eps5, eps6, eps4, eps2, eps8, eps10
      USE elph2,         ONLY : nkf, nqtotf, nktotf, xkf
      USE elph2,         ONLY : etf, chw
      USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3
      USE epwcom,        ONLY : ethrdg_plrn
      USE epwcom,        ONLY : adapt_ethrdg_plrn, init_ethrdg_plrn, nethrdg_plrn
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp_world,      ONLY : world_comm
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum, mp_bcast, mp_size, mp_max
      USE mp_bands, ONLY : intra_bgrp_comm, inter_bgrp_comm, mp_start_bands
      USE mp_bands_util, ONLY : intra_bgrp_comm_ => intra_bgrp_comm
      USE mp_bands_util, ONLY : inter_bgrp_comm_ => inter_bgrp_comm

      IMPLICIT NONE

      REAL(DP),    INTENT(OUT) :: estmteRt(:)
      COMPLEX(DP), INTENT(OUT) :: eigVec(:, :)

      COMPLEX(DP), ALLOCATABLE :: Identity(:, :)
      REAL(DP)    :: rtemp, xxk(3), shift(3)
      REAL(DP)    :: ethrdg_init
      COMPLEX(DP) :: ctemp
      INTEGER     :: iq, inu, ik, ikk, ikq, ik_global, iplrn, ikpg, icount
      INTEGER     :: ibnd, jbnd, itemp, jtemp
      INTEGER     :: indexkn1, indexkn2
      INTEGER :: ithr, nthr

      INTEGER :: npw, npwx, dav_iter, notcnv, btype(nstate_plrn), nhpsi

      INTEGER,      ALLOCATABLE :: iwork(:), ifail(:)
      REAL(DP),     ALLOCATABLE :: rwork(:)
      COMPLEX(DP),  ALLOCATABLE :: work(:)
      COMPLEX(DP),  ALLOCATABLE :: psi(:, :)
      REAL(DP)                  :: ethrdg



      npw = nkf * nbnd_plrn

      npwx = npw
      CALL mp_max(npwx, inter_pool_comm)

      ALLOCATE(psi(1:npwx, 1:nstate_plrn))

      ! JLB: Option for adaptive threshold
      IF (adapt_ethrdg_plrn) THEN
         ethrdg_init = init_ethrdg_plrn
         nthr = nethrdg_plrn
         IF(ionode) THEN 
            WRITE(stdout, "(a)") "     Adaptive threshold on iterative diagonalization activated:"
            WRITE(stdout, "(a)") "     threshold, # iterations, eigenvalue(Ry)"
         END IF 
      ELSE
         nthr = 1
      END IF
      !
      DO ithr = 1, nthr

         psi = czero
         btype(1:nstate_plrn) = 1
         !
         IF (adapt_ethrdg_plrn) THEN
            ethrdg = 10**(LOG10(ethrdg_init) + (ithr-1)*(LOG10(ethrdg_plrn) - LOG10(ethrdg_init))/(nthr-1))
         ELSE
            ethrdg = ethrdg_plrn
         END IF
         !
         ! split eigVector (nqtotf) into parallel pieces psi (nkf), contains corresponding part with Hpsi
         DO ik = 1, nkf
            ik_global = ikqLocal2Global(ik, nktotf)
            DO ibnd = 1, nbnd_plrn
               indexkn1 = (ik-1)*nbnd_plrn + ibnd
               indexkn2 = (ik_global-1)*nbnd_plrn + ibnd
               psi(indexkn1, 1:nstate_plrn) = eigVec(indexkn2, 1:nstate_plrn)
            END DO
         END DO
         ! inter_bgrp_comm should be some non-existing number,
         ! to make the nodes in bgrp equal to 1
         ! intra_bgrp_comm is parallel PW in pwscf
         ! but here it should be parallel K.
         ! Save them before change them
         itemp = intra_bgrp_comm_
         jtemp = inter_bgrp_comm_
         !call mp_start_bands(1, world_comm)
         !
         intra_bgrp_comm_ = inter_pool_comm
         inter_bgrp_comm_ = inter_bgrp_comm
         !
         !write(stdout, *) "test communicator", inter_pool_comm, inter_bgrp_comm, intra_bgrp_comm
         !write(stdout, *) mp_size(inter_pool_comm)
         !write(stdout, *) "past test inter_pool communicator"
         !write(stdout, *) mp_size(inter_bgrp_comm)
         !write(stdout, *) "past test inter_bgrp communicator"
         !write(stdout, *) mp_size(intra_bgrp_comm)
         !write(stdout, *) "past test intra_bgrp communicator" 4*
         CALL start_clock('cegterg_prln')
         CALL cegterg( h_psi_plrn, s_psi_plrn, .false., g_psi_plrn, &
            npw, npwx, nstate_plrn, nstate_plrn*10, 1, psi, ethrdg, &
            estmteRt, btype, notcnv, .false., dav_iter, nhpsi)
         CALL start_clock('cegterg_prln')
         IF(adapt_ethrdg_plrn .AND. ionode) WRITE(stdout, "(a, E14.6, I6, E14.6)") "   ", ethrdg, dav_iter, estmteRt
         IF(notcnv>0 .AND. ionode) WRITE(stdout, "(a)") "   WARNING: Some eigenvalues not converged, & 
         check initialization, ethrdg_plrn or try adapt_ethrdg_plrn"
         !
         intra_bgrp_comm_ = itemp
         inter_bgrp_comm_ = jtemp
         !
         !CALL gatherVector()
         eigVec = czero
         DO ik = 1, nkf
            ik_global = ikqLocal2Global(ik, nktotf)
            DO ibnd = 1, nbnd_plrn
               indexkn1 = (ik-1)*nbnd_plrn + ibnd
               indexkn2 = (ik_global-1)*nbnd_plrn + ibnd
               eigVec(indexkn2, 1:nstate_plrn) = psi(indexkn1, 1:nstate_plrn)
            END DO
         END DO
         CALL mp_sum(eigVec, inter_pool_comm)

      END DO
      !
      DEALLOCATE(psi)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   ! Write ionic positions and displacements
   SUBROUTINE write_plrn_dtau_xsf(dtau, nqf1, nqf2, nqf3, filename, species)
      USE constants_epw, ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
      !USE elph2,         ONLY : nkf, nqtotf
      USE epwcom,        ONLY : nstate_plrn!, nqf1, nqf2, nqf3
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp_world,      ONLY : world_comm
      USE mp,            ONLY : mp_sum, mp_bcast
      USE modes,         ONLY : nmodes
      USE ions_base,     ONLY : nat, amass, ityp, tau, atm, nsp, na, ntypx
      USE cell_base,     ONLY : at, alat

      IMPLICIT NONE

      COMPLEX(DP), INTENT(IN) :: dtau(:, :)
      INTEGER, INTENT(IN) :: nqf1, nqf2, nqf3
      CHARACTER(LEN=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN), OPTIONAL :: species(50)

      REAL(DP)    :: rtemp, cell(3, 3), shift(1:3)
      INTEGER     :: wan_func_file, indexkn1, nbnd_out, nat_all, nptotf, nqf_s(1:3)
      INTEGER     :: ix, iy, iz, iRp_local, iRp, iatm, idir, iatm_all, iatm_sp, ika
      INTEGER     :: iq, inu, ik, ikq, ik_global, ibnd, iplrn, ikpg, isp, Rp_vec(1:3)

      INTEGER,  ALLOCATABLE   :: elements(:)
      REAL(DP), ALLOCATABLE   :: atoms(:,:), displacements(:,:)
      !CHARACTER(LEN=3), ALLOCATABLE  :: elements( : )


      ! total number of atoms is
      ! (number of atoms in unit cell) x (number of cells in the supercell)
      nptotf =  nqf1 * nqf2 * nqf3
      nqf_s  = (/nqf1, nqf2, nqf3/)
      nat_all = nat * nptotf

      ALLOCATE(atoms(3, nat_all))
      ALLOCATE(elements(nat_all))
      ALLOCATE(displacements(3, nat_all))

      atoms = zero
      elements = 0
      displacements = zero

      cell(1:3, 1) = at(1:3, 1) * nqf1
      cell(1:3, 2) = at(1:3, 2) * nqf2
      cell(1:3, 3) = at(1:3, 3) * nqf3


      !print *, "nsp, ityp(1) ", nsp, ityp(1) ! nsp = 0, cannot use it!
      !print *, "at 1:", at(1:3, 1)
      !print *, "at 2:", at(1:3, 2)
      !print *, "at 3:", at(1:3, 3)
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
                  END IF
                  atoms(1:3, iatm_all) = tau(1:3, iatm) + shift(1:3)
                  displacements(1:3, iatm_all) = REAL(dtau(iRp, ika:ika+2))
               END IF
            END DO
         END DO
      END DO


      !displacements = zero

      cell = cell * alat
      atoms = atoms * alat

      CALL write_xsf_file(filename, cell*bohr2ang, elements, atoms*bohr2ang, displacements*bohr2ang)

      DEALLOCATE(atoms, elements, displacements)
   END SUBROUTINE
   !-----------------------------------------------------------------------
   ! JLB: Write ionic positions and displacements for transformed supercell
   SUBROUTINE scell_write_plrn_dtau_xsf(dtau, nqtotf_p, nRp_p, Rp_p, as_p, filename, species)
      USE constants_epw, ONLY : czero, ryd2ev, ryd2mev, zero, bohr2ang
      !USE elph2,         ONLY : nkf, nqtotf
      USE epwcom,        ONLY : nstate_plrn!, nqf1, nqf2, nqf3
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp_world,      ONLY : world_comm
      USE mp,            ONLY : mp_sum, mp_bcast
      USE modes,         ONLY : nmodes
      USE ions_base,     ONLY : nat, amass, ityp, tau, atm, nsp, na, ntypx
      USE cell_base,     ONLY : at, alat

      IMPLICIT NONE

      COMPLEX(DP), INTENT(IN) :: dtau(:, :)
      INTEGER, INTENT(IN) :: nqtotf_p, nRp_p, Rp_p(:,:)
      REAL(DP), INTENT(IN) :: as_p(3,3)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN), OPTIONAL :: species(50)

      REAL(DP)    :: rtemp, cell(3, 3), shift(1:3)
      INTEGER     :: wan_func_file, indexkn1, nbnd_out, nat_all, nqf_s(1:3)
      INTEGER     :: ix, iy, iz, iRp_local, iRp, iatm, idir, iatm_all, iatm_sp, ika
      INTEGER     :: iq, inu, ik, ikq, ik_global, ibnd, iplrn, ikpg, isp, Rp_vec(1:3)

      INTEGER,  ALLOCATABLE   :: elements(:)
      REAL(DP), ALLOCATABLE   :: atoms(:,:), displacements(:,:)
      !CHARACTER(LEN=3), ALLOCATABLE  :: elements( : )


      ! total number of atoms is
      ! (number of atoms in unit cell) x (number of cells in the supercell)
      nat_all = nat * nRp_p

      ALLOCATE(atoms(3, nat_all))
      ALLOCATE(elements(nat_all))
      ALLOCATE(displacements(3, nat_all))

      atoms = zero
      elements = 0
      displacements = zero

      cell(1:3, 1) = as_p(1, 1:3)
      cell(1:3, 2) = as_p(2, 1:3)
      cell(1:3, 3) = as_p(3, 1:3)

      !print *, "nsp, ityp(1) ", nsp, ityp(1) ! nsp = 0, cannot use it!
      !print *, "at 1:", at(1:3, 1)
      !print *, "at 2:", at(1:3, 2)
      !print *, "at 3:", at(1:3, 3)
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
                  END IF
                  atoms(1:3, iatm_all) = tau(1:3, iatm) + shift(1:3)
                  displacements(1:3, iatm_all) = REAL(dtau(iRp, ika:ika+2))
               END IF
            END DO
         END DO
      END DO


      !displacements = zero

      cell = cell * alat
      atoms = atoms * alat

      CALL write_xsf_file(filename, cell*bohr2ang, elements, atoms*bohr2ang, displacements*bohr2ang)

      DEALLOCATE(atoms, elements, displacements)
   END SUBROUTINE
   !-----------------------------------------------------------------------
   SUBROUTINE write_xsf_file(filename, cell, elements, atoms, forces, data_cube)

      IMPLICIT NONE

      INTEGER, INTENT(IN)              :: elements(:)
      REAL(dp), INTENT(IN)             :: cell(3, 3), atoms(:, :)
      REAL(dp), INTENT(IN) , OPTIONAL  :: forces(:, :), data_cube(:, :, :)
      CHARACTER(LEN=*), INTENT(IN):: filename

      REAL(DP)    :: rtemp
      INTEGER     :: file_unit, indexkn1, nbnd_out
      INTEGER     :: iq, inu, ik, ikq, ik_global, ibnd, iplrn, ikpg
      INTEGER     :: ix, iy, iz, nx, ny, nz, iatm, natm, shapeTemp(3)

      file_unit = 602
      natm = UBOUND(elements, DIM=1)

      OPEN (UNIT=file_unit, FILE=TRIM(filename), FORM='formatted', STATUS='unknown')

      WRITE (file_unit, *) '#'
      WRITE (file_unit, *) '# Generated by the EPW polaron code'
      WRITE (file_unit, *) '#'
      WRITE (file_unit, *) '#'
      WRITE (file_unit, *) 'CRYSTAL'
      WRITE (file_unit, *) 'PRIMVEC'
      WRITE (file_unit, '(3f12.7)') cell(1:3, 1)
      WRITE (file_unit, '(3f12.7)') cell(1:3, 2)
      WRITE (file_unit, '(3f12.7)') cell(1:3, 3)

      WRITE (file_unit, *) 'PRIMCOORD'
      ! The second number is always 1 for PRIMCOORD coordinates,
      ! according to http://www.xcrysden.org/doc/XSF.html
      WRITE (file_unit, *)  natm, 1

      DO iatm = 1, natm
         IF(PRESENT(forces)) THEN
            WRITE(file_unit,'(I3, 3x, 3f15.9, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm), forces(1:3, iatm)
         ELSE
            WRITE(file_unit,'(I3, 3x, 3f15.9)') elements(iatm), atoms(1:3, iatm)
         END IF
      END DO

      IF(PRESENT(data_cube)) THEN
         shapeTemp = SHAPE(data_cube)
         WRITE (file_unit, '(/)')
         WRITE (file_unit, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
         WRITE (file_unit, '(3i6)') SHAPE(data_cube)
         WRITE (file_unit, '(3f12.6)') 0.0, 0.0, 0.0
         WRITE (file_unit, '(3f12.7)') cell(1:3, 1)
         WRITE (file_unit, '(3f12.7)') cell(1:3, 2)
         WRITE (file_unit, '(3f12.7)') cell(1:3, 3)
         ! TODO: data cube is probably to large to take in the same way of lattice information
         ! May be usefull and implemented in the furture
         WRITE (file_unit, *) (((data_cube(ix, iy, iz), ix = 1, shapeTemp(1)), &
            iy = 1, shapeTemp(2)), iz = 1, shapeTemp(3))
         WRITE (file_unit, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
      END IF
      CLOSE (file_unit)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE write_plrn_wf(eigvec_wan, filename, etf_all)
      USE constants_epw, ONLY : czero, ryd2ev, ryd2mev
      USE elph2,         ONLY : nkf, nqtotf, nktotf
      USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp_world,      ONLY : world_comm
      USE mp,            ONLY : mp_sum, mp_bcast
      USE epwcom,        ONLY : nbndsub, scell_mat_plrn

      IMPLICIT NONE
      COMPLEX(DP), INTENT(IN) :: eigvec_wan(:, :)
      REAL(DP), INTENT(IN), OPTIONAL :: etf_all(:, :)
      CHARACTER(LEN=*), INTENT(IN) :: filename

      REAL(DP)    :: rtemp
      INTEGER     :: wan_func_file, indexkn1, nbnd_out
      INTEGER     :: iq, inu, ik, ikq, ik_global, ibnd, iplrn, ikpg

      IF(PRESENT(etf_all)) THEN
         nbnd_out = nbnd_plrn
      ELSE
         nbnd_out = nbndsub
      END IF
      wan_func_file = 602

      OPEN(UNIT = wan_func_file, FILE = TRIM(filename))
      IF (scell_mat_plrn) THEN
         WRITE(wan_func_file, '(a, 3I10)') 'Scell', nktotf, nbndsub, nstate_plrn
      ELSE
         WRITE(wan_func_file, '(6I10)') nkf1, nkf2, nkf3, nktotf, nbndsub, nstate_plrn
      END IF

      DO ik = 1, nktotf
         DO ibnd = 1, nbnd_out
            DO iplrn = 1, nstate_plrn
               indexkn1 = (ik-1)*nbnd_out + ibnd
               IF(PRESENT(etf_all)) THEN
                  WRITE(wan_func_file, '(2I5, 4f15.7)') ik, ibnd, etf_all(select_bands_plrn(ibnd), ik)*ryd2ev, &
                     eigvec_wan(indexkn1, iplrn), ABS(eigvec_wan(indexkn1, iplrn))
               ELSE
                  WRITE(wan_func_file, '(2f15.7)') eigvec_wan(indexkn1, iplrn)
               END IF
            END DO
         END DO
      END DO
      CLOSE(wan_func_file)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   SUBROUTINE read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, filename, scell, etf_all)
      USE constants_epw, ONLY : czero
      USE elph2,         ONLY : nkf, nqtotf, nktotf
      USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp_world,      ONLY : world_comm
      USE mp,            ONLY : mp_sum, mp_bcast

      IMPLICIT NONE
      COMPLEX(DP), ALLOCATABLE, INTENT(INOUT) :: eigvec_wan(:, :)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      LOGICAL, INTENT(IN), OPTIONAL  :: scell
      REAL(DP), INTENT(IN), OPTIONAL :: etf_all(:, :)
      INTEGER, INTENT(OUT) :: nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p

      REAL(DP)    :: rtemp
      INTEGER     :: wan_func_file, nPlrn_p, indexkn1
      INTEGER     :: iq, inu, ik, ikq, ik_global, ibnd, iplrn, ikpg
      !JLB (dummy variables read from file)
      INTEGER     :: i1, i2
      REAL(DP)    :: r1
      CHARACTER(LEN=5) :: dmmy

      IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)

      IF(ionode) THEN
         wan_func_file = 602
         OPEN(UNIT = wan_func_file, FILE = TRIM(filename))

         IF (PRESENT(scell) .AND. scell) THEN
            READ(wan_func_file, '(a, 3I10)') dmmy, nktotf_p, nbndsub_p, nPlrn_p
            ! nkf1_p, nkf2_p, nkf3_p should never be called if scell=.true.
            ! Just assigning an arbitrary value
            nkf1_p = 0
            nkf2_p = 0
            nkf3_p = 0
         ELSE
            READ(wan_func_file, '(6I10)') nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, nPlrn_p
            IF(nkf1_p*nkf2_p*nkf3_p /= nktotf_p) THEN
               CALL errore("read_plrn_wf",filename//'Not generated from the uniform grid!', 1)
            END IF
         END IF

         ALLOCATE(eigvec_wan(nbndsub_p * nktotf_p, nPlrn_p))

         eigvec_wan = czero
         DO ik = 1, nktotf_p
            DO ibnd = 1, nbndsub_p
               DO iplrn = 1, nPlrn_p
                  indexkn1 = (ik-1)*nbndsub_p + ibnd
                  !JLB
                  IF(PRESENT(etf_all)) THEN
                     READ(wan_func_file, '(2I5, 3f15.7)') i1, i2, r1, eigvec_wan(indexkn1, iplrn)
                  ELSE
                     READ(wan_func_file, '(2f15.7)') eigvec_wan(indexkn1, iplrn)
                  END IF
                  !JLB
                  !READ(wan_func_file, '(2f15.7)') eigvec_wan(indexkn1, iplrn)
               END DO
            END DO
         END DO
         CLOSE(wan_func_file)
      END IF
      CALL mp_bcast (nkf1_p,  meta_ionode_id, world_comm)
      CALL mp_bcast (nkf2_p,  meta_ionode_id, world_comm)
      CALL mp_bcast (nkf3_p,  meta_ionode_id, world_comm)
      CALL mp_bcast (nktotf_p,meta_ionode_id, world_comm)
      CALL mp_bcast (nPlrn_p, meta_ionode_id, world_comm)
      CALL mp_bcast (nbndsub_p, meta_ionode_id, world_comm)
      IF(.NOT. ALLOCATED(eigvec_wan)) THEN
         ALLOCATE(eigvec_wan(nbndsub_p * nktotf_p, nPlrn_p))
         eigvec_wan = czero
      END IF
      CALL mp_bcast (eigvec_wan, meta_ionode_id, world_comm)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   ! Write Bmat and phonon frequency to filename
   SUBROUTINE write_plrn_bmat(eigvec_wan, filename, etf_all)
      USE constants_epw, ONLY : czero, ryd2ev, ryd2mev
      USE elph2,         ONLY : nkf, nqtotf
      USE epwcom,        ONLY : nstate_plrn, nqf1, nqf2, nqf3
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp_world,      ONLY : world_comm
      USE mp,            ONLY : mp_sum, mp_bcast
      USE modes,         ONLY : nmodes
      !JLB
      USE ions_base,     ONLY : amass, ityp
      USE epwcom,        ONLY : scell_mat_plrn


      IMPLICIT NONE
      COMPLEX(DP), INTENT(IN) :: eigvec_wan(:, :)
      REAL(DP), INTENT(IN), OPTIONAL :: etf_all(:, :)
      CHARACTER(LEN=*), INTENT(IN) :: filename

      REAL(DP)    :: rtemp
      INTEGER     :: wan_func_file, indexkn1, nbnd_out
      INTEGER     :: iq, inu, ik, ikq, ik_global, ibnd, iplrn, ikpg
      !!JLB
      !INTEGER     :: ina, ialpha
      !REAL(dp)    :: cm(3) ! change of center of mass

      !IF(ionode) THEN

      IF(PRESENT(etf_all)) THEN
         nbnd_out = nmodes
      ELSE
         nbnd_out = nmodes
      END IF
      wan_func_file = 602

      OPEN(UNIT = wan_func_file, FILE = TRIM(filename))
      IF (scell_mat_plrn) THEN
         WRITE(wan_func_file, '(a, 2I10)') 'Scell', nqtotf, nmodes
      ELSE
         WRITE(wan_func_file, '(5I10)') nqf1, nqf2, nqf3, nqtotf, nmodes
      END IF

      !cm=0.d0
      !ialpha=0
      DO ik = 1, nqtotf
         DO ibnd = 1, nmodes ! p
            IF(PRESENT(etf_all)) THEN ! \kappa, \alpha
               !WRITE(wan_func_file, '(2I5, 4f15.7)') ik, ibnd, etf_all(ibnd, ik)*ryd2mev, &
               !   eigvec_wan(ik, ibnd), ABS(eigvec_wan(ik, ibnd))
               !!JLB: Changed format for improved accuracy
               WRITE(wan_func_file, '(2I5, 4ES18.10)') ik, ibnd, etf_all(ibnd, ik)*ryd2mev, &
                  eigvec_wan(ik, ibnd), ABS(eigvec_wan(ik, ibnd))
            ELSE
               !!JLB
               !ina = (ibnd - 1) / 3 + 1
               !ialpha = ialpha+1
               !cm(ialpha) = cm(ialpha) + amass(ityp(ina))*REAL(eigvec_wan(ik, ibnd),dp)
               !IF (ialpha==3) ialpha=0
               !!JLB
               !JLB: Changed format for improved accuracy
               WRITE(wan_func_file, '(2ES18.10)') eigvec_wan(ik, ibnd)
            END IF
         END DO
      END DO
      CLOSE(wan_func_file)
      !END IF
   END SUBROUTINE
   !-----------------------------------------------------------------------
   ! Read dtau from filename
   SUBROUTINE read_plrn_dtau(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nmodes_p, filename, scell, etf_all)
      USE constants_epw, ONLY : czero
      USE elph2,         ONLY : nkf, nqtotf, nktotf
      USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp_world,      ONLY : world_comm
      USE mp,            ONLY : mp_sum, mp_bcast
      USE modes,         ONLY : nmodes

      IMPLICIT NONE
      COMPLEX(DP), ALLOCATABLE, INTENT(INOUT) :: eigvec_wan(:, :)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      LOGICAL, INTENT(IN), OPTIONAL :: scell
      REAL(DP), INTENT(IN), OPTIONAL :: etf_all(:, :) !JLB
      INTEGER, INTENT(OUT) :: nkf1_p, nkf2_p, nkf3_p, nktotf_p, nmodes_p

      REAL(DP)    :: rtemp
      INTEGER     :: wan_func_file, nPlrn_p
      INTEGER     :: iq, inu, ik, ikq, ik_global, ibnd, iplrn, ikpg, iatm, icount
      !JLB (dummy variables read from file)
      INTEGER     :: i1, i2
      REAL(DP)    :: r1
      CHARACTER(LEN=5) :: dmmy

      IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)

      IF(ionode) THEN
         wan_func_file = 602
         OPEN(UNIT = wan_func_file, FILE = TRIM(filename))

         IF (PRESENT(scell) .AND. scell) THEN
            READ(wan_func_file, '(a, 2I10)') dmmy, nktotf_p, nmodes_p
            ! nkf1_p, nkf2_p, nkf3_p should never be called if scell=.true.
            ! Just assigning an arbitrary value
            nkf1_p = 0
            nkf2_p = 0
            nkf3_p = 0
         ELSE
            READ(wan_func_file, '(5I10)') nkf1_p, nkf2_p, nkf3_p, nktotf_p, nmodes_p
            IF(nkf1_p*nkf2_p*nkf3_p /= nktotf_p) THEN
               CALL errore('read_plrn_dtau', filename//'Not generated from the uniform grid!', 1)
            END IF
         END IF

         IF(nmodes /= nmodes_p) THEN
            CALL errore('read_plrn_dtau', "Number of phonon modes are different with last run", 1)
         END IF

         ALLOCATE(eigvec_wan(nktotf_p, nmodes_p))

         eigvec_wan = czero
         DO icount = 1, nktotf_p
            DO iatm = 1, nmodes_p
               IF(PRESENT(etf_all)) THEN
                  !READ(wan_func_file, '(2I5, 3f15.7)') i1, i2, r1, eigvec_wan(icount, iatm)
                  !JLB: Changed format for improved accuracy
                  READ(wan_func_file, '(2I5, 3ES18.10)') i1, i2, r1, eigvec_wan(icount, iatm)
               ELSE
                  !READ(wan_func_file, '(2f15.7)') eigvec_wan(icount, iatm)
                  !JLB: Changed format for improved accuracy
                  READ(wan_func_file, '(2ES18.10)') eigvec_wan(icount, iatm)
               END IF
               !READ(wan_func_file, '(2f15.7)') eigvec_wan(icount, iatm)
            END DO
         END DO
         CLOSE(wan_func_file)
      END IF
      CALL mp_bcast (nkf1_p,  meta_ionode_id, world_comm)
      CALL mp_bcast (nkf2_p,  meta_ionode_id, world_comm)
      CALL mp_bcast (nkf3_p,  meta_ionode_id, world_comm)
      CALL mp_bcast (nktotf_p,meta_ionode_id, world_comm)
      CALL mp_bcast (nmodes_p, meta_ionode_id, world_comm)
      IF(.NOT. ALLOCATED(eigvec_wan)) THEN
         ALLOCATE(eigvec_wan(nktotf_p, nmodes_p))
         eigvec_wan = czero
      END IF
      CALL mp_bcast (eigvec_wan, meta_ionode_id, world_comm)
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   !! Fourier transform from eigVecIn to eigVecOut
   !! ttype is 'Bloch2Wan' or 'Wan2Bloch'
   !! Parallel version, each pool calculates its own k point set (nkf),
   !! then the mp_sum is used to sum over different pools.
   !! require the correct initialization of Rp_array
   SUBROUTINE plrn_eigvec_tran(ttype, t_rev, eigVecIn, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
         nrr_k, ndegen_k, irvec_r, dims, eigVecOut, ip_center)
      USE constants_epw, ONLY : czero, twopi, ci, cone, two
      USE elph2,         ONLY : nkf, xkf, etf, chw, nktotf
      USE epwcom,        ONLY : nstate_plrn, nbndsub, time_rev_A_plrn
      USE epwcom,        ONLY : nkf1, nkf2, nkf3
      USE wan2bloch,     ONLY : hamwan2bloch !!=> hamwan2bloch_old
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum
      USE mp_world,      ONLY : mpime
      USE io_global,     ONLY : stdout

      IMPLICIT NONE

      COMPLEX(DP), INTENT(IN) :: eigvecIn(:, :)
      LOGICAL,     INTENT(IN) :: t_rev
      INTEGER,     INTENT(IN) :: nkf1_p, nkf2_p, nkf3_p, nbndsub_p
      INTEGER,     INTENT(IN) :: nrr_k, dims, ndegen_k(:,:,:) ! ! Added for polaron calculations by Chao Lian.
      REAL(DP),    INTENT(IN) :: irvec_r(3, nrr_k)
      CHARACTER(LEN=9), INTENT(IN) :: ttype
      INTEGER, INTENT(IN), OPTIONAL :: ip_center(1:3)
      !JLB
      !LOGICAL, OPTIONAL, INTENT(IN) :: readpol

      COMPLEX(DP), INTENT(OUT) :: eigVecOut(:, :)
      COMPLEX(KIND=dp) :: expTable(3)
      INTEGER :: idir

      REAL(DP)    :: rtemp, xxk(3), shift(3), etf_tmp(nbndsub)
      REAL(DP) :: phi, maxreal !JLB
      COMPLEX(DP) :: ctemp, cufkk(nbndsub, nbndsub), cfac(nrr_k, dims, dims), cufkk_k(nbndsub, nbndsub, nktotf), phase !JLB
      COMPLEX(DP), ALLOCATABLE :: cufkkg ( :, :, :)
      INTEGER     :: itype, iq, inu, ik, ikk, ikq, ik_global, iplrn, ikpg, icount
      INTEGER     :: ibnd, jbnd, ix, iy, iz, indexkn1, indexkn2, i_vec(3), center_shift(1:3), nkf_p(3)
      LOGICAL     :: is_mirror

      nkf_p(1:3) = (/nkf1_p, nkf2_p, nkf3_p/)
      IF(nbndsub_p /= nbndsub) CALL errore('plrnwfwan2bloch','Different bands included in last calculation!', 1)
      IF(ttype == 'Bloch2Wan') THEN
         itype =  1
      ELSE IF (ttype == 'Wan2Bloch') THEN
         itype = -1
      ELSE
         CALL errore('plrn_eigvec_tran', 'Illegal translate form; should be Bloch2Wan or Wan2Bloch!', 1)
      END IF

      IF(PRESENT(ip_center)) THEN
         center_shift(1:3) = nkf_p/2 - ip_center
      ELSE
         center_shift(1:3) = 0
      END IF
      !! itype =  1 : Bloch2Wan: A_{mp} =  \frac{1}{N_p} \sum_{nk}A_{nk} \exp\left(ik\cdot R_p\right)U^\dagger_{mnk}
      !! itype = -1 : Wan2Bloch: A_{nk} = \sum_{mp}A_{mp}\exp(-ik\cdot R_p) U_{mnk}
      !! ibnd -> m, jbnd -> n
      !! R_p from 1 to nkf1/2/3_p, note that loop in the sequence of ix, iy, and iz,
      !! This sequence need to be consistent every time transpose between eigvec_wann and eigvec
      eigVecOut = czero
      DO ik = 1, nkf
         xxk = xkf(1:3, 2 * ik - 1)
         expTable(1:3) = EXP( twopi * ci * xxk(1:3) )
         ik_global = ikqLocal2Global(ik, nktotf)
         is_mirror = (t_rev .AND. (ik_global > ikpg))

         CALL get_cfac(xxk, nrr_k, ndegen_k, irvec_r, dims, cfac)

         CALL hamwan2bloch ( nbndsub, nrr_k, cufkk(1:nbndsub, 1:nbndsub), &
            etf_tmp, chw, cfac, dims, is_mirror)

         IF(itype == 1) cufkk(1:nbndsub, 1:nbndsub) = CONJG(TRANSPOSE(cufkk(1:nbndsub, 1:nbndsub)))
         DO iplrn = 1, nstate_plrn
            !icount = 0
            !! loop over all Wannier position p
            IF (nkf1_p == 0 .or. nkf2_p == 0 .or. nkf3_p == 0) THEN
               CALL errore('plrn_eigvec_tran','Wrong k grid, use nkf1/2/3 to give k grid!', 1)
            END IF
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
                              eigVecIn(indexkn2, iplrn) * ctemp/nktotf * cufkk(ibnd, select_bands_plrn(jbnd)) !JLB: Conjugate transpose taken above!
                        CASE(-1) ! Wan2Bloch !
                           eigVecOut(indexkn2, iplrn) = eigVecOut(indexkn2, iplrn) + &
                              eigVecIn(indexkn1, iplrn) * CONJG(ctemp) * cufkk(select_bands_plrn(jbnd), ibnd) !JLB
                     END SELECT
                  END DO ! jbnd
               END DO ! ibnd
            END DO
         END DO !iplrn
      END DO ! ik
      ! MPI sum due to the loop ik is within local k set
      CALL mp_sum( eigVecOut, inter_pool_comm )
   END SUBROUTINE
   !
   !-----------------------------------------------------------------------
   !! JLB: Fourier transform for non-diagonal supercells
   SUBROUTINE scell_plrn_eigvec_tran(ttype, t_rev, eigVecIn, nktotf_p, nRp_p, Rp_p, nbndsub_p, &
         nrr_k, ndegen_k, irvec_r, dims, eigVecOut)
      USE constants_epw, ONLY : czero, twopi, ci, cone, two
      USE elph2,         ONLY : nkf, xkf, etf, chw, nktotf
      USE epwcom,        ONLY : nstate_plrn, nbndsub, time_rev_A_plrn
      USE epwcom,        ONLY : scell_mat_plrn
      USE wan2bloch,     ONLY : hamwan2bloch !!=> hamwan2bloch_old
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum
      USE mp_world,      ONLY : mpime
      USE io_global,     ONLY : stdout, ionode

      IMPLICIT NONE

      COMPLEX(DP), INTENT(IN) :: eigvecIn(:, :)
      LOGICAL,     INTENT(IN) :: t_rev
      INTEGER,     INTENT(IN) :: nktotf_p, nRp_p, Rp_p(:,:), nbndsub_p
      INTEGER,     INTENT(IN) :: nrr_k, dims, ndegen_k(:,:,:) ! ! Added for polaron calculations by Chao Lian.
      REAL(DP),    INTENT(IN) :: irvec_r(3, nrr_k)
      CHARACTER(LEN=9), INTENT(IN) :: ttype

      COMPLEX(DP), INTENT(OUT) :: eigVecOut(:, :)
      COMPLEX(KIND=dp) :: expTable(3)
      INTEGER :: idir

      REAL(DP)    :: rtemp, xxk(3), shift(3), etf_tmp(nbndsub)
      REAL(DP) :: phi, maxreal !JLB
      COMPLEX(DP) :: ctemp, cufkk(nbndsub, nbndsub), cfac(nrr_k, dims, dims), cufkk_k(nbndsub, nbndsub, nktotf), phase !JLB
      COMPLEX(DP), ALLOCATABLE :: cufkkg ( :, :, :)
      INTEGER     :: itype, iq, inu, ik, ikk, ikq, ik_global, iplrn, ikpg, iRp
      INTEGER     :: ibnd, jbnd, ix, iy, iz, indexkn1, indexkn2
      LOGICAL     :: is_mirror

      IF(nbndsub_p /= nbndsub) CALL errore('plrnwfwan2bloch','Different bands included in last calculation!',1)
      IF(ttype == 'Bloch2Wan') THEN
         itype =  1
      ELSE IF (ttype == 'Wan2Bloch') THEN
         itype = -1
      ELSE
         CALL errore('plrn_eigvec_tran', 'Illegal translate form; should be Bloch2Wan or Wan2Bloch!', 1)
      END IF

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
                              eigVecIn(indexkn2, iplrn) * ctemp/nktotf * cufkk(ibnd, select_bands_plrn(jbnd))
                        CASE(-1) ! Wan2Bloch !
                           eigVecOut(indexkn2, iplrn) = eigVecOut(indexkn2, iplrn) + &
                              eigVecIn(indexkn1, iplrn) * CONJG(ctemp) * cufkk(select_bands_plrn(jbnd), ibnd)
                     END SELECT
                  END DO ! jbnd
               END DO ! ibnd
            END DO
         END DO !iplrn
      END DO ! ik
      ! MPI sum due to the loop ik is within local k set
      CALL mp_sum( eigVecOut, inter_pool_comm )
      !
   END SUBROUTINE
   !!JLB
   !-----------------------------------------------------------------------
   !! Interpolate Ank and write to Ank.band.plrn,
   !! especially used to visualize phonon contribution to polaron in band-mode
   SUBROUTINE interp_plrn_wf(nrr_k, ndegen_k, irvec_r, dims)
      USE constants_epw, ONLY : zero, ryd2ev, czero
      USE io_global,     ONLY : stdout, ionode
      USE epwcom,        ONLY : type_plrn, nbndsub, nstate_plrn
      USE elph2,         ONLY : nktotf, etf

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: nrr_k, dims, ndegen_k(:,:,:) ! ! Added for polaron calculations by Chao Lian.
      REAL(DP), INTENT (IN) :: irvec_r(3, nrr_k)

      COMPLEX(DP), ALLOCATABLE :: eigvec_wan(:, :)
      REAL(DP),   ALLOCATABLE :: dtau_r(:, :)
      INTEGER :: iRp, Rp_vec(3), i_center(2), ip_center(3)

      INTEGER  :: nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p
      INTEGER  :: ik_bm, band_pos, ierr
      REAL(DP) :: efermi

      IF(ionode) WRITE(stdout, "(5x, a)") "Start of interpolation of electronic band structure."
      IF(.NOT. ALLOCATED(etf_all)) THEN
         CALL errore('interp_plrn_wf','etf_all should be correctly prepared before calling interp_plrn_wf', 1)
      END IF

      !!FIXME: When selecting band in solving polaron, nbndsub_p should be changed when output
      CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbndsub_p, 'Amp.plrn')

      i_center = MAXLOC(ABS(eigvec_wan))

      ip_center = index_Rp(i_center(1)/nbndsub_p + 1, (/nkf1_p, nkf2_p, nkf3_p/))
      WRITE(stdout, '(5x, a, i8, 3i5)') "The largest Amp ", i_center(1), ip_center

      CALL plrn_eigvec_tran('Wan2Bloch', .false., eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nbndsub_p, &
         nrr_k, ndegen_k, irvec_r, dims, eigVec, ip_center)

      CALL write_plrn_wf(eigVec, 'Ank.band.plrn', etf_all)

      IF(ALLOCATED(eigvec_wan)) DEALLOCATE(eigvec_wan)
   END SUBROUTINE
   !
   !!-----------------------------------------------------------------------
   !! Interpolate bmat and write to Bmat.band.plrn,
   !! especially used to visualize phonon contribution to polaron in band-mode
   SUBROUTINE interp_plrn_bq(nrr_q, ndegen_q, irvec_q, rws, nrws)
      USE elph2,         ONLY : xqf, wf, nqtotf
      USE modes,         ONLY : nmodes
      USE constants_epw, ONLY : czero
      USE io_global,     ONLY : stdout, ionode

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: nrr_q, ndegen_q(:,:,:) ! ! Added for polaron calculations by Chao Lian.
      INTEGER, INTENT (IN) :: irvec_q(3, nrr_q)
      INTEGER,  INTENT (IN) :: nrws
      REAL(DP), INTENT (IN) :: rws(:, :)

      INTEGER :: nqf1_p, nqf2_p, nqf3_p, nqtotf_p, nmodes_p, ierr
      INTEGER :: iRp, ina, Rp_vec(3), i_center(2), ip_center(3)

      COMPLEX(DP), ALLOCATABLE :: Bmat(:,:)
      COMPLEX(DP), ALLOCATABLE :: dtau(:, :)
      REAL(DP),    ALLOCATABLE :: dtau_r(:, :)

      CALL read_plrn_dtau(dtau, nqf1_p, nqf2_p, nqf3_p, nqtotf_p, nmodes_p, 'dtau.plrn')

      ALLOCATE(dtau_r(nqtotf_p, nmodes/3), STAT = ierr)
      IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating Bmat', 1)
      dtau_r = czero
      DO iRp = 1, nqtotf_p
         DO ina = 1, nmodes/3 ! ika -> kappa alpha
            dtau_r(iRp, ina) = NORM2(REAL(dtau(iRp, (ina-1)*3+1:ina*3)))
         END DO
      END DO
      i_center = MAXLOC(ABS(dtau_r))
      ip_center = index_Rp(i_center(1), (/nqf1_p, nqf2_p, nqf3_p/))

      ALLOCATE(Bmat(nqtotf, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('interp_plrn_bq', 'Error allocating Bmat', 1)
      Bmat = czero

      CALL plrn_bmat_tran('Dtau2Bmat', .false., dtau, nqf1_p, nqf2_p, nqf3_p, &
         nrr_q, ndegen_q, irvec_q, rws, nrws, Bmat, ip_center)

      IF(ionode) CALL write_plrn_bmat(Bmat, 'Bmat.band.plrn', wf)

      DEALLOCATE(dtau)
      DEALLOCATE(Bmat)
      DEALLOCATE(dtau_r)
   END SUBROUTINE
   !
   !!-----------------------------------------------------------------------
   !! Fourier transform between Bmat and dtau
   !! Dtau2Bmat : B_{q\nu} = -1/N_p\sum_{\kappa\alpha p}C_{q\kappa \nu}\Delta\tau_{\kappa\alpha p}  e_{\kappa\alpha\nu}(q)\exp(iqR_p)
   !! Bmat2Dtau : \Delta \tau_{\kappa\alpha p} = -\sum_{q\nu} 1/(C_{q\kappa \nu}) B^*_{q\nu} e_{\kappa\alpha,\nu}(q) \exp(iqR_p)
   !! C_{q\kappa \nu} = N_p\left(\frac{M_k\omega_{q\nu}}{2\hbar}\right)^{\frac{1}{2}} = N_p(M_k)^{\frac{1}{2}}D_{q\nu}
   !! D_{q \nu} = \left(\frac{\omega_{q\nu}}{2\hbar}\right)^{\frac{1}{2}}
   SUBROUTINE plrn_bmat_tran(ttype, t_rev, mat_in, nqf1_p, nqf2_p, nqf3_p, &
         nrr_q, ndegen_q, irvec_q, rws, nrws, mat_out, ip_center)
      USE elph2,         ONLY : xqf, nqtotf, nkf, wf
      USE modes,         ONLY : nmodes
      USE constants_epw, ONLY : eps8, czero, one, two, twopi, zero, ci, cone
      USE ions_base,     ONLY : amass, ityp
      USE wan2bloch,     ONLY : dynwan2bloch, dynifc2blochf
      USE epwcom,        ONLY : lifc
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum
      USE division,      ONLY : fkbounds
      USE io_global,     ONLY : stdout
      USE epwcom,        ONLY : type_plrn

      IMPLICIT NONE

      CHARACTER(LEN=9), INTENT(IN) :: ttype
      INTEGER, INTENT (IN) :: nqf1_p, nqf2_p, nqf3_p
      LOGICAL, INTENT (IN) :: t_rev
      COMPLEX(DP), INTENT(IN) :: mat_in(:, :)
      INTEGER, INTENT (IN) :: nrr_q, ndegen_q(:,:,:)
      INTEGER, INTENT (IN) :: irvec_q(3, nrr_q)
      INTEGER,  INTENT (IN) :: nrws
      REAL(DP), INTENT (IN) :: rws(:, :)
      COMPLEX(DP), INTENT(OUT) :: mat_out(:, :)
      INTEGER, INTENT(IN), OPTIONAL :: ip_center(1:3)
      !JLB
      !LOGICAL, OPTIONAL, INTENT(IN) :: readpol
      LOGICAL     :: mirror_q
      INTEGER     :: iq, inu, ierr, imu, iatm, idir, itype, ika, ip_start, ip_end, iRp, nqf_p(1:3)
      INTEGER     :: ix, iy, iz, ina, nqtotf_p, iqpg, nqf, nptotf, start_modes, Rp_vec(1:3), center_shift(1:3)
      COMPLEX(DP) :: dtemp, shift(3), expTable(3), uf(nmodes, nmodes)
      REAL(DP)    :: xxq(3), xxq_r(3), ctemp, w2(nmodes)!, wf(:,:)
      !JLB
      INTEGER :: jnu, ndegen(nmodes), imode, jmode

      nptotf = nqf1_p * nqf2_p * nqf3_p
      nqf_p(1:3) = (/nqf1_p, nqf2_p, nqf3_p/)

      IF(nptotf <= 0) CALL errore('plrn_eigvec_tran', 'Use correct .plrn file with nqf1_p \= 0!', 1)
      IF(ttype == 'Bmat2Dtau') THEN
         itype =  1
      ELSE IF (ttype == 'Dtau2Bmat') THEN
         itype = -1
      ELSE
         CALL errore('plrn_eigvec_tran', 'Illegal translation form; should be Bmat2Dtau or Dtau2Bmat!', 1)
      END IF

      uf = czero
      w2 = zero
      wf = zero

      mat_out = czero

      CALL fkbounds(nptotf, ip_start, ip_end)

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
            END IF
         END IF
         expTable(1:3) = EXP( twopi * ci * xxq(1:3) )

         ! Get phonon eigenmode and eigenfrequencies
         IF (.NOT. lifc) THEN
            ! Incompatible bugs found 9/4/2020 originated from the latest EPW changes.
            ! parallel q is not working any more due to mp_sum in rgd_blk
            CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, mirror_q)
         ELSE
            CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, mirror_q)
         ENDIF

         DO inu = 1, nmodes
            IF (w2(inu) > -eps8) THEN
               wf(inu, iq) =  DSQRT(ABS(w2(inu)))
            ELSE
               !wf(inu, iq) = -DSQRT(ABS(w2(inu)))
               wf(inu, iq) = 0.d0
            ENDIF
         END DO

         IF(PRESENT(ip_center)) THEN
            center_shift(1:3) = nqf_p/2 - ip_center
            !write(stdout, *) center_shift, ip_center
         ELSE
            center_shift(1:3) = 0
         END IF
         ! For mirror q, calculate the time-symmetric q' and get uf from q'
         ! e_{\kappa\alpha\nu}(-q)= e^*_{\kappa\alpha\nu}(q)
         !!IF(t_rev .and. iq > iqpg) uf = CONJG(uf) !transpose
         start_modes = 1
         DO inu = start_modes, nmodes ! inu -> nu
            IF (wf(inu, iq) < eps8) CYCLE !JLB - cycle zero and imaginary frequency modes
            DO ika = 1, nmodes ! ika -> kappa alpha
               ina = (ika - 1) / 3 + 1
               ctemp = DSQRT(two/(wf(inu, iq) * amass(ityp(ina))))
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
                     mat_out(iRp, ika) = mat_out(iRp, ika) -  cone/REAL(nptotf, dp) * dtemp * ctemp &
                        * (-type_plrn) * CONJG(mat_in(iq, inu))
                  ELSE IF(itype == -1) THEN
                     !  B_{q\nu} = \frac{1}{N_p} \sum_{\kappa\alpha p} D_{\kappa \alpha\nu, p}(q) C_{q}\nu \Delta\tau_{\kappa\alpha p}
                     mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * CONJG(mat_in(iRp, ika)) !JLB: dtau should be real but just in case
                     !mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * mat_in(iRp, ika)
                  END IF
               END DO
            END DO
         END DO
      END DO
      ! sum all the cell index ip
      CALL mp_sum(mat_out, inter_pool_comm )

   END SUBROUTINE
   !!-----------------------------------------------------------------------
   !! JLB: Fourier transform between Bmat and dtau for non-diagonal supercells
   SUBROUTINE scell_plrn_bmat_tran(ttype, t_rev, mat_in, nqtotf_p, nRp_p, Rp_p, &
         nrr_q, ndegen_q, irvec_q, rws, nrws, mat_out)
      USE elph2,         ONLY : xqf, nqtotf, nkf, wf
      USE modes,         ONLY : nmodes
      USE constants_epw, ONLY : eps8, czero, one, two, twopi, zero, ci, cone
      USE ions_base,     ONLY : amass, ityp
      USE wan2bloch,     ONLY : dynwan2bloch, dynifc2blochf
      USE epwcom,        ONLY : lifc
      USE mp_global,     ONLY : inter_pool_comm
      USE mp,            ONLY : mp_sum
      USE division,      ONLY : fkbounds
      USE io_global,     ONLY : stdout, ionode
      USE epwcom,        ONLY : type_plrn

      IMPLICIT NONE

      CHARACTER(LEN=9), INTENT(IN) :: ttype
      INTEGER, INTENT (IN) :: nqtotf_p, nRp_p, Rp_p(:,:)
      LOGICAL, INTENT (IN) :: t_rev
      COMPLEX(DP), INTENT(IN) :: mat_in(:, :)
      INTEGER, INTENT (IN) :: nrr_q, ndegen_q(:,:,:)
      INTEGER, INTENT (IN) :: irvec_q(3, nrr_q)
      INTEGER,  INTENT (IN) :: nrws
      REAL(DP), INTENT (IN) :: rws(:, :)
      COMPLEX(DP), INTENT(OUT) :: mat_out(:, :)
      LOGICAL     :: mirror_q
      INTEGER     :: iq, inu, ierr, imu, iatm, idir, itype, ika, ip_start, ip_end, iRp
      INTEGER     :: ix, iy, iz, ina, iqpg, nqf, start_modes
      COMPLEX(DP) :: dtemp, uf(nmodes, nmodes)
      REAL(DP)    :: xxq(3), xxq_r(3), ctemp, w2(nmodes)!, wf(:,:)
      INTEGER     :: jnu, ndegen(nmodes), imode, jmode

      IF(ttype == 'Bmat2Dtau') THEN
         itype =  1
      ELSE IF (ttype == 'Dtau2Bmat') THEN
         itype = -1
      ELSE
         CALL errore('plrn_eigvec_tran', 'Illegal translation form; should be Bmat2Dtau or Dtau2Bmat!', 1)
      END IF

      uf = czero
      w2 = zero
      wf = zero

      mat_out = czero

      CALL fkbounds(nRp_p, ip_start, ip_end)

      DO iq = 1, nqtotf_p ! iq -> q
         xxq = xqf(1:3, iq)
         xxq_r = xxq(1:3)
         mirror_q = .false.
         ! if we need to force the time-rev symmetry, we have to ensure that the phase of uf is fixed
         ! i.e. uf = uf*(-q)
         IF (t_rev) THEN
            IF (is_mirror_q (iq)) THEN
               xxq_r = xqf(1:3, kpg_map(iq))
               mirror_q = .true.
            END IF
         END IF

         ! Get phonon eigenmode and eigenfrequencies
         IF (.NOT. lifc) THEN
            ! Incompatible bugs found 9/4/2020 originated from the latest EPW changes.
            ! parallel q is not working any more due to mp_sum in rgd_blk
            CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, mirror_q)
         ELSE
            CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, mirror_q)
         ENDIF

         DO inu = 1, nmodes
            IF (w2(inu) > -eps8) THEN
               wf(inu, iq) =  DSQRT(ABS(w2(inu)))
            ELSE
               !wf(inu, iq) = -DSQRT(ABS(w2(inu)))
               wf(inu, iq) = 0.d0
            ENDIF
         END DO

         ! For mirror q, calculate the time-symmetric q' and get uf from q'
         ! e_{\kappa\alpha\nu}(-q)= e^*_{\kappa\alpha\nu}(q)
         !!IF(t_rev .and. iq > iqpg) uf = CONJG(uf) !transpose
         DO inu = 1, nmodes ! inu -> nu
            IF (wf(inu, iq) < eps8) CYCLE !JLB - cycle zero and imaginary frequency modes
            DO ika = 1, nmodes ! ika -> kappa alpha
               ina = (ika - 1) / 3 + 1
               ctemp = DSQRT(two/(wf(inu, iq) * amass(ityp(ina))))
               !
               !DO iRp = 1, nRp
               DO iRp = ip_start, ip_end
                  ! D_{\kappa\alpha\nu,p}(q) = e_{\kappa\alpha,\nu}(q) \exp(iq\cdot R_p)
                  dtemp = uf(ika, inu) * EXP( twopi * ci * DOT_PRODUCT(xxq(1:3), Rp_p(1:3, iRp)) )
                  IF (itype == 1) THEN ! Bqv -> dtau
                     ! \Delta \tau_{\kappa\alpha p} = -\frac{1}{N_p} \sum_{q\nu} C_{\kappa\nu q} D_{\kappa\alpha\nu q}  B^*_{q\nu}
                     ! Dtau(iRp, ika) = Dtau(iRp, ika) + conjg(B(iq, inu)) * ctemp * dtemp
                     mat_out(iRp, ika) = mat_out(iRp, ika) -  cone/REAL(nRp, dp) * dtemp * ctemp &
                        * (-type_plrn) * CONJG(mat_in(iq, inu))
                  ELSE IF(itype == -1) THEN
                     !  B_{q\nu} = \frac{1}{N_p} \sum_{\kappa\alpha p} D_{\kappa \alpha\nu, p}(q) C_{q}\nu \Delta\tau_{\kappa\alpha p}
                     mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * CONJG(mat_in(iRp, ika)) !JLB: dtau should be real but just in case
                     !mat_out(iq, inu) = mat_out(iq, inu) - (-type_plrn) * dtemp/ctemp * mat_in(iRp, ika)
                  END IF
               END DO
            END DO
         END DO
      END DO
      ! sum all the cell index ip
      CALL mp_sum(mat_out, inter_pool_comm )
      !
   END SUBROUTINE
   !-----------------------------------------------------------------------
   SUBROUTINE calc_den_of_state(eigVec, Bmat)
      USE epwcom,        ONLY : nDOS_plrn, edos_max_plrn, edos_min_plrn, edos_sigma_plrn
      USE epwcom,        ONLY : pdos_max_plrn, pdos_min_plrn, pdos_sigma_plrn
      USE epwcom,        ONLY : nstate_plrn, nkf1, nkf2, nkf3
      USE elph2,         ONLY : nkf, nqtotf, nktotf, xkf
      USE elph2,         ONLY : ibndmin, ibndmax, wf
      USE modes,         ONLY : nmodes

      USE constants_epw, ONLY : ryd2mev, czero, one, ryd2ev, two, zero, cone
      USE constants_epw, ONLY : pi, ci, twopi, eps6, eps8, eps5

      IMPLICIT NONE

      INTEGER     :: idos
      INTEGER     :: iq, inu, ik, ikk, ikq, ik_global, iplrn, ikpg, icount
      INTEGER     :: ibnd, jbnd, ix, iy, iz
      INTEGER     :: dos_file, indexkn1

      COMPLEX(DP), INTENT(IN) :: eigVec(:, :), Bmat(:, :)

      REAL(DP), ALLOCATABLE :: rmat_tmp(:, :), rvec_tmp(:)
      REAL(DP), ALLOCATABLE :: edos(:), pdos(:), edos_all(:), pdos_all(:), e_grid(:), p_grid(:)
      REAL(DP) :: temp

      !Calculating DOS
      ALLOCATE(rvec_tmp(nDOS_plrn))
      ALLOCATE(e_grid(nDOS_plrn))
      ALLOCATE(edos(nDOS_plrn))
      ALLOCATE(edos_all(nDOS_plrn))
      temp = MAXVAL(etf_all) * ryd2ev + one
      IF(edos_max_plrn < temp) edos_max_plrn = temp
      temp = MINVAL(etf_all) * ryd2ev - one
      IF(edos_min_plrn > temp) edos_min_plrn = temp

      e_grid = zero
      DO idos = 1, nDOS_plrn
         e_grid(idos) = edos_min_plrn + idos*(edos_max_plrn - edos_min_plrn)/(nDOS_plrn)
      END DO

      edos = zero
      edos_all = zero
      DO ik = 1, nktotf
         DO ibnd = 1, nbnd_plrn
            ! TODO : iplrn
            DO iplrn = 1, 1
               CALL cal_f_delta(e_grid - (etf_all(select_bands_plrn(ibnd), ik) * ryd2ev), &
                  edos_sigma_plrn, rvec_tmp)
               indexkn1 = (ik-1)*nbnd_plrn + ibnd
               edos = edos + (ABS(eigVec(indexkn1, iplrn))**2)*rvec_tmp
               edos_all = edos_all + rvec_tmp
            END DO
         END DO
      END DO

      ALLOCATE(p_grid(nDOS_plrn))
      ALLOCATE(pdos(nDOS_plrn))
      ALLOCATE(pdos_all(nDOS_plrn))

      temp = MAXVAL(wf) * ryd2mev + 10.0_dp
      IF(pdos_max_plrn < temp) pdos_max_plrn = temp
      temp = MINVAL(wf) * ryd2mev - 10.0_dp
      IF(pdos_min_plrn > temp) pdos_min_plrn = temp

      p_grid = zero
      DO idos = 1, nDOS_plrn
         p_grid(idos) = pdos_min_plrn + idos*(pdos_max_plrn - pdos_min_plrn)/(nDOS_plrn)
      END DO

      pdos = zero
      pdos_all = zero
      DO iq = 1, nqtotf
         DO inu = 1, nmodes
            CALL cal_f_delta(p_grid - wf(inu, iq) * ryd2mev, pdos_sigma_plrn, rvec_tmp)
            pdos = pdos + (ABS(Bmat(iq, inu))**2)*rvec_tmp
            pdos_all = pdos_all + rvec_tmp
         END DO
      END DO

      dos_file = 601
      OPEN(UNIT = dos_file, FILE = 'dos.plrn')
      WRITE(dos_file, '(/2x, a/)') '#energy(ev)  A^2   edos  energy(mev)  B^2  pdos'
      DO idos = 1, nDOS_plrn
         WRITE(dos_file, '(6f15.7)') e_grid(idos), edos(idos), &
            edos_all(idos), p_grid(idos), pdos(idos), pdos_all(idos)
      END DO
      CLOSE(dos_file)

      DEALLOCATE(e_grid, edos, edos_all)
      DEALLOCATE(p_grid, pdos, pdos_all)
      DEALLOCATE(rvec_tmp)
   END SUBROUTINE
   !
   SUBROUTINE write_real_space_wavefunction()
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

      IMPLICIT NONE

      INTEGER :: nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbnd_plrn_p
      INTEGER :: nqf_p(3), nqtotf_p, nmodes_p
      INTEGER :: ibnd, jbnd, iline, nAtoms, idir, file_unit, igrid, itemp, indexkn1
      INTEGER :: nxx, nyy, nzz, ipx, ipy, ipz, ip_min, ip_max, ig_vec(1:3)
      INTEGER :: iscx, iscy, iscz, ie, i_species, ivec, iRp
      INTEGER :: n_grid(3), grid_start(3), grid_end(3), n_grid_super(3), r_grid_vec(3)
      INTEGER :: rpc(1:3), ipc(1:3), shift_start(1:3), species(50), Rp_vec(1:3), shift(1:3), ishift
      COMPLEX(dp), ALLOCATABLE :: eigvec_wan(:, :), dtau(:, :), cvec(:)
      COMPLEX(dp) :: Amp, ctemp(1:3), b_vec(1:3)
      CHARACTER(LEN = 60) :: plrn_file
      REAL(dp), ALLOCATABLE :: wann_func(:, :, :, :), rvec(:)
      REAL(dp) :: orig(3), rtempvec(4), cell(3, 3), rtemp(5), tempDen(5,5,5)
      REAL(dp) :: r_cry(3), r_cart(3)

      file_unit = 60512
      ! read Amp.plrn, save eigvec_wan for the latter use
      CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbnd_plrn_p, 'Amp.plrn')
      ! read dtau.plrn, get the displacement.
      CALL read_plrn_dtau(dtau, nqf_p(1), nqf_p(2), nqf_p(3), nqtotf_p, nmodes_p, 'dtau.plrn')
      ! read cube files for the real-space Wannier function Wm(r)
      CALL read_wannier_cube(select_bands_plrn, wann_func, species, &
         n_grid, grid_start, grid_end)

      cell(1:3, 1) = at(1:3, 1) * nqf_p(1) * alat
      cell(1:3, 2) = at(1:3, 2) * nqf_p(2) * alat
      cell(1:3, 3) = at(1:3, 3) * nqf_p(3) * alat

      plrn_file = 'psir_plrn.xsf'
      ! Write the file head including information of structures,
      ! using the same format of
      CALL write_plrn_dtau_xsf(dtau, nqf_p(1), nqf_p(2), nqf_p(3), plrn_file, species)

      orig(1:3) = zero
      n_grid_super(1:3) = nqf_p(1:3) * n_grid(1:3)

      IF(ionode) THEN
         OPEN(UNIT = file_unit, FILE = TRIM(plrn_file), POSITION='APPEND')
         WRITE (file_unit, '(/)')
         WRITE (file_unit, '("BEGIN_BLOCK_DATAGRID_3D",/,"3D_field",/, "BEGIN_DATAGRID_3D_UNKNOWN")')
         WRITE (file_unit, '(3i6)')  n_grid_super / step_wf_grid_plrn
         WRITE (file_unit, '(3f12.6)') zero, zero, zero
         WRITE (file_unit, '(3f12.7)') cell(1:3, 1)*bohr2ang
         WRITE (file_unit, '(3f12.7)') cell(1:3, 2)*bohr2ang
         WRITE (file_unit, '(3f12.7)') cell(1:3, 3)*bohr2ang
      END IF

      b_vec(1:3) = twopi * ci / REAL(n_grid_super(1:3))


      !allocate(cvec(grid_start_min(1): n_grid_super(1) + grid_start_min(1) - 1))
      ALLOCATE(cvec(1:n_grid_super(1)))

      CALL fkbounds(nqtotf_p, ip_min, ip_max)

      !      do nzz = grid_start_min(3), n_grid_super(3) + grid_start_min(3) - 1
      !         do nyy = grid_start_min(2), n_grid_super(2) + grid_start_min(2) - 1
      !            cvec = czero
      !            do nxx = grid_start_min(1), n_grid_super(1) + grid_start_min(1) - 1
      ctemp = czero
      rtemp = czero
      DO nzz = 1, n_grid_super(3), step_wf_grid_plrn
         !IF (MOD(nzz - 1, step_wf_grid_plrn)/=0) CYCLE
         DO nyy = 1, n_grid_super(2), step_wf_grid_plrn
            !IF (MOD(nyy - 1, step_wf_grid_plrn)/=0) CYCLE
            cvec = czero
            DO nxx = 1, n_grid_super(1), step_wf_grid_plrn
               !IF (MOD(nxx - 1, step_wf_grid_plrn)/=0) CYCLE
               !icount = 0
               DO iRp = ip_min, ip_max !-nqf_p(3)/2, (nqf_p(3)+1)/2
                  Rp_vec(1:3) = index_Rp(iRp, nqf_p)
                  rpc(1:3) = Rp_vec(1:3) * n_grid(1:3) !- (n_grid_super(1:3))/2
                  ! To make sure that all the nonzero points are included,
                  ! we need to try from -1 to 1 neighbor supercells
                  DO ishift = 1, 27
                     shift(1:3) = index_shift(ishift)
                     ig_vec(1:3) = (/nxx, nyy, nzz/) - rpc(1:3) + shift(1:3) * n_grid_super(1:3)
                     IF(ALL(ig_vec(1:3) <= grid_end(1:3)) .and. &
                        ALL(ig_vec(1:3) >= grid_start(1:3))) THEN
                        !print *,"igx, igy, igz, nxx, nyy, nzz, ipx, ipy, ipz:",igx, igy, igz, nxx, nyy, nzz, ipx, ipy, ipz
                        DO ibnd = 1, nbndsub !TODO change to nbndsub
                           !cvec(nxx) = cvec(nxx) + wann_func(igx, igy, igz, ibnd)
                           indexkn1 = (iRp - 1) * nbndsub + ibnd
                           !TODO eigvec_wan(indexkn1, 1) should be eigvec_wan(indexkn1, iplrn)
                           cvec(nxx) = cvec(nxx)  +  eigvec_wan(indexkn1, 1) * wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd)
                           !cvec(nxx) = cvec(nxx)  +  wann_func(igx, igy, igz, ibnd) !*  eigvec_wan(indexkn1, 1)
                        END DO !ibnd
                     END IF
                        !END DO ! iscz
                     !END DO ! iscy
                  END DO ! iscx
                     !END DO ! ipz
                  !ND DO ! ipy
               END DO ! ipx
            END DO ! nxx
            CALL mp_sum(cvec, inter_pool_comm)
            IF(ionode) THEN
               !WRITE (file_unit, '(5e13.5)', ADVANCE='yes') REAL(cvec(::step_wf_grid_plrn))
               !JLB: Changed to |\Psi(r)|^{2}, I think it's physically more meaningful
               WRITE (file_unit, '(5e13.5)', ADVANCE='yes') ABS(cvec(::step_wf_grid_plrn))**2
            END IF
            ! Calculate the center of polaron
            ! TODO: not parallel, all the processors are doing the same calculations
            DO nxx = 1, n_grid_super(1)
               r_grid_vec(1:3) = (/nxx - 1, nyy - 1, nzz - 1/)
               ctemp(1:3) = ctemp(1:3) + EXP(b_vec(1:3) * r_grid_vec(1:3)) * (cvec(nxx)**2)
            END DO
            ! End calculating the center of polaron
         END DO
      END DO

      r_cry(1:3) = IMAG(LOG(ctemp(1:3)))/twopi
      ! make crystal coordinates with 0 to 1
      r_cry(1:3) = r_cry - FLOOR(r_cry)
      r_cart(1:3) = zero
      DO idir = 1, 3
         r_cart(1:3) = r_cart(1:3) + r_cry(idir) * cell(1:3, idir)
      END DO


      IF(ionode) THEN
         !WRITE(stdout, "(5x, 'The norm of polaron wavefunction:', 3e11.4)") rtemp
         WRITE(stdout, "(5x, 'The position of polaron:')")
         WRITE(stdout, "(5x, 3f9.4, ' in crystal coordinates')") r_cry(1:3)
         WRITE(stdout, "(5x, 3f9.4, ' in Cartesian coordinates (Angstrom)')") r_cart(1:3)
         WRITE (file_unit, '("END_DATAGRID_3D",/, "END_BLOCK_DATAGRID_3D")')
         CLOSE(file_unit)
         !JLB
         WRITE(stdout, "(5x, '|\Psi(r)|^2 written to file.')")
      END IF
   !DO ipx = 1, nqf_p(1) !-nqf_p(3)/2, (nqf_p(3)+1)/2
      !DO ipy = 1, nqf_p(2)!-nqf_p(3)/2, (nqf_p(3)+1)/2
         !DO ipz = 1, nqf_p(3)!-nqf_p(3)/2, (nqf_p(3)+1)/2
            !icount = icount + 1
               !rpc(1:3) = ((/ipx - 1, ipy - 1, ipz - 1/) - nqf_p(1:3)/2) * n_grid(1:3) !- (n_grid_super(1:3))/2
               !rpc(1:3) = ((/ipx, ipy, ipz/) - nqf_p(1:3)/2) * n_grid(1:3) !- (n_grid_super(1:3))/2
                                 !DO iscx = -1, 1
         !DO iscy = -1, 1
            !DO iscz = -1, 1
!         igx = nxx - rpc(1) + iscx * n_grid_super(1)
!         igy = nyy - rpc(2) + iscy * n_grid_super(2)
!         igz = nzz - rpc(3) + iscz * n_grid_super(3)
      DEALLOCATE(wann_func, cvec, eigvec_wan)
   END SUBROUTINE
   !
   !
   !!JLB: write psir in transformed supercell
   !!     xsf format no longer compatible, 
   !!     as .cube files are written in primitive coords.
   SUBROUTINE scell_write_real_space_wavefunction()

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

      IMPLICIT NONE

      CHARACTER(LEN = 60) :: plrn_file
      INTEGER :: nktotf_p, nkf1_p, nkf2_p, nkf3_p, nbnd_plrn_p, species(50)
      INTEGER :: ibnd, indexkn1, ig_vec(1:3)
      INTEGER :: ir1, ir2, ir3, iRp1, iRp2, n_grid_total
      INTEGER :: ip_min, ip_max
      INTEGER :: n_grid(3), grid_start(3), grid_end(3), n_grid_super(3), r_in_crys_p_sup(3)
      INTEGER :: ishift, shift(3)
      REAL(DP), ALLOCATABLE :: wann_func(:, :, :, :), rvec(:)
      REAL(DP) :: r_in_crys_p(3), r_in_crys_s(3), r_in_cart(3), Rp_in_cart(3)
      REAL(DP) :: p2s(3,3), s2p(3,3)
      COMPLEX(DP) :: cvec
      COMPLEX(DP), ALLOCATABLE :: eigvec_wan(:, :)

      ! Broadcast supercell lattice vectors
      CALL mp_bcast(as, meta_ionode_id, world_comm)

      ! read Amp.plrn, save eigvec_wan
      CALL read_plrn_wf(eigvec_wan, nkf1_p, nkf2_p, nkf3_p, nktotf_p, nbnd_plrn_p, 'Amp.plrn', .true.)
      ! read cube files for the real-space Wannier function Wm(r)
      CALL read_wannier_cube(select_bands_plrn, wann_func, species, &
         n_grid, grid_start, grid_end)
      ! Read list of Rp-s within supercell
      CALL read_Rp_in_S()
      WRITE(stdout, '(a, i12)') "     Number of unit cells within supercell:", nRp

      ! Open file
      !plrn_file = 'psir_plrn.scell.dat'
      plrn_file = 'psir_plrn.scell.csv'
      IF (ionode) THEN
        OPEN (UNIT=iunpsirscell, FILE=TRIM(plrn_file), FORM='formatted', STATUS='unknown')
        !WRITE (iunpsirscell, '(a)') "# x , y , z (Angstrom), |\psi(r)|^2"
        WRITE (iunpsirscell, '(a)') "x , y , z, |\psi(r)|^2"
      END IF

      ! Total number of grid points
      n_grid_total = nRp*n_grid(1)*n_grid(2)*n_grid(3)
      WRITE(stdout, '(a, i12)') "     Total grid points:", n_grid_total
      WRITE(stdout, '(a, i12)') "     Step:", step_wf_grid_plrn

      ! Parallelize iRp
      CALL fkbounds(nRp, ip_min, ip_max)

      ! Matrix to transform from primitive to supercell crystal coordinates
      p2s = matinv3(TRANSPOSE(as))
      p2s = MATMUL(p2s,at)
      ! Supercell to primitive coordinates
      s2p = matinv3(p2s)
      !s2p = matinv3(REAL(scell_mat, DP))


      ! Loop over all the grid points
      DO iRp1 = 1, nRp
         !
         DO ir1 = 1, n_grid(1), step_wf_grid_plrn
            DO ir2 = 1, n_grid(2), step_wf_grid_plrn
               DO ir3 = 1, n_grid(3), step_wf_grid_plrn
                  !
                  r_in_crys_p(1:3) = (/ REAL(ir1-1, DP)/n_grid(1), REAL(ir2-1, DP)/n_grid(2), REAL(ir3-1, DP)/n_grid(3) /) &
                                     + REAL(Rp(1:3, iRp1),DP)
                  !
                  ! Wannier functions stored in (1:ngrid*iRp) list
                  !r_in_crys_p_sup(1:3) = n_grid(1:3) * &
                  !                       ((/ REAL(ir1, DP)/n_grid(1), REAL(ir2, DP)/n_grid(2), REAL(ir3, DP)/n_grid(3) /) + Rp(1:3, iRp1))
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
                       !ig_vec(1:3) = r_in_crys_p_sup(1:3) - Rp(1:3, iRp2) * n_grid(1:3)
                       !
                       IF(ALL(ig_vec(1:3) <= grid_end(1:3)) .and. &
                          ALL(ig_vec(1:3) >= grid_start(1:3))) THEN
                          !
                          DO ibnd = 1, nbndsub ! sum over m
                             !
                             indexkn1 = (iRp2 - 1) * nbndsub + ibnd
                             cvec = cvec + eigvec_wan(indexkn1, 1) * wann_func(ig_vec(1), ig_vec(2), ig_vec(3), ibnd)
                             !
                          END DO !ibnd
                          !
                       END IF
                       !
                     END DO ! ishift
                     !
                  END DO ! iRp2
                  CALL mp_sum(cvec, inter_pool_comm)
                  !
                  ! Write |\psi(r)|^2 data point to file
                  IF (ionode) THEN
                     !WRITE (iunpsirscell, '(3f12.6, E13.5)') r_in_cart(1:3), ABS(cvec)**2
                     WRITE (iunpsirscell, '(f12.6,", ", f12.6,", ", f12.6,", ", E13.5)') r_in_cart(1:3), ABS(cvec)**2
                  END IF
                  !
               END DO ! ir3
            END DO ! ir2
         END DO ! ir1
      END DO ! iRp1
      !
      IF (ionode) THEN
        CLOSE (iunpsirscell)
        WRITE(stdout, "(5x, '|\Psi(r)|^2 written to file.')")
      END IF
      DEALLOCATE(wann_func, eigvec_wan)
      !
   END SUBROUTINE
   !
   ! Read the nth Wannier function from prefix_0000n.cube file
   SUBROUTINE read_wannier_cube(select_bands, wann_func, species, n_grid, grid_start_min, grid_end_max)
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

      IMPLICIT NONE

      INTEGER,  INTENT(IN)  :: select_bands(:)
      REAL(dp), ALLOCATABLE, INTENT(OUT) :: wann_func(:, :, :, :)
      INTEGER,  INTENT(OUT) :: species(50)
      INTEGER,  INTENT(OUT) :: n_grid(3), grid_start_min(3), grid_end_max(3)

      INTEGER :: ibnd_index, ibnd, ie, idir, i_species, nbnd, iline, nAtoms
      INTEGER :: nxx, nyy, nzz, n_len_z
      INTEGER :: grid_start(3), grid_end(3)
      REAL(dp) :: rtempvec(4), norm
      CHARACTER(LEN = 60) :: wancube, temp_str

      nbnd = SIZE(select_bands)
      !print *, "nbnd", nbnd

      ! find the max and min of real space grid of Wannier functions of all Wannier orbitals
      IF(ionode) THEN
         grid_start_min(:) = 100000
         grid_end_max(:) = -100000
         DO ibnd = 1, nbndsub
            WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
            OPEN(UNIT = iun_plot, FILE=TRIM(wancube), FORM='formatted', STATUS='unknown')
            READ(iun_plot, *) temp_str !, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str, temp_str
            !print *, temp_str
            READ(iun_plot, *) n_grid, grid_start, grid_end
            DO idir = 1, 3
               IF(grid_start_min(idir) >= grid_start(idir)) grid_start_min(idir) = grid_start(idir)
               IF(grid_end_max(idir) <= grid_end(idir))   grid_end_max(idir) = grid_end(idir)
            END DO
            CLOSE(iun_plot)
         END DO
         !print *, "     n_grid, grid_start, grid_end", n_grid, grid_start_min, grid_end_max
      END IF

      CALL mp_bcast(n_grid,         meta_ionode_id, world_comm)
      CALL mp_bcast(grid_start_min, meta_ionode_id, world_comm)
      CALL mp_bcast(grid_end_max,   meta_ionode_id, world_comm)
      !CALL mp_bcast(shift_start,    meta_ionode_id, world_comm)

      ! Read the xth Wannier functions from prefix_0000x.cube in ionode
      ! and broadcast to all nodes
      ALLOCATE(wann_func(grid_start_min(1):grid_end_max(1), &
         grid_start_min(2):grid_end_max(2), &
         !grid_start_min(3):grid_end_max(3), nbnd))
         grid_start_min(3):grid_end_max(3), nbndsub))
      wann_func = zero
      species = 0
      IF(ionode) THEN
         DO ibnd = 1, nbndsub
            WRITE(wancube, "(a, '_', i5.5, '.cube')") TRIM(prefix), ibnd
            OPEN(UNIT = iun_plot, FILE=TRIM(wancube), FORM='formatted', STATUS='unknown')
            READ(iun_plot, *) temp_str
            READ(iun_plot, *) n_grid, grid_start, grid_end
            READ(iun_plot, *) nAtoms, rtempvec(1:3)

            DO iline = 1, 3
               READ(iun_plot, '(8A)') temp_str
            END DO
            !read(iun_plot, '(i, 8A)') i_species, temp_str
            !species(1) = i_species
            ie = 1
            DO iline = 1, nAtoms
               READ(iun_plot, '(i4, 4f13.5)') i_species, rtempvec
               !print *, i_species
               IF (iline == 1 ) THEN
                  species(ie) = i_species
                  ie = ie + 1
               ELSE IF (species(ie - 1) /= i_species) THEN
                  species(ie) = i_species
                  ie = ie + 1
               END IF
            END DO
            n_len_z = grid_end(3) - grid_start(3) + 1

            DO nxx = grid_start(1), grid_end(1)
              DO nyy = grid_start(2), grid_end(2)
                DO nzz = grid_start(3), grid_end(3), 6
                  IF (grid_end(3) - nzz < 6) THEN
                     READ(iun_plot, *) wann_func(nxx, nyy, nzz:grid_end(3)-1, ibnd)
                  ELSE
                     READ(iun_plot, '(6E13.5)') wann_func(nxx, nyy, nzz:nzz+5, ibnd)
                  END IF
                ENDDO
              ENDDO
            ENDDO
!            READ(iun_plot, '(6E13.5)') (((wann_func(nxx, nyy, nzz, ibnd),    &
!               nzz = grid_start(3), grid_end(3)), &
!               nyy = grid_start(2), grid_end(2)), &
!               nxx = grid_start(1), grid_end(1))
            CLOSE(iun_plot)
            ! Wannier function is not well normalized
            ! Normalize here will make the calculations with Wannier functions easier
            norm = SUM(wann_func(:, :, :, ibnd) * wann_func(:, :, :, ibnd))
            wann_func(:, :, :, ibnd) = wann_func(:, :, :, ibnd)/SQRT(norm)
         END DO
      END IF
      CALL mp_bcast(wann_func, meta_ionode_id, world_comm)
      CALL mp_bcast(species, meta_ionode_id, world_comm)
      !
   END SUBROUTINE
   !
   SUBROUTINE read_Rp_in_S()
      ! JLB
      !! Read list of Rp unit cell vectors contained on transformed supercell
      !
      USE io_var,        ONLY : iunRpscell
      USE io_global,     ONLY : stdout, ionode, meta_ionode_id
      USE mp,            ONLY : mp_bcast
      USE mp_world,      ONLY : world_comm
      USE elph2,         ONLY : nqtotf
      !
      IMPLICIT NONE
      !
      INTEGER :: iRp, ierr, nRp2
      !
      IF(ionode) THEN
         OPEN(UNIT = iunRpscell, FILE='Rp.scell.plrn', FORM='formatted', STATUS='unknown')
         READ(iunRpscell, *) nRp
         IF (nRp .ne. nqtotf) CALL errore('read_Rp_in_S', 'nRp and nqtotf are not the same!',1)
         CLOSE(UNIT = iunRpscell)
      END IF
      CALL mp_bcast(nRp, meta_ionode_id, world_comm)
      ALLOCATE(Rp(3, nRp), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_Rp_in_S', 'Error allocating Rp', 1)
      Rp = 0
      IF(ionode) THEN
         OPEN(UNIT = iunRpscell, FILE='Rp.scell.plrn', FORM='formatted', STATUS='unknown')
         READ(iunRpscell, *) nRp2
         DO iRp=1,nRp
            READ(iunRpscell, *) Rp(1:3, iRp)
         END DO
         CLOSE(UNIT = iunRpscell)
      END IF
      CALL mp_bcast(Rp, meta_ionode_id, world_comm)
      !
   END SUBROUTINE
   !
   FUNCTION index_Rp(iRp, nqfs)
      USE epwcom, ONLY: nqf1, nqf2, nqf3

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iRp
      INTEGER  :: index_Rp(1:3), nqf_c(1:3)
      INTEGER, INTENT(IN), OPTIONAL  :: nqfs(1:3)


      IF(PRESENT(nqfs)) THEN
         nqf_c(1:3) = nqfs(1:3)
      ELSE
         nqf_c(1:3) = (/nqf1, nqf2, nqf3/)
      END IF

      index_Rp(1) = (iRp - 1)/(nqf_c(2) * nqf_c(3))
      index_Rp(2) = MOD(iRp - 1, nqf_c(2) * nqf_c(3))/nqf_c(3)
      index_Rp(3) = MOD(iRp - 1, nqf_c(3))

      IF(ANY(index_Rp < 0) .or. ANY(index_Rp >= nqf_c)) THEN
         CALL errore('index_Rp','index_Rp not correct!',1)
      END IF
   END FUNCTION
   !
   FUNCTION index_xp(delta_p)
      USE elph2, ONLY : nqtotf
      USE epwcom, ONLY : nqf1, nqf2, nqf3
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: delta_p(1:3)
      INTEGER :: index_xp

      index_xp = delta_p(1) * nqf2 * nqf3 + delta_p(2) * nqf3 + delta_p(3) + 1

      IF(.NOT. ALL(index_Rp(index_xp) == delta_p(1:3))) THEN
         CALL errore('index_xp', 'index_Rp not correct!', 1)
      END IF
   END FUNCTION
   !
   FUNCTION index_shift(ishift)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ishift
      INTEGER  :: index_shift(1:3)

      index_shift(1) = (ishift - 1)/9 - 1
      index_shift(2) = MOD(ishift - 1, 9)/3 - 1
      index_shift(3) = MOD(ishift - 1, 3) - 1

      IF(ANY(index_shift < -1) .or. ANY(index_shift > 1)) THEN
         CALL errore('index_shift', 'index_shift not correct!', 1)
      END IF
   END FUNCTION

END MODULE
