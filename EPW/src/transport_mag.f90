  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  !
  !---------------------------------------------------------------------
  MODULE transport_mag
  !----------------------------------------------------------------------
  !!
  !! This module contains various subroutines when a magnetic field is present
  !!
  !
  CONTAINS
    !-----------------------------------------------------------------------
    SUBROUTINE select_k(inv_tau)
    !-----------------------------------------------------------------------
    !!
    !! This routines selects the k-points that have non-zero inv_tau
    !! We use the following naming convention:
    !!   kpt_ = kpoint coordinate or index
    !!   nkpt_ = number of kpoints
    !!   _bz = full BZ
    !!   _ibz = full IBZ
    !!   _tau = k-point that have inv_tau > eps160
    !!   _fst = k-point within the fsthick (here we consider fsthick including the small buffer)
    !! Dimensions:
    !!   nkpt_ibztau = Number of k-points on the IBZ with non-zero inv_tau
    !!   nkpt_bztau = Number of k-points on the full BZ with non-zero inv_tau
    !!   nkpt_bzfst = Number of k-points within the (extended) fsthick window
    !!
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE input,         ONLY : nstemp
    USE global_var,    ONLY : nbndfst, nktotf, nkpt_ibztau, kpt_ibztau2ibz
    USE ep_constants,  ONLY : eps160
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: inv_tau(nbndfst, nktotf, nstemp)
    !! inv_tau (either inv_tau_all or inv_tau_allcb)
    !
    ! Local variable
    INTEGER :: ik
    !! K-point index
    INTEGER :: itemp
    !! Index over temperature range
    INTEGER :: ierr
    !! Error status
    !
    ! Count the number of k-points with non-zero scattering rates
    ALLOCATE(nkpt_ibztau(nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('select_k', 'Error allocating nkpt_ibztau', 1)
    nkpt_ibztau(:) = 0
    !
    DO itemp = 1, nstemp
      DO ik = 1, nktotf
        ! If any band is non-zero
        IF (SUM(ABS(inv_tau(:, ik, itemp))) > eps160) THEN
          nkpt_ibztau(itemp) = nkpt_ibztau(itemp) + 1
        ENDIF
      ENDDO !ik
    ENDDO !itemp
    !
    ! Create the mapping index for the non-zero inv_tau k-points
    ALLOCATE(kpt_ibztau2ibz(MAXVAL(nkpt_ibztau), nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('select_k', 'Error allocating kpt_ibztau2ibz', 1)
    kpt_ibztau2ibz(:, :) = 0
    nkpt_ibztau(:) = 0
    DO itemp = 1, nstemp
      DO ik = 1, nktotf
        ! If any band is non-zero
        IF (SUM(ABS(inv_tau(:, ik, itemp))) > eps160) THEN
          nkpt_ibztau(itemp) = nkpt_ibztau(itemp) + 1
          kpt_ibztau2ibz(nkpt_ibztau(itemp), itemp) = ik
        ENDIF
      ENDDO !ik
    ENDDO !itemp
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE select_k
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE unfold_k(bztoibz_mat)
    !-----------------------------------------------------------------------
    !!
    !! This routines unfold the nkpt_ibztau k-points from the IBZ to the BZ.
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,       ONLY : DP
    USE input,       ONLY : nstemp, lfast_kmesh, nkf1, nkf2, nkf3
    USE symmetry,    ONLY : nsym => nsym_k
    USE global_var,  ONLY : nktotf, nkpt_ibztau, kpt_ibztau2ibz,       &
                            nkpt_bzfst, kpt_bztau2bz, kpt_bz2bztau,&
                            nkpt_bztau_max, kpt_ibztau2bz
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: bztoibz_mat(nsym, nktotf)
    !! For a given k-point in the IBZ gives the k-point index of all the
    !! k-point in the full BZ that are connected to the current one by symmetry.
    !! nsym + TR is the max number of symmetry
    !
    ! Local variable
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikbz
    !! K-point index in the full BZ
    INTEGER :: itemp
    !! Index over temperature range
    INTEGER :: nb
    !! Index over rotations
    INTEGER :: counter
    !! Counter
    INTEGER :: ierr
    !! Error status
    INTEGER :: nkpt_bztau(nstemp)
    !! Number of kpoints with non-zero inv_tau in full BZ
    !
    ! Count the number of k-points with non-zero inv_tau inside the full BZ
    nkpt_bztau(:) = 0
    !
    ALLOCATE(kpt_ibztau2bz(nsym, MAXVAL(nkpt_ibztau), nstemp), STAT = ierr)
    kpt_ibztau2bz(:,:,:) = 0
    !
    IF (ierr /= 0) CALL errore('unfold_k', 'Error allocating kpt_ibztau2bz', 1)
    ! This is the k-point mapping from full BZ to IBZ of k-points with non-zero inv_tau
    IF (lfast_kmesh) THEN
      ! nkpt_bzfst is the number of k-points in the full BZ within fsthick.
      ALLOCATE(kpt_bz2bztau(nkpt_bzfst, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('unfold_k', 'Error allocating kpt_bz2bztau', 1)
    ELSE
      ALLOCATE(kpt_bz2bztau(nkf1 * nkf2 * nkf3, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('unfold_k', 'Error allocating kpt_bz2bztau', 1)
    ENDIF
    kpt_bz2bztau(:, :) = 0
    !
    DO itemp = 1, nstemp
      DO ik = 1, nkpt_ibztau(itemp)
        ! Loop on the point equivalent by symmetry in the full BZ
        DO nb = 1, nsym
          IF (bztoibz_mat(nb, kpt_ibztau2ibz(ik, itemp)) > 0) THEN
            nkpt_bztau(itemp) = nkpt_bztau(itemp) + 1
            ikbz = bztoibz_mat(nb, kpt_ibztau2ibz(ik, itemp))
            kpt_ibztau2bz(nb, ik, itemp) = ikbz
            kpt_bz2bztau(ikbz, itemp) = nkpt_bztau(itemp)
          ENDIF
        ENDDO ! nb
      ENDDO ! ik
    ENDDO ! itemp
    !
    nkpt_bztau_max = MAXVAL(nkpt_bztau)
    !
    ! We now use the smallest dimension for kpt_bztau2bz
    ALLOCATE(kpt_bztau2bz(nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_k', 'Error allocating kpt_bztau2bz', 1)
    kpt_bztau2bz(:, :) = 0
    !
    DO itemp = 1, nstemp
      counter = 0
      DO ik = 1, nkpt_ibztau(itemp)
        ! Loop on the point equivalent by symmetry in the full BZ
        DO nb = 1, nsym
          IF (bztoibz_mat(nb, kpt_ibztau2ibz(ik, itemp)) > 0) THEN
            counter = counter + 1
            ikbz = bztoibz_mat(nb, kpt_ibztau2ibz(ik, itemp))
            kpt_bztau2bz(counter, itemp) = ikbz
          ENDIF
        ENDDO ! nb
      ENDDO ! ik
    ENDDO ! itemp
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE unfold_k
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE unfold_all(etf_all, vkk_all, inv_tau, f_serta)
    !-----------------------------------------------------------------------
    !!
    !! This routine is used to unfold etf_all, vkk_all, inv_tau and f_serta
    !! from the IBZ (possibly within fsthick) to kpoints within the full BZ but
    !! having non-zero inv_tau.
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,              ONLY : DP
    USE input,              ONLY : nstemp, nkf1, nkf2, nkf3
    USE ep_constants,       ONLY : zero
    USE noncollin_module,   ONLY : noncolin
    USE cell_base,          ONLY : at, bg
    USE symmetry,           ONLY : s => s_k, nsym => nsym_k
    USE global_var,         ONLY : s_bztoibz, nbndfst, nktotf,      &
                                   nkpt_bztau_max, etf_all_b, vkk_all_b, &
                                   wkf_all_b, df_in_b, f_serta_b, f_in_b,    &
                                   f_out_b, inv_tau_b,   &
                                   kpt_bz2bztau, nkpt_ibztau, &
                                   kpt_ibztau2ibz, kpt_ibztau2bz
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity of k
    REAL(KIND = DP), INTENT(in) :: inv_tau(nbndfst, nktotf, nstemp)
    !! inv_tau (either inv_tau_all or inv_tau_allcb)
    REAL(kind = DP), INTENT(in) :: f_serta(3, nbndfst, nktotf, nstemp)
    !! SERTA solution
    !
    ! Local variables
    INTEGER :: ierr
    ! Error status
    INTEGER :: itemp
    !! Index over temperature range
    INTEGER :: nb
    !! Index over rotations
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ik_ibz
    !! K-point index
    INTEGER :: ikbz
    !! K-point index in the BZ
    INTEGER :: ikibz
    !! K-point index in the IBZ
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix (intermediate step)
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in cartesian coordinate

    !
    ALLOCATE(etf_all_b(nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating etf_all_b', 1)
    ALLOCATE(vkk_all_b(3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating vkk_all_b', 1)
    ALLOCATE(wkf_all_b(nkpt_bztau_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating wkf_all_b', 1)
    ALLOCATE(f_serta_b(3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating f_serta_b', 1)
    ALLOCATE(f_in_b(3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating f_in_b', 1)
    ALLOCATE(f_out_b(3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating f_out_b', 1)
    ALLOCATE(inv_tau_b(nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating inv_tau_b', 1)
    ALLOCATE(df_in_b(3, 3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating df_in_cart', 1)
    !
    etf_all_b(:, :, :)    = zero
    vkk_all_b(:, :, :, :) = zero
    wkf_all_b(:)          = zero
    f_serta_b(:, :, :, :) = zero
    f_in_b(:, :, :, :)    = zero
    f_out_b(:, :, :, :)   = zero
    inv_tau_b(:, :, :)    = zero
    df_in_b(:, :, :, :, :) = zero
    !
    IF (noncolin) THEN
      wkf_all_b = 1.0d0 / (nkf1 * nkf2 * nkf3)
    ELSE
      wkf_all_b = 2.0d0 / (nkf1 * nkf2 * nkf3)
    ENDIF
    !
    DO itemp = 1, nstemp
      DO ik_ibz = 1, nkpt_ibztau(itemp)
        ikibz = kpt_ibztau2ibz(ik_ibz, itemp)
        IF(ikibz < 0) CALL errore( 'unfold_all', 'ikibz < 0; it should not happen',1 )
        DO nb = 1, nsym
          ikbz = kpt_ibztau2bz(nb,ik_ibz,itemp)
          IF(ikbz > 0) THEN
            etf_all_b(:, kpt_bz2bztau(ikbz, itemp), itemp) = etf_all(:, ikibz)
            inv_tau_b(:, kpt_bz2bztau(ikbz, itemp), itemp) = inv_tau(:, ikibz, itemp)
            !
            sa(:, :) = DBLE(s(:, :, s_bztoibz(ikbz)))
            sb       = MATMUL(bg, sa)
            sr(:, :) = MATMUL(at, TRANSPOSE(sb))
            sr       = TRANSPOSE(sr)
            !
            DO ibnd = 1, nbndfst
              CALL DGEMV('n', 3, 3, 1.d0, sr, 3, vkk_all(:, ibnd, ikibz), 1, 0.d0, &
                         vkk_all_b(:, ibnd, kpt_bz2bztau(ikbz, itemp), itemp), 1)
              CALL DGEMV('n', 3, 3, 1.d0, sr, 3, f_serta(:, ibnd, ikibz, itemp), 1, 0.d0, &
                         f_serta_b(:, ibnd, kpt_bz2bztau(ikbz, itemp), itemp), 1)
            ENDDO ! ibnd
          ENDIF !ikbz
        ENDDO ! nb
      ENDDO ! ik_ibz
    ENDDO ! itemp
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE unfold_all
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE create_all(etf_all, vkk_all, inv_tau, f_serta)
    !-----------------------------------------------------------------------
    !!
    !! In the case of an homogenous grid without k-point symmetry.
    !! We restrict the kpoints to the one for which the inv_tau is non-zero
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,              ONLY : DP
    USE ep_constants,       ONLY : zero
    USE noncollin_module,   ONLY : noncolin
    USE input,              ONLY : nstemp, nkf1, nkf2, nkf3
    USE global_var,         ONLY : nbndfst, nktotf, nkpt_ibztau,  &
                                   nkpt_bztau_max, etf_all_b, vkk_all_b, &
                                   wkf_all_b, df_in_b, f_serta_b, f_in_b,    &
                                   f_out_b, inv_tau_b, kpt_ibztau2ibz
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies
    REAL(KIND = DP), INTENT(in) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity of k
    REAL(KIND = DP), INTENT(in) :: inv_tau(nbndfst, nktotf, nstemp)
    !! inv_tau (either inv_tau_all or inv_tau_allcb)
    REAL(kind = DP), INTENT(in) :: f_serta(3, nbndfst, nktotf, nstemp)
    !! SERTA solution
    !
    ! Local variables
    INTEGER :: ierr
    ! Error status
    INTEGER :: itemp
    !! Index over temperature range
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikbz
    !! K-point index in the BZ
    !
    nkpt_bztau_max = MAXVAL(nkpt_ibztau)
    !
    ALLOCATE(etf_all_b(nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating etf_all_b', 1)
    ALLOCATE(vkk_all_b(3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating vkk_all_b', 1)
    ALLOCATE(wkf_all_b(nkpt_bztau_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating wkf_all_b', 1)
    ALLOCATE(f_serta_b(3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating f_serta_b', 1)
    ALLOCATE(f_in_b(3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating f_in_b', 1)
    ALLOCATE(f_out_b(3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating f_out_b', 1)
    ALLOCATE(inv_tau_b(nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating inv_tau_b', 1)
    ALLOCATE(df_in_b(3, 3, nbndfst, nkpt_bztau_max, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_all', 'Error allocating df_in_cart', 1)
    !
    etf_all_b(:, :, :)    = zero
    vkk_all_b(:, :, :, :) = zero
    wkf_all_b(:)          = zero
    f_serta_b(:, :, :, :) = zero
    f_in_b(:, :, :, :)    = zero
    f_out_b(:, :, :, :)   = zero
    inv_tau_b(:, :, :)    = zero
    df_in_b(:, :, :, :, :) = zero
    !
    DO itemp = 1, nstemp
      DO ik = 1, nkpt_bztau_max
        ikbz = kpt_ibztau2ibz(ik,itemp)
        etf_all_b(:, ik, itemp) = etf_all(:, ikbz)
        inv_tau_b(:, ik, itemp) = inv_tau(:, ikbz, itemp)
        vkk_all_b(:, :, ik, itemp) = vkk_all(:, :, ikbz)
        f_serta_b(:, :, ik, itemp) = f_serta(:, :, ikbz, itemp)
      ENDDO !ik
    ENDDO !itemp
    !
    IF (noncolin) THEN
      wkf_all_b = 1.0d0 / (nkf1 * nkf2 * nkf3)
    ELSE
      wkf_all_b = 2.0d0 / (nkf1 * nkf2 * nkf3)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE create_all
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE size_indkq(nind, nb_sp, xkf_sp, bztoibz_mat, sparse_q,      &
                          sparse_k, sparse_i, sparse_j, sparse_t, inv_tau, &
                          nkpt_max)
    !-----------------------------------------------------------------------
    !!
    !! This routine create the index map for each ind on that core (nind)
    !! between the (k,q) points associated with that ind index and the
    !! corresponding k (ind1) and k+q (ind2) indexes of the populations f.
    !!
    !! This routine is a dry run with sole purpose to pre-compute nkpt_max
    !! which will then used in create_indkq
    !!
    !! Note: For a given transition probability of index "ind" we may have
    !!       multiple ind1 and ind2 because:
    !!       1) the k-point is a special k-point and needs to be averaged
    !!       2) the k-point is related by symmetry to others and we need to
    !!          computes them explicitely because of B-field removal of
    !!          symmetries.
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE ep_constants,     ONLY : zero, eps160
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg
    USE symmetry,         ONLY : s => s_k, nsym => nsym_k
    USE input,            ONLY : nstemp, lfast_kmesh
    USE global_var,       ONLY : s_bztoibz, nbndfst, nktotf, &
                                 xkf_bz, kpt_bz2bztau, xqf, map_fst
    USE bzgrid,           ONLY : xqf_otf, xkf_otf, kpmq_map
    USE low_lvl,          ONLY : create_interval, bisection
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = 8), INTENT(in) :: nind
    !! Total number of elements per cpu
    INTEGER, INTENT(in) :: nb_sp
    !! Total number of special points
    INTEGER, INTENT(in) :: xkf_sp(97, nb_sp)
    !! Special points indexes and symmetries
    INTEGER, INTENT(in) :: bztoibz_mat(nsym, nktotf)
    !! For a given k-point in the IBZ gives the k-point index of all the
    !! k-point in the full BZ that are connected to the current one by symmetry.
    INTEGER, INTENT(in) :: sparse_q(nind)
    !! Q-point mapping index
    INTEGER, INTENT(in) :: sparse_k(nind)
    !! K-point mapping index
    INTEGER, INTENT(in) :: sparse_i(nind)
    !! Band mapping index
    INTEGER, INTENT(in) :: sparse_j(nind)
    !! Band mapping index
    INTEGER, INTENT(in) :: sparse_t(nind)
    !! Temperature mapping index
    INTEGER, INTENT(inout) :: nkpt_max
    !! Maximum nb of kpoint which have non-zero inv_tau on full BZ corresponding to the ind managed by the current cpu.
    REAL(KIND = DP), INTENT(in) :: inv_tau(nbndfst, nktotf, nstemp)
    !! inv_tau (either inv_tau_all or inv_tau_allcb)
    !
    ! Local variables
    LOGICAL :: special
    !! Is the current k-point a special k-poin
    INTEGER :: ind
    !! Index for sparse matrix written by print_ibte and split on cores
    INTEGER :: ik
    !! K-point index
    INTEGER :: iq
    !! Q-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: jbnd
    !! Local band index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ind1
    !! Index k for that ind
    INTEGER :: ind2
    !! Index k+q for that ind
    INTEGER :: index_sp
    !! K-point index among all the special k-point from [1 ... nb_sp]
    INTEGER :: nb
    !! Index of points in the BZ corresponding to a point in IBZ by a symmetry operation
    INTEGER :: sp
    !! Index of symmetry
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: nkq_abs
    !! Index of the k+q point from the full grid.
    INTEGER :: ierr
    !! Error status
    INTEGER :: n_intval
    !! Number of intervals
    INTEGER, ALLOCATABLE :: val_intval(:)
    !! Value of the first element of each intervals
    INTEGER, ALLOCATABLE :: pos_intval(:)
    !! Position of the first element of each intervals
    REAL(KIND = DP) :: xq(3)
    !! Current q-Cartesian coordinate
    REAL(KIND = DP) :: xk(3)
    !! Current k-Cartesian coordinate
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix (intermediate step)
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in cartesian coordinate
    REAL(KIND = DP) :: sa_sp(3, 3)
    !! Symmetry matrix in crystal for special points
    REAL(KIND = DP) :: sa_tot(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: s_xq(3)
    !! Rotated q-point
    !
    ! Note 1: To find if a k+q point is within the fsthick we need to obtain the mapping
    !         between the index of the point within the fsthick and the index of the point
    !         within the full BZ. This is most efficiently done with bisection.
    ! Note 2: When the number of points within the fshtick window is large, the bissection
    !         is slow. One can speed this up by doing a pre-search since the map_fst is
    !         composed of monotonically increasing numbers (ordered list).
    IF (lfast_kmesh) THEN
      ! We divide map_fst into n_intval intervals
      n_intval = NINT(SQRT(REAL(SIZE(map_fst, 1), KIND = DP)))
      ALLOCATE(val_intval(n_intval), STAT = ierr)
      IF (ierr /= 0) CALL errore('create_indkq', 'Error allocating val_intval', 1)
      ALLOCATE(pos_intval(n_intval), STAT = ierr)
      IF (ierr /= 0) CALL errore('create_indkq', 'Error allocating pos_intval', 1)
      ! We select 1 point every n_interval
      CALL create_interval(SIZE(map_fst, 1), map_fst, n_intval, val_intval, pos_intval)
    ENDIF
    !
    nkpt_max = 0
    DO ind = 1, nind
      iq    = sparse_q(ind)
      ik    = sparse_k(ind)
      ibnd  = sparse_i(ind)
      jbnd  = sparse_j(ind)
      itemp = sparse_t(ind)
      !
      ! Check if k is a special point
      special = .FALSE.
      DO sp = 1, nb_sp
        IF (ik == xkf_sp(1, sp)) THEN
          special = .TRUE.
          index_sp = sp
          EXIT
        ENDIF
      ENDDO
      !
      IF (special) THEN
        DO nb = 1, nsym ! Loops on symmetries.
          IF (bztoibz_mat(nb, ik) > 0) THEN
            ikbz = bztoibz_mat(nb, ik) ! index on the full BZ corresponding to ik in IBZ
            !
            ! Is that k-point such that inv_tau is non-zero ?
            IF (kpt_bz2bztau(ikbz, itemp) > 0) THEN
              IF (lfast_kmesh) THEN
                ! The q-point coordinate is generate on the fly for each q-point
                CALL xqf_otf(iq, xq)
                ! The k-point coordinate is generate on the fly for each k-point
                CALL xkf_otf(map_fst(ikbz), xk)
              ELSE
                xq = xqf(:, iq)
                xk = xkf_bz(:, ikbz)
              ENDIF
              ! Transform the q-point from Crystal to Cartesian coordinate
              CALL cryst_to_cart(1, xq(:), bg, 1)
              !
              ! Do an extra averaging on special k-points liked by symmetry.
              DO sp = 1, nsym
                IF (xkf_sp(sp + 1, index_sp) > 0) THEN
                  sa(:, :)   = DBLE(s(:, :, s_bztoibz(ikbz)))
                  sa_sp(:, :)= DBLE(s(:, :, xkf_sp(sp + 1, index_sp)))
                  sa_tot     = MATMUL(sa, sa_sp)
                  sb         = MATMUL(bg, sa_tot)
                  sr(:, :)   = MATMUL(at, TRANSPOSE(sb))
                  sr         = TRANSPOSE(sr)
                  CALL DGEMV('n', 3, 3, 1.d0, sr, 3, xq(:), 1, 0.d0, s_xq(:), 1)
                  CALL cryst_to_cart(1, s_xq(:), at, -1)
                  CALL kpmq_map(xk, s_xq, 1, nkq_abs)
                  IF (lfast_kmesh) THEN
                    CALL bisection(SIZE(map_fst, 1), map_fst, nkq_abs, n_intval, val_intval, pos_intval)
                    IF (nkq_abs == 0) CYCLE ! The point is not within the fsthick
                  ENDIF ! lfast_kmesh
                  !
                  IF (ABS(inv_tau(ibnd, ik, itemp)) > eps160) THEN
                    ind1 = kpt_bz2bztau(ikbz, itemp)
                    ind2 = kpt_bz2bztau(nkq_abs, itemp)
                    IF (ind1 > 0 .AND. ind2 > 0) THEN
                      nkpt_max = nkpt_max + 1
                    ENDIF
                  ENDIF ! inv_tau
                ENDIF ! xkf_sp
              ENDDO ! sp
            ENDIF ! kpt_bztau2bz
          ENDIF ! bztoibz_mat
        ENDDO ! nb
      ELSE ! Not a special point
        DO nb = 1, nsym ! Loops on symmetries.
          IF (bztoibz_mat(nb, ik) > 0) THEN
            ikbz = bztoibz_mat(nb, ik) ! index on the full BZ corresponding to ik in IBZ
            !
            ! Is that k-point such that inv_tau is non-zero ?
            IF (kpt_bz2bztau(ikbz, itemp) > 0) THEN
              IF (lfast_kmesh) THEN
                ! The q-point coordinate is generate on the fly for each q-point
                CALL xqf_otf(iq, xq)
                ! The k-point coordinate is generate on the fly for each k-point
                CALL xkf_otf(map_fst(ikbz), xk)
              ELSE
                xq = xqf(:, iq)
                xk = xkf_bz(:, ikbz)
              ENDIF
              ! Transform the q-point from Crystal to Cartesian coordinate
              CALL cryst_to_cart(1, xq(:), bg, 1)
              !
              sa(:, :)   = DBLE(s(:, :, s_bztoibz(ikbz)))
              sb         = MATMUL(bg, sa)
              sr(:, :)   = MATMUL(at, TRANSPOSE(sb))
              sr         = TRANSPOSE(sr)
              CALL DGEMV('n', 3, 3, 1.d0, sr, 3, xq(:), 1, 0.d0, s_xq(:), 1)
              CALL cryst_to_cart(1, s_xq(:), at, -1)
              CALL kpmq_map(xk, s_xq, 1, nkq_abs)
              IF (lfast_kmesh) THEN
                CALL bisection(SIZE(map_fst, 1), map_fst, nkq_abs, n_intval, val_intval, pos_intval)
                IF (nkq_abs == 0) CYCLE ! The point is not within the fsthick
              ENDIF ! lfast_kmesh
              !
              IF (ABS(inv_tau(ibnd, ik, itemp)) > eps160) THEN
                ind1 = kpt_bz2bztau(ikbz, itemp)
                ind2 = kpt_bz2bztau(nkq_abs, itemp)
                IF (ind1 > 0 .AND. ind2 > 0) THEN
                  nkpt_max = nkpt_max + 1
                ENDIF
              ENDIF ! inv_tau
            ENDIF ! kpt_bztau2bz
          ENDIF ! bztoibz_mat
        ENDDO ! nb
      ENDIF ! special
    ENDDO ! ind
    !
    WRITE(stdout, '(5x,a,i10)') 'Number of contributing elements for the master core ', nkpt_max
    !
    ! Check that nkpt_max < nind * nsym * nsym
    IF (nkpt_max > nind * nsym * nsym) THEN
      CALL errore('size_indkq', 'nkpt_max cannot be larger than nind * nsym * nsym', 1)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE size_indkq
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE create_indkq(nind, nkpt_max, nb_sp, xkf_sp, bztoibz_mat, sparse_q,      &
                            sparse_k, sparse_i, sparse_j, sparse_t, inv_tau, &
                            nsymk, lsp)
    !-----------------------------------------------------------------------
    !!
    !! This routine create the index map for each ind on that core (nind)
    !! between the (k,q) points associated with that ind index and the
    !! corresponding k (ind1) and k+q (ind2) indexes of the populations f.
    !! Note: For a given transition probability of index "ind" we may have
    !!       multiple ind1 and ind2 because:
    !!       1) the k-point is a special k-point and needs to be averaged
    !!       2) the k-point is related by symmetry to others and we need to
    !!          computes them explicitely because of B-field removal of
    !!          symmetries.
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE ep_constants,     ONLY : zero, eps160
    USE cell_base,        ONLY : at, bg
    USE symmetry,         ONLY : s => s_k, nsym => nsym_k
    USE input,            ONLY : nstemp, lfast_kmesh
    USE global_var,       ONLY : s_bztoibz, nbndfst, nktotf, map_fst, xkf_bz, xqf, &
                                 kpt_bz2bztau, nsym_sp, ind_map
    USE bzgrid,           ONLY : xqf_otf, xkf_otf, kpmq_map
    USE low_lvl,          ONLY : create_interval, bisection
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = 8), INTENT(in) :: nind
    !! Total number of elements per cpu
    INTEGER, INTENT(in) :: nkpt_max
    !! Maximum nb of kpoint which have non-zero inv_tau on full BZ corresponding to the ind managed by the current cpu.
    LOGICAL, INTENT(inout) :: lsp(nind)
    !! At exit, is .true. if the k-point corresponding to the index "ind" is a special point
    INTEGER, INTENT(in) :: nb_sp
    !! Total number of special points
    INTEGER, INTENT(in) :: xkf_sp(97, nb_sp)
    !! Special points indexes and symmetries
    INTEGER, INTENT(in) :: bztoibz_mat(nsym, nktotf)
    !! For a given k-point in the IBZ gives the k-point index of all the
    !! k-point in the full BZ that are connected to the current one by symmetry.
    INTEGER, INTENT(in) :: sparse_q(nind)
    !! Q-point mapping index
    INTEGER, INTENT(in) :: sparse_k(nind)
    !! K-point mapping index
    INTEGER, INTENT(in) :: sparse_i(nind)
    !! Band mapping index
    INTEGER, INTENT(in) :: sparse_j(nind)
    !! Band mapping index
    INTEGER, INTENT(in) :: sparse_t(nind)
    !! Temperature mapping index
    INTEGER, INTENT(inout) :: nsymk(nind)
    !! For a given k-point corresponding to the index "ind", give the number of equivalent kpt by symmetry
    REAL(KIND = DP), INTENT(in) :: inv_tau(nbndfst, nktotf, nstemp)
    !! inv_tau (either inv_tau_all or inv_tau_allcb)
    !
    ! Local variables
    LOGICAL :: special
    !! Is the current k-point a special k-poin
    INTEGER :: ind
    !! Index for sparse matrix written by print_ibte and split on cores
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikpt
    !! Current k-point with max value nkpt_max
    INTEGER :: iq
    !! Q-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: jbnd
    !! Local band index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ind1
    !! Index k for that ind
    INTEGER :: ind2
    !! Index k+q for that ind
    INTEGER :: index_sp
    !! K-point index among all the special k-point from [1 ... nb_sp]
    INTEGER :: nb
    !! Index of points in the BZ corresponding to a point in IBZ by a symmetry operation
    INTEGER :: sp
    !! Index of symmetry
    INTEGER :: isym
    !! Symmetry index
    INTEGER :: jsym
    !! Symmetry index
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: nkq_abs
    !! Index of the k+q point from the full grid.
    INTEGER :: ierr
    !! Error status
    INTEGER :: n_intval
    !! Number of intervals
    INTEGER :: counter_average
    !! Local counter on the number of symmetry related point to a special point.
    INTEGER, ALLOCATABLE :: val_intval(:)
    !! Value of the first element of each intervals
    INTEGER, ALLOCATABLE :: pos_intval(:)
    !! Position of the first element of each intervals
    REAL(KIND = DP) :: xq(3)
    !! Current q-Cartesian coordinate
    REAL(KIND = DP) :: xk(3)
    !! Current k-Cartesian coordinate
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix (intermediate step)
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in cartesian coordinate
    REAL(KIND = DP) :: sa_sp(3, 3)
    !! Symmetry matrix in crystal for special points
    REAL(KIND = DP) :: sa_tot(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: s_xq(3)
    !! Rotated q-point
    !
    ! Note 1: To find if a k+q point is within the fsthick we need to obtain the mapping
    !         between the index of the point within the fsthick and the index of the point
    !         within the full BZ. This is most efficiently done with bisection.
    ! Note 2: When the number of points within the fshtick window is large, the bissection
    !         is slow. One can speed this up by doing a pre-search since the map_fst is
    !         composed of monotonically increasing numbers (ordered list).
    IF (lfast_kmesh) THEN
      ! We divide map_fst into n_intval intervals
      n_intval = NINT(SQRT(REAL(SIZE(map_fst, 1), KIND = DP)))
      ALLOCATE(val_intval(n_intval), STAT = ierr)
      IF (ierr /= 0) CALL errore('create_indkq', 'Error allocating val_intval', 1)
      ALLOCATE(pos_intval(n_intval), STAT = ierr)
      IF (ierr /= 0) CALL errore('create_indkq', 'Error allocating pos_intval', 1)
      ! We select 1 point every n_interval
      CALL create_interval(SIZE(map_fst, 1), map_fst, n_intval, val_intval, pos_intval)
    ENDIF
    !
    ALLOCATE(ind_map(2, nkpt_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('create_indkq', 'Error allocating ind_map', 1)
    ALLOCATE(nsym_sp(nkpt_max), STAT = ierr)
    IF (ierr /= 0) CALL errore('create_indkq', 'Error allocating nsym_sp', 1)
    !
    ikpt = 0
    DO ind = 1, nind
      iq    = sparse_q(ind)
      ik    = sparse_k(ind)
      ibnd  = sparse_i(ind)
      jbnd  = sparse_j(ind)
      itemp = sparse_t(ind)
      !
      ! Check if k is a special point
      special = .FALSE.
      DO sp = 1, nb_sp
        IF (ik == xkf_sp(1, sp)) THEN
          special = .TRUE.
          index_sp = sp
          EXIT
        ENDIF
      ENDDO
      !
      ! Store the list of logical special point for each ind
      lsp(ind) = special
      !
      IF (special) THEN
        isym = 0
        DO nb = 1, nsym ! Loops on symmetries.
          IF (bztoibz_mat(nb, ik) > 0) THEN
            ikbz = bztoibz_mat(nb, ik) ! index on the full BZ corresponding to ik in IBZ
            !
            ! Is that k-point such that inv_tau is non-zero ?
            IF (kpt_bz2bztau(ikbz, itemp) > 0) THEN
              IF (lfast_kmesh) THEN
                ! The q-point coordinate is generate on the fly for each q-point
                CALL xqf_otf(iq, xq)
                ! The k-point coordinate is generate on the fly for each k-point
                CALL xkf_otf(map_fst(ikbz), xk)
              ELSE
                xq = xqf(:, iq)
                xk = xkf_bz(:, ikbz)
              ENDIF
              ! Transform the q-point from Crystal to Cartesian coordinate
              CALL cryst_to_cart(1, xq(:), bg, 1)
              !
              ! Do an extra averaging on special k-points liked by symmetry.
              counter_average = 0
              jsym = 0
              DO sp = 1, nsym
                IF (xkf_sp(sp + 1, index_sp) > 0) THEN
                  counter_average = counter_average + 1
                  sa(:, :)   = DBLE(s(:, :, s_bztoibz(ikbz)))
                  sa_sp(:, :)= DBLE(s(:, :, xkf_sp(sp + 1, index_sp)))
                  sa_tot     = MATMUL(sa, sa_sp)
                  sb         = MATMUL(bg, sa_tot)
                  sr(:, :)   = MATMUL(at, TRANSPOSE(sb))
                  sr         = TRANSPOSE(sr)
                  CALL DGEMV('n', 3, 3, 1.d0, sr, 3, xq(:), 1, 0.d0, s_xq(:), 1)
                  CALL cryst_to_cart(1, s_xq(:), at, -1)
                  CALL kpmq_map(xk, s_xq, 1, nkq_abs)
                  IF (lfast_kmesh) THEN
                    CALL bisection(SIZE(map_fst, 1), map_fst, nkq_abs, n_intval, val_intval, pos_intval)
                    IF (nkq_abs == 0) CYCLE ! The point is not within the fsthick
                  ENDIF ! lfast_kmesh
                  !
                  IF (ABS(inv_tau(ibnd, ik, itemp)) > eps160) THEN
                    ind1 = kpt_bz2bztau(ikbz, itemp)
                    ind2 = kpt_bz2bztau(nkq_abs, itemp)
                    IF (ind1 > 0 .AND. ind2 > 0) THEN
                      isym = isym + 1
                      jsym = jsym + 1
                      ikpt = ikpt + 1
                      ind_map(1, ikpt) = ind1
                      ind_map(2, ikpt) = ind2
                    ENDIF
                  ENDIF ! inv_tau
                ENDIF ! xkf_sp
              ENDDO ! sp
              IF (counter_average > 0) THEN
                DO sp = 1, jsym
                  nsym_sp(ikpt - sp + 1) = counter_average
                ENDDO ! sp
              ENDIF ! counter_average
            ENDIF ! kpt_bztau2bz
          ENDIF ! bztoibz_mat
        ENDDO ! nb
      ELSE ! Not a special point
        isym = 0
        DO nb = 1, nsym ! Loops on symmetries.
          IF (bztoibz_mat(nb, ik) > 0) THEN
            ikbz = bztoibz_mat(nb, ik) ! index on the full BZ corresponding to ik in IBZ
            !
            ! Is that k-point such that inv_tau is non-zero ?
            IF (kpt_bz2bztau(ikbz, itemp) > 0) THEN
              IF (lfast_kmesh) THEN
                ! The q-point coordinate is generate on the fly for each q-point
                CALL xqf_otf(iq, xq)
                ! The k-point coordinate is generate on the fly for each k-point
                CALL xkf_otf(map_fst(ikbz), xk)
              ELSE
                xq = xqf(:, iq)
                xk = xkf_bz(:, ikbz)
              ENDIF
              ! Transform the q-point from Crystal to Cartesian coordinate
              CALL cryst_to_cart(1, xq(:), bg, 1)
              !
              sa(:, :)   = DBLE(s(:, :, s_bztoibz(ikbz)))
              sb         = MATMUL(bg, sa)
              sr(:, :)   = MATMUL(at, TRANSPOSE(sb))
              sr         = TRANSPOSE(sr)
              CALL DGEMV('n', 3, 3, 1.d0, sr, 3, xq(:), 1, 0.d0, s_xq(:), 1)
              CALL cryst_to_cart(1, s_xq(:), at, -1)
              CALL kpmq_map(xk, s_xq, 1, nkq_abs)
              IF (lfast_kmesh) THEN
                CALL bisection(SIZE(map_fst, 1), map_fst, nkq_abs, n_intval, val_intval, pos_intval)
                IF (nkq_abs == 0) CYCLE ! The point is not within the fsthick
              ENDIF ! lfast_kmesh
              !
              IF (ABS(inv_tau(ibnd, ik, itemp)) > eps160) THEN
                ind1 = kpt_bz2bztau(ikbz, itemp)
                ind2 = kpt_bz2bztau(nkq_abs, itemp)
                IF (ind1 > 0 .AND. ind2 > 0) THEN
                  isym = isym + 1
                  ikpt = ikpt + 1
                  ind_map(1, ikpt) = ind1
                  ind_map(2, ikpt) = ind2
                ENDIF
              ENDIF ! inv_tau
            ENDIF ! kpt_bztau2bz
          ENDIF ! bztoibz_mat
        ENDDO ! nb
      ENDIF ! special
      nsymk(ind) = isym
    ENDDO ! ind
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE create_indkq
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE neighbk(map_neigh, cart_der, vec_min)
    !-----------------------------------------------------------------------
    !!
    !! This routines computes and store the 6 nearest neighbours k-point for
    !! all k-points in order to perform finite difference derivative
    !! In addition the routine test if Cartesian derivative is possible or
    !! if crystal derivative has to be used.
    !! In some tilded system with too coarse grids, increment points along the pure Cartesian direction
    !! (x,0,0), (0,y,0) and (0,0,z) cannot be found and one has to use crystal derivative
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE ep_constants,     ONLY : zero, eps6
    USE cell_base,        ONLY : bg
    USE input,            ONLY : nstemp, nkf1, nkf2, nkf3, lfast_kmesh, mp_mesh_k
    USE global_var,       ONLY : map_fst, kpt_bztau2bz, xkf_bz, &
                                 nkpt_bztau_max, kpt_ibztau2ibz
    USE low_lvl,          ONLY : create_interval, bisection
    USE bzgrid,           ONLY : xkf_otf, kpmq_map
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: cart_der
    !! Can we perform Cartesian derivative
    INTEGER, INTENT(inout) :: map_neigh(6, nkpt_bztau_max, nstemp)
    !! Map of the 6 k-point neighbor for k-derivative by finite difference.
    REAL(KIND = DP), INTENT(out) :: vec_min(3)
    !! Norm of the smallest vector in the x, y and z direction
    !
    ! Local variables
    LOGICAL :: found(3)
    !! Found the closest cartesian neighbor in the x, y and z direction
    INTEGER :: xx, yy, zz
    !! Points increment around Gamma in crystal unit
    INTEGER :: x_sign, y_sign, z_sign
    !! Sign of x, y, z (+1 or -1)
    INTEGER :: ik
    !! K-point
    INTEGER :: ikbz
    !! K-points in the full BZ
    INTEGER :: indxp
    !! Index along +x direction
    INTEGER :: indxm
    !! Index along -x direction
    INTEGER :: indyp
    !! Index along +y direction
    INTEGER :: indym
    !! Index along -y direction
    INTEGER :: indzp
    !! Index along +z direction
    INTEGER :: indzm
    !! Index along -z direction
    INTEGER :: ierr
    !! Error status
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: n_intval
    !! Number of intervals
    INTEGER, ALLOCATABLE :: val_intval(:)
    !! Value of the first element of each intervals
    INTEGER, ALLOCATABLE :: pos_intval(:)
    !! Position of the first element of each interval
    REAL(KIND = DP) :: dir_crys(3)
    !! Incremental direction in crystal coordinate
    REAL(KIND = DP) :: dir_cart(3)
    !! Incremental direction in cartesian coordinate
    REAL(KIND = DP) :: xdir_crys_mem(3)
    !! Crystal vector that corresponds to the smallest cartesian vector in the x direction
    REAL(KIND = DP) :: ydir_crys_mem(3)
    !! Crystal vector that corresponds to the smallest cartesian vector in the y direction
    REAL(KIND = DP) :: zdir_crys_mem(3)
    !! Crystal vector that corresponds to the smallest cartesian vector in the z direction
    REAL(KIND = DP) :: xk(3)
    !! Current k-point
    !
    vec_min(:) = 1000.0d0
    found(:) = .FALSE.
    ! First we test if we can do Cartesian derivative (more accurate and symmetric)
    ! Looking for nearest point with pure (x, 0, 0) Cartesian component around Gamma.
    ! We look 4 crystal points around in all directions
    DO xx = -4, 4
      DO yy = -4, 4
        DO zz = -4, 4
          ! We exclude Gamma
          IF (xx == 0 .AND. yy == 0 .AND. zz == 0) CYCLE
          dir_crys = (/REAL(xx) / nkf1, REAL(yy) / nkf2, REAL(zz) / nkf3/)
          CALL DGEMV('n', 3, 3, 1.d0, bg, 3, dir_crys,1 ,0.d0 , dir_cart, 1)
          !
          IF ((ABS(dir_cart(2)) < eps6) .AND. (ABS(dir_cart(3)) < eps6)) THEN
            IF (ABS(dir_cart(1)) < vec_min(1)) THEN
              vec_min(1) = ABS(dir_cart(1))
              xdir_crys_mem(:) = dir_crys(:)
              x_sign = INT(SIGN(1.0d0, dir_cart(1)))
              found(1) = .TRUE.
            ENDIF
          ENDIF ! (x,0,0)
          IF ((ABS(dir_cart(1)) < eps6) .AND. (ABS(dir_cart(3)) < eps6)) THEN
            IF (ABS(dir_cart(2)) < vec_min(2)) THEN
              vec_min(2) = ABS(dir_cart(2))
              ydir_crys_mem(:) = dir_crys(:)
              y_sign = INT(SIGN(1.0d0, dir_cart(2)))
              found(2) = .TRUE.
            ENDIF
          ENDIF ! (0,y,0)
          IF ((ABS(dir_cart(1)) < eps6) .AND. (ABS(dir_cart(2)) < eps6)) THEN
            IF (ABS(dir_cart(3)) < vec_min(3)) THEN
              vec_min(3) = ABS(dir_cart(3))
              zdir_crys_mem(:) = dir_crys(:)
              z_sign = INT(SIGN(1.0d0, dir_cart(3)))
              found(3) = .TRUE.
            ENDIF
          ENDIF ! (0,0,z)
        ENDDO ! zz
      ENDDO ! yy
    ENDDO ! xx
    !
    IF (ANY(found .EQV. .FALSE.)) THEN
      cart_der = .FALSE.
    ENDIF
    !
    IF (lfast_kmesh) THEN
      ! We divide map_fst into n_intval intervals
      n_intval = NINT(SQRT(REAL(SIZE(map_fst, 1), KIND = DP)))
      ALLOCATE(val_intval(n_intval), STAT = ierr)
      IF (ierr /= 0) CALL errore('neighbk', 'Error allocating val_intval', 1)
      ALLOCATE(pos_intval(n_intval), STAT = ierr)
      IF (ierr /= 0) CALL errore('neighbk', 'Error allocating pos_intval', 1)
      ! We select 1 point every n_interval
      CALL create_interval(SIZE(map_fst, 1), map_fst, n_intval, val_intval, pos_intval)
    ENDIF
    !
    ! Now compute and store the nearest neighbor.
    IF (mp_mesh_k) THEN
      DO itemp = 1, nstemp
        DO ik = 1, nkpt_bztau_max
          !
          ikbz = kpt_bztau2bz(ik, itemp)
          IF (ikbz > 0) THEN
            IF (lfast_kmesh) THEN
              CALL xkf_otf(map_fst(ikbz), xk)
              CALL kpmq_map(xk, x_sign * xdir_crys_mem, +1, indxp)
              CALL bisection(SIZE(map_fst, 1), map_fst, indxp, n_intval, val_intval, pos_intval)
              CALL kpmq_map(xk, x_sign * xdir_crys_mem, -1, indxm)
              CALL bisection(SIZE(map_fst, 1), map_fst, indxm, n_intval, val_intval, pos_intval)

              CALL kpmq_map(xk, y_sign * ydir_crys_mem, +1, indyp)
              CALL bisection(SIZE(map_fst, 1), map_fst, indyp, n_intval, val_intval, pos_intval)
              CALL kpmq_map(xk, y_sign * ydir_crys_mem, -1, indym)
              CALL bisection(SIZE(map_fst, 1), map_fst, indym, n_intval, val_intval, pos_intval)

              CALL kpmq_map(xk, z_sign * zdir_crys_mem, +1, indzp)
              CALL bisection(SIZE(map_fst, 1), map_fst, indzp, n_intval, val_intval, pos_intval)
              CALL kpmq_map(xk, z_sign * zdir_crys_mem, -1, indzm)
              CALL bisection(SIZE(map_fst, 1), map_fst, indzm, n_intval, val_intval, pos_intval)
            ELSE ! lfast_kmesh
              xk(:) = xkf_bz(:, ikbz)
              CALL kpmq_map(xk, x_sign * xdir_crys_mem, +1, indxp)
              CALL kpmq_map(xk, x_sign * xdir_crys_mem, -1, indxm)
              CALL kpmq_map(xk, y_sign * ydir_crys_mem, +1, indyp)
              CALL kpmq_map(xk, y_sign * ydir_crys_mem, -1, indym)
              CALL kpmq_map(xk, z_sign * zdir_crys_mem, +1, indzp)
              CALL kpmq_map(xk, z_sign * zdir_crys_mem, -1, indzm)
            ENDIF ! lfast_kmesh
            map_neigh(1, ik, itemp) = indxp
            map_neigh(2, ik, itemp) = indxm
            map_neigh(3, ik, itemp) = indyp
            map_neigh(4, ik, itemp) = indym
            map_neigh(5, ik, itemp) = indzp
            map_neigh(6, ik, itemp) = indzm
          ENDIF ! ikbz
        ENDDO ! ik
      ENDDO ! itemp
    ELSE ! mp_mesh_k
      DO itemp = 1, nstemp
        DO ik = 1, nkpt_bztau_max
          !
          ikbz = kpt_ibztau2ibz(ik, itemp)
          IF (ikbz > 0) THEN
            xk(:) = xkf_bz(:, ikbz)
            CALL kpmq_map(xk, x_sign * xdir_crys_mem, +1, indxp)
            CALL kpmq_map(xk, x_sign * xdir_crys_mem, -1, indxm)
            CALL kpmq_map(xk, y_sign * ydir_crys_mem, +1, indyp)
            CALL kpmq_map(xk, y_sign * ydir_crys_mem, -1, indym)
            CALL kpmq_map(xk, z_sign * zdir_crys_mem, +1, indzp)
            CALL kpmq_map(xk, z_sign * zdir_crys_mem, -1, indzm)
            map_neigh(1, ik, itemp) = indxp
            map_neigh(2, ik, itemp) = indxm
            map_neigh(3, ik, itemp) = indyp
            map_neigh(4, ik, itemp) = indym
            map_neigh(5, ik, itemp) = indzp
            map_neigh(6, ik, itemp) = indzm
          ENDIF ! ikbz
        ENDDO ! ik
      ENDDO ! itemp
    ENDIF ! mp_mesh_k
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE neighbk
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE k_derivative_cart(map_neigh, f_in, vec_min, df_in)
    !-----------------------------------------------------------------------
    !!
    !! Use to compute the Cartesian k-derivative of \partial_E f_nk
    !! General routine that should work for any system
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE input,         ONLY : nstemp, lfast_kmesh, mp_mesh_k
    USE control_flags, ONLY : iverbosity
    USE global_var,    ONLY : nbndfst, nkpt_bztau_max, gtemp, &
                              kpt_bz2bztau
    USE ep_constants,  ONLY : pi, zero, ryd2ev, kelvin2eV
    USE cell_base,     ONLY : alat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: map_neigh(6, nkpt_bztau_max, nstemp)
    !! Map of the 6 k-point neighbor for k-derivative by finite difference.
    REAL(KIND = DP), INTENT(in) :: f_in(3, nbndfst, nkpt_bztau_max, nstemp)
    !! In population for iteration i
    REAL(KIND = DP), INTENT(in) :: vec_min(3)
    !! Norm of the smallest vector in the x, y and z direction
    REAL(KIND = DP), INTENT(inout) :: df_in(3, 3, nbndfst, nkpt_bztau_max, nstemp)
    !! In solution for iteration i, derived
    !! This is ALL the k-points from the full BZ (not only IBZ)
    !
    ! Local variables
    INTEGER :: idex(6)
    !! Index along the 6 directions: idex(1,2,3,..) corresponds to +x, -x, +y, -y, +z and -z direction, respectively
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: ind1
    !! Index on the IBZ
    INTEGER :: ii
    !! Direction index
    REAL(KIND = DP) :: f(3, nbndfst, 6)
    !! Local value of f along the 6 directions
    REAL(KIND = DP) :: sumrule(3, 3)
    !! Variable for sum rule
    REAL(KIND = DP) :: alat_inv
    !! Reciprocal of alat
    !
    alat_inv = (2.0d0 * pi / alat)
    idex(:) = 0
    f(:, :, :) = zero
    IF (mp_mesh_k) THEN
      DO itemp = 1, nstemp
        DO ik = 1, nkpt_bztau_max
          idex(:) = map_neigh(:, ik, itemp)
          !
          DO ii = 1, 6
            IF (idex(ii) == 0) THEN
              IF (.NOT. lfast_kmesh) CALL errore('k_derivative_cart', &
                  'The index cannot be 0 with lfast_kmesh == .FALSE.', 1)
              f(:, :, ii) = zero
            ELSE
              IF (kpt_bz2bztau(idex(ii), itemp) == 0) THEN
                f(:, :, ii) = zero
              ELSE
                ind1        = kpt_bz2bztau(idex(ii), itemp)
                f(:, :, ii) = f_in(:, :, ind1, itemp)
              ENDIF
            ENDIF
          ENDDO
          !
          DO ibnd = 1, nbndfst
            df_in(1, :, ibnd, ik, itemp) = (f(:, ibnd, 1) - f(:, ibnd, 2)) / (2.0d0 * vec_min(1) * alat_inv)
            df_in(2, :, ibnd, ik, itemp) = (f(:, ibnd, 3) - f(:, ibnd, 4)) / (2.0d0 * vec_min(2) * alat_inv)
            df_in(3, :, ibnd, ik, itemp) = (f(:, ibnd, 5) - f(:, ibnd, 6)) / (2.0d0 * vec_min(3) * alat_inv)
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp
    ELSE ! mp_mesh_k
      DO itemp = 1, nstemp
        DO ik = 1, nkpt_bztau_max
          idex(:) = map_neigh(:, ik, itemp)
          DO ii = 1, 6
            IF (idex(ii) > 0) f(:, :, ii) = f_in(:, :, idex(ii), itemp)
          ENDDO
          DO ibnd = 1, nbndfst
            df_in(1, :, ibnd, ik, itemp) = (f(:, ibnd, 1) - f(:, ibnd, 2)) / (2.0d0 * vec_min(1) * alat_inv)
            df_in(2, :, ibnd, ik, itemp) = (f(:, ibnd, 3) - f(:, ibnd, 4)) / (2.0d0 * vec_min(2) * alat_inv)
            df_in(3, :, ibnd, ik, itemp) = (f(:, ibnd, 5) - f(:, ibnd, 6)) / (2.0d0 * vec_min(3) * alat_inv)
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp
    ENDIF ! mp_mesh_k
    !
    ! Checking that population sum rule
    ! Note: We only check the first temperature !
    sumrule(1, 1) = SUM(df_in(1, 1, :, :, 1)) / nkpt_bztau_max
    sumrule(1, 2) = SUM(df_in(1, 2, :, :, 1)) / nkpt_bztau_max
    sumrule(1, 3) = SUM(df_in(1, 3, :, :, 1)) / nkpt_bztau_max
    sumrule(2, 1) = SUM(df_in(2, 1, :, :, 1)) / nkpt_bztau_max
    sumrule(2, 2) = SUM(df_in(2, 2, :, :, 1)) / nkpt_bztau_max
    sumrule(2, 3) = SUM(df_in(2, 3, :, :, 1)) / nkpt_bztau_max
    sumrule(3, 1) = SUM(df_in(3, 1, :, :, 1)) / nkpt_bztau_max
    sumrule(3, 2) = SUM(df_in(3, 2, :, :, 1)) / nkpt_bztau_max
    sumrule(3, 3) = SUM(df_in(3, 3, :, :, 1)) / nkpt_bztau_max
    !
    IF (iverbosity == 4) THEN
      WRITE(stdout, '(5x, a)') ' '
      WRITE(stdout, '(5x, a)') 'Cartesian derivative [more accurate]'
      WRITE(stdout, '(5x, a, 1f10.4, a)') 'Checking the sum rule on df/dk for temperature ', &
                                           gtemp(1) * ryd2ev / kelvin2eV,' K'
      WRITE(stdout, '(5x, a)') '[For cubic materials, in-diagonal terms should go to 0 with increasing grids and fsthick; &
                               &off diagonal terms should be 0 without magnetic fields] '
      WRITE(stdout, '(5x, 1E18.6, 1E18.6, 1E18.6)') sumrule(1, 1), sumrule(2, 1), sumrule(3, 1)
      WRITE(stdout, '(5x, 1E18.6, 1E18.6, 1E18.6)') sumrule(1, 2), sumrule(2, 2), sumrule(3, 2)
      WRITE(stdout, '(5x, 1E18.6, 1E18.6, 1E18.6)') sumrule(1, 3), sumrule(2, 3), sumrule(3, 3)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE k_derivative_cart
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE k_derivative_crys(map_neigh, f_in, df_in)
    !-----------------------------------------------------------------------
    !!
    !! Use to compute the Cartesian k-derivative of \partial_E f_nk
    !! General routine that should work for any system
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE input,         ONLY : nstemp, lfast_kmesh, mp_mesh_k, nkf1, nkf2, nkf3
    USE control_flags, ONLY : iverbosity
    USE global_var,    ONLY : nbndfst, nkpt_bztau_max, gtemp, &
                              kpt_bz2bztau
    USE ep_constants,  ONLY : pi, zero, ryd2ev, kelvin2eV
    USE cell_base,     ONLY : alat, at
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: map_neigh(6, nkpt_bztau_max, nstemp)
    !! Map of the 6 k-point neighbor for k-derivative by finite difference.
    REAL(KIND = DP), INTENT(in) :: f_in(3, nbndfst, nkpt_bztau_max, nstemp)
    !! In population for iteration i
    REAL(KIND = DP), INTENT(inout) :: df_in(3, 3, nbndfst, nkpt_bztau_max, nstemp)
    !! In solution for iteration i, derived
    !! This is ALL the k-points from the full BZ (not only IBZ)
    !
    ! Local variables
    INTEGER :: idex(6)
    !! Index along the 6 directions: idex(1,2,3,..) corresponds to +x, -x, +y, -y, +z and -z direction, respectively
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ik
    !! K-point index
    INTEGER :: i
    !! Cartesian direction
    INTEGER :: ii
    !! Direction index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: ind1
    !! Index on the IBZ
    REAL(KIND = DP) :: f(3, nbndfst, 6)
    !! Local value of f along the 6 directions
    REAL(KIND = DP) :: sumrule(3, 3)
    !! Variable for sum rule
    REAL(KIND = DP) :: alat_inv
    !! Reciprocal of alat
    REAL(KIND = DP) :: dx, dy, dz
    !! Increment in crystal unit
    !
    dx = 1.0 / nkf1
    dy = 1.0 / nkf2
    dz = 1.0 / nkf3
    alat_inv = (2.0d0 * pi / alat)
    idex(:) = 0
    f(:, :, :) = zero
    IF (mp_mesh_k) THEN
      DO itemp = 1, nstemp
        DO ik = 1, nkpt_bztau_max
          idex(:) = map_neigh(:, ik, itemp)
          !
          DO ii = 1, 6
            IF (idex(ii) == 0) THEN
              IF (.NOT. lfast_kmesh) CALL errore('k_derivative_crys', &
                  'The index cannot be 0 with lfast_kmesh == .FALSE.', 1)
              f(:, :, ii) = zero
            ELSE
              IF (kpt_bz2bztau(idex(ii), itemp) == 0) THEN
                f(:, :, ii) = zero
              ELSE
                ind1        = kpt_bz2bztau(idex(ii), itemp)
                f(:, :, ii) = f_in(:, :, ind1, itemp)
              ENDIF
            ENDIF
          ENDDO
          !
          DO ibnd = 1, nbndfst
            df_in(1, :, ibnd, ik, itemp) = (f(:, ibnd, 1) - f(:, ibnd, 2)) / (2.0d0 * dx * alat_inv)
            df_in(2, :, ibnd, ik, itemp) = (f(:, ibnd, 3) - f(:, ibnd, 4)) / (2.0d0 * dy * alat_inv)
            df_in(3, :, ibnd, ik, itemp) = (f(:, ibnd, 5) - f(:, ibnd, 6)) / (2.0d0 * dz * alat_inv)
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp
    ELSE ! mp_mesh_k
      DO itemp = 1, nstemp
        DO ik = 1, nkpt_bztau_max
          idex(:) = map_neigh(:, ik, itemp)
          DO ii = 1, 6
            IF (idex(ii) > 0) f(:, :, ii) = f_in(:, :, idex(ii), itemp)
          ENDDO
          DO ibnd = 1, nbndfst
            DO i = 1, 3
              df_in(1, i, ibnd, ik, itemp) = (f(i, ibnd, 1) - f(i, ibnd, 2)) / (2.0d0 * dx * alat_inv)
              df_in(2, i, ibnd, ik, itemp) = (f(i, ibnd, 3) - f(i, ibnd, 4)) / (2.0d0 * dy * alat_inv)
              df_in(3, i, ibnd, ik, itemp) = (f(i, ibnd, 5) - f(i, ibnd, 6)) / (2.0d0 * dz * alat_inv)
              !
              !Transform derivative with rotation
              df_in(:, i, ibnd, ik, itemp) = MATMUL(at(:, :), df_in(:, i, ibnd, ik, itemp))
            ENDDO ! i
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! itemp
    ENDIF ! mp_mesh_k
    !
    ! Checking that population sum rule
    ! Note: We only check the first temperature !
    sumrule(1, 1) = SUM(df_in(1, 1, :, :, 1)) / nkpt_bztau_max
    sumrule(1, 2) = SUM(df_in(1, 2, :, :, 1)) / nkpt_bztau_max
    sumrule(1, 3) = SUM(df_in(1, 3, :, :, 1)) / nkpt_bztau_max
    sumrule(2, 1) = SUM(df_in(2, 1, :, :, 1)) / nkpt_bztau_max
    sumrule(2, 2) = SUM(df_in(2, 2, :, :, 1)) / nkpt_bztau_max
    sumrule(2, 3) = SUM(df_in(2, 3, :, :, 1)) / nkpt_bztau_max
    sumrule(3, 1) = SUM(df_in(3, 1, :, :, 1)) / nkpt_bztau_max
    sumrule(3, 2) = SUM(df_in(3, 2, :, :, 1)) / nkpt_bztau_max
    sumrule(3, 3) = SUM(df_in(3, 3, :, :, 1)) / nkpt_bztau_max
    !
    IF (iverbosity == 4) THEN
      WRITE(stdout, '(5x, a)') ' '
      WRITE(stdout, '(5x, a)') 'Crystal derivative [less accurate]'
      WRITE(stdout, '(5x, a, 1f10.4, a)') 'Checking the sum rule on df/dk for temperature ', &
                                           gtemp(1) * ryd2ev / kelvin2eV,' K'
      WRITE(stdout, '(5x, a)') '[For cubic materials, in-diagonal terms should go to 0 with increasing grids and fsthick; &
                               &off diagonal terms should be 0 without magnetic fields] '
      WRITE(stdout, '(5x, 1E18.6, 1E18.6, 1E18.6)') sumrule(1, 1), sumrule(2, 1), sumrule(3, 1)
      WRITE(stdout, '(5x, 1E18.6, 1E18.6, 1E18.6)') sumrule(1, 2), sumrule(2, 2), sumrule(3, 2)
      WRITE(stdout, '(5x, 1E18.6, 1E18.6, 1E18.6)') sumrule(1, 3), sumrule(2, 3), sumrule(3, 3)
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE k_derivative_crys
    !-----------------------------------------------------------------------
  !----------------------------------------------------------------------
  END MODULE transport_mag
  !----------------------------------------------------------------------
