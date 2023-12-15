  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE rotate_wavefunction
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines to perform symmetry rotations of
  !! the wavefunctions.
  !!
  !! Partially adapted from elphel2_shuffle.f90
  !!
  !----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  !
  IMPLICIT NONE
  !
  TYPE twodarray
    COMPLEX(DP), ALLOCATABLE :: arr(:, :)
  END TYPE twodarray
  !
  INTEGER :: ntrev
  !! 2 if time-reversal symmetric (time_reversal == .true.), 1 if not.
  LOGICAL, ALLOCATABLE :: exband_rotate(:)
  !! k-point independent list of bands excluded from the calculation
  !! of overlap and projection matrices in W90
  INTEGER, ALLOCATABLE :: degen_grp_inds(:, :)
  !! (nbndep, nkstot) Index of degenerate band groups. Degenerate bands belong
  !! to the same group
  INTEGER, ALLOCATABLE :: degen_grp_ndim(:, :)
  !! (nbndep, nkstot) Dimension of each degenerate band group.
  INTEGER, ALLOCATABLE :: num_degen_grps(:)
  !! (nkstot) Number of groups for each k point.
  INTEGER, ALLOCATABLE :: sym_ktok(:, :, :)
  !! (nkstot, nsym, ntrev) sym_ktok(ik, isym) is the index of S(k(ik))
  INTEGER :: max_num_grps
  !! max(num_degen_grps)
  INTEGER, SAVE :: iq_first_save
  !! Index of the q-point saved in sthmatq_save.
  REAL(DP), ALLOCATABLE :: g0vec_sym(:, :, :, :)
  !! (3, nkstot, nsym, ntrev) G0 in Cartesian representation such that Sk + G0 belongs
  !! to the first BZ.
  COMPLEX(DP), ALLOCATABLE :: emiskv(:, :, :)
  !! (nkstot, nsym, ntrev) exp(-i S(k) v) for fractional translation.
  !
  COMPLEX(DP), ALLOCATABLE :: dmat_all(:, :, :, :, :)
  !! (nbnd, nbnd, nks, nsym, ntrev) <psi_S(k)|psi_k>
  ! TYPE(twodarray), ALLOCATABLE :: dmat(:, :, :)
  !! (max_num_grps, nkstot, nsym) block diagonal matrix <psi_S(k)|psi_k>
  !! dmat(igrp, ik, isym)%arr is square matrix with dimension degen_grp_ndim(igrp, ik)
  COMPLEX(DP), ALLOCATABLE, SAVE :: sthmatq_save(:, :, :, :, :)
  !! sthmatq at irreducible q point.
  !
  CONTAINS
    !
    ! -----------------------------------------------------------------------
    SUBROUTINE setup_rotate_wavefunction(nsym, s, ft, invs)
    !-----------------------------------------------------------------------
    !!
    !! Group energies into degenerate groups and setup degen_grp_inds,
    !! degen_grp_ndim, and num_degen_grps.
    !!
    !! The cutoff for degeneracy, degen_cutoff, is set to eps4.
    !!
    !----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : zero, eps4, eps6, twopi, ci, cone
    USE io_global,     ONLY : ionode_id, meta_ionode
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_bcast
    USE pwcom,         ONLY : nkstot, nbnd
    USE cell_base,     ONLY : at
    USE symm_base,     ONLY : time_reversal
    USE epwcom,        ONLY : filukk, nbndsub
    USE io_var,        ONLY : iunukk
    USE elph2,         ONLY : nbndep
    USE klist_epw,     ONLY : et_all, xk_all
    USE low_lvl,       ONLY : s_crystocart, rotate_cart
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nsym
    !! the number of symmetries of the crystal
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! the symmetry matrices
    REAL(KIND = DP), INTENT(in) :: ft(3, 48)
    !! the fractional traslations in crystal axis
    INTEGER, INTENT(in) :: invs(48)
    !! inverse symmetry matrix
    !
    ! Local variables
    LOGICAL :: dummy_logical
    !! for skipping unnecessary data in ukk file
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ik
    !! k point index
    INTEGER :: jk
    !! k point index
    INTEGER :: isym
    !! symmetry operation index
    INTEGER :: itrev
    !! 1 for symmetry operation without time reversal, 2 with time versal
    INTEGER :: ibnd_all
    !! band index
    INTEGER :: ibnd
    !! band index
    INTEGER :: jbnd
    !! band index
    INTEGER :: igrp
    !! degenerate group index
    INTEGER :: idir
    !! Cartesian direction index
    INTEGER :: ierr
    !! Error index
    REAL(DP) :: rdotk
    !! phase factor for fractional translation
    REAL(DP) :: et_prev
    !! band energy of previous band
    REAL(DP) :: degen_cutoff = eps4
    !! Cutoff for degeneracy
    REAL(DP) :: ft_cart(3)
    !! fractional translation in Cartesian coordinate
    REAL(DP) :: xk_diff(3)
    !! k vector
    REAL(DP) :: xk_diff_crys(3)
    !! k vector in crystal coordinate
    REAL(DP) :: sxkvec(3)
    !! Sk vector
    COMPLEX(DP) :: dummy_complex
    !! for skipping unnecessary data in ukk file
    !
    ntrev = 1
    IF (time_reversal) ntrev = 2
    !
    ALLOCATE(degen_grp_inds(nbndep, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error allocating degen_grp_inds', 1)
    ALLOCATE(degen_grp_ndim(nbndep, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error allocating degen_grp_ndim', 1)
    ALLOCATE(num_degen_grps(nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error allocating num_degen_grps', 1)
    ALLOCATE(sym_ktok(nkstot, nsym, ntrev), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error allocating sym_ktok', 1)
    ALLOCATE(exband_rotate(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error allocating exband_rotate', 1)
    ALLOCATE(g0vec_sym(3, nkstot, nsym, ntrev), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error allocating g0vec_sym', 1)
    ALLOCATE(emiskv(nkstot, nsym, ntrev), STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error allocating emiskv', 1)
    !
    degen_grp_inds(:, :) = 0
    degen_grp_ndim(:, :) = 0
    num_degen_grps(:) = 0
    g0vec_sym(:, :, :, :) = zero
    emiskv(:, :, :) = cone
    !
    ! Compute sym_ktok: index of S(k) in xk_all
    !
    DO ik = 1, nkstot
      !
      DO itrev = 1, ntrev
        !
        DO isym = 1, nsym
          ! sxkvec S(k) = (k*s) if itrev == 1, -(k*s) if itrev == 2
          ! rotate_cart applies s(invs(isym)).
          ! To apply s, one needs MATMUL(k, s), but rotate_cart computes MATMUL(s, k).
          CALL rotate_cart(xk_all(:, ik), s(:, :, invs(isym)), sxkvec)
          IF (itrev == 2) sxkvec = -sxkvec
          !
          ! Find jk such that xk_all(:, jk) == sxkvec + G (mod G)
          !
          DO jk = 1, nkstot
            xk_diff = xk_all(:, jk) - sxkvec
            !
            xk_diff_crys(:) = zero
            DO idir = 1, 3
              xk_diff_crys(:) = xk_diff_crys(:) + at(idir, :) * xk_diff(idir)
            ENDDO
            !
            ! Check xk_all(:, jk) == S(k) + G0
            !
            IF ( ABS(REAL(NINT(xk_diff_crys(1)), DP) - xk_diff_crys(1)) < eps6 .AND. &
                 ABS(REAL(NINT(xk_diff_crys(2)), DP) - xk_diff_crys(2)) < eps6 .AND. &
                 ABS(REAL(NINT(xk_diff_crys(3)), DP) - xk_diff_crys(3)) < eps6 ) THEN
              sym_ktok(ik, isym, itrev) = jk
              !
              g0vec_sym(:, ik, isym, itrev) = xk_diff(:)
              EXIT
            ENDIF
            !
            ! If match is not found, raise error.
            !
            IF (jk == nkstot) THEN
              CALL errore('setup_rotate_wavefunction', &
                  'k2 such that k2 = S(k) + G not found', 1)
            ENDIF
            !
          ENDDO
          !
          ! Setup phase factor emiskv = exp(-i S(k) v).
          ! The actual translation is v = -ft. (see ruotaijk.f90)
          ! So the overall sign is +: emiskv = exp(i S(k) ft).
          !
          IF (ft(1, isym)**2 + ft(2, isym)**2 + ft(3, isym)**2 > 1.0d-8) THEN
            ft_cart(:) = ft(:, isym)
            CALL cryst_to_cart(1, ft_cart, at, +1)
            !
            rdotk = SUM(sxkvec(1:3) * ft_cart(1:3))
            emiskv(ik, isym, itrev) = EXP(ci * twopi * rdotk)
          ELSE
            emiskv(ik, isym, itrev) = cone
          ENDIF
          !
        ENDDO ! isym
      ENDDO ! itrev
    ENDDO ! ik
    !
    ! Read ukk file to setup exband_rotate
    !
    IF (meta_ionode) THEN
      !
      OPEN(iunukk, FILE = filukk, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('setup_rotate_wavefunction', 'error opening ukk file', iunukk)
      !
      ! Skip unnecessary data
      READ(iunukk, *) ibnd, jbnd
      DO ik = 1, nkstot * nbndep * nbndsub
         READ(iunukk, *) dummy_complex
      ENDDO
      DO ik = 1, nkstot * nbndep
        READ(iunukk, *) dummy_logical
      ENDDO
      !
      DO ibnd = 1, nbnd
        READ(iunukk, *) exband_rotate(ibnd)
      ENDDO
      !
      CLOSE(iunukk)
    ENDIF ! meta_ionode
    !
    CALL mp_bcast(exband_rotate, ionode_id, inter_pool_comm)
    !
    ! Check band energy of symmetry equivalent k points are equal
    !
    DO ik = 1, nkstot
      DO itrev = 1, ntrev
        DO isym = 1, nsym
          jk = sym_ktok(ik, isym, itrev)
          DO ibnd = 1, nbnd
            IF (exband_rotate(ibnd)) CYCLE
            IF ( ABS(et_all(ibnd, ik) - et_all(ibnd, jk)) > eps4 ) THEN
              CALL errore('setup_rotate_wavefunction', 'symmetry equivalent &
                &k points have different energy. Use unfold_ksym = .false.', 1)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! isym
      ENDDO ! itrev
    ENDDO ! ik
    !
    ! Loop over k points and bands to find degeneracy
    !
    DO ik = 1, nkstot
      igrp = 1
      ibnd = 0
      DO ibnd_all = 1, nbnd
        !
        IF (exband_rotate(ibnd_all)) CYCLE
        ibnd = ibnd + 1
        !
        IF (ibnd == 1) THEN
          degen_grp_inds(ibnd, ik) = igrp
          degen_grp_ndim(igrp, ik) = degen_grp_ndim(igrp, ik) + 1
          et_prev = et_all(ibnd_all, ik)
          !
        ELSE ! ibnd /= 1
          ! IF not degenerate, increase group index
          IF ( ABS(et_all(ibnd_all, ik) - et_prev) > degen_cutoff ) THEN
            igrp = igrp + 1
          ENDIF
          !
          degen_grp_inds(ibnd, ik) = igrp
          degen_grp_ndim(igrp, ik) = degen_grp_ndim(igrp, ik) + 1
          et_prev = et_all(ibnd_all, ik)
          !
        ENDIF
        !
      ENDDO ! ibnd_all
      !
      num_degen_grps(ik) = igrp
      !
    ENDDO ! ik
    !
    ! Check degen_grp_inds of symmetry equivalent k points are equal
    !
    DO ik = 1, nkstot
      DO itrev = 1, ntrev
        DO isym = 1, nsym
          jk = sym_ktok(ik, isym, itrev)
          IF ( .NOT. ALL(degen_grp_inds(:, ik) == degen_grp_inds(:, jk)) ) THEN
            CALL errore('setup_rotate_wavefunction', 'symmetry equivalent &
              &k points have different degen_grp_inds. Use unfold_ksym = .false.', 1)
          ENDIF
        ENDDO ! isym
      ENDDO ! itrev
    ENDDO ! ik
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE setup_rotate_wavefunction
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE calc_rotation_gauge(nsym, s, invs, gmapsym, eigv, lwin)
    !-----------------------------------------------------------------------
    !! Compute wavefunction gauge matrix <psi(Sk)|S|psi(k)>
    !!
    !! Adapted from elphel2_shuffle.f90
    !!
    !! Currently, developmental stage: compute full matrix dmat_all
    !! Later, TODO: compute only block-diagonal part btw degenerate states
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : zero, czero, cone, eps6
    USE mp,            ONLY : mp_barrier
    USE mp_global,     ONLY : my_pool_id, inter_pool_comm, world_comm
    USE io_global,     ONLY : stdout, meta_ionode
    USE io_files,      ONLY : diropn, prefix
    USE wavefunctions, ONLY : evc
    USE fft_base,      ONLY : dffts
    USE pwcom,         ONLY : nbnd, nks, nkstot
    USE cell_base,     ONLY : at, bg
    USE wvfct,         ONLY : npwx
    USE noncollin_module, ONLY : npol, noncolin
    USE units_lr,      ONLY : iuwfc, lrwfc
    USE io_epw,        ONLY : readwfc, readgmap
    USE io_var,        ONLY : iudmat
    USE elph2,         ONLY : nbndep, gmap, ngk_all, igk_k_all, ibndstart, &
                              ng0vec, g0vec_all_r, ngxxf
    USE division,      ONLY : fkbounds
    USE kfold,         ONLY : ktokpmq, createkmap
    USE klist_epw,     ONLY : xk_loc, xk_all
    USE low_lvl,       ONLY : fractrasl, rotate_cart, s_crystocart
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: lwin(nbndep, nks)
    !! Bands at k within outer energy window
    INTEGER, INTENT(in) :: nsym
    !! the number of symmetries of the crystal
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! the symmetry matrices
    INTEGER, INTENT(in) :: invs(48)
    !! inverse symmetry matrix
    INTEGER, INTENT(in) :: gmapsym(ngxxf, 48)
    !! Correspondence G -> S(G)
    COMPLEX(KIND = DP), INTENT(in) :: eigv(ngxxf, 48)
    !! $e^{iGv}$ for 1 ... nsym (v the fractional translation)
    !
    LOGICAL :: exst
    !! logical variable to check file exists
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: isym
    !! symmetry operation index
    INTEGER :: itrev
    !! 1 for symmetry operation without time reversal, 2 with time versal
    INTEGER :: ig
    !! plane wave index
    INTEGER :: ig0
    !! plane wave index
    INTEGER :: shift_sym
    !! index of G0 in g0vec_all_r such that Sk + G0 is in the first BZ
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: ib_win_fst
    !! First band inside the outer window
    INTEGER :: ib_win_lst
    !! First band inside the outer window
    INTEGER :: imap
    !! Index in gmap for G-sphere igkq translation with G_0
    INTEGER :: ik
    !! k point index
    INTEGER :: ik_gl
    !! global k point index
    INTEGER :: isk_gl
    !! global k point index of Sk
    INTEGER :: isk_gl_
    !! global k point index of Sk
    INTEGER :: ipooltmp
    !! Index of pool for k
    INTEGER :: ipool
    !! Index of pool for Sk
    INTEGER :: npw
    !! Number of k+G-vectors inside 'ecut sphere'
    INTEGER :: npwsk
    !! Number of Sk+G-vectors inside 'ecut sphere'
    INTEGER :: nsk
    !! Index of Sk-point in the pool
    INTEGER :: irec
    !! Record number
    INTEGER :: lrdmat
    !! Size of dmat file
    INTEGER :: iusymk
    !! Unit for reading and writing symk file
    INTEGER :: ierr
    !! Error index
    INTEGER :: igsk_tmp(npwx)
    !! Correspondence Sk+G <-> G
    INTEGER, ALLOCATABLE :: igk(:)
    !! Index for k+G
    INTEGER, ALLOCATABLE :: igsk(:)
    !! Index for k+q+G
    REAL(KIND = DP) :: xk_diff(3)
    !! k-point vector S(k) - k
    REAL(KIND = DP) :: s_cart(3, 3)
    !! Symmetry matrix in Cartesian basis.
    COMPLEX(KIND = DP) :: su2(2, 2)
    !! SU2 rotation matrix to act WFs
    COMPLEX(KIND = DP) :: su2_no_dagger(2, 2)
    !! SU2 that comes out of find_u, we need the dagger of this, which we assign to su2
    COMPLEX(KIND = DP), ALLOCATABLE :: dmat_temp(:, :)
    !! Temporary storage for dmat_all
    COMPLEX(KIND = DP), ALLOCATABLE :: evsk(:, :)
    !! wavefunction at S(k)
    COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:, :)
    !! Auxillary wavefunction
    COMPLEX(KIND = DP), ALLOCATABLE :: aux2(:, :)
    !! Auxillary wavefunction
    COMPLEX(KIND = DP), ALLOCATABLE :: aux3(:, :)
    !! Auxillary wavefunction
    COMPLEX(KIND = DP) :: aux4(npwx * npol, nbnd)
    !! Rotated psi_m,Sk WF by SU(2)
    COMPLEX(KIND = DP) :: aux5(npwx * npol, nbnd)
    !! Rotated psi_n,k WF by SU(2)
    !
    WRITE(stdout, '(5x,a)') 'Compute wavefunction overlap <psi(Sk)|S|psi(k)>'
    !
    ALLOCATE(evsk(npwx * npol, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating evsk', 1)
    ALLOCATE(aux1(dffts%nnr, npol), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating aux1', 1)
    ALLOCATE(aux2(dffts%nnr, npol), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating aux2', 1)
    ALLOCATE(aux3(npwx * npol, nbndep), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating aux3', 1)
    !
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)
    !
    ! close all sequential files in order to re-open them as direct access
    ! close all .wfc files in order to prepare shuffled read
    !
    CLOSE(iuwfc, STATUS = 'keep')
    ! never remove this barrier
    CALL mp_barrier(inter_pool_comm)
    !
    ! First, do necessary allocations
    !
    ! BEGIN FOR_LATER_USE (compute dmat only between degenerate states)
    ! max_num_grps = MAXVAL(num_degen_grps)
    ! ALLOCATE(dmat(max_num_grps, nkstot, nsym), STAT = ierr)
    ! IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error allocating dmat', 1)
    ! DO igrp = 1, max_num_grps
    !   DO ik = 1, nkstot
    !     ndim = degen_grp_ndim(igrp, ik)
    !     !
    !     DO isym = 1, nsym
    !       ALLOCATE(dmat(igrp, ik, isym)%arr(ndim, ndim), STAT = ierr)
    !       IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', &
    !           'Error allocating dmat%arr', 1)
    !     ENDDO
    !     !
    !   ENDDO
    ! ENDDO
    ! END FOR_LATER_USE
    !
    ALLOCATE(dmat_all(nbndep, nbndep, nks, nsym, ntrev))
    !
    DO ik = 1, nks
      !
      WRITE(stdout, '(5x,a,I8)') 'Computing dmat, ik = ', ik
      !
      ! Global k point index
      !
      ik_gl = ik + lower_bnd - 1
      !
      ! Read wavefunction psi(k)
      !
      ipooltmp = my_pool_id + 1
      CALL readwfc(ipooltmp, ik, evc)
      !
      DO itrev = 1, ntrev
        DO isym = 1, nsym
          !
          isk_gl = sym_ktok(ik_gl, isym, itrev)
          !
          ! Find k point and pool index of Sk
          xk_diff(:) = xk_all(:, isk_gl) - xk_loc(:, ik)
          CALL ktokpmq(xk_loc(:, ik), xk_diff, +1, ipool, nsk, isk_gl_)
          !
          if (isk_gl_ /= isk_gl) then
            CALL errore('calc_rotation_gauge', &
                'isk_gl from ktokpmq not equal to isk_gl from sym_ktok', 1)
          endif
          !
          ! Read wavefunction psi(Sk) for isk_gl = sym_ktok(ik_gl, isym)
          !
          CALL readwfc(ipool, nsk, evsk)
          !
          ! ---------------------------------------------------------------------
          ! Compute overlap
          ! Adopted from elphel2_shuffle.f90
          ! Removed parts which are related to the nonlocal projectors
          ! (Will not work for USPPs...)
          ! ---------------------------------------------------------------------
          !
          ! --------------------------------------------------
          !  Calculate SO(3) in cartesian from s_axis_to_cart
          !  then calculate SU(2) transformation matrix from
          !  find_u.
          ! --------------------------------------------------
          aux4(:, :) = czero
          aux5(:, :) = czero
          s_cart(:, :) = zero
          su2_no_dagger(:, :) = czero
          su2(:, :) = czero
          CALL s_crystocart(s(:,:,isym), s_cart, at, bg)
          CALL find_u(s_cart, su2_no_dagger)
          su2 = TRANSPOSE(CONJG(su2_no_dagger))
          !
          ! JML: in contrast to elphel2_shuffle.f90, we use su2 without dagger
          ! JML TODO: WHY?
          su2 = su2_no_dagger
          !
          ! Now we define the igk and igkq from the global igk_k_all
          !
          npw   = ngk_all(ik_gl)
          npwsk = ngk_all(isk_gl)
          !
          ALLOCATE(igk(npw), STAT = ierr)
          IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating igk', 1)
          ALLOCATE(igsk(npwsk), STAT = ierr)
          IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating igsk', 1)
          !
          igk = igk_k_all(1:npw, ik_gl)
          igsk = igk_k_all(1:npwsk, isk_gl)
          !
          !------------------------------------------------------------
          ! Rotate evc --> axu5 by SU2^{dagger}
          ! We transform only psi(k), not psi(Sk).
          !------------------------------------------------------------
          aux4 = evsk
          IF (noncolin) THEN
            IF (itrev == 1) THEN
              DO ibnd = 1, nbnd
                DO ig = 1, npw
                  aux5(ig, ibnd) = su2(1, 1) * evc(ig, ibnd) + su2(1, 2) * evc(ig + npwx, ibnd)
                  aux5(ig + npwx, ibnd) = su2(2, 1) * evc(ig, ibnd) + su2(2, 2) * evc(ig + npwx, ibnd)
                ENDDO
              ENDDO
            ELSE ! itrev == 2: time reversal
              DO ibnd = 1, nbnd
                DO ig = 1, npw
                  aux5(ig + npwx, ibnd) = - ( su2(1, 1) * evc(ig, ibnd) + su2(1, 2) * evc(ig + npwx, ibnd) )
                  aux5(ig, ibnd) = su2(2, 1) * evc(ig, ibnd) + su2(2, 2) * evc(ig + npwx, ibnd)
                ENDDO
              ENDDO
            ENDIF ! itrev
          ELSE
            aux5 = evc
          ENDIF
          !
          ! --------------------------------------------------
          !   Fourier translation of the G-sphere igsk
          ! --------------------------------------------------
          !
          !  Translate by G_0 the G-sphere where evsk is defined,
          !  none of the G-points are lost.
          !
          DO ig0 = 1, ng0vec
            IF ( (ABS(g0vec_sym(1, ik_gl, isym, itrev) - g0vec_all_r(1, ig0)) < eps6) .AND. &
                 (ABS(g0vec_sym(2, ik_gl, isym, itrev) - g0vec_all_r(2, ig0)) < eps6) .AND. &
                 (ABS(g0vec_sym(3, ik_gl, isym, itrev) - g0vec_all_r(3, ig0)) < eps6) ) THEN
              shift_sym = ig0
              EXIT
            ENDIF
            IF (ig0 == ng0vec) THEN
              CALL errore('calc_rotation_gauge', 'cannot find the folding vector in the list', 1)
            ENDIF
          ENDDO
          !
          DO ig = 1, npwsk
            imap = ng0vec * (igsk(ig) - 1) + shift_sym
            igsk_tmp(ig) = gmap(imap)
          ENDDO
          igsk = igsk_tmp
          !
          ! In elphel2_shuffle, both psi(k+q) and psi(k) are rotated, while here
          ! we only rotate only psi(k), not psi(Sk).
          ! So, we apply fractrasl and gmapsym only to aux5 (evc), not aux4 (evsk)
          !
          ! ---------------------------------------------------------------------
          ! phase factor arising from fractional traslations
          ! ---------------------------------------------------------------------
          !
          ! In elphel2_shuffle, we compute {S|v}^{-1}, but here we need {S|v}.
          ! So we use invs(isym) instead of isym.
          !
          ! {S|v}|G+k> = |S(G+k)> * exp(-i S(G+k) v)
          !
          ! gmapsym(:, invs(isym)) computes |S(G+k)>.
          ! fractrasl with invs(isym) computes exp(-i S(G) v) = exp(i G (-S^-1(v))).
          ! So, we need an additional phase factor exp(-i S(k) v) = emiskv.
          ! This is done after calculating the overlap.
          !
          CALL fractrasl(npw, igk, aux5, eigv(:, invs(isym)), cone)
          !
          ! ---------------------------------------------------------------------
          ! wave function rotation to generate matrix elements for the star of q
          ! ---------------------------------------------------------------------
          !
          ! ps. don't use npwx instead of npw, npwsk since the unused elements
          ! may be large and blow up gmapsym (personal experience)
          !
          igk(1:npw) = gmapsym(igk(1:npw), invs(isym))
          !
          ! --------------------------------------------------
          !   Calculation of the matrix element
          ! --------------------------------------------------
          !
          ! This invfft and fwfft calls are the most time-consuming step.
          !
          aux3 = czero
          DO ibnd = 1, nbndep
            aux1 = czero
            aux2 = czero
            CALL invfft_wave(npw, igk, aux5(:, ibnd + ibndstart - 1), aux1)
            !
            IF (itrev == 1) THEN
              aux2 = aux1
            ELSE ! itrev == 2: time revesal
              aux2 = CONJG(aux1)
            ENDIF
            !
            CALL fwfft_wave(npwsk, igsk, aux3(:, ibnd), aux2)
          ENDDO
          !
          ! Calculate overlap between aux4 (psi(Sk)) and aux3 (S * psi(k))
          ! TODO: use ZGEMM
          !
          DO ibnd = 1, nbndep
            DO jbnd = 1, nbndep
              dmat_all(jbnd, ibnd, ik, isym, itrev) &
              = dot_product(aux4(1:npwsk, jbnd + ibndstart - 1), aux3(1:npwsk, ibnd))
              IF (noncolin) THEN
                dmat_all(jbnd, ibnd, ik, isym, itrev) &
                = dmat_all(jbnd, ibnd, ik, isym, itrev) &
                + dot_product(aux4(npwx+1:npwx+npwsk, jbnd + ibndstart - 1), &
                              aux3(npwx+1:npwx+npwsk, ibnd))
              ENDIF
            ENDDO
          ENDDO
          !
          ! Multiply phase factor emiskv for fractional translation.
          !
          dmat_all(:, :, ik, isym, itrev) = dmat_all(:, :, ik, isym, itrev) * emiskv(ik_gl, isym, itrev)
          !
          DEALLOCATE(igk, STAT = ierr)
          IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating igk', 1)
          DEALLOCATE(igsk, STAT = ierr)
          IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating igsk', 1)
          !
        ENDDO ! isym
      ENDDO ! itrev
    ENDDO ! ik
    !
    ! In rotate_epmat, epmatq(iq_first) is slimmed down to the first
    ! ndimwin(ikq), ndimwin(ik) states. Here, we slim down dmat.
    ! This means that states lower than the lower bound of the outer window
    ! are removed.
    !
    ALLOCATE(dmat_temp(nbndep, nbndep), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating dmat_temp', 1)
    !
    DO ik = 1, nks
      !
      ! Find index of first and last band inside the outer window
      !
      ib_win_fst = -1
      ib_win_lst = -1
      !
      DO ibnd = 1, nbndep
        IF (lwin(ibnd, ik)) THEN
          ib_win_fst = ibnd
          EXIT
        ENDIF
      ENDDO
      !
      DO ibnd = nbndep, 1, -1
        IF (lwin(ibnd, ik)) THEN
          ib_win_lst = ibnd
          EXIT
        ENDIF
      ENDDO
      !
      IF (ib_win_fst == -1) CALL errore('calc_rotation_gauge', 'ib_win_fst not found', 1)
      IF (ib_win_lst == -1) CALL errore('calc_rotation_gauge', 'ib_win_lst not found', 1)
      !
      DO itrev = 1, ntrev
        DO isym = 1, nsym
          !
          ! Set dmat_temp to identity
          !
          dmat_temp = czero
          DO ibnd = 1, nbndep
            dmat_temp(ibnd, ibnd) = cone
          ENDDO
          !
          DO jbnd = ib_win_fst, ib_win_lst
            DO ibnd = ib_win_fst, ib_win_lst
              dmat_temp(ibnd - ib_win_fst + 1, jbnd - ib_win_fst + 1) = &
              dmat_all(ibnd, jbnd, ik, isym, itrev)
            ENDDO
          ENDDO
          !
          dmat_all(:, :, ik, isym, itrev) = dmat_temp
          !
        ENDDO ! isym
      ENDDO ! itrev
    ENDDO ! ik
    !
    ! Write dmat to outdir/prefix.dmat file.
    !
    lrdmat = 2 * nbndep * nbndep * nks
    CALL diropn(iudmat, 'dmat', lrdmat, exst)
    !
    DO itrev = 1, ntrev
      DO isym = 1, nsym
        irec = isym + (itrev - 1) * nsym
        CALL davcio(dmat_all(:, :, :, isym, itrev), lrdmat, iudmat, irec, +1)
      ENDDO
    ENDDO
    CLOSE(iudmat)
    !
    ! Write sym_ktok to prefix.symk file. sym_ktok contain all k-points.
    !
    IF (meta_ionode) THEN
      OPEN(NEWUNIT=iusymk, FILE = TRIM(prefix) // '.symk', FORM = 'formatted')
      WRITE(iusymk, '(3i10)') nkstot, nsym, ntrev
      !
      DO ik = 1, nkstot
        DO itrev = 1, ntrev
          DO isym = 1, nsym
            WRITE(iusymk, '(4i8)') ik, isym, itrev, sym_ktok(ik, isym, itrev)
          ENDDO
        ENDDO
      ENDDO
      !
      CLOSE(iusymk)
    ENDIF ! meta_ionode
    !
    !  restore original configuration of files
    !
    CALL diropn(iuwfc, 'wfc', lrwfc, exst)
    ! never remove this barrier - > insures that wfcs are restored to each pool before moving on
    CALL mp_barrier(world_comm)
    !
    DEALLOCATE(dmat_temp, STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating dmat_temp', 1)
    DEALLOCATE(evsk, STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating evsk', 1)
    DEALLOCATE(aux1, STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating aux1', 1)
    DEALLOCATE(aux3, STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating aux3', 1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE calc_rotation_gauge
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE write_qmap(iq, iq_irr, iq_first, isym, isym1, s, ft, timerev)
    !--------------------------------------------------------------------------
    !!
    !! Write symmetry information of q-point to prefix.qmap file.
    !! This file is read in coarse Bloch to Wannier transformation of the upper
    !! Fan matrix.
    !!
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : meta_ionode
    USE io_var,        ONLY : iuqmap
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: timerev
    !! True if we are using time reversal
    INTEGER, INTENT(in) :: iq
    !! Current q-point index.
    INTEGER, INTENT(in) :: iq_irr
    !! Current ireducible q-point
    INTEGER, INTENT(in) :: iq_first
    !! Index of current ireducible q-point in the full q-point list.
    INTEGER, INTENT(in) :: isym
    !! The symmetry which generates the current q in the star
    INTEGER, INTENT(in) :: isym1
    !! The symmetry index of isym before shuffling.
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! The symmetry operation for the eigenmodes
    REAL(KIND = DP), INTENT(in) :: ft(3, 48)
    !! the fractional traslations in crystal axis
    !
    INTEGER :: i
    !
    IF (meta_ionode) THEN
      WRITE(iuqmap, '(3i8, i4, L2)') iq, iq_irr, iq_first, isym1, timerev
      DO i = 1, 3
        WRITE(iuqmap, '(3i8)') s(1:3, i, isym)
      ENDDO
      WRITE(iuqmap, '(3f15.8)') ft(1:3, isym)
    ENDIF
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE write_qmap
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE read_qmap(iq, iq_irr, iq_first, isym, isym1, timerev, s, ft, nsym)
    !--------------------------------------------------------------------------
    !!
    !! Read symmetry information of q-point from prefix.qmap file.
    !!
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : world_comm
    USE io_global,     ONLY : meta_ionode, meta_ionode_id
    USE io_files,      ONLY : prefix
    USE io_var,        ONLY : iuqmap
    USE constants_epw, ONLY : eps4
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Current q-point index.
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! The symmetry operation for the eigenmodes
    REAL(KIND = DP), INTENT(in) :: ft(3, 48)
    !! the fractional traslations in crystal axis
    INTEGER, INTENT(in) :: nsym
    !! Number of symmetry operations
    LOGICAL, INTENT(out) :: timerev
    !! True if we are using time reversal
    INTEGER, INTENT(out) :: iq_irr
    !! Current ireducible q-point
    INTEGER, INTENT(out) :: iq_first
    !! Index of current ireducible q-point in the full q-point list.
    INTEGER, INTENT(out) :: isym
    !! The symmetry which generates the current q in the star
    INTEGER, INTENT(out) :: isym1
    !! The symmetry index of isym before shuffling.
    !
    INTEGER :: i
    !! Counter for Cartesian direction
    INTEGER :: iq_
    !! Counter for q-points
    INTEGER :: iq_dummy
    !! Counter for q-points
    INTEGER :: isym_
    !! Counter for symmetry
    INTEGER :: s_isym(3, 3)
    !! Symmetry matrix read from file
    REAL(KIND = DP) :: ft_isym(3)
    !! Fractional translation read from file
    !
    IF (meta_ionode) THEN
      !
      OPEN(iuqmap, FILE = TRIM(prefix) // '.qmap', FORM = 'formatted')
      !
      ! Read iq-th data from file.
      !
      DO iq_dummy = 1, iq
        READ(iuqmap, '(3i8, i4, L2)') iq_, iq_irr, iq_first, isym1, timerev
        DO i = 1, 3
          READ(iuqmap, '(3i8)') s_isym(1:3, i)
        ENDDO
        READ(iuqmap, '(3f15.8)') ft_isym(1:3)
      ENDDO ! iq_dummy
      !
      ! Find isym such that s(:, :, isym) == s_isym and ft(:, isym) = ft_isym
      !
      isym = -1
      !
      DO isym_ = 1, nsym
        IF ( ALL( s(:, :, isym_) == s_isym ) .AND. &
             ALL( ABS(ft(:, isym_) - ft_isym) < eps4 ) ) THEN
          isym = isym_
          EXIT
        ENDIF
      ENDDO
      !
      CLOSE(iuqmap)
      !
    ENDIF ! meta_ionode
    !
    CALL mp_bcast(iq_irr, meta_ionode_id, world_comm)
    CALL mp_bcast(iq_first, meta_ionode_id, world_comm)
    CALL mp_bcast(isym1, meta_ionode_id, world_comm)
    CALL mp_bcast(isym, meta_ionode_id, world_comm)
    CALL mp_bcast(timerev, meta_ionode_id, world_comm)
    !
    IF (isym == -1) CALL errore('read_qmap', 'isym match not found', 1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE read_qmap
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE read_sym_ktok(time_reversal, nsym, nkstot)
    !--------------------------------------------------------------------------
    !!
    !! Read sym_ktok from prefix.symk file.
    !!
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : world_comm
    USE io_global,     ONLY : meta_ionode, meta_ionode_id
    USE io_files,      ONLY : prefix
    USE constants_epw, ONLY : eps4
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: time_reversal
    !! .true. if time-reversal symmetry is present, .false. otherwise.
    INTEGER, INTENT(in) :: nsym
    !! Number of symmetry operations
    INTEGER, INTENT(in) :: nkstot
    !! Total number of k-points
    !
    INTEGER :: iusymk
    !! Unit for reading and writing symk file
    INTEGER :: nkstot_, nsym_, ntrev_
    !! Number read from file
    INTEGER :: ik, itrev, isym, ik_, isym_, itrev_
    !! Counters
    INTEGER :: ierr
    !! Error index
    !
    ntrev = 1
    IF (time_reversal) ntrev = 2
    !
    ALLOCATE(sym_ktok(nkstot, nsym, ntrev), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_sym_ktok', 'Error allocating sym_ktok', 1)
    !
    IF (meta_ionode) THEN
      !
      OPEN(NEWUNIT=iusymk, FILE = TRIM(prefix) // '.symk', FORM = 'formatted')
      !
      READ(iusymk, '(3i10)') nkstot_, nsym_, ntrev_
      IF (nkstot_ /= nkstot) CALL errore('read_sym_ktok', &
          'nkstot value different in file.', 1)
      IF (nsym_ /= nsym) CALL errore('read_sym_ktok', &
          'nsym value different in file.', 1)
      IF (ntrev_ /= ntrev) CALL errore('read_sym_ktok', &
          'ntrev value different in file.', 1)
      !
      DO ik = 1, nkstot
        DO itrev = 1, ntrev
          DO isym = 1, nsym
            READ(iusymk, '(4i8)') ik_, isym_, itrev_, sym_ktok(ik, isym, itrev)
            !
            IF (ik_ /= ik) CALL errore('read_sym_ktok', &
              'ik from file not equal to ik from loop counter', 1)
            IF (itrev_ /= itrev) CALL errore('read_sym_ktok', &
              'itrev from file not equal to itrev from loop counter', 1)
            IF (isym_ /= isym) CALL errore('read_sym_ktok', &
              'isym from file not equal to isym from loop counter', 1)
          ENDDO
        ENDDO
      ENDDO
      !
      CLOSE(iusymk)
      !
    ENDIF ! meta_ionode
    !
    CALL mp_bcast(sym_ktok, meta_ionode_id, world_comm)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE read_sym_ktok
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE unfold_sthmat(iq, nqc, xqc, lwin, lwinq)
    !--------------------------------------------------------------------------
    !!
    !! If iq is an irreducible q-point, read sthmatq from file in read_sthmat.
    !! Else, unfold sthmatq using gauge matrix dmat.
    !!
    !! sthmatq(ik, kmode, lmode, iq)
    !! = CONJG(rotmat(kmode, imode)) * rotmat(lmode, jmode)
    !! * dmat(k) * sthmat_irr(ik, imode, jmode) * dmat(k)^\dagger
    !!
    !! IF timerev == .TRUE.  : sthmat_irr = sthmatq(iq_first)
    !! IF timerev == .FALSE. : sthmat_irr = CONJG(sthmatq(iq_first))
    !!
    !! rotmat((sna, jdir), (na, idir)) = s_cart(jdir, idir) * phase(na)
    !!
    !! phase(na) = exp( i * 2pi * (xq * rtau(:, isym, na) - xq * ft) )
    !! imode = sna + idir
    !! jmode = na + jdir
    !! sna = irt(isym, na)
    !!
    !!--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : czero, cone, twopi
    USE io_files,      ONLY : diropn
    USE pwcom,         ONLY : nks
    USE ions_base,     ONLY : nat, tau
    USE cell_base,     ONLY : at, bg
    USE symm_base,     ONLY : s, nsym, ft, irt
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : nbndep, sthmatq
    USE kfold,         ONLY : ktokpmq
    USE low_lvl,       ONLY : s_crystocart
    USE read_ahc_files, ONLY : read_sthmat
    USE io_var,        ONLY : iudmat
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: lwin(nbndep, nks)
    !! Outer windows at k+q
    LOGICAL, INTENT(in) :: lwinq(nbndep, nks)
    !! Outer windows at k+q
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    INTEGER, INTENT(in) :: nqc
    !! number of qpoints
    REAL(KIND = DP), INTENT(in) :: xqc(3, nqc)
    !! q-point coordinates (cartesian in units of 2piba)
    !
    LOGICAL :: exst
    !! logical variable to check file exists
    LOGICAL :: timerev
    !! True if using time reversal
    INTEGER :: iq_irr
    !! Irreducible q-point index
    INTEGER :: iq_first
    !! Index of irreducible q-point in the full q-point list
    INTEGER :: isym
    !! The symmetry which generates the current q in the star
    INTEGER :: isym1
    !! The symmetry index of isym before shuffling.
    INTEGER :: na
    !! Index of atom
    INTEGER :: sna
    !! Index of rotated atom
    INTEGER :: imode
    !! Mode index
    INTEGER :: jmode
    !! Mode index
    INTEGER :: lmode
    !! Mode index
    INTEGER :: ik
    !! k point index
    INTEGER :: itrev
    !! 1 for symmetry operation without time reversal, 2 with time versal
    INTEGER :: idir
    !! Cartesian direction index
    INTEGER :: ndim
    !! Dimension of the matrix
    INTEGER :: lrdmat
    !! Size of dmat file
    INTEGER :: irec
    !! Record number
    REAL(KIND = DP) :: xq(3)
    !! The q-point vector
    REAL(KIND = DP) :: xq_irr(3)
    !! The irreducible q-point vector
    REAL(DP) :: ft_cart(3)
    !! fractional translation in Cartesian coordinate
    REAL(KIND = DP) :: arg
    !! Temporary variable for computing phase
    REAL(KIND = DP) :: s_cart(3, 3)
    !! Symmetry matrix in cartesian coordinate
    REAL(KIND = DP) :: rtau_loc(3, 48, nat)
    !! coordinates of direct translations.
    COMPLEX(KIND = DP) :: phase(nat)
    !! Phase factor for symmetry operation
    COMPLEX(KIND = DP) :: rotmat(nmodes, nmodes)
    !! Rotation matrix for epmat
    COMPLEX(KIND = DP) :: rotmat_conjg(nmodes, nmodes)
    !! Rotation matrix for epmat
    COMPLEX(KIND = DP) :: sthmat_irr(nbndep, nbndep, nks, nmodes, nmodes)
    !! epmatq at iq_first
    COMPLEX(KIND = DP) :: sthtmp1(nbndep, nbndep)
    !! Temporary variable for epmat
    COMPLEX(KIND = DP) :: sthtmp2(nbndep, nbndep, nks, nmodes, nmodes)
    !! Temporary variable for epmat
    COMPLEX(KIND = DP) :: sthtmp3(nbndep, nbndep, nks, nmodes)
    !! Temporary variable for epmat
    COMPLEX(KIND = DP) :: dmat(nbndep, nbndep, nks)
    !! Temporary variable for epmat
    !
    CALL start_clock('unfold_sthmat')
    !
    CALL read_qmap(iq, iq_irr, iq_first, isym, isym1, timerev, s, ft, nsym)
    xq = xqc(:, iq)
    xq_irr = xqc(:, iq_first)
    !
    itrev = 1
    IF (timerev) itrev = 2
    !
    ! If iq == iq_first, read sthmatq at irreducible q-point from file.
    !
    IF (iq == iq_first) THEN
      !
      CALL read_sthmat(iq_irr, lwin, lwinq)
      !
      ! For irreducible q-point, unfolding is not needed. Save sthmatq in
      ! sthmatq_save for use at next non-irreducible q-points and return.
      !
      iq_first_save = iq_first
      sthmatq_save = sthmatq
      !
      RETURN
      !
    ENDIF
    !
    ! If iq /= iq_first, rotate sthmatq_save to upmatq(iq).
    ! Check sthmatq_save contains the right irreducible q-point.
    IF (iq_first_save /= iq_first) CALL errore('unfold_sthmat', &
        'iq_first_save /= iq_first')
    sthmatq = sthmatq_save
    !
    ! Setup symmetry operation. Compute rtau, s_cart, ft_cart.
    ! We use rtau_loc because rtau is deallocated in elphon_shuffle_wrap
    ! if epbwrite is true.
    !
    CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau_loc, nat)
    !
    CALL s_crystocart(s(:, :, isym), s_cart, at, bg)
    !
    ft_cart(:) = ft(:, isym)
    CALL cryst_to_cart(1, ft_cart, at, +1)
    !
    ! TODO: Open dmat only once
    !
    lrdmat = 2 * nbndep * nbndep * nks
    CALL diropn(iudmat, 'dmat', lrdmat, exst)
    !
    ! Compute phase factor
    !
    DO na = 1, nat
      arg = twopi * ( SUM( xq(1:3) * (rtau_loc(1:3, isym, na) - ft_cart) ) )
      phase(na) = CMPLX(COS(arg), SIN(arg), KIND = DP)
    ENDDO
    !
    ! Compute rotation matrix
    !
    rotmat = czero
    !
    DO imode = 1, nmodes
      !
      ! imode = 3 * (na - 1) + idir
      !
      na = (imode - 1) / 3 + 1
      idir = imode - 3 * (na - 1)
      !
      sna = irt(isym, na)
      !
      rotmat(3 * (sna - 1) + 1:3*sna, imode) = s_cart(1:3, idir) * phase(na)
      !
    ENDDO
    !
    rotmat_conjg = CONJG(rotmat)
    !
    ! Apply gauge matrix to epmatq
    ! sthtmp2 = dmat(ik) * sthmat_irr * dmat(ik)^\dagger
    ! timerev == .TRUE.  : sthmat_irr = CONJG(sthmatq)
    ! timerev == .FALSE. : sthmat_irr = sthmatq
    ! (Here, sthmatq is read from file and no operation is done yet.)
    !
    sthmat_irr = sthmatq
    IF (timerev) sthmat_irr = CONJG(sthmat_irr)
    !
    ! Read dmat for isym and itrev
    !
    irec = isym1 + (itrev - 1) * nsym
    CALL davcio(dmat, lrdmat, iudmat, irec, -1)
    !
    DO jmode = 1, nmodes
      DO imode = 1, nmodes
        DO ik = 1, nks
          !
          CALL ZGEMM('N', 'N', nbndep, nbndep, nbndep, &
            cone, dmat(1, 1, ik), nbndep, &
            sthmat_irr(1, 1, ik, imode, jmode), nbndep, &
            czero, sthtmp1, nbndep)
          !
          CALL ZGEMM('N', 'C', nbndep, nbndep, nbndep, &
            cone, sthtmp1, nbndep, dmat(1, 1, ik), nbndep, &
            czero, sthtmp2(1, 1, ik, imode, jmode), nbndep)
          !
        ENDDO ! ik
      ENDDO ! imode
    ENDDO ! jmode
    !
    ! Apply rotation matrix to sthmatq:
    !
    ! sthmatq(kmode, lmode)
    ! = CONJG(rotmat(kmode, imode)) * rotmat(lmode, jmode) * sthtmp2(imode, jmode)
    !
    sthmatq(:, :, :, :, :) = czero
    !
    DO jmode = 1, nmodes
      !
      ndim = nbndep * nbndep * nks
      !
      ! sthtmp3(:, kmode) = sum_imode sthtmp2(:, imode, jmode) * CONJG(rotmat(kmode, imode))
      !
      CALL ZGEMM('N', 'T', ndim, nmodes, nmodes, &
        cone, sthtmp2(1, 1, 1, 1, jmode), ndim, rotmat_conjg, nmodes,&
        czero, sthtmp3, ndim)
      !
      ! sthmatq(:, 1:kmode, lmode) += rotmat(lmode, jmode) * sthtmp3(:, 1:kmode)
      !
      DO lmode = 1, nmodes
        sthmatq(:, :, :, :, lmode) = sthmatq(:, :, :, :, lmode) &
        + rotmat(lmode, jmode) * sthtmp3(:, :, :, :)
      ENDDO ! lmode
      !
    ENDDO ! jmode
    !
    CLOSE(iudmat)
    !
    CALL stop_clock('unfold_sthmat')
    !
    ! Here, sthmatq(:, :, ik, :, :) belongs to S*k, not k.
    ! We shuffle the k point indices to order them properly.
    !
    CALL start_clock('shuffle_sthmat')
    CALL shuffle_sthmat(isym1, itrev)
    CALL stop_clock('shuffle_sthmat')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE unfold_sthmat
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE shuffle_sthmat(isym1, itrev)
    !--------------------------------------------------------------------------
    !! Here, redistribute the k-point indices of sthmatq using sym_ktok.
    !! Before shuffling, sthmatq(ik) corresponds to Sk, not k.
    !! After shuffling, sthmatq(ik) corresponds to k.
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE mp,            ONLY : mp_sum
    USE mp_pools,      ONLY : inter_pool_comm
    USE constants_epw, ONLY : czero
    USE pwcom,         ONLY : nks, nkstot
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : nbndep, sthmatq
    USE division,      ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: isym1
    !! The symmetry index of isym before shuffling.
    INTEGER, INTENT(in):: itrev
    !! 1 for symmetry operation without time reversal, 2 with time versal
    !
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: ik
    !! k point index
    INTEGER :: ik_gl
    !! k point index
    INTEGER :: isk
    !! k point index
    COMPLEX(KIND = DP) :: sthmat_temp(nbndep, nbndep, nkstot, nmodes, nmodes)
    !! Temporary storage of sthmat
    !
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)
    !
    sthmat_temp = czero
    DO ik = 1, nks
      ik_gl = ik + lower_bnd - 1
      isk = sym_ktok(ik_gl, isym1, itrev)
      sthmat_temp(:, :, isk, :, :) = sthmatq(:, :, ik, :, :)
    ENDDO
    !
    CALL mp_sum(sthmat_temp, inter_pool_comm)
    !
    sthmatq(:, :, :, :, :) = czero
    !
    DO ik = 1, nks
      ik_gl = ik + lower_bnd - 1
      sthmatq(:, :, ik, :, :) = sthmat_temp(:, :, ik_gl, :, :)
    ENDDO
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE shuffle_sthmat
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE rotate_wfn_deallocate()
    !--------------------------------------------------------------------------
    !!
    !! Deallocate variables in this module
    !!
    !--------------------------------------------------------------------------
    INTEGER :: ierr
    !! Error index
    !
    DEALLOCATE(degen_grp_inds, STAT = ierr)
    IF (ierr /= 0) CALL errore('rotate_wfn_deallocate', 'Error deallocating degen_grp_inds', 1)
    DEALLOCATE(degen_grp_ndim, STAT = ierr)
    IF (ierr /= 0) CALL errore('rotate_wfn_deallocate', 'Error deallocating degen_grp_ndim', 1)
    DEALLOCATE(num_degen_grps, STAT = ierr)
    IF (ierr /= 0) CALL errore('rotate_wfn_deallocate', 'Error deallocating num_degen_grps', 1)
    DEALLOCATE(sym_ktok, STAT = ierr)
    IF (ierr /= 0) CALL errore('rotate_wfn_deallocate', 'Error deallocating sym_ktok', 1)
    DEALLOCATE(exband_rotate, STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error deallocating exband_rotate', 1)
    DEALLOCATE(g0vec_sym, STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error deallocating g0vec_sym', 1)
    DEALLOCATE(emiskv, STAT = ierr)
    IF (ierr /= 0) CALL errore('setup_rotate_wavefunction', 'Error deallocating emiskv', 1)
    ! DEALLOCATE(dmat, STAT = ierr)
    ! IF (ierr /= 0) CALL errore('rotate_wfn_deallocate', 'Error deallocating dmat', 1)
    DEALLOCATE(dmat_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('rotate_wfn_deallocate', 'Error deallocating dmat_all', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE rotate_wfn_deallocate
    !-----------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------------
  END MODULE rotate_wavefunction
  !-----------------------------------------------------------------------------

