  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE symmetry
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines to perform symmetry rotations.
  !!
  !----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  !
  IMPLICIT NONE
  !
  TYPE twodarray
    COMPLEX(KIND = DP), ALLOCATABLE :: arr(:, :)
  END TYPE twodarray
  !
  INTEGER :: s_k(3, 3, 96)
  !! the symmetry operations for the k-points, including TR
  !! IMPORTANT NOTE : s_k differs from s in symm_base.
  !!   1. For non-magnetic systems, it doubles the symmetries by including time reversal operations
  !!   2. It includes a -1 factor for operators including time reversal.
  !!   It should be only applied for time-reversal odd vectors (e.g. wavevectors, current).
  !!   See SUBROUTINE kpoints_time_reversal_init.
  INTEGER :: nsym_k
  !! Number of symmetry operations for the k-points, the maximum is 48*TR
  INTEGER :: t_rev_k(96)
  !! A temporary array for kpoint_grid. Should be cleanup in future.
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
  REAL(KIND = DP), ALLOCATABLE :: g0vec_sym(:, :, :, :)
  !! (3, nkstot, nsym, ntrev) G0 in Cartesian representation such that Sk + G0 belongs
  !! to the first BZ.
  COMPLEX(KIND = DP), ALLOCATABLE :: emiskv(:, :, :)
  !! (nkstot, nsym, ntrev) exp(-i S(k) v) for fractional translation.
  COMPLEX(KIND = DP), ALLOCATABLE :: dmat_all(:, :, :, :, :)
  !! (nbnd, nbnd, nks, nsym, ntrev) <psi_S(k)|psi_k>
  ! TYPE(twodarray), ALLOCATABLE :: dmat(:, :, :)
  !! (max_num_grps, nkstot, nsym) block diagonal matrix <psi_S(k)|psi_k>
  !! dmat(igrp, ik, isym)%arr is square matrix with dimension degen_grp_ndim(igrp, ik)
  COMPLEX(KIND = DP), ALLOCATABLE, SAVE :: sthmatq_save(:, :, :, :, :)
  !! sthmatq at irreducible q point.
  !
  CONTAINS
    !
    !--------------------------------------------------------------------------
    SUBROUTINE rotate_eigenm(iq_first, nqc, isym, s, invs, irt, rtau, xq, cz1, cz2, is_mq)
    !--------------------------------------------------------------------------
    !!
    !!  Here:
    !!
    !!  1). determine the unitary matrix which gives the transformed eigenmodes
    !!     from the first q in the star to the present q
    !!
    !!  2). rotate eigenmodes from the first q in the star (cz1) to the current
    !!      q (cz2)
    !!
    !!  In rotate_epmat.f90:
    !!
    !!  3). rotate the electron-phonon matrix from the cartesian representation
    !!     of the first qpoint of the star to the eigenmode representation
    !!     (using cz1).
    !!
    !!  4). rotate the electron-phonon matrix from the eigenmode representation
    !!     to the cartesian representation of the qpoint iq in the star (with cz2).
    !!
    !!     Step 3 requires using the rotated eigenmodes to set the phases
    !!     (the diagonalizers of the rotated dyn mat will not
    !!     work because of the gages in deg. subspaces and the phases).
    !!
    !!  W.r.t. the standard method of q2qstar_ph.f90, here I construct the
    !!  unitary matrix Gamma defined in Maradudin & Vosko, RMP 1968.
    !!  In this way all rotations are just zgemm operations and we throw away
    !!  all nasty indices.
    !!
    !!  SP - Sept 2019: Cleaning of the routine.
    !
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : dynq
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : cone, czero, twopi, rydcm1, eps10, cmm12meV
    USE control_flags, ONLY : iverbosity
    USE cell_base,     ONLY : at, bg
    USE ions_base,     ONLY : nat, amass, nat, ityp
    USE rigid,         ONLY : cdiagh2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq_first
    !! Originating q point of the star
    INTEGER, INTENT(in) :: nqc
    !! Current q-point index
    INTEGER, INTENT(in) :: isym
    !! The symmetry under consideration
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! The symmetry operation for the eigenmodes
    INTEGER, INTENT(in) :: irt(48, nat)
    !! The rotated of each atom
    INTEGER, INTENT(in) :: invs(48)
    !! The index of the inverse operation
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! The rotated q vector
    REAL(KIND = DP), INTENT(in) :: rtau(3, 48, nat)
    !! The relative position of the rotated atom to the original one
    COMPLEX(KIND = DP), INTENT(inout) :: cz1(nmodes, nmodes)
    !! The eigenvectors for the first q in the star
    COMPLEX(KIND = DP), INTENT(inout) :: cz2(nmodes, nmodes)
    !! The rotated eigenvectors, for the current q in the star
    LOGICAL, OPTIONAL, INTENT(in) :: is_mq
    !! Whether we are in the -q list
    !
    ! Local variables
    INTEGER :: info
    !!
    INTEGER :: imode
    !! Mode index
    INTEGER :: jmode
    !! Mode index
    INTEGER :: nu
    !! Mode index
    INTEGER :: mu
    !! Mode index
    INTEGER :: sna
    !!
    INTEGER :: ism1
    !!
    INTEGER :: na
    !! Atom index
    INTEGER :: nb
    !! Atom index
    REAL(KIND = DP) :: arg
    !!
    REAL(KIND = DP) :: massfac
    !!
    REAL(KIND = DP) :: w1(nmodes)
    !! Phonon freq
    REAL(KIND = DP) :: w2(nmodes)
    !! Phonon freq
    REAL(KIND = DP) :: wtmp(nmodes)
    !!
    REAL(KIND = DP) :: scart(3, 3)
    !!
    COMPLEX(KIND = DP) :: gamma(nmodes, nmodes)
    !! The Gamma matrix for the symmetry operation on the dyn mat
    COMPLEX(KIND = DP) :: cfac
    !!
    COMPLEX(KIND = DP) :: dyn1(nmodes, nmodes)
    !! Dynamical matrix
    COMPLEX(KIND = DP) :: dyn2(nmodes, nmodes)
    !! Dynamical matrix
    !
    ! ------------------------------------------------------------------
    ! diagonalize dynq(iq_first) --> w1, cz1
    ! ------------------------------------------------------------------
    !
    ! first divide by the square root of masses
    !
    DO na = 1, nat
      DO nb = 1, nat
        massfac = 1.d0 / DSQRT(amass(ityp(na)) * amass(ityp(nb)))
        dyn1(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
        dynq(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, iq_first) * massfac
      ENDDO
    ENDDO
    !
    DO jmode = 1, nmodes
      DO imode = 1, jmode
        dyn1(imode, jmode) = &
          (dyn1(imode, jmode) + CONJG(dyn1(jmode, imode))) / 2.d0
      ENDDO
    ENDDO
    !
    CALL cdiagh2(nmodes, dyn1, nmodes, w1, cz1)
    !
    IF (iverbosity == 1) THEN
      !
      ! check the frequencies
      DO nu = 1, nmodes
        IF (w1(nu) > 0.d0) THEN
          wtmp(nu) =  DSQRT(ABS(w1(nu)))
        ELSE
          wtmp(nu) = -DSQRT(ABS(w1(nu)))
        ENDIF
      ENDDO
      WRITE(stdout, '(5x,"Frequencies of the matrix for the first q in the star (cm^{-1})")')
      WRITE(stdout, '(6(2x,f10.5))') (wtmp(nu) * rydcm1, nu = 1, nmodes)
      !
    ENDIF
    !
    ! here dyn1 is dynq(:,:,iq_first) after dividing by the masses
    ! (the true dyn matrix D_q)
    !
    ! -----------------------------------------------------------------------
    ! the matrix gamma (Maradudin & Vosko, RMP, eq. 2.37)
    ! -----------------------------------------------------------------------
    !
    ! I have built the matrix by following step-by-step q2star_ph.f90 and
    ! rotate_and_add_dyn.f90
    !
    ism1 = invs(isym)
    !
    ! the symmetry matrix in cartesian coordinates
    ! (so that we avoid going back and forth with the dynmat)
    ! note the presence of both at and bg in the transform!
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
      cfac = dcmplx(COS(arg), -SIN(arg))
      !
      ! the submatrix (sna,na) contains the rotation scart
      !
      gamma(3 * (sna - 1) + 1:3 * sna, 3 * (na - 1) + 1:3 * na) = cfac * scart
      !
    ENDDO
    !
    !  possibly run some consistency checks
    !
    IF (iverbosity == 1) THEN
      !
      !  D_{Sq} = gamma * D_q * gamma^\dagger (Maradudin & Vosko, RMP, eq. 3.5)
      !
      CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma, &
             nmodes, dyn1 , nmodes, czero , dyn2, nmodes)
      CALL ZGEMM('n', 'c', nmodes, nmodes, nmodes, cone, dyn2, &
             nmodes, gamma, nmodes, czero, dyn1, nmodes)
      !
      DO jmode = 1, nmodes
        DO imode = 1, jmode
          dyn1(imode, jmode) = &
            (dyn1(imode, jmode) + CONJG(dyn1(jmode, imode))) / 2.d0
        ENDDO
      ENDDO
      !
      CALL cdiagh2(nmodes, dyn1, nmodes, w2, cz2)
      !
      ! Check the frequencies
      !
      DO nu = 1, nmodes
        IF (w2(nu) > 0.d0) THEN
          wtmp(nu) =  DSQRT(ABS(w2(nu)))
        ELSE
          wtmp(nu) = -DSQRT(ABS(w2(nu)))
        ENDIF
      ENDDO
      WRITE(stdout, '(5x,"Frequencies of the matrix for the first q in the star (meV)")')
      WRITE(stdout, '(6(2x,f10.5))') (wtmp(nu) * rydcm1 * cmm12meV, nu = 1, nmodes)
      !
    ENDIF
    !
    !
    ! -----------------------------------------------------------------------
    ! transformation of the eigenvectors: e_{Sq} = gamma * e_q  (MV eq. 2.36)
    ! -----------------------------------------------------------------------
    !
    ! ZD: fix the eigenvector issues with -q
    IF(PRESENT(is_mq)) THEN
      IF( is_mq .EQV. .TRUE.) THEN
        CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma, &
         nmodes, CONJG(cz1) , nmodes, czero , cz2, nmodes)
        ! cz2 = CONJG(cz2)! We are just imposing e(-q) = e^*(q)
      ELSE
        CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma, &
        nmodes, cz1 , nmodes, czero , cz2, nmodes)
      ENDIF
    ELSE
      CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, gamma, &
      nmodes, cz1 , nmodes, czero , cz2, nmodes)
    ENDIF
    ! ZD
    !
    !   now check that cz2 diagonalizes dyn2 (after dividing by the masses)
    !
    DO na = 1, nat
      DO nb = 1, nat
        massfac = 1.d0 / DSQRT(amass(ityp(na)) * amass(ityp(nb)))
        dyn2(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb) = &
          dynq(3 * (na - 1) + 1:3 * na, 3 * (nb - 1) + 1:3 * nb, nqc) * massfac
      ENDDO
    ENDDO
    !
    ! the rotated matrix and the one read from file
    IF (iverbosity == 1) WRITE(stdout, '(2f15.10)') dyn2 - dyn1
    !
    ! here I have checked that the matrix rotated with gamma
    ! is perfectly equal to the one read from file for this q in the star
    !
    CALL ZGEMM('c', 'n', nmodes, nmodes, nmodes, cone, cz2, nmodes, dyn2, nmodes, czero, dyn1, nmodes)
    CALL ZGEMM('n', 'n', nmodes, nmodes, nmodes, cone, dyn1, nmodes, cz2, nmodes, czero, dyn2, nmodes)
    !
    DO nu = 1, nmodes
      w2(nu) = ABS(dyn2(nu, nu))
      DO mu = 1, nmodes
        IF (mu /= nu .AND. ABS(dyn2(mu, nu)) > eps10) CALL errore('rotate_eigenm', 'problem with rotated eigenmodes', 0)
      ENDDO
    ENDDO
    !
    IF (iverbosity == 1) THEN
      !
      ! A simple check on the frequencies
      !
      DO nu = 1, nmodes
        IF (w2(nu) > 0.d0 ) then
          wtmp(nu) =  DSQRT(ABS(w2(nu)))
        ELSE
          wtmp(nu) = -DSQRT(ABS(w2(nu)))
        ENDIF
      ENDDO
      WRITE(stdout, '(5x,"Frequencies of the matrix for the current q in the star (cm^{-1})")')
      WRITE(stdout, '(6(2x,f10.5))') (wtmp(nu) * rydcm1, nu = 1, nmodes)
      !
    ENDIF
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE rotate_eigenm
    !--------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE rotate_epmat(cz1, cz2, xq, iq, lwin, lwinq, exband, is_mq)
    !---------------------------------------------------------------------------
    !!
    !! 1). rotate the electron-phonon matrix from the cartesian representation
    !!    of the first qpoint of the star to the eigenmode representation
    !!    (using cz1).
    !!
    !! 2). rotate the electron-phonon matrix from the eigenmode representation
    !!     to the cartesian representation of the qpoint iq (with cz2).
    !!
    !! SP - Sep. 2019: Cleaning.
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE input,         ONLY : exciton
    USE global_var,    ONLY : epmatq, nbndep, epf17
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : cone, czero, one, ryd2mev, eps8
    USE pwcom,         ONLY : nbnd, nks
    USE ions_base,     ONLY : amass, ityp
    USE mp_global,     ONLY : my_pool_id
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: lwin(nbndep, nks)
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: lwinq(nbndep, nks)
    !! Bands at k+q within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
    INTEGER, INTENT(in) :: iq
    !!  Current qpoint
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !  Rotated q vector
    COMPLEX(KIND = DP), INTENT(inout) :: cz1(nmodes, nmodes)
    !! eigenvectors for the first q in the star
    COMPLEX(KIND = DP), INTENT(inout) :: cz2(nmodes, nmodes)
    !!  Rotated eigenvectors for the current q in the star
    LOGICAL, OPTIONAL, INTENT(in) :: is_mq
    !! Whether we are in the -q list
    !
    ! Local variables
    INTEGER :: mu
    !! Counter on phonon branches
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: i
    !! Counter on band index
    INTEGER :: j
    !! Counter on band index
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: massfac
    !! square root of mass
    COMPLEX(KIND = DP) :: eptmp(nmodes)
    !! temporary e-p matrix elements
    COMPLEX(KIND = DP) :: cz_tmp(nmodes, nmodes)
    !! temporary variables
    COMPLEX(KIND = DP) :: cz2t(nmodes, nmodes)
    !! temporary variables
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatq_opt(:, :, :, :)
    !! e-p matrix elements in the outer window
    CHARACTER(LEN = 20) :: tp, rank_char, filename ! ZD: test
    !
    ALLOCATE(epmatq_opt(nbndep, nbndep, nks, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('rotate_epmat', 'Error allocating epmatq_opt', 1)
    !
    ! the mass factors:
    !  1/sqrt(M) for the  direct transform
    !  DSQRT(M)   for the inverse transform
    !
    ! if we set cz1 = cz2 here and we calculate below
    ! cz1 * cz2 we find the identity
    !
    cz2t = cz2
    !
    DO mu = 1, nmodes
      na = (mu - 1) / 3 + 1
      massfac = DSQRT(amass(ityp(na)))
      cz1(mu, :) = cz1(mu, :) / massfac
      cz2(mu, :) = cz2(mu, :) * massfac
      cz2t(mu, :) = cz2t(mu, :) / massfac
    ENDDO
    !
    ! the inverse transform also requires the hermitian conjugate
    !
    cz_tmp = CONJG(TRANSPOSE(cz2))
    cz2 = cz_tmp
    !
    ! slim down to the first ndimwin(ikq), ndimwin(ik) states within the outer window
    !
    epmatq_opt = czero
    DO ik = 1, nks
      jbnd = 0
      DO j = 1, nbndep
        IF (lwin(j, ik)) THEN
          jbnd = jbnd + 1
          ibnd = 0
          DO i = 1, nbndep
            IF (lwinq(i, ik)) THEN
              ibnd = ibnd + 1
              epmatq_opt(ibnd, jbnd, ik, :) = epmatq(i, j, ik, :, iq)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !
    ! ep_mode(j) = cfac * sum_i ep_cart(i) * u(i,j)
    !
    epmatq(:, :, :, :, iq) = czero
    DO ik = 1, nks
      DO jbnd = 1, nbndep
        DO ibnd = 1, nbndep
          !
          ! bring e-p matrix from the cartesian representation of the
          ! first q in the star to the corresponding eigenmode representation
          !
          ! ZD: fix t.r.s-related problem
          IF(PRESENT(is_mq)) THEN
            IF(is_mq .EQV. .TRUE.) THEN
              CALL ZGEMV('t', nmodes, nmodes, cone, CONJG(cz1), nmodes,  &
                        epmatq_opt(ibnd, jbnd, ik, :), 1, czero, eptmp, 1)
            ELSE
              CALL ZGEMV('t', nmodes, nmodes, cone, cz1, nmodes,  &
              epmatq_opt(ibnd, jbnd, ik, :), 1, czero, eptmp, 1)
            ENDIF
          ELSE
            CALL ZGEMV('t', nmodes, nmodes, cone, cz1, nmodes,  &
            epmatq_opt(ibnd, jbnd, ik, :), 1, czero, eptmp, 1)
          ENDIF
          !
          ! ZD: Store the coarse-grid e-ph matrix for ex-ph calculations
          IF(exciton) THEN
            epf17(ibnd, jbnd, ik, :) = eptmp(:)
          ENDIF
          !
          ! rotate epmat in the cartesian representation for this q in the star
          !
          CALL ZGEMV('t', nmodes, nmodes, cone, cz2, nmodes, &
                    eptmp, 1, czero, epmatq(ibnd, jbnd, ik, :, iq), 1)
        ENDDO
      ENDDO
    ENDDO
    !
    DEALLOCATE(epmatq_opt, STAT = ierr)
    IF (ierr /= 0) CALL errore('rotate_epmat', 'Error deallocating epmatq_opt', 1)
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE rotate_epmat
    !---------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE star_q2(xq, at, bg, nsym, s, invs, t_rev, nq, sxq, isq, imq, verbosity, sym_smallq)
    !-----------------------------------------------------------------------
    !!
    !! Generate the star of q vectors that are equivalent to the input one
    !! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
    !! i.e. not those of the small-qroup of q only
    !!
    !! SP - Sept. 2019 - Cleaning
    USE io_global,  ONLY : stdout
    USE kinds,      ONLY : DP
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: verbosity
    !! if true prints several messages.
    INTEGER, INTENT(in) :: nsym
    !! nsym matrices of symmetry operations
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! Symmetry operations
    INTEGER, INTENT(in) :: invs(48)
    !! list of inverse operation indices
    INTEGER, INTENT(in) :: t_rev(48)
    !! Time-reveral sym
    INTEGER, INTENT(out) :: nq
    !! degeneracy of the star of q
    INTEGER, INTENT(out) :: isq(48)
    !! index of q in the star for a given sym
    INTEGER, INTENT(out) :: imq
    !! index of -q in the star (0 if not present)
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! q vector
    REAL(KIND = DP), INTENT(in) :: at(3, 3)
    !! direct lattice vectors
    REAL(KIND = DP), INTENT(in) :: bg(3, 3)
    !! reciprocal lattice vectors
    REAL(KIND = DP), INTENT(inout) :: sxq(3, 48)
    !! list of vectors in the star of q
    !
    ! Local variables
    LOGICAL, EXTERNAL :: eqvect
    !! function used to compare two vectors
    INTEGER :: sym_smallq(48)
    !!
    INTEGER :: nsq(48)
    !! number of symmetry ops. of bravais lattice
    INTEGER :: isym
    !! counters on symmetry ops.
    INTEGER :: ism1
    !! index of inverse of isym
    INTEGER :: iq
    !! q-counter
    INTEGER :: i
    !! Counter
    REAL(KIND = DP) :: saq(3, 48)
    !! auxiliary list of q (crystal coordinates)
    REAL(KIND = DP) :: aq(3)
    !! input q in crystal coordinates
    REAL(KIND = DP) :: raq(3)
    !! rotated q in crystal coordinates
    REAL(KIND = DP) :: zero(3)
    !! a zero vector: used in eqvect
    REAL(KIND = DP), PARAMETER :: accep = 1.e-5_dp
    !! Tolerence
    !
    zero(:) = 0.d0
    saq(:, :) = 0.0d0
    !
    ! go to  crystal coordinates
    !
    DO i = 1, 3
      aq(i) = xq(1) * at(1,i) + xq(2) * at(2,i) + xq(3) * at(3,i)
    ENDDO
    !
    ! create the list of rotated q
    !
    DO i = 1, 48
      nsq(i) = 0
      isq(i) = 0
    ENDDO
    nq = 0
    DO isym = 1, nsym
      ism1 = invs(isym)
      DO i = 1, 3
        raq(i) = s(i, 1, ism1) * aq(1) &
               + s(i, 2, ism1) * aq(2) &
               + s(i, 3, ism1) * aq(3)
      ENDDO
      DO i = 1, 3
        sxq(i, 48) = bg(i, 1) * raq(1) &
                   + bg(i, 2) * raq(2) &
                   + bg(i, 3) * raq(3)
      ENDDO
      IF (t_rev(isym) == 1) raq = -raq
      !
      DO iq = 1, nq
        IF (eqvect(raq, saq(1, iq), zero, accep)) THEN
          isq(isym) = iq
          nsq(iq) = nsq(iq) + 1
        ENDIF
      ENDDO
      IF (isq(isym) == 0) THEN
        nq = nq + 1
        nsq(nq) = 1
        isq(isym) = nq
        saq(:, nq) = raq(:)
        sym_smallq(nq) = isym
        DO i = 1, 3
          sxq(i, nq) = bg(i, 1) * saq(1, nq) &
                     + bg(i, 2) * saq(2, nq) &
                     + bg(i, 3) * saq(3, nq)
        ENDDO
      ENDIF
    ENDDO
    !
    ! SP: Now determine which q in the star is obtained
    !
    ! Set imq index if needed and check star degeneracy
    !
    raq(:) = -aq(:)
    imq = 0
    DO iq = 1, nq
      IF (eqvect(raq, saq(1, iq), zero, accep)) imq = iq
      IF (nsq(iq) * nq /= nsym) CALL errore('star_q', 'wrong degeneracy', iq)
    ENDDO
    !
    ! writes star of q
    !
    IF (verbosity) THEN
      WRITE(stdout, *)
      WRITE(stdout, '(5x,a,i4)') 'Number of q in the star = ', nq
      WRITE(stdout, '(5x,a)') 'List of q in the star:'
      WRITE(stdout, '(7x,i4,3f14.9)') (iq, (sxq(i, iq), i = 1, 3), iq = 1, nq)
      IF (imq == 0) THEN
        WRITE(stdout, '(5x,a)') 'In addition there is the -q list: '
        WRITE(stdout, '(7x,i4,3f14.9)') (iq, (-sxq(i, iq), i = 1, 3), iq = 1, nq)
      ENDIF
    ENDIF
    RETURN
    !---------------------------------------------------------------------------
    END SUBROUTINE star_q2
    !---------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gmap_sym(nsym, s, ft, gmapsym, eigv)
    !-----------------------------------------------------------------------
    !!
    !! For every G vector, find S(G) for all the symmetry operations
    !! of the crystal. Construct the matrix
    !! eigv(ig,isym) = $e^{i G v(S)}$ where v(S) is the (possible)
    !! fractional translation associated with the symmetry operation
    !!
    !! No parallelization on G-vecs at the moment
    !! (actually this is done on the global array, but in elphel2.f90
    !! every processor has just a chunk of the array, I may need some
    !! communication)
    !!
    !----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : twopi, ci, cone
    USE gvect,         ONLY : mill
    USE global_var,    ONLY : mapg, ngxxf
    USE fft_base,      ONLY : dffts
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nsym
    !! the number of symmetries of the crystal
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! the symmetry matrices
    INTEGER, INTENT(out) :: gmapsym(ngxxf, nsym)
    !! the map S(G) = gmapsym (G,S) 1...nsym
    REAL(KIND = DP), INTENT(in) :: ft(3, 48)
    !! the fractional traslations in crystal axis
    COMPLEX(KIND = DP), INTENT(out) :: eigv(ngxxf, nsym)
    !! e^{ iGv} for 1...nsym
    !
    ! Local variables
    INTEGER :: ig
    !! Counter on the G-vector
    INTEGER :: i
    !! Index of the rotated G-vector
    INTEGER :: j
    !! Index of the rotated G-vector
    INTEGER :: k
    !! Index of the rotated G-vector
    INTEGER :: isym
    !! Counter on the symmetry
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: rdotk
    !! $$\mathbf{r}\cdot\mathbf{k}
    !
    i = (dffts%nr1 - 1) / 2
    j = (dffts%nr2 - 1) / 2
    k = (dffts%nr3 - 1) / 2
    ALLOCATE(mapg(-i:i, -j:j, -k:k), STAT = ierr)
    IF (ierr /= 0) CALL errore('gmap_sym', 'Error allocating mapg', 1)
    mapg = 0
    DO ig = 1, ngxxf
      mapg(mill(1, ig), mill(2, ig), mill(3, ig)) = ig
    ENDDO
    !
    ! Loop on the symmetries of the crystal
    !
    DO isym = 1, nsym
      !
      ! Loop on the G vectors
      !
      DO ig = 1, ngxxf
        !
        ! The rotated G-vector
        !
        i = s(1, 1, isym) * mill(1, ig) + s(1, 2, isym) * mill(2, ig) + s(1, 3, isym) * mill(3, ig)
        j = s(2, 1, isym) * mill(1, ig) + s(2, 2, isym) * mill(2, ig) + s(2, 3, isym) * mill(3, ig)
        k = s(3, 1, isym) * mill(1, ig) + s(3, 2, isym) * mill(2, ig) + s(3, 3, isym) * mill(3, ig)
        gmapsym(ig, isym) = mapg(i, j, k)
        !
        ! now the phase factors e^{iGv}
        !
        IF (ft(1, isym)**2 + ft(2, isym)**2 + ft(3, isym)**2 > 1.0d-8) THEN
          !
          rdotk = DBLE(mill(1, ig)) * ft(1, isym) &
                + DBLE(mill(2, ig)) * ft(2, isym) &
                + DBLE(mill(3, ig)) * ft(3, isym)
          !
          ! the actual translation is -v (have a look at ruota_ijk.f90)
          !
          eigv(ig, isym) = EXP(-ci * twopi * rdotk)
          !
        ELSE
          eigv(ig, isym) = cone
        ENDIF
      ENDDO ! ig
      !
    ENDDO
    !
    DEALLOCATE(mapg, STAT = ierr)
    IF (ierr /= 0) CALL errore('gmap_sym', 'Error deallocating mapg', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gmap_sym
    !-----------------------------------------------------------------------
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
    USE ep_constants,  ONLY : zero, eps4, eps6, twopi, ci, cone
    USE io_global,     ONLY : ionode_id, meta_ionode
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_bcast
    USE pwcom,         ONLY : nkstot, nbnd
    USE cell_base,     ONLY : at
    USE symm_base,     ONLY : time_reversal
    USE input,         ONLY : filukk, nbndsub
    USE io_var,        ONLY : iunukk
    USE global_var,    ONLY : nbndep
    USE input,         ONLY : et_all, xk_all
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
      DO ibnd = 1, nbndep
        READ(iunukk, *) jbnd
      ENDDO
      !
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
    USE ep_constants,  ONLY : zero, czero, cone, eps6
    USE mp,            ONLY : mp_barrier
    USE mp_global,     ONLY : my_pool_id, inter_pool_comm, world_comm
    USE io_global,     ONLY : stdout, meta_ionode
    USE io_files,      ONLY : diropn, prefix
    USE wavefunctions, ONLY : evc
    USE fft_base,      ONLY : dffts
    USE fft_wave,      ONLY : invfft_wave, fwfft_wave
    USE pwcom,         ONLY : nbnd, nks, nkstot
    USE cell_base,     ONLY : at, bg
    USE wvfct,         ONLY : npwx
    USE noncollin_module, ONLY : npol, noncolin
    USE units_lr,      ONLY : iuwfc, lrwfc
    USE io,            ONLY : readwfc, readgmap
    USE io_var,        ONLY : iudmat, iusymk
    USE global_var,    ONLY : nbndep, gmap, ngk_all, igk_k_all, ibndkept, &
                              ng0vec, g0vec_all_r, ngxxf
    USE parallelism,   ONLY : fkbounds
    USE kfold,         ONLY : ktokpmq, createkmap
    USE input,         ONLY : xk_loc, xk_all
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
    INTEGER :: hbnd
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
    COMPLEX(KIND = DP), ALLOCATABLE :: aux4(:, :)
    !! Rotated psi_m,Sk WF by SU(2)
    COMPLEX(KIND = DP), ALLOCATABLE :: aux5(:, :)
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
    ALLOCATE(aux4(npwx * npol, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating aux4', 1)
    ALLOCATE(aux5(npwx * npol, nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error allocating aux5', 1)
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
          CALL fractrasl(nbnd, npw, igk, aux5, eigv(:, invs(isym)), cone)
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
            jbnd = ibndkept(ibnd)
            CALL invfft_wave(npwx, npw, igk, aux5(:, jbnd), aux1)
            !
            IF (itrev == 1) THEN
              aux2 = aux1
            ELSE ! itrev == 2: time revesal
              aux2 = CONJG(aux1)
            ENDIF
            !
            CALL fwfft_wave(npwx, npwsk, igsk, aux3(:, ibnd), aux2)
          ENDDO
          !
          ! Calculate overlap between aux4 (psi(Sk)) and aux3 (S * psi(k))
          ! TODO: use ZGEMM
          !
          DO ibnd = 1, nbndep
            DO jbnd = 1, nbndep
              hbnd = ibndkept(jbnd)
              dmat_all(jbnd, ibnd, ik, isym, itrev) &
              = DOT_PRODUCT(aux4(1:npwsk, hbnd), aux3(1:npwsk, ibnd))
              IF (noncolin) THEN
                dmat_all(jbnd, ibnd, ik, isym, itrev) &
                = dmat_all(jbnd, ibnd, ik, isym, itrev) &
                + DOT_PRODUCT(aux4(npwx + 1:npwx + npwsk, hbnd), &
                              aux3(npwx + 1:npwx + npwsk, ibnd))
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
      OPEN(UNIT=iusymk, FILE = TRIM(prefix) // '.symk', FORM = 'formatted')
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
    DEALLOCATE(aux2, STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating aux2', 1)
    DEALLOCATE(aux3, STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating aux3', 1)
    DEALLOCATE(aux4, STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating aux4', 1)
    DEALLOCATE(aux5, STAT = ierr)
    IF (ierr /= 0) CALL errore('calc_rotation_gauge', 'Error deallocating aux5', 1)
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
    USE ep_constants,  ONLY : eps4
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
    USE ep_constants,  ONLY : eps4
    USE io_var,        ONLY : iusymk
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
      OPEN(UNIT=iusymk, FILE = TRIM(prefix) // '.symk', FORM = 'formatted')
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
    USE ep_constants,  ONLY : czero, cone, twopi
    USE io_files,      ONLY : diropn
    USE pwcom,         ONLY : nks
    USE ions_base,     ONLY : nat, tau
    USE cell_base,     ONLY : at, bg
    USE symm_base,     ONLY : s, nsym, ft, irt
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : nbndep, sthmatq
    USE kfold,         ONLY : ktokpmq
    USE low_lvl,       ONLY : s_crystocart
    USE io_ahc,        ONLY : read_sthmat
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
    INTEGER :: ierr
    !! Error status
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
    COMPLEX(KIND = DP), ALLOCATABLE :: sthmat_irr(:, :, :, :, :)
    !! epmatq at iq_first
    COMPLEX(KIND = DP), ALLOCATABLE :: sthtmp1(:, :)
    !! Temporary variable for epmat
    COMPLEX(KIND = DP), ALLOCATABLE :: sthtmp2(:, :, :, :, :)
    !! Temporary variable for epmat
    COMPLEX(KIND = DP), ALLOCATABLE :: sthtmp3(:, :, :, :)
    !! Temporary variable for epmat
    COMPLEX(KIND = DP), ALLOCATABLE :: dmat(:, :, :)
    !! Temporary variable for epmat
    !
    CALL start_clock('unfold_sthmat')
    !
    ALLOCATE(sthmat_irr(nbndep, nbndep, nks, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error allocating sthmat_irr', 1)
    ALLOCATE(sthtmp1(nbndep, nbndep), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error allocating sthtmp1', 1)
    ALLOCATE(sthtmp2(nbndep, nbndep, nks, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error allocating sthtmp2', 1)
    ALLOCATE(sthtmp3(nbndep, nbndep, nks, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error allocating sthtmp3', 1)
    ALLOCATE(dmat(nbndep, nbndep, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error allocating dmat', 1)
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
        'iq_first_save /= iq_first', 1)
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
    DEALLOCATE(sthmat_irr, STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error deallocating sthmat_irr', 1)
    DEALLOCATE(sthtmp1, STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error deallocating sthtmp1', 1)
    DEALLOCATE(sthtmp2, STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error deallocating sthtmp2', 1)
    DEALLOCATE(sthtmp3, STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error deallocating sthtmp3', 1)
    DEALLOCATE(dmat, STAT = ierr)
    IF (ierr /= 0) CALL errore('unfold_sthmat', 'Error deallocating dmat', 1)
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
    USE ep_constants,  ONLY : czero
    USE pwcom,         ONLY : nks, nkstot
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : nbndep, sthmatq
    USE parallelism,   ONLY : fkbounds
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
    INTEGER :: ierr
    !! Error status
    COMPLEX(KIND = DP), ALLOCATABLE :: sthmat_temp(:, :, :, :, :)
    !! Temporary storage of sthmat
    !
    ALLOCATE(sthmat_temp(nbndep, nbndep, nkstot, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('shuffle_sthmat', 'Error allocating sthmat_temp', 1)
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
    DEALLOCATE(sthmat_temp, STAT = ierr)
    IF (ierr /= 0) CALL errore('shuffle_sthmat', 'Error deallocating sthmat_temp', 1)
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
    !--------------------------------------------------------------------------
    SUBROUTINE kpoints_time_reversal_init()
    !--------------------------------------------------------------------------
    !!
    !!  If needed, add time_reversal symmetry into k points grid,
    !!  then, the symmetry of k is (Crystal symmetry)*TR
    !!
    !--------------------------------------------------------------------------
    !
    USE input,              ONLY : scattering
    USE symm_base,          ONLY : s, nsym, time_reversal, t_rev
    USE noncollin_module,   ONLY : noncolin, domag
    USE io_global,          ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER :: isym
    !
    LOGICAL :: tr_epw = .FALSE.
    !
    IF (scattering) tr_epw = time_reversal
    !
    t_rev_k = 0
    !
    IF (tr_epw) THEN
      WRITE(stdout,'(5x, a)') "Add time reversal symmetry into k grid."
      nsym_k = nsym * 2
    ELSE
      nsym_k = nsym
    ENDIF
    !
    ! Copy symmetries from s to s_k, multiply -1 if symmetry involves time reversal.
    !
    IF (domag .AND. noncolin) THEN
      ! Noncollinear magnetism
      DO isym = 1, nsym
        IF (t_rev(isym) == 0) THEN
          s_k(:, :, isym) = s(:, :, isym)
        ELSE
          s_k(:, :, isym) = -s(:, :, isym)
        ENDIF
      ENDDO
    ELSE
      DO isym = 1, nsym
        s_k(:, :, isym) = s(:, :, isym)
      ENDDO
      !
      IF (tr_epw) THEN
        ! Rotation + time reversal applied to k vector, so multiply -1 sign.
        DO isym = 1, nsym
          s_k(:, :, isym + nsym) = -s(:, :, isym)
        ENDDO
      ENDIF
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kpoints_time_reversal_init
    !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  END MODULE symmetry
  !-----------------------------------------------------------------------------
