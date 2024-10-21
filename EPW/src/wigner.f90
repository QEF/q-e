  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2019 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE wigner
  !----------------------------------------------------------------------
  !!
  !! This module contains all the routines linked creation of Wigner-Seitz cell
  !!
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: NSEARCH = 3
  !! Number of Born-von Karman supercell lattice vectors to search for the Wigner-Seitz vector.
  !! [-NSEARCH, NSEARCH]^3 = (2 * NSEARCH + 1)^3 points are searched in total.
  REAL(KIND = DP), PARAMETER :: WS_DIST_TOL = 1.0E-6_DP
  !! Tolerance for finding degeneracies in the Wigner-Seitz selection.
  !! Distance differences below WS_DIST_TOL is considered equal.
  !
  CONTAINS
    !-----------------------------------------------------------------
    SUBROUTINE wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, &
                                 irvec_k,  irvec_q,  irvec_g,  &
                                 ndegen_k, ndegen_q, ndegen_g, &
                                 wslen_k,  wslen_q,  wslen_g,  &
                                 w_centers, dims, tau, dims2)
    !-----------------------------------------------------------------
    !!
    !! Aug  2020 - SP : Optimal WS set for g (wigner_seitzgk and wigner_seitzgq)
    !! June 2018 - SP - CV : Optimal WS set for el and ph (wigner_seitzk and wigner_seitzq)
    !!
    !! This routine wrap the call to three Wigner-Seitz routines:
    !!   wigner_seitzk : Creates a set of WS vectors for each pair of Wannier centers r_n - r_m
    !!                   On exiting, ndegen_k contains the degeneracies of each pairs
    !!                   of wannier centers while irvec_k contains the minimal communal sets of WS vectors.
    !!                   Used for electronic properties
    !!   wigner_seitzq : Creates a set of WS vectors for each pair of atoms tau(nb)-tau(na)
    !!                   On exiting, ndegen_q contains the degeneracies of each pairs
    !!                   of atoms while irvec_q contains the minimal communal sets of WS vectors.
    !!                   Used for phonon properties
    !!   wigner_seitzg : Creates a set of WS vector for each atoms tau(na).
    !!                   On exiting, ndegen_g contains the degeneracies of each atoms
    !!                   while irvec_g contains the minimal communal sets of WS vectors.
    !!                   Used for electron-phonon properties.
    !!
    !! Note 1: ndegen_k, ndegen_q and ndegen_g might contains 0 weigths.
    !! Note 2: No sorting of vectors is needed anymore
    !! Note 3: The dimension 20*nkc1*nkc2*nkc3 should be safe enough.
    !! Note 4: The Wigner-Seitz construction in EPW was done by constructing a cell
    !!         centred unit cell. This is fine for electronic properties (this is what is done in wannier90).
    !!         However for phonon or electron-phonon properties, one can have issues when the cell
    !!         is tilted for example.
    !!         The proper way is to construct a set of WS vectors centred on pairs of atoms (phonons)
    !!         or atoms (el-ph).
    !!         In the matdyn code, a FT grid is constructed with weights centred on pairs of atoms
    !!         and zeros everywhere else.
    !!         EPW now reproduced exactly the results of matdyn for the interpolated phonons at a
    !!         lower computation cost. Indeed we minimize the number of zeros by keeping the union
    !!         of values between all the cells.
    !!         In both cases this is very fast anyway but is important for el-ph properties.
    !!
    !-----------------------------------------------------------------
    USE kinds,     ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nkc1
    !! size of the uniform k mesh
    INTEGER, INTENT(in) :: nkc2
    !! size of the uniform k mesh
    INTEGER, INTENT(in) :: nkc3
    !! size of the uniform k mesh
    INTEGER, INTENT(in) :: nqc1
    !! size of the uniform q mesh
    INTEGER, INTENT(in) :: nqc2
    !! size of the uniform q mesh
    INTEGER, INTENT(in) :: nqc3
    !! size of the uniform q mesh
    INTEGER, INTENT(in) :: dims
    !! Number of bands in the Wannier space
    INTEGER, INTENT(in) :: dims2
    !! Number of atoms
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_k(:, :)
    !! Integer components of the ir-th Wigner-Seitz grid point in the basis
    !! of the lattice vectors for electrons
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_q(:, :)
    !! Integer components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_g(:, :)
    !! Integer components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_k(:, :, :)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid that depend on
    !! Wannier centers $R + r_n - r_m$
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_q(:, :, :)
    !! Wigner-Seitz weights for the phonon grid that depend on
    !! atomic positions $R + \tau(nb) - \tau(na)$
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_g(:, :, :)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    !! atomic positions $R + \tau(na) - r_m$
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_k(:)
    !! real-space length for electrons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_q(:)
    !! real-space length for phonons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_g(:)
    !! real-space length for electron-phonons, in units of alat
    REAL(KIND = DP), INTENT(in) :: w_centers(3, dims)
    !! Wannier centers used for the creation of electron and el-ph WS cells
    REAL(KIND = DP), INTENT(in) :: tau(3, dims2)
    !! Atomic position in the unit cell.
    !
    ! Work Variables
    INTEGER :: ir
    !! Index for WS vectors
    INTEGER :: na
    !! Atom index
    INTEGER :: nrr_k
    !! maximum number of WS vectors for the electrons
    INTEGER :: nrr_q
    !! maximum number of WS vectors for the phonons
    INTEGER :: nrr_g
    !! maximum number of WS vectors for the electron-phonon
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: irvec_kk(:, :)
    !! local integer components of the ir-th Wigner-Seitz grid point
    !! in the basis of the lattice vectors for electrons
    INTEGER, ALLOCATABLE :: irvec_qq(:, :)
    !! local integer components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, ALLOCATABLE :: irvec_gg(:, :)
    !! local integer components of the ir-th Wigner-Seitz grid point for electron-phonons
    !! We use nkc1 instead of nqc1 because the k-grid is always larger or equal to q-grid.
    INTEGER, ALLOCATABLE :: ndegen_kk(:, :, :)
    !! local Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, ALLOCATABLE :: ndegen_qq(:, :, :)
    !! local Wigner-Seitz number of degenerescence (weights) for the phonons grid
    INTEGER, ALLOCATABLE :: ndegen_gg(:, :, :)
    !! local Wigner-Seitz number of degenerescence (weights) for the electron-phonons grid
    REAL(KIND = DP), ALLOCATABLE :: wslen_kk(:)
    !! local real-space length for electrons, in units of alat
    REAL(KIND = DP), ALLOCATABLE :: wslen_qq(:)
    !! local real-space length for phonons, in units of alat
    REAL(KIND = DP), ALLOCATABLE :: wslen_gg(:)
    !! local real-space length for electron-phonon, in units of alat
    !
    ALLOCATE(irvec_kk(3, 20 * nkc1 * nkc2 * nkc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating irvec_kk', 1)
    ALLOCATE(irvec_qq(3, 20 * nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating irvec_qq', 1)
    ALLOCATE(irvec_gg(3, 20 * nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating irvec_gg', 1)
    ALLOCATE(ndegen_kk(20 * nkc1 * nkc2 * nkc3, dims, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating ndegen_kk', 1)
    ALLOCATE(ndegen_qq(20 * nqc1 * nqc2 * nqc3, dims2, dims2), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating ndegen_qq', 1)
    ALLOCATE(ndegen_gg(20 * nqc1 * nqc2 * nqc3, dims2, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating ndegen_gg', 1)
    ALLOCATE(wslen_kk(20 * nkc1 * nkc2 * nkc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating wslen_kk', 1)
    ALLOCATE(wslen_qq(20 * nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating wslen_qq', 1)
    ALLOCATE(wslen_gg(20 * nqc1 * nqc2 * nqc3), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating wslen_gg', 1)
    !
    ! Check the bounds
    IF (nqc1 > nkc1 .OR. nqc2 > nkc2 .OR. nqc3 > nkc3 ) call errore &
       ('wigner_seitz_wrap', ' the phonon grid should be smaller than electron grid', 1)
    !
    CALL start_clock('wigner_seitz')
    !
    ! Now generated the un-sorted points for the electrons, phonons and electron-phonon
    !
    ! If dims > 1, it includes the position of Wannier-Centers
    CALL wigner_seitzkq(nkc1, nkc2, nkc3, irvec_kk, ndegen_kk, wslen_kk, nrr_k, w_centers, dims, .TRUE.)
    CALL wigner_seitzkq(nqc1, nqc2, nqc3, irvec_qq, ndegen_qq, wslen_qq, nrr_q, tau, dims2, .TRUE.)
    CALL wigner_seitzg(nqc1, nqc2, nqc3, irvec_gg, ndegen_gg, wslen_gg, nrr_g, w_centers, tau, dims, dims2, .TRUE.)
    !
    ALLOCATE(irvec_k(3, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating irvec_k', 1)
    ALLOCATE(irvec_q(3, nrr_q), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating irvec_q', 1)
    ALLOCATE(irvec_g(3, nrr_g), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating irvec_g', 1)
    ALLOCATE(ndegen_k(nrr_k, dims, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating ndegen_k', 1)
    ALLOCATE(ndegen_q(nrr_q, dims2, dims2), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating ndegen_q', 1)
    ALLOCATE(ndegen_g(dims, nrr_g, dims2), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating ndegen_g', 1)
    ALLOCATE(wslen_k(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating wslen_k', 1)
    ALLOCATE(wslen_q(nrr_q), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating wslen_q', 1)
    ALLOCATE(wslen_g(nrr_g), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error allocating wslen_g', 1)
    !
    ! Create vectors with correct size.
    DO ir = 1, nrr_k
      ndegen_k(ir, :, :) = ndegen_kk(ir, :, :)
      irvec_k(:, ir)     = irvec_kk(:, ir)
      wslen_k(ir)        = wslen_kk(ir)
    ENDDO
    DO ir = 1, nrr_q
      ndegen_q(ir, :, :) = ndegen_qq(ir, :, :)
      irvec_q(:, ir)     = irvec_qq(:, ir)
      wslen_q(ir)        = wslen_qq(ir)
    ENDDO
    DO ir = 1, nrr_g
      DO na = 1, dims2
        ndegen_g(:, ir, na) = ndegen_gg(ir, na, :)
      ENDDO
      irvec_g(:, ir)     = irvec_gg(:, ir)
      wslen_g(ir)        = wslen_gg(ir)
    ENDDO
    !
    DEALLOCATE(irvec_kk, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating irvec_kk', 1)
    DEALLOCATE(irvec_qq, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating irvec_qq', 1)
    DEALLOCATE(irvec_gg, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating irvec_gg', 1)
    DEALLOCATE(ndegen_kk, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating ndegen_kk', 1)
    DEALLOCATE(ndegen_qq, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating ndegen_qq', 1)
    DEALLOCATE(ndegen_gg, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating ndegen_gg', 1)
    DEALLOCATE(wslen_kk, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating wslen_kk', 1)
    DEALLOCATE(wslen_qq, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating wslen_qq', 1)
    DEALLOCATE(wslen_gg, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitz_wrap', 'Error deallocating wslen_gg', 1)
    !
    CALL stop_clock('wigner_seitz')
    !
    !-----------------------------------------------------------------------------
    END SUBROUTINE wigner_seitz_wrap
    !-----------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------------
    SUBROUTINE wigner_seitzkq(nc1, nc2, nc3, irvec, ndegen, wslen, nrr, shift, dims, sort)
    !-----------------------------------------------------------------------------
    !!
    !! Calculates a grid of points that fall inside of (and eventually
    !! on the surface of) the Wigner-Seitz supercell centered on the
    !! origin of the Bravais lattice with primitive translations
    !! nkc1*a_1+nkc2*a_2+nkc3*a_3
    !!
    !-----------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE ep_constants,  ONLY : eps6
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nc1
    !! size of the uniform k mesh
    INTEGER, INTENT(in) :: nc2
    !! size of the uniform k mesh
    INTEGER, INTENT(in) :: nc3
    !! size of the uniform k mesh
    INTEGER, INTENT(in) :: dims
    !! dim of shift(3,dims)
    INTEGER, INTENT(out) :: irvec(3, 20 * nc1 * nc2 * nc3)
    !! INTEGER components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors
    INTEGER, INTENT(out) :: ndegen(20 * nc1 * nc2 * nc3, dims, dims)
    !! Number of degeneracies
    INTEGER, INTENT(out) :: nrr
    !! number of Wigner-Seitz grid points
    REAL(KIND = DP), INTENT(out) :: wslen(20 * nc1 * nc2 * nc3)
    !! real-space length, in units of alat
    REAL(KIND = DP), INTENT(in) :: shift(3, dims)
    !! Wannier centers (electron WS cell) or atomic positions (lattice WS cell)
    LOGICAL, INTENT(in) :: sort  ! DEBUG
    !! If true, sort the WS vectors to reproduce previous behavior. This should not be needed,
    !! but some tests broke so we keep it for now.
    !
    ! Local variables
    LOGICAL :: found
    !! True if the vector has been found
    INTEGER :: n1, n2, n3
    !! Index for the larger supercell
    INTEGER :: i1, i2, i3
    !! Index to compute |r-R| distance
    INTEGER :: i, ir, irtot, iw, iw2
    !! Iterative index
    INTEGER :: ipol, jpol
    !! Cartesian direction
    INTEGER :: cnt
    !! Count for the Wigner-Seitz degeneracy
    INTEGER :: ierr
    !! Error status
    INTEGER :: nrr_tmp(dims, dims)
    !! Temporary WS matrix
    REAL(KIND = DP) :: adot(3, 3)
    !! Dot product between lattice vector
    REAL(KIND = DP) :: mindist
    !! Minimum distance
    REAL(KIND = DP) :: tot, tot2
    !! Sum of all the weigths
    REAL(KIND = DP) :: ndiff(3)
    !! Distances. Must be now real because of fractional atoms
    REAL(KIND = DP), ALLOCATABLE :: dist(:)
    !! Contains the distance squared |r-R|^2
    INTEGER, ALLOCATABLE :: irvec_tmp(:, :, :, :)
    !! Temp
    INTEGER, ALLOCATABLE :: ndegen_tmp(:, :, :)
    !! Number of degenerate points
    INTEGER :: n  ! DEBUG
    !! Temporary variable for sorting
    INTEGER, ALLOCATABLE :: tmp_ind(:)  ! DEBUG
    !! Temporary index array for sorting
    REAL(KIND = DP), ALLOCATABLE :: tmp_ind_for_sort(:)  ! DEBUG
    !! Temporary variable for sorting
    !
    ALLOCATE(dist((2 * NSEARCH + 1)**3), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error allocating dist', 1)
    ALLOCATE(irvec_tmp(3, 20 * nc1 * nc2 * nc3, dims, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error allocating irvec_tmp', 1)
    ALLOCATE(ndegen_tmp(20 * nc1 * nc2 * nc3, dims, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error allocating ndegen_tmp', 1)
    !
    DO ipol = 1, 3
     DO jpol = 1, 3
       adot(ipol, jpol) = DOT_PRODUCT(at(:, ipol), at(:, jpol))
     ENDDO
    ENDDO
    !
    CALL cryst_to_cart(dims, shift(:, :), bg, -1)
    !
    nrr_tmp(:, :) = 0
    ndegen_tmp(:, :, :) = 0
    irvec_tmp(:, :, :, :) = 0
    !
    DO iw = 1, dims
      DO iw2 = 1, dims
        !
        DO n1 = 0, nc1 - 1
          DO n2 = 0, nc2 - 1
            DO n3 = 0, nc3 - 1
              !
              ! For each iw, iw2 (Wannier functions or atoms), and R = (n1, n2, n3) (lattice
              ! vector in the [0, nc)^3 Born-von Karman supercell), find the mapping of
              ! R + shift(:, iw2) - shift(:, iw) to the Wigner-Seitz cell.
              !
              ! We iterate over the supercell lattice vectors T = (i1 * nc1, i2 * nc2, i3 * nc3)
              ! and find T's such that T + R + shift(:, iw2) - shift(:, iw) is inside the Wigner-Seitz cell.
              ! To do so, we compute dist(T) = |T + R + shift(:, iw2) - shift(:, iw)| for a set of T vectors
              ! in [-NSEARCH, NSEARCH]^3. T vectors with the smallest dist(T) are selected. If multiple T's
              ! give the same dist(T) (up to tolerance eps6) they are all selected.
              !
              ! nrr           : The number of unique T + R vectors.
              ! irvec(3, nrr) : The set of all selected T + R vectors.
              ! ndegen(nrr)   : The number of T vectors with the same minimal dist(T).
              !
              ! Compute dist(T) for all T in [-NSEARCH, NSEARCH]^3
              !
              i = 0
              dist(:) = 0.d0
              DO i1 = -NSEARCH, NSEARCH
                DO i2 = -NSEARCH, NSEARCH
                  DO i3 = -NSEARCH, NSEARCH
                    i = i + 1
                    !
                    ! Calculate dist(T) = |R + T + shift(iw2) - shift(iw)|^2
                    ! R = (n1, n2, n3) * at
                    ! T = (i1 * nc1, i2 * nc2, i3 * nc3) * at
                    !
                    ndiff(1) = n1 + i1 * nc1 + shift(1, iw2) - shift(1, iw)
                    ndiff(2) = n2 + i2 * nc2 + shift(2, iw2) - shift(2, iw)
                    ndiff(3) = n3 + i3 * nc3 + shift(3, iw2) - shift(3, iw)
                    DO ipol = 1, 3
                      DO jpol = 1, 3
                        dist(i) = dist(i) + ndiff(ipol) * adot(ipol, jpol) * ndiff(jpol)
                      ENDDO
                    ENDDO
                    !
                  ENDDO ! i3
                ENDDO ! i2
              ENDDO ! i1
              !IF(iw==iw2) print*,'iw dist ',iw,dist(63)
              !
              mindist = MINVAL(dist)
              !
              ! Count the number of vectors T with the (same) smallest dist(T).
              ! This is the degeneracy.
              !
              cnt = 0
              DO i = 1, SIZE(dist, 1)
                IF (dist(i) < mindist + WS_DIST_TOL) cnt = cnt + 1
              ENDDO
              !
              ! Find T's with the (same) smallest dist(T) and store R + T to irvec_tmp.
              !
              i = 0
              DO i1 = -NSEARCH, NSEARCH
                DO i2 = -NSEARCH, NSEARCH
                  DO i3 = -NSEARCH, NSEARCH
                    i = i + 1
                    !
                    IF (dist(i) < mindist + WS_DIST_TOL) THEN
                      !
                      nrr_tmp(iw, iw2) = nrr_tmp(iw, iw2) + 1
                      ndegen_tmp(nrr_tmp(iw, iw2), iw, iw2) = cnt
                      irvec_tmp(:, nrr_tmp(iw, iw2), iw, iw2) = (/ n1 + i1 * nc1, n2 + i2 * nc2, n3 + i3 * nc3 /)
                      !
                    ENDIF
                    !
                  ENDDO ! i3
                ENDDO ! i2
              ENDDO ! i1
              !
            ENDDO ! n3
          ENDDO ! n2
        ENDDO ! n1
      ENDDO ! iw2
    ENDDO ! iw
    !
    ! DEBUG: Sort irvec_tmp to exactly reproduce the previous behavior
    ! (The order should not matter in principle, but seems to affect transport routines.)
    !
    IF (sort) THEN
      !
      n = MAXVAL(nrr_tmp)
      ALLOCATE(tmp_ind(n), STAT = ierr)
      IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error allocating tmp_ind', 1)
      ALLOCATE(tmp_ind_for_sort(n), STAT = ierr)
      IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error allocating tmp_ind_for_sort', 1)
      !
      DO iw = 1, dims
        DO iw2 = 1, dims
          n = nrr_tmp(iw, iw2)
          !
          ndegen(1:n, 1, 1) = ndegen_tmp(1:n, iw, iw2)
          irvec(:, 1:n) = irvec_tmp(:, 1:n, iw, iw2)
          !
          tmp_ind_for_sort(1:n) = irvec(1, 1:n) * 1.d6 + irvec(2, 1:n) * 1.d3 + irvec(3, 1:n) * 1.d0
          tmp_ind(1:n) = (/ (i, i = 1, n) /)
          CALL hpsort(n, tmp_ind_for_sort(1:n), tmp_ind(1:n))
          !
          DO ir = 1, n
            ndegen_tmp(ir, iw, iw2) = ndegen(tmp_ind(ir), 1, 1)
            irvec_tmp(:, ir, iw, iw2) = irvec(:, tmp_ind(ir))
          ENDDO
        ENDDO
      ENDDO
      !
      DEALLOCATE(tmp_ind, STAT = ierr)
      IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error deallocating tmp_ind', 1)
      DEALLOCATE(tmp_ind_for_sort, STAT = ierr)
      IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error deallocating tmp_ind_for_sort', 1)
      !
    ENDIF
    !
    ! Make a unique common list of Wigner-Seitz R vectors for all iw and iw2.
    !
    nrr = nrr_tmp(1, 1)
    irvec(:, :) = irvec_tmp(:, :, 1, 1)
    DO iw = 1, dims
      DO iw2 = 1, dims
        DO ir = 1, nrr_tmp(iw, iw2)
          found = .FALSE.
          DO irtot = 1, nrr
            IF (ALL(irvec_tmp(:, ir, iw, iw2) == irvec(:, irtot))) THEN
              found = .TRUE.
            ENDIF
          ENDDO !nrr
          IF (.NOT.  found) THEN
            nrr = nrr + 1
            irvec(:, nrr) = irvec_tmp(:, ir, iw, iw2)
          ENDIF
        ENDDO ! ir
      ENDDO ! nb
    ENDDO ! na
    !
    ! Creates a pair of atoms depended degeneracy array but with a number of WS
    ! vectors per pair that is equal to the global set. Populate with zero weights
    ! the one that are not part of that pair set.
    ndegen(:, :, :) = 0
    DO iw = 1, dims
      DO iw2 = 1, dims
        DO ir = 1, nrr_tmp(iw, iw2)
          DO irtot = 1, nrr
            IF (ALL(irvec(:, irtot) == irvec_tmp(:, ir, iw, iw2))) THEN
              ndegen(irtot, iw, iw2) = ndegen_tmp(ir, iw, iw2)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    DO iw = 1, dims
      DO iw2 = 1, dims
        tot = 0.d0
        tot2 = 0.d0
        DO i = 1, nrr
          IF (ndegen(i, iw, iw2) > 0) THEN
            tot2 = tot2 + 1.d0 / DBLE(ndegen(i, iw, iw2))
          ENDIF
        ENDDO
        DO i = 1, nrr_tmp(iw, iw2)
          tot = tot + 1.d0 / DBLE(ndegen_tmp(i, iw, iw2))
        ENDDO
        !
        IF (ABS(tot - DBLE(nc1*nc2*nc3)) > eps6) CALL errore &
         ('wigner_seitzkq',' weights do not add up to nc1*nc2*nc3', 1)
        IF (ABS(tot - tot2) > eps6) CALL errore &
         ('wigner_seitzkq', ' weigths of pair of atoms is not equal to global weights', 1)
      ENDDO
    ENDDO
    !
    IF (nrr > 20 * nc1 * nc2 * nc3) CALL errore &
      ('wigner_seitzkq', 'too many WS points, try to increase the bound 20*nc1*nc2*nc3', 1)
    !
    wslen(:) = 0.d0
    DO i = 1, nrr
      DO ipol = 1, 3
        DO jpol = 1, 3
          wslen(i) = wslen(i) + DBLE(irvec(ipol, i)) * adot(ipol, jpol) * DBLE(irvec(jpol, i))
        ENDDO
      ENDDO
      wslen(i) = DSQRT(wslen(i))
    ENDDO
    !
    CALL cryst_to_cart(dims, shift(:, :), at, 1)
    !
    DEALLOCATE(dist, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error deallocating dist', 1)
    DEALLOCATE(irvec_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error deallocating irvec_tmp', 1)
    DEALLOCATE(ndegen_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error deallocating ndegen_tmp', 1)
    !
    !-----------------------------------------------------------------------------
    END SUBROUTINE wigner_seitzkq
    !-----------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------
    SUBROUTINE wigner_seitzg(nc1, nc2, nc3, irvec, ndegen, wslen, nrr, w_centers, tau, dims, dims2, sort)
    !-----------------------------------------------------------------
    !!
    !! SP - Aug 2020
    !! Calculates a grid of points that fall inside of (and eventually
    !! on the surface of) the Wigner-Seitz supercell centered on each atom - wannier center pairs.
    !! Follows Eq. 66 of PRB 55, 10355 (1997).
    !! We are part of the WS if $R_b + \tau_{\kappa} - r_n$ is inside the supercell.
    !!
    !-----------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE ep_constants,  ONLY : eps6
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nc1
    !! size of the uniform q mesh
    INTEGER, INTENT(in) :: nc2
    !! size of the uniform q mesh
    INTEGER, INTENT(in) :: nc3
    !! size of the uniform q mesh
    INTEGER, INTENT(in) :: dims
    !! dims is either nbndsub or 1 depending on use_ws
    INTEGER, INTENT(in) :: dims2
    !! dims2 is either the number of atoms or 1 depending on use_ws
    INTEGER, INTENT(out) :: irvec(3, 20 * nc1 * nc2 * nc3)
    !! Integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors
    INTEGER, INTENT(out) :: ndegen(20 * nc1 * nc2 * nc3, dims2, dims)
    !! Number of degeneracies
    INTEGER, INTENT(out) :: nrr
    !! Number of Wigner-Seitz grid points
    REAL(KIND = DP), INTENT(in) :: w_centers(3, dims)
    !! Wannier centers
    REAL(KIND = DP), INTENT(in) :: tau(3, dims2)
    !! Atomic positions
    REAL(KIND = DP), INTENT(out) :: wslen(20 * nc1 * nc2 * nc3)
    !! real-space length, in units of alat
    LOGICAL, INTENT(in) :: sort  ! DEBUG
    !! If true, sort the WS vectors to reproduce previous behavior. This should not be needed,
    !! but some tests broke so we keep it for now.
    !
    ! Local variables
    LOGICAL :: found
    !! True if the vector has been foun
    INTEGER :: n1, n2, n3
    !! Index for the larger supercell
    INTEGER :: i1, i2, i3
    !! Index to compute |r-R| distance
    INTEGER :: na
    !! Atom index
    INTEGER :: iw
    !! Wannier function index
    INTEGER :: i
    !! Iterative index
    INTEGER :: ir
    !! Iterative index on the pair of atoms
    INTEGER :: irtot
    !! Iterative index on the combined pair of atoms
    INTEGER :: ipol, jpol
    !! Cartesian direction
    INTEGER :: cnt
    !! Count for the Wigner-Seitz degeneracy
    INTEGER :: ierr
    !! Error status
    INTEGER :: nrr_tmp(dims2, dims)
    !! Temporary array that contains the max number of WS vectors
    !! for a pair of atoms.
    REAL(KIND = DP) :: adot(3, 3)
    !! Dot product between lattice vector
    REAL(KIND = DP) :: mindist
    !! Minimum distance
    REAL(KIND = DP) :: tot, tot2
    !! Sum of all the weigths
    REAL(KIND = DP) :: ndiff(3)
    !! Distances. Must be now real because of fractional atoms
    REAL(KIND = DP), ALLOCATABLE :: dist(:)
    !! Contains the distance squared |r-R|^2
    INTEGER, ALLOCATABLE :: irvec_tmp(:, :, :, :)
    !! Temporary WS vectors for each atoms
    INTEGER, ALLOCATABLE :: ndegen_tmp(:, :, :)
    !! Temporary WS vectors weigths for each atoms
    INTEGER :: n  ! DEBUG
    !! Temporary variable for sorting
    INTEGER, ALLOCATABLE :: tmp_ind(:)  ! DEBUG
    !! Temporary index array for sorting
    REAL(KIND = DP), ALLOCATABLE :: tmp_ind_for_sort(:)  ! DEBUG
    !! Temporary variable for sorting
    !
    ALLOCATE(dist((2 * NSEARCH + 1)**3), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzg', 'Error allocating dist', 1)
    ALLOCATE(irvec_tmp(3, 20 * nc1 * nc2 * nc3, dims2, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzg', 'Error allocating irvec_tmp', 1)
    ALLOCATE(ndegen_tmp(20 * nc1 * nc2 * nc3, dims2, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzg', 'Error allocating ndegen_tmp', 1)
    !
    DO ipol = 1, 3
      DO jpol = 1, 3
        adot(ipol, jpol) = DOT_PRODUCT(at(:, ipol), at(:, jpol))
      ENDDO
    ENDDO
    !
    CALL cryst_to_cart(dims2, tau(:, :), bg, -1)
    CALL cryst_to_cart(dims, w_centers(:, :), bg, -1)
    !
    ! Loop over grid points r on a unit cell that is 8 times larger than a
    ! primitive supercell. In the end nrr contains the total number of grids
    ! points that have been found in the Wigner-Seitz cell
    !
    nrr_tmp(:, :) = 0
    irvec_tmp(:, :, :, :) = 0
    ndegen_tmp(:, :, :) = 0
    !
    DO iw = 1, dims
      DO na = 1, dims2
        DO n1 = 0, nc1 - 1
          DO n2 = 0, nc2 - 1
            DO n3 = 0, nc3 - 1
              !
              ! For R = (n1, n2, n3), compute dist(T) = |T + R + tau(:, iw2) - w_centers(:, iw)|
              ! for all T = (i1 * nc1, i2 * nc2, i3 * nc3) and find T's with the smallest dist(T).
              ! See subroutine wigner_seitzkq for further details.
              !
              i = 0
              dist(:) = 0.d0
              DO i1 = -NSEARCH, NSEARCH
                DO i2 = -NSEARCH, NSEARCH
                  DO i3 = -NSEARCH, NSEARCH
                    i = i + 1
                    !
                    ! Calculate dist(T)
                    ndiff(1) = n1 + i1 * nc1 + tau(1, na) - w_centers(1, iw)
                    ndiff(2) = n2 + i2 * nc2 + tau(2, na) - w_centers(2, iw)
                    ndiff(3) = n3 + i3 * nc3 + tau(3, na) - w_centers(3, iw)
                    DO ipol = 1, 3
                      DO jpol = 1, 3
                        dist(i) = dist(i) + DBLE(ndiff(ipol)) * adot(ipol, jpol) * DBLE(ndiff(jpol))
                      ENDDO
                    ENDDO
                    !
                  ENDDO
                ENDDO
              ENDDO
              !
              mindist = MINVAL(dist)
              !
              ! Count the number of vectors T with the (same) smallest dist(T).
              ! This is the degeneracy.
              !
              cnt = 0
              DO i = 1, SIZE(dist, 1)
                IF (dist(i) < mindist + WS_DIST_TOL) cnt = cnt + 1
              ENDDO
              !
              ! Find T's with the (same) smallest dist(T) and store R + T to irvec_tmp.
              !
              i = 0
              DO i1 = -NSEARCH, NSEARCH
                DO i2 = -NSEARCH, NSEARCH
                  DO i3 = -NSEARCH, NSEARCH
                    i = i + 1
                    !
                    IF (dist(i) < mindist + WS_DIST_TOL) THEN
                      !
                      nrr_tmp(na, iw) = nrr_tmp(na, iw) + 1
                      ndegen_tmp(nrr_tmp(na, iw), na, iw) = cnt
                      irvec_tmp(:, nrr_tmp(na, iw), na, iw) = (/ n1 + i1 * nc1, n2 + i2 * nc2, n3 + i3 * nc3 /)
                      !
                    ENDIF
                    !
                  ENDDO ! i3
                ENDDO ! i2
              ENDDO ! i1
              !
            ENDDO ! n3
          ENDDO ! n2
        ENDDO ! n3
      ENDDO ! na
    ENDDO ! iw
    !
    ! DEBUG: Sort irvec_tmp to exactly reproduce the previous behavior
    ! (The order should not matter in principle, but seems to affect transport routines.)
    !
    IF (sort) THEN
      !
      n = MAXVAL(nrr_tmp)
      ALLOCATE(tmp_ind(n), STAT = ierr)
      IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error allocating tmp_ind', 1)
      ALLOCATE(tmp_ind_for_sort(n), STAT = ierr)
      IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error allocating tmp_ind_for_sort', 1)
      !
      DO iw = 1, dims
        DO na = 1, dims2
          n = nrr_tmp(na, iw)
          !
          ndegen(1:n, 1, 1) = ndegen_tmp(1:n, na, iw)
          irvec(:, 1:n) = irvec_tmp(:, 1:n, na, iw)
          !
          tmp_ind_for_sort(1:n) = irvec(1, 1:n) * 1.d6 + irvec(2, 1:n) * 1.d3 + irvec(3, 1:n) * 1.d0
          tmp_ind(1:n) = (/ (i, i = 1, n) /)
          CALL hpsort(n, tmp_ind_for_sort(1:n), tmp_ind(1:n))
          !
          DO ir = 1, n
            ndegen_tmp(ir, na, iw) = ndegen(tmp_ind(ir), 1, 1)
            irvec_tmp(:, ir, na, iw) = irvec(:, tmp_ind(ir))
          ENDDO
        ENDDO
      ENDDO
      !
      DEALLOCATE(tmp_ind, STAT = ierr)
      IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error deallocating tmp_ind', 1)
      DEALLOCATE(tmp_ind_for_sort, STAT = ierr)
      IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error deallocating tmp_ind_for_sort', 1)
      !
    ENDIF
    !
    ! Now creates a global set of WS vectors from all the atoms pair.
    ! Also remove the duplicated ones.
    nrr = nrr_tmp(1, 1)
    irvec(:, :) = irvec_tmp(:, :, 1, 1)
    DO iw = 1, dims
      DO na = 1, dims2
        DO ir = 1, nrr_tmp(na, iw)
          found = .FALSE.
          DO irtot = 1, nrr
            IF (ALL(irvec_tmp(:, ir, na, iw) == irvec(:, irtot))) THEN
              found = .TRUE.
            ENDIF
          ENDDO !nrr
          IF (.NOT.  found) THEN
            nrr = nrr + 1
            irvec(:, nrr) = irvec_tmp(:, ir, na, iw)
          ENDIF
        ENDDO ! ir
      ENDDO ! na
    ENDDO ! iw
    !
    ! Creates a pair of atoms-dependent degeneracy array but with a number of WS
    ! vectors per pair that is equal to the global set. Populate with zero weights
    ! the one that are not part of that pair set.
    ndegen(:, :, :) = 0
    DO iw = 1, dims
      DO na = 1, dims2
        DO ir = 1, nrr_tmp(na, iw)
          DO irtot = 1, nrr
            IF (ALL(irvec(:, irtot) == irvec_tmp(:, ir, na, iw))) THEN
              ndegen(irtot, na, iw) = ndegen_tmp(ir, na, iw)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    DO iw = 1, dims
      DO na = 1, dims2
        tot = 0.d0
        tot2 = 0.d0
        DO i = 1, nrr
          IF (ndegen(i, na, iw) > 0) THEN
            tot2 = tot2 + 1.d0 / DBLE(ndegen(i, na, iw))
          ENDIF
        ENDDO
        DO i = 1, nrr_tmp(na, iw)
          tot = tot + 1.d0 / DBLE(ndegen_tmp(i, na, iw))
        ENDDO
        !
        IF (ABS(tot - DBLE(nc1 * nc2 * nc3)) > eps6) CALL errore &
           ('wigner_seitzg', ' weights do not add up to nqc1*nqc2*nqc3', 1)
        IF (ABS(tot - tot2) > eps6) CALL errore &
           ('wigner_seitzg', ' weigths of pair of atoms is not equal to global weights', 1)
      ENDDO
    ENDDO
    !
    IF (nrr > 20 * nc1 * nc2 * nc3) CALL errore('wigner_seitzg', 'too many WS points, try to increase the bound 20*nc1*nc2*nc3', 1)
    !
    ! Now we compute the WS distances
    wslen(:) = 0.d0
    DO i = 1, nrr
      DO ipol = 1, 3
        DO jpol = 1, 3
          wslen(i) = wslen(i) + DBLE(irvec(ipol, i)) * adot(ipol, jpol) * DBLE(irvec(jpol, i))
        ENDDO
      ENDDO
      wslen(i) = DSQRT(wslen(i))
    ENDDO
    !
    CALL cryst_to_cart(dims2, tau(:, :), at, 1)
    CALL cryst_to_cart(dims, w_centers(:, :), at, 1)
    !
    DEALLOCATE(dist, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzg', 'Error deallocating dist', 1)
    DEALLOCATE(irvec_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzg', 'Error deallocating irvec_tmp', 1)
    DEALLOCATE(ndegen_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzg', 'Error deallocating ndegen_tmp', 1)
    !
    !------------------------------------------------------------------------------------------
    END SUBROUTINE wigner_seitzg
    !------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE wigner_divide_ndegen(array, npol, nbnd, nrr, nmodes, ndegen, dims)
    !---------------------------------------------------------------------------
    !!
    !! Divide an array in the real-space Wannier representation by the degeneracy
    !! factor ndegen. This subroutine must be called once and only once before
    !! calling Wannier-to-Bloch routines.
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : czero
    USE input,         ONLY : use_ws
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: npol
    !! number of polarizations (leading dimension of array)
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nrr
    !! Number of Wigner-Size points
    INTEGER, INTENT(in) :: nmodes
    !! number of phonon modes
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier functions if use_ws == .TRUE. Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: ndegen(nrr, dims, dims)
    !! Degeneracy factor for each R
    COMPLEX(KIND = DP), INTENT(inout) :: array(npol, nbnd, nbnd, nrr, nmodes)
    !! Matrix in Wannier representation
    !
    ! Local variables
    INTEGER :: iw
    !! Band index
    INTEGER :: iw2
    !! Band index
    INTEGER :: ir
    !! WS vectors for electrons.
    INTEGER :: imode
    !! Counter on phonon modes
    !
    IF (use_ws) THEN
      !
      DO imode = 1, nmodes
        DO iw2 = 1, dims
          DO iw = 1, dims
            DO ir = 1, nrr
              !
              IF (ndegen(ir, iw, iw2) > 0) THEN
                array(:, iw, iw2, ir, imode) &
                = array(:, iw, iw2, ir, imode) / REAL(ndegen(ir, iw, iw2), KIND = DP)
              ELSE ! ndegen == 0
                array(:, iw, iw2, ir, imode) = czero
              ENDIF
              !
            ENDDO ! ir
          ENDDO ! iw
        ENDDO ! iw2
      ENDDO ! imode
      !
    ELSE ! .NOT. use_ws
      !
      DO imode = 1, nmodes
        DO ir = 1, nrr
          array(:, :, :, ir, imode) &
          = array(:, :, :, ir, imode) / REAL(ndegen(ir, 1, 1), KIND = DP)
        ENDDO ! ir
      ENDDO ! imode
      !
    ENDIF ! use_ws
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE wigner_divide_ndegen
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE wigner_divide_ndegen_epmat(array, nbnd, nrr_k, nrr_g, nmodes, ndegen, dims, dims2)
    !---------------------------------------------------------------------------
    !!
    !! Divide an array in the real-space Wannier representation by the degeneracy
    !! factor ndegen. This subroutine must be called once and only once before
    !! calling Wannier-to-Bloch routines.
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : czero
    USE input,         ONLY : use_ws
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nrr_k
    !! Number of Wigner-Size points
    INTEGER, INTENT(in) :: nrr_g
    !! Number of Wigner-Size points
    INTEGER, INTENT(in) :: nmodes
    !! number of modes
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier functions if use_ws == .TRUE. Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: dims2
    !! Is equal to the number of atoms if use_ws == .TRUE. Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: ndegen(dims, nrr_g, dims2)
    !! Degeneracy factor for each R
    COMPLEX(KIND = DP), INTENT(inout) :: array(nbnd, nbnd, nrr_k, nmodes, nrr_g)
    !! Matrix in Wannier representation
    !
    ! Local variables
    INTEGER :: iw
    !! Band index
    INTEGER :: iw2
    !! Band index
    INTEGER :: ir_k, ir_g
    !! WS vectors for electrons.
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: iatm
    !! Counter on atoms
    !
    IF (use_ws) THEN
      !
      DO ir_g = 1, nrr_g
        DO imode = 1, nmodes
          iatm = (imode - 1) / 3 + 1
          DO ir_k = 1, nrr_k
            DO iw2 = 1, nbnd
              DO iw = 1, nbnd
                !
                IF (ndegen(iw, ir_g, iatm) > 0) THEN
                  array(iw, iw2, ir_k, imode, ir_g) &
                  = array(iw, iw2, ir_k, imode, ir_g) / REAL(ndegen(iw, ir_g, iatm), KIND = DP)
                ELSE ! ndegen == 0
                  array(iw, iw2, ir_k, imode, ir_g) = czero
                ENDIF
                !
              ENDDO ! iw
            ENDDO ! iw2
          ENDDO ! ir_k
        ENDDO ! imode
      ENDDO ! ir_g
      !
    ELSE ! .NOT. use_ws
      !
      DO ir_g = 1, nrr_g
        array(:, :, :, :, ir_g) &
        = array(:, :, :, :, ir_g) / REAL(ndegen(1, ir_g, 1), KIND = DP)
      ENDDO ! ir
      !
    ENDIF ! use_ws
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE wigner_divide_ndegen_epmat
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE wigner_divide_ndegen_epmat_dist(array, nbnd, nrr_k, nrr_g, nmodes, ndegen, dims, dims2)
    !---------------------------------------------------------------------------
    !!
    !! Divide an array in the real-space Wannier representation by the degeneracy
    !! factor ndegen. This subroutine must be called once and only once before
    !! calling Wannier-to-Bloch routines.
    !!
    !! JML TODO: Merge with wigner_divide_ndegen_epmat
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : czero
    USE input,         ONLY : use_ws
    USE parallelism,   ONLY : para_bounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT(in) :: nrr_k
    !! Number of Wigner-Size points
    INTEGER, INTENT(in) :: nrr_g
    !! Number of Wigner-Size points
    INTEGER, INTENT(in) :: nmodes
    !! number of modes
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier functions if use_ws == .TRUE. Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: dims2
    !! Is equal to the number of atoms if use_ws == .TRUE. Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: ndegen(dims, nrr_g, dims2)
    !! Degeneracy factor for each R
    COMPLEX(KIND = DP), DIMENSION(:, :, :, :), INTENT(inout) :: array
    !! Matrix in Wannier representation. Size (nbnd, nbnd, nrr_k, ir_stop - ir_start + 1)
    !
    ! Local variables
    INTEGER :: iw
    !! Band index
    INTEGER :: iw2
    !! Band index
    INTEGER :: irn
    !! Combined WS and mode index
    INTEGER :: irn_loc
    !! Combined WS and mode index in the local
    INTEGER :: ir_start, ir_stop
    !! locally start and end points of combined WS and mode index
    INTEGER :: ir_k, ir_g
    !! WS vectors for electrons.
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: iatm
    !! Counter on atoms
    !
    CALL para_bounds(ir_start, ir_stop, nrr_g * nmodes)
    !
    IF (use_ws) THEN
      !
      DO irn = ir_start, ir_stop
        irn_loc = irn - ir_start + 1
        ir_g = (irn - 1)/nmodes + 1
        imode = MOD((irn - 1), nmodes) + 1
        iatm = (imode - 1) / 3 + 1
        DO ir_k = 1, nrr_k
          DO iw2 = 1, nbnd
            DO iw = 1, nbnd
              !
              IF (ndegen(iw, ir_g, iatm) > 0) THEN
                array(iw, iw2, ir_k, irn_loc) &
                = array(iw, iw2, ir_k, irn_loc) / REAL(ndegen(iw, ir_g, iatm), KIND = DP)
              ELSE ! ndegen == 0
                array(iw, iw2, ir_k, irn_loc) = czero
              ENDIF
              !
            ENDDO ! iw
          ENDDO ! iw2
        ENDDO ! ir_k
      ENDDO ! irn
      !
    ELSE ! .NOT. use_ws
      !
      DO irn = ir_start, ir_stop
        irn_loc = irn - ir_start + 1
        ir_g = (irn - 1)/nmodes + 1
        imode = MOD((irn - 1), nmodes) + 1
        array(:, :, :, irn_loc) &
          = array(:, :, :, irn_loc) / REAL(ndegen(1, ir_g, 1), KIND = DP)
      ENDDO ! irn
      !
    ENDIF ! use_ws
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE wigner_divide_ndegen_epmat_dist
    !---------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------
  END MODULE wigner
  !-------------------------------------------------------------------------------------------
