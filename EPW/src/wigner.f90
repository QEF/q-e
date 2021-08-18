  !
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
  IMPLICIT NONE
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
    INTEGER, ALLOCATABLE, INTENT(out) :: irvec_k(:, :)
    !! Integer components of the ir-th Wigner-Seitz grid point in the basis
    !! of the lattice vectors for electrons
    INTEGER, ALLOCATABLE, INTENT(out) :: irvec_q(:, :)
    !! Integer components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, ALLOCATABLE, INTENT(out) :: irvec_g(:, :)
    !! Integer components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, ALLOCATABLE, INTENT(out) :: ndegen_k(:, :, :)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid that depend on
    !! Wannier centers $R + r_n - r_m$
    INTEGER, ALLOCATABLE, INTENT(out) :: ndegen_q(:, :, :)
    !! Wigner-Seitz weights for the phonon grid that depend on
    !! atomic positions $R + \tau(nb) - \tau(na)$
    INTEGER, ALLOCATABLE, INTENT(out) :: ndegen_g(:, :, :)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    !! atomic positions $R + \tau(na) - r_m$
    REAL(KIND = DP), ALLOCATABLE, INTENT(out) :: wslen_k(:)
    !! real-space length for electrons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(out) :: wslen_q(:)
    !! real-space length for phonons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(out) :: wslen_g(:)
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
    INTEGER :: irvec_kk (3, 20 * nkc1 * nkc2 * nkc3)
    !! local integer components of the ir-th Wigner-Seitz grid point
    !! in the basis of the lattice vectors for electrons
    INTEGER :: irvec_qq(3, 20 * nqc1 * nqc2 * nqc3)
    !! local integer components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER :: irvec_gg(3, 20 * nqc1 * nqc2 * nqc3)
    !! local integer components of the ir-th Wigner-Seitz grid point for electron-phonons
    !! We use nkc1 instead of nqc1 because the k-grid is always larger or equal to q-grid.
    INTEGER :: ndegen_kk(20 * nkc1 * nkc2 * nkc3, dims, dims)
    !! local Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER :: ndegen_qq(20 * nqc1 * nqc2 * nqc3, dims2, dims2)
    !! local Wigner-Seitz number of degenerescence (weights) for the phonons grid
    INTEGER :: ndegen_gg(20 * nqc1 * nqc2 * nqc3, dims2, dims)
    !! local Wigner-Seitz number of degenerescence (weights) for the electron-phonons grid
    REAL(KIND = DP) :: wslen_kk(20 * nkc1 * nkc2 * nkc3)
    !! local real-space length for electrons, in units of alat
    REAL(KIND = DP) :: wslen_qq(20 * nqc1 * nqc2 * nqc3)
    !! local real-space length for phonons, in units of alat
    REAL(KIND = DP) :: wslen_gg(20 * nqc1 * nqc2 * nqc3)
    !! local real-space length for electron-phonon, in units of alat
    !
    !  Check the bounds
    IF (nqc1 > nkc1 .OR. nqc2 > nkc2 .OR. nqc3 > nkc3 ) call errore &
       ('wigner_seitz_wrap', ' the phonon grid should be smaller than electron grid', 1)
    !
    ! Now generated the un-sorted points for the electrons, phonons and electron-phonon
    !
    ! If dims > 1, it includes the position of Wannier-Centers
    CALL wigner_seitzkq(nkc1, nkc2, nkc3, irvec_kk, ndegen_kk, wslen_kk, nrr_k, w_centers, dims)
    CALL wigner_seitzkq(nqc1, nqc2, nqc3, irvec_qq, ndegen_qq, wslen_qq, nrr_q, tau, dims2)
    CALL wigner_seitzg(nqc1, nqc2, nqc3, irvec_gg, ndegen_gg, wslen_gg, nrr_g, w_centers, tau, dims, dims2)
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
    !-----------------------------------------------------------------------------
    END SUBROUTINE wigner_seitz_wrap
    !-----------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------------
    SUBROUTINE wigner_seitzkq(nc1, nc2, nc3, irvec, ndegen, wslen, nrr, shift, dims)
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
    USE constants_epw, ONLY : eps6
    USE low_lvl,       ONLY : hpsort_eps_epw
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
    INTEGER, ALLOCATABLE :: ind(:)
    !! Index of sorting
    INTEGER :: nind
    !! The metric tensor
    INTEGER :: ierr
    !! Error status
    INTEGER :: nrr_tmp(dims, dims)
    !! Temporary WS matrix
    INTEGER :: irvec_tmp(3, 20 * nc1 * nc2 * nc3, dims, dims)
    !! Temp
    INTEGER :: ndegen_tmp(20 * nc1 * nc2 * nc3, dims, dims)
    !! Number of degenerate points
    REAL(KIND = DP) :: adot(3, 3)
    !! Dot product between lattice vector
    REAL(KIND = DP) :: dist(125)
    !! Contains the distance squared |r-R|^2
    REAL(KIND = DP) :: mindist
    !! Minimum distance
    REAL(KIND = DP) :: tot, tot2
    !! Sum of all the weigths
    REAL(KIND = DP) :: ndiff(3)
    !! Distances. Must be now real because of fractional atoms
    !
    nind = 20 * nc1 * nc2 * nc3
    IF (nind < 125) THEN
      nind = 125
    ENDIF
    ALLOCATE(ind(nind), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error allocating ind', 1)
    !
    DO ipol = 1, 3
     DO jpol = 1, 3
       adot(ipol, jpol) = DOT_PRODUCT(at(:, ipol), at(:, jpol))
     ENDDO
    ENDDO
    !
    CALL cryst_to_cart(dims, shift(:, :), bg, -1)
    !
    ndegen_tmp(:, :, :) = 0
    irvec_tmp(:, :, :, :) = 0
    DO iw = 1, dims
      DO iw2 = 1, dims
        nrr_tmp(iw, iw2) = 0
        DO n1 = -2 * nc1, 2 * nc1
          DO n2 = -2 * nc2, 2 * nc2
            DO n3 = -2 * nc3, 2 * nc3
              ! Loop over the 5^3 = 125 points R. R=0 corresponds to i1=i2=i3=2, or icnt=63
              i = 0
              dist(:) = 0.d0
              DO i1 = -2, 2
                DO i2 = -2, 2
                  DO i3 = -2, 2
                    i = i + 1
                    ! Calculate distance squared |r-R|^2
                    ndiff(1) = n1 - i1 * nc1 + shift(1, iw2) - shift(1, iw)
                    ndiff(2) = n2 - i2 * nc2 + shift(2, iw2) - shift(2, iw)
                    ndiff(3) = n3 - i3 * nc3 + shift(3, iw2) - shift(3, iw)
                    !ndiff(1) = n1 - i1*nc1 + (shift(1,iw2) + shift(1,iw)) / 2.d0
                    !ndiff(2) = n2 - i2*nc2 + (shift(2,iw2) + shift(2,iw)) / 2.d0
                    !ndiff(3) = n3 - i3*nc3 + (shift(3,iw2) + shift(3,iw)) / 2.d0
                    DO ipol = 1, 3
                      DO jpol = 1, 3
                        dist(i) = dist(i) + DBLE(ndiff(ipol)) * adot(ipol, jpol) * DBLE(ndiff(jpol))
                      ENDDO
                    ENDDO
                    !
                  ENDDO ! i3
                ENDDO ! i2
              ENDDO ! i1
              !IF(iw==iw2) print*,'iw dist ',iw,dist(63)
              !
              ! Sort the 125 vectors R by increasing value of |r-R|^2
              ind(1) = 0 ! required for hpsort_eps (see the subroutine)
              CALL hpsort_eps_epw(125, dist, ind, eps6)
              !
              ! Find all the vectors R with the (same) smallest |r-R|^2;
              ! if R=0 is one of them, then the current point r belongs to
              ! Wignez-Seitz cell => set found to true
              !
              found = .FALSE.
              i = 1
              mindist = dist(1)
              DO WHILE (ABS(dist(i) - mindist) < eps6 .AND. i < 125)
                IF (ind(i) == 63) found = .TRUE.
                i = i + 1
              ENDDO
              !
              IF (found) THEN
                nrr_tmp(iw, iw2) = nrr_tmp(iw, iw2) + 1
                ndegen_tmp(nrr_tmp(iw, iw2), iw, iw2) = i - 1
                irvec_tmp(:, nrr_tmp(iw, iw2), iw, iw2) = (/n1, n2, n3/)
              ENDIF
            ENDDO ! n3
          ENDDO ! n2
        ENDDO ! n1
      ENDDO ! iw2
    ENDDO ! iw
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
    DEALLOCATE(ind, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzkq', 'Error deallocating ind', 1)
    !
    !-----------------------------------------------------------------------------
    END SUBROUTINE wigner_seitzkq
    !-----------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------
    SUBROUTINE wigner_seitzg(nc1, nc2, nc3, irvec, ndegen, wslen, nrr, w_centers, tau, dims, dims2)
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
    USE constants_epw, ONLY : eps6
    USE low_lvl,       ONLY : hpsort_eps_epw
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
    !
    ! Local variables
    LOGICAL :: found
    !! True if the vector has been foun
    INTEGER :: n1, n2, n3
    !! Index for the larger supercell
    INTEGER :: i1, i2, i3
    !! Index to compute |r-R| distance
    INTEGER :: na, iw, iw2
    !! Atom index
    INTEGER :: i
    !! Iterative index
    INTEGER :: ir
    !! Iterative index on the pair of atoms
    INTEGER :: irtot
    !! Iterative index on the combined pair of atoms
    INTEGER :: ipol, jpol
    !! Cartesian direction
    INTEGER, ALLOCATABLE :: ind(:)
    !! Index of sorting
    INTEGER :: nind
    !! The metric tensor
    INTEGER :: ierr
    !! Error status
    INTEGER :: nrr_tmp(dims2, dims)
    !! Temporary array that contains the max number of WS vectors
    !! for a pair of atoms.
    INTEGER :: irvec_tmp(3, 20 * nc1 * nc2 * nc3, dims2, dims)
    !! Temporary WS vectors for each atoms
    INTEGER :: ndegen_tmp(20 * nc1 * nc2 * nc3, dims2, dims)
    !! Temporary WS vectors weigths for each atoms
    REAL(KIND = DP) :: adot(3, 3)
    !! Dot product between lattice vector
    REAL(KIND = DP) :: dist(125)
    !! Contains the distance squared |r-R|^2
    REAL(KIND = DP) :: mindist
    !! Minimum distance
    REAL(KIND = DP) :: tot, tot2
    !! Sum of all the weigths
    REAL(KIND = DP) :: ndiff(3)
    !! Distances. Must be now real because of fractional atoms
    !
    nind = 20 * nc1 * nc2 * nc3
    IF (nind < 125) THEN
      nind = 125
    ENDIF
    ALLOCATE(ind(nind), STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzg', 'Error allocating ind', 1)
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
    DO iw = 1, dims
      DO na = 1, dims2
        DO n1 = -2 * nc1, 2 * nc1
          DO n2 = -2 * nc2, 2 * nc2
            DO n3 = -2 * nc3, 2 * nc3
              !
              ! Loop over the 5^3 = 125 points R. R=0 corresponds to i1=i2=i3=2, or icnt=63
              i = 0
              dist(:) = 0.d0
              DO i1 = -2, 2
                DO i2 = -2, 2
                  DO i3 = -2, 2
                    i = i + 1
                    !
                    ! Calculate distance squared |r-R|^2
                    ndiff(1) = n1 - i1 * nc1 + tau(1, na) - w_centers(1, iw)
                    ndiff(2) = n2 - i2 * nc2 + tau(2, na) - w_centers(2, iw)
                    ndiff(3) = n3 - i3 * nc3 + tau(3, na) - w_centers(3, iw)
                    !ndiff(1) = n1 - i1 * nc1 - tau(1, na) - w_centers(1, iw)
                    !ndiff(2) = n2 - i2 * nc2 - tau(2, na) - w_centers(2, iw)
                    !ndiff(3) = n3 - i3 * nc3 - tau(3, na) - w_centers(3, iw)
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
              ! Sort the 125 vectors R by increasing value of |r-R|^2
              ind(1) = 0 ! required for hpsort_eps (see the subroutine)
              CALL hpsort_eps_epw(125, dist, ind, eps6)
              !
              ! Find all the vectors R with the (same) smallest |r-R|^2;
              ! if R=0 is one of them, then the current point r belongs to
              ! Wignez-Seitz cell => set found to true
              found = .FALSE.
              i = 1
              mindist = dist(1)
              DO WHILE (ABS(dist(i) - mindist) < eps6 .AND. i < 125)
                IF (ind(i) == 63) found = .TRUE.
                i = i + 1
              ENDDO
              !
              IF (found) THEN
                nrr_tmp(na, iw) = nrr_tmp(na, iw) + 1
                ndegen_tmp(nrr_tmp(na, iw), na, iw) = i - 1
                irvec_tmp(:, nrr_tmp(na, iw), na, iw) = (/n1, n2, n3/)
              ENDIF
            ENDDO ! n3
          ENDDO ! n2
        ENDDO ! n3
      ENDDO ! na
    ENDDO ! iw
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
    DEALLOCATE(ind, STAT = ierr)
    IF (ierr /= 0) CALL errore('wigner_seitzg', 'Error deallocating ind', 1)
    !
    !------------------------------------------------------------------------------------------
    END SUBROUTINE wigner_seitzg
    !------------------------------------------------------------------------------------------
    !
    ! -----------------------------------------------------------------------------------------
    SUBROUTINE backtows(v, ws, rws, nrwsx, nrws)
    ! -----------------------------------------------------------------------------------------
    !! 03/2019 - F. Macheda and S. Ponce
    !! This subroutines takes a vector "v" like a k-point and bring it back to the
    !! first Wigner-Seitz cell.
    ! -----------------------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : bg
    USE constants_epw, ONLY : zero, eps8
    USE low_lvl,       ONLY : find_minimum
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: v(3)
    !! Input vector
    INTEGER, INTENT(in) :: nrws
    !! Number of WS vectors
    INTEGER, INTENT(in) :: nrwsx
    !! Maximum number of WS vector
    REAL(KIND = DP), INTENT(in) :: rws(4, nrwsx)
    !! List of WS vectors
    REAL(KIND = DP), INTENT(out) :: ws(3)
    !! Vector back into the first WS cell.
    !
    ! Local variables
    INTEGER :: L, M, N
    !! l,m,n directions
    INTEGER :: i
    !! Index on number of WS vector
    INTEGER :: minpos
    !! Minimum position
    INTEGER, PARAMETER :: nn = 3
    !! number of neighbours
    REAL(KIND = DP) :: dist
    !! Distance
    REAL(KIND = DP) :: distances(nrws)
    !! Dsitance array
    !
    ws(:) = zero
    distances(:) = zero
    !
    alpha : DO L = -nn, nn
      DO M = -nn, nn
        DO N = -nn, nn
          ws(:) = v(:) + L * bg(:, 1) + M * bg(:, 2) + N * bg(:, 3)
          DO i = 1, nrws
            distances(i) = DSQRT((ws(1) - rws(2, i))**2 + (ws(2) - rws(3, i))**2 + (ws(3) - rws(4, i))**2)
          ENDDO
          minpos = find_minimum(distances, nrws)
          dist = DSQRT(ws(1)**2 + ws(2)**2 + ws(3)**2)
          IF (dist < distances(minpos) + eps8) EXIT alpha
        ENDDO ! N
      ENDDO ! M
    ENDDO alpha ! L
    !
    !-----------------------------------------------------------------------------------------
    END SUBROUTINE backtoWS
    !-----------------------------------------------------------------------------------------
    !
  !-------------------------------------------------------------------------------------------
  END MODULE wigner
  !-------------------------------------------------------------------------------------------
