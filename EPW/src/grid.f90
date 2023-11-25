  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE grid
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines related to k-point or q-point grid
  !! generation as well as selection of k/q points.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE loadkmesh_para()
    !-----------------------------------------------------------------------
    !!
    !! Load fine k mesh and distribute among pools
    !!
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode_id, stdout
    USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
    USE mp,        ONLY : mp_bcast, mp_sum
    USE mp_world,  ONLY : mpime
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : filkf, nkf1, nkf2, nkf3, iterative_bte, &
                          rand_k, rand_nk, mp_mesh_k, system_2d, eig_read, vme, &
                          scell_mat_plrn, scell_mat, as, bs
    USE elph2,     ONLY : nkqtotf, nkqf, xkf, wkf, nkf, xkfd, deltaq, &
                          xkf_irr, wkf_irr, bztoibz, s_bztoibz
    USE cell_base, ONLY : at, bg
    USE symm_base, ONLY : s, t_rev, time_reversal, nsym
    USE io_var,    ONLY : iunkf, iunRpscell, iunkgridscell
    USE low_lvl,   ONLY : init_random_seed, matinv3
    USE constants_epw, ONLY : eps4, eps8
    USE noncollin_module, ONLY : noncolin
# if defined(__MPI)
    USE parallel_include, ONLY : MPI_INTEGER2
# endif
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 10) :: coordinate_type
    !! filkf coordinate type (crystal or cartesian)
    LOGICAL, EXTERNAL :: imatches
    !! Regex matching text.
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: idir
    !! Crystal direction (G-vector)
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: i, j, k
    !! Counter on the k-point index along nkf1, nkf2, nkf3
    INTEGER :: rest
    !! rest from the division of nr of q-points over pools
    INTEGER :: ierr
    !! Error status
    INTEGER :: iRp1, iRp2, iRp3, Rpmax, nRp
    !! Number of unit cells within supercell
    INTEGER :: Rp_crys_p(3)
    !! Unit cell vectors in primitive crystal coordinates
    INTEGER, ALLOCATABLE :: Rp(:, :)
    !! List of unit cell vectors within supercell in primitive crystal coords
    INTEGER :: iGs1, iGs2, iGs3, Gsmax, nGs
    !! Number of supercell G-vectors within primitive reciprocal unit cell
    INTEGER :: Gs_crys_s(3)
    !! Supercell G-vectors in supercell reciprocal coordinates
    REAL(KIND = DP) :: ap(3, 3), bp(3, 3)
    !! Auxiliary definitions of real and reciprocal primitive cell vector matrix
    REAL(KIND = DP), ALLOCATABLE :: xkf_(:, :)
    !! coordinates k-points
    REAL(KIND = DP), ALLOCATABLE :: wkf_(:)
    !! weights k-points
    REAL(KIND = DP) :: scell_mat_b(3, 3)
    !! Reciprocal lattice transformation matrix
    REAL(KIND = DP) :: p2s(3, 3), bs2p(3, 3)
    !! Transformation matrix from primitive to supercell crystal coordinates
    REAL(KIND = DP) :: Rp_crys_s(3)
    !! Unit cell vectors in supercell crystal coordinates
    REAL(KIND = DP) :: Gs_crys_p(3)
    !! Supercell G-vectors in primitive crystal coordinates
    REAL(KIND = DP), ALLOCATABLE :: Gs(:, :)
    !! Supercell G-vectors within primitive reciprocal unit cell
    !
    IF (mpime == ionode_id) THEN
      IF (filkf /= '') THEN ! load from file
        !
        WRITE(stdout, *) '    Using k-mesh file: ', TRIM(filkf)
        OPEN(UNIT = iunkf, FILE = filkf, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('loadkmesh_para', 'opening file ' // filkf, ABS(ios))
        READ(iunkf, *) nkqtotf, coordinate_type
        ! Default
        IF (TRIM(coordinate_type) .EQ. ' ') coordinate_type = 'crystal'
        IF (.NOT. imatches("crystal", coordinate_type) .AND. .NOT. imatches("cartesian", coordinate_type)) &
           CALL errore('loadkmesh_para', 'ERROR: Specify either crystal or cartesian coordinates in the filkf file', 1)
        !
        ALLOCATE(xkf_(3, 2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_', 1)
        ALLOCATE(wkf_(2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
        !
        DO ik = 1, nkqtotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          READ(iunkf, *) xkf_(:, ikk ), wkf_(ikk)
          !
          ! SP: This is so we can input a weight of 1 to random file
          !     This way you can feed the same file for the k and q grid
          wkf_(ikk) = wkf_(ikk) * 2.d0
          xkf_(:, ikq) = xkf_(:, ikk)
          wkf_(ikq) = 0.d0
          !
        ENDDO
        CLOSE(iunkf)
        !
        ! redefine nkqtotf to include the k+q points
        !
        nkqtotf = 2 * nkqtotf
        !
        IF (imatches("cartesian", coordinate_type)) THEN
          CALL cryst_to_cart(nkqtotf, xkf_, at, -1)
        ENDIF
        !
      !JLB
      ELSEIF (scell_mat_plrn) THEN
        !
        WRITE(stdout, '(a)') ' '
        WRITE(stdout, '(a)') '     Supercell transformation activated (k), as=S*at'
        WRITE(stdout, '(a,3i4)') '     S(1, 1:3): ', scell_mat(1, 1:3)
        WRITE(stdout, '(a,3i4)') '     S(2, 1:3): ', scell_mat(2, 1:3)
        WRITE(stdout, '(a,3i4)') '     S(3, 1:3): ', scell_mat(3, 1:3)
        !
        ap = TRANSPOSE(at)
        as = MATMUL(scell_mat,ap)
        !
        WRITE(stdout, '(a)') '     Transformed lattice vectors (alat units):'
        WRITE(stdout, '(a,3f12.6)') '     as(1, 1:3): ', as(1, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     as(2, 1:3): ', as(2, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     as(3, 1:3): ', as(3, 1:3)
        !
        scell_mat_b = matinv3(REAL(scell_mat, DP))
        scell_mat_b = TRANSPOSE(scell_mat_b)
        !
        WRITE(stdout, '(a)') '     Reciprocal lattice transformation matrix, Sbar = (S^{-1})^{t}:'
        WRITE(stdout, '(a,3f12.6)') '     Sbar(1, 1:3): ', scell_mat_b(1, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     Sbar(2, 1:3): ', scell_mat_b(2, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     Sbar(3, 1:3): ', scell_mat_b(3, 1:3)
        !
        bp = TRANSPOSE(bg)
        bs = MATMUL(scell_mat_b, bp)
        !
        WRITE(stdout, '(a)') '     Transformed reciprocal lattice vectors (2pi/alat units):'
        WRITE(stdout, '(a,3f12.6)') '     bs(1, 1:3): ', bs(1, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     bs(2, 1:3): ', bs(2, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     bs(3, 1:3): ', bs(3, 1:3)
        WRITE(stdout, '(a)') '  '
        !
        ! Define transformation matrix from primitive crystal coordinates
        ! to supercell crystal coordinates Rp_crys_s = ((a_s)^{T})^{-1} (a_p)^{T} Rp_crys_p
        p2s = matinv3(TRANSPOSE(as))
        p2s = MATMUL(p2s,TRANSPOSE(ap))
        !
        ! Find how many unit cells are contained within the supercell
        Rpmax = 5*MAXVAL(scell_mat) ! This should be large enough to find all
        ALLOCATE(Rp(3, Rpmax**3), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating Rp', 1)
        Rp  = 0
        nRp = 0
        DO iRp1 = -Rpmax, Rpmax
          DO iRp2 = -Rpmax, Rpmax
            DO iRp3 = -Rpmax, Rpmax
              Rp_crys_p = (/iRp1, iRp2, iRp3/)
              Rp_crys_s = MATMUL(p2s, Rp_crys_p)
              ! Unit cell within supercell if Rp\in(0,1) in supercell crystal coordinates
              IF (ALL(Rp_crys_s > -eps8) .AND. ALL(Rp_crys_s < 1.d0-eps8)) THEN
                nRp = nRp + 1
                Rp(1:3, nRp) = Rp_crys_p
              END IF
            END DO
          END DO
        END DO
        WRITE(stdout, '(a, 3i6)') '     Number of unit cells within supercell:', nRp
        !
        ! Write Rp-s in supercell to file
        IF (mpime == ionode_id) THEN
          OPEN(UNIT = iunRpscell, FILE = 'Rp.scell.plrn', ACTION = 'write')
          WRITE(iunRpscell, *) nRp
          DO iRp1 = 1, nRp
            WRITE(iunRpscell, *) Rp(1:3, iRp1)
          END DO
          CLOSE(iunRpscell)
        ENDIF
        !
        IF (ALLOCATED(Rp)) DEALLOCATE(Rp)
        !
        ! Define transformation matrix from reciprocal supercell crystal coordinates
        ! to reciprocal primitive crystal coordinates
        ! Gs_crys_p = ((bp)^{T})^{-1} (bs)^{T} Gs_crys_s
        bs2p = matinv3(TRANSPOSE(bp))
        bs2p = MATMUL(bs2p, TRANSPOSE(bs))
        !
        ! Find how many k-points lie within primitive reciprocal cell
        Gsmax = Rpmax ! This should be large enough to find all
        ALLOCATE(Gs(3, Gsmax**3), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating Gs', 1)
        Gs = 0.d0
        nGs = 0
        DO iGs1 = -Gsmax, Gsmax
          DO iGs2 = -Gsmax, Gsmax
            DO iGs3 = -Gsmax, Gsmax
              Gs_crys_s = (/iGs1, iGs2, iGs3/)
              Gs_crys_p = MATMUL(bs2p, Gs_crys_s)
              ! Gs within primitive reciprocal unit cell if Gs\in(0,1) in crys_p coords.
              IF (ALL(Gs_crys_p > -eps8) .AND. ALL(Gs_crys_p < 1.d0-eps8)) THEN
                nGs = nGs + 1
                Gs(1:3, nGs) = Gs_crys_p
              END IF
            END DO
          END DO
        END DO
        WRITE(stdout, '(a, 3i6)') '     Number of k-points needed:', nGs
        !DO iGs1 = 1, nGs
        !  WRITE(stdout, '(3f12.6)') Gs(1:3, iGs1)
        !END DO
        !
        ! Write Gs-s within unit cell BZ (k-grid) to file
        IF (mpime == ionode_id) THEN
          OPEN(UNIT = iunkgridscell, FILE = 'kgrid.scell.plrn', ACTION = 'write')
          WRITE(iunkgridscell, *) nGs
          DO iGs1 = 1, nGs
            WRITE(iunkgridscell, '(3f12.6)') Gs(1:3, iGs1)
          END DO
          CLOSE(iunkgridscell)
        ENDIF
        !
        ! Save list of needed k-points
        nkqtotf = nGs
        ALLOCATE(xkf_(3, 2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_', 1)
        ALLOCATE(wkf_(2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
        !
        DO ik = 1, nkqtotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          xkf_(:, ikk) = Gs(1:3, ik)
          wkf_(ikk) = 1.d0 ! weight not important for polaron
          !
          xkf_(:, ikq) = xkf_(:, ikk)
          wkf_(ikq) = 0.d0
          !
        ENDDO
        !
        ! redefine nkqtotf to include the k+q points
        nkqtotf = 2 * nkqtotf
        !
        IF (ALLOCATED(Gs)) DEALLOCATE(Gs)
        !
      !JLB
      ELSEIF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_k) THEN
          ! get size of the mp_mesh in the irr wedge
          WRITE(stdout, '(a,3i5)') '     Using uniform MP k-mesh: ', nkf1, nkf2, nkf3
          !
          ! The call to this routine computes the IBZ point xkf_irr, wkf_irr and
          ! returns the number of irr points nkqtotf
          ! xkf_irr and wkf_irr are allocated inside with dimension nkqtotf
          ! xkf_irr is in crystal coordinate
          CALL kpoint_grid_epw(nsym, time_reversal, s, t_rev, nkf1, nkf2, nkf3, nkqtotf)
          !
          ALLOCATE(xkf_(3, 2 * nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_', 1)
          ALLOCATE(wkf_(2 * nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
          !
          DO ik = 1, nkqtotf
            ikk = 2 * ik - 1
            ikq = ikk + 1
            xkf_(:, ikk) = xkf_irr(:, ik)
            xkf_(:, ikq) = xkf_irr(:, ik)
            wkf_(ikk)   = 2.d0 * wkf_irr(ik)
            wkf_(ikq)   = 0.d0
          ENDDO
          !
          IF (iterative_bte) THEN
            ! Fold the points in the region [0-1] from the region -0.5,0.5
            DO ik = 1, 2 * nkqtotf
              DO idir= 1, 3
                IF (xkf_(idir, ik) < 0.0d0) THEN
                  xkf_(idir, ik) = xkf_(idir, ik) + 1.0d0
                ENDIF
              ENDDO
            ENDDO
          ENDIF
          !
          ! redefine nkqtotf to include the k+q points
          !
          nkqtotf = 2 * nkqtotf
          !
        ELSE ! mp_mesh_k
          !
          WRITE(stdout, '(a,3i5)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          !
          nkqtotf = 2 * nkf1 * nkf2 * nkf3
          ALLOCATE(xkf_(3, nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_ ', 1)
          ALLOCATE(wkf_(nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
          wkf_(:) = 0.d0
          DO ik = 1, nkf1 * nkf2 * nkf3
            wkf_(2 * ik - 1) = 2.d0 / DBLE(nkqtotf / 2)
          ENDDO
          DO i = 1, nkf1
            DO j = 1, nkf2
              DO k = 1, nkf3
                ik = (i - 1) * nkf2 * nkf3 + (j - 1) * nkf3 + k
                ikk = 2 * ik - 1
                ikq = ikk + 1
                xkf_(1, ikk) = DBLE(i - 1) / DBLE(nkf1)
                xkf_(2, ikk) = DBLE(j - 1) / DBLE(nkf2)
                xkf_(3, ikk) = DBLE(k - 1) / DBLE(nkf3)
                xkf_(1, ikq) = xkf_(1, ikk)
                xkf_(2, ikq) = xkf_(2, ikk)
                xkf_(3, ikq) = xkf_(3, ikk)
              ENDDO
            ENDDO
          ENDDO
        ENDIF !mp_mesh_k
        !
      ELSEIF (rand_k) THEN  ! random points
        !
        WRITE(stdout, *) '     Using random k-mesh: ', rand_nk
        !
        nkqtotf = rand_nk
        ALLOCATE(xkf_(3, 2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_', 1)
        ALLOCATE(wkf_(2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
        !
        CALL init_random_seed()
        !
        DO ik = 1, nkqtotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          wkf_(ikk) = 2.d0 / DBLE(nkqtotf)
          wkf_(ikq) = 0.d0
          !
          IF (system_2d == 'no') THEN
            CALL random_number(xkf_(:, ikk))
          ELSE
            CALL random_number(xkf_(1:2, ikk))
            xkf_(3, ikk) = 0.d0
          ENDIF
          xkf_(:, ikq) = xkf_(:, ikk)
        ENDDO
        !
        ! redefine nkqtotf to include the k+q points
        !
        nkqtotf = 2 * nkqtotf
        !
      ELSE ! don't know how to get grid
        CALL errore('loadkmesh_para', "Cannot load fine k points", 1)
      ENDIF
    ENDIF
    !
#if defined(__MPI)
    CALL mp_bcast(nkqtotf, ionode_id, inter_pool_comm)
    !
    !  scatter the k points of the fine mesh across the pools
    !
    nkqf = 2 * (nkqtotf / (2 * npool))
    rest = (nkqtotf - nkqf * npool) / 2
    IF (my_pool_id < rest) THEN
      nkqf = nkqf + 2
      lower_bnd = my_pool_id * nkqf + 1
      upper_bnd = lower_bnd + nkqf - 1
    ELSE
      lower_bnd = rest * (nkqf + 2) + (my_pool_id - rest) * nkqf + 1
      upper_bnd = lower_bnd + nkqf - 1
    ENDIF
    !
    nkf = nkqf / 2
    IF (mpime /= ionode_id) THEN
      ALLOCATE(xkf_(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_', 1)
      ALLOCATE(wkf_(nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
      IF (mp_mesh_k) THEN
        ALLOCATE(bztoibz(nkf1 * nkf2 * nkf3), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating bztoibz', 1)
        ALLOCATE(s_bztoibz(nkf1 * nkf2 * nkf3), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating bztoibz', 1)
      ENDIF
    ENDIF
    CALL mp_bcast(xkf_, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkf_, ionode_id, inter_pool_comm)
    IF (mp_mesh_k) THEN
      CALL mp_bcast(bztoibz, ionode_id, inter_pool_comm)
      CALL mp_bcast(s_bztoibz, ionode_id, inter_pool_comm)
      ! The mp_bcast wrapper cannot broadcast integer*2 kinds.
      !CALL MPI_TYPE_CREATE_F90_INTEGER(SIK2, int2type, ierr)
      !CALL MPI_BCAST(s_bztoibz, nkf1 * nkf2 * nkf3, int2type, ionode_id, inter_pool_comm)
      !CALL MPI_BCAST(s_bztoibz, nkf1 * nkf2 * nkf3, MPI_INTEGER2, ionode_id, inter_pool_comm)
    ENDIF
    !DO ik=1, nkf1 * nkf2 * nkf3
    !  IF (mpime == 0) write(900, *) ik, bztoibz(ik), s_bztoibz(ik)
    !  IF (mpime == 0) FLUSH(900)
    !  IF (mpime == 1) write(901, *) ik, bztoibz(ik), s_bztoibz(ik)
    !  IF (mpime == 1) FLUSH(901)
    !ENDDO
    !
#else
    ! In serial the definitions are much easier
    !
    nkqf = nkqtotf
    nkf = nkqf / 2
    lower_bnd = 1
    upper_bnd = nkqf
    !
#endif
    !
    ! Assign the weights and vectors to the correct bounds
    !
    ALLOCATE(xkf(3, nkqf), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf', 1)
    ALLOCATE(wkf(nkqf), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf', 1)
    xkf(:, :) = xkf_(:, lower_bnd:upper_bnd)
    !
    ! KMB: set coordinates of displaced vectors for indabs
    IF (vme == 'wannier' .AND. eig_read) THEN
      ALLOCATE(xkfd(3, nkqf, 6), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkfd', 1)
      deltaq = 0.001d0
      DO ik = 1, nkqf
        !--bring the k point to cartesian coordinates
        CALL cryst_to_cart(1, xkf(:, ik), bg, 1)
        xkfd(:, ik, 1) = xkf(:, ik) + (/deltaq, 0.d0, 0.d0/)
        xkfd(:, ik, 2) = xkf(:, ik) - (/deltaq, 0.d0, 0.d0/)
        xkfd(:, ik, 3) = xkf(:, ik) + (/0.d0, deltaq, 0.d0/)
        xkfd(:, ik, 4) = xkf(:, ik) - (/0.d0, deltaq, 0.d0/)
        xkfd(:, ik, 5) = xkf(:, ik) + (/0.d0, 0.d0, deltaq/)
        xkfd(:, ik, 6) = xkf(:, ik) - (/0.d0, 0.d0, deltaq/)
        !  bring the k point to crystal coordinates
        CALL cryst_to_cart(1, xkf(:, ik), at, -1)
        DO i = 1, 6
          CALL cryst_to_cart(1, xkfd(:, ik, i), at, -1)
        ENDDO
      ENDDO
    ENDIF
    IF (noncolin) THEN
      wkf(:) = wkf_(lower_bnd:upper_bnd) / 2.d0
    ELSE
      wkf(:) = wkf_(lower_bnd:upper_bnd)
    ENDIF
    !
    IF (ABS(SUM(wkf_ (:)) - 2.d0) > eps4) &
      WRITE(stdout, '(5x,"WARNING: k-point weigths do not add up to 1 [loadkmesh_para]")')
    !
    WRITE(stdout, '(5x,"Size of k point mesh for interpolation: ",i10)') nkqtotf
    WRITE(stdout, '(5x,"Max number of k points per pool:",7x,i10)') nkqf
    !
    DEALLOCATE(xkf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error deallocating xkf_', 1)
    DEALLOCATE(wkf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error deallocating wkf_', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE loadkmesh_para
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE loadkmesh_serial()
    !-----------------------------------------------------------------------
    !!
    !! Load fine k mesh in sequential
    !!
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode_id, stdout
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast
    USE mp_world,  ONLY : mpime
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : filkf, nkf1, nkf2, nkf3, &
                          rand_k, rand_nk, mp_mesh_k, system_2d, eig_read, vme
    USE elph2,     ONLY : xkf, wkf, nkqtotf, nkf, nkqf, xkfd, deltaq
    USE cell_base, ONLY : at, bg
    USE symm_base, ONLY : s, t_rev, time_reversal, nsym
    USE io_var,    ONLY : iunkf
    USE low_lvl,   ONLY : init_random_seed
    USE constants_epw, ONLY : eps4
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 10) :: coordinate_type
    !! filkf coordinate type (crystal or cartesian)
    LOGICAL, EXTERNAL  :: imatches
    !! Regex matching text.
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: i, j, k
    !! Counter on the k-point index along nkf1, nkf2, nkf3
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP), ALLOCATABLE :: xkf_tmp(:, :)
    !! coordinates k-points
    REAL(KIND = DP), ALLOCATABLE :: wkf_tmp(:)
    !! weights k-points
    !
    IF (mpime == ionode_id) THEN
      IF (filkf /= '') THEN ! load from file
        !
        ! Each pool gets its own copy from the action=read statement
        !
        WRITE(stdout, *) '     Using k-mesh file: ', TRIM(filkf)
        OPEN(UNIT = iunkf, FILE = filkf, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('loadkmesh_serial', 'opening file ' // filkf, ABS(ios))
        READ(iunkf, *) nkqtotf, coordinate_type
        ! Default
        IF (TRIM(coordinate_type) .EQ. ' ') coordinate_type = 'crystal'
        IF (.NOT. imatches("crystal", coordinate_type) .AND. .NOT. imatches("cartesian", coordinate_type)) &
          CALL errore('loadkmesh_serial', 'ERROR: Specify either crystal or cartesian coordinates in the filkf file', 1)
        !
        ALLOCATE(xkf(3, 2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkf', 1)
        ALLOCATE(wkf(2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating wkf', 1)
        DO ik = 1, nkqtotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          READ(iunkf, *) xkf(:, ikk ), wkf(ikk)
          !
          ! SP: This is so we can input a weight of 1 to random file
          !     This way you can feed the same file for the k and q grid
          wkf(ikk) = wkf(ikk) * 2.d0
          xkf(:, ikq) = xkf(:, ikk)
          wkf(ikq) = 0.d0
          !
        ENDDO
        CLOSE(iunkf)
        !
        ! redefine nkqtotf to include the k+q points
        !
        nkqtotf = 2 * nkqtotf
        !
        IF (imatches("cartesian", coordinate_type)) THEN
          CALL cryst_to_cart(nkqtotf, xkf, at, -1)
        ENDIF
        !
      ELSEIF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_k) THEN
          ! get size of the mp_mesh in the irr wedge
          WRITE(stdout, '(a,3i5)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          !
          ALLOCATE(xkf(3, 2 * nkf1 * nkf2 * nkf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkf', 1)
          ALLOCATE(wkf(2 * nkf1 * nkf2 * nkf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating wkf', 1)
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nsym, time_reversal, .FALSE., s, t_rev, bg, nkf1 * nkf2 * nkf3, &
               0, 0, 0, nkf1, nkf2, nkf3, nkqtotf, xkf, wkf)
          DEALLOCATE(xkf, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error deallocating xkf', 1)
          DEALLOCATE(wkf, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error deallocating wkf', 1)
          ALLOCATE(xkf(3, 2 * nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkf', 1)
          ALLOCATE(wkf(2 * nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating wkf', 1)
          ALLOCATE(xkf_tmp(3, nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkf_tmp', 1)
          ALLOCATE(wkf_tmp(nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating wkf_tmp', 1)
          CALL kpoint_grid(nsym, time_reversal, .FALSE., s, t_rev, bg, nkf1 * nkf2 * nkf3, &
               0, 0, 0, nkf1, nkf2, nkf3, nkqtotf, xkf_tmp, wkf_tmp)
          !
          ! assign to k and k+q for xkf and wkf
          !
          DO ik = 1, nkqtotf
            ikk = 2 * ik - 1
            ikq = ikk + 1
            xkf(:, ikk) = xkf_tmp(:, ik)
            xkf(:, ikq) = xkf_tmp(:, ik)
            wkf(ikk)   = 2.d0 * wkf_tmp(ik)
            wkf(ikq)   = 0.d0
          ENDDO
          DEALLOCATE(xkf_tmp, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error deallocating xkf_tmp', 1)
          DEALLOCATE(wkf_tmp, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error deallocating wkf_tmp', 1)
          !
          ! bring the k point to crystal coordinates
          CALL cryst_to_cart(2 * nkqtotf, xkf, at, -1)
          !
          ! redefine nkqtotf to include the k+q points
          !
          nkqtotf = 2 * nkqtotf
        ELSE
          WRITE (stdout, '(a,3i5)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          !
          nkqtotf = 2 * nkf1 * nkf2 * nkf3
          ALLOCATE(xkf(3, nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkf', 1)
          ALLOCATE(wkf(nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating wkf', 1)
          wkf(:) = 0.d0
          DO ik = 1, nkf1 * nkf2 * nkf3
            wkf(2 * ik - 1) = 2.d0 / DBLE(nkqtotf / 2)
          ENDDO
          DO i = 1, nkf1
            DO j = 1, nkf2
              DO k = 1, nkf3
                ik = (i - 1) * nkf2 * nkf3 + (j - 1) * nkf3 + k
                ikk = 2 * ik - 1
                ikq = ikk + 1
                xkf(1, ikk) = DBLE(i - 1) / DBLE(nkf1)
                xkf(2, ikk) = DBLE(j - 1) / DBLE(nkf2)
                xkf(3, ikk) = DBLE(k - 1) / DBLE(nkf3)
                xkf(1, ikq) = xkf(1, ikk)
                xkf(2, ikq) = xkf(2, ikk)
                xkf(3, ikq) = xkf(3, ikk)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSEIF (rand_k) THEN  ! random points
        WRITE (stdout, *) '    Using random k-mesh: ', rand_nk
        !
        nkqtotf = rand_nk
        ALLOCATE(xkf(3, 2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkf', 1)
        ALLOCATE(wkf(2 * nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating wkf', 1)
        !
        CALL init_random_seed()
        !
        DO ik = 1, nkqtotf
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          wkf(ikk) = 2.d0 / DBLE(nkqtotf)
          wkf(ikq) = 0.d0
          !
          IF (system_2d == 'no') THEN
            CALL random_number(xkf(:, ikk))
          ELSE
            CALL random_number(xkf(1:2, ikk))
            xkf(3, ikk) = 0.d0
          ENDIF
          !
          xkf(:, ikq) = xkf(:, ikk)
          !
        ENDDO
        !
        ! redefine nkqtotf to include the k+q points
        !
        nkqtotf = 2 * nkqtotf
        !
      ELSE ! don't know how to get grid
        CALL errore('loadkmesh_serial', "Cannot load fine k points", 1)
      ENDIF
      !
      ! Serial
      nkf = nkqtotf / 2
      nkqf = nkqtotf
      !
    ENDIF
    !
    CALL mp_bcast(nkf, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkqf, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkqtotf, ionode_id, inter_pool_comm)
    IF (mpime /= ionode_id) THEN
      ALLOCATE(xkf(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkf', 1)
      ALLOCATE(wkf(nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating wkf', 1)
    ENDIF
    CALL mp_bcast(xkf, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkf, ionode_id, inter_pool_comm)
    !
    ! KMB: set coordinates of displaced vectors - indabs
    IF (vme == 'wannier' .AND. eig_read) THEN
      ALLOCATE(xkfd(3, nkqf, 6), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkfd', 1)
      deltaq = 0.001d0
      DO ik = 1, nkqf
        ! Bring the k point to cartesian coordinates
        CALL cryst_to_cart(1, xkf(:, ik), bg, 1)
        xkfd(:, ik, 1) = xkf(:, ik) + (/deltaq, 0.d0, 0.d0/)
        xkfd(:, ik, 2) = xkf(:, ik) - (/deltaq, 0.d0, 0.d0/)
        xkfd(:, ik, 3) = xkf(:, ik) + (/0.d0, deltaq, 0.d0/)
        xkfd(:, ik, 4) = xkf(:, ik) - (/0.d0, deltaq, 0.d0/)
        xkfd(:, ik, 5) = xkf(:, ik) + (/0.d0, 0.d0, deltaq/)
        xkfd(:, ik, 6) = xkf(:, ik) - (/0.d0, 0.d0, deltaq/)
        ! Bring the k point to crystal coordinates
        CALL cryst_to_cart(1, xkf(:, ik), at, -1)
        DO i = 1, 6
          CALL cryst_to_cart(1, xkfd(:, ik, i), at, -1)
        ENDDO
      ENDDO
    ENDIF
    IF (ABS(SUM(wkf) - 2.d0) > eps4) &
      WRITE(stdout,'(5x,"WARNING: k-point weigths do not add up to 1 [loadkmesh_serial]")')
    !
    WRITE(stdout, '(5x,"Size of k point mesh for interpolation: ",i10)') nkqtotf
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE loadkmesh_serial
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE loadkmesh_fullBZ()
    !-----------------------------------------------------------------------
    !!
    !!  Create a k-mesh for fine grid without symmetries on the full grid
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP, sgl
    USE epwcom,        ONLY : nkf1, nkf2, nkf3
    USE elph2,         ONLY : xkf_bz
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: ik
    !! K-point index
    INTEGER :: i, j, k
    !! K-point grid dim
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(xkf_bz(3, nkf1 * nkf2 * nkf3), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fullBZ', 'Error allocating xkf_bz', 1)
    xkf_bz(:, :) = zero
    !
    IF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0)) THEN
      DO i = 1, nkf1
        DO j = 1, nkf2
          DO k = 1, nkf3
            ik = (i - 1) * nkf2 * nkf3 + (j - 1) * nkf3 + k
            xkf_bz(1, ik) = DBLE(i - 1) / DBLE(nkf1)
            xkf_bz(2, ik) = DBLE(j - 1) / DBLE(nkf2)
            xkf_bz(3, ik) = DBLE(k - 1) / DBLE(nkf3)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !-----------------------------------------------------------------------
    END SUBROUTINE loadkmesh_fullBZ
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE loadkmesh_fst(nrr_k, dims, irvec_k, ndegen_k, nelec)
    !-----------------------------------------------------------------------
    !!
    !!  Load fine k mesh and distribute among pools
    !!  We only load the k-points that fall within the fsthick
    !!  This is useful to reduce computational cost in mobility calculations for example
    !!  Only support homogeneous grids.
    !!
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : ionode_id, stdout
    USE mp_global,        ONLY : inter_pool_comm, my_pool_id, npool
    USE mp,               ONLY : mp_bcast, mp_sum, mp_barrier
    USE epwcom,           ONLY : nkf1, nkf2, nkf3, iterative_bte
    USE elph2,            ONLY : wkf_fst, xkf_fst, nkqf, xkf, wkf, nkf, nkqtotf
    USE symm_base,        ONLY : s, t_rev, nsym
    USE constants_epw,    ONLY : byte2Mb, eps4, zero
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Coordinates of real space vector for electrons
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: nelec
    !! Number of electrons
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: idir
    !! Crystal direction (G-vector)
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: rest
    !! rest from the division of nr of q-points over pools
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP), ALLOCATABLE :: xkf_(:, :)
    !! coordinates k-points
    REAL(KIND = DP), ALLOCATABLE :: wkf_(:)
    !! weights k-points
    !
    ! This routine select the k-points with eigenvalues within the fsthick and
    ! then create a bztoibz mapping of those points and their symmetry operation s_bztoibz
    ! xkf_fst and wkf_fst are allocated inside
    CALL kpoint_grid_fst(nsym, s, t_rev, nrr_k, dims, &
                         irvec_k, ndegen_k, nkf1, nkf2, nkf3, nkqtotf, nelec)
    !
    ALLOCATE(xkf_(3, 2 * nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fst', 'Error allocating xkf_', 1)
    ALLOCATE(wkf_(2 * nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fst', 'Error allocating wkf_', 1)
    !
    xkf_(:,:) = zero
    DO ik = 1, nkqtotf
      ikk = 2 * ik - 1
      ikq = ikk + 1
      xkf_(:,ikk)   = xkf_fst(:, ik)
      xkf_(:,ikq)   = xkf_fst(:, ik)
      wkf_(ikk)   = 2.d0 * wkf_fst(ik)
      wkf_(ikq)   = 0.d0
    ENDDO
    DEALLOCATE(xkf_fst, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fst', 'Error deallocating wkf_fst', 1)
    DEALLOCATE(wkf_fst, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fst', 'Error deallocating wkf_fst', 1)
    !
    IF (iterative_bte) THEN
      ! Fold the points in the region [0-1] from the region -0.5,0.5
      DO ik = 1, 2 * nkqtotf
        DO idir = 1, 3
          IF (xkf_(idir, ik) < 0.0d0 ) THEN
            xkf_(idir, ik) = xkf_(idir, ik) + 1.0d0
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    ! redefine nkqtotf to include the k+q points
    !
    nkqtotf = 2 * nkqtotf
    !
#if defined(__MPI)
    CALL mp_bcast(nkqtotf, ionode_id, inter_pool_comm)
    !
    ! scatter the k points of the fine mesh across the pools
    !
    nkqf = 2 * (nkqtotf / (2 * npool))
    rest = (nkqtotf - nkqf * npool) / 2
    IF (my_pool_id < rest) THEN
      nkqf = nkqf + 2
      lower_bnd = my_pool_id * nkqf + 1
      upper_bnd = lower_bnd + nkqf - 1
    ELSE
      lower_bnd = rest * (nkqf + 2) + (my_pool_id - rest) * nkqf + 1
      upper_bnd = lower_bnd + nkqf - 1
    ENDIF
    !
    nkf = nkqf / 2
    CALL mp_bcast(xkf_, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkf_, ionode_id, inter_pool_comm)
    !
#else
    !
    nkqf = nkqtotf
    nkf = nkqf / 2
    lower_bnd = 1
    upper_bnd = nkqf
    !
#endif
    !
    ! Assign the weights and vectors to the correct bounds
    !
    ALLOCATE(xkf(3, nkqf), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fst', 'Error allocating xkf', 1)
    ALLOCATE(wkf(nkqf), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fst', 'Error allocating wkf', 1)
    xkf(:,:) = xkf_(:, lower_bnd:upper_bnd)
    !
    IF (noncolin) THEN
      wkf(:) = wkf_(lower_bnd:upper_bnd) / 2.d0
    ELSE
      wkf(:) = wkf_(lower_bnd:upper_bnd)
    ENDIF
    !
    WRITE(stdout, '(5x,"Size of k point mesh for interpolation: ",i10)' ) nkqtotf
    WRITE(stdout, '(5x,"Max number of k points per pool:",7x,i10)' ) nkqf
    !
    DEALLOCATE(xkf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fst', 'Error deallocating xkf_', 1)
    DEALLOCATE(wkf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadkmesh_fst', 'Error deallocating wkf_', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE loadkmesh_fst
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kpoint_grid_epw(n_sym, time_reversal, s, t_rev, nkc1, nkc2, nkc3, n_irr)
    !-----------------------------------------------------------------------
    !!
    !!  Automatic generation of a uniform grid of k-points with symmetry.
    !!  Routine copied from PW/src/kpoint_grid.f90.
    !!  We had to duplicate because the bztoibz array was deallocated and is needed in  EPW
    !!  This routine is sequential. For parallelized routine, see kpoint_grid_fst
    !!
    USE kinds,            ONLY : DP
    USE division,         ONLY : fkbounds
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE constants_epw,    ONLY : eps6, zero
    USE elph2,            ONLY : bztoibz, s_bztoibz, xkf_irr, wkf_irr
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_INTEGER
#endif
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: time_reversal
    !! True if time reversal
    INTEGER, INTENT(in) :: n_sym
    !! Number of symmetry
    INTEGER, INTENT(in) :: nkc1, nkc2, nkc3
    !! k-grid
    INTEGER, INTENT(in) :: t_rev(48)
    !! Time-reversal sym
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! Symmetry matrice.
    INTEGER, INTENT(inout) :: n_irr
    !! Number of IBZ points
    !
    ! Local variables
    LOGICAL :: in_the_list
    !! Is the current point in the list
    INTEGER(KIND = 2) :: s_save(nkc1 * nkc2 * nkc3)
    !! Temporary symmetry matrix
    INTEGER :: nkr
    !! Total number of points
    INTEGER :: i, j, k
    !! Index on grid size
    INTEGER(KIND = 2) :: ns
    !! Index on symmetry operations
    INTEGER :: n
    !! Global k-point index
    INTEGER :: nk
    !! Equivalent point
    INTEGER :: equiv(nkc1 * nkc2 * nkc3)
    !! Equivalent k-points
    INTEGER :: ik
    !! K-point index
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: fact
    !! Normalization factor
    REAL(KIND = DP) :: xkr(3)
    !! Current point
    REAL(KIND = DP) :: xx, yy, zz
    !! Current point coordinate
    REAL(KIND = DP):: xkg(3, nkc1 * nkc2 * nkc3)
    !! Current point
    REAL(KIND = DP) :: wkk(nkc1 * nkc2 * nkc3)
    !! Weight of the k-point
    !
    nkr = nkc1 * nkc2 * nkc3
    ALLOCATE(bztoibz(nkr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error allocating bztoibz', 1)
    !! Symmetry matrix that links an point to its IBZ friend.
    ALLOCATE(s_bztoibz(nkr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error allocating s_bztoibz', 1)
    !
    equiv(:) = 0
    s_save(:) = 0
    !
    DO i = 1, nkc1
      DO j = 1, nkc2
        DO k = 1, nkc3
          ! this is nothing but consecutive ordering
          n = (k - 1) + ( j- 1 ) * nkc3 + (i - 1) * nkc2 * nkc3 + 1
          ! xkg are the components of the complete grid in crystal axis
          xkg(1, n) = DBLE(i - 1) / nkc1
          xkg(2, n) = DBLE(j - 1) / nkc2
          xkg(3, n) = DBLE(k - 1) / nkc3
        ENDDO
      ENDDO
    ENDDO
    !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
    !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
    DO nk = 1, nkr
      equiv(nk) = nk
    ENDDO
    !
    DO nk = 1, nkr
      !  check if this k-point has already been found equivalent to another
      IF (equiv(nk) == nk) THEN
        wkk(nk) = 1.0d0
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
        DO ns = 1, n_sym
          DO i = 1, 3
            xkr(i) = s(i, 1, ns) * xkg(1, nk) &
                   + s(i, 2, ns) * xkg(2, nk) &
                   + s(i, 3, ns) * xkg(3, nk)
            xkr(i) = xkr(i) - NINT(xkr(i))
          ENDDO
          IF(t_rev(ns) == 1) xkr = -xkr
          xx = xkr(1) * nkc1
          yy = xkr(2) * nkc2
          zz = xkr(3) * nkc3
          in_the_list = ABS(xx - NINT(xx)) <= eps6 .AND. &
                        ABS(yy - NINT(yy)) <= eps6 .AND. &
                        ABS(zz - NINT(zz)) <= eps6
          IF (in_the_list) THEN
            i = MOD(NINT(xkr(1) * nkc1 + 2 * nkc1), nkc1) + 1
            j = MOD(NINT(xkr(2) * nkc2 + 2 * nkc2), nkc2) + 1
            k = MOD(NINT(xkr(3) * nkc3 + 2 * nkc3), nkc3) + 1
            n = (k - 1) + (j - 1) * nkc3 + (i - 1) * nkc2 * nkc3 + 1
            IF (n > nk .AND. equiv(n) == n) THEN
              equiv(n) = nk
              wkk(nk) = wkk(nk) + 1.0d0
              s_save(n) = ns
            ELSE
              IF (equiv(n) /= nk .OR. n < nk) CALL errore('kpoint_grid_epw', &
                 'something wrong in the checking algorithm', 1)
            ENDIF
          ENDIF
           ! We comment this out because n_sym == 48 and not 24.
           !IF (time_reversal) THEN
           !   xx =-xkr(1) * nkc1
           !   yy =-xkr(2) * nkc2
           !   zz =-xkr(3) * nkc3
           !   in_the_list = ABS(xx - NINT(xx)) <= eps6 .AND. &
           !                 ABS(yy - NINT(yy)) <= eps6 .AND. &
           !                 ABS(zz - NINT(zz)) <= eps6
           !   IF (in_the_list) THEN
           !     i = MOD(NINT(-xkr(1) * nkc1  + 2 * nkc1), nkc1) + 1
           !     j = MOD(NINT(-xkr(2) * nkc2  + 2 * nkc2), nkc2) + 1
           !     k = MOD(NINT(-xkr(3) * nkc3  + 2 * nkc3), nkc3) + 1
           !     n = (k - 1) + (j - 1) * nkc3 + (i - 1) * nkc2 * nkc3 + 1
           !     IF (n > nk .AND. equiv(n) == n) THEN
           !       equiv(n) = nk
           !       wkk(nk) = wkk(nk) + 1.0d0
           !       s_save(n) = -ns
           !     ELSE
           !       IF (equiv(n) /= nk .OR. n < nk) CALL errore('kpoint_grid', &
           !       'something wrong in the checking algorithm', 2)
           !     ENDIF
           !   ENDIF
           !ENDIF
        ENDDO
      ENDIF
    ENDDO
    !
    n_irr = 0
    ! Now do the symmetry mapping.
    DO nk = 1, nkr
      bztoibz(nk) = equiv(nk)
    ENDDO
    DO nk = 1, nkr
      ! If its an irreducible point
      IF (equiv(nk) == nk) THEN
        n_irr = n_irr + 1
        bztoibz(nk) = n_irr
        DO ik = nk, nkr
          IF (equiv(ik) == nk) THEN
            bztoibz(ik) = n_irr
          ENDIF
        ENDDO
        ! Then you have the identity matrix
        s_bztoibz(nk) = 1
      ELSE
        s_bztoibz(nk) = s_save(nk)
      ENDIF
    ENDDO
    !
    ! Create the IBZ k-point list and weights
    ! Note: xkf_irr is in crystal coordinate.
    ALLOCATE(xkf_irr(3, n_irr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error allocating xkf_irr', 1)
    ALLOCATE(wkf_irr(n_irr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error allocating wkf_irr', 1)
    xkf_irr(:, :) = zero
    wkf_irr(:)    = zero
    n_irr = 0
    fact = zero
    DO nk = 1, nkr
      IF (equiv(nk) == nk) THEN
        n_irr = n_irr + 1
        wkf_irr(n_irr) = wkk(nk)
        fact = fact + wkf_irr(n_irr)
        ! bring back into to the first BZ
        xkf_irr(:, n_irr) = xkg(:, nk) - NINT(xkg(:, nk))
      ENDIF
    ENDDO
    !
    ! Normalize weights to one
    DO nk = 1, n_irr
      wkf_irr(nk) = wkf_irr(nk) / fact
    ENDDO
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE kpoint_grid_epw
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kpoint_grid_fst(n_sym, s, t_rev, nrr_k, dims, &
                              irvec_k, ndegen_k, nkf1, nkf2, nkf3, n_irr, nelec)
    !-----------------------------------------------------------------------
    !!
    !!  Automatic generation of a uniform fine grid of k-points in the IBZ
    !!  parallelized over k-points.
    !!  Only points within the fsthick are kept.
    !!  bztoibz and s_bztoibz are allocated and computed here with dimension
    !!  nkpt_bzfst = number of point in the full BZ within the fsthick.
    !!
    USE kinds,            ONLY : DP
    USE division,         ONLY : fkbounds
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_world,         ONLY : mpime, world_comm
    USE mp_global,        ONLY : my_pool_id, npool
    USE io_global,        ONLY : stdout
    USE epwcom,           ONLY : fsthick, fermi_energy, nbndsub, scissor
    USE constants_epw,    ONLY : zero, twopi, ci, eps6, eps2, ryd2ev, czero
    USE elph2,            ONLY : chw, wkf_fst, xkf_fst, s_bztoibz, bztoibz, map_fst, &
                                 nkpt_bzfst, nbndskip
    USE wan2bloch,        ONLY : hamwan2bloch
    USE wigner,           ONLY : wigner_seitz_wrap, backtoWS
    USE noncollin_module, ONLY : noncolin
    USE constants_epw,    ONLY : one, two, eps8
# if defined(__MPI)
    USE parallel_include, ONLY : MPI_INTEGER, MPI_SUM, MPI_IN_PLACE, MPI_INTEGER2
# endif
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: n_sym
    !! Number of Bravais symmetry
    INTEGER, INTENT(in) :: nkf1, nkf2, nkf3
    !! Fine k-point grid
    INTEGER, INTENT(in) :: t_rev(48)
    !! Time-reversal symmetry operation
    INTEGER, INTENT(in) :: s(3,3,48)
    !! Symmetry matrix of the crystal
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Coordinates of real space vector for electrons
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: nelec
    !! Number of electrons
    INTEGER, INTENT(out) :: n_irr
    !! Number of IBZ k-points
    !
    ! Local variables
    LOGICAL :: in_the_list
    !! .true. if point is in the list
    LOGICAL :: found
    !! Is the reconstructed k-point part of the original set.
    LOGICAL :: low_core
    !! Are you the lowest rank core with that symmetric k-point
    INTEGER :: nkft
    !! Total number of fine k-points
    INTEGER :: lower_bnd
    !! Lower bound for the k-depend index among the mpi pools
    INTEGER :: upper_bnd
    !! Upper bound for the k-depend index among the mpi pools
    INTEGER :: i, j, k
    !! Index of x,y,z k-points
    INTEGER :: nkpt
    !! Number of k-points per core
    INTEGER :: nkpt_tmp
    !! Number of k-points up to that core
    INTEGER :: nk_inside(npool)
    !! Number of k-points inside the fsthick window
    INTEGER :: ik, jk
    !! K-point index
    INTEGER :: n
    !! Id of the point
    INTEGER(KIND = 2) :: ns
    !! Counter on symmetries
    INTEGER :: iw
    !! Counter on WS dimension
    INTEGER :: iw2
    !! Counter on WS dimension
    INTEGER :: ir
    !! Counter on real-space variable
    INTEGER :: icbm
    !! Index of the CBM
    INTEGER :: pos
    !! Position of the minimum in a vector
    INTEGER :: ierr
    !! Error status
    INTEGER :: nb
    !! Rotation index
    INTEGER :: n_check
    !! Number of full BZ points within the strickt fst
    INTEGER int2type
    !! 2 byte integer type MPI
    INTEGER :: ks(n_sym)
    !! Position of k-point equal by symmetry on the full BZ
    INTEGER :: ks_in(n_sym)
    !! Position of k-point equal by symmetry within the fsthick.
    INTEGER :: val(n_sym)
    !! Minimal value of the equivalent k-point
    INTEGER, ALLOCATABLE :: equiv(:)
    !! k-point equivalence to find IBZ per core
    INTEGER, ALLOCATABLE :: equiv_loc(:)
    !! Local equiv on the full grid of k-points
    INTEGER, ALLOCATABLE :: map_tmp(:)
    !! Temporary map per core inside fsthick
    INTEGER, ALLOCATABLE :: map_para(:)
    !! map of the full BZ homogeneous grid
    INTEGER, ALLOCATABLE :: wkf_in(:)
    !! Global k-point weights of the full BZ inside [fsthick * 1.1] per core
    INTEGER(KIND = 2), ALLOCATABLE :: s_save(:)
    !! Save the rotation index
    REAL(KIND = DP) :: etf(nbndsub)
    !! Eigen-energies for a given k-point
    REAL(KIND = DP) :: xkr(3)
    !! Rotated current k-point
    REAL(KIND = DP) :: xx, yy, zz
    !! Current k-points
    REAL(KIND = DP) :: rdotk(nrr_k)
    !! $r\cdot k$
    REAL(KIND = DP) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP) :: nelec_aux
    !! Temporary nelec, used if etf_mem == 3
    REAL(KIND = DP) :: xkf_rot(3)
    !! Current k-point coordinate rotated with symmetry
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP), ALLOCATABLE :: xkf_para(:,:)
    !! part of the full BZ homogeneous grid
    REAL(KIND = DP), ALLOCATABLE :: xkf_tmp(:,:)
    !! Temporary k-point per core inside fsthick
    REAL(KIND = DP), ALLOCATABLE :: xkf_in(:,:)
    !! Global k-point coordinate of the full BZ inside [fsthick * 1.1]
    COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP) :: cfac(nrr_k, dims, dims)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    !
    ! Total number of k-points in the full BZ
    nkft = nkf1 * nkf2 * nkf3
    nk_inside(:) = 0
    cfac(:, :, :) = czero
    !
    ! Split the total points among cores
    CALL fkbounds(nkft, lower_bnd, upper_bnd)
    nkpt = upper_bnd - lower_bnd + 1
    !
    ! 1) First we find all the points within the fsthick in the full BZ
    ALLOCATE(xkf_para(3, nkpt), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating xkf_para', 1)
    ALLOCATE(map_para(nkpt), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating map_para', 1)
    ALLOCATE(xkf_tmp(3, nkpt), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating xkf_tmp', 1)
    ALLOCATE(map_tmp(nkpt), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating map_tmp', 1)
    xkf_para(:,:) = zero
    xkf_tmp(:,:)  = zero
    map_para(:)   = 0
    map_tmp(:)    = 0
    !
    ! Create a parallelized k-point grids containing all the points in the full BZ.
    DO i = 1, nkf1
      DO j = 1, nkf2
        DO k = 1, nkf3
          ! this is nothing but consecutive ordering
          n = (k - 1) + (j - 1) * nkf3 + (i - 1) * nkf2 * nkf3 + 1
          IF ((n >= lower_bnd) .AND. (n <= upper_bnd)) THEN
            !  xkg are the components of the complete grid in crystal axis
            xkf_para(1, n - lower_bnd + 1) = DBLE(i - 1) / nkf1
            xkf_para(2, n - lower_bnd + 1) = DBLE(j - 1) / nkf2
            xkf_para(3, n - lower_bnd + 1) = DBLE(k - 1) / nkf3
            map_para(n - lower_bnd + 1)    = n
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    ! Compute Hamiltonian : Wannier -> Bloch
    ! We select the k-points for which the eigenenergy is within the fsthick
    IF (nbndskip > 0) THEN
      IF (noncolin) THEN
        nelec_aux = nelec - one * nbndskip
      ELSE
        nelec_aux = nelec - two * nbndskip
      ENDIF
    ELSE
      nelec_aux = nelec
    ENDIF
    !
    icbm = 1
    IF (ABS(scissor) > eps6) THEN
      IF (noncolin) THEN
        icbm = FLOOR(nelec_aux / 1.0d0) + 1
      ELSE
        icbm = FLOOR(nelec_aux / 2.0d0) + 1
      ENDIF
    ENDIF
    ! This is simply because dgemv take only real number (not integer)
    irvec_r = REAL(irvec_k, KIND = DP)
    DO ik = 1, nkpt
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkf_para(:, ik), 1, 0.0_DP, rdotk, 1 )
      DO iw = 1, dims
        DO iw2 = 1, dims
          DO ir = 1, nrr_k
            IF (ndegen_k(ir, iw2, iw) > 0) cfac(ir, iw2, iw) = EXP(ci * rdotk(ir)) / ndegen_k(ir, iw2, iw)
          ENDDO
        ENDDO
      ENDDO
      CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:), chw, cfac, dims)
      !
      ! Apply scissor shift
      etf(icbm:nbndsub) = etf(icbm:nbndsub) + scissor
      !
      ! We take a slightly bigger fsthick as some point do not fully respect crystal symmetry.
      IF (MINVAL(ABS(etf(:) - fermi_energy)) < fsthick * 1.2) THEN
        nk_inside(my_pool_id + 1)             = nk_inside(my_pool_id + 1) + 1
        xkf_tmp(:, nk_inside(my_pool_id + 1)) = xkf_para(:, ik)
        map_tmp(nk_inside(my_pool_id + 1))    = map_para(ik)
      ENDIF
    ENDDO ! ik
    !
    CALL mp_sum(nk_inside, world_comm)
    !
    ! Total number of points inside the fsthick
    nkpt_bzfst = SUM(nk_inside)
    !
    WRITE(stdout, '(5x,a,i9)') 'Number of k-points inside fsthick * 1.2 in the full BZ: ', nkpt_bzfst
    !
    ! Total k-point array with all the kpoints inside fsthick
    ALLOCATE(xkf_in(3, nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating xkf_in', 1)
    ALLOCATE(wkf_in(nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating wkf_in', 1)
    ALLOCATE(map_fst(nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating map_fst', 1)
    ALLOCATE(s_save(nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating s_save', 1)
    xkf_in(:, :) = zero
    wkf_in(:)    = 0
    map_fst(:)   = 0
    s_save(:)    = 0
    !
    ! Number of points up to the current core
    IF (my_pool_id == 0) THEN
      nkpt_tmp = 0
    ELSE
      nkpt_tmp = SUM(nk_inside(1:my_pool_id))
    ENDIF
    !
    ! We create a global list containg all the k-points inside the fsthick (xfk_in)
    ! as well as a map between the full BZ k-grid and the reduced inside fsthick.
    ! Every cores fill it in parallel
    DO ik = 1, nk_inside(my_pool_id + 1)
      xkf_in(:, nkpt_tmp + ik) = xkf_tmp(:, ik)
      map_fst(nkpt_tmp + ik)   = map_tmp(ik)
    ENDDO ! ik
    !
    ! Now merge everything accross cores
    CALL mp_sum(xkf_in, world_comm)
    CALL mp_sum(map_fst, world_comm)
    !
    DEALLOCATE(xkf_para, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating xkf_para', 1)
    DEALLOCATE(xkf_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating xkf_tmp', 1)
    DEALLOCATE(map_para, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating map_para', 1)
    DEALLOCATE(map_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating map_tmp', 1)
    !
    ! 2) We reduce those points to the IBZ using symmetry
    !
    ! equiv(ik) =ik : k-point ik is not equivalent to any previous k-point
    ! equiv(ik)!=ik : k-point ik is equivalent to k-point equiv(ik)
    ALLOCATE(equiv(nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating equiv', 1)
    ALLOCATE(equiv_loc(nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating equiv_loc', 1)
    equiv(:) = 0
    equiv_loc(:) = 0
    !
    ! Parallelize the the k-point inside the fsthick
    CALL fkbounds(nkpt_bzfst, lower_bnd, upper_bnd)
    !
    ! Number of k-points on that core
    nkpt = upper_bnd - lower_bnd + 1
    !
    DO ik = 1, nkpt
      equiv(ik + lower_bnd - 1) = ik + lower_bnd - 1
    ENDDO
    DO ik = 1, nkpt_bzfst
      equiv_loc(ik) = ik
    ENDDO
    !
    DO ik = 1, nkpt
      ! Check if this k-point has already been found equivalent to another
      IF (equiv(ik + lower_bnd - 1) == ik + lower_bnd - 1) THEN
        wkf_in(ik + lower_bnd - 1) = 1
        ! Position of the k-points equivalent by symmetry to the current kpoint
        ks(:) = 0
        ks_in(:) = 0
        DO ns = 1, n_sym
          DO i = 1, 3
            xkr(i) = s(i, 1, ns) * xkf_in(1, ik + lower_bnd - 1) &
                   + s(i, 2, ns) * xkf_in(2, ik + lower_bnd - 1) &
                   + s(i, 3, ns) * xkf_in(3, ik + lower_bnd - 1)
            xkr(i) = xkr(i) - NINT(xkr(i))
          ENDDO
          IF(t_rev(ns) == 1) xkr = -xkr
          xx = xkr(1) * nkf1
          yy = xkr(2) * nkf2
          zz = xkr(3) * nkf3
          in_the_list = ABS(xx - NINT(xx)) <= eps6 .AND. &
                        ABS(yy - NINT(yy)) <= eps6 .AND. &
                        ABS(zz - NINT(zz)) <= eps6
          IF (in_the_list) THEN
            i = MOD(NINT(xkr(1) * nkf1 + 2 * nkf1), nkf1) + 1
            j = MOD(NINT(xkr(2) * nkf2 + 2 * nkf2), nkf2) + 1
            k = MOD(NINT(xkr(3) * nkf3 + 2 * nkf3), nkf3) + 1
            n = (k - 1) + (j - 1) * nkf3 + (i - 1) * nkf2 * nkf3 + 1
            !
            pos = MINLOC(ABS(map_fst - n), 1)
            !
            ks(ns) = n ! Position in the full BZ
            ks_in(ns) = pos ! Position in the nkpt_bzfst subset.
            val(ns) = ABS(map_fst(pos) - n) ! If val is not 0, this means the point is not within fsthick
            !
          ENDIF ! in_the_list
        ENDDO ! n_sym
        !
        low_core = .TRUE.
        DO ns = 1, n_sym
          ! Not the lowest core with that set of equiv. k-points ==> nullify that current k-point position
          ! Note: There is a specific case where we need to keep the point.
          !       If the current k-point has symmetric friend that are outside
          !       the scope of the current core but also outside the fsthick. We
          !       need to have found the point with val(ns) == 0
          IF (ks_in(ns) < lower_bnd .AND. val(ns) == 0) THEN
            equiv(ik + lower_bnd - 1) = 0
            wkf_in(ik + lower_bnd - 1) = 0
            low_core = .FALSE.
            EXIT ! exit the loop
          ENDIF
        ENDDO
        !
        ! If you are the lowest core
        IF (low_core) THEN
          DO ns = 1, n_sym
            IF (ks(ns) > map_fst(ik + lower_bnd - 1) .AND. equiv_loc(ks_in(ns)) == ks_in(ns) .AND. val(ns) == 0) THEN
              equiv_loc(ks_in(ns)) = ik + lower_bnd - 1
              equiv(ks_in(ns)) = ik + lower_bnd - 1
              s_save(ks_in(ns)) = ns
              wkf_in(ik + lower_bnd - 1) = wkf_in(ik + lower_bnd - 1) + 1
            ENDIF
          ENDDO
        ENDIF
        !
      ENDIF ! equiv
      !
    ENDDO ! ik
    !
    CALL mp_sum(equiv, world_comm)
    CALL mp_sum(wkf_in, world_comm)
# if defined(__MPI)
    !CALL MPI_TYPE_CREATE_F90_INTEGER(SIK2, int2type, ierr)
    !CALL MPI_ALLreduce(MPI_IN_PLACE, s_save, nkpt_bzfst, int2type, MPI_SUM, world_comm, ierr)
    CALL MPI_ALLreduce(MPI_IN_PLACE, s_save, nkpt_bzfst, MPI_INTEGER2, MPI_SUM, world_comm, ierr)
#endif
    !
    DEALLOCATE(equiv_loc, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating equiv_loc', 1)
    !
    ! Check that the sum of the weights correctly give the total number of point inside fsthick in the full BZ
    IF (nkpt_bzfst /= SUM(wkf_in)) THEN
      WRITE(stdout,'(5x,a,i9)') 'Reconstituded number of points inside the fsthick in the full BZ from weights ', SUM(wkf_in)
      CALL errore('kpoint_grid_fst', 'The weights do not sum correctly to the number of points.', 1)
    ENDIF
    !
    ALLOCATE(bztoibz(nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating bztoibz', 1)
    ALLOCATE(s_bztoibz(nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating s_bztoibz', 1)
    ALLOCATE(map_tmp(nkpt_bzfst), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating map_tmp', 1)
    bztoibz(:) = 0
    s_bztoibz(:) = 0
    map_tmp(:) = 0
    !
    DO ik = 1, nkpt_bzfst
      bztoibz(ik) = equiv(ik)
    ENDDO
    !
    ! Number of IBZ point within the fsthick * 1.2
    n_irr = 0
    DO ik = 1, nkpt_bzfst
      IF (equiv(ik) == ik) THEN
        n_irr = n_irr + 1
        bztoibz(ik) = n_irr
        DO jk = ik, nkpt_bzfst
          IF (equiv(jk) == ik) THEN
            bztoibz(jk) = n_irr
          ENDIF
        ENDDO ! jk
        map_tmp(n_irr) = ik
      ENDIF ! equiv(ik) == ik
    ENDDO
    !
    ! Now do the symmetry mapping.
    DO ik = 1, nkpt_bzfst
      ! If its an irreducible point
      IF (equiv(ik) == ik) THEN
        ! Then you have the identity matrix
        s_bztoibz(ik) = 1
      ELSE
        s_bztoibz(ik) = s_save(ik)
      ENDIF
    ENDDO
    !
    ! 3) Find irreducible k points and weights
    !
    ALLOCATE(xkf_fst(3, n_irr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating xkf_fst', 1)
    ALLOCATE(wkf_fst(n_irr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating wkf_fst', 1)
    xkf_fst(:, :) = zero
    wkf_fst(:) = zero
    !
    DO ik = 1, n_irr
      xkf_fst(:, ik) = xkf_in(:, map_tmp(ik))
      wkf_fst(ik)    = REAL(wkf_in(map_tmp(ik)), KIND = DP)
    ENDDO
    wkf_fst(:) = wkf_fst(:) / (nkf1 * nkf2 * nkf3)
    !
    DEALLOCATE(wkf_in, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating wkf_in', 1)
    DEALLOCATE(equiv, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating equiv', 1)
    DEALLOCATE(map_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating map_tmp', 1)
    !
    !DBSP
    !DO ik = 1, n_irr
    !  print*,'ik xkf_fst wkf ',ik,  xkf_fst(:, ik), wkf_fst(ik), map_fst(ik)
    !ENDDO
    !
    ! 4) Check that fsthick * 1.2 was enough to take all the symmetry equivalent points
    !
    ! First we take only the IBZ points that are within the strick fsthick for cheking
    n_check = 0
    ALLOCATE(xkf_tmp(3, n_irr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error allocating xkf_tmp', 1)
    xkf_tmp(:, :) = zero
    !
    DO ik = 1, n_irr
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkf_fst(:, ik), 1, 0.0_DP, rdotk, 1 )
      DO iw = 1, dims
        DO iw2 = 1, dims
          DO ir = 1, nrr_k
            IF (ndegen_k(ir, iw2, iw) > 0) cfac(ir, iw2, iw) = EXP(ci * rdotk(ir)) / ndegen_k(ir, iw2, iw)
          ENDDO
        ENDDO
      ENDDO
      CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:), chw, cfac, dims)
      !
      ! Apply scissor shift
      etf(icbm:nbndsub) = etf(icbm:nbndsub) + scissor
      !
      IF (MINVAL(ABS(etf(:) - fermi_energy)) < fsthick) THEN
        n_check = n_check + 1
        xkf_tmp(:, n_check) = xkf_fst(:, ik)
      ENDIF
    ENDDO
    !
    ! Split the total n_check points strictly inside fsthick among cores
    CALL fkbounds(n_check, lower_bnd, upper_bnd)
    nkpt = upper_bnd - lower_bnd + 1
    !
    ! Use symmetries to reconstruct the BZ from IBZ and check that all points were in xkf_in
    DO ik = 1, nkpt
      DO nb = 1, n_sym
        ! Note that s is in crystal
        sa(:, :) = DBLE(s(:, :, nb))
        xkf_rot = MATMUL(sa, xkf_tmp(:, ik + lower_bnd - 1))
        !
        ! Check that the point xkf_rot is part of the orginal xkf_in
        found = .FALSE.
        DO jk = 1, nkpt_bzfst
          IF ((ABS(xkf_rot(1) - xkf_in(1, jk) - NINT(xkf_rot(1) - xkf_in(1, jk))) < eps8) .AND. &
              (ABS(xkf_rot(2) - xkf_in(2, jk) - NINT(xkf_rot(2) - xkf_in(2, jk))) < eps8) .AND. &
              (ABS(xkf_rot(3) - xkf_in(3, jk) - NINT(xkf_rot(3) - xkf_in(3, jk))) < eps8)) THEN
             found = .TRUE.
             EXIT
          ENDIF
        ENDDO
        !
        IF (found .eqv. .FALSE.) CALL errore('kpoint_grid_fst', 'K-point not found. Increase fsthick factor 1.2', 1)
      ENDDO ! nb
    ENDDO ! ik
    !
    DEALLOCATE(xkf_in, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating xkf_in', 1)
    DEALLOCATE(xkf_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_fst', 'Error deallocating xkf_in', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kpoint_grid_fst
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE loadqmesh_para()
    !-----------------------------------------------------------------------
    !!
    !!  Load fine q mesh and distribute among pools
    !!
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
    USE mp,        ONLY : mp_bcast, mp_sum
    USE mp_world,  ONLY : mpime
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE epwcom,    ONLY : filqf, nqf1, nqf2, nqf3, plselfen, specfun_pl, &
                          rand_q, rand_nq, mp_mesh_q, system_2d, lscreen
    USE elph2,     ONLY : xqf, wqf, nqf, nqtotf
    USE cell_base, ONLY : at, bg
    USE symm_base, ONLY : s, t_rev, time_reversal, nsym
    USE io_var,    ONLY : iunqf
    USE noncollin_module, ONLY : noncolin
    USE constants_epw, ONLY : eps4
    USE low_lvl,   ONLY : init_random_seed
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 10) :: coordinate_type
    !! filkf coordinate type (crystal or cartesian)
    LOGICAL, EXTERNAL  :: imatches
    !! Regex matching text.
    INTEGER :: iq
    !! Q-point index
    INTEGER :: lower_bnd
    !! Lower-bound of the core
    INTEGER :: upper_bnd
    !! Lower-bound of the core
    INTEGER :: i, j, k
    !! Directions
    INTEGER :: ios
    !! Status of the reading of the file
    INTEGER :: rest
    !! Remaining of cores numbers
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP), ALLOCATABLE :: xqf_(:, :)
    !! Temporary q-point
    REAL(KIND = DP), ALLOCATABLE :: wqf_(:)
    !! Temporary weight of q-point
    !
    IF (mpime == ionode_id) THEN
      IF (filqf /= '') THEN ! load from file
        !
        WRITE(stdout, *) '    Using q-mesh file: ', TRIM(filqf)
        IF (lscreen) WRITE(stdout, *) '     WARNING: if lscreen=.TRUE., q-mesh needs to be [-0.5:0.5] (crystal)'
        OPEN(UNIT = iunqf, FILE = filqf, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('loadkmesh_para', 'Opening file ' // filqf, ABS(ios))
        READ(iunqf, *) nqtotf, coordinate_type
        IF (TRIM(coordinate_type) .EQ. ' ') coordinate_type = 'crystal'
        IF (.NOT. imatches("crystal", coordinate_type) .AND. .NOT. imatches("cartesian", coordinate_type)) &
          CALL errore('loadqmesh_para', 'ERROR: Specify either crystal or cartesian coordinates in the filqf file', 1)
        !
        ALLOCATE(xqf_(3, nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf_', 1)
        ALLOCATE(wqf_(nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating wqf_', 1)
        !
        DO iq = 1, nqtotf
          !
          READ(iunqf, *) xqf_(:, iq), wqf_(iq)
          !
        ENDDO
        CLOSE(iunqf)
        IF (imatches("cartesian", coordinate_type)) THEN
          CALL cryst_to_cart(nqtotf, xqf_, at, -1)
        ENDIF
        !
      ELSEIF ((nqf1 /= 0) .AND. (nqf2 /= 0) .AND. (nqf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_q) THEN
          IF (lscreen) CALL errore('loadqmesh', 'If lscreen = .TRUE. do not use mp_mesh_q', 1)
          ! get size of the mp_mesh in the irr wedge
          WRITE(stdout, '(a,3i4)') '     Using uniform MP q-mesh: ', nqf1, nqf2, nqf3
          !
          ALLOCATE(xqf_(3, nqf1 * nqf2 * nqf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf_ ', 1)
          ALLOCATE(wqf_(nqf1 * nqf2 * nqf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating wqf_', 1)
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nsym, time_reversal, .FALSE., s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf_, wqf_)
          DEALLOCATE(xqf_, wqf_, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error deallocating xqf_, wqf_', 1)
          ALLOCATE(xqf_ (3, nqtotf), wqf_(nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf_ (3, nqtotf), wqf_', 1)
          CALL kpoint_grid(nsym, time_reversal, .FALSE., s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf_, wqf_)
          !
          ! bring the k point to crystal coordinates
          CALL cryst_to_cart(nqtotf, xqf_, at, -1)
          !
        ELSE
          !
          WRITE (stdout, '(a,3i5)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
          !
          nqtotf =  nqf1 * nqf2 * nqf3
          ALLOCATE(xqf_(3, nqtotf), wqf_(nqtotf) , STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf_(3, nqtotf), wqf_', 1)
          wqf_(:) = 0.d0
          DO iq = 1, nqf1 * nqf2 * nqf3
            wqf_(iq) = 1.d0 / (DBLE(nqtotf))
          ENDDO
          DO i = 1, nqf1
            DO j = 1, nqf2
              DO k = 1, nqf3
                iq = (i - 1) * nqf2 * nqf3 + (j - 1) * nqf3 + k
                xqf_(1, iq) = DBLE(i - 1) / DBLE(nqf1)
                xqf_(2, iq) = DBLE(j - 1) / DBLE(nqf2)
                xqf_(3, iq) = DBLE(k - 1) / DBLE(nqf3)
              ENDDO
            ENDDO
          ENDDO
          IF (lscreen .OR. specfun_pl .OR. plselfen) xqf_(:, :) = xqf_(:, :) - 0.5d0
          !
        ENDIF
        !
      ELSEIF (rand_q) THEN  ! random points
        !
        WRITE(stdout, *) '    Using random q-mesh: ', rand_nq
        !
        nqtotf = rand_nq
        ALLOCATE(xqf_ (3, nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf_ ', 1)
        ALLOCATE(wqf_(nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating wqf_', 1)
        !
        CALL init_random_seed()
        !
        DO iq = 1, nqtotf
          !
          !
          wqf_(iq) = 1.d0 / DBLE(nqtotf)
          !
          IF (system_2d == 'no') THEN
            CALL random_number(xqf_(:, iq))
            IF (lscreen .OR. specfun_pl .OR. plselfen) xqf_(:, iq) = xqf_(:, iq) - 0.5d0
          ELSE
            CALL random_number(xqf_(1:2, iq))
            IF (lscreen .OR. specfun_pl .OR. plselfen) xqf_(1:2, iq) = xqf_(1:2, iq) - 0.5d0
            xqf_(3, iq) = 0.d0
          ENDIF
          !
        ENDDO
        !
      ELSE ! don't know how to get grid
        CALL errore('loadqmesh_para', "Cannot load fine q points", 1)
      ENDIF
    ENDIF
    !
#if defined(__MPI)
    CALL mp_bcast (nqtotf, ionode_id, inter_pool_comm)
    !
    !  scatter the q points of the fine mesh across the pools
    !
    nqf = (nqtotf / npool)
    rest = (nqtotf - nqf * npool) / 2
    IF (my_pool_id < rest) THEN
      nqf = nqf + 2
      lower_bnd = my_pool_id * nqf + 1
      upper_bnd = lower_bnd + nqf - 1
    ELSE
      lower_bnd = rest * (nqf + 2) + (my_pool_id - rest) * nqf + 1
      upper_bnd = lower_bnd + nqf - 1
    ENDIF
    !
    IF (mpime /= ionode_id) THEN
      ALLOCATE(xqf_(3, nqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf_', 1)
      ALLOCATE(wqf_(nqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating wqf_', 1)
    ENDIF
    CALL mp_bcast(xqf_, ionode_id, inter_pool_comm)
    CALL mp_bcast(wqf_, ionode_id, inter_pool_comm)
    !
#else
    !
    ! In serial the definitions are much easier
    !
    nqf = nqtotf
    lower_bnd = 1
    upper_bnd = nqf
    !
#endif
    !
    !  Assign the weights and vectors to the correct bounds
    !
    ALLOCATE(xqf(3, nqf), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf', 1)
    ALLOCATE(wqf(nqf), STAT = ierr)
    IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating wqf', 1)
    xqf(:, :) = xqf_ (:, lower_bnd:upper_bnd)
    IF (noncolin) THEN
      wqf(:) = wqf_(lower_bnd:upper_bnd) / 2.d0
    ELSE
      wqf(:) = wqf_(lower_bnd:upper_bnd)
    ENDIF
    !
    IF (ABS(SUM(wqf_(:)) - 1.d0) > eps4) &
       WRITE(stdout,'(5x,"WARNING: q-point weigths do not add up to 1 [loadqmesh_para]")')
    !
    WRITE(stdout, '(5x,"Size of q point mesh for interpolation: ",i10)') nqtotf
    WRITE(stdout, '(5x,"Max number of q points per pool:",7x,i10)') nqf
    !
    DEALLOCATE(xqf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error deallocating xqf_', 1)
    DEALLOCATE(wqf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error deallocating wqf_', 1)
    !-----------------------------------------------------------------------
    END SUBROUTINE loadqmesh_para
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE loadqmesh_serial()
    !-----------------------------------------------------------------------
    !!
    !!  Load fine q mesh in sequential
    !!
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast
    USE mp_world,  ONLY : mpime
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE epwcom,    ONLY : filqf, nqf1, nqf2, nqf3, &
                          rand_q, rand_nq, mp_mesh_q, system_2d, lscreen, &
                          plselfen, specfun_pl, &
                          scell_mat_plrn, scell_mat, as, bs
    USE elph2,     ONLY : xqf, wqf, nqtotf, nqf
    USE cell_base, ONLY : at, bg
    USE symm_base, ONLY : s, t_rev, time_reversal, nsym
    USE io_var,    ONLY : iunqf
    USE low_lvl,   ONLY : init_random_seed, matinv3
    USE constants_epw, ONLY : eps4, eps8
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 10) :: coordinate_type
    !! filqf coordinate type (crystal or cartesian)
    LOGICAL, EXTERNAL  :: imatches
    !! Regex matching text.
    INTEGER :: iq
    !! Q-index
    INTEGER :: i, j, k
    !! Directions
    INTEGER :: ios
    !! Status integer
    INTEGER :: ierr
    !! Error status
    INTEGER :: iRp1, iRp2, iRp3, Rpmax, nRp
    !! Number of unit cells within supercell
    INTEGER :: Rp_crys_p(3)
    !! Unit cell vectors in primitive crystal coordinates
    INTEGER, ALLOCATABLE :: Rp(:, :)
    !! List of unit cell vectors within supercell in primitive crystal coords
    INTEGER :: iGs1, iGs2, iGs3, Gsmax, nGs
    !! Number of supercell G-vectors within primitive reciprocal unit cell
    INTEGER :: Gs_crys_s(3)
    !! Supercell G-vectors in supercell reciprocal coordinates
    REAL(KIND = DP) :: ap(3, 3), bp(3, 3)
    !! Auxiliary definitions of real and reciprocal primitive cell vector matrix
    REAL(KIND = DP) :: scell_mat_b(3, 3)
    !! Reciprocal lattice transformation matrix
    REAL(KIND = DP) :: p2s(3, 3), bs2p(3, 3)
    !! Transformation matrix from primitive to supercell crystal coordinates
    REAL(KIND = DP) :: Rp_crys_s(3)
    !! Unit cell vectors in supercell crystal coordinates
    REAL(KIND = DP) :: Gs_crys_p(3)
    !! Supercell G-vectors in primitive crystal coordinates
    REAL(KIND = DP), ALLOCATABLE :: Gs(:, :)
    !! Supercell G-vectors within primitive reciprocal unit cell
    !
    IF (mpime == ionode_id) THEN
      IF (filqf /= '') THEN ! load from file
        !
        ! Each pool gets its own copy from the action=read statement
        !
        WRITE(stdout, *) '    Using q-mesh file: ', TRIM(filqf)
        IF (lscreen) WRITE(stdout, *) '     WARNING: if lscreen=.TRUE., q-mesh needs to be [-0.5:0.5] (crystal)'
        OPEN(UNIT = iunqf, FILE = filqf, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('loadqmesh_serial', 'opening file ' // filqf, ABS(ios))
        READ(iunqf, *) nqtotf, coordinate_type
        IF (TRIM(coordinate_type) .EQ. ' ') coordinate_type = 'crystal'
        IF (.NOT. imatches("crystal", coordinate_type) .AND. .NOT. imatches("cartesian", coordinate_type)) &
          CALL errore('loadqmesh_serial', 'ERROR: Specify either crystal or cartesian coordinates in the filqf file', 1)
        !
        ALLOCATE(xqf(3, nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf', 1)
        ALLOCATE(wqf(nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
        !
        DO iq = 1, nqtotf
          !
          READ (iunqf, *) xqf(:, iq), wqf(iq)
          !
        ENDDO
        CLOSE(iunqf)
        IF (imatches("cartesian", coordinate_type)) THEN
          CALL cryst_to_cart(nqtotf, xqf, at, -1)
        ENDIF
        !
      !JLB
      ELSEIF (scell_mat_plrn) THEN
        !
        WRITE(stdout, '(a)') ' '
        WRITE(stdout, '(a)') '     Supercell transformation activated (q), as=S*at'
        WRITE(stdout, '(a,3i4)') '     S(1, 1:3): ', scell_mat(1, 1:3)
        WRITE(stdout, '(a,3i4)') '     S(2, 1:3): ', scell_mat(2, 1:3)
        WRITE(stdout, '(a,3i4)') '     S(3, 1:3): ', scell_mat(3, 1:3)
        !
        ap = TRANSPOSE(at)
        as = MATMUL(scell_mat,ap)
        !
        WRITE(stdout, '(a)') '     Transformed lattice vectors (alat units):'
        WRITE(stdout, '(a,3f12.6)') '     as(1, 1:3): ', as(1, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     as(2, 1:3): ', as(2, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     as(3, 1:3): ', as(3, 1:3)
        !
        scell_mat_b = matinv3(REAL(scell_mat, DP))
        scell_mat_b = TRANSPOSE(scell_mat_b)
        !
        WRITE(stdout, '(a)') '     Reciprocal lattice transformation matrix, Sbar = (S^{-1})^{t}:'
        WRITE(stdout, '(a,3f12.6)') '     Sbar(1, 1:3): ', scell_mat_b(1, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     Sbar(2, 1:3): ', scell_mat_b(2, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     Sbar(3, 1:3): ', scell_mat_b(3, 1:3)
        !
        bp = TRANSPOSE(bg)
        bs = MATMUL(scell_mat_b, bp)
        !
        WRITE(stdout, '(a)') '     Transformed reciprocal lattice vectors (2pi/alat units):'
        WRITE(stdout, '(a,3f12.6)') '     bs(1, 1:3): ', bs(1, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     bs(2, 1:3): ', bs(2, 1:3)
        WRITE(stdout, '(a,3f12.6)') '     bs(3, 1:3): ', bs(3, 1:3)
        WRITE(stdout, '(a)') '  '
        !
        ! Define transformation matrix from primitive crystal coordinates
        ! to supercell crystal coordinates Rp_crys_s = ((a_s)^{T})^{-1} (a_p)^{T} Rp_crys_p
        p2s = matinv3(TRANSPOSE(as))
        p2s = MATMUL(p2s,TRANSPOSE(ap))
        !
        ! Find how many unit cells are contained within the supercell
        Rpmax = 5*MAXVAL(scell_mat) ! This should be large enough to find all
        ALLOCATE(Rp(3, Rpmax**3), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating Rp', 1)
        Rp  = 0
        nRp = 0
        DO iRp1 = -Rpmax, Rpmax
          DO iRp2 = -Rpmax, Rpmax
            DO iRp3 = -Rpmax, Rpmax
              Rp_crys_p = (/iRp1, iRp2, iRp3/)
              Rp_crys_s = MATMUL(p2s, Rp_crys_p)
              ! Unit cell within supercell if Rp\in(0,1) in supercell crystal coordinates
              IF (ALL(Rp_crys_s > -eps8) .AND. ALL(Rp_crys_s < 1.d0-eps8)) THEN
                nRp = nRp + 1
                Rp(1:3, nRp) = Rp_crys_p
              END IF
            END DO
          END DO
        END DO
        WRITE(stdout, '(a, 3i6)') '     Number of unit cells within supercell:', nRp
        !
        IF (ALLOCATED(Rp)) DEALLOCATE(Rp)
        !
        ! Define transformation matrix from reciprocal supercell crystal coordinates
        ! to reciprocal primitive crystal coordinates
        ! Gs_crys_p = ((bp)^{T})^{-1} (bs)^{T} Gs_crys_s
        bs2p = matinv3(TRANSPOSE(bp))
        bs2p = MATMUL(bs2p, TRANSPOSE(bs))
        !
        ! Find how many q-points lie within primitive reciprocal cell
        Gsmax = Rpmax ! This should be large enough to find all
        ALLOCATE(Gs(3, Gsmax**3), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating Gs', 1)
        Gs = 0.d0
        nGs = 0
        DO iGs1 = -Gsmax, Gsmax
          DO iGs2 = -Gsmax, Gsmax
            DO iGs3 = -Gsmax, Gsmax
              Gs_crys_s = (/iGs1, iGs2, iGs3/)
              Gs_crys_p = MATMUL(bs2p, Gs_crys_s)
              ! Gs within primitive reciprocal unit cell if Gs\in(0,1) in crys_p coords.
              IF (ALL(Gs_crys_p > -eps8) .AND. ALL(Gs_crys_p < 1.d0-eps8)) THEN
                nGs = nGs + 1
                Gs(1:3, nGs) = Gs_crys_p
              END IF
            END DO
          END DO
        END DO
        WRITE(stdout, '(a, 3i6)') '     Number of q-points needed:', nGs
        !
        ! Save list of needed q-points
        nqtotf = nGs
        ALLOCATE(xqf(3, nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf', 1)
        ALLOCATE(wqf(nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
        !
        DO iq = 1, nqtotf
          !
          xqf(:, iq) = Gs(1:3, iq)
          wqf(iq) = 1.d0 ! weight not important for polaron
          !
        ENDDO
        !
        IF (ALLOCATED(Gs)) DEALLOCATE(Gs)
        !
      !JLB
      ELSEIF ((nqf1 /= 0) .AND. (nqf2 /= 0) .AND. (nqf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_q) THEN
          IF (lscreen) CALL errore ('loadqmesh', 'If lscreen=.TRUE. do not use mp_mesh_q',1)
          ! get size of the mp_mesh in the irr wedge
          WRITE (stdout, '(a,3i5)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
          !
          ALLOCATE(xqf(3, nqf1 * nqf2 * nqf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf', 1)
          ALLOCATE(wqf(nqf1 * nqf2 * nqf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nsym, time_reversal, .FALSE., s, t_rev, bg, nqf1 * nqf2 * nqf3, &
               0, 0, 0, nqf1, nqf2, nqf3, nqtotf, xqf, wqf)
          DEALLOCATE(xqf, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error deallocating xqf', 1)
          DEALLOCATE(wqf, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error deallocating wqf', 1)
          ALLOCATE(xqf(3, nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf', 1)
          ALLOCATE(wqf(nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
          CALL kpoint_grid(nsym, time_reversal, .FALSE., s, t_rev, bg, nqf1 * nqf2 * nqf3, &
               0,0,0, nqf1, nqf2, nqf3, nqtotf, xqf, wqf)
          !
          ! bring xqf in crystal coordinates
          CALL cryst_to_cart(nqtotf, xqf, at, -1)
          !
        ELSE
          ! currently no offset.
          ! q's are in crystal coordinates in xqf
          WRITE (stdout, '(a,3i5)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
          !
          nqtotf = nqf1 * nqf2 * nqf3
          ALLOCATE(xqf (3, nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf ', 1)
          ALLOCATE(wqf(nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
          wqf(:) = 1.d0 / (DBLE(nqtotf))
          DO i = 1, nqf1
            DO j = 1, nqf2
              DO k = 1, nqf3
                iq = (i - 1) * nqf2 * nqf3 + (j - 1) * nqf3 + k
                xqf(1, iq) = DBLE(i - 1) / DBLE(nqf1)
                xqf(2, iq) = DBLE(j - 1) / DBLE(nqf2)
                xqf(3, iq) = DBLE(k - 1) / DBLE(nqf3)
              ENDDO
            ENDDO
          ENDDO
          IF (lscreen .OR. specfun_pl .OR. plselfen) THEN
            xqf(:, :) = xqf(:, :) - 0.5d0
            WRITE (stdout, '(a)') '     Shifting q mesh to [-0.5:0.5['
          ENDIF
          !
        ENDIF
      ELSEIF (rand_q) THEN  ! random points
        ! random grid
        WRITE (stdout, *) '    Using random q-mesh: ', rand_nq
        !
        nqtotf = rand_nq
        ALLOCATE(xqf(3, nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf', 1)
        ALLOCATE(wqf(nqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
        wqf(:) = 1.d0 / (DBLE(nqtotf))
        !
        CALL init_random_seed()
        !
        DO iq = 1, nqtotf
          !
          IF (system_2d == 'no') THEN
            CALL random_number(xqf(:, iq))
            IF (lscreen .OR. specfun_pl .OR. plselfen) xqf(:, iq) = xqf(:, iq) - 0.5d0
          ELSE
            CALL random_number(xqf(1:2, iq))
            IF (lscreen .OR. specfun_pl .OR. plselfen) xqf(1:2, iq) = xqf(1:2, iq) - 0.5d0
            xqf(3, iq) = 0.d0
          ENDIF
          !
        ENDDO
        IF (lscreen .OR. specfun_pl .OR. plselfen) WRITE(stdout, '(a)') '    Shifting q mesh to [-0.5:0.5['
        !
      ELSE ! don't know how to get grid
        CALL errore('loadqmesh_serial', "Cannot load fine q points", 1)
      ENDIF
      !
      ! Since serial
      nqf = nqtotf
    ENDIF
    !
    CALL mp_bcast(nqf, ionode_id, inter_pool_comm)
    CALL mp_bcast(nqtotf, ionode_id, inter_pool_comm)
    IF (mpime /= ionode_id) THEN
      ALLOCATE(xqf(3, nqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf', 1)
      ALLOCATE(wqf(nqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
    ENDIF
    CALL mp_bcast(xqf, ionode_id, inter_pool_comm)
    CALL mp_bcast(wqf, ionode_id, inter_pool_comm)
    !
    IF (ABS(SUM(wqf) - 1.d0) > eps4) &
      WRITE(stdout,'(5x,"WARNING: q-point weigths do not add up to 1 [loadqmesh_serial]")')
    !
    WRITE(stdout, '(5x,"Size of q point mesh for interpolation: ",i10)' ) nqtotf
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE loadqmesh_serial
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE xqf_otf(iq, xxq)
    !-----------------------------------------------------------------------
    !
    !! This routine computes the q-point coordinate on the fly.
    !! Indeed for very large grids, having all the points in memory is a bottlneck.
    !
    !-----------------------------------------------------------------------
    USE kinds,   ONLY : DP
    USE epwcom,  ONLY : nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! qvectors index
    REAL(KIND = DP), INTENT(inout) :: xxq(3)
    !! Output qvector, in crystal coordinates
    !
    ! Local variables
    INTEGER :: iq1
    !! x-crystal coordinate
    INTEGER :: iq2
    !! y-crystal coordinate
    INTEGER :: iq3
    !! z-crystal coordinate
    !
    ! Integer division from iq = iq3 + iq2 * nqf3 + iq1 * nkqf2 * nkqf3 + 1
    iq1 = (iq - 1) / (nqf2 * nqf3)
    iq2 = ((iq - 1) / nqf3) - iq1 * nqf2
    iq3 = (iq - 1) - iq1 * nqf2 * nqf3 - iq2 * nqf3
    !
    xxq(1) = REAL(iq1, KIND = DP) / REAL(nqf1, KIND = DP)
    xxq(2) = REAL(iq2, KIND = DP) / REAL(nqf2, KIND = DP)
    xxq(3) = REAL(iq3, KIND = DP) / REAL(nqf3, KIND = DP)
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE xqf_otf
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE xkf_otf(ik, xxk)
    !-----------------------------------------------------------------------
    !
    !! This routine computes the k-point coordinate on the fly.
    !! Indeed for very large grids, having all the points in memory is a bottlneck.
    !
    !-----------------------------------------------------------------------
    USE kinds,   ONLY : DP
    USE epwcom,  ONLY : nkf1, nkf2, nkf3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! qvectors index
    REAL(KIND = DP), INTENT(inout) :: xxk(3)
    !! Output qvector, in crystal coordinates
    !
    ! Local variables
    INTEGER :: ik1
    !! x-crystal coordinate
    INTEGER :: ik2
    !! y-crystal coordinate
    INTEGER :: ik3
    !! z-crystal coordinate
    !
    ! Integer division from ik = ik3 + ik2 * nkf3 + ik1 * nkkf2 * nkkf3 + 1
    ik1 = (ik - 1) / (nkf2 * nkf3)
    ik2 = ((ik - 1) / nkf3) - ik1 * nkf2
    ik3 = (ik - 1) - ik1 * nkf2 * nkf3 - ik2 * nkf3
    !
    xxk(1) = REAL(ik1, KIND = DP) / REAL(nkf1, KIND = DP)
    xxk(2) = REAL(ik2, KIND = DP) / REAL(nkf2, KIND = DP)
    xxk(3) = REAL(ik3, KIND = DP) / REAL(nkf3, KIND = DP)
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE xkf_otf
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE qwindow_wrap(totq, nrr_k, dims, ndegen_k, irvec_r, cufkk, cufkq)
    !-----------------------------------------------------------------------
    !!
    !! Author: S. Ponc\'e
    !! This routine determines which q-points falls within the fsthick windows
    !! Store the result in the selecq.fmt file
    !! If the file exists, automatically restart from the file
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, plrn, &
                              scell_mat_plrn, nbndsub, selecqread
    USE elph2,         ONLY : selecq, nqf
    USE mp_world,      ONLY : mpime, world_comm
    USE io_global,     ONLY : ionode_id, stdout
    USE mp,            ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: totq
    !! Total number of q-points computed
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(in) :: dims
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(out) :: cufkk(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), INTENT(out) :: cufkq(nbndsub, nbndsub)
    !! the same, for points k+q
    !
    ! Local
    LOGICAL :: homogeneous
    !! Check if the grids are homogeneous and commensurate
    LOGICAL :: exst
    !! If the file exist
    INTEGER :: ierr
    !! Error status
    INTEGER :: iq
    !! Counter on fine q-point grid
    !
    !
    homogeneous = .FALSE.
    IF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0) .AND. &
        (nqf1 /= 0) .AND. (nqf2 /= 0) .AND. (nqf3 /= 0)) THEN
      IF ((MOD(nkf1, nqf1) == 0) .AND. (MOD(nkf2, nqf2) == 0) .AND. (MOD(nkf3, nqf3) == 0)) THEN
        homogeneous = .TRUE.
      ENDIF
    ELSE
      homogeneous = .FALSE.
    ENDIF
    !
    totq = 0
    !
    IF (plrn .OR. scell_mat_plrn) THEN
      ! For polaron calculations, all the q points have to be included
      totq = nqf
      ALLOCATE(selecq(nqf), STAT = ierr)
      IF (ierr /= 0) CALL errore('qwindow_wrap', 'Error allocating selecq', 1)
      DO iq = 1, nqf
        selecq(iq) = iq
      ENDDO
      !
    ELSE
      ! Check if the file has been pre-computed
      IF (mpime == ionode_id) THEN
        INQUIRE(FILE = 'selecq.fmt', EXIST = exst)
      ENDIF
      CALL mp_bcast(exst, ionode_id, world_comm)
      !
      IF (exst) THEN
        IF (selecqread) THEN
          WRITE(stdout, '(5x,a)')' '
          WRITE(stdout, '(5x,a)')'Reading selecq.fmt file. '
          CALL qwindow(exst, nrr_k, dims, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
        ELSE
          WRITE(stdout, '(5x,a)')' '
          WRITE(stdout, '(5x,a)')'A selecq.fmt file was found but re-created because selecqread == .FALSE. '
          CALL qwindow(.FALSE., nrr_k, dims, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
        ENDIF
      ELSE ! exst
        IF (selecqread) THEN
          CALL errore( 'qwindow_wrap', 'Variable selecqread == .TRUE. but file selecq.fmt not found.',1 )
        ELSE
          CALL qwindow(exst, nrr_k, dims, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
        ENDIF
      ENDIF
      !
      WRITE(stdout, '(5x,a,i8,a)')'We only need to compute ', totq, ' q-points'
      WRITE(stdout, '(5x,a)')' '
      !
    ENDIF ! plrn .OR. scell_mat_plrn
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE qwindow_wrap
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE qwindow(exst, nrr_k, dims, irvec_r, ndegen_k, &
                       cufkk, cufkq, homogeneous)
    !-----------------------------------------------------------------------
    !!
    !! Author: S. Ponc\'e
    !! This routine pre-computes the q-points that falls within the fstichk.
    !! If at least 1 k-point is such that at least one k+q eigenenergy falls
    !! within the user-defined fstichk, then the q-point is taken.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : nqf, xqf, xkf, chw, nkf, nqtotf, &
                              map_rebal, nktotf, bztoibz, map_fst, totq, &
                              selecq
    USE io_global,     ONLY : ionode_id, stdout
    USE io_var,        ONLY : iunselecq
    USE mp_global,     ONLY : npool, world_comm, my_pool_id
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_sum, mp_bcast
    USE constants_epw, ONLY : twopi, ci, zero, eps6, ryd2ev, czero
    USE epwcom,        ONLY : nbndsub, fsthick, use_ws, mp_mesh_k, nkf1, nkf2, &
                              nkf3, iterative_bte, restart_step, scissor,      &
                              ephwrite, etf_mem
    USE noncollin_module, ONLY : noncolin
    USE pwcom,         ONLY : ef, nelec
    USE wan2bloch,     ONLY : hamwan2bloch
    USE poolgathering, ONLY : poolgather
    USE low_lvl,       ONLY : create_interval, bisection
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: exst
    !! If the file exist
    LOGICAL, INTENT(in) :: homogeneous
    !! Check if the grids are homogeneous and commensurate
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(in) :: dims
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(out) :: cufkk(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), INTENT(out) :: cufkq(nbndsub, nbndsub)
    !! the same, for points k+q
    !
    ! Local variable
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: iq
    !! Counter on coarse q-point grid
    INTEGER :: ik, ikk, ikl
    !! Counter on coarse k-point grid
    INTEGER :: icbm
    !! Index for the CBM
    INTEGER :: ibnd
    !! Band index
    INTEGER :: found(npool)
    !! Indicate if a q-point was found within the window
    INTEGER :: iw
    !! Counter on bands when use_ws == .TRUE.
    INTEGER :: iw2
    !! Counter on bands when use_ws == .TRUE.
    INTEGER :: ir
    !! Counter for WS loop
    INTEGER :: nqtot
    !! Total number of q-point for verifications
    INTEGER :: ind1
    !! Index of the k point from the full grid.
    INTEGER :: ind2
    !! Index of the k+q point from the full grid.
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: ierr
    !! Error status
    INTEGER :: nkloc
    !! number of k-point selected on that cpu
    INTEGER :: kmap(nkf)
    !! k-point that are selected for that cpu
    INTEGER :: n_intval
    !! Number of intervals
    INTEGER, ALLOCATABLE :: bztoibz_tmp(:)
    !! Temporary mapping when etf_mem < 3
    INTEGER, ALLOCATABLE :: selecq_tmp(:)
    !! Temporary list of selected q-points
    INTEGER, ALLOCATABLE :: val_intval(:)
    !! Value of the first element of each intervals
    INTEGER, ALLOCATABLE :: pos_intval(:)
    !! Position of the first element of each intervals
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    REAL(KIND = DP) :: xkk(3)
    !! Current k-point on the fine grid
    REAL(KIND = DP) :: xkq(3)
    !! Current k-point on the fine grid
    REAL(KIND = DP) :: rdotk(nrr_k)
    !! $r\cdot k$
    REAL(KIND = DP) :: rdotk2(nrr_k)
    !! $r\cdot k$
    REAL(KIND = DP) :: etf_loc(nbndsub, nkf)
    !! Eigen-energies all full k-grid.
    REAL(KIND = DP) :: etf_locq(nbndsub, nkf)
    !! Eigen-energies all full k-grid.
    REAL(KIND = DP) :: etf_all(nbndsub, nktotf)
    !! Eigen-energies all full k-grid.
    REAL(KIND = DP) :: etf_tmp(nbndsub)
    !! Temporary Eigen-energies at a give k-point
    COMPLEX(KIND = DP) :: cfac(nrr_k, dims, dims)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    COMPLEX(KIND = DP) :: cfacq(nrr_k, dims, dims)
    !! Used to store $e^{2\pi r \cdot k+q}$ exponential
    !
    rdotk(:)  = 0
    rdotk2(:) = 0
    cfac(:, :, :)  = czero
    cfacq(:, :, :) = czero
    !
    IF (exst) THEN
      IF (mpime == ionode_id) THEN
        OPEN(UNIT = iunselecq, FILE = 'selecq.fmt', STATUS = 'old', IOSTAT = ios)
        READ(iunselecq, *) totq
        ALLOCATE(selecq(totq), STAT = ierr)
        IF (ierr /= 0) CALL errore('qwindow', 'Error allocating selecq', 1)
        selecq(:) = 0
        READ(iunselecq, *) nqtot
        READ(iunselecq, *) selecq(:)
        CLOSE(iunselecq)
      ENDIF
      CALL mp_bcast(totq, ionode_id, world_comm)
      IF (mpime /= ionode_id) THEN
        ALLOCATE(selecq(totq), STAT = ierr)
        IF (ierr /= 0) CALL errore('qwindow', 'Error allocating selecq', 1)
      ENDIF
      CALL mp_bcast(nqtot , ionode_id, world_comm)
      CALL mp_bcast(selecq, ionode_id, world_comm)
      IF (nqtot /= nqtotf) THEN
        CALL errore('qwindow', 'Cannot read from selecq.fmt, the q-point grid or &
          & fsthick window are different from read one. Remove the selecq.fmt file and restart.', 1 )
      ENDIF
      !
      IF (homogeneous) THEN
        ! In case of k-point symmetry
        IF (mp_mesh_k .AND. etf_mem < 3) THEN
          IF (iterative_bte .OR. ephwrite) THEN
            ALLOCATE(bztoibz_tmp(nkf1 * nkf2 * nkf3), STAT = ierr)
            IF (ierr /= 0) CALL errore('qwindow', 'Error allocating bztoibz_tmp', 1)
            bztoibz_tmp(:) = 0
            DO ikbz = 1, nkf1 * nkf2 * nkf3
              bztoibz_tmp(ikbz) = map_rebal(bztoibz(ikbz))
            ENDDO
            bztoibz(:) = bztoibz_tmp(:)
            DEALLOCATE(bztoibz_tmp, STAT = ierr)
            IF (ierr /= 0) CALL errore('qwindow', 'Error deallocating bztoibz_tmp', 1)
          ENDIF
        ENDIF ! mp_mesh_k
      ENDIF ! homogeneous
      !
    ELSE
      ALLOCATE(selecq_tmp(nqf), STAT = ierr)
      IF (ierr /= 0) CALL errore('qwindow', 'Error allocating selecq_tmp', 1)
      selecq_tmp(:) = 0
      etf_loc(:, :)  = zero
      etf_locq(:, :) = zero
      etf_all(:, :) = zero
      !
      IF (homogeneous) THEN
        ! First store eigen energies on full grid.
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          xkk = xkf(:, ikk)
          CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
          IF (use_ws) THEN
            DO iw = 1, dims
              DO iw2 = 1, dims
                DO ir = 1, nrr_k
                  IF (ndegen_k(ir, iw2, iw) > 0) THEN
                    cfac(ir, iw2, iw)  = EXP(ci * rdotk(ir)) / ndegen_k(ir, iw2, iw)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ELSE
            cfac(:, 1, 1) = EXP(ci * rdotk(:)) / ndegen_k(:, 1, 1)
          ENDIF
          CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf_loc(:, ik), chw, cfac, dims)
        ENDDO
        CALL poolgather(nbndsub, nktotf, nkf, etf_loc, etf_all)
        !
        ! In case of k-point symmetry
        IF (mp_mesh_k .AND. etf_mem < 3) THEN
          IF (iterative_bte .OR. ephwrite) THEN
            ALLOCATE(bztoibz_tmp(nkf1 * nkf2 * nkf3), STAT = ierr)
            IF (ierr /= 0) CALL errore('qwindow', 'Error allocating bztoibz_tmp', 1)
            bztoibz_tmp(:) = 0
            DO ikbz = 1, nkf1 * nkf2 * nkf3
              bztoibz_tmp(ikbz) = map_rebal(bztoibz(ikbz))
            ENDDO
            bztoibz(:) = bztoibz_tmp(:)
            DEALLOCATE(bztoibz_tmp, STAT = ierr)
            IF (ierr /= 0) CALL errore('qwindow', 'Error deallocating bztoibz_tmp', 1)
          ENDIF
        ENDIF ! mp_mesh_k
        !
        ! Apply a scissor shift to CBM if required by user
        ! The shift is apply to k and k+q
        IF (ABS(scissor) > eps6) THEN
          IF (noncolin) THEN
            icbm = FLOOR(nelec / 1.0d0) + 1
          ELSE
            icbm = FLOOR(nelec / 2.0d0) + 1
          ENDIF
          !
          DO ik = 1, nktotf
            DO ibnd = icbm, nbndsub
              etf_all(ibnd, ik) = etf_all(ibnd, ik) + scissor
            ENDDO
          ENDDO
        ENDIF
        !
        ! Note 1: To find if a k+q point is within the fsthick we need to obtain the mapping
        !         between the index of the point within the fsthick and the index of the point
        !         within the full BZ. This is most efficiently done with bisection.
        ! Note 2: When the number of points within the fshtick window is large, the bissection
        !         is slow. One can speed this up by doing a pre-search since the map_fst is
        !         composed of monotonically increasing numbers (ordered list).
        IF (etf_mem == 3) THEN
          ! We divide map_fst into n_intval intervals
          n_intval = NINT(SQRT(REAL(SIZE(map_fst, 1), KIND = DP)))
          ALLOCATE(val_intval(n_intval), STAT = ierr)
          IF (ierr /= 0) CALL errore('qwindow', 'Error allocating val_intval', 1)
          ALLOCATE(pos_intval(n_intval), STAT = ierr)
          IF (ierr /= 0) CALL errore('qwindow', 'Error allocating pos_intval', 1)
          ! We select 1 point every n_interval
          CALL create_interval(SIZE(map_fst, 1), map_fst, n_intval, val_intval, pos_intval)
        ENDIF
        !
        DO iq = 1, nqf
          ! xqf has to be in crystal coordinate
          IF (etf_mem == 3) THEN
            ! The q-point coordinate is generate on the fly for each q-point
            CALL xqf_otf(iq, xxq)
          ELSE
            xxq = xqf(:, iq)
          ENDIF
          !
          found(:) = 0
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            xkk = xkf(:, ikk)
            xkq = xkk + xxq
            !
            CALL kpmq_map(xkk, (/0d0, 0d0, 0d0/), 1, ind1)
            CALL kpmq_map(xkk, xxq, 1, ind2)
            IF (ind1 == 0 .OR. ind2 == 0) CALL errore ('qwindow', 'ind1 or ind2 cannot be 0', 1)
            !
            IF (etf_mem == 3) THEN
              ! Bisection method to find the index on the grid of the points inside fsthick
              ! from the index on the full BZ grid.
              CALL bisection(SIZE(map_fst, 1), map_fst, ind1, n_intval, val_intval, pos_intval)
              IF (ind1 == 0) CYCLE
              CALL bisection(SIZE(map_fst, 1), map_fst, ind2, n_intval, val_intval, pos_intval)
              IF (ind2 == 0) CYCLE
            ENDIF
            !
            ! Use k-point symmetry
            IF (mp_mesh_k) THEN
              IF ((MINVAL(ABS(etf_all(:, bztoibz(ind1)) - ef)) < fsthick) .AND. &
                  (MINVAL(ABS(etf_all(:, bztoibz(ind2)) - ef)) < fsthick)) THEN
                found(my_pool_id + 1) = 1
                EXIT ! exit the loop
              ENDIF
            ELSE
              IF ((MINVAL(ABS(etf_all(:, ind1) - ef)) < fsthick) .AND. &
                  (MINVAL(ABS(etf_all(:, ind2) - ef)) < fsthick)) THEN
                found(my_pool_id + 1) = 1
                EXIT ! exit the loop
              ENDIF
            ENDIF
            !
          ENDDO ! k-loop
          ! If found on any k-point from the pools
          CALL mp_sum(found, world_comm)
          !
          IF (SUM(found) > 0) THEN
            totq = totq + 1
            selecq_tmp(totq) = iq
            !
            IF (MOD(totq, restart_step) == 0) THEN
              WRITE(stdout, '(5x,a,i15,i15)')'Number selected, total', totq, iq
            ENDIF
          ENDIF
        ENDDO ! iq
      ELSE ! homogeneous
        !
        ! Apply a scissor shift to CBM if required by user
        ! The shift is apply to k and k+q
        IF (ABS(scissor) > eps6) THEN
          IF (noncolin) THEN
            icbm = FLOOR(nelec / 1.0d0) + 1
          ELSE
            icbm = FLOOR(nelec / 2.0d0) + 1
          ENDIF
        ENDIF
        !
        ! First compute the k-points eigenenergies for efficiency reasons
        nkloc = 0
        kmap(:) = 0
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          xkk = xkf(:, ikk)
          CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
          IF (use_ws) THEN
            DO iw = 1, dims
              DO iw2 = 1, dims
                DO ir = 1, nrr_k
                  IF (ndegen_k(ir, iw2, iw) > 0) THEN
                    cfac(ir, iw2, iw)  = EXP(ci * rdotk(ir)) / ndegen_k(ir, iw2, iw)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ELSE
            cfac(:, 1, 1)  = EXP(ci * rdotk(:)) / ndegen_k(:, 1, 1)
          ENDIF
          CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf_tmp(:), chw, cfac, dims)
          !
          IF (ABS(scissor) > eps6) THEN
            DO ibnd = icbm, nbndsub
              etf_tmp(ibnd) = etf_tmp(ibnd) + scissor
            ENDDO
          ENDIF
          ! Check for the k-points in this pool
          IF (MINVAL(ABS(etf_tmp(:) - ef)) < fsthick) THEN
            nkloc = nkloc + 1
            kmap(nkloc) = ik
          ENDIF
        ENDDO ! k-points
        !
        ! Now compute the q-loop doing WS separately for efficiency
        IF (use_ws) THEN
          DO iq = 1, nqf
            IF (etf_mem == 3) THEN
              ! The q-point coordinate is generate on the fly for each q-point
              CALL xqf_otf(iq, xxq)
            ELSE
              xxq = xqf(:, iq)
            ENDIF
            etf_tmp(:) = zero
            found(:) = 0
            DO ikl = 1, nkloc
              ik = kmap(ikl)
              ikk = 2 * ik - 1
              xkk = xkf(:, ikk)
              xkq = xkk + xxq
              CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1)
              DO iw = 1, dims
                DO iw2 = 1, dims
                  DO ir = 1, nrr_k
                    IF (ndegen_k(ir, iw2, iw) > 0) THEN
                      cfacq(ir, iw2, iw) = EXP(ci * rdotk(ir)) / ndegen_k(ir, iw2, iw)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
              CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf_tmp(:), chw, cfacq, dims)
              !
              IF (ABS(scissor) > eps6) THEN
                DO ibnd = icbm, nbndsub
                  etf_tmp(ibnd) = etf_tmp(ibnd) + scissor
                ENDDO
              ENDIF
              !
              IF (MINVAL(ABS(etf_tmp(:) - ef)) < fsthick) THEN
                found(my_pool_id + 1) = 1
                EXIT ! exit the loop
              ENDIF
            ENDDO ! ik
            ! If found on any k-point from the pools
            CALL mp_sum(found, world_comm)
            IF (SUM(found(:)) > 0) THEN
              totq = totq + 1
              selecq_tmp(totq) = iq
              IF (MOD(totq, restart_step) == 0) THEN
                WRITE(stdout, '(5x,a,i12,i12)') 'Number selected, total', totq, iq
              ENDIF
            ENDIF
          ENDDO ! iq
        ELSE ! use_ws
          DO iq = 1, nqf
            IF (etf_mem == 3) THEN
              ! The q-point coordinate is generate on the fly for each q-point
              CALL xqf_otf(iq, xxq)
            ELSE
              xxq = xqf(:, iq)
            ENDIF
            etf_tmp(:) = zero
            found(:) = 0
            DO ikl = 1, nkloc
              ik = kmap(ikl)
              ikk = 2 * ik - 1
              xkk = xkf(:, ikk)
              xkq = xkk + xxq
              CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1)
              cfacq(:, 1, 1) = EXP(ci * rdotk(:)) / ndegen_k(:, 1, 1)
              CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf_tmp(:), chw, cfacq, dims)
              !
              IF (ABS(scissor) > eps6) THEN
                DO ibnd = icbm, nbndsub
                  etf_tmp(ibnd) = etf_tmp(ibnd) + scissor
                ENDDO
              ENDIF
              !
              IF (MINVAL(ABS(etf_tmp(:) - ef) ) < fsthick) THEN
                found(my_pool_id + 1) = 1
                EXIT ! exit the loop
              ENDIF
            ENDDO ! ik
            ! If found on any k-point from the pools
            CALL mp_sum(found, world_comm)
            IF (SUM(found(:)) > 0) THEN
              totq = totq + 1
              selecq_tmp(totq) = iq
              IF (MOD(totq, restart_step) == 0) THEN
                WRITE(stdout, '(5x,a,i12,i12)')'Number selected, total', totq, iq
              ENDIF
            ENDIF
          ENDDO ! iq
        ENDIF ! use_ws
      ENDIF ! homogeneous
      !
      ALLOCATE(selecq(totq), STAT = ierr)
      IF (ierr /= 0) CALL errore('qwindow', 'Error allocating selecq', 1)
      selecq(:) = selecq_tmp(1:totq)
      DEALLOCATE(selecq_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('qwindow', 'Error deallocating selecq_tmp', 1)
      !
      IF (mpime == ionode_id) THEN
        OPEN(UNIT = iunselecq, FILE = 'selecq.fmt', ACTION = 'write')
        WRITE(iunselecq, *) totq    ! Selected number of q-points
        WRITE(iunselecq, *) nqtotf  ! Total number of q-points
        WRITE(iunselecq, *) selecq(1:totq)
        CLOSE(iunselecq)
      ENDIF
      !
    ENDIF ! exst
    !-----------------------------------------------------------------------
    END SUBROUTINE qwindow
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE load_rebal()
    !-----------------------------------------------------------------------
    !!
    !! Author: S. Ponc\'e
    !! Routine used to rebalance the load on k-points.
    !! At the moment this routines is only called in the case of IBTE
    !! using k-point symmetry and an homogeneous grid.
    !! The k-point that are within the fshick are equally splitted among cores.
    !!
    !-----------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : etf, nkf, nkqtotf, xkf, wkf, etf, map_rebal, map_rebal_inv, &
                              lower_bnd, nktotf
    USE epwcom,        ONLY : fsthick, nbndsub, mp_mesh_k
    USE pwcom,         ONLY : ef
    USE mp_global,     ONLY : my_pool_id, npool
    USE constants_epw, ONLY : zero
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm, mp_bcast
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    INTEGER :: ipool
    !! Pool loop
    INTEGER :: ik
    !! Total k-point index
    INTEGER :: ikk
    !! K-point every 2 index
    INTEGER :: ikq
    !! Q-point every 2 index
    INTEGER :: ikpt
    !! Primary k-point index
    INTEGER :: ikpt2
    !! Secondary k-point index
    INTEGER :: rest
    !! Rest of the points
    INTEGER :: tot
    !! Total number of k-point (quotient)
    INTEGER :: counter
    !! temp variable
    INTEGER :: ierr
    !! Error status
    INTEGER :: kpt_in(nkqtotf)
    !! K-points that are within the fshick windows
    INTEGER :: kpt_out(nkqtotf)
    !! K-points that are outside of the fshick windows
    INTEGER :: map_rebal_inv_tmp(nktotf)
    !! Temporary inverse map between the initial ordering of k-point and the rebalanced one
    REAL(KIND = DP) :: xkf_all(3, nkqtotf)
    !! Collect k-point coordinate (and k+q) from all pools in parallel case
    REAL(KIND = DP) :: wkf_all(nkqtotf)
    !! Collect k-point coordinate (and k+q) from all pools in parallel case
    REAL(KIND = DP) :: etf_all(nbndsub, nkqtotf)
    !! Collect all the eigenenergies
    !
    ! Gather all the k-point coordinate and all the weights from all the pools
    ! Gather also all the eigeneneriges
    xkf_all(:, :) = zero
    wkf_all(:) = zero
    etf_all(:, :) = zero
    !
#if defined(__MPI)
    CALL poolgather2(1, nkqtotf, 2 * nkf, wkf, wkf_all )
    CALL poolgather2(3, nkqtotf, 2 * nkf, xkf, xkf_all )
    CALL poolgather2(nbndsub, nkqtotf, 2 * nkf, etf, etf_all )
#else
    xkf_all = xkf
    wkf_all = wkf
    etf_all = etf
#endif
    !
    ALLOCATE(map_rebal(nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('load_rebal', 'Error allocating map_rebal', 1)
    ALLOCATE(map_rebal_inv(nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('load_rebal', 'Error allocating map_rebal_inv', 1)
    !
    kpt_in(:) = 0
    kpt_out(:) = 0
    !
    !The sorting is done by master only
    IF (my_pool_id == 0) THEN
      !
      ikpt = 0
      ikpt2 = 0
      !
      DO ik = 1, nktotf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        IF (MINVAL(ABS(etf_all(:, ikk) - ef)) < fsthick) THEN
          ikpt = ikpt + 1
          kpt_in(ikpt) = ikk
        ELSE
          ikpt2 = ikpt2 + 1
          kpt_out(ikpt2) = ikk
        ENDIF ! fsthick
        !
      ENDDO ! ik
      ! In case of IBZ create a map of old to new points
      IF (mp_mesh_k) THEN
        map_rebal(:) = 0
        map_rebal_inv(:) = 0
        DO ik = 1, ikpt
          map_rebal((kpt_in(ik) + 1) / 2) = ik
          map_rebal_inv(ik) = (kpt_in(ik) + 1) / 2
        ENDDO
        DO ik = 1, ikpt2
          map_rebal((kpt_out(ik) + 1) /2) =  ikpt + ik
          map_rebal_inv(ikpt + ik) = (kpt_out(ik) + 1) /2
        ENDDO
      ENDIF
    ENDIF ! my_pool_id
    !
    CALL mp_bcast(kpt_in, ionode_id, inter_pool_comm)
    CALL mp_bcast(kpt_out, ionode_id, inter_pool_comm)
    CALL mp_bcast(ikpt, ionode_id, inter_pool_comm)
    CALL mp_bcast(ikpt2, ionode_id, inter_pool_comm)
    CALL mp_bcast(map_rebal, ionode_id, inter_pool_comm)
    CALL mp_bcast(map_rebal_inv, ionode_id, inter_pool_comm)
    !
    ! map_rebal contains an array with all the k-point in the window and then
    ! all the k-points out of the window.
    ! We then split those k-points such that the first core has the first k-point,
    ! the second core has the second k-point etc
    !
    tot = (nkqtotf / (2 * npool))         ! quotient
    rest = (nktotf - tot * npool)         ! reminder
    counter = 0
    DO ipool = 1, npool
      DO ik = 1,  tot
        counter = counter + 1
        map_rebal_inv_tmp(counter) = map_rebal_inv(npool * ik - (npool - ipool))
      ENDDO
      !Do the rest
      IF (ipool <= rest) THEN
        counter = counter + 1
        map_rebal_inv_tmp(counter) = map_rebal_inv(npool * (tot + 1) - (npool - ipool))
      ENDIF
    ENDDO
    map_rebal_inv(:) = map_rebal_inv_tmp(:)
    !
    ! Now recontruct map_rebal so that it is the inverse mapping as map_rebal_inv
    DO ik = 1, nktotf
      map_rebal(map_rebal_inv(ik)) = ik
    ENDDO
    !
    ! We then assign this new order to the local k-points and weights on each cpu
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      xkf(:, ikk) = xkf_all(:, 2 * map_rebal_inv(ik + lower_bnd - 1) - 1)
      xkf(:, ikk + 1) = xkf_all(:, 2 * map_rebal_inv(ik + lower_bnd - 1))
      wkf(ikk) = wkf_all(2 * map_rebal_inv(ik + lower_bnd - 1) - 1)
      wkf(ikk + 1) = wkf_all(2 * map_rebal_inv(ik + lower_bnd - 1))
    ENDDO
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE load_rebal
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE symm_mapping(nind, bztoibz_mat, xkf_all, sparse_q, sparse_k)
    !-----------------------------------------------------------------------
    !!
    !! For a given k-point in the IBZ gives the k-point index
    !! of all the k-point in the full BZ that are connected to the current
    !! one by symmetry. nsym + TR is the max number of symmetry
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE symm_base,     ONLY : nsym
    USE elph2,         ONLY : bztoibz, nktotf, ixkqf_tr, s_bztoibz_full, xqf, &
                              nkpt_bzfst, map_fst, s_bztoibz, map_rebal
    USE epwcom,        ONLY : etf_mem, nkf1, nkf2, nkf3, epmatkqread
    USE low_lvl,       ONLY : create_interval, bisection
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = 8), INTENT(in) :: nind
    !! Total number of elements per cpu
    INTEGER, INTENT(inout) :: bztoibz_mat(nsym, nktotf)
    !! For a given k-point in the IBZ gives gives the index of all the kpt in BZ connected by symmetry
    INTEGER, INTENT(in) :: sparse_q(nind)
    !! Q-point mapping index
    INTEGER, INTENT(in) :: sparse_k(nind)
    !! K-point mapping index
    REAL(KIND = DP), INTENT(in) :: xkf_all(3, nktotf)
    !! All the k-points (just k-points, not k and k+q)
    !
    ! Local variables
    INTEGER :: ik
    !! K-point index
    INTEGER :: iq
    !! Q-point index
    INTEGER(KIND = 8) :: ind
    !! Index for sparse matrix
    INTEGER :: ikbz
    !! Index on the full BZ
    INTEGER :: nkq_abs
    !! Index of the k+q point from the full grid.
    INTEGER :: n_intval
    !! Number of intervals
    INTEGER :: ierr
    !! Error index
    INTEGER :: n_sym(nktotf)
    !! Temporary matrix used to count how many symmetry for that k-point
    INTEGER :: bztoibz_tmp(nkf1 * nkf2 * nkf3)
    !! Temporary mapping
    INTEGER, ALLOCATABLE :: val_intval(:)
    !! Value of the first element of each intervals
    INTEGER, ALLOCATABLE :: pos_intval(:)
    !! Position of the first element of each intervals
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    !
    n_sym(:) = 0
    !
    IF (etf_mem == 3) THEN
      !
      DO ikbz = 1, nkpt_bzfst
        ik = bztoibz(ikbz)
        n_sym(ik) = n_sym(ik) + 1
        bztoibz_mat(n_sym(ik), ik) = ikbz
      ENDDO
      !
      ! We divide map_fst into n_intval intervals
      n_intval = NINT(SQRT(REAL(SIZE(map_fst, 1), KIND = DP)))
      ALLOCATE(val_intval(n_intval), STAT = ierr)
      IF (ierr /= 0) CALL errore('symm_mapping', 'Error allocating val_intval', 1)
      ALLOCATE(pos_intval(n_intval), STAT = ierr)
      IF (ierr /= 0) CALL errore('symm_mapping', 'Error allocating pos_intval', 1)
      ! We select 1 point every n_interval
      CALL create_interval(SIZE(map_fst, 1), map_fst, n_intval, val_intval, pos_intval)
      DO ind = 1, nind
        iq = sparse_q(ind)
        ik = sparse_k(ind)
        ! The q-point coordinate is generate on the fly for each q-point
        CALL xqf_otf(iq, xxq)
        !
        CALL kpmq_map(xkf_all(:, ik), xxq, +1, nkq_abs)
        !
        CALL bisection(SIZE(map_fst, 1), map_fst, nkq_abs, n_intval, val_intval, pos_intval)
        ! k + q cannot fall outside the points inside fsthick
        IF (nkq_abs == 0) CALL errore('ibte', 'Error in mapping the vectors', 1)
        !
        s_bztoibz_full(ind) = s_bztoibz(nkq_abs)
        ixkqf_tr(ind) = bztoibz(nkq_abs)
      ENDDO
    ELSE
      ! This call is required because for a epmatkqread restart because then
      ! qwindow is not called and therefore the map_rebal is not applied
      IF (epmatkqread) THEN
        bztoibz_tmp(:) = 0
        DO ikbz = 1, nkf1 * nkf2 * nkf3
          bztoibz_tmp(ikbz) = map_rebal(bztoibz(ikbz))
        ENDDO
        bztoibz(:) = bztoibz_tmp(:)
      ENDIF ! epmatkqread
      !
      ! Now create the mapping matrix
      DO ikbz = 1, nkf1 * nkf2 * nkf3
        ik = bztoibz(ikbz)
        n_sym(ik) = n_sym(ik) + 1
        bztoibz_mat(n_sym(ik), ik) = ikbz
      ENDDO
      !
      DO ind = 1, nind
        iq = sparse_q(ind)
        ik = sparse_k(ind)
        !
        CALL kpmq_map(xkf_all(:, ik), xqf(:, iq), +1, nkq_abs)
        s_bztoibz_full(ind) = s_bztoibz(nkq_abs)
        ixkqf_tr(ind) = bztoibz(nkq_abs)
      ENDDO
    ENDIF ! etf_mem == 3
    !
    WRITE(stdout, '(5x,"Symmetry mapping finished")')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE symm_mapping
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE special_points(nb_sp, xkf_all, xkf_sp)
    !-----------------------------------------------------------------------
    !!
    !! This routine determines  the special k-points that are sent to
    !! themselves as a result of symmetry.
    !! e.g. the point [1 1 1] is sent to itself by the symmetry that exchanges
    !!      the x and y coordinates.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE symm_base,     ONLY : s, nsym
    USE elph2,         ONLY : nkf, nktotf
    USE constants_epw, ONLY : eps6, zero
    USE wigner,        ONLY : backtoWS
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : world_comm
    USE division,      ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: nb_sp
    !! Number of special points
    INTEGER, INTENT(out), ALLOCATABLE :: xkf_sp(:, :)
    !! List of special k-points. The first index is the kpt index and the other
    REAL(KIND = DP), INTENT(in) :: xkf_all(3, nktotf)
    !! All the k-points (just k-points, not k and k+q)
    !
    ! Local variables
    LOGICAL :: sym_found
    !! Logical for IF statement
    INTEGER :: ik
    !! K-point variable
    INTEGER :: nb
    !! Symmetry index
    INTEGER :: counter
    !! Counter on the number of symmetries
    INTEGER :: nrws
    !! Maximum number of WS vectors
    INTEGER :: lower_bnd
    !! Lower bounds index after k para
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: n
    !! Loop index
    INTEGER :: m
    !! Loop index
    INTEGER :: l
    !! Loop index
    INTEGER :: counter_n
    !! Counter for special points on the border
    INTEGER :: ierr
    !! Error status
    INTEGER :: xkt_sp(48, nktotf)
    !! Temp list of special k-points
    INTEGER, PARAMETER :: nrwsx = 200
    !! Variable for WS folding
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix (intermediate step)
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in cartesian coordinate
    REAL(KIND = DP) :: xk(3)
    !! Current k-point coordinate
    REAL(KIND = DP) :: s_xk(3)
    !! Rotated k-point
    REAL(KIND = DP) :: ws(3)
    !! Wigner-Seitz vector
    REAL(KIND = DP) :: rws(4, nrwsx)
    !! Real WS vectors
    REAL(KIND = DP) :: s_xk_border(3, 27)
    !! Look for special points on the border
    !
    ! Split the k-point across cores
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    rws(:, :) = zero
    CALL wsinit(rws, nrwsx, nrws, bg)
    xkt_sp(:, :) = 0
    !
    nb_sp = 0
    DO ik = 1, nkf
      counter = 0
      ! We could skip nb==1 to avoid identity symmetry
      DO nb = 1, nsym
        sa(:, :) = DBLE(s(:, :, nb))
        sb       = MATMUL(bg, sa)
        sr(:, :) = MATMUL(at, TRANSPOSE(sb))
        sr       = TRANSPOSE(sr)
        xk = xkf_all(:, ik + lower_bnd - 1)
        CALL cryst_to_cart(1, xk, bg, 1)
        CALL backtows(xk, ws, rws, nrwsx, nrws)
        xk = ws
        CALL DGEMV('n', 3, 3, 1.d0, sr, 3, xk, 1 ,0.d0 , S_xk, 1)
        !
        ! Needed for border_points
        counter_n = 0
        DO n = -1, 1
          DO m = -1, 1
            DO l = -1, 1
              counter_n = counter_n + 1
              s_xk_border(:, counter_n) = s_xk(:) + REAL(n, KIND = DP) * bg(:, 1) &
                         + REAL(m, KIND = DP) * bg (:, 2) + REAL(l, KIND = DP) * bg (:, 3)
            ENDDO
          ENDDO
        ENDDO
        sym_found = .FALSE.
        !
        IF (DOT_PRODUCT(xk - S_xk, xk - S_xk) < eps6) THEN
          counter = counter + 1
          xkt_sp(counter, ik + lower_bnd - 1) = nb
          sym_found = .TRUE.
        ENDIF
        ! Now check if the symmetry was not found because the point is on border
        IF (.NOT. sym_found) THEN
          DO counter_n = 1, 27
            IF (DOT_PRODUCT(xk(:) - s_xk_border(:, counter_n), xk(:) - s_xk_border(:, counter_n)) < eps6) THEN
              counter = counter + 1
              xkt_sp(counter, ik + lower_bnd - 1) = nb
            ENDIF
          ENDDO
        ENDIF
        !
      ENDDO ! nb
      IF (counter > 1) THEN
        nb_sp = nb_sp + 1
      ENDIF
      !
    ENDDO ! ik
    !
    ! Gather from all cores
    CALL mp_sum(xkt_sp, world_comm)
    CALL mp_sum(nb_sp, world_comm)
    !
    !! 48 symmetries + 1 index for the index of kpt
    ALLOCATE(xkf_sp(49, nb_sp), STAT = ierr)
    IF (ierr /= 0) CALL errore('special_points', 'Error allocating xkf_sp', 1)
    xkf_sp(:, :) = 0
    !
    counter = 0
    DO ik = 1, nktotf
      IF (xkt_sp(2, ik) > 0) THEN
        counter = counter + 1
        xkf_sp(1, counter) = ik
        xkf_sp(2:49, counter) = xkt_sp(:, ik)
      ENDIF
    ENDDO ! ik
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE special_points
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kpmq_map(xk, xq, sign1, nkq)
    !-----------------------------------------------------------------------
    !!
    !! this routine finds the index of k+q or k-q point on the fine k-mesh
    !!
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : nkf1, nkf2, nkf3
    USE constants_epw, ONLY : eps5
    USE mp,        ONLY : mp_bcast, mp_barrier
    USE kfold,     ONLY : backtoBZ
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: sign1
    !! +1 for searching k+q, -1 for k-q
    INTEGER, INTENT(out) :: nkq
    !! the index of k+sign*q
    !
    REAL(KIND = DP), INTENT(in) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! coordinates of q points
    !
    ! Local variables
    LOGICAL :: in_the_list
    !! Check if k point is in the list
    !
    REAL(KIND = DP) :: xx, yy, zz
    !! Temporary variables
    REAL(KIND = DP) :: xxk(3)
    !! k + (sign1) * q
    !
    xxk(:) = xk(:) + DBLE(sign1) * xq(:)
    xx = xxk(1) * nkf1
    yy = xxk(2) * nkf2
    zz = xxk(3) * nkf3
    in_the_list = ABS(xx - NINT(xx)) <= eps5 .AND. &
                  ABS(yy - NINT(yy)) <= eps5 .AND. &
                  ABS(zz - NINT(zz)) <= eps5
    IF (.NOT. in_the_list) CALL errore('kpmq_map', 'k+q does not fall on k-grid', 1)
    !
    !  find the index of this k+q or k-q in the k-grid
    !  make sure xx, yy, zz are in the 1st BZ
    !
    CALL backtoBZ(xx, yy, zz, nkf1, nkf2, nkf3)
    !
    ! since k- and q- meshes are commensurate, nkq can be easily found
    !
    nkq = NINT(xx) * nkf2 * nkf3 + NINT(yy) * nkf3 + NINT(zz) + 1
    !
    !  Now nkq represents the index of k+sign*q on the fine k-grid.
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kpmq_map
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE k_avg(f_out, vkk_all, nb_sp, xkf_sp)
    !-----------------------------------------------------------------------
    !!
    !! This routines enforces symmetry.
    !! Averages points which leaves the k-point unchanged by symmetry
    !!   e.g. k=[1,1,1] and q=[1,0,0] with the symmetry that change x and y gives
    !!        k=[1,1,1] and q=[0,1,0].
    !!
    !! Samuel Ponce & Francesco Macheda
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nstemp
    USE elph2,         ONLY : nkf, nbndfst, nktotf
    USE cell_base,     ONLY : bg, at
    USE constants_epw, ONLY : eps6, zero
    USE symm_base,     ONLY : s, nsym
    USE division,      ONLY : fkbounds
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : world_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nb_sp
    !! Lenght of xkf_sp
    INTEGER, INTENT(in) :: xkf_sp(49, nb_sp)
    !! Special points indexes and symmetries
    REAL(KIND = DP), INTENT(inout) :: f_out(3, nbndfst, nktotf, nstemp)
    !! In solution for iteration i
    REAL(KIND = DP), INTENT(inout) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity of k
    !
    ! Local variables
    LOGICAL :: special_map(nkf)
    !! Special mapping
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Band index
    INTEGER :: lower_bnd
    !! Lower bounds index after k para
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: sp
    !! Local index
    INTEGER ::nb
    !! Local index
    INTEGER :: counter_average
    !! Local counter
    INTEGER :: index_sp(nkf)
    !! Index of special points
    REAL(KIND = DP) :: sa(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sb(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: sr(3, 3)
    !! Symmetry matrix in crystal
    REAL(KIND = DP) :: s_vkk(3, nbndfst)
    !! Rotated vector
    REAL(KIND = DP) :: s_f_out(3, nbndfst)
    !! Rotated vector
    REAL(KIND = DP) :: tmp_vkk(3, nbndfst)
    !! Temporary vector
    REAL(KIND = DP) :: tmp_f_out(3, nbndfst)
    !! Temporary vector
    REAL(KIND = DP) :: f_out_loc(3, nbndfst, nktotf, nstemp)
    !! Local f_out where the k-points have been spread
    REAL(KIND = DP) :: vkk_all_loc(3, nbndfst, nktotf)
    !! Local velocity where the k-points have been spread
    !
    ! Split the k-point across cores
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    f_out_loc(:, :, :, :) = zero
    vkk_all_loc(:, :, :) = zero
    special_map(:) = .FALSE.
    index_sp(:) = 0
    DO ik = 1, nkf
      DO sp = 1, nb_sp
        IF (ik + lower_bnd - 1 == xkf_sp(1, sp)) THEN
          special_map(ik) = .TRUE.
          index_sp(ik) = sp
        ENDIF
      ENDDO
    ENDDO ! ik
    !
    DO itemp = 1, nstemp
      DO ik = 1, nkf
        IF (special_map(ik)) THEN
          counter_average = 0
          tmp_vkk = zero
          tmp_f_out = zero
          DO nb = 1, nsym
            IF (index_sp(ik) > 0) THEN
              IF (xkf_sp(nb + 1, index_sp(ik)) > 0) THEN
                counter_average = counter_average + 1
                sa(:, :) = DBLE(s(:, :, xkf_sp(nb + 1, index_sp(ik))))
                sb       = MATMUL(bg, sa)
                sr(:, :) = MATMUL(at, TRANSPOSE(sb))
                sr       = TRANSPOSE(sr)
                DO ibnd = 1, nbndfst
                  CALL DGEMV('n', 3, 3, 1.d0, sr, 3, vkk_all(:, ibnd, ik + lower_bnd - 1), 1, 0.d0, s_vkk(:, ibnd), 1)
                  CALL DGEMV('n', 3, 3, 1.d0, sr, 3, f_out(:, ibnd, ik + lower_bnd - 1, itemp), 1, 0.d0, s_f_out(:, ibnd), 1)
                  tmp_vkk(:, ibnd) = tmp_vkk(:, ibnd) + s_vkk(:, ibnd)
                  tmp_f_out(:, ibnd) = tmp_f_out(:, ibnd) + s_f_out(:, ibnd)
                ENDDO ! ibnd
              ENDIF
            ENDIF
          ENDDO ! sp
          DO ibnd = 1, nbndfst
            vkk_all_loc(:, ibnd, ik + lower_bnd - 1) = tmp_vkk(:, ibnd) / DBLE(counter_average)
            f_out_loc(:, ibnd, ik + lower_bnd - 1, itemp) = tmp_f_out(:, ibnd) / DBLE(counter_average)
          ENDDO
          !
        ELSE ! not a special point
          DO ibnd = 1, nbndfst
            vkk_all_loc(:, ibnd, ik + lower_bnd - 1) = vkk_all(:, ibnd, ik + lower_bnd - 1)
            f_out_loc(:, ibnd, ik + lower_bnd - 1, itemp) = f_out(:, ibnd, ik + lower_bnd - 1, itemp)
          ENDDO
        ENDIF ! special
      ENDDO! ik
    ENDDO! itemp
    !
    ! Gather from all cores
    CALL mp_sum(vkk_all_loc, world_comm)
    CALL mp_sum(f_out_loc, world_comm)
    !
    f_out = f_out_loc
    vkk_all = vkk_all_loc
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE k_avg
    !-----------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------
  END MODULE grid
  !-----------------------------------------------------------------------
