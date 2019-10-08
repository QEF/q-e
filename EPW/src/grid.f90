  !
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
                          rand_k, rand_nk, mp_mesh_k, system_2d, eig_read, vme
    USE elph2,     ONLY : nkqtotf, nkqf, xkf, wkf, nkf, xkfd, deltaq
    USE cell_base, ONLY : at, bg
    USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
    USE io_var,    ONLY : iunkf
    USE low_lvl,   ONLY : init_random_seed
    USE constants_epw, ONLY : eps4
    USE noncollin_module, ONLY : noncolin
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
    REAL(KIND = DP), ALLOCATABLE :: xkf_(:, :)
    !! coordinates k-points
    REAL(KIND = DP), ALLOCATABLE :: xkf_tmp(:, :)
    !! Temporary k-point
    REAL(KIND = DP), ALLOCATABLE :: xkfval(:, :)
    !! K-point
    REAL(KIND = DP), ALLOCATABLE :: wkf_(:)
    !! weights k-points
    REAL(KIND = DP), ALLOCATABLE :: wkf_tmp(:)
    !! Temporary weights
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
        IF (imatches("cartesian", coordinate_type)) THEN
          CALL cryst_to_cart(nkqtotf, xkf_, at, -1)
        ENDIF
        !
        ! redefine nkqtotf to include the k+q points
        !
        nkqtotf = 2 * nkqtotf
        !
      ELSEIF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_k) THEN
          ! get size of the mp_mesh in the irr wedge 
          WRITE(stdout,'(a,3i4)') '     Using uniform MP k-mesh: ', nkf1, nkf2, nkf3
          call set_sym_bl()
          !
          ALLOCATE(xkf_(3, 2 * nkf1 * nkf2 * nkf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_', 1)
          ALLOCATE(wkf_(2 * nkf1 * nkf2 * nkf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nrot, time_reversal, .FALSE., s, t_rev, bg, nkf1 * nkf2 * nkf3, &
               0, 0, 0, nkf1, nkf2, nkf3, nkqtotf, xkf_, wkf_)
          DEALLOCATE(xkf_, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error deallocating xkf_', 1)
          DEALLOCATE(wkf_, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error deallocating wkf_', 1)
          ALLOCATE(xkf_(3, 2 * nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_', 1)
          ALLOCATE(wkf_(2 * nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
          ALLOCATE(xkf_tmp(3, nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_tmp', 1)
          ALLOCATE(wkf_tmp(nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_tmp', 1)
          ALLOCATE(xkfval(3, 2 * nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkfval', 1)
          xkf_(:, :) = 0.0d0
          xkfval(:, :) = 0.0d0
          CALL kpoint_grid(nrot, time_reversal, .FALSE., s, t_rev, bg, nkf1 * nkf2 * nkf3, &
               0, 0, 0, nkf1, nkf2, nkf3, nkqtotf, xkf_tmp, wkf_tmp)
          !  
          ! assign to k and k+q for xkf and wkf 
          ! 
          ! SP: The variable xkfval is a duplication. However, it allows to avoid some strange 
          !     memory allocation issue. FIXME
          DO ik = 1, nkqtotf
            ikk = 2 * ik - 1
            ikq = ikk + 1
            xkf_(:, ikk)   = xkf_tmp(:, ik)
            xkf_(:, ikq)   = xkf_tmp(:, ik)
            xkfval(:, ikk) = xkf_tmp(:, ik)
            xkfval(:, ikq) = xkf_tmp(:, ik)
            wkf_(ikk)   = 2.d0 * wkf_tmp(ik)
            wkf_(ikq)   = 0.d0
          ENDDO
          DEALLOCATE(xkf_tmp, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error deallocating xkf_tmp', 1)
          DEALLOCATE(wkf_tmp, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error deallocating wkf_tmp', 1)
          !       
          ! bring the k point to crystal coordinates       
          CALL cryst_to_cart(2 * nkqtotf, xkfval, at, -1)
          xkf_(:, :) = xkfval(:, :)
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
          DEALLOCATE(xkfval, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error deallocating xkfval', 1)
          !
        ELSE ! mp_mesh_k
          !
          WRITE(stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          !
          nkqtotf = 2 * nkf1 * nkf2 * nkf3
          ALLOCATE(xkf_ (3, nkqtotf), STAT = ierr)
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
          IF (system_2d) THEN
            CALL random_number(xkf_(1:2, ikk))
            xkf_(3, ikk) = 0.d0
          ELSE
            CALL random_number(xkf_(:, ikk))
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
    ENDIF
    CALL mp_bcast(xkf_, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkf_, ionode_id, inter_pool_comm)
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
    IF (vme .AND. eig_read) THEN
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
    USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
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
        IF (imatches("cartesian", coordinate_type)) THEN
          CALL cryst_to_cart(nkqtotf, xkf, at, -1)
        ENDIF
        !
        ! redefine nkqtotf to include the k+q points
        !
        nkqtotf = 2 * nkqtotf
        !
        ! bring xkf in crystal coordinates
        ! CALL cryst_to_cart(nkqtotf, xkf, at, -1)
        !
      ELSEIF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_k) THEN
          ! get size of the mp_mesh in the irr wedge 
          WRITE(stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          CALL set_sym_bl()
          !                                         
          ALLOCATE(xkf(3, 2 * nkf1 * nkf2 * nkf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating xkf', 1)
          ALLOCATE(wkf(2 * nkf1 * nkf2 * nkf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadkmesh_serial', 'Error allocating wkf', 1)
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nrot, time_reversal, s, t_rev, bg, nkf1 * nkf2 * nkf3, &
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
          CALL kpoint_grid(nrot, time_reversal, s, t_rev, bg, nkf1 * nkf2 * nkf3, &
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
          WRITE (stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
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
          IF (system_2d) THEN
            CALL random_number(xkf(1:2, ikk))
            xkf(3, ikk) = 0.d0
          ELSE
            CALL random_number(xkf(:, ikk))
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
    IF (vme .AND. eig_read) THEN
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
    SUBROUTINE loadkmesh_fullBZ(xkf_bz, nktotbz)
    !-----------------------------------------------------------------------
    !!
    !!  Create a k-mesh for fine grid without symmetries on the full grid
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP, sgl
    USE epwcom,    ONLY : nkf1, nkf2, nkf3
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nktotbz
    !! Total number of k-points 
    REAL(KIND = DP), INTENT(inout) :: xkf_bz(3, nktotbz)
    !! Return the grid on full BZ
    !
    INTEGER :: ik
    !! K-point index
    INTEGER :: i, j, k
    !! K-point grid dim
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
    SUBROUTINE kpoint_grid_epw(nrot, time_reversal, skip_equivalence, s, t_rev, &
                               nkc1, nkc2, nkc3, BZtoIBZ, s_BZtoIBZ)
    !-----------------------------------------------------------------------
    !!
    !!  Automatic generation of a uniform grid of k-points with symmetry. 
    !!  Routine copied from PW/src/kpoint_grid.f90.
    !!  We had to duplicate because the BZtoIBZ array was deallocated and is needed in
    !!  EPW 
    !!
    USE kinds,            ONLY : DP
    USE division,         ONLY : fkbounds
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : world_comm, my_pool_id, npool, inter_pool_comm
    USE kinds_epw,        ONLY : SIK2
    USE constants_epw,    ONLY : eps6
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_INTEGER
#endif
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrot
    !! Number of Bravais symmetry
    INTEGER, INTENT(in) :: nkc1, nkc2, nkc3
    !! k-grid 
    INTEGER, INTENT(in) :: t_rev(48)
    !! Time-reversal sym
    INTEGER, INTENT(in) :: s(3, 3, 48)
    !! Symmetry matrice. 
    INTEGER(SIK2), INTENT(inout) :: s_BZtoIBZ(nkc1 * nkc2 * nkc3)
    !! Symeetry matrix that links an point to its IBZ friend.
    INTEGER, INTENT(inout) :: BZtoIBZ(nkc1 * nkc2 * nkc3)
    !! Number of rotation
    LOGICAL, INTENT(in) :: time_reversal
    !! True if time reversal
    LOGICAL, INTENT(in) :: skip_equivalence
    !! True if equivalent point
    !
    ! Local variables
    LOGICAL :: in_the_list
    !! Is the current point in the list
    INTEGER(SIK2) :: s_save(nkc1 * nkc2 * nkc3)
    !! Temporary symmetry matrix
    INTEGER :: nkr
    !! Total number of points
    INTEGER :: i, j, k
    !! Index on grid size
    INTEGER(SIK2) :: ns
    !! Index on symmetry operations
    INTEGER :: n
    !! Global k-point index
    INTEGER :: nk
    !! Equivalent point
    INTEGER :: equiv(nkc1 * nkc2 * nkc3)
    !! Equivalent k-points
    INTEGER :: ik
    !! K-point index 
    INTEGER :: lower_bnd
    !! K-point paralelization (lower-bound index)
    INTEGER :: upper_bnd
    !! K-point paralelization (upper-bound index) 
    INTEGER :: cumul_nks
    !! Sum of points
    INTEGER :: BZtoIBZ_tmp(nkc1 * nkc2 * nkc3)
    !! Temporrary BZtoIBZ map
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: nkspar(:)
    !! Number of irr points (IBZ)
    REAL(KIND = DP) :: xkr(3)
    !! Current point
    REAL(KIND = DP) :: xx, yy, zz
    !! Current point coordinate
    REAL(KIND = DP), ALLOCATABLE :: xkg(:, :)
    !! Current point
    REAL(KIND = DP), ALLOCATABLE :: wkk(:)
    !! Weight of the k-point
    !
    nkr = nkc1 * nkc2 * nkc3
    ALLOCATE(nkspar(npool), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error allocating nkspar', 1)
    ALLOCATE(xkg(3, nkr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error allocating xkg', 1)
    ALLOCATE(wkk(nkr), STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error allocating wkk', 1)
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
    IF (skip_equivalence) THEN
      CALL infomsg('kpoint_grid_epw', 'ATTENTION: skip check of k-points equivalence')
      wkk = 1.d0
    ELSE
      DO nk = 1, nkr
        !  check if this k-point has already been found equivalent to another
        IF (equiv(nk) == nk) THEN
          wkk(nk) = 1.0d0
          !  check if there are equivalent k-point to this in the list
          !  (excepted those previously found to be equivalent to another)
          !  check both k and -k
          DO ns = 1, nrot
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
    !        IF (time_reversal) THEN
    !           xx =-xkr(1)*nkc1 
    !           yy =-xkr(2)*nkc2 
    !           zz =-xkr(3)*nkc3 
    !           in_the_list=ABS(xx-NINT(xx))<=eps.AND.ABS(yy-NINT(yy))<=eps &
    !                                              .AND. ABS(zz-NINT(zz))<=eps
    !           IF (in_the_list) THEN
    !              i = mod ( nint (-xkr(1)*nkc1  + 2*nkc1), nkc1 ) + 1
    !              j = mod ( nint (-xkr(2)*nkc2  + 2*nkc2), nkc2 ) + 1
    !              k = mod ( nint (-xkr(3)*nkc3  + 2*nkc3), nkc3 ) + 1
    !              n = (k-1) + (j-1)*nkc3 + (i-1)*nkc2*nkc3 + 1
    !              IF (n>nk .AND. equiv(n)==n) THEN
    !                 equiv(n) = nk
    !                 wkk(nk)=wkk(nk)+1.0d0
    !                 s_save(:,:,n) = -s(:,:,ns)
    !              ELSE
    !                 IF (equiv(n)/=nk.OR.n<nk) CALL errore('kpoint_grid', &
    !                 'something wrong in the checking algorithm',2)
    !              ENDIF
    !           ENDIF
    !        ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    ! 
    !  count irreducible points and order them
    nkspar(:) = 0
    DO nk = 1, nkr
      BZtoIBZ(nk) = equiv(nk)
    ENDDO
    !
    CALL fkbounds(nkr, lower_bnd, upper_bnd)
    DO nk = lower_bnd, upper_bnd
      IF (equiv(nk) == nk) THEN
        nkspar(my_pool_id + 1) = nkspar(my_pool_id + 1) + 1
        IF (nkspar(my_pool_id + 1) > nkr) CALL errore('kpoint_grid_epw', 'Too many k-points', 1)
        BZtoIBZ(nk) = nkspar(my_pool_id + 1)
        ! Change all the one above
        DO ik = nk, nkr
          IF (equiv(ik) == nk) THEN
            BZtoIBZ(ik) = nkspar(my_pool_id + 1)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    !Now recompose the vector with the right order 
    CALL mp_sum(nkspar, world_comm)
    cumul_nks = 0
    IF (my_pool_id > 0) THEN
      DO i = 1, my_pool_id
        cumul_nks = cumul_nks + nkspar(i)
      ENDDO
    ENDIF
    DO ik = 1, nkr
      IF((BZtoIBZ(ik) > nkspar(my_pool_id + 1)) .OR. (ik < lower_bnd)) THEN
        BZtoIBZ (ik) = 0
      ELSE
        BZtoIBZ(ik) = BZtoIBZ(ik) + cumul_nks
      ENDIF
    ENDDO
    BZtoIBZ_tmp(:) = 0
    DO i = 1, npool
      IF (my_pool_id + 1 == i) THEN
        DO ik = 1, nkr
          IF (BZtoIBZ_tmp(ik) == 0) THEN
            BZtoIBZ_tmp(ik) = BZtoIBZ(ik)
          ENDIF
        ENDDO
      ENDIF
      CALL mp_bcast(BZtoIBZ_tmp, i - 1, inter_pool_comm)
    ENDDO
    !
    BZtoIBZ = BZtoIBZ_tmp
    !
    ! Now do the symmetry mapping. 
    DO nk = 1, nkr
      ! If its an irreducible point 
      IF (equiv(nk) == nk) THEN
        ! Then you have the identity matrix
        s_BZtoIBZ(nk) = 1
      ELSE
        s_BZtoIBZ(nk) = s_save(nk)  
      ENDIF
    ENDDO
    ! 
    DEALLOCATE(xkg, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error deallocating xkg', 1)
    DEALLOCATE(wkk, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error deallocating wkk', 1)
    DEALLOCATE(nkspar, STAT = ierr)
    IF (ierr /= 0) CALL errore('kpoint_grid_epw', 'Error deallocating nkspar', 1)
    ! 
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE kpoint_grid_epw
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
    USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
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
          call set_sym_bl()
          !
          ALLOCATE(xqf_(3, nqf1 * nqf2 * nqf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf_ ', 1)
          ALLOCATE(wqf_(nqf1 * nqf2 * nqf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating wqf_', 1)
          ! the result of this call is just nkqtotf
          CALL kpoint_grid( nrot, time_reversal, .FALSE., s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf_, wqf_)
          DEALLOCATE(xqf_, wqf_, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error deallocating xqf_, wqf_', 1)
          ALLOCATE(xqf_ (3, nqtotf), wqf_(nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_para', 'Error allocating xqf_ (3, nqtotf), wqf_', 1)
          CALL kpoint_grid( nrot, time_reversal, .FALSE., s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf_, wqf_)
          !  
          ! bring the k point to crystal coordinates       
          CALL cryst_to_cart(nqtotf, xqf_, at, -1)
          !
        ELSE
          !
          WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
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
          IF (system_2d) THEN
            CALL random_number(xqf_(1:2, iq))
            IF (lscreen .OR. specfun_pl .OR. plselfen) xqf_(1:2, iq) = xqf_(1:2, iq) - 0.5d0
            xqf_(3, iq) = 0.d0
          ELSE
            CALL random_number(xqf_(:, iq))
            IF (lscreen .OR. specfun_pl .OR. plselfen) xqf_(:, iq) = xqf_(:, iq) - 0.5d0
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
    USE io_global, ONLY : stdout
    USE epwcom,    ONLY : filqf, nqf1, nqf2, nqf3, &
                          rand_q, rand_nq, mp_mesh_q, system_2d, lscreen, &
                          plselfen, specfun_pl
    USE elph2,     ONLY : xqf, wqf, nqtotf, nqf
    USE cell_base, ONLY : at, bg
    USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
    USE io_var,    ONLY : iunqf
    USE low_lvl,   ONLY : init_random_seed
    USE constants_epw, ONLY : eps4
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
      ELSEIF ((nqf1 /= 0) .AND. (nqf2 /= 0) .AND. (nqf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_q) THEN
          IF (lscreen) CALL errore ('loadqmesh', 'If lscreen=.TRUE. do not use mp_mesh_q',1)
          ! get size of the mp_mesh in the irr wedge 
          WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
          call set_sym_bl()
          !                                         
          ALLOCATE(xqf(3, nqf1 * nqf2 * nqf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf', 1)
          ALLOCATE(wqf(nqf1 * nqf2 * nqf3), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nrot, time_reversal, s, t_rev, bg, nqf1 * nqf2 * nqf3, &
               0, 0, 0, nqf1, nqf2, nqf3, nqtotf, xqf, wqf)
          DEALLOCATE(xqf, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error deallocating xqf', 1)
          DEALLOCATE(wqf, STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error deallocating wqf', 1)
          ALLOCATE(xqf(3, nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating xqf', 1)
          ALLOCATE(wqf(nqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('loadqmesh_serial', 'Error allocating wqf', 1)
          CALL kpoint_grid(nrot, time_reversal, s, t_rev, bg, nqf1 * nqf2 * nqf3, &
               0,0,0, nqf1, nqf2, nqf3, nqtotf, xqf, wqf)
          !
          ! bring xqf in crystal coordinates       
          CALL cryst_to_cart(nqtotf, xqf, at, -1)
          !
        ELSE
          ! currently no offset.  
          ! q's are in crystal coordinates in xqf
          WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
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
          IF (system_2d) THEN
            CALL random_number(xqf(1:2, iq))
            IF (lscreen .OR. specfun_pl .OR. plselfen) xqf(1:2, iq) = xqf(1:2, iq) - 0.5d0
            xqf(3, iq) = 0.d0
          ELSE
            CALL random_number(xqf(:, iq))
            IF (lscreen .OR. specfun_pl .OR. plselfen) xqf(:, iq) = xqf(:, iq) - 0.5d0
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
    SUBROUTINE qwindow(exst, nrr_k, dims, totq, selecq, irvec_r, ndegen_k, & 
                       cufkk, cufkq, homogeneous)
    !-----------------------------------------------------------------------
    !!
    !! This routine pre-computes the q-points that falls within the fstichk.
    !! If at least 1 k-point is such that at least one k+q eigenenergy falls 
    !! within the user-defined fstichk, then the q-point is taken.  
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : nqf, xqf, xkf, chw, nkf, nqtotf, nkqtotf, &
                              map_rebal, nktotf
    USE io_global,     ONLY : ionode_id, stdout
    USE io_var,        ONLY : iunselecq
    USE mp_global,     ONLY : npool, world_comm, my_pool_id
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_sum, mp_bcast
    USE constants_epw, ONLY : twopi, ci, zero, eps6, ryd2ev, czero
    USE epwcom,        ONLY : nbndsub, fsthick, use_ws, mp_mesh_k, nkf1, nkf2, &
                              nkf3, iterative_bte, restart_freq, scissor
    USE noncollin_module, ONLY : noncolin
    USE pwcom,         ONLY : ef, nelec 
    USE cell_base,     ONLY : bg
    USE symm_base,     ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
    USE wan2bloch,     ONLY : hamwan2bloch
    USE io_eliashberg, ONLY : kpmq_map
    USE kinds_epw,     ONLY : SIK2
    USE poolgathering, ONLY : poolgather
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
    INTEGER, INTENT(inout) :: totq
    !! Total number of q-points inside fsthick
    INTEGER, ALLOCATABLE, INTENT(out) :: selecq(:)
    !! List of selected q-points
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
    INTEGER :: BZtoIBZ_tmp(nkf1 * nkf2 * nkf3)
    !! Temporary mapping
    INTEGER :: BZtoIBZ(nkf1 * nkf2 * nkf3)
    !! BZ to IBZ mapping
    INTEGER(SIK2) :: s_BZtoIBZ(nkf1 * nkf2 * nkf3)
    !! symmetry 
    INTEGER :: nkloc
    !! number of k-point selected on that cpu 
    INTEGER :: kmap(nkf)
    !! k-point that are selected for that cpu
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
        READ(iunselecq,*) totq
        ALLOCATE(selecq(totq), STAT = ierr)
        IF (ierr /= 0) CALL errore('qwindow', 'Error allocating selecq', 1)
        selecq(:) = 0
        READ(iunselecq,*) nqtot
        READ(iunselecq,*) selecq(:)
        CLOSE(iunselecq)
      ENDIF
      CALL mp_bcast(totq, ionode_id, world_comm)
      IF (mpime /= ionode_id) ALLOCATE(selecq(totq))
      CALL mp_bcast(nqtot , ionode_id, world_comm)
      CALL mp_bcast(selecq, ionode_id, world_comm)
      IF (nqtot /= nqtotf) THEN
        CALL errore('qwindow', 'Cannot read from selecq.fmt, the q-point grid or &
          & fsthick window are different from read one. Remove the selecq.fmt file and restart.', 1 )
      ENDIF
      !  
    ELSE
      ALLOCATE(selecq(nqf), STAT = ierr)
      IF (ierr /= 0) CALL errore('qwindow', 'Error allocating selecq', 1)
      selecq(:) = 0 
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
                  IF (ndegen_k(ir,iw2,iw) > 0) THEN
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
        IF (mp_mesh_k) THEN
          BZtoIBZ(:) = 0
          s_BZtoIBZ(:) = 0
          ! 
          CALL set_sym_bl()
          !
          ! What we get from this call is BZtoIBZ
          CALL kpoint_grid_epw(nrot, time_reversal, .FALSE., s, t_rev, nkf1, nkf2, nkf3, BZtoIBZ, s_BZtoIBZ)
          ! 
          IF (iterative_bte) THEN
            BZtoIBZ_tmp(:) = 0
            DO ikbz = 1, nkf1 * nkf2 * nkf3
              BZtoIBZ_tmp(ikbz) = map_rebal(BZtoIBZ(ikbz))
            ENDDO
            BZtoIBZ(:) = BZtoIBZ_tmp(:)
          ENDIF
          ! 
          ! 
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
        DO iq = 1, nqf
          xxq = xqf(:, iq)
          ! 
          found(:) = 0
          DO ik = 1, nkf 
            ikk = 2 * ik - 1
            xkk = xkf(:, ikk)
            xkq = xkk + xxq
            !  
            CALL kpmq_map(xkk, (/0d0,0d0,0d0/), 1, ind1) 
            CALL kpmq_map(xkk, xxq, 1, ind2) 
            ! 
            ! Use k-point symmetry
            IF (mp_mesh_k) THEN
              IF (((MINVAL(ABS(etf_all(:, BZtoIBZ(ind1)) - ef)) < fsthick) .AND. &
                    (MINVAL(ABS(etf_all(:, BZtoIBZ(ind2)) - ef)) < fsthick))) THEN
                found(my_pool_id + 1) = 1
                EXIT ! exit the loop 
              ENDIF
            ELSE
              IF (((MINVAL(ABS(etf_all(:, ind1) - ef)) < fsthick) .AND. &
                    (MINVAL(ABS(etf_all(:, ind2) - ef)) < fsthick))) THEN
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
            selecq(totq) = iq
            ! 
            IF (MOD(totq, restart_freq) == 0) THEN
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
            xxq = xqf(:, iq)
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
              selecq(totq) = iq
              IF (MOD(totq, restart_freq) == 0) THEN
                WRITE(stdout, '(5x,a,i12,i12)') 'Number selected, total', totq, iq
              ENDIF
            ENDIF
          ENDDO ! iq
        ELSE ! use_ws
          DO iq = 1, nqf
            xxq = xqf(:, iq)
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
              selecq(totq) = iq
              IF (MOD(totq, restart_freq) == 0) THEN
                WRITE(stdout, '(5x,a,i12,i12)')'Number selected, total', totq, iq
              ENDIF
            ENDIF
          ENDDO ! iq            
        ENDIF ! use_ws
      ENDIF ! homogeneous
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
    !! Total number of k-point
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
    tot = (nkqtotf / (2 * npool))
    rest = (nktotf - tot * npool)
    ! 
    DO ipool = 1, npool
      DO ik = 1,  tot
        map_rebal_inv_tmp(ik + (ipool - 1) * tot) = map_rebal_inv(npool * ik - (npool - ipool))
      ENDDO
    ENDDO
    ! Do the rest
    DO ik = 1, rest
      map_rebal_inv_tmp(ik + npool * tot) = map_rebal_inv(npool * tot + ik)
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
    USE symm_base,     ONLY : s, nrot
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
      DO nb = 1, nrot
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
    SUBROUTINE k_avg(F_out, vkk_all, nb_sp, xkf_sp)
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
    USE symm_base,     ONLY : s, nrot
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
    !! Local F_out where the k-points have been spread
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
      DO sp = 1,nb_sp
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
          tmp_F_out = zero
          DO nb = 1,nrot
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
    CALL mp_sum(F_out_loc, world_comm)
    ! 
    F_out = F_out_loc
    vkk_all = vkk_all_loc
    !  
    !-----------------------------------------------------------------------
    END SUBROUTINE k_avg
    !-----------------------------------------------------------------------
    ! 
  !-----------------------------------------------------------------------
  END MODULE grid
  !-----------------------------------------------------------------------
