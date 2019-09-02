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
    SUBROUTINE loadkmesh_para
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
    USE io_epw,    ONLY : iunkf
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    !
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
    !! Error number
    !
    REAL(KIND = DP), ALLOCATABLE :: xkf_(:, :), xkf_tmp(:, :), xkfval(:, :)
    !! coordinates k-points
    REAL(KIND = DP), ALLOCATABLE :: wkf_(:), wkf_tmp(:)
    !! weights k-points
    !
    IF (mpime == ionode_id) THEN
      IF (filkf /= '') THEN ! load from file (crystal coordinates)
        !
        WRITE(stdout, *) '    Using k-mesh file: ', TRIM(filkf)
        OPEN(UNIT = iunkf, FILE = filkf, STATUS = 'old', FORM = 'formatted', ERR = ierr, IOSTAT = ios)
        IF (ierr /= 0)ALL errore('loadkmesh_para', 'opening file ' // filkf, ABS(ios))
        READ(iunkf, *) nkqtotf 
        !
        ALLOCATE(xkf_(3, 2 * nkqtotf))
        ALLOCATE(wkf_(2 * nkqtotf))
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
          !
          !  bring the k point to crystal coordinates
          ! CALL cryst_to_cart( 1, xkf_ (:,ikk), at, -1)
          !
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
      ELSEIF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_k) THEN
          ! get size of the mp_mesh in the irr wedge 
          WRITE(stdout,'(a,3i4)') '     Using uniform MP k-mesh: ', nkf1, nkf2, nkf3
          call set_sym_bl()
          !
          ALLOCATE(xkf_(3, 2 * nkf1 * nkf2 * nkf3))
          ALLOCATE(wkf_(2 * nkf1 * nkf2 * nkf3))
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nrot, time_reversal, .FALSE., s, t_rev, bg, nkf1 * nkf2 * nkf3, &
               0, 0, 0, nkf1, nkf2, nkf3, nkqtotf, xkf_, wkf_)
          DEALLOCATE(xkf_)
          DEALLOCATE(wkf_)
          ALLOCATE(xkf_(3, 2 * nkqtotf))
          ALLOCATE(wkf_(2 * nkqtotf)) 
          ALLOCATE(xkf_tmp(3, nkqtotf))
          ALLOCATE(wkf_tmp(nkqtotf))
          ALLOCATE(xkfval(3, 2 * nkqtotf))   
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
          DEALLOCATE(xkf_tmp)
          DEALLOCATE(wkf_tmp)
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
          DEALLOCATE(xkfval)
          !
        ELSE ! mp_mesh_k
          !
          WRITE (stdout,'(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          !
          nkqtotf = 2 * nkf1 * nkf2 * nkf3
          ALLOCATE(xkf_ (3, nkqtotf))
          ALLOCATE(wkf_(nkqtotf))
          wkf_(:) = 0.d0
          DO ik = 1, nkf1 * nkf2 * nkf3
            wkf_(2 * ik - 1) = 2.d0 / (DBLE(nktotf))
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
        ALLOCATE(xkf_(3, 2 * nkqtotf))
        ALLOCATE(wkf_(2 * nkqtotf))
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
    nkqf = 2 * (nktotf / npool)
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
    IF (.NOT. ALLOCATED(xkf_)) ALLOCATE(xkf_(3, nkqtotf))
    IF (.NOT. ALLOCATED(wkf_)) ALLOCATE(wkf_(nkqtotf))
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
    ALLOCATE(xkf(3, nkqf))
    ALLOCATE(wkf(nkqf))
    xkf(:, :) = xkf_(:, lower_bnd:upper_bnd)
    ! 
    ! KMB: set coordinates of displaced vectors for indabs
    IF (vme .AND. eig_read) THEN
      ALLOCATE(xkfd(3, nkqf, 6)) 
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
    IF (ABS(sum (wkf_ (:)) - 2.d0) > 1.d-4 ) &
      WRITE(stdout, '(5x,"WARNING: k-point weigths do not add up to 1 [loadkmesh_para]")')
    !
    WRITE(stdout, '(5x,"Size of k point mesh for interpolation: ",i10)') nkqtotf 
    WRITE(stdout, '(5x,"Max number of k points per pool:",7x,i10)') nkqf 
    !
    DEALLOCATE(xkf_)
    DEALLOCATE(wkf_)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE loadkmesh_para
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE loadkmesh_serial
    !-----------------------------------------------------------------------
    !!
    !!  Load fine k mesh in sequential
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
    USE io_epw,    ONLY : iunkf
    USE constants_epw, ONLY : eps4
    !
    IMPLICIT NONE
    !
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
    !! Error number. 
    REAL(KIND = DP), ALLOCATABLE :: xkf_tmp(:, :)
    !! coordinates k-points
    REAL(KIND = DP), ALLOCATABLE :: wkf_tmp(:)
    !! weights k-points
    !
    IF (mpime == ionode_id) THEN
      IF (filkf /= '') THEN ! load from file (crystal coordinates)
        !
        ! Each pool gets its own copy from the action=read statement
        !
        WRITE(stdout, *) '     Using k-mesh file: ', TRIM(filkf)
        OPEN(UNIT = iunkf, FILE = filkf, STATUS = 'old', FORM = 'formatted', ERR = ierr, IOSTAT = ios)
        IF (ierr /= 0) CALL errore('loadkmesh_serial', 'opening file '//filkf, ABS(ios))
        READ(iunkf, *) nkqtotf
        ALLOCATE(xkf(3, 2 * nkqtotf))
        ALLOCATE(wkf(2 * nkqtotf))
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
          !
          !  bring the k point to crystal coordinates
          ! CALL cryst_to_cart( 1, xkf_ (:,ikk), at, -1)
          !
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
        ! bring xkf in crystal coordinates
        ! CALL cryst_to_cart(nkqtotf, xkf, at, -1)
        !
      ELSEIF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_k) THEN
          ! get size of the mp_mesh in the irr wedge 
          WRITE (stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          CALL set_sym_bl()
          !                                         
          ALLOCATE(xkf(3, 2 * nkf1 * nkf2 * nkf3))
          ALLOCATE(wkf(2 * nkf1 * nkf2 * nkf3))
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nrot, time_reversal, s, t_rev, bg, nkf1 * nkf2 * nkf3, &
               0, 0, 0, nkf1, nkf2, nkf3, nkqtotf, xkf, wkf)
          DEALLOCATE(xkf)
          DEALLOCATE(wkf) 
          ALLOCATE(xkf(3, 2 * nkqtotf))
          ALLOCATE(wkf(2 * nkqtotf))
          ALLOCATE(xkf_tmp(3, nkqtotf))
          ALLOCATE(wkf_tmp(nkqtotf))
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
          DEALLOCATE(xkf_tmp)
          DEALLOCATE(wkf_tmp)
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
          ALLOCATE(xkf(3, nkqtotf))
          ALLOCATE(wkf(nkqtotf))
          wkf(:) = 0.d0
          DO ik = 1, nkf1 * nkf2 * nkf3
            wkf(2 * ik - 1) = 2.d0 / (DBLE(nktotf))
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
        ALLOCATE(xkf(3, 2 * nkqtotf))
        ALLOCATE(wkf(2 * nkqtotf))
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
      nkf = nktotf
      nkqf = nkqtotf
      !
    ENDIF
    !
    CALL mp_bcast(nkf, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkqf, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkqtotf, ionode_id, inter_pool_comm)
    IF (mpime /= ionode_id) THEN
      ALLOCATE(xkf(3, nkqtotf))
      ALLOCATE(wkf(nkqtotf))
    ENDIF
    CALL mp_bcast(xkf, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkf, ionode_id, inter_pool_comm)
    !
    ! KMB: set coordinates of displaced vectors - indabs
    IF (vme .AND. eig_read) THEN
      ALLOCATE(xkfd(3, nkqf, 6)) 
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
    USE epwcom,    ONLY : filqf, nkf1, nkf2, nkf3
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nktotbz
    !! Total number of k-points 
    REAL(KIND = DP), INTENT(inout) :: xkf_bz(3, nktotbz)
    !! Return the grid on full BZ
    !
    INTEGER :: ik, i, j, k, ios
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
    !-----------------------------------------------------------------------
    SUBROUTINE loadqmesh_para
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
    USE io_epw,    ONLY : iunqf
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    !
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
    INTEGER :: ierr
    !! Error number.
    INTEGER :: rest
    !! Remaining of cores numbers
    REAL(KIND = DP), ALLOCATABLE :: xqf_(:, :)
    !! Temporary q-point
    REAL(KIND = DP), ALLOCATABLE :: wqf_(:)
    !! Temporary weight of q-point
    !
    IF (mpime == ionode_id) THEN
      IF (filqf /= '') THEN ! load from file (crystal coordinates)
        !
        WRITE(stdout, *) '    Using q-mesh file: ', TRIM(filqf)
        IF (lscreen) WRITE(stdout, *) '     WARNING: if lscreen=.TRUE., q-mesh needs to be [-0.5:0.5] (crystal)' 
        OPEN(UNIT = iunqf, FILE = filqf, STATUS = 'old', FORM = 'formatted', ERR = ierr, IOSTAT = ios)
        IF (ierr /= 0) CALL errore('loadkmesh_para', 'Opening file ' // filqf, ABS(ios))
        READ(iunqf, *) nqtotf
        !
        ALLOCATE(xqf_(3, nqtotf))
        ALLOCATE(wqf_(nqtotf))
        !
        DO iq = 1, nqtotf
          !
          READ(iunqf, *) xqf_(:, iq), wqf_(iq)
          !
        ENDDO
        CLOSE(iunqf)
        !
      ELSEIF ((nqf1 /= 0) .AND. (nqf2 /= 0) .AND. (nqf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_q) THEN
          IF (lscreen) CALL errore('loadqmesh', 'If lscreen = .TRUE. do not use mp_mesh_q', 1)
          ! get size of the mp_mesh in the irr wedge 
          WRITE(stdout, '(a,3i4)') '     Using uniform MP q-mesh: ', nqf1, nqf2, nqf3
          call set_sym_bl()
          !
          ALLOCATE(xqf_ (3, nqf1 * nqf2 * nqf3))
          ALLOCATE(wqf_(nqf1 * nqf2 * nqf3))
          ! the result of this call is just nkqtotf
          CALL kpoint_grid ( nrot, time_reversal, .FALSE., s, t_rev, bg, nqf1*nqf2*nqf3, &
               0,0,0, nqf1,nqf2,nqf3, nqtotf, xqf_, wqf_)
          DEALLOCATE(xqf_, wqf_)
          ALLOCATE(xqf_ (3, nqtotf), wqf_(nqtotf)) 
          CALL kpoint_grid ( nrot, time_reversal, .FALSE., s, t_rev, bg, nqf1*nqf2*nqf3, &
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
          ALLOCATE(xqf_(3, nqtotf), wqf_(nqtotf) )
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
        ALLOCATE(xqf_ (3, nqtotf))
        ALLOCATE(wqf_(nqtotf))
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
      ALLOCATE(xqf_(3, nqtotf))
      ALLOCATE(wqf_(nqtotf))
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
    ALLOCATE(xqf(3, nqf))
    ALLOCATE(wqf(nqf))
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
    DEALLOCATE(xqf_)
    DEALLOCATE(wqf_)
    !-----------------------------------------------------------------------
    END SUBROUTINE loadqmesh_para
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE loadqmesh_serial
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
    USE io_epw,    ONLY : iunqf
    USE constants_epw, ONLY : eps4
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iq
    !! Q-index 
    INTEGER :: i, j, k
    !! Directions
    INTEGER :: ios
    !! Status integer
    INTEGER :: ierr
    !! Error number. 
    !
    IF (mpime == ionode_id) THEN
      IF (filqf /= '') THEN ! load from file (crystal coordinates)
        !
        ! Each pool gets its own copy from the action=read statement
        !
        WRITE(stdout, *) '    Using q-mesh file: ', TRIM(filqf)
        IF (lscreen) WRITE(stdout, *) '     WARNING: if lscreen=.TRUE., q-mesh needs to be [-0.5:0.5] (crystal)'
        OPEN(UNIT = iunqf, FILE = filqf, STATUS = 'old', FORM = 'formatted', ERR = ierr, IOSTAT = ios)
        IF (ierr /= 0) CALL errore('loadqmesh_serial', 'opening file ' // filqf, ABS(ios))
        READ(iunqf, *) nqtotf
        ALLOCATE(xqf(3, nqtotf))
        ALLOCATE(wqf(nqtotf))
        DO iq = 1, nqtotf
          READ (iunqf, *) xqf(:, iq), wqf(iq)
        ENDDO
        CLOSE(iunqf)
        !
        ! bring xqf in crystal coordinates
        ! CALL cryst_to_cart(nqtotf, xqf, at, -1)
        !
      ELSEIF ((nqf1 /= 0) .AND. (nqf2 /= 0) .AND. (nqf3 /= 0)) THEN ! generate grid
        IF (mp_mesh_q) THEN
          IF (lscreen) CALL errore ('loadqmesh', 'If lscreen=.TRUE. do not use mp_mesh_q',1)
          ! get size of the mp_mesh in the irr wedge 
          WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
          call set_sym_bl()
          !                                         
          ALLOCATE(xqf(3, nqf1 * nqf2 * nqf3))
          ALOCATE(wqf(nqf1 * nqf2 * nqf3))
          ! the result of this call is just nkqtotf
          CALL kpoint_grid(nrot, time_reversal, s, t_rev, bg, nqf1 * nqf2 * nqf3, &
               0, 0, 0, nqf1, nqf2, nqf3, nqtotf, xqf, wqf)
          DEALLOCATE(xqf)
          DEALLOCATE(wqf) 
          ALLOCATE(xqf(3, nqtotf))
          ALLOCATE(wqf(nqtotf)) 
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
          ALLOCATE(xqf (3, nqtotf))
          ALLOCATE(wqf(nqtotf))
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
        ALLOCATE(xqf(3, nqtotf))
        ALLOCATE(wqf(nqtotf))
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
      ALLOCATE(xqf(3, nqtotf))
      ALLOCATE(wqf(nqtotf))
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
    USE elph2,         ONLY : nqf, xqf, xkf, chw, etf, nkf, nqtotf, nkqtotf, &
                              map_rebal
    USE io_global,     ONLY : ionode_id, stdout
    USE io_epw,        ONLY : iunselecq
    USE mp_global,     ONLY : npool, inter_pool_comm, world_comm, my_pool_id
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
    INTEGER :: ik, ikk, ikq, ikl
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
    INTEGER :: nkqtotf_tmp
    !! Temporary k-q points.
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: BZtoIBZ_tmp(nkf1 * nkf2 * nkf3)
    !! Temporary mapping
    INTEGER :: BZtoIBZ(nkf1 * nkf2 * nkf3)
    !! BZ to IBZ mapping
    INTEGER :: s_BZtoIBZ(3, 3, nkf1 * nkf2 * nkf3)
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
    REAL(KIND = DP) :: xkf_tmp (3, nkqtotf)
    !! Temporary k-point coordinate (dummy variable)
    REAL(KIND = DP) :: wkf_tmp(nkqtotf)
    !! Temporary k-weights (dummy variable)
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
        READ (iunselecq,*) totq
        ALLOCATE(selecq(totq))
        selecq(:) = 0
        READ (iunselecq,*) nqtot
        READ (iunselecq,*) selecq(:)
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
      ALLOCATE(selecq(nqf))
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
          CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
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
        CALL poolgather(nbndsub, nktotf, nkf, etf_loc, etf_all )
        ! 
        ! In case of k-point symmetry
        IF (mp_mesh_k) THEN
          BZtoIBZ(:) = 0
          s_BZtoIBZ(:, :, :) = 0
          ! 
          IF (mpime == ionode_id) THEN
            ! 
            CALL set_sym_bl()
            !
            ! What we get from this call is BZtoIBZ
            CALL kpoint_grid_epw(nrot, time_reversal, .FALSE., s, t_rev, bg, nkf1 * nkf2 * nkf3, &
                       nkf1, nkf2, nkf3, nkqtotf_tmp, xkf_tmp, wkf_tmp, BZtoIBZ, s_BZtoIBZ)
            ! 
            IF (iterative_bte) THEN
              BZtoIBZ_tmp(:) = 0
              DO ikbz = 1, nkf1 * nkf2 * nkf3
                BZtoIBZ_tmp(ikbz) = map_rebal(BZtoIBZ(ikbz))
              ENDDO
              BZtoIBZ(:) = BZtoIBZ_tmp(:)
            ENDIF
            ! 
          ENDIF ! mpime
          CALL mp_bcast(BZtoIBZ, ionode_id, inter_pool_comm)
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
              WRITE(stdout,'(5x,a,i12,i12)')'Number selected, total', totq, iq
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
          CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
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
              CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1)
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
              CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1)
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
                WRITE(stdout,i '(5x,a,i12,i12)')'Number selected, total', totq, iq
              ENDIF
            ENDIF
          ENDDO ! iq            
        ENDIF ! use_ws
      ENDIF ! homogeneous
      !  
      IF (mpime == ionode_id) THEN
        OPEN(UNIT = iunselecq, FILE = 'selecq.fmt', ACTION = 'write')
        WRITE (iunselecq,*) totq    ! Selected number of q-points
        WRITE (iunselecq,*) nqtotf  ! Total number of q-points 
        WRITE (iunselecq,*) selecq(1:totq)
        CLOSE(iunselecq)
      ENDIF
      ! 
    ENDIF ! exst
    !----------------------------------------------------------------------- 
    END SUBROUTINE qwindow
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE load_rebal
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
    USE io_global,     ONLY : stdout
    USE elph2,         ONLY : etf, nkf, nkqtotf, xkf, wkf, etf, map_rebal, map_rebal_inv
    USE epwcom,        ONLY : fsthick, nbndsub, mp_mesh_k
    USE pwcom,         ONLY : ef
    USE mp_global,     ONLY : my_pool_id, npool
    USE constants_epw, ONLY : zero
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm, mp_bcast
    USE transportcom,  ONLY : lower_bnd
    !
    IMPLICIT NONE
    !  
    INTEGER :: pool_index(npool) 
    !! Index of the current pool
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
    INTEGER :: kpt_in(nkqtotf)
    !! K-points that are within the fshick windows
    INTEGER :: kpt_out(nkqtotf)
    !! K-points that are outside of the fshick windows
    INTEGER :: map_rebal_tmp(nktotf)
    !! Temporary map between the initial ordering of k-point and the rebalanced one
    INTEGER :: map_rebal_inv_tmp(nktotf)
    !! Temporary inverse map between the initial ordering of k-point and the rebalanced one
    !
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
    ALLOCATE(map_rebal(nktotf))
    ALLOCATE(map_rebal_inv(nktotf))
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
    rest = ((nktotf) - tot * npool)
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
  END MODULE grid
  !-----------------------------------------------------------------------
