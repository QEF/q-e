  !----------------------------------------------------------------------------
  MODULE phph
  !----------------------------------------------------------------------------
  !!
  !! This module contains routines related to phonon-phonon coupling
  !! Authored by Yiming Pan and Fabio Caruso
  !! Todo : In a serial running, the calling of MPI subroutines need
  !!        to be repalced.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !--------------------------------------------------------------------------   
    SUBROUTINE psi2_calc()
    !--------------------------------------------------------------------------
    !! This subroutine reads the thrid order force constants produced with 
    !! thirdorder.py and calculate the ph-ph matrix elements on the dense
    !! q grid, according to Eq.(9) of Comp. Phys. Commun. 185, 17471758 (2014)
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout, ionode,ionode_id
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast, mp_gather, &
                                 mp_max, mp_min
    USE ions_base,        ONLY : nat, amass, ityp, tau 
    USE cell_base,        ONLY : at, bg, alat, omega
    USE input,            ONLY : nqf1, nqf2, nqf3, degaussq, phwmin_tdbe,      &
                                 dg_tdbe, dt_tdbe, dph_tdbe
    USE ep_constants,     ONLY : kelvin2Ry, twopi, hbar, ryd2ev, rydcm1,  &
                                 ryd2mev 
    USE constants,        ONLY : pi
    USE global_var,       ONLY : wf,xqf,nqtotf,adapt_smearing
    USE tdbe_common,      ONLY : nifc, ind_atms, ifc3,rvecs, uf_all,      &
                                 vnuq_all, ps2sec,dt_in_ps, nind_p,       &
                                 nind_m, ind_p, ind_m, psi2_p, psi2_m
    USE modes,            ONLY : nmodes
    USE parallelism,      ONLY : fkbounds, fkbounds2
    USE io_var,           ONLY : iunphindx
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE,MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION,          &
                                 MPI_STATUS_IGNORE, MPI_INTEGER, MPI_IN_PLACE,   &
                                 MPI_MODE_RDONLY, MPI_INTEGER8, MPI_SUM,         &
                                 MPI_MODE_DELETE_ON_CLOSE
#endif
    ! 
    !   
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios 
    !! io status
    LOGICAL :: exst
    !! Find if a file exists.
    INTEGER, ALLOCATABLE :: indxp_tmp(:, :)
    !! 
    INTEGER, ALLOCATABLE :: indxm_tmp(:, :)
    !!     
    INTEGER :: i
    !! direction index in x axis
    INTEGER :: j
    !! direction index in y axis
    INTEGER :: k
    !! direction index in z axis
    INTEGER :: ii 
    !!
    INTEGER :: iq1
    !! q-index for 1st phonon
    INTEGER :: iq2
    !! q-index for 2nd phonon
    INTEGER :: iq3p 
    !! q-index for 3rd phonon
    INTEGER :: iq3m 
    !! q-index for 3rd phonon
    INTEGER :: nu1
    !! branch index for 1st phonon
    INTEGER :: nu2 
    !! branch index for 2nd phonon
    INTEGER :: nu3 
    !! branch index for 3rd phonon
    INTEGER :: iqx 
    !! direction index in x axis for phonon
    INTEGER :: iqy
    !! direction index in y axis for phonon
    INTEGER :: iqz 
    !! direction index in z axis for phonon
    INTEGER :: ipool 
    !! index for cpu
    INTEGER :: idim
    !! index for dimension
    INTEGER(KIND = 8) :: indx3p_p(npool) 
    !! Number of ph-ph matrix elements that needs to be
    INTEGER(KIND = 8) :: indx3p_m(npool) 
    !! Number of ph-ph matrix elements that needs to be
    !! calculated in each pool
    INTEGER(KIND = 8) :: ind1
    !! Counter on phonon modes
    INTEGER(KIND = 8) :: ind2
    !! Counter on phonon modes
    INTEGER(KIND = 8) :: indxall_p
    !! Total number of ph-ph matrix elements
    INTEGER(KIND = 8) :: indxall_m
    !! Total number of ph-ph matrix elements
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND) :: indxall_p_
    !! Total number of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: indxall_m_
    !! Total number of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: displ1
    !! displacements in order to balance 
    !! the ph-ph matrix elements between the pools
    INTEGER(KIND = MPI_OFFSET_KIND) :: displ2
    !! displacements in order to balance 
    !! the ph-ph matrix elements between the pools
    INTEGER(KIND = MPI_OFFSET_KIND) :: nindx
    !! number of elements written to file
    INTEGER(KIND = MPI_OFFSET_KIND) :: nindxi 
    !! number of elements read from file
    INTEGER(KIND = MPI_OFFSET_KIND) :: lowerplus
    !! lower bound for parallelization of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: upperplus
    !! upper bound for parallelization of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: lowerminus
    !! lower bound for parallelization of ph-ph matrix elements
    INTEGER(KIND = MPI_OFFSET_KIND) :: upperminus
    !! upper bound for parallelization of ph-ph matrix elements
#endif
    INTEGER :: lower_qbnd, upper_qbnd
    !! lower/upper bound for q parallelization
    REAL(KIND = DP) :: wqnu1
    !! frequency of 1st phonon
    REAL(KIND = DP) :: wqnu2 
    !! frequency of 2nd phonon
    REAL(KIND = DP) :: wqnu3p
    !! frequency of 3rd phonon
    REAL(KIND = DP) :: wqnu3m
    !! frequency of 3rd phonon
    REAL(KIND = DP) :: delta
    !! delta function
    REAL(KIND = DP) :: etap
    !! smearing parameter
    REAL(KIND = DP) :: etam 
    !! smearing parameter
    REAL(KIND = DP) :: q1mq2(3)
    !! q3 = q1 - q2
    REAL(KIND = DP) :: q1pq2(3)
    !! q3 = q1 + q2
    REAL(KIND = DP) :: eta1max 
    !! smearing parameter
    REAL(KIND = DP) :: eta1min
    !! smearing parameter
    REAL(KIND = DP) :: eta2max
    !! smearing parameter
    REAL(KIND = DP) :: eta2min
    !! smearing parameter
    REAL(KIND = DP) :: psi2
    !! ph-ph matrix element
    REAL(KIND = DP) :: w_sum
    !! sum of three phonon frequency
    REAL(KIND = DP) :: inv_eta
    !! inverse of the broadening
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    CHARACTER(LEN = 256) :: filename
    !! filename for mpi_write
    CHARACTER(LEN = 64) :: idim_ch
    !! Name of the files 
    REAL(KIND = DP) :: topsm1
    !! from ps^-1 to the inverse of Rydberg unit of time 
    !
    topsm1 = (ryd2ev * ps2sec) / hbar
    dt_in_ps = dph_tdbe * dt_tdbe / 1000.d0 
    !! dt_tdbe in fs, so devided by 1000 to ps
    !
    CALL start_clock('3-phonon')
    !
    WRITE(stdout, '(5x, a/)') "Phonon-phonon interaction is included, calculate ph-ph coupling matrix elements."
    WRITE(stdout, '(5x, a)') REPEAT('=',67)    
    IF (dg_tdbe)    WRITE(stdout, '(/5x, a)') 'Warning : double grid for phonon-phonon &
                                               interaction is not available, use the   & 
                                               grid nqf1, nqf2 and nqf3 to calculate   &
                                               phonon_phonon interaction' 
    ! 
    WRITE(stdout, '(/5x,a)') "Start running psi2_calc"
    CALL read_ifc3()
    WRITE(stdout, '(/5x,a)') '3rd order force constants have been read from FORCE_CONSTANTS_3RD'
    ! WRITE(stdout, *) 'ityp = ',ityp
    ! WRITE(stdout, *) 'amass = ', amass
    ! WRITE(stdout, *) 'nmodes=',nmodes,"nqtotf=",nqtotf
    CALL fkbounds(nqtotf, lower_qbnd, upper_qbnd)
    ! 
    indxall_p = 0
    indxall_m = 0
    indx3p_p = 0 
    indx3p_m = 0
    !
    !! Here we use atomic unit
    IF (adapt_smearing) THEN
      WRITE(stdout,'(/5x,a)') 'Use adaptive smearing in 3-phonon process'
      eta1min = 1000.0_dp
      eta2min = 1000.0_dp
      eta1max = 0.0_dp
      eta2max = 0.0_dp
    ENDIF
    !
    phwmin_tdbe = phwmin_tdbe / ryd2mev ! to Ryd 
    ! 
    WRITE(stdout, '(/5x, a)') REPEAT('=',67) 
    WRITE(stdout, '(/5x, a/)') 'Step 1: Now count the number of allowed 3-phonon scattering processes'
    DO ii = 1, 2
      !! ii = 1 : count and allocate
      !! ii = 2 : save the index while counting
      ! IF (ii == 2 ) THEN
        ind1 = 0
        ind2 = 0
      ! ENDIF !
      DO iq1 = lower_qbnd, upper_qbnd
        DO nu1 = 1, nmodes
          wqnu1 = wf(nu1, iq1)
          IF (wqnu1 < phwmin_tdbe) CYCLE
          DO iq2 = 1, nqtotf
            DO nu2 = 1, nmodes
              wqnu2 = wf(nu2, iq2)
              IF (wqnu2 < phwmin_tdbe) CYCLE
              q1pq2(:) = xqf(:, iq1) + xqf(:, iq2)
              iqx = MODULO(INT(q1pq2(1) * nqf1), nqf1) + 1
              iqy = MODULO(INT(q1pq2(2) * nqf2), nqf2) + 1
              iqz = MODULO(INT(q1pq2(3) * nqf3), nqf3) + 1
              iq3p = (iqx - 1) * nqf2 * nqf3 + &
                     (iqy - 1) * nqf3 + iqz
              q1mq2(:) = xqf(:, iq1) - xqf(:, iq2)
              iqx = MODULO(INT(q1mq2(1) * nqf1), nqf1) + 1
              iqy = MODULO(INT(q1mq2(2) * nqf2), nqf2) + 1
              iqz = MODULO(INT(q1mq2(3) * nqf3), nqf3) + 1
              iq3m = (iqx - 1) * nqf2 * nqf3 + &
                     (iqy - 1) * nqf3 + iqz
              DO nu3 = 1, nmodes
                wqnu3p = wf(nu3, iq3p)
                wqnu3m = wf(nu3, iq3m)
                IF (adapt_smearing) THEN
                  etap =  eta(&
                            vnuq_all(:, nu2, iq2) -&
                            vnuq_all(:, nu3, iq3p))
                  etam =  eta(&
                            vnuq_all(:, nu2, iq2)-&
                            vnuq_all(:, nu3, iq3m))
                  IF (ii == 1) THEN
                    IF (etap < eta1min) eta1min = etap
                    IF (etap > eta1max) eta1max = etap
                    IF (etam < eta2min) eta2min = etam
                    IF (etam > eta2max) eta2max = etam
                  ENDIF ! ii < 2
                ElSE
                  etap = degaussq
                  etam = degaussq
                ENDIF ! adapt_smearing
                IF (ABS(wqnu1 + wqnu2 - wqnu3p) <= (2.d0 * etap).AND. &
                  wqnu3p > phwmin_tdbe)  THEN
                  IF (ii == 1) THEN
                    indx3p_p(my_pool_id + 1)  = indx3p_p(my_pool_id + 1) + 1
                  ELSE ! ii == 2
                    ind1 = ind1 + 1
                    indxp_tmp(ind1, 1) = iq1
                    indxp_tmp(ind1, 2) = nu1
                    indxp_tmp(ind1, 3) = iq2
                    indxp_tmp(ind1, 4) = nu2
                    indxp_tmp(ind1, 5) = iq3p
                    indxp_tmp(ind1, 6) = nu3
                  ENDIF !! ii == 1 or 2
                ENDIF  !  |wqnu1 + wqnu2 - wqnu3p| < 2 * etap                             
                !
                IF (ABS(wqnu1 - wqnu2 - wqnu3m) <= (2.d0 * etam).AND. &
                  wqnu3m > phwmin_tdbe) THEN 
                  IF (ii == 1) THEN
                    indx3p_m(my_pool_id + 1) = indx3p_m(my_pool_id + 1) + 1
                  ELSE ! ii == 2
                    ind2 = ind2 + 1
                    indxm_tmp(ind2, 1) = iq1
                    indxm_tmp(ind2, 2) = nu1
                    indxm_tmp(ind2, 3) = iq2
                    indxm_tmp(ind2, 4) = nu2
                    indxm_tmp(ind2, 5) = iq3m
                    indxm_tmp(ind2, 6) = nu3                  
                  ENDIF ! ii
                ENDIF ! |wqnu1 - wqnu2 - wqnu3m| < 2 * etam  
              ENDDO ! nu3
            ENDDO ! nu2
          ENDDO ! iq2
        ENDDO ! nu1
      ENDDO ! iq1
      !
      CALL mp_barrier(inter_pool_comm)
#if defined(__MPI)
      IF (ii == 1) THEN
        !! indx3p_p and indx3p_m are kind = 8, can be accepted by mpi_allreduce
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, indx3p_p, npool, MPI_INTEGER8, MPI_SUM, inter_pool_comm, ierr)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, indx3p_m, npool, MPI_INTEGER8, MPI_SUM, inter_pool_comm, ierr)
        !! indx3p_p and indx3p_m are kind = 8, can not be accepted by mp_sum
        ! CALL mp_sum(indx3p_p, inter_pool_comm)
        ! CALL mp_sum(indx3p_m, inter_pool_comm)
!        CALL MPI_ALLREDUCE(indx3p_p,  indx3p_p_all,  npool, MPI_INTEGER8, MPI_SUM, inter_pool_comm, ierr)
!        CALL MPI_ALLREDUCE(indx3p_m, indx3p_m_all, npool, MPI_INTEGER8, MPI_SUM, inter_pool_comm, ierr)
!        indx3p_p  = indx3p_p_all
!        indx3p_m = indx3p_m_all
!         WRITE(stdout,*) 'indx3p_p = ', indx3p_p(:)
!         WRITE(stdout,*) 'indx3p_m = ', indx3p_m(:)
        CALL mp_barrier(inter_pool_comm)
        IF (adapt_smearing) THEN
          CALL mp_max(eta1max, inter_pool_comm)
          CALL mp_min(eta1min, inter_pool_comm)
          CALL mp_max(eta2max, inter_pool_comm)
          CALL mp_min(eta2min, inter_pool_comm)
          WRITE(stdout,'(5x,a,ES20.10,a)') 'eta1max = ', eta1max * rydcm1, ' cm^-1'
          WRITE(stdout,'(5x,a,ES20.10,a)') 'eta2max = ', eta2max * rydcm1, ' cm^-1'
          WRITE(stdout,'(5x,a,ES20.10,a)') 'eta1min = ', eta1min * rydcm1, ' cm^-1'
          WRITE(stdout,'(5x,a,ES20.10,a)') 'eta2min = ', eta2min * rydcm1, ' cm^-1'
        ENDIF
        !
        indxall_p = SUM(indx3p_p)
        indxall_m = SUM(indx3p_m)
        CALL mp_barrier(inter_pool_comm)
        !! indxall_p_ and indxall_p_ are MPI_OFFSET_KIND
        indxall_p_ = INT(indxall_p, KIND = MPI_OFFSET_KIND)
        indxall_m_ = INT(indxall_m, KIND = MPI_OFFSET_KIND)
        WRITE(stdout,'(5x,a,i20)')  "indxall_p = ", indxall_p
        WRITE(stdout,'(5x,a,i20)')  "indxall_m = ", indxall_m 
        WRITE(stdout,'(/5x,a/)') 'Now find the indexes of three phonon scattering '
        !WRITE(stdout,*) 'Check the memory'
        ALLOCATE(indxp_tmp(indx3p_p(my_pool_id + 1), 6), STAT = ierr)
        IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating indxp_tmp', 1)
        ALLOCATE(indxm_tmp(indx3p_m(my_pool_id + 1),6), STAT = ierr)
        IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating indxm_tmp', 1)
        WRITE(stdout,'(5x,a)') 'Step 2: indexing according to the 1st q point'
        !
        indxp_tmp = 0
        indxm_tmp = 0
        ! ind1 = 0
        ! ind2 = 0
        !
        CALL mp_barrier(inter_pool_comm)
      ENDIF ! ii == 1
    ENDDO ! ii
    !
    WRITE(stdout, '(/5x,a/)') 'Step 3: distribute 3-phonon indxes'
    !
    CALL fkbounds2(indxall_p_,  lowerplus, upperplus)
    CALL fkbounds2(indxall_m_, lowerminus, upperminus)
    !   
    nind_p = upperplus - lowerplus + 1
    nind_m = upperminus - lowerminus + 1
    !
    ! WRITE(stdout, *) 'nind_p =', nind_p
    ! WRITE(stdout, *) 'nind_m =', nind_m
    ALLOCATE(ind_p(nind_p, 6), STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating ind_p', 1)
    ALLOCATE(ind_m(nind_m, 6), STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating ind_m', 1)
    !
    ind_p = 0
    ind_m = 0
    !
    DO ii = 1, 2 
      !! ii == 1 : redistribute indxp_tmp
      !! ii == 2 : redistribute indxm_tmp
      IF (ii == 1) THEN
        nindx =  INT(indx3p_p(my_pool_id + 1), KIND = MPI_OFFSET_KIND)
        nindxi = INT(nind_p, KIND = MPI_OFFSET_KIND)
      ELSE
        nindx =  INT(indx3p_m(my_pool_id + 1), KIND = MPI_OFFSET_KIND)
        nindxi = INT(nind_m, KIND = MPI_OFFSET_KIND)
      ENDIF !! ii == 1 or 2
      CALL mp_barrier(inter_pool_comm)
      !
      displ1 = 0_MPI_OFFSET_KIND
      displ2 = 0_MPI_OFFSET_KIND
      !
      IF (ii == 1) THEN
        IF (my_pool_id > 0) THEN
          DO ipool = 1, my_pool_id
            displ1 = displ1 + INT(indx3p_p(ipool), KIND = MPI_OFFSET_KIND)
          ENDDO
        ENDIF
        !
        displ1 = displ1 * 4_MPI_OFFSET_KIND
        displ2 = (lowerplus - 1_MPI_OFFSET_KIND) * 4_MPI_OFFSET_KIND
        !
        WRITE(stdout,'(7x,a)') 'Now distribute indx_plus  evenly across the pools'
      ELSE ! ii == 2
        IF (my_pool_id > 0) THEN
          DO ipool = 1, my_pool_id
            displ1 = displ1 + INT(indx3p_m(ipool), KIND = MPI_OFFSET_KIND)
          ENDDO
        ENDIF
        !
        displ1 = displ1 * 4_MPI_OFFSET_KIND
        displ2 = (lowerminus - 1_MPI_OFFSET_KIND) * 4_MPI_OFFSET_KIND    
        !
        WRITE(stdout,'(7x,a)') 'Now distribute indx_minus evenly across the pools'  
      ENDIF
      CALL mp_barrier(inter_pool_comm)
      !
      IF (ii == 1) THEN    
        DO idim = 1, 6
          WRITE(idim_ch, "(I0)") idim
          filename = 'indph_p'// TRIM(idim_ch)
          ! open file
          CALL MPI_FILE_OPEN(inter_pool_comm, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_OPEN', 1)
          ! Write data to file
          CALL MPI_FILE_WRITE_AT(iunphindx, displ1, indxp_tmp(:,idim), nindx, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_WRITE_AT', 1)
          ! Close file
          CALL MPI_FILE_CLOSE(iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_CLOSE', 1)
          ! Wait for all processes to finish writing
          CALL mp_barrier(inter_pool_comm)
          ! Open file for reading
          CALL MPI_FILE_OPEN(inter_pool_comm, filename, &
               MPI_MODE_RDONLY + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_OPEN', 1)
          ! Read data from file
          CALL MPI_FILE_READ_AT(iunphindx, displ2, ind_p(:,idim), nindxi, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_READ_AT', 1)
          ! Close file
          CALL MPI_FILE_CLOSE(iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_CLOSE', 1)
          !
          CALL mp_barrier(inter_pool_comm)
        ENDDO  ! idim
      ELSE !! ii == 2
        DO idim = 1, 6
          WRITE(idim_ch, "(I0)") idim
          filename = 'indph_m'// TRIM(idim_ch)
          ! open file
          CALL MPI_FILE_OPEN(inter_pool_comm, filename, &
               MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_OPEN', 1)
          ! Write data to file
          CALL MPI_FILE_WRITE_AT(iunphindx, displ1, indxm_tmp(:,idim), nindx, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_WRITE_AT', 1)
          ! Close file
          CALL MPI_FILE_CLOSE(iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_CLOSE', 1)
          ! Wait for all processes to finish writing
          CALL mp_barrier(inter_pool_comm)
          ! Open file for reading
          CALL MPI_FILE_OPEN(inter_pool_comm, filename, MPI_MODE_RDONLY + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_OPEN', 1)
          ! Read data from file
          CALL MPI_FILE_READ_AT(iunphindx, displ2, ind_m(:,idim), nindxi, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_READ_AT', 1)
          ! Close file
          CALL MPI_FILE_CLOSE(iunphindx, ierr)
          IF (ierr /= 0) CALL errore('psi2_calc', 'error in MPI_FILE_CLOSE', 1)
          !
          CALL mp_barrier(inter_pool_comm)
        ENDDO  ! idim
      ENDIF
#endif
    ENDDO !!   
    !
    ! WRITE(stdout, *) 'ind'
    ! DO idim = 1, 6
    ! WRITE(stdout, *) ind_m(10, idim)
    ! WRITE(stdout, *) ind_p(10, idim)
    ! ENDDO
    ! WRITE(stdout, *) 'indxp'
    ! DO idim = 1, 6
    ! WRITE(stdout, *) indxm_tmp(10, idim)
    ! WRITE(stdout, *) indxp_tmp(10, idim)
    ! ENDDO
    DEALLOCATE(indxp_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error deallocating indxp_tmp', 1)
    DEALLOCATE(indxm_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error deallocating indxm_tmp', 1)
    !
    ALLOCATE(psi2_p(nind_p), STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating psi2_p', 1)
    ALLOCATE(psi2_m(nind_m), STAT = ierr)
    IF (ierr /= 0) CALL errore('psi2_calc', 'Error allocating psi2_m', 1)
    psi2_p  = 0.d0
    psi2_m = 0.d0
    WRITE(stdout,'(/5x, a/)') 'Step 4: Calculate ph-ph matrix elements.'
    DO ind1 = 1, nind_p
      iq1  = ind_p(ind1, 1)
      nu1  = ind_p(ind1, 2)
      iq2  = ind_p(ind1, 3)
      nu2  = ind_p(ind1, 4)
      iq3p = ind_p(ind1, 5)
      nu3  = ind_p(ind1, 6)
      CALL psi_interp(nu1, nu2, nu3, iq1, iq2, iq3p, 1.0_dp, psi2)
      IF (adapt_smearing) THEN
        etap =  eta(          &
                  vnuq_all(:, nu2, iq2) -          &
                  vnuq_all(:, nu3, iq3p))
      ElSE
        etap = degaussq
      ENDIF    
      wqnu1  = wf(nu1, iq1)
      wqnu2  = wf(nu2, iq2)
      wqnu3p = wf(nu3, iq3p)
      inv_eta = 1.d0 / etap
      w_sum = wqnu1 + wqnu2 - wqnu3p
      delta = w0gauss(w_sum * inv_eta, 0) *  inv_eta
      psi2_p(ind1) = pi * psi2 * delta / nqtotf 
      psi2_p(ind1) = psi2_p(ind1) * topsm1
      !! Now the unit of psi2_p is ps^-1
    ENDDO
    DO ind2 = 1, nind_m
      iq1  = ind_m(ind2, 1)
      nu1  = ind_m(ind2, 2)
      iq2  = ind_m(ind2, 3)
      nu2  = ind_m(ind2, 4)
      iq3m = ind_m(ind2, 5)
      nu3  = ind_m(ind2, 6)
      CALL psi_interp(nu1, nu2, nu3, iq1, iq2, iq3m, -1.0_dp, psi2)
      IF (adapt_smearing) THEN
        etam=eta(&
               vnuq_all(:, nu2, iq2)-&
               vnuq_all(:, nu3, iq3m))
      ElSE
        etam = degaussq
      ENDIF 
      !
      wqnu1  = wf(nu1, iq1)
      wqnu2  = wf(nu2, iq2)
      wqnu3m = wf(nu3, iq3m)
      inv_eta = 1.d0 / etam
      w_sum = wqnu1 - wqnu2 - wqnu3m
      delta = w0gauss(w_sum * inv_eta, 0) *  inv_eta
      psi2_m(ind2) = pi  * psi2 * delta / nqtotf
      psi2_m(ind2) = psi2_m(ind2) * topsm1
      !! Now the unit of psi2_m is ps^-1
    ENDDO
    !
    ! WRITE(stdout, *) psi2_m(1:10)
    ! WRITE(stdout, *) psi2_p(1:10)
    WRITE(stdout, '(/5x, a)') REPEAT('=',67)
    WRITE(stdout, '(/5x, a/)') '3-phonon scattering matrix elements have been calculated'
    WRITE(stdout, '(5x, a)') REPEAT('=',67)
    CALL mp_barrier(inter_pool_comm)
    !
    RETURN
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE psi2_calc
    !--------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------- 
    SUBROUTINE iphph_calc()
    !-------------------------------------------------------------------------- 
    !! Calculate phonon-phonon collision integral. More details are given by 
    !! the Eq. (S3) of J. Phys. Chem. Lett. 12, 1734.
    !-------------------------------------------------------------------------- 
    !     
    USE kinds,            ONLY : DP
    USE global_var,       ONLY : nqtotf
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE modes,            ONLY : nmodes
    USE tdbe_common,      ONLY : nphon_pre, iph_ph, a_tab, b_tab,c_tab, &
                                 psi2_p, psi2_m, ind_p, ind_m, nstg,    &
                                 dt_in_ps, nind_p, nind_m                       
    !     
    IMPLICIT NONE
    ! 
    INTEGER :: iq1, iq2, iq3
    !! q-indexes of the 3 phonons involved in the ph-ph interaction
    INTEGER :: nu1, nu2, nu3
    !! phonon branch indexes of the 3 phonons involved in the ph-ph interaction
    INTEGER :: istg, pp
    !! Counter on the stages of runge-kutta solver
    INTEGER(KIND = 8) :: indx
    !! Counter on ph-ph coupling matrix element
    REAL(KIND = DP) :: nph1, nph2, nph3
    !! Nonequilibrium phonon distribution of the three phonons involved
    REAL(KIND = DP) :: facp
    !! factor arise from products of phonon distribution
    ! 
    iph_ph = 0.d0
    DO istg = 1, nstg
      DO indx = 1, nind_p
        iq1 = ind_p(indx, 1)
        nu1 = ind_p(indx, 2)
        iq2 = ind_p(indx, 3)
        nu2 = ind_p(indx, 4)
        iq3 = ind_p(indx, 5)
        nu3 = ind_p(indx, 6)
        nph1 = nphon_pre(nu1, iq1)
        nph2 = nphon_pre(nu2, iq2)
        nph3 = nphon_pre(nu3, iq3)
        IF (istg > 1) THEN
          DO pp = 1, istg -1
            nph1 = nph1 + a_tab(istg, pp) * iph_ph(nu1, iq1, pp) * dt_in_ps
            nph2 = nph2 + a_tab(istg, pp) * iph_ph(nu2, iq2, pp) * dt_in_ps
            nph3 = nph3 + a_tab(istg, pp) * iph_ph(nu3, iq3, pp) * dt_in_ps
          ENDDO ! pp
        ENDIF ! istg
!! In Eq. (S3) facp = (1.d0 + nph1) * (1.d0 + nph2) * &
!!      nph3 - nph1 * nph2* (1.d0 + nph3)
!!      It is simplified to the expression below
        facp = - nph1 * (nph2 - nph3) + nph3 * (nph2 + 1.d0)
        iph_ph(nu1, iq1, istg) = iph_ph(nu1, iq1, istg) + psi2_p(indx) * facp
      ENDDO
      !
      DO indx = 1, nind_m
        iq1 = ind_m(indx, 1)
        nu1 = ind_m(indx, 2)
        iq2 = ind_m(indx, 3)
        nu2 = ind_m(indx, 4)
        iq3 = ind_m(indx, 5)
        nu3 = ind_m(indx, 6)
        nph1 = nphon_pre(nu1, iq1)
        nph2 = nphon_pre(nu2, iq2)
        nph3 = nphon_pre(nu3, iq3)
        IF (istg > 1) THEN
          DO pp = 1, istg -1
            nph1 = nph1 + a_tab(istg, pp) * iph_ph(nu1, iq1, pp) * dt_in_ps
            nph2 = nph2 + a_tab(istg, pp) * iph_ph(nu2, iq2, pp) * dt_in_ps
            nph3 = nph3 + a_tab(istg, pp) * iph_ph(nu3, iq3, pp) * dt_in_ps
          ENDDO ! pp
        ENDIF ! istg
!! In Eq. (S3) facp = (1.d0 + nph1) * nph2 * &
!!      nph3 - nph1 * (1.d0 + nph2) * (1.d0 + nph3)
!!      It is simplified to the expression below
        facp = nph2 * nph3 - nph1 * (nph2 + nph3 + 1.d0)
        iph_ph(nu1, iq1, istg) = iph_ph(nu1, iq1, istg) + psi2_m(indx) * facp * 0.5d0
      ENDDO
      CALL mp_sum(iph_ph(:, :, istg), inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
    ENDDO ! istg
    !
    RETURN
    !
    !-------------------------------------------------------------------------- 
    END SUBROUTINE iphph_calc
    !-------------------------------------------------------------------------- 
    !  
    !-------------------------------------------------------------------------- 
    SUBROUTINE inv_tau_rta(temp_min, temp_max, ntemps)
    !-------------------------------------------------------------------------- 
    !! calculate inverse phonon lifetime (linewidth) under RTA, the details of the
    !! expression can be found from many literatures, for example, Eq.(6-8) of 
    !! Comp. Phys. Commun. 185, 17471758 (2014).
    !--------------------------------------------------------------------------     
    ! 
    USE kinds,            ONLY : DP
    USE mp_world,         ONLY : mpime
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE global_var,       ONLY : nqtotf, wf
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE modes,            ONLY : nmodes
    USE tdbe_common,      ONLY : ind_p, ind_m, nind_p, nind_m,      &
                                 psi2_p, psi2_m
    USE ep_constants,     ONLY : kelvin2Ry
    !     
    IMPLICIT NONE
    !
    REAL(KIND=DP),INTENT(in)  :: temp_min, temp_max
    !! lower and upper bound for temperature, in Kelvin
    INTEGER ,INTENT(in) :: ntemps
    !! Number of point in the interval [tmin, tmax]
    REAL(KIND = DP)  :: inv_tau_ph(ntemps, nqtotf,nmodes)
    !! inverse ph-ph lifetime
    REAL(KIND = DP) :: nph2, nph3
    !! Thermal phonon population 
    REAL(KIND = DP) :: temp_in_ry  
    !! Temperature in kelvin (k_B * T)
    INTEGER :: iq1, iq2, iq3
    !! q index for the 3 phonons involved
    INTEGER :: nu1, nu2, nu3
    !! Phonon branch indexes for the 3 phonons involved
    INTEGER(KIND = 8) :: indx
    !! Counter on ph-ph matrix elements
    INTEGER :: itemp
    !! Counter on temperature
    character(len=4) :: temp_in_ch
    !! temperature in string
    character(len=32) ::fname
    !! file name 
    !
    inv_tau_ph = 0.d0
    DO itemp = 1, ntemps
      temp_in_ry = (temp_min + (temp_max - temp_min) / (ntemps - 1)  &
                 * (itemp - 1)) * kelvin2Ry ! 6.33363E-6 is kelvin2Ry
      DO indx = 1, nind_p
        iq1 = ind_p(indx, 1)
        nu1 = ind_p(indx, 2)
        iq2 = ind_p(indx, 3)
        nu2 = ind_p(indx, 4)
        iq3 = ind_p(indx, 5)
        nu3 = ind_p(indx, 6)
        nph2 = 1.0 / (EXP(wf(nu2, iq2) / temp_in_ry) - 1.0)
        ! Bose-Einstein distribution
        nph3 = 1.0 / (EXP(wf(nu3, iq3) / temp_in_ry) - 1.0)
        ! Bose-Einstein distribution
        inv_tau_ph(itemp, iq1, nu1) = inv_tau_ph(itemp, iq1, nu1) +  &
                                      psi2_p(indx) * (nph2 - nph3)
      ENDDO
      !
      DO indx = 1, nind_m
        iq1 = ind_m(indx, 1)
        nu1 = ind_m(indx, 2)
        iq2 = ind_m(indx, 3)
        nu2 = ind_m(indx, 4)
        iq3 = ind_m(indx, 5)
        nu3 = ind_m(indx, 6)
        nph2 = 1.0 / (EXP(wf(nu2, iq2) / temp_in_ry) - 1.0)
        ! Bose-Einstein distribution
        nph3 = 1.0 / (EXP(wf(nu3, iq3) / temp_in_ry) - 1.0)
        ! Bose-Einstein distribution
        inv_tau_ph(itemp, iq1, nu1) = inv_tau_ph(itemp, iq1, nu1) +  &
                  psi2_m(indx) * (nph2 + nph3 + 1.0) / 2.0
      ENDDO
      !
    ENDDO
    CALL mp_sum(inv_tau_ph, inter_pool_comm)
    WRITE(stdout, '(/5x,a/)') 'print the inverse life time of phonon due to phonon-phonon interaction.'
    IF (my_pool_id == ionode_id) THEN
      DO itemp = 1, ntemps
        temp_in_ry = (temp_min + (temp_max - temp_min)/(ntemps - 1) * (itemp - 1))
        WRITE(temp_in_ch, "(I0)") NINT(temp_in_ry)
        fname = "inv_tau_rta_" // TRIM(adjustl(temp_in_ch)) // "K"
        OPEN(unit = 222, FILE = fname)
        DO iq1 = 1, nqtotf
          DO nu1 = 1, nmodes
            WRITE(222, '(2E22.14)') wf(nu1, iq1), inv_tau_ph(itemp, iq1, nu1)
          ENDDO
        ENDDO
        CLOSE(222)
      ENDDO
    ENDIF
    !
    RETURN
    !
    !--------------------------------------------------------------------------------
    END SUBROUTINE inv_tau_rta
    !--------------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------------
    SUBROUTINE read_ifc3()
    !--------------------------------------------------------------------------------
    !! Read 3rd order force constants, generated by thirdorder.py
    !--------------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id    
    USE mp,               ONLY : mp_barrier, mp_bcast
    USE tdbe_common,      ONLY : nifc, ifc3, rvecs, ind_atms
    USE io_var,           ONLY : iun3rdfc
    USE ep_constants,     ONLY : bohr2ang, ryd2ev
    USE cell_base,        ONLY : alat, bg
    !
    implicit none
    !
    INTEGER :: i, j, k
    !! direction index
    INTEGER :: m 
    !!
    INTEGER :: ibid, jbid, kbid
    !! 
    INTEGER :: fcbid
    !!
    INTEGER :: fc 
    !! counter of the fc
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! I/O status 
    !
    IF (my_pool_id == ionode_id) THEN
      OPEN(UNIT = iun3rdfc, FILE = 'FORCE_CONSTANTS_3RD', status="old", IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_ifc3', 'error opening file FORCE_CONSTANTS_3RD', 1)
      READ(iun3rdfc, *) nifc
      ALLOCATE(ind_atms(3, nifc), STAT =ierr)
      IF (ierr /= 0) CALL errore('read_ifc3', 'Error allocating ind_atms', 1)
      ALLOCATE(ifc3(3, 3, 3, nifc), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_ifc3', 'Error allocating ifc3', 1)
      ALLOCATE(rvecs(2, 3, nifc), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_ifc3', 'Error allocating rvecs', 1)
      DO fc = 1, nifc
        READ(iun3rdfc, *) fcbid
        IF (fcbid /= fc) CALL errore('read_ifc3', 'Error in reading ifc3', 1)
        READ(iun3rdfc, *) rvecs(1, :, fc)
        READ(iun3rdfc, *) rvecs(2, :, fc)
        READ(iun3rdfc, *) ind_atms(:, fc)
        READ(iun3rdfc, *) (((ibid, jbid, kbid, ifc3(i, j, k, fc), k = 1, 3), &
                           j = 1, 3), i = 1, 3)
      ENDDO
      CLOSE(iun3rdfc)
      !
      rvecs = rvecs / bohr2ang / alat
      ifc3 = ifc3 * bohr2ang**3 / ryd2ev  ! transform to Rydberg atomic unit
      !! Convert back to crystal coordinates
      CALL cryst_to_cart (nifc, rvecs(1, :, :), bg, -1)
      CALL cryst_to_cart (nifc, rvecs(2, :, :), bg, -1)
      !
    ENDIF ! my_pool_id == ionode_id
    CALL mp_bcast(nifc, ionode_id, inter_pool_comm)
    IF (my_pool_id /= ionode_id) THEN
        ALLOCATE(ind_atms(3, nifc), STAT =ierr)
        IF (ierr /= 0) CALL errore('read_ifc3', 'Error allocating ind_atms', 1)
        ALLOCATE(ifc3(3, 3, 3, nifc), STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc3', 'Error allocating ifc3', 1)
        ALLOCATE(rvecs(2, 3, nifc), STAT = ierr)
        IF (ierr /= 0) CALL errore('read_ifc3', 'Error allocating rvecs', 1)
    ENDIF
    CALL mp_bcast(ind_atms, ionode_id, inter_pool_comm)
    CALL mp_bcast(ifc3,     ionode_id, inter_pool_comm)
    CALL mp_bcast(rvecs,    ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    RETURN
    !
    !--------------------------------------------------------------------------------        
    END SUBROUTINE read_ifc3
    !--------------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------------    
    SUBROUTINE deallocate_phph()
    !--------------------------------------------------------------------------------  
    !! deallocate the 3-phonon matrix elements and their indxes.
    !-------------------------------------------------------------------------------- 
    !
    USE tdbe_common,    ONLY : ind_atms, ifc3, rvecs, psi2_p, &
                               psi2_m, ind_p, ind_m, ph_lw
    !/
    IMPLICIT NONE 
    !
    INTEGER :: ierr
    !
    DEALLOCATE(ind_atms, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ind_atms', 1)
    DEALLOCATE(ifc3,    STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ifc3', 1)
    DEALLOCATE(rvecs,    STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating rvecs', 1)
    DEALLOCATE(psi2_p, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating psi2_p', 1)
    DEALLOCATE(psi2_m, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating psi2_m', 1)
    DEALLOCATE(ind_p, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ind_p', 1)
    DEALLOCATE(ind_m, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_phph', 'Error deallocating ind_m', 1)
    !
    !--------------------------------------------------------------------------------  
    END SUBROUTINE deallocate_phph
    !-------------------------------------------------------------------------------- 
    !    
    !--------------------------------------------------------------------------------  
    SUBROUTINE psi_interp(nu1, nu2, nu3, iq1, iq2, iq3, pm, psi2)
    !-------------------------------------------------------------------------------- 
    !! Calculate the modulus square of ph-ph matrix element from 3rd force constants
    !! Expressions can be found in Eq.(9) of Comp. Phys. Commun. 185, 17471758 (2014)
    !-------------------------------------------------------------------------------- 
    USE kinds,            ONLY : DP 
    USE modes,            ONLY : nmodes
    USE tdbe_common,      ONLY : uf_all, ifc3, rvecs, ind_atms, nifc
    USE global_var,       ONLY : wf, xqf
    USE ions_base,        ONLY : nat, amass, ityp, tau
    USE ep_constants,     ONLY : twopi, ci, eps4
    !
    IMPLICIT NONE 
    !   
    INTEGER, INTENT(IN) :: nu1
    !! branch index of 1st phonon
    INTEGER, INTENT(IN) :: nu2
    !! branch index of 2nd phonon
    INTEGER, INTENT(IN) :: nu3
    !! branch index of 3rd phonon
    INTEGER, INTENT(IN) :: iq1
    !! q-index for 1st phonon
    INTEGER, INTENT(IN) :: iq2
    !! q-index for 2nd phonon
    INTEGER, INTENT(IN) :: iq3 
    !! q-index for 3rd phonon
    REAL(KIND = DP), INTENT(IN) :: pm
    !! plus or minus
    REAL(KIND = DP), INTENT(OUT) :: psi2
    !! results of interpolation
    INTEGER :: fc
    !! index of force constants
    INTEGER :: i
    !! direction index
    INTEGER :: j
    !! direction index
    INTEGER :: k
    !! direction index
    REAL(KIND = DP) :: wqnu1
    !! frequency of 1st phonon
    REAL(KIND = DP) :: wqnu2 
    !! frequency of 2nd phonon
    REAL(KIND = DP) :: wqnu3
    !! frequency of 3rd phonon
    REAL(KIND = DP) :: xxq2(3)
    !! coordinate  of 2nd phonon
    REAL(KIND = DP) :: xxq3(3)
    !! coordinate  of 3rd phonon
    REAL(KIND = DP) :: xxr2(3)
    !! coordinate  of 2nd phonon
    REAL(KIND = DP) :: xxr3(3)
    !! coordinate  of 3rd phonon
    COMPLEX(KIND = DP) :: vectq1(nmodes, nmodes)
    !! eigenmode 1
    COMPLEX(KIND = DP) :: vectq2(nmodes, nmodes)
    !! eigenmode 2
    COMPLEX(KIND = DP) :: vectq3(nmodes, nmodes)
    !! eigenmode 3
    COMPLEX(KIND = DP) :: psi_plus
    !! In case of plus
    COMPLEX(KIND = DP) :: psi_minus
    !! In case of minus
    COMPLEX(KIND = DP) :: factor, psi0     
    ! 
    xxq2  = xqf(:, iq2)
    xxq3 =  xqf(:, iq3)
    wqnu1  = wf(nu1, iq1)
    wqnu2  = wf(nu2, iq2)
    wqnu3  = wf(nu3, iq3)
    vectq1 = uf_all(iq1, :, :)
    vectq2 = uf_all(iq2, :, :)
    vectq3 = uf_all(iq3, :, :)
    psi_plus = 0.d0
    psi_minus = 0.d0   
    psi2  = 0.d0         
    IF (ABS(pm - 1.0) < eps4) THEN
      ! 
      DO fc = 1, nifc
        xxr2 = rvecs(1, :,fc)
        xxr3 = rvecs(2, :,fc)
        factor = 1.d0 / SQRT(amass(ityp(ind_atms(1, fc)))  *  &
                             amass(ityp(ind_atms(2, fc)))  *  &
                             amass(ityp(ind_atms(3, fc)))) *  &
                EXP( twopi * ci * DOT_PRODUCT(xxq2, xxr2)) * &
                EXP(-twopi * ci * DOT_PRODUCT(xxq3, xxr3))
        psi0 = CMPLX(0.d0, 0.d0)
        DO i = 1, 3
          DO j = 1, 3
            DO k = 1, 3
              psi0 = psi0 + ifc3(k, j, i, fc) * &
                    vectq1(k + 3 * (ind_atms(1, fc) - 1), nu1) * &
                    vectq2(j + 3 * (ind_atms(2, fc) - 1), nu2) * &
              CONJG(vectq3(i + 3 * (ind_atms(3, fc) - 1), nu3))
            END DO
          END DO
        END DO
        psi_plus = psi_plus + factor * psi0
      END DO
      psi2= ABS(psi_plus) ** 2 /(wqnu1 * wqnu2 * wqnu3)/4.d0
      !
    ELSEIF(ABS(pm + 1.0) < eps4) THEN
      ! 
      DO fc = 1, nifc
        xxr2 = rvecs(1, :, fc)
        xxr3 = rvecs(2, :, fc)
        factor = 1.d0 / SQRT(amass(ityp(ind_atms(1, fc)))  *  &
                             amass(ityp(ind_atms(2, fc)))  *  &
                             amass(ityp(ind_atms(3, fc)))) *  &
                EXP(-twopi * ci * DOT_PRODUCT(xxq2, xxr2)) * &
                EXP(-twopi * ci * DOT_PRODUCT(xxq3, xxr3))
        psi0 = CMPLX(0.d0, 0.d0)
        DO i = 1, 3
          DO j = 1, 3
            DO k = 1, 3
              psi0 = psi0 + ifc3(k, j, i, fc) * &
                    vectq1(k + 3 * (ind_atms(1, fc) - 1), nu1)  * &
              CONJG(vectq2(j + 3 * (ind_atms(2, fc) - 1), nu2)) * &
              CONJG(vectq3(i + 3 * (ind_atms(3, fc) - 1), nu3))
            END DO
          END DO
        END DO
        psi_minus = psi_minus + factor * psi0
      END DO
      psi2= ABS(psi_minus) ** 2 /(wqnu1 * wqnu2 * wqnu3)/4.d0
      !
    ELSE 
      CALL errore('psi_interp', 'pm must be 1.0 or -1.0', 1)
    ENDIF
    !
    RETURN
    !
    !-------------------------------------------------------------------------------- 
    END SUBROUTINE psi_interp
    !-------------------------------------------------------------------------------- 
    !
    !-------------------------------------------------------------------------------- 
    FUNCTION eta(v)
    !-------------------------------------------------------------------------------- 
    !! Return the base broadening (without prefac) for a mode.
    !-------------------------------------------------------------------------------- 
    USE kinds,            ONLY : DP 
    USE cell_base,         ONLY : bg, alat
    USE ep_constants,    ONLY : ryd2mev, twopi 
    USE input,           ONLY : nqf1, nqf2, nqf3
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(IN) :: v(3)
    REAL(KIND = DP):: eta
    REAL(KIND = DP) :: eta_tmp(3)
    REAL(KIND = DP) :: etamin
    INTEGER :: idim
    !
    etamin = 0.05 / ryd2mev 
    !
    eta=0.d0
    eta_tmp(1) = (twopi / alat) * ABS(DOT_PRODUCT(bg(:, 1), v))/ DBLE(nqf1)
    eta_tmp(2) = (twopi / alat) * ABS(DOT_PRODUCT(bg(:, 2), v))/ DBLE(nqf2)
    eta_tmp(3) = (twopi / alat) * ABS(DOT_PRODUCT(bg(:, 3), v))/ DBLE(nqf3)
    DO idim = 1, 3
      eta = eta + eta_tmp(idim) ** 2
    END DO
    eta = 0.5d0 * DSQRT(eta) / SQRT(12.0d0)
    !! if eta < 0.05 mev, set 0.05 mev
    IF (eta < etamin ) THEN
      eta = etamin
    ENDIF
    !
    RETURN
    !
    !-------------------------------------------------------------------------------- 
    END function eta
    !-------------------------------------------------------------------------------- 
    !
    !-------------------------------------------------------------------------------- 
    SUBROUTINE ph_energy()
    !-------------------------------------------------------------------------------- 
    !! Calculate a list of phonon energy for temperatues ranging from temp_ph_tdbe to 
    !! temp_el_tdbe so that one can use the energy as a reference to find out the 
    !! effective phonon temperature
    !-------------------------------------------------------------------------------- 
    !    
    USE kinds,            ONLY : DP 
    USE input,            ONLY : temp_el_tdbe, temp_ph_tdbe, eps_acoustic
    USE global_var,       ONLY : wf,nqtotf
    use modes,            ONLY : nmodes
    USE parallelism,      ONLY : fkbounds
    USE io_global,        ONLY : stdout
    USE ep_constants,     ONLY : kelvin2Ry
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id 
    USE tdbe_common,          ONLY : phtemp, e_latt, istart, ntph 
    !
    IMPLICIT NONE 
    !
    INTEGER :: ierr
    !! Error info
    INTEGER :: itemp
    !! Counter of temperature
    INTEGER :: iq 
    !! Counter of q points 
    INTEGER :: nu
    !! Counter of phonon branch
    INTEGER :: lower_qbnd
    !! lower q index
    INTEGER :: upper_qbnd
    !! upper q index
    REAL(KIND = DP) :: dtemp
    !! temperature difference
    REAL(KIND = DP) :: n_eq
    !! Equilibirum phonon population
    !
    WRITE(stdout,'(5x,a)') 'Phonon-phonon interaction is treated with RTA'
    WRITE(stdout,'(5x,a)') 'Phonon energy is evaluated on a list of temperature in the range of Tph,Tel'
    phtemp =0.d0
    e_latt = 0.d0
    dtemp = ABS(temp_el_tdbe - temp_ph_tdbe) / DBLE(ntph - 1)
    CALL fkbounds(nqtotf, lower_qbnd, upper_qbnd)
    DO itemp = 1, ntph
      phtemp(itemp) = temp_ph_tdbe + dtemp*(itemp -1) !! Kelvin
      DO iq = lower_qbnd, upper_qbnd
        DO nu = 1, nmodes
          IF (wf(nu,iq) > eps_acoustic) THEN
            n_eq = 1.0 / (EXP(wf(nu, iq) / (phtemp(itemp) * kelvin2Ry)) - 1.0)
            e_latt(itemp) = e_latt(itemp) + n_eq * wf(nu,iq) / nqtotf
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    CALL mp_sum(e_latt, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    istart = 1
    !
    RETURN
    !-------------------------------------------------------------------------------- 
    END SUBROUTINE ph_energy
    !-------------------------------------------------------------------------------- 
    !
    !-------------------------------------------------------------------------------- 
    SUBROUTINE read_ph_lw()
    !-------------------------------------------------------------------------------- 
    !! In case the collision integrals of ph-ph interaction are approximated with
    !! relaxation time approximation, the phonon linewidth can be read from external
    !! files.
    !-------------------------------------------------------------------------------- 
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    use modes,            ONLY : nmodes
    USE global_var,       ONLY : wf,nqtotf
    USE io_global,        ONLY : stdout, ionode,ionode_id
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id
    USE tdbe_common,      ONLY : ph_lw 
    !
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! error info
    INTEGER :: ierr
    !! error info
    INTEGER :: iq
    !! index q
    INTEGER :: nu
    !! index of phonon branch
    WRITE(stdout, '(5x,a)') 'Read phonon linewidth from file ph_lw.dat'
    ALLOCATE (ph_lw(nmodes,nqtotf),STAT =ierr)
    ph_lw = 0.d0
    IF (ierr /= 0) CALL errore('read_ph_lw',' Error allcoating ph_lw',1)
    IF (my_pool_id == ionode_id) THEN
      open(unit = 2, file = 'ph_lw.dat',status="old",IOSTAT = ios)
      DO iq = 1, nqtotf
        DO nu = 1, nmodes
          READ(2,*) ph_lw(nu,iq)
        ENDDO
      ENDDO
      close(2)
      !ph_lw = ph_lw /Ry2THz
    ENDIF
    CALL mp_bcast(ph_lw, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !-------------------------------------------------------------------------------- 
    END SUBROUTINE read_ph_lw
    !-------------------------------------------------------------------------------- 
    !
    !-------------------------------------------------------------------------------- 
    SUBROUTINE iphph_rta_calc()
    !-------------------------------------------------------------------------------- 
    !! This routine calculates ph-ph collision integral 
    !! under relaxation time approximation
    !--------------------------------------------------------------------------------     
    USE kinds,            ONLY : DP 
    USE input,            ONLY : eps_acoustic, twrite_tdbe
    USE global_var,       ONLY : wf,nqtotf
    USE modes,            ONLY : nmodes
    USE parallelism,      ONLY : fkbounds
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id 
    USE io_global,        ONLY : stdout
    USE tdbe_common,      ONLY : nphon_pre, iph_ph, a_tab, b_tab,  &
                                 c_tab, nstg, dt_in_ps, e_latt,    &
                                 phtemp, ph_lw, istart, ntph    
    USE ep_constants,     ONLY : kelvin2Ry
    !
    IMPLICIT NONE
    !
    INTEGER :: lower_qbnd
    !! lower q index
    INTEGER :: upper_qbnd
    !! upper q index
    INTEGER :: iq
    !! q index
    INTEGER :: nu 
    !! phonon branch index
    INTEGER :: istart2
    !! staring index for the guess of temperature 
    INTEGER :: itemp
    !!
    INTEGER :: istg, pp 
    !! Counter on Rugge Kutta solver
    REAL(KIND = DP) :: et
    !! Total energy of the nonequilibrium phonon
    REAL(KIND = DP) :: demin
    !! keep record of min of |et - e_latt(itemp)|
    REAL(KIND = DP) :: detmp 
    !! et - e_latt(itemp)
    REAL(KIND = DP) :: temp_it
    !! Effective lattice temperarure  
    REAL(KIND = DP) :: n_eq
    !! Equilibrium phonon population at temperature
    REAL(KIND = DP) :: nph
    !! Nonqquilibrium phonon population
    !! temp_it
    !
    CALL fkbounds(nqtotf, lower_qbnd, upper_qbnd)
    !
    et = 0.d0
    DO iq = 1, nqtotf
      DO nu = 1, nmodes
        if (wf(nu, iq) > eps_acoustic) then
          et = et + wf(nu,iq) * nphon_pre(nu,iq) / nqtotf
        ENDif
      ENDDO ! nu
    ENDDO ! iq

    IF (istart < 6) THEN
      istart2 = 1
    ELSE
      istart2 = istart - 1
    ENDIF
    !
    demin = 1000.d0
    !
    DO itemp = istart2, ntph
      detmp  = et - e_latt(itemp)
      IF (abs(detmp)< demin) THEN
        demin = abs(detmp)
        istart = itemp
      ELSE
        EXIT
      ENDIF
    ENDDO
    !
    temp_it = phtemp(istart)
    WRITE(stdout, '(5x,a,f10.3,a)') 'lattice temperature : ', temp_it, ' K'
    !
    iph_ph = 0.d0
    DO istg = 1, nstg
      DO iq = lower_qbnd, upper_qbnd
        DO nu = 1, nmodes 
          IF (wf(nu, iq) > eps_acoustic) then
            n_eq = 1.d0 /(EXP(wf(nu, iq) / (temp_it * kelvin2Ry)) - 1.d0)
            nph = nphon_pre(nu,iq)
            IF (istg > 1) THEN
              DO pp = 1, istg - 1 
                nph = nph + a_tab(istg, pp) * iph_ph(nu, iq, pp) * dt_in_ps
              ENDDO
            ENDIF
            iph_ph(nu, iq, istg) = (n_eq-nph)*ph_lw(nu,iq)
          ENDIF
        ENDDO ! nu
      ENDDO ! iq 
      CALL mp_sum(iph_ph(:, :, istg), inter_pool_comm)    
      CALL mp_barrier(inter_pool_comm)
    ENDDO ! istg
    !
    RETURN
    !
    !-------------------------------------------------------------------------------- 
    END SUBROUTINE iphph_rta_calc
    !-------------------------------------------------------------------------------- 
  !
  !-------------------------------------------------------------------------------- 
  END MODULE phph
  !-------------------------------------------------------------------------------- 
