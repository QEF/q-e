  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE qdabs
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines related to phonon assisted + direct absorption over
  !! full spectral region. To achieve this it uses the quasi-degenerate
  !  perturbation theory (QDPT)
  !! 14/04/2022 S. Tiwari: Implementation of QD with subroutines from indabs.f90
  !! 26/07/2022 S. Tiwari: Implementation of QD diagonalization with OMP
  !! 20/08/2022 S. Tiwari: Implementation of QD diagonalization using ELPA
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE qdabs_main(iq,mesh,run_quad)
    !-----------------------------------------------------------------------
    !!
    !! Main routine for QD absorption
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE control_flags, ONLY : iverbosity
    USE input,         ONLY : nstemp, fsthick, degaussw,                                  &
                              eps_acoustic, efermi_read, fermi_energy, nq_init,            &
                              vme, omegamin, omegamax, omegastep, len_mesh, meshnum,      &
                              wf_quasi, start_mesh, mode_res, QD_bin, QD_min, do_CHBB
    USE global_var,    ONLY : etf, ibndmin, nkf, epf17, wkf, nqtotf, wf, wqf,             &
                              sigmar_all, efnew, gtemp,nkqtotf,nkqf,                      &
                              omegap, epsilon2_abs, epsilon2_abs_lorenz,                  &
                              vmef, epsilon2_direct, epsilon2_indirect, epsilon2_qdirect, &
                              totf, nbndfst, nktotf, Energy, E_mesh, H_quad, xkf, xqf,    &
                              E_grid, Qa, tot, totcv, Eigenvec, Eigenval, n_q, sum_E,     &
                              epsilon2_qdirect2, tot_calc, tot_calc_DW, index_buf,        &
                              epsilon2_qdirect2_DW, index_buf_pool, tot_pool, totf_pool,  &
                              totcv_pool,size_tot, Eigenvec_alloc, Eigenvec_alloc_pool,   &
                              stop_qdabs, H_temp, size_m, c_ph, c_dir, c_ph_v, c_dir_v,   &
                              r_tot, H_mat, H_ind1, H_ind2, Eigenvec_alloc_write,         &
                              xkf_write, selecq_QD, lastq, startq
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci,         &
                              eps6, czero, eps8, eps4,eps5
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm,npool, my_pool_id, inter_image_comm 
    USE cell_base,     ONLY : omega,bg,at
    USE mp_world,      ONLY : mpime, world_comm
    USE bzgrid,        ONLY : kpmq_map
    USE low_lvl,       ONLY : init_random_seed
    USE stop,          ONLY : stop_epw
    USE mp_images,     ONLY : my_image_id, nimage
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    INTEGER, INTENT(in) :: mesh
    !! Mesh point for quasidegenerate perturbation
    INTEGER, INTENT(in) :: run_quad
    !! Run or build quasidegenerate perturbation
    !
    ! Local variables
    CHARACTER(LEN = 256) :: nameF
    !! Name of the file
    CHARACTER(LEN = 10) :: c
    !! Number of eta values, in string format
    CHARACTER(LEN = 256) :: format_string
    !! Format string
    CHARACTER(LEN = 20) :: tp
    !! Temperature, in string format
    CHARACTER(LEN = 20) :: mn
    !! Mesh number in string format
    LOGICAL :: adaptive_grid
    !! Future implementation of adaptive grid
    INTEGER :: iw
    !! Index for frequency
    INTEGER :: nomega
    !! Number of points on the photon energy axis
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: ierr
    !! Error status
    INTEGER :: ipol
    !! Index for pool
    INTEGER :: i,j,k
    !! Indices for counters
    INTEGER :: l
    !! index for counters
    INTEGER :: t
    !! index for counters
    INTEGER :: pool_id
    !! ID of  current core
    INTEGER :: assign_nq
    !! nq occupation number for MC
    INTEGER :: assign_n
    !! total occupation number for MC
    INTEGER ::  i_r
    !! random number index for MC
    INTEGER :: imode_r
    !! mode index for MC
    INTEGER :: n
    !! Occupation type
    INTEGER :: size_v
    !! Size of temporaty H_quad,
    INTEGER :: image
    !! Counter over images
    INTEGER :: size_image_maxval(nimage)
    !! Size of maximum Hamiltonian
    INTEGER, PARAMETER :: neta = 9
    !! Broadening parameter
    REAL(KIND = DP) :: ef0
    !! Fermi lefel
    REAL(KIND = DP) :: cfac
    !! Factor for oscillator strength (16 pi^2)
    REAL(KIND = DP) :: xt(3)
    !! buffer for k-point
    REAL(KIND = DP) :: weighta
    !! Gaussian weight at Eigenvalue
    REAL(KIND = DP) :: assign_nr
    !! Occupation number variable for MC
    REAL(KIND = DP) :: imode_rr
    !! mode for determining occupation number for MC
    REAL(KIND = DP) :: i_rr
    !! index for MC 
    REAL(KIND = DP) :: start
    !! Internals for timing, variables for MC integration
    REAL(KIND = DP) :: finish
    !! Internals for timing, variables for MC integration
    REAL(KIND = DP) :: r
    !! Random number for MC integration
    REAL(KIND = DP) :: tot_wf_E
    !! total energy for the MC configuration
    REAL(KIND = DP) :: weightf
    !! Weight for final state
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    !
    n = nq_init
    nomega = INT((omegamax - omegamin) / omegastep) + 1
    !
    !
    cfac = 16.d0 * pi**2
    adaptive_grid = .FALSE.
    ! END allocating variables for qdabs
    !
    CALL mp_barrier(inter_image_comm) 
    IF ((iq == lastq) .AND. (mesh == 1) .AND. (run_quad == 1)) THEN
      ! Monte-Carlo method of integration     
      !
      ! print *, 'inside lastq'
      IF (n == -3) THEN
        WRITE(stdout, '(/5x,a)') 'Monte-Carlo method &
                chosen for partition function integration'
        !
        ALLOCATE(n_q(nqtotf, nmodes), STAT = ierr)
        IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating n_q', 1)
        n_q = 0.d0
        tot_wf_E = 0.d0
        IF (my_pool_id == ionode_id) THEN
          assign_nq = 0
          DO WHILE (assign_nq == 0)
            tot_wf_E = 0.d0
            CALL init_random_seed()
            DO i = 1, nqtotf * nmodes
              !
              CALL RANDOM_NUMBER(assign_nr)
              CALL RANDOM_NUMBER(i_rr)
              CALL RANDOM_NUMBER(imode_rr)
              !
              i_r = INT(FLOOR(i_rr * nqtotf + 1.0))
              imode_r = INT(FLOOR(imode_rr * nmodes + 1.0))
              !
              IF (wf(imode_r,i_r) < eps4) THEN
                assign_n = INT(FLOOR(assign_nr * 100))
              ELSE
                assign_n = INT(FLOOR(assign_nr * 20.0 * dexp(-1.0 * wf(imode_r, i_r) / gtemp(1))))
              ENDIF
              tot_wf_E = tot_wf_E - 1.0 * wf(imode_r, i_r) * assign_n
              n_q(i_r, imode_r) = assign_n
              IF (-1.0 * tot_wf_E > 10.0 * gtemp(1)) EXIT
            ENDDO
            CALL RANDOM_NUMBER(r)
            !
            WRITE(stdout,'(/5x,a,4E22.14)')'configuration E,r,gtemp,tot_wf_E',&
                     DEXP((tot_wf_E) * (gtemp(1))**(-1)), r, gtemp(1), tot_wf_E
            !
            IF (DEXP((tot_wf_E) / gtemp(1)) > r) THEN
              !
              assign_nq = 1
              WRITE(stdout,'(/5x,a)')'Found Configuration'
              !
            ENDIF
          ENDDO
          nameF = 'Distribution.dat'
          OPEN(UNIT = iuindabs, FILE = nameF)
          WRITE(iuindabs, '(a)') '# phonon distribution'
          DO i = 1, nqtotf
            DO imode = 1, nmodes
              WRITE(iuindabs, '(I10,E22.14)') INT(n_q(i, imode)), wf(imode, i) * ryd2ev
            ENDDO
          ENDDO
          CLOSE(iuindabs)
        ENDIF
        sum_E = 0.d0
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(n_q, inter_pool_comm)
      ENDIF
      ! Monte-Carlo integration
      ! Building selecq_QD
      !CALL mp_barrier(inter_pool_comm)
      !CALL mp_sum(selecq_QD, inter_pool_comm)
      !CALL mp_barrier(inter_pool_comm)
      !selecq_QD(:, 1) = 0
      !
      !DO i = 1, meshnum
      !  selecq_QD(i, 1) = selecq_QD(i, 1) + 1
      !  selecq_QD(i, 2) = 1
      !  DO j = 2, nqtotf-1
      !    IF (selecq_QD(i, j + 1) > 0) THEN
      !      selecq_QD(i, 1) = selecq_QD(i, 1) + 1
      !      selecq_QD(i, selecq_QD(i, 1) + 1) = j
      !    ENDIF
      !  ENDDO
      !ENDDO
      !CALL mp_barrier(inter_pool_comm)
    ENDIF
    ! Starting to build quasi degnerate states and Hamiltonian
    ! 
    CALL mp_barrier(inter_image_comm)
    IF (mesh > 0) THEN ! This is here in case we implement adaptive grid
      ! Build H_quad which has all the interaction strengths over k and k+q
      !
      IF ((iq == startq) .AND. (mesh > 1) .AND. (run_quad == 1)) THEN
        !
        WRITE(stdout, '(/5x,a,I10)') 'Performing mesh =', mesh
        !
        IF (INT(Energy(mesh, my_pool_id + 1, my_image_id + 1, 1)) > 0) THEN
          !  
          ALLOCATE(H_quad(INT(Energy(mesh, my_pool_id + 1, my_image_id + 1, 1)), 14), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main','Error allocating H_quad', 1)
          H_quad = -1!0.d0
          H_quad(:, 13) = Energy(mesh, my_pool_id + 1, my_image_id + 1, 1)
          !
        ELSE
          ALLOCATE(H_quad(1, 14), STAT=ierr)
          IF (ierr /= 0) CALL errore('qdabs_main','Error allocating H_quad', 1)
          H_quad = -1!0.d0
          H_quad(:, 13) = -1
        ENDIF
        !
        DO image = 1, nimage
           size_image_maxval(image) = INT(MAXVAL(Energy(mesh,:,image,1))) 
        ENDDO
        ALLOCATE(H_temp(INT(MAXVAL(size_image_maxval(:))), 13), STAT=ierr)
        IF (ierr /= 0) CALL errore('qdabs_main','Error allocating H_temp', 1)
        H_temp = 0.d0
        !
        IF (iverbosity == 5) WRITE(stdout,'(/5x,a,3I10)') 'Size of H_temp',nqtotf,INT(MAXVAL(Energy(mesh,:, my_image_id + 1, 1))),&
                                                  INT(Energy(mesh, my_pool_id + 1, my_image_id + 1, 1))
      ENDIF
      CALL mp_barrier(inter_pool_comm)
      CALL mp_barrier(inter_image_comm)
      !   
      IF ((len_mesh > 2) .AND. (run_quad == 1)) THEN
        !
      !  WRITE(stdout, '(a,I10)') 'doing eig', iq
        CALL build_quasi_eig(iq, E_mesh, len_mesh, mesh, n)
        !
      ENDIF
      CALL mp_barrier(inter_pool_comm)
      CALL mp_barrier(inter_image_comm)
 
      ! H_quad accunulated for last q point
      ! Now build H_mat and diagonalize to obtain Eigenvec
      !
      IF ((iq == lastq) .AND. (meshnum > 1) .AND. (run_quad == 1)) THEN
        !
        !mesh=1 is the first mesh where the QD states are calculated
        !
        IF (mesh == 1) THEN
          !
          CALL mp_barrier(inter_image_comm)
          CALL mp_sum(Energy, inter_pool_comm)
          CALL mp_sum(Energy, inter_image_comm) 
          CALL mp_barrier(inter_pool_comm)
          CALL mp_barrier(inter_image_comm)
          !  
          E_grid = 0.d0
          !
          DO image = 1, nimage  
            DO i = 1, npool
              DO j = 1, len_mesh
                E_grid(j) = E_grid(j)+Energy(j, i, nimage, 1)
              ENDDO
            ENDDO
          ENDDO
          CALL mp_barrier(inter_pool_comm)
          CALL mp_bcast(E_grid, ionode_id, inter_pool_comm)
          !
          DO image = 2, nimage  
            DO i = 2, npool
              DO j = 1,len_mesh
                Energy(j, 1, 1, 2) = Energy(j, 1, 1, 2)+Energy(j, i, image, 2)
                Energy(j, 1, 1, 3) = Energy(j, 1, 1, 3)+Energy(j, i, image, 3)
              ENDDO
            ENDDO
          ENDDO
          ! Write Energies which have the JDOS and QD JDOS  
          IF (my_pool_id == ionode_id) THEN
            !
            nameF = 'Energy.dat'
            OPEN(UNIT = iuindabs, FILE = nameF)
            WRITE(iuindabs, '(a)') '# Energy states for quasidegenerate case'
            WRITE(iuindabs, '(a)') '# omega (eV), No of interactions, JDOS, QD-JDOS'
            DO iw = 1, len_mesh
              WRITE(iuindabs, '(4E22.14)') E_mesh(iw) * ryd2ev,  &
                  E_grid(iw), Energy(iw, 1, 1, 2), Energy(iw, 1, 1, 3)
            ENDDO
            CLOSE(iuindabs)
          ENDIF
          CALL mp_barrier(inter_pool_comm)
          WRITE(stdout, '(5x,a)') 'Finished building matrix'
          !
        ENDIF
        ! END writing Energy
        !
        WRITE(stdout, '(5x,a)') 'Finished sorting, Entering eigenvalue calculation'
        IF (iverbosity == 5) WRITE(stdout, '(5x,a,I10)') 'Number of states in mesh = ', INT(E_grid(mesh))
        !
        CALL mp_barrier(inter_pool_comm)
        !
        ! Prepare for Eigenvalue and Eigenvector solution
        !
        ! Variables with "_pool" are local to processors
        IF (E_grid(mesh) > 0) THEN
          ! Solving the Eigenvector problem
          ! Here we allocate H_mat which keeps the interactions only for
          ! size_v*size_m
          ! Maximum value size_v=4 however that's an extreme case and in
          ! practice size_v=3 works out
          ! One should always check if r_tot<size_v*size_m
          !
          size_m=INT(E_grid(mesh))
          ! This is risky because if we assume all states are independent, the 
          ! maximum size should be 4*size(H_mat) but ideally that does not happen
          ! Weird eigenvalues are a result of this problem but there is no
          ! alternative     
          IF (size_m < 500) THEN
            size_v=4
          ELSE
            size_v=3
          ENDIF
          WRITE(stdout, '(/5x,a,I10)') 'QD array rank = ', size_m
          !
          ALLOCATE(Eigenval(size_m), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating Eigenval(size_m)', 1)
          ALLOCATE(index_buf(size_m,6), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating index_buf', 1)
          index_buf = -1

          ALLOCATE(H_mat(INT(size_v*size_m)), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating H(size_m,3)', 1)
          ALLOCATE(H_ind1(INT(size_v*size_m)), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating H(size_m,3)', 1)
          ALLOCATE(H_ind2(INT(size_v*size_m)), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating H(size_m,3)', 1)
          !
          CALL mp_barrier(inter_pool_comm)
          !
          CALL build_quasi_hamiltonian(E_mesh, len_mesh, mesh, n)
          !
          totf= INT(tot-totcv)
          !  
          CALL mp_barrier(inter_pool_comm)
          !
          DEALLOCATE(H_quad, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs', 'Error deallocating H_quad', 1)
          DEALLOCATE(H_temp, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs', 'Error deallocating H_temp', 1)
          !
          IF (iverbosity == 5) THEN
            WRITE(stdout,'(/5x,a,I10,I10,3I10)') 'tot, totf,totcv,size_tot', tot, totf, totcv, size_tot
          ENDIF
          !
          IF (size_tot > 0) THEN
            ALLOCATE(Eigenvec_alloc_pool(size_tot, tot), STAT = ierr)
            IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating Eigenvec_alloc_pool2(size_tot,tot)', 1)
            ALLOCATE(index_buf_pool(size_tot, 6), STAT = ierr)
            IF (ierr /= 0)CALL errore('qdabs_main', 'Error allocating index_buf_pool(tot,5)', 1)
            Eigenvec_alloc_pool = 0.d0
            index_buf_pool = -1
          ELSE
            ALLOCATE(Eigenvec_alloc_pool(2, 2), STAT = ierr)
            IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating Eigenvec_alloc_pool2(tot,tot)', 1)
            ALLOCATE(index_buf_pool(2, 2), STAT = ierr)
            IF (ierr /= 0)CALL errore('qdabs_main', 'Error allocating index_buf_pool(tot,5)', 1)
            Eigenvec_alloc_pool = 0.d0
            index_buf_pool = -1
          ENDIF
          CALL mp_barrier(inter_pool_comm)
          CALL mp_barrier(inter_image_comm)
          !   
          IF (iverbosity == 5) WRITE(stdout,'(/5x,a)') 'Allocated Eigenvec_pool'
          !  
          CALL diag_quasi_hamiltonian(E_mesh, len_mesh, mesh, n)
          !  
          CALL mp_barrier(inter_pool_comm)
          CALL mp_barrier(inter_image_comm) 
          !  
          IF (wf_quasi /= -1 ) THEN
            IF (wf_quasi == 1) THEN
              ALLOCATE(Eigenvec_alloc_write(npool, 30, tot), STAT=ierr)
              IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating Eigenvec_alloc_write(tot)', 1)
              ALLOCATE(xkf_write(npool, 30, tot, 6), STAT=ierr)
              IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating xkf_write(tot)', 1)
              !  
            ENDIF
            CALL write_wf_quasi(mesh, wf_quasi)
            !
          ENDIF
          !   
          CALL mp_barrier(inter_pool_comm)
          !
          DEALLOCATE(H_mat, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating H(size_m,3)', 1)
          DEALLOCATE(H_ind1, STAT=ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating H(size_m,3)', 1)
          DEALLOCATE(H_ind2, STAT=ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating H(size_m,3)', 1)
          !
          ALLOCATE(epsilon2_qdirect2(3, tot, 1 + mode_res), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating epsilon2_qdirect2', 1)
          ALLOCATE(epsilon2_qdirect2_DW(3, tot), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating epsilon2_qdirect2_DW', 1)
          ALLOCATE(c_ph_v(3,tot), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating epsilon2_qdirect2', 1)
          ALLOCATE(c_dir_v(3,tot), STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error allocating epsilon2_qdirect2', 1)
          !
          epsilon2_qdirect2 = 0.d0
          epsilon2_qdirect2_DW = 0.d0
          c_ph_v = 0.d0
          c_dir_v = 0.d0
          !  
          ! Obtain maximum number of elements per pool
          !  
          WRITE(stdout,'(/5x,a)') 'Finished distributing Eigenvectors'
          !
          CALL mp_bcast(Eigenval,ionode_id, inter_pool_comm)
          CALL mp_barrier(inter_pool_comm)
          CALL mp_barrier(inter_image_comm)
          ! Distribution of Eigenvectors finished 
          IF (iverbosity == 5) THEN
            WRITE(stdout,'(a,I10,2I10)') 'Number of total states, phonon states, and direct states', &
                    INT(tot), INT(totf), INT(totcv)
            !
            DO i = 1, totcv
              WRITE(stdout,'(a,I10,I10,3I10)') 'indextotcv',INT(index_buf(i, 4)), INT(index_buf(i, :))
              !
            ENDDO
          ENDIF
        ENDIF
      ENDIF ! Finished building eigenvectors (run_quad .Eq. 1 is finished !
      CALL mp_barrier(inter_pool_comm)
      !
      ! Finished building eigenvectors (first phase of calculations) !
      ! CALL cpu_time(start)
      !
      ! First do simple calculation for energies with no QD states
      !
      IF ((mesh == 1) .AND. (run_quad == 2) .AND. (start_mesh == 0)) THEN
        !
        IF (do_CHBB) THEN
          IF (iq == 1) THEN
            WRITE(stdout,'(/5x,a)') 'You are additionally performing CHBB   &
                                     calculations (expect slowdown) '
          ENDIF
          CALL calc_CHBB(nomega, iq)
        ENDIF
        DO i = 1,len_mesh-1
          IF (E_grid(i) < 1) THEN
            CALL calc_qdirect_DW(i, nomega, iq, n, sum_E, 1, tot_calc_DW)
          ENDIF
        ENDDO
        !
      ENDIF
      !
      ! Now calculate for QD states 
      !
      IF ((run_quad == 2) .AND. (E_grid(mesh) > 0) .AND. (stop_qdabs == 0)) THEN
        !
        ! epsilon2 for each QD state
        !
        CALL calc_qdirect(mesh, nomega, iq, n, sum_E, tot_calc, .False.)
        ! 
        IF ((mesh > 0) .AND. (mesh < len_mesh)) THEN
          ! states which are in the QD bin but not QD
          !
          CALL calc_qdirect_DW(mesh, nomega, iq, n, sum_E, 2, tot_calc_DW)
          ! 
        ENDIF
        !
        ! Sum QD epsilon2 for all q points
        IF (iq == lastq) THEN
          !
          CALL mp_barrier(inter_pool_comm)
          CALL mp_barrier(inter_image_comm)
          !  
          CALL mp_sum(epsilon2_qdirect2_DW, inter_image_comm)
          CALL mp_sum(epsilon2_qdirect2_DW, inter_pool_comm)
          !!!!!!!!!!! For the second order term which is first order in QD
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CALL calc_qdirect(mesh,nomega,iq,n,sum_E,tot_calc,.TRUE.)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !CALL mp_barrier(inter_image_comm)
          !CALL mp_sum(epsilon2_qdirect2, inter_image_comm)
          !CALL mp_sum(c_ph_v, inter_image_comm)
          !CALL mp_sum(c_dir_v, inter_image_comm)
          !CALL mp_sum(tot_calc, inter_image_comm)
          !CALL mp_sum(tot_calc_DW, inter_image_comm)
          !CALL mp_barrier(inter_image_comm)
          CALL mp_barrier(inter_pool_comm)
          CALL mp_sum(epsilon2_qdirect2, inter_pool_comm)
          CALL mp_sum(c_ph_v, inter_pool_comm)
          CALL mp_sum(c_dir_v, inter_pool_comm)
          CALL mp_sum(tot_calc, inter_pool_comm)
          CALL mp_sum(tot_calc_DW, inter_pool_comm)
          CALL mp_barrier(inter_pool_comm)
          !
          IF ((iverbosity == 5) .AND. (ABS(tot-tot_calc_DW) > 0)) THEN
            ! This checks if all the states in the QD bin are found when
            ! summation over k,q is performed
            ! For mp_mesh_k, tot and tot_calc may not be equal because of extra
            ! states included in QD bin |vk-q>|ck>.
            WRITE(stdout,'(5x,a,I10,3I10)') 'found degenerte states, tot_calc,tot,tot_calc_DW', &
                                             tot_calc, tot, tot_calc_DW, tot_pool
            WRITE(stdout,'(5x,a)') 'Normal behavior for mp_mesh_k .EQ. .True.'
            !
          ENDIF
          !
          ! Sum over all QD states !!!!!!!!!!!!!!!!!
          !
          DO i = 1, tot
            !
            IF (iverbosity == 5) THEN
              !   
              WRITE(stdout,'(/5x,a,E22.14)') 'Eigenval', REAL(Eigenval(i))*ryd2ev
            ENDIF
            ! Sum over all photon energies !!!!!!!!!!!!!!!!!1
            DO iw = 1, nomega
              !
              weightf = w0gauss((((E_mesh(mesh)+E_mesh(mesh+1))/2.0) - omegap(iw)) / degaussw, 0) / degaussw
              weighta = w0gauss((REAL(Eigenval(i)) - omegap(iw)) / degaussw, 0) / degaussw
              !
              DO ipol = 1, 3
                epsilon2_abs(ipol, iw, 1, 1) = epsilon2_abs(ipol, iw, 1, 1) + &
                                       (ABS(epsilon2_qdirect2(ipol, i, 1))**2) &
                               * weighta * cfac * (1.0/omegap(iw)**2) * (1.0/omega)
                !
                epsilon2_abs(ipol, iw, 2, 1) = epsilon2_abs(ipol, iw, 2, 1) + &
                                       (ABS(epsilon2_qdirect2(ipol, i, 1))**2) &
                               * weightf * cfac * (1.0/omegap(iw)**2) * (1.0/omega)
                !
                IF (mode_res>0) THEN
                  DO imode = 1, mode_res
                    epsilon2_abs(ipol, iw, neta + imode, 1) = epsilon2_abs(ipol, iw, neta + imode, 1) + &
                                                         (ABS(epsilon2_qdirect2(ipol, i, imode + 1))**2) &
                                                       * weighta * cfac * (1.0/omegap(iw)**2) * (1.0/omega)
                  ENDDO
                ENDIF

                IF (my_pool_id == ionode_id) THEN
                  c_ph(ipol, iw) = c_ph(ipol, iw) + (ABS(c_ph_v(ipol, i))**2) &
                                 * weighta * cfac * (1.0/omegap(iw)**2) * (1.0/omega)
                  c_dir(ipol, iw) = c_dir(ipol, iw) + (ABS(c_dir_v(ipol, i))**2) &
                                 * weighta * cfac * (1.0/omegap(iw)**2) * (1.0/omega)
                !
                ENDIF
              ENDDO ! ipol
            ENDDO ! Photon
          ENDDO ! tot, Eigenval
          CALL mp_barrier(inter_pool_comm)
          CALL mp_barrier(inter_image_comm)
 
          ! QD  
          IF (iverbosity == 5) THEN
            ! To debug distribution of k points !
            l = 0
            DO i = 1,tot_pool
              !
              IF (INT(index_buf_pool(i, 6) - 1) /= my_pool_id) THEN
                !
                PRINT '(a,2I10)', 'pool id', my_pool_id, INT(index_buf_pool(i, 6) - 1)
                !
              ENDIF
              l = l+1
            ENDDO
          ENDIF
          CALL mp_barrier(inter_pool_comm)
          CALL mp_barrier(inter_image_comm)
          CALL mp_sum(tot_pool, inter_image_comm)
          CALL mp_sum(tot_pool, inter_pool_comm)
          IF (iverbosity == 5) WRITE(stdout,'(/5x,a,I10)') 'tot_pool_sum =', tot_pool
          !
          epsilon2_abs(:, :, 3, 1) = epsilon2_abs(:, :, 1, 1)
          ! Set everything to startig values for next meshpoint
          ! calculation
          !
          ! Deallocate everything that belongs to pools which changes
          ! in next mesh calculation 
          DEALLOCATE(Eigenvec_alloc_pool, STAT = ierr)
          IF (ierr /= 0)CALL errore('qdabs_main', 'Error deallocating Eigenvec_alloc_pool2(tot,tot)', 1)
          DEALLOCATE(index_buf_pool, STAT = ierr)
          IF (ierr /= 0)CALL errore('qdabs_main', 'Error deallocating index_buf_pool(tot,5)', 1)
          DEALLOCATE(Eigenval, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating Eigenval(size_m)', 1)
          DEALLOCATE(index_buf, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating index_buf', 1)
          DEALLOCATE(epsilon2_qdirect2, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating epsilon2_qdirect2', 1)
          DEALLOCATE(epsilon2_qdirect2_DW, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating epsilon2_qdirect2_DW', 1)
          DEALLOCATE(c_ph_v, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating epsilon2_qdirect2', 1)
          DEALLOCATE(c_dir_v, STAT = ierr)
          IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating epsilon2_qdirect2', 1)
          IF (wf_quasi == 1 ) THEN
            ! 
            DEALLOCATE(Eigenvec_alloc_write, STAT = ierr)
            IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating Eigenvec_alloc_write(tot)', 1)
            DEALLOCATE(xkf_write, STAT = ierr)
            IF (ierr /= 0) CALL errore('qdabs_main', 'Error deallocating Eigenvec_alloc_write(tot)', 1)
            ! 
          ENDIF
          !
          CALL mp_barrier(inter_image_comm)
          CALL mp_barrier(inter_pool_comm)
          !
          tot_calc = 0
          tot_calc_DW = 0
          tot_pool = 0
          r_tot = 0
          size_tot = 0
          tot = 0
          !
        ENDIF
      ENDIF
      CALL mp_barrier(inter_image_comm)
      CALL mp_barrier(inter_pool_comm)
    ENDIF
    ! Last calculation of mesh write everything
    ! stop_qdabs for future imeplementation to stop calculation at a meshgrid
    !
    IF (((iq == lastq) .AND. (mesh == meshnum) .AND. &
       (run_quad == 2)) .OR. (stop_qdabs == 1)) THEN
      !  
      ! Sum QD and non-QD part
      !
      CALL write_absorption()
      !
      stop_qdabs = 1
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    ! END Write
    !
    IF (stop_qdabs == 1) THEN
      !
      WRITE(stdout,'(/5x,a,I10)') 'Stopping EPW because stop_qdabs encountered on mesh', INT(mesh)
      CALL stop_epw()
      !
      STOP
    ENDIF
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE qdabs_main
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE prepare_qdabs(nrr_q, irvec_q, ndegen_q, rws, nrws)
    !---------------------------------------------------------------------------
    !! This subroutine prepares the global variables whose size remains constant
    !! throughout the calculation
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE control_flags, ONLY : iverbosity
    USE input,         ONLY : nstemp, fsthick, degaussw,                                  &
                              eps_acoustic, efermi_read, fermi_energy, nq_init, vme,       &
                              omegamin, omegamax, omegastep, len_mesh, meshnum, wf_quasi, &
                              start_mesh, mode_res, QD_bin, QD_min, lifc
    USE global_var,    ONLY : nqtotf, wf, wqf,nkqtotf,nkqf, omegap, epsilon2_abs,         &
                              epsilon2_abs_lorenz, epsilon2_direct, epsilon2_indirect,    &
                              epsilon2_qdirect, nktotf, Energy, E_mesh, E_grid,           &
                              epsilon2_qdirect2, epsilon2_qdirect2_DW, stop_qdabs, c_ph,  &
                              c_dir, c_ph_v, c_dir_v, r_tot, tot_calc, tot_calc_DW, gtemp,&
                              xqf, selecq_QD
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci,         &
                              eps6, czero, eps8, eps4,eps5
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    USE cell_base,     ONLY : omega,bg,at
    USE mp_world,      ONLY : mpime, world_comm
    USE bzgrid,        ONLY : kpmq_map
    USE wannier2bloch, ONLY : dynwan2bloch, dynifc2blochf
    USE mp_images,     ONLY : nimage
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)            :: nrr_q
    !! Number of WS points for phonons 
    INTEGER, INTENT(in)            :: irvec_q(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, INTENT(in)            :: ndegen_q(:, :, :)
    !! Wigner-Seitz weights for the phonon grid that 
    !! depend on atomic positions $R + \tau(nb) - \tau(na)$
    REAL(KIND = DP), INTENT(in)    :: rws(0:3, 200)
    !! Number of real-space Wigner-Seitz vectors
    INTEGER, INTENT(in)            :: nrws
    !! Number of real-space Wigner-Seitz
    !
    ! Local variables
    CHARACTER(LEN = 256) :: nameF
    !! Name of the file
    CHARACTER(LEN = 10) :: c
    !! Number of eta values, in string format
    CHARACTER(LEN = 256) :: format_string
    !! Format string
    CHARACTER(LEN = 20) :: tp,mn
    !! Temperature, in string format
    ! Local variables
    INTEGER :: iw
    !! Index for frequency/energy
    INTEGER :: nomega
    !! Number of points on the photon energy axis
    INTEGER :: m
    !! Counter on denominator imaginary broadening values
    INTEGER :: itemp,imode
    !! Counter on temperatures
    INTEGER :: ierr,ipol
    !! Error status
    INTEGER :: n
    !! Occupation type
    INTEGER :: i
    !! Counter on mesh
    INTEGER,PARAMETER :: neta = 9
    !! Broadening parameter
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: qwidth
    !! QDPT bin width
    REAL(KIND = DP) :: w2(nmodes)
    !! phonon energy at iq=1 Gamma point
    REAL(KIND = DP) :: xxq(3)
    !! q-point
    REAL(KIND = DP) :: xxq_r(3)
    !! q-point
    COMPLEX(KIND = DP) :: uf(nmodes, nmodes)
    !! rotation matrix for the q point
    uf = 0.d0
    n = nq_init
    qwidth = 1.0 / len_mesh
    nomega = INT((omegamax - omegamin) / omegastep) + 1
    w2 = zero
    !
    !
    ! Allocating variables for qdabs
    !
    WRITE(stdout, '(/5x,a)') REPEAT('=',67)
    WRITE(stdout, '(5x,"Phonon-assisted absorption with quasidegenerate perturbation")')
    WRITE(stdout, '(5x,a/)') REPEAT('=',67)
    !
    IF (fsthick < 1.d3) WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a)') 'The following temperatures are calculated:'
    DO itemp = 1, nstemp
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Temperature T = ', gtemp(itemp) * ryd2ev, ' eV'
    ENDDO
    WRITE(stdout,'(/5x,a,f10.6,a)' ) 'size for cell = ', (omega), ' per A'
    WRITE(stdout, '(/5x,a,I10)' ) 'nqtotf = ', nqtotf
    WRITE(stdout, '(/5x,a,I10)' ) 'nktotf = ', nktotf
    ! 
    ALLOCATE(omegap(nomega), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating omegap', 1)
    ALLOCATE(epsilon2_abs(3, nomega, neta + mode_res, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating epsilon2_abs', 1)
    ALLOCATE(epsilon2_abs_lorenz(3, nomega, neta, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating epsilon2_abs_lorenz', 1)
    ALLOCATE(epsilon2_direct(3, nomega, neta, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating epsilon2_abs', 1)
    ALLOCATE(epsilon2_indirect(3, nomega, neta, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating epsilon2_abs', 1)
    ALLOCATE(epsilon2_qdirect(3, nomega, neta+mode_res, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating epsilon2_abs', 1)
    ALLOCATE(c_ph(3,nomega), STAT=ierr)
    IF (ierr /= 0) CALL  errore('prepare_qdabs', 'Error allocating c_ph', 1)
    ALLOCATE(c_dir(3,nomega), STAT=ierr)
    IF (ierr /= 0) CALL  errore('prepare_qdabs', 'Error allocating c_dir', 1)
    ALLOCATE(selecq_QD(meshnum, nqtotf + 2), STAT = ierr)
    IF (ierr /= 0) CALL  errore('prepare_qdabs', 'Error allocating selecq_QD', 1) 
    !
#if defined(__ELPA)
    WRITE(stdout,'(/5x,a)') '************************************************************************ '
    WRITE(stdout,'(/5x,a)') '            QDPT calculation accelerated with ELPA Library               '
    WRITE(stdout,'(/5x,a)') 'Refer: https://elpa.mpcdf.mpg.de/ELPA_PUBLICATIONS.html for citation     '
    WRITE(stdout,'(/5x,a)') '************************************************************************ '
#else
    WRITE(stdout,'(/5x,a)') '************************************************************************ '
    WRITE(stdout,'(/5x,a)') '                      ELPA Library not found                             '
    WRITE(stdout,'(/5x,a)') '   Please install EPW with ELPA for faster and denser calculations       '
    WRITE(stdout,'(/5x,a)') '          Visit: https://elpa.mpcdf.mpg.de for installing ELPA           '
    WRITE(stdout,'(/5x,a)') '************************************************************************ '
#endif
    !
    c_ph = 0.d0
    c_dir = 0.d0
    stop_qdabs = 0
    tot_calc = 0
    tot_calc_DW = 0
    epsilon2_abs = 0.d0
    epsilon2_abs_lorenz = 0.d0
    epsilon2_direct = 0.d0
    epsilon2_indirect = 0.d0
    epsilon2_qdirect = 0.d0
    selecq_QD(:,:) = 0
    !
    xxq = xqf(:, 1)
    xxq_r = xxq
    ! ------------------------------------------------------
    ! S.Tiwari, we need phonon energy at Gamma point to determine
    ! optimal QD_bin
    ! dynamical matrix : Wannier -> Bloch
    ! ------------------------------------------------------
    !
    IF (.NOT. lifc) THEN
      CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, .false.)
    ELSE
      CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, .false.)
    ENDIF
    ! 
    DO iw = 1, nomega
      !
      omegap(iw) = omegamin + (iw - 1) * omegastep
      !
    ENDDO
    !
    IF (QD_bin <= 0) THEN
      !
      QD_bin = DSQRT(DSQRT(ABS(w2(nmodes)))*ryd2ev*1000)*30 / 1000
      !
    ENDIF
    IF (len_mesh == 1) THEN
      !
      len_mesh = INT(((omegamax * ryd2ev)-QD_min)/QD_bin)+2
      meshnum = len_mesh - 1
      !
      WRITE(stdout, '(/5x,a)') '**************************************************** '
      WRITE(stdout, '(/5x,a)') '*****  Automatic generation of QDPT meshgrids  ***** '
      WRITE(stdout, '(/5x,a,I10)') ' Number of QDPT meshgrids to be performed: ', meshnum
      WRITE(stdout, '(/5x,a)') '**************************************************** '
      !
    ELSE
      !
      WRITE(stdout, '(/5x,a,I10)') 'Number of QDPT meshgrids to be performed: ', meshnum
      !
    ENDIF
    ! 
    ALLOCATE(Energy(len_mesh, npool, nimage, 3), STAT=ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating Energy', 1)
    ALLOCATE(E_mesh(len_mesh), STAT=ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating Emesh', 1)
    ALLOCATE(E_grid(len_mesh), STAT=ierr)
    IF (ierr /= 0) CALL errore('prepare_qdabs', 'Error allocating Egrid', 1)
    E_mesh = 0.d0
    Energy = 0.d0
    E_grid=0.d0
    !
    DO i = 1,len_mesh
      IF (i == 1) THEN
        E_mesh(i) = 0.d0
      ELSEIF (i == 2) THEN
        E_mesh(i) = QD_min / ryd2ev
      ELSE
        E_mesh(i) = E_mesh(i-1) + QD_bin / ryd2ev
      ENDIF
    ENDDO
    WRITE(stdout, '(5x,a,f10.2,a)') 'Size of QD bins', DBLE(E_mesh(len_mesh)-E_mesh(len_mesh-1)) &
                 * ryd2ev * 1000, 'meV'
    !
    CALL mp_bcast(tot_calc_DW,ionode_id, inter_pool_comm)
    CALL mp_bcast(tot_calc,ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (my_pool_id == ionode_id) THEN
      IF (wf_quasi > 0) CALL system('mkdir quasi_bands')
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE prepare_qdabs
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE write_absorption()
    !--------------------------------------------------------------------------
    !!
    !! Write imaginary dielectric constant and related properties
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE control_flags, ONLY : iverbosity
    USE input,         ONLY : nstemp, omegamin, omegamax, omegastep,                      &
                              start_mesh, mode_res, QD_bin, QD_min, meshnum, do_CHBB
    USE global_var,    ONLY : gtemp, nkqtotf, nkqf, c_dir, c_ph,                          &
                              omegap, epsilon2_abs, epsilon2_abs_lorenz,                  &
                              epsilon2_direct, epsilon2_indirect, epsilon2_qdirect,       &
                              epsilon2_qdirect2, epsilon2_qdirect2_DW, sum_E
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id, my_image_id, inter_image_comm
    USE stop,          ONLY : stop_epw
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 256) :: nameF
    !! Name of the file
    CHARACTER(LEN = 10) :: c
    !! Number of eta values, in string format
    CHARACTER(LEN = 256) :: format_string
    !! Format string
    CHARACTER(LEN = 20) :: tp,mn
    !! Temperature, in string format
    ! Local variables
    INTEGER :: iw
    !! Index for frequency
    INTEGER :: nomega
    !! Number of points on the photon energy axis
    INTEGER :: m
    !! Counter on denominator imaginary broadening values
    INTEGER :: itemp,imode
    !! Counter on temperatures and modes
    INTEGER :: ierr,ipol
    !! Error status and current pool index
    INTEGER,PARAMETER :: neta = 9
    !! Number of broadening parameter for CHBB theory  
    INTEGER :: neta1
    !! Number of modes+broadening
    !
    CALL mp_barrier(inter_pool_comm)
    CALL mp_barrier(inter_image_comm)
    CALL mp_sum(epsilon2_indirect, inter_pool_comm)
    CALL mp_sum(epsilon2_direct, inter_pool_comm)
    CALL mp_sum(epsilon2_qdirect, inter_pool_comm)
    CALL mp_sum(epsilon2_abs_lorenz, inter_pool_comm)
    CALL mp_sum(sum_E, inter_pool_comm)
    CALL mp_sum(c_ph, inter_pool_comm)
    CALL mp_sum(c_dir, inter_pool_comm)
    CALL mp_sum(epsilon2_indirect, inter_image_comm)
    CALL mp_sum(epsilon2_direct, inter_image_comm)
    CALL mp_sum(epsilon2_abs_lorenz, inter_image_comm)
    CALL mp_sum(sum_E, inter_image_comm)
    CALL mp_sum(c_ph, inter_image_comm)
    CALL mp_sum(c_dir, inter_image_comm)
    CALL mp_barrier(inter_image_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ! Sum QD and non-QD part
    epsilon2_abs(:,:,:,:) = epsilon2_abs(:,:,:,:) + epsilon2_qdirect(:,:,:,:)
    CALL mp_sum(epsilon2_qdirect, inter_image_comm)
    CALL mp_sum(epsilon2_abs, inter_image_comm)
    CALL mp_barrier(inter_image_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    nomega = INT((omegamax - omegamin) / omegastep) + 1
    !
    c = 'X'
    WRITE(c,"(i0)") neta
    format_string = "(5x,f15.6," // TRIM(c) // "E22.14)"
    !
    WRITE(stdout, '(5x,a)')
    WRITE(stdout, '(5x,a)') 'Optical-absorption versus energy'
    WRITE(stdout, '(5x,a)')
    ! For test-farm checking purposes, only show m=1
    WRITE(stdout, '(5x,a)') 'Photon energy (eV), Imaginary dielectric function along x,y,z'
    !
    DO iw = 1, nomega
      WRITE(stdout, '(5x,f15.6,3E22.14)') omegap(iw) * ryd2ev,   &
                       (epsilon2_abs(ipol, iw, 1, 1), ipol = 1, 3)
    ENDDO
    WRITE(stdout, '(5x,a)')
    WRITE(stdout, '(5x,a)') 'directionally averaged imaginary dielectric constant for &
                             & temperature X are  reported in the files               &
                             & epsilon2_indabs_X_finalmesh.dat'
    WRITE(stdout, '(5x,a)')
    !
    ! Output to file
    neta1=neta+mode_res
    IF (my_image_id == 0) THEN
      IF (my_pool_id == ionode_id) THEN
        DO itemp = 1,nstemp
          WRITE(c,"(i0)") neta1 + 1
          WRITE(tp,"(f8.1)") gtemp(itemp) * ryd2ev / kelvin2eV
          WRITE(mn,"(I10)") meshnum
          format_string = "("//TRIM(c) // "E22.14)"
          nameF = 'epsilon2_indabs_' // trim(adjustl(tp)) // 'K_' // trim(adjustl(mn)) // '.dat'
          OPEN(UNIT = iuindabs, FILE = nameF)
          WRITE(iuindabs, '(a)') '# Phonon-assisted absorption versus energy using QDPT'
          WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric &
                                   function along x,y,z'
          WRITE(iuindabs, '(a)') '# omega(eV), im(epsilon), im(epsilon(mid-eigs)), im(epsilon_QD), & 
                                  im(epsilon_phonon), im(epsilon_direct), 0, 0, 0, 0, im(epsilon2(modes) '
          DO iw = 1, nomega
            WRITE(iuindabs, format_string) omegap(iw) * ryd2ev,(SUM(epsilon2_abs(:, iw, m, itemp)) / 3.0d0, m = 1, neta1)
          ENDDO
          CLOSE(iuindabs)
          !
          nameF = 'c_dir' // trim(adjustl(tp)) // 'K_' // trim(adjustl(mn)) // '.dat'
          !  
          OPEN(UNIT = iuindabs, FILE = nameF)
          WRITE(iuindabs, '(a)') '# Direct contribution to epsilon2 vs energy'
          WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary  & 
                                  dielectric function along x,y,z'
          DO iw = 1, nomega
            WRITE(iuindabs,'(/5x,2E22.14)' ) omegap(iw) * ryd2ev,(SUM(c_dir(:,iw))/3.0d0)
          ENDDO
          CLOSE(iuindabs)
          !  
          nameF = 'c_ph' // trim(adjustl(tp)) // 'K_' // trim(adjustl(mn)) // '.dat'
          !  
          OPEN(UNIT = iuindabs, FILE = nameF)
          WRITE(iuindabs, '(a)') '# Phonon-assisted contribution to epsilon2 versus energy'
          WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary &
                                  dielectric function along x,y,z'
          DO iw = 1, nomega
            WRITE(iuindabs, '(/5x,2E22.14)') omegap(iw) * ryd2ev, (SUM(c_ph(:,iw))/3.0d0)
          ENDDO
          CLOSE(iuindabs)
          !
          IF (do_CHBB) THEN
            nameF = 'epsilon2_direct_' // trim(adjustl(tp)) //  'K_' // trim(adjustl(mn)) // '.dat'
            !
            OPEN(UNIT = iuindabs, FILE = nameF)
            WRITE(iuindabs, '(a)') '# Direct absorption versus energy'
            WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary & 
                                    dielectric function along x,y,z'
            DO iw = 1, nomega
              WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_direct(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
            ENDDO
            CLOSE(iuindabs)
            !
            nameF = 'epsilon2_indirect_' // trim(adjustl(tp)) //  'K_' // trim(adjustl(mn)) // '.dat'
            !
            OPEN(UNIT = iuindabs, FILE = nameF)
            WRITE(iuindabs, '(a)') '# Phonon-assisted absorption versus energy'
            WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary &
                                    dielectric function along x,y,z'
            DO iw = 1, nomega
              WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_indirect(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
            ENDDO
            CLOSE(iuindabs)
          ENDIF
          !
          IF (iverbosity == 5) THEN
            nameF = 'epsilon2_qdirect_' // trim(adjustl(tp)) // 'K_' // trim(adjustl(mn)) // '.dat'
            !
            OPEN(UNIT = iuindabs, FILE = nameF)
            WRITE(iuindabs, '(a)') '# Phonon-assisted and direct absorption versus energy'
            WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary &
                                    dielectric function along x,y,z'
            DO iw = 1, nomega
              WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_qdirect(:, iw, m, itemp)) / 3.0d0, m = 1, neta1)
            ENDDO
            CLOSE(iuindabs)
            !
            nameF = 'epsilon2_indabs_lorenz' // trim(adjustl(tp)) //  'K_' // trim(adjustl(mn)) // '.dat'
            !
            OPEN(UNIT = iuindabs, FILE = nameF)
            WRITE(iuindabs, '(a)') '# Phonon-assisted absorption versus energy'
            WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary &
                                     dielectric function along x,y,z'
            DO iw = 1, nomega
              WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_abs_lorenz(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
            ENDDO
            CLOSE(iuindabs)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    CALL mp_barrier(inter_image_comm)
    CALL stop_epw()     
    !-----------------------------------------------------------------------
    END SUBROUTINE write_absorption
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION check_element(initial_E, final_Ev, final_Ec, bound_up, bound_dn)
    !-----------------------------------------------------------------------
    !! Check if the element falls within the boundary, for future use of
    !! adaptive grid
    USE kinds,         ONLY : DP
    !
    REAL(KIND = DP), INTENT (in) :: initial_E
    !! Initial energy value for MB state
    REAL(KIND = DP), INTENT (in) :: final_Ec
    !! Final energy for c-c' transitions
    REAL(KIND = DP), INTENT (in) :: final_Ev
    !! Final energy for v-v' transitions
    REAL(KIND = DP), INTENT (in) :: bound_up
    !! Upper bound
    REAL(KIND = DP), INTENT (in) :: bound_dn
    !! Lower bound
    INTEGER :: check_element
    !! Check if element in the meshgrid
    !
    IF (((initial_E > bound_dn) .AND. (initial_E <= bound_up)) &
     .AND. (((final_Ev <= bound_up) .AND. (final_Ev >= bound_dn)) &
     .OR. ((final_Ec <= bound_up) .AND. (final_Ec >= bound_dn)))) THEN
      !
      check_element = 1
      !  
    ELSE
      !  
      check_element = 0
    ENDIF
    !
    !--------------------------------------------------------------------------
    END FUNCTION
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE build_quasi_eig(iq, E_mesh, len_mesh, meshin, n)
    !--------------------------------------------------------------------------
    !! This subroutine builds the QD Hamiltonian H_quad for meshin,
    !! for meshin==1 it  calculates the nuumber of interactions/elements
    !! present in H_quad for all mesh grids
    !! H_quad is different for each pool
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE control_flags, ONLY : iverbosity
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nstemp, fermi_energy, eps_acoustic,               & 
                              degaussw, fsthick,filkf
    USE global_var,    ONLY : etf, ibndmin, nkf,nqtotf, wf, wqf,               &
                              nbndfst, nktotf, nkqtotf, epf17, nkqf, xkf, xqf, &
                              wkf, gtemp, totf_pool, totcv_pool, tot_pool,     &
                              index_buf_pool, H_quad, stop_qdabs, Energy,      &   
                              maxdim, selecq_QD
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_global,     ONLY : my_pool_id, npool, inter_pool_comm,              &
                              inter_image_comm
    USE cell_base,     ONLY : omega
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi,  &
                              ci, eps40, czero, eps8
    USE bzgrid,        ONLY : kpmq_map
    USE mp_images,     ONLY : nimage, my_image_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Current q index
    REAL(KIND = DP), intent(in) :: E_mesh(len_mesh)
    !! Energy grid
    INTEGER, INTENT(in) :: len_mesh
    !! Length of meshgrid
    INTEGER, INTENT(in) :: meshin
    !! Current mesh
    INTEGER, INTENT(in) :: n
    !! Type of approximation for occupation number
    !
    ! Local variables
    !    
    INTEGER :: i,ii,j,c,v,pool_id,ab,l,leng,suc,nk
    !! integers used for do loops and if statements
    INTEGER :: cv,k_init,k_final,check
    !! integers used for c-c' and v-v' transitions and k states
    INTEGER :: imode
    !! phonon mode
    INTEGER :: k
    !! k-point
    INTEGER :: mesh
    !! QD mesh
    INTEGER :: ikk
    !! k index
    INTEGER :: ikq
    !! k+q index
    REAL(KIND = DP) :: enk
    !! k point K.S energy
    REAL(KIND = DP) :: emkq
    !! k+q point energy
    REAL(KIND = DP) :: ekkvb,ekqvb
    !! k, k+q energy for valence band
    REAL(KIND = DP) :: ekkcb,ekqcb
    !! k, k+q energy for conduction band
    REAL(KIND = DP) :: final_E
    !! final state energy
    REAL(KIND = DP) :: initial_E, weightd
    !! Initial state energy and DOS
    REAL(KIND = DP) :: nqv(nmodes),nq
    !! phonon occupation
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! Gaussian weight function 
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: xt(3)
    !! Buffer array
    ! 
    ii = 1
    xt = 0.d0
    cv = 0
    IF (meshin == 1) THEN
      leng = len_mesh-1
    ELSE
      leng = 1
    ENDIF

    DO mesh = 1, leng !len_mesh-1
      ii = maxdim
      DO k = 1, nkf
        !
        ikk = 2 * k-1
        ikq = ikk+1
        k_init = -1
        k_final = -1
        CALL kpmq_map(xkf(:,ikk), xqf(:,iq), 1, k_final)
        CALL kpmq_map(xkf(:,ikk), xt, 1, k_init)
        ! 
        DO c = 1, nbndfst
          IF (etf(ibndmin - 1 + c, ikk) <= fermi_energy) CYCLE
          !
          DO v = 1, nbndfst
            IF (etf(ibndmin - 1 + v, ikk) > fermi_energy) CYCLE
            !
            cv = 20
            i = v
            DO j = 1, nbndfst
              IF (etf(ibndmin - 1 + j, ikk) > fermi_energy) THEN
                cv = 10
                i = c
              ENDIF
              !
              emkq = etf(ibndmin - 1 + j, ikq) - fermi_energy
              ekkcb = etf(ibndmin - 1 + c, ikk) - fermi_energy
              ekkvb = etf(ibndmin - 1 + v, ikk) - fermi_energy
              ekqcb = etf(ibndmin - 1 + c, ikq) - fermi_energy
              ekqvb = etf(ibndmin - 1 + v, ikq) - fermi_energy
              IF ((ABS(ekkcb) > fsthick) .OR. (ABS(ekqcb) > fsthick)) CYCLE
              !
              DO imode = 1, nmodes
                IF (wf(imode, iq) < eps_acoustic) CYCLE
                DO ab = 1, 2
                  nqv(imode) = wgauss(-wf(imode,iq) / gtemp(1), -99)
                  nqv(imode) = nqv(imode) / (one - two * nqv(imode))

                  IF (n == -1) THEN
                    nq = nqv(imode)
                  ELSEIF (n == -2) THEN
                    nq = FLOOR(nqv(imode))
                  ENDIF

                  IF (cv==10) THEN
                    ! c->c' type transition |vk>|ck>   -> |vk>|ck+q>
                    initial_E = ekkcb-ekkvb
                    final_E = ekkcb-ekkvb-ekkcb+emkq
                  ELSE
                    ! v->v' type transtion  |vk-q>|ck> -> |vk>|ck>
                    initial_E = ekqcb-emkq
                    final_E = ekqcb-emkq+emkq-ekkvb
                  ENDIF

                  IF (ab==1) THEN
                    ! Emission
                    final_E = final_E+wf(imode,iq)
                    nq = nq+1
                  ELSE
                    ! Absorption
                    final_E = final_E-wf(imode,iq)
                    nq = nq-1
                  ENDIF
                  !   
                  pool_id = my_pool_id+1
                  !
                  IF (meshin == 1) THEN
                    weightd = (w0gauss((initial_E-E_mesh(mesh)) / degaussw, 0) / degaussw) * &
                      (w0gauss((final_E-E_mesh(mesh)) / degaussw, 0) / degaussw)

                    Energy(mesh,pool_id, my_image_id + 1, 2) = Energy(mesh,pool_id, my_image_id + 1, 2) + &
                                                               weightd * (wkf(ikk) / 2.0) * wqf(iq)* (1.0/omega)
                    weightd = (w0gauss((initial_E-E_mesh(mesh)) / degaussw, 0) / degaussw)

                    Energy(mesh,pool_id, my_image_id + 1, 3) = Energy(mesh,pool_id, my_image_id + 1, 3) + & 
                                                               weightd * (wkf(ikk) / 2.0) * (1.0/omega)
                    !WRITE(stdout,'(a,E22.14)') 'Energy',Energy(mesh,pool_id,2)
                  ENDIF
                  check = -1
                  suc = 0
                  IF (leng == len_mesh-1) THEN
                    IF ((((initial_E > E_mesh(mesh)) .AND. (initial_E <= E_mesh(mesh + 1))) &
                      .AND. ((final_E <= E_mesh(mesh + 1)) .AND. (final_E > E_mesh(mesh))))) THEN
                      suc = 1
                      IF (mesh == meshin) THEN
                        suc = 2
                      ENDIF
                    ENDIF
                  ENDIF
                  ! S. Tiwari: This If statement is experimental for a future
                  ! Future implementation of adaptive grid
                  !
                  IF ((initial_E > E_mesh(meshin)) .AND. (initial_E <= E_mesh(meshin + 1))) THEN
                    IF ((final_E <= (E_mesh(meshin + 1))) .AND. (final_E > (E_mesh(meshin)))) THEN
                      !  
                      IF (check == -1) THEN
                        suc = 2
                      ELSE
                        WRITE(stdout, '(5x,a)') 'Found_match'
                      ENDIF
                      ! 
                    ELSE
                      suc = 0
                    ENDIF
                  ENDIF
                  IF (suc > 0) THEN
                    IF (meshin == 1) THEN
                      Energy(mesh, pool_id, my_image_id + 1, 1) = Energy(mesh, pool_id, my_image_id + 1, 1) + 1
                      !selecq_QD(mesh, iq + 1) = 1
                    ENDIF
                    IF ((meshin > 1) .AND. (suc == 2)) THEN
                      !
                      H_quad(ii, 1) = ekkcb      !ABS(etf(ibndmin-1+i,ikq)-etf(ibndmin-1+i,ikk))+H_quad(iq,1,1)
                      H_quad(ii, 2) = ekkvb
                      H_quad(ii, 3) = ekqcb
                      H_quad(ii, 4) = c
                      H_quad(ii, 5) = j
                      H_quad(ii, 6) = imode
                      !
                      IF ((cv == 10)) THEN
                        H_quad(ii, 7) = epf17(j, c, imode, k) &
                                   * DSQRT(wqf(iq)) * (1.0 / (DSQRT(2 * (wf(imode, iq)))))
                      ELSE
                        H_quad(ii, 7) = 1.0*CONJG(epf17(v, j, imode, k)) &
                                   * DSQRT(wqf(iq)) * (1.0 / (DSQRT(2 * (wf(imode, iq)))))
                      ENDIF

                      H_quad(ii, 8) = v
                      H_quad(ii, 9) = emkq
                      IF (ab == 1) THEN
                        H_quad(ii, 10) = cv+0
                      ELSE
                        H_quad(ii, 10) = cv+1
                      ENDIF
                      H_quad(ii, 11) = k_init
                      H_quad(ii, 12) = k_final
                      H_quad(ii, 14) = iq
                      maxdim = maxdim+1
                      ii = ii+1
                      !
                    ENDIF ! suc == 2 and meshin > 1
                  ENDIF !suc >0 
                ENDDO ! absorption or emission
              ENDDO ! imode 
            ENDDO ! virtual band c or v
          ENDDO ! valence band 
        ENDDO ! conduction band
      ENDDO ! kpoint
    ENDDO! mesh
    !WRITE(stdout, '(a)' ) 'Done quasi'
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE build_quasi_eig
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE build_quasi_hamiltonian(E_mesh, len_mesh, mesh, n)
    !---------------------------------------------------------------------------
    !! This subroutine builds the quasi-degenrate Hamiltonian using H_quad,
    !! H_quad is distributed over pools, each H_quad is first pulled into H_temp
    !! This H_temp(ii,:) is them pulled into H_mat(index) (which keeps interaction) for QD
    !! states H_ind1(index),H_ind2(index). H_ind1(index) and H_ind1(index)
    !! returns the index of initial and final states which are stored in
    !! index_buf. Full description of states of QD states
    !! |vk>|ck+q>|+-nqv> is kept in index_buf which is in single pool. This
    !! index_buf(tot,6) is then distributed to different pools into index_buf_pool
    !! We then decide if we do OMP diagonalization or ELPA depending on the size
    !! of tot and allocated processors.
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iuindabs
    USE io_global,     ONLY : stdout, ionode_id
    USE modes,         ONLY : nmodes
    USE control_flags, ONLY : iverbosity
    USE input,         ONLY : nstemp, fermi_energy, nkf1, nkf2, nkf3,    &
                              mp_mesh_k, filkf
    USE global_var,    ONLY : etf, ibndmin, gtemp, nkf, nqtotf, wf, wqf, & 
                              nbndfst, epf17, nkqtotf, xkf, nkqf,        &
                              bztoibz, H_quad, Energy, stop_qdabs,       &
                              Eigenvec_alloc_pool, Eigenval,do_elpa,     &
                              H_temp, size_m, H_mat, H_ind1, r_tot,      &
                              tot_pool, index_buf_pool, Eigenvec,        &
                              index_buf, tot, totcv, size_tot, H_ind2,   &
                              indtot, indtotf, indtotcv, indtoti, n_q
    USE mp_global,     ONLY : my_pool_id, npool, inter_pool_comm,        &
                              inter_image_comm
    USE mp,            ONLY : mp_bcast,mp_barrier,mp_sum
    USE cell_base,     ONLY : omega, at, bg
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two,      &
                              zero, pi, ci, eps8, eps20, czero, eps40
    USE bzgrid,        ONLY : kpmq_map
    USE mp_images,     ONLY : my_image_id, nimage
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: E_mesh(len_mesh)
    !! The enegy mesh-grid
    INTEGER, INTENT(in) :: len_mesh
    !! Length of the meshgrid
    INTEGER, INTENT(in) :: mesh
    !! Current mesh point for QDPT
    INTEGER, INTENT(in) :: n
    !! Approximation on occupation number
    !
    ! Local variables 
    !   
    INTEGER :: i, j, l, iq, totf
    !! Counter on eigenstates and total states
    INTEGER :: c, v, cv, t, nk
    !! Counter on conduction and valence states
    !! conduction or valence absorption, number of
    !! k states
    INTEGER :: imode, info
    !! mode  and a buffer integer for ZGEEV
    INTEGER :: k, rest
    !! kpoint and distribution of k point
    INTEGER :: ii, poolf, pooli
    !! counter on all states, final and initial state pool location
    INTEGER :: k_init, k_final
    !! Initial k and final k+q point
    INTEGER :: ierr, size_temp, ierr2
    !! Error counters and size of temporary buffer Hamiltonian
    INTEGER :: suc1, suc2, ipool, index1, index2, nktotf
    !! If many-body states found in current states (suc1,suc2)
    !! location of current many-body states in matrices and total number of k
    !! points
    INTEGER :: numc, numv, numq, numci, numvi
    !! final state and initial state positions
    INTEGER :: image
    !! counter over images
    REAL(KIND = DP) :: enk, enkq, emkq, ekkcb, ekkvb, ekqcb, final_e, Num, nq
    !! KS energies, final state position, occupation, and type of occupation
    REAL(KIND = DP) :: initial_E, start, finish, tot_pool_a(npool)
    !! k point, initial state index, total states in the pool
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Bose-Einstein distribution
    COMPLEX(KIND = DP) :: inter
    !! Electron-phonon interaction
    !
    do_elpa = 0
    r_tot   = 0
    H_mat   = czero
    H_ind1  = 0
    H_ind2  = 0
    suc1    = 0
    suc2    = 0
    tot     = 0
    totcv   = 0
    totf    = 0
    index1  = 0
    index2  = 0
    k_init  = 0
    k_final = 0
    nq      = n
    IF (filkf /=' ') THEN
      nktotf = nkqtotf / 2
    ELSE
      nktotf = nkf1 * nkf2 * nkf3
    ENDIF
    cv = 0
    Num = -1.d0
    IF ((my_pool_id == ionode_id) .AND. (my_image_id == 0)) THEN
      ALLOCATE(indtot(size_m, 5), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_quasi', 'Error allocating indtot(size_m,5)', 1)
      !  
      ALLOCATE(indtotf(size_m), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_quasi', 'Error allocating indtotf(size_m)', 1)
      !   
      ALLOCATE(indtotcv(size_m), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_quasi', 'Error allocating indtotcv(size_m)', 1)
      !
      ALLOCATE(indtoti(size_m,3), STAT = ierr)
      IF (ierr /= 0) CALL errore('qdabs', 'Error allocating indtoti(size_m,3)', 1)
      !
      suc1    = 0
      suc2    = 0
      tot     = 0
      totcv   = 0
      totf    = 0
      index1  = 0
      index2  = 0
      k_init  = 0
      k_final = 0
      nq      = n
      nktotf  = nkf1 * nkf2 * nkf3
      IF (filkf /= ' ') THEN
        nktotf = nkqtotf / 2
      ELSE
        nktotf = nkf1 * nkf2 * nkf3
      ENDIF
      !  
      cv = 0
      Num = -1.d0
      !  
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    CALL mp_barrier(inter_image_comm) 
    DO image = 1, nimage
      !H_temp(:, :) = 0.d0
      CALL mp_barrier(inter_image_comm)
      DO ipool = 1, npool
!        H_temp(:, :) = 0.d0
        size_temp = 0
        CALL mp_barrier(inter_pool_comm)
        CALL mp_barrier(inter_image_comm) 
        IF ((ipool == my_pool_id + 1) .AND. (image == my_image_id + 1)) THEN
          size_temp = INT(H_quad(1, 13))
          IF ((size_temp > 0)) THEN
            H_temp(1:size_temp, 1:12) = H_quad(:, 1:12)
            H_temp(1:size_temp, 13) = H_quad(:, 14)
          ENDIF
        ENDIF
        ! S.T: Here for a future implementation
        !
        !CALL mp_barrier(inter_pool_comm)
        !CALL mp_barrier(inter_image_comm)
        !CALL mp_sum(H_temp, inter_pool_comm)
        !CALL mp_sum(H_temp, inter_image_comm)
        !CALL mp_sum(size_temp, inter_pool_comm)
        !CALL mp_sum(size_temp, inter_image_comm)
        !CALL mp_barrier(inter_pool_comm)
        !CALL mp_barrier(inter_image_comm)
        !
        !CALL mp_barrier(inter_image_comm)
        !CALL mp_bcast(H_temp, image - 1, inter_image_comm)
        !CALL mp_bcast(size_temp, image - 1, inter_image_comm)
        !CALL mp_barrier(inter_image_comm)
        !
        CALL mp_barrier(inter_pool_comm)
        CALL mp_bcast(H_temp, ipool - 1, inter_pool_comm)
        CALL mp_bcast(size_temp, ipool - 1, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
    
        IF ((my_pool_id == ionode_id) .AND. (size_temp > 0) .AND. (my_image_id == 0)) THEN
          DO ii =1, size_temp
            IF ((INT(H_temp(ii, 10)) > 0 ) .AND. (INT(H_temp(ii, 10)) < 22 )) THEN
              ekkcb = REAL(H_temp(ii, 1))
              ekkvb = REAL(H_temp(ii, 2))
              ekqcb = REAL(H_temp(ii, 3))
              c = INT(H_temp(ii, 4))
              j = INT(H_temp(ii, 5))
              imode = INT(H_temp(ii, 6))
              inter = H_temp(ii, 7)
              v = INT(H_temp(ii, 8))
              emkq = REAL(H_temp(ii, 9))
              cv = INT(H_temp(ii, 10))
              k_init = INT(H_temp(ii, 11))
              k_final = INT(H_temp(ii, 12))
              iq = INT(H_temp(ii, 13))
              IF (iverbosity == 5) THEN
                !  
                IF ((cv > 21) .OR. (cv < 10)) THEN
                  !
                  WRITE(stdout,'(/5x,a,9I10,5E22.14)') 'Error:cv,c,v,j,k_init,k_final,ipool,&
                           iq,ii,ekkcb,ekkvb,ekmq,ekqcb', &
                          cv, c, v, j, k_init, k_final, ipool, iq, &
                          ii, ekkcb * ryd2ev, ekkvb * ryd2ev, emkq * ryd2ev, &
                          ekqcb * ryd2ev, (ekqcb-ekkvb) * ryd2ev
                ENDIF
              ENDIF
              IF (n == -1) THEN
                !
                nq = wgauss(-wf(imode,iq) / gtemp(1), -99)
                nq = nq / (one - two * nq)
              ELSEIF (n == -2) THEN
                !  
                nq = wgauss(-wf(imode,iq) / gtemp(1), -99)
                nq = FLOOR( nq / (one - two * nq))
              ELSEIF (n == -3) THEN
                !
                nq = n_q(iq, imode)
              ENDIF
              !
              IF (cv == 10) THEN
                !  
                final_e = ekkcb - ekkvb - ekkcb + emkq + wf(imode,iq)
                initial_E = ekkcb - ekkvb
                Num = nq + 1
              ELSEIF (cv == 20) THEN
                !
                final_e = ekqcb - emkq + emkq - ekkvb + wf(imode,iq)
                initial_E = ekqcb - emkq
                Num = nq + 1
              ELSEIF (cv == 11) THEN
                !
                final_e = ekkcb - ekkvb - ekkcb + emkq - wf(imode,iq)
                initial_E = ekkcb - ekkvb
                Num = nq
              ELSEIF (cv == 21) THEN
                !
                final_e = ekqcb - emkq + emkq - ekkvb - wf(imode,iq)
                initial_E = ekqcb - emkq
                Num = nq
              ENDIF
              !
              IF (Num < 1) THEN
                !
                suc1 = -1
                suc2 = -1
              ENDIF
              !    
              inter = inter*SQRT(Num)
              !
              IF (Num == nq+1) THEN
                !
                Num = nq+1
              ELSEIF (Num == nq) THEN
                !  
                Num = nq-1
              ENDIF
              IF ((final_e < E_mesh(mesh)) .OR. (final_e > E_mesh(mesh + 1)))  THEN
                !
                WRITE(stdout,'(/5x,a,3E22.14)')'final_e', final_e*ryd2ev, E_mesh(mesh)*ryd2ev,     &
                                                          E_mesh(mesh+1)*ryd2ev
              ENDIF
              IF ((initial_E < E_mesh(mesh)) .OR. (initial_E > E_mesh(mesh+1)))  THEN
                !
                WRITE(stdout,'(/5x,a,3E22.14)')'initial_E', initial_E*ryd2ev, E_mesh(mesh)*ryd2ev, &
                                                            E_mesh(mesh+1)*ryd2ev
              ENDIF
              !
              IF ((cv == 10) .OR. (cv == 11)) THEN
                !
                numci = (c-1)*nktotf+k_init
                numvi = (v-1)*nktotf+k_init
                numc = (j-1)*nktotf+k_final
                numv = (v-1)*nktotf+k_init
                numq = (imode-1)*nqtotf+iq
                poolf = pool_final(k_init)!  ipool
                pooli = pool_final(k_init)!ipool
                !
                DO l = 1, totcv
                  IF ((numci == INT(indtoti(l,1))) .AND. (numvi == INT(indtoti(l,2)))) THEN
                    IF (Num >= 0) THEN
                       index1 = INT(indtotcv(l))
                       suc1 = 1
                    ENDIF
                  ENDIF
                ENDDO
                DO l = 1, totf
                  IF ((ABS(numc-INT(indtot(l, 1))) < eps8) .AND. &
                      (ABS(numv-INT(indtot(l, 2))) < eps8) .AND. &
                      (ABS(numq-INT(indtot(l, 3))) < eps8) .AND. & 
                      (ABS(Num-indtot(l, 4)) < eps8)) THEN
                    !    
                    suc2 = 1
                    index2 = INT(indtotf(l))
                  ENDIF
                ENDDO
              ELSE
                numci = (c - 1) * nktotf + k_final
                numvi = (j - 1)* nktotf + k_final
                numc = (c - 1) * nktotf + k_final
                numv = (v - 1) * nktotf + k_init
                numq = (imode - 1) * nqtotf + iq
                !
                poolf = pool_final(k_init)
                pooli = pool_final(k_final)
                !
                DO l = 1, totcv
                  IF ((numci == INT(indtoti(l, 1))) .AND. &
                      (numvi == INT(indtoti(l, 2)))) THEN
                    IF (Num >= 0) THEN
                      !    
                      index1 = INT(indtotcv(l))
                      suc1 = 1
                    ENDIF
                  ENDIF
                ENDDO
                DO l = 1, totf
                  !
                  IF ((ABS(numc-INT(indtot(l, 1))) < eps8) .AND. &
                     (ABS(numv-INT(indtot(l, 2))) < eps8) .AND. &
                     (ABS(numq-INT(indtot(l, 3))) < eps8) .AND. &
                     (ABS(Num-indtot(l, 4)) < eps8)) THEN
                    ! 
                    suc2 = 1
                    index2 = INT(indtotf(l))
                  ENDIF
                ENDDO
              ENDIF
     
              IF (suc1 == 1) THEN
                IF (suc2 == 1) THEN
                  r_tot = r_tot + 1
                  H_mat(r_tot) = inter
                  H_ind1(r_tot) = index1
                  H_ind2(r_tot) = index2
                  !
                  r_tot = r_tot + 1
                  H_mat(r_tot) = CONJG(inter)
                  H_ind1(r_tot) = index2
                  H_ind2(r_tot) = index1
                ELSEIF (suc2 == 0) THEN
                  index2 = tot + 1
                  indtotf(totf + 1) = index2
                  r_tot = r_tot + 1
                  H_mat(r_tot) = inter
                  H_ind1(r_tot) = index1
                  H_ind2(r_tot) = index2
                  !
                  r_tot = r_tot + 1
                  H_mat(r_tot) = CONJG(inter)
                  H_ind1(r_tot) = index2
                  H_ind2(r_tot) = index1
                  !
                  r_tot = r_tot + 1
                  H_mat(r_tot) = final_e!CONJG(inter)
                  H_ind1(r_tot) = index2
                  H_ind2(r_tot) = index2
                  !
                  indtot(totf + 1, 1) = numc
                  indtot(totf + 1, 2) = numv
                  indtot(totf + 1, 3) = numq
                  indtot(totf + 1, 4) = Num
                  indtot(totf + 1, 5) = poolf
                  tot = tot + 1
                  totf = totf + 1
                ENDIF
              ELSEIF (suc1 == 0) THEN
                IF (suc2 == 1) THEN
                  index1 = tot + 1
                  indtotcv(totcv + 1) = index1
                  r_tot = r_tot + 1
                  H_mat(r_tot) = inter
                  H_ind1(r_tot) = index1
                  H_ind2(r_tot) = index2
                  !
                  r_tot = r_tot + 1
                  H_mat(r_tot) = CONJG(inter)
                  H_ind1(r_tot) = index2
                  H_ind2(r_tot) = index1
                  !
                  r_tot = r_tot + 1
                  H_mat(r_tot) = initial_E!CONJG(inter)
                  H_ind1(r_tot) = index1
                  H_ind2(r_tot) = index1
                  !
                  indtoti(totcv + 1, 1) = numci
                  indtoti(totcv + 1, 2) = numvi
                  indtoti(totcv + 1, 3) = pooli
                  totcv = totcv + 1
                  tot = tot + 1
                ELSEIF (suc2 == 0) THEN
                  index1 = tot + 1
                  index2 = tot + 2
                  indtotcv(totcv + 1) = index1
                  indtotf(totf + 1) = index2
                  ! 
                  r_tot = r_tot + 1
                  H_mat(r_tot) = inter
                  H_ind1(r_tot) = index1
                  H_ind2(r_tot) = index2
                  !
                  r_tot = r_tot + 1
                  H_mat(r_tot) = CONJG(inter)
                  H_ind1(r_tot) = index2
                  H_ind2(r_tot) = index1
                  !
                  r_tot = r_tot + 1
                  H_mat(r_tot) = initial_E!CONJG(inter)
                  H_ind1(r_tot) = index1
                  H_ind2(r_tot) = index1
                  !
                  r_tot = r_tot + 1
                  H_mat(r_tot) = final_e!CONJG(inter)
                  H_ind1(r_tot) = index2
                  H_ind2(r_tot) = index2
                  !
                  indtot(totf + 1, 1) = numc
                  indtot(totf + 1, 2) = numv
                  indtot(totf + 1, 3) = numq
                  indtot(totf + 1, 4) = Num
                  indtot(totf + 1, 5) = poolf
                  indtoti(totcv + 1, 1) = numci
                  indtoti(totcv + 1, 2) = numvi
                  indtoti(totcv + 1, 3) = pooli
                  tot = tot + 2
                  totcv = totcv + 1
                  totf = totf + 1
                ENDIF
              ENDIF
            ENDIF
            suc1 = 0
            suc2 = 0
            !  
          ENDDO ! size_temp
        ENDIF ! IF this pool has any interactions
        CALL mp_barrier(inter_pool_comm)
      ENDDO ! ipool
      CALL mp_barrier(inter_image_comm)
    ENDDO 
    !
    CALL mp_barrier(inter_pool_comm)
    CALL mp_barrier(inter_image_comm) 
    !
    CALL mp_bcast(r_tot, ionode_id, inter_image_comm)
    CALL mp_bcast(H_mat, ionode_id, inter_image_comm)
    CALL mp_bcast(H_ind1, ionode_id, inter_image_comm)
    CALL mp_bcast(H_ind2, ionode_id, inter_image_comm)
    !
    CALL mp_bcast(r_tot, ionode_id, inter_pool_comm)
    CALL mp_bcast(H_mat, ionode_id, inter_pool_comm)
    CALL mp_bcast(H_ind1, ionode_id, inter_pool_comm)
    CALL mp_bcast(H_ind2, ionode_id, inter_pool_comm)
    !
    IF ((my_pool_id == ionode_id) .AND. (my_image_id == ionode_id)) THEN
      ! Many body Hamiltonian built
      ! Diagonalization of many body Hamiltonian
      !  
      WRITE(stdout, '(/5x,a,I10)') 'matrix rank =', tot
      !  
      IF (tot == 0) THEN
        WRITE(stdout, '(/5x,a)') 'no states'
      ENDIF
      !  
      DO i = 1, totcv
        index_buf(INT(indtotcv(i)), 6)  = indtoti(i, 3)
        index_buf(INT(indtotcv(i)), 5)  = indtotcv(i)
        index_buf(INT(indtotcv(i)), 1:2) = indtoti(i, 1:2)
      ENDDO
      !  
      DO i = 1, totf
        index_buf(INT(indtotf(i)), 6) = indtot(i, 5)
        index_buf(INT(indtotf(i)), 5) = indtotf(i)
        index_buf(INT(indtotf(i)), 1:4) = indtot(i, 1:4)
      ENDDO
      !  
      suc1 = 0
      !  
      IF (suc1 == 0) THEN
        IF (iverbosity == 5) THEN
          WRITE(stdout,'(/5x,a)') 'Matrix is Hermitian'
        ENDIF          
      ENDIF
      ! FIX in future  
      DEALLOCATE(indtot)
      DEALLOCATE(indtotf)
      DEALLOCATE(indtotcv)
      DEALLOCATE(indtoti)
      !  
    ENDIF
    !
    CALL mp_barrier(inter_pool_comm)
    CALL mp_barrier(inter_image_comm)
    !
    CALL mp_bcast(do_elpa, ionode_id, inter_image_comm)
    CALL mp_bcast(tot, ionode_id, inter_image_comm)
    CALL mp_bcast(totf, ionode_id, inter_image_comm)
    CALL mp_bcast(totcv, ionode_id, inter_image_comm)
    CALL mp_bcast(index_buf, ionode_id, inter_image_comm)
    !
    CALL mp_bcast(do_elpa, ionode_id, inter_pool_comm)
    CALL mp_bcast(tot, ionode_id, inter_pool_comm)
    CALL mp_bcast(totf, ionode_id, inter_pool_comm)
    CALL mp_bcast(totcv, ionode_id, inter_pool_comm)
    CALL mp_bcast(index_buf, ionode_id, inter_pool_comm)
    !
    CALL mp_barrier(inter_pool_comm)
    !
    tot_pool_a = 0.d0
    size_tot = 0
    DO i = 1, tot
      !  
      tot_pool_a(INT(index_buf(i, 6))) = tot_pool_a(INT(index_buf(i, 6))) + 1.0
      !  
    ENDDO
    CALL mp_barrier(inter_pool_comm)
    !
    size_tot = INT(MAXVAL(tot_pool_a))
    !
    WRITE(stdout, '(/5x,a,2I10)') 'Size of Eigenvec:', size_tot, tot
    CALL mp_barrier(inter_pool_comm)
    WRITE(stdout, '(/5x,a,I10)') 'Building Hamiltonian with total elements:', r_tot
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE build_quasi_hamiltonian
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE diag_quasi_hamiltonian(E_mesh, len_mesh, mesh, n)
    !---------------------------------------------------------------------------
    !! This subroutine performs the diagonalization of the QD Hamiltonian. We
    !! save QD Hamiltonian as  sparse matrix H_QD(i,j) as H_mat(index)
    !! where H_ind1(index) .EQ. i and H_ind2(index) .EQ. j. With this we do two
    !! types of diagonalization. (1) OMP, where we simply build Hamiltonian and
    !! solve in each core using ZGEEV while for (2) ELPA, we build the Hamiltonian
    !! block in each core. This is performed using SCALAPACK which allows (i,j)
    !! to be converted to (i_loc,j_loc)_proc where proc is the processor and
    !! i_loc, j_loc are the local indices for that processor block. Once we
    !! finish diagonalization, the Eigenvectors are redistributed from processor
    !! local to gobal index, and then leading dimension that belongs to |vk>|ck+q> is
    !! redistributed to their respective pools using
    !! index_buff_pool(index,5):-> poolf
    !    
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iuindabs
    USE io_global,     ONLY : stdout,ionode_id
    USE modes,         ONLY : nmodes
    USE control_flags, ONLY : iverbosity
    USE input,         ONLY : nstemp, fermi_energy, nkf1,                               &
                              nkf2, nkf3, mp_mesh_k
    USE global_var,    ONLY : etf, ibndmin, gtemp, nkf, nqtotf, wf, wqf, nbndfst, epf17,&
                              nkqtotf, xkf, nkqf, bztoibz, Energy, stop_qdabs,          &
                              Eigenvec_alloc_pool, Eigenval, do_elpa, size_m, r_tot,    &
                              tot_pool, index_buf_pool, Eigenvec, index_buf, tot,       &
                              totcv, H_mat, H_ind1, H_ind2, n_q
    USE mp_global,     ONLY : my_pool_id, npool, inter_pool_comm, inter_image_comm
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE cell_base,     ONLY : omega, at, bg
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi,ci, eps8,  &
                              eps20, czero, eps40
    USE bzgrid,        ONLY : kpmq_map
    USE mp_images,     ONLY : my_image_id
#if defined(__ELPA)
    USE elpa
    !
    class(elpa_t), pointer :: elp
    !! Elpa API 
#endif
    REAL(KIND = DP), INTENT(in) :: E_mesh(len_mesh)
    !! Energy mesh-grid for QDPT    
 
    INTEGER, INTENT(in) :: len_mesh
    !! Length of meshgrid for QDPT
    INTEGER, INTENT(in) :: mesh
    !! Current mesh grid
    INTEGER, INTENT(in) :: n
    !! Occupation number approximation
    INTEGER :: i, iq, totx, toty, totn
    !! Counter on index and distribution of arrays over cores
    INTEGER :: t, nk, l
    !! Counter on indices of final_many-body states
    INTEGER :: size_tot
    !! total size of buffer Hamiltonian
    INTEGER :: ierr, size_temp, ierr2
    !! error flags and temporary buffer size
    REAL(KIND = DP) :: start, finish
    !! estimates cpu time taken to finish diagonalization and distribution
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Bose-Einstein distribution
#if defined(__ELPA)
    INTEGER, EXTERNAL                    :: numroc, INDXG2L, INDXG2P, INDXL2G
    !! SCALAPACK variables
    REAL(KIND = c_double), ALLOCATABLE   :: ev(:), RWORK(:)
    !! Eigenvaues and RWORK for ZGEEV
    COMPLEX(KIND = c_double) :: inter
    !! Interaction strength e-p
    COMPLEX(KIND = c_double),ALLOCATABLE :: al(:,:), zl(:,:), WORK(:), &
                                            buff(:,:), H_f(:,:)
    !! Hamiltonian distributed over processors, eigenvectors, ZGEEV Work,
    !! Buffer array, full Hamiltonian
#else
    REAL(KIND = DP), ALLOCATABLE         :: ev(:), RWORK(:)
    !! Eigenvaues and RWORK for ZGEEV
    COMPLEX(KIND = DP),ALLOCATABLE       :: al(:,:), zl(:,:), WORK(:), &
                                            buff(:,:), H_f(:,:)
    !! Hamiltonian distributed over processors, eigenvectors, ZGEEV Work,
    !! Buffer array, full Hamiltonian
#endif
    ! --------------------------------------------------------------------------------!!
    !! All variables are for diagonalization
    CHARACTER(len=8)                   :: task_suffix
    !! ELPA related task 
    INTEGER                            :: success
    !! Success for ELPA diagonalization or not
    INTEGER                            :: nblk, info
    !! Block length for cyclic distribution
    INTEGER                            :: np_rows, np_cols, na_rows, na_cols
    !! Processor distribution in rows and columns
    INTEGER                            :: my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
    !! Present processor row and column, row/column communicators
    INTEGER                            :: my_blacs_ctxt, sc_desc(9), nprow,npcol,mpierr
    !! BLACS contexts, number of processor rows
    INTEGER                            :: j, il, jl, kl, kp, jp, jg, kg
    !! indices, l:local, g:global, p:processor
    INTEGER, ALLOCATABLE               :: buff_ind(:,:)
    !! BUffer array for redistribution
    INTEGER, ALLOCATABLE               :: na_rows_max(:),na_cols_max(:)
    !! maximum number of elements(row/column) in each processor
    INTEGER, ALLOCATABLE               :: loc(:)
    !! location in an array 
    INTEGER                            :: STATUS
    !! Status if diagonalization
    INTEGER, PARAMETER                 :: error_units = 0
    !! Error unit
    INTEGER                            :: size_vec
    !! Size of the leading dimension of the largest eigenvector
    !!---------------------------------------------------------------------------------!!
    nblk=8
    CALL cpu_time(start)
    WRITE(stdout, '(/5x,a)') 'Entered diag_quasi_Hamiltonian'
#if defined(__ELPA)
    size_vec = MINVAL((/2500 * CEILING(npool/200.d0), 6000/))
#else
    size_vec = tot + 1
#endif
    ALLOCATE(Eigenvec(size_vec, size_vec), STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_quasi', 'Error allocating Eigenvec(size_vec,size_vec)', 1)
    IF (iverbosity == 5) WRITE(stdout, '(/5x,a)') ' Allocated Eigenvec'
    !
    Eigenvec = 0.d0
    Eigenval = 0.d0
    !
    IF ((tot > 1) .AND. (tot < size_vec) .AND. &
        (my_pool_id == ionode_id)) THEN

      ALLOCATE(H_f(tot, tot), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error allocating H_f(tot,tot)', 1)
      ALLOCATE(WORK(size_m * 20), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error allocating WORK(20*size_m)', 1)
      ALLOCATE(RWORK(size_m * 2) ,STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error allocating RWORK(size_m)', 1)
      H_f = 0.d0
      DO i = 1, r_tot
        !
        H_f(H_ind1(i), H_ind2(i)) = H_mat(i)
        !
      ENDDO
      !
      CALL ZGEEV('N', 'V', tot, H_f(1:INT(tot),1:INT(tot)), INT(tot), Eigenval(1:INT(tot)),&
           ' ', INT(tot), Eigenvec(1:INT(tot),1:INT(tot)), INT(tot), WORK, INT(size_m*20), &
          RWORK, info)
      !
      IF (info /= 0) WRITE(stdout, '(/5x,a)') 'Failed diagonalization'
      !
      WRITE(stdout, '(/5x,a)') 'LAPACK diagonalization complete'
      IF (iverbosity == 5) THEN
        WRITE(stdout, '(/5x,a,2I10)') 'Finished Eigenvec calculation LAPACK with total  &
                            &interactions: ', r_tot
      ENDIF
      !
      DEALLOCATE(H_f, STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error deallocating H_f', 1)
      DEALLOCATE(WORK)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error deallocating H_f', 1)
      DEALLOCATE(RWORK)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error deallocating H_f', 1)
      !
      do_elpa = 0
      !  
    ELSEIF (tot >= size_vec) THEN
      do_elpa = 1
    ENDIF
    CALL mp_bcast(do_elpa,ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    tot_pool = 0
    t = 1
    IF (do_elpa == 0) THEN
      WRITE(stdout, '(/5x,a)') 'Used LAPACK for diagonalization'
      !
      CALL mp_bcast(Eigenvec, ionode_id, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !  
      DO i = 1, tot
        IF (my_pool_id == INT(index_buf(i, 6)) - 1) THEN
          ! 
          Eigenvec_alloc_pool(t, 1:tot) = Eigenvec(i, 1:tot)
          index_buf_pool(t, 1:4) = index_buf(i, 1:4)
          index_buf_pool(t, 5) = t
          index_buf_pool(t, 6) = index_buf(i, 6)
          t = t + 1
          ! 
        ENDIF
      ENDDO
      CALL mp_barrier(inter_pool_comm)
      !  
      tot_pool = t - 1
      CALL mp_sum(Eigenval, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
#if defined(__ELPA)
    ELSEIF (do_elpa == 1) THEN
      !  
      ! Allocations
      !
      ALLOCATE(na_rows_max(npool), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error allocating na_rows_max(npool)', 1)
      ALLOCATE(na_cols_max(npool), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error allocating na_cols_max(npool)', 1)
      !
      WRITE(stdout, '(/5x,a)'),'Using ELPA for diagonalization'
      DO np_cols = NINT(SQRT(REAL(npool))), 2, -1
        IF (MOD(npool, np_cols) == 0 ) exit
      ENDDO
      ! at the end of the above loop, npools is always divisible by np_cols
      np_rows = npool / np_cols
      totn = tot
      ! initialise BLACS
      my_blacs_ctxt = inter_pool_comm
      !
      CALL BLACS_Gridinit(my_blacs_ctxt, 'C', np_rows, np_cols)
      CALL BLACS_Gridinfo(my_blacs_ctxt, nprow, npcol, my_prow, my_pcol)
      ! 
      !
      IF ((my_pool_id == ionode_id) .AND. (iverbosity == 5)) THEN
        WRITE(stdout,'(/5x,a)'),'| Past BLACS_Gridinfo.'
      ENDIF
      ! determine the neccessary size of the distributed matrices,
      ! we use the scalapack tools routine NUMROC
      !
      na_rows = numroc(totn, nblk, my_prow, 0, np_rows)
      na_cols = numroc(totn, nblk, my_pcol, 0, np_cols)
      na_rows_max = 0
      na_cols_max = 0
      na_rows_max(my_pool_id + 1) = na_rows
      IF (na_cols > 0) THEN
        !
        na_cols_max(my_pool_id + 1) = na_cols
        !
      ENDIF
      CALL mp_barrier(inter_pool_comm)
      CALL mp_sum(na_rows_max, inter_pool_comm)
      CALL mp_sum(na_cols_max, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      totx = INT(maxval(na_rows_max))
      toty = INT(maxval(na_cols_max))
      !
      IF ((my_pool_id == ionode_id) .AND. (iverbosity == 5)) THEN
           WRITE(stdout, '(/5x,a,I10,I10,I10)'), 'na_rows,na_cols,my_prow', na_rows, na_cols, my_prow
           WRITE(stdout,'(/5x,a,I10,I10,I10)'), 'nprow,npcol,my_pcol', nprow, npcol, my_pcol
      ENDIF
      !  
      ! set up the scalapack descriptor for the checks below
      ! For ELPA the following restrictions hold:
      ! - block sizes in both directions must be identical (args 4 a. 5)
      ! - first row and column of the distributed matrix must be on
      !   row/col 0/0 (arg 6 and 7)
      !
      CALL descinit(sc_desc, totn, totn, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info)
      !  
      IF (info .NE. 0) THEN
        WRITE(error_units,*) 'Error in BLACS descinit! info=',info
        WRITE(error_units,*) 'Most likely this happend since you want to use'
        WRITE(error_units,*) 'more MPI tasks than are possible for your'
        WRITE(error_units,*) 'problem size (matrix size and blocksize)!'
        WRITE(error_units,*) 'The blacsgrid can not be set up properly'
        WRITE(error_units,*) 'Try reducing the number of MPI tasks...'
       ! call MPI_ABORT(inter_pool_comm, 1, mpierr)
      ENDIF
      !  
      ! Error is not assigned for these because there are some issus
      ALLOCATE(al (na_rows, na_cols))
      ALLOCATE(zl (na_rows, na_cols))
      ALLOCATE(ev (totn))
      al=(0.d0, 0.d0)
      zl=(0.d0, 0.d0)
      !
      ALLOCATE(buff_ind(npool, toty + 2), STAT = ierr)
      IF (ierr /= 0)CALL errore('diag_quasi', 'Error allocating buff_ind(tot,tot)', 1)
      ALLOCATE(buff(npool, toty), STAT = ierr)
      IF (ierr /= 0)CALL errore('diag_quasi', 'Error allocating buff(tot,tot)', 1)
      buff = 0.d0
      buff_ind = 0
      CALL mp_barrier(inter_pool_comm)
      WRITE(stdout, '(/5x,a)') 'Past initial array allocation'
      !  
      ! Distribute the Hamiltonian to different processors 
      !
      DO i = 1, r_tot
        !    
        jg = H_ind1(i)
        kg = H_ind2(i)
        !
        inter = H_mat(i)
        jl = INDXG2L(jg, nblk, 0, 0, nprow) !FLOOR((j-1)/REAL(nprow))!
        jp = INDXG2P(jg, nblk, 0, 0, nprow)   !mod(j-1,nprow)
        !
        kl = INDXG2L(kg, nblk, 0, 0, npcol)!FLOOR((k-1)/REAL(npcol))
        kp = INDXG2P(kg, nblk, 0, 0, npcol) !mod(k-1,npcol)
        IF ((jp == my_prow) .AND. (kp == my_pcol)) THEN
          !  
          IF ((jl <= na_rows) .AND. (kl <= na_cols)) THEN
            !
            al(jl, kl) = inter
            !
          ENDIF
        ENDIF
      ENDDO
      !
      CALL mp_barrier(inter_pool_comm)
      !
      IF (iverbosity == 5) WRITE(stdout, '(/5x,a)'),'Finished Hamiltonian setup'
      !  
      ! Start ELPA diagonalization 
      !  
      IF (elpa_init(20180501) /= elpa_ok) THEN
        WRITE(stdout, *), "ELPA API version not supported"
        STOP
      ENDIF
      elp => elpa_allocate()    
      !  
      ! set parameters decribing the matrix and it's MPI distribution
      CALL elp%set("na", totn, success)
      CALL elp%set("nev", totn, success)
      CALL elp%set("local_nrows", na_rows, success)
      CALL elp%set("local_ncols", na_cols, success)
      CALL elp%set("nblk", nblk, success)
      CALL elp%set("mpi_comm_parent", inter_pool_comm, success)
      CALL elp%set("process_row", my_prow, success)
      CALL elp%set("process_col", my_pcol, success)
      success = elp%setup()
      CALL elp%set("solver", elpa_solver_2stage, success)
      CALL elp%set("complex_kernel", elpa_solver_2stage, success)
      ! Calculate eigenvalues/eigenvectors
      !
      WRITE(stdout, '(/5x,a)')'| Entering one-step ELPA solver ... '
      !
      !  
      CALL mp_barrier(inter_pool_comm) ! for correct timings only
      CALL elp%eigenvectors(al, ev, zl, success)
      CALL mp_barrier(inter_pool_comm) ! for correct timings only
      !  
      WRITE(stdout, '(/5x,a)') '| One-step ELPA solver complete.'
      !  
      ! Redistribute the Eigenvector's leading dimension n of E(n,:) to different
      ! pools
      ! The pool is already decided when we build the Hamiltonian and resides
      ! in index_buf(:,6)
      ! The  Eigenvector's leading dimension m, E(:,m) are the eigenvectors
      ! which  also need to  be sorted and should remain the same for all
      ! distributed Eigenvectors
      ! The leading dimension m is sorted for each eigenvector and kept in
      ! buf_ind(pool_id,l)
      buff = 0.d0
      buff_ind = 0
      CALL mp_barrier(inter_pool_comm)
      t = 0
      ALLOCATE(loc(tot), STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error allocating loc', 1) 
      loc = 0
      l = 0
      DO jl = 1, totx
        IF (MOD(jl, 500) == 0) THEN
          WRITE(stdout,'(/5x, a, I10, a, I10)'),'Eigenvector re-distributed = ', jl, '/' , totx
        ENDIF
        DO kl = 1, toty
          IF ((jl <= na_rows) .AND. (kl <= na_cols)) THEN
            !
            kg = INDXL2G(kl, nblk, my_pcol, 0, npcol)
            jg = INDXL2G(jl, nblk, my_prow, 0, nprow)
            buff_ind(my_pool_id+1,1) = jg
            IF ((kg > 0) .AND. (kg <= tot)) THEN
              l = l+1
              buff_ind(my_pool_id + 1, l + 2) = kg
              buff(my_pool_id + 1, l) = zl(jl, kl)
            ENDIF
          ENDIF
        ENDDO
        buff_ind(my_pool_id + 1, 2) = l
        !
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(buff_ind, inter_pool_comm)
        CALL mp_sum(buff, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
        !
        DO i = 1, npool
          IF ((buff_ind(i,1) > 0) .AND. (buff_ind(i,1) <= tot)) THEN
            IF (my_pool_id == INT(index_buf(INT(buff_ind(i, 1)), 6))-1) THEN
              IF (ANY(loc == INT(buff_ind(i, 1)))) THEN
                !
                DO l = 1, INT(buff_ind(i, 2))
                  Eigenvec_alloc_pool(FINDLOC(loc, INT(buff_ind(i, 1))),buff_ind(i, l + 2)) = &  
                                                                      buff(i, l)
                ENDDO
              ELSE
                t = t+1
                DO l = 1, INT(buff_ind(i,2))
                  ! 
                  Eigenvec_alloc_pool(t, INT(buff_ind(i, l + 2))) = buff(i, l)
                  ! 
                ENDDO
                index_buf_pool(t, 1:4) = index_buf(INT(buff_ind(i, 1)), 1:4)
                index_buf_pool(t, 5) = t
                loc(t) = INT(buff_ind(i, 1))
                index_buf_pool(t, 6) = index_buf(buff_ind(i, 1), 6)
                !
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        CALL mp_barrier(inter_pool_comm)
        buff_ind = 0
        buff = 0.d0
        l = 0
      ENDDO
      CALL mp_barrier(inter_pool_comm)
      DEALLOCATE(loc, STAT = ierr)
      IF (ierr /= 0) CALL errore('diag_quasi', 'Error deallocating loc', 1) 
      !
      tot_pool = t
      !   
      CALL mp_sum(t, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      IF (iverbosity == 5) WRITE(stdout, '(/5x,a,2I10)') 'Finished building eigenvectors & 
                                with tot and tot_pool', tot, t
      !  
      Eigenval(1:tot) = ev(:)
      !  
      CALL mp_barrier(inter_pool_comm)
      CALL elpa_deallocate(elp)
      CALL elpa_uninit()
      !  
      IF (my_pool_id == ionode_id) THEN
        !
        PRINT '(/5x,a)','| deallocated elpa.'
        !
      ENDIF
      CALL blacs_gridexit(my_blacs_ctxt)
      DEALLOCATE(al)
      DEALLOCATE(zl)
      DEALLOCATE(ev)
      DEALLOCATE(buff_ind)
      DEALLOCATE(buff)
#endif
    ENDIF
    DEALLOCATE(Eigenvec, STAT = ierr)
    IF (ierr /= 0) CALL errore('diag_quasi', 'Error deallocating Eigenvec', 1)
    !
    ! Finished diagonalization of the many-body Hamiltonian
    !
    CALL cpu_time(finish)
    !
    IF (iverbosity == 5) WRITE(stdout,'(/5x,a,E22.14,a)') 'Time taken for QD bin calculation: ', &
                         (finish-start), 'S'
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE diag_quasi_hamiltonian
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE write_wf_quasi(mesh, typ)
    !---------------------------------------------------------------------------
    !! This subroutine is a utility for  writing QD unitatry matrix/Eigenvector.
    !! It should be noted that for dense mesh grids, Eigenvec is large and causes
    !! memory error. Due to that we only save top 30 eigenvectors with largest
    !! magnitude. Writes a file mesh_INT(mesh).dat for typ=1, where the user
    !! needs to provide a kpoints file for interpolation. The k and k+q are
    !! interpolated to the kpoints value. For typ=2, writes the Hamiltonian, for
    !! typ=3 writes the phonon and direct contribution c_ph_v and c_dir_v for the
    !! meshgrid, for typ=4, only eigenvalues obtained for QD states are written,
    !! which allows to compare eigenvalue  renormalization after QD
    !! diagonalization.
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iuindabs,iunkf
    USE io_global,     ONLY : stdout
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nstemp, fermi_energy, fsthick,           &
                              nkf1, nkf2, nkf3, filkf, mp_mesh_k,      &
                              nq_init
    USE global_var,    ONLY : etf, ibndmin, nkf,nqtotf, wf, wqf,       &
                              nbndfst, nkqtotf, epf17, nkqtotf, xkf,   &
                              Eigenval, Eigenvec_alloc_pool, tot,      &
                              tot_pool, H_ind1, H_ind2, H_mat, r_tot,  &
                              c_ph_v, c_dir_v, index_buf_pool,         &
                              size_tot, Eigenvec_alloc_write,          &
                              size_tot, xkf_write, nkqf, gtemp, xqf
    USE mp_global,     ONLY : my_pool_id, npool, inter_pool_comm,      &
                              ionode_id
    USE cell_base,     ONLY : omega
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two,    &
                              zero, pi, ci, eps6, czero, eps20, eps4
    USE bzgrid,        ONLY : kpmq_map
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: mesh
    !! Current meshgrid    
    INTEGER, INTENT(in) :: typ
    !! Type of information to be printed  
    !
    !Local variables
    ! 
    INTEGER :: i, j, k, t, nktotf, numc, numv, k_init, ikk, ierr,      &
               ikq, ik, cbnd, vbnd, check, indexq, nkfinal,            &
               xk, ipool, num, iqq, xk_final, k_final,                 &
               q_final, xq_final, imode, iq, ab, n, tt, l,             &
               check2(tot), X(2), passk, passkq, passq
    !! Integer variables used for counter and storing indices
    INTEGER, ALLOCATABLE :: x_final(:)
    !! Final k point in direct co-ordinate
    CHARACTER(LEN = 256) :: nameF
    !! File name
    CHARACTER(LEN = 20) :: tp
    !! Temperature
    REAL(KIND = DP) :: ekkcb
    !! Eigen energy on the fine grid relative to the Fermi level for conduction
    !bands
    REAL(KIND = DP) :: ekqcb
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for
    !conduction bands
    REAL(KIND = DP) :: ekkvb
    !! Eigen energy on the fine grid relative to the Fermi level for valence
    !bands
    REAL(KIND = DP) :: ekqvb
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for
    !valence bands
    REAL(KIND = DP) :: ekmk
    !! Eigen energy on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: ekmq,nq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: xt(3),ef0,zeropool(npool),nqv(nmodes)
    !! k point buffer array, Fermi_energy, array of zeros, phonon occupation
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Bose-Einstein distribution 
    REAL(KIND = DP),ALLOCATABLE :: xread(:,:)
    !! Array that stores kpoints provided by the user
    COMPLEX(KIND = DP), ALLOCATABLE:: Eigenvec_alloc_pool2(:,:)
    !! Buffer Eigenvector arrray
    !
    nktotf = nkf1 * nkf2 * nkf3
    n = INT(nq_init)
    ef0 = fermi_energy
    xt = 0.d0
    IF (typ == 1) THEN
      Eigenvec_alloc_write = 0.d0
      xkf_write = 0.d0
      OPEN(UNIT = iunkf, FILE = 'kpoints.txt')
      !    
      READ(iunkf,*) num
      ALLOCATE(xread(3, num), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_wf_quasi', 'Error allocating xread in wf_quasi', 1)
      ALLOCATE(x_final(num), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_wf_quasi', 'Error allocating xread in wf_quasi', 1)
      !
      DO i = 1, num
        !
        READ(iunkf,*) xread(:,i)
        !
      ENDDO
      CLOSE(iunkf)
      x_final(:) = 0
      DO i = 1, num
        DO ik = 1, nkf
          ikk = 2 * ik-1
          xt = 0.d0
          k_init = -1
          !  
          IF ((ABS(xkf(1,ikk)-xread(1,i))  <= eps4) .AND. &
              (ABS(xkf(2,ikk)-xread(2,i))  <= eps4) .AND. &
              (ABS(xkf(3,ikk)-xread(3,i))  <= eps4)) THEN
            CALL kpmq_map(xkf(:, ikk), xt, 1, k_init)
            !
            x_final(i) = k_init
            !
          ENDIF
        ENDDO
      ENDDO
      !  
      CALL mp_barrier(inter_pool_comm)
      CALL mp_sum(x_final, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !  
      k_init    = -1
      k_final   = -1
      q_final   = -1
      tt        =  1
      X(1)      = -1
      X(2)      =  1
      check2(:) = 1
      !  
      DO iqq = 1, nqtotf + 1
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          CALL kpmq_map(xkf(:,ikk), xt, 1, k_init)
          !
          IF (iqq <= nqtotf) THEN
            !
            CALL kpmq_map(xkf(:, ikk), xqf(:, iqq), 1, k_final)
            CALL kpmq_map(xkf(:, ikk), xt, 1, k_init)
            CALL kpmq_map(xqf(:, iqq), xt, 1, q_final)
            !
          ELSE
            ! 
            CALL kpmq_map(xkf(:, ikk), xqf(:, nqtotf), 1, k_final)
            !
          ENDIF
          !
          DO cbnd = 1, nbndfst
            !  the energy of the electron at k (relative to Ef)
            ekkcb = etf(ibndmin - 1 + cbnd, ikk) - ef0
            !
            IF (ekkcb <= 0) CYCLE
            DO vbnd = 1, nbndfst
              !
              ekkvb = etf(ibndmin - 1 + vbnd, ikk) - ef0
              IF (ekkvb >= 0) CYCLE
              indexq = 1
              DO imode = 1, nmodes
                DO ab = 1, 2
                  nqv(imode) = wgauss(-wf(imode, iqq) / gtemp(1), -99)
                  nqv(imode) = nqv(imode) / (one - two * nqv(imode))
                  IF (n == -1) THEN
                    nq = nqv(imode)
                  ELSEIF (n == -2) THEN
                    nq = FLOOR(nqv(imode))
                  ENDIF
                  check = -1
                  IF ((iqq == nqtotf + 1) .AND. (imode == nmodes)) THEN
                    DO i = 1, size_tot
                      !  
                      numc = (cbnd - 1) * nktotf + k_init
                      numv = (vbnd - 1) * nktotf + k_init
                      IF ((numc == INT(index_buf_pool(i, 1))) .AND. &
                          (numv == INT(index_buf_pool(i, 2))) .AND. &
                          (INT(index_buf_pool(i, 3)) < 0)) THEN
                        indexq = INT(index_buf_pool(i, 5))
                        check = 1
                      ENDIF
                    ENDDO
                  ELSE
                    check = INT(check_final(cbnd, k_final, k_init, vbnd, iqq, imode, tot_pool, &
                                index_buf_pool(1:tot_pool, :), (nq - (1 * X(ab))), 2))
                    ! 
                    IF (check /= -1) THEN
                      indexq = INT(index_buf_pool(check, 5))
                    ENDIF
                  ENDIF
                  IF ((check > 0)) THEN
                    DO l = 1, tot
                      IF ((check2(l) <= 30) .AND. (ABS(Eigenvec_alloc_pool(indexq, l)) > eps6)) THEN
                        Eigenvec_alloc_write(my_pool_id + 1, check2(l), l) = Eigenvec_alloc_pool(indexq, l)
                        xk = 0
                        xk_final = 0
                        xq_final = 0
                        !
                        DO j = 1, num
                          IF (ABS(x_final(j) - k_init) == 0) THEN
                            xk = j
                          ENDIF
                          IF (ABS(x_final(j) - k_final) == 0) THEN
                            xk_final = j
                          ENDIF    
                          IF (ABS(x_final(j) - q_final) == 0) THEN
                            xq_final = j
                          ENDIF
                        ENDDO
                        !
                        IF (xk == 0) THEN
                          xk = k_init + nktotf
                        ENDIF
                        !
                        IF (xk_final == 0) THEN
                          xk_final = k_final + nktotf
                        ENDIF
                        !
                        IF (xq_final == 0) THEN
                          xq_final = q_final + nktotf
                        ENDIF
                        !
                        xkf_write(my_pool_id + 1, check2(l), l, 1) = xk
                        xkf_write(my_pool_id + 1, check2(l), l, 2) = vbnd
                        xkf_write(my_pool_id + 1, check2(l), l, 3) = cbnd
                        xkf_write(my_pool_id + 1, check2(l), l, 4) = xk_final
                        xkf_write(my_pool_id + 1, check2(l), l, 5) = xq_final
                        xkf_write(my_pool_id + 1, check2(l), l, 6) = imode
                        check2(l) = check2(l) + 1
                        !
                      ENDIF
                    ENDDO
                  ENDIF
                ENDDO!ab
              ENDDO!imode
            ENDDO!vbnd
          ENDDO !cbnd
        ENDDO !ik
        IF (MOD(iqq, 100) == 0) THEN
          WRITE(stdout,'(/5x,a,I10)') 'progression write iq:', iqq
        ENDIF
      ENDDO!iq
      !  
      CALL mp_barrier(inter_pool_comm)
      CALL mp_sum(Eigenvec_alloc_write, inter_pool_comm)
      CALL mp_sum(xkf_write, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !  
      WRITE(stdout,'(a)') 'Done distributing'
      WRITE(tp,"(I10)") mesh
      nameF = './quasi_bands/mesh_' // trim(adjustl(tp)) // '.dat'
      OPEN(UNIT = iuindabs, FILE = nameF)
      WRITE(iuindabs,'(a,I10)') '# Eigenvector for mesh: ',mesh
      WRITE(iuindabs,'(a)') '# N  k  vband  cband  k+q  q  mode  Eigenval  Re(Eigenvec)  Im(Eigenvec)'
      !
      DO ipool = 1, npool
        DO i = 1, 30
          DO t = 1, tot
            IF (ABS(Eigenvec_alloc_write(ipool, i, t)) > eps20) THEN
              !   
              WRITE(iuindabs,'(I10,9E22.14)') t, xkf_write(ipool, i, t, 1),  &
                      xkf_write(ipool, i, t, 2), xkf_write(ipool, i, t, 3),  &
                      xkf_write(ipool, i, t, 4), xkf_write(ipool, i, t, 5),  &
                      xkf_write(ipool, i, t, 6), REAL(Eigenval(t)) * ryd2ev, &
                      Eigenvec_alloc_write(ipool, i, t)
              !   
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      WRITE(stdout,'(a)') 'finished writing'
      CLOSE(iuindabs)
      DEALLOCATE(x_final, STAT = ierr)
      IF (ierr /=0 ) CALL errore('write_wf_quasi', 'Error deallocating xread in wf_quasi', 1)
      DEALLOCATE(xread, STAT = ierr)
      IF (ierr /=0 ) CALL errore('write_wf_quasi', 'Error deallocating xread in wf_quasi', 1)
      !
    ELSEIF (typ == 2) THEN
      WRITE(tp,"(I10)") mesh
      nameF = './quasi_bands/Hamil_' // trim(adjustl(tp)) // '.dat'
      !
      DO j = 1,npool
        IF (my_pool_id == (j - 1)) THEN
          OPEN(UNIT = iuindabs, FILE = nameF)
          DO k = 1, r_tot
            WRITE(iuindabs,'(2I10,2E22.14)') H_ind1(k), H_ind2(k), (H_mat(k)) * ryd2ev
          ENDDO
          CLOSE(iuindabs)
        ENDIF
      ENDDO
      !
    ELSEIF (typ == 3) THEN
      WRITE(tp,"(I10)") mesh
      nameF = './quasi_bands/c_' // trim(adjustl(tp)) // '.dat'
      !  
      OPEN(UNIT = iuindabs, FILE = nameF)
      DO k = 1, tot
        WRITE(iuindabs,'(4E22.14)') SUM(c_dir_v(1:3, k)), SUM(c_ph_v(1:3, k))
      ENDDO
      CLOSE(iuindabs)
      !
    ELSEIF (typ == 4) THEN
      ALLOCATE(Eigenvec_alloc_pool2(tot, tot))
      Eigenvec_alloc_pool2(:, :) = 0.d0
      WRITE(tp,"(I10)") mesh
      nameF = './quasi_bands/Eigenvals_' // trim(adjustl(tp)) // '.dat'
      !
      OPEN(UNIT = iuindabs, FILE = nameF)
      CALL mp_barrier(inter_pool_comm)
      DO i = 1, npool
        Eigenvec_alloc_pool2(:, :) = 0.d0
        IF (my_pool_id + 1 == i) THEN
          Eigenvec_alloc_pool2(1:tot_pool, :) = Eigenvec_alloc_pool(:, :)
        ENDIF
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(Eigenvec_alloc_pool2, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
        IF (my_pool_id == ionode_id) THEN
          DO k = 1, tot
            DO t = 1, tot
              IF (ABS(Eigenvec_alloc_pool2(k, t)) > eps20) THEN
                WRITE(iuindabs, '(I10,3E22.14)') t, REAL(Eigenval(t)) * ryd2ev, &
                      Eigenvec_alloc_pool2(k, t) !Eigenvec_alloc_write(1,1,k,t)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        CALL mp_barrier(inter_pool_comm)
        !
      ENDDO ! ipool
      CLOSE(iuindabs)
      DEALLOCATE(Eigenvec_alloc_pool2)
    ENDIF
    !---------------------------------------------------------------------------
    END SUBROUTINE write_wf_quasi
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE calc_qdirect(mesh, nomega, iq, n, sum_E, tot_calc, last_calc)
    !---------------------------------------------------------------------------
    !! This subroutine calculates the sum over k and k+q points for epsilon2 for
    !! QD states.
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nstemp, fsthick, degaussw,                            &
                              eps_acoustic, efermi_read, fermi_energy,              &
                              vme, omegamin, omegamax, omegastep, len_mesh,         &
                              nkf1,nkf2,nkf3,scissor,DW,filkf,mode_res
    USE global_var,    ONLY : etf, ibndmin, nkf, epf17, wkf, nqtotf, wf, wqf,       &
                              sigmar_all, efnew, gtemp, nkqtotf, nkqf,              &
                              omegap, vmef, nbndfst, xkf, xqf, etf_ks,              &
                              index_buf_pool, totf_pool, totcv_pool, tot, tot_pool, &
                              epsilon2_qdirect2, epsilon2_qdirect2_DW, E_mesh,      &
                              epsilon2_indirect, Eigenval, c_ph_v,                  &
                              c_dir_v, totf, totcv, Eigenvec_alloc_pool, r_tot, n_q
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi,   &
                              ci, eps6, czero, eps8
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    USE cell_base,     ONLY : omega,bg, at
    USE bzgrid,        ONLY : kpmq_map
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: mesh
    !! Current mesh grid
    INTEGER, INTENT(in)  :: nomega
    !! number of photon gridpoints
    INTEGER, INTENT(in)  :: iq
    !! q point
    INTEGER, INTENT(in)  :: n
    !! Occupation number approximation
    REAL(KIND = DP), INTENT(inout) :: sum_E
    !! Total energy summation (future use)
    INTEGER, INTENT(inout) :: tot_calc
    !! total many-body states calculated
    LOGICAL, INTENT(in)    :: last_calc
    !! If this is the q-point calculation
    ! 
    ! Local variables
    !
    INTEGER :: i, pool_id, k_init, k_final, k_final2, check, check2, vec
    !! index for do loops, pool id, initial and final k, check for
    !! successful finiding of many-body state and eigenvectors
    INTEGER :: j,indexq,nk,nktotf
    !! index for do loops, many-body state,
    INTEGER :: suc,ab,numc,numv
    !! 1 if a many-body state is found, +-1 for absorption and emission,
    !! conduction and valence state location
    INTEGER :: k, l, ikk, ikq, ik, cbnd, vbnd, mbnd, imode
    !! index for k point,
    INTEGER :: ipol, iw
    !! Polarization direction
    INTEGER :: m
    !! Counter on denominator imaginary broadening values
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    INTEGER, PARAMETER :: neta = 9
    !! Broadening parameter
    REAL(KIND = DP) :: alpha, final_e
    !! +-1 for asorption or emission and final state energy
    REAL(KIND = DP) :: ekkcb
    !! Eigen energy on the fine grid relative to the Fermi level for conduction
    !bands
    REAL(KIND = DP) :: ekqcb
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for
    !conduction bands
    REAL(KIND = DP) :: ekkvb
    !! Eigen energy on the fine grid relative to the Fermi level for valence
    !bands
    REAL(KIND = DP) :: ekqvb
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for
    !valence bands
    REAL(KIND = DP) :: ekmk
    !! Eigen energy on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: ekmq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: nqv(nmodes)
    !! Phonon frequencies and phonon occupations on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy
    REAL(KIND = DP) :: st, nq
    !! Phonon occupation
    REAL(KIND = DP) :: wgkk, wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$, $f_{nk}(T)$
    REAL(KIND = DP) :: xt(3),xxk(3),xxkq(3)
    !! Buffer k-point, k+G(0), k+q 
    REAL(KIND = DP) :: weighta, weighte, weightd, weightq
    !!- delta function for absorption, emission, direct
    REAL(KIND = DP) :: inv_eptemp0
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: pfac
    !! Occupation prefactors
    REAL(KIND = DP) :: pface
    !! Occupation prefactors
    REAL(KIND = DP) :: cfac,r
    !! Absorption prefactor
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    COMPLEX(KIND = DP) :: vkk(3, nbndfst, nbndfst)
    !! Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: vkq(3, nbndfst, nbndfst)
    !! Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: Aa(3), Ae(3), Ba(3), Be(3), Ca(3), Ce(3), Da(3), De(3), Qa(3), &
                          DWca(3), DWva(3)    
    !! For summation over virutal states
    COMPLEX(KIND = DP) :: epf(nbndfst, nbndfst, nmodes)
    !! Electron-phonon coupling
    !
    ef0 = fermi_energy
    !
    inv_degaussw = 1.0 / degaussw
    !
    ! 300 K
    ! Epsilon2 prefactor for velocity matrix elements, including factor of 2 
    ! for spin (weights for k-points are divided by 2 to be normalized to 1)
    ! C = 8*pi^2*e^2 = 8*pi^2*2 approx 157.9136704
    !
    cfac = 16.d0 * pi**2
    check = -1
    nktotf = nkf1 * nkf2 * nkf3
    itemp = 1
    IF (filkf /= ' ') THEN
      nktotf = nkqtotf / 2
    ELSE
      nktotf = nkf1 * nkf2 * nkf3
    ENDIF
    nq = n
    xt = 0.d0
    !Eig=0.d0
    DO imode = 1, nmodes
      nqv(imode) = wgauss(-wf(imode,iq) / gtemp(itemp), -99)
      nqv(imode) = nqv(imode) / (one - two * nqv(imode))
      IF (n == -1) THEN
        nq = nqv(imode)
      ELSEIF (n == -2) THEN
        nq = FLOOR(nqv(imode))
      ELSEIF (n == -3) THEN
        ! 
        nq = n_q(iq,imode)
      ENDIF
      sum_E = sum_E+(nq*wf(imode,iq))
    ENDDO

    DO itemp = 1, nstemp
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        CALL kpmq_map(xkf(:, ikk), xqf(:, iq), 1, k_final)
        CALL kpmq_map(xkf(:, ikk), xt, 1, k_init)
        DO cbnd = 1, nbndfst
          DO vbnd = 1, nbndfst
            ! vmef is in units of Ryd * bohr
            vkk(:, cbnd, vbnd) = vmef(:, ibndmin - 1 + cbnd, ibndmin - 1 + vbnd, ikk) !* (wkf(ikk)/2.0)
            vkq(:, cbnd, vbnd) = vmef(:, ibndmin - 1 + cbnd, ibndmin - 1 + vbnd, ikq) !* (wkf(ikk)/2.0)
          ENDDO
        ENDDO
        !
        DO cbnd = 1, nbndfst
          !  the energy of the electron at k (relative to Ef)
          ekkcb = etf(ibndmin - 1 + cbnd, ikk) - ef0
          ekqcb = etf(ibndmin - 1 + cbnd, ikq) - ef0
          IF (ekkcb <= zero) CYCLE
          !
          IF (ABS(ekkcb) < fsthick) THEN
           !
            DO vbnd = 1, nbndfst
             !
              ekkvb = etf(ibndmin - 1 + vbnd, ikk) - ef0
              ekqvb = etf(ibndmin - 1 + vbnd, ikq) - ef0
              IF (ekkvb > zero) CYCLE

              DO imode = 1, nmodes
                epf(:, :, imode) = epf17(:, :, imode,ik)
                nqv(imode) = wgauss(-wf(imode,iq) / gtemp(itemp), -99)
                nqv(imode) = nqv(imode) / (one - two * nqv(imode))
                IF (n == -1) THEN
                  nq = nqv(imode)
                ELSEIF (n == -2) THEN
                  nq = FLOOR(nqv(imode))
                ELSEIF (n == -3) THEN
                  nq = n_q(iq,imode)
                ENDIF
                check = -1
                DO ab = 1, 2
                  suc = 0
                  IF (ab == 1) THEN
                    alpha = 1.0
                    check = INT(check_final(cbnd, k_final, k_init, vbnd, iq, imode, &
                             tot_pool, index_buf_pool(1:tot_pool, :), nq + 1, 2))
                    IF (check /= -1) THEN
                      suc = 1
                      indexq = INT(index_buf_pool(check, 5))
                      tot_calc = tot_calc + 1
                    ENDIF
                    !
                  ELSEIF (ab == 2) THEN
                    check = INT(check_final(cbnd, k_final, k_init, vbnd, iq, imode, &
                            tot_pool, index_buf_pool(1:tot_pool, :), nq - 1, 2))
                    alpha = -1.0
                    IF (check /= -1) THEN
                      suc = 1
                      indexq = INT(index_buf_pool(check, 5))
                      tot_calc = tot_calc + 1
                    ENDIF
                    !
                  ENDIF
                  IF (tot == 0) THEN
                    suc = 0
                  ENDIF
                  !
                  IF (wf(imode,iq) < eps_acoustic) CYCLE
                  !
                  Aa = czero
                  Ba = czero
                  Ca = czero
                  Da = czero
                  DWca= czero
                  DWva= czero
                  Ae = czero
                  Be = czero
                  Ce = czero
                  De = czero
                  !
                  DO mbnd = 1, nbndfst
                    !
                    ! The energy of the electron at k (relative to Ef)
                    ekmk = etf(ibndmin - 1 + mbnd, ikk) - ef0
                    ! The energy of the electron at k+q (relative to Ef)
                    ekmq = etf(ibndmin - 1 + mbnd, ikq) - ef0
                    final_e = (E_mesh(mesh) + E_mesh(mesh + 1)) / 2.0
                    !c-c' transition
                    IF ((ekmq > zero) .AND. (ekmk > zero)) THEN
                      !  
                      IF ((mesh > 1) .AND. (mesh < len_mesh-1)) THEN
                        !   
                        IF (((ekmk - ekkvb) < E_mesh(mesh))  .OR. & 
                            ((ekmk - ekkvb) > E_mesh(mesh + 1))) THEN
                          ! 
                          Aa(:) = Aa(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd)         &
                                           /(final_e - (ekmk - ekkvb))
                          DWca(:) = DWca(:) + epf(cbnd, mbnd, imode) * epf(mbnd, cbnd, imode) &
                                           /(final_e - (ekmk - ekkvb))
                          !
                        ENDIF
                        !   
                      ELSEIF (mesh >= len_mesh - 1) THEN
                        !
                        IF ((ekmk - ekkvb) < E_mesh(mesh)) THEN
                          !
                          Aa(:) =Aa(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd)          &
                                 /(final_e - (ekmk - ekkvb))
                          DWca(:) = DWca(:) + epf(cbnd, mbnd, imode) * epf(mbnd, cbnd, imode) &
                                 /(final_e - (ekmk - ekkvb))
                          !
                        ENDIF
                        !
                      ELSEIF (mesh <= 1) THEN
                        !    
                        IF ((ekmk-ekkvb) > E_mesh(mesh + 1)) THEN
                          !
                          Aa(:) = Aa(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd)         &
                                   /(final_e - (ekmk - ekkvb))
                          DWca(:) = DWca(:) + epf(cbnd, mbnd, imode) * epf(mbnd, cbnd, imode) &
                                   /(final_e - (ekmk - ekkvb))
                          !
                        ENDIF
                        !
                      ENDIF
                      Ca(:) = Ca(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode)             &
                                       / (ekkvb - ekmq - wf(imode, iq) * alpha)
                    !v-v' transition
                    ELSEIF ((ekmq < zero) .AND. (ekmk < zero)) THEN
                      !  
                      IF ((mesh > 1) .AND. (mesh < len_mesh - 1)) THEN
                        !
                        IF (((ekqcb - ekmq) < E_mesh(mesh)) .OR. &
                            ((ekqcb - ekmq) > E_mesh(mesh + 1))) THEN
                          !
                          Ba(:) = Ba(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode)         &
                                   /(final_e - (ekqcb - ekmq))
                          DWva(:) = DWva(:) + epf(vbnd, mbnd, imode) * epf(mbnd, vbnd, imode) &
                                   /(final_e - (ekqcb - ekmq))
                          !
                        ENDIF
                      ELSEIF ((mesh >= len_mesh - 1)) THEN
                        !
                        IF ((ekqcb - ekmq) < E_mesh(mesh)) THEN
                          !
                          Ba(:) = Ba(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode)         &
                                  /(final_e - (ekqcb - ekmq))
                          DWva(:) = DWva(:) + epf(vbnd, mbnd, imode) * epf(mbnd, vbnd, imode) &
                                  /(final_e - (ekqcb - ekmq))
                          !  
                        ENDIF
                      ELSEIF (mesh <= 1) THEN
                        !    
                        IF ((ekqcb - ekmq) > E_mesh(mesh + 1)) THEN
                          ! 
                          Ba(:) = Ba(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode)         &
                                  /(final_e - (ekqcb - ekmq))
                          DWva(:) = DWva(:) + epf(vbnd, mbnd, imode) * epf(mbnd, vbnd, imode) &
                                  /(final_e - (ekqcb - ekmq))
                          !  
                        ENDIF
                      ENDIF
                      Da(:) = Da(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd)             &
                                    /(ekmk - ekqcb - wf(imode, iq) * alpha)
                      !
                    ENDIF
                  ENDDO
                  !
                  IF (n < 2) THEN
                    pfac  =  nq + ((1.0 / 2.0) * alpha + (1.0 / 2.0))
                  ELSE
                    !  Experimental  
                    pfac = (nq + 1.0/2 + alpha*(1.0/2))*(EXP(-nq*wf(imode,iq)))
                  ENDIF
                  IF ((suc == 1) .AND. (.NOT. last_calc)) THEN
                    !
                    Ae(:) = Aa(:)
                    Be(:) = Ba(:)
                    Ce(:) = Ca(:)
                    De(:) = Da(:)
                    !
                    DO ipol = 1, 3
                      !  
                      epsilon2_qdirect2(ipol, 1:tot,1) = epsilon2_qdirect2(ipol, 1:tot,1) +         &
                                                   DSQRT(wkf(ikk) / 2.0) * DSQRT(wqf(iq)) *         &
                              (DSQRT(pfac)  * (Ae(ipol) + Be(ipol) + Ce(ipol) + De(ipol)) *         &
                              CONJG(Eigenvec_alloc_pool(indexq, 1:tot)) / DSQRT(2.0 * wf(imode, iq)))
                      IF (imode <= mode_res) THEN
                        !
                        epsilon2_qdirect2(ipol, 1:tot, imode + 1) =  &
                                          epsilon2_qdirect2(ipol, 1:tot, imode + 1)     +           &
                                          DSQRT(wkf(ikk) / 2.0) * DSQRT(wqf(iq))        *           &
                             (DSQRT(pfac) * (Ae(ipol) + Be(ipol) + Ce(ipol) + De(ipol)) *           &
                             CONJG(Eigenvec_alloc_pool(indexq, 1:tot)) / DSQRT(2.0 * wf(imode, iq)))
                      ENDIF
                      !  
                      c_ph_v(ipol, 1:tot) = c_ph_v(ipol, 1:tot) + DSQRT(wkf(ikk) / 2.0)  *          &
                                                         DSQRT(wqf(iq)) * (DSQRT(pfac)   *          &
                                             (Ae(ipol) + Be(ipol) + Ce(ipol) + De(ipol)) *          &
                              CONJG(Eigenvec_alloc_pool(indexq, 1:tot)) / DSQRT(2.0 * wf(imode,iq)))
                      !
                      epsilon2_qdirect2_DW(ipol, 1:tot) = epsilon2_qdirect2_DW(ipol, 1:tot) +       &
                                                               (wkf(ikk) / 2.0) * (wqf(iq)) *       &
                                       ((pfac)  * CONJG(Eigenvec_alloc_pool(indexq, 1:tot)) *       &
                                                                    (Dwca(ipol)+DWva(ipol)) *       &
                                Eigenvec_alloc_pool(indexq, 1:tot) / (2.0 * wf(imode, iq)))
                    ENDDO ! ipol 
                  ENDIF ! last_calc
                ENDDO ! absorption and emission
              ENDDO  ! imode
              suc = 0
              ! Direct part is calculated for the last iq!!!!!!!!!
              IF ((iq == nqtotf) .AND. (last_calc)) THEN
                Qa = czero
                DO i = 1, tot_pool!cv_pool
                  numc = (cbnd - 1) * nktotf + k_init
                  numv = (vbnd - 1) * nktotf + k_init
                  IF ((numc == INT(index_buf_pool(i, 1))) .AND.&
                     (numv == INT(index_buf_pool(i, 2))) .AND.&
                     (INT(index_buf_pool(i, 3)) < 0)) THEN
                    indexq = INT(index_buf_pool(i, 5))
                    suc = 1
                    tot_calc = tot_calc + 1
                  ENDIF
                ENDDO
                !
                IF (suc == 1) THEN
                  IF (DW == 0) THEN
                    !
                    DO ipol = 1, 3
                      !
                      Qa(ipol) = vkk(ipol, cbnd, vbnd)!*CONJG(Eigenvec_alloc_pool(indexq,vec))
                      epsilon2_qdirect2(ipol, 1:tot, 1) = epsilon2_qdirect2(ipol, 1:tot, 1) +   &
                                                         (Qa(ipol)) * DSQRT(wkf(ikk) / 2.0) *   & 
                                                   CONJG(Eigenvec_alloc_pool(indexq, 1:tot))
                      !   
                      c_dir_v(ipol, 1:tot) = c_dir_v(ipol, 1:tot) +  Qa(ipol) *                 &
                           DSQRT(wkf(ikk) / 2.0) * CONJG(Eigenvec_alloc_pool(indexq, 1:tot))
                      !
                    ENDDO ! ipol
                    !
                  ELSE  
                    ! DW term
                    DO vec = 1, tot
                      Qa(:) = vkk(:, cbnd, vbnd) * CONJG(Eigenvec_alloc_pool(indexq, vec))
                      epsilon2_qdirect2(:, vec, 1) = epsilon2_qdirect2(:, vec, 1) +             &
                                                    (Qa(:)) * DSQRT(wkf(ikk)/2.0)
                      !  
                      DO l = 1, tot
                        IF (l /= vec) THEN
                          IF (ABS(REAL(Eigenval(vec) - Eigenval(l))) > 0.0001) THEN
                            epsilon2_qdirect2(:, vec, 1) = epsilon2_qdirect2(:, vec, 1) +         &
                                    (vkk(:, cbnd, vbnd)*Eigenvec_alloc_pool(indexq, l)  *         &
                                                           epsilon2_qdirect2_DW(:, vec) /         &
                               REAL(Eigenval(vec)-Eigenval(l)+(0*0.001))) * DSQRT(wkf(ikk)/2.0)
                          ENDIF
                        ENDIF
                      ENDDO ! Eigenval
                    ENDDO ! Eigenval
                  ENDIF ! DW
                ENDIF ! If state found
              ENDIF ! Last calc
            ENDDO!vbnd
          ENDIF!fsthick
        ENDDO!cbnd
      ENDDO!ik
    ENDDO! itemp
    !---------------------------------------------------------------------------
    END SUBROUTINE calc_qdirect
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE calc_qdirect_DW(mesh, nomega, iq, n, sum_E, calc, tot_calc_DW)
    !---------------------------------------------------------------------------
    !! This subroutine calculates the epsilon2 for  states not in QD bin or
    !! inside QD bin but Eigenvector (1,m)=1.0 which means these sates do not
    !! entangle. In general, it's 0 for very large grids. This subroutine also
    !! calculates the second order contribution to epsilon2 for degenerate states
    !! which is controlled using the keyword DW. Although, this contibution is
    !! below 1% in most of the tested cases.
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nstemp, fsthick, degaussw,                                 &
                              eps_acoustic, efermi_read, fermi_energy,                   &
                              vme, omegamin, omegamax, omegastep, len_mesh,              &
                              nkf1, nkf2, nkf3, scissor,filkf, mode_res
    USE global_var,    ONLY : etf, ibndmin, nkf, epf17, wkf, nqtotf, wf, wqf,            &
                              sigmar_all, efnew, gtemp, nkqtotf, nkqf, etf_ks,           &
                              omegap, vmef, nbndfst, xkf, xqf, index_buf_pool,           &
                              totf_pool, tot_pool, tot, totcv_pool, n_q, E_mesh,         &
                              epsilon2_qdirect, Eigenval, c_ph, c_dir, totf, totcv
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6,  &
                              czero, eps8
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    USE cell_base,     ONLY : omega, bg, at
    USE bzgrid,        ONLY : kpmq_map
    !
    IMPLICIT NONE 
    !
    INTEGER, INTENT(in)  :: mesh
    !! Index of the current meshgrid
    INTEGER, INTENT(in)  :: nomega
    !! number of points for photon grid
    INTEGER, INTENT(in)  :: iq
    !! q-point index
    INTEGER, INTENT(in)  :: n
    !! phonon occupation approximation 
    REAL(KIND = DP), INTENT(inout) :: sum_E
    !! Total energy summation
    INTEGER, INTENT(in)  :: calc
    !! Type of calculation
    INTEGER, INTENT(inout) :: tot_calc_DW
    !! total number of many-body states that are non-degenrate
    !
    ! Local variables
    !
    INTEGER :: i, eig, vec, pool_id, k_init, k_final, check
    !! counter on many-body states, eigenvalue counter, eigenvector counter,
    !! current id of pool,initial and final k, positive if many-body state found
    INTEGER :: j,indexq,nk,nktotf
    !! Counter, index of many-body state and total number of k points
    INTEGER :: suc,ab,numc,numv
    !! 1 if many-body state found, counter for absorption or emission, conduction and
    !! valence state location
    INTEGER :: k, l, ikk, ikq, ik, cbnd, vbnd, mbnd, iw, imode, DW
    !! counter on k, k+q points, conduction, valence, and intermediate points,
    !! photon energym and second-order controbution
    REAL(KIND = DP) :: alpha
    !! +-1 for absorption or emission
    INTEGER :: ipol
    !! Polarization direction
    INTEGER :: m
    !! Counter on denominator imaginary broadening values
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    INTEGER, PARAMETER :: neta = 9
    !! Broadening parameter
    REAL(KIND = DP) :: final_e
    !! Eigen energy on the fine grid relative to the Fermi level for conduction
    REAL(KIND = DP) :: ekkcb
    !! Eigen energy on the fine grid relative to the Fermi level for conduction
    !bands
    REAL(KIND = DP) :: ekqcb
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for
    !conduction bands
    REAL(KIND = DP) :: ekkvb
    !! Eigen energy on the fine grid relative to the Fermi level for valence
    !bands
    REAL(KIND = DP) :: ekqvb
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for
    !valence bands
    REAL(KIND = DP) :: ekmk
    !! Eigen energy on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: ekmq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: nqv(nmodes)
    !! Phonon frequencies and phonon occupations on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy
    REAL(KIND = DP) :: st,nq
    !! Phonon occupation
    REAL(KIND = DP) :: wgkk, wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$, $f_{nk}(T)$
    REAL(KIND = DP) :: weighta, weighte, weightd, weightq
    !!- delta function for absorption, emission, direct
    REAL(KIND = DP) :: inv_eptemp0
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: pfac
    !! Occupation prefactors
    REAL(KIND = DP) :: pface
    !! Occupation prefactors
    REAL(KIND = DP) :: cfac
    !! Absorption prefactor
    REAL(KIND = DP) :: xt(3),xxk(3),xxkq(3)
    !! k point buffer, k+G(0), K+q
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    COMPLEX(KIND = DP) :: vkk(3, nbndfst, nbndfst)
    !!- Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: vkq(3, nbndfst, nbndfst)
    !!- Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: Aa(3), Ae(3), Ba(3), Be(3), Ca(3), Ce(3), Da(3), De(3), Qa(3)
    !! The four terms in conduction and valence band notation for virtual state
    !! summation
    COMPLEX(KIND = DP) :: epf(nbndfst, nbndfst, nmodes)
    !! electron-phonon matrix
    !    
    itemp = 1
    ! QDPT works for only one temperature
    ef0=fermi_energy
    !
    inv_degaussw = 1.0 / degaussw
    !
    ! 300 K
    ! Epsilon2 prefactor for velocity matrix elements, including factor of 2 for spin 
    ! (weights for k-points are divided by 2 to be normalized to 1)
    ! C = 8*pi^2*e^2 = 8*pi^2*2 approx 157.9136704
    !
    cfac = 16.d0 * pi**2
    check=-1
    nktotf=nkf1 * nkf2 * nkf3
    !
    IF (filkf /=' ') THEN
      nktotf = nkqtotf / 2
    ELSE
      nktotf = nkf1 * nkf2 * nkf3
    ENDIF
    nq = n
    xt = 0.d0
    DO imode = 1, nmodes
      nqv(imode) = wgauss(-wf(imode,iq) / gtemp(itemp), -99)
      nqv(imode) = nqv(imode) / (one - two * nqv(imode))
      IF (n == -1) THEN
        nq = nqv(imode)
      ELSEIF (n == -2) THEN
        nq = FLOOR(nqv(imode))
      ELSEIF (n == -3) THEN
        nq = n_q(iq, imode)
      ENDIF
      sum_E = sum_E + (nq * wf(imode, iq))
    ENDDO
    DO DW = 1, 1
      DO itemp =1, nstemp
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          xxk(:) = xkf(:, ikk)
          xxkq(:) = xqf(:, iq)
          !  
          k_init = -1
          k_final = -1
          CALL kpmq_map(xxk, xxkq, 1, k_final)
          CALL kpmq_map(xxk, xt, 1, k_init)
          DO cbnd = 1, nbndfst
            DO vbnd = 1, nbndfst
              ! vmef is in units of Ryd * bohr
              vkk(:, cbnd, vbnd) = vmef(:, ibndmin - 1 + cbnd, ibndmin - 1 + vbnd, ikk) !* (wkf(ikk)/2.0)
              vkq(:, cbnd, vbnd) = vmef(:, ibndmin - 1 + cbnd, ibndmin - 1 + vbnd, ikq) !* (wkf(ikk)/2.0)
            ENDDO
          ENDDO
          !
          DO cbnd = 1, nbndfst
            !  the energy of the electron at k (relative to Ef)
            ekkcb = etf(ibndmin - 1 + cbnd, ikk) - ef0
            ekqcb = etf(ibndmin - 1 + cbnd, ikq) - ef0
            IF ((ekkcb <= zero) .AND. (ekqcb <= zero)) CYCLE
            !
            IF ((ABS(ekkcb) < fsthick) .AND. (ABS(ekqcb) < fsthick)) THEN
            !
              DO vbnd = 1, nbndfst
                ekkvb = etf(ibndmin - 1 + vbnd, ikk) - ef0
                ekqvb = etf(ibndmin - 1 + vbnd, ikq) - ef0
                !
                IF ((ekkvb > zero) .AND. (ekqvb > zero)) CYCLE
                !
                !wgkq = wgauss(-ekqcb * inv_eptemp0, -99)
                DO imode = 1, nmodes
                  epf(:, :, imode) = epf17(:, :, imode,ik)
                  nqv(imode) = wgauss(-wf(imode,iq) / gtemp(itemp), -99)
                  nqv(imode) = nqv(imode) / (one - two * nqv(imode))
                  IF (n == -1) THEN
                    nq = nqv(imode)
                  ELSEIF (n == -2) THEN
                    nq = FLOOR(nqv(imode))
                  ELSEIF (n == -3) THEN
                    nq = n_q(iq,imode)
                  ENDIF
                  IF (ekqcb > ekkvb + wf(nmodes, iq) + omegamax + 6.0 * degaussw) CYCLE
                  DO ab =1, 2
                    suc = 0
                    !
                    IF (ab == 1) THEN
                      alpha = 1.0
                      IF (calc /= 1) THEN
                        check = INT(check_final(cbnd, k_final, k_init, vbnd, iq, imode, &
                                   tot_pool, index_buf_pool(1:tot_pool, 1:5), nq + 1, 2))
                        IF (check /= -1) THEN
                          suc = 1
                          tot_calc_DW = tot_calc_DW + 1
                        ENDIF
                      ENDIF   
                      !  
                    ELSEIF (ab == 2) THEN
                      alpha = -1.0
                      IF (calc /= 1) THEN 
                        check = INT(check_final(cbnd, k_final, k_init, vbnd, iq, imode, &
                                   tot_pool, index_buf_pool(1:tot_pool, 1:5), nq - 1, 2))
                        IF (check /= -1) THEN
                          suc = 1
                          tot_calc_DW = tot_calc_DW + 1
                        ENDIF                        
                      ENDIF
                    ENDIF
                    !
                    IF (((ekqcb - ekkvb + alpha * wf(imode, iq)) < E_mesh(mesh)) .OR. &
                       ((ekqcb - ekkvb + alpha * wf(imode, iq)) > E_mesh(mesh + 1))) CYCLE
                    !
                    IF (wf(imode, iq) > eps_acoustic) THEN
                      !
                      Aa = czero
                      Ba = czero
                      Ca = czero
                      Da = czero
                      Ae = czero
                      Be = czero
                      Ce = czero
                      De = czero
                      !
                      DO mbnd = 1, nbndfst
                        !
                        ! The energy of the electron at k (relative to Ef)
                        ekmk = etf(ibndmin - 1 + mbnd, ikk) - ef0
                        ! The energy of the electron at k+q (relative to Ef)
                        ekmq = etf(ibndmin - 1 + mbnd, ikq) - ef0
                        !    
                        IF (calc == 1) THEN
                          final_e = ekqcb - ekkvb + alpha * wf(imode, iq)
                        ELSEIF (calc == 2) THEN
                          final_e = (E_mesh(mesh) + E_mesh(mesh + 1)) / 2.0
                        ENDIF
                        !  
                        IF ((ekmq > zero) .AND. (ekmk > zero)) THEN
                          !c->c'  
                          IF ((mesh > 2) .AND. (mesh < len_mesh-1)) THEN
                            !
                            IF (((ekmk - ekkvb) < E_mesh(mesh)) .OR. &
                                ((ekmk - ekkvb) > E_mesh(mesh + 1))) THEN
                              Aa(:) = Aa(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd) &
                                         /(final_e - (ekmk - ekkvb))
                            ENDIF
                            !
                          ELSEIF ((mesh >= len_mesh - 1)) THEN
                            !
                            IF ((ekmk-ekkvb) < E_mesh(mesh - 1)) THEN
                              Aa(:) = Aa(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd) &
                                         /(final_e - (ekmk - ekkvb))
                            ENDIF
                          ELSEIF (mesh <= 2) THEN
                            !
                            IF ((ekmk-ekkvb) > E_mesh(mesh + 1)) THEN
                              Aa(:) = Aa(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd) &
                                          /(final_e -(ekmk - ekkvb))
                            ENDIF
                            !
                          ENDIF
                          !v->c'
                          Ca(:) = Ca(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode)     &
                                    /(ekkvb - ekmq - wf(imode, iq) * alpha)
                          !
                        !
                        ELSEIF ((ekmq < zero) .AND. (ekmk < zero)) THEN
                          ! v->v'  
                          IF ((mesh > 2) .AND. (mesh < len_mesh - 1 )) THEN
                            !
                            IF (((ekqcb - ekmq) < E_mesh(mesh)) .OR. &
                                ((ekqcb - ekmq) > E_mesh(mesh + 1))) THEN
                              !
                              Ba(:) = Ba(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode) &
                                           /(final_e - (ekqcb - ekmq))
                              !
                            ENDIF
                          ELSEIF ((mesh >= len_mesh - 1)) THEN
                        
                            IF ((ekqcb - ekmq) < E_mesh(mesh)) THEN
                        
                              Ba(:) = Ba(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode) &
                                          /(final_e - (ekqcb - ekmq))
                              !
                            ENDIF
                          ELSEIF (mesh <= 2) THEN
                            !
                            IF ((ekqcb - ekmq) > E_mesh(mesh + 1)) THEN
                              !   
                              Ba(:) = Ba(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode) &
                                          /(final_e - (ekqcb - ekmq))
                              !
                            ENDIF
                            !
                          ENDIF
                           ! v->c'
                           Da(:) = Da(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd)    &
                                   /(ekmk-ekqcb-wf(imode,iq)*alpha)
                           ! 
                        ENDIF
                        ! 
                      ENDDO ! mbnd
                      !
                      IF (n < 2) THEN
                         pfac = nq + (1.0 / 2.0) * alpha + (1.0 / 2.0)
                      !
                      ELSE
                        pfac = (nq + 1.0 / 2.0 + alpha * (1.0 / 2.0)) * (EXP(-nq * wf(imode, iq)))
                      ENDIF
                      IF (suc == 1) THEN
                        CYCLE
                      ELSEIF (suc == 0) THEN
                        ! 
                        Ae(:) = Aa(:)
                        Be(:) = Ba(:)
                        Ce(:) = Ca(:)
                        De(:) = Da(:)
                        DO iw = 1, nomega
                          !
                          weighta = w0gauss((ekqcb - ekkvb - omegap(iw) + & 
                                alpha * wf(imode, iq)) / degaussw, 0) / degaussw
                          weightd = w0gauss((ekkcb - ekkvb - omegap(iw)) / (1.0*degaussw), 0) / (1.0 * degaussw)
                          !
                          DO ipol = 1, 3
                            IF ((tot == 0) .OR. (calc == 2)) THEN
                              epsilon2_qdirect(ipol, iw, 1, itemp) = epsilon2_qdirect(ipol, iw, 1, itemp) + &
                                                                            (wkf(ikk) / 2.0) * wqf(iq)    * &
                                                                         (((cfac / omegap(iw)**2) * pfac  * &
                                                                       weighta * ABS(Ae(ipol) + Be(ipol)  + &
                                                     Ce(ipol) + De(ipol))**2 / (2 * wf(imode,iq) * omega)))
                              IF (calc /= 2) THEN

                                epsilon2_qdirect(ipol, iw, 2, itemp) = epsilon2_qdirect(ipol, iw, 1, itemp) + &
                                                                                 (wkf(ikk) / 2.0) * wqf(iq) * &
                                                                 (((cfac / omegap(iw)**2) * pfac  * weighta * &
                                                                ABS(Ae(ipol) + Be(ipol)+Ce(ipol)+De(ipol))**2 &
                                                                  / (2 * wf(imode,iq) * omega)))
                              ENDIF
                      
                              IF (imode <= mode_res) THEN
                                epsilon2_qdirect(ipol, iw, neta+imode, itemp) =                     &  
                                              epsilon2_qdirect(ipol, iw, neta + imode, itemp) +     &
                                                                   (wkf(ikk) / 2.0) * wqf(iq) *     &
                                    (((cfac / omegap(iw)**2) * pfac  * weighta * ABS(Ae(ipol) +     &
                                  Be(ipol) + Ce(ipol) + De(ipol))**2 / (2 * wf(imode, iq) * omega)))
                              ENDIF
                              !  
                              c_ph(ipol, iw) = c_ph(ipol, iw) + (wkf(ikk) / 2.0) * wqf(iq) *        &
                                  (((cfac / omegap(iw)**2) * pfac  * weighta * ABS(Ae(ipol) +       &
                                  Be(ipol) + Ce(ipol) + De(ipol))**2 / (2 * wf(imode, iq) * omega)))
                              !          
                            ENDIF
                            epsilon2_qdirect(ipol, iw, 4, itemp) = epsilon2_qdirect(ipol, iw, 4, itemp) + &
                                                                             (wkf(ikk) / 2.0) * wqf(iq) * &
                                                ((cfac / omegap(iw)**2 * pfac  * weighta * ABS(Ae(ipol) + &
                                   Be(ipol) + Ce(ipol) + De(ipol))**2 / (2 * wf(imode,iq) * omega)))
                            ! 
                          ENDDO ! ipol
                        ENDDO ! omega  
                      ENDIF ! success
                    ENDIF ! eps_acoustic
                  ENDDO ! absorption and emission
                ENDDO  ! imode
                ! 
                IF (iq == nqtotf) THEN
                  Qa = czero
                  suc = 0
                  !  
                  DO i = 1, tot_pool
                    numc = (cbnd-1)*nktotf+k_init
                    numv = (vbnd-1)*nktotf+k_init
                    !
                    IF ((numc == INT(index_buf_pool(i,1))) .AND. &
                        (numv == INT(index_buf_pool(i,2))) .AND. &
                        (INT(index_buf_pool(i,3)) < 0)) THEN
                      suc = 1
                      tot_calc_DW = tot_calc_DW+1
                    ENDIF
                    ! 
                  ENDDO
                  !   
                  IF ((suc == 0) .AND. ((ekkcb - ekkvb) > E_mesh(mesh)) .AND. &
                      ((ekkcb - ekkvb) < E_mesh(mesh + 1))) THEN
                    !
                    Qa(:) = vkk(:, cbnd, vbnd)
                    DO iw=1, nomega
                      !   
                      weightq= (w0gauss((ekkcb - ekkvb - omegap(iw)) / degaussw, 0) / degaussw) * &
                               (wkf(ikk) / 2.0)
                      !
                      DO ipol = 1, 3
                        IF ((tot == 0) .OR. (calc == 2)) THEN

                          epsilon2_qdirect(ipol, iw, 1, itemp) = epsilon2_qdirect(ipol, iw, 1, itemp) + &
                                  ((cfac / omegap(iw)**2 ) * weightq * ABS(Qa(ipol))**2 / omega**1)
                          c_dir(ipol, iw) = c_dir(ipol,iw) + ((cfac / omegap(iw)**2 ) * weightq *       & 
                                             ABS(Qa(ipol))**2 / omega**1)
                          !
                        ENDIF
                        epsilon2_qdirect(ipol, iw, 5, itemp) = epsilon2_qdirect(ipol,iw, 5, itemp) + &
                                ((cfac / omegap(iw)**2 ) * weightq * ABS(Qa(ipol))**2 / omega**1)
                      ENDDO !
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO!vbnd
            ENDIF
          ENDDO!cbnd
        ENDDO!ik
      ENDDO! itemp
    ENDDO! DW
    !--------------------------------------------------------------------------
    END SUBROUTINE calc_qdirect_DW
    !--------------------------------------------------------------------------
    ! 
    !--------------------------------------------------------------------------
    FUNCTION check_final(cbnd, k_final, k_init, vbnd, iq, imode, totf, &
                indextot, n, test)
    !--------------------------------------------------------------------------
    !! This function checks if the final state generated using k and k+q loops
    !! is inside QD bin/QD matrix. If it is, it returns the location of the state
    !! |vk>|ck+q> in the row of Eigenvec matrix.
    !
    USE io_global,     ONLY : stdout, ionode_id
    USE global_var,         ONLY : nqtotf, nkqtotf
    USE input,        ONLY : nkf1, nkf2, nkf3, filkf
    USE kinds,         ONLY : DP
    USE ep_constants,     ONLY : eps4, eps8
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: cbnd
    !! Conduction band index
    INTEGER, INTENT(in) :: k_final
    !! k state of final MB state
    INTEGER, INTENT(in) :: k_init
    !! k-state of initial MB state
    INTEGER, INTENT(in) :: vbnd
    !! Valence band index
    INTEGER, INTENT(in) :: iq
    !! q-point index
    INTEGER, INTENT(in) :: imode
    !! mode index
    INTEGER, INTENT(in) :: totf
    !! total final state
    REAL(KIND = DP), INTENT(in) :: indextot(totf,5)
    !! index buffer to find final state index
    REAL(KIND = DP), INTENT(in) :: n
    !! phonon occupation
    INTEGER, INTENT(in) :: test
    !! test case
    !
    ! Local variables
    !
    INTEGER :: l,numc,numv,numq,nktotf
    !! indices for many-body states
    INTEGER :: check_final
    !! final state index in the global array
    !
    !
    nktotf = nkf1 * nkf2 * nkf3
    IF (filkf /= '') THEN
      nktotf = nkqtotf / 2
    ELSE
      nktotf = nkf1 * nkf2 * nkf3
    ENDIF
    check_final = -1
    numc = (cbnd - 1) * nktotf + k_final
    numv = (vbnd - 1) * nktotf + k_init
    numq = (imode - 1) * nqtotf + iq
    !
    DO l = 1, totf
      !  
      IF (test == 2) THEN
        IF ((ABS(numc - INT(indextot(l, 1))) < eps8) .AND. & 
            (ABS(numv - INT(indextot(l, 2))) < eps8) .AND. &
            (ABS(numq - INT(indextot(l, 3))) < eps8) .AND. &
            (ABS(n - (indextot(l, 4))) < eps8)) THEN
          !  
          check_final = l
        ENDIF
      ELSE
        IF ((ABS(numc - INT(indextot(l, 1))) < eps8) .AND. &
            (ABS(numv - INT(indextot(l, 2))) < eps8)) THEN
          check_final = l
        ENDIF
      ENDIF
    ENDDO
    !---------------------------------------------------------------------------
    END FUNCTION check_final
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    FUNCTION pool_final(k_in)
    !---------------------------------------------------------------------------
    !! This function finds the processor/pool where this particular k_in will end
    !! up in case of mp_mesh_k
    USE io_global,     ONLY : stdout, ionode_id
    USE global_var,    ONLY : nkqtotf,nkf,bztoibz
    USE input,         ONLY : mp_mesh_k
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : eps4, eps8
    USE mp_global,     ONLY : npool
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: k_in
    !! input k-point
    INTEGER             :: rest, k, pool_final
    !! final pool index 
    !  
    k = k_in
    pool_final = CEILING(k / REAL(nkf))
    rest=(nkqtotf - (2 * (nkqtotf /(2 * npool))) * npool) / 2
    !
    IF (mp_mesh_k) THEN
      k = bztoibz(k)
      pool_final = CEILING((k) / REAL(nkf))
      !
    ENDIF
    !
    IF (k > (rest) * (nkf)) THEN
      !
      pool_final = CEILING((k - (rest) * nkf)/REAL(((nkqtotf / 2) - rest * nkf) / &
                  (npool - rest))) + rest
      !
    ENDIF
    !--------------------------------------------------------------------------
    END FUNCTION pool_final
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE calc_CHBB(nomega,iq)
    !------------------------------=-------------------------------------------
    !! CHBB theory for indirect and direct absorption but using second
    !! quantization. Ideally, should provide same result as indabs but due to
    !! Fermi-function mismatch and Pauli-blocking, it provides slightly (order of
    !! 1E-6) different result.
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nstemp, fsthick, degaussw,                      &
                              eps_acoustic, efermi_read, fermi_energy,        &
                              vme, omegamin, omegamax, omegastep
    USE global_var,    ONLY : etf, ibndmin, nkf, epf17, wkf, nqtotf, wf, wqf, &
                              sigmar_all, efnew, gtemp, nkqtotf, nkqf,        &
                              omegap, vmef, nbndfst, nktotf, xkf, xqf,        &
                              epsilon2_direct, epsilon2_indirect,             &
                              epsilon2_qdirect, E_mesh
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, &
                              ci, eps6, czero, eps8
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, npool, my_pool_id
    USE cell_base,     ONLY : omega, bg, at
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: nomega
    !! number of photon grid points
    INTEGER, INTENT(in)  :: iq
    !! q point
    !
    ! Local variables
    !
    INTEGER :: ikk,ikq,ik,cbnd,vbnd,mbnd,iw,imode
    !! Index for k-point, k+q point, global-k, conduction valence and empty
    !! bands, photon index, mode index  
    REAL(KIND = DP) :: alpha
    !! +1 or -1
    INTEGER :: ipol
    !! Polarization direction
    INTEGER :: m
    !! Counter on denominator imaginary broadening values
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    INTEGER, PARAMETER :: neta = 9
    !! Broadening parameter
    REAL(KIND = DP) :: ekkcb
    !! Eigen energy on the fine grid relative to the Fermi level for conduction
    !! bands
    REAL(KIND = DP) :: ekqcb
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for
    !conduction bands
    REAL(KIND = DP) :: ekkvb
    !! Eigen energy on the fine grid relative to the Fermi level for valence
    !! bands
    REAL(KIND = DP) :: ekqvb
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for
    !valence bands
    REAL(KIND = DP) :: ekmk
    !! Eigen energy on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: ekmq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: wq(nmodes),nqv(nmodes)
    !! Phonon frequencies and phonon occupations on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wgkk, wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$, $f_{nk}(T)$
    REAL(KIND = DP) :: weighta, weighte, weightd, weightq
    !!- delta function for absorption, emission, direct
    REAL(KIND = DP) :: inv_eptemp0
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: pfac
    !! Occupation prefactors
    REAL(KIND = DP) :: pface
    !! Occupation prefactors
    REAL(KIND = DP) :: cfac
    !! Absorption prefactor
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: eta(neta) = (/ 0.000, 0.002, 0.005, 0.01, 0.02, & 
                      0.05, 0.1, 0.2, 0.5 /) / ryd2eV
    !! Broadening parameter
    COMPLEX(KIND = DP) :: vkk(3, nbndfst, nbndfst)
    !!- Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: vkq(3, nbndfst, nbndfst)
    !!- Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: Aa(3), Ae(3), Ba(3), Be(3), Ca(3), Ce(3),    &
                         Da(3), De(3), xxk(3), xxkq(3), Qa(3)
    !! The four terms in conduction and valence band notation
    COMPLEX(KIND = DP) :: epf(nbndfst, nbndfst, nmodes)
    !! Electron-Phonon matrix
    !
    ef0=fermi_energy
    !
    inv_degaussw = 1.0 / degaussw
    !
    ! Epsilon2 prefactor for velocity matrix elements, including factor of 2 
    ! for spin (weights for k-points are divided by 2 to be normalized to 1)
    ! C = 8*pi^2*e^2 = 8*pi^2*2 approx 157.9136704
    !
    cfac = 16.d0 * pi**2
    !
    DO itemp = 1, nstemp
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        DO imode = 1, nmodes
          !
          ! the phonon frequency at this q and nu
          wq(imode) = wf(imode, iq)
          !
          epf(:, :, imode) = epf17(:, :, imode,ik)
          IF (wq(imode) > eps_acoustic) THEN
            nqv(imode) = wgauss(-wq(imode) / gtemp(itemp), -99)
            nqv(imode) = nqv(imode) / (one - two * nqv(imode))
          ENDIF
        ENDDO
        !
        DO cbnd = 1, nbndfst
          DO vbnd = 1, nbndfst
            vkk(:, cbnd, vbnd) = vmef(:, ibndmin - 1 + cbnd, ibndmin - 1 + vbnd, ikk)
            vkq(:, cbnd, vbnd) = vmef(:, ibndmin - 1 + cbnd, ibndmin - 1 + vbnd, ikq)
          ENDDO
        ENDDO
        !
        DO cbnd = 1, nbndfst
          !  the energy of the electron at k (relative to Ef)
          ekkcb = etf(ibndmin - 1 + cbnd, ikk) - ef0
          ekqcb = etf(ibndmin - 1 + cbnd, ikq) - ef0
          IF ((ekkcb .le. zero) .AND. (ekqcb .le. zero)) CYCLE
          !
          IF (ABS(ekkcb) < fsthick) THEN
            !
            wgkk = wgauss(-ekkcb * inv_eptemp0, -99)
            !
            DO vbnd = 1, nbndfst
              !
              ekkvb = etf(ibndmin - 1 + vbnd, ikk) - ef0
              ekqvb = etf(ibndmin - 1 + vbnd, ikq) - ef0
              IF ((ekkvb > zero) .AND. (ekqvb > zero)) CYCLE
              !
              IF ((ABS(ekqcb) < fsthick) .AND. (ekqcb < ekkvb + wq(nmodes) + omegamax + 6.0 * degaussw)) THEN
                DO imode = 1, nmodes
                  !
                  IF (wq(imode) > eps_acoustic) THEN
                    !
                    DO m = 1, neta
                      !  
                      Aa = czero
                      Ba = czero
                      Ca = czero
                      Da = czero
                      Ae = czero
                      Be = czero
                      Ce = czero
                      De = czero
                      !
                      DO mbnd = 1, nbndfst
                        !
                        ! The energy of the electron at k (relative to Ef)
                        ekmk = etf(ibndmin - 1 + mbnd, ikk) - ef0
                        ! The energy of the electron at k+q (relative to Ef)
                        ekmq = etf(ibndmin - 1 + mbnd, ikq) - ef0
                        !
                        ! c-c' transitions
                        IF ((ekmk > zero) .AND. (ekmq > zero)) THEN
                          Aa(:) = Aa(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd) &
                                  / (ekqcb - ekmk + wq(imode) + ci * eta(m))
                          !  
                          Ae(:) = Ae(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd) &
                                  / (ekqcb - ekmk - wq(imode) + ci * eta(m))
                          !  
                          Ca(:) = Ca(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode) &
                                  / (ekkvb - ekmq - wq(imode) + ci * eta(m))
                          !  
                          Ce(:) = Ce(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode) &
                                  / (ekkvb - ekmq + wq(imode) + ci * eta(m))
                          !
                        ! v-v' transitions
                        ELSEIF ((ekmk < zero) .AND. (ekmq < zero)) THEN
                          Ba(:) = Ba(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode) &
                                  / (ekmq - ekkvb + wq(imode) + ci * eta(m))
                          !  
                          Be(:) = Be(:) + vkq(:, cbnd, mbnd) * epf(mbnd, vbnd, imode) &
                                  / (ekmq - ekkvb - wq(imode) + ci * eta(m))
                          !
                          Da(:) = Da(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd) &
                                  / (ekmk - ekqcb - wq(imode) + ci * eta(m))
                          !
                          De(:) = De(:) + epf(cbnd, mbnd, imode) * vkk(:, mbnd, vbnd) &
                                  / (ekmk - ekqcb + wq(imode) + ci * eta(m))
                          !  
                        ENDIF
                      ENDDO
                      !
                      pface  =  nqv(imode)
                      pfac   = (nqv(imode) + one)
                      !
                      DO iw = 1, nomega
                        !
                        weighte = w0gauss((ekqcb - ekkvb - omegap(iw) - wq(imode)) / degaussw, 0) / degaussw
                        weighta = w0gauss((ekqcb - ekkvb - omegap(iw) + wq(imode)) / degaussw, 0) / degaussw
                        weightd = w0gauss((ekkcb-ekkvb-omegap(iw)) / (degaussw), 0) / (degaussw)
                        DO ipol = 1, 3
                          epsilon2_indirect(ipol, iw, m, itemp) = epsilon2_indirect(ipol, iw, m, itemp) + &
                                                                             (wkf(ikk) / 2.0) * wqf(iq) * &
                                                               ((cfac / omegap(iw)**2 * pfac  * weighta * &
                              ABS(Aa(ipol) + Ba(ipol) + Ca(ipol) + Da(ipol))**2 / (2 * wq(imode) * omega)))
                          !
                          epsilon2_indirect(ipol, iw, m, itemp) = epsilon2_indirect(ipol, iw, m, itemp) + &
                                                                             (wkf(ikk) / 2.0) * wqf(iq) * &
                                                               ((cfac / omegap(iw)**2 * pface * weighte * &
                              ABS(Ae(ipol) + Be(ipol) + Ce(ipol) + De(ipol))**2 / (2 * wq(imode) * omega)))
                          !
                          IF ((iq == nqtotf) .AND. (imode == 1) .AND. (m == 1)) THEN
                            epsilon2_direct(ipol, iw, 1, itemp) = epsilon2_direct(ipol, iw, 1, itemp) + &
                                                                            ((cfac / (omegap(iw)**2)) * &
                                        weightd * ABS(vkk(ipol, cbnd, vbnd))**2 / omega) * (wkf(ikk)/2.0)
                            !
                            epsilon2_direct(ipol, iw, 2, itemp) = epsilon2_direct(ipol, iw, 2, itemp) + &
                              ((cfac) * weightd * ABS(vkk(ipol, cbnd, vbnd))**2 / omega) * (wkf(ikk)/2.0)
                            !
                            epsilon2_direct(ipol, iw, 3, itemp) = epsilon2_direct(ipol,iw, 3, itemp) + &
                              ((cfac) * ABS(vkk(ipol,cbnd,vbnd))**2 / omega) * (wkf(ikk)/2.0)
                            !
                          ENDIF
                        ENDDO  !ipol
                      ENDDO ! nomega
                    ENDDO ! neta
                  ENDIF ! if wq > acoustic
                ENDDO ! imode
              ENDIF ! endif  ekq in fsthick
            ENDDO ! vbnd
          ENDIF  ! endif  ekk in fsthick
        ENDDO ! cbnd
      ENDDO ! ik
    ENDDO ! itemp
    !--------------------------------------------------------------------------
    END SUBROUTINE calc_CHBB
    !--------------------------------------------------------------------------
    !
  !----------------------------------------------------------------------------
  END MODULE qdabs
  !----------------------------------------------------------------------------

