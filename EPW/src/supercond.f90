  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Roxana Margine
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE supercond
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the routines linked with superconductivity using  
  !! the isotropic or anisotropic Eliashberg formalism. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    !                                                                            
    !----------------------------------------------------------------------- 
    SUBROUTINE eliashberg_init()
    !-----------------------------------------------------------------------
    !
    ! This routine initializes the control variables needed to solve the eliashberg 
    ! equations
    !
    USE kinds,           ONLY : DP
    USE io_global,       ONLY : stdout
    USE epwcom,          ONLY : eliashberg, nkf1, nkf2, nkf3, nsiter, &
                                nqf1, nqf2, nqf3, ntempxx, nswi, nstemp, temps, &
                                muc, lreal, lpade, liso, limag, laniso, lacon, &
                                kerwrite, kerread, imag_read, fila2f, wsfc, wscut, & 
                                tempsmin, tempsmax, rand_q, rand_k
    USE constants_epw,   ONLY : kelvin2eV
    USE eliashbergcom,   ONLY : estemp
    !
    IMPLICIT NONE
    !
    INTEGER :: itemp
    !! Counter on temperature values
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: dtemp
    !! Step in temperature
    !
    IF (eliashberg .AND. liso .AND. laniso) & 
      CALL errore('eliashberg_init', 'liso or laniso needs to be true', 1)
    IF (.NOT. eliashberg .AND. liso) &
      CALL errore('eliashberg_init', 'liso requires eliashberg true', 1)
    IF (.NOT. eliashberg .AND. laniso) & 
      CALL errore('eliashberg_init', 'laniso requires eliashberg true', 1)
    IF (laniso .AND. (fila2f /= ' ')) &
      CALL errore('eliashberg_init', 'anisotropic case can not use fila2f', 1)
    IF (eliashberg .AND. lreal .AND. laniso) & 
      CALL errore('eliashberg_init', 'lreal is implemented only for the isotriopic case', 1)
    IF (eliashberg .AND. lreal .AND. limag) &
      CALL errore('eliashberg_init', 'lreal or limag needs to be true', 1)
    IF (eliashberg .AND. lreal .AND. lacon) &
      CALL errore('eliashberg_init', 'lreal or lacon needs to be true', 1)
    IF (eliashberg .AND. lreal .AND.  lpade) &
      CALL errore('eliashberg_init', 'lreal or lpade needs to be true', 1)
    IF (eliashberg .AND. imag_read .AND. .NOT. limag .AND. .NOT. laniso) &
      CALL errore('eliashberg_init', 'imag_read requires limag true and laniso true', 1)
    IF (eliashberg .AND. lpade .AND. .NOT. limag) &                
      CALL errore('eliashberg_init', 'lpade requires limag true', 1)
    IF (eliashberg .AND. lacon .AND. (.NOT. limag .OR. .NOT. lpade)) & 
      CALL errore('eliashberg_init', 'lacon requires both limag and lpade true', 1)
    IF (eliashberg .AND. lreal .AND. (kerread .AND. kerwrite)) & 
      CALL errore('eliashberg_init', 'kerread cannot be used with kerwrite', 1)
    IF (eliashberg .AND. lreal .AND. (.NOT. kerread .AND. .NOT. kerwrite)) &
      CALL errore('eliashberg_init', 'kerread or kerwrite must be true', 1)
    IF (eliashberg .AND. lreal .AND. wsfc > wscut) & 
      CALL errore('eliashberg_init', 'wsfc should be < wscut', 1)
    IF (eliashberg .AND. lreal .AND. wsfc < 0.d0) & 
      CALL errore('eliashberg_init', 'wsfc should be > 0.d0', 1)
    IF (eliashberg .AND. nswi > 0 .AND. .NOT. limag) &
      CALL errore('eliashberg_init', 'nswi requires limag true', 1)
    IF (eliashberg .AND. nswi < 0) & 
      CALL errore('eliashberg_init', 'nswi should be > 0', 1)
    IF (eliashberg .AND. wscut < 0.d0 ) &
      CALL errore('eliashberg_init', 'wscut should be > 0.d0', 1)
    IF (eliashberg .AND. nstemp < 1) & 
      CALL errore('eliashberg_init', 'wrong number of nstemp', 1)
    IF (eliashberg .AND. MAXVAL(temps(:)) > 0.d0 .AND. tempsmin > 0.d0 .AND. tempsmax > 0.d0) &
      CALL errore('eliashberg_init', 'define either (tempsmin and tempsmax) or temps(:)', 1)
    IF (eliashberg .AND. tempsmax < tempsmin) &
      CALL errore('eliashberg_init', 'tempsmax should be greater than tempsmin', 1)
    IF (eliashberg .AND. nsiter < 1) &
      CALL errore('eliashberg_init', 'wrong number of nsiter', 1)
    IF (eliashberg .AND. muc < 0.d0) &
      CALL errore('eliashberg_init', 'muc should be >= 0.d0', 1) 
    IF (eliashberg .AND. (rand_k .OR. rand_q) .AND. (fila2f == ' ')) &
      CALL errore('eliashberg_init', 'eliashberg requires a uniform grid when fila2f is not used', 1)
    IF (eliashberg .AND. (MOD(nkf1, nqf1) /= 0 .OR. MOD(nkf2, nqf2) /= 0 .OR. MOD(nkf3, nqf3) /= 0 ) .AND. (fila2f == ' ')) &
      CALL errore('eliashberg_init', &
                  'eliashberg requires nkf1,nkf2,nkf3 to be multiple of nqf1,nqf2,nqf3 when fila2f is not used', 1)
    !
    DO itemp = 1, ntempxx
      IF (temps(itemp) > 0.d0) THEN
        nstemp = itemp
      ENDIF
    ENDDO
    !
    ALLOCATE(estemp(nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_init', 'Error allocating estemp', 1)
    estemp(:) = 0.d0
    !
    ! go from K to eV
    IF (MAXVAL(temps(:)) > 0.d0) THEN
      DO itemp= 1, nstemp 
        estemp(itemp) = temps(itemp) * kelvin2eV
      ENDDO
    ELSE
      IF (nstemp == 1) THEN
        estemp(1) = tempsmin * kelvin2eV
      ELSE
        dtemp = (tempsmax - tempsmin) * kelvin2eV / DBLE(nstemp - 1)
        DO itemp = 1, nstemp
          estemp(itemp) = tempsmin * kelvin2eV + DBLE(itemp - 1) * dtemp
        ENDDO
      ENDIF
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_init
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_grid()
    !-----------------------------------------------------------------------
    !
    ! This routine initializes the point grids to solve the
    ! eliashberg equations on real and imag axis
    !
    USE kinds,           ONLY : DP
    USE io_global,       ONLY : stdout
    USE epwcom,          ONLY : nqstep, nswi, nswfc, nswc, nstemp, & 
                                lreal, lpade, limag, lacon, wsfc, wscut
    USE constants_epw,   ONLY : pi, eps6
    USE eliashbergcom,   ONLY : estemp, nsw, nsiw, wsphmax
    !
    IMPLICIT NONE
    !
    INTEGER :: itemp
    !! Counter on temperature values
    INTEGER :: ierr
    !! Error status
    !
    IF (lreal) THEN
      !
      IF (ABS(wsfc) < eps6 .OR. ABS(wscut) < eps6 .OR. nswfc == 0 .OR. nswc == 0) THEN 
        wsfc  = 4.d0 * wsphmax
        wscut = 8.d0 * wsphmax
        nswfc = 4 * nqstep
        nswc  = 2 * nqstep
      ENDIF
      nsw = nswfc + nswc  
      WRITE(stdout, '(5x, a7, f12.6, a11, f12.6)') 'wsfc = ', wsfc, '   wscut = ', wscut
      WRITE(stdout, '(5x, a8, i8, a10, i8, a9, i8)') 'nswfc = ', nswfc, '   nswc = ', nswc, & 
                                                 '   nsw = ', nsw 
      IF (nsw == 0) CALL errore('eliashberg_init', 'wrong number of nsw', 1)
      !
    ELSEIF (limag) THEN
      !
      ALLOCATE(nsiw(nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_init', 'Error allocating nsiw', 1)
      nsiw(:) = 0
      IF (nswi > 0) THEN
        nsiw(:) = nswi
      ELSEIF (wscut > 0.d0) THEN
        DO itemp = 1, nstemp
           nsiw(itemp) = int(0.5d0 * (wscut / pi / estemp(itemp) - 1.d0)) + 1
        ENDDO
      ELSEIF (nswi > 0 .AND. wscut > 0.d0) THEN
        nsiw(:) = nswi
        WRITE(stdout,'(5x,a)') 'when nswi > 0, wscut is not used for limag=.TRUE.'
      ENDIF
      !
      IF (ABS(wscut) < eps6) THEN 
        wscut = 10.d0 * wsphmax
      ENDIF
      ! 
      IF (lpade .OR. lacon) THEN
        nsw = nqstep * NINT(wscut / wsphmax)
        IF (nsw == 0) CALL errore('eliashberg_init', 'wrong number of nsw', 1)
      ENDIF
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_grid
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE evaluate_a2f_lambda
    !-----------------------------------------------------------------------
    !
    ! computes the isotropic spectral function a2F(w), total lambda, and 
    ! distribution of lambda
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iua2ffil, iudosfil, iufillambda, iufillambdaFS
    USE io_files,      ONLY : prefix
    USE phcom,         ONLY : nmodes
    USE cell_base,     ONLY : bg
    USE control_flags, ONLY : iverbosity
    USE elph2,         ONLY : nqtotf, wqf, wf
    USE epwcom,        ONLY : fsthick, eps_acustic, nqstep, degaussq, delta_qsmear, nqsmear, & 
                              degaussw, nkf1, nkf2, nkf3
    USE eliashbergcom, ONLY : nkfs, nbndfs, g2, ixkqf, ixqfs, nqfs, w0g, ekfs, ef0, dosef, wsph, &
                              wkfs, dwsph, a2f_iso, ixkff
    USE constants_epw, ONLY : ryd2ev, eps2, zero, eps16
    USE io_global,     ONLY : stdout, ionode_id
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,      ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: name1
    !! file name
    !
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: ibnd, jbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k or q parallelization
    INTEGER :: iwph
    !! Counter over frequncy
    INTEGER :: ismear
    !! Counter on smearings
    INTEGER :: ibin
    !! Counter on bins
    INTEGER :: nbin, nbink
    !! Number of bins
    INTEGER :: i, j, k
    !! Counter on grid points nkf1, nkf2, nkf3
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: dbin, dbink
    !! Step size in nbin/nbink
    REAL(KIND = DP):: l_sum
    !! total e-ph coupling strength
    REAL(KIND = DP):: lambda_eph
    !! total e-ph coupling strength (a2f integration)
    REAL(KIND = DP) :: x1, x2, x3
    !! Cartesian coordinates of grid points nkf1, nkf2, nkf3
    REAL(KIND = DP) :: weight, weightq
    !! factors in lambda_eph and a2f
    REAL(KIND = DP) :: sigma 
    !! smearing in delta function
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: lambda_max(npool)
    !! Max e-ph coupling strength among all pools
    REAL(KIND = DP), ALLOCATABLE :: l_a2f(:)
    !! total e-ph coupling strength for different ismear (a2f integration)
    REAL(KIND = DP), ALLOCATABLE :: lambda_pairs(:)
    !! Histogram lambda_nk,n'k
    REAL(KIND = DP), ALLOCATABLE :: lambda_k_bin(:)
    !! Histogram lambda_nk
    REAL(KIND = DP), ALLOCATABLE :: a2f(:, :)
    !! Eliashberg sperctral function for different ismear
    REAL(KIND = DP), ALLOCATABLE :: a2f_modeproj(:, :)
    !! Eliashberg sperctral function projected over modes for different ismear
    REAL(KIND = DP), ALLOCATABLE :: phdos(:, :)
    !! Phonon density of states for different ismear
    REAL(KIND = DP), ALLOCATABLE :: phdos_modeproj(:, :)
    !! Phonon density of states  projected over modes for different ismear
    REAL(KIND = DP), ALLOCATABLE :: lambda_k(:, :)
    !! anisotropic e-ph coupling strength 
    ! 
    ! This is only a quick fix since the routine was written for parallel execution - FG June 2014
#if defined(__MPI)
    npool = 1
    my_pool_id = 0
#endif
    !
    ! degaussq is read from the input file in meV and converted to Ryd in epw_readin.f90
    ! go from Ryd to eV
    degaussq = degaussq * ryd2ev
    delta_qsmear = delta_qsmear * ryd2ev
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ALLOCATE(a2f_iso(nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating a2f_iso', 1) 
    ALLOCATE(a2f(nqstep, nqsmear), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating a2f', 1) 
    ALLOCATE(a2f_modeproj(nmodes, nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating a2f_modeproj', 1) 
    a2f_iso(:) = zero
    a2f(:, :) = zero
    a2f_modeproj(:, :) = zero
    !
    ! RM - the 0 index in k is required when printing out values of lambda_k 
    ! When the k-point is outside the Fermi shell, ixkff(ik)=0
    ALLOCATE(lambda_k(0:nkfs, nbndfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating lambda_k', 1)
    lambda_k(:, :) = zero
    !
    l_sum = zero
    lambda_max(:) = zero
    DO ismear = 1, nqsmear
      sigma = degaussq + (ismear - 1) * delta_qsmear
      DO ik = lower_bnd, upper_bnd
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN 
            DO iq = 1, nqfs(ik)
              ! iq0 - index of q-point on the full q-mesh
              iq0 = ixqfs(ik, iq)
              DO jbnd = 1, nbndfs
                IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                  weight = wkfs(ik) * wqf(iq) * w0g(ibnd, ik) * w0g(jbnd, ixkqf(ik, iq0))
                  lambda_eph = 0.d0
                  DO imode = 1, nmodes
                    IF (wf(imode, iq0) > eps_acustic) THEN
                      IF (ismear == 1) THEN 
                        lambda_eph = lambda_eph + g2(ik, iq, ibnd, jbnd, imode) / wf(imode, iq0)
                      ENDIF
                      DO iwph = 1, nqstep
                        weightq  = w0gauss((wsph(iwph) - wf(imode, iq0)) / sigma, 0) / sigma
                        a2f(iwph, ismear) = a2f(iwph, ismear) + weight * weightq * g2(ik, iq, ibnd, jbnd, imode)
                        IF (ismear == 1) THEN
                          a2f_modeproj(imode, iwph) = a2f_modeproj(imode, iwph) + &
                                       weight * weightq * g2(ik, iq, ibnd, jbnd, imode)
                        ENDIF
                      ENDDO ! iwph
                    ENDIF ! wf
                  ENDDO ! imode
                  IF (ismear == 1 .AND. lambda_eph > 0.d0) THEN
                    l_sum = l_sum + weight * lambda_eph
                    weight = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0)) 
                    lambda_k(ik, ibnd) = lambda_k(ik, ibnd) + weight * lambda_eph
                    IF (lambda_eph > lambda_max(my_pool_id + 1)) THEN
                      lambda_max(my_pool_id + 1) = lambda_eph
                    ENDIF
                  ENDIF
                ENDIF ! ekq
              ENDDO ! jbnd
            ENDDO ! iq
          ENDIF ! ekk
        ENDDO ! ibnd
      ENDDO ! ik
    ENDDO ! ismear
    !
    a2f(:, :) = 0.5d0 * a2f(:, :) / dosef
    a2f_modeproj(:, :) = 0.5d0 * a2f_modeproj(:, :) / dosef
    l_sum = l_sum / dosef
    lambda_k(:, :) = 2.d0 * lambda_k(:, :)
    lambda_max(:) = 2.d0 * dosef * lambda_max(:)
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(l_sum, inter_pool_comm)
    CALL mp_sum(a2f, inter_pool_comm)
    CALL mp_sum(a2f_modeproj, inter_pool_comm)
    CALL mp_sum(lambda_max, inter_pool_comm)
    CALL mp_sum(lambda_k, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime == ionode_id) THEN
      !
      ALLOCATE(phdos(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating phdos', 1)
      ALLOCATE(phdos_modeproj(nmodes, nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating phdos_modeproj', 1)
      phdos(:, :) = zero
      phdos_modeproj(:, :) = zero
      !
      DO ismear = 1, nqsmear
        sigma = degaussq + (ismear - 1) * delta_qsmear
        DO iq = 1, nqtotf
          DO imode = 1, nmodes
            IF (wf(imode, iq) > eps_acustic) THEN
              DO iwph = 1, nqstep
                weightq  = w0gauss((wsph(iwph) - wf(imode, iq)) / sigma, 0) / sigma
                phdos(iwph, ismear) = phdos(iwph, ismear) + wqf(iq) * weightq
                IF (ismear == 1) THEN
                  phdos_modeproj(imode, iwph) = phdos_modeproj(imode, iwph) + wqf(iq) * weightq
                ENDIF
              ENDDO ! iwph
            ENDIF ! wf
          ENDDO ! imode
        ENDDO ! iq
      ENDDO ! ismear
      !
      name1 = TRIM(prefix) // '.a2f'
      OPEN(UNIT = iua2ffil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iua2ffil)
      name1 = TRIM(prefix) // '.phdos'
      OPEN(UNIT = iudosfil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iudosfil)
      !
      ALLOCATE(l_a2f(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating l_a2f', 1)
      l_a2f(:) = zero
      !
      DO ismear = 1, nqsmear
        IF (ismear == nqsmear) THEN
          WRITE(iua2ffil, '(" w[meV] a2f for ", i4, " smearing values")') ismear
          WRITE(iudosfil, '(" w[meV] phdos[states/meV] for ", i4, " smearing values")') ismear
        ENDIF
        DO iwph = 1, nqstep
          l_a2f(ismear) = l_a2f(ismear) + a2f(iwph, ismear) / wsph(iwph)
          ! wsph in meV (from eV) and phdos in states/meV (from states/eV)
          IF (ismear == nqsmear) THEN
            WRITE(iua2ffil, '(f12.7, 15f12.7)') wsph(iwph) * 1000.d0, a2f(iwph, :)
            WRITE(iudosfil, '(f12.7, 15f15.7)') wsph(iwph) * 1000.d0, phdos(iwph, :)/ 1000.d0
          ENDIF
        ENDDO
        l_a2f(ismear) = 2.d0 * l_a2f(ismear) * dwsph
      ENDDO
      !
      WRITE(iua2ffil, *) "Integrated el-ph coupling"
      WRITE(iua2ffil, '("  #         ", 15f12.7)') l_a2f(:)
      WRITE(iua2ffil, *) "Phonon smearing (meV)" 
      WRITE(iua2ffil, '("  #         ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * 1000.d0, ismear = 1, nqsmear)
      WRITE(iua2ffil, '("Electron smearing (eV)", f12.7)') degaussw
      WRITE(iua2ffil, '("Fermi window (eV)", f12.7)') fsthick
      WRITE(iua2ffil, '("Summed el-ph coupling ", f12.7)') l_sum
      CLOSE(iua2ffil)
      CLOSE(iudosfil)
      !
      a2f_iso(:) = a2f(:, 1)
      name1 = TRIM(prefix) // '.a2f_iso'
      OPEN(UNIT = iua2ffil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iua2ffil)
      WRITE(iua2ffil, '("w[meV] a2f")')
      DO iwph = 1, nqstep
        ! wsph in meV (from eV) and phdos in states/meV (from states/eV)
        WRITE(iua2ffil, '(f12.7, 100f12.7)') wsph(iwph) * 1000.d0, a2f_iso(iwph)
      ENDDO
      CLOSE(iua2ffil)
      !
      name1 = TRIM(prefix) // '.a2f_proj'
      OPEN(UNIT = iua2ffil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iua2ffil)
      name1 = TRIM(prefix) // '.phdos_proj'
      OPEN(UNIT = iudosfil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iudosfil)
      !
      WRITE(iua2ffil, '("w[meV] a2f a2f_modeproj")')
      WRITE(iudosfil, '("w[meV] phdos[states/meV] phdos_modeproj[states/meV]")') 
      DO iwph = 1, nqstep
        ! wsph in meV (from eV) and phdos in states/meV (from states/eV)
        WRITE(iua2ffil, '(f12.7, 100f12.7)') wsph(iwph) * 1000.d0, a2f_iso(iwph), a2f_modeproj(:, iwph)
        WRITE(iudosfil, '(f12.7, 100f15.7)') wsph(iwph) * 1000.d0, phdos(iwph, 1)/1000.d0, phdos_modeproj(:, iwph) / 1000.d0
      ENDDO
      WRITE(iua2ffil, '(a, f18.7, a, f18.7)') 'lambda_int = ', l_a2f(1), '   lambda_sum = ',l_sum
      CLOSE(iua2ffil)
      CLOSE(iudosfil)
      !
      DEALLOCATE(phdos, STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating phdos', 1)
      DEALLOCATE(phdos_modeproj, STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating phdos_modeproj', 1)
      DEALLOCATE(l_a2f, STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating l_a2f', 1)
      !
    ENDIF
    !
    CALL mp_bcast(a2f_iso, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(a2f, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating a2f', 1)
    DEALLOCATE(a2f_modeproj, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating a2f_modeproj', 1)
    !
    nbink = NINT(1.1d0 * MAXVAL(lambda_k(:, :)) / eps2) + 1 
    dbink = 1.1d0 * MAXVAL(lambda_k(:, :)) / DBLE(nbink) 
    !
    ALLOCATE(lambda_k_bin(nbink), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating lambda_k_bin', 1)
    lambda_k_bin(:) = zero
    !
    nbin = 0
    dbin = zero
    IF (iverbosity == 2) THEN
      nbin = NINT(1.1d0 * MAXVAL(lambda_max(:)) / eps2) + 1
      dbin = 1.1d0 * MAXVAL(lambda_max(:)) / DBLE(nbin)
      ALLOCATE(lambda_pairs(nbin), STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating lambda_pairs', 1)
      lambda_pairs(:) = zero
    ENDIF
    ! 
    WRITE(stdout, '(5x, a13, f21.7, a18, f21.7)') 'lambda_max = ', MAXVAL(lambda_max(:)), & 
                                             '   lambda_k_max = ', MAXVAL(lambda_k(:, :))
    WRITE(stdout, '(a)') ' '
    !
    lambda_k(:, :) = zero
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0 ) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                weight = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0)) / dosef
                lambda_eph = zero
                DO imode = 1, nmodes  
                  IF (wf(imode, iq0) > eps_acustic) THEN
                    lambda_eph = lambda_eph + g2(ik, iq, ibnd, jbnd, imode) / wf(imode,iq0) 
                  ENDIF
                ENDDO
                lambda_eph = 2.d0 * lambda_eph * dosef
                lambda_k(ik, ibnd) = lambda_k(ik, ibnd) +  weight * lambda_eph
                IF (iverbosity == 2) THEN
                  ibin = NINT(lambda_eph / dbin) + 1
                  weight =  w0g(ibnd, ik) * w0g(jbnd,ixkqf(ik, iq0))
                  lambda_pairs(ibin) = lambda_pairs(ibin) + weight
                ENDIF
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
          ibin = NINT(lambda_k(ik, ibnd) / dbink) + 1
          weight = w0g(ibnd, ik)
          lambda_k_bin(ibin) = lambda_k_bin(ibin) + weight
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    ! collect contributions from all pools 
    CALL mp_sum(lambda_k, inter_pool_comm)
    IF (iverbosity == 2) THEN  
      CALL mp_sum(lambda_pairs, inter_pool_comm)
    ENDIF
    CALL mp_sum(lambda_k_bin, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime == ionode_id) THEN
      !
      ! SP: Produced if user really wants it 
      IF (iverbosity == 2) THEN
        name1 = TRIM(prefix) // '.lambda_aniso'
        OPEN(UNIT = iufillambda, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iufillambda)
        WRITE(iufillambda, '(2a12, 2a7)') '# Enk-Ef[eV]', '  lambda_nk', '# kpt', '# band'
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              WRITE(iufillambda, '(2f12.7, 2i7)') ekfs(ibnd, ik) - ef0, lambda_k(ik, ibnd), ik, ibnd
            ENDIF
          ENDDO
        ENDDO
        CLOSE(iufillambda)
      ENDIF
      !
      name1 = TRIM(prefix) // '.lambda_k_pairs'
      OPEN(UNIT = iufillambda, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iufillambda)
      WRITE(iufillambda, '(a12, a30)') '# lambda_nk','  \rho(lambda_nk) scaled to 1'
      DO ibin = 1, nbink
        WRITE(iufillambda,'(2f21.7)') dbink * DBLE(ibin), lambda_k_bin(ibin) / MAXVAL(lambda_k_bin(:))
      ENDDO
      CLOSE(iufillambda)
      !
      ! SP: Produced if user really wants it 
      IF (iverbosity == 2) THEN  
        name1 = TRIM(prefix) // '.lambda_pairs'
        OPEN(UNIT = iufillambda, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iufillambda)
        WRITE(iufillambda, '(a12, a30)') "# lambda_nk,n'k'", "  \rho(lambda_nk,n'k') scaled to 1"
        DO ibin = 1, nbin
          WRITE(iufillambda, '(2f21.7)') dbin * DBLE(ibin), lambda_pairs(ibin) / MAXVAL(lambda_pairs(:))
        ENDDO
        CLOSE(iufillambda)
      ENDIF
      !
      ! SP & RM: .cube file for VESTA plotting (only if iverbosity = 2)
      !
      ! RM - If the k-point is outside the Fermi shell,
      ! ixkff(ik)=0 and lambda_k(0,ibnd) = 0.0
      !
      IF (iverbosity == 2) THEN
        !
        DO ibnd = 1, nbndfs
          !
          IF (ibnd < 10) THEN
            WRITE(name1, '(a, a8, i1, a5)') TRIM(prefix), '.lambda_', ibnd, '.cube'
          ELSEIF (ibnd < 100) THEN
            WRITE(name1, '(a, a8, i2, a5)') TRIM(prefix), '.lambda_', ibnd, '.cube'
          ELSEIF( ibnd < 1000) THEN
            WRITE(name1,'(a, a8, i3, a5)') TRIM(prefix), '.lambda_', ibnd, '.cube'
          ELSE 
            CALL errore( 'eliashberg_setup', 'Too many bands ',1)  
          ENDIF  
          !  
          OPEN(iufillambdaFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
          IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iufillambdaFS)
          WRITE(iufillambdaFS, *) 'Cubfile created from EPW calculation'
          WRITE(iufillambdaFS, *) 'lambda'
          WRITE(iufillambdaFS, '(i5, 3f12.6)') 1, 0.0d0, 0.0d0, 0.0d0
          WRITE(iufillambdaFS, '(i5, 3f12.6)') nkf1, (bg(i, 1) / DBLE(nkf1), i = 1, 3)
          WRITE(iufillambdaFS, '(i5, 3f12.6)') nkf2, (bg(i, 2) / DBLE(nkf2), i = 1, 3)
          WRITE(iufillambdaFS, '(i5, 3f12.6)') nkf3, (bg(i, 3) / DBLE(nkf3), i = 1, 3)
          WRITE(iufillambdaFS, '(i5, 4f12.6)') 1, 1.0d0, 0.0d0, 0.0d0, 0.0d0
          WRITE(iufillambdaFS, '(6f12.6)') (lambda_k(ixkff(ik), ibnd), ik = 1, nkf1 * nkf2 * nkf3)
          CLOSE(iufillambdaFS)
        ENDDO
        !
      ENDIF
      !
      ! SP & RM : Write on file the lambda close to the Fermi surface along with 
      ! Cartesian coordinate, band index, energy distance from Fermi level
      ! and lambda value.
      !
      name1 = TRIM(prefix) // '.lambda_FS'
      OPEN(iufillambdaFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iufillambdaFS)
      WRITE(iufillambdaFS,'(a75)') '#               k-point                  Band Enk-Ef [eV]            lambda'
      DO i = 1, nkf1
        DO j = 1, nkf2
          DO k = 1, nkf3
            ik = k + (j - 1) * nkf3 + (i - 1) * nkf2 * nkf3
            IF (ixkff(ik) > 0) THEN
              DO ibnd = 1, nbndfs
                ! SP: Here take a 0.2 eV interval around the FS.
                IF (ABS(ekfs(ibnd, ixkff(ik)) - ef0) < fsthick) THEN
                !IF (ABS(ekfs(ibnd, ixkff(ik)) - ef0) < 0.2) THEN
                  x1 = bg(1, 1) * (i - 1) / nkf1 + bg(1, 2) * (j - 1) / nkf2 + bg(1, 3) * (k - 1) / nkf3
                  x2 = bg(2, 1) * (i - 1) / nkf1 + bg(2, 2) * (j - 1) / nkf2 + bg(2, 3) * (k - 1) / nkf3
                  x3 = bg(3, 1) * (i - 1) / nkf1 + bg(3, 2) * (j - 1) / nkf2 + bg(3, 3) * (k - 1) / nkf3
                  WRITE(iufillambdaFS, '(3f12.6, i8, f12.6, f24.15)') x1, x2, x3, ibnd, &
                                   ekfs(ibnd, ixkff(ik)) - ef0, lambda_k(ixkff(ik), ibnd)
                ENDIF
              ENDDO ! ibnd
            ENDIF
          ENDDO  ! k
        ENDDO ! j
      ENDDO ! i
      CLOSE(iufillambdaFS)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(lambda_k, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating lambda_k', 1)
    DEALLOCATE(lambda_k_bin, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating lambda_k_bin', 1)
    IF (iverbosity == 2) THEN
      DEALLOCATE(lambda_pairs, STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating lambda_pairs', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE evaluate_a2f_lambda
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------- 
    SUBROUTINE estimate_tc_gap()
    !-----------------------------------------------------------------------
    !
    ! This routine estimates the Tc using Allen-Dynes formula and 
    ! the BCS superconducting gap as the initial guess for delta 
    !  
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nqstep, muc, nstemp
    USE eliashbergcom, ONLY : estemp, wsph, dwsph, a2f_iso, gap0
    USE constants_epw, ONLY : kelvin2eV, zero
    USE io_global, ONLY : stdout, ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    !  
    IMPLICIT NONE
    !  
    INTEGER :: iwph
    !! Counter on frequencies
    ! 
    REAL(KIND = DP):: l_a2f
    !! total e-ph coupling strength (a2f integration)
    REAL(KIND = DP):: logavg
    !! logavg phonon frequency 
    REAL(KIND = DP):: tc
    !! superconding critical temperature
    !
    IF (mpime == ionode_id) THEN
      l_a2f  = zero
      logavg = zero
      DO iwph = 1, nqstep
        l_a2f  = l_a2f  + a2f_iso(iwph) / wsph(iwph)
        logavg = logavg + a2f_iso(iwph) * log(wsph(iwph)) / wsph(iwph)
      ENDDO
      l_a2f  = l_a2f  * 2.d0 * dwsph
      logavg = logavg * 2.d0 * dwsph
      logavg = EXP(logavg / l_a2f)
      WRITE(stdout,'(5x,a,f12.7)') 'Electron-phonon coupling strength = ', l_a2f
      WRITE(stdout,'(a)') ' '
      !
      ! Allen-Dynes estimate of Tc
      !
      tc = logavg / 1.2d0 * EXP(-1.04d0 * (1.d0 + l_a2f) &
                                / (l_a2f - muc * (1.d0 + 0.62d0 * l_a2f)))
      !
      ! initial guess for the gap edge using BCS superconducting ratio 3.52
      !
      gap0 = 3.52d0 * tc / 2.d0
      IF (gap0 <= 0.d0) & 
        CALL errore('estimate_tc_gap', 'initial guess for gap edge should be > 0.d0', 1)
      !
      ! tc in K
      !
      tc = tc / kelvin2eV
      WRITE(stdout, '(5x, a, f12.6, a, f10.5)') 'Estimated Allen-Dynes Tc = ', tc, ' K for muc = ', muc
      WRITE(stdout, '(a)') '  '
      WRITE(stdout, '(5x, a, f12.6, a)') 'Estimated BCS superconducting gap = ', gap0, ' eV'
      WRITE(stdout, '(a)') '  '
      !
      IF (estemp(1) / kelvin2eV > tc) THEN
        WRITE(stdout, '(a)') '  '
        WRITE(stdout, '(5x, a)') 'WARNING WARNING WARNING '
        WRITE(stdout, '(a)') '  '
        WRITE(stdout, '(5x, a, f9.3, a, f9.3, a)') 'The code may crash since tempsmin =', & 
                        estemp(1) / kelvin2eV, ' K is larger than Allen-Dynes Tc = ', tc, ' K'
      ELSEIF (estemp(nstemp) / kelvin2eV > tc) THEN
        WRITE(stdout, '(a)') '  '
        WRITE(stdout, '(5x, a)') 'WARNING WARNING WARNING '
        WRITE(stdout, '(a)') '  '
        WRITE(stdout, '(5x, a, f9.3, a, f9.3, a)') 'The code may crash since tempsmax =', & 
                        estemp(nstemp) / kelvin2eV, ' K is larger than Allen-Dynes Tc = ', tc, ' K'
      ENDIF
      !
    ENDIF
    CALL mp_bcast(gap0, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE estimate_tc_gap
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gen_freqgrid_raxis()
    !-----------------------------------------------------------------------
    !!
    !! Automatic generation of the frequency-grid for real-axis calculations.
    !!
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nswfc, nswc, pwc, wsfc, wscut, lunif
    USE eliashbergcom, ONLY : nsw, ws, dws
    USE constants_epw, ONLY : zero
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iw
    !! Counter over frequency
    INTEGER :: ierr
    !! Error status
    !
    ! define a grid ws in 2 step sizes
    ! 1. a fine grid of nswfc points between (0,wsfc)
    ! 2. a rough grid of nswc points between (wsfc,wscut).
    ! above wsfc the gap function varies slowly
    !
    ! nswfc = nr. of grid points between (0,wsfc)
    ! nswc  = nr. of grid points between (wsfc,wscut)
    !
    WRITE(stdout, '(a)') '    '
    WRITE(stdout, '(5x, a, i6, a)') 'Total number of nsw = ', nsw, ' grid-points are divided in:'
    WRITE(stdout, '(5x, a, i6, a, f12.6, a, f12.6)') 'nswfc = ', nswfc, '  from ', 0.0, ' to ', wsfc  
    WRITE(stdout, '(5x, a, i6, a, f12.6, a, f12.6)') 'nswc  = ', nswc,  '  from ', wsfc, ' to ', wscut
    WRITE(stdout, '(a)') '    '
    !
    ALLOCATE(ws(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('gen_freqgrid_raxis', 'Error allocating ws', 1)
    ALLOCATE(dws(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('gen_freqgrid_raxis', 'Error allocating dws', 1)
    ws(:) = zero
    dws(:) = zero
    !
    DO iw = 1, nswfc
      dws(iw) = wsfc / DBLE(nswfc)
      ws(iw) = DBLE(iw) * dws(iw)
    ENDDO
    DO iw = nswfc + 1, nsw
      dws(iw) = (wscut - wsfc) / DBLE(nswc)
      IF (lunif) THEN 
        ws(iw) = wsfc + DBLE(iw) * dws(iw)
      ELSE 
        ! RM this needs to be checked
        ws(iw) = wsfc + DBLE(iw / nswc)**pwc * (wscut - wsfc)
      ENDIF
    ENDDO
    !
    IF (.NOT. lunif) THEN 
      DO iw = nswfc + 1, nsw - 1
        dws(iw) = ws(iw + 1) - ws(iw)
      ENDDO
      dws(nsw) = dws(nsw - 1)
    ENDIF
    !
    DO iw = 1, nsw
      ! end points contribute only half (trapezoidal integration rule)
      IF ((iw == 1) .OR. (iw == nsw)) THEN
        dws(iw) = 0.5d0 * dws(iw)
      ! boundary points contribute half from left and half from right side
      ELSEIF (iw == nswfc) THEN
        dws(iw) = 0.5d0 * (dws(iw) + dws(iw + 1))
      ENDIF
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gen_freqgrid_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gen_freqgrid_iaxis(itemp)
    !-----------------------------------------------------------------------
    !
    ! Automatic generation of the frequency-grid for imaginary-axis calculations.
    !
    !
    ! input
    !
    ! itemp  - temperature point
    !
    USE epwcom,        ONLY : nqstep, lpade, lacon, laniso
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, wsph, dwsph, estemp, wsphmax
    USE constants_epw, ONLY : pi, zero
    USE low_lvl,       ONLY : mem_size_eliashberg
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    INTEGER :: iw
    !! Counter over frequency   
    INTEGER :: n
    !! frequency index - 1
    INTEGER :: imelt
    !! Required allocation of memory
    INTEGER :: ierr
    !! Error status
    !
    ! frequency-grid for imaginary-axis
    ! nsiw(itemp) = nr. of grid points between (0, wscut) 
    !
    IF (laniso) THEN
      ! memory allocated for wsi and ws
      imelt = nsiw(itemp) + nsw 
      CALL mem_size_eliashberg(2, imelt)
    ENDIF
    !
    ALLOCATE(wsi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('gen_freqgrid_iaxis', 'Error allocating wsi', 1)
    !
    wsi(:) = zero
    DO iw = 1, nsiw(itemp)
      n = iw - 1
      wsi(iw) = DBLE(2 * n + 1) * pi * estemp(itemp) 
      !WRITE(*, *) iw, wsi(iw)
    ENDDO
    !
    ! frequency-grid for real-axis ( Pade approximants and analytic continuation)
    !
    IF (lpade .OR. lacon) THEN
      ALLOCATE(ws(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('gen_freqgrid_iaxis', 'Error allocating ws', 1)
      ws(:) = zero
      DO iw = 1, nsw
        IF (iw <= nqstep) THEN 
          ws(iw) = wsph(iw)
        ELSE
          ws(iw) = wsphmax + DBLE(iw - nqstep) * dwsph
        ENDIF
       !WRITE(*, *) iw, ws(iw), wsph(iw)
      ENDDO
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gen_freqgrid_iaxis
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE gamma_acont(omega, omegap, temp, rgammap, rgammam)
    !-----------------------------------------------------------------------
    !!
    !! computes gammam(w,wp)  (notes RM)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    !
    USE kinds, ONLY : DP
    USE constants_epw, ONLY : eps6, zero, one
    ! 
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency w at point iw on the real-axis
    REAL(KIND = DP), INTENT(in) :: omegap
    !! frequency w' at point iwp on the real-axis
    REAL(KIND = DP), INTENT(in) :: temp
    !! temperature in eV
    REAL(KIND = DP), INTENT(out) :: rgammap
    !! - bose_einstein(w') - fermi_dirac(w + w')
    REAL(KIND = DP), INTENT(out) :: rgammam
    !! bose_einstein(w') + fermi_dirac(-w + w')
    ! 
    rgammap = zero
    rgammam = zero
    IF (ABS(temp) < eps6) THEN
      rgammap = zero
      rgammam = one
    ELSEIF (omegap > zero) THEN 
      rgammap = 0.5d0 * (TANH(0.5d0 * (omega + omegap) / temp) &
                         - 1.d0 / TANH(0.5d0 * omegap / temp))
      rgammam = 0.5d0 * (TANH(0.5d0 * (omega - omegap) / temp) &
                         + 1.d0 / TANH(0.5d0 * omegap / temp))
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gamma_acont
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE dos_quasiparticle(itemp)
    !-----------------------------------------------------------------------
    !!
    !! Computes the quasiparticle density of states in the superconducting state
    !!
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iuqdos
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : lreal, limag, liso, laniso, fsthick
    USE eliashbergcom, ONLY : nsw, estemp, dwsph, ws, dws, delta, adelta, & 
                              wkfs, w0g, nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV, zero, ci
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    CHARACTER(LEN = 256) :: fildos
    !! name dos file

    INTEGER :: iw
    !! Counter over frequency real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: degaussw0
    !! variable to shift the freq. by a small imaginary value
    REAL(KIND = DP) :: weight
    !! quasi-dos factor
    REAL(KIND = DP) :: temp
    !! temperature in K
    REAL(KIND = DP), ALLOCATABLE :: dos_qp(:)
    !! superconducting quasi-particle dos
    !
    COMPLEX(KIND = DP) :: omega
    !! frequency
    !
    degaussw0 = 0.0_DP
    IF (lreal) THEN 
      degaussw0 = 1.d0 * dws(1)
    ELSEIF (limag) THEN 
      degaussw0 = 1.d0 * dwsph
    ENDIF
    !
    temp = estemp(itemp) / kelvin2eV
    IF (temp < 10.d0) THEN
      WRITE(fildos, '(a, a8, f4.2)') TRIM(prefix), '.qdos_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
      WRITE(fildos, '(a, a7, f5.2)') TRIM(prefix), '.qdos_0', temp
    ELSEIF (temp >= 100.d0) THEN
      WRITE(fildos, '(a, a6, f6.2)') TRIM(prefix), '.qdos_', temp
    ENDIF
    OPEN(iuqdos, FILE = fildos, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('dos_quasiparticle', 'error opening file ' // fildos, iuqdos)
    !
    ALLOCATE(dos_qp(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('dos_quasiparticle', 'Error allocating dos_qp', 1)
    dos_qp(:) = zero          
    !
    IF (laniso) THEN
      WRITE(iuqdos, '(5a20)') 'w [eV]', 'N_S/N_F'
      DO iw = 1, nsw 
        omega = ws(iw) + ci * degaussw0
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik)
              dos_qp(iw) = dos_qp(iw) + weight & 
                         * REAL(omega / SQRT(omega * omega - adelta(ibnd, ik, iw) * adelta(ibnd, ik, iw))) 
            ENDIF
          ENDDO
        ENDDO
        WRITE(iuqdos, '(2ES20.10)') ws(iw), dos_qp(iw)
      ENDDO
    ELSEIF (liso) THEN 
      WRITE(iuqdos, '(5a20)') 'w [eV]', 'N_S/N_F'
      DO iw = 1, nsw
        omega = ws(iw) + ci * degaussw0
        dos_qp(iw) = dos_qp(iw) + REAL(omega / SQRT(omega * omega - delta(iw) * delta(iw))) 
        WRITE(iuqdos, '(2ES20.10)') ws(iw), dos_qp(iw)
      ENDDO
    ENDIF
    CLOSE(iuqdos)
    !
    DEALLOCATE(dos_qp, STAT = ierr)
    IF (ierr /= 0) CALL errore('dos_quasiparticle', 'Error deallocating dos_qp', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE dos_quasiparticle
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE free_energy(itemp)
    !-----------------------------------------------------------------------
    !!
    !! Computes the free energy difference between the superconducting and 
    !! normal states
    !!
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufe
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : liso, laniso, fsthick
    USE eliashbergcom, ONLY : estemp, wsi, nsiw, adeltai, aznormi, naznormi, &
                              deltai, znormi, nznormi, &
                              wkfs, w0g, nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : pi, kelvin2eV, zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    ! 
    ! Local variables
    CHARACTER(LEN = 256) :: filfe
    !! name dos file

    INTEGER :: iw
    !! Counter over frequency real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ios
    !! IO error message
    !
    REAL(KIND = DP) :: omega
    !! sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: weight
    !! free energy factor
    REAL(KIND = DP) :: temp
    !! temperature in K
    REAL(KIND = DP) :: dFE
    !! free energy difference between supercond and normal states
    !
    temp = estemp(itemp) / kelvin2eV
    IF (temp < 10.d0) THEN
      WRITE(filfe, '(a, a6, f4.2)') TRIM(prefix), '.fe_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
      WRITE(filfe,'(a, a5, f5.2)') TRIM(prefix), '.fe_0', temp
    ELSEIF (temp >= 100.d0) THEN
      WRITE(filfe,'(a, a4, f6.2)') TRIM(prefix), '.fe_', temp
    ENDIF
    OPEN(iufe, FILE = filfe, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('dos_quasiparticle', 'error opening file ' // filfe, iufe)
    !
    dFE = zero
    IF (laniso) THEN
      DO iw = 1, nsiw(itemp)
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik)
              omega = DSQRT(wsi(iw) * wsi(iw) + adeltai(ibnd, ik, iw) * adeltai(ibnd, ik, iw))
              dFE = dFE - weight * (omega - wsi(iw)) & 
                  * (aznormi(ibnd, ik, iw) - naznormi(ibnd, ik, iw) * wsi(iw) / omega)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSEIF (liso) THEN
      DO iw = 1, nsiw(itemp) 
        omega = DSQRT(wsi(iw) * wsi(iw) + deltai(iw) * deltai(iw))
        dFE = dFE - (omega - wsi(iw)) &
            * (znormi(iw) - nznormi(iw) * wsi(iw) / omega)
      ENDDO
    ENDIF
    dFE = dFE * pi * estemp(itemp)
    WRITE(iufe, '(2ES20.10)') temp, dFE
    CLOSE(iufe)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE free_energy
    !-----------------------------------------------------------------------
    !                                                                            
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_iaxis()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by sum_eliashberg_(an)iso_iaxis
    !!
    !----------------------------------------------------------------------
    !
    USE epwcom, ONLY : liso, laniso
    USE eliashbergcom, ONLY : wsi, deltai, znormi, gap, deltaip, nznormi, keri, &
                              adeltai, adeltaip, aznormi, naznormi, agap
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    DEALLOCATE(wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating wsi', 1)
    DEALLOCATE(deltai, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating deltai', 1)
    DEALLOCATE(znormi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating znormi', 1)
    DEALLOCATE(nznormi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating nznormi', 1)
    DEALLOCATE(gap, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating gap', 1)
    !
    IF (liso) THEN
      DEALLOCATE(deltaip, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating deltaip', 1)
      DEALLOCATE(keri, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating keri', 1)
    ENDIF
    !
    IF (laniso) THEN
      DEALLOCATE(adeltai, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating adeltai', 1)
      DEALLOCATE(adeltaip, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating adeltaip', 1)
      DEALLOCATE(aznormi, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating aznormi', 1)
      DEALLOCATE(naznormi, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating naznormi', 1)
      DEALLOCATE(agap, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_iaxis', 'Error deallocating agap', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_eliashberg_iaxis
    !-----------------------------------------------------------------------
    !                         
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_raxis()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by sum_eliashberg_(an)iso_raxis
    !!
    USE epwcom, ONLY : liso, laniso, lreal, limag, lacon
    USE eliashbergcom, ONLY : ws, delta, znorm, deltap, znormp, &
                              adelta, adeltap, aznorm, aznormp, & 
                              dws, fdwp, bewph, kp, km, gp, gm, & 
                              kp, km, dsumi, zsumi
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    DEALLOCATE(ws, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating ws', 1)
    DEALLOCATE(delta, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating delta', 1)
    DEALLOCATE(znorm, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating znorm', 1)
    !
    IF (liso) THEN 
      DEALLOCATE(deltap, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating deltap', 1)
      DEALLOCATE(znormp, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating znormp', 1)
      !
      IF (lreal) THEN
        DEALLOCATE(dws, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating dws', 1)
        DEALLOCATE(fdwp, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating fdwp', 1)
        DEALLOCATE(bewph, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating bewph', 1)
        DEALLOCATE(kp, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating kp', 1)    
        DEALLOCATE(km, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating km', 1)
      ENDIF
      !
      IF (limag .AND. lacon) THEN
        DEALLOCATE(gp, STAT = ierr)           
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating gp', 1)
        DEALLOCATE(gm, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating gm', 1)
        DEALLOCATE(dsumi, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating dsumi', 1)
        DEALLOCATE(zsumi, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating zsumi', 1)
      ENDIF
    ENDIF
    !
    IF (laniso) THEN
      DEALLOCATE(adelta, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating adelta', 1)
      DEALLOCATE(aznorm, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating aznorm', 1)
      IF (lacon) THEN
        DEALLOCATE(adeltap, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating adeltap', 1)
        DEALLOCATE(aznormp, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_eliashberg_raxis', 'Error deallocating aznormp', 1)
      ENDIF
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_eliashberg_raxis
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_iso()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by eliashberg_init and read_a2f
    !!
    USE epwcom,        ONLY : limag
    USE eliashbergcom, ONLY : a2f_iso, wsph, estemp, nsiw
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    DEALLOCATE(estemp, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_iso', 'Error deallocating estemp', 1)
    DEALLOCATE(wsph, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_iso', 'Error deallocating wsph', 1)
    IF (limag) THEN
      DEALLOCATE(nsiw, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_iso', 'Error deallocating nsiw', 1)
    ENDIF
    DEALLOCATE(a2f_iso, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_iso', 'Error deallocating a2f_iso', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_eliashberg_iso
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_aniso()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by read_frequencies,
    !!  read_eigenvalues, read_kqmap, read_ephmat, eliashberg_init, 
    !!  and evaluate_a2f_lambda subroutines 
    !!
    USE epwcom,        ONLY : limag
    USE elph2,         ONLY : wf, wqf, xqf
    USE eliashbergcom, ONLY : ekfs, xkfs, wkfs, g2, a2f_iso, w0g, & 
                              ixkff, ixkqf, ixqfs, nqfs, memlt_pool, & 
                              wsph, estemp, nsiw
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    DEALLOCATE(estemp, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating estemp', 1)
    DEALLOCATE(wsph, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating wsph', 1)
    IF (limag) THEN
      DEALLOCATE(nsiw, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating nsiw', 1)
    ENDIF
    DEALLOCATE(wf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating wf', 1)
    DEALLOCATE(wqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating wqf', 1)
    DEALLOCATE(xqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating xqf', 1)
    DEALLOCATE(ekfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating ekfs', 1)
    DEALLOCATE(xkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating xkfs', 1)
    DEALLOCATE(wkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating wkfs', 1)
    DEALLOCATE(g2, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating g2', 1)
    DEALLOCATE(a2f_iso, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating a2f_iso', 1)
    DEALLOCATE(w0g, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating w0g', 1)
    DEALLOCATE(ixkff, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating ixkff', 1)
    DEALLOCATE(ixkqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating ixkqf', 1)
    DEALLOCATE(ixqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating ixqfs', 1)
    DEALLOCATE(nqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating nqfs', 1)
    DEALLOCATE(memlt_pool, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating memlt_pool', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_eliashberg_aniso
    !-----------------------------------------------------------------------
    !                           
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_elphon()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by read_frequencies,
    !!  read_eigenvalues, read_kqmap, read_ephmat, and evaluate_a2f_lambda 
    !!
    USE epwcom,        ONLY : limag
    USE elph2,         ONLY : wf, wqf, xqf
    USE eliashbergcom, ONLY : ekfs, xkfs, wkfs, g2, a2f_iso, w0g, &
                              ixkff, ixkqf, ixqfs, nqfs, wsph, memlt_pool
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    DEALLOCATE(wsph, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating wsph', 1)
    DEALLOCATE(wf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating wf', 1)
    DEALLOCATE(wqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating wqf', 1)
    DEALLOCATE(xqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating xqf', 1)
    DEALLOCATE(ekfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ekfs', 1)
    DEALLOCATE(xkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating xkfs', 1)
    DEALLOCATE(wkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating wkfs', 1)
    DEALLOCATE(g2, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating g2', 1)
    DEALLOCATE(a2f_iso, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating a2f_iso', 1)
    DEALLOCATE(w0g, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating w0g', 1)
    DEALLOCATE(ixkff, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ixkff', 1)
    DEALLOCATE(ixkqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ixkqf', 1)
    DEALLOCATE(ixqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ixqfs', 1)
    DEALLOCATE(nqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating nqfs', 1)
    DEALLOCATE(memlt_pool, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating memlt_pool', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_eliashberg_elphon
    !-----------------------------------------------------------------------
  !----------------------------------------------------------------------
  END MODULE supercond
  !----------------------------------------------------------------------
