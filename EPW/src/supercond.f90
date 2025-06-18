  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
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
    !!
    !! This routine initializes the control variables needed to solve the eliashberg
    !! equations
    !!
    USE kinds,            ONLY : DP
    USE mp_global,        ONLY : world_comm
    USE io_global,        ONLY : stdout, ionode_id, ionode
    USE mp_world,         ONLY : mpime
    USE mp,               ONLY : mp_barrier, mp_bcast
    USE input,            ONLY : eliashberg, nkf1, nkf2, nkf3, nsiter, &
                                 nqf1, nqf2, nqf3, nswi, muc, lreal, lpade, &
                                 liso, limag, laniso, lacon, kerwrite, kerread, &
                                 imag_read, fila2f, wscut, rand_q, rand_k, &
                                 ep_coupling, tc_linear, tc_linear_solver, &
                                 gridsamp, fbw, positive_matsu, icoulomb, &
                                 eps_cut_ir, fixsym, emax_coulomb, emin_coulomb
    USE supercond_common, ONLY : size_ir, siz_ir, siz_ir_cl
    USE noncollin_module, ONLY : noncolin
    USE ions_base,        ONLY : amass, ityp, nat, tau
    USE cell_base,        ONLY : at, bg
    USE ep_constants,     ONLY : ryd2ev, zero
    USE global_var,       ONLY : gtemp, elph
    USE io_var,           ONLY : crystal
    USE symm_base,        ONLY : s, set_sym_bl, find_sym
    USE ions_base,        ONLY : nat, tau, ityp, amass
    USE cell_base,        ONLY : alat, omega
    USE modes,            ONLY : nmodes
    USE pwcom,            ONLY : nelec
    USE noncollin_module, ONLY : noncolin, m_loc
    USE low_lvl,          ONLY : fix_sym
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! Contains the state of the opened file
    !
    IF (eliashberg .AND. MINVAL(gtemp(:)) <= 0) &
      CALL errore('eliashberg_init', 'eliashberg requires min temperature > 0', 1)
    IF (eliashberg .AND. liso .AND. laniso) &
      CALL errore('eliashberg_init', 'liso or laniso needs to be true', 1)
    IF (.NOT. eliashberg .AND. liso) &
      CALL errore('eliashberg_init', 'liso requires eliashberg true', 1)
    IF (.NOT. eliashberg .AND. laniso) &
      CALL errore('eliashberg_init', 'laniso requires eliashberg true', 1)
    ! SH: on restarting aniso runs, a2f file can be used for initial gap estimate
    !IF (laniso .AND. (fila2f /= ' ')) &
    !  CALL errore('eliashberg_init', 'anisotropic case can not use fila2f', 1)
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
    IF (eliashberg .AND. nswi > 0 .AND. .NOT. limag) &
      CALL errore('eliashberg_init', 'nswi requires limag true', 1)
    IF (eliashberg .AND. nswi < 0) &
      CALL errore('eliashberg_init', 'nswi should be > 0', 1)
    IF (eliashberg .AND. wscut < 0.d0 ) &
      CALL errore('eliashberg_init', 'wscut should be > 0.d0', 1)
    IF (eliashberg .AND. nsiter < 1) &
      CALL errore('eliashberg_init', 'wrong number of nsiter', 1)
    IF (eliashberg .AND. muc < 0.d0) &
      CALL errore('eliashberg_init', 'muc should be >= 0.d0', 1)
    IF (eliashberg .AND. (rand_k .OR. rand_q) .AND. (fila2f == ' ')) &
      CALL errore('eliashberg_init', 'eliashberg requires a uniform grid when fila2f is not used', 1)
    IF (eliashberg .AND. (MOD(nkf1, nqf1) /= 0 .OR. MOD(nkf2, nqf2) /= 0 .OR. MOD(nkf3, nqf3) /= 0 ) .AND. (fila2f == ' ')) &
      CALL errore('eliashberg_init', &
                  'eliashberg requires nkf1,nkf2,nkf3 to be multiple of nqf1,nqf2,nqf3 when fila2f is not used', 1)
    IF (eliashberg .AND. tc_linear .AND. (lreal .OR. laniso)) &
      CALL errore('eliashberg_init', 'tc_linear cannot be used with lreal or laniso true', 1)
    IF (eliashberg .AND. tc_linear .AND. fbw) &
      CALL errore('eliashberg_init', 'tc_linear cannot be used with fbw true', 1)
    IF (eliashberg .AND. fbw .AND. lacon) &
      CALL errore('eliashberg_init', &
                  'Analytic continuation of Eliashberg equations is not implemented for fbw true', 1)
    IF (eliashberg .AND. (gridsamp < -1 .OR. gridsamp > 3)) &
      CALL errore('eliashberg_init', 'gridsamp must be: -1, 0, 1, 2, or 3', 1)
    IF (eliashberg .AND. tc_linear .AND. .NOT. (tc_linear_solver == 'lapack' .OR. tc_linear_solver == 'power')) &
      CALL errore('eliashberg_init', 'tc_linear_solver must be either lapack or power', 1)
    !
    ! HM added
    IF (laniso .AND. (.NOT. limag)) CALL errore('eliashberg_init', 'limag should be true when laniso == true', 1)
    !
    IF ((.NOT. fbw) .OR. liso) THEN
      IF (gridsamp == 2) THEN
        CALL errore('eliashberg_init', &
                    'The spare-ir sampling is only implemented for solving anisotropic FBW Eliashberg equations.', 1)
      ENDIF
      IF (.NOT. positive_matsu) THEN
        CALL errore('eliashberg_init', &
                   &'Do not use positive_matsu = false except &
                   &when solving the anisotropic FBW Eliashberg equations.', 1)
      ENDIF
    ENDIF
    IF (liso .AND. (gridsamp == 3)) THEN
      CALL errore('eliashberg_init', &
                  'The FFT scheme is only implemented for solving anisotropic Eliashberg equations.', 1)
    ENDIF
    !
    IF ((icoulomb > 0) .AND. (.NOT. (gridsamp == 2))) THEN
      CALL errore('eliashberg_init', 'gridsamp should be 2 true when icoulomb > 0.', 1)
    ENDIF
    IF ((icoulomb > 0) .AND.  (emax_coulomb <= zero)) THEN
      CALL errore('eliashberg_init', 'emax_coulomb should be positive.', 1)
    ENDIF
    IF ((icoulomb > 0) .AND.  (emin_coulomb >= zero)) THEN
      CALL errore('eliashberg_init', 'emin_coulomb should be negative.', 1)
    ENDIF
    !
    IF (laniso .AND. (.NOT. limag)) CALL errore('eliashberg_init', 'limag should be true when laniso == true', 1)
    !
    IF (gridsamp == 2) THEN
      IF (eps_cut_ir >= 1.0d0) THEN
        CALL errore('eliashberg_init', 'eps_cut_ir should be smaller than 1.', 1)
      ELSEIF ((eps_cut_ir < 1.0d0) .AND. (eps_cut_ir > zero)) THEN
        WRITE(stdout, '(5x,a,ES13.3,a)') &
              'eps_cut_ir = ', eps_cut_ir, ' will be used for iterative calculations to solve the Eliashberg equations.'
      ELSE
        WRITE(stdout, '(5x,a)') 'All IR coefficients are used in iterative calculations.'
      ENDIF
    ENDIF
    !
    ! HM: size_ir should be 20 but check the value just in case.
    IF (size_ir < 0) &
      CALL errore('eliashberg_init', 'size_ir should be greater than or equal to 0.', 1)
    !
    IF ((.NOT. positive_matsu) .AND. (gridsamp == 1)) THEN
      CALL errore('eliashberg_init', 'The sparse sampling requires positive_matsu == true', 1)
    ENDIF
    !
    IF(gridsamp == 2) THEN
      ! If sparse-ir sampling is used, determine siz_ir and siz_ir_cl.
      siz_ir = size_ir
      siz_ir_cl = size_ir
    ENDIF
    !
    ! Ryd to eV
    gtemp(:) = gtemp(:) * ryd2ev
    !
    IF (.NOT. elph .AND. .NOT. ep_coupling) THEN
      !
      IF (icoulomb < 1) THEN
        ! We need BZ info to write FS files
        IF (ionode) THEN
          !
          OPEN(UNIT = crystal, FILE = 'crystal.fmt', STATUS = 'old', IOSTAT = ios)
          IF (ios /= 0) CALL errore('eliashberg_init', 'error opening crystal.fmt', crystal)
          READ(crystal, *) nat
          READ(crystal, *) !nmodes
          READ(crystal, *) !nelec
          READ(crystal, *) at
          READ(crystal, *) bg
          READ(crystal, *) !omega
          READ(crystal, *) !alat
          ALLOCATE(tau(3, nat), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_init', 'Error allocating tau', 1)
          READ(crystal, *) tau
          READ(crystal, *) amass
          ALLOCATE(ityp(nat), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_init', 'Error allocating ityp', 1)
          READ(crystal, *) ityp
          READ(crystal, *) noncolin
          READ(crystal, *) !w_centers
          READ(crystal, *) !L
          !
          ! no need further
          DEALLOCATE(tau, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_init', 'Error deallocating tau', 1)
          DEALLOCATE(ityp, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_init', 'Error deallocating ityp', 1)
          CLOSE(crystal)
        ENDIF ! mpime == ionode_id
        CALL mp_bcast(noncolin, ionode_id, world_comm)
        CALL mp_bcast(bg, ionode_id, world_comm)
      ELSE
        IF (ionode) THEN
          !
          OPEN(UNIT = crystal, FILE = 'crystal.fmt', STATUS = 'old', IOSTAT = ios)
          IF (ios /= 0) CALL errore('eliashberg_init', 'error opening crystal.fmt', crystal)
          READ(crystal, *) nat
          READ(crystal, *) nmodes
          READ(crystal, *) nelec
          READ(crystal, *) at
          READ(crystal, *) bg
          READ(crystal, *) omega
          READ(crystal, *) alat
          ALLOCATE(tau(3, nat), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_init', 'Error allocating tau', 1)
          READ(crystal, *) tau
          READ(crystal, *) amass
          ALLOCATE(ityp(nat), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_init', 'Error allocating ityp', 1)
          READ(crystal, *) ityp
          READ(crystal, *) noncolin
          READ(crystal, *) !w_centers
          READ(crystal, *) !L
          ! no need further
          CLOSE(crystal)
        ENDIF ! mpime == ionode_id
        CALL mp_bcast(nat      , ionode_id, world_comm)
        CALL mp_bcast(nmodes   , ionode_id, world_comm)
        CALL mp_bcast(nelec    , ionode_id, world_comm)
        CALL mp_bcast(at       , ionode_id, world_comm)
        CALL mp_bcast(bg       , ionode_id, world_comm)
        CALL mp_bcast(omega    , ionode_id, world_comm)
        CALL mp_bcast(alat     , ionode_id, world_comm)
        IF (mpime /= ionode_id) THEN
          ALLOCATE(tau(3, nat), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_init', 'Error allocating tau', 1)
        ENDIF
        CALL mp_bcast(tau      , ionode_id, world_comm)
        CALL mp_bcast(amass    , ionode_id, world_comm)
        IF (mpime /= ionode_id) THEN
          ALLOCATE(ityp(nat), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_init', 'Error allocating ityp', 1)
        ENDIF
        CALL mp_bcast(ityp     , ionode_id, world_comm)
        CALL mp_bcast(noncolin , ionode_id, world_comm)
        !
        ! Initialize symmetries and create the s matrix
        s(:, :, :) = 0 ! Symmetry in crystal axis with dim: 3,3,48
        CALL set_sym_bl()
        !
        ! Setup crystal symmetry
        IF (.NOT. ALLOCATED(m_loc)) ALLOCATE(m_loc(3, nat))
        CALL find_sym(nat, tau, ityp, .FALSE., m_loc)
        IF (fixsym) CALL fix_sym(.FALSE.)
        !
        DEALLOCATE(tau, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_init', 'Error deallocating tau', 1)
        DEALLOCATE(ityp, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_init', 'Error deallocating ityp', 1)
      ENDIF ! icoulomb
    ENDIF ! .not. elph .and. .not. ep_coupling
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_init
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_grid()
    !-----------------------------------------------------------------------
    !!
    !! This routine initializes the point grids to solve the
    !! eliashberg equations on real and imag axis
    !!
    !! SH: Modified for sparse sampling of Matsubara frequencies (Nov 2021).
    !! RM: Removed variables lunif, nswfc, nswc, wsfc, pwc (Nov 2021).
    !
    USE kinds,           ONLY : DP
    USE io_global,       ONLY : stdout
    USE input,           ONLY : nqstep, nswi, nstemp, lreal, lpade, limag, &
                                lacon, wscut, gridsamp, positive_matsu
    USE global_var,      ONLY : gtemp
    USE ep_constants,    ONLY : eps6
    USE ep_constants,    ONLY : pi
    USE supercond_common,       ONLY : nsw, nsiw, wsphmax, nlambda
    USE io_var,          ONLY : iufilmat
    !
    IMPLICIT NONE
    !
    INTEGER :: itemp
    !! Counter on temperature values
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: lambda_ir
    !! 10 ** nlambda
    REAL(KIND = DP) :: wmax
    !! lambda / beta
    !
    IF (lreal) THEN
      IF (ABS(wscut) < eps6) THEN
        wscut = 10.d0 * wsphmax
      ENDIF
      nsw = nqstep * NINT(wscut / wsphmax)
      IF (nsw == 0) CALL errore('eliashberg_grid', 'wrong number of nsw', 1)
    ENDIF
    !
    IF (limag) THEN
      ! SH: Sparse sampling; gridsamp = -1 is a debug/development option
      !       once set, code reads in the Matsubara indices from file.
      IF (gridsamp == -1) THEN
        nswi = 0
        OPEN(UNIT = iufilmat, FILE = 'matsu-freq.in', STATUS = 'old', IOSTAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_grid', 'Error opening matsu-freq.in', iufilmat)
        nswi = -1
        ierr = 0
        DO WHILE (ierr == 0)
          nswi = nswi + 1
          READ(iufilmat, *, IOSTAT = ierr)
        ENDDO
        CLOSE(iufilmat)
        nsiw(:) = nswi
        IF ((2.d0 * nswi + 1) * pi * gtemp(1) > wscut) THEN
          wscut = (2.d0 * DBLE(nswi) + 1.d0) * pi * gtemp(1)
        ENDIF
      ENDIF
      ! SH: Sparse sampling; with gridsamp = 0 a uniform grid will be generated,
      !       while the gridsamp = 1 generates a grid with points distanced exponentially
      !       from one another (so, sparse sampling!)
      IF ((gridsamp == 0) .OR. (gridsamp == 1) .OR. (gridsamp == 3)) THEN
        ALLOCATE(nsiw(nstemp), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_grid', 'Error allocating nsiw', 1)
        nsiw(:) = 0
        IF (nswi > 0) THEN
          nsiw(:) = nswi
        ELSEIF (wscut > 0.d0) THEN
          DO itemp = 1, nstemp
            nsiw(itemp) = INT(0.5d0 * (wscut / pi / gtemp(itemp) - 1.d0)) + 1
          ENDDO
        ELSEIF (nswi > 0 .AND. wscut > 0.d0) THEN
          nsiw(:) = nswi
          WRITE(stdout,'(5x,a)') 'when nswi > 0, wscut is not used for limag=.TRUE.'
        ENDIF
        IF (ABS(wscut) < eps6) THEN
          wscut = 10.d0 * wsphmax
        ENDIF
        IF (.NOT. positive_matsu) THEN
          nsiw(:) = 2 * nsiw(:)
        ENDIF
      ENDIF
    ENDIF
    !
    IF (gridsamp == 2) THEN
      !
      !! allocate nsiw but not determine the value here; 
      !! nsiw will be determined for each itemp
      !
      ALLOCATE(nsiw(nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_grid', 'Error allocating nsiw', 1)!
      !
      nsiw(:) = 0
      !
      lambda_ir = 1E1_DP ** nlambda
      DO itemp = 1, nstemp
        wmax = lambda_ir * gtemp(itemp)
      ENDDO
    ENDIF
    !
    IF (lpade .OR. lacon) THEN
      nsw = nqstep * NINT(wscut / wsphmax)
      IF (nsw == 0) CALL errore('eliashberg_grid', 'wrong number of nsw', 1)
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
    !!
    !! computes the isotropic spectral function a2F(w), total lambda, and
    !! distribution of lambda
    !!
    !! SH: The "phdos" parts are moved to write_phdos subroutine(Nov 2021).
    !! SM: Updated for distributing freq among images
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iua2ffil, iufillambda, iufillambdaFS
    USE io_files,      ONLY : prefix
    USE modes,         ONLY : nmodes
    USE cell_base,     ONLY : bg
    USE control_flags, ONLY : iverbosity
    USE global_var,    ONLY : wqf, wf, bztoibz
    USE input,         ONLY : fsthick, eps_acoustic, nqstep, degaussq, delta_qsmear, nqsmear, &
                              degaussw, nkf1, nkf2, nkf3
    USE supercond_common,     ONLY : nkfs, nbndfs, g2, ixkqf, ixqfs, nqfs, w0g, ekfs, ef0, dosef, wsph, &
                              wkfs, dwsph, ixkff, ekfs_all, nbndfs_all, ibnd_kfs_all_to_kfs, &
                              ixkf
    USE ep_constants,  ONLY : ryd2ev, eps2, zero, eps16, eps5
    USE io_global,     ONLY : ionode_id, ionode
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE parallelism,   ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: name1
    !! file name
    !
    LOGICAL :: old_version
    !! If true, the format of outputs is compatible with old versions.
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ikfs
    !! Counter on k-points: ikfs = ixkf(bztoibz(ik))
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
    INTEGER :: ibinmin
    !! minimum of ibin
    INTEGER :: ibinmax
    !! maximum of ibin
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
    REAL(KIND = DP):: l_a2f_tmp
    !! Temporary variable for total e-ph coupling strength
    REAL(KIND = DP):: lambda_eph
    !! total e-ph coupling strength (a2f integration)
    REAL(KIND = DP) :: x1, x2, x3
    !! Cartesian coordinates of grid points nkf1, nkf2, nkf3
    REAL(KIND = DP) :: weight, weight2, weight3, weightq
    !! factors in lambda_eph and a2f
    REAL(KIND = DP) :: sigma
    !! smearing in delta function
    REAL(KIND = DP) :: rdum
    !! Dummy for real numbers
    REAL(KIND = DP) :: lambda_step
    !! Current lambda
    REAL(KIND = DP) :: smear, smeark
    !! smearing used for w0gauss
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: lambda_max(npool)
    !! Max e-ph coupling strength among all pools
    REAL(KIND = DP), ALLOCATABLE :: l_a2f(:, :)
    !! total e-ph coupling strength for different ismear (a2f integration)
    REAL(KIND = DP), ALLOCATABLE :: lambda_pairs(:)
    !! Histogram lambda_nk,n'k
    REAL(KIND = DP), ALLOCATABLE :: lambda_k_bin(:)
    !! Histogram lambda_nk
    REAL(KIND = DP), ALLOCATABLE :: a2f(:, :)
    !! Eliashberg sperctral function for different ismear
    REAL(KIND = DP), ALLOCATABLE :: a2f_modeproj(:, :)
    !! Eliashberg sperctral function projected over modes for different ismear
    REAL(KIND = DP), ALLOCATABLE :: lambda_k(:, :)
    !! anisotropic e-ph coupling strength
    !
    old_version = .FALSE.
    !
    ! This is only a quick fix since the routine was written for parallel execution - FG June 2014
#if !defined(__MPI)
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
    ALLOCATE(a2f(nqstep, nqsmear), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating a2f', 1)
    ALLOCATE(a2f_modeproj(nmodes, nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating a2f_modeproj', 1)
    a2f(:, :) = zero
    a2f_modeproj(:, :) = zero
    !
    ! RM - the 0 index in k is required when printing out values of lambda_k
    ! When the k-point is outside the Fermi shell, ixkff(ik)=0
    ALLOCATE(lambda_k(0:nkfs, 0:nbndfs), STAT = ierr)
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
                    IF (wf(imode, iq0) > eps_acoustic) THEN
                      IF (ismear == 1) THEN
                        lambda_eph = lambda_eph + g2(ik, iq, ibnd, jbnd, imode) / wf(imode, iq0)
                      ENDIF
                      weight2 = weight * g2(ik, iq, ibnd, jbnd, imode)
                      DO iwph = 1, nqstep
                        weightq  = w0gauss((wsph(iwph) - wf(imode, iq0)) / sigma, 0) / sigma
                        a2f(iwph, ismear) = a2f(iwph, ismear) + weight2 * weightq
                        IF (ismear == 1) THEN
                          a2f_modeproj(imode, iwph) = a2f_modeproj(imode, iwph) + &
                                       weight2 * weightq
                        ENDIF
                      ENDDO ! iwph
                    ENDIF ! wf
                  ENDDO ! imode
                  IF (ismear == 1 .AND. lambda_eph > 0.d0) THEN
                    l_sum = l_sum + weight * lambda_eph
                    weight3 = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0))
                    lambda_k(ik, ibnd) = lambda_k(ik, ibnd) + weight3 * lambda_eph
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
    IF (ionode) THEN
      !
      !
      name1 = TRIM(prefix) // '.a2f'
      OPEN(UNIT = iua2ffil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iua2ffil)
      !
      ALLOCATE(l_a2f(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating l_a2f', 1)
      l_a2f(:, :) = zero
      !
      DO ismear = 1, nqsmear
        IF (ismear == nqsmear) THEN
          WRITE(iua2ffil, '(" w[meV] a2f and integrated 2*a2f/w for ", i4, " smearing values")') ismear
        ENDIF
        l_a2f_tmp = zero
        DO iwph = 1, nqstep
          l_a2f_tmp = l_a2f_tmp + 2.0d0 * (a2f(iwph, ismear) / wsph(iwph)) * dwsph
          l_a2f(iwph, ismear) = l_a2f_tmp
          ! wsph in meV (from eV)
          IF (ismear == nqsmear) THEN
            WRITE(iua2ffil, '(f12.7, 20f12.7)') wsph(iwph) * 1000.d0, a2f(iwph, :), l_a2f(iwph, :)
          ENDIF
        ENDDO
      ENDDO
      !
      WRITE(iua2ffil, *) "Integrated el-ph coupling"
      WRITE(iua2ffil, '("  #         ", 15f12.7)') l_a2f(nqstep, :)
      WRITE(iua2ffil, *) "Phonon smearing (meV)"
      WRITE(iua2ffil, '("  #         ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * 1000.d0, ismear = 1, nqsmear)
      WRITE(iua2ffil, '("Electron smearing (eV)", f12.7)') degaussw
      WRITE(iua2ffil, '("Fermi window (eV)", f12.7)') fsthick
      WRITE(iua2ffil, '("Summed el-ph coupling ", f12.7)') l_sum
      CLOSE(iua2ffil)
      !
      name1 = TRIM(prefix) // '.a2f_proj'
      OPEN(UNIT = iua2ffil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iua2ffil)
      !
      WRITE(iua2ffil, '("w[meV] a2f a2f_modeproj")')
      DO iwph = 1, nqstep
        ! wsph in meV (from eV)
        WRITE(iua2ffil, '(f12.7, 100f12.7)') wsph(iwph) * 1000.d0, a2f(iwph, 1), a2f_modeproj(:, iwph)
      ENDDO
      WRITE(iua2ffil, '(a, f18.7, a, f18.7)') 'lambda_int = ', l_a2f(nqstep, 1), '   lambda_sum = ',l_sum
      CLOSE(iua2ffil)
      !
      DEALLOCATE(l_a2f, STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating l_a2f', 1)
      !
    ENDIF
    !
    DEALLOCATE(a2f, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating a2f', 1)
    DEALLOCATE(a2f_modeproj, STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error deallocating a2f_modeproj', 1)
    !
    IF (.NOT. old_version) THEN
      nbink = 300
      dbink = 1.1d0 * MAXVAL(lambda_k(:, :)) / DBLE(nbink)
      smeark = 5.0d-1 * dbink
    ELSE
      nbink = NINT(1.1d0 * MAXVAL(lambda_k(:, :)) / eps2) + 1
      dbink = 1.1d0 * MAXVAL(lambda_k(:, :)) / DBLE(nbink)
    ENDIF
    !
    ALLOCATE(lambda_k_bin(nbink), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating lambda_k_bin', 1)
    lambda_k_bin(:) = zero
    !
    nbin = 0
    dbin = zero
    IF (iverbosity == 2) THEN
      IF (.NOT. old_version) THEN
        !nbin = 300
        !dbin = 1.1d0 * MAXVAL(lambda_max(:)) / DBLE(nbin)
        nbin = NINT(1.1d0 * MAXVAL(lambda_max(:)) / eps2) + 1
        dbin = 1.1d0 * MAXVAL(lambda_max(:)) / DBLE(nbin)
        smear = 5.0d-1 * dbin
      ELSE
        nbin = NINT(1.1d0 * MAXVAL(lambda_max(:)) / eps2) + 1
        dbin = 1.1d0 * MAXVAL(lambda_max(:)) / DBLE(nbin)
      ENDIF
      ALLOCATE(lambda_pairs(nbin), STAT = ierr)
      IF (ierr /= 0) CALL errore('evaluate_a2f_lambda', 'Error allocating lambda_pairs', 1)
      lambda_pairs(:) = zero
    ENDIF
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
                  IF (wf(imode, iq0) > eps_acoustic) THEN
                    lambda_eph = lambda_eph + g2(ik, iq, ibnd, jbnd, imode) / wf(imode,iq0)
                  ENDIF
                ENDDO
                lambda_eph = 2.d0 * lambda_eph * dosef
                lambda_k(ik, ibnd) = lambda_k(ik, ibnd) +  weight * lambda_eph
                IF (iverbosity == 2) THEN
                  IF (.NOT. old_version) THEN
                    DO ibin = 1, nbin
                      lambda_step = dbin * DBLE(ibin - 1)
                      weight2 = wkfs(ik) * wqf(iq) * w0gauss((lambda_step - lambda_eph) / smear, 0) / smear
                      weight2 = weight2 * w0g(ibnd, ik) * w0g(jbnd,ixkqf(ik, iq0))
                      lambda_pairs(ibin) = lambda_pairs(ibin) + weight2
                    ENDDO
                  ELSE
                    ibin = NINT(lambda_eph / dbin) + 1
                    weight2 =  wkfs(ik) * wqf(iq) * w0g(ibnd, ik) * w0g(jbnd,ixkqf(ik, iq0))
                    lambda_pairs(ibin) = lambda_pairs(ibin) + weight2
                  ENDIF
                ENDIF
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
          IF (.NOT. old_version) THEN
            DO ibin = 1, nbink
              lambda_step = dbink * DBLE(ibin - 1)
              weight3 = wkfs(ik) * w0gauss((lambda_step - lambda_k(ik, ibnd)) / smeark, 0) / smeark * w0g(ibnd, ik)
              lambda_k_bin(ibin) = lambda_k_bin(ibin) + weight3
            ENDDO
          ELSE
            ibin = NINT(lambda_k(ik, ibnd) / dbink) + 1
            weight3 = w0g(ibnd, ik)
            lambda_k_bin(ibin) = lambda_k_bin(ibin) + weight3
          ENDIF
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
    IF (ionode) THEN
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
      IF (.NOT. old_version) THEN
        WRITE(iufillambda, '(a32)') '# distribution = \rho(lambda_nk)'
        WRITE(iufillambda, '(a21,2(2x,a20))') '#        lambda_nk   ', ' dist. (scaled to 1)', '  dist. (not scaled)'
        DO ibin = nbink, 1, -1
          IF ((lambda_k_bin(ibin) / MAXVAL(lambda_k_bin(:))) > eps5) ibinmin = ibin
        ENDDO
        ibinmin = ibinmin - 1
        IF (ibinmin < 1) ibinmin = 1
        DO ibin = 1, nbink
          IF ((lambda_k_bin(ibin) / MAXVAL(lambda_k_bin(:))) > eps5) ibinmax = ibin
        ENDDO
        ibinmax = ibinmax + 1
        IF (ibinmax > nbink) ibinmax = nbink
        DO ibin = ibinmin, ibinmax
          WRITE(iufillambda,'(1x, ES20.10, 2(2x,ES20.10))') dbink * DBLE(ibin - 1), &
                                  lambda_k_bin(ibin) / MAXVAL(lambda_k_bin(:)), lambda_k_bin(ibin)
        ENDDO
      ELSE
        WRITE(iufillambda, '(a12, a30)') '# lambda_nk','  \rho(lambda_nk) scaled to 1'
        DO ibin = 1, nbink
          WRITE(iufillambda,'(2f21.7)') dbink * DBLE(ibin), lambda_k_bin(ibin) / MAXVAL(lambda_k_bin(:))
        ENDDO
      ENDIF
      CLOSE(iufillambda)
      !
      ! SP: Produced if user really wants it
      IF (iverbosity == 2) THEN
        name1 = TRIM(prefix) // '.lambda_pairs'
        OPEN(UNIT = iufillambda, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iufillambda)
        IF (.NOT. old_version) THEN
          WRITE(iufillambda, '(a37)') "# distribution = \rho(lambda_nk,n'k')"
          WRITE(iufillambda, '(a21,2(2x,a20))') "#    lambda_nk,n'k'  ", ' dist. (scaled to 1)', '  dist. (not scaled)'
          DO ibin = nbin, 1, -1
            IF ((lambda_pairs(ibin) / MAXVAL(lambda_pairs(:))) > eps5) ibinmin = ibin
          ENDDO
          ibinmin = ibinmin - 1
          IF (ibinmin < 1) ibinmin = 1
          DO ibin = 1, nbin
            IF ((lambda_pairs(ibin) / MAXVAL(lambda_pairs(:))) > eps5) ibinmax = ibin
          ENDDO
          ibinmax = ibinmax + 1
          IF (ibinmax > nbin) ibinmax = nbin
          DO ibin = ibinmin, ibinmax
            WRITE(iufillambda,'(1x, ES20.10, 2(2x,ES20.10))') dbin * DBLE(ibin - 1), &
                                    lambda_pairs(ibin) / MAXVAL(lambda_pairs(:)), lambda_pairs(ibin)
          ENDDO
        ELSE
          WRITE(iufillambda, '(a12, a30)') "# lambda_nk,n'k'", "  \rho(lambda_nk,n'k') scaled to 1"
          DO ibin = 1, nbin
            WRITE(iufillambda, '(2f21.7)') dbin * DBLE(ibin), lambda_pairs(ibin) / MAXVAL(lambda_pairs(:))
          ENDDO
        ENDIF
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
        DO ibnd = 1, nbndfs_all
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
          !WRITE(iufillambdaFS, '(6f12.6)') (lambda_k(ixkff(ik), ibnd), ik = 1, nkf1 * nkf2 * nkf3)
          i = 1
          DO ik = 1, nkf1 * nkf2 * nkf3
            ikfs = ixkf(bztoibz(ik))
            IF (ikfs == 0) THEN
              rdum = 0.0d0
            ELSE
              rdum = lambda_k(ikfs, ibnd_kfs_all_to_kfs(ibnd, ikfs))
            ENDIF
            IF (i == 6) THEN
              WRITE(iufillambdaFS, '(f12.6)') rdum
            ELSE
              WRITE(iufillambdaFS, '(f12.6)', advance='no') rdum
            ENDIF
            i = i + 1
            IF (i == 7) i = 1
          ENDDO
          CLOSE(iufillambdaFS)
        ENDDO
        ! HP: Write in .frmsf format compatible with fermisurfer program
        WRITE(name1, '(a, a13)') TRIM(prefix), '.lambda.frmsf'
        OPEN(UNIT = iufillambdaFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('evaluate_a2f_lambda', 'error opening file ' // name1, iufillambdaFS)
        !
        WRITE(iufillambdaFS, '(3i5)') nkf1, nkf2, nkf3
        WRITE(iufillambdaFS, '(i5)') 1
        WRITE(iufillambdaFS, '(i5)') nbndfs_all
        WRITE(iufillambdaFS, '(3f12.6)') (bg(i, 1), i = 1, 3)
        WRITE(iufillambdaFS, '(3f12.6)') (bg(i, 2), i = 1, 3)
        WRITE(iufillambdaFS, '(3f12.6)') (bg(i, 3), i = 1, 3)
        ! HM: Outputting a dummy value in the .frmsf file may form fake Fermi surfaces.
        !     To avoid using a dummy value for the states outside of fsthick window, 
        !     use ekfs_all instead of ekfs.
        !WRITE(iufillambdaFS, '(6f12.6)') ((ekfs(ibnd, ixkff(ik)) - ef0, ik = 1, nkf1 * nkf2 * nkf3), ibnd = 1, nbndfs)
        !WRITE(iufillambdaFS, '(6f12.6)') ((lambda_k(ixkff(ik), ibnd), ik = 1, nkf1 * nkf2 * nkf3), ibnd = 1, nbndfs)
        i = 1
        DO ibnd = 1, nbndfs_all
          DO ik = 1, nkf1 * nkf2 * nkf3
            IF (i == 6) THEN
              WRITE(iufillambdaFS, '(f12.6)') ekfs_all(ibnd, bztoibz(ik)) - ef0
            ELSE
              WRITE(iufillambdaFS, '(f12.6)', advance='no') ekfs_all(ibnd, bztoibz(ik)) - ef0
            ENDIF
            i = i + 1
            IF (i == 7) i = 1
          ENDDO
        ENDDO
        i = 1
        DO ibnd = 1, nbndfs_all
          DO ik = 1, nkf1 * nkf2 * nkf3
            ikfs = ixkf(bztoibz(ik))
            IF (ikfs == 0) THEN
              rdum = 0.0d0
            ELSE
              rdum = lambda_k(ikfs, ibnd_kfs_all_to_kfs(ibnd, ikfs))
            ENDIF
            IF (i == 6) THEN
              WRITE(iufillambdaFS, '(f12.6)') rdum
            ELSE
              WRITE(iufillambdaFS, '(f12.6)', advance='no') rdum
            ENDIF
            i = i + 1
            IF (i == 7) i = 1
          ENDDO
        ENDDO
        CLOSE(iufillambdaFS)
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
      !
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
    !!
    !! This routine estimates the Tc using Allen-Dynes formula and
    !! the BCS superconducting gap as the initial guess for delta
    !!
    !! SH: Updated to write the machine learning estimate for Tc, and
    !!       to limit the "gap0" to 6 decimal digits (Nov 2021).
    !!
    !! SM: update to write Allen and Dynes modified McMillan’s semiempirical result, 
    !! and is needed when lambda > 1.5  [Alen and Dyne, PRB 12, 3, 1975]
    !
    USE kinds,            ONLY : DP
    USE input,            ONLY : nqstep, muc, nstemp, icoulomb
    USE global_var,       ONLY : gtemp
    USE supercond_common,        ONLY : wsph, dwsph, a2f_tmp, gap0, spin_fac
    USE ep_constants,     ONLY : kelvin2eV, zero
    USE noncollin_module, ONLY : noncolin
    USE io_global,        ONLY : stdout, ionode_id, ionode
    USE mp_global,        ONLY : inter_pool_comm
    USE mp,               ONLY : mp_bcast, mp_barrier, mp_sum
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
    REAL(KIND = DP):: tcml
    !! Estimated ML Tc
    REAL(KIND = DP):: tc_AD
    !! Estimated AD-modified Tc
    REAL(KIND = DP):: momavg
    !! avg. of 2nd moment of frequencies
    REAL(KIND = DP):: fomega, fmu
    !! coefficients of ML-based Tc
    REAL(KIND = DP):: f1, f2
    !! coefficients of AD Tc
    REAL(KIND = DP):: lambda1, lambda2
    !! Lambda1, Lambda2, ratio in PRB 12, 3, 1975
    REAL(KIND = DP):: muc_local
    !! alternative value for muc: 0.15
    !
    IF (ionode) THEN
      l_a2f  = zero
      logavg = zero
      momavg = zero
      DO iwph = 1, nqstep
        l_a2f  = l_a2f  + a2f_tmp(iwph) / wsph(iwph)
        logavg = logavg + a2f_tmp(iwph) * log(wsph(iwph)) / wsph(iwph)
        momavg = momavg + a2f_tmp(iwph) * wsph(iwph) * dwsph
      ENDDO
      l_a2f  = l_a2f  * 2.d0 * dwsph
      logavg = logavg * 2.d0 * dwsph
      logavg = EXP(logavg / l_a2f)
      momavg = DSQRT(2.d0 * momavg / l_a2f)
      WRITE(stdout,'(5x,a,f12.7)') 'Electron-phonon coupling strength = ', l_a2f
      WRITE(stdout,'(a)') ' '
      !
      IF ((icoulomb > 0) .AND. (muc > 0.15d0)) THEN
        ! HM: In the case of icoulomb > 0, because we use a value relatively larger than mu^*, 
        !     gap0 is estimated as too small if the input muc value is used in the estimations. 
        !     To avoid that, we will use muc = 0.15 in the estimations.
        muc_local = 0.15d0
        WRITE(stdout, '(5x, a/)') 'muc = 0.15 is used in the following estimations'
      ELSE 
        muc_local = muc
      ENDIF
      !
      ! Standard Dyne modified McMillan formula
      ! [W. McMillan, Phys. Rev. 167, 331 (1968), R. C. Dynes, Solid State Commun. 10, 615 (1972).]
      tc = logavg / 1.2d0 * EXP(-1.04d0 * (1.d0 + l_a2f) &
                                  / (l_a2f - muc_local * (1.d0 + 0.62d0 * l_a2f)))
      !
      ! SM: Allen-Dynes modified McMillan formula with strong-coupling corrections 
      ! [Eqs. (34-38) of Allen and Dyne, PRB 12,3 (1975)]
      !
      lambda1 = 2.46d0 * (1.d0 + 3.8d0 * muc_local)
      lambda2 = 1.82d0 * (1.d0 + 6.3d0 * muc_local) * (momavg / logavg)
      ! Strong-coupling correction factor f1
      f1 = (1.d0 + (l_a2f / lambda1)**(3.d0 / 2.d0))**(1.d0 / 3.d0)
      ! Shape correction factor f2
      f2 = 1.d0 + ((momavg / logavg- 1.d0) * l_a2f**2 / (l_a2f**2 + lambda2**2))
      ! AD-Modified Tc with f1 and f2 factors
      tc_AD = (f1 * f2 * logavg / 1.2d0) * EXP(-1.04d0 * (1.d0 + l_a2f) &
                                / (l_a2f - muc_local * (1.d0 + 0.62d0 * l_a2f)))
      !
      ! SH: ML-based estimate of the Tc (Eqns. [6-8]; Xie et. al.; arXiv:2106.05235)
      fomega = 1.92d0 * ((l_a2f + (logavg / momavg) - muc_local**(1.d0 / 3.d0)) / &
                (DSQRT(l_a2f) * EXP(logavg / momavg))) - 0.08d0
      fmu    = (6.86d0 * EXP(-l_a2f / muc_local)) / ((1.d0 / l_a2f) - muc_local - &
                (logavg / momavg)) + 1.d0
      tcml   = (fomega * fmu * logavg / 1.20d0) * EXP(- (1.04d0 * (1 + l_a2f)) / &
                (l_a2f - muc_local * (1 + 0.62d0 * l_a2f)))
      !
      ! initial guess for the gap edge using BCS superconducting ratio 3.52
      !
      ! SH: the "estimated gap" is restricted to 6 decimal places (in meV units)
      !       for consistency with the gap_from_a2f function. That's, by applying
      !       NINT(gap[eV] * 1.d9) / 1.d9
      gap0 = NINT(3.52d0 * tc / 2.d0 * 1.d9) / 1.d9
      IF (gap0 <= 0.d0) &
        CALL errore('estimate_tc_gap', 'initial guess for gap edge should be > 0.d0', 1)
      !
      ! tc in K
      tc    = tc    / kelvin2eV
      tc_AD = tc_AD / kelvin2eV
      tcml  = tcml / kelvin2eV
      !
      WRITE(stdout, '(5x, a, f8.4, a, f8.4)') 'Estimated Tc using McMillan expression = ', tc, ' K for muc = ', muc
      WRITE(stdout, '(a)') '  '
      WRITE(stdout, '(5x, a, f8.4, a)') 'Estimated Tc using Allen-Dynes modified McMillan expression = ', tc_AD, ' K'
      WRITE(stdout, '(a)') '  '
      WRITE(stdout, '(5x, a, f8.4, a)') 'Estimated Tc using SISSO machine learning model = ', tcml , ' K'
      WRITE(stdout, '(a)') '  '
       WRITE(stdout, '(5x, a, f8.4, a)') 'Estimated w_log = ', logavg * 1000.d0, ' meV'
      WRITE(stdout, '(a)') '  '
      WRITE(stdout, '(5x, a, f8.4, a)') 'Estimated BCS superconducting gap using McMillan Tc = ', gap0 * 1000.d0, ' meV'
      WRITE(stdout, '(a)') '  '
      !
      IF (gtemp(1) / kelvin2eV > tc) THEN
        WRITE(stdout, '(a)') '  '
        WRITE(stdout, '(5x, a)') 'WARNING WARNING WARNING '
        WRITE(stdout, '(a)') '  '
        WRITE(stdout, '(5x, a, f9.3, a, f9.3, a)') 'The code may crash since tempsmin =', &
                        gtemp(1) / kelvin2eV, ' K is larger than McMillan Tc = ', tc, ' K'
      ELSEIF (gtemp(nstemp) / kelvin2eV > tc) THEN
        WRITE(stdout, '(a)') '  '
        WRITE(stdout, '(5x, a)') 'WARNING WARNING WARNING '
        WRITE(stdout, '(a)') '  '
        WRITE(stdout, '(5x, a, f9.3, a, f9.3, a)') 'The code may crash since tempsmax =', &
                        gtemp(nstemp) / kelvin2eV, ' K is larger than McMillan Tc = ', tc, ' K'
      ENDIF
      !
      ! HP: muc needs to be divided by 2 to make consistent with the el-ph contribution
      ! in the SOC calculation
      ! HM: Dividing muc by two is necessary only in the anisotropic calculations. 
      ! Directly substituting half the value into muc poses a significant risk.
      ! spin_fac will be multiplied by the Coulomb term only in the anisotropic calculations.
      IF (noncolin) THEN
        spin_fac = 0.5d0
      ELSE
        spin_fac = 1.0d0
      ENDIF
      !
    ENDIF
    CALL mp_bcast(muc, ionode_id, inter_pool_comm)
    CALL mp_bcast(gap0, ionode_id, inter_pool_comm)
    CALL mp_bcast(spin_fac, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE estimate_tc_gap
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE find_a2f()
    !-----------------------------------------------------------------------
    !!
    !! HP: Read the eliashberg spectral function from a2f file if present,
    !!     otherwise evaluate it and estimate the Tc and initial gap value.
    !!     This is to save time once starting the iso/aniso runs (3/2022).
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_files,      ONLY : prefix
    USE input,         ONLY : fila2f, degaussq, delta_qsmear
    USE supercond_common,     ONLY : wsphmax, wsph
    USE ep_constants,  ONLY : ryd2ev
    !
    IMPLICIT NONE
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: ierr
    !! Error status
    !
    IF (fila2f == ' ') WRITE(fila2f, '(a, a4)') TRIM(prefix), '.a2f'
    INQUIRE(FILE = fila2f, EXIST = exst)
    IF (exst) THEN
      WRITE(stdout, '(5x,a)') 'a2f file is found and will be used to estimate initial gap'
      WRITE(stdout, '(a)')    ' '
      !
      ! This is needed because we do not call evaluate_a2f_lambda
      degaussq = degaussq * ryd2ev
      delta_qsmear = delta_qsmear * ryd2ev
    ELSE
      WRITE(stdout, '(5x,a)') &
        'a2f file is not found to estimate initial gap: calculating a2f files'
      WRITE(stdout, '(a)')    ' '
      CALL evaluate_a2f_lambda
    ENDIF
    !
    ! HP: this is requred to make consistient in iso/aniso calculation
    IF (wsphmax > 0.d0) THEN
      ! read_frequencies
      DEALLOCATE(wsph, STAT = ierr)
      IF (ierr /= 0) CALL errore('find_a2f', 'Error deallocating wsph', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE find_a2f
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gen_freqgrid_raxis()
    !-----------------------------------------------------------------------
    !!
    !! Automatic generation of the frequency-grid for real-axis calculations.
    !!
    !
    USE supercond_common,     ONLY : nsw, ws, dwsph
    USE ep_constants,  ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER :: iw
    !! Counter over frequency
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(ws(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('gen_freqgrid_raxis', 'Error allocating ws', 1)
    ws(:) = zero
    !
    DO iw = 1, nsw
      ws(iw) = DBLE(iw) * dwsph
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gen_freqgrid_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gen_freqgrid_iaxis(itemp, ir_obj)
    !-----------------------------------------------------------------------
    !!
    !! Automatic generation of the frequency-grid for imaginary-axis calculations.
    !!
    !! input
    !!
    !! itemp  - temperature point
    !!
    !! SH: Modified for sparse sampling of Matsubara frequencies (Nov 2021).
    !!
    !! =====================================================================
    !! SH: A note about definition of Matsubara indices/frequencies in epw
    !! RM: updated (Jan 2022)
    !!
    !! In epw, the nsiw(itemp) is the cutoff for Matsubara indicies; i.e.,
    !!   the largest positive Matsubara index "n" is nsiw(itemp)-1.
    !!
    !! nsiw(itemp) = INT(0.5d0 * (wscut / pi / gtemp(itemp) - 1.d0)) + 1
    !!
    !! So, for N=nsiw(itemp), indices are:
    !!
    !!     n=    -N, -(N-1), ..., -1, 0, 1, ...,  N-2,  N-1
    !!   and corresponding 2n+1 factor for frequencies are:
    !!     f= -2N+1,  -2N+3, ..., -1, 1, 3, ..., 2N-3, 2N-1
    !!
    !! The actual frequencies are non-symmetric with respect
    !!   to zero, i.e., the (2N-1)*pi*T has no negative counterpart, and
    !!   the rest are (in terms of 2n+1 factor):
    !!
    !!   F[-N] = -F[N-1]; ...; F[-1] = -F[0]
    !!
    !! With using the lambda_negative, the summations in Eliashberg eqns
    !!   in the epw run only over non-negative "n" indices. That's:
    !!
    !!   iw = 1, 2, ..., N
    !!   n  = 0, 1, ..., N-1
    !!   f  = 1, 3, ..., 2N-1
    !!
    !! If fbw = true and positive_matsu = false, the lambda_negative is not used,
    !! the summations in Eliashberg eqns in the epw run over all "n" indices
    !!
    !!   iw =     1,     2,     3, ..., N-1,  N, N+1       2N-1,   2N
    !!   n  =    -N,  -N+1,  -N+2, ...,  -2, -1,   0, ...,  N-2,  N-1
    !!   f  = -2N+1, -2N+3, -2N+5, ...,  -3, -1,   1, ..., 2N-3, 2N-1
    !! NOTE: nsiw(itemp) = 2N after initializing nsiw in eliashberg_grid
    !! =====================================================================
    !!
    !
    USE input,         ONLY : lpade, lacon, laniso, gridsamp, griddens, tc_linear, &
                              positive_matsu
    USE global_var,    ONLY : gtemp
    USE supercond_common,     ONLY : nsw, nsiw, ws, wsi, dwsph, wsn
    USE ep_constants,  ONLY : zero
    USE ep_constants,  ONLY : pi
    USE low_lvl,       ONLY : mem_size_eliashberg
    USE io_var,        ONLY : iufilmat
    USE mp,            ONLY : mp_bcast, mp_barrier
    USE mp_global,     ONLY : inter_pool_comm, inter_image_comm
    USE io_global,     ONLY : stdout, ionode_id, ionode
    USE control_flags, ONLY : iverbosity
    USE mp_world,      ONLY : mpime
    USE sparse_ir,     ONLY : IR
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    TYPE(IR), INTENT(IN), OPTIONAL :: ir_obj
    !! which contains ir objects such as basis functions
    !
    ! Local variables
    CHARACTER(LEN = 13) :: schm
    !! Type of sampling
    INTEGER :: iw
    !! Counter over frequency
    INTEGER :: n
    !! Matsubara frequency index
    INTEGER :: nsiw_half
    !! INT(0.5d0 * (wscut / pi / gtemp(itemp) - 1.d0)) + 1
    INTEGER(8) :: imelt
    !! Required allocation of memory
    INTEGER :: ierr
    !! Error status
    !
    IF (laniso) THEN
      ! memory allocated for wsi, wsn, and ws
      imelt = 2 * nsiw(itemp) + nsw
      CALL mem_size_eliashberg(2, imelt)
    ENDIF
    !
    IF (gridsamp == 2) nsiw(itemp) = ir_obj%nfreq_f
    !
    ALLOCATE(wsi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('gen_freqgrid_iaxis', 'Error allocating wsi', 1)
    !
    wsi(:) = zero
    !
    ! SH: Sparse sampling:
    !       depending on gridsamp value generates a grid of Matsubara indices:
    !       (-1) read from file, (0) uniform grid, (1) sparse sampling.
    !     Note:
    !       - The wsn(nsiw(itemp)) array is introduces to keep the "indices"
    !       - The griddens (default=1.d0) controls grid density; larger values
    !         give denser mesh, and smaller than 1.d0 make the mesh more sparse
    !
    ALLOCATE(wsn(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('gen_freqgrid_iaxis', 'Error allocating wsn', 1)
    !
    IF (ionode) THEN
      !
      ! Initiate sampling scheme names
      IF (gridsamp == -1) THEN
        schm = 'input        '
      ELSEIF (gridsamp == 0) THEN
        schm = 'uniform      '
      ELSEIF (gridsamp == 1) THEN
        schm = 'sparse       '
      ELSEIF (gridsamp == 2) THEN
        schm = 'sparse-ir    '
      ELSEIF (gridsamp == 3) THEN
        schm = 'uniform (FFT)'
      END IF
      !
      ! input sampling
      IF (gridsamp == -1) THEN
        OPEN(UNIT = iufilmat, FILE = 'matsu-freq.in', STATUS = 'old', IOSTAT = ierr)
        IF (ierr /= 0) CALL errore('gen_freqgrid_iaxis', 'Error opening matsu-freq.in', iufilmat)
        DO iw = 1, nsiw(itemp)
          READ(iufilmat,*) n
          wsn(iw) = n
          wsi(iw) = DBLE(2 * n + 1) * pi * gtemp(itemp)
        ENDDO
        CLOSE(iufilmat)
      END IF
      ! uniform sampling
      IF (((gridsamp == 0) .OR. (gridsamp == 3)) .AND. positive_matsu) THEN
        DO iw = 1, nsiw(itemp)
          n = iw - 1
          wsn(iw) = n
          wsi(iw) = DBLE(2 * n + 1) * pi * gtemp(itemp)
        END DO
      END IF
      ! uniform sampling to consider negative freq.
      IF (((gridsamp == 0) .OR. (gridsamp == 3)) .AND. (.NOT.positive_matsu)) THEN 
        nsiw_half = nsiw(itemp) / 2
        DO iw = 1, nsiw(itemp)
          n = iw - nsiw_half - 1
          wsn(iw) = n
          wsi(iw) = DBLE(2 * n + 1) * pi * gtemp(itemp)
        END DO
      END IF
      ! sparse sampling
      IF ((gridsamp == 1) .AND. positive_matsu) THEN
        n  = 0
        iw = 0
        DO WHILE (n < nsiw(itemp))
          iw      = iw + 1
          wsn(iw) = n
          wsi(iw) = DBLE(2 * n + 1) * pi * gtemp(itemp)
          n = n + NINT(EXP(DBLE(n) / DBLE(nsiw(itemp)) / griddens))
        ENDDO
        ! update the number of freqs. to the true value
        nsiw(itemp) = iw
      END IF
      !
      ! sparse-ir sampling
      IF (gridsamp == 2) THEN
        !
        IF (.NOT. PRESENT(ir_obj)) THEN
          CALL errore('gen_freqgrid_iaxis', 'Error: ir_obj is not given despite gridsamp == 2', 1)
        ENDIF
        !
        DO iw = 1, nsiw(itemp)
          wsi(iw) = DBLE(ir_obj%freq_f(iw)) * pi * gtemp(itemp)
          wsn(iw) = (ir_obj%freq_f(iw) - 1) / 2
        ENDDO
      ENDIF
      !
      ! output the indices to "matsu-freq*.out" file, if iverbosity = 2
      IF (iverbosity == 2) CALL write_matsubara_freq(itemp)
      !
      ! print actual number of Matsubara frequencies
      IF (.NOT. tc_linear) &
      WRITE(stdout, '(5x, a, i6, a, i6, a, a, a)') 'Actual number of frequency points (', &
        itemp, ') = ', nsiw(itemp), ' for ', schm, ' sampling'
      !
    ENDIF
    ! this is important! All cores should be aware of the grid.
    CALL mp_bcast(nsiw(itemp), ionode_id, inter_pool_comm)
    CALL mp_bcast(wsi, ionode_id, inter_pool_comm)
    CALL mp_bcast(wsn, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    ! SM: For ditibuting freq. among images, we need to broadcast among images
    CALL mp_bcast(nsiw(itemp), ionode_id, inter_image_comm)
    CALL mp_bcast(wsi, ionode_id, inter_image_comm)
    CALL mp_bcast(wsn, ionode_id, inter_image_comm)
    CALL mp_barrier(inter_image_comm)
    !
    !
    ! frequency-grid for real-axis ( Pade approximants and analytic continuation)
    !
    IF (lpade .OR. lacon) THEN
      ALLOCATE(ws(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('gen_freqgrid_iaxis', 'Error allocating ws', 1)
      ws(:) = zero
      DO iw = 1, nsw
        ws(iw) = DBLE(iw) * dwsph
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
    ! Local variables
    !
    REAL(KIND = DP) :: inv_temp
    !! Invese temperature inv_etemp = 1/temp. Defined for efficiency reason
    REAL(KIND = DP) :: coth
    !! coth = 1/tanh
    !
    inv_temp = 1.d0 / temp
    !
    coth = 1.d0 / TANH(0.5d0 * omegap * inv_temp)
    !
    rgammap = 0.5d0 * (TANH(0.5d0 * (omega + omegap) * inv_temp) - coth)
    rgammam = 0.5d0 * (TANH(0.5d0 * (omega - omegap) * inv_temp) + coth)
    !
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
    USE input,         ONLY : liso, laniso, fsthick
    USE global_var,    ONLY : gtemp
    USE supercond_common,     ONLY : nsw, dwsph, ws, delta, adelta, &
                              wkfs, w0g, nkfs, nbndfs, ef0, ekfs
    USE ep_constants,  ONLY : kelvin2eV, zero, ci
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_world,      ONLY : mpime
    USE parallelism,   ONLY : fkbounds
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
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
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
    degaussw0 = dwsph
    !
    temp = gtemp(itemp) / kelvin2eV
    IF (temp < 10.d0) THEN
      WRITE(fildos, '(a, a8, f4.2)') TRIM(prefix), '.qdos_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
      WRITE(fildos, '(a, a7, f5.2)') TRIM(prefix), '.qdos_0', temp
    ELSEIF (temp >= 100.d0) THEN
      WRITE(fildos, '(a, a6, f6.2)') TRIM(prefix), '.qdos_', temp
    ENDIF
    !
    ALLOCATE(dos_qp(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('dos_quasiparticle', 'Error allocating dos_qp', 1)
    dos_qp(:) = zero
    !
    IF (laniso) THEN
      !
      CALL fkbounds(nkfs, lower_bnd, upper_bnd)
      !
      DO ik = lower_bnd, upper_bnd
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik)
            DO iw = 1, nsw
              omega = ws(iw) + ci * degaussw0
              dos_qp(iw) = dos_qp(iw) + weight &
                         * REAL(omega / SQRT(omega * omega - adelta(iw, ibnd, ik) * adelta(iw, ibnd, ik)))
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      ! collect contributions from all pools
      CALL mp_sum(dos_qp, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
    ELSEIF (liso) THEN
      DO iw = 1, nsw
        omega = ws(iw) + ci * degaussw0
        dos_qp(iw) = dos_qp(iw) + REAL(omega / SQRT(omega * omega - delta(iw) * delta(iw)))
      ENDDO
    ENDIF
    !
    IF (mpime == ionode_id) THEN
      OPEN(iuqdos, FILE = fildos, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('dos_quasiparticle', 'error opening file ' // fildos, iuqdos)
      WRITE(iuqdos, '(5a20)') 'w [eV]', 'N_S/N_F'
      DO iw = 1, nsw
        WRITE(iuqdos, '(2ES20.10)') ws(iw), dos_qp(iw)
      ENDDO
      CLOSE(iuqdos)
    ENDIF
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
    SUBROUTINE crit_temp_solver(adim, a, eigen, niter)
    !-----------------------------------------------------------------------
    !!
    !! SH: Routine for returning largest eigenvalue of a(adim, adim);
    !!       being used to solve linearized Eliashberg equation
    !! HP: updated 5/11/2021
    !
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : zero, two
    USE input,         ONLY : tc_linear_solver, nsiter, conv_thr_iaxis
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: adim
    !! Dimension of input matrix
    INTEGER, INTENT(out) :: niter
    !! Number of iterations of solver
    REAL(KIND = DP), INTENT(in) :: a(adim, adim)
    !! Matrix for solve
    REAL(KIND = DP), INTENT(out) :: eigen
    !! Maximum eigenvalue
    !
    ! Local variables
    CHARACTER(LEN = 10) :: jobvl, jobvr
    !! Eigenvalue problem-related variables
    LOGICAL :: conv
    !! True if calculation is converged
    INTEGER :: n
    !! Dimension of matrix a
    INTEGER :: ierr
    !! Error status
    INTEGER :: iter
    !! Counter on iteration steps
    INTEGER :: lda, ldvl, ldvr, lwork, ix, iy, ic
    !! Eigenvalue problem-related variables
    REAL(KIND = DP) :: alpha, norm
    !! Eigenvalue problem-related variables
    REAL(KIND = DP), ALLOCATABLE :: work(:)
    !! Eigenvalue problem-related variables
    REAL(KIND = DP), ALLOCATABLE :: x(:), y(:)
    !! Eigenvalue problem-related variables
    REAL(KIND = DP), ALLOCATABLE :: wr(:), wi(:), vl(:, :), vr(:, :)
    !! Eigenvalue problem-related variables
    !
    IF (tc_linear_solver == 'lapack') THEN
      n        = adim
      lwork    = 4 * n
      lda      = n
      ldvl     = n
      ldvr     = n
      jobvl    = 'n'
      jobvr    = 'n'
      !
      ! S.tiwari: some allocations to cure ifort overflow
      ALLOCATE(vl(n, n), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error allocating vl', 1)
      ALLOCATE(vr(n, n), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error allocating vr', 1)
      ALLOCATE(wr(n), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error allocating wr', 1)
      ALLOCATE(wi(n), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error allocating wi', 1)
      ALLOCATE(work(lwork), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error allocating work', 1)
      lwork = -1
      !
      wr(:)    = zero
      wi(:)    = zero
      vl(:, :) = zero
      vr(:, :) = zero
      work(:)  = zero
      !
      ! main solver
      CALL DGEEV(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, ierr)
      IF ( ierr /= 0) CALL errore('crit_temp_solver', 'Error eigenvalue solver failed!', 1)
      lwork = MIN(4 * n, INT(work(1)))
      DEALLOCATE(work, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error deallocating work', 1)
      ALLOCATE(work(lwork), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error allocating work', 1)
      CALL DGEEV(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, ierr)
      IF ( ierr /= 0) CALL errore('crit_temp_solver', 'Error eigenvalue solver failed!', 1)
      !
      !! Find the largest eigenvalue
      eigen = wr(1)
      DO ic = 2, n
        IF (eigen < wr(ic))   eigen = wr(ic)
      ENDDO
      niter = 1
      !
      DEALLOCATE(vl, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error deallocating vl', 1)
      DEALLOCATE(vr, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error deallocating vr', 1)
      DEALLOCATE(wr, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error deallocating wr', 1)
      DEALLOCATE(wi, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error deallocating wi', 1)
      DEALLOCATE(work, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error deallocating work',1)
    ENDIF
    !
    IF (tc_linear_solver == 'power') THEN
      !
      ALLOCATE(x(adim), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error allocating x', 1)
      ALLOCATE(y(adim), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error allocating y', 1)
      !
      n    = adim
      norm = 3.d0 * conv_thr_iaxis
      y(:) = 1.d0
      iter = 1
      ! main solver
      conv = .FALSE.
      DO WHILE (.NOT. conv .AND. iter < nsiter)
        norm = SQRT(SUM(y(1:n) ** two))
        x    = y / norm
        y(:) = zero
        DO ix = 1, n
          DO iy = 1, n
            y(ix) = y(ix) + a(ix, iy) * x(iy)
          ENDDO
        ENDDO
        alpha = zero
        DO ix = 1, n
          alpha = alpha + x(ix) * y(ix)
        ENDDO
        x = y - alpha * x
        norm  = SQRT(SUM(x(1:n) ** two))
        iter  = iter + 1
        IF (norm < conv_thr_iaxis) conv = .TRUE.
      ENDDO
      eigen = alpha
      niter = iter
      !
      DEALLOCATE(x, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error deallocating x', 1)
      DEALLOCATE(y, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_solver', 'Error deallocating y', 1)
      !
    ENDIF
    !
    IF (.NOT. conv .AND. iter == nsiter) THEN
      WRITE(stdout, '(/5x, a, i6)') 'Convergence (tc_linear) was not reached in nsiter = ', iter
      WRITE(stdout, '(5x, a)') 'Increase nsiter or reduce conv_thr_iaxis'
      CALL errore('crit_temp_solver', 'Convergence (tc_linear) was not reached', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE crit_temp_solver
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_matsubara_freq(itemp)
    !-----------------------------------------------------------------------
    !!
    !! SH: If iverbosity = 2; this routine writes the matsubara indices to
    !!       matsu-freq*.out file for each temperature (Nov 2021).
    !!
    USE kinds,         ONLY : DP
    USE supercond_common,     ONLY : nsiw, wsn
    USE global_var,    ONLY : gtemp
    USE ep_constants,  ONLY : kelvin2eV, pi
    USE io_var,        ONLY : iufilmat
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256) :: filename
    !! Output file name
    INTEGER           :: itemp
    !! Counter for temperature
    INTEGER           :: iw
    !! Loop variable for frequencies
    INTEGER           :: ierr
    !! Error status
    REAL(KIND = DP)   :: temp
    !! Temperature
    !
    ! Initiate output file names
    temp = gtemp(itemp) / kelvin2eV
    IF (temp < 10.d0) THEN
      WRITE(filename, 104) 'matsu-freq','_00', temp,'.out'
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
      WRITE(filename, 105) 'matsu-freq','_0' , temp,'.out'
    ELSEIF (temp >= 100.d0) THEN
      WRITE(filename, 106) 'matsu-freq','_'  , temp,'.out'
    ENDIF
    ! Write mats-freq*.out file
    OPEN(UNIT = iufilmat, FILE = filename, IOSTAT = ierr)
    IF (ierr /= 0) CALL errore('write_matsubara_freq', 'Error creating mats-freq.out', 1)
    DO iw =1, nsiw(itemp)
      WRITE(iufilmat,*) wsn(iw)
    ENDDO
    CLOSE(iufilmat)
    !
    104 FORMAT(a10, a3, f4.2, a4)
    105 FORMAT(a10, a2, f5.2, a4)
    106 FORMAT(a10, a1, f6.2, a4)
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE write_matsubara_freq
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_elphon()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by eliashberg_init,
    !!  read_frequencies, read_eigenvalues, read_kqmap, read_ephmat,
    !!  and evaluate_a2f_lambda
    !!
    USE input,         ONLY : liso, laniso
    USE global_var,    ONLY : wf, wqf, xqf, gtemp, bztoibz
    USE supercond_common,     ONLY : ekfs, xkfs, wkfs, g2, w0g, &
                              ixkff, ixkqf, ixqfs, nqfs, memlt_pool, &
                              ibnd_kfs_all_to_kfs, ibnd_kfs_to_kfs_all, &
                              ekfs_all, xkfs_all, wkfs_all
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    ! eliashberg_init
    IF (.NOT. liso .AND. .NOT. laniso) THEN
      DEALLOCATE(gtemp, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_eliashberg_aniso', 'Error deallocating gtemp', 1)
    ENDIF
    ! read_frequencies
    DEALLOCATE(wf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating wf', 1)
    DEALLOCATE(wqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating wqf', 1)
    DEALLOCATE(xqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating xqf', 1)
    ! read_eigenvalues
    DEALLOCATE(ibnd_kfs_all_to_kfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ibnd_kfs_all_to_kfs', 1)
    DEALLOCATE(ibnd_kfs_to_kfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ibnd_kfs_to_kfs_all', 1)
    DEALLOCATE(ekfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ekfs', 1)
    DEALLOCATE(xkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating xkfs', 1)
    DEALLOCATE(wkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating wkfs', 1)
    DEALLOCATE(ekfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ekfs_all', 1)
    DEALLOCATE(xkfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating xkfs_all', 1)
    DEALLOCATE(wkfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating wkfs_all', 1)
    DEALLOCATE(w0g, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating w0g', 1)
    ! read_kqmap
    DEALLOCATE(ixkff, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ixkff', 1)
    DEALLOCATE(bztoibz, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating bztoibz', 1)
    DEALLOCATE(ixkqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ixkqf', 1)
    DEALLOCATE(ixqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating ixqfs', 1)
    DEALLOCATE(nqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating nqfs', 1)
    DEALLOCATE(memlt_pool, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating memlt_pool', 1)
    ! read_ephmat
    DEALLOCATE(g2, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_eliashberg_elphon', 'Error deallocating g2', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_eliashberg_elphon
    !-----------------------------------------------------------------------
    !
  !----------------------------------------------------------------------
  END MODULE supercond
  !----------------------------------------------------------------------
