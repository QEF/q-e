  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------!
  MODULE supercond_vertex
  !----------------------------------------------------------------------
  !!
  !! Origionally by S. Mishra -----------------------------------------------------
  !! Modules that contains all the subroutines related to imaga parallized isotropic 
  !! adiabatic Eliashberg theory
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-------------------------------------------------------------------------------!
    SUBROUTINE prepare_a2f(nrr_k, nrr_q, nrr_g, irvec_r, irvec_q, irvec_g, ndegen_k, &
                        ndegen_q, ndegen_g, nrws, rws, dims, dims2)
    !-------------------------------------------------------------------------------!
    !! S. Mishra: This rotuine checks the restart files for a2f calculated using 
    !! image parallelization (April 2025)
    !
    USE kinds,            ONLY : DP
    USE global_var,       ONLY : totq
    USE io_global,        ONLY : stdout
    USE io_files,         ONLY : prefix, tmp_dir, check_tempdir
    USE io_var,           ONLY : iunrestart
    USE input,            ONLY : fila2f
    USE ions_base,        ONLY : nat
    USE mp,               ONLY : mp_bcast, mp_sum
    USE mp_global,        ONLY : npool, my_pool_id, inter_image_comm, &
                                 nimage, my_image_id, inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: dims2
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, INTENT(in) :: nrr_k
    !! number of electronic WS points
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: nrr_g
    !! number of el-ph WS points
    INTEGER, INTENT(in) :: irvec_q(3, nrr_k)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Coordinates of real space vector for electron-phonon
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(in) :: ndegen_q(nrr_q, dims2, dims2)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, nat)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    INTEGER, INTENT(in) :: nrws
    !! Number of Wigner-Size real space vectors when lifc == .TRUE.
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Coordinates of real space vector for electrons
    REAL(KIND = DP), INTENT(in) :: rws(0:3, nrws)
    !! Real space Wigner-Seitz vector when lifc == .TRUE.
    !
    !! Local variables
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !! Convert integer IDs to character strings
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save restart and a2f files
    CHARACTER(LEN = 256) :: fila2f_tmp
    !! Name temprary a2f file
    CHARACTER(LEN = 256) :: fileiq_tmp
    !! Name tempporary files
    LOGICAL :: exst0, exst, exst2, pfs
    !! If the file exist
    INTEGER :: iq_restart
    !! Counter on coarse q-point grid
    INTEGER :: startq, lastq
    !! q-points sinto images
    INTEGER :: startq_old, lastq_old
    !! For restart options startq and lastq
    INTEGER :: totproc, nimage_old, npool_old
    !! Image-paralleilization check
    INTEGER :: ios
    !! INTEGER variable for I/O control
    !
    !
    IF (fila2f == ' ') WRITE(fila2f, '(a, a4)') TRIM(prefix), '.a2f'
    INQUIRE(FILE = fila2f, EXIST = exst0)
    !
    IF (exst0) THEN
      WRITE(stdout,'(a)') "                      "
      WRITE(stdout,'(a)') "     ------------------------------------------------------------------------ "
      WRITE(stdout,'(a)') "     ------------- WARNING -------------    WARNING  --------  WARNING ------ "
      WRITE(stdout,'(a)') "     a2f file is found. No need to recalculate a2f again.                     "
      WRITE(stdout,'(a)') "     If you want to recalculate a2f using different set of parameters.        "
      WRITE(stdout,'(a)') "     first remove .a2f file and then resubmit calculation.                    "
      WRITE(stdout,'(a)') "     ------------------------------------------------------------------------ "
      WRITE(stdout,'(a)') "  "
    ELSE 
      CALL divide(inter_image_comm, totq, startq, lastq)
      iq_restart = startq
      dirname = TRIM(tmp_dir) // 'tmp_restart/' 
      CALL check_tempdir(dirname, exst, pfs)
      fileiq_tmp = TRIM(dirname) // 'restartq_' // &
           & TRIM(int_to_char(my_image_id + 1)) // '_' // &
           & TRIM(int_to_char(my_pool_id + 1)) // '.fmt'
      fila2f_tmp = TRIM(dirname) // 'a2f_tmp_' // &
           & TRIM(int_to_char(my_image_id + 1)) // '_' // &
           & TRIM(int_to_char(my_pool_id + 1))
      totproc = 0
      INQUIRE(FILE = fileiq_tmp, EXIST = exst)
      INQUIRE(FILE = fila2f_tmp, EXIST = exst2)
      IF (exst .AND. exst2) THEN
        exst = .TRUE.
      ELSE
        exst = .FALSE.
      ENDIF
      IF (exst) totproc = 1
      CALL mp_sum(totproc, inter_pool_comm)
      CALL mp_sum(totproc, inter_image_comm)
      IF (totproc == nimage * npool) THEN
        exst = .TRUE.
      ELSE
        exst = .FALSE.
      ENDIF
      !
      IF (exst) THEN
        !OPEN(UNIT = iunrestart, FILE = fileiq_tmp, STATUS = 'old', IOSTAT = ios)
        OPEN(UNIT = iunrestart, FILE = fileiq_tmp, STATUS = 'old', FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
        READ(iunrestart) startq_old
        READ(iunrestart) lastq_old
        READ(iunrestart) iq_restart
        READ(iunrestart) nimage_old
        READ(iunrestart) npool_old
        CLOSE(iunrestart)
        ! DEBUG
        WRITE(stdout, '(5x,a,i8,a, i8)') 'nimage_old = ', nimage_old, ' nimage = ', nimage
        WRITE(stdout, '(5x,a,i8,a, i8)') 'npool_old = ', npool_old, ' npool = ', npool
        !
        IF (nimage .NE. nimage_old) CALL errore('prepare_a2f', 'nimage is different from the previous one. ', 1)
        IF (npool  .NE. npool_old ) CALL errore('prepare_a2f', 'npool is different from the previous one. ', 1)
        IF (startq .NE. startq_old) CALL errore('prepare_a2f', 'startq is wrong.', 1)
        IF (lastq  .NE. lastq_old ) CALL errore('prepare_a2f', 'lastq is wrong. ', 1)
        !
        IF (iq_restart < lastq_old) THEN
          WRITE(stdout, '(5x,a,i8,a)') 'We restart from ', iq_restart, ' q-points'
          iq_restart = iq_restart + 1
          CALL calculate_a2f(iq_restart, nrr_k, nrr_q, nrr_g, irvec_r, irvec_q, irvec_g, &
                ndegen_k, ndegen_q, ndegen_g, nrws, rws, dims, dims2)
          !
        ELSE
          WRITE(stdout, '(5x,a)') 'All q-points are done, no need to restart !!'
        ENDIF
      ELSE  ! no restartq.fmt file present
        CALL calculate_a2f(iq_restart, nrr_k, nrr_q, nrr_g, irvec_r, irvec_q, irvec_g, &
              ndegen_k, ndegen_q, ndegen_g, nrws, rws, dims, dims2)
      ENDIF ! restart
    ENDIF ! a2f present or not
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE prepare_a2f
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE calculate_a2f(iq_restart, nrr_k, nrr_q, nrr_g, irvec_r, irvec_q, &
                  irvec_g, ndegen_k, ndegen_q, ndegen_g, nrws, rws, dims, dims2)
    !-----------------------------------------------------------------------------
    !
    !! Here we calculate a2f using isotropic Eliashberg approximation
    !!
    !!
    USE kinds,            ONLY : DP
    USE pwcom,            ONLY : ef
    USE input,            ONLY : eps_acoustic, nqstep, degaussq, nbndsub, fsthick,   &
                                 longrange_only, ngaussw, degaussw, lifc, mp_mesh_k, &
                                 restart_step, restart, lfast_kmesh
    USE global_var,       ONLY : wqf, wf, xqf, wkf, nkf, nkqf, chw, nktotf,  nqtotf, xkf, &
                                 ibndmin, ibndmax, bztoibz, totq, selecq
    USE bzgrid,           ONLY : kpmq_map, xqf_otf
    USE wannier2bloch,    ONLY : hamwan2bloch, dynwan2bloch, ephwan2blochp, ephwan2bloch, &
                                 dynifc2blochf
    USE wigner,           ONLY : wigner_divide_ndegen
    USE modes,            ONLY : nmodes
    USE cell_base,        ONLY : at, bg
    USE supercond_common, ONLY : wsph, dwsph, wsphmax
    USE io_global,        ONLY : stdout
    USE mp,               ONLY : mp_barrier, mp_sum
    USE mp_global,        ONLY : inter_pool_comm, inter_image_comm, my_pool_id, my_image_id
    USE ions_base,        ONLY : nat, amass, ityp
    USE ep_constants,     ONLY : ryd2ev, ryd2mev, one, two, zero, czero, eps40, ci, twopi, &
                                 cone
    USE parallelism,      ONLY : poolgather
    USE check_stop,       ONLY : check_stop_now
    USE stop,             ONLY : stop_epw
    ! --------------------------------------------------------------------------------------!
    !
    INTEGER, INTENT(IN) :: iq_restart
    !! restart from iq_restart
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: dims2
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, INTENT(in) :: nrr_k
    !! number of electronic WS points
    INTEGER, INTENT(in) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(in) :: nrr_g
    !! number of el-ph WS points
    INTEGER, INTENT(in) :: irvec_q(3, nrr_k)
    !! Coordinates of real space vector for phonons
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Coordinates of real space vector for electron-phonon
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(in) :: ndegen_q(nrr_q, dims2, dims2)
    !! degeneracy of WS points for phonon
    INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, nat)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    INTEGER, INTENT(in) :: nrws
    !! Number of Wigner-Size real space vectors when lifc == .TRUE.
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !! Coordinates of real space vector for electrons
    REAL(KIND = DP), INTENT(in) :: rws(0:3, nrws)
    !! Real space Wigner-Seitz vector when lifc == .TRUE.
    !
    ! Local variable
    INTEGER :: nk, nkq
    !! Index of the k point from the full bzgrids.
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! Counter on coarse k-point bzgrids
    INTEGER :: iq
    !! Counter on coarse q-point bzgrids
    INTEGER :: ikk
    !! Counter on k-point when you have paired k and q
    INTEGER :: ikq
    !! Paired counter so that q is adjacent to its k
    INTEGER :: na
    !! Counter on atom
    INTEGER :: mu, nu, imode
    !! counter on mode
    INTEGER :: ibnd, jbnd
    !! Counter on band indices
    INTEGER :: iqq
    !! current iq-point
    INTEGER :: iwph
    !! Counter over frequencies omega
    INTEGER :: startq, lastq
    !! Lower/upper bound index for q-points
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    REAL(KIND = DP) :: xkk(3)
    !! Current k-point on the fine bzgrids
    REAL(KIND = DP) :: xkq(3)
    !! current k+q point on the fine bzgrids
    REAL(KIND = DP) :: weight
    !! Product of delta functions in energy and dos
    REAL(KIND = DP) :: weightq
    !! Delta function in phonon mode-q
    REAL(KIND = DP) :: ekk, ekq
    !! Eigen energy at k and k+q on the fine bzgrids relative to the Fermi level
    REAL(KIND = DP) :: dosef
    !! Electronic density of states
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: inv_degaussw
    !! Invese dos inv_degaussw = 1/degaussw. Defined for efficiency reason
    REAL(KIND = DP) :: inv_degaussq
    !! Invese dos inv_degaussq = 1/degaussq. Defined for efficiency reason
    REAL(KIND = DP):: g2                  
    !! Product of two g matrices 
    REAL(KIND = DP) :: wq, inv_wq
    !! phonon freq and its inverse on the fine bzgrids
    REAL(KIND = DP):: lambda_eph
    !! total e-ph coupling strength (a2f integration)
    REAL(KIND = DP):: l_sum
    !! total e-ph coupling strength
    REAL(KIND = DP), ALLOCATABLE :: rdotk(:)
    !! $r\cdot k$
    REAL(KIND = DP), ALLOCATABLE :: rdotk2(:)
    !! $r\cdot k+q$
    REAL(KIND = DP), ALLOCATABLE :: w2(:)
    !! Interpolated phonon frequency
    REAL(KIND = DP), ALLOCATABLE :: w0g1(:, :)
    !! Delta function in energy at k and k+q
    REAL(KIND = DP), ALLOCATABLE :: a2f(:) 
    !!for a2f_iso testing
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Eigen-energies of full k-bzgrids from all pools.
    REAL(KIND = DP), ALLOCATABLE :: etf_loc(:, :)
    !! Eigen-energies local full k-bzgrids.
    REAL(KIND = DP), ALLOCATABLE :: etf_tmp(:)
    !! Temporary Eigen-energies at k-bzgrids.
    REAL(KIND = DP), ALLOCATABLE :: wkf_loc(:)
    !! Local weight of k-points in parallel case
    REAL(KIND = DP), ALLOCATABLE :: wkf_all(:)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: a2f_modeproj(:, :)
    !!for a2f_iso mode-projected
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the density of states at the Fermi level
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkk(:, :)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkq(:, :)
    !! the same, for points k+q
    COMPLEX(KIND = DP), ALLOCATABLE :: uf(:, :)
    !! Rotation matrix for phonons
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE :: cfacq(:)
    !! Used to store $e^{2\pi r \cdot k+q}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwef(:, :, :, :)
    !! e-p matrix  in el wannier - fine Bloch phonon grid
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatf(:, :, :)
    !! e-p matrix  in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: rrf_tmp(:, :, :)
    !! Interpolation position matrix elements on the fine mesh (ipol, nbnd, nbnd)
    COMPLEX(KIND = DP), ALLOCATABLE :: tmp(:, :, :)
    !! Temporary overlap at k and k+q
    !!
    !! Allocation
    ALLOCATE(cfac(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating cfac', 1)
    ALLOCATE(cfacq(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating cfacq', 1)
    ALLOCATE(rdotk(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating rdotk', 1)
    ALLOCATE(rdotk2(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating rdotk2', 1)
    ALLOCATE(w2(3 * nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating w2', 1)
    ALLOCATE(w0g1(nbndsub, nkqf), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating w0g1', 1)
    ALLOCATE(etf_loc(nbndsub, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating etf_loc', 1)
    ALLOCATE(etf_tmp(nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating etf_tmp', 1)
    ALLOCATE(etf_all(nbndsub, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating etf_all', 1)
    ALLOCATE(wkf_loc(nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating wkf_loc', 1)
    ALLOCATE(wkf_all(nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating wkf_all', 1)
    ALLOCATE(cufkk(nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating cufkk', 1)
    ALLOCATE(cufkq(nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating cufkq', 1)
    ALLOCATE(uf(nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating uf', 1)
    ALLOCATE(epmatwef(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating epmatwef', 1)
    ALLOCATE(epmatf(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating epmatf', 1)
    ALLOCATE(rrf_tmp(3, nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating rrf', 1)
    ALLOCATE(tmp(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating rrf', 1)
    ALLOCATE(wsph(nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('supercond_vertex', 'Error allocating wsph', 1)
    ALLOCATE(a2f(nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating a2f', 1)
    ALLOCATE(a2f_modeproj(nmodes, nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error allocating a2f_modeproj', 1)
    !
    !------------------------------------------------------------!
    !
    !start the clock
    CALL start_clock('a2f_iso')
    !
    WRITE(stdout,'(/5x,a)') 'Entering adiabatic isotropic Eliashberg'
    !
    !Initialization
    !
    cfac(:)          = czero
    cfacq(:)         = czero
    rdotk(:)         = zero
    rdotk2(:)        = zero
    w2(:)            = zero
    w0g1(:, :)       = zero
    etf_loc(:, :)    = zero
    etf_all(:, :)    = zero
    etf_tmp(:)       = zero
    wkf_all(:)       = zero
    wkf_loc(:)       = zero
    dosef            = zero
    wsph(:)          = zero
    a2f(:)           = zero
    a2f_modeproj(:, :) = zero
    g2               = zero
    l_sum            = zero
    rrf_tmp(:, :, :) = czero
    epmatwef(:, :, :, :)  = czero
    !
    ! Multiplication is faster than division ==> Important if called a lot
    ! in inner loops
    inv_degaussw = one / degaussw
    inv_degaussq = one / degaussq 
    !
    ! After load_rebal we reordered wkf and xkf. So we need to recalculate 
    ! etf again to correctly calculate dosef
    xxq = 0.d0
    !
    ! nkf is the number of kpoints in the pool
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      xkk = xkf(:, ikk)
      wkf_loc(ik) = wkf(ikk)
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
      cfac(:) = EXP(ci * rdotk(:))
      CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf_loc(:, ik), chw, cfac, .FALSE.)
    ENDDO
    CALL poolgather(nbndsub, nktotf, nkf, etf_loc, etf_all)
    CALL poolgather(1, nktotf, nkf, wkf_loc, wkf_all)
    !
    ! Fermi level and corresponding DOS
    !dosef = dos_ef(0, degaussw, ef, etf, wkf, nkqf, nbndsub)
    dosef = dos_ef(ngaussw, degaussw, ef, etf_loc, wkf_loc, nkf, nbndsub)
    ! N(Ef) is the DOS per spin
    dosef = dosef / two
    inv_dos = one / dosef
    !
    WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
    WRITE (stdout, 101) dosef / ryd2ev, ef * ryd2ev
    !
    ! To get the maximum of phonon frequencies, we have to call dynwan2bloch.
    DO iq = 1, nqtotf
      ! xqf has to be in crystal coordinate
      IF (lfast_kmesh) THEN
        ! The q-point coordinate is generate on the fly for each q-point
        CALL xqf_otf(iq, xxq)
      ELSE
        xxq = xqf(:, iq)
      ENDIF
      ! ------------------------------------------------------
      ! dynamical matrix : Wannier -> Bloch
      ! ------------------------------------------------------
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf, w2, .false.)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq, uf, w2, .false.)
      ENDIF
      ! ...then take into account the mass factors and square-root the frequencies...
      DO nu = 1, nmodes
        ! wf are the interpolated eigenfrequencies (omega on fine bzgrids)
        IF (w2(nu) > 0.0d0) THEN
          wf(nu, iq) =  DSQRT(ABS(w2(nu)))         
        ELSE
          wf(nu, iq) = -DSQRT(ABS(w2(nu)))
        ENDIF
      ENDDO !nu
    ENDDO ! iq
    !
    !Maximum phonon freqeuncy in Ry
    wsphmax = 1.1d0 * MAXVAL(wf(:, :)) ! increase by 10%
    WRITE(stdout,'(5x, "nmodes = ", i9)') nmodes
    WRITE(stdout,'(5x, "w_max (meV) = ", f12.7)') wsphmax * ryd2mev
    ! create phonon bzgrids in Ry 
    dwsph = wsphmax / DBLE(nqstep)
    !
    DO iwph = 1, nqstep
      wsph(iwph) = DBLE(iwph) * dwsph
    ENDDO
    !
    ! Distribute the q-points into images
    CALL divide(inter_image_comm, totq, startq, lastq)
    !
    ! Read from the file if iq_restart > 1
    IF (iq_restart .NE. startq) THEN
      CALL read_restart_a2f(a2f, l_sum)
      WRITE(stdout, '(5x, a, i9, a, i9)') 'We read from iq_restart / lastq = ', iq_restart, ' / ', lastq
    ENDIF
    !
    DO iqq = iq_restart, lastq
      WRITE(stdout, '(5x, a, i10, a, i10)' ) 'Progression iq (fine) = ', iqq - startq + 1, ' / ', lastq - startq + 1
      !
      xxq = 0
      iq  = selecq(iqq)
      xxq = xqf(:, iq)
      ! ------------------------------------------------------
      ! dynamical matrix : Wannier -> Bloch
      ! ------------------------------------------------------
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf, w2, .false.)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq, uf, w2, .false.)
      ENDIF
      !
      DO nu = 1, nmodes
        ! wf are the interpolated eigenfrequencies
        IF (w2(nu) > 0.0d0) THEN
          wf(nu, iq) =  DSQRT(ABS(w2(nu)))         
        ELSE
          wf(nu, iq) = -DSQRT(ABS(w2(nu)))
        ENDIF
        DO mu = 1, nmodes
          na = (mu - 1) / 3 + 1
          uf(mu, nu) = uf(mu, nu) / DSQRT(amass(ityp(na)))
        ENDDO !mu
      ENDDO !nu
      !
      ! --------------------------------------------------------------
      ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
      ! --------------------------------------------------------------
      IF (.NOT. longrange_only) THEN
        CALL ephwan2blochp(nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatwef, nbndsub, nrr_k, dims, dims2)
      ENDIF
      CALL wigner_divide_ndegen(epmatwef, 1, nbndsub, nrr_k, nmodes, ndegen_k, dims)
      !
      ! This is a loop over k blocks in the pool (size of the local k-set)
      DO ik = 1, nkf
        ! xkf is assumed to be in crys coord
        ikk = 2 * ik - 1
        ikq = ikk + 1
        ! Compute coordinates at k and k + q
        xkk = xkf(:, ikk)
        xkq = xkk + xxq
        !
        CALL kpmq_map(xkk, (/0d0, 0d0, 0d0/), 1, nk)
        CALL kpmq_map(xkk, xxq, 1, nkq)
        IF (nk == 0 .OR. nkq == 0) CALL errore ('calculate_a2f', 'nk or nkq cannot be 0', 1)
        ! 
        IF (mp_mesh_k) THEN
          nk  = bztoibz(nk)
          nkq = bztoibz(nkq)
        ENDIF
        !
        CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
        cfac(:) = EXP(ci * rdotk(:))
        CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf_tmp(:), chw, cfac, .FALSE.)
        !
        CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk2, 1)
        cfacq(:) = EXP(ci * rdotk2(:))
        CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf_tmp(:), chw, cfacq, .FALSE.)
        !
        ! --------------------------------------------------------------
        ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch p
        ! --------------------------------------------------------------
        !
        CALL cryst_to_cart(1, xxq, bg, 1)
        tmp(:, :, :)    = czero
        epmatf(:, :, :) = czero
        CALL ephwan2bloch(nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, nmodes, cfac, xxq, rrf_tmp)
        CALL cryst_to_cart(1, xxq, at, -1)
        !
        DO ibnd = ibndmin, ibndmax 
          ekk = etf_all(ibnd, nk) - ef
          DO jbnd = ibndmin, ibndmax
            ekq = etf_all(jbnd, nkq) - ef
            IF ((ABS(ekk) < fsthick) .AND. (ABS(ekq) < fsthick)) THEN  
              !delta[E_k - Ef]
              w0g1(ibnd, ikk) = w0gauss(ekk * inv_degaussw, 0) * inv_degaussw
              w0g1(jbnd, ikq) = w0gauss(ekq * inv_degaussw, 0) * inv_degaussw
              weight          = w0g1(ibnd, ikk) * w0g1(jbnd, ikq) * wkf(ikk) * wqf(iq)
              lambda_eph      = zero
              ! Now do the eigenvector rotation: epmatf(j) = sum_i eptmp(i) * uf(i,j)
              tmp(:, :, :) =  epmatf(:, :, :)
              CALL ZGEMM('n', 'n', nbndsub * nbndsub, nmodes, nmodes, cone, tmp, &
                          nbndsub * nbndsub, uf, nmodes, czero, epmatf, nbndsub * nbndsub)
              !
              DO imode = 1, nmodes
                wq = wf(imode, iq)
                inv_wq = one / (two * wq)
                IF (wq > eps_acoustic) THEN 
                  g2 = ABS(epmatf(jbnd, ibnd, imode))**two * inv_wq
                  ! WRITE(my_pool_id + 9001, '(3G12.4, 2f12.4, E12.5)') iq, ik, imode, ryd2ev * ekk, &
                  !                     ryd2ev * ekq, ryd2ev * ryd2ev * g2
                  lambda_eph = lambda_eph + g2 / wq
                  !
                  DO iwph = 1, nqstep 
                    weightq   = w0gauss((wsph(iwph) - wq) * inv_degaussq, 0) * inv_degaussq
                    a2f(iwph) = a2f(iwph) + weight * g2 * weightq
                    a2f_modeproj(imode, iwph) = a2f_modeproj(imode, iwph) + weight * g2 * weightq
                  ENDDO ! iwph
                ENDIF ! eps-acous
              ENDDO !imode
              IF (lambda_eph > 0.d0) l_sum = l_sum + weight * lambda_eph
            ENDIF !etf-i,j
          ENDDO !jbnd
        ENDDO !ibnd
      ENDDO ! ik
      !! S. Mishra: Check if it's time to write restart files 
      IF (restart .AND. (MOD((iqq - startq + 1), restart_step) == 0 .OR. iqq == lastq)) THEN
        WRITE(stdout, '(a)') ' '
        WRITE(stdout, '(5x, a, i10, a, i10)' ) 'Writing a2f_tmp to file at iqq (fine) = ', iqq, '/', lastq
        WRITE(stdout, '(a)') ' '
        ! FIX ME: create restart for mode-projected a2f
        !CALL write_restart_a2f(startq, lastq, iqq, l_sum, a2f, a2f_modeproj)
        CALL write_restart_a2f(startq, lastq, iqq, a2f, l_sum)
      ENDIF ! restart
      IF ( check_stop_now() ) THEN
       !  In this case the code stops inside after writing the files
        CALL stop_epw()
      ENDIF
      !
    ENDDO  ! iqq
    ! collect contributions from all images (sum over q-points and k-points)
    ! We donot need barrier, mp_sum does that by default
    CALL mp_sum(l_sum, inter_pool_comm)
    CALL mp_sum(a2f, inter_pool_comm)
    CALL mp_sum(a2f_modeproj, inter_pool_comm)
    !CALL mp_barrier(inter_pool_comm)
    CALL mp_sum(l_sum, inter_image_comm)
    CALL mp_sum(a2f, inter_image_comm)
    CALL mp_sum(a2f_modeproj, inter_image_comm)
    !CALL mp_barrier(inter_image_comm)
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'DOS =', f10.6, ' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
    !
    a2f(:) = 0.5d0 * a2f(:) * inv_dos
    l_sum  = l_sum * inv_dos
    WRITE(stdout, '(/,5x,a, f10.6)' ) 'eps_acoustic (mev) = ', eps_acoustic * ryd2mev
    WRITE(stdout, '(/,5x,a, f10.6)' ) 'dosef (1/ev) = ', dosef / ryd2ev
    WRITE(stdout, '(/,5x,a, f10.6)' ) 'l_sum = ', l_sum
    WRITE(stdout, '(/,5x,a, f10.6)' ) 'fsthick = ', fsthick * ryd2ev
    !
    !! Now, write a2f to file named prefix.a2f
    CALL write_a2f(a2f, a2f_modeproj, l_sum, dosef)
    !! Remove the temporary restart files to clean the disk space
    IF (restart) CALL delete_restart_file()
    !
    !! Deallocate
    CALL deallocate_adia_a2f(rdotk, rdotk2, w2, w0g1, etf_all, etf_loc, etf_tmp,  wkf_loc, wkf_all, &
                                  wsph, epmatwef, epmatf, rrf_tmp, tmp, cufkk, cufkq, uf, cfac, cfacq)
    
    DEALLOCATE(a2f, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating a2f', 1) 
    DEALLOCATE(a2f_modeproj, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating a2f_modeproj', 1) 
    !
    CALL stop_clock('a2f_iso')
    CALL print_clock('a2f_iso')
    WRITE(stdout, '(a)') ' '
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE calculate_a2f
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_a2f(a2f, a2f_modeproj, l_sum, dosef)
    !-----------------------------------------------------------------------
    !! S. Mishra: Writing a2f to prefix.a2f file 
    !
    USE kinds,            ONLY : DP
    USE input,            ONLY : nqstep, degaussq, degaussw, fsthick, nqsmear
    USE modes,            ONLY : nmodes
    USE io_files,         ONLY : prefix
    USE io_var,           ONLY : iua2ffil
    USE io_global,        ONLY : stdout, ionode
    USE ep_constants,     ONLY : ryd2eV, ryd2mev, zero
    USE io_files,         ONLY : prefix
    USE supercond_common, ONLY : wsph, dwsph
    !
    REAL(DP), INTENT(in) :: a2f(nqstep)
    !! isotropic a2f
    REAL(DP), INTENT(in) :: a2f_modeproj(nmodes, nqstep)
    !! isotropic a2f
    REAL(DP), INTENT(in) :: l_sum
    !! total e-ph coupling strength
    REAL(DP), INTENT(in) :: dosef
    !! DOS
    !
    !! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name temporary a2f_vertex file
    INTEGER :: iwph
    !! Counter over frequencies omega
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: l_a2f_tmp
    !! Temporary lambda
    REAL(KIND = DP), ALLOCATABLE :: l_a2f(:)
    !! total e-ph coupling strength (a2f integration)
    !
    CALL start_clock('write_a2f')
    ! SM: Only ionode writes to file
    IF (ionode) THEN 
      !
      name1 = TRIM(prefix) // '.a2f'
      OPEN(UNIT = iua2ffil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('write_a2f', 'error opening file ' // name1, iua2ffil)
      !
      ALLOCATE(l_a2f(nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_a2f', 'Error allocating l_a2f', 1)
      l_a2f(:) = zero
      !
      WRITE(iua2ffil, '(" w[meV] a2f and integrated 2*a2f/w for ", i4, " smearing values")') nqsmear
      !
      l_a2f_tmp = zero
      !
      DO iwph = 1, nqstep
        !l_a2f_tmp = l_a2f_tmp + 2.0d0 * (a2f(iwph) / wsph(iwph)) * dwsph
        l_a2f_tmp = l_a2f_tmp + 2.0d0 * (a2f(iwph) / DBLE(iwph))
        l_a2f(iwph) = l_a2f_tmp
        WRITE(iua2ffil, '(f12.7, 20f12.7)') wsph(iwph) * ryd2meV, a2f(iwph), l_a2f(iwph)
      ENDDO !iwph
      !
      WRITE(iua2ffil, *) "Integrated el-ph coupling"
      WRITE(iua2ffil, '("  #         ", f12.7)') l_a2f(nqstep)
      WRITE(iua2ffil, *) "Phonon smearing (meV)"
      WRITE(iua2ffil, '("  #         ", f12.7)') (degaussq * ryd2meV)
      WRITE(iua2ffil, '("Electron smearing (eV)", f12.7)') degaussw * ryd2eV
      WRITE(iua2ffil, '("Fermi window (eV)", f12.7)') fsthick * ryd2eV
      WRITE(iua2ffil, '("DOS (eV)", f12.7)') dosef / ryd2eV
      WRITE(iua2ffil, '("Summed el-ph coupling ", f12.7)') l_sum 
      CLOSE(iua2ffil)
      !
      name1 = TRIM(prefix) // '.a2f_proj'
      OPEN(UNIT = iua2ffil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('write_a2f', 'error opening file ' // name1, iua2ffil)
      !
      WRITE(iua2ffil, '("w[meV] a2f a2f_modeproj")')
      DO iwph = 1, nqstep
        ! wsph in meV (from eV)
        WRITE(iua2ffil, '(f12.7, 100f12.7)') wsph(iwph) * ryd2meV, a2f(iwph), a2f_modeproj(:, iwph)
      ENDDO
      WRITE(iua2ffil, '(a, f18.7, a, f18.7)') 'lambda_int = ', l_a2f(nqstep), '   lambda_sum = ', l_sum
      CLOSE(iua2ffil)
      !
      DEALLOCATE(l_a2f, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_a2f', 'Error deallocating l_a2f', 1)  
    ENDIF
    !
    CALL stop_clock('write_a2f')
    CALL print_clock('write_a2f')
    WRITE(stdout, '(a)') ' '
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE write_a2f
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_restart_a2f(startq, lastq, iqq, a2f, l_sum)
    !-----------------------------------------------------------------------
    !! Write restartq.fmt file and temporary restart files for 
    !! a2f after each restart_intervals
    !
    USE kinds,            ONLY : DP
    USE input,            ONLY : nqstep
    USE io_files,         ONLY : prefix, tmp_dir
    USE io_var,           ONLY : iua2ffil, iunrestart
    USE mp_global,        ONLY : my_pool_id, npool
    USE mp_images,        ONLY : nimage, my_image_id
    !
    INTEGER, INTENT(in) :: startq
    !! total q+p-points 
    INTEGER, INTENT(in) :: lastq
    !! total q+p-points 
    INTEGER, INTENT(in) :: iqq
    !! total q+p-points 
    REAL(DP), INTENT(inout) :: a2f(nqstep)
    !! a2f file
    !REAL(DP), INTENT(inout) :: a2f_modeproj(nmodes, nqstep)
    !! a2f mode-projected file
    REAL(DP), INTENT(inout) :: l_sum
    !! total e-ph coupling strength non-adiabatic
    !
    !Local variables
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !! Convert integer IDs to character strings
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save files
    CHARACTER(LEN = 256) :: fila2f_tmp
    !! Name temprary a2f_vertex file
    CHARACTER(LEN = 256) :: fileiq_tmp
    !! Name eigenvalue file
    INTEGER :: iwph
    !! Counter over frequencies omega
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ierr
    !! Error status
    !
    dirname = TRIM(tmp_dir) // 'tmp_restart/' 
    !CALL create_directory(dirname)
    !
    fileiq_tmp = TRIM(dirname) // 'restartq_' // &
           & TRIM(int_to_char(my_image_id + 1)) // '_' // &
           & TRIM(int_to_char(my_pool_id + 1)) // '.fmt'
    OPEN(UNIT = iunrestart, FILE = fileiq_tmp, STATUS = 'unknown', FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
    !OPEN(UNIT = iunrestart, FILE = fileiq_tmp, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('write_restart_files', 'Error opening file restartq.fmt', 1)
    WRITE(iunrestart) startq
    WRITE(iunrestart) lastq
    WRITE(iunrestart) iqq
    WRITE(iunrestart) nimage
    WRITE(iunrestart) npool
    CLOSE(iunrestart)
    !
    ! Writing a2f_tmp to file
    fila2f_tmp = TRIM(dirname) // 'a2f_tmp_' // &
           & TRIM(int_to_char(my_image_id + 1)) // '_' // &
           & TRIM(int_to_char(my_pool_id + 1))
    OPEN(UNIT = iua2ffil, FILE = fila2f_tmp, STATUS = 'unknown', FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
    IF (ios /= 0) CALL errore('write_restart_files', 'Error opening file prefix.a2f_tmp', 1) 
    !WRITE(iua2ffil, '(E15.6E3)') l_sum
    WRITE(iua2ffil) l_sum
    !! loop over w_ph
    DO iwph = 1, nqstep
      WRITE(iua2ffil) a2f(iwph)
      !WRITE(iua2ffil) a2f_modeproj(:, iwph)
    ENDDO
    !
    CLOSE(iua2ffil)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE write_restart_a2f
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_restart_a2f(a2f, l_sum)
    !-----------------------------------------------------------------------
    ! SM: This is for reading restartq.fmt file and a2f_tmp files for restart
    ! 
    USE kinds,            ONLY : DP
    USE input,            ONLY : nqstep
    USE io_files,         ONLY : prefix, tmp_dir
    USE io_var,           ONLY : iua2ffil
    USE mp_global,        ONLY : my_pool_id, my_image_id
    USE io_global,        ONLY : stdout, ionode_id
    USE ep_constants,     ONLY : zero
    !
    REAL(DP), INTENT(inout) :: a2f(nqstep)
    !! a2f file
    !REAL(DP), INTENT(inout) :: a2f_modeproj(nmodes, nqstep)
    !! a2f mode-projected file
    REAL(DP), INTENT(inout) :: l_sum
    !! total adiabatic e-ph coupling strength 
    !
    !! Local variables
    CHARACTER(LEN = 256) :: fila2f_tmp
    !! Name temprary a2f_vertex file
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !! Convert integer IDs to character strings
    INTEGER :: iwph
    !! Counter over frequencies omega
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ierr
    !! Error status
    !
    dirname = TRIM(tmp_dir) // 'tmp_restart/' 
    fila2f_tmp = TRIM(dirname) // 'a2f_tmp_' // &
           & TRIM(int_to_char(my_image_id + 1)) // '_' // &
           & TRIM(int_to_char(my_pool_id + 1))
    !OPEN(UNIT = iua2ffil, FILE = fila2f_tmp, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    OPEN(UNIT = iua2ffil, FILE = fila2f_tmp, STATUS = 'unknown', FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
    IF (ios /= 0) CALL errore('read_restart_a2f', 'Error opening file prefix.a2f_tmp', 1) 
    READ(iua2ffil) l_sum
    !! loop over w_ph
    DO iwph = 1, nqstep
      !READ(iua2ffil, '(E15.6E3)') a2f(iwph)
      READ(iua2ffil) a2f(iwph)
    ENDDO
    CLOSE(iua2ffil)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_restart_a2f
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE delete_restart_file()
    !-----------------------------------------------------------------------
    !
    !! SM: This is for removing restartq.fmt file and a2f_vertex_tmp files
    !
    USE io_files,         ONLY : prefix, tmp_dir, check_tempdir, delete_if_present
    USE mp_global,        ONLY : my_pool_id, my_image_id
    USE io_global,        ONLY : stdout
    USE input,            ONLY : a2f_iso
    USE io_var,           ONLY : iua2ffil, iunrestart
    !
    !Local variables
    CHARACTER(LEN = 256) :: fileiq_tmp
    !! Temporary restart_qp files
    CHARACTER(LEN = 256) :: fila2f_tmp
    !! Name temporary a2f_vertex files
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !! Convert integer IDs to character strings
    LOGICAL :: exst, exst2, pfs
    !! Does the file exist
    INTEGER :: ios
    !! INTEGER variable for I/O control
    !
    !! delete the temporary file if it exists.
    dirname = TRIM(tmp_dir) // 'tmp_restart/' 
    CALL check_tempdir(dirname, exst, pfs)
    !
    IF (a2f_iso) THEN 
      WRITE(stdout, '(/,5x, a)' ) 'Removing temporary a2f and restart files.'
      fileiq_tmp = TRIM(dirname) // 'restartq_' // &
           & TRIM(int_to_char(my_image_id + 1)) // '_' // &
           & TRIM(int_to_char(my_pool_id + 1)) // '.fmt'
      fila2f_tmp = TRIM(dirname) // 'a2f_tmp_' // &
           & TRIM(int_to_char(my_image_id + 1)) // '_' // &
           & TRIM(int_to_char(my_pool_id + 1))
      ! Check if files exist
      INQUIRE(FILE = fileiq_tmp, EXIST = exst)
      INQUIRE(FILE = fila2f_tmp, EXIST = exst2)
      !! Delete the files if they exist
      IF (exst) THEN
        CALL delete_if_present(fileiq_tmp, .TRUE.)
      ENDIF
      IF (exst2) THEN
        CALL delete_if_present(fila2f_tmp, .TRUE.)
      ENDIF
      WRITE(stdout, '(/5x, a/)') 'Finish removing a2f_tmp and restart_tmp files'
    ENDIF
    WRITE(stdout, '(a)') ' '
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE delete_restart_file
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE deallocate_adia_a2f(rdotk, rdotk2, w2, w0g1, etf_all, etf_loc, etf_tmp,  wkf_loc, &
                      wkf_all, wsph, epmatwef, epmatf, rrf_tmp, tmp, cufkk, cufkq, uf, cfac, cfacq)
    !-----------------------------------------------------------------------
    !! Deallocate variables at the end of fine grid interpolation
    !!
    USE kinds,            ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: rdotk(:)
    !! $r\cdot k$
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: rdotk2(:)
    !! $r\cdot k+q$
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: w2(:)
    !! Interpolated phonon frequency
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: w0g1(:, :)
    !! Delta function in energy
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: etf_all(:, :)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: etf_loc(:, :)
    !! Eigen-energies local full k-bzgrids.
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: etf_tmp(:)
    !! Temporary Eigen-energies at k-bzgrids.
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wkf_loc(:)
    !! Local weight of k-points in parallel case
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wkf_all(:)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wsph(:)
    !! phonon freq. grid
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: epmatwef(:, :, :, :)
    !! e-p matrix  in el wannier - fine Bloch phonon grid
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: epmatf(:, :, :)
    !! e-p matrix  in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: rrf_tmp(:, :, :)
    !! carrier-ionized impurity matrix in smooth Bloch basis
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: tmp(:, :, :)
    !! Temporary overlap at k and k+q
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: cufkk(:, :)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: cufkq(:, :)
    !! the same, for points k+q
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: uf(:, :)
    !! Rotation matrix for phonons
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: cfac(:)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE, INTENT(inout) :: cfacq(:)
    !! Used to store $e^{2\pi r \cdot k+q}$ exponential
    !
    ! Local variables
    INTEGER :: ierr
    !! Error status
    !
    DEALLOCATE(rdotk, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating rdotk', 1)
    DEALLOCATE(rdotk2, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating rdotk2', 1)
    DEALLOCATE(w2, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating w2', 1)
    DEALLOCATE(w0g1, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating w0g1', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating etf_all', 1)
    DEALLOCATE(etf_loc, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating etf_loc', 1)
    DEALLOCATE(etf_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating etf_tmp', 1)    
    DEALLOCATE(wkf_loc, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating wkf_loc', 1)
    DEALLOCATE(wkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating wkf_all', 1)
    DEALLOCATE(wsph, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating wsph', 1) 
    DEALLOCATE(cufkk, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating cufkk', 1)
    DEALLOCATE(cufkq, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating cufkq', 1)
    DEALLOCATE(uf, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating uf', 1)
    DEALLOCATE(cfac, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating cfac', 1)
    DEALLOCATE(cfacq, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating cfacq', 1)
    DEALLOCATE(epmatwef, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating epmatwef', 1)
    DEALLOCATE(epmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating epmatf', 1)
    DEALLOCATE(rrf_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating rrf_tmp', 1)
    DEALLOCATE(tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('calculate_a2f', 'Error deallocating tmp', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_adia_a2f
    !-----------------------------------------------------------------------
  !--------------------------------------------------------------------------
  END MODULE supercond_vertex
  !--------------------------------------------------------------------------
