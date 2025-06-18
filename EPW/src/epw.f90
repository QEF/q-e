  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !-----------------------------------------------------------------------
  PROGRAM epw
  !-----------------------------------------------------------------------
  !! author: Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !! version: v6.0
  !! license: GNU
  !! summary: EPW main driver
  !!
  !! This is the main EPW driver which sets the phases on the wavefunctions,
  !! calls [[wann_run]] and [[elphon_shuffle_wrap]]
  !!
  !! Initial release inside Quantum ESPRESSO made on the 1st of March 2016.
  !!
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout, ionode
  USE mp_world,         ONLY : mpime
  USE mp_global,        ONLY : mp_startup, ionode_id, nimage
  USE control_flags,    ONLY : gamma_only, use_gpu, iverbosity
  USE input,            ONLY : wannierize, nqc1, nqc2, nqc3
  USE global_version,   ONLY : version_number
  USE input,            ONLY : filukk, eliashberg, ep_coupling, epwread, epbread, &
                               lcumulant, nbndsub, do_tdbe
  USE environment,      ONLY : environment_start
  USE global_var,       ONLY : elph
  USE close,            ONLY : close_final, deallocate_epw, remove_out_files
  USE cumulant,         ONLY : spectral_cumulant
  USE wannierization,   ONLY : wann_run
  USE io,               ONLY : openfilepw, loadbm
  USE stop,             ONLY : stop_epw, ephr_deallocate
  use supercond_driver, ONLY : eliashberg_eqs
  USE ep_constants,     ONLY : zero
  USE wannier,          ONLY : build_wannier
  USE check_stop,       ONLY : check_stop_init
  USE tdbe_driver,      ONLY : tdbe
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN = 12) :: code = 'EPW'
  !! Name of the code
  LOGICAL,EXTERNAL    :: check_gpu_support
  !! Name of the program
  INTEGER :: nqc
  !! Number of qpoints on the uniform grid
  INTEGER :: ierr
  !! Error index when reading/writing a file
  INTEGER :: dims
  !! Dims is either nbndsub if use_ws or 1 if not
  INTEGER :: dims2
  !! Dims is either nat if use_ws or 1 if not
  INTEGER :: nrr_k
  !! number of electronic WS points
  INTEGER :: nrr_q
  !! number of phonon WS points
  INTEGER :: nrr_g
  !! number of el-ph WS points
  INTEGER, ALLOCATABLE :: irvec_k(:, :)
  !! Integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors for electrons
  INTEGER, ALLOCATABLE :: irvec_q(:, :)
  !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
  INTEGER, ALLOCATABLE :: irvec_g(:, :)
  !! INTEGER components of the ir-th Wigner-Seitz grid point for electron-phonon
  INTEGER, ALLOCATABLE :: ndegen_k(:, :, :)
  !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
  INTEGER, ALLOCATABLE :: ndegen_q(:, :, :)
  !! Wigner-Seitz weights for the phonon grid that depend on atomic positions $R + \tau(nb) - \tau(na)$
  INTEGER, ALLOCATABLE :: ndegen_g(:, :, :)
  !! Wigner-Seitz weights for the electron-phonon grid that depend on atomic positions $R - \tau(na)$
  REAL(KIND = DP), ALLOCATABLE :: wslen_k(:)
  !! real-space length for electrons, in units of alat
  REAL(KIND = DP), ALLOCATABLE :: wslen_q(:)
  !! real-space length for phonons, in units of alat
  REAL(KIND = DP), ALLOCATABLE :: wslen_g(:)
  !! real-space length for electron-phonons, in units of alat
  REAL(KIND = DP), ALLOCATABLE :: w_centers(:, :)
  !! Wannier centers
  REAL(KIND = DP), ALLOCATABLE :: xqc(:, :)
  !! The qpoints in the uniform coarse q-point grid.
  !
  version_number = '6.0'
  !
  CALL init_clocks(.TRUE.)
  !
  CALL start_clock('EPW')
  !
  gamma_only = .FALSE.
  use_gpu = check_gpu_support()
  IF(use_gpu) Call errore('EPW', 'EPW with GPU NYI', 1)
  !
  CALL mp_startup(start_images = .TRUE.)
  !
  ! Display the logo
  IF (mpime == ionode_id) THEN
    WRITE(stdout, '(a)') "                                                                                      "
    WRITE(stdout, '(a)') "                                       ``:oss/                                        "
    WRITE(stdout, '(a)') "                           `.+s+.     .+ys--yh+     `./ss+.                           "
    WRITE(stdout, '(a)') "                          -sh//yy+`   +yy   +yy    -+h+-oyy                           "
    WRITE(stdout, '(a)') "                          -yh- .oyy/.-sh.   .syo-.:sy-  /yh                           "
    WRITE(stdout, '(a)') "                 `.-.`    `yh+   -oyyyo.     `/syys:    oys      `.`                  "
    WRITE(stdout, '(a)') "               `/+ssys+-` `sh+      `                   oys`   .:osyo`                "
    WRITE(stdout, '(a)') "               -yh- ./syyooyo`                          .sys+/oyo--yh/                "
    WRITE(stdout, '(a)') "               `yy+    .-:-.                             `-/+/:`  -sh-                "
    WRITE(stdout, '(a)') "                /yh.                                              oys                 "
    WRITE(stdout, '(a)') "          ``..---hho---------`   .---------..`      `.-----.`    -hd+---.             "
    WRITE(stdout, '(a)') "       `./osmNMMMMMMMMMMMMMMMs. +NNMMMMMMMMNNmh+.   yNMMMMMNm-  oNMMMMMNmo++:`        "
    WRITE(stdout, '(a)') "       +sy--/sdMMMhyyyyyyyNMMh- .oyNMMmyyyyyhNMMm+` -yMMMdyyo:` .oyyNMMNhs+syy`       "
    WRITE(stdout, '(a)') "       -yy/   /MMM+.`-+/``mMMy-   `mMMh:`````.dMMN:` `MMMy-`-dhhy```mMMy:``+hs        "
    WRITE(stdout, '(a)') "        -yy+` /MMMo:-mMM+`-oo/.    mMMh:     `dMMN/`  dMMm:`dMMMMy..MMMo-.+yo`        "
    WRITE(stdout, '(a)') "         .sys`/MMMMNNMMMs-         mMMmyooooymMMNo:   oMMM/sMMMMMM++MMN//oh:          "
    WRITE(stdout, '(a)') "          `sh+/MMMhyyMMMs- `-`     mMMMMMMMMMNmy+-`   -MMMhMMMsmMMmdMMd/yy+           "
    WRITE(stdout, '(a)') "    `-/+++oyy-/MMM+.`/hh/.`mNm:`   mMMd+/////:-.`      NMMMMMd/:NMMMMMy:/yyo/:.`      "
    WRITE(stdout, '(a)') "   +os+//:-..-oMMMo:--:::-/MMMo. .-mMMd+---`           hMMMMN+. oMMMMMo. `-+osyso:`   "
    WRITE(stdout, '(a)') "   syo     `mNMMMMMNNNNNNNNMMMo.oNNMMMMMNNNN:`         +MMMMs:`  dMMMN/`     ``:syo   "
    WRITE(stdout, '(a)') "   /yh`     :syyyyyyyyyyyyyyyy+.`+syyyyyyyyo:`         .oyys:`   .oyys:`        +yh   "
    WRITE(stdout, '(a)') "   -yh-        ````````````````    `````````              ``        ``          oys   "
    WRITE(stdout, '(a)') "   -+h/------------------------::::::::://////++++++++++++++++++++++///////::::/yd:   "
    WRITE(stdout, '(a)') "   shdddddddddddddddddddddddddddddhhhhhhhhyyyyyssssssssssssssssyyyyyyyhhhhhhhddddh`   "
    WRITE(stdout, '(a)') "                                                                                      "
    WRITE(stdout, '(a)') " Lee, H., Poncé, S., Bushick, K., Hajinazar, S., Lafuente-Bartolome, J.,Leveillee, J.,"
    WRITE(stdout, '(a)') "    Lian, C., Lihm, J., Macheda, F., Mori, H., Paudyal, H., Sio, W., Tiwari, S.,      "
    WRITE(stdout, '(a)') " Zacharias, M., Zhang, X., Bonini, N., Kioupakis, E., Margine, E.R., and Giustino F., "
    WRITE(stdout, '(a)') "                                                     npj Comput Mater 9, 156 (2023)   "
    WRITE(stdout, '(a)') "                                                                                      "
  ENDIF
  !
  CALL environment_start(code)
  !
  ! Do not move remove_out_files: this must be placed right after environment_start
  !
  IF (iverbosity /= -1) THEN
    CALL remove_out_files()
  ENDIF
  !
  ! Read in the input file
  !
  CALL readin()
  !
  ! Add max_seconds to smoothly stop set by user
  !
  CALL check_stop_init()
  !
  IF (epwread .AND. .NOT. epbread) THEN
    WRITE(stdout,'(a)') "                      "
    WRITE(stdout,'(a)') "     ------------------------------------------------------------------------ "
    WRITE(stdout,'(a)') "                   RESTART - RESTART - RESTART - RESTART                         "
    WRITE(stdout,'(a)') "     Restart is done without reading PWSCF save file.                  "
    WRITE(stdout,'(a)') "     Be aware that some consistency checks are therefore not done.                  "
    WRITE(stdout,'(a)') "     ------------------------------------------------------------------------ "
    WRITE(stdout,'(a)') "                      "
  ELSE
    CALL setups()
  ENDIF
  !
  !  Print run info to stdout
  !
  CALL summaries()
  !
  IF (ep_coupling) THEN
    !
    ! In case of restart with arbitrary number of cores.
    IF (epwread .AND. .NOT. epbread) THEN
      CONTINUE
    ELSE
      CALL openfilepw()
    ENDIF
    !
    CALL print_clock('EPW' )
    !
    IF (epwread .AND. .NOT. epbread) THEN
      CONTINUE
    ELSE
      CALL init(.TRUE.)
    ENDIF
    !
    CALL print_clock('EPW')
    !
    ! Generates the perturbation matrix which fixes the gauge of
    ! the calculated wavefunctions
    ! Currently, matices from setphases_wrap are identity matrices.
    ! Thus, for the moment, calling of setphases_wrap is removed.
    ! CALL setphases_wrap()
    !
    IF (wannierize) THEN
      !
      ! Create U(k, k') localization matrix
      CALL wann_run()
    ELSE
      !
      ! Read Wannier matrix from a previous run
      WRITE(stdout, '(/,5x,a,/,3a,/,5x,a,/)') REPEAT('-',67), '     Using ', &
           TRIM(filukk) , ' from disk', REPEAT('-',67)
      ! When wannierize = .false. loadbm should be called in order to load the information
      ! on band manifold determined in Wannierization step.
      CALL loadbm()
    ENDIF
    !
    IF (elph) THEN
      !
      ! Reading of dvscf in the IBZ, reconstruct the el-ph matrix element g
      ! and unfold in the full BZ. The .epb files are created for restart.
      nqc = 0
      ALLOCATE(xqc(3, nqc1 * nqc2 * nqc3), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw', 'Error allocating xqc', 1)
      ALLOCATE(w_centers(3, nbndsub), STAT = ierr)
      IF (ierr /= 0) CALL errore('epw', 'Error allocating w_centers', 1)
      nqc = 0
      xqc(:, :) = zero
      w_centers(:, :) = zero
      CALL ep_coarse_unfolding(nqc, xqc, w_centers)
      !
      ! Rotate from Bloch to Wannier gauge, remove long-range and
      ! fourier transform to real space. The .epmatwp file is created here for restart.
      CALL build_wannier(nqc, xqc, w_centers, dims, dims2, nrr_k, irvec_k, ndegen_k, wslen_k, &
                         nrr_q, irvec_q, ndegen_q, wslen_q, nrr_g, irvec_g, ndegen_g, wslen_g)
      DEALLOCATE(xqc, STAT = ierr)
      IF (ierr /= 0) CALL errore('epw', 'Error deallocating xqc', 1)
      DEALLOCATE(w_centers, STAT = ierr)
      IF (ierr /= 0) CALL errore('epw', 'Error deallocating w_centers', 1)
      !
      ! Interpolate quantities from real-space Wannier to fine Bloch k/q grids and add long-range.
      CALL use_wannier(dims, dims2, nrr_k, irvec_k, ndegen_k, wslen_k, nrr_q, irvec_q, ndegen_q, &
                       wslen_q, nrr_g, irvec_g, ndegen_g, wslen_g)
      !
      CALL ephr_deallocate(irvec_k, irvec_q, irvec_g, ndegen_k, ndegen_q, ndegen_g, &
                           wslen_k, wslen_q, wslen_g)
      !
    ENDIF
    !
    ! Cleanup of the variables
    CALL clean_pw(.FALSE.)
    CALL deallocate_epw()
    !
    ! Close the files
    CALL close_final()
    !
  ENDIF
  !
  IF (lcumulant .AND. ionode) THEN
    CALL spectral_cumulant()
  ENDIF
  !
  IF (eliashberg) THEN
    CALL eliashberg_eqs()
  ENDIF
  !
  ! ymPan: start TDBE simulation
  IF (do_tdbe) THEN
    CALL tdbe()
    CALL deallocate_epw()
  ENDIF
  !
  ! Print statistics and exit gracefully
  CALL stop_epw()
  !
  STOP
  !
  !-----------------------------------------------------------------------
  END PROGRAM epw
  !-----------------------------------------------------------------------
