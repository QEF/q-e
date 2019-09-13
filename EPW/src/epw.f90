  !                                                                            
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
  !! version: v5.2
  !! license: GNU
  !! summary: EPW main driver 
  !!  
  !! This is the main EPW driver which sets the phases on the wavefunctions,
  !! calls [[wann_run]] and [[elphon_shuffle_wrap]]
  !!
  USE io_global,       ONLY : stdout, ionode
  USE mp,              ONLY : mp_bcast, mp_barrier
  USE mp_world,        ONLY : mpime  
  USE mp_global,       ONLY : mp_startup, ionode_id, mp_global_end
  USE control_flags,   ONLY : gamma_only
  USE control_epw,     ONLY : wannierize
  USE global_version,  ONLY : version_number
  USE epwcom,          ONLY : filukk, eliashberg, ep_coupling, epwread, epbread, cumulant
  USE environment,     ONLY : environment_start
  USE elph2,           ONLY : elph 
  USE close_epw,       ONLY : close_final, deallocate_epw
  USE cum_mod,         ONLY : spectral_cumulant
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN = 12) :: code = 'EPW'
  !! Name of the program
  !
  version_number = '5.2.0'
  !
  CALL init_clocks(.TRUE.)
  !
  CALL start_clock('EPW')
  !
  gamma_only = .FALSE.
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
    WRITE(stdout, '(a)') "  S. Ponce, E. R. Margine, C. Verdi, and F. Giustino,                                 "
    WRITE(stdout, '(a)') "                                                Comput. Phys. Commun. 209, 116 (2016) "
    WRITE(stdout, '(a)') "                                                                                      "
  ENDIF
  !
  CALL environment_start(code)
  !
  ! Read in the input file
  !
  CALL epw_readin
  !
  IF (epwread .AND. .NOT. epbread) THEN
    WRITE(stdout,'(a)') "                      "
    WRITE(stdout,'(a)') "     ------------------------------------------------------------------------ "
    WRITE(stdout,'(a)') "                   RESTART - RESTART - RESTART - RESTART                         "
    WRITE(stdout,'(a)') "     Restart is done without reading PWSCF save file.                  "
    WRITE(stdout,'(a)') "     Be aware that some consistency checks are therefore not done.                  "
    WRITE(stdout,'(a)') "     ------------------------------------------------------------------------ "
    WRITE(stdout,'(a)') "                      "
    CALL epw_setup_restart
  ELSE
    CALL epw_setup
  ENDIF
  !
  !  Print run info to stdout
  !
  CALL epw_summary
  !
  IF (ep_coupling) THEN 
    !
    ! In case of restart with arbitrary number of cores.
    IF (epwread .AND. .NOT. epbread) THEN
      CONTINUE
    ELSE 
      CALL openfilepw
    ENDIF
    !
    CALL print_clock('EPW' )
    !
    IF (epwread .AND. .NOT. epbread) THEN
      CONTINUE      
    ELSE
      CALL epw_init(.TRUE.)
    ENDIF
    !
    CALL print_clock('EPW')
    !
    ! Generates the perturbation matrix which fixes the gauge of 
    ! the calculated wavefunctions
    CALL setphases_wrap
    !
    IF (wannierize) THEN
      !
      ! Create U(k, k') localization matrix 
      CALL wann_run
    ELSE
      !
      ! Read Wannier matrix from a previous run
      WRITE(stdout, '(/,5x,a,/,3a,/,5x,a,/)') REPEAT('-',67), '     Using ', &
           TRIM(filukk) , ' from disk', REPEAT('-',67) 
    ENDIF
    !
    IF (elph) THEN
      !
      CALL elphon_shuffle_wrap()
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
  IF (cumulant .AND. ionode) THEN
    CALL spectral_cumulant()
  ENDIF
  !
  IF (eliashberg) THEN
    CALL eliashberg_eqs()
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
