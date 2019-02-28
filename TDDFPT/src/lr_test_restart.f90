!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
LOGICAL FUNCTION test_restart(test_this)
    !
    ! This function tests whether the restart flag is applicable 
    ! Written by O. B. Malcioglu
    !
    USE lr_variables,     ONLY : n_ipol,LR_polarization,restart,bgz_suffix,eels
    USE io_files,         ONLY : prefix, tmp_dir, nd_nmbr, wfc_dir
    USE mp,               ONLY : mp_bcast, mp_barrier,mp_sum
    USE mp_world,         ONLY : world_comm
    USE io_global,        ONLY : ionode, ionode_id, stdout

    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: test_this
    CHARACTER(len=256) :: tempfile, filename, tmp_dir_saved
    LOGICAL :: exst
    CHARACTER(len=6), EXTERNAL :: int_to_char
    INTEGER :: i, temp_restart
    !
    !test_this= 1 : d0psi files
    !test_this= 2 : Lanczos restart files
    !
    temp_restart = 0
    !
    IF (.not.restart) THEN
       test_restart = .false.
       RETURN
    ENDIF
    !
    test_restart = .true.
    !
    ! d0psi files
    !
    IF (test_this == 1) THEN
       !
       ! Check for parallel i/o files that are in wfc_dir
       !
       tmp_dir_saved = tmp_dir
       !
       IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
       !
       IF ( n_ipol == 1 ) THEN
          !
          filename = trim(prefix)//'.d0psi.'//trim(int_to_char(LR_polarization))
          tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
          INQUIRE (file = tempfile, exist = exst)
          IF (.not. exst) THEN
             temp_restart = 1
          ENDIF
          !
       ELSE
          !
          DO i = 1, n_ipol
             !
             filename = trim(prefix)//'.d0psi.'//trim(int_to_char(i))
             tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
             INQUIRE (file = tempfile, exist = exst)
             IF (.not. exst) THEN
                temp_restart = 1
             ENDIF
             !
          ENDDO
          !
       ENDIF
       !
       tmp_dir = tmp_dir_saved
       !
       IF ( wfc_dir /= 'undefined' ) THEN
          !
          ! Check if these files can be read from outdir instead of wfcdir
          !
          IF ( n_ipol == 1 ) THEN 
             !
             filename = trim(prefix)//'.d0psi.'//trim(int_to_char(LR_polarization))
             tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
             INQUIRE (file = tempfile, exist = exst)
             IF (exst) THEN
                temp_restart = 0
             ENDIF
             !
          ELSE
             !
             DO i = 1, n_ipol
                !
                filename = trim(prefix)//'.d0psi.'//trim(int_to_char(i))
                tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
                INQUIRE (file = tempfile, exist = exst)
                IF (exst) THEN
                   temp_restart = 0
                ENDIF
                !
             ENDDO
             !
          ENDIF
       ENDIF 
       !
    ENDIF ! for test_this = 1
    !
    ! Lanczos restart files
    !
    IF (test_this == 2) THEN
       !
       ! Restart files are always written in outdir
       !
       IF ( n_ipol == 1 ) THEN
          filename = trim(prefix)//'.restart_lanczos.'//trim(int_to_char(LR_polarization))
          tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
       ELSE
          filename = trim(prefix)//'.restart_lanczos.'//trim(int_to_char(LR_polarization))
          tempfile = trim(tmp_dir) // trim(filename)//nd_nmbr
       ENDIF
       !
       INQUIRE (file = tempfile, exist = exst)
       !
       IF (.not. exst) THEN
          temp_restart = 1
       ENDIF
       !
       ! End of parallel file i/o
       !
       IF (eels) THEN
          filename = trim(prefix) // trim(bgz_suffix) // trim("dat")
       ELSE
         IF ( n_ipol == 1 ) THEN
            filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
         ELSE
            filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
         ENDIF
       ENDIF
       !
       tempfile = trim(tmp_dir) // trim(filename)
       !
       INQUIRE (file = tempfile, exist = exst)
       !
       IF (.not. exst) THEN
          temp_restart = 1
       ENDIF
       !
    ENDIF ! for test_this = 2
    !
#if defined(__MPI)
    CALL mp_sum(temp_restart,world_comm)
#endif
    !
    IF (temp_restart > 0 ) THEN
       !
       WRITE(stdout,'(5X,"There are missing files!")')
       !
       IF (test_this==1) WRITE(stdout,'(5X,"d0psi files can not be found, &
                                               & trying to recompansate")')
       IF (test_this==2) WRITE(stdout,'(5X,"Lanczos restart files &
                          & can not be found, starting run from scratch")')
       !
       test_restart = .false.
       !
    ENDIF
    !
    RETURN
    !
END FUNCTION test_restart
!-------------------------------------------------------------------------------
