!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!-----------------------------------------------------------------------
subroutine read_alpha ( )
  !-----------------------------------------------------------------------
  !
  !! Read the screening coefficient from file. 
  !
  USE kinds,                 ONLY : DP
  USE control_kcw,           ONLY : alpha_final, tmp_dir_kcw
  USE io_global,             ONLY : stdout, ionode, ionode_id
  USE mp,                    ONLY : mp_bcast
  USE mp_global,             ONLY : intra_image_comm
  USE io_files,              ONLY : prefix
  !
  ! Local Variable
  !
  IMPLICIT NONE
  !
  INTEGER j, dim
  REAL(DP) :: dum
  LOGICAL :: exst, exst1, exst2
  CHARACTER(256) :: filename
  !
  WRITE(stdout,'(/,5x, "READING SCREENING PARAMETERS", 3x,/)') 
  !
  IF ( ionode ) THEN 
    !
    ! Read the alpha file written by kcw_screen (if it exists)
    INQUIRE (FILE=TRIM(tmp_dir_kcw)//TRIM(prefix)//'.alpha.dat', exist=exst1)
    IF (exst1) filename=TRIM(tmp_dir_kcw)//TRIM(prefix)//'.alpha.dat'
    !
    ! Try with user supplied alphafile (this has the priority)
    INQUIRE (FILE='file_alpharef.txt', exist=exst2)
    IF (exst2) filename='file_alpharef.txt'
    !
    exst = (exst1 .OR. exst2)
    !
    IF (.NOT. exst) THEN 
      !
      WRITE(stdout, '(5X, "WARNING: File with alphas  NOT FOUND.")') 
      WRITE(stdout, '(5X, "WARNING: Going to set all the Screening param to 1")') 
      !
      alpha_final = 1.d0
      GOTO 101
      !
    ENDIF
    !
    WRITE(stdout, '(5X,"INFO: alphas read from:", A)') filename
    OPEN (UNIT = 1001, FILE = filename, FORM = 'formatted', STATUS = 'old' )
    READ(1001,*) dim
    !
  ENDIF
  !
  IF ( ionode ) THEN
    !
    DO j = 1, dim
      !
      READ(1001,*) dum, alpha_final(j), dum
      WRITE(stdout,'(5X, "iwann = ", 1I5, 3x, "alpha = ", 1F15.8)') j, alpha_final(j)
      !
    ENDDO
    !
  ENDIF
  !
  CLOSE(1001)
  !
101 CONTINUE
  !
  CALL mp_bcast( dim, ionode_id, intra_image_comm )
  CALL mp_bcast( alpha_final, ionode_id, intra_image_comm )
  !
END subroutine read_alpha
