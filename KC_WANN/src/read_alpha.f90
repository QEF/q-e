!
! Copyright (C) 2001 PWSCF group
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
  USE control_kc_wann,       ONLY : alpha_final
  USE io_global,             ONLY : stdout, ionode, ionode_id
  USE mp,                    ONLY : mp_bcast
  USE mp_global,             ONLY : intra_image_comm
  !
  ! Local Variable
  !
  IMPLICIT NONE
  !
  INTEGER j, dim
  REAL(DP) :: dum
  LOGICAL :: exst
  !
  WRITE(stdout,'(/,5x, "READING SCREENING PARAMETERS", 3x,/)') 
  !
  IF ( ionode ) THEN 
    !
    INQUIRE (FILE='file_alpharef.txt', exist=exst)
    IF (.NOT. exst) THEN 
      !
      WRITE(stdout, '(5X, "WARNING: File file_alpharef.txt NOT FOUND.")') 
      WRITE(stdout, '(5X, "WARNING: Going to set all the Screening param to 1")') 
      !
      alpha_final = 1.d0
      GOTO 101
      !
    ENDIF
    !
    OPEN (UNIT = 1001, FILE = 'file_alpharef.txt', FORM = 'formatted', STATUS = 'old' )
    READ(1001,*) dim
    !
  ENDIF
  !
  IF ( ionode ) THEN
    !
    DO j = 1, dim
      !
      READ(1001,*) dum, alpha_final(j), dum
      WRITE(stdout,'("iwann = ", 1I5, 3x, "alpha = ", 1F15.8)') j, alpha_final(j)
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
