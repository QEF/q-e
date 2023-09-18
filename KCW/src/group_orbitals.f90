!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
SUBROUTINE group_orbitals ( )
  !----------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE control_kcw,          ONLY : num_wann, l_do_alpha, group_alpha, check_spread, tmp_dir_kcw
  USE control_lr,           ONLY : lrpa
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j
  ! ... counters
  !
  COMPLEX(DP) :: sh(num_wann), sh_i, sh_j
  ! ... The self Hartree
  !
  LOGICAL :: lrpa_save, exst
  !
  l_do_alpha = .TRUE. 
  ! ... if .TRUE. the LR calculation for the orbital needs to be done 
  ! 
  IF (.NOT. check_spread) RETURN 
  ! If no check has to be performed, nothing else to do. RETURN
  !
  WRITE( stdout,'(/,5X,"INFO: Group the orbitals according to the SH ... ",/)')
  !
  INQUIRE( file=TRIM(tmp_dir_kcw)//'sh.txt', exist=exst )
  IF (exst) THEN 
    WRITE( stdout,'(/,5X,"INFO: Reading SH from file ... ",/)')
    OPEN (128, file = TRIM(tmp_dir_kcw)//'sh.txt')
    DO i = 1, num_wann
      READ(128,*) sh_i
      sh(i)=sh_i
      WRITE(stdout,'(5X, "orb, Self hartree ", 1i5, 3x, 1F10.6)') i, REAL(sh_i)
    ENDDO
    CLOSE (128)
  ELSE
    WRITE( stdout,'(/,5X,"INFO: Self-Hartree file NOT FOUND ... ")')
    WRITE( stdout,'(  5X,"INFO: Going to re-compute SH ... ",/)')
    !
    ! ... Skip xc calculation inside bare_pot 
    ! ... FIXME: ibetter and cleaner to pass "l_rpa" to "bare_pot" to decide if xc has to be added or not
    lrpa_save = lrpa
    lrpa = .true.
    !
    DO i = 1, num_wann
      ! ... Compute the Self_hartree for each Wannier 
      !
      group_alpha(i)=i
      ! ... as a default each orbital form a group (no grouping) 
      !
      sh_i = CMPLX(0.D0, 0.D0, kind= DP)
      CALL self_hartree ( i, sh_i)
      WRITE(stdout,'(5X, "orb, Self hartree ", 1i5, 3x, 1F10.6)') i, REAL(sh_i)
        sh(i) = sh_i
      !
    ENDDO
    !
    ! Restore the original value of lrpa
    lrpa = lrpa_save
    !
  ENDIF
  !
  ! Check equivalent orbital based on self-hartree
  DO i = 1, num_wann
   ! 
   l_do_alpha (i) = .TRUE.
   group_alpha(i) = i
   ! ... initializaton 
   ! 
   sh_i = CMPLX(0.D0, 0.D0, kind= DP)
   sh_i = sh(i)
   !
   DO j = 1, i-1
     ! 
     sh_j = sh(j)
     !WRITE(stdout,'(5X, "Self hartree i and j", 2i5, 3x, 2F10.6)') i,j, REAL(sh_i), REAL(sh_j)
     IF ( ABS(sh_j - sh_i) .lt. 1e-4) THEN 
       !
       l_do_alpha(i) = .FALSE. 
       group_alpha(i) = group_alpha(j)
       EXIT
       ! ... Exit as soon as I found a match
       !
     ENDIF
   ENDDO 
 ENDDO
 !
 !
 ! ... Summary of the goruping procedure 
 DO i = 1, num_wann  
   !
   IF (l_do_alpha(i) ) THEN
     WRITE(stdout,'(8X, "iwann=", i5, 3X, "DO_LR =", L)') i, l_do_alpha (i)
   ELSE
     WRITE(stdout,'(8X, "iwann=", i5, 3X, "DO_LR =", L, 3x, "--> " i5)') i, l_do_alpha(i), group_alpha(i)
   ENDIF
   !
 ENDDO
 !
 WRITE( stdout,'(/,5X,"INFO: Group the orbitals according to the SH ... DONE ")')
 !
END SUBROUTINE group_orbitals

