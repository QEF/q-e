!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  SUBROUTINE write_qplot_data(auxdyn)
!
!  This routine writes on output the phonon frequencies in a format that
!  can be read by the program plotband. It is used only if the qplot
!  option is used in the phonon code. It writes in a file called fildyn.plot
!
  USE kinds,      ONLY : DP
  USE constants,  ONLY : ry_to_cmm1
  USE ions_base,  ONLY : nat
  USE disp,       ONLY : nqs, omega_disp, x_q, done_iq
  USE control_ph, ONLY : qplot
  USE el_phon,    ONLY : elph_simple, gamma_disp, el_ph_nsigma
  USE mp_images,  ONLY : nimage
  USE output,     ONLY : fildyn
  USE io_global,  ONLY : ionode

  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(IN) :: auxdyn
  CHARACTER(LEN=256) :: filename
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  REAL(DP) :: w1(3*nat)
  INTEGER :: iunit 
  INTEGER :: n, i, isig, iq

  IF ( .NOT. qplot ) CALL errore('write_qplot_data', &
                                 'called in the wrong case', 1)
!
!  Do not write if there are more than one image
!
  IF (nimage > 1) RETURN
!
!  This routine writes the files only when all q points have been done.
!
  DO iq=1,nqs
     IF (.NOT.done_iq(iq)) RETURN
  ENDDO

  IF (ionode) THEN
     iunit = find_free_unit()
     filename = TRIM(auxdyn) // '.freq'
     OPEN (unit=iunit,file=TRIM(filename),status='unknown',form='formatted')
     WRITE(iunit, '(" &plot nbnd=",i4,", nks=",i4," /")') 3*nat, nqs
     DO n=1, nqs
        WRITE(iunit,'(10x,3f10.6)')  x_q(1,n), x_q(2,n), x_q(3,n)
        DO i=1, 3*nat
           w1(i) = SQRT (ABS (omega_disp (i,n)) ) * ry_to_cmm1
           IF ( omega_disp(i,n) < 0.d0) w1(i) = - w1(i)
        ENDDO
        WRITE(iunit,'(6f10.4)') (w1(i), i=1,3*nat)
     END DO
     CLOSE(unit=iunit)
     IF (elph_simple) THEN
        DO isig=1, el_ph_nsigma
           filename = TRIM(auxdyn) // '.elph.'// int_to_char(isig)
           OPEN (unit=iunit,file=TRIM(filename),status='unknown', &
                                                form='formatted')
           WRITE(iunit, '(" &plot nbnd=",i4,", nks=",i4," /")') 3*nat, nqs
           DO n=1, nqs
              WRITE(iunit,'(10x,3f10.6)')  x_q(1,n), x_q(2,n), x_q(3,n)
              WRITE(iunit,'(6f10.4)') (gamma_disp(i,isig,n), i=1,3*nat)
           END DO
           CLOSE(unit=iunit)
        END DO
     ENDIF 
  ENDIF

  RETURN
  END SUBROUTINE write_qplot_data
