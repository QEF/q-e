!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE check_initial_status(auxdyn)
  !-----------------------------------------------------------------------
  !
  ! This routine checks the initial status of the phonon run and prepares
  ! the control of the dispersion calculation. The grid is determined 
  ! by the following variables:
  ! nqs : the number of q points
  ! x_q : the coordinates of the q points
  ! comp_iq : =1 if this q point is calculated in this run, 0 otherwise
  ! done_iq : =1 if already calculated, 0 otherwise
  ! rep_iq  : for each q point how many irreducible representations
  ! done_rep_iq : =1 if the representation has been already calculated
  ! nfs : the number of imaginary frequencies
  ! The last four variables are calculated only at gamma if fpol is .true.
  ! If recover is true this routine checks also that the control parameters
  ! read in input match the values given in the recover file.
  ! Finally it sets the variable:
  ! done_bands : if true the bands have been already calculated for this q
  !
  USE io_global,       ONLY : stdout
  USE control_flags,   ONLY : modenum
  USE ions_base,       ONLY : nat
  USE io_files,        ONLY : tmp_dir
  USE lsda_mod,        ONLY : nspin
  USE scf,             ONLY : rho
  USE disp,            ONLY : nqs, x_q, comp_iq
  USE qpoint,          ONLY : xq
  USE output,          ONLY : fildyn
  USE control_ph,      ONLY : ldisp, recover, done_bands,  &
                              start_q, last_q, current_iq, tmp_dir_ph, lgamma, &
                              ext_recover, ext_restart
  USE ph_restart,      ONLY : ph_readfile, check_status_run, init_status_run
  USE start_k,         ONLY : nks_start
  USE save_ph,         ONLY : save_ph_input_variables
  USE io_rho_xml,      ONLY : write_rho
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256) :: auxdyn
  INTEGER :: iq, iq_start, ierr
  INTEGER :: iu
  !
  ! Initialize local variables
  !
  tmp_dir=tmp_dir_ph
  !
  ! ... Checking the status of the calculation
  !
  IF (recover) THEN
!
!  check if a recover file exists. In this case the first q point is
!  the current one.
!
     IF (.NOT.ext_recover.AND..NOT.ext_restart) THEN
        iq_start=start_q
        done_bands=.FALSE.
     ELSE
        iq_start=current_iq
     ENDIF
!
!  check which representation files are available on the disk and 
!  sets which q points and representations have been already calculated
!
     CALL check_status_run()
!
! write the information on output
!
     IF (last_q<1.OR.last_q>nqs) last_q=nqs
     IF (iq_start<=last_q.AND.iq_start>0) THEN
        WRITE(stdout, &
            '(5x,i4," /",i4," q-points for this run, from", i3,&
               & " to", i3,":")') last_q-iq_start+1, nqs, iq_start, last_q
        WRITE(stdout, '(5x,"  N       xq(1)         xq(2)         xq(3) " )')
        DO iq = 1, nqs
           WRITE(stdout, '(5x,i3, 3f14.9,l6)') iq, x_q(1,iq), x_q(2,iq), &
                             x_q(3,iq)
        END DO
        WRITE(stdout, *)
     ELSEIF (iq_start>last_q) THEN
        WRITE(stdout, &
            '(5x,"Starting q",i4," larger than total number of q points", i4, &
               & " or of last q ", i3)') iq_start, nqs, last_q
     ELSEIF (iq_start<0) THEN
        CALL errore('check_initial_status','wrong iq_start',1)
     ENDIF
  ENDIF
  !
  !  Create a new directory where the ph variables are saved and copy
  !  the charge density there.
  !
  IF (ldisp.OR..NOT.lgamma.OR.modenum/=0) CALL write_rho( rho, nspin )
  CALL save_ph_input_variables()
  !
  IF (.NOT.recover) THEN
     !
     ! recover file not found or not looked for
     !
     done_bands=.FALSE.
     iq_start=start_q
     IF (ldisp) THEN
        !
        ! ... Calculate the q-points for the dispersion
        !
        CALL q_points()
        IF (last_q<1.or.last_q>nqs) last_q=nqs
        !
     ELSE
        !
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        x_q(:,1)=xq(:)
        !
     END IF
     !
     ! This routine initialize the grid control in order to
     ! calculate all q and all representations.
     !
     CALL init_status_run()
     !
  END IF
  !
  !  Set the q points to calculate. If there is the recover file, start from
  !  the q point of the recover file.
  !
  ALLOCATE(comp_iq(nqs))
  comp_iq=0
  DO iq=iq_start,last_q
     comp_iq(iq)=1
  ENDDO
  !
  auxdyn = fildyn
  !
  IF (nks_start==0) CALL errore('check_initial_status','wrong starting k',1)
  !
  RETURN
  END SUBROUTINE check_initial_status

