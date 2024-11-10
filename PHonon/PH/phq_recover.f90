!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine phq_recover
  !-----------------------------------------------------------------------
  !! This subroutine tests if a xml restart file exists with the
  !! information of where the code stopped and, if appropriate, the
  !! partial dynamical matrix and the partial effective charges.
  !! If (\(\text{rec_code}>2\)) \(\text{done_irr}\), \(\text{comp_irr}\)
  !! info on calculated irreps - overrides initialization in 
  !! \(\texttt{phq_setup}\).  
  !! The xml file is in the directory \(\texttt{_phprefix.phsave}\).  
  !! The xml file contains \(\text{where_rec}\), a string with information of the
  !! point where the calculation stopped.
  !
  !! The \(\text{rec_code_read}\) options:
  !
  !! * -1000: nothing has been read. There is no recover file;
  !! * -40:   stops in \(\texttt{phq_setup}\). Only the displacements
  !!          u have been read from file;
  !! * -30:   stops in \(\texttt{phq_init}\). u and dyn(0) read from file;
  !! * -25:   not active yet. Restart in \(\texttt{solve_e_fpol}\);
  !! * -20:   stops in \(\texttt{solve_e}\). All previous. There should be
  !!          a recover file;
  !! * -10:   stops in \(\texttt{solve_e2}\). \(\text{epsilon} and \(\text{zstareu}\)
  !!          are available if requested. There should be a recover file.
  !! *  2:    stops in \(\texttt{phescf}\). All previous, Raman tensor and elop tensor
  !!          are available if required.
  !! *  10:   stops in \(\texttt{solve_linter}\). All previous. There should be a 
  !!          recover file.
  !! *  20:   stops in \(\texttt{phqscf}\). All previous. \(\text{dyn_rec(irr)}\) and
  !!          \(\text{zstarue0(irr)}\) are available.
  !! *  30:   \(\texttt{dynmatrix}\). All previous. \(\text{dyn}\) and \(\text{zstarue}\)
  !!          are available.
  !
  !! The logic of the phonon code recover is the following.  
  !! The recover variable is read from input and never changed. If it is
  !! false it disables completely the recover.  
  !! The control of the code is given by the arrays:
  !
  !! * \(\text{comp_iq, done_iq}\): for each q point if it has to be calculated or
  !!   if it is already available. These are calculated only once by 
  !!   check_initial_status or read from file by the same routine.
  !! * \(\text{comp_irr, done_irr}\): for each irreducible representation if it has
  !!   to be calculated or if it is already calculated. The latter variables
  !!   are valid only for the current q and are calculated in \(\texttt{phq_setup}\)
  !!   and modified here if something is on the file.
  !
  !! \(\text{epsil, done_epsil, zeu, done_zeu, zue, done_zue, lraman, done_lraman,
  !! elop, done_elop} \ldots\) control the electric field calculations. These are
  !! set by \(\texttt{prepare_q}\), or read from file by \(\texttt{phq_setup}\).
  !
  !! The position where the code stopped is in the variable rec_code_read
  !! defined above. This variable allows to exit from a routine if the quantity
  !! calculated by this routine is already saved on file.
  !! It is the responsibility of the routine (not of the calling code)
  !! to known if it has to make the calculation or just exit because the
  !! value of \(\text{rec_code_read} is too high.
  !
  ! if rec_code_read = (-25), -20, -10, 10
  !    It is expected that an unformatted recover file exists.
  !    The following data are in the
  !    unformatted file and are read by
  !    routines solve_e (-20), solve_e2 (-10), solve_linter (10):
  ! iter, dr2, convt
  !    info on status of linear-response calculation for a given irrep.
  ! dvscfin
  !    self-consistent potential for current iteration and irrep
  ! if (okpaw) dbecsum
  !    the change of the D coefficients calculated so far.
  ! if (okvan) int1, int2, int3
  !    arrays used with US potentials : int1 and int2 calculated in dvanqq,
  !    int3 calculatec in newdq (depends upon self-consistency)
  !
  !! \(\text{rec_code_read}\) is valid only for the first q. For the following q
  !! it is reset to -1000 in \(\texttt{clean_pw_ph}\). So the recover file allows to
  !! restart only the current q. However information on other q could
  !! be available in the directory \(\texttt{phsave}\), so this routine reads the
  !! appropriate files and reset \(\text{comp_irr}\) and \(\text{done_irr}\) if 
  !! appropriate.
  !
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE ph_restart,    ONLY : ph_readfile
  USE control_ph,    ONLY : epsil, all_done, zeu, done_epsil, done_zeu, ext_recover, recover, &
                            zue, trans, current_iq, low_directory_check
  USE wvfct,         ONLY : nbnd
  USE el_phon,       ONLY : el_ph_mat, el_ph_mat_rec, done_elph, elph
  USE efield_mod,    ONLY : zstarue0, zstarue0_rec
  USE partial,       ONLY : comp_irr, done_irr
  USE modes,         ONLY : nirr, npert
  USE ramanm,        ONLY : lraman, elop, done_lraman, done_elop
  USE freq_ph,       ONLY : fpol, done_fpol, done_iu, nfs
  USE grid_irr_iq,   ONLY : comp_irr_iq
  USE dynmat,        ONLY : dyn, dyn_rec

  USE qpoint,        ONLY : nksq
  USE control_lr,    ONLY : lgamma, rec_code_read, where_rec
  !
  implicit none
  !
  integer :: irr, ierr, ierr1, iu, npe, imode0
  ! counter on representations
  ! error code
  logical :: exst
  character(len=256) :: filename

  ierr=0
  IF (recover) THEN
     IF (lgamma) CALL ph_readfile('tensors', 0, 0, ierr1)
     IF (fpol.and.lgamma) THEN
        done_fpol=.TRUE.
        DO iu=1,nfs
           CALL ph_readfile('polarization', 0, iu, ierr1)
           done_fpol=done_fpol.AND.done_iu(iu)
        END DO
     ENDIF
     dyn = (0.0_DP, 0.0_DP )
     done_irr=.FALSE.
     imode0=0
     IF (elph) THEN
        el_ph_mat=(0.0_DP, 0.0_DP)
        done_elph=.FALSE.
     ENDIF
     DO irr=0, nirr
        IF (comp_irr_iq(irr,current_iq).OR..NOT.low_directory_check) THEN
           IF (trans.OR.elph) THEN
              CALL ph_readfile('data_dyn', current_iq, irr, ierr1)
              IF (ierr1 == 0) THEN
                 dyn = dyn + dyn_rec
                 IF (zue.and.irr>0) zstarue0 = zstarue0 + zstarue0_rec
              ENDIF
           END IF
           IF ( elph .and. irr > 0 ) THEN
              npe = npert(irr)
              ALLOCATE(el_ph_mat_rec(nbnd,nbnd,nksq,npe))
              CALL ph_readfile('el_phon', current_iq, irr, ierr1)
              IF (ierr1 == 0) THEN
                 el_ph_mat(:,:,:,imode0+1:imode0+npe) = &
                  el_ph_mat(:,:,:,imode0+1:imode0+npe) + el_ph_mat_rec(:,:,:,:)
              ENDIF
              DEALLOCATE(el_ph_mat_rec)
              imode0=imode0 + npe
           END IF
        ENDIF
     ENDDO

     IF (rec_code_read==-40) THEN
        WRITE( stdout, '(/,4x," Modes are read from file ")')
     ELSEIF (rec_code_read==-25) THEN
        WRITE( stdout, '(/,4x," Restart in Polarization calculation")')
     ELSEIF (rec_code_read==-20) THEN
        WRITE( stdout, '(/,4x," Restart in Electric Field calculation")')
     ELSEIF (rec_code_read==-10) then
        WRITE( stdout, '(/,4x," Restart in Raman calculation")')
     ELSEIF (rec_code_read==2) THEN
        WRITE( stdout, '(/,4x," Restart after Electric Field calculation")')
     ELSEIF (rec_code_read==10.OR.rec_code_read==20) then
        WRITE( stdout, '(/,4x," Restart in Phonon calculation")')
     ELSEIF (rec_code_read==30) then
        WRITE( stdout, '(/,4x," Restart after Phonon calculation")')
     ELSE
        call errore ('phq_recover', 'wrong restart data file', -1)
        ierr=1
     ENDIF
  ENDIF
!
  ext_recover = ext_recover .AND. ierr==0
!
! The case in which everything has been already calculated and we just
! recollect all the results must be treated in a special way (it does
! not require any initialization).
! We check here if everything has been done
!
  all_done=.true.
  DO irr = 0, nirr
     IF ( comp_irr(irr) .AND. .NOT.done_irr(irr)  ) all_done=.false.
  ENDDO
  IF (rec_code_read < 2) THEN
     IF (epsil.AND..NOT.done_epsil) all_done=.FALSE.
     IF (zeu.AND..NOT.done_zeu) all_done=.FALSE.
     IF (lraman.AND..NOT.done_lraman) all_done=.FALSE.
     IF (elop.AND..NOT.done_elop) all_done=.FALSE.
     IF (fpol.AND..NOT.done_fpol) all_done=.FALSE.
  END IF

  RETURN
END SUBROUTINE phq_recover
