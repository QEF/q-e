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
  !
  !    This subroutine tests if a xml restart file exists with the
  !    information of where the code stopped and, if appropriate the
  !    partial dynamical matrix and the partial effective charges.
  !    if (rec_code>2) done_irr, comp_irr
  !    info on calculated irreps - overrides initialization in phq_setup.
  !    The xml file is in the
  !    directory _phprefix.phsave. The xml file contains
  !    where_rec  a string with information of the point where the calculation
  !               stopped
  !    rec_code_read  where_rec     status description
  !
  !  -1000                    Nothing has been read. There is no recover file.
  !  -40         phq_setup    Only the displacements u have been read from file
  !  -30         phq_init     u and dyn(0) read from file
  !  -25                      not active yet. Restart in solve_e_fpol
  !  -20         solve_e      all previous. Stopped within solve_e. There
  !                           should be a recover file.
  !  -10         solve_e2     epsilon and zstareu are available if requested.
  !                           Stopped within solve_e2. There should be a
  !                           recover file.
  !   2          phescf       all previous, raman tenson and elop tensor are
  !                           available if required.
  !   10         solve_linter all previous. Stopped within solve linter.
  !                           There should be a recover file.
  !   20         phqscf       all previous dyn_rec(irr) and zstarue0(irr) are
  !                           available.
  !   30         dynmatrix    all previous, dyn and zstarue are available.
  !
  ! The logic of the phonon code recover is the following:
  ! The recover variable is read from input and never changed. If it is
  ! false it disables completely the recover.
  ! The control of the code is given by the arrays:
  ! comp_iq, done_iq : for each q point if it has to be calculated or
  !                    if it is already available. These are calculated
  !                    only once by check_initial_status or read from file
  !                    by the same routine.
  ! comp_irr, done_irr : for each irreducible representation if it has
  !                      to be calculated or if it is already calculated.
  !                      The latter variables are valid only for the current
  !                      q and are calculated in phq_setup and modified here
  !                      if something is on the file.
  ! epsil, done_epsil, zeu, done_zeu, zue, done_zue, lraman, done_lraman,
  ! elop, done_elop ... control the electric field calculations. These are
  ! set by prepare_q, or read from file by phq_setup.
  !
  ! The position where the code stopped is in the variable rec_code_read
  ! defined above. This variable allows to exit from a routine if the quantity
  ! calculated by this routine is already saved on file.
  ! It is the responsibility of the routine (not of the calling code)
  ! to known if it has to make the calculation or just exit because the
  ! value of rec_code_read is too high.
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
  !    rec_code_read is valid only for the first q. For the following q
  !    it is reset to -1000 in clean_pw_ph. So the recover file allows to
  !    restart only the current q. However information on other q could
  !    be available in the directory phsave, so this routine reads the
  !    appropriate files and reset comp_irr and done_irr if appropriate.
  !
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE ph_restart,    ONLY : ph_readfile
  USE control_ph,    ONLY : epsil, rec_code_read, all_done, where_rec,&
                            zeu, done_epsil, done_zeu, ext_recover, recover, &
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
  USE control_lr,    ONLY : lgamma
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
