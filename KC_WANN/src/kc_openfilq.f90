!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE kc_openfilq()
  !----------------------------------------------------------------------------
  !
  !! This subroutine opens all the files necessary for the LR 
  !! calculation at a given q.
  !
  USE control_flags,   ONLY : io_level
  USE units_ph,        ONLY : iudwf, iubar, &
                              iudrhous, lrdrhous,  &
                              lrdwf, lrbar
  USE units_lr,        ONLY : iuwfc, lrwfc
  USE io_files,        ONLY : tmp_dir, diropn, seqopn
  USE control_ph,      ONLY : ext_recover, &
                              tmp_dir_phq, &
                              all_done
  USE save_ph,         ONLY : tmp_dir_save
  USE wvfct,           ONLY : nbnd, npwx
  USE fft_base,        ONLY : dfftp
  USE uspp,            ONLY : okvan
  USE io_files,        ONLY : prefix
  USE noncollin_module,ONLY : npol, nspin_mag
  USE buffers,         ONLY : open_buffer
  USE input_parameters,ONLY : nk1, nk2, nk3
  USE dfile_autoname,  ONLY : dfile_name

  USE control_lr,      ONLY : lgamma
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst, exst_mem
  ! logical variable to check file exists
  ! logical variable to check file exists in memory
  !
  !
  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfilq', 'wrong prefix', 1)
  !
  !     There are six direct access files to be opened in the tmp area
  !
  !     The file with the wavefunctions. In the lgamma case reads those
  !     written by pw.x. In the other cases those calculated by ph.x
  !
  tmp_dir=tmp_dir_phq
  !
  IF (lgamma.AND.nk1.eq.0.AND.nk2.eq.0.AND.nk3.eq.0) tmp_dir=tmp_dir_save
  !
  iuwfc = 30
  lrwfc = nbnd * npwx * npol
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  IF (.NOT.exst.AND..NOT.exst_mem.and..not.all_done) THEN
     CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  !
  ! From now on all files are written with the _ph prefix
  !
  tmp_dir=tmp_dir_phq
  !
  !    The file with deltaV_{bare} * psi
  !    Used in solvelinter
  !
  iubar = 31
  lrbar = nbnd * npwx * npol
  CALL open_buffer (iubar, 'bar', lrbar, io_level, exst_mem, exst, tmp_dir)
  IF (ext_recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.bar not found', 1)
  !
  !    The file with the solution delta psi
  !    Used in solve_linter
  !
  iudwf = 32
  lrdwf =  nbnd * npwx * npol
  CALL open_buffer (iudwf, 'dwf', lrdwf, io_level, exst_mem, exst, tmp_dir)
  IF (ext_recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.dwf not found', 1)
  !
  !   open a file with the static change of the charge
  !
  IF (okvan) THEN
     iudrhous = 35
     lrdrhous =  dfftp%nnr * nspin_mag
     CALL open_buffer (iudrhous, 'prd', lrdrhous, io_level, exst_mem, exst, tmp_dir)
     IF (ext_recover.AND..NOT.exst) &
        CALL errore ('openfilq','file '//trim(prefix)//'.prd not found', 1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE kc_openfilq
