!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE kcw_openfilq()
  !----------------------------------------------------------------------------
  !
  !! This subroutine opens all the files necessary for the LR 
  !! calculation at a given q.
  !
  USE control_flags,     ONLY : io_level
  USE units_lr,          ONLY : iuwfc, lrwfc, iudwf, lrdwf
  USE io_files,          ONLY : tmp_dir, diropn, seqopn
  USE control_kcw,       ONLY : tmp_dir_kcwq, iudvwfc, lrdvwfc, tmp_dir_save
  USE wvfct,             ONLY : nbnd, npwx
  USE io_files,          ONLY : prefix
  USE noncollin_module,  ONLY : npol
  USE buffers,           ONLY : open_buffer
  USE input_parameters,  ONLY : nk1, nk2, nk3
  USE control_kcw,       ONLY : irr_bz

  USE control_lr,        ONLY : lgamma
  USE noncollin_module,  ONLY : domag, noncolin
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
  tmp_dir=tmp_dir_kcwq
  !
  IF (lgamma.AND.nk1.eq.0.AND.nk2.eq.0.AND.nk3.eq.0) tmp_dir=tmp_dir_save
  IF ( noncolin.AND.domag ) tmp_dir=tmp_dir_kcwq
  IF( irr_bz )              tmp_dir=tmp_dir_kcwq
  !
  iuwfc = 30
  lrwfc = nbnd * npwx * npol
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  IF (.NOT.exst.AND..NOT.exst_mem) THEN
     CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  !
  ! From now on all files are written with the _ph prefix
  !
  tmp_dir=tmp_dir_kcwq
  !
  !    The file with deltaV_{bare} * psi
  !    Used in solvelinter
  !
  iudvwfc = 31
  lrdvwfc = nbnd * npwx * npol
  CALL open_buffer (iudvwfc, 'dvwfc', lrdvwfc, io_level, exst_mem, exst, tmp_dir)
  !
  !    The file with the solution delta psi
  !    Used in solve_linter
  !
  iudwf = 32
  lrdwf =  nbnd * npwx * npol
  CALL open_buffer (iudwf, 'dwf', lrdwf, io_level, exst_mem, exst, tmp_dir)
  !
  RETURN
  !
END SUBROUTINE kcw_openfilq
