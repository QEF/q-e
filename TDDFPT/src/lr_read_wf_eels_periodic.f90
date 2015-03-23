!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_read_wf_eels_periodic()
  !---------------------------------------------------------------------
  !
  ! This subroutine reads and stores the ground state 
  ! wavefunctions from PWscf for use in Lanczos linear 
  ! response calculation with a periodic perturbation (q=G).
  !
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, xk
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, &
                                 & diropn, tmp_dir, wfc_dir 
  USE lr_variables,         ONLY : evc0, size_evc
  USE wvfct,                ONLY : npw, igk, nbnd, g2kin, &
                                 & npwx, ecutwfc
  USE fft_base,             ONLY : dffts
  USE kinds,                ONLY : dp
  USE save_ph,              ONLY : tmp_dir_save
  USE control_ph,           ONLY : tmp_dir_phq
  USE noncollin_module,     ONLY : npol
  USE qpoint,               ONLY : nksq
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ibnd, ig, ir, itmp1,itmp2,itmp3
  LOGICAL :: exst
  CHARACTER(len=256) :: filename, tmp_dir_saved
  !
  size_evc = nksq * nbnd * npwx * npol
  !
  ! Read the unperturbed wavefunctions from the PWscf calculation.
  !
  tmp_dir = tmp_dir_save
  !
  ! This is a parallel read, done in wfc_dir
  !
  !tmp_dir_saved = tmp_dir
  !
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
  !
  CALL diropn ( iunwfc, 'wfc', nwordwfc, exst)
  !
  IF (.not.exst .and. wfc_dir == 'undefined') &
     & CALL errore('lr_read_wfc', TRIM( prefix )//'.wfc'//' not found',1)
  !
  IF (.not.exst .and. wfc_dir /= 'undefined') THEN
     !
     WRITE( stdout, '(/5x,"Attempting to read wfc from outdir instead of wfcdir")' )
     CLOSE( UNIT = iunwfc)
     tmp_dir = tmp_dir_saved
     CALL diropn ( iunwfc, 'wfc', nwordwfc, exst)
     IF (.not.exst) CALL errore('lr_read_wfc', TRIM( prefix )//'.wfc'//' not found',1)
     !
  ENDIF
  !
  DO ik = 1, nksq
     !
     CALL davcio(evc0(:,:,ik),nwordwfc,iunwfc,ik,-1)
     !
  ENDDO
  !
  CLOSE( UNIT = iunwfc)
  !
  !tmp_dir = tmp_dir_saved
  !
  tmp_dir = tmp_dir_phq
  !
  ! End of file reading
  !
  RETURN
  !
END SUBROUTINE lr_read_wf_eels_periodic
