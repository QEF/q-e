!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE openfilq()
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine opens all the files necessary for the phononq
  ! ... calculation.
  !
  USE kinds,           ONLY : DP
  USE control_flags,   ONLY : io_level, modenum
  USE units_ph,        ONLY : iuwfc, iudwf, iubar, iucom, iudvkb3, &
                              iudrhous, iuebar, iudrho, iudyn, iudvscf, &
                              lrwfc, lrdwf, lrbar, lrcom, lrdvkb3, &
                              lrdrhous, lrebar, lrdrho, lint3paw, iuint3paw
  USE io_files,        ONLY : tmp_dir, diropn, seqopn
  USE control_ph,      ONLY : epsil, zue, ext_recover, trans, &
                              tmp_dir_phq, start_irr, last_irr, xmldyn, &
                              all_done
  USE save_ph,         ONLY : tmp_dir_save
  USE ions_base,       ONLY : nat
  USE cell_base,       ONLY : at
  USE output,          ONLY : fildyn, fildvscf
  USE wvfct,           ONLY : nbnd, npwx
  USE fft_base,        ONLY : dfftp, dffts
  USE lsda_mod,        ONLY : nspin
  USE uspp,            ONLY : nkb, okvan
  USE uspp_param,      ONLY : nhm
  USE io_files,        ONLY : prefix
  USE noncollin_module,ONLY : npol, nspin_mag
  USE paw_variables,   ONLY : okpaw
  USE control_flags,   ONLY : twfcollect
  USE mp_bands,        ONLY : me_bgrp
  USE io_global,       ONLY : ionode,stdout
  USE buffers,         ONLY : open_buffer
  USE ramanm,          ONLY : lraman, elop, iuchf, iud2w, iuba2, lrchf, lrd2w, lrba2
  USE acfdtest,        ONLY : acfdt_is_active, acfdt_num_der
  USE input_parameters,ONLY : nk1, nk2, nk3
  USE el_phon,         ONLY : elph, elph_mat, iunwfcwann, lrwfcr
  USE dfile_star,      ONLY : dvscf_star
  USE dfile_autoname,  ONLY : dfile_name

  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : lgamma
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  ! integer variable for I/O control
  CHARACTER (len=256) :: filint, fildvscf_rot
  ! the name of the file
  LOGICAL :: exst, exst_mem
  ! logical variable to check file exists
  ! logical variable to check file exists in memory
  !
  REAL(DP) :: edum(1,1), wdum(1,1)
  INTEGER :: ndr, ierr, iq_dummy
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
 !!!!!!!!!!!!!!!!!!!!!!!! ACFDT TEST !!!!!!!!!!!!!!!!
  IF (acfdt_is_active) THEN
     ! ACFDT -test always the wfc is read/written from/to file in tmp_dir_phq
     IF (.not.acfdt_num_der)  then 
        IF (lgamma.AND.modenum==0) tmp_dir=tmp_dir_save
     ENDIF
  ELSE  
     ! this is the standard treatment
     IF (lgamma.AND.modenum==0.AND.nk1.eq.0.AND.nk2.eq.0.AND.nk3.eq.0) tmp_dir=tmp_dir_save
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!! END OF ACFDT TEST !!!!!!!!!!!!!!!!
  iuwfc = 20
  lrwfc = nbnd * npwx * npol
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  IF (.NOT.exst.AND..NOT.exst_mem.and..not.all_done) THEN
     CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  IF (elph_mat) then
     iunwfcwann=733
     lrwfcr= 2 * dffts%nr1x*dffts%nr2x*dffts%nr3x *npol
     if(ionode) then
        CALL diropn (iunwfcwann, 'wfc_r', lrwfcr, exst, dvscf_star%dir)
        IF (.NOT.exst) THEN
           CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc_r not found in Rotated_DVSCF', 1)
        END IF
     endif
  END IF
  !
  ! From now on all files are written with the _ph prefix
  !
  tmp_dir=tmp_dir_phq
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = nbnd * npwx * npol
  CALL open_buffer (iubar, 'bar', lrbar, io_level, exst_mem, exst, tmp_dir)
  IF (ext_recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.bar not found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf =  nbnd * npwx * npol
  CALL open_buffer (iudwf, 'dwf', lrdwf, io_level, exst_mem, exst, tmp_dir)
  IF (ext_recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.dwf not found', 1)
  !
  !   open a file with the static change of the charge
  !
  IF (okvan) THEN
     iudrhous = 25
     lrdrhous =  dfftp%nnr * nspin_mag
     CALL open_buffer (iudrhous, 'prd', lrdrhous, io_level, exst_mem, exst, tmp_dir)
     IF (ext_recover.AND..NOT.exst) &
        CALL errore ('openfilq','file '//trim(prefix)//'.prd not found', 1)
  ENDIF
  !
  !  Optional file(s) containing Delta\rho (opened and written in solve_e
  !  and solve_linter). Used for third-order calculations.
  !
  iudrho = 23
  lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
  !
  !   a formatted file which contains the dynamical matrix in cartesian
  !   coordinates is opened in the current directory

  !   ... by the first node only, other nodes write on unit 6 (i.e./dev/null
  !   exception: electron-phonon calculation from saved data
  !  (iudyn is read, not written, by all nodes)
  !
  IF ( ( .NOT. ionode ) .AND. (.NOT.elph.OR.trans) ) THEN
     iudyn = 6
     GOTO 400
  ENDIF

  IF (((trans.AND.(start_irr/=0.OR.last_irr/=0)).OR.elph).AND..NOT.xmldyn) THEN
     iudyn = 26
     OPEN (unit=iudyn, file=fildyn, status='unknown', err=100, iostat=ios)
100  CALL errore ('openfilq', 'opening file'//fildyn, ABS (ios) )
     REWIND (iudyn)
  ELSE
     iudyn=0
  ENDIF
  !
  !   An optional file for electron-phonon calculations containing deltaVscf
  !
400 IF (trim(fildvscf).NE.' ') THEN
     iudvscf = 27
     IF ( me_bgrp == 0 ) THEN
        IF (trim(dvscf_star%ext).NE.' ' .and. elph_mat) THEN
           fildvscf_rot = dfile_name(xq, at, TRIM(dvscf_star%ext), &
                   TRIM(dvscf_star%dir)//prefix, &
                   generate=.false., index_q=iq_dummy, equiv=.false. )
              
           WRITE(stdout,'(5x,5a)') "Opening dvscf file '",TRIM(fildvscf_rot), &
                   "' (for reading) in directory '",trim(dvscf_star%dir),"'"
              
           CALL diropn (iudvscf, fildvscf_rot, lrdrho, exst, dvscf_star%dir)
        ELSE
           CALL diropn (iudvscf, fildvscf, lrdrho, exst )
        ENDIF
        IF (okpaw) THEN
           filint=TRIM(fildvscf)//'_paw'
           lint3paw = 2 * nhm * nhm * nat * nspin_mag
           iuint3paw=34
           CALL diropn (iuint3paw, filint, lint3paw, exst)
        ENDIF
        END IF
     END IF
  !
  !    In the USPP case we need two files for the Commutator, the first is
  !    given by filbar and a second which just contains P_c x |psi>,
  !    which is required for the calculation of the Born effective carges
  !
  IF (okvan .AND. (epsil .OR. zue)) THEN
     iucom = 28
     lrcom = nbnd * npwx * npol
     CALL open_buffer (iucom, 'com', lrcom, io_level, exst_mem, exst, tmp_dir)
     IF (ext_recover.AND..NOT.exst) &
         CALL errore ('openfilq', 'file '//trim(prefix)//'.com not found', 1)
  !
  !    In the USPP case we also need a file in  order to store derivatives
  !    of kb projectors
  !
     iudvkb3 = 29
     lrdvkb3 = 2 * npwx * nkb * 3
     CALL diropn (iudvkb3, 'dvkb3', lrdvkb3, exst)
     IF (ext_recover.AND..NOT.exst) &
         CALL errore ('openfilq', 'file '//trim(prefix)//'.dvkb3 not found', 1)
  ENDIF
  IF (epsil .OR. zue) THEN
     iuebar = 30
     lrebar =  nbnd * npwx * npol
     CALL open_buffer (iuebar, 'ebar', lrebar, io_level, exst_mem, exst, tmp_dir)
     IF (ext_recover.AND..NOT.exst) &
        CALL errore ('openfilq','file '//trim(prefix)//'.ebar not found', 1)
  ENDIF
  !
  !    files used by raman calculation
  !
  IF (lraman .OR.elop) THEN
     iuchf = 31
     lrchf = 2 * nbnd * npwx * npol
     CALL diropn (iuchf, 'cwf', lrchf, exst)

     iud2w = 32
     lrd2w = 2 * nbnd * npwx * npol
     CALL diropn (iud2w, 'd2w', lrd2w, exst)

     iuba2 = 33
     lrba2 = 2 * nbnd * npwx * npol
     CALL diropn(iuba2, 'ba2', lrba2, exst)
  ENDIF

  RETURN
  !
END SUBROUTINE openfilq
