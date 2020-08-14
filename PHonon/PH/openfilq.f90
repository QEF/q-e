!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
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
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : io_level, modenum
  USE units_ph,         ONLY : iudwf, iubar, iucom, iudvkb3, &
                              iudrhous, iuebar, iudrho, iudyn, iudvscf, &
                              lrdwf, lrbar, lrcom, lrdvkb3, &
                              lrdrhous, lrebar, lrdrho, lint3paw, iuint3paw, &
                              iundnsscf, iudvpsi, lrdvpsi, iugauge
  USE units_lr,         ONLY : iuwfc, lrwfc
  USE io_files,         ONLY : tmp_dir, diropn, seqopn, nwordwfcU
  USE control_ph,       ONLY : epsil, zue, ext_recover, trans, &
                              tmp_dir_phq, start_irr, last_irr, xmldyn, &
                              all_done, newgrid
  USE save_ph,         ONLY : tmp_dir_save
  USE ions_base,       ONLY : nat
  USE cell_base,       ONLY : at
  USE output,          ONLY : fildyn, fildvscf
  USE wvfct,           ONLY : nbnd, npwx
  USE fft_base,        ONLY : dfftp, dffts
  USE lsda_mod,        ONLY : nspin, lsda
  USE uspp,            ONLY : nkb, okvan
  USE uspp_param,      ONLY : nhm
  USE io_files,        ONLY : prefix
  USE noncollin_module,ONLY : npol, nspin_mag, noncolin
  USE paw_variables,   ONLY : okpaw
  USE mp_bands,        ONLY : me_bgrp
  USE spin_orb,        ONLY : domag
  USE io_global,       ONLY : ionode,stdout
  USE buffers,         ONLY : open_buffer, close_buffer
  USE ramanm,          ONLY : lraman, elop, iuchf, iud2w, iuba2, lrchf, lrd2w, lrba2
  USE acfdtest,        ONLY : acfdt_is_active, acfdt_num_der
  USE el_phon,         ONLY : elph, elph_mat, iunwfcwann, lrwfcr
  USE dfile_star,      ONLY : dvscf_star
  USE dfile_autoname,  ONLY : dfile_name
  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : lgamma
  USE units_lr,        ONLY : iuatwfc, iuatswfc
  USE modes,           ONLY : nmodes
  USE ldaU,            ONLY : lda_plus_u, Hubbard_lmax, nwfcU
  USE ldaU_ph,         ONLY : dnsscf_all_modes
  USE mp_pools,        ONLY : me_pool, root_pool
  USE dvscf_interpolate, ONLY : ldvscf_interpolate, nrbase, nrlocal, &
                                wpot_dir, iunwpot, lrwpot
  USE ahc,              ONLY : elph_ahc, ahc_nbnd_gauge
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  ! integer variable for I/O control
  CHARACTER (len=256) :: filint, fildvscf_rot, filwpot
  ! the name of the file
  INTEGER :: ir, irlocal
  !! Real space unit cell index
  INTEGER :: unf_lrwpot, direct_io_factor
  !! record length for opening wpot file
  REAL(DP) :: dummy
  !! dummy variable for calculating direct_io_factor
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
     IF (lgamma.AND.modenum==0.AND..NOT.newgrid ) tmp_dir=tmp_dir_save
     ! FIXME: why this case?
     IF ( noncolin.AND.domag ) tmp_dir=tmp_dir_phq
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!! END OF ACFDT TEST !!!!!!!!!!!!!!!!
  iuwfc = 20
  lrwfc = nbnd * npwx * npol
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  IF (.NOT.exst.AND..NOT.exst_mem.and..not.all_done) THEN
     CALL close_buffer(iuwfc, 'delete') 
     !FIXME Dirty fix for obscure case
     tmp_dir = tmp_dir_phq
     CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
     IF (.NOT.exst.AND..NOT.exst_mem) CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc not found', 1)
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
  !
  ! Files needed for DFPT+U calculation
  !
  IF (lda_plus_u) THEN   
     !
     nwordwfcU = npwx * nwfcU * npol
     !
     ! The unit iuatwfc contains atomic wfcs at k and k+q
     !    
     iuatwfc = 34
     CALL open_buffer (iuatwfc, 'atwfc', nwordwfcU, io_level, exst_mem, exst, tmp_dir)
     !
     ! The unit iuatswfc contains atomic wfcs * S at k and k+q
     !    
     iuatswfc = 35
     CALL open_buffer (iuatswfc, 'satwfc', nwordwfcU, io_level, exst_mem, exst, tmp_dir)
     !
     ! Open a file to write dnsscf_all_modes
     !
     iundnsscf = 36
     IF (trans .OR. elph) THEN
        !
        ! Open a file
        ! Note: if trans=.true. then dnsscf_all_modes will be written to file (see phqscf)
        !
        IF (ionode) CALL seqopn (iundnsscf, 'dnsscf', 'formatted', exst)
        !
        ! If elph=.true. and trans=.true., then dnsscf (dnsscf_all_modes) is computed and 
        ! kept in memory and hence we can directly use it in elphel.
        ! If elph=.true. and trans=.false. (i.e. phonons must be computed in advance),
        ! then we read dnsscf (dnsscf_all_modes) from file.
        !
        IF (elph .AND. .NOT.trans) THEN ! In this case we read dnsscf_all_modes
           !
           IF (ionode) THEN
             !
             IF (.NOT.exst) &
                CALL errore ('openfilq', 'dnsscf file not found, necessary for el-ph calculation, stopping ', 1)
             ! 
             IF (exst) THEN
                !
                ! Here we allocate and deallocate dnsscf_all_modes just to check that it is read from file properly.
                ! In elphel, dnsscf_all_modes will be allocated and read again (for production purposes).
                !
                ALLOCATE (dnsscf_all_modes (2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, nmodes))
                READ(iundnsscf,*,iostat=ios) dnsscf_all_modes
                REWIND(iundnsscf)
                !
                IF (ios.NE.0) &
                   CALL errore ('openfilq', 'dnsscf file corrupted, necessary for el-ph calculation, stopping ', 1)
                ! 
                IF (ios==0) &
                   WRITE( stdout,*) 'THE DNSSCF MATRIX WAS CORRECTLY READ FROM FILE, NECESSARY FOR ELPH+U'
                !
                DEALLOCATE(dnsscf_all_modes)
                !
             ENDIF
             !
           ENDIF
           !
        ENDIF
        !
     ENDIF
     !
  ENDIF
  !
  ! Files needed for dvscf interpolation
  !
  ! Here, root of each pool read different files. Subroutine diropn is not used
  ! because diropn adds processor id at the end of filename.
  !
  IF (ldvscf_interpolate) THEN
    !
    lrwpot = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
    ! Need to multiply direct_io_factor: See diropn in Modules/io_files.f90
    INQUIRE (IOLENGTH=direct_io_factor) dummy
    unf_lrwpot = direct_io_factor * INT(lrwpot, KIND=KIND(unf_lrwpot))
    !
    ! w_pot files are read by the root of each pool
    !
    IF (me_pool == root_pool) THEN
      !
      DO irlocal = 1, nrlocal
        !
        ir = irlocal + nrbase
        !
#if defined(__MPI)
        WRITE(filwpot, '(a,I0,a)') TRIM(wpot_dir) // TRIM(prefix) // '.' &
                                   // 'wpot.irc', ir, '.dat1'
#else
        WRITE(filwpot, '(a,I0,a)') TRIM(wpot_dir) // TRIM(prefix) // '.' &
                                   // 'wpot.irc', ir, '.dat'
#endif
        !
        ! FIXME: better way to set units?
        !
        iunwpot(irlocal) = 1000 + irlocal
        !
        OPEN(iunwpot(irlocal), FILE=TRIM(filwpot), RECL=unf_lrwpot, &
          FORM='unformatted', STATUS='unknown', ACCESS='direct', IOSTAT=ios)
        !
        IF (ios /= 0) CALL errore('openfilq', &
          'problem opening w_pot file ' // TRIM(filwpot), 1)
        !
      ENDDO ! irlocal
      !
    ENDIF ! root_pool
    !
  ENDIF ! ldvscf_interpolate
  !
  ! elph_ahc
  !
  IF (elph_ahc) THEN
    !
    ! File containing delta V_{SCF} * psi
    !
    iudvpsi = 41
    lrdvpsi = ahc_nbnd_gauge * npwx * npol
    CALL open_buffer(iudvpsi, 'dvpsi', lrdvpsi, io_level, exst_mem, exst, tmp_dir)
    IF (ext_recover .AND. .NOT. exst) &
       CALL errore ('openfilq', 'file '//trim(prefix)//'.dvpsi not found', 1)
    !
    ! File for computing gauge
    iugauge = 42
    IF (me_pool == root_pool) THEN
      CALL seqopn(iugauge, 'wfcgauge', 'unformatted', exst, tmp_dir_save)
    ENDIF
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE openfilq
