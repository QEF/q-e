!-----------------------------------------------------------------------
SUBROUTINE lr_read_wf()
  !---------------------------------------------------------------------
  ! ... reads in and stores the ground state wavefunctions
  ! ... for use in Lanczos linear response calculation
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Xiaochuan Ge (2013) to fix some bugs of virt_read and include  davidson

  USE kinds,                ONLY : dp
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, xk
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, diropn,&
                                 & tmp_dir, wfc_dir 
  USE lr_variables,         ONLY : evc0, sevc0 ,revc0, evc0_virt,&
                                 & sevc0_virt, nbnd_total, becp1_virt, &
                                 & becp1_c_virt, no_hxc
  USE realus,               ONLY : igk_k,npw_k
  USE lr_variables,         ONLY : becp1, becp1_c, test_case_no,&
                                 & size_evc, project
  USE wvfct,                ONLY : npw, igk, nbnd, g2kin, npwx, ecutwfc
  USE control_flags,        ONLY : gamma_only,io_level
  USE gvecs,                ONLY : nls, nlsm
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE becmod,               ONLY : bec_type, becp, calbec
  USE realus,               ONLY : real_space, fft_orbital_gamma,&
                                 & initialisation_level,&
                                 & bfft_orbital_gamma, calbec_rs_gamma,&
                                 & add_vuspsir_gamma, v_loc_psir,&
                                 & s_psir_gamma, real_space_debug  
  USE lr_variables,         ONLY : lr_verbosity, lr_exx, davidson
  USE buffers,              ONLY : get_buffer
  USE exx,                  ONLY : exx_grid_init, exx_div_check, exx_restart
  USE funct,                ONLY : dft_is_hybrid
  USE lr_exx_kernel,        ONLY : lr_exx_revc0_init, lr_exx_alloc
  USE wavefunctions_module, ONLY : evc
  USE buffers,              ONLY : open_buffer
  !
  IMPLICIT NONE
  !
  !
  !
  ! local variables
  INTEGER :: ik, ibnd, ig, itmp1,itmp2,itmp3
  LOGICAL :: exst
  CHARACTER(len=256) :: filename, tmp_dir_saved
  !OBM debug
  real(kind=dp) :: obm_debug
  COMPLEX(kind=dp),EXTERNAL :: lr_dot
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_read_wf>")')
  ENDIF
  call start_clock("read_wf")
  !
  IF ((nbnd_total>nbnd .and. davidson) .OR. project) THEN
     CALL virt_read()
  ELSE
     CALL normal_read()
  ENDIF


  WRITE(stdout,'(5x,"Finished reading wfc.")')

  evc(:,:)=evc0(:,:,1)
  IF ( dft_is_hybrid() ) THEN
     CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
     CALL exx_grid_init()
     CALL exx_div_check()
     CALL exx_restart(.true.)
     IF(.NOT. no_hxc) THEN
        lr_exx =.TRUE.
        CALL lr_exx_alloc()
        DO ik=1,nks
           CALL lr_exx_revc0_init(evc0,ik)
        ENDDO
     ENDIF
  ENDIF

  WRITE(stdout,'(5x,"Finished exx setting.")')
  !

  call stop_clock("read_wf")
  RETURN
  !
  CONTAINS
!-----------------------------------------------------------------------
    SUBROUTINE normal_read()
      !
      ! The usual way of reading wavefunctions
      !
      USE lr_variables,             ONLY : check_all_bands_gamma, &
                                          & check_density_gamma, &
                                          & check_vector_gamma, tg_revc0
      USE wavefunctions_module,     ONLY : psic
      USE realus,                   ONLY : tg_psic
      USE mp_global,                ONLY : me_bgrp
      !
      !
      IMPLICIT NONE
      !
      LOGICAL :: use_tg
      INTEGER :: v_siz, incr, ioff, j

      WRITE( stdout, '(/5x,"Normal read")' )

      use_tg=dffts%have_task_groups
      nwordwfc = nbnd * npwx
      size_evc=npwx*nbnd*nks
      incr = 2
      !
      !   Read in the ground state wavefunctions
      !   This is a parallel read, done in wfc_dir
      !
      tmp_dir_saved = tmp_dir
      !IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
      wfc_dir = tmp_dir
      CALL diropn ( iunwfc, 'wfc', 2*nwordwfc, exst)
      !
      IF (.NOT.exst .AND. wfc_dir == 'undefined') CALL errore(&
           &'lr_read_wfc', TRIM( prefix )//'.wfc'//' not found',1) 
      !
      IF (.NOT.exst .AND. wfc_dir /= 'undefined') THEN
         !
         WRITE( stdout, '(/5x,"Attempting to read wfc from outdir &
              &instead of wfcdir")' ) 
         CLOSE( UNIT = iunwfc)
         !
         tmp_dir = tmp_dir_saved
         !
         CALL diropn ( iunwfc, 'wfc', 2*nwordwfc, exst)
         !
         IF (.NOT.exst) CALL errore('lr_read_wfc', &
              &TRIM( prefix )//'.wfc'//' not found',1)
      ENDIF
      !
      IF (gamma_only) THEN
         WRITE( stdout, '(/5x,"Gamma point algorithm")' )
      ELSE
         WRITE( stdout, '(/5x,"Generalised algorithm !warning")' )
      ENDIF
      !
      DO ik=1,nks
         !
         IF (.NOT. real_space_debug > 0 ) THEN 
            !else done in init_realspace realus
            CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw,&
                 & igk, g2kin ) 
            !
            npw_k(ik) = npw
            !
            igk_k(:,ik) = igk(:)
         ENDIF
         !
         !   Read in the ground state wavefunctions
         !   This is a parallel read, done in wfc_dir
         CALL davcio(evc0(:,:,ik),2*nwordwfc,iunwfc,ik,-1)
         !
      ENDDO
      !
      !
      CLOSE( UNIT = iunwfc)
      !
      ! End of file reading
      !
      tmp_dir = tmp_dir_saved
      !
      !
      ! vkb * evc0 and initialization of sevc0
      !
      !
      IF ( okvan ) THEN
         !
         IF ( gamma_only ) THEN
            !
            ! Following line is to be removed when real space
            ! implementation is complete
            CALL init_us_2(npw,igk_k(:,1),xk(1,1),vkb)
            !
            IF (real_space_debug>0) THEN
               !
               !
               DO ibnd=1,nbnd,2
                  CALL fft_orbital_gamma(evc0(:,:,1),ibnd,nbnd)
                  CALL calbec_rs_gamma(ibnd,nbnd,becp1)
                  becp%r(:,ibnd)=becp1(:,ibnd)
                  IF (ibnd + 1 <= nbnd) becp%r(:,ibnd+1)=becp1(:,ibnd+1)
                  CALL s_psir_gamma(ibnd,nbnd)
                  CALL bfft_orbital_gamma(sevc0(:,:,1),ibnd,nbnd)
               ENDDO
               !
            ELSE
               !
               CALL calbec(npw_k(1),vkb,evc0(:,:,1),becp1)
               !
               becp%r=becp1
               !
               CALL s_psi(npwx, npw_k(1), nbnd, evc0(:,:,1), sevc0(:,:,1))
               !
            ENDIF
         ELSE
            !
            ! K point generalized stuff starts here
            !
            DO ik=1,nks
               !
               CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
               !
               CALL calbec(npw_k(ik),vkb,evc0(:,:,ik),becp1_c(:,:,ik))
               !
               becp%k=becp1_c(:,:,ik)
               !
               CALL s_psi (npwx, npw_k(ik), nbnd, evc0(:,:,ik), &
                    & sevc0(:,:,ik)) 
               !
            ENDDO
            !
         ENDIF
         !
      ELSE
         !
         sevc0=evc0
         !
      ENDIF
      !
      !
      ! Inverse fourier transform of evc0
      !
      !
      IF( dffts%have_task_groups ) THEN
         !
         v_siz =  dffts%tg_nnr * dffts%nogrp
         !
         incr = 2 * dffts%nogrp
         !
         tg_revc0=(0.0d0,0.0d0)
         !
      ELSE
         !
         revc0=(0.0d0,0.0d0)
         !
      ENDIF
      !
      IF ( gamma_only ) THEN
         !
         !
         DO ibnd=1,nbnd,incr
            !
            CALL fft_orbital_gamma ( evc0(:,:,1), ibnd, nbnd)
            !
            IF (dffts%have_task_groups) THEN               
               !
               DO j=1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
                  tg_revc0(j,ibnd,1) = tg_psic(j)
               ENDDO
               !
            ELSE
               !
               revc0(1:dffts%nnr,ibnd,1) = psic(1:dffts%nnr)
               !
            ENDIF
            !
         ENDDO
         !
      ELSE
         !
         DO ik=1,nks
            !
            DO ibnd=1,nbnd
               !
               DO ig=1,npw_k(ik)
                  !
                  revc0(nls(igk_k(ig,ik)),ibnd,ik)=&
                       &evc0(ig,ibnd,ik)
                  !
               ENDDO
               !
               dffts%have_task_groups=.false.
               CALL invfft ('Wave', revc0(:,ibnd,ik), dffts)
               dffts%have_task_groups=use_tg
               !
            ENDDO
            !
         ENDDO
         !
      ENDIF
      ! OBM - Last minute check for real space implementation,
      IF ( real_space_debug > 0 .AND. .NOT. gamma_only ) &
           CALL errore( ' iosys ', ' Linear response calculation ' // &
           & 'real space algorithms with k-points not implemented', 1 )
      !
    END SUBROUTINE normal_read
!-----------------------------------------------------------------------
    SUBROUTINE virt_read()
      ! The modifications to read also the virtual orbitals

      USE control_ph,            ONLY : nbnd_occ
      USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type
      USE lr_variables,          ONLY : davidson
      USE io_global,            ONLY : stdout
      !
      !
      IMPLICIT NONE
      !
      COMPLEX(kind=dp), ALLOCATABLE :: evc_all(:,:,:)
      COMPLEX(kind=dp), ALLOCATABLE :: sevc_all(:,:,:)
      REAL(kind=dp),    ALLOCATABLE :: becp1_all(:,:)
      COMPLEX(kind=dp), ALLOCATABLE :: becp1_c_all(:,:,:)
      COMPLEX(kind=dp), ALLOCATABLE :: revc_all(:,:,:)
      !
      ! Check for task groups
      WRITE( stdout, '(/5x,"Virt read")' )
      IF (dffts%have_task_groups) CALL errore ( 'virt_read', 'Task &
           &groups not supported when there are virtual states in the &
           &input.', 1 )
      !
      ! First pretend everything is normal
      nbnd=nbnd_total
      !
      ALLOCATE(revc_all(dffts%nnr,nbnd,nks))
      ALLOCATE(evc_all(npwx,nbnd,nks))
      ALLOCATE(sevc_all(npwx,nbnd,nks))
      ALLOCATE(evc0_virt(npwx,(nbnd-nbnd_occ(1)),nks))
      allocate(sevc0_virt(npwx,(nbnd-nbnd_occ(1)),nks))
      !
      ! XC: It's kind of messy to deallocate and alloacte again becp; This is mainly due to the
      ! change of nbnd. Later this operation need to be done again when we change nbnd back to nbnd_occ
      ! Again, I suggest please, in the future try the best not to change the meaning of global variables
      CALL deallocate_bec_type(becp)
      CALL allocate_bec_type ( nkb, nbnd, becp )

      IF (nkb > 0) THEN
         IF(gamma_only) THEN
            ALLOCATE(becp1_all(nkb,nbnd))
            becp1_all(:,:)=0.0d0
         ELSE
            ALLOCATE(becp1_c_all(nkb,nbnd,nks))
            becp1_c_all(:,:,:)=(0.0d0,0.0d0)
         ENDIF
      ENDIF

      nwordwfc = nbnd * npwx
      size_evc=npwx*nbnd_occ(1)*nks

      !   Read in the ground state wavefunctions
      !   This is a parallel read, done in wfc_dir
      tmp_dir_saved = tmp_dir
      !IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
      wfc_dir = tmp_dir
      CALL diropn ( iunwfc, 'wfc', 2*nwordwfc, exst)
      IF (.NOT.exst .AND. wfc_dir == 'undefined') CALL errore(&
           &'lr_read_wfc', TRIM( prefix )//'.wfc'//' not found',1)
      IF (.NOT.exst .AND. wfc_dir /= 'undefined') THEN
         WRITE( stdout, '(/5x,"Attempting to read from outdir &
              &instead of wfcdir")' )
         CLOSE( UNIT = iunwfc)
         tmp_dir = tmp_dir_saved
         CALL diropn ( iunwfc, 'wfc', 2*nwordwfc, exst)
         IF (.NOT.exst) CALL errore('lr_read_wfc', &
              &TRIM( prefix )//'.wfc'//' not found',1)
      ENDIF
      IF (gamma_only) THEN
         WRITE( stdout, '(/5x,"Gamma point algorithm")' )
      ELSE
         WRITE( stdout, '(/5x,"Generalised algorithm !warning")' )
      ENDIF
      DO ik=1,nks
         IF (.NOT. real_space_debug > 0 ) THEN 
            !else done in init_realspace realus
            CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw,&
                 & igk, g2kin ) 
            npw_k(ik) = npw
            igk_k(:,ik) = igk(:)
         ENDIF
         !   Read in the ground state wavefunctions
         !   This is a parallel read, done in wfc_dir
         CALL davcio(evc_all(:,:,ik),2*nwordwfc,iunwfc,ik,-1)
         !
      ENDDO
      CLOSE( UNIT = iunwfc)

      ! End of file reading
      tmp_dir = tmp_dir_saved

      ! vkb * evc_all and initialization of sevc_all
      IF ( okvan ) THEN
         IF ( gamma_only ) THEN
            ! Following line is to be removed when real space 
            ! implementation is complete
            CALL init_us_2(npw,igk_k(:,1),xk(1,1),vkb)
            IF (real_space_debug>0) THEN
               DO ibnd=1,nbnd,2
                  CALL fft_orbital_gamma(evc_all(:,:,1),ibnd,nbnd)
                  CALL calbec_rs_gamma(ibnd,nbnd,becp1_all)
                  becp%r(:,ibnd)=becp1_all(:,ibnd)
                  IF (ibnd + 1 <= nbnd) becp%r(:,ibnd+1)=becp1_all(:,ibnd+1)
                  CALL s_psir_gamma(ibnd,nbnd)
                  CALL bfft_orbital_gamma(sevc_all(:,:,1),ibnd,nbnd)
               ENDDO
            ELSE
               CALL calbec(npw_k(1),vkb,evc_all(:,:,1),becp1_all)
               becp%r=becp1_all
               CALL s_psi(npwx, npw_k(1), nbnd, evc_all(:,:,1), &
                    & sevc_all(:,:,1))
            ENDIF
         ELSE
            ! K point generalized stuff starts here
            DO ik=1,nks
               CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
               CALL calbec(npw_k(ik),vkb,evc_all(:,:,ik),&
                    & becp1_c_all(:,:,ik),nbnd)
               becp%k=becp1_c_all(:,:,ik)
               CALL s_psi (npwx, npw_k(ik), nbnd, evc_all(:,:,ik), &
                    & sevc_all(:,:,ik))
            ENDDO
         ENDIF
      ELSE
         sevc_all=evc_all
      ENDIF

      ! Inverse fourier transform of evc_all
      !
      !
      revc_all=(0.0d0,0.0d0)
      nbnd=nbnd_occ(1) ! Very important, otherwise there will be bugs 
      CALL deallocate_bec_type(becp)
      CALL allocate_bec_type ( nkb, nbnd, becp )
      !
      IF ( gamma_only ) THEN
         DO ibnd=1,nbnd,2
            IF (ibnd<nbnd) THEN
               DO ig=1,npw_k(1)
                  revc_all(nls(igk_k(ig,1)),ibnd,1)=evc_all(ig,ibnd,1)&
                       &+(0.0d0,1.0d0)*evc_all(ig,ibnd+1,1)
                  revc_all(nlsm(igk_k(ig,1)),ibnd,1)=&
                       &CONJG(evc_all(ig,ibnd,1)&
                       &-(0.0d0,1.0d0)*evc_all(ig,ibnd+1,1))
               ENDDO
            ELSE
               DO ig=1,npw_k(1)
                  revc_all(nls(igk_k(ig,1)),ibnd,1)=evc_all(ig,ibnd,1)
                  revc_all(nlsm(igk_k(ig,1)),ibnd,1)=&
                       &CONJG(evc_all(ig,ibnd,1))
               ENDDO
            ENDIF
            CALL invfft ('Wave', revc_all(:,ibnd,1), dffts)
         ENDDO
      ELSE
         DO ik=1,nks
            DO ibnd=1,nbnd
               DO ig=1,npw_k(ik)
                  revc_all(nls(igk_k(ig,ik)),ibnd,ik)=&
                       &evc_all(ig,ibnd,ik)
               ENDDO
               CALL invfft ('Wave', revc_all(:,ibnd,ik), dffts)
            ENDDO
         ENDDO
      ENDIF

      ! Now everything goes into right place
      evc0=(0.0d0,0.0d0)
      evc0(:,:,:)=evc_all(:,1:nbnd,:)
      !
      sevc0=(0.0d0,0.0d0)
      sevc0(:,:,:)=sevc_all(:,1:nbnd,:)
      !
      revc0=(0.0d0,0.0d0)
      revc0(:,:,:)=revc_all(:,1:nbnd,:)
      !
      IF (nkb>0) THEN
         IF (gamma_only) THEN
            becp1(:,:)=becp1_all(:,1:nbnd)
            becp%r=0.0d0
            becp%r=becp1
         ELSE
            becp1_c(:,:,:)=becp1_c_all(:,1:nbnd,:)
            becp%k=(0.0d0,0.0d0)
            becp%k=becp1_c(:,:,1)
         ENDIF
      ENDIF
      !
      ! Finally retain the conduction bands if needed for projection
      evc0_virt(:,:,:)=evc_all(:,nbnd+1:nbnd_total,:)
      sevc0_virt(:,:,:)=sevc_all(:,nbnd+1:nbnd_total,:)
      IF (nkb>0) THEN
        IF (gamma_only) THEN
          becp1_virt(:,:)=becp1_all(:,nbnd+1:nbnd_total)
        ELSE
          becp1_c_virt(:,:,:)=becp1_c_all(:,nbnd+1:nbnd_total,:)
        ENDIF
      ENDIF

      IF (nkb>0) THEN
         IF (gamma_only) THEN
            DEALLOCATE(becp1_all)
         ELSE
            DEALLOCATE(becp1_c_all)
         ENDIF
      ENDIF
      DEALLOCATE(evc_all)
      DEALLOCATE(sevc_all)
      DEALLOCATE(revc_all)
      ! OBM - Last minute check for real space implementation,
      IF ( real_space_debug > 0 .and. .not. gamma_only ) &
           CALL errore( ' iosys ', ' Linear response calculation ' // &
           & 'real space algorithms with k-points not implemented', 1 )
      !
    END SUBROUTINE virt_read
!-----------------------------------------------------------------------
END SUBROUTINE lr_read_wf
