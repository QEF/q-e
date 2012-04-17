!-----------------------------------------------------------------------
SUBROUTINE lr_read_wf()
  !---------------------------------------------------------------------
  ! ... reads in and stores the ground state wavefunctions
  ! ... for use in Lanczos linear response calculation
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu (2009)
  !
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, xk
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, diropn, tmp_dir, wfc_dir
  USE lr_variables,         ONLY : evc0, sevc0 ,revc0, evc0_virt, sevc0_virt, nbnd_total, &
                                   becp1_virt,becp1_c_virt
  USE realus,               ONLY : igk_k,npw_k
  USE lr_variables,         ONLY : becp1, becp1_c,test_case_no,size_evc,project
  USE wvfct,                ONLY : npw, igk, nbnd, g2kin, npwx, ecutwfc
  USE control_flags,        ONLY : gamma_only
  !use wavefunctions_module,only : evc
  USE gvecs,              ONLY : nls, nlsm
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE becmod,               ONLY : bec_type, becp, calbec
  !use lr_variables,         only : real_space
  !use real_beta,            only : ccalbecr_gamma,s_psir,fft_orbital_gamma,bfft_orbital_gamma
  USE realus,               ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma, real_space_debug
  USE lr_variables,   ONLY : lr_verbosity
  USE buffers,              ONLY : get_buffer
  USE kinds,                 ONLY : dp



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
  !
  IF (nbnd_total>nbnd .or. project) THEN
   CALL virt_read()
  ELSE
   CALL normal_read()
  ENDIF
  !
  !print *, "evc0",lr_dot(evc0(:,:,1),evc0(:,:,1))
  !print *, "sevc0",lr_dot(sevc0(:,:,1),sevc0(:,:,1))
  !print *, "<evc0|sevc0>",lr_dot(evc0(:,:,1),sevc0(:,:,1))
  !print *, "<revc0>",lr_dot(revc0(:,:,1),revc0(:,:,1))
  !print *, "becp1",lr_dot(becp1(:,:),becp1(:,:))
  RETURN
!!!!
  CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE normal_read()
!
!The usual way of reading wavefunctions
!
USE lr_variables, ONLY: check_all_bands_gamma, check_density_gamma,&
                        check_vector_gamma
  IMPLICIT NONE

  nwordwfc = 2 * nbnd * npwx
  size_evc=npwx*nbnd*nks
  !
  !   Read in the ground state wavefunctions
  !   This is a parallel read, done in wfc_dir
  tmp_dir_saved = tmp_dir
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
  CALL diropn ( iunwfc, 'wfc', nwordwfc, exst)
  !
  IF (.not.exst .and. wfc_dir == 'undefined') CALL errore('lr_read_wfc', trim( prefix )//'.wfc'//' not found',1)
  !
  IF (.not.exst .and. wfc_dir /= 'undefined') THEN
    WRITE( stdout, '(/5x,"Attempting to read wfc from outdir instead of wfcdir")' )
    CLOSE( UNIT = iunwfc)
    tmp_dir = tmp_dir_saved
    CALL diropn ( iunwfc, 'wfc', nwordwfc, exst)
    IF (.not.exst) CALL errore('lr_read_wfc', trim( prefix )//'.wfc'//' not found',1)
  ENDIF
  IF (gamma_only) THEN
   WRITE( stdout, '(/5x,"Gamma point algorithm")' )
  ELSE
   WRITE( stdout, '(/5x,"Generalised algorithm !warning")' )
  ENDIF

  DO ik=1,nks
     !
     IF (.not. real_space_debug > 0 ) THEN !else done in init_realspace realus
       CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
       !
       npw_k(ik) = npw
       !
       igk_k(:,ik) = igk(:)
     ENDIF
     !
     CALL davcio(evc0(:,:,ik),nwordwfc,iunwfc,ik,-1)
     !
  ENDDO
  !
  !
  CLOSE( UNIT = iunwfc)
  ! End of file reading
  tmp_dir = tmp_dir_saved
  !print * , "evc0 ",evc0(1:3,1,1)
  !
  !
  ! vkb * evc0 and initialization of sevc0
  !
  !
  IF ( okvan ) THEN
     !
     IF ( gamma_only ) THEN
        !
        ! Following line is to be removed when real space implementation is complete
        CALL init_us_2(npw,igk_k(:,1),xk(1,1),vkb)
        !
        IF (real_space_debug>0) THEN
        !
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
          !rbecp=becp1
          !print *,rbecp
          !
          IF (test_case_no == 1) THEN
           WRITE(stdout,'(/5x,"Test Case 1, dumping Real space calculated rbecp and sevc0",1x)')
           filename=trim(prefix) // "-rbecp-rs.dump"
           OPEN(UNIT=47,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')
           WRITE(unit=47,FMT='("#RBECP SIZE :",i6," number of beta fs",i6," bands",i6)') size(becp%r)&
                                                            ,size(becp%r,1),size(becp%r,2)

           DO itmp2=1, size(becp%r,2)
             WRITE(unit=47,FMT='("#Band no",i3)') itmp2
             DO itmp1=1, size(becp%r,1)
              WRITE(unit=47,FMT=*) becp%r(itmp1,itmp2)
             ENDDO
            ENDDO
           CLOSE(47)
           filename=trim(prefix) // "-sevc0-rs.dump"
           OPEN(UNIT=48,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')
           WRITE(unit=48,FMT='("#SEVC0 SIZE :",i6," NPW ",i6," BANDS ",i6," DIM3",i6)') size(sevc0), &
                                        size(sevc0,1), size(sevc0,2), size(sevc0,3)
           DO itmp2=1, size(sevc0,2)
             WRITE(unit=48,FMT='("#Band no",i3)') itmp2
             DO itmp1=1, size(sevc0,1)
               WRITE(unit=48,FMT='(i6,2x,e21.15, 2x, e21.15,2x)') itmp1, dble(sevc0(itmp1,itmp2,1)), aimag(sevc0(itmp1,itmp2,1))
             ENDDO
            ENDDO

           CLOSE(48)
          ENDIF
          !print *, becp1-rbecp
          !
          ! makedo part until spsi is in place - obsolote
          ! call s_psi(npwx, npw_k(1), nbnd, evc0(:,:,1), sevc0(:,:,1))
        ELSE
           !
           !call pw_gemm('Y',nkb,nbnd,npw_k(1),vkb,npwx,evc0,npwx,becp1,nkb)
           CALL calbec(npw_k(1),vkb,evc0(:,:,1),becp1)
           !
           becp%r=becp1
           !
           CALL s_psi(npwx, npw_k(1), nbnd, evc0(:,:,1), sevc0(:,:,1))
           ! Test case
           IF (test_case_no == 1) THEN
            WRITE(stdout,'(/5x,"Test Case 1, dumping Fourier space calculated rbecp and sevc0",1x)')
            filename=trim(prefix) // "-rbecp.dump"
            OPEN(UNIT=47,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')
            WRITE(unit=47,FMT='("#RBECP SIZE :",i6," number of beta fs",i6," bands",i6)') size(becp%r)&
                                                            ,size(becp%r,1),size(becp%r,2)
            DO itmp2=1, size(becp%r,2)
             WRITE(unit=47,FMT='("#Band no",i3)') itmp2
             DO itmp1=1, size(becp%r,1)
              WRITE(unit=47,FMT=*) becp%r(itmp1,itmp2)
             ENDDO
            ENDDO
            CLOSE(47)
            filename=trim(prefix) // "-sevc0.dump"
            OPEN(UNIT=48,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')
            WRITE(unit=48,FMT='("#SEVC0 SIZE :",i6," NPW ",i6," BANDS ",i6," DIM3",i6)') size(sevc0), &
                                        size(sevc0,1), size(sevc0,2), size(sevc0,3)
            DO itmp2=1, size(sevc0,2)
              WRITE(unit=48,FMT='("#Band no",i3)') itmp2
              DO itmp1=1, size(sevc0,1)
                WRITE(unit=48,FMT='(i6,2x,e21.15, 2x, e21.15,2x)') itmp1, dble(sevc0(itmp1,itmp2,1)), aimag(sevc0(itmp1,itmp2,1))
              ENDDO
             ENDDO
             CLOSE(48)
           ENDIF

           !
        ENDIF
     ELSE
        !
        ! K point generalized stuff starts here
        DO ik=1,nks
           !
           CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
           !
           !call ccalbec(nkb,npwx,npw_k(ik),nbnd,becp1_c(:,:,ik),vkb,evc0(:,:,ik))
           CALL calbec(npw_k(ik),vkb,evc0(:,:,ik),becp1_c(:,:,ik))
           !
           becp%k=becp1_c(:,:,ik)
           !
           CALL s_psi (npwx, npw_k(ik), nbnd, evc0(:,:,ik), sevc0(:,:,ik))
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
  revc0=(0.0d0,0.0d0)
  !
  IF ( gamma_only ) THEN
     !
     DO ibnd=1,nbnd,2
        !
        IF (ibnd<nbnd) THEN
           !
           DO ig=1,npw_k(1)
              !
              revc0(nls(igk_k(ig,1)),ibnd,1)=evc0(ig,ibnd,1)+&
                     (0.0d0,1.0d0)*evc0(ig,ibnd+1,1)
              !
              revc0(nlsm(igk_k(ig,1)),ibnd,1)=conjg(evc0(ig,ibnd,1)-&
                     (0.0d0,1.0d0)*evc0(ig,ibnd+1,1))
              !
           ENDDO
           !
        ELSE
           !
           DO ig=1,npw_k(1)
              !
              revc0(nls(igk_k(ig,1)),ibnd,1)=evc0(ig,ibnd,1)
              !
              revc0(nlsm(igk_k(ig,1)),ibnd,1)=conjg(evc0(ig,ibnd,1))
              !
           ENDDO
           !
        ENDIF
        !
        CALL invfft ('Wave', revc0(:,ibnd,1), dffts)
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
              revc0(nls(igk_k(ig,ik)),ibnd,ik)=evc0(ig,ibnd,ik)
              !
           ENDDO
           !
           CALL invfft ('Wave', revc0(:,ibnd,ik), dffts)
           !
        ENDDO
        !
     ENDDO
     !
  ENDIF

  ! OBM - Last minute check for real space implementation,
  IF ( real_space_debug > 0 .and. .not. gamma_only ) &
       CALL errore( ' iosys ', ' Linear response calculation ' // &
       & 'real space algorithms with k-points not implemented', 1 )
  !
END SUBROUTINE normal_read
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE virt_read()
!
!The modifications to read also the virtual orbitals
!
USE control_ph,            ONLY : nbnd_occ
USE fft_base,     ONLY : dfftp
USE lr_variables, ONLY: check_all_bands_gamma, check_density_gamma,&
                        check_vector_gamma
  IMPLICIT NONE
  COMPLEX(kind=dp), ALLOCATABLE :: evc_all(:,:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: sevc_all(:,:,:)
  real(kind=dp), ALLOCATABLE :: becp1_all(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: becp1_c_all(:,:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: revc_all(:,:,:)

  !First pretend everything is normal
  nbnd=nbnd_total
  !
  ALLOCATE(revc_all(dfftp%nnr,nbnd,nks))
  ALLOCATE(evc_all(npwx,nbnd,nks))
  ALLOCATE(sevc_all(npwx,nbnd,nks))
  IF (nkb > 0) THEN
    IF(gamma_only) THEN
       ALLOCATE(becp1_all(nkb,nbnd))
       becp1_all(:,:)=0.0d0
    ELSE
       ALLOCATE(becp1_c_all(nkb,nbnd,nks))
       becp1_c_all(:,:,:)=(0.0d0,0.0d0)
    ENDIF
  ENDIF


  nwordwfc = 2 * nbnd * npwx
  size_evc=npwx*nbnd_occ(1)*nks
  !
  !   Read in the ground state wavefunctions
  !   This is a parallel read, done in wfc_dir
  tmp_dir_saved = tmp_dir
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
  CALL diropn ( iunwfc, 'wfc', nwordwfc, exst)
  !
  IF (.not.exst .and. wfc_dir == 'undefined') CALL errore('lr_read_wfc', trim( prefix )//'.wfc'//' not found',1)
  !
  IF (.not.exst .and. wfc_dir /= 'undefined') THEN
    WRITE( stdout, '(/5x,"Attempting to read from outdir instead of wfcdir")' )
    CLOSE( UNIT = iunwfc)
    tmp_dir = tmp_dir_saved
    CALL diropn ( iunwfc, 'wfc', nwordwfc, exst)
    IF (.not.exst) CALL errore('lr_read_wfc', trim( prefix )//'.wfc'//' not found',1)
  ENDIF
  !
  IF (gamma_only) THEN
   WRITE( stdout, '(/5x,"Gamma point algorithm")' )
  ELSE
   WRITE( stdout, '(/5x,"Generalised algorithm !warning")' )
  ENDIF

  DO ik=1,nks
     !
     IF (.not. real_space_debug > 0 ) THEN !else done in init_realspace realus
       CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
       !
       npw_k(ik) = npw
       !
       igk_k(:,ik) = igk(:)
     ENDIF
     !
     !   Read in the ground state wavefunctions
     !   This is a parallel read, done in wfc_dir
     CALL davcio(evc_all(:,:,ik),nwordwfc,iunwfc,ik,-1)
     !
  ENDDO
  !
  !
  CLOSE( UNIT = iunwfc)
  ! End of file reading
  tmp_dir = tmp_dir_saved
  !print * , "evc_all ",evc_all(1:3,1,1)
  !
  !
  ! vkb * evc_all and initialization of sevc_all
  !
  !
  IF ( okvan ) THEN
     !
     IF ( gamma_only ) THEN
        !
        ! Following line is to be removed when real space implementation is complete
        CALL init_us_2(npw,igk_k(:,1),xk(1,1),vkb)
        !
        IF (real_space_debug>0) THEN
        !
         !
         !
          DO ibnd=1,nbnd,2
             CALL fft_orbital_gamma(evc_all(:,:,1),ibnd,nbnd)
             CALL calbec_rs_gamma(ibnd,nbnd,becp1_all)
             becp%r(:,ibnd)=becp1_all(:,ibnd)
             IF (ibnd + 1 <= nbnd) becp%r(:,ibnd+1)=becp1_all(:,ibnd+1)
             CALL s_psir_gamma(ibnd,nbnd)
             CALL bfft_orbital_gamma(sevc_all(:,:,1),ibnd,nbnd)
          ENDDO
        ELSE
           !
           CALL calbec(npw_k(1),vkb,evc_all(:,:,1),becp1_all)
           !
           becp%r=becp1_all
           !
           CALL s_psi(npwx, npw_k(1), nbnd, evc_all(:,:,1), sevc_all(:,:,1))
           !
        ENDIF
     ELSE
        !
        ! K point generalized stuff starts here
        DO ik=1,nks
           !
           CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
           !
           CALL calbec(npw_k(ik),vkb,evc_all(:,:,ik),becp1_c_all(:,:,ik),nbnd)
           !
           becp%k=becp1_c_all(:,:,ik)
           !
           CALL s_psi (npwx, npw_k(ik), nbnd, evc_all(:,:,ik), sevc_all(:,:,ik))
           !
        ENDDO
        !
     ENDIF
     !
  ELSE
     !
     sevc_all=evc_all
     !
  ENDIF
  !
  !
  ! Inverse fourier transform of evc_all
  !
  !
  revc_all=(0.0d0,0.0d0)
  !
  IF ( gamma_only ) THEN
     !
     DO ibnd=1,nbnd,2
        !
        IF (ibnd<nbnd) THEN
           !
           DO ig=1,npw_k(1)
              !
              revc_all(nls(igk_k(ig,1)),ibnd,1)=evc_all(ig,ibnd,1)+&
                     (0.0d0,1.0d0)*evc_all(ig,ibnd+1,1)
              !
              revc_all(nlsm(igk_k(ig,1)),ibnd,1)=conjg(evc_all(ig,ibnd,1)-&
                     (0.0d0,1.0d0)*evc_all(ig,ibnd+1,1))
              !
           ENDDO
           !
        ELSE
           !
           DO ig=1,npw_k(1)
              !
              revc_all(nls(igk_k(ig,1)),ibnd,1)=evc_all(ig,ibnd,1)
              !
              revc_all(nlsm(igk_k(ig,1)),ibnd,1)=conjg(evc_all(ig,ibnd,1))
              !
           ENDDO
           !
        ENDIF
        !
        CALL invfft ('Wave', revc_all(:,ibnd,1), dffts)
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
              revc_all(nls(igk_k(ig,ik)),ibnd,ik)=evc_all(ig,ibnd,ik)
              !
           ENDDO
           !
           CALL invfft ('Wave', revc_all(:,ibnd,ik), dffts)
           !
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  !now everything goes into right place
  !
  nbnd=nbnd_occ(1)
  !
  evc0=(0.0d0,0.0d0)
  evc0(:,:,:)=evc_all(:,1:nbnd,:)
  sevc0=(0.0d0,0.0d0)
  sevc0(:,:,:)=sevc_all(:,1:nbnd,:)
  revc0=(0.0d0,0.0d0)
  revc0(:,:,:)=revc_all(:,1:nbnd,:)
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
  IF (project) THEN
   evc0_virt(:,:,:)=evc_all(:,nbnd+1:nbnd_total,:)
   !sevc0_virt(:,:,:)=sevc_all(:,nbnd+1:nbnd_total,:)
   IF (nkb>0) THEN
   IF (gamma_only) THEN
    becp1_virt(:,:)=becp1_all(:,nbnd+1:nbnd_total)
   ELSE
    becp1_c_virt(:,:,:)=becp1_c_all(:,nbnd+1:nbnd_total,:)
   ENDIF
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
