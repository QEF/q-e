!-----------------------------------------------------------------------
subroutine lr_read_wf()
  !---------------------------------------------------------------------
  ! ... reads in and stores the ground state wavefunctions 
  ! ... for use in Lanczos linear response calculation
  !---------------------------------------------------------------------
  !
  ! OBM :
  ! 050608 Modified for calbec interface in v4.0 (w evcx->evcx(:,:,ik or 1)
  !        gamma_only correction
#include "f_defs.h"
  !
  use io_global,            only : stdout
  use klist,                only : nks, xk 
  use cell_base,            only : tpiba2
  use gvect,                only : ngm, g, ecutwfc
  use io_files,             only : nwordwfc, iunwfc, prefix
  use lr_variables,         only : evc0, sevc0 ,revc0
  use realus,               only : igk_k,npw_k
  use lr_variables,         only : becp1, becp1_c,test_case_no,size_evc
  use wvfct,                only : npw, igk, nbnd, g2kin, npwx
  use control_flags,        only : gamma_only
  !use wavefunctions_module, only : evc
  use gsmooth,              only : nr1s, nr2s, nr3s,nrx1s, nrx2s, nrx3s, nls, nlsm
  use uspp,                 only : vkb, nkb, okvan
  use becmod,               only : bec_type, becp, calbec
  !use lr_variables,         only : real_space
  !use real_beta,            only : ccalbecr_gamma,s_psir,fft_orbital_gamma,bfft_orbital_gamma
  USE realus,               ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma, real_space_debug
  USE lr_variables,   ONLY : lr_verbosity 
  USE buffers,              ONLY : get_buffer 
  USE kinds,                 ONLY : dp

  

  !
  implicit none
  !
  !
  
  !
  ! local variables
  integer :: ik, ibnd, ig, itmp1,itmp2,itmp3
  logical :: exst
  character(len=256) :: filename
  real(kind=dp) :: obm_debug
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_read_wf>")')
  endif
  !
  nwordwfc = 2 * nbnd * npwx
  size_evc=npwx*nbnd*nks
  !
  call diropn ( iunwfc, 'wfc', nwordwfc, exst)
  !
  if (.not.exst) call errore('lr_read_wfc', TRIM( prefix )//'.wfc'//' not found',1)
  !
  if (gamma_only) then 
   WRITE( stdout, '(/5x,"Gamma point algorithm")' )
  else
   call errore('lr_read_wfc', 'k-point algorithm is not tested yet',1)
   WRITE( stdout, '(/5x,"Generalised algorithm !warning")' ) 
  endif

  do ik=1,nks
     !
     if (.not. real_space_debug > 0 ) then !else done in init_realspace realus
       CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
       !
       npw_k(ik) = npw
       !
       igk_k(:,ik) = igk(:)
     endif
     !
     !   Read in the ground state wavefunctions
     !
     call davcio(evc0(:,:,ik),nwordwfc,iunwfc,ik,-1)
     !
  enddo
  !
  !
  CLOSE( UNIT = iunwfc)
  !print * , "evc0 ",evc0(1:3,1,1) 
  !
  !
  ! vkb * evc0 and initialization of sevc0
  !
  !
  if ( okvan ) then
     !
     if ( gamma_only ) then
        !
        ! Following line is to be removed when real space implementation is complete
        call init_us_2(npw,igk_k(:,1),xk(1,1),vkb)
        !
        if (real_space_debug>0) then
        ! 
         !
         !
          do ibnd=1,nbnd,2
             call fft_orbital_gamma(evc0(:,:,1),ibnd,nbnd)
             call calbec_rs_gamma(ibnd,nbnd,becp1)
             becp%r(:,ibnd)=becp1(:,ibnd)
             if (ibnd + 1 .le. nbnd) becp%r(:,ibnd+1)=becp1(:,ibnd+1)
             call s_psir_gamma(ibnd,nbnd)
             call bfft_orbital_gamma(sevc0(:,:,1),ibnd,nbnd)
          enddo
          !rbecp=becp1
          !print *,rbecp
          !
          if (test_case_no .eq. 1) then
           write(stdout,'(/5x,"Test Case 1, dumping Real space calculated rbecp and sevc0",1x)')
           filename=trim(prefix) // "-rbecp-rs.dump"
           OPEN(UNIT=47,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')
           write(unit=47,FMT='("#RBECP SIZE :",i6," number of beta fs",i6," bands",i6)') size(becp%r)&
                                                            ,size(becp%r,1),size(becp%r,2) 
          
           do itmp2=1, SIZE(becp%r,2)
             write(unit=47,FMT='("#Band no",i3)') itmp2
             do itmp1=1, SIZE(becp%r,1)
              write(unit=47,FMT=*) becp%r(itmp1,itmp2)
             enddo
            enddo
           close(47)
           filename=trim(prefix) // "-sevc0-rs.dump"
           OPEN(UNIT=48,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')
           write(unit=48,FMT='("#SEVC0 SIZE :",i6," NPW ",i6," BANDS ",i6," DIM3",i6)') size(sevc0), &
                                        size(sevc0,1), size(sevc0,2), size(sevc0,3)
           do itmp2=1, SIZE(sevc0,2)
             write(unit=48,FMT='("#Band no",i3)') itmp2
             do itmp1=1, SIZE(sevc0,1)
               write(unit=48,FMT='(i6,2x,e21.15, 2x, e21.15,2x)') itmp1, DBLE(sevc0(itmp1,itmp2,1)), AIMAG(sevc0(itmp1,itmp2,1))
             enddo
            enddo
           
           close(48)
          endif
          !print *, becp1-rbecp
          !
          ! makedo part until spsi is in place - obsolote
          ! call s_psi(npwx, npw_k(1), nbnd, evc0(:,:,1), sevc0(:,:,1))
        else
           !
           !call pw_gemm('Y',nkb,nbnd,npw_k(1),vkb,npwx,evc0,npwx,becp1,nkb) 
           call calbec(npw_k(1),vkb,evc0(:,:,1),becp1)
           !
           becp%r=becp1
           !
           call s_psi(npwx, npw_k(1), nbnd, evc0(:,:,1), sevc0(:,:,1))
           ! Test case
           if (test_case_no .eq. 1) then 
            write(stdout,'(/5x,"Test Case 1, dumping Fourier space calculated rbecp and sevc0",1x)')
            filename=trim(prefix) // "-rbecp.dump"
            OPEN(UNIT=47,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')
            write(unit=47,FMT='("#RBECP SIZE :",i6," number of beta fs",i6," bands",i6)') size(becp%r)&
                                                            ,size(becp%r,1),size(becp%r,2)  
            do itmp2=1, SIZE(becp%r,2)
             write(unit=47,FMT='("#Band no",i3)') itmp2
             do itmp1=1, SIZE(becp%r,1)
              write(unit=47,FMT=*) becp%r(itmp1,itmp2)
             enddo
            enddo
            close(47)
            filename=trim(prefix) // "-sevc0.dump"
            OPEN(UNIT=48,FILE=filename,STATUS='NEW',ACCESS = 'SEQUENTIAL')  
            write(unit=48,FMT='("#SEVC0 SIZE :",i6," NPW ",i6," BANDS ",i6," DIM3",i6)') size(sevc0), &
                                        size(sevc0,1), size(sevc0,2), size(sevc0,3)
            do itmp2=1, SIZE(sevc0,2)
              write(unit=48,FMT='("#Band no",i3)') itmp2
              do itmp1=1, SIZE(sevc0,1)
                write(unit=48,FMT='(i6,2x,e21.15, 2x, e21.15,2x)') itmp1, DBLE(sevc0(itmp1,itmp2,1)), AIMAG(sevc0(itmp1,itmp2,1))
              enddo
             enddo
             close(48)
           endif

           !
        endif
     else
        !
        ! K point generalized stuff starts here
        do ik=1,nks
           !
           call init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
           !
           !call ccalbec(nkb,npwx,npw_k(ik),nbnd,becp1_c(:,:,ik),vkb,evc0(:,:,ik))
           call calbec(npw_k(ik),vkb,evc0(:,:,ik),becp1_c(:,:,ik))
           !
           becp%k=becp1_c(:,:,ik)
           !
           call s_psi (npwx, npw_k(ik), nbnd, evc0(:,:,ik), sevc0(:,:,ik))
           !
        end do
        !
     end if
     !
  else
     !
     sevc0=evc0
     !
  end if
  !
  !
  ! Inverse fourier transform of evc0
  !
  !
  revc0=(0.0d0,0.0d0)
  !
  if ( gamma_only ) then
     !
     do ibnd=1,nbnd,2
        !
        if (ibnd<nbnd) then
           !
           do ig=1,npw_k(1)
              !
              revc0(nls(igk_k(ig,1)),ibnd,1)=evc0(ig,ibnd,1)+&
                     (0.0d0,1.0d0)*evc0(ig,ibnd+1,1)
              !
              revc0(nlsm(igk_k(ig,1)),ibnd,1)=conjg(evc0(ig,ibnd,1)-&
                     (0.0d0,1.0d0)*evc0(ig,ibnd+1,1))
              !
           end do
           !
        else
           !
           do ig=1,npw_k(1)
              !
              revc0(nls(igk_k(ig,1)),ibnd,1)=evc0(ig,ibnd,1)
              !
              revc0(nlsm(igk_k(ig,1)),ibnd,1)=conjg(evc0(ig,ibnd,1))
              !
           enddo
           !
        endif
        !
        call cft3s(revc0(:,ibnd,1),nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,2)
        !
     end do
     !
  else
     !
     do ik=1,nks
        !
        do ibnd=1,nbnd
           !
           do ig=1,npw_k(ik)
              !
              revc0(nls(igk_k(ig,ik)),ibnd,ik)=evc0(ig,ibnd,ik)
              !
           enddo
           !
           call cft3s(revc0(:,ibnd,ik),nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,2)
           !
        end do
        !
     end do
     !
  end if
  !
  !print * , "evc0 ",evc0(1:3,1,1) 
  !
  !
  !OBM!!! debug---
  !CALL lr_normalise( evc0(:,:,1), obm_debug)
  !print *, "norm of evc0 ",obm_debug
  !OBM!!! debug---
  ! OBM - Last minute check for real space implementation,
  IF ( real_space_debug > 0 .and. .NOT. gamma_only ) &
       CALL errore( ' iosys ', ' Linear response calculation ' // &
       & 'real space algorithms with k-points not implemented', 1 )
  !

  return
end subroutine lr_read_wf
!-----------------------------------------------------------------------
