! Copyrigh(C) 2005-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE exx_std
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only, tqr
  USE io_global,            ONLY : stdout
  USE exx_base,             ONLY : dfftt, exxbuff, locbuff, exxmat, locmat, evc0, &
                                   xkq_collect, index_xk, npwt, exxalfa, eps_occ, &
                                   ibnd_start, ibnd_end, ibnd_buff_start, ibnd_buff_end
  USE klist,                ONLY : nks, igk_k
  USE noncollin_module,     ONLY : noncolin, npol
  !
 CONTAINS
#define _CX(A)  CMPLX(A,0._dp,kind=DP)
#define _CY(A)  CMPLX(0._dp,-A,kind=DP)
  !
  !------------------------------------------------------------------------
  SUBROUTINE exxinit_std( DoLoc )
    !------------------------------------------------------------------------
    !! This subroutine is run before the first H_psi() of each iteration. 
    !! It saves the wavefunctions for the right density matrix, in real space.
    !
    USE wavefunctions,        ONLY : evc, psic
    USE io_files,             ONLY : nwordwfc, iunwfc
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE klist,                ONLY : ngk, nkstot, xk, wk
    USE symm_base,            ONLY : nsym, s, sr
    USE mp_pools,             ONLY : npool, nproc_pool, me_pool, inter_pool_comm
    USE mp_exx,               ONLY : me_egrp, negrp, init_index_over_band,  &
                                     my_egrp_id, inter_egrp_comm,           &
                                     intra_egrp_comm, iexx_start, iexx_end, &
                                     all_start, all_end
    USE mp,                   ONLY : mp_sum, mp_bcast
    USE xc_lib,               ONLY : xclib_get_exx_fraction, start_exx,          &
                                     get_screening_parameter, get_gau_parameter, &
                                     exx_is_active
    USE scatter_mod,          ONLY : gather_grid, scatter_grid
    USE fft_interfaces,       ONLY : invfft
    USE uspp,                 ONLY : nkb, vkb, okvan
    USE paw_variables,        ONLY : okpaw
    USE mp_orthopools,        ONLY : intra_orthopool_comm
    USE exx_base,             ONLY : nkqs, index_sym,  &
                                     exx_set_symm, rir, working_pool, exxdiv, &
                                     erfc_scrlen, gau_scrlen, exx_divergence, &
                                     x_nbnd_occ, nbndproj, local_thr, d_spin
#if defined(__CUDA)
    USE device_memcpy_m,      ONLY : dev_memset
    USE device_fbuff_m,       ONLY : dev_buf
#endif
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: DoLoc
    !! TRUE:  Real Array locbuff(ir, nbnd, nkqs);  
    !! FALSE: Complex Array exxbuff(ir, nbnd/2, nkqs).
    !
    ! ... local variables
    !
    INTEGER :: ik, ibnd, i, j, k, ir, isym, ikq, ig, ierr
    INTEGER :: ibnd_loop_start
    INTEGER :: ipol, jpol
    REAL(DP), ALLOCATABLE :: occ(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: temppsic
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) temppsic
#endif
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    COMPLEX(DP),POINTER     :: psic_nc_d(:,:)
#if defined(__CUDA)
    attributes(DEVICE)      :: psic_nc_d
#endif
    COMPLEX(DP),ALLOCATABLE :: psic_exx(:)
    INTEGER :: nxxs, nrxxs
#if defined(__MPI)
    COMPLEX(DP),ALLOCATABLE  :: temppsic_all(:), psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    INTEGER :: npw, current_ik
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    INTEGER :: h_ibnd
    INTEGER :: evc_offset
    !
    IF ( Doloc ) THEN
        WRITE(stdout,'(/,5X,"Using localization algorithm with threshold: ",&
                & D10.2)') local_thr
        ! IF (.NOT.gamma_only) CALL errore('exxinit','SCDM with K-points NYI',1)
        IF (okvan .OR. okpaw) CALL errore( 'exxinit','SCDM with USPP/PAW not &
                                           &implemented', 1 )
    ENDIF 
    !
    ! Note that nxxs is not the same as nrxxs in parallel case
    nxxs = dfftt%nr1x * dfftt%nr2x * dfftt%nr3x
    nrxxs = dfftt%nnr
    !
#if defined(__MPI)
    IF (noncolin) THEN
       ALLOCATE( psic_all_nc(nxxs,npol), temppsic_all_nc(nxxs,npol) )
    ELSEIF ( .NOT. gamma_only ) THEN
       ALLOCATE( psic_all(nxxs), temppsic_all(nxxs) )
    ENDIF
#endif
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs, npol), psic_nc(nrxxs, npol) )
    ELSEIF ( .NOT. gamma_only ) THEN
       ALLOCATE( temppsic(nrxxs) )
    ENDIF
    !
    CALL divide( inter_egrp_comm, x_nbnd_occ, ibnd_start, ibnd_end )
    !
    IF ( gamma_only ) THEN
        ibnd_buff_start = ibnd_start/2
        IF(mod(ibnd_start,2)==1) ibnd_buff_start = ibnd_buff_start +1
        !
        ibnd_buff_end = ibnd_end/2
        IF(mod(ibnd_end,2)==1) ibnd_buff_end = ibnd_buff_end +1
    ELSE
        ibnd_buff_start = ibnd_start
        ibnd_buff_end   = ibnd_end
    ENDIF
    !
    IF (DoLoc) THEN
      !
      IF (gamma_only) THEN
        IF (.NOT. ALLOCATED(locbuff)) ALLOCATE( locbuff(nrxxs*npol,nbnd,nks) )
        IF (.NOT. ALLOCATED(locmat))  ALLOCATE( locmat(nbnd,nbnd,nks) )
        locbuff = 0.0d0
        locmat = 0.0d0
      ELSE 
        IF (.NOT. ALLOCATED(exxbuff)) ALLOCATE( exxbuff(nrxxs*npol,nbnd,nkqs) )
        IF (.NOT. ALLOCATED(exxmat) ) ALLOCATE( exxmat(nbnd,nkqs,nbnd,nks) )
        exxbuff = (0.0d0, 0.0d0)
        exxmat = 0.0d0
      ENDIF
      !
      IF (.NOT. ALLOCATED(evc0)) then 
        ALLOCATE( evc0(npwx*npol,nbndproj,nks) )
        evc0 = (0.0d0,0.0d0)
      ENDIF
      !
    ELSE
      !
      IF (.NOT. ALLOCATED(exxbuff)) THEN
        IF (gamma_only) THEN
           ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nks))
        ELSE
           ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nkqs))
        END IF
      END IF
      !
    ENDIF
    !
    !assign buffer
    IF(DoLoc) THEN
      locbuff = 0.0_DP
    ELSE
      exxbuff = (0.0_DP,0.0_DP)
    ENDIF
    !
    ! ... This is parallelized over pools. Each pool computes only its k-points
    !
    KPOINTS_LOOP : &
    DO ik = 1, nks
       !
       IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
       !
       ! ik         = index of k-point in this pool
       ! current_ik = index of k-point over all pools
       !
       current_ik = global_kpoint_index( nkstot, ik )
       !
       IF_GAMMA_ONLY : &
       IF (gamma_only) THEN
          !
          h_ibnd = ibnd_start/2
          !
          IF(mod(ibnd_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start-1
          ELSE
             ibnd_loop_start=ibnd_start
          ENDIF
          ! 
          DO ibnd = ibnd_loop_start, ibnd_end, 2
             h_ibnd = h_ibnd + 1
             !
             psic(:) = ( 0._DP, 0._DP )
             !
             IF ( ibnd < ibnd_end ) THEN
                DO ig = 1, npwt
                   psic(dfftt%nl(ig))  = evc(ig,ibnd) + ( 0._dp, 1._dp ) * evc(ig,ibnd+1)
                   psic(dfftt%nlm(ig)) = conjg( evc(ig,ibnd) ) + ( 0._dp, 1._dp ) * conjg( evc(ig,ibnd+1) )
                ENDDO
             ELSE
                DO ig = 1, npwt
                   psic(dfftt%nl(ig))  = evc(ig,ibnd)
                   psic(dfftt%nlm(ig)) = CONJG( evc(ig,ibnd) ) 
                ENDDO
             ENDIF
             !
             CALL invfft( 'Wave', psic, dfftt )
             !
             IF (DoLoc) THEN
               locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+1,ik) = DBLE(  psic_exx(1:nrxxs) )
               IF (ibnd-ibnd_loop_start+evc_offset+2 <= nbnd) &
                  locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+2,ik) = AIMAG( psic_exx(1:nrxxs) )
             ELSE
               exxbuff(1:nrxxs,h_ibnd,ik)=psic(1:nrxxs) 
             ENDIF
             !
          ENDDO
          !
       ELSE IF_GAMMA_ONLY
          !
          npw = ngk (ik)
          IBND_LOOP_K : &
          DO ibnd = ibnd_start, ibnd_end 
             !
             IF (noncolin) THEN
                temppsic_nc(:,:) = ( 0._DP, 0._DP )
                temppsic_nc(dfftt%nl(igk_k(1:npw,ik)),1) = evc(1:npw,ibnd)
                CALL invfft( 'Wave', temppsic_nc(:,1), dfftt )
                temppsic_nc(dfftt%nl(igk_k(1:npw,ik)),2) = evc(npwx+1:npwx+npw,ibnd)
                CALL invfft( 'Wave', temppsic_nc(:,2), dfftt )
             ELSE
                temppsic(:) = ( 0._DP, 0._DP )
                temppsic(dfftt%nl(igk_k(1:npw,ik))) = evc(1:npw,ibnd)
                CALL invfft( 'Wave', temppsic, dfftt )
             ENDIF
             !
             DO ikq = 1, nkqs
                !
                IF (index_xk(ikq) /= current_ik) CYCLE
                isym = ABS(index_sym(ikq) )
                !
                IF (noncolin) THEN ! noncolinear
#if defined(__MPI)
                   DO ipol = 1, npol
                      CALL gather_grid( dfftt, temppsic_nc(:,ipol), temppsic_all_nc(:,ipol) )
                   ENDDO
                   !
                   IF ( me_egrp == 0 ) THEN
                      psic_all_nc(:,:) = (0.0_DP, 0.0_DP)
                      DO ipol = 1, npol
                         DO jpol = 1, npol
                            psic_all_nc(:,ipol) = psic_all_nc(:,ipol) + &
                                          CONJG(d_spin(jpol,ipol,isym)) * &
                                          temppsic_all_nc(rir(:,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDIF
                   !
                   DO ipol = 1, npol
                      CALL scatter_grid( dfftt, psic_all_nc(:,ipol), psic_nc(:,ipol) )
                   ENDDO
#else
                   DO ipol = 1, npol
                      psic_nc(:,ipol) = (0._DP,0._DP)
                      DO jpol = 1, npol
                         psic_nc(:,ipol) = psic_nc(:,ipol) + CONJG(d_spin(jpol,ipol,isym))* &
                                            temppsic_nc(rir(:,isym),jpol)
                      ENDDO
                   ENDDO
#endif
                   !
!civn: not sure about these two cases are really different in this new(old) implementation
                   IF (index_sym(ikq) > 0 ) THEN
                      ! sym. op. without time reversal: normal case
                      exxbuff(1:nrxxs,ibnd,ikq) = psic_nc(:,1)
                      exxbuff(nrxxs+1:nrxxs+nrxxs,ibnd,ikq) = psic_nc(:,2)
                   ELSE
                      ! sym. op. with time reversal: spin 1->2*, 2->-1*
                      exxbuff(1:nrxxs,ibnd,ikq) = CONJG(psic_nc(:,2))
                      exxbuff(nrxxs+1:nrxxs+nrxxs,ibnd,ikq) = -CONJG(psic_nc(:,1))
                   ENDIF
                ELSE ! noncolinear
#if defined(__MPI)
                   CALL gather_grid( dfftt, temppsic, temppsic_all )
                   IF ( me_egrp == 0 ) psic_all(:) = temppsic_all(rir(:,isym))
                   CALL scatter_grid( dfftt, psic_all, psic )
#else
                   psic(:) = temppsic(rir(:,isym))
#endif
                   IF (index_sym(ikq) < 0 ) psic(1:nrxxs) = CONJG(psic(1:nrxxs))
                   exxbuff(1:nrxxs,ibnd,ikq) = psic(1:nrxxs)
                   !
                ENDIF ! noncolinear
                !
             ENDDO
             !
          ENDDO&
          IBND_LOOP_K
          !
       ENDIF&
       IF_GAMMA_ONLY
    ENDDO&
    KPOINTS_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc, psic_nc )
#if defined(__MPI)
       DEALLOCATE( temppsic_all_nc, psic_all_nc )
#endif
    ELSE IF ( .NOT. gamma_only ) THEN
       DEALLOCATE( temppsic )
#if defined(__MPI)
       DEALLOCATE( temppsic_all, psic_all )
#endif
    ENDIF
    !
    ! Each wavefunction in exxbuff is computed by a single pool, collect among 
    ! pools in a smart way (i.e. without doing all-to-all sum and bcast)
    ! See also the initialization of working_pool in exx_mp_init
    ! Note that in Gamma-only LSDA can be parallelized over two pools, and there
    ! is no need to communicate anything: each pools deals with its own spin
    !
    IF ( .NOT. gamma_only ) THEN
       DO ikq = 1, nkqs
         CALL mp_bcast( exxbuff(:,:,ikq), working_pool(ikq), intra_orthopool_comm ) 
       ENDDO
    ENDIF
    !
  END SUBROUTINE exxinit_std
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_std_gamma(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... Gamma-specific version of vexx
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, current_k
    USE klist,          ONLY : xk, nks, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, nbgrp
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : npwt, gt, nqs, index_xkq, x_occupation, g2_convolution
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: RESULT(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL  :: global_kpoint_index
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    INTEGER :: ngmt
    !
    ngmt = dfftt%ngm
    ALLOCATE( fac(ngmt) )
    nrxxs= dfftt%nnr
    !
    ALLOCATE( RESULT(nrxxs), temppsic_dble(nrxxs), temppsic_aimag(nrxxs) )
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=0.0_DP
       psi =0.0_DP
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    ENDIF
    !
    ! Here the loops start
    !
    INTERNAL_LOOP_ON_Q : &
    DO iq=1,nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution(ngmt, gt, xkp, xkq, fac)
       IF ( okvan .and..not.tqr ) CALL qvan_init (ngmt, xkq, xkp)
       !
       LOOP_ON_PSI_BANDS : &
       DO im = 1,m !for each band of psi (the k cycle is outside band)
          IF(okvan) deexx(:) = 0.0_DP
          !
          RESULT = 0.0_DP
          !
          l_fft_doubleband = .false.
          l_fft_singleband = .false.
          !
          IF ( mod(im,2)==1 .and. (im+1)<=m ) l_fft_doubleband = .true.
          IF ( mod(im,2)==1 .and. im==m )     l_fft_singleband = .true.
          !
          IF( l_fft_doubleband ) THEN
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, npwt
                RESULT( dfftt%nl(ig) )  =       psi(ig, im) + (0._DP,1._DP) * psi(ig, im+1)
                RESULT( dfftt%nlm(ig) ) = conjg(psi(ig, im) - (0._DP,1._DP) * psi(ig, im+1))
             ENDDO
!$omp end parallel do
          ENDIF
          !
          IF( l_fft_singleband ) THEN
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, npwt
                RESULT( dfftt%nl(ig) )  =       psi(ig,im)
                RESULT( dfftt%nlm(ig) ) = conjg(psi(ig,im))
             ENDDO
!$omp end parallel do
          ENDIF
          !
          IF( l_fft_doubleband.or.l_fft_singleband) THEN
             CALL invfft ('Wave', RESULT, dfftt)
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                temppsic_dble(ir)  = dble ( RESULT(ir) )
                temppsic_aimag(ir) = aimag( RESULT(ir) )
             ENDDO
!$omp end parallel do
          ENDIF
          !
          RESULT = 0.0_DP
          !
          h_ibnd = ibnd_start/2
          IF(mod(ibnd_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start-1
          ELSE
             ibnd_loop_start=ibnd_start
          ENDIF
          !
          IBND_LOOP_GAM : &
          DO ibnd=ibnd_loop_start,ibnd_end, 2 !for each band of psi
             !
             h_ibnd = h_ibnd + 1
             IF( ibnd < ibnd_start ) THEN
                x1 = 0.0_DP
             ELSE
                x1 = x_occupation(ibnd,  ik)
             ENDIF
             IF( ibnd == ibnd_end) THEN
                x2 = 0.0_DP
             ELSE
                x2 = x_occupation(ibnd+1,  ik)
             ENDIF
             IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
             !
             ! calculate rho in real space. Gamma tricks are used.
             ! temppsic is real; tempphic contains one band in the real part,
             ! another one in the imaginary part; the same applies to rhoc
             !
             IF( mod(im,2) == 0 ) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = exxbuff(ir,h_ibnd,current_k) * temppsic_aimag(ir) / omega
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = exxbuff(ir,h_ibnd,current_k) * temppsic_dble(ir) / omega
                ENDDO
!$omp end parallel do
             ENDIF
             !
             ! bring rho to G-space
             !
             !   >>>> add augmentation in REAL SPACE here
             IF(okvan .and. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL addusxx_r(rhoc,_CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,im)))
                IF(ibnd<ibnd_end) &
                CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),_CX(becpsi%r(:,im)))
             ENDIF
             !
             CALL fwfft ('Rho', rhoc, dfftt)
             !   >>>> add augmentation in G SPACE here
             IF(okvan .and. .not. tqr) THEN
                ! contribution from one band added to real (in real space) part of rhoc
                IF(ibnd>=ibnd_start) &
                   CALL addusxx_g(dfftt, rhoc, xkq,  xkp, 'r', &
                   becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,im) )
                ! contribution from following band added to imaginary (in real space) part of rhoc
                IF(ibnd<ibnd_end) &
                   CALL addusxx_g(dfftt, rhoc, xkq,  xkp, 'i', &
                   becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,im) )
             ENDIF
             !   >>>> charge density done
             !
             vc = 0._DP
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, ngmt
                !
                vc(dfftt%nl(ig))  = fac(ig) * rhoc(dfftt%nl(ig))
                vc(dfftt%nlm(ig)) = fac(ig) * rhoc(dfftt%nlm(ig))
                !
             ENDDO
!$omp end parallel do
             !
             !   >>>>  compute <psi|H_fock G SPACE here
             IF(okvan .and. .not. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL newdxx_g(dfftt, vc, xkq, xkp, 'r', deexx, &
                              becphi_r=x1*becxx(ikq)%r(:,ibnd))
                IF(ibnd<ibnd_end) &
                CALL newdxx_g(dfftt, vc, xkq, xkp, 'i', deexx, &
                              becphi_r=x2*becxx(ikq)%r(:,ibnd+1))
             ENDIF
             !
             !brings back v in real space
             CALL invfft ('Rho', vc, dfftt)
             !
             !   >>>>  compute <psi|H_fock REAL SPACE here
             IF(okvan .and. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL newdxx_r(dfftt, vc, _CX(x1*becxx(ikq)%r(:,ibnd)), deexx)
                IF(ibnd<ibnd_end) &
                CALL newdxx_r(dfftt, vc, _CY(x2*becxx(ikq)%r(:,ibnd+1)), deexx)
             ENDIF
             !
             IF(okpaw) THEN
                IF(ibnd>=ibnd_start) &
                CALL PAW_newdxx(x1/nqs, _CX(becxx(ikq)%r(:,ibnd)),&
                                        _CX(becpsi%r(:,im)), deexx)
                IF(ibnd<ibnd_end) &
                CALL PAW_newdxx(x2/nqs, _CX(becxx(ikq)%r(:,ibnd+1)), &
                                        _CX(becpsi%r(:,im)), deexx)
             ENDIF
             !
             ! accumulates over bands and k points
             !
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                RESULT(ir) = RESULT(ir)+x1* dble(vc(ir))* dble(exxbuff(ir,h_ibnd,current_k))&
                                       +x2*aimag(vc(ir))*aimag(exxbuff(ir,h_ibnd,current_k))
             ENDDO
!$omp end parallel do
             !
          ENDDO &
          IBND_LOOP_GAM
          !
          IF(okvan) THEN
             CALL mp_sum(deexx,intra_bgrp_comm)
             CALL mp_sum(deexx,inter_bgrp_comm)
          ENDIF
          !
          CALL mp_sum( RESULT(1:nrxxs), inter_bgrp_comm)
          !
          ! brings back result in G-space
          !
          CALL fwfft( 'Wave' , RESULT, dfftt )
          !
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*RESULT(dfftt%nl(ig))
          ENDDO
!$omp end parallel do
          ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
          IF(okvan) CALL add_nlxx_pot (lda, hpsi(:,im), xkp, n, &
                           igk_k(1,current_k), deexx, eps_occ, exxalfa)
       ENDDO &
       LOOP_ON_PSI_BANDS
       IF ( okvan .and..not.tqr ) CALL qvan_clean ()
       !
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !
    DEALLOCATE( RESULT, temppsic_dble, temppsic_aimag)
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_std_gamma
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_std_k(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... generic, k-point version of vexx
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, current_k
    USE klist,          ONLY : xk, nks, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, nbgrp
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : gt, nqs, index_xkq, x_occupation, g2_convolution
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: temppsic(:), RESULT(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:),result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: result_g(:), result_nc_g(:,:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ngmt
    !
    ngmt = dfftt%ngm 
    ALLOCATE( fac(ngmt) )
    nrxxs= dfftt%nnr
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol), result_nc(nrxxs,npol) )
       ALLOCATE( result_nc_g(n,npol) )
    ELSE
       ALLOCATE( temppsic(nrxxs), RESULT(nrxxs) )
       ALLOCATE( result_g(n) )
    ENDIF
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=0.0_DP
       psi =0.0_DP
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    ENDIF
    !
    LOOP_ON_PSI_BANDS : &
    DO im = 1,m !for each band of psi (the k cycle is outside band)
       IF(okvan) deexx = 0._DP
       !
       IF (noncolin) THEN
          temppsic_nc = 0._DP
       ELSE
          temppsic    = 0._DP
       ENDIF
       !
       IF (noncolin) THEN
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(dfftt%nl(igk_k(ig,current_k)),1) = psi(ig,im)
          ENDDO
!$omp end parallel do
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(dfftt%nl(igk_k(ig,current_k)),2) = psi(npwx+ig,im)
          ENDDO
!$omp end parallel do
          !
          CALL invfft ('Wave', temppsic_nc(:,1), dfftt)
          CALL invfft ('Wave', temppsic_nc(:,2), dfftt)
          !
       ELSE
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic( dfftt%nl(igk_k(ig,current_k)) ) = psi(ig,im)
          ENDDO
!$omp end parallel do
          CALL invfft ('Wave', temppsic, dfftt)
          !
       ENDIF
       !
       IF (noncolin) THEN
          result_nc = 0.0_DP
       ELSE
          RESULT    = 0.0_DP
       ENDIF
       !
       INTERNAL_LOOP_ON_Q : &
       DO iq=1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          xkq  = xkq_collect(:,ikq)
          !
          ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
          CALL g2_convolution(ngmt, gt, xkp, xkq, fac)
          IF ( okvan .and..not.tqr ) CALL qvan_init (ngmt, xkq, xkp)
          !
          IBND_LOOP_K : &
          DO ibnd=ibnd_start,ibnd_end !for each band of psi
             !
             IF ( abs(x_occupation(ibnd,ik)) < eps_occ) CYCLE IBND_LOOP_K
             !
             !loads the phi from file
             !
             !   >>>> calculate rho in real space
             IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = ( conjg(exxbuff(ir,ibnd,ikq))*temppsic_nc(ir,1) + &
                                conjg(exxbuff(nrxxs+ir,ibnd,ikq))*temppsic_nc(ir,2) )/omega
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir)=conjg(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
                ENDDO
!$omp end parallel do
             ENDIF
             !   >>>> add augmentation in REAL space HERE
             IF(okvan .and. tqr) THEN ! augment the "charge" in real space
                CALL addusxx_r( rhoc, becxx(ikq)%k(:,ibnd), becpsi%k(:,im))
             ENDIF
             !
             !   >>>> brings it to G-space
             CALL fwfft('Rho', rhoc, dfftt)
             !
             !   >>>> add augmentation in G space HERE
             IF(okvan .and. .not. tqr) THEN
                CALL addusxx_g(dfftt, rhoc, xkq, xkp, 'c', &
                   becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,im))
             ENDIF
             !   >>>> charge done
             !
             vc = 0._DP
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, ngmt
                vc(dfftt%nl(ig)) = fac(ig) * rhoc(dfftt%nl(ig)) * &
                                             x_occupation(ibnd,ik) / nqs
             ENDDO
!$omp end parallel do
             !
             ! Add ultrasoft contribution (RECIPROCAL SPACE)
             ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
             IF(okvan .and. .not. tqr) THEN
                CALL newdxx_g(dfftt, vc, xkq, xkp, 'c', deexx, &
                              becphi_c=becxx(ikq)%k(:,ibnd))
             ENDIF
             !
             !brings back v in real space
             CALL invfft ('Rho', vc, dfftt)
             !
             ! Add ultrasoft contribution (REAL SPACE)
             IF(okvan .and. tqr) CALL newdxx_r(dfftt,vc, becxx(ikq)%k(:,ibnd), deexx)
             !
             ! Add PAW one-center contribution
             IF(okpaw) THEN
                CALL PAW_newdxx(x_occupation(ibnd,ik)/nqs, becxx(ikq)%k(:,ibnd), becpsi%k(:,im), deexx)
             ENDIF
             !
             !accumulates over bands and k points
             !
             IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result_nc(ir,1)= result_nc(ir,1) + vc(ir) * exxbuff(ir,ibnd,ikq)
                ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result_nc(ir,2)= result_nc(ir,2) + vc(ir) * exxbuff(ir+nrxxs,ibnd,ikq)
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   RESULT(ir) = RESULT(ir) + vc(ir)*exxbuff(ir,ibnd,ikq)
                ENDDO
!$omp end parallel do
             ENDIF
             !
          ENDDO &
          IBND_LOOP_K
          IF ( okvan .and..not.tqr ) CALL qvan_clean ()
          !
       ENDDO &
       INTERNAL_LOOP_ON_Q
       !
       IF(okvan) THEN
         CALL mp_sum(deexx,intra_bgrp_comm)
         CALL mp_sum(deexx,inter_bgrp_comm)
       ENDIF
       !
       ! bring result back to G-space
       !
       IF (noncolin) THEN
          !
          CALL fwfft ('Wave', result_nc(:,1), dfftt)
          CALL fwfft ('Wave', result_nc(:,2), dfftt)
          !
          !communicate result
          DO ig = 1, n
             result_nc_g(ig,1:npol) = result_nc(dfftt%nl(igk_k(ig,current_k)),1:npol)
          ENDDO
          CALL mp_sum( result_nc_g(1:n,1:npol), inter_bgrp_comm)
          !
          !adds it to hpsi
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)    = hpsi(ig,im)     - exxalfa*result_nc_g(ig,1)
          ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(lda+ig,im)= hpsi(lda+ig,im) - exxalfa*result_nc_g(ig,2)
          ENDDO
!$omp end parallel do
          !
       ELSE
          !
          CALL fwfft ('Wave', RESULT, dfftt)
          !
          !communicate result
          DO ig = 1, n
             result_g(ig) = RESULT(dfftt%nl(igk_k(ig,current_k)))
          ENDDO
          CALL mp_sum( result_g(1:n), inter_bgrp_comm)
          !
          !adds it to hpsi
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*result_g(ig)
          ENDDO
!$omp end parallel do
       ENDIF
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) CALL add_nlxx_pot(lda, hpsi(:,im), xkp, n, igk_k(:,current_k),&
                                       deexx, eps_occ, exxalfa)
       !
    ENDDO &
    LOOP_ON_PSI_BANDS
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, result_nc, result_nc_g )
    ELSE
       DEALLOCATE(temppsic, RESULT, result_g )
    ENDIF
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx)
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_std_k
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_stress_std()
    !-----------------------------------------------------------------------
    !
    ! This is Eq.(10) of PRB 73, 125120 (2006).
    !
    USE constants,            ONLY : fpi, e2, pi, tpi
    USE io_files,             ONLY : iunwfc, nwordwfc
    USE buffers,              ONLY : get_buffer
    USE cell_base,            ONLY : alat, omega, bg, at, tpiba
    USE symm_base,            ONLY : nsym, s
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE wavefunctions,        ONLY : evc
    USE klist,                ONLY : xk, ngk, nks, igk_k
    USE lsda_mod,             ONLY : lsda, current_spin, isk
    USE gvect,                ONLY : g
    USE mp_pools,             ONLY : npool, inter_pool_comm
    USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE uspp,                 ONLY : okvan
    USE exx_base,             ONLY : exxdiv, erfc_scrlen, gau_scrlen, grid_factor, eps, &
                                     nq1, nq2, nq3, nqs, on_double_grid, use_coulomb_vcut_spheric, &
                                     use_coulomb_vcut_ws, vcut, x_gamma_extrapolation, yukawa, index_xkq, &
                                     index_sym, x_occupation
    USE coulomb_vcut_module,  ONLY : vcut_get,  vcut_spheric_get
    !
    ! ---- local variables -------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! local variables
    REAL(DP)   :: exx_stress_std(3,3), exx_stress_(3,3)
    !
    COMPLEX(DP),ALLOCATABLE :: tempphic(:), temppsic(:), RESULT(:)
    COMPLEX(DP),ALLOCATABLE :: tempphic_nc(:,:), temppsic_nc(:,:), &
                               result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: rhoc(:)
    REAL(DP),    ALLOCATABLE :: fac(:), fac_tens(:,:,:), fac_stress(:)
    INTEGER  :: npw, jbnd, ibnd, ik, ikk, ig, ir, ikq, iq, isym
    INTEGER  :: h_ibnd, nqi, iqi, beta, nrxxs, ngm
    INTEGER  :: ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    ! temp array for vcut_spheric
    REAL(DP) :: delta(3,3)
    INTEGER :: ngmt
    !
    IF (npool>1) CALL errore('exx_stress1','stress not available with pools',1)
    IF (noncolin) CALL errore('exx_stress1','noncolinear stress not implemented',1)
    IF (okvan) CALL infomsg('exx_stress1','USPP stress not tested')
    !
    nrxxs = dfftt%nnr
    ngmt = dfftt%ngm
    ngm   = ngmt
    delta = reshape( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    ALLOCATE( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    ALLOCATE( fac_tens(3,3,ngm), fac_stress(ngm) )
    !
    nqi=nqs
    !
    ! loop over k-points
    DO ikk = 1, nks
        current_k = ikk
        IF (lsda) current_spin = isk(ikk)
        npw = ngk(ikk)

        IF (nks > 1) &
            CALL get_buffer(evc, nwordwfc, iunwfc, ikk)

        ! loop over bands
        DO jbnd = 1, nbnd
            !
            temppsic(:) = ( 0._dp, 0._dp )
!$omp parallel do default(shared), private(ig)
            DO ig = 1, npw
                temppsic(dfftt%nl(igk_k(ig,ikk))) = evc(ig,jbnd)
            ENDDO
!$omp end parallel do
            !
            IF(gamma_only) THEN
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw
                    temppsic(dfftt%nlm(igk_k(ig,ikk))) = conjg(evc(ig,jbnd))
                ENDDO
!$omp end parallel do
            ENDIF

            CALL invfft ('Wave', temppsic, dfftt)

            DO iqi = 1, nqi
                !
                iq=iqi
                !
                ikq  = index_xkq(current_k,iq)
                ik   = index_xk(ikq)
                isym = abs(index_sym(ikq))

                ! FIXME: use cryst_to_cart and company as above..
                xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
                IF (index_sym(ikq) < 0) xk_cryst = -xk_cryst
                sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                         s(:,2,isym)*xk_cryst(2) + &
                         s(:,3,isym)*xk_cryst(3)
                xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

                !CALL start_clock ('exxen2_ngmloop')

!$omp parallel do default(shared), private(ig, beta, q, qq, on_double_grid, x)
                DO ig = 1, ngm
                  q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                  q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                  q(3)= xk(3,current_k) - xkq(3) + g(3,ig)

                  q = q * tpiba
                  qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )

                  DO beta = 1, 3
                      fac_tens(1:3,beta,ig) = q(1:3)*q(beta)
                  ENDDO

                  IF (x_gamma_extrapolation) THEN
                      on_double_grid = .true.
                      x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                  ELSE
                      on_double_grid = .false.
                  ENDIF

                  IF (use_coulomb_vcut_ws) THEN
                      fac(ig) = vcut_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSEIF ( use_coulomb_vcut_spheric ) THEN
                      fac(ig) = vcut_spheric_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSEIF (gau_scrlen > 0) THEN
                      fac(ig)=e2*((pi/gau_scrlen)**(1.5d0))* &
                            exp(-qq/4.d0/gau_scrlen) * grid_factor
                      fac_stress(ig) =  e2*2.d0/4.d0/gau_scrlen * &
                            exp(-qq/4.d0/gau_scrlen) *((pi/gau_scrlen)**(1.5d0))* &
                                                                     grid_factor
                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSEIF (qq > 1.d-8) THEN
                      IF ( erfc_scrlen > 0 ) THEN
                        fac(ig)=e2*fpi/qq*(1._dp-exp(-qq/4.d0/erfc_scrlen**2)) * grid_factor
                        fac_stress(ig) = -e2*fpi * 2.d0/qq**2 * ( &
                            (1._dp+qq/4.d0/erfc_scrlen**2)*exp(-qq/4.d0/erfc_scrlen**2) - 1._dp) * &
                            grid_factor
                      ELSE
                        fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
                        fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2 * grid_factor
                      ENDIF

                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSE
                      fac(ig)= -exxdiv ! or rather something else (see f.gygi)
                      fac_stress(ig) = 0._dp  ! or -exxdiv_stress (not yet implemented)
                      IF ( yukawa> 0._dp .and. .not. x_gamma_extrapolation) THEN
                        fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
                        fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2
                      ENDIF
                      IF (erfc_scrlen > 0._dp .and. .not. x_gamma_extrapolation) THEN
                        fac(ig) = e2*fpi / (4.d0*erfc_scrlen**2)
                        fac_stress(ig) = e2*fpi / (8.d0*erfc_scrlen**4)
                      ENDIF
                  ENDIF
                ENDDO
!$omp end parallel do
                !CALL stop_clock ('exxen2_ngmloop')

                IF (gamma_only) THEN
                    !
                    h_ibnd = ibnd_start/2
                    !
                    IF(mod(ibnd_start,2)==0) THEN
                      h_ibnd=h_ibnd-1
                      ibnd_loop_start=ibnd_start-1
                    ELSE
                      ibnd_loop_start=ibnd_start
                    ENDIF
                    !
                    DO ibnd = ibnd_loop_start, ibnd_end, 2     !for each band of psi
                        !
                        h_ibnd = h_ibnd + 1
                        !
                        IF( ibnd < ibnd_start ) THEN
                            x1 = 0._dp
                        ELSE
                            x1 = x_occupation(ibnd,  ik)
                        ENDIF

                        IF( ibnd == ibnd_end) THEN
                            x2 = 0._dp
                        ELSE
                            x2 = x_occupation(ibnd+1,  ik)
                        ENDIF
                        IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
                        !
                        ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                        DO ir = 1, nrxxs
                            tempphic(ir) = exxbuff(ir,h_ibnd,ikq)
                            rhoc(ir)     = conjg(tempphic(ir))*temppsic(ir) / omega
                        ENDDO
!$omp end parallel do
                        ! bring it to G-space
                        CALL fwfft ('Rho', rhoc, dfftt)

                        vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                        DO ig = 1, ngm
                            !
                            vc(:,:) = vc(:,:) + x1 * 0.25_dp * &
                                      abs( rhoc(dfftt%nl(ig)) + &
                                      conjg(rhoc(dfftt%nlm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                            vc(:,:) = vc(:,:) + x2 * 0.25_dp * &
                                      abs( rhoc(dfftt%nl(ig)) - &
                                      conjg(rhoc(dfftt%nlm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                        ENDDO
!$omp end parallel do
                        vc = vc / nqs / 4.d0
                        exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                    ENDDO

                ELSE

                    DO ibnd = ibnd_start, ibnd_end    !for each band of psi
                      !
                      IF ( abs(x_occupation(ibnd,ik)) < 1.d-6) CYCLE
                      !
                      ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                          tempphic(ir) = exxbuff(ir,ibnd,ikq)
                          rhoc(ir)     = conjg(tempphic(ir))*temppsic(ir) / omega
                      ENDDO
!$omp end parallel do

                      ! bring it to G-space
                      CALL fwfft ('Rho', rhoc, dfftt)

                      vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm
                          vc(:,:) = vc(:,:) + rhoc(dfftt%nl(ig))  * &
                                        conjg(rhoc(dfftt%nl(ig)))* &
                                    (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                      ENDDO
!$omp end parallel do
                      vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)

                    ENDDO

                ENDIF ! gamma or k-points

            ENDDO ! iqi
        ENDDO ! jbnd
    ENDDO ! ikk

    DEALLOCATE(tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
    !
    CALL mp_sum( exx_stress_, intra_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    exx_stress_std = exx_stress_
    !
    !-----------------------------------------------------------------------
  END FUNCTION exx_stress_std
  !-----------------------------------------------------------------------
  !
END MODULE exx_std
