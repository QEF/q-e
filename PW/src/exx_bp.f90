! Copyrigh(C) 2005-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE exx_bp
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only, use_gpu, many_fft, tqr
  USE exx_base,             ONLY : dfftt , exxbuff, exxbuff_d, npwt, x_nbnd_occ, &
                                   gt, ggt, gcutmt, gkcut, gstart_t, ngmt_g, &
                                   eps_occ, exxalfa, x_occupation, x_occupation_d, &
                                   ibnd_start, ibnd_end, ibnd_buff_start, ibnd_buff_end
  USE noncollin_module,     ONLY : noncolin, npol
  USE cell_base,            ONLY : at, bg, tpiba2
  USE gvecw,                ONLY : ecutwfc
  USE symm_base,            ONLY : fft_fact
  USE mp_exx,               ONLY : negrp, inter_egrp_comm, intra_egrp_comm, nproc_egrp
  USE klist,                ONLY : nks, xk
  USE mp,                   ONLY : mp_sum
  !
  INTEGER, ALLOCATABLE :: ig_l2gt(:), millt(:,:)
  !
  REAL(DP), ALLOCATABLE :: coulomb_fac(:,:,:)
  !! the Coulomb factor is reused between iterations
  !
  LOGICAL, ALLOCATABLE :: coulomb_done(:,:)
  !! list of which Coulomb factors have been calculated already
  !
 CONTAINS
#define _CX(A)  CMPLX(A,0._dp,kind=DP)
#define _CY(A)  CMPLX(0._dp,-A,kind=DP)
  !
  !-----------------------------------------------------------------------
  SUBROUTINE g2_convolution_all( ngm, g, xk, xkq, iq, current_k )
    !-----------------------------------------------------------------------
    !! Wrapper for g2_convolution.
    !
    USE kinds,     ONLY : DP
    USE klist,     ONLY : nks
    USE exx_base,  ONLY : g2_convolution, nqs
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: ngm
    !! Number of G vectors
    REAL(DP), INTENT(IN) :: g(3,ngm)
    !! Cartesian components of G vectors
    REAL(DP), INTENT(IN) :: xk(3)
    !! current k vector
    REAL(DP), INTENT(IN) :: xkq(3)
    !! current q vector
    INTEGER, INTENT(IN) :: current_k
    !! current k-point index
    INTEGER, INTENT(IN) :: iq
    !! q-grid point index
    !
    ! ... Check if coulomb_fac has been allocated
    IF( .NOT. ALLOCATED( coulomb_fac ) ) ALLOCATE( coulomb_fac(ngm,nqs,nks) )
    !
    ! ... Check if coulomb_done has been allocated
    IF( .NOT. ALLOCATED( coulomb_done) ) THEN
       ALLOCATE( coulomb_done(nqs,nks) )
       coulomb_done = .FALSE.
    ENDIF
    !
    ! ... return if this k and k' already computed, otherwise compute it
    IF ( coulomb_done(iq,current_k) ) RETURN
    !
    CALL g2_convolution( ngm, g, xk, xkq, coulomb_fac(:,iq,current_k) )
    !
    coulomb_done(iq,current_k) = .TRUE.
    !
  END SUBROUTINE g2_convolution_all
  !
  !------------------------------------------------------------------------
  SUBROUTINE set_dfftt_grid_bp( )
    !------------------------------------------------------------------------
    USE exx_bp_utils,         ONLY : smap_exx
    USE command_line_options, ONLY : nmany_, pencil_decomposition_
    USE mp_bands,             ONLY : nyfft
    USE symm_base,            ONLY : fft_fact
    USE fft_types,            ONLY : fft_type_init
    USE recvec_subs,          ONLY : ggen
    !
    IMPLICIT NONE
       LOGICAL :: lpara
       INTEGER :: ngmt
    INTEGER, EXTERNAL :: n_plane_waves
    !
    WRITE( 6, "(5X,'Exchange parallelized over bands (',i4,' band groups)')" ) &
           negrp
    lpara = ( nproc_egrp > 1 )
    CALL fft_type_init( dfftt, smap_exx, "rho", gamma_only, lpara,     &
                        intra_egrp_comm, at, bg, gcutmt, gcutmt/gkcut, &
                        fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_,  &
                        use_pd=pencil_decomposition_ )
    ngmt = dfftt%ngm
    ngmt_g = ngmt
    CALL mp_sum( ngmt_g, intra_egrp_comm )
    ALLOCATE( gt(3,dfftt%ngm) )
    ALLOCATE( ggt(dfftt%ngm)  )
    ALLOCATE( millt(3,dfftt%ngm) )
    ALLOCATE( ig_l2gt(dfftt%ngm) )
    !
    CALL ggen( dfftt, gamma_only, at, bg, gcutmt, ngmt_g, ngmt, &
               gt, ggt, millt, ig_l2gt, gstart_t )
    !
    DEALLOCATE( ig_l2gt )
    DEALLOCATE( millt )
    npwt = n_plane_waves( ecutwfc/tpiba2, nks, xk, gt, ngmt )
    !
    RETURN
    !
  END SUBROUTINE set_dfftt_grid_bp
  !
  !------------------------------------------------------------------------
  SUBROUTINE exxinit_bp( )
    !------------------------------------------------------------------------
    !! This subroutine is run before the first H_psi() of each iteration. 
    !! It saves the wavefunctions for the right density matrix, in real space.
    !
    USE wavefunctions,        ONLY : psic, evc
    USE io_files,             ONLY : iunwfc_exx
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx
    USE klist,                ONLY : ngk, nks, nkstot
    USE mp_exx,               ONLY : me_egrp, init_index_over_band,  &
                                     inter_egrp_comm,           &
                                     iexx_start, iexx_end, &
                                     all_start, all_end
    USE mp,                   ONLY : mp_bcast
    USE scatter_mod,          ONLY : gather_grid, scatter_grid
    USE fft_interfaces,       ONLY : invfft
    USE mp_orthopools,        ONLY : intra_orthopool_comm
    USE exx_base,             ONLY : nkqs, index_xk, index_sym,  &
                                     rir, working_pool, d_spin
    USE exx_bp_utils,         ONLY : change_data_structure, nwordwfc_exx, &
                                     igk_exx, evc_exx, transform_evc_to_exx
#if defined(__CUDA)
    USE device_memcpy_m,      ONLY : dev_memset
    USE device_fbuff_m,       ONLY : dev_buf
#endif
    !
    IMPLICIT NONE
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
    INTEGER :: nxxs, nrxxs
#if defined(__MPI)
    COMPLEX(DP),ALLOCATABLE  :: temppsic_all(:), psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    INTEGER :: npw, current_ik
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ibnd_start_new, ibnd_end_new, max_buff_bands_per_egrp
    INTEGER :: ibnd_exx, evc_offset
    !
    !$acc update device(evc)
    CALL transform_evc_to_exx( 2 )
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
    CALL init_index_over_band( inter_egrp_comm, nbnd, nbnd )
    !
    ! ... this will cause exxbuff to be calculated for every band
    ibnd_start_new = iexx_start
    ibnd_end_new = iexx_end
    !
    IF ( gamma_only ) THEN
        ibnd_buff_start = (ibnd_start_new+1)/2
        ibnd_buff_end   = (ibnd_end_new+1)/2
        max_buff_bands_per_egrp = MAXVAL((all_end(:)+1)/2-(all_start(:)+1)/2)+1
    ELSE
        ibnd_buff_start = ibnd_start_new
        ibnd_buff_end   = ibnd_end_new
        max_buff_bands_per_egrp = MAXVAL(all_end(:)-all_start(:))+1
    ENDIF
    !
    IF (.NOT. ALLOCATED(exxbuff)) THEN
       IF (gamma_only) THEN
          ALLOCATE( exxbuff(nrxxs*npol,ibnd_buff_start:ibnd_buff_start + &
                                        max_buff_bands_per_egrp-1,nkqs) ) ! THIS WORKS as for k
       ELSE
          ALLOCATE( exxbuff(nrxxs*npol,ibnd_buff_start:ibnd_buff_start + &
                                        max_buff_bands_per_egrp-1,nkqs) )
       ENDIF
    ENDIF
    !
    IF (.not. allocated(exxbuff_d) .and. use_gpu) THEN
       IF (gamma_only) THEN
          ALLOCATE( exxbuff_d(nrxxs*npol, ibnd_buff_start:ibnd_buff_start+max_buff_bands_per_egrp-1, nks))
       ELSE
          ALLOCATE( exxbuff_d(nrxxs*npol, ibnd_buff_start:ibnd_buff_start+max_buff_bands_per_egrp-1, nkqs))
       END IF
    ENDIF
    !
    IF (use_gpu) THEN
#if defined (__CUDA)
       ! NB: the array bounds are not passed to the subroutine.
       !
       ! See https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/269311
       !
       ! NB: TO BE CORRECTED WITH THE NEW DeviceXlib LIBRARY that dues internal slicing right!
       CALL dev_memset(exxbuff_d, (0.0_DP,0.0_DP), &
                                 (/ 1,nrxxs*npol/), 1, &
                                 (/ ibnd_buff_start, ibnd_buff_end /), ibnd_buff_start, &
                                 (/ 1,SIZE(exxbuff_d,3)/), 1)
#endif
    ELSE
       !$omp parallel do collapse(3) default(shared) firstprivate(npol,nrxxs,nkqs, &
       !$omp                ibnd_buff_start,ibnd_buff_end) private(ir,ibnd,ikq,ipol)
       DO ikq = 1, SIZE(exxbuff,3) 
          DO ibnd = ibnd_buff_start, ibnd_buff_end
             DO ir = 1, nrxxs*npol
                exxbuff(ir,ibnd,ikq) = (0.0_DP,0.0_DP)
             ENDDO
          ENDDO
       ENDDO
       ! the above loops will replaced with the following line soon
       !CALL threaded_memset(exxbuff, 0.0_DP, nrxxs*npol*SIZE(exxbuff,2)*nkqs*2)
    ENDIF
    !
    ! ... This is parallelized over pools. Each pool computes only its k-points
    !
    KPOINTS_LOOP : &
    DO ik = 1, nks
       !
       IF ( nks > 1 ) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ik )
       !
       ! ik         = index of k-point in this pool
       ! current_ik = index of k-point over all pools
       !
       current_ik = global_kpoint_index( nkstot, ik )
       !
       IF_GAMMA_ONLY : &
       IF (gamma_only) THEN
          !
          IF (MOD(iexx_start,2) == 0) THEN
             ibnd_loop_start = iexx_start-1
          ELSE
             ibnd_loop_start = iexx_start
          ENDIF
          !
          evc_offset = 0
          DO ibnd = ibnd_loop_start, iexx_end, 2
             !
             psic(:) = ( 0._DP, 0._DP )
             !
             IF ( ibnd < iexx_end ) THEN
                IF ( ibnd == ibnd_loop_start .AND. MOD(iexx_start,2) == 0 ) THEN
                   DO ig = 1, npwt
                      psic(dfftt%nl(ig))  = ( 0._DP, 1._DP )*evc_exx(ig,1)
                      psic(dfftt%nlm(ig)) = ( 0._DP, 1._DP )*CONJG(evc_exx(ig,1))
                   ENDDO
                   evc_offset = -1
                ELSE
                   DO ig = 1, npwt
                      psic(dfftt%nl(ig))  = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) &
                           + ( 0._DP, 1._DP ) * evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2)
                      psic(dfftt%nlm(ig)) = CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) ) &
                           + ( 0._DP, 1._DP ) * CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2) )
                   ENDDO
                ENDIF
             ELSE
                DO ig=1,npwt
                   psic(dfftt%nl (ig)) = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1)
                   psic(dfftt%nlm(ig)) = CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) )
                ENDDO
             ENDIF
             !
             CALL invfft( 'Wave', psic, dfftt )
             !
             exxbuff(1:nrxxs,(ibnd+1)/2,current_ik)=psic(1:nrxxs) 
             !
          ENDDO
          !
       ELSE IF_GAMMA_ONLY
          !
          npw = ngk (ik)
          IBND_LOOP_K : &
          DO ibnd = iexx_start, iexx_end
             !
             ibnd_exx = ibnd
             IF (noncolin) THEN
!$omp parallel do default(shared) private(ir) firstprivate(nrxxs)
                DO ir = 1, nrxxs
                   temppsic_nc(ir,1) = ( 0._DP, 0._DP )
                   temppsic_nc(ir,2) = ( 0._DP, 0._DP )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig = 1, npw
                   temppsic_nc(dfftt%nl(igk_exx(ig,ik)),1) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft( 'Wave', temppsic_nc(:,1), dfftt )
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx,npwx)
                DO ig = 1, npw
                   temppsic_nc(dfftt%nl(igk_exx(ig,ik)),2) = evc_exx(ig+npwx,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft( 'Wave', temppsic_nc(:,2), dfftt )
             ELSE
!$omp parallel do default(shared) private(ir) firstprivate(nrxxs)
                DO ir = 1, nrxxs
                   temppsic(ir) = ( 0._DP, 0._DP )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig = 1, npw
                   temppsic(dfftt%nl(igk_exx(ig,ik))) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
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
!$omp parallel do collapse(2)
                      DO ipol = 1, npol
                         DO ir = 1, nxxs
                            psic_all_nc(ir,ipol) = (0.0_DP, 0.0_DP)
                            DO jpol = 1, npol
                               psic_all_nc(ir,ipol) = psic_all_nc(ir,ipol) + &
                                             CONJG(d_spin(jpol,ipol,isym)) * &
                                             temppsic_all_nc(rir(ir,isym),jpol)
                            ENDDO
                         ENDDO
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   DO ipol = 1, npol
                      CALL scatter_grid( dfftt, psic_all_nc(:,ipol), psic_nc(:,ipol) )
                   ENDDO
#else
!$omp parallel do collapse(2)
                   DO ipol = 1, npol
                      DO ir = 1, nxxs
                         psic_nc(ir,ipol) = (0._DP,0._DP)
                         DO jpol = 1, npol
                            psic_nc(ir,ipol) = psic_nc(ir,ipol) + CONJG(d_spin(jpol,ipol,isym))* &
                                               temppsic_nc(rir(ir,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDDO
!$omp end parallel do
#endif
                   !
#if defined (__CUDA)
                   IF (use_gpu) CALL dev_buf%lock_buffer(psic_nc_d, (/nrxxs, npol/), ierr)
                   IF (use_gpu) psic_nc_d = psic_nc
#endif
                   !
                   IF (index_sym(ikq) > 0 ) THEN
                      IF (use_gpu) THEN
                         associate(exxbuff=>exxbuff_d, psic_nc=>psic_nc_d)
                         ! sym. op. without time reversal: normal case
                         !$cuf kernel do 
                         DO ir=1,nrxxs
                            exxbuff(ir,ibnd,ikq)=psic_nc(ir,1)
                            exxbuff(ir+nrxxs,ibnd,ikq)=psic_nc(ir,2)
                         ENDDO
                         end associate
                      ELSE
                         ! sym. op. without time reversal: normal case
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,ikq)
                      DO ir = 1, nrxxs
                         exxbuff(ir,ibnd,ikq) = psic_nc(ir,1)
                         exxbuff(ir+nrxxs,ibnd,ikq) = psic_nc(ir,2)
                      ENDDO
!$omp end parallel do
                      END IF
                   ELSE
                      ! sym. op. with time reversal: spin 1->2*, 2->-1*
                      IF (use_gpu) THEN
                         associate(exxbuff=>exxbuff_d, psic_nc=>psic_nc_d)
                         ! sym. op. with time reversal: spin 1->2*, 2->-1*
                         !$cuf kernel do 
                         DO ir=1,nrxxs
                            exxbuff(ir,ibnd,ikq)=CONJG(psic_nc(ir,2))
                            exxbuff(ir+nrxxs,ibnd,ikq)=-CONJG(psic_nc(ir,1))
                         ENDDO
                         end associate
                      ELSE
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,ikq)
                      DO ir = 1, nrxxs
                         exxbuff(ir,ibnd,ikq) = CONJG(psic_nc(ir,2))
                         exxbuff(ir+nrxxs,ibnd,ikq) = -CONJG(psic_nc(ir,1))
                      ENDDO
!$omp end parallel do
                      ENDIF
                   ENDIF
#if defined(__CUDA)
                IF (use_gpu) CALL dev_buf%release_buffer(psic_nc_d, ierr)
                IF (use_gpu) exxbuff = exxbuff_d
#endif
                ELSE ! noncolinear
#if defined(__MPI)
                   CALL gather_grid( dfftt, temppsic, temppsic_all )
                   IF ( me_egrp == 0 ) THEN
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                      DO ir = 1, nxxs
                         psic_all(ir) = temppsic_all(rir(ir,isym))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   CALL scatter_grid( dfftt, psic_all, psic )
#else
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                   DO ir = 1, nrxxs
                      psic(ir) = temppsic(rir(ir,isym))
                   ENDDO
!$omp end parallel do
#endif
!$omp parallel do default(shared) private(ir) firstprivate(isym,ibnd,ikq)
                   DO ir = 1, nrxxs
                      IF (index_sym(ikq) < 0 ) THEN
                         psic(ir) = CONJG(psic(ir))
                      ENDIF
                      exxbuff(ir,ibnd,ikq) = psic(ir)
                   ENDDO
!$omp end parallel do
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
    CALL change_data_structure( .FALSE. )
    !
    RETURN
    !
  END SUBROUTINE exxinit_bp
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_bp( lda, n, m, psi, hpsi, becpsi )
    !-----------------------------------------------------------------------
    !! Wrapper routine computing V_x\psi, V_x = exchange potential. 
    !! Calls generic version vexx_k or Gamma-specific one vexx_gamma.
    !
    USE becmod,         ONLY : bec_type
    USE mp_exx,         ONLY : negrp, inter_egrp_comm, init_index_over_band
    USE wvfct,          ONLY : nbnd
    USE exx_bp_utils,   ONLY : transform_psi_to_exx, transform_hpsi_to_local, &
                               psi_exx, hpsi_exx
    !
    IMPLICIT NONE
    !
    INTEGER :: lda
    !! input: leading dimension of arrays psi and hpsi
    INTEGER :: n
    !! input: true dimension of psi and hpsi
    INTEGER :: m
    !! input: number of states psi
    COMPLEX(DP) :: psi(lda*npol,m)
    !! input: m wavefunctions
    COMPLEX(DP) :: hpsi(lda*npol,m)
    !! output: V_x*psi
    TYPE(bec_type), OPTIONAL :: becpsi
    !! input: <beta|psi>, optional but needed for US and PAW case
    !
    IF (negrp > 1) THEN
       CALL init_index_over_band( inter_egrp_comm, nbnd, m )
       !
       ! ... transform psi to the EXX data structure
       CALL transform_psi_to_exx( lda, n, m, psi )
    ENDIF
    !
    ! ... calculate the EXX contribution to hpsi
    !
    IF ( gamma_only ) THEN
       IF (negrp == 1)THEN
          IF (.not. use_gpu) CALL vexx_bp_gamma( lda, n, m, psi, hpsi, becpsi )
          IF (      use_gpu) CALL vexx_bp_gamma_gpu( lda, n, m, psi, hpsi, becpsi )
       ELSE
          IF (.not. use_gpu) CALL vexx_bp_gamma( lda, n, m, psi_exx, hpsi_exx, becpsi )
          IF (      use_gpu) CALL vexx_bp_gamma_gpu( lda, n, m, psi_exx, hpsi_exx, becpsi )
       ENDIF
    ELSE
       IF (negrp == 1)THEN
          IF (.not. use_gpu) CALL vexx_bp_k( lda, n, m, psi, hpsi, becpsi )
          IF (      use_gpu) CALL vexx_bp_k_gpu( lda, n, m, psi, hpsi, becpsi )
       ELSE
          IF (.not. use_gpu) CALL vexx_bp_k( lda, n, m, psi_exx, hpsi_exx, becpsi )
          IF (      use_gpu) CALL vexx_bp_k_gpu( lda, n, m, psi_exx, hpsi_exx, becpsi )
       ENDIF
    ENDIF
    !
    IF (negrp > 1) THEN
       !
       ! ... transform hpsi to the local data structure
       !
       CALL transform_hpsi_to_local(lda,n,m,hpsi)
       !
    ENDIF
    !
  END SUBROUTINE vexx_bp
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_bp_gamma( lda, n, m, psi, hpsi, becpsi )
    !-----------------------------------------------------------------------
    !! Gamma-specific version of vexx.
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, current_k, nbnd
    USE klist,          ONLY : xk, nks, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,         ONLY : inter_egrp_comm, my_egrp_id, &
                               intra_egrp_comm, me_egrp, &
                               negrp, max_pairs, egrp_pairs, ibands, nibands, &
                               iexx_istart, iexx_iend, &
                               all_start, all_end, iexx_start, jblock, max_ibands
    USE mp,             ONLY : mp_sum, mp_barrier, mp_circular_shift_left
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : nqs, index_xkq, index_xk, xkq_collect
    USE exx_bp_utils,   ONLY : result_sum, igk_exx
    !
    IMPLICIT NONE
    !
    INTEGER :: lda
    !! input: leading dimension of arrays psi and hpsi
    INTEGER :: n
    !! input: true dimension of psi and hpsi
    INTEGER :: m
    !! input: number of states psi
    COMPLEX(DP) :: psi(lda*npol,max_ibands)
    !! input: m wavefunctions
    COMPLEX(DP) :: hpsi(lda*npol,max_ibands)
    !! output: V_x*psi
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !! input: <beta|psi>, optional but needed for US and PAW case
    !
    ! ... local variables
    !
    COMPLEX(DP), ALLOCATABLE :: result(:,:)
    REAL(DP), ALLOCATABLE :: temppsic_dble (:)
    REAL(DP), ALLOCATABLE :: temppsic_aimag(:)
    !
    COMPLEX(DP), ALLOCATABLE :: vc(:), deexx(:,:)
    INTEGER :: ibnd, ik, im , ikq, iq, ipol
    INTEGER :: ir, ig
    INTEGER :: current_ik
    INTEGER :: ibnd_loop_start
    INTEGER :: nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ialloc
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: iproc, nproc_egrp, ii, ipair
    INTEGER :: jbnd, jstart, jend
    ! ... scratch space for fft of psi and rho
    COMPLEX(DP), ALLOCATABLE :: psi_rhoc_work(:)
    INTEGER :: jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    INTEGER :: exxbuff_index
    INTEGER :: ending_im
    !
    ialloc = nibands(my_egrp_id+1)
    nrxxs = dfftt%nnr
    !
    !ALLOCATE( result(nrxxs), temppsic_DBLE(nrxxs), temppsic_aimag(nrxxs) )
    ALLOCATE( result(nrxxs,ialloc), temppsic_DBLE(nrxxs) )
    ALLOCATE( temppsic_aimag(nrxxs) )
    ALLOCATE( psi_rhoc_work(nrxxs) )
    !
    ALLOCATE( vc(nrxxs) )
    IF (okvan) ALLOCATE( deexx(nkb,ialloc) )
    !
    current_ik = global_kpoint_index( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ALLOCATE( big_result(n,m) )
    big_result = 0.0_DP
    result = 0.0_DP
    !
    DO ii = 1, nibands(my_egrp_id+1)
       IF (okvan) deexx(:,ii) = 0.0_DP
    ENDDO
    !
    ! Here the loops start
    !
    INTERNAL_LOOP_ON_Q : &
    DO iq = 1, nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, current_k )
       IF ( okvan .AND..NOT.tqr ) CALL qvan_init( dfftt%ngm, xkq, xkp )
       !
       DO iegrp = 1, negrp
          !
          ! compute the id of group whose data is currently worked on
          wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
          !
          jblock_start = all_start(wegrp)
          jblock_end   = all_end(wegrp)
          !
          LOOP_ON_PSI_BANDS : &
          DO ii = 1,  nibands(my_egrp_id+1)
             !
             ibnd = ibands(ii,my_egrp_id+1)
             !
             IF (ibnd==0 .OR. ibnd>m) CYCLE
             !
             IF (MOD(ii,2) == 1) THEN
                !
                psi_rhoc_work = (0._DP,0._DP)
                !
                IF ((ii+1) <= MIN(m,nibands(my_egrp_id+1))) THEN
                   ! deal with double bands
!$omp parallel do  default(shared), private(ig)
                   DO ig = 1, npwt
                      psi_rhoc_work( dfftt%nl(ig) )  =       psi(ig, ii) + (0._DP,1._DP) * psi(ig, ii+1)
                      psi_rhoc_work( dfftt%nlm(ig) ) = CONJG(psi(ig, ii) - (0._DP,1._DP) * psi(ig, ii+1))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                IF ( ii == MIN(m,nibands(my_egrp_id+1)) ) THEN
                   ! deal with a single last band
!$omp parallel do  default(shared), private(ig)
                   DO ig = 1, npwt
                      psi_rhoc_work( dfftt%nl(ig) )  =       psi(ig,ii)
                      psi_rhoc_work( dfftt%nlm(ig) ) = CONJG(psi(ig,ii))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                CALL invfft( 'Wave', psi_rhoc_work, dfftt )
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   temppsic_DBLE(ir)  = DBLE( psi_rhoc_work(ir) )
                   temppsic_aimag(ir) = AIMAG( psi_rhoc_work(ir) )
                ENDDO
!$omp end parallel do
                !
             ENDIF
             !
             !
             !determine which j-bands to calculate
             jstart = 0
             jend = 0
             DO ipair = 1, max_pairs
                IF (egrp_pairs(1,ipair,my_egrp_id+1) == ibnd) THEN
                   IF (jstart == 0)THEN
                      jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                      jend = jstart
                   ELSE
                      jend = egrp_pairs(2,ipair,my_egrp_id+1)
                   ENDIF
                ENDIF
             ENDDO
             !
             jstart = MAX(jstart,jblock_start)
             jend = MIN(jend,jblock_end)
             !
             IF (MOD(jstart,2) ==0 ) THEN
                ibnd_loop_start = jstart-1
             ELSE
                ibnd_loop_start = jstart
             ENDIF
             !
             IBND_LOOP_GAM : &
             DO jbnd = ibnd_loop_start, jend, 2 !for each band of psi
                !
                exxbuff_index = (jbnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2
                !
                IF ( jbnd < jstart ) THEN
                   x1 = 0.0_DP
                ELSE
                   x1 = x_occupation(jbnd,ik)
                ENDIF
                IF ( jbnd == jend) THEN
                   x2 = 0.0_DP
                ELSE
                   x2 = x_occupation(jbnd+1,ik)
                ENDIF
                IF ( ABS(x1) < eps_occ .AND. ABS(x2) < eps_occ ) CYCLE
                !
                ! ... calculate rho in real space. Gamma tricks are used.
                ! temppsic is real; tempphic contains one band in the real part,
                ! another one in the imaginary part; the same applies to rhoc
                !
                IF ( MOD(ii,2) == 0 ) THEN
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      psi_rhoc_work(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_aimag(ir) / omega
                   ENDDO
!$omp end parallel do
                ELSE
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      psi_rhoc_work(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_DBLE(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                ! ... bring rho to G-space
                !
                !   >>> add augmentation in REAL SPACE here
                IF (okvan .AND. tqr) THEN
                   IF (jbnd >= jstart) &
                        CALL addusxx_r( psi_rhoc_work, &
                       _CX(becxx(ikq)%r(:,jbnd)), _CX(becpsi%r(:,ibnd)))
                   IF (jbnd < jend) &
                        CALL addusxx_r( psi_rhoc_work, &
                       _CY(becxx(ikq)%r(:,jbnd+1)),_CX(becpsi%r(:,ibnd)))
                ENDIF
                !
                CALL fwfft( 'Rho', psi_rhoc_work, dfftt )
                !   >>>> add augmentation in G SPACE here
                IF (okvan .AND. .NOT. tqr) THEN
                   ! ... contribution from one band added to real (in real space) part of rhoc
                   IF (jbnd >= jstart) &
                        CALL addusxx_g( dfftt, psi_rhoc_work, xkq,  xkp, 'r', &
                        becphi_r=becxx(ikq)%r(:,jbnd), becpsi_r=becpsi%r(:,ibnd) )
                   ! ... contribution from following band added to imaginary (in real space) part of rhoc
                   IF (jbnd < jend) &
                        CALL addusxx_g( dfftt, psi_rhoc_work, xkq,  xkp, 'i', &
                        becphi_r=becxx(ikq)%r(:,jbnd+1), becpsi_r=becpsi%r(:,ibnd) )
                ENDIF
                !   >>>> charge density done
                !
                vc = 0._DP
                !
!$omp parallel do default(shared), private(ig)
                DO ig = 1, dfftt%ngm
                   !
                   vc(dfftt%nl(ig))  = coulomb_fac(ig,iq,current_k) * psi_rhoc_work(dfftt%nl(ig))
                   vc(dfftt%nlm(ig)) = coulomb_fac(ig,iq,current_k) * psi_rhoc_work(dfftt%nlm(ig))
                   !
                ENDDO
!$omp end parallel do
                !
                !   >>>>  compute <psi|H_fock G SPACE here
                IF (okvan .AND. .NOT. tqr) THEN
                   IF (jbnd >= jstart) &
                        CALL newdxx_g( dfftt, vc, xkq, xkp, 'r', deexx(:,ii), &
                           becphi_r=x1*becxx(ikq)%r(:,jbnd) )
                   IF (jbnd<jend) &
                        CALL newdxx_g( dfftt, vc, xkq, xkp, 'i', deexx(:,ii), &
                            becphi_r=x2*becxx(ikq)%r(:,jbnd+1) )
                ENDIF
                !
                !brings back v in real space
                CALL invfft( 'Rho', vc, dfftt )
                !
                !   >>>>  compute <psi|H_fock REAL SPACE here
                IF (okvan .AND. tqr) THEN
                   IF (jbnd >= jstart) &
                        CALL newdxx_r( dfftt,vc, _CX(x1*becxx(ikq)%r(:,jbnd)), deexx(:,ii) )
                   IF (jbnd < jend) &
                        CALL newdxx_r( dfftt,vc, _CY(x2*becxx(ikq)%r(:,jbnd+1)), deexx(:,ii) )
                ENDIF
                !
                IF (okpaw) THEN
                   IF (jbnd >= jstart) &
                        CALL PAW_newdxx( x1/nqs, _CX(becxx(ikq)%r(:,jbnd)),&
                                                _CX(becpsi%r(:,ibnd)), deexx(:,ii) )
                   IF (jbnd < jend) &
                        CALL PAW_newdxx( x2/nqs, _CX(becxx(ikq)%r(:,jbnd+1)),&
                                                _CX(becpsi%r(:,ibnd)), deexx(:,ii) )
                ENDIF
                !
                ! ... accumulates over bands and k points
                !
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result(ir,ii) = result(ir,ii) &
                                 + x1* DBLE(vc(ir))* DBLE(exxbuff(ir,exxbuff_index,ikq)) &
                                 + x2*AIMAG(vc(ir))*AIMAG(exxbuff(ir,exxbuff_index,ikq))
                ENDDO
!$omp end parallel do
                !
             ENDDO &
             IBND_LOOP_GAM
             !
          ENDDO &
          LOOP_ON_PSI_BANDS
          !
          ! get the next nbnd/negrp data
          IF (negrp > 1) CALL mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
          !
       ENDDO ! iegrp
       IF (okvan .AND. .NOT.tqr) CALL qvan_clean()
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !
    DO ii = 1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd == 0 .OR. ibnd>m) CYCLE
       !
       IF (okvan) THEN
          CALL mp_sum( deexx(:,ii), intra_egrp_comm )
       ENDIF
       !
       ! ... brings back result in G-space
       !
       CALL fwfft( 'Wave' , result(:,ii), dfftt )
       !
       ! ... communicate result
       DO ig = 1, n
          big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa*result(dfftt%nl(igk_exx(ig,current_k)),ii)
       ENDDO
       !
       ! ... add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF (okvan) CALL add_nlxx_pot( lda, big_result(:,ibnd), xkp, n, &
                                     igk_exx(1,current_k), deexx(:,ii), eps_occ, exxalfa )
    ENDDO
    !
    CALL result_sum( n*npol, m, big_result )
    !
    IF (iexx_istart(my_egrp_id+1) > 0) THEN
       IF (negrp == 1) THEN
          ending_im = m
       ELSE
          ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
       END IF
       DO im = 1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
           DO ig = 1, n
              hpsi(ig,im) = hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
           ENDDO
!$omp end parallel do
       ENDDO
    ENDIF
    !
    DEALLOCATE( big_result )
    DEALLOCATE( result, temppsic_dble, temppsic_aimag )
    DEALLOCATE( psi_rhoc_work )
    DEALLOCATE( vc )
    IF (okvan) DEALLOCATE( deexx )
    !
  END SUBROUTINE vexx_bp_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_bp_gamma_gpu(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... Gamma-specific version of vexx
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, current_k, nbnd
    USE klist,          ONLY : xk, nks, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,         ONLY : inter_egrp_comm, my_egrp_id, &
                               intra_egrp_comm, me_egrp, &
                               negrp, max_pairs, egrp_pairs, ibands, nibands, &
                               iexx_istart, iexx_iend, &
                               all_start, all_end, iexx_start, jblock, max_ibands
    USE mp,             ONLY : mp_sum, mp_barrier, mp_circular_shift_left
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : nqs, index_xkq, index_xk, xkq_collect
    USE exx_bp_utils,   ONLY : result_sum, igk_exx, igk_exx_d
#if defined(__CUDA)
    USE device_memcpy_m, ONLY : dev_memset
#endif
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,max_ibands)
    COMPLEX(DP)              :: hpsi(lda*npol,max_ibands)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP), ALLOCATABLE :: psi_d(:,:)
#if defined(__CUDA)
    attributes(DEVICE)       :: psi_d
#endif
    COMPLEX(DP),ALLOCATABLE :: result_d(:,:)
#if defined(__CUDA)
    attributes(DEVICE)       :: result_d
#endif
    REAL(DP),ALLOCATABLE :: temppsic_dble_d (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag_d(:)
#if defined(__CUDA)
    attributes(DEVICE)   :: temppsic_dble_d, temppsic_aimag_d
#endif
    !
    COMPLEX(DP),ALLOCATABLE :: vc(:), deexx(:,:), vc_d(:)
    REAL(DP),   ALLOCATABLE :: fac_d(:)
#if defined(__CUDA)
    attributes(DEVICE)   :: vc_d, fac_d
#endif
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ialloc
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    COMPLEX(DP), ALLOCATABLE :: big_result_d(:,:)
#if defined(__CUDA)
    attributes(DEVICE) :: big_result_d
    attributes(PINNED) :: big_result
#endif
    INTEGER :: iproc, nproc_egrp, ii, ipair
    INTEGER :: jbnd, jstart, jend
    ! scratch space for fft of psi and rho
    COMPLEX(DP), ALLOCATABLE :: psi_rhoc_work(:)
    COMPLEX(DP), ALLOCATABLE :: psi_rhoc_work_d(:)
#if defined(__CUDA)
    attributes(DEVICE) :: psi_rhoc_work_d
#endif
    INTEGER :: jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    INTEGER :: exxbuff_index
    INTEGER :: ending_im
    !hack around PGI bug
    INTEGER, POINTER :: dfftt__nl(:)
    INTEGER, POINTER :: dfftt__nlm(:)
#if defined(__CUDA)
    attributes(DEVICE) :: dfftt__nl
    attributes(DEVICE) :: dfftt__nlm
#endif
    !
    ! CUDA Sync
    dfftt__nl=>dfftt%nl_d
    dfftt__nlm=>dfftt%nlm_d
    ALLOCATE(psi_d, source=psi)
    !
    !initial copy of exxbuff
    exxbuff_d = exxbuff
    !
    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac_d(dfftt%ngm) )
    nrxxs= dfftt%nnr
    !
    ALLOCATE( result_d(nrxxs,ialloc))
    ALLOCATE( temppsic_dble_d(nrxxs) )
    ALLOCATE( temppsic_aimag_d(nrxxs) )
    !
    ALLOCATE( psi_rhoc_work(nrxxs) )
    ALLOCATE( psi_rhoc_work_d(nrxxs) )
    !
    ALLOCATE( vc(nrxxs))
    ALLOCATE( vc_d(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb,ialloc))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    allocate(big_result(n,m))
    allocate(big_result_d(n,m))
    big_result = 0.0_DP
#if defined(__CUDA)
    CALL dev_memset(big_result_d,  (0.0_DP, 0.0_DP))
    CALL dev_memset(result_d,  (0.0_DP, 0.0_DP))
#endif
    !
    DO ii=1, nibands(my_egrp_id+1)
       IF(okvan) deexx(:,ii) = 0.0_DP
    END DO
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
       CALL g2_convolution_all(dfftt%ngm, gt, xkp, xkq, iq, current_k)
       IF ( okvan .and..not.tqr ) CALL qvan_init (dfftt%ngm, xkq, xkp)
       !
       ! copy coulomb_fac to device
       fac_d(:) = coulomb_fac(:,iq,current_k)
       !
       DO iegrp=1, negrp
          !
          ! compute the id of group whose data is currently worked on
          wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
          !
          jblock_start = all_start(wegrp)
          jblock_end   = all_end(wegrp)
          !
          LOOP_ON_PSI_BANDS : &
          DO ii = 1,  nibands(my_egrp_id+1)
             !
             ibnd = ibands(ii,my_egrp_id+1)
             !
             IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
             !
             IF ( mod(ii,2)==1 ) THEN
                !
#if defined(__CUDA)
                CALL dev_memset(psi_rhoc_work_d, (0.0_DP,0.0_DP), (/ 1, nrxxs /), 1 )
#endif
                !
                IF ( (ii+1)<=min(m,nibands(my_egrp_id+1)) ) THEN
                   ! deal with double bands
!$cuf kernel do
                   DO ig = 1, npwt
                      psi_rhoc_work_d( dfftt__nl(ig) )  =       psi_d(ig, ii) + (0._DP,1._DP) * psi_d(ig, ii+1)
                      psi_rhoc_work_d( dfftt__nlm(ig) ) = conjg(psi_d(ig, ii) - (0._DP,1._DP) * psi_d(ig, ii+1))
                   ENDDO

                ENDIF
                !
                IF ( ii==min(m,nibands(my_egrp_id+1)) ) THEN
                   ! deal with a single last band
!$cuf kernel do
                   DO ig = 1, npwt
                      psi_rhoc_work_d( dfftt__nl(ig) )  =       psi_d(ig,ii)
                      psi_rhoc_work_d( dfftt__nlm(ig) ) = conjg(psi_d(ig,ii))
                   ENDDO
                   !
                ENDIF
                !
                CALL invfft ('Wave', psi_rhoc_work_d, dfftt)
!$cuf kernel do
                DO ir = 1, nrxxs
                   temppsic_dble_d(ir)  = dble ( psi_rhoc_work_d(ir) )
                   temppsic_aimag_d(ir) = aimag( psi_rhoc_work_d(ir) )
                ENDDO
                !
             ENDIF
             !
             !
             !determine which j-bands to calculate
             jstart = 0
             jend = 0
             DO ipair=1, max_pairs
                IF(egrp_pairs(1,ipair,my_egrp_id+1).eq.ibnd)THEN
                   IF(jstart.eq.0)THEN
                      jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                      jend = jstart
                   ELSE
                      jend = egrp_pairs(2,ipair,my_egrp_id+1)
                   END IF
                END IF
             END DO
             !
             jstart = max(jstart,jblock_start)
             jend = min(jend,jblock_end)
             !
             IF(mod(jstart,2)==0) THEN
                ibnd_loop_start=jstart-1
             ELSE
                ibnd_loop_start=jstart
             ENDIF
             !
             IBND_LOOP_GAM : &
             DO jbnd=ibnd_loop_start,jend, 2 !for each band of psi
                !
                exxbuff_index = (jbnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2
                !
                IF( jbnd < jstart ) THEN
                   x1 = 0.0_DP
                ELSE
                   x1 = x_occupation(jbnd,  ik)
                ENDIF
                IF( jbnd == jend) THEN
                   x2 = 0.0_DP
                ELSE
                   x2 = x_occupation(jbnd+1,  ik)
                ENDIF
                IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
                !
                ! calculate rho in real space. Gamma tricks are used.
                ! temppsic is real; tempphic contains one band in the real part,
                ! another one in the imaginary part; the same applies to rhoc
                !
                IF( mod(ii,2) == 0 ) THEN
!$cuf kernel do
                   DO ir = 1, nrxxs
                      psi_rhoc_work_d(ir) = exxbuff_d(ir,exxbuff_index,ikq) * temppsic_aimag_d(ir) / omega
                   ENDDO

                ELSE
!$cuf kernel do
                   DO ir = 1, nrxxs
                      psi_rhoc_work_d(ir) = exxbuff_d(ir,exxbuff_index,ikq) * temppsic_dble_d(ir) / omega
                   ENDDO

                ENDIF
                !
                ! bring rho to G-space
                !
                !   >>>> add augmentation in REAL SPACE here
                IF(okvan .and. tqr) THEN
                   psi_rhoc_work = psi_rhoc_work_d
                   IF(jbnd>=jstart) &
                        CALL addusxx_r(psi_rhoc_work, &
                       _CX(becxx(ikq)%r(:,jbnd)), _CX(becpsi%r(:,ibnd)))
                   IF(jbnd<jend) &
                        CALL addusxx_r(psi_rhoc_work, &
                       _CY(becxx(ikq)%r(:,jbnd+1)),_CX(becpsi%r(:,ibnd)))
                   psi_rhoc_work_d = psi_rhoc_work
                ENDIF
                !
                CALL fwfft ('Rho', psi_rhoc_work_d, dfftt)
                !   >>>> add augmentation in G SPACE here
                IF(okvan .and. .not. tqr) THEN
                   psi_rhoc_work = psi_rhoc_work_d
                   ! contribution from one band added to real (in real space) part of rhoc
                   IF(jbnd>=jstart) &
                        CALL addusxx_g(dfftt, psi_rhoc_work, xkq,  xkp, 'r', &
                        becphi_r=becxx(ikq)%r(:,jbnd), becpsi_r=becpsi%r(:,ibnd) )
                   ! contribution from following band added to imaginary (in real space) part of rhoc
                   IF(jbnd<jend) &
                        CALL addusxx_g(dfftt, psi_rhoc_work, xkq,  xkp, 'i', &
                        becphi_r=becxx(ikq)%r(:,jbnd+1), becpsi_r=becpsi%r(:,ibnd) )
                   psi_rhoc_work_d = psi_rhoc_work 
                ENDIF
                !   >>>> charge density done
                !
                vc_d = 0._DP
                !
!$cuf kernel do
                DO ig = 1, dfftt%ngm
                   !
                   vc_d(dfftt__nl(ig))  = fac_d(ig) * psi_rhoc_work_d(dfftt__nl(ig))
                   vc_d(dfftt__nlm(ig)) = fac_d(ig) * psi_rhoc_work_d(dfftt__nlm(ig))
                   !
                ENDDO
                !
                !   >>>>  compute <psi|H_fock G SPACE here
                IF(okvan .and. .not. tqr) THEN
                   vc = vc_d
                   IF(jbnd>=jstart) &
                        CALL newdxx_g(dfftt, vc, xkq, xkp, 'r', deexx(:,ii), &
                           becphi_r=x1*becxx(ikq)%r(:,jbnd))
                   IF(jbnd<jend) &
                        CALL newdxx_g(dfftt, vc, xkq, xkp, 'i', deexx(:,ii), &
                            becphi_r=x2*becxx(ikq)%r(:,jbnd+1))
                ENDIF
                !
                !brings back v in real space
                CALL invfft ('Rho', vc_d, dfftt)
                !
                !   >>>>  compute <psi|H_fock REAL SPACE here
                IF(okvan .and. tqr) THEN
                   vc = vc_d
                   IF(jbnd>=jstart) &
                        CALL newdxx_r(dfftt,vc, _CX(x1*becxx(ikq)%r(:,jbnd)), deexx(:,ii))
                   IF(jbnd<jend) &
                        CALL newdxx_r(dfftt,vc, _CY(x2*becxx(ikq)%r(:,jbnd+1)), deexx(:,ii))
                ENDIF
                !
                IF(okpaw) THEN
                   IF(jbnd>=jstart) &
                        CALL PAW_newdxx(x1/nqs, _CX(becxx(ikq)%r(:,jbnd)),&
                                                _CX(becpsi%r(:,ibnd)), deexx(:,ii))
                   IF(jbnd<jend) &
                        CALL PAW_newdxx(x2/nqs, _CX(becxx(ikq)%r(:,jbnd+1)),&
                                                _CX(becpsi%r(:,ibnd)), deexx(:,ii))
                ENDIF
                !
                ! accumulates over bands and k points
                !
!$cuf kernel do
                DO ir = 1, nrxxs
                   result_d(ir,ii) = result_d(ir,ii) &
                                 + x1* dble(vc_d(ir))* dble(exxbuff_d(ir,exxbuff_index,ikq)) &
                                 + x2*aimag(vc_d(ir))*aimag(exxbuff_d(ir,exxbuff_index,ikq))
                ENDDO
                !
             ENDDO &
             IBND_LOOP_GAM
             !
          ENDDO &
          LOOP_ON_PSI_BANDS
          !
          ! get the next nbnd/negrp data
          IF (negrp>1) call mp_circular_shift_left( exxbuff_d(:,:,ikq), me_egrp, inter_egrp_comm )
          !
       ENDDO ! iegrp
       IF ( okvan .and..not.tqr ) CALL qvan_clean ()
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !
    DO ii=1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
       !
       IF(okvan) THEN
          CALL mp_sum(deexx(:,ii),intra_egrp_comm)
       ENDIF
       !
       !
       ! brings back result in G-space
       !
       CALL fwfft( 'Wave' , result_d(:,ii), dfftt )
       !communicate result
       !$cuf kernel do
       DO ig = 1, n
          big_result_d(ig,ibnd) = big_result_d(ig,ibnd) - exxalfa*result_d(dfftt__nl(igk_exx_d(ig,current_k)),ii)
       END DO
       big_result(:,ibnd) = big_result_d(:,ibnd)
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) CALL add_nlxx_pot (lda, big_result(:,ibnd), xkp, n, &
            igk_exx(1,current_k), deexx(:,ii), eps_occ, exxalfa)
    END DO
    !
    CALL result_sum(n*npol, m, big_result)
    IF (iexx_istart(my_egrp_id+1).gt.0) THEN
       IF (negrp == 1) then
          ending_im = m
       ELSE
          ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
       END IF
       DO im=1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
           DO ig = 1, n
              hpsi(ig,im)=hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
           ENDDO
!$omp end parallel do
       END DO
    END IF
    !
    DEALLOCATE(big_result)
    DEALLOCATE(big_result_d)
    DEALLOCATE(result_d)
    DEALLOCATE(temppsic_dble_d)
    DEALLOCATE(temppsic_aimag_d)
    DEALLOCATE(psi_rhoc_work_d)
    DEALLOCATE(psi_d)
    DEALLOCATE(vc)
    DEALLOCATE(vc_d)
    DEALLOCATE(fac_d)
    IF(okvan) DEALLOCATE(deexx)
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_bp_gamma_gpu
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_bp_k( lda, n, m, psi, hpsi, becpsi )
    !-----------------------------------------------------------------------
    !! Generic, k-point version of vexx.
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, current_k, nbnd
    USE klist,          ONLY : xk, nks, nkstot
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,         ONLY : inter_egrp_comm, my_egrp_id, negrp, &
                               intra_egrp_comm, me_egrp, &
                               max_pairs, egrp_pairs, ibands, nibands, &
                               max_ibands, iexx_istart, iexx_iend, &
                               all_start, all_end, iexx_start, jblock
    USE mp,             ONLY : mp_sum, mp_barrier, mp_circular_shift_left
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : nqs, xkq_collect, index_xkq, index_xk
    USE exx_bp_utils,   ONLY : result_sum, igk_exx
    USE io_global,      ONLY : stdout
    !
    !
    IMPLICIT NONE
    !
    INTEGER :: lda
    !! input: leading dimension of arrays psi and hpsi
    INTEGER :: n
    !! input: true dimension of psi and hpsi
    INTEGER :: m
    !! input: number of states psi
    COMPLEX(DP) :: psi(lda*npol,max_ibands)
    !! input: m wavefunctions
    COMPLEX(DP) :: hpsi(lda*npol,max_ibands)
    !! output: V_x*psi
    TYPE(bec_type), OPTIONAL :: becpsi  ! or call a calbec(...psi) instead
    !! input: <beta|psi>, optional but needed for US and PAW case
    !
    ! ... local variables
    !
    COMPLEX(DP),ALLOCATABLE :: temppsic(:,:), result(:,:)
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: result
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) result
#endif
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:,:),result_nc(:,:,:)
    INTEGER :: request_send, request_recv
    !
    COMPLEX(DP),ALLOCATABLE :: deexx(:,:)
    COMPLEX(DP),ALLOCATABLE,TARGET :: rhoc(:,:), vc(:,:)
#if defined(__USE_MANY_FFT)
    COMPLEX(DP),POINTER :: prhoc(:), pvc(:)
#endif
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: rhoc, vc
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) rhoc, vc
#endif
    REAL(DP),   ALLOCATABLE :: fac(:), facb(:)
    INTEGER :: ibnd, ik, im , ikq, iq, ipol
    INTEGER :: ir, ig, ir_start, ir_end
    INTEGER :: irt, nrt, nblock
    INTEGER :: current_ik
    INTEGER :: ibnd_loop_start
    INTEGER :: nrxxs
    REAL(DP) :: x1, x2, xkp(3), omega_inv, nqs_inv
    REAL(DP) :: xkq(3)
    INTEGER, EXTERNAL :: global_kpoint_index
    DOUBLE PRECISION :: max, tempx
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: ir_out, ipair, jbnd
    INTEGER :: ii, jstart, jend, jcount, jind
    INTEGER :: ialloc, ending_im
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    !
    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac(dfftt%ngm) )
    nrxxs= dfftt%nnr
    ALLOCATE( facb(nrxxs) )
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol,ialloc), result_nc(nrxxs,npol,ialloc) )
    ELSE
       ALLOCATE( temppsic(nrxxs,ialloc), result(nrxxs,ialloc) )
    ENDIF
    !
    IF (okvan) ALLOCATE( deexx(nkb,ialloc) )
    !
    current_ik = global_kpoint_index( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ALLOCATE( big_result(n*npol,m) )
    big_result = 0.0_DP
    !
    ! allocate arrays for rhoc and vc
    ALLOCATE( rhoc(nrxxs,jblock), vc(nrxxs,jblock) )
#if defined(__USE_MANY_FFT)
    prhoc(1:nrxxs*jblock) => rhoc(:,:)
    pvc(1:nrxxs*jblock) => vc(:,:)
#endif
    !
    DO ii = 1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd==0 .OR. ibnd>m) CYCLE
       IF (okvan) deexx(:,ii) = 0._DP
       !
       IF (noncolin) THEN
          temppsic_nc(:,:,ii) = 0._DP
       ELSE
!$omp parallel do  default(shared), private(ir) firstprivate(nrxxs)
          DO ir = 1, nrxxs
             temppsic(ir,ii) = 0._DP
          ENDDO
       ENDIF
       !
       IF (noncolin) THEN
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(dfftt%nl(igk_exx(ig,current_k)),1,ii) = psi(ig,ii)
             temppsic_nc(dfftt%nl(igk_exx(ig,current_k)),2,ii) = psi(npwx+ig,ii)
          ENDDO
!$omp end parallel do
          !
          CALL invfft( 'Wave', temppsic_nc(:,1,ii), dfftt )
          CALL invfft( 'Wave', temppsic_nc(:,2,ii), dfftt )
          !
       ELSE
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic( dfftt%nl(igk_exx(ig,current_k)), ii ) = psi(ig,ii)
          ENDDO
!$omp end parallel do
          !
          CALL invfft( 'Wave', temppsic(:,ii), dfftt )
          !
       ENDIF
       !
       IF (noncolin) THEN
!$omp parallel do default(shared) firstprivate(nrxxs) private(ir)
          DO ir = 1, nrxxs
             result_nc(ir,1,ii) = 0.0_DP
             result_nc(ir,2,ii) = 0.0_DP
          ENDDO
       ELSE
!$omp parallel do default(shared) firstprivate(nrxxs) private(ir)
          DO ir = 1, nrxxs
             result(ir,ii) = 0.0_DP
          ENDDO
       ENDIF
       !
    ENDDO
    !
    !precompute these guys
    omega_inv = 1.0 / omega
    nqs_inv = 1.0 / nqs
    !
    !------------------------------------------------------------------------!
    ! Beginning of main loop
    !------------------------------------------------------------------------!
    DO iq = 1, nqs
       !
       ikq = index_xkq(current_ik,iq)
       ik  = index_xk(ikq)
       xkq = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, current_k )
       !
! JRD - below not threaded
       facb = 0D0
       DO ig = 1, dfftt%ngm
          facb(dfftt%nl(ig)) = coulomb_fac(ig,iq,current_k)
       ENDDO
       !
       IF ( okvan .AND..NOT.tqr ) CALL qvan_init( dfftt%ngm, xkq, xkp )
       !
       DO iegrp = 1, negrp
          !
          ! compute the id of group whose data is currently worked on
          wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
          njt = (all_end(wegrp)-all_start(wegrp)+jblock)/jblock
          !
          DO ijt = 1, njt
             !
             jblock_start = (ijt - 1) * jblock + all_start(wegrp)
             jblock_end = MIN(jblock_start+jblock-1,all_end(wegrp))
             !
             DO ii = 1, nibands(my_egrp_id+1)
                !
                ibnd = ibands(ii,my_egrp_id+1)
                !
                IF (ibnd==0 .OR. ibnd>m) CYCLE
                !
                !determine which j-bands to calculate
                jstart = 0
                jend = 0
                !
                DO ipair = 1, max_pairs
                   IF (egrp_pairs(1,ipair,my_egrp_id+1) == ibnd) THEN
                      IF (jstart == 0) THEN
                         jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                         jend = jstart
                      ELSE
                         jend = egrp_pairs(2,ipair,my_egrp_id+1)
                      ENDIF
                   ENDIF
                ENDDO
                !
                jstart = MAX( jstart, jblock_start )
                jend = MIN( jend, jblock_end )
                !
                !how many iters
                jcount = jend-jstart+1
                IF (jcount <= 0) CYCLE
                !
                !----------------------------------------------------------------------!
                !INNER LOOP START
                !----------------------------------------------------------------------!
                !
                nblock = 2048
                nrt = (nrxxs+nblock-1)/nblock
                !
!$omp parallel do collapse(2) private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd = jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = MIN(ir_start+nblock-1, nrxxs)
                      IF (noncolin) THEN
                         DO ir = ir_start, ir_end
                            rhoc(ir,jbnd-jstart+1) = ( CONJG(exxbuff(ir,jbnd-all_start(wegrp)+ &
                                                       iexx_start,ikq))*temppsic_nc(ir,1,ii) + &
                                                       CONJG(exxbuff(nrxxs+ir,jbnd-all_start(wegrp)+ &
                                                       iexx_start,ikq))*temppsic_nc(ir,2,ii) )/omega
                         ENDDO
                      ELSE
!DIR$ vector nontemporal (rhoc)
                         DO ir = ir_start, ir_end
                            rhoc(ir,jbnd-jstart+1) = CONJG(exxbuff(ir,jbnd-all_start(wegrp)+ &
                                                     iexx_start,ikq))*temppsic(ir,ii)*omega_inv
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
!$omp end parallel do
                !
                !   ... add augmentation in REAL space HERE
                IF (okvan .AND. tqr) THEN ! augment the "charge" in real space
                   DO jbnd = jstart, jend
                      CALL addusxx_r(rhoc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd), becpsi%k(:,ibnd))
                   ENDDO
                ENDIF
                !
                !   ... brings it to G-space
#if defined(__USE_MANY_FFT)
                CALL fwfft( 'Rho', prhoc, dfftt, howmany=jcount )
#else
                DO jbnd=jstart, jend
                   CALL fwfft( 'Rho', rhoc(:,jbnd-jstart+1), dfftt )
                ENDDO
#endif
                !
                !   ... add augmentation in G space HERE
                IF (okvan .AND. .NOT. tqr) THEN
                   DO jbnd = jstart, jend
                      CALL addusxx_g( dfftt, rhoc(:,jbnd-jstart+1), xkq, xkp, &
                      'c', becphi_c=becxx(ikq)%k(:,jbnd),becpsi_c=becpsi%k(:,ibnd) )
                   ENDDO
                ENDIF
                !   ... charge done
                !
!call start_collection()
!$omp parallel do collapse(2) private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd = jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = MIN(ir_start+nblock-1,nrxxs)
!DIR$ vector nontemporal (vc)
                      DO ir = ir_start, ir_end
                         vc(ir,jbnd-jstart+1) = facb(ir) * rhoc(ir,jbnd-jstart+1)*&
                                                x_occupation(jbnd,ik) * nqs_inv
                      ENDDO
                   ENDDO
                ENDDO
!$omp end parallel do
!call stop_collection()
                !
                ! Add ultrasoft contribution (RECIPROCAL SPACE)
                ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
                IF (okvan .AND. .NOT. tqr) THEN
                   DO jbnd=jstart, jend
                      CALL newdxx_g( dfftt, vc(:,jbnd-jstart+1), xkq, xkp, 'c',&
                                     deexx(:,ii), becphi_c=becxx(ikq)%k(:,jbnd) )
                   ENDDO
                ENDIF
                !
                !brings back v in real space
#if defined(__USE_MANY_FFT)
                !fft many
                CALL invfft( 'Rho', pvc, dfftt, howmany=jcount )
#else
                DO jbnd = jstart, jend
                   CALL invfft( 'Rho', vc(:,jbnd-jstart+1), dfftt )
                ENDDO
#endif
                !
                ! Add ultrasoft contribution (REAL SPACE)
                IF (okvan .AND. tqr) THEN
                   DO jbnd = jstart, jend
                      CALL newdxx_r( dfftt, vc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd),deexx(:,ii) )
                   ENDDO
                ENDIF
                !
                ! ... Add PAW one-center contribution
                IF (okpaw) THEN
                   DO jbnd = jstart, jend
                      CALL PAW_newdxx( x_occupation(jbnd,ik)/nqs, becxx(ikq)%k(:,jbnd), &
                                       becpsi%k(:,ibnd), deexx(:,ii) )
                   ENDDO
                ENDIF
                !
                ! ... accumulates over bands and k points
                !
!call start_collection()
!$omp parallel do private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd = jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = MIN(ir_start+nblock-1, nrxxs)
                      IF (noncolin) THEN
                         DO ir = ir_start, ir_end
                            result_nc(ir,1,ii) = result_nc(ir,1,ii) + vc(ir,jbnd-jstart+1) * &
                                                 exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq)
                            result_nc(ir,2,ii) = result_nc(ir,2,ii) + vc(ir,jbnd-jstart+1) * &
                                                 exxbuff(ir+nrxxs,jbnd-all_start(wegrp)+iexx_start,ikq)
                         ENDDO
                      ELSE
!!dir$ vector nontemporal (result)
                         DO ir = ir_start, ir_end
                            result(ir,ii) = result(ir,ii) + vc(ir,jbnd-jstart+1)* &
                                            exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq)
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
!$omp end parallel do
!call stop_collection()
                !
                !----------------------------------------------------------------------!
                !INNER LOOP END
                !----------------------------------------------------------------------!
                !
             ENDDO !I-LOOP
          ENDDO !IJT
          !
          ! get the next nbnd/negrp data
          IF (negrp > 1) CALL mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
          !
       ENDDO !iegrp
       !
       IF ( okvan .AND..NOT.tqr ) CALL qvan_clean()
       !
    ENDDO
    !
    !
    !
    DO ii = 1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd==0 .OR. ibnd>m) CYCLE
       !
       IF (okvan) THEN
          CALL mp_sum( deexx(:,ii), intra_egrp_comm )
       ENDIF
       !
       IF (noncolin) THEN
          !brings back result in G-space
          CALL fwfft( 'Wave', result_nc(:,1,ii), dfftt )
          CALL fwfft( 'Wave', result_nc(:,2,ii), dfftt )
          !
          DO ig = 1, n
             big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa* &
                                   result_nc(dfftt%nl(igk_exx(ig,current_k)),1,ii)
             big_result(n+ig,ibnd) = big_result(n+ig,ibnd) - exxalfa* &
                                     result_nc(dfftt%nl(igk_exx(ig,current_k)),2,ii)
          ENDDO
       ELSE
          !
          CALL fwfft( 'Wave', result(:,ii), dfftt )
          !
          DO ig = 1, n
             big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa* &
                                   result(dfftt%nl(igk_exx(ig,current_k)),ii)
          ENDDO
       ENDIF
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF (okvan) CALL add_nlxx_pot( lda, big_result(:,ibnd), xkp, n, igk_exx(:,current_k), &
                                    deexx(:,ii), eps_occ, exxalfa )
       !
    ENDDO
    !
    !deallocate temporary arrays
    DEALLOCATE( rhoc, vc )
    !
    !sum result
    CALL result_sum( n*npol, m, big_result )
    !
    IF (iexx_istart(my_egrp_id+1) > 0) THEN
       !
       IF (negrp == 1) THEN
          ending_im = m
       ELSE
          ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
       ENDIF
       !
       IF (noncolin) THEN
          DO im = 1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(ig,im) = hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(lda+ig,im) = hpsi(lda+ig,im) + big_result(n+ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
          ENDDO
       ELSE
          DO im = 1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(ig,im) = hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
          ENDDO
       ENDIF
    ENDIF
    !
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc, result_nc )
    ELSE
       DEALLOCATE( temppsic, result )
    ENDIF
    !
    DEALLOCATE( big_result )
    DEALLOCATE( fac, facb )
    IF (okvan) DEALLOCATE( deexx )
    !
  END SUBROUTINE vexx_bp_k
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_bp_k_gpu(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... generic, k-point version of vexx
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, current_k, nbnd
    USE klist,          ONLY : xk, nks, nkstot
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,         ONLY : inter_egrp_comm, my_egrp_id, negrp, &
                               intra_egrp_comm, me_egrp, &
                               max_pairs, egrp_pairs, ibands, nibands, &
                               max_ibands, iexx_istart, iexx_iend, &
                               all_start, all_end, iexx_start, jblock
    USE mp,             ONLY : mp_sum, mp_barrier, mp_circular_shift_left
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : nqs, xkq_collect, index_xkq, index_xk
    USE exx_bp_utils,   ONLY : result_sum, igk_exx, igk_exx_d
    !CUDA stuff
    USE mp_exx,         ONLY : iexx_istart_d
    USE io_global,      ONLY : stdout
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,max_ibands)
    COMPLEX(DP)              :: hpsi(lda*npol,max_ibands)
#if defined(__CUDA)
    attributes(DEVICE) :: psi_d, hpsi_d
#endif
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: psi_d(:,:)
    COMPLEX(DP),ALLOCATABLE :: hpsi_d(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_d(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc_d(:,:,:)
    COMPLEX(DP),ALLOCATABLE :: result_d(:,:), result_nc_d(:,:,:)
#if defined(__CUDA)
    attributes(DEVICE) :: temppsic_d, temppsic_nc_d, result_d, result_nc_d
#endif
    INTEGER          :: request_send, request_recv
    !
    COMPLEX(DP),ALLOCATABLE :: deexx(:,:)
    COMPLEX(DP),ALLOCATABLE,TARGET :: rhoc(:,:), vc(:,:)
    COMPLEX(DP),ALLOCATABLE,TARGET :: rhoc_d(:,:), vc_d(:,:)
    COMPLEX(DP),POINTER :: prhoc_d(:), pvc_d(:)
#if defined(__CUDA)
    attributes(DEVICE) :: rhoc_d, vc_d, prhoc_d, pvc_d
#endif
    REAL(DP), ALLOCATABLE :: fac(:), facb(:)
    REAL(DP), ALLOCATABLE :: facb_d(:)
#if defined(__CUDA)
    attributes(DEVICE) :: facb_d
#endif
    INTEGER  :: ibnd, ik, im , ikq, iq, ipol
    INTEGER  :: ir, ig, ir_start, ir_end
    INTEGER  :: irt, nrt, nblock
    INTEGER  :: current_ik
    INTEGER  :: ibnd_loop_start
    INTEGER  :: nrxxs
    REAL(DP) :: x1, x2, xkp(3), omega_inv, nqs_inv
    REAL(DP) :: xkq(3)
    INTEGER, EXTERNAL :: global_kpoint_index
    DOUBLE PRECISION :: max, tempx
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    COMPLEX(DP), ALLOCATABLE :: big_result_d(:,:)
#if defined(__CUDA)
    attributes(DEVICE) :: big_result_d
#endif
    INTEGER :: ir_out, ipair, jbnd
    INTEGER :: ii, jstart, jend, jcount, jind, jcurr
    INTEGER :: ialloc, ending_im
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    INTEGER :: all_start_tmp
    !hack around PGI bug
    INTEGER, POINTER :: dfftt__nl(:)
#if defined(__CUDA)
    attributes(DEVICE) :: dfftt__nl
#endif
    !
    dfftt__nl=>dfftt%nl_d
    !
    CALL start_clock( 'vexx_k_setup' )

    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac(dfftt%ngm) )
    nrxxs= dfftt%nnr
    ALLOCATE( facb(nrxxs) )

    ALLOCATE( psi_d, source=psi )
    ALLOCATE( hpsi_d, source=hpsi )
    ALLOCATE( facb_d(nrxxs) )

    !initial copy of exxbuff
    exxbuff_d = exxbuff
    !
    IF (noncolin) THEN
       ALLOCATE( result_nc_d(nrxxs,npol,ialloc) )

       !temppsic_d knows where it is
       ALLOCATE( temppsic_nc_d(nrxxs,npol,ialloc) )
    ELSE
       ALLOCATE( result_d(nrxxs,ialloc) )

       !temppsic_d knows
       ALLOCATE( temppsic_d(nrxxs,ialloc) )
    ENDIF
    !
    IF(okvan) ALLOCATE(deexx(nkb,ialloc))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    allocate(big_result(n*npol,m))
    big_result = 0.0_DP
    allocate(big_result_d(n*npol,m))
    big_result_d = 0.0_DP
    !
    !allocate arrays for rhoc and vc
    ALLOCATE(rhoc_d(nrxxs,jblock), vc_d(nrxxs,jblock))
    ALLOCATE(rhoc(nrxxs,jblock), vc(nrxxs,jblock))
    
    !
    
    DO ii=1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
       !
       IF(okvan) deexx(:,ii) = 0._DP
       !
       IF (noncolin) THEN
          temppsic_nc_d(:,:,ii) = 0._DP
       ELSE
          temppsic_d(:,ii) = 0._DP
       END IF
       !
       IF (noncolin) THEN
          !$cuf kernel do (1)
          DO ig = 1, n
             temppsic_nc_d(dfftt__nl(igk_exx_d(ig,current_k)),1,ii) = psi_d(ig,ii)
             temppsic_nc_d(dfftt__nl(igk_exx_d(ig,current_k)),2,ii) = psi_d(npwx+ig,ii)
          ENDDO
          CALL invfft ('Wave', temppsic_nc_d(:,1,ii), dfftt)
          CALL invfft ('Wave', temppsic_nc_d(:,2,ii), dfftt)
       ELSE
          !$cuf kernel do (1)
          DO ig = 1, n
             temppsic_d( dfftt__nl(igk_exx_d(ig,current_k)), ii ) = psi_d(ig,ii)
          ENDDO
          CALL invfft ('Wave', temppsic_d(:,ii), dfftt)
       END IF
    END DO

    IF (noncolin) THEN
       result_nc_d = 0.0_DP
    ELSE
       result_d = 0.0_DP
    ENDIF

    ! no longer need psi_d
    DEALLOCATE(psi_d)

    !
    !precompute these guys
    omega_inv = 1.0 / omega
    nqs_inv = 1.0 / nqs
    !
    CALL stop_clock( 'vexx_k_setup' )
    CALL start_clock( 'vexx_k_main' )
    !------------------------------------------------------------------------!
    ! Beginning of main loop
    !------------------------------------------------------------------------!
    vexxmain: DO iq=1, nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution_all(dfftt%ngm, gt, xkp, xkq, iq, current_k)
       !
! JRD - below not threaded
       facb = 0D0
       DO ig = 1, dfftt%ngm
          facb(dfftt%nl(ig)) = coulomb_fac(ig,iq,current_k)
       ENDDO
       facb_d = facb
       !
       IF ( okvan .and..not.tqr ) CALL qvan_init (dfftt%ngm, xkq, xkp)
       !
       DO iegrp=1, negrp
          !
          ! compute the id of group whose data is currently worked on
          wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
          njt = (all_end(wegrp)-all_start(wegrp)+jblock)/jblock
          !
          DO ijt=1, njt
             !
             jblock_start = (ijt - 1) * jblock + all_start(wegrp)
             jblock_end = min(jblock_start+jblock-1,all_end(wegrp))
             !
             DO ii=1, nibands(my_egrp_id+1)
                !
                ibnd = ibands(ii,my_egrp_id+1)
                !
                IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
                !
                !determine which j-bands to calculate
                jstart = 0
                jend = 0
                DO ipair=1, max_pairs
                   IF(egrp_pairs(1,ipair,my_egrp_id+1).eq.ibnd)THEN
                      IF(jstart.eq.0)THEN
                         jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                         jend = jstart
                      ELSE
                         jend = egrp_pairs(2,ipair,my_egrp_id+1)
                      END IF
                   END IF
                END DO
                !
                jstart = max(jstart,jblock_start)
                jend = min(jend,jblock_end)
                !
                !how many iters
                jcount=jend-jstart+1
                if(jcount<=0) cycle
                !
                !----------------------------------------------------------------------!
                !INNER LOOP START
                !----------------------------------------------------------------------!
                !
                nblock=2048
                nrt = (nrxxs+nblock-1)/nblock
                !
associate(rhoc=>rhoc_d, exxbuff=>exxbuff_d)
                all_start_tmp=all_start(wegrp)
                !$cuf kernel do (2)
                DO jbnd=jstart, jend
                   DO ir = 1, nrxxs

                     IF (noncolin) THEN
                       rhoc(ir,jbnd-jstart+1) = &
                       (conjg(exxbuff(ir,jbnd-all_start_tmp+iexx_start,ikq))*temppsic_nc_d(ir,1,ii) +&
                       conjg(exxbuff(nrxxs+ir,jbnd-all_start_tmp+iexx_start,ikq))*temppsic_nc_d(ir,2,ii)) * omega_inv
                     ELSE

                       rhoc(ir,jbnd-jstart+1) = &
                       conjg(exxbuff(ir,jbnd-all_start_tmp+iexx_start,ikq))*temppsic_d(ir,ii)* omega_inv
                     ENDIF

                   ENDDO
                ENDDO
end associate
                !
                !   >>>> add augmentation in REAL space HERE
                IF(okvan .and. tqr) THEN ! augment the "charge" in real space
                   DO jbnd=jstart, jend
                      CALL addusxx_r(rhoc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd), becpsi%k(:,ibnd))
                   ENDDO
                ENDIF
                !
                !   >>>> brings it to G-space
                !
                DO jbnd=jstart, jend, many_fft
                  jcurr = min(many_fft, jend-jbnd+1)
                  prhoc_d(1:nrxxs*jcurr) => rhoc_d(:,jbnd-jstart+1:jbnd-jstart+jcurr)
                  CALL fwfft ('Rho', prhoc_d, dfftt, howmany=jcurr)
                ENDDO
                !
                !   >>>> add augmentation in G space HERE
                IF(okvan .and. .not. tqr) THEN
                   rhoc = rhoc_d
                   DO jbnd=jstart, jend
                      CALL addusxx_g(dfftt, rhoc(:,jbnd-jstart+1), xkq, xkp, &
                      'c', becphi_c=becxx(ikq)%k(:,jbnd),becpsi_c=becpsi%k(:,ibnd))
                   ENDDO
                   rhoc_d = rhoc
                ENDIF
                !   >>>> charge done
                !
associate(vc=>vc_d, facb=>facb_d, rhoc=>rhoc_d, x_occupation=>x_occupation_d)
                !$cuf kernel do (2)
                DO jbnd=jstart, jend
                   DO ir = 1, nrxxs
                         vc(ir,jbnd-jstart+1) = facb(ir) * rhoc(ir,jbnd-jstart+1)*&
                                                x_occupation(jbnd,ik) * nqs_inv
                   ENDDO
                ENDDO
end associate
                !
                ! Add ultrasoft contribution (RECIPROCAL SPACE)
                ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
                IF(okvan .and. .not. tqr) THEN
                   vc = vc_d
                   DO jbnd=jstart, jend
                      CALL newdxx_g(dfftt, vc(:,jbnd-jstart+1), xkq, xkp, 'c',&
                                    deexx(:,ii), becphi_c=becxx(ikq)%k(:,jbnd))
                   ENDDO
                   vc_d = vc
                ENDIF
                !
                !brings back v in real space
                DO jbnd=jstart, jend, many_fft
                  jcurr = min(many_fft, jend-jbnd+1)
                  pvc_d(1:nrxxs*jcurr) => vc_d(:,jbnd-jstart+1:jbnd-jstart+jcurr)
                  CALL invfft ('Rho', pvc_d, dfftt, howmany=jcurr)
                ENDDO
                !
                ! Add ultrasoft contribution (REAL SPACE)
                IF(okvan .and. tqr) THEN
                   vc = vc_d
                   DO jbnd=jstart, jend
                      CALL newdxx_r(dfftt, vc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd),deexx(:,ii))
                   ENDDO
                   vc_d = vc
                ENDIF
                !
                ! Add PAW one-center contribution
                IF(okpaw) THEN
                   vc = vc_d
                   DO jbnd=jstart, jend
                      CALL PAW_newdxx(x_occupation(jbnd,ik)/nqs, becxx(ikq)%k(:,jbnd), becpsi%k(:,ibnd), deexx(:,ii))
                   ENDDO
                   vc_d = vc
                ENDIF
                !
                !accumulates over bands and k points
                !

associate(exxbuff=>exxbuff_d, vc=>vc_d)
                all_start_tmp=all_start(wegrp)
                DO jbnd=jstart, jend
                   !$cuf kernel do (1)
                   DO ir = 1, nrxxs
                      IF (noncolin) THEN
                         result_nc_d(ir,1,ii) = result_nc_d(ir,1,ii) &
                              + vc(ir,jbnd-jstart+1) * exxbuff(ir,jbnd-all_start_tmp+iexx_start,ikq)
                         result_nc_d(ir,2,ii) = result_nc_d(ir,2,ii) &
                              + vc(ir,jbnd-jstart+1) * exxbuff(ir+nrxxs,jbnd-all_start_tmp+iexx_start,ikq)
                      ELSE
                         result_d(ir,ii) = result_d(ir,ii) &
                              + vc(ir,jbnd-jstart+1)*exxbuff(ir,jbnd-all_start_tmp+iexx_start,ikq)
                      ENDIF
                   ENDDO
                ENDDO
end associate
                !
                !----------------------------------------------------------------------!
                !INNER LOOP END
                !----------------------------------------------------------------------!
                !
             END DO !I-LOOP
          END DO !IJT
          !
          ! get the next nbnd/negrp data
          IF (negrp>1) THEN
             call mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
             exxbuff_d = exxbuff
          ENDIF
          !
       END DO !iegrp
       !
       IF ( okvan .and..not.tqr ) CALL qvan_clean ()
    END DO vexxmain

!move this down to after the vexx_k_fin

    CALL stop_clock( 'vexx_k_main' )
    CALL start_clock( 'vexx_k_fin' )
    !
    !
    !
    DO ii=1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
       !
       IF(okvan) THEN
          CALL mp_sum(deexx(:,ii),intra_egrp_comm)
       ENDIF
       !
       !big_result_d=big_result !already initialized along with the big_result=1.0D0
       IF (noncolin) THEN
          !brings back result in G-space
          CALL fwfft ('Wave', result_nc_d(:,1,ii), dfftt)
          CALL fwfft ('Wave', result_nc_d(:,2,ii), dfftt)
          !$cuf kernel do (1)
          DO ig = 1, n
             big_result_d(ig,ibnd) = big_result_d(ig,ibnd) - exxalfa*result_nc_d(dfftt__nl(igk_exx_d(ig,current_k)),1,ii)
             big_result_d(n+ig,ibnd) = big_result_d(n+ig,ibnd) - exxalfa*result_nc_d(dfftt__nl(igk_exx_d(ig,current_k)),2,ii)
          ENDDO
       ELSE
          !
          CALL fwfft ('Wave', result_d(:,ii), dfftt)
          !$cuf kernel do (1)
          DO ig = 1, n
             big_result_d(ig,ibnd) = big_result_d(ig,ibnd) - exxalfa*result_d(dfftt__nl(igk_exx_d(ig,current_k)),ii)
          ENDDO
       ENDIF
       big_result(:,ibnd) = big_result_d(:,ibnd)

       IF(okvan) CALL add_nlxx_pot (lda, big_result(:,ibnd), xkp, n, igk_exx(:,current_k),&
            deexx(:,ii), eps_occ, exxalfa)
       !
    END DO

    ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
    !deallocate temporary arrays
    DEALLOCATE(rhoc, vc)
    IF (noncolin) THEN
       DEALLOCATE( result_nc_d )
    ELSE
       DEALLOCATE( result_d )
    ENDIF

    !dealloc stuff
    DEALLOCATE(rhoc_d, vc_d)
    !
    !sum result
    CALL result_sum(n*npol, m, big_result)
    big_result_d = big_result
    IF (iexx_istart(my_egrp_id+1).gt.0) THEN
       IF (negrp == 1) then
          ending_im = m
       ELSE
          ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
       END IF

       !iexx_istart_d=iexx_istart
       IF(noncolin) THEN
          !$cuf kernel do (2)
          DO im=1, ending_im
             DO ig = 1, n
                hpsi_d(ig,im) = hpsi_d(ig,im) + big_result_d(ig,im+iexx_istart_d(my_egrp_id+1)-1)
                hpsi_d(lda+ig,im) = hpsi_d(lda+ig,im) + big_result_d(n+ig,im+iexx_istart_d(my_egrp_id+1)-1)
             ENDDO
          END DO
       ELSE
          !$cuf kernel do (2)
          DO im=1, ending_im
             DO ig = 1, n
                hpsi_d(ig,im) = hpsi_d(ig,im) + big_result_d(ig,im+iexx_istart_d(my_egrp_id+1)-1)
             ENDDO
          ENDDO
       END IF
    END IF
    hpsi=hpsi_d

    !these need to be deallocated anyhow
    DEALLOCATE(big_result)

    DEALLOCATE(fac, facb )

    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc_d)
    ELSE
       DEALLOCATE(temppsic_d)
    ENDIF

    IF(okvan) DEALLOCATE( deexx)

    DEALLOCATE(big_result_d)
    DEALLOCATE(facb_d)
    DEALLOCATE(hpsi_d)
    !
    CALL stop_clock( 'vexx_k_fin' )
    !
    !------------------------------------------------------------------------
  END SUBROUTINE vexx_bp_k_gpu
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy_bp_gamma()
    !-----------------------------------------------------------------------
    !
    USE constants,               ONLY : fpi, e2, pi
    USE io_files,                ONLY : iunwfc_exx, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g
    USE wvfct,                   ONLY : nbnd, npwx, wg
    USE wavefunctions,           ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_bands,                ONLY : intra_bgrp_comm
    USE mp_exx,                  ONLY : inter_egrp_comm, my_egrp_id, negrp, &
                                        intra_egrp_comm, me_egrp, &
                                        max_pairs, egrp_pairs, ibands, nibands, &
                                        iexx_istart, iexx_iend, &
                                        all_start, all_end, iexx_start, &
                                        init_index_over_band, jblock
    USE mp,                      ONLY : mp_sum, mp_circular_shift_left
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, &
                                        deallocate_bec_type, calbec
    USE paw_variables,           ONLY : okpaw
    USE paw_exx,                 ONLY : PAW_xx_energy
    USE us_exx,                  ONLY : bexg_merge, becxx, addusxx_g, &
                                        addusxx_r, qvan_init, qvan_clean
    USE exx_base,                ONLY : nqs, xkq_collect, index_xkq, index_xk
    USE exx_bp_utils,            ONLY : igk_exx, change_data_structure, &
                                        transform_evc_to_exx, nwordwfc_exx, &
                                        evc_exx
    USE uspp_init,               ONLY : init_us_2
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy_bp_gamma
    !
    ! ... local variables
    !
    REAL(DP) :: energy
    COMPLEX(DP), ALLOCATABLE :: temppsic(:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:)
    REAL(DP),    ALLOCATABLE :: fac(:)
    COMPLEX(DP), ALLOCATABLE :: vkb_exx(:,:)
    INTEGER  :: jbnd, ibnd, ik, ikk, ig, ikq, iq, ir
    INTEGER  :: nrxxs, current_ik, ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    INTEGER :: npw
    INTEGER :: istart, iend, ipair, ii, ialloc
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: exxbuff_index
    INTEGER :: calbec_start, calbec_end
    INTEGER :: intra_bgrp_comm_
    INTEGER :: iegrp, wegrp
    INTEGER :: ibnd_start
    !
    CALL init_index_over_band( inter_egrp_comm, nbnd, nbnd )
    !
    CALL transform_evc_to_exx( 0 )
    !
    ialloc = nibands(my_egrp_id+1)
    !
    nrxxs = dfftt%nnr
    ALLOCATE( fac(dfftt%ngm) )
    !
    ALLOCATE( temppsic(nrxxs), temppsic_DBLE(nrxxs), temppsic_aimag(nrxxs) )
    ALLOCATE( rhoc(nrxxs) )
    ALLOCATE( vkb_exx(npwx,nkb) )
    !
    energy = 0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi )
    !
    IKK_LOOP : &
    DO ikk = 1, nks
       current_ik = global_kpoint_index( nkstot, ikk )
       xkp = xk(:,ikk)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ikk )
       !
       ! ... prepare the |beta> function at k+q
       CALL init_us_2( npw, igk_exx(:,ikk), xkp, vkb_exx )
       !
       ! ... compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       calbec_start = ibands(1,my_egrp_id+1)
       calbec_end = ibands(nibands(my_egrp_id+1),my_egrp_id+1)
       !
       intra_bgrp_comm_ = intra_bgrp_comm
       intra_bgrp_comm = intra_egrp_comm
       !
       CALL calbec( npw, vkb_exx, evc_exx, becpsi, nibands(my_egrp_id+1) )
       !
       intra_bgrp_comm = intra_bgrp_comm_
       !
       IQ_LOOP : &
       DO iq = 1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          !
          xkq = xkq_collect(:,ikq)
          !
          CALL g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, current_ik )
          !
          fac = coulomb_fac(:,iq,current_ik)
          fac(gstart_t:) = 2 * coulomb_fac(gstart_t:,iq,current_ik)
          !
          IF ( okvan .AND. .NOT.tqr ) CALL qvan_init( dfftt%ngm, xkq, xkp )
          !
          DO iegrp = 1, negrp
             !
             ! ... compute the id of group whose data is currently worked on
             wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
             !
             jblock_start = all_start(wegrp)
             jblock_end   = all_end(wegrp)
             !
             JBND_LOOP : &
             DO ii = 1, nibands(my_egrp_id+1)
                !
                jbnd = ibands(ii,my_egrp_id+1)
                !
                IF (jbnd==0 .OR. jbnd>nbnd) CYCLE
                !
                IF ( MOD(ii,2)==1 ) THEN
                   !
                   temppsic = (0._DP,0._DP)
                   !
                   IF ( (ii+1) <= nibands(my_egrp_id+1) ) THEN
                      ! deal with double bands
!$omp parallel do  default(shared), private(ig)
                      DO ig = 1, npwt
                         temppsic( dfftt%nl(ig) )  = &
                              evc_exx(ig,ii) + (0._DP,1._DP) * evc_exx(ig,ii+1)
                         temppsic( dfftt%nlm(ig) ) = &
                              CONJG(evc_exx(ig,ii) - (0._DP,1._DP) * evc_exx(ig,ii+1))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   IF (ii == nibands(my_egrp_id+1)) THEN
                      ! deal with a single last band
!$omp parallel do  default(shared), private(ig)
                      DO ig = 1, npwt
                         temppsic( dfftt%nl(ig) ) = evc_exx(ig,ii)
                         temppsic( dfftt%nlm(ig) ) = CONJG(evc_exx(ig,ii))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   CALL invfft( 'Wave', temppsic, dfftt )
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      temppsic_DBLE(ir) = DBLE( temppsic(ir) )
                      temppsic_aimag(ir) = AIMAG( temppsic(ir) )
                   ENDDO
!$omp end parallel do
                   !
                ENDIF
                !
                !determine which j-bands to calculate
                istart = 0
                iend = 0
                !
                DO ipair = 1, max_pairs
                   IF (egrp_pairs(1,ipair,my_egrp_id+1) == jbnd) THEN
                      IF (istart == 0) THEN
                         istart = egrp_pairs(2,ipair,my_egrp_id+1)
                         iend = istart
                      ELSE
                         iend = egrp_pairs(2,ipair,my_egrp_id+1)
                      ENDIF
                   ENDIF
                ENDDO
                !
                istart = MAX(istart,jblock_start)
                iend = MIN(iend,jblock_end)
                !
                IF (MOD(istart,2) == 0) THEN
                   ibnd_loop_start = istart-1
                ELSE
                   ibnd_loop_start = istart
                ENDIF
                !
                IBND_LOOP_GAM : &
                DO ibnd = ibnd_loop_start, iend, 2       !for each band of psi
                   !
                   exxbuff_index = (ibnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2
                   !
                   IF ( ibnd < istart ) THEN
                      x1 = 0.0_DP
                   ELSE
                      x1 = x_occupation(ibnd,ik)
                   ENDIF
                   !
                   IF ( ibnd < iend ) THEN
                      x2 = x_occupation(ibnd+1,ik)
                   ELSE
                      x2 = 0.0_DP
                   ENDIF
                   ! calculate rho in real space. Gamma tricks are used.
                   ! temppsic is real; tempphic contains band 1 in the real part,
                   ! band 2 in the imaginary part; the same applies to rhoc
                   !
                   IF ( MOD(ii,2) == 0 ) THEN
                      rhoc = 0.0_DP
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         rhoc(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_aimag(ir) / omega
                      ENDDO
!$omp end parallel do
                   ELSE
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         rhoc(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_DBLE(ir) / omega
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   IF (okvan .AND. tqr) THEN
                      IF (ibnd >= istart) &
                           CALL addusxx_r( rhoc,_CX(becxx(ikq)%r(:,ibnd)), &
                                          _CX(becpsi%r(:,jbnd)) )
                      IF (ibnd<iend) &
                           CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)), &
                                          _CX(becpsi%r(:,jbnd)))
                   ENDIF
                   !
                   ! bring rhoc to G-space
                   CALL fwfft( 'Rho', rhoc, dfftt )
                   !
                   IF (okvan .AND. .NOT.tqr) THEN
                      IF (ibnd >= istart) &
                           CALL addusxx_g( dfftt, rhoc, xkq, xkp, 'r', &
                           becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,jbnd-calbec_start+1) )
                      IF (ibnd < iend) &
                           CALL addusxx_g( dfftt, rhoc, xkq, xkp, 'i', &
                           becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,jbnd-calbec_start+1) )
                   ENDIF
                   !
                   vc = 0.0_DP
!$omp parallel do  default(shared), private(ig),  reduction(+:vc)
                   DO ig = 1, dfftt%ngm
                      !
                      ! The real part of rhoc contains the contribution from band ibnd
                      ! The imaginary part    contains the contribution from band ibnd+1
                      !
                      vc = vc + fac(ig) * ( x1 * &
                           ABS( rhoc(dfftt%nl(ig)) + CONJG(rhoc(dfftt%nlm(ig))) )**2 &
                                 +x2 * &
                           ABS( rhoc(dfftt%nl(ig)) - CONJG(rhoc(dfftt%nlm(ig))) )**2 )
                   ENDDO
!$omp end parallel do
                   !
                   vc = vc * omega * 0.25_DP / nqs
                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                   !
                   IF (okpaw) THEN
                      IF (ibnd >= ibnd_start) &
                           energy = energy + exxalfa*wg(jbnd,ikk)*&
                           x1 * PAW_xx_energy(_CX(becxx(ikq)%r(:,ibnd)),_CX(becpsi%r(:,jbnd)) )
                      IF (ibnd < ibnd_end) &
                           energy = energy + exxalfa*wg(jbnd,ikk)*&
                           x2 * PAW_xx_energy(_CX(becxx(ikq)%r(:,ibnd+1)), _CX(becpsi%r(:,jbnd)) )
                   ENDIF
                   !
                ENDDO &
                IBND_LOOP_GAM
             ENDDO &
             JBND_LOOP
             !
             ! get the next nbnd/negrp data
             IF (negrp > 1) CALL mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
             !
          ENDDO ! iegrp
          IF ( okvan .AND. .NOT.tqr ) CALL qvan_clean( )
          !
       ENDDO &
       IQ_LOOP
    ENDDO &
    IKK_LOOP
    !
    DEALLOCATE( temppsic, temppsic_dble, temppsic_aimag )
    !
    DEALLOCATE( rhoc, fac )
    CALL deallocate_bec_type( becpsi )
    !
    CALL mp_sum( energy, inter_egrp_comm )
    CALL mp_sum( energy, intra_egrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy_bp_gamma = energy
    !
    CALL change_data_structure( .FALSE. )
    !
  END FUNCTION  exxenergy_bp_gamma
  !
  !----------------------------------------------------------------------
  SUBROUTINE compute_becpsi( npw_, igk_, q_, evc_exx, becpsi_k )
  !----------------------------------------------------------------------
  !! Calculates becpsi_k = <vkb|evc_exx> - FIXME: untested
  !
  USE kinds,         ONLY : DP
  USE wvfct,         ONLY : npwx, nbnd
  USE uspp,          ONLY : nkb
  USE uspp_param,    ONLY : lmaxkb
  USE becmod,        ONLY : calbec
  USE mp_exx,        ONLY : ibands, nibands, my_egrp_id
  USE uspp_init,     ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw_
  !! number of PWs
  INTEGER, INTENT(IN) :: igk_(npw_)
  !! indices of G in the list of q+G vectors
  REAL(DP), INTENT(IN) :: q_(3)
  !! q vector (2pi/a units)
  COMPLEX(DP), INTENT(IN) :: evc_exx(npwx,nibands(my_egrp_id+1))
  !! wavefunctions from the PW set to exx
  COMPLEX(DP), INTENT(OUT) :: becpsi_k(nkb,nibands(my_egrp_id+1))
  !! <beta|psi> for k points
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE :: vkb_(:,:) !beta functions (npw_ <= npwx)
  INTEGER :: istart, iend
  !
  IF (lmaxkb < 0) RETURN
  !
  istart = ibands(1,my_egrp_id+1)
  iend = ibands(nibands(my_egrp_id+1),my_egrp_id+1)
  !
  write(6,*) 'WARNING: compute_becpsi UNTESTED'
  ALLOCATE( vkb_(npwx,nkb) )
  !
  CALL init_us_2( npw_, igk_, q_, vkb_ )
  !
  CALL calbec( npw_, vkb_, evc_exx, becpsi_k, nibands(my_egrp_id+1) )
  !
  DEALLOCATE( vkb_ )
  !
  RETURN
  !
  END SUBROUTINE compute_becpsi
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy_bp_k()
    !-----------------------------------------------------------------------
    !
    USE constants,               ONLY : fpi, e2, pi
    USE io_files,                ONLY : iunwfc_exx, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g
    USE wvfct,                   ONLY : nbnd, npwx, wg
    USE wavefunctions,           ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_exx,                  ONLY : inter_egrp_comm, my_egrp_id, negrp, &
                                        intra_egrp_comm, me_egrp, &
                                        max_pairs, egrp_pairs, ibands, nibands, &
                                        iexx_istart, iexx_iend, &
                                        all_start, all_end, iexx_start, &
                                        init_index_over_band, jblock
    USE mp_bands,                ONLY : intra_bgrp_comm
    USE mp,                      ONLY : mp_sum, mp_circular_shift_left
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, &
                                        deallocate_bec_type, calbec
    USE paw_variables,           ONLY : okpaw
    USE paw_exx,                 ONLY : PAW_xx_energy
    USE us_exx,                  ONLY : bexg_merge, becxx, addusxx_g, &
                                        addusxx_r, qvan_init, qvan_clean
    USE exx_base,                ONLY : nqs, xkq_collect, index_xkq, index_xk
    USE exx_bp_utils,            ONLY : change_data_structure, &
                                        transform_evc_to_exx, nwordwfc_exx, &
                                        igk_exx, evc_exx
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy_bp_k
    !
    ! ... local variables
    !
    REAL(DP) :: energy
    COMPLEX(DP), ALLOCATABLE :: temppsic(:,:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_nc(:,:,:)
    COMPLEX(DP), ALLOCATABLE,TARGET :: rhoc(:,:)
#if defined(__USE_MANY_FFT)
    COMPLEX(DP), POINTER :: prhoc(:)
#endif
    REAL(DP),    ALLOCATABLE :: fac(:)
    INTEGER  :: npw, jbnd, ibnd, ibnd_inner_start, ibnd_inner_end, ibnd_inner_count, &
                ik, ikk, ig, ikq, iq, ir
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start, nblock, nrt, irt, &
                ir_start, ir_end
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc, omega_inv
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER :: intra_bgrp_comm_
    INTEGER :: ii, ialloc, jstart, jend, ipair
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    !
    CALL init_index_over_band( inter_egrp_comm, nbnd, nbnd )
    !
    CALL transform_evc_to_exx( 0 )
    !
    ialloc = nibands(my_egrp_id+1)
    !
    nrxxs = dfftt%nnr
    ALLOCATE( fac(dfftt%ngm) )
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol,ialloc) )
    ELSE
       ALLOCATE( temppsic(nrxxs,ialloc) )
    ENDIF
    !
    energy = 0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi )
    !
    !precompute that stuff
    omega_inv = 1.0/omega
    !
    IKK_LOOP : &
    DO ikk = 1, nks
       !
       current_ik = global_kpoint_index ( nkstot, ikk )
       xkp = xk(:,ikk)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk(ikk)
       IF ( nks > 1 ) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ikk )
       !
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       intra_bgrp_comm_ = intra_bgrp_comm
       intra_bgrp_comm = intra_egrp_comm
       !
       IF (okvan .OR. okpaw) THEN
          !! FIXME: can be replaced by a call to init_us_2 + calbec
          CALL compute_becpsi( npw, igk_exx(:,ikk), xkp, evc_exx, &
                               becpsi%k(:,ibands(1,my_egrp_id+1)) )
       ENDIF
       !
       intra_bgrp_comm = intra_bgrp_comm_
       !
       ! ... precompute temppsic
       !
       IF (noncolin) THEN
          temppsic_nc = 0.0_DP
       ELSE
          temppsic = 0.0_DP
       ENDIF
       !
       DO ii = 1, nibands(my_egrp_id+1)
          !
          jbnd = ibands(ii,my_egrp_id+1)
          !
          IF (jbnd == 0 .OR. jbnd > nbnd) CYCLE
          !
          !IF ( abs(wg(jbnd,ikk)) < eps_occ) CYCLE
          !
          IF (noncolin) THEN
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic_nc(dfftt%nl(igk_exx(ig,ikk)),1,ii) = evc_exx(ig,ii)
                temppsic_nc(dfftt%nl(igk_exx(ig,ikk)),2,ii) = evc_exx(npwx+ig,ii)
             ENDDO
!$omp end parallel do
             !
             CALL invfft( 'Wave', temppsic_nc(:,1,ii), dfftt )
             CALL invfft( 'Wave', temppsic_nc(:,2,ii), dfftt )
             !
          ELSE
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic(dfftt%nl(igk_exx(ig,ikk)),ii) = evc_exx(ig,ii)
             ENDDO
!$omp end parallel do
             !
             CALL invfft( 'Wave', temppsic(:,ii), dfftt )
             !
          ENDIF
       ENDDO
       !
       IQ_LOOP : &
       DO iq = 1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          !
          xkq = xkq_collect(:,ikq)
          !
          CALL g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, ikk )
          IF ( okvan .AND..NOT.tqr ) CALL qvan_init( dfftt%ngm, xkq, xkp )
          !
          DO iegrp = 1, negrp
             !
             ! ... compute the id of group whose data is currently worked on
             wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
             njt = (all_end(wegrp)-all_start(wegrp)+jblock)/jblock
             !
             IJT_LOOP : &
             DO ijt = 1, njt
                !
                jblock_start = (ijt - 1) * jblock + all_start(wegrp)
                jblock_end = MIN(jblock_start+jblock-1,all_end(wegrp))
                !
                JBND_LOOP : &
                DO ii = 1, nibands(my_egrp_id+1)
                   !
                   jbnd = ibands(ii,my_egrp_id+1)
                   !
                   IF (jbnd==0 .OR. jbnd>nbnd) CYCLE
                   !
                   !determine which j-bands to calculate
                   jstart = 0
                   jend = 0
                   !
                   DO ipair = 1, max_pairs
                      IF (egrp_pairs(1,ipair,my_egrp_id+1) == jbnd) THEN
                         IF (jstart == 0) THEN
                            jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                            jend = jstart
                         ELSE
                            jend = egrp_pairs(2,ipair,my_egrp_id+1)
                         ENDIF
                      ENDIF
                   ENDDO
                   !
                   !these variables prepare for inner band parallelism
                   jstart = MAX(jstart,jblock_start)
                   jend = MIN(jend,jblock_end)
                   ibnd_inner_start = jstart
                   ibnd_inner_end = jend
                   ibnd_inner_count = jend-jstart+1
                   !
                   !allocate arrays
                   ALLOCATE( rhoc(nrxxs,ibnd_inner_count) )
#if defined(__USE_MANY_FFT)
                   prhoc(1:nrxxs*ibnd_inner_count) => rhoc
#endif 
                   !calculate rho in real space
                   nblock = 2048
                   nrt = (nrxxs+nblock-1) / nblock
!$omp parallel do collapse(2) private(ir_start,ir_end)
                   DO irt = 1, nrt
                      DO ibnd = ibnd_inner_start, ibnd_inner_end
                         ir_start = (irt - 1) * nblock + 1
                         ir_end = MIN(ir_start+nblock-1,nrxxs)
                         IF (noncolin) THEN
                            DO ir = ir_start, ir_end
                               rhoc(ir,ibnd-ibnd_inner_start+1) = &
                                 ( CONJG(exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)) * &
                                 temppsic_nc(ir,1,ii) + &
                                 CONJG(exxbuff(nrxxs+ir,ibnd-all_start(wegrp)+iexx_start,ikq)) * &
                                 temppsic_nc(ir,2,ii) ) * omega_inv
                            ENDDO
                         ELSE
                            DO ir = ir_start, ir_end
                               rhoc(ir,ibnd-ibnd_inner_start+1) = omega_inv * &
                                 CONJG(exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)) * &
                                 temppsic(ir,ii)
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
!$omp end parallel do
                   !
                   ! augment the "charge" in real space
                   IF (okvan .AND. tqr) THEN
!$omp parallel do default(shared) private(ibnd) firstprivate(ibnd_inner_start,ibnd_inner_end)
                      DO ibnd = ibnd_inner_start, ibnd_inner_end
                         CALL addusxx_r( rhoc(:,ibnd-ibnd_inner_start+1), &
                                        becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   ! bring rhoc to G-space
#if defined(__USE_MANY_FFT)
                   CALL fwfft ('Rho', prhoc, dfftt, howmany=ibnd_inner_count)
#else
                   DO ibnd = ibnd_inner_start, ibnd_inner_end
                      CALL fwfft('Rho', rhoc(:,ibnd-ibnd_inner_start+1), dfftt)
                   ENDDO
#endif
                   ! augment the "charge" in G space
                   IF (okvan .AND. .NOT. tqr) THEN
                      DO ibnd = ibnd_inner_start, ibnd_inner_end
                         CALL addusxx_g(dfftt, rhoc(:,ibnd-ibnd_inner_start+1), &
                              xkq, xkp, 'c', becphi_c=becxx(ikq)%k(:,ibnd),     &
                              becpsi_c=becpsi%k(:,jbnd))
                      ENDDO
                   ENDIF
                   !
!$omp parallel do reduction(+:energy) private(vc)
                   DO ibnd = ibnd_inner_start, ibnd_inner_end
                      vc=0.0_DP
                      DO ig=1,dfftt%ngm
                         vc = vc + coulomb_fac(ig,iq,ikk) * &
                             DBLE(rhoc(dfftt%nl(ig),ibnd-ibnd_inner_start+1) *&
                             CONJG(rhoc(dfftt%nl(ig),ibnd-ibnd_inner_start+1)))
                      ENDDO
                      vc = vc * omega * x_occupation(ibnd,ik) / nqs
                      energy = energy - exxalfa * vc * wg(jbnd,ikk)
                      !
                      IF (okpaw) THEN
                         energy = energy +exxalfa*x_occupation(ibnd,ik)/nqs*wg(jbnd,ikk) &
                              *PAW_xx_energy(becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                      ENDIF
                   ENDDO
!$omp end parallel do
                   !
                   !deallocate memory
                   DEALLOCATE( rhoc )
                ENDDO &
                JBND_LOOP
                !
             ENDDO&
             IJT_LOOP
             ! get the next nbnd/negrp data
             IF (negrp > 1) call mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
             !
          END DO !iegrp
          !
          IF ( okvan .AND. .NOT.tqr ) CALL qvan_clean()
       ENDDO &
       IQ_LOOP
       !
    ENDDO &
    IKK_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc )
    ELSE
       DEALLOCATE( temppsic )
    ENDIF
    !
    DEALLOCATE( fac )
    !
    CALL deallocate_bec_type( becpsi )
    !
    CALL mp_sum( energy, inter_egrp_comm )
    CALL mp_sum( energy, intra_egrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy_bp_k = energy
    CALL change_data_structure( .FALSE. )
    !
  END FUNCTION  exxenergy_bp_k
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_stress_bp()
    !-----------------------------------------------------------------------
    !! This is Eq.(10) of PRB 73, 125120 (2006).
    !
    USE constants,            ONLY : fpi, e2, pi, tpi
    USE io_files,             ONLY : iunwfc_exx, nwordwfc
    USE buffers,              ONLY : get_buffer
    USE cell_base,            ONLY : alat, omega, bg, at, tpiba
    USE symm_base,            ONLY : nsym, s
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE wavefunctions,        ONLY : evc
    USE klist,                ONLY : xk, ngk, nks
    USE lsda_mod,             ONLY : lsda, current_spin, isk
    USE gvect,                ONLY : g
    USE mp_pools,             ONLY : npool, inter_pool_comm
    USE mp_exx,               ONLY : inter_egrp_comm, intra_egrp_comm, &
                                     ibands, nibands, my_egrp_id, jblock, &
                                     egrp_pairs, max_pairs, negrp, me_egrp, &
                                     all_start, all_end, iexx_start
    USE mp,                   ONLY : mp_sum, mp_circular_shift_left
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE uspp,                 ONLY : okvan
    !
    USE exx_base,             ONLY : nq1, nq2, nq3, nqs, eps, exxdiv,       &
                                     x_gamma_extrapolation, on_double_grid, &
                                     grid_factor, yukawa, erfc_scrlen,      &
                                     use_coulomb_vcut_ws, use_coulomb_vcut_spheric, &
                                     gau_scrlen, vcut, index_xkq, index_xk, index_sym
    USE exx_bp_utils,         ONLY : change_data_structure, transform_evc_to_exx, &
                                     g_exx, igk_exx, nwordwfc_exx, evc_exx
    USE coulomb_vcut_module,  ONLY : vcut_get,  vcut_spheric_get
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    REAL(DP) :: exx_stress_bp(3,3), exx_stress_(3,3)
    !
    COMPLEX(DP),ALLOCATABLE :: tempphic(:), temppsic(:), result(:)
    COMPLEX(DP),ALLOCATABLE :: tempphic_nc(:,:), temppsic_nc(:,:), &
                               result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: rhoc(:)
    REAL(DP),   ALLOCATABLE :: fac(:), fac_tens(:,:,:), fac_stress(:)
    INTEGER  :: npw, jbnd, ibnd, ik, ikk, ig, ir, ikq, iq, isym
    INTEGER  :: nqi, iqi, beta, nrxxs, ngm
    INTEGER  :: ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    REAL(DP) :: delta(3,3)
    INTEGER :: jstart, jend, ii, ipair, jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    INTEGER :: exxbuff_index
    !
    CALL transform_evc_to_exx( 0 )
    !
    IF (npool>1) CALL errore( 'exx_stress2', 'stress not available with pools', 1 )
    IF (noncolin) CALL errore( 'exx_stress2', 'noncolinear stress not implemented', 1 )
    IF (okvan) CALL infomsg( 'exx_stress2', 'USPP stress not tested' )
    !
    nrxxs = dfftt%nnr
    ngm   = dfftt%ngm
    delta = RESHAPE( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    !
    ALLOCATE( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    ALLOCATE( fac_tens(3,3,ngm), fac_stress(ngm) )
    !
    nqi = nqs
    !
    ! ... loop over k-points
    DO ikk = 1, nks
       current_k = ikk
       IF (lsda) current_spin = isk(ikk)
       npw = ngk(ikk)
       !
       IF (nks > 1) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ikk )
       !
       DO iqi = 1, nqi
          !
          iq = iqi
          !
          ikq  = index_xkq(current_k,iq)
          ik   = index_xk(ikq)
          isym = ABS(index_sym(ikq))      
          !      

          ! FIXME: use cryst_to_cart and company as above..      
          xk_cryst(:) = at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)      
          IF (index_sym(ikq) < 0) xk_cryst = -xk_cryst      
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &      
                   s(:,2,isym)*xk_cryst(2) + &      
                   s(:,3,isym)*xk_cryst(3)      
          xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)      
          !      
          !CALL start_clock ('exxen2_ngmloop')      
          !      
!$omp parallel do default(shared), private(ig, beta, q, qq, on_double_grid, x)
          DO ig = 1, ngm      
             IF (negrp == 1) THEN      
                q(1) = xk(1,current_k) - xkq(1) + g(1,ig)      
                q(2) = xk(2,current_k) - xkq(2) + g(2,ig)      
                q(3) = xk(3,current_k) - xkq(3) + g(3,ig)      
             ELSE      
                q(1) = xk(1,current_k) - xkq(1) + g_exx(1,ig)      
                q(2) = xk(2,current_k) - xkq(2) + g_exx(2,ig)      
                q(3) = xk(3,current_k) - xkq(3) + g_exx(3,ig)      
             ENDIF      
             !      
             q = q * tpiba      
             qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )      
             !      
             DO beta = 1, 3      
                fac_tens(1:3,beta,ig) = q(1:3)*q(beta)      
             ENDDO      
             !      
             IF (x_gamma_extrapolation) THEN      
                on_double_grid = .TRUE.      
                x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1      
                on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)      
                x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2      
                on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)      
                x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3      
                on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)      
             ELSE      
                on_double_grid = .FALSE.      
             ENDIF      
             !      
             IF (use_coulomb_vcut_ws) THEN      
                fac(ig) = vcut_get(vcut, q)      
                fac_stress(ig) = 0._dp   ! not implemented      
                IF (gamma_only .AND. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)      
                !      
             ELSEIF ( use_coulomb_vcut_spheric ) THEN      
                fac(ig) = vcut_spheric_get(vcut, q)      
                fac_stress(ig) = 0._dp   ! not implemented      
                IF (gamma_only .AND. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)      
                !      
             ELSEIF (gau_scrlen > 0) THEN      
                fac(ig) = e2*((pi/gau_scrlen)**(1.5d0))* &       
                          EXP(-qq/4.d0/gau_scrlen) * grid_factor       
                fac_stress(ig) =  e2*2.d0/4.d0/gau_scrlen  * &       
                                  EXP(-qq/4.d0/gau_scrlen) * &
                                  ((pi/gau_scrlen)**(1.5d0))*grid_factor       
                IF (gamma_only) fac(ig) = 2.d0 * fac(ig)       
                IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)       
                IF (on_double_grid) fac(ig) = 0._dp       
                IF (on_double_grid) fac_stress(ig) = 0._dp       
                !       
             ELSEIF (qq > 1.d-8) THEN      
                IF ( erfc_scrlen > 0 ) THEN       
                  fac(ig)=e2*fpi/qq*(1._dp-EXP(-qq/4.d0/erfc_scrlen**2)) * grid_factor       
                  fac_stress(ig) = -e2*fpi * 2.d0/qq**2 * ( &       
                      (1._dp+qq/4.d0/erfc_scrlen**2)*EXP(-qq/4.d0/erfc_scrlen**2) - 1._dp) * &       
                      grid_factor       
                ELSE       
                  fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor       
                  fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2 * grid_factor       
                ENDIF       
                !       
                IF (gamma_only) fac(ig) = 2.d0 * fac(ig)       
                IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)       
                IF (on_double_grid) fac(ig) = 0._dp       
                IF (on_double_grid) fac_stress(ig) = 0._dp       
                !      
             ELSE 
                ! 
                fac(ig) = -exxdiv ! or rather something else (see f.gygi)       
                fac_stress(ig) = 0._dp  ! or -exxdiv_stress (not yet implemented)       
                IF ( yukawa> 0._dp .AND. .NOT. x_gamma_extrapolation) THEN       
                  fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )       
                  fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2       
                ENDIF       
                IF (erfc_scrlen > 0._dp .AND. .NOT. x_gamma_extrapolation) THEN       
                  fac(ig) = e2*fpi / (4.d0*erfc_scrlen**2)       
                  fac_stress(ig) = e2*fpi / (8.d0*erfc_scrlen**4)       
                ENDIF 
                !
             ENDIF
             !
          ENDDO      
!$omp end parallel do
          !CALL stop_clock ('exxen2_ngmloop')      
          DO iegrp = 1, negrp      
             !      
             ! compute the id of group whose data is currently worked on      
             wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1      
             !      
             jblock_start = all_start(wegrp)      
             jblock_end = all_end(wegrp)      
             !      
             ! loop over bands      
             DO ii = 1, nibands(my_egrp_id+1)      
                !      
                jbnd = ibands(ii,my_egrp_id+1)      
                !      
                IF (jbnd==0 .OR. jbnd>nbnd) CYCLE      
                !      
                !determine which j-bands to calculate      
                jstart = 0      
                jend = 0      
                DO ipair=1, max_pairs      
                   IF (egrp_pairs(1,ipair,my_egrp_id+1).eq.jbnd)THEN      
                      IF (jstart == 0)THEN      
                         jstart = egrp_pairs(2,ipair,my_egrp_id+1)      
                         jend = jstart      
                      ELSE      
                         jend = egrp_pairs(2,ipair,my_egrp_id+1)      
                      ENDIF      
                   ENDIF      
                ENDDO      
                !      
                jstart = MAX(jstart,jblock_start)      
                jend = MIN(jend,jblock_end)      
                !      
                temppsic(:) = ( 0._dp, 0._dp )      
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw      
                   temppsic(dfftt%nl(igk_exx(ig,ikk))) = evc_exx(ig,ii)      
                ENDDO      
!$omp end parallel do
                !      
                IF (gamma_only) THEN      
!$omp parallel do default(shared), private(ig)
                   DO ig = 1, npw      
                      temppsic(dfftt%nlm(igk_exx(ig,ikk))) = CONJG(evc_exx(ig,ii))      
                   ENDDO      
!$omp end parallel do
                ENDIF      
                !      
                CALL invfft( 'Wave', temppsic, dfftt )      
                !      
                IF (gamma_only) THEN      
                   !      
                   IF (MOD(jstart,2) == 0) THEN      
                      ibnd_loop_start = jstart-1      
                   ELSE      
                      ibnd_loop_start = jstart      
                   ENDIF      
                   !      
                   DO ibnd = ibnd_loop_start, jend, 2     !for each band of psi      
                      !      
                      exxbuff_index = (ibnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2      
                      !      
                      IF ( ibnd < jstart ) THEN      
                         x1 = 0._dp      
                      ELSE      
                         x1 = x_occupation(ibnd,ik)      
                      ENDIF      
                      !      
                      IF ( ibnd == jend) THEN      
                         x2 = 0._dp      
                      ELSE      
                         x2 = x_occupation(ibnd+1,ik)      
                      ENDIF      
                      !      
                      IF ( ABS(x1) < eps_occ .AND. ABS(x2) < eps_occ ) CYCLE      
                      !      
                      ! calculate rho in real space      
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs      
                         tempphic(ir) = exxbuff(ir,exxbuff_index,ikq)      
                         rhoc(ir) = CONJG(tempphic(ir))*temppsic(ir) / omega      
                      ENDDO      
!$omp end parallel do
                      ! bring it to G-space      
                      CALL fwfft( 'Rho', rhoc, dfftt )      
                      !      
                      vc = 0._dp      
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm      
                         !      
                         vc(:,:) = vc(:,:) + x1 * 0.25_dp * &      
                                   ABS( rhoc(dfftt%nl(ig)) + &      
                                   CONJG(rhoc(dfftt%nlm(ig))))**2 * &      
                                   (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))      
                         vc(:,:) = vc(:,:) + x2 * 0.25_dp * &      
                                   ABS( rhoc(dfftt%nl(ig)) - &      
                                   CONJG(rhoc(dfftt%nlm(ig))))**2 * &      
                                   (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))      
                      ENDDO      
!$omp end parallel do
                      vc = vc / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                   ENDDO
                   !
                ELSE
                   !
                   DO ibnd = jstart, jend    !for each band of psi
                      !
                      IF ( ABS(x_occupation(ibnd,ik)) < 1.d-6) CYCLE
                      !
                      ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         tempphic(ir) = exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)
                         rhoc(ir) = CONJG(tempphic(ir))*temppsic(ir) / omega
                      ENDDO
!$omp end parallel do
                      !
                      ! bring it to G-space
                      CALL fwfft( 'Rho', rhoc, dfftt )
                      !
                      vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm
                         vc(:,:) = vc(:,:) + rhoc(dfftt%nl(ig))  * &
                                   CONJG(rhoc(dfftt%nl(ig))) *     &
                                   (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - &
                                   delta(:,:)*fac(ig))
                      ENDDO
!$omp end parallel do
                      !
                      vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                      !
                   ENDDO
                   !
                ENDIF ! gamma or k-points
                !
             ENDDO ! jbnd
             !
             ! get the next nbnd/negrp data
             IF (negrp > 1) CALL mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, &
                                                         inter_egrp_comm )
             !
          ENDDO ! iegrp
          !
       ENDDO ! iqi
       !
    ENDDO ! ikk
    !
    DEALLOCATE( tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
    !
    CALL mp_sum( exx_stress_, intra_egrp_comm )
    CALL mp_sum( exx_stress_, inter_egrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    !
    exx_stress_bp = exx_stress_
    !
    CALL change_data_structure( .FALSE. )
    !
  END FUNCTION exx_stress_bp
  !
END MODULE exx_bp
