MODULE oscdft_functions_gpu
#if defined (__OSCDFT)
   USE kinds,                    ONLY : DP
   USE io_global,                ONLY : ionode, ionode_id, stdout
   USE mp,                       ONLY : mp_bcast
   USE mp_images,                ONLY : intra_image_comm
   USE oscdft_enums
   USE oscdft_occupations,       ONLY : oscdft_new_ns,&
                                        oscdft_get_occupation_numbers
   USE oscdft_context,           ONLY : oscdft_context_type, oscdft_ns_type
   USE oscdft_input,             ONLY : oscdft_input_type
   USE oscdft_indices,           ONLY : oscdft_indices_type,&
                                        oscdft_orbital_indices_type,&
                                        oscdft_constr_indices_type
   USE oscdft_wavefunction,      ONLY : oscdft_wavefunction_type,&
                                        check_bec_type_unallocated_gpu
   USE oscdft_wavefunction_subs, ONLY : oscdft_get_buffer

   PRIVATE
   PUBLIC oscdft_h_diag_gpu, oscdft_h_psi_gpu

   CONTAINS
      SUBROUTINE oscdft_h_diag_gpu(ctx)
         USE lsda_mod,         ONLY : isk, current_spin
         USE klist,            ONLY : ngk
         USE klist,            ONLY : nks
         USE buffers,          ONLY : get_buffer
         USE wvfct,            ONLY : current_k, npwx
         USE noncollin_module, ONLY : npol
         USE g_psi_mod_gpum,   ONLY : h_diag_d, using_h_diag_d

         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),           POINTER       :: inp
         TYPE(oscdft_indices_type),         POINTER       :: idx
         TYPE(oscdft_ns_type),              POINTER       :: nst
         TYPE(oscdft_wavefunction_type),    POINTER       :: wfcO

         REAL(DP)                                         :: alpha
         INTEGER                                          :: curr_dim, h, k, i, ik,&
                                                             iconstr, ioscdft, oidx,&
                                                             isum, osum, h_off, k_off,&
                                                             npw
         COMPLEX(DP), POINTER :: wfcO_wfc(:,:)

         inp => ctx%inp
         idx => ctx%idx
         nst => ctx%nst
         wfcO => ctx%wfcO

         IF (.NOT.ctx%initialized) RETURN
         IF (ctx%warming_up) RETURN
         IF (idx%nconstr == 0) RETURN
         IF (.NOT.ANY(inp%spin_index(idx%iconstr2ioscdft(:)) == current_spin)) RETURN

         CALL start_clock_gpu("oscdft_hdiag")
         npw = ngk(current_k)

         IF (.NOT.ctx%wfc_allocated) THEN
            CALL errore("oscdft_h_psi", "wfc not allocated", 1)
         END IF

         ik = current_k
         CALL using_h_diag_d(1)
         IF (isk(ik) == current_spin) THEN
            IF (nks > 1) THEN
               CALL oscdft_get_buffer(wfcO, ik)
            END IF
            wfcO_wfc => wfcO%wfc
            !$acc data copyin(wfcO_wfc)

            DO iconstr=1,idx%nconstr
               ioscdft = idx%iconstr2ioscdft(iconstr)
               oidx = inp%occup_index(ioscdft)
               IF (inp%start_index(ioscdft) > ctx%global_start_index) CYCLE
               IF (inp%spin_index(ioscdft) == current_spin) THEN
                  curr_dim = idx%ns_dim(ioscdft)
                  DO h=1,curr_dim
                     h_off = wfcO%get_offset(idx%constr, h, iconstr, -1)
                     DO k=1,curr_dim
                        k_off = wfcO%get_offset(idx%constr, k, iconstr, -1)
                        IF (oidx == OCCUP_TRACE) THEN
                           alpha = MERGE(ctx%multipliers(iconstr), 0.D0, h == k)
                        ELSE IF (oidx == OCCUP_SUM) THEN
                           osum = inp%occup_index_sum(2,ioscdft)
                           alpha = ctx%multipliers(iconstr)*&
                                   ctx%nst%occup_eigvects(h,osum,iconstr)*&
                                   ctx%nst%occup_eigvects(k,osum,iconstr)
                        ELSE
                           alpha = ctx%multipliers(iconstr)*&
                                   ctx%nst%occup_eigvects(h,oidx,iconstr)*&
                                   ctx%nst%occup_eigvects(k,oidx,iconstr)
                        END IF
                        !$acc data present(h_diag_d)
                        !$acc data present(wfcO_wfc)
                        !$acc parallel loop
                        DO i=1,npw
                           h_diag_d(i,1) = h_diag_d(i,1)+&
                                         alpha*DBLE(wfcO_wfc(i,k_off)*CONJG(wfcO_wfc(i,h_off)))
                        END DO
                        !$acc end data
                        !$acc end data
                     END DO
                  END DO
               END IF
            END DO
            !$acc end data
         END IF
         CALL stop_clock_gpu("oscdft_hdiag")
      END SUBROUTINE oscdft_h_diag_gpu

      SUBROUTINE oscdft_h_psi_gpu(ctx, lda, n, m, psi_d, hpsi_d)
         USE noncollin_module,    ONLY : npol
         USE klist,               ONLY : nks, nelup, neldw
         USE buffers,             ONLY : get_buffer
         USE lsda_mod,            ONLY : isk, current_spin
         USE io_files,            ONLY : nwordwfc
         USE gvect,               ONLY : gstart
         USE control_flags,       ONLY : gamma_only, offload_type
         USE wvfct,               ONLY : btype, current_k, nbnd, wg
         USE mp_bands,            ONLY : intra_bgrp_comm
         USE mp,                  ONLY : mp_barrier
         USE becmod,              ONLY : bec_type, allocate_bec_type_acc,&
                                         deallocate_bec_type_acc
#if defined(__CUDA)
         USE becmod,              ONLY : calbec_cuf
#endif
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),        POINTER          :: inp
         TYPE(oscdft_indices_type),      POINTER          :: idx
         TYPE(oscdft_wavefunction_type), POINTER          :: wfcO
         INTEGER,                   INTENT(IN)            :: lda, n, m
         COMPLEX(DP),               INTENT(IN)            :: psi_d(lda*npol,m)
         COMPLEX(DP),               INTENT(INOUT)         :: hpsi_d(lda*npol,m)
         INTEGER                                          :: iconstr, ioscdft,&
                                                             ik, h, k, curr_dim,&
                                                             ibnd, i, oidx, isum, osum, jsum,&
                                                             h_off, k_off
         REAL(DP)                                         :: occ_const
         COMPLEX(DP)                                      :: alpha

         COMPLEX(DP),    POINTER :: wfcO_wfc(:,:)
         TYPE(bec_type), TARGET  :: proj
         REAL(DP),       POINTER :: proj_r(:,:)
         COMPLEX(DP),    POINTER :: proj_k(:,:)
#if defined(__CUDA)
         attributes(DEVICE) :: psi_d, hpsi_d
#endif

         inp => ctx%inp
         idx => ctx%idx
         wfcO => ctx%wfcO

         IF (.NOT.ctx%initialized) RETURN
         IF (ctx%warming_up) RETURN
         IF (idx%nconstr == 0) RETURN


         CALL start_clock_gpu("oscdft_hpsi")
         ! sum_hk u_h^I u_k^I |phi_k> <phi_h|psi>
         ik = current_k
         CALL allocate_bec_type_acc(m, wfcO%n, proj)
         IF (gamma_only) THEN
            proj_r => proj%r
         ELSE
            proj_k => proj%k
         END IF
         IF (nks > 1) CALL oscdft_get_buffer(wfcO, ik)
         wfcO_wfc => wfcO%wfc
         !$acc data copyin(wfcO_wfc(:,:))

         ! gets <psi|phi_h>, then call ZGERC where this term is turned to its complex conjugate (<phi_h|psi>)
#if defined(__CUDA)
!civn: remove evc_d and use calbec instead
         CALL calbec_cuf(offload_type, n, psi_d, wfcO_wfc, proj)
#endif

         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            oidx = inp%occup_index(ioscdft)
            IF (inp%start_index(ioscdft) > ctx%global_start_index) CYCLE
            IF ((inp%constraint_applied(ioscdft) /= CONSTR_FALSE).AND.&
                (inp%spin_index(ioscdft) == current_spin)) THEN
               curr_dim = idx%ns_dim(ioscdft)
               DO k=1,curr_dim
                  k_off = wfcO%get_offset(idx%constr, k, iconstr, -1)
                  DO h=1,curr_dim
                     h_off = wfcO%get_offset(idx%constr, h, iconstr, -1)
                     IF (oidx == OCCUP_TRACE) THEN
                        occ_const = MERGE(1.D0, 0.D0, h == k)
                     ELSE IF (oidx == OCCUP_SUM) THEN
                        osum = inp%occup_index_sum(2,ioscdft)
                        occ_const = ctx%nst%occup_eigvects(h,osum,iconstr)*&
                                    ctx%nst%occup_eigvects(k,osum,iconstr)
                     ELSE
                        occ_const = ctx%nst%occup_eigvects(h,oidx,iconstr)*&
                                    ctx%nst%occup_eigvects(k,oidx,iconstr)
                     ENDIF
                     alpha = CMPLX(occ_const * ctx%multipliers(iconstr), 0.D0, kind=DP)
                     IF (gamma_only) THEN
                        !$acc data present(wfcO_wfc(:,:), hpsi_d)
                        !$acc host_data use_device(wfcO_wfc(:,:), proj_r)
                        CALL MYDGER(2*n,m,DBLE(alpha),&
                                        wfcO_wfc(1,h_off),1,&
                                        proj_r(1:m,k_off),1,&
                                        hpsi_d,2*lda)
                        !$acc end host_data
                        !$acc end data
                     ELSE
                        !$acc data present(wfcO_wfc(:,:), hpsi_d)
                        !$acc host_data use_device(wfcO_wfc(:,:), proj_k)
                        CALL MYZGERC(n,m,alpha,&
                                     wfcO_wfc(1,h_off),1,&
                                     proj_k(1:m,k_off),1,&
                                     hpsi_d,lda)
                        !$acc end host_data
                        !$acc end data
                     END IF
                  END DO
               END DO
            END IF
         END DO
         !$acc end data

         CALL deallocate_bec_type_acc(proj)
         CALL mp_barrier(intra_bgrp_comm)
         CALL stop_clock_gpu("oscdft_hpsi")
      END SUBROUTINE oscdft_h_psi_gpu
#endif
END MODULE oscdft_functions_gpu
