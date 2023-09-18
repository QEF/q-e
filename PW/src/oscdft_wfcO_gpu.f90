MODULE oscdft_wfcO_gpu
#if defined (__OSCDFT)
   USE kinds, ONLY : DP
   USE oscdft_wavefunction, ONLY : oscdft_wavefunction_type
   USE oscdft_indices,      ONLY : oscdft_indices_type,&
                                   oscdft_constr_indices_type,&
                                   oscdft_orbital_indices_type
   USE oscdft_input,        ONLY : oscdft_input_type
   USE oscdft_context,      ONLY : oscdft_context_type

   PRIVATE
   PUBLIC orthoOwfc_gpu
   CONTAINS
      SUBROUTINE orthoOwfc_gpu(ctx)
         USE buffers,                  ONLY : get_buffer, save_buffer
         USE klist,                    ONLY : nks, xk, ngk, igk_k
         USE wvfct,                    ONLY : npwx
         USE uspp,                     ONLY : nkb, vkb
         USE uspp_param,               ONLY : upf
         USE becmod,                   ONLY : allocate_bec_type,&
                                              deallocate_bec_type,&
                                              becp, calbec, bec_type
         USE becmod_gpum,              ONLY : becp_d
         USE becmod_subs_gpum,         ONLY : using_becp_auto, using_becp_d_auto, calbec_gpu
         USE noncollin_module,         ONLY : npol
         USE basis,                    ONLY : natomwfc
         USE ions_base,                ONLY : nat, ityp
         USE symm_base,                ONLY : nsym
         USE uspp_init,                ONLY : init_us_2
         USE oscdft_wavefunction,      ONLY : check_bec_type_unallocated
         USE oscdft_wavefunction_subs, ONLY : oscdft_fill_wavefunctions,&
                                              oscdft_poolrecover_overlap
         USE oscdft_wavefunction_subs_gpu,&
                                       ONLY : oscdft_ortho_swfc_gpu,&
                                              oscdft_get_overlap_gpu,&
                                              oscdft_lowdin_ortho_overlap_gpu
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         INTEGER                                          :: ik, npw
         TYPE(oscdft_input_type),           POINTER       :: inp
         TYPE(oscdft_indices_type),         POINTER       :: idx
         TYPE(oscdft_orbital_indices_type), POINTER       :: orbs
         TYPE(oscdft_constr_indices_type),  POINTER       :: constr
         TYPE(oscdft_wavefunction_type),    POINTER       :: wfcO, wfcS, wfcF

         COMPLEX(DP), ALLOCATABLE                         :: wfcatom(:,:), swfcatom(:,:)

         IF (.NOT.ctx%wfc_allocated) THEN
            CALL errore("oscdft_orthoOwfc", "call oscdft_init_wfcO first", 1)
         END IF
         CALL start_clock("oscdft_orthoOwfc_gpu")
         ! WRITE(stdout, *) "orthoOwfc"
         inp    => ctx%inp
         idx    => ctx%idx
         orbs   => idx%orbs
         constr => idx%constr

         wfcO => ctx%wfcO
         wfcS => ctx%wfcS
         wfcF => ctx%forces%wfcO

         ctx%nst%eigvects_set = .false.
         ALLOCATE(wfcatom(npwx*npol,natomwfc), swfcatom(npwx*npol,natomwfc))
         !$acc data create(wfcatom, swfcatom)
         CALL check_bec_type_unallocated(becp)
         CALL allocate_bec_type(nkb, natomwfc, becp)
         CALL using_becp_auto(2)

         DO ik=1,nks
            !$acc host_data use_device(wfcatom, swfcatom)
            CALL atomic_wfc_gpu(ik, wfcatom)
            npw = ngk(ik)
            CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb, .true.)
            CALL using_becp_d_auto(2)
            !$acc data present(vkb(:,:))
            !$acc host_data use_device(vkb)
            CALL calbec_gpu(npw, vkb, wfcatom, becp_d)
            !$acc end host_data
            !$acc end data
            CALL s_psi_gpu(npwx, npw, natomwfc, wfcatom, swfcatom)
            !$acc end host_data

            IF (inp%orthogonalize_swfc) THEN
               CALL oscdft_ortho_swfc_gpu(npwx, npw, natomwfc, wfcatom, swfcatom, .false.)
            ELSE IF (inp%normalize_swfc) THEN
               CALL oscdft_ortho_swfc_gpu(npwx, npw, natomwfc, wfcatom, swfcatom, .true.)
            END IF

            !$acc update host(wfcatom, swfcatom)
            CALL oscdft_fill_wavefunctions(idx, wfcO, ik, wfcatom, swfcatom)
            CALL oscdft_fill_wavefunctions(idx, wfcS, ik, wfcatom, swfcatom)
            CALL oscdft_fill_wavefunctions(idx, wfcF, ik, wfcatom, swfcatom)

            IF (inp%orthogonalize_ns) THEN
               CALL oscdft_get_overlap_gpu(ctx, ik, wfcS%source, wfcatom, swfcatom)
            END IF
         END DO
         CALL deallocate_bec_type(becp)
         CALL using_becp_auto(2)
         !$acc end data
         DEALLOCATE(wfcatom, swfcatom)

         IF (inp%orthogonalize_ns) THEN
            CALL oscdft_poolrecover_overlap(idx)
            CALL oscdft_lowdin_ortho_overlap_gpu(ctx)
         END IF

         CALL stop_clock("oscdft_orthoOwfc_gpu")
      END SUBROUTINE orthoOwfc_gpu
#endif
END MODULE oscdft_wfcO_gpu
