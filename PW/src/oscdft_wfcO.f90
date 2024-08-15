MODULE oscdft_wfcO
#if defined (__OSCDFT)
   USE kinds,               ONLY : DP
   USE io_global,           ONLY : stdout
   USE oscdft_wavefunction, ONLY : oscdft_wavefunction_type
   USE oscdft_indices,      ONLY : oscdft_indices_type,&
                                   oscdft_constr_indices_type,&
                                   oscdft_orbital_indices_type
   USE oscdft_input,        ONLY : oscdft_input_type
   USE oscdft_context,      ONLY : oscdft_context_type

   PRIVATE
   PUBLIC oscdft_init_wfcO, orthoOwfc
   CONTAINS
      SUBROUTINE oscdft_init_wfcO(ctx)
         USE uspp_param,               ONLY : upf
         USE ions_base,                ONLY : ityp, nat
         USE wvfct,                    ONLY : npwx
         USE control_flags,            ONLY : io_level, restart, gamma_only, use_gpu
         USE noncollin_module,         ONLY : npol
         USE symm_base,                ONLY : nsym
         USE io_files,                 ONLY : restart_dir
         USE oscdft_occupations,       ONLY : oscdft_get_occupation_numbers
         USE oscdft_wavefunction_subs, ONLY : oscdft_init_wavefunctions,&
                                              oscdft_debug_print_wavefunctions
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         LOGICAL                                          :: exst, in_wfcO, in_wfcF, in_wfcS,&
                                                             test
         TYPE(oscdft_input_type),           POINTER       :: inp
         TYPE(oscdft_indices_type),         POINTER       :: idx
         TYPE(oscdft_constr_indices_type),  POINTER       :: constr
         TYPE(oscdft_orbital_indices_type), POINTER       :: orbs
         TYPE(oscdft_wavefunction_type),    POINTER       :: wfcO, wfcF, wfcS

         INTEGER :: iunO, iunS, iunF, iunat, iunsat

         inp  => ctx%inp
         idx  => ctx%idx
         constr => idx%constr
         orbs   => idx%orbs

         wfcO => ctx%wfcO
         wfcS => ctx%wfcS
         wfcF => ctx%forces%wfcO

         IF (.NOT. ctx%wfc_allocated) THEN
            iunO   = 40
            iunS   = 41
            iunF   = 42
            iunat  = 0
            iunsat = 0

            IF (inp%debug_print) THEN
               WRITE(stdout, *) ""
               WRITE(stdout, 200) "OSCDFT_INIT_WFCO"
            END IF

            CALL oscdft_init_wavefunctions(ctx, wfcO, iunO, "wfcO", .true., .false., .true.)
            CALL oscdft_init_wavefunctions(ctx, wfcS, iunS, "wfcS", .false., .true., .true.)
            CALL oscdft_init_wavefunctions(ctx, wfcF, iunF, "wfcF", .true., .true., .false.)

            IF (inp%debug_print) THEN
               CALL oscdft_debug_print_wavefunctions(wfcO, "wfcO", "S*atwfc to apply constraint")
               CALL oscdft_debug_print_wavefunctions(wfcS, "wfcS", "S*atwfc to calc ns")
               CALL oscdft_debug_print_wavefunctions(wfcF, "wfcF", "atwfc to calc forces")
            END IF

            IF (restart) THEN
               CALL oscdft_get_occupation_numbers(ctx, .false.)
            END IF
            ctx%wfc_allocated = .true.
         END IF

         CALL orthoOwfc(ctx)

         100 FORMAT("OSCDFT DEBUG: ", A, ": ", I6)
         101 FORMAT("OSCDFT DEBUG: ", A, ": ", ES14.7)
         102 FORMAT("OSCDFT DEBUG: ", A, ": ", *(I6, :, " "))
         103 FORMAT("OSCDFT DEBUG: ", A, ": ", *(ES14.7, :, " "))
         200 FORMAT("OSCDFT DEBUG: ", A)
         300 FORMAT("OSCDFT DEBUG: ", A, ": ", I8)
      END SUBROUTINE oscdft_init_wfcO

      SUBROUTINE orthoOwfc(ctx)
         USE control_flags,            ONLY : offload_type, use_gpu
         USE buffers,                  ONLY : get_buffer, save_buffer
         USE klist,                    ONLY : nks, xk, ngk, igk_k
         USE wvfct,                    ONLY : npwx
         USE uspp,                     ONLY : nkb, vkb
         USE uspp_param,               ONLY : upf
         USE becmod,                   ONLY : allocate_bec_type_acc,&
                                              deallocate_bec_type_acc,&
                                              becp, calbec, bec_type
         USE noncollin_module,         ONLY : npol
         USE basis,                    ONLY : natomwfc
         USE ions_base,                ONLY : nat, ityp
         USE symm_base,                ONLY : nsym
         USE uspp_init,                ONLY : init_us_2
         USE oscdft_wavefunction,      ONLY : check_bec_type_unallocated
         USE oscdft_wavefunction_subs, ONLY : oscdft_fill_wavefunctions,&
                                              oscdft_poolrecover_overlap,&
                                              oscdft_ortho_swfc,&
                                              oscdft_get_overlap,&
                                              oscdft_lowdin_ortho_overlap
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
         CALL start_clock("oscdft_wfcO")
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
         !$acc data create(wfcatom, swfcatom) present(vkb)

         CALL check_bec_type_unallocated(becp)
         CALL allocate_bec_type_acc(nkb, natomwfc, becp)

         DO ik=1,nks
            CALL atomic_wfc(ik, wfcatom)
            npw = ngk(ik)
            CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
            CALL calbec(offload_type, npw, vkb, wfcatom, becp)
            IF (use_gpu) THEN
               CALL s_psi_acc(npwx, npw, natomwfc, wfcatom, swfcatom)
            ELSE
               CALL s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)
            END IF

            IF (inp%orthogonalize_swfc) THEN
               CALL oscdft_ortho_swfc(npwx, npw, natomwfc, wfcatom, swfcatom, .false.)
            ELSE IF (inp%normalize_swfc) THEN
               CALL oscdft_ortho_swfc(npwx, npw, natomwfc, wfcatom, swfcatom, .true.)
            END IF

            !$acc update host(wfcatom, swfcatom)
            CALL oscdft_fill_wavefunctions(idx, wfcO, ik, wfcatom, swfcatom)
            CALL oscdft_fill_wavefunctions(idx, wfcS, ik, wfcatom, swfcatom)
            CALL oscdft_fill_wavefunctions(idx, wfcF, ik, wfcatom, swfcatom)

            IF (inp%orthogonalize_ns) THEN
               CALL oscdft_get_overlap(ctx, ik, wfcS%source, wfcatom, swfcatom)
            END IF
         END DO
         CALL deallocate_bec_type_acc(becp)
         !$acc end data
         DEALLOCATE(wfcatom, swfcatom)

         IF (inp%orthogonalize_ns) THEN
            CALL oscdft_poolrecover_overlap(idx)
            CALL oscdft_lowdin_ortho_overlap(ctx)
            ! CALL oscdft_write_overlap(ctx)
         END IF

         CALL stop_clock("oscdft_wfcO")
      END SUBROUTINE orthoOwfc
#endif
END MODULE oscdft_wfcO
