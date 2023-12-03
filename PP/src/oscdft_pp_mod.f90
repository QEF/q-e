MODULE oscdft_pp_mod
#if defined (__OSCDFT)
   USE kinds,               ONLY : DP
   USE io_global,           ONLY : stdout, ionode
   USE oscdft_wavefunction, ONLY : oscdft_wavefunction_type
   USE oscdft_input,        ONLY : oscdft_input_type
   USE oscdft_indices,      ONLY : oscdft_indices_type, oscdft_orbital_indices_type
   USE oscdft_context,      ONLY : oscdft_ns_type

   PRIVATE
   PUBLIC oscdft_file_type,&
          oscdft_read_file,&
          oscdft_pp_run,&
          oscdft_close_files,&
          oscdft_pp_print_clock,&
          oscdft_copy_wavefunctions

   TYPE oscdft_file_type
      CHARACTER(LEN=256)             :: prefix, outdir
      TYPE(oscdft_wavefunction_type) :: wfc, wfcS
      TYPE(oscdft_input_type)        :: inp
      TYPE(oscdft_indices_type)      :: idx
      TYPE(oscdft_ns_type)           :: nst
   END TYPE oscdft_file_type

   CONTAINS
      SUBROUTINE oscdft_copy_wavefunctions(extension_dest, wfc)
         USE io_files,         ONLY : iunwfc, nwordwfc, restart_dir
         USE buffers,          ONLY : open_buffer, save_buffer
         USE wvfct,            ONLY : nbnd, npwx
         USE noncollin_module, ONLY : npol
         USE klist,            ONLY : nks
         USE control_flags,    ONLY : io_level
         USE wavefunctions,    ONLY : evc
         USE pw_restart_new,   ONLY : read_collected_wfc
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)                  :: extension_dest
         TYPE(oscdft_wavefunction_type), INTENT(INOUT) :: wfc
         LOGICAL                                       :: exst
         INTEGER                                       :: ik, ierr

         WRITE(stdout, 100) extension_dest
         wfc%n = nbnd
         wfc%nword = nbnd * npwx * npol
         CALL open_buffer(wfc%iun, extension_dest, wfc%nword, io_level, exst)
         ALLOCATE(wfc%wfc(npwx * npol, nbnd))
         WRITE(stdout, 101) restart_dir()
         DO ik=1,nks
            CALL read_collected_wfc(restart_dir(), ik, evc)
            IF (nks > 1) THEN
               CALL save_buffer(evc, wfc%nword, wfc%iun, ik)
            ELSE
               wfc%wfc = evc
            END IF
         END DO
         100 FORMAT("OSCDFT: Copying wfc to ", A)
         101 FORMAT("OSCDFT: Reading wfc from ", A)
      END SUBROUTINE oscdft_copy_wavefunctions

      SUBROUTINE init_wfcS(inp, idx, wfcS, iun, extensionS)
         USE basis,                    ONLY : natomwfc, swfcatom
         USE control_flags,            ONLY : gamma_only, io_level
         USE noncollin_module,         ONLY : npol
         USE wvfct,                    ONLY : npwx
         USE buffers,                  ONLY : open_buffer, get_buffer, save_buffer
         USE ions_base,                ONLY : nat, ityp
         USE uspp_param,               ONLY : upf
         USE uspp_init,                ONLY : init_us_2
         USE mp,                       ONLY : mp_barrier
         USE mp_pools,                 ONLY : inter_pool_comm, npool, my_pool_id, me_pool
         USE klist,                    ONLY : nks, xk, ngk, igk_k
         USE uspp,                     ONLY : nkb, vkb
         USE becmod,                   ONLY : allocate_bec_type,&
                                              deallocate_bec_type,&
                                              becp, calbec, bec_type
         USE symm_base,                ONLY : nsym
         USE oscdft_wavefunction_subs, ONLY : oscdft_ortho_swfc, oscdft_get_overlap,&
                                              oscdft_write_overlap,&
                                              oscdft_poolrecover_overlap,&
                                              oscdft_lowdin_ortho_overlap,&
                                              oscdft_init_wavefunctions,&
                                              oscdft_debug_print_wavefunctions,&
                                              oscdft_fill_wavefunctions
         IMPLICIT NONE
         TYPE(oscdft_input_type),        INTENT(INOUT)         :: inp
         TYPE(oscdft_indices_type),      INTENT(INOUT), TARGET :: idx
         TYPE(oscdft_wavefunction_type), INTENT(INOUT)         :: wfcS
         INTEGER,                        INTENT(IN)            :: iun
         CHARACTER(LEN=*),               INTENT(IN)            :: extensionS

         TYPE(oscdft_orbital_indices_type), POINTER            :: orbs

         INTEGER                  :: ioscdft,&
                                     n, l, npw, n1, n2, m1, m2,&
                                     counter, na, nt, iwfc, isym,&
                                     ik, iorb, m
         LOGICAL                  :: in_wfcS, test, exst
         COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)

         orbs => idx%orbs

         CALL start_clock("oscdft_wfcO")

         CALL oscdft_init_wavefunctions(inp, idx, wfcS, iun, extensionS, .false., .true., .true.)
         IF (inp%debug_print) THEN
            CALL oscdft_debug_print_wavefunctions(wfcS, "wfcS", "S*atwfc to calc ns")
         END IF

         IF (wfcS%n > 0) THEN
            ALLOCATE(wfcatom(npwx*npol,natomwfc), swfcatom(npwx*npol,natomwfc))
            CALL allocate_bec_type(nkb, natomwfc, becp)
            DO ik=1,nks
               CALL atomic_wfc(ik, wfcatom)
               npw = ngk(ik)
               CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
               CALL calbec(npw, vkb, wfcatom, becp)
               CALL s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)

               IF (inp%orthogonalize_swfc) THEN
                  CALL oscdft_ortho_swfc(npwx, npw, natomwfc, wfcatom, swfcatom, .false.)
               ELSE IF (inp%normalize_swfc) THEN
                  CALL oscdft_ortho_swfc(npwx, npw, natomwfc, wfcatom, swfcatom, .true.)
               END IF

               CALL oscdft_fill_wavefunctions(idx, wfcS, ik, wfcatom, swfcatom)

               IF (inp%orthogonalize_ns) THEN
                  CALL oscdft_get_overlap(inp, idx, ik, wfcS%source, wfcatom, swfcatom)
               END IF
            END DO
            CALL deallocate_bec_type(becp)
            DEALLOCATE(wfcatom, swfcatom)
         END IF
         IF (inp%orthogonalize_ns) THEN
            CALL oscdft_poolrecover_overlap(idx)
            CALL oscdft_lowdin_ortho_overlap(inp, idx)
            ! CALL oscdft_write_overlap(inp, idx)
         END IF

         CALL stop_clock("oscdft_wfcO")
         100 FORMAT("OSCDFT DEBUG: ", A, ": ", I6)
         102 FORMAT("OSCDFT DEBUG: ", A, ": ", *(I6, :, " "))
      END SUBROUTINE init_wfcS

      SUBROUTINE oscdft_read_file(f)
         USE io_files,         ONLY : nwordwfc, prefix, tmp_dir, wfc_dir
         USE noncollin_module, ONLY : npol
         USE wvfct,            ONLY : nbnd, npwx
         USE input_parameters, ONLY : nat
         USE lsda_mod,         ONLY : nspin
         USE oscdft_context,   ONLY : oscdft_init_indices
         USE symm_base,        ONLY : d1, d2, d3, nsym
         USE oscdft_input,     ONLY : oscdft_read_input
         IMPLICIT NONE

         TYPE(oscdft_file_type), TARGET, INTENT(INOUT) :: f
         INTEGER :: isym, row
         LOGICAL :: needwf

         CALL start_clock("oscdft_init")
         IF (ionode) THEN
            WRITE(stdout, 100) "Reading data from directory:", TRIM(f%outdir)
         END IF
         tmp_dir = f%outdir
         prefix  = f%prefix
         wfc_dir = f%outdir

         needwf = .true.
         CALL read_file_new(needwf)
         IF (.NOT.needwf) THEN
            CALL errore("oscdft_read_file", "wavefunction not collected", 1)
         END IF

         f%wfc%iun  = 42
         nwordwfc = nbnd * npwx * npol
         CALL oscdft_copy_wavefunctions("wfc", f%wfc)

         CALL oscdft_read_input(f%inp)
         CALL oscdft_init_indices(f%idx, f%inp)
         CALL d_matrix(d1, d2, d3)

         ! IF (f%inp%debug_print) THEN
         !    DO isym=1,nsym
         !       WRITE(stdout, 601)
         !       DO row=1,3
         !          WRITE(stdout, 600) "d1", isym, d1(row,:,isym)
         !       END DO
         !       WRITE(stdout, *)
         !       DO row=1,5
         !          WRITE(stdout, 600) "d2", isym, d2(row,:,isym)
         !       END DO
         !       WRITE(stdout, *)
         !       DO row=1,7
         !          WRITE(stdout, 600) "d3", isym, d3(row,:,isym)
         !       END DO
         !    END DO
         !    WRITE(stdout, 601)
         ! END IF

         CALL stop_clock("oscdft_init")

         CALL init_wfcS(f%inp, f%idx, f%wfcS, 40, "wfcS")

         100 FORMAT("OSCDFT: ", A, A)
         600 FORMAT("OSCDFT DEBUG: ", A, "(", I5, "): ", *(F8.5, " "))
         601 FORMAT("=======================================================================================")
      END SUBROUTINE oscdft_read_file

      SUBROUTINE oscdft_pp_run(f)
         USE oscdft_context,     ONLY : oscdft_alloc_nst
         USE oscdft_occupations, ONLY : oscdft_new_ns, oscdft_get_occupation_numbers
         IMPLICIT NONE

         TYPE(oscdft_file_type), INTENT(INOUT) :: f

         CALL oscdft_alloc_nst(f%nst, f%idx%max_ns_dim,&
                               f%idx%nconstr, f%inp%noscdft)
         CALL oscdft_new_ns(f%inp, f%idx, f%wfcS, f%nst, f%wfc)
         CALL oscdft_get_occupation_numbers(f%inp, f%idx, f%nst, .false.)
      END SUBROUTINE oscdft_pp_run

      SUBROUTINE oscdft_close_files(f)
         USE oscdft_wavefunction_subs, ONLY : oscdft_destroy_wavefunctions
         IMPLICIT NONE

         TYPE(oscdft_file_type), INTENT(INOUT) :: f

         CALL oscdft_destroy_wavefunctions(f%wfcS, "DELETE")
         CALL oscdft_destroy_wavefunctions(f%wfc, "DELETE")
      END SUBROUTINE oscdft_close_files

      SUBROUTINE oscdft_pp_print_clock
         IMPLICIT NONE

         WRITE(stdout, '(/,5X,"OSCDFT_PP routines")')
         CALL print_clock("oscdft_init")
         CALL print_clock("oscdft_wfcO")
         CALL print_clock("oscdft_ns")
      END SUBROUTINE oscdft_pp_print_clock
#endif
END MODULE oscdft_pp_mod
