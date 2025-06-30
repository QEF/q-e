MODULE oscdft_base
#if defined (__OSCDFT)
   USE oscdft_context, ONLY : oscdft_context_type
   USE oscdft_input,   ONLY : oscdft_input_type

   PRIVATE

   TYPE(oscdft_context_type), SAVE :: oscdft_ctx
   PUBLIC :: oscdft_ctx, print_oscdft_clocks

   CONTAINS
      SUBROUTINE print_oscdft_clocks(ctx) ! TODO: GPU
         USE io_global,     ONLY : stdout
         USE control_flags, ONLY : use_gpu
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER       :: inp

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         ! TODO: label limited to 12 chars
         WRITE(stdout, '(/,5X,"OSCDFT routines")')
         !                 oscdft_xxxxx  12 chars
         CALL print_clock("oscdft_init")
         CALL print_clock("oscdft_wfcO")
         CALL print_clock("oscdft_ns")
         CALL print_clock("oscdft_hdiag")
         CALL print_clock("oscdft_hpsi")
         CALL print_clock("oscdft_force")
      END SUBROUTINE print_oscdft_clocks
#endif
END MODULE oscdft_base
