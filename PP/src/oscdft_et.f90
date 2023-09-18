PROGRAM oscdft_et
#if defined (__OSCDFT)
   USE mp,                       ONLY : mp_bcast
   USE mp_world,                 ONLY : world_comm
   USE mp_global,                ONLY : mp_startup
   USE io_global,                ONLY : stdout, ionode, ionode_id
   USE io_files,                 ONLY : tmp_dir
   USE plugin_flags,             ONLY : use_oscdft
   USE environment,              ONLY : environment_start, environment_end
   USE oscdft_et_mod,            ONLY : oscdft_file_type,&
                                        oscdft_wavefunction_type,&
                                        oscdft_read_file, get_et,&
                                        oscdft_close_file,&
                                        print_debug, oscdft_et_print_clock,&
                                        oscdft_nelup_neldw_from_input
   USE oscdft_wavefunction_subs, ONLY : oscdft_destroy_wavefunctions
   USE kinds,                    ONLY : DP
   USE klist,                    ONLY : nelup, neldw

   IMPLICIT NONE

   CHARACTER(LEN=256)             :: initial_prefix,&
                                     final_prefix,&
                                     initial_dir,&
                                     final_dir
   INTEGER                        :: ios, ierr, idx
   TYPE(oscdft_file_type)         :: initial_type, final_type
   COMPLEX(DP)                    :: idotf, fdoti, idoti, fdotf
   LOGICAL                        :: print_eigvect, print_matrix
   REAL(DP)                       :: nelup_inp, neldw_inp
   NAMELIST / oscdft_et_namelist / initial_prefix,&
                                   final_prefix,&
                                   initial_dir,&
                                   final_dir,&
                                   print_matrix,&
                                   print_eigvect,&
                                   print_debug,&
                                   nelup_inp, neldw_inp
#if defined(__MPI)
   CALL mp_startup()
#endif
   CALL environment_start("OSCDFT_ET")

   use_oscdft = .true.
   print_matrix = .false.
   print_eigvect = .false.
   nelup_inp = 0.D0
   neldw_inp = 0.D0

   ios = 0
   IF (ionode) THEN
      CALL input_from_file()
      READ(5, oscdft_et_namelist, iostat=ios)
   ENDIF

   CALL mp_bcast(ios, ionode_id, world_comm)
   IF (ios.NE.0) CALL errore("oscdft_et", "reading namelist error", abs(ios))

   CALL mp_bcast(initial_prefix, ionode_id, world_comm)
   CALL mp_bcast(final_prefix,   ionode_id, world_comm)
   CALL mp_bcast(initial_dir,    ionode_id, world_comm)
   CALL mp_bcast(final_dir,      ionode_id, world_comm)
   CALL mp_bcast(print_matrix,   ionode_id, world_comm)
   CALL mp_bcast(print_eigvect,  ionode_id, world_comm)
   CALL mp_bcast(print_debug,    ionode_id, world_comm)
   CALL mp_bcast(nelup_inp,      ionode_id, world_comm)
   CALL mp_bcast(neldw_inp,      ionode_id, world_comm)

   initial_type%prefix      = initial_prefix
   final_type%prefix        = final_prefix
   initial_type%dir         = TRIM(initial_dir)
   final_type%dir           = TRIM(final_dir)

   idx = LEN_TRIM(initial_dir)
   IF (initial_type%dir(idx:idx) /= '/') THEN
      initial_type%dir(idx+1:idx+1) = '/'
   END IF
   idx = LEN_TRIM(final_dir)
   IF (final_type%dir(idx:idx) /= '/') THEN
      final_type%dir(idx+1:idx+1) = '/'
   END IF

   CALL oscdft_read_file(initial_type, final_type, ierr)

   initial_type%inp%print_occup_matrix   = print_matrix
   initial_type%inp%print_occup_eigvects = print_eigvect
   final_type%inp%print_occup_matrix     = print_matrix
   final_type%inp%print_occup_eigvects   = print_eigvect

   CALL oscdft_nelup_neldw_from_input(nelup_inp, neldw_inp)
   CALL get_et(initial_type, final_type)

   CALL oscdft_close_file(initial_type)
   CALL oscdft_close_file(final_type)

   CALL oscdft_et_print_clock()
   CALL environment_end("OSCDFT_ET")
   CALL stop_pp
#endif
END PROGRAM oscdft_et
