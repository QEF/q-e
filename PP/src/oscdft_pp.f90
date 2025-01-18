PROGRAM oscdft_pp
#if defined (__OSCDFT)
   USE mp,                  ONLY : mp_bcast
   USE mp_world,            ONLY : world_comm
   USE mp_global,           ONLY : mp_startup
   USE io_global,           ONLY : stdout, ionode, ionode_id
   USE environment,         ONLY : environment_start, environment_end
   USE oscdft_pp_mod,       ONLY : oscdft_file_type,&
                                   oscdft_read_file,&
                                   oscdft_pp_run,&
                                   oscdft_close_files,&
                                   oscdft_pp_print_clock
   IMPLICIT NONE

   CHARACTER(LEN=256) :: prefix, outdir
   LOGICAL            :: print_matrix, print_eigvect

   INTEGER :: ios, ierr
   TYPE(oscdft_file_type) :: f

   NAMELIST / oscdft_pp_namelist / prefix, outdir

#if defined (__MPI)
   CALL mp_startup
#endif

   CALL environment_start("OSCDFT_PP")

   print_matrix = .false.
   print_eigvect = .false.

   ios = 0
   IF (ionode) THEN
      CALL input_from_file
      READ(5, oscdft_pp_namelist, iostat=ios)
   END IF

   CALL mp_bcast(ios, ionode_id, world_comm)
   IF (ios.NE.0) CALL errore("oscdft_pp", "reading namelist error", abs(ios))

   CALL mp_bcast(prefix,        ionode_id, world_comm)
   CALL mp_bcast(outdir,        ionode_id, world_comm)
   CALL mp_bcast(print_matrix,  ionode_id, world_comm)
   CALL mp_bcast(print_eigvect, ionode_id, world_comm)

   f%prefix = TRIM(prefix)
   f%outdir = TRIM(outdir)

   CALL oscdft_read_file(f)
   CALL oscdft_pp_run(f)
   CALL oscdft_close_files(f)

   CALL oscdft_pp_print_clock()
   CALL environment_end("OSCDFT_PP")
   CALL stop_pp
#endif
END PROGRAM oscdft_pp
