MODULE oscdft_input
#if defined (__OSCDFT)
   USE kinds,       ONLY : DP
   USE parser,      ONLY : field_count, read_line, get_field, parse_unit
   USE io_global,   ONLY : ionode, ionode_id, stdout
   USE mp,          ONLY : mp_bcast
   USE mp_images,   ONLY : intra_image_comm
   USE oscdft_enums

   PRIVATE
   PUBLIC oscdft_input_type, oscdft_read_input

   TYPE oscdft_input_type
      LOGICAL               :: print_occup_matrix,&
                               print_occup_eigvects,&
                               has_max_multiplier,&
                               has_min_multiplier,&
                               get_ground_state_first,&
                               skip_forces,&
                               hpsi_sum_band,&
                               sum_band,&
                               rhoiter,&
                               orthogonalize_swfc,&
                               orthogonalize_ns,&
                               normalize_swfc,&
                               debug_print
      INTEGER               :: noscdft,&
                               warm_up_niter,&
                               convergence_type,& ! CONV_*
                               optimization_method,& ! OPT_*
                               array_convergence_func,& ! CONV_FUNC_*
                               iteration_type,&
                               maxiter,&
                               miniter,&
                               swapping_technique,&
                               test_exit_rho_iter,&
                               test_exit_oscdft_iter
      REAL(DP)              :: min_conv_thr,&
                               max_conv_thr,&
                               final_conv_thr,&
                               conv_thr_multiplier,&
                               min_gamma_n,&
                               max_multiplier,&
                               min_multiplier,&
                               total_energy_threshold
      LOGICAL, ALLOCATABLE  :: print_occup(:) ! ioscdft
      INTEGER, ALLOCATABLE  :: constraint_applied(:),& ! ioscdft; CONSTR_*
                               spin_index(:),& ! ioscdft
                               occup_index(:),& ! ioscdft
                               occup_index_sum(:,:),& ! ioscdft
                               start_index(:) ! ioscdft
      REAL(DP), ALLOCATABLE :: target_occup(:),& ! ioscdft
                               initial_multipliers(:),& ! ioscdft
                               debug_step(:),& ! ioscdft
                               gamma_val(:) ! ioscdft
      CHARACTER(LEN=256), ALLOCATABLE :: orbital_desc(:)
   END TYPE oscdft_input_type

   CONTAINS
      SUBROUTINE oscdft_read_input(inp, filename_inp)
         IMPLICIT NONE

         TYPE(oscdft_input_type),      INTENT(INOUT) :: inp
         CHARACTER(LEN=256), OPTIONAL, INTENT(IN)    :: filename_inp

         INTEGER, EXTERNAL :: find_free_unit
         CHARACTER(LEN=256) :: filename

         INTEGER :: iun
         LOGICAL :: ext

         IF (PRESENT(filename_inp)) THEN
            filename = TRIM(filename_inp)
         ELSE
            filename = 'oscdft.in'
         END IF
         iun = find_free_unit()
         INQUIRE(file=TRIM(filename), exist=ext)
         IF (.NOT.ext) CALL errore("read_oscdft", "missing " // TRIM(filename) // " input file", 1)
         OPEN(unit=iun, file=TRIM(filename), status="old")

         CALL read_namelist(inp, iun)
         CALL alloc_inp(inp)
         CALL read_cards(inp, iun)

         CLOSE(iun)
      END SUBROUTINE oscdft_read_input

      SUBROUTINE capitalize(string)
         IMPLICIT NONE

         CHARACTER(LEN=*), INTENT(INOUT) :: string
         CHARACTER(LEN=1), EXTERNAL      :: capital
         INTEGER                         :: idx

         DO idx=1,LEN(string)
            string(idx:idx) = capital(string(idx:idx))
         ENDDO
      END SUBROUTINE capitalize

      SUBROUTINE read_namelist(inp, iun)
         IMPLICIT NONE

         TYPE(oscdft_input_type), INTENT(INOUT) :: inp
         INTEGER,                 INTENT(IN)    :: iun

         INTEGER :: ios, idx, error_type

         LOGICAL            :: print_occupation_matrix,&
                               print_occupation_eigenvectors,&
                               has_max_multiplier,&
                               has_min_multiplier,&
                               get_ground_state_first,&
                               skip_forces,&
                               sum_band,&
                               rhoiter,&
                               orthogonalize_swfc,&
                               orthogonalize_ns,&
                               normalize_swfc,&
                               debug_print
         INTEGER            :: n_oscdft,&
                               warm_up_niter,&
                               iteration_type,&
                               maxiter,&
                               miniter,&
                               test_exit_rho_iter,&
                               test_exit_oscdft_iter
         CHARACTER(LEN=256) :: convergence_type,&
                               optimization_method,&
                               array_convergence_func,&
                               swapping_technique
         REAL(DP)           :: min_conv_thr,&
                               max_conv_thr,&
                               final_conv_thr,&
                               conv_thr_multiplier,&
                               min_gamma_n,&
                               min_multiplier,&
                               max_multiplier,&
                               total_energy_threshold
         NAMELIST / oscdft / n_oscdft,&
                             warm_up_niter,&
                             convergence_type,&
                             iteration_type,&
                             optimization_method,&
                             array_convergence_func,&
                             min_conv_thr,&
                             max_conv_thr,&
                             final_conv_thr,&
                             conv_thr_multiplier,&
                             print_occupation_matrix,&
                             print_occupation_eigenvectors,&
                             min_gamma_n,&
                             has_max_multiplier,&
                             max_multiplier,&
                             has_min_multiplier,&
                             min_multiplier,&
                             get_ground_state_first,&
                             skip_forces,&
                             maxiter,&
                             miniter,&
                             swapping_technique,&
                             sum_band,&
                             total_energy_threshold,&
                             rhoiter,&
                             debug_print,&
                             test_exit_rho_iter,&
                             test_exit_oscdft_iter,&
                             orthogonalize_swfc,&
                             orthogonalize_ns,&
                             normalize_swfc


         ! defaults
         n_oscdft                        = -1
         iteration_type                  = -1
         warm_up_niter                   = 0
         convergence_type                = "GRADIENT"
         optimization_method             = "GRADIENT DESCENT"
         array_convergence_func          = "MAXVAL"
         min_conv_thr                    = 1.D-3
         max_conv_thr                    = 1.D-1
         final_conv_thr                  = -1.D0
         conv_thr_multiplier             = 0.5D0
         print_occupation_matrix         = .false.
         print_occupation_eigenvectors   = .false.
         min_gamma_n                     = 1.D0
         has_max_multiplier              = .false.
         max_multiplier                  = 0.D0
         has_min_multiplier              = .false.
         min_multiplier                  = 0.D0
         get_ground_state_first          = .FALSE.
         skip_forces                     = .false.
         maxiter                         = 99999
         miniter                         = 0
         swapping_technique              = "NONE"
         sum_band                        = .false.
         total_energy_threshold          = 0.D0
         rhoiter                         = .false.
         test_exit_rho_iter              = -1
         test_exit_oscdft_iter           = -1
         orthogonalize_swfc              = .false.
         orthogonalize_ns                = .false.
         normalize_swfc                  = .false.

         ios = 0
         error_type = 0
         IF (ionode) THEN
            READ(iun, oscdft, iostat=ios)
            IF (ios == 0) THEN
               inp%print_occup_matrix     = print_occupation_matrix
               inp%print_occup_eigvects   = print_occupation_eigenvectors
               inp%noscdft                = n_oscdft
               inp%warm_up_niter          = warm_up_niter
               inp%min_conv_thr           = min_conv_thr
               inp%max_conv_thr           = max_conv_thr
               inp%final_conv_thr         = final_conv_thr
               inp%conv_thr_multiplier    = conv_thr_multiplier
               inp%min_gamma_n            = min_gamma_n
               inp%has_max_multiplier     = has_max_multiplier
               inp%max_multiplier         = max_multiplier
               inp%has_min_multiplier     = has_min_multiplier
               inp%min_multiplier         = min_multiplier
               inp%iteration_type         = iteration_type
               inp%get_ground_state_first = get_ground_state_first
               inp%skip_forces            = skip_forces
               inp%maxiter                = maxiter
               inp%miniter                = miniter
               inp%sum_band               = sum_band
               inp%rhoiter                = rhoiter
               inp%total_energy_threshold = total_energy_threshold
               inp%debug_print            = debug_print
               inp%test_exit_rho_iter     = test_exit_rho_iter
               inp%test_exit_oscdft_iter  = test_exit_oscdft_iter
               inp%orthogonalize_swfc     = orthogonalize_swfc
               inp%orthogonalize_ns       = orthogonalize_ns
               inp%normalize_swfc         = normalize_swfc

               CALL capitalize(convergence_type)
               SELECT CASE (TRIM(convergence_type))
                  CASE ("0", "MULTIPLIERS")
                     inp%convergence_type = CONV_MULTIPLIERS
                  CASE ("1", "GRADIENT")
                     inp%convergence_type = CONV_GRADIENT
                  CASE ("2", "ENERGY")
                     inp%convergence_type = CONV_ENERGY
                  CASE ("-1", "ALWAYS FALSE")
                     inp%convergence_type = CONV_ALWAYS_FALSE
                  CASE ("-2", "ALWAYS TRUE")
                     inp%convergence_type = CONV_ALWAYS_TRUE
                  CASE DEFAULT
                     error_type = 1
               END SELECT

               CALL capitalize(optimization_method)
               SELECT CASE (TRIM(optimization_method))
                  CASE ("0", "CONSTANT GRADIENT DESCENT", "GRADIENT DESCENT")
                     inp%optimization_method = OPT_GRADIENT_DESCENT
                  CASE ("1", "GRADIENT DESCENT2")
                     inp%optimization_method = OPT_GRADIENT_DESCENT2
                  CASE DEFAULT
                     error_type = 2
               END SELECT

               CALL capitalize(array_convergence_func)
               SELECT CASE (TRIM(array_convergence_func))
                  CASE ("0", "MAXVAL")
                     inp%array_convergence_func = CONV_FUNC_MAXVAL
                  CASE ("1", "NORM", "2-NORM")
                     inp%array_convergence_func = CONV_FUNC_NORM
                  CASE ("2", "NORM AVERAGE", "RMS")
                     inp%array_convergence_func = CONV_FUNC_NORM_AVERAGE
                  CASE DEFAULT
                     error_type = 3
               END SELECT

               CALL capitalize(swapping_technique)
               SELECT CASE (TRIM(swapping_technique))
                  CASE ("0", "NONE")
                     inp%swapping_technique = OSCDFT_NONE
                  CASE ("2", "PERMUTE")
                     inp%swapping_technique = OSCDFT_PERMUTE
                  CASE DEFAULT
                     error_type = 4
               END SELECT

               IF (inp%noscdft <= 0) error_type = 5
               IF (inp%iteration_type /= 0 .AND. inp%iteration_type /= 1) THEN
                  error_type = 6
               END IF
            END IF
         END IF
         CALL mp_bcast(ios,        ionode_id, intra_image_comm)
         CALL mp_bcast(error_type, ionode_id, intra_image_comm)

         IF (ios/=0) CALL errore("oscdft_read_namelist", "oscdft.in IO error", ABS(ios))
         SELECT CASE (error_type)
            CASE (1)
               CALL errore("oscdft_read_namelist", "invalid convergence_type", 1)
            CASE (2)
               CALL errore("oscdft_read_namelist", "invalid optimization_method", 1)
            CASE (3)
               CALL errore("oscdft read_namelist", "invalid array_convergence_func", 1)
            CASE (4)
               CALL errore("oscdft_read_namelist", "invalid swapping_technique", 1)
            CASE (5)
               CALL errore("oscdft_read_namelist", "n_oscdft <= 0", 1)
            CASE (6)
               CALL errore("oscdft_read_namelist", "iteration_type invalid", 1)
         END SELECT

         CALL bcast_inp(inp)

         IF (inp%noscdft < 0) CALL errore("oscdft_read_namelist", "n_oscdft must be >= 0", 1)
      END SUBROUTINE read_namelist

      SUBROUTINE bcast_inp(inp)
         IMPLICIT NONE

         TYPE(oscdft_input_type), INTENT(INOUT) :: inp

         CALL mp_bcast(inp%print_occup_matrix,     ionode_id, intra_image_comm)
         CALL mp_bcast(inp%print_occup_eigvects,   ionode_id, intra_image_comm)
         CALL mp_bcast(inp%noscdft,                ionode_id, intra_image_comm)
         CALL mp_bcast(inp%warm_up_niter,          ionode_id, intra_image_comm)
         CALL mp_bcast(inp%convergence_type,       ionode_id, intra_image_comm)
         CALL mp_bcast(inp%optimization_method,    ionode_id, intra_image_comm)
         CALL mp_bcast(inp%array_convergence_func, ionode_id, intra_image_comm)
         CALL mp_bcast(inp%min_conv_thr,           ionode_id, intra_image_comm)
         CALL mp_bcast(inp%max_conv_thr,           ionode_id, intra_image_comm)
         CALL mp_bcast(inp%final_conv_thr,         ionode_id, intra_image_comm)
         CALL mp_bcast(inp%conv_thr_multiplier,    ionode_id, intra_image_comm)
         CALL mp_bcast(inp%min_gamma_n,            ionode_id, intra_image_comm)
         CALL mp_bcast(inp%has_max_multiplier,     ionode_id, intra_image_comm)
         CALL mp_bcast(inp%max_multiplier,         ionode_id, intra_image_comm)
         CALL mp_bcast(inp%has_min_multiplier,     ionode_id, intra_image_comm)
         CALL mp_bcast(inp%min_multiplier,         ionode_id, intra_image_comm)
         CALL mp_bcast(inp%iteration_type,         ionode_id, intra_image_comm)
         CALL mp_bcast(inp%get_ground_state_first, ionode_id, intra_image_comm)
         CALL mp_bcast(inp%skip_forces,            ionode_id, intra_image_comm)
         CALL mp_bcast(inp%maxiter,                ionode_id, intra_image_comm)
         CALL mp_bcast(inp%miniter,                ionode_id, intra_image_comm)
         CALL mp_bcast(inp%swapping_technique,     ionode_id, intra_image_comm)
         CALL mp_bcast(inp%sum_band,               ionode_id, intra_image_comm)
         CALL mp_bcast(inp%rhoiter,                ionode_id, intra_image_comm)
         CALL mp_bcast(inp%total_energy_threshold, ionode_id, intra_image_comm)
         CALL mp_bcast(inp%debug_print,            ionode_id, intra_image_comm)
         CALL mp_bcast(inp%test_exit_rho_iter,     ionode_id, intra_image_comm)
         CALL mp_bcast(inp%test_exit_oscdft_iter,  ionode_id, intra_image_comm)
         CALL mp_bcast(inp%orthogonalize_swfc,     ionode_id, intra_image_comm)
         CALL mp_bcast(inp%orthogonalize_ns,       ionode_id, intra_image_comm)
         CALL mp_bcast(inp%normalize_swfc,         ionode_id, intra_image_comm)
      END SUBROUTINE bcast_inp

      SUBROUTINE alloc_inp(inp)
         IMPLICIT NONE

         TYPE(oscdft_input_type), INTENT(INOUT) :: inp

         ALLOCATE(inp%constraint_applied(inp%noscdft),&
            inp%spin_index(inp%noscdft),&
            inp%occup_index(inp%noscdft),&
            inp%occup_index_sum(20,inp%noscdft),&
            inp%target_occup(inp%noscdft),&
            inp%initial_multipliers(inp%noscdft),&
            inp%print_occup(inp%noscdft),&
            inp%debug_step(inp%noscdft),&
            inp%start_index(inp%noscdft),&
            inp%orbital_desc(inp%noscdft),&
            inp%gamma_val(inp%noscdft))

         inp%constraint_applied   = 0
         inp%spin_index           = 0
         inp%occup_index          = 0
         inp%occup_index_sum(:,:) = 0
         inp%target_occup         = 0.D0
         inp%initial_multipliers  = 0.D0
         inp%print_occup          = .false.
         inp%debug_step           = 0.D0
         inp%start_index          = 1
         inp%gamma_val(:)         = 1.D0
      END SUBROUTINE alloc_inp

      SUBROUTINE read_cards(inp, unit)
         IMPLICIT NONE

         TYPE(oscdft_input_type), INTENT(INOUT) :: inp
         INTEGER,                 INTENT(IN)    :: unit
         CHARACTER(LEN=256)                     :: input_line
         CHARACTER(LEN=80)                      :: card
         LOGICAL                                :: tend,&
                                                   done_target,&
                                                   done_gamma
         INTEGER                                :: i, n_total, ierr

         parse_unit = unit
         100 CALL read_line(input_line, end_of_file=tend)
         IF (tend) GOTO 120
         ! ignore commands
         IF ((input_line.EQ." ").OR.(input_line(1:1).EQ."#").OR.(input_line(1:1).EQ."!")) GOTO 100
         READ(input_line, *) card
         CALL capitalize(input_line)

         done_target = .false.
         done_gamma  = .false.

         SELECT CASE (TRIM(card))
            CASE ("TARGET_OCCUPATION_NUMBER", "TARGET OCCUPATION NUMBER", "TARGET_OCCUPATION_NUMBERS", "TARGET OCCUPATION NUMBERS")
               IF (done_target) CALL errore(&
                  "card_target_occupation_numbers",&
                  "TARGET_OCCUPATION_NUMBERS two occurences found",&
                  1)
               CALL card_target(inp, input_line)
               done_target = .true.
            CASE ("GAMMA_VAL", "GAMMA VAL")
               IF (inp%optimization_method == OPT_GRADIENT_DESCENT2) THEN
                  IF (done_gamma) CALL errore(&
                     "card_target_occupation_numbers",&
                     "TARGET_OCCUPATION_NUMBERS two occurences found",&
                     1)
                  CALL card_gamma_val(inp, input_line)
                  done_gamma = .true.
               ELSE
                  CALL errore(&
                     "oscdft_read_cards",&
                     "GAMMA_VAL needs optimization_methods == 'GRADIENT DESCENT2'",&
                     1)
               END IF
            CASE DEFAULT
               IF (ionode) WRITE(stdout, '(A)') 'Warning: card '//TRIM(input_line)//' ignored'
         END SELECT

         GOTO 100
         120 CONTINUE

         CALL mp_bcast(inp%constraint_applied,  ionode_id, intra_image_comm)
         CALL mp_bcast(inp%spin_index,          ionode_id, intra_image_comm)
         CALL mp_bcast(inp%occup_index,         ionode_id, intra_image_comm)
         CALL mp_bcast(inp%target_occup,        ionode_id, intra_image_comm)
         CALL mp_bcast(inp%initial_multipliers, ionode_id, intra_image_comm)
         CALL mp_bcast(inp%print_occup,         ionode_id, intra_image_comm)
         CALL mp_bcast(inp%debug_step,          ionode_id, intra_image_comm)
         CALL mp_bcast(inp%start_index,         ionode_id, intra_image_comm)
         CALL mp_bcast(inp%occup_index_sum,     ionode_id, intra_image_comm)
         CALL mp_bcast(inp%gamma_val,           ionode_id, intra_image_comm)
         CALL mp_bcast(inp%orbital_desc,        ionode_id, intra_image_comm)
      END SUBROUTINE read_cards

      SUBROUTINE card_target(inp, input_line)
         USE clib_wrappers, ONLY : feval_infix
         IMPLICIT NONE

         TYPE(oscdft_input_type), INTENT(INOUT) :: inp
         CHARACTER(LEN=256),      INTENT(INOUT) :: input_line
         CHARACTER(LEN=256)                     :: field_str
         LOGICAL                                :: tend
         INTEGER                                :: ioscdft, ie, ierr, nfield,&
                                                   print_ifield, ifield, isum

         DO ioscdft=1,inp%noscdft
            CALL read_line(input_line, end_of_file=tend)
            IF (tend) CALL errore("card_target_occupation_numbers",&
                                  "EOF while reading TARGET_OCCUAPTION_NUMBERS",&
                                  ioscdft)

            CALL field_count(nfield, input_line)

            ! constraint_applied
            ifield = 1
            CALL get_field(ifield, field_str, input_line)
            ifield = ifield + 1
            CALL capitalize(field_str)
            SELECT CASE (TRIM(field_str))
               CASE ("FALSE", ".FALSE.", "0", "F")
                  inp%constraint_applied(ioscdft) = CONSTR_FALSE
               CASE ("TRUE", ".TRUE.", "1", "T")
                  inp%constraint_applied(ioscdft) = CONSTR_TRUE
               CASE ("LE")
                  inp%constraint_applied(ioscdft) = CONSTR_LE
               CASE ("GE")
                  inp%constraint_applied(ioscdft) = CONSTR_GE
               CASE ("LE2")
                  inp%constraint_applied(ioscdft) = CONSTR_LE2
               CASE ("GE2")
                  inp%constraint_applied(ioscdft) = CONSTR_GE2
               CASE ("LE3")
                  inp%constraint_applied(ioscdft) = CONSTR_LE3
               CASE ("GE3")
                  inp%constraint_applied(ioscdft) = CONSTR_GE3
               CASE ("DEBUG", "D", "D1")
                  inp%constraint_applied(ioscdft) = CONSTR_D1
               CASE DEFAULT
                  CALL errore("card_target_occupation_numbers",&
                              "constraint_applied: unrecognized field", ioscdft)
            END SELECT
            IF (inp%constraint_applied(ioscdft)/=CONSTR_FALSE.AND.nfield<5) THEN
               CALL errore("card_target_occupation_numbers",&
                           "not enough fields for constraint_applied = .TRUE.",&
                           ioscdft)
            END IF

            ! spin_index
            CALL get_field(ifield, field_str, input_line)
            ifield = ifield + 1
            CALL capitalize(field_str)
            IF (TRIM(field_str) == "UP") THEN
               inp%spin_index(ioscdft) = 1
            ELSE IF (TRIM(field_str) == "DW" .OR. TRIM(field_str) == "DOWN") THEN
               inp%spin_index(ioscdft) = 2
            ELSE
               inp%spin_index(ioscdft) = INT(feval_infix(ierr, field_str))
            END IF
            IF ((inp%spin_index(ioscdft)/=1).AND.(inp%spin_index(ioscdft)/=2)) THEN
               CALL errore("card_target_occupation_numbers",&
                           "spin_index invalid; 1 for SPIN UP, 2 for SPIN DOWN",&
                           ioscdft)
            END IF

            ! something like
            ! [T/F] [spin] 3(2s2p)1(1s) ...
            IF (nfield.GE.ifield) THEN
               CALL get_field(ifield, field_str, input_line)
               ifield = ifield + 1
               field_str = TRIM(field_str)
               CALL capitalize(field_str)
               inp%orbital_desc(ioscdft) = field_str
            ELSE
               IF (inp%constraint_applied(ioscdft)/=CONSTR_FALSE) THEN
                  CALL errore("card_target_occupation_numbers",&
                              "incomplete cards",&
                              ioscdft)
               END IF
            END IF

            IF (inp%constraint_applied(ioscdft)/=CONSTR_FALSE) THEN
               CALL get_field(ifield, field_str, input_line)
               ifield = ifield + 1
               CALL capitalize(field_str)
               SELECT CASE (TRIM(field_str))
                  CASE ("TRACE")
                     inp%occup_index(ioscdft) = OCCUP_TRACE
                  CASE ("SUM")
                     inp%occup_index(ioscdft) = OCCUP_SUM
                  CASE DEFAULT
                     inp%occup_index(ioscdft) = INT(feval_infix(ierr, field_str))
               END SELECT

               IF (inp%occup_index(ioscdft).EQ.OCCUP_SUM) THEN
                  CALL get_field(ifield, field_str, input_line)
                  ifield = ifield + 1
                  inp%occup_index_sum(1,ioscdft) = INT(feval_infix(ierr, field_str))

                  DO isum=1,inp%occup_index_sum(1,ioscdft)
                     CALL get_field(ifield, field_str, input_line)
                     ifield = ifield + 1
                     inp%occup_index_sum(isum+1,ioscdft) = INT(feval_infix(ierr, field_str))
                  END DO
               END IF

               ! target
               CALL get_field(ifield, field_str, input_line)
               ifield = ifield + 1
               inp%target_occup(ioscdft) = feval_infix(ierr, field_str)

               CALL get_field(ifield, field_str, input_line)
               ifield = ifield + 1
               inp%initial_multipliers(ioscdft) = feval_infix(ierr, field_str)
            END IF

            print_ifield = ifield
            IF (nfield>=print_ifield) THEN
               CALL get_field(print_ifield, field_str, input_line)
               CALL capitalize(field_str)
               SELECT CASE (TRIM(field_str))
                  CASE ("TRUE", ".TRUE.", "1", "T")
                     inp%print_occup(ioscdft) = .true.
                  CASE ("FALSE", ".FALSE.", "0", "F")
                     inp%print_occup(ioscdft) = .false.
                  CASE DEFAULT
                     inp%print_occup(ioscdft) = .true.
               END SELECT
               ifield = ifield + 1
            ELSE
               inp%print_occup(ioscdft) = .true.
            END IF

            IF (nfield >= ifield) THEN
               CALL get_field(ifield, field_str, input_line)
               inp%start_index(ioscdft) = INT(feval_infix(ierr, field_str))
            ELSE
               inp%start_index(ioscdft) = 1
            ENDIF

            IF ((nfield>=ifield).AND.(inp%constraint_applied(ioscdft)>=CONSTR_D1)) THEN
               CALL get_field(ifield, field_str, input_line)
               ifield = ifield + 1
               inp%debug_step(ioscdft) = feval_infix(ierr, field_str)
            ENDIF
         END DO
      END SUBROUTINE card_target

      SUBROUTINE card_gamma_val(inp, input_line)
         USE clib_wrappers, ONLY : feval_infix
         IMPLICIT NONE

         TYPE(oscdft_input_type), INTENT(INOUT) :: inp
         CHARACTER(LEN=256),      INTENT(INOUT) :: input_line
         CHARACTER(LEN=256)                     :: field_str
         LOGICAL                                :: tend
         INTEGER                                :: ioscdft, nfield, ierr

         DO ioscdft=1,inp%noscdft
            CALL read_line(input_line, end_of_file=tend)
            IF (tend) CALL errore("card_gamma_val",&
                                  "EOF while reading GAMMA_VAL",&
                                  ioscdft)

            CALL field_count(nfield, input_line)
            IF (nfield/=1) CALL errore("card_gamma_val", "expected 1 field", ioscdft)


            CALL get_field(1, field_str, input_line)
            inp%gamma_val(ioscdft) = feval_infix(ierr, field_str)
         END DO
      END SUBROUTINE card_gamma_val
#endif
   END MODULE oscdft_input
