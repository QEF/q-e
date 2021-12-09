!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE ioneb()
  !-----------------------------------------------------------------------------
  !
  ! ... Copy neb-specific variables from path_input_variables into modules,
  ! ... Variables that have the same name in input file and in the modules
  ! ... are renamed (the logic is the same as in routine "iosys")
  ! ... This is done so that it is possible to use a different input parser
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : autoev, eV_to_kelvin
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : tmp_dir
  USE fcp_module,    ONLY : fcp_check
  !
  USE path_variables, ONLY : lsteep_des, lquick_min, &
                             lbroyden, lbroyden2, llangevin, &
                             lneb, lsmd, restart
  ! renamed variables
  USE path_variables, ONLY : nstep_path_      => nstep_path, &    
                             ds_              => ds, &
                             use_masses_      => use_masses, &
                             CI_scheme_       => CI_scheme, &
                             fixed_tan_       => fixed_tan, &
                             use_freezing_    => use_freezing, &
                             k_max_           => k_max, &
                             k_min_           => k_min, &
                             num_of_images_   => num_of_images, &
                             first_last_opt_  => first_last_opt, &
                             temp_req_        => temp_req, &
                             path_thr_        => path_thr
  !
  USE fcp_variables, ONLY : lfcp_linmin, lfcp_newton, lfcp_coupled
  ! renamed variables
  USE fcp_variables, ONLY : lfcp_ => lfcp, &
                            fcp_mu_ => fcp_mu, &
                            fcp_thr_ => fcp_thr, &
                            fcp_ndiis_ => fcp_ndiis, &
                            fcp_rdiis_ => fcp_rdiis, &
                            fcp_max_volt_ => fcp_max_volt
  !
  USE path_input_parameters_module, ONLY : restart_mode, nstep_path,   &
                               string_method, num_of_images, path_thr, &
                               CI_scheme, opt_scheme, use_masses,      &
                               first_last_opt, temp_req, k_max, k_min, &
                               ds, use_freezing, fixed_tan,            &
                               lfcp, fcp_mu, fcp_thr, fcp_scheme,      &
                               fcp_ndiis, fcp_rdiis, fcp_max_volt
  !
  IMPLICIT NONE
  !
  INTEGER  :: ia, image, nt
  REAL(DP) :: theta, phi
  !
  !
  SELECT CASE(trim( string_method ))
     !
  CASE( 'neb' )
     !
     lneb  = .true.
     !
  CASE( 'smd' )
     !
     lsmd  = .true.
     !
  CASE DEFAULT
     !
     CALL errore( 'ioneb', 'string_method ' // &
                & trim( string_method ) // ' not implemented', 1 )
     !
  END SELECT
  !
  SELECT CASE( trim( restart_mode ) )
  CASE( 'from_scratch' )
     !
     restart        = .false.
     !
  CASE( 'restart' )
     !
     IF ( lneb .or. lsmd ) THEN
        !
        ! ... "path" specific
        !
        restart = .true.
        !
     ENDIF
     !
  CASE DEFAULT
     !
     CALL errore( 'ioneb', &
                & 'unknown restart_mode ' // trim( restart_mode ), 1 )
     !
  END SELECT
  !
!
! check to move after call to iosys in PW
!
!     IF( io_level < 0) CALL errore ( 'ioneb', &
!                       'NEB, SMD do not work with "disk_io" set to "none"', 1)
     !
     !
     IF ( num_of_images < 2 ) &
        CALL errore( 'ioneb', 'string_method=' // trim( string_method ) // &
                   & ': num_of_images must be at least 2', 1 )
     !
     IF ( ( CI_scheme /= "no-CI"  ) .and. &
          ( CI_scheme /= "auto"   ) .and. &
          ( CI_scheme /= "manual" ) ) THEN
        !
        CALL errore( 'ioneb', 'string_method=' // trim( string_method ) // &
                   & ': unknown CI_scheme', 1 )
        !
     ENDIF
     !
     ! ... initialization of logical variables
     !
     lsteep_des  = .false.
     lquick_min  = .false.
     lbroyden    = .false.
     lbroyden2   = .false.
     !
     SELECT CASE( opt_scheme )
     CASE( "sd" )
        !
        lsteep_des = .true.
        !
     CASE( "quick-min" )
        !
        lquick_min = .true.
        !
     CASE( "broyden" )
        !
        lbroyden     = .true.
        !
     CASE( "broyden2" )
        !
        lbroyden2    = .true.
        !
     CASE( "langevin" )
        !
        llangevin = .true.
        !
        IF ( lneb ) &
           CALL errore( 'iosys','string_method=' // trim( string_method ) // &
                      & ': langevin dynamics not implemented', 1 )
        !
        temp_req = temp_req / ( eV_to_kelvin * autoev )
        !
        IF ( temp_req <= 0.D0 ) &
           CALL errore( 'iosys','string_method=' // trim( string_method ) // &
                      & ': tepm_req has not been set', 1 )
        !
        IF ( use_freezing ) &
           WRITE( UNIT = stdout, &
                  FMT = '(5X,"warning: freezing cannot be used in langevin")' )
        !
        use_freezing = .false.
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys','string_method=' // trim( string_method ) // &
                   & ': unknown opt_scheme', 1 )
        !
     END SELECT
     !
  !
  ! ... "path"-optimization variables
  !
  nstep_path_     = nstep_path
  ds_             = ds
  num_of_images_  = num_of_images
  first_last_opt_ = first_last_opt
  use_masses_     = use_masses
  use_freezing_   = use_freezing
  temp_req_       = temp_req
  path_thr_       = path_thr
  CI_scheme_      = CI_scheme
  k_max_          = k_max
  k_min_          = k_min
  fixed_tan_      = fixed_tan
  !
  ! ... resolve fcp_scheme
  !
  lfcp_linmin  = .FALSE.
  lfcp_newton  = .FALSE.
  lfcp_coupled = .FALSE.
  !
  SELECT CASE( fcp_scheme )
  CASE( "lm", "line-min", "line-minimization", "line-minimisation" )
     !
     fcp_scheme = "lm"
     lfcp_linmin = .true.
     !
  CASE( "newton" )
     !
     fcp_scheme = "newton"
     lfcp_newton = .TRUE.
     !
  CASE( "couple", "coupled", "coupling" )
     !
     fcp_scheme = "coupled"
     lfcp_coupled = .TRUE.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys','string_method=' // trim( string_method ) // &
                & ': unknown fcp_scheme', 1 )
     !
  END SELECT
  !
  ! ... "FCP"-optimization variables
  !
  lfcp_         = lfcp
  fcp_mu_       = fcp_mu / autoev
  fcp_thr_      = fcp_thr
  fcp_ndiis_    = fcp_ndiis
  fcp_rdiis_    = fcp_rdiis
  fcp_max_volt_ = fcp_max_volt / autoev
  !
  IF ( lfcp_ ) CALL fcp_check( .TRUE. )
  !
  CALL verify_neb_tmpdir( tmp_dir )
  !
  RETURN
  !
END SUBROUTINE ioneb
!
!-----------------------------------------------------------------------
SUBROUTINE verify_neb_tmpdir( tmp_dir )
  !-----------------------------------------------------------------------
  !
  USE clib_wrappers,    ONLY : f_mkdir
  USE path_input_parameters_module, ONLY : restart_mode
  USE io_files,         ONLY : prefix, check_tempdir, delete_if_present
  USE path_variables,   ONLY : num_of_images
  USE mp_world,         ONLY : world_comm, mpime, nproc
  USE io_global,        ONLY : meta_ionode
  USE mp,               ONLY : mp_barrier
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(inout) :: tmp_dir
  !
  INTEGER             :: ios, image, proc, nofi
  LOGICAL             :: exst, parallelfs
  CHARACTER (len=256) :: file_path, filename
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  !
  file_path = trim( tmp_dir ) // trim( prefix )
  !
  !
  IF ( restart_mode == 'from_scratch' ) THEN
     !
     ! ... let us try to create the scratch directory
     !
     CALL check_tempdir ( tmp_dir, exst, parallelfs )
     !
  ENDIF
  !
  !
  ! ... if starting from scratch all temporary files are removed
  ! ... from tmp_dir ( only by the master node )
  !
  IF ( meta_ionode ) THEN
     !
     ! ... files needed by parallelization among images are removed
     !
     CALL delete_if_present( trim( file_path ) // '.newimage' )
     !
     ! ... file containing the broyden's history
     !
     IF ( restart_mode == 'from_scratch' ) THEN
        !
        CALL delete_if_present( trim( file_path ) // '.broyden' )
        !
     ENDIF
     !
  ENDIF ! end if ionode
  !
  nofi = num_of_images
  !
  DO image = 1, nofi
     !
     file_path = trim( tmp_dir ) // trim( prefix ) //"_" // &
                 trim( int_to_char( image ) ) // '/'
     !
     CALL check_tempdir ( file_path, exst, parallelfs )
     !
     ! ... if starting from scratch all temporary files are removed
     ! ... from tmp_dir ( by all the cpus in sequence )
     !
     IF ( restart_mode == 'from_scratch' ) THEN
        !
        DO proc = 0, nproc - 1
           !
           IF ( proc == mpime ) THEN
              !
              ! ... extrapolation file is removed
              !
              CALL delete_if_present( trim( file_path ) // &
                                    & trim( prefix ) // '.update' )
              !
              ! ... standard output of the self-consistency is removed
              !
              CALL delete_if_present( trim( file_path ) // 'PW.out' )
              !
           ENDIF
           !
           CALL mp_barrier( world_comm )
           !
        ENDDO
        !
     ENDIF ! end restart_mode
     !
  ENDDO ! end do image
  !
  RETURN
  !
END SUBROUTINE verify_neb_tmpdir
