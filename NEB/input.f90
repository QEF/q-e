!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE ioneb(xmlinput,attr)
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine reads input data from standard input ( unit 5 )
  ! ...  Use "-input filename" to read input from file "filename":
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  ---------------------------------------------------------------
  !
  ! ...  access the modules renaming the variables that have the same name
  ! ...  as the input parameters, this is required in order to use a code
  ! ...  independent input parser
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : autoev, eV_to_kelvin, pi, rytoev, &
                            uakbar, amconv, bohr_radius_angs, eps8
  USE io_global,     ONLY : stdout, ionode
  !
  USE cell_base,     ONLY : at, bg, alat, omega, &
                            celldm_ => celldm, &
                            ibrav_  => ibrav, &
                            init_dofree
  !
  USE ions_base,     ONLY : if_pos, ityp, tau, extfor, &
                            ntyp_ => nsp, &
                            nat_  => nat, &
                            amass, tau_format
  !
  !
  USE io_files,      ONLY : tmp_dir 
  !
  USE force_mod,     ONLY : lforce, lstres, force
  !
  USE control_flags, ONLY : &
                            io_level, &
                            lscf, &
                            nstep, &
                            lpath, lneb,   &
                            lsmd,                    &
                            restart
  !
  USE path_variables, ONLY : nstep_path, lsteep_des, lquick_min, &
                             lbroyden, lbroyden2, &
                             llangevin, &
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
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : title, calculation, verbosity, restart_mode, &
                               iprint, tstress, tprnfor, dt, outdir, &
                               wfcdir, prefix, etot_conv_thr, forc_conv_thr, &
                               wf_collect
                               
  USE input_parameters, ONLY : &
                               num_of_images, path_thr, CI_scheme, opt_scheme, &
                               use_masses, first_last_opt, temp_req, k_max,    &
                               k_min, ds, use_freezing, fixed_tan
  !
  !
  ! ... "path" specific
  !
  USE input_parameters, ONLY : pos, full_phs_path_flag
  !
  USE read_namelists_module, ONLY : read_namelists, sm_not_set
  !
  USE read_xml_module,       ONLY : read_xml
  USE iotk_module,           ONLY : iotk_open_read, iotk_close_read,iotk_attlenx
  !
  IMPLICIT NONE
  !
  CHARACTER (len=iotk_attlenx), intent(in) :: attr
  LOGICAL, intent(in) :: xmlinput
  !
  INTEGER  :: ia, image, nt
  REAL(DP) :: theta, phi
  INTEGER  :: iiarg, nargs, iargc, ierr
  CHARACTER (len=50) :: arg
  !
  !
 !
  ! ... all namelists are read
  !
!  IF ( xmlinput ) THEN
!     CALL read_xml ('PW', 1 , attr = attr )
!  ELSE
write(0,*) "before read_namelist"
!     CALL read_namelists( 'SM' )
write(0,*) "after read_namelist"
!  ENDIF
  !
  !
!-----------------------------
! devono andare dopo il call a iosys
!---------------------------
SELECT CASE(trim( calculation ))
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
!     CALL errore( 'iosys', 'calculation ' // &
!                & trim( calculation ) // ' not implemented', 1 )
     !
  END SELECT
  !
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
        restart = .false.
        !
     ENDIF
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', &
                & 'unknown restart_mode ' // trim( restart_mode ), 1 )
     !
  END SELECT
  !
  IF ( lpath ) THEN
     !
write(0,*) "if lpath"
     IF( io_level < 0) CALL errore ( 'iosys', &
                       'NEB, SMD do not work with "disk_io" set to "none"', 1)
     !
     nstep_path = nstep
     !
     IF ( num_of_images < 2 ) &
        CALL errore( 'iosys', 'calculation=' // trim( calculation ) // &
                   & ': num_of_images must be at least 2', 1 )
     !
     IF ( ( CI_scheme /= "no-CI"  ) .and. &
          ( CI_scheme /= "auto"   ) .and. &
          ( CI_scheme /= "manual" ) ) THEN
        !
        CALL errore( 'iosys', 'calculation=' // trim( calculation ) // &
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
write(0,*) "case sd"
        lsteep_des = .true.
        !
     CASE( "quick-min" )
write(0,*) "case quick-min"
        !
        lquick_min = .true.
        !
     CASE( "broyden" )
        !
write(0,*) "case broyden"
        lbroyden     = .true.
        !
     CASE( "broyden2" )
write(0,*) "case broyden2"
        !
        lbroyden2    = .true.
        !
     CASE( "langevin" )
write(0,*) "case langevin"
        !
        llangevin = .true.
        !
        IF ( lneb ) &
           CALL errore( 'iosys','calculation=' // trim( calculation ) // &
                      & ': langevin dynamics not implemented', 1 )
        !
        temp_req = temp_req / ( eV_to_kelvin * autoev )
        !
        IF ( temp_req <= 0.D0 ) &
           CALL errore( 'iosys','calculation=' // trim( calculation ) // &
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
        CALL errore( 'iosys','calculation=' // trim( calculation ) // &
                   & ': unknown opt_scheme', 1 )
        !
     END SELECT
     !
  ENDIF
  !
  ! ... "path"-optimization variables
  !
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
  !
  ! ... read following cards
  !
!  ALLOCATE( ityp( nat_ ) )
!  ALLOCATE( tau(    3, nat_ ) )
!  ALLOCATE( force(  3, nat_ ) )
!  ALLOCATE( if_pos( 3, nat_ ) )
!  ALLOCATE( extfor( 3, nat_ ) )
!  IF ( tfixed_occ ) THEN
!     IF ( nspin_ == 4 ) THEN
!        ALLOCATE( f_inp( nbnd_, 1 ) )
!     ELSE
!        ALLOCATE( f_inp( nbnd_, nspin_ ) )
!     ENDIF
!  ENDIF
  !
!  IF ( tefield ) ALLOCATE( forcefield( 3, nat_ ) )
  !
!write(0,*) "before read cards pw"
!  CALL read_cards_pw ( psfile, tau_format, xmlinput )
!write(0,*) "after read cards pw"
  !
  ! ... set up atomic positions and crystal lattice
  !
  !
     !
     ! ... "path" optimizations specific
     !
     DO image = 1, num_of_images_
        !
        tau = reshape( pos(1:3*nat_,image), (/ 3 , nat_ /) )
        !
        CALL convert_tau ( tau_format, nat_, tau)
        !
        ! ... note that this positions array is in Bohr
        !
        pos(1:3*nat_,image) = reshape( tau, (/ 3 * nat_ /) ) * alat
        !
     ENDDO
     !
  !
!  CALL verify_tmpdir( tmp_dir )
write(0,*) "before verify neb dir"
  CALL verify_neb_tmpdir( tmp_dir )
write(0,*) "after verify neb dir"
  !
! uuuu questo bisogna vedere
     !
  RETURN
  !
END SUBROUTINE ioneb
!
!-----------------------------------------------------------------------
SUBROUTINE verify_neb_tmpdir( tmp_dir )
  !-----------------------------------------------------------------------
  !
  USE wrappers,         ONLY : f_mkdir
  USE input_parameters, ONLY : restart_mode
  USE control_flags,    ONLY : lpath, lbands
  USE io_files,         ONLY : prefix, xmlpun, &
                               delete_if_present, check_writable
  USE pw_restart,       ONLY : pw_readfile
  USE path_variables,   ONLY : num_of_images
  USE mp_global,        ONLY : mpime, nproc, nimage
  USE io_global,        ONLY : ionode
  USE mp,               ONLY : mp_barrier
  USE xml_io_base,      ONLY : copy_file
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(inout) :: tmp_dir
  !
  INTEGER             :: ios, image, proc, nofi
  LOGICAL             :: exst
  CHARACTER (len=256) :: file_path, filename
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  !
  file_path = trim( tmp_dir ) // trim( prefix )
  !
  ! .... check: is the directory writable?
  !
  ios = check_writable ( tmp_dir, mpime )
  !
  IF ( ios /= 0 ) THEN
     !
     ! ... no luck: directory non existent or non writable
     !
     IF ( restart_mode == 'from_scratch' ) THEN
        !
        ! ... let us try to create the scratch directory
        !
        CALL parallel_mkdir ( tmp_dir )
        !
     ELSE
        !
        CALL errore( 'outdir: ', trim( tmp_dir ) // &
                             & ' non existent or non writable', 1 )
        !
     ENDIF
     !
  ELSE
write(0,*) "verify neb dir ok"
     !
     ! ... if starting from scratch all temporary files are removed
     ! ... from tmp_dir ( only by the master node )
     !
     IF ( ionode ) THEN
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
        CALL parallel_mkdir ( file_path )
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
              CALL mp_barrier()
              !
           ENDDO
           !
        ENDIF ! end restart_mode
        !
     ENDDO ! end do image
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE verify_neb_tmpdir
