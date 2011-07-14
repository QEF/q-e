!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_readin
  !-----------------------------------------------------------------------
  !
  !    This routine reads the control variables from standard input (unit 5).
  !    A second routine read_file reads the variables saved to file
  !    by the self-consistent program.
  !
  ! Modified by Osman Baris Malcioglu (2009)
#include "f_defs.h"

  USE lr_variables
  USE kinds,               ONLY : DP
  USE io_files,            ONLY : tmp_dir, prefix,wfc_dir
  USE lsda_mod,            ONLY : current_spin, nspin
  USE control_flags,       ONLY : twfcollect,use_para_diag
  USE scf,                 ONLY : vltot, v, vrs, vnew, &
                                   & destroy_scf_type
  USE fft_base,            ONLY : dfftp
  USE gvecs,               ONLY : doublegrid
  USE wvfct,               ONLY : nbnd, et, wg
  USE lsda_mod,            ONLY : isk
  USE ener,                ONLY : ef
  USE io_global,           ONLY : ionode, ionode_id
  USE klist,               ONLY : nks, wk, nelec
  USE fixed_occ,           ONLY : tfixed_occ
  USE input_parameters,    ONLY : degauss, nosym,wfcdir,outdir
  USE ktetra,              ONLY : ltetra
  USE realus,              ONLY : real_space, real_space_debug, &
                                   init_realspace_vars, qpointlist, &
                                   betapointlist,read_rs_status,newd_r
  USE funct,               ONLY : dft_is_meta
  USE io_global,           ONLY : stdout
  USE control_flags,       ONLY : tqr
  USE iotk_module
  USE charg_resp,          ONLY : w_T_prefix, omeg, w_T_npol, epsil
  USE mp,        ONLY : mp_bcast,mp_barrier
  USE mp_global, ONLY : my_pool_id, intra_image_comm, intra_pool_comm
  USE io_global, ONLY : ionode, ionode_id
  USE DFUNCT,         ONLY : newd
  USE vlocal,         ONLY : strf
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256) :: beta_gamma_z_prefix
          ! fine control of beta_gamma_z file
  CHARACTER(len=80) :: disk_io
          ! Specify the amount of I/O activities
  INTEGER :: ios, iunout,ierr,ipol
  LOGICAL :: auto_rs
  !
  NAMELIST / lr_input / restart, restart_step ,lr_verbosity, prefix, outdir, test_case_no, wfcdir,disk_io
  NAMELIST / lr_control / itermax, ipol, ltammd, real_space, real_space_debug, charge_response, tqr, auto_rs, no_hxc,n_ipol,project
  NAMELIST / lr_post / omeg, beta_gamma_z_prefix, w_T_npol, plot_type, epsil,itermax_int
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_readin>")')
  ENDIF
  auto_rs = .true.
#ifdef __PARA
  IF (ionode) THEN
#endif
  !
  !   Set default values for variables in namelist
  !
  itermax = 500
  restart = .false.
  restart_step = itermax+1
  lr_verbosity = 1
  prefix = 'pwscf'
  disk_io = 'default'
  ltammd = .false.
  ipol = 1
  n_ipol=1
  no_hxc = .false.
  !broadening = 0.0d0
  real_space = .false.
  real_space_debug = 0
  charge_response = 0
  test_case_no = 0
  tqr = .false.
  auto_rs = .true.
  beta_gamma_z_prefix = 'undefined'
  omeg=0.0
  epsil=0.0
  w_T_npol=1
  plot_type=1
 project=.false.
  !terminator="no"
  !grid_coarsening = 1
  !
  !   Reading the namelist lr_input
  !
  CALL input_from_file( )
  !
  READ (5, lr_input, err = 200, iostat = ios)
  200 CALL errore ('lr_readin', 'reading lr_input namelist', abs (ios) )
  !
  !
  !   Reading the namelist lr_control
  !
  READ (5, lr_control, err = 201, iostat = ios)
  201 CALL errore ('lr_readin', 'reading lr_control namelist', abs (ios) )
  !
  !
  !   Reading the namelist lr_post
  !
  IF (charge_response == 1) THEN
   READ (5, lr_post, err = 202, iostat = ios)
   202 CALL errore ('lr_readin', 'reading lr_post namelist', abs (ios) )
   bgz_suffix = trim ( "-stage2.beta_gamma_z." )
   WRITE(stdout,'(/5x,"Prefix of current run is appended by -stage2")')
   IF ( beta_gamma_z_prefix  == 'undefined' ) THEN
    beta_gamma_z_prefix=trim(prefix)
   ENDIF
  ELSE
   bgz_suffix = trim ( ".beta_gamma_z." )
  ENDIF
  !
  ! The status of the real space flags should be read manually
  !
  ! Do not mess with already present wfc structure
  twfcollect = .false.
  !
  outdir = trimcheck(outdir)
  tmp_dir = outdir
  !
  !
  IF ( .not. trim( wfcdir ) == 'undefined' ) THEN
     !
     wfc_dir = trimcheck ( wfcdir )
     !
     !CALL verify_tmpdir( wfc_dir ) ! do not call, since this may erase files
     !
  ENDIF
  !
  !Charge response mode 1 is the "do Lanczos chains twice, conserve memory" scheme
  !
  IF (charge_response == 1 .and. omeg == 0.D0) &
   CALL errore ('lr_readin', 'omeg must be defined for charge response mode 1', 1 )
  IF ( project .and. charge_response /= 1) &
   CALL errore ('lr_readin', 'projection is possible only in charge response mode 1', 1 )

  w_T_prefix = trim( tmp_dir ) // trim( beta_gamma_z_prefix ) // ".beta_gamma_z."
  !

  ierr = 0
  !
  IF ( ipol==4 ) THEN
     !
     n_ipol = 3
     LR_polarization=1
     !
  ELSE
     !
     LR_polarization=ipol
     !
  ENDIF
  IF (itermax_int < itermax) itermax_int=itermax
!
! Limited disk_io support: currently only one setting is supported
!
SELECT CASE( trim( disk_io ) )
  CASE ( 'reduced' )
     !
     lr_io_level = -1
     restart  = .false.
     !
  CASE DEFAULT
     !
     lr_io_level = 1
END SELECT

#ifdef __PARA
  ENDIF
  CALL bcast_lr_input
  CALL mp_bcast(auto_rs, ionode_id)
#endif
  !
  !print *, "post broad"
  !print *, "rs_status"
  outdir = trim( tmp_dir ) // trim( prefix ) // '.save'
  IF (auto_rs) CALL read_rs_status( outdir, ierr )
  IF (real_space) real_space_debug=99
  IF (real_space_debug > 0) real_space=.true.
  IF (lr_verbosity > 1) THEN
      WRITE(stdout,'(5x,"Status of real space flags: TQR=", L5 ,"  REAL_SPACE=", L5)') tqr, real_space
   ENDIF
  !print *, "rs_status-ended"
  !   Now PWSCF XML file will be read, and various initialisations will be done
  !
 !print *, "newd"
 !
  !print *, "read_file"
  !call mp_barrier()
  CALL read_file()

  DEALLOCATE( strf )
  CALL destroy_scf_type(vnew)
  !
  !
  !

  !
  !
  !   Re-initialize all needed quantities from the scf run
  !
  !print *, "init_us_1"
  current_spin=1
  !

  CALL init_us_1 ( )
  !
  IF (tqr) THEN
   CALL newd_r()
  ELSE
   CALL newd() !OBM: this is for the ground charge density
  ENDIF
  !
  IF (dft_is_meta()) &
   CALL errore( ' iosys ', ' Meta DFT ' // &
       & 'not implemented yet', 1 )

   IF ( real_space_debug > 0 ) THEN
   WRITE(stdout,'(/5x,"Real space implementation V.1 D190908",1x)')
    !  !OBM - correct parellism issues
    CALL init_realspace_vars()
    CALL betapointlist()
    WRITE(stdout,'(5X,"Real space initialisation completed")')
  ENDIF

  !
  !print *, "set_vrs"
  CALL set_vrs ( vrs, vltot, v%of_r, 0, 0, dfftp%nnr, nspin, doublegrid )

  DEALLOCATE( vltot )
  CALL destroy_scf_type(v)

  CALL iweights( nks, wk, nbnd, nelec, et, ef, wg, 0, isk)
  !
  !
  IF ( charge_response == 2 ) CALL lr_set_boxes_density()
  !
  !   Checking
  !
  IF ( (ltetra .or. tfixed_occ .or. (degauss /= 0.D0)) ) &
       CALL errore( ' iosys ', ' Linear response calculation ' // &
       & 'not implemented for non-insulating systems', 1 )
  IF ( .not. nosym ) &
       CALL errore( ' iosys ', ' Linear response calculation ' // &
       & 'not implemented with symmetry', 1 )
  !
  !
  !Scalapack related stuff, 
  !
#ifdef __PARA
  use_para_diag = .true.
  CALL check_para_diag( nbnd )
#else
  use_para_diag = .FALSE.
#endif

  RETURN
END SUBROUTINE lr_readin
