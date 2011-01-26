!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine lr_readin
  !-----------------------------------------------------------------------
  !
  !    This routine reads the control variables from standard input (unit 5).
  !    A second routine read_file reads the variables saved to file
  !    by the self-consistent program.
  !
  ! Modified by Osman Baris Malcioglu (2009)
#include "f_defs.h"

  use lr_variables
  USE kinds,               only : DP
  use io_files,            only : tmp_dir, prefix,trimcheck,wfc_dir
  use lsda_mod,            only : current_spin, nspin
  use control_flags,       only : twfcollect
  USE scf,                 ONLY : vltot, v, vrs
  USE grid_dimensions,     ONLY : nrxx
  USE gvecs,             ONLY : doublegrid
  use wvfct,               only : nbnd, et, wg
  use lsda_mod,            only : isk
  use ener,                only : ef
  USE io_global,           ONLY : ionode, ionode_id
  use klist,               only : nks, wk, nelec
  use fixed_occ,           only : tfixed_occ
  use input_parameters,    only : degauss, nosym,wfcdir,outdir
  use ktetra,              only : ltetra
  USE realus,              ONLY : real_space, real_space_debug, &
                                   init_realspace_vars, qpointlist, &
                                   betapointlist,read_rs_status,newd_r
  USE funct,               only : dft_is_meta
  USE io_global,           ONLY : stdout
  USE control_flags,       ONLY : tqr
  USE iotk_module
  USE charg_resp,          ONLY : w_T_prefix, omeg, w_T_npol, epsil
  USE mp,        ONLY : mp_bcast,mp_barrier
  USE mp_global, ONLY : my_pool_id, intra_image_comm, intra_pool_comm
  USE io_global, ONLY : ionode, ionode_id
  USE DFUNCT,         ONLY : newd
  implicit none
  !
  character(len=256) :: beta_gamma_z_prefix
  integer :: ios, iunout,ierr,ipol
  logical :: auto_rs
  !
  namelist / lr_input / restart, restart_step ,lr_verbosity, prefix, outdir, test_case_no, wfcdir
  namelist / lr_control / itermax, ipol, ltammd, real_space, real_space_debug, charge_response, tqr, auto_rs, no_hxc,n_ipol,project
  namelist / lr_post / omeg, beta_gamma_z_prefix, w_T_npol, plot_type, epsil,itermax_int
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_readin>")')
  endif
  auto_rs = .true.
#ifdef __PARA
  if (ionode) then
#endif
  !
  !   Set default values for variables in namelist
  !
  itermax = 500
  restart = .false.
  restart_step = itermax+1
  lr_verbosity = 1
  prefix = 'pwscf'
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
  call input_from_file( )
  !
  read (5, lr_input, err = 200, iostat = ios)
  200 call errore ('lr_readin', 'reading lr_input namelist', abs (ios) )
  ! 
  !
  !   Reading the namelist lr_control
  !
  read (5, lr_control, err = 201, iostat = ios)
  201 call errore ('lr_readin', 'reading lr_control namelist', abs (ios) )
  !
  !
  !   Reading the namelist lr_post
  !
  if (charge_response == 1) then
   read (5, lr_post, err = 202, iostat = ios)
   202 call errore ('lr_readin', 'reading lr_post namelist', abs (ios) )
   bgz_suffix = TRIM ( "-stage2.beta_gamma_z." )
   write(stdout,'(/5x,"Prefix of current run is appended by -stage2")')
   IF ( beta_gamma_z_prefix  == 'undefined' ) then
    beta_gamma_z_prefix=trim(prefix)
   endif
  else
   bgz_suffix = TRIM ( ".beta_gamma_z." )
  endif
  !
  ! The status of the real space flags should be read manually
  ! 
  ! Do not mess with already present wfc structure
  twfcollect = .FALSE.
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
  if (charge_response == 1 .and. omeg == 0.D0) &
   call errore ('lr_readin', 'omeg must be defined for charge response mode 1', 1 )
  if ( project .and. charge_response /= 1) &
   call errore ('lr_readin', 'projection is possible only in charge response mode 1', 1 )

  w_T_prefix = TRIM( tmp_dir ) // TRIM( beta_gamma_z_prefix ) // ".beta_gamma_z."
  !

  ierr = 0
  ! 
  if ( ipol==4 ) then
     !
     n_ipol = 3
     LR_polarization=1 
     !
  else
     !
     LR_polarization=ipol
     !
  end if
  if (itermax_int < itermax) itermax_int=itermax
#ifdef __PARA
  end if
  call bcast_lr_input
  call mp_bcast(auto_rs, ionode_id)
#endif
  !
  !print *, "post broad"
  !print *, "rs_status"
  outdir = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  if (auto_rs) call read_rs_status( outdir, ierr ) 
  if (real_space) real_space_debug=99
  if (real_space_debug > 0) real_space=.true. 
  If (lr_verbosity > 1) THEN
      WRITE(stdout,'(5x,"Status of real space flags: TQR=", L5 ,"  REAL_SPACE=", L5)') tqr, real_space
   endif
  !print *, "rs_status-ended"

  !   Here we finished the reading of the input file.
  !   Now PWSCF XML file will be read, and various initialisations will be done
  !
 !print *, "newd" 
 !
  !print *, "read_file"
  !call mp_barrier()
  call read_file() 
  !
  !
  !   Re-initialize all needed quantities from the scf run
  !
  !print *, "init_us_1"
  current_spin=1
  !

  call init_us_1 ( )
  !
  if (tqr) then
   call newd_r()
  else
   call newd() !OBM: this is for the ground charge density
  endif
  !
  if (dft_is_meta()) &
   CALL errore( ' iosys ', ' Meta DFT ' // &
       & 'not implemented yet', 1 )

   if ( real_space_debug > 0 ) then 
   write(stdout,'(/5x,"Real space implementation V.1 D190908",1x)')
    !  !OBM - correct parellism issues
    call init_realspace_vars()
    call betapointlist()
    write(stdout,'(5X,"Real space initialisation completed")')
  endif

  !
  !print *, "set_vrs"
  call set_vrs ( vrs, vltot, v%of_r, 0, 0, nrxx, nspin, doublegrid )
  !
  call iweights( nks, wk, nbnd, nelec, et, ef, wg, 0, isk)
  !
  !
  if ( charge_response == 2 ) call lr_set_boxes_density()
  !
  !   Checking 
  !
  IF ( (ltetra .OR. tfixed_occ .OR. (degauss /= 0.D0)) ) &
       CALL errore( ' iosys ', ' Linear response calculation ' // &
       & 'not implemented for non-insulating systems', 1 )
  IF ( .NOT. nosym ) &
       CALL errore( ' iosys ', ' Linear response calculation ' // &
       & 'not implemented with symmetry', 1 )
  !
  return
end subroutine lr_readin
