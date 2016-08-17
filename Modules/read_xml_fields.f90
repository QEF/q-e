!
!----------------------------------------------------------------!
!   This module handles the reading of fields in xml inputs      !
!                                                                !
!   written by Simone Ziraldo (08/2010)                          !
!----------------------------------------------------------------!
!
!---------------------------------------------
! TB
! included monopole related stuff, search 'TB'
!---------------------------------------------
!
MODULE read_xml_fields_module
  !
  !
  USE io_global, ONLY : xmlinputunit => qestdin
  USE kinds,     ONLY : DP
  USE input_parameters 
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: read_xml_fields, clean_str
  !
  ! ... temporary variable needed to rebuild the old input_dft variable
  CHARACTER (len = 5) :: exchange, exchange_grad_corr, correlation, &
       correlation_grad_corr, xc_specials
  !
CONTAINS
  !
  !
  !
  !----------------------------------------------------------!
  !    This subroutine does a loop over all fields and       !
  !    sets the parameters that reads in these nodes         !
  !    note: the current implementation doesn't require      !
  !    the fields name or the number of fields               !
  !----------------------------------------------------------!
  SUBROUTINE read_xml_fields ()
    !
    !  
    USE iotk_module, ONLY : iotk_scan_begin, iotk_scan_end, iotk_scan_attr, iotk_attlenx
    USE iotk_unit_interf, ONLY : iotk_rewind
    !
    !
    IMPLICIT NONE
    !
    !
    INTEGER :: ierr, direction1, direction2
    CHARACTER(len = iotk_attlenx) :: attr
    CHARACTER(len = 30) :: name
    CHARACTER(len = 30) :: field
    !
    !
    ! ... the scanning must start from the beginning of the root node
    !
    CALL iotk_rewind( xmlinputunit )
    !
    ! ... default values for strings
    exchange = 'none'
    exchange_grad_corr = 'none'
    correlation = 'none'
    correlation_grad_corr = 'none'
    xc_specials = 'none'
    !
    ! ... fields loop
    !
    DO
       !
       call iotk_scan_begin( xmlinputunit, 'field', attr, direction = direction1, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error scanning begin of field &
           &node', ABS( ierr ) )
       !
       IF ( direction1 == -1 ) THEN
          !
          ! ... the scanning changes direction -> no more fields
          call iotk_scan_end( xmlinputunit, 'field' )
          !
          EXIT
          !
       END IF
       !
       call iotk_scan_attr(attr, 'name', field, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error getting the name of field', &
             ABS( ierr ) )
       !
       ! ... parameters loop
       !
       DO
          !
          CALL iotk_scan_begin(xmlinputunit, 'parameter', attr, direction = direction2, ierr = ierr )
          IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error scanning begin of parameter &
              &node inside '//trim(field)//' field', ABS( ierr ) )
          !
          IF ( direction2 == -1 ) THEN
             !
             ! ... the scanning changes direction -> no more parameters
             CALL iotk_scan_end( xmlinputunit, 'parameter', ierr = ierr)
             !
             EXIT
          END IF
          !
          CALL iotk_scan_attr( attr, 'name', name, ierr = ierr )
          IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error scanning the name of a PARAMETER &
              &inside '//trim(field)//' field', ABS( ierr ) )
          !
          !
          ! ... association string -> name of variable and reading of its value
          CALL read_parameter( name )
          !
          !
          CALL iotk_scan_end( xmlinputunit, 'parameter', ierr = ierr )
          IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error scanning end of '//name//' PARAMETER &
               &inside '//trim(field)//' field', ABS( ierr ) )
          !
       END DO
       !
       !
       call iotk_scan_end( xmlinputunit, 'field', ierr = ierr )
       IF (ierr /= 0) CALL errore( 'read_xml_fields', 'error scanning end of '//field//' field', 1)
       !
    END DO
    !
    ! ... reconstruction of input_dft variable ( parameter used in the old input format )
    !
    ! ... if one of the parameter is set
    IF ( (trim(exchange) /= 'none') .or. (trim(exchange_grad_corr) /= 'none')  &
         .or. (trim(correlation) /= 'none') .or.  (trim(correlation_grad_corr) /= 'none') ) THEN
       !
       ! ... all the parameter must be set
       IF ( (trim(exchange) /= 'none') .and. (trim(exchange_grad_corr) /= 'none')  &
            .and. (trim(correlation) /= 'none') .and.  (trim(correlation_grad_corr) /= 'none') ) THEN

          input_dft = trim(exchange)//'-'//trim(exchange_grad_corr)//'-'&
               //trim(correlation)//'-'//trim(correlation_grad_corr)
       ELSE
          ! ... error: at least one parameter is not set
          CALL errore( 'read_xml_fields', 'all the parameters exchange, exchange_grad_corr, &
               &correlation and correlation_grad_corr must be set', 1 )
          !
       ENDIF
    ELSE
       IF  (trim(xc_specials) /= 'none') input_dft = trim(xc_specials)
    END IF
    !
    RETURN
    !
    !
  END SUBROUTINE read_xml_fields
  !
  !
  !
  !--------------------------------------------------------------!
  !    This routine takes the parameter name as an input and     !
  !    with it reads the correspondig parameter                  !
  !--------------------------------------------------------------!
  SUBROUTINE read_parameter ( name )
    !
    USE iotk_module,ONLY : iotk_scan_dat_inside
    !
    IMPLICIT NONE
    !
    !
    CHARACTER ( len = * ), INTENT(IN) :: name
    !
    INTEGER :: ierr
    !
    ! ... temporary buffers needed for the reading of strings
    !
    CHARACTER ( len = 256 ) :: tmpstr
    !
    ierr = 0
    !
    ! ... list and reading of al the possible parameters
    !
    SELECT CASE (name(1:len_trim(name)))
       !
       !
    CASE ( 'abivol' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, abivol, ierr = ierr )
       !
    CASE ( 'adapt' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, adapt, ierr = ierr )
       !
    CASE ( 'ampre' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ampre, ierr = ierr )
       !
    CASE ( 'assume_isolated' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       assume_isolated = clean_str(tmpstr)
       !
    CASE ( 'bfgs_ndim' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, bfgs_ndim, ierr = ierr )
       !
    CASE ( 'calwf' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, calwf, ierr = ierr )
       !
    CASE ( 'cell_damping' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, cell_damping, ierr = ierr )
       !
    CASE ( 'cell_dofree' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       cell_dofree = clean_str(tmpstr)
       !
    CASE ( 'cell_dynamics' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       cell_dynamics = clean_str(tmpstr)
       !
    CASE ( 'cell_factor' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, cell_factor, ierr = ierr )
       !
    CASE ( 'cell_parameters' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       cell_parameters = clean_str(tmpstr)
       !
    CASE ( 'cell_temperature' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       cell_temperature = clean_str(tmpstr)
       !
    CASE ( 'cell_velocities' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       cell_velocities = clean_str(tmpstr)
       !
    CASE ( 'constrained_magnetization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       constrained_magnetization = clean_str(tmpstr)
       !
    CASE ( 'conv_thr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, conv_thr, ierr = ierr )
       !
    CASE ( 'correlation' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       correlation = clean_str(tmpstr)
       !
    CASE ( 'correlation_grad_corr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       correlation_grad_corr = clean_str(tmpstr)
       !
    CASE ( 'degauss' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, degauss, ierr = ierr )
       !
    CASE ( 'delta_t' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, delta_t, ierr = ierr )
       !
    CASE ( 'diago_cg_maxiter' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, diago_cg_maxiter, ierr = ierr )
       !
    CASE ( 'diago_david_ndim' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, diago_david_ndim, ierr = ierr )
       !
    CASE ( 'diago_full_acc' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, diago_full_acc, ierr = ierr )
       !
    CASE ( 'diago_thr_init' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, diago_thr_init, ierr = ierr )
       !
    CASE ( 'diagonalization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       diagonalization = clean_str(tmpstr)
       !
    CASE ( 'dipfield' )
       CALL iotk_scan_dat_inside( xmlinputunit, dipfield, ierr = ierr )
       !
    CASE ( 'disk_io' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       disk_io = clean_str(tmpstr)
       !
    CASE ( 'dthr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, dthr, ierr = ierr )
       !
    CASE ( 'dt' )
       CALL iotk_scan_dat_inside( xmlinputunit, dt, ierr = ierr )
       !
    CASE ( 'ecfixed' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ecfixed, ierr = ierr )
       !
    CASE ( 'ecutrho' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ecutrho, ierr = ierr )
       !
    CASE ( 'ecutwfc' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ecutwfc, ierr = ierr )
       !
    CASE ( 'edir' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, edir, ierr = ierr )
       !
    CASE ( 'efield' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, efield, ierr = ierr )
       !
    CASE ( 'efield_cart' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, efield_cart, ierr = ierr )
       !
    CASE ( 'efield_phase' )
       CALL iotk_scan_dat_inside( xmlinputunit, efield_phase, ierr = ierr )
       !
    CASE ( 'efx0' )
       CALL iotk_scan_dat_inside( xmlinputunit, efx0, ierr = ierr )
       !
    CASE ( 'efx1' )
       CALL iotk_scan_dat_inside( xmlinputunit, efx1, ierr = ierr )
       !
    CASE ( 'efy0' )
       CALL iotk_scan_dat_inside( xmlinputunit, efy0, ierr = ierr )
       !
    CASE ( 'efy1' )
       CALL iotk_scan_dat_inside( xmlinputunit, efy1, ierr = ierr )
       !
    CASE ( 'efz0' )
       CALL iotk_scan_dat_inside( xmlinputunit, efz0, ierr = ierr )
       !
    CASE ( 'efz1' )
       CALL iotk_scan_dat_inside( xmlinputunit, efz1, ierr = ierr )
       !
    CASE ( 'electron_damping' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, electron_damping, ierr = ierr )
       !
    CASE ( 'electron_dynamics' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       electron_dynamics = clean_str(tmpstr)
       !
    CASE ( 'electron_temperature' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       electron_temperature = clean_str(tmpstr)
       !
    CASE ( 'electron_velocities' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       electron_velocities = clean_str(tmpstr)
       !
    CASE ( 'eamp' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, eamp, ierr = ierr )
       !
    CASE ( 'ekin_conv_thr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ekin_conv_thr, ierr = ierr )
       !
    CASE ( 'ekincw' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ekincw, ierr = ierr )
       !
    CASE ( 'electron_maxstep' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, electron_maxstep, ierr = ierr )
       !
    CASE ( 'scf_must_converge' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, scf_must_converge, ierr = ierr )
       !
    CASE ( 'emass' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, emass, ierr = ierr )
       !
    CASE ( 'emass_cutoff' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, emass_cutoff, ierr = ierr )
       !
    CASE ( 'emaxpos' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, emaxpos, ierr = ierr )
       !
    CASE ( 'eopreg' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, eopreg, ierr = ierr )
       !
    CASE ( 'epol' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, epol, ierr = ierr )
       !
    CASE ( 'etot_conv_thr' )
       CALL iotk_scan_dat_inside( xmlinputunit, etot_conv_thr, ierr = ierr )
       !
    CASE ( 'exchange' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       exchange = clean_str(tmpstr)
       !
    CASE ( 'exchange_grad_corr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       exchange_grad_corr = clean_str(tmpstr)
       !
    CASE ( 'fixed_magnetization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, fixed_magnetization, ierr = ierr )
       !
    CASE ( 'fnosee' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, fnosee, ierr = ierr )
       !
    CASE ( 'fnoseh' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, fnoseh, ierr = ierr )
       !
    CASE ( 'fnosep' )
       CALL iotk_scan_dat_inside( xmlinputunit, fnosep, ierr = ierr )
       !
    CASE ( 'forc_conv_thr' )
       CALL iotk_scan_dat_inside( xmlinputunit, forc_conv_thr, ierr = ierr )
       !
    CASE ( 'force_symmorphic' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, force_symmorphic, ierr = ierr )
       !
    CASE ( 'gdir' )
       CALL iotk_scan_dat_inside( xmlinputunit, gdir, ierr = ierr )
       !
    CASE ( 'grease' )
       CALL iotk_scan_dat_inside( xmlinputunit, grease, ierr = ierr )
       !
    CASE ( 'greash' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, greash, ierr = ierr )
       !
    CASE ( 'greasp' )
       CALL iotk_scan_dat_inside( xmlinputunit, greasp, ierr = ierr )
       !
    CASE ( 'iprint' )
       CALL iotk_scan_dat_inside( xmlinputunit, iprint, ierr = ierr )
       !
    CASE ( 'ion_damping' )
       CALL iotk_scan_dat_inside( xmlinputunit, ion_damping, ierr = ierr )
       !
    CASE ( 'ion_dynamics' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       ion_dynamics = clean_str(tmpstr)
       !
    CASE ( 'ion_nstepe' )
       CALL iotk_scan_dat_inside( xmlinputunit, ion_nstepe, ierr = ierr )
       !
    CASE ( 'ion_positions' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       ion_positions = clean_str(tmpstr)
       !
    CASE ( 'ion_temperature' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       ion_temperature = clean_str(tmpstr)
       !
    CASE ( 'ion_velocities' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       ion_velocities = clean_str(tmpstr)
       !
    CASE ( 'isave' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, isave, ierr = ierr )
       !
    CASE ( 'la2F' )
       CALL iotk_scan_dat_inside( xmlinputunit, la2F, ierr = ierr )
       !
    CASE ( 'lambda' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, lambda, ierr = ierr )
       !
    CASE ( 'lambda_cold' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, lambda_cold, ierr = ierr )
       !
    CASE ( 'lberry' )
       CALL iotk_scan_dat_inside( xmlinputunit, lberry, ierr = ierr )
       !
    CASE ( 'lda_plus_u' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, lda_plus_u, ierr = ierr )
       !
    CASE ( 'lda_plus_u_kind' )
       CALL iotk_scan_dat_inside( xmlinputunit, lda_plus_u_kind, ierr = ierr )
       !
    CASE ( 'lelfield' )
       CALL iotk_scan_dat_inside( xmlinputunit, lelfield, ierr = ierr )
       !
    CASE ( 'lorbm' )
       CALL iotk_scan_dat_inside( xmlinputunit, lorbm, ierr = ierr )
       !
    CASE ( 'lkpoint_dir' )
       CALL iotk_scan_dat_inside( xmlinputunit, lkpoint_dir, ierr = ierr )
       !
    CASE ( 'london' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, london, ierr = ierr )
       !
    CASE ( 'london_rcut' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, london_rcut, ierr = ierr )
       !
    CASE ( 'london_s6' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, london_s6, ierr = ierr )
       !
    CASE ( 'lspinorb' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, lspinorb, ierr = ierr )
       !
    CASE ( 'lforcet' )
       CALL iotk_scan_dat_inside( xmlinputunit, lforcet, ierr = ierr )
       !
    CASE ( 'max_seconds' )
       CALL iotk_scan_dat_inside( xmlinputunit, max_seconds, ierr = ierr )
       !
    CASE ( 'maxiter' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, maxiter, ierr = ierr )
       !
    CASE ( 'maxwfdt' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, maxwfdt, ierr = ierr )
       !
    CASE ( 'mixing_beta' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, mixing_beta, ierr = ierr )
       !
    CASE ( 'mixing_fixed_ns' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, mixing_fixed_ns, ierr = ierr )
       !
    CASE ( 'mixing_mode' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       mixing_mode = clean_str(tmpstr)
       !
    CASE ( 'mixing_ndim' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, mixing_ndim, ierr = ierr )
       !
       ! TB - start
    CASE ( 'monopole' )
       CALL iotk_scan_dat_inside( xmlinputunit, monopole, ierr = ierr )
       !
    CASE ( 'zmon' )
       CALL iotk_scan_dat_inside( xmlinputunit, zmon, ierr = ierr )
       !
    CASE ( 'relaxz' )
       CALL iotk_scan_dat_inside( xmlinputunit, relaxz, ierr = ierr )
       !
    CASE ( 'block' )
       CALL iotk_scan_dat_inside( xmlinputunit, block, ierr = ierr )
       !
    CASE ( 'block_1' )
       CALL iotk_scan_dat_inside( xmlinputunit, block_1, ierr = ierr )
       !
    CASE ( 'block_2' )
       CALL iotk_scan_dat_inside( xmlinputunit, block_2, ierr = ierr )
       !
    CASE ( 'block_height' )
       CALL iotk_scan_dat_inside( xmlinputunit, block_height, ierr = ierr )
       ! TB - end
       !
    CASE ( 'n_inner' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, n_inner, ierr = ierr )
       !
    CASE ( 'nberrycyc' )
       CALL iotk_scan_dat_inside( xmlinputunit, nberrycyc, ierr = ierr )
       !
    CASE ( 'nbnd' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nbnd, ierr = ierr )
       !
    CASE ( 'ndega' )
       CALL iotk_scan_dat_inside( xmlinputunit, ndega, ierr = ierr )
       !
    CASE ( 'ndr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ndr, ierr = ierr )
       !
    CASE ( 'ndw' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ndw, ierr = ierr )
       !
    CASE ( 'nhpcl' )
       CALL iotk_scan_dat_inside( xmlinputunit, nhpcl, ierr = ierr )
       !
    CASE ( 'nhptyp' )
       CALL iotk_scan_dat_inside( xmlinputunit, nhptyp, ierr = ierr )
       !
    CASE ( 'niter_cold_restart' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, niter_cold_restart, ierr = ierr )
       !
    CASE ( 'nit' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nit, ierr = ierr )
       !
    CASE ( 'niter_cg_restart' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, niter_cg_restart, ierr = ierr )
       !
    CASE ( 'noinv' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, noinv, ierr = ierr )
       !
    CASE ( 'noncolin' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, noncolin, ierr = ierr )
       !
    CASE ( 'nosym_evc' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nosym_evc, ierr = ierr )
       !
    CASE ( 'nosym' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nosym, ierr = ierr )
       !
    CASE ( 'nppstr' )
       CALL iotk_scan_dat_inside( xmlinputunit, nppstr, ierr = ierr )
       !
    CASE ( 'nr1' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr1, ierr = ierr )
       !
    CASE ( 'nr1b' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr1b, ierr = ierr )
       !
   CASE ( 'nr1s' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr1s, ierr = ierr )
       !
    CASE ( 'nr2' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr2, ierr = ierr )
       !
    CASE ( 'nr2b' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr2b, ierr = ierr )
       !
    CASE ( 'nr2s' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr2s, ierr = ierr )
       !
    CASE ( 'nr3' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr3, ierr = ierr )
       !
    CASE ( 'nr3b' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr3b, ierr = ierr )
       !
    CASE ( 'nr3s' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr3s, ierr = ierr )
       !
    CASE ( 'nraise' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nraise, ierr = ierr )
       !
    CASE ( 'nsd' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nsd, ierr = ierr )
       !
    CASE ( 'nspin' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nspin, ierr = ierr )
       !
    CASE ( 'nstep' )
       CALL iotk_scan_dat_inside( xmlinputunit, nstep, ierr = ierr )
       !
    CASE ( 'nsteps' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nsteps, ierr = ierr )
       !
    CASE ( 'nwf' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nwf, ierr = ierr )
       !
    CASE ( 'occupations' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       occupations = clean_str(tmpstr)
       !
    CASE ( 'ortho_eps' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ortho_eps, ierr = ierr )
       !
    CASE ( 'ortho_max' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ortho_max, ierr = ierr )
       !
    CASE ( 'orthogonalization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       orthogonalization = clean_str(tmpstr)
       !
    CASE ( 'outdir' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       outdir = clean_str(tmpstr)
       !
    CASE ( 'P_ext' )
       CALL iotk_scan_dat_inside( xmlinputunit, P_ext, ierr = ierr )
       !
    CASE ( 'P_fin' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, P_fin, ierr = ierr )
       !
    CASE ( 'P_in' )
       CALL iotk_scan_dat_inside( xmlinputunit, P_in, ierr = ierr )
       !
    CASE ( 'passop' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, passop, ierr = ierr )
       !
    CASE ( 'pot_extrapolation' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       pot_extrapolation = clean_str(tmpstr)
       !
    CASE ( 'press' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, press, ierr = ierr )
       !
    CASE ( 'press_conv_thr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, press_conv_thr, ierr = ierr )
       !
    CASE ( 'pseudo_dir' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       pseudo_dir = clean_str(tmpstr)
       !
    CASE ( 'pvar' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, pvar, ierr = ierr )
       !
    CASE ( 'q2sigma' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, q2sigma, ierr = ierr )
       !
    CASE ( 'qcutz' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, qcutz, ierr = ierr )
       !
    CASE ( 'refold_pos' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, refold_pos, ierr = ierr )
       !
    CASE ( 'report' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, report, ierr = ierr )
       !
    CASE ( 'remove_rigid_rot' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, remove_rigid_rot, ierr = ierr )
       !
    CASE ( 'restart_mode' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       restart_mode = clean_str(tmpstr)
       !
    CASE ( 'rho_thr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, rho_thr, ierr = ierr )
       !
    CASE ( 'saverho' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, saverho, ierr = ierr )
       !
    CASE ( 'smearing' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       smearing = clean_str(tmpstr)
       !
    CASE ( 'startingpot' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       startingpot = clean_str(tmpstr)
       !
    CASE ( 'startingwfc' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       startingwfc = clean_str(tmpstr)
       !
    CASE ( 'Surf_t' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, Surf_t, ierr = ierr )
       !
    CASE ( 'sw_len' )
       CALL iotk_scan_dat_inside( xmlinputunit, sw_len, ierr = ierr )
       !
    CASE ( 'tabps' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tabps, ierr = ierr )
       !
    CASE ( 'tcg' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tcg, ierr = ierr )
       !
    CASE ( 'tefield' )
       CALL iotk_scan_dat_inside( xmlinputunit, tefield, ierr = ierr )
       !
    CASE ( 'temph' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, temph, ierr = ierr )
       !
    CASE ( 'tempw' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tempw, ierr = ierr )
       !
    CASE ( 'tolp' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tolp, ierr = ierr )
       !
    CASE ( 'tot_charge' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tot_charge, ierr = ierr )
       !
    CASE ( 'tot_magnetization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tot_magnetization, ierr = ierr )
       !
    CASE ( 'tolw' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tolw, ierr = ierr )
       !
    CASE ( 'tprnfor' )
       CALL iotk_scan_dat_inside( xmlinputunit, tprnfor, ierr = ierr )
       !
    CASE ( 'tqr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tqr, ierr = ierr )
       !
    CASE ( 'tq_smoothing' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tq_smoothing, ierr = ierr )
       !
    CASE ( 'tbeta_smoothing' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tbeta_smoothing, ierr = ierr )
       !
    CASE ( 'trust_radius_ini' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, trust_radius_ini, ierr = ierr )
       !
    CASE ( 'trust_radius_max' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, trust_radius_max, ierr = ierr )
       !
    CASE ( 'trust_radius_min' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, trust_radius_min, ierr = ierr )
       !
    CASE ( 'tstress' )
       CALL iotk_scan_dat_inside( xmlinputunit, tstress, ierr = ierr )
       !
    CASE ( 'U_projection_type' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       U_projection_type = clean_str(tmpstr)
       !
    CASE ( 'upscale' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, upscale, ierr = ierr )
       !
    CASE ( 'verbosity' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       verbosity = clean_str(tmpstr)
       !
    CASE ( 'w_1' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, w_1, ierr = ierr )
       !
    CASE ( 'w_2' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, w_2, ierr = ierr )
       !
    CASE ( 'wf_collect' )
       CALL iotk_scan_dat_inside( xmlinputunit, wf_collect, ierr = ierr )
       !
    CASE ( 'wf_efield' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, wf_efield, ierr = ierr )
       !
    CASE ( 'wf_friction' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, wf_friction, ierr = ierr )
       !
    CASE ( 'wf_q' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, wf_q, ierr = ierr )
       !
    CASE ( 'wf_switch' )
       CALL iotk_scan_dat_inside( xmlinputunit, wf_switch, ierr = ierr )
       !
    CASE ( 'wfc_extrapolation' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       wfc_extrapolation = clean_str(tmpstr)
       !
    CASE ( 'wfcdir' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       wfcdir = clean_str(tmpstr)
       !
    CASE ( 'wfdt' )
       CALL iotk_scan_dat_inside( xmlinputunit, wfdt, ierr = ierr )
       !
    CASE ( 'wffort' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, wffort, ierr = ierr )
       !
    CASE ( 'wfsd' )
       CALL iotk_scan_dat_inside( xmlinputunit, wfsd, ierr = ierr )
       !
    CASE ( 'wmass' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, wmass, ierr = ierr )
       !
    CASE ( 'writev' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, writev, ierr = ierr )
       !
    CASE ( 'xc_specials' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       xc_specials = clean_str(tmpstr)
       ! ... up to here
       !
    CASE ( 'xmloutput' )
       CALL iotk_scan_dat_inside( xmlinputunit, xmloutput, ierr = ierr )
       !
    case ( 'vdw_table_name' )
      CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
      vdw_table_name = clean_str(tmpstr)
       !
    case ( 'input_dft' )
      CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
      input_dft = clean_str(tmpstr)
 
    CASE default
       !
       CALL errore( 'read_parameter', 'no parameter with name '//trim( name ), 1 )
       !
       !
    END SELECT
    !
    IF (ierr/=0) THEN
       CALL errore( 'read_parameter', 'problem reading parameter '//trim( name ), 1 )
    END IF
    !
    !
    RETURN
    !
    !
  END SUBROUTINE read_parameter
  !
  !
  !---------------------------------------------------------!
  !  Function that eliminate the tab characters and adjust  !
  !  to the left side the string                            !
  !---------------------------------------------------------!
  FUNCTION clean_str( string )
    !
    !
    IMPLICIT NONE
    !
    !
    CHARACTER (len = *) :: string
    CHARACTER (len = len( string ) ) :: clean_str
    INTEGER :: i
    !
    do i = 1, len( string )
       !
       if ( ichar( string(i:i) ) == 9 ) then
          clean_str(i:i)=' '
       else
          clean_str(i:i)=string(i:i)
       end if
       !
    end do
    !
    clean_str = adjustl( clean_str )
    !
    !
  END FUNCTION clean_str
  !
  !
  !
END MODULE read_xml_fields_module
