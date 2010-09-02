!
!----------------------------------------------------------------!
!   This module handles the reading of fields in xml inputs      !
!                                                                !
!   written by Simone Ziraldo (08/2010)                          !
!----------------------------------------------------------------!
MODULE read_xml_fields_module
  !
  !
  USE io_files, ONLY : xmlinputunit
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
       call iotk_scan_begin( xmlinputunit, 'FIELD', attr, direction = direction1, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error scanning begin of FIELD &
           &node', ABS( ierr ) )
       !
       IF ( direction1 == -1 ) THEN
          !
          ! ... the scanning changes direction -> no more fields
          call iotk_scan_end( xmlinputunit, 'FIELD' )
          !
          EXIT
          !
       END IF
       !
       call iotk_scan_attr(attr, 'name', field, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error getting the name of FIELD', &
             ABS( ierr ) )
       !
       ! ... parameters loop
       !
       DO
          !
          CALL iotk_scan_begin(xmlinputunit, 'PARAMETER', attr, direction = direction2, ierr = ierr )
          IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error scanning begin of PARAMETER &
              &node inside '//trim(field)//' field', ABS( ierr ) )
          !
          IF ( direction2 == -1 ) THEN
             !
             ! ... the scanning changes direction -> no more parameters
             CALL iotk_scan_end( xmlinputunit, 'PARAMETER', ierr = ierr)
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
          CALL iotk_scan_end( xmlinputunit, 'PARAMETER', ierr = ierr )
          IF ( ierr /= 0 ) CALL errore ( 'read_xml_fields', 'error scanning end of '//name//' PARAMETER &
               &inside '//trim(field)//' field', ABS( ierr ) )
          !
       END DO
       !
       !
       call iotk_scan_end( xmlinputunit, 'FIELD', ierr = ierr )
       IF (ierr /= 0) CALL errore( 'read_xml_fields', 'error scanning end of '//field//' field', 1)
       !
    END DO
    !
    ! ... reconstruction of input_dft variable ( parameter used in the old input format )
    !
    ! ... if one of the parameter is setted
    IF ( (trim(exchange) /= 'none') .or. (trim(exchange_grad_corr) /= 'none')  &
         .or. (trim(correlation) /= 'none') .or.  (trim(correlation_grad_corr) /= 'none') ) THEN
       !
       ! ... all the parameter must be setted
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
    ! ... temporary buffers needed for the reading
    !
    CHARACTER ( len = 256 ) :: tmpstr
    REAL( DP ), ALLOCATABLE, DIMENSION(:) :: tmparray
    !
    ierr = 0
    !
    ! ... list and reading of al the possible parameters
    !
    SELECT CASE (name(1:len_trim(name)))
       !
       !
    CASE ( 'verbosity' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       verbosity = clean_str(tmpstr)
       !
    CASE ( 'xmloutput' )
       CALL iotk_scan_dat_inside( xmlinputunit, xmloutput, ierr = ierr )
       !
    CASE ( 'la2F' )
       CALL iotk_scan_dat_inside( xmlinputunit, la2F, ierr = ierr )
       !
    CASE ( 'restart_mode' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       restart_mode = clean_str(tmpstr)
       !
    CASE ( 'wf_collect' )
       CALL iotk_scan_dat_inside( xmlinputunit, wf_collect, ierr = ierr )
       !
    CASE ( 'nstep' )
       CALL iotk_scan_dat_inside( xmlinputunit, nstep, ierr = ierr )
       !
    CASE ( 'iprint' )
       CALL iotk_scan_dat_inside( xmlinputunit, iprint, ierr = ierr )
       !
    CASE ( 'tstress' )
       CALL iotk_scan_dat_inside( xmlinputunit, tstress, ierr = ierr )
       !
    CASE ( 'tprnfor' )
       CALL iotk_scan_dat_inside( xmlinputunit, tprnfor, ierr = ierr )
       !
    CASE ( 'dt' )
       CALL iotk_scan_dat_inside( xmlinputunit, dt, ierr = ierr )
       !
    CASE ( 'outdir' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       outdir = clean_str(tmpstr)
       !
    CASE ( 'wfcdir' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       wfcdir = clean_str(tmpstr)
       !
    CASE ( 'lkpoint_dir' )
       CALL iotk_scan_dat_inside( xmlinputunit, lkpoint_dir, ierr = ierr )
       !
    CASE ( 'max_seconds' )
       CALL iotk_scan_dat_inside( xmlinputunit, max_seconds, ierr = ierr )
       !
    CASE ( 'etot_conv_thr' )
       CALL iotk_scan_dat_inside( xmlinputunit, etot_conv_thr, ierr = ierr )
       !
    CASE ( 'forc_conv_thr' )
       CALL iotk_scan_dat_inside( xmlinputunit, forc_conv_thr, ierr = ierr )
       !
    CASE ( 'disk_io' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       disk_io = clean_str(tmpstr)
       !
    CASE ( 'pseudo_dir' )
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       pseudo_dir = clean_str(tmpstr)
       !
    CASE ( 'tefield' )
       CALL iotk_scan_dat_inside( xmlinputunit, tefield, ierr = ierr )
       !
    CASE ( 'dipfield' )
       CALL iotk_scan_dat_inside( xmlinputunit, dipfield, ierr = ierr )
       !
    CASE ( 'lelfield' )
       CALL iotk_scan_dat_inside( xmlinputunit, lelfield, ierr = ierr )
       !
    CASE ( 'nberrycyc' )
       CALL iotk_scan_dat_inside( xmlinputunit, nberrycyc, ierr = ierr )
       !
    CASE ( 'lberry' )
       CALL iotk_scan_dat_inside( xmlinputunit, lberry, ierr = ierr )
       !
    CASE ( 'gdir' )
       CALL iotk_scan_dat_inside( xmlinputunit, gdir, ierr = ierr )
       !
    CASE ( 'nppstr' )
       CALL iotk_scan_dat_inside( xmlinputunit, nppstr, ierr = ierr )
       !
    CASE ( 'nbnd' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nbnd, ierr = ierr )
       !
    CASE ( 'tot_charge' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tot_charge, ierr = ierr )
       !
    CASE ( 'tot_magnetization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tot_magnetization, ierr = ierr )
       !
    CASE ( 'starting_magnetization' ) 
       IF (ntyp > 0) THEN
          ALLOCATE(tmparray(ntyp))
          CALL iotk_scan_dat_inside( xmlinputunit, tmparray, ierr = ierr )
          starting_magnetization(1:ntyp) = tmparray(1:ntyp)
          DEALLOCATE(tmparray)
       ELSE
          ierr = 1
       END IF
       !
    CASE ( 'ecutwfc' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ecutwfc, ierr = ierr )
       !
    CASE ( 'ecutrho' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ecutrho, ierr = ierr )
       !
    CASE ( 'nr1' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr1, ierr = ierr )
       !
    CASE ( 'nr2' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr2, ierr = ierr )
       !
    CASE ( 'nr3' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr3, ierr = ierr )
       !
    CASE ( 'nr1s' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr1s, ierr = ierr )
       !
    CASE ( 'nr2s' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr2s, ierr = ierr )
       !
    CASE ( 'nr3s' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nr3s, ierr = ierr )
       !
    CASE ( 'nosym' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nosym, ierr = ierr )
       !
    CASE ( 'nosym_evc' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nosym_evc, ierr = ierr )
       !
    CASE ( 'noinv' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, noinv, ierr = ierr )
       !
    CASE ( 'force_symmorphic' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, force_symmorphic, ierr = ierr )
       !
    CASE ( 'occupations' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       occupations = clean_str(tmpstr)
       !
    CASE ( 'degauss' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, degauss, ierr = ierr )
       !
    CASE ( 'smearing' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       smearing = clean_str(tmpstr)
       !
    CASE ( 'nspin' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nspin, ierr = ierr )
       !
    CASE ( 'noncolin' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, noncolin, ierr = ierr )
       !
    CASE ( 'ecfixed' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ecfixed, ierr = ierr )
       !
    CASE ( 'qcutz' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, qcutz, ierr = ierr )
       !
    CASE ( 'q2sigma' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, q2sigma, ierr = ierr )
       !
    CASE ( 'exchange' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       exchange = clean_str(tmpstr)
       !
       ! ... new variables that take the place of input_dft
    CASE ( 'exchange_grad_corr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       exchange_grad_corr = clean_str(tmpstr)
       !
    CASE ( 'correlation' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       correlation = clean_str(tmpstr)
       !
    CASE ( 'correlation_grad_corr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       correlation_grad_corr = clean_str(tmpstr)
       !
    CASE ( 'xc_specials' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       xc_specials = clean_str(tmpstr)
       ! ... up to here
       !
    CASE ( 'lda_plus_u' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, lda_plus_u, ierr = ierr )
       !
    CASE ( 'Hubbard_alpha' ) 
       IF (ntyp > 0) THEN
          ALLOCATE(tmparray(ntyp))
          CALL iotk_scan_dat_inside( xmlinputunit, tmparray, ierr = ierr )
          Hubbard_alpha(1:ntyp) = tmparray(1:ntyp)
          DEALLOCATE(tmparray)
       ELSE
          ierr = 1
       END IF
       !
    CASE ( 'Hubbard_U' ) 
       IF (ntyp > 0) THEN
          ALLOCATE(tmparray(ntyp))
          CALL iotk_scan_dat_inside( xmlinputunit, tmparray, ierr = ierr )
          Hubbard_U(1:ntyp) = tmparray(1:ntyp)
          DEALLOCATE(tmparray)
       ELSE
          ierr = 1
       END IF
       !
    CASE ( 'starting_ns_eigenvalue' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, starting_ns_eigenvalue, ierr = ierr )
       !
    CASE ( 'U_projection_type' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       U_projection_type = clean_str(tmpstr)
       !
    CASE ( 'edir' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, edir, ierr = ierr )
       !
    CASE ( 'emaxpos' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, emaxpos, ierr = ierr )
       !
    CASE ( 'eopreg' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, eopreg, ierr = ierr )
       !
    CASE ( 'eamp' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, eamp, ierr = ierr )
       !
    CASE ( 'angle1' )
       IF (ntyp > 0) THEN
          ALLOCATE(tmparray(ntyp))
          CALL iotk_scan_dat_inside( xmlinputunit, tmparray, ierr = ierr )
          angle1(1:ntyp) = tmparray(1:ntyp)
          DEALLOCATE(tmparray)
       ELSE
          ierr = 1
       END IF
       !
    CASE ( 'angle2' ) 
       IF (ntyp > 0) THEN
          ALLOCATE(tmparray(ntyp))
          CALL iotk_scan_dat_inside( xmlinputunit, tmparray, ierr = ierr )
          angle2(1:ntyp) = tmparray(1:ntyp)
          DEALLOCATE(tmparray)
       ELSE
          ierr = 1
       END IF
       !
    CASE ( 'constrained_magnetization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       constrained_magnetization = clean_str(tmpstr)
       !
    CASE ( 'fixed_magnetization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, fixed_magnetization, ierr = ierr )
       !
    CASE ( 'lambda' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, lambda, ierr = ierr )
       !
    CASE ( 'report' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, report, ierr = ierr )
       !
    CASE ( 'lspinorb' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, lspinorb, ierr = ierr )
       !
    CASE ( 'assume_isolated' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       assume_isolated = clean_str(tmpstr)
       !
    CASE ( 'london' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, london, ierr = ierr )
       !
    CASE ( 'london_s6' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, london_s6, ierr = ierr )
       !
    CASE ( 'london_rcut' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, london_rcut, ierr = ierr )
       !
    CASE ( 'cell_dynamics' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       cell_dynamics = clean_str(tmpstr)
       !
    CASE ( 'press' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, press, ierr = ierr )
       !
    CASE ( 'wmass' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, wmass, ierr = ierr )
       !
    CASE ( 'cell_factor' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, cell_factor, ierr = ierr )
       !
    CASE ( 'press_conv_thr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, press_conv_thr, ierr = ierr )
       !
    CASE ( 'cell_dofree' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       cell_dofree = clean_str(tmpstr)
       !
    CASE ( 'ecutcoarse' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ecutcoarse, ierr = ierr )
       !
    CASE ( 'mixing_charge_compensation' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, mixing_charge_compensation, ierr = ierr )
       !
    CASE ( 'n_charge_compensation' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, n_charge_compensation, ierr = ierr )
       !
    CASE ( 'comp_thr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, comp_thr, ierr = ierr )
       !
    CASE ( 'nlev' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nlev, ierr = ierr )
       !
    CASE ( 'electron_maxstep' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, electron_maxstep, ierr = ierr )
       !
    CASE ( 'conv_thr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, conv_thr, ierr = ierr )
       !
    CASE ( 'mixing_mode' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       mixing_mode = clean_str(tmpstr)
       !
    CASE ( 'mixing_beta' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, mixing_beta, ierr = ierr )
          !
    CASE ( 'mixing_ndim' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, mixing_ndim, ierr = ierr )
       !
    CASE ( 'mixing_fixed_ns' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, mixing_fixed_ns, ierr = ierr )
       !
    CASE ( 'diagonalization' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       diagonalization = clean_str(tmpstr)
    CASE ( 'diago_thr_init' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, diago_thr_init, ierr = ierr )
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
    CASE ( 'efield' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, efield, ierr = ierr )
       !
    CASE ( 'efield_cart' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, efield_cart, ierr = ierr )
       !
    CASE ( 'startingpot' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, startingpot, ierr = ierr )
       !
    CASE ( 'startingwfc' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, startingwfc, ierr = ierr )
       !
    CASE ( 'tqr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tqr, ierr = ierr )
       !
    CASE ( 'ion_dynamics' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       ion_dynamics = clean_str(tmpstr)
       !
    CASE ( 'ion_positions' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       ion_positions = clean_str(tmpstr)
       !
    CASE ( 'phase_space' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       phase_space = clean_str(tmpstr)
       !
    CASE ( 'pot_extrapolation' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       pot_extrapolation = clean_str(tmpstr)
       !
    CASE ( 'wfc_extrapolation' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       wfc_extrapolation = clean_str(tmpstr)
       !
    CASE ( 'remove_rigid_rot' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, remove_rigid_rot, ierr = ierr )
       !
    CASE ( 'ion_temperature' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       ion_temperature = clean_str(tmpstr)
       !
    CASE ( 'tempw' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tempw, ierr = ierr )
       !
    CASE ( 'tolp' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tolp, ierr = ierr )
       !
    CASE ( 'delta_t' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, delta_t, ierr = ierr )
       !
    CASE ( 'nraise' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, nraise, ierr = ierr )
       !
    CASE ( 'refold_pos' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, refold_pos, ierr = ierr )
       !
    CASE ( 'upscale' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, upscale, ierr = ierr )
       !
    CASE ( 'bfgs_ndim' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, bfgs_ndim, ierr = ierr )
       !
    CASE ( 'trust_radius_max' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, trust_radius_max, ierr = ierr )
       !
    CASE ( 'trust_radius_min' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, trust_radius_min, ierr = ierr )
       !
    CASE ( 'trust_radius_ini' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, trust_radius_ini, ierr = ierr )
       !
    CASE ( 'w_1' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, w_1, ierr = ierr )
       !
    CASE ( 'w_2' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, w_2, ierr = ierr )
       !
    CASE ( 'opt_scheme' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       opt_scheme = clean_str(tmpstr)
       !
    CASE ( 'CI_scheme' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, tmpstr, ierr = ierr )
       CI_scheme = clean_str(tmpstr)
       !
    CASE ( 'first_last_opt' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, first_last_opt, ierr = ierr )
       !
    CASE ( 'temp_req' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, temp_req, ierr = ierr )
       !
    CASE ( 'ds' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, ds, ierr = ierr )
       !
    CASE ( 'k_max' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, k_max, ierr = ierr )
       !
    CASE ( 'k_min' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, k_min, ierr = ierr )
       !
    CASE ( 'path_thr' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, path_thr, ierr = ierr )
       !
    CASE ( 'use_masses' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, use_masses, ierr = ierr )
       !
    CASE ( 'use_freezing' ) 
       CALL iotk_scan_dat_inside( xmlinputunit, use_freezing, ierr = ierr )
       !
       !
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
