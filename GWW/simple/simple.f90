!-----------------------------------------------------------------------
program simple
  !!-----------------------------------------------------------------------
  !! 
  !! input:  namelist "&inputsimple", with variables
  !!   prefix       prefix of input files saved by program pwscf
  !!   outdir       temporary directory where files resides
  !!    
  use io_files,  ONLY : prefix, tmp_dir
  use io_files,  ONLY : psfile, pseudo_dir
  use io_global, ONLY : stdout, ionode, ionode_id
  USE mp_global, ONLY: mp_startup
  USE mp_pools, ONLY : kunit
  use mp_world, ONLY: mpime, world_comm 
  USE environment,   ONLY: environment_start
  USE mp, ONLY : mp_bcast
  use ldaU, ONLY : lda_plus_u
  use scf, only : vrs, vltot, v, kedtau
  USE fft_base,             ONLY : dfftp
  use pwcom, only :  nspin
  use uspp, ONLY : okvan
  use realus, ONLY : generate_qpointlist
  USE io_files, ONLY : seqopn
  USE wannier_gw, ONLY : num_nbndv 
  USE gvect, ONLY : ngm
  USE gvecs, ONLY : doublegrid
  USE constants, ONLY : BOHR_RADIUS_SI
  USE cell_base,        ONLY : alat
  USE input_simple, ONLY : num_val,  num_cond, s_bands, deallocate_simple,s_product,&
                       &l_truncated_coulomb,truncation_radius, allocate_simple, &
                       &nkpoints,numpw,l_debug,n_debug,w_type,epsm,lambdam, &
                       &prefix_small , nonlocal_commutator , calc_mode
  USE pwcom, ONLY : igk_k
  !
  IMPLICIT NONE
  CHARACTER(len=9) :: code = 'SIMPLE'
  INTEGER :: ios, kunittmp
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(len=200) :: pp_file
  LOGICAL :: uspp_spsi, ascii, single_file, raw
  LOGICAL :: exst
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(len=256) :: prefix_nc 
  INTEGER :: nglobals
  LOGICAL :: l_dipols
  CHARACTER(LEN=256) :: outdir
  !
  !
  NAMELIST /inputsimple/ prefix,outdir,num_nbndv,num_val,num_cond,s_bands, &
       & s_product,l_truncated_coulomb,truncation_radius,nkpoints,numpw,&
       & l_debug,n_debug,w_type,epsm,lambdam,prefix_small, nonlocal_commutator,&
       & calc_mode
  !
  CALL mp_startup ( )
  CALL environment_start ( code )
  !
  CALL start_clock('simple')
  !
  prefix='export'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  IF ( ionode ) THEN
     CALL input_from_file ( )
     READ(5,inputsimple,IOSTAT=ios)
     IF (ios /= 0) CALL errore ('SIMPLE', 'reading inputsimple namelist', ABS(ios) )
  ENDIF
  !
  tmp_dir = trimcheck( outdir )
  CALL mp_bcast( outdir, ionode_id, world_comm  )
  CALL mp_bcast( tmp_dir, ionode_id , world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm  )
  CALL mp_bcast( num_nbndv, ionode_id , world_comm )
  CALL mp_bcast(num_val, ionode_id, world_comm )
  CALL mp_bcast(num_cond, ionode_id, world_comm )
  CALL mp_bcast(s_bands, ionode_id, world_comm )
  CALL mp_bcast(s_product, ionode_id, world_comm)
  CALL mp_bcast(l_truncated_coulomb, ionode_id, world_comm)
  CALL mp_bcast(truncation_radius, ionode_id, world_comm )
  CALL mp_bcast(nkpoints, ionode_id, world_comm)
  CALL mp_bcast(numpw, ionode_id, world_comm)
  CALL mp_bcast(l_debug, ionode_id, world_comm)
  CALL mp_bcast(n_debug,  ionode_id, world_comm)
  CALL mp_bcast(w_type,  ionode_id, world_comm)
  CALL mp_bcast(epsm,  ionode_id, world_comm)
  CALL mp_bcast(lambdam,  ionode_id, world_comm)
  CALL mp_bcast(prefix_small,  ionode_id, world_comm)
  CALL mp_bcast(nonlocal_commutator,  ionode_id, world_comm)
  CALL mp_bcast(calc_mode,  ionode_id, world_comm)
  !
  IF (.not.l_truncated_coulomb .and. numpw/=0) numpw=numpw+1
  !
  CALL read_file
  !
  CALL openfile_school
  !
#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif
  !
  pp_file= ' '
  uspp_spsi = .FALSE.
  ascii = .FALSE.
  single_file = .FALSE.
  raw = .FALSE.
  !
  CALL read_export(pp_file,kunittmp,uspp_spsi, ascii, single_file, raw)
  !
  CALL summary()
  !
  CALL print_ks_energies()
  !
  CALL hinit0()
  !
  IF (lda_plus_u) THEN
    CALL init_ns()
  ENDIF
  !
  CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  !
  IF ( okvan) CALL generate_qpointlist()
  !
  CALL allocate_simple
  !
  CALL wfc_basis
  !
  IF (calc_mode==0) THEN  ! BSE calculation
      WRITE(stdout,*) '***********************************************************'
      WRITE(stdout,*) 'Preparing for subsequent BSE calculation (for simple_bse.x)'
      WRITE(stdout,*) '***********************************************************'
      CALL product_basis
      CALL v_product
      CALL epe
  ELSEIF(calc_mode==1) then ! IP calculation
      WRITE(stdout,*) '*********************************************************'
      WRITE(stdout,*) 'Preparing for subsequent IP calculation (for simple_ip.x)'
      WRITE(stdout,*) '*********************************************************'
      CALL khamiltonian
  ENDIF
  CALL deallocate_simple
  !
  CALL stop_clock('simple')
  CALL print_clock('optimal_basis')
  CALL print_clock('simple')
  !
  CALL stop_pp
  STOP
  !
END PROGRAM simple



