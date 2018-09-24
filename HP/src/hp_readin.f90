!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE hp_readin()
  !----------------------------------------------------------------------------
  !
  !  This routine reads the input parameters for HP and reads
  !  the data produced by PWscf.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : meta_ionode, meta_ionode_id
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE input_parameters, ONLY : max_seconds
  USE io_files,         ONLY : tmp_dir, prefix, create_directory
  USE control_flags,    ONLY : iverbosity
  USE control_lr,       ONLY : ethr_nscf, lrpa
  USE ldaU_hp,          ONLY : postproc_only, conv_thr_chi, thresh_init,               &
                               at_equiv_criterium, skip_atom, skip_type,               &
                               merge_type, background, collect_chi, tmp_dir_hp,        &
                               perturb_only_atom, sum_pertq, determine_num_pert_only,  &
                               skip_equivalence_q, tmp_dir_save, niter_max,            &
                               disable_type_analysis, docc_thr, num_neigh, lmin, rmax, &
                               nmix, nq1, nq2, nq3, alpha_mix, start_q, last_q, maxter
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER             :: ios, &    ! integer variable for I/O control
                         iter      ! counter on iterations
  CHARACTER (LEN=256) :: outdir
  CHARACTER(LEN=6)    :: int_to_char
  !
  NAMELIST / INPUTHP / prefix, outdir, nq1, nq2, nq3, skip_equivalence_q,             &
                         conv_thr_chi, skip_atom, skip_type, merge_type, iverbosity,  &
                         background, thresh_init, at_equiv_criterium, postproc_only,  &
                         niter_max, alpha_mix, nmix, collect_chi, perturb_only_atom,  &
                         start_q, last_q, sum_pertq, ethr_nscf, num_neigh, lmin,      &
                         determine_num_pert_only, disable_type_analysis, docc_thr,    &
                         max_seconds, rmax
  !
  ! Note: meta_ionode is a single processor that reads the input
  !       Data read from input is subsequently broadcast to all processors
  !       from meta_ionode_id (using the default communicator world_comm)
  !
  IF (meta_ionode) CALL input_from_file ()
  !
  ! Set default values for variables in namelist
  !
  prefix             = 'pwscf'
  conv_thr_chi       = 1.D-5
  thresh_init        = 1.D-14
  ethr_nscf          = 1.D-11
  docc_thr           = 5.D-5
  rmax               = 100.D0
  skip_atom(:)       = .FALSE.
  skip_type(:)       = .FALSE.
  perturb_only_atom(:)    = .FALSE.
  skip_equivalence_q      = .FALSE.
  determine_num_pert_only = .FALSE.
  disable_type_analysis   = .FALSE.
  merge_type(:)      = 0
  at_equiv_criterium = 1
  background         = 'no'
  postproc_only      = .FALSE.
  collect_chi        = .FALSE.
  sum_pertq          = .FALSE.
  num_neigh          = 6
  lmin               = 2
  nq1                = 1
  nq2                = 1
  nq3                = 1
  start_q            = 1
  last_q             = -1
  iverbosity         = 1
  niter_max          = 100     
  alpha_mix(:)       = 0.D0   
  alpha_mix(1)       = 0.3D0  
  nmix               = 4     
  max_seconds        = 1.E+7_DP
  lrpa               = .FALSE.   ! Needed in dv_of_drho
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  ! Reading the namelist inputhp
  !
  IF (meta_ionode) READ( 5, INPUTHP, IOSTAT = ios )
  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore( 'hp_readin', 'reading inputhp namelist', ABS( ios ) )
  !
  ! Setup tmp_dir
  !
  tmp_dir = trimcheck (outdir)
  !
  ! Broadcast input parameters over all processors
  !
  CALL hp_bcast_input()
  !
  ! Here we finished the reading of the input file.
  ! Now allocate space for pwscf variables, read and check them.
  !
  tmp_dir_save = tmp_dir
  tmp_dir_hp = TRIM (tmp_dir) // 'HP' //'/'
  CALL create_directory(tmp_dir_hp)
  !
  ! Read various data produced by PWscf.
  ! In particular, read the unperturbed occupation matrices
  ! via calling the routine read_rho.
  ! read_file calls init_at_1 which initializes tab_at.
  !
  CALL read_file()
  !
  ! Compute the trace of the occupation matrix and then compute
  ! the occupations and the magnetization
  !
  CALL hp_ns_trace()
  !
  IF (postproc_only) THEN
     tmp_dir = tmp_dir_hp
  ELSE
     tmp_dir = tmp_dir_save
  ENDIF
  !
  ! Make sure all the features used in the PWscf calculation 
  ! are actually supported by HP.
  !
  CALL input_sanity()
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE input_sanity()
  !-------------------------------------------------------------------------- 
  ! 
  ! This subroutine aims to gather all of the input sanity checks
  ! (features enabled in PWscf which are unsupported in HP).
  !
  USE klist,            ONLY : lgauss, ltetra, two_fermi_energies
  USE control_flags,    ONLY : gamma_only, tqr
  USE fixed_occ,        ONLY : tfixed_occ
  USE cellmd,           ONLY : lmovecell
  USE noncollin_module, ONLY : i_cons, noncolin
  USE control_flags,    ONLY : twfcollect
  USE mp_global,        ONLY : nproc_pool, nproc_pool_file, &
                               nproc_image_file, nproc_image, nproc_bgrp_file
  USE ldaU,             ONLY : lda_plus_u, U_projection, lda_plus_u_kind
                               !is_hubbard_back
  !
  IMPLICIT NONE
  !
  IF (conv_thr_chi <= 0.D0) CALL errore ('hp_readin', ' Wrong conv_thr_chi ',1)
  !
  IF (nq1.LE.0 .OR. nq2.LE.0 .OR. nq3.LE.0) &
       CALL errore('hp_readin','nq1, nq2, and nq3 must be greater than 0',1)
  !
  IF (start_q <= 0 ) CALL errore('hp_readin', ' Wrong start_q ',1)
  !
  IF (collect_chi .AND. ANY(perturb_only_atom(:))) &
     CALL errore ('hp_readin', 'collect_chi and perturb_only_atom are not allowed to be true together', 1)
  !
  !IF (ANY(is_hubbard_back(:))) &
  !   CALL errore ('hp_readin', 'Calculation of Hubbard parameters with the background is not implemented', 1)
  !
  IF (sum_pertq .AND. .NOT.ANY(perturb_only_atom(:))) &
     CALL errore ('hp_readin', 'sum_pertq can be set to .true. only if perturb_only_atom is .true. for some atom', 1)
  !
  IF (niter_max.LT.1 .OR. niter_max.GT.maxter) &
     CALL errore ('hp_readin', ' Wrong niter_max ', 1)
  !
  DO iter = 1, niter_max
     IF ( alpha_mix(iter).LT.0.D0 .OR. alpha_mix(iter).GT.1.D0 ) &
        CALL errore ('hp_readin', ' Wrong alpha_mix ', iter)
  ENDDO
  !
  IF (num_neigh.LT.1) CALL errore('hp_readin','Not allowed value of num_neigh',1)
  !
  IF (lmin.LT.0 .OR. lmin.GT.3) CALL errore('hp_readin','Not allowed value of lmin',1)
  !
  IF (nmix.LT.1 .OR. nmix.GT.5) &
     & CALL errore ('hp_readin', ' Wrong nmix ', 1) 
  !
  IF (ltetra) CALL errore ('hp_readin', 'HP with tetrahedra is not supported', 1)
  !
  IF (gamma_only) CALL errore('hp_readin',&
     & 'Cannot start from pw.x data file using Gamma-point tricks',1)
  !
  IF (.NOT.lda_plus_u) CALL errore('hp_readin',&
     & 'The HP code can be used only when lda_plus_u=.true.',1)
  !
  IF (lda_plus_u_kind/=0) CALL errore("hp_readin", &
     & ' The HP code supports only lda_plus_u_kind=0',1)
  !
  IF (noncolin) CALL errore('hp_readin','Noncolliner case is not supported',1)
  !
  IF (lmovecell) CALL errore('hp_readin','The HP code is not working after vc-relax',1)
  !
  IF (nproc_image /= nproc_image_file .and. .not. twfcollect) CALL errore('hp_readin', &
     & 'pw.x run with a different number of processors. Use wf_collect=.true.',1)
  !
  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  CALL errore('hp_readin', &
     & 'pw.x run with a different number of pools. Use wf_collect=.true.',1)
  !
  IF (nproc_bgrp_file /= nproc_pool_file) CALL errore('hp_readin', &
     & 'band parallelization is not implemented in HP',1)
  !
  IF (i_cons /= 0) CALL errore('hp_readin',&
     & 'The HP code with constrained magnetization is not yet available',1)
  !
  IF (two_fermi_energies .AND. (ltetra .OR. lgauss)) CALL errore('hp_readin', &
     & 'The HP code with two Fermi energies is not available for metals',1)
  !
  IF (tqr) CALL errore('hp_readin',&
     & 'The HP code with Q in real space is not supported',1)
  !
  IF (tfixed_occ) CALL errore('hp_readin', &
     & 'The HP code with arbitrary occupations not tested',1)
  !
  RETURN
  !
END SUBROUTINE input_sanity
  !
END SUBROUTINE hp_readin
