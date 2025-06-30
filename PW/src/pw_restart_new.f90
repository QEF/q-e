!
! Copyright (C) 2016-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE pw_restart_new
!----------------------------------------------------------------------------
  !
  ! ... New PWscf I/O using xml schema and (optionally) hdf5 binaries
  ! ... Parallel execution: the xml file is written by one processor only
  ! ... ("ionode_id"), read by all processors ;
  ! ... the wavefunction files are written / read by one processor per pool,
  ! ... collected on / distributed to all other processors in pool
  !
  USE kinds, ONLY: dp
  USE qes_types_module, ONLY : output_type, parallel_info_type, &
       general_info_type, input_type, gateInfo_type, dipoleOutput_type, &
       BerryPhaseOutput_type, hybrid_type, vdw_type, dftU_type, smearing_type, sawtoothEnergy_type
  USE qes_write_module, ONLY: qes_write
  USE qes_reset_module, ONLY: qes_reset 
  USE qes_bcast_module,ONLY : qes_bcast
  USE qexsd_module, ONLY: qexsd_xf, qexsd_openschema, qexsd_closeschema, &
       qexsd_readschema
  USE qexsd_input,  ONLY : qexsd_input_obj, qexsd_init_k_points_ibz, &
       qexsd_init_occupations, qexsd_init_smearing,qexsd_init_twochem
  USE qexsd_init,   ONLY: qexsd_init_convergence_info, qexsd_init_algorithmic_info,    & 
                          qexsd_init_atomic_species, qexsd_init_atomic_structure,      &
                          qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
                          qexsd_init_magnetization,qexsd_init_band_structure,          &
                          qexsd_init_dipole_info, qexsd_init_total_energy,             &
                          qexsd_init_vdw, qexsd_init_forces, qexsd_init_stress,        &
                          qexsd_init_outputElectricField, qexsd_init_outputPBC,        &
                          qexsd_init_gate_info, qexsd_init_hybrid,  qexsd_init_dftU,   &
                          qexsd_init_rism3d, qexsd_init_rismlaue, qexsd_init_esm,      &
                          qexsd_init_sawtooth_info, qexsd_occ_obj, qexsd_bp_obj,       & 
                          qexsd_start_k_obj
  USE qexsd_copy,      ONLY : qexsd_copy_parallel_info, &
       qexsd_copy_algorithmic_info, qexsd_copy_atomic_species, &
       qexsd_copy_atomic_structure, qexsd_copy_symmetry, &
       qexsd_copy_basis_set, qexsd_copy_dft, qexsd_copy_efield, &
       qexsd_copy_band_structure, qexsd_copy_magnetization, &
       qexsd_copy_kpoints, qexsd_copy_rism3d, qexsd_copy_rismlaue, &
       qexsd_copy_esm, qexsd_copy_twochem
  USE io_global, ONLY : ionode, ionode_id
  USE io_files,  ONLY : iunpun, xmlfile
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  REAL(DP),ALLOCATABLE       :: local_charges(:), local_mag(:,:) 
  PRIVATE
  PUBLIC :: pw_write_schema, write_collected_wfc
  PUBLIC :: read_xml_file, read_collected_wfc
  !
  CONTAINS
    !------------------------------------------------------------------------
    SUBROUTINE pw_write_schema( only_init, wf_collect )
      !------------------------------------------------------------------------
      !
      ! only_init  = T  write only variables that are known after the 
      !                 initial steps of initialization (e.g. structure)
      !            = F  write the complete xml file
      ! wf_collect = T  if final wavefunctions in portable format are written,
      !              F  if wavefunctions are either not written or are written
      !                 in binary non-portable form (for checkpointing)
      !                 NB: wavefunctions are not written here in any case
      !
      USE control_flags,        ONLY : istep, conv_elec, conv_ions, &
                                       lscf, scf_error, n_scf_steps, &
                                       tqr, tq_smoothing, tbeta_smoothing, &
                                       gamma_only, noinv, smallmem, &
                                       lforce=> tprnfor, tstress, &
                                       mbd_vdw, llondon, lxdm, ts_vdw
      USE constants,            ONLY : e2  
      USE realus,               ONLY : real_space
      USE uspp,                 ONLY : okvan
      USE paw_variables,        ONLY : okpaw
      USE uspp_param,           ONLY : upf
      USE cell_base,            ONLY : at, bg, alat, ibrav
      USE ions_base,            ONLY : nsp, ityp, atm, nat, tau, zv, amass
      USE noncollin_module,     ONLY : noncolin, npol, colin_mag
      USE io_files,             ONLY : psfile, molfile, pseudo_dir
      USE klist,                ONLY : nks, nkstot, xk, ngk, wk, &
                                       lgauss, ngauss, smearing, degauss, nelec, &
                                       two_fermi_energies, nelup, neldw, tot_charge, ltetra,&
                                       degauss_cond,nelec_cond
      USE start_k,              ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       nks_start, xk_start, wk_start
      USE gvect,                ONLY : ngm, ngm_g, g
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE gvecs,                ONLY : ngms_g, dual
      USE fft_base,             ONLY : dffts
      USE wvfct,                ONLY : npwx, et, wg, nbnd,nbnd_cond
      USE ener,                 ONLY : ef, ef_up, ef_dw, vtxc, etxc, ewld, etot, &
                                       ehart, eband, demet, edftd3, elondon, exdm,&
                                       esol, vsol,ef_cond
      USE tsvdw_module,         ONLY : EtsvdW
      USE libmbd_interface,     ONLY : EmbdvdW
      USE gvecw,                ONLY : ecutwfc
      USE fixed_occ,            ONLY : tfixed_occ, f_inp
      USE ktetra,               ONLY : tetra_type
      USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_projectors, &
                                       Hubbard_lmax, Hubbard_l, Hubbard_n, Hubbard_U, Hubbard_Um, Hubbard_Um_nc, & 
                                       Hubbard_J, Hubbard_n2, Hubbard_n3, Hubbard_l2, Hubbard_l3, Hubbard_V,     & 
                                       Hubbard_occ, Hubbard_alpha, Hubbard_alpha_back, nsg, order_um, Hubbard_J0,&
                                       Hubbard_beta, Hubbard_U2, is_hubbard, is_hubbard_back, backall, neighood, &
                                       nsg
      USE symm_base,            ONLY : nrot, nsym, invsym, s, ft, irt, &
                                       t_rev, sname, time_reversal, no_t_rev,&
                                       spacegroup
      USE lsda_mod,             ONLY : nspin, isk, lsda, starting_magnetization, magtot, & 
                                       absmag, local_charges, local_mag 
      USE noncollin_module,     ONLY : angle1, angle2, i_cons, mcons, bfield, &
                                       magtot_nc, lambda, domag, lspinorb
      USE funct,                ONLY : get_dft_short, get_nonlocc_name, dft_is_nonlocc
      
      USE scf,                  ONLY : rho
      USE force_mod,            ONLY : sumfor, force, sigma
      USE extfield,             ONLY : tefield, dipfield, edir, etotefield, &
                                       emaxpos, eopreg, eamp, el_dipole, ion_dipole,&
                                       gate, zgate, relaxz, block, block_1,&
                                       block_2, block_height, etotgatefield ! TB
      USE mp,                   ONLY : mp_sum
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE xc_lib,               ONLY : xclib_dft_is, get_gau_parameter, &
                                       get_screening_parameter, xclib_get_exx_fraction, exx_is_active
      USE exx_base,             ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut
      USE exx,                  ONLY : ecutfock, local_thr 
      USE london_module,        ONLY : scal6, lon_rcut, c6_i
      USE xdm_module,           ONLY : xdm_a1=>a1i, xdm_a2=>a2i
      USE tsvdw_module,         ONLY : vdw_isolated, vdw_econv_thr
      USE input_parameters,     ONLY : verbosity, calculation, ion_dynamics, starting_ns_eigenvalue, &
                                       vdw_corr, london, k_points, assume_isolated, &  
                                       dftd3_threebody, dftd3_version
      USE bp,                   ONLY : lelfield, lberry, el_pol, ion_pol
      !
      USE rap_point_group,      ONLY : elem, nelem, name_class
      USE rap_point_group_so,   ONLY : elem_so, nelem_so, name_class_so
      USE rap_point_group_is,   ONLY : sname_is
      USE bfgs_module,          ONLY : bfgs_get_n_iter
      USE fcp_module,           ONLY : lfcp, fcp_mu
      USE control_flags,        ONLY : ldftd3, do_makov_payne 
      USE Coul_cut_2D,          ONLY : do_cutoff_2D 
      USE esm,                  ONLY : do_comp_esm, esm_nfit, esm_w, esm_a, esm_bc, esm_efield  
      USE martyna_tuckerman,    ONLY : do_comp_mt 
      USE run_info,             ONLY : title
      !
      USE qexsd_module,         ONLY : qexsd_add_all_clocks 
      USE solvmol,              ONLY : nsolV, solVs
      USE rism3d_facade,        ONLY : lrism3d, ecutsolv, qsol, laue_nfit, expand_r, expand_l, &
                                       starting_r, starting_l, buffer_r, buffer_ru, buffer_rv, &
                                       buffer_l, buffer_lu, buffer_lv, both_hands, &
                                       ireference, rism3d_is_laue
      USE two_chem,             ONLY : twochem
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: only_init, wf_collect
      !
      ! Structure containing the output tag
      !
      TYPE(output_type)   :: output_obj
      !
      ! Loop counters and other internal auxiliary variables 
      !
      INTEGER    :: is, viz, na1, na2, nt1, m1, m2 
      !
      ! Auxiliary variables used to format arguments for xml file
      !
      CHARACTER(LEN=37)     :: dft_name
      CHARACTER(LEN=8)      :: smearing_loc
      CHARACTER(LEN=8), EXTERNAL :: schema_smearing
      CHARACTER(LEN=20)     :: occupations
      CHARACTER(LEN=20), EXTERNAL :: schema_occupations
      CHARACTER(LEN=20)     :: pbc_label
      INTEGER               :: npwx_g
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      LOGICAL                  :: dft_is_vdw 
      !
      ! Variables related to RISM output 
      !
      INTEGER               :: isolV
      CHARACTER(LEN=10), ALLOCATABLE :: slabel(:)
      REAL(DP), ALLOCATABLE :: solvrho1(:)
      REAL(DP), ALLOCATABLE :: solvrho2(:)
      ! 
      !
      INTEGER               :: iclass, isym, ielem
      CHARACTER(LEN=15)     :: symop_2_class(48)
      LOGICAL               :: scf_has_converged 
      LOGICAL               :: opt_conv_ispresent
      LOGICAL               :: empirical_vdw
      INTEGER               :: n_opt_steps
      INTEGER               :: n_scf_steps_
      REAL(DP),PARAMETER    :: Ry_to_Ha = 1 / e2 
      !
      ! Auxiliary structures containing optional variables
      !
      TYPE(gateInfo_type)          :: gate_info_opt 
      TYPE(dipoleOutput_type)      :: dipol_opt  
      TYPE(BerryPhaseOutput_type)  :: bp_obj_opt
      TYPE(hybrid_type)            :: hybrid_obj_opt
      TYPE(vdW_type)               :: vdw_obj_opt
      TYPE(dftU_type)              :: dftU_obj_opt
      TYPE(smearing_type)          :: smear_obj_opt 
      TYPE(sawtoothEnergy_type)    :: sawtooth_obj
      !
      ! Copies of optional variables (*_tg) and pointers to them (*_pt)
      ! Pointers are nullified to signal that there is no such variable
      !
      REAL(DP), TARGET :: homo_tg, lumo_tg, ef_tg, degauss_tg, demet_tg, &
           efield_corr_tg, potstat_corr_tg, gatefield_corr_tg, &
           dispersion_energy_tg, &
           london_rcut_tg, london_s6_tg, ts_vdw_econv_thr_tg, &
           xdm_a1_tg, xdm_a2_tg, ecutvcut_tg, scr_par_tg, loc_thr_tg  
      REAL(DP), POINTER :: homo_pt, lumo_pt, ef_pt, degauss_pt, demet_pt, &
           efield_corr_pt, potstat_corr_pt, gatefield_corr_pt, &
           vdw_term_pt, &
           london_rcut_pt, london_s6_pt, ts_vdw_econv_thr_pt, &
           xdm_a1_pt, xdm_a2_pt, ecutvcut_pt, scr_par_pt, loc_thr_pt
      INTEGER, TARGET  :: dftd3_version_tg
      INTEGER,POINTER  :: dftd3_version_pt
      LOGICAL, TARGET  :: dftd3_threebody_tg, ts_vdw_isolated_tg
      LOGICAL, POINTER :: dftd3_threebody_pt, ts_vdw_isolated_pt
      CHARACTER(LEN=20),TARGET   :: vdw_corr_tg, dft_nonlocc_tg
      CHARACTER(LEN=20),POINTER  :: vdw_corr_pt, non_local_term_pt
      ! orphaned pointers?
      REAL(DP), POINTER:: ts_thr_pt
      LOGICAL, POINTER :: ts_isol_pt
      !
      ! Arrays of optional variables
      ! If not allocated, there is no such variable
      !
      REAL(DP), ALLOCATABLE :: london_c6_(:), bp_el_pol(:), bp_ion_pol(:), &
           U_opt(:), Um_opt(:,:,:), J0_opt(:), alpha_opt(:), J_opt(:,:), beta_opt(:), &
           U2_opt(:), alpha_back_opt(:), ef_updw(:), nsg_(:,:,:,:)
      INTEGER,ALLOCATABLE :: n_opt(:), l_opt(:), l2_opt(:), l3_opt(:), n2_opt(:), n3_opt(:)
      LOGICAL, ALLOCATABLE :: backall_opt(:) 
      !
      !
      NULLIFY (homo_pt, lumo_pt, ef_pt, degauss_pt, demet_pt, &
           efield_corr_pt, potstat_corr_pt, gatefield_corr_pt, &
           vdw_term_pt, &
           london_rcut_pt, london_s6_pt, ts_vdw_econv_thr_pt, &
           xdm_a1_pt, xdm_a2_pt, ecutvcut_pt, scr_par_pt, loc_thr_pt )
      NULLIFY (dftd3_version_pt)
      NULLIFY (dftd3_threebody_pt, ts_vdw_isolated_pt)
      NULLIFY (vdw_corr_pt, non_local_term_pt)
      NULLIFY (ts_thr_pt, ts_isol_pt)
      !
      ! Global PW dimensions need to be properly computed, reducing across MPI tasks
      ! If local PW dimensions are not available, set to 0
      !
      ALLOCATE( ngk_g( nkstot ) )
      ngk_g(:) = 0
      IF ( ALLOCATED (ngk) ) THEN
         ngk_g(1:nks) = ngk(:)
         CALL mp_sum( ngk_g(1:nks), intra_bgrp_comm )
         CALL ipoolrecover( ngk_g, 1, nkstot, nks )
      END IF
      ! BEWARE: only the first pool has ngk_g for all k-points
      !
      ! ... compute the maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! XML descriptor
      ! 
      IF ( ionode ) THEN  
         !
         ! ... here we init the variables and finally write them to file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         output_obj%tagname="output"
         output_obj%lwrite = .TRUE.
         output_obj%lread  = .TRUE.
         !
!-------------------------------------------------------------------------------
! ... CONVERGENCE_INFO
!-------------------------------------------------------------------------------
         opt_conv_ispresent = .FALSE.
         SELECT CASE (TRIM( calculation )) 
            CASE ( "relax","vc-relax" )
                opt_conv_ispresent = .TRUE.
                IF (TRIM( ion_dynamics) == 'bfgs' ) THEN 
                    n_opt_steps = bfgs_get_n_iter('bfgs_iter ') 
                ELSE 
                    n_opt_steps = istep 
                END IF 
                scf_has_converged = conv_elec 
                n_scf_steps_ = n_scf_steps
            CASE ("nscf", "bands" )
                n_opt_steps = 0
                scf_has_converged = .FALSE. 
                n_scf_steps_ = 1
            CASE default
                n_opt_steps  = 0 
                scf_has_converged = conv_elec 
                n_scf_steps_ = n_scf_steps
         END SELECT
         ! 
         IF ( opt_conv_ispresent) THEN
             call qexsd_init_convergence_info(output_obj%convergence_info,   &
                        SCF_HAS_CONVERGED = scf_has_converged, &
                        N_SCF_STEPS = n_scf_steps_, SCF_ERROR=scf_error/e2,&
                        OPTIMIZATION_HAS_CONVERGED = conv_ions, &
                        N_OPT_STEPS = n_opt_steps, GRAD_NORM = sumfor, WF_COLLECTED=wf_collect)
         ELSE
             call qexsd_init_convergence_info(output_obj%convergence_info,   &
                        SCF_HAS_CONVERGED = scf_has_converged, &
                        N_SCF_STEPS = n_scf_steps_, SCF_ERROR=scf_error/e2, WF_COLLECTED=wf_collect)
         END IF
         output_obj%convergence_info_ispresent = .TRUE.
         !
            
!-------------------------------------------------------------------------------
! ... ALGORITHMIC_INFO
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_algorithmic_info(output_obj%algorithmic_info, &
              REAL_SPACE_BETA = real_space, REAL_SPACE_Q=tqr , USPP=okvan, PAW=okpaw)
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_SPECIES
!-------------------------------------------------------------------------------
         !
         ! while amass's are always present, starting_mag should not be passed
         ! for nspin==1 or contrained magnetization calculations
         !
         IF (noncolin) THEN
            CALL qexsd_init_atomic_species(output_obj%atomic_species, nsp, atm, psfile, &
                 amass, STARTING_MAGNETIZATION = starting_magnetization, &
                 ANGLE1=angle1, ANGLE2=angle2)
         ELSE IF (nspin==2) THEN 
            CALL qexsd_init_atomic_species(output_obj%atomic_species, nsp, atm, psfile, &
                 amass, STARTING_MAGNETIZATION=starting_magnetization)
         ELSE 
            CALL qexsd_init_atomic_species(output_obj%atomic_species, nsp, atm,psfile, &
                 amass)
         END IF
         output_obj%atomic_species%pseudo_dir = TRIM(pseudo_dir)
         output_obj%atomic_species%pseudo_dir_ispresent = .TRUE.
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_STRUCTURE
!-------------------------------------------------------------------------------
         !         
         CALL qexsd_init_atomic_structure(output_obj%atomic_structure, nsp, atm, ityp, &
              nat, alat*tau, alat, alat*at(:,1), alat*at(:,2), alat*at(:,3), ibrav, natomwfc)
         !
!-------------------------------------------------------------------------------
! ... SYMMETRIES
!-------------------------------------------------------------------------------
         !
         symop_2_class="not found"
         IF (TRIM (verbosity) == 'medium' .OR. TRIM(verbosity) == 'high') THEN
            IF ( noncolin )  THEN 
               symmetries_so_loop:DO isym = 1, nrot 
                  classes_so_loop:DO iclass = 1, 24
                     elements_so_loop:DO ielem=1, nelem_so(iclass)
                        IF ( elem_so(ielem,iclass) == isym) THEN 
                           symop_2_class(isym) = name_class_so(iclass)
                           EXIT symmetries_so_loop
                        END IF
                     END DO elements_so_loop 
                     END DO classes_so_loop
               END DO symmetries_so_loop
            !
            ELSE
               symmetries_loop:DO isym = 1, nrot
                  classes_loop:DO iclass = 1, 12
                     elements_loop:DO ielem=1, nelem (iclass)
                        ! if the time-reversal in collinear systems is not detected
                        IF (colin_mag <= 1) THEN 
                           IF ( elem(ielem,iclass) == isym) THEN
                              symop_2_class(isym) = name_class(iclass)
                              EXIT classes_loop
                           END IF
                        ! if the time-reversal in non-collinear systems is detected
                        ELSE ! IF (colin_mag == 2)
                           IF( sname_is(elem(ielem,iclass)) == sname(isym) ) THEN
                              symop_2_class(isym) = name_class(iclass)
                              EXIT classes_loop
                           END IF
                        END IF
                     END DO elements_loop
                  END DO classes_loop
               END DO symmetries_loop
            END IF
         END IF
         CALL qexsd_init_symmetries(output_obj%symmetries, spacegroup, &
              nsym, nrot, s, ft, sname, t_rev, nat, irt, &
              symop_2_class(1:nrot), verbosity, noncolin, colin_mag)
         output_obj%symmetries_ispresent=.TRUE. 
         !
!-------------------------------------------------------------------------------
! ... BASIS SET
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_basis_set(output_obj%basis_set, gamma_only, ecutwfc/e2, ecutwfc*dual/e2, &
              dfftp%nr1, dfftp%nr2, dfftp%nr3, dffts%nr1, dffts%nr2, dffts%nr3, &
              .FALSE., dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, ngms_g, npwx_g, &
              bg(:,1), bg(:,2), bg(:,3) )
         !
!-------------------------------------------------------------------------------
! ... DFT
!-------------------------------------------------------------------------------
         !
         IF (xclib_dft_is('hybrid') ) THEN 
            IF (get_screening_parameter() > 0.0_DP) THEN
               scr_par_tg = get_screening_parameter() 
               scr_par_pt => scr_par_tg
            END IF 
            IF (ecutvcut > 0.0_DP) THEN 
               ecutvcut_tg = ecutvcut/e2 
               ecutvcut_pt => ecutvcut_tg
            END IF 
            IF ( local_thr > 0._DP) THEN 
               loc_thr_tg = local_thr 
               loc_thr_pt => loc_thr_tg 
            END IF 
            CALL qexsd_init_hybrid(hybrid_obj_opt, DFT_IS_HYBRID = .TRUE., NQ1 = nq1 , NQ2 = nq2, NQ3 =nq3, & 
                                   ECUTFOCK = ecutfock/e2, &
                                   EXX_FRACTION = xclib_get_exx_fraction(), SCREENING_PARAMETER = scr_par_pt, &
                                   EXXDIV_TREATMENT = exxdiv_treatment, X_GAMMA_EXTRAPOLATION = x_gamma_extrapolation,&
                                   ECUTVCUT = ecutvcut_pt, LOCAL_THR = loc_thr_pt )
         ELSE 
            hybrid_obj_opt%lwrite=.false. 
         END IF 

         empirical_vdw = (llondon .OR. ldftd3 .OR. lxdm .OR. ts_vdw .OR. mbd_vdw )
         IF ( dft_is_nonlocc() .OR. empirical_vdw ) THEN 
            IF ( empirical_vdw) THEN
               vdw_corr_tg = TRIM(vdw_corr)
               vdw_corr_pt => vdw_corr_tg
               IF (llondon ) THEN
                  dispersion_energy_tg = elondon/e2
                  london_s6_tg = scal6
                  london_s6_pt => london_s6_tg
                  london_rcut_tg = lon_rcut
                  london_rcut_pt => london_rcut_tg
                  IF (ANY( c6_i(1:nsp) .NE. -1._DP )) THEN
                     ALLOCATE (london_c6_(nsp))
                     london_c6_(1:nsp) = c6_i(1:nsp)
                  END IF
                  !
               ELSE IF ( lxdm ) THEN
                  dispersion_energy_tg = exdm/e2
                  xdm_a1_tg = xdm_a1
                  xdm_a1_pt => xdm_a1_tg
                  xdm_a2_tg = xdm_a2
                  xdm_a2_pt => xdm_a2_tg
                  !
               ELSE IF ( ldftd3) THEN
                  dispersion_energy_tg = edftd3/e2
                  dftd3_version_tg = dftd3_version
                  dftd3_version_pt => dftd3_version_tg
                  dftd3_threebody_tg = dftd3_threebody
                  dftd3_threebody_pt => dftd3_threebody_tg
               ELSE IF ( ts_vdw ) THEN
                  dispersion_energy_tg = 2._DP * EtsvdW/e2
                  ts_vdw_isolated_tg = vdw_isolated
                  ts_vdw_isolated_pt => ts_vdw_isolated_tg
                  ts_vdw_econv_thr_tg = vdw_econv_thr
                  ts_vdw_econv_thr_pt => ts_vdw_econv_thr_tg
               ELSE IF ( mbd_vdw ) THEN
                  dispersion_energy_tg = 2._DP * EmbdvdW/e2 - 2._DP * EtsvdW/e2 !avoiding double-counting
               END IF
               vdw_term_pt => dispersion_energy_tg
            ELSE
               vdw_corr_tg = 'none'
               vdw_corr_pt => vdw_corr_tg
               dft_nonlocc_tg = TRIM(get_nonlocc_name())
               non_local_term_pt => dft_nonlocc_tg
            END IF
            !
            CALL qexsd_init_vdw(vdw_obj_opt, non_local_term_pt, vdw_corr_pt, vdw_term_pt, &
                 ts_thr_pt, ts_isol_pt, london_s6_pt, LONDON_C6 = london_c6_, &
                 LONDON_RCUT =   london_rcut_pt, XDM_A1 = xdm_a1_pt, XDM_A2 = xdm_a2_pt,&
                 DFTD3_VERSION = dftd3_version_pt, DFTD3_THREEBODY = dftd3_threebody_pt)
            !
         ELSE 
            vdw_obj_opt%lwrite=.false. 
         END IF
         IF ( lda_plus_u ) THEN   
            CALL check_and_allocate_real(U_opt, Hubbard_U, fac = Ry_to_Ha)
            CALL check_and_allocate_real(J0_opt, Hubbard_J0 , fac = Ry_to_Ha) 
            CALL check_and_allocate_real(alpha_opt, Hubbard_alpha, fac = Ry_to_Ha) 
            CALL check_and_allocate_real(beta_opt, Hubbard_beta, fac = Ry_to_Ha) 
            CALL check_and_allocate_real(U2_opt, Hubbard_U2, fac = Ry_to_Ha )
            CALL check_and_allocate_real(alpha_back_opt, Hubbard_alpha_back, fac = Ry_to_Ha)
            CALL check_and_allocate_integer(n_opt, Hubbard_n)
            CALL check_and_allocate_integer(l_opt, Hubbard_l)
            CALL check_and_allocate_integer(n2_opt, Hubbard_n2)
            CALL check_and_allocate_integer(l2_opt, Hubbard_l2)
            CALL check_and_allocate_integer(l3_opt, Hubbard_l3)
            CALL check_and_allocate_integer(n3_opt, Hubbard_n3)
            CALL check_and_allocate_logical(backall_opt, backall)
            IF ( ANY(Hubbard_J(:,1:nsp) /= 0.0_DP)) THEN
               ALLOCATE (J_opt(3,nsp)) 
               J_opt(:, 1:nsp) = Hubbard_J(:, 1:nsp)  
            END IF
            IF (ANY(Hubbard_Um(1:2*Hubbard_lmax+1,1:MIN(nspin,2),1:nsp)/=0.0_dp)) THEN
              ALLOCATE (Um_opt(2*Hubbard_lmax+1,MIN(nspin,2),nsp))
              Um_opt(1:2*Hubbard_lmax+1,1:MIN(nspin,2),1:nsp) = &
                      Hubbard_Um(1:2*Hubbard_lmax+1,1:MIN(nspin,2),1:nsp) * Ry_to_Ha  
            ELSE IF (ANY(Hubbard_Um_nc(1:4*Hubbard_lmax+2,1:nsp)/=0.0_dp)) THEN 
              ALLOCATE(Um_opt(4*Hubbard_lmax+2,1,1:nsp))
              Um_opt(1:4*Hubbard_lmax+2,1,1:nsp) = Hubbard_Um_nc(1:4*Hubbard_lmax+2,1:nsp) * Ry_to_Ha 
            END IF 
            IF (lda_plus_u_kind==2) THEN
               ALLOCATE (nsg_(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
               nsg_ = 0.0d0
               DO na1 = 1, nat
                  nt1 = ityp(na1)
                  IF (is_hubbard(nt1)) THEN
                     DO viz = 1, neighood(na1)%num_neigh
                        na2 = neighood(na1)%neigh(viz)
                        IF (na2==na1) THEN
                           DO is = 1, nspin
                              DO m1 = 1, 2*Hubbard_l(nt1)+1
                                 DO m2 = 1, 2*Hubbard_l(nt1)+1
                                    nsg_(m1,m2,is,na1) = DBLE(nsg(m1,m2,viz,na1,is)) 
                                 ENDDO
                              ENDDO
                           ENDDO   
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
            !
            ! Currently rho%nsb is not written/read to/from XML 
            !
            CALL qexsd_init_dftU (dftU_obj_opt, NSP = nsp, PSD = upf(1:nsp)%psd, SPECIES = atm(1:nsp),                & 
                    ITYP = ityp(1:nat), IS_HUBBARD = is_hubbard, IS_HUBBARD_BACK = is_hubbard_back, BACKALL = backall,& 
                    HUBB_OCC = Hubbard_occ, HUBB_n2 = n2_opt, HUBB_L2 = l2_opt, HUBB_L3 = l3_opt, NONCOLIN = noncolin,& 
                    HUBB_N3 = n3_opt, LDA_PLUS_U_KIND = lda_plus_u_kind, U_PROJECTION_TYPE = Hubbard_projectors,      &
                    U =U_opt, Um = Um_opt, U2 = U2_opt, J0 = J0_opt, J = J_opt, n = n_opt, l = l_opt,                 &
                    Hubbard_V = Hubbard_V *Ry_to_Ha, alpha = alpha_opt, beta = beta_opt, alpha_back = alpha_back_opt, & 
                    starting_ns = starting_ns_eigenvalue, Hub_ns = rho%ns, Hub_ns_nc = rho%ns_nc, order_um = order_um, & 
                    Hub_nsg = nsg_)
            !
            IF (ALLOCATED(J_opt)) DEALLOCATE(J_opt)
            IF (ALLOCATED(nsg_))  DEALLOCATE(nsg_)
         ELSE 
           dftU_obj_opt%lwrite=.false. 
         END IF 
         dft_name = get_dft_short()
         !
         CALL qexsd_init_dft  (output_obj%dft, dft_name, hybrid_obj_opt, vdw_obj_opt, dftU_obj_opt)
         CALL qes_reset(hybrid_obj_opt) 
         CALL qes_reset(vdw_obj_opt) 
         CALL qes_reset(dftU_obj_opt) 
         !
!-------------------------------------------------------------------------------
! ... PERIODIC BOUNDARY CONDITIONS 
!-------------------------------------------------------------------------------
         !
         IF (ANY([do_makov_payne, do_comp_mt, do_comp_esm, do_cutoff_2D]))  THEN
            output_obj%boundary_conditions_ispresent=.TRUE.
            IF (do_makov_payne) THEN 
               pbc_label = 'makov_payne' 
            ELSE IF ( do_comp_mt) THEN 
               pbc_label = 'martyna_tuckerman' 
            ELSE IF ( do_comp_esm) THEN 
               pbc_label = 'esm' 
            ELSE IF ( do_cutoff_2D) THEN 
               pbc_label = '2D'
            ELSE 
               CALL errore ('pw_restart_new.f90: ', 'internal error line 470', 1) 
            END IF 
            CALL qexsd_init_outputPBC(output_obj%boundary_conditions, TRIM(pbc_label) )  
            IF (do_comp_esm) THEN
               CALL qexsd_init_esm(output_obj%boundary_conditions%esm, esm_bc, esm_nfit, esm_w, esm_efield, esm_a)
               output_obj%boundary_conditions%esm_ispresent = .TRUE. 
            END IF  
         ENDIF
         !
!-------------------------------------------------------------------------------
! ... MAGNETIZATION
!-------------------------------------------------------------------------------
         !
         output_obj%magnetization_ispresent = .TRUE.  
         IF (noncolin) THEN
           CALL qexsd_init_magnetization(output_obj%magnetization, lsda, noncolin, lspinorb, TOTAL_MAG_NC = magtot_nc,&
             ABSOLUTE_MAG = absmag, ATM = upf(1:nsp)%psd, ITYP = ityp, DO_MAGNETIZATION = domag, & 
             SITE_MAG = local_mag, SITE_CHARGES = local_charges )
         ELSE IF (lsda) THEN 
           CALL qexsd_init_magnetization(output_obj%magnetization, lsda, noncolin, lspinorb, TOTAL_MAG = magtot, &
                ABSOLUTE_MAG = absmag, ATM = upf(1:nsp)%psd, ITYP = ityp, SITE_MAG_POL = local_mag, & 
                SITE_CHARGES = local_charges) 
         ELSE 
           CALL qexsd_init_magnetization(output_obj%magnetization, lsda, noncolin, lspinorb, ABSOLUTE_MAG = 0._DP, &
                ATM = upf(1:nsp)%psd, ITYP = ityp ) 
         END IF 
         !

!--------------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------------
         !
         ! skip if not yet computed
         !
         IF ( only_init ) GO TO 10
         !
         IF ( .NOT. ( lgauss .OR. ltetra )) THEN 
            CALL get_homo_lumo( homo_tg, lumo_tg)
            homo_tg = homo_tg/e2
            homo_pt => homo_tg
            IF ( lumo_tg .LT. 1.d+6 ) THEN
                lumo_tg = lumo_tg/e2
                lumo_pt => lumo_tg
            END IF
         END IF
         IF (nks_start == 0 .AND. nk1*nk2*nk3 > 0 ) THEN 
            CALL qexsd_init_k_points_ibz(qexsd_start_k_obj, "automatic", calculation, &
                 nk1, nk2, nk3, k1, k2, k3, nks_start, alat, at(:,1), .TRUE.)
         ELSE
            CALL qexsd_init_k_points_ibz(qexsd_start_k_obj, k_points, calculation, &
                 nk1, nk2, nk3, k1, k2, k3, nks_start, alat, at(:,1), .TRUE., xk_start, wk_start)
         END IF
         qexsd_start_k_obj%tagname = 'starting_kpoints'
         occupations = schema_occupations( lgauss, ltetra, tetra_type, &
                    tfixed_occ )
         IF ( TRIM (qexsd_input_obj%tagname) == 'input') THEN 
            qexsd_occ_obj = qexsd_input_obj%bands%occupations
         ELSE 
            CALL qexsd_init_occupations ( qexsd_occ_obj, occupations, nspin)
         END IF 
         qexsd_occ_obj%tagname = 'occupations_kind' 
         IF ( two_fermi_energies ) THEN
            ALLOCATE ( ef_updw (2) )
            IF (TRIM(occupations) == 'fixed') THEN  
               ef_updw(1)  = MAXVAL(et(INT(nelup),1:nkstot/2))/e2
               ef_updw(2)  = MAXVAL(et(INT(neldw),nkstot/2+1:nkstot))/e2 
            ELSE 
               ef_updw = [ef_up/e2, ef_dw/e2]
            END IF
         ELSE
            ! The Fermi energy is written also for insulators because it can
            ! be useful for further postprocessing, especially of bands
            ! (for an insulator the Fermi energy is equal to the HOMO/VBMax)
            ef_tg = ef/e2
            ef_pt => ef_tg
         END IF

         IF ( lgauss ) THEN
            IF (TRIM(qexsd_input_obj%tagname) == 'input') THEN 
               smear_obj_opt = qexsd_input_obj%bands%smearing
               IF (twochem) THEN 
               END IF 
            ELSE
               smearing_loc = schema_smearing( smearing )
               CALL qexsd_init_smearing(smear_obj_opt, smearing_loc, degauss/e2)

            END IF  
         ELSE 
            smear_obj_opt%lwrite=.false.  
         END IF 
         !
         CALL qexsd_init_band_structure(  output_obj%band_structure,lsda,noncolin,lspinorb, nelec, natomwfc, &
                                 et, wg, nkstot, xk, ngk_g, wk, SMEARING = smear_obj_opt,  &
                                 STARTING_KPOINTS = qexsd_start_k_obj, OCCUPATIONS_KIND = qexsd_occ_obj, &
                                 WF_COLLECTED = wf_collect, NBND = nbnd, FERMI_ENERGY = ef_pt, EF_UPDW = ef_updw, &
                                 HOMO = homo_pt, LUMO = lumo_pt )
         ! 
         IF (lgauss)  CALL qes_reset (smear_obj_opt)
         CALL qes_reset (qexsd_start_k_obj)
         CALL qes_reset (qexsd_occ_obj)
         !
!------------------------------------------------------------------------------------------
!... TWO CHEM STUFF 
!------------------------------------------------------------------------------------------
        output_obj%two_chem_ispresent = twochem 
        IF (twochem) THEN 
          IF (TRIM(qexsd_input_obj%tagname) == 'input') THEN 
             output_obj%two_chem_ispresent = .TRUE. 
             output_obj%two_chem = qexsd_input_obj%twoch_ 
             output_obj%two_chem%ef_cond_ispresent=.TRUE. 
             output_obj%two_chem%tagname="two_chem" 
             output_obj%two_chem%ef_cond = ef_cond 
          ELSE 
             CALL qexsd_init_twochem(output_obj%two_chem,"two_chem", twochem, nbnd_cond, &
                                                    nelec_cond, degauss_cond, ef_cond)             
          END IF 
        END IF 
!-------------------------------------------------------------------------------------------
! ... TOTAL ENERGY
!-------------------------------------------------------------------------------------------
         !
         IF ( degauss > 0.0d0 ) THEN
            !
            degauss_tg = degauss/e2
            degauss_pt => degauss_tg
            !
            demet_tg = demet/e2
            demet_pt => demet_tg
         END IF
         IF ( tefield ) THEN 
            efield_corr_tg =  etotefield/e2
            efield_corr_pt => efield_corr_tg
         END IF
         IF (lfcp ) THEN
            potstat_corr_tg = fcp_mu * tot_charge / e2
            potstat_corr_pt => potstat_corr_tg
            output_obj%FCP_tot_charge_ispresent = .TRUE.
            output_obj%FCP_tot_charge = tot_charge
            output_obj%FCP_force_ispresent = .TRUE.
            output_obj%FCP_force = (fcp_mu - ef) / e2
         END IF
         IF ( gate) THEN
            gatefield_corr_tg = etotgatefield/e2
            gatefield_corr_pt => gatefield_corr_tg
         END IF

         IF ( lrism3d) THEN 
            CALL qexsd_init_total_energy(output_obj%total_energy, etot/e2, eband/e2, ehart/e2, vtxc/e2, & 
                                       etxc/e2, ewld/e2, degauss_pt, demet_pt, efield_corr_pt, potstat_corr_pt, &
                                       gatefield_corr_pt, DISPERSION_CONTRIBUTION = vdw_term_pt, ESOL = esol/e2, & 
                                       VSOL = vsol/e2 )
         ELSE 
            CALL  qexsd_init_total_energy(output_obj%total_energy, etot/e2, eband/e2, ehart/e2, vtxc/e2, &
                                       etxc/e2, ewld/e2, degauss_pt, demet_pt, efield_corr_pt, potstat_corr_pt,&
                                       gatefield_corr_pt, DISPERSION_CONTRIBUTION = vdw_term_pt) 
         END IF 
         !
!---------------------------------------------------------------------------------------------
! ... FORCES
!----------------------------------------------------------------------------------------------
         !
         IF ( lforce .and. conv_elec ) THEN 
            output_obj%forces_ispresent = .TRUE.
            CALL qexsd_init_forces(output_obj%forces,nat,force,lforce)
         ELSE 
            output_obj%forces_ispresent = .FALSE.
            output_obj%forces%lwrite = .FALSE.  
         END IF 
         !
!------------------------------------------------------------------------------------------------
! ... STRESS 
!------------------------------------------------------------------------------------------------
         IF ( tstress .and. conv_elec ) THEN
            output_obj%stress_ispresent=.TRUE.
            CALL qexsd_init_stress(output_obj%stress, sigma, tstress ) 
         ELSE 
            output_obj%stress_ispresent=.FALSE.
            output_obj%stress%lwrite=.FALSE.
         END IF
!-------------------------------------------------------------------------------------------------
! ... ELECTRIC FIELD
!-------------------------------------------------------------------------------------------------
         output_obj%electric_field_ispresent = ( gate .OR. lelfield .OR. lberry .OR. tefield ) 

         IF ( gate ) THEN 
            CALL qexsd_init_gate_info(gate_info_opt,"gateInfo", etotgatefield/e2, zgate, nelec, &
                   alat, at, bg, zv, ityp)
         ELSE 
            gate_info_opt%lwrite=.false. 
         END IF
         IF ( lelfield ) THEN
            ALLOCATE (bp_el_pol(2), bp_ion_pol(3) )
            bp_el_pol = el_pol 
            bp_ion_pol(1:3) = ion_pol(1:3)
         END IF
         IF (tefield) THEN 
           CALL qexsd_init_sawtooth_info(sawtooth_obj, efield_corr_tg, edir, eamp, emaxpos, eopreg)  
           IF (dipfield) THEN 
              CALL qexsd_init_dipole_info(dipol_opt, el_dipole, ion_dipole, edir, eamp, &
                                  emaxpos, eopreg )  
           ELSE 
             dipol_opt%lwrite=.false. 
           END IF
         ELSE 
           sawtooth_obj%lwrite = .false.
         END IF 
         qexsd_bp_obj%lwrite= lberry 
         IF (output_obj%electric_field_ispresent) &
            CALL qexsd_init_outputElectricField(output_obj%electric_field, lelfield, tefield, dipfield, &
                 lberry, BP_OBJ = qexsd_bp_obj, EL_POL = bp_el_pol, ION_POL = bp_ion_pol,          &
                 GATEINFO = gate_info_opt, SAWTOOTH_OBJ=sawtooth_obj, DIPOLE_OBJ =  dipol_opt)
         !
         CALL qes_reset (gate_info_opt)
         CALL qes_reset (dipol_opt)
         CALL qes_reset( sawtooth_obj) 


!------------------------------------------------------------------------------------------------
! ... 3D-RISM
!------------------------------------------------------------------------------------------------
         !
         IF ( lrism3d ) THEN
            !
            IF ( nsolV > 0 ) THEN
               ALLOCATE( slabel( nsolV ) )
               ALLOCATE( solvrho1( nsolV ) )
               ALLOCATE( solvrho2( nsolV ) )
               DO isolV = 1, nsolV
                  slabel(isolV)   = solVs(isolV)%name
                  solvrho1(isolV) = solVs(isolV)%density
                  solvrho2(isolV) = solVs(isolV)%subdensity
               END DO
            ELSE
               ALLOCATE( solvrho1( 1 ) )
               ALLOCATE( solvrho2( 1 ) )
            END IF
            !
            output_obj%rism3d_ispresent = .TRUE.
            CALL qexsd_init_rism3d(output_obj%rism3d, nsolV, slabel, molfile, solvrho1, solvrho2, ecutsolv/e2)
            !
            DEALLOCATE( slabel )
            DEALLOCATE( solvrho1 )
            DEALLOCATE( solvrho2 )
            !
            IF ( rism3d_is_laue() ) THEN
               output_obj%rismlaue_ispresent = .TRUE.
               CALL qexsd_init_rismlaue(output_obj%rismlaue, both_hands, laue_nfit, ireference, qsol, &
                                        starting_r, expand_r, buffer_r, buffer_ru, buffer_rv, &
                                        starting_l, expand_l, buffer_l, buffer_lu, buffer_lv)
               !
            ELSE
               output_obj%rismlaue_ispresent = .FALSE.
               output_obj%rismlaue%lwrite    = .FALSE.
            END IF
            !
         ELSE
            output_obj%rism3d_ispresent = .FALSE.
            output_obj%rism3d%lwrite    = .FALSE.
            !
            output_obj%rismlaue_ispresent = .FALSE.
            output_obj%rismlaue%lwrite    = .FALSE.
         END IF
         !
!-------------------------------------------------------------------------------
! ... CLOCKS
         CALL qexsd_add_all_clocks()
!-------------------------------------------------------------------------------
! ... ACTUAL WRITING
!-------------------------------------------------------------------------------
 10      CONTINUE
         !
         CALL qexsd_openschema( xmlfile(), iunpun, 'PWSCF', title )
         CALL qes_write (qexsd_xf,output_obj)
         CALL qes_reset (output_obj) 
         CALL qexsd_closeschema()
         !
!-------------------------------------------------------------------------------
         !
      END IF
      DEALLOCATE (ngk_g)
      !
      RETURN
       !
    CONTAINS
       SUBROUTINE check_and_allocate_real(alloc, mydata, fac)
          IMPLICIT NONE
          REAL(DP),ALLOCATABLE  :: alloc(:) 
          REAL(DP)              :: mydata(:)
          REAL(DP), OPTIONAL    :: fac   
          IF ( ANY(mydata(1:nsp) /= 0.0_DP)) THEN 
             ALLOCATE(alloc(nsp)) 
             alloc(1:nsp) = mydata(1:nsp) 
             IF (PRESENT(fac))  alloc = alloc * fac 
          END IF 
          RETURN
       END SUBROUTINE check_and_allocate_real 
       !
       SUBROUTINE check_and_allocate_integer(alloc, mydata)
          IMPLICIT NONE
          INTEGER,ALLOCATABLE  :: alloc(:)
          INTEGER              :: mydata(:)
          IF ( ANY(mydata(1:nsp) /= -1)) THEN
             ALLOCATE(alloc(nsp))
             alloc(1:nsp) = mydata(1:nsp)
          END IF
          RETURN
       END SUBROUTINE check_and_allocate_integer
       !
       SUBROUTINE check_and_allocate_logical(alloc, mydata)
          IMPLICIT NONE
          LOGICAL,ALLOCATABLE  :: alloc(:)
          LOGICAL              :: mydata(:)
          IF ( ANY(mydata(1:nsp))) THEN
             ALLOCATE(alloc(nsp))
             alloc(1:nsp) = mydata(1:nsp)
          END IF
          RETURN
       END SUBROUTINE check_and_allocate_logical
       !
    END SUBROUTINE pw_write_schema
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_collected_wfc( )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_sum, mp_max
      USE io_base,              ONLY : write_wfc
      USE io_files,             ONLY : restart_dir, iunwfc, nwordwfc
      USE cell_base,            ONLY : tpiba, alat, bg
      USE control_flags,        ONLY : gamma_only, smallmem
      USE gvect,                ONLY : ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol
      USE buffers,              ONLY : get_buffer
      USE wavefunctions,        ONLY : evc
      USE exx,                  ONLY : xi, nbndproj
      USE xc_lib,               ONLY : exx_is_active
      USE klist,                ONLY : nks, nkstot, xk, ngk, igk_k
      USE gvect,                ONLY : ngm, g, mill
      USE fft_base,             ONLY : dfftp
      USE wvfct,                ONLY : npwx, et, wg, nbnd
      USE lsda_mod,             ONLY : nspin, isk, lsda
      USE mp_pools,             ONLY : intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : me_bgrp, root_bgrp, intra_bgrp_comm, &
                                       root_bgrp_id, my_bgrp_id
      USE clib_wrappers,        ONLY : f_mkdir_safe
      !
      IMPLICIT NONE
      !
      INTEGER               :: ios, ig, ngg, ipol, ispin
      INTEGER               :: ik, ik_g, ike, iks, npw_g
      INTEGER, EXTERNAL     :: global_kpoint_index
      INTEGER,  ALLOCATABLE :: ngk_g(:), mill_k(:,:)
      INTEGER,  ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:)
      CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
      CHARACTER(LEN=256)    :: dirname
      CHARACTER(LEN=320)    :: filename, filenameace
      !
      dirname = restart_dir ()
      !
      ! ... check that restart_dir exists on all processors that write
      ! ... wavefunctions; create one if restart_dir is not found. This
      ! ... is needed for k-point parallelization, in the case of non-parallel
      ! ... scratch file systems, that are not visible to all processors
      !
      IF ( my_bgrp_id == root_bgrp_id .AND. me_bgrp == root_bgrp ) THEN
         ios = f_mkdir_safe( TRIM(dirname) )
      END IF
      !
      ! ... write wavefunctions and k+G vectors
      !
      iks = global_kpoint_index (nkstot, 1)
      ike = iks + nks - 1
      !
      ! ... ngk_g: global number of k+G vectors
      !
      ALLOCATE( ngk_g( nks ) )
      ngk_g(1:nks) = ngk(1:nks)
      CALL mp_sum( ngk_g, intra_bgrp_comm)
      !
      ! ... The igk_l2g array yields the correspondence between the
      ! ... local k+G index and the global G index
      !
      ALLOCATE ( igk_l2g( npwx ) )
      !
      ! ... the igk_l2g_kdip local-to-global map yields the correspondence
      ! ... between the global order of k+G and the local index for k+G.
      !
      ALLOCATE ( igk_l2g_kdip( npwx ) )
      !
      ALLOCATE ( mill_k( 3, npwx ) )
      !
      k_points_loop: DO ik = 1, nks
         !
         ! ik_g is the index of k-point ik in the global list
         !
         ik_g = ik + iks - 1
         !
         ! ... Compute the igk_l2g array from previously computed arrays
         ! ... igk_k (k+G indices) and ig_l2g (local to global G index map)
         !
         igk_l2g = 0
         DO ig = 1, ngk (ik)
            igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
         END DO
         !
         ! ... npw_g is the maximum G vector index among all processors
         !
         npw_g = MAXVAL( igk_l2g(1:ngk(ik)) )
         CALL mp_max( npw_g, intra_pool_comm )
         !
         igk_l2g_kdip = 0
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik), igk_l2g, &
                              igk_l2g_kdip )
         !
         ! ... mill_k(:,i) contains Miller indices for (k+G)_i
         !
         DO ig = 1, ngk (ik)
            mill_k(:,ig) = mill(:,igk_k(ig,ik))
         END DO
         !
         ! ... read wavefunctions - do not read if already in memory (nsk==1)
         !
         IF ( nks > 1 ) CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
         !
         IF ( nspin == 2 ) THEN
            !
            ! ... LSDA: spin mapped to k-points, isk(ik) tracks up and down spin
            !
            ik_g = MOD ( ik_g-1, nkstot/2 ) + 1 
            ispin = isk(ik)
            filename = TRIM(dirname) // 'wfc' // updw(ispin) // &
                 & TRIM(int_to_char(ik_g))
            !
            if(exx_is_active()) filenameace = TRIM(dirname) // 'ace' // updw(ispin) // &
                 & TRIM(int_to_char(ik_g))
            !
         ELSE
            !
            ispin = 1
            filename = TRIM(dirname) // 'wfc' // TRIM(int_to_char(ik_g))
            !
            if(exx_is_active()) filenameace = TRIM(dirname) // 'ace' // TRIM(int_to_char(ik_g))
            !
         ENDIF
         !
         ! ... Only the first band group of each pool writes
         ! ... No warranty it works for more than one band group
         !
         IF ( my_bgrp_id == root_bgrp_id ) CALL write_wfc( iunpun, &
              filename, root_bgrp, intra_bgrp_comm, ik_g, tpiba*xk(:,ik), &
              ispin, nspin, evc, npw_g, gamma_only, nbnd, &
              igk_l2g_kdip(:), ngk(ik), tpiba*bg(:,1), tpiba*bg(:,2), &
              tpiba*bg(:,3), mill_k, 1.D0 )
         !
         IF ( (my_bgrp_id == root_bgrp_id) .and. exx_is_active() .and. allocated(xi)) then
              CALL write_wfc( iunpun, &
              filenameace, root_bgrp, intra_bgrp_comm, ik_g, tpiba*xk(:,ik), &
              ispin, nspin, xi(:,:,ik), npw_g, gamma_only, nbnd, &
              igk_l2g_kdip(:), ngk(ik), tpiba*bg(:,1), tpiba*bg(:,2), &
              tpiba*bg(:,3), mill_k, 1.D0 )
         END IF 
         !
      END DO k_points_loop
      !
      DEALLOCATE ( mill_k )
      DEALLOCATE ( igk_l2g_kdip )
      DEALLOCATE ( igk_l2g )
      DEALLOCATE ( ngk_g )
      !
      RETURN
      !
    END SUBROUTINE write_collected_wfc
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gk_l2gmap_kdip( npw_g, ngk_g, ngk, igk_l2g, igk_l2g_kdip, igwk )
      !-----------------------------------------------------------------------
      !
      ! ... This subroutine maps local G+k index to the global G vector index
      ! ... the mapping is used to collect wavefunctions subsets distributed
      ! ... across processors.
      ! ... This map is used to obtained the G+k grids related to each kpt
      !
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE mp,                   ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      ! ... Here the dummy variables
      !
      INTEGER, INTENT(IN)  :: npw_g, ngk_g, ngk
      INTEGER, INTENT(IN)  :: igk_l2g(ngk)
      INTEGER, INTENT(OUT) :: igk_l2g_kdip(ngk)
      INTEGER, OPTIONAL, INTENT(OUT) :: igwk(ngk_g)
      !
      INTEGER, ALLOCATABLE :: igwk_(:), itmp(:), igwk_lup(:)
      INTEGER              :: ig, ig_, ngg
      !
      !
      ALLOCATE( itmp( npw_g ) )
      ALLOCATE( igwk_( ngk_g ) )
      !
      itmp(:)  = 0
      igwk_(:) = 0
      !
      DO ig = 1, ngk
         itmp(igk_l2g(ig)) = igk_l2g(ig)
      END DO
      !
      CALL mp_sum( itmp, intra_bgrp_comm )
      !
      ngg = 0
      DO ig = 1, npw_g
         !
         IF ( itmp(ig) == ig ) THEN
            !
            ngg = ngg + 1
            igwk_(ngg) = ig
            !
         END IF
         !
      END DO
      !
      IF ( ngg /= ngk_g ) &
         CALL errore( 'gk_l2gmap_kdip', 'unexpected dimension in ngg', 1 )
      !
      IF ( PRESENT( igwk ) ) THEN
         !
         igwk(1:ngk_g) = igwk_(1:ngk_g)
         !
      END IF
      !
      ALLOCATE( igwk_lup( npw_g ) )
      !
!$omp parallel private(ig_, ig)
!$omp workshare
      igwk_lup = 0
!$omp end workshare
!$omp do
      DO ig_ = 1, ngk_g
         igwk_lup(igwk_(ig_)) = ig_
      END DO
!$omp end do
!$omp do
      DO ig = 1, ngk
         igk_l2g_kdip(ig) = igwk_lup(igk_l2g(ig))
      END DO
!$omp end do
!$omp end parallel
      !
      DEALLOCATE( igwk_lup )
      !
      DEALLOCATE( itmp, igwk_ )
      !
      RETURN
      !
    END SUBROUTINE gk_l2gmap_kdip
    !
    !--------------------------------------------------------------------------
    SUBROUTINE read_xml_file ( wfc_is_collected )
      !------------------------------------------------------------------------
      !
      ! ... This routine allocates space for all quantities already computed
      ! ... in the pwscf program and reads them from the data file.
      ! ... All quantities that are initialized in subroutine "setup" when
      ! ... starting from scratch should be initialized here when restarting
      !
      USE kinds,           ONLY : dp
      USE constants,       ONLY : e2  
      USE gvect,           ONLY : ngm_g, ecutrho
      USE gvecs,           ONLY : ngms_g, dual
      USE gvecw,           ONLY : ecutwfc
      USE fft_base,        ONLY : dfftp, dffts
      USE io_global,       ONLY : stdout
      USE io_files,        ONLY : psfile, pseudo_dir, pseudo_dir_cur, &
           restart_dir, molfile
      USE mp_global,       ONLY : nproc_file, nproc_pool_file, &
           nproc_image_file, ntask_groups_file, &
           nproc_bgrp_file, nproc_ortho_file
      USE ions_base,       ONLY : nat, nsp, ityp, amass, atm, tau, extfor
      USE cell_base,       ONLY : alat, at, bg, ibrav, celldm, omega
      USE fixed_occ,       ONLY : tfixed_occ
      USE force_mod,       ONLY : force
      USE klist,           ONLY : nks, nkstot, xk, wk, tot_magnetization, &
           nelec, nelup, neldw, smearing, degauss, ngauss, lgauss, ltetra,&
           two_fermi_energies,degauss_cond,nelec_cond
      USE ktetra,          ONLY : ntetra, tetra_type
      USE start_k,         ONLY : nks_start, xk_start, wk_start, &
           nk1, nk2, nk3, k1, k2, k3
      USE ener,            ONLY : ef, ef_up, ef_dw,ef_cond
      USE electrons_base,  ONLY : nupdwn, set_nelup_neldw
      USE wvfct,           ONLY : npwx, nbnd, et, wg,nbnd_cond
      USE extfield,        ONLY : forcefield, forcegate, tefield, dipfield, &
           edir, emaxpos, eopreg, eamp, el_dipole, ion_dipole, gate, zgate, &
           relaxz, block, block_1, block_2, block_height
      USE symm_base,       ONLY : nrot, nsym, invsym, s, ft, irt, t_rev, &
           sname, inverse_s, s_axis_to_cart, spacegroup, &
           time_reversal, no_t_rev, nosym, checkallsym
      USE ldaU,            ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_lmax, Hubbard_lmax_back, &
                                  Hubbard_n, Hubbard_l, Hubbard_n2, Hubbard_l2, Hubbard_n3, Hubbard_l3, backall, &
                                  Hubbard_U, Hubbard_U2, Hubbard_J, Hubbard_V, Hubbard_alpha, Hubbard_occ, &
                                  Hubbard_alpha_back, Hubbard_J0, Hubbard_beta, Hubbard_projectors, Hubbard_Um, & 
                                  Hubbard_Um_nc, apply_u, Hubbard_alpha_m, Hubbard_alpha_m_nc
      USE funct,           ONLY : enforce_input_dft, get_dft_short
      USE xc_lib,          ONLY : start_exx, exx_is_active,xclib_dft_is,      &
                                  set_screening_parameter, set_gau_parameter, &
                                  xclib_set_exx_fraction, stop_exx, start_exx  
      USE london_module,   ONLY : scal6, lon_rcut, in_C6
      USE tsvdw_module,    ONLY : vdw_isolated
      USE exx_base,        ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
           exxdiv_treatment, yukawa, ecutvcut
      USE exx,             ONLY : ecutfock, local_thr
      USE control_flags,   ONLY : noinv, gamma_only, tqr, llondon, ldftd3, &
           lxdm, ts_vdw, mbd_vdw, do_makov_payne 
      USE Coul_cut_2D,     ONLY : do_cutoff_2D
      USE esm,             ONLY : do_comp_esm, esm_bc, esm_nfit, esm_w, esm_efield, esm_a
      USE martyna_tuckerman,ONLY: do_comp_mt 
      USE noncollin_module,ONLY : noncolin, npol, angle1, angle2, bfield, &
              nspin_lsda, nspin_gga, nspin_mag, domag, lspinorb, colin_mag
      USE lsda_mod,        ONLY : nspin, isk, lsda, starting_magnetization,&
           current_spin
      USE realus,          ONLY : real_space
      USE basis,           ONLY : natomwfc
      USE uspp,            ONLY : okvan
      USE paw_variables,   ONLY : okpaw
      !
      USE solvmol,         ONLY : nsolV, solVs
      USE rism3d_facade,   ONLY : lrism3d, ecutsolv, qsol, laue_nfit, expand_r, expand_l, &
                                  starting_r, starting_l, buffer_r, buffer_ru, buffer_rv, &
                                  buffer_l, buffer_lu, buffer_lv, both_hands, ireference, &
                                  rism3d_set_laue
      !
      USE mp_images,       ONLY : intra_image_comm
      USE mp,              ONLY : mp_bcast
      USE dftd3_qe,        ONLY : dftd3_in, dftd3, dftd3_xc 
      USE dftd3_api,       ONLY : dftd3_init, dftd3_set_functional 
      USE tsvdw_module,    ONLY : vdw_econv_thr
      USE london_module,   ONLY : init_london
      USE xdm_module,      ONLY : init_xdm
      USE input_parameters,ONLY : verbosity, calculation, ion_dynamics, starting_ns_eigenvalue, &
                                       vdw_corr, london, k_points, assume_isolated, &  
                                       occupations, dftd3_threebody, dftd3_version
      USE two_chem,        ONLY : twochem
      !
      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: wfc_is_collected
      !
      INTEGER  :: i, is, ik, ierr, dum1,dum2,dum3
      LOGICAL  :: magnetic_sym, lvalid_input
      CHARACTER(LEN=37)  :: dft_name
      CHARACTER(LEN=256) ::dft_
      INTEGER           :: npwx_g, llmax, ntmax
      CHARACTER(LEN=320):: filename
      REAL(dp) :: exx_fraction, screening_parameter
      TYPE (output_type)        :: output_obj 
      TYPE (parallel_info_type) :: parinfo_obj
      TYPE (general_info_type ) :: geninfo_obj
      TYPE (input_type)         :: input_obj
      !
      !
      filename = xmlfile ( )
      !
      IF (ionode) CALL qexsd_readschema ( filename, &
           ierr, output_obj, parinfo_obj, geninfo_obj, input_obj)
      CALL mp_bcast(ierr, ionode_id, intra_image_comm)
      IF ( ierr > 0 ) THEN
         CALL errore ( 'read_xml_file', 'fatal error reading xml file', ierr ) 
      ELSE IF ( ierr < 0 ) THEN
         input_obj%tagname = "not_read"
         ! ierr = -1 means that input_obj was not read: do not broadcast it
      ELSE
         CALL qes_bcast(input_obj, ionode_id, intra_image_comm)
      END IF
      CALL qes_bcast(output_obj, ionode_id, intra_image_comm)
      CALL qes_bcast(parinfo_obj, ionode_id, intra_image_comm)
      CALL qes_bcast(geninfo_obj, ionode_id, intra_image_comm) 
      !
      ! ... Now read all needed variables from xml objects
      !
      wfc_is_collected = output_obj%convergence_info%wf_collected 
      lvalid_input = (TRIM(input_obj%tagname) == "input")
      !
      CALL qexsd_copy_parallel_info (parinfo_obj, nproc_file, &
           nproc_pool_file, nproc_image_file, ntask_groups_file, &
           nproc_bgrp_file, nproc_ortho_file)
      !
      pseudo_dir_cur = restart_dir ( )
      CALL qexsd_copy_atomic_species ( output_obj%atomic_species, &
           nsp, atm, amass, starting_magnetization, angle1, angle2, &
           psfile, pseudo_dir ) 
      IF ( pseudo_dir == ' ' ) pseudo_dir=pseudo_dir_cur
      !! Atomic structure section
      !! tau and ityp are allocated inside qexsd_copy_atomic_structure
      !
      CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp, &
           atm, nat, tau, ityp, alat, at(:,1), at(:,2), at(:,3), ibrav, natomwfc )
      !
      !! More initializations needed for atomic structure:
      !! bring atomic positions and crystal axis into "alat" units;
      !! recalculate celldm; compute cell volume, reciprocal lattice vectors
      !
      at = at / alat
      tau(:,1:nat) = tau(:,1:nat)/alat  
      CALL at2celldm (ibrav,alat,at(:,1),at(:,2),at(:,3),celldm)
      CALL volume (alat,at(:,1),at(:,2),at(:,3),omega)
      !!
      !! Basis set section
      CALL qexsd_copy_basis_set ( output_obj%basis_set, gamma_only, ecutwfc,&
           ecutrho, dffts%nr1,dffts%nr2,dffts%nr3, dfftp%nr1,dfftp%nr2,dfftp%nr3, &
           dum1,dum2,dum3, ngm_g, ngms_g, npwx_g, bg(:,1), bg(:,2), bg(:,3) )
      ecutwfc = ecutwfc*e2
      ecutrho = ecutrho*e2
      dual = ecutrho/ecutwfc
      ! FIXME: next line ensures exact consistency between reciprocal and
      ! direct lattice vectors, preventing weird phonon symmetry errors
      ! (due to lousy algorithms, extraordinarily sensitive to tiny errors)
      CALL recips ( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
      !!
      !! DFT section
      CALL qexsd_copy_dft ( output_obj%dft, nsp, atm, &
           dft_name, nq1, nq2, nq3, ecutfock, exx_fraction, screening_parameter, &
           exxdiv_treatment, x_gamma_extrapolation, ecutvcut, local_thr, &
           lda_plus_u, apply_u,lda_plus_u_kind, Hubbard_projectors, Hubbard_n, Hubbard_l, Hubbard_lmax, Hubbard_occ,&
           Hubbard_n2, Hubbard_l2, Hubbard_n3, Hubbard_l3, backall, Hubbard_lmax_back, Hubbard_alpha_back, &
           Hubbard_U, Hubbard_Um, Hubbard_U2, Hubbard_J0, Hubbard_alpha, Hubbard_alpha_m, Hubbard_beta, Hubbard_J, Hubbard_V, &
           vdw_corr, dftd3_version, dftd3_threebody, scal6, lon_rcut, vdw_isolated )
      Hubbard_alpha_back = Hubbard_alpha_back * e2 
      Hubbard_alpha      = Hubbard_alpha      * e2
      Hubbard_beta       = Hubbard_beta       * e2 
      Hubbard_U          = Hubbard_U          * e2 
      Hubbard_Um         = Hubbard_Um         * e2
      Hubbard_Um_nc      = 0.0_DP 
      Hubbard_alpha_m   = Hubbard_alpha_m     * e2 
      Hubbard_alpha_m_nc = 0.0_DP 
      IF (noncolin) THEN
        llmax = SIZE(Hubbard_Um_nc,1)/2   
        ntmax = SIZE(Hubbard_Um_nc,2) 
        Hubbard_Um_nc(:,:) = RESHAPE(Hubbard_Um(:,:,:),[2*llmax,ntmax]) 
        Hubbard_alpha_m_nc(:,:) = RESHAPE(Hubbard_alpha_m(:,:,:),[2*llmax, ntmax])   
        Hubbard_Um = 0.0_DP 
        Hubbard_alpha_m = 0.0_DP 
      END IF
      Hubbard_U2         = Hubbard_U2         * e2 
      Hubbard_V          = Hubbard_V          * e2 
      Hubbard_J0         = Hubbard_J0         * e2 
      Hubbard_J          = Hubbard_J          * e2  
      !! More DFT initializations

      CALL set_vdw_corr ( vdw_corr, llondon, ldftd3, ts_vdw, mbd_vdw, lxdm )
      !FIXME this maybe should be done directly in set_vdw_corr 
      CALL enforce_input_dft ( dft_name, .TRUE. )
      vdw_econv_thr   = input_obj%dft%vdW%ts_vdw_econv_thr
      IF ( lxdm )   CALL init_xdm ( )
      IF ( llondon) CALL init_london ( )
      IF ( ldftd3)  CALL dftd3_iosys ()
      IF ( xclib_dft_is('hybrid') ) THEN
         ecutvcut = ecutvcut*e2
         ecutfock = ecutfock*e2
         CALL xclib_set_exx_fraction( exx_fraction ) 
         CALL set_screening_parameter( screening_parameter )
         CALL start_exx ()
      END IF
      !! Band structure section
      !! et and wg are allocated inside qexsd_copy_band_structure
      CALL qexsd_copy_band_structure( output_obj%band_structure, lsda, &
           nkstot, isk, nbnd, nupdwn(1), nupdwn(2), nelec, xk, &
           wk, wg, ef, ef_up, ef_dw, et )
      ! convert to Ry
      ef = ef*e2
      ef_up = ef_up*e2
      ef_dw = ef_dw*e2
      two_fermi_energies = ( ef_up /= 0.0_dp ) .AND. ( ef_dw /= 0.0_dp )
      et(:,:) = et(:,:)*e2
      !
      ! ... until pools are activated, the local number of k-points nks
      ! ... should be equal to the global number nkstot - k-points are replicated
      !
      nks = nkstot
      !!
      !! Magnetization section
      CALL qexsd_copy_magnetization ( output_obj%magnetization, lsda, noncolin,&
           lspinorb, domag, tot_magnetization )
      !
      bfield = 0.d0
      CALL set_spin_vars( lsda, noncolin, domag, &
           npol, nspin, nspin_lsda, nspin_mag, nspin_gga, current_spin )
      !! Information for generating k-points and occupations
      CALL qexsd_copy_kpoints( output_obj%band_structure, &
           nks_start, xk_start, wk_start, nk1, nk2, nk3, k1, k2, k3, &
           occupations, smearing, degauss )
      degauss = degauss * e2 
      !
      CALL set_occupations( occupations, smearing, degauss, &
           tfixed_occ, ltetra, tetra_type, lgauss, ngauss )
      !! Information for twochem case
      IF (output_obj%two_chem_ispresent) THEN 
        CALL qexsd_copy_twochem(output_obj%two_chem, twochem, nbnd_cond, nelec_cond,degauss_cond,ef_cond)
      ELSE 
        twochem = .FALSE.
      END IF 
      !
      IF (ltetra) ntetra = 6* nk1 * nk2 * nk3 
      IF ( lsda ) &
           CALL set_nelup_neldw(tot_magnetization, nelec, nelup, neldw) 
      !! Symmetry section
      ALLOCATE ( irt(48,nat) )
      IF ( lvalid_input ) THEN 
         CALL qexsd_copy_symmetry ( output_obj%symmetries, &
              spacegroup, nsym, nrot, s, ft, sname, t_rev, invsym, irt, &
              noinv, nosym, no_t_rev, colin_mag, input_obj%symmetry_flags )
         IF (input_obj%electric_field_ispresent) & 
           CALL qexsd_copy_efield ( input_obj%electric_field, &
              tefield, dipfield, edir, emaxpos, eopreg, eamp, &
              gate, zgate, block, block_1, block_2, block_height, relaxz )
         
      ELSE 
         CALL qexsd_copy_symmetry ( output_obj%symmetries, &
              spacegroup, nsym, nrot, s, ft, sname, t_rev, invsym, irt, &
              noinv, nosym, no_t_rev ,colin_mag)
      ENDIF
      !! More initialization needed for symmetry
      magnetic_sym = noncolin .AND. domag
      time_reversal = (.NOT.magnetic_sym) .AND. (.NOT.noinv) 
      CALL inverse_s()
      CALL s_axis_to_cart()
      !! symmetry check - FIXME: must be done in a more consistent way 
      !! IF (nat > 0) CALL checkallsym( nat, tau, ityp)
      !! Algorithmic info
      IF (output_obj%boundary_conditions_ispresent) THEN 
         do_makov_payne = (output_obj%boundary_conditions%assume_isolated == "makov_payne")
         do_comp_mt     = (output_obj%boundary_conditions%assume_isolated == "martyna_tuckerman")
         do_comp_esm    = (output_obj%boundary_conditions%assume_isolated == "esm")
         do_cutoff_2D   = (output_obj%boundary_conditions%assume_isolated == "2D")
      ELSE
         do_makov_payne= .FALSE.
         do_comp_mt    = .FALSE.
         do_comp_esm   = .FALSE.
         do_cutoff_2D  = .FALSE.
      END IF
      IF (do_comp_esm)  CALL qexsd_copy_esm(output_obj%boundary_conditions, esm_bc, esm_nfit, esm_w, esm_efield, esm_a) 
      CALL qexsd_copy_algorithmic_info ( output_obj%algorithmic_info, &
           real_space, tqr, okvan, okpaw )
      !
      !! 3D-RISM
      IF ( output_obj%rism3d_ispresent ) THEN
         lrism3d = .TRUE.
         CALL qexsd_copy_rism3d ( output_obj%rism3d, pseudo_dir, nsolV, solVs, molfile, ecutsolv )
         ecutsolv = ecutsolv * e2
      ELSE
         lrism3d  = .FALSE.
         nsolV    = 0
         ecutsolv = 0.0_DP
      END IF
      !
      !! Laue-RISM
      IF ( output_obj%rismlaue_ispresent ) THEN
         CALL rism3d_set_laue()
         CALL qexsd_copy_rismlaue ( output_obj%rismlaue, both_hands, laue_nfit, ireference, qsol, &
                                    starting_r, expand_r, buffer_r, buffer_ru, buffer_rv, &
                                    starting_l, expand_l, buffer_l, buffer_lu, buffer_lv )
      ELSE
         both_hands = .FALSE.
         laue_nfit  = 0
         expand_r   = -1.0_DP
         expand_l   = -1.0_DP
         starting_r = 0.0_DP
         starting_l = 0.0_DP
      END IF
      !
      ! ... xml data no longer needed, can be discarded
      !
      CALL qes_reset  ( output_obj )
      CALL qes_reset  ( geninfo_obj )
      CALL qes_reset  ( parinfo_obj )
      IF ( TRIM(input_obj%tagname) == "input") CALL qes_reset ( input_obj) 
      !
      ! END OF READING VARIABLES FROM XML DATA FILE
      !
      ALLOCATE( force ( 3, nat ) )
      ALLOCATE( extfor( 3, nat ) )
      IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
      IF ( gate ) ALLOCATE( forcegate( 3, nat ) )
      !
    END SUBROUTINE read_xml_file
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_collected_wfc ( dirname, ik, arr, label_, ierr_ )
      !------------------------------------------------------------------------
      !
      ! ... reads from directory "dirname" (new file format) for k-point "ik"
      ! ... wavefunctions from collected format into distributed array "arr"
      !
      USE control_flags,        ONLY : gamma_only
      USE lsda_mod,             ONLY : nspin, isk
      USE klist,                ONLY : nkstot, nks, ngk, igk_k
      USE wvfct,                ONLY : npwx, nbnd
      USE gvect,                ONLY : ig_l2g
      USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm
      USE mp_pools,             ONLY : me_pool, root_pool, intra_pool_comm
      USE mp,                   ONLY : mp_sum, mp_max
      USE io_base,              ONLY : read_wfc
      USE xc_lib,               ONLY : exx_is_active
      USE exx,                  ONLY : nbndproj
      USE io_global,            ONLY : stdout
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      INTEGER, INTENT(IN) :: ik
      COMPLEX(dp), INTENT(OUT) :: arr(:,:)
      CHARACTER(LEN=3), OPTIONAL, INTENT(IN) :: label_
      INTEGER, OPTIONAL, INTENT(OUT)  :: ierr_
      !
      CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
      CHARACTER(LEN=320)   :: filename, msg 
      CHARACTER(LEN=3)     :: label 
      LOGICAL              :: read_ace
      INTEGER              :: i, ik_g, ig
      INTEGER              :: npol_, nbnd_
      INTEGER              :: ike, iks, ngk_g, npw_g, ispin
      INTEGER, EXTERNAL    :: global_kpoint_index
      INTEGER, ALLOCATABLE :: mill_k(:,:)
      INTEGER, ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:)
      LOGICAL              :: ionode_k
      REAL(DP)             :: scalef, xk_(3), b1(3), b2(3), b3(3)
      !
      ! ... decide whether to read wfc or ace
      !
      if(present(label_)) then 
         label = label_
         if(label.eq."ace") then 
            if(.not.exx_is_active()) CALL errore ('pw_restart-read_collected_wfc',&
                 "ace but not exx_is_active", 1 ) 
            read_ace = .true.
         else if(label.eq."wfc") then
            read_ace = .false.
         else
            CALL errore ('pw_restart - read_collected_wfc', "wrong label", 1 )
         end if
      else
         label = "wfc"
         read_ace = .false.
      end if
      !
      ! ... the root processor of each pool reads
      !
      ionode_k = (me_pool == root_pool)
      !
      iks = global_kpoint_index (nkstot, 1)
      ike = iks + nks - 1
      !
      ! ik_g: index of k-point ik in the global list
      !
      ik_g = ik + iks - 1
      !
      ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
      !
      ALLOCATE ( igk_l2g_kdip( npwx ) )
      !
      ! ... The igk_l2g array yields the correspondence between the
      ! ... local k+G index and the global G index - requires arrays
      ! ... igk_k (k+G indices) and ig_l2g (local to global G index map)
      !
      ALLOCATE ( igk_l2g( npwx ) )
      igk_l2g = 0
      DO ig = 1, ngk(ik)
         igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
      END DO
      !
      ! ... npw_g: the maximum G vector index among all processors
      ! ... ngk_g: global number of k+G vectors for all k points
      !
      npw_g = MAXVAL( igk_l2g(1:ngk(ik)) )
      CALL mp_max( npw_g, intra_pool_comm )
      ngk_g = ngk(ik)
      CALL mp_sum( ngk_g, intra_bgrp_comm)
      !
      ! ... now compute the igk_l2g_kdip local-to-global map
      !
      igk_l2g_kdip = 0
      CALL gk_l2gmap_kdip( npw_g, ngk_g, ngk(ik), igk_l2g, &
           igk_l2g_kdip )
      DEALLOCATE ( igk_l2g )
      !
      IF ( nspin == 2 ) THEN
         !
         ! ... LSDA: spin mapped to k-points, isk(ik) tracks up and down spin
         !
         ik_g = MOD ( ik_g-1, nkstot/2 ) + 1 
         ispin = isk(ik)
         filename = TRIM(dirname) // label // updw(ispin) // &
              & TRIM(int_to_char(ik_g))
         !
      ELSE
         !
         filename = TRIM(dirname) // label // TRIM(int_to_char(ik_g))
         !
      ENDIF
      !
      ! ... Miller indices are read from file (but not used)
      !
      ALLOCATE( mill_k ( 3,npwx ) )
      !
      arr = (0.0_DP, 0.0_DP)
      !
      CALL read_wfc( iunpun, filename, root_bgrp, intra_bgrp_comm, &
           ik_g, xk_, ispin, npol_, arr, npw_g, gamma_only, nbnd_, &
           igk_l2g_kdip(:), ngk(ik), b1, b2, b3, mill_k, scalef, ierr_ )
      !
      DEALLOCATE ( mill_k )
      DEALLOCATE ( igk_l2g_kdip )
      !
      IF ( PRESENT (ierr_) ) THEN
         IF ( ierr_ /= 0 ) RETURN
      END IF
      !
      ! ... here one should check for consistency between what is read
      ! ... and what is expected
      !
      IF(read_ace) THEN
        !
        WRITE(stdout, '(5X,A,I8,A)') 'ACE potential read for ', nbnd_, ' bands'
        nbndproj = nbnd_
        !
      ELSE IF ( nbnd_ < nbnd .and..not. read_ace) THEN
        !
        WRITE (msg,'("The number of bands for this run is",I6,", but only",&
             & I6," bands were read from file")')  nbnd, nbnd_  
        CALL errore ('pw_restart - read_collected_wfc', msg, 1 )
        !
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_collected_wfc
    !
    !------------------------------------------------------------------------
  END MODULE pw_restart_new
