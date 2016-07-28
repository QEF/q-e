!
! Copyright (C) 2005-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE pw_restart
  !----------------------------------------------------------------------------
  !
  ! ... this module contains methods to read and write data produced by PWscf
  !
  ! ... originally written by Carlo Sbraccia  (2005)
  !
  !
  USE iotk_module
  !
#ifdef __XSD
  USE qes_module
  USE qexsd_module, ONLY: qexsd_init_schema, qexsd_openschema, qexsd_closeschema, &
                          qexsd_init_convergence_info, qexsd_init_algorithmic_info, & 
                          qexsd_init_atomic_species, qexsd_init_atomic_structure, &
                          qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
                          qexsd_init_magnetization,qexsd_init_band_structure,  &
                          qexsd_init_total_energy,qexsd_init_forces,qexsd_init_stress, &
                          qexsd_init_outputElectricField
            
#endif
  USE qexml_module, ONLY: qexml_init,qexml_openfile, qexml_closefile, &
                          qexml_write_header, qexml_write_control ,   &
                          qexml_write_cell, qexml_write_moving_cell,  &
                          qexml_write_ions, qexml_write_symmetry,     &
                          qexml_write_efield, qexml_write_planewaves, &
                          qexml_write_spin, qexml_write_magnetization, &
                          qexml_write_xc, qexml_write_exx, qexml_write_occ, &
                          qexml_write_bz,qexml_write_para, qexml_write_bands_info, &
                          qexml_write_bands_pw, qexml_write_esm, qexml_wfc_filename, &
                          default_fmt_version => qexml_default_version, &
                          qexml_kpoint_dirname, &
                          qexml_read_header, qexml_read_cell, qexml_read_moving_cell, &
                          qexml_read_planewaves, qexml_read_ions, qexml_read_spin, &
                          qexml_read_magnetization, qexml_read_xc, qexml_read_occ, qexml_read_bz, &
                          qexml_read_bands_info, qexml_read_bands_pw, qexml_read_symmetry, &
                          qexml_read_efield, qexml_read_para, qexml_read_exx, qexml_read_esm
  !
  USE xml_io_base, ONLY : rho_binary,read_wfc, write_wfc, create_directory
  !
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2, PI
  !
#ifdef __XSD
  USE io_files,  ONLY : iunpun_xsd, xmlpun_schema
#endif
USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun, delete_if_present, &
                        qexml_version, qexml_version_init, pseudo_dir
  !
  USE io_global, ONLY : ionode, ionode_id
  USE mp_images, ONLY : intra_image_comm
  USE mp_pools,  ONLY : my_pool_id
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_bcast, mp_sum, mp_max
  USE parser,    ONLY : version_compare
  !
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), external :: trimcheck
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: pw_writefile, pw_readfile
#ifdef __XSD
  ! 
  PUBLIC :: pw_write_schema, pw_readschema_file, init_vars_from_schema, read_collected_to_evc
#endif
  !
  INTEGER, PRIVATE :: iunout
  !
  LOGICAL :: lcell_read   = .FALSE., &
             lpw_read     = .FALSE., &
             lions_read   = .FALSE., &
             lspin_read   = .FALSE., &
             lstarting_mag_read   = .FALSE., &
             lxc_read     = .FALSE., &
             locc_read    = .FALSE., &
             lbz_read     = .FALSE., &
             lbs_read     = .FALSE., &
             lefield_read = .FALSE., &
             lwfc_read    = .FALSE., &
             lsymm_read   = .FALSE.
  !
  !
  CONTAINS
  !
#ifdef __XSD
    !------------------------------------------------------------------------
    SUBROUTINE pw_write_schema( what )
      !------------------------------------------------------------------------
      !
      USE control_flags,        ONLY : istep, twfcollect, conv_ions, &
                                       lscf, lkpoint_dir, gamma_only, &
                                       tqr, tq_smoothing, tbeta_smoothing, &
                                       noinv, do_makov_payne, smallmem, &
                                       llondon, lxdm, ts_vdw, scf_error, n_scf_steps
      USE realus,               ONLY : real_space
      USE uspp,                 ONLY : okvan
      USE paw_variables,        ONLY : okpaw
      USE uspp_param,           ONLY : upf
      USE global_version,       ONLY : version_number
      USE cell_base,            ONLY : at, bg, alat, tpiba, tpiba2, &
                                       ibrav, celldm
      USE gvect,                ONLY : ig_l2g
      USE ions_base,            ONLY : nsp, ityp, atm, nat, tau, if_pos
      USE noncollin_module,     ONLY : noncolin, npol
      USE io_files,             ONLY : nwordwfc, iunwfc, psfile
      USE buffers,              ONLY : get_buffer
      USE wavefunctions_module, ONLY : evc
      USE klist,                ONLY : nks, nkstot, xk, ngk, wk, qnorm, &
                                       lgauss, ngauss, degauss, nelec, &
                                       two_fermi_energies, nelup, neldw
      USE start_k,              ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       nks_start, xk_start, wk_start
      USE ktetra,               ONLY : ntetra, tetra, ltetra
      USE gvect,                ONLY : ngm, ngm_g, g, mill
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE gvecs,                ONLY : ngms_g, dual
      USE fft_base,             ONLY : dffts
      USE wvfct,                ONLY : npwx, et, wg, nbnd
      USE ener,                 ONLY : ef, ef_up, ef_dw, vtxc, etxc, ewld, etot, &
                                       ehart, eband, demet 
      USE gvecw,                ONLY : ecutwfc
      USE fixed_occ,            ONLY : tfixed_occ, f_inp
      USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, U_projection, &
                                       Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, &
                                       Hubbard_alpha, Hubbard_J0, Hubbard_beta,&
                                       is_hubbard
      USE spin_orb,             ONLY : lspinorb, domag
      USE symm_base,            ONLY : nrot, nsym, invsym, s, ft, irt, &
                                       t_rev, sname, time_reversal, no_t_rev
      USE lsda_mod,             ONLY : nspin, isk, lsda, starting_magnetization, magtot, absmag
      USE noncollin_module,     ONLY : angle1, angle2, i_cons, mcons, bfield, magtot_nc, &
                                       lambda
      USE ions_base,            ONLY : amass
      USE funct,                ONLY : get_dft_name, get_inlc, get_nonlocc_name, dft_is_nonlocc
      USE kernel_table,         ONLY : vdw_table_name
      USE scf,                  ONLY : rho
      USE force_mod,            ONLY : lforce, sumfor, force, sigma, lstres
      USE extfield,             ONLY : tefield, dipfield, edir, etotefield, &
                                       emaxpos, eopreg, eamp
      USE io_rho_xml,           ONLY : write_rho
      USE mp_world,             ONLY : nproc
      USE mp_images,            ONLY : nproc_image
      USE mp_pools,             ONLY : kunit, nproc_pool, me_pool, root_pool, &
                                       intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : nproc_bgrp, me_bgrp, root_bgrp, &
                                       intra_bgrp_comm, nbgrp, ntask_groups
      USE mp_diag,              ONLY : nproc_ortho
      USE funct,                ONLY : get_exx_fraction, dft_is_hybrid, &
                                       get_gau_parameter, &
                                       get_screening_parameter, exx_is_active
      USE exx,                  ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut, ecutfock
      USE cellmd,               ONLY : lmovecell, cell_factor 
      USE martyna_tuckerman,    ONLY : do_comp_mt
      USE esm,                  ONLY : do_comp_esm, esm_nfit, esm_efield, esm_w, &
                                       esm_a, esm_bc
      USE london_module,        ONLY : scal6, lon_rcut, in_c6
      USE xdm_module,           ONLY : xdm_a1=>a1i, xdm_a2=>a2i
      USE tsvdw_module,         ONLY : vdw_isolated, vdw_econv_thr
      USE input_parameters,     ONLY : space_group, verbosity, calculation, ion_dynamics, starting_ns_eigenvalue, &
                                       vdw_corr, london, input_parameters_occupations => occupations
      USE bp,                   ONLY : lelfield, lberry, bp_mod_el_pol => el_pol, bp_mod_ion_pol => ion_pol
      !
      USE rap_point_group,      ONLY : elem, nelem, name_class
      USE rap_point_group_so,   ONLY : elem_so, nelem_so, name_class_so
      USE bfgs_module,          ONLY : bfgs_get_n_iter
      USE qexsd_module,         ONLY : qexsd_dipol_obj, qexsd_bp_obj
     
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      !
      CHARACTER(15)         :: subname="pw_write_schema"
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname, filename
      CHARACTER(LEN=80)     :: vdw_corr_
      INTEGER               :: i, ig, ik, ngg, ierr, ipol, num_k_points
      INTEGER               :: npool, nkbl, nkl, nkr, npwx_g
      INTEGER               :: ike, iks, npw_g, ispin, inlc
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      INTEGER,  ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:), mill_g(:,:)
      LOGICAL               :: lwfc, lrho, lxsd, occupations_are_fixed
      CHARACTER(iotk_attlenx)  :: attr
      INTEGER                  :: iclass, isym, ielem
      CHARACTER(LEN=15)        :: symop_2_class(48)
      LOGICAL                  :: opt_conv_ispresent
      INTEGER                  :: n_opt_steps, h_band
      REAL(DP)                 :: h_energy
      !
      TYPE(output_type) :: output
      
      !
      ! PW dimensions need to be properly computed 
      ! reducing across MPI tasks
      !
      ALLOCATE( ngk_g( nkstot ) )
      !
      ngk_g(1:nks) = ngk(:)
      CALL mp_sum( ngk_g(1:nks), intra_bgrp_comm )
      CALL ipoolrecover( ngk_g, 1, nkstot, nks )
      !CALL mp_bcast( ngk_g, root_bgrp, intra_bgrp_comm )
      !CALL mp_bcast( ngk_g, root_bgrp, inter_bgrp_comm )

      !
      ! ... compute the maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      !DEALLOCATE( ngk_g )
      ! do not deallocate ngk_g here, please, I need it for band_structure_init
      ! P. Delugas 
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      ! 
      IF (tefield .AND. dipfield ) THEN 
          CALL init_dipole_info(qexsd_dipol_obj, rho%of_r)   
          qexsd_dipol_obj%tagname = "dipoleInfo"
      END IF
      ! 
      !
      ! XML descriptor
      ! 
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      CALL qexsd_init_schema( iunpun_xsd )
      !
      !
      IF ( ionode ) THEN  
         !
         ! ... here we init the variables and finally write them to file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL qexsd_openschema(TRIM( dirname ) // '/' // TRIM( xmlpun_schema ))
         output%tagname="output"
         !
!-------------------------------------------------------------------------------
! ... CONVERGENCE_INFO
!-------------------------------------------------------------------------------
         SELECT CASE (TRIM( calculation )) 
            CASE ( "relax","vc-relax" ,"md")
                opt_conv_ispresent = .TRUE.
                IF (TRIM( ion_dynamics) == 'bfgs' ) THEN 
                    n_opt_steps = bfgs_get_n_iter('bfgs_iter ') 
                ELSE 
                    n_opt_steps = istep 
                END IF 
            CASE default
                opt_conv_ispresent = .FALSE.
                n_opt_steps        = 0 
         END SELECT
         ! 
         call qexsd_init_convergence_info(output%convergence_info, &
                                          n_scf_steps=n_scf_steps, scf_error=scf_error, &
                                          opt_conv_ispresent=lforce, &
                                          n_opt_steps=n_opt_steps, grad_norm=sumfor )
         !
!-------------------------------------------------------------------------------
! ... ALGORITHMIC_INFO
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_algorithmic_info(output%algorithmic_info, &
                                          real_space_q=real_space, uspp=okvan, paw=okpaw)
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_SPECIES
!-------------------------------------------------------------------------------
         !
         ! while amass's are always present, starting_mag should not be passed
         ! for nspin==1 or contrained magnetization calculations
         !
         IF (noncolin) THEN
             CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm, psfile, &
                                        amass, STARTING_MAGNETIZATION = starting_magnetization, &
                                        ANGLE1=angle1, ANGLE2=angle2)
         ELSE IF (nspin==2) THEN 
             CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm, psfile, &
                                           amass, STARTING_MAGNETIZATION=starting_magnetization)
         ELSE 
             CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm,psfile, &
                                          amass)
         END IF
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_STRUCTURE
!-------------------------------------------------------------------------------
         !         
         CALL qexsd_init_atomic_structure(output%atomic_structure, nsp, atm, ityp, &
                       nat, tau, 'Bohr', alat, alat*at(:,1), alat*at(:,2), alat*at(:,3), ibrav)
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
                        IF ( elem(ielem,iclass) == isym) THEN
                           symop_2_class(isym) = name_class(iclass)
                           EXIT classes_loop
                        END IF
                     END DO elements_loop
                  END DO classes_loop
               END DO symmetries_loop
            END IF
         END IF
         CALL qexsd_init_symmetries(output%symmetries, nsym, nrot, space_group, &
                                    s, ft, sname, t_rev, nat, irt,symop_2_class(1:nrot), verbosity, &
                                    noncolin)
         !
!-------------------------------------------------------------------------------
! ... BASIS SET
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_basis_set(output%basis_set, gamma_only, ecutwfc/e2, ecutwfc*dual/e2, &
                                   dfftp%nr1, dfftp%nr2, dfftp%nr3, dffts%nr1, dffts%nr2, dffts%nr3, &
                                   .FALSE., dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, ngms_g, npwx_g, &
                                   bg(:,1), bg(:,2), bg(:,3) )
         !
!-------------------------------------------------------------------------------
! ... DFT
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         inlc = get_inlc()
         !
         IF ( lda_plus_u .AND. noncolin) CALL errore(subname,"LDA+U and non-collinear case not implemented in qexsd",10)
         !
         vdw_corr_ = vdw_corr
         IF ( london ) vdw_corr_ = 'grimme-d2'
         CALL qexsd_init_dft(output%dft, dft_name, .TRUE., &
                             dft_is_hybrid(), nq1, nq2, nq3, ecutfock, &
                             get_exx_fraction(), get_screening_parameter(), exxdiv_treatment, &
                             x_gamma_extrapolation, ecutvcut, lda_plus_u, lda_plus_u_kind, 2*Hubbard_lmax+1, &
                             nspin, nsp, 2*Hubbard_lmax+1, nat, atm, ityp, Hubbard_U, Hubbard_J0, Hubbard_alpha, &
                             Hubbard_beta, Hubbard_J, starting_ns_eigenvalue, rho%ns, rho%ns_nc, U_projection, &
                             dft_is_nonlocc(), TRIM(vdw_corr_), TRIM ( get_nonlocc_name()), scal6, in_c6, lon_rcut, xdm_a1, xdm_a2,&
                             vdw_econv_thr, vdw_isolated, is_hubbard, upf(1:nsp)%psd)
         !
!-------------------------------------------------------------------------------
! ... MAGNETIZATION
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_magnetization(output%magnetization, lsda, noncolin, lspinorb, &
                                       magtot, magtot_nc, absmag, domag )
         !

!--------------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------------
         !
         IF (TRIM(input_parameters_occupations) == 'fixed') THEN 
            occupations_are_fixed = .TRUE. 
            IF ( nspin == 1 ) THEN 
               h_band = NINT ( nelec/2.d0) 
            ELSE 
               h_band = NINT ( nelec ) 
            END IF  
            h_energy =MAXVAL (et(h_band, 1:nkstot))
         ELSE 
            occupations_are_fixed = .FALSE. 
            h_energy  = ef 
         END IF 
         CALL  qexsd_init_band_structure(output%band_structure,lsda,noncolin,lspinorb, &
                                         nbnd,nelec, natomwfc, occupations_are_fixed, & 
                                         h_energy,two_fermi_energies, [ef_up,ef_dw], et,wg,nkstot,xk,ngk_g,wk)
         !
!-------------------------------------------------------------------------------------------
! ... TOTAL ENERGY
!-------------------------------------------------------------------------------------------
         !
         IF (tefield) THEN
            CALL  qexsd_init_total_energy(output%total_energy,etot,eband,ehart,vtxc,etxc, &
                                       ewld,degauss,demet, etotefield)
         ELSE 
            CALL  qexsd_init_total_energy(output%total_energy,etot,eband,ehart,vtxc,etxc, &
                                       ewld,degauss,demet)
         END IF
         !
!---------------------------------------------------------------------------------------------
! ... FORCES
!----------------------------------------------------------------------------------------------
         !
         IF ( lforce ) THEN 
            output%forces_ispresent = .TRUE.
            CALL qexsd_init_forces(output%forces,nat,force,lforce)
         ELSE 
            output%forces_ispresent = .FALSE.
            output%forces%lwrite = .FALSE.  
         END IF 
         !
!------------------------------------------------------------------------------------------------
! ... STRESS 
!------------------------------------------------------------------------------------------------
         IF ( lstres) THEN
            output%stress_ispresent=.TRUE.
            CALL qexsd_init_stress(output%stress, sigma, lstres ) 
         ELSE 
            output%stress_ispresent=.FALSE.
            output%stress%lwrite=.FALSE.
         END IF
!-------------------------------------------------------------------------------------------------
! ... ELECTRIC FIELD
!-------------------------------------------------------------------------------------------------
         IF ( lelfield ) THEN
            output%electric_field_ispresent = .TRUE. 
            CALL qexsd_init_outputElectricField(output%electric_field, lelfield, tefield, dipfield, &
                                                lberry, el_pol = bp_mod_el_pol, ion_pol = bp_mod_ion_pol) 
         ELSE IF ( lberry ) THEN 
            output%electric_field_ispresent = .TRUE.
            CALL qexsd_init_outputElectricField(output%electric_field, lelfield, tefield, dipfield, & 
                                                lberry, bp_obj=qexsd_bp_obj) 
         ELSE IF ( tefield .AND. dipfield  ) THEN 
            output%electric_field_ispresent = .TRUE.
            CALL  qexsd_init_outputElectricField(output%electric_field, lelfield, tefield, dipfield, &
                                                  lberry, dipole_obj = qexsd_dipol_obj )                     
         ELSE 
            output%electric_field_ispresent = .FALSE.
         ENDIF
!------------------------------------------------------------------------------------------------
! ... ACTUAL WRITING
!-------------------------------------------------------------------------------
         !
         CALL qes_write_output(iunpun_xsd,output)
         CALL qes_reset_output(output)
         !
!-------------------------------------------------------------------------------
! ... CLOSING
!-------------------------------------------------------------------------------
         !
         CALL qexsd_closeschema()
         !
      END IF
      DEALLOCATE (ngk_g)
      !
      RETURN
       !
    END SUBROUTINE pw_write_schema
    !
#endif
    !------------------------------------------------------------------------
    SUBROUTINE pw_writefile( what )
      !------------------------------------------------------------------------
      !
      USE control_flags,        ONLY : twfcollect, conv_ions, &
                                       lscf, lkpoint_dir, gamma_only, &
                                       tqr, tq_smoothing, tbeta_smoothing, &
                                       noinv, do_makov_payne, smallmem, &
                                       llondon, lxdm, ts_vdw 
      USE realus,               ONLY : real_space
      USE global_version,       ONLY : version_number
      USE cell_base,            ONLY : at, bg, alat, tpiba, tpiba2, &
                                       ibrav, celldm
      USE gvect,                ONLY : ig_l2g
      USE ions_base,            ONLY : nsp, ityp, atm, nat, tau, if_pos
      USE noncollin_module,     ONLY : noncolin, npol
      USE io_files,             ONLY : nwordwfc, iunwfc, psfile
      USE buffers,              ONLY : get_buffer
      USE wavefunctions_module, ONLY : evc
      USE klist,                ONLY : nks, nkstot, xk, ngk, igk_k, wk, qnorm, &
                                       lgauss, ngauss, degauss, nelec, &
                                       two_fermi_energies, nelup, neldw
      USE start_k,              ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       nks_start, xk_start, wk_start
      USE ktetra,               ONLY : ntetra, tetra, ltetra
      USE gvect,                ONLY : ngm, ngm_g, g, mill
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE gvecs,                ONLY : ngms_g, dual
      USE fft_base,             ONLY : dffts
      USE wvfct,                ONLY : npw, npwx, et, wg, nbnd
      USE gvecw,                ONLY : ecutwfc
      USE ener,                 ONLY : ef, ef_up, ef_dw
      USE fixed_occ,            ONLY : tfixed_occ, f_inp
      USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, U_projection, &
                                       Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, &
                                       Hubbard_alpha, Hubbard_J0, Hubbard_beta
      USE spin_orb,             ONLY : lspinorb, domag, lforcet
      USE symm_base,            ONLY : nrot, nsym, invsym, s, ft, irt, &
                                       t_rev, sname, time_reversal, no_t_rev
      USE lsda_mod,             ONLY : nspin, isk, lsda, starting_magnetization
      USE noncollin_module,     ONLY : angle1, angle2, i_cons, mcons, bfield, &
                                       lambda
      USE ions_base,            ONLY : amass
      USE funct,                ONLY : get_dft_name, get_inlc
      USE kernel_table,         ONLY : vdw_table_name
      USE scf,                  ONLY : rho
      USE extfield,             ONLY : tefield, dipfield, edir, &
                                       emaxpos, eopreg, eamp
      USE io_rho_xml,           ONLY : write_rho
      USE mp_world,             ONLY : nproc
      USE mp_images,            ONLY : nproc_image
      USE mp_pools,             ONLY : kunit, nproc_pool, me_pool, root_pool, &
                                       intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : nproc_bgrp, me_bgrp, root_bgrp, &
                                       intra_bgrp_comm, nbgrp, ntask_groups
      USE mp_diag,              ONLY : nproc_ortho
      USE funct,                ONLY : get_exx_fraction, dft_is_hybrid, &
                                       get_gau_parameter, &
                                       get_screening_parameter, exx_is_active
      USE exx,                  ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut, ecutfock
      USE cellmd,               ONLY : lmovecell, cell_factor 
      USE martyna_tuckerman,    ONLY : do_comp_mt
      USE esm,                  ONLY : do_comp_esm, esm_nfit, esm_efield, esm_w, &
                                       esm_a, esm_bc
      USE acfdt_ener,           ONLY : acfdt_in_pw 
      USE london_module,        ONLY : scal6, lon_rcut
      USE tsvdw_module,         ONLY : vdw_isolated

      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      !
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname, filename
      INTEGER               :: i, ig, ik, ngg, ierr, ipol, num_k_points
      INTEGER               :: npool, nkbl, nkl, nkr, npwx_g
      INTEGER               :: ike, iks, npw_g, ispin, inlc
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      INTEGER,  ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:), mill_g(:,:)
      LOGICAL               :: lwfc, lrho, lxsd
      CHARACTER(iotk_attlenx)  :: attr
      !
      !
      SELECT CASE( what )
      CASE( "all" )
         !
         ! ... do not overwrite the scf charge density with a non-scf one
         ! ... (except in the 'force theorem' calculation of MAE where the
         ! ...  charge density differs from the one read from disk)
         !
         lrho  = lscf .OR. lforcet
         lwfc  = twfcollect
         !
      CASE( "config" )
         ! 
         ! ... write just the xml data file, not the charge density and the wavefunctions
         !
         lwfc  = .FALSE.
         lrho  = .FALSE.
         !
      CASE DEFAULT
         !
         CALL errore( 'pw_writefile', 'unexpected case: '//TRIM(what), 1 )
         ! 
      END SELECT
      !
      IF ( ionode ) THEN
         !
         ! ... look for an empty unit (only ionode needs it)
         !
         CALL iotk_free_unit( iunout, ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'pw_writefile ', &
                   'no free units to write wavefunctions', ierr )
      !
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      ! ... create the main restart directory
      !
      CALL create_directory( dirname )
      !
      ! ... create the k-points subdirectories
      !
      IF ( nspin == 2 ) THEN
         num_k_points = nkstot / 2
      ELSE
         num_k_points = nkstot
      END IF
      !
      IF (lkpoint_dir) THEN
         !
         DO i = 1, num_k_points
            !
            CALL create_directory( qexml_kpoint_dirname( dirname, i ) )
            !
         END DO
         !
      END IF
      !
      IF ( nkstot > 0 ) THEN
         !
         ! ... find out the number of pools
         !
         npool = nproc_image / nproc_pool
         !
         ! ... find out number of k points blocks
         !
         nkbl = nkstot / kunit
         !
         ! ... k points per pool
         !
         nkl = kunit * ( nkbl / npool )
         !
         ! ... find out the reminder
         !
         nkr = ( nkstot - nkl * npool ) / kunit
         !
         ! ... Assign the reminder to the first nkr pools
         !
         IF ( my_pool_id < nkr ) nkl = nkl + kunit
         !
         ! ... find out the index of the first k point in this pool
         !
         iks = nkl*my_pool_id + 1
         !
         IF ( my_pool_id >= nkr ) iks = iks + nkr*kunit
         !
         ! ... find out the index of the last k point in this pool
         !
         ike = iks + nkl - 1
         !
      END IF
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      !
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      !
      ! ... collect all G-vectors across processors within the pools
      !
      ALLOCATE( mill_g( 3, ngm_g ) )
      !
      mill_g = 0
      !
      DO ig = 1, ngm
         !
         mill_g(1,ig_l2g(ig)) = mill(1,ig)
         mill_g(2,ig_l2g(ig)) = mill(2,ig)
         mill_g(3,ig_l2g(ig)) = mill(3,ig)
         !
      END DO
      !
      CALL mp_sum( mill_g, intra_bgrp_comm )
      !
      ! ... build the igk_l2g array, yielding the correspondence between
      ! ... the local k+G index and the global G index - see also ig_l2g
      ! ... igk_l2g is build from arrays igk, previously stored in hinit0
      ! ... Beware: for variable-cell case, one has to use starting G and 
      ! ... k+G vectors
      !
      ALLOCATE ( igk_l2g( npwx, nks ) )
      igk_l2g = 0
      !
      DO ik = 1, nks
         npw = ngk (ik)
         CALL gk_l2gmap( ngm, ig_l2g(1), npw, igk_k(1,ik), igk_l2g(1,ik) )
      END DO
      !
      ! ... compute the global number of G+k vectors for each k point
      !
      ALLOCATE( ngk_g( nkstot ) )
      !
      ngk_g = 0
      ngk_g(iks:ike) = ngk(1:nks)
      !
      CALL mp_sum( ngk_g, inter_pool_comm)
      CALL mp_sum( ngk_g, intra_pool_comm)
      !
      ngk_g = ngk_g / nbgrp
      !
      ! ... compute the maximum G vector index among all G+k an processors
      !
      npw_g = MAXVAL( igk_l2g(:,:) )
      !
      CALL mp_max( npw_g, inter_pool_comm )
      CALL mp_max( npw_g, intra_pool_comm )
      !
      ! ... compute the maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! ... define a further l2g map to write gkvectors and wfc coherently
      !
      ALLOCATE ( igk_l2g_kdip( npwx_g, nks ) )
      !
      igk_l2g_kdip = 0
      !
      DO ik = iks, ike
         !
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik-iks+1), &
                              igk_l2g(1,ik-iks+1), igk_l2g_kdip(1,ik-iks+1) )
      END DO
      !
      IF ( ionode ) THEN
         !
         ! ... open XML descriptor
         !
         CALL qexml_init( iunpun )
         CALL qexml_openfile( TRIM( dirname ) // '/' // TRIM( xmlpun ), &
                              'write', BINARY = .FALSE., IERR = ierr  )
         !
         IF (.NOT.(lkpoint_dir)) &
            CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // &
                    & TRIM( xmlpun )//'.eig', BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'pw_writefile ', &
                   'cannot open restart file for writing', ierr )
      !
      IF ( ionode ) THEN  
         !
         ! ... here we start writing the punch-file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_header( "PWSCF", TRIM(version_number) )
         !
!-------------------------------------------------------------------------------
! ... CONTROL 
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_control( PP_CHECK_FLAG=conv_ions, LKPOINT_DIR=lkpoint_dir, &
                            Q_REAL_SPACE=tqr, TQ_SMOOTHING=tq_smoothing, &
                            BETA_REAL_SPACE=real_space, TBETA_SMOOTHING=tbeta_smoothing )
         !
!-------------------------------------------------------------------------------
! ... CELL
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_cell( ibrav, celldm, alat, &
                          at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3), &
                          "Bohr","Bohr","2 pi / a", &
                          do_makov_payne, do_comp_mt, do_comp_esm )
         !
         IF (lmovecell) CALL qexml_write_moving_cell(lmovecell, cell_factor)
         !
!-------------------------------------------------------------------------------
! ... IONS
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_ions( nsp, nat, atm, ityp, psfile, &
                       pseudo_dir, amass, 'a.m.u.', tau, 'Bohr', if_pos, dirname, alat )
         !
!-------------------------------------------------------------------------------
! ... SYMMETRIES
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_symmetry( ibrav, nrot, nsym, invsym, noinv, &
                              time_reversal, no_t_rev, ft, s, sname, "Crystal", irt,  &
                              nat, t_rev )
         !
!-------------------------------------------------------------------------------
! ... ELECTRIC FIELD
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp) 
         !
!
!-------------------------------------------------------------------------------
! ... PLANE_WAVES
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_planewaves( ecutwfc/e2, ecutwfc*dual/e2, npwx_g, gamma_only, &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, &
                                dffts%nr1, dffts%nr2, dffts%nr3, ngms_g, dfftp%nr1, &
                                dfftp%nr2, dfftp%nr3, mill_g, lwfc,'Hartree' )
         !
!-------------------------------------------------------------------------------
! ... SPIN
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_spin( lsda, noncolin, npol, lspinorb, domag )
         !
         CALL qexml_write_magnetization(starting_magnetization, &
                                  angle1*180.0_DP/PI , angle2*180.0_DP/PI, nsp, &
                                  two_fermi_energies, i_cons, mcons, bfield, &
                                  ef_up/e2, ef_dw/e2, nelup, neldw, lambda,'Hartree')
         !
!-------------------------------------------------------------------------------
! ... EXCHANGE_CORRELATION
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         inlc = get_inlc()
         !
         CALL qexml_write_xc( DFT = dft_name, NSP = nsp, LDA_PLUS_U = lda_plus_u,      &
                        LDA_PLUS_U_KIND = lda_plus_u_kind, U_PROJECTION = U_projection, &
                        HUBBARD_LMAX = Hubbard_lmax, HUBBARD_L = Hubbard_l, &
                        HUBBARD_U = Hubbard_U, HUBBARD_J = Hubbard_J, &
                        HUBBARD_J0 = Hubbard_J0, HUBBARD_BETA = Hubbard_beta, &
                        HUBBARD_ALPHA = Hubbard_alpha, &
                        INLC = inlc, VDW_TABLE_NAME = vdw_table_name, &
                        PSEUDO_DIR = pseudo_dir, DIRNAME = dirname, &
                        ACFDT_IN_PW = acfdt_in_pw, &
                        LLONDON = llondon, LONDON_S6 = scal6,         &
                        LONDON_RCUT = lon_rcut, LXDM = lxdm,          &
                        TS_VDW = ts_vdw, VDW_ISOLATED = vdw_isolated )


         IF ( dft_is_hybrid() ) CALL qexml_write_exx &
                       ( x_gamma_extrapolation, nq1, nq2, nq3, &
                         exxdiv_treatment, yukawa, ecutvcut, &
                         get_exx_fraction(), get_gau_parameter(), &
                         get_screening_parameter(), exx_is_active(), ecutfock )
         !
!-------------------------------------------------------------------------------
! ... ESM
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_esm( esm_nfit, esm_efield, esm_w, esm_a, esm_bc )
         !
!-------------------------------------------------------------------------------
! ... OCCUPATIONS
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_occ( LGAUSS = lgauss, NGAUSS = ngauss, &
                         DEGAUSS = degauss/e2,DEGAUSS_UNITS='Hartree', LTETRA = ltetra, NTETRA = ntetra, &
                         TETRA = tetra, TFIXED_OCC = tfixed_occ, LSDA = lsda, &
                         NSTATES_UP = nbnd, NSTATES_DW = nbnd, INPUT_OCC = f_inp )
         !
!-------------------------------------------------------------------------------
! ... BRILLOUIN_ZONE
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3, &
                        '2 pi / a',qnorm, nks_start, xk_start, wk_start )
         !
!-------------------------------------------------------------------------------
! ... PARALLELISM
!-------------------------------------------------------------------------------
         !
         !
         CALL qexml_write_para( kunit, nproc, nproc_pool, nproc_image, &
                                ntask_groups, nproc_bgrp, nproc_ortho )
         !
!-------------------------------------------------------------------------------
! ... CHARGE DENSITY
!-------------------------------------------------------------------------------
         !
         !
         filename = "./charge-density.dat"
         IF ( .NOT. rho_binary ) filename = "./charge-density.xml"
         !
         CALL iotk_link( iunpun, "CHARGE-DENSITY", TRIM(filename), &
                                  CREATE=.FALSE., BINARY=.TRUE. )
         !
!-------------------------------------------------------------------------------
! ... BAND_STRUCTURE_INFO
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_bands_info(  num_k_points, natomwfc, &
                                       nbnd, nbnd, nbnd, &
                                       nspin, nelec,NINT(nelup),NINT(neldw), &
                                       "Hartree", "2 pi / a", &
                                       ef=ef/e2, two_fermi_energies=two_fermi_energies ,&
                                       ef_up=ef_up/e2, ef_down=ef_dw/e2, noncolin=noncolin )
         !
!-------------------------------------------------------------------------------
! ... EIGENVALUES
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_bands_pw( nbnd, num_k_points, nspin, xk, wk, wg,et/e2,"Hartree" , lkpoint_dir , iunout , dirname )
         !
         !
         IF (.NOT. lkpoint_dir ) CALL iotk_close_write( iunout )
         !
!-------------------------------------------------------------------------------
! ... EIGENVECTORS
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "EIGENVECTORS" )
         !
         CALL iotk_write_dat  ( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwx_g )
         !
      END IF
      !
      k_points_loop2: DO ik = 1, num_k_points
         !
         IF ( ionode ) THEN
            !
            CALL iotk_write_begin( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            ! ... G+K vectors
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_GK-VECTORS", ngk_g(ik) )
            !
            IF ( lwfc ) THEN
               !
               filename = qexml_wfc_filename( ".", 'gkvectors', ik, DIR=lkpoint_dir )
               !
               CALL iotk_link( iunpun, "GK-VECTORS", &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               filename = qexml_wfc_filename( dirname, 'gkvectors', ik, &
                                         DIR=lkpoint_dir )
            END IF
            !
         END IF
         !
         IF ( lwfc ) THEN
            !
            IF ( .NOT. smallmem ) CALL write_gk( iunout, ik, filename )
            !
            CALL write_this_wfc ( iunout, ik )
            !
         END IF
         !
         IF ( ionode ) THEN
            !
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
            !
         END IF
         !
      END DO k_points_loop2
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_end( iunpun, "EIGENVECTORS" )
         !
         CALL qexml_closefile( 'write', IERR=ierr)
         !
         !
         CALL delete_if_present( TRIM( dirname ) // '/' // TRIM( xmlpun ) // '.bck' )
         !
      END IF
      !
      DEALLOCATE ( igk_l2g )
      DEALLOCATE ( igk_l2g_kdip )
      !
!-------------------------------------------------------------------------------
! ... CHARGE-DENSITY FILES
!-------------------------------------------------------------------------------
      !
      ! ... also writes rho%ns if lda+U and rho%bec if PAW
      !
      IF ( lrho ) CALL write_rho( rho, nspin )
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( mill_g )
      DEALLOCATE( ngk_g )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'pw_writefile ', 'cannot save history', ierr )
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE write_gk( iun, ik, filename )
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          INTEGER,            INTENT(IN) :: iun, ik
          CHARACTER(LEN=256), INTENT(IN) :: filename
          !
          INTEGER, ALLOCATABLE :: igwk(:,:)
          INTEGER, ALLOCATABLE :: itmp(:)
          !
          !
          ALLOCATE( igwk( npwx_g, nkstot ) )
          !
          igwk(:,ik) = 0
          !
          ALLOCATE( itmp( npw_g ) )
          !
          itmp = 0
          !
          IF ( ik >= iks .AND. ik <= ike ) THEN
             !
             DO ig = 1, ngk(ik-iks+1)
                !
                itmp(igk_l2g(ig,ik-iks+1)) = igk_l2g(ig,ik-iks+1)
                !
             END DO
             !
          END IF
          !
          CALL mp_sum( itmp, inter_pool_comm )
          CALL mp_sum( itmp, intra_pool_comm )
          !
          ngg = 0
          !
          DO ig = 1, npw_g
             !
             if ( itmp(ig) == ig ) THEN
                !
                ngg = ngg + 1
                !
                igwk(ngg,ik) = ig
                !
             END IF
             !
          END DO
          !
          DEALLOCATE( itmp )
          !
          IF ( ionode ) THEN
             !
             CALL iotk_open_write( iun, FILE = TRIM( filename ), &
                                   ROOT="GK-VECTORS", BINARY = .TRUE. )
             !
             CALL iotk_write_dat( iun, "NUMBER_OF_GK-VECTORS", ngk_g(ik) )
             CALL iotk_write_dat( iun, "MAX_NUMBER_OF_GK-VECTORS", npwx_g )
             CALL iotk_write_dat( iun, "GAMMA_ONLY", gamma_only )
             !
             CALL iotk_write_attr ( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
             CALL iotk_write_dat( iun, "K-POINT_COORDS", xk(:,ik), ATTR = attr )
             !
             CALL iotk_write_dat( iun, "INDEX", igwk(1:ngk_g(ik),ik) )
             CALL iotk_write_dat( iun, "GRID", mill_g(1:3,igwk(1:ngk_g(ik),ik)), &
                                  COLUMNS = 3 )
             !
             CALL iotk_close_write( iun )
             !
          END IF
          !
          DEALLOCATE( igwk )
          !
        END SUBROUTINE write_gk
        !
        !--------------------------------------------------------------------
        SUBROUTINE write_this_wfc ( iun, ik )
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN) :: iun, ik
          CHARACTER(LEN=256)  :: filename
          INTEGER :: ispin,ik_eff
          !
          ! ... wavefunctions - do not read if already in memory (nsk==1)
          ! ...                 read only if on this pool (iks <= ik <= ike )
          !
          IF ( ( nks > 1 ) .AND. ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
              CALL get_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
          END IF
          !
          IF ( nspin == 2 ) THEN
             !
             ! ... beware: with pools, isk(ik) has the correct value for 
             ! ... all k-points only on first pool (ionode proc is ok)
             !
             ispin = isk(ik)
             !
             IF ( ionode ) THEN
                !
                filename = qexml_wfc_filename( ".", 'evc', ik, ispin, &
                                                     DIR=lkpoint_dir )
                !
                CALL iotk_link( iunpun, "WFC" // TRIM( iotk_index (ispin) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                !
                filename = qexml_wfc_filename( dirname, 'evc', ik, ispin, & 
                                           DIR=lkpoint_dir )
                !
             END IF
             !
             CALL write_wfc( iunout, ik, nkstot, kunit, ispin, nspin, &
                             evc, npw_g, gamma_only, nbnd, igk_l2g_kdip(:,ik-iks+1),   &
                             ngk(ik-iks+1), filename, 1.D0, &
                             ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
             !
             ik_eff = ik + num_k_points
             !
             ispin = isk(ik_eff)
             !
             ! ... LSDA: now read minority wavefunctions (if not already
             ! ... in memory and if they are on this pool)
             !
             IF ( ( nks > 1 ) .AND. ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
                !
                CALL get_buffer ( evc, nwordwfc, iunwfc, (ik_eff-iks+1) )
                !
             END IF
             !
             IF ( ionode ) THEN
                !
                filename = qexml_wfc_filename( ".", 'evc', ik, ispin, &
                                                     DIR=lkpoint_dir )
                !
                CALL iotk_link( iunpun, "WFC"//TRIM( iotk_index( ispin ) ), &
                                filename, CREATE = .FALSE., BINARY = .TRUE. )
                !
                filename = qexml_wfc_filename( dirname, 'evc', ik, ispin, &
                            DIR=lkpoint_dir )
                !
             END IF
             !
             CALL write_wfc( iunout, ik_eff, nkstot, kunit, ispin, nspin, &
                             evc, npw_g, gamma_only, nbnd, igk_l2g_kdip(:,ik_eff-iks+1), &
                             ngk(ik_eff-iks+1), filename, 1.D0, &
                             ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
             !
          ELSE
             !
             IF ( noncolin ) THEN
                !
                DO ipol = 1, npol
                   !
                   IF ( ionode ) THEN
                      !
                      filename = qexml_wfc_filename( ".", 'evc', ik, ipol, &
                                              DIR=lkpoint_dir )
                      !
                      CALL iotk_link(iunpun,"WFC"//TRIM(iotk_index(ipol)), &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
                      !
                      filename = qexml_wfc_filename( dirname, 'evc', ik, ipol, &
                             DIR=lkpoint_dir)
                      !
                   END IF
                   !
                   ! TEMP  spin-up and spin-down spinor components are written
                   ! TEMP  to different files, like in LSDA - not a smart way
                   !
                   nkl=(ipol-1)*npwx+1
                   nkr= ipol   *npwx
                   CALL write_wfc( iunout, ik, nkstot, kunit, ipol, npol,   &
                                   evc(nkl:nkr,:), npw_g, gamma_only, nbnd, &
                                   igk_l2g_kdip(:,ik-iks+1), ngk(ik-iks+1), &
                                   filename, 1.D0, &
                                   ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
                   !
                END DO
                !
             ELSE
                !
                ispin = 1
                !
                IF ( ionode ) THEN
                   !
                   filename = qexml_wfc_filename( ".", 'evc', ik, DIR=lkpoint_dir )
                   !
                   CALL iotk_link( iunpun, "WFC", filename, &
                                   CREATE = .FALSE., BINARY = .TRUE. )
                   !
                   filename =qexml_wfc_filename( dirname, 'evc', ik, &
                                                      DIR=lkpoint_dir )
                   !
                END IF
                !
                CALL write_wfc( iunout, ik, nkstot, kunit, ispin, nspin, &
                                evc, npw_g, gamma_only, nbnd,            &
                                igk_l2g_kdip(:,ik-iks+1),                &
                                ngk(ik-iks+1), filename, 1.D0, &
                                ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
                !
             END IF
             !
          END IF
          !
       END SUBROUTINE write_this_wfc
       !
    END SUBROUTINE pw_writefile
    !
    !------------------------------------------------------------------------
    SUBROUTINE pw_readfile( what, ierr )
      !------------------------------------------------------------------------
      !
      USE io_rho_xml,    ONLY : read_rho
      USE scf,           ONLY : rho
      USE lsda_mod,      ONLY : nspin
      USE mp_bands,      ONLY : intra_bgrp_comm
      USE mp,            ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: what
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256) :: dirname
      CHARACTER(LEN=80)  :: errmsg
      LOGICAL            :: lcell, lpw, lions, lspin, linit_mag, &
                            lxc, locc, lbz, lbs, lwfc, lheader,          &
                            lsymm, lrho, lefield, ldim, &
                            lef, lexx, lesm
      !
      INTEGER            :: tmp
      !
      ierr = 0
      !
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      ! ... look for an empty unit
      !
      CALL iotk_free_unit( iunout, ierr )
      !
      CALL errore( 'pw_readfile', &
                   'no free units to read wavefunctions', ierr )
      !
      lheader = .NOT. qexml_version_init
      !
      ldim    = .FALSE.
      lcell   = .FALSE.
      lpw     = .FALSE.
      lions   = .FALSE.
      lspin   = .FALSE.
      linit_mag = .FALSE.
      lxc     = .FALSE.
      locc    = .FALSE.
      lbz     = .FALSE.
      lbs     = .FALSE.
      lwfc    = .FALSE.
      lsymm   = .FALSE.
      lrho    = .FALSE.
      lefield = .FALSE.
      lef     = .FALSE.
      lexx    = .FALSE.
      lesm    = .FALSE.
      !
      SELECT CASE( what )
      CASE( 'header' )
         !
         lheader = .TRUE.
         !
      CASE( 'dim' )
         !
         ldim = .TRUE.
         lbz  = .TRUE.
         !
      CASE( 'pseudo' )
         !
         lions = .TRUE.
         !
      CASE( 'config' )
         !
         lcell = .TRUE.
         lions = .TRUE.
         !
      CASE( 'wave' )
         !
         lpw   = .TRUE.
         lwfc  = .TRUE.
         !
      CASE( 'nowavenobs' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag   = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         !
      CASE( 'nowave' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag   = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lbs     = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         !
      CASE( 'all' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag  = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lbs     = .TRUE.
         lwfc    = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         lrho    = .TRUE.
         !
      CASE( 'reset' )
         !
         lcell_read   = .FALSE.
         lpw_read     = .FALSE.
         lions_read   = .FALSE.
         lspin_read   = .FALSE.
         lstarting_mag_read   = .FALSE.
         lxc_read     = .FALSE.
         locc_read    = .FALSE.
         lbz_read     = .FALSE.
         lbs_read     = .FALSE.
         lwfc_read    = .FALSE.
         lsymm_read   = .FALSE.
         lefield_read = .FALSE.
         !
      CASE( 'ef' )
         !
         lef        = .TRUE.
         !
      CASE( 'exx' )
         !
         lexx       = .TRUE.
         !
      CASE( 'esm' )
         !
         lesm       = .TRUE.
         !
      CASE DEFAULT
         !
         CALL errore( 'pw_readfile', 'unknown case '//TRIM(what), 1 )
         !
      END SELECT
      !
      IF ( .NOT. lheader .AND. .NOT. qexml_version_init) &
         CALL errore( 'pw_readfile', 'qexml version not set', 71 )
      !
      IF (  ionode ) THEN
         !
         CALL qexml_init( iunpun )
         CALL qexml_openfile( TRIM( dirname ) // '/' // TRIM( xmlpun ), &
                              'read', BINARY = .FALSE., IERR = ierr  )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr /=0 ) THEN
         errmsg='error opening xml data file'
         GOTO 100
      END IF
      !
      IF ( lheader ) THEN
         !
         CALL read_header( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading header of xml data file'
            GOTO 100
         END IF
         !
      ENDIF
      !
      IF ( ldim ) THEN
         !
         CALL read_dim( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading dimensions in xml data file'
            GOTO 100
         END IF
         !
      ENDIF
      !
      IF ( lcell ) THEN
         !
         CALL read_cell( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading cell info in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lpw ) THEN
         !
         CALL read_planewaves( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading plane-wave info in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lions ) THEN
         !
         CALL read_ions( dirname, ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading info on ions in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lspin ) THEN
         !
         CALL read_spin( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading spin in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF (linit_mag) THEN
         !
         CALL read_magnetization( ierr ) 
         IF ( ierr > 0 ) THEN
            errmsg='error reading magnetization in xml data file'
            GOTO 100
         END IF
        !
      ENDIF
      IF ( lxc ) THEN
         !
         CALL read_xc( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading XC functional in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( locc ) THEN
         !
         CALL read_occupations( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading occupation numbers in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lbz ) THEN
         !
         CALL read_brillouin_zone( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading Brillouin Zone in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lbs ) THEN
         !
         CALL read_band_structure( dirname, ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading band structure in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lwfc ) THEN
         !
         CALL read_wavefunctions( dirname, ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading wavefunctions in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lsymm ) THEN
         !
         CALL read_symmetry( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading symmetry in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lefield ) THEN
         !
         CALL read_efield( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading electric fields in xml data file'
            GOTO 100
         END IF
         !
      END IF

      IF ( lrho ) THEN
         !
         ! ... to read the charge-density we use the routine from io_rho_xml 
         ! ... it also reads ns for ldaU and becsum for PAW
         !
         CALL read_rho( rho, nspin )
         !
      END IF

      IF ( lef ) THEN
         !
         CALL read_ef( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading Fermi energy and number of electrons in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lexx ) THEN
         !
         CALL read_exx( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading hybrid functional in xml data file'
            GOTO 100
         END IF
         !
      END IF
      IF ( lesm ) THEN
         !
         CALL read_esm( ierr )
         IF ( ierr > 0 ) THEN
            errmsg='error reading ESM restart data in xml data file'
            GOTO 100
         END IF
         !
      END IF
      !
      IF (ionode) THEN
         !
         CALL qexml_closefile( 'read', IERR=ierr)
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr > 0 ) THEN
         errmsg='error closing xml data file'
         GOTO 100
      END IF
      !

      RETURN
      !
      ! uncomment to continue execution after an error occurs
      ! 100 IF (ionode) THEN
      !        CALL qexml_closefile( 'read', IERR=tmp)
      !     ENDIF
      !     RETURN
      ! comment to continue execution after an error occurs


100   CALL errore('pw_readfile',TRIM(errmsg),ierr)


      !
    END SUBROUTINE pw_readfile
    !
#if defined  __XSD
     !------------------------------------------------------------------------
     SUBROUTINE pw_readschema_file(ierr, restart_output, restart_input, restart_parallel_info, restart_general_info)
      !------------------------------------------------------------------------
      USE qes_types_module,     ONLY : input_type, output_type, general_info_type, parallel_info_type    
      USE qexsd_reader_module,  ONLY : qexsd_get_output, qexsd_get_input, qexsd_get_general_info, &
                                       qexsd_get_parallel_info
      !
      USE qes_libs_module,      ONLY : qes_write_input, qes_write_output, qes_write_parallel_info, &
                                       qes_write_general_info
      IMPLICIT NONE 
      ! 
      INTEGER                                            :: ierr, iotk_err  
      TYPE( output_type ),OPTIONAL,      INTENT(OUT)   :: restart_output
      TYPE(input_type),OPTIONAL,         INTENT(OUT)   :: restart_input
      TYPE(parallel_info_type),OPTIONAL, INTENT(OUT)   :: restart_parallel_info
      TYPE(general_info_type ),OPTIONAL, INTENT(OUT)   :: restart_general_info
      ! 
      LOGICAL                                   :: found
      CHARACTER(LEN=80)                         :: errmsg = "" 
      CHARACTER(LEN=256)                        :: dirname
      !  
      ! 
      ierr = 0 
      ! 
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      CALL iotk_free_unit( iunpun_xsd, iotk_err )
      !
      CALL errore( 'pw_readschema_file', &
                   'no free units to read xsd output', iotk_err )
      CALL qexsd_init_schema( iunpun_xsd )
      CALL iotk_open_read( iunpun_xsd, TRIM(dirname)//'/'//TRIM(xmlpun_schema))
      !
      IF ( PRESENT ( restart_general_info ) ) THEN 
         CALL qexsd_get_general_info ( iunpun_xsd, restart_general_info , found)
         IF (.NOT. found ) ierr = ierr + 1
         IF ( ierr /=0 ) THEN
            errmsg='error header of xml data file'
            GOTO 100
         END IF
         ! CALL qes_write_general_info( 82, restart_general_info) 
      END IF 
      ! 
      IF ( PRESENT ( restart_parallel_info ) ) THEN 
         CALL qexsd_get_parallel_info ( iunpun_xsd, restart_parallel_info, found ) 
         IF ( .NOT. found ) THEN 
            ierr = ierr + 1  
            errmsg='error parallel_info  of xsd data file' 
            GOTO 100
         END IF
         ! CALL qes_write_parallel_info ( 82, restart_parallel_info )
      END IF  
      ! 
      !
      IF ( PRESENT ( restart_input ) ) THEN   
         CALL qexsd_get_input ( iunpun_xsd, restart_input, found ) 
         IF ( .NOT. found ) THEN 
            ierr = ierr + 1  
            errmsg='error input of xsd data file' 
            GOTO 100
         END IF
         ! CALL qes_write_input( 82, restart_input )  
      END IF 
      ! 
      IF ( PRESENT ( restart_output ) ) THEN 
         CALL qexsd_get_output ( iunpun_xsd, restart_output, found ) 
         IF ( .NOT. found ) THEN 
            ierr = ierr + 1 
            errmsg = 'error output of xsd data file' 
            GOTO 100 
         END IF 
         ! CALL qes_write_output ( 82, restart_output ) 
      END IF 
      CALL iotk_close_read (iunpun_xsd)
      RETURN
 100  CALL errore('pw_readschemafile',TRIM(errmsg),ierr)
    END SUBROUTINE pw_readschema_file
    !  
    !------------------------------------------------------------------------
    SUBROUTINE init_vars_from_schema( what, ierr, output_obj, input_obj, par_info, gen_info )
      !------------------------------------------------------------------------
      !
      USE control_flags,        ONLY : twfcollect
      USE io_rho_xml,           ONLY : read_rho
      USE scf,                  ONLY : rho
      USE lsda_mod,             ONLY : nspin
      USE mp_world,             ONLY : mpime
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE mp,                   ONLY : mp_sum, mp_barrier
      USE qes_types_module,     ONLY : input_type, output_type, general_info_type, parallel_info_type    
!      !
      IMPLICIT NONE
!      !
      CHARACTER(LEN=*), INTENT(IN)           :: what
      TYPE ( output_type), INTENT(IN)        :: output_obj
      TYPE ( input_type ), INTENT(IN)        :: input_obj
      TYPE ( parallel_info_type), INTENT(IN) :: par_info
      TYPE ( general_info_type ), INTENT(IN) :: gen_info
      INTEGER,INTENT (OUT)                   :: ierr 
      !
      CHARACTER(LEN=256) :: dirname
      LOGICAL            :: lcell, lpw, lions, lspin, linit_mag, &
                            lxc, locc, lbz, lbs, lwfc, lheader,          &
                            lsymm, lrho, lefield, ldim, &
                            lef, lexx, lesm
      !
      LOGICAL            :: need_qexml, found, electric_field_ispresent
      INTEGER            :: tmp, iotk_err 
      
      !    
      !
      ierr = 0 
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      !
      !
      !
      ldim    = .FALSE.
      lcell   = .FALSE.
      lpw     = .FALSE.
      lions   = .FALSE.
      lspin   = .FALSE.
      linit_mag = .FALSE.
      lxc     = .FALSE.
      locc    = .FALSE.
      lbz     = .FALSE.
      lbs     = .FALSE.
      lwfc    = .FALSE.
      lsymm   = .FALSE.
      lrho    = .FALSE.
      lefield = .FALSE.
      lef     = .FALSE.
      lexx    = .FALSE.
      lesm    = .FALSE.
      !
     
         
      SELECT CASE( what )
      CASE( 'header' )
         !
         lheader = .TRUE.
         need_qexml = .TRUE.
         !
      CASE ( 'wf_collect' ) 
         ! 
         twfcollect = input_obj%control_variables%wf_collect
         !
      CASE( 'dim' )
         !
         ldim =       .TRUE.
         lbz  =       .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'pseudo' )
         !
         lions = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'config' )
         !
         lcell = .TRUE.
         lions = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'rho' )
         !
         lrho  = .TRUE.
         !
      CASE( 'wave' )
         !
         lpw   = .TRUE.
         lwfc  = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'nowavenobs' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag   = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         need_qexml = .TRUE.
                                                                  
      CASE( 'nowave' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag   = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lbs     = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'all' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag  = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lbs     = .TRUE.
         lwfc    = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         lrho    = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'reset' )
         !
         lcell_read   = .FALSE.
         lpw_read     = .FALSE.
         lions_read   = .FALSE.
         lspin_read   = .FALSE.
         lstarting_mag_read   = .FALSE.
         lxc_read     = .FALSE.
         locc_read    = .FALSE.
         lbz_read     = .FALSE.
         lbs_read     = .FALSE.
         lwfc_read    = .FALSE.
         lsymm_read   = .FALSE.
         lefield_read = .FALSE.
         !
      CASE( 'ef' )
         !
         lef        = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'exx' )
         !
         lexx       = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'esm' )
         !
         lesm       = .TRUE.
         need_qexml = .TRUE.
         !
      END SELECT
      !
   
      electric_field_ispresent = input_obj%electric_field_ispresent 
      !
      IF ( lheader ) THEN 
         CALL readschema_header( gen_info )
      END IF 
      IF ( ldim ) THEN
         !         ! 

         ! 
         CALL readschema_dim(par_info, input_obj%atomic_species, output_obj%atomic_structure, output_obj%symmetries, &
                             output_obj%basis_set, output_obj%band_structure, input_obj) 
                                                                                                           
      ENDIF
      !
      IF ( lcell .AND. ( .NOT. lcell_read) ) THEN
         CALL readschema_cell( output_obj%atomic_structure,  input_obj )
      END IF
      !
      IF ( lpw ) THEN
         twfcollect = input_obj%control_variables%wf_collect
         CALL readschema_planewaves( output_obj%basis_set) 
      END IF
      IF ( lions ) THEN
         CALL readschema_ions( output_obj%atomic_structure, output_obj%atomic_species, input_obj, dirname)
      END IF
      IF ( lspin ) THEN

         CALL readschema_spin( output_obj%magnetization )
      END IF
      IF (linit_mag) THEN
         CALL readschema_magnetization (  output_obj%band_structure,  input_obj%atomic_species, input_obj)
      END IF
      IF ( lxc ) THEN
         CALL readschema_xc (  input_obj%atomic_species, output_obj%dft )
      END IF
      IF ( locc ) THEN
         CALL readschema_occupations( input_obj, output_obj%band_structure )
      END IF
      IF ( lbz ) THEN
         CALL readschema_brillouin_zone( input_obj%k_points_IBZ, input_obj%bands%occupations, &
                                         output_obj%symmetries,  output_obj%band_structure )
      END IF
      IF ( lbs ) THEN
         CALL readschema_band_structure( output_obj%band_structure )
      END IF
      IF ( lwfc ) THEN
         !
         IF (input_obj%control_variables%wf_collect) THEN 
            twfcollect =  input_obj%control_variables%wf_collect 
            CALL read_collected_to_evc(dirname, ierr ) 
         END IF
      END IF
      IF ( lsymm ) THEN
         CALL readschema_symmetry ( output_obj%symmetries, output_obj%basis_set, input_obj%symmetry_flags)
      END IF
      IF ( lefield ) THEN
         IF ( electric_field_ispresent ) THEN 
             CALL readschema_efield( input_obj%electric_field )
         ELSE 
             CALL readschema_efield()
         END IF
      END IF
      IF ( lrho ) THEN
         !
         ! ... to read the charge-density we use the routine from io_rho_xml 
         ! ... it also reads ns for ldaU and becsum for PAW
         !
         CALL read_rho( rho, nspin )
         !
      END IF
      IF ( lef ) THEN
               CALL readschema_ef ( output_obj%band_structure) 
         !
      END IF
      !
      IF ( lexx .AND. input_obj%dft%hybrid_ispresent  ) CALL readschema_exx ( input_obj%dft%hybrid )
      !
      IF ( lesm .AND. input_obj%boundary_conditions_ispresent ) THEN 
         IF ( input_obj%boundary_conditions%esm_ispresent ) &
                        CALL readschema_esm ( input_obj%boundary_conditions%esm) 
      END IF 
      !
      !
      !
      RETURN
      !
      ! uncomment to continue execution after an error occurs
      ! 100 IF (ionode .AND. need_qexml) THEN
      !        CALL qexml_closefile( 'read', IERR=tmp)
      !     ENDIF
      !     RETURN
      ! comment to continue execution after an error occurs
      !
    END SUBROUTINE init_vars_from_schema
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_header( ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the format version of the current xml datafile
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      !
      ierr = 0
      !
      IF ( qexml_version_init ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_header( FORMAT_VERSION = qexml_version, ierr = ierr )
         !
         qexml_version_init = .TRUE.
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr /=0 ) RETURN
      !
      CALL mp_bcast( qexml_version,       ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init,  ionode_id, intra_image_comm )
      !
      !
    END SUBROUTINE read_header
    !
#ifdef __XSD
    !-------------------------------------------------------------------------------
    SUBROUTINE readschema_header (gen_info_obj) 
    !-------------------------------------------------------------------------------
    ! 
    USE io_files,            ONLY: qexsd_fmt, qexsd_version, qexsd_init
    USE qes_types_module,    ONLY: general_info_type
    IMPLICIT NONE 
    !
    TYPE (general_info_type ),INTENT(IN)  :: gen_info_obj   
    ! 
    IF ( qexsd_init ) RETURN 
    qexsd_fmt = TRIM (gen_info_obj%xml_format%NAME)
    qexsd_version = TRIM ( gen_info_obj%xml_format%VERSION)
    qexsd_init = .TRUE. 
    END SUBROUTINE readschema_header 
    ! 
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_dim( ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine collects array dimensions from various sections
      ! ... plus with some other variables needed for array allocation 
      !
      USE ions_base,        ONLY : nat, nsp
      USE symm_base,        ONLY : nsym
      USE gvect,            ONLY : ngm_g, ecutrho
      USE fft_base,         ONLY : dfftp
      USE gvecs,            ONLY : ngms_g, dual
      USE fft_base,         ONLY : dffts
      USE lsda_mod,         ONLY : lsda
      USE noncollin_module, ONLY : noncolin
      USE ktetra,           ONLY : ntetra
      USE klist,            ONLY : nkstot, nelec
      USE wvfct,            ONLY : nbnd, npwx
      USE gvecw,            ONLY : ecutwfc
      USE control_flags,    ONLY : gamma_only
      USE mp_pools,         ONLY : kunit
      USE mp_global,        ONLY : nproc_file, nproc_pool_file, &
                                   nproc_image_file, ntask_groups_file, &
                                   nproc_bgrp_file, nproc_ortho_file
      !
      IMPLICIT NONE
      !
      !CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER  :: npwx_
      LOGICAL  :: found, found2
      CHARACTER(iotk_attlenx)  :: attr
      !
      !
      ! ... first the entire CELL section is read
      ! ... 
      ierr=0
      !
      CALL read_cell( ierr )
      IF ( ierr /= 0) GOTO 100
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_ions( NAT=nat, NSP=nsp, IERR=ierr)
         IF ( ierr /= 0) GOTO 100
         !
         CALL qexml_read_symmetry(NSYM=nsym, FOUND=found, IERR=ierr)
         IF ( ierr /= 0) GOTO 100
         !
         IF ( .NOT. found ) THEN
            !
            nsym = 1
            !
         ENDIF
         !
         CALL qexml_read_planewaves(  ECUTWFC=ecutwfc, ECUTRHO=ecutrho, NPWX=npwx_,GAMMA_ONLY=gamma_only, &
                                      NR1  = dfftp%nr1, NR2  = dfftp%nr2, NR3  = dfftp%nr3, NGM  = ngm_g, &
                                      NR1S = dffts%nr1, NR2S = dffts%nr2, NR3S = dffts%nr3, NGMS = ngms_g, IERR=ierr )
         IF ( ierr /= 0) GOTO 100
         !
         ecutwfc = ecutwfc * e2
         ecutrho = ecutrho * e2
         !
         dual = ecutrho / ecutwfc
         !
         CALL qexml_read_spin( LSDA = lsda, NONCOLIN = noncolin, IERR=ierr )
         IF ( ierr /= 0) GOTO 100
         !
         CALL qexml_read_occ( NTETRA = ntetra, IERR=ierr )
         IF ( ierr /= 0) GOTO 100
         !
         CALL qexml_read_bz( NUM_K_POINTS= nkstot, IERR =  ierr )
         IF ( ierr /= 0) GOTO 100
         !
         IF ( lsda ) nkstot = nkstot * 2
         !
         CALL qexml_read_bands_info( NBND=nbnd, NELEC=nelec, IERR=ierr )
         IF ( ierr /= 0) GOTO 100
         !
         CALL qexml_read_para( KUNIT=kunit, NPROC=nproc_file, NPROC_POOL=nproc_pool_file, &
              NPROC_IMAGE=nproc_image_file, NTASK_GROUPS = ntask_groups_file, &
              NPROC_BGRP=nproc_bgrp_file, NPROC_ORTHO=nproc_ortho_file, FOUND=found, IERR=ierr )
         IF ( ierr /= 0) GOTO 100
         !
         IF ( .NOT. found ) THEN
            !
            kunit = 1
            nproc_file=1
            nproc_pool_file=1
            nproc_image_file=1
            ntask_groups_file=1
            nproc_bgrp_file=1
            nproc_ortho_file=1
            !
         ENDIF
         !
      END IF
      !
100   CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      CALL mp_bcast( nat,        ionode_id, intra_image_comm )
      CALL mp_bcast( nsp,        ionode_id, intra_image_comm )
      CALL mp_bcast( nsym,       ionode_id, intra_image_comm )
      CALL mp_bcast( ecutwfc,    ionode_id, intra_image_comm )
      CALL mp_bcast( ecutrho,    ionode_id, intra_image_comm )
      CALL mp_bcast( dual,       ionode_id, intra_image_comm )
      CALL mp_bcast( npwx_,      ionode_id, intra_image_comm )
      CALL mp_bcast( gamma_only, ionode_id, intra_image_comm )
      CALL mp_bcast( dfftp%nr1,  ionode_id, intra_image_comm )
      CALL mp_bcast( dfftp%nr2,  ionode_id, intra_image_comm )
      CALL mp_bcast( dfftp%nr3,  ionode_id, intra_image_comm )
      CALL mp_bcast( ngm_g,      ionode_id, intra_image_comm )
      CALL mp_bcast( dffts%nr1,  ionode_id, intra_image_comm )
      CALL mp_bcast( dffts%nr2,  ionode_id, intra_image_comm )
      CALL mp_bcast( dffts%nr3,  ionode_id, intra_image_comm )
      CALL mp_bcast( ngms_g,     ionode_id, intra_image_comm )
      CALL mp_bcast( lsda,       ionode_id, intra_image_comm )
      CALL mp_bcast( noncolin,   ionode_id, intra_image_comm )
      CALL mp_bcast( ntetra,     ionode_id, intra_image_comm )
      CALL mp_bcast( nkstot,     ionode_id, intra_image_comm )
      CALL mp_bcast( nelec,      ionode_id, intra_image_comm )
      CALL mp_bcast( nbnd,       ionode_id, intra_image_comm )
      CALL mp_bcast( kunit,      ionode_id, intra_image_comm )
      CALL mp_bcast( nproc_file, ionode_id, intra_image_comm )
      CALL mp_bcast( nproc_pool_file,    ionode_id, intra_image_comm )
      CALL mp_bcast( nproc_image_file,   ionode_id, intra_image_comm )
      CALL mp_bcast( ntask_groups_file,  ionode_id, intra_image_comm )
      CALL mp_bcast( nproc_bgrp_file,    ionode_id, intra_image_comm )
      CALL mp_bcast( nproc_ortho_file,   ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_dim
    ! 
#ifdef __XSD
    !--------------------------------------------------------------------------
    SUBROUTINE readschema_dim(par_info_obj, atomic_species, atomic_structure, symmetries, basis_set, &
                              band_structure, input_obj) 
    ! 
    USE ions_base,        ONLY : nat, nsp
    USE symm_base,        ONLY : nsym
    USE gvect,            ONLY : ngm_g, ecutrho
    USE fft_base,         ONLY : dfftp
    USE gvecs,            ONLY : ngms_g, dual
    USE fft_base,         ONLY : dffts
    USE lsda_mod,         ONLY : lsda
    USE noncollin_module, ONLY : noncolin
    USE ktetra,           ONLY : ntetra
    USE klist,            ONLY : nkstot, nelec
    USE wvfct,            ONLY : nbnd, npwx
    USE gvecw,            ONLY : ecutwfc
    USE control_flags,    ONLY : gamma_only
    USE mp_pools,         ONLY : kunit
    USE mp_global,        ONLY : nproc_file, nproc_pool_file, &
                                 nproc_image_file, ntask_groups_file, &
                                 nproc_bgrp_file, nproc_ortho_file
    !
    USE qes_types_module, ONLY : parallel_info_type, atomic_species_type, atomic_structure_type, &
                                 symmetries_type, basis_set_type, band_structure_type, input_type  
    IMPLICIT NONE 
    !
    TYPE ( parallel_info_type ),INTENT(IN)     :: par_info_obj
    TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_species
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    TYPE ( symmetries_type ),INTENT(IN)        :: symmetries
    TYPE ( basis_set_type ),INTENT(IN)         :: basis_set
    TYPE ( band_structure_type ),INTENT(IN)    :: band_structure 
    TYPE ( input_type ),INTENT(IN)             :: input_obj 
    ! 
    INTEGER                                    :: npwx_
    CALL readschema_cell ( atomic_structure, input_obj ) 
    ! 
    !---------------------------------------------------------------------
    !                                       PARALLEL  DIM 
    !----------------------------------------------------------------------
    nproc_file = par_info_obj%nprocs
    nproc_pool_file = nproc_file/par_info_obj%npool
    nproc_image_file = nproc_file 
    ntask_groups_file = par_info_obj%ntasks
    nproc_bgrp_file = nproc_image_file / par_info_obj%npool / par_info_obj%nbgrp 
    nproc_ortho_file = par_info_obj%ndiag
    !---------------------------------------------------------------------------
    !                                      ATOMS AND SPECIES 
    !--------------------------------------------------------------------------
    nsp = atomic_species%ntyp
    nat = atomic_structure%nat 
    !                                         SIMMETRIES 
    nsym = symmetries%nsym
    !-----------------------------------------------------------------------------
    !                                          BASIS SET 
    !-----------------------------------------------------------------------------
    ecutwfc = basis_set%ecutwfc*e2
    ecutrho = basis_set%ecutrho*e2
    dual = ecutrho/ecutwfc
    npwx_ = basis_set%npwx
    gamma_only= basis_set%gamma_only
    dfftp%nr1 = basis_set%fft_grid%nr1
    dfftp%nr2 = basis_set%fft_grid%nr2          
    dfftp%nr3 = basis_set%fft_grid%nr3
    dffts%nr1 = basis_set%fft_smooth%nr1
    dffts%nr2 = basis_set%fft_smooth%nr2
    dffts%nr3 = basis_set%fft_smooth%nr3
    ngm_g     = basis_set%ngm
    ngms_g    = basis_set%ngms
    !-------------------------------------------------------------------------
    !                                    BAND STRUCTURE  
    !-------------------------------------------------------------------------
    lsda  =    band_structure%lsda
    noncolin = band_structure%noncolin
    nelec =    band_structure%nelec
    nkstot =   band_structure%nks 
    IF ( lsda ) nkstot = nkstot * 2 
    nbnd = band_structure%nbnd
    IF ( input_obj%k_points_IBZ%monkhorst_pack_ispresent ) THEN 
       IF ( TRIM( input_obj%bands%occupations%occupations) == 'tetrahedra' ) THEN 
          ntetra = 6* input_obj%k_points_IBZ%monkhorst_pack%nk1* &
                      input_obj%k_points_IBZ%monkhorst_pack%nk2* &
                      input_obj%k_points_IBZ%monkhorst_pack%nk3 
       END IF 
    ELSE 
       ntetra = 0 
    END IF 
  
    ! 
    END SUBROUTINE readschema_dim
    !
#endif
    !--------------------------------------------------------------------------
    SUBROUTINE read_cell( ierr )
      !------------------------------------------------------------------------
      !
      USE run_info,          ONLY : title
      USE cell_base,         ONLY : ibrav, alat, at, bg, celldm
      USE cell_base,         ONLY : tpiba, tpiba2, omega
      USE cellmd,            ONLY : lmovecell, cell_factor
      USE control_flags,     ONLY : do_makov_payne
      USE martyna_tuckerman, ONLY : do_comp_mt
      USE esm,               ONLY : do_comp_esm
      
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=80) :: bravais_lattice, es_corr
      !
      !
      ierr = 0
      IF ( lcell_read ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_cell( BRAVAIS_LATTICE=bravais_lattice,CELLDM=celldm, ALAT=alat, &
              A1=at(:,1), A2=at(:,2), A3=at(:,3), ES_CORR=es_corr, IERR=ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         SELECT CASE ( TRIM(es_corr))
         CASE ("Makov-Payne")
            do_makov_payne = .true.
            do_comp_mt     = .false.
            do_comp_esm    = .false. 
         CASE ("Martyna-Tuckerman")
            do_makov_payne = .false.
            do_comp_mt     = .true.
            do_comp_esm    = .false.
         CASE ("ESM")
            do_makov_payne = .false.
            do_comp_mt     = .false.
            do_comp_esm    = .true.
         CASE ("None")
            do_makov_payne = .false.
            do_comp_mt     = .false.
            do_comp_esm    = .false.
         END SELECT
         !
         SELECT CASE ( TRIM(bravais_lattice) )
         CASE( "free" )
            ibrav = 0
         CASE( "cubic P (sc)" )
            ibrav = 1
         CASE( "cubic F (fcc)" )
            ibrav = 2
         CASE( "cubic I (bcc)" )
            ibrav = 3
         CASE( "Hexagonal and Trigonal P" )
            ibrav = 4
         CASE( "Trigonal R" )
            ibrav = 5 
         CASE( "Tetragonal P (st)" )
            ibrav = 6
         CASE( "Tetragonal I (bct)" )
            ibrav = 7
         CASE( "Orthorhombic P" )
            ibrav = 8
         CASE( "Orthorhombic base-centered(bco)" )
            ibrav = 9
         CASE( "Orthorhombic face-centered" )
            ibrav = 10
         CASE( "Orthorhombic body-centered" )
            ibrav = 11
         CASE( "Monoclinic P" )
            ibrav = 12
         CASE( "Monoclinic base-centered" )
            ibrav = 13
         CASE( "Triclinic P" )
            ibrav = 14
         CASE DEFAULT
            ibrav = 0
         END SELECT
         !
         ! ... some internal variables
         !
         tpiba  = 2.D0 * pi / alat
         tpiba2 = tpiba**2 
         !
         ! ... to alat units
         !
         at(:,:) = at(:,:) / alat
         !
         CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
         !
         ! ... Generate the reciprocal lattice vectors
         !
         CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
         !
         CALL qexml_read_moving_cell(lmovecell, cell_factor, ierr)
         !
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr > 0 ) RETURN
      !
      CALL mp_bcast( ibrav,     ionode_id, intra_image_comm )
      CALL mp_bcast( alat,      ionode_id, intra_image_comm )
      CALL mp_bcast( celldm,    ionode_id, intra_image_comm )
      CALL mp_bcast( tpiba,     ionode_id, intra_image_comm )
      CALL mp_bcast( tpiba2,    ionode_id, intra_image_comm )
      CALL mp_bcast( omega,     ionode_id, intra_image_comm )
      CALL mp_bcast( at,        ionode_id, intra_image_comm )
      CALL mp_bcast( bg,        ionode_id, intra_image_comm )
      CALL mp_bcast( do_makov_payne, ionode_id, intra_image_comm )
      CALL mp_bcast( do_comp_mt,     ionode_id, intra_image_comm )
      CALL mp_bcast( do_comp_esm,    ionode_id, intra_image_comm )
      CALL mp_bcast( lmovecell, ionode_id, intra_image_comm )
      IF (lmovecell) THEN
         CALL mp_bcast( cell_factor,  ionode_id, intra_image_comm )
      ELSE
         cell_factor=1.0_DP
      END IF
      !
      title = ' '
      !
      lcell_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_cell
    !
#ifdef  __XSD
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_cell(atomic_structure, input_obj  )
    !-----------------------------------------------------------------------
    ! 
    USE run_info,          ONLY : title
    USE cell_base,         ONLY : ibrav, alat, at, bg, celldm
    USE cell_base,         ONLY : tpiba, tpiba2, omega
    USE cellmd,            ONLY : lmovecell, cell_factor
    USE control_flags,     ONLY : do_makov_payne
    USE martyna_tuckerman, ONLY : do_comp_mt
    USE esm,               ONLY : do_comp_esm
    USE qes_types_module,  ONLY : input_type, atomic_structure_type
    !
    IMPLICIT NONE 
    ! 
    TYPE ( atomic_structure_type )            :: atomic_structure 
    TYPE ( input_type )                       :: input_obj
    !
    IF ( lcell_read ) RETURN 
    alat = atomic_structure%alat 
    IF ( atomic_structure%bravais_index_ispresent ) THEN 
       ibrav = atomic_structure%bravais_index 
    ELSE 
       ibrav = 0
    END IF
    celldm = 0.d0
    at(:,1) = atomic_structure%cell%a1
    at(:,2) = atomic_structure%cell%a2
    at(:,3) = atomic_structure%cell%a3
    SELECT CASE  (ibrav ) 
       CASE (1:3) 
          celldm(1) = alat
          celldm(2:6) = 0.d0
       CASE (4) 
          celldm(1) = alat
          celldm(2) = 0.d0
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3)))/alat
          celldm(4:6) = 0.d0
       CASE (5) 
          celldm(1)= alat
          celldm(2:3) = 0.d0
          celldm(4) = DOT_PRODUCT(at(:,1),at(:,2))/(alat**2)
          celldm(5:6) = 0.d0
       CASE (6) 
          celldm(1)= alat 
          celldm(3)= SQRT( DOT_PRODUCT(at(:,2),at(:,2)))/alat
          celldm(2)= 1.d0
          celldm(4:6) = 0.d0
       CASE (7) 
          celldm(1) = alat
          celldm(3) = at(3,3) 
          celldm(2)=0.d0
          celldm(4:6) = 0.d0
       CASE (8)
          celldm(1) = alat
          celldm(2) = SQRT( DOT_PRODUCT (at(:,2),at(:,2)))/alat
          celldm(3) = SQRT( DOT_PRODUCT (at(:,3),at(:,3)))/alat 
          celldm(4:6) = 0.d0
       CASE (9) 
          celldm(1) = alat
          celldm(2) = ABS ( at(2,1)/at(1,1))
          celldm(3) = ABS ( at(3,3)/2.d0/at(1,1))
          celldm(4:6) = 0.d0 
       CASE (10) 
          celldm(1) = alat
          celldm(2) = ABS ( at(2,2)/at(2,1))
          celldm(3) = ABS ( at(3,1)/at(1,1))
          celldm(4:6) = 0.d0
       CASE (11) 
          celldm(1) = alat
          celldm(2) = ABS(at(2,1)/at(1,1))
          celldm(3) = ABS(at(3,1)/at(1,1))
          celldm(4:6) = 0.d0
       CASE (12) 
          celldm(1) = alat 
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(4) = DOT_PRODUCT(at(:,1),at(:,2))/&
                      SQRT(DOT_PRODUCT(at(:,1),at(:,1))*DOT_PRODUCT(at(:,2),at(:,2)))
          celldm(5) =  DOT_PRODUCT(at(:,1),at(:,3))/&
                   SQRT(DOT_PRODUCT(at(:,1),at(:,1))*DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(6) = 0.d0
       CASE (13) 
          celldm(1) = alat
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2)))/(2.d0*at(1,1))
          celldm(3) = ABS (at(3,3)/at(1,3))
          celldm(4) = ATAN(at(2,2)/at(1,2))
          celldm(5:6) = 0.d0
       CASE (14) 
          celldm(1) = alat 
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(4) = DOT_PRODUCT(at(:,3),at(:,2))/SQRT(DOT_PRODUCT(at(:,2),at(:,2))*&
                                                   DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(5) = DOT_PRODUCT(at(:,3),at(:,1))/SQRT(DOT_PRODUCT(at(:,1),at(:,1))*&
                                                   DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(6) = DOT_PRODUCT(at(:,1),at(:,2))/SQRT(DOT_PRODUCT(at(:,2),at(:,2))*&
                                                   DOT_PRODUCT(at(:,1),at(:,1)))
       CASE  default  
          celldm(1) = 1.d0
          IF (alat .GT. 0.d0 ) celldm(1) = alat
          celldm (2:6) = 0.d0
    END SELECT 
    tpiba = 2.d0*PI/alat
    tpiba2= tpiba**2
    omega = ABS (at(1,1)*at(2,2)*at(3,3)+at(1,2)*at(2,3)*at(3,1)+at(1,3)*at(2,1)*at(3,2)-&
                 at(3,1)*at(2,2)*at(1,3)-at(3,2)*at(2,3)*at(1,1)-at(3,3)*at(2,1)*at(1,2))
    at=at / alat
    CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
    IF ( input_obj%boundary_conditions_ispresent )  THEN 
       SELECT CASE ( TRIM( input_obj%boundary_conditions%assume_isolated ))
         CASE ("makov-payne")
            do_makov_payne = .true.
            do_comp_mt     = .false.
            do_comp_esm    = .false. 
         CASE ("martyna-tuckerman")
            do_makov_payne = .false.
            do_comp_mt     = .true.
            do_comp_esm    = .false.
         CASE ("esm")
            do_makov_payne = .false.
            do_comp_mt     = .false.
            do_comp_esm    = .true.
         CASE ("none")
            do_makov_payne = .false.
            do_comp_mt     = .false.
            do_comp_esm    = .false.
       END SELECT 
    END IF         
    !
    title = TRIM(input_obj%control_variables%title)
    SELECT CASE ( TRIM ( input_obj%control_variables%calculation))
       CASE ('vc-relax', 'vc-md' ) 
           lmovecell = .TRUE. 
           IF ( input_obj%cell_control%cell_factor_ispresent ) THEN 
              cell_factor = input_obj%cell_control%cell_factor
           ELSE 
              cell_factor = 1.d0
           END IF 
       CASE default  
           lmovecell = .FALSE. 
    END SELECT 
    lcell_read = .TRUE.  
    !
    END SUBROUTINE readschema_cell
    ! 
#endif
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_ions( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE ions_base, ONLY : nat, nsp, ityp, amass, atm, tau, if_pos
      USE cell_base, ONLY : alat
      USE io_files,  ONLY : psfile, pseudo_dir, pseudo_dir_cur
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: i
      LOGICAL :: exst
      !
      ierr = 0
      IF ( lions_read ) RETURN
      !
      IF ( .NOT. lcell_read ) &
         CALL errore( 'read_ions', 'read cell first', 1 )
      !
      ! this is where PP files should be read from
      !
      pseudo_dir_cur = trimcheck ( dirname ) 
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_ions( NSP=nsp, NAT=nat, ATM=atm, ITYP=ityp, &
                               PSFILE=psfile, AMASS=amass, &
                               TAU=tau, IF_POS=if_pos, PSEUDO_DIR=pseudo_dir, &
                               IERR=ierr )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         DO i = 1, nat
            !
            tau(:,i) = tau(:,i) / alat
            !
         END DO
         !
      END IF
      !
      CALL mp_bcast( nat,    ionode_id, intra_image_comm )
      CALL mp_bcast( nsp,    ionode_id, intra_image_comm )
      CALL mp_bcast( atm,    ionode_id, intra_image_comm )
      CALL mp_bcast( amass,  ionode_id, intra_image_comm )
      CALL mp_bcast( psfile, ionode_id, intra_image_comm )
      CALL mp_bcast( pseudo_dir, ionode_id, intra_image_comm )
      CALL mp_bcast( ityp,   ionode_id, intra_image_comm )
      CALL mp_bcast( tau,    ionode_id, intra_image_comm )
      CALL mp_bcast( if_pos, ionode_id, intra_image_comm )
      !
      lions_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_ions
    !
#ifdef __XSD
    !------------------------------------------------------------------------
    SUBROUTINE readschema_ions( atomic_structure, atomic_species, input_obj, dirname ) 
    !------------------------------------------------------------------------
    ! 
    USE ions_base, ONLY : nat, nsp, ityp, amass, atm, tau, if_pos
    USE cell_base, ONLY : alat
    USE io_files,  ONLY : psfile, pseudo_dir, pseudo_dir_cur
    USE qes_types_module, ONLY: atomic_structure_type, atomic_species_type, input_type 
    ! 
    IMPLICIT NONE 
    ! 
    TYPE ( atomic_structure_type ),INTENT(IN) :: atomic_structure
    TYPE ( atomic_species_type ),INTENT(IN)        :: atomic_species 
    TYPE ( input_type ),INTENT(IN)            :: input_obj 
    CHARACTER(LEN=*), INTENT(IN)              :: dirname
    ! 
    INTEGER                                   :: iat, isp, idx
    CHARACTER(LEN = 3 ),ALLOCATABLE           :: symbols(:) 
    ! 
    IF ( lions_read ) RETURN 
    nat = atomic_structure%nat
    nsp = atomic_species%ntyp
    ALLOCATE ( symbols(nat) ) 
    DO isp = 1, nsp 
       amass(isp) = 0.d0 
       IF (atomic_species%species(isp)%mass_ispresent) amass(isp) = atomic_species%species(isp)%mass
       atm(isp) = TRIM ( atomic_species%species(isp)%name )
       psfile(isp) = TRIM ( atomic_species%species(isp)%pseudo_file) 
    END DO 
    !
    loop_on_atoms:DO iat = 1, nat
       idx = atomic_structure%atomic_positions%atom(iat)%index
       tau(:,idx) = atomic_structure%atomic_positions%atom(iat)%atom 
       symbols(idx)  = TRIM ( atomic_structure%atomic_positions%atom(idx)%name ) 
       loop_on_species:DO isp = 1, nsp
          IF ( TRIM(symbols(idx)) == TRIM (atm(isp))) THEN 
             ityp(iat) = isp 
             exit loop_on_species
          END IF 
       END  DO loop_on_species
    END DO loop_on_atoms
    !
    IF (ALLOCATED(if_pos) ) DEALLOCATE ( if_pos) 
    ALLOCATE (if_pos(3,nat) )
    IF ( input_obj%free_positions_ispresent ) THEN   
       if_pos = input_obj%free_positions%int_mat
    ELSE 
       if_pos = 1
    END IF 
    ! 
    IF ( .NOT. lcell_read .AND. ( atomic_structure%alat_ispresent )) alat = atomic_structure%alat 
    ! 
    pseudo_dir = TRIM(input_obj%control_variables%pseudo_dir)//'/'
    pseudo_dir_cur = TRIM ( dirname)//'/'  
    ! 
    lions_read = .TRUE.
    END SUBROUTINE readschema_ions
    !  
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_symmetry( ierr )
      !------------------------------------------------------------------------
      !
      USE symm_base,       ONLY : nrot, nsym, invsym, s, ft,ftau, irt, t_rev, &
                                  sname, sr, invs, inverse_s, s_axis_to_cart, &
                                  time_reversal, no_t_rev
      USE control_flags,   ONLY : noinv
      USE fft_base,        ONLY : dfftp
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      CHARACTER(iotk_attlenx)  :: attr
      !
      INTEGER  :: i
      LOGICAL  :: found
      !
      ierr = 0
      IF ( lsymm_read ) RETURN
      !
      IF ( .NOT. lpw_read ) &
         CALL errore( 'read_symmetry', 'read planewaves first', 1 )
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_symmetry( NSYM=nsym, NROT=nrot, INVSYM=invsym, NOINV=noinv, &
              TIME_REVERSAL=time_reversal, NO_T_REV=no_t_rev, &
              TRASL=ft, S=s, SNAME=sname, T_REV=t_rev, &
              IRT=irt, FOUND=found, IERR=ierr )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         IF ( .NOT. found ) THEN
            !
            nsym = 1
            s(:,:,nsym) = 0
            s(1,1,nsym) = 1
            s(2,2,nsym) = 1
            s(3,3,nsym) = 1
            sr(:,:,nsym) = DBLE(s(:,:,nsym))
            ftau(:,nsym)= 0
            ft  (:,nsym)= 0.0_DP
            sname(nsym) = 'identity'
            do i = 1, SIZE( irt, 2 )
               irt(nsym,i) = i
            end do
            invsym = .FALSE.
            noinv=.FALSE.
            t_rev(nsym) = 0
            invs(1)=1
            time_reversal=.TRUE.
            no_t_rev=.FALSE.
            !
         ELSE
            !
            DO i = 1, nsym
               !
               ftau(1,i) = NINT( ft(1,i)*DBLE( dfftp%nr1 ) )
               ftau(2,i) = NINT( ft(2,i)*DBLE( dfftp%nr2 ) )
               ftau(3,i) = NINT( ft(3,i)*DBLE( dfftp%nr3 ) )
               !
            END DO
            !
            ! indices of inverse operations and matrices in cartesian axis
            ! are not saved to disk (maybe they should), are recalculated here 
            !
            CALL inverse_s ()
            CALL s_axis_to_cart ()
            !
         END IF
         !
         !
      END IF
      !
      CALL mp_bcast( nsym,   ionode_id, intra_image_comm )
      CALL mp_bcast( nrot,   ionode_id, intra_image_comm )
      CALL mp_bcast( invsym, ionode_id, intra_image_comm )
      CALL mp_bcast( noinv,  ionode_id, intra_image_comm )
      CALL mp_bcast( time_reversal,  ionode_id, intra_image_comm )
      CALL mp_bcast( no_t_rev,       ionode_id, intra_image_comm )
      CALL mp_bcast( s,      ionode_id, intra_image_comm )
      CALL mp_bcast( ftau,   ionode_id, intra_image_comm )
      CALL mp_bcast( ft,     ionode_id, intra_image_comm )
      CALL mp_bcast( sname,  ionode_id, intra_image_comm )
      CALL mp_bcast( irt,    ionode_id, intra_image_comm )
      CALL mp_bcast( t_rev,  ionode_id, intra_image_comm )
      CALL mp_bcast( invs,   ionode_id, intra_image_comm )
      CALL mp_bcast( sr,     ionode_id, intra_image_comm )
      !
      lsymm_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_symmetry
    !
#ifdef __XSD
    !------------------------------------------------------------------------
    SUBROUTINE readschema_symmetry ( symms_obj, basis_obj  , symm_flags_obj) 
    !------------------------------------------------------------------------
      ! 
      USE symm_base,       ONLY : nrot, nsym, invsym, s, ft,ftau, irt, t_rev, &
                                 sname, sr, invs, inverse_s, s_axis_to_cart, &
                                 time_reversal, no_t_rev
      USE control_flags,   ONLY : noinv
      USE fft_base,        ONLY : dfftp
      USE qes_types_module,ONLY : symmetries_type, symmetry_type, basis_type
      ! 
      IMPLICIT NONE   
      ! 
      TYPE ( symmetries_type )               :: symms_obj 
      TYPE ( basis_set_type )                :: basis_obj
      TYPE ( symmetry_flags_type )           :: symm_flags_obj
      INTEGER                                :: isym 
      ! 
      IF (lsymm_read ) RETURN 
      nrot = symms_obj%nrot 
      nsym = symms_obj%nsym
      ! 
      noinv = symm_flags_obj%noinv 
      no_t_rev = symm_flags_obj%no_t_rev 
      !  
      invsym = .FALSE. 
      DO isym = 1, nrot
        s(:,:,isym) = symms_obj%symmetry(isym)%rotation%mat 
        sname(isym) = TRIM ( symms_obj%symmetry(isym)%info%name )  
        IF ( (TRIM(sname(isym)) == "inversion") .AND. (isym .LE. nsym) ) invsym = .TRUE.
        IF ( symms_obj%symmetry(isym)%fractional_translation_ispresent .AND. (isym .LE. nsym) ) THEN
           ft(1:3,isym)  =  symms_obj%symmetry(isym)%fractional_translation(1:3) 
           ftau(1,isym) = NINT( ft(1,isym)*DBLE( basis_obj%fft_grid%nr1 ) )
           ftau(2,isym) = NINT( ft(2,isym)*DBLE( basis_obj%fft_grid%nr2 ) )
           ftau(3,isym) = NINT( ft(3,isym)*DBLE( basis_obj%fft_grid%nr3 ) )
        END IF
        IF ( symms_obj%symmetry(isym)%info%time_reversal_ispresent ) THEN  
           IF (symms_obj%symmetry(isym)%info%time_reversal) THEN 
              t_rev( isym ) = 1
           ELSE
              t_rev( isym ) = 0 
           END IF 
        END IF
        IF ( symms_obj%symmetry(isym)%equivalent_atoms_ispresent .AND. (isym .LE. nsym) )   &
             irt(isym,:) = symms_obj%symmetry(isym)%equivalent_atoms%index_list(:)
      END DO
      CALL inverse_s()
      CALL s_axis_to_cart() 
      lsymm_read = .TRUE. 
    END SUBROUTINE readschema_symmetry 
    !
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_efield( ierr )
      !----------------------------------------------------------------------
      !
      USE extfield, ONLY : tefield, dipfield, edir, emaxpos, eopreg, eamp
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      LOGICAL                       :: found
      !
      ierr = 0
      IF ( lefield_read ) RETURN
      !
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_efield( TEFIELD=tefield, DIPFIELD=dipfield, EDIR=edir, &
                                 EMAXPOS=emaxpos, EOPREG=eopreg, EAMP=eamp, &
                                 FOUND=found, IERR=ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( (ionode).AND.(.NOT.found) ) THEN
         !
         tefield  = .FALSE.
         dipfield = .FALSE.
         !
      END IF
      !
      CALL mp_bcast( tefield,  ionode_id, intra_image_comm )
      CALL mp_bcast( dipfield, ionode_id, intra_image_comm )
      CALL mp_bcast( edir,     ionode_id, intra_image_comm )
      CALL mp_bcast( emaxpos,  ionode_id, intra_image_comm )
      CALL mp_bcast( eopreg,   ionode_id, intra_image_comm )
      CALL mp_bcast( eamp,     ionode_id, intra_image_comm )
      !
      lefield_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_efield
    !
#ifdef __XSD
    !---------------------------------------------------------------------------
    SUBROUTINE readschema_efield( efield_obj  ) 
    !---------------------------------------------------------------------------
      !       
      USE extfield, ONLY : tefield, dipfield, edir, emaxpos, eopreg, eamp
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( electric_field_type),OPTIONAL, INTENT(IN)    :: efield_obj
      ! 
      !
      IF (PRESENT (efield_obj) .AND. (TRIM(efield_obj%electric_potential) == 'sawtooth_potential')) THEN 
         tefield = .TRUE. 
         IF ( efield_obj%dipole_correction_ispresent ) THEN 
            dipfield = efield_obj%dipole_correction
         ELSE 
            dipfield = .FALSE. 
         END IF 
         IF ( efield_obj%electric_field_direction_ispresent ) THEN 
            edir = efield_obj%electric_field_direction
         ELSE 
            edir = 3 
         END IF
         IF ( efield_obj%potential_max_position_ispresent ) THEN 
            edir = efield_obj%electric_field_direction
         ELSE 
            emaxpos = 3 
         END IF 
         IF ( efield_obj%potential_decrease_width_ispresent ) THEN 
            eopreg = efield_obj%potential_decrease_width
         ELSE 
            eopreg = 1.d-1
         END IF 
         IF ( efield_obj%electric_field_amplitude_ispresent ) THEN 
            eamp = efield_obj%electric_field_amplitude
         ELSE 
            eamp = 1.d-3
         END IF
     ELSE
         tefield = .FALSE. 
         dipfield = .FALSE. 
     END IF
  END SUBROUTINE readschema_efield  
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_planewaves( ierr )
      !------------------------------------------------------------------------
      !
      USE gvect,           ONLY : ngm_g, ecutrho
      USE gvecs,           ONLY : ngms_g, dual
      USE gvecw,           ONLY : ecutwfc
      USE fft_base,        ONLY : dfftp
      USE fft_base,        ONLY : dffts
      USE wvfct,           ONLY : npwx
      USE control_flags,   ONLY : gamma_only
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER  :: npwx_
      !
      ierr = 0
      IF ( lpw_read ) RETURN
      !
      !
      IF ( ionode ) CALL qexml_read_planewaves( ECUTWFC=ecutwfc, ECUTRHO=ecutrho, NPWX=npwx_, &
                                     GAMMA_ONLY=gamma_only, &
                                     NR1 = dfftp%nr1, NR2 = dfftp%nr2, NR3 = dfftp%nr3, NGM=ngm_g, &
                                     NR1S= dffts%nr1, NR2S= dffts%nr2, NR3S= dffts%nr3, &
                                     NGMS=ngms_g, IERR=ierr )

      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         ecutwfc = ecutwfc * e2
         ecutrho = ecutrho * e2
         !
         dual = ecutrho / ecutwfc
         !
      END IF
      !
      CALL mp_bcast( ecutwfc,    ionode_id, intra_image_comm )
      CALL mp_bcast( ecutrho,    ionode_id, intra_image_comm )
      CALL mp_bcast( dual,       ionode_id, intra_image_comm )
      CALL mp_bcast( npwx_,      ionode_id, intra_image_comm )
      CALL mp_bcast( gamma_only, ionode_id, intra_image_comm )
      CALL mp_bcast( dfftp%nr1,        ionode_id, intra_image_comm )
      CALL mp_bcast( dfftp%nr2,        ionode_id, intra_image_comm )
      CALL mp_bcast( dfftp%nr3,        ionode_id, intra_image_comm )
      CALL mp_bcast( ngm_g,      ionode_id, intra_image_comm )
      CALL mp_bcast( dffts%nr1,       ionode_id, intra_image_comm )
      CALL mp_bcast( dffts%nr2,       ionode_id, intra_image_comm )
      CALL mp_bcast( dffts%nr3,       ionode_id, intra_image_comm )
      CALL mp_bcast( ngms_g,     ionode_id, intra_image_comm )
      !
      lpw_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_planewaves  
    !
#ifdef __XSD
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_planewaves ( basis_set_obj ) 
    !-----------------------------------------------------------------------
    ! 
    USE gvect,           ONLY : ngm_g, ecutrho
    USE gvecs,           ONLY : ngms_g, dual
    USE gvecw,           ONLY : ecutwfc
    USE fft_base,        ONLY : dfftp
    USE fft_base,        ONLY : dffts
    USE wvfct,           ONLY : npwx
    USE control_flags,   ONLY : gamma_only
    USE qes_types_module,ONLY : basis_set_type
    ! 
    IMPLICIT NONE 
    ! 
    TYPE ( basis_set_type )              :: basis_set_obj  
    !
    IF ( lpw_read ) RETURN 
    ecutwfc = basis_set_obj%ecutwfc*e2
    ecutrho = basis_set_obj%ecutrho*e2
    dual = ecutrho/ecutwfc
    !npwx = basis_set_obj%npwx
    gamma_only= basis_set_obj%gamma_only
    dfftp%nr1 = basis_set_obj%fft_grid%nr1
    dfftp%nr2 = basis_set_obj%fft_grid%nr2          
    dfftp%nr3 = basis_set_obj%fft_grid%nr3
    dffts%nr1 = basis_set_obj%fft_smooth%nr1
    dffts%nr2 = basis_set_obj%fft_smooth%nr2
    dffts%nr3 = basis_set_obj%fft_smooth%nr3
    ngm_g     = basis_set_obj%ngm
    ngms_g    = basis_set_obj%ngms
    !
    lpw_read = .TRUE.
    END SUBROUTINE readschema_planewaves 
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_spin( ierr )
      !------------------------------------------------------------------------
      !
      USE spin_orb,         ONLY : lspinorb, domag
      USE lsda_mod,         ONLY : nspin, lsda
      USE noncollin_module, ONLY : noncolin, npol
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      !
      LOGICAL :: found
      !
      ierr = 0
      IF ( lspin_read ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_spin( lsda, noncolin, npol, lspinorb, domag, ierr )
         !
         IF ( lsda ) THEN
            !
            nspin = 2
            !
         ELSE IF ( noncolin ) THEN
            !
            nspin = 4
            !
         ELSE
            !
            nspin = 1
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( lsda,     ionode_id, intra_image_comm )
      CALL mp_bcast( nspin,    ionode_id, intra_image_comm )
      CALL mp_bcast( noncolin, ionode_id, intra_image_comm )
      CALL mp_bcast( npol,     ionode_id, intra_image_comm )
      CALL mp_bcast( lspinorb, ionode_id, intra_image_comm )
      CALL mp_bcast( domag,    ionode_id, intra_image_comm )
      !
      lspin_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_spin
    !
#ifdef __XSD
    !--------------------------------------------------------------------------
    SUBROUTINE readschema_spin( magnetization_obj) 
    !--------------------------------------------------------------------------
      ! 
      USE spin_orb,         ONLY : lspinorb, domag
      USE lsda_mod,         ONLY : nspin, lsda
      USE noncollin_module, ONLY : noncolin, npol
      USE qes_types_module, ONLY : magnetization_type
      USE symm_base,        ONLY : time_reversal
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( magnetization_type ),INTENT(IN)         :: magnetization_obj
      ! 
      IF ( lspin_read ) RETURN
      lspinorb = magnetization_obj%spinorbit 
      domag =   magnetization_obj%do_magnetization 
      lsda  =   magnetization_obj%lsda
      noncolin = magnetization_obj%noncolin  
      IF ( noncolin .AND. domag ) time_reversal = .FALSE.
      IF ( lsda ) THEN  
        nspin = 2
        npol = 1
      ELSE IF (noncolin ) THEN 
        nspin = 4
        npol = 2
      ELSE 
        nspin =1
        npol = 1 
      END IF
      ! 
      lspin_read = .TRUE.
    END SUBROUTINE readschema_spin 
    !
#endif
    !--------------------------------------------------------------------------
    SUBROUTINE read_magnetization( ierr )
      !------------------------------------------------------------------------
      !
      USE klist,            ONLY : two_fermi_energies, nelup, neldw
      USE ener,             ONLY : ef_up, ef_dw
      USE lsda_mod,         ONLY : starting_magnetization
      USE noncollin_module, ONLY : angle1, angle2, i_cons, mcons, bfield, &
                                   lambda
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      !
      LOGICAL :: found
      INTEGER :: i, nsp
      !
      ierr = 0
      IF ( lstarting_mag_read ) RETURN
      !
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_magnetization(STARTING_MAGNETIZATION=starting_magnetization, &
              ANGLE1=angle1, ANGLE2=angle2, &
              TWO_FERMI_ENERGIES=two_fermi_energies, I_CONS=i_cons, MCONS=mcons, &
              BFIELD=bfield, EF_UP=ef_up, EF_DW=ef_dw, NELUP=nelup, NELDW=neldw, &
              LAMBDA=lambda, FOUND=found, IERR= ierr)
         !
         angle1(:)=angle1(:)*PI/180.d0
         angle2(:)=angle2(:)*PI/180.d0
         !
         IF (two_fermi_energies) THEN
            !
            ef_up = ef_up * e2 
            ef_dw = ef_dw * e2
            !
         ENDIF
         !
      END IF
      !
      CALL mp_bcast( found,  ionode_id, intra_image_comm )
      !
      IF( found ) THEN
         !
         CALL mp_bcast( starting_magnetization,  ionode_id, intra_image_comm )
         CALL mp_bcast( angle1,                  ionode_id, intra_image_comm )
         CALL mp_bcast( angle2,                  ionode_id, intra_image_comm )
         CALL mp_bcast( two_fermi_energies,      ionode_id, intra_image_comm )
         CALL mp_bcast( i_cons,                  ionode_id, intra_image_comm )
         CALL mp_bcast( mcons,                   ionode_id, intra_image_comm )
         CALL mp_bcast( bfield,                  ionode_id, intra_image_comm )
         CALL mp_bcast( nelup,                   ionode_id, intra_image_comm )
         CALL mp_bcast( neldw,                   ionode_id, intra_image_comm )
         CALL mp_bcast( ef_up,                   ionode_id, intra_image_comm )
         CALL mp_bcast( ef_dw,                   ionode_id, intra_image_comm )
         CALL mp_bcast( lambda,                  ionode_id, intra_image_comm )
         !
      ENDIF
      !
      lstarting_mag_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_magnetization
    !
#ifdef __XSD
    !------------------------------------------------------------------------
    SUBROUTINE readschema_magnetization( band_structure_obj, atomic_specs_obj, input_obj) 
      !----------------------------------------------------------------------
      ! 
      USE klist,            ONLY : two_fermi_energies, nelup, neldw
      USE ener,             ONLY : ef_up, ef_dw
      USE lsda_mod,         ONLY : starting_magnetization
      USE noncollin_module, ONLY : angle1, angle2, i_cons, mcons, bfield, &
                                   lambda
      USE electrons_base,   ONLY : set_nelup_neldw
      USE qes_types_module, ONLY : band_structure_type, atomic_species_type, input_type 
      !
      IMPLICIT NONE 
      ! 
      TYPE ( band_structure_type ),INTENT(IN)    :: band_structure_obj
      TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_specs_obj
      TYPE ( input_type ), INTENT(IN)             :: input_obj
      !  
      REAL(DP)                   :: tot_mag_, nelec_, theta, phi, fixed_magnetization(3) 
      INTEGER                    :: nsp_, isp
      !
      IF ( lstarting_mag_read ) RETURN
      bfield = 0.d0
      nelec_ = band_structure_obj%nelec
      two_fermi_energies = band_structure_obj%two_fermi_energies_ispresent
      IF (two_fermi_energies) THEN 
         ef_up = band_structure_obj%two_fermi_energies(1)
         ef_dw = band_structure_obj%two_fermi_energies(2)
         IF ( input_obj%bands%tot_magnetization_ispresent )  THEN 
            tot_mag_ = input_obj%bands%tot_magnetization
         ELSE
            tot_mag_ = 0.d0
         END IF 
         CALL set_nelup_neldw( tot_mag_,  nelec_, nelup, neldw ) 
      END IF 
      nsp_ = atomic_specs_obj%ntyp
      !
      i_cons = 0
      IF (input_obj%spin_constraints_ispresent ) THEN
         lambda = input_obj%spin_constraints%lagrange_multiplier 
         SELECT CASE ( TRIM (input_obj%spin_constraints%spin_constraints ) )
            CASE ( 'atomic') 
               i_cons = 1 
            CASE ( 'atomic_direction' )
               i_cons =  2 
            CASE ( 'total' )
               i_cons = 3 
            CASE ( 'total_direction' ) 
               i_cons = 6 
         END SELECT  
         IF ( input_obj%spin_constraints%target_magnetization_ispresent ) THEN
            fixed_magnetization = input_obj%spin_constraints%target_magnetization 
            SELECT CASE ( i_cons) 
               CASE ( 3 ) 
                  mcons(1,1) = fixed_magnetization(1)
                  mcons(2,1) = fixed_magnetization(2)
                  mcons(3,1) = fixed_magnetization(3)
               CASE ( 6) 
                  mcons(3,1) = fixed_magnetization(3)
               CASE DEFAULT
                  CONTINUE
            END SELECT
         END IF
         !
      END IF
      DO isp = 1, nsp_
         IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent) THEN
             starting_magnetization(isp) = atomic_specs_obj%species(isp)%starting_magnetization      
         END IF                                                                                      
         !                                                                                           
         IF ( band_structure_obj%noncolin ) THEN                                                         
            IF (    atomic_specs_obj%species(isp)%spin_teta_ispresent ) THEN 
               theta = atomic_specs_obj%species(isp)%spin_teta 
               angle1(isp) = theta 
            END IF                                                                  
            IF ( atomic_specs_obj%species(isp)%spin_phi_ispresent ) THEN                  
               phi = atomic_specs_obj%species(isp)%spin_phi
               angle2(isp) = phi
            END IF                                                                     
               !                                                                                     
            IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent .AND. &
                                                                              i_cons == 1 ) THEN 
                !            
                mcons(1,isp) = starting_magnetization(isp) * sin( theta ) * cos( phi )
                mcons(2,isp) = starting_magnetization(isp) * sin( theta ) * sin( phi )
                mcons(3,isp) = starting_magnetization(isp) * cos( theta )
            ELSE IF ( i_cons == 2) THEN  
                mcons(3,isp) = cos(theta) 
            END IF
         ELSE IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent .AND. &
                                                                                  i_cons == 1 ) THEN 
            mcons(1,isp) = starting_magnetization(isp)                                               
         END IF                                                                                      
      END DO   
      !
      lstarting_mag_read = .TRUE.
    END SUBROUTINE readschema_magnetization
    ! 
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_xc( ierr )
      !------------------------------------------------------------------------
      !
      USE ions_base, ONLY : nsp
      USE funct,     ONLY : enforce_input_dft
      USE ldaU,      ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_lmax, &
                            Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_alpha, &
                            Hubbard_J0, Hubbard_beta, U_projection
      USE kernel_table, ONLY : vdw_table_name
      USE acfdt_ener,   ONLY : acfdt_in_pw
      USE control_flags,ONLY : llondon, lxdm, ts_vdw
      USE london_module,ONLY : scal6, lon_rcut
      USE tsvdw_module, ONLY : vdw_isolated
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=20) :: dft_name
      INTEGER           :: nsp_, inlc
      LOGICAL           :: nomsg = .true.
      !
      ierr = 0
      IF ( lxc_read ) RETURN
      !
      IF ( .NOT. lions_read ) &
         CALL errore( 'read_xc', 'read ions first', 1 )
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_xc( dft_name, lda_plus_u, lda_plus_u_kind, U_projection,&
                             Hubbard_lmax, Hubbard_l, nsp_, Hubbard_U, Hubbard_J, &
                             Hubbard_J0, Hubbard_alpha, Hubbard_beta, &
                             inlc, vdw_table_name,  acfdt_in_pw, llondon, scal6, &
                             lon_rcut, lxdm, ts_vdw, vdw_isolated, ierr )
         !
      END IF
      !
      CALL mp_bcast( dft_name,   ionode_id, intra_image_comm )
      CALL mp_bcast( lda_plus_u, ionode_id, intra_image_comm )
      CALL mp_bcast( inlc, ionode_id, intra_image_comm )
      CALL mp_bcast( llondon,    ionode_id, intra_image_comm )
      CALL mp_bcast( lxdm,       ionode_id, intra_image_comm )
      CALL mp_bcast( ts_vdw,     ionode_id, intra_image_comm )
      !
      IF ( lda_plus_u ) THEN
         !
         CALL mp_bcast( lda_plus_u_kind, ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_lmax,  ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_l ,    ionode_id, intra_image_comm )
         CALL mp_bcast( U_projection,  ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_U,     ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_J,     ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_J0,    ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_alpha, ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_beta,  ionode_id, intra_image_comm )
         !
      END IF

      IF ( inlc == 1 .OR. inlc == 2 ) THEN
         CALL mp_bcast( vdw_table_name,  ionode_id, intra_image_comm )
      END IF
      !
      IF ( llondon ) THEN
         CALL mp_bcast( scal6, ionode_id, intra_image_comm )
         CALL mp_bcast( lon_rcut, ionode_id, intra_image_comm )
      END IF
      !
      IF ( ts_vdw ) THEN
         CALL mp_bcast( vdw_isolated, ionode_id, intra_image_comm )
      END IF
      !
      ! SCF EXX/RPA
      !
      CALL mp_bcast( acfdt_in_pw, ionode_id, intra_image_comm )
      !
      IF (acfdt_in_pw) dft_name = 'NOX-NOC'

      ! discard any further attempt to set a different dft
      CALL enforce_input_dft( dft_name, nomsg )
      !
      lxc_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_xc
    !
#ifdef __XSD
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_xc ( atomic_specs, dft_obj ) 
    !-----------------------------------------------------------------------
      ! 
      USE ions_base, ONLY : nsp
      USE funct,     ONLY : enforce_input_dft
      USE ldaU,      ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_lmax, &
                            Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_alpha, &
                            Hubbard_J0, Hubbard_beta, U_projection
      USE kernel_table,     ONLY : vdw_table_name
      USE control_flags,    ONLY : llondon, lxdm, ts_vdw
      USE london_module,    ONLY : scal6, lon_rcut
      USE tsvdw_module,     ONLY : vdw_isolated
      USE qes_types_module, ONLY : atomic_species_type, dft_type
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( atomic_species_type )    :: atomic_specs
      TYPE ( dft_type            )    :: dft_obj
      INTEGER                         :: ihub, nsp_, inlc, isp
      ! 
      CHARACTER(LEN = 256 )           :: label
      CHARACTER(LEN = 20  )           :: dft_name
      CHARACTER(LEN =  3 )            :: symbol
      ! 
      IF ( lxc_read ) RETURN
      dft_name = TRIM(dft_obj%functional) 
      CALL enforce_input_dft ( dft_name, .TRUE. ) 
      lda_plus_u = dft_obj%dftU_ispresent 
      IF ( lda_plus_u ) THEN 
         lda_plus_u_kind = dft_obj%dftU%lda_plus_u_kind
         U_projection = TRIM ( dft_obj%dftU%U_projection_type )
         Hubbard_l =-1 
         IF ( dft_obj%dftU%Hubbard_U_ispresent) THEN 
            loop_on_hubbardU:DO ihub =1, dft_obj%dftU%ndim_Hubbard_U
               symbol = TRIM(dft_obj%dftU%Hubbard_U(ihub)%specie)
               label  = TRIM(dft_obj%dftU%Hubbard_U(ihub)%label ) 
               loop_on_speciesU:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_U(isp) = dft_obj%dftU%Hubbard_U(ihub)%HubbardCommon
                     SELECT CASE ( TRIM (label))
                        CASE ( '1s', '2s', '3s', '4s', '5s', '6s', '7s' ) 
                            Hubbard_l(isp) = 0 
                        CASE ( '2p', '3p', '4p', '5p', '6p' ) 
                            Hubbard_l(isp) = 1 
                        CASE ( '3d', '4d', '5d' ) 
                            Hubbard_l( isp ) = 2 
                        CASE ( '4f', '5f' )  
                            Hubbard_l(isp ) = 3
                        CASE  default 
                            CALL errore ("pw_readschema:", "unrecognized label for Hubbard "//label,&
                                          1 ) 
                     END SELECT   
                     EXIT loop_on_speciesU
                  END IF 
                END DO loop_on_speciesU
            END DO loop_on_hubbardU
         END IF 
         
         IF ( dft_obj%dftU%Hubbard_J0_ispresent ) THEN 
           loop_on_hubbardj0:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J0
               symbol = TRIM(dft_obj%dftU%Hubbard_J0(ihub)%specie)
               loop_on_speciesj0:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN
                     Hubbard_J0(isp) = dft_obj%dftU%Hubbard_J0(ihub)%HubbardCommon
                     EXIT loop_on_speciesj0
                  END IF
               END DO loop_on_speciesj0
            END DO loop_on_hubbardj0
         END IF
         IF ( dft_obj%dftU%Hubbard_alpha_ispresent) THEN 
            loop_on_hubbardAlpha:DO ihub =1, dft_obj%dftU%ndim_Hubbard_alpha
               symbol = TRIM(dft_obj%dftU%Hubbard_alpha(ihub)%specie)
               loop_on_speciesAlpha:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_alpha(isp) = dft_obj%dftU%Hubbard_alpha(ihub)%HubbardCommon
                     EXIT loop_on_speciesAlpha
                  END IF 
                END DO loop_on_speciesAlpha
            END DO loop_on_hubbardAlpha
         END IF 
         IF ( dft_obj%dftU%Hubbard_beta_ispresent) THEN 
            loop_on_hubbardBeta:DO ihub =1, dft_obj%dftU%ndim_Hubbard_beta
               symbol = TRIM(dft_obj%dftU%Hubbard_beta(ihub)%specie)
               loop_on_speciesBeta:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_beta(isp) = dft_obj%dftU%Hubbard_beta(ihub)%HubbardCommon
                     EXIT loop_on_speciesBeta
                  END IF 
                END DO loop_on_speciesBeta
            END DO loop_on_hubbardBeta
         END IF 
         IF ( dft_obj%dftU%Hubbard_J_ispresent) THEN 
            loop_on_hubbardJ:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J
               symbol = TRIM(dft_obj%dftU%Hubbard_J(ihub)%specie)
               loop_on_speciesJ:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_J(:,isp) = dft_obj%dftU%Hubbard_J(ihub)%HubbardJ
                     EXIT loop_on_speciesJ
                  END IF 
                END DO loop_on_speciesJ
            END DO loop_on_hubbardJ
         END IF 
         Hubbard_lmax = MAXVAL( Hubbard_l(1:nsp_) )
      END IF 
      ! 
      IF ( dft_obj%vdW_ispresent ) THEN 
           SELECT CASE( TRIM( dft_obj%vdW%vdw_corr ) )
    !
              CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
                !
                llondon= .TRUE.
                ts_vdw= .FALSE.
                lxdm   = .FALSE.
                !
              CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
                !
                llondon= .FALSE.
                ts_vdw= .TRUE.
                lxdm   = .FALSE.
                !
              CASE( 'XDM', 'xdm' )
                 !
                llondon= .FALSE.
                ts_vdw= .FALSE.
                lxdm   = .TRUE.
                !
              CASE DEFAULT
                !
                llondon= .FALSE.
                ts_vdw = .FALSE.
                lxdm   = .FALSE.
                !
           END SELECT
           SELECT CASE ( TRIM (dft_obj%vdW%non_local_term))
             CASE ('vdw1')  
                inlc = 1
             CASE ('vdw2') 
                inlc = 2
             CASE ('vv10' ) 
                inlc = 3 
             CASE ( 'vdW-DF-x') 
                inlc = 4
             CASE ( 'vdW-DF-y')
                inlc = 5
             CASE ( 'vdW-DF-z')
                inlc = 6
             CASE default 
                inlc = 0 
          END SELECT
          IF (inlc == 0 ) THEN 
             vdw_table_name = ' '
          ELSE IF ( inlc == 3 ) THEN 
             vdw_table_name = 'rVV10_kernel_table'
          ELSE
             vdw_table_name = 'vdW_kernel_table'
          END IF 
          IF (dft_obj%vdW%london_s6_ispresent ) THEN 
             scal6 = dft_obj%vdW%london_s6
          END IF 
          IF ( dft_obj%vdW%london_rcut_ispresent ) THEN 
             lon_rcut = dft_obj%vdW%london_rcut
          END IF 
          IF (dft_obj%vdW%ts_vdW_isolated_ispresent ) THEN 
             vdW_isolated = dft_obj%vdW%ts_vdW_isolated
          END IF 
      END IF 
      !         
      lxc_read = .TRUE.
    END SUBROUTINE readschema_xc
    !  
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_brillouin_zone( ierr )
      !------------------------------------------------------------------------
      !
      USE lsda_mod, ONLY : lsda
      USE klist,    ONLY : nkstot, xk, wk, qnorm
      USE start_k,    ONLY : nks_start, xk_start, wk_start, &
                              nk1, nk2, nk3, k1, k2, k3
      USE symm_base,   ONLY : nrot, s, sname
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      CHARACTER(iotk_attlenx)  :: attr
      !
      INTEGER :: i, ik, num_k_points
      LOGICAL :: found
      !
      ierr = 0
      IF ( lbz_read ) RETURN
      !
      !
      IF ( ionode ) THEN
         !
         ! xk_start and wk_start are ALLOCATABLE inside the function
         CALL qexml_read_bz( NUM_K_POINTS=num_k_points, XK=xk, WK=wk, K1=k1, K2=k2, K3=k3, &
              NK1=nk1, NK2=nk2, NK3=nk3, &
              NKS_START=nks_start, XK_START=xk_start, WK_START=wk_start, QNORM=qnorm, IERR=ierr )
         !
         nkstot = num_k_points
         !
         IF ( lsda ) nkstot = num_k_points * 2
         !
         DO ik = 1, num_k_points
            !
            IF ( lsda ) THEN
               !
               xk(:,ik+num_k_points) = xk(:,ik)
               !
               wk(ik+num_k_points) = wk(ik)
               !
            END IF
            !
         END DO
         !
      END IF
      !
      CALL mp_bcast( nkstot, ionode_id, intra_image_comm )
      CALL mp_bcast( xk,     ionode_id, intra_image_comm )
      CALL mp_bcast( wk,     ionode_id, intra_image_comm )
      CALL mp_bcast( nk1, ionode_id, intra_image_comm )
      CALL mp_bcast( nk2, ionode_id, intra_image_comm )
      CALL mp_bcast( nk3, ionode_id, intra_image_comm )
      CALL mp_bcast( k1, ionode_id, intra_image_comm )
      CALL mp_bcast( k2, ionode_id, intra_image_comm )
      CALL mp_bcast( k3, ionode_id, intra_image_comm )
      CALL mp_bcast( qnorm, ionode_id, intra_image_comm)

      CALL mp_bcast( nks_start, ionode_id, intra_image_comm )
      IF (nks_start>0.and..NOT.ionode) THEN
         IF (.NOT.ALLOCATED(xk_start)) ALLOCATE(xk_start(3,nks_start))
         IF (.NOT.ALLOCATED(wk_start)) ALLOCATE(wk_start(nks_start))
      ENDIF
      IF (nks_start>0) THEN
         CALL mp_bcast( xk_start, ionode_id, intra_image_comm )
         CALL mp_bcast( wk_start, ionode_id, intra_image_comm )
      ENDIF
      CALL mp_bcast(  nrot, ionode_id, intra_image_comm )
      CALL mp_bcast(     s, ionode_id, intra_image_comm )
      CALL mp_bcast( sname, ionode_id, intra_image_comm )
      !
      lbz_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_brillouin_zone
    !
#ifdef __XSD
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE readschema_brillouin_zone( k_pointIBZ_obj , occupations_obj, symmetries_obj, band_struct_obj )
    !-----------------------------------------------------------------------------------------------------
       !
       USE lsda_mod, ONLY : lsda
       USE klist,    ONLY : nkstot, xk, wk, qnorm
       USE start_k,  ONLY : nks_start, xk_start, wk_start, &
                              nk1, nk2, nk3, k1, k2, k3
       USE symm_base,ONLY : nrot, s, sname
       USE qes_types_module, ONLY : k_points_IBZ_type, occupations_type, symmetries_type, band_structure_type
       !
       IMPLICIT NONE
       !
       TYPE ( k_points_IBZ_type ),  INTENT(IN)    :: k_pointIBZ_obj
       TYPE ( occupations_type ),   INTENT(IN)    :: occupations_obj
       TYPE ( symmetries_type ),    INTENT(IN)    :: symmetries_obj 
       TYPE ( band_structure_type ),INTENT(IN)    :: band_struct_obj 
       INTEGER                                    :: ik, isym, nks_
       ! 
       IF ( lbz_read ) RETURN 
       nks_ = band_struct_obj%nks
       nkstot = nks_
       IF ( band_struct_obj%lsda ) nkstot = nkstot * 2  
       ! 
       ! 
       DO ik = 1, nks_
          xk(:,ik) = band_struct_obj%ks_energies(ik)%k_point%k_point(:) 
       END DO 
       !   
       IF ( k_pointIBZ_obj%monkhorst_pack_ispresent ) THEN 
          nks_start = 0 
          nk1 = k_pointIBZ_obj%monkhorst_pack%nk1 
          nk2 = k_pointIBZ_obj%monkhorst_pack%nk2
          nk3 = k_pointIBZ_obj%monkhorst_pack%nk3 
           k1 = k_pointIBZ_obj%monkhorst_pack%k1
           k2 = k_pointIBZ_obj%monkhorst_pack%k2
           k3 = k_pointIBZ_obj%monkhorst_pack%k3
       ELSE IF (k_pointIBZ_obj%nk_ispresent .AND. &
                k_pointIBZ_obj%k_point_ispresent ) THEN 
           nks_start = k_pointIBZ_obj%nk
           ALLOCATE (xk_start(3,nks_start), wk_start(nks_start))
           DO ik =1, nks_start
               xk_start(:,ik) = k_pointIBZ_obj%k_point(ik)%k_point(:) 
               IF ( k_pointIBZ_obj%k_point(ik)%weight_ispresent) THEN 
                  wk_start(ik) = k_pointIBZ_obj%k_point(ik)%weight 
               ELSE 
                  wk_start(ik) = 0.d0
               END IF 
           END DO
       ELSE 
           CALL errore ("pw_readschema: ", &
                        " no information found for initializing brillouin zone information", 1)
       END IF  
       ! 
       IF ( .NOT. lsymm_read  ) THEN 
          nrot = symmetries_obj%nrot
          DO isym =1, symmetries_obj%ndim_symmetry
             s(:,:,isym)     = symmetries_obj%symmetry(isym)%rotation%mat
             sname(isym) = TRIM ( symmetries_obj%symmetry(isym)%info%name) 
         END DO 
       END IF    
       lbz_read = .TRUE.
    END SUBROUTINE readschema_brillouin_zone     
#endif            
    !------------------------------------------------------------------------
    SUBROUTINE read_occupations( ierr )
      !------------------------------------------------------------------------
      !
      USE lsda_mod,       ONLY : lsda, nspin
      USE fixed_occ,      ONLY : tfixed_occ, f_inp
      USE ktetra,         ONLY : ntetra, tetra, ltetra
      USE klist,          ONLY : lgauss, ngauss, degauss, smearing
      USE electrons_base, ONLY : nupdwn 
      USE wvfct,          ONLY : nbnd
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      CHARACTER(iotk_attlenx)  :: attr
      !
      INTEGER :: i
      LOGICAL :: found
      !
      ierr = 0
      IF ( locc_read ) RETURN
      !
      IF ( ionode ) THEN
         !
         ! necessary to don't send nbnd and nspin as input in read_occ
         IF ( .NOT. ALLOCATED( f_inp ) ) THEN
            !
            IF ( nspin == 4 ) THEN
               ALLOCATE( f_inp( nbnd, 1 ) )
            ELSE
               ALLOCATE( f_inp( nbnd, nspin ) )
            ENDIF
            !
         ENDIF
         !
         f_inp( :, :) = 0.0d0
         !
         CALL qexml_read_occ( LGAUSS=lgauss, NGAUSS=ngauss, DEGAUSS=degauss, &
                               LTETRA=ltetra, NTETRA=ntetra, TETRA=tetra, TFIXED_OCC=tfixed_occ, &
                               NSTATES_UP=nupdwn(1), NSTATES_DW=nupdwn(2), INPUT_OCC=f_inp, IERR=ierr )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         IF ( lgauss ) THEN
            !
            SELECT CASE (ngauss )
            CASE (0)
               smearing = 'gaussian'
            CASE (1)
               smearing = 'Methfessel-Paxton'
            CASE (-1)
               smearing = 'Marzari-Vanderbilt'
            CASE (-99)
               smearing = 'Fermi-Dirac'
            CASE DEFAULT
               CALL errore('read_occupations',&
                    'wrong smearing index', abs(1000+ngauss) )
            END SELECT
            !
            degauss = degauss * e2
            !
         ELSE
            !
            ngauss = 0
            degauss = 0.d0
            !
         END IF
         !
         IF ( .NOT. ltetra ) THEN
            !
            ntetra = 0
            !
         END IF
         !
         IF ( .NOT. tfixed_occ ) THEN
            !
            DEALLOCATE( f_inp )
            !
         ENDIF
         !
         !
      END IF
      !
      CALL mp_bcast( lgauss, ionode_id, intra_image_comm )
      !
      IF ( lgauss ) THEN
         !
         CALL mp_bcast( ngauss,  ionode_id, intra_image_comm )
         CALL mp_bcast( degauss, ionode_id, intra_image_comm )
         CALL mp_bcast( smearing, ionode_id, intra_image_comm )
         !
      END IF
      !
      CALL mp_bcast( ltetra, ionode_id, intra_image_comm )
      !
      IF ( ltetra ) THEN
         !
         CALL mp_bcast( ntetra, ionode_id, intra_image_comm )
         CALL mp_bcast( tetra,  ionode_id, intra_image_comm )
         !
      END IF
      !
      CALL mp_bcast( tfixed_occ, ionode_id, intra_image_comm )
      !
      IF ( tfixed_occ ) THEN
         !
         CALL mp_bcast( nupdwn, ionode_id, intra_image_comm )
         !
         IF ( .NOT. ALLOCATED( f_inp ) ) THEN
            !
            IF ( nspin == 4 ) THEN
               ALLOCATE( f_inp( nbnd, 1 ) )
            ELSE
               ALLOCATE( f_inp( nbnd, nspin ) )
            END IF
            !
         ENDIF
         !
         CALL mp_bcast( f_inp, ionode_id, intra_image_comm )
         !
      ENDIF
      !
      locc_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_occupations
    !
#ifdef __XSD
    !--------------------------------------------------------------------------------------------------
    SUBROUTINE readschema_occupations( input_obj, band_struct_obj ) 
      !------------------------------------------------------------------------------------------------
      ! 
      USE lsda_mod,         ONLY : lsda, nspin
      USE fixed_occ,        ONLY : tfixed_occ, f_inp
      USE ktetra,           ONLY : ntetra, ltetra
      USE klist,            ONLY : lgauss, ngauss, degauss, smearing
      USE electrons_base,   ONLY : nupdwn 
      USE wvfct,            ONLY : nbnd
      USE input_parameters, ONLY : input_parameters_occupations => occupations
      USE qes_types_module, ONLY : input_type, band_structure_type
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( input_type ),INTENT(IN)              :: input_obj
      TYPE ( band_structure_type ),INTENT(IN)     :: band_struct_obj 
      INTEGER                                     :: ispin, nk1, nk2, nk3, aux_dim1, aux_dim2 
      ! 
      IF ( locc_read ) RETURN 
      lsda= band_struct_obj%lsda
      nbnd = band_struct_obj%nbnd
      IF ( band_struct_obj%nbnd_up_ispresent ) nupdwn(1) = band_struct_obj%nbnd_up
      IF ( band_struct_obj%nbnd_dw_ispresent ) nupdwn(2) = band_struct_obj%nbnd_dw 
      IF ( lsda )  THEN 
         nspin = 2  
      ELSE IF ( band_struct_obj%noncolin) THEN 
         nspin = 4 
      ELSE 
         nspin = 1 
      END IF 
      !
      lgauss = .FALSE. 
      ltetra = .FALSE. 
      ngauss = 0
      input_parameters_occupations = TRIM ( input_obj%bands%occupations%occupations ) 
      IF (TRIM(input_obj%bands%occupations%occupations) == 'tetrahedra' ) THEN 
        ltetra = .TRUE. 
        nk1 = input_obj%k_points_IBZ%monkhorst_pack%nk1
        nk2 = input_obj%k_points_IBZ%monkhorst_pack%nk2
        nk3 = input_obj%k_points_IBZ%monkhorst_pack%nk3
        ntetra = 6* nk1 * nk2 * nk3 
      ELSE IF ( TRIM (input_obj%bands%occupations%occupations) == 'smearing') THEN 
        lgauss = .TRUE.  
        degauss = input_obj%bands%smearing%degauss
        SELECT CASE ( TRIM( input_obj%bands%smearing%smearing ) )
           CASE ( 'gaussian', 'gauss', 'Gaussian', 'Gauss' )
             ngauss = 0
             smearing  = 'gaussian'
           CASE ( 'methfessel-paxton', 'm-p', 'mp', 'Methfessel-Paxton', 'M-P', 'MP' )
             ngauss = 1
             smearing = 'Methfessel-Paxton'
           CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'Marzari-Vanderbilt', 'M-V', 'MV')
             ngauss = -1
             smearing  = 'Marzari-Vanderbilt'
           CASE ( 'fermi-dirac', 'f-d', 'fd', 'Fermi-Dirac', 'F-D', 'FD')
             ngauss = -99
             smearing = 'Fermi-Dirac'
        END SELECT
      ELSE IF ( TRIM (input_obj%bands%occupations%occupations) == 'from_input' .AND. & 
                input_obj%bands%inputOccupations_ispresent ) THEN 
           tfixed_occ = .TRUE.
           IF ( .NOT. ALLOCATED(f_inp))  &
              aux_dim2 = input_obj%bands%ndim_inputOccupations
              aux_dim1 = MAXVAL(input_obj%bands%inputOccupations(1:aux_dim2)%ndim_vec)
              ALLOCATE (f_inp( aux_dim1, aux_dim2))
           DO ispin = 1, input_obj%bands%ndim_inputOccupations
              f_inp(:,ispin) = input_obj%bands%inputOccupations(ispin)%vec
           END DO
      END IF       
     !
     locc_read = .TRUE.
    END SUBROUTINE readschema_occupations
 !
 !---------------------------------------------------------------------------
#endif 
        
 
    !------------------------------------------------------------------------
    SUBROUTINE read_band_structure( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE control_flags, ONLY : lkpoint_dir
      USE basis,    ONLY : natomwfc
      USE lsda_mod, ONLY : lsda, isk
      USE klist,    ONLY : nkstot, wk, nelec
      USE wvfct,    ONLY : et, wg, nbnd
      USE ener,     ONLY : ef, ef_up, ef_dw
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: ik, ik_eff, num_k_points
      LOGICAL :: found, two_fermi_energies_
      CHARACTER(LEN=256) :: filename
      !
      ierr = 0
      IF ( lbs_read ) RETURN
      !
      IF ( .NOT. lspin_read ) &
         CALL errore( 'read_band_structure', 'read spin first', 1 )
      IF ( .NOT. lbz_read ) &
         CALL errore( 'read_band_structure', 'read band_structure first', 1 )
      !
      !
      IF ( ionode ) THEN
         ! we don't need to read nspin, noncolin
         CALL qexml_read_bands_info( NBND=nbnd, NUM_K_POINTS=num_k_points, NATOMWFC=natomwfc, &
                                     NELEC=nelec, EF=ef, TWO_FERMI_ENERGIES=two_fermi_energies_, &
                                     EF_UP=ef_up, EF_DW=ef_dw, IERR=ierr )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         IF ( .NOT. two_fermi_energies_) THEN
            ef = ef * e2
         ELSE
            ef = 0.d0
            ef_up = ef_up * e2
            ef_dw = ef_dw * e2
         END IF
         !
      END IF
      !
      num_k_points = nkstot
      !
      IF ( lsda ) num_k_points = nkstot / 2
      !
      IF ( ionode ) THEN
         !
         IF (.NOT.lkpoint_dir) filename = TRIM( dirname ) // '/' // TRIM( xmlpun )//'.eig'
         !
         CALL qexml_read_bands_pw( num_k_points, nbnd, nkstot, lsda, lkpoint_dir, filename , ISK=isk, ET=et, WG=wg , IERR=ierr)
         !
         et(:,:) = et(:,:) * e2
         !
         FORALL( ik = 1:nkstot ) wg(:,ik) = wg(:,ik)*wk(ik)
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      CALL mp_bcast( nelec,    ionode_id, intra_image_comm )
      CALL mp_bcast( natomwfc, ionode_id, intra_image_comm )
      CALL mp_bcast( nbnd,     ionode_id, intra_image_comm )
      CALL mp_bcast( isk,      ionode_id, intra_image_comm )
      CALL mp_bcast( et,       ionode_id, intra_image_comm )
      CALL mp_bcast( wg,       ionode_id, intra_image_comm )
      CALL mp_bcast( ef,       ionode_id, intra_image_comm )
      !
      lbs_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_band_structure
    !
#ifdef __XSD
    !------------------------------------------------------------------------
    SUBROUTINE readschema_band_structure( band_struct_obj )
      !------------------------------------------------------------------------
      !
      USE control_flags, ONLY : lkpoint_dir
      USE constants,     ONLY : e2
      USE basis,    ONLY : natomwfc
      USE lsda_mod, ONLY : lsda, isk
      USE klist,    ONLY : nkstot, wk, nelec
      USE wvfct,    ONLY : et, wg, nbnd
      USE ener,     ONLY : ef, ef_up, ef_dw
      USE qes_types_module, ONLY : band_structure_type
      !
      IMPLICIT NONE
      TYPE ( band_structure_type)         :: band_struct_obj
      INTEGER                             :: ik, nbnd_, nbnd_up_, nbnd_dw_
      ! 
      lkpoint_dir = .TRUE.  ! TO BE DISCUSSED 
      lsda = band_struct_obj%lsda
      nkstot = band_struct_obj%nks 
      IF ( lsda) THEN 
        nkstot = nkstot * 2 
        isk(1:nkstot/2) = 1
        isk(nkstot/2+1:nkstot) = 2 
      ELSE 
        isk(1:nkstot)   = 1 
      END IF 
      ! 
      nelec = band_struct_obj%nelec
      nbnd  = band_struct_obj%nbnd 
      natomwfc = band_struct_obj%num_of_atomic_wfc
      IF ( band_struct_obj%fermi_energy_ispresent) THEN 
         ef = band_struct_obj%fermi_energy*e2 
      ELSE IF ( band_struct_obj%two_fermi_energies_ispresent ) THEN 
         ef = 0.d0 
         ef_up = band_struct_obj%two_fermi_energies(1)*e2
         ef_dw = band_struct_obj%two_fermi_energies(2)*e2
      ELSE 
         ef = 0.d0
         ef_up = 0.d0
         ef_dw = 0.d0
      END IF 
      DO ik =1, band_struct_obj%ndim_ks_energies
         IF ( band_struct_obj%lsda) THEN
            IF ( band_struct_obj%nbnd_up_ispresent .AND. band_struct_obj%nbnd_dw_ispresent) THEN
               nbnd_up_ = band_struct_obj%nbnd_up
               nbnd_dw_ = band_struct_obj%nbnd_dw 
            ELSE IF ( band_struct_obj%nbnd_up_ispresent ) THEN 
               nbnd_up_ = band_struct_obj%nbnd_up
               nbnd_dw_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues - nbnd_up_
            ELSE IF ( band_struct_obj%nbnd_dw_ispresent ) THEN 
               nbnd_dw_ = band_struct_obj%nbnd_dw
               nbnd_up_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues - nbnd_dw_ 
            ELSE 
               nbnd_up_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues/2  
               nbnd_dw_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues/2
            END IF
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            wk( ik + band_struct_obj%ndim_ks_energies ) = wk(ik) 
            et(1:nbnd_up_,ik) = band_struct_obj%ks_energies(ik)%eigenvalues(1:nbnd_up_)*e2
            et(1:nbnd_dw_,ik+band_struct_obj%ndim_ks_energies) =  &
                             band_struct_obj%ks_energies(ik)%eigenvalues(nbnd_up_+1:nbnd_up_+nbnd_dw_)*e2
            wg(1:nbnd_up_,ik) = band_struct_obj%ks_energies(ik)%occupations(1:nbnd_up_)*wk(ik)
            wg(1:nbnd_dw_,ik+band_struct_obj%ndim_ks_energies) =  &
                             band_struct_obj%ks_energies(ik)%occupations(nbnd_up_+1:nbnd_up_+nbnd_dw_)*wk(ik)
         ELSE 
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            nbnd_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues
            et (1:nbnd_,ik) = band_struct_obj%ks_energies(ik)%eigenvalues(1:nbnd_)*e2
            wg (1:nbnd_,ik) = band_struct_obj%ks_energies(ik)%occupations(1:nbnd_)*wk(ik)
         END IF  
      END DO 
    END SUBROUTINE readschema_band_structure 
    ! 
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_wavefunctions( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... This routines reads wavefunctions from the new file format and
      ! ... writes them into the old format
      !
      USE control_flags,        ONLY : twfcollect, lkpoint_dir
      USE cell_base,            ONLY : tpiba2
      USE lsda_mod,             ONLY : nspin, isk
      USE klist,                ONLY : nkstot, wk, nks, xk, ngk
      USE wvfct,                ONLY : npw, npwx, et, wg, nbnd
      USE gvecw,                ONLY : ecutwfc
      USE wavefunctions_module, ONLY : evc
      USE io_files,             ONLY : nwordwfc, iunwfc
      USE buffers,              ONLY : save_buffer
      USE gvect,                ONLY : ngm, ngm_g, g, ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol
      USE mp_images,            ONLY : nproc_image, intra_image_comm
      USE mp_pools,             ONLY : kunit, nproc_pool, me_pool, root_pool, &
                                       intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : me_bgrp, nbgrp, root_bgrp, &
                                       intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256)   :: filename
      INTEGER              :: ik, ipol, ik_eff, num_k_points
      INTEGER, ALLOCATABLE :: kisort(:)
      INTEGER              :: npool, nkbl, nkl, nkr, npwx_g
      INTEGER              :: nupdwn(2), ike, iks, npw_g, ispin
      INTEGER, ALLOCATABLE :: ngk_g(:)
      INTEGER, ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:)
      LOGICAL              :: opnd
      REAL(DP),ALLOCATABLE :: gk(:)
      REAL(DP)             :: scalef

      !
      ! The ierr output var is actually not given any value
      ! except this initialization
      !
      ierr = 0
      !
      IF ( iunwfc > 0 ) THEN
         !
         INQUIRE( UNIT = iunwfc, OPENED = opnd )
         !
         IF ( .NOT. opnd ) CALL errore( 'read_wavefunctions', &
                    & 'wavefunctions unit (iunwfc) is not opened', 1 )
      END IF
      !
      IF ( nkstot > 0 ) THEN
         !
         ! ... find out the number of pools
         !
         npool = nproc_image / nproc_pool
         !
         ! ... find out number of k points blocks
         !
         nkbl = nkstot / kunit
         !
         !  k points per pool
         !
         nkl = kunit * ( nkbl / npool )
         !
         ! ... find out the reminder
         !
         nkr = ( nkstot - nkl * npool ) / kunit
         !
         ! ... Assign the reminder to the first nkr pools
         !
         IF ( my_pool_id < nkr ) nkl = nkl + kunit
         !
         ! ... find out the index of the first k point in this pool
         !
         iks = nkl * my_pool_id + 1
         !
         IF ( my_pool_id >= nkr ) iks = iks + nkr * kunit
         !
         ! ... find out the index of the last k point in this pool
         !
         ike = iks + nkl - 1
         !
      END IF
      !
      ! ... find out the global number of G vectors: ngm_g  
      !
      ngm_g = ngm
      !
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      !
      ! ... build the igk_l2g array, yielding the correspondence between
      ! ... the local k+G index and the global G index - see also ig_l2g
      !
      ALLOCATE ( igk_l2g( npwx, nks ) )
      igk_l2g = 0
      !
      ALLOCATE( kisort( npwx ), gk(npwx) )
      !
      DO ik = 1, nks
         !
         kisort = 0
         npw    = npwx
         !
         CALL gk_sort( xk(1,ik+iks-1), ngm, g, &
                       ecutwfc/tpiba2, npw, kisort(1), gk )
         !
         CALL gk_l2gmap( ngm, ig_l2g(1), npw, kisort(1), igk_l2g(1,ik) )
         !
         ngk(ik) = npw
         !
      END DO
      !
      DEALLOCATE( gk, kisort )
      !
      ! ... compute the global number of G+k vectors for each k point
      !
      ALLOCATE( ngk_g( nkstot ) )
      !
      ngk_g = 0
      ngk_g(iks:ike) = ngk(1:nks)
      !
      CALL mp_sum( ngk_g, inter_pool_comm )
      CALL mp_sum( ngk_g, intra_pool_comm )
      ngk_g = ngk_g / nbgrp
      !
      ! ... compute the Maximum G vector index among all G+k an processors
      !
      npw_g = MAXVAL( igk_l2g(:,:) )
      !
      CALL mp_max( npw_g, inter_pool_comm )
      CALL mp_max( npw_g, intra_pool_comm )

      !
      ! ... compute the Maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! 
      ! ... define a further l2g map to read gkvectors and wfc coherently 
      ! 
      ALLOCATE( igk_l2g_kdip( npwx_g, nks ) )
      igk_l2g_kdip = 0
      !
      DO ik = iks, ike
         !
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik-iks+1), &
                              igk_l2g(1,ik-iks+1), igk_l2g_kdip(1,ik-iks+1) )
      END DO
      !
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "EIGENVECTORS" )
         !
      END IF
      !
      num_k_points = nkstot
      !
      IF ( nspin == 2 ) num_k_points = nkstot / 2
      !
      k_points_loop: DO ik = 1, num_k_points
         !
         IF ( ionode ) THEN
            !
            CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            IF ( nspin == 2 .OR. noncolin ) THEN
                !
                CALL iotk_scan_begin( iunpun, "WFC.1", FOUND = twfcollect  )
                IF ( twfcollect ) CALL iotk_scan_end( iunpun, "WFC.1" )
                !
            ELSE
                !
                CALL iotk_scan_begin( iunpun, "WFC", FOUND = twfcollect  )
                IF ( twfcollect ) CALL iotk_scan_end( iunpun, "WFC" )
                !
            ENDIF
            !
         END IF
         !
         CALL mp_bcast( twfcollect, ionode_id, intra_image_comm )
         !
         IF ( .NOT. twfcollect ) THEN
            !
            IF ( ionode ) THEN
               !
               CALL iotk_scan_end( iunpun, &
                                   "K-POINT" // TRIM( iotk_index( ik ) ) )
               !
            END IF
            !
            EXIT k_points_loop
            !
         END IF
         !
         IF ( nspin == 2 ) THEN
            !
            ispin = 1 
            evc=(0.0_DP, 0.0_DP)
            !
            ! ... no need to read isk here: they are read from band structure
            ! ... and correctly distributed across pools in read_file
            !!! isk(ik) = 1
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, ispin, &
                                  DIR=lkpoint_dir ) )
               !
            END IF
            !
            CALL read_wfc( iunout, ik, nkstot, kunit, ispin, nspin,      &
                           evc, npw_g, nbnd, igk_l2g_kdip(:,ik-iks+1),   &
                           ngk(ik-iks+1), filename, scalef, &
                           ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
            !
            IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
               !
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
               !
            END IF
            !
            ispin = 2
            ik_eff = ik + num_k_points
            evc=(0.0_DP, 0.0_DP)
            !
            ! ... no need to read isk here (see above why)
            !isk(ik_eff) = 2
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, ispin, &
                                DIR=lkpoint_dir ) )
               !
            END IF
            !
            CALL read_wfc( iunout, ik_eff, nkstot, kunit, ispin, nspin,      &
                           evc, npw_g, nbnd, igk_l2g_kdip(:,ik_eff-iks+1),   &
                           ngk(ik_eff-iks+1), filename, scalef, &
                           ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
            !
            IF ( ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
               !
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik_eff-iks+1) )
               !
            END IF
            !
         ELSE
            !
            ! ... no need to read isk here (see above why)
            !isk(ik) = 1
            !
            evc=(0.0_DP, 0.0_DP)
            IF ( noncolin ) THEN
               !
               DO ipol = 1, npol
                  !
                  IF ( ionode ) THEN
                     !
                     filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, ipol, &
                                         DIR=lkpoint_dir ) )
                     !
                  END IF
                  !
                  !!! TEMP
                  nkl=(ipol-1)*npwx+1
                  nkr= ipol   *npwx
                  CALL read_wfc( iunout, ik, nkstot, kunit, ispin,          &
                                 npol, evc(nkl:nkr,:), npw_g, nbnd,         &
                                 igk_l2g_kdip(:,ik-iks+1), ngk(ik-iks+1),   &
                                 filename, scalef, & 
                                 ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
                  !
               END DO
               !
            ELSE
               !
               IF ( ionode ) THEN
                  !
                  filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, &
                                         DIR=lkpoint_dir ) )
                  !
               END IF
               !
               ! workaround for pot parallelization ( Viet Nguyen / SdG )
               ! -pot parallelization uses mp_image communicators
               ! note that ionode must be also reset in the similar way 
               ! to image parallelization
               CALL read_wfc( iunout, ik, nkstot, kunit, ispin, nspin,         &
                              evc, npw_g, nbnd, igk_l2g_kdip(:,ik-iks+1),      &
                              ngk(ik-iks+1), filename, scalef, &
                              ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
               !
            END IF
            !
            IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
               !
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
               !
               ! the following two line can be used to debug read_wfc
               ! WRITE(200+10*ik+me_pool,fmt="(2D18.10)") evc
               ! CLOSE(200+10*ik+me_pool )
               !
            END IF
            !
         END IF
         !
         IF ( ionode ) THEN
            !
            CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
         END IF
         !
      END DO k_points_loop
      !
      DEALLOCATE ( igk_l2g )
      DEALLOCATE ( igk_l2g_kdip )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( iunpun, "EIGENVECTORS" )
         !
         !CALL iotk_close_read( iunpun )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_wavefunctions
    !
#ifdef __XSD
    !------------------------------------------------------------------------
    SUBROUTINE read_collected_to_evc( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... This routines reads wavefunctions from the new file format and
      ! ... writes them into the old format
      !
      USE control_flags,        ONLY : twfcollect, lkpoint_dir
      USE cell_base,            ONLY : tpiba2
      USE lsda_mod,             ONLY : nspin, isk
      USE klist,                ONLY : nkstot, wk, nks, xk, ngk
      USE wvfct,                ONLY : npw, npwx, g2kin, et, wg, nbnd
      USE gvecw,                ONLY : ecutwfc
      USE wavefunctions_module, ONLY : evc
      USE io_files,             ONLY : nwordwfc, iunwfc
      USE buffers,              ONLY : save_buffer
      USE gvect,                ONLY : ngm, ngm_g, g, ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol
      USE mp_images,            ONLY : nproc_image, intra_image_comm
      USE mp_pools,             ONLY : kunit, nproc_pool, me_pool, root_pool, &
                                       intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : me_bgrp, nbgrp, root_bgrp, &
                                       intra_bgrp_comm
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256)   :: filename
      INTEGER              :: ik, ipol, ik_eff, num_k_points
      INTEGER, ALLOCATABLE :: kisort(:)
      INTEGER              :: npool, nkbl, nkl, nkr, npwx_g
      INTEGER              :: nupdwn(2), ike, iks, npw_g, ispin
      INTEGER, ALLOCATABLE :: ngk_g(:)
      INTEGER, ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:)
      LOGICAL              :: opnd
      REAL(DP)             :: scalef
      REAL(DP),ALLOCATABLE :: gkin_aux(:)

      !
      ! The ierr output var is actually not given any value
      ! except this initialization
      !
      ierr = 0
      !
      IF ( iunwfc > 0 ) THEN
         !
         INQUIRE( UNIT = iunwfc, OPENED = opnd )
         !
         IF ( .NOT. opnd ) CALL errore( 'read_wavefunctions', &
                    & 'wavefunctions unit (iunwfc) is not opened', 1 )
      END IF
      ! 
      IF ( nkstot > 0 ) THEN
         !
         ! ... find out the number of pools
         !
         npool = nproc_image / nproc_pool
         !
         ! ... find out number of k points blocks
         !
         nkbl = nkstot / kunit
         !
         !  k points per pool
         !
         nkl = kunit * ( nkbl / npool )
         !
         ! ... find out the reminder
         !
         nkr = ( nkstot - nkl * npool ) / kunit
         !
         ! ... Assign the reminder to the first nkr pools
         !
         IF ( my_pool_id < nkr ) nkl = nkl + kunit
         !
         ! ... find out the index of the first k point in this pool
         !
         iks = nkl * my_pool_id + 1
         !
         IF ( my_pool_id >= nkr ) iks = iks + nkr * kunit
         !
         ! ... find out the index of the last k point in this pool
         !
         ike = iks + nkl - 1
         !
      END IF
      !
      ! ... find out the global number of G vectors: ngm_g  
      !
      ngm_g = ngm
      !
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      !
      ! ... build the igk_l2g array, yielding the correspondence between
      ! ... the local k+G index and the global G index - see also ig_l2g
      !
      ALLOCATE ( igk_l2g( npwx, nks ) ,gkin_aux(ngm))
      igk_l2g = 0
      !
      ALLOCATE( kisort( npwx ) )
      !
      DO ik = 1, nks
         !
         kisort = 0
         npw    = npwx
         !
         CALL gk_sort( xk(1,ik+iks-1), ngm, g, &
                       ecutwfc/tpiba2, npw, kisort(1), gkin_aux )
         !
         CALL gk_l2gmap( ngm, ig_l2g(1), npw, kisort(1), igk_l2g(1,ik) )
         !
         ngk(ik) = npw
         !
      END DO
      DEALLOCATE (gkin_aux)
      !
      DEALLOCATE( kisort )
      !
      ! ... compute the global number of G+k vectors for each k point
      !
      ALLOCATE( ngk_g( nkstot ) )
      !
      ngk_g = 0
      ngk_g(iks:ike) = ngk(1:nks)
      !
      CALL mp_sum( ngk_g, inter_pool_comm )
      CALL mp_sum( ngk_g, intra_pool_comm )
      ngk_g = ngk_g / nbgrp
      !
      ! ... compute the Maximum G vector index among all G+k an processors
      !
      npw_g = MAXVAL( igk_l2g(:,:) )
      !
      CALL mp_max( npw_g, inter_pool_comm )
      CALL mp_max( npw_g, intra_pool_comm )

      !
      ! ... compute the Maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! 
      ! ... define a further l2g map to read gkvectors and wfc coherently 
      ! 
      ALLOCATE( igk_l2g_kdip( npwx_g, nks ) )
      igk_l2g_kdip = 0
      !
      DO ik = iks, ike
         !
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik-iks+1), &
                              igk_l2g(1,ik-iks+1), igk_l2g_kdip(1,ik-iks+1) )
      END DO
      !
      !
      !
      num_k_points = nkstot
      !
      IF ( nspin == 2 ) num_k_points = nkstot / 2
      !
      IF ( .NOT. twfcollect ) RETURN 
      k_points_loop: DO ik = 1, num_k_points
         !
         IF ( nspin == 2 ) THEN
            !
            ispin = 1 
            evc=(0.0_DP, 0.0_DP)
            !
            ! ... no need to read isk here: they are read from band structure
            ! ... and correctly distributed across pools in read_file
            !!! isk(ik) = 1
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, ispin, &
                                  DIR=lkpoint_dir ) )
               !
            END IF
            !
            CALL read_wfc( iunout, ik, nkstot, kunit, ispin, nspin,      &
                           evc, npw_g, nbnd, igk_l2g_kdip(:,ik-iks+1),   &
                           ngk(ik-iks+1), filename, scalef, &
                           ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
            !
            IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
               !
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
               !
            END IF
            !
            ispin = 2
            ik_eff = ik + num_k_points
            evc=(0.0_DP, 0.0_DP)
            !
            ! ... no need to read isk here (see above why)
            !isk(ik_eff) = 2
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, ispin, &
                                DIR=lkpoint_dir ) )
               !
            END IF
            !
            CALL read_wfc( iunout, ik_eff, nkstot, kunit, ispin, nspin,      &
                           evc, npw_g, nbnd, igk_l2g_kdip(:,ik_eff-iks+1),   &
                           ngk(ik_eff-iks+1), filename, scalef, &
                           ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
            !
            IF ( ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
               !
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik_eff-iks+1) )
               !
            END IF
            !
         ELSE
            !
            ! ... no need to read isk here (see above why)
            !isk(ik) = 1
            !
            evc=(0.0_DP, 0.0_DP)
            IF ( noncolin ) THEN
               !
               DO ipol = 1, npol
                  !
                  IF ( ionode ) THEN
                     !
                     filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, ipol, &
                                         DIR=lkpoint_dir ) )
                     !
                  END IF
                  !
                  !!! TEMP
                  nkl=(ipol-1)*npwx+1
                  nkr= ipol   *npwx
                  CALL read_wfc( iunout, ik, nkstot, kunit, ispin,          &
                                 npol, evc(nkl:nkr,:), npw_g, nbnd,         &
                                 igk_l2g_kdip(:,ik-iks+1), ngk(ik-iks+1),   &
                                 filename, scalef, & 
                                 ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
                  !
               END DO
               !
            ELSE
               !
               IF ( ionode ) THEN
                  !
                  filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, &
                                         DIR=lkpoint_dir ) )
                  !
               END IF
               !
               CALL read_wfc( iunout, ik, nkstot, kunit, ispin, nspin,         &
                              evc, npw_g, nbnd, igk_l2g_kdip(:,ik-iks+1),      &
                              ngk(ik-iks+1), filename, scalef, &
                              ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
               !
            END IF
            !
            IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
            END IF
            !
         END IF
         !
      END DO k_points_loop
      !
      DEALLOCATE ( igk_l2g )
      DEALLOCATE ( igk_l2g_kdip )
      !
      RETURN
      !

    END SUBROUTINE read_collected_to_evc
    !
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_ef( ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the Fermi energy and the number of electrons
      !
      USE ener,  ONLY : ef, ef_up, ef_dw
      USE klist, ONLY : two_fermi_energies, nelec
      !
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ierr
      !
      ! ... then selected tags are read from the other sections
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_bands_info( EF = ef, EF_UP=ef_up, EF_DW=ef_dw, &
            TWO_FERMI_ENERGIES=two_fermi_energies, NELEC=nelec, IERR=ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr > 0 ) RETURN
      !
      IF (ionode) THEN
         !
         IF (.NOT. two_fermi_energies) THEN
            ef = ef * e2
            ef_up = 0.d0
            ef_dw = 0.d0
         ELSE
            ef = 0.d0
            ef_up = ef_up * e2
            ef_dw = ef_dw * e2
         END IF
         !
      END IF
      !
      CALL mp_bcast( two_fermi_energies, ionode_id, intra_image_comm )
      CALL mp_bcast( ef, ionode_id, intra_image_comm )
      CALL mp_bcast( ef_up, ionode_id, intra_image_comm )
      CALL mp_bcast( ef_dw, ionode_id, intra_image_comm )
      CALL mp_bcast( nelec, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_ef
    !
#ifdef __XSD
    !----------------------------------------------------------------------------------------
    SUBROUTINE readschema_ef ( band_struct_obj )
    !----------------------------------------------------------------------------------------
       !
       USE ener,  ONLY            : ef, ef_up, ef_dw
       USE klist, ONLY            : two_fermi_energies, nelec
       USE qes_types_module, ONLY : band_structure_type 
       ! 
       IMPLICIT NONE 
       ! 
       TYPE ( band_structure_type ),INTENT(IN)      :: band_struct_obj 
       ! 
       two_fermi_energies = band_struct_obj%two_fermi_energies_ispresent 
       nelec = band_struct_obj%nelec
       IF ( two_fermi_energies) THEN 
          ef_up = band_struct_obj%two_fermi_energies(1) 
          ef_dw = band_struct_obj%two_fermi_energies(2)
       ELSE IF ( band_struct_obj%fermi_energy_ispresent ) THEN 
          ef = band_struct_obj%fermi_energy
       END IF 
    END SUBROUTINE readschema_ef 
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_exx( ierr )
      !------------------------------------------------------------------------
      !
      ! ... read EXX variables
      !
      USE funct,                ONLY : set_exx_fraction, set_screening_parameter, &
                                       set_gau_parameter, enforce_input_dft, start_exx
      USE exx,                  ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut, ecutfock
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      REAL(DP) :: exx_fraction, screening_parameter, gau_parameter
      LOGICAL :: exx_is_active, found
      !
      IF ( ionode ) THEN
         CALL qexml_read_exx( X_GAMMA_EXTRAPOLATION=x_gamma_extrapolation, &
              NQX1=nq1, NQX2=nq2, NQX3=nq3, EXXDIV_TREATMENT=exxdiv_treatment, &
              YUKAWA = yukawa, ECUTVCUT=ecutvcut, EXX_FRACTION=exx_fraction, &
              SCREENING_PARAMETER=screening_parameter, GAU_PARAMETER=gau_parameter, &
              EXX_IS_ACTIVE=exx_is_active, ECUTFOCK=ecutfock, FOUND=found, IERR=ierr )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr > 0 ) RETURN
      !
      CALL mp_bcast( found, ionode_id, intra_image_comm )
      !
      IF ( .NOT. found ) RETURN
      !
      CALL mp_bcast( x_gamma_extrapolation, ionode_id, intra_image_comm )
      CALL mp_bcast( nq1, ionode_id, intra_image_comm )
      CALL mp_bcast( nq2, ionode_id, intra_image_comm )
      CALL mp_bcast( nq3, ionode_id, intra_image_comm )
      CALL mp_bcast( exxdiv_treatment, ionode_id, intra_image_comm )
      CALL mp_bcast( yukawa, ionode_id, intra_image_comm )
      CALL mp_bcast( ecutvcut, ionode_id, intra_image_comm )
      CALL mp_bcast( exx_fraction, ionode_id, intra_image_comm )
      CALL mp_bcast( screening_parameter, ionode_id, intra_image_comm )
      CALL mp_bcast( gau_parameter, ionode_id, intra_image_comm )
      CALL mp_bcast( exx_is_active, ionode_id, intra_image_comm )
      CALL mp_bcast( ecutfock, ionode_id, intra_image_comm )
      !
      CALL set_exx_fraction(exx_fraction)
      CALL set_screening_parameter(screening_parameter)
      CALL set_gau_parameter(gau_parameter)
      IF (exx_is_active) CALL start_exx( ) 
      !
      RETURN
      !
    END SUBROUTINE read_exx
    !
#ifdef __XSD
    !------------------------------------------------------------------------
    SUBROUTINE readschema_exx ( hybrid_obj) 
    !------------------------------------------------------------------------
      ! 
      USE funct,                ONLY : set_exx_fraction, set_screening_parameter, &
                                      set_gau_parameter, enforce_input_dft, start_exx
      USE exx,                  ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut, ecutfock
      ! 
      USE  qes_types_module,   ONLY : hybrid_type 
      IMPLICIT NONE
      ! 
      TYPE ( hybrid_type), INTENT(IN)               :: hybrid_obj 
      ! 
      x_gamma_extrapolation = hybrid_obj%x_gamma_extrapolation
      nq1 = hybrid_obj%qpoint_grid%nqx1
      nq2 = hybrid_obj%qpoint_grid%nqx2
      nq3 = hybrid_obj%qpoint_grid%nqx3
      CALL set_exx_fraction( hybrid_obj%exx_fraction) 
      CALL set_screening_parameter ( hybrid_obj%screening_parameter) 
      ecutvcut = hybrid_obj%ecutvcut
      ecutfock = hybrid_obj%ecutfock
      CALL start_exx() 
    END SUBROUTINE  readschema_exx 
#endif    
    !------------------------------------------------------------------------
    SUBROUTINE read_esm( ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads only nelec and ef
      !
      USE esm, ONLY : esm_nfit, esm_efield, esm_w, esm_a, esm_bc
      !
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ierr
      !
      ! ... then selected tags are read from the other sections
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_esm( ESM_NFIT = esm_nfit, ESM_EFIELD = esm_efield, &
           ESM_W = esm_w, ESM_A = esm_a, ESM_BC = esm_bc, IERR=ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr > 0 ) RETURN
      !
      CALL mp_bcast( esm_nfit,    ionode_id, intra_image_comm )
      CALL mp_bcast( esm_efield,  ionode_id, intra_image_comm )
      CALL mp_bcast( esm_w,       ionode_id, intra_image_comm )
      CALL mp_bcast( esm_a,       ionode_id, intra_image_comm )
      CALL mp_bcast( esm_bc,      ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_esm
    !
#ifdef __XSD
    !-----------------------------------------------------------------------------------
    SUBROUTINE readschema_esm ( esm_obj ) 
    !-----------------------------------------------------------------------------------
       ! 
       USE esm, ONLY : esm_nfit, esm_efield, esm_w, esm_a, esm_bc  
       USE qes_types_module, ONLY : esm_type
       ! 
       IMPLICIT NONE 
       ! 
       TYPE ( esm_type ), INTENT(IN)    :: esm_obj 
       ! 
       esm_nfit =   esm_obj%nfit
       esm_efield = esm_obj%efield
       esm_w      = esm_obj%w 
       esm_bc     = esm_obj%bc 
       esm_a      = 0.d0 
    END SUBROUTINE readschema_esm 
    !   
#endif
    !------------------------------------------------------------------------
    SUBROUTINE read_( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... this is a template for a "read section" subroutine
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: idum
      !
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "" )
         !
         CALL iotk_scan_end( iunpun, "" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( idum, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_
    !
    !----------------------------------------------------------------------------
    SUBROUTINE gk_l2gmap( ngm, ig_l2g, ngk, igk, igk_l2g )
      !----------------------------------------------------------------------------
      !
      ! ... This subroutine maps local G+k index to the global G vector index
      ! ... the mapping is used to collect wavefunctions subsets distributed
      ! ... across processors.
      ! ... Written by Carlo Cavazzoni
      !
      IMPLICIT NONE
      !
      ! ... Here the dummy variables
      !
      INTEGER, INTENT(IN)  :: ngm, ngk, igk(ngk), ig_l2g(ngm)
      INTEGER, INTENT(OUT) :: igk_l2g(ngk)
      INTEGER              :: ig
      !
      ! ... input: mapping between local and global G vector index
      !
      DO ig = 1, ngk
         !
         igk_l2g(ig) = ig_l2g(igk(ig))
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE gk_l2gmap
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
      IMPLICIT NONE
      !
      ! ... Here the dummy variables
      !
      INTEGER,           INTENT(IN)  :: npw_g, ngk_g, ngk
      INTEGER,           INTENT(IN)  :: igk_l2g(ngk)
      INTEGER, OPTIONAL, INTENT(OUT) :: igwk(ngk_g), igk_l2g_kdip(ngk)
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
      !
      DO ig = 1, ngk
         !
         itmp(igk_l2g(ig)) = igk_l2g(ig)
         !
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
            !
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
      IF ( PRESENT( igk_l2g_kdip ) ) THEN
         !
         ALLOCATE( igwk_lup( npw_g ) )
         !
!$omp parallel private(ig_, ig)
!$omp workshare
         igwk_lup = 0
!$omp end workshare
!$omp do
         do ig_ = 1, ngk_g
            igwk_lup(igwk_(ig_)) = ig_
         end do
!$omp end do
!$omp do
         do ig = 1, ngk
            igk_l2g_kdip(ig) = igwk_lup(igk_l2g(ig))
         end do
!$omp end do
!$omp end parallel
         !
         DEALLOCATE( igwk_lup )

      END IF
      !
      DEALLOCATE( itmp, igwk_ )
      !
      RETURN
      !
    END SUBROUTINE gk_l2gmap_kdip
    !
END MODULE pw_restart
