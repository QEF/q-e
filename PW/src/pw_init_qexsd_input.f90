! Copyright (C) 2002-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
  !--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE pw_init_qexsd_input(obj,obj_tagname)
  !--------------------------------------------------------------------------------------------------------------------
  !  This routine builds an XML input file, taking the values from the variables
  !  contained in input_parameters MODULE. To work correctly it must be called before
  !  the iosys routine deallocates the input parameters. As the data contained
  !  in XML input are organized differently from those in provided by namelist
  !  input, for some values we need some information from other PW modules. 
  !---------------------------------------------------------------------------
  ! first version march 2016
  !---------------------------------------------------------------------------
  USE input_parameters,  ONLY:  title, calculation, restart_mode, prefix, pseudo_dir, outdir, tstress, tprnfor,       &
                                wf_collect, disk_io, max_seconds, conv_thr, etot_conv_thr, forc_conv_thr,             &
                                press, press_conv_thr,verbosity, iprint, ntyp,                                        &
                                atm => atom_label, psfile => atom_pfile, amass => atom_mass, starting_magnetization,  &
                                angle1, angle2, ip_nat => nat, ip_nspin => nspin, ip_ityp => sp_pos, ip_tau => rd_pos,&
                                ip_atomic_positions => atomic_positions, lspinorb, ip_nqx1 => nqx1, ip_nqx2 => nqx2,  &
                                ip_nqx3 => nqx3, ip_ecutfock => ecutfock, ip_ecutvcut => ecutvcut, localization_thr,  &
                                screening_parameter, exx_fraction, x_gamma_extrapolation, exxdiv_treatment,           &
                                ip_lda_plus_u=>lda_plus_u, ip_lda_plus_u_kind => lda_plus_u_kind,                     & 
                                ip_hubbard_u => hubbard_u, ip_hubbard_j0 => hubbard_j0,                               &
                                ip_hubbard_beta => hubbard_beta, ip_hubbard_alpha => hubbard_alpha,                   &
                                ip_hubbard_j => hubbard_j,  starting_ns_eigenvalue, u_projection_type,                &
                                vdw_corr, london, london_s6, london_c6, london_rcut, london_c6, xdm_a1, xdm_a2,       &
                                ts_vdw_econv_thr, ts_vdw_isolated, dftd3_threebody,dftd3_version,                   &
                                ip_noncolin => noncolin, ip_spinorbit => lspinorb,                                    &
                                nbnd, smearing, degauss, ip_occupations=>occupations, tot_charge, tot_magnetization,  &
                                ip_k_points => k_points, ecutwfc, ip_ecutrho => ecutrho, ip_nr1 => nr1, ip_nr2=>nr2,  &
                                ip_nr3 => nr3, ip_nr1s => nr1s, ip_nr2s => nr2s,ip_nr3s => nr3s, ip_nr1b=>nr1b,       &
                                ip_nr2b=>nr2b, ip_nr3b => nr3b,                                                       &
                                ip_diagonalization=>diagonalization, mixing_mode, mixing_beta,                        &
                                mixing_ndim, tqr, tq_smoothing, tbeta_smoothing, electron_maxstep,                    &
                                diago_thr_init, diago_full_acc,                                                       & 
                                diago_cg_maxiter, diago_ppcg_maxiter, diago_david_ndim,                               &
                                nk1, nk2, nk3, k1, k2, k3, nkstot, ip_xk => xk, ip_wk => wk,                          &
                                ion_dynamics, upscale, remove_rigid_rot, refold_pos, pot_extrapolation,               &
                                wfc_extrapolation, ion_temperature, tempw, tolp, delta_t, nraise, ip_dt => dt,        &
                                bfgs_ndim, trust_radius_min, trust_radius_max, trust_radius_ini, w_1, w_2,            &
                                cell_dynamics, wmass, cell_dofree, cell_factor,                                       &
                                ip_nosym => nosym, ip_noinv => noinv, ip_nosym_evc => nosym_evc,                      & 
                                ip_no_t_rev => no_t_rev, ip_force_symmorphic => force_symmorphic,                     &
                                ip_use_all_frac=>use_all_frac, assume_isolated, esm_bc, esm_w, esm_nfit, esm_efield,  & 
                                ip_lfcpopt => lfcpopt, ip_fcp_mu => fcp_mu,                                           &
                                ecfixed, qcutz, q2sigma,                                                              &    
                                tforces, rd_for,                                                                      &
                                rd_if_pos,                                                                               &
                                tionvel, rd_vel,                                                                      &
                                tefield, lelfield, dipfield, edir, emaxpos, eamp, eopreg, efield, efield_cart,gdir,   &
                                lberry,nppstr,nberrycyc,                                                              &
                                nconstr_inp, nc_fields, constr_type_inp, constr_target_inp, constr_inp, tconstr,      &
                                constr_tol_inp, constrained_magnetization, lambda, fixed_magnetization, input_dft,    &
                                tf_inp, ip_ibrav => ibrav,                                                            &
                                gate, zgate, relaxz, block, block_1, block_2, block_height, real_space
!
  USE fixed_occ,         ONLY:  f_inp               
                                
!
  USE kinds,             ONLY:   DP
  USE parameters,        ONLY:   ntypx
  USE constants,         ONLY:   e2,bohr_radius_angs
  USE ions_base,         ONLY:   iob_tau=>tau
  USE cell_base,         ONLY:   cb_at => at, cb_alat => alat, cb_iforceh => iforceh
  USE funct,             ONLY:   get_dft_is_hybrid => dft_is_hybrid, get_inlc,        &
                                 get_dft_is_nonlocc => dft_is_nonlocc, get_nonlocc_name, get_dft_short
  USE uspp_param,        ONLY:   upf
  USE control_flags,     ONLY:   cf_nstep => nstep 
  USE qes_types_module
  USE qes_libs_module
  USE qexsd_module,      ONLY: qexsd_init_atomic_species, qexsd_init_atomic_structure, qexsd_init_dft, &
                               qexsd_init_hybrid, qexsd_init_vdw, qexsd_init_dftU
  USE qexsd_input  
  IMPLICIT NONE
  ! 
  TYPE (input_type),INTENT(OUT)            ::   obj
  CHARACTER(len=*),INTENT(IN)              ::   obj_tagname
  !
  CHARACTER(80)                            ::   dft_name, diagonalization
  CHARACTER(256)                           ::   tagname
  REAL(DP),ALLOCATABLE                     ::   tau(:,:)
  REAL(DP)                                 ::   alat, a1(3), a2(3), a3(3), gamma_xk(3,1), gamma_wk(1)
  INTEGER                                  ::   inlc,nt
  LOGICAL                                  ::   lsda,dft_is_hybrid,dft_is_nonlocc,is_hubbard(ntypx)=.FALSE., ibrav_lattice
  INTEGER                                  ::   Hubbard_l=0,hublmax=0
  INTEGER                                  ::   iexch, icorr, igcx, igcc, imeta, my_vec(6) 
  INTEGER,EXTERNAL                         ::   set_hubbard_l
  INTEGER                                  ::   lung,l 
  CHARACTER,EXTERNAL                       ::   capital
  CHARACTER(len=20)                        ::   dft_shortname
  CHARACTER(len=25)                        ::   dft_longname
  CHARACTER(LEN=80),TARGET                 ::  vdw_corr_, vdw_nonlocc_
  CHARACTER(LEN=80),POINTER                ::  vdw_corr_pointer => NULL(), vdw_nonlocc_pt=>NULL()
  LOGICAL,TARGET                           ::  gate_tgt, block_tgt, relaxz_tgt
  LOGICAL,POINTER                          ::  gate_ptr, block_ptr, relaxz_ptr
  REAL(DP),TARGET                          ::  block_1_tgt, block_2_tgt, block_height_tgt, zgate_tgt
  REAL(DP),POINTER                         ::  block_1_ptr, block_2_ptr, block_height_ptr, zgate_ptr
  TYPE(hybrid_type),POINTER                ::  hybrid_
  TYPE(dftU_type),POINTER                  ::  dftU_
  TYPE(vdW_type),POINTER                   ::  vdW_
  REAL(DP),TARGET                          ::  xdm_a1_, xdm_a2_, lond_s6_, lond_rcut_, ts_vdw_econv_thr_,&
                                               scr_par_, exx_frc_, ecutvcut_, ecut_fock_, loc_thr_     
  REAL(DP),POINTER                         ::  xdm_a1_pt=>NULL(), xdm_a2_pt=>NULL(), lond_s6_pt=>NULL(), &
                                               lond_rcut_pt=>NULL(), ts_vdw_econv_thr_pt=>NULL(),& 
                                               ecut_fock_opt=>NULL(), scr_par_opt=>NULL(), exx_frc_opt=>NULL(), &
                                               ecutvcut_opt=>NULL(), loc_thr_p => NULL()  
  LOGICAL,TARGET                           ::  empirical_vdw, ts_vdw_isolated_, dftd3_threebody_
  LOGICAL,POINTER                          ::  ts_vdw_isolated_pt=>NULL(), dftd3_threebody_pt=>NULL()
  INTEGER,TARGET                           :: dftd3_version_, spin_ns, nbnd_tg, nq1_tg, nq2_tg, nq3_tg  
  INTEGER,POINTER                          :: dftd3_version_pt=>NULL(), nbnd_pt => NULL(), nq1_pt=>NULL(),&
                                              nq2_pt=>NULL(), nq3_pt=>NULL()  
  REAL(DP),ALLOCATABLE                     :: london_c6_(:), hubbard_U_(:), hubbard_alpha_(:), hubbard_J_(:,:),&
                                              hubbard_J0_(:), hubbard_beta_(:),starting_ns_(:,:,:)
  CHARACTER(LEN=3),ALLOCATABLE             :: species_(:)
  INTEGER, POINTER                         :: nr_1,nr_2, nr_3, nrs_1, nrs_2, nrs_3, nrb_1, nrb_2, nrb_3 
  INTEGER,ALLOCATABLE                      :: nr_(:), nrs_(:), nrb_(:) 
  !
  ! 
  NULLIFY (gate_ptr, block_ptr, relaxz_ptr, block_1_ptr, block_2_ptr, block_height_ptr, zgate_ptr, dftU_, vdW_, hybrid_)
  NULLIFY (nr_1,nr_2,nr_3, nrs_1, nrs_2, nrs_3, nrb_1, nrb_2, nrb_3) 

  obj%tagname=TRIM(obj_tagname)
  IF ( ABS(ip_ibrav)  .GT. 0 ) THEN  
     ibrav_lattice = .TRUE. 
  ELSE
     ibrav_lattice = .FALSE. 
  END IF
  !
  !------------------------------------------------------------------------------------------------------------------------
  !                                                 CONTROL VARIABLES ELEMENT
  !------------------------------------------------------------------------------------------------------------------------
  CALL qexsd_init_control_variables(obj%control_variables,title=title,calculation=calculation,                         &
                                    restart_mode=restart_mode,prefix=prefix,pseudo_dir=pseudo_dir,outdir=outdir,       &
                                    stress=tstress,forces=tprnfor, wf_collect=wf_collect,disk_io=disk_io,              &
                                    max_seconds=max_seconds,etot_conv_thr=etot_conv_thr/e2,forc_conv_thr=forc_conv_thr/e2,   &
                                    press_conv_thr=press_conv_thr,verbosity=verbosity,iprint=iprint, NSTEP = cf_nstep )
  !------------------------------------------------------------------------------------------------------------------------
  !                                                 ATOMIC SPECIES                                                      
  !------------------------------------------------------------------------------------------------------------------------
  IF ( ip_noncolin ) THEN 
     CALL qexsd_init_atomic_species(obj%atomic_species, ntyp,atm, psfile, amass, starting_magnetization, angle1, angle2)
  ELSE IF (ip_nspin == 1 ) THEN
     CALL qexsd_init_atomic_species(obj%atomic_species, ntyp,atm, psfile, amass)
  ELSE IF (ip_nspin == 2 ) THEN
     CALL qexsd_init_atomic_species(obj%atomic_species, ntyp,atm, psfile, amass, starting_magnetization)
  END IF
  !------------------------------------------------------------------------------------------------------------------------
  !                                                 ATOMIC STRUCTURE
  !------------------------------------------------------------------------------------------------------------------------
  ALLOCATE (tau(3, ip_nat))
  alat = cb_alat                      !*bohr_radius_angs
  a1 = cb_at(:,1)*alat
  a2 = cb_at(:,2)*alat
  a3 = cb_at(:,3)*alat
  tau(1:3,1:ip_nat) = iob_tau(1:3,1:ip_nat)*alat
  !
  IF ( ibrav_lattice ) THEN 
     CALL qexsd_init_atomic_structure (obj%atomic_structure, ntyp, atm, ip_ityp, ip_nat, tau, &
                                       ALAT = alat, a1 = a1, a2 = a2, a3 = a3 , ibrav = ip_ibrav )
  ELSE 
     CALL qexsd_init_atomic_structure (obj%atomic_structure, ntyp, atm, ip_ityp, ip_nat, tau,     &
                                    alat = sqrt(sum(a1(1:3)*a1(1:3))), A1 = a1, A2 = a2, A3 = a3 , IBRAV = 0 )
  END IF 
  DEALLOCATE ( tau ) 
  ! 
  !--------------------------------------------------------------------------------------------------------------------------
  !                                                   DFT ELEMENT
  !---------------------------------------------------------------------------------------------------------------------------
  IF ( TRIM(input_dft) .NE. "none" ) THEN 
     dft_name=TRIM(input_dft) 
  ELSE 
     dft_shortname = get_dft_short()        
     dft_name=TRIM(dft_shortname)
  END IF

  dft_is_hybrid=get_dft_is_hybrid()
  IF ( dft_is_hybrid) THEN
     ALLOCATE(hybrid_)
     IF (screening_parameter > 0.0_DP) THEN 
        scr_par_ = screening_parameter 
        scr_par_opt => scr_par_ 
     END IF 
     IF ( exx_fraction > 0.0_DP) THEN 
        exx_frc_ = exx_fraction 
        exx_frc_opt => exx_frc_ 
     END IF 
     IF ( ip_ecutfock > 0.0_DP) THEN 
        ecut_fock_ = ip_ecutvcut/e2 
        ecut_fock_opt => ecut_fock_
     END IF 
     IF ( ip_ecutvcut > 0.0_DP) THEN 
        ecutvcut_ = ip_ecutvcut/e2 
        ecutvcut_opt => ecutvcut_ 
     END IF 
     IF (ANY([ip_nqx1, ip_nqx2, ip_nqx3] /= 0)) THEN 
        nq1_tg = ip_nqx1 
        nq2_tg = ip_nqx2 
        nq3_tg = ip_nqx3
        nq1_pt => nq1_tg
        nq2_pt => nq2_tg
        nq3_pt => nq3_tg
     END IF 
     IF (localization_thr .GT. 0._DP) THEN 
        loc_thr_ = localization_thr
        loc_thr_p => loc_thr_ 
     END IF 
     CALL qexsd_init_hybrid(hybrid_, dft_is_hybrid, NQ1 = ip_nqx1, NQ2= ip_nqx2, NQ3=ip_nqx3,&
                            ECUTFOCK = ecut_fock_opt, EXX_FRACTION = exx_frc_opt,          &
                            SCREENING_PARAMETER = scr_par_opt,  EXXDIV_TREATMENT = exxdiv_treatment,&
                            X_GAMMA_EXTRAPOLATION = x_gamma_extrapolation, ECUTVCUT = ecutvcut_opt, &
                            LOCAL_THR = loc_thr_p )
  END IF
  dft_is_nonlocc=get_dft_is_nonlocc()
  vdw_corr_ = vdw_corr
  IF (london) vdw_corr_ = 'grimme-d2'
  empirical_vdw = .NOT. ( TRIM(vdw_corr_)  == 'none')
  IF (empirical_vdw .OR. dft_is_nonlocc) THEN
    ALLOCATE (vdW_)
    IF ( empirical_vdw ) THEN
        vdw_corr_pointer => vdw_corr_
        SELECT CASE ( TRIM(vdw_corr_))
            CASE ('grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d')
                IF ( london_s6 .NE. 0.75_DP .OR. london_rcut .NE. 200._DP ) THEN
                    lond_s6_ = london_s6
                    lond_s6_pt => lond_s6_
                    lond_rcut_ = london_rcut
                    lond_rcut_pt => lond_rcut_
                END IF
                IF (ANY( london_c6(1:ntyp) .NE. -1._DP )) THEN
                    ALLOCATE (london_c6_(ntyp), species_(ntyp))
                    london_c6_(1:ntyp) = london_c6(1:ntyp)
                    species_(1:ntyp)  = atm(1:ntyp)
                END IF
            CASE ('TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler')
                ts_vdw_isolated_ = ts_vdw_isolated
                ts_vdw_isolated_pt => ts_vdw_isolated_
                ts_vdw_econv_thr_ = ts_vdw_econv_thr
                ts_vdw_econv_thr_pt => ts_vdw_econv_thr_
            CASE ('XDM' , 'xdm')
                xdm_a1_ = xdm_a1
                xdm_a1_pt => xdm_a1_
                xdm_a2_ = xdm_a2
                xdm_a2_pt => xdm_a2_
            CASE ('grimme-d3', 'Grimme-D3', 'DFT-D3', 'dft-d3')
                dftd3_version_ = dftd3_version
                dftd3_version_pt => dftd3_version_
                dftd3_threebody_ = dftd3_threebody
                dftd3_threebody_pt => dftd3_threebody_
        END SELECT
     END IF
     IF (dft_is_nonlocc) THEN
         vdw_nonlocc_ = TRIM(get_nonlocc_name())
         vdw_nonlocc_pt => vdw_nonlocc_
     END IF
     CALL qexsd_init_vdw(vdW_, NON_LOCAL_TERM = vdw_nonlocc_pt, VDW_CORR = vdw_corr_pointer, &
                             TS_THR = ts_vdw_econv_thr_pt, TS_ISOL = ts_vdw_isolated_pt, &
                             LONDON_S6 = lond_s6_pt, LONDON_RCUT = lond_rcut_pt, SPECIES = species_, &
                             XDM_A1 = xdm_a1_pt, XDM_A2 = xdm_a2_pt, DFTD3_VERSION = dftd3_version_pt, &
                             DFTD3_THREEBODY = dftd3_threebody_pt)
  END IF
  !
  IF (ip_lda_plus_u) THEN
     ALLOCATE (dftU_)
     !
     DO nt = 1, ntyp
       !
       is_hubbard(nt) = ip_Hubbard_U(nt)/= 0.0_dp .OR. &
                            ip_Hubbard_alpha(nt) /= 0.0_dp .OR. &
                            ip_Hubbard_J0(nt) /= 0.0_dp .OR. &
                            ip_Hubbard_beta(nt)/= 0.0_dp .OR. &
                            ANY(ip_Hubbard_J(:,nt) /= 0.0_DP)
        IF (is_hubbard(nt)) hublmax = MAX (hublmax, set_hubbard_l(upf(nt)%psd))
        !
     END DO
     IF ( ANY(ip_hubbard_u(1:ntyp) /=0.0_DP)) THEN
        ALLOCATE(hubbard_U_(ntyp))
        hubbard_U_(1:ntyp) = ip_hubbard_u(1:ntyp)
     END IF
     IF (ANY (ip_hubbard_J0 /=0.0_DP)) THEN
        ALLOCATE(hubbard_J0_(ntyp))
        hubbard_J0_ (1:ntyp) = ip_hubbard_J0(1:ntyp)
     END IF
     IF (ANY (ip_hubbard_alpha /=0.0_DP)) THEN
        ALLOCATE(hubbard_alpha_(ntyp))
        hubbard_alpha_ (1:ntyp) = ip_hubbard_alpha(1:ntyp)
     END IF
     IF (ANY (ip_hubbard_beta /=0.0_DP)) THEN
        ALLOCATE(hubbard_beta_(ntyp))
        hubbard_beta_ (1:ntyp) = ip_hubbard_beta(1:ntyp)
     END IF
     IF (ANY (ip_hubbard_J(:,1:ntyp) /=0.0_DP )) THEN
        ALLOCATE(hubbard_J_(3,ntyp))
        hubbard_J_(1:3,1:ntyp) = ip_hubbard_J(1:3,1:ntyp)
      END IF
      IF (ANY(starting_ns_eigenvalue /= -1.0_DP)) THEN
         IF (lsda) THEN
            spin_ns = 2
         ELSE
            spin_ns = 1
         END IF
         ALLOCATE (starting_ns_(2*hublmax+1, spin_ns, ntyp))
         starting_ns_          (1:2*hublmax+1, 1:spin_ns, 1:ntyp) = &
         starting_ns_eigenvalue(1:2*hublmax+1, 1:spin_ns, 1:ntyp)
      END IF
      CALL qexsd_init_dftU(dftU_, NSP = ntyp, PSD = upf(1:ntyp)%psd, SPECIES = atm(1:ntyp), ITYP = ip_ityp(1:ntyp), &
                           IS_HUBBARD = is_hubbard(1:ntyp), LDA_PLUS_U_KIND = ip_lda_plus_u_kind,                   &
                           U_PROJECTION_TYPE=u_projection_type, U=hubbard_U_, J0=hubbard_J0_, NONCOLIN=ip_noncolin, &
                           ALPHA = hubbard_alpha_, BETA = hubbard_beta_, J = hubbard_J_, STARTING_NS = starting_ns_ )
  END IF
  CALL qexsd_init_dft(obj%dft, TRIM(dft_name), hybrid_, vdW_, dftU_)
  IF (ASSOCIATED(hybrid_)) THEN
    CALL qes_reset(hybrid_)
    DEALLOCATE(hybrid_)
  END IF
  IF (ASSOCIATED(vdW_)) THEN
    CALL qes_reset(vdW_)
    DEALLOCATE(vdW_)
  END IF
  IF (ASSOCIATED(dftU_)) THEN
    CALL qes_reset(dftU_)
    DEALLOCATE(dftU_)
  END IF

  !
  !------------------------------------------------------------------------------------------------------------------------
  !                                                   SPIN ELEMENT
  !-------------------------------------------------------------------------------------------------------------------------
  IF (ip_nspin == 2) THEN
     lsda=.TRUE.
  ELSE
     lsda=.FALSE.
  END IF
  CALL qexsd_init_spin(obj%spin, lsda, ip_noncolin, ip_spinorbit)
  !-------------------------------------------------------------------------------------------------------------------------
  !                                                    BANDS ELEMENT
  !-------------------------------------------------------------------------------------------------------------------------
  IF (nbnd /= 0) THEN
     nbnd_tg = nbnd 
     nbnd_pt => nbnd_tg
  END IF 
  IF (tf_inp) THEN
     SELECT CASE (ip_nspin) 
        CASE (2)  
           CALL qexsd_init_bands(obj%bands, nbnd_pt, smearing, degauss, ip_occupations, tot_charge, ip_nspin, &
                                          input_occupations=f_inp(:,1),input_occupations_minority=f_inp(:,2))
        CASE default
           CALL qexsd_init_bands(obj%bands, nbnd_pt, smearing, degauss, ip_occupations, tot_charge, ip_nspin, &
                                                                                input_occupations=f_inp(:,1) )
     END SELECT    
  ELSE 
     IF ( tot_magnetization .LT. 0 ) THEN 
        CALL qexsd_init_bands(obj%bands, nbnd_pt, smearing, degauss, ip_occupations, tot_charge, ip_nspin)
     ELSE
        CALL qexsd_init_bands(obj%bands, nbnd_pt, smearing, degauss, ip_occupations, tot_charge, ip_nspin, &
                              TOT_MAG  = tot_magnetization)
     END IF
  END IF 
  !----------------------------------------------------------------------------------------------------------------------------
  !                                                    BASIS ELEMENT
  !---------------------------------------------------------------------------------------------------------------------------
  IF (ANY([ip_nr1,ip_nr2,ip_nr3] /=0)) THEN 
     ALLOCATE (nr_(3)) 
     nr_ = [ip_nr1,ip_nr2,ip_nr3]
  END IF 
  IF (ANY([ip_nr1s,ip_nr2s,ip_nr3s] /=0)) THEN 
     ALLOCATE (nrs_(3))
     nrs_ = [ip_nr1s,ip_nr2s,ip_nr3s]
  END IF 
  IF (ANY([ip_nr1b,ip_nr2b,ip_nr3b] /=0)) THEN 
     ALLOCATE(nrb_(3)) 
     nrb_ = [ip_nr1b,ip_nr2b,ip_nr3b]
  END IF 

  CALL qexsd_init_basis(obj%basis, ip_k_points, ecutwfc/e2, ip_ecutrho/e2, nr_ , nrs_, nrb_ ) 
  !-----------------------------------------------------------------------------------------------------------------------------
  !                                                    ELECTRON CONTROL
  !------------------------------------------------------------------------------------------------------------------------------
  IF (TRIM(ip_diagonalization) == 'david') THEN 
     diagonalization = 'davidson'
  ELSE 
    diagonalization = ip_diagonalization
  END IF
  CALL qexsd_init_electron_control(obj%electron_control, diagonalization, mixing_mode, mixing_beta, conv_thr/e2,         &
                                   mixing_ndim, electron_maxstep, tqr, real_space, tq_smoothing, tbeta_smoothing, diago_thr_init, &
                                   diago_full_acc, diago_cg_maxiter,  diago_ppcg_maxiter, diago_david_ndim )
  !--------------------------------------------------------------------------------------------------------------------------------
  !                                                   K POINTS IBZ ELEMENT
  !------------------------------------------------------------------------------------------------------------------------------ 
  IF (TRIM(ip_k_points) .EQ. 'gamma' ) THEN 
      gamma_xk(:,1)=[0._DP, 0._DP, 0._DP]
      gamma_wk(1)=1._DP
      CALL qexsd_init_k_points_ibz( obj%k_points_ibz, ip_k_points, calculation, nk1, nk2, nk3, k1, k2, k3, 1,         &
                                    gamma_xk, gamma_wk ,alat,a1,ibrav_lattice) 

  ELSE 
     CALL qexsd_init_k_points_ibz(obj%k_points_ibz, ip_k_points, calculation, nk1, nk2, nk3, k1, k2, k3, nkstot,      &
                                   ip_xk, ip_wk,alat,a1, ibrav_lattice)

  END IF
  !--------------------------------------------------------------------------------------------------------------------------------
  !                                                       ION CONTROL ELEMENT
  !--------------------------------------------------------------------------------------------------------------------------------
  CALL qexsd_init_ion_control(obj%ion_control, ion_dynamics, upscale, remove_rigid_rot, refold_pos,                   &
                              pot_extrapolation, wfc_extrapolation, ion_temperature, tempw, tolp, delta_t, nraise,    &
                              ip_dt, bfgs_ndim, trust_radius_min, trust_radius_max, trust_radius_ini, w_1, w_2)
  !--------------------------------------------------------------------------------------------------------------------------------
  !                                                        CELL CONTROL ELEMENT
  !-------------------------------------------------------------------------------------------------------------------------------
  CALL qexsd_init_cell_control(obj%cell_control, cell_dynamics, press, wmass, cell_factor, cell_dofree, cb_iforceh)
  !---------------------------------------------------------------------------------------------------------------------------------
  !                                SYMMETRY FLAGS
  !------------------------------------------------------------------------------------------------------------------------ 
  obj%symmetry_flags_ispresent = .TRUE.
  CALL qexsd_init_symmetry_flags(obj%symmetry_flags, ip_nosym,ip_nosym_evc, ip_noinv, ip_no_t_rev,                    & 
                                 ip_force_symmorphic, ip_use_all_frac)       
  !------------------------------------------------------------------------------------------------------------------------
  !                              BOUNDARY CONDITIONS
  !---------------------------------------------------------------------------------------------------------------------------- 
  IF (TRIM( assume_isolated ) .EQ. "none" ) THEN
     obj%boundary_conditions_ispresent=.FALSE.
  ELSE 
     obj%boundary_conditions_ispresent = .TRUE.
     IF ( TRIM ( assume_isolated) .EQ. "esm") THEN 
        SELECT CASE (TRIM(esm_bc)) 
          CASE ('pbc', 'bc1' ) 
             CALL qexsd_init_boundary_conditions(obj%boundary_conditions, assume_isolated, esm_bc,&
                                                 ESM_NFIT = esm_nfit, ESM_W = esm_w,ESM_EFIELD = esm_efield)
          CASE ('bc2', 'bc3' ) 
            CALL qexsd_init_boundary_conditions(obj%boundary_conditions, assume_isolated, esm_bc, &
                                                 FCP_OPT = ip_lfcpopt, FCP_MU = ip_fcp_mu, &
                                                 ESM_NFIT = esm_nfit, ESM_W = esm_w,ESM_EFIELD = esm_efield)
        END SELECT 
     ELSE 
        CALL qexsd_init_boundary_conditions(obj%boundary_conditions, assume_isolated) 
     END IF 
  END IF
  !----------------------------------------------------------------------------------------------------------------------------
  !                                                              EKIN FUNCTIONAL 
  !-------------------------------------------------------------------------------------------------------------------------------  
  IF (ecfixed .GT. 1.d-3) THEN
     obj%ekin_functional_ispresent = .TRUE.
     CALL qexsd_init_ekin_functional ( obj%ekin_functional, ecfixed, qcutz, q2sigma)
  ELSE 
     obj%ekin_functional_ispresent = .FALSE.
  END IF
  !-----------------------------------------------------------------------------------------------------------------------------
  !                                                         EXTERNAL FORCES 
  !------------------------------------------------------------------------------------------------------------------------------
  IF ( tforces ) THEN
      obj%external_atomic_forces_ispresent = .TRUE.
      CALL qexsd_init_external_atomic_forces (obj%external_atomic_forces, rd_for,ip_nat)
  ELSE
      obj%external_atomic_forces_ispresent= .FALSE.
  END IF
  !-------------------------------------------------------------------------------------------------------------------------------
  !                                            FREE POSITIONS
  !---------------------------------------------------------------------------------------------------------------------------- 
  IF ( TRIM(calculation) .NE. "scf" .AND. TRIM(calculation) .NE. "nscf" .AND. &
                                           TRIM(calculation) .NE. "bands") THEN
      obj%free_positions_ispresent=.TRUE.
      CALL qexsd_init_free_positions( obj%free_positions, rd_if_pos, ip_nat)
  ELSE
      obj%free_positions_ispresent = .FALSE.
  END IF
  !----------------------------------------------------------------------------------------------------------------------------
  !                                  STARTING IONIC VELOCITIES 
  !-----------------------------------------------------------------------------------------------------------------------------
  IF (tionvel) THEN
     obj%starting_atomic_velocities_ispresent=.TRUE.
     CALL qexsd_init_starting_atomic_velocities(obj%starting_atomic_velocities,tionvel,rd_vel,ip_nat)
  ELSE
     obj%starting_atomic_velocities_ispresent=.FALSE.
  END IF
  !-------------------------------------------------------------------------------------------------------------------------------
  !                                ELECTRIC FIELD
  !--------------------------------------------------------------------------------------------------------------------------- 
  IF (tefield .OR. lelfield .OR. lberry .or. gate ) THEN 
     obj%electric_field_ispresent=.TRUE.
     IF ( gate ) THEN 
         gate_tgt = gate
         gate_ptr => gate_tgt
         zgate_tgt = zgate
         zgate_ptr => zgate_tgt
         block_tgt = block
         block_ptr => block_tgt
         block_1_tgt = block_1
         block_1_ptr => block_1_tgt
         block_2_tgt = block_2
         block_2_ptr => block_2_tgt
         block_height_tgt = block_height
         block_height_ptr => block_height_tgt
         relaxz_tgt = relaxz 
         relaxz_ptr => relaxz_tgt
     END IF
     CALL qexsd_init_electric_field_input(obj%electric_field, tefield, dipfield, lelfield, lberry,       &
                              edir, gdir, emaxpos, eopreg, eamp, efield, efield_cart, nberrycyc, nppstr, &
                              GATE = gate_ptr, ZGATE = zgate_ptr, RELAXZ = relaxz_ptr, BLOCK = block_ptr,&
                              BLOCK_1 = block_1_ptr, BLOCK_2 = block_2_ptr, BLOCK_HEIGHT = block_height_ptr)
  ELSE
     obj%electric_field_ispresent=.FALSE.
  END IF
  !-----------------------------------------------------------------------------------------------------------------------
  !                                     ATOMIC CONSTRAINTS
  !------------------------------------------------------------------------------------------------------------------------ 
  IF (tconstr) THEN
     obj%atomic_constraints_ispresent=.TRUE.
     CALL qexsd_init_atomic_constraints( obj%atomic_constraints, ion_dynamics, tconstr, nconstr_inp,constr_type_inp,  &
                                         constr_tol_inp, constr_target_inp, constr_inp)
  ELSE 
     obj%atomic_constraints_ispresent=.FALSE.
  END IF
  !-----------------------------------------------------------------------------------------------------------------------------
  !                                               SPIN CONSTRAINTS
  !------------------------------------------------------------------------------------------------------------------------------
  
  SELECT CASE (TRIM( constrained_magnetization ))
 
     CASE ("total","total direction") 
          obj%spin_constraints_ispresent=.TRUE.
          CALL qexsd_init_spin_constraints(obj%spin_constraints, constrained_magnetization,lambda,&
                                          fixed_magnetization)
     CASE ("atomic", "atomic direction")
          obj%spin_constraints_ispresent=.TRUE.
          CALL qexsd_init_spin_constraints(obj%spin_constraints, constrained_magnetization, lambda )
     CASE default 
          obj%spin_constraints_ispresent=.FALSE.
  END SELECT
  
  
  obj%lread=.TRUE.
  obj%lwrite=.TRUE.
  ! 
  !
  END SUBROUTINE pw_init_qexsd_input
  !
