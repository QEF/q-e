!
! Copyright (C) 2016-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE read_file()
  !----------------------------------------------------------------------------
  !
  ! Read data produced by pw.x or cp.x - new xml file and binary files
  ! Wrapper routine for backwards compatibility
  !
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : nwordwfc, iunwfc, wfc_dir, tmp_dir, restart_dir
  USE buffers,              ONLY : open_buffer, close_buffer
  USE wvfct,                ONLY : nbnd, npwx
  USE noncollin_module,     ONLY : npol
  USE pw_restart_new,       ONLY : read_collected_to_evc
  USE control_flags,        ONLY : io_level
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : gcutw
  USE klist,                ONLY : nkstot, nks, xk, wk
  USE lsda_mod,             ONLY : isk
  USE wvfct,                ONLY : nbnd, et, wg
  !
  IMPLICIT NONE 
  INTEGER :: ierr
  LOGICAL :: exst, wfc_is_collected
  CHARACTER( LEN=256 )  :: dirname
  !
  !
  ierr = 0 
  !
  ! ... Read the contents of the xml data file
  !
  dirname = restart_dir( )
  WRITE( stdout, '(/,5x,A,/,5x,A)') &
     'Reading data from directory:', TRIM( dirname )
  !
  CALL read_xml_file ( wfc_is_collected )
  !
  ! ... more initializations: pseudopotentials / G-vectors / FFT arrays /
  ! ... charge density / potential / ... , but not KS orbitals
  !
  CALL post_xml_init ( )
  !
  ! ... initialization of KS orbitals
  !
  ! ... distribute across pools k-points and related variables.
  ! ... nks is defined by the following routine as the number 
  ! ... of k-points in the current pool
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  CALL poolscatter( nbnd, nkstot, et, nks, et )
  CALL poolscatter( nbnd, nkstot, wg, nks, wg )
  !
  ! ... allocate_wfc_k also computes no. of plane waves and k+G indices
  ! ... FIXME: the latter should be read from file, not recomputed
  !
  CALL allocate_wfc_k()
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  ! ... io_level = 1 so that a real file is opened
  !
  wfc_dir = tmp_dir
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  !
  ! ... read wavefunctions in collected format, writes them to file
  ! ... FIXME: likely not a great idea
  !
  IF ( wfc_is_collected ) CALL read_collected_to_evc(dirname) 
  !
  CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE read_file
!
!----------------------------------------------------------------------------
SUBROUTINE read_xml_file ( wfc_is_collected )
  !----------------------------------------------------------------------------
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
  USE io_global,       ONLY : stdout, ionode, ionode_id
  USE io_files,        ONLY : psfile, pseudo_dir, pseudo_dir_cur, &
                              restart_dir, xmlfile
  USE mp_global,       ONLY : nproc_file, nproc_pool_file, &
                              nproc_image_file, ntask_groups_file, &
                              nproc_bgrp_file, nproc_ortho_file
  USE ions_base,       ONLY : nat, nsp, ityp, amass, atm, tau, extfor
  USE cell_base,       ONLY : alat, at, bg, ibrav, celldm, omega
  USE force_mod,       ONLY : force
  USE klist,           ONLY : nks, nkstot, xk, wk, tot_magnetization, &
       nelec, nelup, neldw, smearing, degauss, ngauss, lgauss, ltetra
  USE ktetra,          ONLY : ntetra, tetra_type
  USE start_k,         ONLY : nks_start, xk_start, wk_start, &
       nk1, nk2, nk3, k1, k2, k3
  USE ener,            ONLY : ef, ef_up, ef_dw
  USE electrons_base,  ONLY : nupdwn, set_nelup_neldw
  USE wvfct,           ONLY : npwx, nbnd, et, wg
  USE extfield,        ONLY : forcefield, forcegate, tefield, dipfield, &
       edir, emaxpos, eopreg, eamp, el_dipole, ion_dipole, gate, zgate, &
       relaxz, block, block_1, block_2, block_height
  USE symm_base,       ONLY : nrot, nsym, invsym, s, ft, irt, t_rev, &
                              sname, inverse_s, s_axis_to_cart, &
                              time_reversal, no_t_rev, nosym, checkallsym
  USE ldaU,            ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_lmax, &
                              Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_alpha, &
                              Hubbard_J0, Hubbard_beta, U_projection
  USE funct,           ONLY : set_exx_fraction, set_screening_parameter, &
                              set_gau_parameter, enforce_input_dft,  &
                              start_exx, dft_is_hybrid
  USE london_module,   ONLY : scal6, lon_rcut, in_C6
  USE tsvdw_module,    ONLY : vdw_isolated
  USE exx_base,        ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                              exxdiv_treatment, yukawa, ecutvcut
  USE exx,             ONLY : ecutfock, local_thr
  USE control_flags,   ONLY : noinv, gamma_only, tqr, llondon, ldftd3, &
       lxdm, ts_vdw
  USE Coul_cut_2D,     ONLY : do_cutoff_2D
  USE noncollin_module,ONLY : noncolin, npol, angle1, angle2, bfield, &
       nspin_lsda, nspin_gga, nspin_mag
  USE spin_orb,        ONLY : domag, lspinorb
  USE lsda_mod,        ONLY : nspin, isk, lsda, starting_magnetization,&
       current_spin
  USE realus,          ONLY : real_space
  USE basis,           ONLY : natomwfc
  USE uspp,            ONLY : okvan
  USE paw_variables,   ONLY : okpaw
  !
  USE qes_types_module,ONLY : output_type, parallel_info_type, &
       general_info_type, input_type
  USE qes_libs_module, ONLY : qes_reset
  USE qexsd_module,    ONLY : qexsd_readschema
  USE qexsd_copy,      ONLY : qexsd_copy_parallel_info, &
       qexsd_copy_algorithmic_info, qexsd_copy_atomic_species, &
       qexsd_copy_atomic_structure, qexsd_copy_symmetry, &
       qexsd_copy_basis_set, qexsd_copy_dft, qexsd_copy_efield, &
       qexsd_copy_band_structure, qexsd_copy_magnetization, &
       qexsd_copy_kpoints
  USE qes_bcast_module,ONLY : qes_bcast
  USE mp_images,       ONLY : intra_image_comm
  USE mp,              ONLY : mp_bcast
  !
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: wfc_is_collected
  !
  INTEGER  :: i, is, ik, ierr, dum1,dum2,dum3
  LOGICAL  :: magnetic_sym, lvalid_input, lfixed
  CHARACTER(LEN=20) :: dft_name, vdw_corr, occupations
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
  IF (ionode) THEN
     ierr = qexsd_readschema ( filename, output_obj, parinfo_obj, geninfo_obj, input_obj)
  END IF
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF ( ierr > 0 ) CALL errore ( 'read_xml_file', 'fatal error reading xml file', ierr ) 
  CALL qes_bcast(output_obj, ionode_id, intra_image_comm)
  CALL qes_bcast(parinfo_obj, ionode_id, intra_image_comm)
  CALL qes_bcast(geninfo_obj, ionode_id, intra_image_comm) 
  CALL qes_bcast(input_obj, ionode_id, intra_image_comm)
  !
  ! ... Now read all needed variables from xml objects
  !
  wfc_is_collected = output_obj%band_structure%wf_collected
  lvalid_input = (TRIM(input_obj%tagname) == "input")
  !
  CALL qexsd_copy_parallel_info (parinfo_obj, nproc_file, &
       nproc_pool_file, nproc_image_file, ntask_groups_file, &
       nproc_bgrp_file, nproc_ortho_file)
  !
  pseudo_dir_cur = restart_dir ( )
  CALL qexsd_copy_atomic_species ( output_obj%atomic_species, &
       nsp, atm, amass, angle1, angle2, starting_magnetization, &
       psfile, pseudo_dir ) 
  IF ( pseudo_dir == ' ' ) pseudo_dir=pseudo_dir_cur
  !! Atomic structure section
  !! tau and ityp are allocated inside qexsd_copy_atomic_structure
  !
  CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp, &
       atm, nat, tau, ityp, alat, at(:,1), at(:,2), at(:,3), ibrav )
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
       dum1,dum2,dum3, ngm_g, ngms_g, npwx, bg(:,1), bg(:,2), bg(:,3) )
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
       lda_plus_U, lda_plus_U_kind, U_projection, Hubbard_l, Hubbard_lmax, &
       Hubbard_U, Hubbard_J0, Hubbard_alpha, Hubbard_beta, Hubbard_J, &
       vdw_corr, scal6, lon_rcut, vdw_isolated )
  !! More DFT initializations
  CALL set_vdw_corr ( vdw_corr, llondon, ldftd3, ts_vdw, lxdm )
  CALL enforce_input_dft ( dft_name, .TRUE. )
  IF ( dft_is_hybrid() ) THEN
     ecutvcut=ecutvcut*e2
     ecutfock=ecutfock*e2
     CALL set_exx_fraction( exx_fraction ) 
     CALL set_screening_parameter ( screening_parameter )
     CALL start_exx ()
  END IF
  !! Band structure section
  !! et and wg are allocated inside qexsd_copy_band_structure
  CALL qexsd_copy_band_structure( output_obj%band_structure, lsda, &
       nkstot, isk, natomwfc, nbnd, nupdwn(1), nupdwn(2), nelec, xk, &
       wk, wg, ef, ef_up, ef_dw, et )
  ! convert to Ry
  ef = ef*e2
  ef_up = ef_up*e2
  ef_dw = ef_dw*e2
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
  CALL set_spin_vars( lsda, noncolin, lspinorb, domag, &
         npol, nspin, nspin_lsda, nspin_mag, nspin_gga, current_spin )
  !! Information for generating k-points and occupations
  CALL qexsd_copy_kpoints( output_obj%band_structure, &
       nks_start, xk_start, wk_start, nk1, nk2, nk3, k1, k2, k3, &
       occupations, smearing, degauss )
  !
  CALL set_occupations( occupations, smearing, degauss, &
       lfixed, ltetra, tetra_type, lgauss, ngauss )
  IF (ltetra) ntetra = 6* nk1 * nk2 * nk3 
  IF (lfixed) CALL errore('read_file','bad occupancies',1)
  ! FIXME: is this really needed? do we use nelup and neldw?
  IF ( lfixed .AND. lsda ) &
       CALL set_nelup_neldw(tot_magnetization, nelec, nelup, neldw) 
  !! Symmetry section
  ALLOCATE ( irt(48,nat) )
  IF ( lvalid_input ) THEN 
     CALL qexsd_copy_symmetry ( output_obj%symmetries, &
          nsym, nrot, s, ft, sname, t_rev, invsym, irt, &
          noinv, nosym, no_t_rev, input_obj%symmetry_flags )

     CALL qexsd_copy_efield ( input_obj%electric_field, &
          tefield, dipfield, edir, emaxpos, eopreg, eamp, &
          gate, zgate, block, block_1, block_2, block_height, relaxz )

  ELSE 
     CALL qexsd_copy_symmetry ( output_obj%symmetries, &
          nsym, nrot, s, ft, sname, t_rev, invsym, irt, &
          noinv, nosym, no_t_rev )
  ENDIF
  !! More initialization needed for symmetry
  magnetic_sym = noncolin .AND. domag
  time_reversal = (.NOT.magnetic_sym) .AND. (.NOT.noinv) 
  CALL inverse_s()
  CALL s_axis_to_cart()
  !! symmetry check - FIXME: is this needed?
  IF (nat > 0) CALL checkallsym( nat, tau, ityp)
  !! Algorithmic info
  do_cutoff_2D = (output_obj%boundary_conditions%assume_isolated == "2D")
  CALL qexsd_copy_algorithmic_info ( output_obj%algorithmic_info, &
       real_space, tqr, okvan, okpaw )
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
!----------------------------------------------------------------------------
SUBROUTINE post_xml_init (  )
  !----------------------------------------------------------------------------
  !
  ! ... Various initializations needed to start a calculation:
  ! ... pseudopotentials, G vectors, FFT arrays, rho, potential
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE uspp_param,           ONLY : upf
  USE read_pseudo_mod,      ONLY : readpp
  USE uspp,                 ONLY : becsum
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE paw_onecenter,        ONLY : paw_potential
  USE dfunct,               ONLY : newd
  USE funct,                ONLY : get_inlc, get_dft_name
  USE ldaU,                 ONLY : lda_plus_u, eth, init_lda_plus_u, U_projection
  USE esm,                  ONLY : do_comp_esm, esm_init
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_fact 
  USE ions_base,            ONLY : nat, nsp, tau, ityp
  USE recvec_subs,          ONLY : ggen, ggens
  USE gvect,                ONLY : gg, ngm, g, gcutm, mill, ngm_g, ig_l2g, &
                                   eigts1, eigts2, eigts3, gstart, gshells
  USE gvecs,                ONLY : ngms, gcutms 
  USE fft_rho,              ONLY : rho_g2r
  USE fft_base,             ONLY : dfftp, dffts
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE io_rho_xml,           ONLY : read_scf
  USE vlocal,               ONLY : strf
  USE control_flags,        ONLY : gamma_only
  USE control_flags,        ONLY : ts_vdw, tqr, tq_smoothing, tbeta_smoothing
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE wvfct,                ONLY : nbnd, nbndx, et, wg
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin
  USE spin_orb,             ONLY : lspinorb
  USE cell_base,            ONLY : at, bg, set_h_ainv
  USE symm_base,            ONLY : d1, d2, d3
  USE realus,               ONLY : betapointlist, generate_qpointlist, &
                                   init_realspace_vars,real_space
  !
  IMPLICIT NONE
  !
  INTEGER  :: inlc
  REAL(DP) :: ehart, etxc, vtxc, etotefield, charge
  CHARACTER(LEN=20) :: dft_name
  !
  ! ... set G cutoffs and cell factor (FIXME: from setup.f90?)
  !
  CALL set_gcut()
  if (cell_factor == 0.d0) cell_factor = 1.D0
  nbndx = nbnd
  !
  ! ... read pseudopotentials
  ! ... the following call prevents readpp from setting dft from PP files
  !
  dft_name = get_dft_name ()
  CALL readpp ( dft_name )
  !
  ! ... misc PP initialization (from setup.f90)
  !
  okpaw = ANY ( upf(1:nsp)%tpawp )
  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !! average_pp must be called before init_lda_plus_u
  IF ( lda_plus_u ) CALL init_lda_plus_u ( upf(1:nsp)%psd, noncolin )
  !
  ! ... allocate memory for G- and R-space fft arrays (from init_run.f90)
  !
  CALL pre_init()
  ! NB: data_structure uses k-points to compute gkcut
  CALL data_structure ( gamma_only )
  CALL allocate_fft()
  CALL ggen ( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart ) 
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms ) 
  CALL gshells ( lmovecell ) 
  !
  IF (do_comp_esm) CALL esm_init()
  IF (do_cutoff_2D) CALL cutoff_fact()
  !
  ! ... allocate the potentials
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
     CALL d_matrix(d1,d2,d3)
  ENDIF
  !
  ! ... read the charge density in G-space
  !
  CALL read_scf( rho, nspin, gamma_only )
  !
  ! ... bring the charge density to real space
  !
  CALL rho_g2r ( dfftp, rho%of_g, rho%of_r )
  !
  ! ... re-compute the local part of the pseudopotential vltot and
  ! ... the core correction charge (if any) - from hinit0.f90
  !
  CALL init_vloc()
  IF (tbeta_smoothing) CALL init_us_b0()
  IF (tq_smoothing) CALL init_us_0()
  CALL init_us_1()
  IF ( lda_plus_U .AND. ( U_projection == 'pseudo' ) ) CALL init_q_aeps()
  CALL init_at_1()
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2,&
                   dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  CALL setlocal()
  CALL set_rhoc()
  !
  ! ... for real-space PP's
  !
  IF ( tqr ) CALL generate_qpointlist()
  IF (real_space ) THEN
     CALL betapointlist()
     CALL init_realspace_vars()
     WRITE (stdout,'(5X,"Real space initialisation completed")')    
  ENDIF
  !
  ! ... recalculate the potential - FIXME: couldn't make ts-vdw work
  !
  IF ( ts_vdw) THEN
     ! CALL tsvdw_initialize()
     ! CALL set_h_ainv()
     CALL infomsg('read_file_new','*** vdW-TS term will be missing in potential ***')
     ts_vdw = .false.
  END IF
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
       ehart, etxc, vtxc, eth, etotefield, charge, v )
  !
  ! ... More PAW and USPP initializations
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  CALL newd()
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_gcut()
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi, eps8
      USE cell_base, ONLY : alat, tpiba, tpiba2
      USE gvect,     ONLY : ecutrho, gcutm
      USE gvecs,     ONLY : gcutms, dual, doublegrid
      USE gvecw,     ONLY : gcutw, ecutwfc
      !
      !
      ! ... Set the units in real and reciprocal space
      !
      tpiba  = 2.D0 * pi / alat
      tpiba2 = tpiba**2
      !
      ! ... Compute the cut-off of the G vectors
      !
      gcutw =        ecutwfc / tpiba2
      gcutm = dual * ecutwfc / tpiba2
      ecutrho=dual * ecutwfc
      !
      doublegrid = ( dual > 4.0_dp + eps8 )
      IF ( doublegrid ) THEN
         gcutms = 4.D0 * ecutwfc / tpiba2
      ELSE
         gcutms = gcutm
      END IF
      !
    END SUBROUTINE set_gcut
    !
  END SUBROUTINE post_xml_init
