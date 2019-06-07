!
! Copyright (C) 2016 Quantum ESPRESSO group
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
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, postfix, &
       tmp_dir, wfc_dir
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
  dirname = TRIM( tmp_dir ) // TRIM( prefix ) // postfix
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
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE ions_base,            ONLY : nat, ityp, tau, extfor
  USE force_mod,            ONLY : force
  USE klist,                ONLY : nks, nkstot
  USE wvfct,                ONLY : nbnd, et, wg
  USE symm_base,            ONLY : irt
  USE extfield,             ONLY : forcefield, tefield, gate, forcegate
  USE io_files,             ONLY : tmp_dir, prefix, postfix
  USE pw_restart_new,       ONLY : pw_read_schema, readschema_dim, &
       readschema_cell, readschema_ions, readschema_planewaves, &
       readschema_spin, readschema_magnetization, readschema_xc, &
       readschema_occupations, readschema_brillouin_zone, &
       readschema_band_structure, readschema_symmetry, readschema_efield, &
       readschema_outputPBC, readschema_exx, readschema_algo
  USE qes_types_module,     ONLY : output_type, parallel_info_type, &
       general_info_type, input_type
  USE qes_libs_module,      ONLY : qes_reset
#if defined(__BEOWULF)
  USE qes_bcast_module,     ONLY : qes_bcast
  USE mp_images,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_bcast
#endif
  !
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: wfc_is_collected

  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym, ierr
  CHARACTER(LEN=256) :: dirname
  LOGICAL            :: lvalid_input
  TYPE ( output_type)                   :: output_obj 
  TYPE (parallel_info_type)             :: parinfo_obj
  TYPE (general_info_type )             :: geninfo_obj
  TYPE (input_type)                     :: input_obj
  !
  !
#if defined(__BEOWULF)
   IF (ionode) THEN
      CALL pw_read_schema ( ierr, output_obj, parinfo_obj, geninfo_obj, input_obj)
   END IF
   CALL mp_bcast(ierr,intra_image_comm)
   IF ( ierr /= 0 ) CALL errore ( 'read_schema', 'unable to read xml file', abs(ierr) ) 
   CALL qes_bcast(output_obj, ionode_id, intra_image_comm)
   CALL qes_bcast(parinfo_obj, ionode_id, intra_image_comm)
   CALL qes_bcast(geninfo_obj, ionode_id, intra_image_comm) 
   CALL qes_bcast(input_obj, ionode_id, intra_image_comm)
#else
  CALL pw_read_schema ( ierr, output_obj, parinfo_obj, geninfo_obj, input_obj)
  IF ( ierr /= 0 ) CALL errore ( 'read_schema', 'unable to read xml file', abs(ierr) ) 
#endif
  wfc_is_collected = output_obj%band_structure%wf_collected
  !
  ! ... here we read the variables that dimension the system
  !
  CALL readschema_cell( output_obj%atomic_structure ) 
  CALL readschema_dim ( parinfo_obj, output_obj%atomic_species, &
       output_obj%atomic_structure, output_obj%symmetries, &
       output_obj%basis_set, output_obj%band_structure ) 
  !
  ! ... until pools are activated, the local number of k-points nks
  ! ... should be equal to the global number nkstot - k-points are replicated
  !
  nks = nkstot
  !
  ! ... allocate space for arrays to be read in this routine
  !
  ! ... atomic positions, forces, symmetries
  !
  IF ( nat < 0 ) CALL errore( 'read_xml_file', 'wrong number of atoms', 1 )
  ALLOCATE( ityp( nat ) )
  ALLOCATE( tau( 3, nat ) )
  ALLOCATE( force ( 3, nat ) )
  ALLOCATE( extfor( 3, nat ) )
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  IF ( gate ) ALLOCATE( forcegate( 3, nat ) )
  ALLOCATE( irt( 48, nat ) )
  !
  ! ... eigenvalues, weights
  !
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  !
  ! ... here we read all the variables defining the system
  !
  lvalid_input = (TRIM(input_obj%tagname) == "input")
  !
  dirname = TRIM( tmp_dir ) // TRIM( prefix ) // postfix ! FIXME
  CALL readschema_ions( output_obj%atomic_structure, output_obj%atomic_species, dirname)
  CALL readschema_planewaves( output_obj%basis_set) 
  CALL readschema_spin( output_obj%magnetization )
  CALL readschema_magnetization (  output_obj%band_structure,  &
       output_obj%atomic_species, output_obj%magnetization )
  CALL readschema_xc (  output_obj%atomic_species, output_obj%dft )
  CALL readschema_occupations( output_obj%band_structure )
  CALL readschema_brillouin_zone( output_obj%symmetries,  output_obj%band_structure )
  CALL readschema_band_structure( output_obj%band_structure )
  IF ( lvalid_input ) THEN 
     CALL readschema_symmetry (  output_obj%symmetries, output_obj%basis_set, input_obj%symmetry_flags )
     CALL readschema_efield ( input_obj%electric_field )
  ELSE 
     CALL readschema_symmetry( output_obj%symmetries,output_obj%basis_set) 
  ENDIF
  CALL readschema_outputPBC ( output_obj%boundary_conditions)
  IF ( output_obj%dft%hybrid_ispresent  ) THEN
     CALL readschema_exx ( output_obj%dft%hybrid )
  END IF
  CALL readschema_algo(output_obj%algorithmic_info )  
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
  USE noncollin_module,     ONLY : noncolin
  USE spin_orb,             ONLY : lspinorb
  USE funct,                ONLY : get_inlc, get_dft_name
  USE kernel_table,         ONLY : initialize_kernel_table
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
  USE cell_base,            ONLY : at, bg, set_h_ainv
  USE symm_base,            ONLY : d1, d2, d3, checkallsym
  USE realus,               ONLY : betapointlist, generate_qpointlist, &
                                   init_realspace_vars,real_space
  !
  IMPLICIT NONE
  !
  INTEGER  :: inlc
  REAL(DP) :: ehart, etxc, vtxc, etotefield, charge
  CHARACTER(LEN=20) :: dft_name
  !
  ! ... check on symmetry (FIXME: is this needed?)
  !
  IF (nat > 0) CALL checkallsym( nat, tau, ityp)
  !
  ! ... set spin variables, G cutoffs, cell factor (FIXME: from setup.f90?)
  !
  CALL set_gcut()
  if (cell_factor == 0.d0) cell_factor = 1.D0
  CALL set_spin_vars ( )
  nbndx = nbnd
  !
  ! ... read pseudopotentials
  ! ... the following call prevents readpp from setting dft from PP files
  !
  dft_name = get_dft_name ()
  CALL readpp ( dft_name )
  !
  ! ... read the vdw kernel table if needed
  !
  inlc = get_inlc()
  IF (inlc > 0 ) CALL initialize_kernel_table(inlc)
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
    !------------------------------------------------------------------------
    SUBROUTINE set_spin_vars( )
      !------------------------------------------------------------------------
      !
      !  Set various spin-related variables
      !
      USE noncollin_module, ONLY : nspin_lsda, nspin_mag, nspin_gga
      USE spin_orb,  ONLY : domag
      USE lsda_mod, ONLY : nspin, current_spin
      !
      IF (nspin /= 2) current_spin = 1
      !
      nspin_mag  = nspin
      nspin_lsda = nspin
      nspin_gga  = nspin
      IF (nspin==4) THEN
        nspin_lsda=1
        IF (domag) THEN
           nspin_gga=2
        ELSE
           nspin_gga=1
           nspin_mag=1
        ENDIF
      ENDIF
      !
    END SUBROUTINE set_spin_vars
    !
  END SUBROUTINE post_xml_init
