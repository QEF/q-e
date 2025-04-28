!
! Copyright (C) 2016-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE read_file()
  !----------------------------------------------------------------------------
  !
  ! Wrapper routine, for backwards compatibility: reads the xml file,
  ! then reads the wavefunctions in "collected" format and writes them
  ! into "distributed" format, forcing write to file (not to buffer).
  ! NOT TO BE USED IN NEW CODE. Use "read_file_new" instead.
  !
  USE io_global,        ONLY : stdout
  USE control_flags,    ONLY : io_level
  USE buffers,          ONLY : open_buffer, close_buffer, save_buffer
  USE io_files,         ONLY : nwordwfc, iunwfc, restart_dir
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : npol
  USE klist,            ONLY : nks
  USE wavefunctions,    ONLY : evc
  USE pw_restart_new,   ONLY : read_collected_wfc
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  LOGICAL :: exst, wfc_is_collected
  !
  wfc_is_collected = .true.
  CALL read_file_new( wfc_is_collected )
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  ! ... io_level = 1 so that a real file is opened
  !
  nwordwfc = nbnd*npwx*npol
  IF ( io_level /= 0 ) io_level = 1
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  !
  ! ... read wavefunctions in collected format, write them to file
  !
  IF ( wfc_is_collected ) THEN
     !
     WRITE( stdout, '(5x,A)') &
          'Reading collected, re-writing distributed wavefunctions'
     DO ik = 1, nks
        CALL read_collected_wfc ( restart_dir(), ik, evc )
        CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     END DO
     !
  ELSE
     WRITE( stdout, '(5x,A)') &
          'read_file: Wavefunctions in collected format not available'
  END IF
  !
  IF ( io_level /= 0 ) CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE read_file
!
!----------------------------------------------------------------------------
SUBROUTINE read_file_ph( needwf_ph )
  !----------------------------------------------------------------------------
  !
  ! Wrapper routine, for compatibility with the phonon code: as "read_file",
  ! but pool parallelization is done just after the reading of the xml file,
  ! before reading the wavefunction files. To be used ONLY for codes that 
  ! can split processors into pools at run-time depending upon the number
  ! of k-points (unless the number of pools is explicitly specified) 
  !
  USE io_global,        ONLY : stdout
  USE control_flags,    ONLY : io_level
  USE buffers,          ONLY : open_buffer, close_buffer, save_buffer
  USE io_files,         ONLY : nwordwfc, iunwfc, wfc_dir, tmp_dir, restart_dir
  USE wvfct,            ONLY : nbnd, npwx, et, wg
  USE noncollin_module, ONLY : npol
  USE klist,            ONLY : nkstot, nks, xk, wk
  USE lsda_mod,         ONLY : isk
  USE wavefunctions,    ONLY : evc
  USE pw_restart_new,   ONLY : read_collected_wfc
  USE fft_base,         ONLY : dffts
  !
  USE pw_restart_new,   ONLY : read_xml_file
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  LOGICAL :: exst, wfc_is_collected
  LOGICAL, INTENT(IN) :: needwf_ph
  !
  WRITE( stdout, '(/,5x,A)') &
       'Reading xml data from directory:', TRIM( restart_dir() )
  !
  ! ... Read the contents of the xml data file
  !
  CALL read_xml_file ( wfc_is_collected )
  !
  ! Guess parallelization on the basis of data from the scf calculation
  ! Must be done before post_xml_init is called
  !
  CALL setup_para ( dffts%nr3, nkstot, nbnd )
  !
  ! ... more initializations: pseudopotentials / G-vectors / FFT arrays /
  ! ... charge density / potential / ... , but not KS orbitals
  !
  CALL post_xml_init ( )
  !
  ! ... distribute across pools k-points and related variables.
  ! ... nks is defined by the following routine as the number 
  ! ... of k-points in the current pool
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  CALL poolscatter( nbnd, nkstot, et, nks, et )
  !$acc update device(et)
  CALL poolscatter( nbnd, nkstot, wg, nks, wg )
  !
  ! ... allocate_wfc_k also computes no. of plane waves and k+G indices
  ! ... FIXME: the latter should be read from file, not recomputed
  !
  CALL allocate_wfc_k()
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  !
  wfc_dir = tmp_dir
  !
  IF ( wfc_is_collected ) THEN
     !
     nwordwfc = nbnd*npwx*npol
     CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
     !
     ! ... read wavefunctions in collected format, write them to file or buffer
     !
     WRITE( stdout, '(5x,A)') &
          'Reading collected, re-writing distributed wavefunctions in '//TRIM(wfc_dir)
     DO ik = 1, nks
        CALL read_collected_wfc ( restart_dir(), ik, evc )
        CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     END DO
     !
  ELSE
     !
     IF ( needwf_ph ) THEN
        CALL errore ('read_file_ph',' Wavefunctions in collected format not available',1)
     ELSE
        WRITE( stdout, '(5x,A)') 'read_file_ph: Wavefunctions in collected format not needed'
     ENDIF
     !
  END IF
  !
  IF ( io_level /= 0 ) CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE read_file_ph
!
!----------------------------------------------------------------------------
SUBROUTINE read_file_new ( needwf )
  !----------------------------------------------------------------------------
  !
  ! Reads xml data file produced by pw.x or cp.x;
  ! performs initializations related to the contents of the xml file;
  ! if needwf=.t. performs wavefunction-related initialization as well.
  ! Does not actually read wfcs. Returns in "needwf" info on the wfc file
  !
  USE io_global,      ONLY : stdout
  USE io_files,       ONLY : nwordwfc, iunwfc, wfc_dir, tmp_dir, restart_dir
  USE gvecw,          ONLY : gcutw
  USE klist,          ONLY : nkstot, nks, xk, wk
  USE lsda_mod,       ONLY : isk
  USE wvfct,          ONLY : nbnd, et, wg
  USE pw_restart_new, ONLY : read_xml_file
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(INOUT) :: needwf
  !
  LOGICAL :: wfc_is_collected
  !
  WRITE( stdout, '(/,5x,A)') &
       'Reading xml data from directory:', TRIM( restart_dir() )
  !
  ! ... Read the contents of the xml data file
  !
  CALL read_xml_file ( wfc_is_collected )
  !
  ! ... more initializations: pseudopotentials / G-vectors / FFT arrays /
  ! ... charge density / potential / ... , but not KS orbitals
  !
  CALL post_xml_init ( )
  !
  IF ( needwf ) THEN
     IF ( .NOT. wfc_is_collected ) WRITE( stdout, '(5x,A)') &
          'read_file_new: Wavefunctions not in collected format?!?'
     !
     ! ... initialization of KS orbitals
     !
     wfc_dir = tmp_dir
     !
     ! ... distribute across pools k-points and related variables.
     ! ... nks is defined by the following routine as the number 
     ! ... of k-points in the current pool
     !
     CALL divide_et_impera( nkstot, xk, wk, isk, nks )
     CALL poolscatter( nbnd, nkstot, et, nks, et )
     !$acc update device(et)
     CALL poolscatter( nbnd, nkstot, wg, nks, wg )
     !
     ! ... allocate_wfc_k also computes no. of plane waves and k+G indices
     ! ... FIXME: the latter should be read from file, not recomputed
     !
     CALL allocate_wfc_k()
     !
  END IF
  needwf = wfc_is_collected
  !
END SUBROUTINE read_file_new
!----------------------------------------------------------------------------
SUBROUTINE post_xml_init (  )
  !----------------------------------------------------------------------------
  !
  ! ... Various initializations needed to start a calculation:
  ! ... pseudopotentials, G vectors, FFT arrays, rho, potential
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE uspp_param,           ONLY : upf, nhm, nsp
  USE read_pseudo_mod,      ONLY : readpp
  USE uspp,                 ONLY : becsum, allocate_uspp
  USE paw_variables,        ONLY : okpaw, ddd_paw
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE paw_onecenter,        ONLY : paw_potential
  USE dfunct,               ONLY : newd
  USE funct,                ONLY : get_dft_name
  USE ldaU,                 ONLY : lda_plus_u, eth, init_hubbard, Hubbard_projectors, &
                                   lda_plus_u_kind
  USE esm,                  ONLY : do_comp_esm, esm_init
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_fact 
  USE ions_base,            ONLY : nat, nsp, tau, ityp
  USE cell_base,            ONLY : omega
  USE recvec_subs,          ONLY : ggen, ggens
  USE gvect,                ONLY : ecutrho, gg, ngm, g, gcutm, mill, ngm_g, &
                                   ig_l2g, eigts1, eigts2, eigts3, gstart, gshells
  USE gvecs,                ONLY : ngms, gcutms 
  USE gvecw,                ONLY : ecutwfc
  USE fft_rho,              ONLY : rho_g2r
  USE fft_base,             ONLY : dfftp, dffts
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE io_rho_xml,           ONLY : read_scf
  USE vlocal,               ONLY : strf
  USE control_flags,        ONLY : gamma_only, use_gpu
  USE control_flags,        ONLY : ts_vdw, tqr, tq_smoothing, tbeta_smoothing
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE wvfct,                ONLY : nbnd, nbndx, et, wg
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, lspinorb, domag
  USE cell_base,            ONLY : at, bg, set_h_ainv
  USE symm_base,            ONLY : d1, d2, d3
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE realus,               ONLY : betapointlist, generate_qpointlist, &
                                   init_realspace_vars,real_space
  USE solvmol,              ONLY : nsolV, solVs
  USE read_solv_module,     ONLY : read_solvents
  USE rism_module,          ONLY : rism_tobe_alive, rism_pot3d
  USE rism3d_facade,        ONLY : lrism3d, rism3d_initialize, rism3d_read_to_restart
  USE xc_lib,               ONLY : xclib_dft_is_libxc, xclib_init_libxc
  USE atwfc_mod,            ONLY : init_tab_atwfc
  USE beta_mod,             ONLY : init_tab_beta
  USE klist,                ONLY : qnorm
  !
  USE tsvdw_module,         ONLY : tsvdw_initialize
  USE xc_lib,               ONLY : xclib_dft_is
  IMPLICIT NONE
  !
  REAL(DP) :: ehart, etxc, vtxc, etotefield, charge, qmax
  CHARACTER(LEN=37) :: dft_name
  INTEGER :: ierr
  !
  ! ... initialize Libxc if needed
  !
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_init_libxc( nspin, domag )
  !
  ! ... set G cutoffs and cell factor (FIXME: from setup.f90?)
  !
  CALL set_gcut()
  if (cell_factor == 0.d0) cell_factor = 1.D0
  nbndx = nbnd
  !
  ! ... activate 3D-RISM
  !
  IF ( lrism3d ) THEN
     CALL rism_tobe_alive()
  END IF
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
  !! average_pp must be called before init_hubbard
  IF ( lda_plus_u ) THEN
     CALL init_hubbard ( upf(1:nsp)%psd, nspin, noncolin )
  ENDIF
  !
  ! ... allocate memory for G- and R-space fft arrays (from init_run.f90)
  !
  CALL pre_init()
  ! NB: data_structure uses k-points to compute gkcut
  CALL data_structure ( gamma_only )
  CALL allocate_fft()
  CALL ggen ( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart ) 
  !$acc update device(mill, g, gg)
  !
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms ) 
  CALL gshells ( lmovecell )
  !
  IF (do_comp_esm) CALL esm_init()
  IF (do_cutoff_2D) CALL cutoff_fact()
  !
  ! ... allocate the potentials
  !
  call allocate_uspp(use_gpu,noncolin,lspinorb,tqr,nhm,nsp,nat,nspin)
  CALL allocate_locpot()
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
  ! ... bring tau to real space
  IF  ( xclib_dft_is('meta') ) THEN
     CALL rho_g2r (dfftp, rho%kin_g, rho%kin_r)
  ENDIF
  !
  ! ... re-compute the local part of the pseudopotential vltot and
  ! ... the core correction charge (if any) - from hinit0.f90
  !
  CALL init_vloc()
  IF (tbeta_smoothing) CALL init_us_b0(ecutwfc,intra_bgrp_comm)
  IF (tq_smoothing) CALL init_us_0(ecutrho,intra_bgrp_comm)
  !
  ! qmax is the maximum |q+G|, for all G needed by the charge density
  !
  qmax = (qnorm+sqrt(ecutrho))*cell_factor
  CALL init_us_1(nat, ityp, omega, qmax, intra_bgrp_comm)
  !
  ! fill interpolation table for beta functions 
  ! qmax is the maximum |q+G|, for all G needed by the wavefunctions
  !
  qmax = (qnorm + sqrt(ecutwfc))*cell_factor
  CALL init_tab_beta ( qmax, omega, intra_bgrp_comm, ierr )
  !
  IF ( lda_plus_u .AND. ( Hubbard_projectors == 'pseudo' ) ) CALL init_q_aeps()
  !
  CALL init_tab_atwfc( qmax, omega, intra_bgrp_comm, ierr )
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2,&
                   dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  !$acc update device(eigts1(:,:), eigts2(:,:), eigts3(:,:))
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
  ! ... read info needed for 3D-RISM
  !
  IF ( lrism3d ) THEN
     CALL read_solvents( without_density=.TRUE. )
     CALL rism3d_initialize()
     CALL rism3d_read_to_restart()
  END IF
  !
  ! ... recalculate the potential - FIXME: couldn't make ts-vdw work
  !
  IF ( ts_vdw) THEN
      CALL tsvdw_initialize()
      CALL set_h_ainv()
     !CALL infomsg('read_file_new','*** vdW-TS term will be missing in potential ***')
     !ts_vdw = .false.
  END IF
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
       ehart, etxc, vtxc, eth, etotefield, charge, v )
  !
  ! ... recalculate the solvation potential (3D-RISM)
  !
  IF ( lrism3d ) THEN
     CALL rism_pot3d(rho%of_g(:, 1), v%of_r)
  END IF
  !
  ! ... More PAW and USPP initializations
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_paw)
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
