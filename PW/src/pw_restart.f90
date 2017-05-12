!
! Copyright (C) 2005-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! TB
! included monopole related variables, search for 'TB'
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
MODULE pw_restart
  !----------------------------------------------------------------------------
  !
  ! ... this module contains methods to read and write data produced by PWscf
  !
  ! ... originally written by Carlo Sbraccia  (2005)
  !
#if defined(__OLDXML)
  !
  USE iotk_module
  !
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
  PUBLIC :: pw_readfile, pw_writefile, pp_check_file
  PUBLIC :: gk_l2gmap, gk_l2gmap_kdip
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
                                       two_fermi_energies, nelup, neldw, ltetra
      USE start_k,              ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       nks_start, xk_start, wk_start
      USE ktetra,               ONLY : ntetra, tetra, tetra_type
      USE klist,                ONLY : ltetra
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
                                       emaxpos, eopreg, eamp, & !TB
                                       monopole, zmon, block, block_1, &
                                       block_2, block_height, relaxz
      USE io_rho_xml,           ONLY : write_scf
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
      INTEGER               :: nkl, nkr, npwx_g
      INTEGER               :: ike, iks, npw_g, ispin, inlc
      INTEGER, EXTERNAL     :: global_kpoint_index
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
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save/'
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
      iks =  global_kpoint_index (nkstot, 1  )
      ike = iks + nks - 1
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
         CALL qexml_write_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp, &
                                  monopole, zmon, relaxz, block, block_1, block_2,&
                                  block_height ) 
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
                         TETRA_TYPE = tetra_type, TETRA = tetra, TFIXED_OCC = tfixed_occ, LSDA = lsda, &
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
      IF ( lrho ) CALL write_scf( rho, nspin )
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
#if defined __HDF5
          USE hdf5_qe,   ONLY :  prepare_for_writing_final, write_gkhdf5, &
                                 gk_hdf5_write, h5fclose_f
          USE io_files,  ONLY : tmp_dir
#endif

          IMPLICIT NONE
          !
          INTEGER,            INTENT(IN) :: iun, ik
          CHARACTER(LEN=256), INTENT(IN) :: filename
          !
          INTEGER, ALLOCATABLE :: igwk(:,:)
          INTEGER, ALLOCATABLE :: itmp(:)
          CHARACTER(LEN = 256) :: filename_hdf5
          INTEGER              :: ierr 
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
#if defined __HDF5
             filename_hdf5=trim(tmp_dir) //"gk.hdf5"
             CALL prepare_for_writing_final(gk_hdf5_write,inter_pool_comm,filename_hdf5,ik)
             CALL write_gkhdf5(gk_hdf5_write,xk(:,ik),igwk(1:ngk_g(ik),ik), &
                              mill_g(1:3,igwk(1:ngk_g(ik),ik)),ik)
             CALL h5fclose_f(gk_hdf5_write%file_id, ierr)
#else

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
#endif
          END IF
          !
          DEALLOCATE( igwk )
          !
        END SUBROUTINE write_gk
        !
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
      USE io_rho_xml,    ONLY : read_scf
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
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save/'
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
         CALL read_scf( rho, nspin )
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
      USE ktetra,           ONLY : ntetra, tetra_type
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
         CALL qexml_read_occ( NTETRA = ntetra, TETRA_TYPE = tetra_type, IERR=ierr )
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
    !------------------------------------------------------------------------
    SUBROUTINE read_efield( ierr )
      !----------------------------------------------------------------------
      !
      USE extfield, ONLY : tefield, dipfield, edir, emaxpos, eopreg, eamp, & !TB
                           monopole, zmon, relaxz, block, block_1, block_2, block_height
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
                                 MONOPOLE=monopole, ZMON=zmon, RELAXZ=relaxz, & !TB
                                 BLOCK=block, BLOCK_1=block_1, BLOCK_2=block_2, &
                                 BLOCK_HEIGHT=block_height, FOUND=found, IERR=ierr )
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
         monopole = .FALSE.
         !
      END IF
      !
      CALL mp_bcast( tefield,  ionode_id, intra_image_comm )
      CALL mp_bcast( dipfield, ionode_id, intra_image_comm )
      CALL mp_bcast( edir,     ionode_id, intra_image_comm )
      CALL mp_bcast( emaxpos,  ionode_id, intra_image_comm )
      CALL mp_bcast( eopreg,   ionode_id, intra_image_comm )
      CALL mp_bcast( eamp,     ionode_id, intra_image_comm )
      CALL mp_bcast( monopole, ionode_id, intra_image_comm )
      CALL mp_bcast( zmon,     ionode_id, intra_image_comm )
      CALL mp_bcast( relaxz,   ionode_id, intra_image_comm )
      CALL mp_bcast( block,    ionode_id, intra_image_comm )
      CALL mp_bcast( block_1,  ionode_id, intra_image_comm )
      CALL mp_bcast( block_2,  ionode_id, intra_image_comm )
      CALL mp_bcast( block_height, ionode_id, intra_image_comm )
      !
      lefield_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_efield
    !
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

      IF ( inlc > 0 ) THEN
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
    !------------------------------------------------------------------------
    SUBROUTINE read_occupations( ierr )
      !------------------------------------------------------------------------
      !
      USE lsda_mod,       ONLY : lsda, nspin
      USE fixed_occ,      ONLY : tfixed_occ, f_inp
      USE ktetra,         ONLY : ntetra, tetra, tetra_type
      USE klist,          ONLY : ltetra, lgauss, ngauss, degauss, smearing
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
                               LTETRA=ltetra, NTETRA=ntetra, TETRA=tetra, TETRA_TYPE= tetra_type, TFIXED_OCC=tfixed_occ, &
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
         CALL mp_bcast( tetra_type, ionode_id, intra_image_comm )
         if(tetra_type == 0) CALL mp_bcast( tetra,  ionode_id, intra_image_comm )
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
#if defined __HDF5
      USE hdf5_qe,              ONLY : evc_hdf5_write, read_attributes_hdf5, &
                                       prepare_for_reading_final
      USE mp_pools,             ONLY : inter_pool_comm
#endif

      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256)   :: filename
      INTEGER              :: ik, ipol, ik_eff, num_k_points
      INTEGER, ALLOCATABLE :: kisort(:)
      INTEGER              :: nkl, nkr, npwx_g
      INTEGER              :: nupdwn(2), ike, iks, npw_g, ispin
      INTEGER, EXTERNAL    :: global_kpoint_index
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
      iks = global_kpoint_index (nkstot, 1)
      ike = iks + nks - 1
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
    !------------------------------------------------------------------------
    FUNCTION pp_check_file()
      !------------------------------------------------------------------------
      !
      USE io_global,         ONLY : ionode, ionode_id
      USE mp_images,         ONLY : intra_image_comm
      USE control_flags,     ONLY : lkpoint_dir, tqr, tq_smoothing, tbeta_smoothing
      !
      IMPLICIT NONE
      !
      LOGICAL            :: pp_check_file
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: ierr
      LOGICAL            :: lval, found, back_compat
      !
      !
      dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save/'
      filename = TRIM( dirname ) // TRIM( xmlpun )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = filename, IERR = ierr )
      !
      CALL mp_bcast ( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'pp_check_file', 'file ' // &
                 & TRIM( dirname ) // ' not found', ierr )
      !
      ! set a flag for back compatibility (before fmt v1.4.0)
      !
      back_compat = .FALSE.
      !
      IF ( TRIM( version_compare( qexml_version, "1.4.0" )) == "older") &
         back_compat = .TRUE.
      !
      IF ( ionode ) THEN
         !
         IF ( .NOT. back_compat ) THEN
             !
             CALL iotk_scan_begin( iunpun, "CONTROL" ) 
             !
         ENDIF
         !
         CALL iotk_scan_dat( iunpun, "PP_CHECK_FLAG", lval, FOUND = found)
         !
         IF ( .NOT. found ) lval = .FALSE. 
         !
         CALL iotk_scan_dat( iunpun, "LKPOINT_DIR", lkpoint_dir, FOUND = found)
         !
         IF ( .NOT. found ) lkpoint_dir = .TRUE. 
         !
         CALL iotk_scan_dat( iunpun, "Q_REAL_SPACE", tqr, FOUND = found)
         !
         IF ( .NOT. found ) tqr = .FALSE. 
         !
         CALL iotk_scan_dat( iunpun, "TQ_SMOOTHING", tq_smoothing, FOUND = found)
         !
         IF ( .NOT. found ) tq_smoothing = .FALSE. 
         !
         CALL iotk_scan_dat( iunpun, "TBETA_SMOOTHING", tbeta_smoothing, FOUND = found)
         !
         IF ( .NOT. found ) tbeta_smoothing = .FALSE. 
         !
         IF ( .NOT. back_compat ) THEN
             !
             CALL iotk_scan_end( iunpun, "CONTROL" ) 
             !
         ENDIF
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( lval, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( lkpoint_dir, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( tqr, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( tq_smoothing, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( tbeta_smoothing, ionode_id, intra_image_comm )
      !
      pp_check_file = lval
      !
      RETURN
      !
    END FUNCTION pp_check_file
    !
#endif
    !
END MODULE pw_restart

