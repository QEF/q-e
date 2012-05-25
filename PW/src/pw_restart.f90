!
! Copyright (C) 2005-2008 Quantum ESPRESSO group
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
  ! ... written by Carlo Sbraccia (2005)
  !
  USE iotk_module
  USE xml_io_base
  USE xml_io_base, ONLY : default_fmt_version => fmt_version  
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun, delete_if_present, &
                        qexml_version, qexml_version_init, pseudo_dir
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global, ONLY : my_pool_id, intra_image_comm, intra_pool_comm, &
                        my_bgrp_id, intra_image_comm, intra_bgrp_comm, mpime
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
  ! variables to describe qexml current version
  ! and back compatibility
  !
  LOGICAL :: qexml_version_before_1_4_0 = .FALSE.
  !
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE pw_writefile( what )
      !------------------------------------------------------------------------
      !
      USE control_flags,        ONLY : istep, twfcollect, conv_ions, &
                                       lscf, lkpoint_dir, gamma_only, &
                                       tqr, noinv, do_makov_payne 
      USE realus,               ONLY : real_space
      USE global_version,       ONLY : version_number
      USE cell_base,            ONLY : at, bg, alat, tpiba, tpiba2, &
                                       ibrav, celldm
      USE gvect,   ONLY : ig_l2g
      USE ions_base,            ONLY : nsp, ityp, atm, nat, tau, if_pos
      USE noncollin_module,     ONLY : noncolin, npol
      USE io_files,             ONLY : nwordwfc, iunwfc, iunigk, psfile
      USE buffers,              ONLY : get_buffer
      USE wavefunctions_module, ONLY : evc
      USE klist,                ONLY : nks, nkstot, xk, ngk, wk, qnorm, &
                                       lgauss, ngauss, degauss, nelec, &
                                       two_fermi_energies, nelup, neldw
      USE start_k,             ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       nks_start, xk_start, wk_start
      USE ktetra,               ONLY : ntetra, tetra, ltetra
      USE gvect,                ONLY : ngm, ngm_g, g, mill
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE gvecs,              ONLY : ngms_g, dual
      USE fft_base,           ONLY : dffts
      USE wvfct,                ONLY : npw, npwx, g2kin, et, wg, &
                                       igk, nbnd, ecutwfc
      USE ener,                 ONLY : ef, ef_up, ef_dw
      USE fixed_occ,            ONLY : tfixed_occ, f_inp
      USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_lmax,  &
                                       Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_alpha
      USE spin_orb,             ONLY : lspinorb, domag
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
      USE mp_global,            ONLY : kunit, nproc, nproc_pool, me_pool, &
                                       nproc_image, nproc_bgrp, me_bgrp, &
                                       root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm, &
                                       root_bgrp, intra_bgrp_comm, inter_bgrp_comm, nbgrp
      USE funct,                ONLY : get_exx_fraction, dft_is_hybrid, &
                                       get_screening_parameter, exx_is_active
      USE exx,                  ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut
      USE cellmd,               ONLY : lmovecell, cell_factor 
      !
      USE martyna_tuckerman, ONLY: do_comp_mt
      USE esm,               ONLY: do_comp_esm
      
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      !
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname, filename
      INTEGER               :: i, ig, ik, ngg, ierr, ipol, ik_eff, num_k_points
      INTEGER               :: npool, nkbl, nkl, nkr, npwx_g
      INTEGER               :: ike, iks, npw_g, ispin, inlc
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      INTEGER,  ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:), mill_g(:,:)
      LOGICAL               :: lwfc
      REAL(DP), ALLOCATABLE :: raux(:)
      !
      !
      lwfc  = .FALSE.
      !
      SELECT CASE( what )
      CASE( "all" )
         !
         lwfc  = twfcollect
         !
      CASE DEFAULT
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
            CALL create_directory( kpoint_dir( dirname, i ) )
            !
         END DO
         !
      END IF
      !
      IF ( nkstot > 0 ) THEN
         !
         ! ... find out the number of pools
         !
         npool = nproc / nproc_pool
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
      !
      igk_l2g = 0
      !
      IF ( nks > 1 ) REWIND( iunigk )
      !
      DO ik = 1, nks
         !
         npw = ngk (ik)
         IF ( nks > 1 ) READ( iunigk ) igk
         !
         CALL gk_l2gmap( ngm, ig_l2g(1), npw, igk(1), igk_l2g(1,ik) )
         !
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
         CALL iotk_open_write( iunpun, FILE = TRIM( dirname ) // '/' // &
                             & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
         IF (.NOT.(lkpoint_dir)) &
            CALL iotk_open_write( iunout, FILE = TRIM( dirname ) // '/' // &
                    & TRIM( xmlpun )//'.eig', BINARY = .FALSE., IERR = ierr )
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
         CALL write_header( "PWSCF", TRIM(version_number) )
         !
!-------------------------------------------------------------------------------
! ... CONTROL 
!-------------------------------------------------------------------------------
         !
         CALL write_control( PP_CHECK_FLAG=conv_ions, LKPOINT_DIR=lkpoint_dir, &
                             Q_REAL_SPACE=tqr, BETA_REAL_SPACE=real_space )
         !
!-------------------------------------------------------------------------------
! ... CELL
!-------------------------------------------------------------------------------
         !
         CALL write_cell( ibrav, celldm, alat, &
                          at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3), &
                          do_makov_payne, do_comp_mt, do_comp_esm )
         IF (lmovecell) CALL write_moving_cell(lmovecell, cell_factor)
         !
!-------------------------------------------------------------------------------
! ... IONS
!-------------------------------------------------------------------------------
         !
         CALL write_ions( nsp, nat, atm, ityp, psfile, &
                          pseudo_dir, amass, tau, if_pos, dirname, alat )
         !
!-------------------------------------------------------------------------------
! ... SYMMETRIES
!-------------------------------------------------------------------------------
         !
         CALL write_symmetry( ibrav, nrot, nsym, invsym, noinv, &
                              time_reversal, no_t_rev, ft, s, sname, irt,  &
                              nat, t_rev )
         !
!-------------------------------------------------------------------------------
! ... ELECTRIC FIELD
!-------------------------------------------------------------------------------
         !
         CALL write_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp) 
         !
!
!-------------------------------------------------------------------------------
! ... PLANE_WAVES
!-------------------------------------------------------------------------------
         !
         CALL write_planewaves( ecutwfc, dual, npwx_g, gamma_only, dfftp%nr1, dfftp%nr2, &
                                dfftp%nr3, ngm_g, dffts%nr1, dffts%nr2, dffts%nr3, ngms_g, dfftp%nr1, &
                                dfftp%nr2, dfftp%nr3, mill_g, lwfc )
         !
!-------------------------------------------------------------------------------
! ... SPIN
!-------------------------------------------------------------------------------
         !
         CALL write_spin( lsda, noncolin, npol, lspinorb, domag )
         !
         CALL write_magnetization(starting_magnetization, angle1, angle2, nsp, &
                                  two_fermi_energies, i_cons, mcons, bfield, &
                                  ef_up/e2, ef_dw/e2, nelup, neldw, lambda)
         !
!-------------------------------------------------------------------------------
! ... EXCHANGE_CORRELATION
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         inlc = get_inlc()
         !
         CALL write_xc( DFT = dft_name, NSP = nsp, LDA_PLUS_U = lda_plus_u,                  &
                        LDA_PLUS_U_KIND = lda_plus_u_kind, HUBBARD_LMAX = Hubbard_lmax,      &
                        HUBBARD_L = Hubbard_l, HUBBARD_U = Hubbard_U, HUBBARD_J = Hubbard_J, &
                        HUBBARD_ALPHA = Hubbard_alpha, INLC = inlc, VDW_TABLE_NAME = vdw_table_name, &
                        PSEUDO_DIR = pseudo_dir, DIRNAME = dirname)
         IF ( dft_is_hybrid() ) CALL write_exx &
                         ( x_gamma_extrapolation, nq1, nq2, nq3, &
                         exxdiv_treatment, yukawa, ecutvcut, &
                         get_exx_fraction(), &
                         get_screening_parameter(), exx_is_active() )
         !
!-------------------------------------------------------------------------------
! ... OCCUPATIONS
!-------------------------------------------------------------------------------
         !
         CALL write_occ( LGAUSS = lgauss, NGAUSS = ngauss, &
                         DEGAUSS = degauss, LTETRA = ltetra, NTETRA = ntetra, &
                         TETRA = tetra, TFIXED_OCC = tfixed_occ, LSDA = lsda, &
                         NSTATES_UP = nbnd, NSTATES_DOWN = nbnd, F_INP = f_inp )
         !
!-------------------------------------------------------------------------------
! ... BRILLOUIN_ZONE
!-------------------------------------------------------------------------------
         !
         CALL write_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3, &
                        qnorm, nks_start, xk_start, wk_start )
         !
!-------------------------------------------------------------------------------
! ... PARALLELISM
!-------------------------------------------------------------------------------
         !
         CALL write_para( kunit, nproc, nproc_pool, nproc_image ) 
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
         CALL iotk_write_begin( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_write_dat  ( iunpun, "NUMBER_OF_K-POINTS", num_k_points )
         !
         CALL iotk_write_dat  ( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin )
         !
         CALL iotk_write_dat  ( iunpun, "NON-COLINEAR_CALCULATION", noncolin )
         !
         CALL iotk_write_dat  ( iunpun, "NUMBER_OF_ATOMIC_WFC", natomwfc )
         !
         IF ( nspin == 2 ) THEN
            ! Compatibility with CP
            CALL iotk_write_attr( attr, "UP", nbnd, FIRST = .TRUE. )
            CALL iotk_write_attr( attr, "DW", nbnd )
            CALL iotk_write_dat( iunpun, &
                                 "NUMBER_OF_BANDS", nbnd, ATTR = attr )
            CALL iotk_write_attr( attr, "UP", NINT(nelup), FIRST = .TRUE. )
            CALL iotk_write_attr( attr, "DW", NINT(neldw) )
            CALL iotk_write_dat( iunpun, &
                                 "NUMBER_OF_ELECTRONS", nelec, ATTR = attr )
         ELSE
         !
            CALL iotk_write_dat  ( iunpun, "NUMBER_OF_BANDS", nbnd )
            CALL iotk_write_dat  ( iunpun, "NUMBER_OF_ELECTRONS", nelec )
         END IF
         !
         CALL iotk_write_attr ( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
         CALL iotk_write_empty( iunpun, "UNITS_FOR_K-POINTS", ATTR = attr )
         !
         CALL iotk_write_attr ( attr, "UNITS", "Hartree", FIRST = .TRUE. )
         CALL iotk_write_empty( iunpun, "UNITS_FOR_ENERGIES", ATTR = attr )
         !
         ! Fermi energy units in Hartree
         !
         IF (two_fermi_energies) THEN
             !
             CALL iotk_write_dat(iunpun,"TWO_FERMI_ENERGIES",two_fermi_energies)
             CALL iotk_write_dat( iunpun, "ELECTRONS_UP", nelup )
             CALL iotk_write_dat( iunpun, "ELECTRONS_DOWN", neldw )
             CALL iotk_write_dat( iunpun, "FERMI_ENERGY_UP", ef_up / e2 )
             CALL iotk_write_dat( iunpun, "FERMI_ENERGY_DOWN", ef_dw / e2 )
             !
         ELSE
             !
             CALL iotk_write_dat( iunpun, "FERMI_ENERGY", ef / e2)
             !
         ENDIF
         !
         CALL iotk_write_end  ( iunpun, "BAND_STRUCTURE_INFO" )
         !
         !
         CALL iotk_write_begin( iunpun, "EIGENVALUES" )
         !
      END IF
      !
      ALLOCATE( raux( nbnd) )
      !
!-------------------------------------------------------------------------------
! ... EIGENVALUES
!-------------------------------------------------------------------------------
      !
      k_points_loop1: DO ik = 1, num_k_points
         !
         IF ( ionode ) THEN
            !
            CALL iotk_write_begin( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            CALL iotk_write_dat( iunpun, "K-POINT_COORDS", xk(:,ik), COLUMNS=3 )
            !            
            CALL iotk_write_dat( iunpun, "WEIGHT", wk(ik) )
            !
            !
            IF ( nspin == 2 ) THEN
               !
               ispin = 1
               !
               IF (lkpoint_dir) THEN
                  filename=wfc_filename(".",'eigenval1', ik, EXTENSION='xml',&
                                      DIR=lkpoint_dir )
                  !
                  CALL iotk_link( iunpun, "DATAFILE.1", &
                               filename, CREATE = .FALSE., BINARY = .FALSE. )
               ELSE
                  CALL iotk_write_begin( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_UP" )
               ENDIF
               !
               IF ( wk(ik) == 0.D0 ) THEN
                  !
                  raux = wg(:,ik)
                  !
               ELSE
                  !
                  raux = wg(:,ik) / wk(ik)
                  !
               END IF
               !
               IF (lkpoint_dir) THEN
                  filename = wfc_filename( dirname, 'eigenval1', ik, &
                                   EXTENSION='xml',  DIR=lkpoint_dir )
                  !
                   CALL write_eig( iunout, filename, nbnd, et(:, ik) / e2, &
                     "Hartree", OCC = raux(:), IK=ik, ISPIN=ispin )
               ELSE
                   filename=' '
                   CALL write_eig( iunout, filename, nbnd, et(:, ik) / e2, &
                     "Hartree", OCC = raux(:), IK=ik, ISPIN=ispin,  &
                                LKPOINT_DIR=.FALSE. )
               ENDIF
               !
               ispin = 2
               !
               ik_eff = ik + num_k_points
               !
               IF (lkpoint_dir) THEN
                  filename = wfc_filename( ".", 'eigenval2', ik, &
                          EXTENSION='xml',  DIR=lkpoint_dir )
                  !
                  CALL iotk_link( iunpun, "DATAFILE.2", &
                               filename, CREATE = .FALSE., BINARY = .FALSE. )
               ELSE
                  CALL iotk_write_end( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_UP" )
                  CALL iotk_write_begin( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_DW" )
               ENDIF
               !
               IF ( wk(ik_eff) == 0.D0 ) THEN
                  !
                  raux = wg(:,ik_eff)
                  !
               ELSE
                  !
                  raux = wg(:,ik_eff) / wk(ik_eff)
                  !
               END IF
               !
               IF (lkpoint_dir) THEN
                  filename = wfc_filename( dirname, 'eigenval2', ik, &
                             EXTENSION = 'xml',  DIR=lkpoint_dir )
                  !
                  CALL write_eig( iunout, filename, nbnd, et(:, ik_eff) / e2, &
                               "Hartree", OCC = raux(:), IK = ik, ISPIN = ispin)
               ELSE
                  filename=' '
                  CALL write_eig( iunout, filename, nbnd, et(:, ik_eff) / e2, &
                               "Hartree", OCC = raux(:), IK = ik, &
                                         ISPIN = ispin, LKPOINT_DIR=.false.)
                  CALL iotk_write_end( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_DW" )

               ENDIF
               !
            ELSE
               !
               IF (lkpoint_dir) THEN
                  filename = wfc_filename( ".", 'eigenval', ik, &
                          EXTENSION='xml',  DIR=lkpoint_dir )
                  !
                  CALL iotk_link( iunpun, "DATAFILE", &
                               filename, CREATE = .FALSE., BINARY = .FALSE. )
               ELSE
                  CALL iotk_write_begin( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) ) )
               ENDIF
               !
               IF ( wk(ik) == 0.D0 ) THEN
                  !
                  raux(:) = wg(:,ik)
                  !
               ELSE
                  !
                  raux(:) = wg(:,ik) / wk(ik)
                  !
               END IF
               !
               IF (lkpoint_dir) THEN
                  filename = wfc_filename( dirname, 'eigenval', ik, &
                             EXTENSION='xml',  DIR=lkpoint_dir )
                  !
                  CALL write_eig( iunout, filename, nbnd, et(:, ik) / e2, &
                               "Hartree", OCC = raux(:), IK = ik )
               ELSE
                  filename=' '
                  CALL write_eig( iunout, filename, nbnd, et(:, ik) / e2, &
                               "Hartree", OCC = raux(:), IK = ik, &
                                 LKPOINT_DIR=.false. )
                  CALL iotk_write_end( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) ) )
               ENDIF
                  !
            END IF
            !
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
         END IF
         !
      ENDDO k_points_loop1
      !
      !
      IF (.NOT.lkpoint_dir.AND.ionode) CALL iotk_close_write( iunout )
      !
      DEALLOCATE ( raux )
      ! 
      ! 
      IF ( ionode ) THEN 
         !
         CALL iotk_write_end( iunpun, "EIGENVALUES" )
         !
         CALL iotk_write_begin( iunpun, "EIGENVECTORS" )
         !
         CALL iotk_write_dat  ( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwx_g )
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... EIGENVECTORS
!-------------------------------------------------------------------------------
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
               filename = wfc_filename( ".", 'gkvectors', ik, DIR=lkpoint_dir )
               !
               CALL iotk_link( iunpun, "GK-VECTORS", &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               filename = wfc_filename( dirname, 'gkvectors', ik, &
                                         DIR=lkpoint_dir )
            END IF
            !
         END IF
         !
         IF ( lwfc ) THEN
           !
           CALL write_gk( iunout, ik, filename )
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
         CALL iotk_close_write( iunpun )
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
      ! ... do not overwrite the scf charge density with a non-scf one
      ! ... also writes rho%ns if lda+U and rho%bec if PAW
      !
      IF ( lscf ) CALL write_rho( rho, nspin )
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( mill_g )
      DEALLOCATE( ngk_g )
      !
      CALL save_history( dirname, istep )
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
          !
          ! ... wavefunctions
          !
          IF ( nspin == 2 ) THEN
             !
             ! ... beware: with pools, this is correct only on ionode
             !
             ispin = isk(ik)
             !
             IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
                !
                CALL get_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
                !
             END IF
             !
             IF ( ionode ) THEN
                !
                filename = wfc_filename( ".", 'evc', ik, ispin, &
                                                     DIR=lkpoint_dir )
                !
                CALL iotk_link( iunpun, "WFC" // TRIM( iotk_index (ispin) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                !
                filename = wfc_filename( dirname, 'evc', ik, ispin, & 
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
             IF ( ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
                !
                CALL get_buffer ( evc, nwordwfc, iunwfc, (ik_eff-iks+1) )
                !
             END IF
             !
             IF ( ionode ) THEN
                !
                filename = wfc_filename( ".", 'evc', ik, ispin, &
                                                     DIR=lkpoint_dir )
                !
                CALL iotk_link( iunpun, "WFC"//TRIM( iotk_index( ispin ) ), &
                                filename, CREATE = .FALSE., BINARY = .TRUE. )
                !
                filename = wfc_filename( dirname, 'evc', ik, ispin, &
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
             IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
                !
                CALL get_buffer( evc, nwordwfc, iunwfc, (ik-iks+1) )
                !
             END IF
             !
             IF ( noncolin ) THEN
                !
                DO ipol = 1, npol
                   !
                   IF ( ionode ) THEN
                      !
                      filename = wfc_filename( ".", 'evc', ik, ipol, &
                                              DIR=lkpoint_dir )
                      !
                      CALL iotk_link(iunpun,"WFC"//TRIM(iotk_index(ipol)), &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
                      !
                      filename = wfc_filename( dirname, 'evc', ik, ipol, &
                             DIR=lkpoint_dir)
                      !
                   END IF
                   !
                   !!! TEMP
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
                   filename = wfc_filename( ".", 'evc', ik, DIR=lkpoint_dir )
                   !
                   CALL iotk_link( iunpun, "WFC", filename, &
                                   CREATE = .FALSE., BINARY = .TRUE. )
                   !
                   filename = wfc_filename( dirname, 'evc', ik, &
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
      USE mp_global,     ONLY : intra_pool_comm, intra_bgrp_comm
      USE mp,            ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: what
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256) :: dirname
      LOGICAL            :: lcell, lpw, lions, lspin, linit_mag, &
                            lxc, locc, lbz, lbs, lwfc, lheader,          &
                            lsymm, lrho, lefield
      !
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
      !
      SELECT CASE( what )
      CASE( 'header' )
         !
         lheader = .TRUE.
         !
      CASE( 'dim' )
         !
         CALL read_dim( dirname, ierr )
         !
         lbz = .TRUE.
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
      CASE( 'rho' )
         !
         lrho  = .TRUE.
         !
      CASE( 'wave' )
         !
         lpw   = .TRUE.
         lwfc  = .TRUE.
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
         CALL read_ef( dirname, ierr )
         RETURN
      CASE( 'exx' )
         CALL read_exx(dirname, ierr)
         RETURN
         !
      END SELECT
      !
      IF ( .NOT. lheader .AND. .NOT. qexml_version_init) &
         CALL errore( 'pw_readfile', 'qexml version not set', 71 )
      !
      !
      IF ( lheader ) THEN 
         !
         CALL read_header( dirname, ierr )
         !
         ! to be as safe as possible
         !
         IF ( ierr /= 0 ) THEN
             !
             qexml_version = TRIM( default_fmt_version )
             qexml_version_init = .TRUE.
             !
         ENDIF
         !
         ierr = 0 
         !
      ENDIF
      !
      IF ( lcell ) THEN
         !
         CALL read_cell( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lpw ) THEN
         !
         CALL read_planewaves( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lions ) THEN
         !
         CALL read_ions( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lspin ) THEN
         !
         CALL read_spin( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF (linit_mag) THEN
         !
         CALL read_magnetization( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      ENDIF
      IF ( lxc ) THEN
         !
         CALL read_xc( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( locc ) THEN
         !
         CALL read_occupations( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lbz ) THEN
         !
         CALL read_brillouin_zone( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lbs ) THEN
         !
         CALL read_band_structure( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lwfc ) THEN
         !
         CALL read_wavefunctions( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lsymm ) THEN
         !
         CALL read_symmetry( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lefield ) THEN
         !
         CALL read_efield( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      IF ( lrho ) THEN
         !
         ! ... to read the charge-density we use the routine from io_rho_xml 
         ! ... it also reads ns for ldaU and becsum for PAW
         CALL read_rho( rho, nspin )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE pw_readfile
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_header( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the format version of the current xml datafile
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr

      ierr = 0
      IF ( qexml_version_init ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr /=0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "HEADER" )
         !
         CALL iotk_scan_empty( iunpun, "FORMAT", ATTR=attr )
         !
         CALL iotk_scan_attr( attr, "VERSION", qexml_version )
         !
         qexml_version_init = .TRUE.
         !
         CALL iotk_scan_end( iunpun, "HEADER" )
         !
         !
         CALL iotk_close_read( iunpun )
         !
      ENDIF
      !
      CALL mp_bcast( qexml_version,       ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init,  ionode_id, intra_image_comm )
      
      !
      ! init logical variables for versioning
      !
      qexml_version_before_1_4_0 = .FALSE.
      !
      IF ( TRIM( version_compare( qexml_version, "1.4.0" )) == "older" ) &
         qexml_version_before_1_4_0 = .TRUE.
      !
    END SUBROUTINE read_header
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_dim( dirname, ierr )
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
      USE wvfct,            ONLY : nbnd, npwx, ecutwfc
      USE control_flags,    ONLY : gamma_only
      USE mp_global,        ONLY : kunit, nproc_file, nproc_pool_file, &
                                   nproc_image_file
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER  :: npwx_
      LOGICAL  :: found, found2
      !
      !
      ! ... first the entire CELL section is read
      !
      CALL read_cell( dirname, ierr )
      !
      IF ( ierr > 0 ) RETURN
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
      ! ... then selected tags are read from the other sections
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "IONS" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMS", nat )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPECIES", nsp )         
         !
         CALL iotk_scan_end( iunpun, "IONS" )
         !
         CALL iotk_scan_begin( iunpun, "SYMMETRIES", FOUND = found )
         !
         IF ( .NOT. found ) THEN
            !
            nsym = 1
            !
         ELSE
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_SYMMETRIES", nsym )
            !
            CALL iotk_scan_end( iunpun, "SYMMETRIES" )
            !
         END IF
         !
         CALL iotk_scan_begin( iunpun, "PLANE_WAVES" )
         !
         CALL iotk_scan_dat( iunpun, "WFC_CUTOFF", ecutwfc )
         !
         CALL iotk_scan_dat( iunpun, "RHO_CUTOFF", ecutrho )
         !
         ecutwfc = ecutwfc * e2
         ecutrho = ecutrho * e2
         !
         dual = ecutrho / ecutwfc
         !
         CALL iotk_scan_dat( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwx_ )
         !
         CALL iotk_scan_dat( iunpun, "GAMMA_ONLY", gamma_only )
         !
         CALL iotk_scan_empty( iunpun, "FFT_GRID", attr )
         CALL iotk_scan_attr( attr, "nr1", dfftp%nr1 )
         CALL iotk_scan_attr( attr, "nr2", dfftp%nr2 )
         CALL iotk_scan_attr( attr, "nr3", dfftp%nr3 )
         !
         CALL iotk_scan_dat( iunpun, "GVECT_NUMBER", ngm_g )
         !
         CALL iotk_scan_empty( iunpun, "SMOOTH_FFT_GRID", attr )
         CALL iotk_scan_attr( attr, "nr1s", dffts%nr1 )
         CALL iotk_scan_attr( attr, "nr2s", dffts%nr2 )
         CALL iotk_scan_attr( attr, "nr3s", dffts%nr3 )
         !
         CALL iotk_scan_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_g )
         !
         CALL iotk_scan_end( iunpun, "PLANE_WAVES" )
         !
         CALL iotk_scan_begin( iunpun, "SPIN" )
         !
         CALL iotk_scan_dat( iunpun, "LSDA", lsda )
         !
         CALL iotk_scan_dat( iunpun, "NON-COLINEAR_CALCULATION", &
              noncolin, FOUND = found )
         !
         IF ( .NOT. found ) noncolin = .FALSE.
         !
         CALL iotk_scan_end( iunpun, "SPIN" )
         !
         CALL iotk_scan_begin( iunpun, "OCCUPATIONS" )
         !
         CALL iotk_scan_dat( iunpun, &
                             "NUMBER_OF_TETRAHEDRA", ntetra, DEFAULT = 0 )
         !
         CALL iotk_scan_end( iunpun, "OCCUPATIONS" )
         !
         CALL iotk_scan_begin( iunpun, "BRILLOUIN_ZONE" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_K-POINTS", nkstot )
         !
         IF ( lsda ) nkstot = nkstot * 2
         !
         CALL iotk_scan_end( iunpun, "BRILLOUIN_ZONE" )
         !
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd )
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_scan_begin( iunpun, "PARALLELISM", FOUND = found )
         !
         IF ( .NOT. found ) THEN
            !
            kunit = 1
            nproc_file=1
            nproc_pool_file=1
            !
         ELSE
            !
            CALL iotk_scan_dat( iunpun, &
                 "GRANULARITY_OF_K-POINTS_DISTRIBUTION", kunit )
            !
            CALL iotk_scan_dat( iunpun, &
                              "NUMBER_OF_PROCESSORS", nproc_file, &
                               FOUND=found2 )
            IF ( .NOT. found2) nproc_file=1 ! compatibility
            !
            CALL iotk_scan_dat( iunpun, &
                            "NUMBER_OF_PROCESSORS_PER_POOL", nproc_pool_file,&
                               FOUND=found2 )
            IF ( .NOT. found2) nproc_pool_file=1 ! compatibility

            CALL iotk_scan_dat( iunpun, &
                            "NUMBER_OF_PROCESSORS_PER_IMAGE", nproc_image_file,&
                               FOUND=found2 )
            IF ( .NOT. found2) nproc_image_file=1 ! compatibility
            !
            CALL iotk_scan_end( iunpun, "PARALLELISM" )
            !
         END IF
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( nat,        ionode_id, intra_image_comm )
      CALL mp_bcast( nsp,        ionode_id, intra_image_comm )
      CALL mp_bcast( nsym,       ionode_id, intra_image_comm )
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
      CALL mp_bcast( lsda,       ionode_id, intra_image_comm )
      CALL mp_bcast( noncolin,   ionode_id, intra_image_comm )
      CALL mp_bcast( ntetra,     ionode_id, intra_image_comm )
      CALL mp_bcast( nkstot,     ionode_id, intra_image_comm )
      CALL mp_bcast( nelec,      ionode_id, intra_image_comm )
      CALL mp_bcast( nbnd,       ionode_id, intra_image_comm )
      CALL mp_bcast( kunit,      ionode_id, intra_image_comm )
      CALL mp_bcast( nproc_file, ionode_id, intra_image_comm )
      CALL mp_bcast( nproc_pool_file, ionode_id, intra_image_comm )
      CALL mp_bcast( nproc_image_file, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_dim
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_cell( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi
      USE run_info,  ONLY: title
      USE cell_base, ONLY : ibrav, alat, at, bg, celldm
      USE cell_base, ONLY : tpiba, tpiba2, omega
      USE cellmd,    ONLY : lmovecell, cell_factor
      USE control_flags, ONLY : do_makov_payne
      USE martyna_tuckerman, ONLY: do_comp_mt
      USE esm,           ONLY: do_comp_esm
      
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      LOGICAL :: found
      CHARACTER(LEN=80) :: bravais_lattice, es_corr
      !
      !
      ierr = 0
      IF ( lcell_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "CELL" )
         !
         CALL iotk_scan_dat( iunpun, "NON-PERIODIC_CELL_CORRECTION", &
                             es_corr, FOUND=found )
         IF ( .NOT. found ) es_corr="None"
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
         CALL iotk_scan_dat( iunpun, &
                             "BRAVAIS_LATTICE", bravais_lattice )
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
         CALL iotk_scan_dat( iunpun, "LATTICE_PARAMETER", alat )
         !
         ! ... some internal variables
         !
         tpiba  = 2.D0 * pi / alat
         tpiba2 = tpiba**2 
         !
         CALL iotk_scan_dat( iunpun, "CELL_DIMENSIONS", celldm )
         !
         CALL iotk_scan_begin( iunpun, "DIRECT_LATTICE_VECTORS" )
         CALL iotk_scan_dat(   iunpun, "a1", at(:,1) )
         CALL iotk_scan_dat(   iunpun, "a2", at(:,2) )
         CALL iotk_scan_dat(   iunpun, "a3", at(:,3) )
         CALL iotk_scan_end(   iunpun, "DIRECT_LATTICE_VECTORS" )
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
         CALL iotk_scan_end( iunpun, "CELL" )
         !
         CALL iotk_scan_begin( iunpun, "MOVING_CELL", found=lmovecell )
         IF (lmovecell) THEN
            CALL iotk_scan_dat( iunpun, "CELL_FACTOR", cell_factor)
            CALL iotk_scan_end( iunpun, "MOVING_CELL"  )
         END IF
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
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
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      ! this is where PP files should be read from
      !
      pseudo_dir_cur = trimcheck ( dirname ) 
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "IONS" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMS", nat )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
         !
         DO i = 1, nsp
            !
            IF ( qexml_version_before_1_4_0 ) THEN
               !
               CALL iotk_scan_dat( iunpun, "ATOM_TYPE", atm(i) )
               CALL iotk_scan_dat( iunpun, TRIM( atm(i) ) // "_MASS", amass(i) )
               CALL iotk_scan_dat( iunpun, &
                                   "PSEUDO_FOR_" // TRIM( atm(i) ), psfile(i) )
               !
            ELSE
               !
               ! current version
               !
               CALL iotk_scan_begin( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
               !
               CALL iotk_scan_dat( iunpun, "ATOM_TYPE", atm(i) )
               CALL iotk_scan_dat( iunpun, "MASS", amass(i) )
               CALL iotk_scan_dat( iunpun, "PSEUDO", psfile(i) )
               !
               CALL iotk_scan_end( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
               !
            ENDIF
            !
         ENDDO
         !
         ! this is the original location of PP files
         !
         CALL iotk_scan_dat( iunpun, "PSEUDO_DIR", pseudo_dir )
         !
      ENDIF
      !
      IF ( ionode ) THEN
         !
         DO i = 1, nat
            !
            CALL iotk_scan_empty( iunpun, &
                                  "ATOM" // TRIM( iotk_index(i) ), attr )
            !
            CALL iotk_scan_attr( attr, "INDEX",  ityp(i) )
            CALL iotk_scan_attr( attr, "tau",    tau(:,i) )
            CALL iotk_scan_attr( attr, "if_pos", if_pos(:,i) )
            !
            tau(:,i) = tau(:,i) / alat
            !
         END DO
         !
         CALL iotk_scan_end( iunpun, "IONS" )
         !
         CALL iotk_close_read( iunpun )
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
    SUBROUTINE read_symmetry( dirname, ierr )
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
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER  :: i, nat_
      LOGICAL  :: found
      !
      ierr = 0
      IF ( lsymm_read ) RETURN
      !
      IF ( .NOT. lpw_read ) &
         CALL errore( 'read_symmetry', 'read planewaves first', 1 )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "SYMMETRIES", FOUND = found )
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
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_SYMMETRIES", nsym )
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_BRAVAIS_SYMMETRIES", &
                             nrot, FOUND = found )
            IF (.NOT. found) nrot = nsym
            !
            CALL iotk_scan_dat( iunpun, "INVERSION_SYMMETRY", invsym )
            CALL iotk_scan_dat( iunpun, "DO_NOT_USE_TIME_REVERSAL", &
                             noinv, FOUND = found )
            IF (.NOT. found) noinv = .FALSE.
            CALL iotk_scan_dat( iunpun, "TIME_REVERSAL_FLAG", &
                             time_reversal, FOUND = found )
            IF (.NOT. found) time_reversal = .TRUE.
            CALL iotk_scan_dat( iunpun, "NO_TIME_REV_OPERATIONS", &
                             no_t_rev, FOUND = found )
            IF (.NOT. found) no_t_rev = .FALSE.
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMS", nat_ )
            !
            DO i = 1, nsym
               !
               CALL iotk_scan_begin( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
               !
               CALL iotk_scan_empty( iunpun, "INFO", ATTR = attr )
               CALL iotk_scan_attr( attr, "NAME",  sname(i) )
               CALL iotk_scan_attr( attr, "T_REV", t_rev(i) )
               !
               CALL iotk_scan_dat( iunpun, "ROTATION", s(:,:,i) )
               CALL iotk_scan_dat( iunpun, "FRACTIONAL_TRANSLATION", ft(:,i) )
               CALL iotk_scan_dat( iunpun, "EQUIVALENT_IONS", irt(i,1:nat_) )
               !
               ftau(1,i) = NINT( ft(1,i)*DBLE( dfftp%nr1 ) )
               ftau(2,i) = NINT( ft(2,i)*DBLE( dfftp%nr2 ) )
               ftau(3,i) = NINT( ft(3,i)*DBLE( dfftp%nr3 ) )
               !
               CALL iotk_scan_end( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
               !
            END DO
            !
            ! indices of inverse operations and matrices in cartesian axis
            ! are not saved to disk (maybe they should), are recalculated here 
            !
            CALL inverse_s ()
            CALL s_axis_to_cart ()
            !
            DO i = nsym+1, nrot
               !
               CALL iotk_scan_begin( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
               CALL iotk_scan_empty( iunpun, "INFO", ATTR = attr )
               CALL iotk_scan_attr( attr, "NAME",  sname(i) )
               CALL iotk_scan_dat( iunpun, "ROTATION", s(:,:,i) )
               CALL iotk_scan_end( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
               !
            END DO
            !
            CALL iotk_scan_end( iunpun, "SYMMETRIES" )
            !
         END IF
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( nsym,   ionode_id, intra_image_comm )
      CALL mp_bcast( nrot,   ionode_id, intra_image_comm )
      CALL mp_bcast( invsym, ionode_id, intra_image_comm )
      CALL mp_bcast( noinv,  ionode_id, intra_image_comm )
      CALL mp_bcast( time_reversal,  ionode_id, intra_image_comm )
      CALL mp_bcast( no_t_rev,  ionode_id, intra_image_comm )
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
    SUBROUTINE read_efield( dirname, ierr )
      !----------------------------------------------------------------------
      !
      USE extfield, ONLY : tefield, dipfield, edir, emaxpos, eopreg, eamp
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      LOGICAL                       :: found
      !
      ierr = 0
      IF ( lefield_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "ELECTRIC_FIELD", FOUND = found )
         !
         IF ( found ) THEN
            !
            CALL iotk_scan_dat( iunpun, "HAS_ELECTRIC_FIELD", tefield )
            !
            CALL iotk_scan_dat( iunpun, "HAS_DIPOLE_CORRECTION", dipfield )
            !
            CALL iotk_scan_dat( iunpun, "FIELD_DIRECTION", edir )
            !
            CALL iotk_scan_dat( iunpun, "MAXIMUM_POSITION", emaxpos )
            !
            CALL iotk_scan_dat( iunpun, "INVERSE_REGION", eopreg )
            !
            CALL iotk_scan_dat( iunpun, "FIELD_AMPLITUDE", eamp )
            !
            CALL iotk_scan_end( iunpun, "ELECTRIC_FIELD" )
            !
         ELSE
            !
            tefield  = .FALSE.
            dipfield = .FALSE.
            !
         END IF
         !
         CALL iotk_close_read( iunpun )
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
    !------------------------------------------------------------------------
    SUBROUTINE read_planewaves( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE gvect,           ONLY : ngm_g, ecutrho
      USE gvecs,           ONLY : ngms_g, dual
      USE fft_base,        ONLY : dfftp
      USE fft_base,        ONLY : dffts
      USE wvfct,           ONLY : npwx, g2kin, ecutwfc
      USE control_flags,   ONLY : gamma_only
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER  :: npwx_
      !
      ierr = 0
      IF ( lpw_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "PLANE_WAVES" )
         !
         CALL iotk_scan_dat( iunpun, "WFC_CUTOFF", ecutwfc )
         !
         CALL iotk_scan_dat( iunpun, "RHO_CUTOFF", ecutrho )
         !
         ecutwfc = ecutwfc * e2
         ecutrho = ecutrho * e2
         !
         dual = ecutrho / ecutwfc
         !
         CALL iotk_scan_dat( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwx_ )
         !
         CALL iotk_scan_dat( iunpun, "GAMMA_ONLY", gamma_only )
         !
         CALL iotk_scan_empty( iunpun, "FFT_GRID", attr )
         CALL iotk_scan_attr( attr, "nr1", dfftp%nr1 )
         CALL iotk_scan_attr( attr, "nr2", dfftp%nr2 )
         CALL iotk_scan_attr( attr, "nr3", dfftp%nr3 )
         !
         CALL iotk_scan_dat( iunpun, "GVECT_NUMBER", ngm_g )
         !
         CALL iotk_scan_empty( iunpun, "SMOOTH_FFT_GRID", attr )
         CALL iotk_scan_attr( attr, "nr1s", dffts%nr1 )
         CALL iotk_scan_attr( attr, "nr2s", dffts%nr2 )
         CALL iotk_scan_attr( attr, "nr3s", dffts%nr3 )
         !
         CALL iotk_scan_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_g )
         !
         CALL iotk_scan_end( iunpun, "PLANE_WAVES" )
         !
         CALL iotk_close_read( iunpun )
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
    SUBROUTINE read_spin( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE spin_orb,         ONLY : lspinorb, domag
      USE lsda_mod,         ONLY : nspin, lsda
      USE noncollin_module, ONLY : noncolin, npol
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      LOGICAL :: found
      !
      ierr = 0
      IF ( lspin_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "SPIN" )
         !
         CALL iotk_scan_dat( iunpun, "LSDA", lsda )
         !
         CALL iotk_scan_dat( iunpun, "NON-COLINEAR_CALCULATION", &
              noncolin, FOUND = found )
         IF ( .not. found ) noncolin = .FALSE.
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
         IF ( noncolin ) THEN
            !
            CALL iotk_scan_dat( iunpun, "SPINOR_DIM", npol )
            !
         ELSE
            !
            npol = 1
            !
         END IF
         !
         CALL iotk_scan_dat( iunpun, "SPIN-ORBIT_CALCULATION", &
              lspinorb, FOUND = found )
         IF ( .NOT. found ) lspinorb = .FALSE.
         !
         CALL iotk_scan_dat( iunpun, "SPIN-ORBIT_DOMAG", domag, &
               FOUND = found )
         IF ( .NOT. found ) domag = .FALSE.
         !
         CALL iotk_scan_end( iunpun, "SPIN" )
         !
         CALL iotk_close_read( iunpun )
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
    SUBROUTINE read_magnetization( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE constants,        ONLY : PI
      USE klist,            ONLY : two_fermi_energies, nelup, neldw
      USE ener,             ONLY : ef_up, ef_dw
      USE lsda_mod,         ONLY : starting_magnetization
      USE noncollin_module, ONLY : angle1, angle2, i_cons, mcons, bfield, &
                                   lambda
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      LOGICAL :: found
      INTEGER :: i, nsp
      !
      ierr = 0
      IF ( lstarting_mag_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "MAGNETIZATION_INIT", FOUND = found )
         !
         IF( found ) THEN
            !
            CALL iotk_scan_dat(iunpun,"CONSTRAINT_MAG", i_cons)

            CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
            !
            DO i=1,nsp
               !
               CALL iotk_scan_begin( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
               !
               CALL iotk_scan_dat( iunpun, "STARTING_MAGNETIZATION", &
                         starting_magnetization(i) )
               CALL iotk_scan_dat( iunpun, "ANGLE1", angle1(i) )
               CALL iotk_scan_dat( iunpun, "ANGLE2", angle2(i) )
               !
               angle1(i)=angle1(i)*PI/180.d0
               angle2(i)=angle2(i)*PI/180.d0
               !
               IF (i_cons==1.OR.i_cons==2) THEN
                  !
                  CALL iotk_scan_dat( iunpun, "CONSTRANT_1", mcons(1,i) )
                  CALL iotk_scan_dat( iunpun, "CONSTRANT_2", mcons(2,i) )
                  CALL iotk_scan_dat( iunpun, "CONSTRANT_3", mcons(3,i) )
                  !
               ENDIF
               !
               CALL iotk_scan_end( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
               !
            ENDDO

            IF (i_cons==3) THEN
               !
               CALL iotk_scan_dat( iunpun, "FIXED_MAGNETIZATION_1", mcons(1,1) )
               CALL iotk_scan_dat( iunpun, "FIXED_MAGNETIZATION_2", mcons(2,1) )
               CALL iotk_scan_dat( iunpun, "FIXED_MAGNETIZATION_3", mcons(3,1) )
               !
            ELSE IF (i_cons==4) THEN
               ! 
               CALL iotk_scan_dat( iunpun, "MAGNETIC_FIELD_1", bfield(1) )
               CALL iotk_scan_dat( iunpun, "MAGNETIC_FIELD_2", bfield(2) )
               CALL iotk_scan_dat( iunpun, "MAGNETIC_FIELD_3", bfield(3) )
               !
            ENDIF
            !
            CALL iotk_scan_dat(iunpun,"TWO_FERMI_ENERGIES", &
                 two_fermi_energies, FOUND = found)
            IF ( .not. found ) two_fermi_energies=.FALSE.
            !
            IF (two_fermi_energies) THEN
                !
                CALL iotk_scan_dat( iunpun, "FIXED_MAGNETIZATION", mcons(3,1) )
                CALL iotk_scan_dat( iunpun, "ELECTRONS_UP", nelup )
                CALL iotk_scan_dat( iunpun, "ELECTRONS_DOWN", neldw )
                CALL iotk_scan_dat( iunpun, "FERMI_ENERGY_UP", ef_up )
                CALL iotk_scan_dat( iunpun, "FERMI_ENERGY_DOWN", ef_dw )
                !
                ef_up = ef_up * e2 
                ef_dw = ef_dw * e2 
                !
            ENDIF
            !
            IF (i_cons>0) CALL iotk_scan_dat(iunpun,"LAMBDA",lambda)
            !
            CALL iotk_scan_end( iunpun, "MAGNETIZATION_INIT" )
            ! 
         END IF
         !
         CALL iotk_close_read( iunpun )
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
    SUBROUTINE read_xc( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE ions_base, ONLY : nsp
      USE funct,     ONLY : enforce_input_dft
      USE ldaU,      ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_lmax, &
                            Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_alpha
      USE kernel_table, ONLY : vdw_table_name
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=20) :: dft_name
      INTEGER           :: nsp_, inlc
      LOGICAL           :: found, nomsg = .true.
      !
      ierr = 0
      IF ( lxc_read ) RETURN
      !
      IF ( .NOT. lions_read ) &
         CALL errore( 'read_xc', 'read ions first', 1 )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "EXCHANGE_CORRELATION" )
         !
         CALL iotk_scan_dat( iunpun, "DFT", dft_name )
         !
         CALL iotk_scan_dat( iunpun, "LDA_PLUS_U_CALCULATION", lda_plus_u, &
                             FOUND = found )
         IF ( .NOT. found ) lda_plus_u = .FALSE.
         !
         IF ( lda_plus_u ) THEN
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPECIES", nsp_ )
            !
            CALL iotk_scan_dat( iunpun, "LDA_PLUS_U_KIND", lda_plus_u_kind )
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_LMAX", Hubbard_lmax )
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_L", Hubbard_l(1:nsp_) )
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_U", Hubbard_U(1:nsp_) )
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_J", Hubbard_J(1:3,1:nsp_) )
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_ALPHA", Hubbard_alpha(1:nsp_) )
            !
         END IF

         !
         ! Vdw DF 
         !
         CALL iotk_scan_dat( iunpun, "NON_LOCAL_DF", inlc, FOUND = found )
         IF ( found ) THEN 
            !
            IF ( inlc == 1 .OR. inlc == 2 ) THEN
               !
               CALL iotk_scan_dat( iunpun, "VDW_KERNEL_NAME", vdw_table_name )
               !
            END IF
         ELSE
            inlc = 0
         ENDIF
         !
         CALL iotk_scan_end( iunpun, "EXCHANGE_CORRELATION" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( dft_name,   ionode_id, intra_image_comm )
      CALL mp_bcast( lda_plus_u, ionode_id, intra_image_comm )
      CALL mp_bcast( inlc, ionode_id, intra_image_comm )
      !
      IF ( lda_plus_u ) THEN
         !
         CALL mp_bcast( lda_plus_u_kind, ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_lmax,  ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_l ,    ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_U,     ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_J,     ionode_id, intra_image_comm )
         CALL mp_bcast( Hubbard_alpha, ionode_id, intra_image_comm )
         !
      END IF

      IF ( inlc == 1 .OR. inlc == 2 ) THEN
         CALL mp_bcast( vdw_table_name,  ionode_id, intra_image_comm )
      END IF

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
    SUBROUTINE read_brillouin_zone( dirname, ierr )
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
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: i, ik, num_k_points
      LOGICAL :: found
      !
      ierr = 0
      IF ( lbz_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "BRILLOUIN_ZONE" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_K-POINTS", num_k_points )
        
         !
         nkstot = num_k_points
         !
         IF ( lsda ) nkstot = num_k_points * 2
         !
         CALL iotk_scan_empty( iunpun, "MONKHORST_PACK_GRID", attr )
         CALL iotk_scan_attr( attr, "nk1", nk1 )
         CALL iotk_scan_attr( attr, "nk2", nk2 )
         CALL iotk_scan_attr( attr, "nk3", nk3 )   
         CALL iotk_scan_empty( iunpun, "MONKHORST_PACK_OFFSET", attr )
         CALL iotk_scan_attr( attr, "k1", k1 )
         CALL iotk_scan_attr( attr, "k2", k2 )
         CALL iotk_scan_attr( attr, "k3", k3 )   
         !
         DO ik = 1, num_k_points
            !
            CALL iotk_scan_empty( iunpun, "K-POINT" // &
                                & TRIM( iotk_index( ik ) ), attr )
            !
            CALL iotk_scan_attr( attr, "XYZ", xk(:,ik) )
            !
            CALL iotk_scan_attr( attr, "WEIGHT", wk(ik) )
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
         CALL iotk_scan_dat( iunpun, "STARTING_K-POINTS", nks_start, &
                             FOUND = found )
         IF (.NOT. found) nks_start=0
         IF (nks_start > 0 ) THEN
            IF (.NOT.ALLOCATED(xk_start)) ALLOCATE(xk_start(3,nks_start))
            IF (.NOT.ALLOCATED(wk_start)) ALLOCATE(wk_start(nks_start))
         END IF
         DO ik = 1, nks_start
            !
            CALL iotk_scan_empty( iunpun, "K-POINT_START" // &
                                & TRIM( iotk_index( ik ) ), attr )
            !
            CALL iotk_scan_attr( attr, "XYZ", xk_start(:,ik) )
            !
            CALL iotk_scan_attr( attr, "WEIGHT", wk_start(ik) )
            !
         END DO
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_BRAVAIS_SYMMETRIES", &
                             nrot, FOUND = found )
         IF (.NOT. found) THEN
            nrot=0
         ELSE IF (nrot > 0 .AND. nrot < 49 ) THEN
            DO i = 1, nrot
               CALL iotk_scan_begin( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
               !
               CALL iotk_scan_empty( iunpun, "INFO", ATTR = attr )
               CALL iotk_scan_attr( attr, "NAME",  sname(i) )
               CALL iotk_scan_dat( iunpun, "ROTATION", s(:,:,i) )
               CALL iotk_scan_end( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
               !
            END DO
         ELSE
            CALL errore ( 'pw_writefile ', & 
                'incorrect number of symmetries for lattice', nrot )
         END IF
         !
         CALL iotk_scan_dat( iunpun, "NORM-OF-Q", qnorm, FOUND = found )
         IF (.not. found) qnorm=0.0_DP

         CALL iotk_scan_end( iunpun, "BRILLOUIN_ZONE" )
         CALL iotk_close_read( iunpun )
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
    SUBROUTINE read_occupations( dirname, ierr )
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
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: i
      LOGICAL :: found
      !
      ierr = 0
      IF ( locc_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "OCCUPATIONS" )
         !
         CALL iotk_scan_dat( iunpun, "SMEARING_METHOD", lgauss, &
              FOUND = found )
         IF ( .NOT. found ) lgauss = .FALSE.
         !
         IF ( lgauss ) THEN
            !
            CALL iotk_scan_dat( iunpun, "SMEARING_TYPE", ngauss )
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
            CALL iotk_scan_dat( iunpun, "SMEARING_PARAMETER", degauss )
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
         CALL iotk_scan_dat( iunpun, "TETRAHEDRON_METHOD", ltetra, &
              FOUND = found )
         IF ( .NOT. found ) ltetra = .FALSE.
         !
         IF ( ltetra ) THEN
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_TETRAHEDRA", ntetra )
            !
            DO i = 1, ntetra
               !
               CALL iotk_scan_dat( iunpun, "TETRAHEDRON" // &
                                 & iotk_index( i ), tetra(1:4,i) )
               !
            END DO
            !
         ELSE
            !
            ntetra = 0
            !
         END IF
         !
         CALL iotk_scan_dat( iunpun, "FIXED_OCCUPATIONS", tfixed_occ, &
              FOUND = found )
         IF ( .NOT. found ) tfixed_occ = .FALSE.
         !
         IF ( tfixed_occ ) THEN
            !
            CALL iotk_scan_empty( iunpun, "INFO", ATTR=attr, FOUND=found )
            !
            IF ( .NOT. found ) THEN
               !
               nupdwn(1:2) = nbnd
               !
            ELSE
               !
               IF ( qexml_version_before_1_4_0 ) THEN
                  !
                  CALL iotk_scan_attr( attr, "nelup", nupdwn(1) )
                  CALL iotk_scan_attr( attr, "neldw", nupdwn(2) )
                  !
               ELSE
                  !
                  ! current version
                  !
                  CALL iotk_scan_attr( attr, "nstates_up", nupdwn(1) )
                  CALL iotk_scan_attr( attr, "nstates_down", nupdwn(2) )
                  !
               ENDIF
               !
            ENDIF
            !
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
            CALL iotk_scan_dat( iunpun, "INPUT_OCC_UP", f_inp(1:nupdwn(1),1) )
            !
            IF ( lsda ) THEN
               CALL iotk_scan_dat( iunpun, "INPUT_OCC_DOWN", f_inp(1:nupdwn(2),2) )
            ENDIF
            !
         END IF
         !
         CALL iotk_scan_end( iunpun, "OCCUPATIONS" )
         !
         CALL iotk_close_read( iunpun )
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
      INTEGER :: ik, ik_eff, num_k_points, iunout
      LOGICAL :: found, two_fermi_energies_
      !
      ierr = 0
      IF ( lbs_read ) RETURN
      !
      IF ( .NOT. lspin_read ) &
         CALL errore( 'read_band_structure', 'read spin first', 1 )
      IF ( .NOT. lbz_read ) &
         CALL errore( 'read_band_structure', 'read band_structure first', 1 )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF (.NOT.lkpoint_dir) THEN
         !
         IF ( ionode ) &
            CALL iotk_open_read( iunout, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun )//'.eig', IERR = ierr )
         !
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         !
         IF ( ierr > 0 ) RETURN
         !
      END IF
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMIC_WFC", natomwfc, &
              FOUND = found )
         IF ( .NOT. found ) natomwfc = 0
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd )
         !
         CALL iotk_scan_dat( iunpun, "FERMI_ENERGY", ef, FOUND = found )
         !
         IF ( found ) THEN
            ef = ef * e2
         ELSE
            ef = 0.d0
         END IF
         !
         CALL iotk_scan_dat( iunpun, "TWO_FERMI_ENERGIES", &
                 two_fermi_energies_, FOUND = found)
         IF ( .not. found ) two_fermi_energies_=.FALSE.
         !
         IF ( two_fermi_energies_ ) THEN
             !
             CALL iotk_scan_dat( iunpun, "FERMI_ENERGY_UP", ef_up )
             CALL iotk_scan_dat( iunpun, "FERMI_ENERGY_DOWN", ef_dw )
             !
             ef_up = ef_up * e2
             ef_dw = ef_dw * e2
             !
         ENDIF
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE_INFO" )
         !
      END IF
      !
      num_k_points = nkstot
      !
      IF ( lsda ) num_k_points = nkstot / 2
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "EIGENVALUES" )
         !
         k_points_loop: DO ik = 1, num_k_points
            !
            CALL iotk_scan_begin( iunpun, &
                                  "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            IF ( lsda ) THEN
               !
               isk(ik) = 1
               !
               IF (lkpoint_dir) THEN
                  CALL iotk_scan_begin(iunpun, "DATAFILE"//TRIM(iotk_index(1)) &
                                             , FOUND = found)
                  IF (.NOT. found ) GO TO 10 ! workaround: PW-CP compatibility
                  CALL iotk_scan_dat  ( iunpun, "EIGENVALUES", et(:,ik)  )
                  CALL iotk_scan_dat  ( iunpun, "OCCUPATIONS", wg(:,ik) )
                  CALL iotk_scan_end(iunpun, "DATAFILE"//TRIM(iotk_index(1)) )
               ELSE
                  CALL iotk_scan_begin( iunout, &
                  "DATA_EIG"//TRIM( iotk_index(ik) )//"_SPIN_UP", FOUND=found )
                  IF (.NOT. found ) GO TO 10 ! workaround: PW-CP compatibility
                  CALL iotk_scan_dat  ( iunout, "EIGENVALUES", et(:,ik)  )
                  CALL iotk_scan_dat  ( iunout, "OCCUPATIONS", wg(:,ik) )
                  CALL iotk_scan_end( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_UP")
               ENDIF
               !
  10           CONTINUE
               ik_eff = ik + num_k_points
               isk(ik_eff) = 2
               !
               IF (lkpoint_dir) THEN
                  CALL iotk_scan_begin(iunpun,"DATAFILE"//TRIM(iotk_index(2)) &
                                             , FOUND = found)
                  IF (.NOT. found ) GO TO 20 ! workaround: PW-CP compatibility
                  CALL iotk_scan_dat  ( iunpun, "EIGENVALUES", et(:,ik_eff) )
                  CALL iotk_scan_dat  ( iunpun, "OCCUPATIONS", wg(:,ik_eff) )
                  CALL iotk_scan_end( iunpun, "DATAFILE"//TRIM(iotk_index(2)) )
               ELSE
                  CALL iotk_scan_begin( iunout, &
                  "DATA_EIG"//TRIM( iotk_index(ik) )//"_SPIN_DW", FOUND=found )
                  IF (.NOT. found ) GO TO 20 ! workaround: PW-CP compatibility
                  CALL iotk_scan_dat  ( iunout, "EIGENVALUES", et(:,ik_eff) )
                  CALL iotk_scan_dat  ( iunout, "OCCUPATIONS", wg(:,ik_eff) )
                  CALL iotk_scan_end( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_DW")
               ENDIF
  20           CONTINUE
               !
            ELSE
               !
               isk(ik) = 1
               !
               IF (lkpoint_dir) THEN
                  CALL iotk_scan_begin( iunpun, "DATAFILE" , FOUND = found)
                  IF (.NOT. found ) GO TO 15 ! workaround: PW-CP compatibility
                  CALL iotk_scan_dat  ( iunpun, "EIGENVALUES", et(:,ik) )
                  CALL iotk_scan_dat  ( iunpun, "OCCUPATIONS", wg(:,ik) )
                  CALL iotk_scan_end  ( iunpun, "DATAFILE" )
               ELSE
                  CALL iotk_scan_begin( iunout, &
                  "DATA_EIG"//TRIM( iotk_index(ik) ), FOUND = found )
                  IF (.NOT. found ) GO TO 15 ! workaround: PW-CP compatibility
                  CALL iotk_scan_dat  ( iunout, "EIGENVALUES", et(:,ik) )
                  CALL iotk_scan_dat  ( iunout, "OCCUPATIONS", wg(:,ik) )
                  CALL iotk_scan_end( iunout, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) ))
               ENDIF
  15           CONTINUE
               !
            END IF
            !
            CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
         END DO k_points_loop
         !
         et(:,:) = et(:,:) * e2
         !
         FORALL( ik = 1:nkstot ) wg(:,ik) = wg(:,ik)*wk(ik)
         !
         CALL iotk_scan_end( iunpun, "EIGENVALUES" )
         !
         CALL iotk_close_read( iunpun )
         !
         IF (.NOT.lkpoint_dir) CALL iotk_close_read( iunout )
         !
      END IF
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
      USE wvfct,                ONLY : npw, npwx, g2kin, et, wg, nbnd, ecutwfc
      USE wavefunctions_module, ONLY : evc
      USE io_files,             ONLY : nwordwfc, iunwfc
      USE buffers,              ONLY : save_buffer
      USE gvect,                ONLY : ngm, ngm_g, g, ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol
      USE mp_global,            ONLY : kunit, nproc, nproc_pool, me_pool, me_bgrp, nbgrp, &
                                       root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm, &
                                       root_bgrp, intra_bgrp_comm, inter_bgrp_comm
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
      !
      !
      IF ( iunwfc > 0 ) THEN
         !
         INQUIRE( UNIT = iunwfc, OPENED = opnd )
         !
         IF ( .NOT. opnd ) CALL errore( 'read_wavefunctions', &
                    & 'wavefunctions unit (iunwfc) is not opened', 1 )
      END IF
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
      IF ( nkstot > 0 ) THEN
         !
         ! ... find out the number of pools
         !
         npool = nproc / nproc_pool
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
      ALLOCATE( kisort( npwx ) )
      !
      DO ik = 1, nks
         !
         kisort = 0
         npw    = npwx
         !
         CALL gk_sort( xk(1,ik+iks-1), ngm, g, &
                       ecutwfc/tpiba2, npw, kisort(1), g2kin )
         !
         CALL gk_l2gmap( ngm, ig_l2g(1), npw, kisort(1), igk_l2g(1,ik) )
         !
         ngk(ik) = npw
         !
      END DO
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
               filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin, &
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
               filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin, &
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
                     filename = TRIM( wfc_filename( dirname, 'evc', ik, ipol, &
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
                  filename = TRIM( wfc_filename( dirname, 'evc', ik, &
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
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_wavefunctions
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_ef( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads only the Fermi energy
      !
      USE ener, ONLY : ef
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      LOGICAL :: found
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
      ! ... then selected tags are read from the other sections
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE_INFO" )
         !
         CALL iotk_scan_dat( iunpun, "FERMI_ENERGY", ef, FOUND = found )
         !
         IF (found) THEN
            ef = ef * e2
         ELSE
            ef = 0.d0
         END IF
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE_INFO" )
      END IF
      !
      CALL mp_bcast( ef, ionode_id, intra_image_comm )
      !
      IF ( ionode ) CALL iotk_close_read( iunpun )
      !
      RETURN
      !
    END SUBROUTINE read_ef
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_exx( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... read EXX variables
      !
      USE funct,                ONLY : set_exx_fraction, set_screening_parameter, &
                                       enforce_input_dft, start_exx
      USE exx,                  ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      REAL(DP) :: exx_fraction, screening_parameter
      LOGICAL :: exx_is_active, found
      !
      IF ( ionode ) THEN
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      END IF
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr > 0 ) RETURN
      IF ( ionode ) THEN
         CALL iotk_scan_begin( iunpun, "EXACT_EXCHANGE", FOUND = found )
      END IF
      CALL mp_bcast( found, ionode_id, intra_image_comm )
      IF ( ionode ) THEN
         IF ( found ) THEN
            CALL iotk_scan_dat(iunpun, "x_gamma_extrapolation", x_gamma_extrapolation)
            CALL iotk_scan_dat(iunpun, "nqx1", nq1)
            CALL iotk_scan_dat(iunpun, "nqx2", nq2)
            CALL iotk_scan_dat(iunpun, "nqx3", nq3)
            CALL iotk_scan_dat(iunpun, "exxdiv_treatment", exxdiv_treatment)
            CALL iotk_scan_dat(iunpun, "yukawa", yukawa)
            CALL iotk_scan_dat(iunpun, "ecutvcut", ecutvcut)
            CALL iotk_scan_dat(iunpun, "exx_fraction", exx_fraction)
            CALL iotk_scan_dat(iunpun, "screening_parameter", screening_parameter)
            CALL iotk_scan_dat(iunpun, "exx_is_active", exx_is_active)
            CALL iotk_scan_end( iunpun, "EXACT_EXCHANGE" )
         END IF
         CALL iotk_close_read( iunpun )
      END IF
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
      CALL mp_bcast( exx_is_active, ionode_id, intra_image_comm )
      !
      CALL set_exx_fraction(exx_fraction)
      CALL set_screening_parameter(screening_parameter)
      IF (exx_is_active) CALL start_exx( ) 
      !
      RETURN
      !
    END SUBROUTINE read_exx
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
         CALL errore( 'igk_l2g_kdip', 'unexpected dimension in ngg', 1 )
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
