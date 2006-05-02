!
! Copyright (C) 2005-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: pw_writefile, pw_readfile
  !
  INTEGER, PRIVATE :: iunout
  !
  LOGICAL :: lcell_read = .FALSE., &
             lpw_read   = .FALSE., &
             lions_read = .FALSE., &
             lspin_read = .FALSE., &
             lxc_read   = .FALSE., &
             locc_read  = .FALSE., &
             lbz_read   = .FALSE., &
             lbs_read   = .FALSE., &
             lwfc_read  = .FALSE., &
             lsymm_read = .FALSE.
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE pw_writefile( what )
      !------------------------------------------------------------------------
      !
      USE control_flags,        ONLY : istep, modenum, twfcollect
      USE cell_base,            ONLY : at, bg, alat, tpiba, tpiba2, &
                                       ibrav, symm_type, celldm
      USE reciprocal_vectors,   ONLY : ig_l2g
      USE ions_base,            ONLY : nsp, ityp, atm, nat, tau, if_pos
      USE noncollin_module,     ONLY : noncolin, npol
      USE io_files,             ONLY : nwordwfc, iunwfc, psfile, pseudo_dir
      USE wavefunctions_module, ONLY : evc, evc_nc
      USE klist,                ONLY : nks, nkstot, xk, ngk, wk, &
                                       lgauss, ngauss, degauss, nelec, xqq
      USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, ngm, ngm_g, &
                                       g, ig1, ig2, ig3, ecutwfc, dual
      USE basis,                ONLY : natomwfc
      USE gsmooth,              ONLY : nr1s, nr2s, nr3s, ngms_g
      USE ktetra,               ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       ntetra, tetra, ltetra
      USE wvfct,                ONLY : gamma_only, npw, npwx, g2kin, et, wg, &
                                       igk_l2g, nbnd
      USE ener,                 ONLY : ef
      USE fixed_occ,            ONLY : tfixed_occ, f_inp
      USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, &
                                       Hubbard_U, Hubbard_alpha
      USE spin_orb,             ONLY : lspinorb, domag
      USE symme,                ONLY : nsym, invsym, s, ftau, irt, t_rev
      USE char,                 ONLY : sname
      USE lsda_mod,             ONLY : nspin, isk, lsda
      USE ions_base,            ONLY : amass
      USE funct,                ONLY : get_dft_name
      USE scf,                  ONLY : rho
      USE sticks,               ONLY : dfftp
      USE mp_global,            ONLY : kunit, nproc, nproc_pool, me_pool
      USE mp_global,            ONLY : my_pool_id, intra_image_comm, &
                                       intra_pool_comm, inter_pool_comm
      USE mp,                   ONLY : mp_sum, mp_max
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      !
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname, filename, file_pseudo, rho_file_base
      CHARACTER(LEN=80)     :: bravais_lattice
      INTEGER               :: i, ig, ik, ngg, ierr, ipol, ik_eff, num_k_points
      !!!INTEGER,  ALLOCATABLE :: kisort(:)
      INTEGER               :: npool, nkbl, nkl, nkr, npwx_g
      INTEGER               :: ike, iks, npw_g, ispin
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      INTEGER,  ALLOCATABLE :: itmp(:,:)
      LOGICAL               :: lgvec, lwfc
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      !
      !
      lgvec = .FALSE.
      lwfc  = .FALSE.
      !
      SELECT CASE( what )
      CASE( "all" )
         !
         lgvec = .TRUE.
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
      DO i = 1, nkstot
         !
         CALL create_directory( kpoint_dir( dirname, i ) )
         !
      END DO
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
      CALL mp_sum( ngm_g, intra_pool_comm )
      !
      ! ... collect all G-vectors across processors within the pools
      !
      ALLOCATE( itmp( 3, ngm_g ) )
      !
      itmp = 0
      !
      DO ig = 1, ngm
         !
         itmp(1,ig_l2g(ig)) = ig1(ig)
         itmp(2,ig_l2g(ig)) = ig2(ig)
         itmp(3,ig_l2g(ig)) = ig3(ig)
         !
      END DO
      !
      CALL mp_sum( itmp, intra_pool_comm )
      !
      ! ... build the G+k array indexes - BEWARE:
      !!! The following operations are commented out because igk_l2g
      !!! is already calculated in the initialization (routine hinit0)
      !!! It should be done here, but there is a problem with variable
      !!! cell: the starting unit cell and G-vectors should be used
      !
      !!!ALLOCATE( kisort( npwx ) )
      !
      !!!DO ik = 1, nks
         !
         !!!kisort = 0
         !!!npw    = npwx
         !
         !!!CALL gk_sort( xk(1,ik+iks-1), ngm, g, &
         !!!              ecutwfc/tpiba2, npw, kisort(1), g2kin )
         !
         !!!CALL gk_l2gmap( ngm, ig_l2g(1), npw, kisort(1), igk_l2g(1,ik) )
         !
         !!!ngk(ik) = npw
         !
      !!!END DO
      !
      !!!DEALLOCATE( kisort )
      !
      ! ... compute the global number of G+k vectors for each k point
      !
      ALLOCATE( ngk_g( nkstot ) )
      !
      ngk_g = 0
      ngk_g(iks:ike) = ngk(1:nks)
      !
      CALL mp_sum( ngk_g )
      !
      ! ... compute the maximum G vector index among all G+k an processors
      !
      npw_g = MAXVAL( igk_l2g(:,:) )
      !
      CALL mp_max( npw_g )
      !
      ! ... compute the maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      IF ( ionode ) THEN
         !
         ! ... open XML descriptor
         !
         CALL iotk_open_write( iunpun, FILE = TRIM( dirname ) // '/' // &
                             & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      CALL errore( 'pw_writefile ', &
                   'cannot open restart file for writing', ierr )
      !
      IF ( ionode ) THEN  
         !
         ! ... here we start writing the punch-file
         !
!-------------------------------------------------------------------------------
! ... CELL
!-------------------------------------------------------------------------------
         !
         CALL write_cell( ibrav, symm_type, celldm, alat, &
                          at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3) )
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
         CALL write_symmetry( ibrav, symm_type, nsym, &
                          invsym, nr1, nr2, nr3, ftau, s, sname, irt, t_rev )
         !
!-------------------------------------------------------------------------------
! ... PLANE_WAVES
!-------------------------------------------------------------------------------
         !
         CALL write_planewaves( ecutwfc, dual, npwx, gamma_only, nr1, nr2, &
                                nr3, ngm_g, nr1s, nr2s, nr3s, ngms_g, nr1, &
                                nr2, nr3, itmp, lgvec )
         !
!-------------------------------------------------------------------------------
! ... SPIN
!-------------------------------------------------------------------------------
         !
         CALL write_spin( lsda, noncolin, npol, lspinorb, domag )
         !
!-------------------------------------------------------------------------------
! ... EXCHANGE_CORRELATION
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         !
         CALL write_xc( DFT = dft_name, NSP = nsp, LDA_PLUS_U = lda_plus_u, &
                        HUBBARD_LMAX = Hubbard_lmax, HUBBARD_L = Hubbard_l, &
                        HUBBARD_U = Hubbard_U, HUBBARD_ALPHA = Hubbard_alpha )
         !
!-------------------------------------------------------------------------------
! ... OCCUPATIONS
!-------------------------------------------------------------------------------
         !
         CALL write_occ( LGAUSS = lgauss, NGAUSS = ngauss, &
                         DEGAUSS = degauss, LTETRA = ltetra, NTETRA = ntetra, &
                         TETRA = tetra, TFIXED_OCC = tfixed_occ, LSDA = lsda, &
                         NELUP = nbnd, NELDW = nbnd, F_INP = f_inp )
      END IF
      !
      num_k_points = nkstot
      !
      IF ( nspin == 2 ) num_k_points = nkstot / 2
      !
      IF ( ionode ) THEN
         !
!-------------------------------------------------------------------------------
! ... BRILLOUIN_ZONE
!-------------------------------------------------------------------------------
         !
         CALL write_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3 )
         !
!-------------------------------------------------------------------------------
! ... PHONON
!-------------------------------------------------------------------------------         
         !
         CALL write_phonon( modenum, xqq )
         !
!-------------------------------------------------------------------------------
! ... PARALLELISM
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "PARALLELISM" )
         !
         CALL iotk_write_dat( iunpun, &
                              "GRANULARITY_OF_K-POINTS_DISTRIBUTION", kunit )
         !
         CALL iotk_write_end( iunpun, "PARALLELISM" )
         !
!-------------------------------------------------------------------------------
! ... CHARGE-DENSITY
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "CHARGE-DENSITY" )
         !
      END IF
      !
      rho_file_base = 'charge-density'
      !
      IF ( ionode ) &
         CALL iotk_link( iunpun, "RHO", rho_file_base, &
                         CREATE = .FALSE., BINARY = .TRUE. )
      !
      rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )      
      !
      IF ( nspin == 1 ) THEN
         !
         CALL write_rho_xml( rho_file_base, rho(:,1), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
      ELSE IF ( nspin == 2 ) THEN
         !
         ALLOCATE( rhoaux( SIZE( rho, 1 ) ) )
         !
         rhoaux = rho(:,1) + rho(:,2) 
         !
         CALL write_rho_xml( rho_file_base, rhoaux, &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         rho_file_base = 'spin-polarization'
         !
         IF ( ionode ) &
            CALL iotk_link( iunpun, "SPIN_POL", rho_file_base, &
                            CREATE = .FALSE., BINARY = .TRUE. )
         !
         rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
         !
         rhoaux = rho(:,1) - rho(:,2) 
         !
         CALL write_rho_xml( rho_file_base, rhoaux, &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         DEALLOCATE( rhoaux )
         !
      ELSE IF ( nspin == 4 ) THEN
         !
         CALL write_rho_xml( rho_file_base, rho(:,1), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         rho_file_base = 'magnetization.x'
         !
         IF ( ionode ) &
            CALL iotk_link( iunpun, "MAG_X", rho_file_base, &
                            CREATE = .FALSE., BINARY = .TRUE. )
         !
         rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
         !
         CALL write_rho_xml( rho_file_base, rho(:,2), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         rho_file_base = 'magnetization.y'
         !
         IF ( ionode ) &
            CALL iotk_link( iunpun, "MAG_Y", rho_file_base, &
                            CREATE = .FALSE., BINARY = .FALSE. )
         !
         rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
         !
         CALL write_rho_xml( rho_file_base, rho(:,3), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         !
         rho_file_base = 'magnetization.z'
         !
         IF ( ionode ) &
            CALL iotk_link( iunpun, "MAG_Z", rho_file_base, &
                            CREATE = .FALSE., BINARY = .TRUE. )
         !
         rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
         !
         CALL write_rho_xml( rho_file_base, rho(:,4), &
                             nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
      END IF
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_end( iunpun, "CHARGE-DENSITY" )
         !
!-------------------------------------------------------------------------------
! ... BAND_STRUCTURE
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_ATOMIC_WFC", natomwfc )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_BANDS", nbnd )
         !
         CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
         CALL iotk_write_dat( iunpun, "FERMI_ENERGY", ef / e2, ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin )
         !
         CALL iotk_write_begin( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
         !
      END IF
      !
      k_points_loop: DO ik = 1, num_k_points
         !
         IF ( ionode ) THEN
            !
            CALL iotk_write_begin( iunpun, &
                                   "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
            CALL iotk_write_dat( iunpun, &
                                 "K-POINT_COORDS", xk(:,ik), ATTR = attr )
            !            
            CALL iotk_write_dat( iunpun, "WEIGHT", wk(ik) )
            !
            IF ( nspin == 2 ) THEN
               !
               CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
               CALL iotk_write_dat( iunpun, "ET.1", et(:,ik) / e2, ATTR = attr )
               !
               IF ( wk(ik) == 0.D0 ) THEN
                  !
                  CALL iotk_write_dat( iunpun, "OCC.1", wg(:,ik) )
                  !
               ELSE
                  !
                  CALL iotk_write_dat( iunpun, "OCC.1", wg(:,ik) / wk(ik) )
                  !
               END IF
               !
               ik_eff = ik + num_k_points
               !
               CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
               CALL iotk_write_dat( iunpun, "ET.2", &
                                    et(:,ik_eff) / e2, ATTR = attr  )
               !
               IF ( wk(ik_eff) == 0.D0 ) THEN
                  !
                  CALL iotk_write_dat( iunpun, "OCC.2", wg(:,ik_eff) )
                  !
               ELSE
                  !
                  CALL iotk_write_dat( iunpun, "OCC.2", &
                                       wg(:,ik_eff) / wk(ik_eff) )
                  !
               END IF
               !
            ELSE
               !
               CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
               CALL iotk_write_dat( iunpun, "ET.1", et(:,ik) / e2, ATTR = attr )
               !
               IF ( wk(ik) == 0.D0 ) THEN
                  !
                  CALL iotk_write_dat( iunpun, "OCC.1", wg(:,ik) )
                  !
               ELSE
                  !
                  CALL iotk_write_dat( iunpun, "OCC.1", wg(:,ik) / wk(ik) )
                  !
               END IF
               !
            END IF
            !
            IF ( lgvec ) THEN
               !
               ! ... G+K vectors
               !
               filename = wfc_filename( ".", 'gkvectors', ik )
               !
               CALL iotk_link( iunpun, "GK-VECTORS", &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               filename = wfc_filename( dirname, 'gkvectors', ik )
               !
            END IF
            !
         END IF
         !
         IF ( lgvec ) CALL write_gk( iunout, ik, filename )
         !
         IF ( lwfc ) THEN
            !
            ! ... wavefunctions
            !
            IF ( nspin == 2 ) THEN
               !
               ispin = isk(ik)
               !
               IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
                  !
                  CALL davcio( evc, nwordwfc, iunwfc, (ik-iks+1), -1 )
                  !
               END IF
               !
               IF ( ionode ) THEN
                  !
                  filename = wfc_filename( ".", 'evc', ik, ispin )
                  !
                  CALL iotk_link( iunpun, "WFC", filename, &
                                  CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  filename = wfc_filename( dirname, 'evc', ik, ispin )
                  !
               END IF
               !
               CALL write_wfc( iunout, ik, nkstot, kunit, ispin, nspin, &
                               evc, npw_g, nbnd, igk_l2g(:,ik-iks+1),   &
                               ngk(ik-iks+1), filename, 1.D0 )
               !
               ik_eff = ik + num_k_points
               !
               ispin = isk(ik_eff)
               !
               IF ( ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
                  !
                  CALL davcio( evc, nwordwfc, iunwfc, (ik_eff-iks+1), -1 )
                  !
               END IF
               !
               IF ( ionode ) THEN
                  !
                  filename = wfc_filename( ".", 'evc', ik, ispin )
                  !
                  CALL iotk_link( iunpun, "WFC", filename, &
                                  CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  filename = wfc_filename( dirname, 'evc', ik, ispin )
                  !
               END IF
               !
               CALL write_wfc( iunout, ik_eff, nkstot, kunit, ispin, nspin, &
                               evc, npw_g, nbnd, igk_l2g(:,ik_eff-iks+1),   &
                               ngk(ik_eff-iks+1), filename, 1.D0 )
               !
            ELSE
               !
               IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
                  !
                  IF ( noncolin ) THEN
                     !
                     CALL davcio( evc_nc, nwordwfc, iunwfc, (ik-iks+1), -1 )
                     !
                  ELSE
                     !
                     CALL davcio( evc, nwordwfc, iunwfc, (ik-iks+1), -1 )
                     !
                  END IF
                  !
               END IF
               !
               IF ( noncolin ) THEN
                  !
                  DO ipol = 1, npol
                     !
                     IF ( ionode ) THEN
                        !
                        filename = wfc_filename( ".", 'evc', ik, ipol )
                        !
                        CALL iotk_link( iunpun, "WFC", filename, &
                                        CREATE = .FALSE., BINARY = .TRUE. )
                        !
                        filename = wfc_filename( dirname, 'evc', ik, ipol )
                        !
                     END IF
                     !
                     CALL write_wfc( iunout, ik, nkstot, kunit, ispin, nspin, &
                                     evc_nc(:,ipol,:), npw_g, nbnd,           &
                                     igk_l2g(:,ik-iks+1), ngk(ik-iks+1),      &
                                     filename, 1.D0 )
                     !
                  END DO
                  !
               ELSE
                  !
                  IF ( ionode ) THEN
                     !
                     filename = wfc_filename( ".", 'evc', ik )
                     !
                     CALL iotk_link( iunpun, "WFC", filename, &
                                     CREATE = .FALSE., BINARY = .TRUE. )
                     !
                     filename = wfc_filename( dirname, 'evc', ik )
                     !
                  END IF
                  !
                  CALL write_wfc( iunout, ik, nkstot, kunit, ispin, nspin, &
                                  evc, npw_g, nbnd, igk_l2g(:,ik-iks+1),   &
                                  ngk(ik-iks+1), filename, 1.D0 )
                  !
               END IF
               !
            END IF
            !
         END IF
         !
         IF ( ionode ) THEN
            !
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
            !
         END IF
         !
      END DO k_points_loop
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_end( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
         !
         CALL iotk_write_end( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_close_write( iunpun )
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( itmp )
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
          INTEGER, ALLOCATABLE :: itmp1(:)
          !
          !
          ALLOCATE( igwk( npwx_g, nkstot ) )
          !
          igwk(:,ik) = 0
          !
          ALLOCATE( itmp1( npw_g ) )
          !
          itmp1 = 0
          !
          IF ( ik >= iks .AND. ik <= ike ) THEN
             !
             DO ig = 1, ngk(ik-iks+1)
                !
                itmp1(igk_l2g(ig,ik-iks+1)) = igk_l2g(ig,ik-iks+1)
                !
             END DO
             !
          END IF
          !
          CALL mp_sum( itmp1 )
          !
          ngg = 0
          !
          DO ig = 1, npw_g
             !
             if ( itmp1(ig) == ig ) THEN
                !
                ngg = ngg + 1
                !
                igwk(ngg,ik) = ig
                !
             END IF
             !
          END DO
          !
          DEALLOCATE( itmp1 )
          !
          IF ( ionode ) THEN
             !
             CALL iotk_open_write( iun, &
                                   FILE = TRIM( filename ), BINARY = .TRUE. )
             !
             CALL iotk_write_begin( iun,"K-POINT" // iotk_index( ik ), attr )
             !
             CALL iotk_write_dat( iun, "INDEX", igwk(1:ngk_g(ik),ik) )
             CALL iotk_write_dat( iun, "GRID", itmp(1:3,igwk(1:ngk_g(ik),ik)), &
                                  COLUMNS = 3 )
             !
             CALL iotk_write_end( iun, "K-POINT" // iotk_index( ik ) )
             !
             CALL iotk_close_write( iun )
             !
          END IF
          !
          DEALLOCATE( igwk )
          !
        END SUBROUTINE write_gk
        !
    END SUBROUTINE pw_writefile
    !
    !------------------------------------------------------------------------
    SUBROUTINE pw_readfile( what, ierr )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: what
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256) :: dirname      
      LOGICAL            :: lexist, lcell, lpw, lions, lspin, &
                            lxc, locc, lbz, lbs, lwfc, &
                            lsymm, lph, lrho
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
      CALL errore( 'pw_readfile ', &
                   'no free units to read wavefunctions', ierr )      
      !
      lcell = .FALSE.
      lpw   = .FALSE.
      lions = .FALSE.
      lspin = .FALSE.
      lxc   = .FALSE.
      locc  = .FALSE.
      lbz   = .FALSE.
      lbs   = .FALSE.
      lwfc  = .FALSE.
      lsymm = .FALSE.
      lph   = .FALSE.
      lrho  = .FALSE.
      !
      SELECT CASE( what )
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
         lcell = .TRUE.
         lpw   = .TRUE.
         lions = .TRUE.
         lspin = .TRUE.
         lxc   = .TRUE.
         locc  = .TRUE.
         lbz   = .TRUE.
         lbs   = .TRUE.
         lsymm = .TRUE.
         lph   = .TRUE.
         !
      CASE( 'all' )
         !
         lcell = .TRUE.
         lpw   = .TRUE.
         lions = .TRUE.
         lspin = .TRUE.
         lxc   = .TRUE.
         locc  = .TRUE.
         lbz   = .TRUE.
         lbs   = .TRUE.
         lwfc  = .TRUE.
         lsymm = .TRUE.
         lph   = .TRUE.
         lrho  = .TRUE.
         !
      CASE( 'reset' )
         !
         lcell_read = .FALSE.
         lpw_read   = .FALSE.
         lions_read = .FALSE.
         lspin_read = .FALSE.
         lxc_read   = .FALSE.
         locc_read  = .FALSE.
         lbz_read   = .FALSE.
         lbs_read   = .FALSE.
         lwfc_read  = .FALSE.
         lsymm_read = .FALSE.
         !
      END SELECT
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
      IF ( lph ) THEN
         !
         CALL read_phonon( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      !
      IF ( lrho ) THEN
         !
         CALL read_rho( dirname, ierr )
         IF ( ierr > 0 ) RETURN
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE pw_readfile
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_dim( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine collects array dimensions from various sections
      ! ... plus with some other variables needed for array allocation 
      !
      USE ions_base,        ONLY : nat, nsp
      USE symme,            ONLY : nsym
      USE gvect,            ONLY : nr1, nr2, nr3, ngm_g, ecutwfc, dual
      USE gsmooth,          ONLY : nr1s, nr2s, nr3s, ngms_g
      USE lsda_mod,         ONLY : lsda
      USE noncollin_module, ONLY : noncolin
      USE ktetra,           ONLY : ntetra
      USE klist,            ONLY : nkstot, nelec
      USE wvfct,            ONLY : nbnd, npwx, gamma_only
      USE mp_global,        ONLY : kunit
      
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      REAL(DP) :: ecutrho
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
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id )
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
         CALL iotk_scan_begin( iunpun, "SYMMETRIES" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_SYMMETRIES", nsym )
         !
         CALL iotk_scan_end( iunpun, "SYMMETRIES" )
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
         CALL iotk_scan_dat( iunpun, "MAX_NPW", npwx )
         !
         CALL iotk_scan_dat( iunpun, "GAMMA_ONLY", gamma_only )
         !
         CALL iotk_scan_empty( iunpun, "FFT_GRID", attr )
         CALL iotk_scan_attr( attr, "nr1", nr1 )
         CALL iotk_scan_attr( attr, "nr2", nr2 )
         CALL iotk_scan_attr( attr, "nr3", nr3 )         
         !
         CALL iotk_scan_dat( iunpun, "GVECT_NUMBER", ngm_g )
         !
         CALL iotk_scan_empty( iunpun, "SMOOTH_FFT_GRID", attr )
         CALL iotk_scan_attr( attr, "nr1s", nr1s )
         CALL iotk_scan_attr( attr, "nr2s", nr2s )
         CALL iotk_scan_attr( attr, "nr3s", nr3s )
         !
         CALL iotk_scan_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_g )
         !
         CALL iotk_scan_end( iunpun, "PLANE_WAVES" )
         !
         CALL iotk_scan_begin( iunpun, "SPIN" )
         !
         CALL iotk_scan_dat( iunpun, "LSDA", lsda )
         !
         CALL iotk_scan_dat( iunpun, "NON-COLINEAR_CALCULATION", noncolin )
         !
         CALL iotk_scan_end( iunpun, "SPIN" )
         !
         CALL iotk_scan_begin( iunpun, "OCCUPATIONS" )
         !
         CALL iotk_scan_dat( iunpun, &
                             "NUMBER_OF_TETRAHEDRA", ntetra, DEFAULT = 1 )
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
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd )
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_scan_begin( iunpun, "PARALLELISM" )
         !
         CALL iotk_scan_dat( iunpun, &
                             "GRANULARITY_OF_K-POINTS_DISTRIBUTION", kunit )
         !
         CALL iotk_scan_end( iunpun, "PARALLELISM" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( nat,        ionode_id )
      CALL mp_bcast( nsp,        ionode_id )
      CALL mp_bcast( nsym,       ionode_id )
      CALL mp_bcast( ecutwfc,    ionode_id )
      CALL mp_bcast( dual,       ionode_id )
      CALL mp_bcast( npwx,       ionode_id )
      CALL mp_bcast( gamma_only, ionode_id )
      CALL mp_bcast( nr1,        ionode_id )
      CALL mp_bcast( nr2,        ionode_id )
      CALL mp_bcast( nr3,        ionode_id )
      CALL mp_bcast( ngm_g,      ionode_id )
      CALL mp_bcast( nr1s,       ionode_id )
      CALL mp_bcast( nr2s,       ionode_id )
      CALL mp_bcast( nr3s,       ionode_id )
      CALL mp_bcast( ngms_g,     ionode_id )
      CALL mp_bcast( lsda,       ionode_id )
      CALL mp_bcast( noncolin,   ionode_id )
      CALL mp_bcast( ntetra,     ionode_id )
      CALL mp_bcast( nkstot,     ionode_id )
      CALL mp_bcast( nelec,      ionode_id )
      CALL mp_bcast( nbnd,       ionode_id )
      CALL mp_bcast( kunit,      ionode_id )
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
      USE char,      ONLY : title, crystal
      USE cell_base, ONLY : ibrav, alat, symm_type, at, bg, celldm
      USE cell_base, ONLY : tpiba, tpiba2, omega
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=80) :: bravais_lattice
      !
      !
      IF ( lcell_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "CELL" )
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
         IF ( TRIM( bravais_lattice ) == "Trigonal R" .OR. &
              TRIM( bravais_lattice ) == "Hexagonal and Trigonal P" ) THEN
            !
            symm_type = 'hexagonal'
            !
         ELSE
            !
            symm_type = 'cubic'
            !
         END IF
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
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( ibrav,     ionode_id )
      CALL mp_bcast( symm_type, ionode_id )
      CALL mp_bcast( alat,      ionode_id )
      CALL mp_bcast( celldm,    ionode_id )
      CALL mp_bcast( tpiba,     ionode_id )
      CALL mp_bcast( tpiba2,    ionode_id )
      CALL mp_bcast( omega,     ionode_id )
      CALL mp_bcast( at,        ionode_id )
      CALL mp_bcast( bg,        ionode_id )
      !
      ! ... crystal is always set to empty string (see PW/input.f90)
      !
      crystal = ' '
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
      USE io_files,  ONLY : psfile, pseudo_dir
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: i
      !
      !
      IF ( lions_read ) RETURN
      !
      IF ( .NOT. lcell_read ) &
         CALL errore( 'read_ions', 'read cell first', 1 )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      IF ( ierr > 0 ) RETURN
      !
      pseudo_dir = TRIM( dirname )
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
            CALL iotk_scan_dat( iunpun, "ATOM_TYPE", atm(i) )
            CALL iotk_scan_dat( iunpun, TRIM( atm(i) ) // "_MASS", amass(i) )
            CALL iotk_scan_dat( iunpun, &
                                "PSEUDO_FOR_" // TRIM( atm(i) ), psfile(i) )
            !
         END DO
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
         END DO
         !
         tau(:,:) = tau(:,:) / alat
         !
         CALL iotk_scan_end( iunpun, "IONS" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( nat,    ionode_id )
      CALL mp_bcast( nsp,    ionode_id )
      CALL mp_bcast( atm,    ionode_id )
      CALL mp_bcast( amass,  ionode_id )
      CALL mp_bcast( psfile, ionode_id )
      CALL mp_bcast( ityp,   ionode_id )
      CALL mp_bcast( tau,    ionode_id )
      CALL mp_bcast( if_pos, ionode_id )
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
      USE symme, ONLY : nsym, invsym, s, ftau, irt, t_rev
      USE char,  ONLY : sname
      USE gvect, ONLY : nr1, nr2, nr3
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER  :: i
      REAL(DP) :: tmp(3)
      !
      !
      IF ( lsymm_read ) RETURN
      !
      IF ( .NOT. lpw_read ) &
         CALL errore( 'read_symmetry', 'read planewaves first', 1 )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "SYMMETRIES" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_SYMMETRIES", nsym )
         !
         CALL iotk_scan_dat( iunpun, "INVERSION_SYMMETRY", invsym )
         !
         DO i = 1, nsym
            !
            CALL iotk_scan_empty( iunpun, &
                                  "SYMM" // TRIM( iotk_index( i ) ), attr )
            !
            CALL iotk_scan_attr( attr, "ROT",        s(:,:,i) )
            CALL iotk_scan_attr( attr, "T_REV",      t_rev(i) )
            CALL iotk_scan_attr( attr, "FRAC_TRANS", tmp(:) )
            CALL iotk_scan_attr( attr, "NAME",       sname(i) )
            CALL iotk_scan_attr( attr, "EQ_IONS",    irt(i,:) )
            !
            ftau(1,i) = tmp(1) * DBLE( nr1 )
            ftau(2,i) = tmp(2) * DBLE( nr2 )
            ftau(3,i) = tmp(3) * DBLE( nr3 )
            !
         END DO         
         !
         CALL iotk_scan_end( iunpun, "SYMMETRIES" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( nsym,   ionode_id )
      CALL mp_bcast( invsym, ionode_id )
      CALL mp_bcast( s,      ionode_id )
      CALL mp_bcast( ftau,   ionode_id )
      CALL mp_bcast( sname,  ionode_id )
      CALL mp_bcast( irt,    ionode_id )
      !
      lsymm_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_symmetry
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_planewaves( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE gvect,   ONLY : nr1, nr2, nr3, ngm_g, ecutwfc, dual
      USE gsmooth, ONLY : nr1s, nr2s, nr3s, ngms_g
      USE wvfct,   ONLY : gamma_only, npwx, g2kin
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      REAL(DP) :: ecutrho
      !
      !
      IF ( lpw_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
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
         CALL iotk_scan_dat( iunpun, "MAX_NPW", npwx )
         !
         CALL iotk_scan_dat( iunpun, "GAMMA_ONLY", gamma_only )
         !
         CALL iotk_scan_empty( iunpun, "FFT_GRID", attr )
         CALL iotk_scan_attr( attr, "nr1", nr1 )
         CALL iotk_scan_attr( attr, "nr2", nr2 )
         CALL iotk_scan_attr( attr, "nr3", nr3 )         
         !
         CALL iotk_scan_dat( iunpun, "GVECT_NUMBER", ngm_g )
         !
         CALL iotk_scan_empty( iunpun, "SMOOTH_FFT_GRID", attr )
         CALL iotk_scan_attr( attr, "nr1s", nr1s )
         CALL iotk_scan_attr( attr, "nr2s", nr2s )
         CALL iotk_scan_attr( attr, "nr3s", nr3s )
         !
         CALL iotk_scan_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_g )
         !
         CALL iotk_scan_end( iunpun, "PLANE_WAVES" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( ecutwfc,    ionode_id )
      CALL mp_bcast( dual,       ionode_id )
      CALL mp_bcast( npwx,       ionode_id )
      CALL mp_bcast( gamma_only, ionode_id )
      CALL mp_bcast( nr1,        ionode_id )
      CALL mp_bcast( nr2,        ionode_id )
      CALL mp_bcast( nr3,        ionode_id )
      CALL mp_bcast( ngm_g,      ionode_id )
      CALL mp_bcast( nr1s,       ionode_id )
      CALL mp_bcast( nr2s,       ionode_id )
      CALL mp_bcast( nr3s,       ionode_id )
      CALL mp_bcast( ngms_g,     ionode_id )
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
      !
      IF ( lspin_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "SPIN" )
         !
         CALL iotk_scan_dat( iunpun, "LSDA", lsda )
         !
         CALL iotk_scan_dat( iunpun, "NON-COLINEAR_CALCULATION", noncolin )
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
         CALL iotk_scan_dat( iunpun, "SPIN-ORBIT_CALCULATION", lspinorb )
         CALL iotk_scan_dat( iunpun, "SPIN-ORBIT_DOMAG", domag )
         !
         CALL iotk_scan_end( iunpun, "SPIN" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( lsda,     ionode_id )
      CALL mp_bcast( nspin,    ionode_id )
      CALL mp_bcast( noncolin, ionode_id )
      CALL mp_bcast( npol,     ionode_id )
      CALL mp_bcast( lspinorb, ionode_id )
      !
      lspin_read = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE read_spin
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_xc( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE ions_base, ONLY : nsp
      USE funct,     ONLY : set_dft_from_name
      USE ldaU,      ONLY : lda_plus_u, Hubbard_lmax, &
                            Hubbard_l, Hubbard_U, Hubbard_alpha
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=20) :: dft_name
      !
      !
      IF ( lxc_read ) RETURN
      !
      IF ( .NOT. lions_read ) &
         CALL errore( 'read_xc', 'read ions first', 1 )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "EXCHANGE_CORRELATION" )
         !
         CALL iotk_scan_dat( iunpun, "DFT", dft_name )
         !
         CALL iotk_scan_dat( iunpun, "LDA_PLUS_U_CALCULATION", lda_plus_u )
         !
         IF ( lda_plus_u ) THEN
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_LMAX", Hubbard_lmax )
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_L", &
                                Hubbard_l(1:Hubbard_lmax) )
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_U", Hubbard_U(1:nsp) )
            !
            CALL iotk_scan_dat( iunpun, "HUBBARD_ALPHA", Hubbard_alpha(1:nsp) )
            !
         END IF         
         !
         CALL iotk_scan_end( iunpun, "EXCHANGE_CORRELATION" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( dft_name,   ionode_id )
      CALL mp_bcast( lda_plus_u, ionode_id )
      !
      IF ( lda_plus_u ) THEN
         !
         CALL mp_bcast( Hubbard_lmax,  ionode_id )
         CALL mp_bcast( Hubbard_l ,    ionode_id )
         CALL mp_bcast( Hubbard_U,     ionode_id )
         CALL mp_bcast( Hubbard_alpha, ionode_id )
         !
      END IF
      !
      CALL set_dft_from_name( dft_name )
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
      USE klist,    ONLY : nkstot, xk, wk
      USE ktetra,   ONLY : nk1, nk2, nk3, k1, k2, k3
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: ik, num_k_points
      !
      !
      IF ( lbz_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
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
         !
         CALL iotk_scan_end( iunpun, "BRILLOUIN_ZONE" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( nkstot, ionode_id )
      CALL mp_bcast( xk,     ionode_id )
      CALL mp_bcast( wk,     ionode_id )
      CALL mp_bcast( nk1, ionode_id )
      CALL mp_bcast( nk2, ionode_id )
      CALL mp_bcast( nk3, ionode_id )
      CALL mp_bcast( k1, ionode_id )
      CALL mp_bcast( k2, ionode_id )
      CALL mp_bcast( k3, ionode_id )
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
      USE lsda_mod,  ONLY : lsda, nspin
      USE fixed_occ, ONLY : tfixed_occ, f_inp
      USE ktetra,    ONLY : ntetra, tetra, ltetra
      USE klist,     ONLY : lgauss, ngauss, degauss
      USE wvfct,     ONLY : nbnd
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: i
      !
      !
      IF ( locc_read ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "OCCUPATIONS" )
         !
         CALL iotk_scan_dat( iunpun, "SMEARING_METHOD", lgauss )
         !
         IF ( lgauss ) THEN
            !
            CALL iotk_scan_dat( iunpun, "SMEARING_TYPE", ngauss )
            !
            CALL iotk_scan_dat( iunpun, "SMEARING_PARAMETER", degauss )
            !
            degauss = degauss * e2
            !
         END IF
         !
         CALL iotk_scan_dat( iunpun, "TETRAHEDRON_METHOD", ltetra )
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
         END IF
         !
         CALL iotk_scan_dat( iunpun, "FIXED_OCCUPATIONS", tfixed_occ )
         !
         IF ( tfixed_occ ) THEN
            !
            IF ( .NOT. ALLOCATED( f_inp ) ) THEN
               IF ( nspin == 4 ) THEN
                  ALLOCATE( f_inp( nbnd, 1 ) )
               ELSE
                  ALLOCATE( f_inp( nbnd, nspin ) )
               END IF
            END IF
            !
            CALL iotk_scan_dat( iunpun, "INPUT_OCC_UP", f_inp(:,1) )
            !
            IF ( lsda ) &
               CALL iotk_scan_dat( iunpun, "INPUT_OCC_DOWN", f_inp(:,2) )
            !
         END IF         
         !
         CALL iotk_scan_end( iunpun, "OCCUPATIONS" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( lgauss, ionode_id )
      !
      IF ( lgauss ) THEN
         !
         CALL mp_bcast( ngauss,  ionode_id )
         CALL mp_bcast( degauss, ionode_id )
         !
      END IF
      !
      CALL mp_bcast( ltetra, ionode_id )
      !
      IF ( ltetra ) THEN
         !
         CALL mp_bcast( ntetra, ionode_id )
         CALL mp_bcast( tetra,  ionode_id )
         !
      END IF
      !
      CALL mp_bcast( tfixed_occ, ionode_id )
      !
      IF ( tfixed_occ ) CALL mp_bcast( f_inp, ionode_id )
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
      USE basis,    ONLY : natomwfc
      USE lsda_mod, ONLY : lsda, isk
      USE klist,    ONLY : nkstot, wk, nelec
      USE wvfct,    ONLY : et, wg, nbnd
      USE ener,     ONLY : ef
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER :: ik, ik_eff, num_k_points
      !
      !
      IF ( lbs_read ) RETURN
      !
      IF ( .NOT. lspin_read ) &
         CALL errore( 'read_band_structure', 'read spin first', 1 )
      IF ( .NOT. lbz_read ) &
         CALL errore( 'read_band_structure', 'read band_structure first', 1 )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMIC_WFC", natomwfc )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd )
         !
         CALL iotk_scan_dat( iunpun, "FERMI_ENERGY", ef )
         !
         ef = ef * e2
         !
         num_k_points = nkstot
         !
         IF ( lsda ) num_k_points = nkstot / 2
         !
         CALL iotk_scan_begin( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
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
               CALL iotk_scan_dat( iunpun, "ET.1", et(:,ik)  )
               !
               CALL iotk_scan_dat( iunpun, "OCC.1", wg(:,ik) )
               !
               ik_eff = ik + num_k_points
               !
               isk(ik_eff) = 2
               !
               CALL iotk_scan_dat( iunpun, "ET.2", et(:,ik_eff) )
               !
               CALL iotk_scan_dat( iunpun, "OCC.2", wg(:,ik_eff) )
               !
            ELSE
               !
               isk(ik) = 1
               !
               CALL iotk_scan_dat( iunpun, "ET.1", et(:,ik) )
               !
               CALL iotk_scan_dat( iunpun, "OCC.1", wg(:,ik) )
               !
            END IF
            !
            CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
         END DO k_points_loop
         !
         et(:,:) = et(:,:) * e2
         !
         FORALL( ik = 1:nkstot ) wg(:,ik) = wg(:,ik) * wk(ik)
         !
         CALL iotk_scan_end( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( nelec,    ionode_id )
      CALL mp_bcast( natomwfc, ionode_id )
      CALL mp_bcast( nbnd,     ionode_id )
      CALL mp_bcast( isk,      ionode_id )
      CALL mp_bcast( et,       ionode_id )
      CALL mp_bcast( wg,       ionode_id )
      CALL mp_bcast( ef,       ionode_id )
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
      USE control_flags,        ONLY : twfcollect
      USE cell_base,            ONLY : tpiba2
      USE lsda_mod,             ONLY : nspin, isk
      USE klist,                ONLY : nkstot, wk, nelec, nks, xk, ngk
      USE wvfct,                ONLY : npw, npwx, g2kin, et, wg, &
                                       igk_l2g, nbnd
      USE wavefunctions_module, ONLY : evc, evc_nc
      USE reciprocal_vectors,   ONLY : ig_l2g
      USE io_files,             ONLY : nwordwfc, iunwfc
      USE gvect,                ONLY : ngm, ngm_g, ig1, ig2, ig3, g, ecutwfc
      USE noncollin_module,     ONLY : noncolin, npol                             
      USE mp_global,            ONLY : kunit, nproc, nproc_pool
      USE mp_global,            ONLY : my_pool_id, &
                                       intra_pool_comm, inter_pool_comm
      USE mp,                   ONLY : mp_sum, mp_max
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256)   :: filename
      INTEGER              :: i, ig, ik, ipol, ik_eff, num_k_points
      INTEGER, ALLOCATABLE :: kisort(:)
      INTEGER              :: npool, nkbl, nkl, nkr, npwx_g
      INTEGER              :: ike, iks, npw_g, ispin
      INTEGER, ALLOCATABLE :: ngk_g(:)
      INTEGER, ALLOCATABLE :: itmp(:,:)
      LOGICAL              :: exst, opnd
      REAL(DP)             :: scalef
      !
      !
      INQUIRE( UNIT = iunwfc, OPENED = opnd )
      !
      IF ( .NOT. opnd ) &
         CALL errore( 'read_wavefunctions', &
                    & 'wavefunctions unit (iunwfc) is not opened', 1 )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id )
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
      CALL mp_sum( ngm_g, intra_pool_comm )      
      !
      ! ... collect all G vectors across processors within the pools
      !
      ALLOCATE( itmp( 3, ngm_g ) )
      !
      itmp = 0
      !
      DO ig = 1, ngm
         !
         itmp(1,ig_l2g(ig)) = ig1(ig)
         itmp(2,ig_l2g(ig)) = ig2(ig)
         itmp(3,ig_l2g(ig)) = ig3(ig)
         !
      END DO
      !
      CALL mp_sum( itmp, intra_pool_comm )
      !
      ! ... build the G+k array indexes
      !
      ALLOCATE( kisort( npwx ) )
      !
      DO ik = 1, nkstot
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
      CALL mp_sum( ngk_g )
      !
      ! ... compute the Maximum G vector index among all G+k an processors
      !
      npw_g = MAXVAL( igk_l2g(:,:) )
      !
      CALL mp_max( npw_g )
      !
      ! ... compute the Maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g( 1:nkstot ) )
      !      
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_scan_begin( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
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
            CALL iotk_scan_begin( iunpun, &
                                  "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            CALL iotk_scan_begin( iunpun, "WFC", FOUND = twfcollect  )
            !
            IF ( twfcollect ) CALL iotk_scan_end( iunpun, "WFC" )
            !
         END IF
         !
         CALL mp_bcast( twfcollect, ionode_id )
         !
         IF ( .NOT. twfcollect ) THEN
            !
            IF ( ionode ) THEN
               !
               CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
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
            isk(ik) = 1
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin ) )
               !
            END IF
            !
            !
            CALL read_wfc( iunout, ik, nkstot, kunit, ispin, nspin, &
                 evc, npw_g, nbnd, igk_l2g(:,ik-iks+1),   &
                 ngk(ik-iks+1), filename, scalef )
            !
            IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
               !
               CALL davcio( evc, nwordwfc, iunwfc, (ik-iks+1), + 1 )
               !
            END IF
            !
            ispin = 2
            ik_eff = ik + num_k_points
            isk(ik_eff) = 2
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin ) )
               !
            END IF
            !
            CALL read_wfc( iunout, ik_eff, nkstot, kunit, ispin, nspin, &
                           evc, npw_g, nbnd, igk_l2g(:,ik_eff-iks+1),   &
                           ngk(ik_eff-iks+1), filename, scalef )
            !
            IF ( ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
               !
               CALL davcio( evc, nwordwfc, iunwfc, (ik_eff-iks+1), + 1 )
               !
            END IF
            !
         ELSE
            !
            isk(ik) = 1
            !
            IF ( noncolin ) THEN
               !
               DO ipol = 1, npol
                  !
                  IF ( ionode ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evc', ik, ipol ) )
                     !
                  END IF
                  !
                  CALL read_wfc( iunout, ik, nkstot, kunit, ispin,     &
                                 nspin, evc_nc(:,ipol,:), npw_g, nbnd, &
                                 igk_l2g(:,ik-iks+1), ngk(ik-iks+1),   &
                                 filename, scalef )
                  !
               END DO
               !
            ELSE
               !
               IF ( ionode ) THEN
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc', ik ) )
                  !
               END IF
               !
               CALL read_wfc( iunout, ik, nkstot, kunit, ispin, nspin,    &
                              evc, npw_g, nbnd, igk_l2g(:,ik-iks+1),      &
                              ngk(ik-iks+1), filename, scalef )
               !
            END IF
            !
            IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
               !
               IF ( noncolin ) THEN
                  !
                  CALL davcio( evc_nc, nwordwfc, iunwfc, (ik-iks+1), + 1 )
                  !
               ELSE
                  !
                  CALL davcio( evc, nwordwfc, iunwfc, (ik-iks+1), + 1 )
                  !
               END IF
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
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE" )
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
    SUBROUTINE read_phonon( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE control_flags, ONLY : modenum
      USE klist,         ONLY : xqq
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "PHONON" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_MODES", modenum )
         !
         CALL iotk_scan_dat( iunpun, "Q-POINT", xqq(:) )
         !
         CALL iotk_scan_end( iunpun, "PHONON" )
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( modenum, ionode_id )
      CALL mp_bcast( xqq,     ionode_id )      
      !
      RETURN
      !
    END SUBROUTINE read_phonon
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rho( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2
      USE lsda_mod,  ONLY : nspin
      USE scf,       ONLY : rho
      USE sticks,    ONLY : dfftp
      USE mp_global, ONLY : me_pool, nproc_pool
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256) :: rho_file_base
      !
      !
      rho_file_base = TRIM( dirname ) // '/charge-density'
      !
      CALL read_rho_xml( rho_file_base, rho(:,1), &
           nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
      !
      IF ( nspin == 2 ) THEN
         !
         rho_file_base = TRIM( dirname ) // '/spin-polarization'
         !
         CALL read_rho_xml( rho_file_base, rho(:,2), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         rho(:,1) = ( rho(:,1) + rho(:,2) ) / 2.d0  ! ... this is rho_up ...
         !
         rho(:,2) = rho(:,1) - rho(:,2)             ! ... and this is rho_dw
         !
      ELSE IF ( nspin == 4 ) THEN
         !
         rho_file_base = TRIM( dirname ) // '/magnetization.x'
         !
         CALL read_rho_xml( rho_file_base, rho(:,2), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         rho_file_base = TRIM( dirname ) // '/magnetization.y'
         !
         CALL read_rho_xml( rho_file_base, rho(:,3), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
         rho_file_base = TRIM( dirname ) // '/magnetization.z'
         !
         CALL read_rho_xml( rho_file_base, rho(:,4), &
                            nr1, nr2, nr3, nrx1, nrx2, dfftp%ipp, dfftp%npp )
         !
      END IF
      !
      ierr = 0
      !
    END SUBROUTINE read_rho
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
                            & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id )
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
      CALL mp_bcast( idum, ionode_id )
      !
      RETURN
      !
    END SUBROUTINE read_
    !
END MODULE pw_restart
