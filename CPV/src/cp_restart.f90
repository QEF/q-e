!
! Copyright (C) 2005-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE cp_restart
  !-----------------------------------------------------------------------------
  !
  ! ... This module contains subroutines to write and read data required to
  ! ... restart a calculation from the disk
  !
  !
  USE iotk_module
#if defined (__OLDXML)
  USE qexml_module,     ONLY : qexml_init,qexml_openfile, qexml_closefile, &
                          qexml_write_header,qexml_write_control, qexml_write_cell, &
                          qexml_write_ions, qexml_write_planewaves, qexml_write_spin, &
                          qexml_write_xc, qexml_write_occ, qexml_write_bz, qexml_write_para, &
                          qexml_write_bands_info,qexml_write_bands_cp,qexml_write_status_cp, &
                          qexml_kpoint_dirname, qexml_read_header, qexml_read_status_cp, &
                          qexml_read_ions, qexml_read_spin, qexml_read_occ, &
                          qexml_read_bands_info, qexml_read_bands_cp, &
                          qexml_wfc_filename, qexml_restart_dirname
  USE io_files,  ONLY : xmlpun, qexml_version, qexml_version_init
  USE xml_io_base, ONLY  : write_wfc
#endif
  USE xml_io_base,     ONLY  : read_wfc, write_rho, read_print_counter, create_directory
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE io_files,  ONLY : prefix, iunpun
  USE mp,        ONLY : mp_bcast
  USE parser,    ONLY : version_compare
  USE matrix_inversion
  !
  IMPLICIT NONE
  !
  SAVE
  !
#if defined (__OLDXML)
  PRIVATE :: read_cell
#endif
  !
  INTEGER, PRIVATE :: iunout
  !
  ! variables to describe qexml current version
  ! and back compatibility
  !
  LOGICAL, PRIVATE :: qexml_version_before_1_4_0 = .FALSE.
  !
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_writefile( ndw, ascii, nfi, simtime, acc, nk, xk,          &
                             wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,   &
                             taui, cdmi, stau0, svel0, staum, svelm, force,  &
                             vnhp, xnhp0, xnhpm, nhpcl, nhpdim, occ0, occm,  &
                             lambda0,lambdam, xnhe0, xnhem, vnhe, ekincm,    &
                             et, rho, c02, cm2, ctot, iupdwn, nupdwn,        &
                             iupdwn_tot, nupdwn_tot, wfc, mat_z ) ! BS added wfc
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : gamma_only, force_pairing, trhow, &
                                           tksw, twfcollect, do_makov_payne, &
                                           smallmem, llondon, lxdm, ts_vdw,  &
                                           tfor, tpre
      USE control_flags,            ONLY : lwfpbe0nscf, lwfnscf, lwf ! Lingzhu Kong
      USE constants,                ONLY : e2
      USE dener,                    ONLY : detot
      USE io_files,                 ONLY : psfile, pseudo_dir, iunwfc, &
                                           nwordwfc, tmp_dir, diropn
      USE mp_images,                ONLY : intra_image_comm, me_image, &
                                           nproc_image
      USE mp_pools,                 ONLY : nproc_pool, intra_pool_comm, root_pool, inter_pool_comm
      USE mp_bands,                 ONLY : me_bgrp, nproc_bgrp, &
                                           my_bgrp_id, intra_bgrp_comm, &
                                           inter_bgrp_comm, root_bgrp, &
                                           ntask_groups
      USE mp_diag,                  ONLY : nproc_ortho
      USE mp_world,                 ONLY : world_comm, nproc
      USE run_info,                 ONLY : title
      USE gvect,                    ONLY : ngm, ngm_g
      USE gvecs,                    ONLY : ngms_g, ecuts, dual
      USE gvecw,                    ONLY : ngw, ngw_g, ecutwfc
      USE gvect,                    ONLY : ig_l2g, mill
      USE electrons_base,           ONLY : nspin, nelt, nel, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, s_to_r, ainv ! BS added ainv
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           amass, iforce, ind_bck
      USE funct,                    ONLY : get_dft_name, get_inlc, &
           dft_is_hybrid, get_exx_fraction, get_screening_parameter, &
           dft_is_nonlocc, get_nonlocc_name
      USE ldaU_cp,                  ONLY : lda_plus_U, ns, ldmx,Hubbard_l, &
                                           Hubbard_lmax, Hubbard_U
      USE energies,                 ONLY : enthal, ekin, eht, esr, eself, &
                                           epseu, enl, exc, vave
      USE mp,                       ONLY : mp_sum, mp_barrier
      USE fft_base,                 ONLY : dfftp, dffts, dfftb
      USE uspp_param,               ONLY : n_atom_wfc, upf
      USE global_version,           ONLY : version_number
      USE cp_main_variables,        ONLY : descla
      USE cp_interfaces,            ONLY : collect_lambda, collect_zmat
      USE kernel_table,             ONLY : vdw_table_name, kernel_file_name
      USE london_module,            ONLY : scal6, lon_rcut, in_c6
      USE tsvdw_module,             ONLY : vdw_isolated, vdw_econv_thr
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN) :: ndw          !
      LOGICAL,               INTENT(IN) :: ascii        !
      INTEGER,               INTENT(IN) :: nfi          ! index of the current step
      REAL(DP),              INTENT(IN) :: simtime      ! simulated time
      REAL(DP),              INTENT(IN) :: acc(:)       !  
      INTEGER,               INTENT(IN) :: nk           ! number of kpoints
      REAL(DP),              INTENT(IN) :: xk(:,:)      ! k-points coordinates 
      REAL(DP),              INTENT(IN) :: wk(:)        ! k-points weights
      REAL(DP),              INTENT(IN) :: ht(3,3)      ! 
      REAL(DP),              INTENT(IN) :: htm(3,3)     ! 
      REAL(DP),              INTENT(IN) :: htvel(3,3)   ! 
      REAL(DP),              INTENT(IN) :: gvel(3,3)    ! 
      REAL(DP),              INTENT(IN) :: xnhh0(3,3)   ! 
      REAL(DP),              INTENT(IN) :: xnhhm(3,3)   ! 
      REAL(DP),              INTENT(IN) :: vnhh(3,3)    ! 
      REAL(DP),              INTENT(IN) :: taui(:,:)    ! 
      REAL(DP),              INTENT(IN) :: cdmi(:)      ! 
      REAL(DP),              INTENT(IN) :: stau0(:,:)   ! 
      REAL(DP),              INTENT(IN) :: svel0(:,:)   ! 
      REAL(DP),              INTENT(IN) :: staum(:,:)   ! 
      REAL(DP),              INTENT(IN) :: svelm(:,:)   ! 
      REAL(DP),              INTENT(IN) :: force(:,:)   ! 
      REAL(DP),              INTENT(IN) :: xnhp0(:)     ! 
      REAL(DP),              INTENT(IN) :: xnhpm(:)     ! 
      REAL(DP),              INTENT(IN) :: vnhp(:)      ! 
      INTEGER,               INTENT(IN) :: nhpcl        ! 
      INTEGER,               INTENT(IN) :: nhpdim       ! 
      REAL(DP),              INTENT(IN) :: occ0(:)      !  occupations of electronic states
      REAL(DP),              INTENT(IN) :: occm(:)      ! 
      REAL(DP),              INTENT(IN) :: lambda0(:,:,:) ! 
      REAL(DP),              INTENT(IN) :: lambdam(:,:,:) ! 
      REAL(DP),              INTENT(IN) :: xnhe0        ! 
      REAL(DP),              INTENT(IN) :: xnhem        ! 
      REAL(DP),              INTENT(IN) :: vnhe         ! 
      REAL(DP),              INTENT(IN) :: ekincm       ! 
      REAL(DP),              INTENT(IN) :: et(:,:)      !  eigenvalues
      REAL(DP),              INTENT(IN) :: rho(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: c02(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: cm2(:,:)     ! 
      COMPLEX(DP),           INTENT(IN) :: ctot(:,:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: nupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn_tot(:)! 
      INTEGER,               INTENT(IN) :: nupdwn_tot(:)! 
      REAL(DP),              INTENT(IN) :: wfc(:,:)     ! BS 
      REAL(DP),    OPTIONAL, INTENT(IN) :: mat_z(:,:,:) ! 
      !
      LOGICAL               :: write_charge_density
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname
      CHARACTER(LEN=320)    :: filename, sourcefile
      CHARACTER(LEN=4)      :: cspin
      INTEGER               :: kunit, ib, ik_eff
      INTEGER               :: k1, k2, k3
      INTEGER               :: nk1, nk2, nk3
      INTEGER               :: j, i, iss, ig, nspin_wfc, iss_wfc
      INTEGER               :: is, ia, isa, ik, ierr
      INTEGER,  ALLOCATABLE :: ftmp(:,:)
      INTEGER,  ALLOCATABLE :: ityp(:)
      REAL(DP), ALLOCATABLE :: tau(:,:)
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      REAL(DP)              :: omega, htm1(3,3), h(3,3)
      REAL(DP)              :: a1(3), a2(3), a3(3)
      REAL(DP)              :: b1(3), b2(3), b3(3)
      REAL(DP)              :: nelec
      REAL(DP)              :: scalef
      LOGICAL               :: lsda
      REAL(DP)              :: s0, s1, cclock
      INTEGER               :: nbnd_tot
      INTEGER               :: natomwfc, nbnd_
      REAL(DP), ALLOCATABLE :: mrepl(:,:)
      CHARACTER(LEN=256)    :: tmp_dir_save
      LOGICAL               :: exst
      INTEGER               :: inlc
      CHARACTER(iotk_attlenx)  :: attr
      REAL(DP), ALLOCATABLE :: temp_vec(:), wfc_temp(:,:) ! BS 
      !
      k1  = 0
      k2  = 0
      k3  = 0
      nk1 = 0
      nk2 = 0
      nk3 = 0
      !
      ! ... subroutine body
      !
      write_charge_density = trhow
      !
      IF( nspin > 1 .AND. .NOT. force_pairing ) THEN
         !
         !  check if the array storing wave functions is large enought
         !
         IF( SIZE( c02, 2 ) < ( iupdwn( 2 ) + nupdwn(1) - 1 ) ) &
            CALL errore('cp_writefile',' wrong wave functions dimension ', 1 )
         !
      END IF
      !
      IF(  nupdwn_tot(1) < nupdwn(1) ) &
         CALL errore( " writefile ", " wrong number of states ", 1 )
      !
      nbnd_    = nupdwn(1) 
      nbnd_tot = MAX( nupdwn(1), nupdwn_tot(1) )
      nelec = nelt
      !
      ! ... Cell related variables
      ! ... Dirty trick to avoid bogus complaints because ht in intent(in)
      !
      h = ht
      CALL invmat( 3, h, htm1, omega )
      h = TRANSPOSE( ht )
      !
      a1 = ht(1,:)/alat
      a2 = ht(2,:)/alat
      a3 = ht(3,:)/alat
      !
      CALL recips( a1, a2, a3, b1, b2, b3 )
      !
      ! ... Compute array ityp, and tau
      !
      ALLOCATE( ityp( nat ) )
      ALLOCATE( tau( 3, nat ) )
      !
      isa = 0
      !
      DO is = 1, nsp
         !
         DO ia = 1, na(is)
            !
            isa = isa + 1
            ityp(isa) = is
            !
         END DO
         !
      END DO
      !
      natomwfc =  n_atom_wfc ( nat, ityp ) 
      !
      CALL s_to_r( stau0, tau, na, nsp, h )
      !   
      lsda = ( nspin == 2 )
      !
      ALLOCATE( ftmp( nbnd_tot , nspin ) )
      !
      ftmp = 0.0d0
      !
      DO iss = 1, nspin
         !
         ftmp( 1:nupdwn(iss), iss ) = occ0( iupdwn(iss) : iupdwn(iss) + nupdwn(iss) - 1 )
         !
      END DO
      !
#if !defined (__OLDXML)
      !
      CALL errore('cp_writefile','called in wrong case',1)
      !
#else
      !
      ! ... Some ( CP/FPMD ) default values
      !
      IF ( nspin == 2 ) THEN
         kunit = 2
      ELSE
         kunit = 1
      END IF
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
      CALL errore( 'cp_writefile', &
                   'no free units to write wavefunctions', ierr )
      !
      dirname = qexml_restart_dirname( tmp_dir, prefix, ndw )
      !
      ! ... Create main restart directory
      !
      CALL create_directory( dirname )
      !
      ! ... Create k-points subdirectories
      ! ... note: in FPMD and CP k-points are not distributed to processors
      !
      DO i = 1, nk
         !
         CALL create_directory( qexml_kpoint_dirname( dirname, i ) )
         !
      END DO
      IF ( ionode ) THEN
         !
         ! ... Open XML descriptor
         !
         WRITE( stdout, '(/,3X,"writing restart file: ",A)' ) TRIM( dirname )
         !
         CALL qexml_init( iunpun )
         CALL qexml_openfile( TRIM( dirname ) // TRIM( xmlpun ), &
                             & 'write', BINARY = .FALSE., IERR = ierr  )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_writefile ', 'cannot open restart file for writing', ierr )
      !
      s0 = cclock()
      !
      IF ( ionode ) THEN

!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_header( "CP", TRIM(version_number) )
         !
!-------------------------------------------------------------------------------
! ... this flag is used to check if the file can be used for post-processing
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_control( PP_CHECK_FLAG=.TRUE. )
         !
!-------------------------------------------------------------------------------
! ... STATUS
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_status_cp( nfi,simtime,"pico-seconds",TRIM(title), &
                                  ekin, eht, esr, eself, epseu, enl, exc, vave, enthal, &   
                                  'Hartree' )
         !
!-------------------------------------------------------------------------------
! ... CELL
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_cell( ibrav, celldm, alat, a1, a2, a3, b1, b2, b3, &
                          "Bohr","Bohr","2 pi / a", &
                          do_makov_payne, .FALSE., .FALSE. )
         !
!-------------------------------------------------------------------------------
! ... IONS
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_ions( nsp, nat, atm, ityp(ind_bck(:)), &
                          psfile, pseudo_dir, amass, 'a.m.u.', tau(:,ind_bck(:)), &
                          'Bohr', iforce(:,ind_bck(:)), dirname, 1.D0 )
         !
!-------------------------------------------------------------------------------
! ... PLANE_WAVES
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_planewaves( ecutwfc/e2, ecutwfc*dual/e2, ngw_g, gamma_only, &
              dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, & 
              dffts%nr1, dffts%nr2, dffts%nr3, ngms_g,&
              dfftb%nr1, dfftb%nr2, dfftb%nr3, mill, .FALSE.,'Hartree' )
         !
!-------------------------------------------------------------------------------
! ... SPIN
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_spin( lsda, .FALSE., 1, .FALSE., .TRUE. )
         !
!-------------------------------------------------------------------------------
! ... EXCHANGE_CORRELATION
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         inlc = get_inlc()
         !
         CALL qexml_write_xc( DFT = dft_name, NSP = nsp, LDA_PLUS_U = lda_plus_u, &
                        HUBBARD_LMAX = Hubbard_lmax,      &
                        HUBBARD_L = Hubbard_l, HUBBARD_U = Hubbard_U, &
                        INLC = inlc, VDW_TABLE_NAME = vdw_table_name, &
                        PSEUDO_DIR = pseudo_dir, DIRNAME = dirname,   &
                        LLONDON = llondon, LONDON_S6 = scal6,         &
                        LONDON_RCUT = lon_rcut, LXDM = lxdm,          &
                        TS_VDW = ts_vdw, VDW_ISOLATED = vdw_isolated )
         !
!-------------------------------------------------------------------------------
! ... OCCUPATIONS
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_occ( LGAUSS = .FALSE., LTETRA = .FALSE., &
                         TFIXED_OCC = .TRUE., LSDA = lsda, NSTATES_UP = nupdwn_tot(1), &
                         NSTATES_DW = nupdwn_tot(2), INPUT_OCC = DBLE( ftmp ) )
         !
!-------------------------------------------------------------------------------
! ... BRILLOUIN_ZONE
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_bz( nk, xk, wk, k1, k2, k3, nk1, nk2, nk3, &
                        '2 pi / a',0.0_DP )
         !
!-------------------------------------------------------------------------------
! ... PARALLELISM
!-------------------------------------------------------------------------------
         !
         !
         CALL qexml_write_para( kunit, nproc, nproc_pool, nproc_image, &
                                ntask_groups, nproc_bgrp, nproc_ortho )
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... CHARGE-DENSITY
!-------------------------------------------------------------------------------
      !
      IF (write_charge_density) CALL write_rho ( dirname, rho, nspin )
      !
!-------------------------------------------------------------------------------
! ... LDA+U OCCUPATIONS (compatibility with PWscf)
!-------------------------------------------------------------------------------
      !
      IF ( lda_plus_u ) THEN
         !
         IF ( ionode ) THEN
            filename = dirname // '/occup'
            OPEN (UNIT=iunout,FILE=filename,FORM ='formatted',STATUS='unknown')
            WRITE( iunout, * , iostat = ierr) ns
         END IF
         CALL mp_bcast( ierr, ionode_id, intra_image_comm )
         IF ( ierr/=0 ) CALL errore('cp_writefile', 'Writing ldaU ns', 1)
         IF ( ionode ) THEN
            CLOSE( UNIT = iunout, STATUS = 'KEEP' )
         ENDIF
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... TIMESTEPS
!-------------------------------------------------------------------------------
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_attr( attr, "nt", 2, FIRST = .TRUE. )
         !
         CALL iotk_write_begin( iunpun, "TIMESTEPS", attr )
         !
         ! ... STEP0
         !
         CALL iotk_write_begin( iunpun, "STEP0" )
         !
         CALL iotk_write_dat( iunpun, "ACCUMULATORS", acc )
         !
         CALL iotk_write_begin( iunpun, "IONS_POSITIONS" )
         CALL iotk_write_dat(   iunpun, "stau",  stau0(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "svel",  svel0(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "taui",  taui(1:3,1:nat),    COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "cdmi",  cdmi(1:3),          COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "force", force(1:3,1:nat),   COLUMNS=3 )
         CALL iotk_write_end(   iunpun, "IONS_POSITIONS" )
         !
         CALL iotk_write_begin( iunpun, "IONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "nhpcl", nhpcl )
         CALL iotk_write_dat(   iunpun, "nhpdim", nhpdim )
         CALL iotk_write_dat(   iunpun, "xnhp",  xnhp0(1:nhpcl*nhpdim) )
         CALL iotk_write_dat(   iunpun, "vnhp",  vnhp(1:nhpcl*nhpdim) )
         CALL iotk_write_end(   iunpun, "IONS_NOSE" )
         !
         CALL iotk_write_dat( iunpun, "ekincm", ekincm )
         !
         CALL iotk_write_begin( iunpun, "ELECTRONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhe", xnhe0 )
         CALL iotk_write_dat(   iunpun, "vnhe", vnhe )
         CALL iotk_write_end(   iunpun, "ELECTRONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "CELL_PARAMETERS" )
         CALL iotk_write_dat(   iunpun, "ht",    ht )
         CALL iotk_write_dat(   iunpun, "htvel", htvel )
         CALL iotk_write_dat(   iunpun, "gvel",  gvel )
         CALL iotk_write_end(   iunpun, "CELL_PARAMETERS" )
         !
         CALL iotk_write_begin( iunpun, "CELL_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhh", xnhh0 )
         CALL iotk_write_dat(   iunpun, "vnhh", vnhh )
         CALL iotk_write_end(   iunpun, "CELL_NOSE" )
         !
         CALL iotk_write_end( iunpun, "STEP0" )
         !
         ! ... STEPM
         !
         CALL iotk_write_begin( iunpun, "STEPM" )
         !
         CALL iotk_write_begin( iunpun, "IONS_POSITIONS" )
         CALL iotk_write_dat(   iunpun, "stau", staum(1:3,1:nat),  COLUMNS=3 )
         CALL iotk_write_dat(   iunpun, "svel", svelm(1:3,1:nat),  COLUMNS=3 )
         CALL iotk_write_end(   iunpun, "IONS_POSITIONS" )
         !
         CALL iotk_write_begin( iunpun, "IONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "nhpcl", nhpcl )
         CALL iotk_write_dat(   iunpun, "nhpdim", nhpdim )
         CALL iotk_write_dat(   iunpun, "xnhp",  xnhpm(1:nhpcl*nhpdim) )
         CALL iotk_write_end(   iunpun, "IONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "ELECTRONS_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhe", xnhem )
         CALL iotk_write_end(   iunpun, "ELECTRONS_NOSE" )
         !
         CALL iotk_write_begin( iunpun, "CELL_PARAMETERS" )
         CALL iotk_write_dat(   iunpun, "ht",    htm )
         CALL iotk_write_end(   iunpun, "CELL_PARAMETERS" )
         !
         CALL iotk_write_begin( iunpun, "CELL_NOSE" )
         CALL iotk_write_dat(   iunpun, "xnhh", xnhhm )
         CALL iotk_write_end(   iunpun, "CELL_NOSE" )
         !
         CALL iotk_write_end( iunpun, "STEPM" )
         !
         CALL iotk_write_end( iunpun, "TIMESTEPS" )
         !
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... BAND_STRUCTURE_INFO
!-------------------------------------------------------------------------------
      !
      IF ( ionode ) THEN
         ! 
         CALL qexml_write_bands_info(  nk, natomwfc, &
                                       nbnd_tot, nupdwn_tot(1),nupdwn_tot(2),&
                                       nspin, nelec, nel(1), nel(2), &
                                       "Hartree", "2 pi / a")
         !
!-------------------------------------------------------------------------------
! ... EIGENVALUES
!-------------------------------------------------------------------------------
         !
         CALL qexml_write_bands_cp( nbnd_tot, nk, nspin, iupdwn, nupdwn, xk, wk, et, tksw, &
              occ0, occm, "Hartree", "2 pi / a", iunout ,dirname )
         !
         !
         CALL iotk_write_begin( iunpun, "EIGENVECTORS" )
         !
         CALL iotk_write_dat  ( iunpun, "MAX_NUMBER_OF_GK-VECTORS", ngw_g )
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... EIGENVECTORS
!-------------------------------------------------------------------------------
      !
      ! ... Beware: omega may be negative if axis are left-handed!
      !
      scalef = 1.D0 / SQRT( ABS (omega) )
      !
      k_points_loop2: DO ik = 1, nk

         IF( ionode ) THEN

            CALL iotk_write_begin( iunpun, "K-POINT" // TRIM( iotk_index( ik ) ) )
            !
            ! ... G+K vectors
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_GK-VECTORS", ngw_g )
            !
            !
            filename = TRIM( qexml_wfc_filename( ".", 'gkvectors', ik ) )
            !
            CALL iotk_link( iunpun, "GK-VECTORS", filename, CREATE = .FALSE., BINARY = .TRUE. )
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'gkvectors', ik ) )
            !
         END IF
         !
         IF( .NOT. smallmem ) THEN
            CALL write_gk( iunout, ik, filename )
         END IF
         !
         DO iss = 1, nspin
            ! 
            ik_eff = ik + ( iss - 1 ) * nk
            ! 
            iss_wfc = iss
            if( force_pairing ) iss_wfc = 1   ! only the WF for the first spin is allocated
            !
            IF( tksw ) THEN 
               ! 
               !   Save additional WF, 
               !   orthogonal KS states to be used for post processing and PW
               ! 
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( qexml_wfc_filename( ".", 'evc', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( qexml_wfc_filename( ".", 'evc', ik, iss ) )
                     !
                  END IF
                  !
                  IF( nspin == 2 ) THEN
                     CALL iotk_link( iunpun, "WFC" // TRIM( iotk_index (iss) ), &
                                     filename, CREATE = .FALSE., BINARY = .TRUE. )
                  ELSE
                     CALL iotk_link( iunpun, "WFC", filename, CREATE = .FALSE., BINARY = .TRUE. )
                  END IF
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( qexml_wfc_filename( dirname, 'evc', ik, iss ) )
                  !
                  END IF
                  !
               END IF
               !
               ib = iupdwn_tot( iss_wfc )
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,        &
                            ctot( :, ib : ib + nbnd_tot - 1 ), ngw_g, gamma_only,&
                            nbnd_tot, ig_l2g, ngw, filename, scalef, &
                            ionode, root_bgrp, intra_bgrp_comm, inter_bgrp_comm, intra_pool_comm )
               !
            END IF
            !
            IF( twfcollect ) THEN
               !
               !  Save wave function at time t
               !
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( qexml_wfc_filename( ".", 'evc0', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( qexml_wfc_filename( ".", 'evc0', ik, iss ) )
                     !
                  END IF
                  !
                  CALL iotk_link( iunpun, "WFC0" // TRIM( iotk_index (iss) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( qexml_wfc_filename( dirname, 'evc0', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( qexml_wfc_filename( dirname, 'evc0', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,     &
                               c02( :, ib : ib + nbnd_ - 1 ), ngw_g, gamma_only, &
                               nbnd_, ig_l2g, ngw, filename, scalef, &
                               ionode, root_bgrp, intra_bgrp_comm, inter_bgrp_comm, intra_pool_comm )
               !
               !  Save wave function at time t - dt
               !
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( qexml_wfc_filename( ".", 'evcm', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( qexml_wfc_filename( ".", 'evcm', ik, iss ) )
                     !
                  END IF
                  !
                  CALL iotk_link( iunpun, "WFCM" // TRIM( iotk_index (iss) ), &
                                  filename, CREATE = .FALSE., BINARY = .TRUE. )
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( qexml_wfc_filename( dirname, 'evcm', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( qexml_wfc_filename( dirname, 'evcm', ik, iss ) )
                     !
                  END IF
                  !
               END IF
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, iss, nspin,     &
                               cm2( :, ib : ib + nbnd_ - 1 ), ngw_g, gamma_only, &
                               nbnd_, ig_l2g, ngw, filename, scalef, &
                               ionode, root_bgrp, intra_bgrp_comm, inter_bgrp_comm, intra_pool_comm )
               !
            END IF
            !
            cspin = iotk_index( iss )
            !
            ! ... write matrix lambda to file
            !
            ALLOCATE( mrepl( nudx, nudx ) )
            !
            CALL collect_lambda( mrepl, lambda0(:,:,iss), descla(iss) )
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( qexml_wfc_filename( ".", 'lambda0', ik, iss ) )
               !
               CALL iotk_link( iunpun, "LAMBDA0" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .TRUE. )
               !
               CALL iotk_write_dat( iunpun, &
                                    "LAMBDA0" // TRIM( cspin ), mrepl )
               !=============================================================
               !  Lingzhu Kong
               IF ( lwfpbe0nscf .or. lwfnscf ) THEN
                  OPEN(60,file='cp_lambda.dat',status='unknown',form='formatted') 
                  DO j = 1, nudx
                     write(60, '(8f15.8)')(mrepl(i,j),i=1,nudx)
                  ENDDO
                  CLOSE(60)
               ENDIF
               !=============================================================
               !
               !
            END IF
            !
            CALL collect_lambda( mrepl, lambdam(:,:,iss), descla(iss) )
            !
            IF ( ionode ) THEN
               !
               filename = TRIM( qexml_wfc_filename( ".", 'lambdam', ik, iss ) )
               !
               CALL iotk_link( iunpun, "LAMBDAM" // TRIM( cspin ), &
                               filename, CREATE = .TRUE., BINARY = .TRUE. )
               !
               CALL iotk_write_dat( iunpun, &
                                    "LAMBDAM" // TRIM( cspin ), mrepl )
               !
            END IF
            !
            IF( PRESENT( mat_z ) ) THEN
               !
               CALL collect_zmat( mrepl, mat_z(:,:,iss), descla(iss) )
               !
               IF ( ionode ) THEN
                  !
                  filename = TRIM( qexml_wfc_filename( ".", 'mat_z', ik, iss ) )
                  !
                  CALL iotk_link( iunpun, "MAT_Z" // TRIM( cspin ), &
                                  filename, CREATE = .TRUE., BINARY = .TRUE. )
                  !
                  CALL iotk_write_dat( iunpun, "MAT_Z" // TRIM( cspin ), mrepl )
                  !
               END IF
               !
            END IF
            !
            DEALLOCATE( mrepl )
            !
         END DO
         !
         IF ( ionode ) &
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop2
      !
      IF ( ionode ) CALL iotk_write_end( iunpun, "EIGENVECTORS" )
      !
      !-------------------------------------------------------------------------
      ! BS : Wannier centers
      IF ( lwf ) THEN
         !
         IF ( ionode ) THEN
            CALL iotk_write_begin( iunpun, "WANNIER_CENTERS" )
            !
            ALLOCATE(temp_vec(3)) 
            !
            DO iss = 1, nspin
               !
               ALLOCATE( wfc_temp(3,nupdwn(iss)) )
               !
               temp_vec=0.0_DP
               wfc_temp=0.0_DP
               !
               j = 1 !wfc_temp count
               DO i = iupdwn(iss), iupdwn(iss) + nupdwn(iss) -1 
                     !
                     temp_vec(:) = MATMUL( ainv(:,:), wfc(:,i) )
                     !
                     temp_vec(:) = temp_vec(:) - floor (temp_vec(:))
                     !
                     temp_vec(:) = MATMUL( h(:,:), temp_vec(:) )
                     !
                     wfc_temp(:, j) = temp_vec(:)
                     j = j + 1
                     !
               END DO
               !
               CALL iotk_write_dat(   iunpun, "wanniercentres" // TRIM( iotk_index(iss)),  &
                                                                     & wfc_temp(1:3,1:nupdwn(iss)),  COLUMNS=3 )
               !
               DEALLOCATE(wfc_temp)
               !
            ENDDO
            !
            DEALLOCATE(temp_vec)
            !
            CALL iotk_write_end(   iunpun, "WANNIER_CENTERS" )
            !
         ENDIF
         !
      END IF
      ! BS : Wannier centers
      !-------------------------------------------------------------------------
      !
      IF ( ionode ) THEN
         !
         CALL qexml_closefile( 'write', IERR=ierr)
         !
      ENDIF
      !
      call mp_barrier( world_comm )
      !
      IF( (.NOT. twfcollect) .AND. (my_bgrp_id == 0) ) THEN
         !
         tmp_dir_save = tmp_dir
         tmp_dir = TRIM( qexml_kpoint_dirname( tmp_dir, 1 ) )
         !
         iunwfc = 10
         nwordwfc = SIZE( c02 )
         !
         CALL diropn ( iunwfc, 'wfc', 2*nwordwfc, exst )

         CALL davcio ( c02, 2*nwordwfc, iunwfc, 1, +1 )  ! save wave funct
         CALL davcio ( cm2, 2*nwordwfc, iunwfc, 2, +1 )  ! save wave funct
         !
         CLOSE( UNIT = iunwfc, STATUS = 'KEEP' )
         tmp_dir = tmp_dir_save
         !
      END IF
#endif 
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( ftmp )
      DEALLOCATE( tau  )
      DEALLOCATE( ityp )
      !
      s1 = cclock() 
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, &
                '(3X,"restart file written in ",F8.3," sec.",/)' ) ( s1 - s0 )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE cp_writefile
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_readfile( ndr, ascii, nfi, simtime, acc, nk, xk,   &
                            wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,     &
                            taui, cdmi, stau0, svel0, staum, svelm, force,    &
                            vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ0, occm,      &
                            lambda0, lambdam, b1, b2, b3, xnhe0, xnhem, vnhe, &
                            ekincm, c02, cm2, wfc, mat_z ) ! added wfc
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : gamma_only, force_pairing, iverbosity, twfcollect, lwf ! BS added lwf
      USE io_files,                 ONLY : iunpun, xmlpun, iunwfc, nwordwfc, &
                                           tmp_dir, diropn
      USE run_info,                 ONLY : title
      USE gvect,                    ONLY : ngm
      USE gvecw,                    ONLY : ngw, ngw_g
      USE electrons_base,           ONLY : nspin, nbnd, nelt, nel, &
                                           nupdwn, iupdwn, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, s_to_r, r_to_s
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           sort_tau, ityp, ions_cofmass
      USE gvect,       ONLY : ig_l2g, mill
      USE cp_main_variables,        ONLY : nprint_nfi, descla
      USE cp_interfaces,            ONLY : distribute_lambda, distribute_zmat
      USE mp,                       ONLY : mp_sum, mp_bcast
      USE mp_global,                ONLY : intra_image_comm, my_bgrp_id
      USE mp_global,                ONLY : root_bgrp, intra_bgrp_comm, inter_bgrp_comm, intra_pool_comm
      USE parameters,               ONLY : ntypx
      USE constants,                ONLY : eps8, angstrom_au, pi
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)    :: ndr          !  I/O unit number
      LOGICAL,               INTENT(IN)    :: ascii        !
      INTEGER,               INTENT(INOUT) :: nfi          ! index of the current step
      REAL(DP),              INTENT(INOUT) :: simtime      ! simulated time
      REAL(DP),              INTENT(INOUT) :: acc(:)       !
      INTEGER,               INTENT(IN)    :: nk           ! number of kpoints
      REAL(DP),              INTENT(INOUT) :: xk(:,:)      ! k-points coordinates
      REAL(DP),              INTENT(INOUT) :: wk(:)        ! k-points weights
      REAL(DP),              INTENT(INOUT) :: ht(3,3)      !
      REAL(DP),              INTENT(INOUT) :: htm(3,3)     !
      REAL(DP),              INTENT(INOUT) :: htvel(3,3)   !
      REAL(DP),              INTENT(INOUT) :: gvel(3,3)    !
      REAL(DP),              INTENT(INOUT) :: xnhh0(3,3)   !
      REAL(DP),              INTENT(INOUT) :: xnhhm(3,3)   !
      REAL(DP),              INTENT(INOUT) :: vnhh(3,3)    !
      REAL(DP),              INTENT(INOUT) :: taui(:,:)    !
      REAL(DP),              INTENT(INOUT) :: cdmi(:)      !
      REAL(DP),              INTENT(INOUT) :: stau0(:,:)   !
      REAL(DP),              INTENT(INOUT) :: svel0(:,:)   !
      REAL(DP),              INTENT(INOUT) :: staum(:,:)   !
      REAL(DP),              INTENT(INOUT) :: svelm(:,:)   !
      REAL(DP),              INTENT(INOUT) :: force(:,:)   ! 
      REAL(DP),              INTENT(INOUT) :: xnhp0(:)     !      
      REAL(DP),              INTENT(INOUT) :: xnhpm(:)     ! 
      REAL(DP),              INTENT(INOUT) :: vnhp(:)      !  
      INTEGER,               INTENT(INOUT) :: nhpcl        !  
      INTEGER,               INTENT(INOUT) :: nhpdim       !  
      REAL(DP),              INTENT(INOUT) :: occ0(:)      ! occupations
      REAL(DP),              INTENT(INOUT) :: occm(:)      !
      REAL(DP),              INTENT(INOUT) :: lambda0(:,:,:) !
      REAL(DP),              INTENT(INOUT) :: lambdam(:,:,:) !
      REAL(DP),              INTENT(INOUT) :: b1(3)        !
      REAL(DP),              INTENT(INOUT) :: b2(3)        !
      REAL(DP),              INTENT(INOUT) :: b3(3)        !
      REAL(DP),              INTENT(INOUT) :: xnhe0        !
      REAL(DP),              INTENT(INOUT) :: xnhem        !
      REAL(DP),              INTENT(INOUT) :: vnhe         !  
      REAL(DP),              INTENT(INOUT) :: ekincm       !  
      COMPLEX(DP),           INTENT(INOUT) :: c02(:,:)     ! 
      COMPLEX(DP),           INTENT(INOUT) :: cm2(:,:)     ! 
      REAL(DP),              INTENT(INOUT) :: wfc(:,:)     ! BS 
      REAL(DP),    OPTIONAL, INTENT(INOUT) :: mat_z(:,:,:) ! 
      !
      CHARACTER(LEN=256)   :: dirname, kdirname, filename
      CHARACTER(LEN=5)     :: kindex
      CHARACTER(LEN=4)     :: cspin
      INTEGER              :: strlen
      INTEGER              :: kunit
      INTEGER              :: k1, k2, k3
      INTEGER              :: nk1, nk2, nk3
      INTEGER              :: i, j, iss, ig, nspin_wfc, ierr, ik
      REAL(DP)             :: omega, htm1(3,3), hinv(3,3), scalef
      LOGICAL              :: found
      !
      ! ... variables read for testing pourposes
      !
      INTEGER               :: ibrav_
      CHARACTER(LEN=3)      :: atm_(ntypx)
      INTEGER               :: nat_, nsp_, na_
      INTEGER               :: nk_, ik_, nt_
      LOGICAL               :: gamma_only_ , lsda_
      REAL(DP)              :: alat_, a1_(3), a2_(3), a3_(3)
      REAL(DP)              :: zv_ 
      REAL(DP)              :: celldm_(6)
      INTEGER               :: iss_, nspin_, ngwt_, nbnd_ , nbnd_tot
      INTEGER               :: nstates_up_ , nstates_dw_ , ntmp, nel_(2)
      REAL(DP)              :: nelec_ 
      REAL(DP)              :: scalef_
      REAL(DP)              :: wk_
      INTEGER               :: nhpcl_, nhpdim_ 
      INTEGER               :: ib, nb
      INTEGER               :: ik_eff
      REAL(DP)              :: amass_(ntypx)
      INTEGER,  ALLOCATABLE :: ityp_(:) 
      INTEGER,  ALLOCATABLE :: isrt_(:) 
      REAL(DP), ALLOCATABLE :: tau_(:,:) 
      REAL(DP), ALLOCATABLE :: occ_(:) 
      INTEGER,  ALLOCATABLE :: if_pos_(:,:) 
      CHARACTER(LEN=256)    :: psfile_(ntypx)
      CHARACTER(LEN=80)     :: pos_unit
      REAL(DP)              :: s1, s0, cclock
      REAL(DP), ALLOCATABLE :: mrepl(:,:) 
      LOGICAL               :: exst, exist_wfc 
      CHARACTER(LEN=256)    :: tmp_dir_save
      INTEGER               :: io_bgrp_id
      CHARACTER(iotk_attlenx)  :: attr
      !
      ! ... look for an empty unit
      !
      CALL iotk_free_unit( iunout, ierr )
      CALL errore( 'cp_readfile', &
                   'no free units to read wavefunctions', ierr )
      !
#if !defined (__OLDXML)
      CALL errore('cp_readfile','called in wrong case',1)
#else
      kunit = 1
      found = .FALSE.
      exist_wfc = .FALSE.
      !
      dirname = qexml_restart_dirname( tmp_dir, prefix, ndr )
      !
      ! ... Open XML descriptor
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, '(/,3X,"reading restart file: ",A)' ) TRIM( dirname )
         !
         CALL qexml_init( iunpun )
         !
         CALL qexml_openfile( TRIM( dirname ) // TRIM( xmlpun ), &
                              'read', BINARY = .FALSE., IERR = ierr  )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot open restart file for reading', ierr )
      !
      s0 = cclock()
      !
      IF ( ionode ) THEN
         !
         qexml_version = " "
         !
         CALL qexml_read_header( FORMAT_VERSION = qexml_version, ierr = ierr )
         !
         qexml_version_init = .TRUE.
        
         !
         ! init logical variables for versioning
         !
         qexml_version_before_1_4_0 = .FALSE.
         !
         IF ( TRIM( version_compare( qexml_version, "1.4.0" )) == "older" ) &
            qexml_version_before_1_4_0 = .TRUE.
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'error reading the header', ierr )
      !
      CALL mp_bcast( qexml_version,               ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init,          ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_before_1_4_0 , ionode_id, intra_image_comm )
      !
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_status_cp( NFI=nfi,SIMTIME=simtime,TITLE=title, &
                                    FOUND=found, IERR=ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'error reading CP status', ierr )
      !
      IF ( ionode ) THEN
         !
         CALL qexml_closefile( 'read', IERR=ierr)
         !
      ENDIF
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // TRIM( xmlpun )
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), IERR = ierr )
         !
      END IF
      !
      !
      ! ... Read cell and positions
      !
      ALLOCATE( tau_( 3, nat ) )
      ALLOCATE( if_pos_( 3, nat ) )
      ALLOCATE( ityp_( nat ) )
      !
      IF ( ionode ) THEN
         !
         CALL read_cell( ibrav_, celldm_, alat_, a1_, a2_, a3_, b1, b2, b3 )
         !
         CALL recips( a1_, a2_, a3_, b1, b2, b3 )
         !
      END IF
      !
      IF ( ionode ) THEN
         !
         CALL qexml_read_ions( NSP = nsp_, NAT = nat_, ATM = atm_, ITYP = ityp_, &
                               PSFILE = psfile_,AMASS =  amass_, &
                               TAU = tau_, TAU_UNITS = pos_unit, IF_POS = if_pos_, IERR = ierr )
         !
         !
         !
         IF ( ierr == 0 ) THEN
            !
            IF( nsp_ /= nsp .OR. nat_ /= nat ) ierr = 2
            !
            DO i = 1, nat
               !
               IF ( ityp_(i) /= ityp(i) ) ierr = 3
               !
            END DO
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot read positions from restart file', ierr )
      !
      !  Read SPIN infos
      !
      lsda_ = ( nspin == 2 )
      !
      IF( ionode ) THEN
         !
         CALL qexml_read_spin( LSDA = lsda_, IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot read spins from restart file', ierr )
      !
      CALL mp_bcast( lsda_ , ionode_id, intra_image_comm )
      !
      IF( lsda_ .AND. nspin == 1 ) &
         CALL errore( 'cp_readfile', 'LSDA restart file with a spinless run', ierr )
      !
      !  Read Occupations infos
      !
      nstates_up_ = nupdwn( 1 )
      nstates_dw_ = nupdwn( 2 )
      !
      IF( ionode ) THEN
         !
         CALL qexml_read_occ( NSTATES_UP = nstates_up_, NSTATES_DW = nstates_dw_ , IERR = ierr)
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot read occupations from restart file', ierr )
      !
      CALL mp_bcast( nstates_up_ , ionode_id, intra_image_comm )
      CALL mp_bcast( nstates_dw_ , ionode_id, intra_image_comm )
      !
      IF( lsda_ ) THEN
         IF( ( nstates_up_ /=  nupdwn( 1 ) ) .OR. ( nstates_dw_ /=  nupdwn( 2 ) ) ) &
            CALL errore( 'cp_readfile', 'inconsistent number of spin states', ierr )
      END IF

      ! ... read MD timesteps variables
      !
      IF ( ionode ) &
         CALL iotk_scan_begin( iunpun, "TIMESTEPS", attr, FOUND = found )
      ! 
      ierr = 0
      ! 
      IF ( ionode .AND. found ) THEN
         !
         CALL iotk_scan_attr( attr, "nt", nt_ )
         !
         IF ( nt_ > 0 ) THEN
            !
            CALL iotk_scan_begin( iunpun, "STEP0" )
            !
            CALL iotk_scan_dat( iunpun, "ACCUMULATORS", acc )
            !
            CALL iotk_scan_begin( iunpun,"IONS_POSITIONS" )
            CALL iotk_scan_dat(   iunpun, "stau",  stau0(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "svel",  svel0(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "taui",  taui(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "cdmi",  cdmi(1:3) )
            CALL iotk_scan_dat(   iunpun, "force", force(1:3,1:nat) )
            CALL iotk_scan_end(   iunpun, "IONS_POSITIONS" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "nhpcl", nhpcl_ )
            CALL iotk_scan_dat(   iunpun, "nhpdim", nhpdim_ )
            !
            IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
               !
               CALL iotk_scan_dat( iunpun, "xnhp", xnhp0(1:nhpcl*nhpdim) )
               CALL iotk_scan_dat( iunpun, "vnhp", vnhp(1:nhpcl*nhpdim) )
               !
            ELSE
               !
               xnhp0(1:nhpcl*nhpdim) = 0.D0
               vnhp(1:nhpcl*nhpdim)  = 0.D0
               !
            END IF
            !
            CALL iotk_scan_end(   iunpun, "IONS_NOSE" )
            !
            CALL iotk_scan_dat( iunpun, "ekincm", ekincm )
            !
            CALL iotk_scan_begin( iunpun, "ELECTRONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhe", xnhe0 )
            CALL iotk_scan_dat(   iunpun, "vnhe", vnhe )
            CALL iotk_scan_end(   iunpun, "ELECTRONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
            CALL iotk_scan_dat(   iunpun, "ht",    ht )
            CALL iotk_scan_dat(   iunpun, "htvel", htvel )
            CALL iotk_scan_dat(   iunpun, "gvel",  gvel )
            CALL iotk_scan_end(   iunpun, "CELL_PARAMETERS" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhh", xnhh0 )
            CALL iotk_scan_dat(   iunpun, "vnhh", vnhh )
            CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
            !
            CALL iotk_scan_end( iunpun, "STEP0" )
            !
         ELSE
            !
            ierr = 40
            !
            GOTO 100
            !
         END IF
         !
         IF ( nt_ > 1 ) THEN
            !
            CALL iotk_scan_begin( iunpun, "STEPM" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_POSITIONS" )
            CALL iotk_scan_dat(   iunpun, "stau", staum(1:3,1:nat) )
            CALL iotk_scan_dat(   iunpun, "svel", svelm(1:3,1:nat) )
            CALL iotk_scan_end(   iunpun, "IONS_POSITIONS" )
            !
            CALL iotk_scan_begin( iunpun, "IONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "nhpcl", nhpcl_ )
            CALL iotk_scan_dat(   iunpun, "nhpdim", nhpdim_ )
            !
            IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
               !
               CALL iotk_scan_dat( iunpun, "xnhp",  xnhpm(1:nhpcl*nhpdim) )
               !
            ELSE
               !
               xnhpm(1:nhpcl*nhpdim) = 0.D0
               !
            END IF
            !
            CALL iotk_scan_end(   iunpun,"IONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "ELECTRONS_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhe", xnhem )
            CALL iotk_scan_end(   iunpun, "ELECTRONS_NOSE" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
            CALL iotk_scan_dat(   iunpun, "ht", htm )
            CALL iotk_scan_end(   iunpun, "CELL_PARAMETERS" )
            !
            CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
            CALL iotk_scan_dat(   iunpun, "xnhh", xnhhm )
            CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
            !
            CALL iotk_scan_end( iunpun, "STEPM" )
            !
         END IF
         !
         CALL iotk_scan_end( iunpun, "TIMESTEPS" )
         !
      ELSE IF ( ionode ) THEN
         !
         ! ... MD time steps not found, try to recover from CELL and POSITIONS
         ! 
         acc = 0.D0
         ! 
         ALLOCATE( isrt_( nat ) )
         !
         SELECT CASE( TRIM( pos_unit ) )
         CASE( "alat" )
            !
            tau_ = tau_ * alat_
            !
         CASE( "Angstrom" )
            !
            tau_ = tau_ * angstrom_au
            !
         CASE DEFAULT
            !
         END SELECT
         !
         CALL sort_tau( taui, isrt_ , tau_ , ityp_ , nat_ , nsp_ )
         ! 
         ht(1,:) = a1_
         ht(2,:) = a2_
         ht(3,:) = a3_
         !
         CALL invmat( 3, ht, htm1, omega )
         !
         hinv = TRANSPOSE( htm1 )
         !
         CALL r_to_s( taui, stau0, na, nsp, hinv )
         !
         CALL ions_cofmass( taui, amass_ , na, nsp, cdmi )
         !
         staum = stau0
         svel0 = 0.D0
         svelm = 0.D0
         force = 0.D0
         !
         htm   = ht
         htvel = 0.D0
         gvel  = 0.D0
         xnhh0 = 0.D0
         vnhh  = 0.D0
         xnhhm = 0.D0
         !
         xnhe0 = 0.D0
         xnhem = 0.D0
         vnhe  = 0.D0
         !
         ekincm = 0.D0
         !
         xnhp0 = 0.D0
         xnhpm = 0.D0
         vnhp  = 0.D0
         !
         DEALLOCATE( isrt_  )
         !
      END IF
      !
 100  CONTINUE
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF( ierr /= 0 ) THEN
         CALL mp_bcast( attr, ionode_id, intra_image_comm )
         CALL errore( 'cp_readfile ', TRIM( attr ), ierr )
      END IF
      !
      DEALLOCATE( tau_  )
      DEALLOCATE( if_pos_ )
      DEALLOCATE( ityp_ )
      !
      ! ... compute the scale factor
      !
      IF ( ionode ) CALL invmat( 3, ht, htm1, omega )
      !
      CALL mp_bcast( omega, ionode_id, intra_image_comm )
      !
      ! ... Beware: omega may be negative if axis are left-handed!
      !
      scalef = 1.D0 / SQRT( ABS( omega ) )
      !
      ! ... band Structure
      !
      IF ( ionode ) THEN
         !
         ierr = 0
         !
         CALL qexml_read_bands_info( NBND = NBND_TOT, NSPIN = nspin_, NELEC = nelec_, &
                                     NEL_UP = nel_(1), NEL_DOWN = nel_(2) , IERR = ierr)
      ENDIF

      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      CALL errore( 'cp_readfile ', 'error reading bands info', ierr )

      IF ( ionode) THEN
         !
         IF ( nspin_ /= nspin ) THEN
            attr = "spin do not match"
            ierr = 31
            GOTO 90
         END IF
         !
         IF ( nspin == 2 ) THEN
            !
            IF ( ( nel(1) /= nel_(1) ) .OR. ( nel(2) /= nel_(2) ) .OR. ( NINT( nelec_ ) /= nelt ) ) THEN
               attr = "electrons do not match"
               ierr = 33
               GOTO 90
            END IF
            !
         ELSE
            !
            IF ( NINT( nelec_ ) /= nelt ) THEN
               attr = "electrons do not match"
               ierr = 33
               GOTO 90
            END IF
            !
         END IF
         !
         nbnd_ = nbnd_tot
         !
         IF ( nbnd_ < nupdwn(1) ) THEN
            attr = "nbnd do not match"
            ierr = 32
            GOTO 90
         END IF
         !
      END IF
      !
 90   CONTINUE
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF( ierr /= 0 ) THEN
         CALL mp_bcast( attr, ionode_id, intra_image_comm )
         CALL errore( 'cp_readfile ', TRIM( attr ), ierr )
      END IF
      !
      IF( ionode ) THEN
         !
         CALL qexml_read_bands_cp( nk, nbnd_tot, nudx , nspin, iupdwn, &
      nupdwn, occ0, occm, ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile', &
                   'cannot read bands from restart file', ierr )
      !
      IF ( ionode ) THEN
         CALL iotk_scan_begin( iunpun, "EIGENVECTORS" )
      END IF
      !
      k_points_loop2: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         END IF
         !
         DO iss = 1, nspin
            !
            IF ( ionode ) THEN
               !
               CALL iotk_scan_begin( iunpun, "WFC0" // TRIM( iotk_index (iss) ), FOUND = found )
               !
               filename = "WFC0" // TRIM( iotk_index (iss) )
               !
               IF( .NOT. found ) THEN
                  !
                  IF( nspin == 2 ) THEN
                     CALL iotk_scan_begin( iunpun, "WFC" // TRIM( iotk_index (iss) ), FOUND = found )
                     filename = "WFC" // TRIM( iotk_index (iss) )
                  ELSE
                     CALL iotk_scan_begin( iunpun, "WFC", FOUND = found )
                     filename = "WFC"
                  END IF
                  !
               END IF
               !
            END IF
            !
            CALL mp_bcast( found, ionode_id, intra_image_comm )
            !
            IF ( iss == 1 ) THEN
               IF( found ) THEN
                  exist_wfc = .TRUE.
               END IF
            ELSE
               IF( exist_wfc .AND. .NOT. found ) THEN
                  CALL errore( " readfile ", " second spin component of wave functions not found! ", 1 )
               END IF
            END IF
            !
            IF( exist_wfc ) THEN
               !
               IF( .NOT. ( iss > 1 .AND. force_pairing ) ) THEN
                  !
                  ! Only WF with spin 1 are needed when force_pairing is active
                  !
                  ib = iupdwn(iss)
                  nb = nupdwn(iss)
                  !
                  ! filename is not needed we are following the link!
                  !
                  CALL read_wfc( iunpun, ik_eff , nk, kunit, iss_, nspin_, &
                                 c02( :, ib:ib+nb-1 ), ngwt_, nbnd_, ig_l2g, ngw, &
                                 filename, scalef_, &
                                 ionode, root_bgrp, intra_bgrp_comm, &
                                 inter_bgrp_comm, intra_pool_comm, .TRUE. )

                               !ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
                               !ionode, root_bgrp, intra_bgrp_comm, inter_bgrp_comm, intra_pool_comm )

                  !
               END IF
               !
               IF ( ionode ) &
                  CALL iotk_scan_end( iunpun, TRIM(filename) )
               !
               IF ( ionode ) THEN 
                  !
                  CALL iotk_scan_begin( iunpun, "WFCM" // TRIM( iotk_index (iss) ), FOUND = found )
                  !
                  filename = "WFCM" // TRIM( iotk_index (iss) )
                  !
               END IF
               !
               CALL mp_bcast( found, ionode_id, intra_image_comm )
               !
               IF( found ) THEN
                  !
                  IF( .NOT. ( iss > 1 .AND. force_pairing ) ) THEN
                     !
                     ! Only WF with spin 1 are needed when force_pairing is active
                     !
                     ib = iupdwn(iss)
                     nb = nupdwn(iss)
                     !
                     CALL read_wfc( iunpun, ik_eff, nk, kunit, iss_, nspin_, &
                                    cm2( :, ib:ib+nb-1 ), ngwt_, nbnd_, ig_l2g, ngw, &
                                    filename, scalef_ , &
                                    ionode, root_bgrp, intra_bgrp_comm, &
                                    inter_bgrp_comm, intra_pool_comm, .TRUE. )
                     !
                  END IF
                  !
                  IF ( ionode ) &
                     CALL iotk_scan_end( iunpun, TRIM( filename ) )
                  !
               ELSE
                  !
                  cm2 = c02
                  !
               END IF
               !
            END IF
            !
         END DO
         ! 
         !  here the I/O group send wfc to other groups
         !
         io_bgrp_id = 0
         IF( ionode ) io_bgrp_id = my_bgrp_id
         CALL mp_sum( io_bgrp_id, inter_bgrp_comm )
         CALL mp_sum( io_bgrp_id, intra_bgrp_comm )
         !
         CALL mp_bcast( cm2, io_bgrp_id, inter_bgrp_comm )
         CALL mp_bcast( c02, io_bgrp_id, inter_bgrp_comm )
         !
         DO iss = 1, nspin
            !
            ! ... read matrix lambda to file
            !
            cspin = iotk_index( iss )
            !
            ALLOCATE( mrepl( nudx, nudx ) )
            !
            IF( ionode ) THEN
               CALL iotk_scan_dat( iunpun, "LAMBDA0" // TRIM( cspin ), mrepl, FOUND = found )
               IF( .NOT. found ) THEN
                  WRITE( stdout, * ) 'WARNING lambda0 not read from restart file'
                  mrepl = 0.0d0
               END IF
            END IF

            CALL mp_bcast( mrepl, ionode_id, intra_image_comm )

            CALL distribute_lambda( mrepl, lambda0(:,:,iss), descla(iss) )

            IF( ionode ) THEN
               CALL iotk_scan_dat( iunpun, "LAMBDAM" // TRIM( cspin ), mrepl, FOUND = found )
               IF( .NOT. found ) THEN
                  WRITE( stdout, * ) 'WARNING lambdam not read from restart file'
                  mrepl = 0.0d0
               END IF
            END IF
            ! 
            CALL mp_bcast( mrepl, ionode_id, intra_image_comm )

            CALL distribute_lambda( mrepl, lambdam(:,:,iss), descla(iss) )
            !
            IF ( PRESENT( mat_z ) ) THEN
               !
               IF( ionode ) THEN
                  CALL iotk_scan_dat( iunpun, "MAT_Z" // TRIM( iotk_index( iss ) ), mrepl, FOUND = found )
                  IF( .NOT. found ) THEN
                     WRITE( stdout, * ) 'WARNING mat_z not read from restart file'
                     mrepl = 0.0d0
                  END IF
               END IF

               CALL mp_bcast( mrepl, ionode_id, intra_image_comm )

               CALL distribute_zmat( mrepl, mat_z(:,:,iss), descla(iss) )
               !
            END IF
            !
            DEALLOCATE( mrepl )
            !
         END DO
         !
         IF ( ionode ) CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop2
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( iunpun, "EIGENVECTORS" )
         !
      END IF
      !
      !-------------------------------------------------------------------------
      ! BS : Wannier centers
      IF ( ionode ) THEN
         !
         IF ( lwf ) THEN
            !
            CALL iotk_scan_begin( iunpun, "WANNIER_CENTERS" , FOUND=found)
            !
            IF (found) THEN
               !
               DO iss = 1, nspin
                  !
                  CALL iotk_scan_dat(   iunpun, "wanniercentres" // TRIM(iotk_index(iss)), &
                                                   &wfc(1:3, iupdwn(iss):iupdwn(iss) + nupdwn(iss) -1 ) )
                  !
               ENDDO
               !
            ELSE
               WRITE( stdout, * ) 'WARNING wannier centers not read from restart file:'
            ENDIF
            !
            CALL iotk_scan_end(   iunpun, "WANNIER_CENTERS" )
         END IF
         !
      END IF
      ! BS : Wannier centers
      !-------------------------------------------------------------------------
      !
      CALL mp_bcast( qexml_version,      ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( nfi,     ionode_id, intra_image_comm )
      CALL mp_bcast( simtime, ionode_id, intra_image_comm )
      CALL mp_bcast( title,   ionode_id, intra_image_comm )
      CALL mp_bcast( acc,     ionode_id, intra_image_comm )
      !
      CALL mp_bcast( ht,    ionode_id, intra_image_comm )
      CALL mp_bcast( htm,   ionode_id, intra_image_comm )
      CALL mp_bcast( htvel, ionode_id, intra_image_comm )
      CALL mp_bcast( gvel,  ionode_id, intra_image_comm )
      CALL mp_bcast( xnhh0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhhm, ionode_id, intra_image_comm )
      CALL mp_bcast( vnhh,  ionode_id, intra_image_comm )
      CALL mp_bcast( b1,    ionode_id, intra_image_comm )
      CALL mp_bcast( b2,    ionode_id, intra_image_comm )
      CALL mp_bcast( b3,    ionode_id, intra_image_comm )
      !
      CALL mp_bcast( stau0, ionode_id, intra_image_comm )
      CALL mp_bcast( svel0, ionode_id, intra_image_comm )
      CALL mp_bcast( staum, ionode_id, intra_image_comm )
      CALL mp_bcast( svelm, ionode_id, intra_image_comm )
      CALL mp_bcast( taui,  ionode_id, intra_image_comm )
      CALL mp_bcast( force, ionode_id, intra_image_comm )
      CALL mp_bcast( cdmi,  ionode_id, intra_image_comm )
      CALL mp_bcast( xnhp0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhpm, ionode_id, intra_image_comm ) 
      CALL mp_bcast( vnhp,  ionode_id, intra_image_comm )
      !
      CALL mp_bcast( xnhe0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhem, ionode_id, intra_image_comm )
      CALL mp_bcast( vnhe,  ionode_id, intra_image_comm )
      !
      CALL mp_bcast( kunit, ionode_id, intra_image_comm )

      CALL mp_bcast( occ0, ionode_id, intra_image_comm )
      CALL mp_bcast( occm, ionode_id, intra_image_comm )
      !
      !-------------------------------------------------------------------------
      ! BS : Wannier centers
      IF ( lwf ) CALL mp_bcast( wfc, ionode_id, intra_image_comm )
      !-------------------------------------------------------------------------
      !
      IF ( PRESENT( mat_z ) ) &
         CALL mp_bcast( mat_z(:,:,:), ionode_id, intra_image_comm )
      !
      IF ( ionode ) &
         CALL iotk_close_read( iunpun )
      !
      IF( .NOT. exist_wfc ) THEN
         !
         IF( my_bgrp_id == 0 ) THEN
            !
            tmp_dir_save = tmp_dir
            tmp_dir = TRIM( qexml_kpoint_dirname( tmp_dir, 1 ) )
            !
            iunwfc = 10
            nwordwfc = SIZE( c02 )
            !
            CALL diropn ( iunwfc, 'wfc', 2*nwordwfc, exst )
     
            IF ( exst ) THEN
               CALL davcio ( c02, 2*nwordwfc, iunwfc, 1, -1 )  ! read wave funct
               CALL davcio ( cm2, 2*nwordwfc, iunwfc, 2, -1 )  ! read wave funct
               CLOSE( UNIT = iunwfc, STATUS = 'KEEP' )
            ELSE
               CLOSE( UNIT = iunwfc, STATUS = 'DELETE' )
               CALL errore( ' readfile ' , ' no wave function found! ' , 1 )
            END IF

            tmp_dir = tmp_dir_save
            !
         END IF

         CALL mp_bcast( c02, 0, inter_bgrp_comm )
         CALL mp_bcast( cm2, 0, inter_bgrp_comm )

      END IF
      !
      s1 = cclock()
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, &
                '(3X,"restart file read in ",F8.3," sec.",/)' )  ( s1 - s0 )
         !
      END IF
      !
      if ( nprint_nfi == -2 ) then
         CALL read_print_counter( nprint_nfi, tmp_dir, ndr )
         IF( iverbosity > 1 ) write( stdout,*) 'nprint_nfi= ',nprint_nfi
      endif
#endif
      !
      RETURN
      !
    END SUBROUTINE cp_readfile
    ! 
    !------------------------------------------------------------------------
    SUBROUTINE cp_read_wfc( ndr, tmp_dir, ik, nk, iss, nspin, c2, tag )
      !------------------------------------------------------------------------
      !
      USE electrons_base,     ONLY : iupdwn, nupdwn
      USE gvecw,              ONLY : ngw
      USE io_global,          ONLY : ionode
      USE mp_global,          ONLY : root_bgrp, intra_bgrp_comm, inter_bgrp_comm, intra_pool_comm, my_bgrp_id
      USE mp,                 ONLY : mp_bcast, mp_sum
      USE gvect,              ONLY : ig_l2g
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)  :: ndr
      CHARACTER(LEN=*),      INTENT(IN)  :: tmp_dir
      INTEGER,               INTENT(IN)  :: ik, iss, nk, nspin
      CHARACTER,             INTENT(IN)  :: tag
      COMPLEX(DP), OPTIONAL, INTENT(OUT) :: c2(:,:)
      !
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: ik_eff, ib, nb, kunit, iss_, nspin_, ngwt_, nbnd_
      INTEGER            :: io_bgrp_id
      REAL(DP)           :: scalef
      !
#if !defined (__OLDXML)
      CALL errore('cp_read_wfc','called in the wrong case',1)
#else
      kunit = 1
      !
      ik_eff = ik + ( iss - 1 ) * nk
      !
      dirname = qexml_restart_dirname( tmp_dir, prefix, ndr )
      !
      IF ( tag /= 'm' ) THEN
         !
         IF ( nspin == 1 ) THEN
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'evc0', ik ) )
            !
         ELSE
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'evc0', ik, iss ) )
            !
         END IF
         !
      ELSE
         !
         IF ( nspin == 1 ) THEN
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'evcm', ik ) )
            !
         ELSE
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'evcm', ik, iss ) )
            !
         END IF
         !
      END IF
      !
      ib = iupdwn(iss)
      nb = nupdwn(iss)
      !
      CALL read_wfc( iunout, ik_eff, nk, kunit, iss_, nspin_, &
                     c2(:,ib:ib+nb-1), ngwt_, nbnd_, ig_l2g, ngw,  &
                     filename, scalef, & 
                     ionode, root_bgrp, intra_bgrp_comm, &
                     inter_bgrp_comm, intra_pool_comm )
      ! 
      !  here the I/O group send wfc to other groups
      !
      io_bgrp_id = 0
      IF( ionode ) io_bgrp_id = my_bgrp_id
      CALL mp_sum( io_bgrp_id, inter_bgrp_comm )
      CALL mp_sum( io_bgrp_id, intra_bgrp_comm )
      CALL mp_bcast( c2, io_bgrp_id, inter_bgrp_comm )
#endif
      !
      RETURN
      !
    END SUBROUTINE cp_read_wfc
    !
!==============================================================================
!Modified from cp_read_wfc to read valence states for nscf calculations
!Lingzhu Kong
    !------------------------------------------------------------------------
    SUBROUTINE cp_read_wfc_Kong( ndr, tmp_dir, ik, nk, iss, nspin, c2, tag )
      !------------------------------------------------------------------------
      !
      USE kinds,              ONLY : DP
      USE gvecw,              ONLY : ngw
      USE gvect,              ONLY : ig_l2g
      USE wannier_base,       ONLY : vnbsp        
      USE mp_global,          ONLY : intra_bgrp_comm, inter_bgrp_comm, &
                                     root_bgrp, intra_pool_comm
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)  :: ndr
      CHARACTER(LEN=*),      INTENT(IN)  :: tmp_dir
      INTEGER,               INTENT(IN)  :: ik, iss, nk, nspin
      CHARACTER,             INTENT(IN)  :: tag
      COMPLEX(DP), OPTIONAL, INTENT(OUT) :: c2(:,:)
      !
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: ik_eff, ib, nb, kunit, iss_, nspin_, ngwt_, nbnd_
      REAL(DP)           :: scalef
      !
#if !defined (__OLDXML)
      CALL errore('cp_read_wfc_Kong','XSD under development',1)
#else
      kunit = 1
      !
      ik_eff = ik + ( iss - 1 ) * nk
      !
      dirname = qexml_restart_dirname( tmp_dir, prefix, ndr )
      !
      IF ( tag /= 'm' ) THEN
         !
         IF ( nspin == 1 ) THEN
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'evc0', ik ) )
            !
         ELSE
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'evc0', ik, iss ) )
            !
         END IF
         !
      ELSE
         !
         IF ( nspin == 1 ) THEN
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'evcm', ik ) )
            !
         ELSE
            !
            filename = TRIM( qexml_wfc_filename( dirname, 'evcm', ik, iss ) )
            !
         END IF
         !
      END IF
      !
      ib = 1
      nb = vnbsp             
      !
      print *,'before read_wfc me'
      CALL read_wfc( iunout, ik_eff, nk, kunit, iss_, nspin_, &
                     c2(:,ib:ib+nb-1), ngwt_, nbnd_, ig_l2g, ngw,  &
                     filename, scalef, &
                     ionode, root_bgrp, intra_bgrp_comm, &
                     inter_bgrp_comm, intra_pool_comm )

      !
#endif
      RETURN
      !
    END SUBROUTINE cp_read_wfc_Kong
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_read_cell( ndr, tmp_dir, ascii, ht, &
                             htm, htvel, gvel, xnhh0, xnhhm, vnhh )
      !------------------------------------------------------------------------
      !
      USE io_files,  ONLY : iunpun, xmlpun
      USE mp_global, ONLY : intra_image_comm
      USE mp,        ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN)    :: ndr
      CHARACTER(LEN=*), INTENT(IN)    :: tmp_dir
      LOGICAL,          INTENT(IN)    :: ascii
      REAL(DP),         INTENT(INOUT) :: ht(3,3)
      REAL(DP),         INTENT(INOUT) :: htm(3,3)
      REAL(DP),         INTENT(INOUT) :: htvel(3,3)
      REAL(DP),         INTENT(INOUT) :: gvel(3,3)
      REAL(DP),         INTENT(INOUT) :: xnhh0(3,3)
      REAL(DP),         INTENT(INOUT) :: xnhhm(3,3)
      REAL(DP),         INTENT(INOUT) :: vnhh(3,3)
      !
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: strlen
      INTEGER            :: i, ierr, nt_
      LOGICAL            :: found
      !
      ! ... variables read for testing pourposes
      !
      INTEGER          :: ibrav_
      REAL(DP)         :: alat_
      REAL(DP)         :: celldm_(6)
      REAL(DP)         :: a1_(3), a2_(3), a3_(3)
      REAL(DP)         :: b1_(3), b2_(3), b3_(3)
      CHARACTER(iotk_attlenx)  :: attr
      !
#if !defined (__OLDXML)
      CALL errore('cp_read_cell','XSD under development',1)
#else
      !
      dirname = qexml_restart_dirname( tmp_dir, prefix, ndr ) 
      !
      filename = TRIM( dirname ) // TRIM( xmlpun )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                              BINARY = .FALSE., ROOT = attr, IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_read_cell', &
                   'cannot open restart file for reading: ' // TRIM(filename), &
                   ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "TIMESTEPS", attr, FOUND = found )
         !
         IF ( found ) THEN
            !
            CALL iotk_scan_attr( attr, "nt", nt_ )
            !
            IF ( nt_ > 0 ) THEN
               !
               CALL iotk_scan_begin( iunpun, "STEP0" )
               !
               CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
               CALL iotk_scan_dat(   iunpun, "ht",    ht )
               CALL iotk_scan_dat(   iunpun, "htvel", htvel )
               CALL iotk_scan_dat(   iunpun, "gvel",  gvel, &
                                     FOUND = found, IERR = ierr )
               !
               IF ( .NOT. found ) gvel = 0.D0
               !
               CALL iotk_scan_end( iunpun, "CELL_PARAMETERS" )
               !
               CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
               CALL iotk_scan_dat(   iunpun, "xnhh", xnhh0 )
               CALL iotk_scan_dat(   iunpun, "vnhh", vnhh )
               CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
               !
               CALL iotk_scan_end( iunpun, "STEP0" )
               !
            ELSE
               !
               ierr = 40
               !
               GOTO 100
               !
            END IF
            !
            IF( nt_ > 1 ) THEN
               !
               CALL iotk_scan_begin(iunpun,"STEPM")
               !
               CALL iotk_scan_begin( iunpun, "CELL_PARAMETERS" )
               CALL iotk_scan_dat(   iunpun, "ht", htm)
               CALL iotk_scan_end(   iunpun, "CELL_PARAMETERS" )
               !
               CALL iotk_scan_begin( iunpun, "CELL_NOSE" )
               CALL iotk_scan_dat(   iunpun, "xnhh", xnhhm )
               CALL iotk_scan_end(   iunpun, "CELL_NOSE" )
               !
               CALL iotk_scan_end( iunpun, "STEPM" )
               !
            END IF
            !
            CALL iotk_scan_end( iunpun, "TIMESTEPS" )
            !
         ELSE
            !
            ! ... MD steps have not been found, try to restart from cell data
            !
            CALL read_cell( ibrav_, celldm_, alat_, a1_,a2_,a3_, b1_, b2_, b3_ )
            !
            ht(1,:) = a1_
            ht(2,:) = a2_
            ht(3,:) = a3_
            !
            htm   = ht
            htvel = 0.D0
            gvel  = 0.D0
            xnhh0 = 0.D0
            vnhh  = 0.D0
            xnhhm = 0.D0
            !
         END IF
         !
      END IF
      !
 100  CONTINUE
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      CALL mp_bcast( attr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_read_cell ', attr, ierr )
      !
      CALL mp_bcast( ht,    ionode_id, intra_image_comm )
      CALL mp_bcast( htm,   ionode_id, intra_image_comm )
      CALL mp_bcast( htvel, ionode_id, intra_image_comm )
      CALL mp_bcast( gvel,  ionode_id, intra_image_comm )
      CALL mp_bcast( xnhh0, ionode_id, intra_image_comm )
      CALL mp_bcast( xnhhm, ionode_id, intra_image_comm )
      CALL mp_bcast( vnhh,  ionode_id, intra_image_comm )
      !
      IF ( ionode ) CALL iotk_close_read( iunpun )
#endif
      !
      RETURN
      !
    END SUBROUTINE cp_read_cell
    !
#if defined (__OLDXML)
    !----------------------------------------------------------------------------
    SUBROUTINE read_cell( ibrav, celldm, alat, a1, a2, a3, b1, b2, b3 )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(OUT) :: ibrav
      REAL(DP),         INTENT(OUT) :: celldm(6), alat
      REAL(DP),         INTENT(OUT) :: a1(3), a2(3), a3(3)
      REAL(DP),         INTENT(OUT) :: b1(3), b2(3), b3(3)
      !
      CHARACTER(LEN=256) :: bravais_lattice
      !
      !
      CALL iotk_scan_begin( iunpun, "CELL" )
      !
      CALL iotk_scan_dat( iunpun, "BRAVAIS_LATTICE", bravais_lattice )
      !
      SELECT CASE ( TRIM( bravais_lattice ) )
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
      END SELECT
      !
      CALL iotk_scan_dat( iunpun, "LATTICE_PARAMETER", alat )
      CALL iotk_scan_dat( iunpun, "CELL_DIMENSIONS", celldm(1:6) )
      !
      CALL iotk_scan_begin( iunpun, "DIRECT_LATTICE_VECTORS" )
      CALL iotk_scan_dat(   iunpun, "a1", a1 )
      CALL iotk_scan_dat(   iunpun, "a2", a2 )
      CALL iotk_scan_dat(   iunpun, "a3", a3 )
      CALL iotk_scan_end(   iunpun, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_scan_begin( iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      CALL iotk_scan_dat(   iunpun, "b1", b1 )
      CALL iotk_scan_dat(   iunpun, "b2", b2 )
      CALL iotk_scan_dat(   iunpun, "b3", b3 )
      CALL iotk_scan_end(   iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      !
      CALL iotk_scan_end( iunpun, "CELL" )
      !
      RETURN
      !
    END SUBROUTINE
    !
    !
    !----------------------------------------------------------------------------
    SUBROUTINE write_gk( iun, ik, filename )
      !----------------------------------------------------------------------------
       !
       USE gvecw,                    ONLY : ngw, ngw_g
       USE gvect,                    ONLY : ngm, ngm_g
       USE control_flags,            ONLY : gamma_only
       USE gvect,                    ONLY : ig_l2g, mill
       USE mp,                       ONLY : mp_sum
       USE mp_global,                ONLY : intra_bgrp_comm
       USE io_global,                ONLY : ionode
       !
       IMPLICIT NONE
       !
       INTEGER,            INTENT(IN) :: iun, ik
       CHARACTER(LEN=256), INTENT(IN) :: filename
       !
       INTEGER, ALLOCATABLE :: igwk(:)
       INTEGER, ALLOCATABLE :: itmp1(:)
       INTEGER, ALLOCATABLE :: mill_g(:,:)
       INTEGER  :: npwx_g, npw_g, ig, ngg
       REAL(DP) :: xk(3)
       CHARACTER(iotk_attlenx)  :: attr

       ! ... Collect G vectors
       !   
       ALLOCATE( mill_g( 3, ngm_g ) )
       !
       mill_g = 0
       !
       mill_g(:,ig_l2g(1:ngm)) = mill(:,1:ngm)
       !
       CALL mp_sum( mill_g, intra_bgrp_comm )
       !

       xk     = 0.0d0
       npwx_g = ngw_g
       npw_g  = ngw_g

       ALLOCATE( igwk( npwx_g ) )
       ! 
       igwk = 0
       !
       ALLOCATE( itmp1( npw_g ) )
       !
       itmp1 = 0
       !
       !
       DO ig = 1, ngw
          !
          itmp1( ig_l2g( ig ) ) = ig_l2g( ig )
          !
       END DO
       !
       CALL mp_sum( itmp1, intra_bgrp_comm )
       !
       ngg = 0
       !
       DO ig = 1, npw_g
          !
          IF ( itmp1(ig) == ig ) THEN
             !
             ngg = ngg + 1
             !
             igwk( ngg ) = ig
             !
          END IF
          !
       END DO

       DEALLOCATE( itmp1 )
       !
       IF ( ionode ) THEN
          !
          CALL iotk_open_write( iun, FILE = TRIM( filename ), &
                                ROOT="GK-VECTORS", BINARY = .TRUE. )
          !
          CALL iotk_write_dat( iun, "NUMBER_OF_GK-VECTORS", npw_g )
          CALL iotk_write_dat( iun, "MAX_NUMBER_OF_GK-VECTORS", npwx_g )
          CALL iotk_write_dat( iun, "GAMMA_ONLY", gamma_only )
          !
          CALL iotk_write_attr ( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
          CALL iotk_write_dat( iun, "K-POINT_COORDS", xk(:), ATTR = attr )
          !
          CALL iotk_write_dat( iun, "INDEX", igwk( 1:npw_g ) )
          CALL iotk_write_dat( iun, "GRID", mill_g(1:3, igwk(1:npw_g)), COLUMNS = 3 )
          !
          CALL iotk_close_write( iun )
          !
       END IF
       !
       DEALLOCATE( igwk )
       DEALLOCATE( mill_g )

       RETURN

    END SUBROUTINE write_gk
#endif
    !
END MODULE cp_restart
