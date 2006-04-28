!
! Copyright (C) 2005 Quantum-ESPRESSO group
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
  USE iotk_module
  USE xml_io_base
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE io_files,  ONLY : prefix, iunpun, xmlpun
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE :: read_cell
  !
  INTEGER, PRIVATE :: iunout
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_writefile( ndw, scradir, ascii, nfi, simtime, acc, nk, xk, &
                             wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,   &
                             taui, cdmi, stau0, svel0, staum, svelm, force,  &
                             vnhp, xnhp0, xnhpm,nhpcl,nhpdim, occ0, occm,    &
                             lambda0,lambdam, xnhe0, xnhem, vnhe, ekincm,    &
                             et, rho, c04, cm4, c02, cm2, mat_z )
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : gamma_only, force_pairing, reduce_io
      USE io_files,                 ONLY : psfile, pseudo_dir
      USE mp_global,                ONLY : intra_image_comm, me_image, nproc_image
      USE printout_base,            ONLY : title
      USE grid_dimensions,          ONLY : nr1, nr2, nr3, nr1x, nr2x, nr3l
      USE smooth_grid_dimensions,   ONLY : nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
      USE gvecp,                    ONLY : ngm, ngmt, ecutp, gcutp
      USE gvecs,                    ONLY : ngs, ngst, ecuts, gcuts, dual
      USE gvecw,                    ONLY : ngw, ngwt, ecutw, gcutw
      USE reciprocal_vectors,       ONLY : ig_l2g, mill_l
      USE electrons_base,           ONLY : nspin, nbnd, nbsp, nelt, nel, &
                                           nupdwn, iupdwn, f, nudx
      USE cell_base,                ONLY : ibrav, alat, celldm, &
                                           symm_type, s_to_r
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           pmass, amass, iforce, ind_bck
      USE funct,                    ONLY : get_dft_name
      USE energies,                 ONLY : enthal, ekin, eht, esr, eself, &
                                           epseu, enl, exc, vave
      USE mp_global,                ONLY : nproc, mpime
      USE mp,                       ONLY : mp_sum
      USE parameters,               ONLY : nhclm
      USE fft_base,                 ONLY : dfftp
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN) :: ndw          !
      CHARACTER(LEN=*),      INTENT(IN) :: scradir      !
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
      REAL(DP),              INTENT(IN) :: occ0(:,:,:)  ! 
      REAL(DP),              INTENT(IN) :: occm(:,:,:)  ! 
      REAL(DP),              INTENT(IN) :: lambda0(:,:,:) ! 
      REAL(DP),              INTENT(IN) :: lambdam(:,:,:) ! 
      REAL(DP),              INTENT(IN) :: xnhe0        ! 
      REAL(DP),              INTENT(IN) :: xnhem        ! 
      REAL(DP),              INTENT(IN) :: vnhe         ! 
      REAL(DP),              INTENT(IN) :: ekincm       ! 
      REAL(DP),              INTENT(IN) :: et(:,:,:)    ! 
      REAL(DP),              INTENT(IN) :: rho(:,:)     ! 
      COMPLEX(DP), OPTIONAL, INTENT(IN) :: c04(:,:,:,:) ! 
      COMPLEX(DP), OPTIONAL, INTENT(IN) :: cm4(:,:,:,:) ! 
      COMPLEX(DP), OPTIONAL, INTENT(IN) :: c02(:,:)     ! 
      COMPLEX(DP), OPTIONAL, INTENT(IN) :: cm2(:,:)     ! 
      REAL(DP),    OPTIONAL, INTENT(IN) :: mat_z(:,:,:) ! 
      !
      LOGICAL               :: write_charge_density
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname, filename, rho_file_base
      CHARACTER(LEN=4)      :: cspin
      INTEGER               :: kunit, ib, ik_eff
      INTEGER               :: k1, k2, k3
      INTEGER               :: nk1, nk2, nk3
      INTEGER               :: j, i, ispin, ig, nspin_wfc, iss_wfc
      INTEGER               :: is, ia, isa, iss, ise, ik, ierr
      INTEGER,  ALLOCATABLE :: mill(:,:)
      INTEGER,  ALLOCATABLE :: ftmp(:,:)
      INTEGER,  ALLOCATABLE :: ityp(:)
      REAL(DP), ALLOCATABLE :: tau(:,:)
      REAL(DP), ALLOCATABLE :: rhosum(:)
      REAL(DP)              :: omega, htm1(3,3), h(3,3)
      REAL(DP)              :: a1(3), a2(3), a3(3)
      REAL(DP)              :: b1(3), b2(3), b3(3)
      REAL(DP)              :: nelec
      REAL(DP)              :: scalef
      LOGICAL               :: lsda
      REAL(DP)              :: s0, s1, cclock
      !
      !
      write_charge_density = .NOT.reduce_io
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
      dirname = restart_dir( scradir, ndw )
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
         CALL create_directory( kpoint_dir( dirname, i ) )
         !
      END DO
      !
      ! ... Some ( CP/FPMD ) default values
      !
      IF ( nspin == 2 ) THEN
         !
         kunit = 2
         !
      ELSE
         !
         kunit = 1
         !
      END IF
      !
      k1  = 0
      k2  = 0
      k3  = 0
      nk1 = 1
      nk2 = 1
      nk3 = 1
      !
      ! ... Compute Cell related variables
      !
      h = TRANSPOSE( ht )
      !
      CALL invmat( 3, ht, htm1, omega )
      !
      a1 = ht(1,:)
      a2 = ht(2,:)
      a3 = ht(3,:)
      !
      CALL recips( a1, a2, a3, b1, b2, b3 )
      !
      scalef = 1.D0 / SQRT( omega )
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
      CALL s_to_r( stau0, tau, na, nsp, h )
      !   
      ! ... Collect G vectors
      !   
      ALLOCATE( mill( 3, ngmt ) )
      !
      mill = 0
      !
      mill(:,ig_l2g(1:ngm)) = mill_l(:,1:ngm)
      !
      CALL mp_sum( mill, intra_image_comm )
      !
      lsda = ( nspin == 2 )
      !
      ALLOCATE( ftmp( nudx, nspin ) )
      !
      DO ispin = 1, nspin
         !
         iss = iupdwn(ispin)
         !
         ise = iss + nupdwn(ispin) - 1
         !
         ftmp(1:nupdwn(ispin),ispin) = f(iss:ise)
         !
      END DO
      !
      IF ( ionode ) THEN
         !
         ! ... Open XML descriptor
         !
         WRITE( stdout, '(/,3X,"writing restart file: ",A)' ) TRIM( dirname )
         !
         CALL iotk_open_write( iunpun, FILE = TRIM( dirname ) // '/' // &
                             & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_writefile ', &
                   'cannot open restart file for writing', ierr )
      !
      s0 = cclock() 
      !
      IF ( ionode ) THEN
         !
!-------------------------------------------------------------------------------
! ... STATUS
!-------------------------------------------------------------------------------
         !
         CALL iotk_write_begin( iunpun, "STATUS" )
         !
         CALL iotk_write_attr( attr, "ITERATION", nfi, FIRST = .TRUE. )
         CALL iotk_write_empty(iunpun, "STEP", attr )
         !
         CALL iotk_write_attr( attr, "UNITS", "pico-seconds", FIRST = .TRUE. ) 
         CALL iotk_write_dat( iunpun, "TIME", simtime, ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "TITLE", TRIM( title ) )
         !
         CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
         CALL iotk_write_dat( iunpun, "KINETIC_ENERGY", ekin,   ATTR = attr )
         CALL iotk_write_dat( iunpun, "HARTREE_ENERGY", eht,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "EWALD_TERM",     esr,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "GAUSS_SELFINT",  eself,  ATTR = attr )
         CALL iotk_write_dat( iunpun, "LPSP_ENERGY",    epseu,  ATTR = attr )
         CALL iotk_write_dat( iunpun, "NLPSP_ENERGY",   enl,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "EXC_ENERGY",     exc,    ATTR = attr )
         CALL iotk_write_dat( iunpun, "AVERAGE_POT",    vave,   ATTR = attr )
         CALL iotk_write_dat( iunpun, "ENTHALPY",       enthal, ATTR = attr )
         !
         CALL iotk_write_end( iunpun, "STATUS" )      
         !
!-------------------------------------------------------------------------------
! ... CELL
!-------------------------------------------------------------------------------
         !
         CALL write_cell( ibrav, symm_type, &
                          celldm, alat, a1, a2, a3, b1, b2, b3 )
         !
!-------------------------------------------------------------------------------
! ... IONS
!-------------------------------------------------------------------------------
         !
         CALL write_ions( nsp, nat, atm, ityp(ind_bck(:)), &
                          psfile, pseudo_dir, amass, tau(:,ind_bck(:)), &
                          iforce(:,ind_bck(:)), dirname, 1.D0 )
         !
!-------------------------------------------------------------------------------
! ... PLANE_WAVES
!-------------------------------------------------------------------------------
         !
         CALL write_planewaves( ecutw, dual, ngwt, gamma_only, nr1, nr2, &
                                nr3, ngmt, nr1s, nr2s, nr3s, ngst, nr1b, &
                                nr2b, nr3b, mill, .FALSE. )
         !
!-------------------------------------------------------------------------------
! ... SPIN
!-------------------------------------------------------------------------------
         !
         CALL write_spin( lsda, .FALSE., 1, .FALSE., .TRUE. )
         !
!-------------------------------------------------------------------------------
! ... EXCHANGE_CORRELATION
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         CALL write_xc( DFT = dft_name, NSP = nsp, LDA_PLUS_U = .FALSE. )
         !
!-------------------------------------------------------------------------------
! ... OCCUPATIONS
!-------------------------------------------------------------------------------
         !
         CALL write_occ( LGAUSS = .FALSE., LTETRA = .FALSE., &
                         TFIXED_OCC = .TRUE., LSDA = lsda, NELUP = nupdwn(1), &
                         NELDW = nupdwn(2), F_INP = DBLE( ftmp ) )
         !
!-------------------------------------------------------------------------------
! ... BRILLOUIN_ZONE
!-------------------------------------------------------------------------------
         !
         CALL write_bz( nk, xk, wk, k1, k2, k3, nk1, nk2, nk3 )
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
         IF (write_charge_density) CALL iotk_write_begin( iunpun, "CHARGE-DENSITY" )
         !
      END IF
      !
      IF (write_charge_density) then
         !
         rho_file_base = 'charge-density'
         !
         IF ( ionode )&
              CALL iotk_link( iunpun, "RHO", rho_file_base, &
              CREATE = .FALSE., BINARY = .FALSE. )
         !
         rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
         !
         IF ( nspin == 1 ) THEN
            !
            CALL write_rho_xml( rho_file_base, &
                                rho(:,1), nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
         ELSE IF ( nspin == 2 ) THEN
            !
            ALLOCATE( rhosum( SIZE( rho, 1 ) ) )
            !
            rhosum = rho(:,1) + rho(:,2) 
            !
            CALL write_rho_xml( rho_file_base, &
                                rhosum, nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
            DEALLOCATE( rhosum )
            !
            rho_file_base = 'charge-density-up'
            !
            IF ( ionode ) &
                 CALL iotk_link( iunpun, "RHO_UP", rho_file_base, &
                 CREATE = .FALSE., BINARY = .FALSE. )
            !
            rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
            !
            CALL write_rho_xml( rho_file_base, &
                                rho(:,1), nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
            rho_file_base = 'charge-density-dw'
            !
            IF ( ionode ) &
                 CALL iotk_link( iunpun, "RHO_DW", rho_file_base, &
                 CREATE = .FALSE., BINARY = .FALSE. )
            !
            rho_file_base = TRIM( dirname ) // '/' // TRIM( rho_file_base )
            !
            CALL write_rho_xml( rho_file_base, &
                                rho(:,2), nr1, nr2, nr3, nr1x, nr2x, dfftp%ipp, dfftp%npp )
            !
         END IF
         !
      END IF ! write_charge_density
      !
      IF ( ionode ) THEN
         !
         if (write_charge_density) CALL iotk_write_end( iunpun, "CHARGE-DENSITY" )
         !
!-------------------------------------------------------------------------------
! ... TIMESTEPS
!-------------------------------------------------------------------------------
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
         CALL iotk_write_dat(   iunpun, "stau",  stau0(1:3,1:nat) )
         CALL iotk_write_dat(   iunpun, "svel",  svel0(1:3,1:nat) )
         CALL iotk_write_dat(   iunpun, "taui",  taui(1:3,1:nat) )
         CALL iotk_write_dat(   iunpun, "cdmi",  cdmi(1:3) )
         CALL iotk_write_dat(   iunpun, "force", force(1:3,1:nat) )
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
         CALL iotk_write_dat(   iunpun, "stau", staum(1:3,1:nat) )
         CALL iotk_write_dat(   iunpun, "svel", svelm(1:3,1:nat) )
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
!-------------------------------------------------------------------------------
! ... BAND_STRUCTURE
!-------------------------------------------------------------------------------
         ! 
         CALL iotk_write_begin( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin )
         !
         nelec = nelt
         !
         IF ( nspin == 2 ) THEN
            !
            CALL iotk_write_attr( attr, "UP", nel(1), FIRST = .TRUE. )
            CALL iotk_write_attr( attr, "DW", nel(2) )
            CALL iotk_write_dat( iunpun, &
                                 "NUMBER_OF_ELECTRONS", nelec, ATTR = attr )
            !
            CALL iotk_write_attr( attr, "UP", nupdwn(1), FIRST = .TRUE. )
            CALL iotk_write_attr( attr, "DW", nupdwn(2) )
            CALL iotk_write_dat( iunpun, &
                                 "NUMBER_OF_BANDS", nbnd, ATTR = attr )
            !
         ELSE
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec )
            !
            CALL iotk_write_dat( iunpun, "NUMBER_OF_BANDS", nbnd )
            !
         END IF
         !
         CALL iotk_write_begin( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
         !
      END IF
      !
      k_points_loop: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            !
            CALL iotk_write_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
            !
            CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
            CALL iotk_write_dat( iunpun, &
                                 "K-POINT_COORDS", xk(:,ik), ATTR = attr )
            !
            CALL iotk_write_dat( iunpun, "WEIGHT", wk(ik) )
            !
            DO ispin = 1, nspin
               !
               cspin = iotk_index( ispin )
               !
               CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
               CALL iotk_write_dat( iunpun, "ET" // TRIM( cspin ), &
                                    et(:,ik,ispin), ATTR = attr  )
               !
               CALL iotk_write_dat( iunpun, &
                                    "OCC0" // TRIM( cspin ), occ0(:,ik,ispin) )
               CALL iotk_write_dat( iunpun, &
                                    "OCCM" // TRIM( cspin ), occm(:,ik,ispin) )
               !
            END DO
            !
            ! ... G+K vectors
            !
            filename = TRIM( wfc_filename( ".", 'gkvectors', ik ) )
            !
            CALL iotk_link( iunpun, "gkvectors", &
                            filename, CREATE = .FALSE., BINARY = .TRUE. )
            !
            filename = TRIM( wfc_filename( dirname, 'gkvectors', ik ) )
            !
         END IF
         !
         DO ispin = 1, nspin
            ! 
            ik_eff = ik + ( ispin - 1 ) * nk
            ! 
            IF ( ionode ) THEN
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( ".", 'evc', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( ".", 'evc', ik, ispin ) )
                  !
               END IF
               !
               CALL iotk_link( iunpun, "wfc", &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin ) )
                  !
               END IF
               !
            END IF
            !
            iss_wfc = ispin
            if( force_pairing ) iss_wfc = 1   ! only the WF for the first spin is allocated
            !
            IF ( PRESENT( c04 ) ) THEN
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, ispin, nspin,   &
                               c04(:,:,ik,iss_wfc), ngwt, nbnd, ig_l2g, &
                               ngw, filename, scalef )
               !
            ELSE IF ( PRESENT( c02 ) ) THEN
               !
               ib = iupdwn( iss_wfc )
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, ispin, nspin, &
                               c02(:,ib:), ngwt, nbnd, ig_l2g, ngw, &
                               filename, scalef )
               !
            END IF
            !
            IF ( ionode ) THEN
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( ".", 'evcm', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( ".", 'evcm', ik, ispin ) )
                  !
               END IF
               !
               CALL iotk_link( iunpun, "wfcm", &
                               filename, CREATE = .FALSE., BINARY = .TRUE. )
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( dirname, 'evcm', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( dirname, 'evcm', ik, ispin ) )
                  !
               END IF
               !
            END IF
            !
            IF ( PRESENT( cm4 ) ) THEN
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, ispin, nspin,   &
                              cm4(:,:,ik,iss_wfc), ngwt, nbnd, ig_l2g,  &
                              ngw, filename, scalef )
               !
            ELSE IF ( PRESENT( c02 ) ) THEN
               !
               ib = iupdwn(iss_wfc)
               !
               CALL write_wfc( iunout, ik_eff, nk*nspin, kunit, ispin, nspin, &
                               cm2(:,ib:SIZE(cm2,2)), ngwt, nbnd, ig_l2g, ngw, &
                               filename, scalef )
               !
            END IF
            !
         END DO
         !
         IF ( ionode ) &
            CALL iotk_write_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop
      !
      DO ispin = 1, nspin
         ! 
         cspin = iotk_index( ispin )
         !
         IF ( ionode .AND. PRESENT( mat_z ) ) THEN
            ! 
            filename = 'mat_z' // cspin
            !
            CALL iotk_link( iunpun, "mat_z" // TRIM( cspin ), &
                            filename, CREATE = .TRUE., BINARY = .TRUE. )
            !
            CALL iotk_write_dat( iunpun, &
                                 "mat_z" // TRIM( cspin ), mat_z(:,:,ispin) )
            !
         END IF
         !
      END DO
      !
      IF( ionode ) THEN
        !  
        CALL iotk_write_end( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
        !     
        CALL iotk_write_end( iunpun, "BAND_STRUCTURE" )
        ! 
      END IF
      !
      ! ... write matrix lambda to file
      !
      filename = TRIM( dirname ) // '/lambda.dat'
      !
      IF ( ionode ) THEN
         !
         OPEN( UNIT = 10, FILE = TRIM( filename ), &
               STATUS = 'UNKNOWN', FORM = 'UNFORMATTED' )
         !
         WRITE( 10 ) lambda0
         WRITE( 10 ) lambdam
         !
         CLOSE( UNIT = 10 )
         !
      END IF
      !
      IF ( ionode ) CALL iotk_close_write( iunpun )
      !
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( ftmp )
      DEALLOCATE( tau  )
      DEALLOCATE( ityp )
      DEALLOCATE( mill )
      !
      CALL save_history( dirname, nfi )
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
    SUBROUTINE cp_readfile( ndr, scradir, ascii, nfi, simtime, acc, nk, xk,   &
                            wk, ht, htm, htvel, gvel, xnhh0, xnhhm, vnhh,     &
                            taui, cdmi, stau0, svel0, staum, svelm, force,    &
                            vnhp, xnhp0, xnhpm, nhpcl,nhpdim,occ0, occm,      &
                            lambda0, lambdam, b1, b2, b3, xnhe0, xnhem, vnhe, &
                            ekincm, c04, cm4, c02, cm2, mat_z )
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : gamma_only, force_pairing
      USE io_files,                 ONLY : iunpun, xmlpun
      USE printout_base,            ONLY : title
      USE grid_dimensions,          ONLY : nr1, nr2, nr3
      USE smooth_grid_dimensions,   ONLY : nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
      USE gvecp,                    ONLY : ngm, ngmt, ecutp
      USE gvecs,                    ONLY : ngs, ngst
      USE gvecw,                    ONLY : ngw, ngwt, ecutw
      USE electrons_base,           ONLY : nspin, nbnd, nelt, nel, &
                                           nupdwn, iupdwn
      USE cell_base,                ONLY : ibrav, alat, celldm, symm_type, &
                                           s_to_r, r_to_s
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, pmass, &
                                           sort_tau, ityp, ions_cofmass
      USE reciprocal_vectors,       ONLY : ig_l2g, mill_l
      USE cp_main_variables,        ONLY : nprint_nfi
      USE mp,                       ONLY : mp_sum
      USE mp_global,                ONLY : intra_image_comm
      USE parameters,               ONLY : nhclm, ntypx
      USE constants,                ONLY : eps8, angstrom_au
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)    :: ndr          !  I/O unit number
      CHARACTER(LEN=*),      INTENT(IN)    :: scradir      !
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
      REAL(DP),              INTENT(INOUT) :: occ0(:,:,:)  !
      REAL(DP),              INTENT(INOUT) :: occm(:,:,:)  !
      REAL(DP),              INTENT(INOUT) :: lambda0(:,:,:) !
      REAL(DP),              INTENT(INOUT) :: lambdam(:,:,:) !
      REAL(DP),              INTENT(INOUT) :: b1(3)        !
      REAL(DP),              INTENT(INOUT) :: b2(3)        !
      REAL(DP),              INTENT(INOUT) :: b3(3)        !
      REAL(DP),              INTENT(INOUT) :: xnhe0        !
      REAL(DP),              INTENT(INOUT) :: xnhem        !
      REAL(DP),              INTENT(INOUT) :: vnhe         !  
      REAL(DP),              INTENT(INOUT) :: ekincm       !  
      COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: c04(:,:,:,:) ! 
      COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: cm4(:,:,:,:) ! 
      COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: c02(:,:)     ! 
      COMPLEX(DP), OPTIONAL, INTENT(INOUT) :: cm2(:,:)     ! 
      REAL(DP),    OPTIONAL, INTENT(INOUT) :: mat_z(:,:,:) ! 
      !
      CHARACTER(LEN=256)   :: dirname, kdirname, filename
      CHARACTER(LEN=5)     :: kindex
      CHARACTER(LEN=4)     :: cspin
      INTEGER              :: strlen
      INTEGER              :: kunit
      INTEGER              :: k1, k2, k3
      INTEGER              :: nk1, nk2, nk3
      INTEGER              :: i, j, ispin, ig, nspin_wfc, ierr, ik
      REAL(DP)             :: omega, htm1(3,3), hinv(3,3), scalef
      LOGICAL              :: found
      LOGICAL              :: tread_cm
      INTEGER, ALLOCATABLE :: mill(:,:)
      !
      ! ... variables read for testing pourposes
      !
      INTEGER               :: ibrav_
      CHARACTER(LEN=9)      :: symm_type_
      CHARACTER(LEN=3)      :: atm_(ntypx)
      INTEGER               :: nat_, nsp_, na_
      INTEGER               :: nk_, ik_, nt_
      LOGICAL               :: gamma_only_ 
      REAL(DP)              :: alat_, a1_(3), a2_(3), a3_(3)
      REAL(DP)              :: pmass_, zv_ 
      REAL(DP)              :: celldm_(6)
      INTEGER               :: ispin_, nspin_, ngwt_, nbnd_ 
      REAL(DP)              :: nelec_
      REAL(DP)              :: scalef_
      REAL(DP)              :: wk_
      INTEGER               :: nhpcl_, nhpdim_ 
      INTEGER               :: ib
      INTEGER               :: ik_eff
      REAL(DP)              :: amass_(ntypx)
      INTEGER,  ALLOCATABLE :: ityp_(:) 
      INTEGER,  ALLOCATABLE :: isrt_(:) 
      REAL(DP), ALLOCATABLE :: tau_(:,:) 
      INTEGER,  ALLOCATABLE :: if_pos_(:,:) 
      CHARACTER(LEN=256)    :: psfile_(ntypx)
      CHARACTER(LEN=80)     :: pos_unit
      REAL(DP)              :: s1, s0, cclock
      !
      ! ... look for an empty unit
      !
      CALL iotk_free_unit( iunout, ierr )
      !
      CALL errore( 'cp_readfile', &
                   'no free units to read wavefunctions', ierr )
      !
      kunit = 1
      found = .FALSE.
      !
      dirname = restart_dir( scradir, ndr )
      !
      ! ... Open XML descriptor
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/' // TRIM( xmlpun )
         !
         WRITE( stdout, '(/,3X,"reading restart file: ",A)' ) TRIM( dirname )
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                              BINARY = .FALSE., ROOT = attr, IERR = ierr )
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
         CALL iotk_scan_begin( iunpun, "STATUS", FOUND = found )
         !
         IF ( found ) THEN
            !
            CALL iotk_scan_empty( iunpun, "STEP", attr )
            CALL iotk_scan_attr( attr, "ITERATION", nfi )
            CALL iotk_scan_dat( iunpun, "TIME", simtime )
            CALL iotk_scan_dat( iunpun, "TITLE", title )
            CALL iotk_scan_end( iunpun, "STATUS" )
            !
         END IF
         !
      END IF
      !
      ! ... Read cell and positions
      !
      ALLOCATE( tau_( 3, nat ) )
      ALLOCATE( if_pos_( 3, nat ) )
      ALLOCATE( ityp_( nat ) )
      !
      IF ( ionode ) &
         CALL read_cell( ibrav_, symm_type_, celldm_, &
                         alat_, a1_, a2_, a3_, b1, b2, b3 )
      !
      IF ( ionode ) THEN
         !
         CALL read_ions( nsp_, nat_, atm_, ityp_, &
                         psfile_, amass_, tau_, if_pos_, pos_unit, ierr )
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
      ! ... read MD timesteps variables
      !
      IF ( ionode ) &
         CALL iotk_scan_begin( iunpun, "TIMESTEPS", attr, FOUND = found )
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
      scalef = 1.D0 / SQRT( ABS( omega ) )
      !
      ! ... band Structure
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin_ )
         !
         IF ( nspin_ == 2 ) THEN
            !
            CALL iotk_scan_dat( iunpun, &
                                "NUMBER_OF_ELECTRONS", nelec_, ATTR = attr )
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd_, ATTR = attr )
            !
         ELSE
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec_ )
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd_ )
            !
         END IF
         !
         IF ( ( nspin_ /= nspin ) .OR. &
              ( nbnd_ /= nbnd ) .OR. ( NINT( nelec_ ) /= nelt ) ) THEN 
            !
            attr = "electron, bands or spin do not match"
            ierr = 30
            !
            GOTO 100
            !
         END IF
         ! 
         CALL iotk_scan_begin( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
         !
      END IF
      !
      k_points_loop: DO ik = 1, nk
         !
         IF ( ionode ) THEN
            !
            CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
            !
            CALL iotk_scan_dat( iunpun, "WEIGHT", wk_ )
            !
         END IF
         !
         DO ispin = 1, nspin
            !
            cspin = iotk_index( ispin )
            !
            tread_cm = .TRUE.
            !
            ik_eff = ik + ( ispin - 1 ) * nk
            !
            IF ( ionode ) THEN
               !
               CALL iotk_scan_dat( iunpun, &
                                   "OCC0" // TRIM( cspin ), occ0(:,ik,ispin),  &
                                   FOUND = found, IERR = ierr )
               !
               IF ( .NOT. found ) THEN
                  !
                  CALL iotk_scan_dat( iunpun, &
                                      "OCC" // TRIM( cspin ), occ0(:,ik,ispin) )
                  !
               END IF
               !
               CALL iotk_scan_dat( iunpun, &
                                   "OCCM" // TRIM( cspin ), occm(:,ik,ispin ), &
                                   FOUND = found, IERR = ierr )
               !
               occ0(:,ik,ispin) = occ0(:,ik,ispin) * wk_
               !
               IF ( .NOT. found ) THEN
                  !
                  occm(:,ik,ispin) = occ0(:,ik,ispin)
                  !
                  tread_cm = .FALSE.
                  !
               END IF
               !
            END IF
            !
            CALL mp_bcast( tread_cm, ionode_id, intra_image_comm )
            !
            IF ( ionode ) THEN
               !
               IF ( nspin == 1 ) THEN
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc', ik ) )
                  !
               ELSE
                  !
                  filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin ) )
                  !
               END IF
               !
            END IF
            !
            IF( .NOT. ( ispin > 1 .AND. force_pairing ) ) THEN
               !
               ! Only WF with spin 1 are needed when force_pairing is active
               !
               IF ( PRESENT( c04 ) ) THEN
                  !
                  CALL read_wfc( iunout, ik_eff , nk, kunit, ispin_, nspin_, &
                                 c04(:,:,ik,ispin), ngwt_, nbnd_, ig_l2g, &
                                 ngw, filename, scalef_ )
                  !
               ELSE IF ( PRESENT( c02 ) ) THEN
                  !
                  ib = iupdwn(ispin)
                  !
                  CALL read_wfc( iunout, ik_eff , nk, kunit, ispin_, nspin_, &
                                 c02(:,ib:), ngwt_, nbnd_, ig_l2g, ngw, &
                                 filename, scalef_ )
                  !
               END IF
               !
            END IF
            !
            IF ( tread_cm ) THEN 
               !
               IF ( ionode ) THEN
                  !
                  IF ( nspin == 1 ) THEN
                     !
                     filename = TRIM( wfc_filename( dirname, 'evcm', ik ) )
                     !
                  ELSE
                     !
                     filename = TRIM( wfc_filename( dirname, &
                                                    'evcm', ik, ispin ) )
                     !
                  END IF
                  !
               END IF
               !
               IF( .NOT. ( ispin > 1 .AND. force_pairing ) ) THEN
                  !
                  ! Only WF with spin 1 are needed when force_pairing is active
                  !
                  IF ( PRESENT( cm4 ) ) THEN
                     !
                     CALL read_wfc( iunout, ik_eff, nk, kunit, ispin_, nspin_,   &
                                    cm4(:,:,ik,ispin), ngwt_, nbnd_, ig_l2g, &
                                    ngw, filename, scalef_ )
                     !
                  ELSE IF( PRESENT( cm2 ) ) THEN
                     !
                     ib = iupdwn(ispin)
                     !
                     CALL read_wfc( iunout, ik_eff, nk, kunit, ispin_, nspin_, &
                                    cm2(:,ib:), ngwt_, nbnd_, ig_l2g, ngw, &
                                    filename, scalef_ )
                     !
                  END IF
                  !
               END IF
               !
            ELSE
               !
               IF ( PRESENT( cm4 ) ) THEN
                  !
                  cm4 = c04
                  !
               ELSE IF( PRESENT( cm2 ) ) THEN
                  !
                  cm2 = c02
                  !
               END IF
               !
            END IF
            !
         END DO
         !
         IF ( ionode ) &
            CALL iotk_scan_end( iunpun, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO k_points_loop
      !
      DO ispin = 1, nspin
         !
         IF ( ionode .AND. PRESENT( mat_z ) ) &
            CALL iotk_scan_dat( iunpun, "mat_z" // &
                              & TRIM( iotk_index( ispin ) ), mat_z(:,:,ispin) )
         !
      END DO
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( iunpun, "EIGENVALUES_AND_EIGENVECTORS" )
         !
         CALL iotk_scan_end( iunpun, "BAND_STRUCTURE" )
         !
      END IF
      !
 100  CONTINUE
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      CALL mp_bcast( attr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_readfile ', TRIM( attr ), ierr )
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

      CALL mp_bcast( occ0(:,:,:), ionode_id, intra_image_comm )
      CALL mp_bcast( occm(:,:,:), ionode_id, intra_image_comm )
      !
      IF ( PRESENT( mat_z ) ) &
         CALL mp_bcast( mat_z(:,:,:), ionode_id, intra_image_comm )
      !
      IF ( ionode ) &
         CALL iotk_close_read( iunpun )
      !
      ! ... read matrix lambda to file
      !
      filename = TRIM( dirname ) // '/lambda.dat'
      !
      IF ( ionode ) THEN
         !
         INQUIRE( file = TRIM( filename ), EXIST = found )
         !
         IF ( found ) THEN
            !
            OPEN( UNIT = 10, FILE = TRIM( filename ), &
                  STATUS = 'OLD', FORM = 'UNFORMATTED' )
            !
            READ( 10 ) lambda0
            READ( 10 ) lambdam
            !
            CLOSE( UNIT = 10 )
            !
         ELSE
            !
         END IF
         !
      END IF
      !
      DO ispin = 1, nspin
         DO ib = 1, SIZE( lambda0, 2 )
            CALL mp_bcast( lambda0( :, ib, ispin), ionode_id, intra_image_comm )
            CALL mp_bcast( lambdam( :, ib, ispin), ionode_id, intra_image_comm )
         END DO
      END DO
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
      if (ionode.and.(nprint_nfi.eq.-2)) then
         write( stdout,*) 'nprint_nfi= ',nprint_nfi
         CALL read_print_counter( nprint_nfi, scradir, ndr )
         write( stdout,*) 'nprint_nfi= ',nprint_nfi
      endif
      CALL mp_bcast( nprint_nfi, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE cp_readfile
    ! 
    !------------------------------------------------------------------------
    SUBROUTINE cp_read_wfc( ndr, scradir, ik, nk, ispin, nspin, c2, c4, tag )
      !------------------------------------------------------------------------
      !
      USE electrons_base,     ONLY : nbnd, iupdwn
      USE reciprocal_vectors, ONLY : ngwt, ngw, ig_l2g
      !
      IMPLICIT NONE
      !
      INTEGER,               INTENT(IN)  :: ndr
      CHARACTER(LEN=*),      INTENT(IN)  :: scradir
      INTEGER,               INTENT(IN)  :: ik, ispin, nk, nspin
      CHARACTER,             INTENT(IN)  :: tag
      COMPLEX(DP), OPTIONAL, INTENT(OUT) :: c2(:,:)
      COMPLEX(DP), OPTIONAL, INTENT(OUT) :: c4(:,:,:,:)
      !
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: ik_eff, ib, kunit, ispin_, nspin_, ngwt_, nbnd_
      REAL(DP)           :: scalef
      !
      kunit = 1
      !
      ik_eff = ik + ( ispin - 1 ) * nk
      !
      dirname = restart_dir( scradir, ndr )
      !
      IF ( tag /= 'm' ) THEN
         !
         IF ( nspin == 1 ) THEN
            !
            filename = TRIM( wfc_filename( dirname, 'evc', ik ) )
            !
         ELSE
            !
            filename = TRIM( wfc_filename( dirname, 'evc', ik, ispin ) )
            !
         END IF
         !
      ELSE
         !
         IF ( nspin == 1 ) THEN
            !
            filename = TRIM( wfc_filename( dirname, 'evcm', ik ) )
            !
         ELSE
            !
            filename = TRIM( wfc_filename( dirname, 'evcm', ik, ispin ) )
            !
         END IF
         !
      END IF
      !
      IF ( PRESENT( c4 ) ) THEN
         !
         CALL read_wfc( iunout, ik_eff, nk, kunit, ispin_, nspin_,  &
                        c4(:,:,ik,ispin), ngwt_, nbnd_, ig_l2g, &
                        ngw, filename, scalef )
         !
      ELSE IF ( PRESENT( c2 ) ) THEN
         !
         ib = iupdwn(ispin)
         !
         CALL read_wfc( iunout, ik_eff, nk, kunit, ispin_, nspin_, &
                        c2(:,ib:), ngwt_, nbnd_, ig_l2g, ngw,  &
                        filename, scalef )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE cp_read_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE cp_read_cell( ndr, scradir, ascii, ht, &
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
      CHARACTER(LEN=*), INTENT(IN)    :: scradir
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
      CHARACTER(LEN=9) :: symm_type_
      !
      !
      dirname = restart_dir( scradir, ndr ) 
      !
      filename = TRIM( dirname ) // '/' // TRIM( xmlpun )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                              BINARY = .FALSE., ROOT = attr, IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'cp_read_cell', &
                   'cannot open restart file for reading', ierr )
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
            CALL read_cell( ibrav_, symm_type_, celldm_, &
                            alat_, a1_, a2_, a3_, b1_, b2_, b3_ )
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
      !
      RETURN
      !
    END SUBROUTINE cp_read_cell
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_cell( ibrav, symm_type, &
                          celldm, alat, a1, a2, a3, b1, b2, b3 )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(OUT) :: ibrav
      CHARACTER(LEN=*), INTENT(OUT) :: symm_type
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
      IF ( ibrav == 0 ) &
         CALL iotk_scan_dat( iunpun, "CELL_SYMMETRY", symm_type )
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
    !------------------------------------------------------------------------
    SUBROUTINE read_ions( nsp, nat, atm, ityp, psfile, &
                          amass, tau, if_pos, pos_unit, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,            INTENT(OUT) :: nsp, nat
      CHARACTER(LEN=3),   INTENT(OUT) :: atm(:)
      INTEGER,            INTENT(OUT) :: ityp(:)
      CHARACTER(LEN=256), INTENT(OUT) :: psfile(:)
      REAL(DP),           INTENT(OUT) :: amass(:)
      REAL(DP),           INTENT(OUT) :: tau(:,:)
      INTEGER,            INTENT(OUT) :: if_pos(:,:)
      INTEGER,            INTENT(OUT) :: ierr
      CHARACTER(LEN=*),   INTENT(OUT) :: pos_unit
      !
      LOGICAL          :: found
      INTEGER          :: i
      CHARACTER(LEN=3) :: lab
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunpun, "IONS", FOUND = found )
      !
      IF ( .NOT. found ) THEN
         !
         ierr = 1
         !
         RETURN
         !
      END IF
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMS",   nat )
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
      !
      IF ( nsp > SIZE( atm ) .OR. nat > SIZE( ityp ) ) THEN
         !
         ierr = 10
         !
         CALL iotk_scan_end( iunpun, "IONS" )
         !
         RETURN
         !
      END IF
      !
      DO i = 1, nsp
         !
         CALL iotk_scan_dat( iunpun, "ATOM_TYPE", atm(i) )
         !
         CALL iotk_scan_dat( iunpun, &
                             TRIM( atm(i) )//"_MASS", amass(i), ATTR = attr )
         !
         CALL iotk_scan_dat( iunpun, &
                             "PSEUDO_FOR_" // TRIM( atm(i) ), psfile(i) )
         !
      END DO
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_ATOMIC_POSITIONS", attr )
      CALL iotk_scan_attr( attr, "UNITS", pos_unit  )
      !
      DO i = 1, nat
         !
         CALL iotk_scan_empty( iunpun, "ATOM" // TRIM( iotk_index( i ) ), attr )
         CALL iotk_scan_attr( attr, "SPECIES", lab )
         CALL iotk_scan_attr( attr, "INDEX",   ityp(i) )
         CALL iotk_scan_attr( attr, "tau",     tau(:,i) )
         CALL iotk_scan_attr( attr, "if_pos",  if_pos(:,i) )
         !
      END DO
      !
      CALL iotk_scan_end( iunpun, "IONS" )
      !
      RETURN
      !
    END SUBROUTINE read_ions
    !
END MODULE cp_restart
