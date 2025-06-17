!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE cp_restart_new
  !-----------------------------------------------------------------------------
  !! This module contains subroutines to write and read data required to
  !! restart a calculation from the disk.
  !! Important notice:
  !
  !! * only one processor writes (the one for which ionode=.TRUE.);
  !! * all processors read the xml file;
  !! * one processor per band group reads the wavefunctions distributes 
  !!   them within their band group;
  !! * lambda matrices are read by one processor, broadcast to all others.
  !
  USE kinds,     ONLY : DP
  !
  USE qes_types_module 
  USE qes_libs_module  
  USE io_files,  ONLY : iunpun, xmlfile, restart_dir, &
       qexsd_fmt, qexsd_version, create_directory
  USE io_base,   ONLY : write_wfc, read_wfc, write_rhog
  !
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE matrix_inversion
#if defined (__fox)
    USE FoX_wxml
    USE FoX_dom,     ONLY : Node, parseFile, item, getElementsByTagname, &
                            hasAttribute, extractDataAttribute, &
                            extractDataContent, destroy
#else
    USE     wxml
    USE     dom,     ONLY : Node, parseFile, item, getElementsByTagname, &
                            hasAttribute, extractDataAttribute, &
                            extractDataContent, destroy
#endif
  !
  IMPLICIT NONE
  !
  SAVE
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
                             iupdwn_tot, nupdwn_tot, wfc, dt )
      !------------------------------------------------------------------------
      !
      USE control_flags,            ONLY : gamma_only, force_pairing, trhow, &
                                           tksw, do_makov_payne, smallmem,   &
                                           mbd_vdw,                           &
                                           llondon, lxdm, ts_vdw, tfor, tpre
      USE control_flags,            ONLY : lwfpbe0nscf, lwfnscf, lwf ! Lingzhu Kong
      USE constants,                ONLY : e2
      USE parameters,               ONLY : ntypx
      USE dener,                    ONLY : detot
      USE io_files,                 ONLY : psfile, pseudo_dir
      USE mp_images,                ONLY : intra_image_comm, me_image, &
                                           nproc_image
      USE mp_bands,                 ONLY : my_bgrp_id, intra_bgrp_comm, &
                                           root_bgrp, root_bgrp_id
      USE run_info,                 ONLY : title
      USE gvect,                    ONLY : ngm, ngm_g, ecutrho
      USE gvecs,                    ONLY : ngms_g, ecuts
      USE gvecw,                    ONLY : ngw, ngw_g, ecutwfc
      USE gvect,                    ONLY : ig_l2g, mill
      USE electrons_base,           ONLY : nspin, nelt, nel, nudx
      USE cell_base,                ONLY : ibrav, alat, tpiba, s_to_r
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           amass, iforce, ityp 
      USE funct,                    ONLY : get_dft_name, dft_is_nonlocc, get_nonlocc_name
      USE xc_lib,                   ONLY : xclib_dft_is, xclib_get_exx_fraction, &
                                           get_screening_parameter
      USE ldaU_cp,                  ONLY : lda_plus_U, ns, Hubbard_l, Hubbard_n, &
                                           Hubbard_lmax, Hubbard_U
      USE energies,                 ONLY : enthal, ekin, eht, esr, eself, &
                                           epseu, enl, exc, vave
      USE mp,                       ONLY : mp_sum, mp_barrier
      USE fft_base,                 ONLY : dfftp, dffts, dfftb
      USE fft_rho,                  ONLY : rho_r2g
      USE upf_ions,                 ONLY : n_atom_wfc
      USE uspp_param,               ONLY : upf
      USE london_module,            ONLY : scal6, lon_rcut, in_c6
      USE tsvdw_module,             ONLY : vdw_isolated, vdw_econv_thr
      USE clib_wrappers,            ONLY : f_copy
      USE uspp,                     ONLY : okvan
      USE input_parameters,         ONLY : vdw_corr, starting_ns_eigenvalue
      USE qexsd_init, ONLY: qexsd_init_convergence_info, qexsd_init_algorithmic_info,  & 
                          qexsd_init_atomic_species, qexsd_init_atomic_structure,      &
                          qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
                          qexsd_init_magnetization,qexsd_init_band_structure, &
                          qexsd_init_dipole_info, qexsd_init_total_energy,    &
                          qexsd_init_forces,qexsd_init_stress,                &
                          qexsd_init_outputElectricField, qexsd_init_vdw,     &
                          qexsd_init_hybrid, qexsd_init_dftU 
      USE qexsd_input, ONLY: qexsd_init_k_points_ibz
      USE qexsd_module, ONLY: qexsd_openschema, qexsd_closeschema, qexsd_xf,  &
                              qexsd_add_all_clocks
      !
      INTEGER,               INTENT(IN) :: ndw          !
      LOGICAL,               INTENT(IN) :: ascii        !
      INTEGER,               INTENT(IN) :: nfi          ! index of the current step
      REAL(DP),              INTENT(IN) :: simtime      ! simulated time
      REAL(DP),              INTENT(IN) :: dt           ! integration timestep
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
      INTEGER,TARGET,        INTENT(IN) :: nupdwn(:)    ! 
      INTEGER,               INTENT(IN) :: iupdwn_tot(:)! 
      INTEGER,               INTENT(IN) :: nupdwn_tot(:)! 
      REAL(DP),              INTENT(IN) :: wfc(:,:)     ! BS 
      !
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname
      CHARACTER(LEN=320)    :: filename, sourcefile
      INTEGER               :: ik_eff
      INTEGER               :: k1, k2, k3
      INTEGER               :: nk1, nk2, nk3
      INTEGER               :: j, i, iss, ig, nspin_wfc, iss_wfc
      INTEGER               :: is, ia, isa, ik, ierr
      REAL(DP), ALLOCATABLE :: ftmp(:,:)
      REAL(DP), ALLOCATABLE :: tau(:,:)
      COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
      REAL(DP)              :: omega, htm1(3,3), h(3,3)
      REAL(DP)              :: stress(3,3)
      REAL(DP)              :: a1(3), a2(3), a3(3)
      REAL(DP)              :: b1(3), b2(3), b3(3)
      REAL(DP)              :: wk_(2), nelec
      REAL(DP)              :: scalef
      LOGICAL               :: lsda
      REAL(DP)              :: s0, s1
      INTEGER               :: natomwfc, nbnd_, nb, ib
      REAL(DP), ALLOCATABLE :: mrepl(:,:)
      LOGICAL               :: exst
      TYPE(output_type) :: output_obj
      LOGICAL :: is_hubbard(ntypx), empirical_vdw  
      TYPE(occupations_type)       :: bands_occu 
      TYPE(k_points_IBZ_type)      :: k_points_IBZ 
      CHARACTER(LEN=6), EXTERNAL   :: int_to_char
      TYPE (vdW_type),POINTER      :: vdW_ 
      TYPE (dftU_type),POINTER     :: dftU_ 
      TYPE (hybrid_type),POINTER   :: hybrid_
      REAL(DP),ALLOCATABLE         :: london_c6_(:) 
      CHARACTER(LEN=3),ALLOCATABLE :: species_(:) 
      REAL(DP),TARGET              :: lond_rcut_, lond_s6_, ts_vdw_econv_thr_      
      REAL(DP),POINTER             :: london_s6_pt, lonrcut_opt, ts_thr_opt 
      INTEGER,POINTER              :: nbnd_pt, nbnd_up_pt, nbnd_dw_pt 
      CHARACTER(LEN=20),TARGET     :: non_locc_, vdw_corr_ 
      CHARACTER(LEN=20),POINTER    :: non_locc_opt, vdw_corr_opt 
      LOGICAL,POINTER              :: ts_isol_opt 
      LOGICAL,TARGET               :: ts_vdW_isolated_ 
      !
      ! ... subroutine body
      !
      NULLIFY( london_s6_pt, lonrcut_opt, ts_thr_opt, nbnd_pt, nbnd_up_pt, nbnd_dw_pt) 
      NULLIFY ( vdW_, dftU_, hybrid_, non_locc_opt, vdw_corr_opt, ts_isol_opt )
      CALL start_clock('restart')
      !
      IF( force_pairing ) &
            CALL errore('cp_writefile',' force pairing not implemented', 1 )
      IF( tksw ) &
            CALL infomsg('cp_writefile',' Kohn-Sham states not written' )
      !
      lsda = ( nspin == 2 )
      IF( lsda ) THEN
         !
         !  check if the array storing wave functions is large enought
         !
         IF( SIZE( c02, 2 ) < ( iupdwn( 2 ) + nupdwn(1) - 1 ) ) &
            CALL errore('cp_writefile',' wrong dimensions for wave functions', 1 )
         !
      END IF
      !
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
      ! ... Compute array tau
      !
      ALLOCATE( tau( 3, nat ) )
      !
      natomwfc =  n_atom_wfc ( nat, ityp ) 
      !
      CALL s_to_r( stau0, tau, nat, h )
      !
      nbnd_    = nupdwn(1) 
      ALLOCATE( ftmp( nbnd_ , nspin ) )
      ftmp = 0.0d0
      DO iss = 1, nspin
         ftmp( 1:nupdwn(iss), iss ) = occ0( iupdwn(iss) : iupdwn(iss) + nupdwn(iss) - 1 )
      END DO
      !
      ! XML descriptor
      ! 
      dirname=restart_dir(ndw)
      WRITE( stdout, '(/,3X,"writing restart file (with schema): ",A)' ) &
             TRIM(dirname)
      !
      CALL create_directory( TRIM(dirname) )
      !
      IF ( ionode ) THEN
         !
         ! ... here we init the variables and finally write them to file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         output_obj%tagname="output"
         output_obj%lwrite = .TRUE.
!-------------------------------------------------------------------------------
! ... CONVERGENCE_INFO - TO BE VERIFIED   
!-------------------------------------------------------------------------------
!! @note set lwrite to false for this  element P. Delugas 
         CALL qexsd_init_convergence_info(output_obj%convergence_info, &
              scf_has_converged = .FALSE., n_scf_steps=0, scf_error=0.0_dp, &
              n_opt_steps=0, grad_norm=0.0_dp )
         output_obj%convergence_info%lwrite = .FALSE. 
         !
!-------------------------------------------------------------------------------
! ... ALGORITHMIC_INFO
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_algorithmic_info(output_obj%algorithmic_info, &
              real_space_beta=.FALSE., real_space_q=.FALSE., uspp=okvan, paw=.FALSE.)
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_SPECIES
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_atomic_species(output_obj%atomic_species, nsp, atm,&
                 psfile, amass)
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_STRUCTURE
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_atomic_structure(output_obj%atomic_structure, nsp, atm, ityp, &
              nat, tau(:,:), alat, alat*a1(:), alat*a2(:), alat*a3(:), ibrav)
         !
!-------------------------------------------------------------------------------
! ... BASIS SET
!-------------------------------------------------------------------------------
         CALL qexsd_init_basis_set(output_obj%basis_set,gamma_only, ecutwfc/e2, ecutrho/e2, &
              dfftp%nr1, dfftp%nr2, dfftp%nr3, dffts%nr1, dffts%nr2, dffts%nr3, &
              .FALSE., dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, ngms_g, ngw_g, &
              b1(:), b2(:), b3(:) )
!-------------------------------------------------------------------------------
! ... XC FUNCTIONAL
!-------------------------------------------------------------------------------
        dft_name = get_dft_name()
        IF ( lda_plus_U) THEN
           ALLOCATE (dftU_) 
           is_hubbard(:) = (Hubbard_U(:) > 0.0_dp)
           CALL qexsd_init_dftU(OBJ = dftU_, NSP = nsp, PSD = upf(1:nsp)%psd, SPECIES = atm(1:nsp), &
                                ITYP = ityp, IS_HUBBARD  = is_hubbard, LDA_PLUS_U_KIND = 0,         &
                                U_PROJECTION_TYPE = 'atomic', U = Hubbard_U, STARTING_NS = starting_ns_eigenvalue) 
        END IF
        !
        IF (xclib_dft_is('hybrid'))  THEN 
           ALLOCATE (hybrid_) 
           CALL qexsd_init_hybrid(OBJ = hybrid_, DFT_IS_HYBRID = .TRUE. , ECUTFOCK = ecutwfc, &
                                 EXX_FRACTION = xclib_get_exx_fraction(), SCREENING_PARAMETER = get_screening_parameter(),&
                                 EXXDIV_TREATMENT = 'none',  X_GAMMA_EXTRAPOLATION = .FALSE.) 
        END IF 
        empirical_vdW = ( TRIM(vdw_corr) /= 'none' )  
        IF ( empirical_vdW .OR. dft_is_nonlocc() ) THEN 
           ALLOCATE (vdw_)
           IF (empirical_vdw) THEN 
              vdw_corr_ = TRIM (vdw_corr) 
              vdw_corr_opt => vdw_corr_ 
              SELECT CASE(TRIM (vdw_corr_)) 
                CASE ( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d') 
                    lond_s6_ = scal6
                    london_s6_pt => lond_s6_
                    lond_rcut_ = lon_rcut
                    lonrcut_opt => lond_rcut_
                    IF (ANY( in_c6(1:nsp) .NE. -1._DP )) THEN
                       ALLOCATE (london_c6_(nsp), species_(nsp))
                       london_c6_(1:nsp) = in_c6(1:nsp)
                       species_(1:nsp)  = atm(1:nsp)
                   END IF
                CASE ( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler') 
                    ts_vdw_isolated_ = vdw_isolated
                    ts_isol_opt => ts_vdw_isolated_
                    ts_vdw_econv_thr_ = vdw_econv_thr
                    ts_thr_opt => ts_vdw_econv_thr_
              END SELECT 
           ELSE
              vdw_corr_ = 'none'
              vdw_corr_opt => vdw_corr_ 
           END IF 
           IF ( dft_is_nonlocc() ) THEN 
              non_locc_ = TRIM ( get_nonlocc_name()) 
              non_locc_opt => non_locc_ 
           END IF 
           CALL qexsd_init_vdw(vdW_, NON_LOCAL_TERM = non_locc_opt, VDW_CORR = vdw_corr_opt, &
                                  TS_THR = ts_thr_opt, TS_ISOL = ts_isol_opt,& 
                                  LONDON_S6 = london_s6_pt, LONDON_C6 = london_c6_, LONDON_RCUT = lonrcut_opt,& 
                                  SPECIES = species_ )
        END IF    

        CALL qexsd_init_dft(output_obj%dft, dft_name, hybrid_, vdW_, dftU_)
        IF (ASSOCIATED(dftU_)) THEN
           CALL qes_reset(dftU_) 
           DEALLOCATE(dftU_) 
        END IF 
        IF (ASSOCIATED(vdW_) ) THEN 
           CALL qes_reset(vdW_) 
           DEALLOCATE(vdW_) 
        END IF 
        IF ( ASSOCIATED(hybrid_)) THEN 
           CALL qes_reset(hybrid_) 
           DEALLOCATE(hybrid_) 
        END IF 

!-------------------------------------------------------------------------------
! ... MAGNETIZATION
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_magnetization(output_obj%magnetization, LSDA = lsda, NONCOLIN = .false.,&
              SPINORBIT = .false., ABSOLUTE_MAG = 0.d0, ATM = atm, ITYP = ityp )
         output_obj%magnetization_ispresent = lsda
         !
!-------------------------------------------------------------------------------
! ... TOTAL ENERGY
!-------------------------------------------------------------------------------
         CALL  qexsd_init_total_energy(output_obj%total_energy, ETOT = enthal , &
                              EHART = eht, VTXC = vave, ETXC = exc )
!-------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------
         ! TEMP
         IF (lsda) THEN
            wk_ = 1.0_dp
         ELSE
            wk_ = 2.0_dp
         END IF
         CALL qexsd_init_k_points_ibz( k_points_ibz, 'Gamma', &
              'CP',1,1,1,0,0,0,1,alat,a1,.false.,xk,wk_) 
         bands_occu%tagname="occupations_kind"
         bands_occu%lread=.false.
         bands_occu%lwrite=.true.
         bands_occu%spin_ispresent=lsda
         bands_occu%occupations="fixed"
         ! TEMP
         IF (lsda) THEN
            nbnd_up_pt => nupdwn(1) 
            nbnd_dw_pt => nupdwn(2)
         ELSE 
            nbnd_pt => nupdwn(1) 
         END IF 
         CALL qexsd_init_band_structure( OBJ = output_obj%band_structure, LSDA = lsda, NONCOLIN = .FALSE., &
                    LSPINORB=.FALSE., NELEC = nelec, N_WFC_AT = natomwfc,  ET=et, WG = ftmp , NKS = nspin ,&
                    XK = xk , NGK=[ngw_g], WK=wk_, STARTING_KPOINTS= k_points_IBZ, OCCUPATIONS_KIND= bands_occu,& 
                    WF_COLLECTED = .TRUE., NBND = nbnd_pt, NBND_UP = nbnd_up_pt, NBND_DW = nbnd_dw_pt )
         CALL qes_reset (bands_occu)
         CALL qes_reset (k_points_IBZ)
!-------------------------------------------------------------------------------
! ... FORCES
!-------------------------------------------------------------------------------
         !
         output_obj%forces_ispresent=tfor
         CALL qexsd_init_forces(output_obj%forces,nat,force,tfor)
         !
!-------------------------------------------------------------------------------
! ... STRESS - TO BE VERIFIED
!-------------------------------------------------------------------------------
         output_obj%stress_ispresent=tpre
         ! FIXME: may be wrong or incomplete
         IF ( tpre) stress = -MATMUL( detot, ht ) / omega
         CALL qexsd_init_stress(output_obj%stress, stress, tpre )
!------------------------------------------------------------------------------
!--- TIMING 
!------------------------------------------------------------------------------
        CALL qexsd_add_all_clocks() 

!-------------------------------------------------------------------------------
! ... non existent or not implemented fields
!-------------------------------------------------------------------------------
         output_obj%symmetries_ispresent=.false.
         output_obj%electric_field_ispresent=.false.
!-------------------------------------------------------------------------------
! ... ACTUAL WRITING
!-------------------------------------------------------------------------------
         !
         filename = xmlfile(ndw)
         CALL qexsd_openschema( filename, iunpun, 'CPV', title)
         CALL qes_write (qexsd_xf, output_obj)
         CALL qes_reset (output_obj)
         !
         ! CP-SPECIFIC CELL variables
         !
         CALL cp_writecp( qexsd_xf, nfi, simtime, ekin, eht, esr, eself, &
              epseu, enl, exc, vave, enthal, acc, stau0, svel0, taui, cdmi,&
              force, nhpcl, nhpdim, xnhp0, vnhp, ekincm, xnhe0, vnhe, ht,&
              htvel, gvel, xnhh0, vnhh, staum, svelm, xnhpm, xnhem, htm, xnhhm, dt)
         !
         ! Wannier function centers
         !
         IF ( lwf ) CALL cp_writecenters ( qexsd_xf, h, wfc)
         !
         CALL qexsd_closeschema()
         !
!-------------------------------------------------------------------------------
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... WRITE WAVEFUNCTIONS AND LAMBDA MATRICES
!-------------------------------------------------------------------------------
      DO iss = 1, nspin
         !
         ik_eff = iss
         ib = iupdwn(iss)
         nb = nupdwn(iss)
         !
         IF ( my_bgrp_id == root_bgrp_id ) THEN
            !
            ! wfc collected and written by the root processor of the first
            ! band group of each pool/image - no warranty it works for nbgrp >1
            !
            filename = TRIM(dirname) // 'wfc' // TRIM(int_to_char(ik_eff))
            CALL write_wfc( iunpun, filename, root_bgrp, intra_bgrp_comm, &
                 ik_eff, xk(:,1), iss, nspin, c02(:,ib:ib+nb-1), ngw_g, &
                 gamma_only, nb, ig_l2g, ngw,  &
                 tpiba*b1, tpiba*b2, tpiba*b3, mill, scalef )
            ! wavefunctions at time t-dt
            filename = TRIM(dirname) // 'wfcm' // TRIM(int_to_char(ik_eff))
            CALL write_wfc( iunpun, filename, root_bgrp, intra_bgrp_comm, &
                 ik_eff, xk(:,1), iss, nspin, cm2(:,ib:ib+nb-1), ngw_g, &
                 gamma_only, nb, ig_l2g, ngw,  &
                 tpiba*b1, tpiba*b2, tpiba*b3, mill, scalef )
         END IF
         !
         ! matrix of orthogonality constrains lambda at time t
         filename = TRIM(dirname) // 'lambda' // TRIM(int_to_char(ik_eff))
         CALL cp_write_lambda( filename, iunpun, iss, nspin, nudx, &
              lambda0(:,:,iss), ierr )
         ! matrix of orthogonality constrains lambda at time t-dt
         filename = TRIM(dirname) // 'lambdam' // TRIM(int_to_char(ik_eff))
         CALL cp_write_lambda( filename, iunpun, iss, nspin, nudx, &
              lambdam(:,:,iss), ierr )
         !
      END DO
!-------------------------------------------------------------------------------
! ... WRITE PSEUDOPOTENTIALS
!-------------------------------------------------------------------------------
     !
     ! ... copy pseudopotential files into the .save directory
     !
     IF ( ionode ) THEN
        DO is = 1, nsp
           sourcefile= TRIM(pseudo_dir)//psfile(is)
           filename  = TRIM(dirname)//psfile(is)
           IF ( TRIM(sourcefile) /= TRIM(filename) ) &
                ierr = f_copy(sourcefile, filename)
        END DO
     END IF
     !
!-------------------------------------------------------------------------------
! ... CHARGE DENSITY
!-------------------------------------------------------------------------------
      !
     IF (trhow) THEN
        ! Workaround: input rho in real space, bring it to reciprocal space
        ! To be reconsidered once the old I/O is gone
        ALLOCATE ( rhog(ngm, nspin) )
        CALL rho_r2g (dfftp,rho, rhog)
        filename = TRIM(dirname) // 'charge-density' 
        ! Only the first band group collects and writes
        
        IF ( my_bgrp_id == root_bgrp_id ) THEN
           !
           !^^ ... TEMPORARY FIX (newlsda) ...
           IF ( lsda ) THEN
              rhog(:,1) = rhog(:,1) + rhog(:,2) 
              rhog(:,2) = rhog(:,1) - rhog(:,2)*2._dp
           ENDIF
           !^^.......................
           !      
           CALL write_rhog &
                ( filename, root_bgrp, intra_bgrp_comm, &
                tpiba*b1, tpiba*b2, tpiba*b3, gamma_only, &
                mill, ig_l2g, rhog )
        ENDIF
        !
        DEALLOCATE ( rhog )
     END IF
     !
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( ftmp )
      DEALLOCATE( tau  )
      !
      CALL stop_clock('restart')
      CALL print_clock('restart')
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
                            ekincm, c02, cm2, wfc, dt )
      !------------------------------------------------------------------------
      !! Read XML file.
      !
      USE control_flags,            ONLY : gamma_only, force_pairing, llondon,&
                                           ts_vdw, mbd_vdw, lxdm, iverbosity, lwf
      USE run_info,                 ONLY : title
      USE gvect,                    ONLY : ngm
      USE gvecw,                    ONLY : ngw, ngw_g
      USE electrons_base,           ONLY : nspin, nbnd, nupdwn, iupdwn, nudx
      USE cell_base,                ONLY : ibrav, alat, s_to_r, r_to_s
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           ityp, ions_cofmass
      USE gvect,       ONLY : ig_l2g, mill
      USE cp_main_variables,        ONLY : nprint_nfi
      USE ldaU_cp,                  ONLY : lda_plus_U, ns, Hubbard_l, Hubbard_n, &
                                           Hubbard_lmax, Hubbard_U
      USE mp,                       ONLY : mp_sum, mp_bcast
      USE mp_global,                ONLY : nproc_file, nproc_pool_file, &
                                           nproc_image_file, ntask_groups_file,&
                                           nproc_bgrp_file, nproc_ortho_file
      USE mp_pools,                 ONLY : root_pool, intra_pool_comm
      USE parameters,               ONLY : ntypx
      USE constants,                ONLY : eps8, angstrom_au, pi
      USE qes_types_module,         ONLY : output_type, parallel_info_type, &
                                           general_info_type
      USE qes_read_module,          ONLY : qes_read
      USE london_module,            ONLY : scal6, lon_rcut, in_c6
      USE tsvdw_module,             ONLY : vdw_isolated, vdw_econv_thr
      USE qexsd_copy, ONLY:  qexsd_copy_geninfo, qexsd_copy_parallel_info, &
           qexsd_copy_atomic_species, qexsd_copy_atomic_structure, &
           qexsd_copy_basis_set, qexsd_copy_dft, qexsd_copy_band_structure
      !
      INTEGER,               INTENT(IN)    :: ndr          !  I/O unit number
      LOGICAL,               INTENT(IN)    :: ascii        !
      INTEGER,               INTENT(INOUT) :: nfi          ! index of the current step
      REAL(DP),              INTENT(INOUT) :: simtime      ! simulated time
      REAL(DP),              INTENT(INOUT) :: dt           ! integration timestep
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
      !
      CHARACTER(LEN=256)   :: dirname, filename
      INTEGER              :: strlen
      INTEGER              :: k1, k2, k3
      INTEGER              :: nk1, nk2, nk3
      INTEGER              :: i, j, iss, ig, nspin_wfc, ierr, ik
      REAL(DP)             :: omega, htm1(3,3), hinv(3,3), scalef
      LOGICAL              :: found
      !
      ! ... variables read for testing purposes
      !
      INTEGER               :: ibrav_
      CHARACTER(LEN=3)      :: atm_(ntypx)
      INTEGER               :: nat_, nsp_, na_
      INTEGER               :: nk_, isk_(2), nt_, natomwfc
      LOGICAL               :: gamma_only_ , lsda_
      REAL(DP)              :: alat_, a1_(3), a2_(3), a3_(3)
      REAL(DP)              :: zv_ 
      REAL(DP)              :: ecutwfc_, ecutrho_
      INTEGER               :: nr1,nr2,nr3,nr1s,nr2s,nr3s,nr1b,nr2b,nr3b
      INTEGER               :: ngm_g, ngms_g, npw_g 
      INTEGER               :: iss_, nspin_, ngwt_, nbnd_
      INTEGER               :: nbnd_up, nbnd_dw, ntmp
      REAL(DP)              :: nelec_, ef, ef_up, ef_dw
      REAL(DP)              :: scalef_
      REAL(DP)              :: wk_(2)
      INTEGER               :: ib, nb
      REAL(DP)              :: amass_(ntypx)
      INTEGER,  ALLOCATABLE :: ityp_(:) 
      REAL(DP), ALLOCATABLE :: tau_(:,:) 
      REAL(DP), ALLOCATABLE :: occ_(:,:), et_(:,:)
      CHARACTER(LEN=80)     :: pos_unit
      REAL(DP), ALLOCATABLE :: mrepl(:,:) 
      LOGICAL               :: md_found, exist_wfc 
      INTEGER               :: io_bgrp_id
      TYPE ( output_type)   :: output_obj 
      TYPE (parallel_info_type) :: parinfo_obj
      TYPE (general_info_type ) :: geninfo_obj 
      TYPE (Node),POINTER       :: root, nodePointer
      CHARACTER(LEN=20) :: dft_name, vdw_corr
      CHARACTER(LEN=32) :: exxdiv_treatment, Hubbard_projectors
      LOGICAL :: ldftd3
      INTEGER :: nq1, nq2, nq3, lda_plus_U_kind
      REAL(dp):: exx_fraction, screening_parameter, ecutfock, ecutvcut,local_thr
      LOGICAL :: x_gamma_extrapolation
      REAL(dp):: hubbard_dum(3,nsp), hubba_dum(nsp), hubba_dum_dum(1,1,1) 
      LOGICAL :: backall_dum(nsp)
      INTEGER :: hub_l2_dum(nsp), hub_l3_dum(nsp), hub_lmax_back_dum  
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      INTEGER  :: int_dum 
      LOGICAL  :: bool_dum
      !
      !
      dirname = restart_dir(ndr)
      filename= xmlfile(ndr)
      INQUIRE ( file=filename, exist=found )
      IF (.NOT. found ) &
         CALL errore ('cp_readfile', 'xml data file not found', 1)
      !
      ! read XML file into "root" object
      !
      root => parseFile (TRIM(filename))
      !
      ! copy from "root" object into geninfo, parinfo, output objs
      !
      nodePointer => item (getElementsByTagname (root, "general_info"),0)
      ierr = 0 
      IF (ASSOCIATED(nodePointer)) THEN 
         CALL qes_read(nodePointer, geninfo_obj)
      ELSE 
         ierr = ierr + 1 
      END IF
      !
      nodePointer => item (getElementsByTagname (root, "parallel_info"),0) 
      IF (ASSOCIATED(nodePointer)) THEN 
         CALL qes_read(nodePointer, parinfo_obj)
      ELSE 
         ierr = ierr + 10
      END IF
      !
      nodePointer => item (getElementsByTagname (root, "output"),0)
      IF (ASSOCIATED(nodePointer)) THEN 
         CALL qes_read(nodePointer, output_obj)
      ELSE
         ierr = ierr + 101
      END IF
      IF ( ierr > 100) CALL errore ('cp_readfile', 'missing data in file', ierr)
      !
      ! copy CP-specific MD information directly into variables
      !
      CALL cp_readcp ( root, nat, nfi, simtime, acc, stau0, svel0, taui,  &
           cdmi, force, nhpcl, nhpdim, xnhp0, vnhp, ekincm, xnhe0, vnhe, ht,&
           htvel, gvel, xnhh0, vnhh, staum, svelm, xnhpm, xnhem, htm, xnhhm,&
           dt, ierr )
      md_found = ( ierr == 0 )
      IF ( ierr > 0 ) CALL errore ('cp_readcp','bad CP section read',ierr)
      !
      ierr = 0
      !
      ! copy Wannier function centers information directly into variables
      !
      IF ( lwf ) CALL cp_readcenters ( root, wfc)
      !
      ierr = 0
      !   
      CALL destroy (root) 
      !
      ! objects filled, now get variables from objects
      !
      CALL qexsd_copy_geninfo (geninfo_obj, qexsd_fmt, qexsd_version) 
      !
      CALL  qexsd_copy_parallel_info (parinfo_obj, nproc_file, &
            nproc_pool_file, nproc_image_file, ntask_groups_file, &
            nproc_bgrp_file, nproc_ortho_file)
      !
      CALL qexsd_copy_atomic_species (output_obj%atomic_species, nsp_, &
              atm, amass_ )
      IF ( nsp_ /= nsp ) CALL errore ('cp_readfile', 'wrong nsp read', 1)

      ALLOCATE ( tau_(3,nat), ityp_(nat) )
      CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp, &
           atm, nat_, tau_, ityp_, alat_, a1_, a2_, a3_, ibrav_, int_dum )
      IF ( nat_ /= nat ) CALL errore ('cp_readfile', 'wrong nat read', 1)
      !
      CALL recips( a1_, a2_, a3_, b1, b2, b3 )
      IF ( .not.md_found ) THEN
         ! cell not read from CP section: use cell read from xml file
         ht(1,:) = a1_
         ht(2,:) = a2_
         ht(3,:) = a3_
         !
         CALL invmat( 3, ht, htm1, omega )
         hinv = TRANSPOSE( htm1 )
         ! atomic positions not read from CP section: use those from xml file
         ! stau0 contains "scaled" atomic positions (that is, in crystal axis)
         CALL r_to_s( taui, stau0, nat_, hinv )
         CALL ions_cofmass( taui, amass_ , nat_, ityp_, cdmi )
      END IF
      !
      DEALLOCATE ( tau_, ityp_ )
      
      CALL qexsd_copy_basis_set ( output_obj%basis_set, gamma_only_, ecutwfc_,&
           ecutrho_, nr1s, nr2s, nr3s, nr1, nr2, nr3, nr1b, nr2b, nr3b, &
           ngm_g, ngms_g, npw_g, b1, b2, b3 )

      CALL qexsd_copy_dft ( output_obj%dft, nsp, atm, dft_name, &
           nq1, nq2, nq3, ecutfock, exx_fraction, screening_parameter, &
           exxdiv_treatment, x_gamma_extrapolation, ecutvcut, local_thr, &
           lda_plus_U, bool_dum,lda_plus_U_kind, Hubbard_projectors, Hubbard_n, Hubbard_l, Hubbard_lmax,Hubbard_dum, &
           hub_l2_dum, hub_l2_dum, hub_l2_dum, hub_l2_dum, backall_dum, hub_lmax_back_dum, hubba_dum, & 
           Hubbard_U, hubba_dum_dum,hubba_dum, Hubbard_dum(1,:), Hubbard_dum(2,:), hubba_dum_dum, Hubbard_dum(3,:), &
           HUBBARD_J = Hubbard_dum, HUBBARD_V = hubba_dum_dum, VDW_CORR = vdw_corr, dftd3_version = int_dum, & 
           dftd3_3body = bool_dum, SCAL6 = scal6,    & 
           LON_RCUT =lon_rcut, VDW_ISOLATED = vdw_isolated )
      CALL set_vdw_corr (vdw_corr, llondon, ldftd3, ts_vdw, mbd_vdw, lxdm )
      IF ( ldftd3 ) CALL errore('cp_readfile','DFT-D3 not implemented',1)
      !
      lsda_ = output_obj%magnetization_ispresent .AND. output_obj%magnetization%lsda
      IF ( lsda_ .AND. (nspin /= 2) ) CALL errore('cp_readfile','wrong spin',1)
      !
      nbnd_ = nupdwn(1)
      ALLOCATE( occ_(nbnd_, nspin), et_(nbnd_, nspin) )
      CALL qexsd_copy_band_structure( output_obj%band_structure, lsda_, &
           nk_, isk_, nbnd, nbnd_up, nbnd_dw, nelec_, xk, &
           wk_, occ_, ef, ef_up, ef_dw, et_ )
      ! FIXME: in the call, the same array is passed as both occ0 and occm!
      DO iss = 1, nspin
         ib = iupdwn(iss)
         nb = nupdwn(iss)
         occ0(ib:ib+nb-1) = occ_(1:nb,iss)
      END DO
      occm(:) = occ0(:)
      DEALLOCATE (occ_, et_)
      !
      DO iss = 1, nspin
         ib = iupdwn(iss)
         nb = nupdwn(iss)
         CALL cp_read_wfc( ndr, 1, 1, iss, nspin, c02, ' ' )
         CALL cp_read_wfc( ndr, 1, 1, iss, nspin, cm2, 'm', ierr )
         IF ( ierr /= 0) THEN
            CALL infomsg('cp_readfile','wfc at t-dt not found')
            cm2(:,ib:ib+nb-1) = c02(:,ib:ib+nb-1)
         END IF
         ! matrix of orthogonality constrains lambda at time t
         filename = TRIM(dirname) // 'lambda' // TRIM(int_to_char(iss))
         CALL cp_read_lambda( filename, iunpun, iss, nspin, nudx, &
              lambda0(:,:,iss), ierr )
         IF ( ierr /= 0 ) THEN
            CALL infomsg('cp_readfile','lambda not found')
            lambda0 =0.0_dp
            lambdam =0.0_dp
         ELSE
            ! matrix of orthogonality constrains lambda at time t-dt
            filename = TRIM(dirname) // 'lambdam' // TRIM(int_to_char(iss))
            CALL cp_read_lambda( filename, iunpun, iss, nspin, nudx, &
                 lambdam(:,:,iss), ierr )
         END IF
      END DO
      !
      RETURN
      !
    END SUBROUTINE cp_readfile
    !
  !------------------------------------------------------------------------
  SUBROUTINE cp_writecp( xf, nfi, simtime, &
       ekin, eht, esr, eself, epseu, enl, exc, vave, enthal, &
       acc, stau0, svel0, taui, cdmi, force, nhpcl, nhpdim, &
       xnhp0, vnhp, ekincm, xnhe0, vnhe, ht, htvel, gvel, xnhh0, vnhh,      &
       staum, svelm, xnhpm, xnhem, htm, xnhhm, dt) !
    !------------------------------------------------------------------------
    !! Cell related variables, CP-specific.
    !
    USE ions_base, ONLY: nat
    !
    TYPE(xmlf_t),  INTENT(INOUT) :: xf
    INTEGER,  INTENT(IN) :: nfi          ! index of the current step
    REAL(DP), INTENT(IN) :: simtime      ! simulated time
    REAL(DP), INTENT(IN) :: ekin, eht, esr, eself, epseu, enl, exc, vave, &
                            enthal, ekincm  ! energy terms
    REAL(DP), INTENT(IN) :: acc(:)       !  
    REAL(DP), INTENT(IN) :: stau0(:,:)
    REAL(DP), INTENT(IN) :: svel0(:,:)
    REAL(DP), INTENT(IN) :: taui(:,:)
    REAL(DP), INTENT(IN) :: cdmi(:)
    REAL(DP), INTENT(IN) :: force(:,:)
    INTEGER,  INTENT(IN) :: nhpcl
    INTEGER,  INTENT(IN) :: nhpdim
    REAL(DP), INTENT(IN) :: xnhp0(:)
    REAL(DP), INTENT(IN) :: vnhp(:)
    REAL(DP), INTENT(IN) :: xnhe0
    REAL(DP), INTENT(IN) :: vnhe
    REAL(DP), INTENT(IN) :: ht(3,3)
    REAL(DP), INTENT(IN) :: htvel(3,3)
    REAL(DP), INTENT(IN) :: gvel(3,3)
    REAL(DP), INTENT(IN) :: xnhh0(3,3)
    REAL(DP), INTENT(IN) :: vnhh(3,3)
    REAL(DP), INTENT(IN) :: staum(:,:)
    REAL(DP), INTENT(IN) :: svelm(:,:)
    REAL(DP), INTENT(IN) :: xnhpm(:)
    REAL(DP), INTENT(IN) :: xnhem
    REAL(DP), INTENT(IN) :: htm(3,3)
    REAL(DP), INTENT(IN) :: xnhhm(3,3)
    REAL(DP), INTENT(IN) :: dt !last simulation integration timestep used
    !
    !
    IF ( ionode ) THEN
!-------------------------------------------------------------------------------
! ... STATUS
!-------------------------------------------------------------------------------
       !
       CALL xml_NewElement (xf, "STATUS")
       !
       CALL xml_NewElement ( xf, "STEP")
       CALL xml_addAttribute (xf, "ITERATION", nfi) 
       CALL xml_endElement (xf, "STEP") 
       !
       CALL xml_NewElement ( xf, "TIME")
       CALL xml_addAttribute( xf, "UNITS", "pico-seconds")
       CALL xml_addCharacters( xf, simtime )
       CALL xml_EndElement(xf, "TIME")
       !
       CALL xml_NewElement ( xf, "DT")
       CALL xml_addAttribute( xf, "UNITS", "Hartree-time-atomic-units")
       CALL xml_addCharacters( xf, dt )
       CALL xml_EndElement(xf, "DT")
       !
       CALL xml_NewElement ( xf, "TITLE")
       CALL xml_addCharacters ( xf, "temporary title")
       CALL xml_EndElement (xf, "TITLE")
       !
       CALL xml_NewElement( xf, "KINETIC_ENERGY")
       CALL xml_addAttribute( xf, "UNITS", 'Hartree') 
       CALL xml_addCharacters ( xf, ekin)
       CALL xml_EndElement(xf, "KINETIC_ENERGY")
       CALL xml_NewElement( xf, "HARTREE_ENERGY")
       CALL xml_addCharacters ( xf, eht)
       CALL xml_EndElement (xf, "HARTREE_ENERGY")
       CALL xml_NewElement( xf, "EWALD_TERM")
       CALL xml_addCharacters ( xf, esr)
       CALL xml_EndElement (xf, "EWALD_TERM")
       CALL xml_NewElement( xf, "GAUSS_SELFINT")
       CALL xml_addCharacters ( xf, eself)
       CALL xml_EndElement (xf, "GAUSS_SELFINT")
       CALL xml_NewElement( xf, "LPSP_ENERGY")
       CALL xml_addCharacters ( xf, epseu)
       CALL xml_EndElement (xf, "LPSP_ENERGY")
       CALL xml_NewElement( xf, "NLPSP_ENERGY")
       CALL xml_addCharacters ( xf, enl)
       CALL xml_EndElement (xf, "NLPSP_ENERGY")
       CALL xml_NewElement( xf, "EXC_ENERGY")
       CALL xml_addCharacters ( xf, exc)
       CALL xml_EndElement (xf, "EXC_ENERGY")
       CALL xml_NewElement( xf, "AVERAGE_POT")
       CALL xml_addCharacters ( xf, vave)
       CALL xml_EndElement (xf, "AVERAGE_POT")
       CALL xml_NewElement( xf, "ENTHALPY")
       CALL xml_addCharacters ( xf, enthal)
       CALL xml_EndElement (xf, "ENTHALPY")
       !
       CALL xml_endElement( xf, "STATUS" )
       !
!-------------------------------------------------------------------------------
! ... TIMESTEPS
!-------------------------------------------------------------------------------
       !
       CALL xml_NewElement ( xf, "TIMESTEPS")
       CALL xml_addAttribute ( xf, "nt", 2)
       !
       ! ... STEP0
       !
       CALL xml_NewElement ( xf, "STEP0")
       !
       CALL xml_NewElement( xf, "ACCUMULATORS")
       CALL xml_addCharacters ( xf, acc)
       CALL xml_EndElement (xf, "ACCUMULATORS")
       !
       CALL xml_NewElement ( xf, "IONS_POSITIONS" )
       CALL xml_NewElement( xf, "stau")
       CALL xml_addCharacters ( xf, stau0(1:3,1:nat))
       CALL xml_EndElement (xf, "stau")
       CALL xml_NewElement( xf, "svel")
       CALL xml_addCharacters ( xf, svel0(1:3,1:nat) )
       CALL xml_EndElement (xf, "svel")
       CALL xml_NewElement( xf, "taui")
       CALL xml_addCharacters ( xf, taui(1:3,1:nat) )
       CALL xml_EndElement (xf, "taui")
       CALL xml_NewElement( xf, "cdmi")
       CALL xml_addCharacters ( xf, cdmi(1:3) )
       CALL xml_EndElement (xf, "cdmi")
       CALL xml_NewElement( xf, "force")
       CALL xml_addCharacters ( xf, force(1:3, 1:nat) )
       CALL xml_EndElement (xf, "force")
       CALL xml_EndElement( xf, "IONS_POSITIONS" )
       !
       CALL xml_newElement( xf, "IONS_NOSE" )
       CALL xml_NewElement( xf, "nhpcl")
       CALL xml_addCharacters ( xf, nhpcl)
       CALL xml_EndElement (xf, "nhpcl")
       CALL xml_NewElement( xf, "nhpdim")
       CALL xml_addCharacters ( xf, nhpdim)
       CALL xml_EndElement (xf, "nhpdim")
       CALL xml_NewElement( xf, "xnhp")
       CALL xml_addCharacters ( xf, xnhp0(1:nhpcl*nhpdim))
       CALL xml_EndElement (xf, "xnhp")
       CALL xml_NewElement( xf, "vnhp")
       CALL xml_addCharacters ( xf, vnhp(1:nhpcl*nhpdim) )
       CALL xml_EndElement (xf, "vnhp")
       CALL xml_EndElement( xf , "IONS_NOSE" )
       !
       CALL xml_NewElement( xf, "ekincm")
       CALL xml_addCharacters ( xf, ekincm)
       CALL xml_EndElement (xf, "ekincm")

       !
       CALL xml_NewElement ( xf, "ELECTRONS_NOSE" )
       CALL xml_NewElement( xf, "xnhe")
       CALL xml_addCharacters ( xf, xnhe0)
       CALL xml_EndElement (xf, "xnhe")
       CALL xml_NewElement( xf, "vnhe")
       CALL xml_addCharacters ( xf, vnhe)
       CALL xml_EndElement (xf, "vnhe")
       CALL xml_EndElement (  xf, "ELECTRONS_NOSE" )
       !
       CALL xml_NewElement( xf, "CELL_PARAMETERS" )
       CALL xml_NewElement( xf, "ht")
       CALL xml_addCharacters ( xf, ht)
       CALL xml_EndElement (xf, "ht")
       CALL xml_NewElement( xf, "htvel")
       CALL xml_addCharacters ( xf, htvel)
       CALL xml_EndElement (xf, "htvel")
       CALL xml_NewElement( xf, "gvel")
       CALL xml_addCharacters ( xf, gvel)
       CALL xml_EndElement (xf, "gvel")
       CALL xml_EndElement( xf, "CELL_PARAMETERS" )
       !
       CALL xml_NewElement( xf, "CELL_NOSE" )
       CALL xml_NewElement( xf, "xnhh")
       CALL xml_addCharacters ( xf, xnhh0)
       CALL xml_EndElement (xf, "xnhh")
       CALL xml_NewElement( xf, "vnhh")
       CALL xml_addCharacters ( xf, vnhh)
       CALL xml_EndElement (xf, "vnhh")
       CALL xml_EndElement(   xf, "CELL_NOSE" )
       !
       CALL xml_EndElement( xf, "STEP0" )
       !
       ! ... STEPM
       !
       CALL xml_NewElement ( xf, "STEPM" )
       !
       CALL xml_NewElement( xf, "IONS_POSITIONS" )
       CALL xml_NewElement( xf, "stau")
       CALL xml_addCharacters ( xf, staum(1:3, 1:nat) )
       CALL xml_EndElement (xf, "stau")
       CALL xml_NewElement( xf, "svel")
       CALL xml_addCharacters ( xf, svelm(1:3, 1:nat) )
       CALL xml_EndElement (xf, "svel")
       CALL xml_EndElement( xf, "IONS_POSITIONS" )
       !
       CALL xml_NewElement( xf, "IONS_NOSE" )
       CALL xml_NewElement( xf, "nhpcl")
       CALL xml_addCharacters ( xf, nhpcl)
       CALL xml_EndElement (xf, "nhpcl")
       CALL xml_NewElement( xf, "nhpdim")
       CALL xml_addCharacters ( xf, nhpdim)
       CALL xml_EndElement (xf, "nhpdim")
       CALL xml_NewElement( xf, "xnhp")
       CALL xml_addCharacters ( xf, xnhpm(1:nhpcl*nhpdim) )
       CALL xml_EndElement (xf, "xnhp")
       CALL xml_EndElement( xf, "IONS_NOSE" )
       !
       CALL xml_NewElement( xf, "ELECTRONS_NOSE" )
       CALL xml_NewElement( xf, "xnhe")
       CALL xml_addCharacters ( xf, xnhem)
       CALL xml_EndElement (xf, "xnhe")
       CALL xml_EndElement( xf, "ELECTRONS_NOSE" )
       !
       CALL xml_NewElement( xf, "CELL_PARAMETERS" )
       CALL xml_NewElement( xf, "ht")
       CALL xml_addCharacters ( xf, htm)
       CALL xml_EndElement (xf, "ht")
       CALL xml_EndElement( xf, "CELL_PARAMETERS" )
       !
       CALL xml_NewElement( xf, "CELL_NOSE" )
       CALL xml_NewElement( xf, "xnhh")
       CALL xml_addCharacters ( xf, xnhhm)
       CALL xml_EndElement (xf, "xnhh")
       CALL xml_EndElement( xf, "CELL_NOSE" )
       !
       CALL xml_EndElement( xf, "STEPM" )
       !
       CALL xml_EndElement( xf, "TIMESTEPS" )
       !
    ENDIF
    !
  END SUBROUTINE cp_writecp
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_writecenters( xf, h, wfc )
    !------------------------------------------------------------------------
    !! Write Wannier centers.
    !
    USE kinds, ONLY : dp
    USE io_global, ONLY : ionode
    USE cell_base, ONLY : ainv ! what is this? what is the relation with h?
    !
    REAL(DP), INTENT(IN) :: h(:,:), wfc(:,:)
    TYPE(xmlf_t), INTENT(inout) :: xf
    !
    INTEGER :: i, nbnd
    REAL(DP) :: temp_vec(3)
    REAL(DP), ALLOCATABLE :: centers(:,:)
    !
    IF ( ionode ) THEN
       !
       nbnd = SIZE (wfc, 2)
       ALLOCATE ( centers(3,nbnd) )
       CALL xml_NewElement( xf, "WANNIER_CENTERS" )
       !
       temp_vec=0.0_DP
       centers =0.0_DP
       !
       DO i = 1, nbnd
          !
          temp_vec(:) = MATMUL( ainv(:,:), wfc(:,i) )
          temp_vec(:) = temp_vec(:) - floor (temp_vec(:))
          centers(:,i) = MATMUL( h, temp_vec(:) )
          !
       END DO
       !
       
       CALL xml_NewElement( xf, "wanniercentres")
       CALL xml_addNewLine(xf) 
       DO i = 1, nbnd
           CALL xml_addCharacters( xf, centers(1:3,i) )
           CALL xml_addNewLine(xf)
       END DO 
       CALL xml_EndElement (xf, "wanniercentres")
       !
       DEALLOCATE ( centers )
       CALL xml_EndElement( xf, "WANNIER_CENTERS" )
       !
    END IF
    !
  END SUBROUTINE cp_writecenters
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_read_wfc( ndr, ik, nk, iss, nspin, c2, tag, ierr )
    !------------------------------------------------------------------------
    !! Wrapper, and ugly hack, for old \(\texttt{cp_read_wfc}\) called in 
    !! \(\texttt{restart.f90}\).  
    !! If \(\text{ierr}\) is present, returns (\text{ierr}=-1\) if file not
    !! found, 0 otherwise.
    !
    USE mp_bands,           ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE electrons_base,     ONLY : iupdwn, nupdwn
    USE gvecw,              ONLY : ngw, ngw_g
    USE gvect,              ONLY : ig_l2g
    !
    INTEGER,               INTENT(IN)  :: ndr
    INTEGER,               INTENT(IN)  :: ik, iss, nk, nspin
    CHARACTER,             INTENT(IN)  :: tag
    COMPLEX(DP),           INTENT(OUT) :: c2(:,:)
    INTEGER, OPTIONAL,     INTENT(OUT) :: ierr
    !
    INTEGER            :: ib, nb, nbnd, is_, npol
    INTEGER,ALLOCATABLE:: mill_k(:,:)
    CHARACTER(LEN=320) :: filename
    REAL(DP)           :: scalef, xk(3), b1(3), b2(3), b3(3)
    LOGICAL            :: gamma_only
    CHARACTER(LEN=6), EXTERNAL   :: int_to_char
    !
    filename = restart_dir(ndr)
    IF ( tag == 'm' ) THEN
       filename = TRIM(filename) //'wfcm'// TRIM(int_to_char(iss))
    ELSE
       filename = TRIM(filename) //'wfc'// TRIM(int_to_char(iss))
    END IF
    ib = iupdwn(iss)
    nb = nupdwn(iss)
    ! next two lines workaround for bogus complaint due to intent(in)
    is_= iss
    ALLOCATE ( mill_k(3,ngw) )
    !
    ! the first processor of each "band group" reads the wave function,
    ! distributes it to the other processors in the same band group
    !
    CALL read_wfc( iunpun, filename, root_bgrp, intra_bgrp_comm, &
         is_, xk, is_, npol, c2(:,ib:ib+nb-1), ngw_g, gamma_only,&
         nbnd, ig_l2g, ngw, b1,b2,b3, mill_k, scalef, ierr )
    !
    ! Add here checks on consistency of what has been read
    !
    DEALLOCATE ( mill_k)
    !
  END SUBROUTINE cp_read_wfc
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_readcp ( root, nat, nfi, simtime, acc, stau0, svel0, taui,&
       cdmi, force, nhpcl, nhpdim, xnhp0, vnhp, ekincm, xnhe0, vnhe, ht, &
       htvel, gvel, xnhh0, vnhh, staum, svelm, xnhpm, xnhem, htm, xnhhm, &
       dt, ierr )
    !
    !------------------------------------------------------------------------
    !! Cell related variables, CP-specific:
    !
    !! * ierr = -2: nothing found;
    !! * ierr = -1: MD status found, no info on timesteps;
    !! * ierr =  0: MD status and timestep info read;
    !! * ierr =  1: error reading MD status;
    !! * ierr =  2: error reading timestep info.
    !
    !
    TYPE(Node),POINTER,INTENT(IN) :: root
    INTEGER,  INTENT(IN) :: nat
    INTEGER,  INTENT(out) :: nfi
    REAL(DP), INTENT(out) :: simtime
    REAL(DP), INTENT(out) :: ekincm
    REAL(DP), INTENT(out) :: acc(:)
    REAL(DP), INTENT(out) :: stau0(:,:)
    REAL(DP), INTENT(out) :: svel0(:,:)
    REAL(DP), INTENT(out) :: taui(:,:)
    REAL(DP), INTENT(out) :: cdmi(:)
    REAL(DP), INTENT(out) :: force(:,:)
    INTEGER,  INTENT(inout) :: nhpcl
    INTEGER,  INTENT(inout) :: nhpdim
    REAL(DP), INTENT(out) :: xnhp0(:)
    REAL(DP), INTENT(out) :: vnhp(:)
    REAL(DP), INTENT(out) :: xnhe0
    REAL(DP), INTENT(out) :: vnhe
    REAL(DP), INTENT(out) :: ht(3,3)
    REAL(DP), INTENT(out) :: htvel(3,3)
    REAL(DP), INTENT(out) :: gvel(3,3)
    REAL(DP), INTENT(out) :: xnhh0(3,3)
    REAL(DP), INTENT(out) :: vnhh(3,3)
    REAL(DP), INTENT(out) :: staum(:,:)
    REAL(DP), INTENT(out) :: svelm(:,:)
    REAL(DP), INTENT(out) :: xnhpm(:)
    REAL(DP), INTENT(out) :: xnhem
    REAL(DP), INTENT(out) :: htm(3,3)
    REAL(DP), INTENT(out) :: xnhhm(3,3)
    REAL(DP), INTENT(out) :: dt
    INTEGER,  INTENT(out) :: ierr
    !
    LOGICAL :: found
    INTEGER :: nt_, nhpcl_, nhpdim_
    TYPE(Node), POINTER      :: n1Pointer, n2Pointer, n3Pointer, n4Pointer
    !
    ! ... read MD status
    !
    ierr = -2
    n1Pointer => item( getElementsByTagname( root, "STATUS"), 0)
    found = ASSOCIATED( n1Pointer)
    IF ( .NOT.found ) RETURN
    !
    ierr = 1
    n2Pointer => item( getElementsByTagname( n1Pointer, "STEP"), 0)
    found = ASSOCIATED( n2Pointer)
    IF ( .NOT.found ) RETURN
    !
    found = hasAttribute( n2Pointer, "ITERATION") 
    IF (.NOT. found ) RETURN
    CALL extractDataAttribute( n2Pointer, "ITERATION", nfi)
    !
    n2Pointer => item( getElementsByTagname( n1Pointer, "TIME"), 0)
    found = ASSOCIATED(n2Pointer)
    IF ( .NOT.found ) RETURN
    CALL extractDataContent( n2Pointer, simtime)
    !
    n2Pointer => item( getElementsByTagname( n1Pointer, "DT"), 0)
    found = ASSOCIATED(n2Pointer)
    IF ( .NOT.found ) THEN
       dt = -1.0_dp
    ELSE
       CALL extractDataContent( n2Pointer, dt)
    END IF
    !
    ! ... read MD timesteps variables
    !
    n1Pointer => item( getElementsByTagname( root, "TIMESTEPS"), 0)
    found = ASSOCIATED( n1Pointer)
    ! 
    IF ( found ) THEN
       !
       ierr = 0
       !
       CALL extractDataAttribute( n1Pointer, "nt", nt_)
       !
       IF ( nt_ > 0 ) THEN
          !
          n2Pointer => item( getElementsByTagname( n1Pointer, "STEP0"), 0)
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "ACCUMULATORS"), 0)
          CALL extractDataContent(n3Pointer, acc)
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "IONS_POSITIONS"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "stau"), 0)
          CALL extractDataContent(n4Pointer, stau0(1:3,1:nat) )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "svel"), 0)
          CALL extractDataContent(n4Pointer, svel0(1:3,1:nat) )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "taui"), 0)
          CALL extractDataContent(n4Pointer, taui(1:3,1:nat) )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "cdmi"), 0)
          CALL extractDataContent(n4Pointer, cdmi(1:3) )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "force"), 0)
          CALL extractDataContent(n4Pointer, force(1:3,1:nat) )
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "IONS_NOSE"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "nhpcl"), 0)
          CALL extractDataContent(n4Pointer, nhpcl_ )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "nhpdim"), 0)
          CALL extractDataContent(n4Pointer, nhpdim_ )
          IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
             n4Pointer =>item( getElementsByTagname( n3Pointer, "xnhp"), 0)
             CALL extractDataContent(n4Pointer, xnhp0(1:nhpcl*nhpdim ) )
             n4Pointer =>item( getElementsByTagname( n3Pointer, "vnhp"), 0)
             CALL extractDataContent(n4Pointer, vnhp(1:nhpcl*nhpdim ) )
          ELSE
             xnhp0(1:nhpcl*nhpdim) = 0.D0
             vnhp(1:nhpcl*nhpdim)  = 0.D0
          END IF
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "ekincm"), 0)
          CALL extractDataContent(n3Pointer, ekincm) 
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "ELECTRONS_NOSE"), 0) 
          n4Pointer =>item( getElementsByTagname( n3Pointer, "xnhe"), 0)
          CALL extractDataContent(n4Pointer, xnhe0 )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "vnhe"), 0)
          CALL extractDataContent(n4Pointer, vnhe )
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "CELL_PARAMETERS"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "ht"), 0)
          CALL extractDataContent(n4Pointer, ht )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "htvel"), 0)
          CALL extractDataContent(n4Pointer, htvel )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "gvel"), 0)
          CALL extractDataContent(n4Pointer, gvel )
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "CELL_NOSE"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "xnhh"), 0)
          CALL extractDataContent(n4Pointer, xnhh0 )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "vnhh"), 0)
          CALL extractDataContent(n4Pointer, vnhh )
          !
          !
       ELSE
          !
          ierr = 2
          RETURN
          !
       END IF
       !
       IF ( nt_ > 1 ) THEN
          !
          n2Pointer => item( getElementsByTagname( n1Pointer, "STEPM"), 0)
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "IONS_POSITIONS"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "stau"), 0)
          CALL extractDataContent(n4Pointer, staum(1:3,1:nat) )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "svel"), 0)
          CALL extractDataContent(n4Pointer, svelm(1:3,1:nat) )
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "IONS_NOSE"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "nhpcl"), 0)
          CALL extractDataContent(n4Pointer, nhpcl_ )
          n4Pointer =>item( getElementsByTagname( n3Pointer, "nhpdim"), 0)
          CALL extractDataContent(n4Pointer, nhpdim_ )
          !
          IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
             n4Pointer =>item( getElementsByTagname( n3Pointer, "xnhp"), 0)
             CALL extractDataContent(n4Pointer, xnhpm(1:nhpcl*nhpdim))
          ELSE
             xnhpm(1:nhpcl*nhpdim) = 0.D0
          END IF
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "ELECTRONS_NOSE"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "xnhe"), 0)
          CALL extractDataContent(n4Pointer, xnhem )
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "CELL_PARAMETERS"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "ht"), 0)
          CALL extractDataContent(n4Pointer, htm )
          !
          n3Pointer => item( getElementsByTagname( n2Pointer, "CELL_NOSE"), 0)
          n4Pointer =>item( getElementsByTagname( n3Pointer, "xnhh"), 0)
          CALL extractDataContent(n4Pointer, xnhhm )
          !
          !
       END IF
       !
       !
    ELSE
       !
       ierr = -1
       !
       ! ... MD time steps not found, try to recover from CELL and POSITIONS
       ! 
       acc = 0.D0
       ! 
       staum = stau0
       svel0 = 0.D0
       svelm = 0.D0
       force = 0.D0
       !
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
    END IF
    !
  END SUBROUTINE cp_readcp
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_readcenters( root, wfc )
    !------------------------------------------------------------------------
    !! Read Wannier centers.
    !
    USE kinds, ONLY : dp
    USE io_global, ONLY : stdout
    !
    TYPE(Node), POINTER, INTENT(IN)  :: root
    REAL(DP), INTENT(OUT):: wfc(:,:)
    !
    INTEGER :: nbnd, ierr
    LOGICAL :: found
    TYPE(Node), POINTER      :: n1Pointer, n2Pointer, n3Pointer
    !
    nbnd = SIZE (wfc, 2)
    n1Pointer => item( getElementsByTagname(root, "WANNIER_CENTERS"),0) 
    found = ASSOCIATED ( n1Pointer ) 
    !
    IF (found) THEN
       !
       n2Pointer => item( getElementsByTagname ( n1Pointer, "wanniercentres"),0)
       IF (ASSOCIATED (n2Pointer) ) THEN 
          CALL extractDataContent ( n2Pointer, wfc(1:3,1:nbnd), IOSTAT = ierr )
       ELSE
          ierr = 210
       END IF 
       IF ( ierr > 0 ) CALL errore ('cp_readcenters', &
            'error reading Wannier centers',ierr)
       !
    ELSE
       !
       CALL infomsg('cp_readcenters', &
            'Wannier centers not found in restart file')
       wfc(:,:)= 0.0_dp
       !
    END IF
    !
  END SUBROUTINE cp_readcenters
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_read_cell( ndr, ascii, ht, &
                           htm, htvel, gvel, xnhh0, xnhhm, vnhh )
    !------------------------------------------------------------------------
    !
    USE parameters,  ONLY : ntypx
    USE ions_base,   ONLY : nat
    USE qexsd_copy,  ONLY : qexsd_copy_atomic_structure
    USE qes_read_module, ONLY : qes_read
    !
    INTEGER,          INTENT(IN)    :: ndr
    LOGICAL,          INTENT(IN)    :: ascii
    REAL(DP),         INTENT(INOUT) :: ht(3,3)
    REAL(DP),         INTENT(INOUT) :: htm(3,3)
    REAL(DP),         INTENT(INOUT) :: htvel(3,3)
    REAL(DP),         INTENT(INOUT) :: gvel(3,3)
    REAL(DP),         INTENT(INOUT) :: xnhh0(3,3)
    REAL(DP),         INTENT(INOUT) :: xnhhm(3,3)
    REAL(DP),         INTENT(INOUT) :: vnhh(3,3)
    !
    CHARACTER(LEN=320) :: filename
    INTEGER            :: strlen
    INTEGER            :: i, ierr, nt_
    LOGICAL            :: found
    !
    ! ... variables read for testing pourposes
    !
    INTEGER          :: ibrav_
    INTEGER          :: nat_
    INTEGER          :: nsp_
    INTEGER, ALLOCATABLE  :: ityp_(:) 
    REAL(DP)         :: alat_
    REAL(DP)         :: a1_(3), a2_(3), a3_(3)
    REAL(DP)         :: b1_(3), b2_(3), b3_(3)
    REAL(DP), ALLOCATABLE :: tau_(:,:) 
    INTEGER          :: natomwfc_
    CHARACTER(LEN=6) :: atm_(ntypx)
    TYPE(output_type) :: output_obj
    TYPE(Node),POINTER :: root, simpleNode, timestepsNode, cellNode, stepNode
    !
    !
    filename = xmlfile (ndr)
    INQUIRE ( file=filename, exist=found )
    IF (.NOT. found ) &
         CALL errore ('cp_read_cell', 'xml data file not found', 1)
    !
    root => parseFile(filename) 
    !
    timestepsNode => item(getElementsByTagname(root, "TIMESTEPS"),0)
    found = ASSOCIATED(timestepsNode)
    !
    ierr = 0
    IF ( found ) THEN
       !
       CALL extractDataAttribute(timestepsNode, "nt", nt_)
       !
       IF ( nt_ > 0 ) THEN
          !
          stepNode => item ( getElementsByTagname(timestepsNode,"STEP0"),0)
          !
          cellNode => item ( getElementsByTagname(stepNode, "CELL_PARAMETERS"),0)
          !
          simpleNode => item( getElementsByTagname(cellNode, "ht"),0)
          CALL extractDataContent(simpleNode,ht)
          !
          simpleNode => item( getElementsByTagname(cellNode, "htvel"),0)
          CALL extractDataContent(simpleNode,htvel)
          !
          simpleNode => item( getElementsByTagname(cellNode, "gvel"),0)
          IF (ASSOCIATED(simpleNode)) THEN 
              CALL extractDataContent(simpleNode,gvel)
          ELSE 
              gvel = 0.d0
          END IF
          !
          !
          cellNode => item ( getElementsByTagname(stepNode, "CELL_NOSE"),0)
          simpleNode => item( getElementsByTagname(cellNode, "xnhh"),0)
          CALL extractDataContent(simpleNode,xnhh0)
          !
          simpleNode => item( getElementsByTagname(cellNode, "vnhh"),0)
          CALL extractDataContent(simpleNode,vnhh) 
          !
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
          stepNode => item ( getElementsByTagname(timestepsNode, "STEPM"),0)
          !
          cellNode => item ( getElementsByTagname(stepNode, "CELL_PARAMETERS"),0)
          !
          simpleNode => item( getElementsByTagname(cellNode, "ht"),0)
          CALL extractDataContent(simpleNode,htm)
          !
          cellNode => item ( getElementsByTagname(stepNode, "CELL_NOSE"),0)
          !
          simpleNode => item( getElementsByTagname(cellNode, "xnhh"),0)
          CALL extractDataContent(simpleNode,xnhhm)
          !
       END IF
       !
       !
    ELSE
       !
       ! ... MD steps have not been found, try to restart from cell data
       !
       simpleNode  => item ( getElementsByTagname(root, "output"),0)
       CALL qes_read(simpleNode, output_obj) 
       !
       CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp_, &
            atm_, nat_, tau_, ityp_, alat_, a1_, a2_, a3_, ibrav_ , natomwfc_)
       IF ( nat_ /= nat ) CALL errore ('cp_readfile', 'wrong nat read', 1)
       CALL qes_reset (output_obj)
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
    CALL destroy (root)
    !
100 CALL errore( 'cp_read_cell ', 'error reading MD steps', ierr )
    !
  END SUBROUTINE cp_read_cell

  !------------------------------------------------------------------------
  SUBROUTINE cp_write_lambda( filename, iunpun, iss, nspin, nudx, &
       lambda, ierr )
    !------------------------------------------------------------------------
    !! Collect and write matrix lambda to file.
    !
    USE kinds, ONLY : dp
    USE mp, ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    USE io_global, ONLY : ionode, ionode_id
    USE cp_main_variables, ONLY : idesc
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: iunpun, iss, nspin, nudx
    REAL(dp), INTENT(in) :: lambda(:,:)
    INTEGER, INTENT(out) :: ierr
    !
    include 'laxlib.fh'
    !
    REAL(dp), ALLOCATABLE :: mrepl(:,:)
    !
    IF ( ionode ) OPEN( unit=iunpun, file =TRIM(filename), &
         status='unknown', form='unformatted', iostat=ierr)
    CALL mp_bcast (ierr, ionode_id, intra_image_comm )
    IF ( ierr /= 0 ) RETURN
    !
    ALLOCATE( mrepl( nudx, nudx ) )
    CALL collect_lambda( mrepl, lambda, idesc(:,iss) )
    !
    IF ( ionode ) THEN
       WRITE (iunpun, iostat=ierr) mrepl
       CLOSE( unit=iunpun, status='keep')
    END IF
    CALL mp_bcast (ierr, ionode_id, intra_image_comm )
    DEALLOCATE( mrepl )
    !
  END SUBROUTINE cp_write_lambda
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_read_lambda( filename, iunpun, iss, nspin, nudx, &
             lambda, ierr )
    !------------------------------------------------------------------------
    !! Read matrix lambda from file, distribute it.
    !
    USE kinds, ONLY : dp
    USE mp, ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    USE io_global, ONLY : ionode, ionode_id
    USE cp_main_variables, ONLY : idesc
    !
    include 'laxlib.fh'

    CHARACTER(LEN=*), INTENT(in) :: filename
    INTEGER, INTENT(in) :: iunpun, iss, nspin, nudx
    REAL(dp), INTENT(out) :: lambda(:,:)
    INTEGER, INTENT(out) :: ierr
    !
    LOGICAL :: exst
    REAL(dp), ALLOCATABLE :: mrepl(:,:)
    !
    ierr =0
    IF (ionode) INQUIRE( file =TRIM(filename), exist=exst )
    CALL mp_bcast (exst, ionode_id, intra_image_comm )
    IF (.NOT. exst) THEN
       ierr =-1
       RETURN
    END IF
    !
    ALLOCATE( mrepl( nudx, nudx ) )
    IF (ionode) THEN
       OPEN( unit=iunpun, file =TRIM(filename), status='old', &
            form='unformatted')
       READ (iunpun, iostat=ierr) mrepl
       CLOSE( unit=iunpun, status='keep')
    END IF
    CALL mp_bcast( mrepl, ionode_id, intra_image_comm )
    CALL distribute_lambda( mrepl, lambda, idesc(:,iss) )
    DEALLOCATE( mrepl )
    !
  END SUBROUTINE cp_read_lambda
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_write_zmat( ndw, mat_z, ierr )
    !------------------------------------------------------------------------
    !! Collect and write matrix z to file.
    !
    USE kinds, ONLY : dp
    USE mp, ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    USE io_global, ONLY : ionode, ionode_id
    USE cp_main_variables, ONLY : idesc
    USE electrons_base,ONLY: nspin, nudx
    !
    include 'laxlib.fh'

    REAL(dp), INTENT(in) :: mat_z(:,:,:)
    INTEGER, INTENT(in)  :: ndw
    INTEGER, INTENT(out) :: ierr
    !
    CHARACTER(LEN=256)    :: dirname
    CHARACTER(LEN=320)    :: filename
    INTEGER               :: iss
    REAL(dp), ALLOCATABLE :: mrepl(:,:)
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    dirname = restart_dir(ndw) 
    !
    ALLOCATE( mrepl( nudx, nudx ) )
    !
    DO iss = 1, nspin
       !
       CALL collect_zmat( mrepl, mat_z(:,:,iss), idesc(:,iss) )
       !
       filename = TRIM(dirname) // 'mat_z' // TRIM(int_to_char(iss))
       !
       IF ( ionode ) THEN
          OPEN( unit=iunpun, file =TRIM(filename), &
              status='unknown', form='unformatted', iostat=ierr)
          WRITE (iunpun, iostat=ierr) mrepl
          CLOSE( unit=iunpun, status='keep')
       END IF
       !
       CALL mp_bcast (ierr, ionode_id, intra_image_comm )
       !
    END DO
    !
    DEALLOCATE( mrepl )
    !
  END SUBROUTINE cp_write_zmat  
  !------------------------------------------------------------------------
  SUBROUTINE cp_read_zmat( ndr, mat_z, ierr )
    !------------------------------------------------------------------------
    !! Read from file and distribute matrix z.
    !
    USE kinds, ONLY : dp
    USE mp, ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    USE io_global, ONLY : ionode, ionode_id
    USE cp_main_variables, ONLY : idesc
    USE electrons_base,ONLY: nspin, nudx
    !
    include 'laxlib.fh'

    REAL(dp), INTENT(out) :: mat_z(:,:,:)
    INTEGER, INTENT(in)  :: ndr
    INTEGER, INTENT(out) :: ierr
    !
    CHARACTER(LEN=256)    :: dirname
    CHARACTER(LEN=320)    :: filename
    INTEGER               :: iss
    REAL(dp), ALLOCATABLE :: mrepl(:,:)
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    dirname = restart_dir(ndr) 
    !
    ALLOCATE( mrepl( nudx, nudx ) )
    !
    DO iss = 1, nspin
       !
       filename = TRIM(dirname) // 'mat_z' // TRIM(int_to_char(iss))
       IF ( ionode ) THEN
          OPEN( unit=iunpun, file =TRIM(filename), &
            status='old', form='unformatted', iostat=ierr)
          READ (iunpun, iostat=ierr) mrepl
          CLOSE( unit=iunpun, status='keep')
       END IF
       CALL mp_bcast (ierr, ionode_id, intra_image_comm )
       !
       CALL distribute_zmat( mrepl, mat_z(:,:,iss), idesc(:,iss) )
       !
    END DO
    !
    DEALLOCATE( mrepl )
    ! not sure about the following line
    ! CALL mp_bcast( mat_z(:,:,:), ionode_id, intra_image_comm )
    !
  END SUBROUTINE cp_read_zmat
  !------------------------------------------------------------------------
  !
END MODULE cp_restart_new
