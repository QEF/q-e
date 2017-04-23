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
  !
  ! ... This module contains subroutines to write and read data required to
  ! ... restart a calculation from the disk
  !
#if !defined(__OLDXML)
  !
  USE iotk_module
  USE qes_module
  USE qexsd_input, ONLY: qexsd_init_k_points_ibz
  USE qexsd_module, ONLY: qexsd_init_schema, qexsd_openschema, qexsd_closeschema,      &
                          qexsd_init_convergence_info, qexsd_init_algorithmic_info,    & 
                          qexsd_init_atomic_species, qexsd_init_atomic_structure,      &
                          qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
                          qexsd_init_magnetization,qexsd_init_band_structure,          &
                          qexsd_init_dipole_info, qexsd_init_total_energy,             &
                          qexsd_init_forces,qexsd_init_stress,                         &
                          qexsd_init_outputElectricField, input_obj => qexsd_input_obj
  USE io_files,  ONLY : iunpun, xmlpun_schema, prefix, tmp_dir, qexsd_fmt,&
       qexsd_version
  USE io_base,   ONLY : write_wfc, read_wfc
  USE xml_io_base,     ONLY  : write_rho_xml,read_print_counter, create_directory
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE parser,    ONLY : version_compare
  USE matrix_inversion
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
      USE cell_base,                ONLY : ibrav, alat, s_to_r, ainv ! BS added ainv
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           amass, iforce, ind_bck
      USE funct,                    ONLY : get_dft_name, get_inlc, &
           dft_is_hybrid, get_exx_fraction, get_screening_parameter, &
           dft_is_nonlocc, get_nonlocc_name
      USE ldaU_cp,                  ONLY : lda_plus_U, ns, Hubbard_l, &
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
      USE wrappers,                 ONLY : f_copy
      USE uspp,                     ONLY : okvan
      USE input_parameters,         ONLY : vdw_corr, london, starting_ns_eigenvalue
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
      INTEGER               :: kunit, ik_eff
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
      INTEGER               :: natomwfc, nbnd_, nb, ib
      REAL(DP), ALLOCATABLE :: mrepl(:,:)
      CHARACTER(LEN=256)    :: tmp_dir_save
      LOGICAL               :: exst
      INTEGER               :: inlc
      REAL(DP), ALLOCATABLE :: temp_vec(:), wfc_temp(:,:) ! BS 
      TYPE(output_type) :: output_obj
      LOGICAL :: is_hubbard(nsp)
      REAL(dp):: hubbard_dum(3,nsp)
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
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
      ! XML descriptor
      ! 
      WRITE(dirname,'(A,A,"_",I2,".save/")') TRIM(tmp_dir), TRIM(prefix), ndw
      CALL create_directory( TRIM(dirname) )
      !
      CALL qexsd_init_schema( iunpun )
      !
      IF ( ionode ) THEN
         !
         ! ... here we init the variables and finally write them to file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL qexsd_openschema(TRIM( dirname ) // TRIM( xmlpun_schema ))
         output_obj%tagname="output"
         output_obj%lwrite = .TRUE.
!-------------------------------------------------------------------------------
! ... CP-SPECIFIC CELL variables
!-------------------------------------------------------------------------------
         !
         CALL cp_writecp( iunpun, nfi, simtime, ekin, eht, esr, eself, &
              epseu, enl, exc, vave, enthal, acc, stau0, svel0, taui, cdmi,&
              force, nhpcl, nhpdim, xnhp0, vnhp, ekincm, xnhe0, vnhe, ht,&
              htvel, gvel, xnhh0, vnhh, staum, svelm, xnhpm, xnhem, htm, xnhhm)
         !
!-------------------------------------------------------------------------------
! ... CONVERGENCE_INFO - TO BE VERIFIED
!-------------------------------------------------------------------------------
!
         CALL qexsd_init_convergence_info(output_obj%convergence_info, &
              n_scf_steps=0, scf_error=0.0_dp, &
              opt_conv_ispresent=.FALSE., &
              n_opt_steps=0, grad_norm=0.0_dp )
         !
!-------------------------------------------------------------------------------
! ... ALGORITHMIC_INFO
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_algorithmic_info(output_obj%algorithmic_info, &
              real_space_q=.FALSE., uspp=okvan, paw=.FALSE.)
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
              nat, tau(:,ind_bck(:)), alat, alat*a1(:), alat*a2(:), alat*a3(:), ibrav)
         !
!-------------------------------------------------------------------------------
! ... SYMMETRIES
!-------------------------------------------------------------------------------
         output_obj%symmetries%lwrite=.false.
!-------------------------------------------------------------------------------
! ... BASIS SET
!-------------------------------------------------------------------------------
         CALL qexsd_init_basis_set(output_obj%basis_set,gamma_only, ecutwfc/e2, ecutwfc*dual/e2, &
              dfftp%nr1, dfftp%nr2, dfftp%nr3, dffts%nr1, dffts%nr2, dffts%nr3, &
              .FALSE., dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, ngms_g, ngw_g, &
              b1(:), b2(:), b3(:) )
!-------------------------------------------------------------------------------
! ... XC FUNCTIONAL
!-------------------------------------------------------------------------------
         dft_name = get_dft_name()
         is_hubbard(:) = (Hubbard_U(:) > 0.0_dp)
         hubbard_dum(:,:)= 0.0_dp
         CALL qexsd_init_dft(output_obj%dft, dft_name, .true., dft_is_hybrid(), &
              0, 0, 0, ecutwfc, get_exx_fraction(), get_screening_parameter(),&
              'none', .false., 0.0_dp, &
              dft_is_nonlocc(), TRIM(vdw_corr), &
              TRIM ( get_nonlocc_name()), scal6, in_c6, lon_rcut, 0.0_dp, &
              0.0_dp, vdw_econv_thr, vdw_isolated, &
              lda_plus_u, 0, 2*Hubbard_lmax+1, .false.,&
              nspin, nsp, nat, atm, ityp, Hubbard_U,&
              Hubbard_dum(1,:), Hubbard_dum(2,:), Hubbard_dum(3,:),Hubbard_dum,&
              starting_ns_eigenvalue, 'atomic', is_hubbard, upf(1:nsp)%psd, ns )
!-------------------------------------------------------------------------------
! ... MAGNETIZATION
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_magnetization(output_obj%magnetization, lsda, .false.,&
              .false., 0.0_dp, [0.0_dp,0.0_dp, 0.0_dp], 0.0_dp, .false.)
         !
!-------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------
         CALL  qexsd_init_total_energy(output_obj%total_energy,enthal, 0.0_dp, eht,&
              vave, exc, 0.0_dp, 0.0_dp, 0.0_dp)
!-------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------
         ! TEMP
         CALL qexsd_init_k_points_ibz( input_obj%k_points_ibz, 'Gamma', &
              'CP',nk1,nk2,nk3,k1,k2,k3,1,xk,wk,alat,a1,.false.) 
         input_obj%bands%occupations%tagname="occupations"
         input_obj%bands%occupations%lread=.false.
         input_obj%bands%occupations%lwrite=.true.
         input_obj%bands%occupations%spin_ispresent=.false.
         input_obj%bands%occupations%occupations="fixed"
         ! TEMP
         CALL  qexsd_init_band_structure(output_obj%band_structure,lsda, .false., &
              .false., nbnd_, nelec, natomwfc, .true., 0.0_dp , .false., & 
              [0.0_dp,0.0_dp], et,  DBLE( ftmp ), nspin, xk, [ngw_g], wk,  &
              STARTING_KPOINTS = input_obj%k_points_IBZ, &
              OCCUPATION_KIND = input_obj%bands%occupations, &
              WF_COLLECTED = twfcollect)
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
         ! may be wrong or incomplete
         IF ( tpre) h = -MATMUL( detot, ht ) / omega
         CALL qexsd_init_stress(output_obj%stress, h, tpre ) 
!-------------------------------------------------------------------------------
! ... ACTUAL WRITING
!-------------------------------------------------------------------------------
         !
         CALL qes_write_output(iunpun,output_obj)
         CALL qes_reset_output(output_obj)
         !
!-------------------------------------------------------------------------------
! ... CLOSING
!-------------------------------------------------------------------------------
         !
         CALL qexsd_closeschema()
         !
      END IF
      !
!-------------------------------------------------------------------------------
! ... WRITE WFC
!-------------------------------------------------------------------------------
      DO iss = 1, nspin
         !
         ik_eff = iss
         filename = TRIM(dirname) // 'wfc' // TRIM(int_to_char(ik_eff))
         ib = iupdwn(iss)
         nb = nupdwn(iss)
         CALL write_wfc( iunpun, ik_eff, nk, iss, nspin, &
              c02(:,ib:ib+nb-1), ngw_g, gamma_only, nb, ig_l2g, ngw,  &
              filename, scalef, ionode, root_pool, intra_pool_comm )
         !
      END DO
!-------------------------------------------------------------------------------
! ... WRITE PSEUDOPOTENTIALS
!-------------------------------------------------------------------------------
     !
     ! ... copy pseudopotential files into the .save directory
     !
     DO is = 1, nsp
        sourcefile= TRIM(pseudo_dir)//psfile(is)
        filename  = TRIM(dirname)//psfile(is)
        IF ( TRIM(sourcefile) /= TRIM(filename) ) &
             ierr = f_copy(sourcefile, filename)
     END DO
     inlc = get_inlc()
     IF ( inlc > 0 ) THEN 
        sourcefile= TRIM(kernel_file_name)
        filename = TRIM(dirname)//TRIM(vdw_table_name)
        IF ( TRIM(sourcefile) /= TRIM(filename) ) & 
           ierr = f_copy(sourcefile, filename)
     END IF  
     !
!-------------------------------------------------------------------------------
! ... CHARGE DENSITY
!-------------------------------------------------------------------------------
      !
      IF (write_charge_density) then
         !
         filename = TRIM( dirname ) // 'charge-density'
         !
         IF ( nspin == 1 ) THEN
            !
            CALL write_rho_xml( filename, rho(:,1), &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, &
                                dfftp%ipp, dfftp%npp, ionode, intra_bgrp_comm, inter_bgrp_comm )
            !
         ELSE IF ( nspin == 2 ) THEN
            !
            ALLOCATE( rhoaux( SIZE( rho, 1 ) ) )
            !
            rhoaux = rho(:,1) + rho(:,2) 
            !
            CALL write_rho_xml( filename, rhoaux, &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, &
                                dfftp%ipp, dfftp%npp, ionode, intra_bgrp_comm, inter_bgrp_comm )
            !
            filename = TRIM( dirname ) // 'spin-polarization'
            !
            rhoaux = rho(:,1) - rho(:,2) 
            !
            CALL write_rho_xml( filename, rhoaux, &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, &
                                dfftp%ipp, dfftp%npp, ionode, intra_bgrp_comm, inter_bgrp_comm )
            !
            DEALLOCATE( rhoaux )
            !
         END IF
         !
      END IF ! write_charge_density

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
      USE control_flags,            ONLY : gamma_only, force_pairing, llondon,&
                                           ts_vdw, lxdm, iverbosity, twfcollect, lwf
      USE io_files,                 ONLY : iunpun, xmlpun, iunwfc, nwordwfc, &
                                           tmp_dir, diropn
      USE run_info,                 ONLY : title
      USE gvect,                    ONLY : ngm
      USE gvecw,                    ONLY : ngw, ngw_g
      USE electrons_base,           ONLY : nspin, nbnd, nelt, nel, &
                                           nupdwn, iupdwn, nudx
      USE cell_base,                ONLY : ibrav, alat, s_to_r, r_to_s
      USE ions_base,                ONLY : nsp, nat, na, atm, zv, &
                                           sort_tau, ityp, ions_cofmass
      USE gvect,       ONLY : ig_l2g, mill
      USE cp_main_variables,        ONLY : nprint_nfi, descla
      USE cp_interfaces,            ONLY : distribute_lambda, distribute_zmat
      USE ldaU_cp,                  ONLY : lda_plus_U, ns, Hubbard_l, &
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
      USE qexsd_reader_module,      ONLY : qexsd_get_general_info, &
           qexsd_get_parallel_info, qexsd_get_output
      USE kernel_table,             ONLY : vdw_table_name
      USE london_module,            ONLY : scal6, lon_rcut, in_c6
      USE tsvdw_module,             ONLY : vdw_isolated, vdw_econv_thr
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
      INTEGER               :: iss_, nspin_, ngwt_, nbnd_ , nbnd_tot
      INTEGER               :: nstates_up_ , nstates_dw_ , ntmp, nel_(2)
      REAL(DP)              :: nelec_, ef, ef_up, ef_dw
      REAL(DP)              :: scalef_
      REAL(DP)              :: wk_(2)
      INTEGER               :: ib, nb
      INTEGER               :: ik_eff
      REAL(DP)              :: amass_(ntypx)
      INTEGER,  ALLOCATABLE :: ityp_(:) 
      INTEGER,  ALLOCATABLE :: isrt_(:) 
      REAL(DP), ALLOCATABLE :: tau_(:,:) 
      REAL(DP), ALLOCATABLE :: occ_(:,:), et_(:,:)
      INTEGER,  ALLOCATABLE :: if_pos_(:,:) 
      CHARACTER(LEN=256)    :: psfile_(ntypx)
      CHARACTER(LEN=80)     :: pos_unit
      REAL(DP)              :: s1, s0, cclock
      REAL(DP), ALLOCATABLE :: mrepl(:,:) 
      LOGICAL               :: md_found, exist_wfc 
      CHARACTER(LEN=256)    :: tmp_dir_save
      INTEGER               :: io_bgrp_id
      TYPE ( output_type)   :: output_obj 
      TYPE (parallel_info_type) :: parinfo_obj
      TYPE (general_info_type ) :: geninfo_obj 
      CHARACTER(LEN=20) :: dft_name
      CHARACTER(LEN=32) :: exxdiv_treatment, U_projection
      CHARACTER(LEN=256):: vdw_corr
      INTEGER :: nq1, nq2, nq3, lda_plus_U_kind, inlc
      REAL(dp):: ecutfock, exx_fraction, screening_parameter, ecutvcut
      LOGICAL :: x_gamma_extrapolation
      REAL(dp):: hubbard_dum(3,nsp)
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ! ... look for an empty unit
      !
      CALL iotk_free_unit( iunpun, ierr )
      CALL errore( 'cp_readfile', &
                   'no free units to read wavefunctions', ierr )
      !
      CALL qexsd_init_schema( iunpun )
      !
      WRITE(dirname,'(A,A,"_",I2,".save/")') TRIM(tmp_dir), TRIM(prefix), ndr
      filename = TRIM( dirname ) // TRIM( xmlpun_schema )
      INQUIRE ( file=filename, exist=found )
      IF (.NOT. found ) &
         CALL errore ('cp_readfile', 'xml data file not found', 1)
      !
      CALL iotk_open_read( iunpun, TRIM(filename) )
      !
      CALL cp_readcp ( iunpun, nat, nfi, simtime, acc, stau0, svel0, taui,  &
           cdmi, force, nhpcl, nhpdim, xnhp0, vnhp, ekincm, xnhe0, vnhe, ht,&
           htvel, gvel, xnhh0, vnhh, staum, svelm, xnhpm, xnhem, htm, xnhhm,&
           ierr )
      md_found = ( ierr == 0 )
      IF ( ierr > 0 ) CALL errore ('cp_readcp','bad CP section read',ierr)
      !
      CALL qexsd_get_general_info ( iunpun, geninfo_obj, found)
      IF ( .NOT. found ) THEN
         ierr = ierr + 1
      ELSE
         CALL qexsd_copy_general_info (geninfo_obj, qexsd_fmt, qexsd_version) 
      END IF
      !
      CALL qexsd_get_parallel_info ( iunpun, parinfo_obj, found ) 
      IF ( .NOT. found ) THEN
         ierr = ierr + 10
      ELSE
         CALL  qexsd_copy_parallel_info (parinfo_obj, nproc_file, &
              nproc_pool_file, nproc_image_file, ntask_groups_file, &
              nproc_bgrp_file, nproc_ortho_file)
      END IF
      !
      CALL qexsd_get_output ( iunpun, output_obj, found ) 
      IF ( .NOT. found ) ierr = ierr + 101
      IF ( ierr > 100) CALL errore ('cp_readfile', 'missing data in file', ierr)
      !
      CALL qexsd_copy_atomic_species (output_obj%atomic_species, nsp_, atm, &
           psfile_, amass_)
      IF ( nsp_ /= nsp ) CALL errore ('cp_readfile', 'wrong nsp read', 1)

      ALLOCATE ( tau_(3,nat), ityp_(nat), isrt_(nat) )
      CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp, &
           atm, nat_, tau_, ityp_, alat_, a1_, a2_, a3_, ibrav_ )
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
         ! reorder atomic positions according to CP (il-)logic (output in taui)
         CALL sort_tau( taui, isrt_ , tau_ , ityp_ , nat_ , nsp_ )
         ! stau0 contains "scaled" atomic positions (that is, in crystal axis)
         CALL r_to_s( taui, stau0, na, nsp, hinv )
         CALL ions_cofmass( taui, amass_ , na, nsp, cdmi )
      END IF
      !
      DEALLOCATE ( tau_, ityp_, isrt_ )
      
      CALL qexsd_copy_basis_set ( output_obj%basis_set, gamma_only_, ecutwfc_,&
           ecutrho_, nr1s, nr2s, nr3s, nr1, nr2, nr3, nr1b, nr2b, nr3b, &
           ngm_g, ngms_g, npw_g, b1, b2, b3 )

      CALL qexsd_copy_dft ( output_obj%dft, nsp, atm, dft_name, &
           nq1, nq2, nq3, ecutfock, exx_fraction, screening_parameter, &
           exxdiv_treatment, x_gamma_extrapolation, ecutvcut, &
           lda_plus_U, lda_plus_U_kind, U_projection, Hubbard_l, Hubbard_lmax,&
           Hubbard_U, Hubbard_dum(1,:), Hubbard_dum(2,:), Hubbard_dum(3,:), &
           Hubbard_dum, &
           vdw_corr,  llondon, ts_vdw, lxdm, inlc, vdw_table_name, scal6, &
           lon_rcut, vdw_isolated)
      !
      lsda_ = output_obj%magnetization%lsda
      IF ( lsda_ .AND. (nspin /= 2) ) CALL errore('cp_readfile','wrong spin',1)
      !
      nbnd_ = nupdwn(1)
      ALLOCATE( occ_(nbnd_, nspin), et_(nbnd_, nspin) )
      CALL qexsd_copy_band_structure( output_obj%band_structure, lsda_, nk_, &
           isk_, natomwfc, nbnd_, nelec_, wk_, occ_, ef, ef_up, ef_dw, et_ )
      ! FIXME: in the call, the same array is passed as both occ0 and occm!
      DO iss = 1, nspin
         ib = iupdwn(iss)
         nb = nupdwn(iss)
         occ0(ib:ib+nb-1) = occ_(1:nb,iss)
      END DO
      occm(:) = occ0(:)
      DEALLOCATE (occ_, et_)
      !
      CALL iotk_close_read (iunpun)
      !
      DO iss = 1, nspin
         CALL cp_read_wfc( ndr, tmp_dir, 1, 1, iss, nspin, c02, ' ' )
      END DO
      !
      RETURN
      !
    END SUBROUTINE cp_readfile
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_general_info (geninfo_obj, qexsd_fmt, qexsd_version) 
    !-------------------------------------------------------------------------------
    ! 
    USE qes_types_module,    ONLY: general_info_type
    !
    IMPLICIT NONE 
    !
    CHARACTER(LEN=*), INTENT(OUT) :: qexsd_fmt, qexsd_version
    TYPE (general_info_type ),INTENT(IN)  :: geninfo_obj   
    !
    qexsd_fmt = TRIM (geninfo_obj%xml_format%NAME)
    qexsd_version = TRIM ( geninfo_obj%xml_format%VERSION)
    !
  END SUBROUTINE qexsd_copy_general_info
  !
  !---------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_parallel_info (parinfo_obj, nproc_file, &
       nproc_pool_file, nproc_image_file, ntask_groups_file, &
       nproc_bgrp_file, nproc_ortho_file)
    !---------------------------------------------------------------------------    !
    USE qes_types_module, ONLY : parallel_info_type
    !
    IMPLICIT NONE 
    !
    TYPE ( parallel_info_type ),INTENT(IN)     :: parinfo_obj
    INTEGER, INTENT(OUT) :: nproc_file, nproc_pool_file, &
                            nproc_image_file, ntask_groups_file, &
                            nproc_bgrp_file, nproc_ortho_file
    ! 
    nproc_file = parinfo_obj%nprocs
    nproc_pool_file = nproc_file/parinfo_obj%npool
    nproc_image_file = nproc_file 
    ntask_groups_file = parinfo_obj%ntasks
    nproc_bgrp_file = nproc_image_file / parinfo_obj%npool / parinfo_obj%nbgrp 
    nproc_ortho_file = parinfo_obj%ndiag
    !
  END SUBROUTINE qexsd_copy_parallel_info  
  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_atomic_species (atomic_species, nsp, atm, psfile, amass)
    !---------------------------------------------------------------------------    !
    USE qes_types_module, ONLY : atomic_species_type
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_species
    INTEGER, INTENT(out) :: nsp
    CHARACTER(LEN=*), INTENT(out) :: atm(:), psfile(:)
    REAL(dp), INTENT(out) :: amass(:)
    !
    INTEGER :: isp
    !
    nsp = atomic_species%ntyp
    DO isp = 1, nsp 
       amass(isp) = 0.d0 
       IF (atomic_species%species(isp)%mass_ispresent) &
            amass(isp) = atomic_species%species(isp)%mass
       atm(isp) = TRIM ( atomic_species%species(isp)%name )
       psfile(isp) = TRIM ( atomic_species%species(isp)%pseudo_file) 
    END DO
    !
  END SUBROUTINE qexsd_copy_atomic_species

  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_atomic_structure (atomic_structure, nsp, atm, &
       nat, tau, ityp, alat, a1, a2, a3, ibrav )
  !--------------------------------------------------------------------------
    
    USE qes_types_module, ONLY : atomic_structure_type
    USE constants,        ONLY : pi
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    INTEGER, INTENT(in) :: nsp 
    CHARACTER(LEN = 3), INTENT(in) :: atm(:)
    !
    INTEGER, INTENT(out)  :: nat, ibrav, ityp(:)
    REAL(dp), INTENT(out) :: alat, a1(:), a2(:), a3(:), tau(:,:)
    !
    CHARACTER(LEN=3), ALLOCATABLE :: symbols(:)
    INTEGER :: iat, idx, isp
    !
    nat = atomic_structure%nat 
    alat = atomic_structure%alat 
    IF ( atomic_structure%bravais_index_ispresent ) THEN 
       ibrav = atomic_structure%bravais_index 
    ELSE 
       ibrav = 0
    END IF
    ALLOCATE ( symbols(nat) )
    loop_on_atoms:DO iat = 1, nat
       idx = atomic_structure%atomic_positions%atom(iat)%index
       tau(:,idx) = atomic_structure%atomic_positions%atom(iat)%atom 
       symbols(idx)  = TRIM ( atomic_structure%atomic_positions%atom(idx)%name)
       loop_on_species:DO isp = 1, nsp
          IF ( TRIM(symbols(idx)) == TRIM (atm(isp))) THEN 
             ityp(iat) = isp 
             exit loop_on_species
          END IF 
       END  DO loop_on_species
    END DO loop_on_atoms
    DEALLOCATE (symbols)
    IF ( atomic_structure%alat_ispresent ) alat = atomic_structure%alat 
    a1(:) = atomic_structure%cell%a1
    a2(:) = atomic_structure%cell%a2
    a3(:) = atomic_structure%cell%a3

  END SUBROUTINE qexsd_copy_atomic_structure

  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_basis_set ( basis_set, gamma_only, ecutwfc, ecutrho, &
       nr1s, nr2s, nr3s, nr1, nr2, nr3, nr1b, nr2b, nr3b, &
       ngm_g, ngms_g, npw_g, b1, b2, b3 )
  !--------------------------------------------------------------------------   
    !
    USE qes_types_module, ONLY : basis_set_type
    !
    IMPLICIT NONE 
    TYPE ( basis_set_type ),INTENT(IN)         :: basis_set
    LOGICAL, INTENT(out)  :: gamma_only
    INTEGER, INTENT(out)  :: ngm_g, ngms_g, npw_g
    INTEGER, INTENT(out)  :: nr1s, nr2s, nr3s, nr1, nr2, nr3
    INTEGER, INTENT(inout)  :: nr1b, nr2b, nr3b
    REAL(dp), INTENT(out) :: ecutwfc, ecutrho, b1(:), b2(:), b3(:)
    ! 
    ecutwfc = basis_set%ecutwfc
    ecutrho = basis_set%ecutrho
    gamma_only= basis_set%gamma_only
    nr1 = basis_set%fft_grid%nr1
    nr2 = basis_set%fft_grid%nr2          
    nr3 = basis_set%fft_grid%nr3
    nr1s= basis_set%fft_smooth%nr1
    nr2s= basis_set%fft_smooth%nr2
    nr3s= basis_set%fft_smooth%nr3
    IF ( basis_set%fft_box_ispresent ) THEN
       nr1b = basis_set%fft_box%nr1
       nr2b = basis_set%fft_box%nr2
       nr3b = basis_set%fft_box%nr3
    END IF
    ngm_g     = basis_set%ngm
    ngms_g    = basis_set%ngms
    npw_g     = basis_set%npwx
    !
    b1 =  basis_set%reciprocal_lattice%b1
    b2 =  basis_set%reciprocal_lattice%b2
    b3 =  basis_set%reciprocal_lattice%b3
    !
  END SUBROUTINE qexsd_copy_basis_set
  !
  !-----------------------------------------------------------------------
  SUBROUTINE qexsd_copy_dft ( dft_obj, nsp, atm, &
       dft_name, nq1, nq2, nq3, ecutfock, exx_fraction, screening_parameter, &
       exxdiv_treatment, x_gamma_extrapolation, ecutvcut, &
       lda_plus_U, lda_plus_U_kind, U_projection, Hubbard_l, Hubbard_lmax, &
       Hubbard_U, Hubbard_J0, Hubbard_alpha, Hubbard_beta, Hubbard_J, &
       vdw_corr,  llondon, ts_vdw, lxdm, inlc, vdw_table_name, scal6, &
       lon_rcut, vdw_isolated)
    !-------------------------------------------------------------------
    ! 
    USE qes_types_module, ONLY : dft_type
    !
    IMPLICIT NONE 
    TYPE ( dft_type ),INTENT(in) :: dft_obj
    INTEGER, INTENT(in)          :: nsp 
    CHARACTER(LEN=*), INTENT(in) :: atm(nsp)
    ! 
    CHARACTER(LEN=*), INTENT(out) :: dft_name
    ! Variables that may or may not be present should be intent(inout)
    ! so that they do not forget their default value (if any)
    CHARACTER(LEN=*), INTENT(inout) :: exxdiv_treatment
    REAL(dp), INTENT(inout) :: ecutfock, exx_fraction, screening_parameter, &
         ecutvcut
    INTEGER, INTENT(inout) :: nq1, nq2, nq3
    LOGICAL, INTENT(inout) :: x_gamma_extrapolation
    !
    LOGICAL, INTENT(out) :: lda_plus_U
    INTEGER, INTENT(inout) :: lda_plus_U_kind, Hubbard_lmax
    CHARACTER(LEN=*), INTENT(inout) :: U_projection
    INTEGER, INTENT(inout) :: Hubbard_l(:)
    REAL(dp), INTENT(inout) :: Hubbard_U(:), Hubbard_J0(:), Hubbard_J(:,:), &
         Hubbard_alpha(:), Hubbard_beta(:)
    !
    CHARACTER(LEN=256), INTENT(out) :: vdw_corr
    CHARACTER(LEN=256), INTENT(inout) :: vdw_table_name
    LOGICAL, INTENT(out) :: llondon, ts_vdw, lxdm
    INTEGER, INTENT(inout):: inlc
    REAL(dp), INTENT(inout) :: scal6, lon_rcut
    LOGICAL, INTENT(inout) :: vdw_isolated
    !
    CHARACTER(LEN=256 ) :: label
    CHARACTER(LEN=3 )   :: symbol
    INTEGER :: ihub, isp
    !
    dft_name = TRIM(dft_obj%functional)
    IF ( dft_obj%hybrid_ispresent ) THEN
       nq1 = dft_obj%hybrid%qpoint_grid%nqx1
       nq2 = dft_obj%hybrid%qpoint_grid%nqx2
       nq3 = dft_obj%hybrid%qpoint_grid%nqx3
       ecutfock = dft_obj%hybrid%ecutfock
       exx_fraction = dft_obj%hybrid%exx_fraction
       screening_parameter = dft_obj%hybrid%screening_parameter
       exxdiv_treatment = dft_obj%hybrid%exxdiv_treatment
       x_gamma_extrapolation = dft_obj%hybrid%x_gamma_extrapolation
       ecutvcut = dft_obj%hybrid%ecutvcut
    END IF
    !
    lda_plus_u = dft_obj%dftU_ispresent 
    IF ( lda_plus_u ) THEN 
       lda_plus_u_kind = dft_obj%dftU%lda_plus_u_kind
       U_projection = TRIM ( dft_obj%dftU%U_projection_type )
       Hubbard_l =-1 
       IF ( dft_obj%dftU%Hubbard_U_ispresent) THEN 
          loop_on_hubbardU:DO ihub =1, dft_obj%dftU%ndim_Hubbard_U
             symbol = TRIM(dft_obj%dftU%Hubbard_U(ihub)%specie)
             label  = TRIM(dft_obj%dftU%Hubbard_U(ihub)%label ) 
             loop_on_speciesU:DO isp = 1, nsp
                IF ( TRIM(symbol) == TRIM ( atm(isp) ) ) THEN 
                     Hubbard_U(isp) = dft_obj%dftU%Hubbard_U(ihub)%HubbardCommon
                     SELECT CASE ( TRIM (label))
                     CASE ( '1s', '2s', '3s', '4s', '5s', '6s', '7s' ) 
                        Hubbard_l(isp) = 0 
                     CASE ( '2p', '3p', '4p', '5p', '6p' ) 
                        Hubbard_l(isp) = 1 
                     CASE ( '3d', '4d', '5d' ) 
                        Hubbard_l( isp ) = 2 
                     CASE ( '4f', '5f' )  
                        Hubbard_l(isp ) = 3
                     CASE  default 
                        IF (Hubbard_U(isp)/=0) &
                             CALL errore ("pw_readschema:", "unrecognized label for Hubbard "//label, 1 ) 
                     END SELECT
                     EXIT loop_on_speciesU
                  END IF 
                END DO loop_on_speciesU
            END DO loop_on_hubbardU
         END IF 
         
         IF ( dft_obj%dftU%Hubbard_J0_ispresent ) THEN 
            loop_on_hubbardj0:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J0
               symbol = TRIM(dft_obj%dftU%Hubbard_J0(ihub)%specie)
               loop_on_speciesj0:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN
                     Hubbard_J0(isp) = dft_obj%dftU%Hubbard_J0(ihub)%HubbardCommon
                     EXIT loop_on_speciesj0
                  END IF
               END DO loop_on_speciesj0
            END DO loop_on_hubbardj0
         END IF
         IF ( dft_obj%dftU%Hubbard_alpha_ispresent) THEN 
            loop_on_hubbardAlpha:DO ihub =1, dft_obj%dftU%ndim_Hubbard_alpha
               symbol = TRIM(dft_obj%dftU%Hubbard_alpha(ihub)%specie)
               loop_on_speciesAlpha:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN 
                     Hubbard_alpha(isp) = dft_obj%dftU%Hubbard_alpha(ihub)%HubbardCommon
                     EXIT loop_on_speciesAlpha
                  END IF
               END DO loop_on_speciesAlpha
            END DO loop_on_hubbardAlpha
         END IF
         IF ( dft_obj%dftU%Hubbard_beta_ispresent) THEN 
            loop_on_hubbardBeta:DO ihub =1, dft_obj%dftU%ndim_Hubbard_beta
               symbol = TRIM(dft_obj%dftU%Hubbard_beta(ihub)%specie)
               loop_on_speciesBeta:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN 
                     Hubbard_beta(isp) = dft_obj%dftU%Hubbard_beta(ihub)%HubbardCommon
                     EXIT loop_on_speciesBeta
                  END IF
               END DO loop_on_speciesBeta
            END DO loop_on_hubbardBeta
         END IF
         IF ( dft_obj%dftU%Hubbard_J_ispresent) THEN 
            loop_on_hubbardJ:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J
               symbol = TRIM(dft_obj%dftU%Hubbard_J(ihub)%specie)
               loop_on_speciesJ:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN 
                     Hubbard_J(:,isp) = dft_obj%dftU%Hubbard_J(ihub)%HubbardJ
                     EXIT loop_on_speciesJ
                  END IF
               END DO loop_on_speciesJ
            END DO loop_on_hubbardJ
         END IF
         Hubbard_lmax = MAXVAL( Hubbard_l(1:nsp) )
      END IF

      IF ( dft_obj%vdW_ispresent ) THEN 
         vdw_corr = TRIM( dft_obj%vdW%vdw_corr ) 
      ELSE
         vdw_corr = ''
      END IF
      SELECT CASE( TRIM( dft_obj%vdW%vdw_corr ) )
         !
      CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
         !
         llondon= .TRUE.
         ts_vdw= .FALSE.
         lxdm   = .FALSE.
         !
      CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
         !
         llondon= .FALSE.
         ts_vdw= .TRUE.
         lxdm   = .FALSE.
         !
      CASE( 'XDM', 'xdm' )
         !
         llondon= .FALSE.
         ts_vdw= .FALSE.
         lxdm   = .TRUE.
         !
      CASE DEFAULT
         !
         llondon= .FALSE.
         ts_vdw = .FALSE.
         lxdm   = .FALSE.
         !
      END SELECT
      IF ( dft_obj%vdW_ispresent ) THEN 
         SELECT CASE ( TRIM (dft_obj%vdW%non_local_term))
         CASE ('vdw1')  
            inlc = 1
         CASE ('vdw2') 
            inlc = 2
         CASE ('vv10' ) 
            inlc = 3 
         CASE ( 'vdW-DF-x') 
            inlc = 4
         CASE ( 'vdW-DF-y')
            inlc = 5
         CASE ( 'vdW-DF-z')
            inlc = 6
         CASE default 
            inlc = 0 
         END SELECT
         IF (inlc == 0 ) THEN 
            vdw_table_name = ' '
         ELSE IF ( inlc == 3 ) THEN 
            vdw_table_name = 'rVV10_kernel_table'
         ELSE
            vdw_table_name = 'vdW_kernel_table'
         END IF
         IF (dft_obj%vdW%london_s6_ispresent ) THEN 
            scal6 = dft_obj%vdW%london_s6
         END IF
         IF ( dft_obj%vdW%london_rcut_ispresent ) THEN 
            lon_rcut = dft_obj%vdW%london_rcut
         END IF
         IF (dft_obj%vdW%ts_vdW_isolated_ispresent ) THEN 
            vdW_isolated = dft_obj%vdW%ts_vdW_isolated
         END IF
      END IF
   
    END SUBROUTINE qexsd_copy_dft
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_band_structure( band_struct_obj, lsda, nkstot, &
         isk, natomwfc, nbnd, nelec, wk, wg, ef, ef_up, ef_dw, et )
      !------------------------------------------------------------------------
      !
      USE qes_types_module, ONLY : band_structure_type
      !
      IMPLICIT NONE
      TYPE ( band_structure_type)         :: band_struct_obj
      LOGICAL, INTENT(out) :: lsda
      INTEGER, INTENT(out) :: nkstot, natomwfc, nbnd, isk(:)
      REAL(dp), INTENT(out):: nelec, wk(:), wg(:,:)
      REAL(dp), INTENT(out):: ef, ef_up, ef_dw, et(:,:)
      !
      INTEGER :: ik, nbnd_, nbnd_up_, nbnd_dw_
      ! 
      lsda = band_struct_obj%lsda
      nkstot = band_struct_obj%nks 
      IF ( lsda) THEN 
         nkstot = nkstot * 2 
         isk(1:nkstot/2) = 1
         isk(nkstot/2+1:nkstot) = 2 
      ELSE 
         isk(1:nkstot)   = 1 
      END IF
      ! 
      nelec = band_struct_obj%nelec
      nbnd  = band_struct_obj%nbnd 
      natomwfc = band_struct_obj%num_of_atomic_wfc
      IF ( band_struct_obj%fermi_energy_ispresent) THEN 
         ef = band_struct_obj%fermi_energy
         ef_up = 0.d0
         ef_dw = 0.d0
      ELSE IF ( band_struct_obj%two_fermi_energies_ispresent ) THEN 
         ef = 0.d0 
         ef_up = band_struct_obj%two_fermi_energies(1)
         ef_dw = band_struct_obj%two_fermi_energies(2)
      ELSE 
         ef = 0.d0
         ef_up = 0.d0
         ef_dw = 0.d0
      END IF
      DO ik =1, band_struct_obj%ndim_ks_energies
         IF ( band_struct_obj%lsda) THEN
            IF ( band_struct_obj%nbnd_up_ispresent .AND. band_struct_obj%nbnd_dw_ispresent) THEN
               nbnd_up_ = band_struct_obj%nbnd_up
               nbnd_dw_ = band_struct_obj%nbnd_dw 
            ELSE IF ( band_struct_obj%nbnd_up_ispresent ) THEN 
               nbnd_up_ = band_struct_obj%nbnd_up
               nbnd_dw_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues - nbnd_up_
            ELSE IF ( band_struct_obj%nbnd_dw_ispresent ) THEN 
               nbnd_dw_ = band_struct_obj%nbnd_dw
               nbnd_up_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues - nbnd_dw_ 
            ELSE 
               nbnd_up_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues/2  
               nbnd_dw_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues/2
            END IF
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            wk( ik + band_struct_obj%ndim_ks_energies ) = wk(ik) 
            et(1:nbnd_up_,ik) = band_struct_obj%ks_energies(ik)%eigenvalues(1:nbnd_up_)
            et(1:nbnd_dw_,ik+band_struct_obj%ndim_ks_energies) =  &
                 band_struct_obj%ks_energies(ik)%eigenvalues(nbnd_up_+1:nbnd_up_+nbnd_dw_)
            wg(1:nbnd_up_,ik) = band_struct_obj%ks_energies(ik)%occupations(1:nbnd_up_)*wk(ik)
            wg(1:nbnd_dw_,ik+band_struct_obj%ndim_ks_energies) =  &
                 band_struct_obj%ks_energies(ik)%occupations(nbnd_up_+1:nbnd_up_+nbnd_dw_)*wk(ik)
         ELSE 
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            nbnd_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues
            et (1:nbnd_,ik) = band_struct_obj%ks_energies(ik)%eigenvalues(1:nbnd_)
            wg (1:nbnd_,ik) = band_struct_obj%ks_energies(ik)%occupations(1:nbnd_)*wk(ik)
         END IF
      END DO
    END SUBROUTINE qexsd_copy_band_structure
  !------------------------------------------------------------------------
  SUBROUTINE cp_writecp( iunpun, nfi, simtime, &
       ekin, eht, esr, eself, epseu, enl, exc, vave, enthal, &
       acc, stau0, svel0, taui, cdmi, force, nhpcl, nhpdim, &
       xnhp0, vnhp, ekincm, xnhe0, vnhe, ht, htvel, gvel, xnhh0, vnhh,      &
       staum, svelm, xnhpm, xnhem, htm, xnhhm) !
    !------------------------------------------------------------------------
    ! ... Cell related variables, CP-specific
    !
    USE iotk_module
    USE ions_base, ONLY: nat
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: iunpun
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
    !
    CHARACTER(iotk_attlenx)  :: attr
    !
    IF ( ionode ) THEN
!-------------------------------------------------------------------------------
! ... STATUS
!-------------------------------------------------------------------------------
       !
       CALL iotk_write_begin( iunpun, "STATUS" )
       !
       CALL iotk_write_attr( attr, "ITERATION", nfi, FIRST = .TRUE. )
       CALL iotk_write_empty( iunpun, "STEP", attr )
       !
       CALL iotk_write_attr( attr, "UNITS", "pico-seconds", FIRST = .TRUE. )
       CALL iotk_write_dat( iunpun, "TIME", simtime, ATTR = attr )
       !
       CALL iotk_write_dat( iunpun, "TITLE", 'temorary title' )
       !
       CALL iotk_write_attr( attr, "UNITS", 'Hartree', FIRST = .TRUE. )
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
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE cp_writecp
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_read_wfc( ndr, tmp_dir, ik, nk, iss, nspin, c2, tag )
    !------------------------------------------------------------------------
    !
    ! Wrapper for old cp_read_wfc
    !
    USE io_global,          ONLY : ionode
    USE io_files,           ONLY : prefix, iunpun
    USE mp_global,          ONLY : root_pool, intra_pool_comm
    USE electrons_base,     ONLY : iupdwn, nupdwn
    USE gvecw,              ONLY : ngw, ngw_g
    USE gvect,              ONLY : ig_l2g
    !
    IMPLICIT NONE
    !
    INTEGER,               INTENT(IN)  :: ndr
    CHARACTER(LEN=*),      INTENT(IN)  :: tmp_dir
    INTEGER,               INTENT(IN)  :: ik, iss, nk, nspin
    CHARACTER,             INTENT(IN)  :: tag
    COMPLEX(DP),           INTENT(OUT) :: c2(:,:)
    !
    INTEGER            :: ib, nb, nbnd, is_, ns_
    CHARACTER(LEN=320) :: filename
    REAL(DP)           :: scalef
    !
print *, 'entering cp_read_wfc_new'
    IF ( tag == 'm' ) THEN
       WRITE(filename,'(A,A,"_",I2,".save/wfcm",I1)') &
            TRIM(tmp_dir), TRIM(prefix), ndr, iss
    ELSE
       WRITE(filename,'(A,A,"_",I2,".save/wfc",I1)') &
            TRIM(tmp_dir), TRIM(prefix), ndr, iss
    END IF
    ib = iupdwn(iss)
    nb = nupdwn(iss)
    ! next two lines workaround for bogus complaint due to intent(in)
    is_= iss
    ns_= nspin
    CALL read_wfc( iunpun, is_, nk, is_, ns_, &
         c2(:,ib:ib+nb-1), ngw_g, nbnd, ig_l2g, ngw,  &
         filename, scalef, ionode, root_pool, intra_pool_comm )
    !
  END SUBROUTINE cp_read_wfc
  !
  !------------------------------------------------------------------------
  SUBROUTINE cp_readcp ( iunpun, nat, nfi, simtime, acc, stau0, svel0, taui,&
       cdmi, force, nhpcl, nhpdim, xnhp0, vnhp, ekincm, xnhe0, vnhe, ht, &
       htvel, gvel, xnhh0, vnhh, staum, svelm, xnhpm, xnhem, htm, xnhhm, &
       ierr )
    !
    !------------------------------------------------------------------------
    ! ... Cell related variables, CP-specific
    ! ... ierr = -2: nothing found
    ! ... ierr = -1: MD status found, no info on timesteps
    ! ... ierr =  0: MD status and timestep info read
    ! ... ierr =  1: error reading MD status
    ! ... ierr =  2: error reading timestep info
    !
    USE iotk_module
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: iunpun
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
    INTEGER,  INTENT(out) :: ierr
    !
    LOGICAL :: found
    INTEGER :: nt_, nhpcl_, nhpdim_
    CHARACTER(iotk_attlenx)  :: attr
    !
    ! ... read MD status
    !
    ierr = -2
    CALL iotk_scan_begin( iunpun, "STATUS", attr, FOUND = found )
    IF ( .NOT.found ) RETURN
    !
    ierr = 1
    CALL iotk_scan_empty( iunpun, "STEP", ATTR = attr, FOUND = found )
    IF ( .NOT.found ) RETURN
    !
    CALL iotk_scan_attr( attr, "ITERATION", nfi, FOUND = found )
    IF ( .NOT.found ) RETURN
    !
    CALL iotk_scan_dat( iunpun, "TIME", simtime, ATTR = attr, FOUND = found  )
    IF ( .NOT.found ) RETURN
    !
    CALL iotk_scan_end( iunpun, "STATUS", IERR=ierr )
    IF ( ierr /= 0 ) RETURN
    !
    ! ... read MD timesteps variables
    !
    CALL iotk_scan_begin( iunpun, "TIMESTEPS", attr, FOUND = found )
    ! 
    IF ( found ) THEN
       !
       ierr = 0
       !
       CALL iotk_scan_attr( attr, "nt", nt_ )
       !
       IF ( nt_ > 0 ) THEN
          !
          CALL iotk_scan_begin( iunpun, "STEP0" )
          !
          CALL iotk_scan_dat ( iunpun, "ACCUMULATORS", acc )
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
          IF ( nhpcl_ == nhpcl .AND. nhpdim_ == nhpdim ) THEN
             CALL iotk_scan_dat( iunpun, "xnhp", xnhp0(1:nhpcl*nhpdim) )
             CALL iotk_scan_dat( iunpun, "vnhp", vnhp(1:nhpcl*nhpdim) )
          ELSE
             xnhp0(1:nhpcl*nhpdim) = 0.D0
             vnhp(1:nhpcl*nhpdim)  = 0.D0
          END IF
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
          ierr = 2
          RETURN
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
             CALL iotk_scan_dat( iunpun, "xnhp",  xnhpm(1:nhpcl*nhpdim) )
          ELSE
             xnhpm(1:nhpcl*nhpdim) = 0.D0
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
  SUBROUTINE cp_read_cell( ndr, tmp_dir, ascii, ht, &
                           htm, htvel, gvel, xnhh0, xnhhm, vnhh )
    !------------------------------------------------------------------------
    !
    USE parameters,  ONLY : ntypx
    USE ions_base,   ONLY : nat
    USE qexsd_reader_module,  ONLY : qexsd_get_output
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
    INTEGER          :: nat_
    INTEGER          :: nsp_
    INTEGER          :: ityp_(nat) 
    REAL(DP)         :: alat_
    REAL(DP)         :: a1_(3), a2_(3), a3_(3)
    REAL(DP)         :: b1_(3), b2_(3), b3_(3)
    REAL(DP)         :: tau_(3,nat) 
    CHARACTER(LEN=3) :: atm_(ntypx)
    CHARACTER(iotk_attlenx)  :: attr
    TYPE(output_type) :: output_obj
    !
    ! ... look for an empty unit
    !
    CALL iotk_free_unit( iunpun, ierr )
    CALL errore( 'cp_read_cell', 'no free units ', ierr )
    !
    CALL qexsd_init_schema( iunpun )
    !
    WRITE(dirname,'(A,A,"_",I2,".save/")') TRIM(tmp_dir), TRIM(prefix), ndr
    filename = TRIM( dirname ) // TRIM( xmlpun_schema )
    INQUIRE ( file=filename, exist=found )
    IF (.NOT. found ) &
         CALL errore ('cp_read_cell', 'xml data file not found', 1)
    !
    CALL iotk_open_read( iunpun, TRIM(filename) )
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
       CALL qexsd_get_output ( iunpun, output_obj, found ) 
       CALL qexsd_copy_atomic_structure (output_obj%atomic_structure, nsp_, &
            atm_, nat_, tau_, ityp_, alat_, a1_, a2_, a3_, ibrav_ )
       IF ( nat_ /= nat ) CALL errore ('cp_readfile', 'wrong nat read', 1)
       CALL qes_reset_output(output_obj)
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
100 CALL errore( 'cp_read_cell ', attr, ierr )
    !
  END SUBROUTINE cp_read_cell
#endif
  !
END MODULE cp_restart_new
