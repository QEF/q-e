!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE init_run()
  !----------------------------------------------------------------------------
  !
  ! ... this routine initialise the cp code and allocates (calling the
  ! ... appropriate routines) the memory
  !
  USE kinds,                    ONLY : DP
  USE control_flags,            ONLY : nbeg, nomore, lwf, iverbosity, iprint, &
                                       ndr, ndw, tfor, tprnfor, tpre, ts_vdw, &
                                       force_pairing
  USE cp_electronic_mass,       ONLY : emass, emass_cutoff
  USE ions_base,                ONLY : na, nax, nat, nsp, iforce, amass, cdms
  USE ions_positions,           ONLY : tau0, taum, taup, taus, tausm, tausp, &
                                       vels, velsm, velsp, fion, fionm
  USE gvecw,                    ONLY : ngw, ngw_g, g2kin, g2kin_init
  USE smallbox_gvec,            ONLY : ngb
  USE gvecs,                    ONLY : ngms
  USE gvect,                    ONLY : ngm, gstart, gg
  USE fft_base,                 ONLY : dfftp, dffts
  USE electrons_base,           ONLY : nspin, nbsp, nbspx, nupdwn, f
  USE uspp,                     ONLY : nkb, vkb, deeq, becsum,nkbus
  USE core,                     ONLY : rhoc
  USE wavefunctions_module,     ONLY : c0_bgrp, cm_bgrp, phi_bgrp
  USE ensemble_dft,             ONLY : tens, z0t
  USE cg_module,                ONLY : tcg
  USE electrons_base,           ONLY : nudx, nbnd
  USE efield_module,            ONLY : tefield, tefield2
  USE uspp_param,               ONLY : nhm
  USE ions_nose,                ONLY : xnhp0, xnhpm, vnhp, nhpcl, nhpdim
  USE cell_base,                ONLY : h, hold, hnew, velh, tpiba2, ibrav, &
                                       alat, celldm, at, bg
  USE cp_main_variables,        ONLY : lambda, lambdam, lambdap, ema0bg, &
                                       sfac, eigr, taub, &
                                       irb, eigrb, rhog, rhos, rhor,     &
                                       acc, acc_this_run, wfill, &
                                       edft, nfi, vpot, ht0, htm, iprint_stdout
  USE cp_main_variables,        ONLY : allocate_mainvar, descla
  USE energies,                 ONLY : eself, enl, ekin, etot, enthal, ekincm
  USE dener,                    ONLY : detot
  USE time_step,                ONLY : dt2, delt, tps
  USE electrons_nose,           ONLY : xnhe0, xnhem, vnhe
  USE electrons_base,           ONLY : nbspx_bgrp
  USE cell_nose,                ONLY : xnhh0, xnhhm, vnhh
  USE funct,                    ONLY : dft_is_meta, dft_is_hybrid
  USE metagga,                  ONLY : crosstaus, dkedtaus, gradwfc
  !
  USE efcalc,                   ONLY : clear_nbeg
  USE local_pseudo,             ONLY : allocate_local_pseudo
  USE cp_electronic_mass,       ONLY : emass_precond
  USE wannier_subroutines,      ONLY : wannier_startup
  USE cp_interfaces,            ONLY : readfile
  USE ions_base,                ONLY : ions_cofmass
  USE ensemble_dft,             ONLY : id_matrix_init, allocate_ensemble_dft, h_matrix_init
  USE efield_module,            ONLY : allocate_efield, allocate_efield2
  USE cg_module,                ONLY : allocate_cg
  USE wannier_module,           ONLY : allocate_wannier  
  USE io_files,                 ONLY : tmp_dir, prefix
  USE io_global,                ONLY : ionode, stdout
  USE printout_base,            ONLY : printout_base_init
  USE wave_types,               ONLY : wave_descriptor_info
  USE xml_io_base,              ONLY : restart_dir, create_directory
  USE orthogonalize_base,       ONLY : mesure_diag_perf, mesure_mmul_perf
  USE ions_base,                ONLY : ions_reference_positions, cdmi
  USE mp_bands,                 ONLY : nbgrp
  USE mp,                       ONLY : mp_barrier
  USE wrappers
  USE ldaU_cp
  USE control_flags,            ONLY : lwfpbe0nscf         ! exx_wf related 
  USE wavefunctions_module,     ONLY : cv0                 ! exx_wf related
  USE wannier_base,             ONLY : vnbsp               ! exx_wf related
  USE cp_restart,               ONLY : cp_read_wfc_Kong    ! exx_wf related
  USE input_parameters,         ONLY : ref_cell
  USE cell_base,                ONLY : ref_tpiba2, init_tpiba2
  USE tsvdw_module,             ONLY : tsvdw_initialize
  USE exx_module,               ONLY : exx_initialize
  !
  IMPLICIT NONE
  !
  INTEGER            :: i
  CHARACTER(LEN=256) :: dirname
  REAL(DP)           :: a1(3), a2(3), a3(3)
  LOGICAL            :: ftest
  !
  !
  CALL start_clock( 'initialize' )
  !
  ! ... initialize directories
  !
  IF( nbeg < 0 ) THEN
     CALL create_directory( tmp_dir )
  END IF 
  !
  CALL plugin_initialization()

  IF( nbgrp > 1 .AND. force_pairing ) &
     CALL errore( ' init_run ', ' force_pairing with parallelization over bands not implemented yet ', 1 )
  !
  CALL printout_base_init( tmp_dir, prefix )
  !
  dirname = restart_dir( tmp_dir, ndw )
  !
  ! ... Create main restart directory
  !
  CALL create_directory( dirname )
  !
  ! ... initialize g-vectors, fft grids 
  ! ... The number of g-vectors are based on the input celldm!
  !
  CALL init_dimensions()
  !
  ! ... initialization of plugin variables and arrays
  !
  CALL plugin_init_base() 
  ! 
  ! ... initialize atomic positions and cell
  !
  CALL init_geometry()
  !
  ! ... mesure performances of parallel routines
  !
  CALL mesure_mmul_perf( nudx )
  !
  CALL mesure_diag_perf( nudx )
  !
  IF ( lwf ) CALL clear_nbeg( nbeg )
  !
  !=======================================================================
  !     allocate and initialize local and nonlocal potentials
  !=======================================================================
  !
  CALL allocate_local_pseudo( ngms, nsp )
  !
  CALL nlinit()
  !
  !=======================================================================
  !     allocation of all arrays not already allocated in init and nlinit
  !=======================================================================
  !
  CALL allocate_mainvar( ngw, ngw_g, ngb, ngms, ngm, dfftp%nr1,dfftp%nr2,dfftp%nr3, dfftp%nr1x, &
                         dfftp%nr2x, dfftp%npl, dfftp%nnr, dffts%nnr, nat, nax, nsp,   &
                         nspin, nbsp, nbspx, nupdwn, nkb, gstart, nudx, &
                         tpre, nbspx_bgrp )
  ! 
  !=======================================================================
  !     Initialization of the TS-vdW code (RAD)
  !=======================================================================
  !
  IF (ts_vdw) CALL tsvdw_initialize()
  !
  !=======================================================================
  !     Initialization of the exact exchange code (exx_module)
  !=======================================================================
  !exx_wf related
  IF ( dft_is_hybrid() .AND. lwf ) THEN
    !
    CALL exx_initialize()
    !
  END IF
  !
  !  initialize wave functions descriptors and allocate wf
  !
  IF(lwfpbe0nscf) ALLOCATE(cv0( ngw, vnbsp ) )   ! Lingzhu Kong
  ALLOCATE( c0_bgrp( ngw, nbspx ) )
  ALLOCATE( cm_bgrp( ngw, nbspx ) )
  ALLOCATE( phi_bgrp( ngw, nbspx ) )
  !
  IF ( iverbosity > 1 ) THEN
     !
     CALL wave_descriptor_info( wfill, 'wfill', stdout )
     !
  END IF
  !
  ! Depending on the verbosity set the frequency of
  ! verbose information to stdout
  !
  IF( iverbosity < 0 ) iprint_stdout = 100 * iprint
  IF( iverbosity ==0 .OR. iverbosity == 1 ) iprint_stdout = 10 * iprint
  IF( iverbosity > 1 ) iprint_stdout = iprint
  !
  acc          = 0.D0
  acc_this_run = 0.D0
  !
  edft%ent  = 0.D0
  edft%esr  = 0.D0
  edft%evdw = 0.D0
  edft%ekin = 0.D0
  edft%enl  = 0.D0
  edft%etot = 0.D0
  !
  ALLOCATE( becsum(  nhm*(nhm+1)/2, nat, nspin ) )
  ALLOCATE( deeq( nhm, nhm, nat, nspin ) )
  !
  ALLOCATE( vkb( ngw, nkb ) )
  !
  IF ( dft_is_meta() .AND. tens ) &
     CALL errore( ' init_run ', 'ensemble_dft not implemented for metaGGA', 1 )
  !
  IF ( dft_is_meta() .AND. nbgrp > 1 ) &
     CALL errore( ' init_run ', 'band parallelization not implemented for metaGGA', 1 )
  !
  IF ( dft_is_meta() .AND. tpre ) THEN
     !
     ALLOCATE( crosstaus( dffts%nnr, 6, nspin ) )
     ALLOCATE( dkedtaus(  dffts%nnr, 3, 3, nspin ) )
     ALLOCATE( gradwfc(   dffts%nnr, 3 ) )
     !
  END IF
  !
  IF ( lwf ) THEN
     IF( nbgrp > 1 ) &
        CALL errore( ' init_run ', ' wannier with band parallelization not implemented ', 1 )
     CALL allocate_wannier( nbsp, dffts%nnr, nspin, ngm )
  END IF
  !
  IF ( tens .OR. tcg ) THEN
     IF( nbgrp > 1 ) &
        CALL errore( ' init_run ', ' ensemble_dft with band parallelization not implemented ', 1 )
     CALL allocate_ensemble_dft( nkb, nbsp, ngw, nudx, nspin, nbspx, &
                                 dffts%nnr, nat, descla )
  END IF
  !
  IF ( tcg ) THEN 
     CALL allocate_cg( ngw, nbspx,nkbus )
  END IF
  !
  IF ( tefield ) THEN
     IF( nbgrp > 1 ) &
        CALL errore( ' init_run ', ' efield with band paralleliztion not implemented ', 1 )
     CALL allocate_efield( ngw, ngw_g, nbspx, nhm, nax, nsp )
  END IF
  IF ( tefield2 ) THEN
     IF( nbgrp > 1 ) &
        CALL errore( ' init_run ', ' efield with band paralleliztion not implemented ', 1 )
     CALL allocate_efield2( ngw, nbspx, nhm, nax, nsp )
  END IF
  !
  IF ( ALLOCATED( deeq ) ) deeq(:,:,:,:) = 0.D0
  !
  IF ( ALLOCATED( lambda  ) ) lambda  = 0.D0
  IF ( ALLOCATED( lambdam ) ) lambdam = 0.D0
  !
  taum  = tau0
  taup  = 0.D0
  tausm = taus
  tausp = 0.D0
  vels  = 0.D0
  velsm = 0.D0
  velsp = 0.D0
  !
  hnew = h
  !
  IF(lwfpbe0nscf) cv0=( 0.D0, 0.D0 )    ! Lingzhu Kong
  cm_bgrp  = ( 0.D0, 0.D0 )
  c0_bgrp  = ( 0.D0, 0.D0 )
  phi_bgrp = ( 0.D0, 0.D0 )
  !
  IF ( tens ) then
     CALL id_matrix_init( descla, nspin )
     CALL h_matrix_init( descla, nspin )
  ENDIF
  !
  a1(:)=at(:,1)*alat; a2(:)=at(:,2)*alat; a3(:)=at(:,3)*alat 
  IF ( lwf ) CALL wannier_startup( ibrav, alat, a1, a2, a3, &
                                   bg(:,1), bg(:,2), bg(:,3) )
  !
  ! ... Calculate: ema0bg = ecutmass /  MAX( 1.0d0, (2pi/alat)^2 * |G|^2 )
  !
  IF ( ref_cell ) THEN
    WRITE( stdout,'(/,3X,"Reference cell parameters are used in electron mass preconditioning")' )
    WRITE( stdout,'(3X,"ref_tpiba2=",F14.8)' ) ref_tpiba2
    CALL g2kin_init( gg, ref_tpiba2 )
    CALL emass_precond( ema0bg, g2kin, ngw, ref_tpiba2, emass_cutoff ) 
    WRITE( stdout,'(3X,"current_tpiba2=",F14.8)' ) tpiba2
    CALL g2kin_init( gg, tpiba2 )
  ELSE
    WRITE( stdout,'(/,3X,"Cell parameters from input file are used in electron mass preconditioning")' )
    WRITE( stdout,'(3X,"init_tpiba2=",F14.8)' ) init_tpiba2
    CALL g2kin_init( gg, init_tpiba2 )
    CALL emass_precond( ema0bg, g2kin, ngw, init_tpiba2, emass_cutoff ) 
    !WRITE( stdout,'(3X,"current_tpiba2=",F14.8)' ) tpiba2 !BS : DEBUG
    CALL g2kin_init( gg, tpiba2 )
  END IF
  !
  CALL print_legend( )
  !
  CALL ldaU_init()
  !
  IF ( nbeg < 0 ) THEN
     !
     !======================================================================
     !     Initialize from scratch nbeg = -1
     !======================================================================
     !
     nfi = 0
     !
     CALL from_scratch( )
     !
  ELSE
     !
     !======================================================================
     !     nbeg = 0, nbeg = 1
     !======================================================================
     !
     !======================================================================
     ! Kong, read the valence orbitals
     IF(lwfpbe0nscf) THEN
        CALL cp_read_wfc_Kong( 36, tmp_dir, 1, 1, 1, 1, cv0, 'v' )
     ENDIF
     !======================================================================
     i = 1  
     CALL start_clock( 'init_readfile' )
     CALL readfile( i, h, hold, nfi, c0_bgrp, cm_bgrp, taus,   &
                    tausm, vels, velsm, acc, lambda, lambdam, xnhe0, xnhem, &
                    vnhe, xnhp0, xnhpm, vnhp,nhpcl,nhpdim,ekincm, xnhh0, xnhhm,&
                    vnhh, velh, fion, tps, z0t, f )
     !
     CALL from_restart( )
     !
     CALL stop_clock( 'init_readfile' )
     !
  END IF
  !
  !=======================================================================
  !     restart with new averages and nfi=0
  !=======================================================================
  !
  ! ... reset some variables if nbeg < 0 
  ! ... ( new simulation or step counter reset to 0 )
  !
  IF ( nbeg <= 0 ) THEN
     !
     acc = 0.D0
     nfi = 0
     !
  END IF
  !
  IF ( .NOT. tfor .AND. .NOT. tprnfor ) fion(:,:) = 0.D0
  !
  nomore = nomore + nfi
  !
  !  Set center of mass for scaled coordinates
  !
  CALL ions_cofmass( taus, amass, na, nsp, cdms )
  !
  IF ( nbeg <= 0 .OR. lwf ) THEN
     !
     CALL ions_reference_positions( tau0 ) ! BS: screws up msd calculation for lwf ...
     !
  END IF
  !
  CALL stop_clock( 'initialize' )
  !
  RETURN
  !
END SUBROUTINE init_run
