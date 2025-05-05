!
! Copyright (C) 2001-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE init_run()
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE klist,              ONLY : nkstot
  USE start_k,            ONLY : nks_start, nk1, nk2, nk3, k1, k2, k3
  USE symme,              ONLY : sym_rho_init
  USE wvfct,              ONLY : nbnd, et, wg, btype
  USE control_flags,      ONLY : lmd, gamma_only, smallmem, ts_vdw, mbd_vdw, &
                                 lforce => tprnfor, tstress, tqr, use_gpu
  USE gvect,              ONLY : g, gg, mill, gcutm, ig_l2g, ngm, ngm_g, &
                                 gshells, gstart ! to be communicated to the Solvers if gamma_only
  USE gvecs,              ONLY : gcutms, ngms
  USE ions_base,          ONLY : nat, nsp
  USE cell_base,          ONLY : alat, at, bg, set_h_ainv, omega
  USE cellmd,             ONLY : lmovecell
  USE dynamics_module,    ONLY : allocate_dyn_vars
  USE paw_variables,      ONLY : okpaw
  USE paw_init,           ONLY : paw_init_onecenter, allocate_paw_internals
#if defined(__MPI)
  USE paw_init,           ONLY : paw_post_init
#endif
  USE bp,                 ONLY : allocate_bp_efield, bp_global_map
  USE fft_base,           ONLY : dfftp, dffts
  USE recvec_subs,        ONLY : ggen, ggens
  USE wannier_new,        ONLY : use_wannier    
  USE dfunct,             ONLY : newd
  USE martyna_tuckerman,  ONLY : do_comp_mt
  USE esm,                ONLY : do_comp_esm, esm_init
  USE tsvdw_module,       ONLY : tsvdw_initialize
  USE libmbd_interface,   ONLY : init_mbd
  USE Coul_cut_2D,        ONLY : do_cutoff_2D, cutoff_fact 
  USE two_chem,           ONLY : init_twochem, twochem
  USE lsda_mod,           ONLY : nspin
  USE noncollin_module,   ONLY : domag, noncolin, lspinorb
  USE xc_lib,             ONLY : xclib_dft_is_libxc, xclib_init_libxc, &
                                 xclib_dft_is, xclib_set_finite_size_volume, &
                                 dft_has_finite_size_correction
  !
  USE rism_module,        ONLY : lrism, rism_alloc3d
  USE extffield,          ONLY : init_extffield
  USE control_flags,      ONLY : scissor
  USE sci_mod,            ONLY : allocate_scissor
  USE uspp_param,         ONLY : nhm
  USE uspp,               ONLY : allocate_uspp
  !
#if defined (__ENVIRON)
  USE plugin_flags,        ONLY : use_environ
  USE environ_base_module, ONLY : init_environ_base
#endif
  !
  IMPLICIT NONE
  INTEGER :: ierr
  !
#if defined (__ENVIRON)
  REAL(DP) :: at_scaled(3, 3)
  REAL(DP) :: gcutm_scaled
#endif
  !
  CALL start_clock( 'init_run' )
  !
  ! ... calculate limits of some indices, used in subsequent allocations
  !
  CALL pre_init()
  !
  ! ... determine the data structure for fft arrays
  !
  CALL data_structure( gamma_only )
  !
  ! ... print a summary and a memory estimate before starting allocating
  !
  CALL summary()
  CALL memory_report()
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft()
  !
  ! ... generate reciprocal-lattice vectors and fft indices
  !
  CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart, no_global_sort = smallmem )
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms )
  !$acc update device(mill, g, gg)
  !
  IF (gamma_only) THEN
     ! ... Solvers need to know gstart
     call export_gstart_2_solvers(gstart)
  END IF
  !
  IF (do_comp_esm) CALL esm_init(.NOT. lrism)
  !
  ! ... setup the 2D cutoff factor
  !
  IF (do_cutoff_2D) CALL cutoff_fact()
  !
  ! ... setup two chemical potentials calculation
  !
  IF (twochem) CALL init_twochem()
  !
  CALL gshells ( lmovecell )
  !
  ! ... variable initialization for parallel symmetrization
  !
  CALL sym_rho_init ( gamma_only )
  !
  ! ... allocate memory for all other arrays (potentials, wavefunctions etc)
  !
  call allocate_uspp(use_gpu,noncolin,lspinorb,tqr,nhm,nsp,nat,nspin)
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
  ENDIF
  CALL allocate_locpot()
  CALL allocate_bp_efield()
  CALL bp_global_map()
  !
  IF (lrism) CALL rism_alloc3d()
  !
  call plugin_initbase()
#if defined (__LEGACY_PLUGINS)
  CALL plugin_initbase()
#endif 
#if defined (__ENVIRON)
  IF (use_environ) THEN
    IF (alat < 1.D-8) CALL errore('init_run', "Wrong alat", 1)
    at_scaled = at * alat
    gcutm_scaled = gcutm / alat**2
    call init_environ_base(at_scaled, gcutm_scaled, do_comp_mt)
  END IF
#endif
  !
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ), btype( nbnd, nkstot ) )
  !$acc enter data create(et)
  !
  et(:,:) = 0.D0
  !$acc kernels
  et(:,:) = 0.D0
  !$acc end kernels
  wg(:,:) = 0.D0
  btype(:,:) = 1
  !
  IF (ts_vdw .or. mbd_vdw) THEN
     CALL tsvdw_initialize()
     CALL set_h_ainv()
  END IF
  IF (mbd_vdw) THEN
     CALL init_mbd( nks_start, nk1, nk2, nk3, k1, k2, k3, lforce, tstress )
  END IF
  !
  CALL allocate_wfc_k()
  CALL openfil()
  !
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_init_libxc( nspin, domag )
  IF (dft_has_finite_size_correction()) &
       CALL xclib_set_finite_size_volume(REAL(omega*nk1*nk2*nk3))
  IF ( xclib_dft_is('hybrid') ) THEN
     IF ( lmovecell ) CALL infomsg('iosys', &
          'Variable cell and hybrid XC little tested')
     CALL aceinit0()
  END IF
  !
  CALL hinit0()
  !
  CALL potinit()
  !
  CALL newd()
  !
  CALL wfcinit()
  !
  IF(use_wannier) CALL wannier_init()
  !
#if defined(__MPI)
  ! Cleanup PAW arrays that are only used for init
  IF (okpaw) CALL paw_post_init() ! only parallel!
#endif
  !
  IF ( lmd ) CALL allocate_dyn_vars()
  !
  CALL stop_clock( 'init_run' )
  !
  RETURN
  !
END SUBROUTINE init_run
  !
!----------------------------------------------------------------------------
SUBROUTINE pre_init()
  !----------------------------------------------------------------------------
  !
  USE ions_base,        ONLY : nat, ityp
  USE uspp_param,       ONLY : upf, nh, init_uspp_dims
  USE uspp,             ONLY : nkb, nkbus
  IMPLICIT NONE
  INTEGER :: na, nt
  !
  ! init_uspp_dims allocates and fills nh
  !                calculates nhm, nbetam, nwfcm, lmaxkb, lmaxq
  !
  CALL init_uspp_dims( )
  !
  ! calculate the number of beta functions of the solid
  !
  nkb = 0
  nkbus = 0
  do na = 1, nat
     nt = ityp(na)
     nkb = nkb + nh (nt)
     if (upf(nt)%tvanp) nkbus = nkbus + nh (nt)
  enddo
  !
END SUBROUTINE pre_init
