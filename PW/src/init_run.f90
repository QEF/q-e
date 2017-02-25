!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE init_run()
  !----------------------------------------------------------------------------
  !
  USE klist,              ONLY : nkstot
  USE symme,              ONLY : sym_rho_init
  USE wvfct,              ONLY : nbnd, et, wg, btype
  USE control_flags,      ONLY : lmd, gamma_only, smallmem, ts_vdw
  USE cell_base,          ONLY : at, bg, set_h_ainv
  USE cellmd,             ONLY : lmovecell
  USE dynamics_module,    ONLY : allocate_dyn_vars
  USE paw_variables,      ONLY : okpaw
  USE paw_init,           ONLY : paw_init_onecenter, allocate_paw_internals
#if defined(__MPI)
  USE paw_init,           ONLY : paw_post_init
#endif
  USE bp,                 ONLY : allocate_bp_efield, bp_global_map
  USE fft_base,           ONLY : dffts, dtgs
  USE funct,              ONLY : dft_is_hybrid
  USE recvec_subs,        ONLY : ggen
  USE wannier_new,        ONLY : use_wannier    
  USE dfunct,             ONLY : newd
  USE esm,                ONLY : do_comp_esm, esm_init
  USE mp_bands,           ONLY : intra_bgrp_comm, inter_bgrp_comm, nbgrp, root_bgrp_id
  USE mp,                 ONLY : mp_bcast
  USE tsvdw_module,       ONLY : tsvdw_initialize
  USE wavefunctions_module, ONLY : evc
#if defined(__HDF5)
  USE hdf5_qe, ONLY : initialize_hdf5
#endif
  !
  IMPLICIT NONE
  !
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
  IF( smallmem ) THEN
     CALL ggen( gamma_only, at, bg, intra_bgrp_comm, no_global_sort = .TRUE. )
  ELSE
     CALL ggen( gamma_only, at, bg )
  END IF
  !
  IF (do_comp_esm) CALL esm_init()
  CALL gshells ( lmovecell )
  !
  ! ... variable initialization for parallel symmetrization
  !
  CALL sym_rho_init (gamma_only )
  !
  ! ... allocate memory for all other arrays (potentials, wavefunctions etc)
  !
  CALL allocate_nlpot()
  IF (okpaw) THEN
    CALL allocate_paw_internals()
    CALL paw_init_onecenter()
  ENDIF
  CALL allocate_locpot()
  CALL allocate_wfc()
  CALL allocate_bp_efield()
  CALL bp_global_map()
  !
  call plugin_initbase()
  !
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ), btype( nbnd, nkstot ) )
  !
  et(:,:) = 0.D0
  wg(:,:) = 0.D0
  !
  btype(:,:) = 1
  !
  IF (ts_vdw) THEN
    CALL tsvdw_initialize()
    CALL set_h_ainv()
  END IF
  !
  CALL openfil()
  !
  CALL hinit0()
  !
  CALL potinit()
  !
  CALL newd()
#if defined(__HDF5)
  ! calls h5open_f mandatory in any application using hdf5
  CALL initialize_hdf5()
#endif 
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
  IF( nbgrp > 1 ) THEN
     CALL mp_bcast( evc, root_bgrp_id, inter_bgrp_comm )
  ENDIF
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
  USE ions_base,        ONLY : nat, nsp, ityp
  USE uspp_param,       ONLY : upf, lmaxkb, nh, nhm, nbetam
  USE uspp,             ONLY : nkb, nkbus
  IMPLICIT NONE
  INTEGER :: na, nt, nb
  !
  !     calculate the number of beta functions for each atomic type
  !
  lmaxkb = - 1
  DO nt = 1, nsp
     !
     nh (nt) = 0
     !
     ! do not add any beta projector if pseudo in 1/r fmt (AF)
     IF ( upf(nt)%tcoulombp ) CYCLE 
     !
     DO nb = 1, upf(nt)%nbeta
        nh (nt) = nh (nt) + 2 * upf(nt)%lll(nb) + 1
        lmaxkb = MAX (lmaxkb, upf(nt)%lll(nb) )
     ENDDO
     !
  ENDDO
  !
  ! calculate the maximum number of beta functions
  !
  nhm = MAXVAL (nh (1:nsp))
  nbetam = MAXVAL (upf(:)%nbeta)
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


END SUBROUTINE pre_init
