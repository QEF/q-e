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
  USE control_flags,      ONLY : lmd, gamma_only
  USE cell_base,          ONLY : at, bg
  USE cellmd,             ONLY : lmovecell
  USE dynamics_module,    ONLY : allocate_dyn_vars
  USE paw_variables,      ONLY : okpaw
  USE paw_init,           ONLY : paw_init_onecenter, allocate_paw_internals
#ifdef __PARA
  USE paw_init,           ONLY : paw_post_init
#endif
  USE bp,                 ONLY : lberry, lelfield
#ifdef __SOLVENT
  USE fft_base,           ONLY : dfftp
  USE solvent_base,       ONLY : do_solvent
#endif
  USE recvec_subs,        ONLY : ggen
! DCC
!  USE grid_dimensions,    ONLY : nr1x, nr2x, nr3x, nr1, nr2, nr3
!  USE gvect,              ONLY : ecutwfc
!  USE ee_mod,             ONLY : do_comp, do_coarse
! Wannier_ac
  USE wannier_new,        ONLY : use_wannier    
  USE dfunct,                 only : newd
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
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft()
  !
  ! ... generate reciprocal-lattice vectors and fft indices
  !
  CALL ggen ( gamma_only, at, bg )
  CALL gshells ( lmovecell )
  !
  ! ... variable initialization for parallel symmetrization
  !
  CALL sym_rho_init (gamma_only )
  !
  CALL summary()
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
  IF( lberry .or. lelfield) call bp_global_map()
#ifdef __SOLVENT
  IF ( do_solvent ) CALL solvent_initbase( dfftp%nnr )
#endif
! DCC
  ! ... Initializes EE variables
  !
!  IF ( do_comp ) CALL init_ee(nr1x,nr2x,nr3x)
  !
!  IF ( do_coarse )  THEN
!    CALL ggen_coarse()
!    CALL data_structure_coarse( gamma_only, nr1,nr2,nr3, ecutwfc )
!  END IF

  CALL memory_report()
  !
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ), btype( nbnd, nkstot ) )
  !
  et(:,:) = 0.D0
  wg(:,:) = 0.D0
  !
  btype(:,:) = 1
  !
  CALL openfil()
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
#ifdef __PARA
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
