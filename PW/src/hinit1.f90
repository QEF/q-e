!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE hinit1()
  !----------------------------------------------------------------------------
  !! Atomic configuration dependent hamiltonian initialization,
  !! potential, wavefunctions for Hubbard U.  
  !! Important note: it does not recompute structure factors and core charge,
  !! they must be computed before this routine is called.
  !
  USE ions_base,           ONLY : nat, nsp, ityp, tau
  USE cell_base,           ONLY : at, bg, omega, tpiba2
  USE fft_base,            ONLY : dfftp
  USE gvecs,               ONLY : doublegrid
  USE ldaU,                ONLY : lda_plus_u
  USE lsda_mod,            ONLY : nspin
  USE noncollin_module,    ONLY : report
  USE scf,                 ONLY : vrs, vltot, v, kedtau
  USE control_flags,       ONLY : tqr, use_gpu
  USE realus,              ONLY : generate_qpointlist, betapointlist, &
                                  init_realspace_vars, real_space
  USE wannier_new,         ONLY : use_wannier
  USE martyna_tuckerman,   ONLY : tag_wg_corr_as_obsolete
  USE scf,                 ONLY : rho
  USE paw_variables,       ONLY : okpaw, ddd_paw
  USE paw_onecenter,       ONLY : paw_potential
  USE paw_symmetry,        ONLY : paw_symmetrize_ddd
  USE dfunct,              ONLY : newd
  !
  USE scf_gpum,      ONLY : using_vrs
  
  !
  IMPLICIT NONE
  !
  ! these routines can be used to patch quantities that are dependent
  ! on the ions and cell parameters
  !
  CALL plugin_init_ions()
  CALL plugin_init_cell()
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()
  !
  ! ... more position-dependent initializations
  !
  IF ( tqr ) CALL generate_qpointlist()
  !
  IF ( real_space ) THEN
     CALL betapointlist()
     CALL init_realspace_vars()
  ENDIF
  !
  IF ( report /= 0 ) CALL make_pointlists( )
  !
  CALL tag_wg_corr_as_obsolete
  !
  ! ... plugin contribution to local potential
  !
  CALL plugin_scf_potential( rho, .FALSE., -1.d0, vltot )
  !
  ! ... define the total local potential (external+scf)
  !
  CALL using_vrs(1)
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, &
                doublegrid )
  !
  ! ... update the D matrix and the PAW coefficients
  !
  IF ( okpaw ) THEN
     CALL compute_becsum( 1 )
     CALL PAW_potential( rho%bec, ddd_paw )
     CALL PAW_symmetrize_ddd( ddd_paw )
  ENDIF
  ! 
  CALL newd()
  !
  ! ... and recalculate the products of the S with the atomic wfcs used 
  ! ... in LDA+U calculations
  !
  IF (.NOT. use_gpu) THEN
    IF ( lda_plus_u  ) CALL orthoUwfc() 
    IF ( use_wannier ) CALL orthoatwfc( .TRUE. )
  ELSE
    IF ( lda_plus_u  ) CALL orthoUwfc_gpu() 
    IF ( use_wannier ) CALL orthoatwfc_gpu( .TRUE. )
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE hinit1

