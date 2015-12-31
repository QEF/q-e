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
  !
  ! ... Atomic configuration dependent hamiltonian initialization
  ! ... Important note: does not recompute structure factors,
  ! ... they must be computed before this routine is called
  !
  USE ions_base,     ONLY : nat, nsp, ityp, tau
  USE cell_base,     ONLY : at, bg, omega, tpiba2
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm, g
  USE gvecs,         ONLY : doublegrid
  USE ldaU,          ONLY : lda_plus_u
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : vrs, vltot, v, kedtau
  USE control_flags, ONLY : tqr
  USE realus,        ONLY : generate_qpointlist
  USE wannier_new,   ONLY : use_wannier
  USE martyna_tuckerman, ONLY : tag_wg_corr_as_obsolete
  USE scf,           ONLY : rho
  USE paw_variables, ONLY : okpaw, ddd_paw
  USE paw_onecenter, ONLY : paw_potential
  USE paw_symmetry,  ONLY : paw_symmetrize_ddd
  USE dfunct,        ONLY : newd
  !
  IMPLICIT NONE
  !
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()
  !
  ! these routines can be used to patch quantities that are dependent
  ! on the ions and cell parameters
  !
  CALL plugin_init_ions()
  CALL plugin_init_cell()
  !
  ! ... plugin contribution to local potential
  !
  CALL plugin_scf_potential(rho,.FALSE.,-1.d0)
  !
  ! ... define the total local potential (external+scf)
  !
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  !
  IF ( tqr ) CALL generate_qpointlist()
  !
  ! ... update the D matrix and the PAW coefficients
  !
  IF (okpaw) THEN
     CALL compute_becsum(1)
     CALL PAW_potential(rho%bec, ddd_paw)
     CALL PAW_symmetrize_ddd(ddd_paw)
  ENDIF
  ! 
  CALL newd()
  !
  ! ... and recalculate the products of the S with the atomic wfcs used 
  ! ... in LDA+U calculations
  !
  IF ( lda_plus_u ) CALL orthoUwfc () 
  IF ( use_wannier ) CALL orthoatwfc( .true. )
  !
  call tag_wg_corr_as_obsolete
  !
  RETURN
  !
END SUBROUTINE hinit1

