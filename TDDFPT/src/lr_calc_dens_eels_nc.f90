!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE lr_calc_dens_eels_nc (drhoscf, dpsi)
  !-------------------------------------------------------------------------
  !
  ! Calculates response charge-density from linear-response
  ! orbitals and ground-state orbitals (noncollinear case).
  ! See Eq.(36) in B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
  !
  ! It does:
  ! 1) Calculates the response charge-density
  ! 2) Calculates the contribution to the charge-density due to US PP's
  ! 3) Sums up the normal and ultrasoft terms
  ! 4) Symmetrizes the total response charge-density
  !
  ! Written by Iurii Timrov (2013)
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE gvect,                 ONLY : ngm, g
  USE fft_base,              ONLY : dffts, dfftp
  USE klist,                 ONLY : xk, wk
  USE lr_variables,          ONLY : evc0, lr_periodic
  USE wvfct,                 ONLY : nbnd, npwx, npw, igk, g2kin
  USE gvecw,                 ONLY : gcutw
  USE qpoint,                ONLY : npwq, igkq, nksq, ikks, ikqs
  USE units_ph,              ONLY : lrwfc, iuwfc
  USE wavefunctions_module,  ONLY : evc
  USE noncollin_module,      ONLY : npol, nspin_mag
  USE uspp_param,            ONLY : nhm
  USE uspp,                  ONLY : okvan, vkb
  USE lsda_mod,              ONLY : nspin
  USE mp_global,             ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: dpsi(npwx*npol,nbnd,nksq)
  ! input: the perturbed wavefunctions
  COMPLEX(DP), INTENT(out) :: drhoscf(dfftp%nnr,nspin_mag)
  ! input/output: the accumulated change of the charge density
  !
  COMPLEX(DP), ALLOCATABLE :: drhoscfh(:,:), dbecsum(:,:,:), dbecsum_nc(:,:,:,:)
  ! the accumulated change to the charge density on the smooth mesh
  ! the ultrasoft term
  INTEGER :: ik, ikk, ikq, is
  REAL(DP) :: weight
  ! the effective weight of the k point
  !
  CALL start_clock ('lr_calc_dens')
  !
  ALLOCATE (drhoscfh(dffts%nnr,nspin_mag))
  drhoscfh(:,:) = (0.0d0, 0.0d0)
  !
  IF (okvan) THEN
     ALLOCATE (dbecsum((nhm*(nhm + 1))/2,nat,nspin_mag))
     ALLOCATE (dbecsum_nc(nhm,nhm,nat,nspin))
     dbecsum(:,:,:)      = (0.0d0, 0.0d0)
     dbecsum_nc(:,:,:,:) = (0.0d0, 0.0d0)
  ENDIF
  !
  DO ik = 1, nksq
     !
     IF (lr_periodic) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = ikks(ik)
        ikq = ikqs(ik)
     ENDIF
     !
     ! Determination of npw, igk, and npwq, igkq;
     ! g2kin is used here as work space.
     !
     CALL gk_sort( xk(1,ikk), ngm, g, gcutw, npw,  igk,  g2kin )
     CALL gk_sort( xk(1,ikq), ngm, g, gcutw, npwq, igkq, g2kin )
     !
     ! Read the unperturbed wavefuctions evc(k)
     !
     IF (lr_periodic) THEN
        evc(:,:) = evc0(:,:,ik)
     ELSE
        IF (nksq > 1) CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
     ENDIF
     !
     ! The weight of the k point
     !
     weight = wk(ikk)
     !
     ! Calculate beta-functions vkb at point k+q
     ! 
     IF (okvan) CALL init_us_2 (npwq, igkq, xk(1,ikq), vkb)
     !
     ! Calculation of the response charge density
     !
     CALL incdrhoscf_nc(drhoscfh(:,:), weight, ik, dbecsum_nc(:,:,:,:), dpsi(:,:,ik))
     !
  ENDDO 
  !
  ! Interpolate the charge density from the smooth mesh
  ! to a thicker mesh (if doublegrid=.true.)
  ! drhoscfh -> drhoscf
  !
  DO is = 1, nspin_mag
     CALL cinterpolate(drhoscf(1,is), drhoscfh(1,is), 1)
  ENDDO
  !
  IF (okvan) THEN
     !
     ! The calculation of dbecsum is distributed across processors (see addusdbec_nc).
     ! Sum over processors the contributions coming from each slice of bands.
     !
#ifdef __MPI
     CALL mp_sum (dbecsum_nc, intra_bgrp_comm)
#endif
     !
     ! In the noncolinear, spin-orbit case rotate dbecsum
     !
     CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 1)
     !
     ! Calculate the total response charge-density
     ! (sum up the normal and ultrasoft terms)
     !
     CALL lr_addusddens (drhoscf, dbecsum)
     !
  ENDIF
  !
  ! Reduce the charge density response across pools.
  !
#ifdef __MPI
  CALL mp_sum (drhoscf, inter_pool_comm)
#endif
  !
  ! Symmetrization of the charge density response
  ! wrt the small group of q.
  !
#ifdef __MPI
  CALL lr_psym_eels(drhoscf)
#else
  CALL lr_sym_eels(drhoscf)
#endif
  !
  DEALLOCATE (drhoscfh)
  !
  IF (okvan) THEN
     DEALLOCATE (dbecsum)
     DEALLOCATE (dbecsum_nc)
  ENDIF
  !
  CALL stop_clock ('lr_calc_dens')
  !
  RETURN
  !
END SUBROUTINE lr_calc_dens_eels_nc
