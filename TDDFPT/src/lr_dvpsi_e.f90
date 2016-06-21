!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE lr_dvpsi_e(ik,ipol,dvpsi)
  !----------------------------------------------------------------------
  !
  ! On output: dvpsi contains P_c^+ x | psi_ik > in crystal axis
  !            (projected on at(*,ipol) )
  !
  ! dvpsi is computed here (vkb and evc must be set) 
  ! and it is written on file in the routine lr_solve_e.
  !
  ! See Ref.[1] : J. Tobik and A. Dal Corso, JCP 120, 9934 (2004)
  ! for the details of the theory implemented in this routine.
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Rebased wrt PHONON routines. S J Binnie (2011)
  ! Modified by Iurii Timrov (2016)
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba2, at
  USE ions_base,            ONLY : ntyp => nsp
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : xk, ngk, igk_k
  USE wvfct,                ONLY : npwx, nbnd, g2kin, et
  USE wavefunctions_module, ONLY : evc
  USE gvect,                ONLY : g
  USE noncollin_module,     ONLY : noncolin, npol
  USE becmod,               ONLY : allocate_bec_type, calbec, becp, &
                                   & deallocate_bec_type,  bec_type
  USE uspp,                 ONLY : okvan, nkb, vkb
  USE uspp_param,           ONLY : nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE control_lr,           ONLY : nbnd_occ
  USE lr_variables,         ONLY : lr_verbosity, sevc0
  USE io_global,            ONLY : stdout
  USE lrus,                 ONLY : dpqq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ipol, ik
  COMPLEX(kind=dp), INTENT(out) :: dvpsi(npwx*npol,nbnd)
  !
  ! Local variables
  !  
  INTEGER :: ig, ibnd, lter, npw
  ! counters
  ! npw: number of plane-waves at point ik
  REAL(kind=dp) :: atnorm
  COMPLEX(kind=dp),ALLOCATABLE :: d0psi(:,:)
  REAL(DP), ALLOCATABLE  :: h_diag (:,:)
  ! diagonal part of h_scf
  real(DP) ::   anorm
  ! preconditioning cut-off
  REAL(DP), PARAMETER :: thresh = 1.0e-5_DP
  ! the desired convergence of linter
  LOGICAL :: conv_root
  ! true if convergence has been achieved
  COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
  !
  TYPE(bec_type) :: becp1 
  TYPE(bec_type) :: becp2 
  !
  EXTERNAL ch_psi_all, cg_psi
  !
  CALL start_clock ('lr_dvpsi_e')
  !
  conv_root = .TRUE.
  !
  ALLOCATE ( d0psi(npwx*npol,nbnd) )
  d0psi = (0.d0, 0.d0)
  dvpsi = (0.d0, 0.d0)
  !
  npw = ngk(ik)
  !
  CALL allocate_bec_type ( nkb, nbnd, becp1 )
  !
  CALL calbec ( npw, vkb, evc, becp1 )
  !
  CALL allocate_bec_type ( nkb, nbnd, becp2 )
  !
  CALL commutator_Hx_psi (ik, nbnd_occ(ik), becp1, becp2, ipol, d0psi )
  !
  IF (okvan) CALL calbec ( npw, vkb, evc, becp, nbnd)
  !
  ! Orthogonalize d0psi to the valence subspace. Apply P_c^+
  !
  CALL orthogonalize(d0psi, evc, ik, ik, sevc0(:,:,ik), npw, .true.)
  d0psi = -d0psi
  !
  ! Calculate the kinetic energy g2kin: (k+G)^2
  !
  CALL g2_kin(ik)
  !
  ! Calculate the preconditioning matrix h_diag used by cgsolve_all
  !
  ALLOCATE ( h_diag(npwx*npol, nbnd) )
  CALL h_prec(ik, evc, h_diag)
  !
  ! d0psi contains P^+_c [H-eS,x] psi_v for the polarization direction ipol
  ! Now solve the linear systems (H+Q-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
  ! See Eq.(9) in Ref. [1]
  !
  CALL cgsolve_all (ch_psi_all, cg_psi, et (1, ik), d0psi, dvpsi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, nbnd_occ(ik), 1)
  !
  IF (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
       & " lr_dvpsi_e: root not converged ",e10.3)') &
       ik, ibnd, anorm
  !
  FLUSH( stdout )
  DEALLOCATE (h_diag)
  !
  ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
  ! therefore we apply S again, and then subtract the additional term
  ! furthermore we add the term due to dipole of the augmentation charges.
  ! See Eq.(10) in Ref. [1]
  !
  IF (okvan) THEN
     ALLOCATE (spsi ( npwx*npol, nbnd))
     CALL calbec (npw, vkb, dvpsi, becp )
     CALL s_psi(npwx,npw,nbnd,dvpsi,spsi)
     CALL DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
     DEALLOCATE (spsi)
     ALLOCATE (dpqq( nhm, nhm, 3, ntyp))
     CALL compute_qdipol(dpqq)
     CALL qdipol_cryst()
     CALL adddvepsi_us(becp1,becp2,ipol,ik,dvpsi)
     DEALLOCATE (dpqq)
  ENDIF
  !
  IF (okvan) CALL calbec ( npw, vkb, evc, becp, nbnd)
  !
  ! Orthogonalize dvpsi to the valence subspace. Apply P_c^+
  !
  CALL orthogonalize(dvpsi, evc, ik, ik, sevc0(:,:,ik), npw, .true.)
  dvpsi = -dvpsi
  !
  DEALLOCATE (d0psi)
  !
  ! US case: apply the S^{-1} operator
  !
  IF (okvan) THEN
     ALLOCATE (spsi ( npwx*npol, nbnd))
     CALL lr_sm1_psi (.TRUE.,ik,npwx,ngk(ik),nbnd,dvpsi,spsi)
     dvpsi(:,:) = spsi(:,:)
     DEALLOCATE(spsi)
  ENDIF
  !
  ! For some ibrav the crystal axes are not normalized
  ! Here we include the correct normalization
  ! for Lanczos initial wfcs
  !
  atnorm = DSQRT(at(1,ipol)**2 + at(2,ipol)**2 + at(3,ipol)**2)
  !
  dvpsi(:,:) = dvpsi(:,:) / atnorm
  !
  CALL deallocate_bec_type ( becp1 )
  IF (nkb > 0) CALL deallocate_bec_type ( becp2 )
  !
  CALL stop_clock ('lr_dvpsi_e')
  !
  RETURN
  !
END SUBROUTINE lr_dvpsi_e


