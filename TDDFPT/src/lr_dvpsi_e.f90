!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
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
  ! dvpsi is COMPUTED and WRITTEN on file (vkb,evc,igk must be set) 
  ! OBM:                  ^ This is now handled elesewhere
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Rebased wrt PHONON routines. S J Binnie (2011)
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba, at
  USE ions_base,            ONLY : ntyp => nsp
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : npw, npwx, nbnd, igk, g2kin, et
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : noncolin, npol
  USE becmod,               ONLY : allocate_bec_type, calbec, becp, &
                                   & deallocate_bec_type,  bec_type
  USE uspp,                 ONLY : okvan, nkb, vkb
  USE uspp_param,           ONLY : nh, nhm
  USE phus,                 ONLY : dpqq 
  USE control_flags,        ONLY : gamma_only
  USE control_ph,           ONLY : nbnd_occ
  USE realus,               ONLY : npw_k
  USE lr_variables,         ONLY : lr_verbosity
  USE io_global,            ONLY : stdout
  USE qpoint,               ONLY : igkq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ipol, ik
  COMPLEX(kind=dp), INTENT(out) :: dvpsi(npwx*npol,nbnd)
  !
  ! Local variables
  !  
  INTEGER :: ig, ibnd, lter
  ! counters
  REAL(kind=dp) :: atnorm
  COMPLEX(kind=dp),ALLOCATABLE :: d0psi(:,:)
  REAL(DP), ALLOCATABLE  :: h_diag (:,:), eprec(:)
  ! diagonal part of h_scf
  real(DP) ::   anorm
  ! preconditioning cut-off
  REAL(DP), PARAMETER :: thresh = 1.0e-5_DP
  ! the desired convergence of linter
  LOGICAL :: conv_root
  ! true if convergence has been achieved
  COMPLEX(DP), ALLOCATABLE :: work (:,:), spsi(:,:)
  !
  TYPE(bec_type) :: becp1 
  TYPE(bec_type) :: becp2 
  !
  EXTERNAL ch_psi_all, cg_psi
  !
  CALL start_clock ('lr_dvpsi_e')
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<lr_dvpsi_e>")')
  !
  conv_root = .TRUE.
  !
  ALLOCATE(d0psi(npwx*npol,nbnd))
  d0psi = (0.d0, 0.d0)
  dvpsi = (0.d0, 0.d0)
  !
  ALLOCATE (h_diag( npwx*npol, nbnd))
  h_diag = 0.d0
  !
  CALL allocate_bec_type ( nkb, nbnd, becp1 )
  !
  CALL calbec ( npw_k(ik), vkb, evc, becp1 )
  !
  CALL allocate_bec_type ( nkb, nbnd, becp2 )
  !
  CALL commutator_Hx_psi (ik, nbnd_occ(ik), becp1, becp2, ipol, d0psi )
  !
  !    orthogonalize d0psi to the valence subspace: ps = <evc|d0psi>
  !    Apply -P^+_c
  !    NB it uses dvpsi as workspace
  !
  IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd)
  !
  CALL orthogonalize(d0psi, evc, ik, ik, dvpsi,npw_k(ik))
  d0psi = -d0psi
  !
  !   d0psi contains P^+_c [H-eS,x] psi_v for the polarization direction ipol
  !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
  !
  ! eprec is now calculated on the fly for each k point
  !
  ALLOCATE(eprec(nbnd))
  CALL lr_calc_eprec(eprec)
  !
  DO ibnd = 1, nbnd_occ (ik)
     DO ig = 1, npw_k(ik)
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
     ENDDO
     IF (noncolin) THEN
        DO ig = 1, npw_k(ik)
           h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
        ENDDO
     ENDIF
  ENDDO
  !
  ! OBM: upto here, the dvpsi was used as a scratch
  !
  dvpsi(:,:) = (0.d0, 0.d0)
  !
  igkq => igk ! PG: needed by h_psiq, called by ch_psi_all
  !
  CALL cgsolve_all (ch_psi_all, cg_psi, et (1, ik), d0psi, dvpsi, &
       h_diag, npwx, npw_k(ik), thresh, ik, lter, conv_root, anorm, &
       nbnd_occ(ik), 1)
  !
  IF (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
       & " linter: root not converged ",e10.3)') &
       ik, ibnd, anorm
  !
  CALL flush_unit( stdout )
  !
  ! we have now obtained P_c x |psi>.
  ! In the case of USPP this quantity is needed for the Born
  ! effective charges, so we save it to disc
  !
  ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
  ! therefore we apply S again, and then subtract the additional term
  ! furthermore we add the term due to dipole of the augmentation charges.
  !
  IF (okvan) THEN
     !
     ! for effective charges
     !
     ALLOCATE (spsi ( npwx*npol, nbnd))
     CALL calbec (npw_k(ik), vkb, dvpsi, becp )
     CALL s_psi(npwx,npw_k(ik),nbnd,dvpsi,spsi)
     CALL DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
     DEALLOCATE (spsi)
     ALLOCATE (dpqq( nhm, nhm, 3, ntyp))
     CALL compute_qdipol(dpqq)
     CALL qdipol_cryst()
     CALL adddvepsi_us(becp1,becp2,ipol,ik,dvpsi)
     DEALLOCATE (dpqq)
     !
  ENDIF
  !
  ! orthogonalize dvpsi to the valence subspace
  !
  IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd)
  !
  ALLOCATE (work ( npwx, nbnd ) )
  work = evc ! work will be corrupted on exit from orthogonalize
  !
  CALL orthogonalize(dvpsi, evc, ik, ik, work,npw_k(ik))
  dvpsi = -dvpsi
  !
  DEALLOCATE(work)
  DEALLOCATE (h_diag)
  DEALLOCATE (eprec)
  DEALLOCATE (d0psi)
  !
  ! OBM: Addendum to PH dvpsi
  !
  IF (okvan) THEN
     ALLOCATE (spsi ( npwx*npol, nbnd))
     CALL sm1_psi(.TRUE.,ik,npwx,npw_k(ik),nbnd,dvpsi,spsi)
     dvpsi(:,:) = spsi(:,:)
     DEALLOCATE(spsi)
  ENDIF
  !
  ! For some ibrav the crystal axes are not normalized
  ! Here we include the correct normalization
  ! for Lanczos initial wfcs
  !
  atnorm = dsqrt(at(1,ipol)**2+at(2,ipol)**2+at(3,ipol)**2)
  !
  dvpsi(:,:) = dvpsi(:,:)/atnorm
  !
  ! nrec = (ipol - 1)*nksq + ik
  ! call davcio(dvpsi, lrebar, iuebar, nrec, 1)
  ! this_pcxpsi_is_on_file(ik,ipol) = .true.
  !
  CALL deallocate_bec_type ( becp1 )
  IF (nkb > 0) CALL deallocate_bec_type ( becp2 )
  !
  CALL stop_clock ('lr_dvpsi_e')
  !
  RETURN
  !
CONTAINS

  SUBROUTINE lr_calc_eprec(eprec)

    USE kinds,                ONLY : DP
    USE gvect,                ONLY : gstart
    USE wvfct,                ONLY : npw, npwx, nbnd, g2kin
    USE wavefunctions_module, ONLY : evc
    USE realus,               ONLY : npw_k
    USE mp,                   ONLY : mp_sum
    USE mp_global,            ONLY : intra_bgrp_comm

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(INOUT) :: eprec(nbnd)
    COMPLEX(KIND=DP), ALLOCATABLE :: work (:,:)
    REAL(KIND=DP), EXTERNAL :: ddot
    COMPLEX(KIND=DP), EXTERNAL :: ZDOTC
    ! the scalar products
    !
    ALLOCATE (work(npwx,nbnd))
    !
    DO ibnd=1,nbnd
       !
       work = 0.d0
       !
       DO ig = 1,npw_k(ik)
          work(ig,1) = g2kin(ig)*evc(ig,ibnd)
       ENDDO
       !
       IF (gamma_only) THEN
          !
          eprec(ibnd) = 2.0d0*DDOT(2*npw,evc(1,ibnd),1,work,1)
          !
          IF (gstart==2) THEN
             eprec(ibnd) = eprec(ibnd)-DBLE(evc(1,ibnd))*DBLE(work(1,ibnd))
          ENDIF
          !
          eprec(ibnd) = 1.35d0*eprec(ibnd)
          !
       ELSE
          eprec(ibnd) = 1.35d0*ZDOTC(npw_k(ik),evc(1,ibnd),1,work,1)
       ENDIF
       !
    ENDDO
    !
#ifdef __MPI
    CALL mp_sum(eprec, intra_bgrp_comm)
#endif
    !
    DEALLOCATE(work)
    !
    RETURN
    !
  END SUBROUTINE lr_calc_eprec

END SUBROUTINE lr_dvpsi_e


