!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_alloc_init()
  !---------------------------------------------------------------------
  !
  ! This subroutine allocates and initialises linear-response variables.
  !
  USE lr_variables
  USE kinds,                ONLY : dp
  USE ions_base,            ONLY : nat
  USE uspp,                 ONLY : nkb, okvan
  USE uspp_param,           ONLY : nhm
  USE fft_base,             ONLY : dfftp, dffts
  USE klist,                ONLY : nks
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx, nbnd
  USE control_flags,        ONLY : gamma_only
  USE io_global,            ONLY : stdout
  USE charg_resp,           ONLY : w_T, w_T_beta_store, w_T_gamma_store, &
                                 & w_T_zeta_store, w_T_npol,chi
  USE realus,               ONLY : igk_k, npw_k, tg_psic
  USE control_lr,           ONLY : nbnd_occ
  USE noncollin_module,     ONLY : nspin_mag, npol, noncolin
  USE eqv,                  ONLY : dmuxc, evq, dpsi, dvpsi
  USE wavefunctions_module, ONLY : evc
  USE qpoint,               ONLY : igkq, nksq, eigqts
  USE becmod,               ONLY : allocate_bec_type, bec_type, becp
  USE phus,                 ONLY : int3, int3_nc, becp1
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_alloc_init>")')
  ENDIF
  !
  ! Optical case
  !
  IF (.NOT.eels) THEN
     !
     IF (nbnd>nbnd_occ(1)) THEN
        !
        IF (.not. davidson)  WRITE(stdout,'(/,5X,"Warning: There are virtual&
            & states in the input file, trying to disregard in response calculation")')
        nbnd_total = nbnd
        nbnd = nbnd_occ(1)
        !
     ELSE
        !
        IF (davidson) WRITE(stdout,'(/5X,"(XC.G): No virtrual states! To analyse the &
           & property of transitions, it is suggested to calculate some virtual states. ")')
        nbnd_total = nbnd
        !
     ENDIF
     !
  ENDIF
  !
  IF (lr_verbosity > 7) THEN
     WRITE(stdout,'("NPWX=",I15)') npwx
     WRITE(stdout,'("NBND=",I15)') nbnd
     WRITE(stdout,'("NKS=",I15)') nks
     WRITE(stdout,'("NRXX=",I15)') dfftp%nnr
     WRITE(stdout,'("NSPIN_MAG=",I15)') nspin_mag
  ENDIF
  !
  IF (allocated(evc)) THEN
     DEALLOCATE(evc)
     ALLOCATE(evc(npwx*npol,nbnd))
  ENDIF 
  !
  ! Allocate unperturbed and perturbed orbitals
  !
  ALLOCATE(evc0(npwx*npol,nbnd,nksq))
  !
  IF (.NOT.eels) THEN
     ALLOCATE(sevc0(npwx*npol,nbnd,nks))
     sevc0(:,:,:) = (0.0d0,0.0d0)
  ENDIF
  !
  IF (project .or. davidson) THEN
     !
     WRITE(stdout,'(5x,"Allocating ",I5," extra bands for projection")') nbnd_total-nbnd
     ALLOCATE(F(nbnd,(nbnd_total-nbnd),n_ipol))
     ALLOCATE(R(nbnd,(nbnd_total-nbnd),n_ipol))
     ALLOCATE(chi(3,3))
     chi(:,:) = cmplx(0.0d0,0.0d0,dp)
     F(:,:,:) = cmplx(0.0d0,0.0d0,dp)
     R(:,:,:) = cmplx(0.0d0,0.0d0,dp) 
     !
  ENDIF
  !
  ALLOCATE(evc1_old(npwx*npol,nbnd,nksq,2))
  ALLOCATE(evc1(npwx*npol,nbnd,nksq,2))
  ALLOCATE(evc1_new(npwx*npol,nbnd,nksq,2))
  !
  ALLOCATE(sevc1_new(npwx*npol,nbnd,nksq,2))
  sevc1_new(:,:,:,:) = (0.0d0,0.0d0)
  !
  if (.not. pseudo_hermitian) then
    ALLOCATE(sevc1(npwx*npol,nbnd,nksq,2))
    sevc1(:,:,:,:) = (0.0d0,0.0d0)
  endif
  !
  ALLOCATE(d0psi(npwx*npol,nbnd,nksq,n_ipol))
  !
  evc0(:,:,:)       = (0.0d0,0.0d0)
  evc1_old(:,:,:,:) = (0.0d0,0.0d0)
  evc1(:,:,:,:)     = (0.0d0,0.0d0)
  evc1_new(:,:,:,:) = (0.0d0,0.0d0)
  d0psi(:,:,:,:)    = (0.0d0,0.0d0)
  !
  IF (eels) THEN
     ALLOCATE(d0psi2(npwx*npol,nbnd,nksq,n_ipol))
     d0psi2(:,:,:,:) = (0.0d0,0.0d0)
  ENDIF
  !
  ! Allocate the R-space unperturbed orbitals
  !
  IF (dffts%have_task_groups) THEN
     ALLOCATE(tg_revc0(dffts%tg_nnr * dffts%nogrp,nbnd,nksq))
     IF (.NOT. ALLOCATED(tg_psic)) &
          & ALLOCATE( tg_psic(dffts%tg_nnr * dffts%nogrp) )
  ELSE
     IF (.NOT.eels) THEN
        ALLOCATE(revc0(dffts%nnr,nbnd,nksq))
        revc0(:,:,:) = (0.0d0,0.0d0)
     ENDIF
  ENDIF
  !
  ! Optical case: allocate the response charge-density
  !
  IF (.NOT.eels) THEN
     ! 
     IF (gamma_only) THEN
        ALLOCATE(rho_1(dfftp%nnr,nspin_mag))
        rho_1(:,:)=0.0d0
     ELSE
        ALLOCATE(rho_1c(dfftp%nnr,nspin_mag))
        rho_1c(:,:)=(0.0d0,0.0d0)
     ENDIF
     !
  ENDIF
  !
  ! EELS: evq, igkq are allocated and calculated at k+q 
  ! 
  IF (eels) THEN
     !
     ALLOCATE (evq(npwx*npol,nbnd))
     ALLOCATE (dpsi(npwx*npol,nbnd))
     ALLOCATE (dvpsi(npwx*npol,nbnd))
     ALLOCATE (igkq(npwx))
     evq(:,:)   = (0.0d0, 0.0d0)
     dpsi(:,:)  = (0.0d0, 0.0d0)
     dvpsi(:,:) = (0.0d0, 0.0d0)
     !
     IF (okvan) THEN
        !
        ALLOCATE (int3(nhm,nhm,nat,nspin_mag,1))
        int3 = (0.0d0, 0.0d0)
        !
        IF (noncolin) THEN
           ALLOCATE (int3_nc(nhm,nhm,nat,nspin,1))
           int3_nc = (0.0d0, 0.0d0)
        ENDIF
        !
     ENDIF
     !
  ENDIF
  !
  ! Optical case
  !
  IF (charge_response == 1 .AND. .NOT.eels) THEN
     !
     ALLOCATE(w_T_beta_store(itermax_int))
     ALLOCATE(w_T_gamma_store(itermax_int))
     ALLOCATE(w_T_zeta_store(w_T_npol,itermax_int))
     ALLOCATE(w_T(itermax_int))
     w_T_gamma_store(:)=0.0d0
     w_T_beta_store(:)=0.0d0
     w_T_zeta_store(:,:)=cmplx(0.0d0,0.0d0,dp)
     !
  ENDIF
  !
  ! Allocate an array for a derivative of the XC potential
  !
  ALLOCATE(dmuxc (dfftp%nnr, nspin, nspin))
  !
  ! Lanczos coefficients and zeta coefficients
  !
  ALLOCATE(alpha_store(n_ipol,itermax))
  ALLOCATE(beta_store(n_ipol,itermax))
  ALLOCATE(gamma_store(n_ipol,itermax))
  ALLOCATE(zeta_store(n_ipol,n_ipol,itermax))
  alpha_store(:,:)  = 0.0d0
  beta_store(:,:)   = 0.0d0
  gamma_store(:,:)  = 0.0d0
  zeta_store(:,:,:) = (0.0d0,0.0d0)
  !
  IF (gamma_only) THEN
     CALL lr_alloc_init_gamma()
  ELSE
     CALL lr_alloc_init_k()
  ENDIF
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_alloc_init_gamma()
    !
    IF (nkb > 0) THEN
       !
       IF (.not. allocated(becp%r)) CALL allocate_bec_type(nkb,nbnd,becp)
       becp%r(:,:) = 0.0d0
       !
       ALLOCATE(becp_1(nkb,nbnd))
       becp_1(:,:) = 0.0d0
       !
       IF (project .or. davidson) THEN
          ALLOCATE(becp1_virt(nkb,nbnd_total-nbnd))
          becp1_virt(:,:)=0.0d0
       ENDIF
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE lr_alloc_init_gamma
  !
  SUBROUTINE lr_alloc_init_k()
    !
    IF (nkb > 0) THEN
       !
       IF(.not. allocated(becp%k)) CALL allocate_bec_type(nkb,nbnd,becp)
       becp%k(:,:) = (0.0d0,0.0d0)
       !
       IF (.NOT.eels) THEN
          ALLOCATE(becp1_c(nkb,nbnd,nks))
          becp1_c(:,:,:) = (0.0d0,0.0d0)
       ENDIF
       !
       IF (project .or. davidson) THEN
          ALLOCATE(becp1_c_virt(nkb,nbnd_total-nbnd,nks))
          becp1_c_virt(:,:,:) = (0.0d0,0.0d0)
       ENDIF
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE lr_alloc_init_k
  !
END SUBROUTINE lr_alloc_init
!----------------------------------------------------------------------------
