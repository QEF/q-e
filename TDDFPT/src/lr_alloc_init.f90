!-----------------------------------------------------------------------
SUBROUTINE lr_alloc_init()
  !---------------------------------------------------------------------
  ! ... allocates and initialises linear response variables
  !---------------------------------------------------------------------
  ! Modified by Osman Baris Malcioglu in 2009
  !
  USE fft_base,             ONLY : dfftp
  USE fft_base,             ONLY : dffts
  USE klist,                ONLY : nks
  USE lr_variables
  USE uspp,                 ONLY : nkb
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx, nbnd
  USE control_flags,        ONLY : gamma_only
  USE io_global,            ONLY : stdout
  USE charg_resp,           ONLY : w_T, w_T_beta_store, w_T_gamma_store, &
                                  & w_T_zeta_store, w_T_npol,chi
  USE realus,               ONLY : igk_k, npw_k, tg_psic
  USE control_ph,           ONLY : nbnd_occ
  USE noncollin_module,     ONLY : nspin_mag
  USE eqv,                  ONLY : dmuxc
  USE wavefunctions_module, ONLY : evc
  USE kinds,                ONLY : dp
  USE control_ph,           ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  IF (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_alloc_init>")')
  ENDIF
  !
  IF (nbnd>nbnd_occ(1)) THEN
     if(.not. davidson)  WRITE(stdout,'(/,5X,"Warning: There are virtual&
           states in the input file, trying to disregard in response calculation")')
     nbnd_total=nbnd
     nbnd=nbnd_occ(1)
  ELSE
     if(davidson) WRITE(stdout,'(/5X,"(XC.G): No virtrual state!! For analysis the property of transitions, &
               & it is suggested to calculate some virtual states. ")')
     nbnd_total=nbnd
  ENDIF

  IF (lr_verbosity > 7) THEN
   WRITE(stdout,'("NPWX=",I15)') npwx
   WRITE(stdout,'("NBND=",I15)') nbnd
   WRITE(stdout,'("NKS=",I15)') nks
   WRITE(stdout,'("NRXX=",I15)') dfftp%nnr
   WRITE(stdout,'("NSPIN_MAG=",I15)') nspin_mag
  ENDIF

  IF (allocated(evc)) THEN
   DEALLOCATE(evc)
   ALLOCATE(evc(npwx,nbnd))
  ENDIF
  ALLOCATE(evc0(npwx,nbnd,nks))
  ALLOCATE(sevc0(npwx,nbnd,nks))
  IF (project .or. davidson) THEN
   WRITE(stdout,'(5x,"Allocating ",I5," extra bands for projection")') nbnd_total-nbnd
   ALLOCATE(F(nbnd,(nbnd_total-nbnd),n_ipol))
   ALLOCATE(R(nbnd,(nbnd_total-nbnd),n_ipol))
   ALLOCATE(chi(3,3))
   chi(:,:)=cmplx(0.0d0,0.0d0,dp)
   F(:,:,:)=cmplx(0.0d0,0.0d0,dp)
   R(:,:,:)=cmplx(0.0d0,0.0d0,dp)
  ENDIF
  !
  ALLOCATE(evc1_old(npwx,nbnd,nks,2))
  ALLOCATE(evc1(npwx,nbnd,nks,2))
  ALLOCATE(evc1_new(npwx,nbnd,nks,2))
  !if(pseudo_hermitian) then
    ALLOCATE(sevc1_new(npwx,nbnd,nks,2))
    sevc1_new(:,:,:,:)=(0.0d0,0.0d0)
  !else
    ALLOCATE(sevc1(npwx,nbnd,nks,2))
    sevc1(:,:,:,:)=(0.0d0,0.0d0)
  !endif
  ALLOCATE(d0psi(npwx,nbnd,nks,n_ipol))
  !
  IF(dffts%have_task_groups) THEN
     ALLOCATE(tg_revc0(dffts%tg_nnr * dffts%nogrp,nbnd,nks))
     IF(.NOT. ALLOCATED(tg_psic)) &
          & ALLOCATE( tg_psic(dffts%tg_nnr * dffts%nogrp) )
  ELSE
     ALLOCATE(revc0(dffts%nnr,nbnd,nks))
  ENDIF
  !
  IF(gamma_only) THEN
     ALLOCATE(rho_1(dfftp%nnr,nspin_mag))
     rho_1(:,:)=0.0d0
  ELSE
     ALLOCATE(rho_1c(dfftp%nnr,nspin_mag))
     rho_1c(:,:)=(0.0d0,0.0d0)
  ENDIF
  !allocate(rho_tot(dfftp%nnr))
  IF (charge_response == 1 ) THEN
   !allocate(rho_1_tot(dfftp%nnr,nspin_mag)) !Due to broadening this is now done in lr_charg_resp
   !rho_1_tot(:,:)=0.0d0
   !print *,"allocating beta w_t"
   ALLOCATE(w_T_beta_store(itermax_int))
   ALLOCATE(w_T_gamma_store(itermax_int))
   ALLOCATE(w_T_zeta_store(w_T_npol,itermax_int))
   ALLOCATE(w_T(itermax_int))
   w_T_gamma_store(:)=0.0d0
   w_T_beta_store(:)=0.0d0
   w_T_zeta_store(:,:)=cmplx(0.0d0,0.0d0,dp)

  ENDIF
  !if (charge_response /=0) then
  ! allocate(w_T(itermax_int))
  !endif

  ALLOCATE(dmuxc ( dfftp%nnr , nspin , nspin))
  !print *, "dmuxc ALLOCATED",allocated(dmuxc)," SIZE=",size(dmuxc)
  !print *, "nks=",nks

  !allocate (nbnd_occ (nks))
  !

  evc0(:,:,:)=(0.0d0,0.0d0)
  evc1_old(:,:,:,:)=(0.0d0,0.0d0)
  evc1(:,:,:,:)=(0.0d0,0.0d0)
  evc1_new(:,:,:,:)=(0.0d0,0.0d0)
  !rho_tot(:)=0.0d0
  d0psi(:,:,:,:)=(0.0d0,0.0d0)
  !
  ALLOCATE(alpha_store(n_ipol,itermax))
  ALLOCATE(beta_store(n_ipol,itermax))
  ALLOCATE(gamma_store(n_ipol,itermax))
  ALLOCATE(zeta_store(n_ipol,n_ipol,itermax))
  alpha_store(:,:)=0.0d0
  beta_store(:,:)=0.0d0
  gamma_store(:,:)=0.0d0
  zeta_store(:,:,:)=(0.0d0,0.0d0)
  !
  IF(gamma_only) THEN
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
    USE becmod,               ONLY : allocate_bec_type, bec_type, becp
    !
    IF (nkb > 0) THEN
#ifdef __STD_F95
       IF (.not. associated(becp%r)) CALL allocate_bec_type(nkb,nbnd,becp)
#else
       IF (.not. allocated(becp%r)) CALL allocate_bec_type(nkb,nbnd,becp)
#endif
       becp%r(:,:)=0.0d0
       ALLOCATE(becp1(nkb,nbnd))
       becp1(:,:)=0.0d0
       IF (project .or. davidson) THEN
        ALLOCATE(becp1_virt(nkb,nbnd_total-nbnd))
        becp1_virt(:,:)=0.0d0
       ENDIF
    ENDIF
    !
    RETURN
  END SUBROUTINE lr_alloc_init_gamma
  !
  SUBROUTINE lr_alloc_init_k()
    USE becmod,               ONLY : allocate_bec_type, bec_type, becp
    !
    IF (nkb > 0) THEN
#ifdef __STD_F95
       IF(.not. associated(becp%k)) CALL allocate_bec_type(nkb,nbnd,becp)
#else
       IF(.not. allocated(becp%k)) CALL allocate_bec_type(nkb,nbnd,becp)
#endif
       becp%k(:,:)=(0.0d0,0.0d0)
       ALLOCATE(becp1_c(nkb,nbnd,nks))
       becp1_c(:,:,:)=(0.0d0,0.0d0)
       IF (project .or. davidson) THEN
        ALLOCATE(becp1_c_virt(nkb,nbnd_total-nbnd,nks))
        becp1_c_virt(:,:,:)=(0.0d0,0.0d0)
       ENDIF
    ENDIF
    !
    RETURN
  END SUBROUTINE lr_alloc_init_k
  !
END SUBROUTINE lr_alloc_init
!----------------------------------------------------------------------------
