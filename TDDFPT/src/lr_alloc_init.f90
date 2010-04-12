!-----------------------------------------------------------------------
subroutine lr_alloc_init()
  !---------------------------------------------------------------------
  ! ... allocates and initialises linear response variables
  !---------------------------------------------------------------------
  ! Modified by Osman Baris Malcioglu in 2009 
#include "f_defs.h"
  !
  use gvect,             only : nrxx
  use klist,             only : nks
  use lr_variables
  use uspp,              only : nkb
  use lsda_mod,          only : nspin
  use wvfct,             only : npwx, nbnd
  use control_flags,     only : gamma_only
  USE io_global,         ONLY : stdout
  USE charg_resp,        ONLY : w_T, w_T_beta_store, w_T_gamma_store,w_T_zeta_store,w_T_npol,chi
  use realus,            only : igk_k, npw_k
  use control_ph,        ONLY : nbnd_occ
  USE noncollin_module,  ONLY : nspin_mag
  USE eqv,               ONLY : dmuxc
  use wavefunctions_module, only : evc
  use kinds,                only : dp
  !
  implicit none
  !
  If (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_alloc_init>")')
  endif
  if (lr_verbosity > 7) THEN
   WRITE(stdout,'("NPWX=",I15)') npwx
   WRITE(stdout,'("NBND=",I15)') nbnd
   WRITE(stdout,'("NKS=",I15)') nks
   WRITE(stdout,'("NRXX=",I15)') nrxx
   WRITE(stdout,'("NSPIN_MAG=",I15)') nspin_mag
  endif
  !
  if (allocated(evc)) then 
   deallocate(evc)
   allocate(evc(npwx,nbnd))
  endif
  allocate(evc0(npwx,nbnd,nks))
  allocate(sevc0(npwx,nbnd,nks))
  if (project) then
   WRITE(stdout,'(\5x,"Allocating ",I5," extra bands for projection")') nbnd_total-nbnd
   allocate(evc0_virt(npwx,(nbnd_total-nbnd),nks))
   !allocate(sevc0_virt(npwx,(nbnd_total-nbnd),nks))
   allocate(F(nbnd,(nbnd_total-nbnd),n_ipol))
   allocate(R(nbnd,(nbnd_total-nbnd),n_ipol))
   allocate(chi(3,3))
   chi(:,:)=cmplx(0.0d0,0.0d0,dp)
   F(:,:,:)=cmplx(0.0d0,0.0d0,dp)
   R(:,:,:)=cmplx(0.0d0,0.0d0,dp)
  endif
  !
  allocate(evc1_old(npwx,nbnd,nks,2))
  allocate(evc1(npwx,nbnd,nks,2))
  allocate(evc1_new(npwx,nbnd,nks,2))
  allocate(sevc1_new(npwx,nbnd,nks,2))
  allocate(d0psi(npwx,nbnd,nks,n_ipol))
  !
  allocate(revc0(nrxx,nbnd,nks))
  !
  allocate(rho_1(nrxx,nspin_mag))
  rho_1(:,:)=0.0d0
  !allocate(rho_tot(nrxx))
  if (charge_response == 2 ) then 
   !allocate(rho_1_tot(nrxx,nspin_mag)) !Due to broadening this is now done in lr_charg_resp
   !rho_1_tot(:,:)=0.0d0
   !print *,"allocating beta w_t"
   allocate(w_T_beta_store(itermax_int))
   allocate(w_T_gamma_store(itermax_int))
   allocate(w_T_zeta_store(w_T_npol,itermax_int))
   w_T_gamma_store(:)=0.0d0
   w_T_beta_store(:)=0.0d0
   w_T_zeta_store(:,:)=cmplx(0.0d0,0.0d0,dp)

  endif
  if (charge_response /=0) then
   allocate(w_T(itermax_int))
  endif

  allocate(dmuxc ( nrxx , nspin , nspin))
  !print *, "dmuxc ALLOCATED",allocated(dmuxc)," SIZE=",size(dmuxc)
  !print *, "nks=",nks
  
  !allocate (nbnd_occ (nks))
  !
  
  evc0(:,:,:)=(0.0d0,0.0d0)
  evc1_old(:,:,:,:)=(0.0d0,0.0d0)
  evc1(:,:,:,:)=(0.0d0,0.0d0)
  evc1_new(:,:,:,:)=(0.0d0,0.0d0)
  sevc1_new(:,:,:,:)=(0.0d0,0.0d0)
  !rho_tot(:)=0.0d0
  d0psi(:,:,:,:)=(0.0d0,0.0d0)
  !
  allocate(alpha_store(n_ipol,itermax))
  allocate(beta_store(n_ipol,itermax))
  allocate(gamma_store(n_ipol,itermax))
  allocate(zeta_store(n_ipol,n_ipol,itermax))
  alpha_store(:,:)=0.0d0
  beta_store(:,:)=0.0d0
  gamma_store(:,:)=0.0d0
  zeta_store(:,:,:)=(0.0d0,0.0d0)
  !
  if(gamma_only) then
     call lr_alloc_init_gamma()
  else
     call lr_alloc_init_k()
  endif
  !
  return
  !
contains
  !
  subroutine lr_alloc_init_gamma()
    !
    use becmod,               only : allocate_bec_type, bec_type, becp
    !
    if (nkb > 0) then
       !allocate(rbecp(nkb,nbnd))
       if (.not. allocated(becp%r)) call allocate_bec_type(nkb,nbnd,becp)
       becp%r(:,:)=0.0d0
       allocate(becp1(nkb,nbnd))
       becp1(:,:)=0.0d0
       if (project) then
        allocate(becp1_virt(nkb,nbnd_total-nbnd))
        becp1_virt(:,:)=0.0d0
       endif
    endif
    !
    return
  end subroutine lr_alloc_init_gamma
  !
  subroutine lr_alloc_init_k()
    use becmod,               only : allocate_bec_type, bec_type, becp
    !
    if (nkb > 0) then
       if(.not. allocated(becp%k)) call allocate_bec_type(nkb,nbnd,becp) 
       becp%k(:,:)=(0.0d0,0.0d0)
       allocate(becp1_c(nkb,nbnd,nks))
       becp1_c(:,:,:)=(0.0d0,0.0d0)  
       if (project) then
        allocate(becp1_c_virt(nkb,nbnd_total-nbnd,nks))
        becp1_c_virt(:,:,:)=(0.0d0,0.0d0)
       endif
    endif
    !
    return
  end subroutine lr_alloc_init_k
  !
end subroutine lr_alloc_init
!----------------------------------------------------------------------------
