!-----------------------------------------------------------------------
subroutine lr_alloc_init()
  !---------------------------------------------------------------------
  ! ... allocates and initialises linear response variables
  !---------------------------------------------------------------------
  ! OBM
  ! 060508 gamma_only correction
  ! 150608 allocate for becp replaced by allocate_bec (in a crude manner)
#include "f_defs.h"
  !
  use gvect,             only : nrxx
  use klist,             only : nks
  use lr_variables
  use uspp,              only : nkb
  use lsda_mod,          only : nspin
  use wvfct,             only : npwx, nbnd
  use control_flags,     only : gamma_only
  USE io_global,      ONLY : stdout
  USE charg_resp,     ONLY : w_T, w_T_beta_store, w_T_gamma_store
  use realus,            only : igk_k, npw_k
  use control_ph,            ONLY : nbnd_occ
  !
  implicit none
  !
  If (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_alloc_init>")')
  endif
  !
  allocate(evc0(npwx,nbnd,nks))
  allocate(sevc0(npwx,nbnd,nks))
  !
  allocate(evc1_old(npwx,nbnd,nks,2))
  allocate(evc1(npwx,nbnd,nks,2))
  allocate(evc1_new(npwx,nbnd,nks,2))
  allocate(sevc1_new(npwx,nbnd,nks,2))
  allocate(d0psi(npwx,nbnd,nks,n_ipol))
  !
  allocate(revc0(nrxx,nbnd,nks))
  !
  allocate(rho_1(nrxx))
  !allocate(rho_tot(nrxx))
  if (charge_response == 2 ) then 
   allocate(rho_1_tot(nrxx))
   rho_1_tot(:)=0.0d0
   !print *,"allocating beta w_t"
   allocate(w_T_beta_store(itermax))
   allocate(w_T_gamma_store(itermax))
  endif
  if (charge_response /=0 .or. lr_verbosity > 0) then
   allocate(w_T(itermax))
  endif

  allocate(dmuxc ( nrxx , nspin , nspin))
  !print *, "dmuxc ALLOCATED",allocated(dmuxc)," SIZE=",size(dmuxc)
  !print *, "nks=",nks
  
  ! Open shell related
  IF ( .not. ALLOCATED( igk_k ) )    allocate(igk_k(npwx,nks))
  IF ( .not. ALLOCATED( npw_k ) )    allocate(npw_k(nks))
  !allocate (nbnd_occ (nks))
  !
  
  evc0(:,:,:)=(0.0d0,0.0d0)
  evc1_old(:,:,:,:)=(0.0d0,0.0d0)
  evc1(:,:,:,:)=(0.0d0,0.0d0)
  evc1_new(:,:,:,:)=(0.0d0,0.0d0)
  sevc1_new(:,:,:,:)=(0.0d0,0.0d0)
  rho_1(:)=0.0d0
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
    endif
    !
    return
  end subroutine lr_alloc_init_k
  !
end subroutine lr_alloc_init
!----------------------------------------------------------------------------
