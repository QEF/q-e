!-----------------------------------------------------------------------
subroutine lr_dealloc()
  !---------------------------------------------------------------------
  ! ... deallocates all the Lanczos variables
  !---------------------------------------------------------------------
  !
  ! OBM
  ! 050608 gamma_only correction
#include "f_defs.h"
  !
  use lr_variables
  use uspp,              only : nkb
  use control_flags,     only : gamma_only
  use realus,            only : igk_k,npw_k
  USE lr_variables,   ONLY : lr_verbosity
  USE io_global,      ONLY : stdout
  use charg_resp,     ONLY : w_T_beta_store, w_T_gamma_store, w_T

  implicit none
  !
  If (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_dealloc>")')
  endif
  !
  if (allocated(evc0)) deallocate(evc0)
  if (allocated(sevc0)) deallocate(sevc0)
  if (allocated(evc1_old)) deallocate(evc1_old)
  if (allocated(evc1)) deallocate(evc1)
  if (allocated(evc1_new)) deallocate(evc1_new)
  if (allocated(sevc1_new)) deallocate(sevc1_new)
  if (allocated(d0psi)) deallocate(d0psi)
  !
  if (allocated(rho_1)) deallocate(rho_1)
  !if (allocated(rho_tot)) deallocate(rho_tot)
  if (allocated(dmuxc)) deallocate(dmuxc)
  if (allocated(igk_k)) deallocate(igk_k)
  if (allocated(npw_k)) deallocate(npw_k)
  !
  if (allocated(eval1)) deallocate(eval1)
  if (allocated(eval2)) deallocate(eval2)
  if (allocated(vl)) deallocate(vl)
  if (allocated(vr)) deallocate(vr)
  !
  if (allocated(alpha_store)) deallocate(alpha_store)
  if (allocated(beta_store)) deallocate(beta_store) 
  if (allocated(gamma_store)) deallocate(gamma_store)
  if (allocated(zeta_store)) deallocate(zeta_store)
  !
  !Response charge density related
  !
  if (allocated(w_T_beta_store))  deallocate(w_T_beta_store)
  if (allocated(w_T_gamma_store)) deallocate(w_T_gamma_store)
  if (allocated(w_T)) deallocate(w_T)
  if (allocated(rho_1_tot)) deallocate(rho_1_tot)
  !
  if (gamma_only) then
     call lr_dealloc_gamma()
  else
     call lr_dealloc_k()
  endif
  !
  return
  !
contains
  !
  subroutine lr_dealloc_gamma()
    !
    use becmod,               only : bec_type, becp, deallocate_bec_type
    !
    if (nkb > 0) then
       call deallocate_bec_type(becp)
       deallocate(becp1)
    endif
    !
  end subroutine lr_dealloc_gamma
  !
  subroutine lr_dealloc_k()
    !
    use becmod,               only : bec_type, becp, deallocate_bec_type
    !
    if (nkb > 0) then
       call deallocate_bec_type(becp)
       deallocate(becp1_c)
    endif
    !
  end subroutine lr_dealloc_k
  !
end subroutine lr_dealloc
!-----------------------------------------------------------------------
