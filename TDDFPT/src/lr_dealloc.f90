!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------

SUBROUTINE lr_dealloc()
  !---------------------------------------------------------------------
  ! ... deallocates all the Lanczos variables
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu (2009)
  !
  USE lr_variables
  USE uspp,           ONLY : nkb
  USE control_flags,  ONLY : gamma_only
  USE realus,         ONLY : igk_k,npw_k
  USE io_global,      ONLY : stdout
  USE charg_resp,     ONLY : w_T_beta_store, w_T_gamma_store, w_T,&
                           & w_T_zeta_store, chi, rho_1_tot, rho_1_tot_im
  USE eqv,            ONLY : dmuxc
  USE lr_exx_kernel,     ONLY : lr_exx_dealloc

  IMPLICIT NONE
  !
  IF (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_dealloc>")')
  ENDIF
  !
  IF (allocated(evc0)) DEALLOCATE(evc0)
  IF (allocated(sevc0)) DEALLOCATE(sevc0)
  IF (allocated(evc1_old)) DEALLOCATE(evc1_old)
  IF (allocated(evc1)) DEALLOCATE(evc1)
  IF (allocated(evc1_new)) DEALLOCATE(evc1_new)
  IF (allocated(sevc1_new)) DEALLOCATE(sevc1_new)
  IF (allocated(d0psi)) DEALLOCATE(d0psi)
  if (allocated(tg_revc0)) DEALLOCATE( tg_revc0 )
  if (allocated(revc0)) DEALLOCATE( revc0 )

  !
  IF (project) THEN
   DEALLOCATE(evc0_virt)
   !deallocate(sevc0_virt)
   DEALLOCATE(F)
   DEALLOCATE(R)
  ENDIF


  !
  IF (allocated(rho_1)) DEALLOCATE(rho_1)
  !if (allocated(rho_tot)) deallocate(rho_tot)
  IF (allocated(dmuxc)) DEALLOCATE(dmuxc)
  IF (allocated(igk_k)) DEALLOCATE(igk_k)
  IF (allocated(npw_k)) DEALLOCATE(npw_k)
  !
  IF (allocated(eval1)) DEALLOCATE(eval1)
  IF (allocated(eval2)) DEALLOCATE(eval2)
  IF (allocated(vl)) DEALLOCATE(vl)
  IF (allocated(vr)) DEALLOCATE(vr)
  !
  IF (allocated(alpha_store)) DEALLOCATE(alpha_store)
  IF (allocated(beta_store)) DEALLOCATE(beta_store)
  IF (allocated(gamma_store)) DEALLOCATE(gamma_store)
  IF (allocated(zeta_store)) DEALLOCATE(zeta_store)
  !
  !Response charge density related
  !
  IF (allocated(w_T_beta_store))  DEALLOCATE(w_T_beta_store)
  IF (allocated(w_T_gamma_store)) DEALLOCATE(w_T_gamma_store)
  IF (allocated(w_T_zeta_store)) DEALLOCATE(w_T_zeta_store)
  IF (allocated(chi)) DEALLOCATE(chi)
  IF (allocated(w_T)) DEALLOCATE(w_T)
  IF (allocated(rho_1_tot)) DEALLOCATE(rho_1_tot)
  IF (allocated(rho_1_tot_im)) DEALLOCATE(rho_1_tot_im)
  IF (lr_exx) CALL lr_exx_dealloc
  !
  IF (gamma_only) THEN
     CALL lr_dealloc_gamma()
  ELSE
     CALL lr_dealloc_k()
  ENDIF
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_dealloc_gamma()
    !
    USE becmod,               ONLY : bec_type, becp, deallocate_bec_type
    !
    IF (nkb > 0) THEN
       CALL deallocate_bec_type(becp)
       DEALLOCATE(becp1)
       IF (project) THEN
        DEALLOCATE(becp1_virt)
       ENDIF
     ENDIF
    !
  END SUBROUTINE lr_dealloc_gamma
  !
  SUBROUTINE lr_dealloc_k()
    !
    USE becmod,               ONLY : bec_type, becp, deallocate_bec_type
    !
    IF (nkb > 0) THEN
       CALL deallocate_bec_type(becp)
       DEALLOCATE(becp1_c)
       IF (project) THEN
        DEALLOCATE(becp1_c_virt)
       ENDIF
    ENDIF
    !
  END SUBROUTINE lr_dealloc_k
  !
END SUBROUTINE lr_dealloc
!-----------------------------------------------------------------------
