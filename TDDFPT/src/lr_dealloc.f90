!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_dealloc()
  !---------------------------------------------------------------------
  !
  ! This subroutine deallocates all the Lanczos/Davidson variables.
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  !
  USE lr_variables
  USE uspp,           ONLY : nkb
  USE control_flags,  ONLY : gamma_only
  USE realus,         ONLY : tg_psic
  USE klist,          ONLY : ngk, igk_k
  USE io_global,      ONLY : stdout
  USE charg_resp,     ONLY : w_T_beta_store, w_T_gamma_store, w_T,&
                           & w_T_zeta_store, chi, rho_1_tot, rho_1_tot_im
  USE lr_exx_kernel,  ONLY : lr_exx_dealloc
  USE becmod,         ONLY : bec_type, becp, deallocate_bec_type
  USE lrus,           ONLY : int3, int3_nc, becp1, &
                           & bbg, bbk, bbnc
  USE qpoint,         ONLY : ikks, ikqs, eigqts
  USE eqv,            ONLY : dmuxc, evq, dpsi, dvpsi
  USE control_lr,     ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_dealloc>")')
  ENDIF
  !
  IF (allocated(evc0))      DEALLOCATE(evc0)
  IF (allocated(sevc0))     DEALLOCATE(sevc0)
  IF (allocated(evc1_old))  DEALLOCATE(evc1_old)
  IF (allocated(evc1))      DEALLOCATE(evc1)
  IF (allocated(evc1_new))  DEALLOCATE(evc1_new)
  IF (allocated(sevc1_new)) DEALLOCATE(sevc1_new)
  IF (allocated(sevc1))     DEALLOCATE(sevc1)
  IF (allocated(d0psi))     DEALLOCATE(d0psi)
  IF (allocated(d0psi2))    DEALLOCATE(d0psi2)
  IF (allocated(tg_revc0))  DEALLOCATE(tg_revc0)
  IF (allocated(tg_psic))   DEALLOCATE(tg_psic)
  IF (allocated(revc0))     DEALLOCATE(revc0)
  IF (allocated(bbg))       DEALLOCATE(bbg)
  IF (allocated(bbk))       DEALLOCATE(bbk)
  IF (allocated(bbnc))      DEALLOCATE(bbnc)
  !
  IF (project) THEN
     DEALLOCATE(F)
     DEALLOCATE(R)
  ENDIF
  !
  IF (allocated(rho_1))  DEALLOCATE(rho_1)
  IF (allocated(rho_1c)) DEALLOCATE(rho_1c)
  IF (allocated(dmuxc))  DEALLOCATE(dmuxc)
  IF (allocated(igk_k))  DEALLOCATE(igk_k)
  IF (allocated(ngk))    DEALLOCATE(ngk)
  IF (allocated(ikks))   DEALLOCATE(ikks)
  IF (allocated(ikqs))   DEALLOCATE(ikqs)
  !
  ! EELS-related variables
  !
  IF (allocated(dpsi))    DEALLOCATE(dpsi)
  IF (allocated(dvpsi))   DEALLOCATE(dvpsi)
  IF (allocated(eigqts))  DEALLOCATE(eigqts)
  IF (allocated(int3))    DEALLOCATE(int3)
  IF (allocated(int3_nc)) DEALLOCATE(int3_nc)
  !
  IF (eels) THEN
     IF (associated(evq))    DEALLOCATE(evq)
  ELSE
     IF (associated(evq))    NULLIFY(evq)
  ENDIF 
  !
  IF (allocated(becp1)) THEN
     DO ik = 1,size(becp1)
        CALL deallocate_bec_type ( becp1(ik) )
     ENDDO
     DEALLOCATE(becp1)
  ENDIF
  !
  IF (ALLOCATED(nbnd_occ)) DEALLOCATE(nbnd_occ)
  !
  IF (allocated(alpha_store)) DEALLOCATE(alpha_store)
  IF (allocated(beta_store))  DEALLOCATE(beta_store)
  IF (allocated(gamma_store)) DEALLOCATE(gamma_store)
  IF (allocated(zeta_store))  DEALLOCATE(zeta_store)
  !
  ! Response charge density related
  !
  IF (allocated(w_T_beta_store))  DEALLOCATE(w_T_beta_store)
  IF (allocated(w_T_gamma_store)) DEALLOCATE(w_T_gamma_store)
  IF (allocated(w_T_zeta_store))  DEALLOCATE(w_T_zeta_store)
  IF (allocated(chi))             DEALLOCATE(chi)
  IF (allocated(evc0_virt))       DEALLOCATE(evc0_virt)
  IF (allocated(w_T))             DEALLOCATE(w_T)
  IF (allocated(rho_1_tot))       DEALLOCATE(rho_1_tot)
  IF (allocated(rho_1_tot_im))    DEALLOCATE(rho_1_tot_im)
  IF (allocated(cube_save))       DEALLOCATE(cube_save)
  !
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
    IF (nkb > 0) THEN
       !
       CALL deallocate_bec_type(becp)
       !
       DEALLOCATE(becp_1)
       !
       IF (project .or. davidson) THEN
          DEALLOCATE(becp1_virt)
       ENDIF
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE lr_dealloc_gamma
  !
  SUBROUTINE lr_dealloc_k()
    !
    IF (nkb > 0) THEN
       !
       CALL deallocate_bec_type(becp)
       !
       IF (.NOT.eels) DEALLOCATE(becp1_c)
       !
       IF (project .or. davidson) THEN
          DEALLOCATE(becp1_c_virt)
       ENDIF
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE lr_dealloc_k
  !
END SUBROUTINE lr_dealloc
!-----------------------------------------------------------------------
