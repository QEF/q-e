!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .  
!
!---------------------------------------------------------------------------
SUBROUTINE sd0psi()
    !-----------------------------------------------------------------------
    !
    ! This subroutine calculates S * d0psi for US case. 
    ! Input : d0psi - starting Lanczos vector
    ! Output: d0psi = S * d0psi
    !
    ! Modified by Osman Baris Malcioglu (2009)
    ! Modified by Iurii Timrov (EELS extension) (2013)
    !
    USE klist,                ONLY : nks, xk
    USE lr_variables,         ONLY : n_ipol, d0psi, lr_verbosity, eels
    USE uspp,                 ONLY : vkb, nkb, okvan
    USE wvfct,                ONLY : nbnd, npwx
    USE becmod,               ONLY : bec_type, becp, calbec
    USE io_global,            ONLY : stdout
    USE qpoint,               ONLY : nksq
    !
    IMPLICIT NONE
    !
    INTEGER :: ik 
    !
    IF (lr_verbosity > 5) THEN
      WRITE(stdout,'("<sd0psi>")')
    ENDIF
    !
    IF ( nkb==0 .or. (.not.okvan) ) RETURN
    !
    CALL start_clock('sd0psi')
    !
    IF (eels) THEN
       CALL lr_sd0psi_eels()
    ELSE
       CALL lr_sd0psi_optical()
    ENDIF
    !
    CALL stop_clock('sd0psi') 
    !
    RETURN
    !
CONTAINS

SUBROUTINE lr_sd0psi_optical()
    !
    ! Optical case
    !
    USE control_flags,  ONLY : gamma_only
    USE realus,         ONLY : real_space, invfft_orbital_gamma, fwfft_orbital_gamma, &
                             & calbec_rs_gamma, v_loc_psir, s_psir_gamma, &
                             & igk_k, npw_k, real_space_debug
 
    IMPLICIT NONE
    !
    INTEGER :: ip, ibnd
    !
    DO ip=1,n_ipol
       !
       IF (gamma_only) THEN
          !
          IF (real_space_debug>4) THEN
             !
             DO ibnd=1,nbnd,2
                CALL invfft_orbital_gamma(d0psi(:,:,1,ip),ibnd,nbnd)
                CALL calbec_rs_gamma(ibnd,nbnd,becp%r)
                CALL s_psir_gamma(ibnd,nbnd)
                CALL fwfft_orbital_gamma(d0psi(:,:,1,ip),ibnd,nbnd)
             ENDDO
             !
          ELSE
             !
             CALL calbec(npw_k(1),vkb,d0psi(:,:,1,ip),becp)
             CALL s_psi(npwx,npw_k(1),nbnd,d0psi(:,:,1,ip),d0psi(:,:,1,ip))
             !
          ENDIF
          !
       ELSE
         !
         DO ik = 1, nksq
            !
            CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
            CALL calbec(npw_k(ik),vkb,d0psi(:,:,ik,ip),becp)
            CALL s_psi(npwx,npw_k(ik),nbnd,d0psi(:,:,ik,ip),d0psi(:,:,ik,ip))
            !
         ENDDO
         !
       ENDIF
       !
    ENDDO
    !
    RETURN
    !
END SUBROUTINE lr_sd0psi_optical

SUBROUTINE lr_sd0psi_eels()
   !
   ! EELS
   ! Written by I. Timrov (2013)
   !
   USE lr_variables,    ONLY : lr_periodic
   USE qpoint,          ONLY : nksq, npwq, igkq, ikks, ikqs
   USE gvect,           ONLY : ngm, g
   USE wvfct,           ONLY : g2kin
   USE gvecw,           ONLY : gcutw
   USE control_ph,      ONLY : nbnd_occ

   IMPLICIT NONE
   !
   INTEGER :: ikk, ikq
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
      ! Determination of npwq, igkq; g2kin is used here as a workspace.
      !
      CALL gk_sort( xk(1,ikq), ngm, g, gcutw, npwq, igkq, g2kin )
      !
      ! Calculate beta-functions vkb at point k+q
      !
      CALL init_us_2(npwq, igkq, xk(1,ikq), vkb)
      !
      ! Calculate the product of beta-functions vkb with
      ! the response orbitals evc1 : becp%k = <vkb|evc1>
      !
      CALL calbec(npwq, vkb, d0psi(:,:,ik,1), becp, nbnd_occ(ikk))
      !
      ! Apply the S operator
      !
      CALL s_psi(npwx, npwq, nbnd_occ(ikk), d0psi(:,:,ik,1), d0psi(:,:,ik,1))
      !
   ENDDO
   !
   RETURN
   !
END SUBROUTINE lr_sd0psi_eels
    
END SUBROUTINE sd0psi
