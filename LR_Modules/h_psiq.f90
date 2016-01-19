!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE h_psiq (lda, n, m, psi, hpsi, spsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the product of the Hamiltonian
  !     and of the S matrix with a m  wavefunctions  contained
  !     in psi. It first computes the bec matrix of these
  !     wavefunctions and then with the routines hus_1psi and
  !     s_psi computes for each band the required products
  !
  !     Merged with lr_h_psiq June 2011. This function is now used both
  !     in ph.x and turbo_lanczos.x
  !

  USE kinds,  ONLY : DP
  USE wavefunctions_module,  ONLY : psic, psic_nc
  USE becmod, ONLY : bec_type, becp, calbec
  USE noncollin_module, ONLY : noncolin, npol
  USE lsda_mod, ONLY : current_spin
  USE fft_base, ONLY : dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvecs, ONLY: nls
  USE spin_orb, ONLY : domag
  USE scf,    ONLY : vrs
  USE uspp,   ONLY : vkb
  USE wvfct,  ONLY : g2kin, npwx
  USE qpoint, ONLY: igkq
  USE control_flags, ONLY : gamma_only ! Needed only for TDDFPT

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: lda, n, m
  ! input: the leading dimension of the array psi
  ! input: the real dimension of psi
  ! input: the number of psi to compute
  COMPLEX(DP), INTENT(INOUT)  :: psi (lda*npol, m)
  COMPLEX(DP), INTENT(OUT) :: hpsi (lda*npol, m), spsi (lda*npol, m)
  ! input: the functions where to apply H and S
  ! output: H times psi
  ! output: S times psi (Us PP's only)

  !
  !     Here the local variables
  !
  COMPLEX(DP) :: sup, sdwn
  INTEGER :: ibnd
  ! counter on bands
  INTEGER :: j
  ! do loop index


  CALL start_clock ('h_psiq')

  IF (gamma_only) THEN

   CALL h_psiq_gamma()

  ELSE

   CALL h_psiq_k()

  ENDIF

  CALL stop_clock ('h_psiq')
  RETURN

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!k point part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE h_psiq_k()

    IMPLICIT NONE

  CALL start_clock ('init')

  CALL calbec ( n, vkb, psi, becp, m)
  !
  ! Here we apply the kinetic energy (k+G)^2 psi
  !
  hpsi=(0.0_dp,0.0_dp)
  DO ibnd = 1, m
     DO j = 1, n
        hpsi (j, ibnd) = g2kin (j) * psi (j, ibnd)
     ENDDO
  ENDDO
  IF (noncolin) THEN
     DO ibnd = 1, m
        DO j = 1, n
           hpsi (j+lda, ibnd) = g2kin (j) * psi (j+lda, ibnd)
        ENDDO
     ENDDO
  ENDIF
  CALL stop_clock ('init')
  !
  ! the local potential V_Loc psi. First the psi in real space
  !

  DO ibnd = 1, m
     CALL start_clock ('firstfft')
     IF (noncolin) THEN
        psic_nc = (0.0_dp, 0.0_dp)
        DO j = 1, n
           psic_nc(nls(igkq(j)),1) = psi (j, ibnd)
           psic_nc(nls(igkq(j)),2) = psi (j+lda, ibnd)
        ENDDO
        CALL invfft ('Wave', psic_nc(:,1), dffts)
        CALL invfft ('Wave', psic_nc(:,2), dffts)
     ELSE
        psic(:) = (0.0_dp, 0.0_dp)
        DO j = 1, n
           psic (nls(igkq(j))) = psi (j, ibnd)
        ENDDO
        CALL invfft ('Wave', psic, dffts)
     ENDIF
     CALL stop_clock ('firstfft')
     !
     !   and then the product with the potential vrs = (vltot+vr) on the smoo
     !
     IF (noncolin) THEN
        IF (domag) THEN
           DO j=1, dffts%nnr
              sup = psic_nc(j,1) * (vrs(j,1)+vrs(j,4)) + &
                    psic_nc(j,2) * (vrs(j,2)-(0.0_dp,1.0_dp)*vrs(j,3))
              sdwn = psic_nc(j,2) * (vrs(j,1)-vrs(j,4)) + &
                    psic_nc(j,1) * (vrs(j,2)+(0.0_dp,1.0_dp)*vrs(j,3))
              psic_nc(j,1)=sup
              psic_nc(j,2)=sdwn
           ENDDO
        ELSE
           DO j=1, dffts%nnr
              psic_nc(j,1)=psic_nc(j,1) * vrs(j,1)
              psic_nc(j,2)=psic_nc(j,2) * vrs(j,1)
           ENDDO
        ENDIF
     ELSE
        DO j = 1, dffts%nnr
           psic (j) = psic (j) * vrs (j, current_spin)
        ENDDO
     ENDIF
     !
     !   back to reciprocal space
     !
     CALL start_clock ('secondfft')
     IF (noncolin) THEN
        CALL fwfft ('Wave', psic_nc(:,1), dffts)
        CALL fwfft ('Wave', psic_nc(:,2), dffts)
     !
     !   addition to the total product
     !
        DO j = 1, n
           hpsi (j, ibnd) = hpsi (j, ibnd) + psic_nc (nls(igkq(j)), 1)
           hpsi (j+lda, ibnd) = hpsi (j+lda, ibnd) + psic_nc (nls(igkq(j)), 2)
        ENDDO
     ELSE
        CALL fwfft ('Wave', psic, dffts)
     !
     !   addition to the total product
     !
        DO j = 1, n
           hpsi (j, ibnd) = hpsi (j, ibnd) + psic (nls(igkq(j)))
        ENDDO
     ENDIF
     CALL stop_clock ('secondfft')
  ENDDO
  !
  !  Here the product with the non local potential V_NL psi
  !

  CALL add_vuspsi (lda, n, m, hpsi)

  CALL s_psi (lda, n, m, psi, spsi)

END SUBROUTINE h_psiq_k


SUBROUTINE h_psiq_gamma()

    USE becmod, ONLY : becp, calbec
    USE gvect,  ONLY : gstart
    USE realus, ONLY : real_space, invfft_orbital_gamma, &
                       fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                       v_loc_psir, s_psir_gamma, real_space_debug
    USE uspp,                  ONLY : nkb

    IMPLICIT NONE

    CALL start_clock ('init')
    !
    ! Here we apply the kinetic energy (k+G)^2 psi
    !
    IF(gstart==2) psi(1,:)=cmplx(real(psi(1,:),dp),0.0d0,dp)
    !
    DO ibnd=1,m

       DO j=1,n

          hpsi(j,ibnd)=g2kin(j)*psi(j,ibnd)

       ENDDO

    ENDDO

    CALL stop_clock ('init')
      IF (nkb > 0 .and. real_space_debug>2) THEN
        DO ibnd=1,m,2

           CALL invfft_orbital_gamma(psi,ibnd,m,.true.)
          !transform the psi real space, saved in temporary memory

           CALL calbec_rs_gamma(ibnd,m,becp%r)
           !rbecp on psi

           CALL s_psir_gamma(ibnd,m)
           !psi -> spsi

           CALL fwfft_orbital_gamma(spsi,ibnd,m)
           !return back to real space

           CALL invfft_orbital_gamma(hpsi,ibnd,m)
           ! spsi above is now replaced by hpsi

           CALL v_loc_psir(ibnd,m)
           ! hpsi -> hpsi + psi*vrs  (psi read from temporary memory)

           CALL add_vuspsir_gamma(ibnd,m)
           ! hpsi -> hpsi + vusp

           CALL fwfft_orbital_gamma(hpsi,ibnd,m,.true.)
           !transform back hpsi, clear psi in temporary memory

        ENDDO

     ELSE

    CALL vloc_psi_gamma(lda,n,m,psi,vrs(1,current_spin),hpsi)

    IF (noncolin) THEN

        CALL errore ("h_psiq","gamma and noncolin not implemented yet",1)

    ELSE

        CALL calbec ( n, vkb, psi, becp, m)

     ENDIF

     CALL add_vuspsi (lda, n, m, hpsi)
     CALL s_psi (lda, n, m, psi, spsi)

  ENDIF

  END SUBROUTINE h_psiq_gamma

END SUBROUTINE h_psiq
