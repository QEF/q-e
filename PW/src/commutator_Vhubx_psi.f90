!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE commutator_Vhubx_psi(ik, nbnd_calc, vpol, dpsi)
  !-----------------------------------------------------------------------
  !! This routine computes the commutator between the non-local
  !! Hubbard potential and the position operator, applied to \text{psi}
  !! of the current k-point. The result is added to \text{dpsi}.
  !
  ! calculates: [V_hub, r \dot vpol]|psi_nk>
  !
  ! vpol is the polarization vector in Cartesian coordinates.
  ! For crystal coordinate, use vpol = at(:, ipol).
  ! For Cartesian coordinate, use vpol = (1.0, 0.0, 0.0) or other permutations.
  !
  ! Some insights about the formulas here can be found e.g.
  ! in I. Timrov's PhD thesis, Sec. 6.1.3,
  ! https://pastel.archives-ouvertes.fr/pastel-00823758
  !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,          ONLY : DP
  USE io_files,       ONLY : iunhub, iunhub_noS, nwordwfcU
  USE wavefunctions,  ONLY : evc
  USE wvfct,          ONLY : npwx
  USE ions_base,      ONLY : nat, ityp, ntyp => nsp
  USE ldaU,           ONLY : Hubbard_l, Hubbard_U, Hubbard_J0, &
                             is_hubbard, nwfcU, offsetU, oatwfc
  USE uspp,           ONLY : vkb, nkb, okvan
  USE uspp_param,     ONLY : nh, upf
  USE uspp_init,      ONLY : gen_us_dj, gen_us_dy
  USE lsda_mod,       ONLY : lsda, current_spin, isk, nspin
  USE klist,          ONLY : xk, ngk, igk_k
  USE cell_base,      ONLY : tpiba
  USE gvect,          ONLY : g
  USE scf,            ONLY : rho
  USE mp,             ONLY : mp_sum
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE buffers,        ONLY : get_buffer
  USE basis,          ONLY : natomwfc
  USE noncollin_module, ONLY : noncolin, npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k-point index
  INTEGER, INTENT(IN) :: nbnd_calc
  !! Number of bands to calculate [V_hub, x_ipol]|psi_ik>
  REAL(DP), INTENT(IN) :: vpol(3)
  !! Polarization vector in Cartesian coordinates
  COMPLEX(DP), INTENT(OUT) :: dpsi(npwx*npol, nbnd_calc)
  !! Output wavefunction where [V_hub, x_ipol]|psi_ik> is added
  !
  REAL(DP), PARAMETER :: eps = 1.0d-8
  INTEGER     :: na, n ,l, nt, nah, ikb , m, m1, m2, ibnd, ib, ig, jkb, i, &
                 ihubst, ihubst1,  ihubst2, icart, op_spin, npw, offpm, offpmU
  REAL(DP)    :: nsaux, g2k, gk_ig(3)
  REAL(DP), ALLOCATABLE :: xyz(:,:)
  REAL(DP), ALLOCATABLE :: gk_vpol(:)
  !! ipol component of the G/|G| vector, in crystal coordinates
  COMPLEX(DP), ALLOCATABLE :: dkwfcbessel(:,:), dkwfcylmr(:,:), dkwfcatomk(:,:),   &
                 dpqq26(:,:), dpqq38(:,:), dpqq47(:,:), dkvkbbessel(:,:),          &
                 dkvkbylmr(:,:), dkvkb(:,:), aux_1234(:), termi(:,:), trm(:,:),    &
                 wfcatomk(:,:), swfcatomk(:,:), proj1(:,:), proj2(:,:), proj3(:,:)
  CALL start_clock( 'commutator_Vhubx_psi' )
  !
  ! Number of plane waves at point ik
  npw = ngk(ik)
  !
  ALLOCATE (proj1(nbnd_calc,nwfcU))
  ALLOCATE (proj2(nbnd_calc,nwfcU))
  ALLOCATE (proj3(nbnd_calc,nwfcU))
  ALLOCATE (dkwfcbessel(npwx,natomwfc))
  ALLOCATE (dkwfcylmr(npwx,natomwfc))
  ALLOCATE (dkwfcatomk(npwx,nwfcU))
  ALLOCATE (dpqq26(npwx,nwfcU))
  ALLOCATE (dpqq38(npwx,nwfcU))
  ALLOCATE (dpqq47(npwx,nwfcU))
  ALLOCATE (wfcatomk(npwx,nwfcU))
  ALLOCATE (swfcatomk(npwx,nwfcU))
  ALLOCATE (dkvkbbessel(npwx,nkb))
  ALLOCATE (dkvkbylmr(npwx,nkb))
  ALLOCATE (dkvkb(npwx,nkb))
  ALLOCATE (aux_1234(npwx))
  ALLOCATE (termi(npwx,nbnd_calc))
  ALLOCATE (trm(npwx,nbnd_calc))
  ALLOCATE (xyz(3,3))
  ALLOCATE (gk_vpol(npw))
  !
  dpqq26     = (0.d0, 0.d0)
  dpqq38     = (0.d0, 0.d0)
  dpqq47     = (0.d0, 0.d0)
  dkwfcatomk = (0.d0, 0.d0)
  dkvkb      = (0.d0, 0.d0)
  !
  IF (lsda) THEN
    current_spin = isk(ik)
    if (nspin == 2) then
       if (current_spin == 1) then
          op_spin = 2
       else
          op_spin = 1
       endif
    else
       op_spin=1
    endif
  ENDIF
  !
  ! Read the atomic orbitals \phi at k from file (unit iunhub_noS)
  !
  CALL get_buffer (wfcatomk, nwordwfcU, iunhub_noS, ik)
  !
  ! Read S*\phi at k from file (unit iunhub)
  !
  CALL get_buffer (swfcatomk, nwordwfcU, iunhub, ik)
  !
  ! Derivatives w.r.t. k of the atomic wfc
  ! \phi'_(k+G,I,m)_ipol> = exp^-i(k+G)*tau_I * d/dk_ipol[\phi_0(k+G,I,m)]
  ! where \phi_0(k+G,I,m) is the Fourier component of the atomic wfc localized in zero
  !
  ! Derivative of Bessel functions and spherical harmonics wrt to crystal axis ipol
  !
  CALL gen_at_dj (ik, dkwfcbessel)
  CALL gen_at_dy (ik, vpol, dkwfcylmr)
  !
  DO ig = 1, npw
     !
     ! The gk factor is necessary because we do not want the derivative of the Bessel functions
     ! w.r.t. the modulus (calculated in gen_at_dj.f90), but w.r.t.
     ! the cartesian component and then to crystal component ipol
     ! gk_icart= d|k+G|/d(k+G)_icart
     !
     gk_ig = (xk(:,ik) + g(:,igk_k(ig,ik))) * tpiba
     g2k = SUM( gk_ig**2 )
     !
     ! Take the component along the vpol vector
     gk_vpol(ig) = SUM(vpol(:) * gk_ig(:))
     !
     IF (g2k < 1.0d-10) THEN
        gk_vpol(ig) = 0.d0
     ELSE
        gk_vpol(ig) = gk_vpol(ig) / SQRT(g2k)
     ENDIF
     !
  ENDDO
  !
  DO ig = 1, npw
     !
     ! Derivative wrt crystal axis ipol
     ! d|k+G|/d(k+G)_ipol = \sum_{icart} d|k+G|/d(k+G)_icart * at (icart,ipol)
     ! The derivative is done for all the atomic wfc and for each ig
     !
     DO na = 1, nat
         nt = ityp(na)
         IF (is_hubbard(nt)) THEN
            offpmU = offsetU(na)
            offpm  = oatwfc(na)
            DO m1 = 1, 2*Hubbard_l(nt)+1
               dkwfcatomk(ig,offpmU+m1) = dkwfcylmr(ig,offpm+m1)        &
                                        + dkwfcbessel(ig,offpm+m1) * gk_vpol(ig)
            ENDDO
         ENDIF
     ENDDO
     !
  ENDDO
  !
  ! USPP case
  !
  IF (okvan)  THEN
     !
     ! Compute the k derivatives of the beta functions.
     ! Here the derivative of the Bessel functions and the spherical harmonics
     !
     CALL gen_us_dj (ik, dkvkbbessel)
     CALL gen_us_dy (ik, vpol, dkvkbylmr)
     !
     jkb = 0
     DO nt = 1, ntyp
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              DO ikb = 1, nh(nt)
                 jkb = jkb + 1
                 DO ig = 1, npw
                    dkvkb(ig,jkb) = dkvkbylmr(ig,jkb) + dkvkbbessel(ig,jkb) * gk_vpol(ig)
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     !
  ENDIF
  !
  ! Preliminary calculation of various scalar products
  ! The Hubbard terms
  !
  DO nah = 1, nat
     !
     nt = ityp (nah)
     !
     IF (is_hubbard(nt)) THEN
        !
        DO m = 1, 2 * Hubbard_l(nt) + 1
           !
           ihubst = offsetU(nah) + m
           !
           IF (okvan) THEN
              !
              ! vecqqproj for terms 2,3,4,6,7,8
              ! term 6 is the cc of term 2 with m <=> m'
              ! the same holds for 3 and 8, 4 and 7
              ! Note: these are the notations from private notes of A. Floris
              !
              ! FIXME: compute all dpqqNN in vecqqproj, not just one
              ! FIXME: replace dot_product calls with matrix-matrix products
              !
              CALL vecqqproj (npw, vkb, vkb, dkwfcatomk(:,ihubst), dpqq26(:,ihubst))
              CALL vecqqproj (npw, dkvkb, vkb, wfcatomk(:,ihubst), dpqq38(:,ihubst))
              CALL vecqqproj (npw, vkb, dkvkb, wfcatomk(:,ihubst), dpqq47(:,ihubst))
              DO ibnd = 1, nbnd_calc
                 proj3(ibnd,ihubst) = &
                      dot_product (dpqq26(1:npw,ihubst), evc(1:npw,ibnd)) + &
                      dot_product (dpqq47(1:npw,ihubst), evc(1:npw,ibnd)) + &
                      dot_product (dpqq38(1:npw,ihubst), evc(1:npw,ibnd))
              ENDDO
              !
           ENDIF
           !
           DO ibnd = 1, nbnd_calc
              !
              ! Calculate proj (ihubst,ibnd) = < S_{k}\phi_(k,I,m)| psi(ibnd,ik) >
              ! at ihubst (i.e. I, m).
              !
              proj1(ibnd,ihubst) = dot_product (swfcatomk(1:npw,ihubst), evc(1:npw,ibnd))
              proj2(ibnd,ihubst) = dot_product (dkwfcatomk(1:npw,ihubst), evc(1:npw,ibnd))
              !
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_sum(proj1, intra_bgrp_comm)
  CALL mp_sum(proj2, intra_bgrp_comm)
  CALL mp_sum(proj3, intra_bgrp_comm)
  !
  DO nah = 1, nat   ! the Hubbard atom
     !
     nt = ityp (nah)
     !
     IF (is_hubbard(nt)) THEN
            !
            termi = (0.d0, 0.d0)
            !
            DO m1 = 1, 2*Hubbard_l(nt)+1
               !
               ihubst1 = offsetU(nah) + m1
               aux_1234 = (0.d0, 0.d0)
               !
               ! term1 + term2 + term3 + term4
               !
               aux_1234 =  dkwfcatomk(:,ihubst1)
               !
               IF (okvan) THEN
                  aux_1234 = aux_1234 + dpqq26(:,ihubst1) &
                                      + dpqq38(:,ihubst1) &
                                      + dpqq47(:,ihubst1)
               ENDIF
               !
               DO m2 = 1, 2 * Hubbard_l(nt) + 1
                  !
                  ihubst2 = offsetU(nah) + m2
                  !
                  trm = (0.d0, 0.d0)
                  !
                  nsaux = rho%ns(m1, m2, current_spin, nah)
                  !
                  DO ibnd = 1, nbnd_calc
                     trm(:,ibnd) = aux_1234(:) * proj1(ibnd,ihubst2)  + &     ! term_1234
                                   swfcatomk(:,ihubst1) * proj2(ibnd,ihubst2) ! term 5
                  ENDDO
                  !
                  IF (okvan) THEN
                     DO ibnd = 1, nbnd_calc
                        trm(:,ibnd) = trm(:,ibnd) + swfcatomk(:,ihubst1) * &
                                                    proj3(ibnd,ihubst2) ! term_6+7+8
                     ENDDO
                  ENDIF
                  !
                  ! termi (npwx,nbnd), trm (npwx,nbnd), summing for all bands and G vectors
                  !
                  DO ibnd = 1, nbnd_calc
                     IF (m1 == m2) termi(:,ibnd) = termi(:,ibnd) + 0.5d0 * trm(:,ibnd)
                     termi(:,ibnd) = termi(:,ibnd) - nsaux * trm(:,ibnd)
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
            ! We want to have -i d/dk
            !
            DO ibnd = 1, nbnd_calc
               dpsi(:,ibnd) = dpsi(:,ibnd) + (0.d0,-1.d0) * termi(:,ibnd) * &
                                             (Hubbard_U(nt) - Hubbard_J0(nt))
            ENDDO
            !
            !---------------------------------------------------------------
            ! The following is for the J0\=0 case
            !
            IF (nspin.EQ.2 .AND. Hubbard_J0(nt).NE.0.d0) THEN
               !
               termi = (0.d0, 0.d0)
               !
               DO m1 = 1, 2*Hubbard_l(nt)+1
                  !
                  ihubst1 = offsetU(nah) + m1
                  aux_1234 = (0.d0, 0.d0)
                  !
                  ! term1 + term2 + term3 + term4
                  !
                  aux_1234 = dkwfcatomk(:,ihubst1)
                  !
                  IF (okvan) THEN
                     aux_1234 = aux_1234 + dpqq26(:,ihubst1) &
                                         + dpqq38(:,ihubst1) &
                                         + dpqq47(:,ihubst1)
                  ENDIF
                  !
                  DO m2 = 1, 2*Hubbard_l(nt)+1
                     !
                     ihubst2 = offsetU(nah) + m2
                     !
                     trm = (0.d0, 0.d0)
                     !
                     nsaux = rho%ns(m1, m2, op_spin, nah)
                     !
                     DO ibnd = 1, nbnd_calc
                        trm(:,ibnd) = aux_1234(:) * proj1(ibnd,ihubst2)  + & ! term_1234
                                      swfcatomk(:,ihubst1) * proj2(ibnd,ihubst2) ! term 5
                     ENDDO
                     !
                     IF (okvan) THEN
                        DO ibnd = 1, nbnd_calc
                           trm(:,ibnd) = trm(:,ibnd) + swfcatomk(:,ihubst1) &
                                                     * proj3(ibnd,ihubst2)  ! term_6+7+8
                        ENDDO
                     ENDIF
                     !
                     ! termi (npwx,nbnd), trm (npwx,nbnd), summing for all bands and G vectors
                     !
                     DO ibnd = 1, nbnd_calc
                        termi(:,ibnd) = termi(:,ibnd) + nsaux * trm(:,ibnd) ! sign change
                     ENDDO
                     !
                  ENDDO
                  !
               ENDDO
               !
               ! We want to have -i d/dk
               !
               DO ibnd = 1, nbnd_calc
                  dpsi(:,ibnd) = dpsi(:,ibnd) + (0.d0,-1.d0) * termi(:,ibnd) * Hubbard_J0(nt)
               ENDDO
               !
            ENDIF
            !
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)
  DEALLOCATE (proj3)
  DEALLOCATE (dkwfcbessel)
  DEALLOCATE (dkwfcylmr)
  DEALLOCATE (dkwfcatomk)
  DEALLOCATE (dpqq26)
  DEALLOCATE (dpqq38)
  DEALLOCATE (dpqq47)
  DEALLOCATE (wfcatomk)
  DEALLOCATE (swfcatomk)
  DEALLOCATE (dkvkbbessel)
  DEALLOCATE (dkvkbylmr)
  DEALLOCATE (dkvkb)
  DEALLOCATE (aux_1234)
  DEALLOCATE (termi)
  DEALLOCATE (trm)
  DEALLOCATE (xyz)
  DEALLOCATE (gk_vpol)
  !
  CALL stop_clock ('commutator_Vhubx_psi')
  !
  RETURN
  !
END SUBROUTINE commutator_Vhubx_psi

SUBROUTINE vecqqproj (npw, vec1, vec2, vec3, dpqq)
    !
    ! Calculate dpqq (ig) = \sum {na l1 l2} vec1(ig ,na,l1)
    !                * qq(na, l1 ,l2) * < vec2 (na,l2) | vec3 >
    !
    USE kinds,      ONLY : DP
    USE uspp_param, ONLY : nh
    USE ions_base,  ONLY : nat, ityp
    USE uspp,       ONLY : qq_nt, nkb, ofsbeta
    USE wvfct,      ONLY : npwx
    USE mp,         ONLY : mp_sum
    USE mp_bands,   ONLY : intra_bgrp_comm
    !
    IMPLICIT NONE
    !
    ! Index of the displaced atom
    !
    INTEGER, INTENT(IN)      :: npw
    COMPLEX(DP), INTENT(IN)  :: vec1(npwx,nkb)
    COMPLEX(DP), INTENT(IN)  :: vec2(npwx,nkb)
    COMPLEX(DP), INTENT(IN)  :: vec3(npwx)
    COMPLEX(DP), INTENT(OUT) :: dpqq(npwx)
    !
    ! Local variables
    !
    INTEGER     :: na, nt, l1, l2, ig, ibeta1, ibeta2, ibnd
    COMPLEX(DP), ALLOCATABLE :: aux1(:)
    COMPLEX(DP) :: projaux1vec3
    !
    dpqq = (0.d0, 0.d0)
    !
    ALLOCATE (aux1(npwx))
    !
    DO na = 1, nat
       !
       nt = ityp(na)
       !
       DO l1 = 1, nh(nt)
          !
          ibeta1 = ofsbeta(na) + l1
          !
          !  aux1 = \sum_l2 qq_nt(l1,l2,nt) * |vec2(na,l2)>
          !
          aux1 = (0.d0, 0.d0)
          !
          DO l2 = 1, nh(nt)
             ibeta2 = ofsbeta(na) + l2
             aux1(:) = aux1(:) + qq_nt(l1,l2,nt) * vec2(:,ibeta2)
          ENDDO
          !
          projaux1vec3 = dot_product (aux1(1:npw), vec3(1:npw))
          !
          CALL mp_sum(projaux1vec3, intra_bgrp_comm)
          !
          ! Summing on na and l1 for each ig
          !
          dpqq(:) = dpqq(:) + vec1(:,ibeta1) *  projaux1vec3
          !
       ENDDO
       !
    ENDDO
    !
    DEALLOCATE (aux1)
    !
    RETURN
    !
END SUBROUTINE vecqqproj

