!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusstres (sigmanlc)
  !----------------------------------------------------------------------
  !
  !   This routine computes the part of the atomic force which is due
  !   to the dependence of the Q function on the atomic position.
  !   Adds contribution to input sigmanlc, does not sum contributions
  !   from various processors (sum is performed by calling routine)
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega, tpiba
  USE fft_base,   ONLY : dfftp
  USE gvect,      ONLY : ngm, nl, nlm, gg, g, eigts1, eigts2, eigts3, mill
  USE lsda_mod,   ONLY : nspin
  USE scf,        ONLY : v, vltot
  USE uspp,       ONLY : becsum, okvan
  USE uspp_param, ONLY : upf, lmaxq, nh, nhm
  USE control_flags, ONLY : gamma_only
  USE fft_interfaces,ONLY : fwfft
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(inout) :: sigmanlc (3, 3)
  ! the nonlocal stress
  INTEGER :: ig, nt, ih, jh, ijh, ipol, jpol, is, na, nij
  ! counters
  COMPLEX(DP), ALLOCATABLE :: aux(:), aux1(:,:), aux2(:,:), vg(:,:), qgm(:,:)
  ! work space (complex)
  COMPLEX(DP)              :: cfac
  REAL(dp)                 :: fac(3,nspin), sus(3,3)
  ! auxiliary variables
  REAL(DP) , ALLOCATABLE :: qmod(:), ylmk0(:,:), dylmk0(:,:), tbecsum(:,:)
  ! work space (real)
  !
  !
  sus(:,:) = 0.d0
  !
  ALLOCATE ( aux1(ngm,3), aux2(ngm,nspin), qmod(ngm) )
  ALLOCATE ( ylmk0(ngm,lmaxq*lmaxq), dylmk0(ngm,lmaxq*lmaxq) )
  !
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
!$omp parallel do default(shared) private(ig)
  DO ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  ENDDO
!$omp end parallel do
  !
  ! fourier transform of the total effective potential
  !
  ALLOCATE ( vg(ngm,nspin))
  ALLOCATE ( aux(dfftp%nnr) )
  DO is = 1, nspin
     IF ( nspin == 4 .and. is /= 1 ) THEN
        aux(:) = v%of_r(:,is)
     ELSE
        aux(:) = vltot(:) + v%of_r(:,is)
     ENDIF
     CALL fwfft ('Dense', aux, dfftp)
     DO ig = 1, ngm
        vg (ig, is) = aux (nl (ig) )
     ENDDO
  ENDDO
  DEALLOCATE ( aux )
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  ! (no contribution from G=0)
  !
  DO ipol = 1, 3
     CALL dylmr2 (lmaxq * lmaxq, ngm, g, gg, dylmk0, ipol)
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           nij = nh(nt)*(nh(nt)+1)/2
           ALLOCATE (qgm(ngm,nij), tbecsum(nij,nspin) )
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1
                 CALL dqvan2 (ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0, &
                      dylmk0, ipol)
              ENDDO
           ENDDO
           !
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                 !
                 tbecsum(:,:) = becsum(1:nij,na,1:nspin)
                 !
                 CALL dgemm( 'N', 'N', 2*ngm, nspin, nij, 1.0_dp, &
                      qgm, 2*ngm, tbecsum, nij, 0.0_dp, aux2, 2*ngm )
                 !
!$omp parallel do default(shared) private(is, ig)
                 DO is = 1, nspin
                    DO ig = 1, ngm
                       aux2(ig,is) = aux2(ig,is) * CONJG(vg (ig, is))
                    END DO
                 END DO
!$omp end parallel do
!$omp parallel do default(shared) private(ig, cfac)
                 DO ig = 1, ngm
                    cfac = CONJG( eigts1 (mill (1,ig), na) * &
                                  eigts2 (mill (2,ig), na) * &
                                  eigts3 (mill (3,ig), na) )
                     aux1 (ig,1) = cfac * g (1, ig)
                     aux1 (ig,2) = cfac * g (2, ig)
                     aux1 (ig,3) = cfac * g (3, ig)
                ENDDO
!$omp end parallel do
                CALL DGEMM('T','N', 3, nspin, 2*ngm, 1.0_dp, aux1, 2*ngm, &
                           aux2, 2*ngm, 0.0_dp, fac, 3 )    
                DO is = 1, nspin
                   DO jpol = 1, ipol
                      sus (ipol, jpol) = sus (ipol, jpol) - omega * &
                          fac (jpol, is)
                   ENDDO
                ENDDO
              ENDIF
           ENDDO
           DEALLOCATE ( tbecsum, qgm )
        ENDIF
     ENDDO

  ENDDO

  IF (gamma_only) THEN
     sigmanlc(:,:) = sigmanlc(:,:) + 2.0_dp*sus(:,:)
  ELSE
     sigmanlc(:,:) = sigmanlc(:,:) + sus(:,:)
  ENDIF
  DEALLOCATE (ylmk0, dylmk0)
  DEALLOCATE (aux1, aux2, vg, qmod)

  RETURN

END SUBROUTINE addusstres

