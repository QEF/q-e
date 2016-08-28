!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE addusdens(rho)
  !----------------------------------------------------------------------
  !
  USE realus,               ONLY : addusdens_r
  USE control_flags,        ONLY : tqr
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  !
  REAL(kind=dp), INTENT(inout) :: rho(dfftp%nnr,nspin_mag)
  !
  IF ( tqr ) THEN
     CALL addusdens_r(rho,.true.)
  ELSE
     CALL addusdens_g(rho)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE addusdens
!
!----------------------------------------------------------------------
SUBROUTINE addusdens_g(rho)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, nl, nlm, gg, g, &
                                   eigts1, eigts2, eigts3, mill
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : becsum, okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE
  !
  REAL(kind=dp), INTENT(inout) :: rho(dfftp%nnr,nspin_mag)
  !
  !     here the local variables
  !

  INTEGER :: ig, na, nt, ih, jh, ijh, is, nab, nb, nij
  ! counters

  REAL(DP), ALLOCATABLE :: tbecsum(:,:,:)
  ! \sum_kv <\psi_kv|\beta_l><beta_m|\psi_kv> for each species of atoms
  REAL(DP), ALLOCATABLE :: qmod (:), ylmk0 (:,:)
  ! modulus of G, spherical harmonics
  COMPLEX(DP), ALLOCATABLE :: skk(:,:), aux2(:,:)
  ! structure factors, US contribution to rho
  COMPLEX(DP), ALLOCATABLE ::  aux (:,:), qgm(:)
  ! work space for rho(G,nspin), Fourier transform of q

  IF (.not.okvan) RETURN

  CALL start_clock ('addusdens')

  ALLOCATE (aux ( ngm, nspin_mag) )
  ALLOCATE (qmod( ngm), qgm( ngm) )
  ALLOCATE (ylmk0( ngm, lmaxq * lmaxq) )

  aux (:,:) = (0.d0, 0.d0)
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  DO ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  ENDDO
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nij = nh(nt)*(nh(nt)+1)/2
        !
        ! count max number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
        ENDDO
        !
        ALLOCATE ( skk(ngm,nab), tbecsum(nij,nab,nspin_mag), aux2(ngm,nij) )
        !
        nb = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) THEN
              nb = nb + 1
              tbecsum(:,nb,:) = becsum(1:nij,na,1:nspin_mag)
!$omp parallel do default(shared) private(ig)
              DO ig = 1, ngm
                 skk(ig,nb) = eigts1 (mill (1,ig), na) * &
                              eigts2 (mill (2,ig), na) * &
                              eigts3 (mill (3,ig), na)
              ENDDO
!$omp end parallel do
           ENDIF
        ENDDO

        DO is = 1, nspin_mag
           ! sum over atoms
           CALL dgemm( 'N', 'T', 2*ngm, nij, nab, 1.0_dp, skk, 2*ngm,&
                tbecsum(1,1,is), nij, 0.0_dp, aux2, 2*ngm )
           ! sum over lm indices of Q_{lm}
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1
                 CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
!$omp parallel do default(shared) private(ig)
                 DO ig = 1, ngm
                    aux(ig,is) = aux(ig,is) + aux2(ig,ijh)*qgm(ig)
                 ENDDO
!$omp end parallel do
             ENDDO
           ENDDO
        ENDDO
        DEALLOCATE (aux2, tbecsum, skk )
     ENDIF
  ENDDO
  !
  DEALLOCATE (ylmk0)
  DEALLOCATE (qgm, qmod)
  !
  !     convert aux to real space and add to the charge density
  !
#ifdef DEBUG_ADDUSDENS
  CALL start_clock ('addus:fft')
#endif
  DO is = 1, nspin_mag
     psic(:) = (0.d0, 0.d0)
     psic( nl(:) ) = aux(:,is)
     IF (gamma_only) psic( nlm(:) ) = CONJG (aux(:,is))
     CALL invfft ('Dense', psic, dfftp)
     rho(:, is) = rho(:, is) +  DBLE (psic (:) )
  ENDDO
#ifdef DEBUG_ADDUSDENS
  CALL stop_clock ('addus:fft')
#endif
  DEALLOCATE (aux)

  CALL stop_clock ('addusdens')
  RETURN
END SUBROUTINE addusdens_g

