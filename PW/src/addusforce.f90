!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusforce (forcenl)
  !----------------------------------------------------------------------
  !
  USE kinds,        ONLY : dp
  USE ions_base,    ONLY : nat
  USE control_flags,ONLY : tqr
  USE realus,       ONLY : addusforce_r
  !
  IMPLICIT NONE
  REAL(dp), INTENT(INOUT) :: forcenl (3, nat)
  !
  IF ( tqr ) THEN
     CALL addusforce_r (forcenl)
  ELSE
     CALL addusforce_g (forcenl)
  END IF
  !
END SUBROUTINE addusforce
!
!----------------------------------------------------------------------
SUBROUTINE addusforce_g (forcenl)
  !----------------------------------------------------------------------
  !
  !   This routine computes the contribution to atomic forces due
  !   to the dependence of the Q function on the atomic position.
  !   F_j,at= sum_G sum_lm iG_j exp(-iG*R_at) V^*(G) Q_lm(G) becsum(lm,at)
  !   where becsum(lm,at) = sum_i <psi_i|beta_l>w_i<beta_m|psi_i>
  !   On output: the contribution is added to forcenl
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega, tpiba
  USE fft_base,   ONLY : dfftp
  USE gvect,      ONLY : ngm, nl, nlm, gg, g, eigts1, eigts2, eigts3, mill
  USE noncollin_module,   ONLY : nspin_mag
  USE scf,        ONLY : v, vltot
  USE uspp,       ONLY : becsum, okvan
  USE uspp_param, ONLY : upf, lmaxq, nh, nhm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE control_flags, ONLY : gamma_only
  USE fft_interfaces,ONLY : fwfft
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl (3, nat)
  !
  INTEGER :: ig, nt, ih, jh, ijh, nij, ipol, is, na, nb, nab
  REAL(DP) :: fact
  COMPLEX(DP) :: cfac
  ! work space
  COMPLEX(DP), ALLOCATABLE :: aux(:), aux1(:,:,:), vg(:,:), qgm(:,:)
  REAL(DP) , ALLOCATABLE :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:), forceq(:,:)
  !
  IF (.not.okvan) RETURN
  !
  ALLOCATE ( forceq(3,nat) )
  forceq(:,:) = 0.0_dp
  IF (gamma_only) THEN
     fact = 2.d0*omega
  ELSE
     fact = 1.d0*omega
  ENDIF
  !
  ! fourier transform of the total effective potential
  !
  ALLOCATE ( vg(ngm,nspin_mag) )
  ALLOCATE ( aux(dfftp%nnr) )
  DO is = 1, nspin_mag
     IF (nspin_mag==4.and.is/=1) THEN
        aux(:) = v%of_r(:,is)
     ELSE
        aux(:) = vltot (:) + v%of_r (:, is)
     ENDIF
     CALL fwfft ('Dense', aux, dfftp)
     ! Note the factors -i and 2pi/a *units of G) here in V(G) !
     vg (:, is) = aux(nl (:) ) * tpiba * (0.d0, -1.d0)
  ENDDO
  DEALLOCATE (aux)
  !
  ALLOCATE (ylmk0(ngm,lmaxq*lmaxq))
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  !
  ALLOCATE (qmod( ngm))
!$omp parallel do default(shared) private(ig)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  ENDDO
!$omp end parallel do
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        ! qgm contains the Q functions in G space
        !
        nij = nh(nt)*(nh(nt)+1)/2
        ALLOCATE (qgm(ngm,nij))
        ijh = 0
        DO ih = 1, nh (nt)
           DO jh = ih, nh (nt)
              ijh = ijh + 1
              CALL qvan2 (ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0)
           ENDDO
        ENDDO
        !
        ! nab = number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
        ENDDO
        ALLOCATE ( aux1( ngm, na, 3) )
        ALLOCATE ( ddeeq(nij, nab, 3, nspin_mag) )
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                 nb = nb + 1
                 !
                 ! aux1 = product of potential, structure factor and iG
                 !
!$omp parallel do default(shared) private(ig, cfac)
                 do ig = 1, ngm
                    cfac = vg (ig, is) * CONJG(eigts1 (mill(1,ig),na) * &
                                               eigts2 (mill(2,ig),na) * &
                                               eigts3 (mill(3,ig),na) )
                    aux1 (ig, nb, 1) = g (1, ig) * cfac
                    aux1 (ig, nb, 2) = g (2, ig) * cfac
                    aux1 (ig, nb, 3) = g (3, ig) * cfac
                 ENDDO
!$omp end parallel do
                 !
              ENDIF
           ENDDO
           !
           !    ddeeq = dot product of aux1 with the Q functions
           !    No need for special treatment of the G=0 term (is zero)
           !
           DO ipol = 1, 3
              CALL DGEMM( 'C', 'N', nij, nab, 2*ngm, fact, qgm, 2*ngm, &
                   aux1(1,1,ipol), 2*ngm, 0.0_dp, ddeeq(1,1,ipol,is), nij )
           ENDDO
           !
        ENDDO
        !
        DEALLOCATE (aux1)
        DEALLOCATE (qgm)
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                 nb = nb + 1
                 DO ipol = 1, 3
                    DO ijh = 1, nij
                       forceq (ipol, na) = forceq (ipol, na) + &
                            ddeeq (ijh, nb, ipol, is) * becsum (ijh, na, is)
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
        DEALLOCATE ( ddeeq )
     ENDIF
  ENDDO
  !
  CALL mp_sum ( forceq, intra_bgrp_comm )
  !
  forcenl(:,:) = forcenl(:,:) + forceq(:,:)
  !
  DEALLOCATE (qmod)
  DEALLOCATE (ylmk0)
  DEALLOCATE (vg)
  DEALLOCATE (forceq)

  RETURN
END SUBROUTINE addusforce_g
