!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE vhpsi ( lda, n, m, psi, hpsi )
  !-----------------------------------------------------------------------
  !! This routine computes the Hubbard potential applied to the electronic
  !! structure of the current k-point. The result is added to hpsi.
  !! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  USE kinds,            ONLY : dp
  USE ldaU,             ONLY : lda_plus_u_kind, Hubbard_projectors
  USE noncollin_module, ONLY : npol, noncolin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(dp), INTENT(IN) :: psi(lda*npol,m)
  !! the wavefunction
  COMPLEX(dp), INTENT(INOUT) :: hpsi(lda*npol,m)
  !! Hamiltonian dot psi
  !
  IF ( Hubbard_projectors == "pseudo" ) RETURN
  !
  CALL start_clock('vhpsi')
  !
  IF ( lda_plus_u_kind == 0 .OR. lda_plus_u_kind == 1 ) THEN
     !
     CALL vhpsi_U ( lda, n, m, psi, hpsi )
     !
  ELSE IF ( noncolin ) THEN
     !$acc update host(psi, hpsi)
     CALL vhpsi_UV_nc( lda, n, m, psi, hpsi )
     !$acc update device(hpsi)
  ELSE
     !
     CALL vhpsi_UV( lda, n, m, psi, hpsi )
     !
  END IF
  !
  CALL stop_clock('vhpsi')
  !
  RETURN
  !
END SUBROUTINE vhpsi
!
!-----------------------------------------------------------------------
SUBROUTINE vhpsi_UV( lda, np, mps, psi, hpsi )
  !-----------------------------------------------------------------------
  !! This routine applies the Hubbard potential with U_I (=V_II) and V_IJ
  !! to the KS wave functions.  FIXME: not yet ported to GPU
  !! By taking a derivative of the Hubbard energy we obtain the potential
  !! (multiplied by the KS wave function psi_nk):
  !! - \sum_IJ (J\=I) \sum_{m1,m2} V_IJ/2 
  !!       * [ n^IJ_{m1,m2} * |phi^I_m1><phi^J_m2|psi_nk> +
  !!           n^JI_{m2,m1} * |phi^J_m2><phi^I_m1|psi_nk> ]
  !! Using the symmetry with respect to I/J and m1/m2 we can simplify 
  !! the expression for the Hubbard potential above and write it as:
  !! - \sum_IJ (J\=I) \sum_{m1,m2} V_IJ 
  !!       * n^IJ_{m1,m2} * |phi^I_m1><phi^J_m2|psi_nk>
  !! Since in practice the indices I and J are not really equivalent
  !! (due to the way how the generalized Hubbard potential is implemeted), 
  !! instead of the second expression here the first expression is implemented.
  !!
  !! Note: we assume that the phase factor phase_fac at a given k point
  !! has been already computed elsewhere.
  !
  USE kinds,         ONLY : dp
  USE becmod,        ONLY : bec_type, calbec, allocate_bec_type, &
                            deallocate_bec_type
  USE ldaU,          ONLY : Hubbard_l, is_hubbard,   &
                            nwfcU, wfcU, offsetU, lda_plus_u_kind, &
                            is_hubbard_back, Hubbard_l2, offsetU_back, &
                            backall, offsetU_back1
  USE ldaU,          ONLY : ldim_u, neighood, at_sc, phase_fac, Hubbard_V, v_nsg
  USE lsda_mod,      ONLY : current_spin
  USE scf,           ONLY : v
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE control_flags, ONLY : gamma_only
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psi
  COMPLEX(dp), INTENT(IN) :: psi(lda,mps)
  !! the wavefunction
  COMPLEX(dp), INTENT(INOUT) :: hpsi(lda,mps)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  REAL(dp),    ALLOCATABLE :: rtemp(:,:)
  COMPLEX(dp), ALLOCATABLE :: ctemp(:,:), vaux(:,:)
  TYPE(bec_type) :: proj
  !
  COMPLEX(dp) :: phase
  INTEGER :: ldim2, ldimx, ldim1,  m1, m2, equiv_na2, &
             off1, off2, ig, viz, na1, na2, nt1, nt2
  REAL(dp),    ALLOCATABLE :: projauxr(:,:), rvaux(:,:)
  COMPLEX(dp), ALLOCATABLE :: projauxc(:,:), wfcUaux(:,:)
  !
  ! Allocate the array proj = <wfcU|psi>
  CALL allocate_bec_type ( nwfcU,mps, proj )
  !
  !$acc update host(psi, hpsi)
  CALL calbec (np, wfcU, psi, proj)
  ! 
  ! Find the maximum number of magnetic quantum numbers [i.e. MAX(2l+1)]
  !
  ldimx = 0
  DO nt1 = 1, ntyp
     IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
        ldim1 = ldim_u(nt1)
        ldimx = MAX(ldimx,ldim1)
     ENDIF
  ENDDO
  !
  IF (gamma_only) THEN
     ALLOCATE (rtemp(ldimx,mps))
     ALLOCATE (projauxr(ldimx,mps))
     ALLOCATE (rvaux(ldimx,ldimx))
  ELSE
     ALLOCATE (ctemp(ldimx,mps))
     ALLOCATE (projauxc(ldimx,mps))
     ALLOCATE (vaux(ldimx,ldimx))
  ENDIF
  ! 
  ALLOCATE (wfcUaux(np,ldimx))
  !
  DO nt1 = 1, ntyp
     !
     ldim1 = ldim_u(nt1)
     !
     IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
        !
        DO na1 = 1, nat
           !
           IF (ityp(na1).EQ.nt1) THEN
              !
              DO viz = 1, neighood(na1)%num_neigh
                 !
                 na2 = neighood(na1)%neigh(viz)
                 equiv_na2 = at_sc(na2)%at
                 nt2 = ityp(equiv_na2)
                 phase = phase_fac(na2)
                 ldim2 = ldim_u(nt2)
                 !
                 ! Note: Below there is a condition on v_nsg
                 ! because it may be that the user specifies
                 ! Hubbard_alpha for some type for which 
                 ! Hubbard_V was not specified.
                 !
                 IF ( (is_hubbard(nt2).OR.is_hubbard_back(nt2)) .AND. &
                      (Hubbard_V(na1,na2,1).NE.0.d0 .OR. &
                       Hubbard_V(na1,na2,2).NE.0.d0 .OR. &
                       Hubbard_V(na1,na2,3).NE.0.d0 .OR. &
                       Hubbard_V(na1,na2,4).NE.0.d0 .OR. &
                    ANY(v_nsg(:,:,viz,na1,current_spin).NE.0.0d0)) ) THEN
                    !
                    ! Compute the first part of the Hubbard potential, namely:
                    ! - \sum_IJ (J\=I) \sum_{m1,m2} V_IJ/2 
                    !      * n^IJ_{m1,m2} * |phi^I_m1><phi^J_m2|psi_nk> 
                    ! where
                    ! - V_IJ * n^IJ_{m1,m2}   = CONJG(v_nsg)
                    !      <phi^J_m2|Psi_nk>  = proj%r (or proj%k) 
                    !         |phi^I_m1>      = wfcU
                    !
                    wfcUaux(:,:) = (0.0_dp, 0.0_dp)
                    !
                    off1 = offsetU(na1)
                    !
                    DO m1 = 1, ldim_u(nt1)
                       !
                       IF (m1.GT.2*Hubbard_l(nt1)+1) &
                          off1 = offsetU_back(na1) - 2*Hubbard_l(nt1) - 1
                       !
                       IF (backall(nt1) .AND. &
                           m1.GT.(2*Hubbard_l(nt1)+1+2*Hubbard_l2(nt1)+1)) &
                           off1 = offsetU_back1(na1) &
                                 - 2*Hubbard_l(nt1) - 2 - 2*Hubbard_l2(nt1)
                       !
                       DO ig = 1, np
                          wfcUaux(ig,m1) = wfcU(ig,off1+m1)
                       ENDDO
                       !
                    ENDDO 
                    !
                    off2 = offsetU(equiv_na2)
                    !
                    IF (gamma_only) THEN
                       !
                       rvaux(:,:) = 0.0_dp
                       !
                       projauxr(:,:) = 0.0_dp
                       !
                       DO m1 = 1, ldim1
                          !
                          DO m2 = 1, ldim2
                             !
                             rvaux(m2,m1) = DBLE( (v_nsg(m2, m1, viz, na1, current_spin))) * 0.5d0
                             !
                          ENDDO
                          !
                       ENDDO
                       !
                       DO m2 = 1, ldim2
                          !
                          IF (m2.GT.2*Hubbard_l(nt2)+1) &
                             off2 = offsetU_back(equiv_na2) - 2*Hubbard_l(nt2) - 1
                          !
                          IF (backall(nt2) .AND. &
                              m2.GT.(2*Hubbard_l(nt2)+1+2*Hubbard_l2(nt2)+1)) &
                              off2 = offsetU_back1(equiv_na2) &
                                     - 2*Hubbard_l(nt2) - 2 - 2*Hubbard_l2(nt2)
                          !
                          projauxr(m2,:) = DBLE(proj%r(off2+m2,:))
                          !
                       ENDDO
                       !
                       rtemp(:,:) = 0.0_dp
                       !
                       CALL DGEMM ('t','n', ldim1,mps,ldim2, 1.0_dp, &
                            rvaux,ldimx, projauxr,ldimx, 0.0_dp, rtemp, ldimx)
                       !
                       CALL DGEMM ('n','n', 2*np, mps, ldim1, 1.0_dp, &
                            wfcUaux, 2*np, rtemp, ldimx, &
                            1.0_dp, hpsi, 2*lda)
                       !
                    ELSE
                       !
                       vaux(:,:) = (0.0_dp, 0.0_dp)
                       !
                       projauxc(:,:) = (0.0_dp, 0.0_dp)
                       !
                       DO m1 = 1, ldim1
                          !
                          DO m2 = 1, ldim2
                             !
                             vaux(m2,m1) = CONJG( (v_nsg(m2, m1, viz, na1, current_spin))) * 0.5d0
                             !
                          ENDDO
                          !
                       ENDDO
                       !
                       DO m2 = 1, ldim2
                          !
                          IF (m2.GT.2*Hubbard_l(nt2)+1) &
                             off2 = offsetU_back(equiv_na2) - 2*Hubbard_l(nt2) - 1
                          !
                          IF (backall(nt2) .AND. &
                              m2.GT.(2*Hubbard_l(nt2)+1+2*Hubbard_l2(nt2)+1)) &
                              off2 = offsetU_back1(equiv_na2) &
                                     - 2*Hubbard_l(nt2) - 2 - 2*Hubbard_l2(nt2)
                          !
                          projauxc(m2,:) = proj%k(off2+m2,:)
                          !
                       ENDDO
                       !
                       ctemp(:,:) = (0.0_dp,0.0_dp)
                       !
                       CALL ZGEMM ('t','n', ldim1,mps,ldim2, (1.0_dp,0.0_dp), &
                            vaux,ldimx, projauxc,ldimx, (0.0_dp,0.0_dp), ctemp, ldimx)
                       !
                       CALL ZGEMM ('n','n', np, mps, ldim1, phase, &
                            wfcUaux, np, ctemp, ldimx, (1.0_dp,0.0_dp), hpsi, lda)
                       !
                    ENDIF
                    !
                    !
                    ! Compute the second part of the Hubbard potential, namely:
                    ! - \sum_IJ (J\=I) \sum_m1m2 V_IJ/2 * n^JI_m2m1 * |phi^J_m2><phi^I_m1|Psi_nk>
                    ! where
                    ! - V_IJ * n^JI_m2m1   = v_nsg
                    !   <phi^I_m1|Psi_nk>  = proj%r (or proj%k) 
                    !      |phi^J_m2>      = wfcU
                    !
                    wfcUaux(:,:) = (0.0_dp, 0.0_dp)
                    !
                    off2 = offsetU(equiv_na2)
                    ! 
                    DO m2 = 1, ldim_u(nt2)
                       !
                       IF (m2.GT.2*Hubbard_l(nt2)+1) &
                           off2 = offsetU_back(equiv_na2) - 2*Hubbard_l(nt2) - 1
                       !
                       IF (backall(nt2) .AND.  &
                           m2.GT.(2*Hubbard_l(nt2)+1+2*Hubbard_l2(nt2)+1)) &
                           off2 = offsetU_back1(equiv_na2) &
                                  - 2*Hubbard_l(nt2) - 2 - 2*Hubbard_l2(nt2)
                       !
                       DO ig = 1, np
                          wfcUaux(ig,m2) = wfcU(ig,off2+m2)
                       ENDDO
                       !
                    ENDDO 
                    ! 
                    off1 = offsetU(na1)
                    !
                    IF (gamma_only) THEN
                       !
                       projauxr(:,:) = 0.0_dp
                       !
                       DO m1 = 1, ldim1
                          !
                          IF (m1.GT.2*Hubbard_l(nt1)+1) &
                             off1 = offsetU_back(na1) - 2*Hubbard_l(nt1) - 1
                          !
                          IF (backall(nt1) .AND. &
                              m1.GT.(2*Hubbard_l(nt1)+1+2*Hubbard_l2(nt1)+1)) &
                              off1 = offsetU_back1(na1) &
                                     - 2*Hubbard_l(nt1) - 2 - 2*Hubbard_l2(nt1)
                          !
                          projauxr(m1,:) = DBLE(proj%r(off1+m1,:))
                          !
                       ENDDO
                       !
                       rvaux(:,:) = 0.0_dp
                       !
                       DO m1 = 1, ldim1
                          !
                          DO m2 = 1, ldim2
                             !
                             rvaux(m2,m1) = DBLE(v_nsg(m2, m1, viz, na1, current_spin)) * 0.5d0
                             !
                          ENDDO
                          !
                       ENDDO
                       !
                       rtemp(:,:) = 0.0_dp
                       !
                       CALL DGEMM ('n','n', ldim2,mps,ldim1, 1.0_dp, &
                            rvaux,ldimx, projauxr,ldimx, 0.0_dp, rtemp, ldimx)
                       !
                       CALL DGEMM ('n','n', 2*np, mps, ldim2, 1.0_dp, &
                            wfcUaux, 2*np, rtemp, ldimx, &
                            1.0_dp, hpsi, 2*lda)
                       !
                    ELSE
                       !
                       projauxc(:,:) = (0.0_dp,0.0_dp)
                       !
                       do m1 = 1,ldim1
                          !
                          IF (m1.GT.2*Hubbard_l(nt1)+1) &
                             off1 = offsetU_back(na1) - 2*Hubbard_l(nt1) - 1
                          !
                          IF (backall(nt1) .AND. &
                              m1.GT.(2*Hubbard_l(nt1)+1+2*Hubbard_l2(nt1)+1)) &
                              off1 = offsetU_back1(na1) &
                                     - 2*Hubbard_l(nt1) - 2 - 2*Hubbard_l2(nt1)
                          !
                          projauxc(m1,:) = proj%k(off1+m1,:)
                          !
                       end do
                       !
                       vaux(:,:) = (0.0_dp,0.0_dp)
                       !
                       DO m1 = 1,ldim1
                          !
                          DO m2 = 1,ldim2
                             !
                             vaux(m2,m1) = v_nsg(m2, m1, viz, na1, current_spin) * 0.5d0
                             !
                          ENDDO
                          !
                       ENDDO
                       !
                       ctemp(:,:) = (0.0_dp,0.0_dp)
                       !
                       CALL ZGEMM ('n','n', ldim2,mps,ldim1, (1.0_dp,0.0_dp), &
                            vaux,ldimx, projauxc,ldimx, (0.0_dp,0.0_dp), ctemp, ldimx)
                       !
                       CALL ZGEMM ('n','n', np, mps, ldim2, CONJG(phase), &
                            wfcUaux, np, ctemp, ldimx, (1.0_dp,0.0_dp), hpsi, lda)
                       !
                    ENDIF
                    !
                 ENDIF
                 !
              ENDDO ! viz
              !
           ENDIF !nt1 = ityp(na1)
           !
        ENDDO !na1
        !
     ENDIF
     !
  ENDDO
  !
  IF (gamma_only) THEN
     DEALLOCATE (rtemp)
     DEALLOCATE (projauxr)
     DEALLOCATE (rvaux)
  ELSE
     DEALLOCATE (ctemp)
     DEALLOCATE (projauxc)
     DEALLOCATE (vaux)
  ENDIF
  !
  DEALLOCATE (wfcUaux)
  !
  !$acc update device(hpsi)
  !
  CALL deallocate_bec_type (proj)
  !
  RETURN
  !
END SUBROUTINE vhpsi_UV
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
SUBROUTINE vhpsi_UV_nc( lda, np, mps, psi, hpsi )
  !-----------------------------------------------------------------------
  !! Noncollinear version of \(\texttt{vhpsi} routine (A. Smogunov).
  !! Extended to DFT+U+V and reorganised by L. Binci 
  !
  USE kinds,            ONLY: dp
  USE ldaU,             ONLY: Hubbard_lmax, Hubbard_l, is_Hubbard, nwfcU, &
                              wfcU, offsetU, lda_plus_u_kind, ldim_u, v_nsg, &
                              ldim_u, neighood, at_sc, phase_fac, Hubbard_V
  USE scf,              ONLY: v
  USE ions_base,        ONLY: nat, ntyp => nsp, ityp
  USE noncollin_module, ONLY: npol
  USE mp_bands,         ONLY: intra_bgrp_comm
  USE mp,               ONLY: mp_sum
  USE lsda_mod,         ONLY: nspin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psi
  COMPLEX(dp), INTENT(IN) :: psi(lda*npol,mps)
  !! the wavefunction
  COMPLEX(dp), INTENT(INOUT) :: hpsi(lda*npol,mps)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd, na, nwfc, is1, is2, nt, m1, m2
  INTEGER :: ldim2, ldimx, ldim1, equiv_na2, &
              off1, off2, ig, viz, na1, na2, nt1, nt2
  COMPLEX(dp) :: temp, phase
  COMPLEX(dp), ALLOCATABLE :: proj(:,:)
  COMPLEX(dp), ALLOCATABLE :: ctemp(:,:), vaux(:,:)
  COMPLEX(dp), ALLOCATABLE :: projauxc(:,:), wfcUaux(:,:)
  !
  !
  IF ( lda_plus_u_kind == 0 .OR. lda_plus_u_kind == 1 ) CALL errore('vhpsi','incorrectly called',2)
  !
  ALLOCATE( proj(nwfcU, mps) )
  proj(:,:) = (0.0_dp,0.0_dp)
  !
  ! calculate <psi_at | phi_k> 
  CALL ZGEMM ('C', 'N', nwfcU, mps, lda*npol, (1.0_dp, 0.0_dp), wfcU, &
                    lda*npol, psi, lda*npol, (0.0_dp, 0.0_dp),  proj, nwfcU)
  CALL mp_sum ( proj, intra_bgrp_comm )
  !
   !
   ! Find the maximum number of magnetic quantum numbers [i.e. MAX(2l+1)]
   !
   ldimx = 0
   DO nt1 = 1, ntyp
      IF ( is_hubbard(nt1) ) THEN
         ldim1 = ldim_u(nt1)
         ldimx = MAX(ldimx,ldim1)
      ENDIF
   ENDDO
   !
   ALLOCATE (ctemp(ldimx*npol,mps))
   ALLOCATE (projauxc(ldimx*npol,mps))
   ALLOCATE (vaux(ldimx*npol,ldimx*npol))
   !
   ALLOCATE (wfcUaux(lda*npol,ldimx*npol))
   !
   DO nt1 = 1, ntyp
      ldim1 = ldim_u(nt1)
      IF ( is_hubbard(nt1) ) THEN
         DO na1 = 1, nat
            IF (ityp(na1).EQ.nt1) THEN
                     DO viz = 1, neighood(na1)%num_neigh
                        !
                        na2 = neighood(na1)%neigh(viz)
                        equiv_na2 = at_sc(na2)%at
                        nt2 = ityp(equiv_na2)
                        phase = phase_fac(na2)
                        ldim2 = ldim_u(nt2)
                        !
                        IF ( (is_hubbard(nt2) ) .AND. &
                          (Hubbard_V(na1,na2,1).NE.0.d0 .OR. &
                           ANY(v_nsg(:,:,viz,na1,:).NE.0.0d0)) ) THEN
                           !
                           ! Compute the first part of the Hubbard potential, namely:
                           ! - \sum_IJ (J\=I) \sum_{m1,m2} V_IJ/2 
                           !      * n^IJ_{m1,m2} * |phi^I_m1><phi^J_m2|psi_nk> 
                           ! where
                           ! - V_IJ * n^IJ_{m1,m2}   = CONJG(v_nsg)
                           !      <phi^J_m2|Psi_nk>  = proj
                           !         |phi^I_m1>      = wfcU
                           !
                            !  DO is1 = 1, npol
                            !     DO is2 = 1, npol
                           wfcUaux(:,:) = (0.0_dp, 0.0_dp)
                           !
                           off1 = offsetU(na1)
                           DO m1 = 1, ldim1*npol
                              DO ig = 1, lda*npol
                                 wfcUaux(ig,m1) = wfcU(ig,off1+m1)
                              ENDDO
                           ENDDO
                           !write(stdout,*) "phase, wfcUaux(1,1:ldim1)", phase, wfcUaux(1,1:ldim1*npol)
                           !
                           off2 = offsetU(equiv_na2)
                           vaux(:,:) = (0.0_dp, 0.0_dp)
                           projauxc(:,:) = (0.0_dp, 0.0_dp)
                           do is1 =1 , npol
                              do is2 =1 , npol
                                 DO m1 = 1, ldim1
                                    DO m2 = 1, ldim2
                                       vaux(m2+ldim2*(is2-1),m1+ldim1*(is1-1)) = &
                                          CONJG( (v_nsg(m2, m1, viz, na1, npol*(is2-1)+is1))) * 0.5d0
                                    ENDDO
                                 ENDDO
                              enddo 
                           enddo
                           !
                           DO m2 = 1, ldim2*npol
                              projauxc(m2,:) = proj(off2+m2,:)
                           ENDDO
                           ctemp(:,:) = (0.0_dp,0.0_dp)
                           !
                           CALL ZGEMM ('t','n', ldim1*npol,mps,ldim2*npol, (1.0_dp,0.0_dp), &
                                vaux,ldimx*npol, projauxc,ldimx*npol, (0.0_dp,0.0_dp), ctemp, ldimx*npol)
                           !
                           CALL ZGEMM ('n','n', lda*npol, mps, ldim1*npol, phase, &
                                wfcUaux, lda*npol, ctemp, ldimx*npol, (1.0_dp,0.0_dp), hpsi, lda*npol)
                           !
                           ! Compute the second part of the Hubbard potential, namely:
                           ! - \sum_IJ (J\=I) \sum_m1m2 V_IJ/2 * n^JI_m2m1 * |phi^J_m2><phi^I_m1|Psi_nk>
                           ! where
                           ! - V_IJ * n^JI_m2m1   = v_nsg
                           !   <phi^I_m1|Psi_nk>  = proj
                           !      |phi^J_m2>      = wfcU
                           !
                           wfcUaux(:,:) = (0.0_dp, 0.0_dp)
                           off2 = offsetU(equiv_na2)   
                           DO m2 = 1, ldim2*npol
                              DO ig = 1, lda*npol
                                 wfcUaux(ig,m2) = wfcU(ig,off2+m2)
                              ENDDO
                           ENDDO 
                           ! 
                           off1 = offsetU(na1)
                           projauxc(:,:) = (0.0_dp,0.0_dp)
                           do m1 = 1,ldim1*npol
                              projauxc(m1,:) = proj(off1+m1,:)
                           enddo
                           vaux(:,:) = (0.0_dp,0.0_dp)
                           do is1=1,npol
                              do is2=1,npol
                                 DO m1 = 1,ldim1
                                    DO m2 = 1,ldim2
                                       vaux(m2+ldim2*(is2-1),m1+ldim1*(is1-1)) = &
                                            v_nsg(m2, m1, viz, na1, npol*(is2-1)+is1) * 0.5d0
                                    ENDDO
                                 ENDDO
                              enddo
                           enddo
                           !
                           ctemp(:,:) = (0.0_dp,0.0_dp)
                           !
                           CALL ZGEMM ('n','n', ldim2*npol,mps,ldim1*npol, (1.0_dp,0.0_dp), &
                                 vaux,ldimx*npol, projauxc,ldimx*npol, (0.0_dp,0.0_dp), ctemp, ldimx*npol)
                           !
                           CALL ZGEMM ('n','n', lda*npol, mps, ldim2*npol, CONJG(phase), &
                                 wfcUaux, lda*npol, ctemp, ldimx*npol, (1.0_dp,0.0_dp), hpsi, lda*npol)
                           !
                        ENDIF
                     ENDDO 
            ENDIF
         ENDDO
      ENDIF
   !   
   !   
   ENDDO
   !
   DEALLOCATE (ctemp)
   DEALLOCATE (projauxc)
   DEALLOCATE (vaux)
   !
   DEALLOCATE (wfcUaux)
   DEALLOCATE (proj)
   !
   RETURN
   !
  ! --------------------------------------------
END SUBROUTINE vhpsi_UV_nc

