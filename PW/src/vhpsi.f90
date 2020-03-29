!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE vhpsi( ldap, np, mps, psip, hpsi )
  !-----------------------------------------------------------------------
  !! This routine computes the Hubbard potential applied to the electronic
  !! structure of the current k-point. The result is added to hpsi.
  !
  USE kinds,         ONLY : DP
  USE becmod,        ONLY : bec_type, calbec, allocate_bec_type, &
                            deallocate_bec_type
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_Hubbard,   &
                            nwfcU, wfcU, offsetU, lda_plus_u_kind, &
                            is_hubbard_back, Hubbard_l_back, offsetU_back, &
                            backall, offsetU_back1
  USE lsda_mod,      ONLY : current_spin
  USE scf,           ONLY : v
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE control_flags, ONLY : gamma_only
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldap
  !! leading dimension of arrays psip, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psip, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psip
  COMPLEX(DP), INTENT(IN) :: psip(ldap,mps)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(ldap,mps)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  REAL(DP),    ALLOCATABLE :: rtemp(:,:)
  COMPLEX(DP), ALLOCATABLE :: ctemp(:,:), vaux(:,:)
  TYPE(bec_type) :: proj
  !
  CALL start_clock('vhpsi')
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  ! Allocate the array proj
  CALL allocate_bec_type ( nwfcU,mps, proj )
  !
  ! proj = <wfcU|psip>
  CALL calbec (np, wfcU, psip, proj)
  ! 
  IF ( lda_plus_u_kind.EQ.0 .OR. lda_plus_u_kind.EQ.1 ) THEN
     CALL vhpsi_U ()  ! DFT+U
  ELSEIF ( lda_plus_u_kind.EQ.2 ) THEN
     CALL vhpsi_UV () ! DFT+U+V
  ENDIF
  !
  CALL deallocate_bec_type (proj)
  !
  CALL stop_clock('vhpsi')
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE vhpsi_U ()
  !
  ! This routine applies the Hubbard potential with U_I
  ! to the KS wave functions. 
  !
  USE ldaU,      ONLY : ldim_back, ldmx_b, Hubbard_l1_back
  !
  IMPLICIT NONE
  INTEGER :: na, nt, ldim, ldim0
  !
  DO nt = 1, ntyp
     !
     ! Compute the action of the Hubbard potential on the KS wave functions:
     ! V_Hub |psip > = \sum v%ns |wfcU> <wfcU|psip>
     ! where v%ns = U ( delta/2 - rho%ns ) is computed in v_of_rho
     !
     IF ( is_hubbard(nt) ) THEN
        !  
        ldim = 2*Hubbard_l(nt) + 1
        !
        IF (gamma_only) THEN
           ALLOCATE ( rtemp(ldim,mps) )
        ELSE
           ALLOCATE ( ctemp(ldim,mps) )
        ENDIF
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              IF (gamma_only) THEN
                 CALL DGEMM ('n','n', ldim,mps,ldim, 1.0_dp, &
                      v%ns(1,1,current_spin,na),2*Hubbard_lmax+1, &
                      proj%r(offsetU(na)+1,1),nwfcU, 0.0_dp, rtemp, ldim)
                 CALL DGEMM ('n','n', 2*np, mps, ldim, 1.0_dp, &
                      wfcU(1,offsetU(na)+1), 2*ldap, rtemp, ldim, &
                      1.0_dp, hpsi, 2*ldap)
              ELSE
                 !
                 ALLOCATE(vaux(ldim,ldim))
                 !
                 vaux = (0.0_dp,0.0_dp)
                 vaux(:,:) = v%ns(:,:,current_spin,na)
                 !
                 CALL ZGEMM ('n','n', ldim, mps, ldim, (1.0_dp,0.0_dp), &
                      vaux, ldim, proj%k(offsetU(na)+1,1), nwfcU, &
                      (0.0_dp,0.0_dp), ctemp, ldim)
                 !
                 DEALLOCATE(vaux)
                 !
                 CALL ZGEMM ('n','n', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU(1,offsetU(na)+1), ldap, ctemp, ldim, &
                      (1.0_dp,0.0_dp), hpsi, ldap)
                 !
              ENDIF
           ENDIF
        ENDDO
        !
        IF (gamma_only) THEN
           DEALLOCATE ( rtemp )
        ELSE
           DEALLOCATE ( ctemp )
        ENDIF
        !
     ENDIF
     !
     ! If the background is used then compute extra 
     ! contribution to the Hubbard potential
     !
     IF ( is_hubbard_back(nt) ) THEN
        !
        ldim = ldim_back(nt)
        !
        IF (gamma_only) THEN
           ALLOCATE ( rtemp(ldim,mps) )
        ELSE
           ALLOCATE ( ctemp(ldim,mps) )
        ENDIF
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              !
              IF (gamma_only) THEN
                 !
                 ldim = 2*Hubbard_l_back(nt)+1
                 !
                 CALL DGEMM ('n','n', ldim,mps,ldim, 1.0_dp, &
                      v%nsb(1,1,current_spin,na),ldmx_b, &
                      proj%r(offsetU_back(na)+1,1), &
                      nwfcU, 0.0_dp, rtemp, ldim_back(nt))
                 !
                 CALL DGEMM ('n','n', 2*np, mps, ldim, 1.0_dp,   &
                      wfcU(1,offsetU_back(na)+1), 2*ldap, rtemp, &
                      ldim_back(nt), 1.0_dp, hpsi, 2*ldap)
                 !
                 IF (backall(nt)) THEN
                    !
                    ldim  = 2*Hubbard_l1_back(nt)+1
                    ldim0 = 2*Hubbard_l_back(nt)+1
                    !
                    CALL DGEMM ('n','n', ldim,mps,ldim, 1.0_dp,         &
                         v%nsb(ldim0+1,ldim0+1,current_spin,na),        &
                         ldim_back(nt), proj%r(offsetU_back1(na)+1,1),  &
                         nwfcU, 0.0_dp, rtemp, ldim_back(nt))
                    !
                    CALL DGEMM ('n','n', 2*np, mps, ldim, 1.0_dp,       &
                         wfcU(1,offsetU_back1(na)+1), 2*ldap, rtemp,    &
                         ldim_back(nt), 1.0_dp, hpsi, 2*ldap)
                    !
                 ENDIF
                 !
              ELSE
                 !
                 ldim = ldim_back(nt)
                 !
                 ALLOCATE(vaux(ldim,ldim))
                 !
                 vaux = (0.0_dp, 0.0_dp)
                 vaux(:,:) = v%nsb(:,:,current_spin,na)
                 !
                 ldim = 2*Hubbard_l_back(nt)+1
                 !
                 CALL ZGEMM ('n','n', ldim,mps,ldim, (1.0_dp,0.0_dp),   &
                      vaux,ldim_back(nt), proj%k(offsetU_back(na)+1,1), &
                      nwfcU, (0.0_dp,0.0_dp), ctemp, ldim_back(nt))
                 !
                 CALL ZGEMM ('n','n', np, mps, ldim, (1.0_dp,0.0_dp),   &
                      wfcU(1,offsetU_back(na)+1), ldap, ctemp,          &
                      ldim_back(nt), (1.0_dp,0.0_dp), hpsi, ldap)
                 !
                 IF (backall(nt)) THEN
                    !
                    ldim  = 2*Hubbard_l1_back(nt)+1
                    ldim0 = 2*Hubbard_l_back(nt)+1
                    !
                    CALL ZGEMM ('n','n', ldim,mps,ldim,(1.0_dp,0.0_dp), &
                         vaux(ldim0+1,ldim0+1),ldim_back(nt),           &
                         proj%k(offsetU_back1(na)+1,1), nwfcU,          &
                         (0.0_dp,0.0_dp), ctemp, ldim_back(nt))
                    ! 
                    CALL ZGEMM ('n','n', np, mps, ldim, (1.0_dp,0.0_dp),&
                         wfcU(1,offsetU_back1(na)+1), ldap, ctemp,      &
                         ldim_back(nt), (1.0_dp,0.0_dp), hpsi, ldap)
                    !
                 ENDIF
                 !
                 DEALLOCATE(vaux)
                 !
              ENDIF
           ENDIF
        ENDDO
        !
        IF (gamma_only) THEN
           DEALLOCATE ( rtemp )
        ELSE
           DEALLOCATE ( ctemp )
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE vhpsi_U
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
SUBROUTINE vhpsi_UV ()
  !
  ! This routine applies the Hubbard potential with U_I (=V_II) and V_IJ
  ! to the KS wave functions.
  ! By taking a derivative of the Hubbard energy we obtain the potential
  ! (multiplied by the KS wave function psi_nk):
  ! - \sum_IJ (J\=I) \sum_{m1,m2} V_IJ/2 
  !       * [ n^IJ_{m1,m2} * |phi^I_m1><phi^J_m2|psi_nk> +
  !           n^JI_{m2,m1} * |phi^J_m2><phi^I_m1|psi_nk> ]
  ! Using the symmetry with respect to I/J and m1/m2 we can simplify 
  ! the expression for the Hubbard potential above and write it as:
  ! - \sum_IJ (J\=I) \sum_{m1,m2} V_IJ 
  !       * n^IJ_{m1,m2} * |phi^I_m1><phi^J_m2|psi_nk>
  ! Since in practice the indices I and J are not really equivalent
  ! (due to the way how the generalized Hubbard potential is implemeted), 
  ! instead of the second expression here the first expression is implemented.
  !
  ! Note: This routine assumes that the phase factor phase_fac at a given k point
  ! has been already computed elsewhere.
  !
  USE ldaU,      ONLY : ldim_u, neighood, at_sc, phase_fac, Hubbard_V, v_nsg
  !
  IMPLICIT NONE
  COMPLEX(DP) :: phase
  INTEGER :: ldim2, ldimx, ldim1,  m1, m2, equiv_na2, &
             off1, off2, ig, viz, na1, na2, nt1, nt2
  REAL(DP),    ALLOCATABLE :: projauxr(:,:), rvaux(:,:)
  COMPLEX(DP), ALLOCATABLE :: projauxc(:,:), wfcUaux(:,:)
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
                    ! - V_IJ/2 * n^IJ_{m1,m2} = CONJG(v_nsg)
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
                           m1.GT.(2*Hubbard_l(nt1)+1+2*Hubbard_l_back(nt1)+1)) &
                           off1 = offsetU_back1(na1) &
                                 - 2*Hubbard_l(nt1) - 2 - 2*Hubbard_l_back(nt1)
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
                              m2.GT.(2*Hubbard_l(nt2)+1+2*Hubbard_l_back(nt2)+1)) &
                              off2 = offsetU_back1(equiv_na2) &
                                     - 2*Hubbard_l(nt2) - 2 - 2*Hubbard_l_back(nt2)
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
                            1.0_dp, hpsi, 2*ldap)
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
                              m2.GT.(2*Hubbard_l(nt2)+1+2*Hubbard_l_back(nt2)+1)) &
                              off2 = offsetU_back1(equiv_na2) &
                                     - 2*Hubbard_l(nt2) - 2 - 2*Hubbard_l_back(nt2)
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
                            wfcUaux, np, ctemp, ldimx, (1.0_dp,0.0_dp), hpsi, ldap)
                       !
                    ENDIF
                    !
                    !
                    ! Compute the second part of the Hubbard potential, namely:
                    ! - \sum_IJ (J\=I) \sum_m1m2 V_IJ/2 * n^JI_m2m1 * |phi^J_m2><phi^I_m1|Psi_nk>
                    ! where
                    ! - V_IJ/2 * n^JI_m2m1 = v_nsg
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
                           m2.GT.(2*Hubbard_l(nt2)+1+2*Hubbard_l_back(nt2)+1)) &
                           off2 = offsetU_back1(equiv_na2) &
                                  - 2*Hubbard_l(nt2) - 2 - 2*Hubbard_l_back(nt2)
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
                              m1.GT.(2*Hubbard_l(nt1)+1+2*Hubbard_l_back(nt1)+1)) &
                              off1 = offsetU_back1(na1) &
                                     - 2*Hubbard_l(nt1) - 2 - 2*Hubbard_l_back(nt1)
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
                            1.0_dp, hpsi, 2*ldap)
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
                              m1.GT.(2*Hubbard_l(nt1)+1+2*Hubbard_l_back(nt1)+1)) &
                              off1 = offsetU_back1(na1) &
                                     - 2*Hubbard_l(nt1) - 2 - 2*Hubbard_l_back(nt1)
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
                            wfcUaux, np, ctemp, ldimx, (1.0_dp,0.0_dp), hpsi, ldap)
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
  RETURN
  !
END SUBROUTINE vhpsi_UV
!-------------------------------------------------------------------------  
END SUBROUTINE vhpsi
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
SUBROUTINE vhpsi_nc( ldap, np, mps, psip, hpsi )
  !-----------------------------------------------------------------------
  !! Noncollinear version of \(\texttt{vhpsi} routine (A. Smogunov).
  !
  USE kinds,            ONLY: DP
  USE ldaU,             ONLY: Hubbard_lmax, Hubbard_l, is_Hubbard, nwfcU, &
                              wfcU, offsetU
  USE scf,              ONLY: v
  USE ions_base,        ONLY: nat, ntyp => nsp, ityp
  USE noncollin_module, ONLY: npol
  USE mp_bands,         ONLY: intra_bgrp_comm
  USE mp,               ONLY: mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldap
  !! leading dimension of arrays psip, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psip, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psip
  COMPLEX(DP), INTENT(IN) :: psip(ldap*npol,mps)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(ldap*npol,mps)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd, na, nwfc, is1, is2, nt, m1, m2
  COMPLEX(DP) :: temp
  COMPLEX(DP), ALLOCATABLE :: proj(:,:)
  !
  CALL start_clock('vhpsi')
  !
  ALLOCATE( proj(nwfcU, mps) )
  !
!-- FIXME: to be replaced with ZGEMM
! calculate <psi_at | phi_k> 
  DO ibnd = 1, mps
    DO na = 1, nwfcU
       proj(na, ibnd) = dot_product( wfcU(1:ldap*npol, na), psip(1:ldap*npol, ibnd))
    ENDDO
  ENDDO
#if defined(__MPI)
  CALL mp_sum ( proj, intra_bgrp_comm )
#endif
!--

  do ibnd = 1, mps  
    do na = 1, nat  
       nt = ityp (na)  
       if ( is_hubbard(nt) ) then  
          nwfc = 2 * Hubbard_l(nt) + 1

          do is1 = 1, npol
           do m1 = 1, nwfc 
             temp = 0.d0
             do is2 = 1, npol
              do m2 = 1, nwfc  
                temp = temp + v%ns_nc( m1, m2, npol*(is1-1)+is2, na) * &
                              proj(offsetU(na)+(is2-1)*nwfc+m2, ibnd)
              enddo
             enddo
             call zaxpy (ldap*npol, temp, wfcU(1,offsetU(na)+(is1-1)*nwfc+m1),&
                         1, hpsi(1,ibnd),1)
           enddo
          enddo

       endif
    enddo
  enddo

  deallocate (proj)
  CALL stop_clock('vhpsi')

  return
  !
END SUBROUTINE vhpsi_nc

