!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE new_nsg()
  !-----------------------------------------------------------------------
  !! This routine computes the new value for nsgnew (the occupation matrices)
  !! of the DFT+U+V approach.
  !! These quantities are defined as follows: 
  !!
  !! nsg_{I,m1,J,m2,s} = \sum_{k,v}
  !! f_{kv} <\fi^{at}_{I,m1}|\psi_{k,v,s}><\psi_{k,v,s}|\fi^{at}_{J,m2}>
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : Hubbard_l, q_ae, Hubbard_projectors, wfcU, nwfcU,     &
                                   ldim_u, ll, neighood, at_sc, nsgnew, phase_fac, &
                                   max_num_neighbors, Hubbard_l2, backall,     &
                                   offsetU, offsetU_back, offsetU_back1, is_hubbard_back
  USE symm_base,            ONLY : d1, d2, d3, t_rev
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE symm_base,            ONLY : nsym, irt
  USE wvfct,                ONLY : nbnd, npw, npwx, wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE becmod,               ONLY : bec_type, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  USE noncollin_module,     ONLY : npol, colin_mag
#if defined (__OSCDFT)
  USE plugin_flags,         ONLY : use_oscdft
  USE oscdft_base,          ONLY : oscdft_ctx
#endif 
  !
  IMPLICIT NONE
  !
  TYPE (bec_type) :: proj     
  ! proj(nwfcU,nbnd)
  INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, m1, m2, m11, m22, m0, m00, ldim,i_type, is2
  INTEGER :: na1, na2, viz, nt1, nt2, ldim1, ldim2, ldimb, nb1, nb2, viz_b
  INTEGER :: off, off1, off2, off3, eq_na2
  INTEGER :: ldim_std1, ldim_std2
  ! in the nwfcU ordering
  COMPLEX(DP), ALLOCATABLE :: nrg(:,:,:,:,:), nrg_nc(:,:,:,:,:,:)
  COMPLEX(DP) :: phase
  INTEGER, EXTERNAL :: find_viz
  !
  CALL start_clock('new_nsg')
  !
  ldim = 0
  DO nt = 1, ntyp
     ldim = MAX(ldim,ldim_u(nt))
  ENDDO
  !
  ALLOCATE ( nrg (ldim, ldim, max_num_neighbors, nat, nspin) )
  nrg(:,:,:,:,:) = (0.d0, 0.d0)
  nsgnew (:,:,:,:,:) = (0.d0, 0.d0)
  !
  CALL allocate_bec_type ( nwfcU, nbnd, proj ) 
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  ! we start a loop on k points
  !
  DO ik = 1, nks
     !
     IF (lsda) current_spin = isk(ik)
     !
     npw = ngk(ik)
     !
     IF (nks > 1) CALL get_buffer (evc, nwordwfc, iunwfc, ik)
     !
     ! make the projection
     !
     IF ( Hubbard_projectors == 'pseudo' ) THEN
        CALL errore('new_nsg', 'Hubbard_projectors = pseudo is not supported',1)
        !CALL compute_pproj( ik, q_ae, proj )
     ELSE
        IF (nks > 1) CALL get_buffer (wfcU, nwordwfcU, iunhub, ik)
        CALL calbec ( npw, wfcU, evc, proj )
     ENDIF
     !
     ! compute the phase factors for this k and put the result 
     ! in the phase_fac array
     !
     CALL phase_factor(ik)
     !
     ! compute the occupation matrix (nsg_{I,m1,J,m2,s}) of the
     ! atomic orbitals
     !
     DO ibnd = 1, nbnd
        !
        DO na1 = 1, nat  
           ! 
           nt1 = ityp(na1)
           !
           IF (ldim_u(nt1).GT.0) THEN
              !
              ldim1 = ldim_u(nt1)
              !
              DO viz =1, neighood(na1)%num_neigh
                 !
                 na2 = neighood(na1)%neigh(viz)
                 eq_na2 = at_sc(na2)%at
                 nt2 = ityp(eq_na2)
                 ldim2 = ldim_u(nt2)
                 ! 
                 IF (na1.GT.na2) THEN
                    !
                    DO m1 = 1, ldim1
                       DO m2 = 1, ldim2
                          nrg(m2,m1,viz,na1,current_spin) = &
                                CONJG(nrg(m1,m2,find_viz(na2,na1),na2,current_spin))
                       ENDDO
                    ENDDO
                    !
                 ELSE
                    !
                    phase = phase_fac(na2)
                    !
                    DO m1 = 1, ldim1
                       !
                       off1 = offsetU(na1) + m1
                       !
                       IF ( is_hubbard_back(nt1) .AND. m1.GT.2*Hubbard_l(nt1)+1 ) THEN
                          !
                          off1 = offsetU_back(na1) + m1 - 2*Hubbard_l(nt1) - 1
                          !
                          IF ( backall(nt1) .AND. &
                               m1.GT.2*(Hubbard_l(nt1)+Hubbard_l2(nt1)+1)) &
                             off1 = offsetU_back1(na1) + m1 - &
                                    2*(Hubbard_l(nt1)+Hubbard_l2(nt1)+1)
                          !
                       ENDIF
                       !
                       DO m2 = 1, ldim2
                          !
                          off2 = offsetU(eq_na2) + m2
                          !
                          IF ( is_hubbard_back(nt1) .AND. m2.GT.2*Hubbard_l(nt2)+1 ) THEN
                             !
                             off2 = offsetU_back(eq_na2) + m2 - 2*Hubbard_l(nt2) - 1
                             !
                             IF ( backall(nt2) .AND. & 
                                  m2.GT.2*(Hubbard_l(nt2)+Hubbard_l2(nt2)+1)) &
                                off2 = offsetU_back1(eq_na2) + m2 - &
                                       2*(Hubbard_l(nt2)+Hubbard_l2(nt2)+1)
                             !
                          ENDIF
                          !
                          IF (gamma_only) THEN
                             nrg(m2,m1,viz,na1,current_spin) = &
                             nrg(m2,m1,viz,na1,current_spin) + &
                             DBLE( wg(ibnd,ik) * ( proj%r(off1,ibnd) * &
                             CONJG(proj%r(off2,ibnd) * phase) ) )  
                          ELSE
                             nrg(m2,m1,viz,na1,current_spin) = &
                             nrg(m2,m1,viz,na1,current_spin) + &
                             DBLE( wg(ibnd,ik) * ( proj%k(off1,ibnd) * &
                             CONJG(proj%k(off2,ibnd) * phase) ) )
                          ENDIF
                          !
                       ENDDO ! m2
                    ENDDO ! m1
                    !
                 ENDIF
                 !
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     !
  ENDDO ! ik
  !
  CALL deallocate_bec_type (proj) 
  !
  CALL mp_sum( nrg, inter_pool_comm )
  !
  IF (nspin.EQ.1) nrg = 0.5d0 * nrg
  !
  ! impose hermiticity of n_{m1,m2}
  !
  DO ibnd = 1,nbnd
     DO na1 = 1, nat 
        nt1 = ityp(na1)
        IF (ldim_u(nt1).GT.0) THEN
           ldim1 = ldim_u(nt1) 
           DO viz = 1, neighood(na1)%num_neigh
              na2 = neighood(na1)%neigh(viz)
              IF (na1.GT.na2) THEN
                 eq_na2 = at_sc(na2)%at
                 nt2 = ityp (eq_na2)
                 ldim2 = ldim_u(nt2)
                 DO m1 = 1, ldim1
                    DO m2 = 1, ldim2
                       nrg(m2,m1,viz,na1,current_spin) = &
                       ( nrg(m2,m1,viz,na1,current_spin) + &
                       CONJG(nrg(m1,m2,find_viz(na2,na1),na2,current_spin)) )*0.5d0
                       nrg(m1,m2,find_viz(na2,na1),na2,current_spin) =  &
                       CONJG(nrg(m2,m1,viz,na1,current_spin))
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2) .AND. &
                  .NOT.oscdft_ctx%inp%constraint_diag) THEN
     IF (.NOT.oscdft_ctx%conv) THEN
        nsgnew(:,:,:,:,:) = nrg(:,:,:,:,:)
        DEALLOCATE(nrg)
        RETURN
     ENDIF
  ENDIF
#endif
  ! 
  ! symmetrize the quantities nr -> ns
  !
  DO is = 1, nspin
     !
     DO na1 = 1, nat
        !
        nt1 = ityp(na1)
        !
        IF (ldim_u(nt1).GT.0) THEN 
           !
           ldim_std1 = 2*Hubbard_l(nt1)+1
           !
           DO viz = 1, neighood(na1)%num_neigh
              !
              na2 = neighood(na1)%neigh(viz)
              eq_na2 = at_sc(na2)%at
              nt2 = ityp(eq_na2)
              ldim_std2 = 2*Hubbard_l(nt2)+1
              !
              IF (na1.GT.na2) THEN
                 !
                 ! we don't need to compute again
                 DO m1 = 1, ldim_u(nt1)
                    DO m2 = 1, ldim_u(nt2)
                       nsgnew(m2,m1,viz,na1,is) = &
                            CONJG(nsgnew(m1,m2,find_viz(na2,na1),na2,is))
                    ENDDO
                 ENDDO
                 !
              ELSE
                 !
                 DO m1 = 1, ldim_u(nt1)
                    !
                    off  = 1  
                    off2 = 2*Hubbard_l(nt1)+1
                    !
                    IF ( is_hubbard_back(nt1) .AND. m1.GT.ldim_std1 ) THEN
                       !
                       off  = ldim_std1 + 1
                       off2 = ldim_std1 + 2*Hubbard_l2(nt1) + 1 
                       !
                       IF ( backall(nt1) .AND. m1.GT.(ldim_std1+2*Hubbard_l2(nt1)+1) ) THEN
                          off  = ldim_std1 + 2*Hubbard_l2(nt1) + 2
                          off2 = ldim_u(nt1)
                       ENDIF
                       !
                    ENDIF
                    !
                    DO m2 = 1, ldim_u(nt2)
                       !
                       off1 = 1 
                       off3 = 2*Hubbard_l(nt2)+1
                       !
                       IF ( is_hubbard_back(nt2) .AND. m2.GT.ldim_std2 ) THEN
                          !
                          off1 = ldim_std2 + 1
                          off3 = ldim_std2 + 2*Hubbard_l2(nt2) + 1
                          !
                          IF ( backall(nt2) .AND. m2.GT.(ldim_std2+2*Hubbard_l2(nt2)+1) ) THEN
                             off1 = ldim_std2 + 2*Hubbard_l2(nt2) + 2
                             off3 = ldim_u(nt2)
                          ENDIF
                          !
                       ENDIF
                       !
                       ! Perform symmetrization using all available symmetries
                       !
                       DO isym = 1, nsym
                          !
                          CALL symonpair(na1,na2,isym,nb1,nb2)
                          !
                          viz_b = find_viz(nb1,nb2)
                          ! flip spin for time-reversal in collinear case
                          IF ( (colin_mag == 2) .AND. (t_rev(isym) == 1) ) THEN
                             is2 = 3 - is
                          ELSE
                             is2 = is
                          ENDIF
                          ! 
                          DO m0 = off, off2
                             DO m00 = off1, off3
                                IF (ll(m1,nt1).EQ.0 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.0) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        nrg(m00,m0,viz_b,nb1,is2) / nsym
                                ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.0) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d1(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        1.d0 / nsym
                                ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.0) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d2(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        1.d0 / nsym
                                ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.0) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d3(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        1.d0 / nsym

                                ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.1) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d1(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.1) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d1(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d1(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.1) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d2(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d1(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.1) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d3(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d1(m00-off1+1,m2-off1+1,isym) / nsym

                                ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.2) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d2(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.2) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d1(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d2(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.2) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d2(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d2(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.2) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d3(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d2(m00-off1+1,m2-off1+1,isym) / nsym

                                ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.3) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d3(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.3) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d1(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d3(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                  ll(m2,nt2).EQ.3) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d2(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d3(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.3) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d3(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is2) * &
                                        d3(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE
                                   CALL errore ('new_nsg', &
                                        'angular momentum not implemented for at least one type', &
                                        ABS(Hubbard_l(nt1)) )
                                END IF 
                             ENDDO !m00
                          ENDDO !m0
                          !
                       ENDDO !isym
                       !
                    ENDDO !m2
                    !
                 ENDDO  !m1
                 !
              ENDIF !na1 > na2
              ! 
           ENDDO !viz
           !
        ENDIF  !ldim_u > 0
        !
     ENDDO !na1
     !
  ENDDO !is
  !
  DEALLOCATE(nrg)
  !
  CALL stop_clock('new_nsg')
  !
  RETURN
  !
END SUBROUTINE new_nsg 
!
SUBROUTINE new_nsg_nc()
   !-----------------------------------------------------------------------
   !! This routine computes the new value for nsgnew (the occupation matrices)
   !! of the noncollinear DFT+U+V approach.
   !! These quantities are defined as follows: 
   !!
   !! nsg_{I,m1,J,m2,s} = \sum_{k,v}
   !! f_{kv} <\fi^{at}_{I,m1}|\psi_{k,v,s}><\psi_{k,v,s}|\fi^{at}_{J,m2}>
   !
   USE io_global,            ONLY : stdout
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp, ntyp => nsp
   USE klist,                ONLY : nks, ngk
   USE ldaU,                 ONLY : Hubbard_l, q_ae, Hubbard_projectors, wfcU, nwfcU,     &
                                    ldim_u, ll, neighood, at_sc, nsgnew, phase_fac, &
                                    max_num_neighbors, Hubbard_l2, backall,     &
                                    offsetU, offsetU_back, offsetU_back1, is_hubbard_back, &
                                    d_spin_ldau
   USE symm_base,            ONLY : d1, d2, d3
   USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
   USE symm_base,            ONLY : nsym, irt, t_rev
   USE wvfct,                ONLY : nbnd, npw, npwx, wg
   USE control_flags,        ONLY : gamma_only
   USE wavefunctions,        ONLY : evc
   USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
   USE buffers,              ONLY : get_buffer
   USE mp_global,            ONLY : inter_pool_comm
   USE mp,                   ONLY : mp_sum
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE becmod,               ONLY : bec_type, calbec, &
                                    allocate_bec_type, deallocate_bec_type
   USE noncollin_module,     ONLY : npol, noncolin, domag, lspinorb
   !
   IMPLICIT NONE
   !
   TYPE (bec_type) :: proj     
   ! proj(nwfcU,nbnd)
   INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, m1, m2, m11, m22, m3, m4, ldim,i_type, &
              is1, is2, is3, is4
   INTEGER :: na1, na2, viz, nt1, nt2, ldim1, ldim2, ldimb, nb1, nb2, viz_b
   INTEGER :: off, off1, off2, off3, eq_na2
   INTEGER :: ldim_std1, ldim_std2
   ! in the nwfcU ordering
   COMPLEX(DP), ALLOCATABLE :: nrg_nc(:,:,:,:,:,:)
   COMPLEX(DP) :: phase, nrgtmp, sign
   INTEGER, EXTERNAL :: find_viz
   !
   CALL start_clock('new_nsg')
   !
   ldim = 0
   DO nt = 1, ntyp
      ldim = MAX(ldim,ldim_u(nt))
   ENDDO
   !
   ALLOCATE ( nrg_nc (ldim, ldim, max_num_neighbors, nat, npol, npol) )
   nrg_nc(:,:,:,:,:,:) = (0.d0, 0.d0)
   nsgnew (:,:,:,:,:) = (0.d0, 0.d0)
   !
   IF (noncolin) THEN
      ALLOCATE (proj%k (nwfcU, nbnd))      
   ELSE        
      CALL allocate_bec_type ( nwfcU, nbnd, proj )
   ENDIF  
   !
   ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
   !
   ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
   !
   ! we start a loop on k points
   !
   DO ik = 1, nks
      !
      !
      npw = ngk(ik)
      !
      IF (nks > 1) CALL get_buffer (evc, nwordwfc, iunwfc, ik)
      !
      ! make the projection
      !
      IF ( Hubbard_projectors == 'pseudo' ) THEN
         CALL errore('new_nsg', 'Hubbard_projectors = pseudo is not supported',1)
      ELSE
         IF (nks > 1) CALL get_buffer (wfcU, nwordwfcU, iunhub, ik)
         CALL ZGEMM ('C', 'N', nwfcU, nbnd, npwx*npol, (1.0_DP, 0.0_DP), wfcU, &
               npwx*npol, evc, npwx*npol, (0.0_DP, 0.0_DP),  proj%k, nwfcU)
         CALL mp_sum( proj%k( :, 1:nbnd ), intra_bgrp_comm )
      ENDIF
      !
      ! compute the phase factors for this k and put the result 
      ! in the phase_fac array
      !
      CALL phase_factor(ik)
      !
      ! compute the occupation matrix (nsg_{I,m1,J,m2,s}) of the
      ! atomic orbitals
      !
      DO ibnd = 1, nbnd
         !
         DO na1 = 1, nat  
            ! 
            nt1 = ityp(na1)
            !
            IF (ldim_u(nt1).GT.0) THEN
               !
               ldim1 = ldim_u(nt1)
               !
               DO viz =1, neighood(na1)%num_neigh
                  !
                  na2 = neighood(na1)%neigh(viz)
                  eq_na2 = at_sc(na2)%at
                  nt2 = ityp(eq_na2)
                  ldim2 = ldim_u(nt2)
                  ! 
                  IF (na1.GT.na2) THEN
                     !
                     DO m1 = 1, ldim1
                        DO m2 = 1, ldim2
                           DO is1 = 1, npol
                              DO is2 = 1, npol
                                 nrg_nc(m2,m1,viz,na1,is2,is1) = &
                                          CONJG(nrg_nc(m1,m2,find_viz(na2,na1),na2,is1,is2))
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                     !
                  ELSE
                     !
                     phase = phase_fac(na2)
                     !
                     DO m1 = 1, ldim1
                        !
                        off1 = offsetU(na1) + m1
                        !
                        DO m2 = 1, ldim2
                           !
                           off2 = offsetU(eq_na2) + m2
                           !
                           DO is3 = 1, npol
                              DO is4 = 1, npol
                                 IF (domag.or.lspinorb) then
                                    nrg_nc(m2,m1,viz,na1,is4,is3) = &
                                    nrg_nc(m2,m1,viz,na1,is4,is3) + &
                                    dcmplx( wg(ibnd,ik) * ( proj%k(off1+ldim1*(is4-1),ibnd) * &
                                    CONJG(proj%k(off2+ldim2*(is3-1),ibnd) * phase) ) )
                                 ELSE 
                                    nrg_nc(m2,m1,viz,na1,is4,is3) = &
                                    nrg_nc(m2,m1,viz,na1,is4,is3) + &
                                    dble( wg(ibnd,ik) * ( proj%k(off1+ldim1*(is4-1),ibnd) * &
                                    CONJG(proj%k(off2+ldim2*(is3-1),ibnd) * phase) ) )
                                 ENDIF
                              ENDDO
                           ENDDO
                           !
                        ENDDO ! m2
                     ENDDO ! m1
                     !
                  ENDIF
                  !
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      !
   ENDDO ! ik
   !
   CALL deallocate_bec_type (proj) 
   !
   CALL mp_sum( nrg_nc, inter_pool_comm )
   !
   ! impose hermiticity of n_{m1,m2}
   !
   DO ibnd = 1,nbnd
      DO na1 = 1, nat 
         nt1 = ityp(na1)
         IF (ldim_u(nt1).GT.0) THEN
            ldim1 = ldim_u(nt1) 
            DO viz = 1, neighood(na1)%num_neigh
               na2 = neighood(na1)%neigh(viz)
               IF (na1.GT.na2) THEN
                  eq_na2 = at_sc(na2)%at
                  nt2 = ityp (eq_na2)
                  ldim2 = ldim_u(nt2)
                  DO is1 = 1, npol
                     DO is2 = 1, npol
                        DO m1 = 1, ldim1
                           DO m2 = 1, ldim2
                              nrg_nc(m2,m1,viz,na1,is2,is1) = &
                              ( nrg_nc(m2,m1,viz,na1,is2,is1) + &
                              CONJG(nrg_nc(m1,m2,find_viz(na2,na1),na2,is1,is2)) )*0.5d0
                              nrg_nc(m1,m2,find_viz(na2,na1),na2,is1,is2) =  &
                              CONJG(nrg_nc(m2,m1,viz,na1,is2,is1))
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ! 
   ! symmetrize the quantities nr -> ns
   !
   DO na1 = 1, nat
      !
      nt1 = ityp(na1)
      !
      IF (ldim_u(nt1).GT.0) THEN 
         !
         ldim_std1 = 2*Hubbard_l(nt1)+1
         !
         DO viz = 1, neighood(na1)%num_neigh
            !
            na2 = neighood(na1)%neigh(viz)
            eq_na2 = at_sc(na2)%at
            nt2 = ityp(eq_na2)
            ldim_std2 = 2*Hubbard_l(nt2)+1
            !
            IF (na1.GT.na2) THEN
               !
               ! we don't need to compute again
               DO is1 = 1, npol
                  DO is2 = 1, npol
                     is3 = npol*(is2-1) + is1
                     is4 = npol*(is1-1) + is2
                     DO m1 = 1, ldim_u(nt1)
                        DO m2 = 1, ldim_u(nt2)
                           nsgnew(m2,m1,viz,na1,is3) = &
                                 CONJG(nsgnew(m1,m2,find_viz(na2,na1),na2,is4))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
               !
            ELSE
               !
               DO m1 = 1, ldim_u(nt1)
                  !
                  off  = 1
                  off2 = 2*Hubbard_l(nt1)+1
                  !
                  DO m2 = 1, ldim_u(nt2)
                     !
                     off1 = 1
                     off3 = 2*Hubbard_l(nt2)+1
                     !
                     DO is1 = 1, npol
                        DO is2 = 1, npol
                           is = npol*(is1-1) + is2
                           !
                           ! Perform symmetrization using all available symmetries
                           !
                           !nsgnew(m2,m1,viz,na1,is) = nrg_nc(m2,m1,viz,na1,is2,is1)
                           !if(.false.) then
                           DO isym = 1, nsym
                              !
                              sign = (-1)**( Hubbard_l(nt1) + Hubbard_l(nt2) )
                              CALL symonpair(na1,na2,isym,nb1,nb2)
                              !
                              viz_b = find_viz(nb1,nb2)
                              ! 
                              DO m3 = off, off2
                                 DO m4 = off1, off3
                                    DO is3 = 1, npol
                                       DO is4 = 1, npol
                                          nrgtmp = t_rev(isym)*conjg(nrg_nc(m4,m3,viz_b,nb1,is4,is3))+ &
                                                (1-t_rev(isym))*nrg_nc(m4,m3,viz_b,nb1,is4,is3)
                                          IF (ll(m1,nt1).EQ.0 .AND. &
                                              ll(m2,nt2).EQ.0) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) )*& 
                                                      nrgtmp * &
                                                      d_spin_ldau(is2,is4,isym)/nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                                   ll(m2,nt2).EQ.0) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) )*&
                                                      d1(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d_spin_ldau(is2,is4,isym)/nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                                   ll(m2,nt2).EQ.0) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) )*&
                                                      d2(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d_spin_ldau(is2,is4,isym)/nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                                   ll(m2,nt2).EQ.0) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) )*&
                                                      d3(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d_spin_ldau(is2,is4,isym)/nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                                   ll(m2,nt2).EQ.1) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) )* &
                                                      nrgtmp * &
                                                      d_spin_ldau(is2,is4,isym) * &
                                                      d1(m4,m2,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                                   ll(m2,nt2).EQ.1) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d1(m1,m3,isym) * &
                                                      nrgtmp * &
                                                      d1(m2,m4,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                                   ll(m2,nt2).EQ.1) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d2(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d1(m4,m2,isym) *&
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                                   ll(m2,nt2).EQ.1) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d3(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d1(m4,m2,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                                   ll(m2,nt2).EQ.2) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      nrgtmp * &
                                                      d2(m4,m2,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                                   ll(m2,nt2).EQ.2) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d1(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d2(m4,m2,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                                   ll(m2,nt2).EQ.2) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d2(m1,m3,isym) * &
                                                      nrgtmp * &
                                                      d2(m2,m4,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                                   ll(m2,nt2).EQ.2) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d3(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d2(m4,m2,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                                   ll(m2,nt2).EQ.3) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      nrgtmp * &
                                                      d3(m4,m2,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !   
                                          ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                                   ll(m2,nt2).EQ.3) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d1(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d3(m4,m2,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                                   ll(m2,nt2).EQ.3) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d2(m3,m1,isym) * &
                                                      nrgtmp * &
                                                      d3(m4,m2,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                                   ll(m2,nt2).EQ.3) THEN
                                                   nsgnew(m2,m1,viz,na1,is) = &
                                                      nsgnew(m2,m1,viz,na1,is) +  &
                                                      CONJG( d_spin_ldau(is1,is3,isym) ) * &
                                                      d3(m1,m3,isym) * &
                                                      nrgtmp * &
                                                      d3(m2,m4,isym) * &
                                                      d_spin_ldau(is2,is4,isym) / nsym
                                                   !
                                          ELSE
                                             CALL errore ('new_nsg', &
                                                   'angular momentum not implemented for at least one type', &
                                                   ABS(Hubbard_l(nt1)) )
                                          END IF 
                                       ENDDO ! is4
                                       !
                                    ENDDO ! is3
                                    !
                                 ENDDO !m4
                                 !
                              ENDDO !m3
                              !
                           ENDDO !isym
                           !endif
                           !
                        ENDDO !is2
                        !
                     ENDDO !is1
                     !
                  ENDDO !m2
                  !
               ENDDO  !m1
               !
            ENDIF !na1 > na2
            ! 
         ENDDO !viz
         !
      ENDIF  !ldim_u > 0
      !
   ENDDO !na1
   !
   !write(stdout,*) "nsgnew", nsgnew
   DEALLOCATE(nrg_nc)
   !
   CALL stop_clock('new_nsg')
   !
   RETURN
   !
 END SUBROUTINE new_nsg_nc
 ! ------------------------------------------------
 
