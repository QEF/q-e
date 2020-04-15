!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
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
  USE ldaU,                 ONLY : Hubbard_l, q_ae, U_projection, wfcU, nwfcU,     &
                                   ldim_u, ll, neighood, at_sc, nsgnew, phase_fac, &
                                   max_num_neighbors, Hubbard_l_back, backall,     &
                                   offsetU, offsetU_back, offsetU_back1, is_hubbard_back
  USE symm_base,            ONLY : d1, d2, d3
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE symm_base,            ONLY : nsym, irt
  USE wvfct,                ONLY : nbnd, npw, npwx, wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunsat, iunhub
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  !
  IMPLICIT NONE
  !
  TYPE (bec_type) :: proj     
  ! proj(nwfcU,nbnd)
  INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, m1, m2, m11, m22, m0, m00, ldim,i_type
  INTEGER :: na1, na2, viz, nt1, nt2, ldim1, ldim2, ldimb, nb1, nb2, viz_b
  INTEGER :: off, off1, off2, off3, eq_na2
  INTEGER :: ldim_std1, ldim_std2
  ! in the nwfcU ordering
  COMPLEX(DP), ALLOCATABLE :: nrg(:,:,:,:,:)
  COMPLEX(DP) :: phase
  INTEGER, EXTERNAL :: find_viz, type_interaction
  !
  CALL start_clock('new_nsg')
  !
  ldim = 0
  DO nt = 1, ntyp
     ldim = MAX(ldim,ldim_u(nt))
  ENDDO
  !
  ALLOCATE ( nrg (ldim, ldim, max_num_neighbors, nat, nspin) )
  nrg    (:,:,:,:,:) = (0.d0, 0.d0)
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
     IF ( U_projection == 'pseudo' ) THEN
        CALL errore('new_nsg', 'U_projection = pseudo is not supported',1)
        !CALL compute_pproj( ik, q_ae, proj )
     ELSE
        IF (nks > 1) CALL get_buffer (wfcU, nwordwfcU, iunhub, ik)
        CALL calbec ( npw, wfcU, evc, proj )
     END IF
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
                               m1.GT.2*(Hubbard_l(nt1)+Hubbard_l_back(nt1)+1)) &
                             off1 = offsetU_back1(na1) + m1 - &
                                    2*(Hubbard_l(nt1)+Hubbard_l_back(nt1)+1)
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
                                  m2.GT.2*(Hubbard_l(nt2)+Hubbard_l_back(nt2)+1)) &
                                off2 = offsetU_back1(eq_na2) + m2 - &
                                       2*(Hubbard_l(nt2)+Hubbard_l_back(nt2)+1)
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
                       off2 = ldim_std1 + 2*Hubbard_l_back(nt1) + 1 
                       !
                       IF ( backall(nt1) .AND. m1.GT.(ldim_std1+2*Hubbard_l_back(nt1)+1) ) THEN
                          off  = ldim_std1 + 2*Hubbard_l_back(nt1) + 2
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
                          off3 = ldim_std2 + 2*Hubbard_l_back(nt2) + 1
                          !
                          IF ( backall(nt2) .AND. m2.GT.(ldim_std2+2*Hubbard_l_back(nt2)+1) ) THEN
                             off1 = ldim_std2 + 2*Hubbard_l_back(nt2) + 2
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
                          ! 
                          DO m0 = off, off2
                             DO m00 = off1, off3
                                IF (ll(m1,nt1).EQ.0 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.0) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        nrg(m00,m0,viz_b,nb1,is) / nsym
                                ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.0) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d1(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        1.d0 / nsym
                                ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.0) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d2(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        1.d0 / nsym
                                ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.0) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d3(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        1.d0 / nsym

                                ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.1) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d1(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.1) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d1(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d1(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.1) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d2(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d1(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.1) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d3(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d1(m00-off1+1,m2-off1+1,isym) / nsym

                                ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.2) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d2(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.2) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d1(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d2(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.2) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d2(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d2(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.2) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d3(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d2(m00-off1+1,m2-off1+1,isym) / nsym

                                ELSE IF (ll(m1,nt1).EQ.0 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.3) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d3(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.1 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.3) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d1(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d3(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.2 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                  ll(m2,nt2).EQ.3) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d2(m0-off+1,m1-off+1,isym) * &
                                        nrg(m00,m0,viz_b,nb1,is) * &
                                        d3(m00-off1+1,m2-off1+1,isym) / nsym
                                ELSE IF (ll(m1,nt1).EQ.3 .AND. &
                                   m1.ge.off.and.m1.le.off2.and. &
                                   m2.ge.off1.and.m2.le.off3.and. &
                                     ll(m2,nt2).EQ.3) THEN
                                   nsgnew(m2,m1,viz,na1,is) = &
                                        nsgnew(m2,m1,viz,na1,is) +  &
                                        d3(m0-off+1,m1-off+1,isym) * &
                                     nrg(m00,m0,viz_b,nb1,is) * &
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
