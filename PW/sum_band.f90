!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE sum_band()
  !----------------------------------------------------------------------------
  !
  ! ... calculates the symmetrized charge density and sum of occupied
  ! ... eigenvalues.
  ! ... this version works also for metals (gaussian spreading technique)  
  !
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : gamma_only
  USE control_flags,        ONLY : diago_full_acc
  USE cell_base,            ONLY : at, bg, omega
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE ener,                 ONLY : eband, demet, ef, ef_up, ef_dw
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, &
                                   nkstot, wk, xk, nelec, nelup, neldw, &
                                   two_fermi_energies
  USE ktetra,               ONLY : ltetra, ntetra, tetra
  USE ldaU,                 ONLY : lda_plus_U
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : nsym, s, ftau
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param,           ONLY : nh, tvanp, nhm
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, bfield, npol
  USE spin_orb,             ONLY : lspinorb, domag, so, fcoef
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et, btype
  USE mp_global,            ONLY : intra_image_comm, me_image, &
                                   root_image, npool, my_pool_id
  USE mp,                   ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ikb, jkb, ijkb0, ih, jh, ijh, na, np
    ! counters on beta functions, atoms, pseudopotentials  
  INTEGER :: ir, is, ig, ibnd, ik
    ! counter on 3D r points
    ! counter on spin polarizations
    ! counter on g vectors
    ! counter on bands
    ! counter on k points  
  real (DP) demet_up, demet_dw
  !
  !
  CALL start_clock( 'sum_band' )
  !
  becsum(:,:,:) = 0.D0
  rho(:,:)      = 0.D0
  eband         = 0.D0
  demet         = 0.D0
  !
  IF ( .NOT. lgauss .AND. .NOT. ltetra .AND. .NOT. tfixed_occ ) THEN
     !
     ! ... calculate weights for the insulator case
     !
     IF ( two_fermi_energies ) THEN
        !
        CALL iweights( nks, wk, nbnd, nelup, et, ef_up, wg, 1, isk )
        CALL iweights( nks, wk, nbnd, neldw, et, ef_dw, wg, 2, isk )
        !
     ELSE
        !
        CALL iweights( nks, wk, nbnd, nelec, et, ef,    wg, 0, isk )
        !
     END IF
     !
  ELSE IF ( ltetra ) THEN
     !
     ! ... calculate weights for the metallic case
     !
     CALL poolrecover( et, nbnd, nkstot, nks )
     !
     IF ( me_image == root_image ) THEN
        !
        IF (two_fermi_energies) THEN
           !
           CALL tweights( nkstot, nspin, nbnd, nelup, &
                          ntetra, tetra, et, ef_up, wg, 1, isk )
           CALL tweights( nkstot, nspin, nbnd, neldw, &
                          ntetra, tetra, et, ef_dw, wg, 2, isk )
           !
        ELSE
           !
           CALL tweights( nkstot, nspin, nbnd, nelec, &
                          ntetra, tetra, et, ef, wg, 0, isk )
           !
        END IF
       !
     END IF
     !
     CALL poolscatter( nbnd, nkstot, wg, nks, wg )
     !
     CALL mp_bcast( ef, root_image, intra_image_comm )
     !
  ELSE IF ( lgauss ) THEN
     !
     IF ( two_fermi_energies ) THEN
        !
        CALL gweights( nks, wk, nbnd, nelup, degauss, &
                       ngauss, et, ef_up, demet_up, wg, 1, isk )
        CALL gweights( nks, wk, nbnd, neldw, degauss, &
                       ngauss, et, ef_dw, demet_dw, wg, 2, isk )
        !
        demet = demet_up + demet_dw
        !
        bfield(3) = 0.5D0*( ef_up - ef_dw )
        !
     ELSE
        !
        CALL gweights( nks, wk, nbnd, nelec, degauss, &
                       ngauss, et, ef, demet, wg, 0, isk)
     END IF
     !
  ELSE IF ( tfixed_occ ) THEN
     !
     IF ( npool == 1 ) THEN
        !
        wg = f_inp
        !
     ELSE
        !
        wg(:,1) = f_inp(:,my_pool_id+1)
        wg(:,2) = f_inp(:,my_pool_id+1)
        !
     END IF
     !
     ef = - 1.D+20
     !
     DO is = 1, nspin
        !
        DO ibnd = 1, nbnd
           !
           IF ( wg(ibnd,is) > 0.D0 ) ef = MAX( ef, et(ibnd,is) )
           !
        END DO
        !
     END DO
     !
  END IF
  !
  IF ( diago_full_acc ) THEN
     !
     ! ... for diagonalization purposes all the bands are considered occupied
     !
     btype(:,:) = 1
     !
  ELSE
     !
     ! ... for diagonalization purposes a band is considered empty when its 
     ! ... occupation is less than 1.0 %
     !
     btype(:,:) = 1
     !
     FORALL( ik = 1:nks, wk(ik) > 0.D0 ) 
        !
        WHERE( wg(:,ik) / wk(ik) < 0.01D0 ) btype(:,ik) = 0
        !
     END FORALL
     !
  END IF
  !
  ! ... Needed for LDA+U
  !
  IF ( lda_plus_u ) CALL new_ns()  
  !     
  ! ... specific routines are called to sum for each k point the contribution
  ! ... of the wavefunctions to the charge
  !
  IF ( gamma_only ) THEN
     !
     CALL sum_band_gamma()
     !
  ELSE
     !
     CALL sum_band_k()
     !
  END IF    
  !
  ! ... If a double grid is used, interpolate onto the fine grid
  !
  IF ( doublegrid ) THEN
     !
     DO is = 1, nspin
        !
        CALL interpolate( rho(1,is), rho(1,is), 1 )
        !
     END DO
     !
  END IF
  !
  ! ... Here we add the Ultrasoft contribution to the charge
  !
  IF ( okvan ) CALL addusdens()
  !
  IF ( noncolin .AND. .NOT. domag ) rho(:,2:4)=0.D0
  !
  CALL poolreduce( 1, eband )
  CALL poolreduce( 1, demet )
  !
  ! ... symmetrization of the charge density (and local magnetization)
  !
#if defined (__PARA)
  !
  ! ... reduce charge density across pools
  !
  CALL poolreduce( nspin * nrxx, rho )
  !
  IF ( noncolin ) THEN
     !
     CALL psymrho( rho(1,1), nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
     !
     IF ( domag ) &
        CALL psymrho_mag( rho(1,2), nrx1, nrx2, nrx3, &
                          nr1, nr2, nr3, nsym, s, ftau, bg, at )
     !
  ELSE
     !
     DO is = 1, nspin
        !
        CALL psymrho( rho(1,is), nrx1, nrx2, nrx3, &
                      nr1, nr2, nr3, nsym, s, ftau )
        !
     END DO
     !
  END IF
  !
#else
  !
  IF ( noncolin ) THEN
     !
     CALL symrho( rho(1,1), nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
     !
     IF ( domag ) &
        CALL symrho_mag( rho(1,2), nrx1, nrx2, nrx3, &
                         nr1, nr2, nr3, nsym, s, ftau, bg, at )
     !
  ELSE
     !
     DO is = 1, nspin
        !
        CALL symrho( rho(1,is), nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
        !
     END DO
     !
  END IF
  !
#endif
  !  
  CALL stop_clock( 'sum_band' )      
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_gamma()
       !-----------------------------------------------------------------------
       !
       ! ... gamma version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       REAL(DP), ALLOCATABLE :: becp(:,:)
         ! contains <beta|psi>
       !
       !
       ALLOCATE( becp( nkb, nbnd ) )
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       IF ( nks > 1 ) REWIND( iunigk )
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          IF ( nks > 1 ) THEN
             !
             READ( iunigk ) npw, igk
             CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
             !
          END IF
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = 1, nbnd
             !
             ! ... the sum of eband and demet is the integral for  
             ! ... e < ef of e n(e) which reduces for degauss=0 to the sum of 
             ! ... the eigenvalues.
             !
             eband = eband + et(ibnd,ik) * wg(ibnd,ik)
             !
          END DO
          !
          DO ibnd = 1, nbnd, 2
             !
             psic(:) = ( 0.D0, 0.D0 )
             !
             IF ( ibnd < nbnd ) THEN
                !
                ! ... two ffts at the same time
                !
                psic(nls(igk(1:npw)))  = evc(1:npw,ibnd) + &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1)
                psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ibnd) - &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                !
             ELSE
                !
                psic(nls(igk(1:npw)))  = evc(1:npw,ibnd)
                psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ibnd) )
                !
             END IF
             !
             CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
             !
             w1 = wg(ibnd,ik) / omega
             !
             ! ... increment the charge density ...
             !
             IF ( ibnd < nbnd ) THEN
                !
                ! ... two ffts at the same time
                !
                w2 = wg(ibnd+1,ik) / omega
                !
             ELSE
                !
                w2 = w1
                !
             END IF
             !
             DO ir = 1, nrxxs
                !
                rho(ir,current_spin) = rho(ir,current_spin) + &
                                       w1 *  DBLE( psic(ir) )**2 + &
                                       w2 * AIMAG( psic(ir) )**2
                !
             END DO
             !
          END DO
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( .NOT. okvan ) CYCLE k_loop
          !
          IF ( nkb > 0 ) &
             CALL ccalbec( nkb, npwx, npw, nbnd, becp, vkb, evc )
          !
          CALL start_clock( 'becsum' )
          !
          DO ibnd = 1, nbnd
             !
             w1 = wg(ibnd,ik)
             ijkb0 = 0
             !
             DO np = 1, ntyp
                !
                IF ( tvanp(np) ) THEN
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) THEN
                         !
                         ijh = 1
                         !
                         DO ih = 1, nh(np)
                            !
                            ikb = ijkb0 + ih
                            !
                            becsum(ijh,na,current_spin) = &
                                            becsum(ijh,na,current_spin) + &
                                            w1 * becp(ikb,ibnd) * becp(ikb,ibnd)
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + &
                                     w1 * 2.D0 * becp(ikb,ibnd) * becp(jkb,ibnd)
                               !
                               ijh = ijh + 1
                               !
                            END DO
                            !
                         END DO
                         !
                         ijkb0 = ijkb0 + nh(np)
                         !
                      END IF
                      !
                   END DO
                   !
                ELSE
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                      !
                   END DO
                   !
                END IF
                !
             END DO
             !
          END DO
          !
          CALL stop_clock( 'becsum' )
          !
       END DO k_loop
       !
       DEALLOCATE( becp )
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       COMPLEX(DP), ALLOCATABLE :: becp(:,:), becp_nc(:,:,:)
       ! contains <beta|psi>
       !
       COMPLEX(DP), ALLOCATABLE :: be1(:,:), be2(:,:)
       !
       INTEGER :: ipol, kh, kkb, is1, is2
       !
       IF (noncolin) THEN
          ALLOCATE( becp_nc( nkb, npol, nbnd ) )
          IF (lspinorb) ALLOCATE(be1(nhm,2), be2(nhm,2))
       ELSE
          ALLOCATE( becp( nkb, nbnd ) )
       ENDIF
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       IF ( nks > 1 ) REWIND( iunigk )
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          IF ( nks > 1 ) THEN
             !
             READ( iunigk ) npw, igk
             CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
             !
          END IF
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = 1, nbnd
             !
             eband = eband + et(ibnd,ik) * wg(ibnd,ik)
             !
             ! ... the sum of eband and demet is the integral for e < ef of 
             ! ... e n(e) which reduces for degauss=0 to the sum of the 
             ! ... eigenvalues 
             w1 = wg(ibnd,ik) / omega
             IF (noncolin) THEN
                psic_nc = (0.D0,0.D0)
                DO ig = 1, npw
                   psic_nc(nls(igk(ig)),1)=evc(ig     ,ibnd)
                   psic_nc(nls(igk(ig)),2)=evc(ig+npwx,ibnd)
                END DO
                call cft3s (psic_nc(1,1), nr1s, nr2s, nr3s, &
                                          nrx1s,nrx2s,nrx3s, 2)
                call cft3s (psic_nc(1,2), nr1s, nr2s, nr3s, &
                                          nrx1s,nrx2s,nrx3s, 2)
                !
                ! increment the charge density ...
                !
                DO ipol=1,npol
                   DO ir = 1, nrxxs
                      rho (ir, 1) = rho (ir, 1) + &
                      w1*( DBLE(psic_nc(ir,ipol))**2+AIMAG(psic_nc(ir,ipol))**2)
                   END DO
                END DO
                !
                ! In this case, calculate also the three
                ! components of the magnetization (stored in rho(ir,2-4) )
                !
                IF (domag) THEN
                   DO ir = 1,nrxxs
                      rho(ir,2) = rho(ir,2) + w1*2.D0* &
                         (DBLE(psic_nc(ir,1))* DBLE(psic_nc(ir,2)) + &
                         AIMAG(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)))

                      rho(ir,3) = rho(ir,3) + w1*2.D0* &
                         (DBLE(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)) - &
                          DBLE(psic_nc(ir,2))*AIMAG(psic_nc(ir,1)))

                      rho(ir,4) = rho(ir,4) + w1* &
                         (DBLE(psic_nc(ir,1))**2+AIMAG(psic_nc(ir,1))**2 &
                         -DBLE(psic_nc(ir,2))**2-AIMAG(psic_nc(ir,2))**2)
                   END DO
                ELSE
                   rho(:,2:4)=0.d0
                END IF
                !
             ELSE
                !
                psic(:) = ( 0.D0, 0.D0 )
                !
                psic(nls(igk(1:npw))) = evc(1:npw,ibnd)
                !
                CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
                !
                !
                ! ... increment the charge density ...
                !
                DO ir = 1, nrxxs
                   !
                   rho(ir,current_spin) = rho(ir,current_spin) + &
                                            w1 * ( DBLE( psic(ir) )**2 + &
                                                  AIMAG( psic(ir) )**2 )
                   !
                END DO
                !
             END IF
             !
          END DO
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( .NOT. okvan ) CYCLE k_loop
          !
          IF (noncolin) THEN
             IF ( nkb > 0 ) &
                CALL ccalbec_nc( nkb, npwx, npw, npol, nbnd, &
                                                 becp_nc, vkb, evc )
          ELSE
             IF ( nkb > 0 ) &
                CALL ccalbec( nkb, npwx, npw, nbnd, becp, vkb, evc )
          ENDIF
          !
          CALL start_clock( 'becsum' )
          !
          DO ibnd = 1, nbnd
             !
             w1 = wg(ibnd,ik)
             ijkb0 = 0
             !
             DO np = 1, ntyp
                !
                IF ( tvanp(np) ) THEN
                   !
                   DO na = 1, nat
                      !
                      IF (ityp(na)==np) THEN
                         !
                         IF (so(np)) THEN
                            be1=(0.d0,0.d0)
                            be2=(0.d0,0.d0)
                            DO ih = 1, nh(np)
                               ikb = ijkb0 + ih
                               DO kh = 1, nh(np)
                                  IF ((nhtol(kh,np)==nhtol(ih,np)).AND. &
                                  (nhtoj(kh,np)==nhtoj(ih,np)).AND.     &
                                  (indv(kh,np)==indv(ih,np))) THEN
                                     kkb=ijkb0 + kh
                                     DO is1=1,2
                                        DO is2=1,2
                                           be1(ih,is1)=be1(ih,is1)+  &
                                               fcoef(ih,kh,is1,is2,np)*  &
                                                    becp_nc(kkb,is2,ibnd)
                                           be2(ih,is1)=be2(ih,is1)+ &
                                               fcoef(kh,ih,is2,is1,np)* &
                                               CONJG(becp_nc(kkb,is2,ibnd))
                                        END DO
                                     END DO
                                  END IF
                               END DO
                            END DO
                         END IF
                         ijh = 1
                         !
                         DO ih = 1, nh(np)
                            !
                            ikb = ijkb0 + ih
                            !
                            IF (noncolin) THEN
                               !
                               IF (so(np)) THEN
                                  becsum(ijh,na,1)=becsum(ijh,na,1)+ w1*&
                                     (be1(ih,1)*be2(ih,1)+ be1(ih,2)*be2(ih,2))
                                  IF (domag) THEN
                                     becsum(ijh,na,2)=becsum(ijh,na,2)+ w1*&
                                     (be1(ih,2)*be2(ih,1)+ be1(ih,1)*be2(ih,2))
                                     becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                                               w1*(0.d0,-1.d0)*      &  
                                     (be1(ih,2)*be2(ih,1)-be1(ih,1)*be2(ih,2))
                                     becsum(ijh,na,4)=becsum(ijh,na,4)+ w1* &
                                     (be1(ih,1)*be2(ih,1)-be1(ih,2)*be2(ih,2))
                                  ENDIF
                               ELSE
                                  becsum(ijh,na,1) = becsum(ijh,na,1)   &
                                    + w1*( CONJG(becp_nc(ikb,1,ibnd))   &
                                                *becp_nc(ikb,1,ibnd)    &
                                    +      CONJG(becp_nc(ikb,2,ibnd))   &
                                                *becp_nc(ikb,2,ibnd) )
                                  IF (domag) THEN
                                     becsum(ijh,na,2)=becsum(ijh,na,2)  &
                                     + w1*(CONJG(becp_nc(ikb,2,ibnd))   &
                                                 *becp_nc(ikb,1,ibnd)   &
                                     +     CONJG(becp_nc(ikb,1,ibnd))   &
                                                 *becp_nc(ikb,2,ibnd) )
                                     becsum(ijh,na,3)=becsum(ijh,na,3)  &
                                          + w1*2.d0     &
                                      *AIMAG(CONJG(becp_nc(ikb,1,ibnd))* &
                                                   becp_nc(ikb,2,ibnd) )
                                     becsum(ijh,na,4) = becsum(ijh,na,4)    &
                                          + w1*( CONJG(becp_nc(ikb,1,ibnd)) &
                                                      *becp_nc(ikb,1,ibnd)  &
                                          -      CONJG(becp_nc(ikb,2,ibnd)) &
                                                      *becp_nc(ikb,2,ibnd) )
                                  END IF
                               END IF
                            ELSE
                               becsum(ijh,na,current_spin) = &
                                        becsum(ijh,na,current_spin) + &
                                        w1 * DBLE( CONJG( becp(ikb,ibnd) ) * &
                                                          becp(ikb,ibnd) )
                            END IF                       
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               IF (noncolin) THEN
                                  IF (so(np)) THEN
                                     becsum(ijh,na,1)=becsum(ijh,na,1)+ w1*(  &
                                   (be1(jh,1)*be2(ih,1)+be1(jh,2)*be2(ih,2))+ &
                                   (be1(ih,1)*be2(jh,1)+be1(ih,2)*be2(jh,2)))
                                     IF (domag) THEN
                                       becsum(ijh,na,2)=becsum(ijh,na,2)+w1*( &
                                     (be1(jh,2)*be2(ih,1)+be1(jh,1)*be2(ih,2))+&
                                     (be1(ih,2)*be2(jh,1)+be1(ih,1)*be2(jh,2)))
                                       becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                                          w1*(0.d0,-1.d0)*((be1(jh,2)*&
                                          be2(ih,1)-be1(jh,1)*be2(ih,2))+ &
                                         (be1(ih,2)*be2(jh,1)-be1(ih,1)*&
                                                    be2(jh,2)) )
                                       becsum(ijh,na,4)=becsum(ijh,na,4)+ &
                                              w1*((be1(jh,1)*be2(ih,1)- &
                                             be1(jh,2)*be2(ih,2))+  &
                                             (be1(ih,1)*be2(jh,1)-  &
                                              be1(ih,2)*be2(jh,2)) )
                                     END IF
                                  ELSE
                                     becsum(ijh,na,1)= becsum(ijh,na,1)+ &
                                                      w1*2.d0* &
                                     DBLE(CONJG(becp_nc(ikb,1,ibnd))* &
                                                becp_nc(jkb,1,ibnd) + &
                                          CONJG(becp_nc(ikb,2,ibnd))* &
                                                becp_nc(jkb,2,ibnd) )
                                     IF (domag) THEN
                                        becsum(ijh,na,2)=becsum(ijh,na,2)+ &
                                                          w1*2.d0* &
                                           DBLE(CONJG(becp_nc(ikb,2,ibnd))* &
                                                      becp_nc(jkb,1,ibnd) + &
                                                CONJG(becp_nc(ikb,1,ibnd))* &
                                                      becp_nc(jkb,2,ibnd) )
                                        becsum(ijh,na,3)=becsum(ijh,na,3)+ &
                                                       w1*2.d0* &
                                            AIMAG(CONJG(becp_nc(ikb,1,ibnd))* &
                                                        becp_nc(jkb,2,ibnd) + &
                                                  CONJG(becp_nc(ikb,1,ibnd))* &
                                                        becp_nc(jkb,2,ibnd) )
                                        becsum(ijh,na,4)=becsum(ijh,na,4)+ &
                                                       w1*2.d0* &
                                            DBLE(CONJG(becp_nc(ikb,1,ibnd))* &
                                                       becp_nc(jkb,1,ibnd) - &
                                                 CONJG(becp_nc(ikb,2,ibnd))* &
                                                       becp_nc(jkb,2,ibnd) )
                                     END IF
                                  END IF
                               ELSE
                               !
                                   becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + w1 * 2.D0 * &
                                     DBLE( CONJG( becp(ikb,ibnd) ) * &
                                           becp(jkb,ibnd) )
                               ENDIF
                               !            
                               ijh = ijh + 1
                               !
                            END DO
                            !
                         END DO
                         !
                         ijkb0 = ijkb0 + nh(np)
                         !
                      END IF
                      !
                   END DO
                   !
                ELSE
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                      !
                   END DO
                   !
                END IF
                !
             END DO
             !
          END DO
          !
          CALL stop_clock( 'becsum' )
          !
       END DO k_loop
       !
       IF (noncolin) THEN
          DEALLOCATE( becp_nc )
          IF (lspinorb) DEALLOCATE(be1, be2)
       ELSE
          DEALLOCATE( becp )
       ENDIF
       !
       RETURN
       !
     END SUBROUTINE sum_band_k
     !
END SUBROUTINE sum_band
