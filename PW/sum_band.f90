!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE sum_band()
  !----------------------------------------------------------------------------
  !
  ! ... calculates the symmetrized charge density and sum of occupied
  ! ... eigenvalues.
  ! ... this version works also for metals (gaussian spreading technique)  
  !
  USE kinds,           ONLY : DP
  USE wvfct,                ONLY : gamma_only
  USE brilz,                ONLY : omega
  USE basis,                ONLY : nat, ntyp, ityp
  USE ener,                 ONLY : eband, demet, ef
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, &
                                   nkstot, wk, xk, nelec
  USE ktetra,               ONLY : ltetra, ntetra, tetra
  USE ldaU,                 ONLY : lda_plus_U
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : nsym, s, ftau
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE us,                   ONLY : okvan, tvanp, becsum, nh, nkb, vkb
  USE wavefunctions_module, ONLY : evc, psic
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
#ifdef __PARA
  USE para,                 ONLY : me, mypool
#endif
  !
  IMPLICIT NONE
  !
  !
  CALL start_clock( 'sum_band' )
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
  CALL stop_clock( 'sum_band' )      
  !
  RETURN
  !
  CONTAINS
     !
     !--------------------------------------------------------------------------
     SUBROUTINE sum_band_gamma()
       !--------------------------------------------------------------------------
       !
       ! ... gamma version
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
       REAL(KIND=DP) :: w1, w2
         ! weight
       REAL(KIND=DP), ALLOCATABLE :: becp(:,:)
         ! contains <beta|psi>
       !
       !
       ALLOCATE( becp( nkb, nbnd ) )
       !
       becsum(:,:,:) = 0.D0
       rho(:,:)      = 0.D0
       eband         = 0.D0
       demet         = 0.D0
       !
       ! ... calculate weights for the insulator case
       !
       IF ( .NOT. lgauss .AND. .NOT. ltetra .AND. .NOT. tfixed_occ ) THEN
          !
          CALL iweights( nks, wk, nbnd, nelec, et, ef, wg )
          !
          ! ... calculate weights for the metallic case
          !
       ELSE IF ( ltetra ) THEN
          !
#ifdef __PARA
          CALL poolrecover( et, nbnd, nkstot, nks )
          !
          IF ( me == 1 .AND. mypool == 1 ) THEN
             !
#endif
             CALL tweights( nkstot, nspin, nbnd, nelec, ntetra, &
                            tetra, et, ef, wg )
#ifdef __PARA
             !
          END IF
          !
          CALL poolscatter( nbnd, nkstot, wg, nks, wg )
          !
          IF ( me == 1 ) CALL poolbcast( 1, ef )
          !
          CALL broadcast( 1, ef )
          !
#endif
       ELSE IF ( lgauss ) THEN
          !
          CALL gweights( nks, wk, nbnd, nelec, &
                         degauss, ngauss, et, ef, demet, wg )
          !
       ELSE IF ( tfixed_occ ) THEN
          !
          ef = - 1.0D+20
          !
          DO is = 1, nspin
             !
             DO ibnd = 1, nbnd
                !
                wg(ibnd,is) = f_inp(ibnd,is)
                !
                IF ( wg(ibnd,is) > 0.D0 ) ef = MAX( ef, et(ibnd,is) )
                !
             END DO
             !
          END DO
          !
       END IF
       !
       ! ... Needed for LDA+U
       !
       IF ( lda_plus_u ) CALL new_ns()
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
             ! ... the sum of eband and demet is the integral for  
             ! ... e < ef of e n(e) which reduces for degauss=0 to the sum of 
             ! ... the eigenvalues.
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
                                                   w1 * REAL( psic(ir) )**2 + &
                                                   w2 * IMAG( psic(ir) )**2
                !
             END DO
             !
          END DO
          !
          ! ... If we have a US pseudopotential we compute here the sumbec term
          !
          IF ( .NOT. okvan ) CYCLE k_loop
          !
          IF ( nkb > 0 ) &
             CALL pw_gemm( 'Y', nkb, nbnd, npw, &
                           vkb, npwx, evc, npwx, becp, nkb )
          !
          CALL start_clock( 'sumbec' )
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
          CALL stop_clock( 'sumbec' )
          !
       END DO k_loop
       !
       DEALLOCATE( becp )
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
#ifdef __PARA
       CALL poolreduce( 1, eband )
       CALL poolreduce( 1, demet )
#endif
       !
       ! ... symmetrization of the charge density (and local magnetization)
       !
#ifdef __PARA
       !
       ! ... reduce charge density across pools
       !
       CALL poolreduce( nspin * nrxx, rho )
       !
       DO is = 1, nspin
          !
          CALL psymrho( rho(1,is), nrx1, nrx2, nrx3, &
                        nr1, nr2, nr3, nsym, s, ftau )
          !
       END DO
#else
       DO is = 1, nspin
          !
          CALL symrho( rho(1,is), nrx1, nrx2, nrx3, &
                       nr1, nr2, nr3, nsym, s, ftau )
          !
       END DO
#endif
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma
     !
     !
     !--------------------------------------------------------------------
     SUBROUTINE sum_band_k()
       !--------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE becmod,  ONLY : becp
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
       REAL(KIND=DP) :: w1
         ! weight
       !
       !
       becsum(:,:,:) = 0.D0
       rho(:,:)      = 0.D0
       eband         = 0.D0
       demet         = 0.D0
       !
       ! ... calculate weights for the insulator case
       !
       IF ( .NOT. lgauss .AND. .NOT. ltetra .AND. .NOT. tfixed_occ ) THEN
          !
          CALL iweights( nks, wk, nbnd, nelec, et, ef, wg )
          !
          ! ... calculate weights for the metallic case
          !
       ELSE IF ( ltetra ) THEN
#ifdef __PARA
          !
          CALL poolrecover( et, nbnd, nkstot, nks )
          !
          IF ( me == 1 .AND. mypool == 1 ) THEN
             !
#endif
             CALL tweights( nkstot, nspin, nbnd, &
                            nelec, ntetra, tetra, et, ef, wg )
#ifdef __PARA
             !
          ENDIF
          !
          CALL poolscatter( nbnd, nkstot, wg, nks, wg )
          !
          IF ( me == 1 ) CALL poolbcast( 1, ef )
          !
          CALL broadcast( 1, ef )
#endif
          !
       ELSE IF ( lgauss ) THEN
          !
          CALL gweights( nks, wk, nbnd, nelec, &
                         degauss, ngauss, et, ef, demet, wg )
          !
       ELSE IF ( tfixed_occ ) THEN
          !
          ef = - 1.0D+20
          !
          DO is = 1, nspin
             !
             DO ibnd = 1, nbnd
                !
                wg(ibnd,is) = f_inp(ibnd,is)
                !
                IF ( wg(ibnd,is) > 0.D0 ) ef = MAX( ef, et(ibnd,is) )
               !
             END DO
             !
          END DO
          !
       END IF
       !
       ! ... Needed for LDA+U
       !
       IF ( lda_plus_u ) CALL new_ns()
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
             ! ... the factor two is for spin degeneracy
             !
             psic(:) = ( 0.D0, 0.D0 )
             !
             psic(nls(igk(1:npw))) = evc(1:npw,ibnd)
             !
             CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
             !
             w1 = wg(ibnd,ik) / omega
             !
             ! ... increment the charge density ...
             !
             DO ir = 1, nrxxs
                !
                rho(ir,current_spin) = rho(ir,current_spin) + &
                                                 w1 * ( REAL( psic(ir) )**2 + &
                                                        IMAG( psic(ir) )**2 )
                !
             END DO
             !
          END DO
          !
          ! ... If we have a US pseudopotential we compute here the sumbec term
          !
          IF ( .NOT. okvan ) CYCLE k_loop
          !
          IF ( nkb > 0 ) &
             CALL ccalbec( nkb, npwx, npw, nbnd, becp, vkb, evc )
          !
          CALL start_clock( 'sumbec' )
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
                                          w1 * REAL( CONJG( becp(ikb,ibnd) ) * &
                                                     becp(ikb,ibnd) )
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + w1 * 2.D0 * &
                                     REAL( CONJG( becp(ikb,ibnd) ) * &
                                           becp(jkb,ibnd) )
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
          CALL stop_clock( 'sumbec' )
          !
       END DO k_loop
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
#ifdef __PARA
       CALL poolreduce( 1, eband )
       CALL poolreduce( 1, demet )
#endif
       !
       ! ... symmetrization of the charge density (and local magnetization)
       !
#ifdef __PARA
       !
       ! ... reduce charge density across pools
       !
       CALL poolreduce( nspin * nrxx, rho )
       !
       DO is = 1, nspin
          !
          CALL psymrho( rho(1,is), nrx1, nrx2, nrx3, &
                        nr1, nr2, nr3, nsym, s, ftau )
          !
       END DO
#else
       DO is = 1, nspin
          !
          CALL symrho( rho(1,is), nrx1, nrx2, nrx3, &
                       nr1, nr2, nr3, nsym, s, ftau )
          !             
       END DO
#endif
       !
       RETURN
       !
     END SUBROUTINE sum_band_k
     !
END SUBROUTINE sum_band
