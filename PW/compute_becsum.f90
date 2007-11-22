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
SUBROUTINE compute_becsum(iflag)
  !----------------------------------------------------------------------------
  !
  ! ... calculates the becsum term
  ! ... this version works also for metals (gaussian spreading technique)  
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : tpiba2
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk
  USE gvect
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, okvan
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : so
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, g2kin
  USE paw_onecenter,        ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iflag ! if 1 compute also the weights
  !
  ! ... local variables
  !
  INTEGER :: ikb, jkb, ijkb0, ih, jh, ijh, na, np
    ! counters on beta functions, atoms, pseudopotentials  
  INTEGER :: ir, is, ig, ibnd, ik, j
    ! counter on 3D r points
    ! counter on spin polarizations
    ! counter on g vectors
    ! counter on bands
    ! counter on k points  
  !
  !
  CALL start_clock( 'compute_becsum' )
  !
  becsum(:,:,:) = 0.D0
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  !
  IF (iflag==1) CALL weights ( )
  !
  IF (gamma_only) THEN
     CALL compute_becsum_gamma()
  ELSE
     CALL compute_becsum_k()
  ENDIF
  ! ... Needed for PAW: becsum has to be symmetrized so that they reflect a real integral
  ! in k-space, not only on the irreducible zone. For USPP there is no need to do this as
  ! becsums are only used to compute the density, which is symmetrized later.
  !
  IF ( okpaw ) CALL PAW_symmetrize(becsum)
  !
  CALL stop_clock( 'compute_becsum' )      
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_becsum_gamma()
       !-----------------------------------------------------------------------
       !
       ! ... gamma version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
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
          npw = ngk(ik)
          !
          IF ( nks > 1 ) THEN
             !
             READ( iunigk ) igk
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             !
          END IF
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
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
                IF ( upf(np)%tvanp ) THEN
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
     END SUBROUTINE compute_becsum_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_becsum_k()
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
       COMPLEX(DP), ALLOCATABLE :: becsum_nc(:,:,:,:)
       !
       INTEGER :: ipol, kh, kkb, is1, is2, js
       !

       IF (okvan) THEN
          IF (noncolin) THEN
             ALLOCATE(becsum_nc(nhm*(nhm+1)/2,nat,npol,npol))
             becsum_nc=(0.d0, 0.d0)
             ALLOCATE( becp_nc( nkb, npol, nbnd ) )
          ELSE
             ALLOCATE( becp( nkb, nbnd ) )
          END IF
       ENDIF
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
             write(6,*) 'trying to reqinw', iunigk
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          IF ( nks > 1 ) THEN
             !
             write(6,*) 'trying to read', iunigk
             CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
                  igk, g2kin)

!             READ( iunigk ) igk
             write(6,*) 'trying to read', iunwfc
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             !
          END IF
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
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
                IF ( upf(np)%tvanp ) THEN
                   !
                   DO na = 1, nat
                      !
                      IF (ityp(na)==np) THEN
                         !
                         ijh = 1
                         !
                         DO ih = 1, nh(np)
                            !
                            ikb = ijkb0 + ih
                            !
                            IF (noncolin) THEN
                               !
                               DO is=1,npol
                                  !
                                  DO js=1,npol
                                     becsum_nc(ijh,na,is,js) =         &
                                         becsum_nc(ijh,na,is,js)+w1 *  &
                                          CONJG(becp_nc(ikb,is,ibnd)) * &
                                                becp_nc(ikb,js,ibnd)
                                  END DO
                                  !
                               END DO
                               !
                            ELSE
                               !
                               becsum(ijh,na,current_spin) = &
                                        becsum(ijh,na,current_spin) + &
                                        w1 * DBLE( CONJG( becp(ikb,ibnd) ) * &
                                                          becp(ikb,ibnd) )
                               !
                            END IF                       
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               IF (noncolin) THEN
                                  !
                                  DO is=1,npol
                                     !
                                     DO js=1,npol
                                        becsum_nc(ijh,na,is,js) =         &
                                           becsum_nc(ijh,na,is,js) + w1 * &
                                           CONJG(becp_nc(ikb,is,ibnd)) *  &
                                                 becp_nc(jkb,js,ibnd)
                                     END DO
                                     !
                                  END DO
                                  !
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

       IF (noncolin.and.okvan) THEN
          DO np = 1, ntyp
             IF ( upf(np)%tvanp ) THEN
                DO na = 1, nat
                   IF (ityp(na)==np) THEN
                      IF (so(np)) THEN
                         CALL transform_becsum_so(becsum_nc,becsum,na)
                      ELSE
                         CALL transform_becsum_nc(becsum_nc,becsum,na)
                      END IF
                   END IF
                END DO
             END IF
          END DO
       END IF
       !
       IF (okvan) THEN
          IF (noncolin) THEN
             DEALLOCATE( becsum_nc )
             DEALLOCATE( becp_nc )
          ELSE
             DEALLOCATE( becp )
          ENDIF
       END IF
       !
       RETURN
       !
     END SUBROUTINE compute_becsum_k
     !
END SUBROUTINE compute_becsum
