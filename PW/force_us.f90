!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE force_us( forcenl )
  !----------------------------------------------------------------------------
  !
  ! ... nonlocal potential contribution to forces
  ! ... wrapper
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk
  USE gvect,                ONLY : g
  USE uspp,                 ONLY : nkb, vkb, qq, deeq, qq_so, deeq_nc
  USE uspp_param,           ONLY : upf, nh, newpseudo
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE symme,                ONLY : irt, s, nsym
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin
  USE spin_orb,             ONLY : lspinorb
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE buffers,              ONLY : get_buffer
  USE becmod,               ONLY : allocate_bec, deallocate_bec
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! ... the dummy variable
  !
  REAL(DP) :: forcenl(3,nat)
  ! output: the nonlocal contribution
  !
  CALL allocate_bec ( nkb, nbnd )   
  !
  IF ( gamma_only ) THEN
     !
     CALL force_us_gamma( forcenl )
     !
  ELSE
     !
     CALL force_us_k( forcenl )
     !
  END IF  
  !
  CALL deallocate_bec ( )   
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_gamma( forcenl )
       !-----------------------------------------------------------------------
       !
       ! ... calculation at gamma
       !
       USE becmod, ONLY : rbecp, calbec
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       REAL(DP), ALLOCATABLE    :: rdbecp (:,:,:)
       ! auxiliary variable, contains <dbeta|psi>
       COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)
       ! auxiliary variable contains g*|beta>
       REAL(DP) :: ps
       INTEGER       :: ik, ipol, ibnd, ig, ih, jh, na, nt, ikb, jkb, ijkb0
       ! counters
       !
       !
       forcenl(:,:) = 0.D0
       !
       ALLOCATE( rdbecp( nkb, nbnd, 3 ) )    
       ALLOCATE( vkb1(  npwx, nkb ) ) 
       !   
       IF ( nks > 1 ) REWIND iunigk
       !
       ! ... the forces are a sum over the K points and the bands
       !
       DO ik = 1, nks
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk (ik)
          IF ( nks > 1 ) THEN
             READ( iunigk ) igk
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             IF ( nkb > 0 ) &
                CALL init_us_2( npw, igk, xk(1,ik), vkb )
          END IF
          !
          CALL calbec ( npw, vkb, evc, rbecp )
          !
          DO ipol = 1, 3
             DO jkb = 1, nkb
                DO ig = 1, npw
                   vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * g(ipol,igk(ig))
                END DO
             END DO
             !
             CALL calbec ( npw, vkb1, evc, rdbecp(:,:,ipol) )
             !
          END DO
          !
          ijkb0 = 0
          DO nt = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      DO ibnd = 1, nbnd
                         ps = deeq(ih,ih,na,current_spin) - &
                              et(ibnd,ik) * qq(ih,ih,nt)
                         DO ipol = 1, 3
                            forcenl(ipol,na) = forcenl(ipol,na) - &
                                       ps * wg(ibnd,ik) * 2.D0 * tpiba * &
                                       rdbecp(ikb,ibnd,ipol) *rbecp(ikb,ibnd)
                         END DO
                      END DO
                      !
                      IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd = 1, nbnd
                               ps = deeq(ih,jh,na,current_spin) - &
                                    et(ibnd,ik) * qq(ih,jh,nt)
                               DO ipol = 1, 3
                                  forcenl(ipol,na) = forcenl(ipol,na) - &
                                     ps * wg(ibnd,ik) * 2.d0 * tpiba * &
                                     (rdbecp(ikb,ibnd,ipol) *rbecp(jkb,ibnd) + &
                                      rdbecp(jkb,ibnd,ipol) *rbecp(ikb,ibnd) )
                               END DO
                            END DO
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          END DO
       END DO
       !
       ! ... The total D matrix depends on the ionic position via the
       ! ... augmentation part \int V_eff Q dr, the term deriving from the 
       ! ... derivative of Q is added in the routine addusforce
       !
       CALL addusforce( forcenl )
       !
#ifdef __PARA
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( forcenl, inter_pool_comm )
#endif
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the forces. The symmetry matrices are 
       ! ... in the crystal basis so...
       ! ... Transform to crystal axis...
       !
       DO na = 1, nat
          CALL trnvect( forcenl(1,na), at, bg, -1 )
       END DO
       !
       ! ... symmetrize...
       !
       CALL symvect( nat, forcenl, nsym, s, irt )
       !
       ! ... and transform back to cartesian axis
       !
       DO na = 1, nat
          CALL trnvect( forcenl(1,na), at, bg, 1 )
       END DO
       !
       DEALLOCATE( vkb1 )
       DEALLOCATE(rdbecp ) 
       !
       RETURN
       !
     END SUBROUTINE force_us_gamma
     !     
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_k( forcenl )
       !-----------------------------------------------------------------------
       !  
       USE becmod, ONLY : becp, becp_nc, calbec
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       COMPLEX(DP), ALLOCATABLE :: dbecp(:,:,:), dbecp_nc(:,:,:,:)
       ! auxiliary variable contains <beta|psi> and <dbeta|psi>
       COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)
       ! auxiliary variable contains g*|beta>
       COMPLEX(DP) :: psc(2,2), fac
       REAL(DP) :: ps
       INTEGER       :: ik, ipol, ibnd, ig, ih, jh, na, nt, ikb, jkb, ijkb0, &
                        is, js
       ! counters
       !
       !
       forcenl(:,:) = 0.D0
       !
       IF (noncolin) then
          ALLOCATE( dbecp_nc(nkb,npol,nbnd,3) )    
       ELSE
          ALLOCATE( dbecp( nkb, nbnd, 3 ) )    
       ENDIF
       ALLOCATE( vkb1( npwx, nkb ) )   
       ! 
       IF ( nks > 1 ) REWIND iunigk
       !
       ! ... the forces are a sum over the K points and the bands
       !
       DO ik = 1, nks
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          IF ( nks > 1 ) THEN
             READ( iunigk ) igk
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             IF ( nkb > 0 ) &
                CALL init_us_2( npw, igk, xk(1,ik), vkb )
          END IF
          !
          IF (noncolin) THEN
             CALL calbec ( npw, vkb, evc, becp_nc)
          ELSE
             CALL calbec ( npw, vkb, evc, becp)
          ENDIF
          !
          DO ipol = 1, 3
             DO jkb = 1, nkb
                DO ig = 1, npw
                   vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*g(ipol,igk(ig))
                END DO
             END DO
             !
             IF (noncolin) THEN
                IF ( nkb > 0 ) &
                   CALL ZGEMM( 'C', 'N', nkb, nbnd*npol, npw, ( 1.D0, 0.D0 ),&
                            vkb1, npwx, evc, npwx, ( 0.D0, 0.D0 ),    &
                            dbecp_nc(1,1,1,ipol), nkb )
             ELSE
                IF ( nkb > 0 ) &
                   CALL ZGEMM( 'C', 'N', nkb, nbnd, npw, ( 1.D0, 0.D0 ),   &
                            vkb1, npwx, evc, npwx, ( 0.D0, 0.D0 ),      &
                            dbecp(1,1,ipol), nkb )
             END IF
          END DO
          !
          ijkb0 = 0
          DO nt = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      DO ibnd = 1, nbnd
                         IF (noncolin) THEN
                            IF (lspinorb) THEN
                               psc(1,1)=deeq_nc(ih,ih,na,1)-et(ibnd,ik)* &
                                                qq_so(ih,ih,1,nt)
                               psc(1,2)=deeq_nc(ih,ih,na,2)-et(ibnd,ik)* &
                                                qq_so(ih,ih,2,nt)
                               psc(2,1)=deeq_nc(ih,ih,na,3)-et(ibnd,ik)* &
                                                qq_so(ih,ih,3,nt)
                               psc(2,2)=deeq_nc(ih,ih,na,4)-et(ibnd,ik)* &
                                                qq_so(ih,ih,4,nt)
                            ELSE
                               psc(1,1)=deeq_nc(ih,ih,na,1)- &
                                                et(ibnd,ik)*qq(ih,ih,nt)
                               psc(1,2)=deeq_nc(ih,ih,na,2)
                               psc(2,1)=deeq_nc(ih,ih,na,3)
                               psc(2,2)=deeq_nc(ih,ih,na,4)- &
                                                et(ibnd,ik)*qq(ih,ih,nt)
                            END IF 
                            fac=wg(ibnd,ik)
                            DO ipol=1,3
                               DO is=1,npol
                                  DO js=1,npol
                                     forcenl(ipol,na) = forcenl(ipol,na)- &
                                         psc(is,js)*fac*tpiba*( &
                                         CONJG(dbecp_nc(ikb,is,ibnd,ipol))* &
                                         becp_nc(ikb,js,ibnd)+ &
                                         CONJG(becp_nc(ikb,is,ibnd))* &
                                         dbecp_nc(ikb,js,ibnd,ipol) )
                                  END DO
                               END DO
                            END DO
                         ELSE
                            ps = deeq(ih,ih,na,current_spin) - &
                                 et(ibnd,ik) * qq(ih,ih,nt)
                            DO ipol=1,3
                               forcenl(ipol,na) = forcenl(ipol,na) - &
                                        ps * wg(ibnd,ik) * 2.D0 * tpiba * &
                                      DBLE( CONJG( dbecp(ikb,ibnd,ipol) ) * &
                                            becp(ikb,ibnd) )
                            END DO
                         END IF
                      END DO
                      !
                      IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd = 1, nbnd
                               IF (noncolin) THEN
                                  IF (lspinorb) THEN
                                    psc(1,1)=deeq_nc(ih,jh,na,1)-et(ibnd,ik)* &
                                             qq_so(ih,jh,1,nt)
                                    psc(1,2)=deeq_nc(ih,jh,na,2)-et(ibnd,ik)* &
                                              qq_so(ih,jh,2,nt)
                                    psc(2,1)=deeq_nc(ih,jh,na,3)-et(ibnd,ik)* &
                                             qq_so(ih,jh,3,nt)
                                    psc(2,2)=deeq_nc(ih,jh,na,4)-et(ibnd,ik)* &
                                             qq_so(ih,jh,4,nt)
                                  ELSE
                                    psc(1,1)=deeq_nc(ih,jh,na,1) &
                                            -et(ibnd,ik)*qq(ih,jh,nt)
                                    psc(1,2)=deeq_nc(ih,jh,na,2)
                                    psc(2,1)=deeq_nc(ih,jh,na,3)
                                    psc(2,2)=deeq_nc(ih,jh,na,4) &
                                            -et(ibnd,ik)*qq(ih,jh,nt)
                                  END IF
                                  fac=wg(ibnd,ik)
                                  DO ipol=1,3
                                     DO is=1,npol
                                        DO js=1,npol
                                           forcenl(ipol,na) &
                                              =forcenl(ipol,na)- &
                                               psc(is,js)*fac*tpiba*( &
                                          CONJG(dbecp_nc(ikb,is,ibnd,ipol))* &
                                                 becp_nc(jkb,js,ibnd)+ &
                                          CONJG(becp_nc(ikb,is,ibnd))* &
                                                dbecp_nc(jkb,js,ibnd,ipol) + &
                                          CONJG(dbecp_nc(jkb,is,ibnd,ipol))* &
                                                becp_nc(ikb,js,ibnd)+ &
                                          CONJG(becp_nc(jkb,is,ibnd))* &
                                                dbecp_nc(ikb,js,ibnd,ipol) )
                                        END DO
                                     END DO
                                  END DO
                               ELSE
                                  ps = deeq(ih,jh,na,current_spin) - &
                                       et(ibnd,ik) * qq (ih,jh,nt)
                                  DO ipol = 1, 3
                                     forcenl(ipol,na) = forcenl (ipol,na) - &
                                          ps * wg(ibnd,ik) * 2.D0 * tpiba * &
                                       DBLE( CONJG( dbecp(ikb,ibnd,ipol) ) * &
                                             becp(jkb,ibnd) +       &
                                             dbecp(jkb,ibnd,ipol) * &
                                             CONJG( becp(ikb,ibnd) ) )
                                  END DO
                               END IF
                            END DO
                          END DO
                      END IF ! tvanp
                   END DO ! ih = 1, nh(nt)
                   ijkb0 = ijkb0 + nh(nt)
                END IF ! ityp(na) == nt
             END DO ! nat
          END DO ! ntyp
       END DO ! nks
       !
#ifdef __PARA
       CALL mp_sum(  forcenl , intra_pool_comm )
#endif
       !
       DEALLOCATE( vkb1 )
       IF (noncolin) THEN
          DEALLOCATE( dbecp_nc )
       ELSE
          DEALLOCATE( dbecp )
       ENDIF
       !
       ! ... The total D matrix depends on the ionic position via the
       ! ... augmentation part \int V_eff Q dr, the term deriving from the 
       ! ... derivative of Q is added in the routine addusforce
       !
       CALL addusforce( forcenl )
       !
#ifdef __PARA
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( forcenl, inter_pool_comm )
#endif
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the forces. The symmetry matrices are 
       ! ... in the crystal basis so...
       ! ... Transform to crystal axis...
       !
       DO na = 1, nat
          CALL trnvect( forcenl(1,na), at, bg, -1 )
       END DO
       !
       ! ... symmetrize...
       !
       CALL symvect( nat, forcenl, nsym, s, irt )
       !
       ! ... and transform back to cartesian axis
       !
       DO na = 1, nat
          CALL trnvect( forcenl(1,na), at, bg, 1 )
       END DO
       !
       RETURN
       !
     END SUBROUTINE force_us_k
     !     
END SUBROUTINE force_us
