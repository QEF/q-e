
!
! Copyright (C) 2006 QUANTUM-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
subroutine compute_ppsi (ppsi, ppsi_us, ik, ipol, nbnd_occ, current_spin)
  !----------------------------------------------------------------------
  !
  ! On output: ppsi contains P_c^+ p | psi_ik > for the ipol cartesian
  !            coordinate
  !            ppsi_us contains the additional term required for US PP.
  !            See J. Chem. Phys. 120, 9935 (2004) Eq. 10.
  !
  ! (vkb,evc,igk must be set)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba
  USE io_global,            ONLY : stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : et, nbnd, npwx, npw, igk, g2kin
  USE uspp,                 ONLY : nkb, vkb, deeq, qq, qq_so, deeq_nc, okvan
  USE spin_orb,             ONLY : lspinorb
  USE lsda_mod,             ONLY : nspin
  USE gvect,                ONLY : g
  USE klist,                ONLY : xk, nks
  USE noncollin_module,     ONLY : noncolin, npol
  USE becmod,               ONLY : becp, becp_nc
  USE uspp_param,           ONLY : nh, nhm
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ipol, ik, nbnd_occ, current_spin
  !
  COMPLEX(DP) :: ppsi(npwx,npol,nbnd_occ), ppsi_us(npwx,npol,nbnd_occ)
  ! Local variables
  !
  INTEGER :: ig, na, ibnd, jbnd, ikb, jkb, nt, ih, jh, ip, ijkb0
  ! counters

  REAL(DP), ALLOCATABLE  :: gk (:,:)
  ! the derivative of |k+G|
  REAL(DP)  :: vers(3)

  COMPLEX(DP), ALLOCATABLE :: ps2(:,:,:), dvkb (:,:), dvkb1 (:,:),   &
       work (:,:), becp2(:,:), becp2_nc(:,:,:), psc(:,:,:,:), ps(:), &
       ps_nc(:,:), dpqq_so(:,:,:,:,:)

  REAL(DP), ALLOCATABLE :: dpqq(:,:,:,:)

  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  ALLOCATE (work ( npwx, MAX(nkb,1)))
  ALLOCATE (gk ( 3, npwx))    
  IF (nkb > 0) THEN
     IF (noncolin) THEN
        ALLOCATE (becp2_nc (nkb, npol, nbnd))
     ELSE
        ALLOCATE (becp2 (nkb, nbnd))
     END IF

     ALLOCATE (dvkb (npwx, nkb))
     ALLOCATE (dvkb1(npwx, nkb))
     dvkb (:,:) = (0.d0, 0.d0)
     dvkb1(:,:) = (0.d0, 0.d0)
  END IF
  DO ig = 1, npw
     gk (1, ig) = (xk (1, ik) + g (1, igk (ig) ) ) * tpiba
     gk (2, ig) = (xk (2, ik) + g (2, igk (ig) ) ) * tpiba
     gk (3, ig) = (xk (3, ik) + g (3, igk (ig) ) ) * tpiba
     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2
  END DO
  !
  ! this is the kinetic contribution to p :  (k+G)_ipol * psi
  !
  DO ip=1,npol
     DO ibnd = 1, nbnd_occ 
        DO ig = 1, npw
           ppsi(ig,ip,ibnd)=gk(ipol,ig)*evc(ig+npwx*(ip-1),ibnd)
        END DO
     END DO
  END DO

  !
  ! and this is the contribution from nonlocal pseudopotentials
  !
  CALL gen_us_dj (ik, dvkb)
  vers=0.d0
  vers(ipol)=1.d0
  CALL gen_us_dy (ik, vers, dvkb1)
  DO ig = 1, npw
     IF (g2kin (ig) < 1.0d-10) then
        gk (1, ig) = 0.d0
        gk (2, ig) = 0.d0
        gk (3, ig) = 0.d0
     ELSE
        gk (1, ig) = gk (1, ig) / SQRT (g2kin (ig) )
        gk (2, ig) = gk (2, ig) / SQRT (g2kin (ig) )
        gk (3, ig) = gk (3, ig) / SQRT (g2kin (ig) )
     END IF
  END DO

  jkb = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (nt == ityp (na)) THEN
           DO ikb = 1, nh (nt)
              jkb = jkb + 1
              DO ig = 1, npw
                 work (ig,jkb)=dvkb1(ig,jkb)+dvkb(ig,jkb)*gk(ipol,ig) 
              END DO
           END DO
        END IF
     END DO
  END DO
  DEALLOCATE (gk)

  IF (nkb>0) THEN
     IF (noncolin) THEN
        CALL ccalbec_nc(nkb, npwx, npw, npol, nbnd, becp2_nc, work, evc)
     ELSE
        CALL ccalbec(nkb, npwx, npw, nbnd, becp2, work, evc)
     ENDIF
  ENDIF

  ijkb0 = 0
  IF (noncolin) THEN
     ALLOCATE (psc( nkb, 2, nbnd_occ,  2))
     psc=(0.d0,0.d0)
  ELSE
     ALLOCATE (ps2( nkb, nbnd_occ, 2))
     ps2=(0.d0,0.d0)
  END IF
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (nt == ityp (na)) THEN
           DO ih = 1, nh (nt)
              ikb = ijkb0 + ih
              DO jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 DO ibnd = 1, nbnd_occ 
                    IF (noncolin) THEN
                       IF (lspinorb) THEN
                          psc(ikb,1,ibnd,1)=psc(ikb,1,ibnd,1)+(0.d0,-1.d0)* &
                             (becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1)  &
                                 -et(ibnd,ik)*qq_so(ih,jh,1,nt) )+       &
                              becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,2)- &
                                       et(ibnd,ik)* qq_so(ih,jh,2,nt) ) )
                          psc(ikb,2,ibnd,1)=psc(ikb,2,ibnd,1)+(0.d0,-1.d0)*  &
                             (becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,3)  &
                                 -et(ibnd,ik)*qq_so(ih,jh,3,nt) )+       &
                              becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4)- &
                                       et(ibnd,ik)* qq_so(ih,jh,4,nt) ) )
                          psc(ikb,1,ibnd,2)=psc(ikb,1,ibnd,2)+(0.d0,-1.d0)* &
                             (becp_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1)  &
                                 -et(ibnd,ik)*qq_so(ih,jh,1,nt) )+      &
                             becp_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,2)-  &
                                       et(ibnd,ik)* qq_so(ih,jh,2,nt) ) )
                          psc(ikb,2,ibnd,2)=psc(ikb,2,ibnd,2)+(0.d0,-1.d0)*  &
                             (becp_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,3)  &
                                 -et(ibnd,ik)*qq_so(ih,jh,3,nt) )+      &
                             becp_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4)-  &
                                       et(ibnd,ik)* qq_so(ih,jh,4,nt) ) )
                       ELSE
                          psc(ikb,1,ibnd,1)=psc(ikb,1,ibnd,1)+ (0.d0,-1.d0)* &
                              ( becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1) &
                                             -et(ibnd,ik)*qq(ih,jh,nt)) + &
                                becp2_nc(jkb,2,ibnd)*deeq_nc(ih,jh,na,2) )
                          psc(ikb,2,ibnd,1)=psc(ikb,2,ibnd,1)+ (0.d0,-1.d0)* &
                              ( becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4) &
                                             -et(ibnd,ik)*qq(ih,jh,nt))+  &
                                becp2_nc(jkb,1,ibnd)*deeq_nc(ih,jh,na,3) )
                          psc(ikb,1,ibnd,2)=psc(ikb,1,ibnd,2)+ (0.d0,-1.d0)* &
                              ( becp_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1) &
                                             -et(ibnd,ik)*qq(ih,jh,nt))+ &
                                becp_nc(jkb,2,ibnd)*deeq_nc(ih,jh,na,2) )
                          psc(ikb,2,ibnd,2)=psc(ikb,2,ibnd,2)+ (0.d0,-1.d0)* &
                              ( becp_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4) &
                                             -et(ibnd,ik)*qq(ih,jh,nt))+ &
                                becp_nc(jkb,1,ibnd)*deeq_nc(ih,jh,na,3) )
                       END IF
                    ELSE
                       ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                         (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin) &
                         -et(ibnd,ik)*qq(ih,jh,nt))
                       ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp(jkb,ibnd) * &
                         (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin)&
                         -et(ibnd,ik)*qq(ih,jh,nt))
                    ENDIF
                 END DO
              END DO
           END DO
           ijkb0=ijkb0+nh(nt)
        END IF
     END DO
  END DO
  IF (ikb /= nkb .OR. jkb /= nkb) CALL errore ('compute_ppsi', &
                                               'unexpected error',1)

  IF (nkb>0) THEN
     IF (noncolin) THEN
        CALL ZGEMM( 'N', 'N', npwx, nbnd_occ*npol, nkb, &
             (0.d0,0.5d0), vkb, npwx, psc(1,1,1,1), nkb, (1.d0,0.d0), &
              ppsi, npwx )
        CALL ZGEMM( 'N', 'N', npwx, nbnd_occ*npol, nkb, &
             (0.d0,0.5d0), work, npwx, psc(1,1,1,2), nkb, (1.d0,0.d0), &
             ppsi, npwx )
     ELSE
        CALL ZGEMM( 'N', 'N', npw, nbnd_occ, nkb, &
             (0.d0,0.5d0), vkb(1,1), npwx, ps2(1,1,1), nkb, (1.d0,0.0d0), &
             ppsi, npwx )
        CALL ZGEMM( 'N', 'N', npw, nbnd_occ, nkb, &
             (0.d0,0.5d0), work(1,1), npwx, ps2(1,1,2), nkb, (1.d0,0.0d0), &
             ppsi, npwx )
     END IF
  END IF
  IF (noncolin) THEN
     DEALLOCATE (psc)
  ELSE
     DEALLOCATE (ps2)
  END IF
!
!   ppsi contains p - i/2 [x, V_{nl}-eS] psi_v for the ipol polarization
!
!   In the US case there is another term in the matrix element.
!   This term has to be multiplied by the difference of the eigenvalues,
!   so it is calculated separately here and multiplied in the calling
!   routine.

  IF (okvan) THEN
     ppsi_us=(0.d0,0.d0)
     ALLOCATE (dpqq( nhm, nhm, 3, ntyp))
     CALL compute_qdipol(dpqq,ipol)
     IF (noncolin) THEN
        ALLOCATE (ps_nc(nbnd_occ,npol))    
        IF (lspinorb) THEN
           ALLOCATE (dpqq_so( nhm, nhm, nspin, 3, ntyp))
           CALL compute_qdipol_so(dpqq, dpqq_so,ipol)
        END IF
     ELSE
        ALLOCATE (ps(nbnd_occ))    
     ENDIF
     ijkb0 = 0
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (ityp(na)==nt) THEN
              DO ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 IF (noncolin) THEN
                    ps_nc = (0.d0,0.d0)
                 ELSE
                    ps = (0.d0,0.d0)
                 END IF
                 DO jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    DO ibnd=1, nbnd_occ
                       IF (noncolin) THEN
                          DO ip=1,npol
                             IF (lspinorb) THEN
                                ps_nc(ibnd,ip)=ps_nc(ibnd,ip) +          &
                                    (0.d0,1.d0)*(becp2_nc(jkb,1,ibnd)*   &
                                    qq_so(ih,jh,1+(ip-1)*2,nt) +         &
                                    becp2_nc(jkb,2,ibnd) *               &
                                    qq_so(ih,jh,2+(ip-1)*2,nt) )         &
                                  + becp_nc(jkb,1,ibnd)*                 &
                                    dpqq_so(ih,jh,1+(ip-1)*2,ipol,nt)    &
                                  + becp_nc(jkb,2,ibnd)*                 &
                                    dpqq_so(ih,jh,2+(ip-1)*2,ipol,nt)
                             ELSE
                                ps_nc(ibnd,ip)=ps_nc(ibnd,ip)+           &
                                    becp2_nc(jkb,ip,ibnd)*(0.d0,1.d0)*   &
                                    qq(ih,jh,nt)+becp_nc(jkb,ip,ibnd)    &
                                                   *dpqq(ih,jh,ipol,nt)
                             END IF
                          END DO
                       ELSE
                          ps(ibnd) = ps(ibnd) + becp2(jkb,ibnd) *  &
                                (0.d0,1.d0) * qq(ih,jh,nt)   +  &
                                becp(jkb,ibnd) * dpqq(ih,jh,ipol,nt)
                       END IF
                    END DO
                 END DO
                 DO ibnd = 1, nbnd_occ 
                    IF (noncolin) THEN
                       DO ip=1,npol
                          CALL ZAXPY(npw,ps_nc(ibnd,ip),vkb(1,ikb),1,&
                                     ppsi_us(1,ip,ibnd),1)
                       END DO
                    ELSE
                       CALL ZAXPY(npw,ps(ibnd),vkb(1,ikb),1,ppsi_us(1,1,ibnd),1)
                    ENDIF
                 END DO
              END DO
              ijkb0=ijkb0+nh(nt)
           END IF
        END DO
     END DO
     IF (jkb.NE.nkb) CALL errore ('compute_ppsi', 'unexpected error', 1)
     IF (noncolin) THEN
        DEALLOCATE(ps_nc)
     ELSE
        DEALLOCATE(ps)
     END IF
  ENDIF


  IF (nkb > 0) THEN
     DEALLOCATE (dvkb1, dvkb)
     IF (noncolin) THEN
        DEALLOCATE(becp2_nc)
     ELSE
        DEALLOCATE(becp2)
     END IF
  END IF
  DEALLOCATE (work)

  RETURN
END SUBROUTINE compute_ppsi
