!
! Copyright (C) 2002-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdtau(dproj,wfcatom,spsi,alpha,ipol,offset)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the first derivative of the projection
   ! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the atomic displacement
   ! u(alpha,ipol) (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
   ! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
   !
   USE kinds,                ONLY : DP
   USE atom,                 ONLY : nchi, lchi, oc
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE basis,                ONLY : natomwfc
   USE cell_base,            ONLY : tpiba
   USE gvect,                ONLY : g
   USE klist,                ONLY : nks, xk
   USE ldaU,                 ONLY : Hubbard_l, Hubbard_U, Hubbard_alpha
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
   USE uspp,                 ONLY : nkb, vkb, qq
   USE uspp_param,           ONLY : nhm, nh
   USE wavefunctions_module, ONLY : evc
   USE becmod,               ONLY : becp
   
   IMPLICIT NONE
   INTEGER :: &
              alpha,   &! input: the displaced atom
              ipol,    &! input: the component of displacement
              offset    ! input: the offset of the wfcs of the atom "alpha"
   COMPLEX (DP) :: &
           wfcatom(npwx,natomwfc), &! input: the atomic wfc
           spsi(npwx,nbnd),        &! input: S|evc>
           dproj(natomwfc,nbnd)     ! output: the derivative of the projection

   INTEGER :: ig, jkb2, na, m1, ibnd, iwf, nt, ih, jh, ldim

   REAL (DP) :: gvec

   COMPLEX (DP):: ZDOTC

   COMPLEX (DP), ALLOCATABLE :: dwfc(:,:), work(:), dbeta(:), &
                                     betapsi(:,:), dbetapsi(:,:), &
                                     wfatbeta(:,:), wfatdbeta(:,:)
   !      dwfc(npwx,ldim),       ! the derivative of the atomic d wfc
   !      work(npwx),            ! the beta function
   !      dbeta(npwx),           ! the derivative of the beta function
   !      betapsi(nhm,nbnd),     ! <beta|evc>
   !      dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !      wfatbeta(natomwfc,nhm),! <wfc|beta>
   !      wfatdbeta(natomwfc,nhm)! <wfc|dbeta>

   nt = ityp(alpha)

   ldim = 2 * Hubbard_l(nt) + 1

   ALLOCATE ( dwfc(npwx,ldim), work(npwx), dbeta(npwx), betapsi(nhm,nbnd), &
         dbetapsi(nhm,nbnd), wfatbeta(natomwfc,nhm), wfatdbeta(natomwfc,nhm) )

   dproj(:,:) = (0.d0, 0.d0)
   !
   ! At first the derivatives of the atomic wfc and the beta are computed
   !
   IF (Hubbard_U(nt).NE.0.d0.OR.Hubbard_alpha(nt).NE.0.d0) THEN
      DO ig = 1,npw
         gvec = g(ipol,igk(ig)) * tpiba

         ! in the expression of dwfc we don't need (k+G) but just G; k always
         ! multiplies the underived quantity and gives an opposite contribution
         ! in c.c. term because the sign of the imaginary unit.
   
         DO m1 = 1, ldim
            dwfc(ig,m1) = CMPLX(0.d0,-1.d0) * gvec * wfcatom(ig,offset+m1)
         END DO
      END DO

      CALL ZGEMM('C','N',ldim, nbnd, npw, (1.d0,0.d0), &
                  dwfc, npwx, spsi, npwx, (0.d0,0.d0), &
                  dproj(offset+1,1), natomwfc)
   END IF

#ifdef __PARA
   CALL reduce(2*natomwfc*nbnd,dproj)
#endif

   jkb2 = 0
   DO nt=1,ntyp
      DO na=1,nat
         IF ( ityp(na) .EQ. nt ) THEN
            DO ih=1,nh(nt)
               jkb2 = jkb2 + 1
               IF (na.EQ.alpha) THEN
                  DO ig = 1, npw
                     gvec = g(ipol,igk(ig)) * tpiba
                     dbeta(ig) = CMPLX(0.d0,-1.d0) * vkb(ig,jkb2) * gvec
                     work(ig) = vkb(ig,jkb2)
                  END DO
                  DO ibnd=1,nbnd
                     dbetapsi(ih,ibnd)= ZDOTC(npw,dbeta,1,evc(1,ibnd),1)
                     betapsi(ih,ibnd) = becp(jkb2,ibnd)
                  END DO
                  DO iwf=1,natomwfc
                     wfatbeta(iwf,ih) = ZDOTC(npw,wfcatom(1,iwf),1,work,1)
                     wfatdbeta(iwf,ih)= ZDOTC(npw,wfcatom(1,iwf),1,dbeta,1)
                  END DO
               END IF
            END DO
#ifdef __PARA
            CALL reduce(2*nhm*nbnd,dbetapsi)
            CALL reduce(2*natomwfc*nhm,wfatbeta)
            CALL reduce(2*natomwfc*nhm,wfatdbeta)
#endif
            IF (na.EQ.alpha) THEN
               DO ibnd=1,nbnd
                  DO ih=1,nh(nt)
                     DO jh=1,nh(nt)
                        DO iwf=1,natomwfc
                           dproj(iwf,ibnd) = &
                               dproj(iwf,ibnd) + qq(ih,jh,nt) *         &
                               ( wfatdbeta(iwf,ih)*betapsi(jh,ibnd) +   &
                                  wfatbeta(iwf,ih)*dbetapsi(jh,ibnd) )
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         END IF
      END DO
   END DO

   DEALLOCATE ( dwfc, work, dbeta, betapsi, dbetapsi, wfatbeta, wfatdbeta )

   RETURN
END SUBROUTINE dprojdtau
