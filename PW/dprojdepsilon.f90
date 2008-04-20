!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE dprojdepsilon ( ik,dproj,wfcatom,spsi,ipol,jpol )
   !-----------------------------------------------------------------------
   !
   ! This routine computes the first derivative of the projection
   ! <\fi^{at}_{I,m1}|S|\psi_{k,v,s}> with respect to the strain epsilon(i,j)
   ! (we remember that ns_{I,s,m1,m2} = \sum_{k,v}
   ! f_{kv} <\fi^{at}_{I,m1}|S|\psi_{k,v,s}><\psi_{k,v,s}|S|\fi^{at}_{I,m2}>)
   !
   USE kinds,                ONLY : DP
   USE basis,                ONLY : natomwfc
   USE cell_base,            ONLY : tpiba
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp
   USE gvect,                ONLY : g
   USE klist,                ONLY : nks, xk
   USE ldaU,                 ONLY : swfcatom, Hubbard_l, &
                                    Hubbard_U, Hubbard_alpha
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
   USE uspp,                 ONLY : nkb, vkb, qq
   USE uspp_param,           ONLY : upf, nhm, nh
   USE wavefunctions_module, ONLY : evc
   USE becmod,               ONLY : becp
   USE mp_global,            ONLY : intra_pool_comm
   USE mp,                   ONLY : mp_sum

   IMPLICIT NONE
   !
   ! I/O variables first
   !
   INTEGER :: ik, ipol, jpol
   COMPLEX (DP) :: &
           dproj(natomwfc,nbnd),   &! output: the derivative of the projection
           wfcatom(npwx,natomwfc), &! input: the atomic wfc
           spsi(npwx,nbnd)          ! input: S|evc>
   INTEGER :: i, ig, jkb2, lmax_wfc, na, ibnd, iwf, nt, ib, ih,jh, &
              nworddw, nworddb
   REAL (DP) :: xyz(3,3), q, eps, a1, a2
   PARAMETER (eps=1.0d-8)

   COMPLEX (DP) :: ZDOTC

   COMPLEX (DP), ALLOCATABLE :: &
           dwfc(:,:), aux(:,:), work(:), dbeta(:,:), aux1(:,:), &
           betapsi(:,:), dbetapsi(:,:), wfatbeta(:,:), wfatdbeta(:,:)

   !       dwfc(npwx,natomwfc),   ! the derivative of the atomic d wfc
   !       aux(npwx,natomwfc),    ! auxiliary array
   !       work(npwx),            ! the beta function
   !       dbeta(npwx,nkb),       ! the derivative of the beta function
   !       aux1(npwx,nkb),        ! auxiliary array
   !       betapsi(nhm,nbnd),     ! <beta|evc>
   !       dbetapsi(nhm,nbnd),    ! <dbeta|evc>
   !       wfatbeta(natomwfc,nhm),! <wfc|beta>
   !       wfatdbeta(natomwfc,nhm)! <wfc|dbeta>

   REAL (DP), ALLOCATABLE :: gk(:,:), qm1(:)
   !       gk(3,npwx),
   !       qm1(npwx)
   ! xyz are the three unit vectors in the x,y,z directions
   xyz(:,:) = 0.d0
   DO i=1,3
      xyz(i,i) = 1.d0
   END DO

   dproj(:,:) = (0.d0,0.d0)
   !      WRITE( stdout,*) 'dprojde: ik =',ik,' ipol =',ipol,' jpol =',jpol
   !
   ! At first the derivatives of the atomic wfcs: we compute the term
   ! <d\fi^{at}_{I,m1}/d\epsilon(ipol,jpol)|S|\psi_{k,v,s}>
   !
   ALLOCATE ( qm1(npwx), gk(3,npwx) )
   ALLOCATE ( dwfc(npwx,natomwfc), aux(npwx,natomwfc) )

   nworddw = 2*npwx*natomwfc
   nworddb = 2*npwx*nkb

   lmax_wfc = 0
   DO nt=1, ntyp
      lmax_wfc=MAX(lmax_wfc,MAXVAL(upf(nt)%lchi(1:upf(nt)%nwfc)))
   END DO

   ! here the derivative of the Bessel function
      CALL gen_at_dj (ik,natomwfc,lmax_wfc,dwfc)
   ! and here the derivative of the spherical harmonic
      CALL gen_at_dy (ik,natomwfc,lmax_wfc,xyz(1,ipol),aux)

   DO ig = 1,npw
      gk(1,ig) = (xk(1,ik)+g(1,igk(ig)))*tpiba
      gk(2,ig) = (xk(2,ik)+g(2,igk(ig)))*tpiba
      gk(3,ig) = (xk(3,ik)+g(3,igk(ig)))*tpiba
      q = SQRT(gk(1,ig)**2+gk(2,ig)**2+gk(3,ig)**2)
      IF (q.GT.eps) THEN
         qm1(ig)=1.d0/q
      ELSE
         qm1(ig)=0.d0
      END IF
      a1 = -1.d0*gk(jpol,ig)
      a2 = -1.d0*gk(ipol,ig)*gk(jpol,ig)*qm1(ig)
      DO iwf = 1,natomwfc
         dwfc(ig,iwf) = aux(ig,iwf)*a1 + dwfc(ig,iwf)*a2
         IF (ipol.EQ.jpol) dwfc(ig,iwf) = dwfc(ig,iwf) - wfcatom(ig,iwf)*0.5d0
         DO ibnd = 1,nbnd
            dproj(iwf,ibnd) = dproj(iwf,ibnd)+CONJG(dwfc(ig,iwf))*spsi(ig,ibnd)
         END DO
      END DO
   END DO

   a1 = 0.d0
   a2 = 0.d0

#ifdef __PARA
   CALL mp_sum( dproj, intra_pool_comm )
#endif

   DEALLOCATE ( dwfc, aux )
   !
   ! Now the derivatives of the beta functions: we compute the term
   ! <\fi^{at}_{I,m1}|dS/d\epsilon(ipol,jpol)|\psi_{k,v,s}>
   !
   ALLOCATE (dbeta(npwx,nkb), aux1(npwx,nkb), work(npwx), &
             dbetapsi(nhm,nbnd), betapsi(nhm,nbnd), wfatbeta(natomwfc,nhm), &
             wfatdbeta(natomwfc,nhm) )

   ! here the derivative of the Bessel function
      CALL gen_us_dj (ik,dbeta)
   ! and here the derivative of the spherical harmonic
      CALL gen_us_dy (ik,xyz(1,ipol),aux1)

   jkb2 = 0
   DO nt=1,ntyp
      DO na=1,nat
         IF ( ityp(na) .EQ. nt ) THEN
            DO ih=1,nh(nt)
               jkb2 = jkb2 + 1
               DO ig = 1,npw
                  work(ig) = vkb(ig,jkb2)
                  ! now we compute the true dbeta function
                  dbeta(ig,jkb2) = - aux1(ig,jkb2)*gk(jpol,ig) - &
                        dbeta(ig,jkb2) * gk(ipol,ig) * gk(jpol,ig) * qm1(ig)
                  IF (ipol.EQ.jpol) &
                     dbeta(ig,jkb2) = dbeta(ig,jkb2) - work(ig)*0.5d0
               END DO
               DO ibnd = 1,nbnd
                  betapsi(ih,ibnd)= becp(jkb2,ibnd)
                  dbetapsi(ih,ibnd)= ZDOTC(npw,dbeta(1,jkb2),1,evc(1,ibnd),1)
               END DO
               DO iwf=1,natomwfc
                  wfatbeta(iwf,ih) = ZDOTC(npw,wfcatom(1,iwf),1,work,1)
                  wfatdbeta(iwf,ih)= ZDOTC(npw,wfcatom(1,iwf),1,dbeta(1,jkb2),1)
               END DO
            END DO
#ifdef __PARA
            CALL mp_sum( dbetapsi, intra_pool_comm )
            CALL mp_sum( wfatbeta, intra_pool_comm )
            CALL mp_sum( wfatdbeta, intra_pool_comm )
#endif
            DO ibnd = 1,nbnd
               DO ih=1,nh(nt)
                  DO jh = 1,nh(nt)
                     DO iwf=1,natomwfc
                        dproj(iwf,ibnd) = dproj(iwf,ibnd) +               &
                                          qq(ih,jh,nt) *                  &
                               ( wfatdbeta(iwf,ih)*betapsi(jh,ibnd) +     &
                                 wfatbeta(iwf,ih)*dbetapsi(jh,ibnd) )
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END DO
   END DO


   DEALLOCATE (dbeta, aux1, work, dbetapsi, betapsi, wfatbeta, wfatdbeta )
   DEALLOCATE ( qm1, gk )

   RETURN
END SUBROUTINE dprojdepsilon
