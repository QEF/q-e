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
SUBROUTINE dndtau(dns,ldim,alpha,ipol)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
   USE kinds,                ONLY : DP
   USE atom,                 ONLY : nchi, lchi, oc
   USE ions_base,            ONLY : nat, ityp
   USE basis,                ONLY : natomwfc
   USE klist,                ONLY : nks, xk
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE ldaU,                 ONLY : swfcatom, Hubbard_l, &
                                    Hubbard_U, Hubbard_alpha
   USE wavefunctions_module, ONLY : evc
   USE uspp,                 ONLY : nkb, vkb
   USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
   USE becmod,               ONLY : becp
   USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, nwordatwfc, iunat
   
   IMPLICIT NONE

   INTEGER ::               &
              ik,           & ! counter on k points
              ibnd,         & !    "    "  bands
              is,           & !    "    "  spins
              i, na, nt, n, alpha, ipol, counter, m1,m2, l
   INTEGER :: ldim
   REAL (kind=DP) :: &
             dns(ldim,ldim,nspin,nat)
   COMPLEX (kind=DP) :: ZDOTC, c_one, c_zero
   INTEGER, ALLOCATABLE :: offset(:)
   ! offset(nat): offset of d electrons of atom d in the natomwfc ordering
   COMPLEX (kind=DP), ALLOCATABLE :: &
                      proj(:,:), wfcatom(:,:), spsi(:,:),dproj(:,:)
   !                  proj(natomwfc,nbnd), wfcatom(npwx,natomwfc),
   !                  spsi(npwx,nbnd), dproj(natomwfc,nbnd)

   CALL start_clock('dndtau')

   ALLOCATE ( offset(nat) )
   ALLOCATE ( proj(natomwfc,nbnd), dproj(natomwfc,nbnd), &
              spsi(npwx,nbnd), wfcatom(npwx,natomwfc), becp(nkb,nbnd) )

   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   counter = 0
   DO na=1,nat
      offset(na) = 0
      nt=ityp(na)
      DO n=1,nchi(nt)
         IF (oc(n,nt) >= 0.d0) THEN
            l=lchi(n,nt)
            IF (l == Hubbard_l(nt)) offset(na) = counter
            counter = counter + 2 * l + 1
         END IF
      END DO
   END DO

   IF(counter.NE.natomwfc)CALL errore('new_ns','nstart<>counter',1)

   dns(:,:,:,:) = 0.d0
   !
   !    we start a loop on k points
   !
   IF (nks.GT.1) REWIND (iunigk)

   DO ik = 1, nks

      IF (lsda) current_spin = isk(ik)
      !
      ! now we need the first derivative of proj with respect to tau(alpha,ipol)
      !

      IF (nks.GT.1) READ (iunigk) npw, igk

      CALL davcio(evc,nwordwfc,iunwfc,ik,-1)
      CALL davcio(swfcatom,nwordatwfc,iunat,ik,-1)
      c_one= (1.d0, 0.d0)
      c_zero = (0.d0, 0.d0)
      CALL ZGEMM ('C', 'N', natomwfc, nbnd, npw, c_one, &
                            swfcatom, npwx, evc, npwx, c_zero, proj, natomwfc)
#ifdef __PARA
      CALL reduce(2*natomwfc*nbnd,proj)
#endif

      CALL init_us_2 (npw,igk,xk(1,ik),vkb)

      CALL ccalbec(nkb, npwx, npw, nbnd, becp, vkb, evc)

      CALL s_psi  (npwx, npw, nbnd, evc, spsi )

      CALL atomic_wfc( ik, wfcatom )

      dproj(:,:) = (0.d0,0.d0)
      CALL dprojdtau(dproj,wfcatom,spsi,alpha,ipol,offset(alpha))
      !
      ! compute the derivative of occupation numbers (the quantities dn(m1,m2))
      ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
      !
      DO na = 1,nat
         nt = ityp(na)
         IF (Hubbard_U(nt).NE.0.d0.OR.Hubbard_alpha(nt).NE.0.d0) THEN
            DO m1 = 1,ldim
               DO m2 = m1,ldim
                  DO ibnd = 1,nbnd
                     dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                             wg(ibnd,ik) *            &
                                DREAL(  proj(offset(na)+m1,ibnd)  *   &
                                CONJG(dproj(offset(na)+m2,ibnd))  +   &
                                       dproj(offset(na)+m1,ibnd)  *   &
                                CONJG(proj(offset(na)+m2,ibnd)) )
                  END DO
               END DO
            END DO
         END IF
      END DO
   END DO                 ! on k-points

#ifdef __PARA
   CALL poolreduce(ldim*ldim*nspin*nat,dns)
#endif
   !
   ! impose hermiticity of dn_{m1,m2}
   !
   DO na = 1,nat
      DO is = 1,nspin
         DO m1 = 1,ldim
            DO m2 = m1+1,ldim
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            END DO
         END DO
      END DO
   END DO

   DEALLOCATE ( proj, dproj, spsi, wfcatom, becp )
   DEALLOCATE ( offset )

   CALL stop_clock('dndtau')
   RETURN
END SUBROUTINE dndtau

