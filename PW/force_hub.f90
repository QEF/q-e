!
! Copyright (C) 2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
SUBROUTINE force_hub(forceh)
   !----------------------------------------------------------------------
   !
   ! This routine computes the Hubbard contribution to the force. It gives
   ! in output the product (dE_{hub}/dn_{ij}^{alpha})(dn_{ij}^{alpha}
   ! /du(alpha,ipol)) which is the force acting on the atom at tau_{alpha}
   ! (in the unit ceel) along the direction ipol.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE cell_base,            ONLY : at, bg
   USE ldaU,                 ONLY : hubbard_lmax, hubbard_l, hubbard_u, &
                                    hubbard_alpha, ns, U_projection, &
                                    swfcatom
   USE symme,                ONLY : s, nsym, irt
   USE io_files,             ONLY : prefix, iunocc
   USE wvfct,                ONLY : gamma_only, nbnd, npwx, npw, igk
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE mp_global,            ONLY : me_pool, my_pool_id
   USE basis,                ONLY : natomwfc
   USE becmod,               ONLY : becp
   USE uspp,                 ONLY : nkb, vkb
   USE wavefunctions_module, ONLY : evc
   USE klist,                ONLY : nks, xk
   USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, &
                                    iunat, iunsat, nwordatwfc
   USE atom,                 ONLY : nchi, lchi, oc

   IMPLICIT NONE
   REAL (DP) :: forceh(3,nat)  ! output: the Hubbard forces

   COMPLEX (DP), ALLOCATABLE :: &
             proj(:,:), wfcatom(:,:), spsi(:,:)
   !         proj(natomwfc,nbnd), wfcatom(npwx,natomwfc), spsi(npwx,nbnd)
   REAL (DP), ALLOCATABLE :: dns(:,:,:,:)
   !       dns(ldim,ldim,nspin,nat) ! the derivative of the atomic occupations
   INTEGER, ALLOCATABLE :: offset(:)
   ! offset(nat) : offset of d electrons of atom d in the natomwfc ordering

   COMPLEX (DP) :: c_one, c_zero

   INTEGER :: alpha, na, nt, is, m1, m2, ipol, ldim, l, n, ik
   INTEGER :: counter

   LOGICAL ::  exst

   IF (U_projection .NE. "atomic") CALL errore("force_hub", &
                   " forces for this U_projection_type not implemented",1)
   IF (gamma_only) CALL errore('force_huh',&
                   ' LDA+U, forces AND gamma-only not implemented yet',1)

   ldim= 2 * Hubbard_lmax + 1
   ALLOCATE(dns(ldim,ldim,nspin,nat), offset(nat), &
            proj(natomwfc,nbnd), wfcatom(npwx,natomwfc), spsi(npwx,nbnd), &
            becp(nkb,nbnd) )

   forceh(:,:) = 0.d0

   IF ( me_pool == 0 .AND. my_pool_id == 0 ) THEN

      CALL seqopn (iunocc, 'occup', 'formatted', exst)
      READ(iunocc,*) ns
      CLOSE(unit=iunocc,status='keep')

   END IF

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
      CALL davcio(swfcatom,nwordatwfc,iunsat,ik,-1)
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

!      CALL atomic_wfc( ik, wfcatom )
! read atomic wfc instead
      CALL davcio(wfcatom,nwordatwfc,iunat,ik,-1)
   
      DO ipol = 1,3
         DO alpha = 1,nat                 ! the displaced atom
            CALL dndtau_of_k(dns,ldim,offset,proj,wfcatom,spsi,alpha,ipol,ik)
            DO na = 1,nat                 ! the Hubbard atom
               nt = ityp(na)
               IF (Hubbard_U(nt).NE.0.d0.OR. Hubbard_alpha(nt).NE.0.d0) THEN
                  DO is = 1,nspin
                     DO m2 = 1,ldim
                        forceh(ipol,alpha) = forceh(ipol,alpha) -  &
                              Hubbard_U(nt) * 0.5d0           * dns(m2,m2,is,na)
                        DO m1 = 1,ldim
                           forceh(ipol,alpha) = forceh(ipol,alpha) +    &
                              Hubbard_U(nt) * ns(m2,m1,is,na) * dns(m1,m2,is,na)
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO
   END DO

#ifdef __PARA
   CALL poolreduce(3*nat,forceh)
#endif


   DEALLOCATE(dns, offset, proj, wfcatom, spsi, becp)
   
   IF (nspin.EQ.1) forceh(:,:) = 2.d0 * forceh(:,:)
   !
   ! The symmetry matrices are in the crystal basis so...
   ! Transform to crystal axis...
   !
   DO na=1, nat
      CALL trnvect(forceh(1,na),at,bg,-1)
   END DO
   !
   ! ...symmetrize...
   !
   CALL symvect(nat,forceh,nsym,s,irt)
   !
   ! ... and transform back to cartesian axis
   !
   DO na=1, nat
      CALL trnvect(forceh(1,na),at,bg, 1)
   END DO

   RETURN
END SUBROUTINE force_hub
