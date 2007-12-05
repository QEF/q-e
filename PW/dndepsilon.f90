!
! Copyright (C) 2002-2007  Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE dndepsilon ( dns,ldim,ipol,jpol )
   !-----------------------------------------------------------------------
   ! This routine computes the derivative of the ns atomic occupations with
   ! respect to the strain epsilon(ipol,jpol) used to obtain the hubbard
   ! contribution to the internal stres tensor.
   !
   USE kinds,                ONLY : DP
   USE wavefunctions_module, ONLY : evc
   USE ions_base,            ONLY : nat, ityp
   USE basis,                ONLY : natomwfc
   USE klist,                ONLY : nks, xk, ngk
   USE ldaU,                 ONLY : swfcatom, Hubbard_l, &
                                    Hubbard_U, Hubbard_alpha
   USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
   USE uspp,                 ONLY : nkb, vkb
   USE uspp_param,           ONLY : upf
   USE becmod,               ONLY : becp, calbec
   USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, &
                                    iunat, iunsat, nwordatwfc
   USE buffers,              ONLY : get_buffer
   IMPLICIT NONE
   !
   ! I/O variables first
   !
   INTEGER :: ipol, jpol, ldim
   REAL (DP) :: dns(ldim,ldim,nspin,nat)
   !
   ! local variable
   !
   INTEGER :: ik,    & ! counter on k points
              ibnd,  & !    "    "  bands
              is,    & !    "    "  spins
              i, na, nt, n, counter, m1, m2, l
   COMPLEX (DP) :: ZDOTC

   INTEGER, ALLOCATABLE :: offset(:)
   ! offset(nat)  ! offset of d electrons of atom d in the natomwfc ordering
   COMPLEX (DP), ALLOCATABLE :: &
                      proj(:,:), wfcatom(:,:), spsi(:,:), dproj(:,:)
   ! proj(natomwfc,nbnd), wfcatom(npwx,natomwfc),
   ! spsi(npwx,nbnd), dproj(natomwfc,nbnd)

   ALLOCATE (offset(nat), proj(natomwfc,nbnd), wfcatom(npwx,natomwfc),  &
             spsi(npwx,nbnd), dproj(natomwfc,nbnd), becp(nkb,nbnd) )

   !
   ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
   !
   counter = 0
   DO na=1,nat
      offset(na) = 0
      nt=ityp(na)
      DO n=1,upf(nt)%nwfc
         IF (upf(nt)%oc(n) >= 0.d0) THEN
            l=upf(nt)%lchi(n)
            IF (l.EQ.Hubbard_l(nt)) offset(na) = counter
            counter = counter + 2 * l + 1
         END IF
      END DO
   END DO

   IF(counter.NE.natomwfc) CALL errore('new_ns','nstart<>counter',1)

   dns(:,:,:,:) = 0.d0
   !
   !    we start a loop on k points
   !
   IF (nks > 1) REWIND (iunigk)

   DO ik = 1, nks
      IF (lsda) current_spin = isk(ik)
      IF (nks > 1) READ (iunigk) igk
      npw = ngk(ik)
      !
      ! now we need the first derivative of proj with respect to
      ! epsilon(ipol,jpol)
      !
      CALL get_buffer (evc, nwordwfc, iunwfc, ik)
      CALL init_us_2 (npw,igk,xk(1,ik),vkb)
      CALL calbec( npw, vkb, evc, becp )

      CALL s_psi  (npwx, npw, nbnd, evc, spsi )
!      CALL atomic_wfc( ik, wfcatom )
! read atomic wfc instead
      CALL davcio(wfcatom,nwordatwfc,iunat,ik,-1)

      dproj(:,:) = (0.d0,0.d0)

      CALL dprojdepsilon(ik,dproj,wfcatom,spsi,ipol,jpol)

      CALL davcio(swfcatom,nwordatwfc,iunsat,ik,-1)

      DO ibnd = 1, nbnd
         DO i=1,natomwfc
            proj(i,ibnd) = ZDOTC(npw,swfcatom(1,i),1,evc(1,ibnd),1)
         ENDDO
      ENDDO

#ifdef __PARA
       CALL reduce(2*natomwfc*nbnd,proj)
#endif
      !
      ! compute the derivative of the occupation numbers (quantities dn(m1,m2))
      ! of the atomic orbitals. They are real quantities as well as n(m1,m2)
      !
      DO na = 1,nat
         nt = ityp(na)
         IF (Hubbard_U(nt).NE.0.d0.OR.Hubbard_alpha(nt).NE.0.d0) THEN        
            DO m1 = 1, 2 * Hubbard_l(nt) + 1
               DO m2 = m1, 2 * Hubbard_l(nt) + 1
                  DO ibnd = 1,nbnd
                     dns(m1,m2,current_spin,na) = dns(m1,m2,current_spin,na) + &
                                             wg(ibnd,ik) *           &
                               DBLE( proj(offset(na)+m1,ibnd) *      &
                              CONJG(dproj(offset(na)+m2,ibnd) ) +    &
                                    dproj(offset(na)+m1,ibnd)*       &
                              CONJG( proj(offset(na)+m2,ibnd) ) )
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
   ! In nspin.eq.1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dns of one spin component
   !
   IF (nspin.EQ.1) dns = 0.5d0 * dns
   !
   ! impose hermeticity of dn_{m1,m2}
   !
   DO na = 1,nat
      nt = ityp(na)
      DO is = 1,nspin
         DO m1 = 1, 2 * Hubbard_l(nt) + 1
            DO m2 = m1+1, 2 * Hubbard_l(nt) + 1
               dns(m2,m1,is,na) = dns(m1,m2,is,na)
            END DO
         END DO
      END DO
   END DO

   DEALLOCATE (offset, proj, wfcatom, spsi, dproj, becp )

   RETURN
END SUBROUTINE dndepsilon
