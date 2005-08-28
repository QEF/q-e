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
SUBROUTINE dndtau_of_k(dns,ldim,offset,proj,wfcatom,spsi,alpha,ipol,ik)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the derivative of the ns with respect to the ionic
   ! displacement u(alpha,ipol) used to obtain the Hubbard contribution to the
   ! atomic forces.
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE basis,                ONLY : natomwfc
   USE lsda_mod,             ONLY : nspin, current_spin
   USE ldaU,                 ONLY : Hubbard_U, Hubbard_alpha
   USE wvfct,                ONLY : nbnd, npwx, npw, wg
   
   IMPLICIT NONE

   INTEGER ::  alpha, ipol, &
               ik, ibnd, is, na, nt, m1, m2
   INTEGER :: ldim
   REAL (DP) :: &
             dns(ldim,ldim,nspin,nat)
   INTEGER :: offset(nat)
   ! offset(nat): offset of d electrons of atom d in the natomwfc ordering
   COMPLEX (DP) :: &
             proj(natomwfc,nbnd), wfcatom(npwx,natomwfc), spsi(npwx,nbnd)
!   COMPLEX (DP) :: ZDOTC, c_one, c_zero
   COMPLEX (DP), ALLOCATABLE :: &
                      dproj(:,:)
   !                  dproj(natomwfc,nbnd)

   CALL start_clock('dndtau_of_k')

   ALLOCATE ( dproj(natomwfc,nbnd) )

   dns(:,:,:,:) = 0.d0
   !
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
                              DBLE(  proj(offset(na)+m1,ibnd)  *   &
                             CONJG(dproj(offset(na)+m2,ibnd))  +   &
                                    dproj(offset(na)+m1,ibnd)  *   &
                             CONJG(proj(offset(na)+m2,ibnd)) )
               END DO
            END DO
         END DO
      END IF
   END DO
   !
   ! In nspin.eq.1 k-point weight wg is normalized to 2 el/band 
   ! in the whole BZ but we are interested in dns of one spin component
   !
   IF (nspin.EQ.1) dns = 0.5d0 * dns
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

   DEALLOCATE ( dproj ) 

   CALL stop_clock('dndtau_of_k')
   RETURN
END SUBROUTINE dndtau_of_k

