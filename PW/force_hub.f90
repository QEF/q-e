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
   USE kinds,     ONLY : DP
   USE ions_base, ONLY : nat, ityp
   USE cell_base, ONLY : at, bg
   USE ldaU,      ONLY : hubbard_lmax, hubbard_l, hubbard_u, &
                         hubbard_alpha, ns, U_projection
   USE lsda_mod,  ONLY : nspin
   USE symme,     ONLY : s, nsym, irt
   USE io_files,  ONLY : prefix, iunocc
   USE wvfct,     ONLY : gamma_only
   USE mp_global, ONLY : me_pool, my_pool_id

   IMPLICIT NONE
   REAL (kind=DP) :: forceh(3,nat)  ! output: the Hubbard forces

   INTEGER :: alpha, na, nt, is, m1, m2, ipol, ldim

   LOGICAL ::  exst

   REAL (kind=DP), ALLOCATABLE :: dns(:,:,:,:)
   !       dns(ldim,ldim,nspin,nat) ! the derivative of the atomic occupations

   IF (U_projection .NE. "atomic") CALL errore("force_hub", &
                   " forces for this U_projection_type not implemented",1)
   IF (gamma_only) CALL errore('force_huh',&
                   ' LDA+U, forces AND gamma-only not implemented yet',1)
   ldim= 2 * Hubbard_lmax + 1
   ALLOCATE(dns(ldim,ldim,nspin,nat))
   forceh(:,:) = 0.d0
   dns(:,:,:,:) = 0.d0

   IF ( me_pool == 0 .AND. my_pool_id == 0 ) THEN

      CALL seqopn (iunocc, TRIM(prefix)//'.occup', 'formatted', exst)
      READ(iunocc,*) ns
      CLOSE(unit=iunocc,status='keep')

   END IF

   DO ipol = 1,3
      DO alpha = 1,nat                 ! the displaced atom
         CALL dndtau(dns,ldim,alpha,ipol)
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

   DEALLOCATE(dns)
   
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
